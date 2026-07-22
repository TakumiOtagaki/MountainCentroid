"""Topology descriptors for pseudoknot-free RNA secondary structures.

The lossless representation used here has one vertex for every loop region,
including the exterior region, and one weighted edge for every maximal stack
of consecutive base pairs.  Edge weights are the numbers of base pairs in the
stems.  This is a tree for every standard pseudoknot-free dot-bracket string.

The optional ``rag`` mode is a paper-inspired coarse graining based on the
RNA-As-Graphs conventions used by Gopal et al. (2014): isolated base pairs are
ignored and degree-two loops containing at most one unpaired nucleotide are
suppressed.  Pair types are deliberately not filtered, so geometry-only
Mountain Centroid structures can be compared with sequence-valid predictions
using the same mapping.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
import math
from typing import Literal


CoarseGraining = Literal["lossless", "rag"]


@dataclass(frozen=True)
class TreeEdge:
    """A stem joining two loop vertices."""

    left: int
    right: int
    stem_length: int


@dataclass(frozen=True)
class RNATreeGraph:
    """Loop--stem tree derived from one RNA secondary structure."""

    edges: tuple[TreeEdge, ...]
    unpaired_nucleotides: tuple[int, ...]
    exterior_vertex: int
    coarse_graining: CoarseGraining

    @property
    def vertex_count(self) -> int:
        return len(self.unpaired_nucleotides)

    @property
    def stem_count(self) -> int:
        return len(self.edges)

    def adjacency(self) -> tuple[tuple[tuple[int, int], ...], ...]:
        """Return ``(neighbor, stem_length)`` entries for every vertex."""
        neighbors: list[list[tuple[int, int]]] = [
            [] for _ in range(self.vertex_count)
        ]
        for edge in self.edges:
            neighbors[edge.left].append((edge.right, edge.stem_length))
            neighbors[edge.right].append((edge.left, edge.stem_length))
        return tuple(tuple(items) for items in neighbors)

    def degrees(self) -> tuple[int, ...]:
        return tuple(len(items) for items in self.adjacency())


@dataclass(frozen=True)
class TopologyDescriptors:
    """Scalar descriptors computed from one secondary structure."""

    nucleotide_count: int
    base_pair_count: int
    stem_count: int
    multiloop_count: int
    vertex_count: int
    v0: int
    v1: int
    v2: int
    v3: int
    v4: int
    v_ge4: int
    v1_over_v3: float
    v_ge4_over_v3: float
    graph_diameter: int
    graph_radius: int
    average_graph_distance: float
    wiener_index: int
    tree_radius_of_gyration: float
    maximum_ladder_distance: int
    average_ladder_distance: float
    average_stem_length: float
    maximum_stem_length: int
    maximum_mountain_height: int
    mountain_area: int
    average_base_pair_span: float
    maximum_base_pair_span: int
    proximal_pair_fraction_100nt: float
    leaf_fraction: float
    branch_vertex_fraction: float
    degree_entropy: float

    def as_dict(self) -> dict[str, int | float]:
        return dict(self.__dict__)


def _pair_table(structure: str) -> list[int]:
    pair = [-1] * len(structure)
    stack: list[int] = []
    for position, character in enumerate(structure):
        if character == "(":
            stack.append(position)
        elif character == ")":
            if not stack:
                raise ValueError(
                    f"Unbalanced closing bracket at position {position + 1}"
                )
            partner = stack.pop()
            pair[position] = partner
            pair[partner] = position
        elif character != ".":
            raise ValueError(f"Unknown dot-bracket character: {character}")
    if stack:
        positions = ", ".join(str(position + 1) for position in stack)
        raise ValueError(f"Unbalanced opening brackets at positions {positions}")
    return pair


def _lossless_tree(structure: str) -> RNATreeGraph:
    pair = _pair_table(structure)
    unpaired: list[int] = [0]
    edges: list[TreeEdge] = []

    def add_loop(start: int, stop: int, vertex: int) -> None:
        position = start
        loop_unpaired = 0
        while position < stop:
            partner = pair[position]
            if partner < 0:
                loop_unpaired += 1
                position += 1
                continue
            if partner < position:
                raise ValueError("Unexpected closing pair while constructing tree")

            stem_length = 1
            while (
                position + stem_length < partner - stem_length
                and pair[position + stem_length] == partner - stem_length
            ):
                stem_length += 1

            child = len(unpaired)
            unpaired.append(0)
            edges.append(TreeEdge(vertex, child, stem_length))
            add_loop(
                position + stem_length,
                partner - stem_length + 1,
                child,
            )
            position = partner + 1
        unpaired[vertex] = loop_unpaired

    add_loop(0, len(structure), 0)
    return RNATreeGraph(
        edges=tuple(edges),
        unpaired_nucleotides=tuple(unpaired),
        exterior_vertex=0,
        coarse_graining="lossless",
    )


def _rag_coarse_grain(graph: RNATreeGraph) -> RNATreeGraph:
    unpaired = {
        vertex: count for vertex, count in enumerate(graph.unpaired_nucleotides)
    }
    adjacency: dict[int, dict[int, int]] = {
        vertex: {} for vertex in range(graph.vertex_count)
    }
    for edge in graph.edges:
        adjacency[edge.left][edge.right] = edge.stem_length
        adjacency[edge.right][edge.left] = edge.stem_length
    exterior = graph.exterior_vertex

    def contract_edge(left: int, right: int, ignored_pair: bool) -> None:
        nonlocal exterior
        adjacency[left].pop(right)
        adjacency[right].pop(left)
        if ignored_pair:
            unpaired[left] += 2
        unpaired[left] += unpaired[right]
        for neighbor, weight in list(adjacency[right].items()):
            adjacency[neighbor].pop(right)
            adjacency[neighbor][left] = weight
            adjacency[left][neighbor] = weight
        del adjacency[right]
        del unpaired[right]
        if exterior == right:
            exterior = left

    # RAG does not represent isolated single-base-pair stems.
    while True:
        isolated = next(
            (
                (left, right)
                for left, neighbors in adjacency.items()
                for right, weight in neighbors.items()
                if left < right and weight == 1
            ),
            None,
        )
        if isolated is None:
            break
        contract_edge(*isolated, ignored_pair=True)

    # RAG ignores one-residue bulges and internal loops.  Suppression joins
    # the adjacent stacks into one coarse-grained stem.
    while True:
        vertex = next(
            (
                candidate
                for candidate, neighbors in adjacency.items()
                if candidate != exterior
                and len(neighbors) == 2
                and unpaired[candidate] <= 1
            ),
            None,
        )
        if vertex is None:
            break
        (left, left_weight), (right, right_weight) = adjacency[vertex].items()
        adjacency[left].pop(vertex)
        adjacency[right].pop(vertex)
        adjacency[left][right] = left_weight + right_weight
        adjacency[right][left] = left_weight + right_weight
        del adjacency[vertex]
        del unpaired[vertex]

    old_vertices = sorted(adjacency)
    index = {old: new for new, old in enumerate(old_vertices)}
    edges = tuple(
        TreeEdge(index[left], index[right], weight)
        for left in old_vertices
        for right, weight in sorted(adjacency[left].items())
        if left < right
    )
    return RNATreeGraph(
        edges=edges,
        unpaired_nucleotides=tuple(unpaired[old] for old in old_vertices),
        exterior_vertex=index[exterior],
        coarse_graining="rag",
    )


def rna_tree_graph(
    structure: str,
    *,
    coarse_graining: CoarseGraining = "lossless",
) -> RNATreeGraph:
    """Convert standard dot-bracket notation into a loop--stem tree."""
    if coarse_graining not in ("lossless", "rag"):
        raise ValueError(f"Unknown coarse graining: {coarse_graining}")
    graph = _lossless_tree(structure)
    return graph if coarse_graining == "lossless" else _rag_coarse_grain(graph)


def _tree_distance_sums(
    graph: RNATreeGraph,
) -> tuple[int, int]:
    """Return unweighted and stem-length-weighted Wiener indices."""
    vertex_count = graph.vertex_count
    if vertex_count <= 1:
        return 0, 0
    adjacency = graph.adjacency()
    parent = [-1] * vertex_count
    parent_weight = [0] * vertex_count
    order = [graph.exterior_vertex]
    for vertex in order:
        for neighbor, weight in adjacency[vertex]:
            if neighbor == parent[vertex]:
                continue
            parent[neighbor] = vertex
            parent_weight[neighbor] = weight
            order.append(neighbor)
    if len(order) != vertex_count:
        raise ValueError("RNA graph is disconnected")

    subtree_size = [1] * vertex_count
    unweighted = 0
    weighted = 0
    for vertex in reversed(order[1:]):
        size = subtree_size[vertex]
        separated_pairs = size * (vertex_count - size)
        unweighted += separated_pairs
        weighted += separated_pairs * parent_weight[vertex]
        subtree_size[parent[vertex]] += size
    return unweighted, weighted


def _distances_from(
    graph: RNATreeGraph,
    start: int,
    *,
    weighted: bool,
) -> list[int]:
    adjacency = graph.adjacency()
    distances = [-1] * graph.vertex_count
    distances[start] = 0
    stack = [start]
    while stack:
        vertex = stack.pop()
        for neighbor, stem_length in adjacency[vertex]:
            if distances[neighbor] >= 0:
                continue
            increment = stem_length if weighted else 1
            distances[neighbor] = distances[vertex] + increment
            stack.append(neighbor)
    return distances


def _diameter(graph: RNATreeGraph, *, weighted: bool) -> int:
    if graph.vertex_count <= 1:
        return 0
    first = _distances_from(graph, 0, weighted=weighted)
    endpoint = max(range(graph.vertex_count), key=first.__getitem__)
    return max(_distances_from(graph, endpoint, weighted=weighted))


def topology_descriptors(
    structure: str,
    *,
    coarse_graining: CoarseGraining = "lossless",
) -> TopologyDescriptors:
    """Compute topology descriptors for one dot-bracket structure.

    ``tree_radius_of_gyration`` is the Kramers ideal-tree value in units of
    the common edge length ``b``.  ``maximum_ladder_distance`` and
    ``average_ladder_distance`` instead weight each tree edge by its number of
    base pairs.
    """
    graph = rna_tree_graph(structure, coarse_graining=coarse_graining)
    degrees = graph.degrees()
    degree_counts = Counter(degrees)
    vertex_count = graph.vertex_count
    stem_lengths = [edge.stem_length for edge in graph.edges]
    unweighted_wiener, weighted_wiener = _tree_distance_sums(graph)
    pair_table = _pair_table(structure)
    pairs = [
        (left, right)
        for left, right in enumerate(pair_table)
        if right > left
    ]
    pair_spans = [right - left for left, right in pairs]

    height = 0
    mountain_area = 0
    maximum_height = 0
    for position, character in enumerate(structure):
        if character == "(":
            height += 1
        elif character == ")":
            height -= 1
        maximum_height = max(maximum_height, height)
        if position < len(structure) - 1:
            mountain_area += height

    v3 = degree_counts[3]
    v_ge4 = sum(count for degree, count in degree_counts.items() if degree >= 4)
    multiloops = sum(
        count for degree, count in degree_counts.items() if degree >= 3
    )
    pair_count = len(pairs)
    degree_entropy = 0.0
    if vertex_count:
        for count in degree_counts.values():
            probability = count / vertex_count
            degree_entropy -= probability * math.log(probability)

    pair_denominator = vertex_count * (vertex_count - 1) / 2
    return TopologyDescriptors(
        nucleotide_count=len(structure),
        base_pair_count=pair_count,
        stem_count=graph.stem_count,
        multiloop_count=multiloops,
        vertex_count=vertex_count,
        v0=degree_counts[0],
        v1=degree_counts[1],
        v2=degree_counts[2],
        v3=v3,
        v4=degree_counts[4],
        v_ge4=v_ge4,
        v1_over_v3=degree_counts[1] / v3 if v3 else math.nan,
        v_ge4_over_v3=v_ge4 / v3 if v3 else math.nan,
        graph_diameter=_diameter(graph, weighted=False),
        graph_radius=math.ceil(_diameter(graph, weighted=False) / 2),
        average_graph_distance=(
            unweighted_wiener / pair_denominator if pair_denominator else 0.0
        ),
        wiener_index=unweighted_wiener,
        tree_radius_of_gyration=(
            math.sqrt(unweighted_wiener / (vertex_count * vertex_count))
            if vertex_count
            else 0.0
        ),
        maximum_ladder_distance=_diameter(graph, weighted=True),
        average_ladder_distance=(
            weighted_wiener / pair_denominator if pair_denominator else 0.0
        ),
        average_stem_length=(
            sum(stem_lengths) / len(stem_lengths) if stem_lengths else 0.0
        ),
        maximum_stem_length=max(stem_lengths, default=0),
        maximum_mountain_height=maximum_height,
        mountain_area=mountain_area,
        average_base_pair_span=(
            sum(pair_spans) / pair_count if pair_count else 0.0
        ),
        maximum_base_pair_span=max(pair_spans, default=0),
        proximal_pair_fraction_100nt=(
            sum(span <= 100 for span in pair_spans) / pair_count
            if pair_count
            else 0.0
        ),
        leaf_fraction=(degree_counts[1] / vertex_count if vertex_count else 0.0),
        branch_vertex_fraction=(
            multiloops / vertex_count if vertex_count else 0.0
        ),
        degree_entropy=degree_entropy,
    )
