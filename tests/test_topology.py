import math

import pytest

from mountain_centroid.topology import rna_tree_graph, topology_descriptors


def test_single_stem_descriptors() -> None:
    descriptors = topology_descriptors("((...))")

    assert descriptors.nucleotide_count == 7
    assert descriptors.base_pair_count == 2
    assert descriptors.stem_count == 1
    assert descriptors.vertex_count == 2
    assert descriptors.v1 == 2
    assert descriptors.v2 == 0
    assert descriptors.graph_diameter == 1
    assert descriptors.graph_radius == 1
    assert descriptors.average_graph_distance == 1.0
    assert descriptors.tree_radius_of_gyration == 0.5
    assert descriptors.maximum_ladder_distance == 2
    assert descriptors.average_stem_length == 2.0
    assert descriptors.maximum_mountain_height == 2
    assert descriptors.mountain_area == 10


def test_internal_loop_creates_degree_two_vertex() -> None:
    descriptors = topology_descriptors("((..((...))..))")

    assert descriptors.stem_count == 2
    assert descriptors.vertex_count == 3
    assert (descriptors.v1, descriptors.v2, descriptors.v3) == (2, 1, 0)
    assert descriptors.graph_diameter == 2
    assert descriptors.wiener_index == 4
    assert descriptors.average_graph_distance == pytest.approx(4 / 3)
    assert descriptors.tree_radius_of_gyration == pytest.approx(2 / 3)
    assert descriptors.maximum_ladder_distance == 4


def test_multiloop_degree_counts_obey_tree_identity() -> None:
    structure = "((..((...))..((...))..))"
    descriptors = topology_descriptors(structure)

    assert descriptors.stem_count == 3
    assert descriptors.multiloop_count == 1
    assert descriptors.v1 == 3
    assert descriptors.v3 == 1
    assert descriptors.v1_over_v3 == 3.0
    assert descriptors.v_ge4 == 0
    assert descriptors.graph_diameter == 2
    assert descriptors.maximum_ladder_distance == 4


def test_rag_coarse_graining_ignores_one_residue_bulge() -> None:
    structure = "((.((...))))"
    lossless = rna_tree_graph(structure, coarse_graining="lossless")
    rag = rna_tree_graph(structure, coarse_graining="rag")

    assert lossless.vertex_count == 3
    assert [edge.stem_length for edge in lossless.edges] == [2, 2]
    assert rag.vertex_count == 2
    assert len(rag.edges) == 1
    assert rag.edges[0].stem_length == 4


def test_rag_coarse_graining_ignores_isolated_pair() -> None:
    lossless = topology_descriptors("(...)", coarse_graining="lossless")
    rag = topology_descriptors("(...)", coarse_graining="rag")

    assert (lossless.vertex_count, lossless.stem_count) == (2, 1)
    assert (rag.vertex_count, rag.stem_count) == (1, 0)
    assert rag.v0 == 1


def test_unpaired_structure_and_undefined_branch_ratios() -> None:
    descriptors = topology_descriptors("....")

    assert descriptors.vertex_count == 1
    assert descriptors.v0 == 1
    assert descriptors.stem_count == 0
    assert descriptors.graph_diameter == 0
    assert descriptors.tree_radius_of_gyration == 0.0
    assert math.isnan(descriptors.v1_over_v3)


@pytest.mark.parametrize("structure", ["(()", "())", "(..x..)"])
def test_invalid_dot_bracket_is_rejected(structure: str) -> None:
    with pytest.raises(ValueError):
        rna_tree_graph(structure)
