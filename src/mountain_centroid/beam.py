"""Beam-pruned sequence-constrained mountain-profile inference."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence

from .sequence import MIN_HAIRPIN_LENGTH, can_pair, normalise_sequence


@dataclass(frozen=True, slots=True)
class BeamResult:
    """Result of projecting an expected profile onto valid RNA structures."""

    structure: str
    pairs: tuple[tuple[int, int], ...]
    heights: tuple[int, ...]
    squared_error: float


def _bases_available_after_minimum_loop(sequence: str) -> list[int]:
    """Return suffix base masks, used to avoid opening impossible pairs."""
    base_bit = {"A": 1, "C": 2, "G": 4, "U": 8}
    suffix_mask = [0] * (len(sequence) + 1)
    for index in range(len(sequence) - 1, -1, -1):
        suffix_mask[index] = suffix_mask[index + 1] | base_bit[sequence[index]]
    return suffix_mask


def _prune_with_depth_diversity(
    candidates: list[tuple[float, int, int, int]],
    stack_depth: list[int],
    beam_size: int,
) -> list[tuple[float, int, int, int]]:
    """Fill the beam round-robin across reachable mountain depths.

    A purely global prefix-cost beam can discard every deeply nested state long
    before distant pairs become closable. Round-robin pruning across the ten
    best-scoring depth classes preserves the O(B) state budget while retaining
    several alternative stacks at relevant future height ranges.
    """
    ranked = sorted(candidates, key=lambda item: item[:2])
    candidates_at_depth: dict[int, list[tuple[float, int, int, int]]] = {}
    for candidate in ranked:
        depth = stack_depth[candidate[2]]
        candidates_at_depth.setdefault(depth, []).append(candidate)

    depth_order = sorted(
        candidates_at_depth,
        key=lambda depth: candidates_at_depth[depth][0][:2],
    )
    if 0 in depth_order:
        depth_order.remove(0)
        depth_order.insert(0, 0)
    number_of_depths = min(10, beam_size)
    if depth_order and depth_order[0] == 0:
        depth_order = depth_order[:1] + depth_order[1:number_of_depths]
    else:
        depth_order = depth_order[:number_of_depths]

    selected: list[tuple[float, int, int, int]] = []
    rank_within_depth = 0
    while len(selected) < beam_size:
        added = False
        for depth in depth_order:
            group = candidates_at_depth[depth]
            if rank_within_depth >= len(group):
                continue
            selected.append(group[rank_within_depth])
            added = True
            if len(selected) == beam_size:
                break
        if not added:
            break
        rank_within_depth += 1
    selected.sort(key=lambda item: item[:2])
    return selected


def _beam_mountain_centroid_one_direction(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    beam_size: int = 100,
) -> BeamResult:
    """Run the sequence-valid beam search in one scan direction.

    The search scans the sequence from left to right. Its state is the stack of
    currently open base pairs, which makes every emitted structure
    pseudoknot-free. Stack nodes are persistent parent pointers, so opening and
    closing a pair are O(1). Candidate sorting gives O(n B log B) time and
    traceback storage is O(nB). With a fixed beam size B, runtime is O(n) in
    sequence length.

    Beam pruning makes the returned structure an approximation to the
    sequence-constrained Fréchet mean. The implementation is checked against an
    exhaustive oracle on small instances and long-range regression profiles.
    """
    sequence = normalise_sequence(sequence)
    n = len(sequence)
    mu = tuple(float(value) for value in expected_heights)
    if len(mu) != n - 1:
        raise ValueError(
            f"Expected {n - 1} mountain heights for a length-{n} sequence, "
            f"got {len(mu)}"
        )
    if any(not math.isfinite(value) for value in mu):
        raise ValueError("Expected mountain heights must be finite")
    if (
        isinstance(beam_size, bool)
        or not isinstance(beam_size, int)
        or beam_size < 1
    ):
        raise ValueError("beam_size must be a positive integer")

    # Stack 0 is empty. Other nodes store only their top position and parent;
    # no tuple proportional to the nesting depth is copied during transitions.
    stack_parent = [-1]
    stack_top = [-1]
    stack_depth = [0]

    # Trace 0 is empty. Each transition appends one dot-bracket character.
    trace_parent = [-1]
    trace_char = [""]

    # State tuple: (cost, insertion_order, stack_id, trace_id).
    beam: list[tuple[float, int, int, int]] = [(0.0, 0, 0, 0)]
    next_order = 1
    suffix_mask = _bases_available_after_minimum_loop(sequence)
    compatible_mask = {"A": 8, "C": 4, "G": 2 | 8, "U": 1 | 4}

    for index, base in enumerate(sequence):
        position = index + 1
        remaining = n - position
        best_by_stack: dict[int, tuple[float, int, int, int]] = {}

        def add_candidate(
            source_trace: int,
            new_stack: int,
            character: str,
            previous_cost: float,
        ) -> None:
            nonlocal next_order
            depth = stack_depth[new_stack]
            if depth > remaining:
                return
            contribution = 0.0
            if position < n:
                contribution = (depth - mu[index]) ** 2
            trace_id = len(trace_parent)
            trace_parent.append(source_trace)
            trace_char.append(character)
            candidate = (
                previous_cost + contribution,
                next_order,
                new_stack,
                trace_id,
            )
            next_order += 1
            incumbent = best_by_stack.get(new_stack)
            if incumbent is None or candidate[:2] < incumbent[:2]:
                best_by_stack[new_stack] = candidate

        for cost, _, stack_id, trace_id in beam:
            # Prefer an unpaired nucleotide when two paths have exactly the same
            # cost and future state; this supplies a deterministic sparse tie.
            add_candidate(trace_id, stack_id, ".", cost)

            earliest_partner_index = index + MIN_HAIRPIN_LENGTH + 1
            has_future_partner = (
                earliest_partner_index < n
                and suffix_mask[earliest_partner_index]
                & compatible_mask[base]
            )
            if has_future_partner:
                opened_stack = len(stack_parent)
                stack_parent.append(stack_id)
                stack_top.append(index)
                stack_depth.append(stack_depth[stack_id] + 1)
                add_candidate(trace_id, opened_stack, "(", cost)

            if stack_id != 0:
                left_index = stack_top[stack_id]
                hairpin_length = index - left_index - 1
                if (
                    hairpin_length >= MIN_HAIRPIN_LENGTH
                    and can_pair(sequence[left_index], base)
                ):
                    add_candidate(
                        trace_id,
                        stack_parent[stack_id],
                        ")",
                        cost,
                    )

        beam = _prune_with_depth_diversity(
            list(best_by_stack.values()),
            stack_depth,
            beam_size,
        )

    final_state = min(
        (state for state in beam if state[2] == 0),
        key=lambda item: item[:2],
    )
    trace_id = final_state[3]
    characters: list[str] = []
    while trace_id != 0:
        characters.append(trace_char[trace_id])
        trace_id = trace_parent[trace_id]
    structure = "".join(reversed(characters))

    open_positions: list[int] = []
    pairs: list[tuple[int, int]] = []
    heights = [0]
    for position, character in enumerate(structure, start=1):
        if character == "(":
            open_positions.append(position)
        elif character == ")":
            left = open_positions.pop()
            pairs.append((left, position))
        heights.append(len(open_positions))
    pairs.sort()
    squared_error = sum(
        (heights[k] - mu[k - 1]) ** 2 for k in range(1, n)
    )

    return BeamResult(
        structure=structure,
        pairs=tuple(pairs),
        heights=tuple(heights),
        squared_error=float(squared_error),
    )


def beam_mountain_centroid(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    beam_size: int = 100,
) -> BeamResult:
    """Minimize squared mountain loss with bidirectional beam pruning.

    A left-to-right scan can lose a useful set of opening positions before its
    distant closing partners are visible. The equivalent reversed problem has
    the same objective and structural constraints but exposes those partners in
    the opposite order. Running both directions and retaining the lower-loss
    result substantially reduces directional pruning failures while preserving
    O(n B log B) time for a fixed beam size B.

    The returned structure is pseudoknot-free, uses only canonical
    Watson-Crick or GU wobble pairs, and enforces the minimum hairpin length.
    Beam pruning makes it an approximation to the sequence-constrained
    Frechet mean.
    """
    forward = _beam_mountain_centroid_one_direction(
        sequence,
        expected_heights,
        beam_size=beam_size,
    )
    reverse = _beam_mountain_centroid_one_direction(
        sequence[::-1],
        tuple(reversed(expected_heights)),
        beam_size=beam_size,
    )
    reverse_structure = reverse.structure[::-1].translate(
        str.maketrans("()", ")("),
    )
    if reverse.squared_error >= forward.squared_error:
        return forward

    n = len(reverse_structure)
    pairs = tuple(
        sorted((n + 1 - right, n + 1 - left) for left, right in reverse.pairs)
    )
    heights = tuple(reversed(reverse.heights))
    squared_error = sum(
        (heights[k] - float(expected_heights[k - 1])) ** 2
        for k in range(1, n)
    )
    return BeamResult(
        structure=reverse_structure,
        pairs=pairs,
        heights=heights,
        squared_error=float(squared_error),
    )
