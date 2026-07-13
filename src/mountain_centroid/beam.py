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


def beam_mountain_centroid(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    beam_size: int = 100,
) -> BeamResult:
    """Minimize squared mountain loss with sequence-valid beam search.

    The search scans the sequence from left to right. Its state is the stack of
    currently open base pairs, which makes every emitted structure
    pseudoknot-free. Stack nodes are persistent parent pointers, so opening and
    closing a pair are O(1). Candidate sorting gives O(n B log B) time and
    traceback storage is O(nB). With a fixed beam size B, runtime is O(n) in
    sequence length.

    Beam pruning makes the returned structure an approximation to the
    sequence-constrained Fréchet mean. Setting a sufficiently large beam retains
    every prefix state and recovers the exact optimum for small instances.
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

        candidates = sorted(best_by_stack.values(), key=lambda item: item[:2])
        beam = candidates[:beam_size]

        # Keep the all-unpaired prefix as a guaranteed feasible fallback. This
        # reserves at most one beam entry and prevents dead ends near sequence end.
        empty_state = best_by_stack.get(0)
        if empty_state is not None and all(state[2] != 0 for state in beam):
            if beam_size == 1:
                beam = [empty_state]
            else:
                beam[-1] = empty_state
                beam.sort(key=lambda item: item[:2])

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
