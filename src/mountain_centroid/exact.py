"""Exact sequence-constrained mountain-profile projection."""

from __future__ import annotations

from dataclasses import dataclass
from functools import cache
import math
from typing import Sequence

from .sequence import MIN_HAIRPIN_LENGTH, can_pair, normalise_sequence


@dataclass(frozen=True, slots=True)
class ExactDiagnostics:
    """Dynamic-programming work performed for one exact projection."""

    states_evaluated: int
    partner_transitions_evaluated: int
    maximum_external_depth: int

    @property
    def effective_depth_levels(self) -> int:
        """Number of external-depth levels reached, including depth zero."""
        return self.maximum_external_depth + 1


@dataclass(frozen=True, slots=True)
class ExactResult:
    """Exact sequence-constrained Mountain Centroid result."""

    structure: str
    pairs: tuple[tuple[int, int], ...]
    heights: tuple[int, ...]
    squared_error: float
    diagnostics: ExactDiagnostics


def exact_mountain_centroid(
    sequence: str,
    expected_heights: Sequence[float],
) -> ExactResult:
    """Return the exact sequence-valid projection of an expected profile.

    The interval state ``F(i, j, d)`` stores the best cost within ``[i, j]``
    when ``d`` outside pairs enclose the interval.  Only states reachable from
    ``F(0, n - 1, 0)`` are memoized.  If ``D_eff`` external-depth levels are
    reached, the implementation uses O(D_eff n^3) time and O(D_eff n^2)
    memoization space.  These become O(n^4) time and O(n^3) space in the worst
    case because D_eff is O(n).

    Returned pairs are Watson--Crick or GU wobble pairs, satisfy the package's
    minimum hairpin length, and are pseudoknot-free.
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

    legal_partners: tuple[tuple[int, ...], ...] = tuple(
        tuple(
            right
            for right in range(left + MIN_HAIRPIN_LENGTH + 1, n)
            if can_pair(sequence[left], sequence[right])
        )
        for left in range(n)
    )

    states_evaluated = 0
    partner_transitions_evaluated = 0
    maximum_external_depth = 0

    def cut_cost(position: int, depth: int) -> float:
        if position == n - 1:
            return 0.0
        return (depth - mu[position]) ** 2

    def subproblem_cost(left: int, right: int, depth: int) -> float:
        if left > right:
            return 0.0
        return solve(left, right, depth)[0]

    @cache
    def solve(left: int, right: int, depth: int) -> tuple[float, int]:
        # The root recursion guarantees this positional feasibility bound:
        # every outside enclosing pair needs one endpoint on either side.
        if depth > min(left, n - 1 - right):
            raise RuntimeError("Reached an infeasible external-depth state")

        nonlocal states_evaluated
        nonlocal partner_transitions_evaluated
        nonlocal maximum_external_depth
        states_evaluated += 1
        maximum_external_depth = max(maximum_external_depth, depth)

        best_cost = cut_cost(left, depth) + subproblem_cost(
            left + 1,
            right,
            depth,
        )
        best_partner = -1

        for partner in legal_partners[left]:
            if partner > right:
                break
            partner_transitions_evaluated += 1
            candidate = (
                cut_cost(left, depth + 1)
                + subproblem_cost(left + 1, partner - 1, depth + 1)
                + cut_cost(partner, depth)
                + subproblem_cost(partner + 1, right, depth)
            )
            # Strict replacement gives deterministic sparse tie-breaking:
            # unpaired precedes paired, then partners are considered leftmost
            # first.
            if candidate < best_cost:
                best_cost = candidate
                best_partner = partner

        return best_cost, best_partner

    squared_error, _ = solve(0, n - 1, 0)

    characters = ["."] * n
    pairs: list[tuple[int, int]] = []

    def traceback(left: int, right: int, depth: int) -> None:
        if left > right:
            return
        _, partner = solve(left, right, depth)
        if partner < 0:
            traceback(left + 1, right, depth)
            return
        characters[left] = "("
        characters[partner] = ")"
        pairs.append((left + 1, partner + 1))
        traceback(left + 1, partner - 1, depth + 1)
        traceback(partner + 1, right, depth)

    traceback(0, n - 1, 0)
    pairs.sort()

    heights = [0]
    depth = 0
    for character in characters:
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        heights.append(depth)

    return ExactResult(
        structure="".join(characters),
        pairs=tuple(pairs),
        heights=tuple(heights),
        squared_error=float(squared_error),
        diagnostics=ExactDiagnostics(
            states_evaluated=states_evaluated,
            partner_transitions_evaluated=partner_transitions_evaluated,
            maximum_external_depth=maximum_external_depth,
        ),
    )
