"""Exact projection onto unconstrained nonnegative mountain paths."""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Sequence


@dataclass(frozen=True, slots=True)
class RelaxedResult:
    """Exact relaxed Mountain Centroid result."""

    structure: str
    pairs: tuple[tuple[int, int], ...]
    heights: tuple[int, ...]
    squared_error: float


def relaxed_mountain_centroid(expected_heights: Sequence[float]) -> RelaxedResult:
    """Return the exact closest nonnegative unit-step mountain path.

    This relaxation does not enforce sequence complementarity or a minimum
    hairpin length.  For sequence length n, the dynamic program has O(n^2)
    worst-case time and traceback space because the reachable height is O(n).
    More precisely, its time is O(nH), where H is the maximum allowed height.
    """
    mu = tuple(float(value) for value in expected_heights)
    if any(not math.isfinite(value) for value in mu):
        raise ValueError("Expected mountain heights must be finite")
    n = len(mu) + 1
    infinity = float("inf")

    previous = [0.0]
    back: list[list[int]] = [[] for _ in range(n + 1)]
    for position in range(1, n + 1):
        maximum_height = min(position, n - position)
        previous_maximum = min(position - 1, n - position + 1)
        expected = 0.0 if position == n else mu[position - 1]
        current = [infinity] * (maximum_height + 1)
        back[position] = [-1] * (maximum_height + 1)
        for height in range(maximum_height + 1):
            best = (infinity, -1)
            for old_height in (height, height - 1, height + 1):
                if 0 <= old_height <= previous_maximum:
                    candidate = (previous[old_height], old_height)
                    if candidate < best:
                        best = candidate
            contribution = 0.0 if position == n else (height - expected) ** 2
            current[height] = best[0] + contribution
            back[position][height] = best[1]
        previous = current

    heights = [0] * (n + 1)
    for position in range(n, 0, -1):
        heights[position - 1] = back[position][heights[position]]

    characters: list[str] = []
    open_positions: list[int] = []
    pairs: list[tuple[int, int]] = []
    for position in range(1, n + 1):
        difference = heights[position] - heights[position - 1]
        if difference == 1:
            characters.append("(")
            open_positions.append(position)
        elif difference == -1:
            characters.append(")")
            left = open_positions.pop()
            pairs.append((left, position))
        else:
            characters.append(".")

    return RelaxedResult(
        structure="".join(characters),
        pairs=tuple(sorted(pairs)),
        heights=tuple(heights),
        squared_error=float(previous[0]),
    )
