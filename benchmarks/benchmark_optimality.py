#!/usr/bin/env python3
"""Synthetic beam-size agreement benchmark against exhaustive optimization."""

from __future__ import annotations

import argparse
import random
import statistics

from mountain_centroid.beam import beam_mountain_centroid
from mountain_centroid.sequence import MIN_HAIRPIN_LENGTH, can_pair


def exact_squared_error(sequence: str, expected_heights: list[float]) -> float:
    """Exhaust all valid structures; intended only for short benchmarks."""
    n = len(sequence)
    optimum = float("inf")

    def visit(index: int, stack: tuple[int, ...], cost: float) -> None:
        nonlocal optimum
        if cost >= optimum:
            return
        if index == n:
            if not stack:
                optimum = cost
            return

        def contribution(depth: int) -> float:
            if index == n - 1:
                return 0.0
            return (depth - expected_heights[index]) ** 2

        visit(index + 1, stack, cost + contribution(len(stack)))
        if index + MIN_HAIRPIN_LENGTH + 1 < n:
            visit(
                index + 1,
                stack + (index,),
                cost + contribution(len(stack) + 1),
            )
        if stack:
            left = stack[-1]
            if (
                index - left - 1 >= MIN_HAIRPIN_LENGTH
                and can_pair(sequence[left], sequence[index])
            ):
                visit(
                    index + 1,
                    stack[:-1],
                    cost + contribution(len(stack) - 1),
                )

    visit(0, (), 0.0)
    return optimum


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", nargs="+", type=int, default=[8, 10, 12])
    parser.add_argument("--beam-sizes", nargs="+", type=int, default=[5, 10, 25, 100])
    parser.add_argument("--instances", type=int, default=100)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    rng = random.Random(args.seed)
    print(
        "length,beam_size,instances,exact_match_rate,"
        "mean_absolute_gap,max_absolute_gap"
    )
    for length in args.lengths:
        instances = []
        for _ in range(args.instances):
            sequence = "".join(rng.choice("ACGU") for _ in range(length))
            expected_heights = [
                rng.uniform(0.0, min(cut, length - cut, 4))
                for cut in range(1, length)
            ]
            optimum = exact_squared_error(sequence, expected_heights)
            instances.append((sequence, expected_heights, optimum))

        for beam_size in args.beam_sizes:
            gaps = []
            for sequence, expected_heights, optimum in instances:
                result = beam_mountain_centroid(
                    sequence,
                    expected_heights,
                    beam_size=beam_size,
                )
                gaps.append(max(0.0, result.squared_error - optimum))
            exact_matches = sum(gap <= 1e-9 for gap in gaps)
            print(
                f"{length},{beam_size},{len(gaps)},"
                f"{exact_matches / len(gaps):.6f},"
                f"{statistics.mean(gaps):.6f},{max(gaps):.6f}"
            )


if __name__ == "__main__":
    main()
