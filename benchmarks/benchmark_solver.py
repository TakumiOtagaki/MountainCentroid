#!/usr/bin/env python3
"""Reproducible synthetic scaling benchmark for the beam solver."""

from __future__ import annotations

import argparse
import random
import statistics
import time

from mountain_centroid.beam import beam_mountain_centroid


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--lengths", nargs="+", type=int, default=[100, 300, 1000])
    parser.add_argument("--beam-size", type=int, default=100)
    parser.add_argument("--repeats", type=int, default=3)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    if args.repeats < 1:
        raise ValueError("repeats must be positive")

    rng = random.Random(args.seed)
    print("length,beam_size,repeats,median_seconds,min_seconds")
    for length in args.lengths:
        if length < 1:
            raise ValueError("lengths must be positive")
        sequence = "".join(rng.choice("ACGU") for _ in range(length))
        expected_heights = [
            rng.uniform(0.0, min(cut, length - cut, 8))
            for cut in range(1, length)
        ]
        elapsed = []
        for _ in range(args.repeats):
            started = time.perf_counter()
            beam_mountain_centroid(
                sequence,
                expected_heights,
                beam_size=args.beam_size,
            )
            elapsed.append(time.perf_counter() - started)
        print(
            f"{length},{args.beam_size},{args.repeats},"
            f"{statistics.median(elapsed):.6f},{min(elapsed):.6f}"
        )


if __name__ == "__main__":
    main()
