#!/usr/bin/env python3
"""Benchmark the production C++ sequence-constrained route."""

from __future__ import annotations

import argparse
from pathlib import Path
import random
import resource
import statistics
import time

from mountain_centroid.cpp_constrained import (
    cpp_sequence_constrained_mountain_centroid,
    default_cpp_constrained_path,
)


def make_instance(length: int, profile: str, seed: int) -> tuple[str, list[float]]:
    rng = random.Random(seed)
    if profile == "random":
        sequence = "".join(rng.choice("ACGU") for _ in range(length))
    elif profile == "pair_dense":
        sequence = ("GC" * ((length + 1) // 2))[:length]
    else:  # pragma: no cover
        raise ValueError(f"Unknown profile: {profile}")
    expected_heights = [
        rng.uniform(0.0, min(cut, length - cut, 8))
        for cut in range(1, length)
    ]
    return sequence, expected_heights


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lengths",
        nargs="+",
        type=int,
        default=(30, 50, 100, 150, 200, 300),
    )
    parser.add_argument("--instances", type=int, default=3)
    parser.add_argument("--seed", type=int, default=20260717)
    parser.add_argument(
        "--profiles",
        nargs="+",
        choices=("random", "pair_dense"),
        default=("random", "pair_dense"),
    )
    parser.add_argument(
        "--executable", type=Path, default=default_cpp_constrained_path()
    )
    args = parser.parse_args()

    if args.instances < 1:
        raise ValueError("instances must be positive")
    if any(length < 1 for length in args.lengths):
        raise ValueError("lengths must be positive")
    if list(args.lengths) != sorted(args.lengths):
        raise ValueError("lengths must be ascending for cumulative RSS reporting")

    print(
        "profile,length,instance,seconds,states,partner_transitions,"
        "effective_depth_levels,cumulative_child_peak_rss_mib"
    )
    for profile_index, profile in enumerate(args.profiles):
        for length in args.lengths:
            elapsed_values = []
            for instance in range(args.instances):
                instance_seed = (
                    args.seed + profile_index * 1_000_000 + length * 1_000 + instance
                )
                sequence, expected_heights = make_instance(
                    length,
                    profile,
                    instance_seed,
                )
                started = time.perf_counter()
                result = cpp_sequence_constrained_mountain_centroid(
                    sequence,
                    expected_heights,
                    executable=args.executable,
                )
                elapsed = time.perf_counter() - started
                elapsed_values.append(elapsed)
                peak_rss_mib = (
                    resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss / 1024
                )
                diagnostics = result.diagnostics
                print(
                    f"{profile},{length},{instance},{elapsed:.9f},"
                    f"{diagnostics.states_evaluated},"
                    f"{diagnostics.partner_transitions_evaluated},"
                    f"{diagnostics.effective_depth_levels},{peak_rss_mib:.6f}"
                )
            print(
                f"# {profile} length {length}: "
                f"median_seconds={statistics.median(elapsed_values):.9f}, "
                f"max_seconds={max(elapsed_values):.9f}"
            )


if __name__ == "__main__":
    main()
