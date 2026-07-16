#!/usr/bin/env python3
"""Reproducible time and memory benchmark for the exact constrained solver."""

from __future__ import annotations

import argparse
import random
import statistics
import time
import tracemalloc

from mountain_centroid.exact import exact_mountain_centroid


def make_instance(
    length: int,
    profile: str,
    seed: int,
) -> tuple[str, list[float]]:
    """Create one deterministic synthetic sequence and feasible-scale profile."""
    rng = random.Random(seed)
    if profile == "random":
        sequence = "".join(rng.choice("ACGU") for _ in range(length))
    elif profile == "pair_dense":
        sequence = ("GC" * ((length + 1) // 2))[:length]
    else:  # pragma: no cover - argparse restricts the choices
        raise ValueError(f"Unknown profile: {profile}")

    expected_heights = [
        rng.uniform(0.0, min(cut, length - cut, 8))
        for cut in range(1, length)
    ]
    return sequence, expected_heights


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--length", type=int, default=30)
    parser.add_argument("--instances", type=int, default=5)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument(
        "--profiles",
        nargs="+",
        choices=("random", "pair_dense"),
        default=("random", "pair_dense"),
    )
    parser.add_argument(
        "--measure-python-heap",
        action="store_true",
        help="run a separate traced pass for Python allocation peak",
    )
    args = parser.parse_args()

    if args.length < 1:
        raise ValueError("length must be positive")
    if args.instances < 1:
        raise ValueError("instances must be positive")

    print(
        "profile,length,instance,seconds,states,partner_transitions,"
        "effective_depth_levels,python_heap_peak_mib"
    )
    elapsed_by_profile: dict[str, list[float]] = {}
    heap_by_profile: dict[str, list[float]] = {}

    for profile_index, profile in enumerate(args.profiles):
        elapsed_by_profile[profile] = []
        heap_by_profile[profile] = []
        for instance in range(args.instances):
            instance_seed = args.seed + profile_index * 1_000_000 + instance
            sequence, expected_heights = make_instance(
                args.length,
                profile,
                instance_seed,
            )

            started = time.perf_counter()
            result = exact_mountain_centroid(sequence, expected_heights)
            elapsed = time.perf_counter() - started
            elapsed_by_profile[profile].append(elapsed)

            heap_peak_mib: float | None = None
            if args.measure_python_heap:
                tracemalloc.start()
                exact_mountain_centroid(sequence, expected_heights)
                _, heap_peak = tracemalloc.get_traced_memory()
                tracemalloc.stop()
                heap_peak_mib = heap_peak / 1024**2
                heap_by_profile[profile].append(heap_peak_mib)

            diagnostics = result.diagnostics
            heap_field = "" if heap_peak_mib is None else f"{heap_peak_mib:.6f}"
            print(
                f"{profile},{args.length},{instance},{elapsed:.9f},"
                f"{diagnostics.states_evaluated},"
                f"{diagnostics.partner_transitions_evaluated},"
                f"{diagnostics.effective_depth_levels},{heap_field}"
            )

    for profile in args.profiles:
        elapsed = elapsed_by_profile[profile]
        message = (
            f"# {profile}: median_seconds={statistics.median(elapsed):.9f},"
            f" max_seconds={max(elapsed):.9f}"
        )
        if heap_by_profile[profile]:
            heaps = heap_by_profile[profile]
            message += (
                f", median_python_heap_mib={statistics.median(heaps):.6f},"
                f" max_python_heap_mib={max(heaps):.6f}"
            )
        print(message)


if __name__ == "__main__":
    main()
