#!/usr/bin/env python3
"""Benchmark ViennaRNA preprocessing and Mountain Centroid separately."""

from __future__ import annotations

import argparse
import gc
import platform
import random
import time

import RNA

from mountain_centroid.beam import beam_mountain_centroid
from mountain_centroid.bpp_mu import (
    compute_bpp_vienna,
    mountain_expectation_from_bpp,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--lengths",
        nargs="+",
        type=int,
        default=[50, 100, 200, 400, 800],
    )
    parser.add_argument("--instances", type=int, default=5)
    parser.add_argument("--repeats", type=int, default=1)
    parser.add_argument("--beam-size", type=int, default=100)
    parser.add_argument("--temperature", type=float, default=37.0)
    parser.add_argument("--seed", type=int, default=1)
    parser.add_argument("--no-warmup", action="store_true")
    args = parser.parse_args()

    if args.instances < 1 or args.repeats < 1:
        raise ValueError("instances and repeats must be positive")
    if any(length < 1 for length in args.lengths):
        raise ValueError("lengths must be positive")

    if not args.no_warmup:
        warmup_sequence = "GGGAAACCC"
        warmup_bpp = compute_bpp_vienna(
            warmup_sequence,
            temperature=args.temperature,
        )
        warmup_mu = mountain_expectation_from_bpp(warmup_bpp)
        beam_mountain_centroid(
            warmup_sequence,
            warmup_mu,
            beam_size=args.beam_size,
        )

    print(f"# platform={platform.platform()}")
    print(f"# python={platform.python_version()}")
    print(f"# viennarna={getattr(RNA, '__version__', 'unknown')}")
    print(f"# temperature_c={args.temperature}")
    print(f"# master_seed={args.seed}")
    print(
        "length,instance,repeat,sequence_seed,beam_size,"
        "vienna_bpp_seconds,profile_seconds,solver_seconds,total_seconds,pairs"
    )

    master_rng = random.Random(args.seed)
    for length in args.lengths:
        for instance in range(args.instances):
            sequence_seed = master_rng.randrange(2**63)
            sequence_rng = random.Random(sequence_seed)
            sequence = "".join(
                sequence_rng.choice("ACGU") for _ in range(length)
            )

            for repeat in range(args.repeats):
                total_started = time.perf_counter()

                started = time.perf_counter()
                bpp = compute_bpp_vienna(
                    sequence,
                    temperature=args.temperature,
                )
                bpp_seconds = time.perf_counter() - started

                started = time.perf_counter()
                expected_heights = mountain_expectation_from_bpp(bpp)
                profile_seconds = time.perf_counter() - started

                started = time.perf_counter()
                result = beam_mountain_centroid(
                    sequence,
                    expected_heights,
                    beam_size=args.beam_size,
                )
                solver_seconds = time.perf_counter() - started

                total_seconds = time.perf_counter() - total_started
                print(
                    f"{length},{instance},{repeat},{sequence_seed},"
                    f"{args.beam_size},{bpp_seconds:.9f},"
                    f"{profile_seconds:.9f},{solver_seconds:.9f},"
                    f"{total_seconds:.9f},{len(result.pairs)}",
                    flush=True,
                )
                del bpp
                gc.collect()


if __name__ == "__main__":
    main()
