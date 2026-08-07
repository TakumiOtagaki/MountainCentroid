#!/usr/bin/env python3
"""Command-line interface for the single Mountain Centroid estimator."""

from __future__ import annotations

import argparse

from .api import predict


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Sequence-constrained, pseudoknot-free Mountain Centroid "
            "prediction"
        )
    )
    parser.add_argument("--seq", required=True, help="RNA sequence (A/C/G/U)")
    parser.add_argument(
        "--temp",
        type=float,
        default=37.0,
        help=(
            "temperature in °C "
            "(ViennaRNA only; LinearPartition-V is fixed at 37°C)"
        ),
    )
    parser.add_argument(
        "--bpp-backend",
        choices=("vienna", "linearpartition"),
        default="vienna",
        help="BPP backend: ViennaRNA (default) or LinearPartition",
    )
    parser.add_argument(
        "--bpp-beam-size",
        type=int,
        default=100,
        help="LinearPartition beam size [default 100]",
    )
    parser.add_argument(
        "--bpp-cutoff",
        type=float,
        default=0.0,
        help="LinearPartition BPP output cutoff",
    )
    parser.add_argument(
        "--linearpartition-path",
        default=None,
        help="path to the LinearPartition runner script",
    )
    parser.add_argument(
        "--constrained-executable",
        default=None,
        help="path to the compiled sequence-constrained solver",
    )
    args = parser.parse_args()

    prediction = predict(
        args.seq,
        temperature=args.temp,
        bpp_backend=args.bpp_backend,
        bpp_beam_size=args.bpp_beam_size,
        bpp_cutoff=args.bpp_cutoff,
        linearpartition_path=args.linearpartition_path,
        constrained_executable=args.constrained_executable,
    )

    preview = ", ".join(
        f"{value:.3f}" for value in prediction.expected_mountain_heights[:10]
    )
    if len(prediction.expected_mountain_heights) > 10:
        preview += " ..."

    print(f"# sequence (n={len(prediction.sequence)}): {prediction.sequence}")
    print(f"# BPP backend: {prediction.bpp_backend}")
    print("# solver: sequence-constrained interval-depth DP (C++)")
    print(f"# E[h] preview: [{preview}]")
    print(
        "squared_mountain_error =",
        f"{prediction.squared_mountain_error:.6f}",
    )
    print("dot_bracket =", prediction.structure)


if __name__ == "__main__":
    main()
