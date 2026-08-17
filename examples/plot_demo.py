#!/usr/bin/env python3
"""Generate a small demonstration from the public MountainCentroid API."""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Arc

from mountain_centroid import Prediction, predict


DEFAULT_SEQUENCE = "GGGAAAUUCCC"
PROFILE_COLOR = "#009E73"
EXPECTED_COLOR = "#D55E00"


def draw_structure(axis: plt.Axes, prediction: Prediction) -> None:
    """Draw the predicted base pairs above the nucleotide sequence."""
    n = len(prediction.sequence)
    maximum_span = max((right - left for left, right in prediction.pairs), default=1)
    for left, right in prediction.pairs:
        span = right - left
        axis.add_patch(
            Arc(
                ((left + right) / 2, 0),
                width=span,
                height=0.65 * span,
                theta1=0,
                theta2=180,
                color=PROFILE_COLOR,
                linewidth=2,
            )
        )
    axis.axhline(0, color="#303030", linewidth=1)
    for position, nucleotide in enumerate(prediction.sequence, start=1):
        axis.text(position, -0.35, nucleotide, ha="center", va="top", fontsize=10)
    axis.set_xlim(0.5, n + 0.5)
    axis.set_ylim(-0.8, 0.65 * maximum_span / 2 + 0.7)
    axis.set_title(
        f"Mountain Centroid prediction: {prediction.structure}",
        loc="left",
        fontweight="bold",
    )
    axis.set_axis_off()


def draw_profiles(axis: plt.Axes, prediction: Prediction) -> None:
    """Draw the ensemble mean and predicted mountain profiles."""
    n = len(prediction.sequence)
    boundaries = range(1, n)
    axis.plot(
        boundaries,
        prediction.expected_mountain_heights,
        color=EXPECTED_COLOR,
        marker="o",
        markerfacecolor="white",
        linewidth=2,
        label=r"Ensemble mean profile $\mu(k)$",
        zorder=3,
    )
    axis.step(
        range(n + 1),
        prediction.mountain_heights,
        where="post",
        color=PROFILE_COLOR,
        linewidth=2,
        label=r"Mountain Centroid profile $h(k)$",
        zorder=2,
    )
    axis.set_xlim(0, n)
    axis.set_ylim(bottom=0)
    axis.set_xlabel("Nucleotide position $k$")
    axis.set_ylabel("Nesting depth")
    axis.grid(axis="y", color="#E6E6E6", linewidth=0.8)
    axis.spines[["top", "right"]].set_visible(False)
    axis.legend(frameon=False, loc="upper center", ncols=2)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sequence", default=DEFAULT_SEQUENCE)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("docs/assets/mountain_centroid_demo.png"),
    )
    args = parser.parse_args()

    prediction = predict(args.sequence)
    figure, (structure_axis, profile_axis) = plt.subplots(
        2,
        1,
        figsize=(8.0, 5.2),
        gridspec_kw={"height_ratios": (1.05, 1.0)},
        constrained_layout=True,
    )
    draw_structure(structure_axis, prediction)
    draw_profiles(profile_axis, prediction)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(args.output, dpi=200, facecolor="white")
    plt.close(figure)
    print(f"sequence: {prediction.sequence}")
    print(f"structure: {prediction.structure}")
    print(f"squared mountain error: {prediction.squared_mountain_error:.6f}")
    print(f"wrote: {args.output}")


if __name__ == "__main__":
    main()
