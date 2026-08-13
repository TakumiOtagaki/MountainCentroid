"""Python wrapper for the C++ hybrid-objective solver."""

from __future__ import annotations

from dataclasses import dataclass
import math
from pathlib import Path
import subprocess
from typing import Sequence

import numpy as np

from .constrained import ConstrainedDiagnostics
from .formatting import pairs_from_bracket
from .sequence import normalise_sequence


@dataclass(frozen=True, slots=True)
class HybridResult:
    """One hybrid-objective result for a specified weight."""

    alpha: float
    structure: str
    pairs: tuple[tuple[int, int], ...]
    heights: tuple[int, ...]
    squared_mountain_error: float
    centroid_gain: float
    hybrid_objective: float
    solver_seconds: float
    diagnostics: ConstrainedDiagnostics


def default_cpp_hybrid_path() -> Path:
    """Return the solver binary built by ``make hybrid``."""
    return Path(__file__).resolve().parents[2] / "bin" / "hybrid_mountain_centroid"


def cpp_hybrid_mountain_centroid(
    sequence: str,
    expected_heights: Sequence[float],
    bpp: np.ndarray,
    alphas: Sequence[float],
    *,
    executable: str | Path | None = None,
    timeout: float | None = None,
) -> tuple[HybridResult, ...]:
    """Run the C++ solver for one or more hybrid weights."""
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

    probabilities = np.asarray(bpp, dtype=float)
    if probabilities.shape != (n, n):
        raise ValueError(f"bpp must have shape ({n}, {n})")
    upper = probabilities[np.triu_indices(n, 1)]
    if not np.all(np.isfinite(upper)):
        raise ValueError("BPP values must be finite")
    if np.any((upper < 0.0) | (upper > 1.0)):
        raise ValueError("BPP values must lie in [0,1]")

    weights = tuple(float(alpha) for alpha in alphas)
    if not weights:
        raise ValueError("At least one alpha value is required")
    if any(
        not math.isfinite(alpha) or alpha < 0.0 or alpha > 1.0
        for alpha in weights
    ):
        raise ValueError("Alpha values must lie in [0,1]")

    binary = Path(executable) if executable is not None else default_cpp_hybrid_path()
    if not binary.is_file():
        raise FileNotFoundError(
            f"C++ hybrid solver not found at {binary}; run `make hybrid`"
        )
    payload = "\n".join(
        (
            sequence,
            str(len(weights)),
            " ".join(repr(alpha) for alpha in weights),
            " ".join(repr(value) for value in mu),
            " ".join(repr(float(value)) for value in upper),
        )
    )
    completed = subprocess.run(
        [str(binary)],
        input=payload + "\n",
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or "no diagnostic output"
        raise RuntimeError(f"C++ hybrid solver failed: {detail}")

    lines = completed.stdout.splitlines()
    if len(lines) != len(weights):
        raise RuntimeError("C++ hybrid solver returned the wrong number of rows")
    mountain_scale = sum(min(cut, n - cut) ** 2 for cut in range(1, n))
    pair_scale = n // 2
    results = []
    for expected_alpha, line in zip(weights, lines):
        fields = line.split("\t")
        if len(fields) != 7:
            raise RuntimeError("C++ hybrid solver returned malformed output")
        try:
            returned_alpha = float(fields[0])
            scaled_objective = float(fields[2])
            states, transitions, depth_levels = map(int, fields[3:6])
            solver_seconds = float(fields[6])
        except ValueError as error:
            raise RuntimeError(
                "C++ hybrid solver returned malformed diagnostics"
            ) from error
        if not math.isclose(returned_alpha, expected_alpha, abs_tol=1e-15):
            raise RuntimeError("C++ hybrid solver returned alpha values out of order")
        structure = fields[1]
        if len(structure) != n:
            raise RuntimeError(
                "C++ hybrid solver returned a structure of wrong length"
            )
        if depth_levels < 1 or solver_seconds < 0.0:
            raise RuntimeError("C++ hybrid solver returned invalid diagnostics")

        pairs = tuple(pairs_from_bracket(structure))
        heights = [0]
        depth = 0
        for character in structure:
            if character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
            heights.append(depth)
        squared_error = sum(
            (height - expected) ** 2
            for height, expected in zip(heights[1:-1], mu)
        )
        centroid_gain = float(
            sum(
                2.0 * probabilities[left - 1, right - 1] - 1.0
                for left, right in pairs
            )
        )
        hybrid_objective = (
            (1.0 - expected_alpha) * squared_error / mountain_scale
            - expected_alpha * centroid_gain / pair_scale
        )
        expected_scaled = (
            -centroid_gain
            if expected_alpha == 1.0
            else squared_error
            - expected_alpha
            * mountain_scale
            / ((1.0 - expected_alpha) * pair_scale)
            * centroid_gain
        )
        if not math.isclose(
            scaled_objective,
            expected_scaled,
            rel_tol=1e-9,
            abs_tol=1e-9,
        ):
            raise RuntimeError("C++ hybrid solver returned an inconsistent objective")
        results.append(
            HybridResult(
                alpha=expected_alpha,
                structure=structure,
                pairs=pairs,
                heights=tuple(heights),
                squared_mountain_error=squared_error,
                centroid_gain=centroid_gain,
                hybrid_objective=hybrid_objective,
                solver_seconds=solver_seconds,
                diagnostics=ConstrainedDiagnostics(
                    states_evaluated=states,
                    partner_transitions_evaluated=transitions,
                    maximum_external_depth=depth_levels - 1,
                ),
            )
        )
    return tuple(results)
