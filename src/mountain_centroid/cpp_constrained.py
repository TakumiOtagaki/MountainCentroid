"""Python wrapper for the C++ sequence-constrained solver."""

from __future__ import annotations

import math
from pathlib import Path
import subprocess
from typing import Sequence

from .constrained import ConstrainedDiagnostics, ConstrainedResult
from .formatting import pairs_from_bracket
from .sequence import normalise_sequence


def default_cpp_constrained_path() -> Path:
    """Return the solver binary built by ``make constrained``."""
    return (
        Path(__file__).resolve().parents[2]
        / "bin"
        / "sequence_constrained_mountain_centroid"
    )


def cpp_sequence_constrained_mountain_centroid(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    executable: str | Path | None = None,
    timeout: float | None = None,
) -> ConstrainedResult:
    """Run the C++ solver and return the shared constrained-result type."""
    sequence = normalise_sequence(sequence)
    mu = tuple(float(value) for value in expected_heights)
    if len(mu) != len(sequence) - 1:
        raise ValueError(
            f"Expected {len(sequence) - 1} mountain heights for a "
            f"length-{len(sequence)} sequence, got {len(mu)}"
        )
    if any(not math.isfinite(value) for value in mu):
        raise ValueError("Expected mountain heights must be finite")

    binary = (
        Path(executable)
        if executable is not None
        else default_cpp_constrained_path()
    )
    if not binary.is_file():
        raise FileNotFoundError(
            f"C++ constrained solver not found at {binary}; "
            "run `make constrained`"
        )
    payload = sequence + "\n" + " ".join(repr(value) for value in mu) + "\n"
    completed = subprocess.run(
        [str(binary)],
        input=payload,
        text=True,
        capture_output=True,
        check=False,
        timeout=timeout,
    )
    if completed.returncode != 0:
        detail = completed.stderr.strip() or "no diagnostic output"
        raise RuntimeError(f"C++ constrained solver failed: {detail}")

    lines = completed.stdout.splitlines()
    if len(lines) != 3:
        raise RuntimeError("C++ constrained solver returned malformed output")
    structure = lines[0]
    if len(structure) != len(sequence):
        raise RuntimeError(
            "C++ constrained solver returned a structure of wrong length"
        )
    try:
        squared_error = float(lines[1])
        states, transitions, depth_levels = (int(value) for value in lines[2].split())
    except ValueError as error:
        raise RuntimeError(
            "C++ constrained solver returned malformed diagnostics"
        ) from error
    if depth_levels < 1:
        raise RuntimeError("C++ constrained solver returned an invalid depth count")

    pairs = tuple(pairs_from_bracket(structure))
    heights = [0]
    depth = 0
    for character in structure:
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        heights.append(depth)

    return ConstrainedResult(
        structure=structure,
        pairs=pairs,
        heights=tuple(heights),
        squared_error=squared_error,
        diagnostics=ConstrainedDiagnostics(
            states_evaluated=states,
            partner_transitions_evaluated=transitions,
            maximum_external_depth=depth_levels - 1,
        ),
    )
