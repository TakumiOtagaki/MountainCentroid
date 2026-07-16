"""Python wrapper for the standalone C++ exact constrained solver."""

from __future__ import annotations

import math
from pathlib import Path
import subprocess
from typing import Sequence

from .exact import ExactDiagnostics, ExactResult
from .formatting import pairs_from_bracket
from .sequence import normalise_sequence


def default_cpp_exact_path() -> Path:
    """Return the exact solver binary built by ``make exact``."""
    return Path(__file__).resolve().parents[2] / "bin" / "exact_mountain_centroid"


def cpp_exact_mountain_centroid(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    executable: str | Path | None = None,
    timeout: float | None = None,
) -> ExactResult:
    """Run the C++ exact solver and return the shared exact-result type."""
    sequence = normalise_sequence(sequence)
    mu = tuple(float(value) for value in expected_heights)
    if len(mu) != len(sequence) - 1:
        raise ValueError(
            f"Expected {len(sequence) - 1} mountain heights for a "
            f"length-{len(sequence)} sequence, got {len(mu)}"
        )
    if any(not math.isfinite(value) for value in mu):
        raise ValueError("Expected mountain heights must be finite")

    binary = Path(executable) if executable is not None else default_cpp_exact_path()
    if not binary.is_file():
        raise FileNotFoundError(
            f"C++ exact solver not found at {binary}; run `make exact`"
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
        raise RuntimeError(f"C++ exact solver failed: {detail}")

    lines = completed.stdout.splitlines()
    if len(lines) != 3:
        raise RuntimeError("C++ exact solver returned malformed output")
    structure = lines[0]
    if len(structure) != len(sequence):
        raise RuntimeError("C++ exact solver returned a structure of wrong length")
    try:
        squared_error = float(lines[1])
        states, transitions, depth_levels = (int(value) for value in lines[2].split())
    except ValueError as error:
        raise RuntimeError("C++ exact solver returned malformed diagnostics") from error
    if depth_levels < 1:
        raise RuntimeError("C++ exact solver returned an invalid depth count")

    pairs = tuple(pairs_from_bracket(structure))
    heights = [0]
    depth = 0
    for character in structure:
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
        heights.append(depth)

    return ExactResult(
        structure=structure,
        pairs=pairs,
        heights=tuple(heights),
        squared_error=squared_error,
        diagnostics=ExactDiagnostics(
            states_evaluated=states,
            partner_transitions_evaluated=transitions,
            maximum_external_depth=depth_levels - 1,
        ),
    )
