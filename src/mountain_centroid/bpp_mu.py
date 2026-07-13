#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Base-pair probability backends and expected mountain heights."""
from __future__ import annotations

import math
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from typing import List, Literal, Tuple

import numpy as np

try:
    from .utils_tri import tri_to_full
    from .sequence import normalise_sequence
except Exception:
    from utils_tri import tri_to_full
    from sequence import normalise_sequence


Backend = Literal["vienna", "linearpartition"]


def mountain_expectation_from_bpp(bpp: np.ndarray) -> List[float]:
    """Return μ_k = Σ_{i<=k<j} p_ij for cuts k=1..n-1."""
    if bpp.ndim != 2 or bpp.shape[0] != bpp.shape[1]:
        raise ValueError("bpp must be a square matrix")

    n = int(bpp.shape[0])
    diff = [0.0] * (n + 1)
    for i in range(1, n):
        for j in range(i + 1, n + 1):
            pij = float(bpp[i - 1, j - 1])
            diff[i] += pij
            diff[j] -= pij

    mu = [0.0] * max(n - 1, 0)
    acc = 0.0
    for k in range(1, n):
        acc += diff[k]
        mu[k - 1] = acc
    return mu


def compute_bpp_vienna(seq: str, temperature: float = 37.0) -> np.ndarray:
    """Compute an exact thermodynamic BPP matrix with ViennaRNA/RNAlib."""
    import RNA  # ViennaRNA の Python バインディング

    seq = normalise_sequence(seq)
    n = len(seq)
    if n < 2:
        return np.zeros((n, n), dtype=float)

    md = RNA.md()
    md.temperature = float(temperature)

    fc = RNA.fold_compound(seq, md)
    fc.pf()
    bpp_tri = fc.exp_matrices.probs
    return tri_to_full(n, bpp_tri, fc.iindx)


def _find_linearpartition_runner(executable: str | os.PathLike[str] | None) -> Path:
    candidates: list[Path] = []
    if executable is not None:
        candidates.append(Path(executable).expanduser())

    env_path = os.environ.get("MOUNTAIN_CENTROID_LINEARPARTITION")
    if env_path:
        candidates.append(Path(env_path).expanduser())

    # Editable/source checkout: <repo>/src/mountain_centroid/bpp_mu.py
    candidates.append(Path(__file__).resolve().parents[2] / "vendor" / "LinearPartition" / "linearpartition")

    on_path = shutil.which("linearpartition")
    if on_path:
        candidates.append(Path(on_path))

    for candidate in candidates:
        if candidate.is_file():
            return candidate.resolve()

    raise FileNotFoundError(
        "LinearPartition runner not found. Clone submodules recursively and run "
        "`make -C vendor/LinearPartition`, or pass --linearpartition-path."
    )


def compute_bpp_linearpartition(
    seq: str,
    *,
    beam_size: int = 100,
    cutoff: float = 0.0,
    executable: str | os.PathLike[str] | None = None,
) -> np.ndarray:
    """Compute an approximate BPP matrix with LinearPartition-V.

    LinearPartition writes probabilities with five decimal places. Missing
    entries (including pairs pruned by the beam or cutoff) are treated as zero.
    """
    seq = normalise_sequence(seq)
    n = len(seq)
    if n < 2:
        return np.zeros((n, n), dtype=float)
    if beam_size < 0:
        raise ValueError("beam_size must be non-negative")
    if not 0.0 <= cutoff <= 1.0:
        raise ValueError("cutoff must be between 0 and 1")

    runner = _find_linearpartition_runner(executable)
    with tempfile.TemporaryDirectory(prefix="mountain-centroid-lp-") as tmpdir:
        output_path = Path(tmpdir) / "bpp.txt"
        cmd = [
            str(runner),
            "-V",
            "--beamsize",
            str(beam_size),
            "--rewrite",
            str(output_path),
            "--cutoff",
            str(cutoff),
        ]
        completed = subprocess.run(
            cmd,
            input=f"{seq}\n",
            text=True,
            capture_output=True,
            cwd=runner.parent,
            check=False,
        )
        if completed.returncode != 0 or not output_path.is_file():
            details = (completed.stderr or completed.stdout).strip()
            raise RuntimeError(f"LinearPartition failed (exit {completed.returncode}): {details}")

        bpp = np.zeros((n, n), dtype=float)
        for line_number, raw_line in enumerate(output_path.read_text().splitlines(), start=1):
            line = raw_line.strip()
            if not line:
                continue
            fields = line.split()
            if len(fields) != 3:
                raise ValueError(f"Invalid LinearPartition BPP line {line_number}: {raw_line!r}")
            i, j = int(fields[0]), int(fields[1])
            probability = float(fields[2])
            if not (1 <= i < j <= n and 0.0 <= probability <= 1.0):
                raise ValueError(f"Invalid LinearPartition BPP value on line {line_number}: {raw_line!r}")
            bpp[i - 1, j - 1] = probability

    return bpp


def compute_bpp_and_mu(
    seq: str,
    temperature: float = 37.0,
    *,
    backend: Backend = "vienna",
    beam_size: int = 100,
    cutoff: float = 0.0,
    linearpartition_path: str | os.PathLike[str] | None = None,
) -> Tuple[np.ndarray, List[float]]:
    """Compute a BPP matrix and its cut-based expected mountain height."""
    if backend == "vienna":
        bpp = compute_bpp_vienna(seq, temperature=temperature)
    elif backend == "linearpartition":
        if not math.isclose(float(temperature), 37.0, abs_tol=1e-9):
            raise ValueError("LinearPartition-V backend currently supports only its fixed 37°C model")
        bpp = compute_bpp_linearpartition(
            seq,
            beam_size=beam_size,
            cutoff=cutoff,
            executable=linearpartition_path,
        )
    else:
        raise ValueError(f"Unknown BPP backend: {backend}")

    return bpp, mountain_expectation_from_bpp(bpp)
