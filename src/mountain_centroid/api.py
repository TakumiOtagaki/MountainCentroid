"""Small public API for Mountain Centroid prediction."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Sequence

from .bpp_mu import Backend, compute_bpp_and_mu
from .cpp_constrained import cpp_sequence_constrained_mountain_centroid
from .cpp_hybrid import cpp_hybrid_mountain_centroid
from .sequence import normalise_sequence


@dataclass(frozen=True, slots=True)
class Prediction:
    sequence: str
    structure: str
    pairs: tuple[tuple[int, int], ...]
    mountain_heights: tuple[int, ...]
    expected_mountain_heights: tuple[float, ...]
    squared_mountain_error: float
    solver_backend: str
    bpp_backend: str


@dataclass(frozen=True, slots=True)
class HybridPrediction:
    sequence: str
    alpha: float
    structure: str
    pairs: tuple[tuple[int, int], ...]
    mountain_heights: tuple[int, ...]
    expected_mountain_heights: tuple[float, ...]
    squared_mountain_error: float
    centroid_gain: float
    hybrid_objective: float
    solver_backend: str
    bpp_backend: str


def predict_from_profile(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    constrained_executable: str | os.PathLike[str] | None = None,
    bpp_backend: str = "provided",
) -> Prediction:
    """Predict one valid structure from an already computed mean profile."""
    sequence = normalise_sequence(sequence)
    result = cpp_sequence_constrained_mountain_centroid(
        sequence,
        expected_heights,
        executable=constrained_executable,
    )
    return Prediction(
        sequence=sequence,
        structure=result.structure,
        pairs=result.pairs,
        mountain_heights=result.heights,
        expected_mountain_heights=tuple(float(x) for x in expected_heights),
        squared_mountain_error=result.squared_error,
        solver_backend="cpp_interval_depth_dp",
        bpp_backend=bpp_backend,
    )


def predict(
    sequence: str,
    *,
    temperature: float = 37.0,
    bpp_backend: Backend = "vienna",
    bpp_beam_size: int = 100,
    bpp_cutoff: float = 0.0,
    linearpartition_path: str | os.PathLike[str] | None = None,
    constrained_executable: str | os.PathLike[str] | None = None,
) -> Prediction:
    """Compute BPPs and return the sequence-constrained Mountain Centroid."""
    sequence = normalise_sequence(sequence)
    _, expected_heights = compute_bpp_and_mu(
        sequence,
        temperature=temperature,
        backend=bpp_backend,
        beam_size=bpp_beam_size,
        cutoff=bpp_cutoff,
        linearpartition_path=linearpartition_path,
    )
    return predict_from_profile(
        sequence,
        expected_heights,
        bpp_backend=bpp_backend,
        constrained_executable=constrained_executable,
    )


def predict_hybrid_curve(
    sequence: str,
    alphas: Sequence[float],
    *,
    temperature: float = 37.0,
    bpp_backend: Backend = "vienna",
    bpp_beam_size: int = 100,
    bpp_cutoff: float = 0.0,
    linearpartition_path: str | os.PathLike[str] | None = None,
    hybrid_executable: str | os.PathLike[str] | None = None,
) -> tuple[HybridPrediction, ...]:
    """Compute BPPs once and return predictions across hybrid weights."""
    sequence = normalise_sequence(sequence)
    bpp, expected_heights = compute_bpp_and_mu(
        sequence,
        temperature=temperature,
        backend=bpp_backend,
        beam_size=bpp_beam_size,
        cutoff=bpp_cutoff,
        linearpartition_path=linearpartition_path,
    )
    results = cpp_hybrid_mountain_centroid(
        sequence,
        expected_heights,
        bpp,
        alphas,
        executable=hybrid_executable,
    )
    expected_profile = tuple(float(value) for value in expected_heights)
    return tuple(
        HybridPrediction(
            sequence=sequence,
            alpha=result.alpha,
            structure=result.structure,
            pairs=result.pairs,
            mountain_heights=result.heights,
            expected_mountain_heights=expected_profile,
            squared_mountain_error=result.squared_mountain_error,
            centroid_gain=result.centroid_gain,
            hybrid_objective=result.hybrid_objective,
            solver_backend="cpp_hybrid_interval_depth_dp",
            bpp_backend=bpp_backend,
        )
        for result in results
    )


def predict_hybrid(
    sequence: str,
    alpha: float,
    *,
    temperature: float = 37.0,
    bpp_backend: Backend = "vienna",
    bpp_beam_size: int = 100,
    bpp_cutoff: float = 0.0,
    linearpartition_path: str | os.PathLike[str] | None = None,
    hybrid_executable: str | os.PathLike[str] | None = None,
) -> HybridPrediction:
    """Compute one prediction from the hybrid objective."""
    return predict_hybrid_curve(
        sequence,
        (alpha,),
        temperature=temperature,
        bpp_backend=bpp_backend,
        bpp_beam_size=bpp_beam_size,
        bpp_cutoff=bpp_cutoff,
        linearpartition_path=linearpartition_path,
        hybrid_executable=hybrid_executable,
    )[0]
