"""Small public API for Mountain Centroid prediction."""

from __future__ import annotations

from dataclasses import dataclass
import os
from typing import Sequence

from .beam import beam_mountain_centroid
from .bpp_mu import Backend, compute_bpp_and_mu
from .sequence import normalise_sequence


@dataclass(frozen=True, slots=True)
class Prediction:
    sequence: str
    structure: str
    pairs: tuple[tuple[int, int], ...]
    mountain_heights: tuple[int, ...]
    expected_mountain_heights: tuple[float, ...]
    squared_mountain_error: float
    solver_beam_size: int
    bpp_backend: str


def predict_from_profile(
    sequence: str,
    expected_heights: Sequence[float],
    *,
    beam_size: int = 100,
    bpp_backend: str = "provided",
) -> Prediction:
    """Predict one valid structure from an already computed mean profile."""
    sequence = normalise_sequence(sequence)
    result = beam_mountain_centroid(
        sequence,
        expected_heights,
        beam_size=beam_size,
    )
    return Prediction(
        sequence=sequence,
        structure=result.structure,
        pairs=result.pairs,
        mountain_heights=result.heights,
        expected_mountain_heights=tuple(float(x) for x in expected_heights),
        squared_mountain_error=result.squared_error,
        solver_beam_size=beam_size,
        bpp_backend=bpp_backend,
    )


def predict(
    sequence: str,
    *,
    temperature: float = 37.0,
    bpp_backend: Backend = "linearpartition",
    beam_size: int = 100,
    bpp_beam_size: int = 100,
    bpp_cutoff: float = 0.0,
    linearpartition_path: str | os.PathLike[str] | None = None,
) -> Prediction:
    """Compute BPPs and return the beam-pruned Mountain Centroid prediction."""
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
        beam_size=beam_size,
        bpp_backend=bpp_backend,
    )
