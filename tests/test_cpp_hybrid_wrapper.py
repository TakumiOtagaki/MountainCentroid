from pathlib import Path
import subprocess

import numpy as np
import pytest

from mountain_centroid.bpp_mu import mountain_expectation_from_bpp
from mountain_centroid.cpp_hybrid import cpp_hybrid_mountain_centroid


@pytest.fixture(scope="session")
def cpp_hybrid_binary() -> Path:
    repository = Path(__file__).resolve().parents[1]
    subprocess.run(
        ["make", "hybrid"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    return repository / "bin" / "hybrid_mountain_centroid"


def test_cpp_hybrid_wrapper_returns_defined_objective(cpp_hybrid_binary):
    sequence = "GGGAAACCC"
    bpp = np.zeros((len(sequence), len(sequence)), dtype=float)
    bpp[0, 8] = 0.9
    bpp[1, 7] = 0.8
    bpp[2, 6] = 0.7
    mu = mountain_expectation_from_bpp(bpp)
    alphas = (0.0, 0.23076923076923078, 1.0)

    results = cpp_hybrid_mountain_centroid(
        sequence,
        mu,
        bpp,
        alphas,
        executable=cpp_hybrid_binary,
    )

    assert tuple(result.alpha for result in results) == alphas
    mountain_scale = sum(
        min(cut, len(sequence) - cut) ** 2
        for cut in range(1, len(sequence))
    )
    for result in results:
        expected = (
            (1.0 - result.alpha)
            * result.squared_mountain_error
            / mountain_scale
            - result.alpha * result.centroid_gain / 4.0
        )
        assert result.hybrid_objective == pytest.approx(expected)
        assert result.heights[0] == result.heights[-1] == 0
        assert result.solver_seconds >= 0.0


@pytest.mark.parametrize(
    ("alphas", "bpp_value", "message"),
    [
        ((), 0.0, "At least one"),
        ((-0.1,), 0.0, "Alpha values"),
        ((0.5,), 1.1, "BPP values"),
    ],
)
def test_cpp_hybrid_wrapper_rejects_invalid_values(
    cpp_hybrid_binary,
    alphas,
    bpp_value,
    message,
):
    bpp = np.zeros((5, 5), dtype=float)
    bpp[0, 4] = bpp_value
    with pytest.raises(ValueError, match=message):
        cpp_hybrid_mountain_centroid(
            "GAAAU",
            [0.0] * 4,
            bpp,
            alphas,
            executable=cpp_hybrid_binary,
        )


def test_cpp_hybrid_wrapper_rejects_missing_binary(tmp_path):
    with pytest.raises(FileNotFoundError, match="make hybrid"):
        cpp_hybrid_mountain_centroid(
            "GAAAU",
            [0.0] * 4,
            np.zeros((5, 5), dtype=float),
            (0.5,),
            executable=Path(tmp_path) / "missing",
        )
