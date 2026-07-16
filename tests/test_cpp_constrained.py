import random
import subprocess
from pathlib import Path

import pytest

from mountain_centroid.constrained import sequence_constrained_mountain_centroid
from mountain_centroid.cpp_constrained import (
    cpp_sequence_constrained_mountain_centroid,
)


@pytest.fixture(scope="session")
def cpp_constrained_binary() -> Path:
    repository = Path(__file__).resolve().parents[1]
    subprocess.run(
        ["make", "constrained"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    return repository / "bin" / "sequence_constrained_mountain_centroid"


def test_cpp_matches_python_constrained_on_random_instances(
    cpp_constrained_binary,
):
    rng = random.Random(41)
    for _ in range(80):
        n = rng.randint(5, 30)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        mu = [rng.random() * min(cut, n - cut, 6) for cut in range(1, n)]

        python_result = sequence_constrained_mountain_centroid(sequence, mu)
        cpp_result = cpp_sequence_constrained_mountain_centroid(
            sequence,
            mu,
            executable=cpp_constrained_binary,
        )

        assert cpp_result.structure == python_result.structure
        assert cpp_result.pairs == python_result.pairs
        assert cpp_result.heights == python_result.heights
        assert cpp_result.squared_error == pytest.approx(
            python_result.squared_error,
            rel=1e-13,
            abs=1e-13,
        )
        assert cpp_result.diagnostics == python_result.diagnostics


def test_cpp_normalizes_dna_input(cpp_constrained_binary):
    result = cpp_sequence_constrained_mountain_centroid(
        "gaaAt",
        [1.0, 1.0, 1.0, 1.0],
        executable=cpp_constrained_binary,
    )
    assert result.structure == "(...)"


def test_cpp_wrapper_rejects_missing_binary(tmp_path):
    with pytest.raises(FileNotFoundError, match="make constrained"):
        cpp_sequence_constrained_mountain_centroid(
            "GAAAU",
            [1.0] * 4,
            executable=tmp_path / "missing",
        )
