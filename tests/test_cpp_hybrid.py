import random
import subprocess
from pathlib import Path

import numpy as np
import pytest

from mountain_centroid.bpp_mu import mountain_expectation_from_bpp
from mountain_centroid.cpp_constrained import (
    cpp_sequence_constrained_mountain_centroid,
)
from mountain_centroid.formatting import pairs_from_bracket
from mountain_centroid.sequence import MIN_HAIRPIN_LENGTH, can_pair


@pytest.fixture(scope="session")
def cpp_hybrid_binaries() -> tuple[Path, Path]:
    repository = Path(__file__).resolve().parents[1]
    subprocess.run(
        ["make", "constrained", "hybrid"],
        cwd=repository,
        check=True,
        capture_output=True,
        text=True,
    )
    return (
        repository / "bin" / "sequence_constrained_mountain_centroid",
        repository / "bin" / "hybrid_mountain_centroid",
    )


def _run_hybrid(binary, sequence, mu, bpp, alphas):
    upper_triangle = (
        str(float(bpp[left, right]))
        for left in range(len(sequence))
        for right in range(left + 1, len(sequence))
    )
    payload = "\n".join(
        (
            sequence,
            str(len(alphas)),
            " ".join(str(alpha) for alpha in alphas),
            " ".join(str(value) for value in mu),
            " ".join(upper_triangle),
        )
    )
    completed = subprocess.run(
        [str(binary)],
        input=payload + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(completed.stderr.strip())
    rows = []
    for line in completed.stdout.splitlines():
        fields = line.split("\t")
        rows.append(
            {
                "alpha": float(fields[0]),
                "structure": fields[1],
                "objective": float(fields[2]),
            }
        )
    return rows


def _enumerate_structures(sequence):
    structures = []
    n = len(sequence)

    def visit(index, stack, pairs, heights):
        if index == n:
            if not stack:
                structures.append((tuple(pairs), tuple(heights)))
            return
        visit(index + 1, stack, pairs, heights + [len(stack)])
        visit(index + 1, stack + (index,), pairs, heights + [len(stack) + 1])
        if stack:
            left = stack[-1]
            if (
                index - left - 1 >= MIN_HAIRPIN_LENGTH
                and can_pair(sequence[left], sequence[index])
            ):
                visit(
                    index + 1,
                    stack[:-1],
                    pairs + [(left, index)],
                    heights + [len(stack) - 1],
                )

    visit(0, (), [], [])
    return structures


def _scaled_hybrid_objective(pairs, heights, mu, bpp, alpha):
    mountain_loss = sum(
        (height - expected) ** 2
        for height, expected in zip(heights[:-1], mu)
    )
    centroid_gain = sum(2.0 * bpp[left, right] - 1.0 for left, right in pairs)
    if alpha == 1.0:
        return -centroid_gain
    n = len(heights)
    mountain_scale = sum(min(cut, n - cut) ** 2 for cut in range(1, n))
    pair_scale = n // 2
    pair_weight = alpha * mountain_scale / ((1.0 - alpha) * pair_scale)
    return mountain_loss - pair_weight * centroid_gain


def test_cpp_hybrid_matches_exhaustive_objective(cpp_hybrid_binaries):
    _, hybrid_binary = cpp_hybrid_binaries
    rng = random.Random(53)
    alphas = (0.0, 0.001, 0.23076923076923078, 1.0)
    for _ in range(30):
        n = rng.randint(5, 10)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        bpp = np.zeros((n, n), dtype=float)
        for left in range(n):
            for right in range(left + 1, n):
                bpp[left, right] = rng.random()
        mu = mountain_expectation_from_bpp(bpp)
        structures = _enumerate_structures(sequence)
        rows = _run_hybrid(hybrid_binary, sequence, mu, bpp, alphas)

        for row, alpha in zip(rows, alphas):
            optimum = min(
                _scaled_hybrid_objective(pairs, heights, mu, bpp, alpha)
                for pairs, heights in structures
            )
            assert row["alpha"] == pytest.approx(alpha)
            assert row["objective"] == pytest.approx(optimum, abs=1e-11)
            pairs = pairs_from_bracket(row["structure"])
            assert all(
                can_pair(sequence[left - 1], sequence[right - 1])
                and right - left - 1 >= MIN_HAIRPIN_LENGTH
                for left, right in pairs
            )


def test_alpha_zero_matches_existing_cpp_solver(cpp_hybrid_binaries):
    constrained_binary, hybrid_binary = cpp_hybrid_binaries
    rng = random.Random(59)
    for _ in range(30):
        n = rng.randint(5, 30)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        bpp = np.zeros((n, n), dtype=float)
        for left in range(n):
            for right in range(left + 1, n):
                bpp[left, right] = rng.random()
        mu = mountain_expectation_from_bpp(bpp)

        constrained = cpp_sequence_constrained_mountain_centroid(
            sequence,
            mu,
            executable=constrained_binary,
        )
        hybrid = _run_hybrid(hybrid_binary, sequence, mu, bpp, (0.0,))[0]
        assert hybrid["structure"] == constrained.structure
        assert hybrid["objective"] == pytest.approx(constrained.squared_error)


@pytest.mark.parametrize(
    ("alpha", "probability", "message"),
    [
        (-0.1, 0.0, "Alpha must lie in"),
        (1.1, 0.0, "Alpha must lie in"),
        (0.5, 1.1, "invalid BPP value"),
    ],
)
def test_cpp_hybrid_rejects_invalid_input(
    cpp_hybrid_binaries,
    alpha,
    probability,
    message,
):
    _, hybrid_binary = cpp_hybrid_binaries
    sequence = "GAAAU"
    payload = f"{sequence}\n1\n{alpha}\n0 0 0 0\n" + " ".join(
        str(probability) for _ in range(10)
    )
    completed = subprocess.run(
        [str(hybrid_binary)],
        input=payload + "\n",
        text=True,
        capture_output=True,
        check=False,
    )
    assert completed.returncode == 2
    assert message in completed.stderr
