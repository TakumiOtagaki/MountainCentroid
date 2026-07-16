import random

import pytest

from mountain_centroid.constrained import sequence_constrained_mountain_centroid
from mountain_centroid.formatting import pairs_from_bracket
from mountain_centroid.relaxed import relaxed_mountain_centroid
from mountain_centroid.sequence import MIN_HAIRPIN_LENGTH, can_pair

def _exhaustive_squared_error(sequence, mu):
    best = float("inf")
    n = len(sequence)

    def visit(index, stack, heights):
        nonlocal best
        if index == n:
            if stack:
                return
            error = sum(
                (height - expected) ** 2
                for height, expected in zip(heights[:-1], mu)
            )
            best = min(best, error)
            return
        visit(index + 1, stack, heights + [len(stack)])
        if index + MIN_HAIRPIN_LENGTH + 1 < n:
            visit(index + 1, stack + (index,), heights + [len(stack) + 1])
        if stack:
            left = stack[-1]
            if (
                index - left - 1 >= MIN_HAIRPIN_LENGTH
                and can_pair(sequence[left], sequence[index])
            ):
                visit(index + 1, stack[:-1], heights + [len(stack) - 1])

    visit(0, (), [])
    return best


def test_constrained_recovers_nested_sequence_valid_structure_and_objective():
    sequence = "GGGAAACCC"
    mu = [1.0, 2.0, 3.0, 3.0, 3.0, 3.0, 2.0, 1.0]

    result = sequence_constrained_mountain_centroid(sequence, mu)

    assert result.structure == "(((...)))"
    assert result.pairs == ((1, 9), (2, 8), (3, 7))
    assert result.heights == (0, 1, 2, 3, 3, 3, 3, 2, 1, 0)
    assert result.squared_error == 0.0
    assert result.diagnostics.maximum_external_depth == 3


def test_constrained_enforces_pairing_and_hairpin_constraints():
    target = [1.0, 1.0, 1.0, 1.0]

    assert sequence_constrained_mountain_centroid("GAAAU", target).structure == "(...)"
    assert sequence_constrained_mountain_centroid("GAAAG", target).structure == "....."
    assert sequence_constrained_mountain_centroid("GACC", [1.0, 1.0, 1.0]).structure == "...."


def test_constrained_matches_exhaustive_oracle_on_small_random_instances():
    rng = random.Random(29)
    for _ in range(80):
        n = rng.randint(5, 10)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        mu = [rng.random() * min(cut, n - cut, 3) for cut in range(1, n)]

        optimum = _exhaustive_squared_error(sequence, mu)
        result = sequence_constrained_mountain_centroid(sequence, mu)

        assert result.squared_error == pytest.approx(optimum)


def test_constrained_random_outputs_satisfy_structural_invariants():
    rng = random.Random(31)
    for _ in range(40):
        n = rng.randint(5, 35)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        mu = [rng.random() * min(cut, n - cut, 6) for cut in range(1, n)]

        result = sequence_constrained_mountain_centroid(sequence, mu)

        assert pairs_from_bracket(result.structure) == list(result.pairs)
        for left, right in result.pairs:
            assert can_pair(sequence[left - 1], sequence[right - 1])
            assert right - left - 1 >= MIN_HAIRPIN_LENGTH
        assert all(
            not (left_a < left_b < right_a < right_b)
            for left_a, right_a in result.pairs
            for left_b, right_b in result.pairs
        )
        assert result.heights[0] == result.heights[-1] == 0
        assert min(result.heights) == 0
        direct_error = sum(
            (height - expected) ** 2
            for height, expected in zip(result.heights[1:-1], mu)
        )
        assert result.squared_error == pytest.approx(direct_error)


def test_relaxed_objective_bounds_constrained_objective():
    rng = random.Random(37)
    for _ in range(30):
        n = rng.randint(8, 20)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        mu = [rng.random() * min(cut, n - cut, 5) for cut in range(1, n)]

        relaxed = relaxed_mountain_centroid(mu)
        constrained = sequence_constrained_mountain_centroid(sequence, mu)

        assert relaxed.squared_error <= constrained.squared_error + 1e-12


@pytest.mark.parametrize(
    ("sequence", "mu", "message"),
    [
        ("GAAAU", [0.0], "Expected 4"),
        ("GAAAU", [0.0, 0.0, float("nan"), 0.0], "finite"),
        ("GAANU", [0.0] * 4, "unsupported"),
    ],
)
def test_constrained_rejects_invalid_inputs(sequence, mu, message):
    with pytest.raises(ValueError, match=message):
        sequence_constrained_mountain_centroid(sequence, mu)
