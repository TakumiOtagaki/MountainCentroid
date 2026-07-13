import random

import pytest

from mountain_centroid import predict_from_profile
from mountain_centroid.beam import beam_mountain_centroid
from mountain_centroid.formatting import pairs_from_bracket
from mountain_centroid.sequence import MIN_HAIRPIN_LENGTH, can_pair


def _exact_squared_error(sequence, mu):
    best = (float("inf"), None)
    n = len(sequence)

    def visit(index, stack, heights, structure):
        nonlocal best
        if index == n:
            if stack:
                return
            error = sum(
                (height - expected) ** 2
                for height, expected in zip(heights[:-1], mu)
            )
            candidate = (error, "".join(structure))
            if candidate < best:
                best = candidate
            return

        visit(
            index + 1,
            stack,
            heights + [len(stack)],
            structure + ["."],
        )
        if index + MIN_HAIRPIN_LENGTH + 1 < n:
            visit(
                index + 1,
                stack + (index,),
                heights + [len(stack) + 1],
                structure + ["("],
            )
        if stack:
            left = stack[-1]
            if (
                index - left - 1 >= MIN_HAIRPIN_LENGTH
                and can_pair(sequence[left], sequence[index])
            ):
                visit(
                    index + 1,
                    stack[:-1],
                    heights + [len(stack) - 1],
                    structure + [")"],
                )

    visit(0, (), [], [])
    return best


def test_recovers_nested_sequence_valid_structure_and_objective():
    sequence = "GGGAAACCC"
    mu = [1.0, 2.0, 3.0, 3.0, 3.0, 3.0, 2.0, 1.0]

    result = beam_mountain_centroid(sequence, mu, beam_size=100)

    assert result.structure == "(((...)))"
    assert result.pairs == ((1, 9), (2, 8), (3, 7))
    assert result.heights == (0, 1, 2, 3, 3, 3, 3, 2, 1, 0)
    assert result.squared_error == 0.0


def test_enforces_pairing_and_hairpin_constraints():
    target = [1.0, 1.0, 1.0, 1.0]

    assert beam_mountain_centroid("GAAAU", target).structure == "(...)"
    assert beam_mountain_centroid("GAAAG", target).structure == "....."
    assert beam_mountain_centroid("GACC", [1.0, 1.0, 1.0]).structure == "...."


def test_large_beam_matches_exhaustive_oracle_on_small_random_instances():
    rng = random.Random(7)
    for _ in range(30):
        n = rng.randint(5, 9)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        mu = [rng.random() * 2.5 for _ in range(n - 1)]

        exact_error, _ = _exact_squared_error(sequence, mu)
        result = beam_mountain_centroid(sequence, mu, beam_size=100_000)

        assert result.squared_error == pytest.approx(exact_error)


def test_public_profile_api_normalises_sequence_and_reports_direct_error():
    mu = [1.0, 1.0, 1.0, 1.0]
    prediction = predict_from_profile("gaaAt", mu, beam_size=50)

    direct_error = sum(
        (height - expected) ** 2
        for height, expected in zip(
            prediction.mountain_heights[1:-1],
            prediction.expected_mountain_heights,
        )
    )
    assert prediction.sequence == "GAAAU"
    assert prediction.structure == "(...)"
    assert prediction.bpp_backend == "provided"
    assert prediction.squared_mountain_error == direct_error


def test_random_outputs_always_satisfy_structural_invariants():
    rng = random.Random(11)
    for _ in range(50):
        n = rng.randint(5, 60)
        sequence = "".join(rng.choice("ACGU") for _ in range(n))
        mu = [rng.random() * 6 for _ in range(n - 1)]

        result = beam_mountain_centroid(sequence, mu, beam_size=50)

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
        assert all(
            abs(right - left) <= 1
            for left, right in zip(result.heights, result.heights[1:])
        )


@pytest.mark.parametrize("beam_size", [0, -1, True, 1.5])
def test_rejects_invalid_beam_size(beam_size):
    with pytest.raises(ValueError, match="positive integer"):
        beam_mountain_centroid("AAAAA", [0.0] * 4, beam_size=beam_size)
