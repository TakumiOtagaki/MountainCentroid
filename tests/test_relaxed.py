import itertools
import random

import pytest

from mountain_centroid.relaxed import relaxed_mountain_centroid
from mountain_centroid.constrained import sequence_constrained_mountain_centroid


def _exhaustive(mu):
    n = len(mu) + 1
    best = (float("inf"), None)
    for steps in itertools.product((-1, 0, 1), repeat=n):
        heights = [0]
        for step in steps:
            heights.append(heights[-1] + step)
        if heights[-1] != 0 or min(heights) < 0:
            continue
        error = sum((height - expected) ** 2 for height, expected in zip(heights[1:-1], mu))
        best = min(best, (error, tuple(heights)))
    return best[0]


def test_relaxed_recovers_nested_path():
    result = relaxed_mountain_centroid([1, 2, 3, 3, 3, 3, 2, 1])
    assert result.structure == "(((...)))"
    assert result.heights == (0, 1, 2, 3, 3, 3, 3, 2, 1, 0)
    assert result.squared_error == 0.0


def test_relaxed_matches_exhaustive_random_profiles():
    random_generator = random.Random(19)
    for n in range(2, 8):
        for _ in range(10):
            mu = [random_generator.random() * 3 for _ in range(n - 1)]
            result = relaxed_mountain_centroid(mu)
            assert result.squared_error == pytest.approx(_exhaustive(mu))


def test_relaxed_rejects_nonfinite_profile():
    with pytest.raises(ValueError, match="finite"):
        relaxed_mountain_centroid([float("nan")])


def test_relaxed_objective_is_a_lower_bound_for_constrained_prediction():
    sequence = "GGGAAACCC"
    mu = [0.4, 1.3, 2.2, 2.8, 2.7, 2.0, 1.1, 0.3]
    relaxed = relaxed_mountain_centroid(mu)
    constrained = sequence_constrained_mountain_centroid(sequence, mu)
    assert relaxed.squared_error <= constrained.squared_error
