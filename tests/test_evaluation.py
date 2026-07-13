import pytest

from mountain_centroid.evaluation import (
    base_pair_f1,
    mean_squared_mountain_distance,
    mountain_heights,
    squared_mountain_distance,
)


def test_base_pair_f1_uses_exact_pairs():
    assert base_pair_f1("((...))", "((...))") == 1.0
    assert base_pair_f1("(.....)", ".(...).") == 0.0


def test_squared_mountain_distance_is_explicit_and_normalized():
    predicted = "(.....)"
    reference = "......."

    assert mountain_heights(predicted) == (1, 1, 1, 1, 1, 1)
    assert squared_mountain_distance(predicted, reference) == 6.0
    assert mean_squared_mountain_distance(predicted, reference) == 1.0


def test_metrics_reject_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        base_pair_f1("...", "....")
    with pytest.raises(ValueError, match="same length"):
        squared_mountain_distance("...", "....")
