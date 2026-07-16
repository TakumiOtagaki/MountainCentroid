import pytest

from mountain_centroid.evaluation import (
    base_pair_f1,
    mean_squared_mountain_distance,
    mountain_heights,
    normalized_squared_mountain_distance,
    squared_mountain_distance,
)


def test_base_pair_f1_uses_pair_endpoints():
    assert base_pair_f1("((...))", "((...))") == 1.0
    assert base_pair_f1("(.....)", ".(...).") == 0.0


def test_squared_mountain_distance_is_explicit_and_normalized():
    predicted = "(.....)"
    reference = "......."

    assert mountain_heights(predicted) == (1, 1, 1, 1, 1, 1)
    assert squared_mountain_distance(predicted, reference) == 6.0
    assert mean_squared_mountain_distance(predicted, reference) == 1.0
    assert normalized_squared_mountain_distance(predicted, reference) == pytest.approx(
        6.0 / 28.0
    )


def test_max_normalized_mountain_distance_is_bounded():
    assert normalized_squared_mountain_distance("(((.)))", ".......") == 1.0
    assert normalized_squared_mountain_distance(".", ".") == 0.0


def test_metrics_reject_mismatched_lengths():
    with pytest.raises(ValueError, match="same length"):
        base_pair_f1("...", "....")
    with pytest.raises(ValueError, match="same length"):
        squared_mountain_distance("...", "....")
    with pytest.raises(ValueError, match="same length"):
        normalized_squared_mountain_distance("...", "....")
