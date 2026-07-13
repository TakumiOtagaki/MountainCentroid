from mountain_centroid.dp_unpk import dp_nearest_mountain


def test_zero_profile_produces_unpaired_structure():
    heights, pairs, objective = dp_nearest_mountain([0.0, 0.0])

    assert heights == [0, 0, 0, 0]
    assert pairs == []
    assert objective == 0.0


def test_single_arch_profile_is_recovered():
    heights, pairs, _ = dp_nearest_mountain([1.0, 1.0, 0.0])

    assert heights == [0, 1, 1, 0, 0]
    assert pairs == [(1, 3)]

