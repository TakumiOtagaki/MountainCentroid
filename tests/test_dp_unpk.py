from mountain_centroid.dp_unpk import dp_nearest_mountain


def test_zero_profile_produces_unpaired_structure():
    heights, pairs, objective = dp_nearest_mountain([0.0, 0.0])

    assert heights == [0, 0, 0, 0]
    assert pairs == []
    assert objective == 0.0


def test_single_arch_profile_is_recovered():
    mu = [0.8, 1.2, 0.1]
    heights, pairs, objective = dp_nearest_mountain(mu)

    assert heights == [0, 1, 1, 0, 0]
    assert pairs == [(1, 3)]
    direct_squared_error = sum(
        (height - expected) ** 2
        for height, expected in zip(heights[1:-1], mu)
    )
    assert objective == direct_squared_error
