from mountain_centroid.formatting import dot_bracket_from_pairs, pairs_from_bracket


def test_dot_bracket_round_trip():
    pairs = [(1, 8), (2, 5)]
    structure = dot_bracket_from_pairs(8, pairs)

    assert structure == "((..)..)"
    assert pairs_from_bracket(structure) == pairs

