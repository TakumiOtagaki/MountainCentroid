import pytest

from mountain_centroid.formatting import dot_bracket_from_pairs, pairs_from_bracket


def test_dot_bracket_round_trip():
    pairs = [(1, 8), (2, 5)]
    structure = dot_bracket_from_pairs(8, pairs)

    assert structure == "((..)..)"
    assert pairs_from_bracket(structure) == pairs


def test_crossing_pairs_are_rejected():
    with pytest.raises(ValueError, match="crossing"):
        dot_bracket_from_pairs(4, [(1, 3), (2, 4)])
