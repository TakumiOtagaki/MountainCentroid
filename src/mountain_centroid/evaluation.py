"""Evaluation metrics for pseudoknot-free RNA secondary structures."""

from __future__ import annotations

from .formatting import pairs_from_bracket


def base_pair_f1(predicted: str, reference: str) -> float:
    """Return exact-match base-pair F1 for two dot-bracket structures."""
    if len(predicted) != len(reference):
        raise ValueError("Structures must have the same length")
    predicted_pairs = set(pairs_from_bracket(predicted))
    reference_pairs = set(pairs_from_bracket(reference))
    true_positive = len(predicted_pairs & reference_pairs)
    denominator = len(predicted_pairs) + len(reference_pairs)
    return 2.0 * true_positive / denominator if denominator else 0.0


def mountain_heights(structure: str) -> tuple[int, ...]:
    """Return cut heights h(1),...,h(n-1) from standard dot-bracket."""
    depth = 0
    heights: list[int] = []
    for position, character in enumerate(structure):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth < 0:
                raise ValueError(
                    f"Unbalanced closing bracket at position {position + 1}"
                )
        elif character != ".":
            raise ValueError(f"Unknown dot-bracket character: {character}")
        if position < len(structure) - 1:
            heights.append(depth)
    if depth:
        raise ValueError("Unbalanced opening brackets")
    return tuple(heights)


def squared_mountain_distance(predicted: str, reference: str) -> float:
    """Return sum_k (h_pred(k)-h_ref(k))^2 without taking a square root."""
    if len(predicted) != len(reference):
        raise ValueError("Structures must have the same length")
    predicted_heights = mountain_heights(predicted)
    reference_heights = mountain_heights(reference)
    return float(
        sum(
            (predicted_height - reference_height) ** 2
            for predicted_height, reference_height in zip(
                predicted_heights,
                reference_heights,
            )
        )
    )


def mean_squared_mountain_distance(predicted: str, reference: str) -> float:
    """Return squared mountain distance normalized by the number of cuts."""
    number_of_cuts = max(len(reference) - 1, 1)
    return squared_mountain_distance(predicted, reference) / number_of_cuts


def normalized_squared_mountain_distance(predicted: str, reference: str) -> float:
    """Return squared distance divided by the maximum path discrepancy.

    At cut k of a length-n structure, every valid mountain height is bounded by
    min(k, n-k).  Dividing by the sum of these squared bounds gives a
    reference-independent, dimensionless value in [0, 1].
    """
    if len(predicted) != len(reference):
        raise ValueError("Structures must have the same length")
    n = len(reference)
    maximum_squared_distance = sum(
        min(cut, n - cut) ** 2 for cut in range(1, n)
    )
    if maximum_squared_distance == 0:
        return 0.0
    return squared_mountain_distance(predicted, reference) / maximum_squared_distance
