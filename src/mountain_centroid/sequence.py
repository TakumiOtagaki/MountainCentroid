"""RNA sequence normalization and pairing rules used by the estimator."""

from __future__ import annotations


CANONICAL_AND_WOBBLE_PAIRS = frozenset(
    {"AU", "UA", "GC", "CG", "GU", "UG"}
)
MIN_HAIRPIN_LENGTH = 3


def normalise_sequence(sequence: str) -> str:
    """Return an uppercase RNA sequence and reject unsupported symbols."""
    sequence = sequence.upper().replace("T", "U")
    invalid = sorted(set(sequence) - set("ACGU"))
    if invalid:
        raise ValueError(
            f"Sequence contains unsupported symbols: {''.join(invalid)}"
        )
    if not sequence:
        raise ValueError("Sequence must not be empty")
    return sequence


def can_pair(left: str, right: str) -> bool:
    """Return whether two normalized bases form AU, GC, or GU pairs."""
    return left + right in CANONICAL_AND_WOBBLE_PAIRS

