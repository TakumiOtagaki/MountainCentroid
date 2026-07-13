"""Standard pseudoknot-free dot-bracket conversion helpers."""

from __future__ import annotations

from typing import Iterable


def dot_bracket_from_pairs(
    n: int,
    pairs: Iterable[tuple[int, int]],
) -> str:
    """Convert 1-based, noncrossing pairs to standard dot-bracket notation."""
    characters = ["."] * n
    normalized_pairs = sorted(pairs)
    occupied: set[int] = set()
    for left, right in normalized_pairs:
        if not 1 <= left < right <= n:
            raise ValueError(f"Invalid pair ({left}, {right}) for length {n}")
        if left in occupied or right in occupied:
            raise ValueError("Each nucleotide can occur in at most one pair")
        occupied.update((left, right))
        characters[left - 1] = "("
        characters[right - 1] = ")"

    structure = "".join(characters)
    if pairs_from_bracket(structure) != normalized_pairs:
        raise ValueError("Pairs are crossing and cannot use standard dot-bracket")
    return structure


def pairs_from_bracket(structure: str) -> list[tuple[int, int]]:
    """Parse standard pseudoknot-free dot-bracket into 1-based pairs."""
    stack: list[int] = []
    pairs: list[tuple[int, int]] = []
    for position, character in enumerate(structure, start=1):
        if character == ".":
            continue
        if character == "(":
            stack.append(position)
            continue
        if character == ")":
            if not stack:
                raise ValueError(f"Unbalanced closing bracket at position {position}")
            pairs.append((stack.pop(), position))
            continue
        raise ValueError(f"Unknown dot-bracket character: {character}")
    if stack:
        raise ValueError(f"Unbalanced opening brackets at positions {stack}")
    pairs.sort()
    return pairs
