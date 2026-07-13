#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
formatting.py

出力整形（dot-bracket 生成、擬似結び目ありの多括弧表現）。
"""
from __future__ import annotations
from typing import List, Tuple, Dict


def dot_bracket_from_pairs(n: int, pairs: List[Tuple[int, int]]) -> str:
    """非疑似結び目のペアを (i,j) 1-based でもらって dot-bracket を返す。"""
    s = ['.'] * n
    for i, j in pairs:
        s[i - 1] = '('
        s[j - 1] = ')'
    return ''.join(s)


_BRACKETS = [
    ('(', ')'), ('[', ']'), ('{', '}'), ('<', '>'),
    ('A', 'a'), ('B', 'b'), ('C', 'c'), ('D', 'd'),
    ('E', 'e'), ('F', 'f'), ('G', 'g'), ('H', 'h')
]


def bracket_with_pseudoknots(n: int, pairs: List[Tuple[int, int]]) -> str:
    """
    任意交差のペアを複数種の括弧に彩色して dot-bracket 風に表示。
    彩色は貪欲で「同じ色内では交差しない」ように割当。
    色数が足りない場合は 'x'/'X' を使用。
    """
    colors: List[List[Tuple[int, int]]] = []  # 各色のペア集合
    for (i, j) in sorted(pairs):
        placed = False
        for col in colors:
            ok = True
            for (a, b) in col:
                if (a < i < b < j) or (i < a < j < b):
                    ok = False
                    break
            if ok:
                col.append((i, j))
                placed = True
                break
        if not placed:
            colors.append([(i, j)])

    s = ['.'] * n
    for c, col in enumerate(colors):
        left, right = ('x', 'X') if c >= len(_BRACKETS) else _BRACKETS[c]
        for (i, j) in col:
            s[i - 1] = left
            s[j - 1] = right
    return ''.join(s)


def pairs_from_bracket(db: str) -> List[Tuple[int, int]]:
    """
    ドットブラケット（複数括弧種対応）から (i,j) 1-based の塩基対集合を返す。
    未対応の括弧は無視せず例外にする（データ品質担保のため）。
    """
    n = len(db)
    stacks: Dict[str, List[int]] = {}
    left2right = {l: r for (l, r) in _BRACKETS}
    right2left = {r: l for (l, r) in _BRACKETS}
    pairs: List[Tuple[int, int]] = []
    for pos, ch in enumerate(db, start=1):
        if ch == '.':
            continue
        if ch in left2right:
            stacks.setdefault(ch, []).append(pos)
        elif ch in right2left:
            lch = right2left[ch]
            st = stacks.get(lch, [])
            if not st:
                raise ValueError(f"Unbalanced bracket at pos {pos}: {ch}")
            i = st.pop()
            pairs.append((i, pos))
        else:
            raise ValueError(f"Unknown bracket char: {ch}")
    # 未閉括弧のチェック
    for lch, st in stacks.items():
        if st:
            raise ValueError(f"Unbalanced bracket, leftover '{lch}': {st}")
    pairs.sort()
    return pairs
