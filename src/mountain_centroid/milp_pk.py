#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
milp_pk.py

擬似結び目ありの L1-MILP 定式化（PuLP / CBC）。
"""
from __future__ import annotations
from typing import Dict, List, Tuple, Optional


def can_pair(a: str, b: str, allow_gu: bool = True) -> bool:
    a = a.upper(); b = b.upper()
    if (a == 'A' and b == 'U') or (a == 'U' and b == 'A'):
        return True
    if (a == 'C' and b == 'G') or (a == 'G' and b == 'C'):
        return True
    if allow_gu and ((a == 'G' and b == 'U') or (a == 'U' and b == 'G')):
        return True
    return False


def candidate_pairs(seq: str, lmin: int = 3, allow_gu: bool = True, max_span: Optional[int] = None) -> List[Tuple[int, int]]:
    """
    候補ペア集合 A を作る（1-based, i<j）。
    lmin: 最小ヘアピン長（j - i - 1 >= lmin）
    max_span: j - i <= max_span の制限（None なら制限無し）
    """
    seq = seq.upper().replace('T', 'U')
    n = len(seq)
    A: List[Tuple[int, int]] = []
    for i in range(1, n):  # 1..n-1
        for j in range(i + 1, n + 1):
            if (j - i - 1) < lmin:
                continue
            if max_span is not None and (j - i) > max_span:
                continue
            if can_pair(seq[i - 1], seq[j - 1], allow_gu=allow_gu):
                A.append((i, j))
    return A


def milp_pk_l1(mu: List[float], seq: str, lmin: int = 3, allow_gu: bool = True,
               max_span: Optional[int] = None, lam_sparse: float = 0.0) -> Tuple[List[Tuple[int, int]], float]:
    """
    擬似結び目を許す L1-MILP:
      変数 z_ij ∈ {0,1}, y_k ≥0, e_k ≥0
      目的 min Σ e_k + λ Σ z_ij
      s.t. y_k = Σ_{i<=k<j} z_ij
           |y_k - μ_k| ≤ e_k
           各塩基の次数 ≤ 1
    """
    try:
        import pulp
    except ImportError as e:
        raise RuntimeError("PuLP が見つかりません。pip install pulp を実行してください。") from e

    n = len(seq)
    A = candidate_pairs(seq, lmin=lmin, allow_gu=allow_gu, max_span=max_span)

    # カバレッジの前計算
    cover_by_k: Dict[int, List[Tuple[int, int]]] = {k: [] for k in range(1, n)}
    left_of: Dict[int, List[Tuple[int, int]]] = {i: [] for i in range(1, n + 1)}
    right_of: Dict[int, List[Tuple[int, int]]] = {j: [] for j in range(1, n + 1)}
    for (i, j) in A:
        for k in range(i, j):  # i <= k < j
            if 1 <= k <= n - 1:
                cover_by_k[k].append((i, j))
        left_of[i].append((i, j))
        right_of[j].append((i, j))

    # MILP
    prob = pulp.LpProblem("pk_l1", pulp.LpMinimize)

    z = {(i, j): pulp.LpVariable(f"z_{i}_{j}", lowBound=0, upBound=1, cat="Binary") for (i, j) in A}
    y = {k: pulp.LpVariable(f"y_{k}", lowBound=0, cat="Continuous") for k in range(1, n)}
    e = {k: pulp.LpVariable(f"e_{k}", lowBound=0, cat="Continuous") for k in range(1, n)}

    # 高さ定義 & 絶対値線形化
    for k in range(1, n):
        prob += y[k] == pulp.lpSum(z[i, j] for (i, j) in cover_by_k[k])
        prob += y[k] - mu[k - 1] <= e[k]
        prob += mu[k - 1] - y[k] <= e[k]

    # 各塩基は高々1回
    for i in range(1, n + 1):
        prob += pulp.lpSum(z[pair] for pair in left_of[i]) + pulp.lpSum(z[pair] for pair in right_of[i]) <= 1

    # 目的
    obj = pulp.lpSum(e[k] for k in range(1, n))
    if lam_sparse > 0.0:
        obj += lam_sparse * pulp.lpSum(z[i, j] for (i, j) in A)
    prob += obj

    solver = pulp.PULP_CBC_CMD(msg=False)
    prob.solve(solver)

    pairs = [(i, j) for (i, j) in A if pulp.value(z[i, j]) is not None and pulp.value(z[i, j]) > 0.5]
    pairs.sort()
    obj_val = float(pulp.value(prob.objective))
    return pairs, obj_val
