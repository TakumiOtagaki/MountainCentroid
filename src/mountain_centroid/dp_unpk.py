#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
dp_unpk.py

擬似結び目なし（unpk）で μ に最も近い mountain path を L2 厳密 DP で求める。
"""
from __future__ import annotations
from typing import Dict, List, Tuple


def dp_nearest_mountain(mu: List[float]) -> Tuple[List[int], List[Tuple[int, int]], float]:
    """
    L2: min_{t} sum_{k=1}^{n-1} (t_k - mu_k)^2
    s.t. t_0=0, t_n=0, t_k>=0, |t_k - t_{k-1}|<=1
    戻り値:
      - t: 高さ列 t_0..t_n（長さ n+1）
      - pairs: 非疑似結び目のペア集合（1-based）
      - obj: 目的関数値（L2）
    """
    n = len(mu) + 1
    INF = 1e100

    # k ごとの高さ上限（理論上）
    hmax = [0] * (n + 1)
    for k in range(0, n + 1):
        hmax[k] = min(k, n - k)

    # DP テーブル（ローリング）
    prev = [INF] * (hmax[0] + 1)
    prev[0] = 0.0
    # バックポインタ（各 k に辞書で）
    back: List[Dict[int, int]] = [dict() for _ in range(n + 1)]

    for k in range(1, n + 1):
        cur = [INF] * (hmax[k] + 1)
        mu_k = 0.0 if k == n else mu[k - 1]
        for t in range(hmax[k] + 1):
            best_cost = INF
            best_s = None
            # 遷移元 s ∈ {t-1, t, t+1} ∩ [0, hmax[k-1]]
            for s in (t, t - 1, t + 1):
                if 0 <= s <= hmax[k - 1]:
                    # val = prev[s] + (t - mu_k) * (t - mu_k)
                    val = prev[s] + t ** 2 - 2 * t * mu_k 
                    if val < best_cost:
                        best_cost = val
                        best_s = s
            if best_s is not None:
                cur[t] = best_cost
                back[k][t] = best_s
        prev = cur

    # 終点は t_n=0 を強制
    t = [0] * (n + 1)
    t[n] = 0
    for k in range(n, 0, -1):
        t[k - 1] = back[k][t[k]]

    # 括弧復元（非疑似結び目：LIFO）
    stack: List[int] = []
    pairs: List[Tuple[int, int]] = []
    for k in range(1, n + 1):
        dt = t[k] - t[k - 1]
        if dt == +1:
            stack.append(k)
        elif dt == -1:
            if not stack:
                raise RuntimeError("Invalid path: negative pop")
            i = stack.pop()
            pairs.append((i, k))
    pairs.sort()
    obj = sum((t[k] - mu[k - 1]) ** 2 for k in range(1, n))
    return t, pairs, float(obj)
