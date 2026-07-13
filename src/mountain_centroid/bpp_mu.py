#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
bpp_mu.py

ViennaRNA (RNAlib) を用いて base-pair probabilities (bpp) と
ensemble mountain height μ を計算するモジュール。
"""
from __future__ import annotations
from typing import List, Tuple
import numpy as np

try:
    from .utils_tri import tri_to_full
except Exception:
    from utils_tri import tri_to_full


def compute_bpp_and_mu(seq: str, temperature: float = 37.0) -> Tuple[np.ndarray, List[float]]:
    """
    seq: RNA 配列（A/C/G/U）。長さ n。
    戻り値:
      - bpp: NxN 行列 (0-based, 上三角に p_ij, それ以外は 0)
      - mu:  [μ_1, μ_2, ..., μ_{n-1}] where μ_k = sum_{i<=k<j} p_ij
    """
    import RNA  # ViennaRNA の Python バインディング

    seq = seq.upper().replace('T', 'U')
    n = len(seq)
    if n < 2:
        return np.zeros((n, n), dtype=float), [0.0] * (max(n - 1, 0))

    md = RNA.md()
    md.temperature = float(temperature)

    fc = RNA.fold_compound(seq, md)
    fc.pf()
    bpp_tri = fc.exp_matrices.probs
    bpp = tri_to_full(n, bpp_tri, fc.iindx)

    # μ を差分法で効率的に計算
    diff = [0.0] * (n + 1)  # k=0..n（k=1..n-1 を最終利用）
    for i in range(1, n):  # 1..n-1 (1-based)
        for j in range(i + 1, n + 1):  # i+1..n
            pij = bpp[i - 1, j - 1]
            diff[i] += pij
            diff[j] -= pij

    mu = [0.0] * (n - 1)
    acc = 0.0
    for k in range(1, n):
        acc += diff[k]
        mu[k - 1] = acc

    return bpp, mu
