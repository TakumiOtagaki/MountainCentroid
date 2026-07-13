#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
utils_tri.py

ViennaRNA の上三角一次元配列アクセスと復元ユーティリティ。
"""
from __future__ import annotations
from typing import Optional
import numpy as np


def tri_get(arr, i: int, j: int, iindx) -> float:
    """三角配列から (i,j) を取得（1-based, i<j）。C と同じ iindx[i] - j アクセス。
    """
    if i >= j:
        raise ValueError("Require i < j for upper-triangular arrays")
    return arr[iindx[i] - j]


def tri_get_opt(arr, i: int, j: int, iindx) -> Optional[float]:
    """三角配列が None でも動く安全版。None の場合は None を返す。"""
    if arr is None:
        return None
    return tri_get(arr, i, j, iindx)


def tri_to_full(n: int, tri_arr, iindx) -> np.ndarray:
    """ViennaRNAの上三角一次元配列(tri_arr)をNxNのfloat64行列に復元する。
    i,jは1-basedでtri_getを使って取得し、0-basedの[i-1,j-1]に格納。下三角と対角は0。
    """
    out = np.zeros((n, n), dtype=np.float64)
    if tri_arr is None:
        return out
    for i in range(1, n):  # 1..n-1
        for j in range(i + 1, n + 1):  # i+1..n
            out[i - 1, j - 1] = tri_get(tri_arr, i, j, iindx)
    return np.ascontiguousarray(out)
