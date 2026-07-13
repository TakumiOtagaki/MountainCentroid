#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
miqp_unpk.py

擬似結び目なし（unpk）の高さベース MIQP（Gurobi がある環境のみ）。
"""
from __future__ import annotations
from typing import List, Tuple, Optional


def miqp_unpk_height(mu: List[float]) -> Optional[Tuple[List[int], List[Tuple[int, int]], float]]:
    """
    L2 norm の最小化をします
    高さベース MIQP（凸）:
      min sum_k (t_k - mu_k)^2
      s.t. t_0=0, t_n=0, t_k>=0, t_k - t_{k-1} = u_k^+ - u_k^- , u_k^+ + u_k^- <= 1
    Gurobi が見つかれば実行。無ければ None を返す。
    """
    n = len(mu) + 1
    try:
        import gurobipy as gp
        GRB = gp.GRB

        m = gp.Model("unpk_miqp")
        m.setParam("OutputFlag", 0)

        # 変数
        t = m.addVars(range(0, n + 1), vtype=GRB.INTEGER, lb=0, name="t")
        up = m.addVars(range(1, n + 1), vtype=GRB.BINARY, name="up")
        dn = m.addVars(range(1, n + 1), vtype=GRB.BINARY, name="dn")

        # 制約
        m.addConstr(t[0] == 0)
        m.addConstr(t[n] == 0)
        for k in range(1, n + 1):
            m.addConstr(t[k] - t[k - 1] == up[k] - dn[k])
            m.addConstr(up[k] + dn[k] <= 1)
            m.addConstr(t[k] <= min(k, n - k))

        # 目的（二次, 凸）
        obj = gp.QuadExpr()
        for k in range(1, n):  # k=1..n-1 にコスト
            obj += (t[k] - mu[k - 1]) * (t[k] - mu[k - 1])
        m.setObjective(obj, GRB.MINIMIZE)

        m.optimize()
        if m.status != GRB.OPTIMAL:
            return None

        t_sol = [int(round(t[k].X)) for k in range(0, n + 1)]
        # 括弧復元（非疑似結び目）
        stack: List[int] = []
        pairs: List[Tuple[int, int]] = []
        for k in range(1, n + 1):
            dt = t_sol[k] - t_sol[k - 1]
            if dt == +1:
                stack.append(k)
            elif dt == -1:
                i = stack.pop()
                pairs.append((i, k))
        pairs.sort()

        obj_val = sum((t_sol[k] - mu[k - 1]) ** 2 for k in range(1, n))
        return t_sol, pairs, float(obj_val)

    except Exception:
        return None
