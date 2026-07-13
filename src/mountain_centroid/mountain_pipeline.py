#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mountain_pipeline.py

1) ViennaRNA (RNAlib) で base-pair probabilities p_ij を計算し、
   ensemble mountain height μ_k = E[h(k)] を取得。
2) Pseudoknot なし:
   - DP（L2 厳密）で μ に最も近い mountain path を求め、括弧列（dot-bracket）を出力。
   - （オプション）MIQP（高さベース）で同じ目的を解く（Gurobi/CPLEX 等があれば）。
3) Pseudoknot あり:
   - MILP（L1）で μ に最も近い構造（交差許容）を求める。
     候補ペア集合 A は WC/GU と最小ループ長で事前に絞る。

依存:
  - ViennaRNA (pip install ViennaRNA もしくは conda install -c bioconda viennarna)
  - PuLP (pip install pulp)
  - 任意: gurobipy / cplex (ある場合のみ MIQP を実行)

著者メモ:
  - すべて 1-based 的な定義（i<k<j ⇔ cut k を跨ぐ）に合わせ、内部では必要に応じて 0/1-based を変換しています。
  - E[h] は k=1..n-1 の長さ n-1 ベクトルです。
"""

from __future__ import annotations
import argparse
from typing import List, Tuple, Optional

# 相対インポート（パッケージとして実行）/ フォールバック（スクリプト直実行）
try:
    from .bpp_mu import compute_bpp_and_mu
    from .dp_unpk import dp_nearest_mountain
    from .miqp_unpk import miqp_unpk_height
    from .milp_pk import milp_pk_l1
    from .formatting import dot_bracket_from_pairs, bracket_with_pseudoknots
except Exception:  # pragma: no cover - 直実行時のフォールバック
    from bpp_mu import compute_bpp_and_mu
    from dp_unpk import dp_nearest_mountain
    from miqp_unpk import miqp_unpk_height
    from milp_pk import milp_pk_l1
    from formatting import dot_bracket_from_pairs, bracket_with_pseudoknots

# ----------------------------
# ViennaRNA (RNAlib) で E[h] を計算
# ----------------------------

def compute_bpp_and_mu_cli(seq: str, temperature: float = 37.0):
    return compute_bpp_and_mu(seq, temperature)


# ----------------------------
# 擬似結び目なし: DP（L2 厳密）
# ----------------------------

dp_nearest_mountain_cli = dp_nearest_mountain


dot_bracket_from_pairs_cli = dot_bracket_from_pairs


# ----------------------------
# 擬似結び目なし: 整数計画（MIQP, 高さベース）
# ----------------------------

miqp_unpk_height_cli = miqp_unpk_height


# ----------------------------
# 擬似結び目あり: MILP（L1）
# ----------------------------

milp_pk_l1_cli = milp_pk_l1


# ----------------------------
# 出力整形（擬似結び目ありの簡易多括弧表記）
# ----------------------------

bracket_with_pseudoknots_cli = bracket_with_pseudoknots


# ----------------------------
# CLI
# ----------------------------

def main():
    ap = argparse.ArgumentParser(description="E[h]→（非）擬似結び目の構造推定（DP/MIQP/MILP）")
    ap.add_argument("--seq", required=True, help="RNA sequence (A/C/G/U)")
    ap.add_argument("--temp", type=float, default=37.0, help="ViennaRNA temperature (°C)")
    ap.add_argument("--band", type=int, default=None, help="DP のバンド幅（省略で全域）")
    ap.add_argument("--unpk-miqp", action="store_true", help="非疑似結び目 MIQP を実行（Gurobi等がある場合）")
    ap.add_argument("--pk-l1", action="store_true", help="擬似結び目あり L1-MILP を実行")
    ap.add_argument("--lmin", type=int, default=3, help="MILP の最小ヘアピン長（j-i-1 >= lmin）[default 3]")
    # ap.add_argument("--allow-gu", action="store_true", help="MILP の候補に GU を含める")
    allow_gu_default = True
    ap.add_argument("--max-span", type=int, default=None, help="MILP の最大スパン（j-i <= max_span）")
    ap.add_argument("--lam-sparse", type=float, default=0.0, help="MILP のペア数ペナルティ λ（疎にしたいとき）")
    args = ap.parse_args()

    seq = args.seq.upper().replace('T','U')
    n = len(seq)
    print(f"# seq (n={n}): {seq}")

    # 1) E[h]
    bpp, mu = compute_bpp_and_mu_cli(seq, temperature=args.temp)
    print(f"# E[h] length = {len(mu)} (k=1..{n-1})")
    # 先頭数点だけ表示
    preview = ", ".join(f"{x:.3f}" for x in mu[:min(10, len(mu))])
    print(f"# E[h] preview: [{preview}{' ...' if len(mu)>10 else ''}]")

    # 2) 非擬似結び目：DP
    t_dp, pairs_dp, obj_dp = dp_nearest_mountain_cli(mu, band=args.band)
    dbn_dp = dot_bracket_from_pairs_cli(n, pairs_dp)
    print("\n[UNPK][DP L2] objective =", f"{obj_dp:.6f}")
    print("[UNPK][DP L2] dot-bracket:", dbn_dp)

    # 2) 非擬似結び目：MIQP（任意）
    if args.unpk_miqp:
        res = miqp_unpk_height_cli(mu)
        if res is None:
            print("[UNPK][MIQP] solver not found or not optimal -> skipped")
        else:
            t_miqp, pairs_miqp, obj_miqp = res
            dbn_miqp = dot_bracket_from_pairs_cli(n, pairs_miqp)
            print("\n[UNPK][MIQP L2] objective =", f"{obj_miqp:.6f}")
            print("[UNPK][MIQP L2] dot-bracket:", dbn_miqp)

    # 3) 擬似結び目あり：MILP（L1）
    if args.pk_l1:
        pairs_pk, obj_pk = milp_pk_l1_cli(
            mu, seq, lmin=args.lmin, allow_gu=allow_gu_default,
            max_span=args.max_span, lam_sparse=args.lam_sparse
        )
        dbn_pk = bracket_with_pseudoknots_cli(n, pairs_pk)
        print("\n[PK][MILP L1] objective =", f"{obj_pk:.6f}")
        print("[PK][MILP L1] #pairs:", len(pairs_pk))
        print("[PK][MILP L1] pseudo-bracket:", dbn_pk)
        # 先頭いくつかのペアを表示
        print("[PK][MILP L1] sample pairs:", pairs_pk[:min(10, len(pairs_pk))])

if __name__ == "__main__":
    main()