#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
mountain_pipeline.py

1) ViennaRNA または LinearPartition-V で base-pair probabilities p_ij を計算し、
   ensemble mountain height μ_k = E[h(k)] を取得。
2) μ に最も近い pseudoknot-free mountain path を求め、dot-bracket を出力。

依存:
  - ViennaRNA (pip install ViennaRNA もしくは conda install -c bioconda viennarna)
現在の path DP は sequence constraints をまだ実装していない開発用 scaffold。
公開インターフェースは最終的な単一手法を想定し、過去の MIQP/MILP 分岐は持たない。
"""

from __future__ import annotations
import argparse

# 相対インポート（パッケージとして実行）/ フォールバック（スクリプト直実行）
try:
    from .bpp_mu import compute_bpp_and_mu
    from .dp_unpk import dp_nearest_mountain
    from .formatting import dot_bracket_from_pairs
except Exception:  # pragma: no cover - 直実行時のフォールバック
    from bpp_mu import compute_bpp_and_mu
    from dp_unpk import dp_nearest_mountain
    from formatting import dot_bracket_from_pairs

# ----------------------------
# BPP backend で E[h] を計算
# ----------------------------

def compute_bpp_and_mu_cli(seq: str, temperature: float = 37.0, **kwargs):
    return compute_bpp_and_mu(seq, temperature, **kwargs)


# ----------------------------
# 擬似結び目なし: mountain-path DP
# ----------------------------

dp_nearest_mountain_cli = dp_nearest_mountain


dot_bracket_from_pairs_cli = dot_bracket_from_pairs


# ----------------------------
# CLI
# ----------------------------

def main():
    ap = argparse.ArgumentParser(
        description="Expected mountain profile to one pseudoknot-free structure"
    )
    ap.add_argument("--seq", required=True, help="RNA sequence (A/C/G/U)")
    ap.add_argument(
        "--temp",
        type=float,
        default=37.0,
        help="temperature in °C (ViennaRNA only; LinearPartition-V is fixed at 37°C)",
    )
    ap.add_argument(
        "--bpp-backend",
        choices=("vienna", "linearpartition"),
        default="vienna",
        help="BPP backend: exact ViennaRNA (default) or approximate LinearPartition-V",
    )
    ap.add_argument("--beam-size", type=int, default=100, help="LinearPartition beam size [default 100]")
    ap.add_argument("--bpp-cutoff", type=float, default=0.0, help="LinearPartition BPP output cutoff")
    ap.add_argument("--linearpartition-path", default=None, help="Path to LinearPartition runner script")
    args = ap.parse_args()

    seq = args.seq.upper().replace('T','U')
    n = len(seq)
    print(f"# seq (n={n}): {seq}")

    # 1) E[h]
    bpp, mu = compute_bpp_and_mu_cli(
        seq,
        temperature=args.temp,
        backend=args.bpp_backend,
        beam_size=args.beam_size,
        cutoff=args.bpp_cutoff,
        linearpartition_path=args.linearpartition_path,
    )
    print(f"# BPP backend: {args.bpp_backend}")
    print(f"# E[h] length = {len(mu)} (k=1..{n-1})")
    # 先頭数点だけ表示
    preview = ", ".join(f"{x:.3f}" for x in mu[:min(10, len(mu))])
    print(f"# E[h] preview: [{preview}{' ...' if len(mu)>10 else ''}]")

    # 2) 非擬似結び目：DP
    _, pairs_dp, obj_dp = dp_nearest_mountain_cli(mu)
    dbn_dp = dot_bracket_from_pairs_cli(n, pairs_dp)
    print("\n# inference status: relaxed mountain-path scaffold; sequence constraints pending")
    print("squared_mountain_error =", f"{obj_dp:.6f}")
    print("dot_bracket =", dbn_dp)

if __name__ == "__main__":
    main()
