# Public software design

Status: working policy (2026-07-14)

## Two explicit prediction spaces

The public package exposes the same expected squared mountain-loss objective
over two explicitly named prediction spaces. The relaxed projection
projects onto nonnegative unit-step mountain paths without pairability or
minimum-hairpin constraints. Its position-height DP takes O(nH) time and O(nH)
traceback space. The sequence-constrained estimator additionally requires:

- Watson-Crick pairs (AU/UA and GC/CG) plus GU/UG wobble pairs only;
- minimum hairpin length 3;
- no pseudoknots.

LinearPartition-V (fast approximation, public default) and ViennaRNA (optional
partition-function backend) are BPP backends for the same estimator. They are not separate
inference modes.

The former height-only MIQP and pseudoknot-aware L1-MILP branches were research
prototypes. They are intentionally absent from the public CLI. Their last
version remains available in Git history at commit `accde04`.

## Interface target

The stable interface should remain close to:

```text
mountain-centroid --seq SEQUENCE
mountain-centroid --seq SEQUENCE --bpp-beam-size 100
mountain-centroid --seq SEQUENCE --bpp-backend vienna
```

The relaxed solver is available through the public Python function
`relaxed_mountain_centroid`; CLI naming will be finalized with the two-variant
interface. Algorithm-development switches, reference solvers, and paper-only sampling
analyses should not become user-facing modes. The Python constrained DP is a
readable reference implementation; the C++ DP is the production backend.

## Constrained implementation status

The public solver is an interval--external-depth dynamic program. It enforces
sequence-dependent pair admissibility, minimum hairpin length 3, and
pseudoknot-free nesting while minimizing squared mountain-profile error. With
`D_eff` reached external-depth levels it uses `O(D_eff n^3)` time and
`O(D_eff n^2)` memoization space. The C++ backend and Python reference agree
with exhaustive enumeration on short instances.

## Remaining work

Work remaining, in priority order:

1. Complete the RNAstralign comparison against ViennaRNA MFE and the standard
   ViennaRNA centroid, including family-stratified outlier inspection.
2. Implement the Boltzmann-sampling rooted-L2 evaluation in the paper analysis
   code, not as an inference mode.
3. Freeze datasets, versions, seeds, and commands for the TBIO submission.
