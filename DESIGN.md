# Architecture and design

## Objective

For a supplied ensemble mean mountain profile
\(\mu=(\mu_1,\ldots,\mu_{n-1})\), Mountain Centroid minimizes

\[
J(\sigma;\mu)=\sum_{k=1}^{n-1}(h_\sigma(k)-\mu_k)^2
\]

over a specified prediction space. Here \(h_\sigma(k)\) is the number of base
pairs crossing the cut after nucleotide \(k\). BPP inference and optimization
over candidate structures are separate stages.

### Hybrid extension

In the pairability-constrained space, the hybrid solver minimizes

\[
Q_\alpha(\sigma)=(1-\alpha)\frac{J(\sigma;\mu)}{M}
-\alpha\frac{G(\sigma)}{P},
\qquad
G(\sigma)=\sum_{(i,j)\in\sigma}(2p_{ij}-1),
\]

where \(M=\sum_{k=1}^{n-1}\min(k,n-k)^2\), \(P=\lfloor n/2\rfloor\), and
\(0\leq\alpha\leq1\). Thus \(\alpha=0\) is Mountain Centroid and
\(\alpha=1\) is the base-pair centroid (\(\gamma=1\)) endpoint. The pair gain
enters the same interval--external-depth recurrence, so the constrained
solver's complexity is unchanged.

## Prediction spaces

### Geometry-only

The geometry-only dynamic program searches nonnegative paths that:

- start and end at height zero;
- change by \(-1\), \(0\), or \(1\) at each nucleotide;
- stay below the positional bound \(\min(k,n-k)\).

The returned path corresponds to a pseudoknot-free dot-bracket structure, but
it need not satisfy pairability or minimum-hairpin constraints for the input
sequence. The solver uses `O(nH)` time and traceback space, where `H` is the
maximum reachable height.

### Pairability-constrained

The interval--external-depth state `F(i, j, d)` stores the minimum cost within
interval `[i, j]` when `d` outside pairs enclose it. Candidate structures:

- permit AU/UA, GC/CG, and GU/UG pairs;
- require `TURN = 3`;
- exclude pseudoknots.

With `D_eff` reached external-depth levels, the implementation uses
`O(D_eff n^3)` time and `O(D_eff n^2)` memoization space. The corresponding
worst-case bounds are `O(n^4)` time and `O(n^3)` space.

The Python implementation is the readable reference. The C++ implementation is
the production backend and follows the same recurrence and deterministic
tie-breaking rule.

## BPP backend boundary

`bpp_mu.py` converts a BPP matrix into the ensemble mean mountain profile using

\[
\mu_k=\sum_{i\leq k<j}p_{ij}.
\]

ViennaRNA is the default backend; LinearPartition-V is an optional faster
approximation.
LinearPartition beam pruning affects the supplied BPPs and mean profile, not
the Mountain Centroid dynamic programs.

## Public interfaces

The primary command-line interface returns a pairability-constrained
prediction:

```text
mountain-centroid --seq SEQUENCE
mountain-centroid --seq SEQUENCE --bpp-beam-size 100
mountain-centroid --seq SEQUENCE --bpp-backend vienna
mountain-centroid --seq SEQUENCE --alpha 0 --alpha 0.5 --alpha 1
```

The public Python API additionally exposes:

- `predict` and `predict_from_profile`;
- `predict_hybrid` and `predict_hybrid_curve`;
- `relaxed_mountain_centroid`;
- Python and C++ pairability-constrained solvers;
- evaluation and dot-bracket helpers.

Research datasets, manuscript figures, and study-specific statistical analyses
do not belong in this repository.

## Validation

Correctness checks include:

- exhaustive enumeration on small random instances;
- Python/C++ equality of structures, objective values, and diagnostics;
- pairability, `TURN`, balance, and noncrossing invariants;
- direct recomputation of the objective from returned profiles;
- exhaustive hybrid-objective checks, including both endpoints;
- the geometry-only lower-bound relationship;
- LinearPartition and public-API smoke tests.

## Out of scope

The public inference interface does not include pseudoknot prediction,
alternative MIQP/MILP objectives, or study-specific sampling analyses. Earlier
research prototypes remain available in Git history.
