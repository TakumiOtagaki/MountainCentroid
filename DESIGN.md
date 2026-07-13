# Public software design

Status: working policy (2026-07-14)

## One prediction method

The public package will expose one estimator: the sequence-constrained,
pseudoknot-free Fréchet mean under expected squared mountain loss.

The prediction constraints are fixed parts of the method:

- AU, UA, GC, CG, GU, and UG pairs only;
- minimum hairpin length 3;
- no pseudoknots.

ViennaRNA (reference) and LinearPartition-V (fast approximation) are BPP
backends for the same estimator. They are not separate inference modes.

The former height-only MIQP and pseudoknot-aware L1-MILP branches were research
prototypes. They are intentionally absent from the public CLI. Their last
version remains available in Git history at commit `accde04`.

## Interface target

The stable interface should remain close to:

```text
mountain-centroid --seq SEQUENCE
mountain-centroid --seq SEQUENCE --bpp-backend linearpartition --beam-size 100
```

Algorithm-development switches, reference solvers, and paper-only sampling
analyses should not become user-facing modes. If an exact small-instance oracle
is needed for tests, it belongs under `tests/` rather than the installed CLI.

## Implementation status and remaining work

The current dynamic program is an exact projection onto the relaxed space of
nonnegative unit-step mountain paths. It now has a single public route and
reports the actual squared profile error, but it still lacks sequence-dependent
pair admissibility and the minimum hairpin constraint. Consequently it is an
implementation scaffold, not yet the final estimator described in the paper.

Work remaining, in priority order:

1. Replace the relaxed height-path solver with a sequence-valid, pseudoknot-free
   optimizer (or a clearly characterized approximation plus an exact
   small-instance oracle).
2. Add invariant tests for complementarity, hairpin length, noncrossing pairs,
   dot-bracket validity, and equality of the reported objective to the direct
   squared-error calculation.
3. Provide one small Python API around sequence-to-structure prediction; keep
   backend configuration subordinate to it.
4. Add reproducible ViennaRNA-versus-LinearPartition agreement and runtime
   benchmarks.
5. Implement the Boltzmann-sampling rooted-L2 evaluation in the paper analysis
   code, not as an inference mode.
6. Freeze datasets, versions, seeds, and commands for the TBIO submission.
