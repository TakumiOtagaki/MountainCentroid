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
mountain-centroid --seq SEQUENCE --beam-size 200
mountain-centroid --seq SEQUENCE --bpp-backend linearpartition --bpp-beam-size 100
```

Algorithm-development switches, reference solvers, and paper-only sampling
analyses should not become user-facing modes. If an exact small-instance oracle
is needed for tests, it belongs under `tests/` rather than the installed CLI.

## Implementation status

The public solver is a left-to-right beam search over persistent stacks of open
base pairs. It enforces sequence-dependent pair admissibility, minimum hairpin
length 3, and pseudoknot-free nesting while minimizing accumulated squared
mountain-profile error. Candidate generation is O(nB); sorting the beam in the
current Python implementation gives O(n B log B). It is therefore linear in
sequence length for fixed B. The result is approximate whenever pruning drops a
prefix needed by the global optimum.

An exhaustive small-instance oracle lives only in the test suite. It verifies
that the beam implementation recovers the exact constrained optimum when the
beam is large enough; it is not a user-facing inference mode.

## Remaining work

Work remaining, in priority order:

1. Add broader randomized invariant tests and quantify optimality agreement as
   a function of solver beam size on oracle-sized inputs.
2. Add reproducible ViennaRNA-versus-LinearPartition agreement and runtime
   benchmarks.
3. Implement the Boltzmann-sampling rooted-L2 evaluation in the paper analysis
   code, not as an inference mode.
4. Freeze datasets, versions, seeds, and commands for the TBIO submission.
