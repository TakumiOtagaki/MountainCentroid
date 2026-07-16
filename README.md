# MountainCentroid

MountainCentroid exposes two exact prediction spaces for expected squared
mountain loss: a relaxed projection onto mountain paths and a sequence-valid,
pseudoknot-free projection. The exact constrained solver is implemented in
both Python and C++; the Python implementation serves as a readable oracle.

This repository contains only the reusable implementation. Manuscript sources,
figures, and paper-specific outputs live in `MountainCentroidPaper`, which pins
the exact software revision as a Git submodule.

## Provenance

The initial implementation was extracted from
`TakumiOtagaki/ZukerStyleCentroidFold` at commit
`17d39d288ed6900561d93de47a8a1b2b98d2c329`.

## Current scope

- Compute the exact relaxed Mountain Centroid with a position-height dynamic
  program. This variant does not enforce pairability or minimum hairpin length.
- Compute the exact sequence-constrained Mountain Centroid with an
  interval--external-depth dynamic program in `O(D_eff n^3)` time and
  `O(D_eff n^2)` memoization space (`O(n^4)` and `O(n^3)` worst case).
- Compute base-pair probabilities and the expected cut-based mountain height
  with default LinearPartition-V or optional exact ViennaRNA.
- Enforce Watson-Crick pairs (AU/UA and GC/CG) plus GU/UG wobble pairs,
  minimum hairpin length 3, and no pseudoknots during inference.
- Report standard dot-bracket notation and the direct squared profile error.

The older bidirectional beam implementation remains available for engineering
comparisons, but it is not the exact estimator used for the manuscript's
constrained results. LinearPartition beam pruning is a separate upstream
approximation to the BPP-derived expected profile.

The exact solvers are available as the Python APIs
`relaxed_mountain_centroid`, `exact_mountain_centroid`, and
`cpp_exact_mountain_centroid`. The existing command-line prediction route is
retained for compatibility and still exposes the beam implementation.
Pseudoknot prediction and alternative MIQP/MILP objectives are out of scope.

## Installation

```sh
python -m pip install -e .
```

To use the default LinearPartition backend, clone recursively and build the
vendored upstream source:

```sh
git submodule update --init --recursive
make -C vendor/LinearPartition
make exact
```

## Usage

```sh
mountain-centroid --seq ACGUACGUACGU
```

The solver beam can be changed without selecting a different method:

```sh
mountain-centroid --seq ACGUACGUACGU --beam-size 200
```

LinearPartition-V is the default:

```sh
mountain-centroid \
  --seq ACGUACGUACGU \
  --bpp-beam-size 100
```

The optional exact ViennaRNA backend can be selected with:

```sh
mountain-centroid --seq ACGUACGUACGU --bpp-backend vienna
```

LinearPartition uses beam search and therefore approximates the partition
function and BPPs. Backend and beam size must be reported with experimental
results; results from the two backends should not be silently pooled.

Equivalent module invocation:

```sh
python -m mountain_centroid.mountain_pipeline --seq ACGUACGUACGU
```

Solver-only scaling can be measured reproducibly with:

```sh
uv run python benchmarks/benchmark_solver.py \
  --lengths 100 300 1000 3000 --beam-size 100 --repeats 3
```

The production Python-to-C++ exact route can be benchmarked with:

```sh
uv run python benchmarks/benchmark_cpp_exact.py \
  --lengths 30 50 100 150 200 300 --instances 3
```

Approximation quality on synthetic short instances can be checked against
exhaustive optimization with:

```sh
uv run python benchmarks/benchmark_optimality.py
```

An optional exact-backend timing benchmark separates ViennaRNA BPP calculation,
expected-profile construction, and Mountain Centroid inference:

```sh
uv run python benchmarks/benchmark_vienna_pipeline.py \
  --lengths 50 100 200 400 800 --instances 5 --beam-size 100
```

## Repository layout

```text
src/mountain_centroid/   reusable implementation and CLI
tests/                   constraints, exact-oracle, backend, and formatting tests
benchmarks/              reproducible solver scaling benchmark
vendor/LinearPartition/  pinned default BPP backend (Git submodule)
```
