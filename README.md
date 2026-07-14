# MountainCentroid

MountainCentroid predicts the sequence-valid, pseudoknot-free RNA secondary
structure that approximately minimizes expected squared mountain loss.

This repository contains only the reusable implementation. Manuscript sources,
figures, and paper-specific outputs live in `MountainCentroidPaper`, which pins
the exact software revision as a Git submodule.

## Provenance

The initial implementation was extracted from
`TakumiOtagaki/ZukerStyleCentroidFold` at commit
`17d39d288ed6900561d93de47a8a1b2b98d2c329`.

## Current scope

- Compute base-pair probabilities and the expected cut-based mountain height
  with default LinearPartition-V or optional exact ViennaRNA.
- Search valid structures with bidirectional beam-pruned scans and retain the
  lower-loss result.
- Enforce Watson-Crick pairs (AU/UA and GC/CG) plus GU/UG wobble pairs,
  minimum hairpin length 3, and no pseudoknots during inference.
- Report standard dot-bracket notation and the direct squared profile error.

For solver beam size B, each directional scan generates O(nB) candidates, and
beam sorting makes the current Python implementation O(n B log B). With a
fixed B it is linear in sequence length. Persistent stack states make each
open/close transition O(1), avoiding tuple copies proportional to nesting
depth. Beam pruning is an approximation to the mathematical Fréchet mean;
larger beams retain more prefix states. The default is B=100.

The public software intentionally has one prediction route. LinearPartition-V
is the default fast BPP backend; ViennaRNA is an optional exact backend for the
same route. Pseudoknot prediction and alternative MIQP/MILP objectives are out
of scope.

## Installation

```sh
python -m pip install -e .
```

To use the default LinearPartition backend, clone recursively and build the
vendored upstream source:

```sh
git submodule update --init --recursive
make -C vendor/LinearPartition
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
