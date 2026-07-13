# MountainCentroid

MountainCentroid projects the expected mountain profile of an RNA secondary-
structure ensemble onto a discrete, pseudoknot-free mountain path.

This repository contains only the reusable implementation. Manuscript sources,
figures, and paper-specific outputs live in `MountainCentroidPaper`, which pins
the exact software revision as a Git submodule.

## Provenance

The initial implementation was extracted from
`TakumiOtagaki/ZukerStyleCentroidFold` at commit
`17d39d288ed6900561d93de47a8a1b2b98d2c329`.

## Current scope

- Compute base-pair probabilities and the expected cut-based mountain height
  with exact ViennaRNA or approximate LinearPartition-V.
- Find the nearest pseudoknot-free integer mountain path under squared loss.
- Convert the inferred path to dot-bracket notation.
- Provide optional experimental MIQP and pseudoknot-aware L1-MILP solvers.

The path DP is exact for the relaxed mountain-path problem. It currently does
not constrain recovered pairs by nucleotide complementarity or minimum hairpin
length. Sequence-valid constrained inference will be developed and evaluated
before the paper is finalized.

## Installation

```sh
python -m pip install -e .
```

To enable the optional LinearPartition backend, clone recursively and build the
vendored upstream source:

```sh
git submodule update --init --recursive
make -C vendor/LinearPartition
```

## Usage

```sh
mountain-centroid --seq ACGUACGUACGU
```

ViennaRNA remains the default and reference backend. For long sequences,
LinearPartition-V can be selected with:

```sh
mountain-centroid \
  --seq ACGUACGUACGU \
  --bpp-backend linearpartition \
  --beam-size 100
```

LinearPartition uses beam search and therefore approximates the partition
function and BPPs. Backend and beam size must be reported with experimental
results; results from the two backends should not be silently pooled.

Equivalent module invocation:

```sh
python -m mountain_centroid.mountain_pipeline --seq ACGUACGUACGU
```

## Repository layout

```text
src/mountain_centroid/   reusable implementation and CLI
tests/                   unit tests for the path DP and formatting
vendor/LinearPartition/  pinned optional BPP backend (Git submodule)
```
