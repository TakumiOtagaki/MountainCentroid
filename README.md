# MountainCentroid

MountainCentroid selects a single pseudoknot-free RNA secondary structure that
matches an ensemble mean mountain profile under squared mountain distance.
The mean profile can be computed from base-pairing probabilities (BPPs), so the
ensemble does not need to be enumerated.

The package provides two prediction spaces:

- **Geometry-only Mountain Centroid** searches nonnegative unit-step mountain
  paths. It does not enforce nucleotide pairability or a minimum hairpin
  length.
- **Pairability-constrained Mountain Centroid** permits AU/UA, GC/CG, and
  GU/UG pairs, requires at least three nucleotides inside a hairpin
  (`TURN = 3`), and excludes pseudoknots.

Both dynamic programs return a global minimum of their objective over the
corresponding prediction space for the supplied mean profile. The
pairability-constrained solver has a readable Python implementation and a
faster C++ implementation used by the public command-line interface.

## Requirements

- Python 3.10 or later
- a C++17 compiler and `make`
- Git submodules for the default LinearPartition-V backend

ViennaRNA is installed as a Python dependency. LinearPartition-V is pinned as a
Git submodule and built from source.

## Install from a source checkout

Clone recursively, build the two native components, and install the Python
package:

```sh
git clone --recursive https://github.com/TakumiOtagaki/MountainCentroid.git
cd MountainCentroid
make -C vendor/LinearPartition
make constrained
python -m pip install -e .
```

For development with the locked environment:

```sh
uv sync --frozen --extra test
```

The current packaging configuration assumes a source checkout: the
LinearPartition runner and compiled pairability-constrained solver are not
embedded in a platform wheel.

## Command-line use

LinearPartition-V with beam size 100 is the default BPP backend:

```sh
mountain-centroid --seq GGGAAACCC
```

The beam size and BPP output cutoff can be set explicitly:

```sh
mountain-centroid \
  --seq GGGAAACCC \
  --bpp-beam-size 100 \
  --bpp-cutoff 0.0
```

The optional ViennaRNA partition-function backend can be selected with:

```sh
mountain-centroid --seq GGGAAACCC --bpp-backend vienna
```

Equivalent module invocation:

```sh
python -m mountain_centroid.mountain_pipeline --seq GGGAAACCC
```

## Python use

Compute BPPs and a pairability-constrained prediction:

```python
from mountain_centroid import predict

prediction = predict("GGGAAACCC")
print(prediction.structure)
print(prediction.squared_mountain_error)
```

Use an already computed ensemble mean mountain profile:

```python
from mountain_centroid import predict_from_profile

prediction = predict_from_profile(
    "GGGAAACCC",
    [0.9, 1.8, 2.4, 2.4, 2.4, 2.4, 1.5, 0.5],
)
```

The geometry-only solver is available separately:

```python
from mountain_centroid import relaxed_mountain_centroid

result = relaxed_mountain_centroid(
    [0.9, 1.8, 2.4, 2.4, 2.4, 2.4, 1.5, 0.5]
)
```

Pairs returned by the Python APIs use one-based nucleotide indices.

## BPP backends and reproducibility

LinearPartition uses beam pruning and therefore approximates the partition
function and BPPs. The beam size, BPP cutoff, backend revision, and
thermodynamic settings should be reported with experimental results.
LinearPartition and ViennaRNA outputs should not be pooled without explicitly
accounting for the backend difference.

The BPP calculation is upstream of Mountain Centroid optimization. The global
minimum guarantee applies to the supplied mean profile and prediction space; it
does not remove approximation introduced while computing BPPs.

## Complexity

For an RNA of length `n`:

- the geometry-only solver uses `O(nH)` time and traceback space, where `H` is
  the maximum reachable mountain height (`O(n^2)` worst case);
- the pairability-constrained solver uses `O(D_eff n^3)` time and
  `O(D_eff n^2)` memoization space, where `D_eff` is the number of reached
  external-depth levels (`O(n^4)` time and `O(n^3)` space in the worst case).

## Tests and benchmarks

Run the test suite with:

```sh
uv run --frozen --extra test pytest
```

The tests include small-instance exhaustive enumeration, randomized
Python/C++ parity checks, structural invariants, BPP backend checks, and
metric/formatting tests.

Benchmark the production solver with:

```sh
uv run python benchmarks/benchmark_sequence_constrained.py \
  --lengths 30 50 100 150 200 300 --instances 3
```

## Repository layout

```text
src/mountain_centroid/   Python implementation and public API
cpp/                     C++ pairability-constrained solver
tests/                   correctness, parity, backend, and formatting tests
benchmarks/              solver and pipeline benchmarks
vendor/LinearPartition/  pinned default BPP backend (Git submodule)
```

## Provenance

The first implementation was extracted from
`TakumiOtagaki/ZukerStyleCentroidFold` at commit
`17d39d288ed6900561d93de47a8a1b2b98d2c329`. Manuscript analyses, datasets,
and generated figures are maintained separately from this reusable package.

## License

Code and documentation authored for MountainCentroid are available under the
[BSD 3-Clause License](LICENSE). LinearPartition and ViennaRNA are subject to
their respective license terms; see [Third-party software](THIRD_PARTY_LICENSES.md).

## Citation

Please use the citation metadata in [`CITATION.cff`](CITATION.cff). A citation
to the accompanying paper will be added after publication.
