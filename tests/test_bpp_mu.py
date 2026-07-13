from pathlib import Path

import numpy as np
import pytest

from mountain_centroid.bpp_mu import (
    compute_bpp_linearpartition,
    mountain_expectation_from_bpp,
)


def test_mountain_expectation_uses_cut_convention():
    bpp = np.zeros((4, 4), dtype=float)
    bpp[0, 2] = 0.2
    bpp[1, 3] = 0.3

    assert mountain_expectation_from_bpp(bpp) == pytest.approx([0.2, 0.5, 0.3])


def test_linearpartition_backend_smoke():
    runner = Path(__file__).parents[1] / "vendor" / "LinearPartition" / "linearpartition"
    if not runner.is_file() or not (runner.parent / "bin" / "linearpartition_v").is_file():
        pytest.skip("LinearPartition submodule has not been built")

    bpp = compute_bpp_linearpartition("GGGAAACCC", executable=runner, beam_size=100)

    assert bpp.shape == (9, 9)
    assert np.all((0.0 <= bpp) & (bpp <= 1.0))
    assert bpp[0, 8] > 0.5

