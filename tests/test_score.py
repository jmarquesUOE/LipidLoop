"""Tests for the LipiDex-equivalent scoring.

These pin the behaviour that a naive reimplementation would get wrong, so a future
refactor cannot silently drift away from LipiDex.
"""
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.score import dot_product, ppm_diff          # noqa: E402
from lipidloop.filters import Thresholds                    # noqa: E402


def test_identical_spectra_score_1000():
    mz = np.array([239.24, 268.27, 300.29, 337.35])
    it = np.array([508.0, 500.0, 500.0, 999.0])
    assert dot_product(mz, it, mz, it) == pytest.approx(1000.0, rel=1e-9)


def test_disjoint_spectra_score_zero():
    assert dot_product([100.0, 200.0], [999.0, 500.0],
                       [300.0, 400.0], [999.0, 500.0]) == 0.0


def test_reverse_ignores_extra_sample_peaks():
    """A co-isolated contaminant adds sample peaks the library does not predict. The reverse
    score should be unaffected; the forward score should drop."""
    lib_mz, lib_int = [100.0, 200.0], [999.0, 500.0]
    clean = dot_product(lib_mz, lib_int, lib_mz, lib_int, reverse=True)
    dirty = dot_product([100.0, 150.0, 200.0], [999.0, 800.0, 500.0],
                        lib_mz, lib_int, reverse=True)
    assert dirty == pytest.approx(clean, rel=1e-9)

    fwd_clean = dot_product(lib_mz, lib_int, lib_mz, lib_int)
    fwd_dirty = dot_product([100.0, 150.0, 200.0], [999.0, 800.0, 500.0], lib_mz, lib_int)
    assert fwd_dirty < fwd_clean


def test_tolerance_controls_matching():
    lib_mz, lib_int = [100.000], [999.0]
    assert dot_product([100.005], [999.0], lib_mz, lib_int, mz_tol=0.01) > 0
    assert dot_product([100.050], [999.0], lib_mz, lib_int, mz_tol=0.01) == 0.0


def test_ppm_diff_sign_and_magnitude():
    assert ppm_diff(666.6467, 666.64) == pytest.approx(10.05, abs=0.01)
    assert ppm_diff(666.6333, 666.64) < 0


def test_thresholds_reverse_is_stricter():
    t = Thresholds()
    assert t.min_rev_dot_product > t.min_dot_product
    assert t.passes(ppm=5.0, dot=600, rev_dot=750, purity=80)
    assert not t.passes(ppm=5.0, dot=600, rev_dot=650, purity=80)   # rev below 700
    assert not t.passes(ppm=25.0, dot=600, rev_dot=750, purity=80)  # ppm over 20
    assert not t.passes(ppm=5.0, dot=600, rev_dot=750, purity=70)   # purity under 75
