"""Tests for the targeted feature-recovery pass.

The contract that matters: it only ever promotes a peak that an identification already stands
behind, and never one a feature already represents. Those two together are what stop it from
manufacturing features out of noise.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.features import FeatureParams, recover_targeted_features   # noqa: E402

pyopenms = pytest.importorskip("pyopenms")


class Trace:
    """Stands in for a Kernel_MassTrace — only these four methods are used."""

    def __init__(self, mz, rt_seconds, area=1000.0, fwhm=4.0):
        self._mz, self._rt, self._area, self._fwhm = mz, rt_seconds, area, fwhm

    def getCentroidMZ(self): return self._mz
    def getCentroidRT(self): return self._rt
    def computePeakArea(self): return self._area
    def getFWHM(self): return self._fwhm


def feature_map(features=()):
    fm = pyopenms.FeatureMap()
    for mz, rt in features:
        f = pyopenms.Feature()
        f.setMZ(mz)
        f.setRT(rt)
        fm.push_back(f)
    return fm


PARAMS = FeatureParams()


def test_promotes_an_identified_peak_with_no_feature():
    fm = feature_map()
    n = recover_targeted_features([Trace(803.6950, 658.2)], fm, [(803.6950, 658.2)], PARAMS)
    assert n == 1 and fm.size() == 1
    assert fm[0].getMZ() == pytest.approx(803.6950)
    assert fm[0].getMetaValue(b"recovered") == "true"


def test_never_promotes_a_peak_without_an_identification():
    """The whole safeguard: no identification, no feature, however good the peak looks."""
    fm = feature_map()
    assert recover_targeted_features([Trace(803.6950, 658.2, area=1e9)], fm, [], PARAMS) == 0
    assert recover_targeted_features([Trace(803.6950, 658.2)], fm,
                                     [(999.9999, 658.2)], PARAMS) == 0
    assert fm.size() == 0


def test_never_duplicates_a_feature_that_already_exists():
    fm = feature_map([(803.6950, 658.2)])
    assert recover_targeted_features([Trace(803.6950, 658.2)], fm,
                                     [(803.6950, 658.2)], PARAMS) == 0
    assert fm.size() == 1


def test_respects_the_mass_and_time_windows():
    fm = feature_map()
    # 100 ppm away in mass, and 60 s away in time — both outside
    assert recover_targeted_features([Trace(803.7754, 658.2)], fm,
                                     [(803.6950, 658.2)], PARAMS) == 0
    assert recover_targeted_features([Trace(803.6950, 718.2)], fm,
                                     [(803.6950, 658.2)], PARAMS) == 0
    # inside both
    assert recover_targeted_features([Trace(803.6990, 663.0)], fm,
                                     [(803.6950, 658.2)], PARAMS) == 1


def test_carries_the_peak_area_and_width_across():
    fm = feature_map()
    recover_targeted_features([Trace(700.0, 600.0, area=54321.0, fwhm=4.5)], fm,
                              [(700.0, 600.0)], PARAMS)
    assert fm[0].getIntensity() == pytest.approx(54321.0)
    assert fm[0].getWidth() == pytest.approx(4.5)


def test_the_isotope_collision_it_exists_for():
    """SM d41:0 sits within 11 ppm of the M+2 isotope of SM d41:1, so assembly eats it.

    The feature at 801.6915 exists; the peak at 803.6912 has been consumed as its isotope and
    has no feature of its own. An identification at that mass and time brings it back.
    """
    fm = feature_map([(801.6915, 658.8)])
    n = recover_targeted_features([Trace(801.6915, 658.8), Trace(803.6912, 658.8)],
                                  fm, [(803.6950, 658.2)], PARAMS)
    assert n == 1
    assert fm.size() == 2
    assert any(abs(f.getMZ() - 803.6912) < 1e-6 for f in fm)
