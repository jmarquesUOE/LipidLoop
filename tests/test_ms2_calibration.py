"""The MS2 axis measured from calibrant ions: offset, capture, window, and the keep/derive rule."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import pytest

from lipidloop.ms2_calibration import (KEEP_CAPTURE, LADDER, PLAUSIBLE, Ms2Calibration, check_ms2_tolerance,
                                        collect_errors, derive_decision, summarise)


class _Settings:
    def __init__(self, polarity):
        self._p = polarity

    def getPolarity(self):
        return 2 if self._p == "-" else 1


class _Spectrum:
    """The four calls `collect_errors` makes on a pyOpenMS spectrum."""

    def __init__(self, mz, intensity, polarity="+", level=2, filter_string=None):
        self._mz = np.asarray(mz, float); self._it = np.asarray(intensity, float)
        self._pol = polarity; self._level = level; self._fs = filter_string

    def getMSLevel(self):
        return self._level

    def get_peaks(self):
        return self._mz, self._it

    def getInstrumentSettings(self):
        return _Settings(self._pol)

    def metaValueExists(self, key):
        return self._fs is not None and key == "filter string"

    def getMetaValue(self, key):
        return self._fs


def _synthetic(n, offset, sd, polarity="+", seed=0, filter_string=None):
    rng = np.random.default_rng(seed)
    cal = 184.0733 if polarity == "+" else 281.2486
    spectra = []
    for _ in range(n):
        mz = np.append(rng.uniform(100, 900, 20), cal + offset + rng.normal(0, sd))
        it = np.append(rng.uniform(1, 50, 20), 1000.0)
        order = np.argsort(mz)
        spectra.append(_Spectrum(mz[order], it[order], polarity, filter_string=filter_string))
    return spectra


def test_orbitrap_keeps_the_inherited_window():
    errors, n = collect_errors(_synthetic(400, 0.0003, 0.002, filter_string="FTMS + p ESI d Full ms2"), "orbitrap")
    cal = summarise(errors, n)
    assert cal.analyzer == "orbitrap" and cal.n_peaks >= 370
    assert abs(cal.offset_da - 0.0003) < 0.001
    assert cal.capture[0.01] > 0.95 and cal.window in (0.005, 0.01)
    window, msg, plausible = check_ms2_tolerance(0.01, cal)
    assert plausible and "kept" in msg


def test_ion_trap_is_measured_with_its_offset():
    errors, n = collect_errors(_synthetic(600, -0.03, 0.12, filter_string="ITMS + c ESI d Full ms2"), "orbitrap")
    cal = summarise(errors, n)
    assert cal.analyzer == "ion trap"
    assert abs(cal.offset_da + 0.03) < 0.01
    assert cal.capture[0.01] < 0.2 and cal.window in (0.3, 0.5)
    window, msg, plausible = check_ms2_tolerance(0.01, cal)
    assert not plausible and window == cal.window and "OUTSIDE" in msg


def test_a_tof_that_captures_too_little_is_replaced_even_inside_the_band():
    # 88 % within 10 mDa and the rest at 60-100 mDa, the timsTOF pattern
    good = _synthetic(440, 0.0, 0.003, seed=1)
    tail = _synthetic(60, 0.08, 0.01, seed=2)
    errors, n = collect_errors(good + tail, "tof")
    cal = summarise(errors, n)
    assert cal.analyzer == "tof"
    assert cal.capture[0.01] < KEEP_CAPTURE
    window, msg, plausible = check_ms2_tolerance(0.01, cal)
    assert not plausible and window == 0.1 and "captures only" in msg


def test_window_is_clamped_to_the_analyzer_band():
    cal = Ms2Calibration(analyzer="ion trap", n_peaks=500, capture={w: 1.0 for w in LADDER})
    cal.window = max(LADDER[0], PLAUSIBLE["ion trap"][0])
    window, msg, plausible = check_ms2_tolerance(0.5, cal)
    assert plausible and window >= PLAUSIBLE["ion trap"][0]


def test_too_few_calibrant_peaks_changes_nothing():
    errors, n = collect_errors(_synthetic(20, -0.03, 0.12, filter_string="ITMS + c ESI d Full ms2"), "orbitrap")
    cal = summarise(errors, n)
    assert not cal.usable
    window, msg, plausible = check_ms2_tolerance(0.01, cal)
    assert plausible and window == 0.01 and "too few" in msg


def test_negative_mode_uses_the_carboxylate_calibrants():
    errors, n = collect_errors(_synthetic(300, 0.0, 0.002, polarity="-"), "orbitrap")
    cal = summarise(errors, n)
    assert cal.polarity == "-" and cal.n_peaks >= 280 and cal.window in (0.005, 0.01)


# ── the window and the offset are separate decisions ────────────────────────────────────────
# They were once set in one statement, so pinning `ms2_tol` on the target-decoy knee also threw
# away the measured fragment offset — silently, since the measurement still ran and was still
# written to run_config.json with `offset_applied_da: 0.0`. On an ion trap reading 25-47 mDa low
# that is the correction you most want to keep.

def test_offset_mode_keeps_the_configured_window_and_still_corrects_the_drift():
    assert derive_decision("offset", usable=True, plausible=False) == (False, True)
    assert derive_decision("offset", usable=True, plausible=True) == (False, True)


def test_false_applies_neither_which_is_the_documented_report_only_behaviour():
    assert derive_decision(False, usable=True, plausible=False) == (False, False)


def test_true_applies_both():
    assert derive_decision(True, usable=True, plausible=True) == (True, True)


def test_auto_replaces_an_implausible_window_and_takes_the_offset_with_it():
    assert derive_decision("auto", usable=True, plausible=False) == (True, True)
    assert derive_decision("auto", usable=True, plausible=True) == (False, False)


def test_nothing_is_applied_when_the_measurement_is_unusable():
    """Too few calibrant peaks: the configured window stands and no offset is invented."""
    for mode in (True, False, "auto", "offset"):
        assert derive_decision(mode, usable=False, plausible=False) == (False, False)


def test_an_unrecognised_string_raises_instead_of_meaning_true():
    """It used to be truthy, so a typo silently meant `always derive`."""
    with pytest.raises(ValueError, match="expected true, false"):
        derive_decision("yes", usable=True, plausible=True)
