"""Gap filling by re-integration.

The unit under test is not "does it put a number in the hole" — it is whether that number is on
the same scale as the rest of the column, and whether a genuine absence survives. A filled value
on an unknown scale is worse than a gap, because a gap announces itself.
"""
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.gapfill import (GapFilling, _integrate, _peak, _xic,   # noqa: E402
                                count_absences_with_signal, fill_gaps)
from lipidloop.peaks import CompoundGroup, Sample   # noqa: E402

FWHM = 0.09
QUANT = 760.5850


@dataclass
class FakeGroup:
    """Enough of a CompoundGroup for the filler: areas, a quant ion and a retention time."""
    areas: list
    quant_ion: float = QUANT
    retention: float = 5.0
    keep: bool = True


def scans(peak_rt=5.0, height=1000.0, noise=1.0, mz=QUANT, n=61, start=4.7, step=0.01):
    """A gaussian peak on a flat noise floor, as (rt, mz array, intensity array) per scan."""
    out = []
    for i in range(n):
        rt = start + i * step
        signal = height * np.exp(-0.5 * ((rt - peak_rt) / (FWHM / 2.355)) ** 2) + noise
        out.append((rt, np.array([mz - 0.5, mz, mz + 0.5]),
                    np.array([noise, signal, noise])))
    return out


def test_xic_picks_only_the_target_mass():
    s = scans()
    times = np.array([x[0] for x in s])
    axis, signal = _xic(s, times, QUANT, 10.0, 4.8, 5.2)
    assert axis.size > 0
    assert signal.max() > 900          # the peak, not the 1.0 neighbours
    _, off = _xic(s, times, QUANT + 5.0, 10.0, 4.8, 5.2)
    assert off.sum() == 0              # 5 Da away: nothing


def test_a_peak_is_recognised_and_an_empty_window_is_not():
    s = scans()
    times = np.array([x[0] for x in s])
    axis, signal = _xic(s, times, QUANT, 10.0, 4.4, 5.6)
    apex, ratio = _peak(axis, signal, 1.5, 5.0, 0.2)
    assert abs(apex - 5.0) < 0.02 and ratio > 10
    assert _peak(np.array([1.0]), np.array([0.0]), 1.5, 1.0, 0.2) is None


def test_a_real_peak_scores_far_above_noise():
    """What can honestly be asserted here is *separation*, not an absolute threshold.

    The module already got this wrong once by asserting one: apex over the median of the lower
    half of the window called every noise window a peak, filled 100% of one polarity's 7,297 gaps
    and then reported that every whole-group absence carried signal. The threshold that separates
    signal from noise depends on the noise, so it is measured on the data
    (`scripts/gapfill_null.py`) rather than assumed here. This pins the property that makes such
    a measurement possible at all.
    """
    rng = np.random.default_rng(0)
    axis = np.linspace(4.4, 5.6, 121)
    noise_scores = []
    for _ in range(40):
        signal = rng.gamma(shape=2.0, scale=50.0, size=axis.size)
        found = _peak(axis, signal, 0.0, 5.0, 0.2)
        noise_scores.append(found[1] if found else 0.0)

    s = scans(height=5000.0, noise=100.0)
    times = np.array([x[0] for x in s])
    a, sig = _xic(s, times, QUANT, 10.0, 4.4, 5.6)
    peak = _peak(a, sig, 0.0, 5.0, 0.2)
    assert peak is not None
    assert peak[1] > 10 * float(np.median(noise_scores))


def _run(groups, sample_areas_column=0, roles=("sample",), **kw):
    samples = [Sample(file=f"S{i}", role=r) for i, r in enumerate(roles)]
    mzml = {s.file: "unused" for s in samples}
    params = GapFilling(min_calibration_features=2, **kw)
    return samples, mzml, params


def test_the_filled_value_is_calibrated_onto_the_table_scale(monkeypatch):
    """The hazard this module exists for. Raw summed intensity is not pyOpenMS's area; if the
    filler ignores that, the column silently mixes two units. Here the table's areas are 1000x
    the raw sum, and the filled value must come out on the table's scale, not the raw one."""
    import lipidloop.gapfill as gf
    monkeypatch.setattr(gf, "_ms1", lambda path: scans())

    known = [FakeGroup(areas=[0.0, 0.0]) for _ in range(4)]
    raw_sum = None
    for g in known:
        g.areas = [1.0, 1.0]
    # measure what the integrator gives, then set the table areas to 1000x that
    s = scans(); times = np.array([x[0] for x in s])
    axis, signal = _xic(s, times, QUANT, 10.0, 5.0 - 0.2, 5.0 + 0.2)
    raw_sum = _integrate(axis, signal, 5.0, FWHM)
    for g in known:
        g.areas = [raw_sum * 1000.0, raw_sum * 1000.0]
    gap = FakeGroup(areas=[raw_sum * 1000.0, 0.0])

    samples, mzml, params = _run(None, roles=("sample", "sample"))
    report = fill_gaps(known + [gap], samples, {s.file: "x" for s in samples}, params, FWHM)
    assert report.gaps == 1 and report.filled_peak == 1
    assert gap.areas[1] == pytest.approx(raw_sum * 1000.0, rel=0.05)


def test_an_unstable_calibration_refuses_the_injection_rather_than_guessing(monkeypatch):
    """If the ratio between reported and integrated area is all over the place, the relationship
    is not understood and a filled number would be a fabrication. Leave the gap."""
    import lipidloop.gapfill as gf
    monkeypatch.setattr(gf, "_ms1", lambda path: scans())
    # Varying in column 1 — the column being filled, which is the one calibration reads.
    known = [FakeGroup(areas=[1.0, 10.0 ** n]) for n in range(2, 8)]
    gap = FakeGroup(areas=[100.0, 0.0])
    samples, mzml, params = _run(None, roles=("sample", "sample"))
    report = fill_gaps(known + [gap], samples, {s.file: "x" for s in samples}, params, FWHM)
    assert report.refused and gap.areas[1] == 0.0


def test_nothing_there_stays_zero(monkeypatch):
    """A real absence must survive re-integration. This is what makes the presence/absence call
    trustworthy — the filler has to be able to return 'still nothing'."""
    import lipidloop.gapfill as gf
    monkeypatch.setattr(gf, "_ms1", lambda path: scans(height=0.0, noise=0.0))
    known = [FakeGroup(areas=[5.0, 5.0]) for _ in range(4)]
    gap = FakeGroup(areas=[5.0, 0.0])
    samples, mzml, params = _run(None, roles=("sample", "sample"))
    report = fill_gaps(known + [gap], samples, {s.file: "x" for s in samples}, params, FWHM)
    assert gap.areas[1] == 0.0
    assert report.filled_peak == 0


def test_blanks_are_not_filled(monkeypatch):
    """A zero in a blank is information about the blank. Filling it would erase the evidence the
    blank filter runs on."""
    import lipidloop.gapfill as gf
    monkeypatch.setattr(gf, "_ms1", lambda path: scans())
    known = [FakeGroup(areas=[100.0, 100.0]) for _ in range(4)]
    gap = FakeGroup(areas=[100.0, 0.0])
    samples = [Sample(file="S0", role="sample"), Sample(file="B1", role="blank")]
    report = fill_gaps(known + [gap], samples, {s.file: "x" for s in samples},
                       GapFilling(min_calibration_features=2), FWHM)
    assert gap.areas[1] == 0.0 and report.gaps == 0


def test_the_absence_check_counts_groups_that_gained_signal():
    """The number that decides whether the presence filter's whole-group protection stands."""
    cells = {("A",): [0, 1], ("B",): [2, 3]}
    stayed = CompoundGroup(name="a", mw=1.0, retention=1.0, max_area=1.0, has_ms2=True,
                           areas=[5.0, 5.0, 0.0, 0.0])
    gained = CompoundGroup(name="b", mw=1.0, retention=1.0, max_area=1.0, has_ms2=True,
                           areas=[5.0, 5.0, 3.0, 0.0])
    before = {id(stayed): [5.0, 5.0, 0.0, 0.0], id(gained): [5.0, 5.0, 0.0, 0.0]}
    with_signal, tested = count_absences_with_signal([stayed, gained], [], cells, before)
    assert tested == 2 and with_signal == 1
