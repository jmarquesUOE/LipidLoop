"""Filling gaps by going back to the raw file, the way Compound Discoverer does it.

A zero in the table means the feature finder did not detect a peak in that injection. It does not
mean nothing is there. Statistical imputation — half the feature minimum, k-nearest neighbours,
whatever — answers the question by assuming the answer. Re-integration asks the instrument.

Compound Discoverer's `Fill Gaps` node, read out of the facility's own workflows
(`CD_Workflows/20230531_LipiDex_ALIGNED_Pos.cdProcessingWF`, the 25 min Kinetex method these
datasets were acquired with):

    Mass Tolerance            10 ppm
    S/N Threshold             1.5
    Use Real Peak Detection   True

Those are the defaults here, for the obvious reason that a replacement for a pipeline should not
quietly disagree with it about what a filled value means.

## Why this is not just tidier than imputation

It **tests** an assumption the presence filter rests on. `presence.py` protects a group with no
detections at all, on the reasoning that a whole-group zero is a real absence and filling it would
turn a presence/absence result into a weak fold change. That reasoning is sound *if* the zeros are
real. Re-integration is how you find out: if a supposedly absent group turns out to carry a peak
above the S/N threshold, the absence was a detection failure and the on/off call was wrong. The
report says how often that happens, because it is the number that decides whether those calls
survive.

## The scale problem, which is the reason this file is careful

Areas in the table come from pyOpenMS (`FeatureFindingMetabo`, `getIntensity()`). An area computed
here by a different summation is **not guaranteed to be in the same units**, and a table mixing two
scales is worse than a table with gaps — the gaps at least announce themselves.

So the integration is calibrated per injection, against that injection's own detected features:
integrate the raw signal for features whose area is already known, take the median ratio of the
reported area to the integrated one, and apply it. The ratio is reported. If it is unstable
(interquartile spread wide relative to the median) the filling is refused for that injection and
the gaps are left as gaps, because a filled value on an unknown scale is a fabrication.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from statistics import median

import numpy as np


@dataclass
class GapFilling:
    """Re-integrate the raw signal where the feature finder found nothing.

    Defaults are Compound Discoverer's, taken from the facility workflow for this method.

    `rt_window` is how far either side of the group's retention time to look, in minutes. It has
    to absorb the alignment shift — measured at about 5.5 s spread on these methods, and CD's own
    aligner allows 0.5 min — without reaching into a neighbouring peak.

    `max_calibration_spread` refuses the whole injection when the area calibration is too unstable
    to trust. See the module docstring.
    """

    enabled: bool = True
    mass_tolerance_ppm: float = 10.0
    sn_threshold: float = 1.5
    rt_window: float = 0.2
    use_real_peak_detection: bool = True
    max_calibration_spread: float = 0.5
    # Minimum detected features in an injection before its calibration is believed.
    min_calibration_features: int = 50


@dataclass
class GapFillReport:
    gaps: int = 0
    filled_peak: int = 0        # a real peak, clearing S/N
    filled_area: int = 0        # signal present but below S/N — an upper bound, still a number
    empty: int = 0              # nothing there; the zero stands
    refused: list = field(default_factory=list)      # injections whose calibration was unstable
    calibration: dict = field(default_factory=dict)  # injection -> median area ratio
    absences_with_signal: int = 0                    # whole-group zeros that turned out to have a peak
    absences_tested: int = 0

    def __str__(self) -> str:
        out = (f"gap filling: {self.gaps} gaps — {self.filled_peak} filled from a detected peak, "
               f"{self.filled_area} from integrated signal below S/N, {self.empty} left empty")
        if self.calibration:
            out += f"; median area calibration {median(self.calibration.values()):.2f}x"
        if self.refused:
            out += f"; REFUSED {len(self.refused)} injection(s) on unstable calibration"
        if self.absences_tested:
            out += (f"; {self.absences_with_signal}/{self.absences_tested} whole-group absences "
                    f"turned out to carry a peak")
        return out


def _ms1(mzml: str | Path):
    """(retention in minutes, m/z array, intensity array) per MS1 scan, in time order."""
    from pyopenms import MSExperiment, MzMLFile
    experiment = MSExperiment()
    MzMLFile().load(str(mzml), experiment)
    scans = []
    for spectrum in experiment:
        if spectrum.getMSLevel() != 1:
            continue
        mz, intensity = spectrum.get_peaks()
        scans.append((spectrum.getRT() / 60.0, np.asarray(mz), np.asarray(intensity)))
    scans.sort(key=lambda s: s[0])
    return scans


def _xic(scans, times: np.ndarray, target: float, tolerance_ppm: float,
         low: float, high: float) -> tuple[np.ndarray, np.ndarray]:
    """Extracted ion chromatogram over [low, high]: the retention axis and summed intensity."""
    first, last = int(np.searchsorted(times, low)), int(np.searchsorted(times, high))
    if last <= first:
        return np.empty(0), np.empty(0)
    half = target * tolerance_ppm * 1e-6
    axis = np.empty(last - first)
    signal = np.empty(last - first)
    for n, index in enumerate(range(first, last)):
        retention, mz, intensity = scans[index]
        left, right = np.searchsorted(mz, (target - half, target + half))
        axis[n] = retention
        signal[n] = intensity[left:right].sum() if right > left else 0.0
    return axis, signal


def _integrate(axis: np.ndarray, signal: np.ndarray, centre: float, halfwidth: float) -> float:
    """Summed intensity within `halfwidth` of `centre`.

    Deliberately the same procedure used to calibrate against known areas — the calibration is
    only meaningful if the thing being calibrated is measured identically.
    """
    if axis.size == 0:
        return 0.0
    inside = np.abs(axis - centre) <= halfwidth
    return float(signal[inside].sum())


def _peak(axis: np.ndarray, signal: np.ndarray, sn_threshold: float,
          centre: float, halfwidth: float, min_scans: int = 3):
    """(apex retention, signal-to-noise) for a peak at the expected time, or None.

    ⚠ The obvious version of this — apex over the median of the lower half of the window — does
    not work, and failed silently rather than loudly. For a window of pure noise that ratio sits
    around 2-5, so an S/N threshold of 1.5 passes essentially always: on the brain set it filled
    100% of 7,297 negative-mode gaps from a "detected peak" and left none empty, then reported
    that every whole-group absence carried signal. A test that never says no cannot be evidence
    that something is there.

    So: noise is estimated from the **flanks**, outside the peak region, as a robust sigma; the
    apex is taken only within the peak region; and it must be part of a run of at least
    `min_scans` consecutive points above the noise, because a single high scan is a spike and a
    peak has width. That is what `Use Real Peak Detection` means.
    """
    if signal.size == 0 or signal.max() <= 0:
        return None
    core = np.abs(axis - centre) <= halfwidth
    flank = ~core
    if core.sum() == 0 or flank.sum() < 3:
        return None
    baseline = float(np.median(signal[flank]))
    sigma = float(np.median(np.abs(signal[flank] - baseline))) * 1.4826
    if sigma <= 0:
        sigma = max(baseline, 1.0) * 0.1
    inside = np.where(core)[0]
    apex = inside[int(np.argmax(signal[inside]))]
    height = signal[apex] - baseline
    if height <= 0:
        return None
    above = signal > baseline + sigma
    run = 1
    for step in (-1, 1):
        index = apex + step
        while 0 <= index < signal.size and above[index]:
            run += 1
            index += step
    if run < min_scans:
        return None
    return float(axis[apex]), float(height / sigma)


def fill_gaps(groups, samples, mzml_by_sample: dict, params: GapFilling,
              avg_fwhm: float, log=None) -> GapFillReport:
    """Fill zero areas by re-integrating the raw signal. Mutates `group.areas`.

    Only sample and QC injections are filled; a blank's zero is information about the blank.
    """
    say = log or (lambda _: None)
    report = GapFillReport()
    halfwidth = max(avg_fwhm, 0.02)
    targets = [(n, s) for n, s in enumerate(samples) if s.role != "blank"]

    for column, sample in targets:
        path = mzml_by_sample.get(sample.file)
        gaps = [g for g in groups if g.keep and g.quant_ion and g.areas[column] <= 0]
        if path is None or not gaps:
            continue
        scans = _ms1(path)
        if not scans:
            continue
        times = np.array([s[0] for s in scans])

        # ── calibrate against this injection's own detected features
        known = [g for g in groups if g.keep and g.quant_ion and g.areas[column] > 0]
        ratios = []
        for group in known[:: max(1, len(known) // 400)]:
            axis, signal = _xic(scans, times, group.quant_ion, params.mass_tolerance_ppm,
                                group.retention - params.rt_window,
                                group.retention + params.rt_window)
            raw = _integrate(axis, signal, group.retention, halfwidth)
            if raw > 0:
                ratios.append(group.areas[column] / raw)
        if len(ratios) < params.min_calibration_features:
            report.refused.append(sample.file)
            say(f"  gap filling refused for {sample.file}: only {len(ratios)} usable "
                f"calibration features")
            continue
        ratios.sort()
        centre = median(ratios)
        spread = (ratios[int(0.75 * (len(ratios) - 1))] - ratios[int(0.25 * (len(ratios) - 1))])
        if centre <= 0 or spread / centre > params.max_calibration_spread:
            report.refused.append(sample.file)
            say(f"  gap filling refused for {sample.file}: area calibration unstable "
                f"(median {centre:.2f}, IQR {spread:.2f}) — gaps left as gaps")
            continue
        report.calibration[sample.file] = centre

        # ── fill
        for group in gaps:
            report.gaps += 1
            axis, signal = _xic(scans, times, group.quant_ion, params.mass_tolerance_ppm,
                                group.retention - params.rt_window * 3,
                                group.retention + params.rt_window * 3)
            found = _peak(axis, signal, params.sn_threshold, group.retention, params.rt_window)
            if found is None:
                report.empty += 1
                continue
            apex, ratio = found
            real = params.use_real_peak_detection and ratio >= params.sn_threshold
            at = apex if real else group.retention
            area = _integrate(axis, signal, at, halfwidth) * centre
            if area <= 0:
                report.empty += 1
                continue
            group.areas[column] = area
            if real:
                report.filled_peak += 1
            else:
                report.filled_area += 1

    say(str(report))
    return report


def count_absences_with_signal(groups, samples, cells: dict, before: dict) -> tuple[int, int]:
    """How many whole-group absences turned out to carry signal after filling.

    `before` maps id(group) -> the areas as they were. The presence filter's protection of a
    whole-group zero is only correct while those zeros are real; this is the check on it.
    """
    tested = with_signal = 0
    for group in groups:
        if not group.keep:
            continue
        original = before.get(id(group))
        if original is None:
            continue
        for members in cells.values():
            if any(original[i] > 0 for i in members):
                continue
            tested += 1
            if any(group.areas[i] > 0 for i in members):
                with_signal += 1
    return with_signal, tested
