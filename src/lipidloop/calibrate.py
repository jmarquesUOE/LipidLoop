"""Deriving the parameters that describe a run, from that run.

The proposal: run once, measure the data, adjust the thresholds to it, run again. The boundary
that makes it safe, and the reason this module refuses to be general:

    Iterate on parameters that DESCRIBE the measurement.
    Never on parameters that DECIDE what counts as a hit.

Peak width, mass calibration, retention spread and the noise floor are properties of the run and
can be measured. Dot product 500 / reverse 700, fatty-acid purity 75, the blank multiplier and the
minimum feature count are decision rules. Tuning those on the data makes the result unfalsifiable —
"we chose the threshold that gave the most identifications" is not a method — and destroys
comparability with LipiDex, which is the property the whole reimplementation rests on.

That boundary is enforced here rather than documented: `FORBIDDEN` is checked on every derivation
and raises. A convention in a docstring survives exactly as long as the person who wrote it.

## Why the noise floor matters more than the rest

`noise_threshold` is an **absolute intensity**. Intensity scales differ between instruments, so
5,000 on an Orbitrap Fusion Lumos is not 5,000 on a Q Exactive Plus. For a facility pipeline
running on several instruments this is the most likely thing to be silently wrong on the next new
one, and unlike the others it fails quietly: too high and features vanish, too low and mass trace
detection explodes and the run never finishes.

It is derived as a multiple of the file's own measured noise floor, with the multiplier anchored so
that on the instrument the pipeline was validated against it reproduces the validated value. A
derivation that changes the answer on the data it was tuned on is not a calibration, it is a
different pipeline.

## Everything derived is written down

Both the configured and the derived value, with the pass it came from, go into `run_config.json`.
Otherwise "run twice" becomes "run until it looks good", and a result stops being reproducible from
its configuration.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
from pathlib import Path

import numpy as np

# Decision rules. Never derived from the data, whatever a caller asks for.
FORBIDDEN = frozenset({
    "min_dot_product", "min_reverse_dot_product", "min_purity", "dot_product", "reverse_dot",
    "blank_filter.multiplier", "min_feature_count", "min_correlation", "min_cv_improvement",
    "sn_threshold", "min_fraction", "rt_filter_multiplier",
})

# Anchored on the validated runs, by measurement rather than by guess. The noise floor on the
# Lumos brain files is 920 (positive) and 866 (negative) counts, against the `noise_threshold` of
# 5,000 that was fitted against the reference results — so the multiplier that reproduces the
# validated setting on the validated instrument is 5.44 and 5.77. A derivation that changes the
# answer on the data it was tuned against is not a calibration, it is a different pipeline.
NOISE_MULTIPLIER = 5.5

# Nothing derived may leave these ranges. A pathological file should degrade to the configured
# value, not silently reconfigure the run.
BOUNDS = {
    "noise_threshold": (500.0, 100_000.0),
    "chrom_fwhm": (1.0, 120.0),
    "min_fwhm": (0.2, 30.0),
    "max_fwhm": (5.0, 300.0),
    "align_rt_max_difference": (5.0, 300.0),
    "mass_offset_ppm": (-25.0, 25.0),
    "retention_model_max_error": (0.1, 5.0),
}


@dataclass
class Calibration:
    """Which descriptive parameters to derive from a first pass.

    Off by default. A two-pass run costs a full second pass, and a facility pipeline should not
    silently do twice the work — nor silently change its own settings — without being asked.
    """

    enabled: bool = False
    mass_offset: bool = True
    peak_width: bool = True
    alignment: bool = True
    noise_floor: bool = True
    retention_model: bool = True


@dataclass
class Derived:
    """What was measured, what it was, and what it becomes. Written into `run_config.json`."""

    values: dict = field(default_factory=dict)      # parameter -> derived value
    measured: dict = field(default_factory=dict)    # parameter -> what the data showed
    previous: dict = field(default_factory=dict)    # parameter -> what the config said
    clamped: list = field(default_factory=list)
    skipped: dict = field(default_factory=dict)     # parameter -> why nothing was derived

    def __str__(self) -> str:
        if not self.values:
            return "calibration: nothing derived"
        parts = [f"{k} {self.previous.get(k)} -> {v:.4g}" for k, v in sorted(self.values.items())]
        out = "calibration: " + "; ".join(parts)
        if self.clamped:
            out += f" (clamped: {', '.join(self.clamped)})"
        if self.skipped:
            out += f" (skipped: {', '.join(f'{k} — {v}' for k, v in self.skipped.items())})"
        return out


def _guard(name: str) -> None:
    if name in FORBIDDEN:
        raise ValueError(
            f"{name!r} is a decision rule, not a description of the measurement. Deriving it from "
            f"the data makes the result unfalsifiable and breaks comparability with LipiDex. "
            f"See calibrate.FORBIDDEN.")


def _bounded(name: str, value: float, report: Derived) -> float:
    low, high = BOUNDS.get(name, (-np.inf, np.inf))
    if value < low or value > high:
        report.clamped.append(name)
        return float(min(max(value, low), high))
    return float(value)


# ── estimators, each a pure function of what the run measured ───────────────────────────────

def mass_offset_ppm(deltas) -> float | None:
    """The systematic mass error, as the median of the per-match ppm deltas.

    A symmetric +/-10 ppm search window centred on zero is not symmetric around where the ions
    actually are: at a median of -1.78 ppm it is really -8.2/+11.8. Correcting the centre buys
    tolerance on both sides without widening the window, which would cost specificity.
    """
    values = np.asarray([d for d in deltas if d == d], dtype=float)
    return float(np.median(values)) if values.size >= 50 else None


def peak_width_seconds(widths) -> tuple[float, float, float] | None:
    """(chrom_fwhm, min_fwhm, max_fwhm) in seconds, from the measured width distribution.

    Robust percentiles rather than the extremes: one badly integrated feature should not set the
    bound that decides what counts as a peak for the whole run.
    """
    values = np.asarray([w for w in widths if w and w == w], dtype=float)
    values = values[values > 0]
    if values.size < 50:
        return None
    return (float(np.median(values)),
            float(np.percentile(values, 1)),
            float(np.percentile(values, 99)))


def alignment_window_seconds(deviations) -> float | None:
    """How far retention times actually move, as a generous multiple of the observed spread.

    Six times wider than the data needs is not conservative, it is permissive: it lets the
    aligner pair features that are not the same peak.
    """
    values = np.asarray([abs(d) for d in deviations if d == d], dtype=float)
    if values.size < 20:
        return None
    return float(4.0 * np.percentile(values, 95))


def noise_floor(intensities) -> float | None:
    """The file's own noise level: the median of the lowest decile of MS1 peak intensities."""
    values = np.asarray([i for i in intensities if i and i > 0], dtype=float)
    if values.size < 1000:
        return None
    return float(np.median(values[values <= np.percentile(values, 10)]))


def retention_tolerance(residual_sds) -> float | None:
    """The retention-model error limit, from the models' own fitted residuals.

    The 1.0 min floor it replaces is arbitrary — it is neither a property of the chromatography
    nor of the model, and on a fast gradient it is most of the run.
    """
    values = np.asarray([s for s in residual_sds if s and s == s], dtype=float)
    if values.size < 3:
        return None
    return float(3.0 * np.median(values))


# ── putting it together ────────────────────────────────────────────────────────────────────

def derive(config, evidence: dict, params: Calibration | None = None, log=None) -> Derived:
    """Measure the descriptive parameters from a first pass. Returns what changed and why.

    `evidence` carries the raw measurements — `deltas_ppm`, `fwhm_seconds`, `rt_deviations`,
    `ms1_intensities`, `model_residual_sds` — so every estimator is a pure function and can be
    tested against a distribution with a known answer.
    """
    say = log or (lambda _: None)
    params = params or Calibration()
    report = Derived()
    if not params.enabled:
        return report

    def offer(name, value, measured, current, why_skipped):
        _guard(name)
        if value is None:
            report.skipped[name] = why_skipped
            return
        report.measured[name] = measured
        report.previous[name] = current
        report.values[name] = _bounded(name, value, report)

    if params.mass_offset:
        offset = mass_offset_ppm(evidence.get("deltas_ppm", []))
        offer("mass_offset_ppm", offset, offset, getattr(config, "mass_offset_ppm", 0.0),
              "fewer than 50 matched spectra")
    if params.peak_width:
        widths = peak_width_seconds(evidence.get("fwhm_seconds", []))
        if widths is None:
            report.skipped["chrom_fwhm"] = "fewer than 50 detected features"
        else:
            centre, low, high = widths
            offer("chrom_fwhm", centre, centre, config.features.chrom_fwhm, "")
            offer("min_fwhm", low, low, config.features.min_fwhm, "")
            offer("max_fwhm", high, high, config.features.max_fwhm, "")
    if params.alignment:
        window = alignment_window_seconds(evidence.get("rt_deviations", []))
        offer("align_rt_max_difference", window, window,
              config.features.align_rt_max_difference, "fewer than 20 aligned features")
    if params.noise_floor:
        floor = noise_floor(evidence.get("ms1_intensities", []))
        offer("noise_threshold", floor * NOISE_MULTIPLIER if floor else None, floor,
              config.features.noise_threshold, "fewer than 1000 MS1 peaks sampled")
    if params.retention_model:
        limit = retention_tolerance(evidence.get("model_residual_sds", []))
        offer("retention_model_max_error", limit, limit, config.retention_model_max_error,
              "fewer than 3 fitted class models")

    say(str(report))
    return report


def apply(config, report: Derived):
    """A copy of `config` with the derived values in place, and the derivation recorded on it."""
    if not report.values:
        return config
    features = replace(
        config.features,
        **{k: v for k, v in report.values.items()
           if k in {"noise_threshold", "chrom_fwhm", "min_fwhm", "max_fwhm",
                    "align_rt_max_difference"}})
    top = {k: v for k, v in report.values.items()
           if k in {"retention_model_max_error", "mass_offset_ppm"}}
    out = replace(config, features=features, **top)
    out.calibration_report = {
        "derived": dict(report.values),
        "measured": dict(report.measured),
        "previous": dict(report.previous),
        "clamped": list(report.clamped),
        "skipped": dict(report.skipped),
        "pass": 1,
    }
    return out


def evidence_from(output, config, sample_files=3) -> dict:
    """The measurements a calibration pass needs, taken from a finished run.

    Reads a few files rather than all of them for the MS1 noise floor: it is a property of the
    instrument on the day, not of the sample, and sampling three files makes a two-pass run
    affordable enough that people will actually use it.
    """
    result = output.result
    deltas, widths, deviations, residuals = [], [], [], []
    for group in result.compound_groups:
        for compound in group.compounds:
            if compound.fwhm:
                widths.append(compound.fwhm * 60.0)
        for candidate in group.lipid_candidates:
            for lipid in candidate.identifications:
                if lipid.ppm_error == lipid.ppm_error:
                    deltas.append(lipid.ppm_error)
    for sample in output.samples:
        deviations.extend(d[0] * 60.0 for d in sample.rt_deviations)
    for model in (getattr(output, "retention_models", None) or {}).values():
        sd = getattr(model, "residual_sd", None)
        if sd and sd == sd:
            residuals.append(sd)

    intensities = []
    from pathlib import Path as _Path
    mzml_dir = _Path(config.mzml_dir)
    files = sorted(mzml_dir.glob("*.mzML"))[:sample_files] if mzml_dir.exists() else []
    if files:
        from .gapfill import _ms1
        import numpy as _np
        for path in files:
            scans = _ms1(path)
            if scans:
                intensities.extend(_np.concatenate([s[2] for s in scans[::7]]).tolist())
    return {"deltas_ppm": deltas, "fwhm_seconds": widths, "rt_deviations": deviations,
            "ms1_intensities": intensities, "model_residual_sds": residuals}


# Ratio of the configured `noise_threshold` to the measured floor that counts as sane. The
# anchored multiplier is 5.5; outside this band the configured value was fitted for a different
# detector and the run is about to be quietly wrong in one direction or the other.
PLAUSIBLE_NOISE_RATIO = (2.0, 15.0)


def check_noise_threshold(configured: float, floor: float | None) -> tuple[float, str, bool]:
    """(what to use if deriving, a message, whether the configured value is plausible here).

    The third value is what lets the caller adapt rather than warn: a configured threshold inside
    the plausible band belongs to this detector and should stand, one outside it was fitted
    somewhere else and should be replaced.

    The noise floor is the one derivable parameter that needs no second pass — it is measured from
    raw MS1 peaks, which exist before feature detection. So it is always measured, always reported,
    and only applied when asked.

    That split matters because the failure mode is **silence**. An absolute intensity fitted on one
    detector is not wrong in a way anyone notices on the next: too high and features are simply
    absent, too low and the run does not finish and gets blamed on the file. Measuring costs a few
    seconds and converts a silent failure into a loud one, without changing any result that was not
    explicitly asked to change.
    """
    if floor is None or floor <= 0:
        return configured, "noise floor could not be measured — the configured value stands", True
    derived = floor * NOISE_MULTIPLIER
    ratio = configured / floor
    plausible = PLAUSIBLE_NOISE_RATIO[0] <= ratio <= PLAUSIBLE_NOISE_RATIO[1]
    message = (f"MS1 noise floor {floor:,.0f} counts; configured noise_threshold {configured:,.0f} "
               f"is {ratio:.1f}x it (derived would be {derived:,.0f})")
    if not plausible:
        message = ("⚠ " + message + " — OUTSIDE the plausible band "
                   f"{PLAUSIBLE_NOISE_RATIO[0]:g}-{PLAUSIBLE_NOISE_RATIO[1]:g}x. This threshold is "
                   f"an ABSOLUTE INTENSITY and was fitted on a different detector. Set "
                   f"`derive_noise_threshold: true`, or set it by hand for this instrument.")
    return derived, message, plausible


# A DIA isolation window is wide and reused; a DDA one is narrow and follows the precursor.
# Both thresholds are deliberately loose — the point is to catch SWATH-style acquisition, not to
# quibble over a 4 Da DDA window on an older instrument.
DIA_MIN_WIDTH_DA = 5.0

# ⚠ Reuse, not width, is what separates SWATH from DDA — and the original 3.0 was far too low.
#
# A wide isolation window does not make an acquisition DIA. MTBLS5163 fragments through 6.9 Da
# windows and was classified as DIA at 3.0, but it records 715 distinct TRUE precursor masses and
# revisits each window only 6.4x: that is a DDA instrument being generous with its isolation, and
# calling it unsearchable would have discarded 400 correct identifications at R2 0.99.
#
# Measured across the validation set, the two populations do not overlap:
#
#     DDA          ST004797 2.6   ST004651 3.7   MSV000094718 4.1   MSV000095868 4.3
#                  ST003052 6.4   MTBLS5163 6.4
#     SWATH        ST000991 160.9        (40 fixed windows, revisited every cycle)
#
# 20 sits 3x above the highest DDA and 8x below the SWATH case. Real SWATH reuse scales with the
# run length — every window is revisited once per cycle for the whole gradient — so the gap widens
# with file size rather than narrowing.
DIA_MIN_REUSE = 20.0

# Beyond this a "window" is not isolating anything. Waters MSe reports a single 1,150 Da window
# covering the entire mass range, which is all-ion fragmentation wearing an isolation window, and
# its low reuse would otherwise let it pass as DDA.
ALL_ION_MIN_WIDTH_DA = 100.0


def looks_like_dia(mzml_files, sample_files=2, scan_limit=1500) -> str:
    """A message when the files look like DIA, empty when they look like DDA.

    This pipeline identifies by matching one spectrum to one library entry. A DIA spectrum is the
    co-fragmentation of everything in a wide window, so that premise does not hold — and the
    failure is quiet rather than loud. Fed 21 Da SWATH windows, the search loaded 353 candidate
    matches and every one failed the score filters: the right answer, reached with no indication
    that the ACQUISITION MODE was why. A user sees `exit 0` and an empty table.

    Detected from the isolation windows themselves rather than from any metadata field, because
    deposits routinely mislabel this — one states "data dependent acquisition" and ships 147 files
    with no MS2 at all.
    """
    import re
    from collections import Counter
    cv = re.compile(r"<cvParam[^>]*>")
    widths, targets, seen = Counter(), Counter(), 0
    for path in list(mzml_files)[:sample_files]:
        try:
            with Path(path).open(errors="replace") as fh:
                for line in fh:
                    for tag in cv.findall(line):
                        a = dict(re.findall(r'(\w+)="([^"]*)"', tag))
                        name = a.get("name")
                        if name == "isolation window target m/z":
                            targets[round(float(a["value"]), 1)] += 1
                            seen += 1
                        elif name == "isolation window lower offset":
                            widths[round(float(a["value"]), 2)] += 1
                    if seen > scan_limit:
                        break
        except OSError:
            # An unreadable file is skipped; a BUG is not. A bare `except Exception` here
            # swallowed a NameError from a missing import and made the whole check return
            # "looks like DDA" on textbook SWATH — silent, and indistinguishable from a
            # correct negative.
            continue
    if not targets or not widths:
        return ""
    width = max(widths) * 2
    reuse = sum(targets.values()) / len(targets)
    if width >= ALL_ION_MIN_WIDTH_DA:
        return (f"⚠ these files fragment through a {width:.0f} Da window — the whole mass range at "
                f"once. That is all-ion fragmentation, not isolation: nothing was selected, so "
                f"there is no precursor to match a library spectrum to. Feature detection and MS1 "
                f"quantitation are unaffected.")
    if width >= DIA_MIN_WIDTH_DA and reuse >= DIA_MIN_REUSE:
        return (f"⚠ these files look like DIA — {width:.0f} Da isolation windows, "
                f"{len(targets)} of them, each reused {reuse:.0f}x. This pipeline identifies by "
                f"matching ONE spectrum to ONE library entry, which a co-fragmented DIA spectrum "
                f"does not satisfy: expect few or no identifications, and do not read that as the "
                f"sample being empty. Feature detection and quantitation from MS1 are unaffected.")
    return ""


PROFILE_MIN_FRACTION = 0.5     # of sampled spectra, before the file is called profile


def load_centroided(path, log=None):
    """Load an mzML, centroiding it first if the vendor never did. Returns an MSExperiment.

    ⚠ Profile data is not a slow centroid, it is a DIFFERENT measurement, and the pipeline has no
    way to tell without asking. Each chromatographic peak arrives as a dozen raw sampling points
    rather than one fitted m/z, so what the pipeline reads as "the peak" is whichever point
    happened to be tallest in that scan — and that moves with noise.

    Measured on ST004797, across 16,000 strong peaks, tallest-raw-point against fitted centroid:

        median 1.21 ppm    90th pct 10.14 ppm    max 93.84 ppm
        over 10 ppm: 10.2%     over 20 ppm: 3.8%

    The pipeline matches at 10 ppm, so roughly one peak in ten falls outside its own tolerance, and
    the reported m/z jitters scan to scan. Mass traces fragment, identifications go missing
    silently, and mass accuracy cannot be measured at all — the error would be our conversion
    artefact rather than the instrument's.

    None of this announces itself. That deposit ran to healthy-looking feature counts of ~430 per
    file, because the peak-cap threshold hides a 20,308 peaks/scan profile file behind a plausible
    number. The median error of 1.21 ppm is why: most peaks are fine and only the tail does damage.

    Two of eleven studies in the validation set are profile (ST004797 from Sciex Analyst,
    MTBLS2016 at 52,000 peaks/scan), and neither carries a vendor original to re-convert from — so
    this cannot be pushed back onto conversion.
    """
    from pyopenms import MSExperiment, MzMLFile, PeakPickerHiRes

    experiment = MSExperiment()
    MzMLFile().load(str(path), experiment)
    if experiment.empty():
        return experiment

    # pyOpenMS reports it per spectrum: 1 = centroid, 2 = profile, 0 = unknown. Unknown is left
    # alone — picking already-centroided data damages it, so the ambiguous case does nothing.
    profile = sum(1 for s in experiment if s.getType() == 2)
    if profile < PROFILE_MIN_FRACTION * experiment.size():
        return experiment

    picked = MSExperiment()
    PeakPickerHiRes().pickExperiment(experiment, picked, True)
    if log:
        before = np.median([s.size() for s in experiment if s.getMSLevel() == 1] or [0])
        after = np.median([s.size() for s in picked if s.getMSLevel() == 1] or [0])
        log(f"⚠ {Path(path).name} is PROFILE data — the vendor conversion never centroided it. "
            f"Peak-picked on load: {before:,.0f} -> {after:,.0f} MS1 peaks per scan. Left as-is, "
            f"the m/z read from each peak is the tallest raw sampling point, which is over the "
            f"10 ppm matching tolerance for about 1 peak in 10.")
    return picked


def profile_summary(mzml_files, sample_files=4) -> str:
    """A message when the study's files are profile rather than centroided. Empty otherwise.

    `load_centroided` fixes this silently on every load, which is the right behaviour but the wrong
    record: a validation run has to be able to show from its own log what was done to the data.
    Reported once per study rather than once per injection — 224 identical lines is not a record,
    it is noise.
    """
    from pyopenms import MSExperiment, MzMLFile

    profile = checked = 0
    for path in list(mzml_files)[:sample_files]:
        experiment = MSExperiment()
        try:
            MzMLFile().load(str(path), experiment)
        except (RuntimeError, OSError):
            continue
        if experiment.empty():
            continue
        checked += 1
        if sum(1 for s in experiment if s.getType() == 2) >= PROFILE_MIN_FRACTION * experiment.size():
            profile += 1
    if not checked or not profile:
        return ""
    return (f"⚠ this study is PROFILE data ({profile} of {checked} sampled files) — the vendor "
            f"conversion never centroided it. Peak-picking on load. Left as-is the m/z read from "
            f"each peak is the tallest raw sampling point, which on measured data is past the "
            f"10 ppm matching tolerance for about 1 peak in 10.")


def acquisition(path, budget=16 << 20) -> str:
    """What ONE file acquired: "DDA", "DIA", "all-ion", "MS1 only", or "unknown".

    ⚠ Acquisition is a property of the FILE, not of the deposit. Reading one file per dataset and
    generalising produced three wrong verdicts across the validation set, including "MS1 only" for
    a study that identified 591 lipids. Real submissions mix modes:

        ST000991    153 SWATH  +  5 DDA  +  1 MS1-only
        ST003514      DDA x4   +  MS1-only x2   (its record calls this "iterative-MS/MS")
        ST004650      DDA x3   +  MS1-only x3
        ST004503      DDA x3   +  MSe x3

    ST000991 is why this matters beyond bookkeeping. Its 154 SWATH files yielded 21,345 spectral
    hits drawn from just 5 distinct precursor masses — window centres that happen to fall within a
    few mDa of a library entry, matched again in every file. All were discarded downstream, but
    they cost hours of searching. Its real identifications come from the 5 DDA files.
    """
    import re
    level = re.compile(r'name="ms level" value="(\d)"')
    target = re.compile(r'name="isolation window target m/z" value="([\d.]+)"')
    lower = re.compile(r'name="isolation window lower offset" value="([\d.]+)"')
    selected = re.compile(r"selectedIon")

    ms1 = ms2 = sel = 0
    targets: list[float] = []
    widths: list[float] = []
    read = 0
    try:
        with Path(path).open(errors="replace") as fh:
            while read < budget:
                chunk = fh.read(1 << 20)
                if not chunk:
                    break
                read += len(chunk)
                for v in level.findall(chunk):
                    if v == "1":
                        ms1 += 1
                    else:
                        ms2 += 1
                sel += len(selected.findall(chunk))
                targets.extend(float(x) for x in target.findall(chunk))
                widths.extend(float(x) for x in lower.findall(chunk))
    except OSError:
        return "unknown"

    if not ms2:
        return "MS1 only" if ms1 else "unknown"
    if not sel and not targets:
        return "all-ion"
    if widths:
        width = max(widths) * 2
        if width >= ALL_ION_MIN_WIDTH_DA:
            return "all-ion"
        unique = len({round(t, 1) for t in targets}) or 1
        if width >= DIA_MIN_WIDTH_DA and len(targets) / unique >= DIA_MIN_REUSE:
            return "DIA"
    return "DDA"


def searchable(mzml_files) -> tuple[list, dict]:
    """(files worth searching, {reason: [files skipped]}).

    Only DDA can be searched against a spectral library one spectrum at a time. Everything else is
    kept for feature detection and MS1 quantitation, which do not care how the MS2 was acquired —
    the file is excluded from the SEARCH, never from the study.

    `unknown` is searched rather than skipped. A file we cannot classify is not evidence that it is
    unsearchable, and silently dropping data on a failed read is the worse error.
    """
    # ⚠ A single path is iterable, character by character. Passed one, this returned a "kept"
    # list of 84 single letters and an empty skip dict — a confident, entirely wrong answer with
    # no error anywhere. Accept one file or many, since searching one file is a reasonable thing
    # to ask for.
    if isinstance(mzml_files, (str, Path)):
        mzml_files = [mzml_files]
    keep, skip = [], {}
    for f in mzml_files:
        verdict = acquisition(f)
        if verdict in ("DDA", "unknown"):
            keep.append(f)
        else:
            skip.setdefault(verdict, []).append(f)
    return keep, skip


def ms2_absent(mzml_files, sample_files=2) -> str:
    """A message when the files contain no MS2 at all. Empty otherwise.

    The companion to `precursor_selection_missing`, which is deliberately silent on this case: a
    run that never attempted MS2 is a different condition from one that fragmented without
    isolating, and reporting "no precursor selection" for MS1-only data would be wrong.

    Deposits mislabel this routinely. MTBLS2016 is described as data-dependent acquisition and
    ships 147 files, median 789 MB, containing 878 MS1 spectra each and zero MS2. Searched anyway,
    it spent three minutes finding nothing per file before dying on an unrelated truncated
    download — and the emptiness would have looked like a library failure.

    Identification needs MS2. Feature detection, alignment and MS1 quantitation do not, and a
    deposit like this is still perfectly good for those.
    """
    import re
    level = re.compile(r'name="ms level" value="(\d)"')
    seen_ms1 = seen_ms2 = 0
    for path in list(mzml_files)[:sample_files]:
        try:
            with Path(path).open(errors="replace") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), ""):
                    for value in level.findall(chunk):
                        if value == "1":
                            seen_ms1 += 1
                        elif value != "1":
                            seen_ms2 += 1
                    if seen_ms2:
                        return ""          # MS2 exists; nothing to say
        except OSError:
            continue
    if not seen_ms1 or seen_ms2:
        return ""
    return (f"⚠ these files contain NO MS2 — {seen_ms1:,} MS1 spectra sampled and not one at a "
            f"higher level. Nothing can be identified from spectral matching, whatever the deposit "
            f"says it acquired. Feature detection, alignment and MS1 quantitation are unaffected, "
            f"and identifications can still come from the retention model where a class has "
            f"anchors.")


def precursor_selection_missing(mzml_files, sample_files=2, scan_limit=4000) -> str:
    """A message when MS2 exists but nothing was ever isolated. Empty otherwise.

    The DIA check above asks whether the isolation windows are too WIDE. This asks the question
    one step before it: whether there are any at all.

    Found on a Leco Citius LC-HRT deposit whose 1,461 MS2 scans each carry
    `MS:1001880 In-source collision-induced dissociation` and no `selectedIon`, no isolation
    window, no precursor m/z — all-ion fragmentation, where everything eluting at that moment is
    fragmented together. There is nothing to match a library spectrum TO, because nothing was
    isolated. Every identification stage runs, finds nothing, and exits 0.

    Cheap to detect and worth detecting early: that study spent 32 minutes reaching a conclusion
    available in the first second, and the emptiness at the end was indistinguishable from a
    library gap or a threshold set too high.

    ⚠ Silent on files with no MS2 at all — that is a different condition with its own message, and
    conflating them would report "no precursor selection" for an MS1-only run.
    """
    import re
    ms2 = re.compile(r'name="ms level" value="2"')
    selected = re.compile(r"selectedIon|isolation window target")
    allion = re.compile(r"In-source collision-induced dissociation|MS:1001880")

    n_ms2 = n_sel = n_allion = 0
    for path in list(mzml_files)[:sample_files]:
        try:
            with Path(path).open(errors="replace") as fh:
                for chunk in iter(lambda: fh.read(1 << 20), ""):
                    n_ms2 += len(ms2.findall(chunk))
                    n_sel += len(selected.findall(chunk))
                    n_allion += len(allion.findall(chunk))
                    if n_ms2 > scan_limit and n_sel:
                        return ""          # precursors are present; nothing to say
        except OSError:
            continue

    if not n_ms2 or n_sel:
        return ""
    how = ("in-source CID" if n_allion else "no stated activation")
    return (f"⚠ these files carry MS2 but NO PRECURSOR SELECTION — {n_ms2:,} MS2 scans sampled, "
            f"zero with a selected ion or isolation window ({how}). This is all-ion "
            f"fragmentation, not DDA: everything eluting together was fragmented together, so "
            f"there is no precursor to match a library spectrum to and the identification stage "
            f"cannot work on this data. Feature detection and MS1 quantitation are unaffected.")


# How many MS1 peaks per scan a threshold should leave standing.
#
# ⚠ This exists because a multiple of the noise floor is NOT portable. The floor is measured
# correctly on every instrument, but `floor x 5.5` assumes an intensity distribution shaped like an
# Orbitrap's. Waters files carry 10-46x more peaks — a long tail of tiny ones the vendor's
# centroiding leaves in — so on that data the derived threshold of 88 retained 31.85% of peaks,
# 3,814 per scan, and feature detection produced 225,744 compound groups from 106 injections.
#
# Calibrated against the runs that work: Agilent at 715 retains 4.96% (233 peaks/scan), Bruker at
# 550 retains 8.89% (30/scan). A cap of 600 lands the Waters threshold at 307, retaining 5.01% —
# the same fraction as the working Agilent case, reached without assuming anything about units.
MAX_PEAKS_PER_SCAN = 600


def _not_blank(mzml_files, roles: "dict[str, str] | None"):
    """`mzml_files`, blanks excluded when `roles` says which ones they are.

    A blank is not a representative sample by definition -- it is deliberately near-empty of real
    lipid signal. Used as a proxy for "the richest injection" (`threshold_for_feature_survival`,
    by raw file size) it can rank first on disk size alone and still starve at any real threshold,
    which reads as "this batch cannot support a higher threshold" when the actual samples could.
    Caught live on MTBLS5163: `--auto-calibrate` rewrites `BCO MeOH_1` (a blank, and the single
    largest file in the batch by size) through a correction pass, its on-disk size shifts by
    ~12%, and that alone flips the trial's verdict at threshold 550 from "starves" to "survives" --
    the derived threshold then stays at 550 instead of correctly dropping to 275, and every file
    in the batch loses about half its features, calibrated or not. Falls back to every file when
    `roles` is unavailable or filtering would leave nothing, rather than raising.
    """
    files = list(mzml_files)
    if not roles:
        return files
    kept = [f for f in files if roles.get(Path(f).stem) != "blank"]
    return kept or files


def _sampled_intensities(mzml_files, sample_files=3, stride=7, roles: "dict[str, str] | None" = None):
    """(all sampled MS1 peak intensities, number of scans sampled)."""
    import numpy as _np
    from .gapfill import _ms1
    intensities, scans_seen = [], 0
    for path in _not_blank(mzml_files, roles)[:sample_files]:
        scans = _ms1(path)
        if scans:
            taken = scans[::stride]
            scans_seen += len(taken)
            intensities.append(_np.concatenate([s[2] for s in taken]))
    if not intensities:
        return None, 0
    return _np.concatenate(intensities), scans_seen


def measure_noise_floor(mzml_files, sample_files=3, stride=7,
                        roles: "dict[str, str] | None" = None) -> float | None:
    """The MS1 noise floor over a few files. Cheap enough to run on every run."""
    values, _ = _sampled_intensities(mzml_files, sample_files, stride, roles)
    if values is None:
        return None
    return noise_floor(values)


def threshold_for_peak_cap(mzml_files, cap: int = MAX_PEAKS_PER_SCAN,
                           sample_files=3, stride=7,
                           roles: "dict[str, str] | None" = None) -> float | None:
    """The lowest threshold that leaves at most `cap` MS1 peaks per scan.

    Scale-free by construction: it asks how much SURVIVES rather than what the intensity numbers
    mean, so it transfers across vendors whose units and centroiding differ. Returns None when the
    data is already sparser than the cap, in which case it has nothing to say.
    """
    import numpy as _np
    values, scans = _sampled_intensities(mzml_files, sample_files, stride, roles)
    if values is None or scans <= 0:
        return None
    values = values[values > 0]
    keep = cap * scans
    if keep >= values.size:
        return None
    return float(_np.partition(values, values.size - keep)[values.size - keep])

# A file whose richest injection yields fewer features than this is not a sparse sample — it is a
# threshold set above the data. Calibrated against the platforms that work, whose sampled files
# yield 403 (Agilent), 557 (Thermo), 611 (Bruker), 1,524 (Sciex) and 2,331 (Waters) features at
# their working thresholds. The failing case yields 0.
MIN_FEATURES_PER_FILE = 200
FEATURE_TRIAL_FILES = 2
FEATURE_TRIALS = 5


def threshold_for_feature_survival(mzml_files, threshold: float, floor: float,
                                   sample_files: int = FEATURE_TRIAL_FILES,
                                   tries: int = FEATURE_TRIALS,
                                   roles: "dict[str, str] | None" = None) -> float | None:
    """Lower `threshold` until feature detection stops starving. None when it already passes.

    A multiple of the noise floor fails in BOTH directions, and the peak cap only guards one of
    them. On a low-count TOF (Leco Citius: median MS1 peak intensity 37, floor 19) the rule derived
    103, which the monoisotopic peak of a real lipid clears but its 13C isotope — 20-60% as intense
    — does not. pyOpenMS then discards every trace as having no isotope partner, because a feature
    needs two. Measured on that data: 1,065 mass traces found, 0 features built, and 789 features
    at the raw floor. All 153 injections had a median of 2 features and the study died in
    alignment, several stages downstream of the actual cause.

    Rather than predict that from a proxy statistic, this measures the thing itself. Proxies were
    tried and do not separate the cases: retained-peak fraction puts the failing Leco run at 8.52%
    and a healthy Bruker one at 8.24%, and isotope-pair survival calls Waters starved at 5% when it
    detects 2,331 features quite happily. The feature count does separate them, because it IS the
    quantity that matters.

    The threshold is halved toward the floor and the highest surviving value is taken, so a
    detector that never had a problem pays one feature-detection run on one file and keeps its
    threshold unchanged.
    """
    from .features import FeatureParams, detect_features  # local: features imports peaks

    files = [Path(f) for f in mzml_files]
    if not files or threshold <= floor:
        return None
    # Richest injections first. Sampling a blank would read as starvation and drag the threshold
    # down for the whole study -- which is exactly why blanks are excluded by ROLE first, not
    # merely ranked last by size. Size alone is not a safe proxy: MTBLS5163's largest file by disk
    # size is `BCO MeOH_1`, a blank, and a mass-calibration rewrite of it (a pure m/z correction,
    # nothing about its real content) shifted its size by ~12% -- enough to flip this function's
    # own verdict on this one file from "starves at 550" to "survives," which then left the
    # derived threshold at 550 for the whole batch instead of correctly falling to 275, and every
    # file lost about half its features, calibrated or not.
    files = _not_blank(files, roles)
    try:
        files = sorted(files, key=lambda f: f.stat().st_size, reverse=True)[:sample_files]
    except OSError:
        files = files[:sample_files]

    def features_at(value: float) -> int:
        best = 0
        for f in files:
            try:
                found = detect_features(str(f), FeatureParams(noise_threshold=float(value)))
                fm = found[0] if isinstance(found, tuple) else found
                best = max(best, fm.size())
            except (RuntimeError, OSError, ValueError):
                continue
        return best

    if features_at(threshold) >= MIN_FEATURES_PER_FILE:
        return None

    value = threshold
    for _ in range(tries):
        value = max(value / 2.0, float(floor))
        if features_at(value) >= MIN_FEATURES_PER_FILE:
            return value
        if value <= floor:
            break
    # Nothing reached the minimum. The floor is the most permissive honest choice; a genuinely
    # sparse study will still be sparse, but it will not be sparse because of us.
    return float(floor)
