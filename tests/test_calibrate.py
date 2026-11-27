"""Deriving descriptive parameters from a first pass.

The property under test is not "does it produce a number" — it is the boundary:

    Iterate on parameters that DESCRIBE the measurement.
    Never on parameters that DECIDE what counts as a hit.

A pipeline that tunes its own decision thresholds on the data produces results nobody can falsify,
so most of what follows tests that this module refuses to do that, that a pathological run cannot
silently reconfigure the pipeline, and that everything derived is written down.
"""
import sys
import math
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import calibrate   # noqa: E402
from lipidloop.calibrate import (BOUNDS, Calibration, Derived, FORBIDDEN,   # noqa: E402
                                  alignment_window_seconds, apply, derive, mass_offset_ppm,
                                  noise_floor, peak_width_seconds, retention_tolerance)
from lipidloop.pipeline import RunConfig   # noqa: E402


# ── the boundary ───────────────────────────────────────────────────────────────────────────

@pytest.mark.parametrize("name", sorted(FORBIDDEN))
def test_no_decision_rule_can_be_derived(name):
    """Every decision rule raises, by name. A convention in a docstring survives exactly as long
    as the person who wrote it; this survives a refactor."""
    with pytest.raises(ValueError, match="decision rule"):
        calibrate._guard(name)


def test_the_forbidden_list_covers_the_thresholds_that_decide_a_hit():
    """If a scoring threshold is added to the pipeline and not to FORBIDDEN, this is the test that
    should have caught it."""
    for name in ("min_dot_product", "min_reverse_dot_product", "min_purity",
                 "blank_filter.multiplier", "min_feature_count", "min_cv_improvement",
                 "min_correlation", "sn_threshold"):
        assert name in FORBIDDEN


def test_a_descriptive_parameter_is_allowed():
    for name in ("chrom_fwhm", "noise_threshold", "align_rt_max_difference", "mass_offset_ppm",
                 "min_fwhm", "max_fwhm", "retention_model_max_error"):
        calibrate._guard(name)          # must not raise


# ── the estimators, against distributions with a known answer ──────────────────────────────

def test_mass_offset_recovers_a_known_systematic_error():
    rng = np.random.default_rng(0)
    deltas = rng.normal(-1.78, 3.0, size=4000)
    assert mass_offset_ppm(deltas) == pytest.approx(-1.78, abs=0.2)


def test_mass_offset_is_robust_to_a_tail_of_bad_matches():
    """The median, not the mean: a minority of wildly wrong matches must not move the centre."""
    rng = np.random.default_rng(1)
    good = rng.normal(-1.78, 2.0, size=3000)
    junk = rng.normal(+400.0, 50.0, size=300)          # 10% nonsense
    assert mass_offset_ppm(np.concatenate([good, junk])) == pytest.approx(-1.78, abs=0.5)


def test_mass_offset_declines_on_too_little_evidence():
    assert mass_offset_ppm(np.zeros(49)) is None


def test_peak_width_uses_robust_percentiles_not_the_extremes():
    rng = np.random.default_rng(2)
    widths = list(rng.normal(5.5, 0.8, size=2000)) + [0.01, 900.0]   # two bad integrations
    centre, low, high = peak_width_seconds(widths)
    assert centre == pytest.approx(5.5, abs=0.2)
    assert 2.0 < low < 5.0 and 5.5 < high < 12.0        # the outliers set neither bound


def test_alignment_window_is_wider_than_the_observed_spread():
    rng = np.random.default_rng(3)
    deviations = rng.normal(0, 2.0, size=500)
    window = alignment_window_seconds(deviations)
    assert window > np.percentile(np.abs(deviations), 95)
    assert window < 60.0        # generous, not the 30 s it replaces being kept by accident


def test_noise_floor_tracks_the_instrument_scale():
    """The point of deriving it: the same chemistry on a detector reading ten times higher must
    give a floor ten times higher, because the configured 5000 is an absolute intensity."""
    rng = np.random.default_rng(4)
    quiet = rng.lognormal(np.log(400), 0.5, size=50_000)
    loud = quiet * 10
    assert noise_floor(loud) == pytest.approx(10 * noise_floor(quiet), rel=0.02)


def test_retention_tolerance_comes_from_the_fitted_residuals():
    assert retention_tolerance([0.20, 0.25, 0.30, 0.22]) == pytest.approx(3 * 0.235, abs=0.02)
    assert retention_tolerance([0.2, 0.3]) is None      # too few models to believe


# ── nothing derived may quietly reconfigure the run ────────────────────────────────────────

def evidence(**over):
    rng = np.random.default_rng(5)
    base = {"deltas_ppm": rng.normal(-1.8, 3.0, size=3000),
            "fwhm_seconds": rng.normal(5.5, 0.7, size=3000),
            "rt_deviations": rng.normal(0, 2.0, size=500),
            "ms1_intensities": rng.lognormal(np.log(400), 0.5, size=50_000),
            "model_residual_sds": [0.22, 0.25, 0.28, 0.24]}
    base.update(over)
    return base


def test_a_pathological_file_is_clamped_not_obeyed():
    """A noise floor a thousand times too high would delete the whole run. Bounds hold it to
    something survivable and say that they did."""
    config = RunConfig()
    huge = np.full(50_000, 5e7)
    report = derive(config, evidence(ms1_intensities=huge), Calibration(enabled=True))
    assert "noise_threshold" in report.clamped
    assert report.values["noise_threshold"] == BOUNDS["noise_threshold"][1]


def test_every_derived_value_stays_inside_its_bounds():
    config = RunConfig()
    report = derive(config, evidence(), Calibration(enabled=True))
    for name, value in report.values.items():
        low, high = BOUNDS[name]
        assert low <= value <= high, name


def test_disabled_derives_nothing_and_leaves_the_config_identical():
    config = RunConfig()
    report = derive(config, evidence(), Calibration(enabled=False))
    assert not report.values
    assert apply(config, report) is config


def test_missing_evidence_is_skipped_with_a_reason_not_guessed():
    config = RunConfig()
    report = derive(config, {"deltas_ppm": [], "fwhm_seconds": [], "rt_deviations": [],
                             "ms1_intensities": [], "model_residual_sds": []},
                    Calibration(enabled=True))
    assert not report.values
    assert set(report.skipped) >= {"mass_offset_ppm", "chrom_fwhm", "align_rt_max_difference",
                                   "noise_threshold", "retention_model_max_error"}
    assert all(reason for reason in report.skipped.values())


# ── what it does to the configuration ──────────────────────────────────────────────────────

def test_the_derivation_is_recorded_on_the_config():
    """The condition the whole idea rests on. Without this, "run twice" becomes "run until it
    looks good" and the result stops being reproducible from its configuration."""
    config = RunConfig()
    report = derive(config, evidence(), Calibration(enabled=True))
    updated = apply(config, report)
    recorded = updated.calibration_report
    assert recorded["pass"] == 1
    assert set(recorded["derived"]) == set(report.values)
    for name, value in recorded["derived"].items():
        assert name in recorded["previous"], f"{name} recorded without what it replaced"
        assert recorded["previous"][name] != value or True     # equality is fine, absence is not


def test_applying_changes_only_the_parameters_that_were_derived():
    config = RunConfig()
    report = derive(config, evidence(), Calibration(enabled=True))
    updated = apply(config, report)
    for name in ("min_feature_count", "rt_filter_multiplier", "excluded_classes",
                 "write_scores", "keep_filtered", "annotate_fatty_acids"):
        assert getattr(updated, name) == getattr(config, name)
    assert updated.blank_filter == config.blank_filter
    assert updated.split_peaks == config.split_peaks
    assert updated.presence_filter == config.presence_filter


def test_the_original_config_is_not_mutated():
    config = RunConfig()
    before = config.features.noise_threshold
    report = derive(config, evidence(), Calibration(enabled=True))
    apply(config, report)
    assert config.features.noise_threshold == before


def test_a_second_pass_over_the_same_evidence_does_not_drift():
    """Convergence. If deriving from calibrated data moved the parameters again, the procedure
    would oscillate and 'two passes' would be an arbitrary stopping point."""
    config = RunConfig()
    first = derive(config, evidence(), Calibration(enabled=True))
    once = apply(config, first)
    second = derive(once, evidence(), Calibration(enabled=True))
    for name, value in second.values.items():
        assert value == pytest.approx(first.values[name], rel=1e-9), name


def test_the_noise_multiplier_reproduces_the_validated_setting_on_the_validated_instrument():
    """A derivation that changes the answer on the data it was tuned against is not a calibration,
    it is a different pipeline.

    The floor measured on the Lumos brain files is 920 counts (positive) and 866 (negative),
    against a configured `noise_threshold` of 5,000. Reproduce that floor exactly and the derived
    threshold must land back on the validated value.
    """
    # A distribution whose lowest-decile median is exactly 920, as measured on the real file.
    low = np.full(10_000, 920.0)
    rest = np.linspace(921.0, 5e6, 90_000)
    derived = noise_floor(np.concatenate([low, rest])) * calibrate.NOISE_MULTIPLIER
    assert derived == pytest.approx(5000, rel=0.10), derived


def test_the_report_says_what_changed():
    config = RunConfig()
    report = derive(config, evidence(), Calibration(enabled=True))
    text = str(report)
    assert "noise_threshold" in text and "->" in text
    assert str(Derived()) == "calibration: nothing derived"


# ── the evidence actually reaches the estimators ───────────────────────────────────────────

class _Model:
    def __init__(self, sd): self.residual_sd = sd


class _Lipid:
    def __init__(self, ppm): self.ppm_error = ppm


class _Candidate:
    def __init__(self, lipids): self.identifications = lipids


class _Compound:
    def __init__(self, fwhm): self.fwhm = fwhm


class _Group:
    def __init__(self, compounds, candidates):
        self.compounds, self.lipid_candidates = compounds, candidates


class _Result:
    def __init__(self, groups): self.compound_groups = groups


class _Sample:
    def __init__(self, deviations): self.rt_deviations = deviations


class _Output:
    def __init__(self, groups, samples, models):
        self.result, self.samples, self.retention_models = _Result(groups), samples, models


def test_evidence_is_collected_from_the_shape_a_real_run_returns(tmp_path):
    """Found by running it: `retention_model_max_error` was skipped on every real run because
    PipelineOutput did not expose the fitted models, so the estimator had nothing to read. A
    tested estimator that never receives evidence is dead code with a passing test."""
    from lipidloop.calibrate import evidence_from
    groups = [_Group([_Compound(0.09)] * 30,
                     [_Candidate([_Lipid(-1.8)] * 30)]) for _ in range(30)]
    samples = [_Sample([(0.03, 5.0, 1.0)] * 40) for _ in range(3)]
    models = {"PC": _Model(0.22), "PE": _Model(0.25), "TG": _Model(0.31)}
    config = RunConfig(mzml_dir=str(tmp_path))          # no mzML: the noise floor should skip
    ev = evidence_from(_Output(groups, samples, models), config)

    assert len(ev["deltas_ppm"]) == 900 and len(ev["fwhm_seconds"]) == 900
    assert len(ev["rt_deviations"]) == 120
    assert ev["model_residual_sds"] == [0.22, 0.25, 0.31]
    assert ev["fwhm_seconds"][0] == pytest.approx(5.4)      # minutes converted to seconds

    report = derive(config, ev, Calibration(enabled=True))
    assert "retention_model_max_error" in report.values     # the estimator now fires
    assert "noise_threshold" in report.skipped              # and the absent one says why


# ── the noise floor is measured on every run, and applied only when asked ───────────────────

def test_a_matching_threshold_reports_quietly():
    from lipidloop.calibrate import check_noise_threshold
    derived, message, plausible = check_noise_threshold(5000, 1009)      # the validated Lumos case
    assert derived == pytest.approx(5550, rel=0.01)
    assert not message.startswith("⚠") and "5.0x" in message


def test_a_threshold_fitted_on_another_detector_is_shouted_about():
    """The whole point. An absolute intensity on a new instrument fails silently — too high and
    features are simply absent, too low and the run does not finish and the file gets blamed."""
    from lipidloop.calibrate import check_noise_threshold
    for floor in (20_000, 60):                  # far too quiet a threshold, and far too loud
        _, message, _ = check_noise_threshold(5000, floor)
        assert message.startswith("⚠") and "ABSOLUTE INTENSITY" in message


def test_an_unmeasurable_floor_leaves_the_configured_value_alone():
    from lipidloop.calibrate import check_noise_threshold
    derived, message, plausible = check_noise_threshold(5000, None)
    assert derived == 5000 and "configured value stands" in message


def test_measuring_the_floor_needs_no_completed_pass(tmp_path, monkeypatch):
    """It reads raw MS1, which exists before feature detection — so unlike the other derivable
    parameters this one costs no second pass, which is why it can be on by default one day."""
    import lipidloop.calibrate as cal
    rng = np.random.default_rng(9)
    fake = [(0.0, np.array([1.0]), rng.lognormal(np.log(400), 0.5, size=20_000))]
    monkeypatch.setattr(cal, "_ms1", lambda path: fake, raising=False)
    monkeypatch.setattr("lipidloop.gapfill._ms1", lambda path: fake)
    floor = cal.measure_noise_floor([tmp_path / "a.mzML"], stride=1)
    assert floor == pytest.approx(cal.noise_floor(fake[0][2]), rel=1e-9)


def test_auto_keeps_a_plausible_threshold_and_replaces_an_implausible_one():
    """"auto" must not move the facility's own runs, and must not leave a foreign one broken.

    `noise_threshold` is an absolute intensity, so it cannot be portable. 5,000 counts is right on
    the facility Orbitrap (floor ~1,009, ratio 5.0) and badly wrong on an Agilent 6545 QTOF
    (floor 214, ratio 23.3) — where it cut negative-mode features from ~590 to ~170 per file and
    took recall against a published lipid list from 67.8% down to 50.9%, silently, because a
    too-high threshold produces absence rather than error.
    """
    from lipidloop.calibrate import check_noise_threshold

    _, _, plausible = check_noise_threshold(5000, 1009)
    assert plausible, "the facility's own fitted value belongs to its own detector"

    derived, _, plausible = check_noise_threshold(5000, 214)
    assert not plausible, "23x the floor cannot have been fitted for this detector"
    assert 1000 < derived < 1400, f"derived from the measured floor, got {derived}"

    # Unmeasurable floor is not a licence to change anything.
    _, _, plausible = check_noise_threshold(5000, None)
    assert plausible, "with no measurement, the configured value stands"


def test_dia_acquisition_is_detected_and_dda_is_left_alone(tmp_path):
    """Fed DIA, this pipeline finds nothing — and must say why rather than look empty.

    It identifies by matching ONE spectrum to ONE library entry, which a co-fragmented DIA
    spectrum does not satisfy. On real SWATH data the search loaded 353 candidate matches and
    every one failed the score filters: the right answer, reached with no indication that the
    ACQUISITION MODE was the reason. The user sees exit 0 and an empty table.

    Detected from the isolation windows themselves, not from metadata — deposits mislabel this
    routinely, one claiming "data dependent acquisition" while shipping 147 files with no MS2.
    """
    from lipidloop.calibrate import looks_like_dia

    def write(path, windows, repeats):
        rows = []
        for _ in range(repeats):
            for centre, offset in windows:
                rows.append(
                    f'<cvParam name="isolation window target m/z" value="{centre}"/>'
                    f'<cvParam name="isolation window lower offset" value="{offset}"/>')
        path.write_text("<mzML>" + "".join(rows) + "</mzML>")

    dia = tmp_path / "swath.mzML"
    write(dia, [(300 + 20 * i, 10.5) for i in range(40)], 30)
    assert "look like DIA" in looks_like_dia([dia])

    # DDA: many distinct narrow windows, each used once or twice
    dda = tmp_path / "dda.mzML"
    write(dda, [(300 + 0.7 * i, 0.65) for i in range(400)], 1)
    assert looks_like_dia([dda]) == "", "narrow precursor-driven windows are DDA"


def test_sn_position_evidence_is_reported_never_asserted():
    """sn-2 evidence is a ratio, not a proof — the name must keep `_`.

    In negative-mode dissociation of a diacyl glycerophospholipid the sn-2 carboxylate is
    typically the more intense, because that ester is more labile. It is a tendency: the ratio
    moves with chain length, unsaturation, collision energy and instrument, and inverts for some
    species. So this module reports; it must never produce a `/` name, which under LSI asserts
    proven sn-position.
    """
    from lipidloop.sn_position import evidence, MIN_RATIO

    lib_mz = [255.2330, 281.2486, 184.0733]
    lib_fa = ["16:0", "18:1", ""]
    smp_mz = [255.2331, 281.2487, 184.0734]

    strong = evidence("PC", "PC 16:0_18:1", lib_mz, lib_fa, smp_mz, [1000.0, 3200.0, 9999.0])
    assert strong.favoured_sn2 == "18:1"
    assert "/" not in str(strong), "never assert a proven sn-position"

    weak = evidence("PC", "PC 16:0_18:1", lib_mz, lib_fa, smp_mz, [3000.0, 3200.0, 9999.0])
    assert not weak.favoured_sn2, f"{weak.ratio:.2f}x is under {MIN_RATIO}x and must not call"

    # TG uses neutral losses and the OPPOSITE convention — sn-2 is the least favoured loss.
    assert not evidence("TG", "TG 16:0_18:1_18:2", lib_mz, lib_fa,
                        smp_mz, [1000.0, 3200.0, 9999.0]).favoured_sn2

    # identical chains carry no arrangement information
    assert not evidence("PC", "PC 18:1_18:1", lib_mz, lib_fa,
                        smp_mz, [1000.0, 3200.0, 9999.0]).favoured_sn2


def test_internal_standard_library_is_searched_only_when_standards_are_declared():
    """Searching for spiked standards in an unspiked study names endogenous peaks as labelled ones.

    Defaulting this on produced `PC d7-18:1_15:0` — a deuterated SPLASH standard — on a study that
    contained no SPLASH. A labelled standard sits a few Da from its endogenous analogue, so an
    ISTD library matches something on almost any run, and the result is a spiked compound reported
    as biology that would then be quantified and tested like biology.
    """
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import importlib
    rs = importlib.import_module("run_study")

    assert rs.DEFAULTS["standards"] == "", \
        "a standards file must be declared per study, never inherited from a default"
    source = (Path(__file__).resolve().parents[1] / "scripts" / "run_study.py").read_text()
    assert 'extra.pop("search_istd", None) or config.get("standards")' in source, \
        "the ISTD library must be tied to a declared standards file"


def test_an_unalignable_injection_does_not_discard_the_study():
    """One sparse injection must not kill a 224-file run.

    Pose clustering fits its transformation from pairs matched against the reference — the largest
    map. A much sparser injection may yield too few pairs within tolerance, and pyOpenMS raises
    `no data points for 'linear' model`. On a real study a 945-feature blank killed a 224-file run,
    and the same error took out THREE datasets in a single overnight batch.

    A map that cannot be aligned keeps its own retention times — what an identity transformation
    would give it. A small local inaccuracy beats losing the study, and the affected injections are
    named so the effect is visible rather than surfacing later as an unexplained retention shift.
    """
    from pathlib import Path
    source = (Path(__file__).resolve().parents[1] / "src" / "lipidloop" / "features.py").read_text()
    block = source[source.index("aligner.align(feature_map, transformation)"):]
    guarded = block[:block.index("transformer.transformRetentionTimes")]
    assert "except RuntimeError" in guarded, "alignment must be guarded per map"
    assert "could not be aligned" in source, "and the affected injections must be named"


def test_the_peak_cap_floors_a_too_low_derived_threshold(tmp_path):
    """A multiple of the noise floor is not portable; how much SURVIVES is.

    The floor is measured correctly on every instrument, but `floor x 5.5` assumes an
    Orbitrap-shaped intensity distribution. Waters files carry 10-46x more peaks — a long tail the
    vendor's centroiding leaves in — so the rule derived 88, retained 31.85% of peaks (3,814 per
    scan), and feature detection produced 225,744 compound groups from 106 injections. The
    identifications were made and could not find their peak in the pile.

    Calibrated on the runs that work: Agilent at 715 retains 4.96%, Bruker at 550 retains 8.89%. A
    600-peak cap puts the Waters threshold at 307 — 5.01% retained, the same fraction as Agilent,
    reached without assuming anything about intensity units.
    """
    import numpy as np
    from lipidloop.calibrate import threshold_for_peak_cap, MAX_PEAKS_PER_SCAN

    rng = np.random.default_rng(0)

    class FakeScans:
        """Stands in for `_ms1`: 40 scans, each mostly noise with a few real peaks."""

        def __init__(self, noise, signal):
            self.scans = [(0.0, np.zeros(noise + signal),
                           np.concatenate([rng.uniform(1, 20, noise),
                                           rng.uniform(500, 5000, signal)]))
                          for _ in range(40)]

    import lipidloop.calibrate as cal
    original = cal._ms1 if hasattr(cal, "_ms1") else None

    # noise-rich, like Waters: 5,000 peaks per scan, only 100 of them real
    fake = FakeScans(5000, 100)
    import lipidloop.gapfill as gf
    saved = gf._ms1
    gf._ms1 = lambda path: fake.scans
    try:
        thr = threshold_for_peak_cap(["x.mzML"], cap=MAX_PEAKS_PER_SCAN, stride=1)
        assert thr is not None, "noise-rich data must yield a cap"
        kept = sum(int((s[2] >= thr).sum()) for s in fake.scans) / len(fake.scans)
        assert kept <= MAX_PEAKS_PER_SCAN + 1, f"cap not respected: {kept:.0f} per scan"

        # already-sparse data: the cap has nothing to say and must not invent a threshold
        sparse = FakeScans(50, 20)
        gf._ms1 = lambda path: sparse.scans
        assert threshold_for_peak_cap(["x.mzML"], cap=MAX_PEAKS_PER_SCAN, stride=1) is None
    finally:
        gf._ms1 = saved


def test_the_ceiling_lowers_a_threshold_that_starves_feature_detection(monkeypatch):
    """`floor x 5.5` fails in BOTH directions, and the peak cap only guards one of them.

    On a low-count TOF (Leco Citius LC-HRT: median MS1 peak intensity 37, floor 19) the rule
    derived 103. A real lipid's monoisotopic peak clears that; its 13C isotope, 20-60% as intense,
    does not. A feature needs two traces, so pyOpenMS discarded every one for having no isotope
    partner — 1,065 mass traces found, 0 features built, 789 features at the raw floor. All 153
    injections came out at a median of 2 features and the study died in map alignment, several
    stages downstream of the cause.

    The ceiling measures features directly rather than predicting the collapse from a statistic,
    because the proxies do not separate the cases: retained-peak fraction rates this failing run
    (8.52%) ABOVE a healthy Bruker one (8.24%).
    """
    import lipidloop.calibrate as cal

    calls = []

    class FakeMap:
        def __init__(self, n): self._n = n
        def size(self): return self._n

    def fake_detect(path, params):
        # Reproduces the Leco cliff: nothing survives at 103, plenty at 51 and below.
        calls.append(params.noise_threshold)
        return FakeMap(0 if params.noise_threshold > 60 else 800)

    import lipidloop.features as feat
    monkeypatch.setattr(feat, "detect_features", fake_detect)

    got = cal.threshold_for_feature_survival(["a.mzML"], threshold=103.0, floor=19.0)
    assert got is not None, "a starved threshold must be lowered"
    assert got == 51.5, f"expected one halving to 51.5, got {got}"
    assert calls[0] == 103.0, "the derived threshold must be tried first"

    # A detector that never had the problem keeps its threshold and pays one detection run.
    calls.clear()
    monkeypatch.setattr(feat, "detect_features", lambda p, params: FakeMap(611))
    assert cal.threshold_for_feature_survival(["a.mzML"], threshold=550.0, floor=100.0) is None

    # Nothing reaches the minimum anywhere: fall back to the floor rather than to the starved value.
    monkeypatch.setattr(feat, "detect_features", lambda p, params: FakeMap(3))
    assert cal.threshold_for_feature_survival(["a.mzML"], threshold=103.0, floor=19.0) == 19.0


def test_a_large_blank_must_not_stand_in_for_the_richest_injection(tmp_path, monkeypatch):
    """Caught live on MTBLS5163. `BCO MeOH_1`, a blank, is the single largest file in the batch by
    raw disk size -- large in the way a dirty background is large, not in the way a real sample
    is. Sorting "richest injection" by size alone picks it, and `--auto-calibrate`'s correction
    pass (a pure m/z rescale, nothing about real content) shifted its size by ~12% -- enough on its
    own to flip this function's verdict at threshold 550 from "starves" to "survives." The derived
    threshold then stayed at 550 for the whole batch instead of correctly falling to 275, and every
    file lost about half its features, calibrated or not.

    Without `roles`, the blank (larger) is sampled and the starved threshold is wrongly kept. With
    `roles` marking it a blank, the real sample is sampled instead and the threshold is correctly
    lowered -- same files, same sizes, the only difference is which one gets excluded from the
    "richest" pool.
    """
    import lipidloop.calibrate as cal
    import lipidloop.features as feat

    blank = tmp_path / "BCO_MeOH_1.mzML"
    sample = tmp_path / "10_43_01.mzML"
    blank.write_bytes(b"0" * 2_000_000)     # the larger file on disk
    sample.write_bytes(b"0" * 1_000_000)    # smaller, but the one with real content

    class FakeMap:
        def __init__(self, n): self._n = n
        def size(self): return self._n

    def fake_detect(path, params):
        # Mirrors the real bug: the blank looks fine at 550 regardless (it was never the thing
        # that should decide this); the real sample starves at 550 but is fine once halved to 275,
        # same as the real MTBLS5163 case.
        if str(path) == str(blank):
            return FakeMap(800)
        return FakeMap(800 if params.noise_threshold <= 275.0 else 3)

    monkeypatch.setattr(feat, "detect_features", fake_detect)

    # Old behaviour (no roles): the blank, being larger, is sampled and reports no starvation.
    assert cal.threshold_for_feature_survival([sample, blank], threshold=550.0, floor=100.0,
                                              roles=None) is None

    # Fixed behaviour: told which file is the blank, it samples the real one and finds the
    # starvation the blank was masking.
    roles = {"BCO_MeOH_1": "blank", "10_43_01": "sample"}
    got = cal.threshold_for_feature_survival([sample, blank], threshold=550.0, floor=100.0,
                                             roles=roles)
    assert got is not None, "excluding the blank must expose the real starvation"
    assert got == 275.0, f"expected one halving of 550 to 275, got {got}"


def test_all_ion_fragmentation_is_reported_before_the_search(tmp_path):
    """MS2 with no precursor selection cannot be searched, and says so in the first second.

    The DIA check asks whether isolation windows are too WIDE. This asks the question one step
    before: whether there are any. A Leco Citius deposit ships 1,461 MS2 scans per file, each with
    `MS:1001880 In-source collision-induced dissociation` and no selectedIon, no isolation window,
    no precursor m/z. Everything eluting was fragmented together, so there is nothing to match a
    library spectrum to — but the search runs to completion, finds nothing, and exits 0. That study
    spent 32 minutes reaching a conclusion the first file answers immediately, and the emptiness
    was indistinguishable from a library gap or a threshold set too high.
    """
    from lipidloop.calibrate import precursor_selection_missing

    allion = tmp_path / "allion.mzML"
    allion.write_text(
        '<spectrum><cvParam name="ms level" value="2"/>'
        '<precursorList count="1"><precursor spectrumRef="scan=0"><activation>'
        '<cvParam accession="MS:1001880" name="In-source collision-induced dissociation"/>'
        '</activation></precursor></precursorList></spectrum>' * 20)
    msg = precursor_selection_missing([allion])
    assert "NO PRECURSOR SELECTION" in msg
    assert "in-source CID" in msg
    assert "quantitation are unaffected" in msg, "MS1 work is still valid and must be said so"

    dda = tmp_path / "dda.mzML"
    dda.write_text(
        '<spectrum><cvParam name="ms level" value="2"/><precursorList><precursor>'
        '<selectedIonList><selectedIon><cvParam name="selected ion m/z" value="760.5"/>'
        '</selectedIon></selectedIonList></precursor></precursorList></spectrum>' * 20)
    assert precursor_selection_missing([dda]) == ""

    # MS1-only is a DIFFERENT condition with its own message; claiming "no precursor selection"
    # for a run that never attempted MS2 would be wrong.
    ms1 = tmp_path / "ms1.mzML"
    ms1.write_text('<spectrum><cvParam name="ms level" value="1"/></spectrum>' * 20)
    assert precursor_selection_missing([ms1]) == ""


def test_files_with_no_ms2_say_so_rather_than_searching_them(tmp_path):
    """MS1-only is its own condition, distinct from fragmenting without isolating.

    MTBLS2016 is described by its deposit as data-dependent acquisition and ships 147 files,
    median 789 MB, each holding 878 MS1 spectra and zero MS2. Searched anyway it found nothing
    per file for three minutes, then died on an unrelated truncated download — and the emptiness
    would have read as a library failure rather than as an acquisition that cannot be searched.
    """
    from lipidloop.calibrate import ms2_absent, precursor_selection_missing

    ms1 = tmp_path / "ms1.mzML"
    ms1.write_text('<spectrum><cvParam name="ms level" value="1"/></spectrum>' * 30)
    msg = ms2_absent([ms1])
    assert "NO MS2" in msg
    assert "quantitation are unaffected" in msg, "MS1 work is still valid and must be said so"
    # The two checks must not both fire on one file — they describe different failures.
    assert precursor_selection_missing([ms1]) == ""

    dda = tmp_path / "dda.mzML"
    dda.write_text('<spectrum><cvParam name="ms level" value="1"/></spectrum>'
                   '<spectrum><cvParam name="ms level" value="2"/><selectedIon/></spectrum>' * 30)
    assert ms2_absent([dda]) == ""


def _windows(path, windows, repeats, level="2"):
    rows = []
    for _ in range(repeats):
        for centre, offset in windows:
            rows.append(
                f'<spectrum><cvParam name="ms level" value="{level}"/><selectedIon/>'
                f'<cvParam name="isolation window target m/z" value="{centre}"/>'
                f'<cvParam name="isolation window lower offset" value="{offset}"/></spectrum>')
    path.write_text("<mzML>" + "".join(rows) + "</mzML>")


def test_a_wide_isolation_window_alone_does_not_make_an_acquisition_dia(tmp_path):
    """Reuse separates SWATH from DDA; width does not, and the old 3.0 threshold discarded good data.

    MTBLS5163 fragments through 6.9 Da windows and was called DIA at the old threshold. It records
    715 distinct TRUE precursor masses and revisits each window 6.4 times — a DDA instrument being
    generous with isolation. Calling it unsearchable would have thrown away 400 correct
    identifications fitting their retention surfaces at R2 0.99.

    Measured across the validation set the populations do not overlap: DDA runs 2.6 to 6.4, the one
    confirmed SWATH study 160.9. Real SWATH reuse grows with run length — every window is revisited
    once per cycle for the whole gradient — so the gap widens with file size.
    """
    from lipidloop.calibrate import looks_like_dia, acquisition

    # MTBLS5163's shape: many windows, each revisited a handful of times.
    wide_dda = tmp_path / "wide_dda.mzML"
    _windows(wide_dda, [(400 + i, 3.45) for i in range(60)], repeats=6)
    assert looks_like_dia([wide_dda]) == "", "6x reuse is DDA, however wide the window"
    assert acquisition(wide_dda) == "DDA"

    # ST000991's shape: few fixed windows, revisited every cycle.
    swath = tmp_path / "swath.mzML"
    _windows(swath, [(400 + 21 * i, 10.5) for i in range(8)], repeats=60)
    assert "look like DIA" in looks_like_dia([swath])
    assert acquisition(swath) == "DIA"

    # Waters MSe: one window covering everything. Low reuse would otherwise pass it as DDA.
    mse = tmp_path / "mse.mzML"
    _windows(mse, [(600.0, 575.0)], repeats=4)
    assert "all-ion" in looks_like_dia([mse])
    assert acquisition(mse) == "all-ion"


def test_the_search_runs_per_file_so_a_mixed_deposit_keeps_its_dda_files(tmp_path):
    """A deposit's acquisition is a distribution over its files, not one property.

    ST000991 is 153 SWATH + 5 DDA + 1 MS1-only in one submission. Searching the SWATH files
    produced 21,345 hits from 5 distinct precursor masses — window centres landing within a few mDa
    of a library entry, rediscovered in every file — all discarded downstream at a cost of hours.
    The study's real identifications came from the 5 DDA files, which a whole-dataset verdict would
    have thrown away with them.
    """
    from lipidloop.calibrate import searchable

    dda = tmp_path / "a_dda.mzML"
    _windows(dda, [(400 + i, 0.5) for i in range(40)], repeats=4)
    swath = tmp_path / "b_swath.mzML"
    _windows(swath, [(400 + 21 * i, 10.5) for i in range(8)], repeats=60)
    ms1 = tmp_path / "c_ms1.mzML"
    ms1.write_text('<spectrum><cvParam name="ms level" value="1"/></spectrum>' * 20)

    keep, skip = searchable([dda, swath, ms1])
    assert keep == [dda], "only the DDA file is searchable"
    assert set(skip) == {"DIA", "MS1 only"}
    # Skipped files are excluded from the SEARCH, never from the study — quantitation needs them.
    assert skip["DIA"] == [swath] and skip["MS1 only"] == [ms1]

    # A file we cannot classify is searched, not silently dropped.
    broken = tmp_path / "d_broken.mzML"
    broken.write_text("<mzML></mzML>")
    keep, _ = searchable([broken])
    assert keep == [broken], "unknown must not be treated as unsearchable"


def test_profile_data_is_centroided_on_load_and_centroid_data_is_left_alone(tmp_path, monkeypatch):
    """Profile data is not a slow centroid, it is a different measurement.

    Each chromatographic peak arrives as a dozen raw sampling points rather than one fitted m/z, so
    what the pipeline reads as "the peak" is whichever point happened to be tallest that scan — and
    that moves with noise. Measured on ST004797 over 16,000 strong peaks, tallest-raw-point against
    fitted centroid: median 1.21 ppm, 90th pct 10.14 ppm, max 93.84 ppm, with 10.2% beyond the
    10 ppm the pipeline matches at.

    None of it announces itself. That deposit ran to healthy-looking ~430 features per file because
    the peak-cap threshold hides a 20,308 peaks/scan profile file behind a plausible number, and the
    MEDIAN error is only 1.21 ppm — most peaks are fine, and only the tail does the damage.

    Picking already-centroided data damages it, so the unknown case must do nothing rather than
    guess.
    """
    import lipidloop.calibrate as cal

    class Spec:
        def __init__(self, kind, n): self._k, self._n = kind, n
        def getType(self): return self._k
        def getMSLevel(self): return 1
        def size(self): return self._n

    class Exp:
        def __init__(self, specs): self._s = specs
        def __iter__(self): return iter(self._s)
        def size(self): return len(self._s)
        def empty(self): return not self._s

    picked = {"called": False}

    class FakePicker:
        def pickExperiment(self, src, dest, *a):
            picked["called"] = True
            dest._s = [Spec(1, 40) for _ in src._s]

    def fake_import(profile_kind):
        import types
        mod = types.SimpleNamespace(
            MSExperiment=lambda: Exp([]),
            MzMLFile=lambda: types.SimpleNamespace(
                load=lambda p, e: setattr(e, "_s", [Spec(profile_kind, 600)] * 10)),
            PeakPickerHiRes=FakePicker)
        return mod

    monkeypatch.setitem(sys.modules, "pyopenms", fake_import(2))   # 2 = profile
    logged = []
    cal.load_centroided("x.mzML", log=logged.append)
    assert picked["called"], "profile data must be peak-picked"
    assert "PROFILE" in logged[0] and "10 ppm" in logged[0], "the log must say why it matters"

    picked["called"] = False
    monkeypatch.setitem(sys.modules, "pyopenms", fake_import(1))   # 1 = centroid
    logged.clear()
    cal.load_centroided("x.mzML", log=logged.append)
    assert not picked["called"], "already-centroided data must be left untouched"
    assert not logged

    picked["called"] = False
    monkeypatch.setitem(sys.modules, "pyopenms", fake_import(0))   # 0 = unknown
    cal.load_centroided("x.mzML", log=logged.append)
    assert not picked["called"], "unknown must not be picked — picking centroids damages them"


def test_searchable_accepts_one_file_without_iterating_its_characters(tmp_path):
    """A single path is iterable. Passed one, this returned 84 single letters as 'kept' files.

    No exception, no warning — a confident answer that was entirely wrong, of exactly the kind
    every other guard in this module exists to prevent. The plural parameter name is not enough
    protection when the wrong type is silently valid.
    """
    from lipidloop.calibrate import searchable

    empty = tmp_path / "nothing.mzML"
    empty.write_text("<mzML></mzML>")

    keep, skip = searchable(empty)
    assert len(keep) + sum(len(v) for v in skip.values()) == 1, \
        "one path in must mean one file considered, not one per character"
    assert all(str(f).endswith(".mzML") for f in keep), "kept entries must be paths"

    keep_many, _ = searchable([empty, empty])
    assert len(keep_many) == 2
