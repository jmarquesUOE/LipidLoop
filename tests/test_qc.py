"""QC arms, tested where a wrong answer would be silent.

A QC report is read by someone deciding whether to trust their data, so the failure that matters
is not a crash — it is a number that looks reasonable and is wrong.
"""
import csv
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import qc   # noqa: E402


def write_run(tmp_path, columns, areas, names=None, extra=None):
    """A minimal Final_Results.csv with the columns load_run keys off."""
    header = ["Retention Time (min)", "Quant Ion", "Polarity", "Area (max)", "Identification",
              "Lipid Class", "Features Found", "Dot Product", "Identification Source"] + columns
    path = tmp_path / "Final_Results.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for i, row in enumerate(areas):
            meta = (extra or {}).get(i, {})
            writer.writerow([meta.get("rt", 10.0 + i), meta.get("mz", 700.0 + i), "+",
                             max(row), (names or {}).get(i, ""), meta.get("class", "PC"),
                             len(columns), 900, "MS2"] + list(row))
    return path


def test_median_normalise_removes_a_loading_difference(tmp_path):
    """Two injections of identical material at different amounts must become identical."""
    areas = np.array([[100.0, 300.0], [200.0, 600.0], [50.0, 150.0]])
    out = qc.median_normalise(areas)
    assert np.allclose(out[:, 0], out[:, 1])


def test_profile_screen_catches_a_pool_that_is_faithful_in_size_but_not_in_shape(tmp_path):
    """The whole reason the screen is on profile: a full-sized pool can still be wrong."""
    rng = np.random.default_rng(0)
    base = rng.lognormal(10, 1.0, 60)
    columns = [f"Sample_{i:02d}" for i in range(10)] + ["QC_good", "QC_scrambled"]
    areas = np.column_stack([base * rng.uniform(0.9, 1.1, 60) for _ in range(10)]
                            + [base * 1.0, rng.permutation(base)])
    run = qc.load_run(write_run(tmp_path, columns, areas))
    screen = {r["injection"]: r for r in qc.profile_screen(run)}
    assert screen["QC_good"]["rho"] > 0.9
    assert screen["QC_scrambled"]["rho"] < 0.5
    assert screen["QC_scrambled"]["below_sample_band"]
    # and it is not caught by total signal, which is the point
    assert 0.8 < screen["QC_scrambled"]["total"] / screen["QC_good"]["total"] < 1.25


def test_d_ratio_falls_when_the_biology_grows(tmp_path):
    rng = np.random.default_rng(1)
    base = rng.lognormal(10, 1.0, 80)

    def build(sample_spread):
        columns = [f"Sample_{i:02d}" for i in range(12)] + [f"QC_{i:02d}" for i in range(4)]
        areas = np.column_stack(
            [base * rng.lognormal(0, sample_spread, 80) for _ in range(12)]
            + [base * rng.lognormal(0, 0.05, 80) for _ in range(4)])
        return qc.precision(qc.load_run(write_run(tmp_path, columns, areas)))

    tight = build(0.08)
    wide = build(0.60)
    assert wide["median_d_ratio"] < tight["median_d_ratio"]
    assert wide["median_d_ratio"] < 0.5


def test_run_order_finds_a_planted_drift(tmp_path):
    rng = np.random.default_rng(2)
    columns = [f"Sample_{i:02d}" for i in range(20)]
    drifting = np.array([[1000.0 * (1.15 ** j) for j in range(20)] for _ in range(10)])
    flat = rng.lognormal(8, 0.1, (30, 20))
    areas = np.vstack([drifting, flat])
    path = write_run(tmp_path, columns, areas)
    sequence = tmp_path / "seq.csv"
    sequence.write_text("Sample Type,File Name,Sample ID\n"
                        + "".join(f"Unknown,{c},P{i}\n" for i, c in enumerate(columns)))
    run = qc.load_run(path, sequence)
    result = qc.run_order(run)
    # All ten planted features must be found. A couple of the thirty flat ones may cross |rho|
    # 0.5 by chance on twenty points — that is the threshold behaving correctly, not a bug, so
    # the test bounds the false positives instead of forbidding them.
    planted = [qc.spearman(result["injection"], run.areas[i, :20]) for i in range(10)]
    assert all(abs(r) >= 0.5 for r in planted)
    assert 10 <= result["drifting"] <= 13
    assert result["loading_rho"] > 0.9


def test_confounding_separates_a_balanced_design_from_a_confounded_one(tmp_path):
    columns = [f"Sample_{i:02d}" for i in range(20)]
    areas = np.random.default_rng(3).lognormal(8, 0.2, (25, 20))
    path = write_run(tmp_path, columns, areas)
    sequence = tmp_path / "seq.csv"
    sequence.write_text("Sample Type,File Name,Sample ID\n"
                        + "".join(f"Unknown,{c},P{i}\n" for i, c in enumerate(columns)))

    def check(assign):
        meta = tmp_path / "meta.csv"
        with meta.open("w", newline="") as fh:
            writer = csv.writer(fh)
            writer.writerow(["sample", "group"])
            for i, c in enumerate(columns):
                writer.writerow([c, assign(i)])
        return {r["factor"]: r for r in qc.confounding(qc.load_run(path, sequence, meta))}

    alternating = check(lambda i: "A" if i % 2 else "B")
    assert not alternating["group"]["confounded"]

    blocked = check(lambda i: "A" if i < 10 else "B")
    assert blocked["group"]["confounded"], "one group injected first must be caught"


def test_annotation_ambiguity_separates_its_three_cases(tmp_path):
    columns = [f"Sample_{i:02d}" for i in range(8)]
    areas = np.array([[1000.0 + 10 * i for i in range(8)],      # 0 PC 34:1 @ 10.0
                      [990.0 + 10 * i for i in range(8)],       # 1 PC 16:0_18:1, same peak
                      [500.0] * 8,                              # 2 PC 34:1 again, later RT
                      # 3 and 4 share a peak AND track each other, which is what one measurement
                      # reported twice looks like. A pair that does not track is two adjacent
                      # peaks — see test_a_pair_that_does_not_track_is_not_one_measurement.
                      [400.0 + 10 * i for i in range(8)],       # 3 PE 36:2, same peak as 4
                      [800.0 + 21 * i for i in range(8)]])      # 4 SM d34:1, same peak as 3
    extra = {0: {"rt": 10.00, "mz": 760.5851}, 1: {"rt": 10.02, "mz": 760.5851},
             2: {"rt": 14.50, "mz": 760.5851}, 3: {"rt": 12.00, "mz": 744.5543},
             4: {"rt": 12.03, "mz": 744.5543}}
    names = {0: "PC 34:1", 1: "PC 16:0_18:1", 2: "PC 34:1", 3: "PE 36:2", 4: "SM d34:1"}
    run = qc.load_run(write_run(tmp_path, columns, areas, names, extra))
    # The window is one measured peak width; these fixtures sit 0.02-0.03 min apart.
    result = qc.annotation_ambiguity(run, fwhm=0.09)

    repeated = {r["name"]: r for r in result["repeated"]}
    assert repeated["PC 34:1"]["rows"] == 2
    assert repeated["PC 34:1"]["spread"] == pytest.approx(4.5, abs=0.01)

    assert any({p["a"], p["b"]} == {"PC 34:1", "PC 16:0_18:1"} for p in result["summed_pairs"]), \
        "a summed composition beside its own resolved form is not an independent second lipid"
    assert any({p["a"], p["b"]} == {"PE 36:2", "SM d34:1"} for p in result["same_peak"]), \
        "two unrelated names on one peak is the case that inflates a class count"


def test_a_pair_that_does_not_track_is_not_one_measurement(tmp_path):
    """The fix for a real false positive. A fixed 0.2 min window is more than two peak widths on
    this method, so it paired adjacent-but-distinct peaks; every pair it produced correlated
    between -0.20 and +0.53 where one measurement reported twice would correlate at about +1.
    Correlation is now a criterion rather than a printed statistic."""
    columns = [f"Sample_{i:02d}" for i in range(8)]
    areas = np.array([[400.0 + 10 * i for i in range(8)],
                      [800.0 - 21 * i for i in range(8)]])          # anti-correlated
    extra = {0: {"rt": 12.00, "mz": 744.5543}, 1: {"rt": 12.03, "mz": 744.5543}}
    names = {0: "PE 36:2", 1: "SM d34:1"}
    run = qc.load_run(write_run(tmp_path, columns, areas, names, extra))
    result = qc.annotation_ambiguity(run, fwhm=0.09)
    assert not result["same_peak"] and not result["summed_pairs"]
    assert len(result["rejected_pairs"]) == 1
    assert not any({p["a"], p["b"]} == {"PC 34:1", "PC 16:0_18:1"} for p in result["same_peak"])


def test_pool_drift_reports_size_without_correcting_anything(tmp_path):
    columns = [f"Sample_{i:02d}" for i in range(6)] + [f"QC_{i:02d}" for i in range(6)]
    rising = np.array([[1000.0] * 6 + [1000.0 * (1.2 ** j) for j in range(6)]
                       for _ in range(5)])
    flat = np.array([[1000.0] * 12 for _ in range(20)])
    areas = np.vstack([rising, flat])
    sequence = tmp_path / "seq.csv"
    sequence.write_text("Sample Type,File Name,Sample ID\n"
                        + "".join(f"Unknown,{c},P{i}\n" for i, c in enumerate(columns)))
    run = qc.load_run(write_run(tmp_path, columns, areas), sequence)
    drift = qc.pool_drift(run)
    assert drift["systematic"] >= 5
    assert drift["grid"], "thresholds must be reported so one can be chosen later"
    # nothing is altered
    assert run.areas[0, -1] == pytest.approx(1000.0 * 1.2 ** 5)


def test_outliers_catch_different_failures(tmp_path):
    """T-squared finds a sample far from the centre; DModX one the model cannot describe."""
    rng = np.random.default_rng(4)
    base = rng.lognormal(10, 0.6, 120)
    columns = [f"Sample_{i:02d}" for i in range(14)] + ["QC_1", "QC_2", "QC_3"]
    normal = [base * rng.lognormal(0, 0.15, 120) for _ in range(13)]
    extreme = base * rng.lognormal(0, 0.15, 120) * 6.0        # far from the centre
    pools = [base * rng.lognormal(0, 0.03, 120) for _ in range(3)]
    areas = np.column_stack(normal + [extreme] + pools)
    run = qc.load_run(write_run(tmp_path, columns, areas))
    space = qc.score_space(run, ("sample", "qc"))
    flags = qc.outliers(space)
    assert space["columns"][13] == "Sample_13"
    assert 13 in flags["flagged"]


def test_a_single_test_does_not_flag_an_injection(tmp_path):
    """A flag needs corroboration.

    With one test enough, a study of fourteen injections flagged the wrong ones: the two most
    extreme points on the first two components went unflagged while a middling one was flagged,
    each on one test alone. One test crossing a 95% limit is expected by chance at that size, and
    a flag the reader cannot see teaches them to distrust the panel rather than the injection.
    """
    rng = np.random.default_rng(11)
    base = rng.lognormal(10, 0.6, 150)
    columns = [f"Sample_{i:02d}" for i in range(12)] + ["QC_1", "QC_2"]
    normal = [base * rng.lognormal(0, 0.15, 150) for _ in range(11)]
    extreme = base * rng.lognormal(0, 0.15, 150) * 8.0
    pools = [base * rng.lognormal(0, 0.03, 150) for _ in range(2)]
    run = qc.load_run(write_run(tmp_path, columns, np.column_stack(normal + [extreme] + pools)))
    flags = qc.outliers(qc.score_space(run, ("sample", "qc")))

    assert set(flags["flagged"]) == set(flags["corroborated"])
    for i in flags["flagged"]:
        assert len(flags["support"][i]) >= 2, "flagged on a single test"
    for i in flags["weak"]:
        assert i not in flags["flagged"], "a single hit must not flag"


def test_a_column_that_cannot_express_a_contrast_is_not_a_study_factor():
    """Excluded on principle, not by name — the next study will name them differently.

    A real submission produced three metadata columns and the report showed all three as study
    groups: the real one, a replicate index, and `sample_type`, which read "Snap frozen duodenum
    tissue" for every sample and so could never separate anything.
    """
    metadata = {
        "S1": {"group": "A", "sample_type": "tissue", "id": "x1", "replicate": "1"},
        "S2": {"group": "A", "sample_type": "tissue", "id": "x2", "replicate": "2"},
        "S3": {"group": "B", "sample_type": "tissue", "id": "x3", "replicate": "1"},
        "S4": {"group": "B", "sample_type": "tissue", "id": "x4", "replicate": "2"},
    }
    usable = qc.usable_factors(metadata)
    assert "group" in usable
    assert "sample_type" not in usable, "one level cannot express a contrast"
    assert "id" not in usable, "one sample per level is an identifier"
    # a replicate index survives on shape alone — it is a pairing, and needs declaring not guessing
    assert "replicate" in usable


def test_an_unknown_injection_is_not_assumed_to_be_a_sample(tmp_path):
    """The roles map is built from the delivered table, which no longer carries standards columns,
    so it does not mention them at all. Defaulting an unrecognised injection to "sample" let every
    standards injection back into a panel they had just been excluded from — and the count in the
    report was the only sign."""
    import csv as _csv
    path = tmp_path / "spectra.csv"
    with path.open("w", newline="") as fh:
        w = _csv.writer(fh)
        w.writerow(["Sample", "Rank", "Delta m/z (ppm)"])
        for name in ("Sample_01", "Std_Mix_Pos_01", "Blank_Pos_01"):
            for i in range(40):
                w.writerow([name, 1, -1.0 + i * 0.01])
    # roles knows only the sample, exactly as a standards-free delivered table would
    got = qc.mass_accuracy(path, roles={"Sample_01": "sample"})
    plotted = {r["injection"] for r in got["injections"]}
    assert plotted == {"Sample_01"}, f"standards or blanks leaked in: {plotted}"


# ---------------------------------------------------------------------------
# Degenerate studies — the shape that has actually broken runs
#
# Every QC crash so far survived the whole suite: `median_pool_cv` and then
# `missing_identified`, both bare `{}` early returns indexed downstream without a
# guard, both surfacing only on a real overnight batch after the pipeline had
# already written every result file. The tests covered these functions on ordinary
# studies; nothing asked what they return when a study has no samples, no pools or
# a single injection.
#
# These are not edge cases in the deposits: ST000991_DDA is five pool injections and
# nothing else, which is exactly what broke the report.
# ---------------------------------------------------------------------------

def _degenerate(tmp_path, kind):
    """A Run with one degeneracy, built the way `load_run` would see it."""
    shapes = {
        # every injection is a pool — no samples at all (this is ST000991_DDA)
        "pools_only": ["Pool_1", "Pool_2", "Pool_3"],
        # a single injection: no variance, no pairs, nothing to correlate
        "one_injection": ["Sample_1"],
        # samples but nothing to compare them against
        "no_pools": ["Sample_1", "Sample_2", "Sample_3"],
    }
    columns = shapes[kind]
    rng = np.random.default_rng(0)
    areas = rng.lognormal(12, .4, size=(6, len(columns)))
    return qc.load_run(write_run(tmp_path, columns, areas))


@pytest.mark.parametrize("kind", ["pools_only", "one_injection", "no_pools"])
@pytest.mark.parametrize("fn", ["precision", "missingness", "coverage"])
def test_a_degenerate_study_returns_a_dict_that_can_be_rendered(tmp_path, kind, fn):
    """Never a bare {} — for the functions the report indexes WITHOUT a guard.

    That distinction is the whole contract, and getting it wrong in either direction breaks
    something:

      * `precision`, `missingness` and `coverage` are read as `d["key"]`, so `{}` is a KeyError
        and the report dies on the study whose data is most in question.
      * `run_order`, `pool_drift` and `blank_ratio` are read behind `if not d: skip` or
        `if d.get("tested")`. For those `{}` is the documented "nothing to say" signal, and
        making them return a populated dict would defeat the guard and crash the plot on a
        missing inner key instead. They are deliberately excluded here.
    """
    out = getattr(qc, fn)(_degenerate(tmp_path, kind))
    assert isinstance(out, dict)
    assert out, f"{fn} returned a bare dict on a {kind} study"


@pytest.mark.parametrize("fn", ["run_order", "pool_drift", "blank_ratio"])
def test_the_guarded_functions_stay_falsy_rather_than_half_populated(tmp_path, fn):
    """The other half of the contract, pinned so a future 'fix' cannot quietly break the guards.

    Returning a dict carrying only a reason would read as an improvement and would pass every
    consumer's `if not d` check, then fail on `d["log_ratio"]` one line later.
    """
    out = getattr(qc, fn)(_degenerate(tmp_path, "no_pools"))
    assert isinstance(out, dict)
    assert not out or "tested" in out or "log_ratio" in out or "features" in out


def test_acquisition_timestamp_is_found_however_long_the_header(tmp_path):
    """Writers place startTimeStamp at very different depths, and a fixed window misses it.

    Measured across the validation set: 4.1 kB (Thermo), 6.6 kB (Sciex), 9.5 kB and 14.5 kB
    (Waters, Agilent). The reader looked at the first 8 kB, so it found the stamp for three studies
    of five and silently returned "" for the rest — and with no injection order the drift,
    run-order and confounding sections of the report cannot run. They then say "not testable",
    which reads as a property of the data rather than a missing input.
    """
    from lipidloop.qc import _acquisition_time

    stamp = "2022-05-13T11:25:52Z"
    deep = tmp_path / "deep.mzML"
    deep.write_text("<mzML>" + ("<cvParam padding/>" * 1200)
                    + f'<run startTimeStamp="{stamp}">' + "<spectrumList count=\"1\"/></run></mzML>")
    assert len(deep.read_text()) > 20000, "the fixture must be deeper than the old 8 kB window"
    assert _acquisition_time(deep) == stamp

    # A file with no stamp must return "" without reading itself to death.
    none = tmp_path / "none.mzML"
    none.write_text("<mzML><run><spectrumList>" + ("<spectrum index=\"0\"/>" * 5000)
                    + "</spectrumList></run></mzML>")
    assert _acquisition_time(none) == ""
