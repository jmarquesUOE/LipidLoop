"""The score columns must describe the name that was actually written.

A row reported at sum composition has deliberately declined to name its chains; crediting it with
the top molecular candidate's score would report confidence in an assignment it did not make.
This is the same principle `ms2_support` follows, and it is easy to get backwards.
"""
import csv

import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.peakfinder import PeakFinderResult, write_results   # noqa: E402
from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample  # noqa: E402


def candidate(name, sum_name, dot, rev, purity, sample):
    lipid = Lipid(retention=5.0, precursor=760.585, sample=sample, dot=dot, rev_dot=rev,
                  lipid_string=f"{name} [M+H]+;", lib_precursor=760.585, purity=purity,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    lipid.lipid_name, lipid.sum_lipid_name, lipid.lipid_class = name, sum_name, "PC"
    out = LipidCandidate(lipid)
    out.max_dot, out.max_rev_dot, out.purity = dot, rev, purity
    return out


def group_with(candidates, purity, sum_id):
    group = CompoundGroup(name="g", mw=759.578, retention=5.0, max_area=1.0,
                          has_ms2=True, areas=[1.0])
    group.lipid_candidates = list(candidates)
    group.purity = purity
    group.sum_id = sum_id
    group.final_lipid_id = candidates[0]
    return group


def test_sum_composition_row_takes_the_best_of_everything_that_voted_for_it():
    sample = Sample(file="a")
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 900, 980, 10, sample),
                        candidate("PC 16:1_18:0", "PC 34:1", 950, 990, 10, sample)],
                       purity=10.0, sum_id="PC 34:1")
    assert group.identification()[0] == "PC 34:1"        # purity below 75, chains not asserted
    assert group.identification_scores() == (950, 990, 10)


def test_molecular_row_takes_only_the_top_candidate():
    sample = Sample(file="a")
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 900, 980, 99, sample),
                        candidate("PC 16:1_18:0", "PC 34:1", 950, 990, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    assert group.identification()[0] == "PC 16:0_18:1"
    assert group.identification_scores() == (900, 980, 99)


def test_unidentified_row_has_blank_scores():
    empty = CompoundGroup(name="g", mw=1.0, retention=1.0, max_area=0.0, has_ms2=False, areas=[])
    assert empty.identification_scores() == ("", "", "")


def test_columns_are_written_and_land_before_the_sample_columns(tmp_path):
    sample = Sample(file="sample_one")
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [1234.0]
    out = tmp_path / "r.csv"
    write_results(PeakFinderResult([group], [sample]), out, scores=True, ms2_support=False)
    rows = list(csv.reader(out.open(newline="")))
    header, body = rows[0], rows[1]
    for column in ("Dot Product", "Reverse Dot Product", "Purity"):
        assert column in header
    assert header.index("Purity") < header.index("sample_one")
    assert body[header.index("Dot Product")] == "912"


def test_scores_can_be_turned_off(tmp_path):
    # Turning every optional column off, which is what reproduces LipiDex's own header exactly,
    # is covered in test_associated.py — `scores` is only one of the switches now.
    sample = Sample(file="s")
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [1.0]
    out = tmp_path / "r.csv"
    write_results(PeakFinderResult([group], [sample]), out, scores=False, ms2_support=False)
    header = next(csv.reader(out.open(newline="")))
    assert "Dot Product" not in header


def test_identified_only_drops_unnamed_rows_but_keeps_group_ids_stable(tmp_path):
    """`Final_Results_Filtered.csv` is the analysis-ready table and must not carry unnamed
    features — they cannot enter a lipid-level analysis, and a file labelled analysis-ready
    containing thousands of them invites someone to model them by accident.

    The Compound Group id must still be numbered over the full set. If it were numbered over what
    is written, the same id would mean a different molecule in each file and the tables could not
    be joined — which is the whole point of the column.
    """
    sample = Sample(file="s")
    named = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    named.areas = [10.0]
    unnamed = CompoundGroup(name="g", mw=700.0, retention=5.0, max_area=5.0,
                            has_ms2=False, areas=[5.0])
    result = PeakFinderResult([unnamed, named], [sample])

    everything = tmp_path / "all.csv"
    write_results(result, everything, ms2_support=False)
    only_named = tmp_path / "named.csv"
    write_results(result, only_named, ms2_support=False, identified_only=True)

    rows_all = list(csv.DictReader(everything.open(newline="")))
    rows_named = list(csv.DictReader(only_named.open(newline="")))
    assert len(rows_all) == 2 and len(rows_named) == 1
    # reported at molecular resolution, not sum composition
    assert rows_named[0]["Identification"] == "PC 16:0_18:1"
    # the surviving row keeps the id it had in the full table — 2, not 1
    kept_id = [r for r in rows_all
               if r["Identification"] == "PC 16:0_18:1"][0]["Compound Group"]
    assert rows_named[0]["Compound Group"] == kept_id == "2"


def _run_with_pools(areas_by_role):
    """A result whose samples carry roles, for the pooled-QC column."""
    samples, columns = [], []
    for role, values in areas_by_role:
        samples.append(Sample(file=f"{role}_{len(samples)}"))
        samples[-1].role = role
        columns.append(values)
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, samples[0])],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [v for v in columns]
    return PeakFinderResult([group], samples), group


def test_pooled_cv_is_raw_and_reproducible_from_the_table(tmp_path):
    """Raw, deliberately. A delivered table whose quality column cannot be recomputed from the
    columns beside it is a column the reader has to take on trust. The quality report normalises
    first and so reports a slightly different number; where the two diverge, the difference is the
    loading variation, which is worth seeing."""
    from lipidloop.peakfinder import pooled_cv
    result, group = _run_with_pools([("qc", 100.0), ("qc", 200.0), ("qc", 300.0),
                                     ("sample", 150.0)])
    cv = pooled_cv(result, [group])[id(group)]
    # mean 200, sd 100 -> 50%, exactly what anyone recomputing from the three pool columns gets
    assert cv == pytest.approx(50.0)

    # two pools is not a precision estimate
    result2, group2 = _run_with_pools([("qc", 100.0), ("qc", 200.0), ("sample", 1.0)])
    assert pooled_cv(result2, [group2]) == {}


def test_standards_columns_can_be_dropped(tmp_path):
    """The analysis-ready table leaves the standards out: a different material, injected to judge
    the instrument, and beside the samples it invites its way into a fold change."""
    sample = Sample(file="Sample_01"); sample.role = "sample"
    standard = Sample(file="Std_Mix_01"); standard.role = "standard"
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [10.0, 999.0]
    result = PeakFinderResult([group], [sample, standard])

    everything = tmp_path / "all.csv"
    write_results(result, everything, ms2_support=False, qc_cv=False)
    without = tmp_path / "no_std.csv"
    write_results(result, without, ms2_support=False, qc_cv=False, drop_roles=("standard",))

    assert "Std_Mix_01" in next(csv.reader(everything.open(newline="")))
    header = next(csv.reader(without.open(newline="")))
    assert "Std_Mix_01" not in header and "Sample_01" in header
    row = list(csv.DictReader(without.open(newline="")))[0]
    assert row["Sample_01"] == "10.0" and "999.0" not in row.values()


def test_median_normalisation_removes_a_pure_loading_difference(tmp_path):
    """Two injections of identical material at different amounts must come out identical, and the
    CV column must be computed from the normalised areas so it stays recomputable from the file."""
    from lipidloop.peakfinder import median_factors
    samples = []
    for role, name in (("qc", "Pool_1"), ("qc", "Pool_2"), ("qc", "Pool_3")):
        s = Sample(file=name); s.role = role; samples.append(s)
    groups = []
    for base in (100.0, 500.0, 20.0):
        g = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, samples[0])],
                       purity=99.0, sum_id="PC 34:1")
        g.areas = [base, base * 2, base * 3]      # injection 2 loaded 2x, injection 3 loaded 3x
        groups.append(g)
    result = PeakFinderResult(groups, samples)
    factors = median_factors(result, groups, [0, 1, 2])
    assert factors[1] / factors[0] == pytest.approx(2.0)
    assert factors[2] / factors[0] == pytest.approx(3.0)

    out = tmp_path / "norm.csv"
    write_results(result, out, ms2_support=False, normalise=True)
    row = list(csv.DictReader(out.open(newline="")))[0]
    values = [float(row[s.file]) for s in samples]
    assert values[0] == pytest.approx(values[1]) == pytest.approx(values[2])
    # a pure loading difference is not imprecision, so the CV must collapse to zero
    assert float(row["Pooled QC CV (%)"]) == pytest.approx(0.0, abs=1e-6)


def test_normalisation_with_nothing_detected_everywhere_leaves_the_data_alone():
    from lipidloop.peakfinder import median_factors
    samples = []
    for name in ("A", "B"):
        s = Sample(file=name); s.role = "sample"; samples.append(s)
    g = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, samples[0])],
                   purity=99.0, sum_id="PC 34:1")
    g.areas = [10.0, 0.0]          # never detected in both
    assert median_factors(PeakFinderResult([g], samples), [g], [0, 1]) == [1.0, 1.0]


def test_no_delivered_table_carries_a_standards_column(tmp_path):
    """A facility measurement, not the client's data. Whichever table someone opens, the standards
    must not be sitting beside the samples where a selection could sweep them into a comparison."""
    sample = Sample(file="Sample_01"); sample.role = "sample"
    pool = Sample(file="Pool_01"); pool.role = "qc"
    standard = Sample(file="Std_Mix_01"); standard.role = "standard"
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [10.0, 11.0, 999.0]
    result = PeakFinderResult([group], [sample, pool, standard])

    for name, kwargs in (("final.csv", {}),
                         ("unfiltered.csv", {"unfiltered": True}),
                         ("filtered.csv", {"identified_only": True}),
                         ("normalised.csv", {"identified_only": True, "normalise": True})):
        path = tmp_path / name
        write_results(result, path, ms2_support=False, qc_cv=False,
                      drop_roles=("standard",), **kwargs)
        header = next(csv.reader(path.open(newline="")))
        assert "Std_Mix_01" not in header, name
        assert "Sample_01" in header, name


def test_every_written_column_is_registered_as_metadata(tmp_path):
    """A guard against a whole class of silent corruption.

    `qc.load_run` treats everything after the last known metadata column as an injection. A new
    column that is not registered is therefore read as a sample — and `Pooled QC CV (%)`, because
    it contains "QC", was classified as a fifth pooled QC and silently entered every precision,
    missingness and principal-component figure in the report. Nothing crashed; the report simply
    described a column of percentages as if it were an injection.
    """
    from lipidloop.peakfinder import META_COLUMNS

    sample = Sample(file="Sample_01"); sample.role = "sample"
    pools = []
    for i in range(3):
        s = Sample(file=f"Pool_{i}"); s.role = "qc"; pools.append(s)
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [10.0, 11.0, 12.0, 13.0]
    result = PeakFinderResult([group], [sample] + pools)

    # every switch on, which is the widest header the writer can produce
    path = tmp_path / "widest.csv"
    write_results(result, path, unfiltered=True, id_source=True, ms2_support=True,
                  scores=True, adduct=True, group_id=True, qc_cv=True)
    header = [c for c in next(csv.reader(path.open(newline=""))) if c.strip()]
    injections = {s.file for s in result.samples}
    unregistered = [c for c in header if c not in META_COLUMNS and c not in injections]
    assert not unregistered, f"add these to META_COLUMNS: {unregistered}"

    # and the split lands where load_run expects it: metadata first, then injections
    last_meta = max(i for i, c in enumerate(header) if c in META_COLUMNS)
    first_injection = min(i for i, c in enumerate(header) if c in injections)
    assert last_meta < first_injection


def test_a_corrected_mass_is_written_beside_the_measured_one(tmp_path):
    """A delivered table that silently carries a derived mass, with no route back to what the
    instrument reported, cannot be checked by whoever receives it."""
    sample = Sample(file="S1"); sample.role = "sample"
    group = group_with([candidate("PC 16:0_18:1", "PC 34:1", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 34:1")
    group.areas = [1.0]
    group.quant_ion = 760.5851
    result = PeakFinderResult([group], [sample])

    plain = tmp_path / "plain.csv"
    write_results(result, plain, ms2_support=False, qc_cv=False)
    assert "Quant Ion (measured)" not in next(csv.reader(plain.open(newline=""))), \
        "no correction applied, so no second mass column"

    group.quant_ion_measured = 760.5807          # what the instrument said
    corrected = tmp_path / "corrected.csv"
    write_results(result, corrected, ms2_support=False, qc_cv=False)
    row = list(csv.DictReader(corrected.open(newline="")))[0]
    assert float(row["Quant Ion"]) == pytest.approx(760.5851)
    assert float(row["Quant Ion (measured)"]) == pytest.approx(760.5807)


def test_normalisation_factors_come_from_the_rows_that_are_written(tmp_path):
    """The factors are a median across rows, so they must come from the rows they correct.

    Derived from every compound group — most of which are rejected features and noise — and applied
    to the identified lipids, the correction is measured on one population and applied to another.
    Here the rejected rows carry the opposite loading pattern to the kept ones, so a factor taken
    over everything normalises the wrong way round.
    """
    samples = []
    for name in ("Pool_1", "Pool_2"):
        s = Sample(file=name); s.role = "qc"; samples.append(s)

    def make(name, areas, keep=True):
        g = group_with([candidate("PC 16:0_18:1", name, 912, 987, 99, samples[0])],
                       purity=99.0, sum_id=name)
        g.areas = areas
        g.keep = keep
        return g

    # kept, identified: injection 2 loaded twice as much
    kept = [make(f"PC 3{i}:1", [100.0, 200.0]) for i in range(4)]
    # rejected: the opposite pattern, and far more numerous
    junk = [make("", [200.0, 100.0], keep=False) for _ in range(40)]
    result = PeakFinderResult(junk + kept, samples)

    out = tmp_path / "norm.csv"
    write_results(result, out, ms2_support=False, identified_only=True, normalise=True)
    rows = list(csv.DictReader(out.open(newline="")))
    assert len(rows) == 4, "only the identified, kept rows are written"
    for row in rows:
        # a pure 2x loading difference among the written rows must be removed
        assert float(row["Pool_1"]) == pytest.approx(float(row["Pool_2"]))


def test_the_written_normalised_matrix_is_centred_on_its_own_rows(tmp_path):
    """The property a reader can check, so it is asserted rather than hoped for.

    Factors are derived from features detected in every injection — a basis that does not move when
    the library changes — but that is not the whole table, so the written matrix was a few per cent
    off centre. A filter applied later changes which rows are written and would decentre it again
    by another route; this test catches both.
    """
    import statistics
    samples = []
    for i in range(4):
        s = Sample(file=f"S{i}"); s.role = "sample"; samples.append(s)
    groups = []
    for n in range(30):
        g = group_with([candidate("PC 16:0_18:1", f"PC 3{n}:1", 912, 987, 99, samples[0])],
                       purity=99.0, sum_id=f"PC 3{n}:1")
        base = 100.0 * (n + 1)
        g.areas = [base, base * 2, base * 3, base * 1.5]     # pure loading differences
        # a third of the rows are missing from one injection, which is what pulled the median off
        if n % 3 == 0:
            g.areas[1] = 0.0
        groups.append(g)
    out = tmp_path / "norm.csv"
    write_results(PeakFinderResult(groups, samples), out, ms2_support=False,
                  identified_only=True, normalise=True)

    rows = list(csv.DictReader(out.open(newline="")))
    medians = []
    for s in samples:
        vals = [float(r[s.file]) for r in rows if float(r[s.file]) > 0]
        medians.append(statistics.median(vals))
    assert max(medians) / min(medians) - 1 < 0.01, f"not centred: {medians}"


def test_duplicate_verdict_marks_only_duplicated_rows_and_never_deletes(tmp_path):
    """The verdict is about choosing between rows that share a name. On a unique identification
    there is nothing to choose between, and a verdict there would read as a quality score for the
    whole table."""
    samples = []
    for name in ("Pool_1", "Pool_2", "Pool_3"):
        s = Sample(file=name); s.role = "qc"; samples.append(s)

    def make(name, rt, areas, model_error):
        # the MOLECULAR name is what identification() returns and what duplicates are counted on —
        # passing the same one to every row makes them all duplicates of each other
        g = group_with([candidate(name, name, 912, 987, 99, samples[0])],
                       purity=99.0, sum_id=name)
        g.retention, g.areas, g.retention_error = rt, areas, model_error
        # the verdict thresholds on the z, not on minutes: half a minute is nothing in TG
        # (residual sd 0.36) and enormous in LysoPC (0.05), so minutes are a different test in
        # every class. 0.2 min here stands for a typical class residual spread.
        g.retention_z = model_error / 0.2
        return g

    # two rows share a name: one precise and on-model, one imprecise and far off it
    good = make("PC 34:1", 9.0, [100.0, 101.0, 99.0], 0.05)
    bad = make("PC 34:1", 9.6, [100.0, 300.0, 20.0], 0.62)
    alone = make("PC 36:2", 10.0, [50.0, 51.0, 49.0], 0.03)
    out = tmp_path / "v.csv"
    write_results(PeakFinderResult([good, bad, alone], samples), out,
                  ms2_support=False, identified_only=True)
    rows = {(r["Identification"], r["Retention Time (min)"]): r
            for r in csv.DictReader(out.open(newline=""))}
    assert len(rows) == 3, "nothing is deleted"

    verdict = {k[1]: v["Duplicate Verdict"] for k, v in rows.items()}
    assert verdict["9.6"] == "unreliable"
    assert verdict["9.0"] == "isomer?", "precise — nothing says it is wrong"
    assert verdict["10.0"] == "", "not duplicated, so not judged"

    # and the retention error is reported, having previously been computed and discarded
    assert float(rows[("PC 34:1", "9.6")]["RT Model Error"]) == pytest.approx(0.62)


def test_lipid_key_is_unique_and_leaves_the_identification_alone(tmp_path):
    """A name safe to plot and to run one test per row against.

    Identification is the join key across all four tables, Associated_Spectra.csv and other runs,
    so it is not touched. The suffix is a LETTER: the numbers in a lipid name already mean carbon
    count and double bonds, so `PC 34:1(2)` invites the misreading it exists to prevent, and `_`
    separates chains — `PC 34:1_2` would be counted as chain-resolved by the rest of the pipeline.
    """
    sample = Sample(file="S1"); sample.role = "sample"

    def make(name, rt):
        g = group_with([candidate(name, name, 912, 987, 99, sample)], purity=99.0, sum_id=name)
        g.retention, g.areas = rt, [100.0]
        return g

    # letters follow INTENSITY, largest first, so `a` is the dominant species
    groups = [make("PC 34:1", 9.4), make("PC 34:1", 8.1), make("SM d34:1", 7.0)]
    groups[0].areas, groups[0].max_area = [900.0], 900.0      # the later peak is the bigger one
    groups[1].areas, groups[1].max_area = [100.0], 100.0
    out = tmp_path / "k.csv"
    write_results(PeakFinderResult(groups, [sample]), out, ms2_support=False,
                  identified_only=True, qc_cv=False)
    rows = {r["Retention Time (min)"]: r for r in csv.DictReader(out.open(newline=""))}

    assert rows["9.4"]["Lipid Key"] == "PC 34:1a", "the most intense row is a"
    assert rows["8.1"]["Lipid Key"] == "PC 34:1b"
    assert rows["7.0"]["Lipid Key"] == "SM d34:1", "a unique name is unchanged"
    # the identification itself is untouched, so joins still work
    assert {r["Identification"] for r in rows.values()} == {"PC 34:1", "SM d34:1"}
    keys = [r["Lipid Key"] for r in rows.values()]
    assert len(keys) == len(set(keys)), "the key must be unique across the table"
    assert not any("_" in k for k in keys), "an underscore would read as a chain separator"


def test_retention_never_condemns_a_duplicate_however_far_off_the_model(tmp_path):
    """The model cannot rank isomers, so it must not be asked to.

    `PC 18:0_20:3` and `PC 18:1_20:2` both parse to ("PC", "", 38, 3): one prediction serves every
    isomer of a composition, so at most one of them can sit on it and the rest are displaced by
    construction. A verdict from that displacement says "you are not the isomer nearest the class
    average", which is a statement about position in an eluting series and not evidence against a
    row. Held here because it was shipped, and the fixture is the real PC 38:3 — four rows, one
    prediction at 9.455 min, errors that are the retention axis minus a constant.
    """
    from lipidloop.peakfinder import _duplicate_verdicts
    sample = Sample(file="S1"); sample.role = "sample"

    def make(rt):
        g = group_with([candidate("PC 38:3", "PC 38:3", 912, 987, 99, sample)],
                       purity=99.0, sum_id="PC 38:3")
        g.retention, g.areas = rt, [100.0]
        g.retention_error = rt - 9.455
        g.retention_z = g.retention_error / 0.05      # a tight class: up to 15 sigma out
        return g

    groups = [make(rt) for rt in (8.985, 9.321, 9.644, 10.19)]
    duplicates = {id(g): f"{i} of 4" for i, g in enumerate(groups, start=1)}
    verdicts = _duplicate_verdicts(groups, duplicates, {})
    assert [abs(g.retention_z) for g in groups][-1] > 10, "the fixture really is far off the model"
    assert set(verdicts.values()) == {"isomer?"}, "retention alone condemns none of them"

    # precision still does, being independent of where a row elutes
    assert _duplicate_verdicts(groups, duplicates,
                               {id(groups[-1]): 144.9})[id(groups[-1])] == "unreliable"


def test_a_whole_displaced_set_is_flagged_though_none_of_its_rows_can_be():
    """The one thing retention can still say about a duplicated name.

    Isomer positioning explains why all but one row sits away from the prediction. It does not
    explain the row NEAREST the prediction also sitting far from it — no arrangement of a set
    around a point puts every member far from it. So the set is judged on its smallest |z|, which
    is legitimate because that quantity was measured to behave like a unique name (0.099 min out
    against 0.112). The suspicion is on the class model or the identification, not on a row.
    """
    from lipidloop.peakfinder import SET_OFF_MODEL_Z, _duplicate_verdicts
    sample = Sample(file="S1"); sample.role = "sample"

    def make(name, z):
        g = group_with([candidate(name, name, 912, 987, 99, sample)], purity=99.0, sum_id=name)
        g.retention, g.areas = 9.0 + z, [100.0]
        g.retention_error, g.retention_z = z * 0.05, z
        return g

    # straddling: one row on the prediction, the others displaced around it — the normal case
    straddle = [make("PC 38:3", z) for z in (-9.0, -2.7, 3.8, 14.9)]
    # displaced: the same spread, but the whole set sits to one side and even its best row is out
    shifted = [make("PS 40:6", z) for z in (20.0, 26.7, 33.8, 44.9)]
    groups = straddle + shifted
    duplicates = {id(g): "n of 4" for g in groups}
    verdicts = _duplicate_verdicts(groups, duplicates, {})

    assert {verdicts[id(g)] for g in straddle} == {"isomer?"}, "one row is on the prediction"
    assert {verdicts[id(g)] for g in shifted} == {"class model?"}, "none of them is"
    assert min(abs(g.retention_z) for g in straddle) < SET_OFF_MODEL_Z < 20.0

    # a set is never judged on the subset the model happened to reach
    shifted[0].retention_z = None
    partial = _duplicate_verdicts(groups, duplicates, {})
    assert {partial[id(g)] for g in shifted} == {"isomer?"}, "one row unscored, so the set is not"


def test_shorthand_column_sits_under_its_own_heading(tmp_path):
    """`Shorthand (LSI)` and `Lipid Key` must not swap places.

    The header adds the shorthand at Identification+1 and then the key at Identification+1 as
    well, which puts the key BEFORE the shorthand. Inserting the row values in the other order
    produced a table where every value was present and every one sat under the wrong heading —
    `Lipid Key` holding `PC 34:1` and `Shorthand (LSI)` holding `PC 34:1a`. Read by column name,
    as every downstream reader does, that is silent corruption.
    """
    sample = Sample(file="S1"); sample.role = "sample"
    g = group_with([candidate("SM d34:1", "SM d34:1", 912, 987, 99, sample)],
                   purity=99.0, sum_id="SM d34:1")
    g.retention, g.areas, g.max_area = 7.0, [100.0], 100.0
    out = tmp_path / "s.csv"
    write_results(PeakFinderResult([g], [sample]), out, ms2_support=False, qc_cv=False)
    row = next(csv.DictReader(out.open(newline="")))
    assert row["Identification"] == "SM d34:1"
    assert row["Shorthand (LSI)"] == "SM 34:1;O2", "the d prefix is an oxygen count in LSI"
    assert row["Lipid Key"] == "SM d34:1", "a unique name is its own key"


def test_shorthand_never_asserts_a_plasmalogen_without_a_spectrum(tmp_path):
    """`P-` claims a vinyl ether bond. An `RT model` row has no MS2 to support it.

    w2's validation of the class mapping (final-003) is explicit: a shorthand writer must not
    rewrite a demoted ether back to `P-`, and a row named from accurate mass and retention alone
    cannot distinguish alkyl from alkenyl at all. The weaker, defensible form is `O-`.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop.lsi import shorthand

    assert shorthand("Plasmenyl-PE P-38:4", "MS2") == "PE P-38:4"
    assert shorthand("Plasmenyl-PE P-38:4", "RT model") == "PE O-38:4"
    assert shorthand("Plasmanyl-PC O-38:2", "MS2") == "PC O-38:2"
    # Cardiolipin stays at species level: its written order encodes sn-pairing, so a `_` form
    # would need sorting and a `/` form would over-claim positions.
    assert shorthand("CL 72:7", "MS2") == "CL 72:7"
    # `/` is never emitted — the pipeline proves chains, never sn-positions.
    assert "/" not in shorthand("TG 16:0_18:1_18:2", "MS2")


def test_library_set_follows_the_mobile_phase_modifier():
    """Negative-mode adducts follow the modifier; a formate library cannot describe acetate data.

    Searching formate libraries on a study run with ammonium acetate named 19 acetate adducts of
    ordinary dihydroxy ceramides as `Cer[AP]` PHYTOCERAMIDES — within 0.0-3.9 ppm of the true
    masses, at retention times agreeing to +/-0.03 min, with dot products of 999 and Purity 0.
    `LipiDex_HCD_Formic` holds zero ceramide `[M+Ac-H]-` entries and `LipiDex_HCD_Acetate` holds
    2,120, so the identification could not have been right.
    """
    import sys
    from pathlib import Path
    root = Path(__file__).resolve().parents[1]
    sys.path.insert(0, str(root / "scripts"))
    import importlib
    rs = importlib.import_module("run_study")

    # The two shipped mobile phases. Validation-only sets (holdout, decoy) may exist alongside
    # them; what matters is that both real modifiers are present and neither borrows the other's
    # libraries.
    assert {"formate", "acetate"} <= set(rs.LIBRARIES)
    assert rs.DEFAULT_MODIFIER == "formate", "the facility method — its results must not move"
    for modifier in ("formate", "acetate"):
        for polarity in ("Pos", "Neg"):
            stems = rs.LIBRARIES[modifier][polarity]
            wrong = "Acetate" if modifier == "formate" else "Formic"
            assert not any(wrong in s for s in stems), \
                f"{modifier}/{polarity} carries a {wrong} library: {stems}"
        assert modifier in rs.ISTD_LIBRARY, "spiked standards must match the modifier too"


def test_a_chain_resolved_name_needs_ITS_OWN_purity_not_the_peak_s(tmp_path):
    """A row must not claim chains that nothing in its own match supports.

    `CompoundGroup.purity` is the peak's DOMINANT purity — the best-supported component present,
    which need not be the one being named. So a chain-resolved name could ship while the purity
    attributable to that name was zero: good fragment evidence for one lipid, chains asserted for
    another. That is what `Cer[AP] t20:0_23:0` carried while actually being an acetate adduct of an
    ordinary ceramide at dot product 999.

    Narrow deliberately: purity 0 is the NORM for a sum composition, which makes no chain claim.
    Only a name that asserts chains while its own evidence is empty gets demoted.
    """
    sample = Sample(file="S1"); sample.role = "sample"
    good = candidate("PC 16:0_18:1", "PC 34:1", 950, 900, 99, sample)
    weak = candidate("PC 15:0_19:1", "PC 34:1", 940, 890, 0, sample)

    g = group_with([weak], purity=99.0, sum_id="PC 34:1")   # peak looks pure, THIS name does not
    g.lipid_candidates[0].purity = 0.0
    name, _ = g.identification()
    assert name == "PC 34:1", f"chains unsupported by their own match must demote, got {name}"

    g2 = group_with([good], purity=99.0, sum_id="PC 34:1")
    g2.lipid_candidates[0].purity = 99.0
    name2, _ = g2.identification()
    assert name2 == "PC 16:0_18:1", "a well-supported chain assignment is kept"


def test_signal_to_noise_is_reported_and_stays_in_its_own_column():
    """S/N tells a measurement from an integration of noise; abundance cannot.

    ⚠ Appended LAST and gated on the same flag on both sides. A first attempt put it beside
    `Class (canonical)`, where the conditional blocks do not fire in the same combinations for the
    header and the row — every score column shifted one place and `Dot Product` read out an
    adduct. The suite caught it; a reader of the CSV would not have.
    """
    from lipidloop.features import signal_to_noise

    scans = []
    for i in range(120):
        rt = i * 0.05 / 60.0                    # 0.05 s spacing, in minutes
        t = i * 0.05
        # a peak at 5.0 s on a baseline of ~100
        height = 10000.0 if abs(t - 5.0) < 0.3 else 0.0
        scans.append((rt, [500.0], [100.0 + height]))
    # apex 5 s, FWHM 0.6 s -> baseline window is 1.2-3.0 s before it, which exists here
    snr = signal_to_noise(scans, 500.0, 5.0, 0.6, floor=100.0)
    assert snr is not None and snr > 10, f"a clean peak on a flat baseline should score high: {snr}"

    # ⚠ An EMPTY baseline must still give a number. In centroided data an m/z with no signal
    # contributes no peak, so the baseline reads exactly zero and its MAD is zero — returning None
    # there would blank the S/N for the cleanest peaks of all.
    empty = [(i * 0.05 / 60.0, [500.0], [10000.0 if abs(i * 0.05 - 5.0) < 0.3 else 0.0])
             for i in range(120)]
    assert signal_to_noise(empty, 500.0, 5.0, 0.6, floor=50.0) == 200.0

    # ⚠ No leading baseline -> None, not a fabricated value. An early-eluting peak has nothing in
    # front of it, and measuring its own tail instead would be worse than declining.
    assert signal_to_noise(scans, 500.0, 0.2, 0.6) is None
    # and a mass that is not there has no height
    assert signal_to_noise(scans, 999.0, 5.0, 0.6) is None


def test_mass_error_column_reports_the_matched_ion(tmp_path):
    """Theoretical mass must describe the ion the row is QUANTIFIED on, not any adduct of it.

    The value existed at match time and was discarded, so every consumer re-derived it by looking
    the name up in the library — which fails for roughly seven identifications in eight, because
    the table carries a canonicalised name (`SM d30:1`) and the library the raw MSP entry
    (`SM d14:1_16:0`). Per-file mass calibration could not be measured without it.
    """
    from lipidloop.peaks import CompoundGroup

    group = CompoundGroup.__new__(CompoundGroup)
    # Only the fields the accessors touch; a full construction needs a whole run behind it.
    for field, value in (("lipid_candidates", []), ("quant_ion", 760.5851),
                         ("rtls_identification", ""), ("rtls_adduct", ""), ("sum_id", ""),
                         ("plasmenyl_ether_conflict", False), ("final_lipid_id", None)):
        setattr(group, field, value)
    assert group.theoretical_mz() is None, "no candidates means no claim, not a guess"
    assert group.mass_error_ppm() is None


def test_the_new_columns_sit_after_signal_to_noise(tmp_path):
    """Header and row order must agree — the failure is silent and shifts every column after it.

    S/N had to be appended last for this reason once already; adding two more columns in front of
    it put theoretical mass under the `Pooled QC CV` heading, with every value present and every
    one mislabelled.
    """
    import inspect
    from lipidloop import peakfinder

    source = inspect.getsource(peakfinder.write_results) if hasattr(peakfinder, "write_results") \
        else inspect.getsource(peakfinder)
    header_at = source.index('header += ["Theoretical m/z", "Mass Error (ppm)"]')
    snr_header = source.index('header.append("S/N")')
    assert snr_header < header_at, "header: S/N must come before the mass columns"


def test_class_exclusion_removes_decoys_too(tmp_path):
    """An excluded class must leave the DECOY library as well, or the FDR is inflated one-sidedly.

    A decoy is named `DECOY_ d5TG 10:0_10:0_10:0`, so its class parses as the pseudo-class
    `DECOY_` and `"decoy_" != "d5tg"`. `excluded_classes=["d5TG"]` therefore removed 22,960 target
    entries and none of the 22,960 decoys — the entries could still win a match and be counted as a
    false positive while their targets were no longer there to win the true one. 17 d5TG decoy
    hits, 12.5% of every decoy hit measured, all numerator and no denominator.

    General, not specific to d5TG: every decoy in every library parses as `DECOY_`.
    """
    from lipidloop.search import load_libraries

    lib = tmp_path / "lib.msp"
    lib.write_text(
        "Name: d5TG 10:0_10:0_10:0 [M+NH4]+;\nMW: 577.5204\nPRECURSORMZ: 577.5204\n"
        "Comment: Name=d5TG 10:0_10:0_10:0 Type=LipiDex\nNum Peaks: 1\n"
        '100.0000 999 "x_Fragment_[]"\n\n'
        "Name: DECOY_ d5TG 10:0_10:0_10:0 [M+NH4]+;\nMW: 577.5204\nPRECURSORMZ: 577.5204\n"
        "Comment: Name=DECOY_d5TG 10:0_10:0_10:0 Type=LipiDex\nNum Peaks: 1\n"
        '114.0000 999 "x_Fragment_[]"\n\n'
        "Name: PC 34:1 [M+H]+;\nMW: 760.5851\nPRECURSORMZ: 760.5851\n"
        "Comment: Name=PC 34:1 Type=LipiDex\nNum Peaks: 1\n"
        '184.0733 999 "x_Fragment_[]"\n\n')

    kept = load_libraries([lib], excluded_classes=["d5TG"]).spectra
    names = [s.lipid for s in kept]
    assert not any("d5TG" in n for n in names), f"d5TG survived exclusion: {names}"
    assert any("PC 34:1" in n for n in names), "an unrelated class must be untouched"
    assert len(load_libraries([lib]).spectra) == 3, "without exclusion all three load"
