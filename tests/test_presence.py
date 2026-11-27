"""Presence filtering.

The rule is "detected in a fraction of at least one group", not "of all samples", because a lipid
present in one group and absent from another is the strongest result an experiment can produce and
an all-samples rule deletes it. Filling those gaps is gapfill.py's job.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample   # noqa: E402
from lipidloop.presence import (PresenceFilter, apply_presence_filter,   # noqa: E402
                                 presence_cells, read_metadata)


def metadata_file(tmp_path, rows):
    path = tmp_path / "metadata.csv"
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["sample", "age", "genotype"])
        writer.writerows(rows)
    return path


def design(tmp_path):
    """Four groups of two: 3m/12m by WT/ALS, plus a QC and a blank that must be ignored."""
    files, rows = [], []
    for age in ("3m", "12m"):
        for genotype in ("WT", "ALS"):
            for n in range(2):
                name = f"{age}_{genotype}_{n}"
                files.append(Sample(file=name, role="sample"))
                rows.append([name, age, genotype])
    samples = files + [Sample(file="QC1", role="qc"), Sample(file="Blank1", role="blank")]
    return samples, metadata_file(tmp_path, rows)


def group(areas, name="PC 16:0_18:1"):
    out = CompoundGroup(name=name, mw=759.6, retention=5.0, max_area=max(areas + [0.0]),
                        has_ms2=True, areas=list(areas))
    out.quant_ion = 760.585
    lipid = Lipid(retention=5.0, precursor=760.585, sample=Sample(file="x"), dot=900,
                  rev_dot=980, lipid_string=f"{name} [M+H]+;", lib_precursor=760.585, purity=99,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    out.lipid_candidates = [LipidCandidate(lipid)]
    out.final_lipid_id = out.lipid_candidates[0]
    out.purity, out.sum_id = 99.0, name
    return out


def params(tmp_path, **kw):
    """⚠ `min_group_size=2` here, against the production default of 4. These fixtures use cells of
    two so a 70% rule can be exercised in both directions in a readable table; production wants
    n >= 3 before it will call something a group. The floor itself is tested separately."""
    _, meta = design(tmp_path)
    kw.setdefault("min_group_size", 2)
    return PresenceFilter(enabled=True, min_fraction=0.7, metadata=str(meta),
                          group_by=["age", "genotype"], **kw)


# ── the filter ─────────────────────────────────────────────────────────────────────────────

def test_a_lipid_confined_to_one_group_is_kept(tmp_path):
    """The whole point. Present in both members of 3m/WT, absent everywhere else — an all-samples
    rule would drop it at 25% detection, and it is the strongest signal in the experiment."""
    samples, _ = design(tmp_path)
    g = group([10.0, 12.0] + [0.0] * 6 + [0.0, 0.0])
    report = apply_presence_filter([g], samples, params(tmp_path))
    assert g.keep and report.kept == 1 and report.identified_kept == 1


def test_a_lipid_scattered_across_groups_is_dropped(tmp_path):
    """Detected in half of every group: 50% everywhere, 70% nowhere. Same total detection count
    as the case above and the opposite verdict, which is the rule working."""
    samples, _ = design(tmp_path)
    g = group([10.0, 0.0, 11.0, 0.0, 12.0, 0.0, 13.0, 0.0] + [0.0, 0.0])
    report = apply_presence_filter([g], samples, params(tmp_path))
    assert not g.keep and report.dropped == 1
    assert g.filter_reason.startswith("Not detected in 70% of any group")


def test_qc_and_blank_columns_do_not_count_towards_a_group(tmp_path):
    samples, _ = design(tmp_path)
    g = group([10.0] + [0.0] * 7 + [99.0, 99.0])     # one sample, but both QC and blank
    apply_presence_filter([g], samples, params(tmp_path))
    assert not g.keep


def test_already_filtered_rows_are_left_alone(tmp_path):
    samples, _ = design(tmp_path)
    g = group([0.0] * 10)
    g.keep, g.filter_reason = False, "Blank"
    apply_presence_filter([g], samples, params(tmp_path))
    assert g.filter_reason == "Blank"       # not overwritten by the presence reason


def test_without_metadata_it_degrades_to_all_samples_and_says_so(tmp_path):
    samples, _ = design(tmp_path)
    g = group([10.0, 12.0] + [0.0] * 6 + [0.0, 0.0])
    report = apply_presence_filter([g], samples,
                                   PresenceFilter(enabled=True, min_fraction=0.7))
    assert report.ungrouped and not g.keep    # 25% of all samples
    assert "no metadata" in str(report)


def test_metadata_is_read_on_the_first_column(tmp_path):
    _, meta = design(tmp_path)
    table = read_metadata(meta, ["age", "genotype"])
    assert table["3m_WT_0"] == ("3m", "WT")
    assert len(table) == 8


def test_qc_becomes_the_group_when_a_study_has_no_samples():
    """A deposit whose every injection is a pool must still be filtered.

    Before this, `_groups` skipped anything not role `sample`, so an all-QC study produced no cells
    and the filter was skipped entirely — leaving those tables unfiltered while every other study
    in the same comparison had been filtered. Three of the validation deposits are like this
    (MTBKS222_Agilent, ST000991_DDA, ST003077): the pool is the specimen, not a control.
    """
    samples = [Sample(file=f"q{i}", role="qc") for i in range(4)]
    groups = [group([1.0, 1.0, 1.0, 1.0]), group([1.0, 0.0, 0.0, 0.0])]
    report = apply_presence_filter(groups, samples,
                                   PresenceFilter(enabled=True, min_fraction=0.7))
    assert (report.kept, report.dropped) == (1, 1)
    assert "qc=4" in str(report)


def test_qc_as_a_group_can_only_rescue_never_drop():
    """QC is a group like any other, and the rule is a union — so a lipid carried only by the
    samples still passes on the sample cell, and one carried only by the pools passes on the QC
    cell. Neither can delete the other; that is what separates this from the QC rule the module
    warns about, which REQUIRES presence in the pools."""
    samples = ([Sample(file=f"s{i}", role="sample") for i in range(4)]
               + [Sample(file=f"q{i}", role="qc") for i in range(4)])
    par = PresenceFilter(enabled=True, min_fraction=0.7)
    sample_only = [group([1.0] * 4 + [0.0] * 4)]
    assert apply_presence_filter(sample_only, samples, par).kept == 1
    qc_only = [group([0.0] * 4 + [1.0] * 4)]
    assert apply_presence_filter(qc_only, samples, par).kept == 1
    neither = [group([1.0, 0.0, 0.0, 0.0] * 2)]
    assert apply_presence_filter(neither, samples, par).dropped == 1


def test_standards_and_blanks_are_not_groups():
    """A feature in every Std_Mix injection is the standard, not evidence about the study; a blank
    is the thing being filtered against. Neither forms a cell."""
    samples = ([Sample(file=f"s{i}", role="sample") for i in range(4)]
               + [Sample(file=f"std{i}", role="standard") for i in range(4)]
               + [Sample(file=f"b{i}", role="blank") for i in range(4)])
    report = apply_presence_filter([group([0.0] * 4 + [1.0] * 8)], samples,
                                   PresenceFilter(enabled=True, min_fraction=0.7))
    assert report.dropped == 1
    assert "standard" not in str(report) and "blank" not in str(report)


def test_a_cell_below_the_floor_is_not_a_group():
    """1/1 is 100%, so a singleton cell passes anything detected in it once — the fraction rule
    says nothing there. A pooled-tissue atlas with n=1 per tissue would have had the filter
    silently do nothing; now it reports that it found no group large enough and skips.
    """
    samples = [Sample(file=n, role="sample") for n in ("liver", "heart", "kidney")]
    meta = {n: (n,) for n in ("liver", "heart", "kidney")}   # n=1 per tissue
    from lipidloop.presence import _groups, _sized
    assert len(_groups(samples, meta)) == 3
    assert _sized(_groups(samples, meta), 3) == {}
    # And with every cell below the floor the filter skips rather than dropping anything: two
    # injections are not a group, so there is nothing to judge against.
    pair = [Sample(file="a", role="sample"), Sample(file="b", role="sample")]
    report = apply_presence_filter([group([1.0, 0.0])], pair,
                                   PresenceFilter(enabled=True, metadata="", group_by=[]))
    assert report.kept == 0 and report.dropped == 0


def test_two_of_three_survives_at_the_default_threshold():
    """The reason the default is 0.6 rather than 0.7.

    ⚠ At 0.7 a group of three required all three, because 2/3 = 0.667 falls short — so a triplicate
    arm was filtered harder than a twenty-injection one under what looked like the same rule. Most
    of this validation set is small-n, so that quirk was not an edge case, it was the common case.
    """
    samples = [Sample(file=f"s{i}", role="sample") for i in range(3)]
    par = PresenceFilter(enabled=True)                      # default min_fraction
    assert par.min_fraction == 0.6
    two_of_three = [group([1.0, 1.0, 0.0])]
    assert apply_presence_filter(two_of_three, samples, par).kept == 1
    one_of_three = [group([1.0, 0.0, 0.0])]
    assert apply_presence_filter(one_of_three, samples, par).dropped == 1
    # and at the old threshold the same feature would have been dropped
    strict = PresenceFilter(enabled=True, min_fraction=0.7)
    assert apply_presence_filter([group([1.0, 1.0, 0.0])], samples, strict).dropped == 1
