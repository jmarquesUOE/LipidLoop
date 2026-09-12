"""Tests for sample roles and the blank filter.

The filter's whole contract is that it never touches an intensity — it only decides whether a
row survives. Most of these exist to pin that.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.blanks import (BlankFilter, apply_blank_filter, carryover_report,  # noqa: E402
                               infer_role, read_sequence)
from lipidloop.peaks import CompoundGroup, Sample                                  # noqa: E402


def group(areas, name=""):
    g = CompoundGroup(name="", mw=700.0, retention=8.0, max_area=max(areas),
                      has_ms2=False, areas=list(areas))
    if name:
        from lipidloop.peaks import LipidCandidate, Lipid
        lipid = Lipid(retention=8.0, precursor=700.0, sample=Sample(file="x"), dot=900.0,
                      rev_dot=990.0, lipid_string=f"{name} [M+H]+;", lib_precursor=700.0,
                      purity=90, is_lipidex=True, purity_array=[], fragment_masses=[])
        g.lipid_candidates = [LipidCandidate(lipid)]
        g.sum_id = lipid.sum_lipid_name
        g.purity = 90.0
    return g


def samples(roles):
    return [Sample(file=f"f{i}", role=r) for i, r in enumerate(roles)]


ROLES = ["sample", "sample", "sample", "blank"]


def test_role_inference():
    assert infer_role("Blank_01") == "blank"
    assert infer_role("QC_03") == "qc"
    assert infer_role("Pool_02") == "qc"
    assert infer_role("20260602_U5_Heart_12_wkCKD") == "sample"
    assert infer_role("SHAM_1_Heart_12_wk") == "sample"


def test_group_well_above_the_blank_survives():
    groups = [group([1000.0, 1100.0, 900.0, 10.0])]
    report = apply_blank_filter(groups, samples(ROLES), BlankFilter(multiplier=5.0))
    assert groups[0].keep
    assert report.kept == 1 and report.dropped == 0


def test_group_at_blank_level_is_dropped_whole():
    groups = [group([100.0, 110.0, 90.0, 100.0])]
    apply_blank_filter(groups, samples(ROLES), BlankFilter(multiplier=5.0))
    assert not groups[0].keep
    assert groups[0].filter_reason == "Blank"


def test_intensities_are_never_modified():
    """The filter decides; it does not subtract. This is the whole point."""
    original = [1000.0, 1100.0, 900.0, 10.0]
    groups = [group(original)]
    apply_blank_filter(groups, samples(ROLES), BlankFilter())
    assert groups[0].areas == original

    dropped = [100.0, 110.0, 90.0, 100.0]
    groups = [group(dropped)]
    apply_blank_filter(groups, samples(ROLES), BlankFilter())
    assert groups[0].areas == dropped      # dropped, but untouched


def test_multiplier_moves_the_boundary():
    areas = [400.0, 400.0, 400.0, 100.0]        # exactly 4x
    for multiplier, expect_keep in [(3.0, True), (5.0, False)]:
        groups = [group(areas)]
        apply_blank_filter(groups, samples(ROLES), BlankFilter(multiplier=multiplier))
        assert groups[0].keep is expect_keep


def test_max_statistic_is_stricter_than_mean():
    roles = ["sample", "sample", "blank", "blank"]
    areas = [500.0, 500.0, 10.0, 190.0]         # mean blank 100, max blank 190
    kept = {}
    for statistic in ("mean", "max"):
        groups = [group(areas)]
        apply_blank_filter(groups, samples(roles),
                           BlankFilter(multiplier=5.0, statistic=statistic))
        kept[statistic] = groups[0].keep
    assert kept["mean"] and not kept["max"]


def test_a_group_seen_in_too_few_samples_is_dropped():
    groups = [group([1000.0, 0.0, 0.0, 0.0])]
    apply_blank_filter(groups, samples(ROLES), BlankFilter(require_detected_in=2))
    assert not groups[0].keep


def test_no_blanks_means_no_filtering():
    groups = [group([1.0, 1.0, 1.0, 1.0])]
    report = apply_blank_filter(groups, samples(["sample"] * 4), BlankFilter())
    assert groups[0].keep and report.dropped == 0


def test_disabled_filter_keeps_everything():
    groups = [group([100.0, 100.0, 100.0, 100.0])]
    apply_blank_filter(groups, samples(ROLES), BlankFilter(enabled=False))
    assert groups[0].keep


def test_carryover_report_separates_blanks():
    roles = ["sample", "sample", "blank", "blank"]
    groups = [group([1000.0, 1000.0, 5.0, 500.0]),
              group([2000.0, 2000.0, 0.0, 900.0])]
    rows = carryover_report(groups, samples(roles))
    assert len(rows) == 2
    early, late = rows
    assert late["total signal"] > early["total signal"]   # the trailing blank carries more


def test_sequence_file_gives_injection_order(tmp_path):
    path = tmp_path / "seq.csv"
    path.write_text(
        "Bracket Type=4,,,\n"
        "Sample Type,File Name,Sample ID,Path\n"
        "Unknown,Blank_01,Inj_01,D:\\Pos\n"
        "Unknown,QC_01,Inj_03,D:\\Pos\n"
        "Unknown,Blank_02,Inj_51,D:\\Pos\n")
    order = read_sequence(path)
    assert order == {"Blank_01": 1, "QC_01": 3, "Blank_02": 51}


def test_injection_order_comes_from_row_order_not_sample_id(tmp_path):
    """Sample ID is a vial position in two of the three real sequences seen so far.

    Taking its trailing digits gave all five Kiterie blanks injection 8 — they share a vial, so
    they share a position — and the carryover report, whose whole job is to tell a leading blank
    from a trailing one, became meaningless.
    """
    from lipidloop.blanks import read_sequence

    positions = tmp_path / "positions.csv"
    positions.write_text(
        "Bracket Type=4,,,\n"
        "Sample Type,File Name,Sample ID,Position\n"
        "Unknown,Blank_Pos_01,GE8,GE8\n"
        "Unknown,QC_01,RH5,RH5\n"
        "Unknown,Sample_309,RE12,RE12\n"
        "Unknown,Blank_Pos_02,GE8,GE8\n")
    assert read_sequence(positions) == {"Blank_Pos_01": 1, "QC_01": 2,
                                        "Sample_309": 3, "Blank_Pos_02": 4}

    # Where the column really is an injection number, row order agrees with it — so preferring
    # row order costs nothing on the sequence the old behaviour was written against.
    injections = tmp_path / "injections.csv"
    injections.write_text(
        "Sample Type,File Name,Sample ID\n"
        "Unknown,Blank_01,Inj_01\n"
        "Unknown,QC_01,Inj_02\n"
        "Unknown,Sample_A,Inj_03\n")
    assert read_sequence(injections) == {"Blank_01": 1, "QC_01": 2, "Sample_A": 3}


def test_a_file_listed_twice_keeps_its_first_injection(tmp_path):
    """One sequence interleaving both polarities lists each name once per polarity."""
    from lipidloop.blanks import read_sequence

    path = tmp_path / "interleaved.csv"
    path.write_text("Sample Type,File Name,Sample ID\n"
                    "Unknown,Blank_01,A1\nUnknown,Blank_01,A1\nUnknown,QC_01,B2\n")
    assert read_sequence(path) == {"Blank_01": 1, "QC_01": 3}
