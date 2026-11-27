"""Pos/Neg injections paired from the data, not the file names (v1.1 item 20)."""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import pairing, polarity                 # noqa: E402
from lipidloop.pairing import Injection                 # noqa: E402
import qc_report                                         # noqa: E402

strip = polarity._strip_polarity


def _mzml(path, terms, stamp="2026-11-12T17:01:55Z"):
    body = "".join(f'<cvParam accession="{t}"/>' for t in terms)
    path.write_text(f'<mzML><run id="x" startTimeStamp="{stamp}">{body}</run></mzML>')
    return path


def test_scan_polarity_reads_the_scans_and_the_start_time(tmp_path):
    assert pairing.scan_polarity(_mzml(tmp_path / "a.mzML", [pairing.POSITIVE] * 3)) == \
        ("Pos", "2026-11-12T17:01:55Z")
    assert pairing.scan_polarity(_mzml(tmp_path / "b.mzML", [pairing.NEGATIVE]))[0] == "Neg"
    # a polarity-switching file is reported as such, not forced into one mode
    assert pairing.scan_polarity(_mzml(tmp_path / "c.mzML", [pairing.POSITIVE, pairing.NEGATIVE]))[0] == "mixed"
    assert pairing.scan_polarity(tmp_path / "missing.mzML") == ("unknown", "")


def test_sequence_positions_skip_the_bracket_line(tmp_path):
    seq = tmp_path / "sequence.csv"
    seq.write_text("Bracket Type=4,\nSample Type,File Name,Sample ID,Position\n"
                   "Unknown,QC_01_P,QC,RF10\nUnknown,QC_01_N,QC,RF10\nUnknown,EV_1_P,EV,RA1\n")
    assert pairing.sequence_positions(seq) == {"QC_01_P": "RF10", "QC_01_N": "RF10", "EV_1_P": "RA1"}


def _inj(pol, column, vial="", started="", scans=None):
    return Injection(pol, column, vial, "raw header" if vial else "", started, scans or pol)


def test_repeats_of_one_vial_pair_in_acquisition_order():
    """Three pooled-QC injections from one vial: the k-th positive goes with the k-th negative,
    whatever the names say."""
    pos = [_inj("Pos", "QC_01_P", "RF10", "T1"), _inj("Pos", "QC_02_P", "RF10", "T3"),
           _inj("Pos", "EV_1_P", "RA1", "T5")]
    neg = [_inj("Neg", "EV_1_N", "RA1", "T6"), _inj("Neg", "QC_02_N", "RF10", "T4"),
           _inj("Neg", "QC_01_N", "RF10", "T2")]
    pairs = pairing.pair(pos, neg, strip)
    assert [(x.pos.column, x.neg.column, x.column) for x in pairs] == [
        ("QC_01_P", "QC_01_N", "QC_01"), ("QC_02_P", "QC_02_N", "QC_02"), ("EV_1_P", "EV_1_N", "EV_1")]
    assert {x.matched_by for x in pairs} == {"vial (raw header)"}
    assert pairing.disagreements(pairs, strip) == []


def test_vial_wins_over_names_and_the_disagreement_is_reported():
    pos = [_inj("Pos", "S1_Pos", "RA1"), _inj("Pos", "S2_Pos", "RA2")]
    neg = [_inj("Neg", "S2_Neg", "RA1"), _inj("Neg", "S1_Neg", "RA2")]
    pairs = pairing.pair(pos, neg, strip)
    assert [(x.pos.column, x.neg.column) for x in pairs] == [("S1_Pos", "S2_Neg"), ("S2_Pos", "S1_Neg")]
    # the names do not agree, so the column keeps the positive file's name
    assert [x.column for x in pairs] == ["S1_Pos", "S2_Pos"]
    assert len(pairing.disagreements(pairs, strip)) == 2


def test_unequal_injections_per_vial_are_refused():
    pos = [_inj("Pos", "QC_01_P", "RF10"), _inj("Pos", "QC_02_P", "RF10"), _inj("Pos", "A_P", "RA1")]
    neg = [_inj("Neg", "QC_01_N", "RF10"), _inj("Neg", "A_N", "RA1"), _inj("Neg", "B_N", "RA1")]
    with pytest.raises(ValueError, match="did not measure the same vials"):
        pairing.pair(pos, neg, strip)


def test_without_vials_the_names_pair_and_a_mismatch_is_refused():
    pos = [_inj("Pos", "S1_POS"), _inj("Pos", "S2_POS")]
    pairs = pairing.pair(pos, [_inj("Neg", "S2_NEG"), _inj("Neg", "S1_NEG")], strip)
    assert [(x.column, x.neg.column, x.matched_by) for x in pairs] == [
        ("S1", "S1_NEG", "file name"), ("S2", "S2_NEG", "file name")]
    with pytest.raises(ValueError, match="did not measure the same vials"):
        pairing.pair(pos, [_inj("Neg", "S1_NEG"), _inj("Neg", "S9_NEG")], strip)


def test_one_vial_for_every_file_does_not_pair_and_names_are_used():
    """A header that records the same vial for every file says nothing about which is which."""
    pos = [_inj("Pos", "S1_Pos", "1"), _inj("Pos", "S2_Pos", "1"), _inj("Pos", "S3_Pos", "1")]
    neg = [_inj("Neg", "S3_Neg", "1"), _inj("Neg", "S1_Neg", "1"), _inj("Neg", "S2_Neg", "1")]
    assert {x.matched_by for x in pairing.pair(pos, neg, strip)} == {"file name"}
    assert [(x.pos.column, x.neg.column) for x in pairing.pair(pos, neg, strip)] == [
        ("S1_Pos", "S1_Neg"), ("S2_Pos", "S2_Neg"), ("S3_Pos", "S3_Neg")]


def test_a_file_in_the_wrong_folder_is_reported_not_moved():
    injections = [_inj("Pos", "A_P", scans="Neg"), _inj("Neg", "A_N"), _inj("Pos", "B_P", scans="mixed"),
                  _inj("Neg", "B_N", scans="unknown")]
    assert [i.column for i in pairing.wrong_folder(injections)] == ["A_P", "B_P"]


def _tables(tmp_path, pos_cols, neg_cols):
    root = tmp_path / "Results"
    for pol, cols in (("Pos", pos_cols), ("Neg", neg_cols)):
        (root / pol).mkdir(parents=True)
        values = ",".join(str(k + 1) for k in range(len(cols)))
        (root / pol / polarity.NORMALISED).write_text(
            f"Identification,{','.join(cols)}\n{pol} lipid,{values}\n")
    return root


def test_combine_uses_the_pairing_and_reorders_the_negative_columns(tmp_path):
    root = _tables(tmp_path, ["A1_P", "B2_P"], ["B2_N", "A1_N"])
    pos = [_inj("Pos", "A1_P", "RA1"), _inj("Pos", "B2_P", "RA2")]
    neg = [_inj("Neg", "B2_N", "RA2"), _inj("Neg", "A1_N", "RA1")]
    assert polarity.combine(tmp_path, pairs=pairing.pair(pos, neg, strip)) == 2
    lines = (root / polarity.COMBINED).read_text().splitlines()
    assert lines[0].split(",")[-2:] == ["A1", "B2"]
    # Neg values were 1 for B2_N and 2 for A1_N; under columns A1, B2 they read 2, 1
    assert lines[2].split(",")[-2:] == ["2", "1"]


def test_combine_refuses_a_column_the_pairing_does_not_cover(tmp_path):
    _tables(tmp_path, ["A1_P", "B2_P"], ["A1_N", "B2_N"])
    pairs = pairing.pair([_inj("Pos", "A1_P", "RA1")], [_inj("Neg", "A1_N", "RA1")], strip)
    with pytest.raises(ValueError, match="not in the pairing table"):
        polarity.combine(tmp_path, pairs=pairs)


def test_report_shows_the_pairing_and_flags_scans_of_the_other_polarity(tmp_path):
    pairs = pairing.pair([_inj("Pos", "A1_P", "RA1", scans="Neg"), _inj("Pos", "B2_P", "RA2")],
                         [_inj("Neg", "A1_N", "RA1"), _inj("Neg", "B2_N", "RA2")], strip)
    (tmp_path / "Results").mkdir()
    pairing.write(tmp_path / "Results" / pairing.TABLE, pairs)
    html = qc_report.pairing_block(tmp_path)
    assert "2 positive-mode injections" in html and "vial (raw header)" in html
    assert "Scans of another polarity" in html and "<code>A1_P</code> (Neg)" in html
    assert qc_report.pairing_block(tmp_path / "nowhere") == ""
