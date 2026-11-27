"""`Method_Lists/`: the safe, import-ready exclusion and inclusion lists written once per study."""
import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import method_lists as ml  # noqa: E402


def _study(tmp_path: Path, pol: str, msp: Path) -> Path:
    analysis = tmp_path / "Analysis"
    res = analysis / "Results" / pol
    res.mkdir(parents=True)
    (analysis / f"config_{pol}.json").write_text(json.dumps({"libraries": [str(msp)]}))
    sign = "+" if pol == "Pos" else "-"
    # three candidates: one collides with an identified lipid, one with a library-only lipid,
    # one clean
    (res / "Exclusion_List.csv").write_text(
        "m/z,Charge,Polarity,MS2 scans wasted,RT observed (min),In blank,Series,Note\n"
        f"760.5851,1,{sign},900,0.1-15.0,yes,,\"background\"\n"
        f"500.3000,1,{sign},400,0.1-15.0,yes,PEG,\"PEG contamination\"\n"
        f"445.1200,1,{sign},300,0.1-15.0,yes,polysiloxane,\"polysiloxane contamination\"\n")
    (res / "Final_Results.csv").write_text(
        "Compound Group,Retention Time (min),Quant Ion,Polarity,Identification\n"
        f"1,6.1,760.5852,{sign},PC 34:1 [M+H]+\n"
        f"2,6.5,720.5000,{sign},DECOY_ PC 1:0_2:0\n")
    (res / "Inclusion_List.csv").write_text("m/z,Retention (min)\n630.4785,6.23\n594.3779,0.20\n")
    return analysis


def _msp(tmp_path: Path) -> Path:
    p = tmp_path / "lib.msp"
    p.write_text("Name: PC 34:1 [M+H]+\nPRECURSORMZ: 760.5851\nNum Peaks: 1\n184.0733 999\n\n"
                 "Name: PE 20:0_4:0 [M+H]+\nPRECURSORMZ: 500.3010\nNum Peaks: 1\n141.0191 999\n\n"
                 "Name: PE 20:0_4:0 [M-H]-\nPRECURSORMZ: 500.3010\nNum Peaks: 1\n140.0118 999\n\n")
    return p


def test_colliding_entries_are_dropped_and_reported(tmp_path):
    analysis = _study(tmp_path, "Pos", _msp(tmp_path))
    summary = ml.write_method_lists(analysis, say=lambda *_: None)
    out = analysis / "Method_Lists"
    rows = list(csv.DictReader((out / "Exclusion_List_Pos_Thermo.csv").open()))
    assert [r["Mass [m/z]"] for r in rows] == ["445.1200"]
    assert rows[0]["Polarity"] == "Positive" and rows[0]["Start [min]"] == ""
    assert summary["polarities"]["Pos"] == {"exclusion_candidates": 3, "exclusion_kept": 1,
                                            "exclusion_dropped": 2, "inclusion_entries": 2}
    report = (out / "Exclusion_Safety_Pos.md").read_text()
    assert "760.5851" in report and "PC 34:1" in report
    assert "500.3000" in report and "PE 20:0_4:0" in report
    assert (out / "README.md").exists() and (out / "summary.json").exists()


def test_inclusion_list_gets_thermo_layout_and_window(tmp_path):
    analysis = _study(tmp_path, "Pos", _msp(tmp_path))
    ml.write_method_lists(analysis, inclusion_window=0.5, say=lambda *_: None)
    rows = list(csv.DictReader((analysis / "Method_Lists" / "Inclusion_List_Pos_Thermo.csv").open()))
    assert rows[0]["Mass [m/z]"] == "630.4785"
    assert (rows[0]["Start [min]"], rows[0]["End [min]"]) == ("5.73", "6.73")
    assert rows[1]["Start [min]"] == "0.00"          # never negative
    assert rows[0]["Polarity"] == "Positive"


def test_negative_mode_protects_free_fatty_acids(tmp_path):
    analysis = _study(tmp_path, "Neg", _msp(tmp_path))
    res = analysis / "Results" / "Neg"
    (res / "Exclusion_List.csv").write_text(
        "m/z,Charge,Polarity,MS2 scans wasted,RT observed (min),In blank,Series,Note\n"
        "255.2330,1,-,2679,0.1-15.0,yes,,\"background\"\n"          # palmitate [M-H]-
        "445.1200,1,-,300,0.1-15.0,yes,polysiloxane,\"polysiloxane contamination\"\n")
    (res / "Final_Results.csv").write_text(
        "Compound Group,Retention Time (min),Quant Ion,Polarity,Identification\n")
    ml.write_method_lists(analysis, say=lambda *_: None)
    rows = list(csv.DictReader((analysis / "Method_Lists" / "Exclusion_List_Neg_Thermo.csv").open()))
    assert [r["Mass [m/z]"] for r in rows] == ["445.1200"]
    assert "protected class" in (analysis / "Method_Lists" / "Exclusion_Safety_Neg.md").read_text()


def test_decoy_rows_and_decoy_libraries_do_not_protect_anything(tmp_path):
    msp = _msp(tmp_path)
    decoy = tmp_path / "DECOY1_lib.msp"
    decoy.write_text("Name: DECOY_ PC 9:0_9:0 [M+H]+\nPRECURSORMZ: 445.1200\nNum Peaks: 1\n184.0733 999\n\n")
    analysis = _study(tmp_path, "Pos", msp)
    (analysis / "config_Pos.json").write_text(json.dumps({"libraries": [str(msp), str(decoy)]}))
    ml.write_method_lists(analysis, say=lambda *_: None)
    rows = list(csv.DictReader((analysis / "Method_Lists" / "Exclusion_List_Pos_Thermo.csv").open()))
    assert [r["Mass [m/z]"] for r in rows] == ["445.1200"]


def test_missing_polarity_is_skipped(tmp_path):
    analysis = _study(tmp_path, "Pos", _msp(tmp_path))
    summary = ml.write_method_lists(analysis, say=lambda *_: None)
    assert list(summary["polarities"]) == ["Pos"]
    assert not (analysis / "Method_Lists" / "Exclusion_List_Neg_Thermo.csv").exists()
