"""The co-elution filter (v1.1 item 1): a delivered lipid that is really another lipid's ion.

The exact equivalence with the review's frozen scripts (ion family v3 + parent check v1b) is tested on real data
(outside this repository: QE two-by-two and three ground-truth sets, removed rows identical). These tests pin
the rules that decide, and the table handling, without mzML.
"""
import csv
import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import coelution   # noqa: E402
from lipidloop.coelution import CoelutionFilter   # noqa: E402

COLS = ["Compound Group", "Retention Time (min)", "Quant Ion", "Area (max)", "Identification", "Lipid Class",
        "Adduct", "Shorthand (LSI)", "Identification Source", "RT Model Z", "Theoretical m/z", "Mass Error (ppm)",
        "S1", "S2", "S3", "S4"]


def row(group, rt, mz, area, name="", cls="", adduct="", theo="", source="MS2", z=""):
    return {"Compound Group": str(group), "Retention Time (min)": f"{rt}", "Quant Ion": f"{mz}",
            "Area (max)": f"{area}", "Identification": name, "Lipid Class": cls, "Adduct": adduct,
            "Shorthand (LSI)": "", "Identification Source": source if name else "", "RT Model Z": z,
            "Theoretical m/z": f"{theo}" if theo else "", "Mass Error (ppm)": "0.0" if name else "",
            "S1": f"{area}", "S2": f"{area * 0.9}", "S3": f"{area * 1.1}", "S4": f"{area}"}


ISO = coelution.ISO


def test_an_m2_isotope_named_as_the_one_fewer_double_bond_species_is_labelled():
    """The M+2 of PC 36:5 [M+H]+ sits 8.9 mDa from PC 36:4 [M+H]+; named PC 36:4, at the right isotope ratio."""
    a = row(1, 7.30, 780.5538, 1e8, "PC 36:5", "PC", "[M+H]+", 780.5538)
    t = row(2, 7.30, 780.5538 + 2 * ISO, 1.0e7, "PC 36:4", "PC", "[M+H]+", 782.5694)
    labels = coelution.ion_family([a, t], "Pos", ms2=[])
    assert set(labels) == {"2"}
    assert labels["2"][0] == "PC 36:5 [M+H]+" and labels["2"][1].startswith("M+2 isotope")


def test_a_real_species_above_the_isotope_ratio_is_kept():
    a = row(1, 7.30, 780.5538, 1e8, "PC 36:5", "PC", "[M+H]+", 780.5538)
    t = row(2, 7.30, 780.5538 + 2 * ISO, 6.0e7, "PC 36:4", "PC", "[M+H]+", 782.5694)   # 0.6 x parent
    assert coelution.ion_family([a, t], "Pos", ms2=[]) == {}


def test_a_named_sodium_relation_needs_ms2_evidence_of_sodium():
    """Na - H = C2H4 - 3H2 within 2.4 mDa: without the target's own MS2 showing sodium, mass alone never labels."""
    a = row(1, 7.30, 780.5538, 1e8, "PC 36:5", "PC", "[M+H]+", 780.5538)
    neutral = 780.5538 - 1.007276
    t = row(2, 7.30, round(neutral + 22.989218, 4), 1e7, "PC 38:8", "PC", "[M+H]+", 802.5381)
    assert coelution.ion_family([a, t], "Pos", ms2=[]) == {}


def test_a_decoy_is_a_target_never_an_anchor():
    a = row(1, 7.30, 780.5538, 1e8, "DECOY_ PC 36:5", "DECOY_", "[M+H]+")
    t = row(2, 7.30, 780.5538 + 2 * ISO, 1e7, "PC 36:4", "PC", "[M+H]+", 782.5694)
    assert coelution.ion_family([a, t], "Pos", ms2=[]) == {}


def test_the_v1b_apex_rule():
    """PC 34:3 [M+Na]+ on the QE 2x2: r 0.959, apex 0.051 min early, ratio 0.07 -- missed at a fixed 0.05."""
    v1b, v1 = CoelutionFilter(), CoelutionFilter(apex_tol_high_r=0.05)
    assert coelution._explained(0.959, -0.051, 0.07, 0.5, v1b)
    assert not coelution._explained(0.959, -0.051, 0.07, 0.5, v1)
    assert not coelution._explained(0.93, -0.051, 0.07, 0.5, v1b)     # r below 0.95: 0.05 still applies
    assert coelution._explained(0.93, -0.04, 0.07, 0.5, v1b)
    assert not coelution._explained(0.97, -0.08, 0.07, 0.5, v1b)
    assert not coelution._explained(math.nan, 0.0, 0.07, 0.5, v1b)
    assert not coelution._explained(0.99, 0.0, 0.6, 0.5, v1b)          # too strong to be a derivative


def _analysis(tmp_path, rows):
    res = tmp_path / "Results" / "Pos"
    res.mkdir(parents=True)
    for name in (coelution.FINAL, coelution.FILTERED, coelution.NORMALISED):
        with (res / name).open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=COLS)
            w.writeheader()
            w.writerows(rows)
    with (res / coelution.UNFILTERED).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLS + ["Filter Status"])
        w.writeheader()
        w.writerows([dict(r, **{"Filter Status": ""}) for r in rows])
    return res


def test_apply_removes_from_the_delivered_tables_and_labels_the_rest(tmp_path):
    a = row(1, 7.30, 780.5538, 1e8, "PC 36:5", "PC", "[M+H]+", 780.5538)
    t = row(2, 7.30, 780.5538 + 2 * ISO, 1.0e7, "PC 36:4", "PC", "[M+H]+", 782.5694)
    o = row(3, 9.10, 876.8015, 5e7, "TG 52:2", "TG", "[M+NH4]+", 876.8015)
    res = _analysis(tmp_path, [a, t, o])

    summary = coelution.apply(tmp_path, {"Pos": None}, CoelutionFilter())

    read = lambda n: list(csv.DictReader((res / n).open()))
    assert [r["Compound Group"] for r in read(coelution.FILTERED)] == ["1", "3"]
    assert [r["Compound Group"] for r in read(coelution.NORMALISED)] == ["1", "3"]
    final = {r["Compound Group"]: r for r in read(coelution.FINAL)}
    assert len(final) == 3, "nothing is deleted from Final_Results"
    assert final["2"]["Ion Of"] == "PC 36:5 [M+H]+" and final["2"]["Ion Test"] == "ion family"
    assert final["1"]["Ion Of"] == "" and final["3"]["Ion Of"] == ""
    status = {r["Compound Group"]: r["Filter Status"] for r in read(coelution.UNFILTERED)}
    assert status["2"].startswith("co-elution: M+2 isotope") and status["1"] == ""
    assert summary["removed"]["Pos"]["rows"] == 1
    assert json.loads((tmp_path / "coelution_filter.json").read_text())["settings"]["apex_tol_high_r"] == 0.075


def test_disabled_changes_nothing(tmp_path):
    a = row(1, 7.30, 780.5538, 1e8, "PC 36:5", "PC", "[M+H]+", 780.5538)
    t = row(2, 7.30, 780.5538 + 2 * ISO, 1.0e7, "PC 36:4", "PC", "[M+H]+", 782.5694)
    res = _analysis(tmp_path, [a, t])
    before = (res / coelution.FILTERED).read_text()
    assert coelution.apply(tmp_path, {"Pos": None}, CoelutionFilter(enabled=False)) == {}
    assert (res / coelution.FILTERED).read_text() == before


def test_label_columns_never_land_among_the_injections(tmp_path):
    """Readers take everything after the last metadata column as an injection (qc.load_run): appended at
    the end, `Ion Of` was read as a sample and the quality report crashed."""
    from lipidloop.peakfinder import META_COLUMNS
    a = row(1, 7.30, 780.5538, 1e8, "PC 36:5", "PC", "[M+H]+", 780.5538)
    t = row(2, 7.30, 780.5538 + 2 * ISO, 1.0e7, "PC 36:4", "PC", "[M+H]+", 782.5694)
    res = _analysis(tmp_path, [a, t])
    coelution.apply(tmp_path, {"Pos": None}, CoelutionFilter())
    coelution.apply(tmp_path, {"Pos": None}, CoelutionFilter())          # twice: columns not repeated
    columns = next(csv.reader((res / coelution.FINAL).open()))
    start = max(i for i, c in enumerate(columns) if c in META_COLUMNS) + 1
    assert columns[start:] == ["S1", "S2", "S3", "S4"]
    assert all(columns.count(c) == 1 for c in coelution.LABEL_COLUMNS)
