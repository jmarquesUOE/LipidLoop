"""The adduct column, the compound-group key, and `Associated_Spectra.csv`.

`Lipid.java :: toString` drops the adduct on purpose — the identity is the molecule, not the ion.
That is right for the name and leaves the table unable to answer two questions it is regularly
asked: which ion was quantified, and are these two rows at one retention time two adducts of one
molecule. `Na - H` and `+2 carbons, +3 double bonds` differ by 2.4 mDa, so the second question is
not hypothetical.
"""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.associated import COLUMNS, write_associated_spectra   # noqa: E402
from lipidloop.peakfinder import META_COLUMNS, PeakFinderResult, write_results   # noqa: E402
from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample   # noqa: E402
from lipidloop.search import RESULT_COLUMNS   # noqa: E402


def lipid(name, adduct, precursor, sample, ms2_id=1, rank=1, dot=900):
    out = Lipid(retention=5.0, precursor=precursor, sample=sample, dot=dot, rev_dot=980,
                lipid_string=f"{name} {adduct};", lib_precursor=precursor, purity=99,
                is_lipidex=True, purity_array=[], fragment_masses=[],
                ms2_id=ms2_id, rank=rank)
    return out


def group_with(lipids, quant_ion, purity=99.0, sum_id="PC 34:1"):
    group = CompoundGroup(name="g", mw=759.578, retention=5.0, max_area=1.0,
                          has_ms2=True, areas=[1.0])
    group.quant_ion = quant_ion
    candidates = []
    for one in lipids:
        match = next((c for c in candidates if c.lipid_name == one.lipid_name), None)
        if match:
            match.add(one)
        else:
            candidates.append(LipidCandidate(one))
    group.lipid_candidates = candidates
    group.purity = purity
    group.sum_id = sum_id
    group.final_lipid_id = candidates[0]
    return group


# ── the adduct itself ──────────────────────────────────────────────────────────────────────

def test_adduct_is_reported_for_the_identified_row():
    sample = Sample(file="a")
    group = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample)], quant_ion=760.585)
    assert group.identification()[0] == "PC 16:0_18:1"
    assert group.adduct() == "[M+H]+"


def test_the_quantified_ion_comes_first_when_a_molecule_is_seen_as_several():
    """Two ions of one molecule. The row is quantified on one of them, and that is the one a
    reader needs first — the others are supporting evidence, not what `Quant Ion` measures."""
    sample = Sample(file="a")
    group = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample, ms2_id=1),
                        lipid("PC 16:0_18:1", "[M+Na]+", 782.567, sample, ms2_id=2)],
                       quant_ion=782.567)
    assert group.adduct().split("; ")[0] == "[M+Na]+"
    assert set(group.adduct().split("; ")) == {"[M+H]+", "[M+Na]+"}


def test_unidentified_row_has_no_adduct():
    empty = CompoundGroup(name="g", mw=1.0, retention=1.0, max_area=0.0, has_ms2=False, areas=[])
    assert empty.adduct() == ""


# ── the columns ────────────────────────────────────────────────────────────────────────────

def test_adduct_and_group_id_are_written_before_the_sample_columns(tmp_path):
    sample = Sample(file="sample_one")
    group = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample)], quant_ion=760.585)
    group.areas = [1234.0]
    out = tmp_path / "r.csv"
    write_results(PeakFinderResult([group], [sample]), out, ms2_support=False, scores=False)
    header, body = list(csv.reader(out.open(newline="")))[:2]
    assert header.index("Adduct") < header.index("sample_one")
    assert header.index("Compound Group") < header.index("sample_one")
    assert body[header.index("Adduct")] == "[M+H]+"
    assert body[header.index("Compound Group")] == "1"


def test_group_id_is_the_same_in_the_filtered_and_unfiltered_tables(tmp_path):
    """The point of the id: a row dropped from `Final_Results.csv` keeps the number it has in
    `Unfiltered_Results.csv`, so the two can be read together and `Associated_Spectra` can point
    at either. Numbering the filtered table from 1 would silently renumber everything."""
    sample = Sample(file="s")
    kept = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample)], quant_ion=760.585)
    dropped = group_with([lipid("PE 16:0_18:1", "[M+H]+", 718.538, sample)], quant_ion=718.538,
                         sum_id="PE 34:1")
    kept.areas = dropped.areas = [1.0]
    dropped.keep, dropped.filter_reason = False, "Adduct of existing identified peak"
    result = PeakFinderResult([dropped, kept], [sample])

    write_results(result, tmp_path / "final.csv")
    write_results(result, tmp_path / "unfiltered.csv", unfiltered=True)
    final = list(csv.DictReader((tmp_path / "final.csv").open(newline="")))
    unfiltered = list(csv.DictReader((tmp_path / "unfiltered.csv").open(newline="")))

    assert len(final) == 1 and len(unfiltered) == 2
    assert final[0]["Compound Group"] == "2"          # not renumbered to 1
    assert {r["Compound Group"] for r in unfiltered} == {"1", "2"}


def test_every_optional_column_off_gives_exactly_the_lipidex_columns(tmp_path):
    sample = Sample(file="s")
    group = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample)], quant_ion=760.585)
    group.areas = [1.0]
    out = tmp_path / "r.csv"
    write_results(PeakFinderResult([group], [sample]), out, scores=False, ms2_support=False,
                  adduct=False, group_id=False, lipid_key=False, shorthand_column=False)
    header = next(csv.reader(out.open(newline="")))
    assert header == ["Retention Time (min)", "Quant Ion", "Polarity", "Area (max)",
                      "Identification", "Lipid Class", "Features Found", "s", ""]


def test_every_column_the_writer_can_emit_is_declared_metadata(tmp_path):
    """`META_COLUMNS` is how every reader finds where the injections start. If the writer can
    emit a column that is not in it, that column is read as a sample and the areas silently
    gain a phantom injection — which looks like data, not like a bug."""
    sample = Sample(file="an_injection")
    group = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample)], quant_ion=760.585)
    group.areas = [1.0]
    out = tmp_path / "r.csv"
    write_results(PeakFinderResult([group], [sample]), out, unfiltered=True, id_source=True,
                  ms2_support=True, scores=True, adduct=True, group_id=True)
    header = next(csv.reader(out.open(newline="")))
    unaccounted = [c for c in header if c and c not in META_COLUMNS and c != "an_injection"]
    assert unaccounted == []


# ── Associated_Spectra.csv ─────────────────────────────────────────────────────────────────

def _search_file(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        for row in rows:
            writer.writerow({c: row.get(c, "") for c in RESULT_COLUMNS})


def _row(ms2_id, name, precursor, rank=1, dot=900):
    return {"MS2 ID": ms2_id, "Retention Time (min)": 5.0, "Rank": rank, "Identification": name,
            "Precursor Mass": precursor, "Library Mass": precursor, "Delta m/z": -1.2,
            "Dot Product": dot, "Reverse Dot Product": 980, "Purity": 99,
            "Spectral Components": "", "Optimal Polarity": "true", "LipiDex Spectrum": "true",
            "Library": "LipiDex_HCD_Formic.msp", "Potential Fragments": ""}


def test_spectra_are_keyed_to_the_group_they_landed_in(tmp_path):
    sample = Sample(file="s")
    used = lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample, ms2_id=7)
    group = group_with([used], quant_ion=760.585)
    group.areas = [1.0]
    _search_file(tmp_path / "s_Results.csv",
                 [_row(7, "PC 16:0_18:1 [M+H]+;", 760.585),
                  _row(9, "PE 16:0_18:1 [M+H]+;", 718.538)])   # never attached to a group

    written = write_associated_spectra(PeakFinderResult([group], [sample]),
                                       {"s": tmp_path / "s_Results.csv"},
                                       tmp_path / "Associated_Spectra.csv")
    rows = list(csv.DictReader((tmp_path / "Associated_Spectra.csv").open(newline="")))
    assert written == 2 and len(rows) == 2
    by_id = {int(r["MS2 ID"]): r for r in rows}
    assert by_id[7]["Compound Group"] == "1"
    assert by_id[7]["Associated"] == "Associated"
    assert by_id[7]["Adduct"] == "[M+H]+"
    assert by_id[7]["Sample"] == "s"
    # The unattached spectrum is still reported — it was acquired and searched, and a table that
    # only showed the winners could not be used to ask why something is missing.
    assert by_id[9]["Compound Group"] == ""
    assert by_id[9]["Associated"] == ""


def test_all_ranks_are_carried_through(tmp_path):
    sample = Sample(file="s")
    group = group_with([lipid("PC 16:0_18:1", "[M+H]+", 760.585, sample, ms2_id=7)],
                       quant_ion=760.585)
    group.areas = [1.0]
    _search_file(tmp_path / "s_Results.csv",
                 [_row(7, "PC 16:0_18:1 [M+H]+;", 760.585, rank=1, dot=900),
                  _row(7, "PC 16:1_18:0 [M+H]+;", 760.585, rank=2, dot=880)])
    write_associated_spectra(PeakFinderResult([group], [sample]),
                             {"s": tmp_path / "s_Results.csv"},
                             tmp_path / "a.csv")
    rows = list(csv.DictReader((tmp_path / "a.csv").open(newline="")))
    assert [int(r["Rank"]) for r in rows] == [1, 2]
    # Only rank 1 was accepted by the peak finder, and the runner-up says so by having no group.
    assert rows[0]["Compound Group"] == "1" and rows[1]["Compound Group"] == ""


def test_header_is_stable(tmp_path):
    _search_file(tmp_path / "s_Results.csv", [])
    write_associated_spectra(PeakFinderResult([], []), {"s": tmp_path / "s_Results.csv"},
                             tmp_path / "a.csv")
    assert next(csv.reader((tmp_path / "a.csv").open(newline=""))) == COLUMNS
    assert "GaussianScore" not in COLUMNS   # no equivalent here; a constant zero would mislead
