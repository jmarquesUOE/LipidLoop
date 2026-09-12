"""Assigning a class to the polarity that measures it better.

The costly mistakes here are all deletions that nothing downstream would reveal, so the tests are
mostly about what must NOT be dropped.
"""
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import polarity      # noqa: E402


DEFAULT = Path(__file__).resolve().parents[1] / "data/class_polarity.csv"


def test_the_facility_table_loads_and_covers_the_shared_classes():
    table = polarity.load(DEFAULT)
    assert len(table) > 30
    # the ten that appeared in both polarities on the study this was built for
    for name in ("PC", "SM", "LysoPC", "Plasmenyl-PC", "PE", "Plasmenyl-PE",
                 "PI", "PS", "PE-NMe2", "LysoPE"):
        assert name in table, name
    assert table["PC"].polarity == "Pos" and table["PE"].polarity == "Neg"
    # ceramide rows must be marked as resting on the library, not on chemistry
    assert table["Cer[ADS]"].basis == "library"


def test_a_study_file_overrides_the_facility_default(tmp_path):
    override = tmp_path / "class_polarity.csv"
    override.write_text("class,polarity,basis,note\nPC,Neg,study,because this study says so\n")
    table = polarity.load(DEFAULT, override)
    assert table["PC"].polarity == "Neg" and table["PC"].source == "study"
    assert table["PE"].polarity == "Neg" and table["PE"].source == "default"


def test_single_polarity_classes_are_never_eligible():
    """The expensive failure. `Cer[ADS]` is negative-only because the positive ceramide library
    covers phyto bases, not because of chemistry — a rule reaching beyond the shared classes
    deletes it, and it is central to a separate project."""
    table = polarity.load(DEFAULT)
    got = polarity.validate({"PC": 100, "TG": 200}, {"PE": 50, "Cer[ADS]": 18}, table)
    assert got.shared == []
    assert set(got.single) == {"PC", "TG", "PE", "Cer[ADS]"}
    assert got.ok


def test_a_class_in_both_polarities_with_no_assignment_stops_the_run():
    got = polarity.validate({"PC": 10, "Newthing": 4}, {"PC": 2, "Newthing": 6}, {})
    assert got.undecided == ["Newthing", "PC"]
    with pytest.raises(polarity.UndecidedClass) as raised:
        polarity.check(got)
    assert "Newthing" in str(raised.value)


def test_stale_and_moot_entries_are_reported_but_do_not_stop_anything():
    table = polarity.load(DEFAULT)
    got = polarity.validate({"PC": 10}, {"PE": 10}, table)
    assert "CL" in got.stale                 # configured, seen in neither
    assert "PC" in got.moot and "PE" in got.moot   # configured, now one polarity only
    assert got.ok


def write_table(path, rows, columns):
    """A results table. Any extra key a row carries becomes a column, so a test can add a
    retention time without every other test having to know about it."""
    extra = [k for r in rows for k in r
             if k not in ("Identification", "Lipid Class") and k not in columns]
    seen, ordered = set(), []
    for k in extra:
        if k not in seen:
            seen.add(k)
            ordered.append(k)
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["Identification", "Lipid Class"] + ordered + columns)
        w.writeheader()
        w.writerows(rows)


def test_apply_drops_only_the_losing_polarity_and_recentres(tmp_path):
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    for pol, rows in (
        ("Pos", [{"Identification": "PC 34:1", "Lipid Class": "PC", "S1": 100, "S2": 200},
                 {"Identification": "PE 36:2", "Lipid Class": "PE", "S1": 100, "S2": 200},
                 {"Identification": "TG 52:2", "Lipid Class": "TG", "S1": 100, "S2": 200}]),
        ("Neg", [{"Identification": "PC 34:1", "Lipid Class": "PC", "S1": 50, "S2": 100},
                 {"Identification": "PE 36:2", "Lipid Class": "PE", "S1": 50, "S2": 100},
                 {"Identification": "FA 18:1", "Lipid Class": "FA", "S1": 50, "S2": 100}])):
        (root / pol).mkdir(parents=True)
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    result = polarity.apply(tmp_path, polarity.load(DEFAULT))
    assert result["validation"].ok

    def classes(pol, name):
        with (root / pol / name).open(newline="") as fh:
            return [r["Lipid Class"] for r in csv.DictReader(fh)]

    # PC belongs to positive, PE to negative; TG and FA are single-polarity and untouched
    assert classes("Pos", polarity.FILTERED) == ["PC", "TG"]
    assert classes("Neg", polarity.FILTERED) == ["PE", "FA"]

    # and the normalised table is centred on what survived
    with (root / "Pos" / polarity.NORMALISED).open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    import statistics
    med = [statistics.median([float(r[c]) for r in rows]) for c in columns]
    assert max(med) / min(med) - 1 < 0.01


def test_agreement_is_measured_before_the_filter_rewrites_anything(tmp_path):
    """It can only be measured while both tables still contain the shared classes.

    Computed after the rewrite it returns nothing at all, silently — which is what it did on the
    first attempt, next to a docstring saying it had to be computed first.
    """
    root = tmp_path / "Results"
    columns = ["A1", "A2", "B1", "B2"]
    for pol, scale in (("Pos", 1.0), ("Neg", 1.0)):
        (root / pol).mkdir(parents=True)
        rows = []
        for n in range(6):
            # each lipid moves by a different amount, and the two polarities agree on the ordering.
            # A fixture where every fold change is identical has no rank variation, so Spearman is
            # undefined — correct behaviour of the statistic, and a useless test.
            for name, cls, step in ((f"PC 3{n}:1", "PC", n + 1), (f"PE 3{n}:2", "PE", n + 2)):
                rows.append({"Identification": name, "Lipid Class": cls,
                             "A1": 100, "A2": 100,
                             "B1": 100 * (1 + step * scale), "B2": 100 * (1 + step * scale)})
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)
    meta = {}
    for pol in ("Pos", "Neg"):
        path = tmp_path / f"metadata_{pol}.csv"
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["sample", "group"])
            w.writeheader()
            for c in columns:
                w.writerow({"sample": c, "group": "A" if c.startswith("A") else "B"})
        meta[pol] = str(path)

    out = polarity.apply(tmp_path, polarity.load(DEFAULT), metadata=meta)
    got = out["summary"]["agreement"]
    assert got, "agreement must be computed while both tables still hold the shared classes"
    assert got["shared"] == 12 and got["rho"] > 0.9
    # and the filter did run, so the tables no longer share those classes
    assert polarity.identified_classes(root / "Pos" / polarity.FILTERED).keys() == {"PC"}
    assert polarity.identified_classes(root / "Neg" / polarity.FILTERED).keys() == {"PE"}


def test_the_polarity_token_is_removed_wherever_it_sits():
    """It is not always a suffix. This facility's samples end `_Pos` but its pools are
    `Pool_Pos_01`, so a suffix-only rule matched the samples, left the pools mismatched, and the
    two blocks had disjoint columns for four injections."""
    assert polarity._strip_polarity("x_TS_F_Pos", "Pos") == "x_TS_F"
    assert polarity._strip_polarity("Pool_Pos_01", "Pos") == "Pool_01"
    assert polarity._strip_polarity("Pool_Neg_01", "Neg") == "Pool_01"
    # a name that merely starts with the letters is untouched
    assert polarity._strip_polarity("Positive_control", "Pos") == "Positive_control"


def test_combining_refuses_when_the_two_polarities_measured_different_vials(tmp_path):
    """Better to stop than to emit a matrix whose columns are half empty on each side."""
    root = tmp_path / "Results"
    for pol, columns in (("Pos", ["S1_Pos", "S2_Pos"]), ("Neg", ["S1_Neg", "S9_Neg"])):
        (root / pol).mkdir(parents=True)
        rows = [{"Identification": "PC 34:1", "Lipid Class": "PC",
                 columns[0]: 10, columns[1]: 20}]
        write_table(root / pol / polarity.NORMALISED, rows, columns)
    with pytest.raises(ValueError, match="did not measure the same vials"):
        polarity.combine(tmp_path)


def test_applying_the_filter_twice_does_not_erase_what_it_did(tmp_path):
    """The second pass finds nothing to decide, because the classes it removed are no longer in
    both polarities. That is harmless for the tables and destructive for the record: the summary
    the report reads gets overwritten with zeros, and a rebuilt report says the filter did nothing.
    """
    import json
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    for pol in ("Pos", "Neg"):
        (root / pol).mkdir(parents=True)
        rows = [{"Identification": "PC 34:1", "Lipid Class": "PC", "S1": 100, "S2": 100},
                {"Identification": "PE 36:2", "Lipid Class": "PE", "S1": 100, "S2": 100}]
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    table = polarity.load(DEFAULT)
    first = polarity.apply(tmp_path, table)["summary"]
    assert first["shared_classes"] == {"PC": "Pos", "PE": "Neg"}

    second = polarity.apply(tmp_path, table)["summary"]
    assert second["shared_classes"] == {}, "second pass legitimately finds nothing"
    # which is exactly why run_study must not call it again on a report-only rebuild
    assert json.loads((tmp_path / "polarity_filter.json").read_text())["shared_classes"] == {}


def test_chain_names_are_carried_from_the_polarity_that_resolved_them(tmp_path):
    """Quantify on the better ion, name from the better spectrum.

    A choline lipid in positive mode puts almost everything into its head-group fragment and gives
    no usable acyl ions — 0 of 165 PC rows resolved chains in positive on the study this was built
    for — while negative resolves them and is the polarity about to be dropped.
    """
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    (root / "Pos").mkdir(parents=True)
    (root / "Neg").mkdir(parents=True)
    rt = polarity.RT_COLUMN
    pos = [{"Identification": "PC 34:1", "Lipid Class": "PC", rt: 8.00, "S1": 100, "S2": 100},
           {"Identification": "PC 36:2", "Lipid Class": "PC", rt: 8.50, "S1": 100, "S2": 100}]
    neg = [{"Identification": "PC 16:0_18:1", "Lipid Class": "PC", rt: 8.02, "S1": 5, "S2": 5},
           # two resolved names for one sum composition: a real isomer pair, so decline
           {"Identification": "PC 18:0_18:2", "Lipid Class": "PC", rt: 8.51, "S1": 5, "S2": 5},
           {"Identification": "PC 16:0_20:2", "Lipid Class": "PC", rt: 8.52, "S1": 5, "S2": 5}]
    for pol, rows in (("Pos", pos), ("Neg", neg)):
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    out = polarity.apply(tmp_path, polarity.load(DEFAULT))
    assert out["summary"]["chains_carried"] == 1

    with (root / "Pos" / polarity.FILTERED).open(newline="") as fh:
        rows = {r["Identification"]: r for r in csv.DictReader(fh)}
    # ⚠ Keyed on the SUM composition, because the row keeps the name this polarity measured. The
    # molecular name is evidence from another injection and lives in its own column.
    assert "PC 34:1" in rows, "the row keeps the name its own polarity measured"
    assert rows["PC 34:1"]["Chains Resolved"] == "PC 16:0_18:1", "chains recorded beside it"
    assert rows["PC 34:1"]["Chains From"] == "Neg", "and says where they came from"
    assert rows["PC 34:1"]["Chains RT Delta"] == "+0.020", "with how far apart the two peaks were"
    assert rows["PC 34:1"]["S1"] == "100", "quantification stays on the positive ion"
    assert "PC 36:2" in rows, "two isomers share that sum composition — declined, left as sum"
    assert not rows["PC 36:2"]["Chains Resolved"], "and nothing was recorded for it"


def test_chains_are_declined_when_the_two_polarities_disagree_on_retention(tmp_path):
    """⚠ RT confirms the pair. Composition alone would accept a peak eluting a minute away, which
    is a different compound wearing the same sum formula."""
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    rt = polarity.RT_COLUMN
    (root / "Pos").mkdir(parents=True)
    (root / "Neg").mkdir(parents=True)
    pos = [{"Identification": "PC 34:1", "Lipid Class": "PC", rt: 8.00, "S1": 100, "S2": 100}]
    neg = [{"Identification": "PC 16:0_18:1", "Lipid Class": "PC", rt: 11.40, "S1": 5, "S2": 5}]
    for pol, rows in (("Pos", pos), ("Neg", neg)):
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    out = polarity.apply(tmp_path, polarity.load(DEFAULT))
    assert out["summary"]["chains_carried"] == 0

    with (root / "Pos" / polarity.FILTERED).open(newline="") as fh:
        rows = {r["Identification"]: r for r in csv.DictReader(fh)}
    assert not rows["PC 34:1"]["Chains Resolved"], "3.4 min apart is not the same peak"


def test_rt_confirms_but_never_selects_between_candidates(tmp_path):
    """★ The GM3 case, in miniature. In reversed phase +2C and +1DB nearly cancel, so a species and
    its (n-2C, -1DB) partner co-elute BY CHEMISTRY. Here `PC 36:2` sits 0.01 min from the positive
    `PC 34:1` row — closer than `PC 34:1`'s own partner — and must still not be chosen for it,
    because it is a different composition. Composition decides; RT only vetoes."""
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    rt = polarity.RT_COLUMN
    (root / "Pos").mkdir(parents=True)
    (root / "Neg").mkdir(parents=True)
    pos = [{"Identification": "PC 34:1", "Lipid Class": "PC", rt: 8.00, "S1": 100, "S2": 100}]
    neg = [{"Identification": "PC 18:0_18:2", "Lipid Class": "PC", rt: 8.01, "S1": 5, "S2": 5},
           {"Identification": "PC 16:0_18:1", "Lipid Class": "PC", rt: 8.06, "S1": 5, "S2": 5}]
    for pol, rows in (("Pos", pos), ("Neg", neg)):
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    polarity.apply(tmp_path, polarity.load(DEFAULT))
    with (root / "Pos" / polarity.FILTERED).open(newline="") as fh:
        rows = {r["Identification"]: r for r in csv.DictReader(fh)}
    got = rows["PC 34:1"]["Chains Resolved"]
    assert got == "PC 16:0_18:1", f"the 34:1 composition, not the nearer 36:2 peak — got {got!r}"


def test_the_combined_matrix_has_no_repeated_column_name(tmp_path):
    """The result tables already carry a `Polarity` column holding the quant ion's sign. Adding a
    second one of the same name is not a naming quibble: pandas silently renames one, which of the
    two you get depends on the library, and the first check written against it read +/- when it
    meant Pos/Neg."""
    root = tmp_path / "Results"
    columns = ["S1_Pos", "S1_Neg"]
    for pol, column in (("Pos", "S1_Pos"), ("Neg", "S1_Neg")):
        (root / pol).mkdir(parents=True)
        path = root / pol / polarity.NORMALISED
        with path.open("w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["Identification", "Lipid Class", "Polarity", column])
            w.writeheader()
            w.writerow({"Identification": "TG 52:2", "Lipid Class": "TG",
                        "Polarity": "+" if pol == "Pos" else "-", column: 100})
    polarity.combine(tmp_path)
    with (root / polarity.COMBINED).open(newline="") as fh:
        header = next(csv.reader(fh))
    assert len(header) == len(set(header)), f"repeated column: {header}"
    assert "Mode" in header and "Polarity" in header


def test_both_spellings_of_polarity_are_read(tmp_path):
    """`Pos`/`Neg` and `positive`/`negative` both appear in the shipped table, because later rows
    were written by hand. Reading only one vocabulary dropped five rules — BA, Cer[EOS],
    HexCer[EOS], ASM, AHexCer — and every one has a live library, so the first study identifying one
    in both polarities would have refused with rc=3 for no stated reason."""
    path = tmp_path / "t.csv"
    path.write_text("class,polarity,basis,note\n"
                    "PC,Pos,chemistry,\n"
                    "BA,negative,chemistry,\n"
                    "ASM,positive,library,\n")
    table = polarity.load(path)
    assert len(table) == 3
    assert table["BA"].polarity == "Neg"
    assert table["ASM"].polarity == "Pos"


def test_an_unreadable_polarity_raises_rather_than_dropping_the_rule(tmp_path):
    """⚠ The failure direction matters. This module refuses to guess a polarity because a wrong
    guess deletes real lipids from one polarity and nothing downstream shows it. Silently skipping a
    malformed row disables that protection just as effectively, and quietly."""
    path = tmp_path / "t.csv"
    path.write_text("class,polarity,basis,note\nPC,sideways,chemistry,\n")
    with pytest.raises(ValueError, match="sideways"):
        polarity.load(path)


def test_every_row_of_the_shipped_table_is_read():
    """A regression guard on the real file: the count in the log (`polarity filter: N classes`)
    must equal the number of data rows, or rules are being lost."""
    import csv as _csv
    from pathlib import Path as _Path
    src = _Path(__file__).resolve().parents[1] / "data" / "class_polarity.csv"
    with src.open() as fh:
        rows = list(_csv.DictReader(l for l in fh if not l.lstrip().startswith("#")))
    data_rows = [r for r in rows if (r.get("class") or "").strip()]
    assert len(polarity.load(src)) == len(data_rows)


def test_chains_are_taken_from_the_same_polarity_when_it_resolved_them(tmp_path):
    """★ A sum composition and its molecular form can sit on two rows of the SAME table — the
    search resolved chains for one spectrum of the peak and not another. That is one molecule
    named twice, not two findings, and the chains belong beside the sum row.

    Measured across five studies: 70 such pairs, 17 of them co-eluting."""
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    rt = polarity.RT_COLUMN
    (root / "Pos").mkdir(parents=True)
    (root / "Neg").mkdir(parents=True)
    pos = [{"Identification": "PC 34:1", "Lipid Class": "PC", rt: 8.00, "S1": 100, "S2": 100},
           # same peak, chains resolved, in THIS polarity's own table
           {"Identification": "PC 16:0_18:1", "Lipid Class": "PC", rt: 8.03, "S1": 90, "S2": 90}]
    neg = [{"Identification": "PE 36:2", "Lipid Class": "PE", rt: 9.00, "S1": 5, "S2": 5}]
    for pol, rows in (("Pos", pos), ("Neg", neg)):
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    polarity.apply(tmp_path, polarity.load(DEFAULT))
    with (root / "Pos" / polarity.FILTERED).open(newline="") as fh:
        rows = {r["Identification"]: r for r in csv.DictReader(fh)}
    assert rows["PC 34:1"]["Chains Resolved"] == "PC 16:0_18:1"
    assert rows["PC 34:1"]["Chains From"] == "Pos", "from its own table, and says so"
    assert "PC 16:0_18:1" in rows, "the molecular row is annotated, never deleted"


def test_a_chromatographic_isomer_in_the_same_polarity_is_left_alone(tmp_path):
    """⚠ 53 of those 70 pairs sit at DIFFERENT retention times. Those are real isomers — two peaks,
    two entries — and merging them would assert that one compound eluted twice."""
    root = tmp_path / "Results"
    columns = ["S1", "S2"]
    rt = polarity.RT_COLUMN
    (root / "Pos").mkdir(parents=True)
    (root / "Neg").mkdir(parents=True)
    pos = [{"Identification": "PC 34:1", "Lipid Class": "PC", rt: 8.00, "S1": 100, "S2": 100},
           {"Identification": "PC 16:0_18:1", "Lipid Class": "PC", rt: 10.90, "S1": 90, "S2": 90}]
    neg = [{"Identification": "PE 36:2", "Lipid Class": "PE", rt: 9.00, "S1": 5, "S2": 5}]
    for pol, rows in (("Pos", pos), ("Neg", neg)):
        for name in (polarity.FILTERED, polarity.NORMALISED):
            write_table(root / pol / name, [dict(r) for r in rows], columns)

    polarity.apply(tmp_path, polarity.load(DEFAULT))
    with (root / "Pos" / polarity.FILTERED).open(newline="") as fh:
        rows = {r["Identification"]: r for r in csv.DictReader(fh)}
    assert not rows["PC 34:1"]["Chains Resolved"], "2.9 min apart is a different peak"
