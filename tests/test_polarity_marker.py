"""The polarity marker must not outlive the tables it describes."""
import os, time, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[0] / "src"))


def _stale(analysis: Path) -> bool:
    """Mirror of the guard in run_study.py, exercised without running a study."""
    marker = analysis / "polarity_filter.json"
    if not marker.exists():
        return True
    tables = list(analysis.glob("Results/*/Final_Results_Filtered.csv"))
    newest = max((t.stat().st_mtime for t in tables), default=0.0)
    return newest > marker.stat().st_mtime


def test_marker_written_after_the_tables_is_honoured(tmp_path):
    (tmp_path / "Results/Pos").mkdir(parents=True)
    (tmp_path / "Results/Pos/Final_Results_Filtered.csv").write_text("a\n")
    time.sleep(0.01)
    (tmp_path / "polarity_filter.json").write_text("{}")
    assert not _stale(tmp_path), "filter already applied to THESE tables — must not run twice"


def test_a_rebuilt_table_invalidates_the_marker(tmp_path):
    """The bug: a pipeline re-run rewrites the tables unfiltered and leaves the record behind.

    The study then shipped both polarities' shared classes — 853 rows where 689 were delivered.
    """
    (tmp_path / "Results/Pos").mkdir(parents=True)
    (tmp_path / "polarity_filter.json").write_text("{}")
    table = tmp_path / "Results/Pos/Final_Results_Filtered.csv"
    table.write_text("a\n")
    marker = tmp_path / "polarity_filter.json"
    os.utime(table, (marker.stat().st_mtime + 10, marker.stat().st_mtime + 10))
    assert _stale(tmp_path), "tables rebuilt after the marker — the filter must run again"


def test_polarity_agreement_takes_the_median_across_duplicate_rows(tmp_path, monkeypatch):
    """One name, several peaks — not whichever row is last in the file.

    `changes[name] = ...` silently overwrote, so a duplicated lipid was compared using an arbitrary
    chromatographic peak in each polarity, which need not be the same species. On a real study
    `SM d41:2` carried three positive rows at +0.073, +0.141 and -1.752; the last won, and the
    -1.75 was compared against negative mode and reported to the client as an OPPOSITE direction.
    """
    import csv, sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop import polarity

    samples = {"A1": "ctrl", "A2": "ctrl", "B1": "test", "B2": "test"}
    for pol in ("Pos", "Neg"):
        d = tmp_path / "Results" / pol
        d.mkdir(parents=True)
        meta = tmp_path / f"metadata_{pol}.csv"
        with meta.open("w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["sample", "group"])
            for k, v in samples.items(): w.writerow([k, v])
        # two rows agree that the lipid does not move; a third, LAST, says it collapses
        bodies = [["1", "1", "1", "1"], ["1", "1", "1", "1"], ["1", "1", "0.25", "0.25"]]
        if pol == "Neg":
            bodies = [["1", "1", "1.1", "1.1"]]
        with (d / polarity.NORMALISED).open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Identification", "Retention Time (min)", *samples])
            for i, b in enumerate(bodies):
                w.writerow(["SM d41:2", 10.0 + i, *b])
            # the comparison needs three shared lipids before it will report anything
            w.writerow(["PC 34:1", 9.0, "1", "1", "1.2", "1.2"])
            w.writerow(["PE 36:2", 9.5, "1", "1", "0.9", "0.9"])

    out = polarity.agreement(tmp_path, {p: str(tmp_path / f"metadata_{p}.csv")
                                        for p in ("Pos", "Neg")}, "group")
    assert out, "should have produced a comparison"
    assert out["by_lipid"]["SM d41:2"] == "agree", (
        "the median of the three positive rows is ~0, matching negative; last-row-wins gave -2")


def test_sequence_covering_both_polarities_gives_each_injection_its_own_number(tmp_path):
    """One File Name, two injections — the polarity is only in `Path`.

    A real study's sequence carried 152 distinct names across 304 rows: the same name for a
    sample's positive and negative injection. Keeping the first occurrence gave the negative
    injection the positive one's number, and since the acquired files are `<name>_Pos` the bare
    key matched no result column at all — so run order came back empty and every section needing
    it reported "not testable" rather than "the sequence did not match".
    """
    import csv, sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop.blanks import read_sequence

    seq = tmp_path / "seq.csv"
    with seq.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Bracket Type=4"])
        w.writerow(["Sample Type", "File Name", "Path"])
        w.writerow(["Unknown", "S1", r"C:\data\Pos"])
        w.writerow(["Unknown", "S1", r"C:\data\Neg"])
        w.writerow(["Unknown", "S2", r"C:\data\Pos"])

    order = read_sequence(seq)
    assert order["S1_Pos"] == 1, "positive injection is first"
    assert order["S1_Neg"] == 2, "negative is a DIFFERENT injection, not the same number"
    assert order["S1"] == 1, "the bare key keeps its first occurrence, for single-polarity runs"
    assert order["S2_Pos"] == 3


def test_mzml_input_needs_no_conversion_and_no_thermo_parser(tmp_path):
    """A study already in mzML must run without ThermoRawFileParser being installed.

    The pipeline sent every input to the Thermo converter with no format check, and located the
    converter before asking whether anything needed converting — so a study of perfectly good mzML
    failed on a missing Thermo reader it never needed. That is most public data, any collaborator
    sending mzML, and every already-converted dataset in the validation catalogue.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop.convert import Converter, needs_conversion

    assert needs_conversion("x.raw") is True
    assert needs_conversion("x.mzML") is False
    assert needs_conversion("x.mzml") is False, "writers disagree on case; both are converted files"

    mz = tmp_path / "a.mzML"
    mz.write_text("<mzML/>")
    # An executable that does not exist: if convert_all tried to run it, this would raise.
    broken = Converter(tmp_path / "no-such-parser")
    assert broken.convert_all([mz], tmp_path / "out") == [mz], "mzML must pass through untouched"


def test_polarity_suffix_matching_is_case_insensitive():
    """`_POS` and `_pos` name the same vial as `_Pos`.

    Public deposits write `NIST_Iterative_20ev_1_POS`; the facility writes `_Pos`. A
    case-sensitive strip left those columns untouched, so the two polarities harmonised to
    disjoint column sets and `combine` raised "the two runs did not measure the same vials" on
    data where they plainly had.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop.polarity import _strip_polarity

    assert _strip_polarity("NIST_1_POS", "Pos") == "NIST_1"
    assert _strip_polarity("NIST_1_Pos", "Pos") == "NIST_1"
    assert _strip_polarity("NIST_1_pos", "Pos") == "NIST_1"
    assert _strip_polarity("Pool_NEG_01", "Neg") == "Pool_01", "the token is not always a suffix"
    assert _strip_polarity("Positive_control", "Pos") == "Positive_control", "bounded by _ or ends"


def test_a_failed_cross_polarity_join_does_not_discard_a_completed_analysis(tmp_path):
    """Combining the polarities is a convenience, not a precondition.

    `polarity.combine` requires the two runs to have measured the same vials. That is not always
    true: one public study names its polarities `20210303P_...` and `20210303N_...` — a letter
    inside a date token rather than a suffix — and carries 16 positive injections against 12
    negative. Raising there killed a run whose per-polarity results were already written and
    valid, and exited non-zero, which reads as "the analysis failed" when both polarities in fact
    completed.
    """
    import re
    from pathlib import Path
    source = (Path(__file__).resolve().parents[1] / "scripts" / "run_study.py").read_text()
    call = source[source.index("polarity_filter.combine("):]
    guarded = call[:call.index("except polarity_filter.UndecidedClass")]
    assert "except ValueError" in guarded, \
        "combine must be guarded — a join failure cannot discard two completed polarities"
    assert "combined table not written" in guarded, "and it must say what was skipped"


def test_an_empty_delivered_table_is_reported_not_crashed(tmp_path):
    """A study can legitimately filter down to nothing, and the report must say so.

    `load_run` read its columns from `rows[0]`, so a header-only table raised IndexError deep
    inside the report builder — which reads as a crash rather than as "nothing survived the
    filters". Two public datasets hit this in a single batch: a two-injection run cannot satisfy a
    70% presence filter, so every feature is removed and the delivered table is a header alone.
    """
    import csv, sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop import qc

    empty = tmp_path / "Final_Results_Filtered.csv"
    with empty.open("w", newline="") as fh:
        csv.writer(fh).writerow(["Retention Time (min)", "Quant Ion", "Polarity", "Area (max)",
                                 "Identification", "Lipid Class", "Features Found",
                                 "S1_Pos", "S2_Pos"])
    run = qc.load_run(empty)
    assert run.columns == ["S1_Pos", "S2_Pos"], "injections come from the header"
    assert run.rows == []
    assert run.areas.shape == (0, 2), "an empty table still needs a correctly shaped array"


def test_a_sequence_beside_the_raw_files_is_found(tmp_path):
    """A sequence often sits in the polarity folder, not the study root.

    Searching the root alone returned nothing for a study carrying `Pos/Seq_Rand.csv`, so
    injection order went missing and every drift and confounding section of the report would have
    said "not testable" — the same silent failure as the polarity-suffix bug, from a different
    cause.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
    import importlib
    rs = importlib.import_module("run_study")

    study = tmp_path / "study"
    (study / "Pos").mkdir(parents=True)
    (study / "Neg").mkdir()
    seq = study / "Pos" / "Seq_Rand.csv"
    seq.write_text("Sample Type,File Name\nUnknown,S1\n")
    assert rs.find_sequence(study) == str(seq)

    # the root still wins when both exist — it describes the whole study
    root = study / "sequence.csv"
    root.write_text("Sample Type,File Name\nUnknown,S1\n")
    assert rs.find_sequence(study) == str(root)


def test_injection_order_comes_from_the_acquisition_timestamp(tmp_path):
    """The instrument's record of when a run started beats a planned sequence.

    `startTimeStamp` is written from the instrument's own record, so it is what happened rather
    than what was scheduled — and it exists for every file, including a polarity whose sequence
    was never saved. On the skin study, `r1`, `r2` and `r3` of one sample were acquired at 07:07,
    21:55 the PREVIOUS day, and 08:39, so filename order is not run order.

    The file mtime is no substitute: every raw file in that study carries an mtime within ten
    seconds of every other, because they were copied off the instrument in one go.
    """
    import csv, sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop import qc

    mzml = tmp_path / "mzml"
    mzml.mkdir()
    times = {"B": "2024-02-12T20:53:56Z", "A": "2024-02-13T07:07:44Z", "C": "2024-02-13T08:39:50Z"}
    for stem, when in times.items():
        (mzml / f"{stem}.mzML").write_text(f'<mzML><run startTimeStamp="{when}"></run></mzML>')

    table = tmp_path / "Final_Results.csv"
    with table.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Retention Time (min)", "Quant Ion", "Polarity", "Area (max)",
                    "Identification", "Lipid Class", "Features Found", "A", "B", "C"])
        w.writerow([1.0, 100.0, "+", 5.0, "PC 34:1", "PC", 1, 1, 1, 1])

    run = qc.load_run(table, mzml_files=sorted(mzml.glob("*.mzML")))
    assert run.order == {"B": 1, "A": 2, "C": 3}, \
        f"order follows acquisition time, not filename: {run.order}"


def test_the_report_says_where_the_injection_order_came_from():
    """Every drift number rests on the order being right, so its source belongs in the report.

    The two sources answer different questions — the instrument records what was acquired, the
    sequence file what was scheduled — and a disagreement between them means a run was repeated,
    reordered or dropped. Preferring one silently hides that.
    """
    from pathlib import Path
    source = (Path(__file__).resolve().parents[1] / "scripts" / "qc_report.py").read_text()
    assert "Injection order taken from" in source, "the report must name its order source"
    assert "order_disagreement" in source, "and must surface a clash between the two sources"


def test_format_is_detected_from_the_path_not_the_suffix(tmp_path):
    """Waters `.raw` is a DIRECTORY and Thermo `.raw` is a FILE — same suffix, two vendors.

    Dispatching on the suffix alone hands a Waters directory to ThermoRawFileParser, which fails
    in a way that reads as a corrupt file rather than as the wrong reader.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop.convert import detect_format

    thermo = tmp_path / "sample.raw"
    thermo.write_bytes(b"\x01\xa1F\x00i\x00n\x00")
    assert detect_format(thermo) == "thermo"

    waters = tmp_path / "run.raw"
    waters.mkdir()
    (waters / "_FUNC001.DAT").write_bytes(b"\x00")
    assert detect_format(waters) == "waters", "a .raw DIRECTORY is Waters, not Thermo"

    agilent = tmp_path / "run.d"
    (agilent / "AcqData").mkdir(parents=True)
    assert detect_format(agilent) == "agilent"

    bruker = tmp_path / "tims.d"
    bruker.mkdir()
    (bruker / "analysis.tdf").write_bytes(b"\x00")
    assert detect_format(bruker) == "bruker", "Bruker ships Linux readers — handled natively"

    assert detect_format(tmp_path / "x.wiff") == "sciex"
    assert detect_format(tmp_path / "x.mzML") == "mzml"


def test_decoy_classes_do_not_enter_the_polarity_decision():
    """A decoy is not a lipid and has no polarity to assign.

    Decoy entries carry a `DECOY_` prefix, so their class parses as the pseudo-class `DECOY_`. It
    appears in both polarities and is in no class-polarity table, so `check` refused the entire
    run — correct behaviour on an unknown class, wrong input. It killed one dataset in an
    overnight batch.
    """
    import sys
    from pathlib import Path
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop import polarity

    table = {"PC": None}
    v = polarity.validate({"PC": 3, "DECOY_": 9}, {"PC": 2, "DECOY_": 7}, table)
    assert v.undecided == [], f"a decoy must not force a polarity decision: {v.undecided}"

    # a genuinely unknown real class still refuses — the guard must not be weakened
    v2 = polarity.validate({"PC": 3, "Cer[ADS]": 4}, {"PC": 2, "Cer[ADS]": 5}, table)
    assert v2.undecided == ["Cer[ADS]"], "an unknown real class must still stop the run"
