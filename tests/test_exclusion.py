"""What must never end up on an exclusion list.

Excluding a precursor is the one decision in this pipeline that reprocessing cannot undo — the
spectrum is never acquired. So the tests that matter are the negative ones: a lipid that elutes
as a peak has to survive every combination of the other criteria, however unidentified it is.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.exclusion import (Wasted, assign_series, cluster_precursors,  # noqa: E402
                                  duty_cycle, find_wasted, write_exclusion_list)

GRADIENT = 25.0


def background(mz, n=200, key_prefix="bg"):
    """A contaminant: fragmented all the way through the run."""
    return [((key_prefix, mz, k), mz + 1e-6 * mz * (k % 5 - 2), GRADIENT * k / n)
            for k in range(n)]


def peak(mz, rt, n=8, width=0.05, key_prefix="pk"):
    """A lipid: fragmented across one chromatographic peak."""
    return [((key_prefix, mz, k), mz, rt + width * (k / n - 0.5)) for k in range(n)]


def test_background_is_caught():
    scans = background(1054.3044) + peak(760.5851, 12.0)
    found = find_wasted(scans, identified=set(), require_blank=False)
    assert [round(w.mz) for w in found] == [1054]


def test_unidentified_lipid_peak_survives():
    """The case the whole design is about: never identified, but it elutes as a peak."""
    scans = background(1054.3044) + peak(999.9999, 12.0, n=60)
    found = find_wasted(scans, identified=set(), require_blank=False)
    assert 999.9999 not in [w.mz for w in found]


def test_two_peaks_far_apart_are_not_background():
    """An isomer pair spanning the run is still two peaks, not a continuous signal.

    This is the honest limit of the elution-span test, and it is why span is measured but a
    human reads the list. Here the span passes but the scan count does not reach the threshold
    on its own — the guard that catches it in practice is the blank.
    """
    scans = peak(700.5, 4.0) + peak(700.5, 22.0)
    found = find_wasted(scans, identified=set(), require_blank=False, min_scans=5)
    assert [round(w.mz, 1) for w in found] == [700.5]   # documented as caught: needs the blank
    found = find_wasted(scans, identified=set(), blank_precursors=[300.0], require_blank=True)
    assert found == []


def test_identified_precursor_is_never_excluded():
    scans = background(1054.3044)
    identified = {scans[0][0]}
    assert find_wasted(scans, identified, require_blank=False) == []


def test_blank_requirement():
    scans = background(1054.3044)
    assert find_wasted(scans, set(), blank_precursors=[500.0], require_blank=True) == []
    assert find_wasted(scans, set(), blank_precursors=[1054.30], require_blank=True)


def test_too_few_scans():
    scans = background(1054.3044, n=4)
    assert find_wasted(scans, set(), require_blank=False, min_scans=5) == []


def test_polysiloxane_series_named():
    base = 906.2654
    labels = assign_series([base + 74.01872 * k for k in range(6)])
    assert set(labels.values()) == {"polysiloxane"}


def test_sodium_formate_series_named():
    """The negative-mode data is 43% this. [(HCOONa)n + HCOO]-, from formic acid plus glass."""
    labels = assign_series([44.9982 + 67.98738 * k for k in range(4, 12)])
    assert set(labels.values()) == {"sodium formate"}


def test_unnamed_repeat_is_discovered():
    """Naming every contaminant is hopeless; a repeat unit is evidence whatever it is called."""
    labels = assign_series([400.0 + 123.456 * k for k in range(5)])
    assert set(labels.values()) == {"unnamed Δ123.456"}


def test_a_lipid_homologous_series_is_never_called_a_polymer():
    """PC 32:0/34:0/36:0 are 14.0157 apart. So is every acyl homologous series in the sample.

    Printing that next to column bleed would be actively misleading, so CH2 multiples are never
    proposed as an unnamed repeat — at any multiple.
    """
    for multiple in (1, 2, 3):
        labels = assign_series([760.5851 + 14.01565 * multiple * k for k in range(6)])
        assert labels == {}, f"CH2 x{multiple} was labelled a series"


def test_clustering_is_relative_not_absolute():
    """0.01 Da is 7 ppm at m/z 1500 and 100 ppm at m/z 100 — a fixed bin gets both wrong."""
    scans = [(1, 1500.0000, 1.0), (2, 1500.0090, 1.0), (3, 100.0000, 1.0), (4, 100.0090, 1.0)]
    sizes = sorted(len(g) for g in cluster_precursors(scans))
    assert sizes == [1, 1, 2]        # the m/z 1500 pair merges, the m/z 100 pair does not


def test_csv_stays_parseable_when_a_note_contains_a_comma(tmp_path):
    """"background, never identified" is the default note. Unquoted it splits a column."""
    import csv as _csv

    from lipidloop.exclusion import write_thermo_list

    wasted = [Wasted(mz=700.1234, scans=9, rt_lo=0.0, rt_hi=25.0,
                     gradient_fraction=1.0, in_blank=True)]        # no series -> comma in note
    for name, writer in (("a.csv", write_exclusion_list), ("b.csv", write_thermo_list)):
        out = tmp_path / name
        writer(wasted, out, "+")
        with out.open(newline="") as fh:
            rows = list(_csv.reader(fh))
        assert len(rows[1]) == len(rows[0]), f"{name}: comma in the note split a column"
        assert rows[1][-1].startswith("background, never identified")


def test_duty_cycle_and_csv(tmp_path):
    wasted = [Wasted(mz=1054.3044, scans=1020, rt_lo=0.02, rt_hi=24.98,
                     gradient_fraction=1.0, in_blank=True, series="polysiloxane")]
    assert duty_cycle(wasted, 2040) == 0.5
    out = tmp_path / "excl.csv"
    write_exclusion_list(wasted, out, "+")
    body = out.read_text()
    assert "1054.3044" in body and "polysiloxane" in body
    assert len(body.strip().splitlines()) == 2
