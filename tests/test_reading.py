"""Tests for the reader quirks that decide which peaks reach the dot product.

Each of these pins a rule taken from the LipiDex source rather than from the .mgf/.msp format,
because a tidier reading of either format changes the scores.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.mgf import read_mgf          # noqa: E402
from lipidloop.msp import parse_msp         # noqa: E402

MGF = """BEGIN IONS
TITLE=QC.1.1.1 File:"QC.raw", NativeID:"scan=1"
RTINSECONDS=60.5
PEPMASS=500.1234567 372688.2
CHARGE=1+
61.5000 5000.0
100.0000 10000.0
400.0000 8000.0
499.0000 9000.0
60.0000 7000.0
100.0500 1.0
END IONS
BEGIN IONS
TITLE=QC.2.2.1 File:"QC.raw", NativeID:"scan=2"
RTINSECONDS=120.0
PEPMASS=300.0
CHARGE=1-
100.0000 1000.0
END IONS
"""

MSP = """Name: Test-PC 16:0_18:1 [M+H]+;
MW: 500.1234
PRECURSORMZ: 500.1234
Comment: Name=Test OptimalPolarity=true Type=LipiDex
Num Peaks: 4
184.0733 999 "C5H15N1O4P1_Fragment_[]"
100.0000 4 "low_Fragment_[]"
498.5000 500 "near_precursor_Fragment_[]"
250.0000 500 "mid_Fragment_[]"

Name: Other-PE 18:0 [M-H]-;
MW: 400.0
PRECURSORMZ: 400.0
Comment: Name=Other OptimalPolarity=false
Num Peaks: 1
200.0000 999
"""


@pytest.fixture
def files(tmp_path):
    (tmp_path / "s.mgf").write_text(MGF)
    (tmp_path / "l.msp").write_text(MSP)
    return tmp_path


def test_mgf_peak_filters(files):
    spectra = list(read_mgf(files / "s.mgf"))
    kept = dict(zip(spectra[0].mz, spectra[0].intensity))
    assert 100.0 in kept          # ordinary peak
    assert 61.5 in kept           # just above the 61.0 low-mass cutoff
    assert 400.0 in kept
    assert 499.0 not in kept      # within 1.5 Da of the precursor
    assert 60.0 not in kept       # below the low-mass cutoff
    assert 100.05 not in kept     # raw intensity not above 1.0


def test_mgf_precursor_and_rt_rounded_to_three_places(files):
    first = list(read_mgf(files / "s.mgf"))[0]
    assert first.precursor == 500.123
    assert first.retention == 1.008          # 60.5 s / 60


def test_mgf_numbering_counts_every_spectrum(files):
    spectra = list(read_mgf(files / "s.mgf"))
    assert [s.number for s in spectra] == [0, 1]
    assert spectra[1].polarity == "-"


def test_mgf_scaling_drops_peaks_below_half_a_percent(files):
    """Base peak to 999; anything under 5 goes. 5000/10000*999 = 499.5 survives."""
    first = list(read_mgf(files / "s.mgf"))[0]
    first.scale_intensities()
    assert max(first.intensity) == pytest.approx(999.0)
    assert min(first.intensity) >= 5.0


def test_msp_drops_peaks_within_two_daltons_of_precursor(files):
    entry = list(parse_msp(files / "l.msp"))[0]
    assert 498.5 not in entry.mz
    assert 184.0733 in entry.mz


def test_msp_drops_peaks_below_five_after_scaling(files):
    """4/999 of the base peak scales to under 5 and is discarded."""
    entry = list(parse_msp(files / "l.msp"))[0]
    assert 100.0 not in entry.mz
    assert 250.0 in entry.mz


def test_msp_flags_and_polarity_from_the_name(files):
    lipidex, other = list(parse_msp(files / "l.msp"))
    assert lipidex.name == "Test-PC 16:0_18:1 [M+H]+;"
    assert lipidex.is_lipidex and lipidex.optimal_polarity
    assert lipidex.polarity == "positive"
    assert lipidex.annotations[lipidex.mz.index(184.0733)] == "C5H15N1O4P1_Fragment_[]"
    assert not other.is_lipidex and not other.optimal_polarity
    assert other.polarity == "negative"
