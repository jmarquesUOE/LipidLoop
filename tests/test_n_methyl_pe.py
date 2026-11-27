"""PE named as a symmetric N-methyl PE, and PC [M-CH3]- protected by a class model built on itself."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.coelution import ion_family   # noqa: E402
from lipidloop.mgf import SampleSpectrum   # noqa: E402
from lipidloop.msp import LibrarySpectrum   # noqa: E402
from lipidloop.search import Identification, prefer_pe_over_n_methyl   # noqa: E402

NMe2 = "PE-NMe2 16:1_16:1 [M-H]-;"
PE = "PE 16:1_18:1 [M-H]-;"


def hit(name, dot):
    lib = LibrarySpectrum(name=name, precursor_mz=714.5074, mz=[], intensity=[])
    return Identification(None, lib, 0.0, dot, 1000.0)


def order(mz, polarity="-", hits=None):
    ms2 = SampleSpectrum(number=1, precursor=714.5074, retention=7.8, polarity=polarity, mz=mz, intensity=[999.0] * len(mz))
    hits = hits or [hit(NMe2, 900.0), hit(PE, 800.0)]
    return [h.library_spectrum.name for h in prefer_pe_over_n_methyl(hits, ms2, 0.01)]


def test_pe_head_and_no_own_head_moves_the_pe_to_the_top():
    assert order([140.0118, 253.2173, 281.2486]) == [PE, NMe2]
    assert order([196.0380, 253.2173]) == [PE, NMe2]


def test_own_head_keeps_the_score_order():
    assert order([168.0431, 253.2173]) == [NMe2, PE]
    assert order([140.0118, 224.0693, 253.2173]) == [NMe2, PE]


def test_no_head_ion_at_all_keeps_the_score_order():
    assert order([253.2173, 281.2486]) == [NMe2, PE]


def test_positive_mode_and_other_classes_are_untouched():
    assert order([140.0118], polarity="+") == [NMe2, PE]
    other = [hit("PG 16:0_18:1 [M-H]-;", 900.0), hit(PE, 800.0)]
    assert order([140.0118], hits=other)[0].startswith("PG")


def test_a_decoy_is_not_promoted():
    decoys = [hit(NMe2, 900.0), hit("DECOY_ PE 15:1_19:1 [M-H]-;", 800.0)]
    assert order([140.0118], hits=decoys)[0] == NMe2


def _row(group, rt, mz, name, cls, adduct, z=""):
    return {"Compound Group": str(group), "Retention Time (min)": f"{rt}", "Quant Ion": f"{mz}", "Area (max)": "1e7",
            "Identification": name, "Lipid Class": cls, "Adduct": adduct, "Shorthand (LSI)": "", "Identification Source": "MS2",
            "RT Model Z": z, "Theoretical m/z": f"{mz}", "Mass Error (ppm)": "0.0", "S1": "1e7", "S2": "1e7"}


def test_a_class_model_does_not_exempt_pe_nme2_from_the_pc_anchor():
    neutral = 759.5778                                     # PC 34:1
    mz = round(neutral - 15.023475, 4)
    rows = [_row(1, 9.1, mz, "PE-NMe2 16:0_18:1", "PE-NMe2", "[M-H]-", z="0.09")]
    labels = ion_family(rows, "Neg", [], cross=([("PC 34:1", neutral, 9.1)], 0.0))
    assert "[M-CH3]-" in labels["1"][1]


def test_a_class_model_still_exempts_other_classes():
    neutral = 759.5778
    mz = round(neutral - 15.023475, 4)
    rows = [_row(1, 9.1, mz, "PE 35:1", "PE", "[M-H]-", z="0.09")]
    assert not ion_family(rows, "Neg", [], cross=([("PC 34:1", neutral, 9.1)], 0.0))
