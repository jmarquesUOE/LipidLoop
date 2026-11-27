"""Tests for the adduct database and the mass tests built on it.

Every filter in the pairwise sweep reduces to one question — could these two quant ions be the
same neutral molecule — and that question is only answerable with a formula per adduct.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.adducts import (Adduct, check_adduct_known, check_adduct_unknown,  # noqa: E402
                                check_dimer, check_fragment, formula_mass, load_adducts,
                                read_adducts)

PROTON = 1.007825035
SODIUM = 22.9897677

ADDUCTS = [Adduct("[M+H]+", "H1", False, "+", 1),
           Adduct("[M+Na]+", "Na1", False, "+", 1),
           Adduct("[M+NH4]+", "N1H4", False, "+", 1),
           Adduct("[M-H]-", "H-1", True, "-", 1)]


def test_formula_mass_handles_explicit_and_negative_counts():
    assert formula_mass("H1") == pytest.approx(PROTON)
    assert formula_mass("Na1") == pytest.approx(SODIUM)
    assert formula_mass("H-1") == pytest.approx(-PROTON)          # a proton lost
    assert formula_mass("O-1H-1") == pytest.approx(-(15.99491463 + PROTON))  # water less H
    assert formula_mass("H0") == 0.0
    assert formula_mass("") == 0.0


def test_formula_mass_sums_multiple_elements():
    assert formula_mass("N1H4") == pytest.approx(14.003074 + 4 * PROTON)
    assert formula_mass("C2H3O2") == pytest.approx(2 * 12.0 + 3 * PROTON + 2 * 15.99491463)


def test_same_molecule_as_two_adducts_is_recognised():
    """A 760 Da neutral seen as [M+H]+ and [M+Na]+."""
    neutral = 760.0
    assert check_adduct_unknown(ADDUCTS, neutral + PROTON, neutral + SODIUM)


def test_unrelated_masses_are_not_adducts():
    assert not check_adduct_unknown(ADDUCTS, 760.0 + PROTON, 903.4567)


def test_known_adduct_uses_the_identification_to_fix_the_neutral_mass():
    neutral = 760.0
    assert check_adduct_known(ADDUCTS, neutral + PROTON, neutral + SODIUM, "[M+H]+")
    # Wrong adduct named for mass1 -> the neutral mass is wrong -> no match
    assert not check_adduct_known(ADDUCTS, neutral + PROTON, neutral + SODIUM, "[M+Na]+")
    assert not check_adduct_known(ADDUCTS, neutral + PROTON, neutral + SODIUM, None)
    assert not check_adduct_known(ADDUCTS, neutral + PROTON, neutral + SODIUM, "[M+Xx]+")


def test_dimer_detected_in_either_direction():
    neutral = 400.0
    monomer = neutral + PROTON
    dimer = 2 * neutral + PROTON
    assert check_dimer(ADDUCTS, monomer, dimer, "+")
    assert check_dimer(ADDUCTS, dimer, monomer, "+")
    assert not check_dimer(ADDUCTS, monomer, dimer, "-")   # polarity must match


def test_in_source_fragment_matches_a_predicted_fragment():
    assert check_fragment([184.0733, 104.1075], [184.0740])     # within 20 ppm
    assert not check_fragment([184.0733], [200.0])
    assert not check_fragment([], [184.0733])
    assert not check_fragment(None, None)


def test_reads_the_real_adduct_database():
    folders = sorted((Path(__file__).resolve().parents[1] / "data/lipidex_src/libraries").glob("*"))
    if not folders:
        pytest.skip("LipiDex library folders not present")
    adducts = load_adducts(folders)
    names = {a.name for a in adducts}
    assert {"[M+H]+", "[M-H]-", "[M+Na]+", "[M+NH4]+"} <= names
    by_name = {a.name: a for a in adducts}
    assert by_name["[M+H]+"].mass == pytest.approx(PROTON)
    assert by_name["[M-H]-"].mass == pytest.approx(-PROTON)
    assert by_name["[M-H]-"].polarity == "-"
