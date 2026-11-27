"""Spiked standards: the mass arithmetic, and the artefact screen that would delete them.

Both tests here exist because of a real error each. The supplier sheet listed two deuterated
fatty acids at their *unlabelled* masses; and the artefact screen, which is on by default and
needs no configuration, rejects exactly the mass defects that heavy deuteration produces.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.artefacts import plausible_mass, screen           # noqa: E402
from lipidloop.standards import Standard, StandardMix, formula_mass   # noqa: E402


class Spectrum:
    def __init__(self, mz, intensity):
        self.mz, self.intensity = mz, intensity


def test_deuterium_counts_towards_the_mass():
    """C16D31H1O2 is d31-palmitate at 287.43, not palmitic acid at 256.24.

    The sheet this was built from gave 256.2402 — the unlabelled mass — beside the labelled
    formula. Anything reading the mass column instead of the formula looks 31 Da too low.
    """
    assert formula_mass("C16H32O2") == pytest.approx(256.2402, abs=1e-3)
    assert formula_mass("C16D31H1O2") == pytest.approx(287.4348, abs=1e-3)
    assert formula_mass("C18D35H1O2") == pytest.approx(319.4912, abs=1e-3)
    # and the d7 lipids, where the sheet was right
    assert formula_mass("C41H73D7NO8P") == pytest.approx(752.6061, abs=1e-3)


def test_adduct_masses():
    pc = Standard("15:0-18:1-d7-PC", "C41H73D7NO8P", "[M+H]+", "+")
    assert pc.mz == pytest.approx(753.6134, abs=1e-3)
    fa = Standard("Palmitate_d31", "C16D31H1O2", "[M-H]-", "-")
    assert fa.mz == pytest.approx(286.4275, abs=1e-3)


def test_heavily_deuterated_standards_look_impossible_to_the_artefact_screen():
    """The premise of the exemption. d31 and d35 acids carry defects no CHNOPS ion can, which is
    precisely the test the screen uses to find real artefacts."""
    assert not plausible_mass(286.4275)     # d31-palmitate [M-H]-
    assert not plausible_mass(318.4839)     # d35-stearate [M-H]-
    # the d7 species are nowhere near the limit — the failure is selective, and so easy to miss
    assert plausible_mass(753.6134)
    assert plausible_mass(529.3994)


def test_declaring_a_standard_exempts_it_without_disarming_the_screen():
    standards = [Spectrum([286.4275, 500.1], [950.0, 10.0]) for _ in range(20)]
    assert [round(a.mz, 3) for a in screen(standards)] == [286.428], "would be deleted"
    assert screen(standards, protected=[286.4275]) == [], "declared, so kept"

    # a genuine artefact is still caught while the standard is protected
    artefact = [Spectrum([178.295, 500.1], [950.0, 10.0]) for _ in range(20)]
    assert [round(a.mz, 3) for a in screen(artefact, protected=[286.4275])] == [178.295]


def test_mix_loads_from_csv_and_computes_every_mass(tmp_path):
    path = tmp_path / "mix.csv"
    path.write_text("name,formula,adduct,polarity\n"
                    "PC d7,C41H73D7NO8P,[M+H]+,+\n"
                    "Palmitate_d31,C16D31H1O2,[M-H]-,-\n")
    mix = StandardMix.from_csv(path)
    assert len(mix.standards) == 2
    assert len(mix.for_polarity("+")) == 1
    assert mix.protected_mz("-") == pytest.approx([286.4275], abs=1e-3)
    # both polarities when none is named — what the pipeline protects
    assert len(mix.protected_mz()) == 2


def test_offset_check_flags_a_disagreement_and_refuses_when_incomparable():
    """The check exists to catch the case where the identification-derived offset is confidently
    wrong. Its default tolerance is a placeholder: the cross-day difference actually observed on
    this instrument was 2.96 ppm, which sits just inside 3.0, so a same-sequence check wants a
    tighter number than the default and the facility should set it deliberately."""
    from lipidloop.standards import check_offset
    assert check_offset(-4.0, -4.6)["agree"]
    assert not check_offset(-4.0, -12.0)["agree"]
    assert check_offset(-4.0, -12.0)["verdict"] == "DISAGREE"
    assert check_offset(None, -7.0)["verdict"] == "not comparable"
    assert check_offset(-4.0, None)["verdict"] == "not comparable"
    # tolerance is a parameter, not a constant
    assert not check_offset(-4.0, -7.0, tolerance_ppm=1.0)["agree"]


def test_role_comes_from_the_file_name_or_the_sample_type():
    """Two independent rules, because they fail in different situations. The batch generator now
    writes `Std Bracket` as the Xcalibur Sample Type, but sequences already acquired carry
    `Unknown` for the very same injections and are only identifiable by name."""
    from lipidloop.blanks import infer_role
    # named, but the sequence predates the convention
    assert infer_role("Std_Mix_Pos_01", "Unknown") == "standard"
    assert infer_role("Std_mix_Neg_02", "Unknown") == "standard"
    # declared, but generically named
    for declared in ("Std Bracket", "Std_Bracket", "StdBracket", "Standard"):
        assert infer_role("Inj_042", declared) == "standard", declared
    # neither
    assert infer_role("Inj_042", "Unknown") == "sample"
    # the more specific claim wins over standard
    assert infer_role("Blank_Std_01", "Blank") == "blank"
    assert infer_role("Pool_01", "QC") == "qc"


def test_sample_types_are_read_past_the_xcalibur_preamble(tmp_path):
    """Xcalibur writes `Bracket Type=4` above the header row, so a plain DictReader on line 1
    reads the preamble as column names and returns nothing usable."""
    from lipidloop.blanks import read_sample_types
    path = tmp_path / "sequence.csv"
    path.write_text("Bracket Type=4,,,\n"
                    "Sample Type,File Name,Sample ID,Position\n"
                    "Blank,Blank_Pos_01,GE8,GE8\n"
                    "Std Bracket,Std_Mix_Pos_01,RD3,RD3\n")
    types = read_sample_types(path)
    assert types == {"Blank_Pos_01": "Blank", "Std_Mix_Pos_01": "Std Bracket"}


def test_the_cross_check_survives_a_correction_being_applied():
    """The check compares an identification-derived offset with a standards-derived one. Once
    `mass_offset_ppm` is applied the first becomes a *residual* — the precursors were shifted
    before the search — while the standards are measured on the feature's quant ion, which is
    never shifted. Comparing them directly made a correctly calibrated run report a disagreement
    exactly the size of the correction: the check failed precisely when it had worked.
    """
    from lipidloop.standards import check_offset
    applied, residual, standards = -5.77, -0.01, -7.16
    # what the bug did
    assert not check_offset(residual, standards)["agree"]
    # putting both on the instrument's own scale
    assert check_offset(residual + applied, standards)["agree"]


def test_standards_measure_the_instrument_not_the_correction():
    """The fourth time a correction bent the thing meant to check it.

    Once `mass_offset_ppm` began correcting feature masses, measuring the standards on the
    corrected mass meant measuring the residual — the standards-derived offset moved from -7.16 to
    -1.39 the instant the correction was applied. The standards exist to judge the instrument
    independently of the software, so they read the measured mass where one was kept.
    """
    from lipidloop.standards import measured_mz

    class Group:
        quant_ion = 573.3862            # corrected
        quant_ion_measured = 573.3903 - 573.3903 * 5.77e-6 * -1   # what the instrument said
    assert measured_mz(Group()) == Group.quant_ion_measured

    class Uncorrected:
        quant_ion = 573.3862
        quant_ion_measured = None
    assert measured_mz(Uncorrected()) == 573.3862, "no correction applied, so the mass stands"
