"""Tests for automatic artefact screening.

The whole risk of screening automatically is deleting something real, so most of these pin the
things that must NOT be removed.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.artefacts import plausible_mass, screen, strip   # noqa: E402


class Spectrum:
    def __init__(self, mz, intensity):
        self.mz = list(mz)
        self.intensity = list(intensity)


def test_real_lipid_fragments_are_all_possible_masses():
    """Every one of these is a genuine fragment; none may be judged impossible."""
    for mz in (184.0733,   # phosphocholine
               104.1070,   # choline
               264.2686,   # sphingosine
               369.3516,   # cholestadiene
               255.2330,   # 16:0 acyl anion
               290.0881):  # sialic acid
        assert plausible_mass(mz), mz


def test_the_artefact_is_an_impossible_mass():
    """0.290 defect at m/z 178; the richest possible CHNOPS ion there reaches about 0.23."""
    assert not plausible_mass(178.295)
    assert not plausible_mass(178.264)


def test_detects_a_ubiquitous_impossible_base_peak():
    spectra = [Spectrum([178.295, 200.0 + i], [999.0, 100.0]) for i in range(50)]
    found = screen(spectra)
    assert len(found) == 1
    assert found[0].mz == pytest.approx(178.295, abs=0.01)
    assert found[0].fraction == 1.0


def test_never_removes_a_ubiquitous_REAL_fragment():
    """A class fragment can be in every spectrum and the base peak — and must survive.

    This is the failure that would matter: phosphocholine dominating a PC-only run.
    """
    spectra = [Spectrum([184.0733, 200.0 + i], [999.0, 100.0]) for i in range(50)]
    assert screen(spectra) == []


def test_requires_base_peak_dominance_not_just_presence():
    """Present everywhere but always minor: not damaging, so not removed."""
    spectra = [Spectrum([178.295, 200.0 + i], [1.0, 999.0]) for i in range(50)]
    assert screen(spectra) == []


def test_requires_near_universal_presence():
    spectra = ([Spectrum([178.295], [999.0]) for _ in range(40)]
               + [Spectrum([300.0], [999.0]) for _ in range(60)])
    assert screen(spectra) == []


def test_strip_removes_only_the_artefact():
    spectra = [Spectrum([178.295, 184.0733, 264.2686], [999.0, 500.0, 300.0])]
    removed = strip(spectra, screen([Spectrum([178.295, 184.0733], [999.0, 10.0])
                                     for _ in range(50)]))
    assert removed == 1
    assert spectra[0].mz == [184.0733, 264.2686]
    assert spectra[0].intensity == [500.0, 300.0]


def test_strip_with_nothing_found_is_a_no_op():
    spectra = [Spectrum([184.0733], [999.0])]
    assert strip(spectra, []) == 0
    assert spectra[0].mz == [184.0733]
