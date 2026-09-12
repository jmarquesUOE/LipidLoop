"""Tests for the peak finder's building blocks.

The peak finder joins two things that are measured differently — an MS2's retention time and a
feature's aligned retention time — so most of what can go wrong is arithmetic about time and
about names, which is what these cover.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.loess import loess                                     # noqa: E402
from lipidloop.peaks import (GaussianModel, Sample, sum_composition,   # noqa: E402
                              CompoundGroup, Lipid, ppm_diff)
from lipidloop.peakfinder import parse_masses, parse_purity           # noqa: E402


def test_sum_composition_collapses_chains():
    assert sum_composition("PC 16:0_18:1") == "PC 34:1"
    assert sum_composition("TG 18:3_18:2_18:1") == "TG 54:6"
    assert sum_composition("Plasmenyl-PE P-18:0_20:4") == "Plasmenyl-PE P-38:4"
    assert sum_composition("PC[OH] OH-14:1_18:4") == "PC[OH] OH-32:5"


def test_sum_composition_leaves_unparseable_names_alone():
    assert sum_composition("Cholesterol") == "Cholesterol"
    assert sum_composition("SP d18:1") == "SP d18:1"


def test_sum_composition_does_not_truncate_a_chain_it_cannot_fully_parse():
    """★ REGRESSION. `LipiDex2_FAHFA` names its chains `18:1-(O-18:0)`, not this parser's `_`.

    With `re.match` (not `fullmatch`) the leading `18:1` matched and everything after it —
    `-(O-18:0)`, the SECOND chain — was silently dropped, returning `FAHFA 18:1`: a name for a
    563 Da, two-chain, chain-unresolved molecule that is byte-identical to a genuine single free
    fatty acid's sum composition. On Skin_QEplus this is exactly what displaced five real `FA`
    identifications from the delivered table (2026-08-28 displacement test) — not because FAHFA
    out-competed them on evidence (every winning FAHFA match carried Purity 0), but because the
    truncated name looked like the thing it had overwritten.

    A part the regex cannot consume in full must fall through to the documented behaviour: the
    whole name is left alone, not silently shortened to whatever prefix happened to match.
    """
    assert sum_composition("FAHFA 18:1-(O-18:0)") == "FAHFA 18:1-(O-18:0)"
    assert sum_composition("FAHFA 20:0-(O-22:0)") == "FAHFA 20:0-(O-22:0)"
    # A single well-formed chain must still collapse normally — the fix must not make fullmatch
    # reject anything the existing format actually produces.
    assert sum_composition("PC 16:0_18:1") == "PC 34:1"


def test_top_ether_reports_the_neutral_notation():
    """Where plasmenyl and plasmanyl cannot be told apart, neither is asserted."""
    group = CompoundGroup(name="", mw=0.0, retention=0.0, max_area=0.0, has_ms2=False, areas=[])
    group.sum_id = "Plasmenyl-PE P-40:7"
    assert group.top_ether() == "PE O-40:7"


def test_gaussian_peaks_at_the_apex_and_falls_away():
    model = GaussianModel(fwhm=0.07, area=1000.0, apex_rt=10.0)
    assert model.normalized_height(10.0) == pytest.approx(1.0, abs=0.02)
    assert model.normalized_height(10.035) < model.normalized_height(10.0)
    assert model.normalized_height(10.035) == pytest.approx(0.5, abs=0.1)   # half maximum
    assert model.normalized_height(99.0) == 0.001                          # outside the peak


def test_gaussian_height_scales_with_area():
    narrow = GaussianModel(fwhm=0.05, area=1000.0, apex_rt=1.0)
    wide = GaussianModel(fwhm=0.20, area=1000.0, apex_rt=1.0)
    assert narrow.height > wide.height     # same area spread over more time is shorter


def test_loess_follows_a_trend_and_smooths_noise():
    x = [i / 100 for i in range(200)]
    y = [0.5 * v for v in x]
    fitted = loess(x, y)
    assert fitted[10] == pytest.approx(y[10], abs=0.01)
    assert fitted[150] == pytest.approx(y[150], abs=0.01)

    noisy = list(y)
    noisy[100] += 5.0                       # one wild outlier
    smoothed = loess(x, noisy)
    assert abs(smoothed[100] - y[100]) < 2.5   # robustness pulls it back


def test_loess_handles_degenerate_input():
    assert list(loess([], [])) == []
    assert list(loess([1.0, 2.0], [3.0, 4.0])) == [3.0, 4.0]


def test_corrected_rt_maps_through_the_fit():
    sample = Sample(file="a.raw")
    sample.rt_deviations = [(0.1, t / 10, 1.0) for t in range(1, 60)]
    sample.fit()
    corrected = sample.corrected_rt(2.0)
    assert corrected == pytest.approx(1.9, abs=0.05)


def test_corrected_rt_returns_zero_past_the_fitted_range():
    """Past the last fitted point the original returns 0, so an ID simply fails to associate."""
    sample = Sample(file="a.raw")
    sample.rt_deviations = [(0.1, 1.0, 1.0), (0.1, 2.0, 1.0)]
    sample.fit()
    assert sample.corrected_rt(99.0) == 0.0


def test_spectral_components_round_trip():
    parsed = parse_purity("PC 34:1 [M+H]+(82) / PC 34:1 [M+Na]+(18)")
    assert parsed == [("PC 34:1 [M+H]+", 82), ("PC 34:1 [M+Na]+", 18)]
    assert parse_purity("") == []


def test_potential_fragments_parse():
    assert parse_masses("619.4698 | 621.4855 | ") == [619.4698, 621.4855]
    assert parse_masses("") == []


def test_ppm_diff_is_relative_to_the_reference():
    assert ppm_diff(500.001, 500.0) == pytest.approx(2.0, abs=0.01)


def test_lipid_class_is_capitalised():
    """`d5TG` in the library is reported as `D5TG` — parseName upper-cases the first letter."""
    sample = Sample(file="a.raw")
    lipid = Lipid(retention=5.0, precursor=886.789, sample=sample, dot=900.0, rev_dot=990.0,
                  lipid_string="d5TG 16:0_18:1_18:1 [M+NH4]+;", lib_precursor=886.789,
                  purity=90, is_lipidex=True, purity_array=[], fragment_masses=[])
    assert lipid.lipid_class == "D5TG"
    assert lipid.lipid_name == "D5TG 16:0_18:1_18:1"
    assert lipid.sum_lipid_name == "D5TG 52:2"


def test_lipid_name_drops_the_adduct_but_class_survives():
    sample = Sample(file="a.raw")
    lipid = Lipid(retention=5.0, precursor=760.585, sample=sample, dot=900.0, rev_dot=990.0,
                  lipid_string="PC 16:0_18:1 [M+H]+;", lib_precursor=760.5856, purity=90,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    assert lipid.lipid_name == "PC 16:0_18:1"
    assert lipid.sum_lipid_name == "PC 34:1"
    assert lipid.lipid_class == "PC"
    assert lipid.adduct == "[M+H]+"
    assert lipid.polarity == "+"


def test_redundant_groups_merge_on_sum_composition():
    """Two rows reporting the same molecule at the same time are one peak, not two.

    Matched on sum composition rather than the displayed name: a row reported as `PC 34:1` and
    one reported as `PC 16:0_18:1` differ in reporting resolution, not in molecule.
    """
    from lipidloop.peaks import Compound, LipidCandidate

    def with_id(area, retention, name):
        g = CompoundGroup(name="", mw=760.0, retention=retention, max_area=area,
                          has_ms2=True, areas=[area])
        g.avg_fwhm = 0.07
        sample = Sample(file="a.raw")
        lipid = Lipid(retention=retention, precursor=760.585, sample=sample, dot=900.0,
                      rev_dot=990.0, lipid_string=f"{name} [M+H]+;", lib_precursor=760.5856,
                      purity=90, is_lipidex=True, purity_array=[], fragment_masses=[])
        g.compounds = [Compound(mw=760.0, retention=retention, fwhm=0.07, max_mi=1,
                                n_adducts=1, area=area, sample=sample)]
        g.lipid_candidates = [LipidCandidate(lipid)]
        g.final_lipid_id = g.lipid_candidates[0]
        g.sum_id = lipid.sum_lipid_name
        return g

    big = with_id(1000.0, 8.00, "PC 16:0_18:1")
    small = with_id(100.0, 8.02, "PC 34:1")          # same molecule, different resolution
    assert big.final_lipid_id.sum_lipid_name == small.final_lipid_id.sum_lipid_name

    big.merge_compound_group(small)
    assert not small.keep
    assert small.filter_reason == "Redundant Identification"
    assert len(big.compounds) == 2       # evidence pooled, not discarded
    assert big.keep


def test_fwhm_window_divides_rather_than_multiplies():
    """The argument is a divisor. Reading it as a multiplier inverts every window in the sweep."""
    g = CompoundGroup(name="", mw=760.0, retention=8.0, max_area=1.0, has_ms2=True, areas=[])
    g.avg_fwhm = 0.10

    def at(rt):
        return CompoundGroup(name="", mw=760.0, retention=rt, max_area=1.0,
                             has_ms2=True, areas=[])

    # divisor 0.5 -> two peak widths, so 0.15 min away is still inside
    assert g.target_in_fwhm(at(8.15), 0.5)
    assert not g.target_in_fwhm(at(8.25), 0.5)
    # divisor 2.0 -> half a peak width, so 0.15 min away is well outside
    assert not g.target_in_fwhm(at(8.15), 2.0)
    assert g.target_in_fwhm(at(8.04), 2.0)
    # the redundant window, avg_fwhm/0.823, sits between the two
    assert g.target_in_fwhm(at(8.10), 0.823)
    assert not g.target_in_fwhm(at(8.20), 0.823)


def _group_with_candidates(reported_at_sum: bool):
    """A group carrying two molecular candidates that share one sum composition."""
    from lipidloop.peaks import Compound, LipidCandidate

    g = CompoundGroup(name="", mw=760.0, retention=8.0, max_area=1e6, has_ms2=True, areas=[])
    candidates = []
    for molecular, files in (("PC 16:0_18:1", ["a.raw", "b.raw"]), ("PC 14:0_20:1", ["a.raw"])):
        lipids = []
        for f in files:
            lipids.append(Lipid(retention=8.0, precursor=760.585, sample=Sample(file=f),
                                dot=900.0, rev_dot=990.0, lipid_string=f"{molecular} [M+H]+;",
                                lib_precursor=760.5856, purity=90, is_lipidex=True,
                                purity_array=[], fragment_masses=[]))
        candidate = LipidCandidate(lipids[0])
        for extra in lipids[1:]:
            candidate.add(extra)
        candidates.append(candidate)
    g.lipid_candidates = candidates
    g.final_lipid_id = candidates[0]
    g.sum_id = "PC 34:1"
    g.purity = 40.0 if reported_at_sum else 90.0
    return g


def test_ms2_support_counts_spectra_and_files():
    g = _group_with_candidates(reported_at_sum=False)
    assert g.identification()[0] == "PC 16:0_18:1"      # molecular, purity cleared
    assert g.ms2_support() == (2, 2)                    # only that molecule's spectra


def test_ms2_support_for_a_sum_composition_counts_every_contributing_spectrum():
    """The mistake this exists to prevent: a sum-composition row looking unsupported."""
    g = _group_with_candidates(reported_at_sum=True)
    assert g.identification()[0] == "PC 34:1"           # sum, purity below threshold
    assert g.ms2_support() == (3, 2)                    # all three spectra, from two files


def test_ms2_support_is_zero_for_a_retention_model_assignment():
    g = CompoundGroup(name="", mw=760.0, retention=8.0, max_area=1.0, has_ms2=False, areas=[])
    g.rtls_identification = "PC 34:1"
    assert g.identification()[0] == "PC 34:1"
    assert g.identification_source() == "RT model"
    assert g.ms2_support() == (0, 0)


def test_ms2_support_is_zero_when_unidentified():
    g = CompoundGroup(name="", mw=760.0, retention=8.0, max_area=1.0, has_ms2=False, areas=[])
    assert g.ms2_support() == (0, 0)
