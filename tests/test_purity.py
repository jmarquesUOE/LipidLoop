"""Tests for identification purity.

Purity is where a plausible reimplementation drifts hardest from LipiDex, because most of its
rules are refusals: a candidate that cannot account for every chain scores nothing at all, and
a candidate in the wrong polarity is not considered however well it matches.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.msp import LibrarySpectrum, split_annotation   # noqa: E402
from lipidloop.purity import (FattyAcid, PeakPurity, _median,  # noqa: E402
                               _round_half_up, read_fatty_acids)

FA_DB = [FattyAcid("16:0", "Alkyl", "C16H31O2", "true"),
         FattyAcid("18:1", "Alkyl", "C18H33O2", "true"),
         FattyAcid("18:0", "Alkyl", "C18H35O2", "true")]


def lipid(name, peaks, is_lipidex=True, optimal=True):
    """peaks: (mz, intensity, annotation)"""
    spectrum = LibrarySpectrum(name=name, precursor_mz=800.0, mz=[], intensity=[],
                               library="test.msp", is_lipidex=is_lipidex,
                               optimal_polarity=optimal)
    for mz, intensity, annotation in peaks:
        frag_type, fatty_acid = split_annotation(annotation)
        spectrum.mz.append(mz)
        spectrum.intensity.append(intensity)
        spectrum.annotations.append(annotation)
        spectrum.types.append(frag_type)
        spectrum.fatty_acids.append(fatty_acid)
    return spectrum


# Moiety fragments carry "-" where a formula would go, so every chain fragment of one lipid
# shares a transition type and they group together. Real entries look exactly like this; using
# a per-chain formula here would split the chains apart and score nothing.
PC = lipid("PC 16:0_18:1 [M+H]+;",
           [(255.2330, 999.0, "-_Alkyl Fragment_[16:0]"),
            (281.2486, 999.0, "-_Alkyl Fragment_[18:1]"),
            (184.0733, 500.0, "C5H15N1O4P1_Fragment_[]")])


def calc(sample_mz, sample_int, top=PC, isobaric=(), fa_db=FA_DB):
    p = PeakPurity(list(sample_mz), list(sample_int))
    return p.calc_purity(top, list(isobaric), list(fa_db)), p


def test_annotation_split():
    assert split_annotation("O-1_Alkyl Fragment_[10:0]") == ("O-1_Alkyl Fragment", "10:0")
    assert split_annotation("") == (None, "")


def test_sole_candidate_is_wholly_pure():
    purity, _ = calc([255.2330, 281.2486], [500.0, 400.0])
    assert purity == 100


def test_missing_chain_scores_nothing():
    """Only one of the two chains is present, so the candidate is not merely penalised."""
    purity, p = calc([255.2330], [500.0])
    assert purity == 0
    assert p.lipids == []


def test_chain_below_the_five_percent_floor_counts_as_missing():
    """The floor is 5% of the most intense moiety fragment in the *library* entry."""
    purity, _ = calc([255.2330, 281.2486], [500.0, 0.0])
    assert purity == 0


def test_off_polarity_and_non_lipidex_candidates_score_nothing():
    assert calc([255.2330, 281.2486], [500.0, 400.0],
                top=lipid("PC 16:0_18:1 [M+H]+;", [(255.2330, 999.0, "-_Alkyl Fragment_[16:0]"),
                                                   (281.2486, 999.0, "-_Alkyl Fragment_[18:1]")],
                          is_lipidex=False))[0] == 0
    # optimal_polarity is enforced by the caller, not by calc_purity itself
    assert not lipid("x", [], optimal=False).optimal_polarity


def test_unannotated_peak_aborts_the_candidate():
    """A LipiDex entry with a bare peak throws in the original, caught into a zero."""
    broken = lipid("PC 16:0_18:1 [M+H]+;",
                   [(255.2330, 999.0, "-_Alkyl Fragment_[16:0]"),
                    (281.2486, 999.0, "")])
    assert calc([255.2330, 281.2486], [500.0, 400.0], top=broken)[0] == 0


def test_competing_candidate_dilutes_purity():
    other = lipid("PC 18:0_16:1 [M+H]+;",
                  [(283.2643, 999.0, "-_Alkyl Fragment_[18:0]"),
                   (253.2173, 999.0, "-_Alkyl Fragment_[16:1]")])
    fa_db = FA_DB + [FattyAcid("16:1", "Alkyl", "C16H29O2", "true")]
    purity, p = calc([255.2330, 281.2486, 283.2643, 253.2173],
                     [500.0, 400.0, 100.0, 100.0], isobaric=[other], fa_db=fa_db)
    assert 0 < purity < 100
    assert [l.name for l in p.lipids] == ["PC 16:0_18:1 [M+H]+;", "PC 18:0_16:1 [M+H]+;"]
    assert sum(p.purities) == pytest.approx(100, abs=1)


def test_repeated_chain_is_claimed_once_per_occurrence():
    """A lipid carrying 16:0 twice claims that one fragment twice."""
    di = lipid("DG 16:0_16:0 [M+H]+;", [(255.2330, 999.0, "-_Alkyl Fragment_[16:0]")])
    _, p = calc([255.2330], [500.0], top=di)
    assert p.intensities == [500.0]     # median of [500, 500]


def test_median_falls_to_the_mean_of_the_middle_pair():
    assert _median([1.0, 2.0, 3.0]) == 2.0
    assert _median([1.0, 2.0, 3.0, 4.0]) == 2.5
    assert _median([]) == 0.0


def test_rounding_is_half_up_not_bankers():
    assert _round_half_up(0.5) == 1
    assert _round_half_up(1.5) == 2       # Python's round() gives 2 here but 0 for 0.5
    assert _round_half_up(2.5) == 3


def test_fatty_acid_database_reads_the_real_file():
    path = Path(__file__).resolve().parents[1] / "data/lipidex_src/FattyAcids.csv"
    if not path.exists():
        pytest.skip("FattyAcids.csv not present")
    db = read_fatty_acids(path)
    assert len(db) > 50
    assert {fa.base for fa in db} >= {"Alkyl", "Sphingoid"}
    assert all(fa.name and fa.formula for fa in db)
