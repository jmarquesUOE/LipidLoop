"""Two identified rows that are two adducts of one molecule.

The mirror image of split_peaks.py, and the contrast is the point of both. A split peak is one
molecule divided between two rows and is SUMMED; an adduct pair is one molecule and one ghost and
the ghost is REMOVED. Correlation means opposite things in the two cases: adducts of one molecule
must track, while the halves of a split peak are anti-correlated.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.adduct_pairs import (AdductPairFilter, remove_adduct_pairs)   # noqa: E402
from lipidloop.adducts import Adduct   # noqa: E402
from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample   # noqa: E402

FWHM = 0.09
PROTON = 1.007276
SODIUM = 22.989218

ADDUCTS = [Adduct(name="[M+H]+", formula="H1", loss=False, polarity="+", charge=1, mass=PROTON),
           Adduct(name="[M+Na]+", formula="Na1", loss=False, polarity="+", charge=1, mass=SODIUM),
           Adduct(name="[M+NH4]+", formula="N1H4", loss=False, polarity="+", charge=1,
                  mass=18.033823)]


def samples(n=10):
    return [Sample(file=f"S{i}", role="sample") for i in range(n)]


def group(name, lipid_class, adduct, quant, retention, areas, dot=900):
    g = CompoundGroup(name=name, mw=quant - PROTON, retention=retention,
                      max_area=max(areas), has_ms2=True, areas=list(areas))
    g.quant_ion, g.quant_polarity = quant, "+"
    lipid = Lipid(retention=retention, precursor=quant, sample=Sample(file="x"), dot=dot,
                  rev_dot=980, lipid_string=f"{name} {adduct};", lib_precursor=quant, purity=99,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    lipid.lipid_class = lipid_class
    g.lipid_candidates = [LipidCandidate(lipid)]
    g.final_lipid_id = g.lipid_candidates[0]
    g.purity, g.sum_id = 99.0, name
    return g


TRACKING = [10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0]
ALSO = [12.0, 21.0, 33.0, 44.0, 48.0, 61.0, 75.0, 79.0, 95.0, 99.0]
UNRELATED = [90.0, 10.0, 80.0, 20.0, 70.0, 30.0, 60.0, 40.0, 50.0, 55.0]


def real_pair(areas_b=ALSO, name_b="PE 42:4", dot_b=900):
    """`PE 40:1` as [M+H]+, and a row whose mass is that same molecule's sodium adduct — which
    is also, to 3 ppm, protonated PE 42:4. The coincidence this module exists for."""
    neutral = 801.6236
    a = group("PE 40:1", "PE", "[M+H]+", neutral + PROTON, 5.00, TRACKING)
    b = group(name_b, "PE", "[M+Na]+", neutral + SODIUM, 5.02, areas_b, dot=dot_b)
    return [a, b]


def test_a_confirmed_adduct_pair_loses_its_ghost():
    groups = real_pair()
    report = remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.considered == 1 and report.removed == 1
    assert groups[0].keep and not groups[1].keep
    assert groups[1].filter_reason.startswith("Adduct pair: the [M+Na]+ of PE 40:1")


def test_areas_are_never_summed():
    """The distinction from split_peaks.py. Adding a ghost to a real measurement would corrupt
    both, and the two ions have different response factors anyway."""
    groups = real_pair()
    before = list(groups[0].areas)
    remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert groups[0].areas == before


def test_mass_agreement_alone_is_not_enough():
    """412 pairs on the real data are mass-consistent and only 14 track. Without the correlation
    gate this filter would delete real co-eluting lipids by the hundred."""
    groups = real_pair(areas_b=UNRELATED)
    report = remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.considered == 1 and report.removed == 0
    assert report.rejected_by_correlation == 1
    assert all(g.keep for g in groups)


def test_the_canonical_adduct_survives_even_on_a_lower_dot_product():
    """A sodiated identification is the weaker claim on chemistry — the species barely fragments —
    so a dot product earned on few fragments does not outrank one earned on many."""
    groups = real_pair(dot_b=999)
    remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert groups[0].keep and not groups[1].keep


def test_the_same_name_on_a_different_ion_is_handled_here_not_by_split_peaks():
    """Split-peak merging requires the SAME adduct, because it sums and two ions have two
    response factors. So one molecule on two different ions cannot go there — it belongs here,
    where the redundant row is removed rather than added in."""
    groups = real_pair(name_b="PE 40:1")
    report = remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.removed == 1 and report.same_molecule == 1
    assert groups[0].keep and not groups[1].keep


def test_rows_that_do_not_co_elute_are_untouched():
    groups = real_pair()
    groups[1].retention = 6.5
    report = remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.considered == 0 and all(g.keep for g in groups)


def test_two_unrelated_masses_are_not_a_pair():
    a = group("PE 40:1", "PE", "[M+H]+", 802.6309, 5.00, TRACKING)
    b = group("PC 34:1", "PC", "[M+H]+", 760.5851, 5.02, ALSO)
    report = remove_adduct_pairs([a, b], samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.considered == 0 and a.keep and b.keep


def test_disabled_does_nothing():
    groups = real_pair()
    report = remove_adduct_pairs(groups, samples(), AdductPairFilter(enabled=False), FWHM,
                                 adducts=ADDUCTS)
    assert report.removed == 0 and all(g.keep for g in groups)


def test_an_isotope_spacing_is_never_called_an_adduct():
    """Found by auditing a real run: `SM d43:2` was removed against `HexCer[NS] d43:1` on a gap of
    1.0039 — a 13C spacing. `[M]+` and `[M+H]+` differ by 1.00728, only 3.9 mDa away, so an M+1
    isotope that slipped past the isotope filter looks exactly like an adduct relationship."""
    neutral = 800.0
    a = group("HexCer[NS] d43:1", "HexCer", "[M+H]+", neutral + PROTON, 5.00, TRACKING)
    b = group("SM d43:2", "SM", "[M+H]+", neutral + PROTON + 1.0033548, 5.01, ALSO)
    report = remove_adduct_pairs([a, b], samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.rejected_isotope == 1 and report.removed == 0
    assert a.keep and b.keep


def test_the_row_removed_is_the_one_the_masses_accuse():
    """The direction must come from the mass evidence, not from ranking the pair independently.

    `explained` means "second is an adduct of first" and licenses removing second — nothing else.
    Ranking separately and deleting whoever scored worse removed rows on a relationship that did
    not describe them, and left the reason reading "the None of ...".
    """
    groups = real_pair()
    remove_adduct_pairs(groups, samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    dropped = [g for g in groups if not g.keep]
    assert len(dropped) == 1
    assert "None" not in dropped[0].filter_reason
    assert dropped[0].filter_reason.startswith("Adduct pair: the [M+Na]+ of ")


def test_an_accused_row_that_is_the_stronger_one_is_left_alone():
    """If the masses accuse the canonical, better-supported row, the pair is not what it looks
    like. Skipping costs a true positive; removing costs a real lipid."""
    neutral = 801.6236
    # the [M+Na]+ row comes first, so the mass relationship points at the [M+H]+ row
    a = group("PE 42:4", "PE", "[M+Na]+", neutral + SODIUM, 5.00, TRACKING, dot=999)
    a.quant_ion = neutral + SODIUM
    b = group("PE 40:1", "PE", "[M+H]+", neutral + SODIUM - PROTON + PROTON, 5.01, ALSO)
    b.quant_ion = neutral + PROTON
    report = remove_adduct_pairs([a, b], samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.removed == 0 or all("None" not in g.filter_reason for g in (a, b))


def test_one_molecule_on_two_ions_is_removed_but_labelled_differently():
    """`TG 18:2_18:1_20:1` and `TG 56:4` are one molecule at two resolutions, separated by the
    NH4-to-Na spacing of 4.9554. Both names are right and one ROW is redundant — which is not the
    same claim as one name being wrong, so the reason must not say it is."""
    neutral = 900.0
    a = group("TG 56:4", "TG", "[M+NH4]+", neutral + 18.033823, 5.00, TRACKING)
    b = group("TG 18:2_18:1_20:1", "TG", "[M+Na]+", neutral + SODIUM, 5.01, ALSO)
    report = remove_adduct_pairs([a, b], samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.removed == 1 and report.same_molecule == 1
    assert a.keep and not b.keep
    assert "the same molecule on its [M+Na]+" in b.filter_reason


def test_the_same_molecule_on_the_same_ion_is_still_left_to_split_peaks():
    """Compared on sum composition, so a molecular name and its sum collapse. Comparing the
    displayed strings let `TG 18:2_18:1_20:1` past a guard meant to catch `TG 56:4`."""
    neutral = 900.0
    a = group("TG 56:4", "TG", "[M+NH4]+", neutral + 18.033823, 5.00, TRACKING)
    b = group("TG 18:2_18:1_20:1", "TG", "[M+NH4]+", neutral + 18.033823 + 0.0001, 5.01, ALSO)
    report = remove_adduct_pairs([a, b], samples(), AdductPairFilter(), FWHM, adducts=ADDUCTS)
    assert report.considered == 0 and a.keep and b.keep
