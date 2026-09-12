"""Merging one peak that was integrated as two rows.

The rule is not "close and correlated" — correlation is backwards here. A split peak divides a
conserved total and the integration boundary moves between injections, so the two halves are
uncorrelated or anti-correlated; high correlation instead marks two genuinely distinct
co-regulated species, which must be left alone. What decides it is whether summing improves
precision in the pooled QCs, which are the same material in every vial.
"""
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample   # noqa: E402
from lipidloop.split_peaks import SplitPeakMerge, merge_split_peaks   # noqa: E402

FWHM = 0.09


def samples(n_pools=4, n_samples=4):
    out = [Sample(file=f"S{i}", role="sample") for i in range(n_samples)]
    out += [Sample(file=f"QC{i}", role="qc") for i in range(n_pools)]
    return out


def group(name, adduct, retention, quant_ion, areas, sample):
    lipid = Lipid(retention=retention, precursor=quant_ion, sample=sample, dot=900, rev_dot=980,
                  lipid_string=f"{name} {adduct};", lib_precursor=quant_ion, purity=99,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    out = CompoundGroup(name=name, mw=quant_ion - 1.007, retention=retention,
                        max_area=max(areas), has_ms2=True, areas=list(areas))
    out.quant_ion = quant_ion
    out.lipid_candidates = [LipidCandidate(lipid)]
    out.final_lipid_id = out.lipid_candidates[0]
    out.purity = 99.0
    out.sum_id = name
    return out


def split_pair(pools_a, pools_b, name="PC 16:0_18:1", adduct="[M+H]+",
               retention=(5.00, 5.10), quant=(760.5850, 760.5854)):
    """Two rows, four samples then four QCs. Sample areas are held constant so only the QC
    behaviour is under test."""
    sample_areas = [100.0, 100.0, 100.0, 100.0]
    s = samples()
    a = group(name, adduct, retention[0], quant[0], sample_areas + list(pools_a), s[0])
    b = group(name, adduct, retention[1], quant[1], sample_areas + list(pools_b), s[0])
    return [a, b], s


def test_a_split_peak_is_merged_and_the_areas_are_summed():
    """The halves swap which is larger between injections — the total is stable, each half is
    not. Summing must recover the stable total."""
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0])
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.merged == 1
    kept = [g for g in groups if g.keep]
    assert len(kept) == 1
    assert kept[0].areas[4:] == [100.0, 100.0, 100.0, 100.0]
    assert kept[0].areas[:4] == [200.0, 200.0, 200.0, 200.0]   # sample areas summed too
    dropped = next(g for g in groups if not g.keep)
    assert dropped.filter_reason.startswith("Split peak, merged into PC 16:0_18:1")


def test_two_stable_co_eluting_rows_are_left_alone():
    """Both halves already precise. Summing cannot improve on that, so there is no evidence of a
    split and nothing may be merged — this is the isomer case."""
    groups, s = split_pair([100.0, 101.0, 99.0, 100.0], [50.0, 51.0, 49.0, 50.0])
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.merged == 0
    assert all(g.keep for g in groups)


def test_the_larger_row_survives_and_keeps_its_retention():
    groups, s = split_pair([10.0, 90.0, 20.0, 80.0], [90.0, 10.0, 80.0, 20.0])
    groups[1].max_area = 9999.0        # the second row holds more of the peak
    merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    survivor = next(g for g in groups if g.keep)
    assert survivor.retention == 5.10


def test_different_adducts_are_never_summed():
    """`[M+H]+` and `[M+H-H2O]+` are one molecule and two ions with two response factors.
    Summing them adds two scales together, so the adduct must match even though the name does."""
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0])
    groups[1].lipid_candidates[0].identifications[0].adduct = "[M+H-H2O]+"
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.merged == 0 and all(g.keep for g in groups)


def test_different_names_are_never_summed():
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0])
    groups[1].sum_id = "PE 16:0_18:1"
    groups[1].lipid_candidates[0].lipid_name = "PE 16:0_18:1"
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.merged == 0 and all(g.keep for g in groups)


def test_rows_further_apart_than_the_window_are_not_merged():
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0],
                           retention=(5.0, 5.9))       # 10 peak widths apart
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.merged == 0 and report.considered == 0


def test_a_mass_mismatch_blocks_the_merge():
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0],
                           quant=(760.5850, 760.7000))   # 151 ppm apart
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.merged == 0 and report.considered == 0


def test_nothing_is_merged_without_enough_pooled_qcs():
    """The whole criterion is QC precision. With too few QCs there is no evidence, and guessing
    would change quantification on exactly the runs least able to show it had gone wrong."""
    groups, _ = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0])
    two_pools = [Sample(file=f"S{i}", role="sample") for i in range(4)]
    two_pools += [Sample(file=f"QC{i}", role="qc") for i in range(2)]
    for g in groups:
        g.areas = g.areas[:6]
    report = merge_split_peaks(groups, two_pools, SplitPeakMerge(), FWHM)
    assert report.skipped_no_pools and report.merged == 0
    assert all(g.keep for g in groups)


def test_three_fragments_of_one_peak_collapse_to_one_row():
    """An absorbed row must not still be available to absorb a third, or one peak in three parts
    becomes a merged pair plus an orphan."""
    s = samples()
    areas = [[100.0] * 4 + p for p in ([60.0, 5.0, 55.0, 10.0],
                                       [5.0, 60.0, 10.0, 55.0],
                                       [35.0, 35.0, 35.0, 35.0])]
    groups = [group("PC 16:0_18:1", "[M+H]+", 5.0 + 0.05 * i, 760.585 + 0.0001 * i, a, s[0])
              for i, a in enumerate(areas)]
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert len([g for g in groups if g.keep]) == 1
    assert report.merged == 2
    survivor = next(g for g in groups if g.keep)
    assert survivor.areas[4:] == [100.0, 100.0, 100.0, 100.0]


def test_disabled_does_nothing():
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0])
    report = merge_split_peaks(groups, s, SplitPeakMerge(enabled=False), FWHM)
    assert report.merged == 0 and all(g.keep for g in groups)


def test_features_found_counts_the_union():
    """A merged row is present wherever either half was, so `Features Found` must not report
    only the winning half's injections."""
    from lipidloop.peaks import Compound
    groups, s = split_pair([90.0, 10.0, 80.0, 20.0], [10.0, 90.0, 20.0, 80.0])
    groups[0].compounds = [Compound(mw=759.6, retention=5.0, fwhm=FWHM, max_mi=0, n_adducts=1,
                                    area=1.0, sample=s[0])]
    groups[1].compounds = [Compound(mw=759.6, retention=5.1, fwhm=FWHM, max_mi=0, n_adducts=1,
                                    area=1.0, sample=s[1])]
    merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    survivor = next(g for g in groups if g.keep)
    assert len(survivor.compounds) == 2


def test_a_row_undetected_in_the_qcs_does_not_merge_on_no_evidence():
    """Found by running it: an undetected row scores an infinite CV, and `after * factor <= inf`
    is true whatever `after` is, so the pair merged on nothing and landed at 283%. Missing
    evidence must block the merge, not license it — `inf` means unmeasurable, not terrible."""
    groups, s = split_pair([0.0, 0.0, 0.0, 0.0], [10.0, 90.0, 20.0, 80.0])
    report = merge_split_peaks(groups, s, SplitPeakMerge(), FWHM)
    assert report.considered == 1 and report.merged == 0
    assert all(g.keep for g in groups)


def test_a_row_seen_in_too_few_qcs_does_not_merge():
    """Two detections out of four is not a precision estimate. A genuine split is present in the
    QCs on both sides of the cut, so requiring detection is chemistry as well as statistics."""
    groups, s = split_pair([50.0, 60.0, 0.0, 0.0], [10.0, 90.0, 20.0, 80.0])
    report = merge_split_peaks(groups, s, SplitPeakMerge(min_pools=3), FWHM)
    assert report.merged == 0 and all(g.keep for g in groups)


def _systematic_pair(ratio_wanders: bool):
    """Ten injections: six samples then four pooled QCs (that is `samples()`'s order).

    The pools hold a CONSTANT ratio in both variants. That matters: the pools are one homogenate,
    so two independent molecules sit at fixed levels there too, and each row is precise on its own
    — which is what stops the precision test from firing and leaves the ratio across the SAMPLES as
    the only thing separating one divided peak from two co-varying molecules.
    """
    people = samples(n_pools=4, n_samples=6)
    areas_a, areas_b = [], []
    for i, person in enumerate(people):
        true = 1000.0 * (1 + 0.1 * i)
        wander = ratio_wanders and person.role == "sample"
        share = 0.45 * (1 + 0.6 * (i % 3)) if wander else 0.45
        areas_a.append(true * 0.55)
        areas_b.append(true * share)
    a = group("Cer[NDS] d41:2", "[M+FA-H]-", 10.518, 678.6037, areas_a, people[0])
    b = group("Cer[NDS] d41:2", "[M+FA-H]-", 10.629, 678.6036, areas_b, people[0])
    return people, a, b


def test_a_peak_cut_at_the_same_point_in_every_file_is_merged():
    """The other kind of split, which the precision test cannot see.

    A peak whose apex wanders across the boundary is anti-correlated and summing plainly improves
    precision. A peak cut at the SAME point in every file has both halves scaling with the true
    abundance: each is precise on its own, so summing improves nothing measurable, and the pair is
    positively correlated at a near-constant ratio. `Cer[NDS] d41:2` on the study behind this.
    """
    people, a, b = _systematic_pair(ratio_wanders=False)
    report = merge_split_peaks([a, b], people, SplitPeakMerge(min_pools=3), avg_fwhm=0.089)
    assert report.merged == 1 and report.systematic == 1
    assert not b.keep and "merged into" in b.filter_reason
    assert a.areas[0] == pytest.approx(1000.0)


def test_two_molecules_that_merely_co_vary_are_not_merged():
    """The ratio is what separates them. Co-varying molecules reach high correlation easily —
    TG 61:3 on the same run correlates at +0.95 — but their ratio moves, 50.7% there."""
    people, a, b = _systematic_pair(ratio_wanders=True)
    report = merge_split_peaks([a, b], people, SplitPeakMerge(min_pools=3), avg_fwhm=0.089)
    assert report.merged == 0, "a moving ratio is two molecules, not one peak"
    assert a.keep and b.keep
