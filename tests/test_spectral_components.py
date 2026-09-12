"""The per-name purity breakdown, written out."""
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.components import write_spectral_components                   # noqa: E402
from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample     # noqa: E402
from lipidloop.peakfinder import PeakFinderResult                            # noqa: E402


def _group(retention, sum_id, breakdown, weight):
    """A group reported at SUM composition, which is the case this file exists for.

    The breakdown matters most where the peak did not resolve to one molecular name — that is
    where a chain distribution is unavailable from the results table and has to come from the
    fragments. Purity below 75 is what makes `identification()` report the sum.
    """
    sample = Sample(file="S1")
    lipid = Lipid(retention=retention, precursor=760.585, sample=sample, dot=900, rev_dot=980,
                  lipid_string=f"{sum_id} [M+H]+;", lib_precursor=760.585, purity=0,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    lipid.lipid_name = lipid.sum_lipid_name = sum_id
    lipid.lipid_class = sum_id.split(" ", 1)[0]
    candidate = LipidCandidate(lipid)
    candidate.lipid_name, candidate.lipid_class = sum_id, lipid.lipid_class

    group = CompoundGroup(name="g", mw=759.578, retention=retention, max_area=1.0,
                          has_ms2=True, areas=[1.0])
    group.sum_id = sum_id
    group.lipid_candidates = [candidate]
    group.summed_purities = sorted(breakdown, key=lambda p: -p[1])
    group.purity_weight = weight
    group.purity = group.summed_purities[0][1] / weight if weight else 0.0
    return group


def _write(tmp_path, groups):
    result = PeakFinderResult(groups, [Sample(file="S1")])
    path = tmp_path / "components.csv"
    counts = write_spectral_components(result, path)
    return list(csv.DictReader(path.open(newline=""))), counts


def test_the_shares_are_a_distribution_and_not_the_purity_denominator(tmp_path):
    """`purity_weight` is the wrong denominator here, and using it looks almost right.

    It counts the Gaussian score once per ENTRY, so it grows with the number of candidates: the
    shares under it sum to a median of 41% on real data and to less the more candidates a peak
    carries, which makes two peaks incomparable. That is correct for purity, whose question is how
    much belongs to the winner, and wrong for a distribution. Normalised within the peak instead.
    """
    # a denominator of 20 while purity_weight is 8 — the two disagree, and the file must use the
    # first, or a four-candidate peak reads as weaker than a two-candidate one at equal evidence
    group = _group(9.3, "PC 38:3",
                   [("PC 18:0_20:3", 16.0), ("PC 18:1_20:2", 4.0)], weight=8.0)
    rows, (written, peaks) = _write(tmp_path, [group])

    assert (written, peaks) == (2, 1)
    shares = {r["Component"]: float(r["Share of Fragment Signal (%)"]) for r in rows}
    assert shares == {"PC 18:0_20:3": 80.0, "PC 18:1_20:2": 20.0}
    assert sum(shares.values()) == 100.0, "a distribution sums to the whole"
    assert shares["PC 18:0_20:3"] != group.purity, "the purity denominator would give 200"


def test_one_combination_under_two_adducts_is_one_component(tmp_path):
    """Otherwise a molecule competes with itself.

    `TG 18:1_18:1_22:4 [M+NH4]+` and `[M+Na]+` are two library entries for one set of chains. Left
    separate they appear as two isomers at half the share each, which both understates the real
    combination and invents a second one — and on real TG data the same chains routinely appear
    under both adducts.
    """
    group = _group(10.5, "TG 58:6",
                   [("TG 18:1_18:1_22:4 [M+NH4]+", 45.0),
                    ("TG 18:1_18:1_22:4 [M+Na]+", 35.0),
                    ("TG 18:2_18:0_22:4 [M+NH4]+", 20.0)], weight=3.0)
    rows, (written, _) = _write(tmp_path, [group])

    assert written == 2, "three entries, two combinations"
    shares = {r["Component"]: float(r["Share of Fragment Signal (%)"]) for r in rows}
    assert shares == {"TG 18:1_18:1_22:4": 80.0, "TG 18:2_18:0_22:4": 20.0}
    assert rows[0]["Rank"] == "1", "re-ranked after merging, not before"


def test_isomers_are_distinguished_from_co_isolated_other_lipids(tmp_path):
    """Two different questions share one breakdown.

    A component collapsing to the same sum composition is a chain isomer, and the share is a
    statement about how this lipid's chains divide. A component from another class is a statement
    about what else was co-isolated. Summing the two into one distribution would turn a
    co-isolation into a chain composition.
    """
    group = _group(9.3, "PE 36:4",
                   [("PE 16:0_20:4", 120.0), ("PE 18:2_18:2", 60.0), ("PC 34:1", 20.0)],
                   weight=2.0)
    rows, _ = _write(tmp_path, [group])
    isomer = {r["Component"]: r["Isomer"] for r in rows}

    assert isomer["PE 16:0_20:4"] == "yes"
    assert isomer["PE 18:2_18:2"] == "yes", "same sum composition by a different division"
    assert isomer["PC 34:1"] == "no", "another class co-isolated at this mass"
    # the chain distribution is the isomer rows renormalised, and it is the caller's job to say so
    chains = {r["Component"]: float(r["Share of Fragment Signal (%)"])
              for r in rows if r["Isomer"] == "yes"}
    assert sum(chains.values()) < 100.0, "the co-isolated lipid holds the rest"


def test_an_undivided_peak_is_not_written(tmp_path):
    """One component is the whole signal and says nothing.

    Every unambiguous peak would otherwise contribute a row reading 100%, which is what the
    `Purity` column in the results table already says.
    """
    groups = [_group(9.3, "PC 34:1", [("PC 16:0_18:1", 200.0)], weight=2.0),
              _group(9.9, "PC 38:3", [("PC 18:0_20:3", 150.0), ("PC 18:1_20:2", 50.0)],
                     weight=2.0)]
    rows, (written, peaks) = _write(tmp_path, groups)

    assert peaks == 1 and written == 2
    assert {r["Identification"] for r in rows} == {"PC 38:3"}


def test_a_group_with_no_weight_is_skipped_rather_than_divided_by_zero(tmp_path):
    """Purity is 0 by design for LipidBlast and off-polarity entries, and then there is no scale."""
    group = _group(9.3, "Cer[ADS] d42:1",
                   [("Cer[ADS] d18:1_24:0", 0.0), ("Cer[ADS] d18:0_24:1", 0.0)], weight=0.0)
    rows, (written, peaks) = _write(tmp_path, [group])
    assert (rows, written, peaks) == ([], 0, 0)


def test_chain_evidence_asks_whether_the_purity_is_about_THIS_row():
    """A high purity does not mean the purity belongs to the row it is printed on.

    `calc_purity` drops a top hit that scores nothing and then returns `purities[0]`, which by then
    belongs to a different candidate, and `search.py` writes it onto the winner. The delivered data
    contains `Cer[ADS] d17:0_17:0` at **Purity 100** whose entire breakdown reads
    `Cer[NP] t18:0_16:0 (100)` — the phyto ceramide was measured, the dihydro one was named.

    An earlier version of this method keyed on `purity_weight` being non-zero. That is group-wide
    and non-zero for every chain-resolved row in that study, so it returned `fragments` for all 166
    and never once said `library name`. Held here because the broken version looked right and
    passed a test that only exercised an empty breakdown.
    """
    measured = _group(10.5, "Cer[ADS] d18:1_24:0", [("Cer[ADS] d18:1_24:0", 200.0)], weight=2.0)

    # purity 100, and every point of it belongs to a different lipid
    borrowed = _group(9.4, "Cer[ADS] d17:0_17:0", [("Cer[NP] t18:0_16:0", 200.0)], weight=2.0)

    # no breakdown at all is a different failure from a breakdown naming someone else
    ineligible = _group(9.8, "PC 16:0_18:1", [], weight=0.0)
    ineligible.summed_purities = []

    summed = _group(9.1, "PC 34:1", [("PC 16:0_18:1", 150.0)], weight=2.0)
    summed.lipid_candidates[0].lipid_name = "PC 34:1"

    assert measured.chain_evidence() == "fragments"
    assert borrowed.chain_evidence() == "library name", "purity 100, but not this row's purity"
    assert ineligible.chain_evidence() == "not eligible"
    assert summed.chain_evidence() == "", "a sum composition asserts no chains to qualify"

    # the trap: purity alone cannot tell the first two apart
    assert borrowed.purity == measured.purity == 100.0
    assert borrowed.purity_weight and measured.purity_weight, "weight is non-zero for both"


def test_chain_key_sorts_chains_but_keeps_the_class_and_the_hydroxyl_position():
    """Two different questions that a naive string compare gets backwards.

    `PC 16:0_18:1` and `PC 18:1_16:0` are one lipid written two ways — sorting handles it. But
    `PC[OH] OH-16:0_18:2` and `OH-18:2_16:0` are NOT one lipid: the prefix marks which chain carries
    the hydroxyl, and the library never writes a pair in both orders. And `Cer[ADS]` against
    `Cer[NP]` at the same chains is the distinction the whole column exists to make.
    """
    from lipidloop.peaks import _chain_key

    assert _chain_key("PC 16:0_18:1") == _chain_key("PC 18:1_16:0"), "order is not identity"
    assert _chain_key("PC[OH] OH-16:0_18:2") != _chain_key("PC[OH] OH-18:2_16:0"), \
        "the hydroxyl is on a different chain — different molecules"
    assert _chain_key("Cer[ADS] d18:0_16:0") != _chain_key("Cer[NP] t18:0_16:0"), "class matters"
    # the adduct and the trailing semicolon a library name carries must not change identity
    assert _chain_key("Cer[ADS] d17:0_17:0 [M+FA-H]-;") == _chain_key("Cer[ADS] d17:0_17:0")


def test_purity_source_names_the_lipid_the_purity_belongs_to():
    """The cross-class case `Chain Evidence` cannot express.

    `calc_purity` drops a top hit that scores nothing and returns `purities[0]`, which by then
    belongs to a different candidate — so `PS 20:4_22:0` can carry Purity 100 computed from a `PC`.
    `Chain Evidence` only asks whether a component names this row's chains and is blank for sum
    compositions, so it cannot say that. 431 positive and 341 negative matches are affected.
    """
    own = _group(9.4, "Cer[ADS] d18:1_24:0", [("Cer[ADS] d18:1_24:0", 200.0)], weight=2.0)
    borrowed = _group(9.4, "Cer[ADS] d17:0_17:0", [("Cer[NP] t18:0_16:0", 200.0)], weight=2.0)

    assert own.purity_source() == "self"
    assert borrowed.purity_source() == "Cer[NP] t18:0_16:0", "names the lipid it came from"

    # a row that declined to name its chains has not been scored by a different molecule — only by
    # a better-resolved view of its own, so it is still `self`
    summed = _group(9.1, "PC 34:1", [("PC 16:0_18:1", 150.0)], weight=2.0)
    summed.lipid_candidates[0].lipid_name = "PC 34:1"
    assert summed.purity_source() == "self"


def test_a_retention_model_row_cannot_assert_a_proven_vinyl_ether():
    """`P-` claims the alk-1-enyl bond is proven, and an RT-model row has no fragmentation at all.

    Six of 14 `Alkenyl-TG` rows on one study shipped `P-` named by retention time alone. The
    existing ether guard could never reach them: it matches the substrings `Plasmenyl`/`Plasmanyl`,
    which `Alkenyl-TG` does not contain. The class prefix goes with the claim, because `Alkenyl-`
    means alk-1-enyl just as surely as `P-` does.
    """
    sample = Sample(file="S1")
    group = CompoundGroup(name="g", mw=800.0, retention=10.0, max_area=1.0,
                          has_ms2=False, areas=[1.0])
    group.rtls_identification = "Alkenyl-TG P-52:3"

    name, lipid_class = group.identification()
    assert name == "TG O-52:3", "no MS2, so nothing is proven"
    assert lipid_class == "TG", "the class must not keep the claim the name dropped"

    group.rtls_identification = "Plasmanyl-PC O-30:0"
    assert group.identification()[0] == "Plasmanyl-PC O-30:0", "O- is already the weaker claim"


def test_a_retention_model_row_reports_the_ion_it_was_matched_on():
    """The adduct is determined by the match, and was being thrown away.

    `extend_identifications` matches one m/z to one composition and declines anything ambiguous,
    so a surviving hit names exactly one ion — `[M+NH4]+` and `[M+Na]+` differ by ~0.9 Da against
    a 10 ppm window worth 0.008 Da at m/z 800. The candidate tuple dropped it before it reached
    the index, leaving every RT-model row with a blank `Adduct`, which also exempted those rows
    from adduct-pair removal: the check that catches `Na - H` against `+2C +3DB` at 2.4 mDa cannot
    run on a row whose ion is unrecorded.
    """
    group = CompoundGroup(name="g", mw=800.0, retention=10.0, max_area=1.0,
                          has_ms2=False, areas=[1.0])
    group.rtls_identification = "CE 18:2"
    assert group.adduct() == "", "nothing known yet"

    group.rtls_adduct = "[M+NH4]+"
    assert group.adduct() == "[M+NH4]+", "the ion the model matched on"
    assert group.identification()[0] == "CE 18:2", "the NAME still omits the adduct"


def test_the_candidate_tuple_may_or_may_not_carry_an_adduct():
    """Older callers pass a 3-tuple; the index must not break on them."""
    from lipidloop.rtls import RetentionModel, extend_identifications

    model = RetentionModel(lipid_class="CE", prefix="", intercept=0.0, per_carbon=0.5,
                           per_double_bond=-0.3, r2=0.99, residual_sd=0.05, n_used=9, n_total=9)
    groups = []
    for _ in range(2):
        g = CompoundGroup(name="g", mw=666.6, retention=model.predict(18, 2), max_area=1.0,
                          has_ms2=False, areas=[1.0])
        g.quant_ion, g.quant_polarity = 666.6194, "+"
        groups.append(g)

    extend_identifications(groups[:1], {("CE", ""): model},
                           [("CE 18:2", 666.6194, "+")])                      # 3-tuple
    extend_identifications(groups[1:], {("CE", ""): model},
                           [("CE 18:2", 666.6194, "+", "[M+NH4]+")])          # 4-tuple

    assert groups[0].rtls_identification == "CE 18:2" and groups[0].rtls_adduct == ""
    assert groups[1].rtls_identification == "CE 18:2"
    assert groups[1].rtls_adduct == "[M+NH4]+"


def test_the_candidate_tuple_is_built_the_way_the_pipeline_builds_it():
    """Exercises the call site, not just the function it calls.

    `LibrarySpectrum.adduct` is a PROPERTY. `pipeline.py` built its candidate tuple with
    `s.adduct()`, which calls the returned string — `TypeError: 'str' object is not callable`.
    286 tests passed, because every one of them called `extend_identifications` directly with
    hand-made tuples and none went through the line that constructs them. It took down the
    negative half of a five-hour run.
    """
    from lipidloop.msp import LibrarySpectrum
    from lipidloop.peaks import sum_composition

    spectrum = LibrarySpectrum(name="PC 16:0_18:1 [M+H]+;", precursor_mz=760.585, mz=[],
                               intensity=[], annotations=[], library="lib", is_lipidex=True)

    # `polarity` is derived from the adduct, so it is a property too — three of the four fields
    # in this tuple are, which is exactly why calling one of them was easy to get wrong.
    # verbatim the expression in pipeline.py's candidate set comprehension
    candidate = (sum_composition(spectrum.lipid), round(spectrum.precursor_mz, 4),
                 "+" if spectrum.polarity == "positive" else "-", spectrum.adduct)

    assert candidate == ("PC 34:1", 760.585, "+", "[M+H]+")
    assert not callable(spectrum.adduct), "a property, so calling it would fail at runtime"
