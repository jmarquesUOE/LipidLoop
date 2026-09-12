"""The sn-position column: that it is computed where the data exists, and never overclaims."""
import csv
import inspect
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop import sn_position                                        # noqa: E402
from lipidloop.mgf import SampleSpectrum                                 # noqa: E402
from lipidloop.msp import LibrarySpectrum                                # noqa: E402
from lipidloop.search import RESULT_COLUMNS, sn_evidence_cell            # noqa: E402


def _pc() -> LibrarySpectrum:
    """A PC formate entry annotated the way LipiDex writes carboxylate fragments."""
    return LibrarySpectrum(
        name="PC 16:0_18:1 [M+FA-H]-", precursor_mz=804.5755,
        mz=[255.2330, 281.2486, 168.0426], intensity=[999.0, 999.0, 500.0],
        annotations=['FA_Fatty Acid_[16:0]', 'FA_Fatty Acid_[18:1]', 'HG_Head Group_[]'],
        library="LipiDex_HCD_Formic", is_lipidex=True, optimal_polarity=True,
        types=["FA_Fatty Acid", "FA_Fatty Acid", "HG_Head Group"],
        fatty_acids=["16:0", "18:1", ""])


def _spectrum(intensities) -> SampleSpectrum:
    return SampleSpectrum(number=1, precursor=804.5755, retention=10.0, polarity="-",
                          mz=[255.2331, 281.2487, 168.0427], intensity=list(intensities))


# ──────────────────────────────────────────────────────────────────────────────────────────
# The wiring. This is what the first attempt got wrong, and it failed in a way that looked
# like a null result rather than a bug.
# ──────────────────────────────────────────────────────────────────────────────────────────

def test_the_production_row_builder_writes_the_column():
    """⚠ THE REGRESSION THIS FILE EXISTS FOR.

    `pipeline.py` does not call `search.iter_results`; it builds its own row dict inline, and
    `iter_results` serves only the tests and the validation harness. The first attempt wired the
    column into `iter_results` alone, so on real data it never executed and came back empty on all
    430 matches — including declined cases, which should have carried a reason. That was read as
    "nothing qualified" and diagnosed as the wrong hook point, when the truth was that the hook was
    in a function the pipeline never calls.

    Asserting on the source is deliberate: the alternative is a full pipeline run, and the failure
    is structural — one row builder gaining a column the other lacks — so it is visible statically.
    """
    from lipidloop import pipeline

    source = inspect.getsource(pipeline)
    assert 'rows[-1]["sn Evidence"] = sn_evidence_cell(' in source, \
        "pipeline.py must fill the column itself — iter_results is not on the production path"

    # And it must be filled AFTER `Potential Fragments`, whose de-duplication prefix is a join
    # over this dict's values. Ahead of it, a new column changes which masses are suppressed and
    # so changes the in-source-fragment filter downstream.
    assert (source.index('rows[-1]["Potential Fragments"]')
            < source.index('rows[-1]["sn Evidence"]')), \
        "sn Evidence must not land in the Potential Fragments prefix"


def test_both_row_builders_agree_on_the_column_set():
    """One helper, two callers: the column must exist in the declared schema."""
    assert "sn Evidence" in RESULT_COLUMNS
    # Appended, not inserted — `Potential Fragments` must keep its position for the same reason.
    assert RESULT_COLUMNS[-1] == "sn Evidence"
    assert RESULT_COLUMNS.index("Potential Fragments") == len(RESULT_COLUMNS) - 2


def test_the_search_table_column_is_never_blank():
    """A decline carries its reason, which is how "never ran" stays distinguishable from "no call".

    An empty cell is what the broken version produced, and it was indistinguishable from a genuine
    decline. Anything non-empty here fails loudly instead.
    """
    lib, ms2 = _pc(), _spectrum([1000.0, 3200.0, 700.0])
    assert sn_evidence_cell(lib, ms2).strip()

    # A class the rule does not cover still gets a reason, not a blank. Negative adduct on
    # purpose: with a positive one the polarity guard fires first and this would pass without
    # ever reaching the class test.
    cer = _pc()
    cer.name = "Cer[NDS] d18:0_24:1 [M+FA-H]-"
    assert "class not covered" in sn_evidence_cell(cer, ms2)


def test_the_column_round_trips_through_the_search_table():
    """The written form and its reader are a matched pair, and must change together.

    `parse_cell` is what any analysis of the search tables uses to get the call back out — and it is
    what the withdrawn aggregation used. Nothing recomputes the ratio downstream: the sample
    intensities and the library's per-fragment chain annotations do not survive into a
    `_Results.csv` row, so the written cell is the only record.
    """
    cell = sn_evidence_cell(_pc(), _spectrum([1000.0, 3200.0, 700.0]))
    chain, ratio = sn_position.parse_cell(cell)
    assert chain == "18:1"
    assert ratio == 3.2

    # A decline parses as "no call" rather than raising — it is the common case, not an error.
    assert sn_position.parse_cell("class not covered — the sn-2 rule is not established here") \
        == ("", 0.0)
    assert sn_position.parse_cell("") == ("", 0.0)


def test_a_search_table_written_before_this_column_existed_still_reads(tmp_path):
    """Backward compatibility: an old `_Results.csv` has no `sn Evidence` and must still load."""
    old = [c for c in RESULT_COLUMNS if c != "sn Evidence"]
    path = tmp_path / "old_Results.csv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=old)
        w.writeheader()
        w.writerow({c: "" for c in old})

    with path.open(newline="") as fh:
        row = next(csv.DictReader(fh))
    assert row.get("sn Evidence", "") == ""


# ──────────────────────────────────────────────────────────────────────────────────────────
# ⚠ The delivered table must NOT carry this column.
# ──────────────────────────────────────────────────────────────────────────────────────────

def test_the_delivered_table_does_not_carry_the_column(tmp_path):
    """★ Withdrawn from `Final_Results.csv` on purpose — do not let it back without the fix.

    The column was wired end to end and then pulled, because its premise failed validation: on
    delivery-grade negative spectra the call is set by which fatty acids are present rather than by
    how they are arranged. `18:2` loses 14 of 16 pairs, `22:4` and `20:3` lose every one, and
    `18:0` vs `20:4` — the textbook arachidonate PE, and the case the method should be surest
    about — comes out 7-7, an exact coin flip. `22:6` and `22:5` winning nearly everything looks
    canonical but is not evidence: they would win at sn-1 too.

    So a delivered `sn Evidence` value would be read as a position assignment and would not be one.
    The per-spectrum diagnostic in `search/` stays, because it is what made this visible.

    Restoring it needs a measurement that cancels the chain term — a per-fatty-acid response
    baseline from species of known arrangement — not a threshold on the ratio. The aggregation code
    removed with it is in commit 5fe6dea if it is wanted back.
    """
    import csv as _csv
    from lipidloop.peakfinder import META_COLUMNS, PeakFinderResult, write_results
    from lipidloop.peaks import CompoundGroup, Lipid, LipidCandidate, Sample

    sample = Sample(file="sample_one")
    lipid = Lipid(retention=5.0, precursor=804.5755, sample=sample, dot=900.0, rev_dot=980.0,
                  lipid_string="PE 18:0_20:4 [M+FA-H]-", lib_precursor=804.5755, purity=99,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    candidate = LipidCandidate(lipid)
    group = CompoundGroup(name="g", mw=804.5755, retention=5.0, max_area=1.0,
                          has_ms2=True, areas=[1.0])
    group.lipid_candidates = [candidate]
    group.purity = 99.0
    group.sum_id = "PE 38:4"
    group.final_lipid_id = candidate

    out = tmp_path / "Final_Results.csv"
    write_results(PeakFinderResult([group], [sample]), out, scores=True, ms2_support=True)
    header = next(_csv.reader(out.open(newline="")))

    assert "sn Evidence" not in header, "the delivered table must not report sn position"
    assert "Chain Evidence" in header, "the columns beside it must be unaffected"
    assert "sn Evidence" not in META_COLUMNS, "unregister it too, or qc.load_run keeps a ghost"


# ──────────────────────────────────────────────────────────────────────────────────────────
# Polarity. Found by running the wired pipeline on real data, not by reading the module.
# ──────────────────────────────────────────────────────────────────────────────────────────

def test_positive_mode_never_gets_an_sn_call():
    """⚠ The sn-2 rule is negative-mode carboxylate chemistry and must not cross polarities.

    In negative mode both chains leave as carboxylate anions and the sn-2 ester, being more
    labile, gives the more intense one. A protonated phospholipid fragments to the head group
    instead, and its acyl-related ions are neutral losses whose ordering reflects which loss is
    favoured — not ester lability. The module guarded TG for exactly this reason and then let a
    whole polarity through.

    The first wired run made three calls on `PE 18:0_20:4 [M+H]+` before this guard existed.
    """
    lib, ms2 = _pc(), _spectrum([1000.0, 3200.0, 700.0])
    assert sn_position.evidence("PC", "PC 16:0_18:1", lib.mz, lib.fatty_acids,
                                ms2.mz, ms2.intensity, polarity="+").favoured_sn2 == ""
    assert "positive mode" in str(sn_position.evidence(
        "PC", "PC 16:0_18:1", lib.mz, lib.fatty_acids, ms2.mz, ms2.intensity, polarity="positive"))

    # Negative still calls, so the guard is not simply switching the feature off.
    assert sn_position.evidence("PC", "PC 16:0_18:1", lib.mz, lib.fatty_acids,
                                ms2.mz, ms2.intensity, polarity="-").favoured_sn2 == "18:1"


def test_either_polarity_being_positive_declines_the_call():
    """Matching is on precursor mass alone, so the two can disagree — and disagreement declines.

    `search.py` never checks polarity when matching, so a positive spectrum can win a negative
    library entry on a mass coincidence. The sample's polarity is the ion that fragmented; the
    library adduct is the identity claimed. The rule needs both.
    """
    neg_lib = _pc()                                    # "[M+FA-H]-"  -> negative
    pos_lib = _pc()
    pos_lib.name = "PC 16:0_18:1 [M+H]+"

    good = _spectrum([1000.0, 3200.0, 700.0])
    good.polarity = "-"
    assert "favoured sn-2" in sn_evidence_cell(neg_lib, good)

    # positive library entry, negative spectrum
    assert "positive mode" in sn_evidence_cell(pos_lib, good)

    # negative library entry, positive spectrum
    pos_spec = _spectrum([1000.0, 3200.0, 700.0])
    pos_spec.polarity = "+"
    assert "positive mode" in sn_evidence_cell(neg_lib, pos_spec)
