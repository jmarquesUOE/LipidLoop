"""The decoy construction must be wrong in the way it intends, and in no other way.

Two defects motivate this file, both found by measurement rather than by reading:

1. **The chain pattern could not read half the annotations.** It required a digit straight after
   `[`, so `[d18:1]`, `[t16:0]`, `[O-16:0]`, `[P-16:0]` and `[14:0, 16:1]` all failed to match.
   Every sphingoid, ether and cardiolipin entry was therefore skipped by the CH2 construction, and
   — because `head_group_decoy` guards on the same pattern — was handed a head-group decoy the
   module's own docstring forbids.

2. **A decoy can be a correct answer.** `DECOY_ PE-NMe2 15:0_17:1 [M+H]+` built at `--shift 1`
   has the precursor, the formula and the acylium ions of a real `PE 16:0_18:1`. It was 40.3% of
   every decoy hit in the v5 batch — right answers counted as false positives, and outranking and
   deleting the true PE as they went.

⚠ The tests below assert the SUBSTITUTED TEXT, not just that something moved. The failure that
costs most here is silent: dropping the `d` from a shifted `[d18:1]` relabels a sphingoid base as
an acyl chain, which is a different molecule with a different fragment formula, and the file still
parses and still searches.
"""
from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("make_decoy_library",
                                               REPO / "scripts/make_decoy_library.py")
mdl = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(mdl)


# --------------------------------------------------------------------------------------------
# 1. the chain pattern
# --------------------------------------------------------------------------------------------

#: Every annotation dialect that appears in the shipped libraries, with the chain text the
#: pattern must capture. Taken from `grep -oP '_\[[^]]*\]'` over LipiDex_HCD_* and LipidBlast_*.
ANNOTATIONS = [
    ('155.143 100 "O-1_Alkyl Fragment_[10:0]"', "10:0"),
    ('337.3481 999 "-118.0499_Sphingoid Neutral Loss_[d18:0]"', "d18:0"),
    ('256.2282 50 "C2H3O1N1_Phytosphingoid Fragment_[t16:0]"', "t16:0"),
    ('383.3156 999 "N-1H-4_Plasmanyl Neutral Loss_[O-16:0]"', "O-16:0"),
    ('167.1076 999 "-268.0952_Plasmenyl Neutral Loss_[P-16:0]"', "P-16:0"),
    ('591.4031 200 "C3H6O4P1_Cardiolipin DG Fragment_[14:0, 14:0]"', "14:0, 14:0"),
]


@pytest.mark.parametrize("line,chains", ANNOTATIONS)
def test_chain_pattern_reads_every_dialect(line, chains):
    m = mdl.CHAIN.search(line)
    assert m is not None, "the pattern that missed these skipped 14,061 LipidBlast entries"
    assert m.group(1) == chains


def test_unannotated_peak_is_not_a_chain():
    """`[]` must match nothing — it is the head-group/class peak the decoy deliberately keeps."""
    assert mdl.CHAIN.search('144.1019 100 "C7H14O2N1_Fragment_[]"') is None


@pytest.mark.parametrize("chains,shift,expected,n", [
    ("10:0", 3, "13:0", 1),
    ("d18:1", 3, "d21:1", 1),          # ⚠ the `d` must survive: d21:1 is a sphingoid, 21:1 is not
    ("t16:0", 3, "t19:0", 1),
    ("O-16:0", 3, "O-19:0", 1),
    ("P-16:0", 3, "P-19:0", 1),
    ("14:0, 16:1", 3, "17:0, 19:1", 2),
])
def test_shift_chains_carries_the_prefix(chains, shift, expected, n):
    assert mdl.shift_chains(chains, shift) == (expected, n)


def test_sphingoid_entry_gets_a_chain_decoy_and_keeps_its_prefix():
    entry = ["Name: Cer[NP] t16:0_10:0 [M-H]-;", "PRECURSORMZ: 442.3896", "Num Peaks: 1",
             '256.2282 50 "C2H3O1N1_Phytosphingoid Fragment_[t16:0]"']
    out = mdl.shift_entry(entry, 3)
    assert out is not None, "sphingoid entries were skipped entirely by the old pattern"
    assert out[-1] == '298.2751 50 "C2H3O1N1_Phytosphingoid Fragment_[t19:0]"'


def test_two_chain_fragment_moves_by_one_ch2_per_chain():
    """⚠ Both chains move, so the mass must move by 2 * shift * CH2, not by shift * CH2.

    Measured on the cardiolipin DG fragments, the only two-chain annotation shipped:
    `[14:0, 14:0]` 591.4031 vs `[14:0, 16:1]` 617.4188 is 2*CH2 - H2, so each chain's carbons
    count independently. A single-CH2 move would write an m/z its own annotation contradicts.
    """
    entry = ["Name: CL 14:0_14:0_14:0_14:0 [M-H]-;", "PRECURSORMZ: 1239.8392", "Num Peaks: 1",
             '591.4031 200 "C3H6O4P1_Cardiolipin DG Fragment_[14:0, 14:0]"']
    out = mdl.shift_entry(entry, 3)
    mz = float(out[-1].split()[0])
    assert out[-1].endswith('_[17:0, 17:0]"')
    assert mz == pytest.approx(591.4031 + 6 * mdl.CH2, abs=1e-4)


def test_head_group_guard_now_fires_for_sphingoids():
    """The guard reads the same pattern, so fixing one fixed the other.

    ⚠ This is the half of defect 1 that produced a wrong library rather than a missing one:
    4,200 of `DECOY2_HeadGroup_LipiDex_HCD_Formic`'s 14,205 entries DID carry chain fragments and
    should never have been given a head-group decoy.
    """
    entry = ["Name: Cer[NP] t16:0_10:0 [M-H]-;", "PRECURSORMZ: 442.3896", "Num Peaks: 1",
             '256.2282 50 "C2H3O1N1_Phytosphingoid Fragment_[t16:0]"']
    assert mdl.head_group_decoy(entry) is None


def test_head_group_decoy_still_built_for_a_truly_chainless_entry():
    entry = ["Name: SM d18:1_16:0 [M+H]+;", "PRECURSORMZ: 703.5754", "Num Peaks: 1",
             '184.0733 999 "C5H15N1O4P1_Fragment_[]"']
    out = mdl.head_group_decoy(entry)
    assert out is not None
    assert float(out[-1].split()[0]) == pytest.approx(184.0733 + mdl.HEAD_GROUP_SHIFT, abs=1e-4)


# --------------------------------------------------------------------------------------------
# 2. the collision screen
# --------------------------------------------------------------------------------------------

#: The real `PE 16:0_18:1 [M+H]+` as LipidBlast writes it: precursor 718.5387 and the two acylium
#: ions. LipiDex writes the same lipid with glycerol-bearing fragments at 313.2737/339.2894
#: instead, which is why the screen has to look across libraries and not only at the decoy's own.
REAL_PE = ["Name: PE 16:0_18:1 [M+H]+;", "PRECURSORMZ: 718.5387", "Num Peaks: 3",
           '577.519 999 "-141.0191_Neutral Loss_[]"',
           '265.2524 50 "-198.0533_Alkyl Neutral Loss_[16:0]"',
           '239.2368 50 "-198.0533_Alkyl Neutral Loss_[18:1]"']

#: The source of the flagship defect. At +1 CH2 on each of two chains its fragments become the
#: real PE's acylia, because PE-NMe2 is PE plus exactly 2 CH2 in the head group.
PE_NME2 = ["Name: PE-NMe2 15:0_17:1 [M+H]+;", "PRECURSORMZ: 718.5387", "Num Peaks: 7",
           '72.0808 200 "C4H10N1_Fragment_[]"',
           '225.2213 100 "O-1_Alkyl Fragment_[15:0]"',
           '207.2107 50 "H-2O-2_Alkyl Fragment_[15:0]"',
           '549.4877 999 "C-4H-12O-4N-1P-1_Neutral Loss_[]"',
           '170.0577 500 "C4H13N1O4P1_Fragment_[]"',
           '251.2369 100 "O-1_Alkyl Fragment_[17:1]"',
           '233.2264 50 "H-2O-2_Alkyl Fragment_[17:1]"']


def _index(*entries):
    idx = mdl.TargetIndex()
    for e in entries:
        name, precursor, chain_mz, _ = mdl.parse_entry(e)
        idx.bins.setdefault(int(precursor * 100), []).append(
            (name, precursor, chain_mz, mdl.name_composition(name)))
        idx.n += 1
    return idx


def test_shift_1_pe_nme2_is_caught():
    """The whole reason this screen exists. Shift 1 reproduces a real PE; it must not survive."""
    decoy = mdl.shift_entry(PE_NME2, 1)
    _n, precursor, chain_mz, _a = mdl.parse_entry(decoy)
    hit = _index(REAL_PE).collision(precursor, chain_mz, mdl.fragment_composition(decoy))
    assert hit is not None
    assert hit[0] == "PE 16:0_18:1 [M+H]+"


def test_shift_3_pe_nme2_survives():
    """Shift 3 needs a 6 CH2 head-group gap for a two-chain lipid, and no common pair has one."""
    decoy = mdl.shift_entry(PE_NME2, 3)
    _n, precursor, chain_mz, _a = mdl.parse_entry(decoy)
    assert _index(REAL_PE).collision(precursor, chain_mz,
                                     mdl.fragment_composition(decoy)) is None


def test_containment_alone_would_have_passed_the_flagship_defect():
    """⚠ The obvious form of the screen does not work, and this records why.

    "reject when the decoy's chain fragments are a subset of a real target's" reads as the natural
    statement, and it lets `DECOY_ PE-NMe2 15:0_17:1` straight through: the decoy carries four
    chain peaks (each acylium and its water loss) where LipidBlast's real PE lists only the two
    acylia. A library entry is an editor's list, not a spectrum. Hence both directions, plus the
    composition rule.
    """
    decoy = mdl.shift_entry(PE_NME2, 1)
    _n, _p, decoy_frags, _a = mdl.parse_entry(decoy)
    _n2, _p2, real_frags, _a2 = mdl.parse_entry(REAL_PE)
    assert not mdl._contained(decoy_frags, real_frags)     # the rule that would have failed
    assert mdl._contained(real_frags, decoy_frags)         # the rule that catches it


def test_precursor_outside_ms1_tolerance_is_not_a_collision():
    far = list(REAL_PE)
    far[1] = "PRECURSORMZ: 718.6000"        # 62 mDa away: a different bin, a different lipid
    decoy = mdl.shift_entry(PE_NME2, 1)
    _n, precursor, chain_mz, _a = mdl.parse_entry(decoy)
    assert _index(far).collision(precursor, chain_mz, mdl.fragment_composition(decoy)) is None


def test_empty_fragment_set_is_never_a_collision():
    """⚠ Every set contains the empty set. Without this guard the screen would reject every
    head-group decoy on a vacuous subset rather than on evidence."""
    assert _index(REAL_PE).collision(718.5387, [], frozenset()) is None


def test_name_composition_ignores_adduct_and_decoy_prefix():
    assert mdl.name_composition("DECOY_ PE 16:0_18:1 [M+H]+") == frozenset({"16:0", "18:1"})
    assert mdl.name_composition("SM d18:1_16:0 [M+H]+") == frozenset({"d18:1", "16:0"})
    assert mdl.name_composition("Cholesterol [M+H]+") == frozenset()
