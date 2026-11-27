"""Demoting a plasmalogen name (`P-`) to the unresolved ether form (`O-`) must keep the formula.

`P-` does not count the vinyl double bond and `O-` does, so `P-x:y` is `O-x:(y+1)`. The demotion
used to swap the letter only, and every demoted row named a formula 2 H heavier than its ion
(Set 2: the reference `Plasmanyl-PC O-32:1` at m/z 718.58 was delivered as `PC O-32:0`).
"""
import re
import sys
from collections import Counter
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.lsi import shorthand                                              # noqa: E402
from lipidloop.peaks import CompoundGroup, _demote_proven_ether, sum_composition  # noqa: E402
from lipidloop.polarity import _ether_key                                        # noqa: E402

LIBRARY = Path(__file__).resolve().parents[1] / "data" / "libraries" / "LipiDex_HCD_Formic.msp"


def _group(sum_id):
    g = CompoundGroup(name="", mw=0.0, retention=0.0, max_area=0.0, has_ms2=False, areas=[])
    g.sum_id = sum_id
    return g


def test_demotion_adds_the_vinyl_double_bond():
    assert _group("Plasmenyl-PC P-32:0").top_ether() == "PC O-32:1"
    assert _group("Plasmenyl-PE P-16:0_20:4").top_ether() == "PE O-16:1_20:4"
    assert _demote_proven_ether("Alkenyl-TG P-52:3") == "TG O-52:4"
    assert _demote_proven_ether("Alkenyl-DG P-34:1") == "DG O-34:2"
    assert _demote_proven_ether("Plasmanyl-PC O-30:0") == "Plasmanyl-PC O-30:0"


def test_retention_model_shorthand_adds_the_vinyl_double_bond():
    assert shorthand("Plasmenyl-PC P-34:1", "RT model") == "PC O-34:2"
    assert shorthand("Plasmenyl-PE P-16:0_20:4", "RT model") == "PE O-16:1_20:4"
    assert shorthand("Plasmenyl-PC P-34:1", "MS2") == "PC P-34:1", "with MS2 the claim stands"
    assert shorthand("PC O-34:2", "RT model") == "PC O-34:2", "already demoted: no second +1"


def test_a_demoted_row_pairs_with_its_linkage_named_partner():
    """Same formula, other polarity: `_ether_key` must see the demoted `PC O-` row."""
    assert _ether_key("PC O-34:1") == _ether_key("Plasmenyl-PC P-34:0") == _ether_key("Plasmanyl-PC O-34:1")
    assert _ether_key("PC 34:1") is None, "a diacyl PC is not an ether"


def _formula(comment):
    m = re.search(r"Formula=(\S+)", comment)
    return Counter({el: int(n or 1) for el, n in re.findall(r"([A-Z][a-z]?)(\d*)", m.group(1))}) if m else None


@pytest.mark.skipif(not LIBRARY.exists(), reason="library not present")
def test_demoted_names_match_the_library_formula_of_the_alkyl_entry():
    """Against the library itself: every plasmenyl entry, demoted, has the formula of the plasmanyl
    entry of the demoted sum composition wherever the library carries one."""
    entries, name = {}, None
    for line in LIBRARY.open(errors="replace"):
        if line.startswith("Name:"):
            name = line[5:].strip().rstrip(";")
        elif line.startswith("Comment:") and name:
            entries.setdefault(name, _formula(line))
    alkyl = {}
    for n, f in entries.items():
        m = re.match(r"^Plasmanyl-(\S+) (O-\S+) (\[.*)$", n)
        if m and f:
            alkyl[(m.group(1), sum_composition(f"X {m.group(2)}").split(" ", 1)[1], m.group(3))] = f
    checked = 0
    for n, f in entries.items():
        m = re.match(r"^(Plasmenyl-\S+ P-\S+) (\[.*)$", n)
        if not m or not f:
            continue
        demoted = _group(sum_composition(m.group(1))).top_ether()
        head, body = demoted.split(" ", 1)
        partner = alkyl.get((head, body, m.group(2)))
        if partner is None:
            continue
        assert partner == f, f"{n}: demoted to {demoted}, formula differs from the plasmanyl entry"
        checked += 1
    assert checked >= 20, f"only {checked} pairs checked"
