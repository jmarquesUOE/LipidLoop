"""The prenol-quinone entries must agree with the spectra they were built from.

Checked against values the builder was NOT derived from — the m/z the targeted extraction
searched for (`CoQ_check/cq01_xic_pos.py`, `cq03_ubiquinol_check.py`) and the fragment the MS2
screen scored — which is the discipline `test_acylcarnitine_library.py` adopted after a build
that passed its own self-check while every constant fragment sat 0.55 mDa high.

The last test here is the one that matters most, and it is not about this library's arithmetic:
ubiquinol and the ubiquinone M+2 isotope are 8.9 mDa apart, which is INSIDE `search.MS1_TOL`.
If that ever stops being true — a tighter tolerance, a different adduct — the retention
constraint could be relaxed. Until then it cannot be, and this test says so out loud.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LIBRARY = REPO / "data/libraries/Ubiquinones_Positive.msp"
ANCHORS = REPO / "data/libraries/Ubiquinones_RT_anchors.csv"

sys.path.insert(0, str(REPO / "src"))
from lipidloop.search import MS1_TOL              # noqa: E402

#: Precursor m/z as cq01/cq03 computed and extracted them, independent of this builder.
OBSERVED_PRECURSOR = {
    "CoQ10 [M+NH4]+": 880.7177,
    "CoQ10 [M+H]+": 863.6912,
    "CoQ10H2 [M+NH4]+": 882.7334,
    "CoQ9 [M+NH4]+": 812.6551,
}
#: The quinone head group, C10H13O4+ — base peak of every qualifying spectrum in both studies.
HEAD_GROUP = 197.0808
#: The ubiquinone [M+NH4]+ M+2 isotope: 880.7177 + 2 x (13C - 12C).
UBIQUINONE_M2 = 882.7244
#: Measured retention, an in-house Q Exactive Plus study.
RETENTION = {"CoQ10": 21.07, "CoQ10H2": 19.50, "CoQ9": 19.10}

TOL = 0.002


def entries() -> dict[str, tuple[float, list[float], list[int]]]:
    if not LIBRARY.exists():
        pytest.skip(f"{LIBRARY.name} not built — run scripts/build_ubiquinone_library.py")
    out = {}
    for block in LIBRARY.read_text().split("\n\n"):
        m = re.match(r"Name: (CoQ\S+ \[[^\]]+\]\+);", block)
        if not m:
            continue
        pre = float(re.search(r"PRECURSORMZ: ([\d.]+)", block).group(1))
        peaks = [x for x in block.split("\n") if re.match(r"^\d+\.", x)]
        out[m.group(1)] = (pre, [float(x.split()[0]) for x in peaks],
                           [int(x.split()[1]) for x in peaks])
    return out


@pytest.mark.parametrize("name", sorted(OBSERVED_PRECURSOR))
def test_precursor_matches_the_extracted_value(name):
    found = entries()
    assert name in found, f"{name} absent — it is one of the four confirmed ions"
    assert abs(found[name][0] - OBSERVED_PRECURSOR[name]) < TOL


@pytest.mark.parametrize("name", sorted(OBSERVED_PRECURSOR))
def test_quinone_head_group_is_the_base_peak(name):
    """197.0808 was the base peak of every spectrum pooled, so it must be 999 in every entry.

    It is also the only fragment that says "quinone" at all: the rest of the pattern is the
    polyprenyl hydrocarbon ladder, which a sterol gives too.
    """
    pre, mzs, intensities = entries()[name]
    hits = [i for mz, i in zip(mzs, intensities) if abs(mz - HEAD_GROUP) < TOL]
    assert hits == [999], f"{name}: head group {HEAD_GROUP} not present at 999 ({hits})"


def test_no_phosphocholine_fragment():
    """⚠ 184.0733 was MEASURED at median 395 in 7/7 CoQ9 spectra and must not be in the library.

    It is the phosphocholine head group of a co-isolated PC or SM, not a quinone fragment.
    Writing it would make these entries score against every PC spectrum in a run — the failure
    that withdrew `Cer[EOS]`, arriving by a different route.
    """
    for name, (_, mzs, _) in entries().items():
        assert not any(abs(mz - 184.0733) < 0.01 for mz in mzs), f"{name} carries 184.0733"


def test_ubiquinol_cannot_be_separated_from_the_ubiquinone_isotope_by_mass():
    """The whole reason this library ships a retention table and a checker.

    If this test fails because the gap grew past MS1_TOL, the retention constraint is no longer
    load-bearing and the developer notes (PRENOL_QUINONES, not distributed) should be rewritten. It will not fail for that
    reason: the gap is a property of 13C, not of the method.
    """
    quinol = OBSERVED_PRECURSOR["CoQ10H2 [M+NH4]+"]
    gap = abs(quinol - UBIQUINONE_M2)
    assert gap < MS1_TOL, (
        f"{gap * 1000:.1f} mDa now exceeds MS1_TOL {MS1_TOL}; the isotope no longer reaches the "
        f"ubiquinol entry and the retention rule can be reconsidered")
    assert gap / quinol * 1e6 < 15, "the two are within 15 ppm — not resolvable at this peak width"


def test_retention_anchors_cover_every_entry_and_put_the_quinol_first():
    """Every ion in the library needs a measured retention, and the quinol must elute earlier.

    `check_ubiquinol_rt.py` reads this table; a missing row there is a row it cannot judge.
    """
    if not ANCHORS.exists():
        pytest.skip("anchors not built")
    anchors = {(r["name"], r["adduct"]): float(r["retention"])
               for r in csv.DictReader(ANCHORS.open())}
    for name in entries():
        cls, adduct = name.split(" ", 1)
        assert (cls, adduct) in anchors, f"{name} has no measured retention in {ANCHORS.name}"
        assert abs(anchors[(cls, adduct)] - RETENTION[cls]) < 0.3
    assert anchors[("CoQ10H2", "[M+NH4]+")] < anchors[("CoQ10", "[M+NH4]+")] - 0.75, (
        "the reduced form must elute clearly before the oxidised one, or the only test that "
        "separates it from the M+2 isotope has no margin")
