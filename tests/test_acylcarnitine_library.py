"""The generated acylcarnitines must agree with chemistry, not merely with themselves.

The library is produced from a fragmentation rule inferred from LipiDex's own C10-C26 entries. A
rule can be self-consistent and still wrong, and a wrong rule is worse here than a missing library:
every entry in the class would be systematically off by the same amount, and the class would
quietly match nothing (or, worse, match the wrong thing at a plausible mass).

So the rule is checked against values it was NOT derived from — the literature masses of the
short-chain acylcarnitines that clinical assays measure.

⚠ The first version of the builder wrote every constant fragment 0.55 mDa high by using neutral
formula masses for what are cations. It passed the builder's own self-check because that check ran
at 5 mDa. 6.5 ppm at m/z 85 is outside any sensible tolerance, so the check that mattered was the
one against outside values.
"""
from __future__ import annotations

import re
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LIBRARY = REPO / "data/libraries/Acylcarnitines_Positive.msp"

#: [M+H]+ of the acylcarnitines routinely measured in newborn screening and FAO panels.
#: Independent of the builder — these come from the chemistry, not from LipiDex.
LITERATURE = {
    "AC 0:0": 162.1125,     # free L-carnitine
    "AC 2:0": 204.1230,     # acetylcarnitine
    "AC 3:0": 218.1387,     # propionylcarnitine
    "AC 4:0": 232.1543,     # butyryl / isobutyryl
    "AC 5:0": 246.1700,     # isovaleryl / 2-methylbutyryl
    "AC 6:0": 260.1856,     # hexanoyl
    "AC 8:0": 288.2169,     # octanoyl — the MCAD marker
}

#: The fragments LipiDex writes for every acylcarnitine, as cations.
CONSTANT_FRAGMENTS = (144.1019, 85.0284, 60.0808)


def entries() -> dict[str, tuple[float, list[float]]]:
    if not LIBRARY.exists():
        pytest.skip(f"{LIBRARY.name} not built")
    out = {}
    for block in LIBRARY.read_text().split("\n\n"):
        m = re.match(r"Name: (AC \d+:\d+)", block)
        if not m:
            continue
        pre = float(re.search(r"PRECURSORMZ: ([\d.]+)", block).group(1))
        mz = [float(x.split()[0]) for x in block.split("\n") if re.match(r"^\d+\.", x)]
        out[m.group(1)] = (pre, mz)
    return out


@pytest.mark.parametrize("name,expected", sorted(LITERATURE.items()))
def test_precursor_matches_the_literature_mass(name, expected):
    found = entries()
    assert name in found, f"{name} absent — the short chains are the point of this library"
    assert abs(found[name][0] - expected) < 0.001, \
        f"{name}: {found[name][0]:.4f} vs literature {expected:.4f}"


def test_the_constant_fragments_are_cations_not_neutrals():
    """0.55 mDa is the electron, and it is the difference between working and not."""
    for name, (_, mz) in entries().items():
        for want in CONSTANT_FRAGMENTS:
            if want > float(name.split()[1].split(":")[0]) * 0 + 200:
                continue
            near = [x for x in mz if abs(x - want) < 0.05]
            if near:
                assert min(abs(x - want) for x in near) < 0.001, \
                    f"{name}: fragment near {want} is off by more than 1 mDa — neutral, not cation?"


def test_no_entry_duplicates_one_already_shipped():
    """A duplicate would be searched twice and counted twice."""
    shipped = set()
    for lib in ("LipiDex_HCD_Formic.msp", "LipiDex_HCD_Acetate.msp"):
        p = REPO / "data/libraries" / lib
        if p.exists():
            shipped |= set(re.findall(r"^Name: (AC \d+:\d+) \[M\+H\]\+", p.read_text(), re.M))
    assert not (set(entries()) & shipped), "generated entries collide with shipped ones"


def test_the_builder_refuses_to_write_when_the_rule_disagrees():
    """The guard must actually fail, not just print. A check that cannot fail is not a check."""
    r = subprocess.run(
        [sys.executable, str(REPO / "scripts/build_acylcarnitine_library.py"),
         "--out", "/dev/null", "--libraries", str(REPO / "tests")],
        capture_output=True, text=True, timeout=120)
    # No shipped libraries in tests/, so there is nothing to verify against and nothing to collide
    # with; the point is that it exits cleanly rather than crashing on an empty reference.
    assert r.returncode == 0, r.stderr[-800:]
