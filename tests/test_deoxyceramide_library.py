"""The 1-deoxyceramide entries must agree with the spectra they were built from.

Every other library we search is combinatorially generated and validated against itself. This one
exists because four species were confirmed by MS2 in two studies on two instrument vendors
(the developer notes (DEOXYSPHINGOLIPIDS, not distributed)), so it can be checked against something better than its own
arithmetic: the observed fragment masses.

⚠ The first build put every N-acyl amide 1.0078 low, by giving the neutral amide 2n hydrogens
instead of 2n+1. That is close enough to a measured value to read as right in a table and far
enough that not one acyl ion would ever have matched at 0.01 Da — the class would have been
identified on its base ion alone, losing exactly the evidence that names the chain. Hence the
observed values below, which come from the spectra and not from the builder.
"""
from __future__ import annotations

import re
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
LIBRARY = REPO / "data/libraries/DeoxyCeramides_Positive.msp"

#: Observed in an in-house Q Exactive Plus study and MTBLS5163, independent of the builder. Precursor, then the four
#: ions: base - H2O, intact base, N-acyl amide, precursor - H2O.
OBSERVED = {
    "DeoxyCer m18:0_22:0": (608.6340, {"base-H2O": 268.2999, "base": 286.3104,
                                       "acylamide": 340.3574, "prec-H2O": 590.6234}),
    "DeoxyCer m18:0_24:0": (636.6653, {"base-H2O": 268.2999, "base": 286.3104,
                                       "acylamide": 368.3887, "prec-H2O": 618.6547}),
    "DeoxyCer m20:0_24:0": (664.6966, {"base-H2O": 296.3312, "base": 314.3417,
                                       "acylamide": 368.3887, "prec-H2O": 646.6860}),
}

TOL = 0.002


def entries() -> dict[str, tuple[float, list[float], list[str]]]:
    if not LIBRARY.exists():
        pytest.skip(f"{LIBRARY.name} not built — data/ is gitignored")
    out = {}
    for block in LIBRARY.read_text().split("\n\n"):
        m = re.match(r"Name: (DeoxyCer m\d+:\d+_\d+:\d+)", block)
        if not m:
            continue
        pre = float(re.search(r"PRECURSORMZ: ([\d.]+)", block).group(1))
        peaks = [x for x in block.split("\n") if re.match(r"^\d+\.", x)]
        out[m.group(1)] = (pre, [float(x.split()[0]) for x in peaks],
                           [x.split('"')[1] for x in peaks])
    return out


@pytest.mark.parametrize("name", sorted(OBSERVED))
def test_precursor_matches_the_measured_value(name):
    found = entries()
    assert name in found, f"{name} absent — it is one of the confirmed species"
    assert abs(found[name][0] - OBSERVED[name][0]) < TOL


@pytest.mark.parametrize("name", sorted(OBSERVED))
def test_every_measured_ion_is_present(name):
    _, mz, _ = entries()[name]
    for label, want in OBSERVED[name][1].items():
        assert any(abs(m - want) < TOL for m in mz), \
            f"{name}: no peak within {TOL} of the observed {label} at {want:.4f}"


def test_the_acyl_ion_is_the_amide_not_one_hydrogen_short():
    """The regression that motivated this file: 2n hydrogens instead of 2n+1."""
    _, mz, _ = entries()["DeoxyCer m18:0_24:0"]
    assert not any(abs(m - 367.3809) < TOL for m in mz), \
        "the 24:0 acyl ion is 1.0078 low — neutral amide built with 2n hydrogens"


def test_exactly_one_peak_carries_the_sphingoid_chain():
    """Two would make purity.py count two fragments for one declared chain and return 0.

    That is not hypothetical — it is why the phytoceramide entries score Purity 0 beside dot
    products of 999 (the developer notes (CHAIN_RESOLUTION, not distributed)). A library whose every row reports a
    chain-resolved name at Purity 0 is asserting chains it has no evidence for.
    """
    for name, (_, _, types) in entries().items():
        sphingoid = [t for t in types if "Sphingoid Fragment" in t]
        assert len(sphingoid) == 1, \
            f"{name}: {len(sphingoid)} Sphingoid-typed peaks; purity would come back 0"


def test_no_entry_is_isobaric_with_a_canonical_ceramide_water_loss():
    """Saturated only, and this is the reason.

    A 1-deoxyceramide with k double bonds has the formula of the canonical ceramide with k-1
    after that ceramide loses water — the same ion, 0.0 ppm. Such an entry cannot win on
    evidence, only take spectra from a correct assignment. Every name here must be saturated.
    """
    for name in entries():
        chains = name.split(" ", 1)[1].split("_")
        doubles = sum(int(c.split(":")[1]) for c in chains)
        assert doubles == 0, \
            f"{name} carries a double bond — it is an exact isobar of a canonical Cer [M+H-H2O]+"


def test_the_base_ion_is_not_the_dihydro_isobar():
    """268.2999 (m18:0 - H2O), never 266.2842 (d18:0 - 2H2O), which is a real canonical ion."""
    for name, (_, mz, _) in entries().items():
        assert not any(abs(m - 266.2842) < TOL for m in mz), \
            f"{name} writes 266.2842, the doubly-dehydrated dihydro base — that is a canonical ion"
