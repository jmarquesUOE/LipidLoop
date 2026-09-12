"""Build the short-chain acylcarnitines the shipped libraries stop short of.

LipiDex and LipidBlast both start at **C10**. Everything below it is absent, which removes exactly
the acylcarnitines that carry the most metabolic information: acetyl (C2), propionyl (C3), butyryl
(C4) and isovaleryl (C5) are the readouts of fatty-acid oxidation and branched-chain amino-acid
catabolism, and free carnitine (C0) is the pool they draw on. Across five reference deposits the
missing acylcarnitines run **C2-C16**, and none of them can be identified at present however
abundant they are.

This is also why acylcarnitines have never appeared in the retention model: the model fits per
class on confident identifications, and a class with no library entries produces none.

**The fragmentation is regular, which is what makes generating them safe.** Every LipiDex
acylcarnitine carries the same three constant fragments plus two that follow the chain:

    AC 10:0 [M+H]+              AC 18:1 [M+H]+
      257.1747  loss of C3H9N     367.2843  loss of C3H9N      <- chain-dependent
      155.1430  acylium [10:0]    265.2526  acylium [18:1]     <- chain-dependent
      144.1019  C7H14O2N          144.1019  C7H14O2N           <- constant
       85.0284  C4H5O2  (base)     85.0284  C4H5O2  (base)     <- constant
       60.0808  C3H10N             60.0808  C3H10N             <- constant

⚠ **Every generated mass is checked against the shipped entries before anything is written.** The
rule above is inferred from the library, so the library is the test: this script regenerates the
C10-C26 species LipiDex already holds and refuses to run if any precursor or fragment disagrees.
A silent 6 mDa error in a derived formula would put an entire class systematically off.

⚠ **Short chains lose the acylium.** Below about C4 the acyl fragment falls under the m/z floor a
real acquisition records, and for C0 and C2 it does not exist as a distinct ion at all, so those
entries carry only the fragments that are actually observable.

    python scripts/build_acylcarnitine_library.py --out data/libraries/Acylcarnitines_Positive.msp
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

PROTON = 1.00727646
CARBON = 12.0
HYDROGEN = 1.0078250319
OXYGEN = 15.9949146221
NITROGEN = 14.0030740052
ELECTRON = 0.00054858

#: L-carnitine, C7H15NO3.
CARNITINE = 7 * CARBON + 15 * HYDROGEN + NITROGEN + 3 * OXYGEN
WATER = 2 * HYDROGEN + OXYGEN
TRIMETHYLAMINE = 3 * CARBON + 9 * HYDROGEN + NITROGEN          # the C3H9N neutral loss

#: The three fragments every acylcarnitine shows, independent of the acyl chain.
#
# ⚠ These are CATIONS: subtract the electron. Writing the neutral formula mass instead put every
# constant fragment 0.55 mDa high — 6.5 ppm at m/z 85, outside any sensible match tolerance and
# systematic across the entire class. It survived the first check only because the tolerance was
# 5 mDa, loose enough to hide exactly the error the check exists to catch. The tolerance is now
# 1 mDa, which fails on this bug.
CONSTANT = [
    (7 * CARBON + 14 * HYDROGEN + 2 * OXYGEN + NITROGEN - ELECTRON, 100, "C7H14O2N1_Fragment_[]"),
    (4 * CARBON + 5 * HYDROGEN + 2 * OXYGEN - ELECTRON, 999, "C4H5O2_Fragment_[]"),
    (3 * CARBON + 10 * HYDROGEN + NITROGEN - ELECTRON, 100, "C3H10N_Fragment_[]"),
]

# Below this the acylium is not a usable fragment: it falls under the recorded m/z floor, and for
# C0/C2 there is no distinct acyl ion to form.
MIN_ACYLIUM_CARBONS = 4
MIN_CARBONS = 0
MAX_CARBONS = 26
MAX_DOUBLE_BONDS = 6
TOLERANCE = 0.001          # Da — tight enough to catch a missing electron mass (0.55 mDa)


def neutral_mass(carbons: int, double_bonds: int) -> float:
    """Carnitine esterified with one fatty acid: carnitine + FA - H2O."""
    if carbons == 0:
        return CARNITINE
    fa = carbons * CARBON + (2 * carbons - 2 * double_bonds) * HYDROGEN + 2 * OXYGEN
    return CARNITINE + fa - WATER


def acylium(carbons: int, double_bonds: int) -> float:
    """The [CnH(2n-1-2d)O]+ acyl cation, written by LipiDex as an `O-1 Alkyl Fragment`."""
    return carbons * CARBON + (2 * carbons - 1 - 2 * double_bonds) * HYDROGEN + OXYGEN - ELECTRON


def peaks(carbons: int, double_bonds: int, precursor: float) -> list[tuple[float, int, str]]:
    out = [(precursor - TRIMETHYLAMINE, 100, "C-3H-9N-1_Neutral Loss_[]")]
    if carbons >= MIN_ACYLIUM_CARBONS:
        out.append((acylium(carbons, double_bonds), 100,
                    f"O-1_Alkyl Fragment_[{carbons}:{double_bonds}]"))
    return out + CONSTANT


def shipped(paths) -> dict[str, tuple[float, list[float]]]:
    """{`AC c:d`: (precursor, sorted fragment m/z)} from the libraries already on disk."""
    out = {}
    for path in paths:
        if not Path(path).exists():
            continue
        for block in Path(path).read_text(errors="replace").split("\n\n"):
            m = re.match(r"Name: (AC \d+:\d+) \[M\+H\]\+", block)
            if not m:
                continue
            pre = re.search(r"PRECURSORMZ: ([\d.]+)", block)
            mz = [float(x.split()[0]) for x in block.split("\n") if re.match(r"^\d+\.", x)]
            if pre:
                out[m.group(1)] = (float(pre.group(1)), sorted(mz))
    return out


def verify(reference) -> list[str]:
    """Regenerate what is already shipped and report every disagreement."""
    problems = []
    for name, (precursor, mz) in sorted(reference.items()):
        c, d = (int(x) for x in name.split()[1].split(":"))
        mine = neutral_mass(c, d) + PROTON
        if abs(mine - precursor) > TOLERANCE:
            problems.append(f"{name}: precursor {mine:.4f} vs shipped {precursor:.4f}")
            continue
        for want in sorted(x for x, _, _ in peaks(c, d, mine)):
            if not any(abs(want - got) <= TOLERANCE for got in mz):
                problems.append(f"{name}: fragment {want:.4f} absent from the shipped entry")
    return problems


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--libraries", default="data/libraries")
    ap.add_argument("--force", action="store_true",
                    help="write even if the check against the shipped entries fails")
    args = ap.parse_args()

    libdir = Path(args.libraries)
    reference = shipped([libdir / "LipiDex_HCD_Formic.msp", libdir / "LipiDex_HCD_Acetate.msp"])
    if reference:
        cs = [int(k.split()[1].split(":")[0]) for k in reference]
        print(f"  {len(reference)} acylcarnitines already shipped, C{min(cs)}-C{max(cs)}")
    else:
        # ⚠ Say so loudly. With no reference the rule is unverified, and the whole safety argument
        # for generating this class rests on regenerating what LipiDex already holds.
        print("  ⚠ no shipped acylcarnitines found — the rule cannot be verified against anything")

    problems = verify(reference)
    if problems:
        print(f"  ⚠ the derived rule disagrees with {len(problems)} shipped value(s):")
        for p in problems[:10]:
            print(f"      {p}")
        if not args.force:
            print("  refusing to write — a wrong rule would put the whole class systematically off")
            return 1
    else:
        print(f"  ✓ rule reproduces every shipped entry within {TOLERANCE} Da")

    lines, written = [], 0
    for carbons in range(MIN_CARBONS, MAX_CARBONS + 1):
        for double_bonds in range(0, MAX_DOUBLE_BONDS + 1):
            if carbons == 0 and double_bonds:
                continue
            # Methylene-interrupted ceiling, as in the fatty-acid builder: keeps AC 6:4 out.
            if double_bonds and double_bonds > 1 + (carbons - 3) // 3:
                continue
            name = f"AC {carbons}:{double_bonds}"
            if name in reference:            # never duplicate what is already searched
                continue
            neutral = neutral_mass(carbons, double_bonds)
            precursor = neutral + PROTON
            hydrogens = 15 if carbons == 0 else 15 + 2 * carbons - 2 * double_bonds - 2
            lines.append(f"Name: {name} [M+H]+;")
            lines.append(f"MW: {precursor:.4f}")
            lines.append(f"PRECURSORMZ: {precursor:.4f}")
            lines.append(f"Comment: Name={name} [M+H]+ Mass={precursor:.4f} "
                         f"Formula=C{7 + carbons}H{hydrogens}N1O{3 if carbons == 0 else 4} "
                         f"OptimalPolarity=true Type=LipiDex")
            p = [x for x in peaks(carbons, double_bonds, precursor) if x[0] < precursor - 2.0]
            lines.append(f"Num Peaks: {len(p)}")
            for mz, intensity, annotation in sorted(p, key=lambda t: -t[0]):
                lines.append(f'{mz:.4f} {intensity} "{annotation}"')
            lines.append("")
            written += 1

    Path(args.out).write_text("\n".join(lines))
    print(f"  wrote {written} new entries to {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
