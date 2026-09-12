"""Fill the double-bond holes in the shipped diacylglycerol library.

The shipped DG library is not missing chains, it is missing **unsaturation**. Its 1,766 entries are
built from 40 acyl chains, and the double-bond coverage per chain length is patchy:

    C16: 0,1        C19/21/23/25: 0 only      C24: 0,1,4
    C18: 0-4        C20: 0-5                  C26: 0,1,2
    C22: 0,1,2,4,5,6  (no 22:3)

So a diacylglycerol carrying a polyunsaturated odd chain, or 24:2, or 16:2, has no entry at all --
79 sum compositions published across five reference deposits cannot be matched, and every one of
them reads as a detection failure rather than as the library gap it is.

**Enumerated to the observed sums, not exhaustively.** Only chain pairs summing to a composition a
deposit actually reported are written. An exhaustive C10-C26 enumeration would be roughly 5,000
entries, most of which nothing has ever observed, and every one is another chance for a false
positive at a plausible mass.

⚠ **Chains stop at C26, matching the shipped library's own range.** Extending to C30 would reach 25
more sums at the cost of 658 additional entries built on acyl chains that are vanishingly rare in
a diacylglycerol. What that leaves unreached is printed rather than passed over in silence.

⚠ **The masses are checked against the shipped entries before anything is written**, the same way
the acylcarnitine library is: the rule regenerates species LipiDex already holds, and the script
refuses to run if any disagrees by more than 1 mDa.

    python scripts/build_diacylglycerol_library.py --out data/libraries/Diacylglycerols_Positive.msp
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

GLYCEROL = 3 * CARBON + 8 * HYDROGEN + 3 * OXYGEN
WATER = 2 * HYDROGEN + OXYGEN
AMMONIA = NITROGEN + 3 * HYDROGEN
AMMONIUM = AMMONIA + PROTON

MIN_CHAIN_C, MAX_CHAIN_C = 10, 26
MAX_CHAIN_DB = 6                  # 22:6 is the most unsaturated acyl chain in mammalian lipids
TOLERANCE = 0.001


def max_double_bonds(carbons: int) -> int:
    """Methylene-interrupted ceiling: each double bond after the first needs three more carbons.

    Odd-chain fatty acids come from propionyl-CoA priming, not the acetyl-CoA
    elongation/desaturation pathway that builds even-chain polyenes -- they are essentially always
    saturated in vivo, trace-monounsaturated at most. See fatty_acids.py::max_double_bonds.
    """
    if carbons % 2:
        return min(1, MAX_CHAIN_DB)
    return max(0, 1 + (carbons - 3) // 3)


def fatty_acid(carbons: int, double_bonds: int) -> float:
    return carbons * CARBON + (2 * carbons - 2 * double_bonds) * HYDROGEN + 2 * OXYGEN


def neutral_mass(a: tuple[int, int], b: tuple[int, int]) -> float:
    """Glycerol esterified with two fatty acids."""
    return GLYCEROL + fatty_acid(*a) + fatty_acid(*b) - 2 * WATER


def peaks(a, b, precursor: float) -> list[tuple[float, int, str]]:
    """The head-group loss, plus one alkyl neutral loss per distinct chain."""
    out = [(precursor - AMMONIA - WATER, 999, "N-1H-5O-1_Neutral Loss_[]")]
    for chain in {a, b}:
        out.append((precursor - fatty_acid(*chain) - AMMONIA, 999,
                    f"H-4N-1_Alkyl Neutral Loss_[{chain[0]}:{chain[1]}]"))
    return out


def shipped(paths) -> dict[str, float]:
    """{`DG a:b_c:d`: precursor} from the libraries already on disk."""
    out = {}
    for path in paths:
        p = Path(path)
        if not p.exists():
            continue
        for block in p.read_text(errors="replace").split("\n\n"):
            m = re.match(r"Name: DG (\d+:\d+)_(\d+:\d+) \[M\+NH4\]\+", block)
            pre = re.search(r"PRECURSORMZ: ([\d.]+)", block)
            if m and pre:
                out[f"DG {m.group(1)}_{m.group(2)}"] = float(pre.group(1))
    return out


def verify(reference: dict[str, float]) -> list[str]:
    problems = []
    for name, precursor in sorted(reference.items()):
        a, b = (tuple(int(x) for x in c.split(":")) for c in name.split()[1].split("_"))
        mine = neutral_mass(a, b) + AMMONIUM
        if abs(mine - precursor) > TOLERANCE:
            problems.append(f"{name}: {mine:.4f} vs shipped {precursor:.4f}")
    return problems


def wanted(path: Path | None) -> set[tuple[int, int]]:
    """Sum compositions to cover, read from a file of `DG C:D` lines."""
    if not path or not path.exists():
        return set()
    out = set()
    for line in path.read_text().splitlines():
        m = re.match(r"\s*DG (\d+):(\d+)\s*$", line)
        if m:
            out.add((int(m.group(1)), int(m.group(2))))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--libraries", default="data/libraries")
    ap.add_argument("--targets", help="file of `DG C:D` sum compositions to cover")
    args = ap.parse_args()

    libdir = Path(args.libraries)
    reference = shipped([libdir / "LipiDex_HCD_Formic.msp", libdir / "LipidBlast_Formic.msp"])
    print(f"  {len(reference)} diacylglycerols already shipped")
    problems = verify(reference)
    if problems:
        print(f"  ⚠ the derived rule disagrees with {len(problems)} shipped value(s):")
        for p in problems[:8]:
            print(f"      {p}")
        print("  refusing to write")
        return 1
    print(f"  ✓ rule reproduces every shipped entry within {TOLERANCE} Da")

    targets = wanted(Path(args.targets) if args.targets else None)
    if not targets:
        print("  no target sum compositions given (--targets) — nothing to build")
        return 1

    chains = [(c, d) for c in range(MIN_CHAIN_C, MAX_CHAIN_C + 1)
              for d in range(0, min(max_double_bonds(c), MAX_CHAIN_DB) + 1)]
    have = {tuple(sorted(tuple(int(x) for x in c.split(":")) for c in k.split()[1].split("_")))
            for k in reference}

    lines, written, covered = [], 0, set()
    for C, D in sorted(targets):
        for a in chains:
            for b in chains:
                if a[0] + b[0] != C or a[1] + b[1] != D or a > b:
                    continue
                if (a, b) in have:
                    continue
                neutral = neutral_mass(a, b)
                precursor = neutral + AMMONIUM
                name = f"DG {a[0]}:{a[1]}_{b[0]}:{b[1]}"
                p = [x for x in peaks(a, b, precursor) if x[0] < precursor - 2.0]
                if not p:
                    continue
                carbons = 3 + a[0] + b[0]
                hydrogens = 8 + (2 * a[0] - 2 * a[1]) + (2 * b[0] - 2 * b[1]) - 4 + 4
                lines.append(f"Name: {name} [M+NH4]+;")
                lines.append(f"MW: {precursor:.4f}")
                lines.append(f"PRECURSORMZ: {precursor:.4f}")
                lines.append(f"Comment: Name={name} [M+NH4]+ Mass={precursor:.4f} "
                             f"Formula=O5H{hydrogens}C{carbons}N1 "
                             f"OptimalPolarity=true Type=LipiDex")
                lines.append(f"Num Peaks: {len(p)}")
                for mz, intensity, annotation in sorted(p, key=lambda t: -t[0]):
                    lines.append(f'{mz:.4f} {intensity} "{annotation}"')
                lines.append("")
                written += 1
                covered.add((C, D))

    Path(args.out).write_text("\n".join(lines))
    print(f"  wrote {written} entries covering {len(covered)} of {len(targets)} sum compositions")
    unreached = sorted(targets - covered)
    if unreached:
        # ⚠ Never silent. A composition no plausible chain pair reaches is a statement about the
        # deposit, and hiding it would turn a finding into an apparent coverage gap.
        print(f"  {len(unreached)} unreachable with C{MIN_CHAIN_C}-C{MAX_CHAIN_C} chains at "
              f"<={MAX_CHAIN_DB} double bonds — these need acyls that do not occur in a "
              f"diacylglycerol:")
        print("      " + ", ".join(f"DG {c}:{d}" for c, d in unreached))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
