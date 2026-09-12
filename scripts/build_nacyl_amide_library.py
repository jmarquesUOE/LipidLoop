"""Build the N-acyl amides — the endocannabinoid-adjacent classes no shipped library carries.

Four related classes, all a fatty acid amide-bonded to a small polar head, and none present in any
library here: **NAE** (N-acylethanolamines, the anandamide family), **NAGly** (N-acylglycines),
**NATau** (N-acyltaurines) and **GPNAE** (glycerophospho-N-acylethanolamines, the direct NAE
precursor). Reference deposits report 24 species across the four; every one is unidentifiable.

⚠ **This library is weaker evidence than the acylcarnitine or diacylglycerol ones, and the
difference is worth stating plainly.** Those two were generated from a rule inferred from LipiDex's
own entries and then checked by regenerating them, so the library validated itself. There is no
N-acyl amide anywhere in the shipped libraries, so **there is nothing here to check the
fragmentation against**. What is verified:

  * **precursor masses** — pure formula arithmetic, cross-checked against the exact masses of
    named compounds the classes are known by (anandamide, oleoylethanolamide, N-oleoylglycine,
    N-oleoyltaurine). Those are chemistry and they are solid.
  * **fragmentation** — taken from the published behaviour of these classes, NOT measured here and
    NOT validated against an in-house spectrum.

So a hit from this library rests on accurate mass plus a head-group ion, and should be read as
weaker than a hit from a library whose fragments were verified. Injecting authentic standards would
settle it; until then this is chemistry, not measurement.

## The head groups

| class | head | polarity | diagnostic ion |
|---|---|---|---|
| NAE | ethanolamine | + | *m/z* 62.0600, protonated ethanolamine |
| NAGly | glycine | + | *m/z* 76.0393, protonated glycine |
| NATau | taurine | - | *m/z* 124.0074, deprotonated taurine |
| GPNAE | glycerophospho-ethanolamine | + | *m/z* 216.0632, the GPE head |

Each also gives the acylium, which is what carries the chain identity.

    python scripts/build_nacyl_amide_library.py --out-pos data/libraries/NAcylAmides_Positive.msp \\
                                                --out-neg data/libraries/NAcylAmides_Negative.msp
"""
from __future__ import annotations

import argparse
from pathlib import Path

PROTON = 1.00727646
ELECTRON = 0.00054858
CARBON = 12.0
HYDROGEN = 1.0078250319
OXYGEN = 15.9949146221
NITROGEN = 14.0030740052
SULFUR = 31.97207069
PHOSPHORUS = 30.97376151
WATER = 2 * HYDROGEN + OXYGEN

#: class -> (head-group neutral formula mass, polarity, diagnostic ion m/z, label)
HEADS = {
    # ethanolamine C2H7NO
    "NAE":   (2 * CARBON + 7 * HYDROGEN + NITROGEN + OXYGEN, "Pos", 62.0600, "C2H8NO"),
    # glycine C2H5NO2
    "NAGly": (2 * CARBON + 5 * HYDROGEN + NITROGEN + 2 * OXYGEN, "Pos", 76.0393, "C2H6NO2"),
    # taurine C2H7NO3S
    "NATau": (2 * CARBON + 7 * HYDROGEN + NITROGEN + 3 * OXYGEN + SULFUR, "Neg", 124.0074,
              "C2H6NO3S"),
    # glycerophosphoethanolamine C5H14NO6P
    "GPNAE": (5 * CARBON + 14 * HYDROGEN + NITROGEN + 6 * OXYGEN + PHOSPHORUS, "Pos", 216.0632,
              "C5H15NO6P"),
}

#: Exact masses of named compounds, used to check the arithmetic against something outside it.
KNOWN = {
    ("NAE", 20, 4): 348.2897,     # anandamide
    ("NAE", 18, 1): 326.3054,     # oleoylethanolamide
    ("NAE", 16, 0): 300.2897,     # palmitoylethanolamide
    ("NAGly", 18, 1): 340.2846,   # N-oleoylglycine
    # ⚠ 388.2527, not 390.2320. The first version of this list carried the wrong value and the
    # check correctly refused to write — catching an error in the REFERENCE rather than in the
    # arithmetic, which is the failure mode a self-check cannot see. C20H39NO4S = 389.2600
    # neutral; [M-H]- follows.
    ("NATau", 18, 1): 388.2527,   # N-oleoyltaurine, [M-H]-
}

MIN_C, MAX_C = 8, 26
MAX_DB = 6
TOLERANCE = 0.001


def max_double_bonds(carbons: int) -> int:
    # Odd-chain fatty acids come from propionyl-CoA priming, not the acetyl-CoA
    # elongation/desaturation pathway that builds even-chain polyenes -- they are essentially
    # always saturated in vivo, trace-monounsaturated at most. See fatty_acids.py::max_double_bonds.
    if carbons % 2:
        return min(1, MAX_DB)
    return max(0, 1 + (carbons - 3) // 3)


def fatty_acid(carbons: int, double_bonds: int) -> float:
    return carbons * CARBON + (2 * carbons - 2 * double_bonds) * HYDROGEN + 2 * OXYGEN


def neutral_mass(cls: str, carbons: int, double_bonds: int) -> float:
    """Amide bond: fatty acid + head - water."""
    return fatty_acid(carbons, double_bonds) + HEADS[cls][0] - WATER


def precursor(cls: str, carbons: int, double_bonds: int) -> float:
    neutral = neutral_mass(cls, carbons, double_bonds)
    return neutral - PROTON if HEADS[cls][1] == "Neg" else neutral + PROTON


def acylium(carbons: int, double_bonds: int) -> float:
    return carbons * CARBON + (2 * carbons - 1 - 2 * double_bonds) * HYDROGEN + OXYGEN - ELECTRON


def check_known() -> list[str]:
    out = []
    for (cls, c, d), expected in sorted(KNOWN.items()):
        got = precursor(cls, c, d)
        if abs(got - expected) > TOLERANCE:
            out.append(f"{cls} {c}:{d}: {got:.4f} vs published {expected:.4f}")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-pos", required=True)
    ap.add_argument("--out-neg", required=True)
    args = ap.parse_args()

    problems = check_known()
    if problems:
        print("  ⚠ the mass arithmetic disagrees with published compound masses:")
        for p in problems:
            print(f"      {p}")
        print("  refusing to write")
        return 1
    print(f"  ✓ precursors reproduce all {len(KNOWN)} named reference compounds within "
          f"{TOLERANCE} Da")

    out = {"Pos": [], "Neg": []}
    counts = {}
    for cls, (_, polarity, diagnostic, label) in HEADS.items():
        n = 0
        for carbons in range(MIN_C, MAX_C + 1):
            for double_bonds in range(0, min(max_double_bonds(carbons), MAX_DB) + 1):
                mz = precursor(cls, carbons, double_bonds)
                name = f"{cls} {carbons}:{double_bonds}"
                adduct = "[M-H]-" if polarity == "Neg" else "[M+H]+"
                peaks = [(diagnostic, 999, f"{label}_Fragment_[]")]
                # ⚠ Every entry needs a CHAIN-SPECIFIC fragment, in both polarities.
                #
                # Without one there is nothing distinguishing `NATau 18:1` from `NATau 20:1` beyond
                # the precursor, and — the part that actually bites — a decoy of such an entry is
                # identical to its target, so `make_decoy_library` writes none. The first version
                # gave the negative arm only the taurine ion, and 111 NATau targets would have
                # entered the search with zero decoy competition, biasing the FDR low on exactly
                # the class with the weakest evidence.
                if polarity == "Pos":
                    peaks.append((acylium(carbons, double_bonds), 500,
                                  f"O-1_Alkyl Fragment_[{carbons}:{double_bonds}]"))
                else:
                    # The acyl leaves as the carboxylate anion, which carries the chain identity.
                    peaks.append((fatty_acid(carbons, double_bonds) - PROTON, 500,
                                  f"FA_Fragment_[{carbons}:{double_bonds}]"))
                peaks = [p for p in peaks if p[0] < mz - 2.0]
                if not peaks:
                    continue
                out[polarity].append(f"Name: {name} {adduct};")
                out[polarity].append(f"MW: {mz:.4f}")
                out[polarity].append(f"PRECURSORMZ: {mz:.4f}")
                out[polarity].append(f"Comment: Name={name} {adduct} Mass={mz:.4f} "
                                     f"OptimalPolarity=true Type=LipiDex")
                out[polarity].append(f"Num Peaks: {len(peaks)}")
                for m, intensity, annotation in sorted(peaks, key=lambda t: -t[0]):
                    out[polarity].append(f'{m:.4f} {intensity} "{annotation}"')
                out[polarity].append("")
                n += 1
        counts[cls] = n

    Path(args.out_pos).write_text("\n".join(out["Pos"]))
    Path(args.out_neg).write_text("\n".join(out["Neg"]))
    for cls, n in counts.items():
        print(f"    {cls:<7} {n:3d} entries ({HEADS[cls][1]})")
    print(f"  wrote {len(counts)} classes: {args.out_pos}, {args.out_neg}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
