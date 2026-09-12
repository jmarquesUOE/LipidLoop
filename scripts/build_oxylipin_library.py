"""Build a negative-mode oxylipin library and the inclusion list it needs to be any use.

**Read this before adding the library to a run.** On both validation sets the evidence says a
spectral library alone cannot deliver oxylipins, and the reason is acquisition, not identification.

## What the data says

Oxylipins are detectable at MS1 level and almost never fragmented. In skin negative mode:

| | MS1 detection | MS2 spectra across 10 files |
|---|---|---|
| HODE (m/z 295.2279) | 40 of 63 files, area 1.4e7 | **0** |
| HETE/EET (319.2279) | 27 files, area 5.7e6 | **1** |
| 9-HpODE (311.2228) | 2 files | **0** |
| DiHOME (313.2384) | 2 files | 11, all background fragments |

Most of the rest — HEPE, DiHETrE, HDHA, PGF2a, TXB2, RvD1, LXA4, PGA2, 12-HHT — are **not detected
at all**, in either dataset. The ones that are "detected" are often flagged `Adduct of existing
peak`, so even the MS1 evidence is thin.

This is the cholesteryl-ester problem in a more severe form. Oxylipins sit two to three orders of
magnitude below the structural lipids that win every top-N decision, on a run where 46-67% of the
duty cycle already goes to contamination (docs/EXCLUSION_LIST.md). A library cannot identify a
spectrum that was never acquired.

## So this script writes two files

`--out` is the library, and it is **deliberately not in the default set**. Add it only when the
acquisition has changed, or a targeted method is being processed.

`--inclusion-list` is the actionable half: the precursors to force the instrument to fragment. That
is the change that would make the library worth having.

## The fragments, and how far to trust them

Every oxylipin gets the neutral losses that are generic to the class and safe to compute:
`[M-H-H2O]-`, `[M-H-2H2O]-`, `[M-H-CO2]-`, `[M-H-H2O-CO2]-`.

Each hydroxy species also gets **two candidate alpha-cleavage ions** at the carbinol carbon — the
carboxyl-containing fragment that distinguishes 5- from 12- from 15-HETE, which is the whole point
of oxylipin biology. Two, because which one forms depends on a mechanism nothing here can settle:
the alpha,beta-unsaturated carbonyl (full intensity; this is what published values for these
compounds correspond to, and it puts 5-HETE at 115.0395 and 9-HODE at 171.1021) and the
alcohol-retaining form two hydrogens heavier (reduced). Writing both lets a spectrum match on
whichever it actually produces rather than on an assumption.

**Every alpha-cleavage mass here is computed from structure, not measured.** Nothing in these two
datasets could validate them — there are no spectra to validate against. Treat every positional
assignment from this library as provisional until confirmed with authentic standards, which for
oxylipins is the norm rather than a counsel of perfection: the positional isomers co-elute closely
and the biology turns on which one it is.

    python scripts/build_oxylipin_library.py \\
        --out data/libraries/Oxylipins_Negative.msp \\
        --inclusion-list method/Oxylipin_inclusion_list.csv
"""
from __future__ import annotations

import argparse
from pathlib import Path

PROTON = 1.00727646
H2O = 18.0105646
CO2 = 43.9898292
CARBON, HYDROGEN, OXYGEN = 12.0, 1.0078250319, 15.9949146221

# (name, carbons, double bonds, oxygens, hydroxyl position or None, retention hint)
# Hydroxyl position drives the alpha-cleavage ion; None means the class is enumerated without a
# positional claim, and no alpha-cleavage ion is written for it.
OXYLIPINS = [
    # --- linoleic acid derived, C18 ---
    ("9-HODE", 18, 2, 3, 9), ("13-HODE", 18, 2, 3, 13),
    ("9-oxoODE", 18, 2, 3, None), ("13-oxoODE", 18, 2, 3, None),
    ("9-HOTrE", 18, 3, 3, 9), ("13-HOTrE", 18, 3, 3, 13),
    ("9,10-EpOME", 18, 2, 3, None), ("12,13-EpOME", 18, 2, 3, None),
    ("9,10-DiHOME", 18, 1, 4, None), ("12,13-DiHOME", 18, 1, 4, None),
    ("9-HpODE", 18, 2, 4, 9), ("13-HpODE", 18, 2, 4, 13),
    # --- arachidonic acid derived, C20:4 ---
    ("5-HETE", 20, 4, 3, 5), ("8-HETE", 20, 4, 3, 8), ("9-HETE", 20, 4, 3, 9),
    ("11-HETE", 20, 4, 3, 11), ("12-HETE", 20, 4, 3, 12), ("15-HETE", 20, 4, 3, 15),
    ("20-HETE", 20, 4, 3, 20),
    ("5-oxoETE", 20, 4, 3, None), ("15-oxoETE", 20, 4, 3, None),
    ("5,6-EET", 20, 4, 3, None), ("8,9-EET", 20, 4, 3, None),
    ("11,12-EET", 20, 4, 3, None), ("14,15-EET", 20, 4, 3, None),
    ("5,6-DiHETrE", 20, 3, 4, None), ("8,9-DiHETrE", 20, 3, 4, None),
    ("11,12-DiHETrE", 20, 3, 4, None), ("14,15-DiHETrE", 20, 3, 4, None),
    ("LTB4", 20, 4, 4, 5), ("12-HHT", 17, 3, 3, 12),
    # --- prostanoids, C20 ---
    ("PGE2", 20, 2, 5, None), ("PGD2", 20, 2, 5, None), ("PGF2a", 20, 2, 5, None),
    ("PGA2", 20, 3, 4, None), ("PGJ2", 20, 3, 4, None),
    ("8-iso-PGF2a", 20, 2, 5, None), ("6-keto-PGF1a", 20, 1, 6, None),
    ("TXB2", 20, 2, 6, None),
    ("PGE1", 20, 1, 5, None), ("PGF1a", 20, 1, 5, None),
    # --- EPA and DHA derived ---
    ("5-HEPE", 20, 5, 3, 5), ("12-HEPE", 20, 5, 3, 12), ("15-HEPE", 20, 5, 3, 15),
    ("18-HEPE", 20, 5, 3, 18),
    ("4-HDHA", 22, 6, 3, 4), ("14-HDHA", 22, 6, 3, 14), ("17-HDHA", 22, 6, 3, 17),
    ("RvD1", 22, 6, 5, None), ("RvE1", 20, 5, 5, None),
    ("LXA4", 20, 4, 5, None), ("LXB4", 20, 4, 5, None),
    ("Maresin 1", 22, 6, 4, None), ("PD1", 22, 6, 4, None),
]

# Relative intensities for the generic losses. Not measured here — there are no spectra to measure
# — so they are set to describe the ordinary appearance of a hydroxy fatty acid in negative mode
# and are flat enough that the dot product is driven by which peaks are PRESENT rather than by
# their ratios, which is the correct behaviour for a library that has not been calibrated.
INTENSITY_ALPHA = 999
INTENSITY_ALPHA_ALCOHOL = 400
INTENSITY_MINUS_H2O = 500
INTENSITY_MINUS_H2O2 = 200
INTENSITY_MINUS_CO2 = 300
INTENSITY_MINUS_H2O_CO2 = 250

MIN_FRAGMENT_MASS = 61.0
PRECURSOR_CUTOFF = 2.0


def neutral_mass(carbons: int, double_bonds: int, oxygens: int) -> float:
    hydrogens = 2 * carbons - 2 * double_bonds
    return carbons * CARBON + hydrogens * HYDROGEN + oxygens * OXYGEN


def alpha_cleavage(position: int, double_bonds_below: int) -> tuple:
    """The two candidate carboxyl-side fragments at the carbinol carbon, as anions.

    Cleaving the C-C bond next to the hydroxyl-bearing carbon leaves the acid, the chain up to it,
    and the oxygen. Which ion results depends on the mechanism, and **the mechanism is not
    determined by anything in these datasets** — there are no oxylipin spectra here to determine
    it with:

      * retaining the alcohol:            C_n H_(2n-1-2d) O_3
      * with the further loss of H2 (the alpha,beta-unsaturated carbonyl, which is what the
        literature values for these compounds correspond to): C_n H_(2n-3-2d) O_3

    Both are written, the second at full intensity and the first reduced, so a spectrum matches on
    whichever it actually produces rather than on an assumption. That is the correct behaviour for
    a library nothing has calibrated. Confirm with standards before trusting a positional call.
    """
    base = position * CARBON + 3 * OXYGEN
    alcohol = base + (2 * position - 1 - 2 * double_bonds_below) * HYDROGEN
    return alcohol - 2 * HYDROGEN, alcohol


def double_bonds_below(chain_carbons: int, position: int, double_bonds: int) -> int:
    """How many double bonds lie between C1 and the hydroxyl carbon.

    Methylene-interrupted polyenes in these families begin at C5 (arachidonic, EPA), C9 (linoleic,
    alpha-linolenic) or C4 (DHA), so this walks the standard positions rather than guessing. It is
    a structural assumption; a wrong one moves only the alpha-cleavage masses, never the generic
    losses.
    """
    first = {20: 5, 18: 9, 22: 4, 17: 5}.get(chain_carbons, 99)
    return sum(1 for k in range(double_bonds) if first + 3 * k < position)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--inclusion-list")
    ap.add_argument("--rt-start", type=float, default=8.0)
    ap.add_argument("--rt-end", type=float, default=14.0)
    args = ap.parse_args()

    lines, entries = [], []
    for name, carbons, double_bonds, oxygens, position in OXYLIPINS:
        neutral = neutral_mass(carbons, double_bonds, oxygens)
        precursor = neutral - PROTON
        peaks = [
            (precursor - H2O, INTENSITY_MINUS_H2O, "H-2O-1_Neutral Loss_[]"),
            (precursor - 2 * H2O, INTENSITY_MINUS_H2O2, "H-4O-2_Neutral Loss_[]"),
            (precursor - CO2, INTENSITY_MINUS_CO2, "C-1O-2_Neutral Loss_[]"),
            (precursor - H2O - CO2, INTENSITY_MINUS_H2O_CO2, "C-1H-2O-3_Neutral Loss_[]"),
        ]
        if position is not None:
            below = double_bonds_below(carbons, position, double_bonds)
            unsaturated, alcohol = alpha_cleavage(position, below)
            peaks.append((unsaturated, INTENSITY_ALPHA, f"C{position}_Alpha Cleavage_[]"))
            peaks.append((alcohol, INTENSITY_ALPHA_ALCOHOL,
                          f"C{position}_Alpha Cleavage OH_[]"))
        peaks = [p for p in peaks
                 if precursor - p[0] > PRECURSOR_CUTOFF and p[0] > MIN_FRAGMENT_MASS]
        if not peaks:
            continue
        peaks.sort()

        lines.append(f"Name: {name} [M-H]-;")
        lines.append(f"MW: {neutral:.4f}")
        lines.append(f"PRECURSORMZ: {precursor:.4f}")
        lines.append(f"Comment: Name={name} [M-H]- Mass={precursor:.4f} "
                     f"Formula=C{carbons}H{2 * carbons - 2 * double_bonds}O{oxygens} "
                     f"OptimalPolarity=true Type=LipiDex")
        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, intensity, annotation in peaks:
            lines.append(f'{mz:.4f} {intensity} "{annotation}"')
        lines.append("")
        entries.append((name, precursor, carbons, double_bonds, oxygens))

    Path(args.out).write_text("\n".join(lines))
    print(f"{len(entries)} oxylipins written to {args.out}")
    print("  NOT in the default library set — see the module docstring. On both validation sets "
          "these compounds are detected at MS1 and essentially never fragmented.")

    if args.inclusion_list:
        # One row per distinct precursor: the instrument needs masses, not names, and several of
        # these are isomers sharing one.
        by_mass: dict = {}
        for name, precursor, *_ in entries:
            by_mass.setdefault(round(precursor, 4), []).append(name)
        rows = ["m/z,Charge,Polarity,RT start (min),RT end (min),Species,Note"]
        for precursor, names in sorted(by_mass.items()):
            rows.append(f"{precursor:.4f},1,-,{args.rt_start:.1f},{args.rt_end:.1f},"
                        f"{len(names)},\"{'; '.join(names)}\"")
        Path(args.inclusion_list).parent.mkdir(parents=True, exist_ok=True)
        Path(args.inclusion_list).write_text("\n".join(rows) + "\n")
        print(f"\n{len(by_mass)} distinct precursors written to {args.inclusion_list}")
        print("  This is the actionable half. The retention window is a placeholder from where "
              "the C18/C20 acids elute on this method — narrow it with standards before using it.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
