"""Build a negative-mode free-fatty-acid library, at the resolution the fragmentation supports.

No LipiDex library contains a single free fatty acid. Across 91 FA masses the default set has
**zero** negative-mode entries, so a free fatty acid cannot be identified at all, however abundant
it is — and in the skin organoid data they are very abundant indeed: `FA 18:1` is a 1.9e8 peak
present in 63 of 63 files, and every species from C14 to C26 is there, all of them unidentified.

**The hard part is that free fatty acids barely fragment.** Pooling every MS2 taken at a fatty-acid
precursor mass across 16 files, and asking which product ions recur:

| | MS2 spectra | `[M-H-CO2]-` present | relative intensity |
|---|---|---|---|
| FA 16:0 | 711 | **0%** | — |
| FA 18:0 | 561 | **0%** | — |
| FA 18:1 | 305 | **0%** | — |
| FA 24:0 | 21 | **0%** | — |
| FA 20:3 | 11 | **91%** | 999 |
| FA 20:4 | 15 | **100%** | 999 |
| FA 22:6 | 6 | **100%** | 927 |

The carboxylate anion is stable: saturated and monounsaturated acids produce nothing reproducible
at these collision energies, and what looks like a spectrum is co-isolated background — with no
real fragment to normalise against, whatever else is in the isolation window becomes the base peak.
Only polyunsaturated acids, **DB >= 3**, lose CO2, and when they do it is the base peak.

So this library contains **only the species that produce a fragment**, one entry each, one peak
each. Writing entries for the rest would be writing a library of noise, and every one of them would
match any spectrum at that mass.

**The other species are not lost — they are identified by retention instead.** See
`fatty_acids.py`, which names them from accurate mass plus a fitted `RT ~ C + DB` model. That is
the right tool for a class whose members are separated by chromatography and not by fragmentation.

    python scripts/build_fatty_acid_library.py --out data/libraries/FreeFattyAcids_Negative.msp
"""
from __future__ import annotations

import argparse
from pathlib import Path

PROTON = 1.00727646
CARBON = 12.0
HYDROGEN = 1.0078250319
OXYGEN = 15.9949146221
CO2 = 43.98983

# Below this the acid does not lose CO2 at all — measured, see the module docstring.
MIN_DOUBLE_BONDS = 3
MIN_CARBONS = 14
MAX_CARBONS = 30
MAX_DOUBLE_BONDS = 6
# The library reader drops library peaks within this of the precursor. CO2 loss is 44 Da, so it
# always survives; the check is here so a future fragment cannot be added that silently does not.
PRECURSOR_CUTOFF = 2.0


def neutral_mass(carbons: int, double_bonds: int) -> float:
    hydrogens = 2 * carbons - 2 * double_bonds
    return carbons * CARBON + hydrogens * HYDROGEN + 2 * OXYGEN


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-double-bonds", type=int, default=MIN_DOUBLE_BONDS)
    args = ap.parse_args()

    lines, written = [], 0
    for carbons in range(MIN_CARBONS, MAX_CARBONS + 1):
        for double_bonds in range(args.min_double_bonds, MAX_DOUBLE_BONDS + 1):
            # Methylene-interrupted ceiling: each double bond after the first needs three
            # more carbons. Keeps FA 14:6 out of the library, whose mass would otherwise be
            # enumerated and would then claim a real feature.
            if double_bonds > 1 + (carbons - 3) // 3:
                continue
            neutral = neutral_mass(carbons, double_bonds)
            precursor = neutral - PROTON
            fragment = precursor - CO2
            if precursor - fragment <= PRECURSOR_CUTOFF:
                continue
            name = f"FA {carbons}:{double_bonds}"
            lines.append(f"Name: {name} [M-H]-;")
            lines.append(f"MW: {neutral:.4f}")
            lines.append(f"PRECURSORMZ: {precursor:.4f}")
            # No moiety annotation: a free fatty acid IS its own chain, so there is nothing for
            # the purity calculation to resolve. Purity returns 0 and the name is reported as it
            # stands, which is the correct behaviour rather than a limitation.
            lines.append(f"Comment: Name={name} [M-H]- Mass={precursor:.4f} "
                         f"Formula=C{carbons}H{2 * carbons - 2 * double_bonds}O2 "
                         f"OptimalPolarity=true Type=LipiDex")
            lines.append("Num Peaks: 1")
            lines.append(f'{fragment:.4f} 999 "CO2_Loss_Fragment_[]"')
            lines.append("")
            written += 1

    Path(args.out).write_text("\n".join(lines))
    print(f"{written} polyunsaturated fatty acids written to {args.out}")
    print(f"  C{MIN_CARBONS}-C{MAX_CARBONS}, {args.min_double_bonds}-{MAX_DOUBLE_BONDS} double "
          f"bonds, [M-H]- with one fragment: [M-H-CO2]-")
    print("  saturated and monounsaturated acids are deliberately absent — they produce no "
          "reproducible fragment, so they are named by retention instead (fatty_acids.py)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
