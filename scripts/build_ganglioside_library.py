"""Build the negative-mode ganglioside library from first principles, with no third-party input.

This replaces `build_ganglioside_negative_library.py`, which enumerated its species by reading
LipiDex 2's positive-mode ganglioside library. Nothing here is taken from another library: the
glycan compositions are the Svennerholm structures, the ceramide masses are computed from the
elemental formula, and the three product ions are the ones measured in this laboratory's own
spectra (see docs/GANGLIOSIDES.md). Why it is written at sum composition rather than at molecular
resolution is unchanged and is argued in that document.

    python scripts/build_ganglioside_library.py --out data/libraries/Ganglioside_Negative_SumComposition.msp

## What is enumerated

A ganglioside is a ceramide carrying a glycan. The neutral mass is therefore the ceramide mass plus
the condensed residue masses of the sugars, and both halves are enumerated independently:

  * ceramide, as a sum composition, from every dihydroxy sphingoid base (C14 to C22, 0 or 1 double
    bond) combined with every N-acyl chain (C10 to C32, 0 to 6 double bonds), collapsed to total
    carbons and total double bonds. `Cer dC:DB` has the formula C(C)H(2C+1-2DB)NO3.
  * glycan, from the Svennerholm series. `-NANA` means every sialic acid is N-acetylneuraminic
    acid; `-NGNA` means one of them is N-glycolylneuraminic acid, the mixed form, which is what the
    previous library encoded for GD3.

⚠ The previous library, enumerated from LipiDex 2, carried GD2-NGNA at a neutral mass 84.02 Da
below the Svennerholm composition, an error inherited from its source. No identification in this
laboratory's data or in the published corpus was a GD2, so nothing that has been reported depends
on it, but the mass written here is the computed one and differs from the old file.
"""
from __future__ import annotations

import argparse
from pathlib import Path

H = 1.0078250319
C = 12.0
N = 14.0030740
O = 15.9949146
PROTON = 1.00727646

# Condensed (water-lost) residue masses of the monosaccharides that build the glycan.
RESIDUE = {"Hex": 162.052824, "HexNAc": 203.079373, "NeuAc": 291.095417, "NeuGc": 307.090331}

# Svennerholm series: how many of each residue the glycan of each class carries.
GLYCAN = {
    "GM3": {"Hex": 2, "Sia": 1},
    "GM2": {"Hex": 2, "HexNAc": 1, "Sia": 1},
    "GM1": {"Hex": 3, "HexNAc": 1, "Sia": 1},
    "GD3": {"Hex": 2, "Sia": 2},
    "GD2": {"Hex": 2, "HexNAc": 1, "Sia": 2},
    "GD1": {"Hex": 3, "HexNAc": 1, "Sia": 2},
    "GT1": {"Hex": 3, "HexNAc": 1, "Sia": 3},
    "GQ1": {"Hex": 3, "HexNAc": 1, "Sia": 4},
}
# The class and sialic-acid forms the library covers. Kept to the seven the previous library had,
# so that results stay comparable; `--all-classes` adds the all-NeuAc forms of GM1, GD2 and GD3.
DEFAULT_FORMS = [("GM3", "NANA"), ("GM3", "NGNA"), ("GM2", "NANA"), ("GM2", "NGNA"),
                 ("GM1", "NGNA"), ("GD3", "NGNA"), ("GD2", "NGNA")]
EXTRA_FORMS = [("GM1", "NANA"), ("GD3", "NANA"), ("GD2", "NANA")]

# The sphingoid base and N-acyl ranges the enumeration spans.
BASE_CARBONS = range(14, 23)
BASE_DOUBLE_BONDS = (0, 1)
ACYL_CARBONS = range(10, 33)
ACYL_DOUBLE_BONDS = range(0, 7)
MIN_CARBONS, MAX_CARBONS = 24, 45

# Measured in this laboratory's own negative-mode spectra; intensities are the mean relative to the
# base peak over the species confirmed independently in positive mode. See docs/GANGLIOSIDES.md.
NEGATIVE_FRAGMENTS = [
    (290.0881, 999.0, "C11H16NO8_Fragment_[]"),   # [NeuAc - H2O - H]-
    (87.0088, 130.0, "C3H3O3_Fragment_[]"),       # sialic acid ring fragment
    (272.0776, 23.0, "C11H14NO7_Fragment_[]"),    # [NeuAc - 2H2O - H]-
]


def ceramide_mass(carbons: int, double_bonds: int) -> float:
    """`Cer dC:DB` = C(C)H(2C+1-2DB)NO3, the amide of a dihydroxy sphingoid base and a fatty acyl."""
    return carbons * C + (2 * carbons + 1 - 2 * double_bonds) * H + N + 3 * O


def glycan_mass(lipid_class: str, sialic: str) -> float:
    parts = GLYCAN[lipid_class]
    mass = parts.get("Hex", 0) * RESIDUE["Hex"] + parts.get("HexNAc", 0) * RESIDUE["HexNAc"]
    n_sia = parts["Sia"]
    if sialic == "NANA":
        return mass + n_sia * RESIDUE["NeuAc"]
    return mass + (n_sia - 1) * RESIDUE["NeuAc"] + RESIDUE["NeuGc"]


def compositions() -> list[tuple[int, int]]:
    """Every ceramide sum composition reachable from a base and an acyl chain, within the range."""
    found = set()
    for bc in BASE_CARBONS:
        for bd in BASE_DOUBLE_BONDS:
            for ac in ACYL_CARBONS:
                for ad in ACYL_DOUBLE_BONDS:
                    carbons, double_bonds = bc + ac, bd + ad
                    if MIN_CARBONS <= carbons <= MAX_CARBONS:
                        found.add((carbons, double_bonds))
    return sorted(found)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--all-classes", action="store_true",
                    help="also write the all-NeuAc forms of GM1, GD2 and GD3")
    args = ap.parse_args()

    forms = DEFAULT_FORMS + (EXTRA_FORMS if args.all_classes else [])
    grid = compositions()
    lines: list[str] = []
    written = 0
    for lipid_class, sialic in forms:
        head = glycan_mass(lipid_class, sialic)
        for carbons, double_bonds in grid:
            neutral = ceramide_mass(carbons, double_bonds) + head
            precursor = neutral - PROTON
            peaks = [(mz, i, a) for mz, i, a in NEGATIVE_FRAGMENTS if precursor - mz > 2.0]
            if not peaks:
                continue
            name = f"{lipid_class}-{sialic} d{carbons}:{double_bonds}"
            lines.append(f"Name: {name} [M-H]-;")
            lines.append(f"MW: {neutral:.4f}")
            lines.append(f"PRECURSORMZ: {precursor:.4f}")
            lines.append(f"Comment: Name={name} [M-H]- Mass={precursor:.4f} OptimalPolarity=true "
                         f"Type=LipidLoop Source=enumerated")
            lines.append(f"Num Peaks: {len(peaks)}")
            for mz, intensity, annotation in peaks:
                lines.append(f'{mz:.4f} {intensity:.0f} "{annotation}"')
            lines.append("")
            written += 1
    Path(args.out).write_text("\n".join(lines))
    print(f"{written} entries written to {args.out}")
    print(f"  {len(forms)} class forms x {len(grid)} ceramide sum compositions "
          f"(C{MIN_CARBONS}-{MAX_CARBONS}, from bases C{BASE_CARBONS.start}-{BASE_CARBONS.stop-1} "
          f"and acyls C{ACYL_CARBONS.start}-{ACYL_CARBONS.stop-1})")
    print("  enumerated from the Svennerholm structures and computed masses; no third-party library read")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
