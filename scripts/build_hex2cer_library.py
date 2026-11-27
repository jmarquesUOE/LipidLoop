"""Build the negative-mode dihexosylceramide (Hex2Cer) library, reading no third-party file.

    python scripts/build_hex2cer_library.py --out data/libraries/Hex2Cer_Negative_Formate.msp

## Why this library exists

Hex2Cer is the direct product of B4GALT5, the lactosylceramide synthase, and the branch point of
the whole glycosphingolipid tree: everything below it — globosides, and the ganglioside series this
repository already covers — is built on it. Nothing in the shipped libraries names it. LipidBlast
and the LipiDex HCD libraries carry monohexosylceramide only, and the ganglioside library starts at
GM4/GM3, one step further down. So a study could measure the substrate and the products two steps
on, and pass silently over the enzyme's own product.

It was found by the gap it left. Four unnamed features in a facility ceramide-inhibitor study sat
within 3 ppm of Hex2Cer formate adducts, with MS2 acquired and no library to match it against.

## The adduct

**Entries are written as `[M+HCOO]-`, not `[M-H]-`, because that is how the molecule is seen.** A
neutral glycosphingolipid has no acidic proton to lose; in a formate mobile phase it is carried as
the formate adduct. All four species confirmed here are formate adducts, and each fragments first
by losing formic acid to give `[M-H]-` as the base peak. Writing `[M-H]-` entries would put the
precursor 45 Da from anything the instrument selected.

⚠ This ties the library to the mobile phase, exactly as
`reference-negative-adduct-misassignment` warns. It is correct for this facility, which runs
ammonium formate. An acetate method would need `[M+CH3COO]-` at +59.0133 instead.

## What the fragmentation supports, and what it does not

Averaged over the spectra of the four confirmed species, three ions carry the composition:

    [M-H]-          loss of formic acid from the precursor — the base peak
    [M-H-Hex]-      loss of one hexose, 162.0528
    [M-H-2Hex]-     loss of both, landing on the ceramide anion

The second loss is the whole identification. One hexose loss is a monohexosylceramide; **two
sequential losses of 162.0528 arriving at the calculated ceramide anion is what makes this
dihexosyl**, and it was observed directly (d34:1: 860.61 -> 698.55 -> 536.50, against a calculated
Cer d34:1 anion of 536.505).

Intensities are medians relative to the `[M-H]-` base peak, measured over every spectrum of the four
confirmed species in all 39 sample injections of this study (d34:1 and d42:2 contribute 39 spectra
each, d40:1 21, d42:1 16). They are weights for the dot product, not a claim about absolute
response.

⚠ The two hexose losses are not always present: `-Hex` is absent from the d40:1 and d42:1 medians
and `-2Hex` from d42:1, the two species with the fewest and weakest spectra. The library lists both
regardless — a missing peak costs a little dot product, whereas omitting them would discard the only
ions that distinguish this class from a monohexosylceramide.

The rest of the spectrum is the hexose series — 179.0561, 161.0455, 143.0350, 119.0350, 113.0244,
101.0244, 89.0244, 71.0139 — which is intense but says only "a sugar was here". Those ions appear
in any glycosphingolipid spectrum, so weighting them heavily would let a Hex2Cer entry score
against a HexCer or a ganglioside. Only the four most reliable are written, at their measured
intensities, as corroboration rather than evidence.

⚠ **Sum composition, no subclass tag, and the name is `Hex2Cer` rather than `LacCer`.** Two things
are being withheld. The base and acyl chains are not resolved, for the reason argued for the
gangliosides in the developer notes (GANGLIOSIDES, not distributed), so `[NS]` and `[NDS]` cannot be told apart and no tag is
written. And lactosylceramide is only one dihexosylceramide: Gal-Glc-Cer, Gal-Gal-Cer and
Glc-Glc-Cer are isomers that reversed-phase chromatography does not separate and whose spectra here
are identical. `LacCer` would assert a linkage the evidence cannot support, so the name states the
composition. `nomenclature.CANONICAL` already maps `laccer` and `lactosylceramide` onto `Hex2Cer`.

⚠ **No decoy can be built from this library.** `make_decoy_library.py` works by shifting
chain-specific fragments, and every ion here is head-group or precursor derived, annotated `[]`.
Hex2Cer identifications therefore carry the same caveat the gangliosides do: their false-discovery
rate is a ceiling inherited from the rest of the run, not a measurement.
"""
from __future__ import annotations

import argparse
from pathlib import Path

H = 1.0078250319
C = 12.0
N = 14.0030740
O = 15.9949146
PROTON = 1.00727646
ELECTRON = 0.00054858

#: Condensed (water-lost) hexose residue. The same constant the ganglioside builder uses.
HEX = 162.052824
#: [M+HCOO]- = M + HCOOH - H+. Formic acid is 46.005479.
FORMATE_ADDUCT = 46.005479 - PROTON

# The sphingoid base and N-acyl ranges the enumeration spans, and the sum-composition window they
# collapse to. Identical to the ganglioside library, so the two cover the same species grid and a
# Hex2Cer and its GM3 can always be compared.
BASE_CARBONS = range(14, 23)
BASE_DOUBLE_BONDS = (0, 1)
ACYL_CARBONS = range(10, 33)
ACYL_DOUBLE_BONDS = range(0, 7)
MIN_CARBONS, MAX_CARBONS = 24, 45

# Hexose ions, measured in this laboratory's own negative-mode spectra. Intensities are medians
# relative to the [M-H]- base peak at 999, over the spectra of the four confirmed species. Kept to
# the four that were present in every species: the series is diagnostic of a sugar, not of which
# sugar or how many, so it corroborates and must not dominate.
HEXOSE_FRAGMENTS: list[tuple[float, float, str]] = [
    (179.0561, 142.0, "C6H11O6_Fragment_[]"),   # hexose - H
    (161.0455, 90.0, "C6H9O5_Fragment_[]"),     # - H2O
    (101.0244, 243.0, "C4H5O3_Fragment_[]"),    # cross-ring
    (89.0244, 304.0, "C3H5O3_Fragment_[]"),     # cross-ring
]
# Composition-carrying losses, relative to the same base peak. An eighth of it and less, but these
# are the ions that say dihexosyl, so the scoring must still see them.
LOSS_ONE_HEX = 114.0
LOSS_TWO_HEX = 144.0
BASE_PEAK = 999.0


def ceramide_mass(carbons: int, double_bonds: int) -> float:
    """`Cer dC:DB` = C(C)H(2C+1-2DB)NO3, the amide of a dihydroxy sphingoid base and a fatty acyl."""
    return carbons * C + (2 * carbons + 1 - 2 * double_bonds) * H + N + 3 * O


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
    ap.add_argument("--hexoses", type=int, default=2, choices=(2, 3),
                    help="2 for Hex2Cer (measured), 3 for Hex3Cer — Gb3 and the globo/lacto "
                         "trihexosides, whose intensity model is TRANSFERRED, not measured")
    args = ap.parse_args()

    grid = compositions()
    lines: list[str] = []
    written = 0
    n_hex = args.hexoses
    head = n_hex * HEX
    for carbons, double_bonds in grid:
        neutral = ceramide_mass(carbons, double_bonds) + head
        precursor = neutral + FORMATE_ADDUCT
        deprotonated = neutral - PROTON
        peaks = [(deprotonated, BASE_PEAK, "Formate_Loss_Fragment_[]")]
        for k in range(1, n_hex + 1):
            last = k == n_hex
            peaks.append((deprotonated - k * HEX,
                          LOSS_TWO_HEX if last else LOSS_ONE_HEX,
                          "Ceramide_Fragment_[]" if last else f"Hexose_Loss_{k}_Fragment_[]"))
        peaks += HEXOSE_FRAGMENTS
        # Never write a product ion at or above the precursor, and never one below the range the
        # method records.
        peaks = [(mz, i, a) for mz, i, a in peaks if 2.0 < mz < precursor - 2.0 and i > 0]
        if len(peaks) < 3:
            continue
        name = f"Hex{n_hex}Cer d{carbons}:{double_bonds}"
        lines.append(f"Name: {name} [M+HCOO]-;")
        lines.append(f"MW: {neutral:.4f}")
        lines.append(f"PRECURSORMZ: {precursor:.4f}")
        lines.append(f"Comment: Name={name} [M+HCOO]- Mass={precursor:.4f} OptimalPolarity=true "
                     f"Type=LipidLoop Source=enumerated")
        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, intensity, annotation in peaks:
            lines.append(f'{mz:.4f} {intensity:.0f} "{annotation}"')
        lines.append("")
        written += 1
    Path(args.out).write_text("\n".join(lines))
    print(f"{written} entries written to {args.out}")
    print(f"  Hex{n_hex}Cer: {len(grid)} ceramide sum compositions (C{MIN_CARBONS}-{MAX_CARBONS}, from bases "
          f"C{BASE_CARBONS.start}-{BASE_CARBONS.stop - 1} and acyls "
          f"C{ACYL_CARBONS.start}-{ACYL_CARBONS.stop - 1}), written as [M+HCOO]-")
    print("  masses computed from the elemental formula; fragments measured here; "
          "no third-party library read")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
