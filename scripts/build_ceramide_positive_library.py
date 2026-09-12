"""Generate positive-mode entries for the ceramide and glucosylceramide classes that have none.

Of the eight free-ceramide subclasses LipiDex enumerates, **only `Cer[NS]` has positive-mode
entries**. `Cer[AS]`, `Cer[BS]`, `Cer[NDS]`, `Cer[ADS]`, `Cer[BDS]`, `Cer[AP]` and `Cer[NP]` exist
only as `[M-H]-`, so in positive mode they cannot be identified at all — whatever their abundance.

In skin that is the largest single library gap. Screening the unidentified features of the skin
positive run for homologous series (`scripts/find_unclaimed_series.py`) puts a 16-member ladder at
the top, m/z 598.58 to 808.81, **total area 6.2e9**, present in up to 61 of 63 files. Matching its
members against the neutral masses of the negative-mode library identifies them as `Cer[ADS]`,
`Cer[AS]` and `Cer[NP]` — the classes with no positive entries. Their MS2 spectra carry the
sphingoid-base fragments to prove it.

## What positive mode adds, and it is not just coverage

Negative mode cannot separate `Cer[ADS]` from `Cer[NP]`: `Cer[ADS] d18:0_24:0` and
`Cer[NP] t18:0_24:0` are **the same formula and the same exact mass**, and negative-mode spectra
carry nothing that distinguishes them (see docs/FATTY_ACIDS.md's companion note in
results/skin_cer_ads_check.md, where the Cer[ADS]/Cer[AS] call comes down to 1% of the score).

Positive mode does separate them, because the **sphingoid base fragment differs**:

    d18:0 (sphinganine, 2 OH)      302.3054  284.2948  266.2842
    t18:0 (phytosphingosine, 3 OH) 318.3003  300.2897  282.2791

Three masses apart. That is a real gain over the negative-mode assignment, and the reason to build
this rather than treat positive mode as redundant coverage.

**It does not separate everything.** `Cer[AS]` and `Cer[BS]` — alpha- versus beta-hydroxy acyl —
share a formula and a base, so they stay degenerate here exactly as they are in negative mode.

## The fragments

LipiDex's own positive-mode `Cer[NS]` entries carry four peaks — the precursor water loss and
three sphingoid-base ions — and that scheme is reused rather than invented, since it is already
validated to exactness against the reference implementation. Annotations keep LipiDex's format so
the purity calculation can resolve the base the same way it resolves acyl chains elsewhere.

**Two-hydroxyl bases (`d`) keep LipiDex's intensities unchanged**, because this dataset contains
none of them to measure and inventing numbers would be worse than inheriting them.

**Three-hydroxyl bases (`t`) were measured here**, on the very series this library exists to claim
— every member turned out to be a phytoceramide with a long-chain base, t17:0 to t22:0. They get a
**third water loss**, which a two-hydroxyl base cannot produce, and it arrives nearly as intense as
the second. That peak is the discriminator: it is what separates `Cer[NP] t18:0_24:0` from
`Cer[ADS] d18:0_24:0` at identical exact mass.

Bases are parsed from the existing entry names, so the enumeration and the naming are identical to
the negative-mode library — nothing new is invented, only the polarity is added.

    python scripts/build_ceramide_positive_library.py \\
        --out data/libraries/Ceramides_Positive.msp
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.msp import parse_msp          # noqa: E402
from lipidloop.search import DEFAULT_LIBRARIES  # noqa: E402

PROTON = 1.00727646
H2O = 18.0105646
CH2O = 30.0105646
CARBON, HYDROGEN, NITROGEN, OXYGEN = 12.0, 1.0078250319, 14.0030740052, 15.9949146221

# The subclasses with no positive-mode entries anywhere in the default set.
MISSING = ("Cer[AS]", "Cer[BS]", "Cer[NDS]", "Cer[ADS]", "Cer[BDS]", "Cer[AP]", "Cer[NP]",
           # Glucosylceramides are negative-only too. GlcCer[NS] and GlcCer[NDS] carry `d` bases
           # only, so the phyto restriction below removes them anyway — and it should, because
           # `HexCer[NS]` already has 600 positive entries covering the same molecules under a
           # different class name. GlcCer[AP] is the one that is genuinely absent: 215 entries,
           # every one a phyto base, with no HexCer[AP] anywhere to duplicate.
           "GlcCer[NS]", "GlcCer[NDS]", "GlcCer[AP]")

# Sugar losses, taken from LipiDex's own positive `HexCer[NS]` entries rather than assumed: the
# hexose leaves as the free sugar and as the sugar plus water, and both peaks are written.
HEXOSE = 180.0633881
HEXOSE_WATER = HEXOSE + H2O
INTENSITY_HEXOSE_LOSS = 100
_HEXOSYL = ("GlcCer", "GalCer", "HexCer")

# Hydroxyls on the sphingoid base, which is what sets the fragment masses.
BASE_OXYGENS = {"d": 2, "t": 3, "m": 1}

# Two intensity schemes, and the difference between them is the honest part.
#
# PHYTO bases (t, three hydroxyls) were measured here, pooled over six files, on the very series
# this library exists to claim — every member of it turned out to be a phytoceramide with a
# long-chain base, t17:0 to t22:0. A third water loss is available to them and is nearly as
# intense as the second, which is itself diagnostic: a two-hydroxyl base cannot produce it.
#
#     m/z 668.66 (t18:0)   -2H2O 930   -3H2O 898   -H2O 615
#     m/z 724.72 (t22:0)   -3H2O 942   -2H2O 929   -H2O 547
#     m/z 696.69 (t20:0)   -3H2O 507   -H2O  426   -2H2O 413
#
# DIHYDRO and sphingosine bases (d) were not measured here — the skin series contains none — so
# they keep LipiDex's own scheme unchanged, which is already validated to exactness against the
# reference implementation. Inventing intensities for them would be worse than inheriting them.
PHYTO_INTENSITIES = {"-H2O": 600, "-2H2O": 999, "-3H2O": 950, "precursor": 400}
DIHYDRO_INTENSITIES = {"-H2O": 100, "-2H2O": 999, "-H2O-CH2O": 100, "precursor": 50}

# Negative adducts these classes are enumerated under, as the mass added to the neutral. Classes
# differ: free ceramides are [M-H]-, GlcCer[AP] only ever [M+FA-H]-.
FORMIC = 46.00548
ACETIC = 60.02113
NEGATIVE_ADDUCTS = {"[M-H]": 0.0, "[M+FA-H]": FORMIC, "[M+Ac-H]": ACETIC}

PRECURSOR_CUTOFF = 2.0          # the library reader drops peaks closer than this to the precursor
_BASE = re.compile(r"([dtm])(\d+):(\d+)")


def base_mass(prefix: str, carbons: int, double_bonds: int) -> float:
    """Neutral sphingoid base: C_n H_(2n+3-2db) N O_(1..3)."""
    oxygens = BASE_OXYGENS[prefix]
    hydrogens = 2 * carbons + 3 - 2 * double_bonds
    return (carbons * CARBON + hydrogens * HYDROGEN + NITROGEN + oxygens * OXYGEN)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--classes", nargs="*", default=list(MISSING))
    ap.add_argument("--all-bases", action="store_true",
                    help="also emit two-hydroxyl (d) bases — measured to cost 8 correct calls")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]

    # Fingerprint every positive-mode entry that already exists, as (precursor, fragment masses).
    # A new entry with the same fingerprint carries no information the library does not already
    # have, so it cannot help — it can only take matches away from an entry that is already
    # validated. That is not hypothetical: the first build of this library cost 11 correct
    # `Cer[NS]` identifications, each re-labelled `Cer[NDS]`, `Cer[ADS]` or `Cer[NP]` at the same
    # sum composition. LipiDex enumerates these classes over a shared base/acyl grid, so a class
    # label does not constrain the sphingoid base the way the notation suggests, and entries in
    # different classes can be identical in both precursor and fragments. Where that happens the
    # winner is decided by library load order, which is to say by nothing.
    # Never read the file about to be written. Once this library is in DEFAULT_LIBRARIES a naive
    # scan finds its own previous output and calls every entry degenerate with itself — which it
    # silently did, cutting 446 entries to 98 with no error anywhere.
    out_path = Path(args.out).resolve()
    existing: set = set()
    for name in DEFAULT_LIBRARIES:
        path = root / "data/libraries" / f"{name}.msp"
        if not path.exists() or path.resolve() == out_path:
            continue
        for spectrum in parse_msp(path):
            text = spectrum.name.strip()
            if "]" not in text or "-" in text[text.rindex("]") + 1:]:
                continue
            existing.add((round(spectrum.precursor_mz, 3),
                          tuple(sorted(round(m, 3) for m in spectrum.mz))))

    seen: dict[tuple, float] = {}
    for name in DEFAULT_LIBRARIES:
        path = root / "data/libraries" / f"{name}.msp"
        if not path.exists():
            continue
        for spectrum in parse_msp(path):
            if spectrum.lipid_class not in args.classes:
                continue
            text = spectrum.name.strip()
            adduct = text[text.rindex("["):text.rindex("]") + 1]
            # One entry per molecule, recovering the neutral from whichever negative adduct the
            # class happens to be enumerated under. Assuming [M-H]- silently dropped every
            # `GlcCer[AP]`, which exists only as [M+FA-H]-.
            offset = NEGATIVE_ADDUCTS.get(adduct)
            if offset is None:
                continue
            seen.setdefault((spectrum.lipid_class, text.split(" [")[0]),
                            spectrum.precursor_mz + PROTON - offset)

    lines, written, skipped, degenerate, skipped_base = [], 0, 0, 0, 0
    for (lipid_class, molecule), neutral in sorted(seen.items()):
        match = _BASE.search(molecule)
        if match is None:
            skipped += 1
            continue
        prefix, carbons, double_bonds = match.group(1), int(match.group(2)), int(match.group(3))
        # Only three-hydroxyl bases by default, and the reason is measured rather than assumed.
        # A `d`-base entry in these classes carries the same fragments as an existing positive
        # `Cer[NS]` entry with the same base — 85% of them were dropped as exactly degenerate
        # above — and the survivors are not survivors on merit: they win where LipiDex's positive
        # `Cer[NS]` set has an enumeration gap (1,200 positive entries against 2,406 negative), so
        # they are filling holes in a grid rather than matching evidence. Measured cost of
        # emitting them: 8 correct `Cer[NS]` identifications re-labelled `Cer[NDS]`, and recovery
        # against the reference down from 98.5% to 97.1%.
        #
        # `t` bases are different in kind. Nothing else in the library can produce a third water
        # loss, so those entries add a discriminator instead of competing for one.
        if prefix != "t" and not args.all_bases:
            skipped_base += 1
            continue
        base = base_mass(prefix, carbons, double_bonds)
        precursor = neutral + PROTON
        label = f"{prefix}{carbons}:{double_bonds}"

        hexosyl = lipid_class.startswith(_HEXOSYL)
        sugar_peaks = []
        if hexosyl:
            # The base fragments sit on the ceramide left behind, so the sugar loss has to happen
            # before them — these two peaks are what say the molecule carried a hexose at all.
            sugar_peaks = [
                (precursor - HEXOSE, INTENSITY_HEXOSE_LOSS, "C-6H-12O-6_Neutral Loss_[]"),
                (precursor - HEXOSE_WATER, INTENSITY_HEXOSE_LOSS, "C-6H-14O-7_Neutral Loss_[]"),
            ]

        if prefix == "t":
            weights = PHYTO_INTENSITIES
            peaks = [
                (precursor - H2O, weights["precursor"], "H-2O-1_Neutral Loss_[]"),
                (base + PROTON - H2O, weights["-H2O"],
                 f"C3H9O2N1_Sphingoid Fragment_[{label}]"),
                (base + PROTON - 2 * H2O, weights["-2H2O"],
                 f"C3H7O1N1_Sphingoid Fragment_[{label}]"),
                # Only a three-hydroxyl base can lose a third water. This peak is what separates
                # a phytoceramide from the dihydroceramide of identical formula and exact mass.
                (base + PROTON - 3 * H2O, weights["-3H2O"],
                 f"C3H5N1_Sphingoid Fragment_[{label}]"),
            ]
        else:
            weights = DIHYDRO_INTENSITIES
            peaks = [
                (precursor - H2O, weights["precursor"], "H-2O-1_Neutral Loss_[]"),
                (base + PROTON - H2O, weights["-H2O"],
                 f"C3H7O1N1_Sphingoid Fragment_[{label}]"),
                (base + PROTON - 2 * H2O, weights["-2H2O"],
                 f"C3H5N1_Sphingoid Fragment_[{label}]"),
                (base + PROTON - H2O - CH2O, weights["-H2O-CH2O"],
                 f"C2H5N1_Sphingoid Fragment_[{label}]"),
            ]
        peaks += sugar_peaks
        peaks = [p for p in peaks if precursor - p[0] > PRECURSOR_CUTOFF and p[0] > 61.0]
        if len(peaks) < 2:
            skipped += 1
            continue
        peaks.sort()

        fingerprint = (round(precursor, 3), tuple(sorted(round(m, 3) for m, _, _ in peaks)))
        if fingerprint in existing:
            degenerate += 1
            continue
        existing.add(fingerprint)

        lines.append(f"Name: {molecule} [M+H]+;")
        lines.append(f"MW: {neutral:.4f}")
        lines.append(f"PRECURSORMZ: {precursor:.4f}")
        lines.append(f"Comment: Name={molecule} [M+H]+ Mass={precursor:.4f} "
                     f"OptimalPolarity=false Type=LipiDex")
        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, intensity, annotation in peaks:
            lines.append(f'{mz:.4f} {intensity} "{annotation}"')
        lines.append("")
        written += 1

    Path(args.out).write_text("\n".join(lines))
    print(f"{written} positive-mode ceramide entries written to {args.out}")
    if skipped:
        print(f"  {skipped} skipped: no parsable sphingoid base, or too few usable fragments")
    if skipped_base:
        print(f"  {skipped_base} two-hydroxyl (d) bases not emitted — they duplicate the existing "
              f"positive Cer[NS] entries and cost 8 correct calls when tried; --all-bases "
              f"overrides")
    if degenerate:
        print(f"  {degenerate} dropped as degenerate: identical precursor AND fragments to an "
              f"entry the library already has, so they could only steal matches")
    counts: dict[str, int] = {}
    for (lipid_class, _) in seen:
        counts[lipid_class] = counts.get(lipid_class, 0) + 1
    for lipid_class, n in sorted(counts.items()):
        print(f"    {lipid_class:10} {n}")
    print("  OptimalPolarity=false: these are the polarity of last resort for most of these "
          "classes, and the flag governs whether purity is scored at all")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
