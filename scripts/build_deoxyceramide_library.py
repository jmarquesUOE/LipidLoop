"""Generate positive-mode entries for 1-deoxyceramides, the class no library we search contains.

Zero entries across all 308,624 names in the five searched libraries, and no `m`-prefix sphingoid
base anywhere (the developer notes (DEOXYSPHINGOLIPIDS, not distributed)). They are absent by construction, not by abundance —
and they are present in our data: seven species confirmed by MS2 in an in-house Q Exactive Plus study (Thermo Q
Exactive Plus) and `MTBLS5163` (Bruker Q-TOF Compact), two vendors, two labs.

## Why the grid is saturated-only

A 1-deoxyceramide with *k* double bonds has exactly the formula of the canonical ceramide of the
same carbon count with *k-1* after that ceramide loses water, and ceramides lose water freely in
positive mode:

    deoxy 42:1 [M+H]+      = 634.6497
    Cer   42:0 [M+H-H2O]+  = 634.6497      same ion, C42H84NO2+, 0.0 ppm

An unsaturated entry therefore cannot win a match on evidence — it can only take spectra away from
a canonical ceramide that is correctly assigned. The saturated species have no such partner (it
would need a canonical at -1 double bonds), which is why they are the ones that could be confirmed
and the only ones written here. This is the same discipline as the withdrawn `Cer[EOS]` and
`FAHFA` libraries: an entry that cannot be distinguished from something already in the library is
a liability, not coverage.

## The intensities are MEASURED, and that was the blocker

`build_ceramide_positive_library.py` refuses to invent intensities, and for `m` bases there were
none to inherit — a one-hydroxyl base has a genuinely different water-loss ladder, since only one
loss is available to it. Measured here over **154 confirmed spectra** pooled across both studies,
normalised to base peak 999:

    base - H2O        999   present in 154/154     the only water loss a single hydroxyl allows
    N-acyl amide      244   present in 150/154     names the acyl chain
    base [M+H]+        38   present in 108/154
    precursor - H2O    33   present in 104/154
    precursor           0   present in   0/154     never survives

The base ion is `C18H38N+` 268.2999 for `m18:0`. Note it is NOT the 266.2842 that
the developer notes (DEOXYSPHINGOLIPIDS, not distributed) warns about — that isobar is the *unsaturated* deoxy base against the
doubly-dehydrated dihydro base. A saturated `m` base sits 2.0157 away from it, which is what makes
these seven species callable at all.

## Annotation, and why exactly one peak carries the base

`purity.py` picks the most intense annotation type that is neither `_Fragment` nor `_Neutral Loss`,
then requires the number of matched chain fragments of that type to EQUAL the number of chains the
name declares. Writing two `Sphingoid Fragment` peaks for one base makes `found=2` against
`len(lipid_fas)=1`, `_fa_intensities` returns None and **Purity comes back 0** — which is exactly
why the phytoceramide entries score 0 (the developer notes (CHAIN_RESOLUTION, not distributed)). So the base-diagnostic ion is
annotated once, and the intact base is written as a generic `_Fragment_[]` that the purity
calculation ignores. Purity on these entries then means what it should: the base-diagnostic ion
was observed.

The `m` sphingoid bases must exist in `data/lipidex_src/FattyAcids.csv` under base `Sphingoid`
for that to resolve. That file is **gitignored**, like everything under `data/`, so no commit
records the requirement — the build refuses to run without them and `--add-fa-bases` writes them:

    m16:0,Sphingoid,C13H27,true      m19:0,Sphingoid,C16H33,true
    m17:0,Sphingoid,C14H29,true      m20:0,Sphingoid,C17H35,true
    m18:0,Sphingoid,C15H31,true

    python scripts/build_deoxyceramide_library.py \\
        --out data/libraries/DeoxyCeramides_Positive.msp

⚠ Not in DEFAULT_LIBRARIES. Adding a library is not free — measure displacement against the decoys
first (the developer notes (LIBRARY_VERSIONS, not distributed)), as `LipidBlast` taking 278 of 341 mis-attributed matches and the
withdrawn `FAHFA` set both showed.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.msp import parse_msp             # noqa: E402
from lipidloop.search import DEFAULT_LIBRARIES  # noqa: E402

PROTON = 1.00727646
H2O = 18.0105646
CARBON, HYDROGEN, NITROGEN, OXYGEN = 12.0, 1.0078250319, 14.0030740052, 15.9949146221

CLASS = "DeoxyCer"

# Bases: measured range. m18:0 carries most of the signal, m20:0 is the base of the 44:0 species
# in in-house data, and MTBLS5163's odd-chain 41:0 shows a faint m19:0 beside the m18:0 — so the odd bases
# are in on evidence, not on symmetry.
BASES = (16, 17, 18, 19, 20)
# N-acyl: 16:0, 22:0, 23:0, 24:0 and 26:0 are the ones observed. The grid runs a little either
# side of that; it is not extended further because every extra entry is another chance to take a
# match from something correct.
ACYLS = tuple(range(14, 27))

# Measured over 154 confirmed spectra, an in-house Q Exactive Plus study + MTBLS5163. See the module docstring.
INTENSITIES = {"base-H2O": 999, "acylamide": 245, "base": 40, "precursor-H2O": 35}

PRECURSOR_CUTOFF = 2.0      # the library reader drops peaks closer than this to the precursor
MIN_MS2_MASS = 61.0


def base_neutral(carbons: int) -> float:
    """1-deoxy sphingoid base, one hydroxyl: Cn H(2n+3) N O."""
    return carbons * CARBON + (2 * carbons + 3) * HYDROGEN + NITROGEN + OXYGEN


def ceramide_neutral(base_c: int, acyl_c: int) -> float:
    """base + fatty acid - H2O.  C(a+b) H(2a+2b+1) N O2 for the saturated case."""
    n = base_c + acyl_c
    return n * CARBON + (2 * n + 1) * HYDROGEN + NITROGEN + 2 * OXYGEN


def acylamide(acyl_c: int) -> float:
    """The protonated N-acyl amide: neutral fatty amide Cn H(2n+1) N O, plus a proton.

    Measured at 340.3574 for 22:0 and 368.3887 for 24:0 in both studies, where it is the second
    most intense ion and the only one that names the acyl chain.

    ⚠ The neutral amide carries 2n+1 hydrogens, not 2n. Getting that wrong puts every acyl ion
    1.0078 low — close enough to look right beside a measured value and far enough that not one
    of them would ever match at 0.01 Da.
    """
    return acyl_c * CARBON + (2 * acyl_c + 1) * HYDROGEN + NITROGEN + OXYGEN + PROTON


def fa_db_path(root: Path) -> Path:
    return root / "data/lipidex_src/FattyAcids.csv"


def fa_db_has_bases(root: Path, bases) -> list[str]:
    """Which `m` bases are missing from the fatty-acid database, as `Sphingoid` entries."""
    have = set()
    with fa_db_path(root).open(newline="") as fh:
        for row in csv.DictReader(fh):
            if (row.get("Base") or "").strip() == "Sphingoid":
                have.add((row.get("Name") or "").strip())
    return [f"m{b}:0" for b in bases if f"m{b}:0" not in have]


def add_fa_bases(root: Path, bases) -> list[str]:
    """Append the missing `m` bases to the fatty-acid database.

    ⚠ `data/**` is gitignored, so this file is local state that no commit records. Leaving the
    requirement as a sentence in a docstring means a fresh clone builds a library whose every
    entry scores Purity 0 — so the fix is executable rather than written down.

    The `Sphingoid` formula column holds the alkyl tail only (`d18:0` is `C15H31`, carbons minus
    the three-carbon head). A 1-deoxy base has the same tail; the head oxygens differ, and the
    fragment annotation is what carries that.
    """
    path = fa_db_path(root)
    missing = fa_db_has_bases(root, bases)
    if not missing:
        return []
    text = path.read_text()
    if not text.endswith("\n"):
        text += "\n"
    for name in missing:
        carbons = int(name[1:].split(":")[0])
        tail = carbons - 3
        text += f"{name},Sphingoid,C{tail}H{2 * tail + 1},true\n"
    path.write_text(text)
    return missing


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--bases", nargs="*", type=int, default=list(BASES))
    ap.add_argument("--acyls", nargs="*", type=int, default=list(ACYLS))
    ap.add_argument("--unsaturated", action="store_true",
                    help="also emit entries with double bonds — every one of them is an exact "
                         "isobar of a canonical ceramide water-loss ion, so this is for "
                         "measuring the damage, not for shipping")
    ap.add_argument("--add-fa-bases", action="store_true",
                    help="append any missing `m` sphingoid bases to "
                         "data/lipidex_src/FattyAcids.csv (that file is gitignored, so this is "
                         "the only record of what a fresh clone needs)")
    ap.add_argument("--no-check-fa-db", action="store_true")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    out_path = Path(args.out).resolve()

    if args.add_fa_bases:
        added = add_fa_bases(root, args.bases)
        print(f"added to data/lipidex_src/FattyAcids.csv: {', '.join(added) or 'nothing missing'}")

    if not args.no_check_fa_db:
        missing = fa_db_has_bases(root, args.bases)
        if missing:
            print(f"⚠ {len(missing)} base(s) absent from data/lipidex_src/FattyAcids.csv under "
                  f"`Sphingoid`: {', '.join(missing)}")
            print("  Purity will be 0 on every entry using them — the name would be reported "
                  "chain-resolved on no chain evidence, which is the failure "
                  "the developer notes (CHAIN_RESOLUTION, not distributed) is about. Re-run with --add-fa-bases, or pass "
                  "--no-check-fa-db to build anyway.")
            return 1

    # Same guard as build_ceramide_positive_library.py: an entry whose precursor AND fragments
    # already exist carries nothing new and can only take matches away. Never read the file being
    # written, or the library finds its own previous output and calls everything degenerate.
    existing: set = set()
    for name in DEFAULT_LIBRARIES:
        path = root / "data/libraries" / f"{name}.msp"
        if not path.exists() or path.resolve() == out_path:
            continue
        for spectrum in parse_msp(path):
            existing.add((round(spectrum.precursor_mz, 3),
                          tuple(sorted(round(m, 3) for m in spectrum.mz))))

    lines, written, degenerate, thin = [], 0, 0, 0
    combos = [(b, a, 0, 0) for b in args.bases for a in args.acyls]
    if args.unsaturated:
        combos += [(b, a, u, v) for b in args.bases for a in args.acyls
                   for u, v in ((0, 1), (1, 0)) ]

    for base_c, acyl_c, base_db, acyl_db in sorted(combos):
        if base_db or acyl_db:      # unsaturated: shift both the neutral and the base ion
            n = base_c + acyl_c
            neutral = (n * CARBON + (2 * n + 1 - 2 * (base_db + acyl_db)) * HYDROGEN
                       + NITROGEN + 2 * OXYGEN)
        else:
            neutral = ceramide_neutral(base_c, acyl_c)
        precursor = neutral + PROTON
        base_label = f"m{base_c}:{base_db}"
        acyl_label = f"{acyl_c}:{acyl_db}"
        molecule = f"{CLASS} {base_label}_{acyl_label}"

        base_ion = (base_neutral(base_c) - 2 * base_db * HYDROGEN) + PROTON
        peaks = [
            # The one water loss a single hydroxyl allows. Exactly one peak carries the base
            # chain, or the purity calculation's found-vs-declared test fails and returns 0.
            (base_ion - H2O, INTENSITIES["base-H2O"],
             f"C3H5N1_Sphingoid Fragment_[{base_label}]"),
            # Names the acyl chain. Prefix is the delta from the fatty-acid database's own
            # formula (24:0 is C24H47O2 there; the amide ion is C24H50NO+), the same convention
            # as LipiDex's `O-1_Alkyl Fragment_[10:0]`.
            (acylamide(acyl_c) - 2 * acyl_db * HYDROGEN, INTENSITIES["acylamide"],
             f"H3N1O-1_Alkyl Fragment_[{acyl_label}]"),
            # Generic on purpose: a second Sphingoid-typed peak would make the purity
            # calculation count two fragments for one declared chain and return 0.
            (base_ion, INTENSITIES["base"], "C3H7O1N1_Fragment_[]"),
            (precursor - H2O, INTENSITIES["precursor-H2O"], "H-2O-1_Neutral Loss_[]"),
        ]
        peaks = [p for p in peaks
                 if precursor - p[0] > PRECURSOR_CUTOFF and p[0] > MIN_MS2_MASS]
        if len(peaks) < 2:
            thin += 1
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
        # OptimalPolarity=true: these have no acidic group and give nothing in negative mode, so
        # positive is not a fallback for them. The flag governs whether purity is scored at all.
        lines.append(f"Comment: Name={molecule} [M+H]+ Mass={precursor:.4f} "
                     f"OptimalPolarity=true Type=LipiDex")
        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, intensity, annotation in peaks:
            lines.append(f'{mz:.4f} {intensity} "{annotation}"')
        lines.append("")
        written += 1

    Path(args.out).write_text("\n".join(lines))
    print(f"{written} 1-deoxyceramide entries written to {args.out}")
    print(f"  bases m{min(args.bases)}:0-m{max(args.bases)}:0, "
          f"N-acyl {min(args.acyls)}:0-{max(args.acyls)}:0, saturated"
          + (", plus unsaturated (isobaric with canonical water-loss ions)"
             if args.unsaturated else ""))
    if degenerate:
        print(f"  {degenerate} dropped as degenerate with an entry the library already has")
    if thin:
        print(f"  {thin} dropped for too few usable fragments")
    print("  intensities measured on 154 confirmed spectra, an in-house Q Exactive Plus study + MTBLS5163")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
