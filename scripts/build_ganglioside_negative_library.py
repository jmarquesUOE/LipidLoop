"""Build a negative-mode ganglioside library, at the resolution the fragmentation supports.

Neither LipiDex 1 nor LipiDex 2 ships a negative-mode ganglioside library — LipiDex 2's is 4,305
entries, every one `[M+H]+` — so the negative-mode data cannot be searched for them at all. This
generates one from the species LipiDex 2 already enumerates in positive mode.

**It is deliberately built at sum composition, and that is the whole design decision.**

Looking at what these molecules actually produce in negative mode on this instrument, the
richest GM3 spectrum in the dataset contains three ions worth the name:

    290.0881   [NeuAc - H2O - H]-      dominant
    272.0776   [NeuAc - 2H2O - H]-     weak
     87.0088   C3H3O3-                 weak

and nothing else. The glycosidic Y-ions are absent, the ceramide fragment is absent, and the
fatty acyl anion is absent. Every one of those would be needed to say which chains the molecule
carries, so the spectrum says "this is sialylated" and stops.

A library written at molecular resolution would therefore give `GM3 d18:1_24:0` and
`GM3 d20:1_22:0` *identical* spectra, and the dot product could not separate them — the same
single-fragment degeneracy documented for the 184.07 phosphocholine entries in VALIDATION.md,
in a more extreme form. Writing it at sum composition instead makes the library say exactly what
the evidence says: the precursor mass gives the total carbons and double bonds, the fragments
give the head group, and the chain split is not determined.

The entries carry no moiety-fragment annotations, so the purity calculation returns 0 and the
peak finder reports them at sum composition of its own accord. The library and the reporting
agree because both reflect the same absence of evidence.

    python scripts/build_ganglioside_negative_library.py \
        --positive data/libraries/LipiDex2_Ganglioside.msp \
        --out data/libraries/Ganglioside_Negative_SumComposition.msp
"""
from __future__ import annotations

import argparse
import re
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.msp import parse_msp          # noqa: E402
from lipidloop.peaks import sum_composition  # noqa: E402

PROTON = 1.00727646

# Sialic acid fragments, from the observed spectra. Intensities are the mean relative to the
# base peak across the five species independently confirmed in positive mode — a template
# calibrated on confirmed examples, not fitted to the species being searched for.
NEGATIVE_FRAGMENTS = [
    (290.0881, 999.0, "C11H16NO8_Fragment_[]"),     # [NeuAc - H2O - H]-
    (87.0088, 130.0, "C3H3O3_Fragment_[]"),         # sialic acid ring fragment
    (272.0776, 23.0, "C11H14NO7_Fragment_[]"),      # [NeuAc - 2H2O - H]-
]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--positive", required=True, help="LipiDex 2 positive ganglioside .msp")
    ap.add_argument("--out", required=True)
    ap.add_argument("--classes", nargs="*", default=None,
                    help="restrict to these classes, e.g. GM3-NANA (default: all found)")
    args = ap.parse_args()

    # One entry per distinct sum composition. The species list comes from LipiDex 2's own
    # enumeration rather than being invented here; only the resolution is reduced.
    by_sum: dict[tuple[str, float], list[str]] = defaultdict(list)
    for spectrum in parse_msp(args.positive):
        if args.classes and spectrum.lipid_class not in args.classes:
            continue
        molecular = spectrum.lipid
        summed = sum_composition(molecular)
        if summed == molecular and "_" in molecular:
            continue                                  # could not be collapsed; skip
        neutral = spectrum.precursor_mz - PROTON      # entries are all [M+H]+
        by_sum[(summed, round(neutral, 4))].append(molecular)

    lines = []
    for (summed, neutral), members in sorted(by_sum.items()):
        precursor = neutral - PROTON                  # [M-H]-
        peaks = [(mz, inten, ann) for mz, inten, ann in NEGATIVE_FRAGMENTS
                 if precursor - mz > 2.0]
        if not peaks:
            continue
        lines.append(f"Name: {summed} [M-H]-;")
        lines.append(f"MW: {neutral:.4f}")
        lines.append(f"PRECURSORMZ: {precursor:.4f}")
        lines.append(f"Comment: Name={summed} [M-H]- Mass={precursor:.4f} "
                     f"OptimalPolarity=true Type=LipiDex "
                     f"CollapsedFrom={len(members)} molecular species")
        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, inten, annotation in peaks:
            lines.append(f'{mz:.4f} {inten:.0f} "{annotation}"')
        lines.append("")

    Path(args.out).write_text("\n".join(lines))
    collapsed = sum(len(m) for m in by_sum.values())
    print(f"{len(by_sum)} sum compositions written to {args.out}")
    print(f"  collapsed from {collapsed} molecular species in the positive library")
    print(f"  {len(NEGATIVE_FRAGMENTS)} fragments each, all sialic acid — the spectra carry no "
          f"chain information, so neither does the library")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
