"""Build an in-silico library for free sterols, which no LipiDex library covers.

**NOT IN THE DEFAULT LIBRARY SET, AND READ docs/STEROLS.md BEFORE USING IT.** The spectrum below
is a correct cholesterol spectrum. The problem is what else produces it: cholesteryl esters
fragment in-source to the same m/z 369.3516 and the same backbone ions, so on a run with
abundant CE this library matches the CE in-source fragment and reports it as free cholesterol.
That is exactly what happened when it was first tried here — every hit turned out to be CE.

It is usable with an authentic cholesterol standard to fix the retention window, and not
otherwise.

Neither LipiDex version carries free cholesterol. `CE` is there — 40 cholesteryl esters, all
`[M+NH4]+` — but the free sterol is absent from every library in the default set. That gap is
real. Filling it safely is the hard part.

The entry below is built from real spectra: filtering positive-mode data to precursor 369.3516
and requiring the sterol backbone ions gives eleven spectra carrying 54–81 peaks each, and the
fragment masses agree with the backbone series LipiDex already carries inside its CE entries.
The spectrum is right.

**What it identifies may not be.** On the validation set every one of those eleven spectra turned
out to be a cholesteryl ester fragmenting in the source. The 369.3516 chromatogram has exactly
two peaks, at 16.09 and 17.09 min, and they sit on CE 20:4 (16.09) and CE 18:1 (17.08); a third
"cholesterol" row at 16.49 sits on CE 18:2 (16.48). Between 2 and 14 min, where free cholesterol
would elute, the trace is flat background. There was no free cholesterol to find.

So use this only with an authentic cholesterol standard run on the same method, to fix the
retention window and reject anything eluting in the CE region.

    python scripts/build_sterol_library.py --out data/libraries/Sterols_InSilico.msp
"""
from __future__ import annotations

import argparse
from pathlib import Path

PROTON = 1.00727646
WATER = 18.0105646

# Free sterols, as the neutral molecule. Add to this list to extend it; each is emitted as
# [M+H-H2O]+, which is how sterols are normally observed in positive mode.
STEROLS = [
    # name,        neutral monoisotopic mass, formula
    ("Cholesterol", 386.354866, "C27H46O"),
]

# The cholesterol fragment pattern, measured from the confirmed spectra in this dataset and
# normalised to a base peak of 999. Masses agree with the sterol backbone series LipiDex carries
# inside its CE entries.
CHOLESTEROL_FRAGMENTS = [
    (147.1168, 999, "C11H15"), (161.1325, 847, "C12H17"), (135.1168, 633, "C10H15"),
    (95.0855, 588, "C7H11"), (109.1012, 572, "C8H13"), (81.0699, 468, "C6H9"),
    (175.1481, 342, "C13H19"), (149.1325, 333, "C11H17"), (121.1012, 260, "C9H13"),
    (243.2107, 259, "C18H27"), (189.1638, 249, "C14H21"), (215.1794, 237, "C16H23"),
    (203.1794, 190, "C15H23"), (257.2264, 150, "C19H29"), (229.1951, 140, "C17H25"),
]

PRECURSOR_CUTOFF = 2.0    # the library reader drops anything closer than this to the precursor


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lines = []
    for name, neutral, formula in STEROLS:
        precursor = neutral + PROTON - WATER
        peaks = [(mz, i, f) for mz, i, f in CHOLESTEROL_FRAGMENTS
                 if precursor - mz > PRECURSOR_CUTOFF]
        lines.append(f"Name: {name} [M+H-H2O]+;")
        lines.append(f"MW: {neutral:.4f}")
        lines.append(f"PRECURSORMZ: {precursor:.4f}")
        # No moiety annotations: the spectrum carries no chain information, so the purity
        # calculation returns 0 and the peak finder reports the name as it stands.
        lines.append(f"Comment: Name={name} [M+H-H2O]+ Mass={precursor:.4f} "
                     f"Formula={formula} OptimalPolarity=true Type=LipiDex")
        lines.append(f"Num Peaks: {len(peaks)}")
        for mz, intensity, frag_formula in peaks:
            lines.append(f'{mz:.4f} {intensity} "{frag_formula}_Fragment_[]"')
        lines.append("")

    Path(args.out).write_text("\n".join(lines))
    print(f"{len(STEROLS)} sterol(s) written to {args.out}")
    for name, neutral, _ in STEROLS:
        print(f"  {name:14} [M+H-H2O]+ {neutral + PROTON - WATER:9.4f}  "
              f"{len(CHOLESTEROL_FRAGMENTS)} fragments")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
