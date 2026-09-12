"""Find m/z that appear in nearly every MS2 spectrum — instrument artefacts, not chemistry.

A fragment that appears whatever the precursor was is not a fragment. It matters more than it
sounds: every spectrum is normalised to its base peak before scoring, so an artefact that *is*
the base peak scales every genuine fragment down against it, and its unmatched intensity counts
against the forward dot product as well. On this instrument one such artefact costs **36% of
identifications in positive mode and 50% in negative**.

Run this on every new batch rather than reusing a number measured elsewhere. An artefact belongs
to one instrument in one period, not to lipidomics.

**Bin width is wide on purpose.** An artefact's m/z wanders — the one here spreads over about
60 mDa — and a bin narrower than that spread splits it across bins and hides it. With 0.01 bins
the artefact reads as 50% of spectra; with 0.10 bins it reads as 100%, which is the truth.

    python scripts/find_spectral_artefacts.py --mzml data/mzml/val_Pos/QC_01.mzML
"""
from __future__ import annotations

import argparse
import collections
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.mgf import read_mzml_ms2   # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mzml", required=True)
    ap.add_argument("--bin", type=float, default=0.10,
                    help="m/z bin width; wide on purpose, see the module docstring")
    ap.add_argument("--min-fraction", type=float, default=0.3,
                    help="report m/z present in at least this fraction of spectra")
    args = ap.parse_args()

    present = collections.Counter()
    is_base = collections.Counter()
    observed = collections.defaultdict(list)
    total = 0

    for ms2 in read_mzml_ms2(args.mzml):
        if not ms2.mz:
            continue
        total += 1
        top = max(zip(ms2.intensity, ms2.mz))
        is_base[round(top[1] / args.bin)] += 1
        for key in {round(mz / args.bin) for mz in ms2.mz}:
            present[key] += 1
        for mz in ms2.mz:
            observed[round(mz / args.bin)].append(mz)

    print(f"{total} MS2 spectra in {Path(args.mzml).name}\n")
    print(f"{'mean m/z':>10} {'in %':>7} {'base peak %':>12} {'m/z SD (mDa)':>13}   verdict")

    recommended = []
    for key, count in present.most_common(30):
        fraction = count / max(total, 1)
        if fraction < args.min_fraction:
            break
        base_fraction = is_base[key] / max(total, 1)
        values = observed[key]
        spread = 1000 * statistics.stdev(values) if len(values) > 2 else 0.0
        mean_mz = statistics.mean(values)

        if fraction > 0.95 and base_fraction > 0.2:
            verdict = "ARTEFACT — everywhere, usually the base peak"
            recommended.append(round(mean_mz, 2))
        elif fraction > 0.95:
            verdict = "ubiquitous — check it has a plausible formula"
        else:
            verdict = "common"
        print(f"{mean_mz:10.4f} {100*fraction:6.1f}% {100*base_fraction:11.1f}% "
              f"{spread:13.1f}   {verdict}")

    print("\nTwo things mark an artefact rather than a real fragment: it appears whatever the")
    print("precursor was, and its m/z wanders. A real ion holds a stable m/z, and a singly-")
    print("charged organic ion cannot have an arbitrary mass defect — check the formula before")
    print("excluding anything.")
    if recommended:
        print(f'\nadd to the run configuration:    "exclude_mz": {recommended}')
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
