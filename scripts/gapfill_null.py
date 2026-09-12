"""What S/N does an EMPTY position score? The threshold cannot be chosen without knowing.

Compound Discoverer's Fill Gaps ships S/N 1.5, and a threshold is only meaningful against the
noise it is applied to. This probes decoy coordinates — real retention times at m/z offsets where
no lipid can be — and reports the distribution of scores there. Anything at or below that level is
indistinguishable from nothing.

    python scripts/gapfill_null.py --mzml <file.mzML> --results <Final_Results.csv>
"""
import argparse, csv, sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.gapfill import _ms1, _peak, _xic   # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mzml", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--rt-window", type=float, default=0.2)
    ap.add_argument("--ppm", type=float, default=10.0)
    ap.add_argument("--n", type=int, default=400)
    args = ap.parse_args()

    rows = [r for r in csv.DictReader(Path(args.results).open(newline=""))]
    rows = [r for r in rows if r.get("Quant Ion")]
    scans = _ms1(args.mzml)
    times = np.array([s[0] for s in scans])
    rng = np.random.default_rng(0)

    real, decoy, rt_decoy = [], [], []
    picks = rng.choice(len(rows), size=min(args.n, len(rows)), replace=False)
    for i in picks:
        mz = float(rows[i]["Quant Ion"]); rt = float(rows[i]["Retention Time (min)"])
        # Two decoys, because they test different things. The mass decoy (+0.3713 Da, a defect
        # no CHNOPS composition reaches) asks whether the window is specific in m/z. The RT decoy
        # — the SAME mass, 3 min away — asks the question that actually matters for gap filling:
        # does the window pick up unrelated signal at the right mass but the wrong time? An
        # isobaric neighbour or an in-source fragment lives exactly there.
        for target, offset, out in ((mz, 0.0, real), (mz + 0.3713, 0.0, decoy),
                                    (mz, 3.0, rt_decoy)):
            centre = rt + offset
            a, s = _xic(scans, times, target, args.ppm,
                        centre - args.rt_window * 3, centre + args.rt_window * 3)
            found = _peak(a, s, 0.0, centre, args.rt_window)
            out.append(found[1] if found else 0.0)

    real, decoy, rt_decoy = np.array(real), np.array(decoy), np.array(rt_decoy)
    print(f"{len(real)} positions probed in {Path(args.mzml).name}\n")
    print(f"{'percentile':>12} {'at the peak':>14} {'wrong mass':>14} {'wrong time':>14}")
    for q in (50, 75, 90, 95, 99):
        print(f"{q:>11}% {np.percentile(real,q):14.1f} {np.percentile(decoy,q):14.1f} "
              f"{np.percentile(rt_decoy,q):14.1f}")
    print(f"\n  scoring a peak at all:  wrong mass {100*(decoy>0).mean():.1f}%   "
          f"wrong time {100*(rt_decoy>0).mean():.1f}%")
    print(f"\n  {'threshold':>10} {'keeps real':>12} {'wrong mass':>12} {'wrong time':>12}")
    for threshold in (1.5, 3, 5, 10, 20, 50):
        print(f"  {threshold:>9} {100*(real>=threshold).mean():11.1f}% "
              f"{100*(decoy>=threshold).mean():11.1f}% {100*(rt_decoy>=threshold).mean():11.1f}%")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
