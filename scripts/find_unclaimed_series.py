"""Find homologous series among the features nothing identified — the shape of a missing library.

The free-fatty-acid work started from an observation that generalises: a lipid class the libraries
do not cover does not show up as one anonymous peak, it shows up as a **ladder**. Members differ by
CH2 (14.0157) at a steady retention step, and by a double bond (-2.0157) at a smaller negative one.
Nothing else in a chromatogram does that. So the ladders among unidentified features are a direct
readout of which library to build next, ranked by how much signal each would claim.

    python scripts/find_unclaimed_series.py \\
        --results OmicsPilot_Lipidomics/Neg/Final_Results.csv --polarity - --top 15

Reports each series with its members, the retention step per carbon, and the total area it holds.
A series is not an identification — it says a coherent, unclaimed family of compounds is present at
these masses, which is the question worth taking to a library.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

CH2 = 14.01565
TOLERANCE = 0.01
MIN_MEMBERS = 4
MIN_RT_STEP = 0.05      # minutes: members must actually separate, or they are not homologues
MAX_RT_STEP = 3.0


def read_unidentified(path, polarity, min_files):
    with Path(path).open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    out = []
    for row in rows:
        if row.get("Identification", "").strip():
            continue
        if row.get("Filter Status", "").strip():          # already explained as adduct/dimer/...
            continue
        if row.get("Polarity", "").strip() != polarity:
            continue
        if int(row.get("Features Found", 0)) < min_files:
            continue
        out.append((float(row["Quant Ion"]), float(row["Retention Time (min)"]),
                    float(row["Area (max)"]), int(row["Features Found"])))
    out.sort()
    return out


def chain(features, start, tolerance):
    """Walk up from `start` in CH2 steps, taking the most abundant match at each rung."""
    members, probe = [start], start
    while True:
        target = probe[0] + CH2
        candidates = [f for f in features if abs(f[0] - target) <= tolerance]
        if not candidates:
            break
        # Retention must rise: a homologue with one more carbon elutes later, always.
        candidates = [f for f in candidates if f[1] > probe[1]]
        if not candidates:
            break
        probe = max(candidates, key=lambda f: f[2])
        members.append(probe)
    return members


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", required=True)
    ap.add_argument("--polarity", default="-")
    ap.add_argument("--min-files", type=int, default=10)
    ap.add_argument("--min-members", type=int, default=MIN_MEMBERS)
    ap.add_argument("--top", type=int, default=15)
    ap.add_argument("--tolerance", type=float, default=TOLERANCE)
    args = ap.parse_args()

    features = read_unidentified(args.results, args.polarity, args.min_files)
    print(f"{len(features)} unidentified, unfiltered features in {args.polarity} mode "
          f"seen in >= {args.min_files} files\n")

    claimed, series = set(), []
    for start in sorted(features, key=lambda f: -f[2]):
        if start[0] in claimed:
            continue
        members = chain(features, start, args.tolerance)
        if len(members) < args.min_members:
            continue
        steps = [b[1] - a[1] for a, b in zip(members, members[1:])]
        step = sum(steps) / len(steps)
        if not (MIN_RT_STEP <= step <= MAX_RT_STEP):
            continue
        for m in members:
            claimed.add(m[0])
        series.append((members, step))

    series.sort(key=lambda s: -sum(m[2] for m in s[0]))
    print(f"{len(series)} unclaimed homologous series\n")
    for members, step in series[:args.top]:
        total = sum(m[2] for m in members)
        print(f"  {len(members)} members, {members[0][0]:.4f} -> {members[-1][0]:.4f}, "
              f"RT {members[0][1]:.2f} -> {members[-1][1]:.2f} "
              f"({step:+.2f} min per CH2), total area {total:.1e}")
        for mz, rt, area, files in members:
            print(f"      {mz:9.4f}  RT {rt:5.2f}  area {area:8.1e}  in {files} files")
        print()
    if series:
        share = sum(sum(m[2] for m in s[0]) for s in series) / max(
            sum(f[2] for f in features), 1)
        print(f"these series hold {100 * share:.0f}% of the unidentified signal in this polarity")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
