"""Why is a reference molecule missing from our results? Trace it to the stage that lost it.

"We recover 96% of the reference" is only useful if the other 4% can be accounted for. Each
missing molecule fell out at exactly one place, and the four are different problems with different
fixes:

  **never fragmented** — no MS2 anywhere in the batch was identified as it. Nothing downstream
    could have recovered it. If there is an MS1 feature at its mass and retention time, this is a
    data-dependent acquisition problem, not a processing one, and it is what an inclusion list or
    a freed duty cycle would fix.
  **identified, no feature** — spectra were identified but no MS1 feature carried them, so there
    was nothing to quantify. A feature-detection sensitivity problem.
  **filtered** — it reached `Unfiltered_Results.csv` and a named filter dropped it. Recoverable by
    changing a threshold, and `Filter Status` says which one.
  **named differently** — present at the same mass and retention time under another name. A
    reporting-resolution or tie difference, not a loss.

    python scripts/trace_missing.py \\
        --ours OmicsPilot_Lipidomics/Pos/Final_Results.csv \\
        --unfiltered OmicsPilot_Lipidomics/Pos/Unfiltered_Results.csv \\
        --per-file OmicsPilot_Lipidomics/Pos \\
        --reference '.../Pos/CD/Final_ResultsPos.csv' \\
        --report results/missing_trace_Pos.md
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.peaks import sum_composition   # noqa: E402

MZ_TOL = 0.02
RT_TOL = 0.5


def read(path):
    with Path(path).open(newline="", errors="replace") as fh:
        return list(csv.DictReader(fh))


def carbons(name: str) -> "tuple[int, int] | None":
    """Total carbons and double bonds from a sum-composition name, if it has them."""
    match = re.search(r"[dtoOP]?-?(\d+):(\d+)", name)
    return (int(match.group(1)), int(match.group(2))) if match else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--unfiltered")
    ap.add_argument("--per-file", help="folder of <stem>_Results.csv")
    ap.add_argument("--reference", required=True)
    ap.add_argument("--report")
    args = ap.parse_args()

    our_rows = read(args.ours)
    ref_rows = read(args.reference)
    our_sums = {sum_composition(r["Identification"].strip())
                for r in our_rows if r.get("Identification", "").strip()}
    missing = {}
    for row in ref_rows:
        name = row.get("Identification", "").strip()
        if not name:
            continue
        summed = sum_composition(name)
        if summed not in our_sums:
            missing.setdefault(summed, row)
    print(f"{len(missing)} reference molecules missing from {args.ours}")

    # 1. was any spectrum in the batch identified as it?
    identified_anywhere = defaultdict(list)
    if args.per_file:
        for path in glob.glob(os.path.join(args.per_file, "*_Results.csv")):
            if os.path.basename(path) in ("Final_Results.csv", "Unfiltered_Results.csv"):
                continue
            for row in read(path):
                summed = sum_composition(row.get("Identification", "").strip())
                if summed in missing:
                    identified_anywhere[summed].append((os.path.basename(path), row))

    # 2. did it reach the unfiltered table, and under what filter status?
    unfiltered_hit, renamed = {}, {}
    if args.unfiltered:
        rows = read(args.unfiltered)
        for summed, ref in missing.items():
            mz, rt = float(ref["Quant Ion"]), float(ref["Retention Time (min)"])
            for row in rows:
                if abs(float(row["Quant Ion"]) - mz) > MZ_TOL:
                    continue
                if abs(float(row["Retention Time (min)"]) - rt) > RT_TOL:
                    continue
                name = row.get("Identification", "").strip()
                if sum_composition(name) == summed:
                    unfiltered_hit[summed] = row.get("Filter Status", "").strip() or "kept"
                elif name and summed not in renamed:
                    renamed[summed] = (name, row.get("Filter Status", "").strip())

    verdicts = {}
    for summed in missing:
        if summed in unfiltered_hit:
            verdicts[summed] = f"filtered: {unfiltered_hit[summed]}"
        elif identified_anywhere.get(summed):
            verdicts[summed] = "identified, no feature"
        elif summed in renamed:
            verdicts[summed] = f"named differently: {renamed[summed][0]}"
        else:
            verdicts[summed] = "never fragmented"

    tally = Counter(v.split(":")[0] for v in verdicts.values())
    print()
    for verdict, count in tally.most_common():
        print(f"  {verdict:26} {count:4d}")

    # does chain length separate the verdicts?
    print(f"\n{'verdict':26} {'n':>4} {'median C':>9} {'median DB':>10} "
          f"{'median ref area':>16}")
    by_verdict = defaultdict(list)
    for summed, verdict in verdicts.items():
        by_verdict[verdict.split(":")[0]].append(summed)
    for verdict, names in sorted(by_verdict.items(), key=lambda kv: -len(kv[1])):
        chains = [c for n in names if (c := carbons(n))]
        areas = sorted(float(missing[n]["Area (max)"]) for n in names)
        cs = sorted(c for c, _ in chains)
        dbs = sorted(d for _, d in chains)
        print(f"  {verdict:24} {len(names):4d} {cs[len(cs)//2] if cs else 0:9d} "
              f"{dbs[len(dbs)//2] if dbs else 0:10d} {areas[len(areas)//2]:16.1e}")

    # and how do the ones we DO get compare?
    kept = [sum_composition(r["Identification"].strip()) for r in ref_rows
            if r.get("Identification", "").strip()
            and sum_composition(r["Identification"].strip()) in our_sums]
    kept_areas = sorted(float(r["Area (max)"]) for r in ref_rows
                        if sum_composition(r.get("Identification", "").strip()) in our_sums
                        and r.get("Identification", "").strip())
    kept_c = sorted(c for n in kept if (p := carbons(n)) for c in [p[0]])
    if kept_areas:
        print(f"  {'RECOVERED (for comparison)':24} {len(kept):4d} "
              f"{kept_c[len(kept_c)//2] if kept_c else 0:9d} "
              f"{'':10} {kept_areas[len(kept_areas)//2]:16.1e}")

    if args.report:
        write_report(args.report, args=args, missing=missing, verdicts=verdicts,
                     identified_anywhere=identified_anywhere, tally=tally,
                     kept_areas=kept_areas, kept_c=kept_c)
        print(f"\nwritten to {args.report}")
    return 0


def write_report(path, *, args, missing, verdicts, identified_anywhere, tally,
                 kept_areas, kept_c) -> None:
    from datetime import date

    lines = [f"# Where the missing reference molecules were lost", "",
             f"Generated {date.today().isoformat()} by `scripts/trace_missing.py`.", "",
             f"- ours: `{args.ours}`", f"- reference: `{args.reference}`", "",
             f"**{len(missing)} reference molecules are missing.** Each fell out at one stage:",
             "", "| verdict | count | what it means |", "|---|---|---|"]
    meaning = {
        "never fragmented": "no MS2 anywhere in the batch was identified as it — nothing "
                            "downstream could recover it",
        "identified, no feature": "spectra were identified but no MS1 feature carried them",
        "filtered": "reached the unfiltered table and a named filter dropped it",
        "named differently": "present at the same mass and retention time under another name",
    }
    for verdict, count in tally.most_common():
        lines.append(f"| {verdict} | {count} | {meaning.get(verdict, '')} |")

    lines += ["", "## Every one", "",
              "| molecule | reference RT | reference m/z | reference area | verdict |",
              "|---|---|---|---|---|"]
    for summed in sorted(missing, key=lambda s: -float(missing[s]["Area (max)"])):
        ref = missing[summed]
        lines.append(f"| `{summed}` | {float(ref['Retention Time (min)']):.2f} | "
                     f"{float(ref['Quant Ion']):.4f} | {float(ref['Area (max)']):.1e} | "
                     f"{verdicts[summed]} |")

    if kept_areas:
        median_missing = sorted(float(r["Area (max)"]) for r in missing.values())
        lines += ["", "## Abundance", "",
                  "| | n | median reference area |", "|---|---|---|",
                  f"| missing | {len(missing)} | {median_missing[len(median_missing)//2]:.1e} |",
                  f"| recovered | {len(kept_areas)} | {kept_areas[len(kept_areas)//2]:.1e} |", ""]
    Path(path).write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
