"""Compare a finished `Final_Results.csv` against a Compound Discoverer + LipiDex reference.

`validate_experiment.py` runs the pipeline itself and is wired to one dataset. This takes two
result files that already exist, so any run can be scored against any reference without
recomputing it.

    python scripts/compare_to_reference.py \\
        --ours OmicsPilot_Lipidomics/Pos/Final_Results.csv \\
        --reference '.../Validation_2_Skin/Pos/CD/Final_ResultsPos.csv' \\
        --label 'skin organoid, positive'

**Recovery is reported two ways and the difference is not cosmetic.** Matching identification
strings exactly counts `TG 12:0_16:0_20:4` against `TG 48:4` as both a miss and an extra, when it
is one molecule named at two resolutions — the resolution each pipeline reports depends on whether
its fatty-acid purity cleared 75%. Collapsing both sides to sum composition is the honest number.
"""
from __future__ import annotations

import argparse
import csv
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.peaks import sum_composition   # noqa: E402


def identifications(path):
    with Path(path).open(newline="", errors="replace") as fh:
        rows = list(csv.DictReader(fh))
    ids = {r["Identification"].strip() for r in rows if r.get("Identification", "").strip()}
    return rows, ids


def lipid_class(name: str) -> str:
    return name.split()[0] if name.split() else ""


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--ours", required=True)
    ap.add_argument("--reference", required=True)
    ap.add_argument("--label", default="")
    ap.add_argument("--report")
    ap.add_argument("--top-classes", type=int, default=15)
    args = ap.parse_args()

    our_rows, our_ids = identifications(args.ours)
    ref_rows, ref_ids = identifications(args.reference)
    our_sums = {sum_composition(n) for n in our_ids}
    ref_sums = {sum_composition(n) for n in ref_ids}
    shared = our_sums & ref_sums

    lines = [f"# Against the Compound Discoverer + LipiDex reference"
             + (f" — {args.label}" if args.label else ""), "",
             f"- ours: `{args.ours}`", f"- reference: `{args.reference}`", "",
             "| | ours | reference |", "|---|---|---|",
             f"| rows | {len(our_rows):,} | {len(ref_rows):,} |",
             f"| identified rows | {sum(1 for r in our_rows if r.get('Identification','').strip()):,} "
             f"| {sum(1 for r in ref_rows if r.get('Identification','').strip()):,} |",
             f"| distinct lipids (exact name) | {len(our_ids):,} | {len(ref_ids):,} |",
             f"| distinct molecules (sum composition) | {len(our_sums):,} | {len(ref_sums):,} |", "",
             "| | count | of the reference |", "|---|---|---|",
             f"| **shared, same molecule** | **{len(shared):,}** | "
             f"**{100 * len(shared) / max(len(ref_sums), 1):.1f}%** |",
             f"| shared, exact name | {len(our_ids & ref_ids):,} | "
             f"{100 * len(our_ids & ref_ids) / max(len(ref_ids), 1):.1f}% |",
             f"| only ours | {len(our_sums - ref_sums):,} | — |",
             f"| only the reference | {len(ref_sums - our_sums):,} | — |", ""]

    ours_by_class = Counter(lipid_class(n) for n in our_sums)
    ref_by_class = Counter(lipid_class(n) for n in ref_sums)
    shared_by_class = Counter(lipid_class(n) for n in shared)
    lines += ["## By lipid class", "",
              "| class | ours | reference | shared |", "|---|---|---|---|"]
    for name, _ in ref_by_class.most_common(args.top_classes):
        lines.append(f"| {name} | {ours_by_class[name]} | {ref_by_class[name]} | "
                     f"{shared_by_class[name]} |")
    extra = [c for c, _ in ours_by_class.most_common() if c not in ref_by_class][:6]
    if extra:
        lines += ["", "Classes only we report: "
                  + ", ".join(f"{c} ({ours_by_class[c]})" for c in extra), ""]

    missed = sorted(ref_sums - our_sums)
    lines += ["", f"## The {len(missed)} reference molecules we do not report", "",
              ", ".join(f"`{m}`" for m in missed) if missed else "None.", "",
              f"## The {len(our_sums - ref_sums)} molecules only we report", "",
              "**Unverified.** Some will be real lipids the reference missed and some will not; "
              "the reference was produced from the same spectra and shares this pipeline's blind "
              "spots as well as its own.", "",
              ", ".join(f"`{m}`" for m in sorted(our_sums - ref_sums)[:80])
              + (" …" if len(our_sums - ref_sums) > 80 else ""), ""]

    text = "\n".join(lines)
    print(text if not args.report else "\n".join(lines[:28]))
    if args.report:
        Path(args.report).parent.mkdir(parents=True, exist_ok=True)
        Path(args.report).write_text(text)
        print(f"\nwritten to {args.report}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
