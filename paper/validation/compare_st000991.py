"""Old run against new, for ST000991 specifically.

Kept out of the shell script because `score_agreement.py` resolves a dataset to ONE analysis
directory, preferring Analysis_overnight — calling it here would re-score the old run and report
no change whatever the new one did.
"""
import csv
import glob
import sys

sys.path.insert(0, "/home/jair/validation_bench/agreement")
from score_agreement import canon, published  # noqa: E402

pub = {s for _, s in map(canon, published("ST000991")) if s}
for label, d in (("old (all 159 searched)", "Analysis_overnight"),
                 ("new (per-file filter)", "Analysis_perfile")):
    names = []
    for f in glob.glob(f"/home/jair/validation_staged/ST000991/{d}/Results/*/Unfiltered_Results.csv"):
        with open(f, errors="replace") as fh:
            for r in csv.DictReader(fh):
                i = (r.get("Identification") or "").strip()
                if i and not i.upper().startswith("DECOY"):
                    names.append((r.get("Shorthand (LSI)") or i).strip())
    ours = {s for _, s in map(canon, names) if s}
    if not ours:
        print(f"  {label:<26} no results")
        continue
    shared = len(pub & ours)
    print(f"  {label:<26} {len(ours):>4} species  {shared:>4} shared  "
          f"{100 * shared / len(pub):>5.1f}% of theirs  {100 * shared / len(ours):>5.1f}% of ours")
