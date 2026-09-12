"""Write an Xcalibur-style sequence CSV per MTBKS222 arm, from the deposit's own SDRF.

⚠ **The deposit declares its blanks and we were ignoring them.** MTBKS222's SDRF carries a
`Source Name` of `blank sample` and an `Assay Name` naming the file — `Waters_190918_090pp_Pos`.
Without that, `run_study` infers roles from filenames, and `190918_090pp` looks like nothing in
particular, so **every blank was treated as a sample**.

The cost is not cosmetic. On the Waters arm 10 of 18 runs are blanks, and with all of them in one
group the presence filter requires a feature to recur across injections that include solvent. A
real plasma lipid is absent from a blank by definition, so it is filtered out; only 83 of 2,168
features survived, and they survived *because* they also appear in the blanks. That arm's 86.2%
agreement was measured on background.

Roles are read, never guessed: `blank sample` in the SDRF becomes `Blank`, which
`lipidloop.blanks.infer_role` recognises. Everything else stays `Unknown`, which reads as a
sample — the deposit does not distinguish pools, so neither do we.
"""
from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

DEPOSIT = Path("/mnt/datastore/Jair/claudecode/Lipidomics_automation/Validation_Datasets/MTBKS222")
STAGED = Path("/home/jair/validation_staged")


def sdrf_roles() -> dict[str, str]:
    """{file stem: Sample Type} for every run the SDRF describes."""
    out: dict[str, str] = {}
    with (DEPOSIT / "MTBKS222.sdrf.txt").open(errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            assay = (row.get("Assay Name") or "").strip()
            source = (row.get("Source Name") or "").strip().lower()
            if not assay:
                continue
            # The assay name embeds the raw file stem with a polarity suffix appended.
            stem = re.sub(r"_(Pos|Neg)$", "", assay.split("_", 1)[-1] if assay.startswith(
                ("Waters_", "Thermo_", "Agilent_")) else assay)
            for candidate in {assay, stem}:
                out[candidate] = "Blank" if "blank" in source else "Unknown"
    return out


def main() -> int:
    roles = sdrf_roles()
    total = Counter()
    for arm in sorted(STAGED.glob("MTBKS222_*")):
        # ⚠ ONE sequence per arm covering BOTH polarities. `run_study.find_sequence` returns a
        # single file for the whole study, so two per-polarity sequences would leave whichever
        # one it did not pick entirely unmapped — and the blanks in that polarity would go on
        # being treated as samples, which is the bug this exists to fix.
        rows, seen = [], Counter()
        for polarity in ("Pos", "Neg"):
            for f in sorted((arm / polarity).glob("*.mzML")):
                stem = f.stem
                role = roles.get(stem) or roles.get(f"{stem}_{polarity}")
                if role is None:
                    match = [v for k, v in roles.items() if stem in k]
                    role = match[0] if match else "Unknown"
                rows.append({"File Name": stem, "Sample Type": role})
                seen[role] += 1
        if not rows:
            continue
        out = arm / "sequence.csv"
        with out.open("w", newline="") as fh:
            fh.write("Bracket Type=4\n")
            w = csv.DictWriter(fh, fieldnames=["File Name", "Sample Type"])
            w.writeheader()
            w.writerows(rows)
        for stale in arm.glob("sequence_*.csv"):
            stale.unlink()
        print(f"  {arm.name:<26} {len(rows):>3} runs  {dict(seen)}")
        total.update(seen)
    print(f"\n  total across arms: {dict(total)}")
    if not total.get("Blank"):
        print("  ⚠ NO BLANKS FOUND — the SDRF mapping failed; do not re-run on these")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
