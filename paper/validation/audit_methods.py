"""What each deposit ACTUALLY acquired, read from the files rather than from its own description.

Two datasets were queued, run, and found unusable one at a time — ST000987 fragments without
isolating (all-ion, in-source CID) and MTBLS2016 contains no MS2 at all despite being described as
data-dependent. Each cost a run and a diagnosis. The remaining deposits deserve to be checked
before they cost the same, and the check is seconds per dataset against hours per run.

Read from the spectra because deposits are unreliable about this: MTBLS2016's own record says DDA.

Verdicts:
  DDA        MS2 with narrow isolation windows — searchable, this pipeline's premise holds
  DIA        MS2 with wide, repeatedly reused windows — co-fragmented, premise does not hold
  all-ion    MS2 with no isolation at all — nothing was ever selected
  MS1 only   no MS2 — nothing to search
"""
from __future__ import annotations

import json
import re
from pathlib import Path

STAGED = Path("/home/jair/validation_staged")
META = Path("/tmp/claude-1000/-home-jair/1e5bf05d-c2a1-4067-b572-5355f626c8dd/scratchpad/meta")
BUDGET = 120 << 20        # bytes per file; enough for thousands of spectra, bounded on 900 MB files
SAMPLE = 6                # files per polarity

# ⚠ One file per dataset is not enough, and reading only one produced three wrong verdicts.
# ST003514 was called "MS1 only" on a study that identified 591 lipids, because these deposits mix
# acquisitions WITHIN one submission — its own record says "iterative-MS/MS runs were performed",
# meaning survey injections with no MS2 sit beside the DDA ones. A dataset's method is a
# DISTRIBUTION over its files, not a property of the middle one.

_LEVEL = re.compile(r'name="ms level" value="(\d)"')
_TARGET = re.compile(r'name="isolation window target m/z" value="([\d.]+)"')
_LOWER = re.compile(r'name="isolation window lower offset" value="([\d.]+)"')
_SELECTED = re.compile(r"selectedIon")
_ALLION = re.compile(r"MS:1001880|In-source collision-induced dissociation")


def scan(path: Path) -> dict:
    ms1 = ms2 = selected = allion = 0
    targets: list[float] = []
    widths: list[float] = []
    read = 0
    with path.open(errors="replace") as fh:
        while read < BUDGET:
            chunk = fh.read(1 << 20)
            if not chunk:
                break
            read += len(chunk)
            for v in _LEVEL.findall(chunk):
                if v == "1":
                    ms1 += 1
                else:
                    ms2 += 1
            selected += len(_SELECTED.findall(chunk))
            allion += len(_ALLION.findall(chunk))
            targets.extend(float(x) for x in _TARGET.findall(chunk))
            widths.extend(float(x) for x in _LOWER.findall(chunk))
    return {"ms1": ms1, "ms2": ms2, "selected": selected, "allion": allion,
            "targets": targets, "widths": widths}


def verdict(s: dict) -> tuple[str, str]:
    if not s["ms2"]:
        return "MS1 only", "no MS2 in the sampled portion"
    if not s["selected"] and not s["targets"]:
        return "all-ion", ("in-source CID" if s["allion"] else "no isolation stated")
    if s["widths"]:
        width = max(s["widths"]) * 2
        uniq = len(set(round(t, 1) for t in s["targets"])) or 1
        reuse = len(s["targets"]) / uniq
        if width >= 5.0 and reuse >= 3.0:
            return "DIA", f"{width:.0f} Da windows, {uniq} of them, reused {reuse:.0f}x"
        return "DDA", f"{width:.1f} Da isolation"
    return "DDA", "precursors selected"


def claimed(dataset: str) -> str:
    p = META / f"{dataset}.mwtab.txt"
    if p.exists():
        for line in p.open(errors="replace"):
            if line.startswith(("MS:MS_COMMENTS", "AN:ACQUISITION", "MS:INSTRUMENT_NAME")):
                v = line.partition("\t")[2].strip()
                if re.search(r"DDA|data.dependent|DIA|SWATH|MSe|all.ion|iterative", v, re.I):
                    return re.sub(r"\s+", " ", v)[:46]
    j = META / f"{dataset}.summary.json"
    if j.exists():
        try:
            t = json.loads(j.read_text()).get("study_title", "")
            if re.search(r"DDA|DIA|SWATH|MSe", t, re.I):
                return t[:46]
        except (OSError, json.JSONDecodeError):
            pass
    return ""


def spread(files: list[Path], n: int) -> list[Path]:
    """Evenly spaced across the injection list, so a QC block at one end cannot dominate."""
    if len(files) <= n:
        return files
    step = len(files) / n
    return [files[int(i * step)] for i in range(n)]


def main() -> None:
    from collections import Counter
    print(f"{'dataset':<15}{'pol':<5}{'files':>6}   acquisition mix (files per verdict)")
    print("-" * 88)
    for d in sorted(x.name for x in STAGED.iterdir() if x.is_dir()):
        any_files = False
        for pol in ("Pos", "Neg"):
            files = sorted((STAGED / d / pol).glob("*.mzML"))
            if not files:
                continue
            any_files = True
            seen, details = Counter(), {}
            for f in spread(files, SAMPLE):
                try:
                    st = scan(f)
                except OSError:
                    seen["unreadable"] += 1
                    continue
                v, why = verdict(st)
                seen[v] += 1
                details.setdefault(v, why)
            mix = "  ".join(f"{v} x{n} ({details.get(v, '')})" for v, n in seen.most_common())
            usable = "  " if "DDA" in seen else " ✗"
            print(f"{d:<15}{pol:<5}{len(files):>6}{usable} {mix}")
        if not any_files:
            print(f"{d:<15}{'—':<5}{'':>6}   no mzML")
        c = claimed(d)
        if c:
            print(f"{'':<15}{'':<5}{'':>6}   deposit claims: {c}")


if __name__ == "__main__":
    main()
