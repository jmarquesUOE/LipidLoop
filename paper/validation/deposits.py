"""Fetch and cache each deposit's own metadata, so the table can be rebuilt from nothing.

The scripts here used to read a cache under /tmp that a previous session happened to have filled.
That made every number in the validation table depend on a directory with a session lifetime: the
table could be rebuilt on this machine today and nowhere else, ever. For a manuscript that is not a
reproducibility nicety, it is the difference between a result and an anecdote.

So the fetch lives here, the cache is a real directory next to the scripts, and a missing cache is
filled rather than fatal.

Three repositories, three shapes:

  Metabolomics Workbench   REST: summary, analysis, metabolites, and the full mwTab record
  MetaboLights             web service for the study, FTP for the m_*_maf.tsv identification table
  MassIVE                  PROXI, which carries species and keywords and no method at all

⚠ Network failure must never silently produce a DIFFERENT table. A fetch that fails leaves the
cache untouched and the field renders blank with `not stated` as its source — never a plausible
value, never a stale one presented as fresh.
"""
from __future__ import annotations

import csv
import json
import os
import re
import subprocess
from pathlib import Path

CACHE = Path(os.environ.get("VALIDATION_CACHE",
                            Path(__file__).resolve().parent / ".cache"))
TIMEOUT = 120

# Deposits whose files were downloaded by hand rather than by a repository API.
LOCAL_DEPOSITS = Path("/mnt/datastore/Jair/claudecode/Lipidomics_automation/Validation_Datasets")

WORKBENCH = "https://www.metabolomicsworkbench.org/rest/study/study_id"
METABOLIGHTS_WS = "https://www.ebi.ac.uk/metabolights/ws/studies/public/study"
METABOLIGHTS_FTP = "https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public"
MASSIVE = "https://massive.ucsd.edu/ProteoSAFe/proxi/v0.1/datasets?filter="


def _get(url: str, dest: Path) -> bool:
    """Download to `dest` unless it already exists. False when it could not be fetched.

    Deliberately leaves an existing cache file alone on failure: a half-written or empty file would
    be read later as "the deposit says nothing", which is a different and much worse claim than
    "we could not reach the server".
    """
    if dest.exists() and dest.stat().st_size > 0:
        return True
    dest.parent.mkdir(parents=True, exist_ok=True)
    tmp = dest.with_suffix(dest.suffix + ".part")
    try:
        subprocess.run(["curl", "-sS", "--fail", "--max-time", str(TIMEOUT), url, "-o", str(tmp)],
                       check=True, capture_output=True)
    except (subprocess.CalledProcessError, FileNotFoundError, OSError):
        tmp.unlink(missing_ok=True)
        return False
    if tmp.stat().st_size == 0:
        tmp.unlink(missing_ok=True)
        return False
    tmp.rename(dest)
    return True


def ensure(dataset: str) -> None:
    """Fetch whatever this deposit's repository offers. Safe to call repeatedly and offline."""
    if dataset.startswith("ST"):
        _get(f"{WORKBENCH}/{dataset}/summary", CACHE / f"{dataset}.summary.json")
        _get(f"{WORKBENCH}/{dataset}/analysis", CACHE / f"{dataset}.analysis.json")
        _get(f"{WORKBENCH}/{dataset}/metabolites", CACHE / f"{dataset}.mets.json")
        _get(f"{WORKBENCH}/{dataset}/mwtab/txt", CACHE / f"{dataset}.mwtab.txt")
    elif dataset.startswith("MTBLS"):
        _get(f"{METABOLIGHTS_WS}/{dataset}", CACHE / f"{dataset}.ws.json")
        listing = CACHE / f"{dataset}.listing.html"
        if _get(f"{METABOLIGHTS_FTP}/{dataset}/", listing):
            names = set(re.findall(r'href="(m_[^"]*maf\.tsv)"', listing.read_text(errors="replace")))
            for i, name in enumerate(sorted(names)):
                _get(f"{METABOLIGHTS_FTP}/{dataset}/{name}", CACHE / f"{dataset}.{i}.maf.tsv")
    elif dataset.startswith("MSV"):
        _get(f"{MASSIVE}{dataset}", CACHE / f"{dataset}.json")
    elif dataset.startswith("MTBKS"):
        # MetaboLights *Kickstarter* deposits are not served by the study web service, so nothing
        # here is fetched — the MAFs ship with the download and are already on disk.
        #
        # ⚠ This deposit went unscored for weeks reading as "no names deposited", which is the one
        # failure mode a prefix dispatch hides: an unknown prefix matches no branch, ensure()
        # returns silently, and the absence of names is indistinguishable from a deposit that
        # published none. MTBKS222 in fact publishes 16 MAFs — one per vendor per polarity, with
        # retention times — which makes it the most informative reference in the whole set.
        for maf in sorted(LOCAL_DEPOSITS.glob(f"{dataset}/{dataset}.maf.*.txt")):
            arm = maf.name.split(".maf.", 1)[1].removesuffix(".txt")
            dest = CACHE / f"{dataset}.{arm}.maf.tsv"
            if not dest.exists():
                CACHE.mkdir(parents=True, exist_ok=True)
                dest.write_text(maf.read_text(errors="replace"))


def mwtab_fields(dataset: str) -> dict:
    """mwTab as {key: joined distinct values}. Empty when the deposit is not on Workbench."""
    p = CACHE / f"{dataset}.mwtab.txt"
    if not p.exists():
        return {}
    out: dict[str, list[str]] = {}
    for line in p.open(errors="replace"):
        if "\t" not in line:
            continue
        k, _, v = line.partition("\t")
        out.setdefault(k.strip(), []).append(v.strip())
    return {k: " | ".join(dict.fromkeys(v)) for k, v in out.items()}


def summary(dataset: str) -> dict:
    """The deposit's own summary record, list-unwrapped. {} when absent.

    ⚠ Workbench returns a bare object for some studies and a one-element LIST for others, and the
    difference is invisible until a caller does `summary.get(...)` and gets AttributeError on a
    list. `record()` below already unwraps; this did not, and crashed the table builder.
    """
    try:
        d = json.loads((CACHE / f"{dataset}.summary.json").read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    if isinstance(d, list):
        return d[0] if d and isinstance(d[0], dict) else {}
    return d if isinstance(d, dict) else {}


def record(dataset: str, suffix: str) -> dict:
    """Any cached JSON record, list-unwrapped. {} when absent or unparseable."""
    try:
        d = json.loads((CACHE / f"{dataset}.{suffix}").read_text())
    except (OSError, json.JSONDecodeError):
        return {}
    if isinstance(d, list):
        return d[0] if d else {}
    return d


def metabolite_rows(dataset: str) -> list[dict]:
    try:
        d = json.loads((CACHE / f"{dataset}.mets.json").read_text())
    except (OSError, json.JSONDecodeError):
        return []
    rows = list(d.values()) if isinstance(d, dict) else d
    return [r for r in rows if isinstance(r, dict)]


def maf_rows(dataset: str, arm: str = "") -> list[dict]:
    """MAF rows for a deposit, optionally restricted to one instrument arm.

    ⚠ `arm` matters whenever a deposit publishes per-instrument lists. MTBKS222 ran the SAME
    NIST SRM1950 plasma on five vendors and published a MAF for each; only the QExactive Plus arm
    has been processed here. Pooling all five would score our Thermo results against species only
    the Sciex and Bruker instruments reported, deflating agreement for a reason that has nothing
    to do with the software.
    """
    rows = []
    for maf in sorted(CACHE.glob(f"{dataset}*.maf.tsv")):
        if arm and arm.lower() not in maf.name.lower():
            continue
        with maf.open(errors="replace") as fh:
            rows.extend(csv.DictReader(fh, delimiter="\t"))
    return rows


def seed_from(old: Path) -> int:
    """Copy an existing cache in, so a first run does not re-download what is already on disk."""
    if not old.is_dir():
        return 0
    CACHE.mkdir(parents=True, exist_ok=True)
    n = 0
    for f in old.iterdir():
        target = CACHE / f.name
        if f.is_file() and f.stat().st_size > 0 and not target.exists():
            target.write_bytes(f.read_bytes())
            n += 1
    return n


if __name__ == "__main__":
    import sys
    seeded = seed_from(Path("/tmp/claude-1000/-home-jair/"
                            "1e5bf05d-c2a1-4067-b572-5355f626c8dd/scratchpad/meta"))
    if seeded:
        print(f"  seeded {seeded} file(s) from the previous session's cache")
    for d in sys.argv[1:] or sorted(p.name for p in Path("/home/jair/validation_staged").iterdir()
                                    if p.is_dir()):
        ensure(d)
        have = len(list(CACHE.glob(f"{d}.*")))
        print(f"  {d:<15} {have} cached file(s)")
