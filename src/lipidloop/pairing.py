"""Which positive-mode injection is the same vial as which negative-mode one, read from the data (v1.1 item 20).

The combined table needs the two polarities' columns paired. Pairing on file names failed on the
first Lumos sequence named `QC_01_P` / `QC_01_N` (item 19), and a name says only what someone
typed. The instrument records what it did:

* **polarity** in every scan (mzML `MS:1000130` positive / `MS:1000129` negative). Read here to
  CHECK the folder, never to move a file: a file in Pos/ whose scans are negative, or a file with
  both, is reported, and the run goes on. Some public deposits switch polarity inside one file,
  so stopping would break runs that work.
* **vial** in the Thermo raw header (ThermoRawFileParser, header only, about a second per file),
  else the sequence's `Position` column.
* **acquisition time** in the mzML `startTimeStamp`, which orders repeat injections of one vial.

Pairing: same vial, repeats in acquisition order. Without a vial for every injection (mzML-only
public data), or when the vials do not tell injections apart (every file in one vial), the
item-19 name harmonisation is used. Either way a mismatch is refused, never padded.

Only the combined table uses this. Identification and every per-polarity table are untouched.
"""
from __future__ import annotations

import csv
import json
import re
import subprocess
import tempfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

POSITIVE, NEGATIVE = "MS:1000130", "MS:1000129"
HEAD_BYTES = 8_000_000          # ~900 scans of a Lumos file; enough to see both terms in a switching file
TABLE = "Polarity_Pairing.csv"


@dataclass
class Injection:
    polarity: str                   # the folder it came from: Pos or Neg
    column: str                     # its column in the result table (the file stem)
    vial: str = ""
    vial_source: str = ""           # "raw header" or "sequence"
    started: str = ""               # mzML startTimeStamp, ISO 8601
    scans: str = "unknown"          # Pos, Neg, mixed or unknown, from the scans themselves


@dataclass
class Pair:
    column: str                     # the combined table's column
    pos: Injection
    neg: Injection
    matched_by: str


def scan_polarity(mzml: str | Path) -> tuple[str, str]:
    """(polarity of the scans, startTimeStamp) from the head of an mzML file."""
    try:
        with open(mzml, "rb") as fh:
            head = fh.read(HEAD_BYTES).decode("latin-1")
    except OSError:
        return "unknown", ""
    stamp = re.search(r'<run\b[^>]*\bstartTimeStamp="([^"]+)"', head)
    pos, neg = POSITIVE in head, NEGATIVE in head
    scans = "mixed" if pos and neg else "Pos" if pos else "Neg" if neg else "unknown"
    return scans, stamp.group(1) if stamp else ""


def raw_vial(raw: str | Path, parser: str | Path) -> str:
    """Vial position from a Thermo raw header; '' when it cannot be read."""
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "meta.json"
        try:
            subprocess.run([str(parser), "-i", str(raw), "-f", "4", "-m", "0", "-c", str(out)],
                           capture_output=True, text=True, timeout=120)
            meta = json.loads(out.read_text())
        except (OSError, ValueError, subprocess.SubprocessError):
            return ""
    for entry in meta.get("SampleData", []):
        if entry.get("name") == "Vial":
            return str(entry.get("value") or "").strip()
    return ""


def sequence_positions(sequence: str | Path) -> dict[str, str]:
    """File Name -> Position from an Xcalibur sequence export ('Bracket Type=' line first)."""
    try:
        lines = Path(sequence).read_text(encoding="utf-8-sig", errors="replace").splitlines()
    except OSError:
        return {}
    start = next((i for i, line in enumerate(lines) if "File Name" in line), None)
    if start is None:
        return {}
    out = {}
    for row in csv.DictReader(lines[start:]):
        name, position = (row.get("File Name") or "").strip(), (row.get("Position") or "").strip()
        if name and position:
            out[Path(name.replace("\\", "/")).stem] = position
    return out


def _by_name(pos: list[Injection], neg: list[Injection], strip) -> list[Pair]:
    neg_by = {strip(i.column, "Neg"): i for i in neg}
    pos_by = {strip(i.column, "Pos"): i for i in pos}
    if len(neg_by) != len(neg) or len(pos_by) != len(pos) or set(pos_by) != set(neg_by):
        raise ValueError(
            f"sample columns differ between polarities after removing the suffix:\n"
            f"  {sorted(pos_by)}\n  {sorted(neg_by)}\nThe two runs did not measure the same vials.")
    return [Pair(name, i, neg_by[name], "file name") for name, i in pos_by.items()]


def _informative(injections: list[Injection]) -> bool:
    """Vials pair injections only if every one has a vial and they are not all the same vial."""
    vials = [i.vial for i in injections]
    return all(vials) and (len(set(vials)) > 1 or len(injections) <= 2)


def pair(pos: list[Injection], neg: list[Injection], strip) -> list[Pair]:
    """Pair the two polarities' injections, in positive-mode column order.

    `strip` is the name harmoniser (`polarity._strip_polarity`), used for the fallback and to name
    a pair's column. Raises ValueError when the injections cannot all be paired one to one.
    """
    if not _informative(pos + neg):
        return _by_name(pos, neg, strip)
    groups: dict[str, dict[str, list[Injection]]] = defaultdict(lambda: {"Pos": [], "Neg": []})
    for i in pos + neg:
        groups[i.vial][i.polarity].append(i)
    unequal = {v: (len(g["Pos"]), len(g["Neg"])) for v, g in groups.items()
               if len(g["Pos"]) != len(g["Neg"])}
    if unequal:
        detail = ", ".join(f"{v} {p} Pos / {n} Neg" for v, (p, n) in sorted(unequal.items()))
        raise ValueError(f"injections per vial differ between polarities ({detail}). "
                         f"The two runs did not measure the same vials.")
    order = {i.column: k for k, i in enumerate(pos)}
    pairs = []
    for g in groups.values():
        # Repeats of one vial pair in acquisition order; the table's column order breaks ties
        # and stands in when a timestamp is missing.
        for p, n in zip(sorted(g["Pos"], key=lambda i: (i.started, order[i.column])),
                        sorted(g["Neg"], key=lambda i: (i.started, neg.index(i)))):
            sp, sn = strip(p.column, "Pos"), strip(n.column, "Neg")
            source = p.vial_source if p.vial_source == n.vial_source else f"{p.vial_source} / {n.vial_source}"
            pairs.append(Pair(sp if sp == sn else p.column, p, n, f"vial ({source})"))
    pairs.sort(key=lambda x: order[x.pos.column])
    names = [x.column for x in pairs]
    if len(set(names)) != len(names):
        raise ValueError(f"two pairs would share a combined column name: {names}")
    return pairs


def disagreements(pairs: list[Pair], strip) -> list[Pair]:
    """Vial pairs whose file names say otherwise — reported, because one of the two is wrong."""
    return [x for x in pairs if x.matched_by != "file name"
            and strip(x.pos.column, "Pos") != strip(x.neg.column, "Neg")]


def write(path: str | Path, pairs: list[Pair]) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["Column", "Pos file", "Neg file", "Pos vial", "Neg vial", "Pos started",
                    "Neg started", "Pos scans", "Neg scans", "Matched by"])
        for x in pairs:
            w.writerow([x.column, x.pos.column, x.neg.column, x.pos.vial, x.neg.vial,
                        x.pos.started, x.neg.started, x.pos.scans, x.neg.scans, x.matched_by])


def wrong_folder(injections: list[Injection]) -> list[Injection]:
    """Files whose scans are not (only) the polarity of the folder they sit in."""
    return [i for i in injections if i.scans in ("mixed", "Pos", "Neg") and i.scans != i.polarity]
