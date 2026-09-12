"""Blank handling: sample roles and the blank filter.

Neither Compound Discoverer nor LipiDex does this, so nothing here is a reimplementation — it
is a decision, and the reasoning matters more than the code.

**This is a filter, not a subtraction.** Subtracting blank intensities from sample intensities
is the intuitive reading and the wrong operation: it produces negatives you then have to clamp,
it distorts the variance structure every downstream test depends on, and it asserts something
false — that a feature at 60% of blank level is 40% real signal. Contamination and analyte do
not add linearly through one peak near the noise floor.

So the decision is per feature and binary. Either a compound group is credibly above the blanks
and its measured intensities pass through *untouched*, or it is not and the whole row goes.

Two things about blanks that change what the filter means, and which the filter cannot detect
for you:

  * **A blank at the start of a sequence and a blank at the end measure different things.** The
    first measures contamination in the solvent and system; the last measures carryover from
    everything injected before it, which on reversed phase is much larger and is worst for the
    lipids that bind hardest — TG especially. Filtering against a trailing blank at the same
    multiplier will strip real signal. `carryover_report` exists to make that visible before
    you choose which blanks to filter against.
  * **A table whose blanks were already subtracted upstream will be filtered twice**, silently
    deleting real signal. Check the provenance before enabling this.
"""
from __future__ import annotations

import re
from dataclasses import dataclass
from statistics import mean

SAMPLE = "sample"
BLANK = "blank"
QC = "qc"
# A mix of authentic standards, injected to judge the instrument rather than the study. It is not
# a pooled QC: a pool is made from these samples and says whether *this batch* was repeatable,
# while a standard mix is the same material in every batch and every year, so its absolute
# response is comparable across studies and tracks the instrument over time. Both are needed, and
# neither can stand in for the other.
STANDARD = "standard"

_BLANK_PATTERN = re.compile(r"blank|blk", re.I)
_QC_PATTERN = re.compile(r"(^|[_\-])qc([_\-]|\d|$)|pool", re.I)
# `Std_Mix_Pos_01`, `StdMix`, `Standard_02`, `SPLASH`. Checked after blank and QC so a file named
# `Blank_Std_01` is still a blank — the more specific claim wins.
_STANDARD_PATTERN = re.compile(r"std[_\-]?mix|stdmix|(^|[_\-])std([_\-]|\d|$)|standard|splash",
                               re.I)

# Xcalibur `Sample Type`, the instrument's own record of what was injected. Checked ALONGSIDE the
# file name rather than instead of it: the batch generator sets `Std Bracket` going forward, but
# sequences already acquired carry `Unknown` for the same injections and are identified by name.
# Either source may declare a role, so a mix is found whichever way it was recorded.
_BLANK_TYPE = re.compile(r"blank", re.I)
_QC_TYPE = re.compile(r"\bqc\b|quality\s*control", re.I)
_STANDARD_TYPE = re.compile(r"std\s*[_\-]?\s*bracket|standard", re.I)
# A `Sample ID` that declares itself an injection number rather than a vial position.
_INJECTION_ID = re.compile(r"inj[_\-]?\d+", re.I)


def infer_role(filename: str, sample_type: str = "") -> str:
    """A file's role, from its name and its Xcalibur `Sample Type`.

    Either source may declare a role and the more specific claim wins, in the order blank, QC,
    standard, sample. Two independent rules rather than one because they fail in different
    situations: a sequence may be missing or may predate the convention (`Unknown` where a later
    batch writes `Std Bracket`), and a file may be renamed after acquisition. Always overridable
    with `sample_roles` in the run configuration.
    """
    stem = filename.rsplit("/", 1)[-1]
    declared = (sample_type or "").strip()
    if _BLANK_PATTERN.search(stem) or _BLANK_TYPE.search(declared):
        return BLANK
    if _QC_PATTERN.search(stem) or _QC_TYPE.search(declared):
        return QC
    if _STANDARD_PATTERN.search(stem) or _STANDARD_TYPE.search(declared):
        return STANDARD
    return SAMPLE


def infer_roles(filenames, sample_types: "dict[str, str] | None" = None) -> dict[str, str]:
    types = sample_types or {}
    return {name: infer_role(name, types.get(name, "")) for name in filenames}


def read_sample_types(path) -> dict[str, str]:
    """`File Name` -> `Sample Type` from an Xcalibur sequence CSV.

    The header is not on the first line — Xcalibur writes `Bracket Type=4` above it — so the file
    is scanned for the row that names the columns rather than handed to a reader as-is.
    """
    import csv
    from pathlib import Path

    lines = Path(path).read_text(errors="ignore").splitlines()
    start = next((i for i, line in enumerate(lines)
                  if "sample type" in line.lower() and "file name" in line.lower()), None)
    if start is None:
        return {}
    out: dict[str, str] = {}
    for row in csv.DictReader(lines[start:]):
        name = (row.get("File Name") or "").strip()
        if name and name not in out:
            out[name] = (row.get("Sample Type") or "").strip()
    return out


@dataclass(slots=True)
class BlankFilter:
    """Keep a compound group only if it stands clear of the blanks.

    `multiplier` is the fold the samples must exceed the blanks by. 3 is permissive, 5 is the
    usual choice, 10 is strict. `statistic` picks how the blanks are summarised: `mean` is the
    ordinary reading, `max` is stricter and is the right one when a single blank is suspected
    of carryover.
    """

    enabled: bool = True
    multiplier: float = 5.0
    statistic: str = "mean"          # mean | max
    require_detected_in: int = 2     # samples the group must be present in, blanks excluded
    drop_blank_columns: bool = False  # remove blank columns from the written table

    def blank_level(self, values: list[float]) -> float:
        if not values:
            return 0.0
        return max(values) if self.statistic == "max" else mean(values)


@dataclass(slots=True)
class BlankFilterReport:
    kept: int = 0
    dropped: int = 0
    dropped_identified: int = 0
    dropped_names: list[str] = None

    def __post_init__(self):
        if self.dropped_names is None:
            self.dropped_names = []

    def __str__(self) -> str:
        return (f"blank filter: kept {self.kept}, dropped {self.dropped} "
                f"({self.dropped_identified} of them identified)")


def apply_blank_filter(groups, samples, config: BlankFilter) -> BlankFilterReport:
    """Mark compound groups that fail the blank test. Intensities are never modified."""
    report = BlankFilterReport()
    if not config.enabled:
        report.kept = len(groups)
        return report

    sample_columns = [i for i, s in enumerate(samples) if s.role == SAMPLE]
    blank_columns = [i for i, s in enumerate(samples) if s.role == BLANK]
    if not blank_columns or not sample_columns:
        report.kept = len(groups)
        return report

    for group in groups:
        if not group.keep:
            continue
        areas = group.areas
        sample_values = [areas[i] for i in sample_columns if i < len(areas)]
        blank_values = [areas[i] for i in blank_columns if i < len(areas)]

        detected = sum(1 for v in sample_values if v > 0)
        blank_level = config.blank_level(blank_values)
        sample_level = mean(sample_values) if sample_values else 0.0

        passes = detected >= config.require_detected_in and (
            blank_level <= 0.0 or sample_level >= config.multiplier * blank_level)

        if passes:
            report.kept += 1
        else:
            group.keep = False
            group.filter_reason = "Blank"
            report.dropped += 1
            name = group.identification()[0]
            if name:
                report.dropped_identified += 1
                report.dropped_names.append(name)
    return report


def carryover_report(groups, samples, order: dict[str, int] | None = None) -> list[dict]:
    """Compare blanks against each other, so a trailing blank's carryover is visible.

    Returns one row per blank with its total signal and how much of it sits in groups that are
    also strong in samples — which is what carryover looks like, as opposed to contamination.
    """
    blanks = [(i, s) for i, s in enumerate(samples) if s.role == BLANK]
    sample_columns = [i for i, s in enumerate(samples) if s.role == SAMPLE]
    rows = []
    for index, blank in blanks:
        total = 0.0
        shared = 0.0
        detected = 0
        for group in groups:
            if index >= len(group.areas):
                continue
            value = group.areas[index]
            if value <= 0:
                continue
            detected += 1
            total += value
            sample_values = [group.areas[i] for i in sample_columns if i < len(group.areas)]
            if sample_values and mean(sample_values) > value:
                shared += value
        rows.append({
            "file": blank.file,
            "injection": (order or {}).get(blank.file, None),
            "groups detected": detected,
            "total signal": total,
            "fraction shared with samples": (shared / total) if total else 0.0,
        })
    return rows


def read_sequence(path) -> dict[str, int]:
    """Injection order from an Xcalibur sequence CSV: file name -> injection number.

    The order is what tells contamination from carryover, and it is not recoverable from the
    files themselves.

    **The order comes from the row order**, because that is what an Xcalibur sequence means: the
    instrument runs it top to bottom. An earlier version took the trailing digits of `Sample ID`
    instead, which was right for exactly one of the three sequences seen so far:

        CKD rat heart   Sample ID = Inj_01, Inj_02, Inj_51    injection numbers, correct
        skin organoid   Sample ID = RA1, RA2, RA3             VIAL POSITIONS
        Kiterie Lumos   Sample ID = GE8, RH5, RE12            VIAL POSITIONS

    On the Kiterie sequence that gave all five blanks injection 8 — they share a vial, so they
    share a position — and the carryover report became meaningless.

    So `Sample ID` is used only when every value declares itself an injection number (`Inj_01`),
    which also keeps a hand-trimmed subset of a sequence working; otherwise row order wins.

    ⚠ **A name may cover BOTH polarities.** One study's sequence carried 152 distinct `File Name`
    values across 304 rows: the same name for the positive and negative injection of a sample, with
    only the `Path` column telling them apart, while the acquired files are `<name>_Pos.raw` and
    `<name>_Neg.raw`. Keeping the first occurrence then gave the negative injection the positive
    one's number, and — worse — nothing matched the result columns at all, so injection order came
    back empty and every section that needs it went blank without saying why.

    So a polarity-qualified key is emitted alongside the bare one whenever `Path` disambiguates.
    The bare key keeps its first occurrence, as before, for sequences that name one polarity only.
    """
    import csv
    from pathlib import Path

    order: dict[str, int] = {}
    with Path(path).open(newline="", errors="replace") as fh:
        rows = list(csv.reader(fh))
    header_index = None
    for i, row in enumerate(rows):
        if "File Name" in row:
            header_index = i
            break
    if header_index is None:
        return order

    header = rows[header_index]
    name_column = header.index("File Name")
    id_column = header.index("Sample ID") if "Sample ID" in header else None
    path_column = header.index("Path") if "Path" in header else None

    body = [row for row in rows[header_index + 1:]
            if len(row) > name_column and row[name_column].strip()]

    # Use `Sample ID` only when it declares itself an injection number. `Inj_01` says what it is;
    # `GE8` is a vial position wearing the same column. Anything else falls back to row order,
    # which is what an Xcalibur sequence means anyway.
    explicit = None
    if id_column is not None:
        values = [row[id_column].strip() for row in body if len(row) > id_column]
        if values and all(_INJECTION_ID.fullmatch(v) for v in values):
            explicit = True

    for injection, row in enumerate(body, start=1):
        name = row[name_column].strip()
        if explicit:
            digits = re.findall(r"\d+", row[id_column])
            if digits:
                order.setdefault(name, int(digits[-1]))
                continue
        order.setdefault(name, injection)
        if path_column is not None and len(row) > path_column:
            where = row[path_column].lower()
            polarity = "Neg" if "neg" in where else ("Pos" if "pos" in where else "")
            if polarity:
                order.setdefault(f"{name}_{polarity}", injection)
    return order
