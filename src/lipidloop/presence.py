"""Presence filtering — deciding which features are measured often enough to analyse.

A feature detected in a handful of injections is not a measurement, it is a detection. Filtering
on presence is how that gets decided, and the parameter that matters is not the threshold but the
**group the threshold is applied within**. On the brain set at 70%, the identified rows kept are:

    group = genotype        (2 x 20)   534
    group = age x genotype  (4 x 10)   585
    group = age x sex x geno (8 x 5)   626

a wider spread than moving the threshold from 50% to 100% inside any one grouping. The rule is
"detected in at least `min_fraction` of one group", not "of all samples", because a lipid present
in one group and absent from another is the strongest result the experiment can produce and an
all-samples rule deletes it.

**Why not filter on the QCs instead.** Requiring presence in every pooled QC injection looks
stricter and cleaner — 2.8% residual missingness against 13.0%. It is also structurally biased
against the result you are looking for. A pooled QC mixes every sample, so a lipid confined to one
group of ten is diluted tenfold in the pool; on this data the rows that pass the group rule and
fail the QC rule are 4-6x fainter than average, and a quarter of them are present in over 80% of
one group and entirely absent from another. The QC rule does not drop those for being unreliable.
It drops them for being group-specific.

## What happens to the gaps

Nothing, here. The rows this keeps still carry zeros, and filling them is `gapfill.py`'s job — by
re-integrating the raw signal, not by substituting a statistic.

That division matters for one claim in particular. This filter deliberately keeps a lipid present
in one group and absent from another, and the obvious next step is to protect that whole-group
zero from being filled. But "protect it" assumes the absence is real rather than a detection
failure, and only going back to the raw file can tell those apart. So the assumption is not made
here: `gapfill.py` measures those positions and reports how many turned out to carry a peak.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path


@dataclass
class PresenceFilter:
    """Keep a feature detected in at least `min_fraction` of one group.

    `group_by` names columns in `metadata`, a CSV whose first column is the injection file name.
    With no metadata the samples form a single group and the rule degrades to "detected in
    `min_fraction` of all samples", which is stated in the log rather than done silently.

    Filling the gaps is a separate stage — see `gapfill.py`. It re-integrates the raw signal
    rather than substituting a statistic, which is also what tests the assumption below.
    """

    enabled: bool = False
    #: ⚠ 0.6, not the 0.7 this was written with. The change is about what the threshold MEANS at
    #: small n, which is most of this validation set. At 0.7 a group of three needs all three,
    #: because 2/3 is 0.667 and misses — a triplicate arm was being filtered harder than a
    #: twenty-injection one under the same nominal rule. At 0.6 a triplicate needs two of three,
    #: which is what "detected in most of the group" was always meant to say.
    #:
    #: It also removes a special case: ST000991_DDA is one pooled QC across three precursor windows,
    #: where a lipid at m/z 900 can appear in at most 3 of 5 files. 0.6 is the highest threshold no
    #: window-limited lipid can fail on coverage alone, and it was being set per-study to get there.
    min_fraction: float = 0.6
    # A group is a group at n >= 3. Below that the fraction rule barely constrains anything: at
    # n=1 it is vacuous (1/1 is 100%, so a single detection passes) and at n=2 it is a coin toss.
    # Cells smaller than this take no part; if that leaves none, the filter reports it and skips
    # rather than passing everything on a technicality.
    #
    # ⚠ This floor and `min_fraction` interact, and the interaction is why the default is 0.6. At
    # 0.7 a group of three needed ALL THREE (2/3 = 0.667 misses), so a triplicate arm was filtered
    # harder than a twenty-injection one under the same nominal rule. At 0.6 it needs two of three.
    #
    # ⚠ Two validation deposits still fall below it and are meant to: ST003077 is a pooled-tissue
    # atlas with n=1 per tissue, and ST000991_DDA has 2/2/1 injections across three precursor
    # windows. Neither has repeated observations to filter on, and the log says so.
    min_group_size: int = 3
    metadata: str = ""
    group_by: list = field(default_factory=list)


@dataclass
class PresenceReport:
    kept: int = 0
    dropped: int = 0
    identified_kept: int = 0
    identified_dropped: int = 0
    groups: dict = field(default_factory=dict)
    ungrouped: bool = False

    def __str__(self) -> str:
        shape = ", ".join(f"{k}={v}" for k, v in sorted(self.groups.items()))
        where = (f"groups: {shape}" if not self.ungrouped
                 else f"no metadata — grouped by role: {shape}")
        return (f"presence filter: kept {self.kept}, dropped {self.dropped} "
                f"({self.identified_dropped} of them identified); {where}")


def read_metadata(path: str | Path, group_by: list) -> dict[str, tuple]:
    """`{injection file: group key}` from a CSV whose first column names the injection."""
    rows = list(csv.DictReader(Path(path).open(newline="")))
    if not rows:
        return {}
    key_column = list(rows[0])[0]
    return {r[key_column].strip(): tuple(r.get(g, "").strip() for g in group_by)
            for r in rows if r.get(key_column)}


def _groups(samples, metadata: dict[str, tuple]) -> dict[tuple, list[int]]:
    """The cells the filter works in: one per group, where a group is whatever the study says.

    A group is a biological arm, a lab in a ring trial, a tissue in an atlas, a precursor window in
    a DDA series — or the pooled QCs. They are all just groups, and the rule is a union: a feature
    survives if it was detected reliably in ANY ONE of them.

    ⚠ **Adding QC as a group can only ever rescue a feature, never drop one**, because the rule is
    "detected in at least one group". That is the opposite of the QC rule this module warns about at
    the top, which REQUIRES presence in the pooled QCs and does drop group-specific lipids for being
    group-specific. Requiring is biased; offering is not. What QC-as-a-group buys is the study whose
    every injection is a pool — a NIST replicate series, a pooled-tissue atlas, a DDA window series
    — which previously produced no cells at all, so the filter was skipped and those tables alone
    carried unfiltered noise while every other study in the same comparison had been filtered.

    Blanks and standard/IS injections are not groups: a blank is the thing being filtered against,
    and a feature present in every Std_Mix injection is the standard, not evidence about the study.

    With no metadata the ROLE is the group, so pools and samples form separate cells rather than one
    mixed cell — which is what makes the all-pool study work without any metadata at all.
    """
    out: dict[tuple, list[int]] = {}
    for index, sample in enumerate(samples):
        if sample.role in ("blank", "standard", "istd"):
            continue
        key = metadata.get(sample.file)
        if key is None:
            key = (sample.role or "sample",)
        out.setdefault(key, []).append(index)
    return out


def _sized(cells: dict, minimum: int) -> dict:
    """Drop cells too small for the fraction rule to mean anything. If that leaves none, the caller
    reports the filter as skipped rather than passing everything on a technicality."""
    return {k: v for k, v in cells.items() if len(v) >= max(1, minimum)}


def presence_cells(samples, params: "PresenceFilter") -> dict[tuple, list[int]]:
    """The sample groups the filter works in. Shared with gap filling, which has to ask the same
    question about the same cells — two different groupings would make the absence check
    meaningless."""
    metadata: dict[str, tuple] = {}
    if params.metadata and params.group_by:
        metadata = read_metadata(params.metadata, params.group_by)
    return _sized(_groups(samples, metadata), params.min_group_size)


def apply_presence_filter(groups_of_compounds, samples, params: PresenceFilter,
                          log=None) -> PresenceReport:
    """Drop any compound group not detected in `min_fraction` of at least one sample group."""
    say = log or (lambda _: None)
    report = PresenceReport()
    metadata: dict[str, tuple] = {}
    if params.metadata and params.group_by:
        metadata = read_metadata(params.metadata, params.group_by)
    cells = _sized(_groups(samples, metadata), params.min_group_size)
    report.ungrouped = not metadata
    report.groups = {"/".join(k) if k != ("all",) else "all": len(v) for k, v in cells.items()}
    if not cells:
        say(f"presence filter: no group of at least {params.min_group_size} injections, skipped")
        return report

    for group in groups_of_compounds:
        if not group.keep:
            continue
        identified = bool(group.identification()[0])
        passes = False
        for members in cells.values():
            detected = sum(1 for i in members if group.areas[i] > 0)
            if members and detected / len(members) >= params.min_fraction - 1e-9:
                passes = True
                break
        if passes:
            report.kept += 1
            report.identified_kept += identified
        else:
            group.keep = False
            group.filter_reason = (f"Not detected in {params.min_fraction:.0%} of any group")
            report.dropped += 1
            report.identified_dropped += identified
    say(str(report))
    return report


