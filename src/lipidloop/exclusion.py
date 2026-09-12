"""Find precursors that spend MS2 duty cycle and never return an identification.

A DDA run has a fixed budget. These files fit about nine MS2 into a 0.87 s cycle, which is why
there are only ~1.2 MS1 scans per second and roughly five points across a 4.7 s peak. Every MS2
spent on column bleed is one not spent on a lipid, and it is paid for twice — once in the
identification that was not attempted, once in the MS1 sampling that makes the peak thinner.

On the validation data that budget is being spent badly. Over 24 files per polarity, counting only
precursors never identified in *any* file, eluting across more than half the gradient, and equally
present in the blanks:

    positive    285 precursors    43.0% of all MS2    polysiloxane 24.5%, PEG 2.8%
    negative    343 precursors    67.0% of all MS2    sodium formate 44.6%, PTFE/PFPE 7.4%

The split is almost perfectly complementary — siloxane dominates positive mode and is invisible in
negative, salt clusters dominate negative and are a rounding error in positive. Put these on the
instrument exclusion list and that budget goes back to the sample. See docs/EXCLUSION_LIST.md.

## Why this is safe, and where it stops being safe

Excluding a precursor is irreversible in a way that filtering afterwards is not: the spectrum is
never acquired, so a mistake here cannot be undone in reprocessing. "Never identified" alone is a
bad criterion — most unidentified spectra are real ions the library does not cover, and excluding
them would freeze the method's blind spots in place.

What separates contamination from an unidentified lipid is **chromatography**. A lipid elutes as a
peak, a few seconds wide. Background is present at every retention time because it enters the
source continuously. So the load-bearing criterion is elution span, not identification failure,
and the identification test is only there to make sure we are not excluding something the pipeline
is currently getting right.

Four conditions, all required:

1. **Never identified**, in any file of the batch.
2. **Fragmented at least `MIN_SCANS` times** — below that it is not worth a list entry.
3. **Spans more than `MIN_GRADIENT_FRACTION` of the run.** The one that carries the argument.
4. **Also fragmented in a blank**, when blanks are given. Contamination does not need the sample.

Homologous-series membership is reported but deliberately *not* required: it is corroboration for
a human reading the list, and requiring it would miss single-compound contaminants.

The list is a recommendation to put into the acquisition method. Nothing here changes processing,
and the script writes a file rather than editing anything.
"""
from __future__ import annotations

import bisect
import statistics
from dataclasses import dataclass

MIN_SCANS = 5                    # fewer than this is not worth an exclusion-list entry
MIN_GRADIENT_FRACTION = 0.5      # must be present across at least half the run
CLUSTER_PPM = 10.0               # precursors this close are the same ion
BLANK_TOLERANCE = 0.01

# Repeat units of the polymers that contaminate every LC-MS lab. Reported as corroboration;
# membership is never a criterion, because a contaminant does not have to be a polymer.
SERIES = {
    "polysiloxane": 74.01872,    # C2H6OSi — column bleed, PDMS from tubing and septa
    "PEG": 44.02621,             # C2H4O — detergents, plasticiser, plastic labware
    "PPG": 58.04186,             # C3H6O — antifoam, mould release
    "PTFE/PFPE": 99.99361,       # C2F4 — fluoropolymer tubing, pump oil; negative mode mostly
    "PFPE": 65.99146,            # CF2O — the other perfluoropolyether repeat
    # Salt clusters, not polymers, but they behave identically: continuous, blank-present,
    # unidentifiable. Sodium formate is what a formic-acid mobile phase makes out of sodium
    # leached from glassware, and on this data it is the single largest consumer of negative-mode
    # duty cycle. Worth fixing at the bench as well as excluding.
    "sodium formate": 67.98738,  # CHO2Na — [(HCOONa)n + HCOO]-
    "sodium acetate": 82.00303,  # C2H3O2Na — the same thing on an acetate method
}
SERIES_TOLERANCE = 0.01
MIN_SERIES_MEMBERS = 3

# Unnamed series are discovered as well. Naming every contaminant a lab can have is hopeless —
# the negative-mode data here is dominated by a Δ 67.987 series that matches none of the above —
# and the argument for exclusion never rested on knowing what the compound is. A repeat unit
# holding across several ions is evidence of a synthetic polymer whatever it turns out to be.
UNNAMED_RANGE = (20.0, 200.0)
MIN_UNNAMED_MEMBERS = 4          # stricter than a named repeat: nothing corroborates it

# A repeat of CH2 is what a LIPID series looks like — PC 32:0, 34:0, 36:0 are 14.0157 apart, and
# so are the members of every acyl-chain homologous series in the sample. Discovering that as a
# "polymer" and printing it next to column bleed would be actively misleading, so multiples of
# CH2 are never proposed as an unnamed repeat. None of the named units is one: PEG 44.026 sits
# between 3xCH2 (42.047) and 4xCH2 (56.063).
CH2 = 14.01565
CH2_GUARD = 0.02


@dataclass(slots=True)
class Wasted:
    """One precursor that is costing duty cycle."""

    mz: float
    scans: int                  # MS2 spent on it across the batch
    rt_lo: float
    rt_hi: float
    gradient_fraction: float    # share of the run it is present across
    in_blank: bool
    series: str = ""            # e.g. "polysiloxane" — corroboration, not a criterion

    def __str__(self) -> str:
        tags = [f"{self.scans} MS2",
                f"RT {self.rt_lo:.1f}-{self.rt_hi:.1f} ({100 * self.gradient_fraction:.0f}% of run)"]
        if self.in_blank:
            tags.append("in blank")
        if self.series:
            tags.append(self.series)
        return f"m/z {self.mz:.4f}: " + ", ".join(tags)


def cluster_precursors(scans, ppm: float = CLUSTER_PPM):
    """Group `(key, mz, retention)` triples into distinct precursors.

    The instrument reports the same ion with a few ppm of scatter between scans, so a raw m/z is
    not an identity. Sorted sweep with a relative gap, rather than fixed bins: at m/z 1500 a
    0.01 bin is 7 ppm and would split the siloxane series, and at m/z 100 it is 100 ppm and would
    merge two different ions.
    """
    ordered = sorted(scans, key=lambda s: s[1])
    if not ordered:
        return []
    groups, current = [], [ordered[0]]
    for scan in ordered[1:]:
        if (scan[1] - current[-1][1]) / max(scan[1], 1.0) * 1e6 <= ppm:
            current.append(scan)
        else:
            groups.append(current)
            current = [scan]
    groups.append(current)
    return groups


def assign_series(mz_values: list[float], tolerance: float = SERIES_TOLERANCE,
                  min_members: int = MIN_SERIES_MEMBERS,
                  discover: bool = True) -> dict[float, str]:
    """Label each m/z with the polymer series it belongs to, where it belongs to one.

    A member needs at least `min_members` in its chain — two ions separated by 74.019 is a
    coincidence, six of them is column bleed. Named repeats are matched first; with `discover`
    the leftovers are searched for repeats nobody named, which is what the negative-mode data
    needs.
    """
    ordered = sorted(mz_values)
    labels: dict[float, str] = {}
    for name, repeat in SERIES.items():
        _chain_series(ordered, labels, repeat, name, tolerance, min_members)

    while discover:
        remaining = [mz for mz in ordered if mz not in labels]
        repeat = _commonest_repeat(remaining, tolerance)
        if repeat is None:
            break
        before = len(labels)
        _chain_series(ordered, labels, repeat, f"unnamed Δ{repeat:.3f}", tolerance,
                      MIN_UNNAMED_MEMBERS)
        if len(labels) == before:
            break
    return labels


def _chain_series(ordered, labels, repeat, name, tolerance, min_members) -> None:
    remaining = [mz for mz in ordered if mz not in labels]
    for start in remaining:
        if start in labels:
            continue
        chain, probe = [start], start
        while True:
            nxt = _nearest(remaining, probe + repeat, tolerance)
            if nxt is None or nxt in labels:
                break
            chain.append(nxt)
            probe = nxt
        if len(chain) >= min_members:
            for mz in chain:
                labels[mz] = name


def _commonest_repeat(mz_values: list[float], tolerance: float) -> "float | None":
    """The spacing that recurs most among unlabelled candidates, if any recurs enough.

    Pairwise rather than adjacent: a polymer series is interleaved with everything else in a
    sorted m/z list, so its members are rarely neighbours.
    """
    low, high = UNNAMED_RANGE
    tally: dict[float, list[float]] = {}
    for index, first in enumerate(mz_values):
        for second in mz_values[index + 1:]:
            delta = second - first
            if low <= delta <= high:
                tally.setdefault(round(delta, 2), []).append(delta)
    for key in [k for k in tally if _is_ch2_multiple(k)]:
        del tally[key]
    if not tally:
        return None
    best = max(tally.values(), key=len)
    if len(best) < MIN_UNNAMED_MEMBERS - 1:
        return None
    return statistics.mean(best)


def _is_ch2_multiple(delta: float, tolerance: float = CH2_GUARD) -> bool:
    """True for 14.016, 28.031, 42.047 ... — a lipid homologous series, not a polymer."""
    return abs(delta - round(delta / CH2) * CH2) <= tolerance


def _nearest(ordered: list[float], target: float, tolerance: float) -> "float | None":
    index = bisect.bisect_left(ordered, target - tolerance)
    best, best_error = None, tolerance
    while index < len(ordered) and ordered[index] <= target + tolerance:
        error = abs(ordered[index] - target)
        if error <= best_error:
            best, best_error = ordered[index], error
        index += 1
    return best


def find_wasted(scans, identified: set, blank_precursors: "list[float] | None" = None,
                min_scans: int = MIN_SCANS,
                min_gradient_fraction: float = MIN_GRADIENT_FRACTION,
                require_blank: bool = True) -> list[Wasted]:
    """The exclusion-list candidates.

    `scans` is `(key, precursor_mz, retention_minutes)` for every MS2 in the batch; `identified`
    is the set of keys that produced a library match. Pass `blank_precursors` to require that a
    candidate also fires in a blank — the strongest single piece of evidence available without
    re-acquiring anything.
    """
    scans = list(scans)
    if not scans:
        return []
    retentions = [s[2] for s in scans]
    gradient = max(retentions) - min(retentions)
    if gradient <= 0:
        return []

    blank = sorted(blank_precursors) if blank_precursors else []

    found = []
    for group in cluster_precursors(scans):
        if len(group) < min_scans:
            continue
        if any(s[0] in identified for s in group):
            continue
        rts = [s[2] for s in group]
        fraction = (max(rts) - min(rts)) / gradient
        if fraction < min_gradient_fraction:
            continue
        mz = statistics.mean(s[1] for s in group)
        seen_in_blank = _in_blank(blank, mz)
        if require_blank and blank and not seen_in_blank:
            continue
        found.append(Wasted(mz=mz, scans=len(group), rt_lo=min(rts), rt_hi=max(rts),
                            gradient_fraction=fraction, in_blank=seen_in_blank))

    for candidate, label in assign_series([w.mz for w in found]).items():
        for w in found:
            if w.mz == candidate:
                w.series = label
    found.sort(key=lambda w: -w.scans)
    return found


def _in_blank(blank: list[float], mz: float, tolerance: float = BLANK_TOLERANCE) -> bool:
    if not blank:
        return False
    index = bisect.bisect_left(blank, mz - tolerance)
    return index < len(blank) and blank[index] <= mz + tolerance


def duty_cycle(wasted: list[Wasted], total_ms2: int) -> float:
    """Share of the MS2 budget the candidates are consuming."""
    return sum(w.scans for w in wasted) / max(total_ms2, 1)


def write_exclusion_list(wasted: list[Wasted], path, polarity: str, charge: int = 1) -> None:
    """Write the list in the same shape as `method/CE_inclusion_list.csv`.

    Readable form, carrying the evidence for every entry: how much duty cycle it cost, where it
    eluted, whether the blank had it, which series it belongs to. This is the one to review.
    `write_thermo_list` writes the same rows in the instrument's import format.

    No retention window: these are present throughout, so a window would only create a hole in
    the exclusion where the contaminant is still being fragmented.
    """
    from pathlib import Path

    lines = ["m/z,Charge,Polarity,MS2 scans wasted,RT observed (min),In blank,Series,Note"]
    for w in wasted:
        note = f"{w.series} contamination" if w.series else "background, never identified"
        lines.append(f"{w.mz:.4f},{charge},{polarity},{w.scans},"
                     f"{w.rt_lo:.2f}-{w.rt_hi:.2f},{'yes' if w.in_blank else 'no'},"
                     f'{w.series or ""},"{note}"')
    Path(path).write_text("\n".join(lines) + "\n")


def write_thermo_list(wasted: list[Wasted], path, polarity: str, charge: int = 1,
                      start: "float | None" = None, end: "float | None" = None) -> None:
    """The same rows in the column layout Xcalibur's mass-list import expects.

    Empty `Start`/`End` mean the whole run, which is what a continuously-present contaminant
    needs — a window would leave a hole where the contaminant is fragmented again. Pass `start`
    and `end` only to restrict an entry deliberately.

    How wide an entry really is comes from the method's exclusion mass tolerance, not from this
    file: the instrument skips any precursor within that tolerance of a listed mass, so a real
    lipid close enough is skipped too. That is the one failure mode this format cannot express —
    run `scripts/check_exclusion_safety.py` against the libraries and the identifications before
    importing. (The isolation window is a separate thing and a milder one: a lipid inside it is
    still selected on its own m/z and still fragmented, just co-isolated with the contaminant,
    which was happening before the list existed.)
    """
    from pathlib import Path

    sign = "Positive" if polarity.strip() in {"+", "Positive", "positive"} else "Negative"
    lines = ["Mass [m/z],Formula [M],Species,CS [z],Polarity,Start [min],End [min],"
             "(N)CE,MSX ID,Comment"]
    for w in wasted:
        comment = f"{w.series} contamination" if w.series else "background, never identified"
        lines.append(f"{w.mz:.4f},,,{charge},{sign},"
                     f"{'' if start is None else f'{start:.2f}'},"
                     f"{'' if end is None else f'{end:.2f}'},,,"
                     f'"{comment} ({w.scans} MS2)"')
    Path(path).write_text("\n".join(lines) + "\n")


def series_breakdown(wasted: list[Wasted], total_ms2: int) -> list[tuple]:
    """`(series, precursors, scans, share)` per series, largest first."""
    grouped: dict[str, list[Wasted]] = {}
    for w in wasted:
        grouped.setdefault(w.series or "unassigned", []).append(w)
    rows = [(name, len(members), sum(m.scans for m in members),
             sum(m.scans for m in members) / max(total_ms2, 1))
            for name, members in grouped.items()]
    rows.sort(key=lambda r: -r[2])
    return rows
