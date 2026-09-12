"""Merging one chromatographic peak that was integrated as two compound groups.

A split peak is one molecule whose feature was cut in two, so its signal is divided between two
rows. It is not two isomers, and the difference is measurable rather than a matter of judgement:

    same identification, same adduct, median abs delta m/z  0.0006     37 pairs
    same identification, different adduct                   4.9549      3 pairs

Three orders of magnitude. The second group is one molecule seen as two ions — `[M+H]+` against
`[M+H-H2O]+` differing by 18.0105, water to four decimal places — and those must NOT be summed,
because two ions of one molecule respond differently and adding them adds two scales together.
Only the first group is a split.

**Correlation is the wrong test, and backwards.** The obvious rule — merge rows that correlate
highly — merges nothing: of 36 same-name close pairs on the Kiterie positive run, none reached
rho 0.8 and the median was +0.16, with 9 outright negative. A split peak divides a conserved
total, and the integration boundary moves between injections, so when one half gets more the
other gets less. High correlation is instead the signature of two genuinely distinct co-regulated
species — exactly the pairs that must be left alone.

**The test used here is precision in the pooled QCs.** Every QC vial holds the same material, so
variance between them is technical by construction. If summing two rows makes that variance fall,
the division between them was technical. That is not a proxy for the question, it is the question:

    median QC CV of the worse half   65.4%
    median QC CV of the sum          13.6%

With fewer than `min_pools` QC injections there is no such evidence and nothing is merged — the
peak finder's own retention-only redundancy filter is left to it, and the log says so. Guessing
in the absence of QCs would silently change quantification on exactly the runs least able to
show it had gone wrong.

Runs after the peak finder rather than inside it: the finder reproduces LipiDex, and this is an
addition on top of it, like the blank filter.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from statistics import mean, stdev


@dataclass
class SplitPeakMerge:
    """When two rows are one peak.

    `max_peak_widths` is in units of the run's mean FWHM, so it travels between gradients — a
    fixed 0.2 min window is two peak widths on one method and half a peak width on another.

    `max_ppm` guards the quant ion. With the identification and the adduct already required to
    match it is nearly redundant, which is the point: it costs nothing and it stops a coincidence
    of naming from summing two masses that are not the same ion.

    `min_cv_improvement` is the factor by which summing must beat the worse half in the pooled
    QCs. The observed split pairs improve by about 4.8x and the co-eluting-but-distinct pairs by
    about 1.1x, so 1.5 sits in open ground rather than on a boundary that needs tuning.

    A second signature is tested when the precision test declines, because the two kinds of split
    look nothing alike. The precision test was built for a peak whose apex WANDERS across the
    integration boundary from file to file: one half gains what the other loses, the pair is
    ANTI-correlated, and summing them plainly improves precision. A peak cut at the SAME point in
    every file behaves oppositely — both halves scale with the true abundance, so each is precise
    on its own, the pair is positively correlated at a near-constant ratio, and summing improves
    nothing measurable. `Cer[NDS] d41:2` on the study this was added for: 0.111 min apart against a
    0.089 min peak width, same formate adduct, both halves in all fourteen injections, rho +0.83,
    ratio 1.21 varying by 18.6%, and CVs of 13.8% and 11.3% that summing could not better.

    `systematic_ratio_cv` is the load-bearing one. Two molecules that merely co-vary reach high
    correlation easily — `TG 61:3` on the same run correlates at +0.95 — but their ratio moves
    (50.7% there). A ratio held to within a quarter across samples spanning two biological groups
    is what one peak divided by a fixed boundary looks like, and not what two molecules do unless
    they are locked together.
    """

    enabled: bool = True
    max_peak_widths: float = 2.0
    max_ppm: float = 10.0
    min_cv_improvement: float = 1.5
    min_pools: int = 3
    systematic: bool = True
    systematic_rho: float = 0.8
    systematic_ratio_cv: float = 25.0
    systematic_min_injections: int = 6


@dataclass
class SplitPeakReport:
    merged: int = 0
    systematic: int = 0
    considered: int = 0
    skipped_no_pools: bool = False
    pools: int = 0
    examples: list[tuple[str, float, float]] = field(default_factory=list)

    def __str__(self) -> str:
        if self.skipped_no_pools:
            return (f"split-peak merge: skipped, {self.pools} pooled QC injections "
                    f"(need the QCs to tell a split from two isomers)")
        out = (f"split-peak merge: {self.merged} of {self.considered} candidate pairs merged "
               f"on {self.pools} pooled QCs")
        if self.systematic:
            out += (f" ({self.systematic} of them on a constant ratio across injections rather "
                    f"than on precision — a peak cut at the same point in every file)")
        if self.examples:
            worst = max(self.examples, key=lambda e: e[1] / e[2] if e[2] > 0 else 0)
            out += f" (best: {worst[0]} QC CV {worst[1]:.0f}% -> {worst[2]:.0f}%)"
        return out


def _cv(values: list[float], min_detected: int) -> float:
    """Coefficient of variation, in percent, or `inf` when there is nothing to measure.

    A row detected in only one or two QC injections has no measurable precision, and `inf` is the
    honest answer — but it must not be read as "very imprecise, so merging can only help". That
    is what `_try_merge` guards: an infinite CV is missing evidence, not bad evidence, and the
    two are opposite in what they license.
    """
    detected = [v for v in values if v is not None and v > 0]
    if len(detected) < min_detected or len(values) < 3:
        return float("inf")
    centre = mean(values)
    if centre <= 0:
        return float("inf")
    return 100.0 * stdev(values) / centre


def _finite(value: float) -> bool:
    return value == value and value not in (float("inf"), float("-inf"))


def _ppm(a: float, b: float) -> float:
    return abs(a - b) / b * 1e6 if b else float("inf")


def merge_split_peaks(groups, samples, params: SplitPeakMerge, avg_fwhm: float,
                      log=None) -> SplitPeakReport:
    """Sum rows that are one peak integrated twice. Mutates `groups` in place.

    The winner keeps its identity and takes the summed areas; the loser is marked not-kept with
    a reason naming what absorbed it, so `Unfiltered_Results.csv` still shows it. Compounds move
    across too, so `Features Found` counts the union rather than the winning half.
    """
    say = log or (lambda _: None)
    report = SplitPeakReport()
    pools = [i for i, s in enumerate(samples) if s.role == "qc"]
    # The systematic signature is measured across every injection carrying sample material: a
    # ratio held over the pools alone is not evidence, because the pools are one homogenate.
    columns = [i for i, s in enumerate(samples) if s.role in ("sample", "qc")]
    report.pools = len(pools)
    if not params.enabled:
        return report
    if len(pools) < params.min_pools:
        report.skipped_no_pools = True
        say(str(report))
        return report

    window = avg_fwhm * params.max_peak_widths
    live = sorted((g for g in groups if g.keep and g.quant_ion is not None),
                  key=lambda g: g.retention)

    # Absorbed groups stay out of later comparisons, so three fragments of one peak collapse
    # into one row rather than into an arbitrary pair plus an orphan.
    for i, first in enumerate(live):
        if not first.keep:
            continue
        for second in live[i + 1:]:
            if second.retention - first.retention > window:
                break
            if not second.keep:
                continue
            if not _same_species(first, second):
                continue
            if _ppm(first.quant_ion, second.quant_ion) > params.max_ppm:
                continue
            report.considered += 1
            merged = _try_merge(first, second, pools, params)
            if merged is None and params.systematic:
                merged = _try_systematic(first, second, columns, params)
                if merged is not None:
                    report.systematic += 1
            if merged is not None:
                report.merged += 1
                report.examples.append(merged)

    say(str(report))
    return report


def _try_systematic(winner, loser, columns: list[int], params: SplitPeakMerge):
    """Merge a peak that was cut at the same point in every file.

    Requires the pair to be present together in most injections, to rise and fall together, and to
    hold a near-constant ratio while doing so. The last is what separates one divided peak from two
    molecules that merely co-vary: co-varying molecules reach high correlation readily, but their
    ratio moves.
    """
    if winner.max_area < loser.max_area:
        winner, loser = loser, winner
    pairs = [(winner.areas[i], loser.areas[i]) for i in columns
             if i < len(winner.areas) and i < len(loser.areas)
             and winner.areas[i] > 0 and loser.areas[i] > 0]
    if len(pairs) < params.systematic_min_injections:
        return None

    import statistics
    a = [p[0] for p in pairs]
    b = [p[1] for p in pairs]
    ratios = [x / y for x, y in pairs]
    mean_ratio = statistics.mean(ratios)
    if mean_ratio <= 0:
        return None
    ratio_cv = 100.0 * statistics.stdev(ratios) / mean_ratio
    if ratio_cv > params.systematic_ratio_cv:
        return None

    mean_a, mean_b = statistics.mean(a), statistics.mean(b)
    num = sum((x - mean_a) * (y - mean_b) for x, y in pairs)
    den = (sum((x - mean_a) ** 2 for x in a) * sum((y - mean_b) ** 2 for y in b)) ** 0.5
    if not den or num / den < params.systematic_rho:
        return None

    name, _ = winner.identification()
    _absorb(winner, loser)
    return (f"{name} (systematic)", ratio_cv, num / den)


def _same_species(a, b) -> bool:
    """Same molecule and same ion. Both are required, and for different reasons.

    The name alone would sum `[M+H]+` onto `[M+H-H2O]+`, which are one molecule but two ions with
    two response factors. The adduct alone would sum two different molecules that happen to share
    an ion type.
    """
    name_a, _ = a.identification()
    name_b, _ = b.identification()
    if not name_a or name_a != name_b:
        return False
    return a.adduct() == b.adduct()


def _try_merge(winner, loser, pools: list[int], params: SplitPeakMerge):
    """Merge if summing improves precision in the pools. Returns (name, before, after) or None.

    Which row survives is decided by area, not by score: both rows carry the same identification
    by construction, so the only thing to choose between them is which holds more of the peak.
    """
    if winner.max_area < loser.max_area:
        winner, loser = loser, winner
    a = [winner.areas[i] for i in pools]
    b = [loser.areas[i] for i in pools]
    summed = [x + y for x, y in zip(a, b)]
    cv_a, cv_b, after = (_cv(a, params.min_pools), _cv(b, params.min_pools),
                         _cv(summed, params.min_pools))
    # Every one of the three must be measurable. Without this, a row undetected across the QCs
    # scores an infinite CV, `after * factor <= inf` is true whatever `after` is, and the pair
    # merges on no evidence at all — which on the Kiterie run produced a "best" merge going from
    # an unmeasurable CV to 283%. Missing evidence must block the merge, not wave it through.
    if not all(map(_finite, (cv_a, cv_b, after))):
        return None
    worse = max(cv_a, cv_b)
    if not (after * params.min_cv_improvement <= worse):
        return None

    _absorb(winner, loser)
    name, _ = winner.identification()
    return (name, worse, after)


def _absorb(winner, loser) -> None:
    """The winner takes the summed areas and the loser's compounds; the loser is kept in the
    unfiltered table with a reason naming what absorbed it, so nothing vanishes silently."""
    winner.areas = [x + y for x, y in zip(winner.areas, loser.areas)]
    winner.max_area = max(winner.areas) if winner.areas else winner.max_area
    winner.compounds.extend(loser.compounds)
    loser.keep = False
    name, _ = winner.identification()
    loser.filter_reason = f"Split peak, merged into {name} at {winner.retention:.3f} min"
