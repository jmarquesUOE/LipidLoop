"""Abundance correlation across samples, as a redundancy test.

An adduct, an in-source fragment and an isotope of the same molecule are not independent
measurements of anything: they are the same peak, split. Whatever the molecule does across the
sample set, all of them do together. So if an unidentified feature co-elutes with an identified
one *and* rises and falls with it across every sample, it is almost certainly a second view of
that molecule rather than a compound in its own right.

LipiDex 2 adds this; LipiDex 1 has nothing like it. It is strictly more information than the
mass-difference tests in `adducts.py`, which ask only whether two quant ions *could* be the same
neutral molecule — a question that has a fair number of coincidental yes answers at 20 ppm.
Correlation asks whether they *behave* like it, and coincidence is much harder across 24 samples.

**Only unidentified features are removed.** That is the safety property that makes this
affordable: an identification can never be lost to it, so the worst case is losing an
unannotated row.

Two ways it can mislead, both worth knowing before raising the threshold:

  * **A homogeneous sample set makes everything correlate.** Injection-to-injection loading
    varies, every lipid tracks it, and with replicates of one material the correlation matrix
    goes to one. The threshold has to be judged against how much genuine biological variation
    the set contains, not chosen from a textbook.
  * **Co-regulated lipids are real.** Two species in the same pathway can correlate at 0.95
    honestly. Co-elution is what stops that being a problem here — they also have to be at the
    same retention time — but the risk is not zero.
"""
from __future__ import annotations

from dataclasses import dataclass

DEFAULT_MIN_CORRELATION = 0.95
DEFAULT_MIN_POINTS = 5


def pearson(xs: list[float], ys: list[float]) -> float:
    """Linear correlation. Returns 0 when either series has no spread."""
    n = len(xs)
    if n < 2:
        return 0.0
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    covariance = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    var_x = sum((x - mean_x) ** 2 for x in xs)
    var_y = sum((y - mean_y) ** 2 for y in ys)
    if var_x <= 0 or var_y <= 0:
        return 0.0
    return covariance / (var_x * var_y) ** 0.5


def _ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__)
    ranks = [0.0] * len(values)
    i = 0
    while i < len(order):
        j = i
        while j + 1 < len(order) and values[order[j + 1]] == values[order[i]]:
            j += 1
        shared = (i + j) / 2.0 + 1.0
        for k in range(i, j + 1):
            ranks[order[k]] = shared
        i = j + 1
    return ranks


def spearman(xs: list[float], ys: list[float]) -> float:
    """Rank correlation — robust to the one enormous sample that dominates a Pearson."""
    if len(xs) < 2:
        return 0.0
    return pearson(_ranks(xs), _ranks(ys))


CORRELATIONS = {"pearson": pearson, "spearman": spearman}


@dataclass(slots=True)
class CorrelationFilter:
    """How correlation is used.

    `guard_mass_relationships` is the well-founded use and is on: before removing a group as an
    adduct, in-source fragment, dimer or isotope of a *specific* partner, require the two to
    track each other. Measured on this data the mass tests are already enriched for it —
    in-source removals have a median r of 0.96 and adduct-of-identified 0.93 against a 0.85
    baseline — so the guard mostly confirms, and refuses the minority that do not behave.

    `enabled` is the standalone LipiDex 2-style filter: remove any unidentified feature that
    co-elutes with and tracks *any* identified one, with no mass relationship required. It is
    off, because measured here it is far too blunt — see docs/CORRELATION.md.
    """

    enabled: bool = False
    guard_mass_relationships: bool = True
    guard_min_correlation: float = 0.5
    min_correlation: float = DEFAULT_MIN_CORRELATION
    correlation_type: str = "pearson"
    min_points: int = DEFAULT_MIN_POINTS
    fwhm_divisor: float = 0.5      # the same co-elution window the adduct sweep uses


def correlate_areas(a: list[float], b: list[float], columns: list[int],
                    method="pearson", min_points: int = DEFAULT_MIN_POINTS) -> float | None:
    """Correlation over the samples where both were detected. None if too few."""
    xs, ys = [], []
    for i in columns:
        if i < len(a) and i < len(b) and a[i] > 0 and b[i] > 0:
            xs.append(a[i])
            ys.append(b[i])
    if len(xs) < min_points:
        return None
    return CORRELATIONS[method](xs, ys)


def filter_correlated_unidentified(groups, samples, config: CorrelationFilter) -> int:
    """Remove unidentified groups that co-elute with, and track, an identified one.

    `groups` must be retention-sorted, as the peak finder keeps them. Returns the number
    removed.
    """
    if not config.enabled:
        return 0

    columns = [i for i, s in enumerate(samples) if s.role != "blank"]
    if len(columns) < config.min_points:
        return 0

    identified = [g for g in groups if g.keep and g.final_lipid_id is not None]
    if not identified:
        return 0

    # Bucket the identified groups by retention so each candidate is a couple of lookups.
    bucket_size = 0.1
    index: dict[int, list] = {}
    for group in identified:
        index.setdefault(int(group.retention / bucket_size), []).append(group)

    removed = 0
    for group in groups:
        if not group.keep or group.final_lipid_id is not None:
            continue
        window = group.avg_fwhm / config.fwhm_divisor if config.fwhm_divisor else 0.0
        if window <= 0:
            continue
        lo = int((group.retention - window) / bucket_size)
        hi = int((group.retention + window) / bucket_size)
        for b in range(lo, hi + 1):
            for other in index.get(b, ()):
                if abs(other.retention - group.retention) >= window:
                    continue
                r = correlate_areas(group.areas, other.areas, columns,
                                    config.correlation_type, config.min_points)
                if r is not None and r >= config.min_correlation:
                    group.keep = False
                    group.filter_reason = "Correlated with identified peak"
                    removed += 1
                    break
            if not group.keep:
                break
    return removed
