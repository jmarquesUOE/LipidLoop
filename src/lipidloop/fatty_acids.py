"""Name free fatty acids from accurate mass and retention, because they will not fragment.

A free fatty acid is the one lipid class that spectral matching cannot reach. The carboxylate
anion is stable: across 16 files, `FA 16:0` produced 711 MS2 spectra and `FA 18:0` 561, and not one
of them contains a reproducible product ion. Only polyunsaturated acids lose CO2, so
`FreeFattyAcids_Negative.msp` covers DB >= 3 and stops there. Everything from `FA 14:0` to
`FA 26:1` is invisible to a dot product no matter how abundant it is — and in the skin organoid
data `FA 18:1` is a 1.9e8 peak found in all 63 files.

**What identifies them instead is the homologous series itself.** One unidentified feature at
m/z 281.2486 could be anything. Thirty-three features sitting on thirty-three exact fatty-acid
masses, whose retention times jointly fit

    RT = intercept + per_carbon x C + per_double_bond x DB

with R2 0.94, a positive carbon term and a negative double-bond term, could not be a coincidence:
that is a chromatographic law being obeyed, and obeying it is the evidence. On the skin negative
data the fit is `RT = 5.08 + 0.327xC - 0.446xDB`, residual sd 0.344 min over C14-C26.

So this is deliberately **not** a per-feature test. Each candidate on its own is only a mass, and
a mass is not an identification. The set is tested as a set, and only if the set behaves like a
homologous series is any member named.

## The gates, and why each is there

`MIN_CANDIDATES` — a handful of masses will always land somewhere; a series needs members.
`MIN_R2` and `MAX_RESIDUAL_SD` — the fit has to be good, or there is no law being obeyed.
`per_carbon > 0` and `per_double_bond < 0` — the signs are the chemistry. A reversed-phase column
retains longer chains and releases unsaturated ones. A "model" with the wrong signs has fitted
noise, and no number of points rescues it.
`MAX_SIGMA` — a member more than this far off the fitted surface is not named, even though its
mass matches.

If the fit fails any gate, **nothing is named at all**. A bad fit is an absence of evidence.

Names are written with `Identification Source` = `RT model`, never as MS2 evidence, so a
downstream reader can always tell which rows were named without a spectrum.
"""
from __future__ import annotations

from dataclasses import dataclass

PROTON = 1.00727646
CARBON = 12.0
HYDROGEN = 1.0078250319
OXYGEN = 15.9949146221

MIN_CARBONS = 12
MAX_CARBONS = 30
MAX_DOUBLE_BONDS = 6
MZ_TOLERANCE = 0.01

MIN_CANDIDATES = 8         # fewer than this is not a series
MIN_R2 = 0.85
MAX_RESIDUAL_SD = 0.75     # minutes
MAX_SIGMA = 3.0            # a member this far off the surface is not named
MIN_FILES = 3              # a candidate seen in fewer files than this does not join the fit
# Retention against carbon number is not straight. Over C14-C30 the step between consecutive
# saturated acids grows monotonically from 0.21 to 0.81 min, so a line leaves systematic residuals
# and a wider naming window than the data deserves: sd 0.530 min linear against 0.393 with a C^2
# term, which tightens the 3-sigma window from 1.63 to 1.18 min. The quadratic is only fitted when
# there are enough points to determine it, and only kept when it actually helps.
MIN_POINTS_FOR_CURVE = 12
CURVE_IMPROVEMENT = 0.9    # the quadratic must cut the residual sd by at least this factor
CURVE_FLOOR = 0.02         # minutes: below this a line already fits, and a curve fits only noise
# Absolute floor on the trimming threshold. Without it, a fit that is very good collapses: the
# residual spread goes to floating-point noise, `3 x spread` becomes ~1e-15, and almost every
# point is trimmed as an outlier from a model it fits exactly. The surviving handful then spans a
# fraction of the carbon range and the sign gates are checked over the wrong interval.
MIN_TRIM_WINDOW = 0.05     # minutes


def neutral_mass(carbons: int, double_bonds: int) -> float:
    return (carbons * CARBON + (2 * carbons - 2 * double_bonds) * HYDROGEN + 2 * OXYGEN)


def deprotonated(carbons: int, double_bonds: int) -> float:
    return neutral_mass(carbons, double_bonds) - PROTON


def max_double_bonds(carbons: int) -> int:
    """How many double bonds a chain of this length can carry.

    Natural polyenes are methylene-interrupted, so each double bond after the first needs three
    more carbons: 1 + (C-3)//3. That gives 20:6 and 22:7 as the ceilings — comfortably above
    arachidonic (20:4), EPA (20:5) and DHA (22:6) — while refusing `FA 14:6`, which would
    otherwise be enumerated and whose mass would then claim a real feature.

    Odd-chain is a separate case, not a smaller version of the even-chain rule. Even-chain
    polyenes come from acetyl-CoA elongation plus desaturation; odd-chain fatty acids come from
    propionyl-CoA priming instead (gut-microbial propionate, or catabolism of branched/odd
    precursors) — a pathway that does not build polyenes. Mammalian odd-chain species are
    essentially always saturated, with monounsaturated forms (margaroleic acid, FA 17:1) at trace
    level and nothing beyond that reported. Without this, `candidate_masses()` enumerates FA 13:3,
    15:4, 19:5 and similar — caught live on the NIST SRM 1950 method-development batch, where FA
    15:4 was a genuine, reproducible, blank-clear peak that is chemically not a real compound. A
    real peak on an impossible mass still needs a name from somewhere, and an over-wide candidate
    grid is happy to offer it one.
    """
    if carbons % 2:
        return min(1, MAX_DOUBLE_BONDS)
    return max(0, min(MAX_DOUBLE_BONDS, 1 + (carbons - 3) // 3))


def candidate_masses(min_carbons: int = MIN_CARBONS, max_carbons: int = MAX_CARBONS,
                     ceiling: int = MAX_DOUBLE_BONDS) -> dict:
    """`{(C, DB): [M-H]- m/z}` for every plausible free fatty acid."""
    out = {}
    for carbons in range(min_carbons, max_carbons + 1):
        for double_bonds in range(0, min(ceiling, max_double_bonds(carbons)) + 1):
            out[(carbons, double_bonds)] = deprotonated(carbons, double_bonds)
    return out


@dataclass(slots=True)
class SeriesModel:
    intercept: float
    per_carbon: float
    per_double_bond: float
    r2: float
    residual_sd: float
    n_used: int
    per_carbon_squared: float = 0.0
    carbon_range: tuple = (0, 0)

    def slope_at(self, carbons: float) -> float:
        """d(retention)/d(carbon) at this chain length."""
        return self.per_carbon + 2.0 * self.per_carbon_squared * carbons

    @property
    def usable(self) -> bool:
        """The gates. The sign checks are the chemistry, and they are checked on the slope.

        With a C^2 term the bare linear coefficient is no longer the per-carbon slope — on the
        skin data the quadratic fit is `-0.441/C +0.0192/C^2`, whose slope at C20 is +0.327,
        entirely correct, while the linear coefficient alone is negative. Testing that coefficient
        rejected a good model and named nothing. What has to be positive is the derivative, across
        the whole range the model was fitted over, and since it is linear in C the two ends decide
        it.
        """
        low, high = self.carbon_range if self.carbon_range != (0, 0) else (1, 1)
        rising = self.slope_at(low) > 0 and self.slope_at(high) > 0
        return (self.n_used >= MIN_CANDIDATES and self.r2 >= MIN_R2
                and self.residual_sd <= MAX_RESIDUAL_SD
                and rising and self.per_double_bond < 0)

    def predict(self, carbons: int, double_bonds: int) -> float:
        return (self.intercept + self.per_carbon * carbons
                + self.per_carbon_squared * carbons * carbons
                + self.per_double_bond * double_bonds)

    def __str__(self) -> str:
        low, high = self.carbon_range
        if self.per_carbon_squared:
            slope = f"{self.slope_at(low):+.3f}..{self.slope_at(high):+.3f}/C (curved)"
        else:
            slope = f"{self.per_carbon:+.3f}/C"
        return (f"FA n={self.n_used} R2={self.r2:.3f} sd={self.residual_sd:.3f} min  "
                f"{slope} {self.per_double_bond:+.3f}/DB"
                f"{'' if self.usable else '   [not usable]'}")


def _fit_once(points, quadratic: bool, trim_sigma: float, iterations: int,
              protected: int = 0):
    """Trimmed least squares. Returns `(SeriesModel, keep mask)` or None."""
    import numpy as np

    carbons = np.array([p[0] for p in points], float)
    doubles = np.array([p[1] for p in points], float)
    retention = np.array([p[2] for p in points], float)
    terms = 4 if quadratic else 3
    keep = np.ones(len(points), bool)

    def design(mask):
        columns = [np.ones(mask.sum()), carbons[mask], doubles[mask]]
        if quadratic:
            columns.append(carbons[mask] ** 2)
        return np.column_stack(columns)

    def predict(coefficients):
        out = coefficients[0] + coefficients[1] * carbons + coefficients[2] * doubles
        return out + coefficients[3] * carbons ** 2 if quadratic else out

    coefficients = None
    for _ in range(iterations + 1):
        matrix = design(keep)
        if matrix.shape[0] < max(MIN_CANDIDATES, terms + 2) \
                or np.linalg.matrix_rank(matrix) < terms:
            return None
        coefficients, *_ = np.linalg.lstsq(matrix, retention[keep], rcond=None)
        residuals = retention - predict(coefficients)
        spread = float(np.std(residuals[keep], ddof=1)) if keep.sum() > 3 else 0.0
        if spread <= 0:
            break
        fresh = keep & (np.abs(residuals) <= max(trim_sigma * spread, MIN_TRIM_WINDOW))
        if protected:
            fresh[:protected] = True      # standards are never trimmed
        if fresh.sum() < max(MIN_CANDIDATES, terms + 2) or fresh.sum() == keep.sum():
            break
        keep = fresh

    residuals = (retention - predict(coefficients))[keep]
    observed = retention[keep]
    total = float(((observed - observed.mean()) ** 2).sum())
    r2 = 1.0 - float((residuals ** 2).sum()) / total if total > 0 else 0.0
    sd = float(np.std(residuals, ddof=terms)) if keep.sum() > terms else 0.0
    used = carbons[keep]
    return SeriesModel(float(coefficients[0]), float(coefficients[1]), float(coefficients[2]),
                       r2, sd, int(keep.sum()),
                       float(coefficients[3]) if quadratic else 0.0,
                       (float(used.min()), float(used.max())))


def fit_series(points, trim_sigma: float = 3.0, iterations: int = 3,
               anchors=()) -> "SeriesModel | None":
    """Fit `RT ~ C + DB`, and `+ C^2` as well when the data can determine it.

    `anchors` are measured standards. They join the fit and are **never trimmed** — a point whose
    identity is known cannot be an outlier from a surface meant to describe it; if it disagrees
    with the surface, the surface is wrong. Observed candidates are trimmed as usual.

    The quadratic is kept only if it cuts the residual sd by a real margin. Curvature over a
    C14-C30 range is real on this method, but it is a property of the gradient rather than of the
    chemistry — the whole series elutes inside a rising IPA ramp — so it does not transfer to
    another method, and a curve fitted to a handful of points extrapolates badly. The residual sd
    is what sets the window inside which a feature gets named, so a spuriously small one is the
    dangerous failure, not a slightly wide one.
    """
    points = list(anchors) + list(points)
    protected = len(list(anchors))
    if len(points) < MIN_CANDIDATES:
        return None
    linear = _fit_once(points, False, trim_sigma, iterations, protected)
    if linear is None or len(points) < MIN_POINTS_FOR_CURVE:
        return linear
    curved = _fit_once(points, True, trim_sigma, iterations, protected)
    if linear.residual_sd <= CURVE_FLOOR:
        return linear      # a line already predicts retention to better than a second
    if curved is None or curved.residual_sd > CURVE_IMPROVEMENT * linear.residual_sd:
        return linear
    return curved


def read_standards(path) -> dict:
    """`{(C, DB): retention}` from a CSV of authentic standards — `name,retention` per row.

    `FA 18:1,10.36`. Blank path or missing file gives an empty dict, which is the ordinary case:
    the surface is then fitted from the observed series alone.
    """
    import csv
    import re
    from pathlib import Path

    if not path or not Path(path).exists():
        return {}
    out = {}
    with Path(path).open(newline="") as fh:
        for row in csv.reader(fh):
            if len(row) < 2:
                continue
            match = re.search(r"(\d+):(\d+)", row[0])
            try:
                retention = float(row[1])
            except ValueError:
                continue          # header line, or a comment
            if match:
                out[(int(match.group(1)), int(match.group(2)))] = retention
    return out


def annotate(groups, polarity: str = "-", mz_tolerance: float = MZ_TOLERANCE,
             min_files: int = MIN_FILES, max_sigma: float = MAX_SIGMA,
             standards: "dict | None" = None, log=None):
    """Name unidentified compound groups that form a fatty-acid homologous series.

    Only groups that are still kept, still unidentified, and of the right polarity are eligible —
    nothing this does can overwrite an identification made from a spectrum.

    `standards` maps `(C, DB)` to a measured retention time from an authentic standard run on the
    same method. When given, those points anchor the surface and are never trimmed out of it,
    while the observed series still extends it into the compositions no standard covers — the
    hybrid model. Ten standards spread across the range pin the shape far better than fifty
    observations that all have to be taken on trust.
    """
    say = log or (lambda _: None)
    masses = candidate_masses()
    standards = standards or {}

    eligible = [g for g in groups
                if g.keep and not g.identification()[0] and g.quant_polarity == polarity]
    matched = []
    for group in eligible:
        for (carbons, double_bonds), mz in masses.items():
            if abs(group.quant_ion - mz) <= mz_tolerance:
                matched.append((carbons, double_bonds, group))
                break
    if not matched:
        return 0

    # One feature per composition for the fit: the most abundant, and only if it is seen in
    # enough files. Fitting on every match would let a scattering of one-off features set the
    # surface that is meant to judge them.
    best: dict = {}
    for carbons, double_bonds, group in matched:
        if len(group.compounds) < min_files:
            continue
        key = (carbons, double_bonds)
        if key not in best or group.max_area > best[key].max_area:
            best[key] = group
    observed = [(c, db, g.retention) for (c, db), g in best.items()]
    anchors = [(c, db, rt) for (c, db), rt in sorted(standards.items())]
    if anchors:
        say(f"  free fatty acids: {len(anchors)} authentic standards anchor the surface")
    model = fit_series(observed, anchors=anchors)
    if model is None or not model.usable:
        say(f"  free fatty acids: {len(best)} candidate masses, "
            f"{'no model could be fitted' if model is None else str(model)} — nothing named")
        return 0
    say(f"  free fatty acids: {model}")

    limit = max_sigma * model.residual_sd
    named = 0
    for carbons, double_bonds, group in matched:
        error = group.retention - model.predict(carbons, double_bonds)
        if abs(error) > limit:
            continue
        group.rtls_identification = f"FA {carbons}:{double_bonds}"
        named += 1
    say(f"  free fatty acids: {named} features named within {limit:.2f} min of the surface")
    return named
