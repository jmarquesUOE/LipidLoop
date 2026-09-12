"""Retention-time modelling within a lipid class — the criterion LipiDex 1 never used.

On a reversed-phase column a lipid's retention is largely set by two numbers: how many acyl
carbons it carries, and how many double bonds. More carbons is more hydrophobic and elutes
later; more double bonds kinks the chain and elutes earlier. Within one class that relationship
is close to linear, so

    retention ~ intercept + per_carbon * carbons + per_double_bond * double_bonds

is a usable model, and an identification that sits far off its class's own surface is probably
wrong. LipiDex 1 does not use this at all — its only retention test on an identification asks
whether the lipid lands in the window where its class elutes, never whether it lands where *that
member* of the class should.

**Provenance.** This is the *effective carbon number* relationship, long established in
reversed-phase lipid chromatography, and LipiDex 2 applies it as "ECN versus retention time"
inside its Degreaser quality-control module (Hutchins et al., *Cell Systems* 2018 for LipiDex 1;
LipiDex 2, 2024). ⚠ An earlier version of this docstring credited it to "RTLS", which is wrong:
RTLS is LipiDex 2's *real-time library search*, an acquisition method that triggers class-targeted
MSn, and has nothing to do with retention modelling.

**One deliberate departure from textbook ECN.** ECN is usually a single number, ECN = C - 2*DB,
which fixes the double-bond coefficient at exactly twice the carbon one. Here the two are fitted
independently, and on real data the ratio is not 2: across the eight best-fitting classes of one
study it ran 1.24 to 1.86, mean 1.50, and differed by class — 1.24 for PE against 1.86 for
Cer[NS]. Forcing the textbook 2.0 would bias every class, and by different amounts.

Fitted on this data the coefficients come out chemically sensible everywhere the fit is good:
+0.14 to +0.42 min per carbon, -0.18 to -0.76 min per double bond.

**The fit is diagnostic in its own right.** A class whose members are correctly identified fits
tightly — Alkenyl-TG R² 0.98, AC 0.94, Cer[NS] 0.92. A class that does not fit is telling you
something: PS comes out at R² 0.19 with a residual spread of three minutes and coefficients of
the wrong sign, which is not a failure of the model but a statement that those identifications
are not a homologous series.

Because the data being fitted contains the outliers being looked for, the fit is iterative:
fit, drop the worst residuals, refit. Without that a handful of bad rows drag the surface onto
themselves and nothing looks like an outlier.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

from .nomenclature import canonical_class

import numpy as np

from .peaks import sum_composition

# `PC 34:1`, `SM d36:1`, `Plasmenyl-PE P-38:4`, `PC[OH] OH-32:5`
_SUM_NAME = re.compile(r"^(\S+)\s+([A-Za-z]*-?)(\d+):(\d+)$")

MIN_POINTS = 6           # below this there is no surface to speak of
TRIM_SIGMA = 3.0         # residuals beyond this are dropped and the model refitted
MAX_ITERATIONS = 3
MIN_R2 = 0.7             # a model worse than this is not used to judge anything
MAX_RESIDUAL_SD = 1.0    # minutes


def parse_sum_name(name: str) -> tuple[str, str, int, int] | None:
    """`SM d36:1` -> ("SM", "d", 36, 1). None if the name is not a sum composition.

    A molecular name is collapsed first: `TG 16:0_18:1_18:2` becomes `TG 52:3` and is judged the
    same as any other TG 52:3. Without this the model silently declines to judge every row reported
    at molecular resolution, which on the skin set is 78 identifications: every TG and DG whose
    fatty-acid purity cleared 75%, and the gangliosides. They were not being protected; they were
    being skipped.

    ⚠ **This is a limitation accepted, not a distinction without a difference**, and an earlier
    version of this docstring claimed the latter — that same-sum rows "elute within seconds of each
    other". They do not. Measured over 131 duplicated names in one positive-mode study, the gap
    between rows sharing a name has a median of 0.41 min; only 2% are within 0.1 min and 36% exceed
    half a minute. `PC 18:0_20:3` and `PC 18:1_20:2` are chromatographically distinct and both
    return ("PC", "", 38, 3).

    So the model has no term that could separate them — true — but the consequence is that it gives
    every isomer of a composition ONE prediction, and all but one of them is displaced from it by
    construction. Judge a row against this model only when its name appears once; see the note above
    `peakfinder.UNRELIABLE_CV` for the verdict this cost.
    """
    text = (name or "").strip()
    match = _SUM_NAME.match(text)
    if not match:
        match = _SUM_NAME.match(sum_composition(text))
    if not match:
        return None
    # ⚠ Canonicalise the class HERE, because every grouping decision downstream reads it.
    #
    # The model fits ONE SURFACE PER CLASS at MIN_POINTS = 6. Our own two libraries write the same
    # molecules as `GlcCer[NS]` and `HexCer[NS]`, so a study searching both split one class into
    # two half-sized ones and could leave BOTH under the minimum — losing the model for
    # hexosylceramides entirely, and with it the ability to extend those identifications or to
    # catch a subclass substitution in them.
    return (canonical_class(match.group(1)), match.group(2),
            int(match.group(3)), int(match.group(4)))


@dataclass(slots=True)
class RetentionModel:
    lipid_class: str
    prefix: str
    intercept: float
    per_carbon: float
    per_double_bond: float
    r2: float
    residual_sd: float
    n_used: int
    n_total: int

    @property
    def key(self) -> tuple[str, str]:
        return (self.lipid_class, self.prefix)

    @property
    def usable(self) -> bool:
        """Whether this model is good enough to judge an identification by.

        A poorly fitting class is not evidence that its members are wrong — it is an absence of
        evidence either way, and must not be used to remove anything.
        """
        return (self.n_used >= MIN_POINTS and self.r2 >= MIN_R2
                and self.residual_sd <= MAX_RESIDUAL_SD
                and self.per_carbon > 0 and self.per_double_bond < 0)

    def predict(self, carbons: int, double_bonds: int) -> float:
        return self.intercept + self.per_carbon * carbons + self.per_double_bond * double_bonds

    def __str__(self) -> str:
        return (f"{self.lipid_class}{' ' + self.prefix if self.prefix else '':<4} "
                f"n={self.n_used}/{self.n_total} R2={self.r2:.3f} "
                f"sd={self.residual_sd:.3f} min  "
                f"{self.per_carbon:+.3f}/C {self.per_double_bond:+.3f}/DB"
                f"{'' if self.usable else '   [not usable]'}")


def _consensus(carbons, doubles, retention):
    """Which points to start the fit from, found by consensus rather than by using all of them.

    ⚠ Seeding the trim loop with every point lets the contamination set the threshold meant to
    catch it. On this data the PS class carried six end-of-gradient rows against sixteen real ones;
    the first least-squares line was dragged towards them, the residual spread came out at 6.5 min,
    nothing was ever trimmed, the class failed to fit — and because an unfittable class is exempt
    from the retention filter, the junk that broke the fit escaped it. **Contamination bought
    immunity from the test designed to catch it.**

    Robustifying the trim scale alone was not enough: at one reliability threshold the class fitted
    at R2 0.942 and at a slightly stricter one it collapsed to R2 -1.070, which is a knife-edge
    rather than a fix. RANSAC finds the largest self-consistent subset instead of trusting a seed,
    and recovers PS at R2 0.880 with none of the wash retained.

    Nor could this be caught without a model: retention runs about +0.25 min per carbon and
    -0.36 per double bond, so extra carbons are offset by extra double bonds and species far apart
    in composition genuinely co-elute. "Too many molecules too close together" flags real
    chromatography.

    Falls back to using every point when scikit-learn is absent or the fit fails, which is the
    previous behaviour.
    """
    everything = np.ones(len(retention), dtype=bool)
    if len(retention) < MIN_POINTS + 2:
        return everything
    try:
        from sklearn.linear_model import LinearRegression, RANSACRegressor
        design = np.column_stack([carbons, doubles])
        model = RANSACRegressor(LinearRegression(), random_state=0).fit(design, retention)
        inliers = np.asarray(model.inlier_mask_, dtype=bool)
        return inliers if inliers.sum() >= MIN_POINTS else everything
    except Exception:                      # noqa: BLE001 - a failed consensus is not fatal
        return everything


def fit_model(lipid_class: str, prefix: str, points: list[tuple[int, int, float]],
              trim_sigma: float = TRIM_SIGMA,
              max_iterations: int = MAX_ITERATIONS) -> RetentionModel | None:
    """Least squares with iterative trimming. `points` are (carbons, double_bonds, retention)."""
    if len(points) < MIN_POINTS:
        return None

    carbons = np.array([p[0] for p in points], dtype=float)
    doubles = np.array([p[1] for p in points], dtype=float)
    retention = np.array([p[2] for p in points], dtype=float)
    keep = _consensus(carbons, doubles, retention)

    coefficients = None
    for _ in range(max_iterations + 1):
        design = np.column_stack([np.ones(keep.sum()), carbons[keep], doubles[keep]])
        if design.shape[0] < MIN_POINTS or np.linalg.matrix_rank(design) < 3:
            return None
        coefficients, *_ = np.linalg.lstsq(design, retention[keep], rcond=None)

        predicted = coefficients[0] + coefficients[1] * carbons + coefficients[2] * doubles
        residuals = retention - predicted
        # ⚠ A ROBUST scale, not the standard deviation. With std the outliers set the very
        # threshold meant to catch them: on this data PS carried five end-of-gradient rows at
        # 24.5 min against a real population at 7-13, the residual spread came out at 6.5 min,
        # `trim_sigma * spread` was therefore enormous, and nothing was ever trimmed. The class
        # then failed to fit (R2 0.020) and — because an unfittable class is exempted — the junk
        # that broke the fit escaped the filter entirely. Contamination bought immunity from the
        # test designed to catch it.
        #
        # With MAD the same class trims 12 of 29 points and fits at R2 0.942, sd 0.218 min.
        centre = np.median(residuals[keep])
        spread = float(np.median(np.abs(residuals[keep] - centre)) * 1.4826) if keep.sum() > 3 \
            else 0.0
        if spread <= 0:
            break
        fresh = keep & (np.abs(residuals) <= trim_sigma * spread)
        if fresh.sum() < MIN_POINTS or fresh.sum() == keep.sum():
            break
        keep = fresh

    predicted = coefficients[0] + coefficients[1] * carbons + coefficients[2] * doubles
    residuals = retention[keep] - predicted[keep]
    total = np.sum((retention[keep] - retention[keep].mean()) ** 2)
    r2 = 1.0 - float(np.sum(residuals ** 2) / total) if total > 0 else 0.0
    residual_sd = float(np.std(residuals, ddof=1)) if keep.sum() > 3 else 0.0

    return RetentionModel(lipid_class=lipid_class, prefix=prefix,
                          intercept=float(coefficients[0]), per_carbon=float(coefficients[1]),
                          per_double_bond=float(coefficients[2]), r2=r2,
                          residual_sd=residual_sd, n_used=int(keep.sum()), n_total=len(points))


def fit_models(observations: list[tuple[str, float]]) -> dict[tuple[str, str], RetentionModel]:
    """Fit one model per (class, chain prefix). `observations` are (sum name, retention)."""
    grouped: dict[tuple[str, str], list[tuple[int, int, float]]] = {}
    for name, retention in observations:
        parsed = parse_sum_name(name)
        if parsed is None:
            continue
        lipid_class, prefix, carbons, doubles = parsed
        grouped.setdefault((lipid_class, prefix), []).append((carbons, doubles, retention))

    models: dict[tuple[str, str], RetentionModel] = {}
    for (lipid_class, prefix), points in grouped.items():
        model = fit_model(lipid_class, prefix, points)
        if model is not None:
            models[(lipid_class, prefix)] = model
    return models


def retention_error(models, name: str, retention: float) -> float | None:
    """How far an identification sits from where its class says it should. None if unjudgeable."""
    parsed = parse_sum_name(name)
    if parsed is None:
        return None
    model = models.get((parsed[0], parsed[1]))
    if model is None or not model.usable:
        return None
    return retention - model.predict(parsed[2], parsed[3])


def extend_identifications(groups, models, candidates, mz_tol_ppm: float = 20.0,
                           sigma: float = 2.0, max_rt_error: float = 0.5) -> int:
    """Name unidentified features from the retention model — "extending" the library's reach.

    A model fitted on the lipids that *were* identified predicts where any other member of the
    class should elute. An unidentified feature whose mass matches such a member, at the
    retention time the model predicts, is a plausible assignment even though no MS2 was ever
    taken of it.

    **This is weaker evidence than an identification and is kept separate from one.** It rests
    on an accurate mass and a retention time, with no fragmentation behind it, so it cannot
    distinguish isomers and cannot tell a lipid from anything else of the same formula eluting
    at the same moment. Assignments land in `rtls_identification`, are reported with their
    source, and are never merged into the MS2-derived identifications.

    Ambiguity is declined rather than guessed: if two candidates both fit, neither is assigned.

    `candidates` are (sum name, m/z, polarity) or (sum name, m/z, polarity, adduct). The adduct
    is worth carrying even though the *name* deliberately omits it: the match is to one m/z, and
    different adducts of one composition are far apart — `[M+NH4]+` and `[M+Na]+` differ by about
    0.9 Da against a 10 ppm window that is 0.008 Da at m/z 800 — so a surviving hit names exactly
    one ion. Dropping it left every retention-model row with a blank `Adduct`, which also exempted
    those rows from adduct-pair removal: the check that catches `Na - H` against `+2C +3DB` at
    2.4 mDa cannot run on a row whose ion is unrecorded.

    Returns the number assigned.
    """
    usable = {k: m for k, m in models.items() if m.usable}
    if not usable:
        return 0

    # Index candidates by integer m/z, keeping only those whose class has a usable model.
    index: dict[int, list[tuple[str, float, str, RetentionModel, int, int]]] = {}
    for candidate in candidates:
        name, mz, polarity = candidate[0], candidate[1], candidate[2]
        adduct = candidate[3] if len(candidate) > 3 else ""
        parsed = parse_sum_name(name)
        if parsed is None:
            continue
        model = usable.get((parsed[0], parsed[1]))
        if model is None:
            continue
        index.setdefault(int(mz), []).append(
            (name, mz, polarity, model, parsed[2], parsed[3], adduct))

    assigned = 0
    for group in groups:
        if not group.keep or group.final_lipid_id is not None or group.quant_ion is None:
            continue
        hits = []
        for key in (int(group.quant_ion) - 1, int(group.quant_ion), int(group.quant_ion) + 1):
            for name, mz, polarity, model, carbons, doubles, adduct in index.get(key, ()):
                if polarity and group.quant_polarity and polarity != group.quant_polarity:
                    continue
                if abs(mz - group.quant_ion) / group.quant_ion * 1e6 > mz_tol_ppm:
                    continue
                error = abs(group.retention - model.predict(carbons, doubles))
                if error <= min(max_rt_error, sigma * model.residual_sd) or \
                        error <= sigma * model.residual_sd <= max_rt_error:
                    hits.append((error, name, adduct))
        if len(hits) != 1:
            continue          # nothing, or ambiguous — decline rather than guess
        group.rtls_identification = hits[0][1]
        group.rtls_adduct = hits[0][2]
        assigned += 1
    return assigned
