"""Quality control for a finished run — the analysis half.

Everything here answers a question asked *before* any biology, and the answers go into one
report (`scripts/qc_report.py`) that is sent to the person whose samples these are. That reader
is not a mass spectrometrist, so every number has to be one they can act on.

Five arms, each generalised from the per-study scripts written for the ALS brain batch:

    precision      how repeatable is the measurement, and does any injection not belong
    run order      is acquisition drift confounded with the biology
    annotation     is one peak being reported as two lipids
    blanks         what is in the background, and is it contamination or carryover
    coverage       what was actually measured, and on what evidence

## Two decisions worth stating

**A pooled QC is screened on its profile, not its total signal.** An injection can be small and
still be a faithful copy of the pool; it can also be full-sized and wrong. Correlating each pool
against the median sample profile after median normalisation is the test that separates those,
and it is what caught a pool at rho 0.785 on the July batch while total signal looked survivable.

**Precision and contrast are different questions.** The CV says whether the instrument repeats.
The ratio of sample spread to QC spread says whether the biology is larger than the noise. A run
can pass the first and fail the second, and then the data is repeatable but cannot answer the
question it was collected for. Both are reported, and the report says which is which.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .blanks import infer_role

MIN_DETECTED_FRACTION = 0.5     # a feature must be seen in this share of the injections used
DRIFT_RHO = 0.5                 # |Spearman| against injection number that counts as drifting
RUN_IN = 5                      # leading injections examined separately for a start-of-run effect


# ── loading ────────────────────────────────────────────────────────────────────────────────

@dataclass(slots=True)
class Run:
    """A finished run: the feature table, the injections, and what each injection is."""

    rows: list
    columns: list
    areas: np.ndarray
    roles: dict
    order: dict = field(default_factory=dict)
    # Where the order came from, and whether the two sources agreed. A drift analysis rests
    # entirely on the order being right, so how it was obtained belongs in the report rather than
    # being an invisible implementation choice.
    order_source: str = ""
    order_disagreement: str = ""
    metadata: dict = field(default_factory=dict)   # column name -> {factor: value}

    def indices(self, *roles) -> list:
        return [i for i, c in enumerate(self.columns) if self.roles[c] in roles]

    @property
    def samples(self) -> list:
        return self.indices("sample")

    @property
    def pools(self) -> list:
        return self.indices("qc")

    @property
    def blank_columns(self) -> list:
        return self.indices("blank")


def usable_factors(metadata: dict) -> list[str]:
    """Metadata columns that could actually be a study group.

    Two are excluded on principle rather than by name, because the next study will use different
    names. A column with **one level** expresses no contrast — `sample_type` was "Snap frozen
    duodenum tissue" for all ten samples and appeared as a study group. A column where **every
    level holds one sample** is an identifier, not a factor.

    What survives is still not guaranteed to be the comparison the submitter cares about; a
    replicate index survives both rules and is a pairing, not a group. That distinction needs
    declaring, not guessing, which is what the study configuration will carry.
    """
    from collections import Counter

    keys = sorted({k for row in metadata.values() for k in row})
    out = []
    for key in keys:
        counts = Counter(row.get(key) for row in metadata.values() if row.get(key))
        if len(counts) < 2:
            continue                       # no contrast
        if max(counts.values()) < 2:
            continue                       # an identifier
        out.append(key)
    return out


_STAMP = re.compile(r'startTimeStamp="([^"]+)"')


def _acquisition_time(path) -> str:
    """The instrument's own record of when this run started, from the mzML header.

    ⚠ **Read until the spectra start, not a fixed 8 kB.** Writers place `startTimeStamp` at very
    different depths — 4.1 kB (Thermo), 6.6 kB (Sciex), 9.5 kB and 14.5 kB (Waters, Agilent) — so
    an 8 kB window found it for three studies of five and silently returned "" for the rest. The
    cost is not a blank field: with no order, the drift, run-order and confounding sections of the
    quality report cannot run, and they report "not testable" rather than saying the input was
    missing.
    """
    try:
        with Path(path).open(errors="replace") as fh:
            head = ""
            while len(head) < 1_000_000:
                chunk = fh.read(65536)
                if not chunk:
                    break
                head += chunk
                match = _STAMP.search(head)
                if match:
                    return match.group(1)
                # The stamp lives in <run>, which precedes the spectra. Once those start it is
                # not coming, and reading a 780 MB file to prove that is not free.
                if "<spectrumList" in head or "<spectrum " in head:
                    break
    except Exception:
        return ""
    return ""


def load_run(results, sequence=None, metadata=None, mzml_files=None) -> Run:
    with Path(results).open(newline="") as fh:
        reader = csv.DictReader(fh)
        # ⚠ Take the columns from the HEADER, not from the first row. A delivered table can be
        # legitimately empty — a two-injection run cannot satisfy a 70% presence filter, and a
        # study whose features are all sparse will filter down to nothing — and `rows[0]` then
        # raises IndexError deep inside the report builder, which reads as a crash rather than as
        # "nothing survived the filters". Two public datasets hit this in one batch.
        columns = list(reader.fieldnames or [])
        rows = list(reader)
    # Everything after the last metadata column is an injection. Anchoring on one named column
    # instead breaks whenever the writer's optional columns change: with `ms2_support` off, an
    # anchor of `Features Found` silently reads `Dot Product` and `Purity` as samples and the
    # numbers still look plausible. Listing the metadata and taking what follows cannot do that.
    from .peakfinder import META_COLUMNS
    start = max((i for i, c in enumerate(columns) if c in META_COLUMNS), default=-1) + 1
    injections = [c for c in columns[start:] if c.strip()]
    # An empty table still needs a correctly SHAPED array: `np.array([])` is 1-D and every
    # downstream `areas[:, i]` then raises. Shape it (0, n_injections) so the emptiness carries
    # through as "no features" rather than as a different error further along.
    areas = (np.array([[float(r[c] or 0) for c in injections] for r in rows], dtype=float)
             if rows else np.empty((0, len(injections)), dtype=float))

    order = {}
    # ⚠ Prefer the ACQUISITION TIMESTAMP over the sequence file where it can be read.
    #
    # `startTimeStamp` in the mzML is written by the converter from the instrument's own record of
    # when the run started. It is what actually happened, whereas a sequence file is what was
    # planned — and the two diverge whenever a run is repeated, reordered or dropped. It also
    # exists for every file, including polarities whose sequence was never saved: this study has
    # `Pos/Seq_Rand.csv` and nothing for negative at all.
    #
    # The file's mtime is NOT a substitute. On this study every raw file carries an mtime within
    # ten seconds of every other, because they were copied off the instrument in one go.
    if not order and mzml_files:
        stamps = {}
        for path in mzml_files:
            stamp = _acquisition_time(path)
            if stamp:
                stamps[Path(path).stem] = stamp
        if len(stamps) >= max(2, 0.8 * len(mzml_files)):
            for n, stem in enumerate(sorted(stamps, key=lambda k: stamps[k]), start=1):
                for column in injections:
                    if column == stem or re.sub(r"_(Pos|Neg)$", "", column, flags=re.I) == stem:
                        order[column] = n

    order_source = "acquisition timestamp" if order else ""

    # ⚠ When BOTH sources exist, compare them. A disagreement is not a nuisance — it means the
    # sequence and the instrument record different things, which happens when a run is repeated,
    # reordered or dropped mid-sequence. Silently preferring one hides a fact the analyst needs.
    disagreement = ""
    if order and sequence:
        from .blanks import read_sequence
        planned = read_sequence(sequence)
        pairs = []
        for column, n in order.items():
            stem = re.sub(r"_(Pos|Neg)$", "", column, flags=re.I)
            p = planned.get(column, planned.get(stem))
            if p:
                pairs.append((n, p))
        if len(pairs) >= 3:
            ranks = sorted(pairs)
            moved = sum(1 for i, (_, p) in enumerate(ranks, start=1)
                        if p != sorted(x[1] for x in ranks)[i - 1])
            if moved:
                disagreement = (f"{moved} of {len(pairs)} injections sit in a different position "
                                f"in the sequence file than the instrument recorded")

    if not order and sequence:
        from .blanks import read_sequence
        raw_order = read_sequence(sequence)
        # ⚠ Resolve against the columns that actually exist. A sequence may name a sample without
        # its polarity (`AAAA0161_HT`) while the injection is `AAAA0161_HT_Pos` — matching on the
        # bare key alone returns nothing, and every order-dependent section silently reports
        # "not testable" rather than "the sequence did not match".
        for column in injections:
            if column in raw_order:
                order[column] = raw_order[column]
                continue
            stem = re.sub(r"_(Pos|Neg)$", "", column, flags=re.I)
            if stem in raw_order:
                order[column] = raw_order[stem]
        if order:
            order_source = "sequence file"
        if injections and not order:
            print(f"  ⚠ sequence file matched none of the {len(injections)} injections — "
                  f"run order is unavailable, so drift and confounding cannot be tested")

    factors = {}
    if metadata and Path(metadata).exists():
        with Path(metadata).open(newline="") as fh:
            for record in csv.DictReader(fh):
                key = next((v for k, v in record.items()
                            if k.lower() in {"file", "filename", "sample", "injection"}), None)
                if key:
                    factors[key.strip()] = {k: v for k, v in record.items()
                                            if k.lower() not in {"file", "filename", "sample",
                                                                 "injection"}}
    return Run(rows=rows, columns=injections, areas=areas,
               roles={c: infer_role(c) for c in injections}, order=order,
               order_source=order_source, order_disagreement=disagreement, metadata=factors)


# ── shared numerics ────────────────────────────────────────────────────────────────────────

def median_normalise(areas):
    """Divide each injection by its median non-zero area.

    Removes differences in how much material reached the detector, which is what a CV computed
    on raw areas would otherwise be measuring.
    """
    out = areas.astype(float).copy()
    for j in range(out.shape[1]):
        column = out[:, j]
        present = column[column > 0]
        if present.size:
            out[:, j] = column / np.median(present)
    return out


def pareto_log(matrix):
    """log10(x+1), mean-centre, divide by sqrt(sd). Columns with no spread are dropped."""
    logged = np.log10(matrix + 1.0)
    centred = logged - logged.mean(axis=0, keepdims=True)
    spread = logged.std(axis=0, ddof=1)
    usable = spread > 1e-9
    return centred[:, usable] / np.sqrt(spread[usable]), usable


def pca(matrix, components=3):
    u, s, vt = np.linalg.svd(matrix, full_matrices=False)
    k = min(components, len(s))
    variance = (s ** 2) / max(matrix.shape[0] - 1, 1)
    return u[:, :k] * s[:k], variance[:k] / variance.sum(), vt[:k]


def spearman(x, y) -> float:
    from scipy.stats import rankdata
    rx, ry = rankdata(x) - np.mean(rankdata(x)), rankdata(y) - np.mean(rankdata(y))
    denominator = np.sqrt((rx ** 2).sum() * (ry ** 2).sum())
    return float((rx * ry).sum() / denominator) if denominator else 0.0


def _feature_subset(areas, index, min_fraction=MIN_DETECTED_FRACTION):
    sub = areas[:, index]
    detected = (sub > 0).sum(axis=1) / max(len(index), 1)
    return sub[detected >= min_fraction], detected >= min_fraction


def sample_feature_mask(run: "Run", min_fraction=MIN_DETECTED_FRACTION):
    """Which features are analysable, decided on the SAMPLE columns alone.

    ⚠ Deciding it per PCA instead — "detected in half the injections *used here*" — silently
    changes the feature set between panels, and not in a neutral direction. The pooled QCs contain
    nearly everything, so adding eight QC columns hands every sparse feature eight free detections:
    a feature seen in 17 of 40 samples fails the samples-only bar and passes the samples+QC bar at
    17+8 of 48. The samples+QC panel then contains features the samples-only panel correctly
    rejected, and the reader reads the difference as a data problem.

    One mask, computed once, applied to every panel — so the panels differ only in which
    injections are drawn, which is what a reader already assumes.
    """
    index = run.samples
    if not index:
        return np.ones(run.areas.shape[0], dtype=bool)
    detected = (run.areas[:, index] > 0).sum(axis=1) / len(index)
    return detected >= min_fraction


def identified_mask(run: "Run"):
    """Rows carrying a name. The PCA is worth seeing both ways — see `score_space`."""
    return np.array([bool(r.get("Identification", "").strip()) for r in run.rows], dtype=bool)


# ── arm 1: precision ───────────────────────────────────────────────────────────────────────

def profile_screen(run: Run) -> list:
    """Each injection's Spearman against the median sample profile, median-normalised.

    The screen that judges a pool on whether it *looks like* the samples rather than on how much
    of it arrived. Returns one record per injection, with the range the samples themselves span
    so a pool can be read against its own reference.
    """
    normalised = median_normalise(run.areas)
    samples = run.samples
    if not samples:
        return []
    reference = np.median(normalised[:, samples], axis=1)
    usable = reference > 0
    out = []
    for j, column in enumerate(run.columns):
        values = normalised[usable, j]
        out.append({"injection": column, "role": run.roles[column],
                    "order": run.order.get(column),
                    "rho": spearman(reference[usable], values),
                    "total": float(run.areas[:, j].sum()),
                    "detected": int((run.areas[:, j] > 0).sum())})
    sample_rhos = [r["rho"] for r in out if r["role"] == "sample"]
    band = (min(sample_rhos), max(sample_rhos)) if sample_rhos else (0.0, 1.0)
    for record in out:
        record["below_sample_band"] = record["rho"] < band[0]
    return out


def precision(run: Run, mask=None) -> dict:
    """CV on the pools, CV on the samples, and the D-ratio between them.

    `mask` restricts the rows. Computed over every quantified feature these describe the
    instrument; computed over the identified subset they describe **the data the reader will
    actually analyse**, which is usually the number they mean. The two differ, because
    unidentified features are on average fainter and therefore noisier, so the report gives both
    and says which is which rather than quoting one and calling it precision.
    """
    pools, samples = run.pools, run.samples
    def unusable(reason: str) -> dict:
        """The full key set with no numbers in it.

        ⚠ Every early return here must go through this. Returning a bare {} makes callers
        KeyError, and it does so only on studies that are already in trouble — the report for a
        run with zero usable features is exactly when you need it to render and say why. That is
        how `KeyError: 'median_pool_cv'` killed a 32-minute run at the final step, after the fix
        for the OTHER early return had already been made.
        """
        return {"usable": False, "reason": reason,
                "features": 0, "total_features": int(run.areas.shape[0]),
                "pool_cv": np.array([]), "sample_cv": np.array([]), "d_ratio": np.array([]),
                "median_pool_cv": float("nan"), "median_sample_cv": float("nan"),
                "median_d_ratio": float("nan"), "under_20": float("nan"),
                "under_30": float("nan"), "d_under_half": float("nan")}

    if len(pools) < 2 or not samples:
        # Not every study runs a pooled QC series — the skin validation set has blanks and none.
        return unusable(f"{len(pools)} pooled QC injection(s): precision needs at least 2")
    normalised = median_normalise(run.areas)
    if mask is not None:
        mask = np.asarray(mask, dtype=bool)
        if mask.sum() < 3:
            return unusable(f"{int(mask.sum())} feature(s) left after filtering: "
                            f"precision needs at least 3")
        normalised = normalised[mask]
    pool_values, sample_values = normalised[:, pools], normalised[:, samples]
    present = (pool_values > 0).all(axis=1)

    def cv(matrix):
        mean = matrix.mean(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            return 100 * matrix.std(axis=1, ddof=1) / np.where(mean > 0, mean, np.nan)

    pool_cv = cv(pool_values[present])
    sample_cv = cv(sample_values[present])
    with np.errstate(invalid="ignore", divide="ignore"):
        sample_sd = sample_values[present].std(axis=1, ddof=1)
        d_ratio = pool_values[present].std(axis=1, ddof=1) / np.where(sample_sd > 0, sample_sd,
                                                                     np.nan)
    finite = np.isfinite(pool_cv)
    d_finite = np.isfinite(d_ratio)
    return {"usable": True, "reason": "",
            "features": int(present.sum()),
            "total_features": int(normalised.shape[0]),
            "pool_cv": pool_cv[finite], "sample_cv": sample_cv[np.isfinite(sample_cv)],
            "d_ratio": d_ratio[d_finite],
            "median_pool_cv": float(np.median(pool_cv[finite])),
            "median_sample_cv": float(np.median(sample_cv[np.isfinite(sample_cv)])),
            "median_d_ratio": float(np.median(d_ratio[d_finite])),
            "under_20": float((pool_cv[finite] < 20).mean()),
            "under_30": float((pool_cv[finite] < 30).mean()),
            "d_under_half": float((d_ratio[d_finite] < 0.5).mean())}


def score_space(run: Run, keep_roles, mask=None) -> dict | None:
    """PCA over the injections whose role is in `keep_roles`.

    `mask` fixes the feature set from outside, so the same features are used in every panel. Pass
    `sample_feature_mask(run)` for the measurement-level view and `... & identified_mask(run)` for
    the view closest to what will actually be analysed; both are worth showing, because an
    unidentified feature is still a real measurement but will never be reported.
    """
    index = run.indices(*keep_roles)
    if len(index) < 3:
        return None
    if mask is None:
        sub, mask = _feature_subset(run.areas, index)
    else:
        sub = run.areas[np.asarray(mask), :][:, index]
    if sub.shape[0] < 3:
        return None
    scaled, _ = pareto_log(sub.T)
    scores, explained, loadings = pca(scaled)
    return {"columns": [run.columns[i] for i in index],
            "roles": [run.roles[run.columns[i]] for i in index],
            "scores": scores, "explained": explained, "loadings": loadings,
            "scaled": scaled, "features": int(sub.shape[0])}


def contrast(space: dict) -> dict | None:
    """Sample spread against pool spread — whether the biology is bigger than the noise."""
    if space is None:
        return None
    pools = [i for i, r in enumerate(space["roles"]) if r == "qc"]
    samples = [i for i, r in enumerate(space["roles"]) if r == "sample"]
    if not pools or not samples:
        return None
    centre = space["scores"][pools].mean(axis=0)
    pool_spread = float(np.linalg.norm(space["scores"][pools] - centre, axis=1).mean())
    sample_spread = float(np.linalg.norm(space["scores"][samples] - centre, axis=1).mean())
    return {"pool_spread": pool_spread, "sample_spread": sample_spread,
            "ratio": sample_spread / pool_spread if pool_spread else float("inf")}


def outliers(space: dict) -> dict:
    """Hotelling's T-squared and DModX.

    T-squared finds an injection far from the centre *of the model*; DModX finds one the model
    does not describe, which can sit innocently near the centre of a score plot. They catch
    different failures and a report that runs only one of them will miss the other.
    """
    from scipy.stats import f as f_dist

    scores, scaled = space["scores"], space["scaled"]
    n, k = scores.shape
    variance = scores.var(axis=0, ddof=1)
    variance[variance <= 0] = np.inf
    t2 = (scores ** 2 / variance).sum(axis=1)
    limit = (k * (n - 1) / (n - k)) * f_dist.ppf(0.95, k, n - k) if n > k + 1 else np.inf

    # ── robust T-squared: the same distance, but centred and scaled on medians.
    # Classical T² uses the sample covariance, which the outliers themselves inflate, so two or
    # more outliers mask each other — exactly the situation a QC report is looking for. Median and
    # MAD cannot be inflated by a minority, so the distance stays honest as the count rises.
    centre = np.median(scores, axis=0)
    spread = np.median(np.abs(scores - centre), axis=0) * 1.4826
    spread[spread <= 0] = np.inf
    robust = np.sqrt((((scores - centre) / spread) ** 2).sum(axis=1))
    r_med = np.median(robust)
    r_mad = np.median(np.abs(robust - r_med)) * 1.4826
    robust_limit = r_med + 3 * r_mad if r_mad > 0 else np.inf

    # ── distance to the k nearest injections: no covariance model at all, so it fails differently
    # from both T² tests. An injection alone in its region of the space is far from its neighbours
    # however the spread is estimated.
    k = max(2, min(5, n - 2))
    gaps = np.sqrt(((scores[:, None, :] - scores[None, :, :]) ** 2).sum(axis=2))
    np.fill_diagonal(gaps, np.inf)
    neighbour = np.sort(gaps, axis=1)[:, :k].mean(axis=1)
    n_med = np.median(neighbour)
    n_mad = np.median(np.abs(neighbour - n_med)) * 1.4826
    neighbour_limit = n_med + 3 * n_mad if n_mad > 0 else np.inf

    # ── isolation forest, on the SCORES rather than the features.
    # It isolates a point by random axis-parallel splits, so on a raw matrix of several hundred
    # features against forty samples a random split rarely touches whatever makes a sample
    # unusual and the score approaches noise. On the first components it is well within range,
    # and it matches how the other tests are computed. scikit-learn's implementation is used
    # rather than a local one: this number goes into a report that has to be defended.
    forest = np.zeros(n, dtype=bool)
    forest_score = np.zeros(n)
    try:
        from sklearn.ensemble import IsolationForest
        model = IsolationForest(n_estimators=400, random_state=0, contamination="auto")
        model.fit(scores)
        forest_score = -model.score_samples(scores)      # larger is more anomalous
        f_med = np.median(forest_score)
        f_mad = np.median(np.abs(forest_score - f_med)) * 1.4826
        forest = forest_score > (f_med + 3 * f_mad if f_mad > 0 else np.inf)
    except ImportError:
        pass

    residual = scaled - scores @ space["loadings"]
    dmodx = np.sqrt((residual ** 2).sum(axis=1) / max(scaled.shape[1] - k, 1))
    # A distance-to-model is one-sided and skewed, so a robust cut beats a normal-theory limit.
    median = np.median(dmodx)
    mad = np.median(np.abs(dmodx - median)) / 0.6745
    dmodx_limit = median + 3 * mad if mad > 0 else np.inf
    # Five tests, and what matters is how many agree. **A single hit does not flag an injection.**
    # One test crossing a 95% limit is expected by chance in a run of this size, and on a study of
    # fourteen injections that produced exactly the wrong answer: the two most extreme points on
    # the first two components went unflagged while a middling one was flagged, each on one test
    # alone. Flags a reader cannot see are worse than no flags, because they teach the reader to
    # distrust the panel rather than the injection. Single hits are still returned, as `weak`, for
    # anyone who wants them.
    votes = {
        "T\u00b2": t2 > limit,
        "robust T\u00b2": robust > robust_limit,
        "DModX": dmodx > dmodx_limit,
        "neighbour distance": neighbour > neighbour_limit,
        "isolation forest": forest,
    }
    support = {i: [name for name, hit in votes.items() if hit[i]] for i in range(n)}
    return {"t2": t2, "t2_limit": float(limit), "dmodx": dmodx,
            "dmodx_limit": float(dmodx_limit),
            "robust_t2": robust, "robust_limit": float(robust_limit),
            "neighbour": neighbour, "neighbour_limit": float(neighbour_limit),
            "forest": forest_score,
            "support": support,
            "weak": sorted(i for i in range(n) if len(support[i]) == 1),
            "flagged": sorted(i for i in range(n) if len(support[i]) >= 2),
            "corroborated": sorted(i for i in range(n) if len(support[i]) >= 2)}


# ── arm 2: run order ───────────────────────────────────────────────────────────────────────

def run_order(run: Run, rho_cut=DRIFT_RHO, run_in=RUN_IN) -> dict:
    """Per-feature drift against injection number, with and without the first injections.

    A start-of-run effect and a drift across the whole run need different remedies — the first
    is fixed by discarding a few injections, the second by a covariate — so they are separated
    rather than averaged into one number.
    """
    index = [i for i, c in enumerate(run.columns)
             if run.roles[c] == "sample" and c in run.order]
    if len(index) < 8:
        return {}
    injection = np.array([run.order[run.columns[i]] for i in index], float)
    keep_late = injection > np.sort(injection)[min(run_in, len(injection) - 1)]

    per_feature, per_feature_late = [], []
    for row in range(run.areas.shape[0]):
        values = run.areas[row, index]
        if (values > 0).sum() < 0.8 * len(index):
            continue
        per_feature.append(spearman(injection, values))
        if keep_late.sum() >= 8:
            per_feature_late.append(spearman(injection[keep_late], values[keep_late]))
    per_feature = np.array(per_feature)
    per_feature_late = np.array(per_feature_late) if per_feature_late else per_feature

    loading = run.areas[:, index].sum(axis=0)
    return {"tested": len(per_feature), "rho": per_feature, "rho_late": per_feature_late,
            "drifting": int((np.abs(per_feature) >= rho_cut).sum()),
            "drifting_late": int((np.abs(per_feature_late) >= rho_cut).sum()),
            "loading_rho": spearman(injection, loading),
            "loading_rho_late": spearman(injection[keep_late], loading[keep_late])
            if keep_late.sum() >= 8 else float("nan"),
            "injection": injection, "loading": loading,
            "columns": [run.columns[i] for i in index]}


def confounding(run: Run) -> list:
    """Is any biological factor confounded with when the sample was injected?

    Drift only matters if it lines up with the comparison being made. A run that drifts badly but
    is properly randomised is recoverable; a run that drifts mildly while one group sits at the
    start is not.
    """
    from scipy.stats import kruskal, mannwhitneyu

    index = [i for i, c in enumerate(run.columns)
             if run.roles[c] == "sample" and c in run.order and c in run.metadata]
    if len(index) < 8:
        return []
    out = []
    factors = usable_factors({run.columns[i]: run.metadata[run.columns[i]] for i in index})
    for factor in factors:
        groups = {}
        for i in index:
            value = (run.metadata.get(run.columns[i], {}).get(factor) or "").strip()
            if value:
                groups.setdefault(value, []).append(run.order[run.columns[i]])
        groups = {k: v for k, v in groups.items() if len(v) >= 3}
        if len(groups) < 2:
            continue
        arrays = [np.array(v, float) for v in groups.values()]
        try:
            if len(arrays) == 2:
                p = float(mannwhitneyu(arrays[0], arrays[1], alternative="two-sided").pvalue)
            else:
                p = float(kruskal(*arrays).pvalue)
        except ValueError:
            continue
        out.append({"factor": factor,
                    "levels": {k: (len(v), float(np.median(v))) for k, v in groups.items()},
                    "p": p, "confounded": p < 0.05})
    return out


def pool_drift(run: Run) -> dict:
    """How much each lipid moves across the pooled QCs — the input a QC-RSC decision needs.

    QC-RSC corrects a run by fitting each feature's trend through the pooled QCs and dividing it
    out. Whether to do that is a threshold decision this facility has not yet made, so nothing is
    corrected here: what is reported is the size of the trend that a correction would remove, at
    several thresholds, so the threshold can be chosen against real numbers rather than a default.

    Two measures, because they disagree in a useful way. The Spearman against pool order says
    whether the movement is *systematic*; the first-to-last fold says whether it is *large*. A
    feature can drift monotonically by 5% (systematic, harmless) or jump twofold at random
    (large, not a drift).
    """
    pools = [(run.order.get(run.columns[i], i), i) for i in run.pools]
    pools.sort()
    if len(pools) < 4:
        return {}
    index = [i for _, i in pools]
    position = np.arange(len(index), dtype=float)
    normalised = median_normalise(run.areas)[:, index]
    usable = (normalised > 0).all(axis=1)
    values = normalised[usable]
    rho = np.array([spearman(position, row) for row in values])
    fold = values[:, -1] / np.where(values[:, 0] > 0, values[:, 0], np.nan)
    with np.errstate(invalid="ignore"):
        log_fold = np.log2(fold)
    finite = np.isfinite(log_fold)
    grid = []
    for cut in (0.4, 0.6, 0.8):
        for lf in (0.5, 1.0):
            grid.append({"rho": cut, "log2_fold": lf,
                         "features": int(((np.abs(rho) >= cut) &
                                          finite & (np.abs(log_fold) >= lf)).sum())})
    return {"pools": [run.columns[i] for i in index],
            "order": [o for o, _ in pools], "features": int(usable.sum()),
            "rho": rho, "log_fold": log_fold[finite],
            "median_abs_rho": float(np.median(np.abs(rho))),
            "median_abs_log_fold": float(np.median(np.abs(log_fold[finite]))),
            "systematic": int((np.abs(rho) >= 0.6).sum()),
            "grid": grid,
            "totals": [float(run.areas[:, i].sum()) for i in index]}


# ── arm 3: annotation ambiguity ────────────────────────────────────────────────────────────

# Strips a TRAILING ADDUCT only — `[M+H]+`, `[M+FA-H]-`, `[M-2H]2-`.
# ⚠ It must not be `\s*\[.*$`. That also ate the subclass bracket that is part of the molecule's
# name: `Cer[NS] d34:1` became `Cer`, `GlcCer[AP] t42:0` became `GlcCer`, `PC[OH] 36:2` became
# `PC`. Every sphingolipid and hydroxy subclass therefore collapsed onto its class, which inflated
# the repeated-name count and made the report's worked example print a class where a molecule
# should be.
_ADDUCT = re.compile(r"\s*\[M[^\]]*\][+-]\d*\s*;?\s*$")


def annotation_ambiguity(run: Run, rt_window=None, mz_window=0.01,
                         fwhm=None, min_correlation=0.8) -> dict:
    """Where one measurement is being reported as more than one lipid.

    Three distinct failures, and they are not interchangeable:

      * **one name, several retention times** — either real isomers, which must stay separate,
        or one compound split across features, which must not;
      * **several names, one peak** — the same measurement counted twice, which doubles the
        apparent evidence in any class or enrichment summary;
      * **a summed composition alongside its own resolved form** — `PC 34:1` and
        `PC 16:0_18:1` at the same peak are one molecule reported at two resolutions.

    ⚠ **The retention window must be scaled to the measured peak width, and correlation must be a
    criterion rather than a printed decoration.** A fixed 0.2 min window is more than two peak
    widths on this method (FWHM 0.077 min negative, 0.091 positive), and every pair it produced was
    1.1-1.8 widths apart with areas correlating between -0.20 and +0.53 — where one measurement
    reported twice would correlate at about +1. They were adjacent but distinct peaks, and all 15
    were artefacts of the window.

    ⚠ **And the conceptual limit stands even when a pair is real:** a sum composition can
    correspond to several molecular species, so pairing `PC 34:1` with `PC 16:0_18:1` assumes that
    resolved form is the right one of the possibilities. The pairing can say two rows are one
    measurement. It can never confirm the chain assignment.
    """
    # Default to one peak width. The peak finder measures FWHM per run, and it differs between
    # methods and between polarities, so a constant in minutes cannot be right for both.
    if rt_window is None:
        rt_window = fwhm if fwhm else 0.09
    identified = [(i, r) for i, r in enumerate(run.rows) if r.get("Identification", "").strip()]
    by_name: dict = {}
    for i, row in identified:
        name = row["Identification"].strip()
        by_name.setdefault(_ADDUCT.sub("", name), []).append(
            (i, float(row["Retention Time (min)"]), float(row["Quant Ion"]), name))

    repeated = []
    for name, entries in sorted(by_name.items()):
        if len(entries) < 2:
            continue
        retentions = sorted(e[1] for e in entries)
        repeated.append({"name": name, "rows": len(entries),
                         "retentions": retentions,
                         "spread": retentions[-1] - retentions[0],
                         "areas": [float(run.rows[e[0]]["Area (max)"]) for e in entries]})

    order = sorted(identified, key=lambda t: (float(t[1]["Retention Time (min)"]),
                                              float(t[1]["Quant Ion"])))
    same_peak, summed_pairs, rejected = [], [], []
    for a in range(len(order)):
        ia, ra = order[a]
        for b in range(a + 1, len(order)):
            ib, rb = order[b]
            if float(rb["Retention Time (min)"]) - float(ra["Retention Time (min)"]) > rt_window:
                break
            if abs(float(rb["Quant Ion"]) - float(ra["Quant Ion"])) > mz_window:
                continue
            na = _ADDUCT.sub("", ra["Identification"].strip())
            nb = _ADDUCT.sub("", rb["Identification"].strip())
            if na == nb:
                continue
            r = _correlation(run, ia, ib)
            record = {"a": na, "b": nb,
                      "rt": (float(ra["Retention Time (min)"]),
                             float(rb["Retention Time (min)"])),
                      "mz": float(ra["Quant Ion"]),
                      "gap_widths": ((float(rb["Retention Time (min)"])
                                      - float(ra["Retention Time (min)"])) / rt_window),
                      "correlation": r}
            # Two rows are one measurement only if they behave like one. Without this the list is
            # a list of neighbours.
            if r is None or r < min_correlation:
                rejected.append(record)
                continue
            if _is_resolved_form(na, nb) or _is_resolved_form(nb, na):
                summed_pairs.append(record)
            else:
                same_peak.append(record)
    # The surplus is the number that matters and the one a reader cannot derive: rows occupied by
    # names that appear more than once, minus one row each for the names themselves.
    multi = [r for r in repeated if r["rows"] > 1]
    occupied = sum(r["rows"] for r in multi)
    return {"identified_rows": len(identified),
            "distinct_names": len(by_name),
            "repeated_names": len(multi),
            "rows_occupied": occupied,
            "surplus": occupied - len(multi),
            "single_row_names": len(by_name) - len(multi),
            "rt_window": rt_window,
            "repeated": sorted(repeated, key=lambda r: -r["rows"]),
            "same_peak": sorted(same_peak, key=lambda r: -(r["correlation"] or 0)),
            "summed_pairs": summed_pairs,
            "rejected_pairs": rejected}


def _is_resolved_form(summed: str, molecular: str) -> bool:
    """`PC 34:1` against `PC 16:0_18:1` — the same molecule at two resolutions."""
    from .peaks import sum_composition
    return "_" in molecular and "_" not in summed and sum_composition(molecular) == summed


def _correlation(run: Run, i: int, j: int) -> float | None:
    samples = run.samples
    if len(samples) < 4:
        return None
    a, b = run.areas[i, samples], run.areas[j, samples]
    if (a > 0).sum() < 4 or (b > 0).sum() < 4:
        return None
    return spearman(a, b)


# ── arm 4: blanks ──────────────────────────────────────────────────────────────────────────

def blank_report(run: Run) -> list:
    """Per blank: how much is in it, and how much of that also sits in the samples."""
    samples = run.samples
    if not samples or not run.blank_columns:
        return []
    sample_mean = run.areas[:, samples].mean(axis=1)
    out = []
    for j in run.blank_columns:
        values = run.areas[:, j]
        present = values > 0
        total = float(values[present].sum())
        shared = float(values[present & (sample_mean > values)].sum())
        out.append({"injection": run.columns[j], "order": run.order.get(run.columns[j]),
                    "features": int(present.sum()), "total": total,
                    "shared": (shared / total) if total else 0.0})
    return sorted(out, key=lambda r: (r["order"] is None, r["order"]))


# ── arm 5: coverage ────────────────────────────────────────────────────────────────────────

def coverage(run: Run) -> dict:
    """What was measured, and on what evidence."""
    from .peaks import sum_composition

    identified = [r for r in run.rows if r.get("Identification", "").strip()]
    classes: dict = {}
    for row in identified:
        classes.setdefault(row.get("Lipid Class", "").strip() or "?", set()).add(
            sum_composition(row["Identification"].strip()))
    sources: dict = {}
    for row in identified:
        sources[row.get("Identification Source", "").strip() or "MS2"] = \
            sources.get(row.get("Identification Source", "").strip() or "MS2", 0) + 1
    dots = [int(r["Dot Product"]) for r in identified if r.get("Dot Product", "").strip()]
    return {"rows": len(run.rows), "identified": len(identified),
            "molecules": len({sum_composition(r["Identification"].strip())
                              for r in identified}),
            "classes": {k: len(v) for k, v in sorted(classes.items(),
                                                     key=lambda kv: -len(kv[1]))},
            "sources": sources,
            "dot_median": float(np.median(dots)) if dots else None,
            "dot_low": int(sum(1 for d in dots if d < 700))}


# ── arm 6: missing values ──────────────────────────────────────────────────────────────────

def missingness(run: Run) -> dict:
    """How much is missing, and — the part that decides the treatment — why.

    Counting is not enough. The three mechanisms have different remedies, and picking the wrong
    one biases the result rather than merely weakening it:

        MCAR   missing independently of everything     any imputation is unbiased; rare in MS
        MAR    explained by something OBSERVED         correctable by conditioning on it
        MNAR   depends on the unobserved value itself  left-censored; mean/median biases upward

    So four questions are asked in order: does missingness track abundance (MNAR), does it track a
    technical covariate (MAR, correctable), does it track a biological factor (the dangerous one —
    a comparison then partly compares detection rather than biology), and what does a formal MCAR
    test say.

    The headline needs framing or it reads as a disaster: an untargeted table is mostly sparse
    unidentified features, and the number a reader should carry is the identified subset.
    """
    samples = run.samples
    if not samples:
        # ⚠ The full key set with no numbers in it, never a bare {}.
        #
        # This is the THIRD time this shape has bitten: `precision` had two such returns and only
        # one had been fixed, and the report then died on `missing_identified` in an overnight
        # batch — after the pipeline had written every result file, so the data was complete and
        # the run was still marked FAIL. A report for a study with no usable samples is exactly
        # when it must render and say why.
        return {"usable": False, "reason": "no sample injections", "injections": 0, "features": 0,
                "missing_all": float("nan"), "missing_identified": float("nan"),
                "complete_all": 0, "complete_identified": 0, "identified": 0}
    areas = run.areas[:, samples]
    present = areas > 0
    named = identified_mask(run)

    out = {
        "injections": len(samples),
        "features": int(areas.shape[0]),
        "missing_all": float(1 - present.mean()),
        "missing_identified": float(1 - present[named].mean()) if named.any() else float("nan"),
        "complete_all": int((present.all(axis=1)).sum()),
        "complete_identified": int((present[named].all(axis=1)).sum()) if named.any() else 0,
        "identified": int(named.sum()),
        "sparse": int((present.mean(axis=1) < 0.25).sum()),
    }

    # 1 ─ against abundance. The signature of left-censoring, and the dominant mechanism in MS.
    mean_when_seen = np.array([row[row > 0].mean() if (row > 0).any() else 0.0 for row in areas])
    usable = mean_when_seen > 0
    quartiles = []
    if usable.sum() >= 8:
        edges = np.quantile(mean_when_seen[usable], [0, 0.25, 0.5, 0.75, 1.0])
        for lo, hi in zip(edges[:-1], edges[1:]):
            band = usable & (mean_when_seen >= lo) & (mean_when_seen <= hi)
            if band.any():
                quartiles.append(float(1 - present[band].mean()))
    out["by_abundance"] = quartiles
    out["abundance_gap"] = (quartiles[0] - quartiles[-1]) if len(quartiles) == 4 else float("nan")

    # 2 ─ against technical covariates. A relationship here is MAR, and correctable.
    per_injection = 1 - present.mean(axis=0)
    out["per_injection"] = {run.columns[c]: float(v) for c, v in zip(samples, per_injection)}
    if run.order:
        seq = [run.order.get(run.columns[c], 0) for c in samples]
        if len(set(seq)) > 2:
            out["vs_order"] = spearman(seq, per_injection)
    totals = areas.sum(axis=0)
    out["vs_loading"] = spearman(totals, per_injection)

    # 3 ─ against the biological factors. Must be reported even when negative: if one group is
    # systematically less complete, a downstream comparison is partly comparing detection.
    factors = {}
    for column in {k for c in samples for k in run.metadata.get(run.columns[c], {})}:
        groups = {}
        for position, c in enumerate(samples):
            value = run.metadata.get(run.columns[c], {}).get(column)
            if value:
                groups.setdefault(value, []).append(per_injection[position])
        if len(groups) >= 2 and all(len(v) >= 3 for v in groups.values()):
            levels = {k: float(np.mean(v)) for k, v in groups.items()}
            factors[column] = {"levels": levels,
                               "spread": max(levels.values()) - min(levels.values()),
                               "p": _kruskal([groups[k] for k in groups])}
    out["by_factor"] = factors

    # 4 ─ a formal MCAR test, reported with its caveat rather than as a verdict.
    out["mcar_p"] = _little_like(present)
    out["mechanism"] = _mechanism(out)
    return out


def _kruskal(groups) -> float:
    try:
        from scipy.stats import kruskal
        return float(kruskal(*groups).pvalue)
    except Exception:
        return float("nan")


def _little_like(present) -> float:
    """A chi-square on whether missingness is homogeneous across injections.

    Not Little's test proper — that needs the values as well as the pattern — and on a table this
    size it is near-certain to reject whatever the mechanism is. Reported for completeness, with
    the caveat that rejection alone does not separate MAR from MNAR.
    """
    try:
        from scipy.stats import chi2_contingency
        counts = np.vstack([(~present).sum(axis=0), present.sum(axis=0)])
        counts = counts[:, counts.sum(axis=0) > 0]
        if counts.shape[1] < 2 or (counts.sum(axis=1) == 0).any():
            return float("nan")
        return float(chi2_contingency(counts).pvalue)
    except Exception:
        return float("nan")


def _mechanism(summary: dict) -> str:
    """The dominant mechanism, in one word, with the evidence behind it in the report."""
    gap = summary.get("abundance_gap", float("nan"))
    if gap == gap and gap > 0.15:
        return "MNAR"
    order = abs(summary.get("vs_order", 0.0) or 0.0)
    if order > 0.5:
        return "MAR"
    return "inconclusive"


def blank_ratio(run: Run, multiplier: float = 3.0) -> dict:
    """How far the samples stand above the blanks, per feature.

    "Shared with the blanks" is not a useful statement on its own: a feature 1000x higher in the
    samples than in a blank is background that does not matter, and one 5x higher is background
    that does. Both are "shared". What decides it is the RATIO, and its distribution is the thing
    worth showing.

    Returns log10(mean sample / blank level) for every feature detected in a blank, so the
    multiplier the filter uses can be drawn on the same axis as the data it acts on.
    """
    samples, blanks = run.samples, run.blank_columns
    if not samples or not blanks:
        return {}
    blank_level = run.areas[:, blanks].mean(axis=1)
    sample_mean = run.areas[:, samples].mean(axis=1)
    seen = (blank_level > 0) & (sample_mean > 0)
    if not seen.any():
        return {}
    ratio = sample_mean[seen] / blank_level[seen]
    named = identified_mask(run)[seen]
    return {"log_ratio": np.log10(ratio), "identified": named,
            "features": int(seen.sum()),
            "below": float((ratio < multiplier).mean()),
            "below_identified": (float((ratio[named] < multiplier).mean())
                                 if named.any() else float("nan")),
            "median": float(np.median(ratio)),
            "multiplier": float(multiplier)}


def mass_accuracy(spectra_path, order: "dict[str, int] | None" = None,
                  roles: "dict[str, str] | None" = None) -> dict:
    """Instrument mass accuracy per injection, from the rank-1 library matches.

    Uses the identifications rather than the spiked standards because there are tens of thousands
    of them against a handful of standards, which is what makes a per-injection figure possible at
    all. The standards are the independent check on this number, not its source — see
    `standards.check_offset`.

    A systematic offset is calibration and is corrected with `mass_offset_ppm`; a *trend across
    injections* is the mass axis moving during the run, which no single constant can fix, and is
    the reason this is reported per injection rather than as one number.

    **Samples and pooled QCs.** A blank has almost no analyte, and its matches are too few and too
    poorly centroided to compare with a sample's. Standards are left out too: the mix is deuterated
    and largely absent from the library, so it produces few matches and appears in one polarity but
    not the other — an asymmetry in what the library can name, arriving on the page as though it
    were an asymmetry in the instrument. The standards' own mass accuracy is reported against their
    known masses in `Standards_Performance_<polarity>.csv`, which is a better measurement of them
    than a library match is.
    """
    import csv
    from collections import defaultdict
    from pathlib import Path

    if not spectra_path or not Path(spectra_path).exists():
        return {}
    from .blanks import infer_role

    wanted = {"sample", "qc"}
    seen_role: dict[str, str] = {}
    per = defaultdict(list)
    for row in csv.DictReader(open(spectra_path)):
        try:
            if int(row["Rank"]) != 1:
                continue
            # ⚠ Fall back to inferring the role, never to assuming "sample". The delivered table
            # no longer carries standards columns, so the roles map built from it does not mention
            # them at all — and a default of "sample" then let every standards injection straight
            # back into a panel they had just been excluded from.
            if row["Sample"] not in seen_role:
                declared = (roles or {}).get(row["Sample"])
                seen_role[row["Sample"]] = declared or infer_role(row["Sample"])
            if seen_role[row["Sample"]] not in wanted:
                continue
            per[row["Sample"]].append(float(row["Delta m/z (ppm)"]))
        except (KeyError, ValueError):
            continue
    if not per:
        return {}

    injections = []
    for name, deltas in per.items():
        if len(deltas) < 20:          # too few to give a stable median
            continue
        injections.append({"injection": name, "n": len(deltas),
                           "median_ppm": float(np.median(deltas)),
                           "iqr_ppm": float(np.percentile(deltas, 75) - np.percentile(deltas, 25)),
                           "order": (order or {}).get(name)})
    if not injections:
        return {}
    injections.sort(key=lambda r: (r["order"] is None, r["order"], r["injection"]))
    medians = np.array([r["median_ppm"] for r in injections])

    out = {"injections": injections,
           "median_ppm": float(np.median(medians)),
           "spread_ppm": float(medians.max() - medians.min()),
           "n_matches": int(sum(r["n"] for r in injections))}

    placed = [r for r in injections if r["order"] is not None]
    if len(placed) >= 6:
        x = np.array([r["order"] for r in placed], dtype=float)
        y = np.array([r["median_ppm"] for r in placed])
        out["drift_ppm"] = float(np.polyfit(x, y, 1)[0] * (x.max() - x.min()))
        out["rho"] = spearman(x, y)
        out["first_ppm"], out["last_ppm"] = float(y[0]), float(y[-1])
    return out
