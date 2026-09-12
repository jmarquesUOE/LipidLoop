"""Run the peak finder over a quantified feature table plus per-file identifications.

`CDPeakFinder.java :: runQuantitation`, in order. The feature table can come from Compound
Discoverer's export or from `features.py`; both are normalised to the same three-level
compound group / compound / feature model before anything here runs, so this stage never
learns which produced it.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from .adducts import (Adduct, check_adduct_known, check_adduct_unknown, check_dimer,
                      check_fragment, load_adducts)
from .correlate import CorrelationFilter, correlate_areas, filter_correlated_unidentified
from .rtls import fit_models, parse_sum_name, retention_error
from .lsi import shorthand
from .nomenclature import canonical_class
from .peaks import (ISOTOPE_SPACING, MIN_DOT_PRODUCT, MIN_FA_PURITY, MIN_ID_NUM,
                    MIN_REV_DOT_PRODUCT, MIN_RT_MULTIPLIER, UNALIGNED_PPM, Compound,
                    CompoundGroup, Feature, Lipid, Sample, ppm_diff)

ADDUCT_TOLERANCE = 0.01

# Two constants quoted from the original rather than chosen. 0.823 is a FWHM *divisor* — two
# rows reporting the same lipid are one peak within avg_fwhm/0.823, about 1.2 peak widths.
# 0.6745 rescales a median absolute deviation to a standard deviation for a normal distribution.
REDUNDANT_FWHM_DIVISOR = 0.823
MAD_TO_SIGMA = 0.6745

# Classes whose spectra are far more informative in one polarity than the other. LipiDex takes
# this from the library's own OptimalPolarity flag, which the result file carries per row.
PREFERRED_POLARITY_COLUMN = "Optimal Polarity"


@dataclass(slots=True)
class PeakFinderResult:
    compound_groups: list[CompoundGroup]
    samples: list[Sample]

    @property
    def kept(self) -> list[CompoundGroup]:
        return [g for g in self.compound_groups if g.keep]

    @property
    def identified(self) -> list[CompoundGroup]:
        return [g for g in self.kept if g.final_lipid_id is not None]


class PeakFinder:
    def __init__(self, compound_groups: list[CompoundGroup], samples: list[Sample],
                 rt_filter: bool = True, rt_filter_multiplier: float = MIN_RT_MULTIPLIER,
                 adduct_filtering: bool = True, in_source_filtering: bool = True,
                 min_feature_count: int = MIN_ID_NUM, adducts: list[Adduct] | None = None,
                 retention_model_filter: bool = True,
                 retention_model_sigma: float = 3.0,
                 retention_model_max_error: float = 1.0,
                 correlation: CorrelationFilter | None = None):
        self.compound_groups = compound_groups
        self.samples = samples
        self.rt_filter = rt_filter
        self.rt_filter_multiplier = rt_filter_multiplier
        self.adduct_filtering = adduct_filtering
        self.in_source_filtering = in_source_filtering
        self.min_feature_count = min_feature_count
        # Without an adduct database the sweep cannot ask whether two quant ions are the same
        # neutral molecule, and degenerates to comparing raw mass differences.
        self.adducts = adducts or []
        self.retention_model_filter = retention_model_filter
        self.retention_model_sigma = retention_model_sigma
        self.retention_model_max_error = retention_model_max_error
        self.retention_models: dict = {}
        self.correlation = correlation or CorrelationFilter()
        self._correlation_columns = [i for i, s in enumerate(samples) if s.role != "blank"]
        self.avg_fwhm = 0.0
        self.imported_lipids: list[Lipid] = []
        self.unassigned_lipids: list[Lipid] = []

    # ── pipeline

    def run(self, result_files: dict[str, Path], log=None) -> PeakFinderResult:
        say = log or (lambda _: None)

        self.compound_groups.sort(key=lambda g: g.retention)
        for group in self.compound_groups:
            group.calc_fwhm()
            group.find_quant_ion()
        self.avg_fwhm = (sum(g.avg_fwhm for g in self.compound_groups)
                         / max(len(self.compound_groups), 1))
        say(f"{len(self.compound_groups)} compound groups, mean FWHM {self.avg_fwhm:.3f} min")

        self._calculate_rt_deviations()
        self._fill_rt_gaps()
        self._find_isobaric_features()

        self.imported_lipids = self._load_ids(result_files)
        say(f"{len(self.imported_lipids)} identifications loaded")

        associated = self._associate_ids()
        passing = sum(1 for l in self.imported_lipids
                      if l.dot > MIN_DOT_PRODUCT and l.rev_dot > MIN_REV_DOT_PRODUCT)
        say(f"{passing} identifications pass the score filters, "
            f"{associated} associated to features")

        self._score_and_filter_by_dot_product()
        for group in self.compound_groups:
            group.filter_sum_ids()
            group.calculate_weighted_purity()
            group.best_lipid_id()

        self._filter_results()
        if self.rt_filter:
            self._check_class_rt_distribution(self.rt_filter_multiplier)
        if self.retention_model_filter:
            self._filter_by_retention_model(say)

        # Last, because it needs the final identifications: an unidentified feature that
        # co-elutes with an identified one and tracks it across every sample is the same peak
        # seen twice. Only unidentified groups can be removed, so no identification is at risk.
        removed = filter_correlated_unidentified(self.compound_groups, self.samples,
                                                 self.correlation)
        if removed and log:
            say(f"{removed} unidentified features removed as correlated with an identified peak")

        for group in self.compound_groups:
            for compound in group.compounds:
                compound.sample.compounds_found += 1
                compound.sample.fwhms.append(compound.fwhm)

        result = PeakFinderResult(self.compound_groups, self.samples)
        say(f"{len(result.kept)} groups kept, {len(result.identified)} identified")
        return result

    # ── retention-time alignment

    def _calculate_rt_deviations(self):
        """Deviation between each feature's measured and aligned retention, then LOESS per sample."""
        for group in self.compound_groups:
            for compound in group.compounds:
                for feature in compound.features:
                    if feature.real_retention != 0.0:
                        feature.sample.rt_deviations.append(
                            (feature.real_retention - feature.retention,
                             feature.retention, feature.area))
        for sample in self.samples:
            sample.fit()

    def _fill_rt_gaps(self):
        """Features the unaligned table never matched get their measured retention interpolated."""
        for group in self.compound_groups:
            for compound in group.compounds:
                for feature in compound.features:
                    if feature.real_retention == 0.0:
                        offset = feature.sample.corrected_rt(feature.retention)
                        feature.real_retention = (feature.retention
                                                  if offset == 0.0
                                                  else feature.retention + (feature.retention - offset))

    def _find_isobaric_features(self):
        """Neighbouring features of the same mass and polarity compete to own an MS2.

        Without this a spectrum lands on whichever isobaric peak is checked first rather than
        the one actually eluting at that moment.
        """
        by_sample: dict[str, list[Feature]] = {}
        for group in self.compound_groups:
            for compound in group.compounds:
                for feature in compound.features:
                    by_sample.setdefault(feature.sample.file, []).append(feature)

        for features in by_sample.values():
            features.sort(key=lambda f: f.mass)
            for i, feature in enumerate(features):
                for j in range(i + 1, len(features)):
                    other = features[j]
                    if ppm_diff(other.mass, feature.mass) >= UNALIGNED_PPM:
                        break
                    if other.polarity != feature.polarity:
                        continue
                    if abs(other.retention - feature.retention) < self.avg_fwhm * 6.0:
                        feature.isobaric_neighbours.append(other)
                        other.isobaric_neighbours.append(feature)
        for group in self.compound_groups:
            group.no_coeluting_peaks = not any(
                f.isobaric_neighbours for c in group.compounds for f in c.features)

    # ── identifications

    def _load_ids(self, result_files: dict[str, Path]) -> list[Lipid]:
        lipids: list[Lipid] = []
        by_file = {s.file: s for s in self.samples}
        for sample_name, path in result_files.items():
            sample = by_file.get(sample_name)
            if sample is None:
                continue
            with Path(path).open(newline="") as fh:
                for row in csv.DictReader(fh):
                    try:
                        lipid = Lipid(
                            retention=float(row["Retention Time (min)"]),
                            precursor=float(row["Precursor Mass"]),
                            sample=sample,
                            dot=float(row["Dot Product"]),
                            rev_dot=float(row["Reverse Dot Product"]),
                            lipid_string=row["Identification"],
                            lib_precursor=float(row["Library Mass"]),
                            purity=int(row["Purity"] or 0),
                            is_lipidex=row["LipiDex Spectrum"].strip().lower() == "true",
                            purity_array=parse_purity(row["Spectral Components"]),
                            fragment_masses=parse_masses(row["Potential Fragments"]),
                            # Optional, so a hand-made or third-party `_Results.csv` still loads;
                            # only `Associated_Spectra` needs them.
                            ms2_id=int(row.get("MS2 ID") or 0),
                            rank=int(row.get("Rank") or 0))
                    except (ValueError, KeyError):
                        continue
                    lipid.preferred_polarity = row[PREFERRED_POLARITY_COLUMN].strip().lower() == "true"
                    sample.purities.append(lipid.purity)
                    lipids.append(lipid)
        return lipids

    def _associate_ids(self) -> int:
        """Attach each identification to the feature it was fragmented from."""
        associated = 0
        index = RetentionIndex(self.compound_groups, self.avg_fwhm)
        for lipid in self.imported_lipids:
            if not (lipid.dot > MIN_DOT_PRODUCT and lipid.rev_dot > MIN_REV_DOT_PRODUCT):
                continue
            found = False
            for group in index.near(lipid.corrected_retention):
                if lipid.polarity == "+" and not group.positive_feature:
                    continue
                if lipid.polarity == "-" and not group.negative_feature:
                    continue
                if group.is_unique_quant_ion(lipid.precursor):
                    continue
                for compound in group.compounds:
                    if abs(lipid.corrected_retention - compound.retention) > 1.0:
                        break
                    for feature in compound.features:
                        if feature.accepts(lipid, group.no_coeluting_peaks):
                            feature.lipids.append(lipid)
                            group.add_sum_id(lipid.sum_lipid_name)
                            feature.sample.lipids += 1
                            found = True
                            break
                    if found:
                        break
                if found:
                    break
            if found:
                associated += 1
            else:
                self.unassigned_lipids.append(lipid)
        return associated

    def _score_and_filter_by_dot_product(self):
        """Drop candidates below the score thresholds and score where on the peak each sat.

        `CDFeature :: filterIDByDotProduct` removes only what falls *below* the thresholds — it
        does not keep only the best-scoring candidate. Every survivor goes on to vote in the
        sum-composition majority and to contribute to the weighted purity, which is the point:
        a peak with several credible identities should come out impure and be reported at sum
        composition, and it cannot if all but one candidate has been thrown away first.
        """
        for group in self.compound_groups:
            for compound in group.compounds:
                for feature in compound.features:
                    for lipid in feature.lipids:
                        if lipid.dot < MIN_DOT_PRODUCT or lipid.rev_dot < MIN_REV_DOT_PRODUCT:
                            lipid.keep = False
                            continue
                        lipid.gaussian_score = feature.peak_model.normalized_height(
                            lipid.corrected_retention)

    # ── filtering

    def _filter_results(self):
        """Remove isotopes, adducts, in-source fragments, dimers and redundant identifications.

        Walks backwards from each group over the earlier groups still within 2 FWHM, which is
        why the list has to be retention-sorted before this runs.
        """
        groups = self.compound_groups
        for i in range(1, len(groups)):
            current = groups[i]
            j = i - 1
            while j >= 0 and current.target_in_fwhm(groups[j], 2.0):
                other = groups[j]
                self._filter_pair(current, other)
                j -= 1

    def _tracks(self, a: CompoundGroup, b: CompoundGroup) -> bool:
        """Do two groups rise and fall together across the samples?

        Used to confirm a mass relationship before acting on it. Two ions of one molecule are
        the same peak split, so they cannot disagree across samples; two coincidentally
        adduct-spaced compounds easily can. Returns True when there is not enough data to
        judge — the mass evidence then stands on its own rather than being overruled by silence.
        """
        if not self.correlation.guard_mass_relationships:
            return True
        r = correlate_areas(a.areas, b.areas, self._correlation_columns,
                            self.correlation.correlation_type, self.correlation.min_points)
        return True if r is None else r >= self.correlation.guard_min_correlation

    def _filter_pair(self, current: CompoundGroup, other: CompoundGroup):
        """The pairwise sweep, in `CDPeakFinder :: filterResults` order.

        Isotope, then adducts against an identified peak, then adducts between two unidentified
        peaks, then in-source fragments, then dimers, then redundant identifications. Order
        matters: an unidentified peak removed as an adduct is no longer available to be removed
        as a dimer, and the reason recorded is the first one that applied.
        """
        if current.quant_ion is None or other.quant_ion is None:
            return
        close = current.target_in_fwhm(other, 0.5)
        gap = abs(current.quant_ion - other.quant_ion)

        # ── M+1 isotope of an unidentified neighbour
        if (abs(ISOTOPE_SPACING - gap) < ADDUCT_TOLERANCE and close
                and current.quant_polarity == other.quant_polarity):
            heavier = current if current.quant_ion > other.quant_ion else other
            if heavier.final_lipid_id is None and heavier.keep and self._tracks(current, other):
                heavier.keep = False
                heavier.filter_reason = "M+1 isotope"

        if self.adduct_filtering and close and current.keep and other.keep:
            # ── adduct of an identified peak: the identification fixes one neutral mass
            current_adduct = _adduct_of(current)
            other_adduct = _adduct_of(other)
            if (check_adduct_known(self.adducts, current.quant_ion, other.quant_ion,
                                   current_adduct)
                    or check_adduct_known(self.adducts, other.quant_ion, current.quant_ion,
                                          other_adduct)):
                if current.final_lipid_id is None and other.final_lipid_id is not None \
                        and self._tracks(current, other):
                    current.keep = False
                    current.filter_reason = "Adduct of existing identified peak"
                elif current.final_lipid_id is not None and other.final_lipid_id is None \
                        and self._tracks(current, other):
                    other.keep = False
                    other.filter_reason = "Adduct of existing identified peak"

            # ── adduct between two unidentified peaks: keep the larger
            if (current.keep and other.keep
                    and current.final_lipid_id is None and other.final_lipid_id is None
                    and check_adduct_unknown(self.adducts, current.quant_ion, other.quant_ion)
                    and self._tracks(current, other)):
                weaker = current if current.max_area < other.max_area else other
                weaker.keep = False
                weaker.filter_reason = "Adduct of existing peak"

        # ── in-source fragment: one group's quant ion is another's predicted fragment.
        # No `keep` guard, deliberately: the original tests this whatever has already happened
        # to the pair, so a group removed as an adduct can still be re-marked here. The reason
        # a row carries is therefore the LAST rule that applied to it, not the first.
        if self.in_source_filtering and close:
            if check_fragment(_fragments_of(current), other.quant_ions) \
                    and current.quant_polarity == other.quant_polarity:
                other.keep = False
                other.filter_reason = "In-source fragment"
            elif check_fragment(_fragments_of(other), current.quant_ions) \
                    and current.quant_polarity == other.quant_polarity:
                current.keep = False
                current.filter_reason = "In-source fragment"

        # ── dimer: remove the unidentified one. No `keep` guard, as above.
        if (self.adduct_filtering and close
                and check_dimer(self.adducts, current.quant_ion, other.quant_ion,
                                current.quant_polarity)
                and self._tracks(current, other)):
            heavier = current if current.quant_ion > other.quant_ion else other
            if heavier.final_lipid_id is None:
                heavier.keep = False
                heavier.filter_reason = "Dimer"

        # ── same lipid reported twice: pool the evidence into the larger peak. Matched on SUM
        # composition, not the displayed name, so `PC 34:1` and `PC 16:0_18:1` collapse.
        if (current.final_lipid_id is not None and other.final_lipid_id is not None
                and current.target_in_fwhm(other, REDUNDANT_FWHM_DIVISOR)
                and current.final_lipid_id.sum_lipid_name
                == other.final_lipid_id.sum_lipid_name):
            # Winner by area, EXCEPT that a group already removed cannot be the winner.
            # A deliberate deviation from the original, which picks purely by area: an earlier
            # rule may already have taken the larger of the two, and merging into it then
            # removes both copies and loses the lipid outright. Ceramides are where this bites
            # here — [M+H]+ and its in-source water loss elute together, one copy goes as an
            # in-source fragment and the other as redundant.
            if current.keep and not other.keep:
                winner, loser = current, other
            elif other.keep and not current.keep:
                winner, loser = other, current
            elif current.max_area > other.max_area:
                winner, loser = current, other
            else:
                winner, loser = other, current
            winner.merge_compound_group(loser)


    def _filter_by_retention_model(self, say=None):
        """Remove identifications that sit far off their own class's retention surface.

        Distinct from `_check_class_rt_distribution`, which only asks whether a lipid falls in
        the window where its class elutes. This asks whether it falls where *that member* of the
        class should, given its carbons and double bonds — a much sharper question, and the one
        LipiDex 1 never asks.

        A class whose model does not fit is left alone entirely. A bad fit is an absence of
        evidence, not evidence of absence, and must not be allowed to remove anything.
        """
        identified = [g for g in self.compound_groups if g.keep and g.final_lipid_id is not None]

        # Fit on the well-detected rows only, then judge everything. A compound found in half
        # the samples is far more likely to be a real peak correctly assigned than a one-off,
        # and fitting on all rows is circular: the duplicated identifications that most need
        # judging are exactly the ones that spoil the fit meant to judge them. On this data PC
        # goes from R2 0.70 and unusable to R2 0.87 purely by fitting on the reliable subset.
        threshold = max(2, len(self.samples) // 2)
        reliable = [g for g in identified if len(g.compounds) >= threshold] or identified
        self.retention_models = fit_models([(g.identification()[0], g.retention)
                                            for g in reliable])
        if not self.retention_models:
            return

        usable = [m for m in self.retention_models.values() if m.usable]
        if say:
            say(f"retention models: {len(usable)} usable of {len(self.retention_models)} classes, "
                f"fitted on {len(reliable)} of {len(identified)} identifications")
            for model in sorted(usable, key=lambda m: -m.n_used)[:8]:
                say(f"    {model}")

        removed = 0
        for group in identified:
            error = retention_error(self.retention_models, group.identification()[0],
                                    group.retention)
            group.retention_error = error
            if error is None:
                continue
            model = self.retention_models[parse_sum_name(group.identification()[0])[:2]]
            if model.residual_sd:
                group.retention_z = error / model.residual_sd
            group.retention_model_r2 = model.r2
            limit = max(self.retention_model_sigma * model.residual_sd,
                        self.retention_model_max_error)
            if abs(error) > limit:
                group.keep = False
                group.filter_reason = "RT model outlier"
                removed += 1
        if say and removed:
            say(f"    {removed} identifications removed as retention-model outliers")

    def _check_class_rt_distribution(self, multiplier: float):
        """Remove identifications that fall outside where their class elutes.

        `CDPeakFinder :: checkClassRTDist`. The window is the class median plus or minus
        `multiplier * MAD / 0.6745` — the median absolute deviation rescaled to a standard
        deviation, so it is robust to the outliers it is looking for. Classes with four or
        fewer identifications are left alone, there being no distribution to speak of.

        Note what this is *not*: it is a clustering test, not an elution-order test. LipiDex
        never uses the fact that within a class retention rises with carbon number and falls
        with double bonds, which is the strongest orthogonal retention criterion available in
        reversed-phase lipidomics. Nothing here changes that; it is only reproduced faithfully.
        """
        by_class: dict[str, list[CompoundGroup]] = {}
        for group in self.compound_groups:
            if group.final_lipid_id is not None:
                by_class.setdefault(group.final_lipid_id.lipid_class, []).append(group)

        for groups in by_class.values():
            if len(groups) <= 4:
                continue
            retentions = [g.retention for g in groups]
            median = _median(retentions)
            mad = _median([abs(r - median) for r in retentions])
            window = abs(multiplier * mad / MAD_TO_SIGMA)
            if window <= 0:
                continue
            for group in groups:
                if group.keep and not (median - window <= group.retention <= median + window):
                    group.keep = False
                    group.filter_reason = "RT out of class range"


def _adduct_of(group: CompoundGroup) -> str | None:
    if group.final_lipid_id is None or not group.final_lipid_id.identifications:
        return None
    return group.final_lipid_id.identifications[0].adduct


def _fragments_of(group: CompoundGroup) -> list[float]:
    if group.final_lipid_id is None or not group.final_lipid_id.identifications:
        return []
    return group.final_lipid_id.identifications[0].fragment_masses


def _median(values: list[float]) -> float:
    if not values:
        return 0.0
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle] + ordered[middle - 1]) / 2.0


class RetentionIndex:
    """Compound groups bucketed by retention time, as `createIndexMap` does it."""

    def __init__(self, groups: list[CompoundGroup], avg_fwhm: float, bucket: float = 0.1):
        self.groups = groups
        self.bucket = bucket
        self.window = max(avg_fwhm * 3.0, 1e-6)
        self.map: dict[int, list[CompoundGroup]] = {}
        for group in groups:
            self.map.setdefault(int(group.retention / bucket), []).append(group)

    def near(self, retention: float):
        lo = int((retention - self.window) / self.bucket)
        hi = int((retention + self.window) / self.bucket)
        for b in range(lo, hi + 1):
            yield from self.map.get(b, ())


def parse_purity(text: str) -> list[tuple[str, int]]:
    """`PC 34:1 [M+H]+(82) / PC 34:1 [M+Na]+(18)` -> [(name, 82), (name, 18)]."""
    out = []
    for part in (text or "").split(" / "):
        part = part.strip()
        if not part or not part.endswith(")") or "(" not in part:
            continue
        name, _, value = part.rpartition("(")
        try:
            out.append((name.strip(), int(value.rstrip(")"))))
        except ValueError:
            continue
    return out


def parse_masses(text: str) -> list[float]:
    out = []
    for part in (text or "").split("|"):
        part = part.strip()
        if part:
            try:
                out.append(float(part))
            except ValueError:
                pass
    return out


RESULT_HEADER = ["Retention Time (min)", "Quant Ion", "Polarity", "Area (max)",
                 "Identification", "Lipid Class", "Features Found"]

# Every column `write_results` can emit before the per-injection areas. Anything in a results
# file that is not one of these is an injection. Kept here, beside the writer, because every
# reader has to agree with the writer about where the sample columns start — and a reader that
# guesses wrong does not fail, it silently reports `Purity` as a sample and the numbers still
# look plausible.
# Every non-injection column this module can write. `qc.load_run` takes everything AFTER the last
# of these to be an injection, so a new column missing from this set is silently read as a sample.
# `Pooled QC CV (%)` did exactly that: it contains "QC", so it was not merely counted as an extra
# injection but classified as a fifth pooled QC, and every precision, missingness and
# principal-component figure was computed over it. `test_every_written_column_is_registered` fails
# if a column is ever added without being listed here.
META_COLUMNS = frozenset(RESULT_HEADER) | {
    "Compound Group", "Filter Status", "Adduct", "Dot Product", "Reverse Dot Product",
    "Purity", "MS2 Spectra", "MS2 Files", "Identification Source", "Pooled QC CV (%)",
    "Quant Ion (measured)", "Duplicate Name", "RT Model Error", "Duplicate Verdict",
    "Chains From", "Lipid Key", "RT Model Z", "RT Model R2", "Chain Evidence", "Purity Source",
    "Shorthand (LSI)", "Class (canonical)", "S/N", "Theoretical m/z", "Mass Error (ppm)"}


def pooled_cv(result: PeakFinderResult, groups,
              scale: "dict[int, float] | None" = None) -> dict[int, float]:
    """Per-feature coefficient of variation across the pooled QC injections, in percent.

    Computed on the areas **exactly as written in this table**, with no normalisation, so anyone
    filtering at 20 or 30% can recompute the number from the columns beside it and get the same
    answer. That reproducibility is the point: a delivered table whose quality column cannot be
    derived from the table is a column you have to take on trust.

    It therefore differs from the precision the quality report shows, which median-normalises
    first and so removes injection-to-injection loading differences. Where loading is steady the
    two nearly agree — on the study this was written against, median 8.8% here against 7.9%
    there — and where they diverge, the difference *is* the loading variation, which is worth
    seeing rather than hiding.

    Needs at least three pools with signal for the feature. Fewer is not a precision estimate, and
    a standard deviation over two points reported to a decimal place invites a filter to be applied
    to it anyway, so those come back absent rather than optimistic.
    """
    from statistics import mean, stdev

    columns = [i for i, s in enumerate(result.samples) if s.role == "qc"]
    if len(columns) < 3:
        return {}
    out: dict[int, float] = {}
    for group in groups:
        values = [group.areas[i] / ((scale or {}).get(i) or 1.0) for i in columns
                  if i < len(group.areas) and group.areas[i] > 0]
        if len(values) < 3:
            continue
        centre = mean(values)
        if centre > 0:
            out[id(group)] = 100.0 * stdev(values) / centre
    return out


def median_factors(result: PeakFinderResult, groups, columns) -> list[float]:
    """Per-injection median-normalisation factors, scaled to average one.

    Derived only from features detected in **every** included injection. A feature missing from one
    column would otherwise contribute a zero to that column's median and drag its factor down, so
    the correction would be largest exactly where detection was worst — which is backwards.

    Returns all-ones when no feature is present everywhere, because a normalisation with nothing to
    derive itself from should leave the data alone rather than invent a scale.
    """
    from statistics import mean, median

    complete = [g for g in groups
                if all(i < len(g.areas) and g.areas[i] > 0 for i in columns)]
    if not complete:
        return [1.0] * len(columns)
    totals = [median([g.areas[i] for g in complete]) for i in columns]
    centre = mean(totals)
    if not centre:
        return [1.0] * len(columns)
    return [t / centre if t else 1.0 for t in totals]


def _lipid_keys(groups) -> dict[int, str]:
    """A name unique within the table, for plotting and for one-row-per-test analyses.

    `Identification` is left exactly as it is: it is the join key across all four tables, across
    `Associated_Spectra.csv` and across runs, and mutating it to gain uniqueness would break joins
    to buy a property that a second column supplies for nothing.

    Duplicated names take a LETTER — `PC 34:1a`, `PC 34:1b` — **ordered by intensity, largest
    first**, so `a` is the dominant species and the letters rank what matters rather than what
    happens to elute first. On `PC 38:3` the four rows span nearly three hundred fold, from
    1.2 billion down to 3.9 million; retention order would have put the smallest of those in the
    middle and said nothing. ⚠ The consequence is that a letter is only stable within a study: if
    the relative intensities shift between runs, so do the letters. Join on `Identification` plus
    retention across runs, never on the key.

    A letter rather
    than a number because the numbers in a lipid name already mean carbon count and double bonds,
    so `PC 34:1(2)` invites exactly the misreading it is meant to prevent; and not an underscore,
    because `_` separates chains, and "does the name contain an underscore" is how the pipeline
    tells a chain-resolved name from a sum composition.

    Every member is suffixed, including the first. With `PC 34:1` beside `PC 34:1b` a reader cannot
    tell whether the bare name is unique or merely the earliest of several; with `a` and `b` the
    presence of a suffix always means the name is shared.

    Unique names are returned unchanged, so the column can be used for every row without thinking.
    """
    order: dict[str, list] = {}
    for group in groups:
        name = (group.identification()[0] or "").strip()
        if name:
            order.setdefault(name, []).append(group)
    out: dict[int, str] = {}
    for name, members in order.items():
        if len(members) < 2:
            out[id(members[0])] = name
            continue
        for index, group in enumerate(sorted(members, key=lambda g: -g.max_area)):
            out[id(group)] = f"{name}{_letter(index)}"
    return out


def _letter(index: int) -> str:
    """a..z, then aa, ab, ... — so a name shared by more than twenty-six rows still resolves."""
    letters = ""
    index += 1
    while index:
        index, remainder = divmod(index - 1, 26)
        letters = chr(ord("a") + remainder) + letters
    return letters


def _duplicate_names(groups) -> dict[int, str]:
    """`n of m` per group, for identifications carried by more than one row.

    Empty when every name is unique, so a table without duplicates gains no column and no reader
    goes looking for a distinction that is not there.
    """
    order: dict[str, list] = {}
    for group in groups:
        name = (group.identification()[0] or "").strip()
        if name:
            order.setdefault(name, []).append(group)
    out: dict[int, str] = {}
    for name, members in order.items():
        if len(members) < 2:
            continue
        for index, group in enumerate(sorted(members, key=lambda g: g.retention), start=1):
            out[id(group)] = f"{index} of {len(members)}"
    return out


# A row is judged unreliable above this pooled-QC CV: the report's own band, Dunn et al. 2011, 30%
# for untargeted profiling.
#
# ⚠ Retention deliberately does NOT contribute a verdict here, though it did until Jair asked
# whether the model could tell `PC 18:0_20:3` from `PC 18:1_20:2`. It cannot — both collapse to
# ("PC", "", 38, 3), so every row of a duplicated name receives the SAME prediction and their
# retention errors are one number minus a constant. At most one isomer can sit on that prediction;
# the rest are displaced *by construction*, which is measurable: within a duplicate set the nearest
# row sits 0.099 min out, in line with unique names at 0.112, and its siblings sit at 0.391. An
# "off-model" verdict on a duplicate therefore meant "you are not the isomer closest to the class
# average" — a statement about position in an eluting series, not evidence against the row.
# `RT Model Z` and `RT Model R2` are still reported, and remain sound for unique names, where
# displacement has no such innocent explanation.
#
# What retention CAN still say about a duplicated name is a statement about the whole set rather
# than about a row within it: if even the row NEAREST the prediction is far from it, no arrangement
# of isomers explains that, and the class model or the identification is wrong. The threshold
# applies to the set's smallest |z| and is the same three sigma, which is legitimate precisely
# because that quantity was measured to behave like a unique name — nearest row 0.099 min out
# against 0.112 for names appearing once.
UNRELIABLE_CV = 30.0
SET_OFF_MODEL_Z = 3.0


def cvs_preview(result, groups, qc_cv: bool) -> dict:
    return pooled_cv(result, groups) if qc_cv else {}


def _duplicate_verdicts(groups, duplicates: dict, cvs: dict) -> dict[int, str]:
    """Why a duplicated row might not be worth keeping — or that nothing says it isn't.

    Only rows sharing a name are judged: for a unique identification there is nothing to choose
    between, and a verdict there would read as a quality score for the whole table.

    Deliberately does NOT delete. Of 165 duplicated names on the study this was written for, this
    marks a minority; the rest are several rows all precise, which is what genuine chromatographic
    isomers sharing a sum-composition label look like. Choosing between those needs a
    chromatographic argument, not a threshold.

    Two criteria, and they answer different questions. `unreliable` is per row and is pooled-QC
    precision. `class model?` is per NAME: it marks every row of a name whose rows are ALL far from
    the prediction, which is not a claim that any one of them is wrong but that the class model or
    the identification is. See `SET_OFF_MODEL_Z` for why the set can be judged when its rows
    individually cannot.
    """
    members: dict[str, list] = {}
    for group in groups:
        if id(group) in duplicates:
            members.setdefault((group.identification()[0] or "").strip(), []).append(group)

    displaced = set()
    for name, rows in members.items():
        scores = [abs(g.retention_z) for g in rows if g.retention_z is not None]
        # every row must have a score: a set judged on the two of five that carry one would be
        # judged on whichever rows the model happened to reach, which is not the whole set
        if len(scores) == len(rows) and scores and min(scores) > SET_OFF_MODEL_Z:
            displaced.add(name)

    out: dict[int, str] = {}
    for group in groups:
        if id(group) not in duplicates:
            continue
        reasons = []
        cv = cvs.get(id(group))
        if cv is not None and cv > UNRELIABLE_CV:
            reasons.append("unreliable")
        if (group.identification()[0] or "").strip() in displaced:
            reasons.append("class model?")
        out[id(group)] = "; ".join(reasons) if reasons else "isomer?"
    return out


def recentre(groups, columns, factors) -> list[float]:
    """Fold a second scaling into `factors` so per-injection medians are equal on these rows.

    Medians are taken over cells with signal: a zero is a non-detection, not a small number, and
    letting it into a median would drag a column down in proportion to how much it failed to
    detect — the opposite of a correction.
    """
    from statistics import mean, median

    scaled = []
    for position, factor in zip(columns, factors):
        values = [g.areas[position] / factor for g in groups
                  if position < len(g.areas) and g.areas[position] > 0]
        scaled.append(median(values) if values else 0.0)
    centre = mean([v for v in scaled if v > 0] or [1.0])
    return [f * (s / centre) if s > 0 else f for f, s in zip(factors, scaled)]


def write_results(result: PeakFinderResult, path: str | Path, unfiltered: bool = False,
                  id_source: bool = False, ms2_support: bool = True,
                  scores: bool = True, adduct: bool = True,
                  group_id: bool = True, identified_only: bool = False,
                  qc_cv: bool = True, drop_roles: "tuple[str, ...]" = (),
                  normalise: bool = False, lipid_key: bool = True,
                  shorthand_column: bool = True) -> None:
    """Write `Final_Results.csv`, or `Unfiltered_Results.csv` with the filter reason column.

    The optional columns all go after `Features Found` — the same place LipiDex itself inserts
    `Filter Status` — so anything locating the sample columns by name is unaffected.

    `scores` adds `Dot Product`, `Reverse Dot Product` and `Purity` for the name as reported. A
    result table without them asks the reader to trust every row equally, and the rows are not
    equal: on this data identified rows run from a dot product near 1000 down to below 500, and
    the difference between those two is the difference between a molecule and a guess. It matters
    most when nothing has been filtered out — see `unfiltered`, where a weak row and a strong one
    sit side by side.

    `ms2_support` adds `MS2 Spectra` and `MS2 Files`: how much fragmentation evidence stands
    behind the identification. `Features Found` counts something quite different — the files
    where the *MS1 feature* was detected — and the two can disagree sharply. A row can be
    quantified in seven files on the strength of a single spectrum. On by default because that
    distinction is easy to misread and expensive to misread.

    `id_source` distinguishes an MS2-derived identification from one assigned by the retention
    model. Off unless retention-model extension is used, because a row named without
    fragmentation behind it must say so.

    `adduct` adds `Adduct`: which ion the row was quantified on. See `CompoundGroup.adduct` for
    why the name alone is not enough.

    `group_id` adds `Compound Group`, a stable identifier for the row. It is the group's position
    in the unfiltered set, so the same group carries the same id in `Final_Results.csv` and in
    `Unfiltered_Results.csv`, and `Associated_Spectra.csv` can point at it. Without it the only
    join key is name plus retention time, which is precisely the thing that is ambiguous when a
    molecule appears on more than one row.

    `lipid_key` adds `Lipid Key`, a name unique within the table for plotting and for
    one-row-per-test analyses. Off gives exactly the columns LipiDex writes.

    `identified_only` drops rows carrying no name. Used for `Final_Results_Filtered.csv`, the
    analysis-ready table: an unnamed feature cannot enter a lipid-level statistical analysis, and
    leaving thousands of them in a file labelled analysis-ready invites someone to model them by
    accident. They are not lost — `Final_Results.csv` keeps every kept row, named or not, and
    `Unfiltered_Results.csv` keeps the rejected ones with the reason.

    Set all of them false for exactly the columns LipiDex writes.
    """
    header = list(RESULT_HEADER)
    if shorthand_column:
        # Beside the name, not instead of it. The pipeline's own name carries what the shorthand
        # cannot — `Cer[NS]` names the LipiDex subclass, `Plasmanyl-PC` says which ether matched —
        # so both are written. See `lsi.py` for the claims this column must never over-state.
        # Optional, because "every column off gives exactly what LipiDex writes" is a promise the
        # test suite holds this writer to, and a new column on by default would break it.
        header.insert(header.index("Identification") + 1, "Shorthand (LSI)")
    if shorthand_column:
        # ⚠ APPENDED, not inserted. The insertions above are index-based and the comment on the
        # row side records what happened last time two of them were ordered wrongly: every value
        # present, every one under the wrong heading. Appending is order-independent.
        #
        # `Lipid Class` keeps whatever the matching library called it, because that is provenance —
        # it says which library produced the hit. This column says which class it IS, with the
        # spelling differences resolved: our own two libraries write the same molecules as
        # `GlcCer[NS]` and `HexCer[NS]`, and grouping on the raw name splits one class in two.
        header.append("Class (canonical)")
    if group_id:
        header.insert(0, "Compound Group")
    if unfiltered:
        header.append("Filter Status")
    if adduct:
        header.append("Adduct")
    if scores:
        header += ["Dot Product", "Reverse Dot Product", "Purity", "Chain Evidence",
                   "Purity Source"]
    if ms2_support:
        header += ["MS2 Spectra", "MS2 Files"]
    if id_source:
        header.append("Identification Source")
    # Written only when a correction was actually applied — an extra column identical to the one
    # beside it invites the reader to look for a difference that is not there.
    measured = any(g.quant_ion_measured for g in result.compound_groups)
    if measured:
        header.append("Quant Ion (measured)")
    # ⚠ WHICH ROWS are written is decided BEFORE the normalisation factors, because the factors
    # are a median across rows and must come from the rows they will be applied to. Deriving them
    # from every compound group — thirteen thousand of them, mostly rejected features and noise —
    # and applying them to eight hundred identified lipids means correcting one population by a
    # factor measured on another.
    #
    # `unfiltered` decides both what is written and whether the reason column appears, because a
    # table containing dropped rows without saying which they are would be a trap.
    groups = result.compound_groups if unfiltered else result.kept
    if identified_only:
        groups = [g for g in groups if (g.identification()[0] or "").strip()]

    # A row is a compound group, not a lipid: 44% of positive rows on the study this was added for
    # shared a name with another row. Three causes need three different answers — a chromatographic
    # isomer, a second adduct, a split peak — so nothing is merged here. The count makes the
    # duplication visible, because the file is called analysis-ready and someone reading that will
    # run one test per row. Counted over the rows WRITTEN: "1 of 4" has to mean four rows in this
    # file, not four somewhere upstream.
    duplicates = _duplicate_names(groups)
    keys = _lipid_keys(groups) if lipid_key else {}
    if keys:
        header.insert(header.index("Identification") + 1, "Lipid Key")
    if duplicates:
        header.append("Duplicate Name")
    # How far a row sits from where its class says its composition should elute. Computed by the
    # peak finder and, until now, thrown away, though it is the only per-row evidence independent
    # of intensity.
    model_errors = any(g.retention_error is not None for g in groups)
    if model_errors:
        # Three numbers, because one does not carry it. The error says how far off in minutes; the
        # z says how far off *for this class*, which is the only version comparable between
        # classes; the R2 says whether the class model is worth believing at all — a z of 3 against
        # an R2 of 0.19 means little, and PS fits at 0.19 here.
        #
        # ⚠ Reported, not acted on, and read them only for rows whose name appears once. The model
        # keys on sum composition, so every row of a duplicated name shares one prediction — see
        # the note above `UNRELIABLE_CV`.
        header += ["RT Model Error", "RT Model Z", "RT Model R2"]
    verdicts = _duplicate_verdicts(groups, duplicates, cvs_preview(result, groups, qc_cv))
    if verdicts:
        header.append("Duplicate Verdict")

    # Columns to omit. The analysis-ready table drops the standards: they are a different material
    # injected to judge the instrument, and leaving them beside the samples invites them into a
    # fold change. They keep their own table, Standards_Performance.csv.
    keep = [i for i, s in enumerate(result.samples) if s.role not in drop_roles]
    factors = (median_factors(result, groups, keep)
               if normalise else [1.0] * len(keep))
    if normalise:
        # ⚠ Re-centre on the rows actually written. The factors above come from features detected
        # in EVERY injection, which is the stable basis — it cannot be pulled about by which
        # features happened to be missing where — but it is not the whole table, so the written
        # matrix is left a few per cent off centre. On this study that residual was 12% in positive.
        #
        # Two steps rather than one because they answer different needs: derive on a set that does
        # not move when the library changes, then make the delivered file self-consistent, so
        # "every sample has the same median on the rows in this file" is true by construction and
        # can be asserted. A filter applied later — polarity, say — changes which rows are written
        # and would otherwise silently decentre the matrix again, by a different route.
        factors = recentre(groups, keep, factors)
        scale = dict(zip(keep, factors))
    else:
        scale = dict(zip(keep, factors))
    # The CV column is always computed from the areas in the *same file*, so it can be recomputed
    # from the columns beside it whichever table the reader opens.
    cvs = pooled_cv(result, groups, scale if normalise else None) if qc_cv else {}
    if cvs:
        header.append("Pooled QC CV (%)")
    # ⚠ LAST, and gated on the same flag as the row side. Height over a robust local baseline sigma for the quantified ion,
    # in the sample where it is strongest — so a reader can tell an area that is a measurement
    # from one that is an integration of noise, which the abundance column cannot show.
    #
    # Appended at the very end rather than beside `Class (canonical)`: the conditional blocks above
    # do not fire in the same combinations on the header and row sides, and a first attempt to put
    # it there shifted every score column one place, so `Dot Product` read out an adduct.
    if shorthand_column:
        header.append("S/N")
    # ⚠ After S/N, for the same reason S/N went after `Class (canonical)`: the conditional blocks
    # above do not fire in the same combinations on the header and row sides, so anything inserted
    # mid-table shifts every score column and `Dot Product` starts reading out an adduct.
    if shorthand_column:
        header += ["Theoretical m/z", "Mass Error (ppm)"]
    header += [result.samples[i].file for i in keep]
    # Numbered over the full set, never over what is written, so `Compound Group` means the same
    # thing in every table and a row can be joined across them.
    ids = {id(g): n for n, g in enumerate(result.compound_groups, start=1)}
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header + [""])
        for group in groups:
            name, lipid_class = group.identification()
            row = [group.retention, group.quant_ion, group.quant_polarity, group.max_area,
                   name, lipid_class, len(group.compounds)]
            # `Identification Source` decides whether a plasmalogen claim is allowed: an
            # `RT model` row has no MS2, so `P-` is demoted to `O-` rather than asserted.
            # ⚠ Insert in the SAME ORDER the header does, or the columns shift silently. The
            # header adds `Shorthand (LSI)` at Identification+1 and then `Lipid Key` at
            # Identification+1 as well, which puts the key BEFORE the shorthand. Inserting the row
            # values the other way round produced a table whose `Lipid Key` column held shorthand
            # and vice versa — every value present, every one under the wrong heading.
            if shorthand_column:
                row.insert(RESULT_HEADER.index("Identification") + 1,
                           shorthand(name, group.identification_source()))
            if keys:
                row.insert(RESULT_HEADER.index("Identification") + 1, keys.get(id(group), ""))
            if group_id:
                row.insert(0, ids.get(id(group), ""))
            if shorthand_column:
                row.append(canonical_class(lipid_class or name))
            if unfiltered:
                row.append(group.filter_reason)
            if adduct:
                row.append(group.adduct())
            if scores:
                row += (list(group.identification_scores())
                        + [group.chain_evidence(), group.purity_source()])
            if ms2_support:
                row += list(group.ms2_support())
            if id_source:
                row.append(group.identification_source())
            if measured:
                row.append(group.quant_ion_measured if group.quant_ion_measured else "")
            if duplicates:
                row.append(duplicates.get(id(group), ""))
            if model_errors:
                row.append(f"{group.retention_error:+.3f}"
                           if group.retention_error is not None else "")
                row.append(f"{group.retention_z:+.2f}" if group.retention_z is not None else "")
                row.append(f"{group.retention_model_r2:.3f}"
                           if group.retention_model_r2 is not None else "")
            if verdicts:
                row.append(verdicts.get(id(group), ""))
            if cvs:
                value = cvs.get(id(group))
                row.append(f"{value:.1f}" if value is not None else "")
            if shorthand_column:
                snrs = [f.snr for c in group.compounds for f in c.features
                        if f.snr is not None and group.quant_ion is not None
                        and abs(f.mass - group.quant_ion) < 0.01]
                row.append(f"{max(snrs):.0f}" if snrs else "")
            # ⚠ AFTER the S/N block, mirroring the header exactly. Placed before it — which is
            # where this first went — the row carries theoretical mass where the header says
            # `Pooled QC CV`, and every column from there rightward reads the wrong heading.
            if shorthand_column:
                theoretical = group.theoretical_mz()
                error = group.mass_error_ppm()
                row.append(f"{theoretical:.4f}" if theoretical is not None else "")
                row.append(f"{error:+.2f}" if error is not None else "")
            areas = [group.areas[i] / (scale.get(i) or 1.0) if i < len(group.areas) else ""
                     for i in keep]
            writer.writerow(row + areas + [""])
