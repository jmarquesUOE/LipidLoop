"""The peak finder: join MS2 identifications to quantified features and write `Final_Results.csv`.

`compound_discoverer/CDPeakFinder.java`, plus `CDCompoundGroup`, `CDCompound`, `CDFeature`,
`peak_finder/Sample`, `Lipid` and `GaussianModel`.

The identification stage answers "what is this spectrum". This stage answers the question the
result file actually reports: *which chromatographic peak* is that lipid, and how much of it is
there in each sample. The two are joined on mass, polarity, sample and retention time, and the
retention-time part is the hard one, because an MS2 is taken at whatever moment the instrument
decided to fragment while the feature's retention time has been aligned across files. LipiDex
recovers the offset by fitting the deviation between each feature's aligned and measured
retention time, per sample, and mapping the MS2 back through it.

Three layers, mirroring what Compound Discoverer exports:

    compound group   one molecule, aligned across samples   -> one output row
      compound       that molecule in one sample
        feature      one adduct ion of it in that sample     <- IDs attach here

An identification attaches to a *feature*, is promoted to the compound group by majority vote
over sum compositions, and is reported at molecular resolution only if the purity clears
`MINFAPURITY`; otherwise the group is named at sum composition.
"""
from __future__ import annotations

import csv
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

from .loess import loess

MAX_PPM_DIFF = 20.0        # CDPeakFinder.MINPPMDIFF
# `Utilities.MINRTMULTIPLIER`, set from the GUI's FWHM spinner, default 2.0. It governs two
# different things: the FWHM multiple inside which an MS2 may be attributed to a feature, and
# the width of the per-class retention window. Not to be confused with CDPeakFinder's *instance*
# field of the same name, set from a separate spinner and passed to checkClassRTDist as an
# argument that method never reads — dead in the original, so it is not reproduced.
MIN_RT_MULTIPLIER = 2.0
MIN_DOT_PRODUCT = 500.0
MIN_REV_DOT_PRODUCT = 700.0
MIN_FA_PURITY = 75.0
MIN_ID_NUM = 1             # compounds per group required to keep it
UNALIGNED_PPM = 5.0        # matchUnalignedFeature
ISOTOPE_SPACING = 1.003


def ppm_diff(mass1: float, mass2: float) -> float:
    return abs(mass1 - mass2) / mass2 * 1e6


_ADDUCT_SUFFIX = re.compile(r"\s*\[M[^\]]*\][+-]?.*$")


def _chain_key(name: str) -> tuple:
    """`Cer[ADS] d17:0_17:0 [M+FA-H]-` -> ("Cer[ADS]", ("17:0", "d17:0")).

    Identity of a chain-resolved lipid for the purpose of asking whether a purity breakdown is
    about *this* row. Chains are sorted because component names permute their order —
    `PC[OH] OH-16:0_18:2` and `OH-18:2_16:0` are one lipid, and comparing the strings says they
    are two. The class is part of the key: `Cer[ADS]` and `Cer[NP]` at the same chains are
    different molecules, and telling them apart is the whole point.
    """
    text = _ADDUCT_SUFFIX.sub("", (name or "").strip()).rstrip(";").strip()
    lipid_class, _, chains = text.partition(" ")
    return (lipid_class, tuple(sorted(chains.split("_"))))


def _class_of(name: str, fallback: str) -> str:
    """The class token of a name, so a demoted name cannot disagree with its own class column."""
    token = (name or "").strip().split(" ", 1)[0]
    return token or fallback


def _demote_proven_ether(name: str) -> str:
    """`Alkenyl-TG P-52:3` -> `TG O-52:3`. `P-` claims a *proven* alk-1-enyl bond.

    Shaped like `top_ether`: the class prefix goes with the claim it encoded, because `Alkenyl-TG`
    means alk-1-enyl just as surely as `P-` does, and keeping it while writing `O-` would rebuild
    the same contradiction one field along. Anything already at `O-` is untouched — that is the
    weaker claim and is what remains when nothing is proven.
    """
    text = (name or "").strip()
    lipid_class, _, rest = text.partition(" ")
    if not rest.startswith("P-"):
        return text
    return f"{lipid_class.split('-')[-1]} O-{rest[2:]}"


def _sum_key(name: str):
    """`PC 16:0_18:1` -> ("PC", 34, 1). None if the chains cannot be read."""
    lipid_class, chains = _chain_key(name)
    carbons = doubles = 0
    for chain in chains:
        found = re.search(r"(\d+):(\d+)", chain)
        if not found:
            return None
        carbons += int(found.group(1))
        doubles += int(found.group(2))
    return (lipid_class, carbons, doubles)


def _round_half_up(x: float) -> int:
    return int(math.floor(x + 0.5))


# ─────────────────────────────────────────── peak shape

class GaussianModel:
    """A Gaussian standing in for the real chromatographic peak.

    Compound Discoverer reports an area and a width but not a trace, so LipiDex rebuilds the
    peak as a Gaussian of that FWHM, sampled every FWHM/12 out to ±3 FWHM, and integrates it by
    trapezium to recover a height. The point of it is `normalized_height`: how far up the peak
    an MS2 was taken. A spectrum fragmented on the apex is better evidence for the feature's
    identity than one fragmented on the tail, and that weight carries into the purity.
    """

    __slots__ = ("fwhm", "area", "apex_rt", "x", "y", "height")

    def __init__(self, fwhm: float, area: float | None, apex_rt: float,
                 height: float | None = None):
        self.fwhm = fwhm
        self.area = area
        self.apex_rt = apex_rt
        self.x: list[float] = []
        self.y: list[float] = []

        sigma = fwhm / 2.3548 if fwhm else 1e-9
        time = -fwhm * 3.0
        step = fwhm / 12.0 if fwhm else 1.0
        while time < fwhm * 3.0:
            self.x.append(time)
            # 2.71828 rather than e: the original hard-codes it, and the height it produces
            # feeds a weighted average, so the difference is not quite nothing.
            self.y.append(2.71828 ** (-(time ** 2 / (2.0 * sigma ** 2))))
            time += step

        trapezoid = sum((self.y[i] + self.y[i + 1]) / 2.0 * (self.x[i + 1] - self.x[i])
                        for i in range(len(self.x) - 1))
        if height is None:
            self.height = (area / trapezoid) if trapezoid else 0.0
        else:
            self.height = height

    def normalized_height(self, time: float) -> float:
        """Height at `time` on a 0-1 scale, or 0.001 outside the modelled window."""
        normalized = time - self.apex_rt
        for i in range(len(self.x) - 1):
            if self.x[i] < normalized < self.x[i + 1]:
                if abs(normalized - self.x[i]) > abs(normalized - self.x[i + 1]):
                    return self.y[i + 1]
                return self.y[i]
        return 0.001


# ─────────────────────────────────────────── data model

@dataclass(slots=True)
class Sample:
    file: str
    cd_id: str = ""
    role: str = "sample"        # sample | blank | qc — see blanks.py
    injection: int | None = None
    rt_deviations: list[tuple[float, float, float]] = field(default_factory=list)
    _rt: list[float] = field(default_factory=list)
    _smoothed: list[float] = field(default_factory=list)
    ppm_errors: list[float] = field(default_factory=list)
    purities: list[int] = field(default_factory=list)
    fwhms: list[float] = field(default_factory=list)
    compounds_found: int = 0
    lipids: int = 0
    max_feature_rt: float = 0.0

    def fit(self) -> None:
        """LOESS the (aligned RT, deviation) pairs — `Sample.java :: loessFit`."""
        if not self.rt_deviations:
            return
        self.rt_deviations.sort(key=lambda d: d[1])
        # The interpolator needs strictly increasing x.
        retentions = []
        last = None
        for _, retention, _ in self.rt_deviations:
            if last is not None and retention <= last:
                retention = last + 0.00001
            retentions.append(retention)
            last = retention
        deviations = [d[0] for d in self.rt_deviations]
        self._rt = retentions
        self._smoothed = list(loess(retentions, deviations))

    def corrected_rt(self, rt: float) -> float:
        """Map a measured retention time into the aligned frame.

        Returns 0.0 past the end of the fitted range, which is what the original does — and
        means a late MS2 simply fails to associate rather than associating wrongly.
        """
        for i, retention in enumerate(self._rt):
            if rt - retention < 0:
                return rt - self._smoothed[i]
        return 0.0


@dataclass(slots=True)
class Feature:
    """One adduct ion in one sample. Identifications attach here."""
    adduct: str
    charge: int
    mw: float
    mass: float          # m/z
    retention: float     # aligned
    fwhm: float
    mi: int
    area: float
    parent_area_percent: float
    sample: Sample
    polarity: str = ""
    #: Height over a robust local baseline sigma; None where no leading baseline exists.
    snr: float | None = None
    real_retention: float = 0.0
    peak_model: GaussianModel | None = None
    lipids: list["Lipid"] = field(default_factory=list)
    isobaric_neighbours: list["Feature"] = field(default_factory=list)

    def __post_init__(self):
        self.polarity = "+" if self.charge > 0 else "-"
        self.peak_model = GaussianModel(self.fwhm, self.area, self.retention)
        self.sample.max_feature_rt = max(self.sample.max_feature_rt, self.retention)

    def most_intense_at(self, rt: float) -> bool:
        """True if no isobaric neighbour is taller here — `CDFeature :: isMostIntensePeak`."""
        mine = self.peak_model.normalized_height(rt) * self.peak_model.height
        return all(n.peak_model.normalized_height(rt) * n.peak_model.height <= mine
                   for n in self.isobaric_neighbours)

    def accepts(self, lipid: "Lipid", no_coeluting: bool) -> bool:
        """`CDFeature :: checkLipid` — polarity, mass, sample, peak dominance, retention."""
        if lipid.polarity != self.polarity:
            return False
        if not ppm_diff(lipid.precursor, self.mass) < MAX_PPM_DIFF:
            return False
        if self.sample.file != lipid.sample.file:
            return False
        if not no_coeluting and not self.most_intense_at(self.sample.corrected_rt(lipid.retention)):
            return False
        window = self.fwhm * MIN_RT_MULTIPLIER
        return self.real_retention - window <= lipid.retention <= self.real_retention + window


@dataclass(slots=True)
class Compound:
    mw: float
    retention: float
    fwhm: float
    max_mi: int
    n_adducts: int
    area: float
    sample: Sample
    features: list[Feature] = field(default_factory=list)


@dataclass(slots=True)
class Lipid:
    """One row of a `_Results.csv`, parsed into the identification the peak finder joins on."""
    retention: float
    precursor: float
    sample: Sample
    dot: float
    rev_dot: float
    lipid_string: str
    lib_precursor: float
    purity: int
    is_lipidex: bool
    purity_array: list[tuple[str, int]]
    fragment_masses: list[float]
    corrected_retention: float = 0.0
    lipid_class: str = ""
    adduct: str = ""
    polarity: str = ""
    lipid_name: str = ""
    sum_lipid_name: str = ""
    ppm_error: float = 0.0
    keep: bool = True
    gaussian_score: float = 0.0
    preferred_polarity: bool = False
    # Where this identification came from in the search output. Carried so `Associated_Spectra`
    # can name the compound group a spectrum ended up in — without it the evidence table can
    # only be joined back on name and retention time, which is the ambiguity it exists to settle.
    ms2_id: int = 0
    rank: int = 0

    def __post_init__(self):
        self.ppm_error = (self.precursor - self.lib_precursor) / self.lib_precursor * 1e6
        self.corrected_retention = self.sample.corrected_rt(self.retention)
        self.sample.ppm_errors.append(self.ppm_error)
        self._parse()

    def _parse(self):
        name = self.lipid_string.strip().rstrip(";").strip()
        match = re.search(r"(\[M[^\]]*\][+-]\d?)", name)
        self.adduct = match.group(1) if match else ""
        body = name[:match.start()].strip() if match else name
        self.lipid_class = body.split(" ", 1)[0] if " " in body else body
        # `Lipid.parseName` upper-cases the first letter of the class, so the library's `d5TG`
        # is reported as `D5TG`. The rebuilt name is class + chains, so this propagates.
        if self.lipid_class:
            self.lipid_class = self.lipid_class[0].upper() + self.lipid_class[1:]
            if " " in body:
                body = self.lipid_class + " " + body.split(" ", 1)[1]
            else:
                body = self.lipid_class
        # `Lipid.java :: toString` drops the adduct: the reported identity is the molecule, and
        # the same molecule seen as [M+H]+ and [M+Na]+ must collapse to one name.
        self.lipid_name = body
        self.polarity = "+" if "]+" in self.adduct else "-"
        self.sum_lipid_name = sum_composition(body)


def sum_composition(body: str) -> str:
    """`PC 16:0_18:1` -> `PC 34:1`. Chains that do not parse are left alone.

    ⚠ **`fullmatch`, not `match`.** `LipiDex2_FAHFA` names its two chains `18:1-(O-18:0)` —
    hyphen-and-parenthesis, not this parser's `_`. `chains.split("_")` never splits it, so the
    whole string reaches the regex as one "chain"; `match` matched the leading `18:1` and silently
    dropped everything after it, returning `FAHFA 18:1` for a 563 Da, two-chain, chain-unresolved
    molecule — a name identical to a genuine single free fatty acid's sum composition. On
    Skin_QEplus this displaced five real `FA` identifications from the delivered table under a
    name that looked like confirmation of the thing it had actually overwritten. `fullmatch`
    requires the whole part to be consumed, so an unrecognised delimiter now falls through to the
    documented behaviour — the name is left alone — instead of being read as a legitimate chain.
    """
    if " " not in body:
        return body
    lipid_class, chains = body.split(" ", 1)
    parts = chains.split("_")
    carbons = doubles = 0
    prefixes = []
    for part in parts:
        # Chain prefixes are letters, optionally hyphenated: d18:1, P-18:0, O-16:0, OH-14:1.
        m = re.fullmatch(r"([A-Za-z]*-?)(\d+):(\d+)", part.strip())
        if not m:
            return body
        if m.group(1):
            prefixes.append(m.group(1))
        carbons += int(m.group(2))
        doubles += int(m.group(3))
    prefix = prefixes[0] if prefixes else ""
    return f"{lipid_class} {prefix}{carbons}:{doubles}"


@dataclass(slots=True)
class CompoundGroup:
    name: str
    mw: float
    retention: float
    max_area: float
    has_ms2: bool
    areas: list[float]
    compounds: list[Compound] = field(default_factory=list)
    quant_ion: float | None = None
    # The mass as the instrument reported it, kept when a systematic correction is applied to
    # `quant_ion`. A delivered table that silently carries a derived mass, with no way back to the
    # measurement, cannot be checked by whoever receives it.
    quant_ion_measured: float | None = None
    quant_polarity: str = ""
    quant_ions: list[float] = field(default_factory=list)
    avg_fwhm: float = 0.0
    max_isotope: int = 0
    sum_ids: list[str] = field(default_factory=list)
    lipid_candidates: list["LipidCandidate"] = field(default_factory=list)
    sum_id: str = ""
    purity: float = 0.0
    final_lipid_id: "LipidCandidate | None" = None
    keep: bool = True
    filter_reason: str = ""
    positive_feature: bool = False
    negative_feature: bool = False
    no_coeluting_peaks: bool = True
    summed_purities: list[tuple[str, float]] = field(default_factory=list)
    # The denominator `summed_purities` is normalised by. Kept because without it the breakdown is
    # a set of unnormalised weighted sums that only the top entry has a scale for: `purity` is
    # top / this, so every other entry is unreadable unless this is retained.
    purity_weight: float = 0.0
    plasmenyl_ether_conflict: bool = False
    # Assigned from the per-class retention model, with no MS2 behind it. Deliberately a
    # separate field from the identification: it is weaker evidence and stays labelled as such.
    rtls_identification: str = ""
    # The ion the retention-model match was made on. Determined by the match — one m/z,
    # one adduct — and previously discarded, leaving those rows outside adduct-pair removal.
    rtls_adduct: str = ""
    retention_error: float | None = None
    # The same error in units of its own class's scatter, and how well that class fits at all.
    # Minutes are not comparable between classes — half a minute is nothing in TG (residual sd
    # 0.36) and enormous in LysoPC (0.05) — so a threshold in minutes is a different test in every
    # class. The z is what the outlier filter already thresholds on.
    retention_z: float | None = None
    retention_model_r2: float | None = None

    # ── setup

    def calc_fwhm(self):
        if self.compounds:
            self.avg_fwhm = sum(c.fwhm for c in self.compounds) / len(self.compounds)

    def find_quant_ion(self):
        """The most abundant distinct ion in the group is what gets quantified."""
        best = 0.0
        for compound in self.compounds:
            for feature in compound.features:
                if self.is_unique_quant_ion(feature.mass):
                    self.quant_ions.append(feature.mass)
                if feature.area > best:
                    best = feature.area
                    self.quant_ion = feature.mass
                    self.quant_polarity = feature.polarity
                if feature.polarity == "+":
                    self.positive_feature = True
                else:
                    self.negative_feature = True

    def is_unique_quant_ion(self, mass: float) -> bool:
        return all(ppm_diff(mass, q) >= MAX_PPM_DIFF for q in self.quant_ions)

    def add_sum_id(self, sum_id: str):
        self.sum_ids.append(sum_id)

    def target_in_fwhm(self, other: "CompoundGroup", divisor: float) -> bool:
        """Within `avg_fwhm / divisor` of another group — `CDCompoundGroup.targetInFWHM`.

        The argument DIVIDES the peak width, it does not multiply it. So the pairwise filters,
        called with 0.5, look two peak widths either side, while the outer walk, called with
        2.0, looks only half a peak width. Reading it as a multiplier inverts every window in
        the sweep and makes the filters four times too narrow.
        """
        return abs(self.retention - other.retention) < self.avg_fwhm / divisor

    # ── identification

    def filter_sum_ids(self):
        """Majority vote over sum compositions; everything else is dropped."""
        if not self.sum_ids:
            return
        counts: dict[str, int] = {}
        for sum_id in self.sum_ids:
            counts[sum_id] = counts.get(sum_id, 0) + 1
        most_common = max(counts, key=lambda k: counts[k])
        self.sum_id = most_common

        by_name: dict[str, LipidCandidate] = {}
        for compound in self.compounds:
            for feature in compound.features:
                for lipid in feature.lipids:
                    if lipid.sum_lipid_name != most_common:
                        lipid.keep = False
                        continue
                    candidate = by_name.get(lipid.lipid_name)
                    if candidate is None:
                        by_name[lipid.lipid_name] = LipidCandidate(lipid)
                    else:
                        candidate.add(lipid)
        self.lipid_candidates = sorted(by_name.values())

    def calculate_weighted_purity(self):
        """Weighted purity of the group — `CDCompoundGroup.calculateWeightedPurity`.

        Each identification carries a per-name purity breakdown (the `Spectral Components`
        column) and a Gaussian score saying how far up the peak that spectrum was taken. Those
        breakdowns are summed per name across every spectrum, weighted by the score, and the
        purity of the group is

            highest per-name weighted sum / total weight

        Two details decide the number, and getting either wrong shifts identifications between
        molecular and sum-composition reporting:

        * It is the **top name's** sum over the total weight, not the mean of the top hits'
          purities. `puritySum` is accumulated in the original and then never read — dead code
          that a reimplementation naturally mistakes for the answer.
        * The total weight counts the Gaussian score **once per entry in the breakdown**, not
          once per spectrum. A spectrum whose signal splits across three candidates therefore
          contributes three times as much denominator as one that is unambiguous, which is what
          makes a genuinely mixed peak fall below the threshold.
        """
        weight_sum = 0.0
        summed: dict[str, int] = {}
        for compound in self.compounds:
            for feature in compound.features:
                for lipid in feature.lipids:
                    if not (lipid.keep and lipid.preferred_polarity):
                        continue
                    for name, value in lipid.purity_array:
                        weight_sum += lipid.gaussian_score
                        # Each contribution is rounded before it is added, not after.
                        summed[name] = summed.get(name, 0) + _round_half_up(
                            value * lipid.gaussian_score)

        self.summed_purities = sorted(summed.items(), key=lambda p: -p[1])
        self.purity_weight = weight_sum
        if self.summed_purities and weight_sum > 0:
            self.purity = self.summed_purities[0][1] / weight_sum
        else:
            self.purity = 0.0
        self._check_plasmenyl_ether()

    def _check_plasmenyl_ether(self):
        """Plasmenyl and plasmanyl forms both present, neither convincing.

        The two differ by a vinyl ether versus an alkyl ether at sn-1, and positive-mode spectra
        frequently cannot tell them apart. Where the evidence splits between them, LipiDex
        declines to choose and reports the neutral `O-` ether notation instead of asserting a
        plasmalogen. Worth knowing before reading a plasmalogen result off this column.
        """
        if not self.summed_purities:
            return
        # Compared against the raw weighted sums, not the normalised purity — as in the
        # original. Verified: `CDCompoundGroup.java:477-508` does the same, testing
        # `summedPurities.get(0).purity` against MINFAPURITY and then raw `pSum`/`eSum`
        # against the same constant. Two quantities, one threshold, deliberately.
        top_name, top_purity = self.summed_purities[0]
        if top_purity > MIN_FA_PURITY:
            return
        if "Plasmenyl" not in top_name and "Plasmanyl" not in top_name:
            return
        plasmenyl = sum(v for n, v in self.summed_purities if "Plasmenyl" in n)
        plasmanyl = sum(v for n, v in self.summed_purities if "Plasmanyl" in n)
        if 0.0 < plasmenyl < MIN_FA_PURITY and 0.0 < plasmanyl < MIN_FA_PURITY:
            self.plasmenyl_ether_conflict = True

    def top_ether(self) -> str:
        """`Plasmenyl-PE P-40:7` -> `PE O-40:7` — MzCompoundGroup.java :: getTopEther."""
        try:
            return (self.sum_id[self.sum_id.index("-") + 1:self.sum_id.index(" ")]
                    + " O" + self.sum_id[self.sum_id.rindex("-"):])
        except ValueError:
            return self.sum_id

    def best_lipid_id(self):
        if self.lipid_candidates:
            self.final_lipid_id = self.lipid_candidates[0]

    def merge_compound_group(self, other: "CompoundGroup") -> None:
        """Absorb another group reporting the same lipid — `CDCompoundGroup.mergeCompoundGroup`.

        The loser's compounds move across and the identification is recomputed over the union,
        so the surviving row is supported by every spectrum that supported either. It is a
        merge rather than a deletion: the evidence is pooled, not discarded.
        """
        self.compounds.extend(other.compounds)
        self.filter_sum_ids()
        self.calculate_weighted_purity()
        self.best_lipid_id()
        other.keep = False
        other.filter_reason = "Redundant Identification"

    # ── output

    def ms2_support(self) -> tuple[int, int]:
        """Spectra and distinct files behind the identification as reported.

        Deliberately counts the support for the name actually written, not for the best
        candidate. A row reported at sum composition is supported by every spectrum that voted
        for that sum composition, however its chains were resolved; a row reported at molecular
        resolution is supported only by the spectra backing that molecule. Counting the latter
        under a sum-composition heading would understate it, which is the mistake that made
        `GM3-NANA d34:1` look like it had no spectra at all when it has ten.

        Returns (0, 0) for a row named by the retention model — that is the point of it having
        no MS2, and `identification_source` says so.
        """
        name, _ = self.identification()
        if not name or not self.lipid_candidates:
            return 0, 0
        if name == self.sum_id or self.plasmenyl_ether_conflict:
            lipids = [l for c in self.lipid_candidates for l in c.identifications]
        else:
            lipids = list(self.lipid_candidates[0].identifications)
        return len(lipids), len({l.sample.file for l in lipids})

    def identification_scores(self) -> tuple:
        """Best dot product, reverse dot product and purity behind the name as reported.

        Scored against the name actually written, on the same principle as `ms2_support`: a row
        reported at sum composition takes the best score among every candidate that collapses to
        that sum, because they all voted for it; a row reported at molecular resolution takes the
        top candidate's. Taking the top candidate's score for a sum-composition row would credit
        the row with a molecular assignment it deliberately did not make.

        Blank for a row named by the retention model — there is no spectrum, and
        `identification_source` says so.
        """
        name, _ = self.identification()
        if not name or not self.lipid_candidates:
            return "", "", ""
        if name == self.sum_id or self.plasmenyl_ether_conflict:
            candidates = self.lipid_candidates
        else:
            candidates = self.lipid_candidates[:1]
        return (round(max(c.max_dot for c in candidates)),
                round(max(c.max_rev_dot for c in candidates)),
                max(c.purity for c in candidates))

    def chain_evidence(self) -> str:
        """Whether the chains in this row's name were measured here, or inherited from a library.

        **A chain-resolved name is not necessarily chain evidence.** `Purity` scores the winning
        identification; it does not gate the name that identification carries. The name comes from
        whichever library entry won the dot product, and library entries are routinely written at
        chain resolution — so a row can read `Cer[ADS] d18:0_24:1` on no chain evidence at all.

        ⚠⚠ **A high purity does not mean the purity belongs to this row.** `calc_purity` drops a top
        hit that scores nothing and then returns `purities[0]`, which is by then a *different*
        candidate's score, and `search.py` writes it onto the winner's row. The clearest case in the
        delivered data is `Cer[ADS] d17:0_17:0` at **Purity 100** whose entire breakdown reads
        `Cer[NP] t18:0_16:0 (100)`. The phyto ceramide was measured; the dihydro one was named.

        So the test is **whether any component names THIS row's lipid**, not whether a component
        exists. An earlier version of this method keyed on `purity_weight` being non-zero, which is
        group-wide and non-zero for every chain-resolved row in that study — it would have returned
        `fragments` for all 166 and never once said `library name`. It was caught in review before
        it ever ran, which is the only reason it is not in a delivered table.

        Chains are compared as a multiset, so `PC 16:0_18:1` and `PC 18:1_16:0` are one lipid.
        ⚠ `PC[OH] OH-16:0_18:2` and `OH-18:2_16:0` are NOT — the prefix marks which chain carries
        the hydroxyl, both exist as separate library entries at the same precursor (832.5704), and
        the sorted key distinguishes them correctly. An earlier version of this docstring gave that
        pair as an example of permutation, contradicting the test beside it.

        Measured on one negative-mode study, chain-resolved rows: **139 `fragments`, 27
        `library name`** (Cer[ADS] 15, Cer[NDS] 12) out of 166. Positive mode: 48 and 0. The
        `library name` rows are **not weak matches** — the lowest dot product among them is 699
        against 511 for rows that did score — so no threshold would catch them and raising one would
        delete good identifications instead.

        Blank where the name is a sum composition: no chains are asserted, so there is nothing to
        qualify, and a value would imply a claim the row declined to make.
        """
        name, _ = self.identification()
        if "_" not in (name or ""):
            return ""
        if not self.summed_purities:
            # No breakdown at all: off-polarity, or nothing scored. Distinguished from a breakdown
            # that names someone else, because they are different failures with different fixes.
            return "not eligible"
        mine = _chain_key(name)
        if any(_chain_key(component) == mine for component, _ in self.summed_purities):
            return "fragments"
        return "library name"

    def identification_source(self) -> str:
        if self.lipid_candidates:
            return "MS2"
        return "RT model" if self.rtls_identification else ""

    def adduct(self) -> str:
        """The ion or ions behind the name as reported, the quantified one first.

        `Lipid.java :: toString` drops the adduct, because the reported identity is the molecule
        and not the ion it happened to be seen as. That is right for the name and wrong for the
        table: with no adduct column a reader cannot tell which ion `Quant Ion` measures, and
        cannot check whether two rows at one retention time are two adducts of one molecule —
        which they can be, since `Na - H` and `+2 carbons, +3 double bonds` differ by 2.4 mDa.
        Reported in its own column so `Identification` stays exactly what LipiDex writes.

        Scoped to the same candidates as `identification()`, on the same principle as
        `identification_scores`: a row reported at sum composition takes the ions of every
        candidate that collapses to that sum, a row reported at molecular resolution only the
        top candidate's.

        Ordered by how near the supporting spectra's precursor sits to `Quant Ion`, so the ion
        the row is actually quantified on leads; ties broken by number of spectra. More than one
        means the molecule was seen as several ions, not that the assignment is uncertain.
        """
        name, _ = self.identification()
        if not name:
            return ""
        if not self.lipid_candidates:
            # A retention-model row has no candidates by construction, but its ion is known: the
            # model matched one m/z, and `extend_identifications` declines anything ambiguous.
            return self.rtls_adduct
        if name == self.sum_id or self.plasmenyl_ether_conflict:
            candidates = self.lipid_candidates
        else:
            candidates = self.lipid_candidates[:1]
        seen: dict[str, tuple[float, int]] = {}
        for candidate in candidates:
            for lipid in candidate.identifications:
                if not lipid.adduct:
                    continue
                gap = (abs(lipid.precursor - self.quant_ion)
                       if self.quant_ion is not None else 0.0)
                nearest, count = seen.get(lipid.adduct, (gap, 0))
                seen[lipid.adduct] = (min(nearest, gap), count + 1)
        return "; ".join(sorted(seen, key=lambda a: (seen[a][0], -seen[a][1])))

    def theoretical_mz(self) -> float | None:
        """The library mass of the ion `quant_ion` measures, or None when nothing was matched.

        This is what makes mass error reportable per row, and calibration measurable per file. It
        was already computed at match time (`LipidIdentification.lib_precursor`) and then thrown
        away, so every downstream consumer had to re-derive it by looking the name up in the
        library — which fails for about seven identifications in eight, because the table carries a
        canonicalised name (`SM d30:1`) and the library the raw MSP entry (`SM d14:1_16:0`).

        The ion is chosen the way `adduct()` chooses it: nearest to the quantified m/z, so the
        value describes the peak the row is actually measured on rather than some other adduct of
        the same molecule.
        """
        name, _ = self.identification()
        if not name or not self.lipid_candidates:
            return None
        candidates = (self.lipid_candidates if name == self.sum_id or self.plasmenyl_ether_conflict
                      else self.lipid_candidates[:1])
        best = None
        for candidate in candidates:
            for lipid in candidate.identifications:
                mz = getattr(lipid, "lib_precursor", None) or getattr(lipid, "precursor", None)
                if not mz:
                    continue
                gap = abs(mz - self.quant_ion) if self.quant_ion is not None else 0.0
                if best is None or gap < best[0]:
                    best = (gap, float(mz))
        return best[1] if best else None

    def mass_error_ppm(self) -> float | None:
        """Observed minus theoretical, in ppm. None when either side is unknown."""
        theoretical = self.theoretical_mz()
        if theoretical is None or not self.quant_ion:
            return None
        return (self.quant_ion - theoretical) / theoretical * 1e6

    def identification(self) -> tuple[str, str]:
        """Name and class, at the resolution the evidence supports — CDCompoundGroup.toString.

        ⚠ **The class is derived from the name wherever the name has been demoted**, rather than
        carried from the winning candidate. Otherwise the two disagree: the name says `PC O-38:6`,
        declining to assert a plasmalogen, while the class still says `Plasmenyl-PC` — so anyone
        filtering by class counts a plasmalogen the name refuses to claim. On one study that is 18
        rows across `Plasmenyl-PC` and `Plasmenyl-PE`.
        """
        if not self.lipid_candidates:
            if self.rtls_identification:
                # ⚠ A retention-model row has NO fragmentation behind it, so it cannot prove an
                # alk-1-enyl bond — and `P-` under the LSI shorthand is exactly that claim. Six of
                # the 14 `Alkenyl-TG` rows on one study asserted a proven vinyl ether on retention
                # time alone. The existing ether guard could never reach them: it matches on the
                # substrings `Plasmenyl`/`Plasmanyl`, which `Alkenyl-TG` does not contain.
                demoted = _demote_proven_ether(self.rtls_identification)
                return demoted, demoted.split(" ", 1)[0]
            return "", ""
        lipid_class = self.lipid_candidates[0].lipid_class
        if self.plasmenyl_ether_conflict:
            return self.top_ether(), _class_of(self.top_ether(), lipid_class)
        if self.purity == 0.0 and "Plasmenyl" in self.sum_id:
            return self.top_ether(), _class_of(self.top_ether(), lipid_class)
        if self.purity < MIN_FA_PURITY:
            return self.sum_id, lipid_class
        # ⚠ `self.purity` is the GROUP's dominant purity — the best-supported component in the
        # peak, which need not be the one about to be named. So a chain-resolved name could ship
        # while the purity attributable to THAT name was zero: the group had good fragment
        # evidence for some lipid, and the row claimed chains for a different one.
        #
        # Measured across three experiments, `Purity = 0` on a chain-resolved name is the single
        # most reliable predictor of a false identification here. It is what `Cer[AP] t20:0_23:0`
        # carried while matching an acetate adduct of an ordinary ceramide at dot product 999, and
        # what the held-out ether phosphatidylcholines carried at dot product 1000.
        #
        # Narrow on purpose. Purity 0 is the NORM for a sum composition — 328 of 509 rows on one
        # study — because a sum name makes no chain claim, and demoting those would destroy most
        # identifications for no gain. Only a name that ASSERTS chains while nothing supports them
        # is a contradiction, and only that is demoted: 6 rows of 82 chain-resolved on that study.
        top = self.lipid_candidates[0]
        if getattr(top, "purity", 0.0) < MIN_FA_PURITY and top.lipid_name != self.sum_id:
            return self.sum_id, lipid_class
        return top.lipid_name, lipid_class

    def purity_source(self) -> str:
        """`self`, or the name of the lipid this row's `Purity` actually belongs to.

        `calc_purity` drops a top hit that scores nothing and then returns `purities[0]`, which by
        then belongs to a different candidate — so a row can carry another lipid's purity, at any
        value including 100. That is LipiDex's behaviour, verified against `PeakPurity.java` down
        to a guard its author wrote and commented out, so it is reported rather than repaired.

        On one study this affects **431 positive and 341 negative** rank-1 matches, and the worst
        cases are cross-class: `PS 20:4_22:0` at Purity 100 from a `PC`. `Chain Evidence` cannot
        express those — it only asks whether a component names this row's chains, and is blank for
        sum compositions. This covers every row.

        A component that collapses to the same sum composition is still `self`: a row that declined
        to name its chains has not been scored by a different molecule, only by a better-resolved
        view of its own.
        """
        if not self.summed_purities:
            return ""
        name, _ = self.identification()
        if not name:
            return ""
        top = _ADDUCT_SUFFIX.sub("", self.summed_purities[0][0].strip()).rstrip(";").strip()
        if top == name or _chain_key(top) == _chain_key(name):
            return "self"
        mine, theirs = _sum_key(name), _sum_key(top)
        if mine is not None and mine == theirs:
            return "self"
        return top


@dataclass(slots=True)
class LipidCandidate:
    """One molecular identity within a compound group, pooled over the spectra supporting it."""
    lipid_name: str = ""
    sum_lipid_name: str = ""
    lipid_class: str = ""
    max_dot: float = 0.0
    max_rev_dot: float = 0.0
    count: int = 1
    purity: int = 0
    identifications: list[Lipid] = field(default_factory=list)
    preferred_polarity: bool = False

    def __init__(self, lipid: Lipid):
        self.lipid_name = lipid.lipid_name
        self.sum_lipid_name = lipid.sum_lipid_name
        self.lipid_class = lipid.lipid_class
        self.max_dot = lipid.dot
        self.max_rev_dot = lipid.rev_dot
        self.count = 1
        self.identifications = [lipid]
        self.purity = lipid.purity
        self.preferred_polarity = lipid.preferred_polarity

    def add(self, lipid: Lipid):
        self.max_dot = max(self.max_dot, lipid.dot)
        self.max_rev_dot = max(self.max_rev_dot, lipid.rev_dot)
        self.identifications.append(lipid)
        # Lipid.compareTo: preferred polarity first, then how far up the peak it was taken.
        self.identifications.sort(key=lambda l: (not l.preferred_polarity, -l.gaussian_score))
        self.count += 1
        self.purity = self.identifications[0].purity
        self.preferred_polarity = self.identifications[0].preferred_polarity

    def __lt__(self, other: "LipidCandidate") -> bool:
        """Preferred polarity wins, then the identification with more spectra behind it."""
        if self.preferred_polarity != other.preferred_polarity:
            return self.preferred_polarity
        return self.count > other.count
