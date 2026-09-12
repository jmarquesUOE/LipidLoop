"""Spiked standards: what to expect, and how the instrument performed against them.

A pooled QC is made from the study and answers "was this batch repeatable". A standard mix is the
same material in every batch and every year, so it answers a different question — "is the
instrument where it was last month" — and its absolute response is comparable across studies in a
way a pool's never is. Both are needed; neither substitutes for the other.

Two jobs here:

  * declare the expected ions, so the artefact screen can be told they are expected. Heavily
    deuterated standards *must* be declared: the screen rejects masses whose defect is impossible
    for a singly-charged CHNOPS ion, deuterium is not CHNOPS, and d31/d35 fatty acids land inside
    the rejection region. See `artefacts.screen`.
  * measure each standard in each injection — response, mass accuracy, retention — and write it
    out, so performance can be tracked between runs rather than only within one.

Masses are computed from formulae here rather than copied from a supplier sheet. That is not
fussiness: the sheet this was built from listed the two deuterated fatty acids at their
*unlabelled* neutral masses, 256.2402 for C16D31H1O2 and 284.2715 for C18D35H1O2, which are wrong
by 31.19 and 35.22 Da. A formula cannot make that mistake.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean, median, stdev

def measured_mz(group) -> "float | None":
    """The mass as the instrument reported it, before any systematic correction.

    Standards exist to judge the instrument independently of the software. Once `mass_offset_ppm`
    began correcting feature masses, measuring the standards on `quant_ion` meant measuring the
    correction rather than the instrument: the standards-derived offset moved from -7.16 to -1.39
    the moment the correction was applied, which is the residual, and the cross-check then compared
    a residual against a residual with the offset added back to one side only.
    """
    return getattr(group, "quant_ion_measured", None) or group.quant_ion


MONOISOTOPIC = {
    "C": 12.0, "H": 1.00782503207, "D": 2.01410177785, "N": 14.0030740048,
    "O": 15.9949146196, "P": 30.97376163, "S": 31.97207100, "Na": 22.98976928,
}
PROTON = 1.007276466
ELECTRON = 0.000548580

# Adduct name -> mass added to the neutral. Deliberately small: these are the ions a lipid
# standard mix actually gives, not a general adduct table.
ADDUCTS = {
    "[M+H]+": PROTON,
    "[M+NH4]+": MONOISOTOPIC["N"] + 4 * MONOISOTOPIC["H"] - ELECTRON,
    "[M+Na]+": MONOISOTOPIC["Na"] - ELECTRON,
    "[M-H]-": -PROTON,
    "[M+HCOO]-": MONOISOTOPIC["C"] + MONOISOTOPIC["H"] + 2 * MONOISOTOPIC["O"] + ELECTRON,
    "[M-CH3]-": -(MONOISOTOPIC["C"] + 3 * MONOISOTOPIC["H"]) + ELECTRON,
}

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def formula_mass(formula: str) -> float:
    """Monoisotopic neutral mass. `D` is deuterium, not an element symbol collision."""
    total = 0.0
    for element, count in _TOKEN.findall(formula.replace(" ", "")):
        if not element:
            continue
        if element not in MONOISOTOPIC:
            raise ValueError(f"unknown element {element!r} in {formula!r}")
        total += MONOISOTOPIC[element] * int(count or 1)
    return total


@dataclass(slots=True)
class Standard:
    name: str
    formula: str
    adduct: str
    polarity: str            # "+" or "-"
    retention: float = 0.0   # optional, filled from observation

    @property
    def neutral_mass(self) -> float:
        return formula_mass(self.formula)

    @property
    def mz(self) -> float:
        if self.adduct not in ADDUCTS:
            raise ValueError(f"unknown adduct {self.adduct!r} for {self.name}")
        return self.neutral_mass + ADDUCTS[self.adduct]


@dataclass
class StandardMix:
    """The spiked mix. `protected_mz` is what the artefact screen must be told about."""

    standards: list[Standard] = field(default_factory=list)

    @classmethod
    def from_csv(cls, path: str | Path) -> "StandardMix":
        """`name,formula,adduct,polarity` — mass is computed, never read."""
        out = []
        for row in csv.DictReader(open(path)):
            if not (row.get("name") or "").strip():
                continue
            out.append(Standard(name=row["name"].strip(), formula=row["formula"].strip(),
                                adduct=row["adduct"].strip(), polarity=row["polarity"].strip()))
        return cls(out)

    def for_polarity(self, polarity: str) -> list[Standard]:
        return [s for s in self.standards if s.polarity == polarity]

    def protected_mz(self, polarity: str = "") -> list[float]:
        chosen = self.for_polarity(polarity) if polarity else self.standards
        return [s.mz for s in chosen]


#: ⚠ **A pre-run scan of raw MS1 cannot tell a spiked study from an unspiked one.** An attempt at
#: `detect_in_files` — search each standard's exact mass across MS1 scans, call the mix present on a
#: quorum — was written and removed, because it does not work:
#:
#:     Skin_QEplus    spiked      1.1-2.1% of MS1 scans contain the mass
#:     ST004797         not spiked  4.2-13.9%
#:
#: The known-spiked study scores LOWEST. In a dense spectrum some peak sits within 10 ppm of almost
#: any mass, so raw-scan presence measures spectral density, not the analyte.
#:
#: What does work is the FEATURE table: the same standards appear there at +0.1 to +2.4 ppm as
#: detected features. So the quorum has to be applied AFTER the run, on features — not before it,
#: on scans. That is why the ISTD library cannot simply be gated on a cheap pre-check.
#:
#: Thresholds, when it is implemented on features: a mix is present when >= 2 standards are each
#: detected in >= half the injections. On Skin_QEplus three standards were each found in 58-60 of
#: 65; on unspiked samples the same test gives 1 of 5 at p = 0.134, the single hit being the d7-TG
#: mass that matches at baseline nearly everywhere.


MIN_STANDARDS = 2
MIN_INJECTION_FRACTION = 0.5


def detect_in_features(mix: StandardMix, polarity: str, groups, samples,
                       ppm: float = 10.0, min_standards: int = MIN_STANDARDS,
                       min_fraction: float = MIN_INJECTION_FRACTION) -> dict:
    """Is this mix actually spiked into the samples? Decided on FEATURES, never on raw scans.

    Unlike `measure`, this looks at ordinary sample injections rather than standards injections —
    the question is whether an internal standard is present *in the samples*, which is what makes
    it usable for mass accuracy and recovery, and is a different fact from a standards vial having
    been run.

    ⚠ **One standard found is the expected background, not evidence.** The d7-TG mass in particular
    matches at baseline in nearly every dataset tested, so a single hit says nothing. Both
    thresholds must hold: at least `min_standards` distinct compounds, each detected in at least
    `min_fraction` of the injections. A real spike is consistent across a run; a coincidence is not.

    Returns {present, standards_found, injections, per_standard}. `per_standard` gives the
    detection count for every compound in the mix, present or absent, so a near-miss is visible
    rather than collapsed into a boolean.
    """
    wanted = mix.for_polarity(polarity)
    columns = [i for i, s in enumerate(samples) if getattr(s, "role", "sample") == "sample"]
    if not wanted or not columns:
        return {"present": False, "standards_found": 0, "injections": len(columns),
                "per_standard": {}}

    seen: dict[str, int] = {}
    for standard in wanted:
        tolerance = standard.mz * ppm / 1e6
        best = 0
        for group in groups:
            mz = measured_mz(group)
            if mz is None or abs(mz - standard.mz) > tolerance:
                continue
            areas = getattr(group, "areas", None)
            if areas is None:
                best = max(best, len(columns))     # no per-injection detail: count it as present
                continue
            best = max(best, sum(1 for i in columns
                                 if i < len(areas) and areas[i] and areas[i] > 0))
        seen[standard.name] = best

    need = max(1, int(round(len(columns) * min_fraction)))
    found = sum(1 for n in seen.values() if n >= need)
    return {"present": found >= min_standards, "standards_found": found,
            "injections": len(columns), "per_standard": seen}


def measure(mix: StandardMix, polarity: str, groups, samples, ppm: float = 10.0,
            rt_tolerance: float = 0.5) -> list[dict]:
    """Response, mass accuracy and retention for each standard in each standards injection.

    Matched against the compound groups the pipeline already built, so this reports the same
    numbers the result tables do rather than a second, differently-derived measurement.
    """
    columns = [i for i, s in enumerate(samples) if s.role == "standard"]
    if not columns:
        return []

    out = []
    for standard in mix.for_polarity(polarity):
        tolerance = standard.mz * ppm / 1e6
        hits = [g for g in groups
                if measured_mz(g) and abs(measured_mz(g) - standard.mz) <= tolerance
                and (not standard.retention
                     or abs(g.retention - standard.retention) <= rt_tolerance)]
        if not hits:
            out.append({"standard": standard.name, "adduct": standard.adduct,
                        "expected_mz": round(standard.mz, 4), "found": False})
            continue
        best = max(hits, key=lambda g: g.max_area)
        areas = [best.areas[i] for i in columns if i < len(best.areas)]
        present = [a for a in areas if a > 0]
        row = {
            "standard": standard.name, "adduct": standard.adduct,
            "expected_mz": round(standard.mz, 4), "found": True,
            "observed_mz": round(measured_mz(best), 4),
            "mass_error_ppm": round(1e6 * (measured_mz(best) - standard.mz) / standard.mz, 2),
            "retention": round(best.retention, 3),
            "injections": len(columns), "detected_in": len(present),
            "median_area": round(median(present), 1) if present else 0.0,
        }
        if len(present) > 1:
            row["cv_percent"] = round(100 * stdev(present) / mean(present), 1)
        out.append(row)
    return out


def write(rows: list[dict], path: str | Path) -> int:
    if not rows:
        return 0
    fields: list[str] = []
    for row in rows:
        for key in row:
            if key not in fields:
                fields.append(key)
    with Path(path).open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    return len(rows)


def offset_ppm(mix: StandardMix, polarity: str, groups, ppm: float = 20.0) -> "float | None":
    """Median mass error over the spiked standards, in ppm.

    This is the *check*, not the calibration. Thousands of rank-1 library matches estimate the
    mass offset far better than a handful of standards ever could, so the offset applied to a run
    should be derived from the identifications. What the standards add is independence: the
    identification-derived offset is circular — it assumes the identifications are right — while a
    spiked standard's mass is known before the run starts. Agreement between the two is evidence
    the calibration is real; disagreement means one of them is wrong, and that is worth stopping
    for.

    Window is deliberately wide (20 ppm): the whole point is to measure an offset that may be
    large, and a tight window centred on zero would reject exactly the error being looked for.
    """
    errors = []
    for standard in mix.for_polarity(polarity):
        tolerance = standard.mz * ppm / 1e6
        hits = [g for g in groups
                if measured_mz(g) and abs(measured_mz(g) - standard.mz) <= tolerance]
        if not hits:
            continue
        best = max(hits, key=lambda g: g.max_area)
        errors.append(1e6 * (measured_mz(best) - standard.mz) / standard.mz)
    if not errors:
        return None
    return median(errors)


def check_offset(derived: "float | None", from_standards: "float | None",
                 tolerance_ppm: float = 3.0) -> dict:
    """Compare the offset derived from identifications against the spiked standards."""
    out = {"derived_ppm": derived, "standards_ppm": from_standards,
           "tolerance_ppm": tolerance_ppm}
    if derived is None or from_standards is None:
        out["verdict"] = "not comparable"
        return out
    out["difference_ppm"] = from_standards - derived
    out["agree"] = abs(out["difference_ppm"]) <= tolerance_ppm
    out["verdict"] = "agree" if out["agree"] else "DISAGREE"
    return out
