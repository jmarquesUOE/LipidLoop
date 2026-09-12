"""Detect instrument artefacts in MS2 spectra, automatically, on every run.

A peak that appears whatever the precursor was is not a fragment. It matters because every
spectrum is normalised to its base peak before scoring: an artefact that *is* the base peak
scales every genuine fragment down against it, and its unmatched intensity counts against the
forward dot product too. On the instrument this was written against, one such peak costs 36% of
identifications in positive mode and 50% in negative.

**This screens on every analysis rather than relying on a configured list.** A hard-coded m/z
would be right for one instrument in one period and wrong everywhere else, which is no use to
anyone running the method in another lab. The detector has to travel with the method.

## Telling an artefact from a real ubiquitous fragment

The danger of screening automatically is deleting something real. Lipidomics has genuinely
common fragments — phosphocholine at 184.0733 is in every PC and SM spectrum — and removing one
would be far worse than the artefact it was meant to catch. Three conditions must hold together:

1. **Present in nearly every spectrum** regardless of precursor. Real class fragments are not:
   phosphocholine appears in 9.9% of positive spectra here, because only some precursors are
   phosphocholines. An artefact appears in 100%.
2. **Usually the base peak.** This is what makes it damaging rather than merely present, and it
   is the condition that justifies acting at all.
3. **No possible formula.** The decisive test. A singly-charged ion built from CHNOPS cannot have
   an arbitrary mass defect: the richest possible in hydrogen is a saturated alkyl cation, whose
   defect runs at about 0.00115 per unit m/z. At m/z 178 nothing real can exceed ~0.23, and the
   artefact here sits at 0.290. Meanwhile every genuine fragment checked — phosphocholine 0.074,
   choline 0.108, sphingosine 0.261, cholestadiene 0.344 — is comfortably inside the bound.

Condition 3 is what makes this safe to run unattended. Frequency alone would eventually delete a
real fragment on a sample type where one class dominates; an impossible mass cannot be real
whatever the sample.

Anything removed is reported, never dropped silently.
"""
from __future__ import annotations

import collections
import statistics
from dataclasses import dataclass

# Mass defect of a saturated alkyl cation per unit m/z — the most hydrogen-rich, and so the
# highest-defect, ion that CHNOPS chemistry can produce. The margin covers calibration error.
MAX_DEFECT_PER_MZ = 0.00115
DEFECT_MARGIN = 0.03

MIN_FRACTION = 0.95        # present in at least this share of spectra
MIN_BASE_FRACTION = 0.20   # and the base peak in at least this share
BIN = 0.10                 # wide: an artefact's m/z wanders, and a narrow bin splits it


def plausible_mass(mz: float, defect_margin: float = DEFECT_MARGIN) -> bool:
    """Could a singly-charged CHNOPS ion have this m/z?"""
    defect = mz - int(mz)
    return defect <= MAX_DEFECT_PER_MZ * mz + defect_margin


@dataclass(slots=True)
class Artefact:
    mz: float
    fraction: float
    base_fraction: float
    spread_mda: float

    def __str__(self) -> str:
        return (f"m/z {self.mz:.3f} in {100*self.fraction:.0f}% of spectra, base peak in "
                f"{100*self.base_fraction:.0f}%, m/z spread {self.spread_mda:.0f} mDa, "
                f"no possible formula")


def screen(spectra, min_fraction: float = MIN_FRACTION,
           min_base_fraction: float = MIN_BASE_FRACTION, bin_width: float = BIN,
           protected: "list[float] | tuple[float, ...]" = ()) -> list[Artefact]:
    """Find artefacts in a set of spectra. `spectra` need only expose `.mz` and `.intensity`.

    `protected` exempts masses you expect to see. **Deuterated standards need this.** The decisive
    test asks whether a singly-charged CHNOPS ion could carry the observed mass defect, and
    deuterium is not in CHNOPS: d31-palmitate at m/z 286.4275 and d35-stearate at 318.4839 have
    defects of 0.43 and 0.48 against limits of 0.36 and 0.40, so a spiked standard mix lands in
    exactly the rejection region the screen was built to catch, and would be deleted per file with
    only a line in the log to say so. The d7 species are unaffected — seven deuteriums is not
    enough to trip it — which makes the failure selective and easy to miss.

    Passing the standards you spiked is therefore not a convenience. Without it the screen is
    wrong about them by construction.
    """
    present = collections.Counter()
    is_base = collections.Counter()
    observed = collections.defaultdict(list)
    total = 0

    for spectrum in spectra:
        if not spectrum.mz:
            continue
        total += 1
        top_intensity, top_mz = max(zip(spectrum.intensity, spectrum.mz))
        is_base[round(top_mz / bin_width)] += 1
        for key in {round(mz / bin_width) for mz in spectrum.mz}:
            present[key] += 1
        for mz in spectrum.mz:
            observed[round(mz / bin_width)].append(mz)

    if total == 0:
        return []

    found = []
    for key, count in present.items():
        fraction = count / total
        if fraction < min_fraction:
            continue
        base_fraction = is_base[key] / total
        if base_fraction < min_base_fraction:
            continue
        values = observed[key]
        mean_mz = statistics.mean(values)
        if plausible_mass(mean_mz):
            continue          # it could be a real ion — leave it alone
        if any(abs(mean_mz - m) <= bin_width for m in protected):
            continue          # expected, and declared: a spiked standard, not an artefact
        found.append(Artefact(mz=mean_mz, fraction=fraction, base_fraction=base_fraction,
                              spread_mda=1000 * statistics.stdev(values) if len(values) > 2 else 0.0))
    found.sort(key=lambda a: -a.fraction)
    return found


def strip(spectra, artefacts, tolerance: float = BIN) -> int:
    """Remove artefact peaks from spectra in place. Returns the number of peaks removed."""
    masses = [a.mz for a in artefacts]
    if not masses:
        return 0
    removed = 0
    for spectrum in spectra:
        if not spectrum.mz:
            continue
        keep = [i for i, mz in enumerate(spectrum.mz)
                if not any(abs(mz - bad) <= tolerance for bad in masses)]
        removed += len(spectrum.mz) - len(keep)
        spectrum.mz = [spectrum.mz[i] for i in keep]
        spectrum.intensity = [spectrum.intensity[i] for i in keep]
    return removed
