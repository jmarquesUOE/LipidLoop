"""`Spectral_Components.csv` — how the chain-bearing signal in a peak divides between candidates.

`Purity` reports one number: the share of the moiety-fragment signal belonging to the winning
name. That number is the top of a *breakdown* — every candidate that explained some of the signal,
with its share — which `calculate_weighted_purity` computes and, until now, discarded. This writes
it out.

**Why it is worth writing.** Asked how a sum composition divides between its chain combinations —
`PC 38:3` as so much `18:0_20:3` and so much `18:1_20:2` — the results table can only answer where
the isomers happened to separate chromatographically *and* each triggered MS2 *and* each won a dot
product on its own peak. On one study that is 8 of 155 chain-resolved compositions, and the
selection runs the wrong way: a minor isomer fails all three conditions, so near-equal splits
survive and lopsided ones are reported as a single peak. The breakdown is a *within-peak*
measurement and needs none of that — co-eluting isomers are exactly the case it describes.

**`Isomer`** marks components collapsing to the same sum composition as the reported name. The
breakdown contains isobaric candidates from other classes as well, which are a different question
(what else is in this peak) from chain composition (how this lipid's chains divide), and mixing
them would make a chain distribution out of a co-isolation.

⚠⚠ **A share is not a molar fraction, and this file must not be read as one.** It is a share of
*fragment intensity*, and fragment intensity is not proportional to abundance: the *sn*-2
carboxylate is preferentially released in glycerophospholipids, and response varies with chain
length and unsaturation. Turning these into composition requires calibration against authentic
standards. The column is named `Share of Fragment Signal (%)` rather than anything shorter for
that reason.

⚠⚠ **It is worth least exactly where the library enumerates most.** A peak's components are
library entries that explained some of its chain-bearing fragments, and for triacylglycerols the
library holds nearly every arithmetic combination — so `TG 58:6` returns 27 components down to
`TG 10:2_24:1_24:1`, which is the combinatorial library re-expressed rather than a measurement.
Two-chain classes with discriminating carboxylate fragments (PE, PI, PG, PS) are the case this
file answers well. **Read the component count as a warning: a peak with 20 of them has not been
resolved into a distribution, it has been spread across a catalogue.**

⚠ Absence is not evidence of absence either. A candidate reaches the breakdown only if *all* of its
chains were found — `purity.py` returns nothing for a candidate whose second chain fragment fell
below the 5% intensity floor, so it is missing from the denominator rather than scored low.
"""
from __future__ import annotations

import csv
import re
from pathlib import Path

from .peakfinder import PeakFinderResult
from .rtls import parse_sum_name

COLUMNS = ["Compound Group", "Retention Time (min)", "Identification", "Lipid Class",
           "Component", "Share of Fragment Signal (%)", "Isomer", "Rank"]

# `TG 18:1_18:1_22:4 [M+NH4]+` and the same chains as `[M+Na]+` are two library entries for ONE
# combination, and left alone they appear as two components splitting one share between them —
# which reads as two isomers at half the abundance each. Merged on the chains.
_ADDUCT = re.compile(r"\s*\[M[^\]]*\][+-]?\s*$")


def _bare(name: str) -> str:
    return _ADDUCT.sub("", (name or "").strip())


def _composition(name: str):
    parsed = parse_sum_name(name)
    return parsed if parsed else None


def write_spectral_components(result: PeakFinderResult, path: str | Path) -> tuple[int, int]:
    """Write the per-name purity breakdown for every group that has one.

    Returns (components written, peaks they came from).

    Numbered over the full set exactly as `write_results` numbers it, so `Compound Group` joins
    against every other table.
    """
    ids = {id(g): n for n, g in enumerate(result.compound_groups, start=1)}
    rows = peaks = 0
    with Path(path).open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(COLUMNS)
        for group in result.compound_groups:
            merged: dict[str, float] = {}
            for component, value in group.summed_purities:
                bare = _bare(component)
                merged[bare] = merged.get(bare, 0.0) + value
            # A single component is the whole of the signal and says nothing: it is what an
            # unambiguous peak looks like, and writing it would triple the file to no purpose.
            if len(merged) < 2:
                continue
            # ⚠ NOT `purity_weight`. That denominator counts the Gaussian score once per ENTRY,
            # so it grows with the number of candidates and the shares under it sum to well under
            # 100 — a median of 41% on one run, and less the more candidates a peak has. That is
            # the right denominator for purity, whose question is "how much belongs to the
            # winner", and the wrong one for a distribution, whose shares must be comparable
            # between peaks. Normalised within the peak instead, so the column sums to 100.
            total = sum(merged.values())
            if total <= 0:
                continue
            name, lipid_class = group.identification()
            winner = _composition(name)
            peaks += 1
            ordered = sorted(merged.items(), key=lambda p: -p[1])
            for rank, (component, value) in enumerate(ordered, start=1):
                other = _composition(component)
                writer.writerow([
                    ids.get(id(group), ""), group.retention, name, lipid_class, component,
                    round(value / total * 100.0, 2),
                    "yes" if winner and other and winner == other else "no",
                    rank,
                ])
                rows += 1
    return rows, peaks
