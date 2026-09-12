"""Two identified rows that are two adducts of one molecule — where at most one name is right.

The peak finder's adduct sweep only ever fires when one member of a pair is **unidentified**:
`checkAdduct` against a named peak, or between two unnamed peaks. There is no branch for two
identified groups, so a pair that both carry names is never compared. That hole is not academic,
because of one coincidence:

    Na - H                        = 21.98194
    +2 carbons, +3 double bonds   = 21.98435   (C2H4 - 3 x H2)
    difference   2.40 mDa   =   3.0 ppm at m/z 800

**The sodium adduct of `PE 40:1` sits 3 ppm from protonated `PE 42:4`** — inside any ordinary
search window, at the same retention time, in a class that supplies both rungs of the ladder. Both
rows get named, both survive, and their neutral masses agree, which means **at most one of the two
names can be true**.

## Why this removes rather than sums

The opposite of `split_peaks.py`, and the distinction is the whole point. A split peak is one
molecule divided between two rows, so the rows are summed. An adduct pair is one molecule and one
**ghost**: adding them would add a real measurement to an artefact, and would also add two ions
with two response factors. The loser has to go.

## Mass agreement is not enough on its own

Measured on the brain positive run, **412** differently-named pairs at one peak have masses
consistent with two adducts of one molecule — but only **14** also track across the samples. With
eight positive adducts the mass test alone is permissive and co-elution inside a lipid class is
ordinary, so mass agreement is a screen, not a verdict.

Correlation is the confirmation here, and note that this is the reverse of `split_peaks.py`, where
correlation is actively misleading. Two adducts of one molecule are two ions of the same eluting
compound, so they must rise and fall together; the adduct ratio drifts with the matrix, which adds
noise but cannot make them disagree. A split peak has no such requirement — its halves are
anti-correlated, because the integration boundary moves. Same statistic, opposite meaning, because
the underlying physical claim is different.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .adducts import Adduct, load_adducts, ppm_diff
from .correlate import correlate_areas
from .peaks import sum_composition

# 13C - 12C. `[M]+` and `[M+H]+` differ by a proton, 1.00728, which sits only 3.9 mDa from this —
# under 5 ppm at m/z 800, inside any sane window. Without this guard an M+1 isotope that survived
# the isotope filter is removed as an "adduct", and the audit found exactly that: `SM d43:2`
# deleted against `HexCer[NS] d43:1` on a gap of 1.0039.
ISOTOPE_SPACING = 1.0033548
ISOTOPE_GUARD_DA = 0.01

# Adducts that a class is normally quantified on. A sodiated or ammoniated species that has to
# outrank one of these is making the weaker claim: it fragments poorly, so its identification
# rests on less evidence even when the dot product happens to be comparable.
CANONICAL = {
    "TG": ("[M+NH4]+",), "DG": ("[M+NH4]+",), "CE": ("[M+NH4]+",),
    "PC": ("[M+H]+",), "PE": ("[M+H]+",), "PS": ("[M+H]+",), "PI": ("[M+H]+",),
    "SM": ("[M+H]+",), "Cer": ("[M+H]+",),
}
DEFAULT_CANONICAL = ("[M+H]+", "[M-H]-")


@dataclass
class AdductPairFilter:
    """Remove the weaker of two identified rows that are two adducts of one molecule.

    `min_correlation` is the confirmation the mass test cannot give. On this data mass agreement
    alone flags 412 positive pairs and only 14 survive tracking, so this is doing most of the work.

    `max_peak_widths` is in units of the run's mean FWHM, as everywhere else, so it travels
    between gradients.
    """

    enabled: bool = True
    max_peak_widths: float = 1.0
    max_ppm: float = 10.0
    min_correlation: float = 0.8
    min_points: int = 8
    correlation_type: str = "spearman"


@dataclass
class AdductPairReport:
    considered: int = 0
    removed: int = 0
    rejected_by_correlation: int = 0
    rejected_isotope: int = 0
    ambiguous: int = 0
    same_molecule: int = 0
    examples: list = field(default_factory=list)

    def __str__(self) -> str:
        out = (f"adduct pairs: {self.removed} removed of {self.considered} mass-consistent pairs "
               f"({self.rejected_by_correlation} rejected for not tracking, "
               f"{self.ambiguous} ambiguous, {self.rejected_isotope} isotope spacings skipped; "
               f"{self.same_molecule} were one molecule on two ions)")
        if self.examples:
            name, winner, adduct = self.examples[0]
            out += f" (e.g. {name} removed as the {adduct} of {winner})"
        return out


def _canonical(lipid_class: str, adduct: str) -> bool:
    return adduct in CANONICAL.get(lipid_class, DEFAULT_CANONICAL)


def _neutral(group, adducts: list[Adduct]) -> float | None:
    """The neutral mass implied by the row's own quant ion and its declared adduct."""
    declared = (group.adduct() or "").split("; ")[0]
    if not declared or group.quant_ion is None:
        return None
    match = next((a for a in adducts if a.name == declared), None)
    return group.quant_ion - match.mass if match else None


def _explains(neutral: float, other, adducts: list[Adduct], max_ppm: float) -> str | None:
    """The adduct by which `other`'s quant ion would be this neutral molecule, if any."""
    if other.quant_ion is None:
        return None
    for adduct in adducts:
        if adduct.polarity != other.quant_polarity:
            continue
        if ppm_diff(other.quant_ion - adduct.mass, neutral) < max_ppm:
            return adduct.name
    return None


def remove_adduct_pairs(groups, samples, params: AdductPairFilter, avg_fwhm: float,
                        adducts: list[Adduct] | None = None, log=None) -> AdductPairReport:
    """Drop the weaker row of each confirmed adduct pair. Mutates `groups`.

    Which row survives, in order: the one whose adduct is canonical for its class, then the higher
    dot product. Canonical first rather than score first because a sodiated identification is the
    weaker claim on chemistry — the species barely fragments — and a dot product earned on few
    fragments is not comparable to one earned on many.
    """
    say = log or (lambda _: None)
    report = AdductPairReport()
    if not params.enabled:
        return report
    table = adducts if adducts is not None else load_adducts()
    if not table:
        say("adduct pairs: no adduct database, skipped")
        return report

    columns = [i for i, s in enumerate(samples) if s.role != "blank"]
    window = avg_fwhm * params.max_peak_widths
    live = sorted((g for g in groups
                   if g.keep and g.quant_ion is not None and g.identification()[0]),
                  key=lambda g: g.retention)

    for i, first in enumerate(live):
        if not first.keep:
            continue
        for second in live[i + 1:]:
            if second.retention - first.retention > window:
                break
            if not second.keep or not first.keep:
                continue
            name_a, class_a = first.identification()
            name_b, class_b = second.identification()
            # ⚠ On SUM COMPOSITION, not the displayed name. `TG 18:2_18:1_20:1` and `TG 56:4` are
            # one molecule reported at two resolutions, and comparing the strings called them
            # different — which is precisely what `sum_composition` exists to prevent.
            same_molecule = sum_composition(name_a) == sum_composition(name_b)
            same_adduct = ((first.adduct() or "").split("; ")[0]
                           == (second.adduct() or "").split("; ")[0])
            if same_molecule and same_adduct:
                continue                      # one peak cut in two: split_peaks.py's business
            if abs(abs(first.quant_ion - second.quant_ion) - ISOTOPE_SPACING) < ISOTOPE_GUARD_DA:
                report.rejected_isotope += 1
                continue
            neutral_a, neutral_b = _neutral(first, table), _neutral(second, table)
            explained = (neutral_a is not None
                         and _explains(neutral_a, second, table, params.max_ppm)) or None
            reversed_ = (neutral_b is not None
                         and _explains(neutral_b, first, table, params.max_ppm)) or None
            if not explained and not reversed_:
                continue
            report.considered += 1
            r = correlate_areas(first.areas, second.areas, columns,
                                params.correlation_type, params.min_points)
            if r is None or r < params.min_correlation:
                report.rejected_by_correlation += 1
                continue

            # ⚠ The direction comes from the MASS EVIDENCE, not from the ranking. `explained`
            # means "second's ion is an adduct of first's molecule", and that is the only claim
            # the masses support — it licenses removing `second` and nothing else. Ranking the
            # pair independently and then removing whichever scored worse deleted rows on a
            # relationship that did not describe them, and left the reason string reading
            # "the None of ...", which is how the audit found it.
            if explained and reversed_:
                winner, loser = _rank(first, class_a, second, class_b)
                adduct_of = explained if loser is second else reversed_
            elif explained:
                winner, loser, adduct_of = first, second, explained
            else:
                winner, loser, adduct_of = second, first, reversed_
            # If the row the masses accuse is also the better-supported one, the pair is not what
            # it looks like. Skipping costs a true positive; removing costs a real lipid.
            if _rank(winner, winner.identification()[1],
                     loser, loser.identification()[1])[0] is loser:
                report.ambiguous += 1
                continue
            loser.keep = False
            # Two genuinely different cases, and the reason has to say which. Different molecules
            # means one of the two NAMES is wrong. The same molecule on two ions means both names
            # are right and one ROW is redundant — nothing was misidentified.
            kind = ("the same molecule on its" if same_molecule else "the")
            loser.filter_reason = (f"Adduct pair: {kind} {adduct_of} of "
                                   f"{winner.identification()[0]} at {winner.retention:.3f} min")
            if same_molecule:
                report.same_molecule += 1
            report.removed += 1
            report.examples.append((loser.identification()[0],
                                    winner.identification()[0], adduct_of))
    say(str(report))
    return report


def _rank(first, class_a: str, second, class_b: str):
    """(winner, loser). Canonical adduct first, then dot product."""
    canonical_a = _canonical(class_a, (first.adduct() or "").split("; ")[0])
    canonical_b = _canonical(class_b, (second.adduct() or "").split("; ")[0])
    if canonical_a != canonical_b:
        return (first, second) if canonical_a else (second, first)
    dot_a = first.identification_scores()[0] or 0
    dot_b = second.identification_scores()[0] or 0
    return (first, second) if dot_a >= dot_b else (second, first)
