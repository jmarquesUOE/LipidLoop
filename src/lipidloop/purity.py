"""Identification purity — `spectrum_searcher/PeakPurity.java`.

A dot product says the spectrum looks like the library entry. It does not say the spectrum
belongs to *only* that lipid. At a given precursor mass several isobaric species can be
co-isolated, and each contributes its own fragments. Purity asks a different question: of the
signal that the candidate identifications between them explain, how much belongs to the winner?

The measure is built from moiety fragments — the peaks that carry chain information, an acyl
or sphingoid fragment or neutral loss, rather than a headgroup peak every member of the class
shares. For each candidate LipiDex chooses the most intense moiety fragment type present, finds
the sample intensity of every fragment of that type, and requires *all* of the lipid's chains to
be found or the candidate scores nothing. The winner's intensity as a percentage of the sum
across candidates is the purity, and 75% is the threshold in `Thresholds`.

Three behaviours are load-bearing and none is obvious:

  * **Only the top-ranked identification gets a purity**, and only if it is a LipiDex spectrum
    in its optimal polarity. Off-polarity entries score 0 — which is why positive-mode PC,
    annotated `OptimalPolarity=false`, comes back 0 throughout.
    ⚠ **The `is_lipidex` gate never fires in this configuration, and two earlier explanations of
    what it does were wrong.** `lib_gen/Lipid.java` writes `Type=LipiDex` as an unconditional
    string literal, so *every* library a LipiDex generator produces carries it — including ones
    built from LipidBlast class definitions. Measured: the flag is true for 100% of search rows in
    both polarities. There is **no by-design zero for LipidBlast**; off-polarity is the only one
    that exists. (The second wrong explanation blamed generic peak annotation leaving
    `transition_type` unset. Also false: 97.8% of LipidBlast entries can set it, because
    `Sphingoid Fragment` has a space before `Fragment` and survives a filter that excludes
    `_Fragment`.)
  * **A missing chain is fatal, not partial.** `faCount != faDB.size()` returns null and the
    candidate contributes nothing, so a lipid whose second chain fragment fell below the 5%
    intensity floor is not merely penalised, it is absent from the denominator.
  * **Shared fragments are charged once.** Isobaric candidates commonly predict the same
    fragment mass; an intensity correction accumulated across candidates stops that one peak
    being counted in full for each of them.

An unannotated peak in a LipiDex spectrum throws in the original and is caught into a 0.0. That
is reproduced rather than repaired: it changes which candidates enter the denominator.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import floor
from pathlib import Path

from .msp import LibrarySpectrum

MIN_PURITY_FOR_FRAGMENTS = 5.0   # PeakPurity.java: minPurity
MOIETY_INTENSITY_FLOOR = 0.05    # findFAIntensities is called with maxFragType * .05
MASS_EQUALITY = 1e-6             # getIndexOf


@dataclass(slots=True)
class FattyAcid:
    name: str
    base: str
    formula: str
    enabled: str


def read_fatty_acids(path: str | Path) -> list[FattyAcid]:
    """Read `src/backup/FattyAcids.csv`.

    LipiDex skips any line containing "Name" — the header — and ignores the Enabled column
    here, so disabled acids still take part in purity.
    """
    out: list[FattyAcid] = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or "Name" in line:
            continue
        split = line.split(",")
        out.append(FattyAcid(split[0], split[1], split[2], split[3] if len(split) > 3 else ""))
    return out


def _median(values: list[float]) -> float:
    """getMedianChainIntensity — even-length falls to the mean of the middle pair."""
    if not values:
        return 0.0
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2 == 1:
        return ordered[middle]
    return (ordered[middle] + ordered[middle - 1]) / 2.0


def _round_half_up(x: float) -> int:
    return int(floor(x + 0.5))


@dataclass(slots=True)
class PeakPurity:
    """One purity calculation, over one sample spectrum and its candidate identifications."""

    sample_mz: list[float]
    sample_intensity: list[float]
    mz_tol: float = 0.01
    lipids: list[LibrarySpectrum] = field(default_factory=list)
    intensities: list[float] = field(default_factory=list)
    purities: list[int] = field(default_factory=list)
    matched_masses: list[float] = field(default_factory=list)
    _top_masses: list[float] = field(default_factory=list)
    _corrections: list[float] = field(default_factory=list)

    def calc_purity(self, top: LibrarySpectrum, isobaric: list[LibrarySpectrum],
                    fa_db: list[FattyAcid]) -> int:
        best = self._purity_intensity(top, fa_db, top_match=True)
        if best > 0.0:
            self.intensities.append(best)
            self.lipids.append(top)
            self.matched_masses.extend(top.mz)

        for candidate in isobaric:
            value = self._purity_intensity(candidate, fa_db, top_match=False)
            if value > 0.0 and self._is_unique(candidate):
                self.intensities.append(value)
                self.lipids.append(candidate)
                if value > MIN_PURITY_FOR_FRAGMENTS:
                    self.matched_masses.extend(candidate.mz)

        if not self.intensities:
            return 0
        total = sum(self.intensities)
        self.purities = [_round_half_up(v / total * 100) for v in self.intensities]
        return self.purities[0]

    def _is_unique(self, candidate: LibrarySpectrum) -> bool:
        """isUniqueLipid — a substring test, so `PC 16:0_18:1` blocks its own sum-composition."""
        return not any(l.name.find(candidate.name) >= 0 for l in self.lipids)

    def _match(self, mass: float) -> float:
        """findTransitionMatch — sample intensity at a library mass, first hit wins."""
        for m, i in zip(self.sample_mz, self.sample_intensity):
            if abs(m - mass) < self.mz_tol:
                return i
        return 0.0

    def _index_of(self, masses: list[float], mass: float) -> int:
        for i, m in enumerate(masses):
            if abs(m - mass) < MASS_EQUALITY:
                return i
        return -1

    def _purity_intensity(self, ls: LibrarySpectrum, fa_db: list[FattyAcid],
                          top_match: bool) -> float:
        if not ls.is_lipidex:
            return 0.0

        # Most intense moiety fragment type. Generic fragments and neutral losses carry no
        # chain information; cardiolipin DG and PUFA fragments are excluded by name.
        transition_type = ""
        max_intensity = 0.0
        for frag_type, intensity in zip(ls.types, ls.intensity):
            if frag_type is None:
                return 0.0          # the NullPointerException the original catches
            if ("_Fragment" not in frag_type and "_Neutral Loss" not in frag_type
                    and "DG Fragment" not in frag_type and "PUFA" not in frag_type):
                if intensity > max_intensity:
                    max_intensity = intensity
                    transition_type = frag_type
        if transition_type == "":
            return 0.0

        moiety = ""
        if " Fragment" in transition_type:
            moiety = transition_type[transition_type.index("_") + 1:].replace(" Fragment", "")
        elif " Neutral Loss" in transition_type:
            moiety = transition_type[transition_type.index("_") + 1:].replace(" Neutral Loss", "")
        lipid_fas = self._parse_fatty_acids(ls, fa_db, moiety) if moiety else []

        matched = self._fa_intensities(ls, lipid_fas, transition_type,
                                       max_intensity * MOIETY_INTENSITY_FLOOR)
        if matched is None:
            return 0.0
        if len(matched) < 3 or top_match:
            return _median(matched)
        return sorted(matched)[-1]

    @staticmethod
    def _parse_fatty_acids(ls: LibrarySpectrum, fa_db: list[FattyAcid],
                           fa_type: str) -> list[FattyAcid]:
        """The chains named in the lipid, as database entries of the given base type.

        `PC 16:0_18:1 [M+H]+;` -> the text between the first and last space -> 16:0, 18:1.
        A chain the database does not carry under this base type simply does not appear, and
        the caller treats a short list as a failed identification.
        """
        name = ls.name
        first, last = name.find(" "), name.rfind(" ")
        if first < 0 or last <= first:
            return []
        return [fa for chain in name[first + 1:last].split("_")
                for fa in fa_db if fa.name == chain and fa.base == fa_type]

    def _fa_intensities(self, ls: LibrarySpectrum, lipid_fas: list[FattyAcid],
                        transition_type: str, min_intensity: float) -> list[float] | None:
        """Sample intensity of every fragment of one moiety type. None if a chain is missing."""
        matched: list[float] = []
        found = 0
        for k, mass in enumerate(ls.mz):
            if ls.types[k] != transition_type:
                continue
            intensity = self._match(mass)
            if intensity <= min_intensity:
                continue
            corrected = self._index_of(self._top_masses, mass)
            if corrected > -1:
                intensity = max(intensity - self._corrections[corrected], 0.0)
            # A lipid carrying the same chain twice claims that fragment twice.
            for _ in range(sum(1 for fa in lipid_fas if fa.name == ls.fatty_acids[k])):
                matched.append(intensity)
                found += 1

        if found != len(lipid_fas):
            return None

        # Record what this candidate has claimed, so a later candidate predicting the same
        # fragment cannot count the same sample peak again.
        for k, mass in enumerate(ls.mz):
            if ls.types[k] != transition_type:
                continue
            correction = _median(matched)
            index = self._index_of(self._top_masses, mass)
            if mass not in self._top_masses:
                self._top_masses.append(mass)
                self._corrections.append(correction)
            else:
                self._corrections[index] = min(self._corrections[index] + correction,
                                               self._match(mass))
        return matched
