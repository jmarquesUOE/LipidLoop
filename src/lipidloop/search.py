"""Search sample MS2 spectra against the libraries, reproducing LipiDex's spectrum searcher.

`SpectrumSearcher.java :: matchLibrarySpectra` bins the library by precursor mass and scores
every entry within `ms1_tol` of a sample precursor. Two behaviours here are worth naming
because they are not what a fresh design would choose, and results differ if you "fix" them:

  * **Polarity is never checked.** Matching is on precursor mass alone, so a positive-mode
    spectrum can be identified against a negative-mode library entry if the masses coincide.
    The `Optimal Polarity` column reports the library entry's own flag; it does not filter.
  * **Ties break by descending library precursor.** `LibrarySpectrum.compareTo` returns 1 when
    the *other* precursor is larger, which sorts the library high-to-low, and the bin scan then
    visits bins low-to-high. Equal dot products therefore resolve to the heavier entry within a
    bin. It reads like a sign slip in the original, but it is what produced the reference
    results, so it is reproduced rather than corrected.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import floor
from pathlib import Path
from typing import Iterable, Sequence

from . import sn_position
from .mgf import SampleSpectrum, read_mgf
from .msp import LibrarySpectrum, parse_msp
from .purity import FattyAcid, PeakPurity
from .score import INT_WEIGHT, MASS_WEIGHT, dot_product

# The default library set, IN LOAD ORDER — the order breaks ties between entries that score
# identically, and this one reproduces the reference run. Append to it rather than reordering.
# Every library here is LipiDex 1's (MIT) or built in this project. The LipiDex 2 ganglioside
# library was removed on 2026-09-13: LipiDex 2 grants no licence, and the negative-mode ganglioside
# library is now enumerated from the Svennerholm structures instead of from it
# (scripts/build_ganglioside_library.py, manuscript/licensing/).
DEFAULT_LIBRARIES = ("LipidBlast_Formic", "LipiDex_HCD_Formic", "LipiDex_HCD_Hydroxy", "Ganglioside_Negative_SumComposition",
                     "Ceramides_Positive")

# Free fatty acids are configured separately (`RunConfig.fatty_acid_libraries`) rather than mixed
# into the set above. They are a different kind of evidence: three quarters of the class produces
# no fragment at all and is named by retention, and the retention surface belongs to one column
# and one gradient. Keeping the list separate is what lets it be swapped for a measured standards
# library without touching the lipid libraries — see docs/FATTY_ACIDS.md.
FATTY_ACID_LIBRARIES = ("FreeFattyAcids_Negative",)

MS1_TOL = 0.01      # SpectrumSearcherGUI.java: ms1TolBox default
MS2_TOL = 0.01      # SpectrumSearcherGUI.java: ms2TolBox default
MAX_RESULTS = 1     # SpectrumSearcherGUI.java: maxSearchResultsSpinner default
MIN_DOT_TO_REPORT = 1.0   # matchLibrarySpectra only records a hit above this


def _round(x: float, places: int = 0) -> float:
    """Java's Math.round: half up, toward positive infinity — not Python's half-to-even."""
    scale = 10.0 ** places
    return floor(x * scale + 0.5) / scale


@dataclass(slots=True)
class Identification:
    spectrum: SampleSpectrum
    library_spectrum: LibrarySpectrum
    delta_ppm: float
    dot: float
    reverse_dot: float


class LibraryIndex:
    """Library spectra binned by precursor mass, as LipiDex bins them."""

    def __init__(self, spectra: Iterable[LibrarySpectrum], ms1_tol: float = MS1_TOL):
        self.ms1_tol = ms1_tol
        self.bin_size = ms1_tol / 10.0
        # Descending precursor, stable — see the class docstring on tie-breaking.
        self.spectra: list[LibrarySpectrum] = sorted(spectra, key=lambda s: -s.precursor_mz)
        for s in self.spectra:
            order = sorted(range(len(s.mz)), key=s.mz.__getitem__)
            s.mz = [s.mz[k] for k in order]
            s.intensity = [s.intensity[k] for k in order]
            s.annotations = [s.annotations[k] for k in order]
            s.types = [s.types[k] for k in order]
            s.fatty_acids = [s.fatty_acids[k] for k in order]

        self.precursors = [s.precursor_mz for s in self.spectra]
        # SpectrumSearcher seeds minMass at 999.0 and only ever lowers it.
        self.min_mass = min([999.0] + self.precursors) if self.precursors else 999.0
        self.max_mass = max([0.0] + self.precursors)

        self.bins: dict[int, list[int]] = {}
        for i, p in enumerate(self.precursors):
            self.bins.setdefault(self._bin(p), []).append(i)

    def _bin(self, mass: float) -> int:
        return int((mass - self.min_mass) / self.bin_size)

    def candidates(self, precursor: float) -> list[int]:
        """Indices within ms1_tol, in LipiDex's scan order (bins ascending)."""
        out: list[int] = []
        for b in range(self._bin(precursor - self.ms1_tol), self._bin(precursor + self.ms1_tol) + 1):
            for i in self.bins.get(b, ()):
                if abs(self.precursors[i] - precursor) < self.ms1_tol:
                    out.append(i)
        return out


def load_libraries(paths: Sequence[str | Path], ms1_tol: float = MS1_TOL,
                   excluded_classes: Sequence[str] = ()) -> LibraryIndex:
    """Load and index the libraries, optionally dropping whole lipid classes.

    Exclusion happens here rather than downstream so an excluded class cannot win a match in
    the first place — it is not in the library at all, and nothing further along has to know.

    The case that motivates it: `LipiDex_HCD_Formic.msp` carries 22,960 `d5TG` entries, a
    combinatorial set of deuterated triacylglycerol standards and about 13% of the whole
    library. Unless d5 standards were actually spiked into the sample, every one of those is a
    mass coincidence waiting to be matched, and a hit is a false positive by construction. If
    you *do* spike SPLASH or similar, take this out — otherwise the standards vanish.
    
    ⚠ **The class test must see past the DECOY_ prefix.** A decoy entry is named
    `DECOY_ d5TG 10:0_10:0_10:0`, so its class parses as the pseudo-class `DECOY_` and
    `"decoy_" != "d5tg"` — meaning `excluded_classes=["d5TG"]` removed 22,960 TARGET entries and
    none of the 22,960 decoys.

    That is exactly backwards for a decoy. Those entries could still win a match and be counted as
    a false positive, while their targets were no longer in the library to win the true one: 17
    d5TG decoy hits, 12.5% of every decoy hit measured, all of them in the FDR numerator against a
    class with no denominator. The fault is general — every decoy in every library parses as
    `DECOY_`, so any excluded class behaved this way.
    """

    excluded = {c.strip().lower() for c in excluded_classes if c.strip()}
    spectra: list[LibrarySpectrum] = []
    dropped = 0
    for p in paths:
        for spectrum in parse_msp(p):
            klass = spectrum.lipid_class.lower()
            if klass.startswith("decoy"):
                # Take the class from the NAME instead, past the prefix, so a decoy is excluded on
                # the same rule as the target it was built from.
                rest = spectrum.lipid.split(None, 1)
                if len(rest) > 1 and rest[0].lower().startswith("decoy"):
                    klass = rest[1].split(None, 1)[0].lower()
            if excluded and klass in excluded:
                dropped += 1
                continue
            spectra.append(spectrum)
    return LibraryIndex(spectra, ms1_tol=ms1_tol)


def search_spectrum(ms2: SampleSpectrum, index: LibraryIndex, ms2_tol: float = MS2_TOL,
                    mass_weight: float = MASS_WEIGHT,
                    int_weight: float = INT_WEIGHT) -> list[Identification]:
    """Score one sample spectrum against the library. Returns hits sorted best first."""
    if not (index.min_mass < ms2.precursor < index.max_mass):
        return []

    ms2.scale_intensities()
    if not ms2.mz:
        return []
    if not ms2.scaled or ms2.mz != sorted(ms2.mz):
        order = sorted(range(len(ms2.mz)), key=ms2.mz.__getitem__)
        ms2.mz = [ms2.mz[k] for k in order]
        ms2.intensity = [ms2.intensity[k] for k in order]

    hits: list[Identification] = []
    for i in index.candidates(ms2.precursor):
        lib = index.spectra[i]
        dot = dot_product(ms2.mz, ms2.intensity, lib.mz, lib.intensity, mz_tol=ms2_tol,
                          reverse=False, mass_weight=mass_weight, int_weight=int_weight,
                          presorted=True)
        if dot <= MIN_DOT_TO_REPORT:
            continue
        rev = dot_product(ms2.mz, ms2.intensity, lib.mz, lib.intensity, mz_tol=ms2_tol,
                          reverse=True, mass_weight=mass_weight, int_weight=int_weight,
                          presorted=True)
        delta = (ms2.precursor - lib.precursor_mz) / lib.precursor_mz * 1e6
        hits.append(Identification(ms2, lib, _round(delta, 4), dot, rev))

    hits.sort(key=lambda h: -h.dot)   # stable: preserves the bin scan order on ties
    return hits


def _jstr(x: float) -> str:
    """Java's Double.toString, near enough for masses in this range."""
    return repr(float(x))


def calc_purity(hits: list[Identification], ms2: SampleSpectrum, fa_db: list[FattyAcid],
                mz_tol: float = MS2_TOL) -> tuple[int, PeakPurity]:
    """Purity of the top identification — `SampleSpectrum.java :: calcPurityAll`.

    The competing candidates are the other identifications at exactly the same library
    precursor and polarity. Anything that is not a LipiDex spectrum in its optimal polarity
    scores 0, and the returned PeakPurity is then empty: no spectral components are written
    for it either.
    """
    calculator = PeakPurity(ms2.mz, ms2.intensity, mz_tol=mz_tol)
    if not hits or hits[0].library_spectrum.name == "":
        return 0, calculator

    top = hits[0].library_spectrum
    isobaric = [h.library_spectrum for h in hits[1:]
                if h.library_spectrum.polarity == top.polarity
                and h.library_spectrum.precursor_mz == top.precursor_mz
                and abs(top.precursor_mz - h.library_spectrum.precursor_mz) < mz_tol]

    if top.is_lipidex and top.optimal_polarity:
        return calculator.calc_purity(top, isobaric, fa_db), calculator
    return 0, calculator


def _potential_fragments(purity: int, calculator: PeakPurity, top: LibrarySpectrum,
                         prefix: str) -> str:
    """The `Potential Fragments` column.

    Above 1 it lists the masses the candidates between them explained; below 1 it falls back to
    the top entry's own predicted fragments. Exactly 1 lands in neither branch and writes
    nothing — a gap in the original, kept because it shows up in real result files.

    The de-duplication is a substring test against the row built so far, so a mass already
    printed anywhere in the row — including inside another number — suppresses it.
    """
    out = ""
    if purity > 1:
        for mass in calculator.matched_masses:
            if _jstr(mass) not in prefix + out:
                out += f"{_jstr(mass)} | "
    elif purity < 1:
        for mass in top.mz:
            out += f"{_jstr(mass)} | "
    return out


def sn_evidence_cell(lib: LibrarySpectrum, ms2: SampleSpectrum,
                     ms2_tol: float = MS2_TOL) -> str:
    """The `sn Evidence` cell for one matched spectrum.

    ⚠ **Every row builder must call this, not `sn_position.evidence` directly.** The production
    pipeline does not go through `iter_results` — `pipeline.py` builds its own row dict inline, and
    `iter_results`/`search_file` serve the tests and the validation harness. The first attempt at
    this column wired only `iter_results`, so on real data the code never ran and the column came
    back empty on all 430 matches *including* the declined cases that should have carried a reason.
    That looked like "nothing qualified" and was diagnosed as a bad hook point; it was a hook in a
    function the pipeline never calls. One helper, called from both, is what keeps that from
    recurring.

    The library entry supplies the name, the class and the per-fragment chain annotations; the
    sample spectrum supplies the intensities. Both are live here and neither survives into the peak
    finder, which reads the search tables back from CSV — so this is the only stage that *can*
    compute it.

    Never blank: a decline carries its reason, which is what distinguishes "class not covered" from
    "fragment not observed" from "the code never ran".

    ⚠ **Both polarities are checked, and either one being positive declines the call.** Matching
    here is on precursor mass alone — this module's own docstring says polarity is never checked —
    so a positive-mode spectrum can match a negative-mode library entry on a mass coincidence. The
    sample's polarity is the ion that was actually fragmented; the library entry's adduct is the
    identity being claimed. The sn-2 rule needs both to be negative, and disagreement between them
    is itself a reason to make no call.
    """
    return str(sn_position.evidence(
        lib.lipid_class, lib.lipid, lib.mz, lib.fatty_acids, ms2.mz, ms2.intensity, tol=ms2_tol,
        polarity="+" if (sn_position._is_positive(ms2.polarity)
                         or sn_position._is_positive(lib.polarity)) else "-"))


RESULT_COLUMNS = ["MS2 ID", "Retention Time (min)", "Rank", "Identification", "Precursor Mass",
                  "Library Mass", "Delta m/z", "Dot Product", "Reverse Dot Product", "Purity",
                  "Spectral Components", "Optimal Polarity", "LipiDex Spectrum", "Library",
                  "Potential Fragments",
                  # Appended, never inserted. `Potential Fragments` de-duplicates masses with a
                  # substring test against the row built so far, and in `pipeline.py` that prefix
                  # is a join over the row dict's values — so a column added ahead of it changes
                  # which masses are suppressed, and the in-source-fragment filter downstream
                  # reads a different set. A new column at the end is inert.
                  "sn Evidence"]


def iter_results(mgf_path: str | Path, index: LibraryIndex, fa_db: list[FattyAcid],
                 ms2_tol: float = MS2_TOL, max_results: int = MAX_RESULTS,
                 min_ms2_mass: float = 61.0):
    """Yield `(rows, hits)` per identified spectrum — rows in LipiDex's `_Results.csv` layout.

    The hits come out alongside so a caller can see the candidates behind a row, which the
    validation needs in order to recognise a tie.

    Only the top-ranked identification carries a purity; lower ranks report 0, as LipiDex does,
    because `calcPurityAll` is only ever called on the first element of the array.
    """
    for ms2 in read_mgf(mgf_path, min_ms2_mass=min_ms2_mass):
        hits = search_spectrum(ms2, index, ms2_tol=ms2_tol)
        if not hits or hits[0].dot <= MIN_DOT_TO_REPORT:
            continue
        purity, calculator = calc_purity(hits, ms2, fa_db, mz_tol=ms2_tol)

        rows: list[dict] = []
        for rank, hit in enumerate(hits[:max_results], start=1):
            lib = hit.library_spectrum
            row_purity = purity if rank == 1 else 0
            components = " / ".join(f"{l.name.replace(';', '')}({p})"
                                    for l, p in zip(calculator.lipids, calculator.purities))
            prefix = (f"{ms2.number},{_jstr(ms2.retention)},{rank},{lib.name},"
                      f"{_jstr(ms2.precursor)},{_jstr(lib.precursor_mz)},"
                      f"{_jstr(hit.delta_ppm)},{int(_round(hit.dot))},"
                      f"{int(_round(hit.reverse_dot))},{row_purity},{components},"
                      f"{str(lib.optimal_polarity).lower()},"
                      f"{str(lib.is_lipidex).lower()},{lib.library},")
            rows.append({
                "MS2 ID": ms2.number,
                "Retention Time (min)": ms2.retention,
                "Rank": rank,
                "Identification": lib.name,
                "Precursor Mass": ms2.precursor,
                "Library Mass": lib.precursor_mz,
                "Delta m/z": hit.delta_ppm,
                "Dot Product": int(_round(hit.dot)),
                "Reverse Dot Product": int(_round(hit.reverse_dot)),
                "Purity": row_purity,
                "Spectral Components": components,
                "Optimal Polarity": str(lib.optimal_polarity).lower(),
                "LipiDex Spectrum": str(lib.is_lipidex).lower(),
                "Library": lib.library,
                "Potential Fragments": _potential_fragments(row_purity, calculator, lib, prefix),
                "sn Evidence": sn_evidence_cell(lib, ms2, ms2_tol=ms2_tol),
            })
        yield rows, hits


def search_file(mgf_path: str | Path, index: LibraryIndex, fa_db: list[FattyAcid],
                ms2_tol: float = MS2_TOL, max_results: int = MAX_RESULTS,
                min_ms2_mass: float = 61.0) -> list[dict]:
    """Search one .mgf and return every result row."""
    return [row for rows, _ in iter_results(mgf_path, index, fa_db, ms2_tol=ms2_tol,
                                            max_results=max_results, min_ms2_mass=min_ms2_mass)
            for row in rows]
