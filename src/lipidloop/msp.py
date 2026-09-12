"""Parse .msp spectral libraries.

The libraries are reused unchanged from LipiDex — LipidBlast and the LipiDex HCD sets.
Format is the NIST .msp text convention as LipiDex writes it:

    Name: Cer[ADS] d18:0_24:0 [M-H]-;
    MW: 666.64
    PRECURSORMZ: 666.64
    Comment: Name=... Mass=... Formula=... OptimalPolarity=true Type=LipiDex
    Num Peaks: 11
    618.6195 250 "-48.0211_Neutral Loss_[]"
    337.3481 999 "-118.0499_Sphingoid Neutral Loss_[d18:0]"
    ...

The annotation in quotes is retained: the fragment *type* is what makes a match
interpretable, and the purity calculation needs to know which peaks are fatty-acid
fragments rather than generic neutral losses.

**Two dialects.** LipiDex 1 writes the adduct inside the name (`Cer[ADS] d18:0_24:0 [M-H]-;`)
and flags in a free-text `Comment`. LipiDex 2 writes bare names and puts everything in its own
field — `PRECURSORTYPE`, `IONMODE`, `COMPOUNDCLASS`, `OPTIMALPOLARITY`, `ISLIPIDEX` — with CRLF
line endings. Both are read here, and a v2 entry is normalised to the v1 shape by appending its
`PRECURSORTYPE` to the name, so nothing downstream needs to know which it came from. Reading a
v2 file as though it were v1 does not fail loudly: it silently yields names with no adduct,
which makes every entry look negative-mode.

Reading otherwise reproduces `SpectrumSearcher.java :: readMSP` rather than the format in general,
because the search results depend on which peaks it chose to keep:

  * a peak is kept only if `precursor - mz > 2.0`, so the precursor region never contributes;
  * intensities are then scaled to the base peak at 999 and anything below 5 is discarded.

Peaks are plain lists, not numpy arrays: a library holds hundreds of thousands of spectra of
a few dozen peaks each, and per-spectrum array construction costs more than it returns.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

_ADDUCT = re.compile(r"(\[M[^\]]*\][+-])")

PRECURSOR_CUTOFF = 2.0        # readMSP: precursor region excluded from the library spectrum
MIN_SCALED_INTENSITY = 5.0    # LibrarySpectrum.java :: scaleIntensities
BASE_PEAK = 999.0


def split_annotation(annotation: str) -> tuple[str | None, str]:
    """Transition type and fatty acid from an annotation — `Transition.java :: addType`.

        "O-1_Alkyl Fragment_[10:0]"  ->  ("O-1_Alkyl Fragment", "10:0")

    The type keeps its formula prefix, because that is what LipiDex compares on when it looks
    for the transitions belonging to one moiety. An unannotated peak has no type at all, which
    matters: purity aborts on the first one it meets.
    """
    if not annotation:
        return None, ""
    fatty_acid = ""
    if "[" in annotation and "]" in annotation:
        fatty_acid = annotation[annotation.rindex("[") + 1:annotation.rindex("]")]
    if "_" not in annotation:
        return None, fatty_acid
    return annotation[:annotation.rindex("_")], fatty_acid


# ⚠ `GlcCer` is an overclaim, and the two shipped libraries disagree about it.
#
# Glucosyl- and galactosylceramide are stereoisomers. Reversed-phase chromatography does not
# resolve them and their MS2 spectra are indistinguishable, so no untargeted RP-LC/MS2 method can
# tell which one it saw. The evidence supports "a hexose", which is what `HexCer` says.
#
# It is also a real functional bug rather than a naming preference: LipidBlast writes `GlcCer[NS]`
# (5,418 entries) and LipiDex writes `HexCer[NS]` (4,800) for the same molecules, and both are
# searched together. The retention model fits ONE SURFACE PER CLASS and needs MIN_POINTS = 6
# members, so splitting one class across two names splits its anchors and can leave both halves
# under the minimum — losing the RT model for hexosylceramides entirely.
#
# Normalised on read rather than by editing the .msp files, so the libraries stay byte-identical
# to their published sources and the provenance survives.
_HEXOSE = re.compile(r"^(?:Glc|Gal|Glu|Hex)Cer")


@dataclass(slots=True)
class LibrarySpectrum:
    name: str
    precursor_mz: float
    mz: list[float]
    intensity: list[float]
    annotations: list[str] = field(default_factory=list)
    library: str = ""
    is_lipidex: bool = False
    optimal_polarity: bool = False
    types: list[str | None] = field(default_factory=list)
    fatty_acids: list[str] = field(default_factory=list)

    @property
    def adduct(self) -> str:
        m = _ADDUCT.search(self.name)
        return m.group(1) if m else ""

    @property
    def lipid(self) -> str:
        """Name without the adduct and trailing punctuation, hexosylceramides normalised."""
        raw = _ADDUCT.sub("", self.name).strip().rstrip(";").strip()
        return _HEXOSE.sub("HexCer", raw, count=1)

    @property
    def lipid_class(self) -> str:
        return self.lipid.split(" ", 1)[0] if " " in self.lipid else self.lipid

    @property
    def polarity(self) -> str:
        """LipiDex decides on the name alone: anything without `]+` is negative."""
        return "positive" if "]+" in self.name else "negative"

    def scale_intensities(self) -> None:
        """Base peak to 999, drop everything under 5 — LibrarySpectrum.java lines 144-156."""
        if not self.intensity:
            return
        top = max(self.intensity)
        if top <= 0:
            return
        keep = [k for k, v in enumerate(self.intensity)
                if (v / top) * BASE_PEAK >= MIN_SCALED_INTENSITY]
        self.mz = [self.mz[k] for k in keep]
        self.intensity = [(self.intensity[k] / top) * BASE_PEAK for k in keep]
        self.annotations = [self.annotations[k] for k in keep]
        self.types = [self.types[k] for k in keep]
        self.fatty_acids = [self.fatty_acids[k] for k in keep]


def parse_msp(path: str | Path, library_name: str | None = None) -> Iterator[LibrarySpectrum]:
    """Stream spectra from an .msp file, applying LipiDex's peak filter and scaling.

    Streams rather than loading: LipiDex_HCD_Formic.msp is 47 MB and LipidBlast 23 MB, and a
    full search may load several at once.
    """
    path = Path(path)
    lib = library_name or path.name
    name = ""
    precursor = 0.0
    entry: LibrarySpectrum | None = None
    is_lipidex = False
    optimal_polarity = False
    peak_start = False
    precursor_type = ""     # v2 only: the adduct, which v1 keeps inside the name

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")   # v2 files are CRLF

            if "Name:" in line or "NAME:" in line:
                if entry is not None and precursor > 0.0:
                    entry.scale_intensities()
                    yield entry
                entry = None
                precursor = 0.0
                is_lipidex = False
                optimal_polarity = False
                peak_start = False
                precursor_type = ""
                # LipiDex drops one character after the colon, not the whitespace run.
                name = line[line.index(":") + 1:][1:]

            if "OptimalPolarity=true" in line or line.upper().startswith("OPTIMALPOLARITY: TRUE"):
                optimal_polarity = True
            if "LipiDex" in line or line.upper().startswith("ISLIPIDEX: TRUE"):
                is_lipidex = True
            if line.upper().startswith("PRECURSORTYPE:"):
                precursor_type = line.split(":", 1)[1].strip()

            if "Num Peaks:" in line:
                peak_start = True
                # Normalise a v2 entry to the v1 shape: the adduct belongs in the name, which
                # is where polarity, the adduct itself and the sum composition are read from.
                full_name = f"{name} {precursor_type}" if precursor_type else name
                entry = LibrarySpectrum(name=full_name, precursor_mz=precursor, mz=[], intensity=[],
                                        annotations=[], library=lib, is_lipidex=is_lipidex,
                                        optimal_polarity=optimal_polarity)

            elif peak_start and "." in line and "Num" not in line and entry is not None:
                split = line.split("\t") if "\t" in line else line.split(" ")
                mz = float(split[0])
                if precursor - mz > PRECURSOR_CUTOFF:
                    annotation = (line[line.index('"') + 1:line.rindex('"')]
                                  if is_lipidex and '"' in line else "")
                    frag_type, fatty_acid = split_annotation(annotation)
                    entry.mz.append(mz)
                    entry.intensity.append(float(split[1]))
                    entry.annotations.append(annotation)
                    entry.types.append(frag_type)
                    entry.fatty_acids.append(fatty_acid)

            if "PRECURSORMZ:" in line:
                precursor = float(line[line.rindex(" ") + 1:])

    if entry is not None and precursor > 0.0:
        entry.scale_intensities()
        yield entry
