"""Read the .mgf that Compound Discoverer exports, the way LipiDex reads it.

Reproduces `SpectrumSearcher.java :: readMGF` and `SampleSpectrum.java`. The filters look
arbitrary written out, and they are, but they decide which peaks reach the dot product and so
they decide the score:

  * `precursor - mz > 1.5` — drops the precursor and its isotope envelope;
  * `mz > min_ms2_mass` — the low-mass cutoff, 61.0 in the LipiDex GUI;
  * `intensity > 1.0` — raw counts, before any scaling;
  * then base peak to 999 and drop anything below 5, at search time.

Two details a rewrite gets wrong by being tidy. The precursor is rounded to three decimals
when stored but the *unrounded* value is what the 1.5 Da filter uses, and spectra are numbered
by counting `END IONS` — including spectra that end up with no usable peaks — because that
number is the `MS2 ID` column downstream.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator

MIN_MS2_MASS = 61.0        # SpectrumSearcherGUI.java: lowMassBox default
PRECURSOR_CUTOFF = 1.5     # readMGF
MIN_RAW_INTENSITY = 1.0    # readMGF
MIN_SCALED_INTENSITY = 5.0  # SampleSpectrum.java :: scaleIntensities
BASE_PEAK = 999.0


def _round3(x: float) -> float:
    """Java's Math.round(x*1000)/1000.0 — half-up, not banker's rounding."""
    from math import floor
    return floor(x * 1000.0 + 0.5) / 1000.0


@dataclass(slots=True)
class SampleSpectrum:
    number: int
    precursor: float
    retention: float
    polarity: str
    mz: list[float] = field(default_factory=list)
    intensity: list[float] = field(default_factory=list)
    scaled: bool = False

    def scale_intensities(self) -> None:
        """Base peak to 999, drop under 5 — SampleSpectrum.java lines 54-66."""
        if self.scaled or not self.intensity:
            self.scaled = True
            return
        top = max(self.intensity)
        if top > 0:
            keep = [k for k, v in enumerate(self.intensity)
                    if (v / top) * BASE_PEAK >= MIN_SCALED_INTENSITY]
            self.mz = [self.mz[k] for k in keep]
            self.intensity = [(self.intensity[k] / top) * BASE_PEAK for k in keep]
        self.scaled = True


def read_mzml_ms2(path: str | Path, min_ms2_mass: float = MIN_MS2_MASS,
                  exclude_mz: "list[float] | None" = None,
                  exclude_tolerance: float = 0.06) -> Iterator[SampleSpectrum]:
    """Same spectra, same filters, straight from the .mzML.

    Lets the pipeline run without an .mgf export in the middle. Numbering counts every MS2
    scan in the file, so `MS2 ID` means the same thing it does when reading an .mgf that was
    exported from the same run — but only if that export also kept every MS2 scan.
    """
    from .calibrate import load_centroided

    experiment = load_centroided(path)

    excluded = list(exclude_mz or ())
    number = 0
    for spectrum in experiment:
        if spectrum.getMSLevel() != 2:
            continue
        precursors = spectrum.getPrecursors()
        if not precursors:
            number += 1
            continue
        precursor = precursors[0].getMZ()
        polarity = "-" if spectrum.getInstrumentSettings().getPolarity() == 2 else "+"
        out = SampleSpectrum(number=number, precursor=_round3(precursor),
                             retention=_round3(spectrum.getRT() / 60.0), polarity=polarity)
        mz_array, intensity_array = spectrum.get_peaks()
        for mz, intensity in zip(mz_array, intensity_array):
            # Instrument artefacts are dropped before scaling, which is the point: every
            # spectrum is normalised to its base peak, so an artefact that IS the base peak
            # scales every real fragment down against it and depresses the dot product.
            if any(abs(mz - bad) <= exclude_tolerance for bad in excluded):
                continue
            if (precursor - mz) > PRECURSOR_CUTOFF and mz > min_ms2_mass \
                    and intensity > MIN_RAW_INTENSITY:
                out.mz.append(float(mz))
                out.intensity.append(float(intensity))
        number += 1
        yield out


def read_mgf(path: str | Path, min_ms2_mass: float = MIN_MS2_MASS) -> Iterator[SampleSpectrum]:
    """Stream spectra from an .mgf, applying LipiDex's peak filters."""
    path = Path(path)
    precursor = 0.0
    retention = 0.0
    polarity = "+"
    spectrum: SampleSpectrum | None = None
    count = 0

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith("#"):
                continue

            if "PEPMASS=" in line:
                body = line[line.index("=") + 1:]
                precursor = float(body[:body.rindex(" ")] if " " in body else body)
            elif "RTINSECONDS" in line:
                body = line[line.index("=") + 1:]
                retention = float(body[:body.rindex(" ")] if " " in body else body) / 60.0
            elif "CHARGE" in line:
                polarity = "-" if "-" in line else "+"
            elif "END IONS" in line:
                count += 1
                if spectrum is not None:
                    yield spectrum
                spectrum = None
                polarity = "+"
                precursor = 0.0
                retention = 0.0
            elif "." in line and "TITLE" not in line:
                mz_str, int_str = line.split(" ")[:2]
                if spectrum is None:
                    spectrum = SampleSpectrum(number=count, precursor=_round3(precursor),
                                              retention=_round3(retention), polarity=polarity)
                mz, inten = float(mz_str), float(int_str)
                if (precursor - mz) > PRECURSOR_CUTOFF and mz > min_ms2_mass \
                        and inten > MIN_RAW_INTENSITY:
                    spectrum.mz.append(mz)
                    spectrum.intensity.append(inten)
