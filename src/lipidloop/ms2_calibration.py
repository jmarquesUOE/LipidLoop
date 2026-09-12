"""The MS2 mass axis, measured from calibrant ions that lipid fragment spectra always carry.

`ms2_tol` is an absolute m/z window, so it belongs to the analyzer that produced the fragment
spectra. An Orbitrap or a TOF holds fragments to a few millidaltons; a linear ion trap to a few
tenths of a dalton, and it reads low by tens of millidaltons besides. The value inherited from
LipiDex, 0.01 Da, is right for the first class and useless for the second, and nothing in a
configuration file can know which class the next study comes from. So the run measures it.

Measured on every study run through the pipeline on 2026-09-12 (three files per polarity, absolute
-mass calibrants only): every Orbitrap and every TOF captured 95 to 100 % of true fragment peaks
within 10 mDa with a median offset under 1 mDa (Set 1, both Set 2 arms, Q Exactive Plus and HF,
Exploris 480, Agilent 6546/6545/6530A, Sciex TripleTOF 6600, X500R, ZenoTOF 7600); the Lumos with
ion-trap MS2 captured 13 % within 10 mDa, 97 % within 300 mDa in positive mode and 94 % within
500 mDa in negative, reading low by 25 and 47 mDa. See docs/MS2_TOLERANCE.md for the table.

## Why calibrant ions and not the library search

Deriving the window from a search measures the search with itself and needs a window to start.
Ions whose exact mass is known and which nearly every lipidomics run produces need neither:

    positive   184.0733  phosphocholine head group (PC, SM, LPC)
               264.2686  sphingosine d18:1 backbone (SM, Cer, HexCer)
    negative   255.2330 / 279.2330 / 281.2486 / 303.2330 / 327.2330  16:0, 18:2, 18:1, 20:4, 22:6 carboxylates

Neutral losses (the 59 Da trimethylamine loss, the 141 Da phosphoethanolamine loss) were tried and
rejected as calibrants: their expected m/z inherits the error of the reported precursor, which is
the isolation centre or a picked isotope, and they read 15 to 20 mDa off on instruments whose
absolute-mass ions read within 1 mDa. Cholesterol's 369.3516 was rejected for the same reason on
a ZenoTOF (+20 mDa where 184 read +0.2).

## The boundary this respects

The window DESCRIBES the measurement; it does not decide what counts as a hit. The dot-product
thresholds stay where they are. As with the noise floor (`calibrate.py`), the measurement is made
on every run and reported, and applied only when the configured value cannot belong to the
analyzer that produced the spectra -- so a validated Orbitrap run is bit-for-bit unchanged.
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

CALIBRANTS = {
    "+": {"PC head 184.0733": 184.0733, "sphingosine 264.2686": 264.2686},
    "-": {"FA 16:0 255.2330": 255.2330, "FA 18:2 279.2330": 279.2330, "FA 18:1 281.2486": 281.2486,
          "FA 20:4 303.2330": 303.2330, "FA 22:6 327.2330": 327.2330},
}
TOP_N = 10                 # the calibrant must be among the N most intense peaks of its spectrum
MIN_RELATIVE = 0.02        # and at least this fraction of the base peak
SEARCH_WINDOW = {"orbitrap": 0.05, "tof": 0.10, "ion trap": 0.70, "unknown": 0.70}
# Windows a search would use, from which the smallest one capturing TARGET_CAPTURE is chosen.
LADDER = (0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7)
TARGET_CAPTURE = 0.97
# Where a configured window can plausibly belong, per analyzer class. Inside the band the
# configured value stands; outside it was set for a different analyzer and is replaced.
PLAUSIBLE = {"orbitrap": (0.005, 0.05), "tof": (0.005, 0.10), "ion trap": (0.2, 0.7), "unknown": (0.005, 0.7)}
MIN_PEAKS = 200            # fewer calibrant peaks than this and nothing is derived
# A configured window inside its analyzer's band is still replaced when it captures fewer true
# fragments than this: the Bruker timsTOF Pro is a TOF, and 0.01 Da captures 88 % of its fragment
# peaks against 95-100 % on every other TOF and Orbitrap measured.
KEEP_CAPTURE = 0.90
SAMPLE_FILES = 3

ANALYZER_TERMS = {"MS:1000484": "orbitrap", "MS:1000079": "orbitrap", "MS:1000084": "tof",
                  "MS:1000264": "ion trap", "MS:1000078": "ion trap", "MS:1000083": "ion trap", "MS:1000082": "ion trap"}
# Vendors whose converters declare no analyzer term: their DDA lipidomics instruments are TOFs.
TOF_VENDORS = ("Bruker", "Waters", "Agilent", "Sciex", "SCIEX", "AB SCIEX")


def file_analyzer(path: str | Path) -> str:
    """The analyzer class the mzML header declares; a vendor name as fallback; else `unknown`.

    A Thermo hybrid declares both an Orbitrap and an ion trap and the per-scan filter string
    decides (`scan_analyzer`). Bruker and Waters converters write no analyzer term at all.
    """
    head = open(path, "rb").read(400_000).decode("utf-8", "ignore")
    found = {ANALYZER_TERMS[k] for k in ANALYZER_TERMS if k in head}
    if "orbitrap" in found:
        return "orbitrap"
    if "tof" in found:
        return "tof"
    if "ion trap" in found:
        return "ion trap"
    if any(v in head for v in TOF_VENDORS):
        return "tof"
    return "unknown"


def scan_analyzer(spectrum, default: str) -> str:
    """Thermo writes the analyzer into every scan's filter string: ITMS is the trap, FTMS the Orbitrap."""
    if spectrum.metaValueExists("filter string"):
        fs = str(spectrum.getMetaValue("filter string"))
        if fs.startswith("ITMS"):
            return "ion trap"
        if fs.startswith("FTMS"):
            return "orbitrap"
    return default


@dataclass
class Ms2Calibration:
    """What the calibrant ions said about this run's MS2 axis."""

    analyzer: str = "unknown"
    n_spectra: int = 0
    n_peaks: int = 0
    polarity: str = ""
    offset_da: float = 0.0                       # median observed minus expected
    capture: dict = field(default_factory=dict)  # window (Da) -> fraction captured, after the offset
    window: float = 0.0                          # smallest ladder window capturing TARGET_CAPTURE, in band
    files: list = field(default_factory=list)

    @property
    def usable(self) -> bool:
        return self.n_peaks >= MIN_PEAKS and self.window > 0

    def to_dict(self) -> dict:
        return {"analyzer": self.analyzer, "n_spectra": self.n_spectra, "n_calibrant_peaks": self.n_peaks,
                "polarity": self.polarity, "offset_da": round(self.offset_da, 5),
                "capture_after_offset": {str(k): round(v, 4) for k, v in self.capture.items()},
                "window_da": self.window, "files": self.files}

    def __str__(self) -> str:
        if not self.n_peaks:
            return "MS2 axis: no calibrant ion found — nothing measured"
        caps = ", ".join(f"{int(k * 1000)} mDa {100 * v:.0f}%" for k, v in self.capture.items()
                         if k in (0.01, 0.02, 0.05, 0.3, 0.5))
        return (f"MS2 axis: {self.analyzer}, {self.n_peaks:,} calibrant peaks in {self.n_spectra:,} spectra, "
                f"offset {1000 * self.offset_da:+.1f} mDa; captured after offset: {caps}; "
                f"window for {int(100 * TARGET_CAPTURE)}%: {self.window} Da")


def collect_errors(spectra, default_analyzer: str = "unknown"):
    """(polarity, analyzer, error_da) per accepted calibrant peak, from an iterable of pyOpenMS spectra."""
    out = []
    n = 0
    for sp in spectra:
        if sp.getMSLevel() != 2:
            continue
        n += 1
        mz, it = sp.get_peaks()
        if len(mz) < 3:
            continue
        pol = "-" if sp.getInstrumentSettings().getPolarity() == 2 else "+"
        an = scan_analyzer(sp, default_analyzer)
        w = SEARCH_WINDOW.get(an, 0.7)
        base = float(it.max())
        top = float(np.sort(it)[-TOP_N]) if len(it) >= TOP_N else float(it.min())
        for expected in CALIBRANTS[pol].values():
            if expected < mz[0] - w or expected > mz[-1] + w:   # outside the scan range, with the search margin
                continue
            sel = np.abs(mz - expected) <= w
            if not sel.any():
                continue
            k = np.where(sel)[0][np.argmax(it[sel])]
            if it[k] < top or it[k] < MIN_RELATIVE * base:
                continue
            out.append((pol, an, float(mz[k] - expected)))
    return out, n


def summarise(errors, n_spectra: int, files=()) -> Ms2Calibration:
    """Offset, capture after the offset, and the window; from `collect_errors` output."""
    cal = Ms2Calibration(n_spectra=n_spectra, files=list(files))
    if not errors:
        return cal
    pols = [e[0] for e in errors]
    cal.polarity = max(set(pols), key=pols.count)
    ans = [e[1] for e in errors]
    cal.analyzer = max(set(ans), key=ans.count)
    err = np.array([e[2] for e in errors])
    cal.n_peaks = int(len(err))
    cal.offset_da = float(np.median(err))
    corrected = np.abs(err - cal.offset_da)
    cal.capture = {w: float(np.mean(corrected <= w)) for w in LADDER}
    lo, hi = PLAUSIBLE.get(cal.analyzer, PLAUSIBLE["unknown"])
    window = next((w for w in LADDER if cal.capture[w] >= TARGET_CAPTURE), LADDER[-1])
    cal.window = float(min(max(window, lo), hi))
    return cal


def measure_ms2_axis(mzml_files, sample_files: int = SAMPLE_FILES, roles=None, log=None) -> Ms2Calibration:
    """Measure on the largest non-blank files. Cheap enough to run on every run."""
    from pyopenms import MSExperiment, MzMLFile

    files = [Path(f) for f in mzml_files]
    if roles:
        samples = [f for f in files if roles.get(f.stem, "sample") != "blank"]
        files = samples or files
    files = [f for f in files if not re.search(r"blank", f.name, re.I)] or files
    try:
        files = sorted(files, key=lambda f: f.stat().st_size, reverse=True)[:sample_files]
    except OSError:
        files = files[:sample_files]
    errors, n = [], 0
    for f in files:
        try:
            exp = MSExperiment(); MzMLFile().load(str(f), exp)
        except Exception as e:                                   # noqa: BLE001 - a bad file must not stop the run
            if log:
                log(f"MS2 axis: {f.name} could not be read ({e}); skipped")
            continue
        e, k = collect_errors(exp, file_analyzer(f))
        errors += e; n += k
    return summarise(errors, n, files=[f.name for f in files])


def check_ms2_tolerance(configured: float, cal: Ms2Calibration) -> tuple[float, str, bool]:
    """(window to use if deriving, message, whether the configured window is plausible here).

    Same contract as `calibrate.check_noise_threshold`: the third value lets `auto` keep a
    configured window that belongs to this analyzer and replace one that cannot.
    """
    if not cal.usable:
        return configured, f"MS2 axis: {cal.n_peaks} calibrant peaks, too few to judge — configured ms2_tol {configured} stands", True
    lo, hi = PLAUSIBLE.get(cal.analyzer, PLAUSIBLE["unknown"])
    in_band = lo <= configured <= hi
    cap = cal.capture.get(configured) if configured in cal.capture else None
    enough = cap is None or cap >= KEEP_CAPTURE
    at = f", which captures {100 * cap:.0f}% of true fragments here" if cap is not None else ""
    if in_band and enough:
        msg = f"ms2_tol {configured} kept — consistent with a {cal.analyzer}{at} (measured window {cal.window} Da)"
    elif in_band:
        msg = (f"⚠ ms2_tol {configured} is in the band for a {cal.analyzer} but captures only {100 * cap:.0f}% of true "
               f"fragments here (below {int(100 * KEEP_CAPTURE)}%). Measured window {cal.window} Da after a "
               f"{1000 * cal.offset_da:+.1f} mDa offset.")
    else:
        msg = (f"⚠ ms2_tol {configured} is OUTSIDE the plausible band {lo}-{hi} Da for a {cal.analyzer}{at}. "
               f"Measured window {cal.window} Da after a {1000 * cal.offset_da:+.1f} mDa offset.")
    return cal.window, msg, in_band and enough
