"""Convert Bruker timsTOF `.d` (TDF) to mzML, without ProteoWizard and without Windows.

The ProteoWizard build on the conversion machine is 3.0.9987 from 2016, which predates timsTOF
support entirely — `[ReaderFail] don't know how to read`. That left 13 files and the 397 published
species of MTBKS222's PASEF arm unreachable.

They do not need it. alphatims reads TDF by opening the SQLite database directly and decompressing
`analysis.tdf_bin` itself; it is a reimplementation rather than a wrapper, so no Bruker SDK, no
Visual C++ runtime, nothing beyond the venv. Verified on these files: 14,050 frames, 32.4M peaks.

**PASEF is not conventional DDA, and the conversion has to respect that.** A frame is 765 mobility
slices, and one MS2 frame carries fragments from SEVERAL precursors, each occupying its own slice
range. So spectra are built PER PRECURSOR — summing only the scans between that precursor's
`ScanNumBegin` and `ScanNumEnd` — rather than by collapsing whole frames, which would merge
co-isolated precursors back together and throw away the separation PASEF exists to provide.

⚠ **This is a transformation, not a format change, and the methods section should say so.** The
mobility dimension is used to separate precursors and is then summed away; the resulting spectra
are cleaner than conventional DDA at the same isolation width, so this arm is not on equal footing
with the other vendors. The per-precursor mobility is written into the mzML so the dimension is at
least recoverable.

    python scripts/tdf_to_mzml.py <file.d> --out <dir>
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

MS1, PASEF_MS2 = 0, 8


def _global_metadata(source: Path) -> dict:
    """The TDF's own GlobalMetadata table, or {} if it cannot be read."""
    import sqlite3
    try:
        con = sqlite3.connect(f"file:{source / 'analysis.tdf'}?mode=ro", uri=True)
        return dict(con.execute("SELECT Key, Value FROM GlobalMetadata"))
    except Exception:
        return {}


def polarity_of(meta: dict) -> str:
    """`positive` / `negative` from the acquisition method name.

    ⚠ **There is no polarity field in a TDF.** It is carried in `MethodName` —
    `LCMS_timsON_pos_default.m` against `LCMS_timsON_neg.m` — and nowhere else in GlobalMetadata.

    The first version of this converter wrote no polarity at all, so all 13 files came out
    undetermined and the stager, which reads polarity from the mzML precisely so it never has to
    guess from a filename, could place none of them. Guessing from the filename here would be the
    same error this pipeline keeps finding in other people's data.
    """
    method = (meta.get("MethodName") or "").lower()
    if "_pos" in method:
        return "positive"
    if "_neg" in method:
        return "negative"
    return ""


def _set_polarity(spectrum, polarity: str) -> None:
    """Write the scan-polarity cvParam the stager and the pipeline both read."""
    from pyopenms import InstrumentSettings, IonSource

    settings = InstrumentSettings()
    settings.setPolarity(IonSource.Polarity.POSITIVE if polarity == "positive"
                         else IonSource.Polarity.NEGATIVE)
    spectrum.setInstrumentSettings(settings)


MERGE_PPM = 60.0


# ⚠ numpy at module scope: `convert` imports it locally, but `_merge` is called from there AND
# tested directly, so it cannot rely on the caller's import.
import numpy as np  # noqa: E402


def _merge(mz, intensity, ppm: float = MERGE_PPM):
    """Sum a precursor's repeat frames onto ONE mass axis, rather than concatenating them.

    ⚠ This is the difference between a usable spectrum and an unusable one. PASEF re-acquires a
    precursor across consecutive frames, and the same ion lands in a slightly different TOF bin each
    time. Concatenating the peak lists therefore writes the SAME fragment several times at slightly
    different m/z, with the same intensity:

        782.5156 (6909)  782.5067 (6909)  782.4713 (6909)  782.4669 (6909)

    one ion, four copies, spread 0.049 Da. LipiDex's forward dot product adds every unmatched sample
    peak to its denominator (`score.py`), so three spurious copies of every real fragment collapse
    it — measured median forward dot 8 against Thermo's 926 — while the reverse score, which ignores
    unmatched peaks, stayed at a perfect 1000. That combination is what the defect looked like from
    the outside.

    ⚠ 60 ppm, not the 10-20 ppm used everywhere else in this pipeline. The spread being merged is
    frame-to-frame jitter of 45-63 ppm, so a tolerance set to the instrument's resolving power does
    not close it — which is why an earlier duplicate check at 20 ppm reported none and the
    hypothesis was dismissed. The window has to exceed the jitter it is there to remove.
    """
    if not len(mz):
        return mz, intensity
    order = np.argsort(mz)
    mz, intensity = mz[order], intensity[order]
    out_mz, out_in = [], []
    start = 0
    for i in range(1, len(mz) + 1):
        if i < len(mz) and (mz[i] - mz[start]) / mz[start] * 1e6 <= ppm:
            continue
        block_mz, block_in = mz[start:i], intensity[start:i]
        total = float(block_in.sum())
        # Intensity-weighted centroid: the strongest observation of the ion should carry the
        # position, not the arithmetic middle of the jitter.
        out_mz.append(float((block_mz * block_in).sum() / total) if total else float(block_mz[0]))
        out_in.append(total)
        start = i
    return np.array(out_mz), np.array(out_in)


def convert(source: Path, target: Path) -> tuple[int, int]:
    """(MS1 spectra, MS2 spectra) written. Skips work already done."""
    import warnings

    import numpy as np
    from pyopenms import MSExperiment, MSSpectrum, MzMLFile, Precursor

    if target.exists() and target.stat().st_size > 0:
        return 0, 0

    meta = _global_metadata(source)
    polarity = polarity_of(meta)
    if not polarity:
        # Loudly, not silently. A file with no polarity cannot be staged, and a run that writes it
        # anyway produces an arm that is quietly missing half its data.
        raise ValueError(f"no polarity in MethodName={meta.get('MethodName')!r}")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        import alphatims.bruker as ab
        data = ab.TimsTOF(str(source))

    frames = data.frames
    experiment = MSExperiment()
    n_ms1 = n_ms2 = 0

    # ── MS1: one spectrum per frame, summed across the whole mobility range ───────────────────
    for _, frame in frames[frames.MsMsType == MS1].iterrows():
        index = int(frame.Id)
        peaks = data[index]
        if peaks is None or not len(peaks):
            continue
        spectrum = MSSpectrum()
        spectrum.setMSLevel(1)
        _set_polarity(spectrum, polarity)
        spectrum.setRT(float(frame.Time))
        mz = peaks["mz_values"].to_numpy(dtype=float)
        intensity = peaks["intensity_values"].to_numpy(dtype=float)
        order = np.argsort(mz)                      # mzML requires ascending m/z
        spectrum.set_peaks((mz[order], intensity[order]))
        experiment.addSpectrum(spectrum)
        n_ms1 += 1

    # ── MS2: one spectrum per PRECURSOR, over that precursor's mobility slice only ────────────
    # ⚠ Sum a precursor's repeat frames rather than writing one spectrum each.
    #
    # PASEF re-acquires the same precursor across consecutive frames to accumulate signal — 2.2
    # frames each here, and measured on a plasma run they sit a median 0.3 s apart, none more than
    # 10 s. They are the same chromatographic moment, so writing them separately throws the
    # accumulation away: median 10 peaks per MS2, with only half the spectra reaching ten peaks at
    # all. That is too thin for a dot product to mean much.
    fragments = data.fragment_frames
    precursors = data.precursors.set_index("Id")
    times = dict(zip(frames.Id.astype(int), frames.Time.astype(float)))

    collected: dict[int, list] = {}
    for _, row in fragments.iterrows():
        collected.setdefault(int(row.Precursor), []).append(row)

    for pid, rows in collected.items():
        if pid not in precursors.index:
            continue
        info = precursors.loc[pid]
        mz_parts, intensity_parts, energies = [], [], []
        for row in rows:
            peaks = data[int(row.Frame), int(row.ScanNumBegin):int(row.ScanNumEnd)]
            if peaks is None or not len(peaks):
                continue
            mz_parts.append(peaks["mz_values"].to_numpy(dtype=float))
            intensity_parts.append(peaks["intensity_values"].to_numpy(dtype=float))
            energies.append(float(row.CollisionEnergy))
        if not mz_parts:
            continue
        mz, intensity = _merge(np.concatenate(mz_parts), np.concatenate(intensity_parts))

        spectrum = MSSpectrum()
        spectrum.setMSLevel(2)
        _set_polarity(spectrum, polarity)
        spectrum.setRT(float(times.get(int(rows[0].Frame), 0.0)))
        spectrum.set_peaks((mz, intensity))

        precursor = Precursor()
        # Monoisotopic where the instrument resolved one; the isolation centre otherwise, which is
        # what was actually targeted.
        mono = float(info.MonoisotopicMz) if info.MonoisotopicMz == info.MonoisotopicMz else 0.0
        precursor.setMZ(mono or float(rows[0].IsolationMz))
        precursor.setCharge(int(info.Charge) if info.Charge == info.Charge else 0)
        width = float(rows[0].IsolationWidth)
        precursor.setIsolationWindowLowerOffset(width / 2)
        precursor.setIsolationWindowUpperOffset(width / 2)
        spectrum.setPrecursors([precursor])
        # ⚠ Keep the mobility. It is the dimension being summed away, and without it nobody
        # reprocessing this file can tell PASEF output from ordinary DDA — the same gap this
        # deposit's own MAFs leave by publishing no CCS at all.
        spectrum.setMetaValue("PASEF precursor id", pid)
        spectrum.setMetaValue("PASEF scan number", float(info.ScanNumber))
        spectrum.setMetaValue("PASEF frames summed", len(mz_parts))
        if energies:
            spectrum.setMetaValue("PASEF collision energy", sum(energies) / len(energies))
        experiment.addSpectrum(spectrum)
        n_ms2 += 1

    # Sample name and acquisition time, both from GlobalMetadata. The timestamp is what injection
    # order is derived from everywhere else in this project; without it this arm would be the one
    # study with no order, having just fixed that for every other.
    # ⚠ MSExperiment has setDateTime directly; there is no setExperimentalSettings, and calling
    # one raises only when the conversion is otherwise complete — thirteen minutes in.
    if meta.get("AcquisitionDateTime"):
        from pyopenms import DateTime
        stamp = DateTime()
        stamp.set(meta["AcquisitionDateTime"].split("+")[0].replace("T", " ")[:19])
        experiment.setDateTime(stamp)
    # ⚠ Write the instrument identity. Without it the converted mzML carries a bare `instrument
    # model` and no vendor, model or serial, so a table built by reading instrument metadata out of
    # the files cannot say what produced this arm — measured, and it is the only arm in the set
    # with no identity at all. GlobalMetadata has the fields; they were simply not carried across.
    from pyopenms import Instrument, Software
    instrument = Instrument()
    instrument.setName(meta.get("InstrumentName") or "Bruker timsTOF")
    instrument.setVendor("Bruker Daltonics")
    if meta.get("InstrumentSerialNumber"):
        instrument.setMetaValue("instrument serial number", meta["InstrumentSerialNumber"])
    experiment.setInstrument(instrument)

    experiment.sortSpectra(True)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(".partial")
    MzMLFile().store(str(tmp), experiment)
    tmp.rename(target)
    return n_ms1, n_ms2


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("sources", nargs="+")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    out = Path(args.out)
    done = 0
    for s in args.sources:
        source = Path(s)
        target = out / f"{source.stem}.mzML"
        try:
            a, b = convert(source, target)
        except Exception as exc:                     # noqa: BLE001 — one bad file must not stop the set
            print(f"  FAIL  {source.name}: {type(exc).__name__}: {exc}", flush=True)
            continue
        size = target.stat().st_size / 1e6 if target.exists() else 0
        print(f"  {source.name:<36} MS1 {a:>5}  MS2 {b:>6}  {size:>6.0f} MB", flush=True)
        done += 1
    print(f"  converted {done} of {len(args.sources)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
