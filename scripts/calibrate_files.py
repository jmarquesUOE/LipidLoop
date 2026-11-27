"""Measure and correct each file's mass axis, before the features are linked.

A fixed 10 ppm linking window assumes every injection sits on the same mass axis. On MTBLS5163
they do not: three QC injections sat at +13.5 ppm against a batch median near zero, so their
features could never join the main consensus group. They formed their own — two features instead
of forty-one — leaving a zero in the main row for every one of those injections. 48 of 336 lipids
were split that way; on ST004797, 84 of 701, stranding 1,365 measurements.

Nothing is wrong with the instrument. A ~15 ppm drift across a 50-injection batch is ordinary
QTOF behaviour, and the damage is entirely ours: we link at a fixed tolerance and never measure
per-file mass error.

⚠ **This is not the existing `mass_offset_ppm`.** That is one scalar per POLARITY, applied AFTER
grouping, to re-centre the identification search window. It cannot help here, because a single
number cannot represent both a batch median near zero and one injection at +16: correcting the run
by +16 would push the forty-one well-behaved injections out of the window instead.

**Two references, because they answer different questions and neither is sufficient alone.**

*Relative* — each file against the batch consensus m/z of confidently identified lipids. This is
what puts every injection on a common axis, which is the thing that fixes linking, and it is
plentiful: ~180 calibrants per file rather than ~25.

*Absolute* — the batch as a whole against theoretical library masses. Consensus alone would make
the files agree with each other while leaving them jointly wrong, so the batch-level bias is
measured separately and folded into every file's correction.

⚠ Library names are NOT the reported names. Results carry a canonicalised form (`SM d30:1`) while
the library key is the raw MSP entry (`SM d14:1_16:0`), so a name lookup matches about one
identification in eight. That is enough for one batch-level number and nowhere near enough for a
per-file estimate — which is the other reason the per-file work is done against consensus.

    python scripts/calibrate_files.py <study> --results <Unfiltered_Results.csv> --out <dir>
"""
from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

MATCH_PPM = 30.0          # window for finding a calibrant peak; wider than the drift being measured
RT_WINDOW = 0.15          # minutes either side
MIN_CALIBRANTS = 20       # below this the median is not trustworthy and the file is left alone
MIN_CORRECTION = 2.0      # ppm; below this, correcting is noise-fitting


def library_masses(libdir: Path) -> dict[tuple[str, str], float]:
    """{(lipid name, adduct): theoretical precursor m/z} across the target libraries."""
    from lipidloop.msp import parse_msp
    out: dict[tuple[str, str], float] = {}
    for p in sorted(libdir.glob("*.msp")):
        if "DECOY" in p.name or "HOLDOUT" in p.name:
            continue
        for s in parse_msp(str(p)):
            name = (s.lipid or "").strip()
            adduct = (getattr(s, "adduct", "") or "").strip()
            if name and s.precursor_mz > 0:
                out.setdefault((name, adduct), s.precursor_mz)
                out.setdefault((name, ""), s.precursor_mz)
    return out


def calibrants(results: Path, masses: dict, min_files: int, min_dot: float,
               absolute: bool = False) -> list[tuple[float, float]]:
    """[(reference m/z, retention time)] from confident identifications.

    `absolute=False` uses the observed consensus m/z — many points, relative alignment only.
    `absolute=True` uses the theoretical library mass — few points, but true.
    """
    out = []
    for r in csv.DictReader(results.open(errors="replace")):
        name = (r.get("Identification") or "").strip()
        if not name or name.upper().startswith("DECOY"):
            continue
        try:
            dot = float(r.get("Dot Product") or 0)
            seen = int(r.get("Features Found") or 0)
            rt = float(r["Retention Time (min)"])
        except (KeyError, ValueError):
            continue
        if dot < min_dot or seen < min_files:
            continue
        if absolute:
            adduct = (r.get("Adduct") or "").strip()
            mz = masses.get((name, adduct)) or masses.get((name, ""))
        else:
            try:
                mz = float(r["Quant Ion"])
            except (KeyError, ValueError):
                mz = None
        if mz:
            out.append((mz, rt))
    return out


def measure(path: Path, refs: list[tuple[float, float]]) -> tuple[float | None, int, float]:
    """(median ppm error, how many calibrants matched, spread)."""
    from pyopenms import MzMLFile, MSExperiment
    experiment = MSExperiment()
    MzMLFile().load(str(path), experiment)
    ms1 = [s for s in experiment if s.getMSLevel() == 1]
    deltas = []
    for mz, rt in refs:
        best = None
        for s in ms1:
            if abs(s.getRT() / 60 - rt) >= RT_WINDOW:
                continue
            peaks, intensities = s.get_peaks()
            for m, i in zip(peaks, intensities):
                if abs(m - mz) / mz * 1e6 < MATCH_PPM and (best is None or i > best[1]):
                    best = (m, i)
        if best:
            # `best[0]` comes from pyOpenMS's `get_peaks()`, a numpy array -- a bare numpy float64
            # here would poison everything downstream that touches this delta (the median, the
            # correction-floor comparison's bool, the spread tuple), none of which is JSON
            # serializable. One cast at the source, rather than one at every consumer.
            deltas.append(float(1e6 * (best[0] - mz) / mz))
    if len(deltas) < MIN_CALIBRANTS:
        return None, len(deltas), 0.0
    return float(statistics.median(deltas)), len(deltas), float(statistics.pstdev(deltas))


MANIFEST = "calibration.json"


def already_calibrated(directory: Path) -> dict:
    """{filename: applied ppm} for a directory that has been calibrated before, else {}."""
    path = Path(directory) / MANIFEST
    if not path.exists():
        return {}
    try:
        return json.loads(path.read_text()).get("applied_ppm", {})
    except (OSError, json.JSONDecodeError):
        return {}


def apply(source: Path, target: Path, ppm: float) -> None:
    """Write a copy with every m/z shifted by -ppm, MS1, MS2 and precursors alike."""
    from pyopenms import MzMLFile, MSExperiment
    experiment = MSExperiment()
    MzMLFile().load(str(source), experiment)
    factor = 1.0 - ppm * 1e-6
    out = MSExperiment()
    for s in experiment:
        mz, intensity = s.get_peaks()
        s.set_peaks((mz * factor, intensity))
        # ⚠ Precursors too. Correcting the spectra and leaving the precursor behind would shift
        # every MS2 away from the peak it was taken from, breaking identification rather than
        # fixing linking.
        pres = s.getPrecursors()
        if pres:
            for p in pres:
                p.setMZ(p.getMZ() * factor)
            s.setPrecursors(pres)
        out.addSpectrum(s)
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(".partial")
    MzMLFile().store(str(tmp), out)
    tmp.rename(target)


def calibrate_study(study: Path, results: Path, out: Path, libraries: Path = Path("data/libraries"),
                    min_dot: float = 900.0, min_files: int = 10, measure_only: bool = False,
                    polarities: tuple = ("Pos", "Neg"), say=print) -> dict:
    """Measure per-file mass drift against `results`' confident identifications and, unless
    `measure_only`, write corrected mzML into `out`.

    Walks `study/<polarity>` for every `polarity` in `polarities` — pass **exactly one** polarity
    (e.g. `("Neg",)`) to calibrate that polarity alone against ITS OWN identifications. Passing
    both when `results` came from only one polarity's search matches Neg files against Pos-derived
    reference masses and gives a wrong correction, silently — confirmed real on `MTBKS222_Waters`
    (combined-invocation Neg spread -12.3 to +5.4 ppm; isolated-Neg-only spread -0.4 to +9.0 ppm,
    a different answer, not just a noisier one). Callers that need a two-polarity study calibrated
    call this once per polarity, both times pointed at the same `out`.

    Returns a result dict — see `main()` for its use as the CLI's own summary — rather than only
    printing, so `run_study.py`'s `--auto-calibrate` orchestration can inspect the outcome
    (whether anything was corrected, the offsets, the spread) without scraping stdout.
    """
    masses = library_masses(Path(libraries))
    refs = calibrants(Path(results), masses, min_files, min_dot)
    abs_refs = calibrants(Path(results), masses, min_files, min_dot, absolute=True)
    say(f"  {len(refs)} relative calibrants (consensus m/z), "
        f"{len(abs_refs)} absolute (library mass)")
    if len(refs) < MIN_CALIBRANTS:
        say("  too few calibrants — refusing to calibrate on this")
        return {"refused": "too few calibrants", "relative_calibrants": len(refs),
                "absolute_calibrants": len(abs_refs)}

    files = sorted(f for p in polarities for f in (study / p).glob("*.mzML"))

    # ⚠ NEVER measure on already-calibrated files, and never calibrate them twice.
    #
    # The first pass exists to measure the mass axis as the instrument left it. Run against
    # corrected files it finds offsets near zero, concludes nothing needs doing, and the real
    # drift becomes invisible — or, if a correction IS applied, it lands on top of the previous
    # one. Neither failure raises anything: the numbers simply come out wrong and confident.
    #
    # This is not hypothetical now that the conversion cache is shared across analyses of a study:
    # calibrated output written into that cache would be picked up by the next first pass as if it
    # were raw.
    prior = {}
    for p in polarities:
        prior = prior or already_calibrated(study / p)
    if prior:
        say(f"  ⚠ these files were already calibrated ({len(prior)} of them, "
            f"median {statistics.median(list(prior.values())):+.1f} ppm applied).")
        say("    Measuring on corrected data would report a drift near zero and hide the real "
            "offset. Point --study at the ORIGINALS.")
        return {"refused": "already calibrated", "prior_ppm": prior}

    # One batch-level absolute bias, measured on a few files and applied to all. Per-file absolute
    # estimates would rest on ~25 points each; the batch estimate rests on all of them together.
    bias = 0.0
    if len(abs_refs) >= MIN_CALIBRANTS:
        sample = files[:: max(1, len(files) // 5)][:5]
        vals = [v for v, n, _ in (measure(f, abs_refs) for f in sample) if v is not None]
        if vals:
            bias = statistics.median(vals)
            say(f"  batch absolute bias {bias:+.1f} ppm (from {len(sample)} files against "
                f"theoretical masses) — folded into every correction")
    else:
        say("  ⚠ too few library-matched identifications for an absolute reference; "
            "files will be aligned to each other but the batch bias stays unknown")
    offsets, per_file, corrected = {}, {}, 0
    for f in files:
        ppm, n, sd = measure(f, refs)
        if ppm is None:
            say(f"  {f.name:<34} only {n} calibrants — left alone")
            per_file[f.name] = {"ppm": None, "n": n, "sd": None, "corrected": False}
            if not measure_only:
                target = out / f.parent.name / f.name
                target.parent.mkdir(parents=True, exist_ok=True)
                if not target.exists():
                    target.symlink_to(f.resolve())
            continue
        ppm += bias
        offsets[f.name] = ppm
        will_correct = abs(ppm) >= MIN_CORRECTION
        mark = "" if will_correct else "  (within noise, not corrected)"
        say(f"  {f.name:<34} {ppm:+6.1f} ppm  sd {sd:4.1f}  n={n}{mark}")
        per_file[f.name] = {"ppm": ppm, "n": n, "sd": sd, "corrected": will_correct}
        if not measure_only and will_correct:
            apply(f, out / f.parent.name / f.name, ppm)
            corrected += 1
        elif not measure_only:
            target = out / f.parent.name / f.name
            target.parent.mkdir(parents=True, exist_ok=True)
            if not target.exists():
                target.symlink_to(f.resolve())
    if offsets and not measure_only:
        # A record beside the output, so a later pass can tell corrected files from raw ones
        # rather than inferring it from a directory name.
        for polarity in polarities:
            d = out / polarity
            if d.exists():
                (d / MANIFEST).write_text(json.dumps(
                    {"source": str(study), "applied_ppm":
                     {k: round(v, 3) for k, v in offsets.items()}}, indent=2))
    spread = None
    if offsets:
        v = list(offsets.values())
        spread = (min(v), max(v))
        say(f"\n  spread across files: {min(v):+.1f} to {max(v):+.1f} ppm "
            f"(range {max(v)-min(v):.1f} ppm)")
        say(f"  {corrected} corrected, {len(files)-corrected} passed through")
    return {"relative_calibrants": len(refs), "absolute_calibrants": len(abs_refs),
           "batch_bias_ppm": bias if len(abs_refs) >= MIN_CALIBRANTS else None,
           "total_files": len(files), "corrected_files": corrected,
           "spread_ppm": spread, "per_file": per_file}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--results", required=True, help="Unfiltered_Results.csv from a first pass")
    ap.add_argument("--out", required=True)
    ap.add_argument("--libraries", default="data/libraries")
    ap.add_argument("--min-dot", type=float, default=900.0)
    ap.add_argument("--min-files", type=int, default=10)
    ap.add_argument("--measure-only", action="store_true")
    args = ap.parse_args()

    result = calibrate_study(Path(args.study), Path(args.results), Path(args.out),
                             Path(args.libraries), args.min_dot, args.min_files,
                             args.measure_only)
    return 1 if result.get("refused") else 0


if __name__ == "__main__":
    raise SystemExit(main())
