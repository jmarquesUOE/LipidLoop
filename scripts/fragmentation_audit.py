"""Were unidentified features fragmented, and if so why did the spectra not match?

The reported `MS2 Spectra` column cannot answer the first question: `ms2_support()` counts spectra
supporting *the name written*, so it is zero for every unnamed row by construction. Asking it gives
a clean 0% that is a property of the column, not of the data. This script goes to the acquired
spectra instead.

Two guards, because the naive version of this measurement is badly wrong:

  * a decoy control. With ~700,000 MS2 spectra per polarity, coincidence alone can manufacture a
    high hit rate, so the same test runs with precursor mass displaced by a few Da and with
    retention displaced by minutes. If the real rate does not stand well clear of the decoys, the
    measurement means nothing.
  * MS2-identified rows as a positive control. They must score 100%; anything less is a bug in the
    matching here rather than a finding.

Produces the numbers behind manuscript §5.2, Tables 9 and 10.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os

import numpy as np

# Displacements for the decoy control. Chosen to be implausible as real adduct or isotope spacings
# while staying in the same mass region, so the decoy tests coincidence and not a different
# chemistry.
MZ_DECOYS = (+3.111, -3.111, +7.377, -7.377)
RT_DECOYS = (+4.0, -4.0)

# Suppressive artefacts of manuscript §4, per polarity.
ARTEFACTS = {"Pos": (178.295,), "Neg": (97.473,)}


def precursors(mzml_dir: str, skip_blanks: bool = True):
    """Every MS2 precursor (m/z, retention in minutes) plus its peak list, sorted by m/z."""
    from pyopenms import MSExperiment, MzMLFile

    out = []
    for path in sorted(glob.glob(os.path.join(mzml_dir, "*.mzML"))):
        if skip_blanks and "Blank" in os.path.basename(path):
            continue
        run = MSExperiment()
        MzMLFile().load(path, run)
        for spectrum in run:
            if spectrum.getMSLevel() != 2:
                continue
            pre = spectrum.getPrecursors()
            if not pre:
                continue
            mzs, intensities = spectrum.get_peaks()
            if len(mzs) == 0:
                continue
            out.append((pre[0].getMZ(), spectrum.getRT() / 60.0, mzs, intensities))
    out.sort(key=lambda s: s[0])
    return out


def rows(results: str):
    for row in csv.DictReader(open(results)):
        try:
            mz = float(row["Quant Ion"])
            rt = float(row["Retention Time (min)"])
        except (TypeError, ValueError):
            continue
        named = bool((row.get("Identification") or "").strip())
        yield mz, rt, named, (row.get("Identification Source") or "").strip()


def matches(pmz, prt, mz, rt, ppm, rt_tol):
    """Indices of spectra acquired on this feature."""
    tol = mz * ppm / 1e6
    lo, hi = np.searchsorted(pmz, mz - tol), np.searchsorted(pmz, mz + tol)
    return [j for j in range(lo, hi) if abs(prt[j] - rt) <= rt_tol]


def audit(mzml_dir: str, results: str, polarity: str, ppm: float, rt_tol: float):
    spectra = precursors(mzml_dir)
    pmz = np.array([s[0] for s in spectra])
    prt = np.array([s[1] for s in spectra])
    print(f"\n=== {polarity}: {len(spectra):,} MS2 spectra acquired")

    buckets: dict[str, list] = {}
    for mz, rt, named, source in rows(results):
        key = f"identified ({source})" if named else "unidentified"
        buckets.setdefault(key, []).append((mz, rt))

    for key in sorted(buckets):
        feats = buckets[key]
        hit = [bool(matches(pmz, prt, mz, rt, ppm, rt_tol)) for mz, rt in feats]
        real = 100 * sum(hit) / len(feats)
        decoy_mz = np.mean([
            100 * np.mean([bool(matches(pmz, prt, mz + d, rt, ppm, rt_tol)) for mz, rt in feats])
            for d in MZ_DECOYS])
        decoy_rt = np.mean([
            100 * np.mean([bool(matches(pmz, prt, mz, rt + d, ppm, rt_tol)) for mz, rt in feats])
            for d in RT_DECOYS])
        print(f"  {key:22s} n={len(feats):5d}  fragmented {real:5.1f}%  "
              f"never {100 - real:5.1f}%   decoy m/z {decoy_mz:4.1f}%  decoy RT {decoy_rt:4.1f}%")

    # Table 10: why the acquired spectra failed.
    for key in sorted(buckets):
        if key.startswith("identified") and "MS2" not in key:
            continue
        tic, peaks, artefact = [], [], []
        for mz, rt in buckets[key]:
            for j in matches(pmz, prt, mz, rt, ppm, rt_tol):
                _, _, mzs, intensities = spectra[j]
                total = intensities.sum()
                if total <= 0:
                    continue
                tic.append(total)
                peaks.append(len(mzs))
                artefact.append(sum(
                    intensities[k] for k in range(len(mzs))
                    if any(abs(mzs[k] - a) < 0.01 for a in ARTEFACTS[polarity])) / total)
        if tic:
            print(f"  {key:22s} spectra={len(tic):6d}  median TIC {np.median(tic):12,.0f}  "
                  f"median peaks {np.median(peaks):4.0f}  "
                  f"artefact {100 * np.median(artefact):5.1f}% of TIC")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mzml", required=True, help="directory of converted mzML")
    parser.add_argument("--results", required=True, help="Final_Results.csv")
    parser.add_argument("--polarity", required=True, choices=("Pos", "Neg"))
    parser.add_argument("--ppm", type=float, default=10.0)
    parser.add_argument("--rt-tol", type=float, default=0.15, help="minutes")
    args = parser.parse_args()
    audit(args.mzml, args.results, args.polarity, args.ppm, args.rt_tol)


if __name__ == "__main__":
    main()
