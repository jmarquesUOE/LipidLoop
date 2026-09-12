"""Measure the MS2 mass error of a run from calibrant ions that lipid spectra always carry.

The MS2 window (`ms2_tol`) is an absolute m/z tolerance, so it belongs to the analyzer that
produced the fragment spectra: an Orbitrap holds fragments to a few millidaltons, a linear ion
trap to a few tenths of a dalton. Rather than derive the window from a library search, which would
measure the search with itself, this measures it from ions whose exact mass is known and which
occur in nearly every lipidomics run:

  positive mode  184.0733  phosphocholine head group (PC, SM, LPC)
                 264.2686  sphingosine d18:1 backbone fragment (SM, Cer, HexCer)
                 369.3516  cholesterol minus water (CE, cholesterol)
                 and neutral losses from the isolated precursor: 59.0735 (trimethylamine, PC and SM
                 in ion-trap CID, where the 184 ion falls below the low-mass cutoff), 183.0660
                 (phosphocholine, PC), 141.0191 (phosphoethanolamine, PE), 17.0265 (ammonia, TG and
                 DG ammonium adducts), 162.0528 (hexose, HexCer)
  negative mode  255.2330  16:0 carboxylate      279.2330  18:2      281.2486  18:1
                 283.2643  18:0                  303.2330  20:4      327.2330  22:6
                 241.0119  inositol phosphate minus water (PI)   152.9958  glycerophosphate minus water
                 and the neutral loss 60.0211 (methyl formate from formate adducts of PC and SM)

For every MS2 spectrum, each calibrant whose expected m/z lies inside the scan range is looked for
as the most intense peak within a search window (0.05 Da for Fourier-transform and TOF analyzers,
0.7 Da for ion traps), and accepted only if that peak is among the ten most intense in the spectrum
and at least 2 % of the base peak. The error is observed minus expected. Per file and polarity the
script reports the number of accepted calibrant peaks, the median error (a systematic offset), and
the fraction of calibrant peaks captured by a window of 5, 10, 20, 50, 100, 300, 500 and 700 mDa,
which is what a search window buys: `ms2_tol` 0.01 captures the `cap_10` fraction of true fragments. The analyzer is read per scan from the
Thermo filter string where present (ITMS / FTMS) and otherwise from the file's instrument
configuration.

    .venv/bin/python scripts/ms2_mass_error.py --out errors.csv  label=/dir/with/mzml [label=/dir ...]
    .venv/bin/python scripts/ms2_mass_error.py --out errors.csv --max-files 3 label=/dir ...
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

POS_ABS = {"PC head 184": 184.0733, "sphingosine 264": 264.2686, "cholesterol 369": 369.3516}
POS_NL = {"NL trimethylamine 59": 59.0735, "NL phosphocholine 183": 183.0660, "NL phosphoethanolamine 141": 141.0191,
          "NL ammonia 17": 17.0265, "NL hexose 162": 162.0528}
NEG_ABS = {"FA 16:0 255": 255.2330, "FA 18:2 279": 279.2330, "FA 18:1 281": 281.2486, "FA 18:0 283": 283.2643,
           "FA 20:4 303": 303.2330, "FA 22:6 327": 327.2330, "PI 241": 241.0119, "glycerophosphate 153": 152.9958}
NEG_NL = {"NL methyl formate 60": 60.0211}
TOP_N = 10
MIN_REL = 0.02
SEARCH = {"ion trap": 0.7, "orbitrap": 0.05, "tof": 0.1, "unknown": 0.7}
# The window statistics use only fragment ions with an absolute mass. Neutral-loss calibrants
# inherit the error of the reported precursor m/z (isolation centre, isotope picked) and are kept
# as diagnostics in the per-calibrant file only.
CORE = {"PC head 184", "sphingosine 264", "FA 16:0 255", "FA 18:2 279", "FA 18:1 281", "FA 20:4 303", "FA 22:6 327"}


def analyzer_of(spectrum, default: str) -> str:
    fs = spectrum.getMetaValue("filter string") if spectrum.metaValueExists("filter string") else None
    if fs:
        fs = str(fs)
        if fs.startswith("ITMS"):
            return "ion trap"
        if fs.startswith("FTMS"):
            return "orbitrap"
    return default


ANALYZER_TERMS = {"MS:1000484": "orbitrap", "MS:1000079": "orbitrap", "MS:1000084": "tof", "MS:1000264": "ion trap",
                  "MS:1000078": "ion trap", "MS:1000083": "ion trap", "MS:1000082": "ion trap"}


def file_analyzers(path) -> str:
    """The analyzer class declared in the mzML instrument configuration (first 400 kB of the file).

    A Thermo hybrid declares both an Orbitrap and an ion trap; the per-scan filter string then
    decides (`analyzer_of`). Everything else declares one analyzer class, and a TOF is the common case.
    """
    head = open(path, "rb").read(400_000).decode("utf-8", "ignore")
    found = [ANALYZER_TERMS[k] for k in ANALYZER_TERMS if k in head]
    if not found:
        return "unknown"
    if "orbitrap" in found:
        return "orbitrap"
    if "tof" in found:
        return "tof"
    return "ion trap"


def measure_file(path: Path):
    import pyopenms as oms
    exp = oms.MSExperiment(); oms.MzMLFile().load(str(path), exp)
    default = file_analyzers(path)
    rows = []
    n_ms2 = 0
    for sp in exp:
        if sp.getMSLevel() != 2:
            continue
        n_ms2 += 1
        mz, it = sp.get_peaks()
        if len(mz) < 3:
            continue
        pol = "-" if sp.getInstrumentSettings().getPolarity() == 2 else "+"
        an = analyzer_of(sp, default)
        w = SEARCH.get(an, 0.7)
        base = it.max(); top = np.sort(it)[-TOP_N] if len(it) >= TOP_N else it.min()
        prec = sp.getPrecursors()[0].getMZ() if sp.getPrecursors() else None
        cal = dict(POS_ABS if pol == "+" else NEG_ABS)
        for name, loss in (POS_NL if pol == "+" else NEG_NL).items():
            if prec:
                cal[name] = prec - loss
        for name, expected in cal.items():
            if expected < mz[0] or expected > mz[-1]:
                continue
            sel = np.abs(mz - expected) <= w
            if not sel.any():
                continue
            k = np.where(sel)[0][np.argmax(it[sel])]
            if it[k] < top or it[k] < MIN_REL * base:
                continue
            rows.append((pol, an, name, expected, mz[k] - expected, it[k] / base))
    return rows, n_ms2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("sets", nargs="+", help="label=/path/to/dir (mzML found recursively) or label=/path/file.mzML")
    ap.add_argument("--out", required=True)
    ap.add_argument("--max-files", type=int, default=4, help="per set and polarity; largest files first")
    args = ap.parse_args()
    out = open(args.out, "w", newline=""); w = csv.writer(out)
    cal_out = open(args.out.replace(".csv", "_by_calibrant.csv"), "w", newline=""); wc = csv.writer(cal_out)
    wc.writerow(["set", "file", "polarity", "analyzer", "calibrant", "n", "median_error_mDa", "cap_10", "cap_20", "cap_50", "cap_300", "cap_500"])
    WINDOWS = [5, 10, 20, 50, 100, 300, 500, 700]
    w.writerow(["set", "file", "polarity", "analyzer", "n_ms2", "n_calibrant_peaks", "n_calibrants_used",
                "median_error_mDa", "mad_mDa", "median_error_ppm"] + [f"cap_{x}" for x in WINDOWS])
    for arg in args.sets:
        label, path = arg.split("=", 1); path = Path(path)
        files = [path] if path.is_file() else sorted(path.rglob("*.mzML"), key=lambda f: -f.stat().st_size)
        samples = [f for f in files if not re.search(r"blank|_qc|qc_|pooled_qc|^qc", f.name, re.I)]
        files = samples or files
        # up to max-files per polarity, judged by file name when the folder is split, else by content
        picked = []
        seen = {"+": 0, "-": 0}
        for f in files:
            hint = "-" if re.search(r"neg|_n_|_n\.|/Neg/", str(f), re.I) else ("+" if re.search(r"pos|_p_|_p\.|/Pos/", str(f), re.I) else "?")
            if hint != "?" and seen[hint] >= args.max_files:
                continue
            picked.append(f)
            if hint != "?":
                seen[hint] += 1
            if hint == "?" and len(picked) >= 2 * args.max_files:
                break
        for f in picked:
            try:
                rows, n_ms2 = measure_file(f)
            except Exception as e:
                print(f"{label} {f.name}: FAILED {e}", file=sys.stderr); continue
            for pol in ("+", "-"):
                allsub = [r for r in rows if r[0] == pol]
                sub = [r for r in allsub if r[2] in CORE]
                if not sub:
                    continue
                an = max({r[1] for r in sub}, key=lambda a: sum(1 for r in sub if r[1] == a))
                err = np.array([r[4] for r in sub]) * 1000.0
                ppm = np.array([r[4] / r[3] * 1e6 for r in sub])
                used = len({r[2] for r in sub})
                med = float(np.median(err)); mad = float(np.median(np.abs(err - med)))
                caps = [float(np.mean(np.abs(err) <= x)) for x in WINDOWS]
                for name in sorted({r[2] for r in allsub}):
                    e = np.array([r[4] for r in allsub if r[2] == name]) * 1000.0
                    wc.writerow([label, f.name, pol, an, name, len(e), f"{np.median(e):.2f}"] + [f"{float(np.mean(np.abs(e) <= x)):.3f}" for x in (10, 20, 50, 300, 500)])
                cal_out.flush()
                w.writerow([label, f.name, pol, an, n_ms2, len(sub), used, f"{med:.2f}", f"{mad:.2f}", f"{np.median(ppm):.1f}"]
                           + [f"{c:.4f}" for c in caps])
                out.flush()
                print(f"{label:22s} {f.name[:30]:30s} {pol} {an:9s} ms2={n_ms2:6d} cal={len(sub):6d} offset={med:+7.2f} mDa  "
                      + "  ".join(f"{x}:{100*c:5.1f}%" for x, c in zip(WINDOWS, caps)))
    out.close()


if __name__ == "__main__":
    main()
