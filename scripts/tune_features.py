"""Tune feature-detection thresholds against the identifications, not against a reference.

Compound Discoverer is a comparator, not ground truth, and optimising to match it would bake
its choices in. The identifications are a better target: every MS2 that clears the score
thresholds is an independent statement that something real elutes at that mass and time, so a
good parameter set finds an MS1 feature for as many of them as it can.

Two numbers per parameter set, and they pull against each other:

  **recall**  fraction of confident identifications that have a feature at their mass and time.
              Raising it is the point.
  **features** how many features the run produced in total. Most are never identified, and a
              setting that doubles this to gain a point of recall has mostly bought noise —
              which then costs time and gives the adduct, blank and in-source filters more to
              chew through.

Neither alone is the answer, so the output is the trade-off, sorted, with a suggested pick:
the setting whose extra features per extra identification is still reasonable. Judge it, do not
read the top line.

The recovery pass is deliberately switched off here — it exists to paper over detection misses,
and would hide exactly what this is trying to measure.

    python scripts/tune_features.py --mzml data/mzml/val_Pos/QC_01.mzML \
        --results outputs/validation_Pos/QC_01_Results.csv
"""
from __future__ import annotations

import argparse
import csv
import itertools
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.features import FeatureParams, _set   # noqa: E402

MATCH_PPM = 20.0        # the tolerance the peak finder associates on
MATCH_RT = 15.0         # seconds either side of the MS2


def load_targets(path: Path) -> list[tuple[float, float]]:
    targets = []
    for row in csv.DictReader(path.open(newline="")):
        try:
            if int(row["Dot Product"]) > 500 and int(row["Reverse Dot Product"]) > 700:
                targets.append((float(row["Precursor Mass"]),
                                float(row["Retention Time (min)"]) * 60.0))
        except (ValueError, KeyError):
            continue
    # One target per distinct (mass, time): a lipid fragmented forty times is one thing to find.
    unique = []
    for mz, rt in sorted(targets):
        if not any(abs(m - mz) / mz * 1e6 < MATCH_PPM and abs(r - rt) < MATCH_RT
                   for m, r in unique):
            unique.append((mz, rt))
    return unique


def load_ms1(mzml: Path):
    from pyopenms import MSExperiment, MzMLFile
    experiment = MSExperiment()
    MzMLFile().load(str(mzml), experiment)
    ms1 = MSExperiment()
    for spectrum in experiment:
        if spectrum.getMSLevel() == 1:
            ms1.addSpectrum(spectrum)
    ms1.sortSpectra(True)
    return ms1


def detect(ms1, params: FeatureParams):
    from pyopenms import (ElutionPeakDetection, FeatureFindingMetabo, FeatureMap,
                          MassTraceDetection)
    mtd = MassTraceDetection()
    p = mtd.getParameters()
    _set(p, {"mass_error_ppm": float(params.mass_error_ppm),
             "noise_threshold_int": float(params.noise_threshold),
             "chrom_peak_snr": float(params.chrom_peak_snr),
             "min_trace_length": float(params.min_trace_length),
             "max_trace_length": float(params.max_trace_length)})
    mtd.setParameters(p)
    traces = []
    mtd.run(ms1, traces, 0)

    epd = ElutionPeakDetection()
    p = epd.getParameters()
    _set(p, {"width_filtering": params.width_filtering,
             "chrom_peak_snr": float(params.chrom_peak_snr),
             "chrom_fwhm": float(params.chrom_fwhm),
             "min_fwhm": float(params.min_fwhm),
             "max_fwhm": float(params.max_fwhm)})
    epd.setParameters(p)
    split = []
    epd.detectPeaks(traces, split)

    ffm = FeatureFindingMetabo()
    p = ffm.getParameters()
    _set(p, {"isotope_filtering_model": params.isotope_filtering_model,
             "remove_single_traces": "true" if params.remove_single_traces else "false",
             "mz_scoring_by_elements": "false", "report_convex_hulls": "false",
             "charge_lower_bound": params.charge_low,
             "charge_upper_bound": params.charge_high})
    ffm.setParameters(p)
    feature_map = FeatureMap()
    chromatograms = []
    ffm.run(split, feature_map, chromatograms)
    return feature_map


def recall(feature_map, targets) -> int:
    index: dict[int, list[tuple[float, float]]] = {}
    for feature in feature_map:
        index.setdefault(int(feature.getMZ()), []).append((feature.getMZ(), feature.getRT()))
    found = 0
    for mz, rt in targets:
        hit = False
        for key in (int(mz) - 1, int(mz), int(mz) + 1):
            for fmz, frt in index.get(key, ()):
                if abs(fmz - mz) / mz * 1e6 < MATCH_PPM and abs(frt - rt) < MATCH_RT:
                    hit = True
                    break
            if hit:
                break
        found += hit
    return found


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mzml", required=True)
    ap.add_argument("--results", required=True)
    ap.add_argument("--noise", nargs="*", type=float, default=[1000, 2000, 5000, 10000])
    ap.add_argument("--snr", nargs="*", type=float, default=[2.0, 3.0])
    ap.add_argument("--min-trace", nargs="*", type=float, default=[2.0, 3.0, 4.0])
    ap.add_argument("--single", nargs="*", default=["true", "false"])
    # ⚠ The one detection setting with no recorded reason. pyOpenMS defaults to
    # `metabolites (5% RMS)` — an isotope-pattern plausibility check — and this pipeline disables
    # it. On one file that carries 18% more features; whether any of them are IDENTIFIED lipids is
    # what this sweep answers, and it is the only way to tell a free specificity gain from a
    # recall cost.
    ap.add_argument("--isotope-model", nargs="*", dest="isotope_models",
                    default=["none", "metabolites (5% RMS)", "metabolites (2% RMS)"])
    args = ap.parse_args()

    targets = load_targets(Path(args.results))
    print(f"{len(targets)} distinct identified peaks to find in {Path(args.mzml).name}\n")
    ms1 = load_ms1(Path(args.mzml))

    rows = []
    grid = list(itertools.product(args.noise, args.snr, args.min_trace, args.single,
                                 args.isotope_models))
    for i, (noise, snr, min_trace, single, iso) in enumerate(grid, 1):
        params = FeatureParams(isotope_filtering_model=iso,
                               noise_threshold=noise, chrom_peak_snr=snr,
                               min_trace_length=min_trace,
                               remove_single_traces=(single == "true"))
        t0 = time.time()
        feature_map = detect(ms1, params)
        found = recall(feature_map, targets)
        rows.append({"noise": noise, "snr": snr, "min_trace": min_trace, "single": single,
                     "isotope_model": iso,
                     "features": feature_map.size(), "found": found,
                     "recall": found / max(len(targets), 1), "seconds": time.time() - t0})
        print(f"  [{i:2d}/{len(grid)}] noise={noise:<6g} snr={snr} min_trace={min_trace} "
              f"single={single:<5} -> {feature_map.size():6d} features, "
              f"recall {100*found/max(len(targets),1):5.1f}%", flush=True)

    print("\n" + "=" * 78)
    print("TRADE-OFF, best recall first")
    print("=" * 78)
    print(f"{'noise':>7} {'snr':>4} {'trace':>6} {'single':>7} {'features':>9} {'recall':>8} "
          f"{'feat/ID':>8} {'sec':>5}")
    for r in sorted(rows, key=lambda r: (-r["recall"], r["features"])):
        print(f"{r['noise']:7g} {r['snr']:4g} {r['min_trace']:6g} {r['single']:>7} "
              f"{r['features']:9d} {100*r['recall']:7.1f}% "
              f"{r['features']/max(r['found'],1):8.1f} {r['seconds']:5.1f}")

    # Pareto front: nothing else gets more recall for fewer features.
    front = [r for r in rows
             if not any(o["recall"] >= r["recall"] and o["features"] < r["features"]
                        for o in rows)]
    front.sort(key=lambda r: r["features"])
    print("\nPARETO FRONT — each row buys recall no cheaper elsewhere")
    previous = None
    for r in front:
        extra = ""
        if previous:
            d_found = r["found"] - previous["found"]
            d_feat = r["features"] - previous["features"]
            extra = (f"   +{d_found} IDs for +{d_feat} features"
                     f" = {d_feat / d_found:.0f} features per ID" if d_found > 0 else "")
        print(f"  noise={r['noise']:<6g} snr={r['snr']} trace={r['min_trace']} "
              f"single={r['single']:<5} {r['features']:6d} features "
              f"recall {100*r['recall']:5.1f}%{extra}")
        previous = r

    print("\nThe last step where features-per-ID is still modest is usually the right pick;"
          "\nbeyond it you are buying noise. Confirm on a held-out file before adopting.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
