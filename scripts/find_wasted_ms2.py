"""Find what the MS2 duty cycle is being wasted on, and write an exclusion list for it.

Answers two questions about an acquisition method, from data already acquired:

  1. Which precursors are being fragmented over and over and never identified?
  2. If those slots went somewhere else, would anything be waiting for them?

The second is the one that decides whether an exclusion list is worth the trouble, and `--simulate`
answers it directly: it detects MS1 features, finds the ones that never received an MS2 at all,
and walks the run cycle by cycle handing each freed slot to the most abundant unfragmented feature
co-eluting at that moment. That is an upper bound — a real instrument re-ranks by intensity and
its dynamic exclusion will not follow this exactly — but it is measured on the actual run rather
than assumed.

    python scripts/find_wasted_ms2.py \
        --mzml data/mzml/val_Pos/QC_0*.mzML data/mzml/val_Pos/2*.mzML \
        --results outputs/validation_Pos --blank data/mzml/val_Pos/Blank_0*.mzML \
        --polarity + --out method/Exclusion_List_Pos.csv \
        --thermo method/Exclusion_List_Pos_Thermo.csv \
        --report results/wasted_ms2_Pos.md --simulate

Three outputs because they have three different readers: `--out` is the annotated list carrying
the evidence per row and is the one to review, `--thermo` is the same rows in Xcalibur's mass-list
import layout, `--report` is the whole analysis in markdown with provenance. Pass the blanks to
`--blank` only — listing them under `--mzml` as well counts a file that is almost nothing but
contamination as a sample, and inflates the duty-cycle share.

Read docs/EXCLUSION_LIST.md before putting anything in an acquisition method. Excluding a
precursor is the one decision in this pipeline that cannot be undone by reprocessing.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.exclusion import (MIN_GRADIENT_FRACTION, MIN_SCANS,   # noqa: E402
                                  duty_cycle, find_wasted, series_breakdown,
                                  write_exclusion_list, write_thermo_list)


def read_ms2(path):
    """`(mz, retention_min)` for every MS2 scan, in acquisition order."""
    from pyopenms import MSExperiment, MzMLFile

    experiment = MSExperiment()
    MzMLFile().load(str(path), experiment)
    out = []
    for spectrum in experiment:
        if spectrum.getMSLevel() != 2:
            continue
        precursors = spectrum.getPrecursors()
        out.append((precursors[0].getMZ() if precursors else 0.0, spectrum.getRT() / 60.0))
    return out


def identified_ids(results_dir, stem) -> set:
    path = Path(results_dir) / f"{stem}_Results.csv"
    if not path.exists():
        return set()
    with path.open(newline="") as fh:
        return {int(row["MS2 ID"]) for row in csv.DictReader(fh)}


def simulate(mzml, contaminant_mz, tolerance=0.01):
    """How many never-fragmented features the freed slots could have reached, in this file."""
    from pyopenms import MSExperiment, MzMLFile

    from lipidloop.features import FeatureParams, detect_features

    experiment = MSExperiment()
    MzMLFile().load(str(mzml), experiment)
    cycles, current = [], None
    for spectrum in experiment:
        if spectrum.getMSLevel() == 1:
            current = {"rt": spectrum.getRT() / 60.0, "mz": []}
            cycles.append(current)
        elif spectrum.getMSLevel() == 2 and current is not None:
            precursors = spectrum.getPrecursors()
            current["mz"].append(precursors[0].getMZ() if precursors else 0.0)

    bad = sorted(contaminant_mz)

    def is_contaminant(mz):
        index = bisect.bisect_left(bad, mz - tolerance)
        return index < len(bad) and bad[index] <= mz + tolerance

    feature_map, _ = detect_features(mzml, FeatureParams())
    features = []
    for feature in feature_map:
        box = feature.getConvexHull().getBoundingBox()
        features.append({"mz": feature.getMZ(), "rt": feature.getRT() / 60.0,
                         "intensity": feature.getIntensity(),
                         "lo": box.minPosition()[0] / 60.0, "hi": box.maxPosition()[0] / 60.0})

    fragmented_at = sorted((mz, cycle["rt"]) for cycle in cycles for mz in cycle["mz"])
    masses = [mz for mz, _ in fragmented_at]

    def was_fragmented(mz, rt, mz_tol=0.02, rt_window=0.25):
        index = bisect.bisect_left(masses, mz - mz_tol)
        while index < len(masses) and masses[index] <= mz + mz_tol:
            if abs(fragmented_at[index][1] - rt) <= rt_window:
                return True
            index += 1
        return False

    missed = [f for f in features if not was_fragmented(f["mz"], f["rt"])]
    missed.sort(key=lambda f: -f["intensity"])

    reached = set()
    for cycle in cycles:
        freed = sum(1 for mz in cycle["mz"] if is_contaminant(mz))
        if not freed:
            continue
        available = [k for k, f in enumerate(missed)
                     if f["lo"] - 0.02 <= cycle["rt"] <= f["hi"] + 0.02 and k not in reached]
        reached.update(available[:freed])
    return features, missed, reached


def write_report(path, *, args, wasted, total, share, simulation, per_file) -> None:
    """The analysis, in full, in a file that outlives the terminal it was printed in.

    Every table the run produced, and enough provenance to say what it was run on — the numbers
    age with the column, the labware and the mobile phase, so a list without its date and its
    input files is a list nobody can safely reuse.
    """
    from datetime import date

    lines = [f"# Wasted MS2 duty cycle — polarity `{args.polarity}`", "",
             f"Generated {date.today().isoformat()} by `scripts/find_wasted_ms2.py`.", "",
             "## What this was run on", "",
             f"- **{len(args.mzml)} files**, {total:,} MS2 scans",
             f"- identifications from `{args.results}`" if args.results else
             "- **no `--results` given** — every precursor counted as unidentified, so the "
             "elution-span test is carrying this alone",
             f"- blanks: {', '.join(Path(b).name for b in args.blank) or 'none'}", "",
             "| file | MS2 | identified |", "|---|---|---|"]
    for stem, count, hits in per_file:
        lines.append(f"| {stem} | {count:,} | {hits:,} |")

    lines += ["", "## Criteria", "",
              "All four must hold. Elution span is the load-bearing one — a lipid elutes as a "
              "peak, background is present at every retention time.", "",
              "1. never identified in **any** file of the batch",
              f"2. fragmented at least **{args.min_scans}** times",
              f"3. spans more than **{100 * args.min_gradient_fraction:.0f}%** of the run",
              "4. also fragmented in a blank" if args.blank else
              "4. blank test **not applied** (no `--blank` given)", "",
              "## Headline", "",
              f"**{len(wasted)} precursors are consuming {sum(w.scans for w in wasted):,} of "
              f"{total:,} MS2 — {100 * share:.1f}% of the duty cycle.**", "",
              "## By series", "",
              "The repeat spacing between members is what names the source.", "",
              "| series | precursors | MS2 | share of duty cycle |", "|---|---|---|---|"]
    for name, count, scans, fraction in series_breakdown(wasted, total):
        lines.append(f"| {name} | {count} | {scans:,} | {100 * fraction:.1f}% |")

    if simulation:
        features, missed, reached, intensities = simulation
        lines += ["", "## Where the freed slots would go", "",
                  f"Simulated on `{Path(args.mzml[0]).name}`: MS1 features are detected, the ones "
                  "that never received an MS2 are found, and the run is walked cycle by cycle "
                  "handing each freed slot to the most abundant unfragmented feature co-eluting "
                  "at that moment.", "",
                  "| | count |", "|---|---|",
                  f"| features detected | {features:,} |",
                  f"| never fragmented | {missed:,} ({100 * missed / max(features, 1):.0f}%) |",
                  f"| reachable with the freed slots | **{reached:,} "
                  f"({100 * reached / max(missed, 1):.0f}%)** |"]
        if intensities:
            lines.append(f"| median intensity of those | {statistics.median(intensities):,.0f} |")
            lines.append(f"| largest | {max(intensities):,.0f} |")
        lines += ["", "**Upper bound.** A real instrument re-ranks by intensity live and applies "
                  "its own dynamic exclusion, so it will not follow this exactly, and some freed "
                  "slots will land on isotopes and adducts already covered. The claim is that the "
                  "budget exists and something worth fragmenting is waiting for it — not that "
                  "identifications rise by this much."]

    lines += ["", "## Every candidate", "",
              "Sorted by duty cycle consumed. This is the table to review before importing "
              "anything: an entry with a plausible lipid formula deserves a manual check.", "",
              "| m/z | MS2 | RT observed (min) | % of run | in blank | series |",
              "|---|---|---|---|---|---|"]
    for w in wasted:
        lines.append(f"| {w.mz:.4f} | {w.scans:,} | {w.rt_lo:.2f}–{w.rt_hi:.2f} | "
                     f"{100 * w.gradient_fraction:.0f}% | {'yes' if w.in_blank else 'no'} | "
                     f"{w.series or ''} |")
    lines += ["", "---", "",
              "See [EXCLUSION_LIST.md](../docs/EXCLUSION_LIST.md) for what these criteria are "
              "chosen against, and where they stop being safe. Excluding a precursor is the one "
              "decision in this pipeline that reprocessing cannot undo.", ""]
    Path(path).write_text("\n".join(lines))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--mzml", nargs="+", required=True)
    ap.add_argument("--results", help="folder of <stem>_Results.csv from a completed run")
    ap.add_argument("--blank", nargs="*", default=[],
                    help="blank .mzML; candidates must also fire here")
    ap.add_argument("--polarity", default="+")
    ap.add_argument("--out", help="annotated CSV, one row per candidate, with the evidence")
    ap.add_argument("--thermo", help="the same rows in Xcalibur mass-list import layout")
    ap.add_argument("--report", help="markdown: every table this run produced, plus provenance")
    ap.add_argument("--simulate", action="store_true",
                    help="detect features and estimate what the freed slots would reach")
    ap.add_argument("--min-scans", type=int, default=MIN_SCANS)
    ap.add_argument("--min-gradient-fraction", type=float, default=MIN_GRADIENT_FRACTION)
    args = ap.parse_args()

    if not args.results:
        print("no --results given: every precursor counts as unidentified, so the elution-span "
              "test is carrying this alone. Treat the output as a shortlist to check by hand.\n")

    scans, identified, total, per_file = [], set(), 0, []
    for path in args.mzml:
        stem = Path(path).stem
        hits = identified_ids(args.results, stem) if args.results else set()
        here = 0
        for number, (mz, rt) in enumerate(read_ms2(path)):
            key = (stem, number)
            scans.append((key, mz, rt))
            if number in hits:
                identified.add(key)
            here += 1
        total += here
        per_file.append((stem, here, len(hits)))
        print(f"  {stem}: {total} MS2 so far, {len(identified)} identified")

    blank_mz = [mz for path in args.blank for mz, _ in read_ms2(path)]
    if args.blank:
        print(f"  {len(args.blank)} blank(s): {len(blank_mz)} MS2")

    wasted = find_wasted(scans, identified, blank_mz or None,
                         min_scans=args.min_scans,
                         min_gradient_fraction=args.min_gradient_fraction)

    share = duty_cycle(wasted, total)
    print(f"\n{len(wasted)} precursors are consuming {sum(w.scans for w in wasted)} of {total} "
          f"MS2 — {100 * share:.1f}% of the duty cycle\n")
    print(f"{'m/z':>10} {'MS2':>6} {'RT range':>14} {'% of run':>9}  {'blank':>5}  series")
    for w in wasted[:30]:
        print(f"{w.mz:10.4f} {w.scans:6d} {w.rt_lo:6.2f}-{w.rt_hi:6.2f} "
              f"{100 * w.gradient_fraction:8.0f}%  {'yes' if w.in_blank else 'no':>5}  {w.series}")
    if len(wasted) > 30:
        print(f"  ... and {len(wasted) - 30} more, all in the output file")

    print("\nby series:")
    for name, count, count_scans, fraction in series_breakdown(wasted, total):
        print(f"  {name:16} {count:4d} precursors, {count_scans:7d} MS2 ({100 * fraction:.1f}%)")

    simulation = None
    if args.simulate:
        print("\nsimulating where the freed slots would go (first file only)")
        features, missed, reached = simulate(args.mzml[0], [w.mz for w in wasted])
        intensities = [missed[k]["intensity"] for k in reached]
        simulation = (len(features), len(missed), len(reached), intensities)
        print(f"  {len(features)} features, {len(missed)} never fragmented "
              f"({100 * len(missed) / max(len(features), 1):.0f}%)")
        print(f"  freed slots could reach {len(reached)} of them "
              f"({100 * len(reached) / max(len(missed), 1):.0f}%)")
        if intensities:
            print(f"  their intensity: median {statistics.median(intensities):,.0f}, "
                  f"max {max(intensities):,.0f}")
        print("  upper bound: a real instrument re-ranks by intensity and its dynamic exclusion "
              "will not follow this exactly.")

    for path, writer in ((args.out, write_exclusion_list), (args.thermo, write_thermo_list)):
        if path:
            Path(path).parent.mkdir(parents=True, exist_ok=True)
            writer(wasted, path, args.polarity)
            print(f"\nwritten to {path}")
    if args.report:
        Path(args.report).parent.mkdir(parents=True, exist_ok=True)
        write_report(args.report, args=args, wasted=wasted, total=total, share=share,
                     simulation=simulation, per_file=per_file)
        print(f"written to {args.report}")
    if args.out or args.thermo:
        print("\nRead docs/EXCLUSION_LIST.md first. Excluding a precursor is the one decision "
              "here that reprocessing cannot undo.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
