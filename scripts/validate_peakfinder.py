"""Validate the peak finder against a reference `Final_Results.csv`.

Feeds Compound Discoverer's own aligned and unaligned exports in, so feature detection is held
constant and any difference is the peak finder's.

    python scripts/validate_peakfinder.py --polarity Pos
"""
from __future__ import annotations

import os
import argparse
import csv
import sys
import time
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.cd_tables import read_aligned, read_unaligned      # noqa: E402
from lipidloop.adducts import load_adducts                       # noqa: E402
from lipidloop.peakfinder import PeakFinder, write_results        # noqa: E402

DATA = Path(os.environ.get("LIPIDLOOP_REFERENCE_DIR", "reference/Test_CD-Lipidex"))  # the CD + LipiDex reference export


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--polarity", default="Pos")
    ap.add_argument("--out", default="outputs")
    ap.add_argument("--min-features", type=int, default=2,
                    help="compounds required per group; the reference run used 2")
    args = ap.parse_args()

    folder = DATA / args.polarity / f"CD_{args.polarity}"
    aligned = folder / f"Aligned_{args.polarity}.csv"
    unaligned = folder / f"Unaligned_{args.polarity}.csv"
    reference = folder / "Final_Results.csv"

    t0 = time.time()
    groups, samples = read_aligned(aligned, min_feature_count=args.min_features)
    print(f"aligned: {len(groups)} compound groups, {len(samples)} samples, "
          f"{sum(len(c.features) for g in groups for c in g.compounds)} features "
          f"({time.time()-t0:.1f}s)")

    adducts = load_adducts()
    print(f"adducts: {len(adducts)} from the library folders")
    finder = PeakFinder(groups, samples, min_feature_count=args.min_features, adducts=adducts)
    for group in groups:
        group.calc_fwhm()
        group.find_quant_ion()
    avg_fwhm = sum(g.avg_fwhm for g in groups) / max(len(groups), 1)
    matched = read_unaligned(unaligned, groups, avg_fwhm)
    print(f"unaligned: {matched} features matched to measured retention times")

    result_files = {s.file: DATA / args.polarity / f"{s.file.split('.raw')[0]}_Results.csv"
                    for s in samples}
    missing = [str(p) for p in result_files.values() if not p.exists()]
    if missing:
        print("missing result files:", missing)
        return 1

    result = finder.run(result_files, log=print)

    out = Path(args.out)
    out.mkdir(exist_ok=True)
    write_results(result, out / f"Final_Results_{args.polarity}.csv")
    write_results(result, out / f"Unfiltered_Results_{args.polarity}.csv", unfiltered=True)

    # ── compare
    ref_rows = list(csv.DictReader(reference.open(newline="")))
    ref_ided = [r for r in ref_rows if r["Identification"]]
    ours = result.kept
    ours_ided = [g for g in ours if g.identification()[0]]

    print(f"\nrows kept:        ours {len(ours):5d}   reference {len(ref_rows):5d}")
    print(f"rows identified:  ours {len(ours_ided):5d}   reference {len(ref_ided):5d}")

    # Match on quant ion + retention time, which is what identifies a row across the two files.
    def key(mz, rt):
        return (round(float(mz), 2), round(float(rt), 1))

    ref_by_key = {}
    for r in ref_rows:
        ref_by_key.setdefault(key(r["Quant Ion"], r["Retention Time (min)"]), []).append(r)

    hit = miss = 0
    name_same = name_diff = 0
    examples = []
    for g in ours_ided:
        k = key(g.quant_ion, g.retention)
        candidates = ref_by_key.get(k, [])
        if not candidates:
            miss += 1
            continue
        hit += 1
        ours_name = g.identification()[0].strip()
        if any(c["Identification"].strip() == ours_name for c in candidates):
            name_same += 1
        else:
            name_diff += 1
            if len(examples) < 10:
                examples.append((g.retention, g.quant_ion, ours_name,
                                 candidates[0]["Identification"]))

    print(f"\nof our {len(ours_ided)} identified rows:")
    print(f"  matched to a reference row by quant ion + RT: {hit}")
    print(f"  no reference row there:                       {miss}")
    print(f"  identification identical:                     {name_same}")
    print(f"  identification differs:                       {name_diff}")

    ref_keys = {key(r["Quant Ion"], r["Retention Time (min)"]) for r in ref_ided}
    our_keys = {key(g.quant_ion, g.retention) for g in ours_ided}
    print(f"\nreference identified rows we did not identify: {len(ref_keys - our_keys)}")
    print(f"identified rows the reference did not:         {len(our_keys - ref_keys)}")

    if examples:
        print("\nfirst identification differences:")
        for rt, mz, a, b in examples:
            print(f"  RT {rt:7.3f} m/z {mz:10.4f}  ours {a!r}  ref {b!r}")

    print("\nfilter reasons:",
          Counter(g.filter_reason for g in result.compound_groups if not g.keep).most_common())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
