"""Compare lipidloop's search against a LipiDex `_Results.csv` for the same .mgf.

This is the validation that matters before any feature detection exists: the .mgf is the same
input LipiDex read, so every difference is ours. Feature detection is not in the loop at all.

    python scripts/validate_qc.py --polarity Pos --file QC_01 [--limit 2000]
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

from lipidloop.purity import read_fatty_acids            # noqa: E402
from lipidloop.search import DEFAULT_LIBRARIES, iter_results, load_libraries  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("LIPIDLOOP_REFERENCE_DIR", "reference/Test_CD-Lipidex"))  # the CD + LipiDex reference export
LIBRARY_DIR = ROOT / "data/libraries"
# Load order is not cosmetic: it breaks ties between entries that score identically, so it has
# to match the order the libraries were ticked in the LipiDex GUI for the reference run.
DEFAULT_ORDER = list(DEFAULT_LIBRARIES)
FATTY_ACIDS = ROOT / "data/lipidex_src/FattyAcids.csv"


def load_reference(path: Path) -> dict[int, dict]:
    with path.open(newline="") as fh:
        return {int(r["MS2 ID"]): r for r in csv.DictReader(fh)}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--polarity", default="Pos")
    ap.add_argument("--file", default="QC_01")
    ap.add_argument("--limit", type=int, default=0, help="stop after N sample spectra")
    ap.add_argument("--show", type=int, default=10, help="how many mismatches to print")
    ap.add_argument("--libraries", nargs="+", default=DEFAULT_ORDER,
                    help="library stems, in load order (which decides tie-breaks)")
    args = ap.parse_args()
    libraries = [LIBRARY_DIR / f"{stem}.msp" for stem in args.libraries]

    mgf = DATA / args.polarity / f"{args.file}.mgf"
    ref_path = DATA / args.polarity / f"{args.file}_Results.csv"
    reference = load_reference(ref_path)

    t0 = time.time()
    index = load_libraries(libraries)
    print(f"library: {len(index.spectra)} spectra, "
          f"{index.min_mass:.4f}-{index.max_mass:.4f} Da, {time.time()-t0:.1f}s")

    t0 = time.time()
    fa_db = read_fatty_acids(FATTY_ACIDS)
    ours: dict[int, dict] = {}
    n = 0
    for rows, hits in iter_results(mgf, index, fa_db):
        row = rows[0]
        top = hits[0].dot
        n += 1
        if args.limit and row["MS2 ID"] >= args.limit:
            break
        ours[row["MS2 ID"]] = {"name": row["Identification"],
                               "dot": row["Dot Product"],
                               "rev": row["Reverse Dot Product"],
                               "delta": row["Delta m/z"],
                               "lib": row["Library"],
                               "precursor": row["Precursor Mass"],
                               "rt": row["Retention Time (min)"],
                               "purity": row["Purity"],
                               "components": row["Spectral Components"],
                               "fragments": row["Potential Fragments"],
                               # Candidates whose score ties the winner to within floating-
                               # point noise. Which one LipiDex returned is not a decision
                               # its algorithm made, so a difference here is not an error.
                               # Keyed by name *and* library: a tie can span two libraries.
                               "tied": {(h.library_spectrum.name.strip(),
                                         h.library_spectrum.library) for h in hits
                                        if abs(h.dot - top) <= 1e-9 * max(abs(top), 1.0)}}
    print(f"searched in {time.time()-t0:.1f}s, {len(ours)} identified")

    scope = set(ours) | {k for k in reference if not args.limit or k <= max(ours, default=0)}
    only_ours = sorted(set(ours) - set(reference))
    only_ref = sorted(k for k in scope if k in reference and k not in ours)
    both = sorted(set(ours) & set(reference))

    print(f"\nreference rows in scope: {len([k for k in scope if k in reference])}")
    print(f"agree on which spectra are identified: {len(both)}")
    print(f"only ours: {len(only_ours)}   only LipiDex: {len(only_ref)}")

    same_name = same_dot = same_rev = same_delta = tied = 0
    same_purity = same_components = same_fragments = 0
    purity_mismatches = []
    mismatches = []
    for k in both:
        a, b = ours[k], reference[k]
        ref_row = (b["Identification"].strip(), b["Library"])
        name_ok = (a["name"].strip(), a["lib"]) == ref_row
        dot_ok = a["dot"] == int(b["Dot Product"])
        rev_ok = a["rev"] == int(b["Reverse Dot Product"])
        delta_ok = abs(a["delta"] - float(b["Delta m/z"])) < 1e-4
        same_name += name_ok
        same_dot += dot_ok
        same_rev += rev_ok
        same_delta += delta_ok

        # Purity only means anything where the identification agrees; comparing it across a
        # tie would be comparing the purity of two different lipids.
        if name_ok:
            purity_ok = a["purity"] == int(b["Purity"])
            components_ok = a["components"] == b["Spectral Components"]
            fragments_ok = a["fragments"] == b["Potential Fragments"]
            same_purity += purity_ok
            same_components += components_ok
            same_fragments += fragments_ok
            if not (purity_ok and components_ok and fragments_ok):
                purity_mismatches.append((k, a, b))

        if not name_ok and ref_row in a["tied"]:
            tied += 1
        elif not (name_ok and dot_ok and rev_ok):
            mismatches.append((k, a, b))

    n_both = max(len(both), 1)
    print(f"identification identical: {same_name}/{len(both)} ({100*same_name/n_both:.2f}%)")
    print(f"  + tied on score:        {tied}  "
          f"(-> {100*(same_name+tied)/n_both:.2f}% accounted for)")
    print(f"dot product identical:    {same_dot}/{len(both)} ({100*same_dot/n_both:.2f}%)")
    print(f"reverse dot identical:    {same_rev}/{len(both)} ({100*same_rev/n_both:.2f}%)")
    print(f"delta m/z identical:      {same_delta}/{len(both)} ({100*same_delta/n_both:.2f}%)")

    n_named = max(same_name, 1)
    print(f"\nof the {same_name} rows where the identification agrees:")
    print(f"purity identical:         {same_purity}/{same_name} ({100*same_purity/n_named:.2f}%)")
    print(f"spectral components:      {same_components}/{same_name} "
          f"({100*same_components/n_named:.2f}%)")
    print(f"potential fragments:      {same_fragments}/{same_name} "
          f"({100*same_fragments/n_named:.2f}%)")

    if purity_mismatches:
        print(f"\nfirst {min(args.show, len(purity_mismatches))} purity mismatches:")
        for k, a, b in purity_mismatches[:args.show]:
            print(f"  MS2 {k} {a['name']!r} prec {a['precursor']}")
            print(f"    ours: purity={a['purity']} components={a['components']!r}")
            print(f"          fragments={a['fragments']!r}")
            print(f"    ref:  purity={b['Purity']} components={b['Spectral Components']!r}")
            print(f"          fragments={b['Potential Fragments']!r}")

    if mismatches:
        print(f"\nfirst {min(args.show, len(mismatches))} mismatches:")
        for k, a, b in mismatches[:args.show]:
            print(f"  MS2 {k} prec {a['precursor']}")
            print(f"    ours: {a['name']!r} dot={a['dot']} rev={a['rev']} [{a['lib']}]")
            print(f"    ref:  {b['Identification']!r} dot={b['Dot Product']} "
                  f"rev={b['Reverse Dot Product']} [{b['Library']}]")
    if only_ref:
        print(f"\nfirst {min(args.show, len(only_ref))} LipiDex-only spectra:")
        for k in only_ref[:args.show]:
            b = reference[k]
            print(f"  MS2 {k} prec {b['Precursor Mass']} {b['Identification']!r} "
                  f"dot={b['Dot Product']} [{b['Library']}]")
        print("  libraries:", Counter(reference[k]["Library"] for k in only_ref).most_common())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
