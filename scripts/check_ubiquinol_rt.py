"""Is a reported ubiquinol really ubiquinol, or the ubiquinone M+2 isotope?

The one check `Ubiquinones_Positive.msp` cannot make for itself, and the reason that library is
documented the way it is:

    CoQ10   [M+NH4]+  M+2 isotope   882.7244       (880.7177 + 2 x 1.0033548)
    CoQ10H2 [M+NH4]+  monoisotopic  882.7334
                                    --------
                                    0.0090 Da  =  10.2 ppm

`search.MS1_TOL` is 0.01 Da absolute, so the isotope is inside the window that picks library
candidates, and the MS2 cannot separate them either — both give the quinone head group 197.0808
as base peak, because it is the same head group. Nothing in the peak finder stops this: its
isotope rule fires at M+1 spacing only, and only on a peak carrying NO identification
(known limitations of the peak finder: isotope promotion, and an isotope rule applied only to unidentified peaks), which an identified ubiquinol row is not.

**They separate on retention and on nothing else.** Measured on an in-house Q Exactive Plus study: ubiquinone apexes at
21.07 and its M+2 isotope with it; ubiquinol apexes at 19.50, about 1.5 min earlier, because the
reduced form carries two hydroxyls where the oxidised form carries two ketones and is therefore
more polar on reversed phase. The isotope contributed ZERO at 19.5 in all three files profiled.

So: a ubiquinol row that sits at the ubiquinone retention is the isotope, whatever its dot
product. Run this on any study that reports one.

    python scripts/check_ubiquinol_rt.py Analysis_*/Results/Pos/Final_Results_Filtered.csv

Exit status is 1 if any ubiquinol row fails, so it can gate a delivery. `--anchors` points at a
retention table other than the shipped one — use it when the LC method is not the facility
gradient, because every number here belongs to that gradient and to no other.
"""
from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
ANCHORS = REPO / "data/libraries/Ubiquinones_RT_anchors.csv"

#: A quinol row must sit at least this far BEFORE the quinone of the same chain length. Measured
#: separation is 1.57 min (21.07 - 19.50); half of it is the tolerance, which is still six times
#: the 0.25 min retention window of either peak.
MIN_SEPARATION = 0.75

OXIDISED = {"CoQ6H2": "CoQ6", "CoQ7H2": "CoQ7", "CoQ8H2": "CoQ8", "CoQ9H2": "CoQ9",
            "CoQ10H2": "CoQ10", "CoQ11H2": "CoQ11"}


def lipid_class(row: dict) -> str:
    """The class as the table reports it, falling back to the first token of the name."""
    name = (row.get("Lipid Class") or "").strip()
    if name:
        return name
    ident = (row.get("Identification") or "").strip()
    return ident.split(" ", 1)[0] if ident else ""


def expected_retention(path: Path) -> dict[str, float]:
    """class -> measured retention, from the anchors table."""
    out: dict[str, float] = {}
    if not path.exists():
        return out
    for row in csv.DictReader(path.open()):
        out.setdefault(row["name"], float(row["retention"]))
    return out


def check(table: Path, anchors: dict[str, float], tolerance: float) -> int:
    rows = list(csv.DictReader(table.open()))
    if not rows:
        print(f"{table}: empty")
        return 0
    quinols = [r for r in rows if lipid_class(r) in OXIDISED]
    quinones = [r for r in rows if lipid_class(r) in set(OXIDISED.values())]
    if not quinols:
        print(f"{table}: no quinol rows — nothing to check "
              f"({len(quinones)} ubiquinone row(s) present)")
        return 0

    failures = 0
    for row in quinols:
        cls = lipid_class(row)
        rt = float(row["Retention Time (min)"])
        partner = OXIDISED[cls]
        # Prefer the ubiquinone measured in THIS run; fall back to the anchor table, which is the
        # facility gradient and may not be this study's.
        observed = [float(r["Retention Time (min)"]) for r in quinones
                    if lipid_class(r) == partner]
        if observed:
            # The NEAREST ubiquinone row is the one to beat: if the quinol is not clearly earlier
            # than even the closest quinone peak, it is that peak's isotope.
            reference, source = min(observed, key=lambda t: abs(t - rt)), "this run"
        elif partner in anchors:
            reference, source = anchors[partner], "anchor table"
        else:
            print(f"  ? {cls} at {rt:.2f} — no {partner} to compare against, in this run or the "
                  f"anchors. UNVERIFIED: do not deliver it as ubiquinol.")
            failures += 1
            continue
        separation = reference - rt
        if separation < tolerance:
            failures += 1
            print(f"  FAIL {cls} at RT {rt:.2f} against {partner} at {reference:.2f} "
                  f"({source}): {separation:+.2f} min, under the {tolerance:.2f} min minimum. "
                  f"This is the {partner} M+2 isotope at 882.7244, not ubiquinol — they are "
                  f"10 ppm apart and MS1_TOL is 0.01 Da.")
        else:
            print(f"  ok   {cls} at RT {rt:.2f}, {separation:.2f} min before {partner} "
                  f"at {reference:.2f} ({source}) — the reduced form elutes earlier, as it must.")
    print(f"{table}: {len(quinols)} quinol row(s), {failures} failing")
    return failures


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("tables", nargs="+", help="Final_Results_Filtered.csv / Unfiltered_Results.csv")
    ap.add_argument("--anchors", default=str(ANCHORS))
    ap.add_argument("--min-separation", type=float, default=MIN_SEPARATION,
                    help=f"minutes the quinol must elute before its quinone (default "
                         f"{MIN_SEPARATION}); raise it for a shallower gradient")
    args = ap.parse_args()

    anchors = expected_retention(Path(args.anchors))
    if not anchors:
        print(f"⚠ no anchors read from {args.anchors} — comparison falls back to this run only")
    failures = sum(check(Path(t), anchors, args.min_separation) for t in args.tables)
    if failures:
        print(f"\n{failures} ubiquinol row(s) could not be distinguished from the ubiquinone M+2 "
              f"isotope. the developer notes (PRENOL_QUINONES, not distributed).")
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
