"""Would dropping the survey scan from 120K to 60K merge anything the library can tell apart?

Orbitrap resolving power falls with mass as R(m) = R200 * sqrt(200/m), so a setting quoted at m/z
200 is worth half as much at m/z 800 -- exactly where lipids sit. Two peaks of comparable intensity
are treated here as resolved when their separation reaches the peak width, dm >= m / R(m).

A pair only matters if BOTH of these hold:
  * it is resolved at 120K and not at 60K. A pair too close for either setting is unaffected by the
    change, and one far enough apart for both is safe. Only the band between them is at risk.
  * the two species can actually be present at the same time. Chromatography resolves most isobars
    for free, so the observed-data mode requires them to co-elute. The window must be about one peak
    width -- a generous window is the trap here. At 0.2 min against a 0.063 min FWHM this screen
    flagged five pairs with apparent biases of 250-350%, every one of which turned out to be
    separated by two to three peak widths and to contribute nothing at all. Co-elution within a
    window several times the peak width is not co-elution.

Isotope peaks are included as species in their own right: the M+2 of a lipid overlapping the
monoisotopic peak of a relative with one fewer double bond is the classic lipidomics interference,
and asking whether resolution fixes it is the point of the exercise.

Produces the numbers behind manuscript §5.3.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
from collections import defaultdict

C13 = 1.0033548

def resolving_power(mz: float, at200: float) -> float:
    return at200 * (200.0 / mz) ** 0.5

def peak_width(mz: float, at200: float) -> float:
    return mz / resolving_power(mz, at200)

def parse_msp(path: str):
    """(name, precursor m/z, polarity) per entry."""
    name = mz = None
    with open(path, errors="ignore") as fh:
        for line in fh:
            if line[:5].lower() == "name:":
                name = line.split(":", 1)[1].strip().rstrip(";")
                mz = None
            elif line[:12].upper() == "PRECURSORMZ:":
                try:
                    mz = float(line.split(":", 1)[1].strip())
                except ValueError:
                    mz = None
            elif line.strip().lower().startswith("num peaks") and name and mz:
                adduct = re.search(r"\[M[^\]]*\]([+-])", name)
                yield name, mz, (adduct.group(1) if adduct else "")
                name = mz = None

def species(paths, polarity, isotopes=2):
    """Every library ion, plus its isotope peaks, as (m/z, label)."""
    out, seen = [], set()
    for path in paths:
        for name, mz, sign in parse_msp(path):
            if polarity and sign != polarity:
                continue
            if (name, round(mz, 4)) in seen:
                continue
            seen.add((name, round(mz, 4)))
            out.append((mz, name, 0))
            for k in range(1, isotopes + 1):
                out.append((mz + k * C13, name, k))
    out.sort()
    return out

def at_risk(items, low_res, high_res, same_name_ok=False):
    """Pairs resolved at high_res but not at low_res."""
    hits = []
    for i in range(len(items)):
        mz_i, name_i, iso_i = items[i]
        width_low = peak_width(mz_i, low_res)
        for j in range(i + 1, len(items)):
            mz_j, name_j, iso_j = items[j]
            delta = mz_j - mz_i
            if delta > width_low:
                break
            if delta <= 0:
                continue
            if not same_name_ok and name_i == name_j:
                continue          # a lipid against its own isotope is not an ambiguity
            if delta > peak_width(mz_i, high_res):
                hits.append((mz_i, name_i, iso_i, name_j, iso_j, delta))
    return hits

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--libraries", nargs="+", required=True)
    ap.add_argument("--polarity", default="+", choices=("+", "-", ""))
    ap.add_argument("--low", type=float, default=60000)
    ap.add_argument("--high", type=float, default=120000)
    ap.add_argument("--observed", default="", help="Final_Results.csv: restrict to lipids seen")
    ap.add_argument("--rt-window", type=float, default=0.07,
                    help="min; should be about one peak FWHM, not several")
    args = ap.parse_args()

    paths = [p for pattern in args.libraries for p in sorted(glob.glob(pattern))]
    items = species(paths, args.polarity)
    print(f"library ions (+{2} isotopes each): {len(items):,} from {len(paths)} files, "
          f"polarity {args.polarity!r}")
    for mz in (400, 800, 1200):
        print(f"  at m/z {mz}: peak width {peak_width(mz, args.high) * 1000:5.1f} mDa at "
              f"{args.high/1000:.0f}K, {peak_width(mz, args.low) * 1000:5.1f} mDa at "
              f"{args.low/1000:.0f}K")

    hits = at_risk(items, args.low, args.high)
    print(f"\nPairs resolved at {args.high/1000:.0f}K but NOT at {args.low/1000:.0f}K, "
          f"whole library: {len(hits):,}")

    if args.observed:
        names, rts = set(), defaultdict(list)
        for row in csv.DictReader(open(args.observed)):
            nm = (row.get("Identification") or "").strip()
            if not nm:
                continue
            names.add(nm)
            try:
                rts[nm].append(float(row["Retention Time (min)"]))
            except (KeyError, ValueError):
                pass

        def base(label):
            return re.sub(r"\s*\[M[^\]]*\][+-]\d*\s*;?\s*$", "", label).strip()

        seen_names = {base(n) for n in names}
        obs = [it for it in items if base(it[1]) in seen_names]
        print(f"\nrestricted to {len(seen_names):,} lipids identified here -> {len(obs):,} ions")
        hits_obs = at_risk(obs, args.low, args.high)
        print(f"  pairs in the at-risk band: {len(hits_obs):,}")

        coeluting = []
        for mz, a, ia, b, ib, delta in hits_obs:
            ra, rb = rts.get(base(a), []), rts.get(base(b), [])
            if ra and rb and min(abs(x - y) for x in ra for y in rb) <= args.rt_window:
                coeluting.append((mz, a, ia, b, ib, delta))
        print(f"  ... of which co-elute within {args.rt_window} min: {len(coeluting):,}")
        for mz, a, ia, b, ib, delta in sorted(coeluting, key=lambda h: -h[5])[:25]:
            ta = f"{a} (M+{ia})" if ia else a
            tb = f"{b} (M+{ib})" if ib else b
            print(f"    m/z {mz:9.4f}  {delta*1000:5.1f} mDa  needs R={mz/delta:,.0f}"
                  f"   {ta}  vs  {tb}")

if __name__ == "__main__":
    main()
