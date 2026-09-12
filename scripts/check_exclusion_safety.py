"""Check an exclusion list against real lipids before it goes near an acquisition method.

An exclusion entry does not remove one m/z. It removes everything the instrument cannot tell
apart from that m/z, and how much that is depends on a method setting — the exclusion mass
tolerance — not on anything the list file can express. So the check has to happen here.

Two things a lipid can collide with, and they are different failures:

  **within the exclusion tolerance** — the lipid matches the list entry and is never fragmented.
    This is the failure that loses an identification, and the one to act on.

  **within the isolation window** (±0.55 Da on this method) — the lipid is still selected on its
    own m/z and still gets an MS2, but the contaminant is co-isolated into it and the spectrum is
    chimeric. Worth knowing, not a reason to drop the entry: that co-isolation is happening today,
    with or without the list.

Checked against two populations, because they answer different questions:

  **identified in this experiment** — lipids we would demonstrably have lost. Concrete.
  **anywhere in the libraries** — lipids we could lose in some future sample. Conservative, and
    the right one for a method that will run on tissue this batch never contained.

    python scripts/check_exclusion_safety.py \\
        --list method/Exclusion_List_Pos.csv --polarity + \\
        --results outputs/validation_Pos/Final_Results.csv \\
        --report results/exclusion_safety_Pos.md --safe-list method/Exclusion_List_Pos_Safe.csv

`--safe-list` writes the list again with every entry that collides inside the exclusion tolerance
removed. Prefer that to editing by hand.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.msp import parse_msp                  # noqa: E402
from lipidloop.search import DEFAULT_LIBRARIES       # noqa: E402

ISOLATION_HALF_WIDTH = 0.55      # read from the mzML: ±0.55 Da on this method
DEFAULT_PPM = 10.0               # a typical exclusion mass tolerance


def adduct_sign(name: str) -> str:
    """`+` or `-` from a library entry's adduct.

    Read after the closing bracket, not by searching for `]-`: that misses every multiply-charged
    entry, since `[M-2H]2-` ends `]2-` and would be filed as positive. It put doubly-charged
    cardiolipins into the positive-mode collision check the first time round.
    """
    if "]" not in name:
        return "+"
    return "-" if "-" in name[name.rindex("]") + 1:] else "+"


def library_lipids(paths, polarity):
    """`(mz, name)` for every library entry in this polarity, sorted by m/z."""
    want = "+" if polarity.strip() in {"+", "Positive", "positive"} else "-"
    out = []
    for path in paths:
        for spectrum in parse_msp(path):
            name = spectrum.name.strip()
            if adduct_sign(name) == want:
                out.append((spectrum.precursor_mz, name))
    out.sort()
    return out


def identified_lipids(results_csv, polarity):
    """`(quant ion m/z, identification)` for identified rows of the matching polarity."""
    want = "+" if polarity.strip() in {"+", "Positive", "positive"} else "-"
    out = []
    with Path(results_csv).open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("Polarity", "").strip() != want or not row.get("Identification", "").strip():
                continue
            out.append((float(row["Quant Ion"]), row["Identification"].strip()))
    out.sort()
    return out


def collisions(mz, lipids, tolerance):
    """Every lipid within `tolerance` of `mz`."""
    masses = [m for m, _ in lipids]
    index = bisect.bisect_left(masses, mz - tolerance)
    out = []
    while index < len(masses) and masses[index] <= mz + tolerance:
        out.append(lipids[index])
        index += 1
    return out


def read_list(path):
    with Path(path).open(newline="") as fh:
        return list(csv.DictReader(fh))


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--list", required=True, help="Exclusion_List_*.csv to check")
    ap.add_argument("--polarity", default="+")
    ap.add_argument("--libraries", nargs="*", default=None)
    ap.add_argument("--results", help="Final_Results.csv — lipids actually identified here")
    ap.add_argument("--no-protect-fatty-acids", action="store_true",
                    help="do not protect free-fatty-acid masses (see the note in main)")
    ap.add_argument("--ppm", type=float, default=DEFAULT_PPM,
                    help="exclusion mass tolerance in the acquisition method")
    ap.add_argument("--isolation", type=float, default=ISOLATION_HALF_WIDTH,
                    help="half the isolation width, Da")
    ap.add_argument("--report")
    ap.add_argument("--safe-list", help="write the list again without the colliding entries")
    ap.add_argument("--safe-thermo", help="the same, in Xcalibur mass-list import layout")
    args = ap.parse_args()

    root = Path(__file__).resolve().parents[1]
    paths = args.libraries or [root / "data/libraries" / f"{n}.msp" for n in DEFAULT_LIBRARIES]
    paths = [p for p in map(Path, paths) if p.exists()]

    entries = read_list(args.list)
    library = library_lipids(paths, args.polarity)
    identified = identified_lipids(args.results, args.polarity) if args.results else []

    # Free-fatty-acid masses are protected whether or not this run identified any, and the reason
    # is a trap this check walked into on real data. The wasted-MS2 screen removes precursors that
    # are never identified, elute across the run and appear in blanks — and a free fatty acid meets
    # all three, because it does not fragment and its mass is also produced in-source by every
    # phospholipid carrying that acyl chain. On the skin set the generated list contained
    # m/z 255.2327 (FA 16:0, 2,679 MS2), 283.2641 (FA 18:0) and 281.2483 (FA 18:1). Excluding those
    # would have permanently prevented ever fragmenting palmitate, stearate or oleate.
    #
    # Checking against `--results` catches it only where the run already names them, which is why
    # it must not be the only guard: the same list generated before the fatty-acid work, or on a
    # tissue where the series is too sparse to name, would have carried them through. A class that
    # cannot be identified by MS2 looks exactly like contamination to a screen built on "never
    # identified".
    if not args.no_protect_fatty_acids and args.polarity.strip() in {"-", "Negative", "negative"}:
        from lipidloop.fatty_acids import candidate_masses
        protected = sorted((mz, f"FA {c}:{db} (protected class)")
                           for (c, db), mz in candidate_masses().items())
        identified = sorted(identified + protected)
        print(f"protecting {len(protected)} free-fatty-acid masses regardless of identification")
    print(f"{len(entries)} exclusion entries, {len(library):,} library lipids "
          f"({len(paths)} libraries), {len(identified)} identified here")

    rows = []
    for entry in entries:
        mz = float(entry["m/z"])
        window = mz * args.ppm / 1e6
        rows.append({
            "mz": mz, "scans": int(entry["MS2 scans wasted"]), "series": entry["Series"],
            "lib_tol": collisions(mz, library, window),
            "id_tol": collisions(mz, identified, window),
            "lib_iso": collisions(mz, library, args.isolation),
            "id_iso": collisions(mz, identified, args.isolation),
        })

    unsafe = [r for r in rows if r["lib_tol"] or r["id_tol"]]
    lost_ids = [r for r in rows if r["id_tol"]]
    chimeric = [r for r in rows if r["id_iso"] and not r["id_tol"]]

    print(f"\nat ±{args.ppm:.0f} ppm (the exclusion tolerance):")
    print(f"  entries colliding with a lipid identified here: {len(lost_ids)}")
    print(f"  entries colliding with any library lipid:       {len(unsafe)}")
    print(f"at ±{args.isolation} Da (the isolation window):")
    print(f"  entries co-isolating with a lipid identified here: "
          f"{len([r for r in rows if r['id_iso']])}")

    for r in sorted(lost_ids, key=lambda r: -r["scans"]):
        names = ", ".join(n for _, n in r["id_tol"][:3])
        print(f"  ! {r['mz']:10.4f} {r['scans']:6d} MS2 {r['series'] or '':16} -> {names}")

    if args.safe_list or args.safe_thermo:
        keep = [e for e, r in zip(entries, rows) if not (r["lib_tol"] or r["id_tol"])]
        if args.safe_list:
            Path(args.safe_list).parent.mkdir(parents=True, exist_ok=True)
            with Path(args.safe_list).open("w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=entries[0].keys())
                writer.writeheader()
                writer.writerows(keep)
        if args.safe_thermo:
            from lipidloop.exclusion import Wasted, write_thermo_list
            write_thermo_list([Wasted(mz=float(e["m/z"]), scans=int(e["MS2 scans wasted"]),
                                      rt_lo=0.0, rt_hi=0.0, gradient_fraction=1.0,
                                      in_blank=e["In blank"] == "yes", series=e["Series"])
                               for e in keep], args.safe_thermo, args.polarity)
        print(f"\n{len(keep)} of {len(entries)} entries kept "
              f"({len(entries) - len(keep)} dropped) -> "
              f"{', '.join(filter(None, (args.safe_list, args.safe_thermo)))}")

    if args.report:
        write_report(args.report, args=args, rows=rows, entries=entries,
                     library=library, identified=identified, paths=paths)
        print(f"written to {args.report}")
    return 0


def write_report(path, *, args, rows, entries, library, identified, paths) -> None:
    from datetime import date

    lost = [r for r in rows if r["id_tol"]]
    lib_hit = [r for r in rows if r["lib_tol"]]
    iso = [r for r in rows if r["id_iso"]]
    wasted_total = sum(r["scans"] for r in rows)

    out = [f"# Exclusion-list safety check — polarity `{args.polarity}`", "",
           f"Generated {date.today().isoformat()} by `scripts/check_exclusion_safety.py`.", "",
           "An exclusion entry does not remove one m/z. It removes everything the instrument "
           "cannot tell apart from it, and how much that is depends on a method setting — the "
           "exclusion mass tolerance — that the list file has no column for.", "",
           "## What was checked", "",
           f"- **{len(entries)} entries** in `{args.list}`",
           f"- **{len(library):,} library lipids** in this polarity, from "
           f"{', '.join(p.stem for p in paths)}",
           f"- **{len(identified)} lipids identified** in `{args.results}`" if args.results
           else "- no `--results` given: only the library check ran",
           f"- exclusion tolerance **±{args.ppm:.0f} ppm**, isolation window "
           f"**±{args.isolation} Da** (read from the mzML)", "",
           "## Result", "",
           "| check | entries | what it would mean |", "|---|---|---|",
           f"| collides with a lipid **identified here**, within the exclusion tolerance | "
           f"**{len(lost)}** | an identification we would lose |",
           f"| collides with **any library lipid**, within the exclusion tolerance | "
           f"{len(lib_hit)} | an identification we could lose in some other sample |",
           f"| co-isolates with a lipid identified here (±{args.isolation} Da) | {len(iso)} | "
           "chimeric spectrum — already happening, with or without the list |", ""]

    if not lost:
        out += ["**No entry collides with anything identified in this experiment.** The list is "
                "safe to import at this tolerance.", ""]
    else:
        out += ["### Entries that would cost an identification", "",
                "| m/z | MS2 wasted | series | lipid |", "|---|---|---|---|"]
        for r in sorted(lost, key=lambda r: -r["scans"]):
            out.append(f"| {r['mz']:.4f} | {r['scans']:,} | {r['series'] or ''} | "
                       f"{', '.join(n for _, n in r['id_tol'])} |")
        out.append("")

    if lib_hit:
        out += ["### Entries colliding with a library lipid not seen here", "",
                "Conservative: these lipids were not identified in this batch, so nothing is lost "
                "today. They matter if the method will run on a tissue this batch did not "
                "contain.", "",
                "| m/z | MS2 wasted | series | library lipid |", "|---|---|---|---|"]
        for r in sorted(lib_hit, key=lambda r: -r["scans"]):
            if r["id_tol"]:
                continue
            names = ", ".join(n for _, n in r["lib_tol"][:4])
            more = f" (+{len(r['lib_tol']) - 4} more)" if len(r["lib_tol"]) > 4 else ""
            out.append(f"| {r['mz']:.4f} | {r['scans']:,} | {r['series'] or ''} | {names}{more} |")
        out.append("")

    safe = [r for r in rows if not (r["lib_tol"] or r["id_tol"])]
    out += ["## If the colliding entries are dropped", "",
            f"- entries kept: **{len(safe)} of {len(entries)}**",
            f"- duty cycle still recovered: **{sum(r['scans'] for r in safe):,} of "
            f"{wasted_total:,} MS2** "
            f"({100 * sum(r['scans'] for r in safe) / max(wasted_total, 1):.1f}% of what the full "
            "list would recover)", "",
            "## The honest reading of the co-isolation column", "",
            f"{len(iso)} entries sit within ±{args.isolation} Da of a lipid identified here. That "
            "is not a reason to drop them. Those contaminants are already being co-isolated into "
            "those lipids' spectra today — excluding the contaminant as a *precursor* does not "
            "put it into the isolation window, it was there anyway. The number is here because it "
            "measures how crowded this m/z space is, which is worth knowing when reading any "
            "spectrum from this method.", ""]
    Path(path).write_text("\n".join(out))


if __name__ == "__main__":
    raise SystemExit(main())
