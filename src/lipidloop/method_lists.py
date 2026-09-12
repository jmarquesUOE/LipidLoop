"""Method-ready mass lists for the next acquisition, written once per study after both polarities.

`pipeline.py` writes, per polarity, the evidence-bearing `Exclusion_List.csv`, its Xcalibur-layout
twin `Exclusion_List_Thermo.csv`, and a bare `Inclusion_List.csv` (m/z, retention). None of those
is what goes into the instrument as-is, for three reasons this module exists to close:

1. An exclusion entry removes everything within the method's exclusion mass tolerance of it, so a
   real lipid close enough to a contaminant mass would never be fragmented again. The list has to
   be checked against the lipids identified in this run and against every library lipid, and the
   colliding entries dropped, BEFORE import. `scripts/check_exclusion_safety.py` did this by hand;
   here it happens on every run.
2. The inclusion list is two numbers per row; Xcalibur wants the same mass-list layout as the
   exclusion list, with a retention window per entry.
3. Both lists, both polarities, plus the safety report belong in ONE folder with unambiguous
   names, so the person building the method does not have to look inside `Results/<pol>/`.

Output, under `<analysis>/Method_Lists/`:

    Exclusion_List_Pos_Thermo.csv     safe list, Xcalibur mass-list import layout
    Exclusion_List_Neg_Thermo.csv
    Inclusion_List_Pos_Thermo.csv     same layout, Start/End = retention ± `inclusion_window`
    Inclusion_List_Neg_Thermo.csv
    Exclusion_Safety_Pos.md           what was dropped and why, per polarity
    Exclusion_Safety_Neg.md
    README.md                         counts, the tolerances assumed, and the one thing the file
                                      cannot express (the method's own exclusion tolerance)

Nothing here changes any result table.
"""
from __future__ import annotations

import bisect
import csv
import json
from dataclasses import dataclass
from pathlib import Path

from .exclusion import Wasted, write_thermo_list
from .msp import parse_msp

ISOLATION_HALF_WIDTH = 0.55     # Da; ±0.55 on the validated method
DEFAULT_PPM = 10.0              # the exclusion mass tolerance a method typically uses
DEFAULT_INCLUSION_WINDOW = 0.5  # min either side of the observed retention

THERMO_HEADER = ("Mass [m/z],Formula [M],Species,CS [z],Polarity,Start [min],End [min],"
                 "(N)CE,MSX ID,Comment")


def adduct_sign(name: str) -> str:
    """`+` or `-` from a library entry's adduct, read after the closing bracket so that
    `[M-2H]2-` is negative and not mistaken for positive."""
    if "]" not in name:
        return "+"
    return "-" if "-" in name[name.rindex("]") + 1:] else "+"


def library_lipids(paths, polarity: str) -> list[tuple[float, str]]:
    """`(m/z, name)` for every library entry in this polarity, sorted by m/z. Decoy libraries are
    skipped: a decoy is not a lipid a method could lose."""
    want = "+" if polarity in {"+", "Pos", "Positive", "positive"} else "-"
    out = []
    for path in paths:
        path = Path(path)
        if not path.exists() or path.name.upper().startswith("DECOY"):
            continue
        for spectrum in parse_msp(path):
            name = spectrum.name.strip()
            if adduct_sign(name) == want:
                out.append((spectrum.precursor_mz, name))
    out.sort()
    return out


def identified_lipids(results_csv, polarity: str) -> list[tuple[float, str]]:
    """`(quant ion m/z, identification)` for identified, non-decoy rows of this polarity."""
    want = "+" if polarity in {"+", "Pos", "Positive", "positive"} else "-"
    out = []
    with Path(results_csv).open(newline="") as fh:
        for row in csv.DictReader(fh):
            name = row.get("Identification", "").strip()
            if row.get("Polarity", "").strip() != want or not name:
                continue
            if name.upper().startswith("DECOY"):
                continue
            out.append((float(row["Quant Ion"]), name))
    out.sort()
    return out


def protected_fatty_acids() -> list[tuple[float, str]]:
    """Free-fatty-acid `[M-H]-` masses, protected in negative mode whether or not this run named
    any. A free fatty acid does not fragment, elutes wherever its parent phospholipids do and shows
    up in blanks, so the wasted-MS2 screen cannot tell it from contamination; excluding palmitate
    would be permanent."""
    from .fatty_acids import candidate_masses
    return sorted((mz, f"FA {c}:{db} (protected class)") for (c, db), mz in candidate_masses().items())


def collisions(mz: float, lipids: list[tuple[float, str]], tolerance: float) -> list[tuple[float, str]]:
    masses = [m for m, _ in lipids]
    i = bisect.bisect_left(masses, mz - tolerance)
    out = []
    while i < len(masses) and masses[i] <= mz + tolerance:
        out.append(lipids[i])
        i += 1
    return out


@dataclass
class Checked:
    entry: dict
    lib_tol: list
    id_tol: list
    id_iso: list

    @property
    def safe(self) -> bool:
        return not (self.lib_tol or self.id_tol)


def check_exclusion(entries: list[dict], library: list, identified: list,
                    ppm: float = DEFAULT_PPM, isolation: float = ISOLATION_HALF_WIDTH) -> list[Checked]:
    out = []
    for e in entries:
        mz = float(e["m/z"])
        window = mz * ppm / 1e6
        out.append(Checked(e, collisions(mz, library, window), collisions(mz, identified, window),
                           collisions(mz, identified, isolation)))
    return out


def _wasted(e: dict) -> Wasted:
    return Wasted(mz=float(e["m/z"]), scans=int(e["MS2 scans wasted"]), rt_lo=0.0, rt_hi=0.0,
                  gradient_fraction=1.0, in_blank=e.get("In blank", "") == "yes",
                  series=e.get("Series", ""))


def write_inclusion_thermo(rows: list[dict], path, polarity: str, window: float,
                           charge: int = 1) -> None:
    """`Inclusion_List.csv` rows (`m/z`, `Retention (min)`) in the Xcalibur mass-list layout, each
    with a retention window so the instrument targets the mass where the feature actually eluted."""
    sign = "Positive" if polarity in {"+", "Pos", "Positive", "positive"} else "Negative"
    lines = [THERMO_HEADER]
    for r in rows:
        rt = float(r["Retention (min)"])
        lines.append(f"{float(r['m/z']):.4f},,,{charge},{sign},{max(rt - window, 0.0):.2f},"
                     f"{rt + window:.2f},,,\"inclusion: MS1 feature never fragmented\"")
    Path(path).write_text("\n".join(lines) + "\n")


def _safety_report(path, polarity: str, checked: list[Checked], ppm: float, isolation: float,
                   n_library: int, n_identified: int, libraries: list) -> None:
    dropped = [c for c in checked if not c.safe]
    kept = [c for c in checked if c.safe]
    chimeric = [c for c in checked if c.id_iso and not c.id_tol]
    lines = [f"# Exclusion list safety, {polarity}", "",
             f"{len(checked)} candidate entries checked at ±{ppm:.0f} ppm (the exclusion mass "
             f"tolerance assumed for the method) against {n_identified} lipids identified in this "
             f"run and {n_library:,} library lipids ({len(libraries)} libraries).", "",
             f"- kept: **{len(kept)}** (written to the Thermo list)",
             f"- dropped, collide with a lipid identified here: "
             f"**{len([c for c in dropped if c.id_tol])}**",
             f"- dropped, collide with some library lipid: "
             f"**{len([c for c in dropped if c.lib_tol and not c.id_tol])}**",
             f"- kept but co-isolating (±{isolation} Da) with an identified lipid: "
             f"{len(chimeric)} (already chimeric today; excluding the contaminant helps, not hurts)",
             ""]
    if dropped:
        lines += ["## Dropped entries", "", "| m/z | MS2 wasted | series | collides with |", "|---|---|---|---|"]
        for c in sorted(dropped, key=lambda c: -int(c.entry["MS2 scans wasted"])):
            names = ", ".join(n for _, n in (c.id_tol or c.lib_tol)[:3])
            lines.append(f"| {float(c.entry['m/z']):.4f} | {c.entry['MS2 scans wasted']} | "
                         f"{c.entry.get('Series', '')} | {names} |")
        lines.append("")
    lines += ["The method's own exclusion tolerance decides how wide each entry really is; if the "
              f"method uses more than ±{ppm:.0f} ppm, re-run this check with that value.", ""]
    Path(path).write_text("\n".join(lines))


def write_method_lists(analysis: Path, ppm: float = DEFAULT_PPM,
                       isolation: float = ISOLATION_HALF_WIDTH,
                       inclusion_window: float = DEFAULT_INCLUSION_WINDOW,
                       say=print) -> dict:
    """Build `<analysis>/Method_Lists/` from whatever polarities `<analysis>/Results/` holds.

    Libraries are the ones the run itself searched (from `config_<pol>.json`), minus decoys, so
    the collision check sees exactly the lipids this configuration could name. Returns a summary
    dict (also written as README.md) so a caller can log it.
    """
    analysis = Path(analysis)
    out = analysis / "Method_Lists"
    out.mkdir(exist_ok=True)
    summary: dict = {"ppm": ppm, "isolation_half_width": isolation,
                     "inclusion_window_min": inclusion_window, "polarities": {}}
    for pol in ("Pos", "Neg"):
        res = analysis / "Results" / pol
        excl = res / "Exclusion_List.csv"
        incl = res / "Inclusion_List.csv"
        if not res.is_dir() or not (excl.exists() or incl.exists()):
            continue
        cfg = analysis / f"config_{pol}.json"
        libraries = json.loads(cfg.read_text()).get("libraries", []) if cfg.exists() else []
        info: dict = {}
        if excl.exists():
            entries = list(csv.DictReader(excl.open(newline="")))
            library = library_lipids(libraries, pol)
            identified = []
            for name in ("Final_Results.csv", "Final_Results_Filtered.csv"):
                if (res / name).exists():
                    identified = identified_lipids(res / name, pol)
                    break
            if pol == "Neg":
                identified = sorted(identified + protected_fatty_acids())
            checked = check_exclusion(entries, library, identified, ppm, isolation)
            kept = [c.entry for c in checked if c.safe]
            write_thermo_list([_wasted(e) for e in kept], out / f"Exclusion_List_{pol}_Thermo.csv",
                              "+" if pol == "Pos" else "-")
            _safety_report(out / f"Exclusion_Safety_{pol}.md", pol, checked, ppm, isolation,
                           len(library), len(identified), libraries)
            info["exclusion_candidates"] = len(entries)
            info["exclusion_kept"] = len(kept)
            info["exclusion_dropped"] = len(entries) - len(kept)
        if incl.exists():
            rows = list(csv.DictReader(incl.open(newline="")))
            write_inclusion_thermo(rows, out / f"Inclusion_List_{pol}_Thermo.csv", pol,
                                   inclusion_window)
            info["inclusion_entries"] = len(rows)
        summary["polarities"][pol] = info
        say(f"method lists {pol}: exclusion {info.get('exclusion_kept', 0)} kept of "
            f"{info.get('exclusion_candidates', 0)} ({info.get('exclusion_dropped', 0)} dropped "
            f"for colliding with a lipid), inclusion {info.get('inclusion_entries', 0)} entries")

    lines = ["# Method lists", "",
             "Mass lists for the next acquisition of this study, in the Xcalibur mass-list import "
             "layout, one file per polarity and purpose. Import the exclusion list into the "
             "method's exclusion (reject) mass list and the inclusion list into its inclusion "
             "list; nothing here has been applied to any instrument by the pipeline itself.", "",
             f"- exclusion entries were checked at ±{ppm:.0f} ppm against every lipid identified "
             f"in this run and every lipid in the searched libraries; colliding entries are "
             f"dropped and listed in `Exclusion_Safety_<pol>.md`. If the method's exclusion "
             f"tolerance is wider than ±{ppm:.0f} ppm, re-check with that value.",
             f"- inclusion entries carry a retention window of ±{inclusion_window} min around "
             f"the observed apex.",
             "- Start/End are empty on exclusion entries: contamination is present for the whole "
             "run and a window would leave a hole.", ""]
    for pol, info in summary["polarities"].items():
        lines.append(f"- {pol}: exclusion {info.get('exclusion_kept', 0)} of "
                     f"{info.get('exclusion_candidates', 0)} kept "
                     f"(`Exclusion_List_{pol}_Thermo.csv`), inclusion "
                     f"{info.get('inclusion_entries', 0)} (`Inclusion_List_{pol}_Thermo.csv`)")
    (out / "README.md").write_text("\n".join(lines) + "\n")
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary
