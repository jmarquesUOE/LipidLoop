"""Figure 5 data: wasted MS2 duty cycle per study and polarity (A), contaminants shared across
vendors (B), and Set 2 without versus with the pipeline's lists (C).

Writes `data/*.csv` beside this script. Sources: every run's own `Duty_Cycle.json` (the corpus in
its 2026-09-04 recalibrated pass, Set 1 in its decoy-matched run, Set 2 in both arms); the
cross-vendor contaminant table from Supporting Information Table S8 (derived from the same runs
by formula matching of recurring unassigned masses); and, for (C), the analysis directories plus
the mzML scan headers of every Set 2 replicate.
"""
from __future__ import annotations

import bisect
import csv
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)
WORK = Path("/mnt/work/lipidomics_work")

CORPUS = {  # dataset -> vendor
    "MTBKS222_Waters": "Waters", "MTBKS222_Thermo": "Thermo", "MTBKS222_Agilent": "Agilent",
    "MTBKS222_IMS_Normal": "Sciex", "MTBKS222_IMS_HighMass": "Sciex", "MTBKS222_BrukerBAF": "Bruker",
    "MTBLS5163": "Bruker", "ST003077": "Thermo", "ST003514": "Agilent", "ST004797": "Sciex",
    "ST004797_oxtg": "Sciex", "MSV000095868": "Agilent", "ST000991_DDA": "Sciex",
    "ST002705": "Bruker", "ST003052": "Thermo",
}
NIST = {"MTBKS222_Waters", "MTBKS222_Thermo", "MTBKS222_Agilent", "MTBKS222_IMS_Normal",
        "MTBKS222_IMS_HighMass", "MTBKS222_BrukerBAF", "ST003514"}
INHOUSE = {
    "Set 1 (Q Exactive)": (WORK / "PDAC_KPC_Set3" / "Analysis_v8_decoy_2026-09-10", "Thermo"),
    "Set 2 without lists (Fusion Lumos)": (WORK / "NIST_v4_2026-09-11" / "Analysis_v8_2026-09-11", "Thermo"),
    "Set 2 with lists (Fusion Lumos)": (WORK / "NIST_v5_2026-09-11" / "Analysis_v8_2026-09-11", "Thermo"),
}
NAMED = ["polysiloxane", "sodium formate", "sodium acetate", "PEG", "PPG", "PTFE/PFPE"]


def wasted_rows() -> list[dict]:
    out = []
    def add(label, analysis, vendor, nist, group):
        for pol in ("Pos", "Neg"):
            f = analysis / "Results" / pol / "Duty_Cycle.json"
            if not f.exists():
                continue
            d = json.loads(f.read_text())
            row = {"dataset": label, "vendor": vendor, "nist": int(nist), "group": group,
                   "polarity": pol, "total_ms2": d["total_ms2"],
                   "wasted": 100 * d["duty_cycle_wasted_fraction"]}
            other = 0.0
            for s in d["by_series"]:
                name = s["series"]
                if name in NAMED or name == "unassigned":
                    row[name] = row.get(name, 0.0) + 100 * s["fraction"]
                else:
                    other += 100 * s["fraction"]
            for name in NAMED + ["unassigned"]:
                row.setdefault(name, 0.0)
            row["other named series"] = other
            out.append(row)
    for ds, vendor in CORPUS.items():
        add(ds, WORK / ds / "Analysis_recal_2026-09-04", vendor, ds in NIST, "corpus")
    for label, (analysis, vendor) in INHOUSE.items():
        add(label, analysis, vendor, label.startswith("Set 2"), "in-house")
    return out


def contaminants() -> list[dict]:
    # Supporting Information Table S8, cross-vendor table (derived from the same runs).
    rows = [
        ("PEG", "series matcher", "Agilent;Bruker;Sciex;Thermo", 95428),
        ("sodium acetate", "series matcher", "Agilent;Bruker;Sciex;Thermo", 49872),
        ("DEHP (phthalate)", "formula match", "Agilent;Bruker;Sciex;Thermo", 16342),
        ("cholesterol in-source fragment", "formula match", "Agilent;Bruker;Sciex;Thermo", 6711),
        ("polysiloxane", "series matcher", "Bruker;Sciex;Thermo", 191958),
        ("sodium formate", "series matcher", "Bruker;Sciex;Thermo", 235807),
        ("free fatty acids", "formula match", "Agilent;Sciex;Thermo", 12244),
        ("DBP, erucamide", "formula match", "Bruker;Sciex;Thermo", 13169),
    ]
    return [{"contaminant": c, "found_by": f, "vendors": v, "wasted_scans": n} for c, f, v, n in rows]


def _ms2_precursors(path: Path) -> list[float]:
    txt = path.read_text(errors="ignore")
    out = []
    for s in re.findall(r"<spectrum [^>]*>(.*?)</spectrum>", txt, flags=re.S):
        if re.search(r'MS:1000511[^>]*value="2"', s):
            m = re.search(r'MS:1000744[^>]*value="([\d.]+)"', s)
            if m:
                out.append(float(m.group(1)))
    return out


def paired() -> list[dict]:
    from lipidloop.peaks import sum_composition
    arms = {"without lists": WORK / "NIST_v4_2026-09-11", "with lists": WORK / "NIST_v5_2026-09-11"}
    lists = WORK / "NIST_v4_2026-09-11" / "Analysis_v8_2026-09-11" / "Method_Lists"
    out = []
    for pol in ("Pos", "Neg"):
        excl = sorted(float(r["Mass [m/z]"]) for r in csv.DictReader((lists / f"Exclusion_List_{pol}_Thermo.csv").open()))
        for arm, w in arms.items():
            a = w / "Analysis_v8_2026-09-11"
            d = json.loads((a / "Results" / pol / "Duty_Cycle.json").read_text())
            rows = [r for r in csv.DictReader((a / "Results" / pol / "Final_Results_Filtered.csv").open())
                    if r.get("Identification", "").strip() and not r["Identification"].upper().startswith("DECOY")]
            ms2_conf = [r for r in rows if r.get("Identification Source", "") == "MS2"]
            # MS2 scans on an excluded mass, all three replicates pooled
            n = hit = 0
            for mz_file in sorted((w / "mzml" / pol).glob(f"NIST_{pol}_0*.mzML")):
                for mz in _ms2_precursors(mz_file):
                    n += 1
                    i = bisect.bisect_left(excl, mz - mz * 10e-6)
                    if i < len(excl) and excl[i] <= mz + mz * 10e-6:
                        hit += 1
            out.append({"polarity": pol, "arm": arm, "wasted": 100 * d["duty_cycle_wasted_fraction"],
                        "on_excluded_mass": 100 * hit / n, "ms2_scans_replicates": n,
                        "ms2_confirmed": len(ms2_conf), "delivered": len(rows),
                        "molecules": len({sum_composition(r["Identification"]) for r in ms2_conf})})
            print(out[-1])
    return out


def write(name: str, rows: list[dict]) -> None:
    with (DATA / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    write("wasted.csv", wasted_rows())
    write("contaminants_by_vendor.csv", contaminants())
    write("paired.csv", paired())
    print("written to", DATA)
