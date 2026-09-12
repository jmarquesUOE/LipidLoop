"""Figure 3 data: per-class retention models on Set 1 (A to D) and the target-versus-decoy
dot-product distributions on Set 1 and Set 2 (E).

Retention panels: every identified, non-decoy row of Set 1's unfiltered output, parsed to
(class, carbons, double bonds, retention) with the pipeline's own parser, and the model refitted
with the pipeline's own `fit_model` on the MS2-identified rows so the drawn line is the line the
pipeline used. Rows named by the retention model without MS2 are carried with their source so the
panel can show them in a different colour at the same size.

Decoy panel: rank-1 matches from each run's `Associated_Spectra.csv`, target versus `DECOY_`.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

RUNS = {
    "Set 1": Path("/mnt/work/lipidomics_work/PDAC_KPC_Set3/Analysis_v8_decoy_2026-09-10"),
    "Set 2": Path("/mnt/work/lipidomics_work/NIST_v4_2026-09-11/Analysis_v8_2026-09-11"),
}
# (set, class as `Class (canonical)`, polarity) for the four retention panels: the two most
# populated positive-mode classes on both in-house sets, so the same model is seen on two
# instruments and two gradients.
CLASSES = [("Set 1", "TG", "Pos"), ("Set 1", "PC", "Pos"), ("Set 2", "TG", "Pos"), ("Set 2", "PC", "Pos")]


def retention_rows() -> tuple[list[dict], list[dict]]:
    from lipidex_py.rtls import fit_model, parse_sum_name
    points, models = [], []
    for label, cls, pol in CLASSES:
        run = RUNS[label]
        rows = list(csv.DictReader((run / "Results" / pol / "Unfiltered_Results.csv").open()))
        fit_pts = []
        for r in rows:
            name = r.get("Identification", "").strip()
            if not name or name.upper().startswith("DECOY") or r.get("Class (canonical)") != cls:
                continue
            parsed = parse_sum_name(name)
            if parsed is None:
                continue
            _, prefix, c, db = parsed
            rt = float(r["Retention Time (min)"])
            src = r.get("Identification Source", "")
            points.append({"set": label, "class": cls, "polarity": pol, "prefix": prefix, "carbons": c,
                           "double_bonds": db, "retention": rt, "source": src,
                           "filter_status": r.get("Filter Status", "")})
            if src == "MS2" and prefix == "":
                fit_pts.append((c, db, rt))
        m = fit_model(cls, "", fit_pts)
        models.append({"set": label, "class": cls, "polarity": pol, "intercept": m.intercept, "per_carbon": m.per_carbon,
                       "per_double_bond": m.per_double_bond, "r2": m.r2, "residual_sd": m.residual_sd,
                       "n_used": m.n_used, "n_total": m.n_total, "usable": int(m.usable)})
        print(label, cls, pol, m)
    return points, models


def decoy_rows() -> list[dict]:
    out = []
    for label, run in RUNS.items():
        for pol in ("Pos", "Neg"):
            for r in csv.DictReader((run / "Results" / pol / "Associated_Spectra.csv").open()):
                if r.get("Rank", "1") not in ("1", ""):
                    continue
                out.append({"set": label, "polarity": pol,
                            "kind": "decoy" if r["Name"].upper().startswith("DECOY") else "target",
                            "dot": int(float(r["Dot Product"]))})
    return out


def write(name, rows):
    with (DATA / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    pts, models = retention_rows()
    write("retention_points.csv", pts)
    write("retention_models.csv", models)
    write("decoy_scores.csv", decoy_rows())
    print("written to", DATA)
