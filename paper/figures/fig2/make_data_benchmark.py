"""Figure 2C and 2D data: Set 2 (first arm, no lists) against the NIST SRM 1950 interlaboratory
consensus (Bowden 2017; reference table built by manuscript/benchmark/build_reference.py).

C: per class family, consensus entries (>= 5 laboratories; the paper's 339) recovered with MS2,
   recovered by the retention model only, or missed; plus the out-of-scope entries no spectral
   library can name (eicosanoids, bile acids, cholesterol, S1P, free fatty acids).
D: PE, the same consensus species measured in both polarities: log10 median area against log10
   consensus concentration, positive and negative mode side by side.
"""
from __future__ import annotations

import csv
import sys
from collections import defaultdict
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "manuscript" / "benchmark"))
import score_nist  # noqa: E402

RUN = Path("/mnt/work/lipidomics_work/NIST_v4_2026-09-11/Analysis_v8_2026-09-11/Results")
DATA = HERE / "data"


def main() -> None:
    ref = score_nist.load_reference()
    sp = score_nist.our_species(RUN)
    have_ms2 = {k for k, v in sp.items() if v["source"] == "MS2"}
    groups = defaultdict(list)
    for r in ref:
        if r["source"] == "Bowden2017" and r["tier"] != "Inf":
            groups[r["group"]].append(r)
    fam = defaultdict(lambda: {"ms2": 0, "rt_model": 0, "missed": 0})
    out_of_scope = 0
    pe = []
    for g, alts in groups.items():
        f = next((a["family"] for a in alts if a["kind"] != "other"), "other")
        if f == "other" or f not in score_nist.IN_SCOPE:
            out_of_scope += 1
            continue
        keys = [a["key"] for a in alts]
        if any(k in have_ms2 for k in keys):
            fam[f]["ms2"] += 1
        elif any(k in sp for k in keys):
            fam[f]["rt_model"] += 1
        else:
            fam[f]["missed"] += 1
        if f == "PE" and alts[0]["tier"] == "Ref":
            k = next((k for k in keys if k in sp), None)
            if k is not None:
                for pol, (area, cv, src, cv_norm) in sp[k]["by_pol"].items():
                    if src == "MS2" and area > 0:
                        pe.append({"consensus": alts[0]["name"], "polarity": pol, "value_nmol_mL": alts[0]["value_nmol_mL"],
                                   "area": area})
    rows = [{"family": f, **v, "n": v["ms2"] + v["rt_model"] + v["missed"]} for f, v in fam.items()]
    rows.sort(key=lambda r: -r["n"])
    rows.append({"family": "out of scope", "ms2": 0, "rt_model": 0, "missed": 0, "n": out_of_scope})
    with (DATA / "benchmark_recall.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["family", "n", "ms2", "rt_model", "missed"]); w.writeheader(); w.writerows(rows)
    with (DATA / "benchmark_pe.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["consensus", "polarity", "value_nmol_mL", "area"]); w.writeheader(); w.writerows(pe)
    print(rows); print(len(pe), "PE points")


if __name__ == "__main__":
    main()
