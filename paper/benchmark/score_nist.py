"""Score a pipeline run on NIST SRM 1950 against the published lists (see build_reference.py).

For each run (a Results/ folder with Pos/ and Neg/ Final_Results.csv):

  * every identified row is reduced to the same key as the reference (family, kind, C, DB) at
    sum composition; a species is `MS2` if any row carrying it was identified from a spectrum,
    else `RT model`; its area is the median over the NIST replicate columns of the polarity with
    the larger median, its CV the replicate CV in that polarity.
  * consensus recall: the fraction of Bowden 2017 consensus entries (per tier, per family) whose
    key, or any key of the isobaric group, is among our species. Reported for MS2-only and for
    MS2 + retention-model species. Entries outside what the libraries can name (eicosanoids, bile
    acids, cholesterol, sphingoid-base phosphates) are listed separately and count against nothing.
  * precision proxy: the fraction of our species, in families any list covers, that appear in at
    least one published list (Bowden all tiers + Quehenberger 2010 + Godzien 2024). A species
    absent from every list is not shown wrong, only unconfirmed; the number is a floor.
  * quantitative agreement, within each family with >= 5 matched Ref-tier entries: Spearman rho
    of log10(median area) against log10(consensus nmol/mL). No response factors are applied, so
    only the within-family rank order is a fair test.
  * dispersion: the median replicate CV of matched species next to the median inter-laboratory
    COD of the same consensus entries.

    .venv/bin/python manuscript/benchmark/score_nist.py v3=/path/Results v4=/path/Results ...
"""
from __future__ import annotations

import csv
import math
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import median

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parents[1] / "src"))
from build_reference import key_of  # noqa: E402
from lipidloop.peaks import sum_composition  # noqa: E402

REF = HERE / "reference" / "srm1950_reference_lists.csv"
OUT = HERE / "results"
META_PREFIX = ("Blank", "QC", "blank", "qc")
# Families a spectral identification can reach. Free fatty acids are out: the carboxylate anion
# gives no diagnostic HCD fragment, the pipeline names them only as a homologous series, and on
# a DDA lipid method that series is not there (Jair, 2026-09-12: "the method is mainly for lipids
# with MS2 anyway"). They are reported with the other out-of-scope entries, never as recall.
IN_SCOPE = {"PC", "PE", "LPC", "LPE", "PI", "PS", "PG", "PA", "LPI", "SM", "Cer", "HexCer", "TG", "DG", "CE", "AC", "CerOH"}


def load_reference():
    rows = list(csv.DictReader(REF.open()))
    for r in rows:
        r["key"] = (r["family"], r["kind"], int(r["carbons"]) if r["carbons"] else None,
                    int(r["double_bonds"]) if r["double_bonds"] else None)
    return rows


def our_species(results: Path) -> dict:
    """key -> {name, source, area, cv, polarity}"""
    species: dict = {}
    for pol in ("Pos", "Neg"):
        path = results / pol / "Final_Results.csv"
        if not path.exists():
            continue
        rows = list(csv.DictReader(path.open()))
        cols = [c for c in rows[0].keys() if c and c not in _META and not c.startswith(META_PREFIX)]
        # per-column median over identified rows: the replicate CV is also computed after dividing
        # each column by it, which removes an injection-volume difference between replicates
        ident = [r for r in rows if r.get("Identification", "").strip()]
        colmed = {c: float(np.median([float(r[c]) for r in ident if r[c] and float(r[c]) > 0] or [1.0])) for c in cols}
        for r in rows:
            name = r.get("Identification", "").strip()
            if not name or name.upper().startswith("DECOY"):
                continue
            key = key_of(sum_composition(name))
            if key[1] == "other":
                key = (key[0], "other", None, None)
            areas = np.array([float(r[c] or 0) for c in cols])
            areas = areas[areas > 0]
            area = float(np.median(areas)) if len(areas) else 0.0
            cv = float(np.std(areas, ddof=1) / np.mean(areas) * 100) if len(areas) > 1 else float("nan")
            normed = np.array([float(r[c] or 0) / colmed[c] for c in cols]); normed = normed[normed > 0]
            cv_norm = float(np.std(normed, ddof=1) / np.mean(normed) * 100) if len(normed) > 1 else float("nan")
            src = "MS2" if r.get("Identification Source", "") == "MS2" else "RT model"
            cur = species.get(key)
            if cur is None or (src == "MS2" and cur["source"] != "MS2") or (src == cur["source"] and area > cur["area"]):
                species[key] = {"name": name, "source": src, "area": area, "cv": cv, "cv_norm": cv_norm, "polarity": pol,
                                "family": key[0] if key[1] != "other" else "other",
                                "by_pol": cur["by_pol"] if cur else {}}
            bp = species[key]["by_pol"]
            if pol not in bp or area > bp[pol][0]:
                bp[pol] = (area, cv, src, cv_norm)
    return species


_META = {"Compound Group", "Retention Time (min)", "Quant Ion", "Polarity", "Area (max)", "Identification", "Lipid Key",
         "Shorthand (LSI)", "Lipid Class", "Features Found", "Class (canonical)", "Adduct", "Dot Product",
         "Reverse Dot Product", "Purity", "Chain Evidence", "Purity Source", "MS2 Spectra", "MS2 Files",
         "Identification Source", "Duplicate Name", "RT Model Error", "RT Model Z", "RT Model R2", "Duplicate Verdict",
         "S/N", "Theoretical m/z", "Mass Error (ppm)", "Filter Status", "Quant Ion (measured)", "Isotope Pattern",
         "Retention Model Error (min)"}


def spearman(x, y):
    from scipy.stats import spearmanr
    r = spearmanr(x, y)
    return float(r.correlation), float(r.pvalue)


def score(label: str, results: Path, ref: list[dict]) -> dict:
    sp = our_species(results)
    have_ms2 = {k for k, v in sp.items() if v["source"] == "MS2"}
    have_all = set(sp)
    out = {"run": label, "species_total": len(have_all), "species_ms2": len(have_ms2)}

    # --- consensus recall (Bowden), by tier and by family, groups recovered by any alternative
    groups: dict = defaultdict(list)
    for r in ref:
        if r["source"] == "Bowden2017":
            groups[r["group"]].append(r)
    recall_rows = []
    per_tier = defaultdict(lambda: {"n": 0, "ms2": 0, "all": 0, "out_of_scope": 0})
    per_family = defaultdict(lambda: {"n": 0, "ms2": 0, "all": 0})
    for g, alts in groups.items():
        tier = alts[0]["tier"]
        keys = [a["key"] for a in alts]
        fam = next((a["family"] for a in alts if a["kind"] != "other"), "other")
        if fam == "other" or fam not in IN_SCOPE:
            per_tier[tier]["out_of_scope"] += 1
            continue
        hit_ms2 = any(k in have_ms2 for k in keys)
        hit_all = any(k in have_all for k in keys)
        per_tier[tier]["n"] += 1; per_tier[tier]["ms2"] += hit_ms2; per_tier[tier]["all"] += hit_all
        if tier != "Inf":
            per_family[fam]["n"] += 1; per_family[fam]["ms2"] += hit_ms2; per_family[fam]["all"] += hit_all
        matched = next((k for k in keys if k in have_all), None)
        recall_rows.append({"run": label, "tier": tier, "consensus": alts[0]["name"], "family": fam,
                            "value_nmol_mL": alts[0]["value_nmol_mL"], "n_labs": alts[0]["n_labs"], "cod_pct": alts[0]["cod_pct"],
                            "found": "MS2" if hit_ms2 else ("RT model" if hit_all else ""),
                            "our_name": sp[matched]["name"] if matched else "", "our_area": sp[matched]["area"] if matched else "",
                            "our_cv_pct": sp[matched]["cv"] if matched else "", "our_key": matched})
    out["recall_by_tier"] = {t: dict(v) for t, v in per_tier.items()}
    out["recall_by_family"] = {f: dict(v) for f, v in per_family.items()}
    both = {"n": per_tier["Ref"]["n"] + per_tier["COD>40"]["n"], "ms2": per_tier["Ref"]["ms2"] + per_tier["COD>40"]["ms2"],
            "all": per_tier["Ref"]["all"] + per_tier["COD>40"]["all"]}
    out["recall_339"] = both

    # --- precision proxy: our species in families any list covers, found in any list
    listed = {r["key"] for r in ref if r["kind"] != "other"}
    fam_listed = {r["family"] for r in ref if r["kind"] != "other"}
    prec = defaultdict(lambda: {"n": 0, "listed": 0})
    for k, v in sp.items():
        if v["family"] in fam_listed:
            bucket = prec[(v["family"], v["source"])]
            bucket["n"] += 1; bucket["listed"] += k in listed
    out["precision"] = {f"{f}|{s}": dict(v) for (f, s), v in prec.items()}
    tot = defaultdict(lambda: {"n": 0, "listed": 0})
    for (f, s), v in prec.items():
        tot[s]["n"] += v["n"]; tot[s]["listed"] += v["listed"]
    out["precision_total"] = {s: dict(v) for s, v in tot.items()}
    out["unlisted"] = sorted(v["name"] for k, v in sp.items() if v["family"] in fam_listed and k not in listed)

    # --- quantitative agreement, Ref tier, within family
    # one polarity per family (the one carrying more matched MS2 entries), so areas share a scale
    quant_pol = defaultdict(lambda: defaultdict(list))
    for r in recall_rows:
        if r["tier"] == "Ref" and r["found"] == "MS2" and r["our_key"] is not None:
            for pol, (area, cv, src, cv_norm) in sp[r["our_key"]]["by_pol"].items():
                if src == "MS2" and area > 0:
                    quant_pol[r["family"]][pol].append((math.log10(float(r["value_nmol_mL"])), math.log10(area),
                                                       cv, float(r["cod_pct"]) if r["cod_pct"] else float("nan"), cv_norm))
    quant, quant_polarity = {}, {}
    for fam, d in quant_pol.items():
        pol = max(d, key=lambda p: len(d[p]))
        quant[fam] = d[pol]; quant_polarity[fam] = pol
    rho = {}
    allx, ally = [], []
    for fam, pts in quant.items():
        if len(pts) >= 5:
            x, y = zip(*[(p[0], p[1]) for p in pts])
            r_, p_ = spearman(x, y)
            rho[fam] = {"n": len(pts), "rho": round(r_, 2), "p": p_, "polarity": quant_polarity[fam]}
        for p in pts:
            allx.append(p[0]); ally.append(p[1])
    out["spearman"] = rho
    out["spearman_pooled"] = {"n": len(allx), "rho": round(spearman(allx, ally)[0], 2)} if len(allx) > 5 else {}
    cvs = [p[2] for pts in quant.values() for p in pts if not math.isnan(p[2])]
    cods = [p[3] for pts in quant.values() for p in pts if not math.isnan(p[3])]
    cvn = [p[4] for pts in quant.values() for p in pts if not math.isnan(p[4])]
    out["dispersion"] = {"median_replicate_cv_pct": round(median(cvs), 1) if cvs else None,
                         "median_replicate_cv_normalised_pct": round(median(cvn), 1) if cvn else None,
                         "median_interlab_cod_pct": round(median(cods), 1) if cods else None, "n": len(cvs)}
    return out, recall_rows, sp


def main(argv):
    ref = load_reference()
    OUT.mkdir(exist_ok=True)
    summary, all_recall = [], []
    for arg in argv:
        label, path = arg.split("=", 1)
        res, recall_rows, sp = score(label, Path(path), ref)
        summary.append(res); all_recall += recall_rows
        with (OUT / f"species_{label}.csv").open("w", newline="") as fh:
            w = csv.writer(fh); w.writerow(["family", "kind", "carbons", "double_bonds", "name", "source", "polarity", "area", "cv_pct", "cv_normalised_pct", "in_any_list"])
            listed = {r["key"] for r in ref if r["kind"] != "other"}
            for k, v in sorted(sp.items(), key=lambda kv: (str(kv[0][0]), str(kv[0][1]), kv[0][2] or 0, kv[0][3] or 0)):
                w.writerow([k[0], k[1], k[2], k[3], v["name"], v["source"], v["polarity"], f"{v['area']:.0f}", f"{v['cv']:.1f}", f"{v['cv_norm']:.1f}", int(k in listed)])
    with (OUT / "consensus_recall_rows.csv").open("w", newline="") as fh:
        for r in all_recall:
            r.pop("our_key", None)
        w = csv.DictWriter(fh, fieldnames=list(all_recall[0].keys())); w.writeheader(); w.writerows(all_recall)
    # ---- report
    lines = ["# NIST SRM 1950 benchmark", "",
             "Consensus = Bowden 2017 (NIST interlaboratory exercise); Ref = >= 5 labs and COD <= 40 % (the LipidQC benchmark set), "
             "COD>40 = >= 5 labs with COD > 40 %, Ref + COD>40 = the paper's 339; Inf = 3-4 labs. Out-of-scope entries "
             "(eicosanoids, bile acids, cholesterol, S1P, free fatty acids) are excluded from the denominators.", ""]
    hdr = "| run | species (MS2 / all) | Ref recall MS2 | Ref recall MS2+RT | 339 recall MS2 | 339 recall MS2+RT | Inf recall MS2+RT | listed (MS2) | listed (RT model) | pooled rho (n) | median rep CV raw / column-normalised / interlab COD |"
    lines += [hdr, "|" + "---|" * 11]
    for s in summary:
        t = s["recall_by_tier"]; r339 = s["recall_339"]; pt = s["precision_total"]
        def pct(a, n): return f"{a}/{n} ({100*a/n:.0f}%)" if n else "-"
        lines.append(f"| {s['run']} | {s['species_ms2']} / {s['species_total']} | {pct(t['Ref']['ms2'], t['Ref']['n'])} | {pct(t['Ref']['all'], t['Ref']['n'])} | "
                     f"{pct(r339['ms2'], r339['n'])} | {pct(r339['all'], r339['n'])} | {pct(t['Inf']['all'], t['Inf']['n'])} | "
                     f"{pct(pt.get('MS2', {}).get('listed', 0), pt.get('MS2', {}).get('n', 0))} | {pct(pt.get('RT model', {}).get('listed', 0), pt.get('RT model', {}).get('n', 0))} | "
                     f"{s['spearman_pooled'].get('rho', '-')} ({s['spearman_pooled'].get('n', 0)}) | {s['dispersion']['median_replicate_cv_pct']} % / {s['dispersion']['median_replicate_cv_normalised_pct']} % / {s['dispersion']['median_interlab_cod_pct']} % (n={s['dispersion']['n']}) |")
    lines += ["", "## Consensus recall by family (Ref + COD>40, MS2 / MS2+RT of n)", ""]
    fams = sorted({f for s in summary for f in s["recall_by_family"]}, key=lambda f: -max(s["recall_by_family"].get(f, {"n": 0})["n"] for s in summary))
    lines.append("| family | n | " + " | ".join(s["run"] for s in summary) + " |"); lines.append("|---|---|" + "---|" * len(summary))
    for f in fams:
        n = max(s["recall_by_family"].get(f, {"n": 0})["n"] for s in summary)
        cells = [f"{s['recall_by_family'].get(f, {'ms2': 0})['ms2']} / {s['recall_by_family'].get(f, {'all': 0})['all']}" for s in summary]
        lines.append(f"| {f} | {n} | " + " | ".join(cells) + " |")
    lines += ["", "## Precision proxy by family (our species found in any published list / our species), MS2 only", ""]
    pf = sorted({k.split("|")[0] for s in summary for k in s["precision"]})
    lines.append("| family | " + " | ".join(s["run"] for s in summary) + " |"); lines.append("|---|" + "---|" * len(summary))
    for f in pf:
        cells = []
        for s in summary:
            v = s["precision"].get(f"{f}|MS2", {"n": 0, "listed": 0})
            cells.append(f"{v['listed']}/{v['n']}" if v["n"] else "-")
        lines.append(f"| {f} | " + " | ".join(cells) + " |")
    lines += ["", "## Within-family Spearman rho, log area vs log consensus (Ref tier, MS2, one polarity per family, families with >= 5 matches)", ""]
    sf = sorted({f for s in summary for f in s["spearman"]})
    lines.append("| family | " + " | ".join(s["run"] for s in summary) + " |"); lines.append("|---|" + "---|" * len(summary))
    for f in sf:
        lines.append(f"| {f} | " + " | ".join(f"{s['spearman'][f]['rho']} (n={s['spearman'][f]['n']}, {s['spearman'][f]['polarity']})" if f in s["spearman"] else "-" for s in summary) + " |")
    lines += ["", "## Out of scope (Bowden entries no spectral library can name)", ""]
    oos = sorted({g["name"] for g in ref if g["source"] == "Bowden2017" and (g["kind"] == "other" or g["family"] not in IN_SCOPE)})
    lines.append(", ".join(oos))
    for s in summary:
        lines += ["", f"## {s['run']}: MS2 species in covered families that no list carries ({len(s['unlisted'])})", "", ", ".join(s["unlisted"])]
    (OUT / "REPORT.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines[:12]))
    print("written", OUT / "REPORT.md")


if __name__ == "__main__":
    main(sys.argv[1:])
