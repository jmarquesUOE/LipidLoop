"""Compare the 25-min (2026) and 30-min (2024) gradients on separation quality.

The confound is that these methods were run on different samples: skin organoid against rat
heart. Anything computed on all identifications therefore measures the lipidome as much as the
chromatography. So the headline comparison here is restricted to the identifications the two
datasets SHARE - same molecules, two methods - and the all-identification numbers are reported
only as context.

Neither the LipiDex export nor the Compound Discoverer export carries peak width, so classical
peak capacity is not recoverable without re-reading the .raw files. What is recoverable from
RT + quant ion:

  A  window usage      how much of the run the lipids actually occupy
  B  co-isolation      how many other identified lipids sit inside the MS2 isolation window at
                       the same retention time - this is what produces chimeric MS2, and
                       chimeric MS2 is what fails the LipiDex 75% purity filter
  C  homologue         seconds of retention gained per CH2 and per double bond within a class;
     selectivity       the direct measure of resolving power for the isomer-adjacent species
                       that lipidomics actually has to separate
  D  class order       whether classes come off as distinct blocks or overlap

Caveat carried into the output: these are two runs two years apart, so column age and
instrument state are inside the comparison as well as the gradient.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

S = Path("/tmp/claude-1000/-home-jair/7a1e4579-31fb-421c-93fb-62aa8acacc72/scratchpad")
OUT = Path(__file__).resolve().parents[1] / "outputs"

# gradient window = first to last time the mobile phase is still changing
# flow rate is read from the .raw header (PumpModule.LoadingPump.Flow.Nominal) and it is
# NOT the same between these methods: 150 uL/min in 2024 against 300 in 2026, exactly 2x.
# Retention time scales as 1/flow, so every time-domain measure - seconds per CH2, peak
# spacing, co-elution counts - is inflated 2x for the 2024 run by flow alone. Everything
# below is therefore computed in ELUTION VOLUME (uL = RT x flow), in which flow cancels.
# Resolution is flow-invariant to first order because peak width scales as 1/flow too.
SETS = {
    "skin_pos": (S / "skin2024/Pos_Final_Results.csv", "2024 / 30 min / 150 uL", 30.0, 150.0),
    "skin_neg": (S / "skin2024/Neg_Final_Results.csv", "2024 / 30 min / 150 uL", 30.0, 150.0),
    "ckd_pos":  (S / "dunja2/Pos_Final_Results.csv",   "2026 / 25 min / 300 uL", 25.0, 300.0),
    "ckd_neg":  (S / "dunja2/Neg_Final_Results.csv",   "2026 / 25 min / 300 uL", 25.0, 300.0),
}
ISOLATION_MZ = 1.0     # +/- Da, typical quadrupole isolation for DDA
COELUTE_UL = 15.0      # +/- uL of eluent, a flow-independent stand-in for a peak width


def parse_chain(ident: str):
    """Total carbons and double bonds summed over all acyl chains in an identification."""
    chains = re.findall(r"(\d+):(\d+)", str(ident))
    if not chains:
        return np.nan, np.nan
    c = sum(int(a) for a, _ in chains)
    d = sum(int(b) for _, b in chains)
    return c, d


def load(path: Path, flow: float) -> pd.DataFrame:
    d = pd.read_csv(path)
    d = d.loc[:, ~d.columns.str.startswith("Unnamed")]
    d = d[d["Identification"].notna()].copy()
    d["rt"] = pd.to_numeric(d["Retention Time (min)"], errors="coerce")
    d["mz"] = pd.to_numeric(d["Quant Ion"], errors="coerce")
    d["cls"] = d["Lipid Class"].astype(str)
    d["ident"] = d["Identification"].astype(str).str.strip()
    d[["nC", "nDB"]] = d["ident"].apply(lambda s: pd.Series(parse_chain(s)))
    d = d.dropna(subset=["rt", "mz"]).reset_index(drop=True)
    d["flow"] = flow
    d["ev"] = d["rt"] * flow          # elution volume, uL - the flow-independent axis
    return d


def coisolation(d: pd.DataFrame) -> pd.Series:
    """Per lipid, how many OTHER identified lipids fall inside the isolation window at the
    same RT. A precursor with neighbours here yields chimeric MS2."""
    ev, mz = d["ev"].to_numpy(), d["mz"].to_numpy()
    near = (np.abs(ev[:, None] - ev[None, :]) <= COELUTE_UL) & \
           (np.abs(mz[:, None] - mz[None, :]) <= ISOLATION_MZ)
    np.fill_diagonal(near, False)
    return pd.Series(near.sum(1), index=d.index)


def selectivity(d: pd.DataFrame, min_n: int = 8) -> pd.DataFrame:
    """Within each class, regress RT on carbon number and double-bond count.
    Slopes are seconds of retention per CH2 and per double bond."""
    rows = []
    for cls, g in d.groupby("cls"):
        g = g.dropna(subset=["nC", "nDB"])
        if len(g) < min_n or g["nC"].nunique() < 3:
            continue
        X = np.column_stack([np.ones(len(g)), g["nC"], g["nDB"]])
        try:
            beta, *_ = np.linalg.lstsq(X, g["ev"].to_numpy(), rcond=None)
        except np.linalg.LinAlgError:
            continue
        pred = X @ beta
        ss_res = float(((g["ev"] - pred) ** 2).sum())
        ss_tot = float(((g["ev"] - g["ev"].mean()) ** 2).sum())
        rows.append({"class": cls, "n": len(g),
                     "uL_per_CH2": beta[1], "uL_per_DB": beta[2],
                     "r2": 1 - ss_res / ss_tot if ss_tot > 0 else np.nan,
                     "ev_span_uL": g["ev"].max() - g["ev"].min()})
    return pd.DataFrame(rows).sort_values("n", ascending=False)


def summarise(d: pd.DataFrame, run_min: float, label: str) -> dict:
    co = coisolation(d)
    lo, hi = d["rt"].quantile([0.01, 0.99])
    elo, ehi = d["ev"].quantile([0.01, 0.99])
    gaps = np.diff(np.sort(d["ev"].to_numpy()))
    return {"set": label, "n": len(d),
            "window_min": hi - lo, "pct_of_run": (hi - lo) / run_min * 100,
            "window_uL": ehi - elo,
            "ids_per_mL": len(d) / (ehi - elo) * 1000,
            "median_gap_uL": float(np.median(gaps)),
            "pct_coisolated": float((co > 0).mean() * 100),
            "median_neighbours": float(co.median()),
            "mean_neighbours": float(co.mean())}


def main() -> None:
    data = {k: load(p, f) for k, (p, _, _, f) in SETS.items()}
    lines = []

    def w(s=""):
        print(s)
        lines.append(s)

    w("=" * 78)
    w("METHOD COMPARISON - 2024 (30 min) against 2026 (25 min)")
    w("=" * 78)
    w()
    w("FLOW RATE DIFFERS: 150 uL/min in 2024, 300 in 2026 (from the .raw headers), so all")
    w("measures below are in elution VOLUME, not time. Column oven also differs, 40 vs 50 C.")
    w("IPA ramp: 2024 10-22 min = 1800 uL;  2026 12-20 min = 2400 uL - so the 2026 ramp is")
    w("volumetrically SHALLOWER despite being shorter in time.")
    w()

    w("-" * 78)
    w("1. ALL IDENTIFICATIONS  (confounded by sample - context only)")
    w("-" * 78)
    rows = [summarise(data[k], SETS[k][2], f"{k}  [{SETS[k][1]}]") for k in SETS]
    allsum = pd.DataFrame(rows)
    w(allsum.to_string(index=False, float_format=lambda x: f"{x:.2f}"))
    w()
    w("  ids_per_mL        identifications per mL of eluent (crowding, flow-independent)")
    w("  median_gap_uL     microlitres of eluent between adjacent identifications")
    w("  pct_coisolated    % of lipids with another identified lipid inside +/-1 Da, +/-15 uL")
    w()

    # ---- 2. paired on shared identifications ------------------------------------
    w("-" * 78)
    w("2. PAIRED ON SHARED IDENTIFICATIONS  (same molecules, two methods)")
    w("-" * 78)
    for pol in ("pos", "neg"):
        a, b = data[f"skin_{pol}"], data[f"ckd_{pol}"]
        # collapse duplicate identifications to the most intense
        ac = a.sort_values("Area (max)", ascending=False).drop_duplicates("ident")
        bc = b.sort_values("Area (max)", ascending=False).drop_duplicates("ident")
        shared = sorted(set(ac["ident"]) & set(bc["ident"]))
        if len(shared) < 20:
            w(f"  {pol}: only {len(shared)} shared identifications, skipped")
            continue
        A = ac.set_index("ident").loc[shared]
        B = bc.set_index("ident").loc[shared]
        sa = summarise(A.reset_index(), 30.0, f"skin {pol} [2024/30min]")
        sb = summarise(B.reset_index(), 25.0, f"ckd  {pol} [2026/25min]")
        w(f"  {pol.upper()}  -  {len(shared)} identifications present in both runs")
        w(pd.DataFrame([sa, sb]).to_string(index=False,
                                           float_format=lambda x: f"{x:.2f}"))
        r = np.corrcoef(A["rt"], B["rt"])[0, 1]
        w(f"    RT correlation between methods r = {r:.3f}  "
          f"(elution ORDER is {'preserved' if r > 0.9 else 'NOT preserved'})")
        w(f"    span: 2024 {A['rt'].max()-A['rt'].min():.1f} min, "
          f"2026 {B['rt'].max()-B['rt'].min():.1f} min")
        w()

    # ---- 3. homologue selectivity ------------------------------------------------
    w("-" * 78)
    w("3. HOMOLOGUE SELECTIVITY  (seconds of RT per CH2 / per double bond, within class)")
    w("-" * 78)
    sel = {}
    for k in SETS:
        s = selectivity(data[k])
        s["set"] = k
        sel[k] = s
    for pol in ("pos", "neg"):
        sk, ck = sel[f"skin_{pol}"], sel[f"ckd_{pol}"]
        m = sk.merge(ck, on="class", suffixes=("_2024", "_2026"))
        if not len(m):
            continue
        w(f"  {pol.upper()}  - classes resolvable in both")
        cols = ["class", "n_2024", "uL_per_CH2_2024", "uL_per_DB_2024", "r2_2024",
                "n_2026", "uL_per_CH2_2026", "uL_per_DB_2026", "r2_2026"]
        w(m[cols].to_string(index=False, float_format=lambda x: f"{x:.2f}"))
        w(f"    median uL/CH2:  2024 {m['uL_per_CH2_2024'].median():.1f}   "
          f"2026 {m['uL_per_CH2_2026'].median():.1f}")
        w(f"    median uL/DB :  2024 {m['uL_per_DB_2024'].median():.1f}   "
          f"2026 {m['uL_per_DB_2026'].median():.1f}")
        w(f"    median r2     :  2024 {m['r2_2024'].median():.3f}   "
          f"2026 {m['r2_2026'].median():.3f}   (fit quality = orderly homologue elution)")
        w()
        m.to_csv(OUT / f"method_selectivity_{pol}.csv", index=False)

    # ---- 4. ceramides specifically ----------------------------------------------
    w("-" * 78)
    w("4. CERAMIDES IN NEGATIVE MODE  (the Cer[ADS] screening question)")
    w("-" * 78)
    for k in ("skin_neg", "ckd_neg"):
        d = data[k]
        cer = d[d["cls"].str.contains("Cer", case=False, na=False)]
        w(f"  {k} [{SETS[k][1]}]  {len(cer)} ceramide identifications")
        for cls, g in cer.groupby("cls"):
            co = coisolation(d).loc[g.index]
            w(f"    {cls:<16} n={len(g):<4} RT {g['rt'].min():5.2f}-{g['rt'].max():5.2f} "
              f"span {g['ev'].max()-g['ev'].min():6.0f} uL   "
              f"co-isolated {100*(co>0).mean():5.1f}%")
        w()

    allsum.to_csv(OUT / "method_comparison_summary.csv", index=False)
    (OUT / "method_comparison.txt").write_text("\n".join(lines))
    print(f"\nwritten: {OUT}/method_comparison.txt (+ csv)")


if __name__ == "__main__":
    main()
