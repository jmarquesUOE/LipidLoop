"""Compare our quantitation against the facility reference, on the same raw files.

Identification has been validated exhaustively — 14,959 matches reproducing LipiDex to the
rounding. **Quantitation never has.** In the reference workflow Compound Discoverer produces the
areas and LipiDex only identifies; this pipeline replaces BOTH halves, and the half that assigns
numbers has never been compared to anything.

Nothing here needs instrument time: both sides already processed the same injections, so for every
molecule they agree on there are two independent area series to compare.

Two questions, and the second is the one that matters:

  * **Do the two agree on how a molecule behaves across samples?** Correlation of the log areas,
    per molecule. This is what a differential analysis actually depends on — an absolute scale
    factor is harmless, a per-molecule disagreement is not.
  * **Is there a systematic scale difference?** The median ratio. A constant factor points at a
    structural difference in what is being integrated rather than at noise.

⚠ **Parsing these two tables is where the errors live, not the statistics.** Four attempts were
needed: `str.isdigit()` silently rejects scientific notation, so half the areas read as missing;
sample columns were detected by file extension, which ours do not carry; and the substitutions that
strip `.raw` and `(F1)` were applied in the wrong order, so nothing matched at all — which
presented as "0 samples in common" rather than as an error. Each failure looked like a result.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
from lipidloop.peaks import sum_composition  # noqa: E402

#: Columns that are not per-sample areas. Anything else is treated as a sample, because our tables
#: carry no file extension and the reference's do — no pattern matches both.
META = {
    "Retention Time (min)", "Quant Ion", "Polarity", "Area (max)", "Area (Max.)", "Identification",
    "Lipid Class", "Features Found", "Dot Product", "Reverse Dot Product", "Purity", "MS2 Spectra",
    "MS2 Files", "Identification Source", "Filter Status", "Adduct", "Checked", "Name",
    "Molecular Weight", "RT [min]", "MS2", "Shorthand (LSI)", "Lipid Key", "Class (canonical)",
    "Compound Group", "Chain Evidence", "Purity Source", "RT Model Error", "RT Model Z",
    "Duplicate Name", "Duplicate Verdict", "Quant Ion (measured)", "Pooled QC CV (%)",
}

MIN_PAIRED = 8      # samples a molecule must be measured in on BOTH sides to be correlated


def number(x) -> float | None:
    """⚠ `float()`, never `str.isdigit()` — an area written 1.2E8 is not a digit string."""
    try:
        return float(str(x).replace(",", ""))
    except (TypeError, ValueError):
        return None


def sample_key(column: str) -> str:
    """`Area: 1_Sample_r1.raw (F1)` and `1_Sample_r1` to the same key.

    ⚠ Order matters: the `(F1)` suffix comes off BEFORE the extension, because the extension is
    not at the end of the string. Reversed, nothing matches and it reads as no shared samples.
    """
    c = re.sub(r"^Area:\s*", "", column).strip()
    c = re.sub(r"\s*\(F\d+\)$", "", c)
    return re.sub(r"\.(raw|mzML|mzXML)$", "", c, flags=re.I)


def load(path: Path) -> dict[str, dict[str, float]]:
    """{sum composition: {sample: area}}. Rows sharing a composition are combined by median."""
    rows = list(csv.DictReader(path.open(errors="replace")))
    if not rows:
        return {}
    cols = [k for k in rows[0] if k and k not in META]
    by: dict[str, list[dict[str, float]]] = {}
    for r in rows:
        name = (r.get("Identification") or "").strip()
        if not name or name.upper().startswith("DECOY"):
            continue
        key = sum_composition(name)
        if key:
            by.setdefault(key, []).append({sample_key(c): (number(r.get(c)) or 0.0) for c in cols})
    return {k: {s: float(np.median([d[s] for d in v])) for s in v[0]} for k, v in by.items()}


def compare(reference: Path, ours: Path, label: str) -> dict | None:
    ref, us = load(reference), load(ours)
    if not ref or not us:
        print(f"  {label}: could not read one side"); return None
    samples = sorted(set(next(iter(ref.values()))) & set(next(iter(us.values()))))
    shared = [k for k in us if k in ref]
    if not samples:
        print(f"  {label}: ⚠ no samples in common — check `sample_key`, not the data"); return None

    r, ratio = [], []
    for k in shared:
        a = np.array([ref[k][s] for s in samples])
        b = np.array([us[k][s] for s in samples])
        m = (a > 0) & (b > 0)
        if m.sum() < MIN_PAIRED:
            continue
        c = np.corrcoef(np.log10(a[m]), np.log10(b[m]))[0, 1]
        if np.isfinite(c):
            r.append(c); ratio.append(float(np.median(b[m] / a[m])))
    r, ratio = np.array(r), np.array(ratio)
    if not r.size:
        print(f"  {label}: no molecule measured in >= {MIN_PAIRED} shared samples"); return None

    print(f"  {label}:  {len(shared)} molecules shared, {len(samples)} samples, {r.size} correlated")
    print(f"     median r {np.median(r):.3f}   r>0.9 {100*np.mean(r>0.9):.1f}%   "
          f"r>0.8 {100*np.mean(r>0.8):.1f}%   r<0.5 {100*np.mean(r<0.5):.1f}%")
    print(f"     area ratio ours/reference: median {np.median(ratio):.2f}  "
          f"IQR {np.percentile(ratio,25):.2f}-{np.percentile(ratio,75):.2f}")
    return {"label": label, "n": int(r.size), "median_r": float(np.median(r)),
            "frac_r_gt_09": float(np.mean(r > 0.9)), "frac_r_lt_05": float(np.mean(r < 0.5)),
            "median_ratio": float(np.median(ratio))}


def main() -> None:
    V = Path("/mnt/datastore/Jair/claudecode/Lipidomics_automation/Validation sets")
    pairs = [
        (V / "Validation_2_Skin/Pos/CD/Final_ResultsPos.csv",
         Path("OmicsPilot_Lipidomics/Pos/Final_Results.csv"), "skin Pos"),
        (V / "Validation_2_Skin/Neg/CD_Neg/Final_ResultsNeg.csv",
         Path("OmicsPilot_Lipidomics/Neg/Final_Results.csv"), "skin Neg"),
        (V / "Validation_1_CKD/Pos/CD/Final_Results.csv",
         Path("OmicsPilot_Lipidomics/CKD_Pos/Final_Results.csv"), "CKD Pos"),
        (V / "Validation_1_CKD/Neg/CD_Neg/Final_Results.csv",
         Path("OmicsPilot_Lipidomics/CKD_Neg/Final_Results.csv"), "CKD Neg"),
    ]
    for ref, ours, label in pairs:
        if ref.exists() and ours.exists():
            compare(ref, ours, label)
        else:
            missing = ref if not ref.exists() else ours
            print(f"  {label}: missing {missing}")


if __name__ == "__main__":
    main()
