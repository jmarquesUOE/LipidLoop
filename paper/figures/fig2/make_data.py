"""Figure 2 data: exactness per file (A), recovery of the reference (B), Set 2 score pairs (C).

Writes `data/*.csv` beside this script so the panel scripts never touch the datastore. Re-run this
when an analysis directory or a reference table changes.

Sources
  Set 1 exactness   scripts/validate_qc.py over every Set 1 file with the reference's three
                    libraries (see manuscript/SUPPORTING_INFORMATION.md, Table S4). The per-file
                    totals are re-derived here from the same inputs.
  Set 2 exactness   LipiDex's Peak Finder associated-spectra table matched to our per-file search
                    output (screen-off run) by file, precursor and retention.
  Recovery          scripts/compare_to_reference.py logic, unfiltered output vs CD Final_Results.
"""
from __future__ import annotations

import bisect
import csv
import glob
import io
import contextlib
import os
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

SET1_REF = Path("/mnt/datastore/Jair/claudecode/Lipidomics_automation/Validation_Datasets/20250423_PDAC_KPC_Set_3")
SET1_RUN = Path("/mnt/work/lipidomics_work/PDAC_KPC_Set3/Analysis_v8_decoy_2026-09-10")
SET2_REF = Path("/mnt/datastore/Jair/claudecode/Lipidomics_automation/Validation sets/NIST/20260911_NIST_v4_Serum")
SET2_RUN = Path("/mnt/work/lipidomics_work/NIST_v4_2026-09-11/Analysis_v8_2026-09-11")
SET2_NOSCREEN = Path("/mnt/work/lipidomics_work/NIST_v4_2026-09-11/Analysis_v8_noscreen_2026-09-11")
REF_LIBRARIES = ["LipidBlast_Formic", "LipiDex_HCD_Formic", "LipiDex_HCD_Hydroxy"]
# Set 1 negative-mode LipiDex exports that are genuine (the other 19 are copies of the positive files)
SET1_NEG_REAL = {"Acat1_KO_C1_01", "Acat1_KO_C1_03"}


def strip(name: str) -> str:
    return re.sub(r"\s*\[.*$", "", name.strip().rstrip(";")).strip()


def set1_exactness() -> list[dict]:
    import validate_qc as vq
    vq.DATA = SET1_REF
    rows = []
    for pol in ("Pos", "Neg"):
        for mgf in sorted((SET1_REF / pol / "mgf").glob("*.mgf")):
            stem = mgf.stem
            if pol == "Neg" and stem not in SET1_NEG_REAL:
                continue
            sys.argv = ["validate_qc.py", "--polarity", f"{pol}/mgf", "--file", stem, "--show", "0",
                        "--libraries", *REF_LIBRARIES]
            buf = io.StringIO()
            with contextlib.redirect_stdout(buf):
                try:
                    vq.main()
                except SystemExit:
                    pass
            t = buf.getvalue()
            g = lambda p: int(re.search(p, t).group(1))
            rows.append({"set": "Set 1", "polarity": pol, "file": stem,
                         "spectra": g(r"reference rows in scope: (\d+)"),
                         "identical": g(r"identification identical: (\d+)/"),
                         "tie": g(r"tied on score:\s+(\d+)"),
                         "dot_identical": g(r"dot product identical:\s+(\d+)/")})
            print("set 1", pol, stem, rows[-1]["spectra"])
    return rows


def set2_exactness() -> tuple[list[dict], list[dict]]:
    summary, pairs = [], []
    for pol in ("Pos", "Neg"):
        ref = list(csv.DictReader((SET2_REF / pol / "CD" / "Associated_Spectra.csv").open(newline="")))
        ours = {}
        for f in glob.glob(str(SET2_NOSCREEN / "Results" / pol / "search" / "*_Results.csv")):
            stem = os.path.basename(f).replace("_Results.csv", "")
            rs = [r for r in csv.DictReader(open(f)) if r["Rank"] == "1"]
            rs.sort(key=lambda r: float(r["Precursor Mass"]))
            ours[stem] = (rs, [float(r["Precursor Mass"]) for r in rs])
        per = {}
        for r in ref:
            stem = r["Sample"].split(" (")[0].replace(".raw", "")
            d = per.setdefault(stem, {"set": "Set 2", "polarity": pol, "file": stem, "spectra": 0,
                                      "identical": 0, "tie": 0, "dot_identical": 0})
            d["spectra"] += 1
            rs, masses = ours[stem]
            prec, rt = float(r["Precursor"]), float(r["Retention"])
            c = [x for x in rs[bisect.bisect_left(masses, prec - 0.003):bisect.bisect_right(masses, prec + 0.003)]
                 if abs(float(x["Retention Time (min)"]) - rt) <= 0.003]
            if not c:
                continue
            c.sort(key=lambda x: abs(float(x["Retention Time (min)"]) - rt))
            o = c[0]
            same_name = strip(o["Identification"]) == strip(r["Name"])
            same_dot = int(float(o["Dot Product"])) == int(float(r["DotProduct"]))
            d["identical"] += same_name
            d["tie"] += (not same_name) and same_dot
            d["dot_identical"] += same_dot
            pairs.append({"polarity": pol, "file": stem, "lipidex_dot": int(float(r["DotProduct"])),
                          "ours_dot": int(float(o["Dot Product"])), "name": strip(r["Name"])})
        summary += list(per.values())
    return summary, pairs


def recovery() -> list[dict]:
    from lipidex_py.peaks import sum_composition
    out = []
    for label, run, ref, std in (("Set 1", SET1_RUN, SET1_REF, "D5TG"), ("Set 2", SET2_RUN, SET2_REF, "D5TG")):
        for pol in ("Pos", "Neg"):
            def names(path):
                return {sum_composition(r["Identification"].strip()) for r in csv.DictReader(open(path))
                        if r.get("Identification", "").strip() and not r["Identification"].upper().startswith("DECOY")}
            ours = names(run / "Results" / pol / "Unfiltered_Results.csv")
            other = names(run / "Results" / ("Neg" if pol == "Pos" else "Pos") / "Unfiltered_Results.csv")
            theirs = names(ref / pol / "CD" / "Final_Results.csv")
            rec = theirs & ours
            missed = theirs - ours
            standard = {n for n in missed if n.startswith(std)}
            in_other = {n for n in missed - standard if n in other}
            absent = missed - standard - in_other
            out.append({"set": label, "polarity": pol, "reference": len(theirs), "recovered": len(rec),
                        "other_polarity": len(in_other), "absent": len(absent), "standard": len(standard),
                        "absent_names": "; ".join(sorted(absent)), "other_names": "; ".join(sorted(in_other))})
            print(label, pol, out[-1])
    return out


def write(name: str, rows: list[dict]) -> None:
    with (DATA / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


if __name__ == "__main__":
    s1 = set1_exactness()
    s2, pairs = set2_exactness()
    write("exactness.csv", s1 + s2)
    write("set2_dot_pairs.csv", pairs)
    write("recovery.csv", recovery())
    print("written to", DATA)
