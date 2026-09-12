"""The validation table: one row per deposit, with everything a reviewer will ask for.

Every column is read from the deposit or from the run, never assembled by hand — the instrument
from the study's own metadata record, the injection count from the staged mzML, the agreement from
`score_agreement.py`. Retyping any of it into a manuscript by hand is where the errors get in.

Two distinctions the table has to preserve, because collapsing either would misrepresent the work:

  * **injections are not samples.** Every count here is FILES, which include pools, blanks and
    technical replicates. Where the deposit states a biological n it is carried separately.
  * **"no names deposited" is not 0% agreement.** Six of these deposits publish `m/z_RT` feature
    tables with no identification column — four of them by design, as data-PROCESSING benchmarks.
    A blank cell and a 0% cell make opposite claims about the software.
"""
from __future__ import annotations

# Run from anywhere: the sibling `deposits` module is next to this file, not on the path.
import sys, pathlib; sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))

import csv
import glob
import json
import re
import sys
from pathlib import Path

# ⚠ Import the SIBLING scorer, never a copy in the output directory.
#
# This line used to read `sys.path.insert(0, "/home/jair/validation_bench/agreement")`, putting the
# output folder ahead of this one, where a 2026-08-21 copy of score_agreement.py still sat. The
# table builder silently used that copy — 164 lines behind, with no Analysis_v3 preference, no
# RefMet and no per-arm MAF support — so it reported MTBLS5163 at 66.1% where the current scorer
# says 71.2%, and showed MTBKS222 as publishing no identifications at all.
#
# The manuscript's own blocking list says no table may mix code vintages. The tool that builds the
# table was mixing them.
from score_agreement import score, ARM  # noqa: E402

import deposits

STAGED = Path("/home/jair/validation_staged")
OUT = Path(__file__).resolve().parent

# ⚠ Species is READ from the deposit, never typed here.
#
# An earlier version of this file hardcoded these, and hardcoded MSV000094718 as Mus musculus. It
# is Bos taurus — bovine liver, NCBI TaxID 9913, stated plainly in the MassIVE PROXI record. A
# wrong species in a validation table is the kind of error a reviewer finds and a reader cannot
# recover from, and it got in because the field looked too small to be worth sourcing.
#
# What remains here is instrument and matrix for repositories with no machine-readable field for
# them, each carrying the page it was read from. Species is not in this table at all.
NON_WORKBENCH = {
    "MTBLS5163":    ("Bruker maXis II QTOF", "Dendritic cells"),
    "MTBLS2016":    ("Bruker maXis QTOF", "Whole organism"),
    "MSV000094718": ("Agilent QTOF", "Liver + solvent blanks"),
    "MSV000095868": ("Agilent QTOF", "Gut microbiome / caecal"),
    "Skin_QEplus":  ("Q Exactive Plus", "Skin organoid"),
    # MetaboLights Kickstarter: NIST SRM 1950 plasma on five vendors. Only the QExactive Plus arm
    # has been processed — the other three staged arms are unconverted, and 18 of the 26 Thermo
    # .raw are 0-byte failed downloads. Repository is MetaboLights, not in-house; the MTBKS prefix
    # simply is not what the ST/MTBLS/MSV test looks for.
    "MTBKS222":     ("Thermo Q Exactive Plus", "Blood (plasma), NIST SRM 1950"),
}

# ⚠ Several deposits are ONE study split across accessions. Counting them separately overstates n.
#
#   ST004503 + ST004626   one benchmark, reversed phase, positive and negative
#   ST004650 + ST004651   the same benchmark on HILIC, positive and negative
#   ST000987 + ST000991   parts V and IX of one nine-platform comparison, same 153 samples
#
# All four of the first group are titled "Untargeted METABOLOMICS", and two of them are HILIC. A
# lipid library against HILIC metabolomics is not a lipidomics validation, whatever it scores, so
# the scope is carried as its own column rather than left for a reader to infer from the title.
STUDY_GROUP = {
    "ST004503": "CSU-benchmark-RP", "ST004626": "CSU-benchmark-RP",
    "ST004650": "CSU-benchmark-HILIC", "ST004651": "CSU-benchmark-HILIC",
    "ST000987": "nine-platform", "ST000991": "nine-platform",
}


def species_of(dataset: str, mwtab: dict, summary: dict) -> tuple[str, str]:
    """(species, where it was read from). Blank when no source states it — never guessed."""
    v = mwtab.get("SU:SUBJECT_SPECIES", "")
    if v:
        return v, "mwtab:SU:SUBJECT_SPECIES"
    v = summary.get("species", "")
    if v:
        return v, "workbench:summary"
    for suffix, src in (("json", "massive:proxi"), ("ws.json", "metabolights:ws")):
        d = deposits.record(dataset, suffix)
        if not d:
            continue
        for group in (d.get("species") or []):
            for entry in (group if isinstance(group, list) else [group]):
                if isinstance(entry, dict) and "scientific name" in str(entry.get("name", "")):
                    return str(entry.get("value") or ""), src
        for org in ((d.get("content") or {}).get("organism") or []):
            n = (org.get("organismName") or "").strip()
            if n and n.lower() != "blank":
                return n, src
    return "", "not stated"


def assay_scope(mwtab: dict, summary: dict) -> str:
    """lipidomics / metabolomics / methods, from the deposit's own title and column chemistry."""
    title = (mwtab.get("ST:STUDY_TITLE", "") or summary.get("study_title", "")).lower()
    chrom = (mwtab.get("CH:CHROMATOGRAPHY_TYPE", "") or "").upper()
    scope = "lipidomics" if "lipidom" in title else (
        "metabolomics" if "metabolom" in title else "unstated")
    if "HILIC" in chrom:
        scope += " (HILIC)"
    return scope


mwtab_fields = deposits.mwtab_fields
summary = deposits.summary


def deposited_rt(dataset: str) -> str:
    """Does the deposit publish a retention time with its identifications?

    This is the column that decides whether an agreement figure can be checked at PEAK level or
    only at name level, and the answer splits by REPOSITORY rather than by submitter. MetaboLights
    MAF files carry a `retention_time` column; the Metabolomics Workbench metabolite records have
    no equivalent field, so every Workbench study here publishes none — 1,066 species in ST004797,
    415 in ST003077, and not one retention time among them.

    Matching on name alone invents agreement: two pipelines can both report `PC 34:1` while
    pointing at different peaks, and without retention there is no way to tell.
    """
    n = total = 0
    for row in deposits.maf_rows(dataset, ARM.get(dataset, "")):
        total += 1
        if (row.get("retention_time") or "").strip():
            n += 1
    for row in deposits.metabolite_rows(dataset):
        total += 1
        if any((row.get(k) or "").strip() for k in ("retention_time", "ri", "RT")):
            n += 1
    if not total:
        return ""
    if n == total:
        return "all"
    return f"{100 * n // total}%" if n else "none"


def rows() -> list[dict]:
    out = []
    for d in sorted(x.name for x in STAGED.iterdir() if x.is_dir()):
        npos = len(glob.glob(str(STAGED / d / "Pos" / "*.mzML")))
        nneg = len(glob.glob(str(STAGED / d / "Neg" / "*.mzML")))
        f, s = mwtab_fields(d), summary(d)
        inst = f.get("MS:INSTRUMENT_NAME", "")
        stype = f.get("CO:SAMPLE_TYPE", "")

        if d in NON_WORKBENCH:
            i2, t2 = NON_WORKBENCH[d]
            inst, stype = inst or i2, stype or t2
        species, species_src = species_of(d, f, s)
        sc = score(d)
        out.append({
            "dataset": d,
            "repository": ("Workbench" if d.startswith("ST") else
                           "MetaboLights" if d.startswith(("MTBLS", "MTBKS")) else
                           "MassIVE" if d.startswith("MSV") else "in-house"),
            "instrument": re.sub(r"\s*(Hybrid Quadrupole-Orbitrap )?Mass Spectrometer$", "", inst),
            "type": f.get("MS:INSTRUMENT_TYPE", "") or ("Orbitrap" if "Exactive" in inst else "QTOF"),
            "pos": npos, "neg": nneg, "injections": npos + nneg,
            "deposit_n": s.get("number_of_samples", ""),
            "sample_type": stype, "species": species, "species_source": species_src,
            "study_group": STUDY_GROUP.get(d, d), "assay_scope": assay_scope(f, s),
            "our_ids": sc.get("ours", 0),
            "their_names": sc.get("published", 0),
            "shared": sc.get("shared_species", ""),
            "of_theirs": round(sc["of_theirs"], 1) if sc.get("status") == "ok" else "",
            "of_ours": round(sc["of_ours"], 1) if sc.get("status") == "ok" else "",
            "status": sc.get("status", ""),
            "deposited_rt": deposited_rt(d),
        })
    return out


def markdown(rs: list[dict]) -> str:
    def agree(r):
        if r["status"] == "no names deposited":
            return "— *no IDs deposited*"
        if r["status"] == "no run":
            return "*not yet run*"
        return f"**{r['of_theirs']}%** ({r['shared']}/{r['their_names']})"

    def rt(r):
        return {"all": "**all**", "none": "none", "": "—"}.get(r["deposited_rt"], r["deposited_rt"])

    lines = ["| dataset | repository | instrument | inj. (+/−) | sample type | our IDs "
             "| agreement with deposit | deposit publishes RT |",
             "|---|---|---|---|---|---|---|---|"]
    for r in sorted(rs, key=lambda x: (x["status"] != "ok", -x["our_ids"])):
        inj = f"{r['injections']} ({r['pos']}/{r['neg']})" if r["injections"] else "—"
        st = (r["sample_type"] or "—")
        if r["species"]:
            st += f" · *{r['species']}*"
        lines.append(f"| {r['dataset']} | {r['repository']} | {r['instrument'] or '—'} | {inj} | "
                     f"{st} | {r['our_ids'] or '—'} | {agree(r)} | {rt(r)} |")
    return "\n".join(lines)


if __name__ == "__main__":
    for d in sorted(x.name for x in STAGED.iterdir() if x.is_dir()):
        deposits.ensure(d)          # fetch anything missing; a failure leaves the field blank
    rs = rows()
    OUT.mkdir(parents=True, exist_ok=True)
    with (OUT / "validation_table.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rs[0].keys()))
        w.writeheader(); w.writerows(rs)
    md = markdown(rs)
    (OUT / "validation_table.md").write_text(md + "\n")
    print(md)
    ok = [r for r in rs if r["status"] == "ok"]
    print(f"\n{len(ok)} deposits with published names; "
          f"mean agreement {sum(r['of_theirs'] for r in ok)/len(ok):.1f}%" if ok else "")
    print(f"total injections processed: {sum(r['injections'] for r in rs)}")
    with_rt = [r for r in rs if r["deposited_rt"] == "all"]
    named = [r for r in rs if r["their_names"]]
    print(f"deposits publishing retention times for every entry: {len(with_rt)} of {len(named)} "
          f"that publish identifications ({', '.join(sorted(r['dataset'] for r in with_rt))})")
