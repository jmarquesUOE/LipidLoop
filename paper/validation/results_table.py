"""The results table: what each deposit reported, what we found, and how confident either is.

One row per dataset. Every number is read from a file — the deposit's own record, or a run's own
output — and where a number is not available the cell says so rather than being left to look like
a zero.

⚠ **Three distinctions this table exists to keep.**

*Detected is not identified.* A compound group is a peak the software grouped across injections; an
identification is a name attached to one. The ratio between them is the interesting quantity, and
collapsing them hides the failure mode where a run detects 51,388 things and names one.

*Injections are not samples, and searched is not staged.* Acquisition is a property of the file:
several deposits mix DDA with SWATH or MS1-only in one submission, so the number of files SEARCHED
is smaller than the number staged, and both differ from the number of animals.

*A missing comparison is not a failed one.* Six deposits publish no identifications at all — four
of them by design, as data-processing benchmarks. Those cells are blank, never 0%.
"""
from __future__ import annotations

import csv
import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import deposits                                            # noqa: E402
from score_agreement import canon, complete, ours, ours_delivered, published, ORDER  # noqa: E402
from lipidloop.calibrate import acquisition               # noqa: E402

STAGED = Path("/home/jair/validation_staged")

PLATFORM = {
    "ST004797": "Sciex ZenoTOF 7600", "ST003514": "Agilent 6545 QTOF",
    "ST000991": "Sciex X500R QTOF", "ST003052": "Thermo Q Exactive HF",
    "MTBLS5163": "Bruker maXis II", "ST002705": "Bruker micrOTOF-Q II",
    "MSV000095868": "Agilent QTOF", "MSV000094718": "Agilent QTOF",
}
# Datasets in the lipidomics DDA set. Others are excluded with a reason, listed under the table.
EXCLUDED = {
    "ST000987": "all-ion (in-source CID), no precursor selection",
    "MTBLS2016": "no MS2 in any of 147 files, despite being described as DDA",
    "MSV000094718": "methods deposit (bovine liver, vendor contaminants), not a study",
    "ST004503": "metabolomics benchmark", "ST004626": "metabolomics benchmark",
    "ST004650": "metabolomics benchmark, HILIC", "ST004651": "metabolomics benchmark, HILIC",
    # ⚠ Held out on a CONVERTER fault, not on its data. Its forward dot product has a median of 8
    # against Thermo's 926 while its reverse dot is 1000 — every library peak present and far more
    # besides, which is a conversion producing profile-like spectra, not a bad instrument. Without
    # this line it appeared in the table with 0 groups, 0 identifications and "0.0% agreement",
    # which reads as a measured result and is not one.
    "MTBKS222_BrukerTDF": "timsTOF conversion unresolved — forward dot median 8 vs Thermo 926",
    "Skin_QEplus": "in-house, held separately",
    # ⚠ Was missing this line while Skin_QEplus had it — both are in-house with no public
    # identifications (validation_table.md already labels this one "in-house"), but only skin was
    # excluded from the public-corpus table. Rat_Heart_Lumos was counted inside "15 studies" in the
    # abstract/§2.3 by omission, double-dipping it as both a primary validation set and a corpus
    # member. It belongs with skin, held out of the public corpus the same way.
    "Rat_Heart_Lumos": "in-house, held separately",
    # Same rat-heart data, anonymised and re-staged (raw filenames stripped of disease-group/PI
    # identifiers) as a third primary validation set, not a public corpus member. Excluded for the
    # same reason as Rat_Heart_Lumos above, under its new name.
    "Rat_Heart": "in-house, held separately (anonymised primary validation set, not a public deposit)",
    # Third primary validation set (mouse pancreatic tumour, positive mode only), not a public
    # deposit — same category as Rat_Heart/Skin_QEplus, must not appear in the public-corpus table.
    # ⚠ Was keyed "PDAC_KPC" — no such directory exists, the real staged name is "PDAC_KPC_Set3"
    # (Table 2's "Set 3"). The dead key never matched, so this in-house set leaked into the public
    # corpus as a 16th study once it got a recognized complete run (2026-09-08). Same double-dip
    # class as the Rat_Heart_Lumos incident below — caught the same way, by checking the table's
    # own row count against its "Fifteen studies" prose after a fresh run.
    "PDAC_KPC_Set3": "in-house, held separately (primary validation set, not a public deposit)",
    # ⚠ NOT AN INDEPENDENT ARM — score_agreement.py's EXCLUDED_ARMS already knew this
    # (MTBKS222_Sciex's ARM comment), but this table's own EXCLUDED never got the matching entry,
    # so RESULTS_TABLE.md was triple-counting the same 24 injections as three separate studies.
    # Verified by reading instrument configuration and file identity out of the mzML, not the
    # folder name. This is the manual "17→15 studies" dedup the abstract has been doing by hand —
    # encoding it here makes it automatic instead of a step someone has to remember.
    "MTBKS222_Sciex": "the same injections as IMS_Normal + IMS_HighMass, staged again under a third name",
}


def run_dir(dataset: str) -> str:
    """The newest COMPLETE run. ⚠ An in-flight run must never be scored.

    `Analysis_ba` was mid-flight when this was first written — positive mode finished, negative
    mode still detecting features — and being newest it was picked, so the table silently described
    half a study. A directory counts as complete only when every polarity present has written its
    `Final_Results.csv`, which the pipeline does last.
    """
    for name in ORDER:
        # ⚠ Directories only. `Path(...) / "Results/*/"` drops the trailing slash, so the glob
        # also matched `QC_report.html` and `Combined_Filtered_Normalised.csv`, and a file has no
        # `Final_Results.csv` inside it — so every complete run was judged incomplete.
        if complete(STAGED / dataset / name):
            return name
    return ""


def run_metrics(dataset: str) -> dict:
    """Compound groups, identifications and decoys, read from the run this table scores."""
    d = run_dir(dataset)
    groups = ids = decoys = 0
    for f in glob.glob(str(STAGED / dataset / d / "Results/*/Unfiltered_Results.csv")) if d else []:
        for row in csv.DictReader(open(f, errors="replace")):
            groups += 1
            name = (row.get("Identification") or "").strip()
            if not name:
                continue
            if name.upper().startswith("DECOY"):
                decoys += 1
            else:
                ids += 1
    return {"analysis": d, "groups": groups, "id_rows": ids, "decoys": decoys}


def files(dataset: str) -> tuple[int, int]:
    """(staged injections, files an acquisition check finds searchable).

    ⚠ Count `.raw` as well as `.mzML`. This globbed only mzML, so ST003077 — staged as 14 Thermo
    `.raw` — returned (0, 0), and the study was written off in EXCLUDED as "not staged". It was
    staged the whole time, and it delivers 545 rows across 27 lipid classes at 0.12% decoy FDR.
    The pipeline converts `.raw` itself; only this table could not see them. A hand-typed exclusion
    then outlived the measurement that prompted it, which is exactly what README.md says this
    file must not do.
    """
    fs = sorted(STAGED.joinpath(dataset).glob("*/*.mzML"))
    if not fs:
        raw = sorted(STAGED.joinpath(dataset).glob("*/*.raw"))
        if raw:
            # `acquisition()` parses mzML, so classify the conversions rather than the .raw. They
            # live in <study>/mzml/<polarity>/ — one per input, written on the first run.
            converted = sorted(STAGED.joinpath(dataset).glob("mzml/*/*.mzML"))
            searchable = sum(1 for f in converted if acquisition(f) in ("DDA", "unknown"))
            return len(raw), searchable
        return 0, 0
    # ⚠ Classify EVERY file, not a sample. Evenly spaced sampling reported 0 searchable files for a
    # study with 5 usable ones, because its DDA injections are 6 of 159 and clustered by name
    # (`QC_DDA_*`). The searchable count is exactly the number a sample is worst at estimating —
    # it is small, and it is not spread through the run.
    return len(fs), sum(1 for f in fs if acquisition(f) in ("DDA", "unknown"))


# Exclusions that are FACTUAL CLAIMS about the data, mapped to the check that proves them. The
# rest of EXCLUDED are editorial judgements ("metabolomics benchmark", "not a study", "held
# separately") which no measurement can make and which are legitimately typed by hand.
#
# ⚠ The distinction is the whole point. "not staged" was a fact, it went stale, and nothing
# noticed for months because a reason string is never re-read. A judgement cannot go stale; a
# fact can, and must therefore be re-checked every time the table is built.
CHECKABLE = {"not staged": lambda d: files(d)[0] == 0}


def audit_exclusions() -> list[str]:
    """Complaints about exclusions whose stated reason is no longer true."""
    problems = []
    for dataset, reason in sorted(EXCLUDED.items()):
        check = CHECKABLE.get(reason)
        if check and not check(dataset):
            problems.append(f"{dataset}: excluded as {reason!r}, but that is no longer true")
    return problems


def main() -> None:
    for complaint in audit_exclusions():
        print(f"⚠ STALE EXCLUSION — {complaint}", file=sys.stderr)
    sets = [d for d in sorted(x.name for x in STAGED.iterdir() if x.is_dir())
            if d not in EXCLUDED]
    rows = []
    for d in sets:
        deposits.ensure(d)
        n_files, n_search = files(d)
        m = run_metrics(d)
        pub = published(d)
        searched, _ = ours(d)                  # audit trail — the recall ceiling, and the FDR base
        mine, _ = ours_delivered(d)            # what an analyst receives
        ps = {s for _, s in map(canon, pub) if s}
        os_ = {s for _, s in map(canon, mine) if s}
        searched_species = {s for _, s in map(canon, searched) if s}
        # ⚠ Agreement is scored on the DELIVERED set. A deposit publishes its deliverable, so
        # scoring our audit trail against their curated list is not like-for-like.
        shared = len(ps & os_)
        rows.append({
            "dataset": d, "platform": PLATFORM.get(d, "—"),
            "files": n_files, "searched": n_search,
            "groups": m["groups"], "id_rows": m["id_rows"],
            "species": len(os_), "searched_species": len(searched_species),
            "decoys": m["decoys"],
            "fdr": 100 * m["decoys"] / max(m["id_rows"] + m["decoys"], 1),
            "theirs": len(ps), "shared": shared,
            "agree": 100 * shared / len(ps) if ps else None,
            "analysis": m["analysis"],
        })

    hdr = ("| dataset | platform | inj. (searched) | compound groups | ID rows "
           "| species delivered | species searched | decoy FDR | deposit species | shared "
           "| agreement |")
    print(hdr)
    print("|" + "---|" * 11)
    for r in sorted(rows, key=lambda x: -(x["agree"] or -1)):
        agree = f"**{r['agree']:.1f}%**" if r["agree"] is not None else "— *none deposited*"
        theirs = r["theirs"] or "—"
        shared = r["shared"] if r["theirs"] else "—"
        print(f"| {r['dataset']} | {r['platform']} | {r['files']} ({r['searched']}) "
              f"| {r['groups']:,} | {r['id_rows']:,} | {r['species']:,} "
              f"| {r['searched_species']:,} "
              f"| {r['fdr']:.2f}% | {theirs} | {shared} | {agree} |")

    # ⚠ Every row must come from the SAME analysis. `run_dir` falls back to the newest COMPLETE
    # run per dataset, so building this while a batch is in flight silently mixes versions —
    # finished studies score v6, unfinished ones score v5 — and the totals then describe a run
    # that never happened. Report it rather than printing a plausible-looking table.
    versions = {r["analysis"] for r in rows if r["analysis"]}
    if len(versions) > 1:
        print(f"⚠ MIXED ANALYSIS VERSIONS — {', '.join(sorted(versions))}. These totals combine "
              f"runs from different code. Rebuild once one batch has finished.", file=sys.stderr)
        for r in sorted(rows, key=lambda x: x["analysis"]):
            print(f"    {r['dataset']:<18} {r['analysis'] or '(no complete run)'}", file=sys.stderr)
    missing = [r["dataset"] for r in rows if not r["analysis"]]
    if missing:
        print(f"⚠ NO COMPLETE RUN for: {', '.join(missing)} — their numbers are zeros, not results.",
              file=sys.stderr)

    tot = lambda k: sum(r[k] for r in rows)
    print(f"\n**Totals** — {len(rows)} studies, {tot('files'):,} injections "
          f"({tot('searched'):,} searchable), {tot('groups'):,} compound groups, "
          f"{tot('id_rows'):,} identification rows, {tot('decoys')} decoys "
          f"({100*tot('decoys')/max(tot('id_rows')+tot('decoys'),1):.2f}% FDR).")
    with_pub = [r for r in rows if r["agree"] is not None]
    if with_pub:
        print(f"{len(with_pub)} deposits publish identifications: "
              f"{sum(r['theirs'] for r in with_pub):,} species, "
              f"agreement {min(r['agree'] for r in with_pub):.1f}–"
              f"{max(r['agree'] for r in with_pub):.1f}%.")
    print("\n**Excluded**, with reason:\n")
    for k, v in sorted(EXCLUDED.items()):
        print(f"- `{k}` — {v}")


if __name__ == "__main__":
    main()
