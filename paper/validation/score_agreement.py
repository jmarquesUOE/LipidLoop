"""Agreement between our identifications and the ones each deposit published itself.

This is deliberately NOT called recall. Recall implies the published list is truth, and on every
deposit examined so far it is not: it is another software's output, made with another library, at
another evidence threshold. Two pipelines disagreeing tells you they disagree — which one is right
needs the spectrum, not the count.

So the number reported is symmetric and its two halves are kept apart:

  * **we found / they published** — how much of their list we reproduce
  * **they published / we found** — how much of ours they carry

⚠ Most deposits publish NO names at all. Four of the studies staged here (the Waters/ADAP-BIG
benchmark set) deposit `m/z_RT` feature tables by design — they exist to benchmark data PROCESSING,
so an identification column would defeat their purpose. Those score as "not applicable", never as
zero agreement, and the distinction has to survive into the table: a blank cell and a 0% cell say
opposite things about the software.

Matching is on the LSI shorthand the pipeline already emits, at two levels — molecular (chains
resolved on both sides) and species (collapsed to one total) — because crediting a species-level
hit as a molecular one claims a distinction we did not make.
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

sys.path.insert(0, "/home/jair/repos/Lipidomics-automation/src")

import deposits
import refmet
from lipidloop.nomenclature import canonical_class

STAGED = Path("/home/jair/validation_staged")

#: Which analysis directory a dataset's numbers come from, newest code first. Module-level so the
#: results table and the scorer cannot disagree about which run they are describing.
# ⚠ Add every new analysis name HERE when one is created. This list went stale: the studies moved
# to Analysis_v5 and then Analysis_v6, neither of which was listed, so `run_dir` returned "" for
# every dataset, `run_metrics` returned zeros, and the results table could not be regenerated at
# all — it had become a snapshot of v4 runs that no longer exist on disk. A table nobody can
# rebuild is a table nobody can check.
ORDER = ("Analysis_recal_2026-09-04", "Analysis_calibrated", "Analysis_v7", "Analysis_v6", "Analysis_v5",
         "Analysis_v4", "Analysis_v3", "Analysis_ba", "Analysis_v2", "Analysis_centroided",
         "Analysis_perfile", "Analysis_overnight", "Analysis_decoy")


def complete(directory: Path) -> bool:
    """⚠ Has every polarity finished? An in-flight run must never be scored.

    A run that is still going has written `Unfiltered_Results.csv` for the polarity it finished and
    nothing for the one it is working on, so a newest-first search picks it and silently describes
    half a study — which is exactly what happened while a bile-acid re-run was mid-flight. The
    pipeline writes `Final_Results.csv` last, so its presence in every polarity is the signal.
    """
    results = directory / "Results"
    if not results.is_dir():
        return False
    polarities = [d for d in results.iterdir() if d.is_dir()]
    return bool(polarities) and all((d / "Final_Results.csv").exists() for d in polarities)

_CHAIN = re.compile(r"(\d+):(\d+)")
_ADDUCT = re.compile(r"[\[\(]M[+-].*$")
# Spiked internal standards. Workbench marks them `refmet_name = "Standard"`, and the raw names
# carry SPLASH/ISTD/(d7) tags. They are deuterated compounds added by the lab, so counting them in
# the denominator asks us to reproduce something we deliberately never report as an endogenous
# identification — it would depress agreement for behaving correctly.
_ISTD = re.compile(r"ISTD|SPLASH|EquiSPLASH|\(d\d+\)|-d\d+\b", re.I)

#: Classes built on a sphingoid base, where a `d`/`t` chain prefix carries an implied oxygen count.
SPHINGOID = {"Cer", "HexCer", "Hex2Cer", "Hex3Cer", "SM", "SHexCer", "AcylHexCer", "AcylSM",
             "GM1", "GM2", "GM3", "GD1a", "GD1b", "GD2", "GD3", "GT1b", "GQ1b", "PE-Cer", "PI-Cer"}

# Class-name synonyms. Two lists naming the same molecule differently is a notation disagreement,
# not a disagreement about chemistry, and scoring it as a miss inflates our apparent error.
#
# ⚠ HexCer is the RIGHT name and the others are overclaims. Glucosyl- and galactosylceramide are
# stereoisomers that reversed-phase chromatography does not separate and whose MS2 spectra are
# indistinguishable, so an untargeted method cannot tell them apart. Anything reporting `GlcCer`
# from RP-LC/MS2 alone is asserting a distinction it has no evidence for; the honest name is the
# hexosyl one. Same argument one sugar up: Hex2Cer, not LacCer.
# ⚠ ONE synonym map, imported from the package — not a second copy here.
#
# There were two. The scorer carried its own 22 entries while `nomenclature.py` carried 28, and
# although they never contradicted each other on a shared key, the package knew `hex1cer` (which is
# how RefMet spells it) and the scorer did not. Since `refmet.py` keys on the PACKAGE map and
# `canon()` keyed on this one, the species comparison and the RefMet comparison were normalising
# differently — in the same table, on the same row.
#
# ⚠ SE is still deliberately NOT mapped to CE: a steryl ester is not always a cholesteryl ester,
# our libraries carry CE only, and the mapping would credit us with sterol esters we cannot
# identify. That gap is real and belongs in the coverage table.


def _class(name: str) -> str:
    """Canonical class name. Subclass tags survive — Cer[ADS] must not merge with Cer[NP]."""
    return canonical_class(name)


def canon(name: str) -> tuple[str, str]:
    """(molecular key, species key) from a lipid name written in any of the usual dialects.

    Deposits write LSI loosely and inconsistently: `PC(16:0/18:1)` with parentheses from
    LipidSearch, `SM 36:3,O2` with a comma where the standard says semicolon, a trailing
    `|PS 18:0_22:6` alternative after a pipe, an adduct glued on the end. None of those
    differences are chemical, so all are normalised away before comparing.
    """
    s = (name or "").strip()
    if not s:
        return "", ""
    s = s.split("|")[0].strip()          # keep the primary of `PS 40:6|PS 18:0_22:6`
    # `PC P-32:1 or PC O-32:2` is the deposit declining to choose between two ethers, NOT a
    # two-chain molecular species. Read literally it yields a fabricated molecular key, which is
    # how 35 species-level entries were being counted as chain-resolved.
    s = re.split(r"\s+or\s+", s, flags=re.I)[0].strip()
    s = _ADDUCT.sub("", s).strip()
    s = s.replace("(", " ").replace(")", " ")
    s = s.replace(",O", ";O").replace("/", "_")
    s = re.sub(r"\s+", " ", s).strip()

    head, _, tail = s.partition(" ")
    if not tail:
        return "", ""
    cls = _class(head.strip().rstrip(":"))
    chains = _CHAIN.findall(tail)
    if not chains:
        return "", ""
    oxy = ";O" + (re.search(r";O(\d*)", tail).group(1) or "") if ";O" in tail else ""
    # ⚠ `Cer d18:1/16:0` and `Cer 34:1;O2` are THE SAME MOLECULE, and without this they never
    # matched. The `d`/`t` prefix is the old dialect for the sphingoid base's hydroxyls — d = 2,
    # t = 3 — which LSI writes explicitly as `;O2` / `;O3`. Our libraries use the prefix form and
    # the deposits use the LSI form, so every sphingolipid in five deposits was scored as a miss:
    # 106 SM, 89 Cer, 37 SHexCer and 33 HexCer, read as an oxidised-species library gap that would
    # have been "filled" by building entries we already hold.
    #
    # Expand the prefix rather than stripping the suffix: `;O3` on a ceramide is a real
    # phytoceramide and must stay distinct from the `;O2` base form.
    if not oxy and cls.split("[")[0] in SPHINGOID:
        base = re.search(r"\b([dte])(\d+):", tail)
        if base:
            oxy = ";O" + {"d": "2", "t": "3", "e": "1"}[base.group(1).lower()]
    ether = "O-" if re.search(r"\bO-", tail) else ("P-" if re.search(r"\bP-", tail) else "")

    parts = sorted(f"{c}:{d}" for c, d in chains)
    molecular = f"{cls} {ether}{'_'.join(parts)}{oxy}" if len(parts) > 1 else ""
    tc = sum(int(c) for c, _ in chains)
    td = sum(int(d) for _, d in chains)
    species = f"{cls} {ether}{tc}:{td}{oxy}"
    return molecular, species


# Deposits that published a separate list per instrument, and which arm a study directory holds.
# Scoring against the pooled list would compare one instrument to five instruments' union.
#
# ⚠ The MAF names are the deposit's, not ours: its Bruker timsTOF list is `BrukerPASEF` and its
# BAF list is `BrukerDDA`, and the Sciex one is `UC_Davis_Sciex`. Mapping our directory names
# straight through would match no MAF at all and report every vendor arm as publishing nothing.
ARM = {
    "MTBKS222":           "Thermo",
    "MTBKS222_Thermo":    "Thermo",
    "MTBKS222_Agilent":   "Agilent",
    "MTBKS222_Waters":    "Waters",
    # ⚠ Not one Sciex arm but three submissions, each with its own MAF and gradient. The RT check
    # is what caught it: against UC_Davis_Sciex the Tsugawa files sat +4.22 min out with 0% inside
    # 0.2 min; against IMS_Normal they land at -0.01 min and 85%.
    "MTBKS222_Sciex":        "UC_Davis_Sciex",
    "MTBKS222_IMS_Normal":   "IMS_Normal",
    "MTBKS222_IMS_HighMass": "IMS_HighMass",
    "MTBKS222_BrukerTDF": "BrukerPASEF",
    "MTBKS222_BrukerBAF": "BrukerDDA",
}

#: ⚠ Arms deliberately left out of the comparison, and why. Kept as a named set rather than a
#: deletion so a missing row is explained rather than merely absent.
EXCLUDED_ARMS = {
    # Its three NIST SRM 1950 replicates were acquired on three DIFFERENT DAYS — 2017-03-14, 03-15
    # and 03-16 — so they are not a triplicate: a feature has to survive three days of instrument
    # drift to appear in all three. It is also the oldest data in the collection by two years. On
    # v5 it returned 21 identifications against 206-778 for the sibling arms, with zero decoy hits,
    # so it cannot even carry an error rate. A row for it would measure its acquisition design, not
    # its instrument.
    "MTBKS222_UCDavis": "3 replicates across 3 separate days (2017-03-14/15/16); oldest in the set",
    # ⚠ NOT AN INDEPENDENT ARM. Its 24 injections ARE MTBKS222_IMS_Normal + MTBKS222_IMS_HighMass —
    # the same files staged twice under a third name. Counting it as a separate vendor arm
    # double-counts those two in every cross-vendor total. Found by reading the instrument
    # configuration and file identity out of the mzML rather than trusting the folder name.
    "MTBKS222_Sciex": "the same injections as IMS_Normal + IMS_HighMass, staged again under a third name",
    # ⚠ SCOPE, not a defect. PASEF is not conventional DDA and this pipeline is built for
    # conventional DDA: one spectrum per precursor, mobility summed away, HCD Orbitrap libraries.
    # After the converter fix the arm still delivers 13 lipids from 6 injections, with 58.7% of
    # peaks at or above the precursor. Handling PASEF properly means carrying mobility through
    # detection, association and scoring — a different project. See docs/TIMSTOF_CONVERTER_DEFECT.md.
    "MTBKS222_BrukerTDF": "PASEF is out of scope for this pipeline — kept as-is, not quotable",
}

#: A study directory that is one arm of a larger deposit -> the deposit its MAFs live under.
DEPOSIT_OF = {k: "MTBKS222" for k in ARM if k.startswith("MTBKS222_")}


def published(dataset: str) -> list[str]:
    """Names the deposit itself published, or [] when it published none."""
    names: list[str] = []
    deposits.ensure(DEPOSIT_OF.get(dataset, dataset))
    if True:
        try:
            for r in deposits.metabolite_rows(dataset):
                ref = (r.get("refmet_name") or "").strip()
                raw = (r.get("metabolite_name") or "").strip()
                if ref == "Standard" or _ISTD.search(raw):
                    continue
                # RefMet leaves `refmet_name` empty when it cannot map the entry; the submitter's
                # own name is then the only claim on offer.
                names.append(ref or raw)
        except (json.JSONDecodeError, AttributeError):
            pass
    for row in deposits.maf_rows(DEPOSIT_OF.get(dataset, dataset), ARM.get(dataset, "")):
        names.append((row.get("metabolite_identification") or "").strip())
    # UNKNOWN_n placeholders are the deposit saying it did NOT identify the feature.
    return [n for n in names if n and not n.upper().startswith("UNKNOWN")]


def _read(paths) -> list[str]:
    found: list[str] = []
    for h in paths:
        with open(h, errors="replace") as fh:
            for row in csv.DictReader(fh):
                ident = (row.get("Identification") or "").strip()
                if not ident or ident.upper().startswith("DECOY"):
                    continue
                found.append((row.get("Shorthand (LSI)") or ident).strip())
    return found


def ours_delivered(dataset: str) -> tuple[list[str], str]:
    """The identifications an analyst actually RECEIVES, from `Final_Results_Filtered.csv`.

    ⚠ Use this, not `ours`, when scoring against a deposit's published list. A deposit publishes
    its deliverable — the list it chose to stand behind. `ours` reads `Unfiltered_Results.csv`,
    which is the audit trail before the blank, presence, adduct, dimer and isotope filters, so
    scoring it against their curated list compares a pre-filter set with a post-filter one and
    gives us more chances to match. Measured across the studies that publish identifications, that
    asymmetry was worth +4 to +24 percentage points of "agreement":

        MTBLS5163        67.8% -> 44.1%
        ST003077         77.1% -> 55.7%
        ST003514         54.6% -> 41.4%
        MTBKS222_Thermo  73.1% -> 68.3%
        ST004797         28.7% -> 24.4%

    `ours` stays for the FDR, where unfiltered is correct: the decoys are unfiltered too, so that
    comparison is already like-for-like and filtering one side would break it.
    """
    for name in ORDER:
        if not complete(STAGED / dataset / name):
            continue
        hits = sorted(glob.glob(str(STAGED / dataset / name / "Results/*/Final_Results_Filtered.csv")))
        if hits:
            return _read(hits), name
    return [], ""


def ours(dataset: str) -> tuple[list[str], str]:
    found: list[str] = []
    # ⚠ Preference is EXPLICIT and newest-code-first, and the chosen directory is returned.
    #
    # This used to prefer `Analysis_overnight`, which silently re-scored pre-fix runs after a
    # re-processing pass: two datasets kept their old numbers while two showed new ones, in one
    # table, with nothing to indicate the difference. A validation table that mixes code versions
    # is worse than one that is merely out of date, because it looks consistent.
    for name in ORDER:
        if not complete(STAGED / dataset / name):
            continue
        hits = sorted(glob.glob(str(STAGED / dataset / name / "Results/*/Unfiltered_Results.csv")))
        if hits:
            return _read(hits), name
    # ST003514's recall work ran out of the benchmark tree rather than the staged tree.
    for root, pat in ((STAGED, f"{dataset}/Analysis_*/Results/*/Unfiltered_Results.csv"),
                      (Path("/home/jair/validation_bench"),
                       f"{dataset.lower()}_recall/Analysis_acetate/Results/*/Unfiltered_Results.csv")):
        hits = sorted(glob.glob(str(root / pat)))
        if hits:
            return _read(hits), Path(hits[0]).parents[1].name
    return found, "none"


def score(dataset: str) -> dict:
    pub, (mine, source) = published(dataset), ours(dataset)
    if not pub:
        return {"dataset": dataset, "published": 0, "ours": len(set(mine)),
                "source": source, "status": "no names deposited"}

    pm = {m for m, _ in map(canon, pub) if m}
    ps = {s for _, s in map(canon, pub) if s}
    om = {m for m, _ in map(canon, mine) if m}
    os_ = {s for _, s in map(canon, mine) if s}
    if not mine:
        return {"dataset": dataset, "published": len(ps), "ours": 0,
                "source": source, "status": "no run"}

    shared_mol = len(pm & om)
    shared_sp = len(ps & os_)

    # A second comparison under RefMet — the one nomenclature both sides can be written in without
    # either of us choosing it. Workbench applies it to every deposit it holds, so agreement here
    # is not an artefact of the synonym table above.
    #
    # ⚠ RefMet resolves `HexCer` to `GalCer`, committing to a stereoisomer that reversed-phase MS2
    # cannot distinguish. We do not adopt its names, only use them as a shared key; the honest
    # level is its `sub_class`, which says HexCer.
    # ⚠ The RefMet figure is computed over the entries RefMet can NAME, which is not all of them.
    # On one deposit that is 66% — so a higher RefMet agreement there is mostly a smaller
    # denominator, not recovered agreement, and the coverage column has to be read with it.
    # Where coverage is complete the two figures barely differ (100% -> no change at all; 98% ->
    # one point), which says nomenclature was not the main source of disagreement.
    #
    # That RefMet declines to name a third of that deposit's list is itself evidence about the
    # list: those are the implausible chain splits, the plant lipids and the exotic conjugates.
    pr = {refmet.refmet_name(x) for x in ps}
    orr = {refmet.refmet_name(x) for x in os_}
    pr.discard(""); orr.discard("")
    shared_refmet = len(pr & orr)
    of_theirs_refmet = 100 * shared_refmet / len(pr) if pr else 0.0
    return {
        "dataset": dataset,
        "published": len(ps),
        "ours": len(os_),
        "shared_species": shared_sp,
        "shared_molecular": shared_mol,
        "of_theirs": 100 * shared_sp / len(ps) if ps else 0.0,
        "of_ours": 100 * shared_sp / len(os_) if os_ else 0.0,
        "refmet_theirs": len(pr), "refmet_ours": len(orr), "refmet_shared": shared_refmet,
        "of_theirs_refmet": of_theirs_refmet,
        "source": source,
        "status": "ok",
    }


if __name__ == "__main__":
    sets = sorted({p.name for p in STAGED.iterdir() if p.is_dir()})
    rows = [score(d) for d in sets]
    hdr = (f"{'dataset':<16}{'theirs':>8}{'ours':>7}{'shared':>8}{'molec':>7}"
           f"{'of theirs':>11}{'of ours':>9}{'RefMet':>9}{'cov':>6}   source")
    print(hdr); print("-" * len(hdr))
    for r in rows:
        if r["status"] != "ok":
            print(f"{r['dataset']:<16}{r['published']:>8}{r['ours']:>7}{'':>8}{'':>7}{'':>11}{'':>9}   {r['status']}")
        else:
            print(f"{r['dataset']:<16}{r['published']:>8}{r['ours']:>7}{r['shared_species']:>8}"
                  f"{r['shared_molecular']:>7}{r['of_theirs']:>10.1f}%{r['of_ours']:>8.1f}%"
                  f"{r['of_theirs_refmet']:>8.1f}%"
                  f"{100 * r['refmet_theirs'] / max(r['published'], 1):>5.0f}%   {r['source']}")
    out = Path("/home/jair/validation_bench/agreement/agreement.csv")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["dataset", "published", "ours", "shared_species",
                                           "shared_molecular", "of_theirs", "of_ours",
                                           "refmet_theirs", "refmet_ours", "refmet_shared",
                                           "of_theirs_refmet", "source", "status"])
        w.writeheader(); w.writerows(rows)
    print(f"\nwrote {out}")
