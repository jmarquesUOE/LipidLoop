"""Co-elution filter: a delivered lipid that is really another lipid's ion (v1.1 item 1).

The peak finder's adduct / isotope / in-source sweep (`peakfinder._filter_pair`) removes the member of
a pair that has NO name. When both members are named nothing ties them together, and when the weak
one is named and the strong one is not, the parent goes and the derivative is delivered. Measured on
the review's ground truth and the Block B runs: M+2 isotopes delivered as one-fewer-double-bond
species, [M+Na]+ / [M+K]+ as two-carbon-longer lipids, PC [M-CH3]- as PE-NMe2, alkylammonium adducts
as +2C TG.

Two independent tests run on the finished tables of both polarities, after the polarity rule. A
delivered row is removed when EITHER labels it -- their union is what was validated:

  ion family     anchor-based. Every named row predicts the ions its molecule must also produce
                 (adducts, class in-source products, dimers, M+1..M+3, cross-polarity choline ions);
                 a co-eluting row at one of them is labelled. Sodium relations need positive MS2
                 evidence of sodium (146.98 / M-59 / M-43); the alkylammonium rule needs a run-level
                 gate (the run must show the adduct). Ported unchanged from the review's frozen v3
                 (ion_family_FROZEN_v3_2026-10-20.py, sha256 a387a02898c4c9c8).
  shape test     row-based. A stronger co-eluting partner at an adduct-swap / water-loss / M+2 /
                 demethylation shift whose XIC has the same SHAPE over the target's peak (r >= 0.9),
                 the same apex, and a derivative-sized area ratio. Ported from the frozen parent check
                 v1 (parent_check_FROZEN_v1_2026-11-04.py, sha256 007063cd5c9e1843) with ONE change
                 (v1b, Jair 2026-11-04): the apex may be up to `apex_tol_high_r` (0.075 min) off when
                 r >= `high_r` (0.95), else `apex_tol` (0.05). On the QE, Na adducts apex 0.01-0.06
                 min before their [M+H]+; at a fixed 0.05 PC 34:3 [M+Na]+ (offset 0.051, r 0.96)
                 survived in one arm and not another. Ground truth unchanged (one extra 'unexplained'
                 row removed on MTBLS394 QE). `apex_tol_high_r = apex_tol` restores v1 exactly.

Rules learned in the review and kept here:
  - MS2 never rescues or vetoes a Na/K relation via m/z 184: the 1 Da isolation window co-isolates
    the [M+H]+ neighbour, so Na adducts of PC show 184 as base peak.
  - a generic water loss from a named diacyl parent is the plasmalogen formula (PC 38:5 - H2O =
    PC P-38:5): only unnamed or hydroxy / lyso / Cer / CE / FA parents lose water.
  - co-elution within 0.05 min is not enough: PDAC ether TG peak ~2.5 s after TG n:(DB-1) and their
    ratio climbs ~10x across the peak (separate compounds), hence the shape test.
  - mzML retention time can be in seconds (LipidMS deposit): the unit is read from the file.
  - decoys are targets, never anchors.

Nothing is deleted from `Final_Results.csv`: it gains `Ion Of / Ion Relation / Ion ppm / Ion dRT /
Ion r / Ion Test`. `Final_Results_Filtered.csv` and the normalised table drop the labelled rows;
`Unfiltered_Results.csv` records why; `coelution_filter.json` holds the counts.
"""
from __future__ import annotations

import csv
import json
import math
import re
import statistics
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, asdict
from pathlib import Path

import numpy as np

POLARITIES = ("Pos", "Neg")
FILTERED = "Final_Results_Filtered.csv"
NORMALISED = "Final_Results_Filtered_median_Normalised.csv"
FINAL = "Final_Results.csv"
UNFILTERED = "Unfiltered_Results.csv"
LABEL_COLUMNS = ["Ion Of", "Ion Relation", "Ion ppm", "Ion dRT", "Ion r", "Ion Test"]


@dataclass
class CoelutionFilter:
    enabled: bool = True
    rt_tol: float = 0.05            # ion family: co-elution window for generic / class relations
    max_ms2_files: int = 12         # ion family: mzML files indexed for the MS2 evidence
    apex_tol: float = 0.05          # shape test
    apex_tol_high_r: float = 0.075  # shape test, when r >= high_r (v1b); = apex_tol for v1
    high_r: float = 0.95


# ════════════════════════════════════════════════════════════════════════════════════════════════
# ion family (frozen v3, unchanged logic)
# ════════════════════════════════════════════════════════════════════════════════════════════════

PPM = 5.0            # after removing the run's median offset
ISO_PPM = 4.0        # 13C spacing vs 2H spacing differ by ~11 ppm at m/z 800; keep below half
ISO_RT = 0.04        # min; isotopes co-elute exactly
ISO_SLACK = 1.5
ISO_RT_WIDE = 0.25   # a saturated parent's apex drifts (PC 40:6 at 9.86 vs its M+2 at 10.05, LipidMS);
                     # allowed only when MASS decides: <=3 ppm from the isotope AND >=4 ppm worse on its own name
ISO = 1.0033548378

# adduct name as the pipeline writes it -> mass shift from the neutral
ADDUCT = {"[M+H]+": 1.007276, "[M+Na]+": 22.989218, "[M+NH4]+": 18.033823, "[M+K]+": 38.963158,
          "[M+H-H2O]+": -17.003288, "[M]+": 0.0, "[M+Li]+": 7.016003,
          "[M-H]-": -1.007276, "[M+FA-H]-": 44.998201, "[M+HCOO]-": 44.998201, "[M+Ac-H]-": 59.013851,
          "[M+CH3COO]-": 59.013851, "[M+Cl]-": 34.969402, "[M-CH3]-": -15.023475, "[M-2H]2-": None}

# ions any molecule may also show, per polarity: (label, shift from neutral, multiplier of M).
#  - formate and acetate adducts differ by exactly CH2, so allowing both links every homologue pair
#    (SM d34:0 <-> SM d33:0). Only the run's own modifier adduct is generated (detected from anchors).
#  - a generic water loss is formula-identical to the plasmalogen (PE 38:5 - H2O = PE P-38:6) and to
#    hydroxy-lipid relatives (PC[OH] 36:5 - H2O = PC 36:6); water loss only for classes that really lose it.
#  - Na - H = C2H4 - 3H2 within 2.4 mDa: sodium relations demand exact co-elution (ISO_RT), not rt_tol.
#  - [M+Na-2H]- is rare and coincidence-prone; dropped.
GENERIC = {
    "+": [("[M+H]+", 1.007276, 1), ("[M+Na]+", 22.989218, 1), ("[M+K]+", 38.963158, 1),
          ("[M+NH4]+", 18.033823, 1),
          ("[2M+H]+", 1.007276, 2), ("[2M+Na]+", 22.989218, 2), ("[2M+NH4]+", 18.033823, 2)],
    "-": [("[M-H]-", -1.007276, 1), ("[M+Cl]-", 34.969402, 1), ("[2M-H]-", -1.007276, 2)],
}
MODIFIER_ADDUCT = {"acetate": ("[M+OAc]-", 59.013851), "formate": ("[M+HCOO]-", 44.998201)}
TIGHT = ("[M+Na]+", "[M+K]+", "[2M+Na]+")          # coincidence-prone: require exact co-elution
CROSS_RT = 0.10      # min, after aligning the two polarity runs
GENERIC_LABELS = {g[0] for v in GENERIC.values() for g in v} | {"[M+OAc]-", "[M+HCOO]-"}
CHOLINE = ("PC", "SM", "LPC", "Plasmenyl-PC", "Plasmanyl-PC", "LysoPC")
CERAMIDE = ("Cer", "Cer[NS]", "Cer[NDS]", "Cer[AS]", "Cer[ADS]", "Cer[AP]", "Cer[NP]", "Cer[EOS]")
GLYCEROL = ("DG", "MG")
SPECIFIC = {
    "-": [(CHOLINE, "[M-CH3]-", -15.023475)],    # + [M+modifier-CH3]- added at run time
    "+": [(("PE", "Plasmenyl-PE", "Plasmanyl-PE", "LPE", "LysoPE"), "[M+H-141.019]+ (PE head loss)", 1.007276 - 141.019094),
          (("PS", "LPS", "LysoPS"), "[M+H-185.009]+ (PS head loss)", 1.007276 - 185.008923),
          (CHOLINE, "[M+H-183.066]+ (phosphocholine loss)", 1.007276 - 183.066044),
          (CERAMIDE, "[M+H-H2O]+", 1.007276 - 18.010565),
          (CERAMIDE, "[M+H-2H2O]+", 1.007276 - 2 * 18.010565),
          (GLYCEROL, "[M+H-H2O]+", 1.007276 - 18.010565)],
}

# Alkylammonium adducts [M+CnH2n+4N]+ of NH4-adduct classes: at every TG apex an ion at +n x 14.01565 Da from
# [M+NH4]+, formula-IDENTICAL to the NH4 adduct of the lipid with n more carbons. Proof it is not the homolog: TG
# standards alone (Lumos) -- the +28 ion's MS2 carries ONLY the parent's [M+H-RCOOH]+ ion and none 28 Da higher; a
# real homolog would also elute ~0.25 min per carbon later. Seen in Set 2 (Lumos 2026) and MTBLS394 QE, absent from
# Set 1 and LipidMS -> solvent/system dependent. Only n = 2; three guards: a run-level gate, the target's own class
# model (RT Model Z > ALKYL_MAX_Z keeps it), and the MS2 homolog veto.
ALKYL_CLASSES = ("TG", "DG", "CE", "MG")
ALKYL = [(f"[M+C{n}H{2 * n + 4}N]+ (alkylammonium, +{n}CH2)", 18.033823 + n * 14.01565, n) for n in (2,)]
ALKYL_GATE = 0.05      # share of named NH4 TG/CE/DG/MG anchors with a co-apex weaker +2CH2 row (TG alone 21%,
ALKYL_GATE_MIN = 3     # Set 2 10%, MTBLS394 QE 45% vs LipidMS 1%, Set 1 1%; off-mass controls 0%)
ALKYL_MAX_Z = -1.0
ALKYL_LABELS = {a[0] for a in ALKYL}
ALKYL_N = {a[0]: a[2] for a in ALKYL}

NA_LIKE = ("[M+Na]+", "[M+K]+", "[2M+Na]+")
MOD_LABELS = ("[M+OAc]-", "[M+HCOO]-")
# classes that form the modifier adduct in Neg. Anionic glycerophospholipids deprotonate instead -- generating
# their formate adduct made exact-formula coincidences (PS [M-H]- = PE(+1 DB) [M+HCOO]-).
FORMS_MOD_ADDUCT = ("PC", "SM", "LPC", "LysoPC", "Plasmenyl-PC", "Plasmanyl-PC", "PC[OH]",
                    "TG", "DG", "MG", "CE", "HexCer", "GlcCer", "GalCer", "SHexCer", "LacCer")
NO_H_CLASSES = ("TG", "DG", "CE", "MG")            # protonated ion barely forms: no [M+H]+ edge from these


def _family(cls):
    return {"LysoPC": "LPC", "LysoPE": "LPE", "LysoPS": "LPS"}.get(cls, cls)


def _binom_ratio(k, n_c):
    """Expected (M+k)/M from carbon alone: Poisson approx with lambda = 0.0107 n_C."""
    lam = 0.0107 * n_c
    return lam ** k / math.factorial(k)


def _pearson(a, b):
    pts = [(x, y) for x, y in zip(a, b) if x > 0 and y > 0]
    if len(pts) < 4:
        return None
    xs, ys = [math.log(p[0]) for p in pts], [math.log(p[1]) for p in pts]
    if statistics.pstdev(xs) == 0 or statistics.pstdev(ys) == 0:
        return None
    return round(statistics.correlation(xs, ys), 3)


def load_ms2(folder, max_files=12, log=None):
    """MS2 spectra (rt min, precursor, m/z array, intensity array) of up to `max_files` mzML, files whose
    name suggests QC / pool / DDA first."""
    if folder is None or not Path(folder).exists():
        return []
    import pyopenms as oms
    files = sorted(Path(folder).glob("*.mzML"),
                   key=lambda f: (not any(k in f.name.upper() for k in ("QC", "POOL", "DDA")), f.name))
    out = []
    for f in files[:max_files]:
        e = oms.MSExperiment()
        opts = oms.PeakFileOptions()
        opts.setMSLevels([2])
        fh = oms.MzMLFile()
        fh.setOptions(opts)
        fh.load(str(f), e)
        for sp in e:
            if sp.getMSLevel() != 2 or not sp.getPrecursors():
                continue
            m, i = sp.get_peaks()
            out.append((sp.getRT() / 60, sp.getPrecursors()[0].getMZ(), m, i))
    if log:
        log(f"  MS2 index {Path(folder).name}: {len(out)} spectra from {min(len(files), max_files)} files")
    return out


def _spectra(ms2, mz, rt, rt_win=0.1, ppm=10.0):
    return [s for s in ms2 if abs(s[0] - rt) <= rt_win and abs(s[1] - mz) / mz * 1e6 <= ppm]


def _sodiated_by_ms2(t, ms2, rt_win=0.1, ppm=10.0):
    """Does the target's OWN MS2 show it is a sodium adduct? True / False, or None with no MS2 or no rule.
    Positive evidence, not the absence of the protonated signature: a co-isolated protonated species can make
    184.07 the base peak of a sodiated precursor's spectrum.
      choline lipids:  146.9814 (C2H5O4PNa+) >= 2% of base, or [M+Na - 59.0735]+ (trimethylamine) >= 5%
      PE:              [M+Na - 43.0422]+ (aziridine) >= 5%"""
    fam = _family(t["Lipid Class"])
    mz, rt = float(t["Quant Ion"]), float(t["Retention Time (min)"])
    spectra = _spectra(ms2, mz, rt, rt_win, ppm)
    if not spectra:
        return None
    if fam in ("PC", "SM", "LPC", "Plasmenyl-PC", "Plasmanyl-PC"):
        tests = [(146.9814, 0.02), (mz - 59.0735, 0.05)]
    elif fam in ("PE", "LPE", "Plasmenyl-PE", "Plasmanyl-PE"):
        tests = [(mz - 43.0422, 0.05)]
    else:
        return None
    votes = []
    for _, _, m, i in spectra:
        if len(i) == 0:
            continue
        base = float(i.max())
        votes.append(any(max((float(ii) for mm, ii in zip(m, i) if abs(mm - x) <= 0.01), default=0.0) >= f * base
                         for x, f in tests))
    if not votes:
        return None
    return sum(votes) > len(votes) / 2


def _own_class_by_ms2(t, pol, ms2, rt_win=0.1, ppm=10.0):
    """True if the target's own Neg MS2 shows its CLAIMED class's signature, so it is not someone else's ion.
      PS: [M-H - 87.0320]- (serine loss) >= 20% of base;  PI: 241.0119 or 259.0224 >= 5% of base.
    Other classes: None. PE-NMe2 vs PC [M-CH3]- is the same anion -- no signature exists."""
    if pol != "Neg" or not t["Identification"]:
        return None
    fam = _family(t["Lipid Class"])
    mz, rt = float(t["Quant Ion"]), float(t["Retention Time (min)"])
    if fam in ("PS", "LPS"):
        tests = [(mz - 87.0320, 0.20)]
    elif fam in ("PI", "LPI", "LysoPI"):
        tests = [(241.0119, 0.05), (259.0224, 0.05)]
    else:
        return None
    votes = []
    for _, _, m, i in _spectra(ms2, mz, rt, rt_win, ppm):
        if len(i) == 0:
            continue
        base = float(i.max())
        votes.append(any(max((float(ii) for mm, ii in zip(m, i) if abs(mm - x) <= 0.01), default=0.0) >= f * base
                         for x, f in tests))
    if not votes:
        return None
    return sum(votes) > len(votes) / 2


def _homolog_by_ms2(t, a, n, ms2, rt_win=0.1, ppm=10.0, tol=0.01, min_shift=0.3):
    """True if the target's own MS2 carries the anchor's fragments shifted by n CH2: a real homolog contributes."""
    def spectra(r):
        mz, rt = float(r["Quant Ion"]), float(r["Retention Time (min)"])
        return [x for x in _spectra(ms2, mz, rt, rt_win, ppm) if len(x[3])]
    ts, as_ = spectra(t), spectra(a)
    if not ts or not as_:
        return None
    _, _, pm, pi = max(as_, key=lambda x: float(x[3].sum()))
    pm = [float(m) for m, i in zip(pm, pi) if i >= 0.02 * pi.max()]
    on = lambda x: any(abs(x - y) <= tol for y in pm)
    votes = []
    for _, _, m, i in ts:
        keep = [(float(mm), float(ii)) for mm, ii in zip(m, i) if ii >= 0.05 * i.max()]
        tot = sum(ii for _, ii in keep)
        shifted = sum(ii for mm, ii in keep if on(mm - n * 14.01565) and not on(mm))
        votes.append(shifted / tot >= min_shift if tot else False)
    return sum(votes) > len(votes) / 2


def _choline_anchors(pos_rows, neg_rows, log=None):
    """Choline lipids named by MS2 in Pos, with neutral mass, and the Neg-minus-Pos RT shift measured on the same
    names seen in both polarities."""
    pos = [r for r in pos_rows if r["Identification"]]
    neg = [r for r in neg_rows if r["Identification"]]
    anchors = []
    for r in pos:
        if _family(r["Lipid Class"]) in ("PC", "SM", "LPC", "Plasmenyl-PC", "Plasmanyl-PC") \
                and r["Identification Source"] == "MS2" and r["Adduct"].split(";")[0].strip() == "[M+H]+":
            anchors.append((r["Identification"], float(r["Theoretical m/z"] or r["Quant Ion"]) - 1.007276,
                            float(r["Retention Time (min)"])))
    key = lambda r: (r["Lipid Class"], r["Shorthand (LSI)"] or r["Identification"])
    pos_rt = {key(r): float(r["Retention Time (min)"]) for r in pos if r["Identification Source"] == "MS2"}
    d = [float(r["Retention Time (min)"]) - pos_rt[key(r)] for r in neg
         if r["Identification Source"] == "MS2" and key(r) in pos_rt]
    d = [x for x in d if abs(x) < 0.5]
    shift = statistics.median(d) if len(d) >= 5 else 0.0
    if log:
        log(f"  cross-polarity: {len(anchors)} Pos choline anchors; Neg-Pos RT shift {shift:+.3f} min "
            f"from {len(d)} shared names")
    return anchors, shift


def ion_family(rows, pol, ms2, rt_tol=0.05, cross=None, log=None) -> dict:
    """Label rows of one polarity's Final_Results that are another row's ion. Returns
    {Compound Group: (Ion Of, Ion Relation, Ion ppm, Ion dRT, Ion r)}."""
    say = log or (lambda _: None)
    rows = [dict(r) for r in rows]
    if not rows:
        return {}
    header = list(rows[0].keys())
    samples = [s for s in header[header.index("Mass Error (ppm)") + 1:] if s]
    sign = "+" if pol == "Pos" else "-"
    errs = sorted(float(r["Mass Error (ppm)"]) for r in rows if r["Identification"] and r["Mass Error (ppm)"])
    offset = errs[len(errs) // 2] if errs else 0.0
    for r in rows:
        r["_mz"] = float(r["Quant Ion"]) / (1 + offset * 1e-6)
        r["_rt"] = float(r["Retention Time (min)"])
        r["_area"] = float(r["Area (max)"] or 0)
        r["_areas"] = [float(r[s] or 0) for s in samples]
    rows_by_rt = sorted(range(len(rows)), key=lambda i: rows[i]["_rt"])
    generic = list(GENERIC[sign])
    specific = list(SPECIFIC[sign])
    minor = set()
    if sign == "-":
        n_ac = sum("[M+Ac-H]-" in r["Adduct"] or "CH3COO" in r["Adduct"] for r in rows if r["Identification"])
        n_fa = sum("[M+FA-H]-" in r["Adduct"] or "HCOO" in r["Adduct"] for r in rows if r["Identification"])
        mod = "acetate" if n_ac >= n_fa else "formate"
        lab, sh = MODIFIER_ADDUCT[mod]
        generic.append((lab, sh, 1))
        specific.append((CHOLINE, lab[:-2] + "-CH3]-", sh - 15.023475))
        # the OTHER modifier's adduct still forms from traces (LipidMS acetate run: formate clusters ~20% of
        # acetate). Exactly CH2 from the main adduct of the homologue: across classes, exact co-elution only.
        minor_lab, minor_sh = MODIFIER_ADDUCT["formate" if mod == "acetate" else "acetate"]
        generic.append((minor_lab, minor_sh, 1))
        minor.add(minor_lab)
        say(f"  {pol}: modifier from anchors = {mod} (acetate {n_ac}, formate {n_fa})")

    alkyl_on = False
    if sign == "+":
        top = [r for r in rows if r["Identification"] and r["Lipid Class"] in ALKYL_CLASSES
               and r["Adduct"].split(";")[0].strip() == "[M+NH4]+" and "DECOY" not in r["Identification"].upper()]

        def has_plus28(a):
            mz = a["_mz"] + 2 * 14.01565
            return any(abs(t["_rt"] - a["_rt"]) <= ISO_RT and abs(t["_mz"] - mz) / mz * 1e6 <= PPM
                       and t["_area"] < a["_area"] for t in rows)
        k = sum(has_plus28(a) for a in top)
        alkyl_on = bool(top) and k >= ALKYL_GATE_MIN and k / len(top) >= ALKYL_GATE
        say(f"  {pol}: alkylammonium gate {k}/{len(top)} NH4 anchors with a co-apex +2CH2 row -> "
            f"{'ON' if alkyl_on else 'off'} (needs >= {ALKYL_GATE_MIN} and >= {ALKYL_GATE:.0%})")

    edges = {}            # target index -> list of (anchor index, relation, ppm, drt)
    for ai, a in enumerate(rows):
        if "DECOY" in a["Identification"].upper():
            continue
        named = bool(a["Identification"])
        # an unnamed row anchors only generic relations, assuming the plain (de)protonated ion
        add = a["Adduct"].split(";")[0].strip() if named else ("[M+H]+" if sign == "+" else "[M-H]-")
        shift = ADDUCT.get(add)
        if shift is None:
            continue
        neutral = a["_mz"] - shift
        cls = a["Lipid Class"]
        forms_mod = cls in FORMS_MOD_ADDUCT or cls.startswith(("Cer", "HexCer", "GlcCer", "GalCer"))
        ions = [(lab, neutral * mult + sh) for lab, sh, mult in generic
                if (named or lab in TIGHT or lab in ("[M+NH4]+", "[M+Cl]-"))
                and not (lab == "[M+H]+" and cls in NO_H_CLASSES)
                and not (lab in MOD_LABELS and not forms_mod)]
        if named:
            ions += [(lab, neutral + sh) for classes, lab, sh in specific if cls in classes]
        if alkyl_on and named and sign == "+" and add == "[M+NH4]+" and cls in ALKYL_CLASSES:
            ions += [(lab, neutral + sh) for lab, sh, _ in ALKYL]
        n_c = 0.065 * neutral
        for ti in rows_by_rt:
            t = rows[ti]
            if ti == ai:
                continue
            weaker = t["_area"] < a["_area"]           # generic relations and isotopes need it; class-specific
                                                       # in-source products have one direction only
            drt = t["_rt"] - a["_rt"]
            if abs(drt) > max(rt_tol, ISO_RT_WIDE):
                continue
            hit = None
            for lab, mz in ions:
                if not weaker and lab in GENERIC_LABELS:
                    continue
                e = (t["_mz"] - mz) / mz * 1e6
                own = add.replace("[M+FA-H]-", "[M+HCOO]-").replace("[M+Ac-H]-", "[M+OAc]-")
                if lab in minor and (abs(drt) > ISO_RT or _family(t["Lipid Class"]) == _family(cls)):
                    continue
                if lab in ALKYL_LABELS and (abs(drt) > ISO_RT or not weaker):
                    continue
                if lab in ALKYL_LABELS and t["Identification"] and t.get("RT Model Z") not in ("", None) \
                        and float(t["RT Model Z"]) > ALKYL_MAX_Z:
                    continue
                if abs(drt) > rt_tol:
                    break
                if abs(e) <= PPM and lab != own and (lab not in TIGHT or abs(drt) <= ISO_RT):
                    if lab in NA_LIKE and t["Identification"]:
                        same_lipid = (_family(t["Lipid Class"]) == _family(cls)
                                      and (t["Shorthand (LSI)"] or t["Identification"]) == (a["Shorthand (LSI)"] or a["Identification"]))
                        if not same_lipid:
                            # mass alone cannot separate Na-H from C2H4-3H2 (3 ppm): the target's own MS2 must
                            # SHOW sodium; no MS2 or a class without a rule -> leave it alone
                            if _sodiated_by_ms2(t, ms2) is not True:
                                continue
                    if lab in ALKYL_LABELS and _homolog_by_ms2(t, a, ALKYL_N[lab], ms2) is True:
                        continue
                    hit = (lab, e)
                    break
            if hit is None and abs(drt) <= ISO_RT_WIDE and weaker:
                for k in (1, 2, 3):
                    target = a["_mz"] + k * ISO
                    e = (t["_mz"] - target) / target * 1e6
                    own_err = None
                    if t["Identification"] and t.get("Theoretical m/z"):
                        th = float(t["Theoretical m/z"])
                        own_err = (t["_mz"] - th) / th * 1e6
                    mass_decides = (abs(e) <= 3.0 and own_err is not None and abs(own_err) >= abs(e) + 4.0)
                    if (abs(drt) > ISO_RT or not named) and not mass_decides:
                        continue                        # an unnamed parent needs the mass test (held-out ST004797)
                    if abs(e) <= ISO_PPM:
                        ratio = t["_area"] / a["_area"] if a["_area"] else 9
                        if ratio <= ISO_SLACK * _binom_ratio(k, n_c):
                            hit = (f"M+{k} isotope (area ratio {ratio:.3f} vs exp {_binom_ratio(k, n_c):.3f}"
                                   + (f"; own name {own_err:+.1f} ppm" if own_err is not None else "") + ")", e)
                        break
            if hit and _own_class_by_ms2(t, pol, ms2) is True:
                hit = None                              # its own MS2 says it is what it is called
            if hit:
                edges.setdefault(ti, []).append((ai, hit[0], round(hit[1], 1), round(drt, 3)))

    # cross-polarity anchors: a choline lipid named in POSITIVE mode at this RT fixes a neutral mass whose
    # negative-mode ions ([M+mod]-, [M-CH3]-, [M+mod-CH3]-, and their M+1..M+3) cannot be other lipids.
    cross_hits = {}
    if cross and sign == "-":
        anchors, shift_rt = cross
        lab_mod = [g for g in generic if g[0] in ("[M+OAc]-", "[M+HCOO]-") and g[0] not in minor][0]
        rel = [(lab_mod[0], lab_mod[1]), ("[M-CH3]-", -15.023475), (lab_mod[0][:-2] + "-CH3]-", lab_mod[1] - 15.023475)]
        for name, neutral, prt in anchors:
            for ti, t in enumerate(rows):
                if abs(t["_rt"] - (prt + shift_rt)) > CROSS_RT or _family(t["Lipid Class"]) in ("PC", "SM", "LPC", "Plasmenyl-PC", "Plasmanyl-PC"):
                    continue
                # a name its own class retention model supports is not overridden by a mass coincidence.
                # Not for PE-NMe/PE-NMe2: rows named so are mostly these same PC ions, so their class
                # model is fitted on them and would vouch for them.
                z = t.get("RT Model Z", "")
                if (t["Identification"] and z not in ("", None) and abs(float(z)) <= 2.0
                        and not t["Lipid Class"].startswith("PE-NMe")):
                    continue
                if _own_class_by_ms2(t, "Neg", ms2) is True:
                    continue
                for lab, sh in rel:
                    # isotopes only of the MAIN Neg ion (the modifier adduct)
                    for k in ((0, 1, 2, 3) if lab == lab_mod[0] else (0,)):
                        mz = neutral + sh + k * ISO
                        e = (t["_mz"] - mz) / mz * 1e6
                        if abs(e) <= PPM and (t["Identification"] or k == 0):
                            cross_hits[ti] = (name, f"{lab}{' M+' + str(k) if k else ''} of the Pos-named lipid", round(e, 1),
                                              round(t["_rt"] - prt - shift_rt, 3))
                            break
                    if ti in cross_hits:
                        break

    # resolve: a row explained by an anchor that is itself unexplained is an ion; mutual -> larger area wins
    explained = {}
    for ti, lst in edges.items():
        for ai, rel, e, drt in sorted(lst, key=lambda x: -rows[x[0]]["_area"]):
            mutual = any(a2 == ti for a2, *_ in edges.get(ai, []))
            if mutual and rows[ai]["_area"] < rows[ti]["_area"]:
                continue
            explained[ti] = (ai, rel, e, drt)
            break
    # an anchor that is itself explained cannot explain others (chains collapse to the root)
    for ti in list(explained):
        ai = explained[ti][0]
        seen = {ti}
        while ai in explained and ai not in seen:
            seen.add(ai)
            ai = explained[ai][0]
        explained[ti] = (ai,) + explained[ti][1:]

    labels = {}
    for i, r in enumerate(rows):
        if i in cross_hits and i not in explained:
            name, rel, e, drt = cross_hits[i]
            labels[r["Compound Group"]] = (name + " (Pos)", rel, e, drt, "")
        elif i in explained:
            ai, rel, e, drt = explained[i]
            a = rows[ai]
            labels[r["Compound Group"]] = (f"{a['Identification'] or 'group ' + a['Compound Group']} {a['Adduct']}",
                                           rel, e, drt, _pearson(r["_areas"], a["_areas"]))
    say(f"  {pol}: ion family offset {offset:+.2f} ppm; {len(labels)} rows labelled as ions "
        f"({sum(1 for r in rows if r['Identification'] and r['Compound Group'] in labels)} of them named)")
    return labels


# ════════════════════════════════════════════════════════════════════════════════════════════════
# shape test (frozen parent check v1, apex rule v1b)
# ════════════════════════════════════════════════════════════════════════════════════════════════

# parent m/z = target m/z + shift
SHIFT = {"Pos": {"Na swap ([M+Na]+ of [M+H]+)": -21.98194, "K swap ([M+K]+ of [M+H]+)": -37.95588,
                 "NH4 swap ([M+NH4]+ of [M+H]+)": -17.02655, "Na swap ([M+Na]+ of [M+NH4]+)": -4.95539,
                 "K swap ([M+K]+ of [M+NH4]+)": -20.92934, "H2O loss": +18.01056, "M+2 of DB+1": -2.01565},
         "Neg": {"formate swap ([M+HCOO]- of [M-H]-)": -46.00548, "Cl swap ([M+Cl]- of [M-H]-)": -35.97613,
                 "demethylation ([M-CH3]- of [M+HCOO]-)": +60.02113, "H2O loss": +18.01056,
                 "M+2 of DB+1": -2.01565}}
H2O_PARENT = re.compile(r"OH|LysoP|LPC|LPE|Cer|CE |FA ")
# the 'parent' is the same species in another adduct when the target's own adduct list has the swapped adduct
SAME = {"Na swap": "Na]", "K swap": "K]", "NH4 swap": "NH4]", "formate swap": "FA-H]", "Cl swap": "Cl]",
        "H2O loss": "H2O]"}
TOL, TOL_M2, DRT_PAIR, DRT_XIC, DRT_MS2 = 5e-6, 12e-6, 0.05, 0.25, 0.15


def _sumcomp(n):
    """TG 16:0_18:1_18:1 -> TG 52:2."""
    n = n.split(" [")[0]
    m = re.match(r"^(\S+) ((?:OH-|[dtmOP]-?)?)([\d:_/;]+)(.*)$", n)
    if not m or "_" not in m.group(3) and "/" not in m.group(3):
        return n
    ch = re.findall(r"(\d+):(\d+)", m.group(3))
    return f"{m.group(1)} {m.group(2)}{sum(int(c) for c, _ in ch)}:{sum(int(d) for _, d in ch)}"


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return 0.0


def _m2_expected(mz):
    """Expected M+2/M for a lipid of this m/z (nC ~ m/z / 17, ~8 O)."""
    nc = mz / 17
    return (nc * 0.0107) ** 2 / 2 + 8 * 0.002


def _scan_file(args):
    """One pass over one mzML: XIC points for each requested (mz, rt window)."""
    from pyteomics import mzml
    path, xic_req = args
    xic = {k: [] for k in xic_req}
    lo_t = min([w[1] for w in xic_req.values()] + [1e9])
    hi_t = max([w[2] for w in xic_req.values()] + [-1])
    with mzml.MzML(path) as rd:
        for s in rd:
            if s.get("ms level") == 2:
                continue
            st = s["scanList"]["scan"][0]["scan start time"]
            t = float(st) / (60 if getattr(st, "unit_info", "minute") == "second" else 1)
            if t < lo_t or t > hi_t:
                continue
            mzs, it = s["m/z array"], s["intensity array"]
            for k, (m, a, b) in xic_req.items():
                if a <= t <= b:
                    lo, hi = np.searchsorted(mzs, [m * (1 - 5e-6), m * (1 + 5e-6)])
                    xic[k].append((t, float(it[lo:hi].max()) if hi > lo else 0.0))
    return path, xic


def _explained(shape_r, apex_offset, ratio, ratio_limit, params: CoelutionFilter) -> bool:
    """The DERIVATIVE test: same shape over the target's peak, same apex, derivative-sized ratio. v1b: the apex may
    be `apex_tol_high_r` off when the shapes correlate at `high_r` or better (NaN shape never passes)."""
    apex_tol = params.apex_tol_high_r if shape_r >= params.high_r else params.apex_tol
    return shape_r >= 0.9 and abs(apex_offset) <= apex_tol and ratio <= ratio_limit


def shape_test(unf, fin, mzml_dir, params: CoelutionFilter, pol, log=None) -> dict:
    """Delivered MS2 rows that are a derivative ion of a stronger co-eluting partner. Returns
    {Compound Group: (partner name, relation, shape r, apex offset, ratio)} for DERIVATIVE rows only.

    The frozen v1 also reported MS2 modified cosine, m/z 184 and the PC/PE isobar flag; none of them decided
    anything and they are not ported."""
    say = log or (lambda _: None)
    fin = [r for r in fin if r["Identification"] and r["Identification Source"] == "MS2"]
    if not fin or not unf:
        return {}
    cols = list(unf[0].keys())
    samples = [c for c in cols[cols.index("Mass Error (ppm)") + 1:] if c and not c.startswith("Blank")]
    samples = [s for s in samples if mzml_dir and (Path(mzml_dir) / f"{s}.mzML").exists()]
    files = [str(Path(mzml_dir) / f"{s}.mzML") for s in samples]
    cand = []
    for r in fin:
        mz, t, a = _num(r["Quant Ion"]), _num(r["Retention Time (min)"]), _num(r["Area (max)"])
        hits = []
        for lab, d in SHIFT[pol].items():
            p, tol = mz + d, TOL_M2 if lab.startswith("M+2") else TOL
            for u in unf:
                if u["Compound Group"] == r["Compound Group"]:
                    continue
                m2, t2, a2 = _num(u["Quant Ion"]), _num(u["Retention Time (min)"]), _num(u["Area (max)"])
                if abs(m2 - p) / p < tol and abs(t2 - t) <= DRT_PAIR and a2 >= a:
                    hits.append((lab, u))
        cand.append(dict(r=r, mz=mz, rt=t, area=a, hits=hits))
    xic_req = {}
    for i, c in enumerate(cand):
        if c["hits"]:
            xic_req[("t", i)] = (c["mz"], c["rt"] - DRT_XIC, c["rt"] + DRT_XIC)
        for j, (lab, u) in enumerate(c["hits"]):
            xic_req[("p", i, j)] = (_num(u["Quant Ion"]), c["rt"] - DRT_XIC, c["rt"] + DRT_XIC)
    if not xic_req or not files:
        say(f"  {pol}: shape test, {len(cand)} delivered, 0 with a partner")
        return {}
    with ProcessPoolExecutor(min(8, len(files))) as ex:
        xic = {f: x for f, x in ex.map(_scan_file, [(f, xic_req) for f in files])}

    out = {}
    for i, c in enumerate(cand):
        r = c["r"]
        name = r["Identification"]
        best = None
        for j, (lab, u) in enumerate(c["hits"]):
            pm = _num(u["Quant Ion"])
            rs, offs = [], []
            for f in files:
                xt, xp = xic[f][("t", i)], xic[f][("p", i, j)]
                if len(xt) < 5:
                    continue
                yt, yp = np.array([v for _, v in xt]), np.array([v for _, v in xp])
                tt = np.array([t for t, _ in xt])
                near = np.abs(tt - c["rt"]) <= DRT_PAIR  # the target's own apex, then its peak +-0.1 min
                if not near.any() or yt[near].max() <= 0:
                    continue
                apex = tt[near][yt[near].argmax()]
                keep = (np.abs(tt - apex) <= 0.1) & ((yt > 0) | (yp > 0))
                if keep.sum() >= 5 and yt[keep].std() > 0 and yp[keep].std() > 0:
                    rs.append(np.corrcoef(yt[keep], yp[keep])[0, 1])
                    offs.append(apex - tt[keep][yp[keep].argmax()])
            ratios = [_num(r[s]) / _num(u[s]) for s in samples if _num(r[s]) > 0 and _num(u[s]) > 0]
            ratio = float(np.median(ratios)) if ratios else c["area"] / _num(u["Area (max)"])
            rmed = float(np.median(rs)) if rs else np.nan
            off = float(np.median(offs)) if offs else np.nan
            lim = 1.5 * _m2_expected(pm) if lab.startswith("M+2") else 0.5
            explained = _explained(rmed, off, ratio, lim, params)
            # a water loss needs a parent that loses water -- unnamed, or hydroxylated / lyso / ceramide /
            # cholesteryl ester. A named diacyl parent minus H2O has the formula of the plasmalogen.
            if lab == "H2O loss" and u["Identification"] and not H2O_PARENT.search(u["Identification"]):
                explained = False
            same = (u["Identification"] and _sumcomp(u["Identification"]) == _sumcomp(name)) or any(
                lab.startswith(k) and a in r["Adduct"] for k, a in SAME.items())
            v = "SAME-SPECIES" if same else "DERIVATIVE" if explained else "PARTNER"
            h = dict(verdict=v, shift=lab, parent=u["Identification"] or f"group {u['Compound Group']} (no ID)",
                     shape_r=rmed, apex_offset=off, ratio=ratio)
            order = {"DERIVATIVE": 0, "SAME-SPECIES": 1, "PARTNER": 2}
            if best is None or (order[v], -np.nan_to_num(rmed)) < (order[best["verdict"]],
                                                                    -np.nan_to_num(best["shape_r"])):
                best = h
        if best and best["verdict"] == "DERIVATIVE":
            out[r["Compound Group"]] = (best["parent"], best["shift"], best["shape_r"], best["apex_offset"],
                                        best["ratio"])
    say(f"  {pol}: shape test, {len(cand)} delivered, {sum(bool(c['hits']) for c in cand)} with a partner, "
        f"{len(out)} derivative")
    return out


# ════════════════════════════════════════════════════════════════════════════════════════════════
# applying both
# ════════════════════════════════════════════════════════════════════════════════════════════════

def _read(path: Path):
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        return [f for f in (reader.fieldnames or []) if f], list(reader)


def _write(path: Path, fields, rows):
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _with_label_columns(fields):
    """The label columns go INSIDE the metadata block, never after the injections.

    ⚠ Readers take every column after the last metadata column to be an injection
    (`qc.load_run`, see `peakfinder.META_COLUMNS`). Appended at the end, `Ion Of` was read as a
    sample and the quality report died on 'group 166 (no ID)' (NIST + TG spike run, 2026-11-04).
    Placed after `Identification Source`, they sit between metadata columns and nothing counts them.
    """
    if all(c in fields for c in LABEL_COLUMNS):
        return list(fields)
    out = [c for c in fields if c not in LABEL_COLUMNS]
    at = (out.index("Identification Source") + 1 if "Identification Source" in out
          else out.index("Mass Error (ppm)") if "Mass Error (ppm)" in out else len(out))
    return out[:at] + LABEL_COLUMNS + out[at:]


def apply(analysis, mzml_dirs: dict, params: CoelutionFilter | None = None, log=None) -> dict:
    """Run both tests on `analysis/Results/<pol>` and remove their union from the delivered tables.

    `mzml_dirs`: {polarity: folder of that polarity's mzML}. Returns the summary also written to
    `analysis/coelution_filter.json`."""
    from .polarity import _normalise, _numeric_columns
    params = params or CoelutionFilter()
    say = log or (lambda _: None)
    root = Path(analysis) / "Results"
    pols = [p for p in POLARITIES if (root / p / FINAL).exists() and (root / p / FILTERED).exists()]
    if not params.enabled or not pols:
        return {}
    final = {p: _read(root / p / FINAL) for p in pols}
    cross = (_choline_anchors(final["Pos"][1], final["Neg"][1], log=say)
             if {"Pos", "Neg"} <= set(pols) else None)
    summary = {"settings": asdict(params), "removed": {}, "labelled_in_final": {}}
    for pol in pols:
        ms2 = load_ms2(mzml_dirs.get(pol), params.max_ms2_files, log=say)
        fields_f, rows_f = final[pol]
        labels = ion_family(rows_f, pol, ms2, params.rt_tol, cross if pol == "Neg" else None, log=say)
        _, unf = _read(root / pol / UNFILTERED)
        fields_d, delivered = _read(root / pol / FILTERED)
        shape = shape_test(unf, delivered, mzml_dirs.get(pol), params, pol, log=say)

        # Final_Results.csv: labelled, never deleted
        for r in rows_f:
            g = r["Compound Group"]
            a, b = labels.get(g), shape.get(g)
            if a:
                r.update(zip(LABEL_COLUMNS, list(a) + ["ion family" + (" + shape" if b else "")]))
            elif b:
                parent, rel, rr, off, ratio = b
                r.update({"Ion Of": parent, "Ion Relation": rel, "Ion ppm": "",
                          "Ion dRT": f"{off:+.3f}", "Ion r": f"{rr:.3f}", "Ion Test": "shape"})
        _write(root / pol / FINAL, _with_label_columns(fields_f), rows_f)

        drop = set(labels) | set(shape)
        why = {}
        for g in drop:
            if g in labels:
                why[g] = f"co-elution: {labels[g][1]} of {labels[g][0]}"
            else:
                why[g] = f"co-elution: {shape[g][1]} of {shape[g][0]} (shape r {shape[g][2]:.2f})"
        removed = 0
        for name in (FILTERED, NORMALISED):
            path = root / pol / name
            if not path.exists():
                continue
            fields, rows = _read(path)
            kept = [r for r in rows if r.get("Compound Group") not in drop]
            if name == FILTERED:
                removed = len(rows) - len(kept)
            elif kept and len(kept) != len(rows):
                _normalise(kept, _numeric_columns(fields))
            _write(path, fields, kept)
        fields_u, unf = _read(root / pol / UNFILTERED)
        if "Filter Status" in fields_u:
            for r in unf:
                g = r.get("Compound Group")
                if g in why and not (r.get("Filter Status") or "").strip():
                    r["Filter Status"] = why[g]
            _write(root / pol / UNFILTERED, fields_u, unf)
        summary["removed"][pol] = {"rows": removed,
                                   "ion family": sum(1 for g in labels if g not in shape and any(
                                       r["Compound Group"] == g for r in delivered)),
                                   "shape": sum(1 for g in shape if g not in labels),
                                   "both": sum(1 for g in shape if g in labels)}
        summary["labelled_in_final"][pol] = len(drop)
        say(f"  {pol}: co-elution filter removed {removed} delivered rows "
            f"(ion family {summary['removed'][pol]['ion family']}, shape {summary['removed'][pol]['shape']}, "
            f"both {summary['removed'][pol]['both']})")
    (Path(analysis) / "coelution_filter.json").write_text(json.dumps(summary, indent=2))
    return summary
