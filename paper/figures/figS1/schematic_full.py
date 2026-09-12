"""Figure S1: the full pipeline schematic, every stage and every decision, in the same conventions
as Figure 1 (drawing code in `manuscript/figures/schematic_lib.py`): a question in every diamond
with a labelled yes AND no exit, sinks for rows that leave the delivered table, repeated input
boxes instead of arrows that span bands, and the geometric check.

    .venv/bin/python manuscript/figures/figS1/schematic_full.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from schematic_lib import CELLS, band, box, dec, edge, render, reset, row, sink, write  # noqa: E402

reset(1240)

# ── 1 conversion ─────────────────────────────────────────────────────────────────────────────
band(20, 330, "1 — CONVERSION   (reader chosen from the file itself, never from its suffix; one mzML cache per study)")
raw, cached, which = row(56, 58, [(220, "vendor data\n.raw  .d  .wiff, Pos/ and Neg/"),
                                  (220, "cached mzML newer\nthan the raw file?"),
                                  (190, "which vendor\nwrote the file?")], "proc", gap=70, x0=80)
CELLS[raw]["kind"] = "input"; dec(cached, 51); dec(which, 51)
trfp, mscv, tims = row(150, 56, [(220, "ThermoRawFileParser\nThermo .raw"),
                                 (230, "msconvert (wine)\nAgilent, Waters, Sciex"),
                                 (220, "alphatims + writer\nBruker timsTOF")], "proc", gap=40, x0=80)
prof = box(CELLS[mscv]["x"] + 30, 244, 170, 68, "profile\nspectra?", "dec")
cent = box(CELLS[prof]["x"] + 230, 250, 150, 56, "centroid once", "proc")
mzml = box(980, 250, 220, 56, ".mzML\ncached per study, atomic write", "file")
edge(raw, cached); edge(cached, which, "no"); edge(cached, mzml, "yes")
edge(which, trfp, "Thermo"); edge(which, mscv, "Agilent, Waters, Sciex"); edge(which, tims, "Bruker")
edge(trfp, prof); edge(mscv, prof); edge(tims, prof)
edge(prof, cent, "yes"); edge(cent, mzml); edge(prof, mzml, "no", via="under")

# ── 2 identification ─────────────────────────────────────────────────────────────────────────
band(370, 330, "2 — IDENTIFICATION   (per injection: every MS2 spectrum against the target AND the decoy libraries)")
ms2, art, srch, offs, pur, chain = row(408, 58, [(120, "MS2 scans\nfrom the mzML"), (190, "artefact present in\n≥ 95% of spectra?"),
                                                 (190, "score every match\nLipiDex algorithm, exact"),
                                                 (150, "apply the\nmass offset"), (170, "fatty-acid\npurity ≥ 75%?"),
                                                 (160, "name at\nchain level")], "proc", gap=30)
dec(art, 403); dec(pur, 403)
strip = box(CELLS[art]["x"] + 5, 520, 180, 50, "strip the artefact peak", "proc")
tgt = box(CELLS[srch]["x"] - 40, 600, 150, 56, "target libraries\nHCD, per modifier", "proc")
dc1 = box(CELLS[srch]["x"] + 125, 600, 160, 56, "decoy: chains shifted\nby whole CH2 units", "proc")
dc2 = box(CELLS[srch]["x"] + 300, 600, 160, 56, "decoy: fragments displaced\nby an impossible mass", "proc")
sumc = box(CELLS[pur]["x"], 520, 170, 50, "name at\nsum composition", "proc")
sres = box(CELLS[chain]["x"], 520, 160, 50, "search/*.csv\nper injection", "file")
edge(mzml, ms2); edge(ms2, art); edge(art, strip, "yes"); edge(strip, srch); edge(art, srch, "no")
edge(tgt, srch); edge(dc1, srch); edge(dc2, srch)
edge(srch, offs); edge(offs, pur); edge(pur, chain, "yes"); edge(pur, sumc, "no"); edge(chain, sres); edge(sumc, sres)

# ── 3 feature detection ──────────────────────────────────────────────────────────────────────
band(720, 250, "3 — FEATURE DETECTION   (+ per-file mass calibration as a second pass)")
mz3, feat, align, link, drift = row(758, 58, [(110, ".mzML"), (170, "detect features\n(pyOpenMS)"), (170, "align retention\n(pose clustering)"),
                                              (170, "link injections\n10 ppm, 30 s"), (190, "per-file drift\n≥ 2 ppm?")], "proc", gap=40)
CELLS[mz3]["kind"] = "input"; dec(drift, 753)
noise = box(CELLS[link]["x"] - 20, 870, 200, 56, "measure the noise floor\n(reported, not applied)", "proc")
recal = box(CELLS[drift]["x"] + 20, 870, 240, 56, "write calibrated mzML,\nprecursors shifted too", "proc")
edge(mz3, feat); edge(feat, align); edge(align, link); edge(link, drift)
edge(drift, recal, "yes"); edge(recal, feat, "second pass", True); edge(drift, noise, "no")

# ── 4 peak finder ────────────────────────────────────────────────────────────────────────────
band(990, 420, "4 — PEAK FINDER   (join identifications to features, then the reproduced filters)")
s4, join, thr, name = row(1028, 58, [(110, "search/*.csv"), (150, "join MS2\n→ features"),
                                     (180, "dot ≥ 500 and\nreverse dot ≥ 700?"), (170, "name the row\nmolecular or sum")], "proc", gap=40, x0=80)
CELLS[s4]["kind"] = "input"; dec(thr, 1023)
k4a = sink(thr, 1135)
rtm, offm = row(1028, 58, [(170, "fit the per-class\nretention model"), (170, "off its class\nmodel?")], "proc", gap=40, x0=CELLS[name]["x"] + 210)
dec(offm, 1023)
k4b = sink(offm, 1135)
add, ins, red, win = row(1225, 58, [(190, "adduct, dimer or isotope\nof another peak?"), (170, "in-source\nfragment?"),
                                    (170, "redundant\nidentification?"), (210, "outside the class RT\nwindow? (off by default)")], "proc", gap=40, x0=80)
for d in (add, ins, red, win):
    dec(d, 1220)
k4c = sink(add, 1330); k4d = sink(ins, 1330); k4e = sink(red, 1330); k4f = sink(win, 1330)
edge(noise, join); edge(s4, join); edge(join, thr); edge(thr, name, "yes"); edge(thr, k4a, "no")
edge(name, rtm); edge(rtm, offm); edge(offm, k4b, "yes"); edge(offm, add, "no")
edge(add, ins, "no"); edge(add, k4c, "yes"); edge(ins, red, "no"); edge(ins, k4d, "yes")
edge(red, win, "no"); edge(red, k4e, "yes"); edge(win, k4f, "yes")

# ── 5 post-filters ───────────────────────────────────────────────────────────────────────────
band(1430, 410, "5 — POST-FILTERS")
split, adpr, blk, ffa = row(1468, 58, [(180, "split-peak merge\n(needs pooled QCs)"), (180, "adduct-pair removal\n(tracking ρ ≥ 0.8)"),
                                       (180, "sample ≥ 3 × blank?"), (190, "free fatty acids\nmass + retention surface")], "proc", gap=45)
dec(blk, 1463)
k5a = sink(blk, 1575)
ext, pres, gap = row(1660, 58, [(200, "retention-model extension\nMSI level 3, labelled"), (200, "detected in ≥ 60%\nof one group?"),
                                (200, "gap filling\nre-integrated from the mzML")], "proc", gap=60)
dec(pres, 1655)
k5b = sink(pres, 1765)
edge(win, split, "no"); edge(split, adpr); edge(adpr, blk); edge(blk, ffa, "yes"); edge(blk, k5a, "no")
edge(ffa, ext); edge(ext, pres); edge(pres, gap, "yes")
edge(pres, k5b, "no")

# ── 6 output tables ──────────────────────────────────────────────────────────────────────────
band(1860, 215, "6 — OUTPUT TABLES   (per polarity)")
unf, fin, filt = row(1896, 66, [(240, "Unfiltered_Results.csv\nevery row + Filter Status"), (215, "Final_Results.csv\nwhat survived"),
                                (240, "Final_Results_Filtered.csv\nanalysis-ready")], "file", gap=50)
norm, assoc, comp = row(1990, 60, [(240, "…_median_Normalised.csv"), (215, "Associated_Spectra.csv\nfrom the search stage"),
                                   (240, "Spectral_Components.csv\nfrom the purity breakdown")], "file", gap=50)
edge(gap, unf, via="target"); edge(gap, fin, via="target"); edge(gap, filt, via="target"); edge(filt, norm)

# ── 7 cross-polarity, report, feedback ───────────────────────────────────────────────────────
band(2095, 330, "7 — CROSS-POLARITY, REPORT AND FEEDBACK TO THE INSTRUMENT   (needs both polarities)")
polq, keep, agree, comb, rep = row(2135, 58, [(180, "class identified in\nboth polarities?"),
                                              (220, "keep the polarity its chemistry\nfavours, carry chain names"),
                                              (150, "cross-polarity\nagreement ρ"),
                                              (250, "Combined_Filtered_Normalised.csv\nthe two blocks on their own scales"),
                                              (130, "QC report")], "proc", gap=28)
dec(polq, 2130); CELLS[comb]["kind"] = "file"; CELLS[rep]["kind"] = "file"
passt = box(CELLS[keep]["x"], 2245, 220, 50, "pass through untouched", "proc")
s7, duty, lists = row(2335, 58, [(110, "search/*.csv"), (250, "duty-cycle diagnosis: never-identified\nprecursors, unfragmented features"),
                                 (250, "exclusion + inclusion lists, collision-\nchecked, import layout → next acquisition")], "proc", gap=40, x0=80)
CELLS[s7]["kind"] = "input"; CELLS[lists]["kind"] = "file"
edge(filt, polq); edge(polq, keep, "yes"); edge(polq, passt, "no"); edge(keep, agree); edge(passt, agree); edge(agree, comb); edge(comb, rep)
edge(s7, duty); edge(duty, lists)

if __name__ == "__main__":
    write(render(), HERE, "figS1", png_width=2400)
