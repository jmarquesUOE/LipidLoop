"""Figure 1: the pipeline as a compact band schematic. Drawing code and conventions live in
`manuscript/figures/schematic_lib.py`; this file is the layout only.

    .venv/bin/python manuscript/figures/fig1/schematic.py
"""
from __future__ import annotations

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
from schematic_lib import CELLS, band, box, dec, edge, render, reset, row, sink, write  # noqa: E402

reset(1240)
VR, VL = 1190, 50

# ── layout ───────────────────────────────────────────────────────────────────────────────────
# Every diamond asks a question and has a labelled yes and a labelled no, each ending in a box.
# No arrow travels further than one band: where a later band consumes an earlier file, that file
# is repeated as a small input box at the start of the band; rows a decision removes end in a
# small sink box under the diamond, all of them the same file, Unfiltered_Results.csv.



# 1 conversion
band(20, 130, "1 — CONVERSION   (reader chosen from the file itself, one mzML cache per study)")
raw, fmt, conv, mzml = row(56, 58, [(190, "vendor data\n.raw  .d  .wiff"), (200, "detect format,\nchoose the reader"),
                                    (300, "convert: ThermoRawFileParser |\nmsconvert | alphatims"),
                                    (190, ".mzML\ncached, atomic write")], "proc", gap=45)
CELLS[raw]["kind"] = "input"; CELLS[mzml]["kind"] = "file"
edge(raw, fmt); edge(fmt, conv); edge(conv, mzml)

# 2 identification (per injection)
band(170, 235, "2 — IDENTIFICATION   (per injection, MS2 against target AND decoy libraries)")
ms2, art, srch, pur, chain = row(210, 58, [(130, "MS2 scans\nfrom the mzML"), (190, "artefact present in\n≥ 95% of spectra?"),
                                          (220, "score every match\nLipiDex algorithm, exact\ntarget + decoy libraries"),
                                          (180, "fatty-acid\npurity ≥ 75%?"),
                                          (170, "name at\nchain level")], "proc", gap=36)
dec(art, 205); dec(pur, 205)
strip = box(CELLS[art]["x"] + 5, 318, 180, 50, "strip the artefact peak", "proc")
sumc = box(CELLS[pur]["x"], 318, 180, 50, "name at\nsum composition", "proc")
sres = box(CELLS[chain]["x"], 318, 170, 50, "search/*.csv\nper injection", "file")
edge(mzml, ms2); edge(ms2, art); edge(art, strip, "yes"); edge(strip, srch); edge(art, srch, "no")
edge(srch, pur); edge(pur, chain, "yes"); edge(pur, sumc, "no"); edge(chain, sres); edge(sumc, sres)

# 3 feature detection and calibration
band(425, 165, "3 — FEATURE DETECTION   (+ per-file mass calibration as a second pass)")
mz3, feat, drift, cal = row(462, 58, [(120, ".mzML"), (280, "detect features, align retention,\nlink injections (pyOpenMS)"),
                                      (190, "per-file drift\n≥ 2 ppm?"),
                                      (270, "write calibrated mzML\nsecond pass on corrected files")], "proc", gap=45)
CELLS[mz3]["kind"] = "input"; dec(drift, 457)
edge(mz3, feat); edge(feat, drift); edge(drift, cal, "yes"); edge(cal, feat, "", True)

# 4 peak finder
band(610, 235, "4 — PEAK FINDER   (join identifications to features, then filter)")
s4, join, thr, name, rtm, offm, pff = row(647, 58, [(110, "search/*.csv"), (130, "join MS2\n→ features"),
                                                    (170, "dot ≥ 500 and\nreverse dot ≥ 700?"),
                                                    (110, "name the row"), (140, "fit per-class\nretention model"),
                                                    (150, "off its class\nmodel?"),
                                                    (200, "adduct, dimer, isotope,\nin-source fragment,\nredundant name (one filter\nstatus each)")], "proc", gap=24)
CELLS[s4]["kind"] = "input"; dec(thr, 642); dec(offm, 642)
k4a = sink(thr, 755); k4b = sink(offm, 755)
edge(drift, join, "no"); edge(s4, join)
edge(join, thr); edge(thr, name, "yes"); edge(thr, k4a, "no"); edge(name, rtm); edge(rtm, offm)
edge(offm, pff, "no"); edge(offm, k4b, "yes")

# 5 post-filters
band(865, 235, "5 — POST-FILTERS")
dup, blk, pres, gap = row(902, 58, [(250, "adduct-pair removal,\nsplit-peak merge"), (180, "sample ≥ 3 × blank?"),
                                    (200, "detected in ≥ 60%\nof one group?"), (210, "gap filling\nre-integrated from mzML")], "proc", gap=50)
dec(blk, 897); dec(pres, 897)
k5a = sink(blk, 1010); k5b = sink(pres, 1010)
edge(pff, dup); edge(dup, blk); edge(blk, pres, "yes"); edge(blk, k5a, "no"); edge(pres, gap, "yes"); edge(pres, k5b, "no")

# 6 output tables (per polarity)
band(1120, 125, "6 — OUTPUT TABLES   (per polarity)")
unf, fin, filt, assoc = row(1155, 58, [(230, "Unfiltered_Results.csv\nevery row + Filter Status"),
                                       (190, "Final_Results.csv\nwhat survived"),
                                       (220, "Final_Results_Filtered.csv\nanalysis-ready"),
                                       (220, "Associated_Spectra.csv\nevidence per row, from the search")], "file", gap=40)
edge(gap, unf, via="target"); edge(gap, fin, via="target"); edge(gap, filt, via="target")

# 7 cross-polarity, report, feedback
band(1265, 235, "7 — CROSS-POLARITY, REPORT AND FEEDBACK TO THE INSTRUMENT")
polq, keep, comb, duty, lists = row(1305, 58, [(180, "class identified in\nboth polarities?"),
                                               (220, "keep the polarity its chemistry\nfavours, carry chain names"),
                                               (180, "combined table\n+ QC report"),
                                               (210, "duty-cycle diagnosis\nnever-identified precursors,\nunfragmented features"),
                                               (220, "exclusion + inclusion lists\ncollision-checked, import layout\n→ loaded into the next acquisition")], "proc", gap=26)
dec(polq, 1300)
passt = box(CELLS[keep]["x"], 1415, 220, 50, "pass through untouched", "proc")
CELLS[comb]["kind"] = "file"; CELLS[lists]["kind"] = "file"
edge(filt, polq); edge(polq, keep, "yes"); edge(polq, passt, "no"); edge(keep, comb); edge(passt, comb)
edge(assoc, duty); edge(duty, lists)




if __name__ == "__main__":
    write(render(), HERE, "fig1")
