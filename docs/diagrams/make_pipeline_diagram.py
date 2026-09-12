"""Generate a draw.io diagram of the pipeline: .raw in, CSVs out, decisions shown.

Editable source for `pipeline_raw_to_csv.html`, which is the one to read — it needs nothing
installed. Kept in step with it deliberately: AGNOSTIC, so no counts, timings or per-study
numbers, because those belong to one analysis and would be wrong for the next. The warnings live
in the HTML beside the figure rather than inside it.
"""
from html import escape

W = []          # cells
_id = [10]

def nid():
    _id[0] += 1
    return f"n{_id[0]}"

def box(x, y, w, h, text, style, parent="1"):
    i = nid()
    W.append(f'<mxCell id="{i}" value="{escape(text)}" style="{style}" vertex="1" parent="{parent}">'
             f'<mxGeometry x="{x}" y="{y}" width="{w}" height="{h}" as="geometry"/></mxCell>')
    return i

def edge(a, b, text="", style="", parent="1"):
    i = nid()
    s = ("edgeStyle=orthogonalEdgeStyle;rounded=0;html=1;fontSize=10;"
         "endArrow=block;endFill=1;strokeColor=#555555;" + style)
    W.append(f'<mxCell id="{i}" value="{escape(text)}" style="{s}" edge="1" parent="{parent}" '
             f'source="{a}" target="{b}"><mxGeometry relative="1" as="geometry"/></mxCell>')

PROC = ("rounded=1;whiteSpace=wrap;html=1;arcSize=12;fillColor=#dae8fc;strokeColor=#6c8ebf;"
        "fontSize=11;verticalAlign=middle;")
DEC  = ("rhombus;whiteSpace=wrap;html=1;fillColor=#fff2cc;strokeColor=#d6b656;fontSize=10;")
FILE = ("shape=note;whiteSpace=wrap;html=1;size=14;fillColor=#d5e8d4;strokeColor=#82b366;"
        "fontSize=10;align=left;spacingLeft=6;verticalAlign=top;")
INPUT= ("shape=parallelogram;perimeter=parallelogramPerimeter;whiteSpace=wrap;html=1;"
        "fixedSize=1;fillColor=#f8cecc;strokeColor=#b85450;fontSize=11;")
NOTE = ("text;html=1;whiteSpace=wrap;fontSize=9;align=left;verticalAlign=top;fontColor=#666666;")
GROUP= ("rounded=0;whiteSpace=wrap;html=1;fillColor=none;strokeColor=#999999;dashed=1;"
        "verticalAlign=top;align=left;spacingLeft=8;fontSize=11;fontStyle=1;fontColor=#666666;")

# ── stage bands
box(20,  20, 1240, 210, "1 — CONVERSION   (reader chosen by detected format, never by suffix)", GROUP)
box(20, 170, 1240, 210, "2 — IDENTIFICATION   (per injection, MS2 vs spectral libraries)", GROUP)
box(20, 400, 1240, 150, "3 — FEATURE DETECTION   (pyOpenMS: detect, align, link)", GROUP)
box(20, 570, 1240, 300, "4 — PEAK FINDER   (join identifications to features, then filter)", GROUP)
box(20, 890, 1240, 250, "5 — POST-FILTERS", GROUP)
box(20,1160, 1240, 230, "6 — OUTPUT TABLES   (per polarity)", GROUP)
box(20,1410, 1240, 250, "7 — CROSS-POLARITY  +  REPORT   (run_study.py, needs BOTH polarities)", GROUP)

# ── 1 conversion
# ⚠ Vendor-aware. Waters `.raw` is a DIRECTORY and Thermo `.raw` is a FILE — same suffix, two
# vendors, two readers — so the reader is chosen from what the path IS.
raw   = box(60,  55, 190, 60, "vendor data\n.raw .d .wiff", INPUT)
fmt   = box(300, 50, 170, 70, "detect\nformat ?", DEC)
trfp  = box(520, 30, 200, 52, "ThermoRawFileParser\nThermo .raw", PROC)
mscv  = box(520, 88, 200, 52, "msconvert (wine)\nAgilent Waters Sciex", PROC)
tims  = box(520,146, 200, 52, "alphatims + writer\nBruker timsTOF", PROC)
mzml  = box(770, 88, 190, 60, ".mzML\ncached per STUDY", FILE)
prof  = box(1000, 50, 170, 70, "profile ?\ncentroid once", DEC)
cache = box(1000,140, 170, 60, "mzML newer\nthan source ?", DEC)
edge(raw, fmt)
edge(fmt, trfp); edge(fmt, mscv); edge(fmt, tims)
edge(trfp, mzml); edge(mscv, mzml); edge(tims, mzml)
edge(mzml, prof); edge(mzml, cache)

# ── 2 identification
ms2   = box(60, 215, 170, 60, "MS2 scans\nread from mzML", PROC)
artef = box(280, 210, 175, 70, "impossible\nmass defect ?", DEC)
srch  = box(500, 215, 190, 60, "library search\nsearch.py + .msp", PROC)
offs  = box(740, 215, 180, 60, "mass_offset_ppm\napplied to search", PROC)
pur   = box(970, 210, 170, 70, "purity ≥ 75 % ?", DEC)
searchcsv = box(500, 300, 190, 60, "search/*.csv\none per injection", FILE)
edge(mzml, ms2); edge(ms2, artef); edge(artef, srch, "no — keep")
edge(artef, ms2, "yes — strip")
edge(srch, offs); edge(offs, pur); edge(srch, searchcsv)

# ── 3 features
feat  = box(60, 445, 190, 60, "feature detection\nfeatures.py", PROC)
align = box(300, 445, 190, 60, "RT alignment\nLOESS per sample", PROC)
link  = box(540, 445, 190, 60, "link across\ninjections", PROC)
noise = box(780, 445, 190, 60, "noise floor\nmeasured, not assumed", PROC)
edge(mzml, feat); edge(feat, align); edge(align, link); edge(link, noise)

# ── 4 peak finder
join  = box(60, 615, 180, 60, "join MS2 → features\npeakfinder.py", PROC)
score = box(280, 610, 175, 70, "dot ≥ 500 and\nrev dot ≥ 700 ?", DEC)
name  = box(500, 615, 200, 60, "name the row\nmolecular OR sum", PROC)
rtmod = box(740, 615, 200, 60, "per-class retention\nmodel  RT ~ C + DB", PROC)
rtout = box(990, 610, 180, 70, "|z| > 3σ from\nits class model ?", DEC)
addsw = box(60, 730, 180, 60, "adduct / dimer /\nisotope sweep", PROC)
insrc = box(280, 730, 175, 60, "in-source fragment\nremoval", PROC)
redun = box(500, 730, 200, 60, "redundant\nidentification", PROC)
rtwin = box(740, 725, 200, 70, "class RT window\n(OFF — clips long chains)", DEC)
edge(noise, join); edge(pur, join); edge(join, score); edge(score, name, "yes")
edge(name, rtmod); edge(rtmod, rtout)
edge(rtout, addsw, "no — keep"); edge(addsw, insrc); edge(insrc, redun); edge(redun, rtwin)

# ── 5 post filters
split = box(60,  935, 190, 60, "split-peak merge\npooled-QC precision", PROC)
adpr  = box(300, 935, 190, 60, "adduct-pair removal\nNa−H vs +2C+3DB", PROC)
blank = box(540, 935, 190, 60, "blank filter\nmean × 3, in ≥ 2", PROC)
ffa   = box(780, 935, 190, 60, "free fatty acids\nRT ~ C + DB surface", PROC)
ext   = box(1020, 935, 190, 60, "retention-model\nextension (no MS2)", PROC)
pres  = box(60, 1040, 190, 70, "detected in ≥ 70 %\nof ANY one group ?", DEC)
gap   = box(300, 1045, 190, 60, "gap filling\nre-integrate the raw", PROC)
edge(rtwin, split, "kept"); edge(split, adpr); edge(adpr, blank); edge(blank, ffa); edge(ffa, ext)
edge(ext, pres); edge(pres, gap, "yes")

# ── 6 outputs
unf  = box(60, 1215, 210, 90, "Unfiltered_Results.csv\nEVERY compound group\n+ Filter Status = why it went", FILE)
fin  = box(300, 1215, 210, 90, "Final_Results.csv\nwhat survived the filters", FILE)
filt = box(540, 1215, 210, 90, "Final_Results_Filtered.csv\nANALYSIS-READY\npresence-filtered + gap-filled", FILE)
norm = box(780, 1215, 210, 90, "…_median_Normalised.csv\nper-sample loading removed", FILE)
assoc= box(1020, 1215, 210, 40, "Associated_Spectra.csv", FILE)
comp = box(1020, 1265, 210, 40, "Spectral_Components.csv", FILE)
edge(pres, unf, "no → recorded"); edge(gap, filt); edge(gap, fin); edge(filt, norm)
edge(searchcsv, assoc); edge(pur, comp)

# ── 7 cross polarity
polq = box(60, 1455, 200, 70, "class seen in\nBOTH polarities ?", DEC)
polf = box(310, 1460, 200, 60, "keep it in the polarity\nthat measures it better", PROC)
agree= box(560, 1460, 200, 60, "agreement ρ\n⚠ BEFORE rewriting", PROC)
comb = box(810, 1455, 200, 70, "Combined_Filtered_\nNormalised.csv\nMode + Scale + Agreement", FILE)
rep  = box(1050, 1460, 180, 60, "QC report\nHTML + PDF", FILE)
edge(filt, polq); edge(polq, polf, "yes"); edge(polq, comb, "no — passes through")
edge(polf, agree); edge(agree, comb); edge(comb, rep); edge(norm, polq)

xml = ('<mxfile host="app.diagrams.net" agent="lipidomics-pipeline">'
       '<diagram id="pipeline" name="Lipidomics pipeline">'
       '<mxGraphModel dx="1400" dy="900" grid="1" gridSize="10" guides="1" tooltips="1" '
       'connect="1" arrows="1" fold="1" page="1" pageScale="1" pageWidth="1300" pageHeight="1700" '
       'math="0" shadow="0"><root>'
       '<mxCell id="0"/><mxCell id="1" parent="0"/>'
       + "".join(W) +
       '</root></mxGraphModel></diagram></mxfile>')

out = str(Path(__file__).resolve().parent / "pipeline_raw_to_csv.drawio")
open(out, "w").write(xml)
print(f"wrote {out}  ({len(W)} cells, {len(xml)} bytes)")
