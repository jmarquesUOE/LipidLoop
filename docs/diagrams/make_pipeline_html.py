"""Standalone HTML+SVG of the pipeline. No install, no network, opens in any browser.

Deliberately AGNOSTIC: no counts, no timings, no per-study numbers. Those belong to one analysis
and would be wrong for the next one — a study of a thousand injections does not take the same time
as a study of twenty. The diagram states what the pipeline DOES; the notes below it state what to
watch for. Neither claims a measurement.
"""
from html import escape

CELLS, EDGES = [], []
CANVAS_W = 1160

def box(x, y, w, h, text, kind, rx=7):
    CELLS.append(dict(x=x, y=y, w=w, h=h, text=text, kind=kind, rx=rx))
    return len(CELLS) - 1

def edge(a, b, label="", dash=False):
    EDGES.append((a, b, label, dash))

def row(y, h, items, kind, gap=70, total=None):
    """Lay a row of boxes out CENTRED on the canvas."""
    widths = [w for w, _ in items]
    span = sum(widths) + gap * (len(items) - 1)
    x = (CANVAS_W - span) / 2
    out = []
    for (w, text) in items:
        out.append(box(x, y, w, h, text, kind))
        x += w + gap
    return out

BAND_H = {}
def band(y, h, title):
    return box(30, y, CANVAS_W - 60, h, title, "band")

# ── 1 conversion
#
# ⚠ Vendor-aware, not Thermo-only. The reader is chosen from what the path IS, never from its
# suffix: Waters `.raw` is a DIRECTORY and Thermo `.raw` is a FILE — same suffix, two vendors, two
# readers. An earlier version of this diagram showed a single ThermoRawFileParser box, which was
# true of one vendor in five.
band(20, 262, "1 — CONVERSION   (reader chosen by detected format, never by suffix)")
raw, fmt = row(52, 58, [(230, "vendor data\n.raw .d .wiff, Pos/ and Neg/"),
                        (180, "detect format")], "mixed", gap=70)
CELLS[raw]["kind"] = "input"; CELLS[fmt]["kind"] = "dec"; CELLS[fmt]["h"] = 68
CELLS[fmt]["y"] = 47
trfp, mscv, tims = row(136, 52, [(230, "ThermoRawFileParser\nThermo .raw"),
                                 (250, "msconvert (wine)\nAgilent, Waters, Sciex"),
                                 (230, "alphatims + writer\nBruker timsTOF")], "proc", gap=40)
# ⚠ Centred, not aligned to the first converter. Left-aligned it sits directly above `MS2 scans`,
# and the drop from here to feature detection then runs straight through that box — which the
# geometric check catches and a reader would simply misread.
mzml = box(325, 210, 250, 58,
           ".mzML\ncached per STUDY, atomic write", "file")
prof = box(CELLS[mzml]["x"] + 290, 205, 180, 68,
           "profile data?\nyes → centroid once", "dec")
edge(raw, fmt); edge(fmt, trfp); edge(fmt, mscv); edge(fmt, tims)
edge(trfp, mzml); edge(mscv, mzml); edge(tims, mzml); edge(mzml, prof)
cache = box(CELLS[prof]["x"] + 250, 205, 180, 68, "mzML newer than source?\nno → convert", "dec")
edge(prof, cache)

# ── 2 identification
#
# ⚠ Both halves of the search are drawn. The false-discovery rate is a headline claim and a reader
# cannot judge it from a single "library search" box: TARGET and DECOY are searched together
# against the same spectra, which is what makes the decoy count a rate rather than a tally.
#
# The two decoy constructions differ and the difference matters — chain-shuffled where the score
# depends on the chains, head-group-displaced where it does not, because a chain shuffle in a class
# whose scored intensity is chain-independent produces a copy of its own target rather than a decoy.
band(304, 268, "2 — IDENTIFICATION   (per injection: MS2 against target AND decoy libraries)")
ms2, artef = row(340, 58, [(180, "MS2 scans"), (190, "artefact screen\nyes → strip")],
                 "proc", gap=110)
CELLS[artef]["kind"] = "dec"; CELLS[artef]["h"] = 68; CELLS[artef]["y"] = 335
# ⚠ Placed explicitly, not auto-centred. Centred, this row starts at x=205 — exactly where the
# mzML → detect-features edge drops — and the edge then runs through `target libraries`.
tgt  = box(280, 424, 210, 54, "target libraries", "proc")
dec1 = box(530, 424, 230, 54, "decoy — chain shuffle", "proc")
dec2 = box(800, 424, 230, 54, "decoy — head group", "proc")
srch, offs, pur = row(506, 58, [(190, "score every match"), (190, "mass offset applied"),
                                (180, "purity ≥ threshold?")], "proc", gap=50)
CELLS[pur]["kind"] = "dec"; CELLS[pur]["h"] = 68; CELLS[pur]["y"] = 501
# ⚠ On the purity row, right of it, not on the library row: at x=940 on the library row it sat
# on top of `decoy — head group` (800-1030) — a real overlap the checker does not test for.
scsv = box(CANVAS_W - 220, 506, 190, 58, "search/*.csv\nper injection", "file")
edge(mzml, ms2); edge(ms2, artef)
edge(artef, tgt, "keep"); edge(artef, dec1); edge(artef, dec2)
edge(tgt, srch); edge(dec1, srch); edge(dec2, srch)
edge(srch, offs); edge(offs, pur); edge(pur, scsv, "yes")

# ── 3 features
# ⚠ Linking is where a drifting mass axis does its damage, so the calibration loop is drawn here
# rather than left implicit. Linking uses a fixed tolerance, and a file whose mass axis sits outside
# it can never join the consensus group — its features form their own, leaving a zero in the main
# row for every injection concerned. The loop measures each file's offset against the standards or
# the batch consensus, writes corrected copies, and links again.
#
# The dashed return edge is the second pass. It is a LOOP, not a step, which is the whole point:
# the offset cannot be measured until the features have been linked once.
band(594, 200, "3 — FEATURE DETECTION   (+ per-file mass calibration → a SECOND PASS)")
feat, align, link, noise = row(626, 58, [(180, "detect features"), (180, "align retention"),
                                         (180, "link injections"), (180, "measure noise floor")], "proc")
edge(mzml, feat); edge(feat, align); edge(align, link); edge(link, noise)
# ⚠ Right of x=325: two edges into the peak finder drop there, and a box in the way is a
# crossing the geometric check rejects.
drift = box(560, 706, 210, 68, "per-file offset\nbeyond tolerance?", "dec")
recal = box(830, 712, 230, 56, "write calibrated mzML\nprecursors too", "proc")
edge(link, drift); edge(drift, recal)
# ⚠ The second pass re-enters at the SEARCH, not at linking. `calibrate_files.apply` shifts MS1,
# MS2 and precursors alike — correcting the spectra and leaving the precursors behind would push
# every MS2 off the peak it was taken from — so once the precursors move, identification changes
# too. Pass two is the whole pipeline over corrected files.
#
# Said in the box rather than drawn as an edge: a back-arrow from here to band 2 spans four bands,
# and the router finds no clear column — it crossed `decoy — head group`, `search/*.csv` and
# `measure noise floor`. A long arrow nobody can follow is worse than a sentence.
CELLS[recal]["text"] = "write calibrated mzML\n→ re-run from the search"

# ── 4 peak finder
band(818, 250, "4 — PEAK FINDER   (join identifications to features, then filter)")
join, score, name = row(852, 58, [(190, "join MS2 → features"), (180, "scores above\nthresholds?"),
                                  (190, "name the row\nmolecular OR sum")], "proc")
CELLS[score]["kind"] = "dec"; CELLS[score]["h"] = 68; CELLS[score]["y"] = 847
rtmod, rtout = row(936, 58, [(210, "fit per-class\nretention model"), (180, "off its class\nmodel?")], "proc")
CELLS[rtout]["kind"] = "dec"; CELLS[rtout]["h"] = 68; CELLS[rtout]["y"] = 931
addsw, insrc, redun, rtwin = row(1018, 58,
    [(170, "adduct / dimer /\nisotope"), (170, "in-source\nfragment"),
     (170, "redundant\nidentification"), (180, "class RT window\n(off by default)")], "proc")
CELLS[rtwin]["kind"] = "dec"; CELLS[rtwin]["h"] = 68; CELLS[rtwin]["y"] = 1013
edge(noise, join); edge(pur, join); edge(join, score); edge(score, name, "yes")
edge(name, rtmod); edge(rtmod, rtout); edge(rtout, addsw, "no")
edge(addsw, insrc); edge(insrc, redun); edge(redun, rtwin)

# ── 5 post filters
band(1090, 210, "5 — POST-FILTERS")
split, adpr, blank = row(1122, 58, [(190, "split-peak merge"), (190, "adduct-pair removal"),
                                   (170, "blank filter")], "proc")
ffa, ext, pres, gap_ = row(1202, 58, [(180, "free fatty acids"), (190, "retention-model\nextension"),
                                     (190, "present in a group?"), (170, "gap filling")], "proc")
CELLS[pres]["kind"] = "dec"; CELLS[pres]["h"] = 68; CELLS[pres]["y"] = 1197
edge(rtwin, split, "kept"); edge(split, adpr); edge(adpr, blank); edge(blank, ffa)
edge(ffa, ext); edge(ext, pres); edge(pres, gap_, "yes")

# ── 6 outputs
band(1312, 208, "6 — OUTPUT TABLES   (per polarity)")
unf, fin, filt = row(1346, 76, [(230, "Unfiltered_Results.csv\nevery group + Filter Status"),
                                (215, "Final_Results.csv\nwhat survived"),
                                (240, "Final_Results_Filtered.csv\nANALYSIS-READY")], "file")
norm, assoc, comp = row(1438, 66, [(240, "…_median_Normalised.csv"),
                                   (215, "Associated_Spectra.csv\nfrom the search stage"),
                                   (230, "Spectral_Components.csv\nfrom the purity breakdown")], "file")
edge(pres, unf, "no", True); edge(gap_, fin); edge(gap_, filt); edge(filt, norm)


# ── 7 cross polarity
band(1534, 150, "7 — CROSS-POLARITY + REPORT   (needs BOTH polarities)")
polq, polf, agree = row(1570, 58, [(200, "class in BOTH\npolarities?"),
                                   (210, "keep the polarity that\nmeasures it better"),
                                   (180, "agreement ρ")], "proc")
CELLS[polq]["kind"] = "dec"; CELLS[polq]["h"] = 68; CELLS[polq]["y"] = 1565
comb, rep = row(1652, 56, [(250, "Combined_Filtered_Normalised.csv"), (170, "QC report")], "file")
edge(norm, polq); edge(polq, polf, "yes")
edge(polq, comb, "no", True); edge(polf, agree); edge(agree, comb); edge(comb, rep)

CANVAS_H = 1734

# ─────────────────────────────────────────────── render
STYLE = {
 "band":  ("none",    "#c4c4c4", "#7a7a7a", 11.5, "start",  True,  "5 4"),
 "proc":  ("#dae8fc", "#6c8ebf", "#17395e", 11,   "middle", False, ""),
 "dec":   ("#fff2cc", "#d6b656", "#6b5300", 10.5, "middle", False, ""),
 "file":  ("#d5e8d4", "#82b366", "#1e4620", 10.5, "middle", False, ""),
 "input": ("#f8cecc", "#b85450", "#6b1f1c", 11,   "middle", False, ""),
 "mixed": ("#dae8fc", "#6c8ebf", "#17395e", 11,   "middle", False, ""),
}
def anchor(c, side):
    x, y, w, h = c["x"], c["y"], c["w"], c["h"]
    return {"t": (x + w/2, y), "b": (x + w/2, y + h),
            "l": (x, y + h/2), "r": (x + w, y + h/2)}[side]


def rows_of(cells):
    """Cluster shapes into rows by centre, so a corridor can clear the TALLEST of them.

    Diamonds are taller than the boxes beside them. Deriving a corridor from one box's bottom
    edge therefore put it straight through the diamond next door — the source of every remaining
    crossing.
    """
    mids = sorted({round((c["y"] + c["h"] / 2) / 12) * 12 for c in cells if c["kind"] != "band"})
    out = []
    for m in mids:
        members = [c for c in cells if abs(c["y"] + c["h"] / 2 - m) < 22 and c["kind"] != "band"]
        if members:
            out.append((min(c["y"] for c in members), max(c["y"] + c["h"] for c in members)))
    merged = []
    for top, bot in out:
        if merged and top <= merged[-1][1] + 4:
            merged[-1] = (merged[-1][0], max(merged[-1][1], bot))
        else:
            merged.append((top, bot))
    return merged


def corridor(y_from, y_to):
    """Clear horizontal band between the two rows those y values fall in."""
    for (t1, b1), (t2, b2) in zip(ROWS, ROWS[1:]):
        if b1 <= max(y_from, y_to) and t2 >= min(y_from, y_to) and b1 < t2:
            if min(y_from, y_to) <= b1 and max(y_from, y_to) >= t2:
                return (b1 + t2) / 2
    return (y_from + y_to) / 2


def route(ca, cb):
    """Orthogonal path that leaves the corridors between rows, never crossing a box.

    Same row and adjacent -> a straight segment in the gap, which is where the label goes.
    Otherwise -> drop into the corridor below the source row, travel there, then enter the
    target from above. The corridor is empty by construction because rows are laid out with
    vertical clearance, so nothing is crossed.
    """
    mid = lambda c: c["y"] + c["h"] / 2
    same_row = abs(mid(ca) - mid(cb)) < 22
    if same_row and cb["x"] > ca["x"]:
        (x1, y1), (x2, y2) = anchor(ca, "r"), anchor(cb, "l")
        top = min(ca["y"], cb["y"])
        above = max([b for _, b in ROWS if b <= top + 2], default=top - 26)
        return [(x1, y1), (x2, y2)], ((x1 + x2) / 2, (above + top) / 2 + 3)
    if same_row:                                   # right-to-left, loop under the row
        (x1, y1), (x2, y2) = anchor(ca, "b"), anchor(cb, "b")
        my = max(ca["y"] + ca["h"], cb["y"] + cb["h"]) + 20
        return [(x1, y1), (x1, my), (x2, my), (x2, y2)], ((x1 + x2) / 2, my + 3)
    down = cb["y"] > ca["y"]
    (x1, y1) = anchor(ca, "b" if down else "t")
    (x2, y2) = anchor(cb, "t" if down else "b")
    my = corridor(ca["y"] + ca["h"] / 2, cb["y"] + cb["h"] / 2)
    if abs(x1 - x2) < 6:
        return [(x1, y1), (x2, y2)], ((x1 + x2) / 2, (y1 + y2) / 2 - 4)
    return [(x1, y1), (x1, my), (x2, my), (x2, y2)], ((x1 + x2) / 2, my + 3)


ROWS = rows_of(CELLS)

svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {CANVAS_W} {CANVAS_H}" '
       f'style="width:100%;height:auto;background:#fff" font-family="system-ui,-apple-system,sans-serif">',
       '<defs><marker id="a" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto">'
       '<path d="M0,0 L0,6 L8,3 z" fill="#7a7a7a"/></marker></defs>']

labels = []
for a, b, text, dash in EDGES:
    pts, (lx, ly) = route(CELLS[a], CELLS[b])
    d = "M " + " L ".join(f"{x} {y}" for x, y in pts)
    svg.append(f'<path d="{d}" fill="none" stroke="#9a9a9a" stroke-width="1.3"'
               + (' stroke-dasharray="5 4"' if dash else '') + ' marker-end="url(#a)"/>')
    if text:
        labels.append((lx, ly, text))

for c in CELLS:
    fill, stroke, fg, fs, align, bold, dash = STYLE[c["kind"]]
    da = f' stroke-dasharray="{dash}"' if dash else ""
    if c["kind"] == "dec":
        cx, cy = c["x"] + c["w"]/2, c["y"] + c["h"]/2
        svg.append(f'<polygon points="{cx},{c["y"]} {c["x"]+c["w"]},{cy} {cx},{c["y"]+c["h"]} '
                   f'{c["x"]},{cy}" fill="{fill}" stroke="{stroke}" stroke-width="1.2"/>')
    else:
        svg.append(f'<rect x="{c["x"]}" y="{c["y"]}" width="{c["w"]}" height="{c["h"]}" '
                   f'rx="{c["rx"]}" fill="{fill}" stroke="{stroke}" stroke-width="1.2"{da}/>')
    lines = c["text"].split("\n")
    if align == "start":
        tx, ty, ta = c["x"] + 12, c["y"] + 17, "start"
    else:
        tx, ta = c["x"] + c["w"]/2, "middle"
        ty = c["y"] + c["h"]/2 - (len(lines) - 1) * fs * 0.60 + fs * 0.36
    weight = " font-weight='600'" if (bold or c["kind"] == "file") else ""
    for i, ln in enumerate(lines):
        svg.append(f'<text x="{tx}" y="{ty + i*fs*1.3}" font-size="{fs}" fill="{fg}" '
                   f'text-anchor="{ta}"{weight}>{escape(ln)}</text>')

# labels last, on a white plate, so a line can never cut through "yes" or "no"
for lx, ly, text in labels:
    w = 7 + len(text) * 5.6
    svg.append(f'<rect x="{lx - w/2}" y="{ly - 10}" width="{w}" height="14" rx="3" '
               f'fill="#ffffff" stroke="none" opacity="0.95"/>')
    svg.append(f'<text x="{lx}" y="{ly}" font-size="10" fill="#6a6a6a" text-anchor="middle" '
               f'font-weight="600">{escape(text)}</text>')
svg.append("</svg>")

WARN = [
 ("Purity decides the resolution of the name, but does not gate it",
  "The purity threshold decides whether a row is reported at molecular or at sum composition. It "
  "scores the <em>winning</em> identification — it is not a gate on the name that identification "
  "carries, and library entries are routinely written at chain resolution. A row can therefore be "
  "chain-resolved on no chain evidence at all, which is what the <code>Chain Evidence</code> column "
  "reports."),
 ("A class with no usable retention model is exempt from the retention check",
  "The model is fitted per class. A class with too few members, or one that is not a homologous "
  "series, produces no model — and every row in it passes the retention check <strong>by never "
  "taking it</strong>. Those rows carry a blank model error, which reads as though nothing was "
  "wrong. <code>RT Model R2</code> says whether the class model is worth believing."),
 ("The class retention window is off by default",
  "It removes long-chain species by construction rather than by evidence, so the per-class "
  "retention model does that job instead."),
 ("Presence is judged per group, never globally",
  "A row is kept if it was detected in enough of <em>any one</em> group. Judged across the whole "
  "study instead, a species present in only one arm — often the interesting one — is deleted for "
  "being absent from the other."),
 ("Gaps are measured, not imputed",
  "Where the peak finder found nothing, the raw file is re-integrated at that mass and retention "
  "time, so a filled value is a measurement at that coordinate rather than a number inferred from "
  "other samples. A group with no detections at all keeps its zeros: presence and absence are never "
  "erased."),
 ("A class in both polarities with no rule stops the run",
  "The assignment table is checked against every run. A class identified in both polarities but "
  "absent from the table halts processing rather than being guessed, because a wrong guess deletes "
  "real lipids and nothing downstream shows it."),
 ("The two polarities are not two measurements of one thing",
  "Different ions with different ionisation efficiencies. They are never averaged or merged; the "
  "combined matrix carries the mode and scale of each block to say so."),
 ("Conversion is cached",
  "An mzML newer than its raw file is reused, so reprocessing skips the conversion stage entirely. "
  "Conversion is written atomically — an interrupted convert must not leave a short file that the "
  "freshness check then trusts."),
]

rows = "".join(f"<tr><td><strong>{escape(t)}</strong><br><span>{b}</span></td></tr>" for t, b in WARN)

html = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Lipidomics pipeline — .raw to delivered CSVs</title>
<style>
 body{{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;margin:0;padding:26px 18px;
      background:#f4f4f6;color:#1b1b1b;line-height:1.6}}
 .wrap{{max-width:1240px;margin:0 auto;background:#fff;padding:28px 34px 38px;
        border-radius:10px;box-shadow:0 1px 4px rgba(0,0,0,.10)}}
 h1{{font-size:20px;margin:0 0 4px}} h2{{font-size:15px;margin:30px 0 10px}}
 p.sub{{margin:0 0 18px;color:#666;font-size:13px}}
 .key{{display:flex;gap:22px;flex-wrap:wrap;margin:0 0 18px;font-size:12px;color:#555;
       justify-content:center}}
 .key span{{display:flex;align-items:center;gap:7px}}
 .sw{{width:15px;height:12px;border-radius:3px;display:inline-block;border:1px solid #aaa}}
 .fig{{overflow-x:auto;border:1px solid #ececec;border-radius:6px;padding:8px;text-align:center}}
 table{{border-collapse:collapse;font-size:13px;width:100%}}
 td{{border-bottom:1px solid #eee;padding:9px 4px;vertical-align:top}}
 td span{{color:#444}}
 code{{background:#f2f2f4;padding:1px 5px;border-radius:3px;font-size:12px}}
</style></head><body><div class="wrap">
<h1>Lipidomics pipeline — <code>.raw</code> to the delivered CSVs</h1>
<p class="sub">Seven bands, in the order the pipeline executes them. Yellow diamonds are decisions;
each one changes what happens to a row. The diagram carries no counts or timings — those belong to
one study and would be wrong for the next.</p>
<div class="key">
 <span><i class="sw" style="background:#f8cecc"></i>input</span>
 <span><i class="sw" style="background:#dae8fc"></i>processing step</span>
 <span><i class="sw" style="background:#fff2cc"></i>decision</span>
 <span><i class="sw" style="background:#d5e8d4"></i>file written</span>
</div>
<div class="fig">{''.join(svg)}</div>

<p class="foot">The reading notes — what each decision costs when it is wrong — are in
<code>PIPELINE_FIGURE_LEGEND.md</code>, kept separate so this page can be used as a
manuscript figure unaltered.</p>
</div></body></html>"""

out = str(Path(__file__).resolve().parent / "pipeline_raw_to_csv.html")
open(out, "w").write(html)
print(f"wrote {out}  ({len(html)//1024} KB, self-contained)")
