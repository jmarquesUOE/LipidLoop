"""Shared drawing code for the band schematics (Figure 1 compact, Figure S1 full).

Conventions, the same in both: bands in execution order; rounded boxes are steps, diamonds are
questions with a labelled yes AND no exit that each end in a box; green boxes are files written;
grey dashed boxes are sinks (rows that leave the delivered table, all `Unfiltered_Results.csv`);
orange boxes are inputs, repeated at the start of a band that reads an earlier file so that no
arrow travels further than one band. Edges are orthogonal and routed through the corridors
between rows; `check()` refuses a label on a shape or an edge through one.
"""
from __future__ import annotations

import re
from html import escape
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import style  # noqa: E402

CELLS, EDGES = [], []
CANVAS_W = 1240
D = 68


def reset(width: int = 1240) -> None:
    global CANVAS_W
    CELLS.clear(); EDGES.clear(); CANVAS_W = width


def box(x, y, w, h, text, kind, rx=7):
    CELLS.append(dict(x=x, y=y, w=w, h=h, text=text, kind=kind, rx=rx))
    return len(CELLS) - 1


def edge(a, b, label="", dash=False, via=None):
    """`via` = x of a clear column for a long drop, or "target" to drop straight to the corridor
    above the target's row (when the source's own corridor holds a second row of boxes)."""
    EDGES.append((a, b, label, dash, via))


def row(y, h, items, kind, gap=60, x0=None):
    widths = [w for w, _ in items]
    span = sum(widths) + gap * (len(items) - 1)
    x = (CANVAS_W - span) / 2 if x0 is None else x0
    out = []
    for (w, text) in items:
        out.append(box(x, y, w, h, text, kind))
        x += w + gap
    return out


def band(y, h, title):
    return box(30, y, CANVAS_W - 60, h, title, "band")


def dec(idx, y=None):
    CELLS[idx]["kind"] = "dec"; CELLS[idx]["h"] = D
    CELLS[idx]["y"] = y if y is not None else CELLS[idx]["y"] - 5


def sink(under, y, text="→ Unfiltered_Results.csv\nwith its filter status"):
    c = CELLS[under]
    return box(c["x"] + c["w"] / 2 - 95, y, 190, 44, text, "sink")


# ── render ───────────────────────────────────────────────────────────────────────────────────
P = style.PASTEL
STYLE = {
 "band":  ("none",         "#c4c4c4", "#7a7a7a", 11.5, "start",  True,  "5 4"),
 "proc":  (P["sky"],       "#6c8ebf", "#17395e", 10.5, "middle", False, ""),
 "dec":   (P["yellow"],    "#d6b656", "#6b5300", 10,   "middle", False, ""),
 "file":  (P["green"],     "#5f9a6e", "#1e4620", 10.5, "middle", False, ""),
 "input": (P["vermillion"],"#b85450", "#6b1f1c", 10.5, "middle", False, ""),
 "note":  ("none",         "none",    "#6a6a6a", 10.5, "start",  False, ""),
 "sink":  ("#f2f2f2",      "#9a9a9a", "#4a4a4a", 9.5,  "middle", False, "4 3"),
}


def anchor(c, side):
    x, y, w, h = c["x"], c["y"], c["w"], c["h"]
    return {"t": (x + w/2, y), "b": (x + w/2, y + h), "l": (x, y + h/2), "r": (x + w, y + h/2)}[side]


def rows_of(cells):
    mids = sorted({round((c["y"] + c["h"] / 2) / 12) * 12 for c in cells if c["kind"] not in ("band", "note")})
    out = []
    for m in mids:
        members = [c for c in cells if abs(c["y"] + c["h"] / 2 - m) < 22 and c["kind"] not in ("band", "note")]
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
    for (t1, b1), (t2, b2) in zip(ROWS, ROWS[1:]):
        if b1 <= max(y_from, y_to) and t2 >= min(y_from, y_to) and b1 < t2:
            if min(y_from, y_to) <= b1 and max(y_from, y_to) >= t2:
                return (b1 + t2) / 2
    return (y_from + y_to) / 2


def route(ca, cb, via=None):
    mid = lambda c: c["y"] + c["h"] / 2
    if via == "under":
        # same row, loop under it (used when a straight segment would cross a box in between)
        (x1, y1), (x2, y2) = anchor(ca, "b"), anchor(cb, "b")
        my = max(ca["y"] + ca["h"], cb["y"] + cb["h"]) + 22
        return [(x1, y1), (x1, my), (x2, my), (x2, y2)], ((x1 + x2) / 2, my + 3)
    if via == "target":
        # straight down out of the source to the corridor just above the target's row, across, in:
        # for a source whose own corridor is occupied by a second row of boxes in the same band
        down = cb["y"] > ca["y"]
        (x1, y1) = anchor(ca, "b" if down else "t")
        (x2, y2) = anchor(cb, "t" if down else "b")
        top = cb["y"]; prv = max([b for _, b in ROWS if b <= top + 2], default=top - 40)
        c2 = (prv + top) / 2
        return [(x1, y1), (x1, c2), (x2, c2), (x2, y2)], ((x1 + x2) / 2, c2 - 4)
    if via is not None:
        # leave the source vertically into the clear corridor beside its row, travel to the margin
        # column, run along it to the corridor beside the target row, come back over the target and
        # enter it vertically: every horizontal segment lies in a corridor, every vertical one in
        # the margin or directly above/below a box, so nothing is crossed
        down = cb["y"] > ca["y"]
        (x1, y1) = anchor(ca, "b" if down else "t")
        (x2, y2) = anchor(cb, "t" if down else "b")
        def corridor_below(c):
            bot = c["y"] + c["h"]
            nxt = min([t for t, _ in ROWS if t >= bot - 2], default=bot + 40)
            return (bot + nxt) / 2
        def corridor_above(c):
            top = c["y"]
            prv = max([b for _, b in ROWS if b <= top + 2], default=top - 40)
            return (prv + top) / 2
        c1 = corridor_below(ca) if down else corridor_above(ca)
        c2 = corridor_above(cb) if down else corridor_below(cb)
        pts = [(x1, y1), (x1, c1), (via, c1), (via, c2), (x2, c2), (x2, y2)]
        return pts, (x1 + 16, y1 + (13 if down else -6))
    same_row = abs(mid(ca) - mid(cb)) < 22
    if same_row and cb["x"] > ca["x"]:
        (x1, y1), (x2, y2) = anchor(ca, "r"), anchor(cb, "l")
        return [(x1, y1), (x2, y2)], ((x1 + x2) / 2, y1 - 7)
    if same_row:
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


ROWS = []


def render() -> str:
    global ROWS
    ROWS = rows_of(CELLS)
    CANVAS_H = max(c["y"] + c["h"] for c in CELLS) + 20
    svg = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {CANVAS_W} {CANVAS_H}" '
           f'width="{CANVAS_W}" height="{CANVAS_H}" font-family="DejaVu Sans, Arial, sans-serif">',
           '<rect x="0" y="0" width="100%" height="100%" fill="#ffffff"/>',
           '<defs><marker id="a" markerWidth="9" markerHeight="9" refX="8" refY="3" orient="auto">'
           '<path d="M0,0 L0,6 L8,3 z" fill="#7a7a7a"/></marker></defs>']
    labels = []
    for a, b, text, dash, via in EDGES:
        pts, (lx, ly) = route(CELLS[a], CELLS[b], via)
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
        elif c["kind"] != "note":
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
    for lx, ly, text in labels:
        w = 7 + len(text) * 5.6
        svg.append(f'<rect x="{lx - w/2}" y="{ly - 10}" width="{w}" height="14" rx="3" fill="#ffffff" opacity="0.95"/>')
        svg.append(f'<text x="{lx}" y="{ly}" font-size="10" fill="#6a6a6a" text-anchor="middle" '
                   f'font-weight="600">{escape(text)}</text>')
    svg.append("</svg>")
    return "".join(svg)


def check(svg: str) -> dict:
    """The geometric check from docs/diagrams/check_diagram.py: labels and edges clear of shapes."""
    rects = [tuple(map(float, m)) for m in re.findall(r'<rect x="([-\d.]+)" y="([-\d.]+)" width="([\d.]+)" height="([\d.]+)" rx="7"', svg)]
    rects = [r for r in rects if r[2] < 400]
    rhombi = []
    for p in re.findall(r'<polygon points="([\d.]+),([\d.]+) ([\d.]+),([\d.]+) ([\d.]+),([\d.]+) ([\d.]+),([\d.]+)"', svg):
        v = list(map(float, p)); xs, ys = v[0::2], v[1::2]
        rhombi.append(((min(xs)+max(xs))/2, (min(ys)+max(ys))/2, (max(xs)-min(xs))/2, (max(ys)-min(ys))/2))
    in_rect = lambda px, py, r: r[0] <= px <= r[0]+r[2] and r[1] <= py <= r[1]+r[3]
    in_rhom = lambda px, py, d: abs(px-d[0])/d[2] + abs(py-d[1])/d[3] <= 1.0
    hits = lambda px, py: any(in_rect(px, py, r) for r in rects) or any(in_rhom(px, py, d) for d in rhombi)
    labels = re.findall(r'<text x="([\d.]+)" y="([\d.]+)" font-size="10" fill="#6a6a6a"[^>]*>([^<]+)</text>', svg)
    bad_labels = [t for x, y, t in labels
                  if any(hits(float(x)-(7+len(t)*5.6)/2+(7+len(t)*5.6)*i/6, float(y)-10+14*j/3) for i in range(7) for j in range(4))]
    bad_edges = 0
    for d in re.findall(r'<path d="(M [^"]+)"', svg):          # edges only; the marker path is "M0,0 ..."
        pts = [tuple(map(float, p.split())) for p in d.replace("M ", "").split(" L ")]
        for (x1, y1), (x2, y2) in zip(pts, pts[1:]):
            n = 40
            inner = [(x1 + (x2-x1)*k/n, y1 + (y2-y1)*k/n) for k in range(3, n-2)]
            if any(hits(px, py) for px, py in inner):
                bad_edges += 1
                print("  crossing segment", (x1, y1), "->", (x2, y2))
    return {"labels on a shape": bad_labels, "edges through a shape": bad_edges}




def write(svg: str, out_dir: Path, stem: str, png_width: int = 2400) -> dict:
    import cairosvg
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"{stem}.svg").write_text(svg)
    cairosvg.svg2pdf(bytestring=svg.encode(), write_to=str(out_dir / f"{stem}.pdf"))
    cairosvg.svg2png(bytestring=svg.encode(), write_to=str(out_dir / f"{stem}.png"), output_width=png_width)
    result = check(svg)
    print(f"written {out_dir / stem}.svg / .pdf / .png   {result}")
    return result
