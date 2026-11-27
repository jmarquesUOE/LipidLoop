"""Geometric check of the diagram. Diamonds are tested as RHOMBI, not bounding boxes —
a label just past a vertex is in clear space, and treating it as a box invents collisions."""
import re, sys
s = open("docs/diagrams/pipeline_raw_to_csv.html").read()
svg = re.search(r"<svg.*?</svg>", s, re.S).group(0)
W, H = map(float, re.search(r'viewBox="0 0 ([\d.]+) ([\d.]+)"', svg).groups())

rects = [(float(a),float(b),float(c),float(d)) for a,b,c,d in
         re.findall(r"<rect x=\"([-\d.]+)\" y=\"([-\d.]+)\" width=\"([\d.]+)\" height=\"([\d.]+)\" rx=\"[67]\"", svg)]
rects = [r for r in rects if r[2] < 400]                       # drop the band frames
rhombi = []
for p in re.findall(r'<polygon points="([\d.]+),([\d.]+) ([\d.]+),([\d.]+) ([\d.]+),([\d.]+) ([\d.]+),([\d.]+)"', svg):
    v = list(map(float, p))
    xs, ys = v[0::2], v[1::2]
    cx, cy = (min(xs)+max(xs))/2, (min(ys)+max(ys))/2
    rhombi.append((cx, cy, (max(xs)-min(xs))/2, (max(ys)-min(ys))/2))

def in_rect(px, py, r):  return r[0] <= px <= r[0]+r[2] and r[1] <= py <= r[1]+r[3]
def in_rhom(px, py, d):  return abs(px-d[0])/d[2] + abs(py-d[1])/d[3] <= 1.0
def hits(px, py):
    return any(in_rect(px,py,r) for r in rects) or any(in_rhom(px,py,d) for d in rhombi)

labels = [(float(x),float(y),t) for x,y,t in re.findall(
    r'<text x="([\d.]+)" y="([\d.]+)" font-size="10" fill="#6a6a6a"[^>]*>([^<]+)</text>', svg)]
bad = []
for lx,ly,t in labels:                                          # sample the label plate
    w = 7+len(t)*5.6
    pts = [(lx-w/2+w*i/6, ly-10+14*j/3) for i in range(7) for j in range(4)]
    if any(hits(px,py) for px,py in pts): bad.append(t)

paths = re.findall(r'<path d="M ([\d.]+) ([\d.]+)([^"]*)"', svg)
cross = []
for x1,y1,rest in paths:
    pts=[(float(x1),float(y1))]+[(float(a),float(b)) for a,b in re.findall(r"L ([\d.]+) ([\d.]+)", rest)]
    for (ax,ay),(bx,by) in zip(pts,pts[1:]):
        n = max(int(abs(bx-ax)+abs(by-ay))//4, 1)
        inner = [(ax+(bx-ax)*k/n, ay+(by-ay)*k/n) for k in range(1, n)]   # exclude endpoints
        if sum(1 for px,py in inner if hits(px,py)) > 2: cross.append((ax,ay,bx,by))

out = [r for r in rects if r[0]<0 or r[1]<0 or r[0]+r[2]>W or r[1]+r[3]>H]
cx=[r[0]+r[2]/2 for r in rects]
print(f"  labels on a shape        : {len(bad)}  {bad if bad else '-'}")
print(f"  edges through a shape    : {len(cross)}  {cross[:3] if cross else '-'}")
print(f"  shapes outside the canvas: {len(out)}  {out if out else '-'}")
print(f"  centring: mean box {sum(cx)/len(cx):.0f} vs canvas {W/2:.0f}")
sys.exit(1 if (bad or cross or out) else 0)
