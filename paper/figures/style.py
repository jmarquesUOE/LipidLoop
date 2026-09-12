"""One style for every manuscript figure.

House rules (Jair): no titles inside the image (legends live in DRAFT.md), panel letters only,
colour-blind-safe palette in pastel tones, colour carries the group on same-size marks, no star
markers, SVG + PDF output, one colour code used consistently across all figures.

The palette is Okabe-Ito blended 40% towards white so it reads as pastel while keeping the hue
separation the original was designed for.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

# Okabe-Ito, pastelised.
def _pastel(hex_colour: str, mix: float = 0.40) -> str:
    r, g, b = (int(hex_colour[i:i + 2], 16) for i in (1, 3, 5))
    r, g, b = (round(c + (255 - c) * mix) for c in (r, g, b))
    return f"#{r:02x}{g:02x}{b:02x}"

_OKABE = {
    "orange": "#E69F00", "sky": "#56B4E9", "green": "#009E73", "yellow": "#F0E442",
    "blue": "#0072B2", "vermillion": "#D55E00", "purple": "#CC79A7", "black": "#000000",
}
PASTEL = {k: _pastel(v) for k, v in _OKABE.items()}
PASTEL["grey"] = "#BFBFBF"
PASTEL["lightgrey"] = "#E3E3E3"

# The colour code. Same meaning in every figure.
COLOUR = {
    # polarity
    "positive": PASTEL["orange"],
    "negative": PASTEL["blue"],
    # in-house sets
    "Set 1": PASTEL["green"],
    "Set 2": PASTEL["purple"],
    # target vs decoy
    "target": PASTEL["sky"],
    "decoy": PASTEL["vermillion"],
    # exactness categories
    "identical": PASTEL["sky"],
    "tie": PASTEL["yellow"],
    "other": PASTEL["vermillion"],
    # recovery categories
    "recovered": PASTEL["sky"],
    "other polarity": PASTEL["yellow"],
    "absent": PASTEL["vermillion"],
    "standard": PASTEL["grey"],
    # acquisition arms
    "without lists": PASTEL["grey"],
    "with lists": PASTEL["green"],
    # vendors
    "Agilent": PASTEL["orange"],
    "Bruker": PASTEL["sky"],
    "Sciex": PASTEL["green"],
    "Thermo": PASTEL["purple"],
    "Waters": PASTEL["vermillion"],
    # NIST SRM 1950 emphasis
    "nist": "#4D4D4D",
    # neutral
    "neutral": PASTEL["grey"],
    "edge": "#4D4D4D",
}

# Contaminant series in the duty-cycle figures, same colour wherever a series appears.
SERIES = {
    "polysiloxane": PASTEL["purple"],
    "sodium formate": PASTEL["blue"],
    "sodium acetate": PASTEL["sky"],
    "PEG": PASTEL["orange"],
    "PPG": PASTEL["yellow"],
    "PTFE/PFPE": PASTEL["green"],
    "other named series": PASTEL["vermillion"],
    "unassigned": PASTEL["lightgrey"],
}
SERIES_ORDER = list(SERIES)

# Double-bond count in retention panels: a sequential, colour-blind-safe ramp (viridis) but
# pastelised the same way, so it sits with the rest.
def double_bond_colours(n: int) -> list[str]:
    cmap = plt.get_cmap("viridis")
    out = []
    for i in range(n):
        r, g, b, _ = cmap(i / max(n - 1, 1))
        out.append(_pastel("#%02x%02x%02x" % (round(r * 255), round(g * 255), round(b * 255)), 0.30))
    return out

SINGLE_COLUMN = 3.33   # inches, ACS
DOUBLE_COLUMN = 7.0

def apply() -> None:
    mpl.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 8,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "legend.fontsize": 7,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.linewidth": 0.6,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "lines.linewidth": 1.0,
        "lines.markersize": 4,
        "legend.frameon": False,
        "svg.fonttype": "none",       # keep text as text in SVG
        "pdf.fonttype": 42,           # embed TrueType, editable text in PDF
        "figure.dpi": 150,
    })

def letter(ax, text: str) -> None:
    """Panel letter at the top-left corner, outside the axes, bold, the only label a panel carries."""
    ax.text(-0.12, 1.04, text, transform=ax.transAxes, fontsize=10, fontweight="bold",
            va="bottom", ha="right")

def save(fig, stem: Path | str) -> list[Path]:
    """Write <stem>.svg and <stem>.pdf beside each other."""
    stem = Path(stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    out = []
    for ext in ("svg", "pdf"):
        p = stem.with_suffix(f".{ext}")
        fig.savefig(p, bbox_inches="tight")
        out.append(p)
    return out
