"""Figure 3E: dot-product distributions of rank-1 target and decoy matches, Set 1 and Set 2, both
polarities, with the delivery threshold (500) marked."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

THRESHOLD = 500


def load():
    rows = list(csv.DictReader((HERE / "data" / "decoy_scores.csv").open()))
    for r in rows:
        r["dot"] = int(r["dot"])
    return rows


def draw(axes):
    """`axes`: four axes, one per (set, polarity)."""
    rows = load()
    bins = np.arange(0, 1001, 25)
    panels = [("Set 1", "Pos"), ("Set 1", "Neg"), ("Set 2", "Pos"), ("Set 2", "Neg")]
    for ax, (s, pol) in zip(axes, panels):
        t = [r["dot"] for r in rows if r["set"] == s and r["polarity"] == pol and r["kind"] == "target"]
        d = [r["dot"] for r in rows if r["set"] == s and r["polarity"] == pol and r["kind"] == "decoy"]
        ax.hist(t, bins=bins, color=style.COLOUR["target"], edgecolor=style.COLOUR["edge"], linewidth=0.2,
                label="target")
        ax.hist(d, bins=bins, color=style.COLOUR["decoy"], edgecolor=style.COLOUR["edge"], linewidth=0.2,
                label="decoy")
        ax.axvline(THRESHOLD, color=style.COLOUR["edge"], linewidth=0.6, linestyle="--")
        ax.set_yscale("log")
        ax.set_xlim(0, 1000)
        sign = "+" if pol == "Pos" else "−"
        ax.text(0.03, 0.95, f"{s} ({sign})\ntarget {len(t):,}\ndecoy {len(d):,}", transform=ax.transAxes,
                va="top", fontsize=6)
        if s == "Set 1" and pol == "Pos":
            ax.legend(loc="upper center", fontsize=6, ncol=2, bbox_to_anchor=(2.3, 1.28), handlelength=1.2)
        ax.set_xlabel("Dot product")
    axes[0].set_ylabel("Rank-1 matches")
    return axes


if __name__ == "__main__":
    style.apply()
    fig, axes = plt.subplots(1, 4, figsize=(style.DOUBLE_COLUMN, 1.9), sharey=True)
    fig.subplots_adjust(wspace=0.25)
    draw(list(axes)); style.letter(axes[0], "E")
    style.save(fig, HERE / "panels" / "panel_e")
