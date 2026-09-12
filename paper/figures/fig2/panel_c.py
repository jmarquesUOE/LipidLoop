"""Figure 2C: LipiDex dot product against the reimplementation's, every Set 2 spectrum matched
between the two (positive and negative mode). All points sit on the diagonal except the internal
standard the pipeline excludes from its library."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "set2_dot_pairs.csv").open()))
    for r in rows:
        r["lipidex_dot"] = int(r["lipidex_dot"])
        r["ours_dot"] = int(r["ours_dot"])
    return rows


def draw(ax):
    rows = load()
    ax.plot([0, 1000], [0, 1000], color=style.COLOUR["edge"], linewidth=0.5, zorder=1)
    for pol, label in (("Pos", "positive"), ("Neg", "negative")):
        sub = [r for r in rows if r["polarity"] == pol]
        ax.scatter([r["lipidex_dot"] for r in sub], [r["ours_dot"] for r in sub], s=9,
                   color=style.COLOUR[label], edgecolors=style.COLOUR["edge"], linewidths=0.2,
                   alpha=0.8, label=f"{label} (n = {len(sub):,})", zorder=2)
    off = [r for r in rows if r["lipidex_dot"] != r["ours_dot"]]
    for r in off:
        ax.annotate(r["name"].split(" ")[0], (r["lipidex_dot"], r["ours_dot"]), fontsize=6,
                    xytext=(6, -8), textcoords="offset points")
    ax.set_xlabel("Dot product, LipiDex")
    ax.set_ylabel("Dot product, reimplementation")
    ax.set_xlim(0, 1020)
    ax.set_ylim(0, 1020)
    ax.set_aspect("equal")
    ax.legend(loc="upper left", fontsize=6.5)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 3.0))
    draw(ax)
    style.letter(ax, "C")
    style.save(fig, HERE / "panels" / "panel_c")
