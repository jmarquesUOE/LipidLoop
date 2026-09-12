"""Figure S3: the decoy fraction among rank-1 spectrum matches as a function of the dot-product
threshold, with the reverse dot product held at its 700 threshold, on Set 1 and Set 2 (first arm)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def draw(ax):
    rows = list(csv.DictReader((HERE / "data" / "decoy_vs_threshold.csv").open()))
    for (s, pol, key, ls) in (("Set 1", "Pos", "Set 1", "-"), ("Set 1", "Neg", "Set 1", "--"), ("Set 2", "Pos", "Set 2", "-"), ("Set 2", "Neg", "Set 2", "--")):
        sub = [r for r in rows if r["set"] == s and r["polarity"] == pol and int(r["threshold"]) > 0]
        ax.plot([int(r["threshold"]) for r in sub], [float(r["decoy_fraction_pct"]) for r in sub], linestyle=ls, marker="o", markersize=3,
                color=style.COLOUR[key], label=f"{s}, {'positive' if pol == 'Pos' else 'negative'}")
    ax.axvline(500, color=style.COLOUR["edge"], linewidth=0.6, linestyle=":")
    ax.set_xlabel("Dot-product threshold (reverse dot ≥ 700)")
    ax.set_ylabel("Decoy fraction of rank-1 matches (%)")
    ax.set_ylim(0, 1.0)
    ax.legend(loc="upper center", fontsize=6, ncol=2, bbox_to_anchor=(0.5, -0.25))
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.6))
    draw(ax)
    for p in style.save(fig, HERE / "figS3"):
        print("written", p)
    import cairosvg
    cairosvg.svg2png(url=str(HERE / "figS3.svg"), write_to=str(HERE / "figS3.png"), output_width=1400)
