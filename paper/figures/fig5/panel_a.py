"""Figure 5A: MS2 duty cycle spent on never-identified precursors, per study and polarity,
stacked by contaminant series; the corpus plus Set 1 and both Set 2 arms."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "wasted.csv").open()))
    for r in rows:
        for k in style.SERIES_ORDER + ["wasted"]:
            r[k] = float(r[k])
        r["nist"] = r["nist"] == "1"
    # in-house first, then corpus by descending positive-mode waste
    inhouse = [r for r in rows if r["group"] == "in-house"]
    corpus = [r for r in rows if r["group"] == "corpus"]
    pos_w = {r["dataset"]: r["wasted"] for r in corpus if r["polarity"] == "Pos"}
    corpus.sort(key=lambda r: (-pos_w.get(r["dataset"], 0), r["dataset"], r["polarity"] != "Pos"))
    return inhouse + corpus


def draw(ax):
    rows = load()
    y = list(range(len(rows)))
    left = [0.0] * len(rows)
    for name in style.SERIES_ORDER:
        vals = [r[name] for r in rows]
        ax.barh(y, vals, left=left, color=style.SERIES[name], edgecolor=style.COLOUR["edge"],
                linewidth=0.25, height=0.8, label=name)
        left = [a + b for a, b in zip(left, vals)]
    for yi, r in zip(y, rows):
        ax.text(r["wasted"] + 0.8, yi, f"{r['wasted']:.0f}%", va="center", fontsize=5.5)
    ax.set_yticks(y)
    labels = []
    for r in rows:
        sign = "+" if r["polarity"] == "Pos" else "−"
        labels.append(f"{r['dataset']} ({sign})")
    ax.set_yticklabels(labels, fontsize=5.5)
    for tick, r in zip(ax.get_yticklabels(), rows):
        if r["nist"]:
            tick.set_color(style.COLOUR["nist"]); tick.set_fontweight("bold")
    n_in = sum(1 for r in rows if r["group"] == "in-house")
    ax.axhline(n_in - 0.5, color=style.COLOUR["edge"], linewidth=0.5, linestyle=":")
    ax.invert_yaxis()
    ax.set_xlabel("MS2 scans on never-identified precursors (%)")
    ax.set_xlim(0, 80)
    ax.legend(loc="upper center", fontsize=6, ncol=4, bbox_to_anchor=(0.45, -0.11), handlelength=1.2, columnspacing=1.0)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN + 0.6, 5.2))
    draw(ax); style.letter(ax, "A")
    style.save(fig, HERE / "panels" / "panel_a")
