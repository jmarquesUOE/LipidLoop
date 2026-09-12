"""Figure 2A: per file, spectra LipiDex identified, split into identical identification, exact
score tie, and other, for Set 1 (23 files) and Set 2 (8 files)."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "exactness.csv").open()))
    for r in rows:
        for k in ("spectra", "identical", "tie", "dot_identical"):
            r[k] = int(r[k])
        r["other"] = r["spectra"] - r["identical"] - r["tie"]
    order = {"Set 1": 0, "Set 2": 1}
    rows.sort(key=lambda r: (order[r["set"]], r["polarity"] != "Pos", r["file"]))
    return rows


def draw(ax):
    rows = load()
    y = list(range(len(rows)))
    left = [0] * len(rows)
    for key in ("identical", "tie", "other"):
        vals = [r[key] for r in rows]
        ax.barh(y, vals, left=left, color=style.COLOUR[key], edgecolor=style.COLOUR["edge"],
                linewidth=0.3, height=0.8, label={"identical": "identical identification",
                                                   "tie": "exact score tie", "other": "different"}[key])
        left = [a + b for a, b in zip(left, vals)]
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['file']} ({'+' if r['polarity'] == 'Pos' else '−'})" for r in rows], fontsize=5.5)
    ax.invert_yaxis()
    ax.set_xlabel("MS2 spectra identified by LipiDex")
    # set separators and labels
    n1 = sum(1 for r in rows if r["set"] == "Set 1")
    ax.axhline(n1 - 0.5, color=style.COLOUR["edge"], linewidth=0.5, linestyle=":")
    xmax = max(r["spectra"] for r in rows)
    ax.text(xmax * 1.02, (n1 - 1) / 2, "Set 1", rotation=270, va="center", ha="left", fontsize=7)
    ax.text(xmax * 1.02, n1 + (len(rows) - n1 - 1) / 2, "Set 2", rotation=270, va="center", ha="left", fontsize=7)
    ax.legend(loc="lower right", fontsize=6.5)
    ax.set_xlim(0, xmax * 1.08)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 4.2))
    draw(ax)
    style.letter(ax, "A")
    style.save(fig, HERE / "panels" / "panel_a")
