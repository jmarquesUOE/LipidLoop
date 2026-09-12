"""Figure 2B: recovery of the reference workflow's molecules, Set 1 and Set 2, both polarities,
as stacked fractions: recovered, found in the other polarity, absent, internal standard."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

KEYS = (("recovered", "recovered"), ("other_polarity", "other polarity"),
        ("absent", "absent"), ("standard", "standard"))
LABEL = {"recovered": "recovered", "other polarity": "reported in the other polarity",
         "absent": "not reported", "standard": "internal standard (excluded)"}


def load():
    rows = list(csv.DictReader((HERE / "data" / "recovery.csv").open()))
    for r in rows:
        for k in ("reference", "recovered", "other_polarity", "absent", "standard"):
            r[k] = int(r[k])
    rows.sort(key=lambda r: (r["set"], r["polarity"] != "Pos"))
    return rows


def draw(ax):
    rows = load()
    x = list(range(len(rows)))
    bottom = [0.0] * len(rows)
    for key, colour in KEYS:
        vals = [100 * r[key] / r["reference"] for r in rows]
        ax.bar(x, vals, bottom=bottom, color=style.COLOUR[colour], edgecolor=style.COLOUR["edge"],
               linewidth=0.3, width=0.7, label=LABEL[colour])
        bottom = [a + b for a, b in zip(bottom, vals)]
    for xi, r in zip(x, rows):
        ax.text(xi, 100 * r["recovered"] / r["reference"] - 4, f"{100 * r['recovered'] / r['reference']:.1f}%",
                ha="center", va="top", fontsize=6.5, color=style.COLOUR["edge"])
        ax.text(xi, 101.5, f"n = {r['reference']}", ha="center", va="bottom", fontsize=6)
    ax.set_xticks(x)
    ax.set_xticklabels([f"{r['set']}\n{'positive' if r['polarity'] == 'Pos' else 'negative'}" for r in rows])
    ax.set_ylabel("Molecules (%)")
    ax.set_ylim(0, 112)
    ax.legend(loc="upper center", fontsize=6, ncol=2, bbox_to_anchor=(0.5, -0.28), handlelength=1.2, columnspacing=1.0)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.8))
    draw(ax)
    style.letter(ax, "B")
    style.save(fig, HERE / "panels" / "panel_b")
