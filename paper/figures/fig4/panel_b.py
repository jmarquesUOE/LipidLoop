"""Figure 4B: the two-pass mass calibration on the three studies whose within-batch drift exceeded
the 2 ppm floor: identifications before and after, with decoy hits printed."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "calibration.csv").open()))
    for r in rows:
        for k in ("before", "after", "decoy_before", "decoy_after"):
            r[k] = int(r[k])
        r["drift_ppm"] = float(r["drift_ppm"])
    return rows


def draw(ax):
    rows = load()
    y = list(range(len(rows)))
    h = 0.38
    ax.barh([yi - h / 2 for yi in y], [r["before"] for r in rows], height=h, color=style.COLOUR["without lists"],
            edgecolor=style.COLOUR["edge"], linewidth=0.3, label="before calibration")
    ax.barh([yi + h / 2 for yi in y], [r["after"] for r in rows], height=h, color=style.COLOUR["target"],
            edgecolor=style.COLOUR["edge"], linewidth=0.3, label="after calibration")
    for yi, r in zip(y, rows):
        ax.text(r["before"] + 6, yi - h / 2, f"{r['before']}  ({r['decoy_before']} decoy)", va="center", fontsize=5.5)
        ax.text(r["after"] + 6, yi + h / 2, f"{r['after']}  ({r['decoy_after']} decoy)", va="center", fontsize=5.5)
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['instrument']} ({'+' if r['polarity'] == 'Pos' else '−'})\n{r['drift_ppm']} ppm drift"
                        for r in rows], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Identified lipids")
    ax.set_xlim(0, max(r["after"] for r in rows) * 1.6)
    ax.legend(loc="upper center", fontsize=6, ncol=2, bbox_to_anchor=(0.5, -0.16), handlelength=1.2)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.6))
    draw(ax); style.letter(ax, "B")
    style.save(fig, HERE / "panels" / "panel_b")
