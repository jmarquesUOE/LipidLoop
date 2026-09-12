"""Figure 5C: Set 2 without versus with the exclusion and inclusion lists generated from its own
first acquisition: wasted duty cycle, MS2 scans on an excluded mass, MS2-confirmed identifications,
per polarity, as paired bars."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

METRICS = [("wasted", "wasted MS2 (%)"), ("on_excluded_mass", "MS2 on an excluded mass (%)"),
           ("ms2_confirmed", "MS2-confirmed lipids")]


def load():
    rows = list(csv.DictReader((HERE / "data" / "paired.csv").open()))
    for r in rows:
        for k in ("wasted", "on_excluded_mass", "ms2_confirmed", "delivered", "molecules"):
            r[k] = float(r[k])
    return rows


def draw(axes):
    """`axes` is a list of three axes, one per metric."""
    rows = load()
    groups = [("Pos", "positive"), ("Neg", "negative")]
    arms = ["without lists", "with lists"]
    w = 0.36
    for ax, (key, label) in zip(axes, METRICS):
        for gi, (pol, pname) in enumerate(groups):
            for ai, arm in enumerate(arms):
                r = next(r for r in rows if r["polarity"] == pol and r["arm"] == arm)
                x = gi + (ai - 0.5) * w
                ax.bar(x, r[key], width=w * 0.95, color=style.COLOUR[arm],
                       edgecolor=style.COLOUR["edge"], linewidth=0.3, label=arm if gi == 0 else None)
                ax.text(x, r[key], f"{r[key]:.0f}", ha="center", va="bottom", fontsize=6)
        ax.set_xticks([0, 1]); ax.set_xticklabels([p for _, p in groups])
        ax.set_ylabel(label)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.15)
    axes[1].legend(loc="upper center", fontsize=6, ncol=2, bbox_to_anchor=(0.5, -0.26), handlelength=1.2)
    return axes


if __name__ == "__main__":
    style.apply()
    fig, axes = plt.subplots(1, 3, figsize=(style.DOUBLE_COLUMN, 2.2))
    fig.subplots_adjust(wspace=0.5)
    draw(list(axes)); style.letter(axes[0], "C")
    style.save(fig, HERE / "panels" / "panel_c")
