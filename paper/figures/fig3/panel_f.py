"""Figure 3F: one target library entry and its decoy, the worked example of the decoy construction.
TG 16:0_18:1_18:1 [M+NH4]+ from LipiDex_HCD_Formic: the precursor and the intensities are kept, every
chain-specific fragment is shifted by whole CH2 units (16:0 -> 19:0, 18:1 -> 21:1), so the decoy asserts
acyl chains that cannot coexist with the precursor mass it shares with the target."""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

# (m/z, intensity, label, label height offset) — the pairs sit 26 m/z apart, so labels are staggered
TARGET = [(577.5190, 999, "−18:1 acyl", 30), (603.5347, 999, "−16:0 acyl", 150), (239.2369, 100, "16:0 acylium", 30), (265.2526, 100, "18:1 acylium", 150)]
DECOY = [(619.5660, 999, "−21:1 acyl", 30), (645.5817, 999, "−19:0 acyl", 150), (281.2838, 100, "19:0 acylium", 30), (307.2995, 100, "21:1 acylium", 150)]
PRECURSOR = 876.802


def draw(ax):
    for mz, it, lab, dy in TARGET:
        ax.vlines(mz, 0, it, color=style.COLOUR["target"], linewidth=1.6)
    for mz, it, lab, dy in DECOY:
        ax.vlines(mz, 0, -it, color=style.COLOUR["decoy"], linewidth=1.6)
    ax.text(590, 1040, "acyl neutral losses\n16:0 and 18:1", ha="center", va="bottom", fontsize=5.5, color=style.COLOUR["edge"])
    ax.text(252, 150, "acylium ions\n16:0 and 18:1", ha="center", va="bottom", fontsize=5.5, color=style.COLOUR["edge"])
    ax.text(632, -1040, "acyl neutral losses\n19:0 and 21:1", ha="center", va="top", fontsize=5.5, color=style.COLOUR["edge"])
    ax.text(294, -150, "acylium ions\n19:0 and 21:1", ha="center", va="top", fontsize=5.5, color=style.COLOUR["edge"])
    ax.axhline(0, color=style.COLOUR["edge"], linewidth=0.6)
    ax.axvline(PRECURSOR, color=style.COLOUR["edge"], linewidth=0.6, linestyle="--")
    ax.text(PRECURSOR - 8, 1500, "shared\nprecursor", ha="right", va="top", fontsize=5.5)
    ax.text(0.02, 0.97, "target", transform=ax.transAxes, va="top", fontsize=6.5, fontweight="bold", color=style.COLOUR["edge"])
    ax.text(0.02, 0.03, "decoy", transform=ax.transAxes, va="bottom", fontsize=6.5, fontweight="bold", color=style.COLOUR["edge"])
    ax.set_xlim(200, 900)
    ax.set_ylim(-1550, 1550)
    ax.set_yticks([-999, -500, 0, 500, 999]); ax.set_yticklabels(["999", "500", "0", "500", "999"])
    ax.set_xlabel("m/z")
    ax.set_ylabel("Library intensity")
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.2))
    draw(ax); style.letter(ax, "F")
    style.save(fig, HERE / "panels" / "panel_f")
