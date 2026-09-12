"""Figure 4D: PE on Set 2 (first arm), the same consensus species measured in both polarities:
log10 median area against log10 consensus concentration, positive and negative mode, with a
least-squares line per polarity and the within-polarity Spearman rho and slope printed (on log-log
axes a proportional response has slope 1). The polarity rule keeps PE from negative mode."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from scipy.stats import spearmanr

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "benchmark_pe.csv").open()))
    for r in rows:
        r["x"] = math.log10(float(r["value_nmol_mL"])); r["y"] = math.log10(float(r["area"]))
    return rows


def draw(ax):
    rows = load()
    import numpy as np
    for pol, key, marker in (("Pos", "positive", "o"), ("Neg", "negative", "s")):
        sub = [r for r in rows if r["polarity"] == pol]
        x = np.array([r["x"] for r in sub]); y = np.array([r["y"] for r in sub])
        rho = spearmanr(x, y).correlation
        # least-squares line on the log-log axes: a proportional response has slope 1
        slope, intercept = np.polyfit(x, y, 1)
        xs = np.array([x.min() - 0.05, x.max() + 0.05])
        ax.plot(xs, slope * xs + intercept, color=style.COLOUR[key], linewidth=1.0, zorder=1)
        ax.scatter(x, y, s=16, marker=marker, color=style.COLOUR[key], edgecolor=style.COLOUR["edge"], linewidth=0.3,
                   zorder=2, label=f"{key} mode, ρ = {rho:.2f}, slope = {slope:.2f} (n = {len(sub)})")
    ax.set_xlabel("Consensus concentration, log10 nmol/mL")
    ax.set_ylabel("Median area, log10")
    ax.legend(loc="upper center", fontsize=6, ncol=1, handlelength=1.0, bbox_to_anchor=(0.5, -0.22))
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.6))
    draw(ax); style.letter(ax, "D")
    style.save(fig, HERE / "panels" / "panel_d")
