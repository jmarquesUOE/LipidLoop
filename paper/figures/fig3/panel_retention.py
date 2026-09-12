"""Figure 3A to D: retention time against total acyl carbons for one lipid class of Set 1, points
coloured by double-bond count, the pipeline's fitted model drawn per double-bond count, and the
features the model named without MS2 shown at the same size in a distinct colour.

`draw(ax, cls, pol)` draws one class; the merge script calls it four times."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

MAX_DB = 15   # one colour scale for every class panel


def load(set_: str, cls: str, pol: str):
    pts = [r for r in csv.DictReader((HERE / "data" / "retention_points.csv").open())
           if r["set"] == set_ and r["class"] == cls and r["polarity"] == pol and r["prefix"] == ""]
    for r in pts:
        r["carbons"] = int(r["carbons"]); r["double_bonds"] = int(r["double_bonds"])
        r["retention"] = float(r["retention"])
    model = next(r for r in csv.DictReader((HERE / "data" / "retention_models.csv").open())
                 if r["set"] == set_ and r["class"] == cls and r["polarity"] == pol)
    for k in ("intercept", "per_carbon", "per_double_bond", "r2", "residual_sd"):
        model[k] = float(model[k])
    return pts, model


def draw(ax, set_: str, cls: str, pol: str):
    pts, m = load(set_, cls, pol)
    dbs = sorted({r["double_bonds"] for r in pts})
    palette = style.double_bond_colours(MAX_DB + 1)
    colours = {db: palette[min(db, MAX_DB)] for db in dbs}
    ms2 = [r for r in pts if r["source"] == "MS2"]
    ext = [r for r in pts if r["source"] != "MS2"]
    c_lo, c_hi = min(r["carbons"] for r in pts), max(r["carbons"] for r in pts)
    for db in dbs:
        sub = [r for r in ms2 if r["double_bonds"] == db]
        if sub:
            ax.scatter([r["carbons"] for r in sub], [r["retention"] for r in sub], s=14,
                       color=colours[db], edgecolors=style.COLOUR["edge"], linewidths=0.3, zorder=3)
        cs = np.array([c for c in range(c_lo, c_hi + 1)
                       if any(r["double_bonds"] == db and abs(r["carbons"] - c) <= 4 for r in pts)])
        if len(cs) >= 2:
            ax.plot(cs, m["intercept"] + m["per_carbon"] * cs + m["per_double_bond"] * db,
                    color=colours[db], linewidth=0.8, zorder=2)
    if ext:
        ax.scatter([r["carbons"] for r in ext], [r["retention"] for r in ext], s=14,
                   facecolors="white", edgecolors=style.COLOUR["decoy"], linewidths=0.8, zorder=4)
    # Focus the axes on the model's own range; points the fit trimmed can sit far off it (a
    # void-volume PC, a late-eluting TG) and would otherwise compress every real point.
    preds = [m["intercept"] + m["per_carbon"] * r["carbons"] + m["per_double_bond"] * r["double_bonds"] for r in pts]
    lo, hi = min(preds) - 0.7, max(preds) + 1.2   # headroom for the two-line label at the top left
    outside = sum(1 for r in pts if not lo <= r["retention"] <= hi)
    ax.set_ylim(lo, hi)
    ax.set_xlabel("Total acyl carbons")
    ax.set_ylabel("Retention time (min)")
    sign = "+" if pol == "Pos" else "−"
    label = f"{set_}, {cls} ({sign})  n = {m['n_used']}/{m['n_total']}  R² = {m['r2']:.2f}"
    if outside:
        label += f"\n{outside} trimmed point{'s' if outside > 1 else ''} outside the axis"
    ax.text(0.03, 0.97, label, transform=ax.transAxes, va="top", ha="left", fontsize=6)
    return dbs, colours


def shared_legend(fig, axes, max_db: int):
    """One colour bar for double-bond count and one marker legend, shared by the four panels."""
    import matplotlib as mpl
    from matplotlib.colors import ListedColormap, BoundaryNorm
    cols = style.double_bond_colours(max_db + 1)
    cmap = ListedColormap(cols)
    norm = BoundaryNorm(np.arange(-0.5, max_db + 1.5, 1), cmap.N)
    sm = mpl.cm.ScalarMappable(cmap=cmap, norm=norm)
    cb = fig.colorbar(sm, ax=axes, orientation="vertical", fraction=0.025, pad=0.02, ticks=range(0, max_db + 1, 2))
    cb.set_label("double bonds", fontsize=7)
    cb.ax.tick_params(labelsize=6)
    handles = [Line2D([], [], marker="o", linestyle="", markersize=4, markerfacecolor=style.COLOUR["neutral"],
                      markeredgecolor=style.COLOUR["edge"], markeredgewidth=0.3, label="identified by MS2"),
               Line2D([], [], marker="o", linestyle="", markersize=4, markerfacecolor="white",
                      markeredgecolor=style.COLOUR["decoy"], markeredgewidth=0.8, label="named by the model, no MS2"),
               Line2D([], [], color=style.COLOUR["edge"], linewidth=0.8, label="fitted model per double-bond count")]
    axes[0].legend(handles=handles, loc="upper center", fontsize=6, ncol=3, bbox_to_anchor=(1.15, 1.32),
                   handletextpad=0.4, columnspacing=1.2)


if __name__ == "__main__":
    style.apply()
    for letter, (set_, cls, pol) in zip("ABCD", [("Set 1", "TG", "Pos"), ("Set 1", "PC", "Pos"), ("Set 2", "TG", "Pos"), ("Set 2", "PC", "Pos")]):
        fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.6))
        draw(ax, set_, cls, pol); style.letter(ax, letter)
        shared_legend(fig, [ax], 15)
        style.save(fig, HERE / "panels" / f"panel_{letter.lower()}")
        plt.close(fig)
