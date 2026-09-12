"""Figure 4A: agreement with each deposit's own published identifications, per study, coloured
by vendor; NIST SRM 1950 acquisitions in bold; decoy rate printed at the bar end; studies that
publish no list shown as labelled empty rows."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Patch

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "corpus.csv").open()))
    for r in rows:
        r["agreement"] = float(r["agreement"]) if r["agreement"] else None
        r["decoy_fdr"] = float(r["decoy_fdr"]); r["nist"] = r["nist"] == "1"
        r["injections"] = int(r["injections"])
    with_list = sorted([r for r in rows if r["agreement"] is not None], key=lambda r: -r["agreement"])
    without = sorted([r for r in rows if r["agreement"] is None], key=lambda r: r["dataset"])
    return with_list + without


def draw(ax):
    rows = load()
    y = list(range(len(rows)))
    for yi, r in zip(y, rows):
        if r["agreement"] is not None:
            ax.barh(yi, r["agreement"], color=style.COLOUR[r["vendor"]], edgecolor=style.COLOUR["edge"],
                    linewidth=0.3, height=0.75)
            ax.text(r["agreement"] + 1, yi, f"{r['agreement']:.0f}%   decoy {r['decoy_fdr']:.2f}%",
                    va="center", fontsize=6)
        else:
            ax.barh(yi, 100, color="white", edgecolor=style.COLOUR["edge"], linewidth=0.3, height=0.75,
                    linestyle=":")
            ax.text(1.5, yi, f"no published list   decoy {r['decoy_fdr']:.2f}%", va="center", fontsize=6,
                    color=style.COLOUR["edge"])
    ax.set_yticks(y)
    ax.set_yticklabels([f"{r['instrument']}, {r['injections']} inj." for r in rows], fontsize=6)
    for tick, r in zip(ax.get_yticklabels(), rows):
        if r["nist"]:
            tick.set_fontweight("bold")
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Agreement with the deposit's published identifications (%)")
    n_with = sum(1 for r in rows if r["agreement"] is not None)
    ax.axhline(n_with - 0.5, color=style.COLOUR["edge"], linewidth=0.5, linestyle=":")
    handles = [Patch(facecolor=style.COLOUR[v], edgecolor=style.COLOUR["edge"], linewidth=0.3, label=v)
               for v in ("Agilent", "Bruker", "Sciex", "Thermo", "Waters")]
    ax.legend(handles=handles, loc="upper center", fontsize=6, ncol=5, bbox_to_anchor=(0.5, -0.20), handlelength=1.2, columnspacing=1.0)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.DOUBLE_COLUMN * 0.6, 3.6))
    draw(ax); style.letter(ax, "A")
    style.save(fig, HERE / "panels" / "panel_a")
