"""Figure 5B: contaminants recurring across vendors, as a presence matrix (contaminant x vendor),
restricted to the four vendors with measurable contamination; cell colour = vendor, cell text
= how the contaminant was found."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402

VENDORS = ["Agilent", "Bruker", "Sciex", "Thermo"]
SHORT = {"cholesterol in-source fragment": "cholesterol fragment", "DEHP (phthalate)": "DEHP"}


def load():
    rows = list(csv.DictReader((HERE / "data" / "contaminants_by_vendor.csv").open()))
    for r in rows:
        r["vendors"] = set(r["vendors"].split(";"))
        r["wasted_scans"] = int(r["wasted_scans"])
    rows.sort(key=lambda r: (-len(r["vendors"]), -r["wasted_scans"]))
    return rows


def draw(ax):
    rows = load()
    for i, r in enumerate(rows):
        for j, v in enumerate(VENDORS):
            present = v in r["vendors"]
            ax.add_patch(Rectangle((j, i), 1, 1, facecolor=style.COLOUR[v] if present else "white",
                                   edgecolor=style.COLOUR["edge"], linewidth=0.4))
        ax.text(len(VENDORS) + 0.15, i + 0.5, f"{r['wasted_scans']:,}", va="center", fontsize=6)
    ax.set_xlim(0, len(VENDORS) + 1.9)
    ax.set_ylim(len(rows), 0)
    ax.set_xticks([j + 0.5 for j in range(len(VENDORS))])
    ax.set_xticklabels(VENDORS, fontsize=6.5, rotation=30, ha="right")
    ax.set_yticks([i + 0.5 for i in range(len(rows))])
    ax.set_yticklabels([f"{SHORT.get(r['contaminant'], r['contaminant'])} ({'formula' if r['found_by'] == 'formula match' else 'series'})"
                        for r in rows], fontsize=6)
    ax.text(len(VENDORS) + 0.15, -0.25, "MS2 scans", fontsize=6, va="bottom")
    ax.tick_params(length=0)
    for sp in ax.spines.values():
        sp.set_visible(False)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 2.4))
    draw(ax); style.letter(ax, "B")
    style.save(fig, HERE / "panels" / "panel_b")
