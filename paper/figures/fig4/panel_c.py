"""Figure 4C: Set 2 (first arm) against the NIST SRM 1950 interlaboratory consensus, per class
family: consensus entries recovered with MS2, recovered by the retention model only, or missed;
the out-of-scope entries no spectral library can name drawn apart, hatched."""
from __future__ import annotations

import csv
import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import style  # noqa: E402


def load():
    rows = list(csv.DictReader((HERE / "data" / "benchmark_recall.csv").open()))
    for r in rows:
        for k in ("n", "ms2", "rt_model", "missed"):
            r[k] = int(r[k])
    return rows


def draw(ax):
    rows = load()
    scope = [r for r in rows if r["family"] != "out of scope"]
    oos = next(r for r in rows if r["family"] == "out of scope")
    y = list(range(len(scope)))
    h = 0.7
    ax.barh(y, [r["ms2"] for r in scope], height=h, color=style.COLOUR["recovered"], edgecolor=style.COLOUR["edge"],
            linewidth=0.3, label="recovered, MS2")
    ax.barh(y, [r["rt_model"] for r in scope], left=[r["ms2"] for r in scope], height=h, color=style.COLOUR["other polarity"],
            edgecolor=style.COLOUR["edge"], linewidth=0.3, label="recovered, retention model only")
    ax.barh(y, [r["missed"] for r in scope], left=[r["ms2"] + r["rt_model"] for r in scope], height=h,
            color=style.COLOUR["neutral"], edgecolor=style.COLOUR["edge"], linewidth=0.3, label="missed")
    yo = len(scope) + 0.6
    ax.barh([yo], [oos["n"]], height=h, color="white", edgecolor=style.COLOUR["edge"], linewidth=0.3, hatch="////",
            label="out of scope (no library)")
    for yi, r in zip(y, scope):
        rec = r["ms2"] + r["rt_model"]
        ax.text(r["n"] + 0.8, yi, f"{rec}/{r['n']}", va="center", fontsize=5.5)
    ax.text(oos["n"] + 0.8, yo, f"{oos['n']}", va="center", fontsize=5.5)
    ax.set_yticks(y + [yo])
    ax.set_yticklabels([r["family"] for r in scope] + ["eicosanoids, bile acids,\ncholesterol, S1P, free FA"], fontsize=6)
    ax.invert_yaxis()
    ax.set_xlabel("Consensus lipids (≥ 5 laboratories)")
    ax.set_xlim(0, max(r["n"] for r in rows) * 1.22)
    total_n = sum(r["n"] for r in scope); total_rec = sum(r["ms2"] + r["rt_model"] for r in scope); total_ms2 = sum(r["ms2"] for r in scope)
    leg = ax.legend(loc="upper center", fontsize=5.5, ncol=2, bbox_to_anchor=(0.45, -0.16), handlelength=1.2, columnspacing=0.8,
                    title=f"{total_ms2} recovered with MS2, {total_rec} in all, of {total_n} in scope", title_fontsize=6)
    return ax


if __name__ == "__main__":
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 3.0))
    draw(ax); style.letter(ax, "C")
    style.save(fig, HERE / "panels" / "panel_c")
