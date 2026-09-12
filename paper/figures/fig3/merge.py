"""Figure 3: per-class retention models (A to D, 2 x 2), the decoy test (E, four axes in one row) and the
decoy construction worked example (F).

    python manuscript/figures/fig3/merge.py      -> fig3/fig3.svg, fig3/fig3.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import style             # noqa: E402
import panel_retention   # noqa: E402
import panel_e           # noqa: E402
import panel_f           # noqa: E402

CLASSES = [("Set 1", "TG", "Pos"), ("Set 1", "PC", "Pos"), ("Set 2", "TG", "Pos"), ("Set 2", "PC", "Pos")]


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(style.DOUBLE_COLUMN, 9.2))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[2.4, 1, 0.85], hspace=0.5)
    top = GridSpecFromSubplotSpec(2, 2, subplot_spec=gs[0], wspace=0.3, hspace=0.4)
    axes_r = []
    for i, ((set_, cls, pol), letter) in enumerate(zip(CLASSES, "ABCD")):
        ax = fig.add_subplot(top[i // 2, i % 2])
        panel_retention.draw(ax, set_, cls, pol); style.letter(ax, letter)
        axes_r.append(ax)
    panel_retention.shared_legend(fig, axes_r, panel_retention.MAX_DB)
    bottom = GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[1], wspace=0.25)
    axes_e = [fig.add_subplot(bottom[0, i]) for i in range(4)]
    for ax in axes_e[1:]:
        ax.sharey(axes_e[0]); ax.tick_params(labelleft=False)
    panel_e.draw(axes_e); style.letter(axes_e[0], "E")
    last = GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[2], width_ratios=[1.3, 1], wspace=0.3)
    ax_f = fig.add_subplot(last[0, 0])
    panel_f.draw(ax_f); style.letter(ax_f, "F")
    for p in style.save(fig, HERE / "fig3"):
        print("written", p)


if __name__ == "__main__":
    main()
