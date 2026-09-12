"""Figure 5: MS2 duty cycle and feedback to the instrument. A (per-study waste, left, full height),
B (cross-vendor contaminants, top right), C (Set 2 without/with lists, bottom right, three axes).

    python manuscript/figures/fig5/merge.py      -> fig5/fig5.svg, fig5/fig5.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import style     # noqa: E402
import panel_a   # noqa: E402
import panel_b   # noqa: E402
import panel_c   # noqa: E402


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(style.DOUBLE_COLUMN, 7.2))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.15, 1], height_ratios=[2.6, 1],
                  wspace=0.6, hspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    sub = GridSpecFromSubplotSpec(1, 3, subplot_spec=gs[1, :], wspace=0.5)
    axes_c = [fig.add_subplot(sub[0, i]) for i in range(3)]
    panel_a.draw(ax_a); style.letter(ax_a, "A")
    panel_b.draw(ax_b); style.letter(ax_b, "B")
    panel_c.draw(axes_c); style.letter(axes_c[0], "C")
    for p in style.save(fig, HERE / "fig5"):
        print("written", p)


if __name__ == "__main__":
    main()
