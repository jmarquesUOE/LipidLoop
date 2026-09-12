"""Figure 4: validation beyond Set 1. A (corpus agreement per study), B (two-pass calibration),
C (Set 2 against the NIST consensus, per class), D (PE against the consensus, both polarities).

    python manuscript/figures/fig4/merge.py      -> fig4/fig4.svg, fig4/fig4.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import style     # noqa: E402
import panel_a   # noqa: E402
import panel_b   # noqa: E402
import panel_c   # noqa: E402
import panel_d   # noqa: E402


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(style.DOUBLE_COLUMN, 7.4))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.25, 1], height_ratios=[1, 1], wspace=0.75, hspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    panel_a.draw(ax_a); style.letter(ax_a, "A")
    panel_b.draw(ax_b); style.letter(ax_b, "B")
    panel_c.draw(ax_c); style.letter(ax_c, "C")
    panel_d.draw(ax_d); style.letter(ax_d, "D")
    for p in style.save(fig, HERE / "fig4"):
        print("written", p)


if __name__ == "__main__":
    main()
