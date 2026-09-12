"""Figure 4: one configuration on five vendors. A (agreement per study, left) and B (two-pass
calibration, right); letters only, legend in DRAFT.md.

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


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(style.DOUBLE_COLUMN, 3.8))
    gs = GridSpec(1, 2, figure=fig, width_ratios=[1.25, 1], wspace=0.75)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    panel_a.draw(ax_a); style.letter(ax_a, "A")
    panel_b.draw(ax_b); style.letter(ax_b, "B")
    for p in style.save(fig, HERE / "fig4"):
        print("written", p)


if __name__ == "__main__":
    main()
