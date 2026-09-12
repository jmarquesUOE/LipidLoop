"""Figure 2: exactness and recovery. Composes panels A (per-file exactness), B (recovery) and
C (Set 2 score pairs) into one double-column figure; letters only, legend in DRAFT.md.

    python manuscript/figures/fig2/merge.py      -> fig2/fig2.svg, fig2/fig2.pdf
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import style          # noqa: E402
import panel_a        # noqa: E402
import panel_b        # noqa: E402
import panel_c        # noqa: E402


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(style.DOUBLE_COLUMN, 5.2))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1.25, 1], height_ratios=[1, 1.05],
                  wspace=0.45, hspace=0.95)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    panel_a.draw(ax_a); style.letter(ax_a, "A")
    panel_b.draw(ax_b); style.letter(ax_b, "B")
    panel_c.draw(ax_c); style.letter(ax_c, "C")
    for p in style.save(fig, HERE / "fig2"):
        print("written", p)


if __name__ == "__main__":
    main()
