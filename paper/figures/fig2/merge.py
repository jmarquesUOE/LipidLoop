"""Figure 2: identification is correct where it can be checked. A (per-file agreement with the
reference workflow), B (recovery of its molecules), C (Set 2 against the NIST SRM 1950 consensus,
per class), D (PE against the consensus in both polarities); letters only, legend in DRAFT.md.

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
import panel_d        # noqa: E402


def main() -> None:
    style.apply()
    fig = plt.figure(figsize=(style.DOUBLE_COLUMN, 8.2))
    gs = GridSpec(3, 2, figure=fig, width_ratios=[1.05, 1], height_ratios=[0.75, 1.45, 1.0],
                  wspace=0.5, hspace=0.75)
    ax_a = fig.add_subplot(gs[:, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 1])
    ax_d = fig.add_subplot(gs[2, 1])
    panel_a.draw(ax_a); style.letter(ax_a, "A")
    panel_b.draw(ax_b); style.letter(ax_b, "B")
    panel_c.draw(ax_c); style.letter(ax_c, "C")
    panel_d.draw(ax_d); style.letter(ax_d, "D")
    for p in style.save(fig, HERE / "fig2"):
        print("written", p)


if __name__ == "__main__":
    main()
