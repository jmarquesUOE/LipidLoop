"""Figure S2: LipiDex's dot product against LipidLoop's for every matched Set 2 spectrum.

    python manuscript/figures/figS2/merge.py      -> figS2/figS2.svg, figS2/figS2.pdf, figS2/figS2.png
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import style    # noqa: E402
import panel_a  # noqa: E402


def main() -> None:
    style.apply()
    fig, ax = plt.subplots(figsize=(style.SINGLE_COLUMN, 3.0))
    panel_a.draw(ax)
    for p in style.save(fig, HERE / "figS2"):
        print("written", p)
    import cairosvg
    cairosvg.svg2png(url=str(HERE / "figS2.svg"), write_to=str(HERE / "figS2.png"), output_width=1400)


if __name__ == "__main__":
    main()
