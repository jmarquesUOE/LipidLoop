"""Write `<analysis>/Method_Lists/` for an analysis that already exists.

`run_study.py` does this on every run; this is for analyses produced before it did, or for
re-checking with a different exclusion tolerance once the method's own value is known.

    python scripts/build_method_lists.py /path/to/Analysis_dir [--ppm 10] [--inclusion-window 0.5]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.method_lists import (DEFAULT_INCLUSION_WINDOW, DEFAULT_PPM,   # noqa: E402
                                     ISOLATION_HALF_WIDTH, write_method_lists)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("analysis", help="analysis directory holding Results/Pos and/or Results/Neg")
    ap.add_argument("--ppm", type=float, default=DEFAULT_PPM,
                    help="exclusion mass tolerance the acquisition method will use")
    ap.add_argument("--isolation", type=float, default=ISOLATION_HALF_WIDTH,
                    help="half the isolation width, Da")
    ap.add_argument("--inclusion-window", type=float, default=DEFAULT_INCLUSION_WINDOW,
                    help="retention window either side of each inclusion entry, min")
    args = ap.parse_args()
    summary = write_method_lists(Path(args.analysis), ppm=args.ppm, isolation=args.isolation,
                                 inclusion_window=args.inclusion_window)
    if not summary["polarities"]:
        print("nothing to do: no Results/<pol>/Exclusion_List.csv or Inclusion_List.csv found",
              file=sys.stderr)
        return 1
    print(f"written to {Path(args.analysis) / 'Method_Lists'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
