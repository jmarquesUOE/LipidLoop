"""The report must render, end to end, on the studies that have actually broken it.

Every QC failure so far passed the unit suite first. `median_pool_cv`, then `missing_identified`
— both bare `{}` early returns, both found by an overnight batch rather than by a test, and both
on the same kind of study: one with no usable samples. The unit tests exercised the functions in
isolation on ordinary inputs, so nothing ever assembled a degenerate run and asked the report to
draw it.

This is the missing layer. It is slow by the standards of the rest of the suite because it builds
figures and a PDF, and that is the point — the crashes lived between the statistics and the
rendering, which is exactly the seam no faster test crosses.

The pools-only case is not hypothetical: ST000991_DDA is five pool injections and nothing else.
"""
from __future__ import annotations

import csv
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[1]

SHAPES = {
    # ST000991_DDA exactly: every injection a pool, no samples to report on
    "pools_only": ["Pool_1", "Pool_2", "Pool_3", "Pool_4", "Pool_5"],
    # nothing to compare samples against, and no order to test drift along
    "no_pools": ["Sample_1", "Sample_2", "Sample_3"],
    # a single injection: no variance, no pairs, no correlation defined
    "one_injection": ["Sample_1"],
}


def _write(path: Path, columns: list[str]) -> None:
    rng = np.random.default_rng(0)
    areas = rng.lognormal(12, .4, size=(40, len(columns)))
    header = ["Retention Time (min)", "Quant Ion", "Polarity", "Area (max)", "Identification",
              "Lipid Class", "Features Found", "Dot Product", "Identification Source"] + columns
    with path.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for i, row in enumerate(areas):
            writer.writerow([10 + i * .1, 700 + i, "+", max(row),
                             "PC 16:0_18:1" if i % 3 else "", "PC",
                             len(columns), 900, "MS2"] + list(row))


@pytest.mark.slow  # registered in pytest.ini
@pytest.mark.parametrize("kind", sorted(SHAPES))
def test_the_report_renders_on_a_degenerate_study(tmp_path, kind):
    for name in ("Final_Results.csv", "Unfiltered_Results.csv"):
        _write(tmp_path / name, SHAPES[kind])

    out = tmp_path / "QC"
    result = subprocess.run(
        [sys.executable, str(REPO / "scripts/qc_report.py"),
         "--results-pos", str(tmp_path / "Final_Results.csv"),
         "--unfiltered-pos", str(tmp_path / "Unfiltered_Results.csv"),
         "--study", f"degenerate {kind}", "--owner", "test", "--out", str(out)],
        capture_output=True, text=True, timeout=900)

    # Report the traceback, not just a non-zero code. An exit number alone sent one of these
    # failures back to an overnight log to be diagnosed the slow way.
    assert result.returncode == 0, f"{kind} failed:\n{result.stderr[-3000:]}"
    assert (out / "QC_report.html").exists(), f"{kind} exited 0 but wrote no HTML"
