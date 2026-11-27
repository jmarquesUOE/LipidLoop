"""LOESS smoothing, as `Sample.java` uses it for retention-time alignment.

LipiDex fits Apache Commons Math's `LoessInterpolator(0.1, 1)` to (retention, deviation) pairs
and reads the smoothed deviation back to map an MS2's measured retention time into the aligned
frame. Reproduced here rather than pulled in from statsmodels or scipy: it is forty lines,
and the bandwidth and robustness-iteration count are part of the algorithm being reproduced,
not tuning knobs.

Local linear regression with tricube weights over the nearest `bandwidth * n` points, then
`robustness_iters` re-weightings by bisquare of the residuals.
"""
from __future__ import annotations

import numpy as np

DEFAULT_BANDWIDTH = 0.1     # LoessInterpolator(0.1, 1)
DEFAULT_ROBUSTNESS = 1


def loess(x, y, bandwidth: float = DEFAULT_BANDWIDTH,
          robustness_iters: int = DEFAULT_ROBUSTNESS) -> np.ndarray:
    """Smoothed y at each x. x must be sorted ascending and strictly increasing."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    n = x.size
    if n == 0:
        return np.array([])
    if n <= 2:
        return y.copy()

    bandwidth_points = max(2, int(bandwidth * n))
    robustness = np.ones(n)
    fitted = np.zeros(n)

    for iteration in range(robustness_iters + 1):
        for i in range(n):
            # Window of the nearest bandwidth_points neighbours, walked outward from i.
            lo = max(0, i - bandwidth_points // 2)
            hi = min(n, lo + bandwidth_points)
            lo = max(0, hi - bandwidth_points)

            xs = x[lo:hi]
            ys = y[lo:hi]
            radius = max(abs(x[i] - xs[0]), abs(xs[-1] - x[i]))
            if radius <= 0:
                fitted[i] = y[i]
                continue

            weights = (1.0 - np.abs((xs - x[i]) / radius) ** 3) ** 3
            weights = np.clip(weights, 0.0, None) * robustness[lo:hi]

            total = weights.sum()
            if total <= 0:
                fitted[i] = y[i]
                continue

            mean_x = (weights * xs).sum() / total
            mean_y = (weights * ys).sum() / total
            dx = xs - mean_x
            variance = (weights * dx * dx).sum()
            if variance <= 1e-12:
                fitted[i] = mean_y
                continue
            slope = (weights * dx * (ys - mean_y)).sum() / variance
            fitted[i] = mean_y + slope * (x[i] - mean_x)

        if iteration < robustness_iters:
            residuals = np.abs(y - fitted)
            median = np.median(residuals)
            if median <= 0:
                break
            scaled = np.clip(residuals / (6.0 * median), 0.0, 1.0)
            robustness = (1.0 - scaled ** 2) ** 2

    return fitted
