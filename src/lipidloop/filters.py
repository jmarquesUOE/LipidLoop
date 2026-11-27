"""Identification filters, values quoted from peak_finder/Utilities.java lines 12-17.

These are hard-coded in LipiDex and not exposed in its GUI. They are exposed here so a run
can record what it used, but the defaults reproduce LipiDex exactly.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class Thresholds:
    max_ppm_diff: float = 20.0        # MAXPPMDIFF
    min_dot_product: float = 500.0    # MINDOTPRODUCT
    min_rev_dot_product: float = 700.0  # MINREVDOTPRODUCT
    min_fa_purity: float = 75.0       # MINFAPURITY
    min_rt_multiplier: float = 0.5    # MINRTMULTIPLIER
    min_id_num: int = 1               # MINIDNUM

    def passes(self, ppm: float, dot: float, rev_dot: float, purity: float) -> bool:
        """The reverse threshold is deliberately stricter than the forward one: a real
        identification must explain the library's predicted fragments even when the spectrum
        carries extra peaks from co-isolation."""
        return (abs(ppm) <= self.max_ppm_diff
                and dot >= self.min_dot_product
                and rev_dot >= self.min_rev_dot_product
                and purity >= self.min_fa_purity)


MASS_OF_ELECTRON = 0.00054858026  # Utilities.java line 32
