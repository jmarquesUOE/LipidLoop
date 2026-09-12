"""Find MS1 features a DDA run detected but never gave an MS2 to.

The complementary problem to `exclusion.py`. Exclusion frees duty cycle; this decides where some
of it should go on the next injection of the same study. A DDA run only fragments what looks
intense enough, ranked at the moment of acquisition — a real lipid can lose that competition by
co-eluting with something brighter, or simply by being dimmer in this exact sample set. Losing
Top-N is not evidence against the ion, only evidence the instrument never got to ask what it was.
So a candidate here is judged on abundance and library plausibility, not on any property of the
compound, and the list is a recommendation for the *next* injection of more of the same study —
nothing here identifies anything in the current run.

The list is a recommendation to put into the acquisition method. Nothing here changes processing,
and the script writes a file rather than editing anything.
"""
from __future__ import annotations

import bisect
import statistics
from dataclasses import dataclass


@dataclass(slots=True)
class Included:
    """One MS1 feature that lost Top-N and matches something the study's own library expects."""

    mz: float
    retention: float
    area: float
    library_names: list[str]    # matching library entries within `ms1_tol`, closest first
    n_matches: int

    def __str__(self) -> str:
        names = "; ".join(self.library_names[:2])
        return f"m/z {self.mz:.4f} @ {self.retention:.2f} min, area {self.area:,.0f}: {names}"


def find_inclusion_candidates(groups, all_scans, index, ms1_tol: float = 0.01,
                              rt_window: "float | None" = None,
                              max_matches: int = 3) -> list["Included"]:
    """The re-injection targets: MS1 features this run's own Top-N never gave an MS2 to.

    Three conditions, all required:

    1. **Low abundance relative to this run's own confirmed lipids** — `max_area` below the
       median of features THIS batch identified by MS2. An absolute intensity cutoff does not
       travel between instruments or methods; the batch's own confirmed population is the one
       abundance scale that is always available and already calibrated to the run, the same
       reasoning `calibrate.py` uses for the noise floor.
    2. **Never fragmented at all**, anywhere in the batch — no MS2 scan landed within `ms1_tol`
       of the feature's mass and within `rt_window` of where it eluted. This is the mirror image
       of `find_wasted`: that function starts from spectra that WERE acquired and never
       identified; this one starts from features that were detected on MS1 and never got a
       spectrum in the first place. `CompoundGroup.has_ms2` looks like the obvious flag for this
       and is deliberately not used — feature detection sets it `False` unconditionally at
       construction (features.py) and nothing in this codebase ever sets it `True`, so it carries
       no real MS2 history; trusting it would flag every feature in every run as a candidate.
    3. **Matches a library entry within `ms1_tol`.** Without this, a low-abundance, unfragmented
       feature is exactly as likely to be noise as a real lipid. The match is what turns "some
       ion was there" into "something this study's own search would have named, if only it had
       been asked."

    `groups` needs `final_lipid_id`, `quant_ion`, `max_area`, `retention` and `keep` — pass
    `PeakFinder.run()`'s compound groups, not the raw feature-detection output, so that
    identification status is known. Groups the presence/adduct/blank filters already rejected
    (`keep=False`) are skipped: a filtered-out feature (an in-source fragment, an adduct pair, a
    blank contaminant) is not a candidate worth spending an instrument's duty cycle re-confirming.

    `rt_window` defaults to six times this batch's mean chromatographic FWHM — the same "close
    enough to be the same peak" multiple `peakfinder.py` already uses for feature matching.
    """
    identified_areas = [g.max_area for g in groups
                        if g.final_lipid_id is not None and g.max_area > 0]
    if not identified_areas:
        return []
    abundance_ceiling = statistics.median(identified_areas)

    if rt_window is None:
        widths = [g.avg_fwhm for g in groups if getattr(g, "avg_fwhm", 0.0) > 0]
        rt_window = statistics.mean(widths) * 6.0 if widths else 0.25

    scans = sorted(all_scans, key=lambda s: s[1])
    scan_mz = [s[1] for s in scans]
    scan_rt = [s[2] for s in scans]

    found = []
    for group in groups:
        if not group.keep or group.final_lipid_id is not None:
            continue
        if group.quant_ion is None or not (0 < group.max_area < abundance_ceiling):
            continue
        if _ever_scanned(scan_mz, scan_rt, group.quant_ion, group.retention, ms1_tol, rt_window):
            continue
        matches = index.candidates(group.quant_ion)
        if not matches:
            continue
        names = [index.spectra[i].name for i in matches[:max_matches]]
        found.append(Included(mz=group.quant_ion, retention=group.retention,
                              area=group.max_area, library_names=names, n_matches=len(matches)))

    found.sort(key=lambda c: c.area)
    return found


def _ever_scanned(scan_mz: list[float], scan_rt: list[float], mz: float, retention: float,
                  ms1_tol: float, rt_window: float) -> bool:
    """Whether any MS2 scan in the batch landed on this feature's mass and elution window."""
    lo = bisect.bisect_left(scan_mz, mz - ms1_tol)
    hi = bisect.bisect_right(scan_mz, mz + ms1_tol)
    return any(abs(scan_rt[i] - retention) <= rt_window for i in range(lo, hi))


def write_inclusion_list(candidates: list[Included], path) -> None:
    """`m/z,Retention (min)` — the re-injection mass list, nothing else.

    Deliberately just the two numbers an acquisition method needs. The evidence behind each row
    (area, library match, why it qualified) is `Included` itself for anything that wants it
    programmatically; the file is the mass list, not the review document.
    """
    from pathlib import Path

    lines = ["m/z,Retention (min)"]
    for c in candidates:
        lines.append(f"{c.mz:.4f},{c.retention:.2f}")
    Path(path).write_text("\n".join(lines) + "\n")
