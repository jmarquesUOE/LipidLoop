"""Read the mobile-phase modifier off the data instead of taking it on trust.

`--modifier` decides which library is searched, and it is currently a hand-entered field sourced
from whatever the deposit happened to record. Three of seventeen curated studies say UNKNOWN, and
because `formate` is also the default, an unevidenced guess is indistinguishable from an evidenced
fact once it reaches the results.

⚠ WHY THIS MATTERS MORE THAN A TIDY METADATA FIELD. The formate and acetate adducts differ by
EXACTLY one CH2:

    [M+HCOO]-    +44.9982
    [M+CH3COO]-  +59.0139
    difference    14.0157   =  CH2 to within 0.000 mDa

So `[PC 34:1 + HCOO]-` and `[PC 33:1 + CH3COO]-` have identical exact mass. This is a degeneracy,
not a near-miss: no mass accuracy separates them, at any resolution. Search the wrong modifier and
every choline lipid in negative mode is silently renamed to its odd-chain neighbour, at full dot
product, with a mass error of zero. Nothing downstream can catch it — which is why the check has
to happen at the front.

HOW IT WORKS. A mobile phase leaves its own fingerprint in the background. Sodium leached from
glassware makes sodium formate clusters on a formic-acid method and sodium acetate clusters on an
acetate one, and `exclusion.py` already models both repeat units because they dominate negative-mode
duty cycle. This module reuses that machinery and simply asks which of the two is present.

⚠ A tie is NOT evidence of absence. On thin sampling a known-formate study returned zero members of
either series, because its base peaks were all real lipids. `Detection.call` is None in that case
and the caller must treat it as "no evidence", never as a contradiction of the declared modifier.

⚠ One study need not have one modifier. MSV000094718 is a comparison OF ammonium salts — formate,
acetate, bicarbonate and fluoride in different sub-experiments. A per-study answer is not
well-defined there, and a mixed result is a true description rather than a failure.

Off by default. It costs a pass over the MS1 scans, and most runs already know their modifier from
the deposit; it exists for data that arrives without provenance.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from .exclusion import SERIES, SERIES_TOLERANCE, assign_series

# The two series that are diagnostic of the mobile phase. Everything else in SERIES is a
# contaminant that says nothing about which salt was used — siloxane comes off the column
# whatever you pump through it.
DIAGNOSTIC = {"sodium formate": "formate", "sodium acetate": "acetate"}

# Sampling bounds. A background cluster series is present across the whole run by definition, so
# it does not need every scan to be found — and reading every MS1 peak of a 224-file study is
# exactly the kind of thing that put 4.2 GB into swap.
MAX_SPECTRA_PER_FILE = 200
TOP_PEAKS_PER_SPECTRUM = 40
MIN_MEMBERS = 3                  # assign_series' own floor: two ions apart is a coincidence
MIN_MARGIN = 2.0                 # the winner must have this many times the loser's members
# ⚠ An absolute floor as well as a ratio, because a ratio cannot see thinness: 6 members against 0
# is an infinite margin and almost no evidence. Measured on real studies, a present series gives
# 177-4479 members; MSV000095868 gave 6 across three files and would otherwise have been called
# with the same confidence as a study carrying 700.
MIN_EVIDENCE = 30


@dataclass(slots=True)
class Detection:
    """What the background says the mobile phase was."""

    call: "str | None" = None            # "formate" | "acetate" | None when there is no evidence
    counts: dict = field(default_factory=dict)
    files: int = 0
    reason: str = ""

    @property
    def confident(self) -> bool:
        return self.call is not None

    def __str__(self) -> str:
        if not self.counts:
            return f"no salt-cluster series found across {self.files} file(s) — no evidence either way"
        seen = ", ".join(f"{k} {v}" for k, v in sorted(self.counts.items(), key=lambda kv: -kv[1]))
        head = self.call or "inconclusive"
        return f"{head} ({seen}; {self.files} file(s)){'; ' + self.reason if self.reason else ''}"

    def disagrees_with(self, declared: str) -> bool:
        """True only when there is positive evidence AGAINST the declared modifier.

        Absence of evidence never disagrees — see the tie note in the module docstring.
        """
        return self.confident and declared in DIAGNOSTIC.values() and self.call != declared


def survey_peaks(path, max_spectra: int = MAX_SPECTRA_PER_FILE,
                 top_peaks: int = TOP_PEAKS_PER_SPECTRUM) -> list[float]:
    """The most intense MS1 peaks from an evenly spaced sample of one file's scans.

    Evenly spaced rather than the first N: a gradient's first scans are all solvent front, and a
    series that is present across the whole run should be sampled across the whole run.
    """
    from pyopenms import MSExperiment, MzMLFile

    experiment = MSExperiment()
    MzMLFile().load(str(path), experiment)
    ms1 = [s for s in experiment if s.getMSLevel() == 1]
    if not ms1:
        return []
    step = max(1, len(ms1) // max_spectra)
    out: list[float] = []
    for spectrum in ms1[::step][:max_spectra]:
        mzs, intensities = spectrum.get_peaks()
        if len(mzs) == 0:
            continue
        order = sorted(range(len(mzs)), key=lambda i: intensities[i], reverse=True)[:top_peaks]
        out.extend(float(mzs[i]) for i in order)
    return out


def detect(peaks_by_file: "list[list[float]]", tolerance: float = SERIES_TOLERANCE) -> Detection:
    """Which diagnostic salt series is present in these peak lists.

    Series membership is decided per file and then pooled. Pooling the raw m/z across files first
    would let three ions from three different runs form a "chain" that never existed in any of
    them — the chain has to be real somewhere before it counts.
    """
    counts: dict[str, int] = {}
    files = 0
    for peaks in peaks_by_file:
        if not peaks:
            continue
        files += 1
        # discover=False: unnamed repeats are the exclusion list's business, not this question's,
        # and letting discovery run here only adds noise to a two-way decision.
        labels = assign_series(sorted(set(round(mz, 4) for mz in peaks)),
                               tolerance=tolerance, min_members=MIN_MEMBERS, discover=False)
        for name in labels.values():
            if name in DIAGNOSTIC:
                counts[name] = counts.get(name, 0) + 1

    if not counts:
        return Detection(None, {}, files, "")
    ranked = sorted(counts.items(), key=lambda kv: -kv[1])
    top, top_n = ranked[0]
    runner_n = ranked[1][1] if len(ranked) > 1 else 0
    if top_n < MIN_EVIDENCE:
        return Detection(None, counts, files,
                         f"only {top_n} members of the winning series — too thin to call "
                         f"(a present series normally gives hundreds)")
    if runner_n and top_n < MIN_MARGIN * runner_n:
        return Detection(None, counts, files,
                         f"both series present at comparable levels (x{top_n / runner_n:.1f}) — "
                         f"the run may genuinely mix modifiers")
    return Detection(DIAGNOSTIC[top], counts, files, "")


def detect_files(paths) -> Detection:
    """`detect` over a list of mzML paths, reading each one once."""
    return detect([survey_peaks(p) for p in paths])
