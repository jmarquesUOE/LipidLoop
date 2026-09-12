# Gap filling

A zero in the table means the feature finder did not detect a peak in that injection. It does not
mean nothing is there. Statistical imputation — half the feature minimum, kNN, QRILC — answers the
question by assuming the answer. Re-integration asks the instrument.

## The parameters are Compound Discoverer's, deliberately

Read out of the facility's own workflow,
`CD_Workflows/20230531_LipiDex_ALIGNED_Pos.cdProcessingWF`, the 25 min Kinetex method these
datasets were acquired with:

    Fill Gaps (GapFillingNode)
        Mass Tolerance            10 ppm
        S/N Threshold             1.5
        Use Real Peak Detection   True

The 2025 Acquity 12 min workflow uses 5 ppm / S/N 3 for the same node. A replacement for a pipeline
should not quietly disagree with it about what a filled value means, so these are the defaults.

For the same reason `Mark Background Compounds` in those workflows sets **`Max. Sample/Blank = 3`**,
and the blank filter here now matches it. Running at 5 made this pipeline *stricter* than the
software it replaces, which is a difference nobody chose.

## It tests an assumption, which is the real reason it is here

`presence.py` deliberately keeps a lipid present in one group and absent from another — the
strongest result an experiment can produce. The tempting next step is to protect that whole-group
zero from ever being filled, on the grounds that filling it would turn presence/absence into a weak
fold change.

That reasoning is only sound **if the zeros are real**. Re-integration is how you find out. If a
supposedly absent group turns out to carry a peak above the S/N threshold, the absence was a
detection failure and the on/off call was wrong. Every run reports the count, because it is the
number that decides whether those calls stand.

## The scale problem, and why the filler can refuse

Areas in the table come from pyOpenMS `FeatureFindingMetabo` (`getIntensity()`). An area computed
here by a different summation is **not guaranteed to be in the same units**, and a table mixing two
scales is worse than a table with gaps — a gap at least announces itself.

So the integration is calibrated per injection against that injection's own detected features:
integrate the raw signal where the area is already known, take the median ratio of reported to
integrated, apply it. Two guards:

| guard | default | what it prevents |
|---|---|---|
| `min_calibration_features` | 50 | calibrating off a handful of features |
| `max_calibration_spread` | 0.5 | accepting a ratio whose IQR is more than half its median |

If either fails the injection is **refused** — its gaps stay gaps and the log says so. A filled
value on an unknown scale is a fabrication, and refusing is the only honest response.

## What gets filled

Only sample and QC injections. A zero in a blank is information about the blank, and filling it
would erase the evidence the blank filter runs on.

Each gap gets one of three outcomes, counted in the report:

| outcome | meaning |
|---|---|
| filled from a detected peak | a local maximum clearing S/N — a real measurement |
| filled from integrated signal | signal present but below S/N — an upper bound, still a number |
| left empty | nothing there; the zero stands, and now it means something |

## Configuration

    "gap_filling": {
      "enabled": true,
      "mass_tolerance_ppm": 10.0,
      "sn_threshold": 1.5,
      "rt_window": 0.2,
      "use_real_peak_detection": true,
      "max_calibration_spread": 0.5,
      "min_calibration_features": 50
    }

`rt_window` is minutes either side of the group's retention time. It has to absorb the alignment
shift — about 5.5 s spread on these methods, and CD's own aligner allows 0.5 min — without reaching
into a neighbouring peak.

Output goes to **`Final_Results_Filtered.csv`**, the analysis-ready table. `Final_Results.csv`
keeps its zeros so the two can be diffed.

## Is it right? Held-out validation

Calibrate on half of an injection's detected features, predict the other half, compare to the
areas already known:

    held-out predicted / actual    median 1.00    IQR 0.93-1.13    93.5% within 2x

and, the part that matters, **no intensity dependence** — across five quintiles from a median area
of 173,000 to 48,000,000 the ratio stays between 0.99 and 1.04. If the integral were inflated by
co-eluting signal in the 10 ppm window, faint features would be the worst case. They are not.

## What the gaps actually are

Gap rate against a row's own median intensity, per quintile:

    faintest  33.2%  ->  24.3%  ->  17.4%  ->  14.1%  ->  5.8%  brightest

So a gap is mostly a detection threshold being crossed, not a random failure — though even the
brightest quintile loses 5.8%. It is **not** a retention-time problem: the filled peak sits a
median 1.5 s from the group's retention time, 85.6% within one peak width and never more than 12 s
away, against a linking tolerance of 30 s.

Ranked inside its own row, a filled value lands at the **21st percentile** of that row's observed
values; 69% fall below the row's median and 53% in its bottom quartile. That is where
below-detection signal belongs — low, overlapping the bottom of the observed distribution.

⚠ **A statistic that misled, recorded so it is not repeated.** The first check asked how often a
filled value exceeded the *smallest observed value in its row*, and 74% did, which looked like
gross overestimation. It is not: the minimum of ~32 draws is an extreme order statistic sitting in
the low tail by construction, so any typical value clears it. Comparing a typical value against an
extreme answers nothing. The percentile rank above is the honest form of the same question.

## What this means for whole-group absences

Every whole-group absence tested carried a peak — 74/74 positive, 45/45 negative. Given the
intensity gradient above, the fair reading is **"present but faint enough to fall below detection
in that group"**, not "detection failed at random". Either way the absence was not real, so those
rows are fold changes rather than presence/absence calls — but the distinction matters when writing
about them.
