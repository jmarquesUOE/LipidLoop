# Deriving parameters from the run

The proposal was: run once, measure the data, adjust the thresholds to it, run again. The boundary
that makes it safe:

> **Iterate on parameters that DESCRIBE the measurement.
> Never on parameters that DECIDE what counts as a hit.**

Peak width, mass calibration, retention spread and the noise floor are properties of the run and
can be measured. Dot product 500 / reverse 700, fatty-acid purity 75, the blank multiplier and the
minimum feature count are decision rules. Tuning those on the data makes the result unfalsifiable —
*"we chose the threshold that gave the most identifications"* is not a method — and destroys
comparability with LipiDex, which is what the whole reimplementation rests on.

**The boundary is enforced, not documented.** `calibrate.FORBIDDEN` is checked on every derivation
and raises. A convention in a docstring survives exactly as long as the person who wrote it, and
there is a test per forbidden name.

## What is derived, and from what

| parameter | derived from | why the configured value was wrong |
|---|---|---|
| `mass_offset_ppm` | median ppm delta over matched spectra | a symmetric ±10 ppm window centred on 0 is really −8.2/+11.8 at a median of −1.78 ppm |
| `chrom_fwhm` | median measured peak width | 5.0 s configured against 5.5 s observed here and 4.7 s on the CKD method; a 60 min gradient would differ far more |
| `min_fwhm` / `max_fwhm` | 1st and 99th percentiles of the width distribution | 1–20 s were guesses |
| `align_rt_max_difference` | 4 × the 95th percentile of observed retention deviation | 30 s against a ~5.5 s spread is not conservative, it is permissive — it lets the aligner pair features that are not the same peak |
| `retention_model_max_error` | 3 × the median fitted residual sd | the 1.0 min floor is arbitrary, and on a fast gradient it is most of the run |
| **`noise_threshold`** | 5.5 × the file's own measured noise floor | **an absolute intensity** — see below |

## The noise floor is measured on EVERY run, whether or not it is used

It is the one derivable parameter that needs no second pass — it comes from raw MS1 peaks, which
exist before feature detection. So it is always measured and always reported, and applied only when
`derive_noise_threshold` says so.

That split is deliberate, because **the failure mode is silence**. A threshold in absolute counts
fitted on one detector is not visibly wrong on the next: set too high, features are simply absent
and nobody knows what is missing; set too low, mass trace detection explodes and the run is blamed
on the file. Measuring costs seconds and turns a silent failure into a loud one:

    Pos: MS1 noise floor 1,009 counts; configured noise_threshold 5,000 is 5.0x it
    Neg: MS1 noise floor   941 counts; configured noise_threshold 5,000 is 5.3x it

Outside a plausible 2–15× band the log carries a warning naming the cause. The measured floor is
written into `run_config.json` on every run, so a result carries the evidence for whether its noise
threshold suited the instrument that produced it — including results produced before anyone thought
to ask.

## The noise floor is the one that matters

`noise_threshold` is an absolute intensity, so 5,000 on an Orbitrap Fusion Lumos is not 5,000 on a
Q Exactive Plus. For a facility pipeline running on several instruments this is the most likely
thing to be silently wrong on the next new one, and it fails quietly in both directions: too high
and features vanish, too low and mass trace detection explodes and the run never finishes.

**The multiplier is anchored by measurement, not chosen.** The floor on the Lumos brain files is
920 counts in positive and 866 in negative, against the `noise_threshold` of 5,000 that was fitted
against the reference results — so the multiplier that reproduces the validated setting on the
validated instrument is 5.44 and 5.77. It is set to **5.5**. A derivation that changes the answer
on the data it was tuned against is not a calibration, it is a different pipeline, and there is a
test pinning exactly that.

## What is under the threshold: measured, on the pooled QCs

The obvious objection to a detection floor is that nobody has looked below it. The pooled QCs are
the right material to look with, because they are the same sample eight times — anything recovered
can be tested for reproducibility instead of merely counted.

Full pipeline on the eight positive pooled QCs, at three thresholds:

| `noise_threshold` | runtime | features | identified rows | **distinct molecules** | median QC CV | in all 8 QCs |
|---|---|---|---|---|---|---|
| 5,000 | 128 s | 3,509 | 586 | **416** | 44.1% | 39.1% |
| 2,500 | 160 s | 4,348 | 591 | **416** | 47.3% | — |
| 1,000 | 218 s | 4,837 | 590 | **418** | 49.1% | — |

**It does not explode** — 1.7× the runtime at a fifth of the threshold, so the practical objection
to lowering it does not hold on these files. What does hold is that almost nothing real comes back.

The features gained, judged on eight injections of one material:

| gained at | new features | present in all 8 QCs | median CV |
|---|---|---|---|
| 2,500 | 2,033 | **16.1%** | **85.3%** |
| 1,000 | 2,577 | **14.0%** | **85.0%** |

against 39.1% and 44.1% for the features already there. **Eighty-six per cent of what a lower
threshold recovers cannot be detected reliably in identical material, at a CV of 85%.** That is
noise being admitted, not sensitivity being gained.

The molecular payoff is two:

    gained  PC 47:5, PC 48:4, PS 38:2, PS 38:4, PS 40:2, SM d37:6
    lost    PC O-38:0, PC O-38:7, PE 42:4, TG 16:1_18:2_18:2

six gained and four lost, for 1,328 extra features and 70% more runtime — and the 145 extra
identified rows are almost entirely additional copies of molecules already in the table, which is
the repeated-row problem of §8 rather than new coverage.

**So 5,000 is well chosen on this instrument, and now for a stated reason rather than by
inheritance.** The justification is reproducibility in identical material, not compute cost. What
does not follow is that 5,000 transfers — it is still an absolute intensity, which is why the floor
is measured on every run.

## Guard rails

- **Bounds.** Nothing derived may leave `calibrate.BOUNDS`; a pathological file is clamped and the
  clamp is reported, rather than silently reconfiguring the run.
- **Evidence thresholds.** Each estimator declines rather than guesses — 50 matched spectra for the
  mass offset, 50 features for peak width, 20 deviations for alignment, 1,000 MS1 peaks for the
  noise floor, 3 fitted models for the retention limit. A skip is recorded with its reason.
- **Robust statistics throughout.** Medians and percentiles, never means or extremes: a tail of
  bad matches must not move the mass calibration, and one badly integrated feature must not set
  the bound that decides what a peak is.
- **Convergence.** Deriving again from calibrated data returns the same values, so two passes is
  not an arbitrary stopping point.

## Everything derived is written down

Both passes write configuration. `run_config.calibrated.json` is what pass 2 ran under;
`run_config.json` carries a `calibration_report` with the derived value, what was measured, what it
replaced, anything clamped and anything skipped.

Without that, *"run twice"* becomes *"run until it looks good"*, and a result stops being
reproducible from its configuration — which is the property the entire pipeline rests on.

## What it actually did on the validated instrument

A full two-pass run on the Kiterie positive set (1,841 s against 914 s for one pass):

    mass_offset_ppm            0.0  ->  -1.00
    chrom_fwhm                 5.0  ->   4.17 s
    min_fwhm / max_fwhm     1 / 20  ->   1.75 / 15.03 s
    align_rt_max_difference   30.0  ->   5.0 s   (CLAMPED: 4 x p95 was 2.8 s)
    noise_threshold           5000  ->   5549    (measured floor 1009)
    retention_model_max_error  skipped on pass 1 — see below

**And it did not improve the result.**

| | rows | identified | distinct | median QC CV | rows over 30% CV |
|---|---|---|---|---|---|
| one pass | 1,671 | 575 | **437** | **20.1%** | **185** |
| two passes | 1,667 | **589** | 433 | 20.7% | 199 |

Fourteen more identified rows, four fewer distinct molecules, and marginally worse precision. That
is the honest outcome and it is what should have been expected: **the configured values were fitted
against this instrument's own reference results**, so on this instrument there is nothing to
recover. `noise_threshold` landing at 5,549 against a configured 5,000 is the calibration agreeing
with the existing setting to within 11%, not finding a fault in it.

**The value is transfer, not improvement.** On a Q Exactive Plus, a 60 min gradient, or any
instrument whose intensity scale differs, the configured absolute values are the ones that would be
silently wrong. Running it here shows the machinery is correct — it derives, clamps, records and
converges — and shows that it costs a little precision when there was nothing to fix. **Enable it
when the method or the instrument changes; leave it off when running a validated method on the
instrument it was validated on.**

⚠ `retention_model_max_error` skipped on the first real run because `PipelineOutput` did not expose
the fitted models, so the estimator had nothing to read. A tested estimator that never receives
evidence is dead code with a passing test; the models are now exposed and the evidence path is
tested against the shape a real run returns.

## Using it

    "calibration": { "enabled": true }

Off by default: a two-pass run costs a full second pass, and a facility pipeline should not
silently do twice the work — nor silently change its own settings — without being asked.
