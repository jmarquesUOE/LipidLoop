# What the filters remove, and the one that is wrong

Every filter in the peak finder was reproduced from LipiDex rather than designed here. Most of
them are deduplication and are plainly right: an M+1 isotope, an adduct of a peak already in the
table, a dimer and an in-source fragment are all the same molecule counted twice. One is not.

## Where identifications go, skin organoid set, 63 files per polarity

| filter status | positive, identified | negative, identified | positive, unidentified |
|---|---|---|---|
| **kept** | 1,191 | 614 | 9,627 |
| **RT out of class range** | **131** | **54** | 0 |
| Redundant Identification | 79 | 11 | 0 |
| RT model outlier | 24 | 3 | 0 |
| Blank | 22 | 2 | 567 |
| In-source fragment | 3 | 33 | 180 |
| Adduct of existing peak | 0 | 0 | 11,469 |
| Dimer | 0 | 0 | 2,033 |
| M+1 isotope | 0 | 0 | 1,616 |

**`RT out of class range` is the single largest destroyer of identifications** — 9.0% of identified
rows in positive mode, 7.5% in negative — and it is the one that is wrong.

## Exactly what it does

`_check_class_rt_distribution`, reproducing `CDPeakFinder :: checkClassRTDist`. Five steps:

1. Collect every compound group that carries an identification (`final_lipid_id is not None`) and
   bucket it by **lipid class**. Unidentified groups are never touched.
2. **Skip any class with four or fewer** identified groups — no distribution to speak of.
3. Take each group's **aligned retention time** and compute the class **median**, then the
   **median absolute deviation** about it.
4. `window = |multiplier × MAD / 0.6745|`. The 0.6745 rescales a MAD to a standard deviation for a
   normal distribution, and the multiplier is **2.0**, so the half-width is **2.965 × MAD** —
   nominally ±2σ.
5. Any group whose retention falls outside `median ± window` has `keep = False` and
   `filter_reason = "RT out of class range"`.

It runs **after** the deduplication filters and **before** the retention model, so a row it removes
never reaches the model that could have vouched for it.

Two things it is not. It is not the retention *model* filter — that is
`_filter_by_retention_model`, which is chemistry-aware and reports `RT model outlier`. And it is
not the LOESS retention correction, which adjusts retention times and removes nothing.

### Where the 2.0 comes from

`Utilities.MINRTMULTIPLIER` is 0.5 in the source and is overwritten at startup from the GUI's
**FWHM** spinner, default 2.0. `CDPeakFinder` has a *second* field of the same name, set from its
own spinner (default 3.5) and passed to `checkClassRTDist` as an argument **the method never
reads** — dead in the original, so not reproduced. So the window is 2.965 × MAD rather than the
5.19 × MAD the visible spinner implies, and no LipiDex user has ever been able to widen it from
that control.

### What that produces on real data

Per class, skin organoid positive mode:

| class | identified | median RT | MAD | window ± | class actually spans | **cut** |
|---|---|---|---|---|---|---|
| PC | 358 | 12.74 | 1.47 | **4.34** | 2.24 – 23.34 | **53** |
| PC[OH] | 104 | 11.38 | 0.63 | **1.88** | 9.93 – 20.50 | **17** |
| Plasmanyl-PC | 176 | 15.01 | 1.79 | 5.30 | 10.13 – 22.51 | 16 |
| DG | 47 | 15.21 | 1.43 | 4.25 | 11.42 – 23.12 | 7 |
| PI | 60 | 11.85 | 0.52 | **1.54** | 10.74 – 14.95 | 5 |
| LysoPC | 33 | 9.74 | 0.57 | 1.70 | 8.90 – 12.38 | 5 |
| SM | 81 | 14.72 | 2.28 | 6.76 | 10.32 – 20.38 | 0 |

`PC[OH]` is the clearest case: the identifications cluster tightly, MAD 0.63, so the window is
±1.88 min — while the class genuinely elutes across 10.6 minutes. Everything past 13.3 min goes.

**The MAD is what makes this self-defeating.** It is robust by design, so the tails cannot widen
the window that is cutting them — the tighter the middle of a class clusters, the more aggressively
its ends are removed, no matter how many real long-chain species are out there.

## Why it is wrong

It keeps identifications within a robust two-sigma window on the *distribution* of retention times
within a lipid class.

In reversed-phase lipidomics that is the wrong shape for the data. Retention inside a class is not
a cluster with a centre and noise — it rises almost linearly with carbon number. On this method
the fitted per-class models run at **0.6 to 0.7 min per carbon**, so a class spanning C30 to C50
spans eight minutes of gradient, legitimately. A two-sigma window on that distribution **clips both
chain-length tails by construction**, and the longer tail is the one anyone reaching for
very-long-chain species cares about.

The module docstring already said the filter "is a clustering test, not an elution-order test".
This is what that costs. In positive mode `PC 42:1` is removed at RT 18.18 with a **dot product of
1000** and an area of **5.5e8** — a perfect spectral match on one of the largest peaks in the run,
discarded for eluting 5.4 min after the median of its class.

**It is not useless.** The same PC window also removes `PC 36:2` and `PC 34:1` at RT 2.2 with areas
of 4e4 — void-volume junk carrying a confident-looking spectrum, which is exactly the case the
filter was written for. The problem is that a symmetric window around a median cannot tell the
early-eluting artefact from the late-eluting long chain, because it is not looking at chain
length at all.

### The evidence

Tracing every reference molecule we failed to report (`scripts/trace_missing.py`):

| | missing | RT out of class range | identified, no feature | never fragmented |
|---|---|---|---|---|
| skin, positive | 25 | **12** | 10 | 2 |
| skin, negative | 11 | **4** | 3 | 2 (+2 blank) |

And the missing molecules are **not** what you would guess:

| | median carbons | median reference area |
|---|---|---|
| missing, skin positive | **43** | 6.9e7 |
| recovered, skin positive | 37 | 1.4e8 |
| missing, CKD positive | **50–57** | up to 4.2e9 |
| recovered, CKD positive | 42 | 6.2e7 |

**Chain length separates them; abundance does not.** Four of the twelve clipped in skin positive
mode are among the largest peaks in the entire reference — `Plasmanyl-PC O-46:1` at 9.2e8,
`PS 42:1` at 9.0e8, `PC 44:2` at 6.0e8, `PE O-40:2` at 4.6e8. These are not a marginal tail being
tidied away; they are big, well-measured lipids removed for eluting where a long chain elutes.

The same pattern holds on the CKD rat heart set, where the missing molecules have a median of 50
carbons against 42 for the recovered ones and include peaks at 4.2e9 and 3.6e9.

### The pipeline already has the right test

`retention_model_filter` fits `retention ≈ intercept + per_carbon·C + per_double_bond·DB` per class
(R² 0.92–0.97 on this data) and removes identifications the model cannot place. It knows that a
long chain elutes late because it is long. It flagged **24 rows in positive mode against the class
window's 131**, and of the 25 missing reference molecules it accounted for exactly one.

The two filters disagree because one models the chemistry and the other models a histogram.

## What to do

Three options, in order of preference. All three were run end to end on the skin set, so the trade is measured rather than argued:

| | default | **`rt_filter: false`** | `keep_filtered: true` |
|---|---|---|---|
| skin positive, recovery | 95.8% | **98.0%** | 98.5% |
| missing | 25 | **12** | 9 |
| rows in the table | 10,818 | **10,927** | 27,633 |
| skin negative, recovery | 96.9% | **98.0%** | 99.4% |
| missing | 11 | **7** | 2 |
| rows | 3,786 | **3,826** | 5,779 |

**1. Turn the class window off and keep the retention model.** `"rt_filter": false`. The surgical
fix: it removes the filter that is wrong and keeps the one that is right. It costs **109 extra rows
in positive mode and 40 in negative** — the identifications the window was clipping, less the ones
the model then judges properly — and buys most of the recovery. This is the setting to use unless
you have a reason not to.

**2. `"keep_filtered": true`.** Nothing is dropped; `Final_Results.csv` carries every compound
group with a `Filter Status` column saying what *would* have been removed. This is what the skin
run in `OmicsPilot_Lipidomics/` currently uses, because it was asked for explicitly. It recovers
another 0.5–1.7 points, and it costs **16,800 extra rows in positive mode**.

**Read the caveat before using that table.** Keeping everything also keeps the rows the
deduplicating filters exist to remove: in positive mode 11,469 adducts, 2,033 dimers and 1,616
isotopes of peaks already present. Those are genuinely the same molecules counted again, and
anything that sums or counts the table without reading `Filter Status` will multiply-count them.
`Filter Status` blank means kept.

**3. Raise `rt_filter_multiplier`** from 2.0. Blunt — it widens the window for every class equally,
including the ones where it is doing useful work.

## The right model for the job — and it is already here

The void-volume PCs *should* go. The question is what removes them without also removing
`PC 42:1` at 18.18. The answer is the retention model the pipeline already fits, used **instead
of** the class window rather than after it.

`PC = -1.73 + 0.477·C - 0.737·DB`, fitted on this run (R² 0.909, residual sd 0.925 min):

| row | observed RT | model predicts | residual | verdict |
|---|---|---|---|---|
| PC 36:2 | 2.24 | 13.97 | **−11.73** | removed |
| PC 34:1 | 2.27 | 13.75 | **−11.48** | removed |
| PC 34:1 | 7.99 | 13.75 | **−5.76** | removed |
| PC 34:2 | 7.87 | 13.02 | **−5.15** | removed |
| PC 34:1 | 12.74 | 13.75 | −1.01 | kept |
| PC 44:2 | 18.06 | 17.79 | +0.27 | kept |
| PC 42:1 | 18.18 | 17.57 | +0.61 | kept |
| PC 43:1 | 18.91 | 18.05 | +0.86 | kept |

Run end to end with `rt_filter: false`, **every identified row below RT 8 is removed as an
`RT model outlier`** — the three at 2.2 and, more usefully, a second cluster at 7.87–7.99 that
looks perfectly good by score (`PC 34:1`, dot 997, area 2.3e6) and is 5.8 min from where a 34:1
PC belongs. The class window happened to catch those too, but only because they fell outside a
window drawn for unrelated reasons.

Everything that survives below RT 9 is chemistry: acylcarnitines `AC 12:0`–`AC 16:1`,
`LysoPC 14:1`, `SP d18:1`. Short, polar, genuinely early.

**Why it works where a window cannot.** The window knows only "PC" and a median. The model knows
this row claims to be a 36:2, and a 36:2 PC elutes at 13.97 min on this gradient. That makes the
same test one-sided in the way chemistry actually is: a lipid can be far too early (it is not what
the name says, or it is an in-source fragment), and being *late* for its composition is equally
diagnostic, but neither has anything to do with where the class median sits.

### What makes it safe to remove data with

- fitted **only on well-detected rows** (present in at least half the samples), so the duplicated
  and one-off rows that most need judging cannot spoil the fit meant to judge them — PC goes from
  R² 0.70 to 0.909 on this alone;
- **iteratively trimmed** at 3σ, up to three passes;
- a class whose model does not fit — fewer than 6 points, R² < 0.7, residual sd > 1 min, or
  coefficients of the wrong sign — **is not used to judge anything**. A bad fit is an absence of
  evidence, not evidence of absence;
- the limit is `max(3 × residual sd, 1.0 min)`, so a very tight class cannot acquire an absurdly
  narrow window.

### Where it is still improvable

- **The fit is least squares inside the trim loop.** A grossly wrong point perturbs the first
  iteration before it is trimmed. Theil–Sen or Huber weighting would be immune from the start.
  Three trim passes cope here, but this is the weak link.
- **Retention against carbon is not quite linear.** The saturated Cer[ADS] ladder shows smooth
  systematic residuals under a straight line. A quadratic term, or an equivalent-carbon-number
  form, would shrink the residual sd and let the threshold tighten — PC's ±2.78 min limit is set
  entirely by its 0.925 min residual sd.
- **A void-volume floor** would be the obvious backstop for classes too small to model, and is
  physics rather than statistics: nothing retained on a C18 column can elute at t₀. It is **not
  implemented, because on this data it would remove nothing** — with `rt_filter: false` no
  identified row survives below RT 8 that should not. Worth adding the day a class with fewer than
  six identifications catches void junk.
- **It is necessary, not sufficient.** An in-source fragment that happens to elute where its
  mis-assigned name predicts will pass. Retention is one orthogonal axis, not a proof.

### One thing that was quietly missing

`parse_sum_name` only accepted sum compositions, so **every row reported at molecular resolution
was silently exempt from the model** — 78 identifications on the skin set, every TG and DG whose
fatty-acid purity cleared 75%, plus the gangliosides. They were not being protected; they were
being skipped. Molecular names are now collapsed first (`TG 16:0_18:1_18:2` → `TG 52:3`), which is
sound because retention depends on total carbons and double bonds, not on how they are split
between chains, and the model has no term that could tell two splits apart.

The effect is small and in the right direction: one more usable class per polarity, two more
outliers caught in positive mode and fifteen more in negative, and negative-mode recovery from
97.7% to **98.0%**. Nothing that was being kept correctly started being removed.

## Reading the scores

`Final_Results.csv` now carries `Dot Product`, `Reverse Dot Product` and `Purity` for the name as
reported (`write_scores`, on by default). They matter most when nothing has been filtered out,
because a weak row and a strong one then sit side by side.

On the skin set the identified rows run from a dot product of **501 to 1000**, median 975
(positive) and 895 (negative), with **195 positive and 116 negative rows below 700**. A row at 520
and a row at 990 are not the same kind of statement, and before these columns existed the table
did not say so.

The scores describe the name that was written, not the best candidate available: a row reported at
sum composition takes the best score among every candidate that collapsed to that sum, because
they all voted for it. Crediting it with the top molecular candidate's score would report
confidence in a chain assignment the row deliberately declined to make.

## The other two loss routes

`RT out of class range` is the fixable one. The rest of the missing molecules are:

- **identified, no feature** (10 positive, 3 negative): spectra were identified but no MS1 feature
  carried them, so there was nothing to quantify. These are late-eluting and long — `CE 26:0` at
  25.5 min, `TG 60:9` at 24.9 min on a 30 min gradient. Feature-detection sensitivity, and
  `max_fwhm` (20 s) is the first parameter to suspect where peaks broaden late.
- **never fragmented** (2 per polarity): no MS2 anywhere in the batch was identified as them.
  Nothing in processing could have recovered these; they are what an inclusion list or a freed
  duty cycle is for — see [EXCLUSION_LIST.md](EXCLUSION_LIST.md).

## ⚠ The trim scale must be robust, or contamination buys immunity (fixed 2026-08-13)

`fit_model` trims outliers iteratively at `trim_sigma * spread`. `spread` was the standard
deviation of the residuals — **which the outliers themselves set**, so the threshold meant to catch
them was inflated by them.

Found by reading a QC report: `PS 40:6` reported at **RT 24.47 min** on a 25 min gradient, where PS
elutes at 8–11. Six identified rows survived past 24 min spanning 0.09 min between six distinct
molecules, `PS 34:1` through `PS 44:10` — ten carbons and nine double bonds apart, eluting 0.02 min
from each other. That is the end-of-gradient wash, not chromatography.

The mechanism is the part worth remembering:

    PS residual spread with those rows in   6.50 min
    trim threshold = 3 x 6.50               19.5 min   -> nothing is ever trimmed
    resulting fit                           R2 0.020   -> class declared NOT USABLE
    an unusable class is exempt             -> the junk that broke the fit escapes the filter

**Contamination bad enough to break the model bought the whole class immunity from the model**, and
the worse the junk the more certain its survival.

**Fix: a robust scale.** `spread` is now the median absolute deviation of the residuals, scaled by
1.4826. On the same PS data it trims 12 of 29 points and fits at **R² 0.942, sd 0.218 min**.

    usable class models   10 of 17  ->  14 of 17
    PS                    R2 0.020  ->  0.942
    PE O-                 R2 0.000  ->  0.904

Every class that changed did so by rejecting junk, not by loosening a standard, and the classes
that already fitted are unaffected.

⚠ **The general lesson:** any threshold computed from the data it is meant to filter must use a
robust estimator, or the thing being detected sets its own detection limit. The same fault appeared
independently in `split_peaks` (an infinite CV read as "very bad" rather than "unmeasurable") and in
`qc.outliers` (classical Hotelling's T² masking multiple outliers, now paired with a robust form).
