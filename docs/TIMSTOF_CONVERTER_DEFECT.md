# The timsTOF arm is broken in the converter, not the data

Measured 2026-08-26 on `sample_plasma_1-3_1_7192` (MTBKS222 BrukerPASEF arm).

## What the funnel looked like

| arm | spectra loaded | pass score | associated with a feature |
|---|---|---|---|
| Thermo | 2,988 | 2,324 (78%) | 2,043 |
| Sciex IMS | 5,570 | 4,159 (75%) | 3,803 |
| **timsTOF** | 8,624 | **294 (3.4%)** | **76** |

Scored against its own MAF that gives 5.6% agreement. ⚠ **That number is meaningless and must not
be quoted** — it measures the converter, not the instrument.

## The tell

| | timsTOF | Thermo |
|---|---|---|
| forward dot, median | **8** | 926 |
| reverse dot, median | **1000** | 995 |
| mass error, median | 4.8 ppm | 2.8 ppm |

Reverse 1000 with forward 8 is diagnostic. `score.py` reproduces LipiDex exactly: `reverse=True`
counts only matched peaks, while the forward score adds every unmatched peak — library or sample —
to the denominator. So the matched peaks align perfectly and almost everything else does not.

## What it is not

Three hypotheses tested and killed, recorded so they are not tried again:

1. **Duplicated peaks from summing frames** — median fraction of adjacent peaks within 20 ppm was
   0.00. ⚠ But see below: the window was too tight and this conclusion was wrong for the wrong
   reason.
2. **MS2 not corresponding to MS1** — 320 of 321 sampled precursors have an MS1 peak within
   20 ppm, and the RT axes align (MS1 0.0-25.0 min, MS2 0.2-25.0 min).
3. **Spectra too sparse or too rich overall** — raw, unbiased, timsTOF carries a median of 23 peaks
   per MS2 against Thermo's 17, over 64,033 spectra against 3,271. Richness is fine.

## What it is

**58.7% of timsTOF peaks are dropped by the `precursor - mz > 1.5` filter, against 5.1% for
Thermo.** They sit at or above the precursor. A singly-charged lipid fragment cannot.

Three unrelated precursors — 821.5329, 769.4473, 686.4392 — return the same top peaks:

    514.9x, 432.9x, 596.8x, 350.98, 268.99, 186.99      spaced ~81.94 apart

A background ion series, identical across precursors that share nothing. The per-precursor mobility
slice is not isolating that precursor's fragments.

And the same nominal ion recurs **within one spectrum** at 514.873 / 514.913 / 514.927 / 514.931 /
514.935 / 514.942 — a spread of 69 mDa, 134 ppm. ⚠ This is why hypothesis 1 above was dismissed
wrongly: the duplicate test used a 20 ppm window, which is 0.01 Da at m/z 515 and cannot see a
0.07 Da smear. The frames being summed are not on a common mass axis, so concatenating their peak
lists produces a smeared cluster where there should be one peak.

## Where to look next

`scripts/tdf_to_mzml.py`, the MS2 loop. Two candidates, not yet separated:

- the scan range. `data[frame, ScanNumBegin:ScanNumEnd]` may not be selecting what the fragment
  frame's precursor table says it is — if the slice is wrong or ignored, every spectrum gets the
  frame's dominant background, which is exactly what is observed.
- the frame summing. Peaks from different frames arriving on different mass axes, then
  concatenated rather than merged onto a common axis.

The check that separates them: take one MS2 frame carrying several precursors and confirm that
different scan ranges return *different* peaks. If they do not, it is the slice.

⚠ Until this is fixed the arm produces 20 identifications from 8,624 spectral matches, and no
timsTOF identification count, agreement figure or FDR from this pipeline is usable.


---

# RESOLVED IN PART — 2026-08-26

## What it actually was

**Frame summing by concatenation.** PASEF re-acquires a precursor across consecutive frames; the
same ion lands in a slightly different TOF bin each time. The converter concatenated the peak lists,
so a real fragment was written 3-4 times at slightly different m/z with identical intensity:

    782.5156 (6909)  782.5067 (6909)  782.4713 (6909)  782.4669 (6909)

⚠ **This was the FIRST hypothesis, dismissed above on a test that could not have detected it.** The
duplicate check used a 20 ppm window; the jitter is 45-63 ppm. The clean 0.00 it returned was
meaningless, and three further hypotheses were pursued before returning here. If you re-test this,
use a window wider than the jitter.

Fixed in `scripts/tdf_to_mzml.py::_merge`, 60 ppm, intensities summed, intensity-weighted centroid.

| measured on sample_plasma_1-3_1_7192, sampled across the whole run | before | after |
|---|---|---|
| peaks within 60 ppm of a neighbour | 0.06 | **0.00** |
| forward dot, median, among hits | 9 | **22** |
| hits reaching dot >= 500 | 6% | **21%** |

## Two earlier claims corrected

- **"3.4% of spectra pass scoring, against 75-78% for the other arms."** The denominator is not
  comparable. timsTOF writes 64,033 MS2 per file against Thermo's 3,271, because PASEF has the duty
  cycle to fragment background ions the Orbitrap never reaches. In absolute terms the hit counts are
  similar — roughly 1,150 against 942 on one file. A low percentage here is partly the instrument
  working as designed.
- **"The MS2 m/z axis is wrong."** It is not. That came from sampling `mz[:40]` — the lowest-m/z
  peaks, which are the noise region. Across all peaks the most common MS2 value in the file is
  184.1, the phosphocholine head group, followed by 760.6 / 786.6 / 758.6 / 810.6.

## What is still wrong

Forward dot 22 against Thermo's 925, with reverse dot still a perfect 1000.

**58.7% of timsTOF peaks sit within 1.5 Da of the precursor or above it** (Thermo: 5.1%). The
`precursor - mz > 1.5` filter therefore strips most of the spectrum, leaving a median of ~5 peaks
against Thermo's 13. A sparse sample matched to a 20-peak library entry scores a perfect reverse dot
— every sample peak is in the library — and a near-zero forward dot, because the unmatched library
peaks fill the denominator. That is precisely the signature still observed.

⚠ So the open question is no longer "why do the spectra score badly" but **"why is so much of each
spectrum in the precursor region"**. Candidates, untested:

- low fragmentation efficiency at the collision energies used, leaving the precursor intact;
- the mobility slice carrying co-eluting ions of similar mobility but higher m/z;
- isotope envelopes of the precursor being retained where the other vendors' centroiding removes them.

**Until that is settled, timsTOF identification counts remain provisional.** The arm is queued so it
produces a row on current code, not because it is on equal footing with the other four vendors.

---

# DECISION, 2026-08-27 — the pipeline is not valid for timsTOF PASEF data

**The arm is kept as it is and is not pursued further.** Not because the converter is unfixable, but
because the pipeline's assumptions do not hold for PASEF, and making them hold is a different
project from the one this software is.

## What was actually established

The converter defect was real and is fixed: frames were concatenated rather than merged, so every
fragment appeared 3–4 times at slightly different m/z. Merging at 60 ppm moved the forward dot from
9 to 22 and raised the share of hits passing `dot >= 500` from 6% to 21%.

**It was not enough, and the remainder is not a converter problem.** After the fix:

| | timsTOF | Thermo |
|---|---|---|
| forward dot, median among hits | **22** | 925 |
| reverse dot | 1000 | 995 |
| peaks within 1.5 Da of the precursor or above it | **58.7%** | 5.1% |
| scored peaks per MS2 after the precursor filter | **~5** | 13 |
| delivered identifications, whole arm | **13 lipids, 6 injections** | 322 |

## Why this is a scope boundary, not a bug

PASEF is not conventional DDA and the pipeline is built for conventional DDA:

- **One spectrum, one precursor.** A PASEF frame carries fragments from many precursors separated
  in the mobility dimension. Our conversion sums that dimension away, so the separation survives
  only because PASEF selected the pairs as separate precursors in the first place. Nothing
  downstream uses mobility.
- **The precursor region dominates.** 58.7% of peaks sit at or above the precursor, so the standard
  `precursor - 1.5 Da` filter strips most of each spectrum. A sparse query against a rich library
  gives a perfect reverse dot and a near-zero forward dot — which is exactly the observed signature.
  Changing that means changing the scoring, which changes every other arm too.
- **The libraries are HCD Orbitrap.** timsTOF CID fragmentation is chemically different, so even a
  perfect converter is matching against spectra from another instrument class.

⚠ Handling PASEF properly means carrying mobility through feature detection, association and
scoring — a CCS axis alongside retention. That is a genuine and worthwhile project (see the
instrument grant work, where mobility resolves 1,756 co-eluting pairs per injection at 99% MS2
cosine below 0.5). It is not a fix to this pipeline.

## What may and may not be said

- ✅ The **duty-cycle and mobility measurements** stand — they are measured from the raw data and do
  not depend on our identification working.
- ✅ The converter reads TDF natively via alphatims, with no Bruker SDK and no Windows. That stands.
- ⚠ **No timsTOF identification count, agreement figure or FDR from this pipeline may be quoted.**
  13 lipids from 6 injections is a statement about scope, not about the instrument.
- ⚠ The arm must not appear in any cross-vendor identification table. Excluded in
  `manuscript/validation/score_agreement.py::EXCLUDED_ARMS`.
