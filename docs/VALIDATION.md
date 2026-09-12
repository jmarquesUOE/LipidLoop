# Validation: does `lipidloop` return what LipiDex returns?

The identification stage is validated. Feature detection is not, because it does not exist yet.

## What was tested

The reference export (`LIPIDLOOP_REFERENCE_DIR`, see `tests/test_validation.py`) holds four QC injections per polarity taken through
Compound Discoverer and then LipiDex, which leaves a per-file `QC_0*_Results.csv` next to each
`.mgf`. Those `.mgf` files are the *same input LipiDex read*, so running `lipidloop` on them
puts nothing between input and result except this code — no conversion, no feature detection,
no alignment. Any difference is ours.

`scripts/validate_qc.py` does the comparison; `tests/test_validation.py` pins it.

    python scripts/validate_qc.py --polarity Pos --file QC_01

Libraries: `LipidBlast_Formic.msp`, `LipiDex_HCD_Formic.msp`, `LipiDex_HCD_Hydroxy.msp` —
the three the reference run actually used, taken from the LipiDex install rather than
re-downloaded, so a mismatch could never be blamed on a different library version.
175,240 spectra after LipiDex's own peak filtering. Fatty acids from
`data/lipidex_src/FattyAcids.csv`.

**Load order matters and is not cosmetic.** It decides ties between entries scoring identically,
and the same lipid appears in more than one library with fragment masses differing in the
fourth decimal. LipidBlast first is what reproduces this reference run; in LipiDex that order
is whatever sequence the libraries were ticked in the GUI, which no result file records.

## Result

| | Pos QC_01 | QC_02 | QC_03 | QC_04 | Neg QC_01 | QC_02 | QC_03 | QC_04 |
|---|---|---|---|---|---|---|---|---|
| identifications | 2273 | 2512 | 2466 | 2483 | 1300 | 1348 | 1306 | 1271 |
| spectra identified by only one of us | 0 | 0 | 0 | 0 | 0 | 0 | 0 | 0 |
| dot product identical | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |
| reverse dot identical | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |
| delta m/z identical | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |
| identification identical | 99.65% | 99.44% | 99.59% | 99.72% | 100% | 100% | 100% | 100% |
| remainder tied on score | 8 | 14 | 10 | 7 | 0 | 0 | 0 | 0 |
| **purity identical** | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |
| **spectral components identical** | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |
| **potential fragments identical** | 100% | 100% | 100% | 100% | 100% | 100% | 100% | 100% |

Across 14,959 identifications: the same set of spectra identified, every dot product, reverse
dot product and mass error identical to the digit, and every one of the 39 differing names
explained below. On the 14,920 rows where the identification agrees, **every remaining column
is identical too** — purity, the spectral-components breakdown with each competing lipid's
percentage, and the potential-fragments list including its de-duplication. The result file is
reproduced in full. Search time is about 5 s per positive-mode file plus 3 s to load the
libraries, against a LipiDex GUI run measured in minutes.

## The 39 that differ are ties, not errors

All 39 are positive-mode, and all are cases where two library entries score *mathematically
equally* and the winner is decided by the last bit of a floating-point number.

The clearest example is Plasmenyl-PC vs Plasmanyl-PC. `Plasmanyl-PC O-18:0_18:3 [M+H]+` and
`Plasmenyl-PC P-18:0_18:2 [M+H]+` are isomers — same formula, same precursor 770.6064 — and in
`LipiDex_HCD_Formic.msp` each has exactly **one** fragment: 184.0733 for one, 184.0728 for the
other, both at intensity 999.

For a single-peak library spectrum the dot product is

    1000 * (s*l)^2 / (S * l^2)  =  1000 * s^2 / S

where `l` is the library peak's weighted intensity. It cancels. The two entries therefore have
identical scores no matter what their fragment masses are, and which one is returned depends on
whether `Math.pow` in Java and `pow` in CPython round a `**0.9` the same way — which they do
not always do. LipiDex's own choice is not consistent either: it picks the plasmenyl form in
some spectra and the plasmanyl form in others at the same mass.

Both members of such a pair are equally supported by the evidence, which is the real point:
this is a limitation of single-fragment library entries, not of either implementation. The
validation therefore accepts any member of the tied set, matched on name *and* source library —
ties can span two libraries, and the same lipid name appears in more than one.

This matters for interpretation elsewhere: an O-/P- annotation resting on a lone 184.07
headgroup fragment carries no information distinguishing the two, whatever the results file says.

The same tie surfaced a second way while validating purity, and it is worth recording because
it looks like a bug and is not. `TG 18:2_18:2_18:2 [M+Na]+` exists in both LipidBlast and
`LipiDex_HCD_Formic` with bitwise-identical dot products and fragment masses of 621.4855 and
621.4853. Only the first of two identically-named candidates survives `isUniqueLipid`, so which
library loads first decides which mass reaches the `Potential Fragments` column — and 621.4855
is then suppressed as a duplicate of a mass already printed, while 621.4853 is not. Three rows
across the eight files hung on it. Load order is therefore part of the run configuration, not a
detail.

---

# The peak finder

Validated separately and to a different standard, because it cannot be exact in the way the
identification stage is.

`scripts/validate_peakfinder.py` feeds **Compound Discoverer's own** aligned and unaligned
exports into our peak finder, so feature detection is held constant and any difference is the
peak finder's. Positive mode, the same four QC injections:

| | ours | reference |
|---|---|---|
| compound groups in | 3,781 | 3,781 |
| rows kept | 2,643 | 2,688 |
| rows identified | **440** | 440 |
| identified rows matching a reference row on quant ion + RT | 439 of 440 | — |
| **of those, identification identical** | **436 (99.3%)** | — |
| reference identifications missed | 2 | — |
| identifications the reference does not have | 2 | — |

So the join itself — which chromatographic peak an MS2 belongs to, at what retention time, with
what area — lands on the same row 417 times out of 465, and names it identically 95% of the
time. That is good enough to trust the stage and not good enough to call it exact.

Every filter category now lands in the right range:

| filter | ours | reference |
|---|---|---|
| Adduct of existing peak | 624 | 566 |
| Adduct of existing identified peak | 168 | 162 |
| In-source fragment | 136 | 146 |
| Redundant identification | 69 | 70 |
| M+1 isotope | 55 | 58 |
| Dimer | 46 | 59 |
| RT out of class range | 33 | 32 |

**Known gaps, in order of size:**

* **38 fewer rows kept than the reference** — 2,650 against 2,688. In aggregate the two remove
  almost the same rows; the split between reason labels differs, and it can, because a pair may
  satisfy more than one rule and the label recorded is the last that applied.
* **Feature width, on our own detection.** Peak widths from pyOpenMS have a median close to
  Compound Discoverer's mean (0.079 vs 0.070 min) but a long tail, mean 0.121. Every retention
  window in the peak finder scales with a group's width, so on our features the windows are
  wider than on CD's and the sweep cuts slightly harder — real-experiment recovery is 92.6%
  where CD's own features give a near-exact match. Clamping the tail was tried and moved the
  mean to 0.105 without changing recovery, so it was not kept.
* **3 identifications differ**: one isobaric candidate chosen differently (`PE 42:9` against
  `PC 37:6`), one ether assignment (`PE 38:5` against `PE O-38:6`), and one row we identify that
  the reference leaves blank.
* **2 reference identifications missed, 2 added.**

---

# End to end, with no Compound Discoverer at all

`scripts/run_pipeline.py` on the four positive QC `.raw` files — conversion, pyOpenMS feature
detection, library search, peak finder — against what Compound Discoverer plus LipiDex produced
from the same raws. **62 seconds** once the `.mzML` files exist, about 150 s including
conversion.

| | ours | CD + LipiDex |
|---|---|---|
| features per file | 5,617–5,968 | ~5,200 |
| rows | 5,517 | 2,688 |
| identified rows | 570 | 440 |
| **distinct lipids** | **442** | **366** |
| shared | **343 — 93.7% of the reference's lipids** | |
| only ours | 99 | |
| only the reference | 23 | |

Per class the two track closely: TG 150/126, PC 92/75, PE 65/55, SM 23/22, PI 10/10.

**93.7% recovery is the honest headline, and the 99 extra lipids are not yet characterised.**
They are either genuine additional identifications — pyOpenMS finds somewhat more features, and
the peak finder keeps more rows — or false positives. Deciding which needs manual inspection of
a sample of them, and until that is done the extras should be treated as unverified rather than
as a sensitivity gain. The 23 missed include `D5TG 52:2`, an internal standard, which is worth
chasing on its own.

---

# On a real experiment

`scripts/validate_experiment.py` runs the whole pipeline over a real sequence — the CKD rat
heart set, 26 files per polarity: 20 samples, 4 QCs spread through the run (injections 3, 13,
27, 49) and 2 blanks (injections 1 and 51). Four identical QC injections cannot test
retention-time alignment, because there is nothing to align; this can.

**It immediately found two defects the QC-only test could not.**

1. `MapAlignmentAlgorithmPoseClustering.align()` only *computes* the transformation — it fills
   the `TransformationDescription` and leaves the map untouched. Applying it needs a second call
   through `MapAlignmentTransformer`. Without that the pipeline ran cleanly, produced plausible
   output, and had done **no alignment at all**. The only symptom was that measured and aligned
   retention times were identical, which identical injections cannot reveal.
2. Even once alignment worked, the measured retention time was being overwritten with the
   aligned one, so the peak finder's retention correction had nothing to fit.

Then the *working* alignment made the spread worse — 6.8 s to 12.5 s — because the
pose-clustering defaults pair features up to **100 s and 0.3 Da** apart. Built for maps
disagreeing by minutes; on high-resolution data with seconds of drift it pairs unrelated
features. Tightened to 30 s / 10 ppm.

**And a third: every feature was labelled positive polarity.** pyOpenMS reports a feature's
charge as a positive magnitude whatever the run, so taking the sign from the charge is
accidentally right in a positive run and catastrophic in a negative one — no negative-mode
identification could attach to any feature, and negative mode returned **zero** identified
compounds while otherwise running perfectly. Polarity now comes from the mzML.

## Result, positive mode

| | |
|---|---|
| RT spread of a compound across files, median | 5.6 s → **5.5 s** after alignment |
| same, 90th percentile | 17.5 s → **17.2 s** |
| alignment moves features by | median 0.8 s, 90th pct 2.6 s |
| QC drift, injection 3 → 49 | **+2.6 s** |
| distinct lipids | **499** vs reference 444 |
| shared, same molecule | **416 — 97.7% of the reference's lipids** |
| shared, exact name | 425 — 95.7% |
| genuinely absent | **10 molecules** |

93.7% is the same figure the QC set gave, which is the point of measuring it twice on
independent data.

## Result, negative mode

| | |
|---|---|
| RT spread, median | 3.6 s → **2.6 s** after alignment |
| same, 90th percentile | 10.2 s → **9.4 s** |
| QC drift, injection 3 → 49 | +1.6 s |
| distinct lipids | **231** vs reference 211 |
| shared, same molecule | **183 — 95.3% of the reference's lipids** |
| shared, exact name | 189 — 89.6% |
| genuinely absent | **9 molecules** |

Alignment helps more in negative mode than positive: median spread 3.6 s to 2.3 s, against a
marginal improvement in positive.

**Recovery is measured two ways, and the difference matters.** Comparing identification strings
exactly counts a molecular-versus-sum reporting difference — `TG 12:0_16:0_20:4` against
`TG 48:4` — as both a miss and an extra, when it is the same molecule named at two resolutions.
Collapsing both sides to sum composition is the honest number, and it is 2-6 points higher.

**The honest reading is that this chromatography barely needs alignment.** 2.6 s of drift across
49 injections on a 25-minute gradient is very stable, and the residual 5.5 s spread is dominated
by apex determination at ~5 MS1 points per peak, not by drift. Alignment now does no harm and a
marginal good; it is not what limits this data. Feature detection is.

## Blanks behave exactly as their positions predict

| blank | injection | groups | fraction of its signal shared with samples |
|---|---|---|---|
| Blank_01 | 1 | 1,001 | **9%** |
| Blank_02 | 51 | 1,001 | **35%** |

The leading blank measures contamination; the trailing one measures carryover, and a third of
its signal sits in compounds that are strong in the samples. **Filtering against the trailing
blank at the same multiplier would strip real lipids.** At 5× against both, 1,881 groups were
dropped, 32 of them identified.

## What is not validated

* **Feature detection per feature.** ~5,900 per file against Compound Discoverer's ~5,200 is
  the right order, but no feature-by-feature comparison has been run. Different algorithms will
  not agree feature-for-feature; recovery of identified compounds, above, is the useful measure.
* **The ~63 extra molecules.** Unverified, not a claimed improvement.
* **Four identified lipids have no feature to attach to.** Diagnosed in
  [FEATURE_DETECTION.md](FEATURE_DETECTION.md): one never traced, one dropped for having no
  isotope partner, and two whose trace is consumed as the M+2 isotope of the homologue with one
  more double bond — a systematic collision, since 2 H and 2 x 13C differ by only 11 ppm.
* **Feature width.** Peak widths are bounded at 20 s FWHM, about four times both the observed
  median (4.7 s) and the `chrom_fwhm` parameter, which drops 3.2% of features and brings the
  mean from 0.121 to 0.102 min against Compound Discoverer's 0.070. It did not move recovery,
  so width is not what limits this.
* **`d5TG` is excluded from the library by default.** It is 22,960 entries in
  `LipiDex_HCD_Formic.msp`, 13% of the whole library, and a deuterated standard set: unless d5
  standards were spiked, a match is a mass coincidence. Two of the reference's lipids are
  D5TG and are therefore deliberately not reported. Turn `excluded_classes` off if you spike
  SPLASH.
* **Which blanks to filter against.** The default filters against both. Given the carryover
  measured above, filtering against the leading blank only — or using `statistic: max` — is
  arguably the better default, but that is a judgement about this dataset rather than a
  general result.
* **Alignment.** Four QC injections are the right test for identification — identical samples,
  so any difference is a bug — and the wrong one for retention-time alignment, which needs real
  samples with real drift. This is the outstanding item before the pipeline should be used on
  a real experiment.
