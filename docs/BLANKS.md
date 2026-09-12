# What is actually in the blanks

Measured on the CKD rat heart validation set — 26 files per polarity, 20 samples, 4 QCs, and
two blanks at injections 1 and 51. `scripts/analyse_blanks.py` reproduces all of it from an
`Unfiltered_Results.csv`. **No filtering decision is made here**; this is the evidence for one.

## Blank signal is a small, mostly separable fraction

| | positive | negative |
|---|---|---|
| mean total signal, samples | 1.41e11 | 1.01e10 |
| mean total signal, blanks | 4.66e9 / 3.46e9 | 9.7e8 / 1.19e9 |
| blanks as a share of sample signal | **~3%** | **~10%** |
| compound groups seen in any blank | 1,310 of 13,278 | 231 of 2,557 |
| **identified** groups seen in any blank | **88 of 832** | **3 of 247** |

## The signal splits cleanly, which is the useful part

Classifying every group present in a blank by its blank-to-sample ratio:

| | positive | median RT | negative |
|---|---|---|---|
| contamination (blank ≥ sample) | 885 | 7.56 min | 173 |
| substantial (20–100% of sample) | 298 | 7.49 min | 52 |
| carryover-like (5–20%) | 29 | 8.87 min | 2 |
| trace (<5%) | 98 | 12.96 min | 4 |

Most blank signal is **contamination, not carryover** — it is at or above sample level, so it is
not analyte that leaked in from a previous injection. It concentrates around 7–8 min, the
phospholipid region.

Carryover is nevertheless visible, and exactly where the injection order predicts:

| blank | injection | contamination | substantial | carryover-like | trace |
|---|---|---|---|---|---|
| Blank_01 | 1 | 866 | 90 | 16 | 29 |
| Blank_02 | 51 | 605 | 290 | 21 | 85 |

The trailing blank has **three times as many groups sitting at 20–100% of sample level**. That
is carryover: real analyte in the wrong vial. It is why the two blanks should not be assumed
interchangeable.

**One prediction did not survive contact with the data.** Carryover was expected to be worst for
TG, which binds hardest to a reversed-phase column. It is not: TG identified in blanks sits at
0.0% of sample level. The carryover here is phospholipid, not neutral lipid.

## The threshold barely matters

| rule | groups dropped | identified dropped |
|---|---|---|
| mean of blanks, 3× | 1,156 | 23 |
| mean of blanks, 5× | 1,183 | 26 |
| mean of blanks, 10× | 1,202 | 27 |
| mean of blanks, 20× | 1,212 | 27 |
| max of blanks, 5× | 1,197 | 27 |
| **Blank_01 only, 5×** | **956** | **21** |

Positive mode. Between 3× and 20× — nearly an order of magnitude — the cost changes by four
identified compounds. The blank signal is effectively bimodal: a group is either clearly
blank-dominated or clearly sample-dominated, and very little sits in between. **The multiplier
is not a sensitive parameter on this data**, so it is not worth agonising over.

In **negative mode the question does not arise**: zero identified lipids are dropped at any
threshold tested, from 3× to 20×.

## What filtering would actually cost, named

26 identified groups fail 5× in positive mode. Nineteen are additional rows for lipids reported
elsewhere in the table, so dropping them loses no lipid. **Seven lipids would disappear
entirely:**

| lipid | RT | blank/sample | in n samples |
|---|---|---|---|
| PE O-36:5 | 8.63 | 425% | 3 |
| Plasmanyl-PC O-38:2 | 8.08 | 129% | 20 |
| PC 31:1 | 8.58 | 129% | 14 |
| PC[OH] OH-31:0 | 9.57 | 126% | 6 |
| PC[OH] OH-38:1 | 7.60 | 109% | 19 |
| TG 57:0 | 18.49 | 26% | 18 |
| PC[OH] OH-42:8 | 7.70 | 23% | 12 |

Five of the seven are **more abundant in the blank than in the samples**. Those are contaminants
that happened to match a library entry, and losing them is a gain. The two worth a second look
are `TG 57:0` and `PC[OH] OH-42:8`, both around a quarter of sample level and detected in most
samples — plausibly real signal with a contaminant contribution.

## A separate flag that came out of this

`PC[OH] OH-37:5` is reported at four retention times: 7.04, 8.26, 16.37 and 23.89 min. One
molecule elutes once. Looking wider, **15 of our 517 identified lipids are reported at retention
times more than 3 minutes apart**, `PC 38:4` at twenty different times spanning 5.2–16.0 min.

The reference does the same thing but far less: 2 lipids rather than 15, and `PC 38:4` at 6
retention times rather than 20. So this is partly inherent to matching a headgroup-dominated
spectrum, and partly **our redundant-identification filter being weaker than LipiDex's** — it
only collapses duplicates within twice the peak width of the strongest, so widely separated
repeats survive.

This is worth knowing before choosing an analysis method: for a lipid appearing in twenty rows,
which row carries the quantity is not defined by the result file.

## Options, not a recommendation

1. **Filter against both blanks at 5×.** Costs 26 identified groups, 7 lipids, 5 of which are
   contaminants. Simple, defensible, the conventional choice.
2. **Filter against the leading blank only.** Costs 21 identified groups. Treats the trailing
   blank as carryover — real analyte — and declines to penalise a lipid for having been abundant
   in the previous injection.
3. **Filter, but flag rather than delete.** Keep every row and add the blank ratio as a column,
   so the decision moves downstream where it can be made per lipid.
4. **Do not filter positive mode at all**, given 88 of 832 identified groups touch a blank and
   only 7 lipids hinge on it; inspect those 7 by hand instead.

Negative mode is unaffected by the choice.
