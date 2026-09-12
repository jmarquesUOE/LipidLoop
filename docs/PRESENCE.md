# Presence filtering

A feature detected in a handful of injections is a detection, not a measurement. Deciding which
rows are analysable is a filter on presence — and the parameter that matters is not the threshold
but **the group the threshold is applied within**.

## The group matters more than the threshold

Identified rows kept, Kiterie brain positive, 891 identified rows, 40 samples in a balanced
2 × 2 × 2:

| group used | 50% | 60% | 70% | 80% | 90% | 100% |
|---|---|---|---|---|---|---|
| genotype (2 × 20) | 603 | 567 | **534** | 498 | 457 | 394 |
| age × genotype (4 × 10) | 664 | 621 | **585** | 547 | 513 | 471 |
| age × sex × genotype (8 × 5) | 713 | 713 | **626** | 626 | 544 | 544 |

At a fixed 70% the grouping moves the answer by 92 rows; moving the threshold from 50% to 100%
inside the middle row moves it by 193. Both matter, and only one of them usually gets stated.

**Note the repeats in the 8 × 5 row.** With five samples per cell there are only six reachable
thresholds, so 50% and 60% are the same filter, as are 70% and 80%, and 90% and 100%. A threshold
tuned to a precision the design cannot express is a decision that was never made.

**Configured here: `age × genotype`, 70%** — the group is the unit being compared, and with ten
samples per group one detection does not decide the outcome.

## Why not filter on the QCs

Requiring presence in all eight pooled QC injections looks stricter and cleaner:

| | rows kept | distinct | residual missing | rows fully complete |
|---|---|---|---|---|
| ≥ 70% of one age × genotype group | 585 | — | ~13% | ~52% |
| present in 100% of QCs | 429 | 357 | **2.8%** | **74.6%** |

It is also **structurally biased against the result the experiment is for**. A pooled QC mixes
every sample, so a lipid confined to one group of ten is diluted tenfold in the pool. Of the 198
positive rows that pass the group rule and fail the QC rule, **25% are present in over 80% of one
group and entirely absent from another** — 49 rows in positive, 27 in negative — and they are
4–6× fainter than the rows that pass. The QC rule does not drop them for being unreliable. It
drops them for being group-specific.

The two rules are also nested rather than complementary: `A AND B` keeps 428 where B alone keeps
429, and `A OR B` keeps 627 where A alone keeps 626. Choosing between them is a stringency
decision, not a choice of criterion.

**Use QC completeness as an annotation, not a filter** — flag which retained rows are
QC-complete so a reader can weight them.

## Honest cost of the group rule

The other 75% of those extra rows are not group-specific; they are faint and patchy. The rule
admits noise, and that is the price of keeping the on/off species.

## The gaps are measured, not imputed

Filling is `gapfill.py`'s job and it does not substitute a statistic — it goes back to the raw file
and re-integrates at the expected mass and retention time, which is what Compound Discoverer's
`Fill Gaps` node does. See [GAP_FILLING.md](GAP_FILLING.md).

That division matters for one claim this filter makes. It deliberately keeps a lipid present in one
group and absent from another, and the tempting next step is to protect that whole-group zero from
being filled. But protecting it **assumes** the absence is real rather than a detection failure,
and only the raw file can tell those apart. So the assumption is not made: those positions are
measured like any other, and the run reports how many turned out to carry a peak.

The filled table is **`Final_Results_Filtered.csv`** — the analysis-ready one.
`Final_Results.csv` keeps its zeros, so the two can be diffed to see precisely what was filled.

## Configuration

    "presence_filter": {
      "enabled": true,
      "min_fraction": 0.7,
      "metadata": "metadata_Pos.csv",
      "group_by": ["age", "genotype"]
    }

`metadata` is a CSV whose first column is the injection file name. With no metadata the samples
form one group and the rule degrades to "detected in `min_fraction` of all samples" — which is
reported in the log rather than done silently.

Applied to **every** compound group, not only the identified ones: how reliably a feature was
detected is a property of the measurement, and an unidentified row surviving on two detections
while a named one is dropped on nine would be indefensible.
