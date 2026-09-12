# Study metadata — ours, and how to write your own

Two different things live here, and the distinction matters:

```
metadata/
  deposit/    metadata derived from each deposit's own records — 32 files
  curated/    metadata we WROTE, because the deposit had none or its own was unusable — 7 files
  modifiers.csv   mobile-phase modifier per study, with the evidence for each
  v5_plan.csv     what each study-polarity run was actually given
```

**`curated/` is a research artefact, not configuration.** Five of these deposits cannot be processed
correctly without it, and the judgement in them is not recoverable from the raw files. They are
published so the analysis can be reproduced and, more importantly, disagreed with.

---

## The format — this is what you need for your own data

A metadata CSV is one row per injection. The **first column is the injection file name** (no
extension); the column named by `--group-by` (default `group`) carries the group.

```csv
sample,group
20240612_Ctrl_1,control
20240612_Ctrl_2,control
20240612_Treat_1,treated
20240612_Treat_2,treated
```

Run it with:

```
python scripts/run_study.py <study> --metadata-pos metadata_Pos.csv \
                                    --metadata-neg metadata_Neg.csv --group-by group
```

Extra columns are allowed and ignored by the pipeline — use them to keep the design readable
(`genotype`, `lab`, `tissue`, `window` all appear below). Only the column you name in `--group-by`
is used.

**Without metadata** every injection forms one group per role, which is stated in the log rather
than done silently. That is a safe default and a poor experiment: the presence filter keeps a
feature detected in most of *one* group, so a lipid real in one arm and absent from another is only
protected if the arms are declared.

⚠ **A group needs n ≥ 3**, and the default threshold is **0.6**. Below n=3 the fraction rule barely
constrains anything — at n=1 it is vacuous, since 1/1 is 100%. Cells smaller than that take no part,
and if that leaves none the filter reports it and skips.

The threshold is 0.6 rather than 0.7 because of what 0.7 means at small n: a group of three needed
**all three**, since 2/3 = 0.667 falls short. A triplicate arm was therefore filtered harder than a
twenty-injection one under what looked like the same rule — and most of this validation set is
small-n, so that was the common case rather than an edge case. At 0.6 a triplicate needs two of
three.

⚠ **Roles are inferred, not guessed at random.** Blanks, QCs and standards are identified from the
file name and from the `Sample Type` column of an Xcalibur sequence, and they are excluded from the
groups: a blank is what you filter against, and a feature in every Std_Mix injection is the
standard. Pooled QCs *do* form a group of their own — see below.

---

## What is in `curated/`, and why each one exists

**`MTBKS222_Agilent_metadata_{Pos,Neg}.csv`** — six replicate injections of one NIST SRM 1950 pool
per polarity, all of which classify as QC. Declared as one group (`pool`) so the filter has cells to
work on.

**`ST000991_DDA_metadata_Pos.csv`** — ⚠ one pooled QC, five injections, **three different precursor
windows** (`mz300-1100` ×2, `mz300-700` ×2, `mz700-1100` ×1). A method comparison, not a biological
design. A lipid at m/z 900 can appear in at most 3 of 5 files; one at m/z 500 reaches 4 of 5. So the
group is all five, and **0.6 is the highest threshold no window-limited lipid can fail on coverage
alone** — which is one of the two reasons 0.6 is now the pipeline default rather than a per-study
override. Grouping by window instead gives 2/2/1, all below the n≥3 floor, and no filter at all. The
window is kept as its own column.

**`ST003077_metadata_{Pos,Neg}.csv`** — seven pooled tissues, n=1 each, treated as **one group of
seven replicates** so the filter has something to work on. ⚠ The cost is stated rather than hidden: a
lipid confined to one or two tissues is detected in 1–2 of 7 and is dropped. That is the trade for
filtering an n=1 design at all. Tissue is kept as its own column.

**`MTBLS2016_metadata_Pos.csv`** — a four-lab ring trial of two genotypes with no metadata at all.
⚠ Grouped by genotype alone the labs are pooled and a lipid only one lab detects is dropped — which
in a ring trial is the result, not noise. The cell is therefore **genotype × lab**, written
pre-combined because `--group-by` takes one column; `genotype` and `lab` are kept separately for the
between-lab question. QCs form their own cells per lab.

**`ST004503_metadata_Pos.repaired.csv`** — the deposit's own `group` column reads
`QC_QC_lcqtof_t3_pos_1_041_54.raw`: a real group prefix with the file name glued on, so every row was
its own cell and the filter could not fire. Removing the file name recovers `QC_QC`,
`sample_sample`, `dda_dda` — exactly the shape its sibling deposits ST004626 and ST004650 ship
correctly. **Repaired, not discarded.**

⚠ The test for "broken" is one-directional: the group *containing* the sample name is the defect.
The reverse — a sample name containing the group, as in `Pooled_liver_pos` against group `liver` — is
a perfectly good design. Testing containment both ways rejected ST003077's real seven-tissue atlas.

---

## `modifiers.csv`

Mobile-phase modifier per study, with the **source and the evidence quoted** — which decides whether
formate or acetate adducts are searched. Three studies are marked `UNKNOWN` because the deposit
carries no structured method; they were run on the formate default and the log says so. A wrong
modifier means the wrong adducts are searched and identifications collapse, so this is recorded
rather than inferred from which setting gave more hits.

## `v5_plan.csv`

What each study-polarity run was actually given: modifier, whether the method is evidenced or
defaulted, which metadata file, the group column, and whether a presence filter was possible at all.
Written before the batch so the choice per study is inspectable, and kept after so it is
recoverable.

---

## The general lesson, if you are adapting this to your own deposits

Half the work in this validation set was not processing — it was establishing what each deposit
actually is. A pooled QC series is not a replicate series. A DDA campaign of 81 injections is not a
group of 81 replicates. An atlas with n=1 per tissue cannot support a presence filter at all.

The pipeline cannot infer any of that, and a wrong grouping fails silently — it deletes real lipids
and reports a clean run. Declaring the design is the part that needs a person.
