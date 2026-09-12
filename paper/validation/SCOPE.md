# Scope: this paper is about DDA

Decided 2026-08-22. Recorded here because the datasets, the code paths and several findings all
depend on it, and because two deposits look like successes until you ask how they were acquired.

## The rule

**Only data-dependent acquisition is processed and reported.** SWATH, other DIA, all-ion / MSe and
MS1-only injections are excluded from identification. This is already enforced in code rather than
by convention: `calibrate.searchable()` classifies EVERY FILE and the search runs only on the DDA
ones (`calibrate.acquisition`, wired at `pipeline.run`).

Excluded files are dropped from the SEARCH only. Feature detection, alignment and MS1
quantitation still run on them — those do not care how MS2 was collected.

## Why, in one paragraph

The pipeline identifies by matching one spectrum to one library entry. A DIA spectrum is the
co-fragmentation of everything in a wide window, so that premise does not hold — and the failure is
silent rather than loud. On ST000991, searching its 154 SWATH files produced **21,345 spectral hits
drawn from just 5 distinct precursor masses**: isolation-window centres that happen to fall within
a few mDa of a library entry, rediscovered in every file. Every one was rejected downstream, so the
cost was hours rather than a wrong answer — but nothing about the run said "this acquisition cannot
be searched".

## What this does NOT mean

It is not a claim that DIA cannot work. A prototype on the same data separated real fragments from
co-eluting decoys by chromatographic correlation at **AUC 0.729, 3x enrichment** — real signal,
too weak to build on without the full parameter sweep and an end-to-end retrieval test. That work
belongs in a second paper, where ST000991 is an unusually good test set: it contains the **same
pooled QC run both DDA (5 injections) and SWATH (15)**, so recovery is measurable rather than
argued. See `/home/jair/validation_bench/dia/PLAN_dia_deconvolution.md`.

## Consequences for the dataset table

| dataset | consequence |
|---|---|
| ST000991 | **5 searchable files of 159.** Agreement is 63.7% either way — the SWATH files contributed nothing — but the evidence base is 5 injections and the paper must say so |
| ST000987 | excluded entirely: all-ion, in-source CID, no precursor selection anywhere |
| MTBLS2016 | excluded entirely: no MS2 at all, despite the deposit describing it as data-dependent |
| ST004503 / ST004626 / ST004650 / ST004651 | contain MSe files; also metabolomics, two of them HILIC |
| ST004797, ST003052, ST002705, MSV000095868, MTBLS5163, ST003514 | unaffected — DDA throughout |

## The findings this scope produces, which are worth reporting

1. **Deposits mix acquisition modes within a single submission**, and no metadata field says so:
   ST000991 is 153 SWATH + 5 DDA + 1 MS1-only; ST003514 is DDA x4 + MS1-only x2; ST004503 is
   DDA x3 + MSe x3. Anyone benchmarking on public data is doing this unknowingly.
2. **Deposits mislabel acquisition.** MTBLS2016 is described as data-dependent and contains zero
   MS2 across 147 files.
3. **Reuse separates SWATH from DDA; width does not.** Measured: DDA 2.6–6.4x window reuse,
   SWATH 160.9x. A 7 Da isolation window is not evidence of DIA — MTBLS5163 fragments through
   6.9 Da windows, records 715 true precursor masses, and is ordinary DDA.
