# Two-pass mass calibration — the process, with values

Every number in the validation tables was produced in one pass, uncalibrated, with a fixed 10 ppm
linking tolerance. That is fine for twelve of the fifteen studies and wrong for three, and the
reason is measurable rather than a matter of judgement.

## Why a second pass is needed at all

A fixed linking tolerance assumes every injection in a batch sits on the same mass axis. Measured
across the validation set, within-batch drift ranges over two orders of magnitude:

| instrument | within-batch spread |
|---|---|
| Agilent QTOF (MTBKS222) | 0.1 ppm |
| Thermo Q Exactive Plus | 0.2 ppm |
| Sciex X500R | 0.5 ppm |
| Agilent QTOF (MSV000094718) | 0.8 ppm |
| Sciex IMS QTOF | 0.9 ppm |
| Thermo Q Exactive HF | 1.0 ppm |
| Agilent 6545 | 1.8 ppm |
| Sciex ZenoTOF 7600 | 2.6 ppm |
| Sciex TripleTOF | 3.4 ppm |
| Sciex IMS (high mass) | 3.5 ppm |
| **Waters Xevo G2** | **9.3 ppm** |
| **Bruker maXis II** | **17.1 ppm** |

⚠ **Analyser class predicts nothing.** Three TOFs are tighter than the Orbitraps; the median TOF
spread is 1.8 ppm against 0.6 ppm for Orbitrap, and the single outlier is one instrument, not a
family. So a preset — "25 ppm for TOF, 10 for Orbitrap" — is either too loose for most instruments
or too tight for some. The spread has to be **measured per batch**.

On MTBLS5163 the consequence is concrete: three QC injections sit at +13.5 ppm against a batch
median near zero, so at a 10 ppm link tolerance their features can never join the main consensus
group. They form their own — two features instead of forty-one — leaving a zero in the main row
for every one of those injections. **48 of 336 lipids split that way; on ST004797, 84 of 701,
stranding 1,365 measurements.**

## The windows actually in force

| stage | window | configurable? |
|---|---|---|
| feature detection | 10 ppm | `features.mass_error_ppm` |
| **feature linking** | **10 ppm** | `features.link_mz_tolerance` ← where the damage happens |
| map alignment | 10 ppm | `features.align_mz_max_difference` |
| **identification search** | **20 ppm** | ⚠ hardcoded `MAX_PPM_DIFF` in `peaks.py` |
| gap filling | 10 ppm | `gap_filling.mass_tolerance_ppm` |
| adduct pairing | 20 ppm | `adducts.py` |

⚠ **The identification search has always been 20 ppm**, inherited from LipiDex as
`CDPeakFinder.MINPPMDIFF`. Widening the search is not what this is about — the search was never
the constraint. **Linking** is, and it runs at half the search window.

Verified against Skin_QEplus, the first study run with per-row theoretical masses: 99.7% of
negative and 99.0% of positive identifications fall inside 20 ppm, nothing beyond 30 ppm.

## The process

### Pass one — measure, do not correct

Run the study exactly as now, with one change:

    features.link_mz_tolerance   10 -> 25 ppm

⚠ **Loosen LINKING only.** Detection stays at 10 ppm and the search is already 20 ppm. A loose
link tolerance in pass one costs little because pass one's output is used only to pick calibrants,
and calibrants are drawn at `dot >= 900` — near-isobar merges a wide window might create do not
survive that filter. The reported numbers never come from pass one.

Its purpose is to produce a consensus in which drifted files are present rather than split off,
because a file whose features never joined a group cannot have its offset measured.

### Measure the per-file offset

    python scripts/calibrate_files.py <study> \
        --results <pass-one>/Results/<polarity>/Unfiltered_Results.csv \
        --out <study>/mzml_calibrated --measure-only

| parameter | value | why |
|---|---|---|
| calibrant selection | `dot >= 900`, seen in >= half the files (cap 5) | ~204 per study on MTBLS5163 |
| match window | 30 ppm | wider than the drift being measured, or it cannot be seen |
| RT window | ±0.15 min | tight enough that a neighbouring peak is not picked |
| minimum calibrants | 20 | below this the median is not trustworthy; the file is left alone |
| minimum correction | 2.0 ppm | below this, correcting is fitting noise |

**Two references, for two different questions.**

*Relative* — each file against the batch consensus m/z. Plentiful (~204 calibrants) and it is what
puts every injection on a common axis, which is what fixes linking.

*Absolute* — the batch against theoretical library masses, folded in as a single bias so the files
are not merely made to agree with each other while staying jointly wrong.

⚠ Until the `Theoretical m/z` column (added 2026-08-26) is present in a study's results, the
absolute reference rests on ~25 name-matched identifications rather than all ~204 — results carry a
canonicalised name (`SM d30:1`), the library the raw MSP entry (`SM d14:1_16:0`). **Re-run once
with the column live and the absolute reference becomes as plentiful as the relative one.**

### Write calibrated files

Same command without `--measure-only`. Every m/z is scaled by `1 - ppm x 1e-6` — **MS1, MS2 and
precursors alike.** Correcting the spectra and leaving the precursors behind would shift every MS2
away from the peak it was taken from, breaking identification while appearing to fix linking.

Files inside ±2 ppm are symlinked rather than rewritten. A `calibration.json` beside the output
records the applied ppm per file.

⚠ **Pass one must never read calibrated files.** Run against corrected data it finds offsets near
zero, concludes nothing needs doing, and the real drift becomes invisible — or a second correction
lands on top of the first. Neither raises anything. The script refuses to run on a directory
carrying a `calibration.json`, and that guard matters now the conversion cache is shared across
analyses of a study.

### Pass two — the reported run

Run on the calibrated files with linking returned to tight:

    features.link_mz_tolerance   25 -> 10 ppm

Tight is correct here **because the drift has been removed rather than accommodated**. A permanently
wide link tolerance would merge genuinely distinct near-isobars — the error this pipeline criticises
in other tools.

## What to check afterwards

Not "did it run" but "did it help", against the pass-one numbers:

| measure | expectation |
|---|---|
| lipids split across two rows | falls — 48 of 336 on MTBLS5163, 84 of 701 on ST004797 |
| stranded measurements | falls — 1,365 on ST004797 |
| identifications | rises slightly, as split rows merge |
| **decoy rate** | ⚠ **may move either way** — merging rows merges their decoy hits too |
| per-file spread on the calibrated files | under 2 ppm, or the correction did not take |

The decoy rate is the one to watch. If calibration merges things that should not merge, that is
where it shows.

## Which studies need it

**Three**: MTBLS5163 (17.1 ppm), MTBKS222_Waters (9.3 ppm), ST004797 (2.6 ppm but 109 injections
and 84 split lipids). The other twelve measured at 0.1–1.8 ppm and would be corrected by less than
the 2 ppm floor.

⚠ **Run all fifteen anyway, once.** Not for calibration — for the `Theoretical m/z` column, which
gives absolute mass accuracy per study and closes the mass-accuracy axis listed as never-run in
`1.md`. Calibrating three studies now would leave a table where three rows have accuracy figures
and twelve do not.
