# The MS2 window, measured

`ms2_tol` is an absolute m/z window and belongs to the analyzer that produced the fragment spectra. `ms2_calibration.py` measures it on every run from ions with a known exact mass that lipid spectra always carry (phosphocholine 184.0733 and sphingosine 264.2686 in positive mode; the 16:0, 18:2, 18:1, 20:4 and 22:6 carboxylates in negative mode). For each calibrant the most intense peak within a search window is taken (0.05 Da Orbitrap, 0.10 Da TOF, 0.70 Da ion trap), accepted if it is among the ten most intense peaks and at least 2 % of the base peak. The median error is the offset; the fraction of calibrant peaks within a window, after subtracting the offset, is what that window captures.

## Every study the pipeline has run, 2026-09-12 (three sample files per polarity, blanks excluded)

| instrument | pol | analyzer | calibrant peaks | offset mDa | ≤5 mDa | ≤10 | ≤20 | ≤50 | ≤100 | ≤300 | ≤500 |
|---|---|---|---|---|---|---|---|---|---|---|---|
| Set 1 (Q Exactive) | + | orbitrap | 2,957 | +0.2 | 98 | 98 | 98 | 100 | 100 | 100 | 100 |
| Set 1 (Q Exactive) | - | orbitrap | 1,515 | -0.1 | 98 | 98 | 99 | 100 | 100 | 100 | 100 |
| Set 2 arm 1 (Fusion Lumos) | + | orbitrap | 3,016 | +0.5 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| Set 2 arm 1 (Fusion Lumos) | - | orbitrap | 1,524 | -1.0 | 97 | 98 | 98 | 100 | 100 | 100 | 100 |
| Set 2 arm 2 (Fusion Lumos) | + | orbitrap | 3,605 | +0.4 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| Set 2 arm 2 (Fusion Lumos) | - | orbitrap | 1,807 | -0.9 | 97 | 98 | 98 | 100 | 100 | 100 | 100 |
| Fusion Lumos, ion-trap MS2 | + | ion trap | 11,600 | -24.6 | 6 | 13 | 26 | 57 | 85 | 97 | 98 |
| Fusion Lumos, ion-trap MS2 | - | ion trap | 6,613 | -46.9 | 4 | 8 | 16 | 37 | 56 | 80 | 94 |
| Q Exactive Plus | + | orbitrap | 1,769 | -0.3 | 98 | 98 | 98 | 100 | 100 | 100 | 100 |
| Q Exactive Plus | - | orbitrap | 990 | -0.3 | 100 | 100 | 100 | 100 | 100 | 100 | 100 |
| Q Exactive HF | + | orbitrap | 2,110 | -0.2 | 97 | 97 | 98 | 100 | 100 | 100 | 100 |
| Q Exactive HF | - | orbitrap | 1,350 | +0.3 | 98 | 99 | 99 | 100 | 100 | 100 | 100 |
| Exploris 480 | + | orbitrap | 14,495 | -0.3 | 99 | 99 | 99 | 100 | 100 | 100 | 100 |
| Exploris 480 | - | orbitrap | 14,060 | -0.2 | 97 | 99 | 99 | 100 | 100 | 100 | 100 |
| Agilent 6546 | + | tof | 2,723 | -0.1 | 99 | 100 | 100 | 100 | 100 | 100 | 100 |
| Agilent 6546 | - | tof | 2,411 | +0.1 | 94 | 95 | 96 | 97 | 100 | 100 | 100 |
| Agilent 6545 | + | tof | 260 | -0.3 | 93 | 94 | 95 | 97 | 100 | 100 | 100 |
| Agilent 6530A | + | tof | 436 | -1.1 | 85 | 95 | 98 | 99 | 100 | 100 | 100 |
| Agilent 6530A | - | tof | 108 | -0.1 | 96 | 99 | 99 | 99 | 100 | 100 | 100 |
| Bruker timsTOF Pro (wide search, see text) | + | unknown | 4,269 | +0.3 | 72 | 73 | 73 | 79 | 82 | 88 | 90 |
| Bruker timsTOF Pro (wide search, see text) | - | unknown | 9,107 | -0.6 | 57 | 59 | 61 | 64 | 68 | 78 | 93 |
| Bruker maXis II | + | unknown | 438 | -0.3 | 95 | 97 | 97 | 99 | 99 | 100 | 100 |
| Bruker micrOTOF-Q II | + | unknown | 1,023 | -1.5 | 87 | 90 | 91 | 91 | 98 | 99 | 100 |
| Bruker micrOTOF-Q II | - | unknown | 712 | -6.6 | 23 | 80 | 90 | 91 | 93 | 95 | 98 |
| Sciex TripleTOF 6600 | + | tof | 3,964 | -0.8 | 99 | 100 | 100 | 100 | 100 | 100 | 100 |
| Sciex TripleTOF 6600 | - | tof | 4,067 | -2.3 | 83 | 96 | 99 | 100 | 100 | 100 | 100 |
| Sciex TripleTOF 6600 (high mass) | + | tof | 4,994 | -0.8 | 99 | 100 | 100 | 100 | 100 | 100 | 100 |
| Sciex TripleTOF 6600 (high mass) | - | tof | 4,599 | -2.0 | 81 | 96 | 99 | 100 | 100 | 100 | 100 |
| Sciex X500R | + | tof | 1,934 | +0.3 | 98 | 99 | 99 | 100 | 100 | 100 | 100 |
| Sciex ZenoTOF 7600 | + | tof | 237 | -0.2 | 97 | 100 | 100 | 100 | 100 | 100 | 100 |
| Sciex ZenoTOF (oxTG) | + | tof | 2,424 | -0.3 | 95 | 99 | 99 | 100 | 100 | 100 | 100 |
| Waters Xevo G2 | + | unknown | 171 | +1.2 | 96 | 98 | 98 | 98 | 98 | 99 | 100 |

Re-measured with the Bruker vendor fallback (searched as a TOF, 0.10 Da): timsTOF Pro positive 3,518 peaks, offset +0.3 mDa, captured 88 % at 10 mDa, 89 % at 20, 95 % at 50, 100 % at 100; negative 6,191 peaks, offset −0.2 mDa, 88 / 90 / 94 / 100 %. The timsTOF is the one TOF whose 0.01 Da window captures fewer than 90 % of true fragments; its measured window is 0.1 Da.

## What the numbers decided

- Every Orbitrap and every other TOF captures 95 to 100 % of true fragments within 10 mDa with an offset under 1 mDa: the inherited 0.01 Da stands and nothing on those runs moves.
- The Fusion Lumos with ion-trap MS2 reads 25 mDa low in positive mode and 47 mDa low in negative, and needs 0.3 Da (positive) to 0.5 Da (negative) to capture 97 %; after the offset is subtracted the window narrows. `derive_ms2_tol: auto` replaces 0.01 there.
- Neutral-loss calibrants (59, 141, 183 Da losses) were tried and rejected: they inherit the error of the reported precursor and read 15 to 20 mDa off where absolute-mass ions read within 1 mDa. Cholesterol 369.3516 was rejected after reading +20 mDa on a ZenoTOF.
- Bruker and Waters mzML declare no analyzer term; a vendor name in the header is taken as a TOF.

## The rule

`auto` (default): the configured window is kept when it lies in the analyzer's plausible band (Orbitrap 0.005-0.05 Da, TOF 0.005-0.10, ion trap 0.2-0.7) AND captures at least 90 % of the calibrant peaks; otherwise the measured window is used and the measured offset is subtracted from every fragment m/z before scoring. `true` always applies the measurement, `false` only reports it. The measurement, the decision and the values used are written to `run_config.json` under `measured_ms2` on every run. The window describes the measurement; the dot-product thresholds that decide a hit are never derived (`calibrate.FORBIDDEN`).

Measurement script for a folder of mzML files without running the pipeline: `scripts/ms2_mass_error.py`.

## Proof of concept, 2026-09-12: the ion-trap batch from the default to the measured window

NIST SRM 1950, Fusion Lumos with parallel ion-trap MS2 (6 replicates + 2 blanks per polarity), same run in every respect except the MS2 window (mass calibration second pass on, matched decoys):

| window | identified rows +/− | MS2-confirmed rows +/− | decoy rows +/− | MS2 species (both polarities) | NIST consensus recall (231 in scope) |
|---|---|---|---|---|---|
| 0.01 Da, the inherited default, derivation off | 453 / 64 | 293 / 58 | 1 of 590 / 0 of 78 | 243 | 154 |
| 0.3 Da measured by the run, 24 mDa offset applied (`auto`) | 730 / 193 | 646 / 183 | 2 of 1,020 / 3 of 237 | 440 | 206 |
| 0.3 Da set by hand, no offset (control) | 728 / 194 | 643 / 183 | 2 of 1,021 / 3 of 240 | 440 | 206 |

The run measured 0.3 Da from 11,643 calibrant peaks without being told the spectra came from a trap; against the default it doubles the MS2-confirmed rows and lifts consensus recall from 154 to 206 of 231 at the same decoy rate. The offset correction on top of the right window changes three rows: a 0.3 Da window already covers a 24 mDa shift. On the Bruker timsTOF Pro arm of the corpus (a TOF capturing 88 % of fragments at 10 mDa) the measured 0.1 Da window changed nothing (164 vs 165 identified rows), so the rule is neutral where the default already works and decisive where it does not.
