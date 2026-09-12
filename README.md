# LipidLoop

An open, command-line pipeline for untargeted LC-MS/MS lipidomics: from vendor raw files to an
analysis-ready lipid table, with a target-decoy error rate on every batch, a per-class retention
model whose calls are labelled, one polarity per class with chain evidence carried across, per-file
mass calibration, and a diagnosis of the acquisition that is fed back to the instrument as exclusion
and inclusion lists. Every setting that changes the answer is written beside the result.

```
vendor raw (.raw / .d / .wiff)
   │  ThermoRawFileParser, msconvert or alphatims, chosen from the file itself
   ▼
mzML (cached once per study)
   │  pyOpenMS feature detection, alignment, linking; two-pass mass calibration
   │  MS2 search against the target AND decoy libraries; artefact screen
   │  peak finder: join identifications to features, retention model, filters
   │  post-filters: adduct pairs, split peaks, blank ratio, presence, gap filling
   ▼
Results/{Pos,Neg}/Final_Results_Filtered.csv        analysis-ready, one row per compound group
Results/{Pos,Neg}/Unfiltered_Results.csv            every row, with its filter status
Results/{Pos,Neg}/Associated_Spectra.csv            every spectral match, keyed to its row
Results/Combined_Filtered_Normalised.csv            one polarity per class, chains carried across
Results/QC_report.html                              the batch, its decoy rate, its acquisition
Method_Lists/                                       exclusion + inclusion lists for the next run
```

Spectral scoring is compatible with LipiDex (Hutchins et al. 2018) and reproduces it to the digit on
the same input, so every difference in output is attributable to what LipidLoop adds. The
algorithm, the constants and the reproduced filters are documented in
[docs/LIPIDEX_ALGORITHM.md](docs/LIPIDEX_ALGORITHM.md); the validation in
[docs/VALIDATION.md](docs/VALIDATION.md).

## Install

Python 3.11 or later, Linux or macOS (Windows works for Thermo files; other vendors need msconvert).

```bash
git clone https://github.com/jmarquesUOE/LipidLoop.git
cd LipidLoop
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

**Spectral libraries** are not in the tree. Download the archive attached to the release and unpack
it into `data/libraries/`:

```bash
curl -L -o lipidloop_libraries.zip \
  https://github.com/jmarquesUOE/LipidLoop/releases/download/v1.0.0/lipidloop_libraries_v1.0.0.zip
unzip lipidloop_libraries.zip -d data/libraries/
```

The archive holds the LipiDex and LipidBlast libraries (MIT), the libraries generated here
(`data/libraries/PROVENANCE_*.md` records where each came from and what was changed) and the
matching decoy libraries. Decoys can be regenerated with `scripts/make_decoy_library.py`.

**Converters.** Thermo `.raw` needs [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser)
(put it in `data/tools/ThermoRawFileParser/` or pass `--set thermo_parser=<path>`); Agilent, Waters
and Sciex need msconvert (ProteoWizard); Bruker `.d` is read natively.

## Run a study

A study is a folder with `Pos/` and/or `Neg/` holding the raw files, one file per injection, blanks
and QCs named as such.

```bash
python scripts/run_study.py /path/to/study --name "my study" --modifier formate \
    --decoy-matched --auto-calibrate
```

`--modifier` selects the library set for the mobile-phase modifier (formate or acetate);
`--decoy-matched` searches the matched decoy libraries and reports the error rate;
`--auto-calibrate` runs the two-pass per-file mass calibration. Every setting can be overridden with
`--set key=value`; the full configuration used is written to `Analysis_*/config_{Pos,Neg}.json`.
`python scripts/run_study.py --help` lists everything.

## What is in the repository

| folder | content |
|---|---|
| `src/lipidloop/` | the package: conversion, feature detection, search, decoys, peak finder, retention model, polarity rule, calibration, artefact screen, duty-cycle diagnosis, method lists, QC report |
| `scripts/` | `run_study.py` (the entry point), library builders, decoy generation, method-list export, calibration, the validation and audit tools |
| `tests/` | 425 tests; `pytest -m "not integration and not slow"` runs in about a minute |
| `docs/` | one document per stage, each stating what the stage does, what it was tested against and what it gets wrong |
| `data/` | the class-to-polarity table, adduct and fatty-acid tables, library provenance |
| `paper/` | every figure and benchmark of the manuscript, regenerated from the scripts and data in the folder |

## Validation, in one paragraph

On 16,435 MS2 spectra the identification stage returns the same identified spectra as LipiDex with
identical scores. End to end it recovers 99.4 to 99.5 % of the reference workflow's molecules in
positive mode and 92.8 to 94.7 % in negative mode on two Orbitrap instruments, and 81 to 86 % of
the NIST SRM 1950 interlaboratory consensus lipids within its scope. On fifteen public studies from
five instrument vendors, 738 injections searched unmodified, it agrees with each deposit's published
identifications at 25 to 79 % at a corpus-wide decoy rate of 0.13 %. Fed back to the instrument, its
duty-cycle diagnosis cut the wasted MS2 budget from 40 % and 63 % to 15 % and 23 % on the same
sample. Details and the data behind each number are in `paper/`.

## Citing

Marques, J. G.; von Kriegsheim, A. LipidLoop: a target-decoy lipidomics pipeline with
retention-time evidence and acquisition feedback. Manuscript in preparation (2026). Software: this repository, release 1.0.0 (`CITATION.cff`).

## Licence

MIT (see `LICENSE`). The LipiDex and LipidBlast spectral libraries are redistributed under their own
MIT licences; the two LipiDex tables in `data/lipidex_src/` likewise.
