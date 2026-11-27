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
   │  post-filters: adduct pairs, split peaks, co-elution, blank ratio, presence, gap filling
   ▼
Results/{Pos,Neg}/Final_Results_Filtered.csv        analysis-ready, one row per compound group
Results/{Pos,Neg}/Unfiltered_Results.csv            every row, with its filter status
Results/{Pos,Neg}/Associated_Spectra.csv            every spectral match, keyed to its row
Results/Combined_Filtered_Normalised.csv            one polarity per class, chains carried across
Results/QC_report.html                              the batch, its decoy rate, its acquisition
Method_Lists/                                       exclusion + inclusion lists for the next run
```

Spectral scoring is compatible with LipiDex (Hutchins et al. 2018) and reproduces its scores to
the last digit on the same input, so every difference in output is attributable to what LipidLoop
adds.

## Quick start

```bash
git clone https://github.com/jmarquesUOE/LipidLoop.git && cd LipidLoop
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
curl -L -o lipidloop_libraries.zip \
  https://github.com/jmarquesUOE/LipidLoop/releases/download/v1.1.0/lipidloop_libraries_v1.1.0.zip
unzip lipidloop_libraries.zip -d data/libraries/
python scripts/run_study.py /path/to/study --name "my study" --modifier formate \
    --decoy-matched --auto-calibrate
```

A study is a folder with `Pos/` and/or `Neg/` holding one raw file per injection.

## Documentation

| | |
|---|---|
| [docs/USER_GUIDE.md](docs/USER_GUIDE.md) | install, converters, study layout, options, outputs, method lists, library builders |
| [docs/METHOD.md](docs/METHOD.md) | what each stage does, in run order |
| [docs/LIMITATIONS.md](docs/LIMITATIONS.md) | what the error estimate does not cover, and other known limits |
| [docs/LIPIDEX_ALGORITHM.md](docs/LIPIDEX_ALGORITHM.md) | the LipiDex scoring logic, as implemented |
| [CHANGELOG.md](CHANGELOG.md) | changes by release |

## What is in the repository

| folder | content |
|---|---|
| `src/lipidloop/` | the package: conversion, feature detection, search, decoys, peak finder, retention model, co-elution, polarity rule, calibration, artefact screen, method lists, QC |
| `scripts/` | `run_study.py` (the entry point), `run_pipeline.py`, `qc_report.py`, library and decoy builders, calibration, method-list tools, converters |
| `tests/` | 526 tests; `pytest -q` runs in under a minute |
| `docs/` | user guide, method, limitations, pipeline diagram |
| `data/` | class-to-polarity table, adduct and fatty-acid tables, standards, library provenance |

## Validation, in one paragraph

On the same fragment-spectrum files LipiDex read, the identification stage annotates the identical
set of 21,392 spectra with identical scores. It recovers 223 of the 310 NIST SRM 1950
interlaboratory consensus lipids within the libraries' scope with fragment-spectrum evidence, and 12
more through the retention model. On fifteen public studies from five instrument vendors (738
injections searched with one configuration) it agrees with each deposit's published
identifications at 25 to 79 %, at a corpus-wide decoy rate of 0.23 % (23 of 9,951 rows). Against
ground truth (authentic standards, a curated spiked serum and a closed-world spike-in experiment),
the decoy estimate measures chain-assignment error as designed. Separately, 4 to 8 % of delivered
rows were isotopes, adducts or in-source ions of another lipid, an error no spectral decoy
represents and that the co-elution step now recognises. Fed back to the instrument as exclusion,
inclusion or both lists, its output raised spectrum-confirmed lipids by 4, 6 and 9 % in a
two-by-two acquisition experiment. Details are in the manuscript and its Supporting Information.

## Citing

Marques, J. G.; von Kriegsheim, A. LipidLoop: a target-decoy lipidomics pipeline with
retention-time evidence and acquisition feedback. Manuscript submitted (2026). Software: this
repository, release 1.1.0 (`CITATION.cff`, `CHANGELOG.md`).

## Licence

The **code** in this repository is MIT. See `LICENSE`.

The **spectral library archive**, distributed as a release asset, is not all MIT. It bundles
material from three sources under different terms, itemised per file in `NOTICE.md`:

* libraries from LipiDex 1.1, verbatim or derived, under the MIT licence of `coongroup/LipiDex`;
* libraries built from the MS-DIAL Tandem Mass Spectral Atlas under **CC BY 4.0**, which requires
  attribution, a statement of changes and the licence URI;
* libraries generated for LipidLoop, under this repository's MIT licence.

`NOTICE.md` is also included inside the archive, because the MIT licence requires its copyright
notice to travel with every copy and CC BY 4.0 requires the attribution to do the same.

No material from LipiDex 2 is distributed; its repository states no licence. The ganglioside library
is built from published structures instead, by `scripts/build_ganglioside_library.py`.

The two LipiDex tables in `data/lipidex_src/` are MIT, from `coongroup/LipiDex`.

## Funding

J.G. Marques was supported by Tenovus Scotland (grant E24-15). Mass spectrometry was performed at the Mass Spectrometry
Facility of the Institute of Genetics and Cancer, University of Edinburgh; the Orbitrap Fusion Lumos was funded by Wellcome
Trust Multiuser Equipment Grant 208402/Z/17/Z.
