# User guide

## Install

Python 3.11 or later, on Linux or macOS.

```bash
git clone https://github.com/jmarquesUOE/LipidLoop.git && cd LipidLoop
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

Install from the checkout (`-e`): the pipeline reads its tables from `data/` next to the source.

**Spectral libraries** are a release asset, not part of the tree:

```bash
curl -L -o lipidloop_libraries.zip \
  https://github.com/jmarquesUOE/LipidLoop/releases/download/v1.1.0/lipidloop_libraries_v1.1.0.zip
unzip lipidloop_libraries.zip -d data/libraries/
```

**Converters.**
- Thermo `.raw`: [ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser), the
  self-contained Linux build, in `data/tools/ThermoRawFileParser/` (or `--set thermo_parser=<path>`).
- Agilent, Waters and Sciex: msconvert (ProteoWizard).
- Bruker `.d`: read natively; `scripts/tdf_to_mzml.py` converts a timsTOF file on its own.
- mzML is accepted as is. Profile mzML can be centroided first with `scripts/precentroid.py`.

Run the tests with `pytest -q` (under a minute; tests that need the library archive skip until it
is unpacked).

## A study folder

```
my_study/
  Pos/   one raw file per injection, positive mode
  Neg/   one raw file per injection, negative mode
  sequence.csv   optional: injection order and Sample Type, as exported from the instrument
```

Blanks and pooled QCs are recognised from the sequence or from the file name (`Blank`, `QC`,
`PooledQC`, an acquisition-mode tag such as `QCMS01`). Either polarity can be absent.

## Run it

```bash
python scripts/run_study.py /path/to/my_study --name "my study" --modifier formate \
    --decoy-matched --auto-calibrate
```

| option | what it does |
|---|---|
| `--modifier formate\|acetate` | picks the library set for the mobile-phase modifier (negative-mode adducts follow it) |
| `--decoy-matched` | searches the decoy libraries built from the same targets and reports the error estimate |
| `--auto-calibrate` | two-pass mass calibration: measure each file's mass offset, correct, detect again |
| `--metadata-pos/--metadata-neg FILE`, `--group-by COLUMN` | biological groups for the presence filter; without them the samples form one group (said in the log and the report), and a template is written for you to fill in |
| `--set key=value` | overrides any setting; the override is written into the configuration |
| `--no-polarity-filter` | delivers both polarities instead of one per class |
| `--no-coelution-filter` | switches off the co-elution step |
| `--class-polarity FILE` | your own class-to-polarity table instead of `data/class_polarity.csv` |
| `--analysis-dir DIR` | where the results go (default `<study>/Analysis_<date>/`) |

`python scripts/run_study.py --help` lists everything. `scripts/run_pipeline.py` runs one
polarity from one configuration file, for scripted use.

## What you get

```
Analysis_<date>/
  config_Pos.json, config_Neg.json      every setting that changes the answer, configured or measured,
                                        and the commit of the code that ran
  run_Pos.log, run_Neg.log
  Results/Pos/Final_Results_Filtered.csv    the delivered table, one row per compound group
  Results/Pos/Unfiltered_Results.csv        every row, with the filter that excluded it
  Results/Pos/Associated_Spectra.csv        every spectral match, keyed to its row
  Results/Neg/...                           the same for negative mode
  Results/Combined_Filtered_Normalised.csv  one polarity per class, chain evidence carried across
  Results/QC_report.html / .pdf             the batch, its error estimate, its acquisition
  Method_Lists/                             exclusion and inclusion lists for the next acquisition
  coelution_filter.json, polarity_filter.json
```

Columns worth knowing in the delivered table:

- `Identification`, `Lipid Class`, `Adduct`, `Quant Ion`: the name, the class, and the ion used
  for the area.
- `Identification Source`: MS2, or the retention model without a spectrum (MSI level 3).
- `Dot Product`, `Reverse Dot Product`, `Purity`, `MS2 Spectra`: the evidence behind the name.
- `Chains Resolved`, `Chains From`, `Chain Evidence`: whether the acyl chains were resolved, and
  from which polarity.
- `RT Model Error`, `RT Model Z`, `RT Model Window`: the retention-model verdict, including whether
  the row lies outside the range the class model was fitted on.
- `Duplicate Name`, `Duplicate Verdict`: rows sharing a name, and what was decided about them.
- decoy rows keep their `DECOY_` name and stay in every table, so the estimate can be recomputed.

Other files: `Spectral_Components.csv` (how the chain-bearing signal of mixed spectra divides),
`coelution_filter.json` and `polarity_filter.json` (what those steps removed and why),
`Polarity_Pairing.csv` (which positive and negative injections are the same vial).

## The error estimate

The QC report gives, per batch and polarity, the fraction of spectrum-supported rows assigned to
decoys. It estimates **chain-assignment error within a correct class and sum composition**. It
does not cover rows that are another lipid's isotope, adduct or in-source ion, or rows named by the
retention model without a spectrum. See [METHOD.md](METHOD.md) and [LIMITATIONS.md](LIMITATIONS.md).

## Feeding the results back to the instrument

`Method_Lists/` holds an exclusion list (contaminant-like ions that used MS2 without ever being
identified) and an inclusion list (lipids identified only weakly or only by retention), in the
instrument's import format. `scripts/build_method_lists.py` rebuilds them from an existing analysis;
`scripts/check_exclusion_safety.py` checks an exclusion list against the lipids it could hide.

## Building libraries

Each `scripts/build_*_library.py` regenerates one library of the archive from formulae and
measured spectra; `scripts/make_decoy_library.py` builds the matching decoy library. The
provenance of every archived library is in `data/libraries/PROVENANCE_*.md` and its licence in
`NOTICE.md`.
