# The pipeline

```
.raw
  │  convert.py      ThermoRawFileParser, self-contained Linux build, no mono
  ▼
.mzML
  ├─ features.py     pyOpenMS: mass traces -> elution peaks -> features -> align -> link
  └─ mgf.py          MS2 scans, read straight out of the same mzML
       │  search.py + purity.py   library search, LipiDex-equivalent scoring and purity
       ▼
  peakfinder.py      join identifications to features, filter, write
       ▼
Final_Results.csv
```

Run it:

    python scripts/run_pipeline.py --raw data/*.raw

`--libraries` defaults to the standard set under `data/libraries/`; pass it only to search
something else.

or with a saved configuration, which is the reproducible form:

    python scripts/run_pipeline.py --config run_config.json

Every run writes its `run_config.json` next to the results. That is not tidiness: **library
load order changes the output** — it breaks ties between entries that score identically, and the
same lipid exists in more than one library with fragment masses differing in the fourth decimal.
LipiDex takes that order from the sequence the libraries were ticked in its GUI and records it
nowhere, which is why an archived LipiDex result cannot be fully reproduced from the file alone.

## Setup

`pip install -r requirements.txt` covers Python. ThermoRawFileParser is a separate download —
the self-contained Linux build from
[compomics/ThermoRawFileParser](https://github.com/compomics/ThermoRawFileParser/releases),
unpacked into `data/tools/ThermoRawFileParser/`. It needs no mono and no .NET runtime.

The libraries and the fatty-acid database come from a LipiDex install: `src/msp_files/*.msp`
and `src/backup/FattyAcids.csv`. The default set, **in load order**, is

    LipidBlast_Formic  LipiDex_HCD_Formic  LipiDex_HCD_Hydroxy  LipiDex2_Ganglioside

The first three are LipiDex 1's; `LipiDex2_Ganglioside.msp` comes from a **LipiDex 2** install
(`SpectralLibraries/`) and is the one library from v2 that adds anything here — see
[LIBRARY_VERSIONS.md](LIBRARY_VERSIONS.md). Append to that list rather than reordering it: load
order breaks ties between entries that score identically. Use the copies from the install that produced any results you
intend to compare against rather than re-downloading, so a difference can never be a library
version.

## Stage notes

**Conversion** is skipped when an `.mzML` already exists and is newer than its `.raw`. About
40 s and 150 MB per file here.

**Feature detection** is the stage with real parameters, and `noise_threshold` is the one that
matters: below about 5,000 on these files, mass trace detection produces hundreds of thousands
of traces and the run does not finish in any useful time. At 5,000 it finds ~5,900 features per
file in about 3 s, against Compound Discoverer's ~5,200.

A caveat this data makes concrete: a DDA run spends its duty cycle on MS2. These files carry
1,838 MS1 scans against 15,019 MS2 — about 1.2 Hz — and with peaks near 0.07 min FWHM that is
roughly five MS1 points across a peak. Thin for any peak picker. It is a property of the
acquisition, not of the software, and `min_trace_length` has to be set against it.

**Identification** is exact against LipiDex; see [VALIDATION.md](VALIDATION.md).

**The peak finder** is where an identification stops being "this spectrum looks like PC 34:1"
and becomes a row with a retention time and an area per sample. The join is on mass, polarity,
sample and retention time, and the retention time is the awkward part: an MS2 is taken whenever
the instrument decided to fragment, while a feature's retention time has been aligned across
files. LipiDex recovers the offset by fitting the deviation between each feature's measured and
aligned retention time, per sample, with a LOESS, and maps the MS2 back through it. That fit is
reproduced in `loess.py` rather than imported, because its bandwidth and robustness-iteration
count are part of the algorithm.

## What a run writes

Delivered folder layout, built by `pipeline.standard_output_dir` and `pipeline.report_dir` so a
new study does not reinvent it and the report cannot drift from where the tables actually are:

    Analysis_<date>/
      Results/Pos/   the files below
      Results/Neg/   the same, other polarity
      Report_QC/     the quality report, both polarities in one document

Inside each polarity folder:

    Final_Results_Filtered.csv the analysis-ready table: IDENTIFIED LIPIDS ONLY, gaps measured
                               by re-integrating the raw MS1 rather than imputed
    Final_Results.csv          one row per compound group, named or not, filters applied, gaps
                               kept as zeros
    Unfiltered_Results.csv     every group including rejected ones, with `Filter Status` giving
                               the reason each was dropped
    Associated_Spectra.csv     one row per spectral match at every rank, keyed to its compound group
    run_config.json            every setting that changes the answer
    run_config.calibrated.json the settings pass 2 ran under, if calibration was enabled
    search/                    the per-injection `_Results.csv` the peak finder reads

`Compound Group` is numbered over the full set and is stable across all four tables, so a row in
any of them joins to the same group in the others.

`Final_Results.csv` always keeps its zeros. The filled table is a second file, never a
replacement, so nothing downstream is handed filled values it did not ask for and the two can be
diffed to see exactly what was filled. The gaps are **measured**, not imputed — see
[GAP_FILLING.md](GAP_FILLING.md). See [SPLIT_PEAKS.md](SPLIT_PEAKS.md) for rows that are one peak
twice, [ADDUCT_PAIRS.md](ADDUCT_PAIRS.md) for rows that are one molecule seen as two ions, and
[PRESENCE.md](PRESENCE.md) for which rows are analysable at all, and
[CALIBRATION.md](CALIBRATION.md) for deriving a run's parameters from the run.

The per-injection files are an intermediate, not a deliverable — at fifty-odd injections a
polarity they buried the four files a reader wants, so they live in `search/`. They are still
written, unchanged, because the peak finder reads them back and `validate_qc.py` compares them
against real LipiDex output.

## What the output columns mean

`Final_Results.csv` is one row per compound group — one molecule, aligned across samples.

| column | meaning |
|---|---|
| Compound Group | stable id for the row; the same number in `Unfiltered_Results.csv` and `Associated_Spectra.csv` |
| Retention Time (min) | aligned retention of the group |
| Quant Ion | the most abundant ion in the group; the one quantified |
| Polarity | polarity of that ion |
| Area (max) | largest area across samples |
| Identification | blank if unidentified; sum composition if purity < 75; molecular species otherwise |
| Lipid Class | class of the top candidate |
| Adduct | which ion the row was quantified on; several, quantified one first, if the molecule was seen as more than one |
| Features Found | how many samples the group was detected in |
| one column per sample | area in that sample |

`Adduct` exists because `Lipid.java :: toString` drops it — the reported identity is the
molecule, not the ion it was seen as. That is right for the name and leaves the table unable to
answer two ordinary questions: which ion `Quant Ion` measures, and whether two rows at one
retention time are two adducts of one molecule. The second is not hypothetical — `Na - H` and
`+2 carbons, +3 double bonds` differ by **2.4 mDa, 3.0 ppm at m/z 800**, so the sodium adduct of
`PE 40:1` sits inside a normal search window of protonated `PE 42:4`.

`Compound Group` numbers every group in the unfiltered set, so a row dropped from
`Final_Results.csv` keeps its number and the tables join exactly. Without it the only key is
name plus retention time, which is precisely what is ambiguous when one molecule appears on
several rows.

### `Associated_Spectra.csv`

One row per spectral match across the whole run, replacing the job of opening fifty per-injection
files and joining them by hand. Columns: `Compound Group`, `Name` (with adduct), `Adduct`,
`Sample`, `MS2 ID`, `Rank`, `Precursor`, `Library Mass`, `Delta m/z (ppm)`, `Retention`,
`Polarity`, `Dot Product`, `Reverse Dot Product`, `Purity`, `Library`, `Associated`.

Two departures from LipiDex's file of the same name. It carries `Compound Group`, so the join
back to the results table is exact rather than a match on name and retention time. And it keeps
matches that were **not** associated with any group — `Associated` is blank for those — because a
table showing only the winners cannot be used to ask why something is missing. On the Kiterie
positive run: 83,350 matches over 53 injections, 35.7% associated, referencing 1,287 groups.

`GaussianScore` is absent. The field exists on `Lipid` because Compound Discoverer supplies it
and this pipeline has no equivalent; a column of zeros would read as a measurement rather than an
absence.

The identification column reports at the resolution the evidence supports. A row reading
`PC 34:1` is not a failure to determine the chains — it is a statement that the fatty-acid
purity did not clear 75%, so the chains are not asserted. `PE O-38:6` rather than a plasmenyl
or plasmanyl name means the evidence split between the two and neither was asserted; that
matters if you are reading plasmalogens off this column.
