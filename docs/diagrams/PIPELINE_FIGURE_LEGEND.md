# Pipeline figure — legend and reading notes

The figure itself is `pipeline_raw_to_csv.html`, kept free of prose so it can be used as a
manuscript figure unaltered. This is what belongs in a caption or a supplementary note.

⚠ **No numbers from any one study.** No counts, no timings, no thresholds — a study of a
thousand injections does not behave like one of twenty, and a figure that quotes either is
wrong for the next reader.

## What each decision costs when it is wrong

### Purity decides the resolution of the name, but does not gate it

The purity threshold decides whether a row is reported at molecular or at sum composition. It scores the *winning* identification — it is not a gate on the name that identification carries, and library entries are routinely written at chain resolution. A row can therefore be chain-resolved on no chain evidence at all, which is what the `Chain Evidence` column reports.

### A class with no usable retention model is exempt from the retention check

The model is fitted per class. A class with too few members, or one that is not a homologous series, produces no model — and every row in it passes the retention check **by never taking it**. Those rows carry a blank model error, which reads as though nothing was wrong. `RT Model R2` says whether the class model is worth believing.

### The class retention window is off by default

It removes long-chain species by construction rather than by evidence, so the per-class retention model does that job instead.

### Presence is judged per group, never globally

A row is kept if it was detected in enough of *any one* group. Judged across the whole study instead, a species present in only one arm — often the interesting one — is deleted for being absent from the other.

### Gaps are measured, not imputed

Where the peak finder found nothing, the raw file is re-integrated at that mass and retention time, so a filled value is a measurement at that coordinate rather than a number inferred from other samples. A group with no detections at all keeps its zeros: presence and absence are never erased.

### A class in both polarities with no rule stops the run

The assignment table is checked against every run. A class identified in both polarities but absent from the table halts processing rather than being guessed, because a wrong guess deletes real lipids and nothing downstream shows it.

### The two polarities are not two measurements of one thing

Different ions with different ionisation efficiencies. They are never averaged or merged; the combined matrix carries the mode and scale of each block to say so.

### Conversion is cached

An mzML newer than its raw file is reused, so reprocessing skips the conversion stage entirely. Conversion is written atomically — an interrupted convert must not leave a short file that the freshness check then trusts.

## The two additions that postdate the first version of this figure

### Mass calibration is a loop, not a step

Linking uses a fixed mass tolerance, and a file whose axis has drifted outside it can never join
the consensus group — its features form their own, leaving a zero in the main row for every
injection concerned. Measured across the validation set, within-batch drift ranges from 0.1 ppm to
17.1 ppm, and analyser class predicts nothing: three time-of-flight instruments are tighter than
the Orbitraps.

⚠ **The second pass re-runs the whole pipeline from the search, not just the linking.** Calibration
shifts MS1, MS2 **and precursors** — correcting the spectra while leaving the precursors behind
would push every MS2 off the peak it came from — so once the precursors move, the identifications
move with them.

⚠ And the loop is the point: **the offset cannot be measured until the features have been linked
once.** A first pass must therefore link loosely enough that drifted files are present to be
measured at all; pass two links tight, because the drift has been removed rather than accommodated.

This is why pre-centroiding matters. ST004797 is 199.5 GB of profile data at ~6 h a pass;
centroided once it is 33 GB at about an hour, which is the difference between two passes being
affordable and not.

### Target and decoy are searched together

The decoy count is a rate only because both are searched against the same spectra. Two
constructions appear because one does not fit every class: **chain shuffle** where the score depends
on the chains, **head-group displacement** where it does not.

⚠ A chain shuffle in a class whose scored intensity is chain-independent produces a *copy* of its
own target rather than a decoy — measured shared intensity reaches 0.97 for HexCer[NS]. And a class
with no decoy at all does not have a low error rate; it has an **uncountable** one.

## Which table to use

| file | |
|---|---|
| `Final_Results_Filtered.csv` | the analysis-ready table: filtered, gaps measured. **Start here.** |
| `Final_Results.csv` | what survived the filters, before the presence filter and gap filling |
| `Unfiltered_Results.csv` | every compound group with its `Filter Status` — why a row is absent |
| `…_median_Normalised.csv` | the filtered table, median-normalised |
| `Associated_Spectra.csv` | every spectral match behind a row |
| `Spectral_Components.csv` | how chain-bearing signal in a peak divided between candidates |
