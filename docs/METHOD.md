# Method

What each stage does, in the order a run executes it. The manuscript's Experimental Section is the
citable description; this page follows it. The pipeline figure is in [diagrams/](diagrams/).

## 1. Conversion

The reader is chosen from the file: ThermoRawFileParser for Thermo `.raw`, msconvert for Agilent,
Waters and Sciex, native reading for Bruker `.d`. One mzML is written per raw file and cached for
the study, so a reprocessing does not convert again.

## 2. Annotation (per injection)

Each MS2 spectrum is first screened for chemically implausible artefact ions, then scored against
the target and decoy libraries together. Scoring is compatible with LipiDex 1
([LIPIDEX_ALGORITHM.md](LIPIDEX_ALGORITHM.md)): MS1 and MS2 tolerances of 0.01 Da, dot product and
reverse dot product against the library spectrum, spectral purity from fatty-acid fragment
intensities with a 75 % threshold for a chain-level name, and the same tie-breaking, including
library load order, which is recorded. Names are hierarchical: class from head-group evidence, sum
composition from the precursor, acyl chains only when a chain-specific fragment is measured;
sn-positions are not reported.

The MS2 mass axis is measured on every run from fragment ions of known mass. When the configured
fragment window is inconsistent with the analyser that produced the spectra, the window and its
offset are re-derived; the decision is written to the configuration.

## 3. Decoys and the error estimate

For each target library a decoy library holds one decoy per target entry. The precursor mass,
head-group fragments and intensities are kept and each chain-specific fragment is shifted by whole
methylene units, so a decoy hit is a chain assignment that cannot coexist with its precursor: a
same-class chain-assignment error. Entries without chain-specific fragments (13 classes, about 11 %
of positive-mode entries) use a second construction that displaces every fragment by an implausible
mass. Targets and decoys compete in one search, one to one, so no multiplier is applied.

The estimate is the number of spectrum-supported rows assigned to decoys divided by all
spectrum-supported rows, per batch and polarity, with a binomial interval when the decoy count is
small. Decoy rows pass through every filter with their names visible. The retention model never
judges a decoy, which makes the post-filter estimate marginally conservative.

## 4. Feature detection, alignment and calibration

pyOpenMS detects, aligns and links features (noise threshold 5,000 counts, S/N 3, trace length at
least 3 s, peak width 1 to 20 s, alignment 30 s and 10 ppm). The MS1 noise floor is measured per
run; a configured threshold outside 2 to 15 times the floor is re-derived from it, bounded so that
no more than 600 peaks per scan remain and at least 200 features survive in the richest injection.
With `--auto-calibrate`, each file's mass offset is measured, the file corrected and detection
repeated.

## 5. Peak finder

Annotations are joined to features and the filters of the reference workflow applied: sample to
blank ratio 3, presence in at least 60 % of one biological group of at least three samples. Zero
areas are re-integrated from the mzML within 10 ppm and 0.2 min. No row is deleted from the audit
output: a row that fails a filter stays in `Unfiltered_Results.csv` with the name of the filter.

## 6. Retention model

Within each class, retention time is fitted to total carbon and double-bond number by least squares,
with iterative trimming beyond three median absolute deviations. A model needs at least six points,
R2 at least 0.7 and a residual SD of at most 1 min; otherwise the class is declined and logged. A row
is flagged when it deviates by more than the larger of 3 class SD and 1 min, and marked as an
extrapolation when it lies outside the class model's training range. A strong outlier supported by
few spectra is removed even inside the 1 min floor. An MS1 feature without a spectrum is named only
when exactly one library mass within tolerance lies on its class model (MSI level 3).

## 7. Rows that may be one lipid

- **Adduct pairs**: two named groups within one peak width whose ions are two adducts of one
  neutral mass (10 ppm) collapse only if they track across samples (Spearman rho >= 0.8); the
  class-canonical adduct is kept and areas are never summed.
- **Split peaks**: two groups with the same name and adduct within two peak widths and 10 ppm are
  summed only when the summed row's CV across pooled QCs improves 1.5-fold on the worse of the two;
  never when their chain assignments differ; no action without pooled QCs.
- **Co-elution**: each named lipid predicts the isotopes, adducts, in-source and other-polarity ions
  it should produce; a delivered row explained by a stronger co-eluting feature (same XIC shape,
  r >= 0.9; same apex within 0.05 min; ratio fitting an adduct or the isotope abundance) is
  recognised as that ion. Alkylammonium adducts (+28.0313 Da from [M+NH4]+) are labelled when the
  run shows the series.
- Other rows sharing a name stay as unresolved candidates, flagged when the retention model
  disagrees.

## 8. One polarity per class

Each class annotated in both polarities is reported from the polarity favoured by its head-group
chemistry (`data/class_polarity.csv`; choline and sphingoid classes from positive, acidic
glycerophospholipids, gangliosides and lyso-PE from negative). A class seen in one polarity passes
untouched. Chain-level evidence from the dropped polarity is carried to the kept row in its own
column. Positive and negative injections of the same vial are paired from the data for the combined
table. `--no-polarity-filter` delivers both polarities.

## 9. Reporting and acquisition feedback

The QC report describes the batch: noise floor, mass and MS2 axes, decoy estimate, retention models,
filter losses, and how the MS2 budget was spent. Exclusion lists (ions fragmented repeatedly, never
identified, present in blanks, eluting across much of the gradient) and inclusion lists (lipids
identified weakly or by retention only) are written in the instrument's import format for the next
acquisition.

## Reproducibility

Every run writes `config_Pos.json` / `config_Neg.json` with every setting that changes the answer,
configured or measured, and the commit of the code that ran, and a per-spectrum table linking every
match to the row it annotated.
