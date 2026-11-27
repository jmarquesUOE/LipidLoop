# Changelog

## 1.1.0 (2026-11-27)

The release the revised manuscript was produced with: its reprocessed results come from this code
(commits e42a1dc to 6be2d31 of the development history; later commits change only the combined
table's polarity pairing and provenance) and the library archive `lipidloop_libraries_v1.1.0.zip`
attached to the release. Item numbers in code comments ("v1.1 item N") refer to the list below.

### Changes made after the 1.0.0 snapshot, before the revision

These changes were made between 1.0.0 and the start of the revision (September to October 2026)
and were not in a public release until now:

* **Polarity filter applies to decoys.** A decoy is filtered by the class it imitates. Before this,
  every decoy skipped the class-polarity rule and stayed in a polarity where its real counterparts
  had been deleted, which inflated the delivered decoy rate.
* **MS2 axis.** The fragment window and the fragment-axis offset are separate decisions
  (`derive_ms2_tol = auto | true | false | offset`). An unrecognised value now raises an error
  instead of being treated as `true`.
* **Chain evidence.** A chain-resolved name has to win the feature-level vote of the spectra
  behind it. This is a no-op unless `chain_evidence` is on, which it is not by default.
* **Default libraries.** The positive set no longer names the LipiDex 2 ganglioside library and
  its decoy (1.0.0's defaults named them although its archive did not contain them). The
  1-deoxyceramide library (saturated only) is added to the positive set. In the negative set the
  extended ganglioside library replaces the standard one; it is a superset that adds GM4 and the
  all-NeuAc forms of GM1, GD2 and GD3.
* **Nomenclature.** Coenzyme Q spellings map to `CoQ9`, `CoQ10` and `CoQ10H2`.
* **Class-polarity table.** New rows for GM3-NGNA, PI[OH] and the prenol quinones.

### New in 1.1.0 (changes made during the revision)

* **Co-elution step** (`coelution.py`; items 1 and 2). Removes a delivered lipid that is an isotope, adduct or
  in-source ion of another lipid at the same retention time, even when both rows are named. This
  combines an ion-family pass with a parent check on peak shape, apex and ratio. The step includes
  the alkylammonium adduct rule.
* **Split peaks** (item 3). A merge is refused when the two halves name different molecular species.
* **Pooled QC** (item 4). A pooled QC is recognised when its name carries an acquisition-mode tag, e.g. `QCMS01`.
* **Retention model** (items 5 and 8). Verdicts made outside the class's training window are flagged. A strong
  outlier on weak MS2 evidence (few spectra) is removed even inside the 1 min floor.
* **Polarity rule** (item 7). Ether species seen in both polarities are paired.
* **Decoys and isomers** (items 12 to 14). Decoy names are kept in `top_ether`: an ether decoy
  used to lose its `DECOY_` prefix there and could be delivered as a lipid. In negative mode a
  PE-NMe or PE-NMe2 top hit with PE head-group ions and none of its own gives way to the best PE.
  PE-NMe gets no retention-time exemption in the co-elution cross-polarity test.
* **Ether names** (item 18). A demoted plasmalogen name keeps its formula: `P-x:y` becomes `O-x:(y+1)`.
* **Combined table** (items 19 and 20). Pos and Neg injections are paired from the data, not from the file name. A
  single-letter polarity token (`_P` / `_N`) is accepted.
* **Provenance** (item 21). Every `config_*.json`, `run_config.json` and the QC report record the commit of
  the code that ran (`+dirty` if it was modified).
* **Libraries** (item 15 and new libraries). `LipiDex_Splash_ISTD_{Formic,Acetate}.msp`: the PS d7 head-group loss
  intensity is corrected. New in the archive: `DeoxyCeramides_Positive` and its decoy,
  `Ganglioside_Negative_SumComposition_Extended`, and `Ubiquinones_Positive` with its retention
  anchors and provenance. The archive now holds 63 files.
* New library builders: `build_deoxyceramide_library.py`, `build_hex2cer_library.py`,
  `build_ubiquinone_library.py`, plus `check_ubiquinol_rt.py`.

### Repository

The repository now holds the code, tests, user documentation and data tables. The development
notes, the validation and audit scripts and the figure code of the first submission are no longer
distributed.

### Known and not changed

Feature finding is not fully deterministic: reruns of the same commit can differ by one or two
rows. Rows in a class without a retention model pass unjudged.

## 1.0.0 (2026-09-17)

First public release, cited by the manuscript as first submitted. Archived at
doi 10.5281/zenodo.22813962; superseded by 1.1.0.
