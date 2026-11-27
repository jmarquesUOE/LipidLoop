# BileAcids_Negative*.msp — where these came from and what was changed

**Source.** MS-DIAL Tandem Mass Spectral Atlas VS69 (negative), Zenodo record 10953284, entries
whose `COMPOUNDCLASS` is `BileAcid` or `BASulfate`.

**Licence.** CC-BY 4.0. Attribution: Tsugawa H. *et al.*, *A lipidome atlas in MS-DIAL 4*, Nature
Biotechnology 2020, doi 10.1038/s41587-020-0531-2. Redistribution and modification are permitted
with attribution and a statement of changes; this file is that statement.

**Changes made** (by `scripts/import_msdial_class.py`, which is the reproducible record):

1. **Filtered** to the two bile-acid classes — 47 entries of 792,757.
2. **Relabelled `ST ` to `BA `.** MS-DIAL names bile acids with the same `ST` prefix it uses for
   sterols and separates them only in `COMPOUNDCLASS`. This pipeline takes the class from the first
   token of the name, so importing verbatim would put bile acids and sterols on ONE retention
   surface. They elute nothing like each other; the fit would be junk and would take the sterols
   with it. `BA` is what MS-DIAL's own result tables use.
3. **Dropped 11 entries with fewer than 3 peaks.** A two-peak theoretical spectrum whose second
   peak is the precursor cannot discriminate anything, and bile acids are the class where that
   matters most — see below.
4. **Split by adduct.** `[M-H]-` (25) forms whatever the modifier and goes in both arms;
   `[M+CH3COO]-` (11) exists only in an acetate mobile phase and is loaded for the acetate arm
   alone. Shipping a modifier-specific adduct to the wrong run is how a library entry gets matched
   to the wrong molecule.

## ⚠ What this library cannot do

These are **theoretical** spectra — the Comment field says "created from the information of SCIEX
5600 and 6600" — not measured ones.

**31 of the original 47 entries share a precursor mass with another entry.** Bile acids are largely
stereoisomers differing in hydroxyl orientation, not in chain length or unsaturation, so mass does
not separate them and neither does a fragmentation model that does not know stereochemistry.
**Retention is the only thing that can.**

That is what the 96-well authentic bile-acid plate is for: 47 of its 96 compounds sit in isobaric
groups, six at 408.2876 alone. Until those are injected and their retentions measured, an
identification from this library should be read as "a bile acid of this composition", not as a
named stereoisomer.
