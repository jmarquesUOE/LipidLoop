# UltraLongCer_EOS / UltraLongHexCer_EOS / AcylSM / AcylHexCer — provenance

**Source.** MS-DIAL Tandem Mass Spectral Atlas VS69 (positive), Zenodo record 10953284,
`COMPOUNDCLASS` in {`Cer_EOS`, `HexCer_EOS`, `ASM`, `AHexCer`}. 28,154 entries.

**Licence.** CC-BY 4.0. Tsugawa H. *et al.*, *A lipidome atlas in MS-DIAL 4*, Nature Biotechnology
2020, doi 10.1038/s41587-020-0531-2. `scripts/import_msdial_class.py` is the statement of changes.

## Why these four

No library shipped here contained a single omega-esterified ceramide, acyl sphingomyelin or acyl
hexosylceramide — **zero entries across all four classes**. These are the skin-barrier lipids, and
the gap is one of chain length as much as class: our libraries stop at 26 carbons on a single
chain, while these reach **64**.

| class | entries | longest single chain | past 26 C |
|---|---|---|---|
| Cer[EOS] | 3,179 | 26–**64** | 98.8% |
| HexCer[EOS] | 3,759 | 28–**64** | 100% |
| ASM | 10,416 | 30–48 | 100% |
| AHexCer | 10,800 | 16–30 | 41% |

## Changes made

1. **Filtered** to the four classes, 28,154 of 1.06 M atlas entries.
2. **Relabelled the class into the name.** MS-DIAL puts the subclass only in `COMPOUNDCLASS`: its
   `Cer_EOS` entries are *named* `Cer …`, and `ASM` entries are named `SM …`. This pipeline takes
   the class from the first token, so importing verbatim would merge ultra-long esterified
   ceramides into ordinary `Cer` and acyl-SM into ordinary `SM` — sharing one retention surface
   with molecules 30 carbons shorter. Written as `Cer[EOS]`, `HexCer[EOS]` and `ASM`.
3. **Rewrote three notation differences**, each of which silently destroyed the sum composition:
   - `;2O` → `;O2` (oxygen count, written the other way round)
   - `/` → `_` between chains. MS-DIAL's slash ASSERTS sn-position here and, worse, was not read
     as a chain separator at all: `Cer[EOS] 14:1;O2/26:1;O2` summed to **`Cer[EOS] 14:1`** instead
     of 40:2 — silently dropping the second chain, and with it the ultra-long tail that is the
     entire reason for importing the class.
   - the parenthesised esterified chain, `(O-14:0)` leading and `(FA 14:0)` trailing, which no
     parser here reads. ⚠ It appears at BOTH ends, so a one-sided fix works for AHexCer and fuses
     ASM into `;O214:0`, which parses as no chain at all.

   After these, **100% of entries yield a parseable sum composition and retention class**, against
   0% for ASM and AHexCer before.

## ⚠ What this library cannot do

Theoretical spectra, positive mode `[M+H]+` only — the source ships no negative-mode entries for
these classes, so they are unreachable in a negative run whatever the chemistry says.

**Not yet validated on real data.** The bile-acid import is the cautionary precedent: 36 entries,
correctly loaded, produced zero identifications because the analytes were not retained by the
chromatography. Having a library is not evidence that the class is measurable. These should be run
against **skin** data, where the classes are abundant and structurally important, before any claim
is made — rat plasma annotates 4 AHexCer and essentially no EOS ceramides, so it cannot test them.
