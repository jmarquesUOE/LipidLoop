# Di- and trihexosylceramides (Hex2Cer, Hex3Cer), negative mode — provenance

Built by `scripts/build_hex2cer_library.py`, which is the reproducible record.

    python scripts/build_hex2cer_library.py --out data/libraries/Hex2Cer_Negative_Formate.msp
    python scripts/build_hex2cer_library.py --hexoses 3 --out data/libraries/Hex3Cer_Negative_Formate.msp

**176 target entries each, ceramide sum compositions C24–C45, all `[M+HCOO]-`. No decoys** — see
below.

⚠ **Hex3Cer's intensity model is TRANSFERRED from Hex2Cer, not measured.** No Hex3Cer spectrum was
available to measure — the class is absent from the data that prompted it (see "First use"). That
weakens any positive identification, which would need manual spectral confirmation before being
believed. It does not weaken a negative: a zero-match result does not depend on the intensity model.

## Why it exists

Hex2Cer is the direct product of B4GALT5, the lactosylceramide synthase, and the branch point of
the glycosphingolipid tree: globosides and the whole ganglioside series are built on it. Nothing in
the shipped libraries names it. LipidBlast and the LipiDex HCD libraries carry monohexosylceramide
only, and `Ganglioside_Negative_SumComposition` starts at GM4/GM3, one step further down. A study
could therefore measure the substrate (HexCer) and the products two steps on (GM3, GM2, GM1) and
pass silently over the enzyme's own product.

It was found by the gap it left. Four unnamed features in a facility ceramide-inhibitor study
(`20251017_Cer_Inhibitor`) sat within 3 ppm of Hex2Cer formate adducts, at plausible retention
times, with MS2 acquired and no library entry to match it against.

## The adduct, and why it is not `[M-H]-`

A neutral glycosphingolipid has no acidic proton to lose. In a formate mobile phase it is carried
as the formate adduct, and every species confirmed here was seen that way. Each fragments first by
losing formic acid to give `[M-H]-` as the base peak. Writing `[M-H]-` entries would put the
precursor 45 Da away from anything the instrument selected, and the library would match nothing.

⚠ **This ties the library to the mobile phase.** It is correct for this facility, which runs
ammonium formate with formic acid. An acetate method needs `[M+CH3COO]-` at +59.0133 instead, and
this file must not be used against one — the failure would be silent. This is the same trap that
named 19 acetate adducts of ordinary ceramides as phytoceramides on ST003514.

## How the identification is justified

Masses are computed from the elemental formula: `Cer dC:DB` is C(C)H(2C+1-2DB)NO3, plus two
condensed hexose residues at 162.052824 each. Against the four features that prompted the library,
the computed `[M+HCOO]-` agrees to **−2.6 to −1.9 ppm**, which is this run's own global mass bias in
that window rather than an error in the library.

Three ions carry the composition, and the second loss is the whole identification:

    [M-H]-          loss of formic acid from the precursor — the base peak
    [M-H-Hex]-      loss of one hexose
    [M-H-2Hex]-     loss of both, landing on the ceramide anion

One hexose loss is a monohexosylceramide. **Two sequential losses of 162.0528 arriving at the
calculated ceramide anion is what makes this dihexosyl**, and it was observed directly — d34:1 runs
860.61 → 698.55 → 536.50 against a calculated Cer d34:1 anion of 536.505, and d42:2 runs
970.72 → 808.66 → 646.62 against 646.622.

Fragment intensities are medians relative to the `[M-H]-` base peak, measured over every spectrum of
the four confirmed species across all 39 sample injections (d34:1 and d42:2 contribute 39 spectra
each, d40:1 21, d42:1 16). Only four of the eight hexose ions are written: that series is intense
but says only "a sugar was here", appears in any glycosphingolipid spectrum, and weighting it
heavily would let a Hex2Cer entry score against a HexCer or a ganglioside.

## ⚠ No decoy, so no false-discovery measurement

`make_decoy_library.py` builds a decoy by shifting chain-specific fragments, which requires a
fragment annotated with a chain. Every ion here is head-group or precursor derived, annotated `[]`,
so the method produces nothing. Hex2Cer identifications carry the same caveat the gangliosides do:
their false-discovery rate is a ceiling inherited from the rest of the run, not a measurement of
this class.

## ⚠ What the name withholds

`Hex2Cer`, not `LacCer`. Lactosylceramide is only one dihexosylceramide — Gal-Glc-Cer, Gal-Gal-Cer
and Glc-Glc-Cer are isomers that reversed-phase chromatography does not separate and whose spectra
here are identical. `LacCer` would assert a linkage the evidence cannot support. `nomenclature.py`
already maps `laccer` and `lactosylceramide` onto `Hex2Cer`, so either dialect reads correctly.

Sum composition, and no subclass tag: the base and acyl chains are not resolved, so `[NS]` and
`[NDS]` cannot be told apart. The argument is the one made for the gangliosides in
the developer notes (GANGLIOSIDES, not distributed).

## First use

`20251017_Cer_Inhibitor`, re-searched 2026-09-15 into `Analysis_2026-09-15_Hex2Cer`. Seven species
identified in negative mode, dot products 958–999 against reverse 1000. Adding the library changed
**nothing else**: positive stayed at 426 identifications with no class altered, and negative went
247 → 254, the seven Hex2Cer being the only difference.

⚠ One of the seven is doubtful. `Hex2Cer d42:1` appears twice, at RT 8.42 and 8.59. The 8.59 row is
sound (−1.9 ppm, dot 997, 11 spectra); the 8.42 row is −6.1 ppm on 2 spectra and co-elutes with
d40:1. The pipeline flagged the pair `isomer?` and did not merge them, because that study has no
pooled QC injections and a split peak cannot be told from two isomers without them.

## Hex3Cer first use — a recorded negative

Built to test whether material accumulating at Hex2Cer was leaving down the globo branch
(A4GALT -> Gb3), which nothing in the facility libraries could see. Searched against
`20251017_Cer_Inhibitor`, all 42 negative injections, into `Analysis_2026-09-15_Hex23Cer`.

**Nothing.** 0 identifications, 0 rows that matched and then failed a filter, and **0 spectral
matches to any of the 176 entries**, associated or orphan. An independent accurate-mass search over
all 869 negative features found nothing within 8 ppm of any Hex3Cer composition either, in a mass
range the run covered (features span m/z 202-1698) and sampled (619 MS2 scans in the Gb3 window
across 8 injections).

So those cells make no Gb3, and the attenuation from Hex2Cer down to the gangliosides is not
diversion into the globo series. Lc3 and the lacto/neolacto branch are covered by the same masses —
at sum composition Gb3 and Lc3 are indistinguishable — and were equally absent.
