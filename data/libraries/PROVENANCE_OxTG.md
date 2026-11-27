# OxTG_Positive.msp — generated, not imported, and not yet enabled

**752 entries, 376 parent sum compositions × 2 oxygen counts, `[M+NH4]+`.**
Built by `scripts/build_oxtg_positive_library.py`, which is the reproducible record.

## Why generated

No vendor ships one. Zero `OxTG` COMPOUNDCLASS entries and zero TG names carrying `;O` across
MS-DIAL's Tandem Mass Spectral Atlas VS69 — 1.06 million spectra, both polarities. MS-DIAL produces
those annotations from rules at search time, so there is nothing to import.

## What constrains it

**Autoxidation requires a bis-allylic CH₂** — the carbon between two double bonds. A chain with
DB ≥ 2 oxidises readily; 18:1 is far slower and 16:0/18:0 do not go at all. Only parents carrying
an oxidisable chain are enumerated. That constraint is the difference between 752 entries and
~60,000, and it is the reason this library **cannot** reproduce the two fully saturated OxTG
annotations found in a public deposit (`TG 16:0_16:0_21:0;O2`, `TG 16:0_18:0_21:0;O2`).

⚠ It *can* generate `TG[OH] OH-54:2` and `OH-53:2`, which that deposit also annotates as
impossible molecular species (`18:1_18:1_18:0;O3`). That is not a failure: at SUM level those
compositions have oxidisable molecular forms (`18:2_18:0_18:0`). The deposit's error is at
molecular resolution, which this library deliberately does not claim.

## Naming — and why not the LSI `;O` form

⚠ `TG 52:4;O` collapses onto ordinary TG. Verified against the live code: `sum_composition` returns
`"TG 52:4"` and `parse_sum_name` judges it on the ordinary TG retention surface, so the entry would
be merged by duplicate-name detection, split-peak summing and adduct-pair removal, and inherit a
surface fitted to molecules 48 Da lighter. The `[OH]`/`OH-` idiom — LipiDex's own, for hydroxylated
phospholipids — keeps its class and its surface: `TG[OH] OH-52:4` → class `TG[OH]`, prefix `OH-`.

## Sum composition, not molecular species

The oxidised chain is well determined: its neutral loss differs from an ordinary fatty acid by
**36 mDa** (`FA 18:2;O` 294.2195 vs `FA 19:2` 294.2559), resolvable at any sensible tolerance. The
two unoxidised chains are determined only as a PAIR, because the surviving diacyl ion gives their
sum. A molecular-resolution library would give `TG 18:1_18:2_18:2;O` and `TG 18:0_18:3_18:2;O`
near-identical spectra and the dot product could not separate them.

## ⚠⚠ Fragment INTENSITIES are provisional — this is the main limitation

Masses are exact. Intensities are **not measured**, and the reason matters: the only public deposit
annotating OxTGs is profile-mode data whose MS2, even after peak-picking, splits one fragment into
several centroids —

    precursor 910.7470:  599.5012 (100%)  599.5082 (77%)  599.5151 (20%)  599.4752 (15%)

all one peak. Any intensity statistic from that is meaningless. The pattern used asserts only what
the fragmentation demands: the acyl-loss ions exist, and the oxidised-chain loss is prominent.
`--intensities` re-parameterises it from clean data without regenerating anything else.

**Treat a match as evidence of composition, not of confidence.**

## ⚠ NOT wired into `run_study.py`

Deliberately. It is unvalidated on real data and its intensities are provisional, so enabling it
would put an in-silico library with an invented intensity model into production runs. The
bile-acid import is the precedent for what "we have the library now" is worth: 36 correct entries,
zero identifications, because the chromatography did not retain the analytes.

Enable it after:

1. **Acquisition.** The oxidised-chain discrimination needs ≥15K MS2 resolution — an ion trap at
   unit resolution cannot separate 36 mDa and destroys the only diagnostic. MS1 wants ≥120K; note
   Orbitrap resolution is quoted at m/z 200 and falls as 1/√(m/z), so 120K @200 is ~56K at m/z 920,
   just under the ~60K needed to resolve the `;O3` isobar from `TG (n+4):(d+4)`.
2. **Decoys.** `DECOY1_OxTG_Positive.msp` is built — 752 entries, none skipped. Run it.
3. **A positive control.** Spike an oxidised standard into an extracted blank. Nothing here shows
   the pipeline can detect an OxTG that is definitely present.
4. **The artefact question.** Oxidised TGs form in the vial. See
   `Projects/Oxidised_TGs/REPORT_ST004797_oxidised_TG.md` — three explanations tested there, all
   failed, and the planned 72 h stability run is what settles it for our own method.
