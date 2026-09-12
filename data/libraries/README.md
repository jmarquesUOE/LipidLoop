# Libraries — and why the mobile-phase modifier chooses them

`run_study.py --modifier {formate,acetate}` picks the set. **Default is formate**, which is the
facility method — Edinburgh runs **ammonium formate with formic acid**. Facility results must not
move, and do not.

**Read the modifier from the deposit, do not infer it.** Workbench mwTab carries it verbatim in
`CH:SOLVENT_A`/`CH:SOLVENT_B`, and MetaboLights ISA-Tab in the chromatography protocol. Inferring
it from identification counts fails exactly where it matters: on a positive-mode-only study the
two library sets differed by 2 identifications out of 215, with **zero** negative-mode rows — the
polarity where the modifier changes anything at all.

Of 14 public deposits whose method could be read, **13 are formate** and one is acetate. The
modifiers found are recorded in `paper/validation/metadata/modifiers.csv`
with the quoted evidence for each.

## The modifier is not a preference

Negative-mode adducts follow the mobile-phase modifier. On a study run with ammonium **acetate**,
searching the **formate** libraries named 19 acetate adducts of ordinary dihydroxy ceramides as
`Cer[AP]` **phytoceramides** — a class that was not present at all:

    Cer[AP] t42:0  was really  Cer 18:0;O2/22:0 [M+CH3COO]-   0.3 ppm, RT within 0.01 min
    Cer[AP] t46:1  was really  Cer 20:1;O2/24:0 [M+CH3COO]-   0.1 ppm

An acetate adduct of a dihydroxy ceramide has almost exactly the mass of a deprotonated
phytoceramide two carbons longer with an extra oxygen. Dot products were 982–999.

| library | ceramide `[M+Ac-H]-` entries |
|---|---|
| `LipiDex_HCD_Formic` | **0** |
| `LipiDex_HCD_Acetate` | **2,120** |

Running the same data with `--modifier acetate`: **`Cer[AP]` 19 → 0**, distinct identifications
397 → 430, and recall against the reference study's curated set 87.5% → **88.9%** (unambiguous
67.5% → **72.7%**).

## The set

| | formate | acetate |
|---|---|---|
| in-silico | `LipidBlast_Formic` | `LipidBlast_Acetate` |
| HCD | `LipiDex_HCD_Formic` | `LipiDex_HCD_Acetate` |
| hydroxy | `LipiDex_HCD_Hydroxy` | (shared — carries both) |
| spiked standards | `LipiDex_Splash_ISTD_Formic` | `LipiDex_Splash_ISTD_Acetate` |
| gangliosides, ceramides | shared | shared |

Source: LipiDex 1.1, `src/msp_files/`.

**The SPLASH ISTD libraries are searched by default.** 21 entries, so they cost nothing, and they
give every run a set of masses known *before* it started — the only non-circular input a mass
calibration can have, and a run-quality check that needs no biology.

## Still missing

`LipiDex_HCD_Plants` (29,302) and `LipiDex_HCD_ULCFA` (1,725) are available and not searched.
ULCFA matters for skin and meibum; Plants is a good **negative control** — searching plant lipids
against human plasma measures the false-positive rate, which nobody publishes.
