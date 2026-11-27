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
modifiers found are recorded in a table that is not distributed
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

## Excluded from redistribution

No LipiDex 2 material is redistributed. Six libraries derived from it are held locally and are not in
the release archive, because LipiDex 2 grants no licence (no licence file, no terms in the README,
nothing in its distributed archive; see `manuscript/licensing/`): `LipiDex2_Ganglioside.msp`,
`LipiDex2_FAHFA.msp`, `LipiDex2_Ornithine.msp`, `DECOY1_LipiDex2_Ganglioside.msp`,
`DECOY2_HeadGroup_LipiDex2_FAHFA.msp` and `HOLDOUT_LipiDex2_Ganglioside.msp`. They are also out of
the default library sets, so a normal run never loads them.

`Ganglioside_Negative_SumComposition.msp` used to be enumerated from the LipiDex 2 positive library.
Since 2026-09-13 it is built by `scripts/build_ganglioside_library.py` from the Svennerholm glycan
structures with computed ceramide masses and fragments measured in this laboratory, so it carries no
third-party content and ships in the archive.

⚠ `DECOY2_HeadGroup_LipiDex2_FAHFA.msp` also carries a defect: none of its 1,089 entries has the
`DECOY_` name prefix, so a hit on it would never be counted as a decoy. Nothing uses it.

## Two ganglioside libraries, and which is which

`Ganglioside_Negative_SumComposition.msp` (1,232 entries) is the library the **manuscript** used.
Seven class forms: GM3 and GM2 as both NANA and NGNA, and GM1, GD3, GD2 as NGNA. Built by
`scripts/build_ganglioside_library.py` with no flags. **Do not regenerate it with other options** —
the paper's Set 1 and Set 2 configurations name this file, and changing it changes those results.

`Ganglioside_Negative_SumComposition_Extended.msp` (2,112 entries) is the **facility** library and
is what `scripts/run_study.py` loads by default. Built with `--all-classes`. It adds the all-NeuAc
forms of GM1, GD3 and GD2, and both forms of **GM4**, the sialylgalactosylceramide with a single
hexose, which is the simplest sialylated glycosphingolipid and the only one here with one hexose.
GM4 is largely a brain and kidney lipid and is scarce in plasma, which is why it is not in the
manuscript set.

Neither library has a decoy counterpart, so ganglioside rows carry no target-decoy estimate.


## Prenol quinones — a category, not a class

`Ubiquinones_Positive.msp` (4 entries, added 2026-10-02) holds coenzyme Q: `CoQ10` as `[M+NH4]+`
and `[M+H]+`, ubiquinol `CoQ10H2` as `[M+NH4]+`, and `CoQ9` as `[M+NH4]+`. Before it, **no library
here contained a single quinone** — LIPID MAPS category PR was missing in full, which is why
coenzyme Q has never appeared in a delivered table despite being a better peak than the
free-cholesterol positive control on our own data.

It is **not** in the default set: a quinone takes no decoy from the CH2 construction, and its
displacement is unmeasured. Search it with `--extra-library Ubiquinones_Positive`.

⚠⚠ **Ubiquinol is 10 ppm from the ubiquinone M+2 isotope** (882.7334 against 882.7244), inside the
0.01 Da `MS1_TOL`, with an identical MS2. They separate on retention alone. Any run that reports
ubiquinol must be checked with `scripts/check_ubiquinol_rt.py`, which reads
`Ubiquinones_RT_anchors.csv`. See `PROVENANCE_Ubiquinones.md` and the developer notes (PRENOL_QUINONES, not distributed).
