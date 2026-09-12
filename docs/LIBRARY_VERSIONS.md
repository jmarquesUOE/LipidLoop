# Does the LipiDex 2 library increase identifications?

Short answer: **no, not from the updated versions of the libraries we already use, and barely
from the new ones.** Measured rather than assumed.

## The three shared libraries are content-identical

`LipiDex2_HCD_Formic.msp` is 56.8 MB against v1's 45.3, which looks like growth and is not.
Comparing on lipid **and** adduct, across all three libraries we search:

| library | v1 entries | v2 entries | only in v2 |
|---|---|---|---|
| LipiDex HCD Formic | 97,082 | 97,082 | 0 |
| LipidBlast Formic (v2: `InSilico_LipidBlast`) | 51,772 | 51,772 | 0 |
| HCD Hydroxy | 15,904 | 15,904 | 0 |
| **total** | **134,529** | **134,529** | **0** |

Every entry is present in both. The 25% file growth is per-entry metadata that v1 kept in a
free-text comment or not at all — `PRECURSORTYPE`, `FORMULA`, `IONMODE`, `COMPOUNDCLASS`,
`OPTIMALPOLARITY`, `ISLIPIDEX`.

A first comparison suggested 9,316 cardiolipins differed. They do not: v1 writes
`CL 14:0_14:0_14:0_14:0 [M-2H]2-` and v2 writes the same lipid with the adduct in its own field
and spelled `[M-2H]-2`. A notation change, not a coverage change.

## The genuinely new libraries add 0.17%

LipiDex 2 ships classes v1 has no library for. Searched over all 26 positive-mode files
alongside the usual three:

| | library spectra | identifications |
|---|---|---|
| v1 three libraries | 152,280 | 24,060 |
| plus FAHFA, Ganglioside, Ornithine | 159,114 | **24,102** |

| new library | entries | identifications |
|---|---|---|
| LipiDex2_FAHFA | 1,089 | **0** |
| LipiDex2_Ornithine | 1,440 | **0** |
| LipiDex2_Ganglioside | 4,305 | **42** |

Nothing was displaced — the three original libraries return exactly 17,286 / 6,358 / 416 in both
runs — so the 42 are purely additive. That is **0.17%** of identifications.

All 42 are GM3-NANA, four species (`d18:1_16:0`, `_18:0`, `_23:0`, `_24:0`), eluting tightly at
8.3, 8.6, 9.7 and 9.9 min. GM3 in heart tissue is plausible and the retention order by chain
length is right, so these are not obviously spurious. They are also not strong: dot products
mostly 510–800 against a 500 floor, on library entries carrying eight peaks. **Gangliosides are
normally analysed in negative mode**, and a positive-mode assignment on a middling score is worth
confirming before it is used.

Ornithine lipids returning nothing is expected — they are bacterial. FAHFA returning nothing in
heart is unsurprising at these concentrations.

## Worth having anyway: the parser now reads both dialects

The v2 format was not readable here, and the failure mode was quiet rather than loud. v1 puts the
adduct inside the name; v2 leaves the name bare and puts the adduct in `PRECURSORTYPE`. Reading a
v2 file with a v1 parser yields names with no adduct, and since polarity is read off the name,
**every entry silently looks negative-mode** — the same class of bug that made our negative mode
return zero identifications. `msp.py` now reads both and normalises v2 to the v1 shape, so
nothing downstream knows the difference.

## The ganglioside library is now in the default set

Added, and measured in place. It is the last entry in `DEFAULT_LIBRARIES`, appended rather than
inserted so the existing load order — which is what reproduces the reference run — is untouched.

| | without | with |
|---|---|---|
| library spectra | 152,280 | 156,585 |
| output rows | 5,084 | 5,081 |
| distinct molecules | 529 | **534** |
| reference molecules recovered | 416 (97.7%) | **416 (97.7%)** |

Purely additive: nothing displaced, five molecules gained, agreement with the reference
unchanged. The five are

| RT | m/z | identification | detected in |
|---|---|---|---|
| 8.25 | 1153.7263 | GM3-NANA d34:1 | 16 files |
| 8.99 | 1209.7888 | GM3-NANA d18:1_20:0 | 7 |
| 9.43 | 1237.8195 | GM3-NANA d18:1_22:0 | 15 |
| 9.68 | 1251.8339 | GM3-NANA d41:1 | 5 |
| 9.94 | 1265.8512 | GM3-NANA d18:1_24:0 | 23 |

Retention rises monotonically with chain length across all five, which is what a real homologous
series does and a coincidence would not. Note the retention model cannot check that here — five
members is below its six-point minimum — so the ordering is a manual observation, not a filter
result.

**The caveat stands.** These sit on dot products of roughly 510–800 against a 500 floor, on
library entries carrying eight peaks, and gangliosides are normally analysed in negative mode.
Treat them as worth confirming rather than as settled.

## Conclusion

Switching to the v2 versions of the libraries we already use gains nothing measurable. The
ganglioside library gains five species and is now on by default. The FAHFA and ornithine
libraries are not relevant to this tissue and are not included.

The real improvements in LipiDex 2 are in the software, not the spectra — see
[LIPIDEX2.md](LIPIDEX2.md).

## Confirming the gangliosides in negative mode

The GM3 calls came from a positive-mode library on middling scores, so they were worth checking
independently. **A library search cannot do it**: `LipiDex2_Ganglioside.msp` is 4,305 entries and
every one is `[M+H]+`. Neither LipidBlast nor the LipiDex HCD libraries carry a single ganglioside
entry in any polarity. There is nothing to search negative-mode data against.

A structural fragment can. GM3 is mono-sialylated, and in negative mode sialic acid gives a
diagnostic dehydrated NeuAc at **m/z 290.088** — no library required.

| species | [M+H]+ | [M−H]− | MS1 scans at [M−H]− | MS2 taken | contain 290.088 |
|---|---|---|---|---|---|
| GM3-NANA d34:1 | 1153.7263 | 1151.7117 | 38 | 6 | **6** |
| GM3-NANA d18:1_20:0 | 1209.7888 | 1207.7742 | 55 | 12 | **6** |
| GM3-NANA d18:1_22:0 | 1237.8195 | 1235.8049 | 72 | 7 | **7** |
| GM3-NANA d41:1 | 1251.8339 | 1249.8193 | 48 | 6 | **6** |
| GM3-NANA d18:1_24:0 | 1265.8512 | 1263.8366 | 118 | 6 | **6** |

All five have real MS1 signal at the deprotonated mass at the right retention time, all five were
fragmented, and **every one of those spectra carries the sialic acid fragment**.

The control is what makes that mean something. Across negative QC_01, a peak within 0.01 of
290.088 appears in **16 of 13,109 MS2 spectra (0.1%)**, and in **14 of 2,719 spectra with
precursor above 1100 m/z (0.5%)**. A fragment present in 0.5% of comparable spectra and in 100%
of these is not a coincidence. **The gangliosides are real.**

One honest caveat: the observed fragment sits at 290.0858–290.0864 against a theoretical 290.0881,
about 7 ppm low. That is worse than this instrument's usual accuracy and is more than one would
like on a diagnostic ion, though low-mass fragments are the least well calibrated part of the
range. It does not change the conclusion, given the specificity.

### There are more gangliosides here than we identify

Sweeping negative mode for that fragment finds **28 distinct precursors carrying it, of which we
identify 5**. The unidentified ones are not scattered — they form the rest of the homologous
series:

| m/z | RT | identified |
|---|---|---|
| 1151.7 | 8.26 | **yes** |
| 1177.7 | 8.31 | no |
| 1179.7 | 8.61 | no |
| 1205.7 | 8.67 | no |
| 1207.8 | 8.99 | **yes** |
| 1233.8 | 9.07 | no |
| 1235.8 | 9.44 | **yes** |
| 1249.8 | 9.68 | **yes** |
| 1261.8 | 9.45 | no |
| 1263.8 | 9.95 | **yes** |
| 1277.8 | 10.19 | no |

The identified members step by 28 Da — two carbons — with retention rising each time. Each
unidentified member sits 2 Da below an identified one and elutes slightly earlier, which is one
additional double bond behaving exactly as it should. This is the retention-model pattern from
[RTLS.md](RTLS.md), visible by eye.

**The limit on ganglioside coverage here is the library, not the data.** The positive-mode
library finds the saturated members that happen to ionise well as `[M+H]+`; the negative-mode
data plainly contains the unsaturated ones too, and there is no library in either LipiDex version
to identify them with.

## `LipiDex_HCD_ULCFA` and `LipiDex2_FAHFA` in negative mode — measured 2026-08-28

The note above ("FAHFA returning nothing in heart is unsurprising") is about a **positive-mode**
survey. Negative is the polarity these species actually ionise in, and it tells a different story
— on Skin_QEplus (Neg, decoy-matched, one library added at a time against the `Analysis_v6`
baseline of 970 identifications / 319 delivered species / 0.51% FDR):

| library | real rank-1 hits | decoy hits | delivered species lost | gained |
|---|---|---|---|---|
| `LipiDex_HCD_ULCFA` (full, 1,725 entries) | 310 | 0 | 0 | 0 |
| `LipiDex_HCD_ULCFA_Deprotonated` (1,426, acetate dropped) | 175 | 0 | 0 | 0 |
| `LipiDex2_FAHFA` (1,089 entries) | 113 | 0 | **5** | 1 |

⚠ **First measurement of this was worthless — `--extra-library` was a silent no-op.** `run_study.py`
accepted the flag, parsed it into `args.extra_libraries`, and never copied it into the config dict.
Three runs reported "identical to baseline in every number," which read as a clean null result and
was actually three runs that never searched the library under test — only its decoy loaded. Fixed
in `0e2ae63`, with a regression test that drives `main()` through `sys.argv` (no earlier test
caught it because every one calls `build_config` directly with a hand-built `extra` dict — exactly
the dict `main()` was failing to assemble) and a shell-level check that greps the written
`run_config.json` for the library stem before trusting any run's numbers.

**ULCFA is clean on this dataset — and inert.** Both versions win real spectra at the search stage
(310 / 175, zero decoys) but none clear `MIN_DOT_PRODUCT`/`MIN_REV_DOT_PRODUCT` to reach the
delivered table: 0 lost, 0 gained, either way. The deprotonated subset's surviving hits carry a
substantially higher median dot (671 against 432) and zero acetate contamination, confirming the
filter works — acetate entries (`[M+FA-H]-` differs from `[M+HCOO]-` by exactly one CH2,
`modifier.py`) accounted for 160 of the full library's 310 hits and are a formate-method false-match
risk by construction. Skin has no genuinely ultra-long negative-mode phospholipid signal at the
confidence this method reaches, so this measurement bounds the library's cost at zero without yet
proving its benefit — the real test is a matrix that needs 30+ C chains.

**FAHFA displaces five real identifications, and does not earn the displacement.** All 113 winning
FAHFA matches carry **`Purity 0`** — the identical no-chain-evidence signature that failed
`Cer[EOS]`'s decoy test (`8c4b34b`). The five lost species — `FA 16:1`, `FA 18:1`, `FA 18:2`,
`FA 20:1`, `FA 22:1` — were all `Identification Source = RT model` in baseline, and the replacement
rows were not what they first appeared to be:

⚠⚠ **The displacement mechanism was a naming bug, not a genuine identification.** `FAHFA
18:1-(O-18:0)` won those spectra, but `sum_composition()` split chains on `_` — FAHFA uses
`-(O-...)` — so the second chain was silently dropped and the delivered row read `FAHFA 18:1`: a
name byte-identical to a real single free fatty acid's sum composition, for a 563 Da,
chain-unresolved molecule. Fixed in `6bfe6d5` (`re.match` → `re.fullmatch`, so an unrecognised
delimiter now falls through to the function's own documented contract — leave the name alone —
instead of truncating to whatever prefix happened to match). General fix: any future library using
a non-`_` chain delimiter would have hit the same silent truncation.

**FAHFA is withdrawn on the same grounds as `Cer[EOS]`: Purity 0 on every winning match, real
displacement of delivered identifications, and no chain annotations to build an informative decoy
from** (confirmed at the raw-file level earlier — 1,089 entries, 3 peaks each, zero annotated).
Left on disk, matched decoy at `DECOY2_HeadGroup_LipiDex2_FAHFA.msp`, not in `DEFAULT_LIBRARIES`.

## `FAMLS_Measured_Negative.msp` — a real, measured FA library (facility use, not the default set)

Added 2026-08-30. 26 entries, real averaged apex-windowed MS2 spectra from authentic IROA-FAMLS
free-fatty-acid standards run on the Lumos — not the templated single-CO2-loss-peak entries
`FreeFattyAcids_Negative.msp` ships (every entry there is one combinatorially-generated fragment,
never measured). Full provenance per entry in its `Comment` field: original standard name,
measured RT, how many scans were averaged.

⚠ **Not in `DEFAULT_LIBRARIES` or `FATTY_ACIDS` (`scripts/run_study.py`) — deliberately.** Adding
it there would pull it into every future run through this script, manuscript validation reruns
included, with no displacement/collision testing done. This is facility material: opt in per run
with `--set fatty_acid_libraries='["data/libraries/FreeFattyAcids_Negative.msp",
"data/libraries/FAMLS_Measured_Negative.msp"]'`, which needs zero changes to any default and
touches nothing for a run that doesn't ask for it.

Before this goes into any default set: the same displacement/decoy-safety work already done for
every other library addition in this project (§4, and the `Cer[EOS]`/ULCFA/FAHFA precedent) --
matched decoys, a collision check against every real target, measured effect on delivered
identifications. None of that has been done for this library yet.

Full facility write-up, the RT model built alongside it, and the RT-only predictions for
compounds without usable MS2 (the IROA-FAMLS standards run, authors' data)
(`FAMLS_Measured_Negative.msp`'s own source, `RT_MODEL.md`, `FAMLS_RT_predicted.md`,
`RT_MS2_summary.md`).
