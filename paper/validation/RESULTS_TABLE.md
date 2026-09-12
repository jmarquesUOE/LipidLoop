# Public corpus: what was downloaded, how it was screened, what was kept

Public lipidomics deposits matching untargeted LC-MS/MS were downloaded from MetaboLights,
Metabolomics Workbench and MassIVE without pre-filtering on platform, tissue or reported outcome.
Every file, not every deposit label, was then screened before searching:

- **Acquisition type**, from precursor-isolation-window reuse measured directly on the mzML: DDA
  fragments each precursor once per cycle (narrow window, low reuse); DIA/SWATH revisits a wide
  window repeatedly (high reuse). Files scored DIA, all-ion or MS1-only are excluded per file, not
  per study -- several deposits mix acquisition types within one submission (e.g. `ST000991`: 153
  SWATH + 5 DDA + 1 MS1-only; only the 5 DDA files are searched).
- **Chromatography**, from the deposit's own method record: reversed-phase kept, HILIC excluded. A
  lipid library searched against HILIC separation is a coverage mismatch, not a lipidomics test --
  this is why four further public deposits, titled as untargeted metabolomics benchmarks and split
  evenly RP/HILIC, are excluded from this corpus entirely (`ST004503`/`ST004626` are RP but are
  metabolomics benchmarks, not lipidomics studies; `ST004650`/`ST004651` are the same benchmark on
  HILIC).
- **Genuine independence**, checked rather than assumed from staging metadata: one deposit's
  combined "Sciex" arm turned out to be its Normal and HighMass arms staged again under a third
  name (below) and was excluded as a duplicate, not counted as an extra study.
- **In-house primary validation sets** (rat heart, skin organoid, PDAC KPC) are held out of this
  corpus entirely, whatever their own identification counts -- they are scored against this lab's
  own Compound Discoverer + LipiDex reference elsewhere (Table 2), not against a public deposit's
  published list, and counting them here would double them into both roles.

Fifteen studies passed all three criteria and are searched below. Injections staged but not
searched, and deposits excluded outright, are listed with reasons at the end of this table.

| dataset | platform | inj. (searched) | compound groups | ID rows | species delivered | species searched | decoy FDR | deposit species | shared | agreement |
|---|---|---|---|---|---|---|---|---|---|---|
| MTBKS222_Waters | — | 17 (17) | 1,517 | 260 | 191 | 202 | 0.00% | 145 | 114 | **78.6%** |
| MTBKS222_Thermo | — | 8 (8) | 1,440 | 467 | 315 | 337 | 0.00% | 309 | 201 | **65.0%** |
| ST003077 | — | 14 (14) | 4,482 | 1,708 | 446 | 832 | 0.12% | 415 | 246 | **59.3%** |
| MTBKS222_IMS_Normal | — | 12 (12) | 3,199 | 711 | 400 | 489 | 0.00% | 503 | 276 | **54.9%** |
| MTBKS222_IMS_HighMass | — | 12 (12) | 3,995 | 675 | 408 | 488 | 0.00% | 497 | 266 | **53.5%** |
| MTBKS222_Agilent | — | 12 (12) | 6,719 | 749 | 423 | 476 | 0.00% | 550 | 281 | **51.1%** |
| MTBLS5163 | Bruker maXis II | 50 (50) | 990 | 359 | 218 | 309 | 0.00% | 59 | 29 | **49.2%** |
| ST003514 | Agilent 6545 QTOF | 30 (20) | 1,547 | 562 | 254 | 346 | 0.00% | 370 | 181 | **48.9%** |
| MTBKS222_BrukerBAF | — | 13 (13) | 3,061 | 206 | 146 | 163 | 0.00% | 261 | 115 | **44.1%** |
| ST004797 | Sciex ZenoTOF 7600 | 224 (224) | 26,353 | 2,401 | 423 | 672 | 0.29% | 1066 | 271 | **25.4%** |
| MSV000095868 | Agilent QTOF | 184 (184) | 8,101 | 113 | 71 | 98 | 0.00% | — | — | — *none deposited* |
| ST000991_DDA | — | 5 (5) | 1,027 | 325 | 228 | 244 | 0.61% | — | — | — *none deposited* |
| ST002705 | Bruker micrOTOF-Q II | 28 (28) | 17,654 | 479 | 254 | 361 | 0.21% | — | — | — *none deposited* |
| ST003052 | Thermo Q Exactive HF | 120 (120) | 7,802 | 275 | 121 | 175 | 0.36% | — | — | — *none deposited* |
| ST004797_oxtg | — | 19 (19) | 2,913 | 617 | 236 | 318 | 0.00% | — | — | — *none deposited* |

**Totals** — 15 studies, 748 injections (738 searchable), 90,800 compound groups, 9,907 identification rows, 13 decoys (0.13% FDR).
10 deposits publish identifications: 4,175 species, agreement 25.4–78.6%.

⚠ **Regenerated 2026-09-08** against `Analysis_recal_2026-09-04` (the `--auto-calibrate` re-run
that fixed the blank-classification bug, `2e64de1`) — supersedes the numbers above from the prior
`Analysis_calibrated`/`Analysis_v7` runs. Two of the three highest-stakes datasets moved more than
1pp: **MTBLS5163 66.1%→49.2%** (-16.9pp) and **MTBKS222_Waters 85.5%→78.6%** (-6.9pp); the third,
ST004797, barely moved (26.5%→25.4%). Everything else shifted ≤1.1pp. Not yet reviewed for
whether these two drops are the bug fix correctly removing inflated matches, or a new regression —
check before quoting either number in prose. `PDAC_KPC_Set3` (Table 2's in-house "Set 3") was
briefly leaking into this table as a 16th study before this regeneration — its exclusion key was
stale (`PDAC_KPC`, no such directory), fixed in `results_table.py`.

**Excluded**, with reason:

- `MSV000094718` — methods deposit (bovine liver, vendor contaminants), not a study
- `MTBKS222_BrukerTDF` — timsTOF conversion unresolved — forward dot median 8 vs Thermo 926
- `MTBKS222_Sciex` — the same injections as IMS_Normal + IMS_HighMass, staged again under a third
  name. Caught by reading instrument configuration and file identity out of the mzML, then scoring
  retention time two ways: against the deposit's own Sciex MAF the 24 files disagree by +4.22 min
  at 0% agreement inside a 0.2 min window; against `IMS_Normal`'s own MAF instead, the same files
  land at -0.01 min and 85%. Not a coincidental overlap -- the deposit's own reference for this arm
  does not describe these files, `IMS_Normal`'s does.
- `MTBLS2016` — no MS2 in any of 147 files, despite being described as DDA
- `PDAC_KPC_Set3` — in-house, held separately (primary validation set, not a public deposit)
- `Rat_Heart` — in-house, held separately (anonymised primary validation set, not a public deposit)
- `Rat_Heart_Lumos` — in-house, held separately
- `ST000987` — all-ion (in-source CID), no precursor selection
- `ST004503` — metabolomics benchmark
- `ST004626` — metabolomics benchmark
- `ST004650` — metabolomics benchmark, HILIC
- `ST004651` — metabolomics benchmark, HILIC
- `Skin_QEplus` — in-house, held separately
