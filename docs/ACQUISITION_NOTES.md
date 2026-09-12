# Acquisition notes — inclusion list, exclusion list, and a spectral artefact

Written from the CKD rat heart validation set (26 files per polarity, 25 min method, 5 mM
ammonium formate / 0.1% formic acid in all three channels). Three findings for the next method:
one you already planned, two you did not.

---

# 1. Cholesteryl esters go on the inclusion list

**Decided.** `method/CE_inclusion_list.csv` — 40 species as `[M+NH4]+` plus the same 40 as
`[M+Na]+` as a backup adduct, over a 14.5–19.5 min window.

## Why it is needed

One CE is identified across 26 files, from a library of 40. The cause is not the library, not
the mobile phase, and only partly biology — see
[CHOLESTERYL_ESTERS.md](CHOLESTERYL_ESTERS.md) for the full diagnosis. In short:

* **Zero MS2 fall within 0.01 Da of any CE library precursor except `CE 20:4`.** They were never
  fragmented. Nothing downstream can recover that.
* `CE 20:4` `[M+NH4]+` peaks at 2.5e6; `CE 18:1` at 7.4e4, `18:2` 6.9e4, `16:0` 2.8e4, `18:0`
  1.2e4. A 1e4 precursor in a window where TG run at 1e9–1e10 never reaches the top-N list.
* Rats esterify arachidonate preferentially, so `CE 20:4` genuinely dominating is expected and
  not an artefact — but the others are real and measurable, just never selected.

## Two things worth pairing with it

**Lower the source / ion-transfer-tube voltage.** The shared in-source fragment 369.3516 peaks at
7.9e6 — *larger than the intact `CE 20:4` precursor it came from*, and 22× the intact `CE 18:1`
adduct. A large part of the CE population is destroyed before mass selection. A short voltage
ramp on a pooled QC, watching the 369.3516-to-`[M+NH4]+` ratio, would show how much comes back.
With the inclusion list forcing selection, a stronger precursor also means a better spectrum.

**Add an exclusion.** A background ion at **668.571** is selected roughly every 0.11 min through
the whole run — around 180 MS2 per file spent on it. It sits 96 ppm from `CE 18:1` (668.6346),
close enough to look like CE in a coarse search and far enough to never match. It turned out to be
one member of a much larger problem: see §3.

---

# 2. A spectral artefact is costing about a third of all identifications

Not planned, and larger than the CE problem.

**m/z 178.29 is present in 100% of MS2 spectra and is the base peak in ~69% of them — in both
polarities.**

It is not a real ion. A singly-charged organic ion cannot have that mass defect — +0.295 at
m/z 178 is beyond any possible CH formula — and its measured m/z wanders by ±20 mDa, which a real
ion does not.

## Why it matters more than it looks

Every spectrum is normalised to its base peak before scoring, in this pipeline and in LipiDex.
When the artefact *is* the base peak, every genuine fragment is scaled down against it, and the
unmatched artefact intensity also counts against the forward dot product. Both push scores below
threshold.

Measured over four files per polarity, searching the same spectra with and without it:

| | identifications | distinct lipids |
|---|---|---|
| positive, as acquired | 4,128 | 649 |
| **positive, artefact removed** | **5,634 (+36%)** | **744 (+15%)** |
| negative, as acquired | 917 | 239 |
| **negative, artefact removed** | **1,376 (+50%)** | **319 (+33%)** |

The median dot barely moves (922 → 929), so this is not inflating scores across the board — it is
rescuing spectra that were being penalised. The gains land in the abundant real classes:
TG +741, PC +322, PE +80, SM +16.

## The full screen: nothing else comes close

Screened every m/z across four files per polarity, ~65,000 positive and ~50,000 negative spectra:

| polarity | m/z | in % of spectra | base peak in | m/z SD |
|---|---|---|---|---|
| positive | **178.295** | **100%** | **69%** | 10.4 mDa |
| negative | **178.264** | **100%** | **67%** | 10.4 mDa |
| negative | 112.985 | 33% | 6% | 4.9 mDa |

**Only one m/z exceeds 95%**, and it is the same feature in both polarities. Nothing else is
above 33%, and the runner-up is the base peak in only 6% of spectra so it barely touches the
scaling.

One detail that matters for screening: **the artefact's m/z wanders over about 60 mDa**, so a
narrow bin splits it and hides it. Binned at 0.01 it reads as 50% of spectra; binned at 0.10 it
reads as 100%, which is the truth. `find_spectral_artefacts.py` therefore defaults to a 0.10 bin,
and reports the m/z spread alongside the frequency — a stable m/z means a real ion, a wandering
one does not.

One caveat kept in view: `PC[OH]` gains disproportionately, 66 → 249, and that is a class already
flagged as suspect elsewhere (reported at multiple retention times, high blank/sample ratio). Some
of the gain may be in identifications that are not trustworthy for other reasons.

## What to do

**On the instrument**, find it. Present in every spectrum at that intensity, it is most likely
electronic or a persistent background in the collision cell rather than the sample — a clean
solvent injection with the same method would confirm within one run. Worth doing before the next
batch: this is 36% of identifications.

**In software: screened automatically, on every run, per file.** A configured m/z would be right
for one instrument in one period and useless to another lab running the same published method, so
the detector travels with the method rather than a list of masses. `artefacts.py` runs on every
analysis by default, and everything it removes is reported in the log:

    artefact in 20260602_Rat_15_SHAM_Heart_4_wk: m/z 178.295 in 100% of spectra,
    base peak in 77%, m/z spread 10 mDa, no possible formula — 27746 peaks removed

`scripts/find_spectral_artefacts.py` runs the same screen standalone if you want to look before
processing. `exclude_mz` still exists for masses that must be removed by hand.

### How it avoids deleting real fragments

This is the risk that matters. Lipidomics has genuinely common fragments, and removing one would
be far worse than the artefact. Three conditions must hold together:

1. **Present in ≥95% of spectra**, whatever the precursor. Real class fragments are not:
   phosphocholine 184.0733 is in 9.9% of positive spectra here, because only some precursors are
   phosphocholines.
2. **Base peak in ≥20%.** This is what makes it damaging rather than merely present.
3. **No possible formula** — the decisive one. A singly-charged CHNOPS ion cannot have an
   arbitrary mass defect: the most hydrogen-rich possibility, a saturated alkyl cation, runs at
   about 0.00115 per unit m/z, so nothing real at m/z 178 exceeds ~0.23. The artefact sits at
   0.290.

| fragment | mass defect | bound | verdict |
|---|---|---|---|
| phosphocholine 184.0733 | 0.074 | 0.232 | possible — kept |
| choline 104.1070 | 0.108 | 0.145 | possible — kept |
| sphingosine 264.2686 | 0.261 | 0.334 | possible — kept |
| cholestadiene 369.3516 | 0.344 | 0.454 | possible — kept |
| **178.295** | **0.290** | **0.226** | **impossible — removed** |

Condition 3 is what makes it safe to run unattended. Frequency alone would eventually delete a
real fragment on a sample type where one class dominates. An impossible mass cannot be real
whatever the sample.

### Effect with screening on, both polarities, 26 files each

| | before | after |
|---|---|---|
| positive, reference molecules recovered | 97.7% | **98.6%** |
| positive, genuinely absent | 10 | **6** |
| negative, reference molecules recovered | 95.3% | **98.4%** |
| negative, genuinely absent | 9 | **3** |

The reference itself was produced from these same artefact-laden spectra, so it suffered the same
suppression — which is why we now recover nearly all of it. **The corollary is that molecules
only we report rose sharply, 87 → 137 in positive and 47 → 101 in negative, and those are not
validated.** Some will be real lipids the artefact was hiding from both pipelines; some will not.

---

---

# 3. Half the MS2 duty cycle is being spent on contamination

The 668.571 ion in §1 was not a one-off — it is one PEG oligomer out of hundreds of contaminant
ions. Screening every precursor across 24 files per polarity (20 samples + 4 QC; the blanks are
used as the blank test, not counted as samples) and keeping only those **never identified in any
file**, **eluting across more than half the gradient**, and **equally present in the blanks**:

| | precursors | MS2 spent on them | share of duty cycle |
|---|---|---|---|
| positive | 285 | 166,497 / 386,881 | **43.0%** |
| negative | 343 | 200,661 / 299,574 | **67.0%** |

and the sources are identifiable from the repeat spacing between members:

| polarity | dominant series | repeat | share | where it comes from |
|---|---|---|---|---|
| + | polysiloxane | 74.019 (C₂H₆OSi) | 24.5% | column bleed, PDMS from tubing and septa |
| − | **sodium formate** | 67.987 (CHO₂Na) | **44.6%** | formic acid + sodium leached from glassware |
| − | PTFE/PFPE | 99.994 (C₂F₄) | 7.4% | fluoropolymer tubing, pump oil |

The two polarities have almost disjoint problems: siloxane dominates positive and is invisible in
negative, salt clusters dominate negative and are a rounding error in positive.

**The sodium formate finding is a bench problem before it is a software one.** Nearly half the
negative-mode duty cycle is going to [(HCOONa)ₙ + HCOO]⁻ clusters running from m/z 316.9451 up past
1500 — the mobile phase reacting with sodium out of glassware. Excluding them recovers the budget;
using plastic and low-sodium water stops them being made.

**Is the freed budget worth anything?** Measured rather than assumed. In positive QC_01, 2,916 of
5,200 MS1 features never receive an MS2 at all, and simulating the reallocation cycle by cycle the
freed slots could reach 2,437 of them (84%); negative, 335 of 376 (89%). Even among the 500 most
abundant positive features, 52 never got fragmented. That is an upper bound — the instrument
re-ranks live and its own dynamic exclusion will not follow the simulation — but the budget exists
and there is plenty waiting for it.

The current lists are in [`method/`](../method/README.md) and the full analysis behind them in
[`results/wasted_ms2_Pos.md`](../results/wasted_ms2_Pos.md) and
[`results/wasted_ms2_Neg.md`](../results/wasted_ms2_Neg.md). They are generated per batch and never
reused — column bleed depends on the column, PEG on the labware, sodium formate on the glassware:

    python scripts/find_wasted_ms2.py --mzml <run>/*.mzML --results <results-dir> \
        --blank <run>/Blank_*.mzML --polarity + \
        --out method/Exclusion_List_Pos.csv --thermo method/Exclusion_List_Pos_Thermo.csv \
        --report results/wasted_ms2_Pos.md --simulate

**Read [EXCLUSION_LIST.md](EXCLUSION_LIST.md) before putting any of it in a method** — excluding a
precursor is the only decision here that reprocessing cannot undo. The collision check has been
run (`scripts/check_exclusion_safety.py`): at a ±10 ppm exclusion tolerance **no entry on either
list collides with a lipid identified in this experiment**, and the `_Safe` variants, which also
drop anything within 10 ppm of any library lipid, still recover 98.4% / 99.7% of the duty cycle.
Use those.

---

# 4. Oxylipins need an inclusion list before they need a library

Asked for as a library; the data says the library is the second thing to build.

| | MS1 detection, skin negative | MS2 spectra across 10 files |
|---|---|---|
| HODE (m/z 295.2279) | 40 of 63 files, area 1.4e7 | **0** |
| HETE/EET (319.2279) | 27 files, area 5.7e6 | **1** |
| 9-HpODE, DiHOME | 2 files each | 0, and 11 of background |

HEPE, DiHETrE, HDHA, PGF2α, TXB2, RvD1, LXA4, PGA2 and 12-HHT are **not detected at all**, in
either dataset. **The oxylipins that are there are essentially never fragmented** — this is the
cholesteryl-ester problem in a more severe form. They sit two to three orders of magnitude below
the structural lipids that win every top-N decision, on a run where 46–67% of the duty cycle
already goes to contamination.

`method/Oxylipin_inclusion_list.csv` — 18 distinct precursors covering 54 species. The retention
window is a placeholder from where the C18/C20 acids elute here and **should be narrowed with
standards**. `data/libraries/Oxylipins_Negative.msp` exists alongside it and is **deliberately not
in the default set**: a library cannot identify a spectrum that was never acquired, and adding it
now would generate matches against two-file features with no fragmentation behind them.

Realistically, oxylipins at these levels want their own extraction and a targeted method. The
inclusion list is what makes them possible on this one.

---

## Summary for the next method

| change | expected gain |
|---|---|
| CE inclusion list (`method/CE_inclusion_list.csv`) | CE identified beyond `CE 20:4` |
| find and remove the 178.295 artefact | ~36% more identifications, everything |
| **exclusion list, negative mode** (`method/Exclusion_List_Neg_Thermo.csv`) | **67% of duty cycle back** — mostly sodium formate clusters |
| **exclusion list, positive mode** (`method/Exclusion_List_Pos_Thermo.csv`) | **43% of duty cycle back** — mostly siloxane column bleed |
| less sodium at the bench (plastic, low-Na water) | stops the negative-mode clusters being made |
| **oxylipin inclusion list** (`method/Oxylipin_inclusion_list.csv`) | 18 precursors currently detected but never fragmented |
| lower source voltage | stronger intact CE precursors, better spectra |

The inclusion and exclusion lists belong in the same method and work together: the exclusion list
frees the budget, the inclusion list decides where some of it goes.
