# Chromatography: are the lipids eluting while the gradient is still moving?

`scripts/plot_gradient_vs_elution.py` overlays where each lipid class elutes on the solvent
programme of the method that produced it. A class sitting on a flat part of the programme is
being separated isocratically, and the gradient is contributing nothing to its resolution.

Each dataset carries its own `Method`, because these runs do not share a programme and
plotting one against another's gradient is misleading.

## Methods on record

Read from the `.raw` header the same way `RawMethodExtractor` does: recover UTF-16LE strings
from the first 96 MB, keep lines carrying a time marker or a `LoadingPump.%B/%C` value.

**`M2026` — 25 min.** `Data - 2026\20260602_Dunja_Lipidomics`
B 15→100% over 1.5–5 min, hold to 12, IPA 12→20 min to 95% C, wash, re-equilibrate.

**`M2024` — 30 min.** `Data - 2024\20240212_Martina_BrownLab`
B 15→100% over 0.5–5 min, hold to 10, **IPA 10→22 min** to 95% C, 99% at 26, then wash.

Same B ramp; the difference is the IPA segment, 12 minutes rather than 8, starting 2 minutes
earlier.

## What the runs show

| run | polarity | identified | elution span | % of run |
|---|---|---|---|---|
| CKD rat heart (M2026) | positive | 642 | 6.6–18.7 | 48% |
| CKD rat heart (M2026) | negative | 238 | 6.8–12.7 | 24% |
| Skin organoid (M2024) | positive | 939 | 9.1–25.2 | 54% |
| Skin organoid (M2024) | negative | 409 | 9.2–19.5 | 34% |

Both methods work: polar lipids elute in the B plateau, neutral lipids (TG, CE) on the IPA
ramp. Neither wastes the run the way the older 2025 ATRX data does, where nothing was
identified past 11 min.

## Flow rate is not the same, and it dominates everything

`PumpModule.LoadingPump.Flow.Nominal` in the headers: **150 uL/min in 2024, 300 in 2026** -
exactly 2x. Column oven differs too, 40 vs 50 C.

This matters more than the gradient shape. Retention time scales as 1/flow, so any time-domain
measure of separation - seconds per CH2, peak spacing, co-elution counts in a fixed time window
- is inflated 2x for the 2024 run by flow alone. Peak width scales as 1/flow as well, so
resolution is flow-invariant to first order and the time-domain numbers are simply misleading.

`scripts/compare_methods.py` therefore works in **elution volume** (uL = RT x flow), in which
flow cancels. Homologue selectivity, uL per CH2, median across classes:

| | 2024 / 150 uL/min | 2026 / 300 uL/min |
|---|---|---|
| positive | 61.7 | 65.3 |
| negative | 66.3 | 61.1 |

The two methods are **equivalent**. An earlier version of this file claimed 2024 resolved
homologues about 2x better; that was the flow difference and nothing else.

In volume terms the 2026 IPA ramp is also the shallower of the two - 12-20 min at 300 uL/min is
2400 uL, against 1800 uL for 10-22 min at 150 - despite being shorter in time.

Where they do differ, 2026 is ahead on the classes that elute during the IPA ramp:

| class | 2024 uL/CH2 | 2026 uL/CH2 |
|---|---|---|
| Cer[NDS] (neg) | 83.8 | **122.8** |
| Cer[NS] (pos) | 87.1 | **128.2** |
| TG (pos) | 31.5 | **76.2** |
| SM (neg) | 78.9 | **99.4** |

Elution order is identical between the methods (RT correlation r = 0.987 on 222 shared
positive-mode identifications), so the column chemistry is behaving the same way in both.

## Where this bears on Cer[ADS]

Negative mode is the only polarity that sees Cer[ADS] - LipidBlast carries it as `[M-H]-` and
`[M+FA-H]-` only - so the negative-mode window is the one that matters.

- **Skin, M2024:** 28 Cer[ADS] over 1352 uL of eluent.
- **CKD, M2026:** 4 Cer[ADS] over 268 uL.

That gap is **sample and library, not chromatography**. Skin is genuinely rich in alpha-hydroxy
ceramides; rat heart is not. Per-CH2 selectivity for the dihydroceramides is if anything better
under M2026.

One observation about M2026 survives the flow correction, because it concerns composition
rather than time: negative-mode lipids elute 6.8-12.7 min while the B plateau runs 5-12 min, so
most of that separation happens at constant mobile phase. The predicted consequence did not
appear, though - neg PE selectivity is 70.3 uL/CH2 under M2024 against 72.6 under M2026. High-
organic isocratic elution is still resolving these homologues perfectly well, so the plateau is
not costing anything measurable and is not a reason to change the method.

**Conclusion: keep the 25-min method.** It matches the 30-min method on separation, beats it on
ceramides and TG, and does it in five fewer minutes. The cost is solvent - 7500 uL per run
against 4500.

## Adding a proposed gradient

Add a `Method` in the script and pass it to `figure()`. The classes are drawn from an existing
`Final_Results.csv`, so you can see where they currently land against a programme you are
considering before committing instrument time to it.
