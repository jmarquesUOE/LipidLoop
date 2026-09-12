# Why only one cholesteryl ester is identified

One CE — `CE 20:4` — is identified across 26 files, from a library of 40. `CE 16:0`, `18:0`,
`18:1` and `18:2` are all visible in MS1 and none is identified. Four candidate explanations,
three of them wrong.

## Not the library

All four are present at the correct precursor masses: `CE 16:0` 642.6189, `CE 18:0` 670.6502,
`CE 18:1` 668.6346, `CE 18:2` 666.6189, each with 19 peaks and 369.3516 as base.

## Not the mobile phase

Checked in the acquisition method rather than assumed — all three channels are
`5 mM Ammonium Formate / 0.1% Formic acid`. Ammonium is present, and `CE 20:4` proves the
`[M+NH4]+` chemistry works.

## Partly real biology

| | [M+NH₄]⁺ apex |
|---|---|
| CE 20:4 | **2.5e6** |
| CE 18:1 | 7.4e4 |
| CE 18:2 | 6.9e4 |
| CE 16:0 | 2.8e4 |
| CE 18:0 | 1.2e4 |

`CE 20:4` is 30–200× the others. Rats esterify arachidonate preferentially — unlike humans,
where CE 18:2 dominates — so an arachidonate-dominated CE profile is expected here and is not an
artefact. The others are genuinely minor. That alone does not explain a *total* absence of
identification, though.

## The actual reason: they were never fragmented

**Zero MS2 spectra fall within 0.01 Da of any of those library precursors.** Only `CE 20:4` has
any — 9 across ten files, all at 16.07–16.12 min, all containing 369.3516.

At first sight DDA looks like it *is* selecting `CE 18:1`: 180 MS2 have a precursor within 0.5 Da
of 668.6346. They are all at **668.571**, which is 64 mDa — **96 ppm** — away, and they fragment
to 178.29 and 89.06 with no 369.35 at all. That is a persistent background ion, selected roughly
every 0.11 min across the entire run, and it has nothing to do with cholesteryl esters.

So a precursor at 1–7e4, in a window where the triacylglycerols run at 1e9–1e10, never reaches the
top-N list.

**In-source fragmentation makes it worse.** The shared fragment 369.3516 peaks at 7.9e6 — larger
than the intact `CE 20:4` precursor it comes from, and 22× the intact `CE 18:1` adduct. A
substantial part of the CE population is being converted to a fragment before mass selection, so
the intact adduct that DDA has to see is depleted by the source.

## How to fix it, in order of value

1. **Targeted inclusion list.** The CE `[M+NH4]+` masses with a 15–19 min retention window. Direct
   and certain: it removes the intensity competition entirely, and the list is nine masses.
2. **Reduce in-source fragmentation** — lower the source or ion-transfer-tube voltage. Some of the
   7.9e6 sitting in the fragment channel comes back to the precursor, where it can be selected.
   Worth a short voltage ramp on a pooled QC, watching the 369.3516-to-`[M+NH4]+` ratio.
3. **Exclusion list for the background.** The 668.571 ion alone takes ~180 selections per file.
   The 178.29 artefact appears in most spectra across this dataset — including the negative-mode
   ganglioside spectra — and is worth chasing as contamination in its own right.
4. **Do not expect a software fix.** With one identified CE there is no retention model to fit
   (`rtls.py` needs six), so identifications cannot be extended from retention. The MS1 features
   are present and could be quantified at MS1 level, but that is an assertion from mass and
   retention alone and should be labelled as such — the `Identification Source` column exists for
   exactly that.

## Do not quantify CE on 369.3516

Every CE gives it, so it is shared across the whole class and cannot be attributed to any one
species. It is also the reason free cholesterol cannot be identified in this data at all — see
[STEROLS.md](STEROLS.md).
