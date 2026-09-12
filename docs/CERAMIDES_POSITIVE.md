# Positive-mode ceramides and glucosylceramides, and the 85% of a library that could only do harm

Of the eight free-ceramide subclasses LipiDex enumerates, **only `Cer[NS]` has positive-mode
entries**. `Cer[AS]`, `Cer[BS]`, `Cer[NDS]`, `Cer[ADS]`, `Cer[BDS]`, `Cer[AP]` and `Cer[NP]` exist
only as `[M-H]-`.

In skin that showed up as the largest unclaimed signal in the whole dataset. Screening the
unidentified positive-mode features for homologous series
(`scripts/find_unclaimed_series.py`) put a 16-member ladder at the top — m/z 598.58 to 808.81,
**total area 6.2×10⁹**, in up to 61 of 63 files — and matching its members against the neutral
masses of the negative-mode library named them `Cer[ADS]`, `Cer[AS]` and `Cer[NP]`.

Reading their spectra settled it more precisely. Every member is a **phytoceramide with a
long-chain base**:

| precursor | sphingoid-base ions | base |
|---|---|---|
| 654.64 | 286.2754 / 268.2646 / 250.2536 | t17:0 |
| 668.66 | 300.2907 / 282.2798 / 264.2691 | t18:0 |
| 696.69 | 328.3219 / 310.3110 / 292.3005 | t20:0 |
| 724.72 | 356.3528 / 338.3424 / 320.3320 | t22:0 |

Three ions 18.011 apart each time — a base losing water three times. Textbook skin.

## Why positive mode is not just more coverage

Negative mode **cannot** separate `Cer[ADS] d18:0_24:0` from `Cer[NP] t18:0_24:0`: same formula,
same exact mass to five decimals, and negative-mode spectra carry nothing that distinguishes them
(`results/skin_cer_ads_check.md` — the related `Cer[ADS]`/`Cer[AS]` call comes down to 1% of the
dot product).

Positive mode separates them, because **only a three-hydroxyl base can lose a third water**:

    Cer[ADS] d18:0_24:0    266.2842   284.2948          (two hydroxyls, two losses)
    Cer[NP]  t18:0_24:0    264.2686   282.2791  300.2897  (three)

That is the reason to build this, and it is why the library ended up covering only part of what it
started as.

## The 85% that had to be dropped

The first build emitted 3,592 entries — every molecule in the seven missing classes. It **cost 11
correct identifications**: `Cer[NS] d31:1` became `Cer[NDS] d31:1`, `Cer[NS] d34:3` became
`Cer[ADS] d34:2`, and recovery against the reference fell from 98.5% to 96.6%.

The cause is that LipiDex enumerates these classes over a **shared base/acyl grid**, so a class
label does not constrain the sphingoid base the way the notation suggests. `Cer[NDS]` entries exist
with unsaturated bases, and their fragments are then identical to a `Cer[NS]` entry of the same
base. Where two entries are identical in precursor *and* fragments, the winner is decided by
library load order — which is to say by nothing.

Two filters, applied in order:

**1. Exact degeneracy.** Any new entry whose (precursor, fragment masses) already exists in a
positive-mode entry is dropped: it carries no information the library lacks, so it can only take
matches from an entry that is already validated. **3,069 of 3,592 went this way.** Recovery
recovered to 97.1%.

**2. Two-hydroxyl bases entirely.** The 523 survivors were not surviving on merit. LipiDex's
positive `Cer[NS]` set has 1,200 entries against 2,406 negative, so it has enumeration *gaps*, and
a surviving `d`-base entry wins by filling a hole in the grid rather than by matching evidence.
Emitting them still cost 8 correct calls. `t` bases are different in kind: nothing else in the
library can produce a third water loss, so they add a discriminator instead of competing for one.

**446 entries survive**, all phyto. `--all-bases` overrides, and the flag exists so the measurement
can be repeated rather than taken on trust.

## Result

| | before | after |
|---|---|---|
| skin positive, recovery | 98.5% (9 missing) | **98.5% (9 missing)** |
| `Cer[NP]` rows | 0 | **33** |
| `Cer[AP]` rows | 0 | **14** |
| `Cer[NS]` rows | 57 | **129** |
| total ceramide area | — | 3.17×10¹⁰ |

**No regression, and a class that was previously invisible in this polarity is now reported.** The
`Cer[NS]` rise is a side effect worth noting: with the phytoceramides claiming their own spectra,
sphingosine ceramides stop being assigned to them.

## Fragments and intensities

LipiDex's own positive `Cer[NS]` scheme is reused — precursor water loss plus three sphingoid-base
ions, with annotations in its format so the purity calculation can resolve the base the same way it
resolves acyl chains elsewhere.

Intensities are split by base type, and the split is the honest part:

- **`t` bases were measured here**, pooled over six files, on the very series this library exists
  to claim: `-2H2O` 999, `-3H2O` 950, `-H2O` 600, precursor `-H2O` 400.
- **`d` bases keep LipiDex's numbers unchanged.** This dataset contains none to measure, and
  inventing them would be worse than inheriting values already validated to exactness. They are
  not emitted by default in any case.

## Glucosylceramides, added the same way and mostly ruled out the same way

`GlcCer[NS]`, `GlcCer[NDS]` and `GlcCer[AP]` are negative-only too, so the same builder covers
them. Two of the three fall out immediately:

- **`GlcCer[NS]` and `GlcCer[NDS]` carry `d` bases only** (1,806 entries each), and `HexCer[NS]`
  already has 600 positive entries covering the same molecules under a different class name. The
  phyto restriction removes them without needing a special case, which is a small vindication of
  the rule.
- **`GlcCer[AP]` is the genuinely absent one**: 215 entries, every one a phyto base, with no
  `HexCer[AP]` anywhere to duplicate. It is the only glucosylceramide class the build emits.

Hexosyl entries get two extra peaks, taken from LipiDex's own positive `HexCer[NS]` entries rather
than assumed — the sugar leaving as the free hexose (180.063) and as hexose plus water (198.074).
The base fragments then sit on the ceramide left behind.

### The result is a clean negative

**Not one `GlcCer[AP]` was identified in either set**, and recovery was unchanged — skin positive
98.5%, rat heart positive 98.8%. The rat heart set reports no phytoceramides at all (38 `Cer[NS]`,
1 `HexCer[NS]`), which is the expected contrast: phyto bases are a skin lipid. That
needed checking rather than shrugging at, because their masses land within the 0.01 Da search
window of abundant phosphatidylcholines — 356 features at `GlcCer[AP]` masses in skin positive, the
largest being `PC 36:2` at 1.6×10¹⁰.

Rescoring every PC-assigned spectrum in one file against the whole library:

| | |
|---|---|
| PC-assigned spectra with a `GlcCer[AP]` candidate at all | 18 |
| median margin, winner over `GlcCer` | **884** points |
| times `GlcCer` came within 50 points | 2 |
| times `GlcCer` won | **0** |

Both near misses are on spectra whose *winner* scores 35 and 37 — far below the dot-product floor
of 500, so neither is reported by anything. **The class is absent from this tissue, not crowded out
of it.** The entries stay: measured risk, not assumed safety.

## Two bugs the glucosylceramide pass exposed

**The builder was reading its own output.** Once `Ceramides_Positive` joined `DEFAULT_LIBRARIES`,
the degeneracy scan found the previous build and declared every entry degenerate with itself,
silently cutting 446 entries to 98 with no error anywhere. A generator that consumes the library
set it writes into has to exclude its own output path, and nothing about the failure looked like a
failure.

**`GlcCer[AP]` was being dropped entirely.** The builder recovered each neutral mass by assuming
`[M-H]-`, and that class is enumerated only as `[M+FA-H]-`. It now reads whichever negative adduct
a class actually uses. The symptom was a silent zero — the class simply never appeared in the
output counts.

## A positive-only library is scoped to positive runs

`Ceramides_Positive` is **not in the negative-mode configurations**. LipiDex matches a library
entry to a spectrum on **precursor mass alone, with no polarity check** — a quirk of the original,
reproduced faithfully (docs/LIPIDEX_ALGORITHM.md) — so an `[M+H]+` entry can in principle score
against a negative-mode spectrum. Excluding a polarity-specific library from the runs of the other
polarity costs nothing and removes the possibility.

**A false alarm on the way there, recorded because the mistake is easy to repeat.** The skin
negative results contain 70 `Cer[AP]` / `Cer[NP]` / `GlcCer[AP]` rows, and they looked at first
like exactly the harm above. They are not: checking the `Library` column shows every one came from
`LipidBlast_Formic.msp`, LipiDex's own **negative-mode** entries — the shipped set has 750 negative
`Cer[AP]`, 750 negative `Cer[NP]` and 215 negative `GlcCer[AP]`. They are ordinary, correct
identifications that predate this work entirely.

**The class name alone does not say which library or which polarity produced a row.** The `Library`
column does, and it is the thing to check before concluding that a new library has broken
something.

Measured effect of the scoping: skin negative moves by one molecule, 99.4% either way.
`RunConfig.libraries` is per-run, so this is configuration rather than code, and the same scoping
already applies in the other direction to `fatty_acid_libraries`.

## Limits

- **`Cer[AS]` and `Cer[BS]`** — α- versus β-hydroxy acyl — share a formula and a base, so they stay
  degenerate here exactly as in negative mode.
- **The base carbon number is determined; the acyl split follows by subtraction.** A row named
  `Cer[NP] t18:0_24:0` rests on the base fragment, not on an acyl fragment.
- **Measured on one tissue.** The `t`-base intensities come from skin, where phytoceramides are
  abundant. Another matrix may fragment them differently, and the numbers should be re-measured
  before being relied on elsewhere.
- **Glucosyl- and hexosyl-ceramides were not covered.** `GlcCer[NS]`, `GlcCer[NDS]` and
  `GlcCer[AP]` are equally negative-only and would need the hexose loss modelled first.
