# Free fatty acids: the class a spectral library cannot reach

No LipiDex library contains a single free fatty acid. Across 91 fatty-acid masses the default set
has **zero** negative-mode entries, so a free fatty acid could not be identified at all — while in
the skin organoid data they are among the largest peaks in the file. `FA 18:1` is a 1.9e8 feature
present in **63 of 63 files**, and every species from C14 to C30 is there, all unidentified.

## Do they fragment? Barely

Pooling every MS2 taken at a fatty-acid precursor across 16 files and asking which product ions
recur, rather than reading one spectrum:

| | MS2 spectra | `[M−H−CO₂]⁻` present | relative intensity |
|---|---|---|---|
| FA 16:0 | 711 | **0%** | — |
| FA 18:0 | 561 | **0%** | — |
| FA 18:1 | 305 | **0%** | — |
| FA 24:0 | 21 | **0%** | — |
| FA 26:1 | 18 | **0%** | — |
| **FA 20:3** | 11 | **91%** | 999 |
| **FA 20:4** | 15 | **100%** | 999 |
| **FA 22:6** | 6 | **100%** | 927 |
| FA 22:5 | 12 | 17% | 984 |

The carboxylate anion is stable. Saturated and monounsaturated acids produce **nothing
reproducible** at these collision energies — and what looks like a spectrum is co-isolated
background, because with no real fragment to normalise against, whatever else is in the isolation
window becomes the base peak. The apparent "fragments" of `FA 16:0` are 170.9862, 186.9283,
214.9931: mass defects that no organic ion can have. They are the polymer and salt background
already documented in [EXCLUSION_LIST.md](EXCLUSION_LIST.md).

**Only polyunsaturated acids lose CO₂**, from about three double bonds up, and when they do it is
the base peak.

## So the library contains only what fragments

`data/libraries/FreeFattyAcids_Negative.msp` — 63 entries, C14–C30, **DB ≥ 3 only**, `[M−H]⁻` with
one peak: `[M−H−CO₂]⁻`. Built by `scripts/build_fatty_acid_library.py`.

Saturated and monounsaturated acids are **deliberately absent**. Writing entries for them would be
writing a library of noise, and each would match any spectrum at its mass. This is the same
decision as the negative-mode ganglioside library: write what the evidence supports and stop.

The enumeration uses a methylene-interrupted ceiling — each double bond after the first needs three
more carbons, `1 + (C−3)//3` — so `FA 20:5` and `FA 22:6` are included and `FA 14:6`, whose mass
would otherwise be free to claim a real feature, is not.

## The rest are named by retention, and that is not a compromise

One unidentified feature at m/z 281.2486 could be anything. **Forty-two features sitting on
forty-two exact fatty-acid masses, whose retention times jointly obey `RT ~ C + DB`, could not be a
coincidence.** A chromatographic law is being obeyed, and obeying it is the evidence.

`fatty_acids.py`, on by default (`annotate_fatty_acids`), negative mode only. On the skin data:

    free fatty acids: FA n=42 R2=0.969 sd=0.390 min  +0.096..+0.710/C (curved) -0.650/DB
    free fatty acids: 48 features named within 1.17 min of the surface

56 rows in all: **48 named from retention, 8 confirmed by MS2** through the CO2-loss library, and
the `Identification Source` column says which is which on every row.

The resulting ladder, C14 to C30, one member per carbon, with the monounsaturated series shifted
earlier throughout:

| | 14:0 | 16:0 | 18:0 | 20:0 | 22:0 | 24:0 | 26:0 | 28:0 | 30:0 |
|---|---|---|---|---|---|---|---|---|---|
| RT | 9.87 | 10.32 | 10.84 | 11.47 | 12.26 | 13.29 | 14.56 | 16.00 | 17.59 |

### It is tested as a set, never per feature

Each candidate on its own is only a mass, and a mass is not an identification. Every gate is on the
*series*:

- `MIN_CANDIDATES` (8) — a handful of masses always land somewhere; a series needs members;
- `MIN_R2` (0.85) and `MAX_RESIDUAL_SD` (0.75 min) — the fit has to be good, or no law is being
  obeyed;
- **retention must rise with carbons and fall with double bonds** — the signs are the chemistry.
  Reversed phase retains longer chains and releases unsaturated ones. A model with the wrong signs
  has fitted noise, and no number of points rescues it. The rise is checked on the *slope* across
  the fitted carbon range, not on a coefficient — see Curvature below for why that distinction
  matters;
- `MAX_SIGMA` (3) — a member further than this off the fitted surface is not named, mass match or
  no mass match;
- a candidate seen in fewer than 3 files may be named but does **not** join the fit, so a scatter
  of one-off features cannot set the surface that judges them.

**If any gate fails, nothing is named at all.**

Names are written with `Identification Source` = `RT model`, so a reader can always tell which rows
have no spectrum behind them.

### The hold-out test

The eight species the library *can* confirm by fragmentation — `FA 20:3`, `20:4`, `20:5`, `22:3`,
`22:4`, `22:5`, `22:6`, `24:4`, with dot products 726–992 — are an independent check. Fitting the
surface on **those eight alone** and predicting the other 48:

| | median error | 90th percentile | max |
|---|---|---|---|
| linear | **0.33 min** | 2.82 | 4.72 |

Eight seeds spanning C20–C24 with 3–6 double bonds predict species from C14 to C30 with a median
error of twenty seconds. **The species identified by fragmentation land on the same surface as the
species named by retention.** That is the validation of the method, and it is the argument for
trusting the 48 that no spectrum could ever have reached.

### The curvature is the gradient, not the chemistry

The whole fatty-acid series elutes **inside the IPA ramp**. On this method (`M2024`) channel C
starts at 10.0 min and reaches 95% at 22; the fatty acids run 9.57 to 17.59 min, so **82% of the
class elutes while the solvent is still changing**, and the rest is pressed up against the moment
it starts.

That shows in the data. The ten species eluting before 10 min are squeezed into a 0.4 min band —
`FA 14:1` 9.57, `FA 20:5` 9.67, `FA 16:2` 9.73, `FA 22:6` 9.80, `FA 18:3` 9.81, `FA 14:0` 9.87,
`FA 20:4` 9.91 — and their residuals are systematically positive for the polyunsaturates
(+0.64, +0.53, +0.39) and negative for the short saturates (−0.29, −0.10). **The polyunsaturates
are being held until the eluent gets strong enough, and they arrive together.** Double-bond
resolution is poor at that end of the class for a chromatographic reason, not an analytical one.

**So the C² term is an empirical patch on a gradient-shaped effect, and it will not transfer to
another method.** A segmented fit was tried, on the expectation of a kink at the IPA onset, and
there is no kink: the best breakpoint lands at C24 (12.78 min), mid-ramp, and buys only 0.02 min of
residual sd over the quadratic. The ramp is gradual, so its effect is gradual.

| model | residual sd | R² |
|---|---|---|
| linear in C | 0.530 | 0.9369 |
| quadratic in C | 0.393 | 0.9660 |
| segmented, break fitted at C24 | 0.371 | 0.9697 |

This is the argument for anchoring on standards rather than refining the polynomial.

### Hybrid: anchoring on authentic standards

`fatty_acid_standards` takes a CSV of measured retention times, `name,retention`:

    FA 16:0,10.32
    FA 18:1,10.36
    FA 20:4,9.91

Those points join the fit and are **never trimmed**. A point whose identity is known cannot be an
outlier from a surface meant to describe it — if a standard disagrees with the surface, the surface
is wrong, so it stays in and drags the fit rather than being quietly discarded to make the fit look
good. The observed series still extends the surface into the compositions no standard covers.

Ten standards spread across C14–C30, including two or three either side of the IPA onset, would pin
the shape far better than fifty observations that all have to be taken on trust — and would say
directly whether the curvature is the gradient or the chemistry.

### The list is separate on purpose

`fatty_acid_libraries` is its own configuration field, not appended to `libraries`. Fatty acids are
a different kind of evidence: three quarters of the class produces no fragment at all, and the
retention surface belongs to one column and one gradient. Keeping the list separate is what lets it
be swapped for a standards-built library without touching the lipid libraries, and what keeps the
load order of those libraries — which breaks score ties — undisturbed.

### Curvature, mechanically

Retention against carbon number is not straight: the step between consecutive saturated acids grows
monotonically from 0.21 min (C14→C15) to 0.81 min (C29→C30). A line leaves systematic residuals and
therefore a wider naming window than the data deserves — **sd 0.530 min linear against 0.393 with a
C² term**, tightening the 3σ window from 1.63 to 1.18 min.

The quadratic is fitted only when there are at least 12 points to determine it, kept only if it cuts
the residual sd by at least 10%, and never fitted at all when the line already predicts retention to
better than 0.02 min. The residual sd sets the window inside which a feature gets named, so a
spuriously small one is the dangerous failure — not a slightly wide one.

**The sign gate had to be rewritten for this, and it caught itself doing the wrong thing.** With a
C² term the bare linear coefficient is no longer the per-carbon slope: the fitted surface here is
`−0.441/C +0.0192/C²`, whose slope at C20 is **+0.327** — entirely correct — while the linear
coefficient alone is negative. Testing that coefficient rejected a good model and named nothing at
all. What must be positive is the *derivative*, checked across the whole carbon range the model was
fitted over; being linear in C, the two ends decide it. The failure was safe in the right
direction, which is the point of the gate.

A second bug surfaced with it: when a fit is very good the residual spread collapses to
floating-point noise, `3 × spread` becomes ~1e-15, and almost every point is trimmed as an outlier
from a model it fits exactly — leaving a handful spanning a fraction of the carbon range, over
which the sign gates are then checked. The trimming threshold now has an absolute floor of
0.05 min.

## The negative control: CKD rat heart

Running the same code on the other validation set is the strongest evidence the gates work, because
there it **declines**:

    free fatty acids: 5 candidate masses, no model could be fitted — nothing named

Rat heart is not skin. Side by side, on features sitting within 0.01 Da of a fatty-acid mass in
negative mode:

| | skin organoid | CKD rat heart |
|---|---|---|
| features on fatty-acid masses | 61 | 22 |
| **distinct compositions** | **57** | **8** |
| compositions with more than one peak | 2 | **5** |
| kept | 56 | 8 |
| already flagged `In-source fragment` | 3 | **10** |
| already flagged `Adduct of existing peak` | 0 | **4** |

In skin the masses come as 57 compositions with one peak each: a homologous series. In rat heart
they come as **8 compositions, five of them at several retention times**, and two thirds of the
features are already flagged as in-source fragments or adducts of peaks that are in the table
anyway. `FA 18:0` alone appears at 7.06, 7.31, 9.76, 10.58 and 13.55 min — that is not a compound,
it is the fatty acyl anion falling off five different phospholipids in the source.

Five eligible candidates is below `MIN_CANDIDATES`, so no surface was fitted and nothing was named.
**A method that names 56 fatty acids in the tissue that has them and none in the tissue that does
not is behaving correctly.** Skin makes free and very-long-chain fatty acids; heart, extracted this
way, does not.

It is also a reminder of what the in-source filter is doing upstream. Without it, those 14 rat
heart features would have been eligible candidates, and 19 candidates across 8 compositions might
well have fitted *something*.

## Limits, plainly

- **Double-bond position and geometry are not determined.** `FA 18:1` is oleate *or* vaccenate; the
  model separates compositions, not isomers, and HCD would not separate them either.
- **The model is method-specific.** Every coefficient here belongs to this column and gradient.
  It is refitted per run, as it must be.
- **Odd-chain polyenes are excluded from `candidate_masses()` entirely.** `max_double_bonds()` caps
  odd carbon counts at 1 double bond — mammalian odd-chain fatty acids come from propionyl-CoA
  priming, not the acetyl-CoA elongation/desaturation pathway that builds even-chain polyenes, so
  species like `FA 13:3`, `15:4`, `19:5` are never real. Found live on the NIST SRM 1950
  method-development batch: `FA 15:4` was a genuine, reproducible, blank-clear peak — real signal,
  impossible chemistry — that the old ceiling would have offered a name.
- **Very-long-chain species still deserve a look before use.** `FA 30:1` at 4.6e7 in 59 files is
  plausible in skin, which genuinely makes very-long-chain fatty acids, but it is named on
  retention alone at the extrapolated end of the surface.
- **An in-source fragment that lands on a fatty-acid mass at the right retention would be named.**
  Retention is one orthogonal axis, not a proof — the same caveat as everywhere else in this
  pipeline, and the reason `Identification Source` exists.
