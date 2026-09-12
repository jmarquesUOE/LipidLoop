# Retention modelling within a lipid class

The criterion LipiDex 1 never uses. Implemented in `rtls.py`, on by default as a filter, off by
default as a way of naming things.

⚠ **Not "RTLS".** An earlier version of this file, like an earlier version of `rtls.py`'s own
docstring, credited this to LipiDex 2's "RTLS" — wrong. RTLS is LipiDex 2's *real-time library
search*, an acquisition-time method that triggers class-targeted MSn; it has nothing to do with
retention modelling. The peak-finder parameter that does the same thing this module does is
**`UseRTLS`** — *"creates a retention time model within each lipid class based on fatty acyl
carbons and fatty acyl unsaturations to assess identification likelihood"* (read from an
installed LipiDex2 v1.1.2 copy's help text, see `docs/LIPIDEX2.md`) — a confusingly similar name
for an unrelated feature. This file keeps calling the module `rtls.py`/"RTLS" as our own internal
shorthand below; it is not a claim about LipiDex 2's naming.

## The model

On a reversed-phase column a lipid's retention is largely set by two numbers — acyl carbons and
double bonds. More carbons is more hydrophobic and elutes later; more double bonds kinks the
chain and elutes earlier. Within one class that is close to linear:

    retention ≈ intercept + per_carbon × carbons + per_double_bond × double_bonds

Fitted on the CKD rat heart set the coefficients come out chemically sensible wherever the fit
is good: **+0.14 to +0.42 min per carbon, −0.18 to −0.76 min per double bond**.

This is a sharper question than LipiDex 1's `checkClassRTDist`, which only asks whether a lipid
falls in the window where its class elutes. This asks whether it falls where *that member* of
the class should.

## Two things that make it work

**Fit iteratively.** The data being fitted contains the outliers being looked for, so a handful
of bad rows drag the surface onto themselves and nothing looks wrong. Fit, drop residuals beyond
3σ, refit.

**Fit on the well-detected rows, then judge everything.** Fitting on all rows is circular: the
duplicated identifications that most need judging are exactly the ones that spoil the fit meant
to judge them. PC on this data fits at R² 0.70 over all rows — below the usability bar, so PC
would never be judged — and at **R² 0.87** when fitted only on compounds found in at least half
the samples. That one change is the difference between the filter doing nothing for PC and
doing most of its work there.

**A class that does not fit is left alone.** A bad fit is an absence of evidence, not evidence of
absence. A model is used only with ≥6 points, R² ≥ 0.7, residual SD ≤ 1 min, and coefficients of
the chemically correct sign.

## The fit is a diagnostic in itself

| class | R² | residual SD | per C | per DB |
|---|---|---|---|---|
| Alkenyl-TG P- | 0.980 | 0.09 min | +0.176 | −0.339 |
| AC | 0.936 | 0.06 | +0.143 | −0.180 |
| Cer[NS] d | 0.924 | 0.32 | +0.415 | −0.762 |
| TG | 0.865 | 0.47 | +0.250 | −0.417 |
| **PS** | **0.186** | **3.17** | **−0.535** | **+1.283** |

PS fits at R² 0.19, with a residual spread of three minutes and coefficients of the **wrong
sign** — retention falling with carbon number and rising with unsaturation. That is not a failure
of the model. It is the model saying those PS identifications are not a homologous series, which
agrees with PS appearing at eleven different retention times from 6.1 to 24.4 min.

## As a filter — on by default

Removes identifications sitting more than `max(3σ, 1 min)` from their class's surface.

| | RTLS off | RTLS on | reference |
|---|---|---|---|
| rows | 4,946 | 4,870 | 4,128 |
| identified rows | 984 | **908** | 642 |
| distinct lipids | 530 | 529 | 444 |
| **lipids reported >3 min apart** | **13** | **3** | 2 |
| reference molecules recovered | 416 | **416** | — |
| molecules only ours | 82 | 81 | — |

76 identified rows removed, **at no cost to reference recovery at all**. `PC 38:4` goes from
13 rows spanning 6.5–11.8 min to **5 rows spanning 8.7–9.7**. Removals concentrate where the
duplication was: PC 51, PE 12, TG 10, SM 2, Cer[NS] 1.

This is what finally closed the duplicate-retention-time problem — the redundant-ID merge could
not, because it only collapses rows within about 1.2 peak widths of each other, and these were
minutes apart.

## As a way of naming features — off by default

`retention_model_extend` predicts where any library member of a modelled class should elute and
names unidentified features that match on accurate mass and predicted retention. Ambiguity is
declined rather than guessed: if two candidates fit, neither is assigned. Assignments are kept
in a separate field and reported in an `Identification Source` column as `RT model` rather than
`MS2`, and the column only appears when the option is on, so the default output keeps LipiDex's
schema exactly.

**Measured, it is not good enough to turn on by default.** On this data it named 186 features
and added 89 distinct molecules — of which **exactly one** was a molecule the reference also
reports. The other 88 are unverifiable from this data.

That is the expected shape of the evidence rather than a bug. An assignment from accurate mass
plus retention has no fragmentation behind it: it cannot distinguish isomers, and it cannot tell
a lipid from anything else of the same formula eluting at the same moment. It is a hypothesis
generator, useful when you want candidates to go back and target with MS2, and it should not be
carried into a quantitative result as though it were an identification.

## What the model cannot do: isomers

`parse_sum_name` collapses a molecular name before judging it, so `TG 16:0_18:1_18:2` is judged as
`TG 52:3`. This is deliberate — without the collapse the model silently declines to judge every row
reported at chain resolution, 78 identifications on the skin set — but it has a consequence that an
earlier version of this document got wrong.

**`PC 18:0_20:3` and `PC 18:1_20:2` both parse to (PC, 38, 3), so the model issues one prediction
for every isomer of a composition.** It cannot rank them, cannot say which is which, and cannot say
that either is implausible.

The earlier justification for collapsing was that rows sharing a sum "elute within seconds of each
other". They do not. Over 131 duplicated names in one positive-mode study the gap between rows
sharing a name has a **median of 0.41 min; only 2% are within 0.1 min and 36% exceed half a
minute.** Collapsing is still right, for the reason above, but it is a limitation accepted rather
than a distinction without a difference.

### The consequence for judging rows

At most one isomer can sit on the single prediction, so the rest are displaced *by construction*.
Measured within duplicate sets:

| | distance from the prediction |
|---|---|
| names appearing once | 0.112 min (median) |
| nearest row of a duplicate set | 0.099 min |
| its siblings | 0.391 min |

So the nearest row of a set behaves like a unique name, and the others do not. **Judge a row against
this model only when its name appears once.** A verdict of "off-model" on a duplicated row meant
only that it was not the isomer closest to the class average; it was withdrawn from
`Duplicate Verdict` for that reason (see `peakfinder.UNRELIABLE_CV`).

### Two scores, and why one is not enough

Rows carry `RT Model Error` in minutes, `RT Model Z` in standard deviations of their own class
model, and `RT Model R2`.

**Minutes are not comparable between classes.** Half a minute is unremarkable in TG, whose model
scatters by 0.36 min, and large in LysoPC at 0.05. A fixed threshold in minutes is therefore a
different test in every class, which is why one in minutes flagged 11.5% of positive rows against
2.9% of negative — an artefact of the classes present, not a finding about the rows.

**R² says whether the class model is worth believing at all.** A deviation of three standard
deviations against an R² of 0.19 means very little, and PS fits at 0.19 on real data with
coefficients of the wrong sign. Reported beside the z so it is not read as though every class were
LysoPC at 0.992.

### The one thing the set can still say

If even the row *nearest* the prediction is far from it, no arrangement of isomers explains that,
and the class model or the identification is wrong. A name whose smallest |z| exceeds
`SET_OFF_MODEL_Z` (3.0) is flagged `class model?` on every one of its rows. The threshold is
legitimate on this quantity and on no other in a duplicate set, precisely because the nearest row
was measured to behave like a unique name.

A set is judged only when *every* row carries a z; otherwise the verdict would be a claim about
whichever rows the model happened to reach, wearing the name of the whole set.

⚠ **Untested.** It fires on nothing in either polarity of the study it was written against — 108
duplicated names in positive, 5 in negative, none flagged, worst set 1.65 σ against a threshold of
3. A criterion that has never fired has also never been shown to fire when it should.
