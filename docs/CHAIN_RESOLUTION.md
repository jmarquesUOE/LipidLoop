# Chain resolution, and what a chain-resolved name is worth

A row reading `PC 34:1` reports a **sum composition**: 34 carbons and one double bond across the
chains, without saying how they are divided. A row reading `PE 18:2_18:2` reports a **molecular
species**: the chains are named, though not their *sn* positions — the underscore is deliberate and
a slash would claim positions the spectra do not support.

This document is about how often the second happens, why, and what it does and does not license.

## It is a property of the polarity, not of the pipeline

How often a table names chains is read as a quality of the software. It is not.

| | rows resolved to chains |
|---|---|
| positive | **6%** |
| negative | **53%** |

The cause is the ion. A choline lipid in negative mode is a formate adduct and fragments through
the loss of formate rather than into carboxylates, so negative resolves **0 of 29 SM, 0 of 15
LysoPC, 0 of 1 PC**. Positive cannot resolve them either. The classes negative *does* resolve — PE
51%, PI 78% — are the ones the polarity rule keeps in negative anyway, while positive's copies of
those classes resolve **0 of 79 PE**.

⇒ **The two polarities resolve chains for disjoint sets of classes, and each already keeps the
classes it resolves.** Carrying molecular names across from the discarded copies was implemented on
the expectation that choline lipids would gain chain detail from negative. They do not. It
transfers ten rows of one class (Plasmenyl-PC) and is harmless, but the premise was wrong.

⚠ Any chain-resolved PC in a historical table came from a **negative-mode spectrum**, not from
better processing.

## ⚠ A chain-resolved name is not necessarily chain evidence

**`Purity` is a score on the winning identification, not a gate on the name it carries.** The name
comes from whichever library entry won the dot product, and library entries are frequently written
at chain resolution. So a row can be reported at chain resolution on no chain evidence at all, and
`Purity = 0` beside a molecular name is exactly that case.

Purity is 0 by design for entries outside their **optimal polarity** (see `purity.py`), so a zero
is not automatically a warning — but it does mean the chains were not established by the fragments
in this run. ⚠ There is **no** by-design zero for LipidBlast, contrary to what this document and
`purity.py` both said until two independent reviews measured it.

**Check `Purity` before reading chains off any row.** In one negative-mode table, 3 of the 8
compositions carrying more than one chain combination were ceramide pairs whose minor member scored
zero purity, with a dot product about 1% behind the major (974 against 964). Those are near-ties in
the library resolved arbitrarily, not two species measured.

## Can the chain distribution be measured? Not from the result table

The natural question — for `PC 38:3`, how much is `18:0_20:3` and how much `18:1_20:2` — cannot be
answered from these tables, and the obstacle is structural rather than one of sample size.

Of 155 negative-mode compositions carrying a chain-resolved name, **8 carry two different
combinations**, all of them pairs. Five are glycerophospholipids with purity 100 on both members,
and they agree with known biology:

| composition | combinations | split |
|---|---|---|
| PI 36:4 | `16:0_20:4` / `18:2_18:2` | 95 / 5 |
| Plasmenyl-PE P-38:4 | `P-18:0_20:4` / `P-16:0_22:4` | 90 / 10 |
| PG 38:5 | `18:1_20:4` / `18:2_20:3` | 79 / 21 |
| PI 38:6 | `16:0_22:6` / `18:2_20:4` | 75 / 25 |
| PE 36:4 | `16:0_20:4` / `18:2_18:2` | 58 / 42 |

`PI 36:4` at 95% `16:0_20:4` is the canonical arachidonate-at-*sn*-2 species, which is
encouraging — and five compositions is not a distribution.

⚠⚠ **The selection runs the wrong way.** To be seen at all, a minor isomer must (1) separate
chromatographically, (2) trigger MS2, and (3) win a dot product. A species present at 5% fails all
three. So near-equal splits are observed and lopsided ones are reported as a single peak — **the
bias is toward the splits that matter least**, which is the opposite of what a prior needs.

### Where the answer actually lives — `Spectral_Components.csv`

Purity is built from **moiety fragments** — the acyl and sphingoid fragments that carry chain
identity — and produces a per-name breakdown, the `Spectral Components` of
`calculate_weighted_purity`. That breakdown is the within-peak chain distribution, and it does not
require chromatographic separation, which is the binding constraint above. It is now written
(`components.py`), one row per candidate per peak:

| column | |
|---|---|
| `Compound Group` | joins to every other table |
| `Component` | the candidate's name, **adduct stripped** |
| `Share of Fragment Signal (%)` | its share of the chain-bearing signal in this peak |
| `Isomer` | `yes` if it collapses to the same sum composition as the reported name |
| `Rank` | 1 is the largest share |

On a four-injection positive QC set: **1,642 components over 286 peaks**, against 8 compositions
obtainable from the results table.

`Isomer` separates two questions that share one breakdown. A component with the same sum
composition is a chain isomer, and its share says how this lipid's chains divide; a component from
another class says what else was co-isolated. Summing the two would turn a co-isolation into a
chain composition.

#### Three traps, all of them load-bearing

**The purity denominator is the wrong one.** `purity_weight` counts the Gaussian score once per
*entry*, so it grows with the number of candidates: shares under it summed to a median of **41%**
and to less the more candidates a peak carried, making two peaks incomparable. That is correct for
purity, whose question is how much belongs to the winner, and wrong for a distribution. Shares are
normalised within the peak instead and sum to 100.

**One combination appears under several adducts.** `TG 18:1_18:1_22:4 [M+NH4]+` and the same chains
as `[M+Na]+` are two library entries for one molecule; left separate they read as two isomers at
half the share each. Merged on the chains, then re-ranked. This removed 289 spurious components
from the run above.

**It is worth least exactly where the library enumerates most.** For triacylglycerols the library
holds nearly every arithmetic combination, so `TG 58:6` returns 27 components down to
`TG 10:2_24:1_24:1` — the combinatorial library re-expressed, not a measurement. **Read the
component count as a warning**: a peak with 20 of them has been spread across a catalogue, not
resolved into a distribution.

#### What it does not yet answer

⚠ **Purity is only computed in a class's optimal polarity**, so a positive-mode run yields
breakdowns for neutral lipids and essentially nothing else — the positive QC set above is 261 TG
peaks, 19 Alkenyl-TG, and no glycerophospholipids at all. **The two-chain classes this file answers
best (PE, PI, PG, PS) are negative-mode, and are not yet demonstrated on real data.**

⚠ **Fragment intensity is not molar ratio.** The *sn*-2 carboxylate is preferentially released in
glycerophospholipids and response varies with chain length and unsaturation, so a share is a ratio
requiring calibration against authentic standards, not a composition.

⚠ **Absence is not evidence of absence.** A candidate reaches the breakdown only if *all* of its
chains were found: `purity.py` returns nothing for a candidate whose second chain fragment fell
below the 5% intensity floor, so it is missing from the denominator rather than scored low.

## Enumeration is not a substitute

Enumerating the arithmetically possible combinations of a sum composition is easy and nearly
useless without a biological constraint. Over the chain set actually used in the shipped library
names — 78 distinct chains, every carbon number 10 to 26 including all odd ones, up to 7 double
bonds — a sum composition yields:

| composition | arithmetic combinations | peaks observed |
|---|---|---|
| PC 34:1 | 11 | — |
| PC 38:3 | 19 | 4 |
| TG 52:3 | **173** | — |

The first candidates for `PC 38:3` include `12:0_26:3` and `14:1_24:2`. **The library's chain list
is itself combinatorial and was never a statement about what exists in mammals**, so it cannot
supply the plausibility filter. That filter is a question for the literature, not for the data.

## Related

- `docs/RTLS.md` — why the retention model cannot rank isomers of one composition
- `docs/SPLIT_PEAKS.md` — when several rows of one name are one molecule rather than isomers
- `manuscript/DRAFT.md` §7.1, §8.6, §11.2

## The `Chain Evidence` column

Every chain-resolved row now says where its chains came from:

| value | meaning |
|---|---|
| `fragments` | the chains were tested against moiety fragments in this run |
| `library name` | inherited from the library entry that won the dot product — **never tested** |
| *(blank)* | the row reports a sum composition, so it asserts no chains to qualify |

The test is **whether any component of the purity breakdown names THIS row's lipid** — not whether a
breakdown exists.

⚠⚠ **A high purity does not mean the purity belongs to the row it is printed on.**
`purity.py::calc_purity` drops a top hit that scores nothing and then returns `purities[0]`, which
by then belongs to a different candidate, and `search.py` writes it onto the winner's row. The
clearest case in the delivered data is **`Cer[ADS] d17:0_17:0` at `Purity = 100`** whose entire
breakdown reads `Cer[NP] t18:0_16:0 (100)`. The phyto ceramide was measured; the dihydro one was
named. Across rank-1 matches this affects **431 positive and 341 negative** spectra.

Measured on one negative-mode study, chain-resolved rows: **139 `fragments`, 27 `library name`**
(Cer[ADS] 15, Cer[NDS] 12) of 166. Positive: 48 and 0.

⚠ **They are not weak matches.** The lowest dot product among them is **699**, against **511** for
rows that scored — so raising a dot-product cutoff removes good identifications and leaves these in
place. It is a provenance problem, and the column is the fix.

⚠ **Three explanations of this were wrong before one was measured, two of them published here.**
First: `purity.py` documents that LipidBlast entries score 0 *by design* via `if not ls.is_lipidex`.
**That gate never fires** — `lib_gen/Lipid.java` writes `Type=LipiDex` as an unconditional string
literal, so every generated library carries it. Second, and written in this file as the correction:
that LipidBlast peaks are annotated generically, leaving `transition_type` unset. **Also false** —
**97.8% of LipidBlast entries (53,965 of 55,169) can set it**, because annotations such as
`C1O1_Sphingoid Fragment_[d14:1]` put a *space* before `Fragment` and so survive a filter that
excludes `_Fragment`. LipidBlast is not merely eligible, it is the dominant source: it wins 278 of
the 341 mis-attributed negative matches.

**The truth is simpler than either story: there is no by-design zero for LipidBlast at all.** The
only one that exists is optimal polarity — all 5,608 rank-1 positive `PC` rows are
`OptimalPolarity=false` with purity 0.

**And `Spectral Components` is not discarded** — it ships per injection in `Results/<pol>/search/*.csv`,
alongside `Optimal Polarity` and `LipiDex Spectrum`. `Spectral_Components.csv` aggregates it to the
compound group, which the per-injection files cannot do; the claim that the data was thrown away was
wrong.

## The `sn Evidence` column

Chain resolution answers *which* chains; it says nothing about *where* they sit. The spectra that
resolved the chains usually also carry that information, and until now nothing read it.

In negative-mode dissociation of a diacyl glycerophospholipid both chains leave as carboxylate
anions, and the **sn-2 carboxylate is typically the more intense** because that ester is the more
labile. `sn_position.py` measures the ratio of the two chains' strongest diagnostic fragments,
using the library's own per-fragment chain annotations, so no fragmentation chemistry is
re-derived.

**It reports; it never asserts.** The name keeps `_`. Writing `/` would claim a proven sn-position
under LSI, and an intensity ratio is not proof — the ratio moves with chain length, unsaturation,
collision energy and instrument, and it inverts for some species. The column sits beside the name
where a reader can weigh it.

**It ships in ONE place — the per-injection search tables — and deliberately not in the delivered
one:**

| table | contents |
|---|---|
| `Results/<pol>/search/*_Results.csv` | one cell per matched spectrum, **never blank** — a decline carries its reason |
| `Final_Results.csv` | **nothing. Withdrawn** — see the validation section below |

It was wired end to end, aggregated onto the delivered row, measured, and then pulled. A delivered
value would be read as a position assignment, and it is not one. The per-spectrum cell stays because
it is a diagnostic — it is what made the failure visible — and it sits beside the `Dot Product` a
reader needs in order to weigh it. The aggregation code is in commit `5fe6dea` if the premise is
ever repaired.

### Two traps, both of which bit

⚠ **The pipeline does not use `search.py::iter_results`.** `pipeline.py` builds its own row dict
inline; `iter_results` and `search_file` serve the tests and the validation harness. The first
attempt at this column wired `iter_results` alone, so on real data **the code never ran** and the
column came back empty on all 430 matches — including declines, which should have carried a reason.
That read as "nothing qualified" and was diagnosed as the wrong hook point, when the hook was
simply in a function the pipeline never calls. Both builders now call one helper,
`search.py::sn_evidence_cell`, and a test asserts the production path fills the column.

⚠ **The rule is negative-mode chemistry and originally crossed polarities.** A protonated
phospholipid does not give carboxylates — it fragments to the head group, and its acyl-related ions
are neutral losses whose intensity ordering reflects which loss is favoured, not ester lability. The
module guarded TG for precisely this reason and then let a whole polarity through: the first wired
run made three calls on **`PE 18:0_20:4 [M+H]+`**. Both the sample's polarity and the library
entry's adduct must be negative, and they can disagree because matching is on precursor mass alone.

⚠ **The column is appended, never inserted.** `Potential Fragments` de-duplicates masses with a
substring test against the row built so far, and in `pipeline.py` that prefix is a join over the row
dict's values. A column added ahead of it changes which masses are written and so changes the
in-source-fragment filter downstream — an evidence column silently editing the delivered table.

### What it cannot do

It cannot be recomputed later. It needs the sample intensities and the library's per-fragment chain
annotations, and neither survives into a `_Results.csv` row — `peaks.py::Lipid` is parsed back from
that CSV. So the search stage is the only place it can be computed, and the written format and
`sn_position.parse_cell` are a matched pair that must change together.

It also inherits the quality of the identification it describes. A call on a dot product of 46
describes an arrangement of chains the row probably does not have. **No dot-product gate is applied
inside this column** — inventing one is the failure mode documented above for `Chain Evidence` — but
the delivered table only ever sees identifications that already cleared `MIN_DOT_PRODUCT` (500) and
`MIN_REV_DOT_PRODUCT` (700) to associate with a compound group.

### ⚠⚠ The premise fails its first validation — do not read the delivered column as evidence yet

Measured on Rat_Heart_Lumos negative mode, three injections, restricted to spectra good enough to
reach the delivered table (`dot > 500`, `reverse dot > 700`): **177 calls, of which 59.3% put the
SATURATED chain at sn-2.**

| | |
|---|---|
| counter-canonical (sn-2 the *less* unsaturated chain) | **105 (59.3%)** |
| canonical (sn-2 the more unsaturated chain) | 72 (40.7%) |

⚠ **Not pseudoreplication.** Counted per distinct species rather than per spectrum — one lipid
contributes many spectra — it is **31 of 54 species (57.4%)**, so repeated counting of a few
species is not what produces it.

### It is carboxylate response, not ester lability, and "the saturated chain wins" is too crude

Collapsed to species and scored as pairwise outcomes, the winner is largely a property of **which
two chains are present**, not of their arrangement:

| chain | outcome across pairs |
|---|---|
| `22:6` | **7 W – 1 L** |
| `22:5` | 4 W – 1 L |
| `18:1` | 8 W – 4 L |
| `20:4` | **7 W – 7 L** |
| `18:0` | **7 W – 7 L** |
| `18:2` | **2 W – 14 L** |
| `22:4` | 0 W – 3 L |
| `20:3` | 0 W – 2 L |

A single fixed ranking of fatty acids predicts **79.6%** of the 54 outcomes. That number is fitted
on the same pairs and so is an upper bound, but the structure is not in doubt.

★★ **The decisive case is `18:0` vs `20:4`: 7–7, an exact coin flip.** `PE 18:0/20:4` is the
textbook arachidonate phospholipid — sn-1 stearate, sn-2 arachidonate — and it is where a working
method should be most confident. It has no opinion. Meanwhile `22:6` and `22:5` win nearly
everything, which reads as canonical but is not evidence of anything: they would win at sn-1 too.

So the ratio tracks **how well each fatty acid gives a carboxylate anion**, and where two chains
have similar response the result is random. Ester lability, if present, is buried under it.

⚠ **An earlier read of this same data said the opposite, and the difference was the filter.** Over
*all* matched spectra the top calls are 20:4 and 22:6 at sn-2, which looks like textbook confirmation.
Restricting to delivery-grade matches inverts the balance. A validation that selects on nothing is
not a validation.

**Outcome: the delivered column was withdrawn** (`Final_Results.csv` no longer carries
`sn Evidence`, and it is unregistered from `META_COLUMNS`); a test fails if it reappears. **What
stands** is the wiring, which is correct and tested, and the per-spectrum column in `search/`, which
is the diagnostic that made the problem visible in the first place.

This follows the precedent set by the `Cer[EOS]` libraries: measured against its own test, failed
it, withdrawn the same day, evidence kept on disk.

**What would fix it** is a comparison that cancels the chain term rather than a threshold on the
ratio — the same molecule's two carboxylates compared against the *same pair* seen in a lipid of
known arrangement, or a per-fatty-acid response normalisation derived from species where the
arrangement is not in doubt. Until then this column measures a ratio, and the ratio is not position.
