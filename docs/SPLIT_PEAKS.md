# Split peaks

One molecule whose chromatographic peak was integrated as two compound groups, so its signal is
divided between two rows. Distinct from an isomer pair, which is two molecules, and from an
adduct pair, which is one molecule seen as two ions. All three look the same in a results table —
several rows at one retention time — and they need three different treatments.

    split peak    one molecule, one ion, one peak cut in two    ->  SUM the rows
    adduct pair   one molecule, two ions                        ->  keep one, do not sum
    isomers       two molecules                                 ->  leave both alone

## Telling them apart

`Adduct` is what makes this decidable. Same-identification pairs within two peak widths, Kiterie
positive run, 53 injections:

| | pairs | median rho | median QC CV, worse half -> summed | median abs delta m/z |
|---|---|---|---|---|
| **same adduct** | 37 | +0.14 | **65.4% -> 13.6%** | **0.0006** |
| different adduct | 3 | +0.48 | 148.8% -> 26.0% | 4.9549 |

Three orders of magnitude between 0.6 mDa and 4.95 Da. The different-adduct pairs are real
chemistry and double as a check on the adduct assignment: `Cer[NS] d41:1` at delta m/z **18.0105**
is `[M+H]+` against `[M+H-H2O]+`, water to four decimal places, and `TG 52:5` at 4.9549 is the
NH4/Na spacing, 4.9554.

## Correlation is the wrong test, and backwards

The obvious rule — merge rows that are close and correlate highly — **merges nothing**. Of 36
same-name close pairs, none reached rho 0.8; the median was **+0.16** and nine were negative.

A split peak divides a conserved total, and the integration boundary moves between injections, so
when one half gets more the other gets less. **Splits are uncorrelated or anti-correlated.** High
correlation instead marks two genuinely distinct co-regulated species — precisely the pairs that
must be left alone. A correlation rule keeps the splits and merges the isomers.

## The test that is used

**Does summing improve precision in the pooled QCs?** Every QC vial holds the same material, so
variance between QC injections is technical by construction. If summing two rows makes that
variance fall, the division between them was technical. It is not a proxy for the question — it
is the question.

    PC 40:3    rho -0.62    QC CV 107% / 60%   ->  20%
    PE 40:4    rho -0.23    QC CV 185% /  9%   ->   7%

## The rule

Two rows are merged when **all** of these hold:

| condition | default | why |
|---|---|---|
| same identification | — | |
| **same adduct** | — | two ions of one molecule have two response factors; summing them adds two scales |
| within `max_peak_widths` of each other | 2.0 | in units of the run's mean FWHM, so it travels between gradients — a fixed 0.2 min window is two peak widths on one method and half of one on another |
| quant ions within `max_ppm` | 10 | nearly redundant given the two above, which is the point: it costs nothing and blocks a naming coincidence |
| summed QC CV beats the worse half by `min_cv_improvement` | 1.5 | splits improve ~4.8x, distinct co-eluting rows ~1.1x, so 1.5 sits in open ground rather than on a boundary needing tuning |

The larger row survives and takes the summed areas; the loser is marked not-kept with
`Split peak, merged into <name> at <rt> min`, so `Unfiltered_Results.csv` still shows it and the
decision is auditable. Compounds move across too, so `Features Found` counts the union rather than
the winning half's injections. An absorbed row is out of later comparisons, so one peak in three
fragments collapses to one row rather than a merged pair plus an orphan.

## With no QCs, nothing is merged

Below `min_pools` (default 3) pooled QC injections the stage does nothing and says so in the log.
The whole criterion is QC precision; without QCs there is no evidence, and guessing would change
quantification on exactly the runs least able to show it had gone wrong. The peak finder's own
retention-only redundancy filter still applies.

## Where it runs

After the peak finder, before the blank filter — so a group is judged against the blanks on the
whole of its signal rather than on whichever half happened to be the larger row. It is deliberately
outside the peak finder: that reproduces LipiDex, and this is an addition on top of it, like the
blank filter.

The existing `merge_compound_group` inside the finder is a different thing. It pools *evidence*
for redundant identifications and sets `keep = False` on the loser, but `areas` is built at group
construction, so it never sums. Under it, the smaller half of a split peak was discarded.

## The adduct case is handled separately, and oppositely

`adduct_pairs.py` removes the weaker of two identified rows that are two adducts of one molecule —
see the mass coincidence above. It is worth holding the two side by side, because the same
statistic points in opposite directions:

| | split peak | adduct pair |
|---|---|---|
| what it is | one molecule, one ion, one peak cut in two | one molecule, two ions, two names |
| the rows | same identification | different identifications |
| correlation | **anti-correlated** — the integration boundary moves, so the halves trade signal | **must track** — two ions of the same eluting compound |
| treatment | **sum** | **remove the ghost** |

Summing an adduct pair would add two ions with two response factors, and add a real measurement to
an artefact. Removing a split peak would throw away half the signal. Getting the two confused is
the failure mode both modules are written to prevent.

## Configuration

    "split_peaks": {
      "enabled": true,
      "max_peak_widths": 2.0,
      "max_ppm": 10.0,
      "min_cv_improvement": 1.5,
      "min_pools": 3
    }

Written into `run_config.json` with everything else, so a merged table carries the rule that
merged it.

## After merging: judging what is left

Rows still sharing a name after merging are isomers, and `Duplicate Verdict` says what, if
anything, argues against one:

| verdict | scope | meaning |
|---|---|---|
| `unreliable` | the row | pooled-QC CV above 30% (Dunn et al. 2011) |
| `class model?` | **the whole name** | every row of it is far from the class retention prediction |
| `isomer?` | the row | nothing argues against it |

⚠ **Retention does not judge individual rows here, and a verdict that did was withdrawn.** The
retention model issues one prediction per sum composition, so all but one isomer is displaced from
it by construction — the verdict meant "you are not the isomer nearest the class average". See
`docs/RTLS.md` for the measurement behind that and for `class model?`, which judges the set rather
than a row.

This leaves precision as the only per-row criterion. The loss is real rather than a tidying-up: the
two agreed only weakly (ρ +0.22), so retention was a genuine second view — it was the wrong view
for this question.

Nothing here deletes. Of 165 duplicated names on the study this was built against, a minority are
marked; the rest are several rows all precise, which is what genuine chromatographic isomers
sharing a sum-composition label look like. Choosing between those needs a chromatographic argument,
not a threshold — and `docs/CHAIN_RESOLUTION.md` covers why the chain names themselves are weaker
evidence than they look.
