# Adduct pairs

Two identified rows at one retention time whose masses say they are two adducts of **one**
molecule. Their neutral masses agree, so at most one of the two names can be true.

## The hole this closes

The peak finder's adduct sweep, faithfully reproduced from LipiDex, only fires when one member of
a pair is unidentified:

| branch | fires when |
|---|---|
| `check_adduct_known` | one identified, the other **not** |
| `check_adduct_unknown` | **both** unidentified |

There is no branch for two identified groups, so a pair that both carry names is never compared.

## Why it happens often enough to matter

    Na - H                        = 21.98194
    +2 carbons, +3 double bonds   = 21.98435   (C2H4 - 3 x H2)
    difference   2.40 mDa   =   3.0 ppm at m/z 800

**The sodium adduct of `PE 40:1` sits 3 ppm from protonated `PE 42:4`** — inside any ordinary
search window, at the same retention time, in a class that supplies both rungs of the ladder.

Negative mode barely suffers: `[M-H]-`, `[M+Ac-H]-` and `[M+Cl]-` spacings do not collide with the
CH2/double-bond ladder the way `Na - H` does. On the brain run the mass screen flags 412 positive
pairs against 20 negative.

## Mass agreement is a screen, not a verdict

Of those 412 positive pairs, **only 14 also track across the samples**. With eight positive adducts
the mass test alone is permissive and co-elution inside a lipid class is ordinary, so acting on
mass alone would delete real lipids by the hundred.

Correlation is the confirmation: two adducts of one molecule are two ions of the same eluting
compound and must rise and fall together. The adduct ratio drifts with matrix sodium, which adds
noise, but cannot make them disagree.

⚠ **This is the reverse of the split-peak rule, and the reversal is not a special case.** There,
correlation is actively misleading — a split peak divides a conserved total and its halves are
anti-correlated. Here it is required. The statistic is the same; the physical claim being tested is
not. See [SPLIT_PEAKS.md](SPLIT_PEAKS.md) for the two side by side.

## Which row survives

1. **Canonical adduct for the class** — `[M+NH4]+` for TG/DG/CE, `[M+H]+` for the phospholipids
   and sphingolipids.
2. **Then dot product.**

Canonical first rather than score first, deliberately. A sodiated species barely fragments, so a
sodiated identification rests on fewer fragments; a dot product earned on three peaks is not
comparable to one earned on twelve, and ranking on the number alone would let the ghost win.

The loser is marked `Adduct pair: the [M+Na]+ of <winner> at <rt> min`, so it stays visible in
`Unfiltered_Results.csv` and the decision can be audited.

## It removes, it never sums

The opposite of a split peak. Adding the two would combine a real measurement with an artefact,
and combine two ions with two response factors. Nothing is added to the winner.

## What it does on real data

Kiterie brain, 53 injections per polarity:

| | mass-consistent pairs | removed | rejected: not tracking | ambiguous | isotope spacings |
|---|---|---|---|---|---|
| positive | 108 | **18** | 89 | 1 | 7 |
| negative | 2 | **0** | 2 | 0 | 6 |

Every removal, audited individually on sum composition:

| chain difference | n | what it is |
|---|---|---|
| +2C, +3DB | 14 | the `Na - H` coincidence — `LysoPE 20:4` as the `[M+Na]+` of `LysoPE 18:1`, `PE 42:4` of `PE 40:1` |
| -2C, -3DB | 1 | the same, in the other direction |
| +5C, +3DB | 1 | the same, across classes: PC and PE headgroups differ by 3 x CH2 = 42.047, so `PE 46:4` = `PE 44:1` = `PC 41:1` + (Na - H) |
| +0C, +2DB | 1 | **water loss** — `Plasmanyl-PC O-42:4` is the `[M+H-H2O]+` of `PC 42:2`, at -18.0030 |
| 0, 0 | 1 | one molecule on two ions |

**16 of 18 are the sodium coincidence**, and every one of the remaining two is a known mass
relationship. Negative removing nothing is the expected control.

⚠ **The water-loss case is worth its own note.** In-source dehydration of a diacyl species lands on
the mass of an ether lipid with one more double bond, so it inflates plasmalogen and plasmanyl
counts specifically — the class this pipeline already warns cannot be assigned reliably in positive
mode (`docs/VALIDATION.md`). This filter catches it only when the two co-elute and track; it is not
a general solution to that problem.

## Two cases, and they are not the same claim

| the two rows are | what it means | reason recorded |
|---|---|---|
| **different molecules** | their neutral masses agree, so at most one NAME is right | `the [M+Na]+ of <winner>` |
| **the same molecule** | both names are right; one ROW is redundant | `the same molecule on its [M+Na]+ of <winner>` |

The second case is one molecule seen as two ions — e.g. `TG 18:2_18:1_20:1` against `TG 56:4`,
separated by the NH4-to-Na spacing of 4.9554. Nothing was misidentified there, and saying so would
be wrong, but the redundant row still has to go: the same molecule must not be quantified twice on
two ions with two response factors.

⚠ **Sameness is decided on sum composition, never on the displayed name.** `TG 18:2_18:1_20:1` and
`TG 56:4` are one molecule at two resolutions, and comparing the strings called them different —
which is exactly what `sum_composition` exists to prevent, and what the peak finder's own
redundancy filter already uses. A pair that is the same molecule on the **same** ion is left to
`split_peaks.py`, which sums it; this module only takes the same molecule on **different** ions,
which cannot be summed.

## Two defects the first run had, found by auditing every removal

Worth recording, because both produced plausible-looking output.

**An isotope spacing read as an adduct.** `SM d43:2` was removed against `HexCer[NS] d43:1` on a
gap of 1.0039 — a 13C spacing, not an adduct. `[M]+` and `[M+H]+` differ by 1.00728, which is 3.9
mDa away and inside 10 ppm at m/z 800, so any M+1 isotope surviving the isotope filter matched. Now
guarded explicitly; it skips 6-7 pairs per polarity.

**The direction came from the ranking, not the masses.** `explained` means "second is an adduct of
first" and licenses removing *second* — nothing else. Ranking the pair independently and deleting
whichever scored worse could delete the other row, on a relationship that did not describe it. The
tell was the reason string reading `the None of ...`. The direction now comes from the mass
evidence, ranking only breaks a symmetric case, and a pair whose accused row is the better
supported one is left alone as **ambiguous** rather than acted on.

## Configuration

    "adduct_pairs": {
      "enabled": true,
      "max_peak_widths": 1.0,
      "max_ppm": 10.0,
      "min_correlation": 0.8,
      "min_points": 8,
      "correlation_type": "spearman"
    }
