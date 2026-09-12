# Abundance correlation as a redundancy test

LipiDex 2 adds a correlation filter; LipiDex 1 has nothing like it. Implemented in
`correlate.py` in two forms, only one of which survived measurement.

The idea is sound. An adduct, an in-source fragment and an isotope of one molecule are not
independent measurements of anything — they are the same peak, split. Whatever the molecule does
across the sample set, all of them do together. That is strictly more information than the mass
tests in `adducts.py`, which ask only whether two quant ions *could* be the same neutral
molecule, a question with a fair number of coincidental yes answers at 20 ppm.

## Form 1: the standalone filter — measured, and left off

LipiDex 2's description: *"Filter redundant unidentified features that have high correlation and
co-elute with identified features."* No mass relationship required — co-elution plus correlation
is the whole test.

Only unidentified features are removed, so no identification can be lost to it. But that safety
property says nothing about whether the rows removed *should* be, and here they should not.

Applying the same rule to the **identified** rows gives a false-positive estimate: those have
their own MS2, so a hit there is the filter calling a real compound redundant.

| threshold | unidentified removed | identified that would be hit |
|---|---|---|
| 0.80 | 1,090 / 3,962 (27.5%) | 666 / 908 (**73.3%**) |
| 0.90 | 668 (16.9%) | 447 (**49.2%**) |
| 0.95 | 394 (9.9%) | 252 (**27.8%**) |
| 0.98 | 199 (5.0%) | 108 (11.9%) |
| 0.99 | 112 (2.8%) | 48 (5.3%) |

At any threshold that removes a useful number of rows, the same rule would condemn a quarter or
more of the compounds we have spectra for. Co-eluting lipids genuinely track each other — they
share injection-to-injection loading, and lipids of one class at one retention time are often
co-regulated — so "correlates with something identified nearby" is close to universal.

The estimate is not perfect: identified rows are more abundant and better measured, and two
identified rows at one retention time may genuinely be an adduct pair. It overstates the rate
somewhat. Not by enough to matter at 27.8%.

**Default: off.**

## Form 2: as a guard on a specific mass relationship — kept

The same measurement points at where correlation *is* informative. Taking each removal the mass
sweep already makes and asking how well it correlates with a co-eluting identified peak:

| removed as | n | median r | ≥0.95 |
|---|---|---|---|
| In-source fragment | 341 | **0.96** | 54.3% |
| Adduct of existing identified peak | 469 | **0.93** | 39.9% |
| *(kept, unidentified — baseline)* | 1,844 | 0.85 | 21.4% |

Both categories where the partner is by construction an *identified* peak come out clearly
enriched against the baseline. The mass tests are mostly right, and correlation agrees with them.

So correlation is used as a **confirmation of a specific pair** rather than a filter in its own
right: before removing a group as an adduct, in-source fragment, dimer or isotope of a particular
partner, the two must track each other. Where there is too little overlap to judge, the mass
evidence stands on its own — silence does not overrule it.

**Default: on, at r ≥ 0.5** — a low bar, deliberately. It is not there to be selective; it is
there to refuse the minority of removals where the two ions demonstrably do not behave like one
molecule.

Cost, measured: it blocks about 214 removals, 3.6% of them, taking output rows from 4,870 to
5,084 against the reference's 4,128. **Reference recovery, identified rows and distinct
molecules are all unchanged** — it touches unidentified rows only.

That cost is accepted on an asymmetry rather than on a measured gain, and the reasoning should be
visible: wrongly removing a peak loses a real compound, wrongly keeping one leaves an extra
unannotated row. The two are not equally bad. There is no ground truth for unidentified rows, so
the benefit cannot be measured directly — if row count matters more to you than that asymmetry,
set `guard_mass_relationships` to false.

## Note for a different experiment

Both forms depend on how much genuine variation the sample set carries. These are 20 biological
samples plus QCs; on a set of replicates of one material every correlation goes to one and the
standalone filter would remove nearly everything, while the guard would block nothing. Neither
default should be assumed to transfer without re-running the two tables above.
