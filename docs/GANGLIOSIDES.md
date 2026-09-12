# A negative-mode ganglioside library

Neither LipiDex version ships one. `LipiDex2_Ganglioside.msp` is 4,305 entries and **every one
is `[M+H]+`**; LipidBlast and the LipiDex HCD libraries contain no ganglioside entry in any
polarity. So negative-mode data — where sialylated species ionise best — could not be searched
for them at all.

`scripts/build_ganglioside_negative_library.py` generates one from the species LipiDex 2 already
enumerates, and the interesting part is what resolution it is written at.

## What the fragmentation actually supports

The richest GM3 negative spectrum in this dataset contains three ions worth the name:

| m/z | intensity | assignment |
|---|---|---|
| 290.0858 | 398,762 | [NeuAc − H₂O − H]⁻ |
| 87.0074 | 52,392 | C₃H₃O₃⁻, sialic acid |
| 272.0759 | 9,145 | [NeuAc − 2H₂O − H]⁻ |

and these, checked for explicitly, are **absent**:

| expected | |
|---|---|
| [M−H−NeuAc]⁻ glycosidic Y-ion | absent |
| [M−H−NeuAc−Hex]⁻ | absent |
| ceramide fragment | absent |
| [24:0 FA−H]⁻ acyl anion | absent |

Every one of those would be needed to say which chains the molecule carries. The spectrum says
*this is sialylated* and stops.

## So the library is written at sum composition

A library at molecular resolution would give `GM3 d18:1_24:0` and `GM3 d20:1_22:0` **identical**
spectra, and the dot product could not separate them — the single-fragment degeneracy documented
for the 184.07 phosphocholine entries in [VALIDATION.md](VALIDATION.md), in a more extreme form.

Written at sum composition it says exactly what the evidence says: the precursor mass gives total
carbons and double bonds, the fragments give the head group, the chain split is not determined.
4,305 molecular species collapse to **931 sum compositions**.

The entries carry no moiety-fragment annotations, so the purity calculation returns 0 and the
peak finder reports them at sum composition of its own accord. Library and reporting agree
because both reflect the same absence of evidence.

The three intensities are the mean relative to base peak across the five species independently
confirmed in positive mode — a template calibrated on confirmed examples, not fitted to the
species being searched for.

## Result

Negative mode, 26 files. Nothing displaced: identifications go from 4,494 to 4,573, the extra 79
all ganglioside, and positive mode is bit-identical with and without it (3,130 either way, none
matching the new library).

**Eight GM3 species, where negative mode previously had none:**

| identification | RT | quant ion | files | MS2 spectra |
|---|---|---|---|---|
| GM3-NANA d34:1 | 8.26 | 1151.6988 | 18 | 5 |
| GM3-NANA d36:1 | 8.61 | 1179.7296 | 14 | 10 |
| GM3-NANA d38:2 | 8.67 | 1205.7445 | 9 | 2 |
| GM3-NANA d38:1 | 9.00 | 1207.7614 | 18 | 13 |
| GM3-NANA d42:2 | 9.44 | 1261.8068 | 12 | 4 |
| GM3-NANA d40:1 | 9.45 | 1235.7928 | 18 | 16 |
| GM3-NANA d41:1 | 9.70 | 1249.8084 | 17 | 12 |
| GM3-NANA d42:1 | 9.97 | 1263.8244 | 20 | 16 |

Three of these — **d36:1, d38:2 and d42:2** — are the previously unidentified precursors found by
sweeping for the sialic acid fragment. They were in the data all along with nothing to identify
them against.

**The retention behaviour is independent confirmation.** Retention rises monotonically with
carbon number across d34:1 → d42:1, and each `:2` species elutes before its `:1` counterpart —
d38:2 at 8.67 against d38:1 at 9.00, d42:2 at 9.44 against d42:1 at 9.97. One more double bond
elutes earlier, which is what the retention model in [RTLS.md](RTLS.md) formalises and which no
part of this library encodes.

## Two caveats

**The weakest two are weak.** `d38:2` rests on 2 spectra and `d42:2` on 4. The `MS2 Spectra`
column carries that.

**The positive-mode molecular names are supported — I was too pessimistic above.** Checked
rather than assumed: positive mode does carry chain-diagnostic fragments, though not the ones the
library leans on hardest.

| | n spectra | sphingoid 264.2686 | ceramide Y-ion | **both** | direct N-acyl ion |
|---|---|---|---|---|---|
| d18:1_16:0 | 13 | 12 | 12 | **12** | **0** |
| d18:1_22:0 | 15 | 7 | 7 | **7** | 3 |
| d18:1_24:0 | 15 | 15 | 15 | **15** | 7 |

The **direct N-acyl fragment is weak or absent** — 0 of 13 spectra for 16:0, and only 3 of 15 for
22:0, at around 0.5% of base peak. On its own that would leave the acyl chain undetermined, which
is what prompted the caveat.

But it is not needed. The **sphingoid ion fixes the long-chain base** (264.2686 = d18:1, strong,
at 458–953 relative intensity) and the **ceramide Y-ion fixes the total ceramide** (520.51 /
604.60 / 632.63, at 353–680). The N-acyl chain then follows by subtraction: total minus
sphingoid. For `d18:1_16:0`, a Y-ion at 520.51 means ceramide d34:1 and a sphingoid ion at 264.27
means d18:1, so the acyl is 16:0 and nothing else. A `d20:1_14:0` would give the same Y-ion but a
sphingoid ion at 292.30, which is not observed.

Both determining ions co-occur in the **same spectrum** 12/13, 7/15 and 15/15 times. So the
molecular assignments stand, on a different pair of ions than the library's own weighting
suggests.

The contrast with negative mode is the real point: positive gives the sphingoid base and the
ceramide, negative gives only the head group. **For structure these species are better served in
positive mode, not negative** — the opposite of the usual advice for gangliosides, which is about
sensitivity rather than structural information.
