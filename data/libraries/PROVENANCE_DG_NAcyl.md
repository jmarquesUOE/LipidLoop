# Diacylglycerols and N-acyl amides — provenance

    python scripts/build_diacylglycerol_library.py --out data/libraries/Diacylglycerols_Positive.msp \
                                                   --targets <file of `DG C:D` lines>
    python scripts/build_nacyl_amide_library.py --out-pos data/libraries/NAcylAmides_Positive.msp \
                                                --out-neg data/libraries/NAcylAmides_Negative.msp
    # then a matching decoy for EACH, via scripts/make_decoy_library.py

| library | entries | decoys | verified against |
|---|---|---|---|
| `Diacylglycerols_Positive` | 342 | 342 | LipiDex's own 1,025 DG entries, ≤1 mDa |
| `NAcylAmides_Positive` (NAE, NAGly, GPNAE) | 333 | 333 | 4 named compounds |
| `NAcylAmides_Negative` (NATau) | 111 | 111 | 1 named compound |

## Diacylglycerols — a double-bond gap, not a chain gap

The shipped DG library holds 1,766 entries built from 40 acyl chains, but its unsaturation coverage
is patchy: C19/21/23/25 appear **only** as saturated, C24 only as 0/1/4, C22 has no 22:3, C16 stops
at 16:1. A diacylglycerol carrying a polyunsaturated odd chain has no entry, so 78 sum compositions
published across five deposits could not be matched — each reading as a detection failure rather
than a library gap.

**Enumerated to the observed sums only.** An exhaustive C10–C26 enumeration would be ~5,000 entries
that nothing has observed, each one a chance for a false positive at a plausible mass.

⚠ **Chains stop at C26, matching the shipped library's range.** 41 of the 78 sums are unreachable
that way and the builder prints every one. They are not a coverage gap: `DG 53:6` through
`DG 61:11` would need C27–C30 acyls, which do not occur in a diacylglycerol, and `DG 72:0` would
need two saturated C36 chains. Their carbon counts sit in **triacylglycerol** range, which is what
an in-source fatty-acid loss from a TG produces — see the in-source fragment section of paper 1 §8.

## N-acyl amides — weaker evidence, deliberately

⚠ **Read a hit from this library differently to one from the DG or acylcarnitine libraries.** Those
were generated from a rule inferred from LipiDex's own entries and then checked by regenerating
them: the library validated itself. **No N-acyl amide exists in any shipped library here**, so
there was nothing to check the fragmentation against.

- **Precursor masses are solid** — formula arithmetic, cross-checked against anandamide,
  oleoylethanolamide, palmitoylethanolamide, N-oleoylglycine and N-oleoyltaurine.
- **Fragmentation is asserted from published behaviour**, not measured here and not validated
  against an in-house spectrum.

So a hit rests on accurate mass plus a head-group ion plus one chain-specific ion. Authentic
standards would settle it; until then this is chemistry, not measurement.

⚠ **The negative arm needed a chain-specific fragment to be searchable at all.** The first version
gave NATau only the taurine ion at *m/z* 124.0074. With no chain-specific peak, nothing
distinguishes `NATau 18:1` from `NATau 20:1` beyond the precursor — and `make_decoy_library`
correctly wrote **zero** decoys, because a decoy of such an entry is identical to its target. That
would have put 111 targets into the search with no decoy competition, biasing the measured FDR low
on precisely the class with the weakest evidence. The acyl carboxylate anion is now included.
