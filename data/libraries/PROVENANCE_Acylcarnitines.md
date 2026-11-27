# Acylcarnitines, positive mode — provenance

Built by `scripts/build_acylcarnitine_library.py`, which is the reproducible record.

    python scripts/build_acylcarnitine_library.py --out data/libraries/Acylcarnitines_Positive.msp
    python scripts/make_decoy_library.py data/libraries/Acylcarnitines_Positive.msp \
                                         data/libraries/DECOY1_Acylcarnitines_Positive.msp

**86 target entries, C0–C26; 81 decoys.** Five targets (C0–C3) have no chain-specific fragment, so
a decoy of them would be identical to the target and is deliberately not written.

## Why it exists

LipiDex and LipidBlast both begin at **C10**. Everything below is absent, which removes precisely
the acylcarnitines that carry metabolic information: acetyl (C2), propionyl (C3), butyryl (C4),
isovaleryl (C5), and free carnitine (C0) as the pool they draw on. Across five reference deposits
the unmatched acylcarnitines run C2–C16.

It is also why acylcarnitines have never appeared in the per-class retention model — that model
fits on confident identifications, and a class with no library produces none.

## How the masses are justified

The fragmentation rule is **inferred from LipiDex's own C10–C26 entries**, so the library is used
as the test: the builder regenerates all 40 shipped species and refuses to write if any precursor
or fragment disagrees by more than **1 mDa**.

⚠ That check alone was not sufficient. The first version wrote every constant fragment 0.55 mDa
high — neutral formula masses for what are cations — and passed a self-check running at 5 mDa.
6.5 ppm at *m/z* 85 would have put the entire class systematically off while looking correct.
`tests/test_acylcarnitine_library.py` therefore checks against **literature masses the rule was not
derived from** (newborn-screening / FAO panel values), which is what actually caught it.

## Known limits

- **Short chains lose the acylium.** Below C4 the acyl fragment falls under the *m/z* floor a real
  acquisition records, and for C0/C2 no distinct acyl ion forms. Those entries carry only the
  constant fragments, so they are identified on precursor plus three shared ions — weaker evidence
  than a long-chain acylcarnitine, and the purity guard should be read accordingly.
- **Isomers are not resolved.** `AC 5:0` is isovaleryl *and* 2-methylbutyryl; `AC 4:0` is butyryl
  *and* isobutyryl. Reversed-phase separates some of these partially and the spectra do not
  distinguish them at all, so the name is the sum composition and nothing more.
- **Unvalidated against authentic standards on this instrument.** The two standards plates include
  carnitines (see the lab notes); until those are injected, these entries are chemistry, not
  measurement.
