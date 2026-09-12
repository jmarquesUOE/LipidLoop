# Free cholesterol — a retracted identification, and why it was wrong

**Conclusion first: there is no detectable free cholesterol in this dataset.** An earlier version
of this document reported two cholesterol peaks. Both were cholesteryl esters fragmenting in the
ion source. The library is no longer in the default set.

## The gap is real

No LipiDex library contains free cholesterol, in either version. `CE` is there — 40 cholesteryl
esters, all `[M+NH4]+` — but the free sterol is absent, so it cannot be identified at all. That
part still stands.

## What went wrong

Cholesterol is seen as `[M+H−H2O]+` at 369.3516. Filtering to that precursor and requiring the
sterol backbone ions gave eleven spectra between 16.1 and 17.1 min, carrying 54–81 peaks and all
fourteen backbone ions, whose fragment masses agree with the series LipiDex already carries
inside its CE entries. Built into a library entry, it identified two compound groups at dot
products of 652–944.

Every check I ran was internally consistent. **None of them tested the alternative explanation:
cholesteryl esters fragment in the source to exactly m/z 369.3516 and exactly those backbone
ions.** A CE in-source fragment is indistinguishable from free cholesterol by mass and by MS2. It
is distinguishable by retention — free cholesterol is far less hydrophobic than a cholesteryl
ester and elutes much earlier — and I never looked.

## The evidence that settles it

The 369.3516 chromatogram over the whole run has exactly two peaks, and each sits on a CE:

| 369.3516 peak | intensity | coincident precursor |
|---|---|---|
| 16.09 | 7.9e6 | **CE 20:4 [M+NH₄]⁺ at 16.09** |
| 17.09 | 1.6e6 | **CE 18:1 [M+NH₄]⁺ at 17.08** |

and the third reported row, at 16.49, sits on **CE 18:2 at 16.48**.

Binned by minute, the trace is flat background (2–7e4) everywhere else. Between 2 and 14 min,
where free cholesterol should elute, the maximum is 3.4e5 — 4% of the CE-region signal and not a
credible peak.

The fragment exceeding its own precursor is the signature: 369.3516 reaches 7.9e6 where CE 20:4
itself reaches 2.5e6, and 1.6e6 where CE 18:1 reaches 7.4e4 — 22-fold. In-source fragmentation of
CE is efficient, so the fragment channel outruns the surviving precursor. Across samples the
17.07 row correlates with CE 20:4 at **r = 0.84**, which is what a fragment does with its parent
and not what two independent lipids do.

## Why the in-source filter did not catch it

The peak finder has an in-source-fragment filter, and it works by asking whether a group's quant
ion appears in the *predicted fragment list of an identified neighbour*. It failed here for a
reason worth recording:

* **Only one CE was identified.** CE 20:4 was; CE 18:1 and CE 18:2 were not, because their
  `[M+NH4]+` precursors are weak — 7.4e4 and 6.9e4 — and never produced a passing MS2. With no
  identified parent, there is nothing for the filter to match against.
* **The one identified CE was too far away.** CE 20:4 sits at 16.09 and the reported row at
  16.49, about 0.4 min apart, well outside the two-peak-width window the filter uses.

So the filter is not broken; it is blind whenever the parent is present but unidentified. That is
common for CE, whose ammonium adducts ionise poorly while their in-source fragment is enormous.

## What would make this safe

An **authentic cholesterol standard** run on the same method, to fix the retention window. That
is the only thing that separates free cholesterol from a CE fragment here, because mass and MS2
cannot. `scripts/build_sterol_library.py` still generates the entry and carries this warning; it
is deliberately not in `DEFAULT_LIBRARIES`.

The same caution applies to any sterol added to it. Sterols sharing a formula are already
indistinguishable from each other by MS2, and sterol esters make it worse by producing the free
sterol's spectrum from a completely different molecule.

## What still holds from the earlier work

* `CE` entries all carry the same spectrum — 369.3516 plus the sterol backbone — because the acyl
  chain contributes nothing. For CE that is sound: the sterol is fixed, so precursor mass gives
  the acyl chain exactly.
* Only **one** CE was identified in 26 files, against a library of 40. The abundant CE species —
  16:0, 18:1, 18:2 — are present in the MS1 data and were not identified. That is a genuine gap,
  and unlike the cholesterol claim it is worth pursuing.
