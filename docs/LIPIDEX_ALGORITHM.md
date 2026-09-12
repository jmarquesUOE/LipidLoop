# LipiDex identification logic, as implemented

Extracted from the LipiDex Java source (`C:\LipiDex\LipiDex\src`, Coon lab) so this pipeline
reproduces its behaviour rather than approximating it. Every number below is quoted from the
source, with the file and symbol it came from. Nothing here is inferred.

Reference: Hutchins, Russell & Coon (2018) *Cell Systems* 6:621. LipiDex is GUI-only —
`App.java` is a nine-line Swing launcher with no command-line entry point, which is why this
pipeline reimplements the logic instead of driving the binary.

---

## 1. Spectral similarity

`spectrum_searcher/SampleSpectrum.java :: calcDotProduct(libArray, mzTol, reverse, massWeight, intWeight)`

Sample and library peak lists are sorted by m/z and merged with a tolerance of `mzTol`. Each
peak ends up in one of three states: matched in both, library-only, or sample-only.

For every peak, with `m` = m/z, `I` = intensity:

    weighted(m, I) = m^massWeight * I^intWeight

Sums are accumulated as:

| peak state | contributes to |
|---|---|
| in both | `numerSum`, `libSum`, `sampleSum` |
| library only | `libSum` (forward only) |
| sample only | `sampleSum` (forward only), **at half intensity** |

    score = 1000 * numerSum^2 / (sampleSum * libSum)

`reverse = true` skips the library-only and sample-only terms entirely, so the reverse dot
product ignores sample peaks the library does not predict. That is what makes it robust to
co-isolated contaminant fragments.

**Weights** (`SpectrumSearcher.java`, lines 43-44):

    intWeight  = 1.2
    massWeight = 0.9

**Note the asymmetry**: for matched peaks the code uses `m^massWeight`, but for unmatched peaks
it uses `massWeight * m`. This is in the original and is reproduced here for fidelity.

**Mass binning** (`SpectrumSearcher.java`, line 67): `massBinSize = ms1Tol / 10.0`

Unmatched sample peaks are halved (`sampleIntensities / 2.0`) before weighting, in the forward
score only.

---

## 2. Identification purity

`spectrum_searcher/PeakPurity.java :: calcPurity(...)`

**This is not spectral noise.** It is the confidence that the top-scoring lipid explains the
spectrum, given every isobaric candidate that could also explain it.

1. For the top match, compute an intensity contribution from its diagnostic fatty-acid
   fragments (`getPurityIntensity`, using the fatty-acid database).
2. Do the same for every isobaric candidate identification that is a unique lipid.
3. Purity = round(top_intensity / sum_of_all_intensities * 100)

`minPurity = 5.0` — candidates above this contribute their mass to in-source fragment
screening. Note it is compared against a raw *intensity*, not against the percentage, despite
the name.

Reimplementing this requires the fatty-acid database, because the contributions are computed
from FA fragment intensities rather than from the full spectrum. It is
`src\backup\FattyAcids.csv` that the searcher reads, not the per-library copies under
`src\libraries\`, and the `Enabled` column is ignored — disabled acids still take part.

Implemented in `lipidloop/purity.py`. The details that decide the number:

* **Only the top-ranked identification is scored**, and only if it is a LipiDex spectrum in its
  optimal polarity. Off-polarity entries return 0 — which is why positive-mode PC, carrying
  `OptimalPolarity=false`, reports 0 throughout a positive run.
  ⚠ The `isLipidex` half of that gate **never fires here**: `lib_gen/Lipid.java` stamps
  `Type=LipiDex` unconditionally on every generated library, LipidBlast-derived ones included, and
  `SpectrumSearcher.java:628` sets the flag on an unscoped `line.contains("LipiDex")`. Off-polarity
  is the only by-design zero.
* **The moiety type is matched exactly, formula prefix and all.** A transition's type is the
  annotation up to its last underscore, so `-_Alkyl Fragment_[16:0]` gives `-_Alkyl Fragment`.
  Moiety fragments carry `-` where a formula would go, which is what lets the chains of one
  lipid group together; the generic fragments carry a real formula and are excluded by name
  (`_Fragment`, `_Neutral Loss`, `DG Fragment`, `PUFA`).
* **A missing chain is fatal.** `faCount != faDB.size()` returns null and the candidate
  contributes nothing at all, so a lipid whose second chain fragment fell below the 5%
  intensity floor is absent from the denominator rather than penalised within it.
* **Repeated chains are claimed once each.** `numOccurences` counts how often the chain appears
  in the lipid's own name, so `DG 16:0_16:0` claims its single 16:0 fragment twice.
* **Shared fragments are corrected across candidates**, keyed on library mass to within 1e-6.
  Isobaric candidates whose fragment masses differ in the fourth decimal — as the same lipid
  does between LipidBlast and the LipiDex libraries — miss each other and both count the peak.
* `getPurityIntensity` returns the median contribution for the top match or for fewer than
  three chains, and the maximum otherwise.
* An unannotated peak in a LipiDex entry throws in the original and is caught into a 0.0.

---

## 3. Filter thresholds

`peak_finder/Utilities.java`, lines 12-17. Hard-coded, with one exception noted below:

| constant | value | meaning |
|---|---|---|
| `MAXPPMDIFF` | 20.0 | maximum precursor mass error, ppm |
| `MINDOTPRODUCT` | 500.0 | forward dot product floor |
| `MINREVDOTPRODUCT` | 700.0 | reverse dot product floor |
| `MINFAPURITY` | 75.0 | minimum fatty-acid purity |
| `MINRTMULTIPLIER` | 0.5 → **2.0** | retention-time window multiplier — see below |
| `MINIDNUM` | 1 | minimum features per compound group to retain |
| `MASSOFELECTRON` | 0.00054858026 | used in adduct mass calculation |

The reverse threshold being higher than the forward one (700 vs 500) is deliberate: a real
identification must explain the library's predicted fragments well, even if the spectrum
carries extra peaks from co-isolation.

**`MINRTMULTIPLIER` is the exception and it is worth knowing about.** The 0.5 in the source is
overwritten at startup from the GUI's *FWHM* spinner, whose default is **2.0** — so 2.0 is the
value that actually runs. It then governs two unrelated things: the FWHM multiple inside which an
MS2 may be attributed to a feature, and the half-width of the per-class retention window in
`checkClassRTDist`.

There is a *second* field of the same name on `CDPeakFinder`, set from its own spinner (default
3.5) and passed to `checkClassRTDist` as an argument **that method never reads**. It is dead in
the original, so it is not reproduced here. The practical consequence: the retention window is
2.965 × MAD wide rather than the 5.19 × MAD the visible spinner implies, and no LipiDex user has
ever been able to widen it from that spinner. See [FILTERS.md](FILTERS.md) for what that costs.

---

## 4. Output columns

`SpectrumSearcher.java` line 409 writes:

    Precursor Mass, Library Mass, Delta m/z, Dot Product, Reverse Dot Product, Purity,
    Spectral Components, Optimal Polarity, LipiDex Spectrum, Library, Potential Fragments

This pipeline emits the same columns so existing downstream analysis reads its output
unchanged.

---

## 5. Libraries

`.msp` format, reused as-is. Counts of the non-canonical ceramide subclasses, which is what
motivated this project:

| library | Cer[ADS] | Cer[NDS] | Cer[AP] |
|---|---|---|---|
| LipidBlast_Formic / _Acetate | 2,408 | 7,224 | 1,290 |
| LipiDex_HCD_Formic | 0 | 0 | 640 |
| LipiDex_HCD_Acetate | 0 | 0 | 960 |

**Cer[ADS] exists only in LipidBlast, and only in negative mode** (`[M-H]-` and `[M+FA-H]-`,
602 entries each). A search that omits LipidBlast cannot return Cer[ADS] at all, regardless of
signal. Library selection must therefore be explicit and recorded in the run config.
