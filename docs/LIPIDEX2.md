# What LipiDex 2 changed

Read from an installed copy of **LipiDex2 v1.1.2** (`C:\Program Files\LipiDex2-v1.1.2`) — the
config files, the spectral libraries, and the .NET assembly metadata and help strings. Not from
the paper, so treat the algorithmic detail as inferred from names and UI text rather than read
off source. The v1 conclusions elsewhere in these docs came from actual Java source and are
firmer.

## It is a rewrite, not a revision

Java/Swing became **C# / .NET WPF**, with EntityFramework and SQLite. The module set is
`SpectrumSearcher`, `PeakFinder`, `LibraryGenerator`, `LibraryForge`, `SpectralAnnotator`,
`LipidQC` / *Degreaser*.

Two input changes matter more than the language:

* **It reads Thermo `.raw` directly** (`ThermoFisher.CommonCore.RawFileReader.dll`). No MGF or
  mzML export step for the MS2 side.
* **It reads Compound Discoverer `.cdResult` databases directly**, including CD's own
  retention-time alignment tables (`FileAlignmentCorrectionItems`, `CorrectionCurveItems`),
  rather than the CSV exports v1 parsed.

The second cuts both ways. It is a better interface to CD, and it is a *tighter* coupling to it:
searching `.cdResult` across the assemblies returns 30-odd hits and MZmine returns none, where
v1 shipped an `mzmine` package alongside `compound_discoverer`. **Compound Discoverer looks
required for the quantitative half of LipiDex 2**, where v1 at least had an alternative. No
command-line entry point was found either, so the headless argument for this project stands.

## Spectrum searching

| | v1 | v2 |
|---|---|---|
| mass tolerance | Da only, 0.01 default | **ppm or Da**, separately for MS1 / MS2 / MSn |
| MS levels | MS2 | **MSn, including MS3** — `MsnOrder`, child spectra, `MSnOrderMustMatch` |
| fragmentation type | ignored | **matched** — `UseAnyFragmentationTypeForMatching`, off by default |
| mass analyser | ignored | **matched** — `UseAnyMassAnalyzerForMatching` |
| retention time | ignored | **searchable window** — `MinRtSearchRange` / `MaxRtSearchRange` |
| dot product threshold | one global value | **per lipid class** — `UseClassDotProducts` |

Two of those are direct answers to problems documented elsewhere in this repo.

**Per-class dot products.** The UI text: *"phospholipids that have their MS2 spectra dominated
by the 184 m/z peak may want to have a lower dot product threshold to increase identification
rate."* That is exactly the single-fragment problem behind the 39 unresolvable ties in
[VALIDATION.md](VALIDATION.md) — v2 does not solve the ambiguity, but it stops one threshold
having to serve both a PC dominated by one headgroup peak and a TG with real chain fragments.

**Fragmentation-type matching**, with the help text *"Because fragmentation type greatly affects
the peaks in a library spectrum, it is strongly recommended to leave this unchecked"* — i.e. by
default an HCD library only matches HCD spectra. v1 matched anything to anything, which is the
same permissiveness as its never checking polarity.

The core scoring looks unchanged: `intWeight`, `massWeight`, forward and reverse dot product,
and the fatty-acid database `FattyAcids.csv` is **byte-identical to v1's**.

## Peak finder

New parameters, with the two substantial ones first.

**`UseRTLS` — a retention-time model per lipid class.** The description: *"Creates a retention
time model within each lipid class based on fatty acyl carbons and fatty acyl unsaturations to
assess identification likelihood."* Plus `FinalIDFromRTLS`, `RTLSRelationshipsFound` and a
*"maximum retention time error from retention time model"*.

This is the elution-order criterion **v1 does not have**. v1's only retention-time test on an
identification is `checkClassRTDist`, a clustering check asking whether a lipid sits in the
window where its class elutes; it never uses the fact that within a class retention rises with
carbon number and falls with double bonds. v2 does, and uses it to score identification
likelihood and to assign an ID outright.

**`CorrelationFilter` — abundance correlation across samples.** *"Filter redundant unidentified
features that have high correlation and co-elute with identified features"*, with a selectable
correlation type (Pearson among them). Our adduct sweep decides that two quant ions are the same
molecule from mass and retention time alone; v2 additionally requires them to rise and fall
together across the sample set, which is much harder to satisfy by coincidence.

Then:

* **Peak Quality Factors** from Compound Discoverer 3.3+, with a minimum weighted-PQF filter —
  *"the weighted score is a combination of the 4 base PQFs"*.
* `UseRtMadFactor` / `MaxRtMadFactor` — the class-RT MAD filter is now a real, toggleable
  parameter. In v1 the GUI spinner for it was **wired to an argument the method never read**
  (see the commit history here); v2 appears to have fixed that.
* `FwhmWindowMultiplier` — one explicit name, where v1 overloaded a single constant for both the
  ID-association window and the class-RT window.
* `GetConnectedRedundantCompoundGroups` — redundancy resolved over connected components rather
  than pairwise, which avoids the both-copies-removed failure we had to patch.
* `WriteGapStatus`, `WritePQFs` — gap-fill provenance in the output.
* Filters kept from v1, now with explicit help text: isotope (*"differ by 1.003 Thomsons, the
  m/z of a neutron"*), homo-dimer, in-source fragmentation, adduct.

## Adduct database, fixed

| v1 | v2 |
|---|---|
| `Name,Formula,Loss,Polarity,Charge` | `Adduct,GainFormula,LossFormula,Charge` |
| negative counts in formulas — `H-1`, `O-1H-1` | separate gain and loss formulas |
| polarity in its own column, unsigned charge | **signed charge** — `-1`, `-2`, `1` |

The signed charge is the same trap that made our negative mode return zero identifications: a
charge magnitude with polarity carried elsewhere is easy to lose. Contents differ slightly too —
v2 adds `[M-H-H2O]-` and drops `[M+K]+` and `[M+C2H6N2]+`.

## LipidQC, the "Degreaser" — entirely new

A QC suite with no v1 counterpart, driven by a metadata CSV that assigns each raw file a role:
**Sample**, **QcRep**, **QcDilution**.

* **%RSD** across replicate injections — *"Include at least 3 Quality Control Replicates"*.
* **Linear dynamic range** from a dilution series, with the sample quantitation distribution
  drawn as its marginal.
* **Peak quality factor distributions**.
* **Log2 quantitation distributions** per sample.
* **Retention time vs m/z modelling**.
* **Library-versus-experimental mirror plots**.
* **Run order and batch** handling, with run order recovered from file metadata so it survives
  copying between machines.

The role system is the same idea as `blanks.py` here, arrived at independently — though note
LipiDex 2's roles are Sample / QcRep / QcDilution, with **no blank role**, so blank filtering
still appears to be the analyst's problem.

## Libraries

| | v1 | v2 |
|---|---|---|
| HCD Formic | 45.3 MB | **56.8 MB** |
| HCD Acetate | 35.4 MB | 44.0 MB |
| LipidBlast Formic | 22.0 MB | 28.0 MB, renamed `InSilico_LipidBlast` |

New: **CID MS2 and CID MS3** libraries, **FAHFA**, **Gangliosides**, **Ornithine lipids**, and an
`UltimateSplash_ISTD` set.

## What is worth taking

In rough order of value to this project:

1. **RTLS.** A per-class retention model on carbons and unsaturations is the strongest orthogonal
   criterion available and we have everything needed to fit one. It would attack the residual
   duplicate-retention-time problem directly.
2. **Correlation filter.** We already carry per-sample areas for every compound group, so
   requiring a correlation before calling two ions the same molecule is a small change to the
   adduct sweep and a large gain in specificity.
3. **Per-class dot product thresholds.** Cheap, and aimed exactly at the classes where our
   single-threshold behaviour is worst.
4. **Fragmentation type and mass analyser in matching.** We ignore both, as v1 did.
5. **ppm tolerances** rather than fixed Da.

Not worth copying: the tighter binding to `.cdResult`. Removing Compound Discoverer is the point
of this project, and on that axis v2 moves the wrong way.
