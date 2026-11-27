# Prenol quinones — ubiquinone and ubiquinol, positive mode — provenance

Built by `scripts/build_ubiquinone_library.py`, which is the reproducible record.

    python scripts/build_ubiquinone_library.py \
        --out data/libraries/Ubiquinones_Positive.msp \
        --anchors data/libraries/Ubiquinones_RT_anchors.csv

**4 target entries, no decoys.** Added 2026-10-02. No third-party file is read: precursor and
fragment masses are computed from elemental formulae, fragment intensities are measured on this
laboratory's own spectra. MIT, as the rest of our own material (`NOTICE.md`).

| entry | precursor | retention | spectra pooled |
|---|---|---|---|
| `CoQ10 [M+NH4]+` | 880.7177 | 21.07 | 12 |
| `CoQ10 [M+H]+` | 863.6912 | 21.07 | 10 |
| `CoQ10H2 [M+NH4]+` (ubiquinol) | 882.7334 | 19.50 | 11 |
| `CoQ9 [M+NH4]+` | 812.6551 | 19.10 | 7 |

Neutral formulae: CoQ10 `C59H90O4` 862.6839, ubiquinol-10 `C59H92O4` 864.6996, CoQ9 `C54H82O4`
794.6213. `--predicted-homologues` adds CoQ6–CoQ11 and their quinols by n x C5H8 (68.0626); they
are predictions, off by default, and only CoQ9, CoQ10 and ubiquinol-10 are confirmed here.

## Why it exists

Prenol lipids are an entire LIPID MAPS category — PR — that no library in the default set covers.
Zero quinone entries anywhere, which is why coenzyme Q has never appeared in a delivered table
from this pipeline or from the Compound Discoverer + LipiDex workflow it replaces; the
an in-house study's 940-lipid matrix has no quinone row for the same reason. The category has
had a colour reserved in `qc_report.py` since that plot was written and nothing has ever filled
it.

It is not absent from the samples. On an in-house Q Exactive Plus study, positive mode, ubiquinone is a
**better peak than the free-cholesterol positive control**: max/median 812 against cholesterol's
45, 222 non-zero scans, present in all ten injections examined, with `[M+H]+` and `[M+Na]+`
agreeing at the same retention.

## Source of the values

Our own data, 2026-10-02:
a facility re-analysis (raw data not distributed)

| script | what it contributed |
|---|---|
| `cq01_xic_pos.py` | targeted MS1 XIC, ±5 ppm, 10 injections: the precursors, the retentions, and the demonstration that the class is there at all |
| `cq03_ubiquinol_check.py` | the ubiquinol-vs-isotope test — co-elution profiles, adduct agreement, MS2 head group |
| `cq04_coq_ms2_peaklist.py` | the fragment intensities written into the library |

Study: an in-house Q Exactive Plus study, Thermo Q Exactive Plus, facility reversed-phase gradient,
ammonium formate. **Every retention here belongs to that gradient and to no other.**

## The intensities are measured; the masses are computed

`cq04` takes every MS2 whose precursor is within 20 ppm of a target, whose retention is inside
that target's window, **and whose base peak is the quinone head group 197.0808**. That last gate
is not fussiness: without it a minority of co-isolated spectra contributed 184.0733
(phosphocholine) and a 292/293/320/321 series at median intensity above 900, which is a
neighbour's spectrum, not a weak CoQ fragment. Each qualifying spectrum is scaled to base peak
999, fragments are binned at 15 ppm, and a fragment is written if it is present in at least 70 %
of the spectra for that adduct, at its median intensity.

Two fragment families, both generic to the class rather than to the homologue:

* **the quinone head** — 197.0808 `C10H13O4+`, the base peak of every qualifying spectrum
  (100 % presence), with its losses 182.0574 (−CH3·), 169.0859 (−CO), 154.0624 and 151.0754;
* **the polyprenyl chain** — the CnH(2n−1)+ / CnH(2n−3)+ hydrocarbon ladder, 67 to 189.

Masses are written from the assigned formula, not from the measurement: the measured low-mass
ions run +10 to +16 ppm on this instrument, about 1 mDa at m/z 81. That is well inside `MS2_TOL`
(0.01 Da) and there is no reason to write the instrument's drift into a library.

⚠ `184.0733` is **deliberately absent** although it was measured at median 395 in 7/7 CoQ9
spectra, and `tests/test_ubiquinone_library.py` enforces that. It is a co-isolated PC/SM head
group; writing it would make these entries score against every PC spectrum in a run.

## ⚠⚠ The trap: ubiquinol is 10 ppm from the ubiquinone M+2 isotope

    CoQ10   [M+NH4]+  M+2 isotope   882.7244       (880.7177 + 2 × 1.0033548)
    CoQ10H2 [M+NH4]+  monoisotopic  882.7334
                                    --------
                                    0.0090 Da  =  10.2 ppm at m/z 882

`search.MS1_TOL` is **0.01 Da, absolute**, so the isotope falls inside the window that selects
library candidates: the ubiquinol entry is a candidate for every MS2 taken on the ubiquinone
isotope. MS2 cannot break the tie — both give 197.0808 as base peak, because it is the same
molecule's head group, and the isolation width is roughly 100× the mass difference. Nothing in
the peak finder catches it either: its isotope rule fires at M+1 spacing and only on a peak with
no identification (known limitations of the peak finder: isotope promotion, and an isotope rule applied only to unidentified peaks), and an identified ubiquinol row is precisely
the case it skips.

**They separate on retention and on nothing else**, measured over three files:

    CoQ10   [M+NH4]+  880.7177   apex 21.05–21.09
    CoQ10   M+2       882.7244   apex 21.05–21.08    contributes ZERO at 19.5
    CoQ10H2 [M+NH4]+  882.7334   apex 19.50–19.60    ~1.5 min earlier
    CoQ10H2 [M+Na]+   887.6888   apex 19.52–19.56    the adduct agrees; an isotope cannot do that

The reduced form carries two hydroxyls where the oxidised form carries two ketones, so it is more
polar and elutes earlier on reversed phase. The direction is chemistry, not an artefact.

The `.msp` format has no retention field, so the constraint is carried in three places instead:

1. `Ubiquinones_RT_anchors.csv` — the measured retention of every ion, machine-readable;
2. `scripts/check_ubiquinol_rt.py` — tests a delivered table against it and names any ubiquinol
   row that co-elutes with ubiquinone as the isotope it is. Exit status 1, so it can gate a
   delivery. **Run it on any study that reports ubiquinol.**
3. `RT=` and `RTWindow=` in each entry's `Comment`, so the number travels with the entry. No
   reader consumes them; they are for whoever opens the file.

## Known limits

- **Not in `DEFAULT_LIBRARIES`.** A library enters the default set only after its displacement
  against matched decoys has been measured (the developer notes (LIBRARY_VERSIONS, not distributed)); `Cer[EOS]` and `FAHFA`
  were both withdrawn on that test. These entries also cannot take the CH2 decoy construction —
  `make_decoy_library.py` shifts chain fragments identified by their `[c:d]` annotation, and a
  quinone has no acyl chain — so a decoy for them needs its own construction first. Searching
  them meanwhile is `--extra-library Ubiquinones_Positive`, which is the mechanism built for a
  library that must not sit in the defaults; see the developer notes (PRENOL_QUINONES, not distributed).
- **Purity is always 0.** There is no acyl chain to count, the same as `Sterols_InSilico`. The
  names are defined structures rather than sum compositions, so nothing is lost by it, but a
  purity filter will drop these entries if one is applied without thinking.
- **No authentic standard.** Everything here is measured on endogenous material. The reduced head
  group 199.0965 `C10H15O4+`, which would be positive evidence for ubiquinol specifically, was
  looked for in all 11 spectra and observed in **none** — the quinol's head oxidises in the
  source, so its spectrum is the quinone's. A ubiquinol standard is what would settle the
  identification formally.
- **The quinol/quinone ratio is not a redox read-out.** This extraction has no antioxidant, so
  most of what the ratio reports is auto-oxidation during extraction (measured 0.077 → 0.19
  across samples of one study). The species are measurable; the redox state is not, without a
  redox-preserving extraction.
- **No `[M+Na]+` entries.** Both forms give sodiated ions in MS1 at the same retention (885.6731,
  887.6888) and they are in the anchors table for targeted work, but not one sodiated MS2 passed
  the screen, so there is nothing measured to write.
- **DMQ10 is not here.** The COQ7 substrate was searched for and is not reliably detected
  (max/median 15.6, retention scattered by 3 min — noise). It needs a targeted method.
