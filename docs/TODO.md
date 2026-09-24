# Known defects to fix

Ordered by how much delivered data they affect. Each one has been seen on real studies, with the evidence recorded.

---

## 1. M+2 isotopes of unsaturated species are promoted as more-saturated species ⚠ largest

**What happens.** One extra double bond is −2.01565 Da; two ¹³C is +2.00671 Da. They differ by only **8.9 mDa (~11 ppm at m/z 800)**, so the M+2 isotope of `PC n:d+1` lands almost on the monoisotopic mass of `PC n:d`. When that isotope peak becomes its own feature and is matched to the library, a **fake, more-saturated species** is created.

**Scale, on `20260921_Chinmay_Lipidomics`:**
- **27 of 214 PC features, carrying 26% of the total PC signal.** Also SM (8 features, 26% of signal) and TG (8, 10%).
- **Positive mode only** in that run; every negative-mode class was clean (PE, PI, PG, PS, plasmenyl-PE: 0 suspects).
- It biases composition, not just single features: PC's double-bond index changed from −0.021 to **+0.051** (sign flip) and the WT PUFA share from 16.1% to 19.6% once they were removed.

**How to detect it (all three agree, and mass accuracy is the strongest):**
1. **Mass accuracy.** Real pairs one double bond apart sit within **±1 ppm** of the 2H expectation (measured on PE in this run). The artefacts sit **−10 to −17 ppm** from it, i.e. far outside instrument performance, and within ~1 mDa of the 2 × ¹³C spacing.
2. **Co-elution.** A genuine extra double bond elutes **~21 s earlier** for PC on this method (PE 19 s, TG 31 s, PI 18 s). The artefacts co-elute to **within 1 s**.
3. **Abundance.** They sit at roughly the theoretical M+2 share (10–18% of the partner).

**Suggested fix.** Before promoting or keeping a feature, test it against every co-eluting feature of the same class with one more double bond: if the mass gap is closer to 2 × ¹³C than to 2H **and** the apex RT difference is far below the class's measured per-double-bond shift, drop it (or flag it) whatever its library match. Do not gate this on whether the feature has an identification — see defect 3.

**Evidence:** `/mnt/datastore/Jair/Projects/JMJD5/Analysis_Multiomics_2026/01_lipidomics/` → `isotope_suspect_features.csv`, `isotope_suspects_by_class.csv`, `UNSATURATION.md`, `pc_annotation_check.py`. Raw-file confirmation pending.

---

## 2. `recover_targeted_features` is fed per-file MS2 targets → fake on/off features

`pipeline.py` builds the recovery target list **per injection** (`targets[Path(mzml).stem] = ...`), so a peak is only rescued in the files where DDA happened to fragment and match it. The same peak is left as an isotope elsewhere and reported as **zero**, which reads as a biological on/off difference.

**Case:** `PA 44:6` (Neg, m/z 803.56, RT 8.73) in `20260921_Chinmay_Lipidomics` was reported 0 in all three WT and ~1.1e7 in all three KO. Raw files showed the peak present in WT at 9.3–9.9e6 — a real ratio of ~1.2×, not on/off. It is the ¹³C isotope of **PC 34:2 [M+HCOO]⁻**; the run's −5 ppm offset put it on the PA 44:6 mass.

**Suggested fix.** Pool the targets across the study, or after recovery re-measure every rescued m/z and RT in all files. **Every on/off feature in past deliveries should be re-checked.**

---

## 3. The M+1 isotope filter skips any peak that carries an identification

`peakfinder.py` removes an M+1 isotope only when `heavier.final_lipid_id is None`. A library-matched peak therefore survives even when it sits exactly one isotope spacing above a larger co-eluting peak and tracks it across samples. This is what lets defects 1 and 2 reach the output.

**Suggested fix.** Apply the isotope test regardless of identification; where it fires on an identified peak, flag it rather than silently keeping it.

---

## 4. Isotopes of a class kept in the other polarity are not recognised

The adduct/isotope filter missed a negative-mode isotope of **PC**, a class this pipeline keeps in positive mode only, so it was never compared against its own parent. Isotope and adduct relationships should be tested against all detected features in that polarity, not only those of classes retained there.
