# Why four identified lipids have no feature

Traced on the CKD rat heart set, positive mode. Four of the fourteen molecules the reference
reports and we do not are cases where **the identification stage found the lipid confidently but
there is no MS1 feature for it to attach to** — `SM d41:0` from 43 spectra at dot 987, and no
compound group within 20 ppm and 0.3 min.

The MS1 signal is not the problem. All four are far above the noise floor:

| | MS1 points | apex | M+1 | noise floor |
|---|---|---|---|---|
| SM d41:0 | 37 | 5,591,052 | 943,073 | 5,000 |
| SM d42:0 | 59 | 11,237,934 | 1,873,954 | 5,000 |
| PC 48:5 | 9 | 16,975 | 9,372 | 5,000 |
| AC 24:2 | 5 | 22,312 | 7,769 | 5,000 |

Following them through the three detection stages gives three different answers.

## AC 24:2 — never traced

`MassTraceDetection` produces no trace at all. Five MS1 points at an apex of 22,000 on a
1.2 Hz duty cycle is close to the floor of what a trace can be built from.

**Recoverable** by dropping `noise_threshold` from 5,000 to 2,000 — at the cost of **11× more
features**, 5,379 to 57,934 per file.

## PC 48:5 — traced, peaked, then discarded for having no isotope

A trace exists and `ElutionPeakDetection` finds a peak at RT 11.86, FWHM 4.3 s. It is then
dropped by `remove_single_traces`, which requires a feature to have more than one isotope
trace. Its M+1 sits at 9,372 — above the floor but evidently not traced consistently enough
across the peak.

**Recoverable** by `remove_single_traces=false` — at the cost of **6× more features**, 5,200 to
30,981.

## SM d41:0 and SM d42:0 — the trace is consumed as somebody else's isotope

This one is not a threshold and cannot be tuned away.

At RT 10.98 there is a feature at m/z 801.6915 carrying five traces:

    801.6915   802.6902   803.6912   804.6944   805.6981

That third trace is where `SM d41:0` elutes. It has been assigned as the **M+2 isotope of the
mono-unsaturated homologue**, and a trace can belong to only one feature.

The arithmetic is why. One double bond is 2.01565 Da; two ¹³C are 2.00671 Da. They differ by
**8.9 mDa, which at m/z 803 is 11 ppm** — so the monoisotopic peak of the saturated species sits
almost exactly on the M+2 isotope of the species with one more double bond. For a lipid class
whose members differ by single double bonds, that collision is systematic rather than unlucky,
which is why it takes out `SM d41:0` and `SM d42:0` together.

Three isotope-scoring settings were tried and none recovers them:

| | features | SM d41:0 | SM d42:0 |
|---|---|---|---|
| baseline | 5,200 | no | no |
| `mz_scoring_13C=true` | 5,064 | no | no |
| `mz_scoring_by_elements=true` | 5,301 | no | no |
| `remove_single_traces=false` | 30,981 | no | no |

A separate feature does exist at m/z 803.7054, but at RT 11.39 — too far from the MS2 at 10.97
to associate. Whether the `SM d41:0` assignment is itself sound is worth a second look: the MS2
precursor at 803.695 falls between the M+2 isotope (803.692) and the saturated species (803.701),
and an isolation window takes both.

## What this costs and what it would cost to fix

The current settings favour precision. Recovering `AC 24:2` and `PC 48:5` means 6–11× more
features per file, which is slower at every downstream stage and gives the adduct, blank and
in-source filters far more to chew on — the extra rows are mostly single-trace noise, which is
what `remove_single_traces` exists to remove.

## The targeted recovery pass

Implemented, in `recover_targeted_features`. Identification runs before feature detection in the
pipeline precisely so that this is possible: after assembly, every elution peak is checked
against the identifications that passed the score thresholds, and any peak that **an MS2 has
already identified** but that **no feature represents** is promoted to a feature of its own,
carrying the trace's own peak area and width.

The safeguard is the ordering, not a threshold. Nothing is recovered without an identification
standing behind it, so the pass cannot manufacture features out of noise; the worst it can do is
promote a peak that really was an isotope, at a mass and retention time where a spectrum matched
a library entry anyway.

It costs about **2.5% more features** — 112 to 162 per file against ~5,000 — where the parameter
route costs 6 to 11 times more.

| | before | after |
|---|---|---|
| positive, molecules recovered | 96.7% | **97.7%** |
| positive, genuinely absent | 14 | **10** |
| negative, molecules recovered | 94.3% | **95.3%** |
| negative, genuinely absent | 11 | **9** |

Brought back: `SM d41:0`, `SM d42:0` — the two the pass exists for — plus `PC 48:5`,
`PE-NMe 38:4` and `Plasmanyl-PC O-30:0`.

`AC 24:2` is still absent and always will be by this route: it is never traced at all, so there
is no elution peak to promote.

**The extra identifications are not validated.** Molecules only we report went from 63 to 80 in
positive mode, concentrated in SM and Cer[NS] — exactly the classes where the isotope collision
happens, which is consistent with recovering real saturated species, and equally consistent with
promoting isotope peaks that a spectrum happened to match. Distinguishing the two needs manual
inspection and has not been done.

## What remains, positive mode

Ten molecules, and no longer one pattern: two deliberate `D5TG` exclusions, three removed by our
own blank and class-retention filters and so pending those decisions, two halves of the O-/P-
ether tie, one below the reverse dot threshold, one never traced, one new.

## Can the thresholds be tuned?

Yes, and `scripts/tune_features.py` does it — but the answer it gives is *don't*.

**The objective matters.** Optimising to match Compound Discoverer would bake its choices in;
it is a comparator, not ground truth. The identifications are the better target: every MS2 that
clears the score thresholds is an independent statement that something real elutes at that mass
and time, so a parameter set is scored by how many of them get a feature. That is self-contained
and measures exactly the failure mode above. The recovery pass is switched off during tuning,
since it exists to paper over detection misses and would hide what is being measured.

Sweeping noise threshold, SNR, minimum trace length and `remove_single_traces` over three files
— positive QC, a positive sample, negative QC — gives a consistent picture:

| setting | features | detection recall |
|---|---|---|
| noise 5000, `remove_single_traces=true` (current) | 5.2k / 3.8k / 1.3k | 74.0% / 81.2% / 84.3% |
| noise 10000, `remove_single_traces=false` | 17.2k / 11.7k / 3.1k | **84.5% / 89.1% / 92.5%** |

`remove_single_traces` is the dominant lever, worth 8–10 points on its own — far more than the
noise floor or the SNR — and raising the noise floor pays for most of the features it costs.
Everything else is nearly flat: dropping the noise floor from 10,000 to 1,000 buys two points of
recall for four times the features.

**And it was tried end to end, and rejected.** Detection recall 74% → 84% produced:

| | current | tuned |
|---|---|---|
| molecules recovered | 97.7% | 98.1% |
| output rows | 4,929 | **11,109** (reference 4,128) |
| runtime | ~430 s | 621 s |

Two more molecules for 6,181 more rows — about 3,000 rows per molecule. The extra features are
overwhelmingly *unidentified*: identified groups rose from 835 to 1,099 while distinct molecules
went from 528 to 541. The targeted recovery pass gets the same identified peaks for **2.5% more
features rather than 170%**, because it asks for the specific peaks that are missing instead of
lowering the bar for everything.

So the tuning is worth running — it is how the `remove_single_traces` effect was found, and it
would matter on data where recovery could not help — but on this data the recovery pass
dominates loosening the thresholds, and the defaults stay where they are. That is the reason
they are where they are, rather than an accident.
