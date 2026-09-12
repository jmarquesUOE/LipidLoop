# What the MS2 duty cycle is being spent on

A DDA run has a fixed budget. These files fit about nine MS2 into a 0.87 s cycle, which is why
there are only ~1.2 MS1 scans per second and roughly five points across a 4.7 s peak. Every MS2
spent on column bleed is one not spent on a lipid, and it is paid for twice — once in the
identification that was not attempted, once in the MS1 sampling that makes the peak thinner.

    python scripts/find_wasted_ms2.py \
        --mzml data/mzml/val_Pos/*.mzML --results outputs/validation_Pos \
        --blank data/mzml/val_Pos/Blank_0*.mzML --polarity + \
        --out method/Exclusion_List_Pos.csv \
        --thermo method/Exclusion_List_Pos_Thermo.csv \
        --report results/wasted_ms2_Pos.md --simulate

Three outputs, because they have three different readers. `--out` is the annotated list, one row
per candidate carrying the evidence — **the one to review**. `--thermo` is the same rows in the
Xcalibur mass-list import layout. `--report` is the whole analysis in markdown, every table plus
the per-file breakdown and enough provenance to say what it was run on.

The current lists live in [`method/`](../method/README.md) and the current reports in
[`results/`](../results/). Both are tracked in git; `data/` and `outputs/` are not.

## What it found here

Counting only precursors that are **never identified in any file**, **elute across more than half
the gradient**, and are **equally present in the blanks**, over 24 files per polarity:

| | files | precursors | MS2 spent on them | share of all MS2 |
|---|---|---|---|---|
| positive | 24 | 285 | 166,497 / 386,881 | **43.0%** |
| negative | 24 | 343 | 200,661 / 299,574 | **67.0%** |

24 files, not 26: the two blanks per polarity are used as the blank test and are not counted as
samples. Counting them as both inflates the duty-cycle share, because a blank is almost nothing
but contamination.

The four QC injections alone give 44.6% and 67.4% — the same answer from a sixth of the data, so a
short pilot is enough to build the list.

Broken down by homologous series, which is how you tell where it is coming from:

| polarity | series | repeat | precursors | share of duty cycle | source |
|---|---|---|---|---|---|
| + | **polysiloxane** | 74.019 (C₂H₆OSi) | 67 | **24.5%** | column bleed, PDMS from tubing and septa |
| + | unassigned | — | 173 | 12.6% | |
| + | PEG | 44.026 (C₂H₄O) | 18 | 2.8% | detergent, plasticiser, plastic labware |
| + | sodium formate | 67.987 (CHO₂Na) | 24 | 2.0% | formic acid + sodium from glass |
| + | PTFE/PFPE | 99.994 (C₂F₄) | 3 | 1.1% | fluoropolymer tubing, pump oil |
| − | **sodium formate** | 67.987 (CHO₂Na) | 214 | **44.6%** | formic acid + sodium from glass |
| − | unassigned | — | 122 | 15.0% | |
| − | PTFE/PFPE | 99.994 (C₂F₄) | 4 | 7.4% | fluoropolymer tubing, pump oil |
| − | polysiloxane | 74.019 (C₂H₆OSi) | 3 | 0.0% | |

The split is almost perfectly complementary: siloxane dominates positive mode and is invisible in
negative, salt clusters dominate negative and are a rounding error in positive. Two different
problems wearing the same costume.

**Nearly half the negative-mode duty cycle goes to one salt-cluster series, and two thirds
to background overall.** The Δ 67.987 series is
[(HCOONa)ₙ + HCOO]⁻ — the mobile phase reacting with sodium from glassware — and its members run
from m/z 316.9451 up past 1500. That one is worth fixing at the bench as well as excluding: less
sodium in, fewer clusters to exclude.

## Would it actually improve identification?

Yes, and the size of the gain can be measured rather than assumed. `--simulate` detects MS1
features, finds the ones that never received an MS2 at all, and walks the run cycle by cycle
handing each freed slot to the most abundant unfragmented feature co-eluting at that moment:

| | features | never fragmented | reachable with the freed slots |
|---|---|---|---|
| positive, QC_01 | 5,200 | 2,916 (56%) | **2,437 (84%)** |
| negative, QC_01 | 1,253 | 376 (30%) | **335 (89%)** |

The starved-feature half of this is not a QC artefact: a rat heart sample from the same batch has
5,025 features with 2,762 (55%) never fragmented, the same picture.

Even among the 500 most abundant positive-mode features, 52 never got an MS2; among the top
1,000, 173. Those are not marginal signals — the median feature the freed slots would reach has an
intensity of 6.4×10⁵, and the largest 1.4×10⁸.

**This is an upper bound.** A real instrument re-ranks by intensity live and applies its own
dynamic exclusion, so it will not follow the simulation exactly, and some freed slots will land on
isotopes, adducts, or other things already covered. The honest claim is that the budget exists and
something worth fragmenting is waiting for it — not that identifications rise by 88%.

There is a second effect that costs nothing to take. Whatever is not spent on MS2 is spent on MS1:
at ~9 MS2 per cycle, removing 45% of them shortens the cycle and lifts the MS1 rate. Five points
across a peak is thin for any peak picker, and this is the cheapest way to get more.

## Why this is safe, and where it stops being safe

Excluding a precursor is the only decision in this pipeline that reprocessing cannot undo. The
artefact screen can be turned off and the data re-searched; a spectrum that was never acquired is
gone. So the criteria are chosen against the failure mode that matters — excluding a real lipid.

**"Never identified" is a bad criterion on its own.** Most unidentified spectra are real ions the
library does not cover, and excluding them would freeze the method's blind spots in place: the
lipid classes we cannot yet name would become the ones we can never discover. Free cholesterol was
missing from every library in the default set until recently; on a "never identified" rule it would
have gone on the exclusion list.

**Chromatography is what separates the two.** A lipid elutes as a peak, a few seconds wide.
Contamination is present at every retention time because it enters the source continuously. That
is the load-bearing criterion — the identification test is only there to make sure we are not
excluding something the pipeline is currently getting right.

Both are backed by the blank. Sample matrix is not required to produce any of these ions.

### The limits, stated plainly

- **A compound that genuinely elutes throughout would be excluded.** Free fatty acids tailing
  across a whole gradient on a badly conditioned column are the realistic case. Look at the list
  before using it; anything with a plausible lipid formula deserves a manual check.
- **The series labels are corroboration, not proof.** They are assigned by repeat spacing, and a
  spacing can coincide. Multiples of CH₂ are never proposed as an unnamed series precisely because
  that is what a *lipid* homologous series looks like — PC 32:0 / 34:0 / 36:0 are 14.0157 apart.
- **`--results` matters.** Without it every precursor counts as unidentified and only the elution
  test is doing any work. The script says so when you omit it.
- **The list ages.** Column bleed depends on the column, PEG on the labware, sodium formate on the
  glassware. Re-run it per batch, exactly like the artefact screen.

## Relationship to the other two lists

Three separate things, easy to confuse:

| | what it does | when | reversible |
|---|---|---|---|
| [artefact screen](../src/lipidloop/artefacts.py) | removes a *fragment* m/z from acquired spectra | processing | yes |
| **exclusion list** (this) | stops a *precursor* being fragmented | acquisition | **no** |
| [CE inclusion list](CHOLESTERYL_ESTERS.md) | forces a precursor to be fragmented | acquisition | n/a |

## The trap the collision check actually caught

Worth stating plainly, because it nearly went the other way. On the skin set the generated
negative-mode list contained **m/z 255.2327 (2,679 MS2), 283.2641 and 281.2483** — palmitate,
stearate and oleate. Excluding those would have permanently prevented ever fragmenting the three
commonest fatty acids in biology.

They qualified because they meet every criterion honestly. A free fatty acid **is** never
identified by MS2 — it does not fragment ([FATTY_ACIDS.md](FATTY_ACIDS.md)). It **does** elute
across the whole run, because its mass is produced in-source by every phospholipid carrying that
acyl chain. And it **is** in the blanks. *"Never identified, everywhere, in the blanks"* describes
contamination and it equally describes a class this pipeline could not yet name.

Checking against `--results` catches it only where the run already names them, and the list built
before the fatty-acid work would have carried them straight through. So free-fatty-acid masses are
now **protected in negative mode whether or not the run identified any** — the mass table, not the
result file, is the guard. On the rat heart set, where the fatty-acid series is too sparse to name
and nothing protected it, `FA 23:0` was on the list until this was added.

**The general lesson: a class you cannot identify looks exactly like contamination to a screen
built on "never identified".** Any future library that closes a gap should be followed by
regenerating and re-checking these lists.

## One thing the file format cannot tell you

**Check the list against real lipids before importing**, at the exclusion mass tolerance your
method uses. An exclusion entry does not remove one m/z — it removes everything the instrument
cannot tell apart from it, and that width is a method setting the list file has no column for.
`scripts/check_exclusion_safety.py` does the check:

    python scripts/check_exclusion_safety.py --list method/Exclusion_List_Pos.csv --polarity + \
        --results outputs/validation_Pos/Final_Results.csv --ppm 10 \
        --report results/exclusion_safety_Pos.md \
        --safe-list method/Exclusion_List_Pos_Safe.csv \
        --safe-thermo method/Exclusion_List_Pos_Safe_Thermo.csv

Two collisions are possible and they are different failures. A lipid **within the exclusion
tolerance** matches the entry and is never fragmented — that is the one that costs an
identification. A lipid **within the isolation window** (±0.55 Da here) still gets selected on its
own m/z and still gets an MS2; the contaminant is merely co-isolated into it. That second one is
worth knowing and is *not* a reason to drop an entry: the co-isolation happens today, with or
without the list, and excluding a precursor does not put it into anyone's isolation window.

At ±10 ppm on this data:

| | entries | collides with a lipid identified here | collides with any library lipid | co-isolates (±0.55 Da) |
|---|---|---|---|---|
| positive | 285 | **0** | 23 | 77 |
| negative | 343 | **0** | 12 | 98 |

**Nothing on either list would cost an identification made in this experiment.** Dropping every
entry within 10 ppm of *any* library lipid as well — species this batch never contained but another
tissue might — still recovers 98.4% (positive) and 99.7% (negative) of the duty cycle the full list
recovers. The conservative list is essentially free, so use it: `Exclusion_List_*_Safe.csv`.

The exclusion and inclusion lists are complementary and belong in the same method: the exclusion
list frees the budget, the inclusion list decides where some of it goes. Cholesteryl esters are the
worked example — only one CE was identified because DDA never picked the others, and the freed
slots are exactly what an inclusion list needs to spend.

## `Method_Lists/`: the import-ready lists, written on every run (added 2026-09-11)

Everything above was a by-hand sequence: generate the list, run `check_exclusion_safety.py`,
convert the inclusion list, collect the files. Jair, 2026-09-11: "1 2 3 should be hardcoded in the
software." So `run_study.py` now calls `lipidloop.method_lists.write_method_lists()` after both
polarities have run and before the report, and writes `<analysis>/Method_Lists/`:

| file | what |
|---|---|
| `Exclusion_List_{Pos,Neg}_Thermo.csv` | the exclusion candidates that survive the collision check, Xcalibur mass-list layout, Start/End empty (whole run) |
| `Inclusion_List_{Pos,Neg}_Thermo.csv` | the inclusion candidates in the same layout, Start/End = apex ± 0.5 min |
| `Exclusion_Safety_{Pos,Neg}.md` | what was dropped and which lipid it collided with |
| `README.md`, `summary.json` | counts and the tolerances assumed |

The check is the same as `check_exclusion_safety.py`: ±10 ppm (the exclusion tolerance a method
typically uses) against every non-decoy lipid identified in `Final_Results.csv` and every lipid
in the libraries the run searched (from `config_<pol>.json`, decoy libraries skipped), with the
free-fatty-acid masses protected in negative mode whether or not the run named any. The one thing
the file still cannot express is the method's own exclusion tolerance: if it is wider than
±10 ppm, re-run `scripts/build_method_lists.py <analysis> --ppm N`. That script also builds the
folder for analyses produced before this existed.
