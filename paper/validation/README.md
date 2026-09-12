# Validation table — how each number in it is produced

Every field is read from the deposit or from a run. Nothing here is typed by hand, and that is a
rule with a reason: this file's predecessor hardcoded MSV000094718 as *Mus musculus* when it is
*Bos taurus* — bovine liver, NCBI TaxID 9913, stated plainly in the MassIVE record. The field
looked too small to be worth sourcing.

| script | what it does |
|---|---|
| `build_table.py` | the table itself — one row per deposit, every species carrying a `species_source` |
| `score_agreement.py` | our identifications against the ones each deposit published, on LSI shorthand |
| `audit_methods.py` | what each deposit ACTUALLY acquired, sampled per file |
| `compare_st000991.py` | one dataset's old run against its new one |

## Three distinctions the table has to preserve

**Injections are not samples.** Every count is FILES, which include pools, blanks and technical
replicates. Where a deposit states a biological n it is carried separately, and the two disagree
often enough that reporting one as the other would be wrong on most rows.

**"No names deposited" is not 0% agreement.** Six deposits publish `m/z_RT` feature tables with no
identification column, four of them by design as data-PROCESSING benchmarks. A blank cell and a 0%
cell make opposite claims about the software.

**Accessions are not studies.** ST004503+ST004626 are one benchmark; ST004650+ST004651 are the same
benchmark on HILIC; ST000987+ST000991 are parts V and IX of one nine-platform comparison sharing
153 samples. Counting accessions overstates n. `assay_scope` carries the related awkwardness: all
four CSU deposits are titled *metabolomics*, two of them HILIC, and a lipid library against HILIC
metabolomics is not a lipidomics validation whatever it scores.

## What agreement means, and does not

It is deliberately not called recall. Recall implies the published list is truth, and it is not —
it is another software's output, made with another library, at another evidence threshold. Two
pipelines disagreeing tells you they disagree.

Spiked internal standards are excluded from the denominator (Workbench marks them
`refmet_name = "Standard"`; raw names carry SPLASH/ISTD/(d7) tags). They are compounds the lab
added, so counting them would penalise us for correctly declining to report them as endogenous.

⚠ Two deposits cannot be compared at all, for reasons that are not our software's: ST000987
fragments without isolating (all-ion, in-source CID) and MTBLS2016 contains no MS2 despite being
described as data-dependent. Their published names come from MS1 mass and retention time. Scoring
against them measures a difference in evidence policy, not in performance.
