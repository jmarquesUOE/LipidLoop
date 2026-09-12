"""Our measured retention times against the ones a deposit published — the first external check.

The per-class retention model has so far only ever been judged against itself: it flags a lipid
whose retention is impossible for its class, and the only confirmation available was internal
consistency. MTBKS222 changes that. Its MAFs carry a retention time for every one of the 589
identifications in the Thermo arm, measured on the same injections we processed, so for the first
time there is an outside answer to compare against.

Two questions, and they are different:

**Do we put the same lipid at the same time?** Same raw files, same chromatography, so a
disagreement is a disagreement about which peak a name belongs to — not about method. This tests
peak assignment, and it is the one place where a name-level agreement score is blind: two pipelines
can "agree" on `PC 34:1` while pointing at different peaks.

**Does the model's own outlier flag predict those disagreements?** If the lipids we flagged as
retention outliers are the same ones the deposit places elsewhere, the flag is measuring something
real. If they are unrelated, it is measuring noise, however internally consistent it looked.

⚠ Matching is on the LSI shorthand, never on the raw name. The deposit writes `Cer-NS d18:1/17:0`
and we write `Cer[NS] 18:1;O2/17:0`; comparing raw strings would score a near-total miss and read
as a chromatography problem.
"""
from __future__ import annotations

import csv
import statistics
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import deposits                                    # noqa: E402
import score_agreement                             # noqa: E402

STAGED = Path("/home/jair/validation_staged")
# Their MAF is one row per adduct, so the same lipid appears several times at one retention time.
AGREE_MIN = 0.2          # minutes; inside this the two are calling the same peak


def _key(name: str) -> str:
    """The sum-composition key, borrowed from the scorer rather than reinvented.

    ⚠ A first pass here matched on LSI shorthand alone and shared only 74 lipids where the
    scorer shares 180. The reason is granularity, not spelling: the deposit writes
    `Cer-NS d18:1/22:0` and we write `Cer 40:1`, so a shorthand comparison drops every class
    where the two sides resolve chains differently — 98 TGs, 32 SMs, every ceramide. Summing to
    a species key is what makes them comparable, and `score_agreement.canon` already does it.
    """
    molecular, species = score_agreement.canon(name)
    return species or molecular


def published_rt(dataset: str, arm: str) -> dict[str, float]:
    out: dict[str, list[float]] = {}
    for row in deposits.maf_rows(dataset, arm):
        name = (row.get("metabolite_identification") or "").strip()
        rt = (row.get("retention_time") or "").strip()
        if not name or not rt or name.upper().startswith("UNKNOWN"):
            continue
        key = _key(name)
        if not key:
            # ⚠ canon() returns "" for a name it cannot parse. Left in, the empty key matches
            # the empty key on the other side and appears in the shared set as a lipid with no
            # name and a real-looking retention difference.
            continue
        try:
            out.setdefault(key, []).append(float(rt))
        except ValueError:
            continue
    return {k: statistics.median(v) for k, v in out.items()}


def ours(dataset: str) -> dict[str, tuple[float, float | None]]:
    """{key: (retention time, |model z|)} across both polarities."""
    out: dict[str, tuple[float, float | None]] = {}
    for path in sorted(STAGED.glob(f"{dataset}/Analysis_v3/Results/*/Final_Results.csv")):
        for row in csv.DictReader(path.open(errors="replace")):
            name = (row.get("Identification") or "").strip()
            if not name or name.upper().startswith("DECOY"):
                continue
            try:
                rt = float(row["Retention Time (min)"])
            except (KeyError, ValueError):
                continue
            z = row.get("RT Model Z") or ""
            try:
                zval = abs(float(z.replace("+", "")))
            except ValueError:
                zval = None
            key = _key(row.get("Shorthand (LSI)") or name)
            if key:
                out.setdefault(key, (rt, zval))
    return out


def main() -> None:
    dataset, arm = "MTBKS222", "Thermo"
    theirs, mine = published_rt(dataset, arm), ours(dataset)
    shared = sorted(set(theirs) & set(mine))

    print(f"  published {len(theirs)} distinct lipids with a retention time")
    print(f"  ours      {len(mine)}")
    print(f"  shared    {len(shared)}")
    if not shared:
        print("\n  nothing shared — check the shorthand mapping before reading anything into this")
        return

    deltas = [(k, mine[k][0] - theirs[k], mine[k][1]) for k in shared]
    absd = sorted(abs(d) for _, d, _ in deltas)
    agree = [d for d in deltas if abs(d[1]) <= AGREE_MIN]

    print(f"\n  |ΔRT| median {statistics.median(absd):.2f} min, "
          f"p90 {absd[int(.9 * len(absd)) - 1]:.2f} min, max {max(absd):.2f} min")
    print(f"  within {AGREE_MIN} min: {len(agree)} of {len(shared)} "
          f"({100 * len(agree) / len(shared):.0f}%) — the same peak, not just the same name")

    # The question the model exists to answer: does its own flag predict the disagreements?
    flagged = [(k, d, z) for k, d, z in deltas if z is not None and z >= 2.0]
    clean = [(k, d, z) for k, d, z in deltas if z is not None and z < 2.0]
    if flagged and clean:
        f_bad = sum(1 for _, d, _ in flagged if abs(d) > AGREE_MIN) / len(flagged)
        c_bad = sum(1 for _, d, _ in clean if abs(d) > AGREE_MIN) / len(clean)
        f_n = sum(1 for _, d, _ in flagged if abs(d) > AGREE_MIN)
        c_n = sum(1 for _, d, _ in clean if abs(d) > AGREE_MIN)
        print(f"\n  model flagged (|z| >= 2): {f_n} of {len(flagged)} disagree with the deposit")
        print(f"  model clean   (|z| <  2): {c_n} of {len(clean)} disagree")
        try:
            from scipy.stats import fisher_exact
            _, pval = fisher_exact([[f_n, len(flagged) - f_n], [c_n, len(clean) - c_n]])
            print(f"  Fisher exact p = {pval:.4f}")
            # ⚠ Report the counts, not only the ratio. This rests on 7 flagged lipids and the
            # effect is carried by 2 of them; an enrichment quoted alone would read as far
            # sturdier than it is.
        except ImportError:
            pass
        if c_bad > 0:
            print(f"  enrichment: {f_bad / c_bad:.1f}x — a flagged lipid is this much more likely "
                  f"to be placed elsewhere by the deposit")
        elif f_bad > 0:
            print("  enrichment: every disagreement is in the flagged set")
    else:
        print("\n  no usable z values — the model column is empty for the shared set")

    worst = sorted(deltas, key=lambda t: -abs(t[1]))[:10]
    print("\n  largest disagreements (ours minus theirs):")
    for k, d, z in worst:
        print(f"    {k:<28} {d:+7.2f} min   ours {mine[k][0]:6.2f}  theirs {theirs[k]:6.2f}"
              f"   z {z if z is not None else '—'}")


if __name__ == "__main__":
    main()
