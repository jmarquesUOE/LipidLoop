"""Assign a lipid class to the polarity that measures it better, and drop it from the other.

The two polarities are not two measurements of one thing. A choline lipid is a protonated ion in
positive and a formate adduct in negative; an acidic phospholipid is the reverse. Reporting both
invites a reader to compare numbers that are not on a common scale, and on this facility's first
paired study a third of the lipids found in both disagreed in sign.

**Only classes a run identifies in BOTH polarities are eligible.** A class seen in one polarity is
passed through untouched whatever the table says — on the study this was written for, that is 29 of
39 classes, 444 positive and 149 negative lipids. A rule applied to every class would have deleted
half the table, including `Cer[ADS]`, which is negative-only purely because the positive ceramide
library covers phyto bases and not dihydro ones.

Runs after both polarities, because deciding which classes are ambiguous needs both inventories,
and the pipeline processes one polarity at a time.

**Nothing is averaged or merged.** Different ionisation efficiencies mean the areas are not on a
common scale, and an average of two correct numbers is still meaningless.
"""
from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path
from statistics import mean, median

POLARITIES = ("Pos", "Neg")

RT_COLUMN = "Retention Time (min)"
#: How far apart the same molecule may elute in the two polarities before a cross-polarity chain
#: assignment is refused. Generous by design — it exists to catch a gross mismatch, not to do
#: chromatography. Measured cross-polarity offset on GM3 is +0.018 min (spread 0.010 across five
#: species, separate injections); median peak FWHM across these studies is 0.110 min, so this is
#: about 1.8 peak widths.
CHAIN_RT_TOLERANCE = 0.20
FILTERED = "Final_Results_Filtered.csv"
NORMALISED = "Final_Results_Filtered_median_Normalised.csv"
UNFILTERED = "Unfiltered_Results.csv"


# Identifications from a decoy library, which exist only to be counted as false positives.
DECOY_PREFIX = "DECOY"


class UndecidedClass(Exception):
    """A class appears in both polarities and nothing says which one to prefer."""


@dataclass
class Assignment:
    polarity: str
    basis: str = ""
    note: str = ""
    source: str = "default"


#: `Pos`/`Neg` are the canonical spellings; the file also contains `positive`/`negative` because
#: later rows were written by hand. Both are accepted — the vocabulary of the file is a fact to be
#: read, not a rule to be enforced at the cost of losing rules.
_SPELLINGS = {"pos": "Pos", "positive": "Pos", "+": "Pos",
              "neg": "Neg", "negative": "Neg", "-": "Neg"}


def _normalise_polarity(value: str | None) -> str:
    return _SPELLINGS.get((value or "").strip().lower(), (value or "").strip())


def load(default: str | Path, override: str | Path | None = None) -> dict[str, Assignment]:
    """Facility defaults, then the study's own file on top.

    A study that disagrees says so in its own file rather than editing the shared one, so the
    facility default stays a statement about chemistry and the study's exception stays visible as
    an exception.
    """
    table: dict[str, Assignment] = {}
    for path, source in ((default, "default"), (override, "study")):
        if not path or not Path(path).exists():
            continue
        with Path(path).open() as fh:
            rows = csv.DictReader(line for line in fh if not line.lstrip().startswith("#"))
            for row in rows:
                name = (row.get("class") or "").strip()
                polarity = _normalise_polarity(row.get("polarity"))
                if not name:
                    continue
                # ⚠ A row that cannot be read is RAISED, not skipped. This module exists to refuse
                # to guess a polarity; silently discarding a rule is the same failure wearing a
                # quieter face. Five rows — BA, Cer[EOS], HexCer[EOS], ASM, AHexCer — were written
                # `positive`/`negative` where the older rows say `Pos`/`Neg`, and were dropped for
                # it: 45 rules loaded from 50. Every one has a live library, so the first study to
                # identify one in both polarities would refuse with rc=3, exactly as Rat_Heart_Lumos did
                # on an undeclared ganglioside. The log said `polarity filter: 45 classes in the
                # table` every run and nobody had a reason to know 50 was the number.
                if polarity not in POLARITIES:
                    raise ValueError(
                        f"{path}: class {name!r} has polarity {row.get('polarity')!r}, which is "
                        f"not one of {sorted(POLARITIES)} nor a spelling of them. Fix the row — "
                        f"dropping it would silently disable the rule.")
                table[name] = Assignment(polarity=polarity,
                                         basis=(row.get("basis") or "").strip(),
                                         note=(row.get("note") or "").strip(),
                                         source=source)
    return table


def identified_classes(path: str | Path) -> dict[str, int]:
    """Class -> identified row count, from a result table."""
    counts: dict[str, int] = {}
    with Path(path).open(newline="") as fh:
        for row in csv.DictReader(fh):
            if not (row.get("Identification") or "").strip():
                continue
            name = (row.get("Lipid Class") or "").strip()
            if name:
                counts[name] = counts.get(name, 0) + 1
    return counts


@dataclass
class Validation:
    """Cases A-D from the requirement. Only A stops the run."""

    undecided: list = field(default_factory=list)   # A: in both polarities, not in the table
    stale: list = field(default_factory=list)       # B: in the table, in neither polarity
    moot: list = field(default_factory=list)        # C: in the table, now only one polarity
    single: list = field(default_factory=list)      # D: one polarity — the normal case
    shared: list = field(default_factory=list)      # in both AND decided

    @property
    def ok(self) -> bool:
        return not self.undecided


def validate(pos: dict[str, int], neg: dict[str, int],
             table: dict[str, Assignment]) -> Validation:
    """Compare the configured table against what this run actually identified.

    The inventory is a property of the LIBRARIES, not of chemistry: extend the positive ceramide
    library and `Cer[ADS]` becomes a both-polarities class and enters this decision for the first
    time. A hard-coded table would not notice, and the failure is silent and in the worst
    direction — it keeps producing output, and the output is missing lipids nobody asked it to drop.
    """
    # ⚠ A decoy is not a lipid and has no polarity to assign. Decoy entries carry a `DECOY_`
    # prefix, so their class parses as the pseudo-class `DECOY_` — which appears in both
    # polarities, is absent from every class-polarity table, and made this function refuse the
    # whole run. Correct behaviour on an unknown class, wrong input: a validation artefact should
    # never have reached a decision about real chemistry. It killed one dataset in a batch.
    seen_pos = {c for c in pos if not c.upper().startswith(DECOY_PREFIX)}
    seen_neg = {c for c in neg if not c.upper().startswith(DECOY_PREFIX)}
    both = seen_pos & seen_neg
    out = Validation()
    out.undecided = sorted(both - set(table))
    out.stale = sorted(set(table) - (seen_pos | seen_neg))
    out.moot = sorted(set(table) & (seen_pos ^ seen_neg))
    out.single = sorted(seen_pos ^ seen_neg)
    out.shared = sorted(both & set(table))
    return out


def check(validation: Validation) -> None:
    """Raise on case A. Guessing is how `Cer[ADS]` gets deleted."""
    if validation.undecided:
        raise UndecidedClass(
            f"{len(validation.undecided)} class(es) appear in both polarities but are not in the "
            f"class-polarity table: {', '.join(validation.undecided)}. Add them to the study's "
            f"class_polarity.csv, or to the facility default if the assignment is general. "
            f"Refusing to choose: a wrong guess deletes real lipids from one polarity and nothing "
            f"downstream would show it.")


# ── applying it ────────────────────────────────────────────────────────────────────────────────

def _numeric_columns(fieldnames) -> list[str]:
    from .peakfinder import META_COLUMNS
    return [c for c in fieldnames if c and c not in META_COLUMNS]


def _normalise(rows: list[dict], columns: list[str]) -> None:
    """Median-normalise in place, then re-centre on these rows.

    Mirrors `peakfinder.median_factors` and `peakfinder.recentre`; `test_polarity.py` pins the two
    against each other. It is done again here because the surviving row set is not known until the
    polarity filter has run, and a matrix normalised before that filter is no longer centred on the
    rows it ends up containing — which is the bug this pipeline has already had once.
    """
    def value(row, column):
        try:
            return float(row[column] or 0)
        except (TypeError, ValueError):
            return 0.0

    complete = [r for r in rows if all(value(r, c) > 0 for c in columns)]
    if complete:
        totals = [median([value(r, c) for r in complete]) for c in columns]
        centre = mean(totals) or 1.0
        factors = [t / centre if t else 1.0 for t in totals]
    else:
        factors = [1.0] * len(columns)

    for row in rows:
        for column, factor in zip(columns, factors):
            row[column] = value(row, column) / factor

    # re-centre on the rows written, over cells with signal
    medians = []
    for column in columns:
        present = [value(r, column) for r in rows if value(r, column) > 0]
        medians.append(median(present) if present else 0.0)
    centre = mean([m for m in medians if m > 0] or [1.0])
    for row in rows:
        for column, m in zip(columns, medians):
            if m > 0:
                row[column] = value(row, column) / (m / centre)


def _spearman(a, b) -> float:
    def rank(v):
        order = sorted(range(len(v)), key=lambda i: v[i])
        out = [0.0] * len(v)
        i = 0
        while i < len(order):
            j = i
            while j + 1 < len(order) and v[order[j + 1]] == v[order[i]]:
                j += 1
            share = (i + j) / 2 + 1
            for k in range(i, j + 1):
                out[order[k]] = share
            i = j + 1
        return out

    x, y = rank(list(a)), rank(list(b))
    n = len(x)
    if n < 3:
        return float("nan")
    mx, my = mean(x), mean(y)
    num = sum((p - mx) * (q - my) for p, q in zip(x, y))
    den = (sum((p - mx) ** 2 for p in x) * sum((q - my) ** 2 for q in y)) ** 0.5
    return num / den if den else float("nan")


def agreement(analysis: str | Path, metadata: dict[str, str], group_by: str = "group") -> dict:
    """How well the two polarities agree on the lipids they both identified.

    A run-level quality figure, computed BEFORE the filter, because afterwards the two tables share
    no classes and the comparison becomes impossible. A run where this collapses is saying
    something about the acquisition rather than about the biology.

    Fold change is taken on the normalised matrix, between the first two levels of the study
    grouping. With no grouping, or fewer than three shared lipids, it returns nothing rather than a
    number that cannot mean anything.
    """
    import math

    root = Path(analysis) / "Results"
    per: dict[str, dict[str, float]] = {}
    for polarity in POLARITIES:
        path = root / polarity / NORMALISED
        meta_path = metadata.get(polarity)
        if not path.exists() or not meta_path or not Path(meta_path).exists():
            return {}
        groups: dict[str, str] = {}
        with Path(meta_path).open(newline="") as fh:
            for row in csv.DictReader(fh):
                name = (row.get("sample") or "").strip()
                if name:
                    groups[name] = (row.get(group_by) or "").strip()
        levels = sorted({v for v in groups.values() if v})
        if len(levels) < 2:
            return {}
        with path.open(newline="") as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            return {}
        columns = [c for c in rows[0] if c in groups]
        changes: dict[str, list] = {}
        for row in rows:
            name = (row.get("Identification") or "").strip()
            if not name:
                continue
            sides = []
            for level in levels[:2]:
                vals = []
                for column in columns:
                    if groups[column] != level:
                        continue
                    try:
                        value = float(row[column] or 0)
                    except (TypeError, ValueError):
                        continue
                    if value > 0:
                        vals.append(math.log2(value))
                sides.append(mean(vals) if vals else None)
            if sides[0] is not None and sides[1] is not None:
                changes.setdefault(name, []).append(sides[1] - sides[0])
        # ⚠ One name, several chromatographic peaks — so take the MEDIAN across them, not
        # whichever row happened to be last in the file. `changes[name] = ...` silently overwrote,
        # which meant a duplicated lipid was compared using an arbitrary peak in each polarity, and
        # need not even have been the same species. `SM d41:2` carries three positive-mode rows at
        # +0.073, +0.141 and -1.752; the last one won, and the -1.75 was then compared against
        # negative mode and reported to the client as an OPPOSITE direction. About a third of the
        # flagged lipids are on more than one row, so the flag was part artefact.
        per[polarity] = {k: median(v) for k, v in changes.items()}

    shared = sorted(set(per["Pos"]) & set(per["Neg"]))
    if len(shared) < 3:
        return {}
    pos = [per["Pos"][k] for k in shared]
    neg = [per["Neg"][k] for k in shared]
    opposite = sum(1 for a, b in zip(pos, neg) if a * b < 0)
    return {"shared": len(shared), "rho": _spearman(pos, neg),
            "opposite": opposite,
            "opposite_fraction": opposite / len(shared),
            "levels": levels[:2],
            # per lipid, so the disagreement survives into the combined table instead of being
            # resolved out of sight by the filter
            "by_lipid": {k: ("agree" if per["Pos"][k] * per["Neg"][k] >= 0 else "OPPOSITE")
                         for k in shared}}


def apply(analysis: str | Path, table: dict[str, Assignment], log=None,
          metadata: "dict[str, str] | None" = None, group_by: str = "group") -> dict:
    """Drop each shared class from the polarity that measures it worse, and rewrite the tables.

    `Unfiltered_Results.csv` and `Final_Results.csv` are left alone: the first is the audit trail
    and the second keeps the full filtered set, which is how every other filter in this pipeline
    behaves. Only the two analysis-ready tables are cut.
    """
    say = log or (lambda _: None)
    root = Path(analysis) / "Results"
    counts = {p: identified_classes(root / p / FILTERED) for p in POLARITIES}
    report = validate(counts["Pos"], counts["Neg"], table)
    check(report)

    # ⚠ BEFORE anything is rewritten. Afterwards the two tables share no classes by construction,
    # so the comparison silently returns nothing — which is what it did when this was computed at
    # the end, next to a docstring saying it had to be computed at the start.
    shared_agreement = agreement(analysis, metadata or {}, group_by)

    # Carry chain detail across before anything is dropped. A choline lipid fragmenting in
    # positive mode puts almost everything into the 184.07 head group and gives no usable acyl
    # ions, so its chains are only readable in negative — where the polarity filter is about to
    # remove it. Measured on the study this was written for: 0 of 165 PC rows chain-resolved in
    # positive against 51% of PE in negative. Quantify on the better ion, name from the better
    # spectrum.
    carried = _carry_chain_names(root, report, table, log=say)

    removed: dict[str, dict[str, int]] = {p: {} for p in POLARITIES}
    for polarity in POLARITIES:
        other = [c for c in report.shared if table[c].polarity != polarity]
        drop = set(other)
        for name in (FILTERED, NORMALISED):
            path = root / polarity / name
            if not path.exists():
                continue
            with path.open(newline="") as fh:
                reader = csv.DictReader(fh)
                fields = [f for f in (reader.fieldnames or []) if f]
                rows = [r for r in reader]
            kept, cut = [], {}
            for row in rows:
                lipid_class = (row.get("Lipid Class") or "").strip()
                if lipid_class in drop:
                    cut[lipid_class] = cut.get(lipid_class, 0) + 1
                else:
                    kept.append(row)
            if name == NORMALISED and kept:
                _normalise(kept, _numeric_columns(fields))
            with path.open("w", newline="") as fh:
                writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
                writer.writeheader()
                writer.writerows(kept)
            if name == FILTERED:
                removed[polarity] = cut
        if removed[polarity]:
            listed = ", ".join(f"{k} {v}" for k, v in sorted(removed[polarity].items()))
            say(f"  {polarity}: removed {sum(removed[polarity].values())} rows — {listed}")

    _annotate_unfiltered(root, report, table, removed)

    say(f"polarity filter: {len(report.shared)} classes decided, "
        f"{len(report.single)} single-polarity classes passed through untouched")

    summary = {
        "chains_carried": carried,
        "shared_classes": {c: table[c].polarity for c in report.shared},
        "basis": {c: table[c].basis for c in report.shared},
        "removed": removed,
        "undecided": report.undecided, "stale": report.stale, "moot": report.moot,
        "single_polarity_classes": report.single,
        "kept_untouched": {p: sum(n for c, n in identified_classes(root / p / FILTERED).items())
                           for p in POLARITIES},
        "agreement": shared_agreement,
    }
    import json
    (Path(analysis) / "polarity_filter.json").write_text(json.dumps(summary, indent=2))
    return {"validation": report, "removed": removed, "table": table, "summary": summary}


def _carry_chain_names(root: Path, report: Validation, table: dict[str, Assignment],
                        log=None, rt_tol: float = CHAIN_RT_TOLERANCE) -> int:
    """Record the chain-resolved name the other polarity found, beside the row it belongs to.

    Same molecule, matched on lipid class and sum composition, so `PC 34:1` in positive gains
    `PC 16:0_18:1` when negative resolved it. Areas are untouched — quantification stays on the ion
    that measures it best.

    ⚠ The name is written to `Chains (other polarity)`, NOT over `Identification`. Overwriting made
    an inference from a different injection read as a measurement in this one, and nothing in the
    row said otherwise except a provenance column nobody joins on. The sum composition is what THIS
    polarity measured; the molecular name is evidence from elsewhere and is labelled as such.

    ★ THREE GATES, and the order matters:

      1. Same lipid class and same SUMMED composition. Necessary and never sufficient on its own.
      2. Retention time within `rt_tol`. CONFIRMS the pair; it must never select between candidates
         — see below.
      3. Exactly one resolved name for that composition. Two is a real isomer pair, and picking one
         would assert chains the data does not choose between.

    ⚠ WHY RT CONFIRMS BUT CANNOT SELECT. In reversed phase, +2 carbons and +1 double bond very
    nearly cancel, so a species and its (n-2C, -1DB) partner co-elute by chemistry rather than by
    accident. Measured on GM3 in Rat_Heart_Lumos: `GM3-NANA d42:2` at 9.439 min and `d40:1` at
    9.452 min, 0.013 min apart, and the per-class RT model puts d42:2 at Z = +0.00 — a perfect fit,
    not a bad annotation. Matching on nearest RT would have paired the positive `d18:1_22:0`
    (sums to 40:1) with `d42:2`. Composition decides; RT only vetoes.

    The tolerance is generous on purpose: the cross-polarity offset measured on those same GM3
    species is +0.018 min (Neg later, spread 0.010 across five species) while median peak FWHM
    across these studies is 0.110 min. `rt_tol` is there to catch a gross mismatch — a different
    peak entirely — not to do chromatography.
    """
    say = log or (lambda _: None)
    from .peaks import sum_composition

    def _rt(row) -> "float | None":
        try:
            return float(row.get(RT_COLUMN) or "")
        except ValueError:
            return None

    carried = declined_rt = 0
    for polarity in POLARITIES:
        if not (root / polarity / FILTERED).exists():
            continue
        # ⚠ TWO DIFFERENT SCOPES, and conflating them cost this function its same-polarity half.
        #
        # `winning` — classes the polarity table decided IN FAVOUR of this polarity — bounds what
        # may be read from the OTHER table, because only a class present in both was ever decided.
        #
        # Reading this polarity's own table is not bounded by it at all. A PC found only in
        # positive is not "shared", has no polarity rule and needs none, and can still carry a sum
        # composition on one row and its molecular form on another. Gating that on `report.shared`
        # meant the same-polarity case never ran for the classes it most applies to.
        winning = [c for c in report.shared if table[c].polarity == polarity]
        other = "Neg" if polarity == "Pos" else "Pos"
        # sum composition -> {name: [(retention time, source polarity)]}. Keep every RT: the same
        # molecular species can appear on more than one row, and a match may use any of them.
        #
        # ★ BOTH POLARITIES ARE SOURCES, this one included. A sum composition and its molecular
        # form can sit on two rows of the SAME table — the search resolved chains for one spectrum
        # of the peak and not another — and that is the same molecule reported twice, not two
        # findings. Measured across five studies: 70 such pairs, of which 17 co-elute (one peak,
        # named twice) and 53 sit at different retention times, which are real chromatographic
        # isomers and must stay apart. The RT gate is what separates the two cases, so the same
        # rule serves both and nothing extra is needed to tell them apart.
        resolved: dict[str, dict[str, list]] = {}
        for source in (other, polarity):
            table_path = root / source / FILTERED
            if not table_path.exists():
                continue
            restrict = winning if source == other else None
            with table_path.open(newline="") as fh:
                for row in csv.DictReader(fh):
                    name = (row.get("Identification") or "").strip()
                    if "_" not in name:
                        continue
                    if restrict is not None and (row.get("Lipid Class") or "").strip() not in restrict:
                        continue
                    (resolved.setdefault(sum_composition(name), {})
                             .setdefault(name, []).append((_rt(row), source)))
        if not resolved:
            continue

        path = root / polarity / FILTERED
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            fields = [f for f in (reader.fieldnames or []) if f]
            rows = list(reader)
        for column in ("Chains Resolved", "Chains From", "Chains RT Delta"):
            if column not in fields:
                fields.insert(fields.index("Identification") + 1, column)

        for row in rows:
            name = (row.get("Identification") or "").strip()
            if not name or "_" in name:
                continue          # already chain-resolved, or nothing identified
            options = resolved.get(name)
            if not options:
                continue
            here = _rt(row)
            best: dict[str, tuple] = {}
            if here is not None:
                # Gate 2. An RT that cannot be read on either side does not veto — absence of a
                # retention time is not evidence the peaks differ, and older tables have no column.
                for n, hits in options.items():
                    usable = [(r, src) for r, src in hits if r is not None]
                    if not usable:
                        best[n] = (None, hits[0][1])
                        continue
                    r, src = min(usable, key=lambda x: abs(x[0] - here))
                    if abs(r - here) <= rt_tol:
                        best[n] = (r, src)
                if not best:
                    declined_rt += 1
                    continue
            else:
                best = {n: hits[0] for n, hits in options.items()}
            if len(best) > 1:
                continue          # a genuine isomer pair — the data does not choose
            chosen, (there, source) = next(iter(best.items()))
            row["Chains Resolved"] = chosen
            row["Chains From"] = source
            if here is not None and there is not None:
                row["Chains RT Delta"] = f"{there - here:+.3f}"
            carried += 1

        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)
    if carried or declined_rt:
        say(f"  chain-resolved names: {carried} recorded in `Chains Resolved`"
            + (f", {declined_rt} declined on retention time" if declined_rt else "")
            + " (the row keeps the name its own polarity measured; `Chains From` says "
              "which table the chains came from)")
    return carried


def _annotate_unfiltered(root: Path, report: Validation, table: dict[str, Assignment],
                         removed: dict) -> None:
    """Record the polarity decision in the audit trail rather than only in the delivered tables.

    A reader must be able to see that `SM` left the negative table and why, without opening the
    configuration. Only rows that were otherwise kept are marked — a row already dropped for being
    a blank or an isotope keeps the reason it was actually dropped.
    """
    for polarity in POLARITIES:
        path = root / polarity / UNFILTERED
        if not path.exists():
            continue
        losing = {c for c in report.shared if table[c].polarity != polarity}
        if not losing:
            continue
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            fields = [f for f in (reader.fieldnames or []) if f]
            rows = list(reader)
        if "Filter Status" not in fields:
            continue
        for row in rows:
            name = (row.get("Lipid Class") or "").strip()
            if (name in losing and not (row.get("Filter Status") or "").strip()
                    and (row.get("Identification") or "").strip()):
                row["Filter Status"] = f"polarity: {name} assigned to {table[name].polarity}"
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)


COMBINED = "Combined_Filtered_Normalised.csv"

# What `Pooled QC CV (%)` means for a row, so a reader can filter without inventing thresholds.
# The bands are the report's: this facility's own 20% target, and the 30% in common use for
# untargeted profiling (Dunn et al. 2011).
def qc_flag(cv: str) -> str:
    try:
        value = float(cv)
    except (TypeError, ValueError):
        return ""
    if value <= 20:
        return "good"
    if value <= 30:
        return "acceptable"
    return "unreliable"


def combine(analysis: str | Path, agreement_by_lipid: "dict | None" = None, log=None) -> int:
    """One analysis-ready matrix from the two polarity tables.

    Every analysis of a two-polarity run otherwise begins by doing the same four things by hand:
    strip the polarity suffix so the sample columns line up, decide what to do about classes
    measured in both, concatenate, and remember not to re-normalise. The last is silent when
    forgotten. The pipeline knows all four, so it emits the table itself.

    ⚠ **The two blocks are NOT on a common scale and are deliberately not renormalised together.**
    Each polarity is already median-normalised within itself; ionisation efficiency differs, so on
    the study this was written for the sample medians were about 7.0e6 in positive and 1.5e6 in
    negative. Renormalising the concatenation would be dominated by whichever polarity contributed
    more rows and would fold two scales into one that means nothing. The `Scale` column says which
    block a row belongs to, because the consequence is not obvious: per-lipid tests are unaffected,
    since every row is compared against itself across samples, but anything comparing rows to each
    other — principal components, clustering, heat maps, class sums — must standardise per row
    first or the positive block wins on magnitude and row count alone.
    """
    say = log or (lambda _: None)
    root = Path(analysis) / "Results"
    from .peakfinder import META_COLUMNS

    blocks, samples = [], None
    for polarity in POLARITIES:
        path = root / polarity / NORMALISED
        if not path.exists():
            continue
        with path.open(newline="") as fh:
            reader = csv.DictReader(fh)
            fields = [f for f in (reader.fieldnames or []) if f]
            rows = list(reader)
        if not rows:
            continue
        injections = [c for c in fields if c not in META_COLUMNS]
        # `..._Pos` and `..._Neg` name the same vial. Harmonised so the blocks stack; a mismatch
        # here would silently produce a matrix with two disjoint sets of columns.
        harmonised = [_strip_polarity(c, polarity) for c in injections]
        if samples is None:
            samples = harmonised
        elif harmonised != samples:
            raise ValueError(
                f"sample columns differ between polarities after removing the suffix:\n"
                f"  {samples}\n  {harmonised}\nThe two runs did not measure the same vials.")
        blocks.append((polarity, fields, injections, rows))

    if not blocks:
        return 0
    annotation = [c for c in blocks[0][1] if c in META_COLUMNS]
    # `Mode`, not `Polarity`: the result tables already carry a `Polarity` column holding the
    # quant ion's sign, and two columns of the same name is not a naming quibble — a reader loading
    # this with pandas gets one silently renamed, and which of the two they end up with depends on
    # the library. It bit this check first: a count of `Polarity` returned +/- rather than Pos/Neg.
    header = ["Mode", "Scale"] + annotation + ["QC Flag"] + samples
    if agreement_by_lipid:
        header.insert(header.index("QC Flag") + 1, "Polarity Agreement")

    out = root / COMBINED
    written = 0
    with out.open("w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(header)
        for polarity, fields, injections, rows in blocks:
            for row in rows:
                line = [polarity, f"median-normalised within {polarity}"]
                line += [row.get(c, "") for c in annotation]
                line.append(qc_flag(row.get("Pooled QC CV (%)", "")))
                if agreement_by_lipid:
                    line.append(agreement_by_lipid.get(
                        (row.get("Identification") or "").strip(), ""))
                line += [row.get(c, "") for c in injections]
                writer.writerow(line)
                written += 1
    say(f"wrote {COMBINED}: {written} lipids over {len(blocks)} polarities, "
        f"{len(samples)} injections — the two blocks are NOT on a common scale")
    return written


def _strip_polarity(column: str, polarity: str) -> str:
    """Remove the polarity token wherever it sits in the name.

    It is not always a suffix. This facility's samples are `..._F_Pos` but its pools are
    `Pool_Pos_01` — the token is in the middle — so a suffix-only rule matched the samples, left
    the pools as `Pool_Pos_01` against `Pool_Neg_01`, and the two blocks then had disjoint columns
    for four injections. The guard caught it; this is the fix.

    Bounded by underscores or the ends of the string, so a sample legitimately called
    `Positive_control` keeps its name.

    ⚠ Case-insensitive. The facility writes `_Pos`, but public deposits write `_POS`, `_pos` and
    `_positive` — a case-sensitive match left those columns untouched, so the two polarities
    harmonised to two disjoint sets and `combine` raised "the two runs did not measure the same
    vials" on data where they plainly had. Matching the case of one lab's file-naming habit is not
    something to encode.
    """
    import re
    return re.sub(rf"(?:^|_){polarity}(?=_|$)", "", column, count=1, flags=re.IGNORECASE)
