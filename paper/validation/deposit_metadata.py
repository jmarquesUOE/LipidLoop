"""Write per-dataset sample roles and groups from whatever each deposit actually publishes.

Every study here has been processed as **one group of samples with no blanks**, because
`run_study` infers roles from filenames and nothing else was supplied. Two things follow, and both
are defects in our numbers rather than missing conveniences:

**Blanks counted as samples.** The presence filter requires a feature to recur across the injection
group. Where solvent blanks sit in that group, a real lipid — absent from a blank by definition —
fails the filter, and what survives is disproportionately what appears in *both*. On MTBKS222's
Waters arm, 9 of 17 runs are blanks and 83 of 2,168 features survived; its 86.2% agreement, the
highest in the whole table, was measured on background.

**Every condition in one group.** In a case/control design, demanding that a feature appear across
all injections penalises exactly the features that differ between groups — the biology.

Nothing here is guessed. Each deposit's own metadata is read, in whatever form it publishes:

| source | studies | gives |
|---|---|---|
| SDRF (`Source Name`, `Assay Name`) | MTBKS222 | blank/sample roles per run |
| mwTab `SUBJECT_SAMPLE_FACTORS` | Workbench `ST*` | factors and raw file names |
| MetaboLights `s_*.txt` | `MTBLS*` | sample sheet |
| Xcalibur sequence already on disk | in-house | roles and injection order |

⚠ **A file that cannot be matched is reported, never silently defaulted.** Defaulting to "sample"
is precisely the bug this exists to fix, and it would look identical to success.

    python manuscript/validation/deposit_metadata.py [--write]
"""
from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
# ⚠ src too: _role delegates to lipidloop.blanks so the metadata cannot drift from the
# classifier that reads it.
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))
import deposits  # noqa: E402

STAGED = Path("/home/jair/validation_staged")
LOCAL_SDRF = deposits.LOCAL_DEPOSITS

#: Roles INFERRED from the data where a deposit publishes no metadata at all.
#
# ⚠ These are a weaker claim than everything else in this file and are recorded as such. Every
# other role here is read from a deposit's own SDRF, mwTab or sample sheet; these are read from
# the run names plus the signal level, and they go into a separate `sequence_inferred.csv` with a
# provenance note rather than being blended into the deposit-sourced sequences.
#
# MTBLS5163: `BCO%20MeOH` and `BCO%20EXT%201-3` are a methanol blank and three extraction blanks —
# `%20` is an encoded space, which is also why `blank|blk` matching never sees them. Each carries
# roughly a third of the median identification count (78-97 against 290), which is what a blank
# carrying background looks like. The name and the data agree; the deposit says nothing.
INFERRED_BLANKS = {
    "MTBLS5163": ("BCO", "run names decode to 'BCO MeOH' and 'BCO EXT 1-3'; "
                         "0.27-0.33x the median identification count"),
}

#: Sample-type strings `lipidloop.blanks.infer_role` recognises.
BLANK, QC, STANDARD, UNKNOWN = "Blank", "QC", "Std Bracket", "Unknown"
# ⚠ Standards were missing from this list while `blanks.infer_role` has supported them all along,
# so a standards injection was written as Unknown and read as a sample — the same failure as
# blanks-as-samples, and it matters more than it looks: the quality report measures mass accuracy
# and recovery against the spiked mix, and a standards run sitting in the sample group corrupts
# the very panel it should be validating.
#
# `Std Bracket` is written rather than `Standard` because that is the Xcalibur string
# `blanks._STANDARD_TYPE` matches.


_STAMP = re.compile(r'startTimeStamp="([^"]+)"')


def acquisition_time(path: Path) -> str:
    """The instrument's own record of when this run started, from the mzML header.

    ⚠ This, not the filename and not the sequence row. Verified on MTBLS5163: ordering 50 files by
    timestamp reproduces the injection numbers embedded in their names exactly, 1 through 50,
    strictly increasing. A deposit that renamed its files, or shipped no sequence at all, still
    carries the truth in the header — and every study here except one does.

    Injection order is not a nicety. Without it the drift, run-order and confounding sections of
    the quality report cannot run at all, and they report "not testable" rather than announcing
    that the input was missing.
    """
    try:
        with Path(path).open(errors="replace") as fh:
            head = ""
            while len(head) < 1_000_000:
                chunk = fh.read(65536)
                if not chunk:
                    break
                head += chunk
                match = _STAMP.search(head)
                if match:
                    return match.group(1)
                if "<spectrumList" in head or "<spectrum " in head:
                    break
    except OSError:
        return ""
    return ""


def _role(*texts: str) -> str:
    r"""The role, decided by the PIPELINE's own classifier rather than a copy of it.

    ⚠ Two implementations of one rule will diverge, and this pair already had. My `\bqc\b` did not
    match `QC_stdmix` — `_` is a word character — while `blanks._QC_PATTERN` reads
    `(^|[_\-])qc([_\-]|\d|$)` and does. So this file wrote "Std Bracket" for a run the pipeline
    then read as a QC. Metadata that disagrees with its reader is worse than no metadata, because
    both look correct in isolation and only the combination is wrong.

    `infer_role` takes a filename and a declared type and returns one of blank/qc/standard/sample.
    Its answer is mapped back to the Xcalibur strings it itself recognises, so what is written here
    round-trips to exactly the same role.
    """
    from lipidloop.blanks import infer_role

    blob = " ".join(t or "" for t in texts if t)
    role = infer_role(blob, blob)
    return {"blank": BLANK, "qc": QC, "standard": STANDARD}.get(role, UNKNOWN)


def from_mwtab(dataset: str) -> dict[str, tuple[str, str]]:
    """{raw file stem: (role, group)} from Workbench's SUBJECT_SAMPLE_FACTORS block."""
    path = deposits.CACHE / f"{dataset}.mwtab.txt"
    if not path.exists():
        return {}
    out: dict[str, tuple[str, str]] = {}
    for line in path.read_text(errors="replace").splitlines():
        if not line.startswith("SUBJECT_SAMPLE_FACTORS"):
            continue
        parts = [p.strip() for p in line.split("\t")]
        if len(parts) < 4:
            continue
        sample, factors = parts[2], parts[3]
        extra = parts[4] if len(parts) > 4 else ""
        # `Name:Value | Name:Value` -> a single group label, joined so a two-factor design does
        # not collapse to one of its factors.
        group = "_".join(v.split(":", 1)[1].strip() for v in factors.split("|")
                         if ":" in v) or "all"
        group = re.sub(r"\s+", "", group)[:40]
        role = _role(sample, factors, extra)
        for m in re.finditer(r"RAW_FILE_NAME[^=]*=([^;]+)", extra):
            stem = Path(m.group(1).strip()).stem
            if stem:
                out[stem] = (role, group)
        if not re.search(r"RAW_FILE_NAME", extra):
            out[sample] = (role, group)
    return out


def from_metabolights(dataset: str) -> dict[str, tuple[str, str]]:
    """{raw file stem: (role, group)} from a MetaboLights s_*.txt sample sheet."""
    out: dict[str, tuple[str, str]] = {}
    for path in sorted(deposits.CACHE.glob(f"{dataset}*.sample.tsv")):
        with path.open(errors="replace") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                name = (row.get("Sample Name") or row.get("Source Name") or "").strip()
                if not name:
                    continue
                factors = " ".join(v for k, v in row.items()
                                   if k and k.startswith("Factor Value") and v)
                out[name] = (_role(name, factors),
                             re.sub(r"\s+", "", factors)[:40] or "all")
    return out


def from_sdrf(dataset: str) -> dict[str, tuple[str, str]]:
    """{raw file stem: (role, group)} from a MetaboLights SDRF sitting with the download."""
    path = LOCAL_SDRF / dataset / f"{dataset}.sdrf.txt"
    if not path.exists():
        return {}
    out: dict[str, tuple[str, str]] = {}
    with path.open(errors="replace") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            assay = (row.get("Assay Name") or "").strip()
            source = (row.get("Source Name") or "").strip()
            name = (row.get("Characteristics[sample_name]") or "").strip()
            if not assay:
                continue
            role = _role(source, name)
            group = re.sub(r"\s+", "", source)[:40] or "all"
            stem = re.sub(r"_(Pos|Neg)$", "", assay)
            for key in {assay, stem, stem.split("_", 1)[-1]}:
                out[key] = (role, group)
    return out


def _key(name: str) -> str:
    """A filename reduced to what is comparable across a deposit and a staged copy.

    ⚠ Separators differ and nothing else does. ST004797 declares `Sample.1.mzML` and the staged
    file is `Sample 1.mzML`; matching literally left all 224 files unmatched and the whole study
    without groups. Dots, spaces, underscores and hyphens all collapse.
    """
    return re.sub(r"[.\s_-]+", "", name).lower()


def resolve(stem: str, table: dict[str, tuple[str, str]]) -> tuple[str, str] | None:
    """Match a staged file to a deposit record: exact, then separator-insensitive, then containment."""
    if stem in table:
        return table[stem]
    norm = {_key(k): v for k, v in table.items()}
    k = _key(stem)
    if k in norm:
        return norm[k]
    for key, value in norm.items():
        if key and (key in k or k in key):
            return value
    return None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--write", action="store_true", help="write the files; otherwise report only")
    args = ap.parse_args()

    studies = sorted(p for p in STAGED.iterdir() if p.is_dir()
                     and any((p / q).is_dir() for q in ("Pos", "Neg")))
    grand = Counter()
    for study in studies:
        name = study.name
        deposit = "MTBKS222" if name.startswith("MTBKS222") else name
        deposits.ensure(deposit)
        table = from_sdrf(deposit) or from_mwtab(deposit) or from_metabolights(deposit)
        files = [f for q in ("Pos", "Neg") for f in sorted((study / q).glob("*.mzML"))]
        if not files:
            continue
        if not table:
            marker = INFERRED_BLANKS.get(deposit)
            if not marker:
                print(f"  {name:<24} — no deposit metadata published")
                grand["no source"] += 1
                continue
            prefix, why = marker
            rows = [{"File Name": f.stem,
                     "Sample Type": BLANK if f.stem.startswith(prefix) else UNKNOWN}
                    for f in files]
            n = sum(1 for r in rows if r["Sample Type"] == BLANK)
            print(f"  {name:<24} {len(files):>3} files  {{'Blank': {n}, 'Unknown': {len(rows)-n}}}"
                  f"  ⚠ INFERRED")
            grand["inferred blank"] += n
            if args.write:
                out = study / "sequence_inferred.csv"
                with out.open("w", newline="") as fh:
                    fh.write(f"# INFERRED, not published by the deposit: {why}\n")
                    fh.write("Bracket Type=4\n")
                    w = csv.DictWriter(fh, fieldnames=["File Name", "Sample Type"])
                    w.writeheader(); w.writerows(rows)
                # `find_sequence` matches *equence*.csv, so this is picked up like any other —
                # the filename and the header comment are what mark it as inferred.
                (study / "sequence.csv").write_text(out.read_text())
            continue

        # Acquisition order first, so the sequence is written in the order the instrument ran.
        stamps = {f: acquisition_time(f) for f in files}
        dated = sorted((v, f) for f, v in stamps.items() if v)
        order = {f: i for i, (_, f) in enumerate(dated, start=1)}

        rows, groups, seen, unmatched = [], {}, Counter(), []
        for f in files:
            hit = resolve(f.stem, table)
            if hit is None:
                unmatched.append(f.stem)
                seen["unmatched"] += 1
                continue
            role, group = hit
            rows.append({"File Name": f.stem, "Sample Type": role,
                         "Injection Order": order.get(f, ""),
                         "Acquisition Time": stamps.get(f, "")})
            groups[f.stem] = group
            seen[role] += 1

        note = f"{dict(seen)}"
        dated_n = sum(1 for f in files if stamps.get(f))
        print(f"  {name:<24} {len(files):>3} files  {note}"
              f"  order {dated_n}/{len(files)}")
        if unmatched:
            # ⚠ Never default these to "sample". That is the bug being fixed, and it would look
            # exactly like success.
            print(f"      ⚠ {len(unmatched)} unmatched, left out: {', '.join(unmatched[:4])}"
                  f"{' …' if len(unmatched) > 4 else ''}")
        grand.update(seen)

        if args.write and rows:
            with (study / "sequence.csv").open("w", newline="") as fh:
                fh.write("Bracket Type=4\n")
                w = csv.DictWriter(fh, fieldnames=["File Name", "Sample Type",
                                                   "Injection Order", "Acquisition Time"])
                w.writeheader()
                # Written in acquisition order; a sequence that lists runs in a different order
                # than they happened is worse than none.
                w.writerows(sorted(rows, key=lambda r: (r["Injection Order"] == "",
                                                        r["Injection Order"])))
            for polarity in ("Pos", "Neg"):
                pol_files = sorted((study / polarity).glob("*.mzML"))
                if not pol_files:
                    continue
                with (study / f"metadata_{polarity}.csv").open("w", newline="") as fh:
                    w = csv.writer(fh); w.writerow(["sample", "group"])
                    for f in pol_files:
                        if f.stem in groups:
                            w.writerow([f.stem, groups[f.stem]])

    print(f"\n  totals: {dict(grand)}")
    print("  " + ("written" if args.write else "dry run — pass --write to write the files"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
