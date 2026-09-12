"""Why is each published lipid we do not report missing?

An agreement percentage says how many were missed, not why, and the two have completely different
consequences. A class our libraries do not contain is a coverage decision we can reverse in an
afternoon; a species present in the library and not found is a detection failure; an annotation
that cannot exist is theirs, not ours. Reported together they average into a number that means
nothing.

Each miss is assigned the FIRST reason that applies, in this order — most decisive first, so a
chemically impossible entry is never also counted as a library gap:

  1. chemically impossible   the composition cannot exist for that class
  2. void                    annotated at or before the column's void time, where nothing is
                             resolved and the MS2 is chemical background
  3. class absent            we hold no library entry for the class at all
  4. species absent          the class is covered, this composition is not
  5. wrong polarity          covered, but only in the polarity this study did not run
  6. not detected            in the library, in the right polarity, and still not found
                             — the only category that is a failure of the software
"""
from __future__ import annotations

import csv
import json
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import deposits                                             # noqa: E402
from score_agreement import STAGED, canon, ours, published  # noqa: E402
from lipidloop.msp import parse_msp                        # noqa: E402
from lipidloop.nomenclature import canonical_class         # noqa: E402

REPO = Path(__file__).resolve().parents[2]
VOID_MIN = 1.0            # minutes; see the caveat printed at the end
# Minimum total carbons a class can physically have: one chain per acyl group, none below 12 C
# in practice, so anything under this cannot be the molecule named.
FLOOR = {"MG": 12, "DG": 24, "TG": 36, "PC": 24, "PE": 24, "PS": 24, "PG": 24, "PI": 24}


def dataset_libraries(dataset: str, analysis: str) -> list[str]:
    """Every `.msp` path actually configured for this dataset's resolved run, both polarities.

    ⚠ NOT `glob(data/libraries/*.msp)`. That counted `LipiDex_HCD_Plants.msp` and
    `OxTG_Positive.msp` as coverage this pipeline holds -- both real files on disk, neither wired
    into any run's config; they appear only as example text in `--extra-library`'s help. A species
    that exists only in one of those was misclassified `NOT DETECTED (in library, right
    polarity)` -- the one category that blames the software -- when the true reason is "never
    searched." Reading `run_config.json`'s own `libraries` list is what `ours()` already
    identified as this dataset's real run; crediting any other set of files as "held" compares
    the deposit against a pipeline that was never actually run.
    """
    paths: set[str] = set()
    for pol in ("Pos", "Neg"):
        cfg = STAGED / dataset / analysis / "Results" / pol / "run_config.json"
        if cfg.exists():
            paths.update(json.loads(cfg.read_text()).get("libraries", []))
    return sorted(paths)


def library_index(msp_paths):
    """(classes held, species held, species by polarity) from exactly these library files."""
    classes, species = set(), set()
    by_pol = defaultdict(set)
    for p in msp_paths:
        if "DECOY" in p or "HOLDOUT" in p:
            continue
        for s in parse_msp(p):
            key = canon(s.lipid)[1]
            if not key:
                continue
            species.add(key)
            # ⚠ Base class, brackets stripped. Our libraries are entirely subclass-tagged
            # (`HexCer[NS]`, `HexCer[AP]`) while deposits write plain `HexCer`, so comparing the
            # tagged form reported a class we hold thousands of entries for as absent.
            classes.add(canonical_class(key).split("[")[0])
            by_pol[s.polarity].add(key)
    return classes, species, by_pol


def impossible(key: str) -> bool:
    m = re.match(r"([A-Za-z0-9\[\]-]+) (\d+):(\d+)", key)
    if not m:
        return False
    cls = canonical_class(m.group(1)).split("[")[0]
    return cls in FLOOR and int(m.group(2)) < FLOOR[cls]


def deposit_rts(dataset: str) -> dict[str, float]:
    """Published retention times, where the deposit gives them."""
    out: dict[str, float] = {}
    for row in deposits.maf_rows(dataset):
        n = (row.get("metabolite_identification") or "").strip()
        rt = (row.get("retention_time") or "").strip()
        if n and rt:
            try:
                out[canon(n)[1]] = float(rt)
            except ValueError:
                pass
    extra = Path("/mnt/datastore/Jair/Projects/Non-Canonical-Ceramides_[ADS]_[NDS]/"
                 "OKC-HET_rat_plasma_multiomics/metadata")
    for f in extra.glob(f"{dataset}_*_annotations.csv"):
        for row in csv.DictReader(f.open(errors="replace")):
            n = (row.get("metabolite_name") or "").strip()
            rt = (row.get("Average Rt(min)") or "").strip()
            if n and rt:
                try:
                    out.setdefault(canon(n)[1], float(rt))
                except ValueError:
                    pass
    return out


def main() -> None:
    # Every deposit that publishes identifications, one row per genuinely independent run.
    # ⚠ Not bare "MTBKS222"/"ST000991" -- neither is a staged directory (only the per-arm folders
    # are), so `ours()` found nothing for either, `src` printed "none", and every published
    # species counted as "missing" against zero local results: not inflated, meaningless. MTBKS222
    # is six real vendor arms (see manuscript §2.3); `MTBKS222_Sciex` is excluded here as it is
    # everywhere else in this repo -- its 24 files ARE IMS_Normal + IMS_HighMass, staged again
    # under a third name (score_agreement.py's ARM comment), not a seventh run to report.
    for d in ("MTBKS222_Waters", "MTBKS222_Thermo", "MTBKS222_Agilent", "MTBKS222_IMS_Normal",
              "MTBKS222_IMS_HighMass", "MTBKS222_BrukerBAF", "MTBLS5163", "ST000991_DDA",
              "ST003514", "ST004797"):
        pub, (mine, src) = published(d), ours(d)
        ps = {s for _, s in map(canon, pub) if s}
        os_ = {s for _, s in map(canon, mine) if s}
        missing = ps - os_
        rts = deposit_rts(d)

        if src == "none":
            print(f"  {d}  —  no complete run found; skipped rather than compared against nothing\n")
            continue

        # ⚠ Per dataset, not once globally. Different studies load different library sets
        # (formate vs acetate, extra libraries) -- crediting one dataset's run with another's
        # coverage would misattribute a miss exactly the way the unwired-library bug did.
        lib_paths = dataset_libraries(d, src)
        classes, species, by_pol = library_index(lib_paths)
        polarities = {p.name for p in (STAGED / d / src / "Results").glob("*")
                      if p.is_dir() and p.name in ("Pos", "Neg")}
        run_pol = {"Pos": "positive", "Neg": "negative"}
        held = set().union(*(by_pol[run_pol[p]] for p in polarities)) if polarities else species

        reasons, examples = Counter(), {}
        for k in missing:
            cls = canonical_class(k).split("[")[0]
            if impossible(k):
                r = "chemically impossible"
            elif k in rts and rts[k] <= VOID_MIN:
                r = "void (<=1 min)"
            elif cls not in classes:
                r = "class absent from our libraries"
            elif k not in species:
                r = "species absent (class covered)"
            elif k not in held:
                r = "wrong polarity for this study"
            else:
                r = "NOT DETECTED (in library, right polarity)"
            reasons[r] += 1
            examples.setdefault(r, k)

        print(f"  {d}  —  {len(missing)} of {len(ps)} published species not reported "
              f"({src}, {'+'.join(sorted(polarities))}; libraries searched: {len(species):,} "
              f"species across {len(classes)} classes)")
        for r, n in reasons.most_common():
            print(f"     {n:>5}  {100*n/len(missing):>5.1f}%  {r:<44} e.g. {examples[r]}")
        print()

    print(f"  ⚠ the void cutoff is {VOID_MIN} min, taken from one deposit's gradient. It should be "
          f"derived per run\n     from the first confirmed lipid, not fixed — on a 6 min method it "
          f"is too generous, on a 25 min\n     method too strict. Deposits that publish no "
          f"retention times cannot be tested for it at all.")


if __name__ == "__main__":
    main()
