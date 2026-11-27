"""`Associated_Spectra.csv` — every spectral match, in one table, keyed to its compound group.

LipiDex writes one `_Results.csv` per injection and, separately, an `Associated_Spectra.csv`
holding the retained matches across the whole run. The per-injection files are the interface the
peak finder reads; this is the one a person reads. With 53 injections per polarity the per-file
layout means the question "what evidence named this row?" cannot be asked without opening 53
files and joining them by hand.

Two things here that LipiDex's own file does not have:

**`Compound Group`.** LipiDex's table carries no key, so the only way back to `Final_Results.csv`
is to match on name plus retention time — which is exactly the join that is ambiguous when one
molecule appears on several rows, and useless when two rows at one retention time carry different
names. The id written here is the group's position in the unfiltered set, the same id
`Final_Results.csv` and `Unfiltered_Results.csv` carry, so the join is exact.

**Every rank, not just the winner.** A rank-1 match with a rank-2 close behind it at a similar dot
product is a different kind of claim from a rank-1 standing alone, and the table should be able to
show the difference. `Rank` makes dropping the rest a one-column filter for anyone who wants
LipiDex's view.

`GaussianScore` is deliberately absent. The field exists on `Lipid` because Compound Discoverer
supplies it, and this pipeline has no equivalent — a column that is always zero reads as a
measurement rather than an absence.
"""
from __future__ import annotations

import csv
import re
from collections import defaultdict
from pathlib import Path

from .peakfinder import PeakFinderResult

COLUMNS = ["Compound Group", "Name", "Adduct", "Sample", "MS2 ID", "Rank",
           "Precursor", "Library Mass", "Delta m/z (ppm)", "Retention", "Polarity",
           "Dot Product", "Reverse Dot Product", "Purity", "Library", "Associated"]


def _links(result: PeakFinderResult) -> dict[tuple[str, int, int], list[int]]:
    """(sample, MS2 id, rank) -> the compound groups that identification ended up in.

    A list rather than a single id: a spectrum is attached per feature, and nothing stops one
    match being accepted by features in two groups. Collapsing that to one id would hide it.
    """
    ids = {id(g): n for n, g in enumerate(result.compound_groups, start=1)}
    out: dict[tuple[str, int, int], list[int]] = defaultdict(list)
    for group in result.compound_groups:
        gid = ids[id(group)]
        for candidate in group.lipid_candidates:
            for lipid in candidate.identifications:
                key = (lipid.sample.file, lipid.ms2_id, lipid.rank)
                if gid not in out[key]:
                    out[key].append(gid)
    return out


def write_associated_spectra(result: PeakFinderResult,
                             result_files: dict[str, Path],
                             path: str | Path) -> int:
    """Write one long-format table over every per-injection search result. Returns the row count.

    Reads the per-injection files back rather than taking the rows in memory, so the table is a
    faithful account of what the peak finder was given — if the two ever disagree, the bug is
    visible here instead of being papered over.
    """
    links = _links(result)
    written = 0
    with Path(path).open("w", newline="") as out:
        writer = csv.DictWriter(out, fieldnames=COLUMNS)
        writer.writeheader()
        for sample, source in result_files.items():
            with Path(source).open(newline="") as fh:
                for row in csv.DictReader(fh):
                    try:
                        ms2_id, rank = int(row["MS2 ID"]), int(row["Rank"])
                    except (ValueError, KeyError):
                        continue
                    name = row["Identification"].strip().rstrip(";").strip()
                    groups = links.get((sample, ms2_id, rank), [])
                    writer.writerow({
                        "Compound Group": "; ".join(str(g) for g in groups),
                        "Name": name,
                        "Adduct": _adduct(name),
                        "Sample": sample,
                        "MS2 ID": ms2_id,
                        "Rank": rank,
                        "Precursor": row["Precursor Mass"],
                        "Library Mass": row["Library Mass"],
                        "Delta m/z (ppm)": row["Delta m/z"],
                        "Retention": row["Retention Time (min)"],
                        "Polarity": "+" if "]+" in name else "-",
                        "Dot Product": row["Dot Product"],
                        "Reverse Dot Product": row["Reverse Dot Product"],
                        "Purity": row["Purity"],
                        "Library": row["Library"],
                        "Associated": "Associated" if groups else "",
                    })
                    written += 1
    return written


def _adduct(name: str) -> str:
    match = re.search(r"(\[M[^\]]*\][+-]\d?)", name)
    return match.group(1) if match else ""
