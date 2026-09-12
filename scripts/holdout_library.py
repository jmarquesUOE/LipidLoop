"""Remove named lipids from a library, so a search can be asked what it does without the answer.

A target-decoy measures whether a wrong answer beats a right one that is present. It is blind to
the case that has produced every real defect found so far: the right answer is **absent**, and the
software takes the best available wrong one rather than declining. That is how 19 acetate adducts
became phytoceramides.

This removes every entry whose lipid name is held out — all adducts, all charge states — so the
re-search has no way to be right about those peaks. What it reports there is a direct measure of
its willingness to assert without evidence.

    python scripts/holdout_library.py IN.msp OUT.msp --names held_out.txt
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path


_CD = re.compile(r"(\d+):(\d+)")


def lipid_of(name_line: str) -> str:
    """The lipid, stripped of adduct and trailing punctuation."""
    v = name_line.split(":", 1)[1].strip() if ":" in name_line else name_line.strip()
    v = re.sub(r"\s*\[[^\]]*\][+-]?\d*\s*;?\s*$", "", v)
    return v.strip().rstrip(";").strip()


def composition(lipid: str) -> str:
    """(class, total carbons, total double bonds) — the key both sides can be compared on.

    ⚠ Holding out by NAME removes almost nothing. The pipeline reports a sum composition
    (`PC 34:1`) whenever it cannot resolve chains, while the library holds molecular species
    (`PC 16:0_18:1`, `PC 15:0_19:1`, ...). Matching literally removed 66 entries for 89 held-out
    lipids — the search would still have found every one of them under a chain-resolved name, and
    the experiment would have measured nothing.

    Collapsing both sides to a composition removes every species that could produce the held-out
    name, which is what "the answer is absent" has to mean.
    """
    head = lipid.split(" ", 1)[0] if " " in lipid else lipid
    body = lipid[len(head):]
    ch = _CD.findall(body)
    if not ch:
        return lipid.upper()
    ether = "O" if re.search(r"\bO-", body) else ("P" if re.search(r"\bP-", body) else "")
    return (f"{head.upper()}{ether}|{sum(int(c) for c, _ in ch)}:"
            f"{sum(int(d) for _, d in ch)}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--names", type=Path, required=True,
                    help="one lipid name per line, without adduct")
    args = ap.parse_args()

    held = {composition(n.strip()) for n in args.names.read_text().splitlines() if n.strip()}
    kept = removed = 0
    entry: list[str] = []
    hit = False

    def flush(fh):
        nonlocal entry, hit, kept, removed
        if entry:
            if hit:
                removed += 1
            else:
                fh.write("\n".join(entry) + "\n\n")
                kept += 1
        entry, hit = [], False

    with args.out.open("w") as fh:
        for line in args.source.open(errors="replace"):
            line = line.rstrip("\n")
            if line.strip() == "":
                flush(fh)
                continue
            if line.lower().startswith("name:") and composition(lipid_of(line)) in held:
                hit = True
            entry.append(line)
        flush(fh)

    print(f"  {args.source.name}: kept {kept:,}, removed {removed:,} "
          f"({len(held)} names held out)")


if __name__ == "__main__":
    main()
