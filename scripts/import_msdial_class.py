"""Extract one MS-DIAL COMPOUNDCLASS into a library this pipeline can search.

MS-DIAL's atlas is MSP and our parser already reads that dialect — `NAME:` plus `PRECURSORTYPE:`
are recombined into the v1 shape internally — so this is extraction and relabelling, not a format
conversion.

⚠ **One relabelling, and it matters.** MS-DIAL files bile acids under the name prefix `ST`, the
same prefix it uses for sterols, and distinguishes them only in `COMPOUNDCLASS`. Our pipeline takes
the class from the first token of the name, so importing them verbatim would put bile acids and
sterols on ONE retention surface. They elute nothing like each other, so the fit would be junk and
would take the sterols down with it. They are written out as `BA`, which is what MS-DIAL's own
result tables call them.

Licence: the MS-DIAL libraries are CC-BY 4.0 (Tsugawa et al.). Redistribution and modification are
permitted with attribution and a statement of changes — this script IS the statement of changes,
and PROVENANCE.md records the attribution.
"""
from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path


def _chains(name: str) -> str:
    """Rewrite MS-DIAL chain notation into this pipeline's: `_` between chains, no parentheses."""
    # ⚠ The parenthesised chain appears at BOTH ends: leading for AHexCer `(O-14:0)16:1…`,
    # trailing for ASM `…30:1;O2(FA 14:0)`. Substituting a trailing separator only works for the
    # first and silently fuses the second into `;O214:0`, which parses as no chain at all.
    # Wrap it in separators and collapse afterwards, which is position-independent.
    out = re.sub(r"\((?:O-|FA )?(\d+:\d+)\)", r"_\1_", name)
    out = out.replace("/", "_")
    out = re.sub(r"_+", "_", out).replace(" _", " ")
    return re.sub(r"_+(\s|$)", r"\1", out)


def blocks(path: Path, compound_class: str):
    text = path.read_text(errors="replace")
    want = re.compile(rf"^COMPOUNDCLASS: ({compound_class})\s*$", re.M)
    for block in re.split(r"\n(?=NAME:)", text):
        if want.search(block):
            yield block


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("source", type=Path, help="MS-DIAL atlas .msp")
    ap.add_argument("out", type=Path)
    ap.add_argument("--class", dest="klass", required=True,
                    help="COMPOUNDCLASS regex, e.g. 'BileAcid|BASulfate'")
    ap.add_argument("--rename", default="",
                    help="rewrite the name's leading token, e.g. ST=BA")
    ap.add_argument("--min-peaks", type=int, default=3,
                    help="drop entries thinner than this; a 2-peak spectrum whose second peak is "
                         "the precursor cannot discriminate between stereoisomers")
    args = ap.parse_args()

    old, _, new = args.rename.partition("=")
    kept, thin = [], 0
    for block in blocks(args.source, args.klass):
        n = int(re.search(r"^Num Peaks: (\d+)", block, re.M).group(1))
        if n < args.min_peaks:
            thin += 1
            continue
        if old and new:
            block = re.sub(rf"^NAME: {re.escape(old)} ", f"NAME: {new} ", block, count=1, flags=re.M)
        # ⚠ Three dialect differences, each of which silently breaks the sum composition.
        #
        # `;2O` vs `;O2` — the oxygen count, written the other way round.
        # `/` vs `_`     — MS-DIAL separates chains with a slash, which here ASSERTS sn-position
        #                  and, worse, is not read as a chain separator at all: `Cer[EOS]
        #                  14:1;O2/26:1;O2` summed to `Cer[EOS] 14:1` instead of 40:2, silently
        #                  dropping the second chain and with it the ultra-long tail that is the
        #                  entire reason for importing this class.
        # `(O-14:0)` / `(FA 14:0)` — the esterified chain in parentheses, which no parser here
        #                  reads, leaving AHexCer with no parseable sum and no retention class.
        block = re.sub(r";(\d)O\b", r";O\1", block)
        block = re.sub(r"^(NAME: .*)$", lambda m: _chains(m.group(1)), block, flags=re.M)
        kept.append(block.strip())

    args.out.write_text("\n\n".join(kept) + "\n")
    print(f"  {len(kept)} entries written to {args.out}")
    print(f"  {thin} dropped for fewer than {args.min_peaks} peaks")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
