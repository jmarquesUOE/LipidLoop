"""Generate an oxidised-triacylglycerol library, constrained by chemistry rather than combinatorics.

No vendor ships one: zero `OxTG` entries and zero TG names carrying `;O` across MS-DIAL's
1.06-million-spectrum atlas, both polarities. MS-DIAL generates those annotations from rules at
search time, so there is nothing to import and the class has to be generated — as the negative-mode
ganglioside library here was.

## What constrains the enumeration

**Autoxidation needs a bis-allylic CH2**, the carbon between two double bonds. A chain with
DB >= 2 oxidises readily; 18:1 is far slower and 16:0/18:0 do not go at all. Only parents carrying
an oxidisable chain are enumerated. That single constraint is the difference between ~500 entries
and ~60,000, and it is also the strongest falsifier available: 5 of the 22 OxTG annotations in one
public deposit carry no oxidisable chain anywhere in the molecule, and this generator cannot
reproduce them.

## Naming — and why NOT the LSI `;O` form

⚠ `TG 52:4;O` collapses onto the ordinary TG of the same C:DB. Verified against the live code:
`sum_composition("TG 50:3;O3")` returns `"TG 50:3"` and `parse_sum_name` judges it on the ordinary
TG retention surface — so such an entry would be merged by duplicate-name detection, split-peak
summing and adduct-pair removal, and would inherit a surface fitted to molecules 48 Da lighter.

The existing `[OH]`/`OH-` idiom, which LipiDex uses for hydroxylated phospholipids, survives all of
that: `TG[OH] OH-52:4` keeps its own class and its own retention surface.

## Sum composition, not molecular species

The oxidised chain's identity is well determined — its neutral loss differs from an ordinary fatty
acid by 36 mDa, resolvable at any sensible tolerance. The two unoxidised chains are determined only
as a PAIR, because the surviving diacyl ion gives their sum. Resolving `TG 18:1_18:2_18:2;O` from
`TG 18:0_18:3_18:2;O` needs both unoxidised acylium ions, and in the spectra examined that region
is at or below noise. A molecular-resolution library would give those two entries near-identical
spectra and the dot product could not separate them — the ganglioside argument, one level up.

## ⚠ Fragment INTENSITIES are provisional

Masses here are exact. Intensities are not measured, and the reason is worth recording: the one
public deposit annotating OxTGs is profile-mode data whose MS2, even after peak-picking, splits a
single fragment into several centroids —

    precursor 910.7470:  599.5012 (100%)  599.5082 (77%)  599.5151 (20%)  599.4752 (15%)

all one peak. Intensity medians from that are meaningless. The pattern used below asserts only what
the fragmentation demands — the acyl-loss ions exist and the oxidised-chain loss is prominent — and
`--intensities` exists so it can be re-parameterised from clean data without regenerating anything
else. **Until then, treat a match on this library as evidence of composition, not of confidence.**
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

H = 1.0078250319
C = 12.0
O = 15.9949146
NH3 = 17.0265491
PROTON = 1.00727646
OXYGEN = O


def fa_neutral(c: int, d: int, ox: int = 0) -> float:
    """Neutral fatty acid CnH(2n-2d)O2, plus `ox` extra oxygens."""
    return C * c + H * (2 * c - 2 * d) + O * 2 + ox * OXYGEN


def tg_neutral(chains, ox: int = 0) -> float:
    """Triacylglycerol from three fatty acids: glycerol + 3 esters, losing 3 water."""
    glycerol = C * 3 + H * 8 + O * 3
    return glycerol + sum(fa_neutral(c, d) for c, d in chains) - 3 * (H * 2 + O) + ox * OXYGEN


def parse_chains(name: str):
    return [(int(a), int(b)) for a, b in re.findall(r"(\d+):(\d+)", name)]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("out", type=Path)
    ap.add_argument("--parents", type=Path, required=True,
                    help="library whose TG molecular species define the parents to oxidise")
    ap.add_argument("--oxygens", default="1,2",
                    help="oxygen counts to generate; 1=hydroxy/keto/epoxide, 2=hydroperoxide")
    ap.add_argument("--min-db", type=int, default=2,
                    help="a chain needs this many double bonds to autoxidise (bis-allylic CH2)")
    ap.add_argument("--intensities", default="999,300,150",
                    help="oxidised-chain loss, ordinary-chain loss, [M+H]+ — PROVISIONAL")
    args = ap.parse_args()

    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from lipidloop.msp import parse_msp

    oxygens = [int(x) for x in args.oxygens.split(",")]
    i_ox, i_ord, i_mh = (int(x) for x in args.intensities.split(","))

    # Parents: TG molecular species carrying at least one oxidisable chain. Collapsed to the SUM
    # the oxidised form will be reported at, so each sum is generated once however many molecular
    # species share it.
    sums: dict[tuple[int, int], list] = {}
    for s in parse_msp(str(args.parents)):
        if not s.lipid.startswith("TG "):
            continue
        chains = parse_chains(s.lipid.split(" ", 1)[1])
        if len(chains) != 3 or max(d for _, d in chains) < args.min_db:
            continue
        key = (sum(c for c, _ in chains), sum(d for _, d in chains))
        sums.setdefault(key, chains)

    blocks = []
    for (tc, td), chains in sorted(sums.items()):
        for ox in oxygens:
            precursor = tg_neutral(chains, ox) + NH3 + PROTON
            peaks = []
            # the oxidised chain leaves as RCOOH + n oxygens, with ammonia
            oxidisable = [(c, d) for c, d in chains if d >= args.min_db]
            for c, d in dict.fromkeys(oxidisable):
                frag = precursor - (fa_neutral(c, d, ox) + NH3)
                # ⚠ The chain token must be a bare `c:d`. `[OH-18:2]` does not parse, and the
                # decoy builder keys on this annotation — an entry it cannot read gets no decoy,
                # so the class would ship with no FDR estimate at all. The oxidation is carried by
                # the fragment TYPE, which is where the other libraries put it.
                peaks.append((frag, i_ox, f"O-{ox}_Oxidised Alkyl Neutral Loss_[{c}:{d}]"))
            # an ordinary chain leaves instead, keeping the oxygen on the diacyl ion
            for c, d in dict.fromkeys(chains):
                frag = precursor - (fa_neutral(c, d, 0) + NH3)
                peaks.append((frag, i_ord, f"O-0_Alkyl Neutral Loss_[{c}:{d}]"))
            peaks.append((precursor - NH3, i_mh, "M-0_Precursor_[]"))
            peaks.sort()

            name = f"TG[OH] OH-{tc}:{td}" + (f";O{ox}" if ox > 1 else ";O")
            blocks.append(
                f"Name: {name} [M+NH4]+;\n"
                f"MW: {precursor:.4f}\n"
                f"PRECURSORMZ: {precursor:.4f}\n"
                # ⚠ `Type=LipiDex` is what makes the parser read the chain annotations at all;
                # without it every annotation is silently dropped and purity, chain evidence and
                # the decoy builder all see an unannotated spectrum.
                f"Comment: Name={name} [M+NH4]+ Mass={precursor:.3f} Type=LipiDex "
                f"generated=build_oxtg_positive_library.py intensities=PROVISIONAL\n"
                f"Num Peaks: {len(peaks)}\n"
                + "".join(f"{m:.4f}\t{i}\t\"{a}\"\n" for m, i, a in peaks))

    args.out.write_text("\n".join(blocks))
    print(f"  {len(blocks)} entries -> {args.out}")
    print(f"  {len(sums)} parent sum compositions carrying a chain with >= {args.min_db} DB")
    print(f"  oxygen counts: {oxygens}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
