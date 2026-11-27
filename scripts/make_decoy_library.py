"""Build a target-decoy library: same class, same precursor, WRONG CHAINS.

A decoy is only useful if it competes. The plants decoy — galactolipids against human plasma —
scored 0 of 591, which bounds one kind of error and cannot touch the kind that actually bites:
naming the wrong SPECIES within a class that is genuinely present. That is what happened when 19
acetate adducts of ordinary ceramides were named as phytoceramides, and a plants library would
never have caught it.

So the decoy is built to be wrong in exactly that way:

  * **precursor m/z kept** — it competes for the same spectra, which is what makes the count a rate
  * **class and head-group fragments kept** — it still looks like that class, so the class-level
    evidence cannot separate it
  * **chain-specific fragments shifted by a whole number of CH2** — the decoy asserts acyl chains
    that cannot coexist with its own precursor mass

LipiDex's MSP annotates which fragments are chain-specific: `"O-1_Alkyl Fragment_[10:0]"` carries
the chain in brackets, while `"C7H14O2N1_Fragment_[]"` does not. That annotation is what makes this
possible without re-deriving the fragmentation chemistry.

A hit on such an entry is a false positive of the most informative kind: the software had the right
class and the wrong molecule, and could not tell.

⚠ **Classes whose fragments are ALL head-group get no decoy from this method**, because there is
no chain-specific peak to shift. Re-measured 2026-08-27 with the fixed `CHAIN` pattern, over the
positive-mode entries of the configured formate Pos pool: **13 classes, 12,300 of 112,497 entries,
10.9%**. (The figure this replaced — "21 classes, 32,315 entries, 18.3%" — was taken with the
broken pattern, which counted every sphingoid and ether entry as uncovered. It is void, and its
method could not be reproduced from anything in the repo, so the two numbers are not comparable
term for term; the method for this one is `_w1/coverage.py` in the rebuild scratch.)
SM is still in the uncovered set in positive mode — all 1,202 positive entries — so the head-group
construction is still needed. A class with no decoy cannot produce a decoy hit, so its
false-discovery rate is not low, it is **uncountable**, and the reported FDR silently assumes the
covered classes speak for it.

`--head-group-shift` builds a decoy for exactly those entries by a different rule: keep the
precursor, and displace every fragment by a mass that **no combination of atoms can produce**. That
gives up the "right class, wrong molecule" property — such a decoy is wrong about the class too —
but it is the only construction available when nothing in the spectrum names a chain, and an
uncountable class is worse than a bluntly-measured one.

⚠ Never apply it to entries that DO have chain-specific fragments; the CH2 decoy is the more
informative measurement there, and mixing the two would make the rate mean two different things at
once.

⚠ **A decoy can be a CORRECT answer, and then it is worse than useless.** `--screen-against`
rejects any decoy that a real library entry at the same precursor cannot be told apart from. See
`TargetIndex.collision` for the three rules and for the measured case that motivates them. Always
pass the whole target pool the study will search — a decoy built from one library can be
indistinguishable from a target in another, and the flagship case is exactly that.

    python scripts/make_decoy_library.py IN.msp OUT.msp --shift 3 \
        --screen-against data/libraries/LipidBlast_Formic.msp ...   # the WHOLE pool
    python scripts/make_decoy_library.py IN.msp OUT.msp --head-group-shift   # the uncovered classes
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

CH2 = 14.015650      # one methylene: the spacing of an acyl homologous series

# ⚠ The chain annotation is NOT always `[16:0]`. LipiDex writes sphingoid bases with a hydroxy-count
# prefix (`[d18:1]`, `[t18:0]`, `[m18:0]`), ether chains with a linkage prefix (`[O-16:0]` alkyl,
# `[P-18:1]` plasmenyl), and cardiolipin's DG fragments with TWO chains in one bracket
# (`[14:0, 16:1]`). The original pattern required a digit immediately after `[`, so it matched none
# of those. That was not a cosmetic miss — it broke the construction in two opposite directions at
# once:
#
#   * those entries were skipped by the CH2 construction, so plasmalogens, ether lipids, ceramides,
#     sphingomyelins and cardiolipins had no chain decoy at all — 14,061 of LipidBlast_Formic's
#     55,169 entries and 4,200 of LipiDex_HCD_Formic's 104,167;
#   * `head_group_decoy()` guards on this same pattern, so the guard never fired for them and they
#     WERE given a head-group decoy — which this module's own docstring forbids.
#     `DECOY2_HeadGroup_LipiDex_HCD_Formic` carried 4,200 such entries out of 14,205.
#
# `[]` (a peak with no chain) must still match nothing, which is what the `\d+:\d+` core enforces.
CHAIN = re.compile(r'_\[((?:[dtmDTM]|[OPop]-)?\d+:\d+(?:,\s*(?:[dtmDTM]|[OPop]-)?\d+:\d+)*)\]"?\s*$')

#: One chain inside that bracket. The prefix is captured so the substitution can carry it through:
#: ⚠ a shifted `[d18:1]` is `[d21:1]`, not `[21:1]`. Dropping the `d` silently relabels a sphingoid
#: base as an acyl chain, which is a different molecule and a different fragment formula.
_ONE_CHAIN = re.compile(r'([dtmDTM]|[OPop]-)?(\d+):(\d+)')


def shift_chains(chains: str, shift: int) -> tuple[str, int]:
    """`"14:0, 16:1"`, 3 -> `("17:0, 19:1", 2)`. Returns the new text and how many chains moved.

    Every chain named in one bracket moves, and the count comes back because the fragment's mass
    must move by `shift * CH2` PER CHAIN. Measured on the cardiolipin DG fragments, which are the
    only two-chain annotation in the shipped libraries: `[14:0, 14:0]` at 591.4031 against
    `[14:0, 16:1]` at 617.4188 is 26.0157 = 2*CH2 - H2, i.e. mass tracks each chain's carbons
    independently. Moving both chains but the mass by only one CH2 would write a fragment whose
    m/z and whose formula disagree.
    """
    moved = 0

    def one(m: "re.Match[str]") -> str:
        nonlocal moved
        moved += 1
        return f"{m.group(1) or ''}{int(m.group(2)) + shift}:{m.group(3)}"

    return _ONE_CHAIN.sub(one, chains), moved


def infer_chain_peaks(entries: "list[tuple[str, list[str]]]", tol: float = 0.01) -> dict[str, set]:
    """Which peaks are chain-specific, for libraries that do not annotate them.

    ⚠ Four libraries carry NO chain annotations at all — AcylHexCer (10,800 entries), AcylSM
    (10,416), UltraLongHexCer_EOS (3,759), UltraLongCer_EOS (3,179) — so `CHAIN` finds nothing and
    the CH2 construction cannot touch them. Together they are 28,283 entries, 13.5% of the formate
    positive pool, with no decoy of any kind. Their false-discovery rate is not low, it is
    uncountable.

    But they are not uninformative: every one of their peak lists is DISTINCT (10,800 of 10,800 for
    AcylHexCer), so a decoy built from them can genuinely compete. The annotation is missing, not
    the information — and the information can be recovered from the library itself.

    **The method.** Take two entries whose composition differs by exactly one CH2 in one chain. A
    peak carrying that chain must appear 14.0157 higher in the heavier entry; a head-group peak must
    appear at the same mass in both. Differencing homologues therefore reconstructs the annotation.

    Measured on AcylHexCer over 56 homologous pairs: 168 peaks move by exactly one CH2, 392 stay.

    Returns {entry name: set of chain-specific m/z}. An entry with no homologue in the library gets
    an empty set and is skipped by the caller rather than guessed at.
    """
    import re as _re
    num = _re.compile(r"(\d+):(\d+)")
    parsed = []
    for name, lines in entries:
        peaks = []
        for line in lines:
            t = line.strip()
            if t and t[0].isdigit():
                try:
                    peaks.append(float(t.split()[0]))
                except (ValueError, IndexError):
                    pass
        parsed.append((name, num.findall(name), peaks))

    out: dict[str, set] = {}
    by_len: dict[int, list] = {}
    for item in parsed:
        by_len.setdefault(len(item[1]), []).append(item)

    for name, comp, peaks in parsed:
        if not comp or not peaks:
            continue
        moved: set = set()
        for other, ocomp, opeaks in by_len.get(len(comp), ()):
            if other == name or not opeaks:
                continue
            diff = [(int(a[0]) - int(b[0]), int(a[1]) - int(b[1])) for a, b in zip(comp, ocomp)]
            # exactly one chain one carbon longer, everything else identical
            if sorted(diff) != sorted([(0, 0)] * (len(diff) - 1) + [(1, 0)]):
                continue
            for m in peaks:
                if any(abs(m - n - CH2) < tol for n in opeaks) and \
                   not any(abs(m - n) < tol for n in opeaks):
                    moved.add(round(m, 4))
            if moved:
                break
        if moved:
            out[name] = moved
    return out


def shift_entry(lines: list[str], shift: int) -> list[str] | None:
    """One MSP entry with its chain fragments moved `shift` methylenes. None if it has none."""
    out, moved = [], 0
    for line in lines:
        if line.lower().startswith("name:"):
            out.append(line.replace("Name:", "Name: DECOY_", 1))
            continue
        if line.startswith(("MW:", "PRECURSORMZ:", "Num Peaks:")):
            out.append(line)                     # precursor untouched — it must still compete
            continue
        if line.startswith("Comment:"):
            out.append(line.replace("Name=", "Name=DECOY_", 1))
            continue
        m = CHAIN.search(line)
        if not m:
            out.append(line)                     # head-group and class fragments stay
            continue
        parts = line.split(None, 2)
        if len(parts) < 2:
            out.append(line)
            continue
        try:
            mz = float(parts[0])
        except ValueError:
            out.append(line)
            continue
        note = parts[2] if len(parts) > 2 else ""
        chains, n_chains = shift_chains(m.group(1), shift)
        # A plain replacement string, not a template: chain text is digits and a letter, but
        # `re.sub` would still read a stray backslash or `\g` in it as a group reference.
        note = CHAIN.sub(lambda _m, c=chains: f'_[{c}]"', note)
        out.append(f"{mz + n_chains * shift * CH2:.4f} {parts[1]} {note}")
        moved += 1
    return out if moved else None


#: A displacement no formula can match: not a multiple of CH2 (14.01565), H2 (2.01565), O
#: (15.99491) or any common neutral loss, and non-integer so it cannot coincide with a nominal
#: mass difference. Large enough to leave the isotope envelope, small enough to stay in the
#: recorded range.
HEAD_GROUP_SHIFT = 6.3721


def head_group_decoy(lines: list[str], shift: float = HEAD_GROUP_SHIFT) -> "list[str] | None":
    """A decoy for an entry with no chain-specific fragment: precursor kept, every peak displaced.

    Returns None when the entry HAS chain annotations — those belong to the CH2 construction, which
    measures a more informative error.
    """
    if any(CHAIN.search(l) for l in lines):
        return None
    out, moved = [], 0
    for line in lines:
        if line.lower().startswith("name:"):
            out.append(line.replace("Name:", "Name: DECOY_", 1)); continue
        if line.startswith(("MW:", "PRECURSORMZ:", "Num Peaks:")):
            out.append(line); continue           # precursor untouched — it must still compete
        if line.startswith("Comment:"):
            out.append(line.replace("Name=", "Name=DECOY_", 1)); continue
        parts = line.split(None, 2)
        try:
            mz = float(parts[0])
        except (ValueError, IndexError):
            out.append(line); continue
        # ⚠ Displace UP, never down: a downward shift can walk a fragment into the noise region
        # below the recorded range, and a decoy peak that cannot be observed cannot compete.
        note = parts[2] if len(parts) > 2 else ""
        out.append(f"{mz + shift:.4f} {parts[1]} {note}".rstrip())
        moved += 1
    return out if moved else None


# ---------------------------------------------------------------------------------------------
# The collision screen
# ---------------------------------------------------------------------------------------------
#
# ⚠ This is the primary fix in this module, and it exists because a decoy can be a CORRECT answer.
#
# `DECOY_ PE-NMe2 15:0_17:1 [M+H]+` was built at `--shift 1` from `PE-NMe2 14:0_16:1`. PE-NMe2 is
# PE plus 2 CH2 in the head group; +1 CH2 on each of two chains puts those 2 CH2 back. The result
# has the same precursor (718.5387), the same molecular formula (C39H77N1O8P1) and the same
# acylium ions (239.2370, 265.2525) as a REAL `PE 16:0_18:1 [M+H]+` — so it is a true description
# of an abundant human lipid wearing the wrong class label. It scored as a false positive while
# being right (40.3% of every decoy hit in the v5 batch), and because ties go to the heavier
# entry it also OUTRANKED and deleted the true PE.
#
# The condition generalises:
#
#     mass(head group A) - mass(head group B) = shift * (number of distinct annotated chains)
#
# so no single `--shift` is safe for every class. At shift 1 with two chains you need a 2 CH2 gap
# (PE/PE-NMe2). At shift 3 with ONE chain you need a 3 CH2 gap — and PC is exactly PE + 3 CH2, so
# the lyso classes are exposed where the diacyl ones are not. Anticipating the next colliding pair
# is a losing game; screening the built decoy against the real library is not, because it catches
# the whole error class without anyone having to name it in advance.
#
# The test is deliberately asymmetric. A decoy is rejected when its chain-bearing fragments are a
# SUBSET of a real target's at the same precursor, not only when they are equal: a subset means
# every chain the decoy asserts is genuinely present in a real molecule that competes for the same
# spectra, so no evidence in the spectrum can tell them apart.

# The searcher's own windows — `search.py :: MS1_TOL / MS2_TOL`, which are LipiDex's GUI defaults.
# ⚠ Use these, not the 20 ppm identification window: the screen has to ask what the SEARCHER can
# separate, and the searcher bins on absolute Da.
MS1_TOL = 0.01
MS2_TOL = 0.01

_NAME_ADDUCT = re.compile(r"\[M[^\]]*\][+-]")


def parse_entry(lines: list[str]) -> tuple[str, float, list[float], list[float]]:
    """Name, precursor, chain-bearing peak m/z, all peak m/z — the fields the screen compares on."""
    name, precursor = "", 0.0
    chain_mz: list[float] = []
    all_mz: list[float] = []
    for line in lines:
        if line.lower().startswith("name:"):
            name = line.split(":", 1)[1].strip().rstrip(";").strip()
            continue
        if line.startswith("PRECURSORMZ:"):
            try:
                precursor = float(line.split(":", 1)[1])
            except ValueError:
                pass
            continue
        if line.startswith(("MW:", "Comment:", "Num Peaks:")):
            continue
        parts = line.split(None, 2)
        try:
            mz = float(parts[0])
        except (ValueError, IndexError):
            continue
        all_mz.append(mz)
        if CHAIN.search(line):
            chain_mz.append(mz)
    return name, precursor, chain_mz, all_mz


def lipid_class(name: str) -> str:
    """`DECOY_ PE-NMe2 15:0_17:1 [M+H]+` -> `PE-NMe2`. Prefix and adduct stripped."""
    n = name.replace("DECOY_", "", 1).strip()
    return n.split(" ", 1)[0] if " " in n else n


def name_composition(name: str) -> frozenset[str]:
    """The chains a NAME asserts: `PE 16:0_18:1 [M+H]+` -> {"16:0", "18:1"}.

    A set, not a list: `TG 16:0_16:0_18:1` and `TG 16:0_18:1` assert the same distinct chains, and
    the fragment annotations cannot tell them apart either.
    """
    core = _NAME_ADDUCT.sub("", name.replace("DECOY_", "", 1)).strip().rstrip(";").strip()
    if " " not in core:
        return frozenset()
    return frozenset(c.strip() for c in core.split(" ", 1)[1].split("_") if ":" in c)


def fragment_composition(lines: list[str]) -> frozenset[str]:
    """The chains an entry's FRAGMENTS assert, after shifting. `[14:0, 16:1]` contributes both."""
    chains: set[str] = set()
    for line in lines:
        m = CHAIN.search(line)
        if m:
            chains.update(c.strip() for c in m.group(1).split(","))
    return frozenset(chains)


def _contained(small: list[float], large: list[float]) -> bool:
    """Every m/z in `small` has a partner in `large` within MS2_TOL."""
    return all(any(abs(a - b) <= MS2_TOL for b in large) for a in small)


class TargetIndex:
    """Every real library entry, binned on precursor, ready to be asked "is this decoy you?".

    Binning is on `int(precursor * 100)` with the neighbouring bins scanned, which covers the
    0.01 Da MS1 window and keeps the lookup O(1) against a 200,000-entry pool.
    """

    def __init__(self) -> None:
        self.bins: dict[int, list[tuple[str, float, list[float], frozenset[str]]]] = {}
        self.n = 0

    def add_file(self, path: Path) -> int:
        added = 0
        for lines in read_entries(path):
            name, precursor, chain_mz, _all_mz = parse_entry(lines)
            if not precursor:
                continue
            self.bins.setdefault(int(precursor * 100), []).append(
                (name, precursor, chain_mz, name_composition(name)))
            added += 1
        self.n += added
        return added

    def _near(self, precursor: float):
        base = int(precursor * 100)
        for b in (base - 1, base, base + 1):
            for row in self.bins.get(b, ()):
                if abs(row[1] - precursor) <= MS1_TOL:
                    yield row

    def collision(self, precursor: float, frag_mz: list[float],
                  composition: frozenset[str] = frozenset()) -> tuple[str, str] | None:
        """`(colliding target name, which rule fired)`, or None.

        Three rules, and all three are needed — this is not belt-and-braces, it is measurement.
        On `DECOY_ PE-NMe2 15:0_17:1 [M+H]+`, the flagship defect, only two of the three fire:

        | rule | fires? | why |
        |---|---|---|
        | `decoy-in-target`  | **no**  | the decoy carries 4 chain peaks (acylium + water loss), LipidBlast's real `PE 16:0_18:1` lists only the 2 acylia |
        | `target-in-decoy`  | yes | those 2 acylia, 265.2524 / 239.2368, ARE 2 of the decoy's 4, to 0.0002 Da |
        | `same-composition` | yes | the shifted fragments assert {16:0, 18:1} and a real `PE 16:0_18:1` sits at the same 718.5387 |

        ⚠ So the containment test alone, in the direction it is natural to write it, would have
        passed the very decoy this screen exists to remove. A library entry is an editor's list of
        the peaks worth naming, not the spectrum — the real PE spectrum does contain the water
        losses, LipiDex just did not write them down. Containment must therefore be tested in BOTH
        directions, and the composition rule carries the argument that does not depend on how
        completely either library was annotated.
        """
        if not frag_mz:
            return None                    # ⚠ vacuously-true subset: never a collision
        for name, _mz, target_frags, target_comp in self._near(precursor):
            if composition and target_comp and composition == target_comp:
                return name, "same-composition"
            if not target_frags:
                continue
            if _contained(frag_mz, target_frags):
                return name, "decoy-in-target"
            if _contained(target_frags, frag_mz):
                return name, "target-in-decoy"
        return None


def read_entries(path: Path):
    """Yield each blank-line-delimited MSP entry as a list of lines."""
    entry: list[str] = []
    for line in path.open(errors="replace"):
        line = line.rstrip("\n")
        if line.strip() == "":
            if entry:
                yield entry
                entry = []
            continue
        entry.append(line)
    if entry:
        yield entry


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("source", type=Path)
    ap.add_argument("out", type=Path)
    ap.add_argument("--head-group-shift", action="store_true",
                    help="build decoys ONLY for entries with no chain-specific fragment, by "
                         "displacing every peak. For the classes the CH2 construction cannot "
                         "reach — SM, CerP, the sulfatides — which otherwise have no decoy at "
                         "all and so an uncountable rather than a low error rate.")
    ap.add_argument("--shift", type=int, default=3,
                    help="methylenes to move chain fragments by. Not 0, and large enough that the "
                         "shifted fragment cannot be the real one within tolerance. Default 3.")
    ap.add_argument("--screen-against", type=Path, nargs="*", default=[], metavar="TARGET.msp",
                    help="reject any decoy whose chain-bearing fragments are a subset of a real "
                         "entry's at the same precursor. ⚠ Pass EVERY target library the study "
                         "searches, not just this decoy's own source: the searcher scores one "
                         "spectrum against the whole pool at once, so a decoy built from "
                         "LipidBlast can still be indistinguishable from a LipiDex target.")
    ap.add_argument("--rejects", type=Path, default=None,
                    help="write the rejected decoys and the real entry each collided with to this "
                         "CSV. The screen is only trustworthy if what it removed can be read.")
    args = ap.parse_args()
    if args.shift == 0:
        ap.error("--shift 0 would copy the target library, not decoy it")

    index = TargetIndex()
    for target in args.screen_against:
        n = index.add_file(target)
        print(f"  indexed {n:,} target entries from {target.name}")

    rejects: list[tuple[str, str, str, str]] = []   # decoy, class, colliding target, rule
    by_class: dict[str, int] = {}
    by_rule: dict[str, int] = {}
    written = skipped = 0

    def emit(fh, entry: list[str]) -> None:
        nonlocal written, skipped
        made = (head_group_decoy(entry) if args.head_group_shift
                else shift_entry(entry, args.shift))
        if not made:
            skipped += 1
            return
        if index.n:
            name, precursor, chain_mz, all_mz = parse_entry(made)
            # ⚠ A head-group decoy has no chain-bearing peak by construction, so screening it on
            # chain fragments would test nothing. Screen it on its whole displaced spectrum
            # instead — that IS what it asserts.
            probe = chain_mz if chain_mz else (all_mz if args.head_group_shift else [])
            hit = index.collision(precursor, probe, fragment_composition(made))
            if hit:
                cls = lipid_class(name)
                rejects.append((name, cls, hit[0], hit[1]))
                by_class[cls] = by_class.get(cls, 0) + 1
                by_rule[hit[1]] = by_rule.get(hit[1], 0) + 1
                return
        fh.write("\n".join(made) + "\n\n")
        written += 1

    with args.out.open("w") as fh:
        for entry in read_entries(args.source):
            emit(fh, entry)

    print(f"  wrote {written:,} decoy entries to {args.out}")
    if args.head_group_shift:
        # ⚠ In this mode the skipped entries are the ones that DO carry chain fragments — the
        # opposite of the default. The default wording here would read as a failure when it is
        # the intended division of labour between the two constructions.
        print(f"  skipped {skipped:,} that carry chain-specific fragments — those belong to the "
              f"CH2 construction, which measures the more informative error")
    else:
        print(f"  skipped {skipped:,} with no chain-specific fragment — a decoy of those would be "
              f"identical to its target and would score as a target")
    if index.n:
        print(f"  rejected {len(rejects):,} by the collision screen "
              f"(indistinguishable from a real target at the same precursor)")
        for rule, n in sorted(by_rule.items(), key=lambda kv: -kv[1]):
            print(f"      rule {rule:18s} {n:6,}")
        for cls, n in sorted(by_class.items(), key=lambda kv: -kv[1]):
            print(f"      {cls:20s} {n:6,}")
    if args.rejects:
        with args.rejects.open("w") as fh:
            fh.write("decoy,class,collides_with,rule\n")
            for name, cls, hit, rule in rejects:
                fh.write(f'"{name}","{cls}","{hit}","{rule}"\n')


if __name__ == "__main__":
    main()
