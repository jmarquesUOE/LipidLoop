"""sn-position evidence from the relative intensity of acyl fragments already in the spectrum.

The pipeline resolves *which* chains a lipid carries — 43% of MS2-backed negative-mode
identifications on a reference study — and reports them as `PC 16:0_18:1`, where `_` says the
chains are known and their positions are not. The spectra that proved the chains usually also
carry information about their arrangement, and nothing currently reads it.

**The chemistry.** In negative-mode collisional dissociation of a diacyl glycerophospholipid, both
acyl chains appear as carboxylate anions, and the one esterified at **sn-2 is typically the more
intense**: the sn-2 ester is more labile, so its carboxylate forms preferentially. The effect is
reported consistently across PC, PE, PS, PG and PI.

**Why this is reported and not asserted.** It is a tendency, not a rule. The ratio moves with chain
length, unsaturation, collision energy and instrument, and it inverts for some species. So this
module returns the measured ratio and which chain it favours, and the caller writes that into its
own column — the identification NAME keeps `_`.

⚠ **Never emit `/` from this.** Under the LSI shorthand `PC 16:0/18:1` asserts sn-1 and sn-2 as
proven. An intensity ratio is probabilistic evidence, and w2's validation of the class mapping was
explicit that the pipeline must not upgrade a claim it cannot support. `_` stays; the evidence
sits beside it where a reader can weigh it.

⚠ **Restricted to classes where the rule is established.** Triacylglycerols are excluded: with
three chains and neutral-loss rather than carboxylate fragments, the intensity ordering reflects
which loss is favoured rather than a clean sn-2 preference, and sn-2 in a TG is the *least*
favoured neutral loss — the opposite convention. Applying the phospholipid rule there would invert
the answer.
"""
from __future__ import annotations

import re
from dataclasses import dataclass

# Classes where the sn-2 carboxylate is reported to dominate in negative mode. Diacyl
# glycerophospholipids and their lyso forms only.
SN2_DOMINANT = {
    "PC", "PE", "PS", "PG", "PI", "PA",
    "Plasmanyl-PC", "Plasmenyl-PC", "Plasmanyl-PE", "Plasmenyl-PE",
    "PC[OH]", "PE[OH]",
}

# Below this the two fragments are not meaningfully different and no call is made. Chosen high
# deliberately: a 1.2x difference is within the run-to-run variation of fragment intensity.
MIN_RATIO = 1.5

_CD = re.compile(r"(\d+):(\d+)")


@dataclass(slots=True)
class SnEvidence:
    """What the fragment intensities say about chain arrangement."""

    chains: dict[str, float]          # chain -> observed intensity of its diagnostic fragment
    favoured_sn2: str = ""            # the chain the ratio puts at sn-2, or "" if undecided
    ratio: float = 0.0                # intensity of the favoured chain over the next
    note: str = ""

    def __str__(self) -> str:
        # ⚠ Surface the REASON when no call is made. Returning "" for every declined case made
        # the column uniformly empty and gave no way to tell "class not covered" from "fragment
        # not observed" from "the code never ran" — which is exactly the debugging position this
        # is currently in.
        if self.favoured_sn2:
            return f"{self.favoured_sn2} favoured sn-2 ({self.ratio:.1f}x)"
        return self.note or "no call"


def _observed(mz: float, sample_mz, sample_intensity, tol: float) -> float:
    """Intensity of the sample peak nearest `mz`, or 0.0 when nothing is within tolerance."""
    best, best_gap = 0.0, tol
    for m, i in zip(sample_mz, sample_intensity):
        gap = abs(m - mz)
        if gap <= best_gap:
            best, best_gap = float(i), gap
    return best


def _is_positive(polarity: str) -> bool:
    """True for a stated positive polarity; False for negative or for nothing stated."""
    p = (polarity or "").strip().lower()
    return p.startswith("+") or p.startswith("pos")


def evidence(lipid_class: str, lipid_name: str, library_mz, library_chains,
             sample_mz, sample_intensity, tol: float = 0.01,
             polarity: str = "") -> SnEvidence:
    """Weigh the acyl-fragment intensities for one matched spectrum.

    `library_chains` is the per-fragment chain annotation (`LibrarySpectrum.fatty_acids`), which is
    what makes this cheap: the library already says which fragment belongs to which chain, so no
    fragmentation chemistry has to be re-derived here.

    `polarity` accepts `+`/`-` or `positive`/`negative`, from either the sample spectrum or the
    library entry's adduct. Callers on the production path pass it; it defaults to unstated, which
    is permissive, because a caller that cannot say should not be silently refused.
    """
    # ⚠ NEGATIVE MODE ONLY, and this guard was missing until real data exposed it.
    #
    # Everything above is carboxylate chemistry: in negative mode both chains leave as carboxylate
    # anions and the sn-2 ester, being the more labile, gives the more intense one. A protonated
    # phospholipid does not do this — it fragments to the head group, and what acyl-related ions
    # appear are neutral losses whose intensity ordering reflects which loss is favoured, not
    # ester lability. Scoring those with the sn-2 rule is the exact error the TG exclusion above
    # exists to prevent, applied to a whole polarity instead of one class.
    #
    # It was not hypothetical: the first wired run made three calls on `PE 18:0_20:4 [M+H]+`,
    # because the class matched, the library entry carried chain annotations, and nothing ever
    # asked which polarity the spectrum came from.
    if _is_positive(polarity):
        return SnEvidence({}, note="positive mode — the sn-2 rule is negative-mode chemistry")

    if lipid_class not in SN2_DOMINANT:
        return SnEvidence({}, note="class not covered — the sn-2 rule is not established here")

    chains = [f"{c}:{d}" for c, d in _CD.findall(lipid_name.partition(" ")[2])]
    if len(chains) != 2:
        # One chain is a lyso species with nothing to compare; three or more is not a class
        # this rule covers, and should have been excluded above.
        return SnEvidence({}, note="needs exactly two resolved chains")
    if chains[0] == chains[1]:
        return SnEvidence({}, note="identical chains — arrangement is unobservable")

    found: dict[str, float] = {}
    for mz, chain in zip(library_mz, library_chains):
        if not chain:
            continue
        key = chain.strip()
        if key not in chains:
            continue
        seen = _observed(float(mz), sample_mz, sample_intensity, tol)
        # A chain can have several diagnostic fragments; take the strongest, as the comparison is
        # between the best evidence for each chain rather than between arbitrary fragments.
        if seen > found.get(key, 0.0):
            found[key] = seen

    if len(found) < 2 or min(found.values()) <= 0:
        return SnEvidence(found, note="only one chain's fragment observed")

    ordered = sorted(found.items(), key=lambda kv: -kv[1])
    ratio = ordered[0][1] / ordered[1][1]
    if ratio < MIN_RATIO:
        return SnEvidence(found, ratio=ratio,
                          note=f"inconclusive ({ratio:.1f}x, under {MIN_RATIO}x)")
    return SnEvidence(found, favoured_sn2=ordered[0][0], ratio=ratio)


# The written form of a call, and its reader, kept together deliberately: the per-injection search
# table is the only carrier between the search stage and the delivered table, so a change to the
# format has to change the parser in the same edit or the delivered column silently goes blank.
_CELL = re.compile(r"^(\S+) favoured sn-2 \(([\d.]+)x\)$")


def parse_cell(cell: str) -> tuple[str, float]:
    """The chain and ratio out of a written cell; `("", 0.0)` for a declined one.

    Declines are stored as prose (`"class not covered — ..."`) rather than as an empty string, so
    this must not treat "unparseable" as an error — a decline is the common case and a legitimate
    one.
    """
    m = _CELL.match((cell or "").strip())
    return (m.group(1), float(m.group(2))) if m else ("", 0.0)
