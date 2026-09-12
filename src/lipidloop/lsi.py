"""Translate the pipeline's lipid names into LSI shorthand (Liebisch et al., J Lipid Res 2020).

Why a column and not a rename: the pipeline's own names carry information the shorthand does not —
`Cer[NS]` names the LipiDex subclass, `Plasmanyl-PC` says which ether the library matched — and
throwing that away to gain interoperability would be a bad trade. Both are written, side by side.

Why it is needed at all: scoring recall against a published list is a name-matching problem, and
that is where the answer is decided. Against the NIST SRM 1950 list, a naive comparison reported
**SM 0/51 and Cer 0/29** on data where 18 SM and 5 Cer had plainly been identified — because their
notation is `SM 18:1;O2/16:0` and ours is `SM d34:1`. An LSI column makes that comparison
mechanical instead of a per-study translation exercise.

## The three rules this must not break

**`/` is a claim, `_` is not.** In LSI, `PC 16:0/18:1` asserts sn-1 and sn-2; `PC 16:0_18:1` says
the chains are known and their positions are not. The pipeline proves chains from fragments, never
positions, so it emits `_` — always. Writing `/` would assert something no MS2 here can support.

**`P-` is a claim too.** A plasmalogen has a vinyl ether bond. Where the pipeline demoted a
`Plasmenyl-` identification to `O-` because the evidence did not separate alkyl from alkenyl, this
must NOT put the `P-` back. w2's validation of the class mapping is explicit about that, and an
`RT model` row — which has no MS2 at all — can never be `P-`.

**Cardiolipin is emitted at species level.** `CL 72:7`, not a chain list. The written order of a CL
name encodes sn-pairing, so sorting it into a `_` form would destroy that, and a `/` form would
over-claim positions. Species level is the only honest option.

Built on `lsi_class_mapping.csv` (w1's survey, validated in `final-003`).
"""
from __future__ import annotations

import re
from functools import lru_cache

# Class translation. Left: what the pipeline writes. Right: the LSI head, with `{}` where the body
# goes. Sphingoid classes take an oxygen designator rather than the d/t letter — `SM d34:1` becomes
# `SM 34:1;O2` — because `d`/`t` count hydroxyls and LSI states them as `;O2`/`;O3`.
BASE_OXYGEN = {"d": ";O2", "t": ";O3", "m": ";O"}

CLASS = {
    "Cer": "Cer", "Cer[NS]": "Cer", "Cer[NDS]": "Cer", "Cer[NP]": "Cer",
    "Cer[AS]": "Cer", "Cer[ADS]": "Cer", "Cer[AP]": "Cer", "Cer[BS]": "Cer",
    "Cer[BDS]": "Cer", "Cer[EOS]": "Cer",
    "HexCer": "HexCer", "HexCer[NS]": "HexCer", "HexCer[NDS]": "HexCer", "HexCer[AP]": "HexCer",
    "GlcCer": "HexCer", "GlcCer[NS]": "HexCer", "GlcCer[NDS]": "HexCer", "GlcCer[AP]": "HexCer",
    "SM": "SM", "CerP": "CerP", "SPB": "SPB", "SPBP": "SPBP",
    "PC": "PC", "PE": "PE", "PS": "PS", "PG": "PG", "PI": "PI", "PA": "PA",
    "LysoPC": "LPC", "LysoPE": "LPE", "LysoPS": "LPS", "LysoPG": "LPG",
    "LysoPI": "LPI", "LysoPA": "LPA",
    "TG": "TG", "DG": "DG", "MG": "MG", "CE": "CE", "FA": "FA", "AC": "CAR",
    "CL": "CL", "MLCL": "MLCL",
}

# Ether classes. The pipeline spells the linkage out in the class; LSI puts `O-` or `P-` on the
# body. `Plasmanyl` is the alkyl ether (`O-`); `Plasmenyl` is the vinyl ether (`P-`).
ETHER = {
    "Plasmanyl-PC": ("PC", "O-"), "Plasmenyl-PC": ("PC", "P-"),
    "Plasmanyl-PE": ("PE", "O-"), "Plasmenyl-PE": ("PE", "P-"),
    "Plasmanyl-PS": ("PS", "O-"), "Plasmenyl-PS": ("PS", "P-"),
    "Plasmanyl-TG": ("TG", "O-"), "Plasmenyl-TG": ("TG", "P-"),
    "Alkanyl-TG": ("TG", "O-"), "Alkenyl-TG": ("TG", "P-"),
    "Alkanyl-DG": ("DG", "O-"), "Alkenyl-DG": ("DG", "P-"),
    "LysoPC O": ("LPC", "O-"), "LysoPE O": ("LPE", "O-"),
}

SPECIES_LEVEL_ONLY = {"CL", "MLCL"}          # see the module docstring

_CHAIN = re.compile(r"(\d+):(\d+)")


@lru_cache(maxsize=4096)
def shorthand(name: str, source: str = "") -> str:
    """LSI shorthand for one pipeline name. Empty string when it cannot be stated honestly.

    `source` is the row's `Identification Source`. An `RT model` row has no MS2 behind it, so it
    can never carry `P-` — with no fragments there is nothing to distinguish a vinyl ether from an
    alkyl one, and asserting the harder claim from the weaker evidence is exactly backwards.
    """
    name = (name or "").strip()
    if not name:
        return ""

    head, _, body = name.partition(" ")
    if not body:
        return ""

    ether = ""
    if head in ETHER:
        lsi_head, ether = ETHER[head]
        if ether == "P-" and source == "RT model":
            # No spectrum: state the weaker, defensible form.
            ether = "O-"
    elif head in CLASS:
        lsi_head = CLASS[head]
    else:
        base = re.sub(r"\[[^\]]*\]", "", head)         # Cer[XYZ] -> Cer, for subclasses not listed
        lsi_head = CLASS.get(base, "")
        if not lsi_head:
            return ""

    # An `O-`/`P-` already written into the body (`PC O-38:2`) is the same claim; do not double it.
    body = body.strip()
    m = re.match(r"^([OP])-\s*", body)
    if m:
        if not ether:
            ether = f"{m.group(1)}-"
        elif ether[0] != m.group(1) and source == "RT model":
            ether = "O-"
        body = body[m.end():]

    oxygen = ""
    # ⚠ A hydroxyl chain prefix must be read BEFORE anything else, and it was not read at all.
    #
    # LipiDex writes hydroxylated glycerophospholipids as `PC[OH] OH-16:0_18:2`. `OH-` is matched
    # by neither the ether pattern above (`^([OP])-` needs O or P immediately before the hyphen)
    # nor the sphingoid pattern below, so it was skipped in silence and the chains parsed as if
    # nothing were attached: `PC[OH] OH-16:0_18:2` exported as `PC 16:0_18:2`. That is a
    # hydroxylated lipid reported as its unmodified parent — 331 rows across the staged validation
    # results — and it makes our output score as a miss against a deposit's `PC …;O` and as a
    # FALSE MATCH against its ordinary `PC`.
    m = re.match(r"^\(?OH\)?(\d*)-\s*", body)
    if m:
        count = int(m.group(1)) if m.group(1) else 1
        oxygen = ";O" if count == 1 else f";O{count}"
        body = body[m.end():]
    else:
        m = re.match(r"^([dtm])(?=\d)", body)
        if m:
            oxygen = BASE_OXYGEN[m.group(1)]
            body = body[1:]

    chains = _CHAIN.findall(body)
    if not chains:
        return ""

    if lsi_head in SPECIES_LEVEL_ONLY or len(chains) == 1:
        total_c = sum(int(c) for c, _ in chains)
        total_d = sum(int(d) for _, d in chains)
        return f"{lsi_head} {ether}{total_c}:{total_d}{oxygen}"

    # ⚠ `_`, never `/`. The pipeline proves chains from fragments and never proves sn-position.
    # The chain ORDER is preserved rather than sorted: for cardiolipin it encodes sn-pairing, and
    # for everything else the written order is what the library asserted.
    joined = "_".join(f"{int(c)}:{int(d)}" for c, d in chains)
    return f"{lsi_head} {ether}{joined}{oxygen}"


def sum_shorthand(name: str, source: str = "") -> str:
    """Species-level form: chains collapsed to one total. `PC 16:0_18:1` -> `PC 34:1`.

    Published lists mix the two levels, so a comparison needs both. Kept separate rather than
    replacing `shorthand`, because collapsing discards the chain evidence the pipeline worked for.
    """
    full = shorthand(name, source)
    if not full:
        return ""
    head, _, body = full.partition(" ")
    chains = _CHAIN.findall(body)
    if len(chains) <= 1:
        return full
    ether = "O-" if "O-" in body else ("P-" if "P-" in body else "")
    oxygen = ""
    m = re.search(r";O\d*$", body)
    if m:
        oxygen = m.group(0)
    return (f"{head} {ether}{sum(int(c) for c, _ in chains)}:"
            f"{sum(int(d) for _, d in chains)}{oxygen}")
