"""One name for one lipid class, whatever the source called it.

Class names are not standardised, and the mismatches are silent — they read as a missing class
rather than as a spelling difference. Four have already cost real work here:

    LysoPA  vs  LPA          LipidBlast against MS-DIAL
    AC      vs  CAR / ACar   our libraries against MS-DIAL
    GlcCer  vs  HexCer       ⚠ our OWN two libraries, for the same molecules
    ST      vs  BA           MS-DIAL files bile acids under its sterol prefix

The third is the expensive one. The retention model fits ONE SURFACE PER CLASS and needs six
members; one class arriving under two names splits its anchors and can leave both halves under the
minimum, losing the model for that class entirely. So this is not cosmetic tidying — a canonical
class is what the retention model must group by.

⚠ **It is worse than a between-tools problem.** A single MS-DIAL-annotated deposit (ST004797)
writes `TG` and `TAG`, `DG` and `DAG`, `HexCer` and `GlcCer`, and `BA` and `ST` — all in one output
table. There is no stable "MS-DIAL convention" to normalise against, which is why this maps to a
canonical form rather than translating between conventions.

## What is deliberately NOT here

**No structure identifiers are invented.** InChIKey, SMILES and LIPID MAPS IDs name ONE defined
structure, and most identifications here are sum compositions — `PC 34:1` is a family, not a
molecule. An identifier can only be carried through when the source library supplied one AND the
identification is molecular-level. LipidBlast and LipiDex supply none at all.

**No SE -> CE mapping.** A steryl ester is not always a cholesteryl ester, and our libraries carry
CE only. Mapping them would credit us with sterol esters we cannot identify; the gap is real and
belongs in the coverage table.

## Provenance

The shorthand this pipeline emits (`lsi.py`) already IS the LIPID MAPS convention — Liebisch et
al. 2020, "Update on LIPID MAPS classification, nomenclature, and shorthand notation for MS-derived
lipid structures", J Lipid Res 61:1539. There is no separate LIPID MAPS naming to convert to.
"""
from __future__ import annotations

import re

#: variant -> canonical. Lower-cased keys are matched case-insensitively.
CANONICAL: dict[str, str] = {
    # hexosylceramides. HexCer is the honest name: glucosyl- and galactosylceramide are
    # stereoisomers that reversed-phase chromatography does not separate and whose MS2 spectra are
    # indistinguishable, so `GlcCer` asserts what the evidence cannot support.
    "glccer": "HexCer", "galcer": "HexCer", "glucer": "HexCer", "hex1cer": "HexCer",
    "glucosylceramide": "HexCer", "galactosylceramide": "HexCer",
    "laccer": "Hex2Cer", "hex2cer": "Hex2Cer", "lactosylceramide": "Hex2Cer",
    # acylcarnitines
    "car": "AC", "acar": "AC", "carnitine": "AC", "carnitinec": "AC", "acylcarnitine": "AC",
    # glycerolipids
    "tag": "TG", "dag": "DG", "mag": "MG", "triacylglycerol": "TG", "diacylglycerol": "DG",
    # lyso forms
    "lysopc": "LPC", "lysope": "LPE", "lysopi": "LPI", "lysops": "LPS",
    "lysopg": "LPG", "lysopa": "LPA", "lysosm": "LSM",
    # cholesteryl esters
    "cholesterylester": "CE", "che": "CE",
    # bile acids. MS-DIAL names them with its sterol prefix and separates them only in
    # COMPOUNDCLASS; on a name-derived class that would put them on the sterols' retention surface,
    # which they resemble not at all.
    "bileacid": "BA", "basulfate": "BA",
}

_SUBCLASS = re.compile(r"\[[^\]]*\]$")
# ⚠ Gangliosides carry their subclass as a SUFFIX, not a bracket, and so were invisible to the
# normalisation every other subclass gets. Our libraries hold 748 `GM3-NANA` and 748 `GM3-NGNA`
# entries; a deposit writes plain `GM3`; nothing matched, and the class read as one we did not
# cover at all — the fifth naming mismatch found this way, after LysoPA/LPA, AC/CAR, GlcCer/HexCer
# and ST/BA.
#
# NANA and NGNA are different molecules — N-acetyl- against N-glycolyl-neuraminic acid, 16 Da
# apart, and rodents make NGNA where humans make only NANA. So the distinction is preserved as a
# bracket tag rather than discarded; it simply stops hiding the base class.
_SIALIC = re.compile(r"^(G[MDTQ]\d[a-c]?)-(NANA|NGNA)$", re.I)
# ⚠ Underscore AND hyphen. MS-DIAL writes `GlcCer_NS`; MetaboLights deposits write `Cer-NS` and
# `HexCer-NS`. Only the underscore was matched, so the hyphenated dialect fell through to
# CANONICAL as a whole unknown token and `Cer-NS` stayed `Cer-NS` — reported as a class we hold no
# library for, while thousands of `Cer[NS]` entries sat in the libraries unmatched.
_MSDIAL_SUBCLASS = re.compile(r"[_-](NS|NDS|AS|ADS|AP|NP|BS|BDS|EOS|HS|HDS)$")


def canonical_class(name: str) -> str:
    """The canonical class for a class token or a full lipid name.

    Subclass tags are preserved, not flattened: `Cer[ADS]` and `Cer[NP]` are different molecules
    from different enzymes and must not merge. MS-DIAL's underscore form (`GlcCer_NS`) is
    rewritten to the bracket form so the two dialects agree.
    """
    token = (name or "").strip().split()[0] if name and " " in name else (name or "").strip()
    if not token:
        return ""
    token = _MSDIAL_SUBCLASS.sub(lambda m: f"[{m.group(1)}]", token)
    token = _SIALIC.sub(lambda m: f"{m.group(1).upper()}[{m.group(2).upper()}]", token)
    tag = _SUBCLASS.search(token)
    head = _SUBCLASS.sub("", token)
    mapped = CANONICAL.get(head.lower(), head)
    return f"{mapped}{tag.group(0)}" if tag else mapped


def canonical_name(name: str) -> str:
    """The full lipid name with only its class token canonicalised. Chains are left untouched."""
    s = (name or "").strip()
    if " " not in s:
        return canonical_class(s)
    head, _, rest = s.partition(" ")
    return f"{canonical_class(head)} {rest}"


def variants_of(canonical: str) -> list[str]:
    """Every spelling known to map to this canonical class, for a report or a lookup table."""
    return sorted({k for k, v in CANONICAL.items() if v == canonical} | {canonical.lower()})
