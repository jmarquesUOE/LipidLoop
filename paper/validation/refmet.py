"""Map lipid names onto RefMet, the Metabolomics Workbench's standardised nomenclature.

Deliberately in the COMPARISON layer, not in the pipeline. RefMet lives behind a web service, and
identification must not depend on an external host being up — a run that names lipids differently
because a server was slow is not reproducible. Here it only ever affects a comparison, and a
failed fetch degrades to "no RefMet name" rather than to a different answer.

## Why it is worth having at all

Class names are not standardised, and a single MS-DIAL-annotated deposit writes `TG` and `TAG`,
`DG` and `DAG`, `HexCer` and `GlcCer`, and `BA` and `ST` — in one output table. So there is no
"their convention" to normalise against. RefMet is the one shared target: Workbench applies it to
every deposit it holds, which makes it the only nomenclature under which our identifications and a
public deposit's are written the same way by construction rather than by our own mapping.

## ⚠ What a RefMet match does and does not mean

RefMet names lipids at the level the evidence supports, mostly in Liebisch shorthand — so
`PC 34:1` maps cleanly and carries `super_class`, `main_class` and `sub_class`. It also carries
`inchi_key` and `exactmass`, but ONLY for entries defined enough to have one: a sum composition is
a family, not a molecule, and its InChIKey field is empty. Never treat a RefMet hit as a structure.
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "src"))

import deposits  # noqa: E402
from lipidloop.nomenclature import canonical_class  # noqa: E402

URL = "https://www.metabolomicsworkbench.org/rest/refmet/all"
CACHE_FILE = "refmet_all.json"

_TABLE: dict[str, dict] | None = None


def _key(name: str) -> str:
    """Lookup key: case-folded, whitespace-collapsed, class token canonicalised.

    The canonicalisation is what makes `GlcCer 42:1;O2` and `HexCer 42:1;O2` reach the same RefMet
    row, and `CAR 16:0` reach `AC 16:0`'s.
    """
    s = re.sub(r"\s+", " ", (name or "").strip())
    if not s:
        return ""
    head, _, rest = s.partition(" ")
    if rest:
        s = f"{canonical_class(head)} {rest}"
    return s.lower()


def table() -> dict[str, dict]:
    """The whole RefMet table, keyed for lookup. Fetched once, then cached on disk."""
    global _TABLE
    if _TABLE is not None:
        return _TABLE
    path = deposits.CACHE / CACHE_FILE
    if not (path.exists() and path.stat().st_size > 0):
        deposits._get(URL, path)
    try:
        raw = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        _TABLE = {}
        return _TABLE
    rows = raw.values() if isinstance(raw, dict) else raw
    _TABLE = {}
    for row in rows:
        if not isinstance(row, dict):
            continue
        name = (row.get("name") or "").strip()
        if name:
            _TABLE.setdefault(_key(name), row)
    return _TABLE


def lookup(name: str) -> dict | None:
    """The RefMet row for a lipid name, or None. Tries the name, then its sum composition."""
    t = table()
    if not t:
        return None
    hit = t.get(_key(name))
    if hit:
        return hit
    # A molecular name may exist in RefMet only at species level: PC 16:0_18:1 -> PC 34:1.
    m = re.match(r"^(\S+)\s+(.*)$", (name or "").strip())
    if not m:
        return None
    chains = re.findall(r"(\d+):(\d+)", m.group(2))
    if len(chains) < 2:
        return None
    oxy = re.search(r";O\d*", m.group(2))
    total = f"{sum(int(c) for c, _ in chains)}:{sum(int(d) for _, d in chains)}"
    return t.get(_key(f"{m.group(1)} {total}{oxy.group(0) if oxy else ''}"))


def refmet_name(name: str) -> str:
    hit = lookup(name)
    return (hit.get("name") or "").strip() if hit else ""


def classification(name: str) -> tuple[str, str, str]:
    """(super_class, main_class, sub_class) — RefMet's hierarchy, or blanks."""
    hit = lookup(name)
    if not hit:
        return "", "", ""
    return (hit.get("super_class") or "", hit.get("main_class") or "",
            hit.get("sub_class") or "")


if __name__ == "__main__":
    t = table()
    print(f"  RefMet entries: {len(t):,}")
    for n in ("PC 34:1", "GlcCer 42:1;O2", "HexCer 42:1;O2", "CAR 16:0", "AC 16:0",
              "TAG 52:2", "PC 16:0_18:1", "Cer 42:1;O2", "BA 24:1;O3"):
        hit = lookup(n)
        sup, main, sub = classification(n)
        print(f"  {n:<18} -> {(hit.get('name') if hit else '(no match)'):<18} {sub or main or sup}")
