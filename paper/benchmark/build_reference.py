"""Build one harmonised table of published NIST SRM 1950 lipid lists, for the benchmark.

Sources (all at sum-composition level, values in nmol/mL where given):

  * Bowden 2017 (J Lipid Res 58:2275, NIST interlaboratory exercise). The 318 consensus values
    distributed inside LipidQC v1.0 (Ulmer 2017, Anal Chem 89:13069; `consensus values` sheet):
    254 lipids reported by >= 5 laboratories with COD <= 40 % (tier `Ref`) and 64 reported by 3
    or 4 laboratories (tier `Inf`). Supplemental Tables S1-S5 of the paper add the 86 lipids with
    >= 5 laboratories but COD > 40 % (tier `COD>40`); Ref + COD>40 is the paper's "339 lipids".
  * Quehenberger 2010 (J Lipid Res 51:3299, LIPID MAPS consortium), Supplemental Tables 1A-6B,
    parsed from the PDF text (`pdftotext -layout`). Only species-level lipids are kept; eicosanoids,
    sterols, dolichols and CoQ are recorded as `other` so they count against nothing.
  * Godzien 2024 (J Lipid Res 65:100671, ST003514), the published endogenous list (presence only,
    no concentrations).

Every entry is reduced to a KEY that our own identifications can be reduced to as well:
(family, kind, carbons, double bonds), where family is a class family (PC, PE, LPC, LPE, PI, PS,
PG, PA, SM, Cer, HexCer, TG, DG, CE, FA, ...), kind is `` for diacyl/acyl or `ether` for
plasmanyl/plasmenyl (a plasmenyl P-C:DB is written as ether C:DB+1, its formula equivalent, so
`PC O-34:1` and `PC P-34:0` are one key, as the consensus tables treat them). Isobaric groups such
as `PC O-34:1/P-34:0/33:1` keep one `group` id with one row per alternative: the group is
recovered when any alternative is identified.

    .venv/bin/python manuscript/benchmark/build_reference.py <raw_dir>

<raw_dir> holds `consensus_values.csv` (exported from LipidQC-v1.0.xlsm), `bowden_table_0..7.csv`
(the docx tables S1-S6), `queh_supp.txt` and `ST003514_published_lipids.csv`.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "reference" / "srm1950_reference_lists.csv"

FAMILY = {
    "FFA": "FA", "FA": "FA", "TAG": "TG", "TG": "TG", "DAG": "DG", "DG": "DG", "CER": "Cer", "CER[NS]": "Cer",
    "CER[NDS]": "Cer", "HEXCER": "HexCer", "HEXCER[NS]": "HexCer", "HEXCER[NDS]": "HexCer", "GLCCER": "HexCer",
    "GLCCER[NDS]": "HexCer", "GLCCER[NS]": "HexCer", "LACCER": "Hex2Cer", "CEROH": "CerOH", "CE": "CE", "CHE": "CE",
    "PC": "PC", "PE": "PE", "PI": "PI", "PS": "PS", "PG": "PG", "PA": "PA", "LPC": "LPC", "LYSOPC": "LPC",
    "LPE": "LPE", "LYSOPE": "LPE", "LPI": "LPI", "LYSOPI": "LPI", "SM": "SM", "PLASMANYL-PC": "PC",
    "PLASMENYL-PC": "PC", "PLASMANYL-PE": "PE", "PLASMENYL-PE": "PE", "ALKANYL-TG": "TG", "ALKENYL-TG": "TG",
    "CAR": "AC", "AC": "AC", "PC O": "PC", "PE O": "PE", "LPC O": "LPC", "LPE O": "LPE", "PC P": "PC", "PE P": "PE",
    "LPC P": "LPC", "LPE P": "LPE", "TG O": "TG", "TG P": "TG",
}

_CHAIN = re.compile(r"^(?P<cls>[A-Za-z][A-Za-z0-9\[\]\-]*?)\s*[\s(]\s*(?P<pre>d|t|O-|P-|OH-)?(?P<c>\d+):(?P<db>\d+)(?P<suf>[ep])?\)?(?P<rest>.*)$")


def key_of(name: str) -> tuple[str, str, int | None, int | None]:
    """One lipid name (any of the dialects above) -> (family, kind, C, DB); unparseable -> other."""
    s = name.strip()
    m = _CHAIN.match(s)
    if not m:
        return (s, "other", None, None)
    cls = m.group("cls").upper().replace("1,2-", "").replace("1,3-", "")
    fam = FAMILY.get(cls)
    if fam is None:
        return (s, "other", None, None)
    c, db = int(m.group("c")), int(m.group("db"))
    pre, suf = m.group("pre") or "", m.group("suf") or ""
    kind = ""
    if pre == "O-" or suf == "e" or cls.startswith("PLASMANYL") or cls == "ALKANYL-TG":
        kind = "ether"
    elif pre == "P-" or suf == "p" or cls.startswith("PLASMENYL") or cls == "ALKENYL-TG":
        kind, db = "ether", db + 1
    if fam == "Cer" and pre == "OH-":
        fam = "CerOH"
    if "[OH]" in cls:
        return (s, "other", None, None)
    return (fam, kind, c, db)


def ref_alternatives(name: str) -> list[str]:
    """`PC O-34:1/P-34:0/33:1` -> [`PC O-34:1`, `PC P-34:0`, `PC 33:1`]; `PE(36:2e)/PE(36:1p)` likewise."""
    parts = [p.strip() for p in name.split("/")]
    if len(parts) == 1:
        return parts
    head = parts[0].split(" ")[0] if " " in parts[0] else parts[0].split("(")[0]
    out = [parts[0]]
    for p in parts[1:]:
        out.append(p if p.upper().startswith(head.upper()) else f"{head} {p}")
    return out


def rows_for(source, tier, name, value, unc, n_labs, cod, group):
    out = []
    for alt in ref_alternatives(name):
        fam, kind, c, db = key_of(alt)
        out.append({"source": source, "tier": tier, "name": name, "alternative": alt, "family": fam, "kind": kind,
                    "carbons": c, "double_bonds": db, "group": group, "value_nmol_mL": value, "uncertainty": unc,
                    "n_labs": n_labs, "cod_pct": cod})
    return out


def bowden(raw: Path) -> list[dict]:
    out = []
    seen = set()
    for r in csv.DictReader((raw / "consensus_values.csv").open()):
        name = r["Measureand"]
        if not name or name.startswith("Total"):
            continue
        unit = r["UnitNumerator"]
        value = float(r["Value"]) * (1e-3 if unit == "pmol" else 1.0)
        unc = float(r["StandardUncertainty"] or 0) * (1e-3 if unit == "pmol" else 1.0)
        tier = "Ref" if r["CertRef"] == "Ref" else "Inf"
        seen.add(name.upper())
        out += rows_for("Bowden2017", tier, name, value, unc, int(float(r["NumLabs"] or 0)), r["CoefficientOfDispersion(%)"], f"B:{name}")
    for i in range(5):                                # docx tables 0-4 = Supplemental Tables S1-S5
        for r in csv.DictReader((raw / f"bowden_table_{i}.csv").open()):
            name = r["Lipid"].strip()
            if not name or name.upper() in seen:
                continue
            unit = r["Units"]
            value = float(r["Consensus Location"]) * (1e-3 if unit == "pmol/mL" else 1.0)
            unc = float(r["Standard Uncertainty"]) * (1e-3 if unit == "pmol/mL" else 1.0)
            seen.add(name.upper())
            out += rows_for("Bowden2017", "COD>40", name, value, unc, int(r["# of Labs"]), r["COD (%)"], f"B:{name}")
    return out


_SECTION = re.compile(r"^Supplemental Table (\d[A-D]?)")
_NUM = r"([\d.]+)\s+([\d.]+)"


def quehenberger(raw: Path) -> list[dict]:
    out, section, base = [], "", ""
    acc: dict[str, list] = {}

    def add(name, mean, sem, tag=""):
        alts = ref_alternatives(name)
        k = "/".join(alts)
        if k in acc:                                   # ω-isomers of one FA, 1,2- and 1,3-DG: summed
            acc[k][1] += mean; acc[k][2] = (acc[k][2] ** 2 + sem ** 2) ** 0.5
        else:
            acc[k] = [name, mean, sem]

    for line in (raw / "queh_supp.txt").read_text().splitlines():
        m = _SECTION.match(line.strip())
        if m:
            section = m.group(1); base = ""; continue
        s = line.strip()
        if not s:
            continue
        if section == "1A":
            m = re.search(r"\s(\d+:\d+)\s*(\([^)]*\))?\s+" + _NUM, s)
            if m and not s.startswith("Chain"):
                add(f"FA {m.group(1)}", float(m.group(3)), float(m.group(4)))
        elif section in ("2A", "5C"):
            m = re.match(r"(TG|CE)\((\d+:\d+)\)\s+" + _NUM, s)
            if m:
                add(f"{m.group(1)} {m.group(2)}", float(m.group(3)), float(m.group(4)))
        elif section in ("2B", "2C"):
            m = re.match(r"1,[23]-DG\(\s*(\d+:\d+)\)\s+" + _NUM, s)
            if m:
                add(f"DG {m.group(1)}", float(m.group(2)), float(m.group(3)))
        elif section == "3A":
            m = re.match(r"((?:L?P[ACEGIS]\(\d+:\d+[ep]?\)/?)+)\d?\s+" + _NUM, s)
            if m:
                add(m.group(1), float(m.group(2)), float(m.group(3)))
        elif section == "4A":
            m = re.match(r"C(\d+:\d+)\s+" + _NUM, s)
            if m:
                add(f"SM d{m.group(1)}", float(m.group(2)), float(m.group(3)))
        elif section == "4B":
            m = re.match(r"d(\d+):(\d)\d?\s+C(\d+):(\d+)\s+" + _NUM, s)
            if m:
                c = int(m.group(1)) + int(m.group(3)); db = int(m.group(2)) + int(m.group(4))
                add(f"HexCer d{c}:{db}", float(m.group(5)), float(m.group(6)))
        elif section == "4C":
            m = re.match(r"d(\d+):(\d)\d?$", s)
            if m:
                base = f"{m.group(1)}:{m.group(2)}"; continue
            m = re.match(r"C(\d+):(\d+)\s+" + _NUM, s)
            if m and base:
                bc, bdb = base.split(":")
                add(f"Cer d{int(bc) + int(m.group(1))}:{int(bdb) + int(m.group(2))}", float(m.group(3)), float(m.group(4)))
    for name, mean, sem in acc.values():
        out += rows_for("Quehenberger2010", "LIPIDMAPS", name, mean, sem, 1, "", f"Q:{name}")
    return out


def _godzien_sum(name: str) -> str:
    """`Cer 18:1;O2/16:0` -> `Cer 34:1`; `PC O-16:0/18:1` -> `PC O-34:1`; `SM 42:4;O2` -> `SM 42:4`."""
    name = re.sub(r";O\d?", "", name)
    head, _, chains = name.partition(" ")
    parts = chains.split("/")
    if len(parts) < 2:
        return name
    c = db = 0; pre = ""
    for part in parts:
        m = re.fullmatch(r"([OP]-)?(\d+):(\d+)", part.strip())
        if not m:
            return name
        pre = pre or (m.group(1) or ""); c += int(m.group(2)); db += int(m.group(3))
    return f"{head} {pre}{c}:{db}"


def godzien(raw: Path) -> list[dict]:
    out = []
    for r in csv.DictReader((raw / "ST003514_published_lipids.csv").open()):
        if r["ORIGIN"] != "endogenous":
            continue
        name = r["SHORTHAND NOTATION"].strip() or r["SYNONYM NAME: SUM COMPOSITION"].strip()
        name = _godzien_sum(name)
        out += rows_for("Godzien2024", "ST003514", name, "", "", 1, "", f"G:{name}")
    return out


def main(raw: Path) -> None:
    rows = bowden(raw) + quehenberger(raw) + godzien(raw)
    with OUT.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    import collections
    for src in ("Bowden2017", "Quehenberger2010", "Godzien2024"):
        sub = [r for r in rows if r["source"] == src]
        groups = {r["group"] for r in sub}
        tiers = collections.Counter(next(x for x in sub if x["group"] == g)["tier"] for g in groups)
        fams = collections.Counter(r["family"] if r["kind"] != "other" else "other" for r in sub)
        print(f"{src}: {len(groups)} entries {dict(tiers)}")
        print("   ", dict(sorted(fams.items(), key=lambda kv: -kv[1])))
        print("    other:", sorted({r['alternative'] for r in sub if r['kind'] == 'other'})[:40])
    print("written", OUT)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
