"""The adduct database, and the mass tests the peak finder's filters are built on.

`lib_gen/Adduct.java` plus the `checkAdduct*` / `checkDimer` / `checkFragment` helpers in
`peak_finder/Utilities.java`.

The database the peak finder actually reads is **`src/peak_finder/Possible_Adducts.csv`**, not
the `Adducts.csv` inside each library folder. They are not the same list: the peak finder's
carries `[M+K]+` and `[M+C2H6N2]+`, which the library files lack, and omits `[M+FA-H]-` and
`[M-2H]2-`, which they have. Potassium adducts are common enough in positive mode that using
the library list instead measurably under-detects adducts and dimers.

This is what makes the adduct sweep work at all. One molecule appears in the feature table
several times over — as `[M+H]+`, `[M+Na]+`, `[M+NH4]+`, as a dimer, and as whatever survives
in-source fragmentation — and every one of those is a separate compound group competing for the
same row in the result file. Without a formula for each adduct there is no way to ask whether
two quant ions could be the same neutral molecule, and the sweep degenerates into comparing
mass differences against a hard-coded list.

The test is always the same shape: strip the adduct's formula mass off each quant ion and see
whether the two neutral masses agree to within 20 ppm.

Reproduced from the semantics rather than transliterated: the original parses a formula by
string surgery over an element table, this uses a regular expression over the same table. The
element masses are the ones LipiDex carries, which are ordinary monoisotopic values.
"""
from __future__ import annotations

import csv
import re
from dataclasses import dataclass
from pathlib import Path

MAX_PPM_DIFF = 20.0

# peak_finder/Utilities.java: elementsArray and masses, in the same order. The Xa/Xb/Xc/Xd
# placeholders stand in for heavy isotopes after `removeHeavyElements` rewrites (2H), (13C),
# (15N) and (18O); they are kept so a labelled formula still resolves.
ELEMENT_MASSES = {
    "H": 1.007825035, "Xa": 2.014101779, "Li": 7.016003, "B": 11.0093055, "C": 12.0,
    "Xb": 13.00335483, "N": 14.003074, "Xc": 15.00010897, "O": 15.99491463, "Xd": 17.9991603,
    "F": 18.99840322, "Na": 22.9897677, "Mg": 23.9850423, "P": 30.973762, "S": 31.9720707,
    "Cl": 34.96885272, "K": 38.9637074, "Ca": 39.9625906, "Cr": 51.9405098, "Mn": 54.9380471,
    "Fe": 55.9349393, "Ni": 57.9353462, "Co": 58.9331976, "Cu": 62.9295989, "Zn": 63.9291448,
    "As": 74.9215942, "Br": 78.9183361, "Se": 79.9165196, "Mo": 97.9054073, "Pd": 105.903478,
    "Ag": 106.905092, "Cd": 113.903357, "I": 126.904473, "Au": 196.966543, "Hg": 201.970617,
}

_HEAVY = {"(2H)": "Xa", "(13C)": "Xb", "(15N)": "Xc", "(18O)": "Xd"}
_TOKEN = re.compile(r"([A-Z][a-z]?)(-?\d*)")


def formula_mass(formula: str) -> float:
    """Monoisotopic mass of a formula, counts explicit and possibly negative.

    `H1` is a proton gained, `H-1` one lost, `O-1H-1` a water loss less a proton. A bare element
    with no count means one. Unknown elements contribute nothing, as in the original.
    """
    if not formula:
        return 0.0
    for heavy, placeholder in _HEAVY.items():
        formula = formula.replace(heavy, placeholder)
    total = 0.0
    for element, count in _TOKEN.findall(formula):
        if element not in ELEMENT_MASSES:
            continue
        total += ELEMENT_MASSES[element] * (int(count) if count not in ("", "-") else 1)
    return total


def ppm_diff(mass1: float, mass2: float) -> float:
    if mass2 == 0:
        return float("inf")
    return abs(mass1 - mass2) / abs(mass2) * 1e6


@dataclass(slots=True)
class Adduct:
    name: str
    formula: str
    loss: bool
    polarity: str
    charge: int
    mass: float = 0.0

    def __post_init__(self):
        self.mass = formula_mass(self.formula)


def read_adducts(path: str | Path) -> list[Adduct]:
    """Read an `Adducts.csv` from a LipiDex library folder."""
    out: list[Adduct] = []
    with Path(path).open(newline="", errors="replace") as fh:
        for row in csv.DictReader(fh):
            name = (row.get("Name") or "").strip()
            if not name:
                continue
            try:
                charge = int(float(row.get("Charge") or 1))
            except ValueError:
                charge = 1
            out.append(Adduct(name=name, formula=(row.get("Formula") or "").strip(),
                              loss=(row.get("Loss") or "").strip().upper() == "TRUE",
                              polarity=(row.get("Polarity") or "").strip(),
                              charge=charge))
    return out


DEFAULT_ADDUCTS = Path(__file__).resolve().parents[2] / "data/lipidex_src/peak_finder/Possible_Adducts.csv"


def load_adducts(sources=None) -> list[Adduct]:
    """Load the adduct database.

    Accepts the peak finder's `Possible_Adducts.csv` directly, or library folders containing an
    `Adducts.csv`, or nothing — in which case the peak finder's own list is used, which is what
    LipiDex does.
    """
    if not sources:
        sources = [DEFAULT_ADDUCTS] if DEFAULT_ADDUCTS.exists() else []
    by_name: dict[str, Adduct] = {}
    for source in sources:
        path = Path(source)
        if path.is_dir():
            path = path / "Adducts.csv"
        if not path.exists():
            continue
        for adduct in read_adducts(path):
            by_name.setdefault(adduct.name, adduct)
    return list(by_name.values())


def check_adduct_known(adducts: list[Adduct], mass1: float, mass2: float,
                       mass1_adduct: str | None) -> bool:
    """Could mass2 be another adduct of the molecule mass1 was identified as?

    The identification fixes which adduct mass1 is, so its neutral mass is known exactly;
    mass2 is then tried against every adduct in the database.
    """
    if not mass1_adduct:
        return False
    reference = next((a for a in adducts if a.name == mass1_adduct), None)
    if reference is None:
        return False
    neutral1 = mass1 - reference.mass
    return any(ppm_diff(mass2 - a.mass, neutral1) < MAX_PPM_DIFF for a in adducts)


def check_adduct_unknown(adducts: list[Adduct], mass1: float, mass2: float) -> bool:
    """Could the two masses be two *different* adducts of one unidentified molecule?"""
    for i, first in enumerate(adducts):
        neutral1 = mass1 - first.mass
        for j, second in enumerate(adducts):
            if i == j:
                continue
            if ppm_diff(mass2 - second.mass, neutral1) < MAX_PPM_DIFF:
                return True
    return False


def check_dimer(adducts: list[Adduct], mass1: float, mass2: float, polarity: str) -> bool:
    """Is one of the two a dimer of the other, sharing an adduct?"""
    for adduct in adducts:
        if adduct.polarity != polarity:
            continue
        if ppm_diff((mass1 - adduct.mass) * 2.0, mass2 - adduct.mass) < MAX_PPM_DIFF:
            return True
        if ppm_diff(mass1 - adduct.mass, (mass2 - adduct.mass) * 2.0) < MAX_PPM_DIFF:
            return True
    return False


def check_fragment(fragment_masses, quant_ions) -> bool:
    """Is another group's quant ion one of this lipid's own predicted fragments?

    That is an in-source fragment: the molecule broke apart in the source, and the piece was
    picked up and quantified as though it were a compound in its own right. The fragment list
    is the `Potential Fragments` column, which is why that column has to be populated for this
    filter to do anything at all.
    """
    for fragment in fragment_masses or ():
        for ion in quant_ions or ():
            if ppm_diff(fragment, ion) < MAX_PPM_DIFF:
                return True
    return False
