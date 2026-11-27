"""Generate positive-mode entries for the prenol quinones — ubiquinone (CoQ) and ubiquinol.

**A whole LIPID MAPS category that no library we search contains.** Prenol lipids are category PR
of the LIPID MAPS classification, and `qc_report.py` has had a colour reserved for them since the
category plot was written; nothing has ever filled it. Zero quinone entries across every library
in the default set, so coenzyme Q has never reached a delivered table from this pipeline, or from
the Compound Discoverer + LipiDex workflow it replaces. It is not absent from the samples. On the
in-house Q Exactive Plus positive-mode data ubiquinone is a **better peak than
the free-cholesterol positive control**: max/median 812 against cholesterol's 45, present in all
ten injections examined.

That is the same kind of gap as the 1-deoxyceramides (the developer notes (DEOXYSPHINGOLIPIDS, not distributed)): absent by
construction, not by abundance.

## What is written, and what is deliberately not

Four entries, all confirmed in our own data by MS2 with the quinone head group as base peak:

    CoQ10   [M+NH4]+   880.7177   RT 21.07    the dominant ion, 12 qualifying spectra
    CoQ10   [M+H]+     863.6912   RT 21.07    10 spectra
    CoQ10H2 [M+NH4]+   882.7334   RT 19.5     ubiquinol, the reduced form, 11 spectra
    CoQ9    [M+NH4]+   812.6551   RT 19.1     7 spectra

`[M+Na]+` is seen in MS1 for both forms at the same retention (885.6731 and 887.6888) and is
**not written**: not one sodiated spectrum passed the MS2 screen, so there is nothing measured to
write and a sodiated fragmentation is not the protonated one. DMQ10, the COQ7 substrate, is **not
written either** — it was looked for and is not reliably there (max/median 15.6, retention
scattered by 3 min, i.e. noise), and an entry for a compound we cannot detect can only take
matches from something else. `--predicted-homologues` adds CoQ6-CoQ11 and their quinols by
n x C5H8; they are predictions, off by default, and the reason is below.

## ⚠⚠ The trap: ubiquinol is 10 ppm from the ubiquinone M+2 isotope

    CoQ10   [M+NH4]+  M+2 isotope   882.7244       (880.7177 + 2 x 1.0033548)
    CoQ10H2 [M+NH4]+  monoisotopic  882.7334
                                    --------
                                    0.0090 Da  =  10.2 ppm at m/z 882

`search.MS1_TOL` is **0.01 Da**, absolute, so the isotope falls INSIDE the window that selects
library candidates for a precursor: the ubiquinol entry below is a candidate for every MS2 taken
on the ubiquinone isotope. MS2 cannot break the tie either — both ions give the same quinone head
group at 197.0808 as base peak, because it is the same molecule's head group. The isolation width
is ~100x the mass difference, so the two hypotheses are never isolated apart.

**They separate on retention and on nothing else.** Measured here, three files, apexes on a common
grid:

    CoQ10   [M+NH4]+  880.7177   apex 21.05-21.09
    CoQ10   M+2       882.7244   apex 21.05-21.08     contributes ZERO at 19.5
    CoQ10H2 [M+NH4]+  882.7334   apex 19.50-19.60     ubiquinol, ~1.5 min earlier
    CoQ10H2 [M+Na]+   887.6888   apex 19.52-19.56     the adduct agrees, which the isotope cannot

The reduced form carries two hydroxyls where the oxidised form carries two ketones, so it is more
polar and elutes EARLIER on reversed phase. The direction is chemistry, not an artefact.

**This library cannot enforce that on its own.** The .msp format has no retention field, and
`peakfinder._filter_pair` removes an isotope only at M+1 spacing and only when the heavier peak
carries no identification (known limitations of the peak finder: isotope promotion, and an isotope rule applied only to unidentified peaks) — an identified ubiquinol row is exactly
the case it skips. So the constraint is carried in three places instead:

  * `data/libraries/Ubiquinones_RT_anchors.csv`, the measured retention of every ion written here;
  * `scripts/check_ubiquinol_rt.py`, which tests a delivered table against it and names any
    ubiquinol row that co-elutes with ubiquinone as the isotope it is;
  * `RT=` and `RTWindow=` in each entry's `Comment`, so the number travels with the entry. Nothing
    downstream reads them — they are for whoever opens the file.

**Run the checker on any study where ubiquinol is reported.** A ubiquinol row at the ubiquinone
retention is the isotope, whatever its dot product.

## The intensities are measured, the masses are computed

Repo rule, and both halves matter. Intensities come from `cq04_coq_ms2_peaklist.py` over ten
injections of an in-house Q Exactive Plus study: every MS2 whose precursor
is within 20 ppm of the target, whose retention is in the target's window, and **whose base peak
is the quinone head** — that last gate matters, without it a minority of co-isolated spectra
contributed 184.0733 (phosphocholine) and a 292/293/320/321 series at median intensity 900+,
which is a neighbour's spectrum and not a weak CoQ fragment. Each spectrum scaled to base peak
999, fragments binned at 15 ppm, median intensity across spectra, kept at >= 70 % presence.

Masses are computed from the assigned formula rather than copied from the measurement, because
the measured low-mass ions run +10 to +16 ppm on this instrument (about 1 mDa at m/z 81 — well
inside `MS2_TOL` 0.01 Da, but there is no reason to write the drift into a library).

    python scripts/build_ubiquinone_library.py --out data/libraries/Ubiquinones_Positive.msp

⚠ **Not in DEFAULT_LIBRARIES**, and that is the project rule rather than timidity: a library goes
into the default set only after its displacement against matched decoys has been measured
(the developer notes (LIBRARY_VERSIONS, not distributed); `Cer[EOS]` and `FAHFA` were both withdrawn on that test). These entries
also cannot take the usual CH2 decoy construction — `make_decoy_library.py` shifts chain fragments
identified by their `[c:d]` annotation, and a quinone has no acyl chain to shift — so a decoy for
them needs its own construction first. Searching them meanwhile is one flag — `--extra-library`,
which appends to the modifier's default set and records the choice in that run's own
configuration. ⚠ Not `--set libraries=[...]`: `run_study.py` rebuilds `config["libraries"]` from
the stem list a few lines after reading the overrides, so that form is silently discarded.

    python scripts/run_study.py STUDY --extra-library Ubiquinones_Positive

the developer notes (PRENOL_QUINONES, not distributed) has the full account.
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path

# Monoisotopic masses, as everywhere else in this repo.
EM = {"C": 12.0, "H": 1.0078250319, "N": 14.0030740052, "O": 15.9949146221, "Na": 22.9897692809}
ELECTRON = 0.00054858
PROTON = EM["H"] - ELECTRON
ISOPRENE = 5 * EM["C"] + 8 * EM["H"]        # 68.0626, the homologue increment
CARBON_13 = 1.0033548                       # 13C - 12C, for the isotope arithmetic in the docs

ADDUCTS = {
    "[M+H]+": PROTON,
    "[M+NH4]+": EM["N"] + 4 * EM["H"] - ELECTRON,
    "[M+Na]+": EM["Na"] - ELECTRON,
}

PRECURSOR_CUTOFF = 2.0      # the library reader drops peaks closer than this to the precursor
MIN_MS2_MASS = 61.0         # below the m/z floor a real acquisition records

_TOKEN = re.compile(r"([A-Z][a-z]?)(\d*)")


def neutral_mass(formula: str) -> float:
    """Monoisotopic mass of a neutral CHNO formula."""
    return sum(EM[element] * (int(count) if count else 1)
               for element, count in _TOKEN.findall(formula) if element)


def cation_mass(formula: str) -> float:
    """Monoisotopic m/z of a singly charged cation written as its formula, e.g. `C10H13O4`.

    ⚠ The electron. A fragment cation is the neutral formula MINUS one electron, and leaving that
    out puts every constant fragment 0.55 mDa high. That is 6.5 ppm at m/z 85 and it is the exact
    mistake `build_acylcarnitine_library.py` made and documented, so it is not repeated here.
    """
    return neutral_mass(formula) - ELECTRON


#: Neutral species. CoQ-n carries n isoprene units; the quinol is the quinone plus 2 H.
#: CoQ10 C59H90O4 862.6839, CoQ9 C54H82O4 794.6213, ubiquinol-10 C59H92O4 864.6996 — the three
#: that `cq01_xic_pos.py` and `cq03_ubiquinol_check.py` searched for and found.
QUINONES = {
    "CoQ9": ("C54H82O4", 9, "ubiquinone-9, oxidised"),
    "CoQ10": ("C59H90O4", 10, "ubiquinone-10, oxidised"),
    "CoQ9H2": ("C54H84O4", 9, "ubiquinol-9, reduced"),
    "CoQ10H2": ("C59H92O4", 10, "ubiquinol-10, reduced"),
}

#: Fragment intensities MEASURED on an in-house Q Exactive Plus study, 10 injections, by
#: `CoQ_check/cq04_coq_ms2_peaklist.py`: median over the qualifying spectra, base peak 999, kept
#: at >= 70 % presence. Masses are computed from these formulae, not from the measurement.
#:
#: Two families, and both are generic to the class rather than to the homologue:
#:   * the QUINONE HEAD — 197.0808 C10H13O4+ and its losses, 182.0574 (-CH3.), 169.0859 (-CO),
#:     154.0624, 151.0754. 197.0808 is the base peak of every qualifying spectrum, 100 % presence.
#:   * the POLYPRENYL CHAIN — the CnH(2n-1)+ / CnH(2n-3)+ hydrocarbon ladder, the same series a
#:     sterol gives and the reason these entries carry no chain annotation (see below).
#:
#: ⚠ `184.0733` was measured at median 395 in 7/7 CoQ9 spectra and is NOT written: it is the
#: phosphocholine head group of a co-isolated PC/SM, and writing it would make this entry match
#: every PC spectrum in the run.
SPECTRA: dict[tuple[str, str], list[tuple[str, int]]] = {
    # 12 qualifying spectra. The dominant ion of the class.
    ("CoQ10", "[M+NH4]+"): [
        ("C10H13O4", 999), ("C6H9", 525), ("C7H11", 401), ("C8H13", 248), ("C8H11", 221),
        ("C9H13", 208), ("C7H9", 193), ("C5H7", 144), ("C9H15", 123), ("C11H17", 110),
    ],
    # 10 qualifying spectra, and the richest pattern — the protonated species fragments further,
    # which is why it carries the head-group losses and the long end of the prenyl ladder.
    ("CoQ10", "[M+H]+"): [
        ("C10H13O4", 999), ("C7H11", 192), ("C6H9", 170), ("C8H13", 148), ("C9H13", 91),
        ("C8H11", 88), ("C11H17", 86), ("C10H15", 85), ("C9H15", 85), ("C7H9", 81),
        ("C5H9", 60), ("C12H19", 56), ("C11H15", 56), ("C12H17", 56), ("C10H13", 53),
        ("C9H11", 51), ("C9H11O2", 48), ("C9H13O3", 48), ("C5H7", 47), ("C13H19", 43),
        ("C6H11", 36), ("C12H15", 33), ("C13H21", 33), ("C8H9", 32), ("C10H17", 32),
        ("C13H17O4", 32), ("C14H21", 31), ("C8H10O3", 31), ("C6H7", 30), ("C14H19", 29),
        ("C11H13", 26), ("C7H7", 26),
    ],
    # 11 qualifying spectra. ⚠ Read the module docstring before using this entry: it is the one
    # that can be fed by the ubiquinone M+2 isotope, and its spectrum cannot tell the difference.
    # Note what is NOT here: the reduced head group 199.0965 (C10H15O4+) was looked for in every
    # one of these spectra and observed in NONE. The quinol's head oxidises in the source, so the
    # spectrum is the quinone's. An authentic ubiquinol standard is what would settle it.
    ("CoQ10H2", "[M+NH4]+"): [
        ("C10H13O4", 999), ("C6H9", 256), ("C7H11", 166), ("C8H13", 112), ("C9H11O2", 88),
        ("C9H13", 88), ("C9H13O3", 87), ("C8H11", 83), ("C7H9", 71), ("C5H7", 64),
        ("C9H15", 61), ("C5H9", 61),
    ],
    # 7 qualifying spectra. Same pattern as CoQ10, as it must be: one isoprene unit shorter
    # changes the precursor and nothing about the head group.
    ("CoQ9", "[M+NH4]+"): [
        ("C10H13O4", 999), ("C6H9", 550), ("C7H11", 366), ("C8H13", 259), ("C9H13", 197),
        ("C5H9", 188), ("C8H11", 187), ("C11H17", 168), ("C10H15", 138),
    ],
}

#: What the predicted homologues inherit, and from which measured entry. A quinol takes the
#: quinone's spectrum because that is what the quinol was measured to give (see above).
PREDICT_FROM = {"quinone": ("CoQ10", "[M+NH4]+"), "quinol": ("CoQ10H2", "[M+NH4]+")}

#: Measured retention, an in-house Q Exactive Plus study, 10 injections (cq01) and 3 (cq03 profiles).
#: ⚠ The reason this file exists at all: ubiquinol and the ubiquinone M+2 isotope are 10 ppm
#: apart and separate on retention alone.
RETENTION = {
    ("CoQ10", "[M+NH4]+"): (21.07, 0.15), ("CoQ10", "[M+H]+"): (21.07, 0.15),
    ("CoQ10", "[M+Na]+"): (21.05, 0.15),
    ("CoQ10H2", "[M+NH4]+"): (19.50, 0.25), ("CoQ10H2", "[M+Na]+"): (19.54, 0.25),
    ("CoQ9", "[M+NH4]+"): (19.10, 0.25),
}


def homologue(base: str, n: int) -> tuple[str, float]:
    """`CoQ10` shifted to n isoprene units: the name and the neutral mass."""
    formula, units, _ = QUINONES[base]
    reduced = base.endswith("H2")
    stem = "CoQ" + str(n) + ("H2" if reduced else "")
    return stem, neutral_mass(formula) + (n - units) * ISOPRENE


def entry(name: str, adduct: str, neutral: float, fragments: list[tuple[str, int]],
          formula: str, note: str, confirmed: bool) -> list[str]:
    """One .msp block, in the LipiDex 1 dialect `msp.parse_msp` reads."""
    precursor = neutral + ADDUCTS[adduct]
    peaks = [(cation_mass(f), i, f) for f, i in fragments]
    peaks = [(mz, i, f) for mz, i, f in peaks
             if precursor - mz > PRECURSOR_CUTOFF and mz >= MIN_MS2_MASS]
    peaks.sort()
    rt, window = RETENTION.get((name, adduct), (None, None))

    lines = [f"Name: {name} {adduct};",
             f"MW: {neutral:.4f}",
             f"PRECURSORMZ: {precursor:.4f}"]
    # ⚠ Keys beyond the five LipiDex writes are ignored by every reader here, which is what makes
    # them safe to add — and `RT=` is the only place the retention constraint travels with the
    # entry itself. Do not write a literal "Name:" or "Num Peaks:" into this line; the parser
    # looks for those substrings anywhere on a line.
    comment = (f"Comment: Name={name} {adduct} Mass={precursor:.4f} Formula={formula} "
               f"OptimalPolarity=true Type=LipiDex Category=Prenol "
               f"Evidence={'measured' if confirmed else 'predicted'}")
    if rt is not None:
        comment += f" RT={rt:.2f} RTWindow={window:.2f}"
    if note:
        comment += f" Note={note}"
    lines.append(comment)
    # No moiety annotations: a quinone has no acyl chain, so there is nothing for the purity
    # calculation to count and it returns 0 — the same as `Sterols_InSilico`. The name is a
    # defined structure, not a sum composition, so Purity 0 costs no information here.
    lines.append(f"Num Peaks: {len(peaks)}")
    for mz, intensity, frag in peaks:
        lines.append(f'{mz:.4f} {intensity} "{frag}_Fragment_[]"')
    lines.append("")
    return lines


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--predicted-homologues", action="store_true",
                    help="also write CoQ6-CoQ11 and their quinols, [M+NH4]+, from the measured "
                         "CoQ10 pattern. PREDICTIONS: only CoQ9, CoQ10 and ubiquinol-10 have been "
                         "confirmed here, and every unconfirmed entry is another chance to take a "
                         "match from something correct.")
    ap.add_argument("--anchors", help="also write the measured retention table to this path")
    args = ap.parse_args()

    lines: list[str] = []
    written: list[tuple[str, str, float, bool]] = []
    for (name, adduct), fragments in SPECTRA.items():
        formula, _, note = QUINONES[name]
        neutral = neutral_mass(formula)
        lines += entry(name, adduct, neutral, fragments, formula, note.replace(" ", "_"), True)
        written.append((name, adduct, neutral + ADDUCTS[adduct], True))

    if args.predicted_homologues:
        for n in range(6, 12):
            for base, kind in (("CoQ10", "quinone"), ("CoQ10H2", "quinol")):
                stem, neutral = homologue(base, n)
                if (stem, "[M+NH4]+") in SPECTRA:
                    continue                       # measured already, do not overwrite it
                src = SPECTRA[PREDICT_FROM[kind]]
                units = 2 * n + (2 if kind == "quinol" else 0)
                formula = f"C{9 + 5 * n}H{10 + 8 * n + (2 if kind == 'quinol' else 0)}O4"
                assert abs(neutral_mass(formula) - neutral) < 1e-6, (formula, neutral, units)
                lines += entry(stem, "[M+NH4]+", neutral, src, formula, "predicted_homologue",
                               False)
                written.append((stem, "[M+NH4]+", neutral + ADDUCTS["[M+NH4]+"], False))

    Path(args.out).write_text("\n".join(lines))
    print(f"{len(written)} entries written to {args.out}")
    for name, adduct, mz, confirmed in written:
        print(f"  {name:10} {adduct:10} {mz:9.4f}  {'measured' if confirmed else 'PREDICTED'}")

    iso = (neutral_mass(QUINONES['CoQ10'][0]) + ADDUCTS['[M+NH4]+']) + 2 * CARBON_13
    quinol = neutral_mass(QUINONES['CoQ10H2'][0]) + ADDUCTS['[M+NH4]+']
    print(f"\n  ⚠ ubiquinol [M+NH4]+ {quinol:.4f} vs ubiquinone M+2 isotope {iso:.4f}: "
          f"{abs(quinol - iso) * 1000:.1f} mDa, {abs(quinol - iso) / quinol * 1e6:.1f} ppm — "
          f"inside MS1_TOL (0.01 Da). Separable by RETENTION only; run "
          f"scripts/check_ubiquinol_rt.py on any study reporting ubiquinol.")

    if args.anchors:
        rows = ["name,adduct,mz,retention,window,evidence"]
        for (name, adduct), (rt, window) in sorted(RETENTION.items()):
            formula = QUINONES[name][0]
            rows.append(f"{name},{adduct},{neutral_mass(formula) + ADDUCTS[adduct]:.4f},"
                        f"{rt:.2f},{window:.2f},measured")
        Path(args.anchors).write_text("\n".join(rows) + "\n")
        print(f"  retention anchors -> {args.anchors}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
