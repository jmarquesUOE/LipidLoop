"""End-to-end check against `Final_Results.csv`-grade output LipiDex itself produced.

The .mgf here is the same file LipiDex read, so nothing but this code stands between input and
result — no feature detection, no alignment. Skipped when the datastore is not mounted, since
the reference data lives there and is far too large to vendor.

Slow (~15 s): loads 175k library spectra and searches ~29k sample spectra. Run it with
`pytest -m integration`, or leave it in the default run and take the wait.
"""
import os
import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.purity import read_fatty_acids                  # noqa: E402
from lipidloop.search import iter_results, load_libraries      # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(os.environ.get("LIPIDLOOP_REFERENCE_DIR", "reference/Test_CD-Lipidex"))
# Load order breaks ties between entries that score identically, so it has to match the order
# the libraries were ticked in the GUI for the reference run — LipidBlast first. Deliberately
# NOT `DEFAULT_LIBRARIES`: this test reproduces one specific LipiDex run, which was made without
# the ganglioside library, so searching a different set would be testing something else.
LIBRARIES = [ROOT / "data/libraries/LipidBlast_Formic.msp",
             ROOT / "data/libraries/LipiDex_HCD_Formic.msp",
             ROOT / "data/libraries/LipiDex_HCD_Hydroxy.msp"]
FATTY_ACIDS = ROOT / "data/lipidex_src/FattyAcids.csv"

pytestmark = [
    pytest.mark.integration,
    pytest.mark.skipif(not all(p.exists() for p in LIBRARIES + [FATTY_ACIDS]),
                       reason="libraries not present under data/"),
    pytest.mark.skipif(not DATA.exists(), reason="datastore not mounted"),
]


@pytest.fixture(scope="module")
def index():
    return load_libraries(LIBRARIES)


def run(polarity, name, index):
    reference = {int(r["MS2 ID"]): r
                 for r in csv.DictReader((DATA / polarity / f"{name}_Results.csv").open(newline=""))}
    fa_db = read_fatty_acids(FATTY_ACIDS)
    ours = {}
    for rows, hits in iter_results(DATA / polarity / f"{name}.mgf", index, fa_db):
        top = hits[0].dot
        ours[rows[0]["MS2 ID"]] = (
            rows[0],
            {(h.library_spectrum.name.strip(), h.library_spectrum.library) for h in hits
             if abs(h.dot - top) <= 1e-9 * max(abs(top), 1.0)})
    return reference, ours


@pytest.mark.parametrize("polarity,name", [("Neg", "QC_01"), ("Pos", "QC_01")])
def test_reproduces_lipidex(polarity, name, index):
    reference, ours = run(polarity, name, index)

    assert set(ours) == set(reference), "different spectra identified"

    for k, (row, tied) in ours.items():
        ref = reference[k]
        assert row["Dot Product"] == int(ref["Dot Product"]), f"dot product, MS2 {k}"
        assert row["Reverse Dot Product"] == int(ref["Reverse Dot Product"]), f"reverse, MS2 {k}"
        assert row["Delta m/z"] == pytest.approx(float(ref["Delta m/z"]), abs=1e-4), f"delta, MS2 {k}"

        # An exact tie is not a decision the algorithm makes: where two library entries score
        # equal, which one comes back is down to the last floating-point bit, and Java and
        # CPython round `pow` differently. Accept any member of the tied set, matched on name
        # *and* library — a tie can span two libraries, and the same lipid appears in both.
        ours_row = (row["Identification"].strip(), row["Library"])
        ref_row = (ref["Identification"].strip(), ref["Library"])
        assert ours_row == ref_row or ref_row in tied, f"identification, MS2 {k}: {ours_row}"

        # Purity is only comparable where the identification agrees — across a tie it would be
        # the purity of a different lipid.
        if ours_row == ref_row:
            assert row["Purity"] == int(ref["Purity"]), f"purity, MS2 {k}"
            assert row["Spectral Components"] == ref["Spectral Components"], f"components, MS2 {k}"
            assert row["Potential Fragments"] == ref["Potential Fragments"], f"fragments, MS2 {k}"
