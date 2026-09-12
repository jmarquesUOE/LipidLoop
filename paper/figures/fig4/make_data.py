"""Figure 4 data: the public corpus (A) and the two-pass calibration on the three studies that
needed it (B).

(A) is parsed from `manuscript/validation/RESULTS_TABLE.md`, the corpus table the manuscript
quotes, so the figure cannot drift from it. (B) is Supporting Information Table S7, quoted.
"""
from __future__ import annotations

import csv
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
DATA = HERE / "data"
DATA.mkdir(exist_ok=True)

VENDOR = {
    "MTBKS222_Waters": "Waters", "MTBKS222_Thermo": "Thermo", "MTBKS222_Agilent": "Agilent",
    "MTBKS222_IMS_Normal": "Sciex", "MTBKS222_IMS_HighMass": "Sciex", "MTBKS222_BrukerBAF": "Bruker",
    "MTBLS5163": "Bruker", "ST003077": "Thermo", "ST003514": "Agilent", "ST004797": "Sciex",
    "ST004797_oxtg": "Sciex", "MSV000095868": "Agilent", "ST000991_DDA": "Sciex",
    "ST002705": "Bruker", "ST003052": "Thermo",
}
INSTRUMENT = {
    "MTBKS222_Waters": "Xevo G2 QTOF", "MTBKS222_Thermo": "Q Exactive Plus", "MTBKS222_Agilent": "6546 QTOF",
    "MTBKS222_IMS_Normal": "TripleTOF 6600", "MTBKS222_IMS_HighMass": "TripleTOF 6600 (high mass)",
    "MTBKS222_BrukerBAF": "timsTOF Pro", "MTBLS5163": "maXis II", "ST003077": "Exploris 480",
    "ST003514": "6545 QTOF", "ST004797": "ZenoTOF 7600", "ST004797_oxtg": "ZenoTOF 7600 (oxTG)",
    "MSV000095868": "6530A QTOF", "ST000991_DDA": "X500R QTOF", "ST002705": "micrOTOF-Q II",
    "ST003052": "Q Exactive HF",
}
NIST = {"MTBKS222_Waters", "MTBKS222_Thermo", "MTBKS222_Agilent", "MTBKS222_IMS_Normal",
        "MTBKS222_IMS_HighMass", "MTBKS222_BrukerBAF", "ST003514"}


def corpus() -> list[dict]:
    text = (ROOT / "manuscript" / "validation" / "RESULTS_TABLE.md").read_text()
    rows = []
    for line in text.splitlines():
        if not line.startswith("| ") or "| dataset |" in line or line.startswith("|---"):
            continue
        cells = [c.strip() for c in line.strip("|").split("|")]
        if len(cells) < 11 or cells[0] not in VENDOR:
            continue
        ds = cells[0]
        agree = re.sub(r"[*%]", "", cells[10]).strip()
        agree = agree if re.fullmatch(r"[\d.]+", agree) else ""
        rows.append({"dataset": ds, "vendor": VENDOR[ds], "instrument": INSTRUMENT[ds],
                     "nist": int(ds in NIST), "injections": int(cells[2].split("(")[1].rstrip(")")),
                     "id_rows": int(cells[4].replace(",", "")), "species": int(cells[5].replace(",", "")),
                     "decoy_fdr": float(cells[7].rstrip("%")),
                     "deposit_species": int(cells[8].replace(",", "")) if cells[8] not in ("—", "") else "",
                     "agreement": float(agree) if agree and agree != "—" else ""})
    assert len(rows) == 15, len(rows)
    return rows


def calibration() -> list[dict]:
    # Supporting Information Table S7.
    return [
        {"study": "MTBLS5163", "instrument": "maXis II", "polarity": "Pos", "drift_ppm": 17.1, "before": 232, "after": 258, "decoy_before": 2, "decoy_after": 0},
        {"study": "MTBKS222_Waters", "instrument": "Xevo G2", "polarity": "Pos", "drift_ppm": 9.3, "before": 199, "after": 211, "decoy_before": 0, "decoy_after": 0},
        {"study": "MTBKS222_Waters", "instrument": "Xevo G2", "polarity": "Neg", "drift_ppm": 9.3, "before": 41, "after": 50, "decoy_before": 0, "decoy_after": 0},
        {"study": "ST004797", "instrument": "ZenoTOF 7600", "polarity": "Pos", "drift_ppm": 2.6, "before": 551, "after": 561, "decoy_before": 1, "decoy_after": 0},
        {"study": "ST004797", "instrument": "ZenoTOF 7600", "polarity": "Neg", "drift_ppm": 2.6, "before": 111, "after": 117, "decoy_before": 0, "decoy_after": 0},
    ]


def write(name, rows):
    with (DATA / name).open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    write("corpus.csv", corpus())
    write("calibration.csv", calibration())
    print("written to", DATA)
