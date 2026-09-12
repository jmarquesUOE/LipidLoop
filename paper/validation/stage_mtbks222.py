"""Split MTBKS222's converted runs into one study per vendor arm, polarity read from the file.

⚠ **One study per vendor, never all five together.** MTBKS222 is the same NIST SRM 1950 plasma on
five instruments. Pooling them into a single run would ask the feature grouper to align peaks
across different chromatography and different mass analysers, which is meaningless — and each
vendor publishes its own MAF, so they have to be scored separately anyway.

⚠ **Polarity comes from the mzML, not the filename.** The five vendors each name it differently —
`Pooled_Neg_`, `QExPlus_LIPn_`, `..._Neg.wiff`, `190918_101nn` — and a bare letter inside a date
token defeats any suffix rule. The cvParam is authoritative and costs one pass of the header.
"""
from __future__ import annotations

import re
import shutil
import sys
from pathlib import Path

STAGED = Path("/home/jair/validation_staged")
SOURCE = STAGED / "MTBKS222/_converted"
RAW = Path("/mnt/datastore/Jair/claudecode/Lipidomics_automation/Validation_Datasets/MTBKS222/raw")


def polarity(path: Path, limit: int = 4000) -> str:
    """`Pos` / `Neg` from the first scan's cvParam, or "" if neither appears."""
    try:
        with path.open(errors="replace") as fh:
            for i, line in enumerate(fh):
                if "MS:1000130" in line or 'name="positive scan"' in line:
                    return "Pos"
                if "MS:1000129" in line or 'name="negative scan"' in line:
                    return "Neg"
                if i > limit:
                    break
    except OSError:
        pass
    return ""


def vendor_of(stem: str) -> str:
    """Which acquisition this mzML came from, by matching it back to the raw directory."""
    for suffix in (".d", ".raw", ".wiff"):
        p = RAW / f"{stem}{suffix}"
        _ = stem
        if not p.exists():
            continue
        if p.is_file():
            if suffix != ".wiff":
                return "Thermo"
            # ⚠ Three different SUBMISSIONS share the .wiff extension, each with its own MAF and
            # its own gradient. Lumping them as "Sciex" scored the Tsugawa files against the UC
            # Davis list and produced a 4.22 min retention offset with 0% agreement — which reads
            # as broken chromatography and is really two labs' data compared to each other.
            if "HighMass" in stem:
                return "IMS_HighMass"
            if "Normal" in stem:
                return "IMS_Normal"
            return "UCDavis"
        names = {c.name.lower() for c in p.iterdir()}
        if "analysis.tdf" in names:
            return "BrukerTDF"
        if "analysis.baf" in names:
            return "BrukerBAF"
        if "acqdata" in names:
            return "Agilent"
        if any(n.startswith("_func") for n in names):
            return "Waters"
    return ""


def main() -> int:
    files = sorted(SOURCE.glob("*.mzML"))
    if not files:
        print("  nothing converted yet")
        return 1
    counts: dict[tuple[str, str], int] = {}
    unplaced = []
    for f in files:
        vendor, pol = vendor_of(f.stem), polarity(f)
        if not vendor or not pol:
            unplaced.append((f.name, vendor or "?", pol or "?"))
            continue
        target = STAGED / f"MTBKS222_{vendor}" / pol
        target.mkdir(parents=True, exist_ok=True)
        link = target / f.name
        if not link.exists():
            link.symlink_to(f.resolve())
        counts[(vendor, pol)] = counts.get((vendor, pol), 0) + 1

    for (vendor, pol), n in sorted(counts.items()):
        print(f"  MTBKS222_{vendor:<10} {pol}  {n:3d} files")
    if unplaced:
        # ⚠ Never silent. A file that lands nowhere is invisible in every downstream count.
        print(f"  ⚠ {len(unplaced)} could not be placed (vendor/polarity):")
        for name, v, p in unplaced[:10]:
            print(f"      {name}  vendor={v} polarity={p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
