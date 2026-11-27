"""Centroid a profile study once, to disk, instead of on every load.

ST004797 is 199.5 GB of profile mzML — ~25,000 sampling points per MS1 scan, 781 MB for a single
blank. The pipeline already centroids it on load (`calibrate.load_centroided`), so all of that is
read, peak-picked in memory and discarded on every run: about six hours per pass, which is why the
study sits last in every batch and why a two-pass mass calibration on it is unattractive.

Doing it once and keeping the result changes nothing about the numbers — the same
`PeakPickerHiRes`, the same settings — and makes repeat runs affordable.

⚠ **The profile originals are never touched.** They are read-only on the shared store and they are
what let the isolation-window, duty-cycle and mass-accuracy work happen at all; centroiding is
one-way. Output goes to a separate directory and the study is re-staged to point at it, so
reverting means re-pointing the symlinks.

    python scripts/precentroid.py <study> --out <dir> [--workers 4]
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def centroid_one(source: Path, target: Path) -> tuple[int, int]:
    """(bytes in, bytes out). Skips work already done."""
    from pyopenms import MzMLFile
    from lipidloop.calibrate import load_centroided

    if target.exists() and target.stat().st_size > 0:
        return source.stat().st_size, target.stat().st_size
    experiment = load_centroided(str(source))
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_suffix(".partial")
    # ⚠ Write to a temporary name and rename. A run interrupted midway otherwise leaves a
    # truncated mzML that looks complete, and the next pass skips it — silently analysing half a
    # file. The rename is atomic on the same filesystem.
    MzMLFile().store(str(tmp), experiment)
    tmp.rename(target)
    return source.stat().st_size, target.stat().st_size


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("study")
    ap.add_argument("--out", required=True)
    ap.add_argument("--min-free-gb", type=int, default=20)
    args = ap.parse_args()

    study, out = Path(args.study), Path(args.out)
    files = sorted(f for p in ("Pos", "Neg") for f in (study / p).glob("*.mzML"))
    if not files:
        print("  nothing to do")
        return 1
    print(f"  {len(files)} files")

    import shutil
    done = bytes_in = bytes_out = 0
    for f in files:
        free = shutil.disk_usage(out.parent if out.parent.exists() else "/").free / 1e9
        if free < args.min_free_gb:
            print(f"  STOP  only {free:.0f} GB free, below the {args.min_free_gb} GB floor")
            break
        polarity = f.parent.name
        try:
            a, b = centroid_one(f, out / polarity / f.name)
        except Exception as exc:                      # noqa: BLE001 — one bad file must not stop the set
            print(f"  FAIL  {f.name}: {type(exc).__name__}: {exc}")
            continue
        bytes_in += a; bytes_out += b; done += 1
        if done % 10 == 0 or done == len(files):
            print(f"    {done}/{len(files)}  {bytes_in/1e9:.1f} GB -> {bytes_out/1e9:.1f} GB "
                  f"({bytes_in/max(bytes_out,1):.1f}x smaller)", flush=True)
    print(f"  done {done}/{len(files)}: {bytes_in/1e9:.1f} GB -> {bytes_out/1e9:.1f} GB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
