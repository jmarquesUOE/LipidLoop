"""Did MS2 scans reach the AGC target, or did the clock run out?

An Orbitrap fills until it reaches the AGC target *or* the maximum injection time, whichever comes
first. A scan whose injection time sits exactly at the method maximum stopped because it ran out of
time, not because it had enough ions — the spectrum is ion-starved, and its quality is set by
however many ions happened to arrive in that window.

Reads injection time per MS2 scan from the mzML, and the method settings from the `.raw` header
(UTF-16, no ThermoRawFileParser needed) so the ceiling is read from the method rather than inferred
from the data.

Deliberately does *not* compare total ion current against the AGC target: TIC here is the sum of
centroided peak intensities in detector units, not a count of ions, and the two are not the same
quantity. The injection time alone answers the question.

Produces the numbers behind manuscript §5.3.
"""
from __future__ import annotations

import argparse
import glob
import os
import re
import xml.etree.ElementTree as ET

import numpy as np

NS = "{http://psi.hupo.org/ms/mzml}"
WANT = ("ms level", "ion injection time", "total ion current", "scan start time")


def method_ms2(raw_path: str, probe: int = 6_000_000) -> dict:
    """MS2 branch of the acquisition method, from the raw header."""
    text = open(raw_path, "rb").read(probe).decode("utf-16-le", errors="ignore").replace("\x00", "")
    text = re.sub(r"[^\x20-\x7e]+", " | ", text)
    start = re.search(r"MSn Level = 2", text)
    if not start:
        return {}
    segment = text[start.start():start.start() + 900]
    keys = ("Maximum Injection Time (ms)", "Absolute AGC Target", "Normalized AGC Target (%)",
            "Orbitrap Resolution", "Isolation Window (m/z)", "HCD Collision Energy/Energies (%)")
    out = {}
    for key in keys:
        found = re.search(re.escape(key) + r" = ([^|]+)", segment)
        if found:
            out[key] = found.group(1).strip()
    return out


def scans(mzml_path: str):
    """(ms level, injection time, TIC, retention) per spectrum.

    Only spectrum elements are cleared. Clearing every element as it ends destroys a spectrum's
    own children before its end event fires, which silently yields nothing at all.
    """
    for _, element in ET.iterparse(mzml_path, events=("end",)):
        if element.tag != NS + "spectrum":
            continue
        got = {}
        for param in element.iter(NS + "cvParam"):
            name = param.get("name")
            if name in WANT:
                got[name] = float(param.get("value"))
        if "ms level" in got:
            yield (got["ms level"], got.get("ion injection time"),
                   got.get("total ion current", 0.0), got.get("scan start time", 0.0))
        element.clear()


def audit(mzml_dir: str, raw_dir: str, pattern: str, tolerance: float = 1e-9):
    files = sorted(glob.glob(os.path.join(mzml_dir, pattern)))
    if not files:
        raise SystemExit(f"no mzML matching {pattern!r} in {mzml_dir}")

    raws = sorted(glob.glob(os.path.join(raw_dir, pattern.replace(".mzML", ".raw"))))
    if raws:
        method = method_ms2(raws[0])
        print("  method (MS2): " + ", ".join(f"{k} = {v}" for k, v in method.items()))
        ceiling = float(method.get("Maximum Injection Time (ms)", "nan"))
    else:
        ceiling = float("nan")

    times, tics, levels, rts = [], [], [], []
    for path in files:
        for level, injection, tic, rt in scans(path):
            levels.append(level)
            rts.append(rt)
            if level == 2 and injection is not None:
                times.append(injection)
                tics.append(tic)
    times, tics = np.array(times), np.array(tics)
    if not np.isfinite(ceiling):
        ceiling = times.max()
        print(f"  no raw header found; taking the observed maximum {ceiling:.1f} ms as the ceiling")

    at = times >= ceiling - tolerance
    print(f"  {len(files)} files, {len(times):,} MS2 scans")
    print(f"  injection time: min {times.min():.1f}  median {np.median(times):.1f}  "
          f"max {times.max():.1f} ms   (method maximum {ceiling:.0f} ms)")
    print(f"  AT THE CEILING: {100 * at.mean():.1f}%   filled early: {100 * (~at).mean():.1f}%")
    if at.any() and (~at).any():
        print(f"  median TIC — at ceiling {np.median(tics[at]):,.0f}   "
              f"filled early {np.median(tics[~at]):,.0f}  "
              f"(ratio {np.median(tics[~at]) / max(np.median(tics[at]), 1):.0f}x)")

    levels, rts = np.array(levels), np.array(rts)
    ms1 = np.where(levels == 1)[0]
    if len(ms1) > 2:
        per = np.diff(ms1) - 1
        cycle = np.diff(rts[ms1]) * 60
        share = np.median(per) * ceiling / 1000 / max(np.median(cycle), 1e-9)
        print(f"  cycle: {np.median(per):.0f} MS2 per MS1, {np.median(cycle):.2f} s; "
              f"MS2 ion accumulation is {100 * share:.0f}% of it")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mzml", required=True)
    parser.add_argument("--raw", default="", help="directory of matching .raw, for method settings")
    parser.add_argument("--pattern", default="*.mzML")
    args = parser.parse_args()
    audit(args.mzml, args.raw or args.mzml, args.pattern)


if __name__ == "__main__":
    main()
