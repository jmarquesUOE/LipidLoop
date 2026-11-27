"""Build CH2 decoys for libraries that do not annotate their chain fragments.

⚠ Four libraries carry NO chain annotations, so `make_decoy_library.py` cannot touch them:

    AcylHexCer_Positive           10,800
    AcylSM_Positive               10,416
    UltraLongHexCer_EOS_Positive   3,759
    UltraLongCer_EOS_Positive      3,179

28,283 entries — 13.5% of the formate positive pool — with no decoy of any kind, so their
false-discovery rate is uncountable rather than low.

They are not uninformative. Every peak list in AcylHexCer is distinct (10,800 of 10,800), so a
decoy built from them competes. What is missing is the annotation, not the information, and the
annotation can be recovered by differencing homologues: two entries one CH2 apart in composition
must have their chain peaks 14.0157 apart and their head-group peaks identical.

    python scripts/build_unannotated_decoys.py --shift 3
"""
from __future__ import annotations
import argparse, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from make_decoy_library import CH2, read_entries          # noqa: E402

LIB = ROOT / "data/libraries"
SOURCES = ["AcylHexCer_Positive", "AcylSM_Positive",
           "UltraLongHexCer_EOS_Positive", "UltraLongCer_EOS_Positive"]
NUM = re.compile(r"(\d+):(\d+)")


def entry_name(lines):
    for l in lines:
        if l.lower().startswith("name:"):
            return l.split(":", 1)[1].strip()
    return ""


def peaks_of(lines):
    out = []
    for l in lines:
        t = l.strip()
        if t and t[0].isdigit():
            try:
                out.append(float(t.split()[0]))
            except (ValueError, IndexError):
                pass
    return out


def build(stem: str, shift: int, tol: float = 0.01):
    entries = [list(e) for e in read_entries(LIB / f"{stem}.msp")]
    meta = [(entry_name(e), NUM.findall(entry_name(e)), peaks_of(e), e) for e in entries]
    by_len = {}
    for m in meta:
        by_len.setdefault(len(m[1]), []).append(m)

    written = skipped = 0
    out_lines = []
    for name, comp, peaks, lines in meta:
        moved = set()
        for oname, ocomp, opeaks, _ in by_len.get(len(comp), ()):
            if oname == name or not opeaks or not peaks:
                continue
            d = [(int(a[0]) - int(b[0]), int(a[1]) - int(b[1])) for a, b in zip(comp, ocomp)]
            if sorted(d) != sorted([(0, 0)] * (len(d) - 1) + [(1, 0)]):
                continue
            for p in peaks:
                if any(abs(p - q - CH2) < tol for q in opeaks) and \
                   not any(abs(p - q) < tol for q in opeaks):
                    moved.add(round(p, 4))
            if moved:
                break
        if not moved:
            skipped += 1          # no homologue: skipped rather than guessed at
            continue
        new = []
        for l in lines:
            t = l.strip()
            # ⚠ read_entries strips newlines. Every branch must put one back, or the whole header
            # collapses onto one line and the file is unparseable.
            if t.lower().startswith("name:"):
                # ⚠ Case matters: these libraries write `NAME:`, not `Name:`. A decoy that does not
                # carry the DECOY_ prefix is counted as a TARGET identification — the FDR does not
                # just lose accuracy, it inverts.
                head, _, rest = t.partition(":")
                new.append(f"{head}: DECOY_{rest.strip()}\n")
            elif t and t[0].isdigit():
                parts = t.split(None, 1)
                mz = float(parts[0])
                # precursor and head-group peaks untouched: the decoy must still look like its
                # class and compete for the same spectra.
                if round(mz, 4) in moved:
                    mz += shift * CH2
                new.append(f"{mz:.4f} {parts[1] if len(parts) > 1 else ''}".rstrip() + "\n")
            else:
                new.append(l if l.endswith("\n") else l + "\n")
        out_lines.extend(new)
        out_lines.append("\n")
        written += 1
    (LIB / f"DECOY1_{stem}.msp").write_text("".join(out_lines))
    return written, skipped


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--shift", type=int, default=3)
    a = ap.parse_args()
    total_w = total_s = 0
    for stem in SOURCES:
        if not (LIB / f"{stem}.msp").exists():
            print(f"  {stem:<32} MISSING"); continue
        w, s = build(stem, a.shift)
        total_w += w; total_s += s
        print(f"  {stem:<32} {w:>6} decoys, {s:>6} skipped (no homologue)", flush=True)
    print(f"\n  {total_w} decoys written, {total_s} entries had no homologue to difference against")
