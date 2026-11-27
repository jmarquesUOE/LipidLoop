"""Get everything into `.mzML`, converting only what needs converting.

Thermo `.raw` goes through ThermoRawFileParser, which ships a self-contained Linux build — no mono
and no .NET runtime to install, which is the whole reason this is a subprocess rather than a
library call.

**Input already in mzML passes straight through.** That sounds obvious and was not the case: the
pipeline sent every input to the Thermo converter with no format check, so a study of perfectly
good mzML failed — and failed on a *missing Thermo converter it never needed*, because the
converter was located before anyone asked whether there was anything to convert. Public data,
collaborators' data and instruments that export mzML directly were all unusable.

Conversion is skipped when an `.mzML` already exists and is newer than its `.raw`. Converting is
the slowest step per file that has no algorithmic content, and the output is deterministic.
"""
from __future__ import annotations

import os
import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path

INDEXED_MZML = "2"   # -f: 0 MGF, 1 mzML, 2 indexed mzML

# Suffixes that need no conversion. Case-insensitive: writers disagree (`.mzML`, `.mzml`), and a
# case-sensitive check would silently send a converted file back to the converter.
READY = {".mzml"}


def needs_conversion(path: str | Path) -> bool:
    """True when this input has to go through a vendor reader before the pipeline can read it."""
    return Path(path).suffix.lower() not in READY


def detect_format(path: str | Path) -> str:
    """What this path actually is: `mzml`, `thermo`, `agilent`, `bruker`, `waters`, `sciex`, ``.

    ⚠ Dispatch on WHAT THE PATH IS, never on its suffix alone. Waters `.raw` is a DIRECTORY and
    Thermo `.raw` is a FILE — the same suffix, two vendors, two readers. Handing a Waters
    directory to ThermoRawFileParser fails in a way that reads as a corrupt file.
    """
    p = Path(path)
    suffix = p.suffix.lower()
    if suffix == ".mzml":
        return "mzml"
    if p.is_dir():
        if (p / "analysis.tdf").exists() or (p / "analysis.baf").exists():
            return "bruker"
        if (p / "AcqData").is_dir():
            return "agilent"
        if any(p.glob("_FUNC*.DAT")):
            return "waters"
        return ""
    if suffix == ".raw":
        return "thermo"
    if suffix in (".wiff", ".wiff2"):
        return "sciex"
    return ""


class MsConvert:
    """ProteoWizard's msconvert, for the vendor formats Thermo's own reader cannot open.

    Agilent, Waters and Sciex readers are Windows libraries; they run under wine inside the
    ProteoWizard container. Bruker is not routed here — it ships Linux readers and is handled
    natively, which is both faster and one less dependency.

    ⚠ The container must NOT be left in udocker's PRoot execution mode. PRoot patches ELF binaries
    in place, which corrupts wine's 32-bit loader and produces `not enough space for ELF headers`
    on every run — including afterwards under other modes, because the patch persists in the
    container. Create the container fresh and set `--execmode=R2` before first use.

    Requires unprivileged user namespaces (`kernel.apparmor_restrict_unprivileged_userns=0`).
    """

    def __init__(self, udocker: str = "udocker", container: str = "pwiz2"):
        self.udocker = udocker
        self.container = container

    def available(self) -> bool:
        return bool(shutil.which(self.udocker)) and Path(
            "/proc/sys/kernel/apparmor_restrict_unprivileged_userns").read_text().strip() == "0"

    def convert(self, source: str | Path, out_dir: str | Path) -> Path:
        source, out_dir = Path(source).resolve(), Path(out_dir).resolve()
        out = out_dir / f"{source.stem}.mzML"
        if out.exists() and out.stat().st_size > 0:
            return out
        out_dir.mkdir(parents=True, exist_ok=True)
        script = (
            "export WINEPREFIX=/wineprefix64 WINEDEBUG=-all HOME=/tmp; "
            f"wine /wineprefix64/drive_c/pwiz/msconvert.exe '/in/{source.name}' "
            "--mzML --64 --zlib -o /out"
        )
        result = subprocess.run(
            [self.udocker, "run",
             f"--volume={source.parent}:/in", f"--volume={out_dir}:/out",
             "--entrypoint=/bin/sh", self.container, "-c", script],
            capture_output=True, text=True)
        if not out.exists():
            raise ConversionError(
                f"msconvert produced nothing for {source.name}:\n"
                f"{result.stdout[-800:]}\n{result.stderr[-800:]}")
        return out


class ConversionError(RuntimeError):
    pass


@dataclass(slots=True)
class Converter:
    """Locates ThermoRawFileParser and runs it."""

    executable: Path

    #: Where the parser is unpacked in a checkout. The directory and the executable inside it
    #: share a name, which is how the release archive is laid out.
    BUNDLED = Path(__file__).resolve().parents[2] / "data/tools/ThermoRawFileParser"

    @classmethod
    def find(cls, root: Path | None = None) -> "Converter":
        """Locate the parser: explicit root, then $THERMORAWFILEPARSER, then bundled, then PATH.

        ⚠ The bundled copy used to be found ONLY when a caller passed `root` explicitly, so a bare
        `Converter.find()` raised even with the parser sitting unpacked in the checkout. Nothing in
        the pipeline passes a root, which meant Thermo `.raw` could not be converted unattended —
        and the way that surfaced was two datasets silently never staged, because the stager reads
        polarity from mzML and there were no mzML to read.
        """
        candidates = []
        if root is not None:
            candidates.append(Path(root) / "ThermoRawFileParser")
        env = os.environ.get("THERMORAWFILEPARSER")
        if env:
            env_path = Path(env)
            # Accept either the executable itself or the directory holding it.
            candidates.append(env_path if env_path.is_file() else env_path / "ThermoRawFileParser")
        candidates.append(cls.BUNDLED / "ThermoRawFileParser")

        for candidate in candidates:
            if candidate.is_file():
                if not os.access(candidate, os.X_OK):
                    # Unpacking loses the executable bit often enough that a "not found" here
                    # would send the reader looking for a missing file that is right there.
                    raise ConversionError(
                        f"ThermoRawFileParser found at {candidate} but is not executable. "
                        f"Run: chmod +x {candidate}")
                return cls(candidate)

        found = shutil.which("ThermoRawFileParser")
        if found:
            return cls(Path(found))
        raise ConversionError(
            f"ThermoRawFileParser not found. Download the self-contained LINUX build from "
            f"https://github.com/compomics/ThermoRawFileParser/releases and unpack it into "
            f"{cls.BUNDLED}/ , or set $THERMORAWFILEPARSER. (The release page also offers osx and "
            f"win archives whose names sort ahead of the linux one — take the linux zip.)")

    def convert(self, raw: str | Path, out_dir: str | Path, force: bool = False) -> Path:
        raw = Path(raw)
        out = Path(out_dir) / f"{raw.stem}.mzML"
        out.parent.mkdir(parents=True, exist_ok=True)

        if out.exists() and not force and out.stat().st_mtime >= raw.stat().st_mtime:
            return out

        # Convert to a temporary name and rename only on success. A conversion that is interrupted
        # — the process killed, the machine rebooted — otherwise leaves a truncated .mzML that is
        # newer than its .raw, so the check above trusts it and every later run reuses it. That
        # failure is silent until something tries to parse it, which may be a different stage on a
        # different day. A rename within one directory is atomic, so the cache only ever holds
        # files that finished.
        # The suffix must stay `.mzML`: ThermoRawFileParser appends `.mzML` to `-b` unless the
        # path already ends in it, so a name like `x.mzML.partial` comes back as
        # `x.mzML.partial.mzML` and the check below then reports a successful conversion as failed.
        partial = out.with_name(f"{out.stem}.partial.mzML")
        partial.unlink(missing_ok=True)
        result = subprocess.run(
            [str(self.executable), "-i", str(raw), "-b", str(partial), "-f", INDEXED_MZML],
            capture_output=True, text=True)
        if result.returncode != 0 or not partial.exists():
            partial.unlink(missing_ok=True)
            raise ConversionError(f"conversion failed for {raw.name}:\n{result.stdout}\n{result.stderr}")
        partial.replace(out)
        return out

    def convert_all(self, raws, out_dir, force: bool = False, log=None) -> list[Path]:
        """Convert what needs it, pass through what does not. Order is preserved."""
        out = []
        for raw in raws:
            if not needs_conversion(raw):
                # Used where it lies. Copying it into `mzml_dir` would double the storage for a
                # large study and give two paths for one file, which is how a stale copy ends up
                # being analysed.
                if log:
                    log(f"already mzML, using as is: {Path(raw).name}")
                out.append(Path(raw))
                continue
            if log:
                log(f"converting {Path(raw).name}")
            out.append(self.convert(raw, out_dir, force=force))
        return out
