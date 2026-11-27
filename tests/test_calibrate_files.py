"""Per-file mass calibration must never run on its own output.

The first pass exists to measure the mass axis as the instrument left it. Run against corrected
files it finds offsets near zero, concludes nothing needs doing, and the real drift becomes
invisible — or, if a correction is applied anyway, it lands on top of the previous one. Neither
failure raises anything: the numbers come out wrong and confident.

This stopped being hypothetical when the conversion cache became shared across analyses of a
study, because calibrated output written into that cache would be read by the next first pass as
if it were raw.
"""
from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
SCRIPT = REPO / "scripts/calibrate_files.py"


def test_refuses_to_calibrate_files_it_has_already_calibrated(tmp_path):
    study = tmp_path / "study"
    (study / "Pos").mkdir(parents=True)
    (study / "Pos" / "a.mzML").write_text("<mzML></mzML>")
    (study / "Pos" / "calibration.json").write_text(
        json.dumps({"source": "earlier", "applied_ppm": {"a.mzML": 17.4}}))

    results = tmp_path / "results.csv"
    results.write_text("Identification,Dot Product,Features Found,Retention Time (min),"
                       "Quant Ion,Adduct\n")
    out = tmp_path / "out"

    r = subprocess.run([sys.executable, str(SCRIPT), str(study), "--results", str(results),
                        "--out", str(out)], capture_output=True, text=True, timeout=300)
    assert r.returncode != 0, "a second calibration must fail, not proceed quietly"
    assert not out.exists() or not list(out.rglob("*.mzML")), "nothing may be written"


def test_the_manifest_records_what_was_applied(tmp_path):
    """Without a record, corrected files are indistinguishable from raw ones except by folder name
    — and folder names are exactly what gets lost when data is moved between disks."""
    import importlib.util
    spec = importlib.util.spec_from_file_location("calibrate_files", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    d = tmp_path / "Pos"
    d.mkdir()
    assert module.already_calibrated(d) == {}, "an unmarked directory must read as uncalibrated"
    (d / "calibration.json").write_text(json.dumps({"applied_ppm": {"x.mzML": -3.2}}))
    assert module.already_calibrated(d) == {"x.mzML": -3.2}


def _load_module():
    import importlib.util
    spec = importlib.util.spec_from_file_location("calibrate_files", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


class _FakeSpectrum:
    """`get_peaks()` returns numpy arrays here, matching real pyOpenMS -- a bare numpy float
    leaking out of `measure()` is exactly what made `Calibration.json` unwritable (`TypeError:
    Object of type bool is not JSON serializable`, from a numpy bool comparison), caught live on
    MTBKS222_Thermo. A fake that hands back plain Python floats would never have caught it."""

    def __init__(self, level, rt_sec, peaks, intensities):
        import numpy as np
        self._level, self._rt = level, rt_sec
        self._peaks, self._int = np.array(peaks, dtype="float64"), np.array(intensities, dtype="float64")

    def getMSLevel(self): return self._level
    def getRT(self): return self._rt
    def get_peaks(self): return (self._peaks, self._int)


class _FakeExperiment:
    def __init__(self): self._specs = []
    def __iter__(self): return iter(self._specs)


def _fake_pyopenms(spectra_by_path):
    """A `pyopenms` stand-in whose `MzMLFile().load(path, experiment)` looks up canned spectra by
    path, so `measure()` can be driven without a real mzML file or a real MS engine."""
    import types

    class FakeMzMLFile:
        def load(self, path, experiment):
            experiment._specs = spectra_by_path.get(str(path), [])

    return types.SimpleNamespace(MzMLFile=FakeMzMLFile, MSExperiment=_FakeExperiment)


def test_calibrate_study_isolates_the_polarity_it_is_given(tmp_path, monkeypatch):
    """The bug this exists to prevent: one invocation walking BOTH Pos and Neg unconditionally
    matches Neg files against Pos-derived reference masses (or vice versa) and reports a wrong
    correction, silently. Confirmed real on MTBKS222_Waters -- combined-invocation Neg spread was
    -12.3 to +5.4 ppm, isolated-Neg-only spread was -0.4 to +9.0 ppm, a different answer, not just
    a noisier one. `polarities=("Neg",)` must touch only `study/Neg`, never `study/Pos`."""
    module = _load_module()

    study = tmp_path / "study"
    (study / "Pos").mkdir(parents=True)
    (study / "Neg").mkdir(parents=True)
    (study / "Pos" / "a.mzML").write_text("<mzML></mzML>")
    (study / "Neg" / "b.mzML").write_text("<mzML></mzML>")

    ref_mz = 760.5851
    # A Pos file whose peak sits far off (would report a large, wrong ppm if ever read) and a Neg
    # file whose peak sits close to the reference -- if isolation works, only the Neg file's
    # (small, correct) offset is ever measured, and Pos is never touched at all.
    pos_path = str((study / "Pos" / "a.mzML").resolve())
    neg_path = str((study / "Neg" / "b.mzML").resolve())
    spectra_by_path = {
        pos_path: [_FakeSpectrum(1, 300.0, [ref_mz * 1.001], [1e6])] * 25,
        neg_path: [_FakeSpectrum(1, 300.0, [ref_mz * (1 + 3e-6)], [1e6])] * 25,
    }
    monkeypatch.setitem(sys.modules, "pyopenms", _fake_pyopenms(spectra_by_path))

    results = tmp_path / "results.csv"
    results.write_text(
        "Identification,Dot Product,Features Found,Retention Time (min),Quant Ion,Adduct\n"
        + "".join(f"Lipid {i},950,20,5.0,{ref_mz},[M-H]-\n" for i in range(25)))

    out = tmp_path / "out"
    result = module.calibrate_study(study, results, out, min_files=1,
                                    polarities=("Neg",), measure_only=True, say=lambda *a: None)

    assert set(result["per_file"]) == {"b.mzML"}, "Pos must never be read when isolated to Neg"
    # ~3 ppm, not the ~1000 ppm a Pos file's peak would give if isolation had failed and Neg's
    # measurement leaked in reference masses matched against the wrong file.
    assert abs(result["per_file"]["b.mzML"]["ppm"] - 3.0) < 1.0
    json.dumps(result)   # Calibration.json's own write path -- must never raise on a numpy type


def test_calibrate_study_refuses_with_too_few_calibrants(tmp_path):
    module = _load_module()
    study = tmp_path / "study"
    (study / "Neg").mkdir(parents=True)
    results = tmp_path / "results.csv"
    results.write_text("Identification,Dot Product,Features Found,Retention Time (min),"
                       "Quant Ion,Adduct\n")
    out = tmp_path / "out"
    result = module.calibrate_study(study, results, out, polarities=("Neg",), say=lambda *a: None)
    assert result["refused"] == "too few calibrants"
    assert not out.exists() or not list(out.rglob("*.mzML"))
