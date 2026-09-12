

def test_find_locates_the_bundled_parser_without_being_told_where(tmp_path, monkeypatch):
    """A bare `Converter.find()` must work when the parser is unpacked in the checkout.

    It used to find the bundled copy ONLY when a caller passed `root`, and nothing in the pipeline
    passes one — so Thermo `.raw` could not be converted unattended even with the parser sitting
    right there. The way that surfaced was two datasets silently absent from every batch for
    months: no conversion meant no mzML, the stager reads polarity from mzML, so it filed every
    file as "Unknown" and staged nothing. No error was raised at any point.
    """
    import os
    from lipidloop.convert import Converter, ConversionError

    monkeypatch.delenv("THERMORAWFILEPARSER", raising=False)

    bundled = tmp_path / "bundled" / "ThermoRawFileParser"
    bundled.mkdir(parents=True)
    exe = bundled / "ThermoRawFileParser"
    exe.write_text("#!/bin/sh\n"); exe.chmod(0o755)
    monkeypatch.setattr(Converter, "BUNDLED", bundled)

    assert Converter.find().executable == exe, "the bundled copy must be found with no arguments"

    # An explicit root still wins over the bundle.
    other = tmp_path / "other"
    other.mkdir()
    other_exe = other / "ThermoRawFileParser"
    other_exe.write_text("#!/bin/sh\n"); other_exe.chmod(0o755)
    assert Converter.find(root=other).executable == other_exe

    # $THERMORAWFILEPARSER accepts the executable itself or the directory holding it.
    monkeypatch.setenv("THERMORAWFILEPARSER", str(other_exe))
    assert Converter.find().executable == other_exe
    monkeypatch.setenv("THERMORAWFILEPARSER", str(other))
    assert Converter.find().executable == other_exe
    monkeypatch.delenv("THERMORAWFILEPARSER")

    # Unpacking loses the executable bit often enough that "not found" would send the reader
    # hunting for a file that is plainly there. Say what is actually wrong.
    exe.chmod(0o644)
    try:
        Converter.find()
    except ConversionError as e:
        assert "not executable" in str(e) and "chmod" in str(e)
    else:
        raise AssertionError("a non-executable parser must not be returned")

    # Nothing anywhere: the message must name the LINUX build, because the release page offers
    # osx-arm64 and osx archives whose names sort ahead of it and installing one fails obscurely.
    exe.unlink()
    monkeypatch.setattr("shutil.which", lambda name: None)
    try:
        Converter.find()
    except ConversionError as e:
        assert "LINUX" in str(e) and str(bundled) in str(e)
    else:
        raise AssertionError("expected ConversionError when no parser exists")


def test_conversion_cache_is_shared_across_analyses(tmp_path):
    """The cache must sit beside the STUDY, or every re-run reconverts from scratch.

    The comment in run_study always said the point was for "the cache to be found again when the
    study is reprocessed", and the code put it under the analysis directory instead — so a new
    --analysis-dir started empty every time. CKD and Skin are 178 Thermo .raw between them; a v3/v4
    pair plus a calibration second pass paid for the same conversion three times over.
    """
    import json
    import sys
    from pathlib import Path as _Path
    sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "scripts"))
    from run_study import build_config

    study = tmp_path / "study"
    (study / "Pos").mkdir(parents=True)
    dirs = []
    for name in ("Analysis_a", "Analysis_b"):
        analysis = study / name
        analysis.mkdir()
        path = build_config(study, analysis, "Pos", study / "Pos", "", "", "group", {})
        dirs.append(json.loads(_Path(path).read_text())["mzml_dir"])

    assert dirs[0] == dirs[1], "two analyses of one study must share one conversion cache"
    assert "Analysis_" not in dirs[0], "the cache must not live inside an analysis directory"


def test_byte_identical_files_are_staged_once(tmp_path):
    """De-duplicating by stem is not enough when a deposit mirrors its own files.

    MSV000094718 stages every vendor injection twice — once under a short name, once under its full
    MassIVE path — different stems, identical bytes. Its negative arm is 9 real injections presented
    as 18. A presence filter requiring "detected in at least 2 samples" is then satisfied by one
    injection counted twice, and every per-file CV treats a copy as an independent replicate.
    """
    import sys
    from pathlib import Path as _Path
    sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "scripts"))
    from run_study import acquisitions

    folder = tmp_path / "Pos"
    folder.mkdir()
    content = "<mzML>" + "x" * 5000 + "</mzML>"
    (folder / "run_1.mzML").write_text(content)
    (folder / "mirrored_path_prefix_run_1.mzML").write_text(content)     # same bytes, other name
    (folder / "run_2.mzML").write_text(content.replace("x", "y"))        # genuinely different

    kept = acquisitions(folder)
    assert len(kept) == 2, f"expected 2 unique acquisitions, got {[p.name for p in kept]}"
    # The deposit's own short name survives, not the mirrored path.
    assert any(p.name == "run_1.mzML" for p in kept), "the shorter name should be the one kept"


def test_same_size_different_content_is_not_deduplicated(tmp_path):
    """Size is only the cheap pre-filter — two files of equal length must both survive."""
    import sys
    from pathlib import Path as _Path
    sys.path.insert(0, str(_Path(__file__).resolve().parents[1] / "scripts"))
    from run_study import acquisitions

    folder = tmp_path / "Neg"
    folder.mkdir()
    (folder / "a.mzML").write_text("<mzML>" + "a" * 4000 + "</mzML>")
    (folder / "b.mzML").write_text("<mzML>" + "b" * 4000 + "</mzML>")
    assert len(acquisitions(folder)) == 2
