"""The facility entry point's discovery rules.

Both tests are about picking the *right* file when several match. Getting these wrong is silent:
the run completes, the report is produced, and the injection order or the grouping is quietly from
the wrong source.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

import run_study      # noqa: E402


def test_randomised_sequence_wins_over_the_draft_it_was_made_from(tmp_path):
    """`SeqRand.csv` sits beside `Seq.csv` and only the randomised one describes what the
    instrument actually ran. Alphabetical order picks `Seq.csv`, which is the wrong answer and
    would silently give every downstream drift and carryover figure the wrong injection order."""
    (tmp_path / "Seq.csv").write_text("x")
    (tmp_path / "SeqRand.csv").write_text("x")
    assert Path(run_study.find_sequence(tmp_path)).name == "SeqRand.csv"


def test_plain_sequence_file_is_found(tmp_path):
    (tmp_path / "sequence.csv").write_text("x")
    assert Path(run_study.find_sequence(tmp_path)).name == "sequence.csv"


def test_the_current_convention_wins_over_the_legacy_names(tmp_path):
    """`sequence.csv` is what the batch generator writes now; SeqRand/Seq are older studies."""
    for name in ("sequence.csv", "SeqRand.csv", "Seq.csv"):
        (tmp_path / name).write_text("x")
    assert Path(run_study.find_sequence(tmp_path)).name == "sequence.csv"


def test_no_sequence_is_not_an_error(tmp_path):
    assert run_study.find_sequence(tmp_path) == ""


def test_only_polarity_folders_holding_raw_files_count(tmp_path):
    (tmp_path / "Pos").mkdir()
    (tmp_path / "Pos" / "a.raw").write_text("x")
    (tmp_path / "Neg").mkdir()          # exists but empty — not a polarity to run
    (tmp_path / "Other").mkdir()
    (tmp_path / "Other" / "b.raw").write_text("x")
    found = run_study.polarity_folders(tmp_path)
    assert set(found) == {"Pos"}


def test_an_interrupted_conversion_leaves_no_file_to_trust(tmp_path, monkeypatch):
    """The failure this prevents is silent and delayed. A killed conversion used to leave a
    truncated .mzML newer than its .raw, so the freshness check trusted it and every later run
    reused it; the error surfaced later, in a different stage, as an XML parse failure."""
    import subprocess
    from lipidloop.convert import ConversionError, Converter

    exe = tmp_path / "parser"
    exe.write_text("#!/bin/sh\nexit 1\n")
    exe.chmod(0o755)
    raw = tmp_path / "sample.raw"
    raw.write_text("x")
    out_dir = tmp_path / "mzml"

    def fake(cmd, **kwargs):
        # write a partial file, then report failure, exactly as an interrupted run would.
        # ThermoRawFileParser appends `.mzML` unless the path already ends in it — reproduced
        # here, because getting the temporary name wrong made a successful conversion look failed.
        target = Path(cmd[cmd.index("-b") + 1])
        assert target.suffix == ".mzML", f"parser would rename {target.name}"
        target.write_text("<mzML>truncated")
        return subprocess.CompletedProcess(cmd, 1, "", "killed")

    monkeypatch.setattr(subprocess, "run", fake)
    with pytest.raises(ConversionError):
        Converter(exe).convert(raw, out_dir)
    assert list(out_dir.glob("*.mzML")) == [], "a partial conversion must not be left behind"
    assert list(out_dir.glob("*.partial")) == []


def test_set_with_a_dotted_key_edits_inside_the_block(tmp_path):
    """`--set presence_filter.min_fraction=0.6` must change that one setting and leave metadata and
    group_by intact. A flat update would write a top-level key of that literal name, which nothing
    reads — the flag would look like it worked and do nothing."""
    import json
    path = run_study.build_config(
        tmp_path, tmp_path, "Pos", tmp_path, sequence="", metadata="", group_by="group",
        extra={"modifier": "formate", "presence_filter.min_fraction": 0.6})
    block = json.loads(Path(path).read_text())["presence_filter"]
    assert block["min_fraction"] == 0.6
    assert block["enabled"] is True and "group_by" in block


def test_extra_library_and_decoy_flags_reach_the_config(tmp_path, monkeypatch):
    """★ REGRESSION. `--extra-library` was accepted by argparse and then went nowhere.

    `args.extra_libraries` was never copied into the `extra` dict `main()` builds, so
    `build_config`'s `extra.pop("extra_libraries", ...)` always saw `None` and the flag was a
    silent no-op — accepted on the command line, recorded in no config, changing no search.

    It was not caught by any existing test because every one of them calls `build_config`
    directly with a hand-built `extra` dict, which is exactly the dict that was never assembled
    correctly by `main()`. This one drives `main()` itself through `sys.argv`, the only path that
    exercises the broken wiring.

    Caught by three ULCFA/FAHFA library-displacement runs (2026-08-28) that searched only the
    DECOY library they were paired with — the real target never loaded — and reported "0 real
    hits, 0 decoy hits, 0 displaced", identical to baseline in every number. That reads as a clean
    null result and was actually three runs that tested nothing.
    """
    study = tmp_path / "study"
    (study / "Pos").mkdir(parents=True)
    (study / "Pos" / "sample.mzML").touch()

    captured: dict = {}

    def fake_build_config(study, analysis, polarity, folder, sequence, metadata, group_by, extra):
        captured.update(extra)
        raise SystemExit(0)   # abort before anything tries to actually run

    monkeypatch.setattr(run_study, "build_config", fake_build_config)
    monkeypatch.setattr(
        sys, "argv",
        ["run_study.py", str(study), "--analysis-dir", str(tmp_path / "out"),
         "--extra-library", "LipiDex_HCD_ULCFA",
         "--extra-library", "LipiDex2_FAHFA",
         "--decoy", "DECOY1_LipiDex_HCD_ULCFA"])

    with pytest.raises(SystemExit):
        run_study.main()

    assert captured.get("extra_libraries") == ["LipiDex_HCD_ULCFA", "LipiDex2_FAHFA"], \
        "--extra-library must reach build_config, repeatable, in the order given"
    assert "DECOY1_LipiDex_HCD_ULCFA" in captured.get("decoy_libraries", []), \
        "--decoy must still work — this regression must not silently break its neighbour too"
