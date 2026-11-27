"""Every run records the commit of the code that ran (v1.1 item 21).

On 2026-11-12 a set of runs was labelled with one commit while the checkout was at the next, and
nothing in the output could tell. The record now comes from the checkout itself.
"""
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.pipeline import RunConfig           # noqa: E402
from lipidloop.version import code_version         # noqa: E402
import qc_report                                    # noqa: E402

version = code_version.__wrapped__      # uncached: each test changes the repository under it


def _git(repo, *args):
    # Identity passed per command: the test must not depend on, or write, any git configuration.
    return subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
                           *args], capture_output=True, text=True, check=True).stdout.strip()


def _repo(tmp_path):
    _git(tmp_path, "init", "-q")
    (tmp_path / "code.py").write_text("x = 1\n")
    _git(tmp_path, "add", "code.py")
    _git(tmp_path, "commit", "-q", "-m", "one")
    return tmp_path


def test_clean_checkout_gives_its_short_commit(tmp_path):
    repo = _repo(tmp_path)
    assert version(repo) == _git(repo, "rev-parse", "--short", "HEAD")


def test_edited_tracked_file_is_marked_dirty(tmp_path):
    repo = _repo(tmp_path)
    (repo / "code.py").write_text("x = 2\n")
    assert version(repo).endswith("+dirty")


def test_untracked_files_do_not_mark_it_dirty(tmp_path):
    """Results and logs written inside a checkout do not change what the code does."""
    repo = _repo(tmp_path)
    (repo / "run.log").write_text("output")
    assert "+dirty" not in version(repo)


def test_outside_a_checkout_says_unknown(tmp_path):
    assert version(tmp_path) == "unknown"


def test_run_config_keeps_the_commit_through_save_and_load(tmp_path):
    config = RunConfig(raw_files=[], libraries=[], output_dir=str(tmp_path))
    config.code_version = "abc1234"
    config.save(tmp_path / "run_config.json")
    assert json.loads((tmp_path / "run_config.json").read_text())["code_version"] == "abc1234"
    assert RunConfig.load(tmp_path / "run_config.json").code_version == "abc1234"


class _Arm:
    def __init__(self, label, results):
        self.label, self.results = label, str(results)


def test_report_states_one_commit_or_each_when_they_differ(tmp_path):
    for pol, v in (("Pos", "abc1234"), ("Neg", "abc1234")):
        (tmp_path / pol).mkdir()
        (tmp_path / pol / "run_config.json").write_text(json.dumps({"code_version": v}))
    arms = [_Arm(p, tmp_path / p / "Final_Results.csv") for p in ("Pos", "Neg")]
    assert qc_report.code_versions_text(arms) == "<code>abc1234</code>"
    (tmp_path / "Neg" / "run_config.json").write_text(json.dumps({"code_version": "def5678"}))
    text = qc_report.code_versions_text(arms)
    assert "Pos <code>abc1234</code>" in text and "Neg <code>def5678</code>" in text


def test_report_says_not_recorded_for_an_old_run(tmp_path):
    (tmp_path / "Pos").mkdir()
    (tmp_path / "Pos" / "run_config.json").write_text(json.dumps({"ms2_tol": 0.3}))
    assert "not recorded" in qc_report.code_versions_text([_Arm("Pos", tmp_path / "Pos" / "x.csv")])
