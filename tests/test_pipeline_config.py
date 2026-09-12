"""RunConfig.load() must tolerate config fields it does not itself declare.

`run_study.py` writes audit-only fields onto a saved config (e.g. `modifier_detected`,
`modifier_evidence` from `--detect-modifier`) so the file explains itself without a second
document. Passing every JSON key straight into the constructor made that a breaking change: a
real `--detect-modifier` run crashed on reload with `TypeError: RunConfig.__init__() got an
unexpected keyword argument 'modifier_detected'` because `load()` blindly spread the whole file.
"""
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.pipeline import RunConfig  # noqa: E402


def test_load_ignores_unknown_top_level_keys(tmp_path):
    path = tmp_path / "config.json"
    path.write_text(json.dumps({
        "mzml_dir": "mzml",
        "output_dir": "out",
        "modifier_detected": "acetate",
        "modifier_evidence": {"sodium formate": 6181, "sodium acetate": 4442},
    }))
    config = RunConfig.load(path)
    assert config.mzml_dir == "mzml"
    assert config.output_dir == "out"


def test_load_still_honours_fields_it_declares(tmp_path):
    """The fix must not turn into "accept anything" -- declared fields still take their value."""
    path = tmp_path / "config.json"
    path.write_text(json.dumps({"output_dir": "custom_out", "ms2_tol": 0.3}))
    config = RunConfig.load(path)
    assert config.output_dir == "custom_out"
    assert config.ms2_tol == 0.3
