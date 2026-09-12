"""Feature detection's tolerance for damaged input.

A study is a folder of files someone downloaded, and downloads are interrupted. What the pipeline
does with a file it cannot read decides whether one bad download costs one injection or the whole
run.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))


def test_a_truncated_mzml_is_skipped_not_fatal(tmp_path):
    """⚠ One bad download used to abort an entire study.

    MTBLS2016 lost a 147-file run to a single mzML that ended mid-base64 inside a `<binary>`
    element, discarding ~15 minutes of completed feature detection — and it had TWO such files, so
    fixing the one named in the traceback would have failed again later in the same run.

    The check is a 200-byte tail read rather than a try/except, because pyOpenMS is not reliably
    survivable on a damaged file: one truncation raises `Parse Error`, which can be caught, and
    another **segfaults the interpreter**, which cannot.
    """
    from lipidloop.features import file_polarity, mzml_is_complete
    bad = tmp_path / "truncated.mzML"
    bad.write_bytes(b'<?xml version="1.0"?><indexedmzML><mzML><binary>AAAA')
    assert mzml_is_complete(bad) is False
    assert file_polarity(bad) == ""          # skipped, not raised, not segfaulted

    good = tmp_path / "whole.mzML"
    good.write_bytes(b'<?xml version="1.0"?><indexedmzML><mzML></mzML></indexedmzML>')
    assert mzml_is_complete(good) is True
