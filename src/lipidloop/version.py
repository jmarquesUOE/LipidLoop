"""Which code produced a run, read from the checkout that is running — never from a label passed in.

⚠ Runs used to be named by hand after the commit someone believed was checked out, and a run
labelled with the wrong commit could not be told apart from the output. The commit is now written into every config and the QC report by
the code itself, so the record cannot drift from what ran.
"""
from __future__ import annotations

import subprocess
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


@lru_cache(maxsize=1)
def code_version(root: str | Path | None = None) -> str:
    """Short commit of the checkout, `+dirty` when tracked files differ from it.

    Untracked files are ignored: a results folder or a log inside the checkout does not change
    what the code does, and counting them would mark every working copy dirty. `unknown` when the
    package is not running from a git checkout (an installed copy, a zip of the source).
    """
    where = str(root or ROOT)
    try:
        commit = subprocess.run(["git", "-C", where, "rev-parse", "--short", "HEAD"],
                                capture_output=True, text=True, timeout=10)
        if commit.returncode != 0 or not commit.stdout.strip():
            return "unknown"
        status = subprocess.run(["git", "-C", where, "status", "--porcelain", "--untracked-files=no"],
                                capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return "unknown"
    dirty = "+dirty" if status.returncode == 0 and status.stdout.strip() else ""
    return commit.stdout.strip() + dirty
