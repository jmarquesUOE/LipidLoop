"""Run the whole pipeline: .raw in, Final_Results.csv out.

    python scripts/run_pipeline.py --config run_config.json
    python scripts/run_pipeline.py --raw 'data/*.raw' --libraries data/libraries/*.msp
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.calibrate import apply, derive, evidence_from
from lipidloop.pipeline import RunConfig, run   # noqa: E402
from lipidloop.search import DEFAULT_LIBRARIES  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", help="JSON run configuration")
    ap.add_argument("--raw", nargs="*", default=[], help=".raw files")
    ap.add_argument("--libraries", nargs="*", default=[],
                    help=".msp libraries IN LOAD ORDER — it breaks score ties. "
                         f"Defaults to {', '.join(DEFAULT_LIBRARIES)} under data/libraries/")
    ap.add_argument("--fatty-acids", default=str(ROOT / "data/lipidex_src/FattyAcids.csv"))
    ap.add_argument("--out", default="outputs/pipeline")
    args = ap.parse_args()

    if args.config:
        config = RunConfig.load(args.config)
    else:
        if not args.raw:
            ap.error("--raw is required without --config")
        libraries = args.libraries or [str(ROOT / f"data/libraries/{stem}.msp")
                                       for stem in DEFAULT_LIBRARIES]
        missing = [p for p in libraries if not Path(p).exists()]
        if missing:
            ap.error("libraries not found: " + ", ".join(missing))
        config = RunConfig(raw_files=[str(p) for p in args.raw],
                           libraries=[str(p) for p in libraries],
                           fatty_acids=args.fatty_acids,
                           output_dir=args.out)

    t0 = time.time()
    say = lambda m: print(f"[{time.time()-t0:7.1f}s] {m}", flush=True)   # noqa: E731
    output = run(config, log=say)

    # A calibrated run is two passes: measure the descriptive parameters from the first, then run
    # again with them. Only parameters that DESCRIBE the measurement are derived — never one that
    # decides what counts as a hit; `calibrate.FORBIDDEN` enforces that rather than trusting it.
    if config.calibration.enabled:
        say("pass 1 complete — measuring the run")
        report = derive(config, evidence_from(output, config), config.calibration, log=say)
        if report.values:
            calibrated = apply(config, report)
            calibrated.calibration.enabled = False      # the second pass must not recurse
            calibrated.save(Path(config.output_dir) / "run_config.calibrated.json")
            say("pass 2 — re-running with the measured parameters")
            output = run(calibrated, log=say)
            output.config.calibration_report = calibrated.calibration_report
            output.config.save(Path(config.output_dir) / "run_config.json")
        else:
            say("nothing could be measured — the single-pass result stands")

    print(f"done in {time.time()-t0:.0f}s -> {output.final_results}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
