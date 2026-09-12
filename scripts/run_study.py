"""Run a whole study: point it at a submission folder, get results and a quality report.

    python scripts/run_study.py /path/to/20260813_Client_Lipidomics \
        --study "Client lipidomics" --owner "Dr Someone"

This is the facility entry point. `run_pipeline.py` processes one polarity from one config;
everything above that — finding the raw files, writing a config per polarity, keeping the settings
that were validated, running both, building the report, laying the folders out the same way every
time — was being done by hand for each study, which is how two studies end up processed
differently without anyone deciding to.

What it does:

    <study>/Pos/*.raw, <study>/Neg/*.raw    discovered
    <study>/sequence.csv                    found if present: injection order AND Sample Type
    <study>/Analysis_<date>/
        config_Pos.json, config_Neg.json    written, and kept, so the run is reproducible
        run_Pos.log, run_Neg.log
        Results/Pos/, Results/Neg/          the tables, plus search/
        Results/QC_report.pdf               both polarities in one document

**Settings are the validated ones and are not invented per study.** They come from `DEFAULTS`
below, which is the configuration the reference datasets were processed with. A study that needs
something different should say so on the command line, so the difference is recorded rather than
discovered later.

**Grouping.** The presence filter keeps a lipid detected in most of *one group*, so that a species
real in one arm and absent in the other survives. That needs to know the groups. Pass
`--metadata FILE --group-by COLUMN`; without it the samples form a single group, which is stated
in the log and in the report rather than done silently. If no metadata is given, a template listing
the sample files is written for you to fill in and re-run.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter
import json
import subprocess
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

MIN_CORRECTION = 2.0   # ppm; mirrors calibrate_files.py's own floor, for the log line only
DEFAULT_MODIFIER = "formate"      # the facility method
LIB = ROOT / "data/libraries"

# The configuration the reference datasets were validated with. Anything not named here takes the
# RunConfig default. Kept in one place so every study in the facility is processed the same way
# unless someone deliberately says otherwise.
DEFAULTS = {
    "fatty_acids": str(ROOT / "data/lipidex_src/FattyAcids.csv"),
    # ⚠ NOT a default. A study that did not spike standards must not be searched for them, and
    # must not have its mass calibration cross-checked against compounds that are not there.
    # `--standards FILE` declares them per study; `run_study` leaves this empty otherwise.
    "standards": "",
    # The class retention-window filter removes long-chain species by construction — see
    # docs/FILTERS.md — so it is off, and the per-class retention model does the job instead.
    "rt_filter": False,
    "retention_model_extend": True,
    "retention_model_extend_ppm": 10.0,
    # Compound Discoverer used 3 on this facility's own workflows. 5 was stricter than anything
    # the facility had ever actually run.
    "blank_filter": {"enabled": True, "multiplier": 3.0, "statistic": "mean",
                     "require_detected_in": 2, "drop_blank_columns": False},
}

# Libraries by MOBILE-PHASE MODIFIER, because negative-mode adducts follow the modifier and a
# library built for one does not describe the other.
#
# ⚠ This is not a preference. On a dataset run with ammonium ACETATE, searching the formate
# libraries named 19 acetate adducts of ordinary ceramides as `Cer[AP]` PHYTOCERAMIDES — matching
# the true masses within 0.0-3.9 ppm at retention times agreeing to +/-0.03 min, with dot products
# of 999 and `Purity = 0`. `LipiDex_HCD_Formic` contains **zero** ceramide `[M+Ac-H]-` entries;
# `LipiDex_HCD_Acetate` contains **2,120**. The identification could not have been right.
#
# The facility runs formate, so that stays the default and its results do not move.
LIBRARIES = {
    "formate": {
        "Pos": ["LipidBlast_Formic", "LipiDex_HCD_Formic", "LipiDex_HCD_Hydroxy",
                "LipiDex2_Ganglioside", "Ceramides_Positive",
                # ⚠ WITHDRAWN 2026-08-27 — "UltraLongCer_EOS_Positive", "UltraLongHexCer_EOS_Positive"
                # FAILED the target-decoy test on the one dataset they exist for. On Skin_QEplus:
                # 20 real EOS matches against 197 decoy hits (90.8% decoys), and the DECOYS SCORED
                # BETTER than the targets (median dot 832 vs 700). A library whose decoys outscore
                # its targets is matching on mass alone. They tripled that study's FDR, 1.76% ->
                # 6.87%, and put 63 decoy EOS rows into the delivered table against 12 real ones —
                # every real one a sum composition at Purity 0 with no chain evidence.
                # These were added to close the stratum-corneum barrier-ceramide gap (every other
                # ceramide library we hold stops at 26 C). The gap is real; these libraries do not
                # close it. Left on disk for investigation — they carry no chain annotations, which
                # is why their decoys had to be built by differencing homologues, and that is the
                # first thing to look at. Re-enable only with a decoy measurement that passes.
                "AcylSM_Positive", "AcylHexCer_Positive", "Acylcarnitines_Positive",
                "Diacylglycerols_Positive", "NAcylAmides_Positive"],
        "Neg": ["LipidBlast_Formic", "LipiDex_HCD_Formic", "LipiDex_HCD_Hydroxy",
                "Ganglioside_Negative_SumComposition", "BileAcids_Negative",
                "NAcylAmides_Negative"],
    },
    "acetate": {
        "Pos": ["LipidBlast_Acetate", "LipiDex_HCD_Acetate", "LipiDex_HCD_Hydroxy",
                "LipiDex2_Ganglioside", "Ceramides_Positive",
                # ⚠ WITHDRAWN 2026-08-27 — "UltraLongCer_EOS_Positive", "UltraLongHexCer_EOS_Positive"
                # FAILED the target-decoy test on the one dataset they exist for. On Skin_QEplus:
                # 20 real EOS matches against 197 decoy hits (90.8% decoys), and the DECOYS SCORED
                # BETTER than the targets (median dot 832 vs 700). A library whose decoys outscore
                # its targets is matching on mass alone. They tripled that study's FDR, 1.76% ->
                # 6.87%, and put 63 decoy EOS rows into the delivered table against 12 real ones —
                # every real one a sum composition at Purity 0 with no chain evidence.
                # These were added to close the stratum-corneum barrier-ceramide gap (every other
                # ceramide library we hold stops at 26 C). The gap is real; these libraries do not
                # close it. Left on disk for investigation — they carry no chain annotations, which
                # is why their decoys had to be built by differencing homologues, and that is the
                # first thing to look at. Re-enable only with a decoy measurement that passes.
                "AcylSM_Positive", "AcylHexCer_Positive", "Acylcarnitines_Positive",
                "Diacylglycerols_Positive", "NAcylAmides_Positive"],
        "Neg": ["LipidBlast_Acetate", "LipiDex_HCD_Acetate", "LipiDex_HCD_Hydroxy",
                "Ganglioside_Negative_SumComposition", "BileAcids_Negative",
                "BileAcids_Negative_Acetate", "NAcylAmides_Negative"],
    },
    # A validation-only set: the acetate libraries with a sample of confidently identified lipids
    # removed. Asks what the search does when the right answer is ABSENT — the failure mode a
    # target-decoy is structurally blind to.
    "acetate-holdout": {
        "Pos": ["HOLDOUT_LipidBlast_Acetate", "HOLDOUT_LipiDex_HCD_Acetate",
                "HOLDOUT_LipiDex_HCD_Hydroxy", "LipiDex2_Ganglioside", "Ceramides_Positive",
                # ⚠ EOS libraries withdrawn here too — see the note in the formate set. This arm
                # must stay comparable with the ones it is a holdout OF, and a library that only
                # this set searches would confound "the right answer was removed" with "a bad
                # library took the match".
                "AcylSM_Positive", "AcylHexCer_Positive", "Acylcarnitines_Positive",
                "Diacylglycerols_Positive", "NAcylAmides_Positive"],
        "Neg": ["HOLDOUT_LipidBlast_Acetate", "HOLDOUT_LipiDex_HCD_Acetate",
                "HOLDOUT_LipiDex_HCD_Hydroxy", "Ganglioside_Negative_SumComposition"],
    },
}

# Spiked internal standards, matched to the modifier for the same reason. Searching them costs
# nothing (21 entries) and gives the run a set of masses known BEFORE it started — the only
# non-circular input a mass calibration can have, and a run-quality check that needs no biology.
ISTD_LIBRARY = {"formate": "LipiDex_Splash_ISTD_Formic",
                "acetate": "LipiDex_Splash_ISTD_Acetate",
                "acetate-holdout": "HOLDOUT_LipiDex_Splash_ISTD_Acetate"}

# Decoys must be drawn from the SAME libraries the targets are, or the comparison is between a
# handicapped target set and a fair decoy one. `--decoy-matched` picks the pair for whichever
# modifier is in force, so a formate study is never scored against acetate decoys.
#
# ⚠ Every target library added needs its decoy added here too. The FDR is targets scored against
# decoys drawn from the same space; growing the target set alone biases the estimate LOW, because
# the new entries can win matches while nothing competes with them. The acylcarnitines added 86
# targets, so their 81 decoys go in beside them.
MATCHED_DECOYS = {
    "formate": ["DECOY1_LipiDex_HCD_Formic", "DECOY1_LipidBlast_Formic",
                "DECOY1_Acylcarnitines_Positive", "DECOY1_Diacylglycerols_Positive",
                "DECOY1_NAcylAmides_Positive", "DECOY1_NAcylAmides_Negative",
                # ⚠ Covers the classes the CH2 construction cannot reach — SM alone contributed 910
                # identifications with no decoy at all, so its error rate was uncountable rather
                # than low. Kept as its own library so its hits stay separable: the construction is
                # new and under test, and the FDR can be recomputed with or without it.
                # Built 2026-08-27. Both carry chain annotations already, so they take the CH2
                # construction — they were simply never built. Together they close ~20,000 of the
                # ~48,000 formate-pool entries that had no decoy of any kind, including the five
                # [OH] classes, which had none in EITHER polarity.
                "DECOY1_LipiDex_HCD_Hydroxy", "DECOY1_LipiDex2_Ganglioside",
                # Built 2026-08-27 by differencing homologues: these four carry NO chain
                # annotations, so the chain peaks were recovered from the library itself — two
                # entries one CH2 apart must have their chain peaks 14.0157 apart and their
                # head-group peaks identical. 25,094 decoys, all verified: precursor preserved,
                # every moved peak shifted by exactly 3xCH2.
                "DECOY1_AcylHexCer_Positive", "DECOY1_AcylSM_Positive",
                # Their targets are withdrawn (see above), so their decoys go too: decoys must be
                # drawn from the same space as the targets or the FDR is measured against a library
                # nothing is being searched from.
                "DECOY2_HeadGroup_LipiDex_HCD_Formic"],
    "acetate": ["DECOY1_LipiDex_HCD_Acetate", "DECOY1_LipidBlast_Acetate",
                "DECOY1_Acylcarnitines_Positive", "DECOY1_Diacylglycerols_Positive",
                "DECOY1_NAcylAmides_Positive", "DECOY1_NAcylAmides_Negative",
                # ⚠ Covers the classes the CH2 construction cannot reach — SM alone contributed 910
                # identifications with no decoy at all, so its error rate was uncountable rather
                # than low. Kept as its own library so its hits stay separable: the construction is
                # new and under test, and the FDR can be recomputed with or without it.
                # Built 2026-08-27. Both carry chain annotations already, so they take the CH2
                # construction — they were simply never built. Together they close ~20,000 of the
                # ~48,000 formate-pool entries that had no decoy of any kind, including the five
                # [OH] classes, which had none in EITHER polarity.
                "DECOY1_LipiDex_HCD_Hydroxy", "DECOY1_LipiDex2_Ganglioside",
                # Built 2026-08-27 by differencing homologues: these four carry NO chain
                # annotations, so the chain peaks were recovered from the library itself — two
                # entries one CH2 apart must have their chain peaks 14.0157 apart and their
                # head-group peaks identical. 25,094 decoys, all verified: precursor preserved,
                # every moved peak shifted by exactly 3xCH2.
                "DECOY1_AcylHexCer_Positive", "DECOY1_AcylSM_Positive",
                # Their targets are withdrawn (see above), so their decoys go too: decoys must be
                # drawn from the same space as the targets or the FDR is measured against a library
                # nothing is being searched from.
                "DECOY2_HeadGroup_LipiDex_HCD_Acetate"],
}
# Free fatty acids are named from mass plus a fitted retention surface, and the library has no
# positive-mode entries for the class, so it is a negative-mode arm only.
FATTY_ACIDS = {"Neg": [str(LIB / "FreeFattyAcids_Negative.msp")], "Pos": []}


# Input the pipeline can start from. Thermo `.raw` is converted; `.mzML` is used as it lies.
# ⚠ Discovery used to look for `*.raw` only, so a folder of perfectly good mzML was invisible and
# the entry point reported "no Pos/ or Neg/ folder with .raw files" — which reads as an empty study
# rather than an unsupported format. Public data, collaborators' data and any instrument exporting
# mzML directly were all unusable through this script.
INPUT_PATTERNS = ("*.raw", "*.mzML", "*.mzml")


def acquisitions(folder: Path) -> list[Path]:
    """Every file in one polarity folder the pipeline can take, de-duplicated and sorted.

    A study holding both `X.raw` and a previously converted `X.mzML` must not process X twice, so
    when a stem appears in more than one format the mzML wins — it is what conversion would have
    produced anyway, and using it skips the work.
    """
    by_stem: dict[str, Path] = {}
    for pattern in INPUT_PATTERNS:
        for path in folder.glob(pattern):
            stem = path.stem
            if stem in by_stem and by_stem[stem].suffix.lower() == ".mzml":
                continue                       # already have the converted form
            by_stem[stem] = path
    return _drop_identical(sorted(by_stem.values(), key=lambda p: p.stem))


def _drop_identical(paths: list[Path]) -> list[Path]:
    """Remove files that are byte-identical to one already kept, whatever they are called.

    ⚠ De-duplicating by STEM is not enough. MSV000094718 stages every vendor injection twice —
    once under a short name and once under its full MassIVE path, `raw_BLE_Vendor1_1_Neg.mzML` and
    `updates_2024-07-17_joshuaroberts_..._BLE_Vendor1_1_Neg.mzML` — different stems, identical
    bytes. Its negative arm is **9 real injections presented as 18**.

    That is not a cosmetic miscount. A presence filter requiring "detected in at least 2 samples"
    is satisfied by one injection counted twice, and every per-file CV treats a copy as an
    independent replicate, understating variance.

    Size first, hash only on collision: hashing every file in a 224-injection study to find
    duplicates that usually do not exist would cost more than it saves.
    """
    import hashlib

    by_size: dict[int, list[Path]] = {}
    for path in paths:
        try:
            by_size.setdefault(path.stat().st_size, []).append(path)
        except OSError:
            by_size.setdefault(-1, []).append(path)

    drop: set[Path] = set()
    for size, group in by_size.items():
        if size < 0 or len(group) < 2:
            continue
        seen: dict[str, Path] = {}
        for path in group:
            try:
                with path.open("rb") as fh:
                    digest = hashlib.md5(fh.read(4_000_000)).hexdigest()
            except OSError:
                continue
            if digest in seen:
                # Keep the shorter name: the deposit's own, rather than the path it was mirrored
                # under. Ties broken alphabetically so the choice is reproducible.
                keep, other = sorted((seen[digest], path), key=lambda p: (len(p.name), p.name))
                seen[digest] = keep
                drop.add(other)
            else:
                seen[digest] = path

    if drop:
        # ⚠ Say what was dropped. A silently smaller file count reads as staging having gone wrong.
        print(f"  ⚠ {len(drop)} duplicate file(s) ignored — byte-identical to another staged file")
        for path in sorted(drop)[:4]:
            print(f"      {path.name}")
    return [p for p in paths if p not in drop]


def polarity_folders(study: Path) -> dict[str, Path]:
    found = {}
    for name in ("Pos", "Neg"):
        folder = study / name
        if folder.is_dir() and acquisitions(folder):
            found[name] = folder
    return found


def find_sequence(study: Path) -> str:
    """The sequence actually run, which is not always the first file that matches.

    `sequence.csv` is the current convention and wins outright. Older studies have a randomised
    sequence sitting beside the un-randomised draft it was made from — `SeqRand.csv` next to
    `Seq.csv` — where only the randomised one describes what the instrument did, and alphabetical
    order picks the wrong one, so the randomised name is tried before the plain one.
    """
    # ⚠ Look inside the polarity folders too. A sequence often sits BESIDE the raw files it
    # describes rather than at the study root — `Pos/Seq_Rand.csv` — and searching the root alone
    # returns nothing, so injection order goes missing and every drift and confounding section of
    # the report says "not testable" when the file was there all along.
    roots = [study] + [study / p for p in ("Pos", "Neg") if (study / p).is_dir()]
    for pattern in ("sequence.csv", "*and*.csv", "Seq*.csv", "*equence*.csv"):
        for root in roots:
            for candidate in sorted(root.glob(pattern)):
                if "seq" in candidate.name.lower():
                    return str(candidate)
    return ""


def find_submission_form(study: Path) -> str:
    """The SSF the submitter emailed in. Named per submission, so matched on shape not name."""
    for candidate in sorted(study.glob("*.xlsx")):
        if candidate.name.startswith("~$"):        # an Excel lock file, not a form
            continue
        return str(candidate)
    return ""


def metadata_from_form(form: str, analysis: Path, folders: dict[str, Path],
                       group_by: str) -> tuple[dict[str, str], object]:
    """Write a metadata CSV per polarity from the submission form.

    The submitter already declared the design on the form. Retyping it into a CSV is how a study
    gets analysed against groups nobody checked.
    """
    # ⚠ The facility layer is OPTIONAL. `facility/` holds everything specific to the Edinburgh IGC
    # — the submission-form reader and the branded report — so the pipeline itself can be run by
    # anyone without it. A missing facility layer is a stated absence, not a crash: the study simply
    # gets a metadata template to fill in, which is the same path a study with no form takes.
    sys.path.insert(0, str(ROOT))
    try:
        from facility import submission as ssf
    except ImportError:
        print("  no facility layer — submission form not read; fill in the metadata template")
        return {}, None

    parsed = ssf.read(form)
    if not parsed or not parsed.samples:
        return {}, None
    written: dict[str, str] = {}
    for polarity, folder in folders.items():
        stems = [p.stem for p in acquisitions(folder)]
        rows = ssf.to_metadata(parsed, stems, group_column=group_by)
        if not rows:
            continue
        path = analysis / f"metadata_{polarity}.csv"
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            writer.writeheader()
            writer.writerows(rows)
        written[polarity] = str(path)
        print(f"  {polarity}: matched {len(rows)} of {len(stems)} files to the form "
              f"(the rest are blanks, pools and standards, identified by role)")
    return written, parsed


def write_metadata_template(study: Path, analysis: Path, folders: dict[str, Path],
                            group_by: str) -> None:
    """A fillable metadata file, so the next run can group properly.

    Sample files only — blanks, pools and standards are identified by role and are not part of a
    biological group.
    """
    from lipidloop.blanks import infer_role, read_sample_types

    sequence = find_sequence(study)
    types = read_sample_types(sequence) if sequence else {}
    for polarity, folder in folders.items():
        path = analysis / f"metadata_{polarity}.template.csv"
        rows = []
        for raw in acquisitions(folder):
            stem = raw.stem
            if infer_role(stem, types.get(stem, "")) != "sample":
                continue
            rows.append({"sample": stem, group_by: ""})
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=["sample", group_by])
            writer.writeheader()
            writer.writerows(rows)
        print(f"  wrote {path.name} — fill in `{group_by}` and re-run with "
              f"--metadata-{polarity.lower()}")


def build_config(study: Path, analysis: Path, polarity: str, folder: Path,
                 sequence: str, metadata: str, group_by: str, extra: dict) -> Path:
    config = dict(DEFAULTS)
    config["raw_files"] = [str(p) for p in acquisitions(folder)]
    modifier = extra.pop("modifier", DEFAULT_MODIFIER)
    stems = list(LIBRARIES[modifier][polarity])
    # ⚠ Extra libraries are OPT-IN per run, never added to the defaults above.
    #
    # A library that is generated rather than measured, or imported and not yet validated, must not
    # reach a production run by being listed once and forgotten — the failure would be silent
    # identifications from an invented intensity model. `--extra-library` makes each such use a
    # deliberate act, recorded in that run's own configuration.
    for stem in extra.pop("extra_libraries", None) or []:
        if stem not in stems:
            stems.append(stem)
    # ⚠ Search the spiked-standard libraries ONLY when standards were actually spiked.
    #
    # Defaulting this on created a false positive on the first unspiked study it met: a peak was
    # named `PC d7-18:1_15:0`, a deuterated SPLASH standard, in a sample containing no SPLASH at
    # all. A labelled standard sits a few Da from its endogenous analogue, so an ISTD library will
    # match something on almost any run — and the result is a spiked compound reported as biology,
    # which would then be quantified and tested like biology.
    #
    # Tied to whether a standards file is declared, because that is the same fact stated once: a
    # study that spiked standards has a standards CSV, and one that did not has none.
    if extra.pop("search_istd", None) or config.get("standards"):
        stems.append(ISTD_LIBRARY[modifier])
    # A library of lipids that CANNOT be in the sample. Every identification from it is a false
    # positive by construction, so it measures the false-positive rate directly — a number
    # lipidomics does not usually report because there is normally no way to know. It changes
    # nothing else: the decoy competes for the same spectra as the real libraries, which is the
    # point. Only ever added when asked for.
    for stem in extra.pop("decoy_libraries", []):
        stems.append(stem)
    config["libraries"] = [str(LIB / f"{stem}.msp") for stem in stems]
    config["fatty_acid_libraries"] = FATTY_ACIDS[polarity]
    config["output_dir"] = str(analysis / "Results" / polarity)
    # ⚠ The conversion cache goes BESIDE THE STUDY, not on the system disk. A 79 GB study of
    # Thermo .raw converts to a similar volume of mzML, and `data/mzml/` had already accumulated
    # 86 GB from previous studies against 77 GB free — the next large run fills the disk and takes
    # the machine with it. Keeping it beside the study also means the cache is found again when
    # the study is reprocessed, rather than being regenerated into a path that depends on the
    # working directory.
    #
    # ⚠ Beside the STUDY, not beside the analysis — which is what the paragraph above always
    # intended and the code did not do. Under `analysis/`, every new `--analysis-dir` starts with
    # an empty cache and reconverts the lot: CKD and Skin are 178 Thermo .raw between them, so a
    # v3/v4 pair, a calibration second pass and any library re-run each paid for the same
    # conversion again. Conversion output is deterministic, so there is nothing to isolate.
    config["mzml_dir"] = str(study / "mzml" / polarity)
    config["sequence_file"] = sequence
    config["presence_filter"] = {"enabled": True, "min_fraction": 0.6,
                                 "metadata": metadata or "",
                                 "group_by": [group_by] if metadata else []}
    # ⚠ Dotted keys address a field INSIDE a nested block — `--set presence_filter.min_fraction=0.6`
    # changes one setting and leaves the rest of the block alone. A flat update would instead add a
    # top-level key literally named "presence_filter.min_fraction", which nothing reads: the flag
    # would appear to work and change nothing. Unknown blocks are still created, so a typo surfaces
    # in the written config rather than vanishing.
    for key, value in extra.items():
        if "." in key:
            head, _, leaf = key.partition(".")
            block = config.get(head)
            if not isinstance(block, dict):
                block = {}
            config[head] = {**block, leaf: value}
        else:
            config[key] = value
    path = analysis / f"config_{polarity}.json"
    path.write_text(json.dumps(config, indent=2))
    return path


def run(config: Path, log: Path) -> int:
    with log.open("w") as fh:
        result = subprocess.run(
            [sys.executable, str(ROOT / "scripts/run_pipeline.py"), "--config", str(config)],
            stdout=fh, stderr=subprocess.STDOUT)
    return result.returncode


def build_report(analysis: Path, folders: dict[str, Path], sequence: str,
                 metadata: dict[str, str], args) -> int:
    command = [sys.executable, str(ROOT / "scripts/qc_report.py"),
               "--out", str(analysis / "Results")]
    for polarity in folders:
        results = analysis / "Results" / polarity
        key = polarity.lower()
        command += [f"--results-{key}", str(results / "Final_Results.csv"),
                    f"--unfiltered-{key}", str(results / "Unfiltered_Results.csv"),
                    f"--spectra-{key}", str(results / "Associated_Spectra.csv")]
        if metadata.get(polarity):
            command += [f"--metadata-{key}", metadata[polarity]]
    if sequence:
        command += ["--sequence", sequence]
    for name in ("study", "owner", "instrument", "method"):
        value = getattr(args, name, "")
        if value:
            command += [f"--{name}", value]
    return subprocess.run(command).returncode


def run_auto_calibration(study: Path, analysis: Path, folders: dict[str, Path], sequence: str,
                         metadata: dict[str, str], group_by: str,
                         extra: dict) -> tuple[dict[str, Path], dict[str, dict]]:
    """Pass one (loose 25ppm linking, measurement only) -> per-polarity drift measurement,
    isolated -> corrected mzML only where warranted. See `docs/MASS_CALIBRATION.md` for the
    process this implements and the thresholds it uses; nothing here is a new number, it is the
    manual procedure documented there made automatic.

    Returns `(folders, calibration)`. `folders` is the polarity -> mzML-directory map the caller's
    MAIN pass should read from: the calibrated directory where correction was applied, the
    study's own original folder otherwise -- so the caller's existing per-polarity loop needs no
    other change to become "pass two". `calibration` is `{polarity: result}` from
    `calibrate_files.calibrate_study`, for the caller to write as `Calibration.json` once the main
    pass's `Results/<polarity>/` exists.

    Doubles the runtime of any study it runs on: pass one is a full pipeline pass, and its
    identifications are never reported, only its `Unfiltered_Results.csv` used to pick calibrants.
    That is why this is opt-in (`--auto-calibrate`) rather than a default like the exclusion and
    inclusion lists, which cost nothing beyond a single ordinary pass.
    """
    from calibrate_files import calibrate_study   # scripts/ is already on sys.path[0]

    pass1 = analysis / "_calib_pass1"
    pass1.mkdir(parents=True, exist_ok=True)
    calibrated_root = analysis / "_calibrated_mzml"
    calibration: dict[str, dict] = {}
    new_folders = dict(folders)

    for polarity, folder in folders.items():
        pass1_extra = {**extra, "features.link_mz_tolerance": 25.0}
        config = build_config(study, pass1, polarity, folder, sequence,
                              metadata.get(polarity, ""), group_by, pass1_extra)
        log = pass1 / f"run_{polarity}.log"
        print(f"\nauto-calibrate: pass one ({polarity}, 25ppm linking, measurement only) -> {log}")
        code = run(config, log)
        print(f"  exit {code}")
        if code != 0:
            print(f"  ⚠ auto-calibrate pass one failed for {polarity}; running {polarity} "
                  f"uncalibrated", file=sys.stderr)
            calibration[polarity] = {"refused": f"pass one failed, exit {code}"}
            continue

        results_csv = pass1 / "Results" / polarity / "Unfiltered_Results.csv"
        n_files = len(acquisitions(folder))
        # `--min-files` defaults to 10 in the manual procedure, which no study under 10 files per
        # polarity can ever satisfy (MTBKS222_Waters: 9 Pos/8 Neg, needed --min-files 4 by hand).
        # Scale it to the study instead of carrying a fixed default that silently refuses small
        # studies.
        min_files = max(3, min(10, n_files // 2))

        # Isolate: `calibrate_study` walks every polarity it is given, so calibrating a
        # two-polarity study in one call cross-contaminates -- Neg files measured against
        # Pos-derived reference masses give a wrong correction, silently (confirmed real on
        # MTBKS222_Waters, see calibrate_files.py's own docstring on this). `isolate` is a
        # throwaway directory holding a symlink to ONLY this polarity's real source folder, under
        # the name `calibrate_study` expects, so a single call here only ever sees one polarity.
        isolate = analysis / "_calib_isolate" / polarity
        isolate.mkdir(parents=True, exist_ok=True)
        link = isolate / polarity
        if not link.exists():
            link.symlink_to(folder.resolve())

        result = calibrate_study(isolate, results_csv, calibrated_root,
                                 min_files=min_files, polarities=(polarity,), say=print)
        calibration[polarity] = result
        if result.get("refused"):
            print(f"  auto-calibrate ({polarity}): {result['refused']} — running uncalibrated")
        elif result.get("corrected_files"):
            print(f"  auto-calibrate ({polarity}): {result['corrected_files']} of "
                  f"{result['total_files']} file(s) corrected, spread {result['spread_ppm']} ppm "
                  f"— pass two will use the calibrated files")
            new_folders[polarity] = calibrated_root / polarity
        else:
            print(f"  auto-calibrate ({polarity}): nothing exceeded the {MIN_CORRECTION} ppm "
                  f"correction floor — running the original files")
    return new_folders, calibration


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("study", help="submission folder holding Pos/ and/or Neg/")
    ap.add_argument("--metadata-pos", default="")
    ap.add_argument("--metadata-neg", default="")
    ap.add_argument("--form", default="",
                    help="submission form (.xlsx); found in the study folder by default")
    ap.add_argument("--group-by", default="group",
                    help="metadata column defining the biological groups (default: group)")
    ap.add_argument("--name", dest="study_name", default="",
                    help="study name for the report (defaults to the folder name)")
    ap.add_argument("--owner", default="")
    ap.add_argument("--instrument", default="Orbitrap Fusion Lumos")
    ap.add_argument("--method", default="Kinetex EVO 150 mm, 25 min")
    ap.add_argument("--modifier", choices=sorted(LIBRARIES), default=DEFAULT_MODIFIER,
                    help="mobile-phase modifier — picks the library set. Negative-mode adducts "
                         "follow the modifier: searching formate libraries on acetate data names "
                         "acetate adducts as different lipids entirely. Default: "
                         f"{DEFAULT_MODIFIER} (the facility method).")
    ap.add_argument("--detect-modifier", action="store_true",
                    help="read the modifier off the background salt clusters and warn if it "
                         "contradicts --modifier. OFF by default: it costs a pass over the MS1 "
                         "scans and most runs already know their modifier from the deposit. Turn "
                         "it on for data that arrives without provenance. It never overrides "
                         "--modifier — it reports, because absence of a series is not evidence "
                         "of the other one.")
    ap.add_argument("--decoy-matched", action="store_true",
                    help="add the decoy libraries built from THIS modifier's targets. Use this "
                         "rather than naming them by hand: a formate study scored against acetate "
                         "decoys compares a handicapped target set with a fair decoy one.")
    ap.add_argument("--extra-library", action="append", default=[], dest="extra_libraries",
                    metavar="STEM",
                    help="add a library by file stem, e.g. OxTG_Positive. Repeatable. For "
                         "generated or unvalidated libraries that must not sit in the defaults.")
    ap.add_argument("--decoy", action="append", default=[], metavar="LIBRARY",
                    help="add a library of lipids that cannot be present, as a negative control. "
                         "Every identification it wins is a false positive by construction. "
                         "`--decoy LipiDex_HCD_Plants` searches galactolipids (MGDG/DGDG/SQDG), "
                         "which occur only in photosynthetic membranes.")
    ap.add_argument("--set", action="append", default=[], metavar="KEY=JSON",
                    help="override any config field in BOTH polarities, e.g. --set rt_filter=true")
    # Per-polarity because this is the setting that genuinely differs between them: on the
    # instrument these were written against, positive sat at -1.0 ppm and negative at -5.8.
    # Measured by the quality report's mass-accuracy section; set here once it is known.
    ap.add_argument("--mass-offset-pos", type=float, default=None,
                    help="systematic mass error to correct, ppm, positive mode")
    ap.add_argument("--mass-offset-neg", type=float, default=None,
                    help="systematic mass error to correct, ppm, negative mode")
    ap.add_argument("--analysis-dir", default="",
                    help="defaults to <study>/Analysis_<today>")
    ap.add_argument("--class-polarity", default="",
                    help="class -> preferred polarity table; defaults to the facility file, with "
                         "<analysis>/class_polarity.csv layered on top when present")
    ap.add_argument("--no-polarity-filter", action="store_true",
                    help="leave both polarities intact; the tables then contain the same class "
                         "measured two ways, on scales that cannot be compared")
    ap.add_argument("--report-only", action="store_true",
                    help="rebuild the report from results already produced")
    ap.add_argument("--auto-calibrate", action="store_true",
                    help="measure per-file mass drift (a throwaway 25ppm-linking pass) and, where "
                         "any file exceeds the 2ppm correction floor, run on corrected mzML "
                         "instead of the originals. See docs/MASS_CALIBRATION.md. Doubles the "
                         "runtime of any study it runs on -- opt in per study, not a default.")
    args = ap.parse_args()

    sys.path.insert(0, str(ROOT / "src"))
    study = Path(args.study).resolve()
    if not study.is_dir():
        ap.error(f"not a directory: {study}")
    folders = polarity_folders(study)
    if not folders:
        ap.error(f"no Pos/ or Neg/ folder with .raw or .mzML files under {study}")

    analysis = Path(args.analysis_dir) if args.analysis_dir else \
        study / f"Analysis_{date.today().isoformat()}"
    analysis.mkdir(parents=True, exist_ok=True)
    sequence = find_sequence(study)

    decoys = list(args.decoy)
    if args.decoy_matched:
        decoys += MATCHED_DECOYS.get(args.modifier, [])
    # ⚠ `--extra-library` was accepted by argparse, populated `args.extra_libraries`, and then
    # went nowhere: nothing ever copied it into `extra`, so `config()`'s
    # `extra.pop("extra_libraries", ...)` always saw None and the flag was a silent no-op. Three
    # ULCFA/FAHFA displacement runs on Skin_QEplus (2026-08-28) searched only the DECOY library
    # they were paired with — the target never loaded — and reported "0 real hits, 0 decoy hits,
    # 0 displaced" identically to baseline, which read as a clean null result and was not one.
    extra: dict = {"modifier": args.modifier, "decoy_libraries": decoys,
                  "extra_libraries": list(args.extra_libraries)}
    for item in args.set:
        key, _, value = item.partition("=")
        try:
            extra[key] = json.loads(value)
        except json.JSONDecodeError:
            extra[key] = value

    metadata = {k: v for k, v in (("Pos", args.metadata_pos), ("Neg", args.metadata_neg)) if v}
    form = args.form or find_submission_form(study)
    print(f"study     {study}")
    print(f"analysis  {analysis}")
    # Say which formats were found, per polarity. A study that is half-converted is normal, and
    # a reader needs to see that rather than a single total that hides it.
    def describe(folder: Path) -> str:
        files = acquisitions(folder)
        kinds = Counter("mzML" if f.suffix.lower() == ".mzml" else f.suffix.lstrip(".").lower()
                        for f in files)
        return " + ".join(f"{n} {k}" for k, n in sorted(kinds.items()))
    counts = ", ".join(f"{name} ({describe(folder)})" for name, folder in folders.items())
    print(f"polarities {counts}")
    print(f"modifier  {args.modifier} — searching the {args.modifier} libraries")
    if args.detect_modifier:
        # ⚠ Reports, never overrides. The formate and acetate adducts differ by exactly CH2, so a
        # wrong modifier renames every negative-mode choline lipid to its odd-chain neighbour with
        # zero mass error — worth a warning at the front, where it is still cheap to act on.
        # Negative mode carries the signal: salt clusters dominate there and are a rounding error
        # in positive, so survey Neg when there is one.
        from lipidloop.modifier import detect_files
        survey_folder = folders.get("Neg") or next(iter(folders.values()), None)
        survey = sorted(Path(survey_folder).glob("*.mzML")) if survey_folder else []
        if not survey:
            print("detect    no mzML to survey — modifier left as declared")
        else:
            found = detect_files(survey)
            print(f"detect    {found}")
            if found.disagrees_with(args.modifier):
                print(f"⚠ WARNING the background says {found.call}, not {args.modifier}. These two "
                      f"adducts differ by exactly CH2, so the wrong one renames negative-mode "
                      f"choline lipids to odd-chain species at zero mass error. Check before "
                      f"trusting this run.")
            extra["modifier_detected"] = found.call or ""
            extra["modifier_evidence"] = found.counts
    if decoys:
        print(f"decoy     {', '.join(decoys)} — every hit is a FALSE POSITIVE by construction")
    print(f"sequence  {sequence or 'none found — injection order and Sample Type unavailable'}")

    parsed = None
    if not metadata and form:
        print(f"form      {Path(form).name}")
        from_form, parsed = metadata_from_form(form, analysis, folders, args.group_by)
        metadata.update(from_form)
        if parsed:
            listed = len(parsed.samples)
            declared = parsed.declared_samples
            if declared and str(declared).strip() != str(listed):
                print(f"  ⚠ the form declares {declared} samples and lists {listed}. "
                      f"dar checks a form against itself before acquisition — this is only a "
                      f"note that the two disagree.")
            print(f"  conditions: {parsed.conditions}")
    if not metadata:
        print("metadata  none — samples form a single group for the presence filter")
        write_metadata_template(study, analysis, folders, args.group_by)

    calibration: dict = {}
    if args.auto_calibrate and not args.report_only:
        folders, calibration = run_auto_calibration(study, analysis, folders, sequence,
                                                     metadata, args.group_by, extra)

    if not args.report_only:
        for polarity, folder in folders.items():
            per_polarity = dict(extra)
            offset = getattr(args, f"mass_offset_{polarity.lower()}")
            if offset is not None:
                per_polarity["mass_offset_ppm"] = offset
            config = build_config(study, analysis, polarity, folder, sequence,
                                  metadata.get(polarity, ""), args.group_by, per_polarity)
            log = analysis / f"run_{polarity}.log"
            print(f"\nrunning {polarity} -> {log}")
            code = run(config, log)
            print(f"  exit {code}")
            if code != 0:
                print(f"  ⚠ {polarity} failed; see {log}", file=sys.stderr)
                return code
            # Written here, not inside `run_auto_calibration`, because `Results/<polarity>/` does
            # not exist until this pass -- the pass whose folder this is (calibrated or original)
            # -- has actually run. `checked: False` distinguishes "never checked" from "checked,
            # clean" for the QC report.
            if args.auto_calibrate:
                cal = calibration.get(polarity, {})
                (analysis / "Results" / polarity / "Calibration.json").write_text(json.dumps(
                    {"checked": True, "used_calibrated_files": bool(cal.get("corrected_files")),
                     **cal}, indent=2))

    # Between the runs and the report, because it needs both inventories to know which classes are
    # ambiguous, and the report must describe the tables as delivered.
    # ⚠ Applied once per set of result tables, and keyed on whether it already has been rather
    # than on which flags were passed. It is not idempotent in the way that matters: run a second
    # time it finds nothing left to decide — the classes it removed are no longer in both
    # polarities — and would overwrite the summary with zeros, so a rebuilt report would state that
    # the filter did nothing, having deleted the evidence that it did.
    #
    # Keying on the record rather than on `--report-only` is what makes both orders work: a report
    # rebuild straight after a fresh pipeline run still applies it, and a second rebuild does not.
    marker = analysis / "polarity_filter.json"
    already = marker.exists()
    # ⚠ ...but only if it still describes the tables that are there. A full pipeline run rewrites
    # `Final_Results_Filtered.csv` from scratch, unfiltered, while leaving this record in place —
    # so the guard refused, and the study shipped BOTH polarities' shared classes: 853 rows where
    # 689 were delivered, and 390 where 312 were. It said so in the log and the tables were wrong
    # anyway, which is the whole problem with a record that outlives what it records.
    if already:
        tables = list(analysis.glob("Results/*/Final_Results_Filtered.csv"))
        newest = max((t.stat().st_mtime for t in tables), default=0.0)
        if newest > marker.stat().st_mtime:
            print("\npolarity filter: the results were rebuilt after it last ran — reapplying")
            marker.unlink()
            already = False
    if not already and not args.no_polarity_filter and len(folders) > 1:
        from lipidloop import polarity as polarity_filter
        default = Path(args.class_polarity) if args.class_polarity else ROOT / "data/class_polarity.csv"
        override = analysis / "class_polarity.csv"
        table = polarity_filter.load(default, override if override.exists() else None)
        print(f"\npolarity filter: {len(table)} classes in the table"
              + (f", {override.name} layered on top" if override.exists() else ""))
        try:
            outcome = polarity_filter.apply(analysis, table, log=print,
                                            metadata=metadata, group_by=args.group_by)
            try:
                polarity_filter.combine(
                    analysis, (outcome["summary"].get("agreement") or {}).get("by_lipid"),
                    log=print)
            except ValueError as exc:
                # ⚠ Combining the two polarities into one table is a CONVENIENCE. It requires the
                # two runs to have measured the same vials, which is not always true and is not
                # this pipeline's business to insist on: one public study names its polarities
                # `20210303P_...` and `20210303N_...` — a letter inside a date token, not a
                # suffix — and has 16 positive injections against 12 negative.
                #
                # Raising here killed a run whose per-polarity results were already written and
                # perfectly valid, and returned a non-zero exit that read as "the analysis
                # failed". Both polarities are complete; only the joined table is missing.
                print(f"\n  ⚠ combined table not written — {exc}", file=sys.stderr)
                print("    Both polarities completed; only the cross-polarity join was skipped.",
                      file=sys.stderr)
        except polarity_filter.UndecidedClass as exc:
            print(f"\n  REFUSED: {exc}", file=sys.stderr)
            return 3

    elif already:
        print("\npolarity filter: already applied to these tables — "
              "delete polarity_filter.json to run it again")

    # Method-ready mass lists, once per study after both polarities: the exclusion list checked
    # against this run's identifications and the searched libraries (colliding entries dropped),
    # the inclusion list in the same Xcalibur layout, both polarities in one folder. Written
    # before the report so a rebuilt report can point at them; never fatal, the tables are done.
    print("\nwriting method lists")
    try:
        from lipidloop.method_lists import write_method_lists
        write_method_lists(analysis, say=lambda msg: print(f"  {msg}"))
    except Exception as exc:                                     # noqa: BLE001
        print(f"  ⚠ method lists not written — {exc}", file=sys.stderr)

    print("\nbuilding the quality report")
    args.study = args.study_name or (parsed.study_name if parsed else "") or study.name
    if not args.owner and parsed and parsed.pi:
        args.owner = parsed.pi
    return build_report(analysis, folders, sequence, metadata, args)


if __name__ == "__main__":
    raise SystemExit(main())
