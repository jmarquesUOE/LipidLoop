"""End to end: `.raw` in, `Final_Results.csv` out, no GUI anywhere.

    raw   -> convert.py    ThermoRawFileParser
    mzML  -> features.py   pyOpenMS: detect, align, link
          -> search.py     library search over the MS2 scans in the same mzML
          -> peakfinder.py join identifications to features, filter, write

Every setting the run depends on lives in one `RunConfig` and is written next to the results.
That is not tidiness for its own sake: library load order alone changes the output, and LipiDex
records it nowhere, so a result file that does not carry its configuration cannot be reproduced.
"""
from __future__ import annotations

import json
from statistics import median
from dataclasses import asdict, dataclass, field, fields, replace
from pathlib import Path

from .adduct_pairs import AdductPairFilter, remove_adduct_pairs
from .artefacts import screen, strip
from .associated import write_associated_spectra
from .components import write_spectral_components
from .blanks import (BlankFilter, apply_blank_filter, carryover_report, infer_roles,
                     read_sample_types, read_sequence)
from .exclusion import duty_cycle, find_wasted, series_breakdown, write_exclusion_list, write_thermo_list
from .inclusion import find_inclusion_candidates, write_inclusion_list
from .calibrate import (MAX_PEAKS_PER_SCAN, MIN_FEATURES_PER_FILE, Calibration,
                        check_noise_threshold, looks_like_dia, measure_noise_floor,
                        ms2_absent, precursor_selection_missing, profile_summary,
                        searchable,
                        threshold_for_feature_survival, threshold_for_peak_cap)
from .correlate import CorrelationFilter
from .convert import Converter, needs_conversion
from .features import FeatureParams, run_feature_detection
from .mgf import read_mzml_ms2
from .peaks import sum_composition
from .peakfinder import PeakFinder, write_results
from .purity import read_fatty_acids
from .fatty_acids import annotate as annotate_fatty_acids
from .fatty_acids import read_standards
from .rtls import extend_identifications
from .gapfill import GapFilling, count_absences_with_signal, fill_gaps
from .presence import PresenceFilter, apply_presence_filter, presence_cells
from .split_peaks import SplitPeakMerge, merge_split_peaks
from .standards import StandardMix, check_offset, measure as measure_standards
from .standards import offset_ppm as standards_offset_ppm
from .standards import write as write_standards
from .adducts import load_adducts
from .search import (_potential_fragments, calc_purity, load_libraries,
                     search_spectrum, sn_evidence_cell, _jstr, _round)


@dataclass
class RunConfig:
    """Everything that changes the answer, in one place."""

    raw_files: list[str] = field(default_factory=list)
    libraries: list[str] = field(default_factory=list)   # ORDER MATTERS: it breaks score ties
    # Free-fatty-acid spectra, kept as their own list rather than appended to `libraries`. Only
    # the polyunsaturates are in it — the rest of the class does not fragment — so this list is
    # the small MS2-confirmable part of a class that is otherwise named by retention. Swap it for
    # a library built from authentic standards when you have one.
    fatty_acid_libraries: list[str] = field(default_factory=list)
    # CSV of authentic-standard retention times, `name,retention` (e.g. `FA 18:1,10.36`). When
    # given, the standards anchor the retention surface and are never trimmed out of the fit,
    # while observed candidates still extend it — the hybrid model. Empty means the surface is
    # fitted from the observed series alone.
    fatty_acid_standards: str = ""
    fatty_acids: str = ""
    # Folders holding Adducts.csv — the adduct sweep is inert without them.
    adduct_libraries: list[str] = field(default_factory=list)
    # Where the result tables land. The delivered convention is
    #
    #     Analysis_<date>/Results/<polarity>/     the four CSVs plus search/
    #     Analysis_<date>/Report_QC/              the quality report
    #
    # so a study folder holds one Analysis_<date> per reprocessing and the results of the two
    # polarities sit side by side under one parent. `standard_output_dir` builds it; a config may
    # still name any directory outright, because reprocessing someone else's layout is normal.
    output_dir: str = "outputs"
    mzml_dir: str = "data/mzml"
    thermo_parser: str = "data/tools/ThermoRawFileParser"

    ms1_tol: float = 0.01
    ms2_tol: float = 0.01
    min_ms2_mass: float = 61.0
    # Artefact handling. `screen_artefacts` detects them per file on every run and is ON:
    # a hard-coded m/z would be right for one instrument in one period and useless to anyone
    # running this method in another lab, so the detector travels with the method rather than
    # a list of masses. Anything it removes is reported. `exclude_mz` adds masses by hand on
    # top, for artefacts that are real ions and so cannot be detected by implausible mass.
    screen_artefacts: bool = True
    exclude_mz: list = field(default_factory=list)
    # MS1->MS2 duty-cycle accounting and exclusion-list candidates, computed from the same MS2
    # scans and identification keys the search loop already reads -- no second mzML parse. ON by
    # default, same reasoning as `screen_artefacts`: a hard-coded exclusion list is right for one
    # method in one period, so this travels with the method instead. Writes
    # `Exclusion_List.csv` (review form, the evidence per candidate) and
    # `Exclusion_List_Thermo.csv` (Xcalibur mass-list import form) beside the results, and a
    # `Duty_Cycle.json` summary the QC report reads. See exclusion.py; nothing is excluded on an
    # instrument by this alone -- run `scripts/check_exclusion_safety.py` before importing a list.
    build_exclusion_list: bool = True
    # The complementary list: MS1 features THIS run detected, low abundance, matching this
    # study's own libraries, that never got an MS2 attempt at all -- forced-fragmentation targets
    # for the study's next injection, not identifications for this one. ON by default for the
    # same reason as `build_exclusion_list`: it is computed from data the pipeline already has
    # (feature areas, quant ions, the library index), so leaving it off would only be throwing
    # away a free-to-compute recommendation. Needs the full feature table, so it runs after peak
    # finding and every filter (stage 6d), not alongside the exclusion-list step, which only has
    # MS2-scan data at the point it runs. Writes `Inclusion_List.csv` (m/z, retention) beside the
    # results. See inclusion.py; nothing is included on an instrument by this alone.
    build_inclusion_list: bool = True
    # CSV of spiked standards: `name,formula,adduct,polarity`. Masses are computed from the
    # formulae, never read from a column — the supplier sheet this was built against listed two
    # deuterated fatty acids at their unlabelled masses, wrong by 31 and 35 Da. Declaring the mix
    # also EXEMPTS it from the artefact screen, which is not optional for deuterated standards:
    # the screen rejects mass defects impossible for CHNOPS, deuterium is not CHNOPS, and d31/d35
    # acids fall inside the rejection region. See standards.py.
    standards: str = ""
    # Lipid classes dropped from the library before searching. d5TG is 13% of
    # LipiDex_HCD_Formic and is a deuterated standard set — remove this entry if you spike
    # d5 standards deliberately. Recorded in run_config.json like everything else.
    excluded_classes: list[str] = field(default_factory=lambda: ["d5TG"])
    min_feature_count: int = 2
    rt_filter: bool = True
    rt_filter_multiplier: float = 2.0
    adduct_filtering: bool = True
    in_source_filtering: bool = True
    # Retention modelling within each lipid class (LipiDex 2 calls it RTLS).
    retention_model_filter: bool = True
    retention_model_sigma: float = 3.0
    retention_model_max_error: float = 1.0
    # Systematic mass error, in ppm, applied to centre the search window where the ions
    # actually are. Derived by `calibrate.py`; 0.0 means uncorrected.
    mass_offset_ppm: float = 0.0
    # Naming unidentified features from the model. OFF by default: an assignment with no MS2
    # behind it is weaker evidence than an identification and should be opted into.
    retention_model_extend: bool = False
    retention_model_extend_sigma: float = 2.0
    # Mass window for naming a feature from the model. The default was 20 ppm, which is
    # four times the instrument's measured accuracy on these runs (median error
    # -1.78 ppm, 78% of matches inside 5 ppm) and admits candidates the mass evidence
    # does not support. 10 ppm is generous against a calibrated Orbitrap.
    retention_model_extend_ppm: float = 10.0
    retention_model_extend_max_error: float = 0.5
    # Adds MS2 Spectra / MS2 Files to the output. Extends LipiDex's schema by two columns;
    # set false for exactly the columns LipiDex writes.
    write_ms2_support: bool = True
    # Dot product, reverse dot product and purity for the name as reported. On by default: a
    # result table without them asks the reader to trust every identified row equally, and rows
    # run from a dot product near 1000 down to below 500.
    write_scores: bool = True
    # Write every compound group to Final_Results.csv, with `Filter Status` saying what would
    # have been dropped, instead of dropping it. Off by default so the headline figures stay
    # comparable with LipiDex's own output — but see docs/FILTERS.md: the class retention-window
    # filter removes abundant long-chain species by construction, so `true` is the honest setting
    # when the long-chain end of a class matters.
    keep_filtered: bool = False
    # Name free fatty acids from mass + a fitted RT ~ C + DB surface. On because no library can
    # reach them — the default set has zero negative-mode entries across 91 fatty-acid masses,
    # and the carboxylate anion does not fragment. Nothing is named unless the whole candidate
    # set behaves like a homologous series; see fatty_acids.py.
    annotate_fatty_acids: bool = True
    # `noise_threshold` is an ABSOLUTE INTENSITY and does not transfer between detectors. The floor
    # is measured on every run and reported whatever this says; setting it true also *applies* the
    # derived value, in the same pass — it needs no second one, being read from raw MS1.
    # "auto" | True | False.
    #
    # ⚠ `noise_threshold` is an ABSOLUTE INTENSITY, so it cannot be portable. 5,000 counts was
    # fitted on the facility Orbitrap; on an Agilent 6545 QTOF the MS1 noise floor is 214 counts,
    # making that default 23x the floor. On a real dataset it cut negative-mode features from 590
    # to 170 per file and took recall against a published lipid list from 67.8% to 50.9% — with no
    # error, because the failure mode of a too-high threshold is simply absence.
    #
    # "auto" keeps the configured value when it sits in the plausible band (2-15x the measured
    # floor) and derives when it does not. So a facility run, where the fitted value is correct by
    # construction, is bit-for-bit unchanged and reproducible; a run on someone else's instrument
    # adapts instead of silently under-reporting. True and False force the behaviour.
    derive_noise_threshold: "bool | str" = "auto"
    correlation: CorrelationFilter = field(default_factory=CorrelationFilter)
    features: FeatureParams = field(default_factory=FeatureParams)

    # Roles are inferred from file names when not given; always overridable, because a
    # blank that is not named like one would otherwise be filtered against itself.
    sample_roles: dict = field(default_factory=dict)
    sequence_file: str = ""      # Xcalibur sequence CSV, for injection order
    blank_filter: BlankFilter = field(default_factory=BlankFilter)
    split_peaks: SplitPeakMerge = field(default_factory=SplitPeakMerge)
    adduct_pairs: AdductPairFilter = field(default_factory=AdductPairFilter)
    presence_filter: PresenceFilter = field(default_factory=PresenceFilter)
    gap_filling: GapFilling = field(default_factory=GapFilling)
    calibration: Calibration = field(default_factory=Calibration)
    # Filled in by a calibrated run: what was measured, what it replaced, which pass.
    calibration_report: dict = field(default_factory=dict)
    # Written by every run, whether or not the threshold was derived from it, so a result carries
    # the evidence for whether its noise threshold suited the instrument that produced it.
    measured_noise_floor: float = 0.0

    @classmethod
    def load(cls, path: str | Path) -> "RunConfig":
        data = json.loads(Path(path).read_text())
        feature_params = FeatureParams(**data.pop("features", {}))
        blank = BlankFilter(**data.pop("blank_filter", {}))
        split = SplitPeakMerge(**data.pop("split_peaks", {}))
        pairs = AdductPairFilter(**data.pop("adduct_pairs", {}))
        presence = PresenceFilter(**data.pop("presence_filter", {}))
        gaps = GapFilling(**data.pop("gap_filling", {}))
        calib = Calibration(**data.pop("calibration", {}))
        correlation = CorrelationFilter(**data.pop("correlation", {}))
        # ⚠ Not every key a config carries configures the run. `run_study.py` also writes
        # audit-only fields onto it (e.g. `modifier_detected`, `modifier_evidence` from
        # --detect-modifier) so a saved config explains itself without a second file. Passing
        # `**data` straight into the constructor made every such addition a breaking change:
        # add one informational field to the writer and every saved config from before it existed
        # (or every field the constructor does not declare) fails to load. Keep only what this
        # dataclass actually declares; the rest was worth writing for a human to read, not for
        # this constructor to consume.
        known = {f.name for f in fields(cls)}
        data = {k: v for k, v in data.items() if k in known}
        return cls(features=feature_params, blank_filter=blank, correlation=correlation,
                   split_peaks=split, adduct_pairs=pairs, presence_filter=presence,
                   gap_filling=gaps, calibration=calib, **data)

    def save(self, path: str | Path) -> None:
        data = asdict(self)
        Path(path).write_text(json.dumps(data, indent=2))


@dataclass
class PipelineOutput:
    """What a run produced, for a caller that wants more than the file path."""
    final_results: Path
    result: object          # PeakFinderResult
    samples: list
    config: "RunConfig"
    injection_order: dict
    # Per-class retention models, exposed because a calibration pass takes its retention tolerance
    # from their fitted residuals. Without this the estimator has nothing to read and skips.
    retention_models: dict = field(default_factory=dict)


def polarity_label(out) -> str:
    """`Pos`/`Neg` taken from the output folder's name, which is what the layout calls it."""
    return Path(out).name or "run"


def standard_output_dir(analysis_dir: str | Path, polarity: str) -> Path:
    """`Analysis_<date>/Results/<polarity>` — the delivered folder convention.

    Kept in one place so a new study does not reinvent it, and so the quality report and the
    pipeline cannot drift apart about where results live. `report_dir` is its counterpart.
    """
    return Path(analysis_dir) / "Results" / polarity


def report_dir(analysis_dir: str | Path) -> Path:
    """`Analysis_<date>/Report_QC` — one report per analysis, covering both polarities."""
    return Path(analysis_dir) / "Report_QC"


def run(config: RunConfig, log=None) -> PipelineOutput:
    say = log or (lambda _: None)
    out = Path(config.output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # 1 ─ conversion
    # ⚠ Locate the converter ONLY if something actually needs converting. `Converter.find` raises
    # when ThermoRawFileParser is absent, so finding it up front made a study of pure mzML fail on
    # a missing Thermo reader it never needed — which is most public data and any collaborator who
    # sends mzML.
    to_convert = [r for r in config.raw_files if needs_conversion(r)]
    if to_convert:
        converter = Converter.find(Path(config.thermo_parser))
        mzml_files = converter.convert_all(config.raw_files, config.mzml_dir, log=say)
    else:
        mzml_files = [Path(r) for r in config.raw_files]
        say(f"{len(mzml_files)} files already mzML — no conversion needed")
    say(f"{len(mzml_files)} mzML files ready")

    # 2 ─ identification.
    #
    # The acquisition is checked BEFORE the search, not after. On an all-ion deposit the search
    # runs to completion, finds nothing, and exits 0 — 32 minutes to reach a conclusion that the
    # first file answers in a second, and an emptiness indistinguishable from a library gap.
    for note in (profile_summary(mzml_files), ms2_absent(mzml_files),
                 precursor_selection_missing(mzml_files)):
        if note:
            say(note)

    # ⚠ Acquisition is a property of the FILE, so the verdict is per file and the search runs only
    # on the ones it can work on. Deposits mix modes within a single submission, and the mix is not
    # visible from any metadata field.
    #
    # ST000991 is the case: 153 SWATH files, 5 DDA, 1 MS1-only. Searching the SWATH files produced
    # 21,345 hits drawn from 5 distinct precursor masses — window centres landing within a few mDa
    # of a library entry, rediscovered in every file. Every one was discarded downstream, so the
    # cost was hours rather than a wrong answer, and the study's real identifications came from the
    # 5 DDA files.
    #
    # Skipped files are excluded from the SEARCH only. They keep their MS1, and feature detection,
    # alignment and quantitation below run on all of them.
    search_files, skipped = searchable(mzml_files)
    if skipped:
        for reason, files in sorted(skipped.items()):
            say(f"⚠ {len(files)} of {len(mzml_files)} file(s) are {reason} and cannot be searched "
                f"against a spectral library — excluded from identification, kept for feature "
                f"detection and quantitation. First: {Path(files[0]).name}")
        if not search_files:
            say("⚠ no file in this study can be searched; identifications can still come from the "
                "retention model where a class has anchors.")

    index = load_libraries([Path(p) for p in config.libraries + config.fatty_acid_libraries],
                           ms1_tol=config.ms1_tol,
                           excluded_classes=config.excluded_classes)
    say(f"library: {len(index.spectra)} spectra"
        + (f" (excluding {', '.join(config.excluded_classes)})" if config.excluded_classes else ""))
    fa_db = read_fatty_acids(config.fatty_acids) if config.fatty_acids else []

    mix = StandardMix.from_csv(config.standards) if config.standards else None
    protected: list[float] = []
    if mix:
        # Every declared mass, both polarities. Guessing the run's polarity from a folder name
        # would be one rename away from silently deleting a standard, and protecting an ion that
        # cannot appear in this polarity costs nothing.
        protected = mix.protected_mz()
        say(f"standards: {len(mix.standards)} declared, exempt from the artefact screen")

    # Roles, moved ahead of identification (originally computed just before feature detection)
    # so the exclusion-list step below knows which files are blanks without a second pass.
    sample_types = read_sample_types(config.sequence_file) if config.sequence_file else {}
    roles = dict(infer_roles([Path(p).stem for p in mzml_files], sample_types))
    roles.update(config.sample_roles)

    result_files: dict[str, Path] = {}
    artefacts_seen: dict[str, list] = {}
    # Every MS2 scan in the batch, `(key, precursor_mz, retention_minutes)`, and which keys were
    # identified -- the exclusion-list step's own input, built as a byproduct of the loop below
    # rather than a second mzML parse.
    all_scans: list[tuple[tuple[str, int], float, float]] = []
    identified_keys: set[tuple[str, int]] = set()
    # Identified peaks, for the targeted feature-recovery pass. Identification runs before
    # feature detection precisely so this is available.
    targets: dict[str, list[tuple[float, float]]] = {}
    for mzml in search_files:
        rows = []
        stem = Path(mzml).stem
        spectra = list(read_mzml_ms2(mzml, min_ms2_mass=config.min_ms2_mass,
                                     exclude_mz=config.exclude_mz))
        if config.screen_artefacts:
            found = screen(spectra, protected=protected)
            if found:
                removed = strip(spectra, found)
                artefacts_seen[stem] = {"artefacts": found, "peaks_removed": removed}
                for artefact in found:
                    say(f"  artefact in {stem}: {artefact} — "
                        f"{removed} peaks removed")
        # A measured precursor carries the instrument's systematic mass error. Correcting it
        # here centres the search window on where the ions actually are, rather than widening the
        # window — which would cost specificity to buy the same tolerance.
        if config.mass_offset_ppm:
            for ms2 in spectra:
                ms2.precursor -= ms2.precursor * config.mass_offset_ppm * 1e-6
        for ms2 in spectra:
            key = (stem, ms2.number)
            all_scans.append((key, ms2.precursor, ms2.retention))
            hits = search_spectrum(ms2, index, ms2_tol=config.ms2_tol)
            if not hits:
                continue
            identified_keys.add(key)
            purity, calculator = calc_purity(hits, ms2, fa_db, mz_tol=config.ms2_tol)
            top = hits[0]
            lib = top.library_spectrum
            components = " / ".join(f"{l.name.replace(';', '')}({p})"
                                    for l, p in zip(calculator.lipids, calculator.purities))
            rows.append({
                "MS2 ID": ms2.number, "Retention Time (min)": ms2.retention, "Rank": 1,
                "Identification": lib.name, "Precursor Mass": ms2.precursor,
                "Library Mass": lib.precursor_mz, "Delta m/z": top.delta_ppm,
                "Dot Product": int(_round(top.dot)),
                "Reverse Dot Product": int(_round(top.reverse_dot)),
                "Purity": purity, "Spectral Components": components,
                "Optimal Polarity": str(lib.optimal_polarity).lower(),
                "LipiDex Spectrum": str(lib.is_lipidex).lower(),
                "Library": lib.library,
            })
            # The in-source-fragment filter matches another group's quant ion against this
            # lipid's predicted fragments, so leaving this column blank disables that filter.
            prefix = ",".join(str(v) for v in rows[-1].values())
            rows[-1]["Potential Fragments"] = _potential_fragments(
                purity, calculator, lib, prefix)
            # ⚠ AFTER the line above, never before. `prefix` is a join over this dict's values, and
            # `_potential_fragments` suppresses any mass already present in it as a substring. A
            # column added earlier would land in that prefix and change which masses are written,
            # which changes the in-source-fragment filter — a new evidence column silently editing
            # the delivered table.
            rows[-1]["sn Evidence"] = sn_evidence_cell(lib, ms2, ms2_tol=config.ms2_tol)
        targets[Path(mzml).stem] = [
            (r["Precursor Mass"], r["Retention Time (min)"] * 60.0) for r in rows
            if r["Dot Product"] > 500 and r["Reverse Dot Product"] > 700]
        # Into `search/`, not the run folder. These are the peak finder's input, one per
        # injection, and at 50-plus files a polarity they buried the four tables a reader
        # actually wants. `Associated_Spectra.csv` is written at the top from exactly these.
        search_dir = out / "search"
        search_dir.mkdir(parents=True, exist_ok=True)
        path = search_dir / f"{Path(mzml).stem}_Results.csv"
        _write_rows(path, rows)
        result_files[Path(mzml).stem] = path
        say(f"{Path(mzml).stem}: {len(rows)} identifications")

    # 1a ─ artefact-screening summary. Written whenever screening ran, whether or not anything was
    # found -- otherwise the QC report cannot tell "checked, clean" from "never checked" (the file
    # predates `screen_artefacts`, or it was turned off). `artefacts_seen` was already computed as
    # a byproduct of the search loop above; this is the only place it is written anywhere.
    if config.screen_artefacts:
        (out / "Artefacts.json").write_text(json.dumps({
            "files_screened": len(search_files),
            "files_with_artefacts": len(artefacts_seen),
            "by_file": {
                stem: {
                    "peaks_removed": info["peaks_removed"],
                    "artefacts": [{"mz": a.mz, "fraction": a.fraction,
                                  "base_fraction": a.base_fraction, "spread_mda": a.spread_mda}
                                 for a in info["artefacts"]],
                } for stem, info in artefacts_seen.items()},
        }, indent=2))
        say(f"wrote Artefacts.json: {len(artefacts_seen)} of {len(search_files)} file(s) "
            f"carried an artefact")

    # 2a ─ MS1->MS2 duty cycle and exclusion-list candidates, from the scans and identification
    # keys just built above. A precursor earns a place on the list only by chromatography (never
    # identified, fragmented enough times to matter, present across most of the gradient, and in
    # a blank when one exists) -- never by "unidentified" alone, which would just as happily flag
    # a real lipid the library does not cover. See exclusion.py for the full reasoning.
    if config.build_exclusion_list and all_scans:
        blank_stems = {stem for stem, role in roles.items() if role == "blank"}
        blank_mz = sorted(mz for (stem, _), mz, _rt in all_scans if stem in blank_stems)
        wasted = find_wasted(all_scans, identified_keys, blank_precursors=blank_mz or None)
        total_ms2 = len(all_scans)
        cycle = duty_cycle(wasted, total_ms2)
        say(f"MS1→MS2: {total_ms2:,} MS2 scans over {len(search_files)} file(s); "
            f"{len(wasted):,} unidentified precursors consuming {cycle:.1%} of the duty cycle")
        pol_sign = "+" if polarity_label(out) == "Pos" else "-"
        write_exclusion_list(wasted, out / "Exclusion_List.csv", pol_sign)
        write_thermo_list(wasted, out / "Exclusion_List_Thermo.csv", pol_sign)
        breakdown = series_breakdown(wasted, total_ms2)
        (out / "Duty_Cycle.json").write_text(json.dumps({
            "total_ms2": total_ms2,
            "files_searched": len(search_files),
            "wasted_precursors": len(wasted),
            "wasted_scans": sum(w.scans for w in wasted),
            "duty_cycle_wasted_fraction": cycle,
            "by_series": [{"series": name, "precursors": count, "scans": scans,
                          "fraction": fraction} for name, count, scans, fraction in breakdown],
        }, indent=2))
        say(f"wrote Exclusion_List.csv / Exclusion_List_Thermo.csv / Duty_Cycle.json: "
            f"{len(wasted)} candidates")

    # 2b ─ the noise floor. Always measured, because the failure mode of an absolute intensity
    # threshold on a new detector is silence: too high and features are simply absent, too low and
    # the run does not finish. Applied only when asked.
    dia = looks_like_dia(mzml_files)
    if dia:
        say(dia)

    floor = measure_noise_floor(mzml_files, roles=roles)
    derived_noise, message, plausible = check_noise_threshold(
        config.features.noise_threshold, floor)
    say(message)

    # ⚠ A multiple of the floor is not portable, and trusting it alone broke a run.
    #
    # The floor is measured correctly everywhere, but `floor x 5.5` assumes an Orbitrap-shaped
    # intensity distribution. On Waters data — 10-46x more peaks, a long tail the vendor's
    # centroiding leaves in — it derived 88, retained a third of all peaks, and feature detection
    # produced 225,744 compound groups from 106 injections, burying the identifications.
    #
    # So the derived value is floored by a scale-free second estimate: the threshold that leaves at
    # most MAX_PEAKS_PER_SCAN standing. Whichever is HIGHER wins, because the failure being guarded
    # is a threshold far too low. On the instruments where the multiplier already worked this
    # changes nothing — the cap is the looser of the two there.
    capped = threshold_for_peak_cap(mzml_files, roles=roles)
    if capped and derived_noise and capped > derived_noise:
        say(f"⚠ the derived threshold {derived_noise:,.0f} would leave more than "
            f"{MAX_PEAKS_PER_SCAN} MS1 peaks per scan standing — raising it to {capped:,.0f}, "
            f"which is where this data's peak density says signal starts. A threshold below that "
            f"turns noise into features.")
        derived_noise = capped

    # ⚠ And the same multiplier fails in the OTHER direction, which the cap above cannot catch.
    #
    # On a low-count TOF (Leco Citius: median MS1 peak intensity 37) `floor x 5.5` derived 103 — a
    # value a real lipid's monoisotopic peak clears but its 13C isotope does not. Since a feature
    # needs two traces, every trace was discarded for having no isotope partner: 1,065 traces
    # found, 0 features built, all 153 injections at a median of 2 features, and the study died in
    # alignment several stages downstream of the cause.
    #
    # Predicting this from a proxy does not work — retained-peak fraction rates the failing run
    # (8.52%) above a healthy Bruker one (8.24%). So the ceiling measures the real quantity: detect
    # features on the richest sampled file and step the threshold down until they stop starving.
    # Costs one detection run on instruments that never had the problem, and changes nothing there.
    if derived_noise:
        survivable = threshold_for_feature_survival(mzml_files, derived_noise, floor, roles=roles)
        if survivable and survivable < derived_noise:
            say(f"⚠ the threshold {derived_noise:,.0f} leaves fewer than {MIN_FEATURES_PER_FILE} "
                f"features in the richest injection — lowering it to {survivable:,.0f}. Above that "
                f"this detector's counts are so low that the 13C isotope falls below the cut, and "
                f"a feature without its isotope partner is discarded.")
            derived_noise = survivable
    mode = config.derive_noise_threshold
    if mode == "auto":
        # Derive only when the configured value cannot belong to this detector. Keeping it when it
        # is plausible is what makes "auto" safe to default: the facility's own runs do not move.
        use_derived = bool(floor) and not plausible
    else:
        use_derived = bool(mode) and bool(floor)
    if use_derived:
        config = replace(config, features=replace(config.features,
                                                  noise_threshold=float(derived_noise)))
        say(f"noise_threshold derived from this run: {derived_noise:,.0f}"
            + ("  (configured value was implausible for this detector)" if mode == "auto" else ""))
    elif mode == "auto" and floor:
        say(f"noise_threshold {config.features.noise_threshold:,.0f} kept — "
            f"consistent with this detector's floor")
    config.measured_noise_floor = float(floor) if floor else 0.0

    # 3 ─ features
    # Sample Type from the sequence, where there is one, alongside the file names. Two sources
    # because they fail in different situations — a sequence may predate the naming convention,
    # a file may be renamed after acquisition — and either may declare a role. (Roles themselves
    # are computed earlier, ahead of identification, for the exclusion-list step.)
    counts = {r: sum(1 for v in roles.values() if v == r) for r in set(roles.values())}
    say(f"roles: {counts}")

    groups, samples = run_feature_detection(mzml_files, config.features,
                                            min_feature_count=config.min_feature_count,
                                            log=say, roles=roles, targets=targets)

    order = read_sequence(config.sequence_file) if config.sequence_file else {}
    for sample in samples:
        sample.injection = order.get(sample.file)

    # 4 ─ peak finder
    adducts = load_adducts(config.adduct_libraries)
    if adducts:
        say(f"adducts: {len(adducts)}")
    finder = PeakFinder(groups, samples, adducts=adducts,
                        retention_model_filter=config.retention_model_filter,
                        retention_model_sigma=config.retention_model_sigma,
                        retention_model_max_error=config.retention_model_max_error,
                        correlation=config.correlation,
                        rt_filter=config.rt_filter,
                        rt_filter_multiplier=config.rt_filter_multiplier,
                        adduct_filtering=config.adduct_filtering,
                        in_source_filtering=config.in_source_filtering,
                        min_feature_count=config.min_feature_count)
    result = finder.run({s.file: result_files[s.file] for s in samples if s.file in result_files},
                        log=say)

    # 4a ─ the same systematic correction, now on the feature masses. Applied here rather than at
    # detection so the peak finder's own comparisons — adduct spacings, isotope spacings — are made
    # on the measured values, where a uniform shift cancels anyway.
    #
    # This is what puts the retention-model extension on the corrected axis: it matches a feature's
    # `quant_ion` against a library mass in a window centred on zero, so an uncorrected offset
    # spends half that window before the search begins. The MS2 route was already corrected, at the
    # precursor, before the library lookup.
    #
    # The measured value is kept on every group and written beside the corrected one.
    if config.mass_offset_ppm:
        shifted = 0
        for group in result.compound_groups:
            if group.quant_ion:
                group.quant_ion_measured = group.quant_ion
                group.quant_ion -= group.quant_ion * config.mass_offset_ppm * 1e-6
                shifted += 1
        say(f"applied {config.mass_offset_ppm:+.2f} ppm to {shifted:,} feature masses; the "
            f"measured values are kept in `Quant Ion (measured)`")

    # 4b ─ split peaks. Before the blank filter, so a group is judged on the whole of its signal
    # rather than on whichever half of it happened to be the larger row.
    if config.split_peaks.enabled:
        merge_split_peaks(result.compound_groups, samples, config.split_peaks,
                          finder.avg_fwhm, log=say)

    # 4c ─ two identified rows that are two adducts of one molecule. After the split merge, so
    # a molecule split across rows is one row before being compared with anything else, and
    # before the blank filter so a ghost is never the thing that survives it.
    if config.adduct_pairs.enabled:
        remove_adduct_pairs(result.compound_groups, samples, config.adduct_pairs,
                            finder.avg_fwhm, adducts=adducts, log=say)

    # 5 ─ blanks. After the peak finder, so a group is judged on its final areas.
    if config.blank_filter.enabled:
        for row in carryover_report(result.kept, samples, order):
            say(f"blank {row['file']} (inj {row['injection']}): "
                f"{row['groups detected']} groups, "
                f"{row['fraction shared with samples']:.0%} of its signal shared with samples")
        report = apply_blank_filter(result.compound_groups, samples, config.blank_filter)
        say(str(report))

    # 6a ─ free fatty acids, which no spectral library can reach. They do not fragment: 711 MS2
    # spectra of FA 16:0 across 16 files contain no reproducible product ion. What identifies
    # them is that thirty-odd features on exact fatty-acid masses jointly obey RT ~ C + DB.
    # After the blank filter so a contaminant series cannot be named, before the general
    # retention-model extension so those names are available to it.
    if config.annotate_fatty_acids:
        annotate_fatty_acids(result.compound_groups, polarity="-",
                             standards=read_standards(config.fatty_acid_standards), log=say)

    # 6b ─ extend identifications from the retention model, if asked
    if config.retention_model_extend and finder.retention_models:
        candidates = {(sum_composition(s.lipid), round(s.precursor_mz, 4),
                       "+" if s.polarity == "positive" else "-", s.adduct)
                      for s in index.spectra}
        assigned = extend_identifications(
            result.compound_groups, finder.retention_models, sorted(candidates),
            mz_tol_ppm=config.retention_model_extend_ppm,
            sigma=config.retention_model_extend_sigma,
            max_rt_error=config.retention_model_extend_max_error)
        say(f"retention model named {assigned} features that had no MS2")

    # 6c ─ presence. Last of the filters, and after the retention-model naming, so a row is
    # judged on the areas it will actually be analysed with. Applied to every compound group
    # rather than only the identified ones: how reliably a feature was detected is a property of
    # the measurement, and an unidentified row surviving on two detections while a named one is
    # dropped on nine would be indefensible.
    if config.presence_filter.enabled:
        apply_presence_filter(result.compound_groups, samples, config.presence_filter, log=say)

    # 6d ─ inclusion-list candidates: real MS1 features this run's own Top-N never gave an MS2
    # to, for this study's NEXT injection rather than this one. Runs here, after every filter and
    # after the retention-model extension (6b), because a feature that step just named is no
    # longer a candidate, and `keep`/`quant_ion` are only final once every earlier stage has had
    # its say. See inclusion.py; nothing here changes this run's identifications.
    if config.build_inclusion_list:
        included = find_inclusion_candidates(result.compound_groups, all_scans, index,
                                             ms1_tol=config.ms1_tol)
        write_inclusion_list(included, out / "Inclusion_List.csv")
        say(f"wrote Inclusion_List.csv: {len(included)} candidates "
            f"for the next injection of this study")

    final = out / "Final_Results.csv"
    # With `keep_filtered`, `Final_Results.csv` carries every compound group and the
    # `Filter Status` column that says what would have been dropped, rather than dropping it.
    # The scores are what make that usable: nothing is hidden, and the reader can see which rows
    # are worth trusting instead of trusting the filter to have been right.
    # `Identification Source` is written whenever any row could have been named without a
    # spectrum. Free-fatty-acid annotation is exactly that case, and a row named from retention
    # alone must say so in the table, not only in the log.
    id_source = config.retention_model_extend or config.annotate_fatty_acids
    # Standards are kept out of every delivered table. They are a facility measurement — the same
    # material injected in every study to judge the instrument — not the client's data, and a
    # column of them sitting beside the samples is one accidental selection away from a fold
    # change. They are reported in full in Standards_Performance.csv.
    write_results(result, final, unfiltered=config.keep_filtered, id_source=id_source,
                  ms2_support=config.write_ms2_support, scores=config.write_scores,
                  drop_roles=("standard",))
    write_results(result, out / "Unfiltered_Results.csv", unfiltered=True, id_source=id_source,
                  drop_roles=("standard",),
                  ms2_support=config.write_ms2_support, scores=config.write_scores)
    # `Final_Results_Filtered.csv` is the analysis-ready table: filtered, and with the gaps
    # measured rather than left blank. A second table, never a replacement — `Final_Results.csv`
    # keeps its zeros, so what was filled can always be recovered by diffing the two.
    if config.presence_filter.enabled and config.gap_filling.enabled:
        # Exactly the groups this table will contain. Filling the whole set instead still wrote
        # the right file — the extra rows are never written — but reported a count an order of
        # magnitude too large, counting rows the reader will never see.
        written = result.compound_groups if config.keep_filtered else result.kept
        original = [list(g.areas) for g in written]
        before = {id(g): a for g, a in zip(written, original)}
        report = fill_gaps(written, samples, {Path(m).stem: m for m in mzml_files},
                           config.gap_filling, finder.avg_fwhm, log=say)
        # The presence filter protects a group with no detections at all, on the reasoning that a
        # whole-group zero is a real absence. That is only true while the zeros are real, and
        # re-integration is the test. Reported because it decides whether those calls stand.
        cells = presence_cells(samples, config.presence_filter)
        report.absences_with_signal, report.absences_tested = count_absences_with_signal(
            written, samples, cells, before)
        if report.absences_tested:
            say(f"  whole-group absences carrying a peak after re-integration: "
                f"{report.absences_with_signal}/{report.absences_tested} — a presence/absence "
                f"call is only as good as the detection behind it")
        write_results(result, out / "Final_Results_Filtered.csv", unfiltered=config.keep_filtered,
                      id_source=id_source, ms2_support=config.write_ms2_support,
                      scores=config.write_scores, identified_only=True,
                      drop_roles=("standard",))
        # The same table, median-normalised. Written beside the unnormalised one rather than
        # instead of it: normalisation is a choice about the analysis, and a delivered table that
        # has silently been scaled cannot be un-scaled by whoever receives it. Blanks come out
        # too — a normalisation factor derived from a column with almost no signal is not a
        # correction, it is noise given authority.
        write_results(result, out / "Final_Results_Filtered_median_Normalised.csv",
                      unfiltered=config.keep_filtered, id_source=id_source,
                      ms2_support=config.write_ms2_support, scores=config.write_scores,
                      identified_only=True, drop_roles=("standard", "blank"), normalise=True)
        say("wrote Final_Results_Filtered_median_Normalised.csv — the same lipids, "
            "median-normalised, blanks excluded")
        for group, areas in zip(written, original):
            group.areas = areas
        named = sum(1 for g in written if (g.identification()[0] or "").strip())
        say(f"wrote Final_Results_Filtered.csv — the analysis-ready table: {named} identified "
            f"lipids of {len(written)} kept rows. Unnamed features are deliberately absent; they "
            f"cannot enter a lipid-level analysis and are kept in Final_Results.csv")

    # Standards performance. Written whenever a mix is declared and standards injections are
    # present — it judges the instrument, not the study, and is the one table comparable with
    # other runs, because the mix is the same material every time.
    if mix:
        polarity = "-" if any(s.quant_polarity == "-" for s in result.kept if s.quant_polarity) \
            else "+"
        # Cross-check the mass calibration. The offset a run should use is derived from the
        # identifications — tens of thousands of matches beat a handful of standards — but that
        # estimate assumes the identifications are right, so it cannot validate itself. A spiked
        # standard's mass is known before the run, which makes it independent evidence rather
        # than a better measurement. Agreement is the point; disagreement means one of them is
        # wrong and should stop the run being trusted.
        deltas = [l.ppm_error for g in result.compound_groups for c in g.lipid_candidates
                  for l in c.identifications if l.ppm_error == l.ppm_error]
        # ⚠ Compare like with like. `mass_offset_ppm` shifts the MS2 precursors before the
        # search, so once a correction is applied the identification-derived error is what is
        # LEFT OVER, while the standards are measured on the feature's quant ion, which is never
        # shifted and so still carries the instrument's full error. Comparing the two directly
        # made a correctly calibrated run report a disagreement exactly the size of the correction
        # — the check failed precisely when it had worked. Adding the applied offset back puts
        # both on the instrument's own scale.
        residual = median(deltas) if deltas else None
        derived = (residual + config.mass_offset_ppm) if residual is not None else None
        from_standards = standards_offset_ppm(mix, polarity, result.compound_groups)
        check = check_offset(derived, from_standards)
        if check["verdict"] == "not comparable":
            say("mass calibration: standards not found in this run, no cross-check possible")
        else:
            applied = (f", after correcting {config.mass_offset_ppm:+.2f} "
                       f"(residual {residual:+.2f})") if config.mass_offset_ppm else ""
            say(f"mass calibration: {derived:+.2f} ppm from {len(deltas):,} identifications"
                f"{applied}, {from_standards:+.2f} ppm from the standards, "
                f"difference {check['difference_ppm']:+.2f} ppm — {check['verdict']}")
            if not check["agree"]:
                say("  ⚠ the two disagree. One of them is wrong: either the identifications are "
                    "not what they claim, or the standards were mis-declared. Do not apply an "
                    "offset until this is resolved.")

        rows = measure_standards(mix, polarity, result.compound_groups, samples)
        if rows:
            # Not into Results/. The standards judge the instrument, not the study, so the table
            # belongs with the facility's own records rather than in the folder that goes to the
            # client. Under the delivered layout that is the analysis root, one level above
            # Results/<polarity>; anywhere else it sits beside the results, because guessing a
            # parent directory in a layout we do not recognise is how a file goes missing.
            standards_dir = out.parent.parent if out.parent.name == "Results" else out
            standards_dir.mkdir(parents=True, exist_ok=True)
            target = standards_dir / f"Standards_Performance_{polarity_label(out)}.csv"
            n = write_standards(rows, target)
            found = sum(1 for r in rows if r.get("found"))
            say(f"wrote {target.name} (facility record, outside Results/): "
                f"{found}/{n} declared standards found")
            for r in rows:
                if not r.get("found"):
                    say(f"  ⚠ standard NOT FOUND: {r['standard']} {r['adduct']} "
                        f"expected m/z {r['expected_mz']}")
        else:
            say("standards declared but no injection carries the `standard` role — "
                "name them Std_Mix_* or set sample_roles in the config")

    spectra_table = out / "Associated_Spectra.csv"
    n_spectra = write_associated_spectra(result, result_files, spectra_table)
    say(f"wrote {spectra_table.name}: {n_spectra} spectral matches over "
        f"{len(result_files)} injections")

    components_table = out / "Spectral_Components.csv"
    n_components, n_divided = write_spectral_components(result, components_table)
    say(f"wrote {components_table.name}: {n_components} components over "
        f"{n_divided} peaks whose chain-bearing signal divides")
    config.save(out / "run_config.json")
    say(f"wrote {final}")
    return PipelineOutput(final_results=final, result=result, samples=samples,
                          config=config, injection_order=order,
                          retention_models=dict(finder.retention_models or {}))


def _write_rows(path: Path, rows: list[dict]) -> None:
    import csv
    from .search import RESULT_COLUMNS
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=RESULT_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)
