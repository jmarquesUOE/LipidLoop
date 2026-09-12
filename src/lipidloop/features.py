"""Feature detection, alignment and linking — the half Compound Discoverer used to do.

pyOpenMS in-process, not MZmine as a subprocess: no batch XML that has to be authored in a GUI
first, no second install, and the parameters live in the same config as the LipiDex constants.

    mass traces          MassTraceDetection      m/z traces through time
      elution peaks      ElutionPeakDetection    split traces into chromatographic peaks
        features         FeatureFindingMetabo    assemble isotope patterns into features
          alignment      MapAlignerPoseClustering  put every file on one retention-time axis
            linking      FeatureGroupingAlgorithmQT  the same feature across files

The output is deliberately the same three-level structure Compound Discoverer exports —
compound group, compound, feature — so `peakfinder.py` cannot tell which produced it.

One caveat this data makes concrete. A DDA run spends most of its duty cycle on MS2: these
files carry 1,838 MS1 scans against 15,019 MS2, about 1.2 Hz, and with peaks around 0.07 min
FWHM that is roughly five MS1 points across a peak. That is thin for any peak picker, and it
is a property of the acquisition rather than of the software — `min_trace_length` has to be set
against it rather than against a textbook value.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as _np

from .peaks import Compound, CompoundGroup, Feature, Sample


@dataclass(slots=True)
class FeatureParams:
    """Detection settings. Defaults tuned against the Compound Discoverer reference run."""

    mass_error_ppm: float = 10.0
    # Tuned with scripts/tune_features.py, which scores a parameter set by how many confident
    # identifications get a feature. That found a much higher-recall setting — noise 10000 with
    # remove_single_traces off lifts detection recall from 74% to 84% — and it was tried and
    # rejected, because end to end it bought 0.4 points of molecule recovery for 2.25x the
    # output rows. The extra features are overwhelmingly unidentified, and the targeted
    # recovery pass already gets the identified peaks for 2.5% more features instead of 170%.
    # Kept here as the reason not to loosen these: see docs/FEATURE_DETECTION.md.
    noise_threshold: float = 5000.0
    chrom_peak_snr: float = 3.0
    min_trace_length: float = 3.0        # seconds; ~5 MS1 points per peak at 1.2 Hz
    max_trace_length: float = -1.0       # capping this costs minutes per file for nothing
    width_filtering: str = "fixed"
    # Peak width bounds in SECONDS, enforced when width_filtering is "fixed". The OpenMS
    # default ceiling of 60 s is a whole minute — on a 25 min gradient whose peaks have a
    # median FWHM of 4.7 s that is not a peak, it is an undeconvolved trace. Left at 60 the
    # widths have a median close to Compound Discoverer's mean but a tail reaching 14x that,
    # and since every retention window in the peak finder scales with a feature's width, one
    # such feature widens its own filter windows. 20 s is about four times both the observed
    # median and `chrom_fwhm`, and drops 3.2% of features.
    chrom_fwhm: float = 5.0
    min_fwhm: float = 1.0
    max_fwhm: float = 20.0
    isotope_filtering_model: str = "none"
    remove_single_traces: bool = True
    charge_low: int = 1
    charge_high: int = 2
    # Linking across files
    link_rt_tolerance: float = 15.0      # seconds
    link_mz_tolerance: float = 10.0      # ppm

    # Alignment. The pose-clustering defaults are built for maps that disagree by minutes and
    # by 0.3 Da; on high-resolution data with a couple of seconds of drift they pair features
    # that have nothing to do with each other and make the spread worse. These are tight
    # enough to only pair genuine counterparts.
    align_enabled: bool = True
    align_rt_max_difference: float = 30.0    # seconds
    align_mz_max_difference: float = 10.0    # ppm
    align_max_shift: float = 60.0            # seconds

    # Targeted recovery of identified peaks that feature assembly consumed as isotopes.
    recover_identified: bool = True
    recovery_ppm: float = 20.0               # same tolerance the peak finder associates on
    recovery_rt_window: float = 15.0         # seconds either side of the MS2

    def to_dict(self) -> dict:
        return {f: getattr(self, f) for f in self.__dataclass_fields__}


def _set(params, values: dict):
    for key, value in values.items():
        key = key.encode() if isinstance(key, str) else key
        if isinstance(value, str):
            value = value.encode()
        params.setValue(key, value)


SNR_NEAR_FWHM = 2.0      # baseline starts this many FWHM before the apex
SNR_FAR_FWHM = 5.0       # and ends here
SNR_MIN_SCANS = 8        # fewer than this and no estimate is made rather than a bad one


def signal_to_noise(scans, mz: float, rt_seconds: float, fwhm: float,
                    tol: float = 0.01, floor: float = 1.0) -> float | None:
    """Peak height over a robust local baseline sigma. None when the baseline cannot be measured.

    ⚠ **Leading side only.** A real chromatographic peak tails, so the region AFTER the apex still
    carries analyte where the region before does not. A symmetric window measures the peak's own
    tail and reports it as noise — which deflates S/N worst for exactly the tailing peaks the
    number is meant to flag.

    ⚠ **MAD, not standard deviation.** A baseline window this wide frequently contains another
    peak. One such peak destroys an SD and leaves a MAD essentially untouched.

    ⚠ **The convention has to travel with the number.** This is height over a robust sigma
    (1.4826 x MAD). The pharmacopoeial definition uses peak-to-peak noise and gives values roughly
    five times lower for the same peak, so an unqualified "S/N > 3" means nothing.

    Returns None near the start of the gradient, where there is no baseline in front of the peak.
    That is honest: an early-eluting peak has no leading baseline, and inventing one from its own
    tail would be worse than declining.
    """
    if fwhm <= 0 or not scans:
        return None
    lo = rt_seconds - SNR_FAR_FWHM * fwhm
    hi = rt_seconds - SNR_NEAR_FWHM * fwhm
    base, height = [], 0.0
    half = 0.5 * fwhm
    for rt, mzs, its in scans:
        t = rt * 60.0
        if not (lo < t < hi) and abs(t - rt_seconds) > half:
            continue
        if len(mzs) == 0:
            value = 0.0
        else:
            j = _np.searchsorted(mzs, mz)
            value = 0.0
            for k in (j - 1, j, j + 1):
                if 0 <= k < len(mzs) and abs(mzs[k] - mz) <= tol:
                    value = max(value, float(its[k]))
        if lo < t < hi:
            base.append(value)
        else:
            height = max(height, value)
    if len(base) < SNR_MIN_SCANS or height <= 0:
        return None
    b = _np.asarray(base)
    sigma = 1.4826 * float(_np.median(_np.abs(b - _np.median(b))))
    # ⚠ A clean extracted chromatogram has an EMPTY baseline — in centroided data an m/z with no
    # signal contributes no peak at all, so every baseline scan reads exactly zero and the MAD is
    # zero with it. Returning None there would leave S/N blank for precisely the cleanest peaks,
    # which is the opposite of useful.
    #
    # Noise cannot be below the detector's own floor, so that is the substitute. The number then
    # means "at least this", which is the honest reading when the local baseline is unmeasurable.
    return float(height / max(sigma, floor))


def detect_features(mzml: str | Path, params: FeatureParams, log=None):
    """One file to a FeatureMap."""
    from pyopenms import (ElutionPeakDetection, FeatureFindingMetabo, FeatureMap,
                          MassTraceDetection, MSExperiment)

    from .calibrate import load_centroided

    experiment = load_centroided(mzml, log=log)

    ms1 = MSExperiment()
    for spectrum in experiment:
        if spectrum.getMSLevel() == 1:
            ms1.addSpectrum(spectrum)
    ms1.sortSpectra(True)

    mtd = MassTraceDetection()
    p = mtd.getParameters()
    _set(p, {"mass_error_ppm": params.mass_error_ppm,
             "noise_threshold_int": params.noise_threshold,
             "chrom_peak_snr": params.chrom_peak_snr,
             "min_trace_length": params.min_trace_length,
             "max_trace_length": params.max_trace_length})
    mtd.setParameters(p)
    traces = []
    mtd.run(ms1, traces, 0)

    epd = ElutionPeakDetection()
    p = epd.getParameters()
    _set(p, {"width_filtering": params.width_filtering,
             "chrom_peak_snr": params.chrom_peak_snr,
             "chrom_fwhm": params.chrom_fwhm,
             "min_fwhm": params.min_fwhm,
             "max_fwhm": params.max_fwhm})
    epd.setParameters(p)
    split_traces = []
    epd.detectPeaks(traces, split_traces)

    ffm = FeatureFindingMetabo()
    p = ffm.getParameters()
    _set(p, {"isotope_filtering_model": params.isotope_filtering_model,
             "remove_single_traces": "true" if params.remove_single_traces else "false",
             "mz_scoring_by_elements": "false",
             "report_convex_hulls": "true",
             "charge_lower_bound": params.charge_low,
             "charge_upper_bound": params.charge_high})
    ffm.setParameters(p)

    feature_map = FeatureMap()
    chromatograms = []
    ffm.run(split_traces, feature_map, chromatograms)
    feature_map.setUniqueIds()
    # ⚠ Computed HERE, where the spectra are already in memory. An earlier version left a module
    # global for a caller to populate and nothing ever did, so the column existed and was empty in
    # every run — a silent failure that a batch would have carried all night.
    # ⚠ Assign back through the INDEX. Iterating a pyOpenMS FeatureMap yields COPIES, so
    # `for f in feature_map: f.setMetaValue(...)` writes to a temporary and is lost without error —
    # which is how this shipped populating 0 of 1,815 features while the estimator itself worked.
    scans = [(sp.getRT() / 60.0, *sp.get_peaks()) for sp in ms1]
    floor = max(params.noise_threshold * 0.2, 1.0)
    for i in range(feature_map.size()):
        feature = feature_map[i]
        width = feature.getWidth() or (feature.getMetaValue("FWHM") or 0.0)
        if width <= 0:
            continue
        value = signal_to_noise(scans, feature.getMZ(), feature.getRT(), float(width), floor=floor)
        if value is not None:
            feature.setMetaValue("SNR", float(value))
            feature_map[i] = feature

    return feature_map, split_traces


def recover_targeted_features(split_traces, feature_map, targets, params, log=None):
    """Build features for confident identifications that feature finding left with none.

    Feature assembly can consume a chromatographic peak as another feature's isotope, and a
    trace can only belong to one feature, so the peak disappears. The systematic case is a
    lipid class whose members differ by a double bond: 2 H is 2.01565 Da and 2 x 13C is
    2.00671, only 8.9 mDa apart, so at m/z 800 the monoisotopic peak of the saturated species
    lands within 11 ppm of the M+2 isotope of the one with one more double bond. Sphingomyelins
    lose whole species to this.

    Lowering thresholds to get those peaks back costs six to eleven times more features, almost
    all of it noise. This instead asks a targeted question: is there an elution peak here that
    **an MS2 spectrum has already identified with a passing score**, and that no feature
    currently represents? If so, promote it. Nothing is recovered without an identification
    standing behind it, so the pass cannot manufacture features out of noise — the worst it can
    do is promote a peak that really was an isotope, at a mass and time where a spectrum
    matched a library entry anyway.

    `targets` are (m/z, retention in seconds) pairs. Returns the number recovered.
    """
    from pyopenms import Feature

    if not targets:
        return 0

    tolerance_ppm = params.recovery_ppm
    window = params.recovery_rt_window

    # Index targets and existing features by integer m/z so each trace is a couple of lookups.
    target_index: dict[int, list[tuple[float, float]]] = {}
    for mz, rt in targets:
        target_index.setdefault(int(mz), []).append((mz, rt))

    covered: dict[int, list[tuple[float, float]]] = {}
    for feature in feature_map:
        covered.setdefault(int(feature.getMZ()), []).append((feature.getMZ(), feature.getRT()))

    def matches(index, mz, rt):
        for key in (int(mz) - 1, int(mz), int(mz) + 1):
            for other_mz, other_rt in index.get(key, ()):
                if abs(other_mz - mz) / mz * 1e6 < tolerance_ppm and abs(other_rt - rt) < window:
                    return True
        return False

    recovered = 0
    for trace in split_traces:
        mz, rt = trace.getCentroidMZ(), trace.getCentroidRT()
        if not matches(target_index, mz, rt):
            continue
        if matches(covered, mz, rt):
            continue

        feature = Feature()
        feature.setMZ(mz)
        feature.setRT(rt)
        feature.setIntensity(trace.computePeakArea())
        feature.setCharge(1)
        feature.setWidth(trace.getFWHM())
        feature.setMetaValue("recovered", "true")
        feature_map.push_back(feature)
        covered.setdefault(int(mz), []).append((mz, rt))
        recovered += 1

    if recovered:
        feature_map.setUniqueIds()
        if log:
            log(f"  recovered {recovered} features from identified peaks with no feature")
    return recovered


def mzml_is_complete(path: str | Path) -> bool:
    """Does this mzML actually end?

    ⚠ Checked BEFORE handing the file to pyOpenMS, because pyOpenMS is not reliably survivable on a
    damaged one. A file truncated mid-`<binary>` raised `Parse Error`, which can be caught; a file
    damaged differently **segfaulted the interpreter**, which cannot. A 200-byte read is the only
    defence that works for both.

    Cheap enough to run over everything: 2,143 staged files took minutes over SMB, and it finds
    every bad file at once rather than one per fifteen-minute run.
    """
    try:
        with open(path, "rb") as fh:
            fh.seek(max(0, Path(path).stat().st_size - 200))
            tail = fh.read()
    except OSError:
        return False
    return b"</indexedmzML>" in tail or b"</mzML>" in tail


def file_polarity(mzml: str | Path) -> str:
    """Polarity of a run, from the mzML itself.

    pyOpenMS reports a feature's charge as a positive magnitude whatever the polarity, so the
    sign cannot be recovered from the feature. Taking it from the charge silently labels every
    feature positive — which is accidentally right in a positive run and means *no* negative
    identification can ever attach to a feature in a negative one.

    ⚠ Returns "" for a file that cannot be parsed, rather than raising. One truncated download —
    an mzML ending mid-base64 inside a `<binary>` element — used to abort an entire study from here,
    discarding every file already processed. MTBLS2016 lost a 147-file run and ~15 minutes of
    completed feature detection to a single bad file, and it had TWO, so fixing the one named in the
    traceback would have failed again later in the same run.

    The caller decides what to do with "". Skipping is not free either — a study that quietly drops
    3 of 147 injections is a different study — so the skip must be logged, not swallowed.
    """
    from pyopenms import MSExperiment, MzMLFile

    if not mzml_is_complete(mzml):
        return ""
    experiment = MSExperiment()
    try:
        MzMLFile().load(str(mzml), experiment)
    except RuntimeError:
        return ""
    for spectrum in experiment:
        polarity = spectrum.getInstrumentSettings().getPolarity()
        if polarity in (1, 2):
            return "+" if polarity == 1 else "-"
    return "+"


def align(feature_maps: list, params: "FeatureParams | None" = None) -> None:
    """Put every file on one retention-time axis, using the largest map as the anchor.

    `align` only *computes* the transformation — it fills the TransformationDescription and
    leaves the map untouched. Applying it takes a second call through MapAlignmentTransformer.
    Miss that and the pipeline runs cleanly, produces plausible output, and has silently done
    no alignment at all; the only symptom is that measured and aligned retention times are
    identical, which four identical QC injections cannot reveal.
    """
    from pyopenms import MapAlignmentAlgorithmPoseClustering, MapAlignmentTransformer

    if len(feature_maps) < 2:
        return
    aligner = MapAlignmentAlgorithmPoseClustering()
    if params is not None:
        p = aligner.getParameters()
        _set(p, {"superimposer:mz_pair_max_distance": params.align_mz_max_difference / 1e6 * 800,
                 "superimposer:max_shift": params.align_max_shift,
                 "pairfinder:distance_RT:max_difference": params.align_rt_max_difference,
                 "pairfinder:distance_MZ:max_difference": params.align_mz_max_difference,
                 "pairfinder:distance_MZ:unit": "ppm"})
        aligner.setParameters(p)
    reference = max(feature_maps, key=lambda m: m.size())
    aligner.setReference(reference)
    transformer = MapAlignmentTransformer()
    failed = []
    for feature_map in feature_maps:
        if feature_map is reference:
            continue
        transformation = _new_transformation()
        try:
            aligner.align(feature_map, transformation)
        except RuntimeError as exc:
            # ⚠ One unalignable injection must not discard the study.
            #
            # Pose clustering fits a transformation from pairs matched against the reference — the
            # LARGEST map. When an injection is much sparser, there may be too few pairs within
            # tolerance and pyOpenMS raises `no data points for 'linear' model`. On a real study
            # that meant a 945-feature blank killed a 224-file run, and the same error took out
            # three datasets in one batch.
            #
            # A map that cannot be aligned keeps its own retention times, which is exactly what an
            # identity transformation would give it. That is a small, local inaccuracy; losing the
            # study is not. The names are reported so a reader can see WHICH injections were left
            # untransformed rather than discovering it as an unexplained retention shift.
            failed.append((getattr(feature_map, "getMetaValue", lambda _: None)("spectra_data"),
                           str(exc)))
            continue
        transformer.transformRetentionTimes(feature_map, transformation, True)
        try:
            feature_map.updateRanges()
        except Exception:
            pass
    if failed:
        names = []
        for meta, _ in failed:
            try:
                names.append(Path(meta[0].decode()).name)
            except Exception:
                names.append("?")
        print(f"  ⚠ {len(failed)} of {len(feature_maps)} injections could not be aligned to the "
              f"reference and keep their own retention times: {', '.join(names[:6])}"
              + (" ..." if len(names) > 6 else "")
              + ". Usually far sparser than the reference — blanks against a rich sample.")


def _new_transformation():
    from pyopenms import TransformationDescription
    return TransformationDescription()


def link(feature_maps: list, params: FeatureParams):
    """One consensus feature per real compound, spanning the files it was seen in."""
    from pyopenms import ConsensusMap, FeatureGroupingAlgorithmQT

    grouper = FeatureGroupingAlgorithmQT()
    p = grouper.getParameters()
    _set(p, {"distance_RT:max_difference": params.link_rt_tolerance,
             "distance_MZ:max_difference": params.link_mz_tolerance,
             "distance_MZ:unit": "ppm"})
    grouper.setParameters(p)

    consensus = ConsensusMap()
    grouper.group(feature_maps, consensus)
    return consensus


def feature_snr(feature_maps) -> list[dict]:
    """{feature unique id: S/N} per map.

    ⚠ The consensus map holds `FeatureHandle`, not `Feature` — it carries position and intensity
    and NO metadata, and has no `metaValueExists` at all. Reading the S/N off the handle raised
    AttributeError and killed every run in a batch. The unique id is the only thing that survives
    the consensus step, so the value has to be looked up by it, exactly as retention already is.
    """
    return [{f.getUniqueId(): float(f.getMetaValue("SNR"))
             for f in fm if f.metaValueExists("SNR")} for fm in feature_maps]


def original_retentions(feature_maps: list) -> list[dict]:
    """Each feature's retention time *before* alignment, keyed by unique id.

    Must be called before `align`, which rewrites the maps in place. Without this the measured
    and aligned retention times are the same number and the peak finder's retention correction
    has nothing to fit — which is exactly the failure a QC-only test cannot show, because
    identical injections have nothing to align.
    """
    return [{f.getUniqueId(): f.getRT() for f in feature_map} for feature_map in feature_maps]


def build_compound_groups(consensus, feature_maps: list, sample_names: list[str],
                          min_feature_count: int = 2, roles: dict[str, str] | None = None,
                          originals: list[dict] | None = None,
                          snrs: list[dict] | None = None,
                          polarities: list[str] | None = None
                          ) -> tuple[list[CompoundGroup], list[Sample]]:
    """Turn a ConsensusMap into the compound group / compound / feature model.

    pyOpenMS has no notion of a "compound" sitting between the consensus feature and the
    per-file feature, so each file's contribution becomes one compound carrying one feature.
    Adduct grouping, which Compound Discoverer does and which is what makes a CD compound hold
    several ions, is left to the peak finder's adduct sweep instead.
    """
    roles = roles or {}
    samples = [Sample(file=name, cd_id=f"F{i+1}", role=roles.get(name, "sample"))
               for i, name in enumerate(sample_names)]
    blank_columns = {i for i, s in enumerate(samples) if s.role == "blank"}
    groups: list[CompoundGroup] = []

    for consensus_feature in consensus:
        handles = list(consensus_feature.getFeatureList())
        # A blank must not satisfy "detected in at least N files": two blanks sharing a
        # contaminant would otherwise create a compound group out of nothing.
        if sum(1 for h in handles if h.getMapIndex() not in blank_columns) < min_feature_count:
            continue

        areas = [0.0] * len(samples)
        compounds: list[Compound] = []
        max_area = 0.0

        for handle in handles:
            index = handle.getMapIndex()
            if index >= len(samples):
                continue
            sample = samples[index]
            rt = handle.getRT() / 60.0
            mz = handle.getMZ()
            area = handle.getIntensity()
            charge = abs(handle.getCharge() or 1)
            if polarities is not None and index < len(polarities) and polarities[index] == "-":
                charge = -charge
            width = getattr(handle, "getWidth", lambda: 0.0)() / 60.0

            areas[index] = area
            max_area = max(max_area, area)

            feature = Feature(adduct="", charge=charge, mw=mz, mass=mz, retention=rt,
                              fwhm=width or 0.05, mi=1, area=area,
                              parent_area_percent=100.0, sample=sample)
            # Measured retention, before alignment moved it. Falls back to the aligned value
            # only when the original is unavailable.
            measured = None
            if originals is not None and index < len(originals):
                measured = originals[index].get(handle.getUniqueId())
            feature.real_retention = (measured / 60.0) if measured is not None else rt
            # ⚠ Measured against the sample's OWN scans. Noise varies across a gradient and
            # between injections, so a global floor cannot stand in for a local baseline — which
            # is the whole reason for computing this per peak rather than reusing the run's
            # derived threshold.
            if snrs is not None and index < len(snrs):
                feature.snr = snrs[index].get(handle.getUniqueId())
            compound = Compound(mw=mz, retention=rt, fwhm=feature.fwhm, max_mi=1,
                                n_adducts=1, area=area, sample=sample)
            compound.features.append(feature)
            compounds.append(compound)

        group = CompoundGroup(name="", mw=consensus_feature.getMZ(),
                              retention=consensus_feature.getRT() / 60.0,
                              max_area=max_area, has_ms2=False, areas=areas)
        group.compounds = compounds
        groups.append(group)

    groups.sort(key=lambda g: g.retention)
    return groups, samples


def run_feature_detection(mzml_files: list[Path], params: FeatureParams | None = None,
                          min_feature_count: int = 2, log=None,
                          roles: dict[str, str] | None = None,
                          targets: dict[str, list[tuple[float, float]]] | None = None
                          ) -> tuple[list[CompoundGroup], list[Sample]]:
    params = params or FeatureParams()
    say = log or (lambda _: None)

    feature_maps = []
    polarities = []
    targets = targets or {}
    unreadable: list[str] = []
    for path in mzml_files:
        polarity = file_polarity(path)
        if not polarity:
            # ⚠ Named, counted, and carried to the end — never silently dropped. A run that
            # processes 145 of 147 files is a different experiment from one that processes 147.
            unreadable.append(Path(path).name)
            say(f"⚠ SKIPPED {Path(path).name}: cannot be parsed — truncated or corrupt mzML. "
                f"Feature detection continues without it.")
            continue
        polarities.append(polarity)
        feature_map, split_traces = detect_features(path, params)
        found = feature_map.size()
        if params.recover_identified:
            recover_targeted_features(split_traces, feature_map,
                                      targets.get(Path(path).stem, []), params, log=say)
        say(f"{Path(path).name}: {found} features"
            + (f" (+{feature_map.size() - found} recovered)" if feature_map.size() > found else ""))
        feature_maps.append(feature_map)

    originals = original_retentions(feature_maps)
    snrs = feature_snr(feature_maps)
    if params.align_enabled:
        align(feature_maps, params)
    shifts = [abs(f.getRT() - originals[i][f.getUniqueId()])
              for i, fm in enumerate(feature_maps) for f in fm
              if f.getUniqueId() in originals[i]]
    if shifts:
        shifts.sort()
        say(f"alignment moved features by a median {shifts[len(shifts)//2]:.1f} s, "
            f"90th percentile {shifts[int(len(shifts)*0.9)]:.1f} s")

    consensus = link(feature_maps, params)
    say(f"{consensus.size()} consensus features across {len(feature_maps)} files")

    if len(set(polarities)) == 1:
        say(f"polarity: {polarities[0]}")
    else:
        say(f"polarity: mixed {sorted(set(polarities))}")

    return build_compound_groups(consensus, feature_maps,
                                 [Path(p).stem for p in mzml_files],
                                 min_feature_count=min_feature_count, roles=roles,
                                 originals=originals, snrs=snrs, polarities=polarities)
