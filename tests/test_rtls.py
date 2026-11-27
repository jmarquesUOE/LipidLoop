"""Tests for retention modelling within a lipid class.

The model's job is to say where a member of a class *should* elute. Most of what can go wrong
is it being trusted when it should not be — on too few points, on a class that is not a
homologous series, or with coefficients pointing the wrong way chemically.
"""
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
from lipidloop.peaks import CompoundGroup                                  # noqa: E402
from lipidloop.rtls import (extend_identifications, fit_model, fit_models,  # noqa: E402
                             parse_sum_name, retention_error)


def series(lipid_class="PC", prefix="", intercept=1.0, per_c=0.25, per_db=-0.35,
           carbons=range(30, 44), doubles=(0, 1, 2)):
    """A clean homologous series obeying the model exactly."""
    return [(f"{lipid_class} {prefix}{c}:{d}", intercept + per_c * c + per_db * d)
            for c in carbons for d in doubles]


def test_parses_sum_names_including_prefixes():
    assert parse_sum_name("PC 34:1") == ("PC", "", 34, 1)
    assert parse_sum_name("SM d36:1") == ("SM", "d", 36, 1)
    assert parse_sum_name("Plasmenyl-PE P-38:4") == ("Plasmenyl-PE", "P-", 38, 4)
    assert parse_sum_name("PC[OH] OH-32:5") == ("PC[OH]", "OH-", 32, 5)
    # A molecular name is collapsed first and judged as its sum. Retention depends on total
    # carbons and double bonds, not on how they are split between chains, and the model has no
    # term that could tell two splits apart. Before this, every row reported at molecular
    # resolution was silently exempt from the retention filter.
    assert parse_sum_name("PC 16:0_18:1") == ("PC", "", 34, 1)
    assert parse_sum_name("TG 16:0_18:1_18:2") == ("TG", "", 52, 3)
    assert parse_sum_name("nonsense") is None
    assert parse_sum_name("") is None


def test_recovers_the_coefficients_it_was_built_from():
    models = fit_models(series())
    model = models[("PC", "")]
    assert model.usable
    assert model.per_carbon == pytest.approx(0.25, abs=0.001)
    assert model.per_double_bond == pytest.approx(-0.35, abs=0.001)
    assert model.r2 == pytest.approx(1.0, abs=1e-6)


def test_too_few_points_gives_no_model():
    assert fit_model("PC", "", [(30, 0, 8.0), (32, 0, 8.5)]) is None


def test_an_outlier_is_trimmed_rather_than_absorbed():
    """Without iterative trimming a wild point drags the surface and hides itself."""
    observations = series()
    observations.append(("PC 40:2", 25.0))          # far off the surface
    models = fit_models(observations)
    model = models[("PC", "")]
    assert model.n_used < model.n_total             # it was dropped
    assert model.per_carbon == pytest.approx(0.25, abs=0.01)   # and did not distort the fit
    error = retention_error(models, "PC 40:2", 25.0)
    assert error is not None and abs(error) > 5     # and is now visible as an outlier


def test_a_class_that_is_not_a_homologous_series_is_not_usable():
    """Wrong-sign coefficients mean the identifications are not what they claim to be."""
    scrambled = [("PS 36:1", 20.0), ("PS 38:2", 6.0), ("PS 40:6", 15.0), ("PS 34:0", 9.0),
                 ("PS 42:4", 7.0), ("PS 36:4", 22.0), ("PS 38:0", 11.0), ("PS 40:1", 8.0)]
    model = fit_models(scrambled)[("PS", "")]
    assert not model.usable
    # and therefore cannot be used to remove anything
    assert retention_error({("PS", ""): model}, "PS 36:1", 20.0) is None


def test_unusable_models_never_judge():
    models = fit_models(series(carbons=range(30, 36), doubles=(0,)))
    for model in models.values():
        if not model.usable:
            assert retention_error(models, "PC 34:0", 99.0) is None


def group(mz, rt, polarity="+"):
    g = CompoundGroup(name="", mw=mz, retention=rt, max_area=1e6, has_ms2=False, areas=[])
    g.quant_ion = mz
    g.quant_polarity = polarity
    return g


def test_extension_names_a_feature_where_the_model_predicts_it():
    models = fit_models(series())
    model = models[("PC", "")]
    predicted = model.predict(36, 1)
    unknown = group(790.5, predicted)
    n = extend_identifications([unknown], models, [("PC 36:1", 790.5, "+")])
    assert n == 1 and unknown.rtls_identification == "PC 36:1"
    assert unknown.identification()[0] == "PC 36:1"
    assert unknown.identification_source() == "RT model"


def test_extension_declines_when_the_retention_is_wrong():
    models = fit_models(series())
    unknown = group(790.5, 2.0)          # nowhere near where PC 36:1 should elute
    assert extend_identifications([unknown], models, [("PC 36:1", 790.5, "+")]) == 0
    assert unknown.rtls_identification == ""


def test_extension_declines_ambiguity_rather_than_guessing():
    models = fit_models(series())
    predicted = models[("PC", "")].predict(36, 1)
    unknown = group(790.5, predicted)
    two = [("PC 36:1", 790.5, "+"), ("PC 36:1", 790.5005, "+")]
    assert extend_identifications([unknown], models, two) == 0


def test_extension_never_overwrites_an_ms2_identification():
    models = fit_models(series())
    identified = group(790.5, models[("PC", "")].predict(36, 1))
    identified.final_lipid_id = object()      # already identified from a spectrum
    assert extend_identifications([identified], models, [("PC 36:1", 790.5, "+")]) == 0


def test_one_class_under_two_spellings_fits_one_surface():
    """The retention model must group by CANONICAL class, or a spelling difference halves it.

    The model fits one surface per class and needs MIN_POINTS = 6 members. Our own two libraries
    write the same molecules as `GlcCer[NS]` (LipidBlast, 5,418 entries) and `HexCer[NS]`
    (LipiDex, 4,800), and a study searches both. Grouping on the raw name split one class into two
    half-sized ones, either of which could fall under the minimum — losing the model for
    hexosylceramides entirely, and with it the ability to extend those identifications or to catch
    a subclass substitution in them.
    """
    from lipidloop.rtls import fit_models, parse_sum_name

    assert parse_sum_name("GlcCer[NS] 42:1")[0] == "HexCer[NS]"
    assert parse_sum_name("HexCer[NS] 42:1")[0] == "HexCer[NS]"
    assert parse_sum_name("CAR 16:0")[0] == "AC"
    assert parse_sum_name("TAG 52:2")[0] == "TG"

    # Subclass tags survive: Cer[ADS] and Cer[NP] are different molecules from different enzymes.
    assert parse_sum_name("Cer[ADS] 42:0")[0] == "Cer[ADS]"
    assert parse_sum_name("Cer[NP] 42:0")[0] == "Cer[NP]"

    points = []
    for i, (c, d) in enumerate([(34, 0), (36, 1), (38, 2), (40, 1), (42, 0),
                                (34, 1), (36, 2), (38, 0), (40, 3), (42, 2)]):
        name = ("GlcCer[NS] %d:%d" if i % 2 else "HexCer[NS] %d:%d") % (c, d)
        points.append((name, 2.0 + 0.20 * c - 0.35 * d))

    models = fit_models(points)
    assert len(models) == 1, f"two spellings must fit one surface, got {list(models)}"
    only = next(iter(models.values()))
    assert only.n_used == 10, "every member must reach the fit, not half of them"


# ── flag rows judged outside the model's own training window (v1.1 item 5) ──────────────────────
# TG 30:0 and 36:0 spiked into NIST SRM 1950 (true) eluted before the plasma-TG model's training points
# and were rejected at z -11.3 / -9.7. Withholding the verdict there also released five wrong off-peak
# rows, so the verdict stands and the row is flagged as judged by extrapolation.

def test_a_model_records_the_window_of_its_own_training_points():
    model = fit_models(series())[("PC", "")]
    rts = [rt for _, rt in series()]
    assert model.rt_min == pytest.approx(min(rts)) and model.rt_max == pytest.approx(max(rts))
    assert model.covers(model.rt_min) and model.covers(model.rt_max)
    assert not model.covers(model.rt_min - 0.01) and not model.covers(model.rt_max + 0.01)


def test_a_trimmed_point_does_not_stretch_the_window():
    points = series() + [("PC 30:0", 2.0)]          # far off its own surface: trimmed
    model = fit_models(points)[("PC", "")]
    assert model.rt_min > 2.0


def _finder(groups, n_samples=4):
    from lipidloop.peakfinder import PeakFinder
    from lipidloop.peaks import Sample
    f = PeakFinder.__new__(PeakFinder)
    f.compound_groups = groups
    f.samples = [Sample(file=f"S{i}", role="sample") for i in range(n_samples)]
    f.retention_model_sigma = 3.0
    f.retention_model_max_error = 1.0
    f.retention_model_weak_z = 5.0
    f.retention_model_weak_spectra = 3
    f.retention_model_weak_min_error = 0.075
    return f


def _named(name, retention, n_compounds=4, spectra=1):
    from lipidloop.peaks import Lipid, LipidCandidate, Sample
    lipids = [Lipid(retention=retention, precursor=800.0, sample=Sample(file=f"S{i % 4}", role="sample"),
                    dot=950, rev_dot=990, lipid_string=f"{name} [M+NH4]+;", lib_precursor=800.0,
                    purity=0, is_lipidex=True, purity_array=[], fragment_masses=[])
              for i in range(spectra)]
    g = CompoundGroup(name=name, mw=800.0, retention=retention, max_area=1.0, has_ms2=True, areas=[1.0])
    g.compounds = [object()] * n_compounds
    g.lipid_candidates = [LipidCandidate(lipids[0])]
    for lipid in lipids[1:]:
        g.lipid_candidates[0].identifications.append(lipid)
    g.final_lipid_id = g.lipid_candidates[0]
    g.sum_id = name
    return g


def test_a_verdict_outside_the_window_stands_and_is_flagged():
    train = [_named(f"TG {c}:{d}", 12.0 + 0.25 * (c - 44) - 0.45 * d)
             for c in range(44, 64, 2) for d in (0, 1, 2, 3) if 12.0 + 0.25 * (c - 44) - 0.45 * d >= 11.9]
    spike = _named("TG 30:0", 5.5)                  # outside the window, 3 min before the extrapolation (8.5)
    off_peak = _named("TG 42:0", 13.20)             # inside the window, 1.7 min off its surface
    good = _named("TG 50:1", 12.0 + 0.25 * 6 - 0.45)
    f = _finder(train + [spike, off_peak, good])
    f._filter_by_retention_model()
    model = f.retention_models[("TG", "")]
    assert not model.covers(5.5) and model.covers(13.2)
    assert not spike.keep and spike.filter_reason == "RT model outlier"
    assert spike.retention_outside_window.endswith("(extrapolated)")
    assert not off_peak.keep and off_peak.retention_outside_window == ""
    assert good.keep and good.retention_outside_window == ""


def _tight_pc(noise=0.06):
    """A PC series fitting to ~0.05 min, so 1 min (the error floor) is ~20 sigma."""
    jitter = [noise, -noise, 0.0]
    return [_named(f"PC {c}:{d}", 2.0 + 0.25 * c - 0.35 * d + jitter[(c + d) % 3], spectra=8)
            for c in range(30, 44) for d in (0, 1, 2)]


def test_a_strong_outlier_on_few_spectra_goes_inside_the_error_floor():
    """v1.1 item 8: QE 2x2 NIST PC 20:4 at z -20 on one spectrum was kept by the 1 min floor."""
    train = _tight_pc()
    pred = 2.0 + 0.25 * 36 - 0.35 * 1
    weak = _named("PC 36:1", pred - 0.8, spectra=1)        # 0.8 min: inside the floor, far in sigma
    strong = _named("PC 36:1", pred - 0.8, spectra=4)      # same place, more spectra: kept
    f = _finder(train + [weak, strong])
    f._filter_by_retention_model()
    assert abs(weak.retention_z) > 5 and abs(weak.retention_error) < 1.0
    assert not weak.keep and weak.filter_reason == "RT model outlier (weak MS2)"
    assert strong.keep
    assert all(g.keep for g in train)


def test_the_weak_ms2_rule_can_be_switched_off():
    train = _tight_pc()
    weak = _named("PC 36:1", 2.0 + 0.25 * 36 - 0.35 - 0.8, spectra=1)
    f = _finder(train + [weak])
    f.retention_model_weak_z = 0
    f._filter_by_retention_model()
    assert weak.keep


def test_a_high_z_within_a_peak_width_is_not_an_outlier():
    """A class fitted to ~0.002 min turns 0.02 min into z -10 (PI 16:1_16:0, QE 2x2 NIST arm B)."""
    train = _tight_pc(noise=0.004)
    pred = 2.0 + 0.25 * 36 - 0.35 * 1
    near = _named("PC 36:1", pred - 0.05, spectra=1)       # many sigma, but under 0.075 min
    far = _named("PC 36:1", pred - 0.10, spectra=1)
    f = _finder(train + [near, far])
    f._filter_by_retention_model()
    assert abs(near.retention_z) > 5 and near.keep
    assert not far.keep and far.filter_reason == "RT model outlier (weak MS2)"
