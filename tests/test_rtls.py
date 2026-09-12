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
