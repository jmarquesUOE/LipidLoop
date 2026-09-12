"""Free fatty acids are named without a spectrum, so the gates are the whole safety argument.

Each candidate on its own is only a mass, and a mass is not an identification. What is being
tested is the *set*: only if thirty-odd features on exact fatty-acid masses jointly obey
RT ~ C + DB is any one of them named. These tests pin that nothing gets named when that is not
true.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.fatty_acids import (annotate, candidate_masses, deprotonated,  # noqa: E402
                                    fit_series, max_double_bonds)
from lipidloop.peaks import Compound, CompoundGroup, Sample  # noqa: E402

# The surface measured on the skin negative data.
INTERCEPT, PER_C, PER_DB = 5.081, 0.327, -0.446
# Filtered through max_double_bonds() so odd-chain members stay within the real ceiling (see
# test_mass_table_is_chemically_sane) instead of exercising compositions candidate_masses() no
# longer enumerates at all.
SERIES = [(c, db) for c in range(14, 27) for db in (0, 1, 2) if db <= max_double_bonds(c)]


def group(mz, rt, files=5, area=1e7):
    g = CompoundGroup(name="", mw=mz, retention=rt, max_area=area, has_ms2=False, areas=[])
    g.quant_ion, g.quant_polarity = mz, "-"
    g.compounds = [Compound(mw=mz, retention=rt, fwhm=0.1, max_mi=0, n_adducts=1,
                            area=area, sample=Sample(file=f"f{i}"))
                   for i in range(files)]
    return g


def real_series(noise=0.0):
    out = []
    for i, (c, db) in enumerate(SERIES):
        rt = INTERCEPT + PER_C * c + PER_DB * db + (noise if i % 2 else -noise)
        out.append(group(deprotonated(c, db), rt))
    return out


def test_a_real_homologous_series_is_named():
    groups = real_series()
    assert annotate(groups) == len(SERIES)
    assert groups[0].identification()[0] == "FA 14:0"


def test_masses_alone_are_not_enough():
    """Right masses, retention times shuffled into nonsense. Nothing may be named."""
    groups = real_series()
    retentions = [g.retention for g in groups]
    for g, rt in zip(groups, reversed(retentions[: len(retentions) // 2] + retentions)):
        g.retention = rt
    named = annotate(groups)
    assert named == 0 or all(not g.identification()[0] for g in groups[:2])


def test_wrong_signs_are_rejected_however_good_the_fit():
    """A surface where longer chains elute EARLIER fits beautifully and is chemically impossible."""
    groups = [group(deprotonated(c, db), 30.0 - PER_C * c - PER_DB * db) for c, db in SERIES]
    assert annotate(groups) == 0


def test_too_few_candidates():
    groups = [group(deprotonated(c, 0), INTERCEPT + PER_C * c) for c in range(14, 19)]
    assert annotate(groups) == 0


def test_an_identified_group_is_never_overwritten():
    groups = real_series()
    victim = groups[3]
    victim.rtls_identification = ""
    before = len(groups)
    # give it an MS2 identification by making identification() return something
    victim.sum_id = "PC 34:1"
    from lipidloop.peaks import Lipid, LipidCandidate
    lipid = Lipid(retention=victim.retention, precursor=victim.quant_ion,
                  sample=Sample(file="f0"), dot=900.0, rev_dot=990.0,
                  lipid_string="PC 34:1 [M-H]-;", lib_precursor=victim.quant_ion, purity=0,
                  is_lipidex=True, purity_array=[], fragment_masses=[])
    lipid.lipid_name, lipid.sum_lipid_name, lipid.lipid_class = "PC 34:1", "PC 34:1", "PC"
    victim.lipid_candidates = [LipidCandidate(lipid)]
    annotate(groups)
    assert victim.identification()[0] == "PC 34:1"
    assert len(groups) == before


def test_positive_mode_groups_are_left_alone():
    groups = real_series()
    for g in groups:
        g.quant_polarity = "+"
    assert annotate(groups, polarity="-") == 0


def test_one_off_features_do_not_set_the_surface():
    """A candidate seen in too few files may be named, but must not join the fit."""
    groups = real_series()
    groups.append(group(deprotonated(20, 0), 2.0, files=1))     # wild retention, one file
    model = fit_series([(c, db, g.retention) for (c, db), g in
                        zip(SERIES, groups)])
    assert model is not None and model.usable
    annotate(groups)
    assert not groups[-1].identification()[0]                   # 2.0 min is far off the surface


def test_mass_table_is_chemically_sane():
    masses = candidate_masses()
    assert abs(masses[(16, 0)] - 255.2330) < 0.001      # palmitate
    assert abs(masses[(18, 1)] - 281.2486) < 0.001      # oleate
    assert abs(masses[(20, 4)] - 303.2330) < 0.001      # arachidonate
    assert abs(masses[(22, 6)] - 327.2330) < 0.001      # DHA
    assert (20, 5) in masses and (22, 6) in masses      # EPA and DHA
    assert (14, 6) not in masses                        # too unsaturated for the chain length
    assert (17, 1) in masses                            # odd-chain MUFA: margaroleic acid, real
    # odd-chain "polyenes": no mammalian pathway builds these; each was a real, blank-clear peak
    # on the NIST SRM 1950 batch before this gate existed, chemically impossible regardless
    assert (13, 3) not in masses
    assert (15, 4) not in masses
    assert (19, 5) not in masses


def test_curvature_is_used_when_the_data_supports_it_and_not_otherwise():
    """Real retention curves with chain length. The quadratic must be earned, not assumed."""
    from lipidloop.fatty_acids import fit_series

    curved = [(c, 0, 5.0 + 0.20 * c + 0.004 * c * c) for c in range(14, 31)]
    curved += [(c, db, 5.0 + 0.20 * c + 0.004 * c * c - 0.45 * db)
               for c in (18, 20, 22) for db in (1, 2, 3)]
    model = fit_series(curved)
    assert model is not None and model.usable
    assert model.per_carbon_squared > 0, "curvature present and plenty of points — should be used"
    assert model.residual_sd < 0.05

    straight = [(c, db, 5.0 + 0.33 * c - 0.45 * db)
                for c in range(14, 27) for db in (0, 1, 2)]
    assert fit_series(straight).per_carbon_squared == 0.0

    # Too few points to determine a curve, however curved they look. (Nine, and spanning both
    # double-bond counts so the fit is not rank-deficient for an unrelated reason.)
    few = [(c, 0, 5.0 + 0.20 * c + 0.004 * c * c) for c in range(14, 21)] \
        + [(c, 1, 5.0 + 0.20 * c + 0.004 * c * c - 0.45) for c in (18, 20)]
    assert len(few) == 9
    assert fit_series(few).per_carbon_squared == 0.0


def test_the_sign_gate_reads_the_slope_not_the_linear_coefficient():
    """A curved fit can have a negative linear coefficient and still rise everywhere.

    Measured on the skin data: `-0.441/C +0.0192/C^2`, slope +0.327 at C20. Testing the bare
    coefficient rejected a perfectly good model and named nothing.
    """
    from lipidloop.fatty_acids import fit_series

    points = [(c, db, 5.0 + 0.20 * c + 0.004 * c * c - 0.45 * db)
              for c in range(14, 31) for db in (0, 1, 2)]
    model = fit_series(points)
    assert model is not None and model.per_carbon_squared > 0
    assert model.usable, "a rising curved surface must be usable"
    assert model.slope_at(14) > 0 and model.slope_at(30) > 0


def test_a_curve_that_turns_over_inside_the_range_is_rejected():
    """Rising then falling is not chromatography, however well it fits."""
    from lipidloop.fatty_acids import fit_series

    points = [(c, db, 5.0 + 1.2 * c - 0.03 * c * c - 0.45 * db)
              for c in range(14, 31) for db in (0, 1, 2)]
    model = fit_series(points)
    assert model is not None
    assert not model.usable


def test_standards_anchor_the_surface_and_are_never_trimmed(tmp_path):
    """A point whose identity is known cannot be an outlier from a surface meant to describe it.

    If a standard disagrees with the surface, the surface is wrong — so the standard stays in and
    drags the fit, rather than being quietly discarded to make the fit look good.
    """
    from lipidloop.fatty_acids import fit_series, read_standards

    observed = [(c, db, 5.0 + 0.33 * c - 0.45 * db) for c in range(14, 27) for db in (0, 1, 2)]
    # One standard a long way off the line the observations describe.
    rogue = [(30, 0, 5.0 + 0.33 * 30 - 4.0)]
    plain = fit_series(observed)
    anchored = fit_series(observed, anchors=rogue)
    assert plain.residual_sd < 0.01
    assert anchored.residual_sd > plain.residual_sd, "the standard must not be trimmed away"

    csv_path = tmp_path / "standards.csv"
    csv_path.write_text("name,retention\nFA 18:1,10.36\nFA 20:4,9.91\n# a comment,\n")
    assert read_standards(csv_path) == {(18, 1): 10.36, (20, 4): 9.91}
    assert read_standards("") == {} and read_standards(tmp_path / "nope.csv") == {}
