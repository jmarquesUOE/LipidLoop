"""The mobile-phase modifier read off the background salt clusters."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.exclusion import SERIES                      # noqa: E402
from lipidloop.modifier import detect                       # noqa: E402

NA_FORMATE = SERIES["sodium formate"]
NA_ACETATE = SERIES["sodium acetate"]


def chain(start, repeat, n):
    """`n` members of a cluster series, which is what the instrument actually sees."""
    return [start + i * repeat for i in range(n)]


def plenty(start, repeat, chains=12, n=6):
    """Enough series members to clear MIN_EVIDENCE.

    Real negative-mode data carries 177-4479 members of a present series, so a test that proves
    a call must supply comparable evidence — otherwise it is really testing the thinness floor.
    Several short chains rather than one long one, which is how clusters actually appear:
    different adducts of the same salt, each with its own ladder.
    """
    out = []
    for i in range(chains):
        out.extend(chain(start + i * 0.37, repeat, n))
    return out


def test_formate_series_is_called_formate():
    assert detect([plenty(158.964, NA_FORMATE)]).call == "formate"


def test_acetate_series_is_called_acetate():
    assert detect([plenty(186.995, NA_ACETATE)]).call == "acetate"


def test_no_series_is_no_evidence_not_a_contradiction():
    """⚠ The failure that matters. A known-formate study returned zero members of either series
    on thin sampling, because its base peaks were all real lipids. That must read as "no
    evidence", never as a reason to override the declared modifier."""
    d = detect([[200.1, 311.4, 590.2, 782.6]])
    assert d.call is None
    assert not d.confident
    assert not d.disagrees_with("formate")
    assert not d.disagrees_with("acetate")


def test_a_mixed_run_refuses_to_choose():
    """MSV000094718 is a comparison OF ammonium salts. Both series present at comparable levels is
    a true description of such a run, and forcing a winner would invent a fact."""
    d = detect([plenty(158.964, NA_FORMATE) + plenty(186.995, NA_ACETATE)])
    assert d.call is None
    assert "mix" in d.reason


def test_a_clear_margin_still_wins_when_both_are_present():
    peaks = plenty(158.964, NA_FORMATE, chains=24) + plenty(186.995, NA_ACETATE, chains=5)
    assert detect([peaks]).call == "formate"


def test_disagreement_is_only_reported_against_a_diagnostic_modifier():
    d = detect([plenty(186.995, NA_ACETATE)])
    assert d.disagrees_with("formate")
    assert not d.disagrees_with("acetate")
    # hydroxy/other libraries are not a formate-vs-acetate question at all
    assert not d.disagrees_with("hydroxy")


def test_two_ions_are_a_coincidence_not_a_series():
    """assign_series' floor, restated here because it is the whole basis of the claim."""
    assert detect([chain(158.964, NA_FORMATE, 2)]).call is None


def test_series_are_decided_per_file_then_pooled():
    """Three ions from three different runs must not form a chain that existed in none of them."""
    per_file = [[158.964], [158.964 + NA_FORMATE], [158.964 + 2 * NA_FORMATE]]
    assert detect(per_file).call is None


def test_the_two_adducts_differ_by_exactly_ch2():
    """The reason this module exists: [M+HCOO]- and [M+CH3COO]- are one CH2 apart, so a wrong
    modifier renames every negative-mode choline lipid to its odd-chain neighbour with zero mass
    error. Guard the arithmetic so nobody 'fixes' the constants apart."""
    formate_adduct, acetate_adduct = 44.99820, 59.01385
    assert abs((acetate_adduct - formate_adduct) - 14.01565) < 0.001


def test_a_thin_series_is_not_a_confident_call():
    """⚠ Found by running this on real data. A ratio cannot see thinness: 6 members against 0 is an
    infinite margin and almost no evidence. Real studies carrying a series give 177-4479 members;
    MSV000095868 gave 6 across three files and was being called with the same confidence."""
    from lipidloop.modifier import MIN_EVIDENCE
    d = detect([chain(158.964, NA_FORMATE, 6)])
    assert d.call is None
    assert "too thin" in d.reason
    # and the same series, present properly, still calls
    assert detect([plenty(158.964, NA_FORMATE)]).call == "formate"
