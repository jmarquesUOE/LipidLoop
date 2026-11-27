"""A chain split is reported only where an ion the rivals do not predict is observed."""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.msp import LibrarySpectrum                      # noqa: E402
from lipidloop.mgf import SampleSpectrum                       # noqa: E402
from lipidloop.search import chain_support                     # noqa: E402


def lib(name, mz, inten=None):
    return LibrarySpectrum(name=name, precursor_mz=700.0, mz=list(mz),
                           intensity=list(inten or [100.0] * len(mz)))


def sample(mz):
    return SampleSpectrum(number=1, precursor=700.0, retention=10.0, polarity="-",
                          mz=list(mz), intensity=[100.0] * len(mz))


SHARED = [200.0, 300.0]          # head group / backbone: identical between the two splits


def test_sum_composition_name_asserts_nothing():
    top = lib("Cer[ADS] d42:0", SHARED)
    assert chain_support(top, [], sample(SHARED))[0] is True


def test_no_rival_in_library_is_not_evidence():
    """The common case: 73% of chain names in a delivered table won by default."""
    top = lib("Cer[ADS] d17:0_25:0", SHARED + [400.0])
    ok, seen, diag = chain_support(top, [], sample(SHARED + [400.0]))
    assert ok is False and seen == 0 and diag == 0


def test_diagnostic_ion_observed_keeps_the_split():
    top = lib("Cer[ADS] d17:0_25:0", SHARED + [400.0])
    rival = lib("Cer[ADS] d18:0_24:0", SHARED + [414.0])      # one CH2 away
    ok, seen, diag = chain_support(top, [rival], sample(SHARED + [400.0]))
    assert ok is True and seen == 1 and diag == 1


def test_diagnostic_ion_absent_collapses_to_sum():
    top = lib("Cer[ADS] d17:0_25:0", SHARED + [400.0])
    rival = lib("Cer[ADS] d18:0_24:0", SHARED + [414.0])
    ok, seen, diag = chain_support(top, [rival], sample(SHARED))   # the 400.0 ion is not there
    assert ok is False and seen == 0 and diag == 1


def test_shared_ions_are_never_diagnostic():
    """Both splits predict these, so observing them says nothing about the split."""
    top = lib("Cer[ADS] d17:0_25:0", SHARED)
    rival = lib("Cer[ADS] d18:0_24:0", SHARED)
    ok, beaten, rivals = chain_support(top, [rival], sample(SHARED))
    assert ok is False and beaten == 0 and rivals == 1


def test_weak_library_peaks_do_not_count():
    top = lib("Cer[ADS] d17:0_25:0", SHARED + [400.0], [100.0, 100.0, 1.0])
    rival = lib("Cer[ADS] d18:0_24:0", SHARED + [414.0])
    ok, beaten, rivals = chain_support(top, [rival], sample(SHARED + [400.0]))
    assert ok is False and beaten == 0 and rivals == 1


def test_library_name_carries_its_adduct():
    """⚠ Regression: `sum_composition` leaves an unparseable name alone, so passing the FULL
    library name made every comparison a silent no-op — columns filled in, nothing collapsed."""
    from lipidloop.search import lipid_body
    from lipidloop.peaks import sum_composition

    body, suffix = lipid_body("Cer[ADS] d17:0_25:0 [M+FA-H]-;")
    assert body == "Cer[ADS] d17:0_25:0"
    assert suffix == " [M+FA-H]-;"
    assert sum_composition(body) == "Cer[ADS] d42:0"
    # and the round trip keeps the adduct the rest of the pipeline reads
    assert sum_composition(body) + suffix == "Cer[ADS] d42:0 [M+FA-H]-;"


def test_gate_sees_through_the_adduct():
    top = lib("Cer[ADS] d17:0_25:0 [M+FA-H]-;", SHARED + [400.0])
    rival = lib("Cer[ADS] d18:0_24:0 [M+FA-H]-;", SHARED + [414.0])
    assert chain_support(top, [rival], sample(SHARED))[0] is False       # ion absent -> collapse
    assert chain_support(top, [rival], sample(SHARED + [400.0]))[0] is True


def test_every_rival_must_be_ruled_out():
    """One rival that cannot be excluded is enough to collapse the name."""
    top = lib("Cer[ADS] d17:0_25:0", SHARED + [400.0])
    beaten = lib("Cer[ADS] d18:0_24:0", SHARED + [414.0])      # 400.0 rules this one out
    unbeaten = lib("Cer[ADS] d19:0_23:0", SHARED + [400.0])    # predicts 400.0 too
    ok, n, total = chain_support(top, [beaten], sample(SHARED + [400.0]))
    assert ok is True and n == 1 and total == 1
    ok, n, total = chain_support(top, [beaten, unbeaten], sample(SHARED + [400.0]))
    assert ok is False and n == 1 and total == 2


def test_crowding_alone_does_not_collapse_a_name():
    """⚠ Regression: the union-of-rivals form collapsed 91% of TG purely because the library
    holds ~76 splits per TG sum composition. Many rivals, each individually ruled out, is fine."""
    top = lib("TG 16:1_17:2_20:5", SHARED + [400.0])
    rivals = [lib(f"TG r{i}", SHARED + [414.0 + i]) for i in range(40)]
    ok, n, total = chain_support(top, rivals, sample(SHARED + [400.0]))
    assert ok is True and n == 40 and total == 40


def _group(names_counts, sum_id):
    """A compound group whose candidates are (name, spectra-count) pairs."""
    from lipidloop.peaks import CompoundGroup, LipidCandidate

    g = CompoundGroup(name="", mw=700.0, retention=10.0, max_area=1.0,
                      has_ms2=True, areas=[1.0])
    g.sum_id = sum_id
    cands = []
    for name, count in names_counts:
        c = LipidCandidate.__new__(LipidCandidate)
        c.lipid_name = name
        c.sum_lipid_name = sum_id
        c.lipid_class = "Cer[ADS]"
        c.max_dot = 900.0
        c.max_rev_dot = 1000.0
        c.count = count
        c.purity = 100
        c.identifications = []
        c.preferred_polarity = False
        cands.append(c)
    g.lipid_candidates = cands
    g.purity = 100.0
    return g


def test_feature_keeps_chains_when_most_spectra_carried_them():
    g = _group([("Cer[ADS] d17:0_25:0", 9), ("Cer[ADS] d42:0", 3)], "Cer[ADS] d42:0")
    assert g.identification()[0] == "Cer[ADS] d17:0_25:0"


def test_feature_demotes_when_most_spectra_failed_the_test():
    """⚠ Regression: a per-spectrum gate alone is inert here. The peak finder names a feature
    from its best candidate, so a minority that kept the chain name carried it through and none
    of the six delivered Cer[ADS] names changed."""
    g = _group([("Cer[ADS] d17:0_25:0", 4), ("Cer[ADS] d42:0", 9)], "Cer[ADS] d42:0")
    assert g.identification()[0] == "Cer[ADS] d42:0"


def test_no_collapsed_sibling_means_no_op():
    """With `chain_evidence` off nothing is ever rewritten, so this rule cannot fire."""
    g = _group([("Cer[ADS] d17:0_25:0", 4)], "Cer[ADS] d42:0")
    assert g.identification()[0] == "Cer[ADS] d17:0_25:0"
