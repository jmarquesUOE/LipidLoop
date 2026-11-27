"""The re-injection candidate list: low abundance, never fragmented, library-plausible.

Nothing here is an identification, so the tests that matter are about what keeps a feature OFF
the list as much as what puts it on -- an already-identified, already-fragmented, already-filtered,
or unmatched feature has no business being a re-injection target.
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from lipidloop.inclusion import Included, find_inclusion_candidates, write_inclusion_list  # noqa: E402


class _FakeGroup:
    """Stands in for `CompoundGroup` -- only the fields `find_inclusion_candidates` reads."""

    def __init__(self, quant_ion, retention, max_area, final_lipid_id=None, keep=True,
                avg_fwhm=0.05):
        self.quant_ion = quant_ion
        self.retention = retention
        self.max_area = max_area
        self.final_lipid_id = final_lipid_id
        self.keep = keep
        self.avg_fwhm = avg_fwhm


class _FakeSpectrum:
    def __init__(self, name):
        self.name = name


class _FakeIndex:
    """Stands in for `LibraryIndex` -- a fixed mass -> spectrum-index mapping, tolerance-free."""

    def __init__(self, entries):
        self.spectra = [_FakeSpectrum(name) for _mz, name in entries]
        self._masses = [mz for mz, _name in entries]

    def candidates(self, precursor, tol=0.01):
        return [i for i, mz in enumerate(self._masses) if abs(mz - precursor) <= tol]


def test_low_abundance_unfragmented_library_match_is_a_candidate():
    identified = [_FakeGroup(quant_ion=760.5851, retention=12.0, max_area=1_000_000,
                             final_lipid_id="PC 34:1")]
    candidate = _FakeGroup(quant_ion=700.5000, retention=8.0, max_area=1_000)
    index = _FakeIndex([(700.5000, "PC 32:0; [M+H]+")])
    found = find_inclusion_candidates(identified + [candidate], all_scans=[], index=index)
    assert [round(c.mz, 4) for c in found] == [700.5]
    assert found[0].library_names == ["PC 32:0; [M+H]+"]


def test_already_identified_group_is_never_a_candidate():
    group = _FakeGroup(quant_ion=700.5, retention=8.0, max_area=1_000, final_lipid_id="PC 32:0")
    index = _FakeIndex([(700.5, "PC 32:0; [M+H]+")])
    assert find_inclusion_candidates([group], all_scans=[], index=index) == []


def test_high_abundance_feature_is_not_a_candidate():
    """Above the median of this batch's own MS2-confirmed lipids -- it isn't losing Top-N."""
    identified = [_FakeGroup(quant_ion=760.5851, retention=12.0, max_area=1_000,
                             final_lipid_id="PC 34:1")]
    bright = _FakeGroup(quant_ion=700.5, retention=8.0, max_area=1_000_000)
    index = _FakeIndex([(700.5, "PC 32:0; [M+H]+")])
    assert find_inclusion_candidates(identified + [bright], all_scans=[], index=index) == []


def test_feature_that_was_fragmented_is_not_a_candidate():
    """The whole point: only features that never got an MS2 attempt belong on this list."""
    identified = [_FakeGroup(quant_ion=760.5851, retention=12.0, max_area=1_000_000,
                             final_lipid_id="PC 34:1")]
    candidate = _FakeGroup(quant_ion=700.5, retention=8.0, max_area=1_000)
    scans = [(("f", 1), 700.5, 8.0)]     # an MS2 scan did land on this mass and retention
    index = _FakeIndex([(700.5, "PC 32:0; [M+H]+")])
    found = find_inclusion_candidates(identified + [candidate], all_scans=scans, index=index)
    assert found == []


def test_no_library_match_is_not_a_candidate():
    identified = [_FakeGroup(quant_ion=760.5851, retention=12.0, max_area=1_000_000,
                             final_lipid_id="PC 34:1")]
    candidate = _FakeGroup(quant_ion=700.5, retention=8.0, max_area=1_000)
    index = _FakeIndex([(999.9, "PC 40:0; [M+H]+")])
    assert find_inclusion_candidates(identified + [candidate], all_scans=[], index=index) == []


def test_filtered_out_group_is_never_a_candidate():
    """`keep=False` means an earlier filter already judged it not real; not worth re-injecting."""
    identified = [_FakeGroup(quant_ion=760.5851, retention=12.0, max_area=1_000_000,
                             final_lipid_id="PC 34:1")]
    candidate = _FakeGroup(quant_ion=700.5, retention=8.0, max_area=1_000, keep=False)
    index = _FakeIndex([(700.5, "PC 32:0; [M+H]+")])
    assert find_inclusion_candidates(identified + [candidate], all_scans=[], index=index) == []


def test_inclusion_csv_writer(tmp_path):
    """Just the mass list -- m/z and retention, nothing else, no Thermo import form."""
    included = [Included(mz=700.5000, retention=8.0, area=1_000.0,
                         library_names=["PC 32:0; [M+H]+"], n_matches=1)]
    out = tmp_path / "Inclusion_List.csv"
    write_inclusion_list(included, out)
    body = out.read_text()
    assert body == "m/z,Retention (min)\n700.5000,8.00\n"
