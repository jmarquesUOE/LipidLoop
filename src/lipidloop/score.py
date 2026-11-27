"""Spectral similarity, reproducing LipiDex exactly.

See docs/LIPIDEX_ALGORITHM.md section 1. Implemented from
`spectrum_searcher/SampleSpectrum.java :: calcDotProduct`, including the two quirks that a
naive reimplementation would smooth over:

  * unmatched sample peaks are halved before weighting, forward score only;
  * matched peaks weight mass as m**massWeight, unmatched peaks as massWeight*m.

Both are in the original. They are reproduced deliberately so scores are comparable with
existing LipiDex results rather than merely similar.
"""
from __future__ import annotations

INT_WEIGHT = 1.2     # SpectrumSearcher.java line 43
MASS_WEIGHT = 0.9    # SpectrumSearcher.java line 44


def dot_product(sample_mz, sample_int, lib_mz, lib_int, mz_tol: float = 0.01,
                reverse: bool = False, mass_weight: float = MASS_WEIGHT,
                int_weight: float = INT_WEIGHT, presorted: bool = False) -> float:
    """LipiDex dot product, scaled to 0-1000.

    reverse=True ignores peaks present in only one of the two spectra, which is what makes
    the reverse score robust to co-isolated contaminants.

    Plain Python lists rather than numpy: a search calls this a few hundred thousand times on
    spectra of a few dozen peaks each, where array construction costs more than the arithmetic
    it saves. `presorted` skips the sort when the caller already keeps peaks in mass order,
    which the library and sample readers both do.
    """
    if presorted:
        s_mz, s_in = sample_mz, sample_int
        l_mz, l_in = lib_mz, lib_int
    else:
        s = sorted(zip(map(float, sample_mz), map(float, sample_int)))
        l = sorted(zip(map(float, lib_mz), map(float, lib_int)))
        s_mz = [p[0] for p in s]
        s_in = [p[1] for p in s]
        l_mz = [p[0] for p in l]
        l_in = [p[1] for p in l]

    n_lib, n_sample = len(l_mz), len(s_mz)

    numer = lib_sum = sample_sum = 0.0
    i = j = 0
    while i < n_lib and j < n_sample:
        diff = s_mz[j] - l_mz[i]
        if abs(diff) > mz_tol:
            if l_mz[i] < s_mz[j]:
                if not reverse:
                    lib_sum += (mass_weight * l_mz[i] * l_in[i] ** int_weight) ** 2
                i += 1
            else:
                if not reverse:
                    sample_sum += (mass_weight * s_mz[j] * (s_in[j] / 2.0) ** int_weight) ** 2
                j += 1
        else:
            lw = l_mz[i] ** mass_weight * l_in[i] ** int_weight
            sw = s_mz[j] ** mass_weight * s_in[j] ** int_weight
            lib_sum += lw ** 2
            sample_sum += sw ** 2
            numer += sw * lw
            i += 1
            j += 1

    if not reverse:
        while i < n_lib:
            lib_sum += (mass_weight * l_mz[i] * l_in[i] ** int_weight) ** 2
            i += 1
        while j < n_sample:
            sample_sum += (mass_weight * s_mz[j] * (s_in[j] / 2.0) ** int_weight) ** 2
            j += 1

    if numer <= 0.0 or sample_sum <= 0.0 or lib_sum <= 0.0:
        return 0.0
    return 1000.0 * numer ** 2 / (sample_sum * lib_sum)


def ppm_diff(observed: float, reference: float) -> float:
    """SampleSpectrum.java :: calcPPMDiff — signed, relative to the reference mass."""
    return (observed - reference) / reference * 1e6
