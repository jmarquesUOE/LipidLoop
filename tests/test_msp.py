

def test_hexosylceramides_are_not_reported_as_glucosyl():
    """`GlcCer` claims a stereochemistry reversed-phase MS2 cannot establish.

    Glucosyl- and galactosylceramide are stereoisomers: RP chromatography does not resolve them
    and their MS2 spectra are indistinguishable. The evidence supports "a hexose".

    It is also a functional bug. LipidBlast writes `GlcCer[NS]` (5,418 entries) and LipiDex writes
    `HexCer[NS]` (4,800) for the same molecules, and a study searches both. The retention model
    fits one surface per class at MIN_POINTS = 6, so one class under two names splits its anchors
    and can drop both halves below the minimum — losing the RT model for the class entirely.
    """
    from lipidloop.msp import LibrarySpectrum

    def cls(name):
        return LibrarySpectrum(name=name, precursor_mz=1.0, mz=[], intensity=[]).lipid_class

    for written in ("GlcCer[NS] 42:1;O2 [M+H]+;", "GalCer[NS] 42:1;O2 [M+H]+;",
                    "HexCer[NS] 42:1;O2 [M+H]+;"):
        assert cls(written) == "HexCer[NS]", written

    # Sulfatide is a different lipid and must survive untouched — the pattern is anchored, so
    # `SHexCer` does not match, and neither does the disaccharide.
    assert cls("SHexCer 42:1;O2 [M-H]-;") == "SHexCer"
    assert cls("Hex2Cer 34:1;O2 [M+H]+;") == "Hex2Cer"
    assert cls("Cer[ADS] 42:0;O3 [M+H]+;") == "Cer[ADS]"
