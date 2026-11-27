

def test_ganglioside_sialic_acid_suffixes_stop_hiding_the_base_class():
    """Gangliosides carry their subclass as a SUFFIX, so they escaped every other normalisation.

    Our libraries hold 748 `GM3-NANA` and 748 `GM3-NGNA` entries. A deposit writes plain `GM3`.
    Nothing matched, and the class read as one we did not cover at all — the fifth mismatch of this
    kind, after LysoPA/LPA, AC/CAR, GlcCer/HexCer and ST/BA.

    NANA and NGNA are genuinely different molecules — N-acetyl- against N-glycolyl-neuraminic acid,
    16 Da apart, and rodents make NGNA where humans make only NANA — so the distinction is kept as
    a bracket tag rather than discarded. It simply stops concealing the base class.
    """
    from lipidloop.nomenclature import canonical_class

    assert canonical_class("GM3-NANA 34:1") == "GM3[NANA]"
    assert canonical_class("GM3-NGNA 34:1") == "GM3[NGNA]"
    assert canonical_class("GD3-NGNA 36:1") == "GD3[NGNA]"
    assert canonical_class("GM3 34:1") == "GM3"
    # The two sialic acids must NOT merge — they differ by an oxygen.
    assert canonical_class("GM3-NANA 34:1") != canonical_class("GM3-NGNA 34:1")
