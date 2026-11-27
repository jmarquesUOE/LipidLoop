

def test_a_hydroxyl_chain_prefix_survives_into_the_shorthand():
    """`OH-` was read by nothing, so a hydroxylated lipid exported as its unmodified parent.

    LipiDex writes hydroxylated glycerophospholipids as `PC[OH] OH-16:0_18:2`. The `OH-` prefix is
    matched by neither the ether pattern (`^([OP])-` needs O or P immediately before the hyphen)
    nor the sphingoid `d`/`t`/`m` pattern, so it was skipped in silence and the chains parsed as if
    nothing were attached.

    The consequence is worse than a missing annotation: the row scores as a MISS against a
    deposit's `PC …;O` and as a FALSE MATCH against its ordinary `PC`, in the same comparison.
    331 rows across the staged validation results were affected.
    """
    from lipidloop.lsi import shorthand

    assert shorthand("PC[OH] OH-16:0_18:2", "MS2") == "PC 16:0_18:2;O"
    assert shorthand("PE[OH] OH-18:0_20:4", "MS2") == "PE 18:0_20:4;O"

    # The ether prefixes must be untouched — `O-` is a linkage, `OH-` is an oxygen.
    assert shorthand("PC O-16:0_18:2", "MS2") == "PC O-16:0_18:2"
    assert shorthand("PC P-16:0_18:2", "MS2") == "PC P-16:0_18:2"
    # And the sphingoid bases still carry their own oxygen count.
    assert shorthand("SM d34:1", "MS2") == "SM 34:1;O2"
    assert shorthand("PC 16:0_18:2", "MS2") == "PC 16:0_18:2"
