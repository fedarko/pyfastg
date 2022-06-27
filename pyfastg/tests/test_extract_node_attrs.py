import pytest
from ..pyfastg import extract_node_attrs


def test_basic():
    attrs = extract_node_attrs("EDGE_3_length_100_cov_28.087'")
    assert attrs == {"name": "3-", "length": 100, "cov": 28.087}

    attrs = extract_node_attrs("EDGE_123_length_160_cov_90")
    assert attrs == {"name": "123+", "length": 160, "cov": 90}

    # Check that ".5" is interpretable by python
    attrs = extract_node_attrs("EDGE_123hi_length_160_cov_.5'")
    assert attrs == {"name": "123hi-", "length": 160, "cov": 0.5}


def test_alphabetical_ids():
    attrs = extract_node_attrs("EDGE_asdf_length_100_cov_29.087'")
    assert attrs == {"name": "asdf-", "length": 100, "cov": 29.087}

    attrs = extract_node_attrs("EDGE_lmao_length_1_cov_3")
    assert attrs == {"name": "lmao+", "length": 1, "cov": 3}


def test_bad_declaration():
    bd = "EDGE_asdf_length_100_cog_29.087'"
    with pytest.raises(ValueError) as exc_info:
        extract_node_attrs(bd)
    assert str(exc_info.value) == (
        "Wasn't able to find all expected info (edge name, length, coverage) "
        'in the declaration "{}". Please remember that pyfastg only '
        "supports SPAdes-dialect FASTG files."
    ).format(bd)

    with pytest.raises(ValueError) as exc_info:
        extract_node_attrs("")
    assert str(exc_info.value) == (
        "Wasn't able to find all expected info (edge name, length, coverage) "
        'in the declaration "". Please remember that pyfastg only '
        "supports SPAdes-dialect FASTG files."
    )


def test_bad_coverages_caught_by_regex():
    # The |||| case could've slipped through in earlier pyfastg versions
    bad_decls = (
        "EDGE_asdf_length_100_cov_",
        "EDGE_asdf_length_100_cov_1e10",
        "EDGE_asdf_length_100_cov_||||",
        "EDGE_asdf_length_100_cov_1,234",
        "EDGE_asdf_length_100_cov_1,234.5",
        "EDGE_asdf_length_100_cov_1 234.5",
    )
    for bd in bad_decls:
        with pytest.raises(ValueError) as exc_info:
            extract_node_attrs(bd)
        assert str(exc_info.value) == (
            "Wasn't able to find all expected info (edge name, length, "
            'coverage) in the declaration "{}". Please remember that pyfastg '
            "only supports SPAdes-dialect FASTG files."
        ).format(bd)


def test_bad_coverage_passes_regex_but_bad_float():
    bad_decls = (
        "EDGE_asdf_length_100_cov_.",
        "EDGE_asdf_length_100_cov_...",
        "EDGE_asdf_length_100_cov_1.2.3",
        "EDGE_asdf_length_100_cov_3....",
        "EDGE_asdf_length_100_cov_....3",
    )
    for bd in bad_decls:
        with pytest.raises(ValueError) as exc_info:
            extract_node_attrs(bd)
        assert "could not convert string to float" in str(exc_info.value)


def test_leading_or_trailing_info_not_allowed():
    bad_decls = (
        ">EDGE_asdf_length_100_cov_29.087'",
        "asdfEDGE_asdf_length_100_cov_29.087'",
        "EDGE_asdf_length_100_cov_29.087'ghij",
        "asdfEDGE_asdf_length_100_cov_29.087'ghij",
        "eEDGE_asdf_length_100_cov_29.087",
        "EDGE_asdf_length_100_cov_29.087_",
        "EDGE_asdf_length_100_cov_29.087_asdf",
        "eEDGE_asdf_length_100_cov_29.087'",
    )
    for bd in bad_decls:
        with pytest.raises(ValueError) as exc_info:
            extract_node_attrs(bd)
        assert str(exc_info.value) == (
            "Wasn't able to find all expected info (edge name, length, "
            'coverage) in the declaration "{}". Please remember that pyfastg '
            "only supports SPAdes-dialect FASTG files."
        ).format(bd)
