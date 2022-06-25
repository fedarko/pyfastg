import pytest
from ..pyfastg import extract_node_attrs


def test_basic():
    attrs = extract_node_attrs("EDGE_3_length_100_cov_28.087'")
    assert attrs == {"name": "3-", "length": 100, "cov": 28.087}

    attrs = extract_node_attrs("EDGE_123_length_160_cov_90")
    assert attrs == {"name": "123+", "length": 160, "cov": 90}


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
