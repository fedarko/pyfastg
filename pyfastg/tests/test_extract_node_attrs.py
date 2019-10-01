from ..pyfastg import extract_node_attrs


def test_basic():
    attrs = extract_node_attrs("EDGE_3_length_100_cov_28.087'")
    assert attrs == {"name": "3-", "length": 100, "coverage": 28.087}

    attrs = extract_node_attrs("EDGE_123_length_160_cov_90")
    assert attrs == {"name": "123+", "length": 160, "coverage": 90}
