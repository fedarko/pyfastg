from ..pyfastg import extract_node_len_cov_rc

def test_basic():
    attrs = extract_node_len_cov_rc("EDGE_3_length_100_cov_28.087'")
    assert attrs == {"name": "3", "length": 100, "coverage": 28.087, "rc": True}

    attrs = extract_node_len_cov_rc("EDGE_123_length_160_cov_90")
    assert attrs == {"name": "123", "length": 160, "coverage": 90, "rc": False}
