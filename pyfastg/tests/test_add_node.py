import pytest
import networkx as nx
from ..pyfastg import add_node_to_digraph


def test_basic():
    g = nx.DiGraph()
    add_node_to_digraph(
        g, {"name": "asdf", "cov": 5.2, "seq": "ATCGCC", "length": 6}
    )
    assert "asdf" in g.nodes
    assert g.nodes["asdf"]["cov"] == 5.2
    assert g.nodes["asdf"]["gc"] == 4 / 6.0
    assert g.nodes["asdf"]["length"] == 6


def test_insufficient_attrs():
    g = nx.DiGraph()

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {})
    assert "name not present for all nodes" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {"name": "123"})
    assert "length not present for all nodes" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {"name": "123", "length": 2})
    assert "cov not present for all nodes" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {"name": "123", "length": 2, "cov": 6.3})
    assert "seq not present for all nodes" in str(exc_info.value)

    # Finally, this should work
    add_node_to_digraph(
        g, {"name": "123", "length": 2, "cov": 6.3, "seq": "AG"}
    )
    assert "123" in g.nodes


def test_length_mismatch():
    g = nx.DiGraph()
    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(
            g, {"name": "asdf", "cov": 5.2, "seq": "A", "length": 6}
        )
    assert (
        "Length given vs. actual seq. length differs for node asdf"
        in str(exc_info.value)
    )
