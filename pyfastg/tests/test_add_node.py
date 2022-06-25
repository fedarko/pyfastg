import pytest
import networkx as nx
from ..pyfastg import add_node_to_digraph


def test_basic():
    def check_asdf(g):
        assert "asdf" in g.nodes
        assert g.nodes["asdf"]["cov"] == 5.2
        assert g.nodes["asdf"]["gc"] == 4 / 6.0
        assert g.nodes["asdf"]["length"] == 6

    def check_ghjasdf(g):
        assert "ghjasd" in g.nodes
        assert g.nodes["ghjasd"]["cov"] == 100
        assert g.nodes["ghjasd"]["gc"] == 1 / 3.0
        assert g.nodes["ghjasd"]["length"] == 3

    g = nx.DiGraph()

    # 1. Add node "asdf" to g
    add_node_to_digraph(
        g, {"name": "asdf", "cov": 5.2, "seq": "ATCGCC", "length": 6}
    )
    check_asdf(g)

    # 2. Add node "ghjasd" to g
    add_node_to_digraph(
        g,
        {
            "name": "ghjasd",
            "cov": 100,
            "seq": "CAT",
            "length": 3,
            "outgoing_node_names": ["asdf", "qwerty", "hamborgar"],
        },
    )
    # This should have added three new nodes (ghjasdf, qwerty, hamborgar)
    # qwerty and hamborgar, however, don't have any attributes (yet)

    # Double-check that asdf's attributes were not somehow lost
    check_asdf(g)
    check_ghjasdf(g)
    assert "qwerty" in g.nodes
    assert "hamborgar" in g.nodes

    assert ("ghjasd", "asdf") in g.edges
    assert ("ghjasd", "qwerty") in g.edges
    assert ("ghjasd", "hamborgar") in g.edges

    # 3. Add node "hamborgar" to g (it's already in there but is "empty")
    add_node_to_digraph(
        g, {"name": "hamborgar", "cov": 33.3, "seq": "AAAA", "length": 4}
    )
    # Again, check that prior nodes' attributes are ok
    check_asdf(g)
    check_ghjasdf(g)

    assert "qwerty" in g.nodes
    assert "hamborgar" in g.nodes
    assert ("ghjasd", "asdf") in g.edges
    assert ("ghjasd", "qwerty") in g.edges
    assert ("ghjasd", "hamborgar") in g.edges
    assert g.nodes["hamborgar"]["cov"] == 33.3
    assert g.nodes["hamborgar"]["gc"] == 0
    assert g.nodes["hamborgar"]["length"] == 4


def test_insufficient_attrs():
    g = nx.DiGraph()

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {})
    assert "name not present for all edges" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {"name": "123"})
    assert "length not present for all edges" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {"name": "123", "length": 2})
    assert "cov not present for all edges" in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        add_node_to_digraph(g, {"name": "123", "length": 2, "cov": 6.3})
    assert "seq not present for all edges" in str(exc_info.value)

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
    assert "Length given vs. actual seq. length differs for edge asdf" in str(
        exc_info.value
    )
