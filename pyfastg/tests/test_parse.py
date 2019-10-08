import pytest
from pyfastg import parse_fastg


def test_parse_medium_assembly_graph():
    digraph = parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
    assert len(digraph.nodes) == 66
    assert len(digraph.edges) == 86
    for i in range(1, 34):
        si = str(i)
        for suffix in ("+", "-"):
            name = si + suffix
            assert name in digraph.nodes
            # As a sanity check, make sure that node 1's attrs were all parsed
            # successfully
            # (Ideally we'd do this automatically for large-ish graphs like
            # this one, but ... that would require writing another FASTG parser
            # here :)
            if i == 1:
                assert digraph.nodes[name]["length"] == 9909
                assert digraph.nodes[name]["cov"] == 6.94721
                assert digraph.nodes[name]["gc"] == 5153 / 9909.0


def test_parse_small_assembly_graph():
    """Tests a simple manually-created assembly graph.

    Also tests an identical graph that happens to have a bunch of blank lines
    around the sequences (shouldn't cause a problem).
    """
    for fn in ("small", "whitespace_in_seq"):
        digraph = parse_fastg("pyfastg/tests/input/{}.fastg".format(fn))
        assert len(digraph.nodes) == 6
        assert len(digraph.edges) == 8
        i2length = {1: 9, 2: 3, 3: 5}
        i2cov = {1: 4.5, 2: 100, 3: 16.5}
        i2gc = {1: 5 / 9.0, 2: 2 / 3.0, 3: 3 / 5.0}
        for i in range(1, 4):
            si = str(i)
            for suffix in ("+", "-"):
                name = si + suffix
                assert name in digraph.nodes
                assert digraph.nodes[name]["cov"] == i2cov[i]
                assert digraph.nodes[name]["length"] == i2length[i]
                assert digraph.nodes[name]["gc"] == i2gc[i]

        valid_edges = (
            ("2+", "1+"),
            ("2+", "3-"),
            ("2+", "3+"),
            ("1+", "3-"),
            ("3-", "2-"),
            ("3+", "1-"),
            ("3+", "2-"),
            ("1-", "2-"),
        )
        for e in valid_edges:
            assert e in digraph.edges


def test_parse_multicolon_assembly_graph():
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/multicolon.fastg")
    assert "multiple ':'s found in line" in str(exc_info.value)


def test_parse_stuff_at_start_assembly_graph():
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/stuff_at_start.fastg")
    assert 'File doesn\'t start with a ">" character' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/whitespace_at_start.fastg")
    assert 'File doesn\'t start with a ">" character' in str(exc_info.value)
