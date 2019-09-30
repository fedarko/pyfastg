from pyfastg import parse_fastg

def test_parse_simple_assembly_graph():
    digraph = parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
    assert len(digraph.nodes) == 66
    assert len(digraph.edges) == 86
    for i in range(1, 34):
        si = str(i)
        assert si + "+" in digraph.nodes
        assert si + "-" in digraph.nodes
