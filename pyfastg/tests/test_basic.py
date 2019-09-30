from pyfastg import parse_fastg

def test_simple_assembly_graph():
    graph = parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
    assert len(graph) == 66
