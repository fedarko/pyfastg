from ..pyfastg import FASTGNode


def test_fastgnode_basic():
    n = FASTGNode("supercoolnode", 100, 98.9, True)
    assert n.name == "supercoolnode-"
    assert n.length == 100
    assert n.cov == 98.9
    assert n.reverse_complemented == True

    n_str = n.__str__().split("\n")
    assert n_str[0] == "Node: supercoolnode-"
    assert n_str[1] == "Length: 100"
    assert n_str[2] == "Coverage: 98.9"
    assert n_str[3] == "Reverse Complemented?: True"
