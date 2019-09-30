import re
import networkx as nx


class FASTGNode(object):
    """Simple class that stores information about a node in the graph."""

    def __init__(self, name, length, cov, reverse_complemented):
        if reverse_complemented:
            self.name = name + "-"
        else:
            self.name = name + "+"
        self.length = length
        self.cov = cov
        self.reverse_complemented = reverse_complemented

    def __str__(self):
        return str(
            "Node: {0}\nLength: {1}\nCoverage: {2}\nReverse Complemented?: {3}"
        ).format(
            self.name, self.length, self.cov, str(self.reverse_complemented)
        )


def extract_node_len_cov_rc(node_declaration):
    """Returns four specific attributes of a FASTG node declaration.

    As an example, extract_node_len_cov_rc("EDGE_3_length_100_cov_28.087'")
    should return {"name": "3", "length": 100, "coverage": 28.087, "rc": True}.

    Returns
    -------

    dict
        A mapping of "name", "length", "coverage", and "rc" to the
        corresponding node attributes. The "name" value should be a str,
        the "length" value should be an int, the "coverage" value should be a
        float, and the "rc" value should be a bool.
    """
    rc = False
    if node_declaration.endswith("'"):
        rc = True
        node_declaration = node_declaration[0:-1]
    p = re.compile(
        r"EDGE_(?P<node>\d*?)_length_(?P<length>\d*?)_cov_(?P<cov>[\d|\.]*)"
    )
    m = p.search(node_declaration)
    return {
        "name": m.group("node"),
        "length": int(m.group("length")),
        "coverage": float(m.group("cov")),
        "rc": rc,
    }


def make_node(declaration):
    attrs = extract_node_len_cov_rc(declaration)
    return FASTGNode(
        attrs["name"], attrs["length"], attrs["coverage"], attrs["rc"]
    )


def parse_fastg(f):
    # This is a 2-D list where each entry in the list corresponds to
    # [starting node name, either None or a comma-separated list of node names]
    # The second entry indicates the nodes that the starting node name has an
    # outgoing edge to.
    # A TODO here is storing node sequence information -- or at least recording
    # GC content, etc.
    node_neighs = []
    with open(f, "r") as graph_file:
        for line in graph_file:
            if line.startswith(">"):
                line_no_sc = line.strip().strip(";")
                colons = sum([1 for x in line_no_sc if x == ":"])
                if colons > 1:
                    raise ValueError(
                        "multiple ':'s found in line, and can only "
                        + "be used to separate nodes from neighbor "
                        + "list\n"
                    )
                elif colons == 0:
                    # orphaned node or terminal node
                    node_neighs.append([line_no_sc, None])
                else:
                    node_neighs.append(line_no_sc.split(":"))

    # Convert node_neighs to a NetworkX DiGraph
    digraph = nx.DiGraph()
    for node, neighs in node_neighs:
        node_obj = make_node(node)
        # We don't bother storing the reverse_complemented FASTGNode attribute
        # since that information is now encoded in the node name (if it ends in
        # "+" it's not reverse-complemented; if it ends in "-" then it is)
        digraph.add_node(
            node_obj.name, length=node_obj.length, coverage=node_obj.cov
        )
        if neighs is not None:
            for outgoing_neighbor_node in neighs.split(","):
                neighbor_node_obj = make_node(outgoing_neighbor_node)
                digraph.add_edge(node_obj.name, neighbor_node_obj.name)

    return digraph
