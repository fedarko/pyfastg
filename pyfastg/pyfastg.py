import re
import networkx as nx
from skbio import DNA


def extract_node_attrs(node_declaration):
    """Returns three specific attributes of a FASTG node declaration.

    As an example, extract_node_attrs("EDGE_3_length_100_cov_28.087'")
    should return {"name": "3-", "length": 100, "cov": 28.087}.

    Returns
    -------

    dict
        A mapping of "name", "length", and "cov" to the corresponding
        corresponding node attributes. The "name" value should be a str,
        the "length" value should be an int, the "cov" value should be a
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
    name = m.group("node")
    if rc:
        name += "-"
    else:
        name += "+"
    return {
        "name": name,
        "length": int(m.group("length")),
        "cov": float(m.group("cov")),
    }


def add_node_to_digraph(digraph, node_attrs):
    for required_attr in ("name", "length", "cov", "seq"):
        if required_attr not in node_attrs:
            raise ValueError(
                "{} not present for all nodes".format(required_attr)
            )
    if len(node_attrs["seq"]) != node_attrs["length"]:
        raise ValueError(
            "Length given vs. actual seq. length differs for node {}".format(
                node_attrs["name"]
            )
        )
    node_gc_content = DNA(node_attrs["seq"]).gc_content()
    digraph.add_node(
        node_attrs["name"],
        length=node_attrs["length"],
        cov=node_attrs["cov"],
        gc=node_gc_content,
    )
    if "outgoing_node_names" in node_attrs:
        for neighbor_node_name in node_attrs["outgoing_node_names"]:
            digraph.add_edge(node_attrs["name"], neighbor_node_name)


def parse_fastg(f):

    digraph = nx.DiGraph()
    curr_node_attrs = {}
    with open(f, "r") as graph_file:
        for line in graph_file:
            stripped_line = line.strip()
            if stripped_line.startswith(">"):
                if len(curr_node_attrs) > 0:
                    add_node_to_digraph(digraph, curr_node_attrs)
                    curr_node_attrs = {}
                line_no_sc = stripped_line.strip(";")
                colons = sum([1 for x in line_no_sc if x == ":"])
                if colons > 1:
                    raise ValueError(
                        "multiple ':'s found in line, and can only "
                        + "be used to separate nodes from neighbor "
                        + "list\n"
                    )
                elif colons == 0:
                    # orphaned node or terminal node
                    curr_node_attrs = extract_node_attrs(line_no_sc)
                else:
                    # This node has at least one outgoing edge
                    line_no_sc_split = line_no_sc.split(":")
                    curr_node_attrs = extract_node_attrs(line_no_sc_split[0])
                    curr_node_attrs["outgoing_node_names"] = []
                    for neighbor_decl in line_no_sc_split[1].split(","):
                        neighbor_name = extract_node_attrs(neighbor_decl)[
                            "name"
                        ]
                        curr_node_attrs["outgoing_node_names"].append(
                            neighbor_name
                        )
            elif len(stripped_line) > 0:
                if "seq" in curr_node_attrs:
                    curr_node_attrs["seq"] += stripped_line
                else:
                    curr_node_attrs["seq"] = stripped_line
    # Account for the last node described at the bottom of the file
    if len(curr_node_attrs) > 0:
        add_node_to_digraph(digraph, curr_node_attrs)

    return digraph
