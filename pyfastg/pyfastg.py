import re
import networkx as nx
from skbio import DNA


def extract_node_attrs(node_declaration):
    """Returns three specific attributes of a FASTG node declaration.

        As an example, extract_node_attrs("EDGE_3_length_100_cov_28.087'")
        should return {"name": "3-", "length": 100, "cov": 28.087}.

        Parameters
        ----------

        node_declaration: str
            A string tentatively representing a FASTG node declaration (even
            though it should start with ">EDGE"). We impose some pretty strict
            criteria on how this declaration can be structured in order to make
            sure that we get the same amount of information from all nodes in
            the assembly graph file.

        Returns
        -------

        dict
            A mapping of "name", "length", and "cov" to the corresponding
            corresponding node attributes. The "name" value should be a str,
            the "length" value should be an int, the "cov" value should be a
            float, and the "rc" value should be a bool.

        Raises
        ------

        If the regular expression we use to retrieve information from a node
        declaration does not match the node_declaration, this will raise a
        ValueError.
    """
    rc = False
    if node_declaration.endswith("'"):
        rc = True
        node_declaration = node_declaration[0:-1]
    p = re.compile(
        r"EDGE_(?P<node>[a-zA-Z\d]+?)_length_(?P<length>\d+?)_cov_(?P<cov>[\d|\.]+)"  # noqa
    )
    m = p.search(node_declaration)
    if m is None:
        raise ValueError(
            'Wasn\'t able to find all info in the node declaration "{}". '
            "This assembly graph is likely formatted in a way that pyfastg "
            "does not support."
        )
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


def check_all_attrs_present(
    attr_mapping, attrs=("name", "length", "cov", "seq")
):
    """Given a dict of attributes, ensures that all attrs are keys in the dict.

        This is used to make sure all nodes in the graph have the information
        we expect nodes to have.

        Parameters
        ----------

        attr_mapping: dict
            A dictionary describing the attributes of a node.

        attrs: collection
            All of the required attributes that we expect to be present as keys
            in attr_mapping. The defaults for this are pretty reasonable.

        Returns
        -------

        None

        Raises
        ------

        If any of the entries in "attrs" are not present in attr_mapping, this
        will raise a ValueError.
    """
    for required_attr in attrs:
        if required_attr not in attr_mapping:
            raise ValueError(
                "{} not present for all nodes".format(required_attr)
            )


def add_node_to_digraph(digraph, node_attrs):
    """Adds a node (and potentially some edges) to an existing DiGraph object.

        Parameters
        ----------

        digraph: networkx.DiGraph
            An existing graph object, to which a new node and potentially some
            edges will be added.

        node_attrs: dict
            A key/value representation of the attributes the new node will
            have.  This includes everything check_all_attrs_present() looks for
            by default, as well as optionally an "outgoing_node_names"
            attribute: if this outgoing_node_names attribute is present as a
            key in node_attrs, then this will add edges from the new node to
            all of the node IDs contained in node_attrs["outgoing_node_names"].

            (Even if the neighbor nodes referred to in these edges have not
            already been added to the graph, adding an edge referring to them
            should add the corresponding nodes to the graph -- and later on,
            when those nodes are explicitly added to the graph, the edges
            incident on them should remain.)

        Returns
        -------

        None

        Raises
        ------

        ValueError: if the length of the node's seq attribute differs from its
        actual length attribute.
    """
    check_all_attrs_present(node_attrs)
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
    """Given a path to a FASTG file, returns a representation of its structure.

        Parameters
        ----------

        f: str
            Path to a (SPAdes-style) FASTG file.

        Returns
        -------

        networkx.DiGraph
            A representation of the structure of the input assembly graph.
            This DiGraph can immediately be used with the NetworkX library.
    """

    digraph = nx.DiGraph()
    curr_node_attrs = {}

    # First off, quickly validate the file on an initial pass through
    # (TODO: extend this to check up front that all IDs and sequences look
    # valid, etc.?)
    with open(f, "r") as graph_file:
        if not graph_file.readline().startswith(">"):
            raise ValueError(
                "File doesn't start with a \">\" character. This doesn't seem "
                "like a FASTG file."
            )

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
    add_node_to_digraph(digraph, curr_node_attrs)

    # Ensure that all nodes referenced in the digraph have the guaranteed
    # attributes
    # This could *not* happen if node A has an edge to node B, but node B is
    # never actually declared in the FASTG file (i.e. we never see a line that
    # goes ">EDGE_B_..."). This scenario would mean the input graph is invalid.
    for node_name in digraph.nodes:
        check_all_attrs_present(
            digraph.nodes[node_name], attrs=("length", "cov", "gc")
        )

    return digraph
