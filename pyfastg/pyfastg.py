import re
import networkx as nx


def extract_node_attrs(declaration):
    """Returns three specific attributes of a FASTG edge declaration.

    As an example, extract_node_attrs("EDGE_3_length_100_cov_28.087'")
    should return {"name": "3-", "length": 100, "cov": 28.087}.

    The function name might be confusing, sorry. This parses a declared edge in
    the FASTG file, but we will convert it to a node in the returned NetworkX
    graph (and arguably "edges" in the FASTG spec are really just nodes,
    anyway).

    Parameters
    ----------
    declaration: str
        A string tentatively representing a FASTG edge declaration.
        We impose some pretty strict criteria on how this can be structured in
        order to make sure that we get the same amount of information from all
        edges in the assembly graph file.

    Returns
    -------
    dict
        A mapping of "name", "length", and "cov" to this edge's attributes.
        The "name" value will be a str, the "length" value will be an int,
        and the "cov" value will be a float.

    Raises
    ------
    ValueError
        If the regular expression we use to retrieve information from a
        declaration does not match our regex, or is explicitly not supported
        by pyfastg in some way.

        If trying to parse the length or coverage values from the declaration
        fails.
    """
    # Quick early checks
    if declaration.startswith("~"):
        raise ValueError(
            "pyfastg does not support the ~ notation described in the FASTG "
            "spec."
        )
    if "[" in declaration:
        raise ValueError(
            "pyfastg does not support the [] notation described in the FASTG "
            "spec."
        )

    rc = False
    if declaration.endswith("'"):
        rc = True
        nonrc_declaration = declaration[:-1]
    else:
        nonrc_declaration = declaration
    # Regarding the (?P<name>...) syntax, see
    # https://docs.python.org/3/library/re.html#index-17. This is a way of
    # saying "capture this group, and also call it 'name'" -- this makes it
    # easy to extract the edge name, length, etc. from the declaration.
    p = re.compile(
        r"^EDGE_(?P<name>[a-zA-Z\d]+)_length_(?P<length>\d+)_cov_(?P<cov>[\d\.]+)$"  # noqa
    )
    m = p.search(nonrc_declaration)
    if m is None:
        raise ValueError(
            (
                "Wasn't able to find all expected info (edge name, length, "
                'coverage) in the declaration "{}". Please remember that '
                "pyfastg only supports SPAdes-dialect FASTG files."
            ).format(declaration)
        )
    name = m.group("name")
    if rc:
        name += "-"
    else:
        name += "+"

    # These will trigger errors if the length or coverage are invalid. In
    # practice, the length should be fine (since we know that it will be a
    # string comprised of just digits).
    slen = int(m.group("length"))
    # However, since the coverage regex allows for the . character (since
    # coverage is usually a float), we could have problematic cases slip
    # through the regex -- for example, "cov_....." or something silly.
    # We could modify the coverage regex to be fancier and only detect numbers
    # with a single ".", but it's easier (and simpler) to just let Python's
    # float() function do the dirty work for us.
    scov = float(m.group("cov"))
    return {"name": name, "length": slen, "cov": scov}


def check_all_attrs_present(
    attr_mapping, attrs=("name", "length", "cov", "seq")
):
    """Given a dict of attributes, ensures that all attrs are keys in the dict.

    This is used to make sure all nodes (representing FASTG edges) in the graph
    have the information we expect them to have.

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
    ValueError
        If any of the entries in "attrs" are not present in attr_mapping.
    """
    for required_attr in attrs:
        if required_attr not in attr_mapping:
            raise ValueError(
                "{} not present for all edges".format(required_attr)
            )


def check_unique(values, err_msg):
    """Checks that a given list contains all unique values.

    Parameters
    ----------
    values: list

    err_msg: str

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If "values" is not unique. (This particular error's message will be
        err_msg.)
    """
    if len(set(values)) < len(values):
        raise ValueError(err_msg)


def validate_seq_and_compute_gc(seq):
    """Validates a sequence and computes its GC content.

    This replaces our previous dependency on scikit-bio.DNA() to do this.

    Parameters
    ----------
    seq: str
        Represents a DNA or RNA sequence.

    Returns
    -------
    gc: float
        GC content. Will be within the range [0, 1].
        If the length of seq is zero, the returned GC content will be zero.

    Raises
    ------
    ValueError
        If the sequence does not only contain characters from the alphabet
        {A, C, G, T, U}.

    References
    ----------
    Zero-length sequences are ok, per the FASTG spec:
    http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf

    Returning a GC content of zero for zero-length sequences is consistent
    with the behavior of scikit-bio and of BioPython:
    http://scikit-bio.org/docs/latest/generated/skbio.sequence.DNA.gc_content.html
    https://biopython.org/docs/dev/api/Bio.SeqUtils.html#Bio.SeqUtils.GC
    """
    nongc_chars = "ATU"
    gc_chars = "GC"
    gc_ct = 0
    seq_len = 0
    for char in seq:
        if char in gc_chars:
            gc_ct += 1
        elif char not in nongc_chars:
            raise ValueError(
                (
                    'Sequence "{}" contains character(s) not in the alphabet '
                    "{{A, C, G, T, U}}."
                ).format(seq)
            )
        seq_len += 1

    if seq_len == 0:
        return 0
    else:
        return gc_ct / seq_len


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
            "Length given vs. actual seq. length differs for edge {}".format(
                node_attrs["name"]
            )
        )
    node_gc_content = validate_seq_and_compute_gc(node_attrs["seq"])
    node_name = node_attrs["name"]
    digraph.add_node(
        node_name,
        length=node_attrs["length"],
        cov=node_attrs["cov"],
        gc=node_gc_content,
    )
    if "outgoing_node_names" in node_attrs:
        outgoing_nodes = node_attrs["outgoing_node_names"]
        check_unique(
            outgoing_nodes,
            "Node {} has duplicate outgoing adjacencies.".format(node_name),
        )
        for neighbor_node_name in node_attrs["outgoing_node_names"]:
            digraph.add_edge(node_name, neighbor_node_name)


def update_and_check_decl(node_name, declaration, nodename2declaration):
    """Updates and checks a dict regarding declaration consistency.

    Here's what's going on. An edge's full "declaration" (e.g.
    "EDGE_1_length_9909_cov_6.94721") contains multiple types of information
    besides the edge's actual identifier: it also describes length and
    coverage. An edge's declaration can be repeated multiple times in the
    FASTG file (whenever another edge has an outgoing adjacency to this edge,
    the declaration gets repeated; this is because, from the FASTG spec's
    perspective, this bulky declaration is really the full name of the edge).

    However, to us, this declaration isn't the edge's name. The edge's name
    (well, really the *node*'s name in our output NetworkX graph) will be based
    on just some of the edge's info -- "EDGE_1_length_9909_cov_6.94721" becomes
    just "1+", with the other information stored as attributes of this node.

    The problem we want to address here is *inconsistency*. What do you do with
    a malformed graph, where we later see a reference to an edge with the
    declaration "EDGE_1_length_9909_cov_7"? According to the FASTG spec, this
    is a completely different edge, since its name is different; but from our
    perspective it's still "1+", right?

    This function addresses this problem. The first time we see an edge's
    declaration (whether it is on the line that this edge is actually declared,
    or if it's just because another edge had an outgoing adjacency to this
    edge), we store the exact declaration we see. Every time we see this edge's
    declaration again, we verify that the declarations are consistent.

    ... This is probably overkill, but better safe than sorry.

    Parameters
    ----------
    node_name: str
        The name of a node in the NetworkX graph. Something like "1+".

    declaration: str
        The full declaration of this node's corresponding edge in the FASTG
        file. Something like "EDGE_1_length_9909_cov_6.942721".

    nodename2declaration: dict
        Maps node names to their seen declarations.

        If node_name is not already a key in this dict, we will add node_name
        to this dict (with declaration as the value to which node_name maps).

        If node_name is already a key in this dict, we will check that
        declaration matches the existing declaration stored for node_name in
        this dict.

    Returns
    -------
    None

    Raises
    ------
    ValueError
        If node_name is in nodename2declaration, but
        nodename2declaration[node_name] != declaration (i.e. the declarations
        are inconsistent).
    """
    if node_name in nodename2declaration:
        existing_declaration = nodename2declaration[node_name]
        if existing_declaration != declaration:
            raise ValueError(
                (
                    "Node {} has inconsistent edge declarations in the FASTG: "
                    'we already saw "{}", but we just saw "{}".'
                ).format(node_name, existing_declaration, declaration)
            )
    else:
        nodename2declaration[node_name] = declaration


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

    # Maps the name of a node in our graph (e.g. "1+") to its declaration in
    # the FASTG file (e.g. "EDGE_1_length_100_cov_20").
    # Since a node's declaration could be repeated multiple times (e.g. other
    # edges have outgoing adjacencies to it), we use this dict to ensure that
    # all instances of a node's declaration are identical. This lets us catch
    # subtle inconsistencies (e.g. let's say we later see a reference to
    # "EDGE_1_length_101_cov_20"; why does this have a different length?).
    nodename2decl = {}
    with open(f, "r") as graph_file:
        for line in graph_file:
            stripped_line = line.strip()
            if stripped_line.startswith(">"):
                if len(curr_node_attrs) > 0:
                    add_node_to_digraph(digraph, curr_node_attrs)
                    curr_node_attrs = {}

                if not stripped_line.endswith(";"):
                    raise ValueError(
                        (
                            'The edge declaration line "{}" must end with a ; '
                            "character. (The sequence for this edge should be "
                            "given on the next line.)"
                        ).format(stripped_line)
                    )

                line_no_sc = stripped_line[:-1]
                colon_ct = line_no_sc.count(":")
                if colon_ct > 1:
                    raise ValueError(
                        (
                            'Multiple ":" characters found in line "{}". An '
                            'edge line can only have exactly 0 or 1 ":" '
                            "characters; pyfastg doesn't support FASTG "
                            '"properties" yet.'
                        ).format(stripped_line)
                    )
                elif colon_ct == 0:
                    # orphaned node or terminal node
                    # The [1:] slices off the starting > character
                    curr_node_decl = line_no_sc[1:]
                    curr_node_attrs = extract_node_attrs(curr_node_decl)
                    update_and_check_decl(
                        curr_node_attrs["name"], curr_node_decl, nodename2decl
                    )
                else:
                    # This node has at least one outgoing edge
                    line_no_sc_split = line_no_sc.split(":")
                    # The [1:] slices off the starting > character
                    curr_node_decl = line_no_sc_split[0][1:]
                    curr_node_attrs = extract_node_attrs(curr_node_decl)
                    update_and_check_decl(
                        curr_node_attrs["name"], curr_node_decl, nodename2decl
                    )

                    # Check outgoing edges
                    curr_node_attrs["outgoing_node_names"] = []
                    for neighbor_decl in line_no_sc_split[1].split(","):
                        neighbor_attrs = extract_node_attrs(neighbor_decl)
                        neighbor_name = neighbor_attrs["name"]
                        update_and_check_decl(
                            neighbor_name, neighbor_decl, nodename2decl
                        )
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
