import pytest
import tempfile
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


def test_parse_multicolon():
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/multicolon.fastg")
    first_bad_line = (
        ">EDGE_1_length_9_cov_4.5:EDGE_3_length_5_cov_16.5':"
        "EDGE_2_length_3_cov_100;"
    )
    assert str(exc_info.value) == (
        'Multiple ":" characters found in line "{}". An edge line can only '
        'have exactly 0 or 1 ":" characters; pyfastg doesn\'t support FASTG '
        '"properties" yet.'
    ).format(first_bad_line)


def test_parse_stuff_at_start():
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/stuff_at_start.fastg")
    assert 'File doesn\'t start with a ">" character' in str(exc_info.value)

    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/whitespace_at_start.fastg")
    assert 'File doesn\'t start with a ">" character' in str(exc_info.value)


def test_parse_no_rc():
    """Tests a graph where 1+ -> 3-, but not 3+ -> 1-, exists.

    Ensures that this second edge isn't automatically created.

    Also, this graph has a case where 2+ exists but not but 2-,
    and where 4- exists but not 4+. Both these cases are ok, and
    we verify that these reverse-complement nodes are also not automatically
    created.
    """
    digraph = parse_fastg("pyfastg/tests/input/norc.fastg")
    assert len(digraph.nodes) == 6
    assert len(digraph.edges) == 2

    assert set(digraph.nodes) == {"1+", "1-", "2+", "3+", "3-", "4-"}
    assert ("1+", "3-") in digraph.edges
    assert ("1-", "2+") in digraph.edges


def test_parse_no_seq():
    """Tests a graph where an edge isn't defined.

    This is kind of like an integration test that verifies that the
    call to check_all_attrs_present() after creating the graph is good.
    """
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/noseq.fastg")
    assert str(exc_info.value) == "length not present for all edges"


def test_parse_bad_len():
    """Tests a graph where the sequence length != the declared length.

    This is kind of like an integration test that verifies that the
    call to check_all_attrs_present() after creating the graph is good.
    """
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/bad_len.fastg")
    assert str(exc_info.value) == (
        "Length given vs. actual seq. length differs for edge 1-"
    )


def test_parse_no_semicolon():
    """Tests a graph where the edge declaration and seq are on the same line.

    Technically allowed by the FASTG spec, but unsupported by us, for the very
    well-thought-out reason of "I don't want to refactor this code."
    """
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/nosc.fastg")
    assert str(exc_info.value) == (
        "The edge declaration line \">EDGE_1_length_9_cov_4.5':"
        'EDGE_2_length_3_cov_100; ATGGGCGAT" must end with a ; character. '
        "(The sequence for this edge should be given on the next line.)"
    )

    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/nosc2.fastg")
    assert str(exc_info.value) == (
        'The edge declaration line ">EDGE_1_length_9_cov_4.5\'" '
        "must end with a ; character. "
        "(The sequence for this edge should be given on the next line.)"
    )


def test_parse_inconsistent_declarations():
    """Tests graphs with inconsistent edge declarations."""
    # i know that having separate files for this -- rather than just mocking
    # the files with tempfiles or something -- is bad practice, but this way is
    # quickest to set up right now...
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/inconsistent_cov.fastg")
    assert str(exc_info.value) == (
        "Node 3- has inconsistent edge declarations in the FASTG: we already "
        'saw "EDGE_3_length_5_cov_16.6\'", but we just saw '
        '"EDGE_3_length_5_cov_16.5\'".'
    )

    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/inconsistent_len.fastg")
    assert str(exc_info.value) == (
        "Node 1+ has inconsistent edge declarations in the FASTG: we already "
        'saw "EDGE_1_length_20_cov_4.5", but we just saw '
        '"EDGE_1_length_9_cov_4.5".'
    )

    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/inconsistent_both.fastg")
    assert str(exc_info.value) == (
        "Node 1+ has inconsistent edge declarations in the FASTG: we already "
        'saw "EDGE_1_length_20_cov_5.2", but we just saw '
        '"EDGE_1_length_9_cov_4.5".'
    )

    # This one is fun -- what if a declaration occurs twice on the same line,
    # but it's inconsistent?
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/inconsistent_on_same_line.fastg")
    assert str(exc_info.value) == (
        "Node 1+ has inconsistent edge declarations in the FASTG: we already "
        'saw "EDGE_1_length_20_cov_5.2", but we just saw '
        '"EDGE_1_length_9_cov_20".'
    )


def test_parse_multi_graph():
    """Tests a graph with multiple copies of a single adjacency."""
    with pytest.raises(ValueError) as exc_info:
        parse_fastg("pyfastg/tests/input/multi_graph.fastg")
    assert str(exc_info.value) == (
        "Node 1+ has duplicate outgoing adjacencies."
    )


def write_tempfile(text):
    # yanked from strainflye/tests/utils_for_testing.py
    fh = tempfile.NamedTemporaryFile(suffix=".fastg")
    with open(fh.name, "w") as f:
        f.write(text)
    return fh


def test_parse_whitespace_btwn_adjs():
    # " , "
    fh = write_tempfile(
        """>EDGE_1_length_9_cov_4.5:EDGE_3_length_5_cov_16.5' , EDGE_2_length_3_cov_1;
ATCGCCCAT
>EDGE_1_length_9_cov_4.5':EDGE_2_length_3_cov_1';
ATGGGCGAT
>EDGE_2_length_3_cov_1:EDGE_1_length_9_cov_4.5;
CGA
>EDGE_2_length_3_cov_1';
TCG
>EDGE_3_length_5_cov_16.5:EDGE_1_length_9_cov_4.5',EDGE_2_length_3_cov_1';
GGATC
>EDGE_3_length_5_cov_16.5':EDGE_2_length_3_cov_1';
GATCC
"""
    )
    with pytest.raises(ValueError) as ei:
        parse_fastg(fh.name)
    assert str(ei.value) == (
        "Wasn't able to find all expected info (edge name, length, coverage) "
        'in the declaration "EDGE_3_length_5_cov_16.5\' ". Please remember '
        "that pyfastg only supports SPAdes-dialect FASTG files."
    )

    # " ,"
    fh = write_tempfile(
        """>EDGE_1_length_9_cov_4.5:EDGE_3_length_5_cov_16.5' ,EDGE_2_length_3_cov_1;
ATCGCCCAT
>EDGE_1_length_9_cov_4.5':EDGE_2_length_3_cov_1';
ATGGGCGAT
>EDGE_2_length_3_cov_1:EDGE_1_length_9_cov_4.5;
CGA
>EDGE_2_length_3_cov_1';
TCG
>EDGE_3_length_5_cov_16.5:EDGE_1_length_9_cov_4.5',EDGE_2_length_3_cov_1';
GGATC
>EDGE_3_length_5_cov_16.5':EDGE_2_length_3_cov_1';
GATCC
"""
    )
    with pytest.raises(ValueError) as ei:
        parse_fastg(fh.name)
    assert str(ei.value) == (
        "Wasn't able to find all expected info (edge name, length, coverage) "
        'in the declaration "EDGE_3_length_5_cov_16.5\' ". Please remember '
        "that pyfastg only supports SPAdes-dialect FASTG files."
    )

    # ", "
    fh = write_tempfile(
        """>EDGE_1_length_9_cov_4.5:EDGE_3_length_5_cov_16.5', EDGE_2_length_3_cov_1;
ATCGCCCAT
>EDGE_1_length_9_cov_4.5':EDGE_2_length_3_cov_1';
ATGGGCGAT
>EDGE_2_length_3_cov_1:EDGE_1_length_9_cov_4.5;
CGA
>EDGE_2_length_3_cov_1';
TCG
>EDGE_3_length_5_cov_16.5:EDGE_1_length_9_cov_4.5',EDGE_2_length_3_cov_1';
GGATC
>EDGE_3_length_5_cov_16.5':EDGE_2_length_3_cov_1';
GATCC
"""
    )
    with pytest.raises(ValueError) as ei:
        parse_fastg(fh.name)
    assert str(ei.value) == (
        "Wasn't able to find all expected info (edge name, length, coverage) "
        'in the declaration " EDGE_2_length_3_cov_1". Please remember '
        "that pyfastg only supports SPAdes-dialect FASTG files."
    )


def test_parse_tilde_before_outgoing_adj_name():
    """These are technically supported in the spec, but we don't allow them."""
    fh = write_tempfile(
        """>EDGE_1_length_9_cov_4.5:EDGE_3_length_5_cov_16.5',~EDGE_2_length_3_cov_1;
ATCGCCCAT
>EDGE_1_length_9_cov_4.5':EDGE_2_length_3_cov_1';
ATGGGCGAT
>EDGE_2_length_3_cov_1:EDGE_1_length_9_cov_4.5;
CGA
>EDGE_2_length_3_cov_1';
TCG
>EDGE_3_length_5_cov_16.5:EDGE_1_length_9_cov_4.5',EDGE_2_length_3_cov_1';
GGATC
>EDGE_3_length_5_cov_16.5':EDGE_2_length_3_cov_1';
GATCC
"""
    )
    with pytest.raises(ValueError) as ei:
        parse_fastg(fh.name)
    assert str(ei.value) == (
        "pyfastg does not support the ~ notation described in the FASTG spec."
    )


def test_parse_brackets_after_outgoing_adj_name():
    """Again, technically supported in the spec, but we don't allow this.

    This is a way of encoding adjacency-specific properties."""
    fh = write_tempfile(
        """>EDGE_1_length_9_cov_4.5:EDGE_3_length_5_cov_16.5'[prop=val],EDGE_2_length_3_cov_1;
ATCGCCCAT
>EDGE_1_length_9_cov_4.5':EDGE_2_length_3_cov_1';
ATGGGCGAT
>EDGE_2_length_3_cov_1:EDGE_1_length_9_cov_4.5;
CGA
>EDGE_2_length_3_cov_1';
TCG
>EDGE_3_length_5_cov_16.5:EDGE_1_length_9_cov_4.5',EDGE_2_length_3_cov_1';
GGATC
>EDGE_3_length_5_cov_16.5':EDGE_2_length_3_cov_1';
GATCC
"""
    )
    with pytest.raises(ValueError) as ei:
        parse_fastg(fh.name)
    assert str(ei.value) == (
        "pyfastg does not support the [] notation described in the FASTG spec."
    )
