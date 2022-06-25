# pyfastg: a minimal Python library for parsing networks from SPAdes FASTG files
[![pyfastg CI](https://github.com/fedarko/pyfastg/actions/workflows/main.yml/badge.svg)](https://github.com/fedarko/pyfastg/actions/workflows/main.yml)
[![Code Coverage](https://codecov.io/gh/fedarko/pyfastg/branch/master/graph/badge.svg)](https://codecov.io/gh/fedarko/pyfastg)
[![PyPI](https://img.shields.io/pypi/v/pyfastg)](https://pypi.org/project/pyfastg)

## The FASTG file format
FASTG is a format to describe genome assemblies, geared toward accurately representing the ambiguity resulting from sequencing limitations, ploidy, or other factors that complicate representation of a seqence as a simple string.  The official spec for the FASTG format can be found [here](http://fastg.sourceforge.net/).

pyfastg parses graphs that follow **a subset of this specification**: in
particular, it is designed to work with files output by the
[SPAdes](http://cab.spbu.ru/software/spades/) family of assemblers.

## The pyfastg library
pyfastg contains `parse_fastg()`, a function that accepts as input a path
to a SPAdes FASTG file. This function parses the structure of the specified
file, returning a [NetworkX](https://networkx.github.io) `DiGraph` object representing
the structure of the graph.

pyfastg is very much in its infancy, so it may be most useful as a starting point.
Pull requests are welcome!

## Note about the graph topology

Version 1.00 of the [FASTG spec](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf) contains
the following sentence (in section 6, page 7):

> Note also that strictly speaking, [the structure described in a FASTG file] is not a graph at all, as we have not specified a notion of vertex. However in many cases one can without ambiguity define vertices and thereby associate a _bona fide_ digraph, and we do so frequently in this document to illustrate concepts.

We do this in pyfastg. **"Edges" in the FASTG file will be represented as nodes
in the NetworkX graph, and "adjacencies" between edges in the FASTG file will
be represented as edges in the NetworkX graph.** As far as we're aware, this is
usually how these graphs are visualized.

### Installation
pyfastg can be installed using pip:

```bash
pip install pyfastg
```

### Quick Example
The second line (which points to one of pyfastg's test assembly graphs)
assumes that you're located in the root directory of the pyfastg repo.

```python
>>> import pyfastg
>>> g = pyfastg.parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
>>> # g is now a NetworkX DiGraph! We can do whatever we want with this object.
>>> # Example: List the sequences in this graph (these are "edges" in the FASTG
>>> # file, but are represented as nodes in g)
>>> g.nodes()
NodeView(('1+', '29-', '1-', '6-', '2+', '26+', '27+', '2-', '3+', '4+', '6+', '7+', '3-', '33-', '9-', '4-', '5+', '5-', '28+', '7-', '8+', '28-', '9+', '8-', '12-', '10+', '12+', '10-', '24-', '32-', '11+', '30-', '11-', '27-', '19-', '13+', '25+', '31-', '13-', '14+', '14-', '26-', '15+', '15-', '23-', '16+', '16-', '17+', '17-', '19+', '18+', '33+', '18-', '20+', '20-', '22+', '21+', '21-', '22-', '23+', '24+', '25-', '29+', '30+', '31+', '32+'))
>>> # Example: Get details for a single sequence (length, coverage, GC-content)
>>> g.nodes["15+"]
{'length': 193, 'cov': 6.93966, 'gc': 0.5492227979274611}
>>> # Example: Get information about the graph's connectivity
>>> import networkx as nx
>>> components = list(nx.weakly_connected_components(g))
>>> for c in components:
...     print(len(c), "nodes")
...     print(c)
...
33 nodes
{'8-', '17-', '15+', '30+', '16+', '26-', '25+', '19+', '7+', '23+', '14-', '18-', '10-', '29-', '20-', '27-', '11-', '5-', '3+', '2-', '12-', '13+', '31-', '6+', '1+', '21-', '24-', '32-', '22+', '28+', '4+', '33-', '9-'}
33 nodes
{'26+', '29+', '18+', '3-', '2+', '8+', '15-', '24+', '9+', '17+', '27+', '28-', '11+', '6-', '20+', '14+', '19-', '13-', '4-', '21+', '5+', '31+', '22-', '12+', '25-', '30-', '10+', '1-', '7-', '32+', '23-', '33+', '16-'}
```

### Required File Format (tl;dr: SPAdes-dialect FASTG files only)
Currently, pyfastg is hardcoded to parse FASTG files created by the SPAdes assembler. Other valid FASTG files that don't follow the pattern used by SPAdes for edge names are not supported.

In particular, each edge in the file must have a name formatted like:

```bash
EDGE_1_length_9909_cov_6.94721
```

The edge ID (here, `1`) can contain the characters `a-z`, `A-Z`, and `0-9`.

The edge length (here, `9909`) can contain the characters `0-9`.

The edge coverage (here, `6.94721`) can contain the characters `0-9` and `.`.

An edge name can optionally end with a `'` character, indicating that
this edge's reverse complement is being referenced.

We assume that each sequence (the line(s) between edge declarations)
consists only of the characters `A`, `C`, `G`, `T`, or `U`. Lowercase characters
or degenerate nucleotides are not allowed; this matches section 15 of version
1.00 of the [FASTG spec](http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf).
(The FASTG spec doesn't explicitly allow for uracil [`U`], but we do anyway to
allow for RNA sequences.)

Leading and trailing whitespace in sequence lines will be ignored, so something
like
```bash
    ATC

 G     
```
is technically valid (however, a line like `ATC G` is not valid since the inner
space, ` `, would be considered part of the sequence).

It is also worth noting that pyfastg **only creates nodes based on the edges
explicitly described in the graph**: if your graph only contains edges
`EDGE_1_...`, `EDGE_2_...`, and `EDGE_3_...`, then
pyfastg will only create nodes 1+, 2+, and 3+, and not the reverse complement
nodes 1-, 2-, 3-, etc.

Similarly, if your graph contains an adjacency from edge `EDGE_1_...` to
`EDGE_2_...'`, then this adjacency will only be represented as a single edge
(1+ → 2-) in pyfastg's output graph. The implied reverse-complement of this
edge (2+ → 1-) will not be automatically created.

### Identified node attributes
Nodes in the returned `DiGraph` (represented in the FASTG file as `EDGE_`s)
contain three attribute fields:

1. `length`: the length of the sequence (represented as a python `int`)
2. `cov`: the coverage of the sequence (represented as a python `float`)
2. `gc`: the GC-content of the sequence (represented as a python `float`)

Furthermore, every node's name will end in `-` if the node is a "reverse
complement" (i.e. if its declaration in the FASTG file ends in a `'` character) and `+` otherwise.

### Dependencies

- [NetworkX](https://networkx.github.io)

### License
pyfastg is licensed under the MIT License. Please see the `LICENSE` file for details.
