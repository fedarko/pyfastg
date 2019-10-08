# pyfastg: a minimal Python library for parsing networks from SPAdes FASTG files
[![Build Status](https://travis-ci.org/fedarko/pyfastg.svg?branch=master)](https://travis-ci.org/fedarko/pyfastg) [![Code Coverage](https://codecov.io/gh/fedarko/pyfastg/branch/master/graph/badge.svg)](https://codecov.io/gh/fedarko/pyfastg)

## The FASTG file format
FASTG is a format to describe genome assemblies, geared toward accurately representing the ambiguity resulting from sequencing limitations, ploidy, or other factors that complicate representation of a seqence as a simple string.  The official spec for the FASTG format can be found [here](http://fastg.sourceforge.net/).

pyfastg parses graphs that follow **a subset of this specification**: in
particular, it is designed to work with files output by the
[SPAdes](http://cab.spbu.ru/software/spades/) family of assemblers.

## pyfastg
pyfastg contains `parse_fastg()`, a function that accepts as input a path
to a SPAdes FASTG file. This function parses the structure of the specified
file, returning a [NetworkX](https://networkx.github.io) `DiGraph` object representing
the structure of the graph.

pyfastg is very much in its infancy, so it may be most useful as a starting point.
Pull requests welcome!

### Quick Example

```python
>>> import pyfastg
>>> g = pyfastg.parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
>>> # g is now a NetworkX DiGraph! We can do whatever we want with this object.
>>> # Example: List the nodes in g
>>> g.nodes()
NodeView(('1+', '29-', '1-', '6-', '2+', '26+', '27+', '2-', '3+', '4+', '6+', '7+', '3-', '33-', '9-', '4-', '5+', '5-', '28+', '7-', '8+', '28-', '9+', '8-', '12-', '10+', '12+', '10-', '24-', '32-', '11+', '30-', '11-', '27-', '19-', '13+', '25+', '31-', '13-', '14+', '14-', '26-', '15+', '15-', '23-', '16+', '16-', '17+', '17-', '19+', '18+', '33+', '18-', '20+', '20-', '22+', '21+', '21-', '22-', '23+', '24+', '25-', '29+', '30+', '31+', '32+'))
>>> # Example: Get details for a single node (length, coverage, and GC-content)
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
Currently, pyfastg is hardcoded to parse FASTG files created by the SPAdes assembler. Other valid FASTG files that don't follow the pattern used by SPAdes for node names are not supported.

In particular, each node in the file must be declared as

```bash
>EDGE_1_length_9909_cov_6.94721
```

The node ID (here, `1`) can contain the characters `a-z`, `A-Z`, and `0-9`.

The node length (here, `9909`) can contain the characters `0-9`.

The node coverage (here, `6.94721`) can contain the characters `0-9` and `.`.

We assume that each node sequence (the line(s) between node declarations)
consists only of valid DNA characters, as determined by
[`skbio.DNA`](http://scikit-bio.org/docs/latest/generated/skbio.sequence.DNA.html).
Leading and trailing whitespace in sequence lines will be ignored, so something
like
```bash
    ATC

 G     
```
is perfectly valid (however, `ATC G` is not since the inner space, ` `, will be
considered part of the sequence).

It is also worth noting that pyfastg **only creates nodes/edges based on those
observed in the graph**: if your graph only contains nodes 1+, 2+, and 3+, then
this won't automatically create reverse complement nodes 1-, 2-, 3-, etc.

### Identified node attributes
Nodes in the returned `DiGraph` (represented in the FASTG file as `EDGE_`s)
contain three attribute fields:

1. `length`: the length of the node (represented as a python `int`)
2. `cov`: the coverage of the node (represented as a python `float`)
2. `gc`: the GC-content of the node's sequence (represented as a python `float`)

Furthermore, every node's name will end in `-` if the node is a "reverse
complement" (i.e. if its declaration in the FASTG file ends in a `'` character) and `+` otherwise.

### Dependencies

- [NetworkX](https://networkx.github.io)
- [scikit-bio](http://scikit-bio.org/)
