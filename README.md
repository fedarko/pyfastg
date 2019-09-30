# pyfastg: a minimal Python library for parsing networks from FASTG files

## The FASTG file format
FASTG is a format to describe genome assemblies, geared toward accuratly representing the ambiguity resulting from sequencing limitations, ploidy, or ther factors that complicate representation of a seqence as a simple string.  The official spec for the fastg format can be found [here](http://fastg.sourceforge.net/).

## pyfastg
This library contains `parse_fastg()`, a function that accepts as input a path
to a SPAdes FASTG file. This function parses the structure of the specified
file, returning a NetworkX `DiGraph` object representing the structure of the
graph.

This library is very much in its infancy, so it may be most useful as a starting point.  Pull requests welcome!

### Quick Example

```python
>>> import pyfastg
>>> g = pyfastg.parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
>>> # g is now a NetworkX DiGraph
>>> g.nodes()
NodeView(('1+', '29+', ...))
```

### Known Limitations
Currently, pyfastg is hardcoded to parse FASTG files created by the SPAdes assembler. Other valid FASTG files that don't follow the pattern used by SPAdes for contig names are not supported.

Due to the way that FASTG handles sequence gaps, ambiguities, etc, we do not attempt to actually parse the sequence that the nodes represent, merely the relationships between the nodes.

### Identified node attributes
Nodes in the returned `DiGraph` (represented in the FASTG file as `EDGE_`s)
contain two attribute fields:

1. `length`: the length of the node (represented as a python `int`)
2. `coverage`: the coverage of the noed (represented as a python `float`)

Furthermore, every node's name will end in `-` if the node is a "reverse
complement", and `+` otherwise.
