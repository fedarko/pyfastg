# pyfastg: a minimal Python library for parsing networks from FASTG files

# Introduction
FASTG is a format to describe genome assemblies, geared toward accuratly representing the ambiguity resulting from sequencing limitations, ploidy, or ther factors that complicate representation of a seqence as a simple string.  The official spec for the fastg format can be found [here](http://fastg.sourceforge.net/).

# Implementation
This library is very much in its infancy, so it may be most useful as a starting point.  Pull requests welcome!

## Known Limitations
Currently, this is hardcoded to parse FASTG files created by the SPAdes assembler. Other valid FASTG files that dont follow the pattern used by SPAdes for contig names are not supported.

Futher, the user is left to use their prefered method of manipulating networks once the data have been parsed.

Due to the way that FASTG handles sequence gaps, abiguities, etc, we do not attempt to actually parse the sequence that the nodes represent, merely the relationships between the nodes.

## Classes
In an effort to simplify the representation, there are two types of classes for nodes: the main *FastgNode* class, and the *FastgNodeNeighbor* class.

### FastgNode
This stores basic information about  the node: the name (usually a number), length of the sequence, a list of neighboring nodes, and whether the sequence is reverse-complimented.

### FastgNodeNeighbor
This stores limited information about the node as its relationship to another node: the name (usually a number) and whether the sequence is reverse-complimented.  This provides enough information to find the full record for the neighbor given a list of all nodes.
