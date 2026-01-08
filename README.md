# pyfastg: a minimal Python library for parsing FASTG files

<!-- Regarding centering images, see https://gist.github.com/DavidWells/7d2e0e1bc78f4ac59a123ddf8b74932d -->
<div align="center">
<a href="https://github.com/fedarko/pyfastg/actions/workflows/main.yml"><img src="https://github.com/fedarko/pyfastg/actions/workflows/main.yml/badge.svg" alt="pyfastg CI" /></a>
<a href="https://codecov.io/gh/fedarko/pyfastg"><img src="https://codecov.io/gh/fedarko/pyfastg/branch/master/graph/badge.svg" alt="Code Coverage" /></a>
<a href="https://pypi.org/project/pyfastg"><img src="https://img.shields.io/pypi/v/pyfastg?color=006dad" alt="PyPI" /></a>
<a href="https://anaconda.org/bioconda/pyfastg"><img src="https://img.shields.io/conda/vn/bioconda/pyfastg.svg?color=43b02a" alt="bioconda" /></a>
</div>

## The FASTG file format
FASTG is a format for describing sequencing assembly graphs. It attempts to
accurately represent the ambiguity resulting from sequencing limitations, ploidy,
or other factors that complicate representation of a seqence as a simple string.

The latest specification for the FASTG format is version 1.00, as of writing;
the original FASTG website is down, but
an archived version of the v1.00 specification is
[accessible here](https://web.archive.org/web/20211209213905/http://fastg.sourceforge.net/FASTG_Spec_v1.00.pdf).
Whenever the rest of this documentation mentions "the FASTG spec," this is in reference
to this version of the specification.

pyfastg is a Python library designed to parse graphs that follow
**a subset of the FASTG spec**. In particular, pyfastg is designed to
work with files output by the [SPAdes](https://ablab.github.io/spades/)
family of assemblers. It also now supports files output by
[MEGAHIT](https://github.com/voutcn/megahit)!

## The pyfastg library
The pyfastg library contains `parse_fastg()`, a function that
takes as input a path to a FASTG file. `parse_fastg()` reads this
FASTG file and returns a [NetworkX](https://networkx.org/)
`DiGraph` object representing the structure of the assembly graph.

pyfastg is useful as a starting point for other applications.
Using this NetworkX `DiGraph` object, we can do whatever we want with the
assembly graph: analyze it, convert it to other formats, visualize it, etc.

### Note about the graph topology

The FASTG spec contains the following sentence (in section 6, page 7):

> Note also that strictly speaking, [the structure described in a FASTG file] is not a graph at all, as we have not specified a notion of vertex. However in many cases one can without ambiguity define vertices and thereby associate a _bona fide_ digraph, and we do so frequently in this document to illustrate concepts.

We use the following approach to get around this problem: **"edges" in the FASTG file will be represented as nodes in the NetworkX graph produced by pyfastg, and "adjacencies" between edges in the FASTG file will be represented as edges in the NetworkX graph produced by pyfastg.**

As far as we're aware, this "conversion" from edges to nodes matches
how FASTG files have often been visualized in the past.

### Installation
pyfastg can be installed using [pip](https://pip.pypa.io/) or [conda](https://conda.io/):

#### Installation using pip

```bash
pip install pyfastg
```

#### Installation using conda

```bash
conda install -c bioconda pyfastg
```

#### Dependencies
As of writing, pyfastg's only direct dependency (which should be installed
automatically when running either of the above installation commands) is
[NetworkX](https://networkx.org/). pyfastg requires a minimum NetworkX
version of 2.

As of writing, pyfastg supports Python 3.6 and up.
pyfastg might be able to work with earlier versions of Python,
but we do not explicitly test against these.

### Quick example: using pyfastg to load and analyze an assembly graph
The second line (which points to one of pyfastg's test assembly graphs)
assumes that you're located in the root directory of the pyfastg repo.

```python
>>> import pyfastg
>>> g = pyfastg.parse_fastg("pyfastg/tests/input/assembly_graph.fastg")
>>> # g is now a NetworkX DiGraph! We can do whatever we want with this object.
>>>
>>> # Example: List the sequences in this graph (these are "edges" in the FASTG
>>> # file, but are represented as nodes in g)
>>> g.nodes()
NodeView(('1+', '29-', '1-', '6-', '2+', '26+', '27+', '2-', '3+', '4+', '6+', '7+', '3-', '33-', '9-', '4-', '5+', '5-', '28+', '7-', '8+', '28-', '9+', '8-', '12-', '10+', '12+', '10-', '24-', '32-', '11+', '30-', '11-', '27-', '19-', '13+', '25+', '31-', '13-', '14+', '14-', '26-', '15+', '15-', '23-', '16+', '16-', '17+', '17-', '19+', '18+', '33+', '18-', '20+', '20-', '22+', '21+', '21-', '22-', '23+', '24+', '25-', '29+', '30+', '31+', '32+'))
>>>
>>> # Example: Get details for a single sequence (length, coverage, GC-content)
>>> g.nodes["15+"]
{'length': 193, 'cov': 6.93966, 'gc': 0.5492227979274611}
>>>
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

### Details about the required input file format (tl;dr: SPAdes-dialect FASTG files only)
Currently, pyfastg is only designed to parse FASTG files created by the SPAdes or MEGAHIT
assemblers. Other valid FASTG files that don't follow the formats used by SPAdes or
MEGAHIT are not supported. (If you would like us to add support for a new assembler's
output, please open an issue!)

#### Edge names

Each sequence in the file should have a name formatted like:

<table align="center">
 <thead>
  <tr>
   <th>SPAdes</td>
   <th>MEGAHIT</td>
  </tr>
 </thead>
 <tbody>
  <tr>
   <td><code>EDGE_1_length_9909_cov_6.94721</code></td>
   <td><code>NODE_1_length_9909_cov_6.9472_ID_1</code></td>
  </tr>
 </tbody>
</table>

In MEGAHIT FASTG files, these sequences are referred to as `NODE`s instead of
`EDGE`s. We will keep saying "edge" throughout this documentation for the sake of
simplicity.

The edge ID (here, `1`) can contain the characters `a-z`, `A-Z`, and `0-9`.
In MEGAHIT files, there are two IDs included in each name -- one after the
`NODE_` and one at the very end of the name (`ID_`). We will only use the
first one.

The edge length (here, `9909`) can contain the characters `0-9`.

The edge coverage (here, `6.94721`) can contain the characters `0-9` and `.`.

An edge name can optionally end with a `'` character, indicating that
this edge is a reverse complement. We will refer to whether or not an edge name
ends with `'` as its _orientation_: an edge that does not end with a `'` has a
`+` orientation, and an edge name that ends with a `'` has a `-` orientation.

All edge names in a FASTG file should be consistent with respect to a given
ID and orientation.
If, in a single FASTG file, pyfastg sees a reference to an edge named
`EDGE_1_length_9909_cov_6.94721` and also a reference to an edge named
`EDGE_1_length_9908_cov_6.95` (with the same ID [`1`]
and orientation [`+`], but a different length and/or coverage)
then it will throw an error.

#### Edge declaration lines

Here, we refer to each line starting with `>` as an _edge declaration_. An
edge's sequence is described in the line(s) following its edge declaration
(until the next edge declaration); additionally, the outgoing adjacencies from
this edge to other edges may be described on this line, if present. For example,
the line

```
>EDGE_1_length_5_cov_10:EDGE_2_length_3_cov_1,EDGE_3_length_6_cov_2.5',EDGE_4_length_8_cov_5.1;
```

indicates that the edge `EDGE_1_length_5_cov_10` has three outgoing adjacencies: to the
edges `EDGE_2_length_3_cov_1`, `EDGE_3_length_6_cov_2.5'`, and `EDGE_4_length_8_cov_5.1`.
This line would thus result in three "edges" being created
in the NetworkX graph produced by pyfastg: (`1+` → `2+`), (`1+` → `3-`), and (`1+` → `4+`).

Each edge declaration must end with a `;` character (after removing trailing
whitespace). Section 15 of the FASTG spec mentions that having a newline
after the semicolon isn't required, but we require it here for the sake of
simplicity.

#### Edge sequences

We assume that each sequence (the line(s) between edge declarations)
consists only of the characters `A`, `C`, `G`, `T`, or `U`. So, more complex
types of strings (e.g. the "stuffed gaps" described in the FASTG spec) are
not allowed in an edge's sequence.

Additionally, lowercase characters or degenerate nucleotides are not allowed;
this matches section 15 of the FASTG spec.
The FASTG spec doesn't explicitly allow for uracil (`U`), but we allow it
anyway in order to support RNA sequences. (`U` and `T` are allowed to be contained
in the same sequence,
[in the unlikely case that this is needed](https://en.wikipedia.org/wiki/Uracil#In_DNA).)

Leading and trailing whitespace in sequence lines will be ignored, as will
blank lines within a sequence. So, something like

```
>EDGE_1_length_4_cov_100;
    ATC

 G     
```

is technically valid: this sequence is read as `ATCG`.
However, the following example:

```
>EDGE_1_length_4_cov_100;
ATC G
```

is not valid and will cause pyfastg to throw an error.
This is because the inner space between
the `C` and the `G` would be read as part of the sequence.

### Details about the output NetworkX graph

#### Node names and attributes
Nodes in the returned `DiGraph` (corresponding to edges in the FASTG file)
will contain three attribute fields:

1. `length`: the length of the sequence (represented as a python `int`)
2. `cov`: the coverage of the sequence (represented as a python `float`)
2. `gc`: the GC-content (in the range [0, 1]) of the sequence (represented as a python `float`)

Each node's name is a python `str` created by concatenating edge IDs and orientations.
For example, `EDGE_1_length_9909_cov_6.94721` will correspond to a node named `1+`.
This naming scheme is analogous to that used by
[Bandage](https://github.com/rrwick/Bandage/wiki/Single-vs-double-node-style).

#### About reverse complements

pyfastg **only creates nodes based on the edges
explicitly described in the FASTG file**. If a file only describes edges
`EDGE_1_length_5_cov_10`, `EDGE_2_length_6_cov_10'`, and `EDGE_3_length_7_cov_15`, then
pyfastg will only create nodes `1+`, `2-`, and `3+`, and not the reverse complement
nodes `1-`, `2+`, `3-`, etc.

Similarly, if a file contains an adjacency from edge `EDGE_1_length_5_cov_10` to
`EDGE_2_length_6_cov_10'`, then this adjacency will only be represented as a single edge
(`1+` → `2-`) in pyfastg's output graph. The implied reverse-complement of this
edge (`2+` → `1-`) will not be created unless the file explicitly
contains an adjacency from `EDGE_2_length_6_cov_10` to `EDGE_1_length_5_cov_10'`.

## Information for pyfastg developers

Pull requests are welcome! If you're interested in developing pyfastg's code,
this section provides some instructions for getting started.

### Setting up a development "environment" for pyfastg

You will probably want to fork this repository and then clone your fork to your
computer. Once you do this, `cd` into the root of the repository and run

```bash
pip install -e .[dev]
```

to install pyfastg in "editable mode." Thanks to the `[dev]` flag, this will also install
pyfastg's development dependencies (see the `extras_require` line in
[`setup.py`](https://github.com/fedarko/pyfastg/blob/master/setup.py) for details).

### Testing, linting, and formatting the code

pyfastg's [`Makefile`](https://github.com/fedarko/pyfastg/blob/master/Makefile)
contains targets that perform these three tasks:

- Run tests: `make test`
- Lint and style-check the code: `make stylecheck`
- Automatically style the code: `make style`

These targets should all be run from the root of the pyfastg repository. They
should hopefully be self-explanatory, but let us know if you have
any questions.

## Test data sources

The test data file located in `pyfastg/tests/input/megahit-example.fastg` was
computed by running `megahit_core contig2fastg` on `k21.contigs.fa`, a FASTA file from
[this GitHub repository](https://github.com/jraysajulga/megahit-contig2fastg-wrapper).

## Changelog
See pyfastg's
[`CHANGELOG.md`](https://github.com/fedarko/pyfastg/blob/master/CHANGELOG.md) file
for information on the changes included with new pyfastg releases.

## License
pyfastg is licensed under the MIT License. Please see pyfastg's
[`LICENSE`](https://github.com/fedarko/pyfastg/blob/master/LICENSE) file for details.

## Contact
The recommended way to get in touch with pyfastg's developers is by
[opening a GitHub issue](https://github.com/fedarko/pyfastg/issues).
