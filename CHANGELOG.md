# pyfastg changelog

## pyfastg v0.1.0 (June 26, 2022)

It's been a while! This release includes a fair amount of changes; the bulk of
these make the parser more strict when dealing with unsupported or malformed
FASTG files.

### New features
- Sequences containing `U` (uracil) are now explicitly allowed. Previously,
  these would have triggered an error from scikit-bio about them not being
  valid DNA characters.

### Documentation
- Improved some of the error messages produced when pyfastg encounters
  malformed graphs.

- Added information to the README (and various error messages) to be
  clearer about the fact that the sequences in the input FASTG file are "edges,"
  but are represented as nodes in the output NetworkX graph.

- Updated the README to clarify pyfastg's behavior when an edge's implied
  reverse-complement (e.g. for 1+ → 2-, this would be 2+ → 1-) does not
  explicitly exist in the FASTG file. (The implied edge is not added if it
  isn't described in the FASTG file.)

- Restructured and improved the README in various small ways (e.g. moved
  the information on pyfastg's dependencies to be alongside the installation
  instructions).
  
### Backward-incompatible changes

N/A. Some of the ways in which we've made validation more strict might be
technically "backward-incompatible," but these sorts of situations were already
unsupported by pyfastg due to not being part of the general SPAdes FASTG dialect.

### Bug fixes
- The error thrown upon seeing an invalid declaration was meant to show the
  problematic declaration, but it just showed the string `{}` (due to Python's
  `.format()` function not being used). This has been fixed -- this error
  message now shows the declaration that caused a problem.

- The `|` character was previously
  [included](https://github.com/fedarko/pyfastg/blob/3b99ba8624339e3efdc87de1d122e5f401b4537a/pyfastg/pyfastg.py#L43)
  in the set of accepted characters our regular expression used when extracting
   an edge's coverage (e.g. the coverage `12|3.5` would have technically been
   accepted by the regular expression).

    - This shouldn't have caused any major problems, since [a few lines later](https://github.com/fedarko/pyfastg/blob/3b99ba8624339e3efdc87de1d122e5f401b4537a/pyfastg/pyfastg.py#L60)
      Python would have thrown an error about such a string not being convertible
      to a float. However, this is now caught further up in this function: the
      regular expression no longer accepts `|`.

- Added checks for rare problems with inconsistent edge declarations.

    - Previously, edge declarations were parsed using a [regular expression](https://github.com/fedarko/pyfastg/blob/3b99ba8624339e3efdc87de1d122e5f401b4537a/pyfastg/pyfastg.py#L43)
  that did not use the `^` or `$` characters (see
  [Python's documentation on regular expressions](https://docs.python.org/3/library/re.html)
  for details).

    - Long story short, this meant that an edge formatted like
      `asdfEDGE_1_length_9_cov_3ghij` would have technically been accepted, even
      though it contained extra information before (`asdf`) and after (`ghij`) the
      relevant parts of the edge name.

    - This could have caused problems in rare cases (e.g. if a graph contained
      both an edge labelled `asdfEDGE_1_length_9_cov_3ghij` and an edge named
      `EDGE_1_length_9_cov_3`, these edges would have been treated as if they
      were the same edge since the ID `1` is the same).

    - To prevent this sort of problem, we now enforce that edge
      names cannot have any extra leading or trailing information.

    - Additionally, we now check each occurrence of an edge's declaration in
      the FASTG file (both on the line that this edge is declared, and on any
      lines where other edges have outgoing adjacencies to this edge). If any
      declarations are inconsistent (e.g. we see `EDGE_1_length_20_cov_5.2` on
      one line, but `EDGE_1_length_20_cov_4.5` on another line), we will raise
      an error about this inconsistency.

- We now detect and throw an error if an edge has multiple outgoing adjacencies
  to the same edge (e.g. we see a line formatted like `>x:y,z,y`).

    - Previously, nothing would have happened in this case (only one copy of
      the edge from `x`'s node to `y`'s node would be created), since the type
      of graph we create (`networkx.DiGraph`) does not support multi-edges.
      Throwing an error in this case seems like the better way to handle this.

- Detect and throw an error if an edge declaration line does not end with `;`.

- Detect and throw clear errors if we see syntax associated with certain
  unsupported things in the FASTG spec (e.g. the `~` or `[]` notations).

### Performance improvements

N/A. The extra layers of validation added in this version might slow down
pyfastg slightly, but (in our opinion) this is worth it.

### Development improvements
- Switched from Travis CI to GitHub Actions ([#3](https://github.com/fedarko/pyfastg/issues/3)).

- Set up the GitHub Actions CI to test against Python 3.6, 3.7, 3.8, 3.9, 3.10,
  and the latest development version of 3.11. Previously, the CI only checked
  against Python 3.6 and 3.7.

- Added development instructions to the README.

### Miscellaneous
- Removed dependency on scikit-bio for validating sequences and computing GC
  content ([#2](https://github.com/fedarko/pyfastg/issues/2)). Now, the only
  non-development dependency of pyfastg is NetworkX. (That being said, we may
  add more dependencies in the future if absolutely necessary.)

## pyfastg v0.0.0 (October 14, 2019)

Initial release.

Note that pyfastg v0.0.0 was released on PyPI on October 7, 2019; there are a
few small changes to pyfastg's documentation that were committed to the GitHub
repository from October 7, 2019 to October 14, 2019. These small changes are
thus not reflected in the PyPI version of the README for v0.0.0.
