# Bioseqlib #

_Victor E. Bazterra_

_Center for High Performance Computing_

_University of Utah (2005)_

# Description #

Bioseqlib is small library in C++ for bioinformatics applications developted at Center for High Performance Computing University of Utah. Originally, it was built as a part of a bigger project providing basics algorithms like memory management, sequence manipulation and dynamic programming. Since then, the library has grown little more including binary tree operations and some algorithms for calculating phylogenetic trees. I/O functions are done by seqio library, a C library designed to read biological sequences from different databases. A series of wrapper functions are used to access seqio capabilities from bioseqlib.

The structure of this library will allow to be a natural container for more complex implementations. It is organized to improve code portability and re-usability. Another feature is that it is well documented using doxygen as documentation generator.


# Bioseqlib features #

  * Basic containers to represent biological sequences and any of their possible alignments
  * Amino acid alphabet including the PAM and BLOSUM score matrices.
  * A general binary tree definition and operation.
  * Algorithms to create a rooted or unrooted trees given the distances between vertices.
  * C++ adaptation of seqio library.
  * Algorithms to calculate the sequence metrics and/or distances.
  * Linear and quadratic space algorithms for optimal global alignment between two sequences using affinitive gaps.

# Installation #

It is necessary to define BIOSEQLIBROOT that should be ponting to the root directory where the library is install.

After that just type 'make' in the BIOSEQLIBROOT directory.

Example:

```
~> mkdir bioseqlib
~> cp bioseqlib-x.x.tar.gz bioseqlib/
~> cd bioseqlib
bioseqlib> export BIOSEQLIBROOT=$PWD (for bash)
```