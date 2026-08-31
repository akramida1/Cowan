# NOTICE

This file records where the parts of this package came from. The terms under
which each part may be used are in `LICENSE.md`.

## Code provenance

This software is a substantially revised, extended, and modernized descendant of
the atomic-structure programs originally developed by Robert D. Cowan at Los
Alamos. The four main programs — RCN (`FOR/RCN36K.F`), RCN2 (`FOR/RCN2K.F`),
RCG (`FOR/rcg11kd.f`) and RCE (`FOR/RCE20KD.F`) — descend directly from his, and
each retains in its header the 1980 Los Alamos / U.S. Government notice under
which Cowan released it. **That notice permits free distribution anywhere in the
world but forbids sale, and it still governs these files.** It is reproduced in
section 1 of `LICENSE.md`.

The present implementation retains parts of the historical Cowan code and much
of the scientific and computational lineage of those programs, but has undergone
extensive redevelopment. Among other changes, static storage and extensive use of
`EQUIVALENCE` have been replaced by dynamic memory management, parallel
processing has been introduced, true *jj* coupling has been added to RCG and RCE,
and numerous algorithms and data structures have been modified for efficiency and
for use on modern computer systems.

That redevelopment happened in two stages, and the distinction matters legally as
well as historically:

- **Through the 2021 release**, the work was done by Alexander Kramida as an
  employee of the National Institute of Standards and Technology, and was
  distributed by NIST from <https://doi.org/10.18434/T4/1502500>. It includes the
  packed binary file format (`PK`/`UNPK`), eigenvector recognition in RCE, the
  Perl utilities now retired to `CODE/perl_legacy/`, and a large number of
  corrections to the original programs. Software developed by NIST employees is
  not subject to copyright protection within the United States; the NIST software
  notice under which those versions were distributed is reproduced in section 2
  of `LICENSE.md`.
- **From 2022 onward**, the work has been done by Alexander Kramida at the
  Department of Astronomy, University of Maryland. It includes the Cray-pointer
  restructuring of memory use, true *jj* coupling, the speed and memory work,
  OpenMP parallelization, the per-configuration term cutoff, the Python
  utilities and their test suite, the port to Linux and macOS, and the
  fundamental-constants work. These portions are offered under the MIT License;
  see section 3 of `LICENSE.md`.

## Other contributors

The `INWOL` and `IQS` options of RCG (columns 77-78 of a configuration card in
`ING11`), and the associated code in `FOR/rcg11kd.f`, were contributed by
Yu. Ralchenko.

## GRASP memory-management routines

Three memory-management routines used by this package,

- `ALLOC`
- `RALLOC`
- `DALLOC`

originated in development versions of the GRASP atomic-structure software and
were provided to Alexander Kramida during collaborative GRASP development at
NIST.

These routines occur in two components of the present package and therefore
appear twice in the source tree, in `FOR/rcg11kd.f` and in `FOR/RCE20KD.F`.

The original source comments identifying their authors and development history
have intentionally been retained. In particular, the routines acknowledge
contributions by Charlotte Froese Fischer, Farid A. Parpia, Jacek Bieroń, and
other GRASP developers as recorded in the individual source files.

These routines predate GRASP2018. They employ the older GRASP2K-style
dynamic-memory mechanism based on Cray pointers and therefore are not copies of
the later Fortran 95 memory-management implementation distributed with GRASP2018.

GRASP2018 and its successors are distributed by the GRASP/COMPAS collaboration
under the MIT License. The inclusion of the three historical routines in this
package should not, however, be interpreted as a claim that their earlier source
text became MIT-licensed solely by virtue of the later GRASP2018 release.

The routines are retained with their original authorship acknowledgments in
recognition of the contributions of the GRASP developers to this software.

## Cowan's documentation

`DOCS/` contains verbatim, unaltered copies of R. D. Cowan's original
documentation for the programs, preserved for historical purposes. It is
redistributed as received and is not modified by this package.

## Bundled runtime libraries

The Windows distribution ships six DLLs from the GCC and MinGW-w64 projects
beside the executables, because OpenMP cannot be statically linked with this
toolchain. They are redistributed unmodified, as built by WinLibs; their
licences are unchanged and are summarised in section 4 of `LICENSE.md`. On Linux
and macOS the corresponding libraries come from the system and nothing is
bundled.

## Attribution

Users of this software are requested to acknowledge its Cowan, NIST, and GRASP
heritage as appropriate in scientific publications.

Once a publication describing the present package is available, users will also
be requested to cite that publication. Citation information will be added to the
repository when available.
