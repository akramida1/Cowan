# NOTICE

## Code provenance

This software is a substantially revised, extended, and modernized descendant of the atomic-structure programs originally developed by Robert D. Cowan at Los Alamos.

The present implementation retains parts of the historical Cowan code and much of the scientific and computational lineage of those programs, but has undergone extensive redevelopment. Among other changes, static storage and extensive use of `EQUIVALENCE` have been replaced by dynamic memory management, parallel processing has been introduced, and numerous algorithms and data structures have been modified for efficiency and for use on modern computer systems.

Intermediate versions of this software were developed and distributed at the National Institute of Standards and Technology (NIST). Those versions were distributed under the NIST software notice reproduced in the `LICENSE` file.

## GRASP memory-management routines

Three memory-management routines used by this package,

- `ALLOC`
- `RALLOC`
- `DALLOC`

originated in development versions of the GRASP atomic-structure software and were provided to Alexander Kramida during collaborative GRASP development at NIST.

These routines occur in two components of the present package and therefore appear twice in the source tree.

The original source comments identifying their authors and development history have intentionally been retained. In particular, the routines acknowledge contributions by Charlotte Froese Fischer, Farid A. Parpia, Jacek Bieroń, and other GRASP developers as recorded in the individual source files.

These routines predate GRASP2018. They employ the older GRASP2K-style dynamic-memory mechanism based on Cray pointers and therefore are not copies of the later Fortran 95 memory-management implementation distributed with GRASP2018.

GRASP2018 and its successors are distributed by the GRASP/COMPAS collaboration under the MIT License. The inclusion of the three historical routines in this package should not, however, be interpreted as a claim that their earlier source text became MIT-licensed solely by virtue of the later GRASP2018 release.

The routines are retained with their original authorship acknowledgments in recognition of the contributions of the GRASP developers to this software.

## Attribution

Users of this software are requested to acknowledge its Cowan, NIST, and GRASP heritage as appropriate in scientific publications.

Once a publication describing the present package is available, users will also be requested to cite that publication. Citation information will be added to the repository when available.