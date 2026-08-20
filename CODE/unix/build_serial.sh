#!/usr/bin/env bash
#
# The non-OpenMP build: gfortran_O3.bat's counterpart.
#
#   ./build_serial.sh [program ...]
#
# Everything the threaded build produces, minus -fopenmp, so the
# executables need no libgomp at run time.  The numbers are the same
# either way -- no parallel region in this package reassociates a sum --
# so this costs only speed.
#
# It is one line because build.sh does the work; the separate name is
# there because the batch files it replaces are separate.

exec "$(cd "$(dirname "$0")" && pwd)/build.sh" --serial "$@"
