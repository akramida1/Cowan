#!/usr/bin/env bash
#
# The debug build: gfortran_dbg.bat's counterpart.
#
#   ./build_debug.sh [program ...]
#
# -O0 -ggdb3 -fbacktrace in place of -O3, and serial, because a debugger
# is far easier to follow one thread at a time and the codes give the
# same numbers threaded or not.
#
# -fno-automatic is kept even here.  It puts locals in static storage;
# drop it and the numbers change, and a debug build would no longer be
# debugging the same code.
#
# The debugger is gdb on Linux.  On macOS gdb is awkward to codesign, and
# lldb -- which ships with the Xcode command line tools -- reads gfortran's
# DWARF well enough for this package.

exec "$(cd "$(dirname "$0")" && pwd)/build.sh" --debug "$@"
