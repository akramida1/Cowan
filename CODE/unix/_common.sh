#!/usr/bin/env bash
#
# Shared by rcn, rcn2, rcg and rce.  This file is sourced, not run.
#
# It carries the two things all four wrappers do: print the program's
# banner, and time the run.  On Windows the second of those is a separate
# batch file, mark_time.bat, because cmd.exe has no functions worth the
# name and a `call` was the only way to share the arithmetic.  A shell has
# functions, so there is no reason for a separate mark_time here and none
# is provided.

# CODE/, the directory the executables and the .py tools live in.  The
# wrappers are one level below it, in CODE/unix/.
#
# `cd ... && pwd` and not `readlink -f`: readlink on macOS is BSD readlink
# and has no -f.  This gives an absolute, symlink-free directory either
# way, which is what RCG's COWAN.CFG needs.
cowan_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)

# ------------------------------------------------------------------ banner

# cowan_banner <box-indent> <line> ...
#
# The box around START is drawn in ASCII.  The batch files use the CP437
# line-drawing characters, which are right on a Windows console and are
# mojibake anywhere else -- a UTF-8 terminal shows them as accented
# letters.  The width and the indentation are the batch files' own, so the
# banner still lands in the same columns.
cowan_banner() {
    local indent=$1 pad line
    shift

    # cls, guarded.  Only when someone is actually watching: with output
    # redirected to a file the escape sequence would be written into it.
    if [ -t 1 ] && command -v clear >/dev/null 2>&1; then
        clear
    fi

    echo "."
    for line in "$@"; do
        echo "$line"
    done

    pad=$(printf '%*s' "$indent" '')
    echo "$pad+-------+"
    echo "$pad| START |"
    echo "$pad+-------+"
}

# ------------------------------------------------------------------ timing

# The batch files print the time to a hundredth of a second, which they
# can because cmd.exe's %time% carries it.  Portable shell does not:
# `date +%N` is a GNU extension and prints a literal N on macOS, and bash
# 5's $EPOCHREALTIME is not in the bash 3.2 that macOS ships.  So this
# counts whole seconds, using bash's own SECONDS, and the runs this timing
# is for are minutes to hours long.
cowan_stamp() {
    date '+%H:%M:%S, %a %m/%d/%Y'
}

# cowan_run <executable> [argument ...]
#
# Runs it, times it, and hands back its exit code -- which mark_time.bat
# does not do, having lost it to the arithmetic before the end.
cowan_run() {
    local exe=$1 rc elapsed hr min sec
    shift

    if [ ! -x "$exe" ]; then
        echo "$(basename "$0"): $(basename "$exe") is not in $cowan_dir" >&2
        echo "  Build it first:  $cowan_dir/unix/build.sh" >&2
        return 127
    fi

    echo "Started at                                   $(cowan_stamp)"
    echo "."

    SECONDS=0
    "$exe" "$@"
    rc=$?
    elapsed=$SECONDS

    echo "."
    hr=$((elapsed / 3600))
    min=$(((elapsed % 3600) / 60))
    sec=$((elapsed % 60))
    echo "Finished at                                  $(cowan_stamp)"
    echo "Time spent:                                  $hr hr $min min $sec sec"

    return $rc
}
