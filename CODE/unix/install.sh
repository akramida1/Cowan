#!/usr/bin/env bash
#
# Install the Cowan package on Linux or macOS.
#
#   ./install.sh                 into ~/.local/lib/cowan, or, run as root,
#                                into /usr/local/lib/cowan
#   ./install.sh --user          force ~/.local/lib/cowan
#   ./install.sh --system        force /usr/local/lib/cowan
#   ./install.sh --prefix DIR    install into DIR
#   ./install.sh --dry-run       print what would be copied, change nothing
#   ./install.sh --rebuild-cfp   rebuild the data decks even if they are there
#
# There is no Windows counterpart: on Windows the author copies CODE\ to
# C:\cowan\code by hand and puts that on PATH.  This script exists because
# Unix can enforce something Windows only lets us promise -- see below.
#
# WHAT IS INSTALLED, AND WHAT IS NOT
#
# The installed tree is CODE/ minus everything that is Windows-only or
# development-only:
#
#   installed    the four executables; CODE/unix/ (the wrappers); the .py
#                tools, toolbox/ and lsjj/; the data files RCG and RCN2
#                read -- FOR072, FOR073, FOR074, SENIOR, LSJJ.DAT --
#                plus ING11.CFP and calc_TP_M1E2.par.  The three FOR07x
#                are the CFP decks and LSJJ.DAT is the LS-to-jj deck;
#                this script builds all four first, see below.
#
#   not          *.bat and mark_time.bat (cmd.exe), *.exe (Windows
#                builds), *.dll (the MinGW runtime; on Unix libgfortran
#                and libgomp come from the system), tests/ and pytest.ini,
#                perl_legacy/, __pycache__/ and .pytest_cache/
#
# THE PERMISSIONS MODEL, WHICH IS THE POINT OF THIS SCRIPT
#
# The package's functional division is that the code directory is never
# written to during a run: every program runs in the working directory and
# writes only there.  On Windows that is a convention held up by the batch
# files.  Here it is made real -- directories 755, data files 644,
# executables and wrappers 755 -- so that a stray OPEN in the code
# directory fails loudly instead of quietly corrupting the package.
#
# It survives contact with the one program that would breach it.  RCG
# opens FOR072/073/074 with STATUS='UNKNOWN', which is a read-write open,
# and libgfortran retries a refused read-write open read-only; so an
# ordinary RCG run reads 644 decks without complaint.  The exception is a
# CFP-generation run, which really does write them -- hence the ordering
# below, and the permission dance in make_cfp.
#
# THE DERIVED DATA DECKS ARE BUILT HERE
#
# FOR072, FOR073, FOR074 and LSJJ.DAT are not part of the distribution.
# The three FOR07x are binary files RCG writes itself, from
# CODE/ING11.CFP, and RCG will not run without them.  LSJJ.DAT is the
# LS-to-jj transformation deck, written by CODE/lsjj/mkdeck.py from the
# same input; RCG opens it only for a run that asks for true jj coupling.
# So this script builds them: if they are not in the source tree it runs
# make_cfp before the copying -- here, from the unpacked tree, which is
# the last moment CODE/ is certain to be writable and the only one at
# which no permission dance is needed.  It takes well under a minute and
# asks nothing.  --rebuild-cfp forces it even when the decks are already
# there.
#
# LSJJ.DAT is the one file that may still be absent afterwards.  Its
# generator needs numpy -- the coefficients are null spaces and
# degenerate eigenspaces of the coupling operator, and another linear
# algebra would return a different, equally valid basis, which would no
# longer be the deck that was validated against Cowan's own coefficients.
# numpy is therefore the one thing in this package that has to be
# pip-installed, and it is optional: make_cfp says so and carries on, the
# install completes, and everything except true jj coupling works.
#
# make_cfp is given a scratch directory to work in.  It runs in the
# current directory -- ING11 and OUTG11 are written there, and an existing
# ING11 is backed up -- and the directory the user happened to invoke
# install.sh from is not ours to write in.
#
# Running make_cfp by hand is still meaningful afterwards, for a rebuild
# against an RCG compiled with a compiler whose unformatted record layout
# differs from gfortran's.  Against an installed tree that means either
# owning the tree or using sudo, and make_cfp says so.
#
# See build.sh for compiling, and README.md for what to do next.

set -u

self=$(basename "$0")

die() {
    echo "$self: $*" >&2
    exit 1
}

here=$(cd "$(dirname "$0")" && pwd)          # CODE/unix
code_dir=$(cd "$here/.." && pwd)             # CODE

# ---------------------------------------------------------------- options

prefix=""
dry=0
rebuild_cfp=0

while [ $# -gt 0 ]; do
    case "$1" in
        --user)    prefix=$HOME/.local/lib/cowan ;;
        --system)  prefix=/usr/local/lib/cowan ;;
        --prefix)  [ $# -ge 2 ] || die "--prefix needs a directory"
                   prefix=$2; shift ;;
        --prefix=*) prefix=${1#--prefix=} ;;
        --dry-run|-n) dry=1 ;;
        --rebuild-cfp) rebuild_cfp=1 ;;
        -h|--help) sed -n '3,11p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *)         die "unknown option $1" ;;
    esac
    shift
done

# The default follows who is running.  Installing to /usr/local as an
# ordinary user would fail at the first mkdir; installing to ~ as root
# would leave a root-owned tree in the invoking user's home.
if [ -z "$prefix" ]; then
    if [ "$(id -u)" = 0 ]; then
        prefix=/usr/local/lib/cowan
    else
        prefix=$HOME/.local/lib/cowan
    fi
fi

# Checked here rather than left to fail later as a missing FOR072.  RCG is
# told where the package is through COWAN.CFG, reads that line with
# FORMAT(A) and cuts it at the first blank, so a directory with a space in
# its name cannot work however carefully the wrapper quotes it.  This is
# the same test the rcg wrapper makes at run time; making it here means
# the user hears about it once, at install, and not on every run.
case "$prefix" in
    *\ *) echo "$self: the install directory must not contain a space:" >&2
          echo "  $prefix" >&2
          echo "  RCG reads that path from COWAN.CFG and cuts it at the first blank." >&2
          exit 1 ;;
esac

case "$prefix" in
    /*) ;;
    *)  prefix=$(pwd)/$prefix ;;
esac

# --------------------------------------------------------- what to copy

# Three lists, one per mode.  Kept explicit rather than derived from a
# find(1) with exclusions: a new Windows-only file dropped into CODE/
# should not silently start being installed.

executables="RCN36K RCN2K rcg11kd RCE20KD"

data="SENIOR ING11.CFP calc_TP_M1E2.par"

# Kept apart from $data because they are the one thing in the source tree
# that may legitimately be absent: they are built below rather than
# distributed.  They are copied with the rest, at the same mode.
#
# LSJJ.DAT is one of them because make_cfp builds it too, and it is the
# only entry in either list that may still be absent after that run: its
# generator needs numpy, which is optional.  The copying loop lets that
# one through; see there.
cfp_decks="FOR072 FOR073 FOR074 LSJJ.DAT"

# The wrappers, including build*.sh.  Installing the build scripts is
# deliberate: cowan_run's "Build it first: .../unix/build.sh" names them,
# and a build.sh run from an installed tree fails with its own clear
# message about COWAN_SRC rather than with a missing file.
wrappers="_common.sh build.sh build_debug.sh build_serial.sh
          rcn rcn2 rcg rce
          add_eav calc_TP_M1E2 conf_list config_purity conv_out
          insert_level_ids make_cfp make_levels_list_prepub
          make_levels_list_pub printout
          read_aa reorder_ing11 s11 scale_param subst_exp_levs
          subst_params_ing11 update11 install.sh"

tools=""
for f in "$code_dir"/*.py; do
    [ -e "$f" ] || continue
    tools="$tools $(basename "$f")"
done

toolbox=""
for f in "$code_dir"/toolbox/*.py; do
    [ -e "$f" ] || continue
    toolbox="$toolbox $(basename "$f")"
done

# The LS-to-jj generator.  Installed for make_cfp's sake and not for the
# user's: nothing in unix/ runs it directly, but make_cfp reaches for
# CODE/lsjj/mkdeck.py beside itself, in the installed tree as in this one.
lsjj=""
for f in "$code_dir"/lsjj/*.py; do
    [ -e "$f" ] || continue
    lsjj="$lsjj $(basename "$f")"
done

# ------------------------------------------------------------- preflight

missing=""
for f in $executables; do
    [ -f "$code_dir/$f" ] || missing="$missing $f"
done
if [ -n "$missing" ]; then
    echo "$self: not built yet -- missing in $code_dir:$missing" >&2
    echo "  Compile first:  $here/build.sh" >&2
    echo "  (The .exe files beside them are the Windows builds and are not" >&2
    echo "  installed; the Unix executables carry no suffix.)" >&2
    exit 1
fi

for f in $data; do
    [ -f "$code_dir/$f" ] || missing="$missing $f"
done
for f in $wrappers; do
    [ -f "$here/$f" ] || missing="$missing unix/$f"
done
[ -z "$missing" ] || die "missing from the source tree:$missing"

[ -n "$toolbox" ] || die "no toolbox/ package in $code_dir; the .py tools import it"
[ -n "$lsjj" ] || die "no lsjj/ package in $code_dir; make_cfp builds LSJJ.DAT with it"

# ---------------------------------------------------------------- copying

# One place where the dry run and the real run part company, so that
# --dry-run exercises every path decision and only the effects are
# suppressed.
run() {
    if [ "$dry" = 1 ]; then
        echo "  would: $*"
    else
        "$@" || die "failed: $*"
    fi
}

# place <source> <destination> <mode>
#
# The chmod u+w is for reinstalling over an install this script has
# already locked down, and over one make_cfp has been through: it puts the
# decks back to a-w on the way out, and cp cannot overwrite a file it may
# not write.  Only the owner can do this, which is the same person who can
# overwrite the file at all, so it adds no permission that was not already
# there.  Silent when it fails -- the cp that follows is the real test and
# gives the better message.
place() {
    if [ "$dry" = 1 ]; then
        echo "  would: install $1 -> $2 ($3)"
        return 0
    fi
    [ -e "$2" ] && chmod u+w "$2" 2>/dev/null
    cp "$1" "$2" || die "could not copy $1 to $2"
    chmod "$3" "$2" || die "could not set mode $3 on $2"
}

echo "Installing into $prefix"
[ "$dry" = 1 ] && echo "(dry run -- nothing will be written)"
echo

# -------------------------------------------------------- the data decks

# Before anything is copied, so that a failure here costs nothing, and so
# that what gets installed is the deck this machine's RCG has just written.
cfp_missing=""
for f in $cfp_decks; do
    [ -f "$code_dir/$f" ] || cfp_missing="$cfp_missing $f"
done

if [ -n "$cfp_missing" ] || [ "$rebuild_cfp" = 1 ]; then
    if [ "$dry" = 1 ]; then
        echo "data decks"
        echo "  would: $here/make_cfp   (in a scratch directory)"
        echo
    else
        [ -w "$code_dir" ] ||
            die "cannot build the data decks: $code_dir is not writable"

        # A directory of our own, removed however make_cfp leaves it: see
        # the header.  mktemp -d and not a fixed name, so that two installs
        # cannot collide over the ING11 each of them writes.
        scratch=$(mktemp -d "${TMPDIR:-/tmp}/cowan-cfp.XXXXXX") ||
            die "could not create a working directory for make_cfp"

        # COWAN_NO_CLEAR keeps make_cfp's screen clear from wiping the
        # install's own output; its banner is wanted, the clear is not.
        ( cd "$scratch" || exit 1; COWAN_NO_CLEAR=1 "$here/make_cfp" )
        rc=$?
        rm -rf "$scratch"
        [ "$rc" -eq 0 ] ||
            die "make_cfp failed (status $rc); the data decks were not built"
        echo
    fi
fi

# 755 on the directories: readable and traversable by everyone, writable
# only by the owner.  That is what makes the read-only guarantee hold for
# other users of a system-wide install.
for d in "" /unix /toolbox /lsjj; do
    run mkdir -p "$prefix$d"
    run chmod 755 "$prefix$d"
done

echo "executables (755)"
for f in $executables; do
    place "$code_dir/$f" "$prefix/$f" 755
done

echo "wrappers (755)"
for f in $wrappers; do
    place "$here/$f" "$prefix/unix/$f" 755
done

# 644 and not 755, although the .py files carry a shebang and are
# executable in the repository: the supported way to run them is the
# wrapper in unix/, which is what the README's command lines name and what
# finds an interpreter.  Installing them non-executable keeps that the
# only entry point instead of leaving two that can drift apart.
echo "python tools (644)"
for f in $tools; do
    place "$code_dir/$f" "$prefix/$f" 644
done
for f in $toolbox; do
    place "$code_dir/toolbox/$f" "$prefix/toolbox/$f" 644
done
for f in $lsjj; do
    place "$code_dir/lsjj/$f" "$prefix/lsjj/$f" 644
done

# The data files last, and 644 last of all, because that is the mode the
# whole arrangement turns on.  644 and not 444: the owner keeps write
# permission so that make_cfp can be run without a chmod when the package
# is installed in a home directory.  Other users get read-only, which is
# what an ordinary RCG run needs -- libgfortran downgrades the refused
# read-write open of STATUS='UNKNOWN' to read-only and carries on.
#
# The one tolerated absence, and it is deliberately narrow: only
# LSJJ.DAT, and only with a line saying so.  Everything else missing here
# is a broken source tree or a make_cfp that did not do its job, and
# place() should die on it as it always has.
echo "data (644)"
lsjj_deck=1
for f in $data $cfp_decks; do
    if [ ! -f "$code_dir/$f" ] && [ "$f" = LSJJ.DAT ]; then
        lsjj_deck=0
        echo "  LSJJ.DAT was not built -- not installed (see the warning above)"
        continue
    fi
    place "$code_dir/$f" "$prefix/$f" 644
done

# ------------------------------------------------------------------ done

echo
if [ "$dry" = 1 ]; then
    echo "Dry run complete.  Nothing was written."
    exit 0
fi

echo "Installed."
echo
echo "Put the wrappers on PATH -- they are what the README's command lines name:"
echo
echo "    export PATH=\"$prefix/unix:\$PATH\""
echo
case "${SHELL:-}" in
    *zsh) echo "Add that line to ~/.zshrc." ;;
    *)    echo "Add that line to ~/.bashrc (bash) or ~/.zshrc (zsh, the macOS default)." ;;
esac
echo "Then run rcn, rcn2, rcg, rce and the tools from a working directory of"
echo "your own.  Nothing is ever written into $prefix during a run."
echo
echo "Two Unix-specific settings, if a large case misbehaves:"
echo
echo "  OMP_STACKSIZE  RCG and RCE are threaded.  A thread's private stack"
echo "                 comes from OMP_STACKSIZE, whose default is a few MB --"
echo "                 smaller than the Windows default.  If a large case"
echo "                 crashes inside a parallel region and a serial build of"
echo "                 the same case does not, raise it:"
echo "                     export OMP_STACKSIZE=512M"
echo "                 and raise the main thread's own limit to match:"
echo "                     ulimit -s 65536"
echo "                 (These are starting values, not measured ones.)"
echo
echo "  .bss size      -fno-automatic puts every local in static storage, so"
echo "                 these executables carry a large .bss.  Should it ever"
echo "                 pass 2 GB -- after raising KMX in rcgpar.f, say -- the"
echo "                 small code model runs out of reach and the link fails"
echo "                 with a relocation error.  Add -mcmodel=medium then:"
echo "                     FFLAGS=-mcmodel=medium $prefix/unix/build.sh"
echo
echo "The data decks FOR072, FOR073, FOR074 and LSJJ.DAT are built on the"
echo "machine that installs the package, not distributed, and they were built"
echo "here.  The three CFP decks are portable as bytes between gfortran"
echo "targets, so they have to be built again only for an RCG compiled with a"
echo "compiler whose unformatted record layout differs from gfortran's."
echo "LSJJ.DAT is not portable in that way -- its coefficients come out of"
echo "LAPACK, and a different LAPACK moves a few dozen of them in the last"
echo "bit -- which is another reason it is built here rather than shipped."
echo "Either way the command is the same:"
echo
echo "    $prefix/unix/make_cfp"
echo
echo "run from a working directory of your own, and as the owner of $prefix."

# Repeated here because this is the text the user actually reads: the
# warning make_cfp printed is some way up the screen by now, and the
# consequence -- one feature missing from an install that otherwise
# succeeded -- is easy to scroll past.
if [ "$lsjj_deck" = 0 ]; then
    echo
    echo "LSJJ.DAT was NOT built, so true jj coupling (column 38 of the ING11"
    echo "main card set to 8) is not available.  Its generator needs numpy,"
    echo "which is the one thing in this package that has to be installed with"
    echo "pip.  Everything else works without it.  To get the deck:"
    echo
    echo "    python3 -m pip install numpy"
    echo "    $prefix/unix/make_cfp"
    echo
    echo "the second command run as the owner of $prefix, from a working"
    echo "directory of your own."
fi
