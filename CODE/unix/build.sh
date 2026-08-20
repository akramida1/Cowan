#!/usr/bin/env bash
#
# Build the Cowan package on Linux or macOS.
#
#   ./build.sh                  all four programs; RCG and RCE threaded
#   ./build.sh --serial         the same, without OpenMP
#   ./build.sh --debug          -O0 -ggdb3, serial, for gdb or lldb
#   ./build.sh rcg11kd RCE20KD  only the programs named
#
# This is the Unix counterpart of gfortran_O3.bat, gfortran_O3_OMP.bat and
# gfortran_dbg.bat.  Those take one source stem per call because that is
# how the author drives them by hand; on Unix there is no reason not to do
# all four in one go, so this is one script with two switches.
#
# The executables land in CODE/, named as the batch files name them minus
# the .exe: RCN36K, RCN2K, rcg11kd, RCE20KD.  gfortran's -o leaves no
# suffix here, and the wrappers and the Python tools look for exactly
# those names beside CODE/.
#
# The flag sets are the ones recorded in HANDOFF.md section 3 and
# CLAUDE.md section 2.  Two of them are load-bearing and are not to be
# relaxed to buy speed:
#
#   -fno-automatic MUST accompany -fopenmp.  Without it -frecursive puts
#   locals on the stack and the numbers change.  This script refuses to
#   assemble a command line that has one and not the other.
#
#   -ffast-math and -march=native are never added.  The first alters the
#   eigensolver's arithmetic; the second buys about 1% and makes the
#   executable unfit to distribute.
#
# -fallow-invalid-boz is passed to every program, though gfortran_O3.bat
# passes it only to the OpenMP build: gcc 13 compiles Cowan's BOZ
# constants without it, later ones do not, and a Unix user is far more
# likely to have gcc 14 or 15.  It downgrades an error to a warning and
# changes no code.
#
# Environment:
#   FC          compiler to use, if the probe below picks the wrong one
#   FFLAGS      extra flags, appended last; screened for the two above
#   COWAN_SRC   the FOR/ directory, if this script has been moved away
#               from the repository

set -eu

die() {
    echo "build.sh: $*" >&2
    exit 1
}

here=$(cd "$(dirname "$0")" && pwd)          # CODE/unix
code_dir=$(cd "$here/.." && pwd)             # CODE
src_dir=${COWAN_SRC:-$code_dir/../FOR}

[ -d "$src_dir" ] || die "no Fortran sources at $src_dir (set COWAN_SRC)"
src_dir=$(cd "$src_dir" && pwd)

# ---------------------------------------------------------------- options

serial=0
debug=0
wanted=""

for arg in "$@"; do
    case "$arg" in
        --serial)  serial=1 ;;
        --debug)   debug=1; serial=1 ;;
        -h|--help) sed -n '3,8p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        -*)        die "unknown option $arg" ;;
        *)         wanted="$wanted $arg" ;;
    esac
done

# ------------------------------------------------------------- compiler

# On Linux the compiler is gfortran.  On macOS, Homebrew's gcc installs
# gfortran-15, gfortran-14 and so on, and often no unsuffixed gfortran at
# all, so a plain `command -v gfortran` is not enough to go on.
if [ -n "${FC:-}" ]; then
    command -v "$FC" >/dev/null 2>&1 || die "FC is set to $FC, which is not on PATH"
else
    FC=""
    for candidate in gfortran gfortran-15 gfortran-14 gfortran-13 gfortran-12 gfortran-11; do
        if command -v "$candidate" >/dev/null 2>&1; then
            FC=$candidate
            break
        fi
    done
fi

if [ -z "$FC" ]; then
    echo "build.sh: no gfortran found on PATH." >&2
    case "$(uname -s)" in
        Darwin) echo "  Install one with:  brew install gcc" >&2
                echo "  It arrives as gfortran-15 (or similar), which this script looks for." >&2 ;;
        *)      echo "  Install one with:  sudo apt install gfortran        (Debian, Ubuntu)" >&2
                echo "                     sudo dnf install gcc-gfortran    (Fedora, RHEL)" >&2 ;;
    esac
    echo "  Or set FC to the compiler you want:  FC=/path/to/gfortran ./build.sh" >&2
    exit 1
fi

# ---------------------------------------------------------------- flags

flags="-std=legacy -fshort-enums -ftracer -fno-backslash"
flags="$flags -fallow-argument-mismatch -fcray-pointer -fallow-invalid-boz"
flags="$flags -fno-automatic"

# -malign-double is an x86 option; on Apple Silicon and on ARM Linux
# gfortran rejects it outright.
case "$(uname -m)" in
    arm64|aarch64) ;;
    *) flags="$flags -malign-double" ;;
esac

if [ "$debug" = 1 ]; then
    # gfortran_dbg.bat's flags.  -fcheck=bounds is deliberately absent:
    # the EQUIVALENCE and dummy-array style this package is built on
    # trips it on correct code.
    flags="$flags -O0 -ggdb3 -fbacktrace"
else
    flags="$flags -O3"
fi

# Anything the user adds goes last, but not these.
for word in ${FFLAGS:-}; do
    case "$word" in
        -ffast-math|-funsafe-math-optimizations|-Ofast)
            die "$word changes the eigensolver's arithmetic; results would not be bit-exact" ;;
        -fautomatic|-frecursive)
            die "$word undoes -fno-automatic; locals must stay in static storage" ;;
    esac
done
flags="$flags ${FFLAGS:-}"

# The weld, asserted rather than assumed: every command line this script
# runs goes through here first.
check_weld() {
    case " $1 " in
        *" -fopenmp "*)
            case " $1 " in
                *" -fno-automatic "*) ;;
                *) die "refusing to build with -fopenmp and without -fno-automatic" ;;
            esac ;;
    esac
}

# ------------------------------------------------------------- programs

# name : source file : threaded
#
# RCN and RCN2 have no parallel regions, so -fopenmp would only cost them
# a runtime dependency.  The source spellings are the ones on disk: the
# batch files get away with %1.F for all four because Windows does not
# distinguish rcg11kd.f from rcg11kd.F, and here it does.
programs="RCN36K:RCN36K.F:0 RCN2K:RCN2K.F:0 rcg11kd:rcg11kd.f:1 RCE20KD:RCE20KD.F:1"

# INCLUDE lines are resolved by gfortran against the source's own
# directory, and are matched case-sensitively here but not on Windows.
# Checking them up front turns a wall of "Cannot open include file" into
# one line that says which file and what it should have been called.
check_includes() {
    local names inc actual bad=0
    names=$(grep -i "^ *include *'" "$1" 2>/dev/null |
            sed "s/.*[Ii][Nn][Cc][Ll][Uu][Dd][Ee] *'\([^']*\)'.*/\1/" |
            sort -u || true)
    for inc in $names; do
        [ -e "$src_dir/$inc" ] && continue
        actual=$(ls "$src_dir" | grep -ix -- "$inc" || true)
        if [ -n "$actual" ]; then
            echo "  $(basename "$1") includes '$inc', but the file is named $actual" >&2
        else
            echo "  $(basename "$1") includes '$inc', which is not in $src_dir" >&2
        fi
        bad=1
    done
    return $bad
}

built=""
for entry in $programs; do
    name=${entry%%:*}
    rest=${entry#*:}
    source=${rest%%:*}
    threaded=${rest#*:}

    if [ -n "$wanted" ]; then
        case " $wanted " in
            *" $name "*) ;;
            *) continue ;;
        esac
    fi

    [ -f "$src_dir/$source" ] || die "$src_dir/$source is missing"

    if ! check_includes "$src_dir/$source"; then
        die "$name cannot be compiled until the INCLUDE names above match the files on disk"
    fi

    this=$flags
    if [ "$threaded" = 1 ] && [ "$serial" = 0 ]; then
        this="$this -fopenmp"
    fi
    check_weld "$this"

    echo "==> $name"
    echo "    $FC $this -o $code_dir/$name $src_dir/$source"
    # gfortran's one expected warning is that -fno-automatic overrides
    # -frecursive.  Anything else on stderr is worth reading, so unlike
    # the batch files this does not divert it to a file.
    if ! (cd "$src_dir" && "$FC" $this -o "$code_dir/$name" "$src_dir/$source"); then
        die "$name failed to compile"
    fi
    built="$built $name"
done

[ -n "$built" ] || die "nothing to build; known programs are RCN36K RCN2K rcg11kd RCE20KD"

echo
echo "Built in $code_dir:$built"
if [ "$serial" = 0 ]; then
    echo "RCG and RCE are threaded; set OMP_NUM_THREADS to choose how many."
    echo "Thread count does not change the numbers, only the time."
fi
