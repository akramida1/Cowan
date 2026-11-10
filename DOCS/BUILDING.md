# Building the Fortran

[Back to the README](../README.md)

## Recompiling the Fortran

### On Windows

The distributed executables were compiled with

    GNU Fortran (MinGW-W64 x86_64-ucrt-posix-seh, built by Brecht
    Sanders, r8) 13.2.0

from the **WinLibs** distribution of GCC. To install it:

1. Download `winlibs-x86_64-posix-seh-gcc-13.2.0-mingw-w64ucrt-11.0.1-r8.zip`
   (216 MB) from the release page
   https://github.com/brechtsanders/winlibs_mingw/releases/tag/13.2.0posix-18.1.5-11.0.1-ucrt-r8
   Its SHA-256 is
   `8332dbad378edc4b3defc46fca76087a980dec28010adeb39dd8f967c2e75dda`.
2. Extract it to `C:\`. The archive holds a single top-level `mingw64\`
   folder, so this gives `C:\mingw64` (762 MB on disk). No installer and
   no administrator rights are involved.
3. Append `C:\mingw64\bin` to the `Path` variable, exactly as you
   appended `COWAN\CODE` in step 2 of the
   [Windows installation](INSTALL.md#windows).
4. Check the result with `gfortran --version`. It must report the
   version quoted above.

The same folder supplies `gdb`, the debugger that `gfortran_dbg.bat`
builds for, and `mingw32-make`.

Earlier releases of this package took gfortran from Strawberry Perl,
whose `C:\Strawberry\c` folder is this same WinLibs build: the compiler
binaries and all six runtime DLLs are byte-for-byte identical. Anyone
who already has Strawberry Perl installed therefore needs to change
nothing. Installing WinLibs alone is simply quicker, and it pins the
compiler version.

`winget install -e --id BrechtSanders.WinLibs.POSIX.UCRT` is a shorter
route but installs the current GCC (16.1.0 at the time of writing) and
cannot be pinned back to 13.2.0. Use it only if you do not need to
reproduce the distributed executables exactly; a different compiler
version is a different optimizer.

The sources are recompiled with the batch file `gfortran_O3.bat` (or
`gfortran_O3_OMP.bat` for a multithreaded build, or `gfortran_dbg.bat`
for debugging with `gdb`), described in [FILES.md](FILES.md).

### On Linux and macOS

There is nothing to recompile *for* -- compiling is the installation, and
`CODE/unix/build.sh` does it; see [INSTALL.md](INSTALL.md#linux-and-macos).
Use it again after changing a source file or an array dimension, naming
the programs you want rebuilt if you do not want all four. `FFLAGS=...`
in the environment appends flags to what it uses, and `FC=...` names a
compiler other than the `gfortran` it finds on the PATH.

### The flags, and which of them matter

The two build systems use the same flag set, and three points in it are
load-bearing rather than a matter of taste:

- **`-fno-automatic` must accompany `-fopenmp`.** Without it, `-fopenmp`
  implies `-frecursive`, which moves local variables to the stack, and
  the computed numbers change. `build.sh` refuses to assemble a command
  line that has one without the other.
- **`-ffast-math` is never added.** It alters the arithmetic of the
  eigensolver.
- **`-march=native` is not used** for anything meant to be distributed.
  It is worth about 1% and makes the executable unfit to run on another
  machine.

Every other flag is ordinary optimization. RCG and RCE are threaded, and
a threaded build is bit-identical to a serial one at any thread count,
because no parallel region reassociates a floating-point sum.


## Compiling for other operating systems

This section is about systems other than Windows, Linux and macOS, or
about compilers other than gfortran. On those three with gfortran,
nothing below has to be done: `gfortran_O3_OMP.bat` and
`CODE/unix/build.sh` build the package as it stands.

**Read this first if you are thinking of another compiler.** Building
with anything other than gfortran will **drop the parallelization**. RCG
and RCE are threaded with OpenMP as gfortran implements it, and the build
that ships is a `-fopenmp` build; a compiler without OpenMP, or one whose
OpenMP the sources have not been checked against, gives a serial
executable. On a small case that costs nothing. On a large one it is the
difference between a run and an overnight run: **expect a slowdown of
more than a factor of ten** on a many-core machine, since the threading
scales with the cores it is given. Weigh that against whatever the other
compiler is being chosen for.

- **How `RCG` finds its data files.** `RCG` has to locate `SENIOR`,
`LSJJ.DAT`, `FOR072`, `FOR073` and `FOR074`, which live with the
executable and not in the working directory. This is the mechanism that
makes the package's code/workspace split work -- the code directory is
never written to during a run, and every program runs in, and writes
only into, the working directory. It is governed by the `MSDOS` flag in
the first lines of `RCG`, where those files are opened: with `MSDOS=1`,
which is how the sources ship, `RCG` reads the package directory from a
one-line `COWAN.CFG` file that the `rcg` wrapper writes into the working
directory just before the run and deletes just after.

	**This design is implemented, and verified, for Windows, Linux and
	macOS built with gfortran, and for no other combination.** The
	`COWAN.CFG` line is read with `FORMAT(A)`, cut at its first blank,
	and then scanned backwards for a backslash so that a Windows path
	can be truncated after the last one; a Unix path contains no
	backslash, so that scan falls through and the path is used whole.
	The wrappers on each system -- `rcg.bat` and `CODE/unix/rcg` --
	write the line in the form each needs. On any other operating
	system, or with a compiler whose file-name handling differs, expect
	to have to do this yourself: set `MSDOS=0` and put the package
	directory, ending in the system's own path separator, in the
	hard-coded `FIL=` assignment in the `ELSE` branch beside it. One
	consequence of a path with a blank in it, on any system, is that
	`RCG` cuts the path there and then fails to find its data files;
	install the package somewhere whose full path has no spaces.

- The other change to look at is the system-dependent timing routine:
`SECONDS()` in `RCN`, `RCN2` and `RCG`, which returns a time in seconds,
and `CLOCK()` in `RCE`, which returns one in minutes. (`RCG` has a
`CLOCK()` as well, but it only divides `SECONDS()` by 60 and is not
machine-dependent.) Several variants are supplied in each; uncomment the
one that suits your system and comment out the code that is active by
putting a 'C' in the first column of its lines. Do this in every program
of the package. If you do not know the timing call for your system,
simply make the routine return zero -- the times are printed for
information only and nothing depends on them.

	On Linux and macOS nothing in fact has to be changed. `RCG` and
	`RCE` already call the Fortran standard intrinsic `SYSTEM_CLOCK`,
	which needs no library and behaves the same on all three systems;
	`RCN` and `RCN2` call `SECNDS` and `ETIME`, which are GNU
	extensions that `gfortran` also provides on Linux and macOS. The
	commented blocks headed *Use this code for Linux or macOS* in the
	four sources give the standard-intrinsic equivalents for anyone
	compiling with something other than `gfortran`, and record the one
	trap: `SYSTEM_CLOCK` must be given `INTEGER*8` arguments, because
	`gfortran` returns millisecond resolution for `INTEGER*4` ones and
	nanosecond resolution (`clock_gettime`) for `INTEGER*8`. Do not
	substitute `CPU_TIME` or `ETIME` in `RCG` or `RCE`: both return CPU
	time summed over the threads, which in an OpenMP build overstates
	the elapsed time by roughly the thread count.

	Compiling may not be enough by itself to run `RCG`. The CFP
	(coefficients of fractional parentage) decks `FOR072`, `FOR073` and
	`FOR074` are unformatted files. They are not distributed -- they are
	built on the installed machine by `RCG` itself -- so a deck already
	built by a gfortran `RCG` will not be read by an `RCG` compiled with
	a compiler whose unformatted record layout differs, and has to be
	built again with the new executable. Build them with a special input
	file, which is supplied: copy
	`ING11.CFP` to `ING11`, delete the first 8 comment lines from the
	copy, and run `RCG`. Those comment lines state the minimum array
	dimensions `RCG` needs to handle the deck. `make_cfp` (`make_cfp.bat`
	on Windows, `CODE/unix/make_cfp` on Linux and macOS) does exactly
	this and puts the results where `RCG` will look for them.

	The same applies to `LSJJ.DAT`, in the same way and for the same
	reason: it too is an unformatted file whose record layout is the
	compiler's, it too is not distributed, and `make_cfp` builds it in
	the same run -- so a deck written under one compiler is refused, not
	misread, by an `RCG` built with another, and the remedy is to run
	`make_cfp` again with the new executable. It is less portable than
	the CFP decks rather than more: its coefficients come out of LAPACK
	through `numpy`, and a different LAPACK build moves a few dozen of
	them in the last bit, so two machines do not produce identical bytes
	even when their record layout agrees. Nothing `RCG` prints moves with
	that. No other program of the package needs anything beyond the
	change to the timing routine described above.

- The dynamic memory allocation functions `ALLOC`, `DALLOC`, and `RALLOC`
in `RCG` and `RCE` (introduced in 2026) are machine-specific. As they
ship, they declare their pointers `INTEGER*8` and call `MALLOC`, which is
right for 64-bit Windows, Linux and macOS with gfortran and needs no
change on any of the three. On other platforms, see the presently
commented-out versions of the machine-specific lines and uncomment the
ones that suit (there are versions for DEC, SUN, Linux and VAX). These
routines need Cray pointers, which is why `-fcray-pointer` is in the
flag set; a compiler without that extension will not build `RCG` or
`RCE` at all.

