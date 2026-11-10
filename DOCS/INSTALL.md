# Installing the Cowan codes

[Back to the README](../README.md)

## Requirements

- **A 64-bit system:** Windows 10 or later, a current Linux distribution,
  or macOS. The distributed executables are x86-64, compiled with gfortran
  on 64-bit Windows 11, and will not run on a 32-bit system. On Linux and
  macOS the programs are compiled locally, again with gfortran and again
  only for 64 bits; the sources build equally well for 64-bit ARM (Apple
  Silicon, aarch64 Linux), but macOS has not been exercised -- see the
  note on testing in the [README](../README.md#platforms-and-testing).
- **Memory and disk:** at least 16 GB of free RAM for the array dimensions
  the distributed executables are built with; on a machine with less,
  recompile with smaller dimensions ([BUILDING.md](BUILDING.md)). A
  calculation may need up to 20 GB of free disk space, depending on the
  number and complexity of the configurations computed together.
- **Python 3**, on the PATH, for the utility codes. They use the standard
  library only, except `calc_TP` (below). They were developed and tested
  on Python 3.14 and use no syntax newer than 3.8. On Windows the `.bat`
  files call the `py` launcher; on Linux and macOS the wrappers look for
  `python3`, falling back to a `python` that is Python 3.
- **`numpy`**, the only Python package the package uses, for two things.
  `make_cfp` builds the `LSJJ.DAT` deck with `CODE/lsjj/mkdeck.py`, which
  needs it; there it is optional -- without it `make_cfp` says so, builds
  the CFP decks and finishes normally, and everything works except true
  *jj* coupling. The algebra cannot be done without it honestly: the
  coefficients are null spaces and degenerate eigenspaces of the coupling
  operator, and another implementation would return a different, equally
  valid basis -- no longer the deck that was validated against Cowan's
  own coefficients. And `calc_TP`, where it is not optional: that program
  pairs a trial's eigenvectors with the unvaried run's, a matrix product
  per *J* block, and it stops with the `pip install numpy` line rather
  than running slowly.
- **A Fortran compiler (gfortran)**: on Linux and macOS, to build the
  package in the first place, but only then -- once built, the programs
  link against the system libgfortran and libgomp and nothing further is
  required to run them. On Windows, only to **recompile** the sources:
  the gfortran runtime libraries are distributed in `CODE\`.

## Windows

(Linux and macOS users should skip to the next section.)

The installation procedure is designed for installation from a
compressed (ZIP) file downloaded in its entirety from the Web.

To install the package on your PC:

0. Install Python 3 from https://www.python.org/downloads/ , selecting
"Add python.exe to PATH" in the installer. That is all the utility codes
need; the four main programs need nothing at all beyond the DLLs shipped
in `CODE\`. If you intend to use true *jj* coupling, add `numpy` as well:

        pip install numpy

	It is the only Python package the package needs. It is optional for
	step 6 below, where `make_cfp` builds the `LSJJ.DAT` deck -- skip it
	and everything else still works, except true *jj* coupling -- but
	`calc_TP` requires it outright.
	Install a Fortran compiler only if you intend to recompile the
	sources -- see [BUILDING.md](BUILDING.md).
1. Extract all files from the COWAN.ZIP archive under directory `Cowan-main` PRESERVING THE
SUB-DIRECTORY STRUCTURE into the directory you choose to host the
package (typically `c:\COWAN`, but any location will do). Control check: 
your ~\COWAN directory must contain subdirectories CODE, FOR, WORK, DOCS,
and the README.md file - exactly like in Cowan-main of the zip archive.
2. Add the path of the Cowan executables to the system `Path` variable.
If you are installing to a network location, you may need to map it to a
drive letter first: the Windows command interpreter (cmd.exe) does not
accept UNC paths (those beginning with `\\`). Ask your system
administrator for help if necessary.

	- In the Windows search box on the taskbar, type `envir`; the best
	match should be 'Edit the system environment variables'. Open it and
	click `Environment Variables`. Select the variable named `Path`, in
	either the user or the system list, and click `Edit`, then `New`, and
	type the path of the COWAN\CODE directory -- for example
	`c:\COWAN\CODE`. See Windows Help for more on changing `Path`.
3. Open a Command Prompt window: type `cmd` in the Windows search box
and click the Command Prompt icon.
4. Make the window at least 120 characters wide. Right-click its title
bar, choose Properties, set Window Size Width (and the font size, colors
and anything else you want), click OK, choose 'Save properties for
future windows with same title' and click OK.
5. Type `PATH` and press Enter. If the COWAN\CODE directory does not
appear in the result, close the window and go back to step 2.
6. Create the binary data decks. They are not part of the distribution --
they are built here, on your machine, and RCG will not run without them.
Change to the COWAN\WORK directory, e.g. `cd c:\COWAN\WORK`, and execute
the make_cfp.bat batch file:

  `make_cfp <Enter>`

	This builds the binary CFP 'decks' (files `FOR072`, `FOR073`, `FOR074`
	in `COWAN\CODE`) and the LS-to-*jj* deck `LSJJ.DAT` beside them. It
	should take < half a minute, and you should receive a message

	`Done. New CFP binary files FOR072, FOR073, and FOR074 are created in "...COWAN\CODE\" directory.`
	`      The LS-to-jj deck LSJJ.DAT was rebuilt there as well.`

**If you do not see the first of those two lines, or any of the files `FOR072`, `FOR073`, `FOR074` has a zero size, RCG will not work! Contact A. Kramida to resolve the problem.**

	The second line is different. If it says instead that `LSJJ.DAT` was
	NOT rebuilt, read the warning above it: almost always it means
	`numpy` is not installed (step 0). RCG will run, and everything will
	work except true *jj* coupling. Install `numpy` and run `make_cfp`
	again whenever you want it.

7. Execute the test task defined by the sample input files `IN36` and
`IN2` in the `COWAN\WORK` directory:

        rcn  <Enter>
        rcn2 <Enter>
        rcg  <Enter>
        rce  <Enter>
	Each of these commands should report a `NORMAL EXIT`. If any of them
	fails, the installation was not successful; write to A. Kramida at
	akramida @ umd.edu.

	After a successful installation, the directory COWAN should contain
	`README.md`, `LICENSE.md`, `NOTICE.md`, and the following
	sub-directories:

	- `CODE` - all executable codes and the data files they need.
	- `DOCS` - this documentation; in `DOCS/Cowan/`, R. D. Cowan's own
					 manuals, unaltered: `RCN_DOC.txt`, `RCG_DOC.txt`,
					 `RCE_DOC.txt` and `readme.cowan.htm`.
	- `FOR`  - the Fortran sources of the main programs.
	- `WORK` - holds the sample input files `IN36` and `IN2`. Running the
					 four programs in turn on them is the installation test.
					 Its sub-directories `AI`, `AI/DIEL`, `Forbid`, `Mo15`,
					 `Pr3`, `W70` and `U90` hold further worked examples,
					 listed in [FILES.md](FILES.md).

	Besides the main programs RCN, RCN2, RCG and RCE, the package includes
	the utility codes described in [UTILITIES.md](UTILITIES.md). They are
	written in Python and need
	Python 3. Each is launched by a wrapper of the same name -- a `.bat`
	file in `CODE/` on Windows, a shell script in `CODE/unix/` on Linux
	and macOS -- which finds the code beside itself, so they can be run
	from any working directory once that directory is on the PATH.


## Linux and macOS

The four main programs are compiled from source here; everything else is
the same package. Five steps, none of them needing administrator rights.

1. **Install gfortran and Python 3.** Any reasonably recent gfortran
   will do; the reference build is GCC 13.

        sudo apt install gfortran python3        # Debian, Ubuntu
        sudo dnf install gcc-gfortran python3    # Fedora, RHEL
        brew install gcc python                  # macOS

   On macOS, gfortran comes from the Homebrew (or MacPorts) GCC package.
   Apple's own command-line tools provide `clang`, which is a C compiler
   and **not** a Fortran compiler; there is no Fortran compiler in Xcode.
   Homebrew installs the command as `gfortran`, and `build.sh` will find
   it on the PATH like any other.

   If you intend to use true *jj* coupling, add `numpy` as well:

        sudo apt install python3-numpy           # Debian, Ubuntu
        sudo dnf install python3-numpy           # Fedora, RHEL
        pip3 install numpy                       # or, in a virtualenv

   It is the only Python package the package needs, it is needed by one
   step only -- the `LSJJ.DAT` deck that `make_cfp` builds in step 4 --
   and it is optional. Skip it and everything else still works.

2. **Unpack the package** wherever you like -- `~/cowan` is a good
   choice. Do not put it in a directory whose name contains a space:
   RCG is told where its data files are through a one-line `COWAN.CFG`
   that it reads with `FORMAT(A)` and cuts at the first blank, so a path
   with a space in it cannot work. The build and install scripts check
   for this and say so.

   The zip GitHub builds from the repository carries no Unix permission
   bits, so nothing in `CODE/unix/` comes out of it executable. Restore
   that before anything else:

        chmod +x CODE/unix/*

   Without it the two commands below are refused as not executable.
   (`bash CODE/unix/build.sh` works around it for a single run.) Only
   the unpacked tree needs this: `install.sh` sets the modes itself on
   everything it installs.

3. **Compile:**

        CODE/unix/build.sh

   That builds all four programs into `CODE/`, named `RCN36K`, `RCN2K`,
   `rcg11kd` and `RCE20KD` -- the same names as the Windows executables
   without the `.exe`. RCG and RCE are built with OpenMP and are
   threaded. It takes well under a minute. `build.sh --serial` builds
   without OpenMP, and `build.sh --debug` builds unoptimized for `gdb`
   or `lldb`; `build.sh rcg11kd RCE20KD` rebuilds only the programs
   named. The flags are not to be adjusted casually -- see
   [BUILDING.md](BUILDING.md) for which of them are load-bearing.

4. **Install:**

        CODE/unix/install.sh

   This copies the package into `~/.local/lib/cowan` (or, run as root,
   into `/usr/local/lib/cowan`; `--prefix DIR` puts it anywhere, and
   `--dry-run` shows what would be copied without copying it). What is
   installed is the four executables, the wrappers, the Python tools and
   the data files the programs read; the `.bat` files, the `.exe` files,
   the Windows DLLs, the test suite and `perl_legacy/` are left behind.

   The point of installing rather than running from the unpacked tree is
   the file permissions it sets: directories and executables 755, data
   files 644. The package's functional division is that **the code
   directory is never written to during a run** -- every program runs in
   the working directory and writes only there. On Windows that is a
   convention the batch files hold up. Here it is enforced by the
   operating system, so that a stray `OPEN` in the code directory fails
   loudly instead of quietly corrupting the package.

   `install.sh` also builds the binary data decks, which is the one thing
   it does before copying anything. `FOR072`, `FOR073` and `FOR074` are
   not part of the distribution; `RCG` writes them itself, from
   `CODE/ING11.CFP`, and it will not run without them. `LSJJ.DAT`, the
   LS-to-*jj* deck, is built in the same run, by `CODE/lsjj/mkdeck.py`
   from the same input. This is the counterpart of step 6 of the Windows
   procedure, done for you: the script runs `CODE/unix/make_cfp` from the
   unpacked tree, which is the last moment `CODE/` is certain to be
   writable, in a scratch directory of its own. It takes well under a
   minute and asks nothing. `--rebuild-cfp` forces it even when the decks
   are already there.

   `LSJJ.DAT` is the one file that may still be missing afterwards: its
   generator needs `numpy`. If `numpy` is not installed, `make_cfp` says
   so, builds the three CFP decks, and the install finishes normally --
   `install.sh` repeats the warning in its closing summary. Everything
   works except true *jj* coupling.

   The three CFP decks are portable as bytes between gfortran targets, so
   a deck built on any of the three supported systems is read unchanged
   on the other two; they have to be built again only for an `RCG`
   compiled with a compiler whose unformatted record layout differs from
   gfortran's. `LSJJ.DAT` is not portable in that way: its coefficients
   come out of LAPACK, and a different LAPACK moves a few dozen of them
   in the last bit. Nothing `RCG` prints moves with them, but it is
   another reason the deck is built where it is read. Either way the
   remedy is the same -- run `make_cfp` by hand, from a working directory
   of your own, and as the owner of the installed tree.

5. **Put the wrappers on the PATH.** `install.sh` prints the exact line;
   it is

        export PATH="$HOME/.local/lib/cowan/unix:$PATH"

   Add it to `~/.bashrc` (bash) or `~/.zshrc` (zsh, the macOS default),
   open a new terminal, and check with `which rcg`.

Then run the installation test exactly as step 7 of the Windows
procedure does, in a working directory of your own holding copies of
`WORK/IN36` and `WORK/IN2`:

        rcn
        rcn2
        rcg
        rce

Each should report `NORMAL EXIT`. **Every other command line in this
documentation works unchanged**, on all three systems: the wrappers carry the
same names as the `.bat` files, they take the same arguments, and they
print the same banners.

One more thing is worth knowing about, though it is not normally needed:

- **Thread stack size.** RCG and RCE are threaded, and a thread's
  private stack comes from `OMP_STACKSIZE`, whose default is a few MB --
  smaller than the Windows default. If a very large case crashes inside
  a parallel region and a serial build of the same case does not, raise
  it with `export OMP_STACKSIZE=512M` and raise the main thread's own
  limit to match with `ulimit -s 65536`.

