# Cowan's atomic-structure codes, modernized for 64-bit systems

A. Kramida, Department of Astronomy, University of Maryland, College Park,
MD, USA -- akramida @ umd.edu

**Version 26.2.2**, November 10, 2026 -- see
[DOCS/CHANGES.md](DOCS/CHANGES.md). Earlier versions were distributed
from <https://doi.org/10.18434/T4/1502500>.

To see which version you have: on Windows, the Details tab of an
executable's Properties in Explorer, or in PowerShell
`(Get-Item CODE\RCE20KD.exe).VersionInfo.FileVersion`; on Linux and macOS,
the `VERSION=` line of `version.txt` in the installed tree.

R. D. Cowan's programs for atomic energy levels, radiative transition
wavelengths and probabilities, electron-impact excitation cross sections
and photoionization cross sections, written by him between 1962 and 1991
and described in R. D. Cowan, *The Theory of Atomic Structure and Spectra*
(University of California Press, Berkeley, 1981). The package holds the
four main programs -- RCN, RCN2, RCG and RCE, as Fortran sources and as
Windows executables -- the data files they read, and Python utility codes
that prepare their input and interpret their output. Several people have
contributed to them since; [NOTICE.md](NOTICE.md) says who. Write to the
author with any question about their operation.

## Platforms and testing

The package runs on **Windows, Linux and macOS**, with the same input
decks, output files, utility codes and command names on all three.
Executables are supplied for 64-bit Windows, with the gfortran runtime
libraries they need; on Linux and macOS the four programs are compiled
from `FOR/` by one command.

The Windows build is the one the author uses daily. The Linux build has
been verified end to end on x86-64 (Ubuntu 24.04, gfortran 13.3): the
Pr III case gives a byte-identical `TAPE2E` on both systems. **macOS has
not been exercised** -- it is expected to work, not tested; please report
what you find.

## Requirements

- A 64-bit system with at least 16 GB of RAM. The executables are built
  for array dimensions that suit most tasks; for larger ones, recompile
  ([DOCS/BUILDING.md](DOCS/BUILDING.md)).
- Python 3 for the utility codes, and `numpy` for true *jj* coupling
  (optional) and for `calc_TP` (required).
- gfortran on Linux and macOS; on Windows only to recompile.

No path in the package may contain a space. Details are in
[DOCS/INSTALL.md](DOCS/INSTALL.md#requirements).

## Installing

**Windows** -- in full in [DOCS/INSTALL.md](DOCS/INSTALL.md#windows):

1. Install Python 3, with "Add python.exe to PATH"; `pip install numpy`.
2. Extract the zip archive, keeping its directories, e.g. to `c:\COWAN`.
3. Add `c:\COWAN\CODE` to the `Path` variable.
4. In `c:\COWAN\WORK`, run `make_cfp`. It builds the binary data decks
   RCG cannot run without.

**Linux and macOS** -- in full in
[DOCS/INSTALL.md](DOCS/INSTALL.md#linux-and-macos):

    chmod +x CODE/unix/*
    CODE/unix/build.sh
    CODE/unix/install.sh

and put on the PATH the line `install.sh` prints at the end.

## Running

From a directory holding the input files `IN36` and `IN2`:

    rcn
    rcn2
    rcg
    rce

Each program reports `NORMAL EXIT`; run in `WORK/`, this is the
installation test. After a least-squares fit with RCE, `update11` carries
the fitted parameters back into RCG's input and runs RCG again. On the
first RCE run of a fit, `rce --ground` puts the calculated ground level
at zero. Worked examples, each with its own `README.txt`, are in the
subdirectories of `WORK/`.

**Comparing with output of versions before 26.2.1:** the physical
constants were brought up to CODATA 2022 in 26.2.1, and some printed
quantities moved slightly -- autoionization rates by +0.05 %. See
[DOCS/CHANGES.md](DOCS/CHANGES.md#the-fundamental-constants).

## Documentation

| Document | Contents |
|---|---|
| [DOCS/INSTALL.md](DOCS/INSTALL.md) | requirements and installation |
| [DOCS/MANUAL_ADDITIONS.md](DOCS/MANUAL_ADDITIONS.md) | what differs from Cowan's manuals: running the programs, file formats, new options |
| [DOCS/UTILITIES.md](DOCS/UTILITIES.md) | the Python utility codes |
| [DOCS/BUILDING.md](DOCS/BUILDING.md) | recompiling; other compilers and operating systems |
| [DOCS/CHANGES.md](DOCS/CHANGES.md) | release notes, program versions, bugs fixed |
| [DOCS/FILES.md](DOCS/FILES.md) | every file in the package |
| [DOCS/Cowan/](DOCS/Cowan/) | R. D. Cowan's original manuals for RCN/RCN2, RCG and RCE, unaltered |

## License

The package is not covered by a single licence: see
[LICENSE.md](LICENSE.md) and [NOTICE.md](NOTICE.md).
