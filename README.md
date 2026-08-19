# COWAN'S ATOMIC-STRUCTURE CODES, MODERNIZED FOR 64-BIT WINDOWS

                            by A. Kramida
      (Adapted also for IBM RS/6000, VAX/VMS, Linux, and other platforms)

               Department of Astronomy, University of Maryland
                           College Park, MD, USA
                    e-mail: akramida @ umd.edu

                 Earlier versions distributed from
                https://doi.org/10.18434/T4/1502500

        Current version release date: August 19, 2026

### CONTENTS

The package consists of a set of computer programs for calculation
of energy levels, radiative transition wavelengths and probabilities,
electron-impact excitation cross sections, and photoionization cross
sections. The programs were originally written by Robert D. Cowan
between 1962 and 1991
(see R. D. Cowan, The Theory of Atomic Structure and Spectra, University
of California Press, Berkeley - Los Angeles - London, 1981).
Several people have contributed to them since. Complete instructions on
how to use the programs are in the documentation files distributed with
this package, in `DOCS\`. Write to A. Kramida at the address above with
any question about their operation.

Besides the four main programs -- RCN, RCN2, RCG and RCE, supplied both
as executables and as Fortran sources -- the package contains a number
of further files: data files the main programs read, and utility codes
that prepare their input and interpret their output. Follow the
installation procedure below to have the package work as intended. Every
file in the package is listed in the addendum.

The executable files supplied here are built for array dimensions that
suit most tasks. If your calculation does not fit them, recompile the
programs and, if necessary, rebuild the data files used by RCG (see the
instructions below). On a PC, use the gfortran compiler from the WinLibs
distribution; the runtime libraries it needs are already in `CODE\`, so
nothing has to be installed to *run* the package. The codes also compile
on other platforms -- see COMPILING FOR OTHER OPERATING SYSTEMS below.

**The utility codes are Python and no longer use Perl at all**, and Perl
is no longer part of the toolchain either. The Perl versions the Python
replaced are kept, unchanged and unused, in `CODE\perl_legacy\`: they are
the reference the Python was accepted against, they are not on the PATH,
and nothing calls them.

---

### CURRENT VERSIONS

| Code                  | Date     |
|-----------------------|----------|
| RCN  ver. 36k         | Jan 2026 |
| RCN2 ver. 2k          | May 2018 |
| RCG  ver. 11kd        | Jan 2026 |
| RCE  ver. 20kd        | Jan 2026 |
| **Python utility codes:** |        |
| printout.py                 | Aug 2026 |
| conv_out.py                 | Aug 2026 |
| add_eav.py                  | Aug 2026 |
| s11.py                      | Aug 2026 |
| subst_params_ing11.py       | Aug 2026 |
| subst_exp_levs.py           | Aug 2026 |
| reorder_ing11.py            | Aug 2026 |
| scale_param.py              | Aug 2026 |
| conf_list.py                | Aug 2026 |
| read_aa.py                  | Aug 2026 |
| calc_TP_M1E2.py             | Aug 2026 |
| config_purity.py            | Aug 2026 |
| make_levels_list_prepub.py  | Aug 2026 |
| make_levels_list_pub.py     | Aug 2026 |
| toolbox\ (shared modules)   | Aug 2026 |

Each of these (except the last four) replaces a Perl program of the same name, kept in
`CODE\perl_legacy\`. Two Perl files have no Python namesake:
`conv_cowan.pl` was `conv_out.pl`'s own library and moved with it, and
`vacair.pl` is now `toolbox\vacair.py`. The `config_purity`
and two `make_levels_list` programs, and the `toolbox` library itself,
are new: they never had a Perl version.

### SYSTEM REQUIREMENTS

This version of the package was compiled with gfortran on 64-bit
Windows 11. The executable files will not run on a 32-bit operating
system.

The array dimensions the distributed executables are built with require
at least 16 GB of RAM. On a machine with less memory, or on another
operating system, recompile the sources with smaller dimensions and the
system-specific changes described below.

### INSTALLATION (Only for 64-bit Windows-based PCs; mainframe and workstation users can skip this section)

The installation procedure is designed for installation from a
compressed (ZIP) file downloaded in its entirety from the Web.

To install the package on your PC:

0. Install Python 3 from https://www.python.org/downloads/ , selecting
"Add python.exe to PATH" in the installer. That is all the utility codes
need; the four main programs need nothing at all beyond the DLLs shipped
in `CODE\`. Install a Fortran compiler only if you intend to recompile
the sources -- see RECOMPILING UNDER WINDOWS below.
1. Extract all files from the COWAN.ZIP archive PRESERVING THE
SUB-DIRECTORY STRUCTURE into the directory you choose to host the
package (typically `c:\COWAN`, but any location will do).
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
6. Change to the COWAN\WORK directory, e.g. `cd c:\COWAN\WORK`.
Execute the make_cfp.bat batch file:

  `make_cfp <Enter>`

	This builds the binary CFP 'decks' (files `FOR072`, `FOR073`, `FOR074` in
	`COWAN\CODE`). It should take < half a minute, and you should receive
	a message 

	`Done. New CFP binary files for072, ... are created in ...COWAN\CODE\ directory`. 

**If you do not see this message, or any of the files `FOR072`, `FOR073`, `FOR074` has a zero size, RCG will not work! Contact A. Kramida to resolve the problem.**

7. Execute the test task defined by the sample input files `IN36` and
`IN2` in the `COWAN\WORK` directory:

        RCN  <Enter>
        RCN2 <Enter>
        RCG  <Enter>
        RCE  <Enter>
	Each of these commands should report a `NORMAL EXIT`. If any of them
	fails, the installation was not successful; write to A. Kramida at
	akramida @ umd.edu.

	After a successful installation, the directory COWAN should contain
	this `README.md` file, `license.txt`, and the following
	sub-directories:

	- `CODE` - all executable codes and the data files they need.
	- `DOCS` - R. D. Cowan's own manuals, unaltered: `RCN_DOC.txt`,
					 `RCG_DOC.txt`, `RCE_DOC.txt` and `readme.cowan.htm`.
	- `FOR`  - the Fortran sources of the main programs.
	- `WORK` - holds the sample input files `IN36` and `IN2`. Running the
					 four programs in turn on them is the installation test.
					 Its sub-directories `AI`, `AI\DIEL`, `Forbid`, `Pr3`,
					 `W70` and `U90` hold further worked examples, listed in
					 the addendum below.

	Besides the main programs RCN, RCN2, RCG and RCE, the package includes
	the utility codes described below. They are written in Python and need
	Python 3. Each is launched by a `.bat` file of the same name, which
	finds the code beside itself, so they can be run from any working
	directory once `COWAN\CODE` is on the PATH.


### RUNNING THE PROGRAMS
Programs are normally run in the following order:
`RCN`, `RCN2`, `RCG`, `RCE`.

On a PC these are the names of the batch files. Type them at the command
prompt from the directory holding the input files. Each program produces
the input for the next, and `RCG` and `RCE` are often run several times
after the first pass.

After a least-squares fit (LSF) with RCE, one normally wants to
recalculate the transition rates by running RCG with the fitted
parameters from the LSF output. In this PC version of the package the
easiest way to do that is to run the utility code:

`update11 2`

(if the `in36` file contains configurations of both parities) or

`update11`

(if there is only one parity in the calculation). Besides transferring the
fitted Slater parameters from the `RCE` output (RCEOUT, to be exact) into the 
`ing11` file and re-running `RCG` with these parameters, this does some 
additional useful things (see more details below in the explanation of the
update11 code).

### ADDITIONS TO THE MANUAL

1. The input/output file names are given in the table below.
		
	| Program | Input file(s)   | Output file(s) | Contents of the output file(s) |
	|---------|-----------------|----------------|--------------------------------|
	| RCN     | IN36    (ASCII) | OUT36   (ASCII)| Diagnostics                    |
	|         |                 | TAPE2N (binary)| Input for RCN2 (wavefunctions) |
	| RCN2    | IN2     (ASCII) | OUT2    (ASCII)| Diagnostics                    |
	|         | TAPE2N (binary) | ING11   (ASCII)| Input for RCG  (can be edited) |
	| RCG     | ING11   (ASCII) | OUTG11  (ASCII)| Diagnostics                    |
	|         |                 | OUTGINE (ASCII)| Input for RCE (can be edited)  |
	|         |                 | TAPE2E (binary)| Input for RCE                  |
	| RCE     | OUTGINE (ASCII) | OUTE    (ASCII)| Diagnostics                    |
	|         | TAPE2E (binary) | LEVELS1 (ASCII)|   -"-                          |
	|         |[RCEINP] (ASCII) | LEVELS2 (ASCII)|   -"-                          |
	|         |                 | LEVELS3 (ASCII)|   -"-                          |
	|         |                 | PARVALS (ASCII)|   -"-                          |
	|         |                 | RCEOUT  (ASCII)    | Can be renamed to and used as RCEINP |
	|         |                 | [RCEINP] (ASCII)   | Vector input file for RCE            |
	|         |                 | [RCEINP.HF] (ASCII)| Copy of RCEINP created first time you run RCE.^*)^|

	^*)^ Note: if you have a `RCEINP.HF` file from older runs with different 
	`IN2` or `IN36` files in the current directory, delete it before running `RCE`,
	so that the newly created `RCEINP.HF` would contain properly scaled (as 
	defined in `IN2`) HFR parameters.

	File `RCEINP` can be created by `RCE` or can be used by `RCE` as an
	input file, depending on the flag `INVECT` in the OUTGINE file (see below).

2. Compared to Cowan's original version, the `RCG` program
has two more input parameters that must be specified in the
"first control card" (usually, the first line in the file `ING11`).
They both are in the format `INTEGER*1` and are specified in
columns 77 and 78 (in place of the disabled parameter `ICTC`; see the
manual for `RCG`). Both default to 0, which gives the printout format R. D.
Cowan designed. If the second parameter, `IQS` (column 78 in
`ING11`), is set to a non-zero value *n* < 10, the quantum numbers of
the *n*-th open subshell appear in the parent term names in the output.
This matters when there are three or more open subshells, where Cowan's
original output left it ambiguous which subshell a number belonged to.
The first of the two parameters, `INWOL`, has no effect.

	If `IQS` equals zero and `INWOL` (column 77 in `ING11`) equals 1, then
	the alternate way of choice of the parent subshell is used. For
	questions, refer to the author of these two parameters, Yu. Ralchenko
	(present e-mail address yralchen @ umd.edu).

3. The parameter `DMIN` (columns 61-65 of the `RCG` input file `ING11`)
originally had to be a non-negative floating-point number less than 1.0,
and it set the minimum line strength retained in the `RCG` transition
list. There is a second use of it that the manual does not mention: if
it is greater than 1.0, only lines with *gA* >= 10^DMIN^ s^-1^ are
included. That option now works with every line sorting; originally it
did not work at all for a wavelength-sorted list. **Set `DMIN` to 0.00
for a forbidden-transition calculation** -- see `WORK\Forbid\README.txt`.

4. The designations of terms have changed slightly. In the original
version of `RCG` two different terms could end up with the same name,
even with the additional (e.g., seniority) quantum numbers in use. To
remove that ambiguity I added a distinguishing letter ('a', 'b', ...)
to such terms. For example, the term names `(3P) 3Da` and `(3P) 3Db` mean
different ^3^D terms: their parent terms have the same name but different
genealogy. For such cases, the user may have to manually trace the
history of the summation of quantum numbers that produced these terms,
using some additional printout in `OUTG11`. This task is automated by the
utility codes `printout.bat` and `conv_out.bat` (see below), which read the
summation history from `genealogy.txt` rather than reconstructing it.

5. Notation of the configuration-interaction parameters has been
changed. The first two symbols are the sequential numbers of the two
interacting configurations (the same as it was in Cowan's original RCG),
but now the letters 'A','B', ..., 'Z', 'a', ..., 'z' are used for
numbers 10, 11, ..., 35, 36, ..., 61, in place of the asterisk (`*`)
that Fortran writes when a value does not fit the format specified
(`I1`). In the latest version of Cowan's `RCG`, the format of
the CI parameters was extended from 8 to 10 characters, eliminating the
need for such encoding. However, I did not implement this extension
in the present version of the package.

6. The alternate-coupling scheme input files for `RCE` (other than
*LS*/*JJ*) are created by `RCG` if the following conditions are met:
	 1. `NLSMAX`, columns 23-25 of the first line of `ING11`, does not
			forbid them. Since the 2026 revision a **blank** field means
			"no limit" and permits couplings 3-7, and only an explicitly
			punched `0` forbids them; 999 also means "no limit" and is
			still accepted. This is the one field in which a blank no
			longer reproduces the pre-2026 behavior, and the change is
			deliberate: `KCPLD(I)=0` is what asks for coupling *I*, and 0
			is what a blank column gives, so a card that left both fields
			blank used to contradict itself.
	 2. Only one zero is present in columns 33-37 of the first line
			of `ING11` (`KCPLD(3..7)`): no output if non-zero.
	 3. If coupling schemes No. 5, 6 or 7 are needed, the number of
			open subshells `NOSUBC` must be not less than 3 (add more
			open inner shell configurations in `IN36` if it is not so!).

7. Two new optional input parameters have been introduced:
`E0MAX(1)` and `E0MAX(2)`, which restrict the list of lines included by `RCG`
in `OUTG11` so that the energy of the first-parity levels does not exceed
`E0MAX(1)` (if it is non-blank), and the energy of the second-parity
levels does not exceed `E0MAX(2)` (if it is non-blank). These parameters
can be specified in the optional **rescale card**, e.g., on the
additional line included at the top of `ING11` with zero in position 5
and floating-point numbers `E0MAX(1)` and `E0MAX(2)` in positions 66-75 and
76-85, respectively. The values are given in units of kK (10^3^ cm^-1^).

8. The parameter `ICRIT2` is disabled in `RCE20k`: its value has no
effect. The program exits after completing the fit, whether the
parameters were linked or unlinked, and never starts another run with
the next value of `LSQM`.

9. *Suspended* execution is not supported in `RCE20k`: intermediate
results cannot be saved after the *n*-th iteration and resumed from
there. The new `RCEOUT`/`RCEINP` files (see below) achieve almost the
same thing; the only cost is one extra initial iteration when a new run
is started from an `RCEOUT` copied into `RCEINP`. The `I216` parameter is
disabled in `RCE`, and its position in `OUTGINE` is taken by `IW6`, which
directs some output to the screen when it is negative.

10. Earlier versions of `RCE` refused to fit when every energy level
was flagged "experimental": the run printed `Iteration No 0` and stopped,
and, worse, wrote `RCEOUT` with every level starred, so feeding it back as
`RCEINP` also discarded the experimental energies. The 2026 revision
removed that restriction. A case in which every level is known is the
best-determined one there is, and `RCE` now fits it like any other,
whatever iteration count `OUTGINE` carries. The only conditions that still
end a run after the zeroth cycle are the two that leave nothing to fit: no
varied parameter at all, or no known level anywhere.

11. `RCE` accepts an alternative input file format. On its first run it
writes a new input file, `RCEINP`; on every run it also writes `RCEOUT`,
which has the same structure. This makes it far easier to change the
input -- to fix or free an energy level or a Slater parameter, for
example. `RCEOUT` can be edited in any text editor and renamed to
`RCEINP`, and it then becomes the input for the following iterations.
These features can be ignored and the old `OUTGINE` file used as before.

	To read `RCEINP` instead of `OUTGINE`, put a non-zero value in the
	25th column of the first line of `OUTGINE` (the parameter `INVECT`,
	"use vector-input file"). Every subsequent line of `OUTGINE` is then
	ignored, except the first control line of each parity section, and
	the corresponding data are read from `RCEINP`. Since the 2026
	revision `INVECT` also selects how levels are re-ordered after each
	iteration:

	| `INVECT` | Input | Level ordering after each iteration |
	|---|---|---|
	| 0 | `OUTGINE`; `RCEINP` is written for later fitting | full eigenvector recognition |
	| 1 | `RCEINP` | full eigenvector recognition |
	| 2 | `RCEINP` | eigenvectors matched only for levels tied to an experimental energy; all others ordered by increasing energy |
	| 3 | `RCEINP` | eigenvectors ignored; every level ordered by increasing energy (Cowan's original sorting) |

	Values 2 and 3 trade accuracy of the level matching for speed and are
	meant for cases too large for full recognition. The fit itself uses
	only the levels that carry an experimental energy, and mode 2 still
	matches those through their eigenvectors.

	`RCE` uses a new algorithm for sorting the eigenvalues. It keeps the
	correspondence between the level flags from one iteration to the next
	far better than Cowan's original options, but it does not guarantee
	the correct flags in every case. Where the coupling scheme is very far
	from pure *LS*, it may be necessary to increase the number of
	eigenvector components used to recognize a level. The default is 5; to
	raise it, increase the parameter `MAXC` in the source and recompile.
	Do the same if `RCE` terminates with 'Failed to identify the
	eigenvalues', which is very unlikely. In practice that message
	indicates diverging iterations, which is usually cured by fixing some
	of the Slater parameters.
	
	Two further new parameters appear on the first line of `OUTGINE`:
	`IORDER` and `ISAMPLE`, in columns 26-30 and 31-35. If `IORDER` is
	non-zero, the sampling algorithm starts reordering the eigenvectors
	only after *n* = `IORDER` iterations; before that they are sorted by
	increasing energy within each *J* group.

	If `ISAMPLE` is non-zero, the eigenvector list from the previous
	iteration is used as the sample for reordering. Otherwise the sample
	is always the list from the zeroth (more precisely, the `IORDER`-th)
	iteration, which comes from `RCEINP` or from `OUTGINE` according to
	`INVECT`.

	The same sampling algorithm assigns the level labels -- the unique
	term name each level carries -- if the parameter `CRIT1` is positive;
	its default is negative. `CRIT1` is in columns 61-65 of the line
	preceding the final '`   -1`' line of each parity section in
	`OUTGINE`, and has been there since version 20 of `RCE`.

**NOTE 1:**

Changing columns 16-30 of the first line of the _SECOND_ parity section
of `OUTGINE` has no effect. Only the very _FIRST_ line of `OUTGINE` is
examined for these options.

**NOTE 2:**

For `RCE` iterations to work properly, the configuration labels in
`IN36` must be no longer than 5 characters and must contain no spaces
and no dashes. Use primes, quotes, or asterisks to distinguish
configurations with different cores.

**Description of `RCEINP` and `RCEOUT` file formats:**

1. The line containing flags dealing with printing modes is the same as
   in `OUTGINE`. Additional features produced by the same print flags:
   a) Parameter `IPRNA` (if zero) makes all non-zero angular coefficients
      for radial integrals (for each "experimental" level) to appear in
      the output (`OUTE`), in much more readable form than it was before.
      (Originally, this parameter caused printing of all the matrices to 
      be diagonalized.)
   b) Parameter `IPRNV` (if zero) makes all non-zero derivatives of energy 
      levels over the Slater parameters to be printed (originally, _ALL_
      the derivative matrices were printed).
   c) Parameter `IPRNSQ` = 2 will suppress all the eigenvectors printout;
      `IPRNSQ` = 1 will suppress (as it was originally) only the output of
      "eigenvector components squared".
2. The iteration control line (same as in `OUTGINE`).
3. The group of lines describing the eigenvectors resembles the format
   of the `LEVELS1` file in *LS* coupling: the energy; an asterisk (`*`)
   for an "unknown" level or a blank for an "experimental" one; the
   calculated (fitted) energy; the *J* value; and then one or more groups
   describing the eigenvector components, each of them the signed
   amplitude (multiplied by 100 and rounded), the configuration name and
   the term name. Up to five components may be present.
4. The group of lines describing the Slater parameters resembles the
   format of the `OUTE` file: the parameter name, the parameter flag, the
   value, the maximum value for that parameter (usually 0.000, meaning no
   maximum), and the denominator if the `DENOM` flag is set. For
   configuration-interaction parameters the names of the two
   configurations are printed for reference. Do not change the columns
   the numbers occupy in the `RCEINP` file that the first run of `RCE`
   creates.

The same structure (1-4) follows for the second parity, if there is one.

### HARDWARE REQUIREMENTS (for PC users)

The package needs a 64-bit x86-64 PC with at least 16 GB of free RAM.
A calculation may need up to 20 GB of free disk space, depending on the
number and complexity of the configurations computed together.

### SOFTWARE REQUIREMENTS (for PC users)
Windows 10 or later.

**Python 3 must be installed** and its `py` launcher must be on the PATH;
that is all the utility codes need. They use the Python standard library
only -- there is nothing to `pip install`. They were developed and tested
on Python 3.14 and use no syntax newer than 3.8.

A Fortran compiler is needed only if you intend to **recompile** the
sources. Nothing else in the package requires one: the gfortran runtime
libraries are distributed in `CODE\`.


### SUMMARY OF THE MAJOR BUGS FIXED
- Multiple problems with memory allocation have been fixed in the 'kd'
versions of `RCG` and `RCE` (Jan 2026) by introducing dynamic memory 
allocation routines from `GRASP2K` and removing most of the 
`EQUIVALENCE` statements. Recompiling the codes with changed array
dimensions is no longer a major problem.
- A couple of small bugs have been fixed in `RCN`, compared to the 
previous official release of 2021: 1) the line that made the code to 
bypass the fine tuning of the convergence-managing block (adjustment 
of `alfmin` and `alfmax`) for ion charge > 1 has been commented out, as
convergence problems were observed for some highly charged ions, such 
as U^28+^ (for the 4f^9^ configuration). 2) Some highly charged complex 
ions require up to 450 iterations of the self-consistent field 
calculation. Thus, the maximum number of iterations has been increased
from 200 to 470, and the corresponding array dimensions have been 
adjusted accordingly.
- Incorrect identification of Rydberg-series configurations in `RCN2`
led to the appearance of non-zero *R~d~*^0^ configuration integrals in some
cases, e.g., when *n*sp^3^ and sp^3^*n*'s configurations are included in the
set, or when some core-excited configurations are present in addition
to a normal Rydberg series. See A. Kramida, Comput. Phys. Commun. **215**,
47 (2017); Erratum: **232**, 266 (2018).
- Multiple bug fixes made in the McGuinness version of `RCN` of 2004
have been incorporated in the present release of `RCN`. This solved poor
convergence problems in some cases, e.g., term-dependent calculations
in Ar I.
- Abnormal termination of `RCG` on large calculations in earlier
versions was caused by a mismatch between the declared sizes of some
`COMMON` blocks, which a bug in the Salford Fortran linker turned into
incorrect memory allocation. All `COMMON` blocks carrying the same name
are now declared with the same size.
- Several dimension checks have been added throughout `RCG`.
- Checks that the range of the term sorting keys is large enough were
added in several places in `RCG`. If it is not, the diagnostic
'Insufficient term sorting key' appears in the output, and `RCG` has to
be recompiled with a wider range (consult A. Kramida).
- An error that produced wrong _JJ_ notations in the output files was
fixed in `RCG` and `RCE`.
- An error that produced unexpected blank lines in the `OUTGINE` file
was fixed in `RCG`.
- Errors that sometimes put control characters into the `ING11` and
`OUTGINE` files were fixed in `RCN2` and `RCG`.
- `RCN2` did not recognize the use of parameters with illegal range
(`IABG`=2, see the description of the input file for `RCG`), and if it was
necessary, the `ING11` file had to be edited by hand. This is now
fixed.
- In some cases Cowan's original `RCG` produced enormous temporary files
(several gigabytes). They are now 10 to 100 times smaller, which speeds
up `RCG` and `RCE` by a factor of 2 to 5.
- In all programs, `real*4` constants and functions (`ABS`, `EXP`, `LOG`
etc.) have been replaced with `real*8` constants and corresponding
functions (`DABS`, `DEXP`, `DLOG` etc.). It is significant for some Fortran
compilers.
- Errors in `RCG`/`RCE` that led to wrong notation of levels in alternate
coupling schemes in `RCE` output files have been fixed.
- `RCE` gave up after the zeroth cycle whenever every energy level of a
configuration set carried an experimental value, and starred every level
of the `RCEOUT` it then wrote, which discarded those energies if the file
was fed back as `RCEINP`. See item 10 of the additions to the manual
above.
- Gaps in the numbering of linked parameter groups in `RCE` (a group
dropped while a fit is being developed leaves a hole) added a column to
the least-squares system that no parameter feeds. The fitted parameters
were unaffected, but the count of free parameters, and hence the
degrees of freedom and every standard deviation derived from them, were
too large. The groups are now renumbered consecutively before they are
counted.

### RECOMPILING UNDER WINDOWS
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
   appended `COWAN\CODE` in installation step 2 above.
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
for debugging with `gdb`), described in the ADDENDUM below.

### COMPILING FOR OTHER OPERATING SYSTEMS
- To compile the sources for an operating system other than Windows,
first remove the MS-DOS-specific code that `RCG` uses to locate its data
files `SENIOR`, `LSJJ.DAT`, `FOR072`, `FOR073` and `FOR074`. It is easy to
find: it is in the first lines of `RCG`, where those files are opened.
Setting the `MSDOS` flag to zero and editing the hard-coded path is
enough.
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

	Compiling is not enough by itself to run `RCG`. The first `RCG` run
	must use a special input file, supplied with the package, to generate
	the CFP (coefficients of fractional parentage) decks -- the files
	`FOR072`, `FOR073` and `FOR074`. Copy `ING11.CFP` to `ING11`, delete
	the first 8 comment lines from the copy, and run `RCG`. Those comment
	lines state the minimum array dimensions `RCG` needs to handle the
	deck. No other program of the package needs anything beyond the
	change to the timing routine described above.

- The dynamic memory allocation functions `ALLOC`, `DALLOC`, and `RALLOC` in 
`RCG` and `RCE` (introduced in 2026) are machine-specific. See the presently 
commented out versions of the machine-specific lines and uncomment them 
according to the hardware platform used (there are versions for DEC, SUN,
Linux, and VAX).

### ADDENDUM: The list of files in the distribution package.

1) Root directory        - This file and the license.

	Contents:

        ```
        README.md          - this file;
        license.txt        - the terms the package is distributed under.
        ```

2) DOCS directory        - R. D. Cowan's own documentation, copied from
                            his directory at LANL in December 1994 and
                            preserved unaltered.

	Contents:

        ```
        RCN_DOC.txt        - manual for the RCN, RCN2 and HF8 programs;
        RCG_DOC.txt        - manual for the RCG program;
        RCE_DOC.txt        - manual for the RCE program;
        readme.cowan.htm   - R. D. Cowan's write-up for his own version
                             of the package, December 1999.
        ```
 
3) CODE directory   - the executable files, the data files they need,
                       and the utility codes.

	Contents:

	a) Main executables:

        ```
        RCN36K.exe   - Single-configuration radial wavefunctions
        RCN2K.exe    - Single-configuration and configuration-interaction
                       radial integrals and transition matrix elements
        rcg11kd.exe  - Angular matrix elements, energy levels, spectra
                       and cross sections
        RCE20KD.exe  - The least-squares fitting program
        ```
	
	b) Main batch files for running the main executables:
 
        ```
        rcn.bat       - batch file for running the RCN program;
        rcn2.bat      - batch file for running the RCN2 program;
        rcg.bat       - batch file for running the RCG program;
        rce.bat       - batch file for running the RCE program;
        make_cfp.bat  - batch file for creating the binary CFP decks.
        mark_time.bat - called by the four above; times the run.
        ```
 
	c) Data files common to all calculations:

        ```
        FOR072        | Binary CFP decks used by RCG. They are not in
        FOR073        | the distribution package; they are created
        FOR074        | during installation (step 6 above).

        SENIOR        - Text file used by RCG.
        LSJJ.DAT      - Binary deck of Cowan's LS coefficients of
                        fractional parentage together with the jj ones,
                        from which RCG builds the true LS -> jj
                        transformation. It is read only by a run that
                        asks for true jj coupling (KCPLD(8) = 8 in
                        column 38 of the ING11 main card), so a missing
                        file is not an error until then. It is
                        unformatted, and machine-specific with it;
                        FOR\LSJJ\mkdeck.py rebuilds it.
        ```

	d) Compilation batch files and the gfortran runtime libraries:

        ```
        gfortran_O3.bat     - Batch file for compiling a source file with
                              the GNU gfortran compiler. Usage:
                              gfortran_O3 <source_file_name>
                              (The name must be given without the extension,
                              which must be .F)
        gfortran_O3_OMP.bat - Same, with OpenMP enabled. RCG and RCE are
                              threaded; the result is bit-identical to the
                              serial build at any thread count, because no
                              parallel region reassociates a sum.
        gfortran_dbg.bat    - Same as gfortran_O3, but with options for
                              debugging with gdb (-O0 -ggdb3, no OpenMP).
        libgfortran-5.dll   | The gfortran runtime, as distributed by
        libgcc_s_seh-1.dll  | WinLibs 13.2.0 r8. All six are needed:
        libquadmath-0.dll   | RCG and RCE are OpenMP builds and import
        libwinpthread-1.dll | libgomp, which in turn imports libdl.
        libgomp-1.dll       |
        libdl.dll           |
        ```
		(These are the libraries the distributed executables were linked
		against. If you rebuild with a different compiler, replace all six
		with that compiler's own, and do not mix the two sets.)

	e) Auxiliary files (rarely needed):

        `ING11.CFP`     - the complete CFP input deck for `RCG`.

	f) Utility files:

        f.1) Python utilities:

                       Each of the .bat files is accompanied by a
                       corresponding .py Python source file with the
                       same name. You can run the .py files directly by

                       py <file_name.py> [parameters]

                       from any directory -- no PYTHONPATH is needed,
                       because every program puts its own folder on the
                       import path before it imports anything. The .bat
                       file does nothing but call the .py of its own
                       name and hand back its exit code, so `conf_list`
                       and `py C:\COWAN\CODE\conf_list.py` are the same
                       run.

                       Every program prints its own usage when it is
                       given arguments it cannot make sense of, and the
                       docstring at the top of each .py file is the full
                       description of what it does and why.

                       You can modify any of the .py files as needed.
                       However, if you find any bugs or introduce
                       new useful options, please notify me, so that
                       your changes and bug fixes could be
                       incorporated in future releases.

        toolbox\      - The shared library the programs are built on: the
                       readers for IN36, ING11, OUTG11, genealogy.txt and
                       the RCE files, the ASD notation formatter, and the
                       pieces the individual programs assemble. Four of
                       its modules will also print what they have read,
                       which is useful when a deck is not doing what you
                       expect:

                       py -m toolbox.in36 [IN36]
                       py -m toolbox.outg11 [OUTG11]
                       py -m toolbox.genealogy
                       py -m toolbox.jjkeys

                       They must be run with `-m` as shown, from a
                       directory that has COWAN\CODE on PYTHONPATH, or
                       run the launcher .bat of the program you want
                       instead. They are library modules, not commands.

        tests\        - The test suite (`py -m pytest` from COWAN\CODE).
                       A large part of it runs the retired Perl programs
                       from perl_legacy\ side by side with the Python and
                       compares the output, which is the evidence that
                       the port is faithful. Those tests read a finished
                       run and skip themselves when the case in WORK\ has
                       not been run -- the cases ship as inputs only. Run
                       one as its README.txt describes and its tests come
                       back.

        perl_legacy\  - The retired Perl utilities, unchanged. Nothing in
                       the package calls them and they are not on the
                       PATH; they are kept as the reference the Python
                       was accepted against, and because they are the
                       record of what the codes used to do.

        add_eav.bat   - Python program used to add a given energy shift
                       to all Eav parameters in the ING11 file.
                       Usage:

                       add_eav <energy_shift_in_kK>

                       The same action can be achieved by including
                       the rescale 'card' in ING11, i.e. by adding
                       at the top of ING11 a line with zero in
                       position 5 and the value of the energy shift
                       in positions 21 through 30.

        update11.bat  - A batch file to be executed after the LSF fitting
                       with RCE to transfer the fitted parameters into ING11,
                       run RCG with those fitted parameters, and run
                       printout.bat (see below).

        s11.bat       - Python program for substituting the parameters
                       from the RCEOUT file (produced by RCE
                       least-squares fitting program) back into the
                       ING11 file (input for RCG program);
                       Usage:

                       s11

                       (with no parameters - for calculations with
                       one parity set or for substituting the
                       fitted parameters into ING11 for the first
                       parity set only);

                       s11 2

                       (for substituting the fitted parameters into
                       ING11 for two parity sets).

        printout.bat  - Python program to print out the level list in
                       a form similar to LEVELS2, but with
                       configuration and term labels including
                       proper parent terms in the order of summation
                       of shells, in any of the coupling schemes RCG
                       and RCE can produce. The eigenvector
                       compositions with parentage are read from the
                       OUTG11 file, while the experimental energies
                       are substituted from the RCEOUT file. Thus,
                       to produce correct results, you must first
                       substitute the fitted parameters from RCEOUT
                       into ING11 (e.g., by running s11 2) and
                       re-run RCG. The program also produces a file
                       with fitted and HF parameter values and ratios
                       LSF/HF. Usage:

                       printout <levels_out_name> [params_out_name]

                       where <levels_out_name> is the name of the
                       output levels file (typically, levels.txt),
                       and [params_out_file] is the optional name
                       of the parameters output file (if omitted,
                       the default is params.txt). The output files
                       are ASCII tab-delimited files best viewed
                       with Excel or any other spreadsheet program.

                       This program takes into account the possible
                       reordering of subshells (see the description of
                       reorder_ing11.bat below). For it to work, the
                       files IN36, ING11, OUTG11, genealogy.txt, OUTE,
                       RCEINP.HF and LEVELS1 must be present in the
                       current directory. genealogy.txt is written by
                       RCG only when IPGEN, columns 83-84 of the ING11
                       main card, is 2. The scaling factors for the HF
                       parameters are read from OUTG11. Therefore,
                       for the ratios LSF/HF to be correct, the
                       first run of RCE, which produces RCEINP.HF,
                       must be made with the sequence of runs
                       RCN/RCN2/RCG/RCE. The files OUTE and RCEINP.HF
                       are skipped if you specify an optional --no-rce
                       parameter in the command line.

                       The Python version is no longer limited to LS
                       and JJ. Its options are:

                       --couplings LS,JJ
                          Which coupling schemes to print, in this
                          order. Any of the eight RCG and RCE can
                          produce may be named.
                       --from LSJK=JK
                          Take a coupling from a sub-folder, for the
                          common case where the run that has it is a
                          second RCE run in a directory of its own.
                          Merging the two level lists used to be hand
                          work.
                       --from LSr=LS@reordered
                          The same, under a name of your own, so that
                          one levels file can hold two blocks of the
                          same coupling -- which is what a reordered
                          run is.
                       --dominant-config
                          Add a conf0 column: the level's leading
                          configuration as RCE determines it, from the
                          full-precision eigenvector. The weight of a
                          configuration in a level is the same in every
                          coupling scheme, so this is the one
                          configuration name that does not depend on
                          the scheme. The f.2 programs below need it.
                       --drop-full-shells 5p,4d
                          Leave the named subshells out of the printed
                          labels where they are closed.
                       --components 5   --min-percent 1.5
                          How much of each eigenvector composition to
                          print.
                       --match eigenvector|energy-order
                          How a level of a second folder is tied to the
                          same level here.

                       Labels are built from `genealogy.txt`, so the
                       middle of a coupling chain is read rather than
                       guessed, and the principal quantum numbers are
                       taken from IN36 through the ING11 slots, so a
                       reordered deck is labeled correctly. Every
                       label carries the full electron count.

        conv_out.bat  - Python program to convert both level and line
                       lists from the output of RCG and RCE. Input
                       parameters:
                       1) Output levels file name (e.g., 'levels.txt').
                          The format of this particular file is 
                          slightly different from a similar file created 
                          by the printout code. Namely, it contains 
                          several additional columns. The most important 
                          one is the unique level number, which is used 
                          to label the lower and upper levels in the 
                          simultaneously created transitions file.
                       2) Output transitions file name (e.g., tr.txt)
                       3) RCE or noRCE
                          (for including or omitting RCE data such
                          as experimental energy levels)
                       4) LTE or DR or AUGER or CM or
                          BF:<time_of_flight_ns>

                          This specifies how to scale the line
                          intensities. The most popular use is LTE
                          (local thermodynamic equilibrium, i.e. a
                          Boltzmann population distribution). DR stands
                          for 'dielectronic recombination'. This option, as well as
                          AUGER requires autoionization rates to be
                          available in two separate files, one for
                          the first parity set and one for the second.
                          If this option is specified, the names of
                          these files must be given at the end of the
                          parameter list (parameters 7 and 8).

                          CM is a special 'Cascade Matrix' option
                          in calculations of Auger electron spectra.

                          BF stands for 'Beam-foil'. For this option,
                          an additional parameter, time of flight
                          in nanoseconds, must be given together
                          with 'BF:'. It means the time required for
                          the emitting ions to cross the field of
                          view of the spectrometer. The intensity
                          is estimated by integrating the line
                          intensity from each of the decaying levels
                          within this 'time of flight' period.
                       5) Effective excitation temperature in 1000
                          cm-1. This parameter is required with all
                          options, but is used only in LTE.
                       6) Max. scaled intensity
                          A number meaning the ratio of max and min
                          intensities of the lines to be included in
                          the transition list. The min. intensity is
                          assigned the value of 1, and the rest are
                          scaled accordingly. Lines with intensity
                          smaller than 0.5 are omitted. All
                          intensities are rounded to an integer
                          number.
                       7,8) Optional parameters: Autoionization
                          rate files for the first and second parity
                          configuration sets. These files are
                          generated by the utility code read_aa.bat,
                          described below: one run of it per parity,
                          made in the directory holding that parity's
                          OUTG11, with its output redirected to the
                          file named here. That working arrangement --
                          two directories, even configurations in one
                          and odd in the other -- is why there are two
                          files and not one.
                        9) Optional parameter: IDEN
                          It is used to generate input files 
                          ENLEV.DAT and TRANS.DAT for the visual
                          line identification code IDEN2 
                          (V. I. Azarov, A. Kramida, and M. Ya. 
                          Vokhmentsev, Comput. Phys. Commun. 225, 
                          149–153 (2018); 
                          DOI:10.1016/j.cpc.2017.12.012). This 
                          will work only if the observed line list
                          DLV.DAT and an auxiliary file numset.dat
                          are present in the working directory
                          (for explanation of those files, see the
                          above article on the IDEN2 code).
                             If the option IDEN is used, it may be 
                          followed by additional optional parameters 
                          (added in the year 2024 release):
                          truncate
                             and/or
                          linear 
                             or 
                          convert2log.
                             The option 'truncate' can be specified if one
                          wants to reduce the number of energy levels in 
                          the ENLEV file by deleting the levels that do
                          not have any transitions from them.
                             The option 'linear' should be specified if a 
                          linear intensity scale is desired for the 
                          calculated line intensities in the TRANS.DAT
                          file (default is a log scale).
                             The option 'convert2log' can be used to 
                          convert the observed intensities in the DLV.DAT
                          file to a log scale matching the default scale 
                          of calculated intensities. A new DLV.DAT file 
                          will then be created. Options 'linear' and 
                          'convert2log' are incompatible.

                       In the transitions output file, with the LTE,
                       BF, and DR options, the wavelengths are given
                       in angstroms. Those between 2000 and 20000
                       angstroms (or, more precisely, those
                       corresponding to wavenumbers between 5000 and
                       50000 cm-1) are given in standard air,
                       otherwise in vacuum. With the AUGER options,
                       the transition energy in eV is printed
                       instead.

                       For each line, the configuration and
                       LS term label are printed for both levels,
                       as well as the Landé g factors, radiative
                       lifetimes, gA and A values, the value of the
                       cancellation factor, and, if autoionization
                       data are included, the corresponding
                       autoionization rates and branching fractions.

                       The positional parameter list above is long
                       and strictly ordered, and a missing argument
                       shifts everything after it. Every setting may
                       therefore be given as a long option instead --
                       --model, --temperature, --max-intensity, --ai,
                       --couplings, --components, --min-percent,
                       --iden, --truncate, --linear, --convert2log --
                       and the two forms may be mixed freely.

                       Better still, every setting may be written into
                       a `conv_out.par` file in the working directory,
                       in the same layout as `calc_TP_M1E2.par`: the
                       value first, then a semicolon, then the key,
                       then prose saying what it is for. Run conv_out
                       with no par file and no arguments and it writes
                       a fully commented sample and stops. The par file
                       supplies the defaults and the command line
                       overrides them, so the historical positional
                       form above keeps working exactly as it did.

                       This program began as a reimplementation of an
                       older Borland Pascal code of the same name, and
                       a few of that program's options have still not
                       been carried over. What it gains in exchange is
                       that it places no limit on the number of levels
                       or transitions, it is portable, and it can write
                       any of the coupling schemes RCG produces rather
                       than LS and JJ alone.

        read_aa.bat
                     - Python program to process the output of RCG
                       for autoionization calculations (see above).
                       Usage:

                       read_aa <cores> <IP> [--set N] [--folder F]

                       <cores> is a comma-separated list of
                       name:energy giving the core states of the
                       ionized atom with their energies above that
                       ion's ground state, in the run's own energy
                       unit, and with no blanks anywhere in it. A
                       name is a regular expression matched against
                       RCG's channel labels, so one column of the
                       table may sum several of them. <IP> is the
                       ionization energy of the atom being
                       calculated, in the same unit, which is what
                       puts the core energies on this ion's scale.
                       The table goes to standard output; redirect
                       it to the file conv_out is to be given.

                       py C:\COWAN\CODE\read_aa.py explains the rest
                       and gives worked examples.

        reorder_ing11.bat
                     - Python program for reordering the subshells in
                       the ING11 file. This is needed if you want to
                       change the default summation order of the
                       subshells. For example, the configuration
                       1s.2s.2p of Li-like spectra is best described
                       by first combining the 2s and 2p subshells to
                       form intermediate LS terms 3P* or 1P*, and
                       then combining this intermediate term with
                       the 2S term of the 1s subshell. The default
                       order of shell summation in RCG is left to
                       right, i.e., 1s is combined with 2s to
                       produce 3S of 1S, then this is combined with
                       the remaining 2p shell to produce the final
                       4P* or 2P* term. Usage:

                       reorder_ing11 <output_filename> <new_order>

                       where <output_filename> is the name of the
                       new ING11 file (the original ING11 must be
                       present in the current directory; it is used
                       as input), and <new_order> is a sequence of
                       digits, e.g., 231 to specify that shells 2
                       and 3 will go first in the new ING11,
                       followed by shell 3. You can recognize the
                       shells in the top section of ING11 and count
                       their numbers from left to right, starting
                       with 1. The program rearranges not only the
                       shell designations, but also the corresponding
                       Slater and CI parameters.

        scale_param.bat
                     - Python program for scaling the parameter values
                       in the RCEINP file by given factor(s).
                       Parameters:
                       1) input_file_name
                       Name of the input file (normally, a saved copy
                       of the RCEINP file).
                       2) output_file_name
                       Name of the new RCEINP file to produce.
                       3) One or more of the following tokens
                       separated by space:
                       F=<F_factor>
                       G=<G_factor>
                       CI=<CI_factor>
                       GROUP<number>=<factor>
                       where the factors given are used to multiply
                       all parameters from the corresponding group
                       of the input file. The GROUP number is the
                       absolute value of the numerical flag used to
                       group parameters in RCEINP, e.g., 51, 52, etc.
                       4) conf_set
                       (1 for the first parity set, 2 for the second
                       parity set, or 12 for both).
                       5) FIX
                       An optional parameter. If specified, all Slater
                       parameters that were 'free' in the input file
                       will be 'fixed' in the output file for the
                       given parity set(s).

        subst_exp_levs.bat
                     - Python program to substitute experimental levels
                       into RCEINP from another file.

                       Usage:

                       subst_exp_levs <input file> <output_file> [params=flags|noflags|Eav|F|G|Z|CI|none]

                       This program is used to substitute experimental
                       level values from one RCEINP (or RCEOUT) file
                       into another, possibly created with a different
                       set of configurations. Typically, you would
                       need it if you made an LSF fitting with a
                       limited set of configurations and then
                       discovered you need more configurations to
                       include. Then, after you run RCN/RCN2/RCG/RCE
                       with an extended configuration set, instead of
                       manually editing again your new RCEINP file,
                       you can use this program to transfer the level
                       values from your previous calculations to the
                       new RCEINP file, which is the <output_file> in
                       the parameter list.

                       Optionally, you can specify a certain action
                       for the parameters from your old calculation.

                       params=flags will set only the flags for the
                       same parameters in the new RCEINP (i.e., make
                       them fixed or free or grouped in the same way
                       as they were in your old RCEINP).

                       params=noflags will transfer only the values
                       of the parameters but not the flags.

                       params=Eav or params=F or params=G or params=Z
                       or params=CI will transfer only the values and
                       flags for the corresponding parameter group
                       (Z stands for ZETA). Both values and flags
                       will be transferred.

                       params=none will skip all parameter values and
                       flags (leave them unchanged in the new RCEINP
                       file).

        subst_params_ing11.bat
                     - Python program to substitute parameter values
                       from one ING11 file into another one. The
                       substitution is done line by line. If the
                       configuration name(s) on the line read from the
                       input file matches those in the output file, the
                       parameter line(s) for these configurations are
                       substituted in place of similar lines in the
                       output file. The input ING11 file is assumed to
                       have both parities, while the output ING11 is
                       supposed to have only one parity set. The main
                       purpose is to substitute fitted parameters to
                       set up ING11 for autoionization calculations.

                       Usage:

                       subst_params_ing11 <input_ING11> <output_ING11> <conf_set>

                       where <conf_set> can be 1 or 2, for conf. set #1
                       or #2 in the <input_ING11>.

        conf_list.bat
                     - Python program to print out the lists of
                       configurations included in the first and second
                       parity sets of IN36.

        toolbox\vacair.py
                     - The module conv_out uses to convert wavelengths
                       between vacuum and standard air. The functions
                       are Lvac() and Lair(); their arguments are the
                       air and the vacuum wavelength in angstroms,
                       respectively. This was vacair.pl.

        calc_TP_M1E2.bat
                     - Python program to run Monte Carlo trials for
                       estimation of uncertainties of M1 and E2 transition
                       rates (see A. Kramida, Atoms 2, 86–122 (2014),
                       https://doi.org/10.3390/atoms2020086).
                       The settings that used to be hard-coded in the
                       Perl source are now in `calc_TP_M1E2.par`, a
                       commented parameter file. Run the program with
                       no par file present and it writes one into the
                       working directory and stops, so that you can
                       read it through and fill it in; `--write-par`
                       writes a fresh one at any time.

                       The RCG run it reads must print the forbidden
                       lines it is asked about, which needs two things
                       on the IN2 main card. IMAG and IQUAD, columns 49
                       and 50, switch the M1 and E2 spectra on. DMIN,
                       columns 61-65, must be **0.00**: a DMIN between
                       1 and 50 means "omit every transition weaker
                       than gA = 10**DMIN", and forbidden lines are
                       almost all weaker than the 5.00 or 10.0 that
                       suits electric-dipole work, so the spectrum
                       comes out empty. See `WORK\Forbid\README.txt`.

        f.2) From a fitted calculation to a published level list:

                       Three programs with no Perl ancestor, run in
                       this order, each handing the next a file. The
                       middle file is the one a person edits:

                       levels.txt
                         |  config_purity
                         v
                       coupling_choice.txt   <- you mark this one
                         |  make_levels_list_prepub
                         v
                       levels_prepub.txt
                         |  make_levels_list_pub
                         v
                       levels_pub.txt

                       They rest on the invariance noted under
                       --dominant-config above: because recoupling is
                       block-diagonal in the configuration, the weight
                       of a configuration in a level is the same in
                       every scheme. Purities of two schemes are
                       therefore comparable, and a composition may
                       quote each configuration in its own scheme
                       without the percentages ceasing to add up.

                       Prepare the input with

                       printout levels.txt --dominant-config
                                --couplings LS,JJ,LSJK --from LSJK=JK
                                --components 10 --min-percent 0.5

                       Parity is the first column of all three
                       outputs: deleting a column for an article table
                       is easy, reconstructing one is not.

        config_purity.bat
                     - Which coupling scheme describes each
                       configuration best. Writes coupling_choice.txt:
                       one row per configuration in IN36 order, one
                       10-character column per coupling scheme found in
                       levels.txt, and in each cell the purity -- the
                       unweighted mean, over the levels of that
                       configuration, of the leading component's
                       percentage. A star marks the chosen scheme, one
                       per row, and the file is written with the star
                       already on the highest purity. That is a
                       suggestion, not a verdict: previous literature
                       on the spectrum, and the schemes used for
                       neighboring isoelectronic ions in the database,
                       routinely outweigh a few percent gained by a
                       rarely used scheme. Move the stars by hand.

                       The file is fixed-width because it is meant to
                       be edited. It is never silently discarded:
                       writing over an existing one needs --force, and
                       --merge updates the numbers while leaving the
                       stars where you put them.

                       Options: --levels, -o, --force, --merge.

        make_levels_list_prepub.bat
                     - One coupling scheme per configuration, applied.
                       For each level it takes from each scheme's block
                       only those components whose configuration was
                       assigned that scheme, concatenates them across
                       schemes, sorts by percentage and prints. A level
                       mixing two configurations may therefore read
                       44% in JJ for one and 35% in LS for the other.
                       Writes levels_prepub.txt.

                       Options: --levels, --choice, -o, --components,
                       --min-percent.

        make_levels_list_pub.bat
                     - Designations unique within each parity and J,
                       which is what a published or database level
                       needs. The name comes from the leading component
                       where the level is reasonably pure; where it is
                       not, from further down the composition; failing
                       that, from letter prefixes (a, b, ...) assigned
                       per configuration and term across all J, in
                       order of increasing energy. Where no meaningful
                       name can be assigned at all, the name is the
                       integer part of the energy. Writes
                       levels_pub.txt, in cm-1.

                       --nlev takes one number per configuration set,
                       in IN36 order, saying how many levels to
                       publish; the default is the position of the
                       highest observed level in each set, and a set
                       with no observed level makes --nlev mandatory.

                       Options: --nlev, --n_comp, -o,
                       --min-name-percent, --letters-max,
                       --leading-min, --max-outside.

4) FOR directory    - Fortran source files of the Cowan programs.
                       Fixed-form Fortran 77, built with -std=legacy.

        RCN36K.F
        RCN2K.F
        rcg11kd.f
        rcgpar.f      - Include file for rcg11kd.f (array dimensions)
        RCE20KD.F
        RCE20KD_PAR.F - Include file for RCE20KD.F (array dimensions)
        LSJJCOM.F     | Include files for the true jj coupling
        LSJJID.F      | added to RCG and RCE in 2026
        LSJJPAR.F     |
        RCEJJ2E.F     - Include file for RCE20KD.F, true jj coupling
        LSJJ\         - The Python that builds CODE\LSJJ.DAT
                        (mkdeck.py) and the scripts that validated the
                        coefficients against Cowan's own. Its README.txt
                        is the long-form account of that work.


5) WORK directory   - Sample directory for Cowan calculation projects,
                       and the examples that come with the package. The
                       two files in WORK itself are the installation
                       test; each sub-directory is a worked example of
                       one thing the codes can do.

        IN36          - sample input file for RCN (for two configurations
                       of Sn VIII);
        IN2           - sample input file for RCN2;

        AI\           - B IV example of an autoionization calculation:
                       IN36 and IN2 only, run the four programs in an
                       empty directory holding a copy of them.

        AI\DIEL\      - the same physics set up by RCN's DIEL option,
                       which writes the continuum configuration cards
                       for you instead of your punching them by hand.
                       Inputs (indiel, in36, in2, gcard) and reference
                       outputs to compare your run against. See its
                       README.txt.

        Forbid\       - Sn VIII 4d7, forbidden lines: an RCG run with
                       the M1 and E2 spectra switched on (IMAG and
                       IQUAD, columns 49-50 of the IN2 main card), and
                       the example for calc_TP_M1E2. IN36 and IN2 only.
                       Its README.txt explains why DMIN, columns 61-65,
                       must be 0.00 for forbidden transitions -- with
                       the usual electric-dipole cut of 5.00 or more,
                       RCG prints the section headings and no lines
                       under them.

        Pr3\          - Pr III, a real least-squares fit, and the worked
                       example of everything downstream of one: the
                       RCEOUT/RCEINP input format, the update11 cycle,
                       the merge of a coupling computed in a second
                       folder (Pr3\JK, the same fit asked for in
                       LS + LSJK instead of LS + JJ), and the three
                       programs of section f.2. It is also the conv_out
                       example, and carries the dlv.dat and numset.dat
                       an IDEN run needs. Inputs only -- the run writes
                       about 412 MB. Its README.txt gives the command
                       sequence step by step, with the figures each
                       step should report.

        W70\          - W LXX, boron-like tungsten: the true jj example
                       (KCPLD(8) = 8, column 38 of the ING11 main
                       card), and with it the true-jj_keys.txt key to
                       the 8-character basis labels and the LEVELS*J
                       files. It is also the example of a run with
                       nothing observed, so there is no fit, every
                       parameter is fixed, and make_levels_list_pub
                       has to be told --nlev. IN2 and IN36 only; the
                       run takes seconds. See its README.txt.

        U90\          - U XC, lithium-like uranium: a small, fast
                       least-squares fit over the K-shell-excited
                       configurations, and the reorder_ing11 example.
                       U90\reordered is the same fit with the ING11
                       subshell slots permuted, and shows what a
                       different coupling order does to the level
                       designations and to the percentages -- the same
                       eigenstates and the same energies, expanded in
                       a different basis. Inputs only. See its
                       README.txt.


                                  August 19, 2026

                                  Alexander Kramida,
                                  e-mail akramida @ umd.edu
