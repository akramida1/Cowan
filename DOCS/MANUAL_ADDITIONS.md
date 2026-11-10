# Additions to Cowan's manuals

[Back to the README](../README.md)

R. D. Cowan's own manuals are in [`Cowan/`](Cowan/):
[RCN_DOC.txt](Cowan/RCN_DOC.txt) for RCN, RCN2 and HF8,
[RCG_DOC.txt](Cowan/RCG_DOC.txt) and [RCE_DOC.txt](Cowan/RCE_DOC.txt).
They are preserved unaltered, so they do not describe what has changed
since. This document does: how the programs are run now, their input and
output files, and the options added to RCG and RCE.

## Running the programs
Programs are normally run in the following order:
`rcn`, `rcn2`, `rcg`, `rce`.

These are the names of the wrappers -- the `.bat` files in `CODE/` on
Windows, the shell scripts in `CODE/unix/` on Linux and macOS. Type them
at the command prompt from the directory holding the input files. The
commands are spelled the same way on every platform, in lower case, here
and in every README of the package; Windows is case-insensitive, so `RCN`
works there as well. Each program produces the input for the next, and
RCG and RCE are often run several times after the first pass.

After a least-squares fit (LSF) with RCE, one normally wants to
recalculate the transition rates by running RCG with the fitted
parameters from the LSF output. The easiest way to do that is to run the
utility code:

`update11`

Nothing has to be typed after it: update11 reads the `in36` file, counts
the parities in it, and substitutes both parity sets when both are there.
(The older form `update11 2` still works and still means both.) Besides
transferring the fitted Slater parameters from the `RCE` output (RCEOUT,
to be exact) into the `ing11` file and re-running `RCG` with these
parameters, this does some additional useful things (see `update11` in
[UTILITIES.md](UTILITIES.md)).


## Additions to the manual

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
for a forbidden-transition calculation** -- see `WORK/Forbid/README.txt`.

	**In an `RCG` transition list, *gA* is the primary reference quantity,
	not log(*gf*).** If you need *S*, *f* or *A*, take them from *gA* and the
	transition energy. The reason is the printed precision, not the physics:
	log(*gf*) is written with three decimals, and half a unit in the last
	decimal of a *base-10 logarithm* is ln(10) x 0.0005 = 1.15x10^-3^, i.e.
	**0.115 % in *gf*** -- already visible in the second significant decimal
	when you convert back. Recovering a line strength from log(*gf*)
	therefore injects about 0.1 % of random-sign scatter that is simply not
	present in *gA*.

	The natural expectation is the opposite, that *gf* should be the better
	of the two, since it comes from *S* with a lower power of the transition
	energy (sigma rather than sigma^3^), so that rounding of the energy
	propagates into it less strongly. That is true as far as it goes, and it
	is beside the point: the dominant error here is not in the energy, it is
	contributed by the log(10) taken at print time, and that step is the one
	routinely overlooked. *gA* is printed in floating form and carries its
	own digits; log(*gf*) does not. This has been checked more than once and
	the conclusion has always been the same.

4. The designations of terms have changed slightly. In the original
version of `RCG` two different terms could end up with the same name,
even with the additional (e.g., seniority) quantum numbers in use. To
remove that ambiguity I added a distinguishing letter ('a', 'b', ...)
to such terms. For example, the term names `(3P) 3Da` and `(3P) 3Db` mean
different ^3^D terms: their parent terms have the same name but different
genealogy. For such cases, the user may have to manually trace the
history of the summation of quantum numbers that produced these terms,
using some additional printout in `OUTG11`. This task is automated by the
utility codes `printout` and `conv_out` (see [UTILITIES.md](UTILITIES.md)), which read the
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

8. Each configuration card of `ING11` now has its own **term cutoff** in
columns 83-85, limiting the number of parent *LS* terms of that
configuration's **first subshell**. It is a per-configuration version of
`NLT11` (columns 26-27 of the main card), which is one number for the
whole run and exempts the last configuration of each parity, and which
therefore cannot serve two configurations whose cores differ -- Dy VIII
with f^6^, f^7^ and f^8^ cores is the case it was written for. **A blank
field means "use `NLT11`"**, which is what every existing deck has, so
nothing changes for them; the field is applied with no exemption when it
is punched. Cowan's card ends at column 80 and `RCN2` still writes it
that way, so old records simply blank-fill -- but note that `RCN2` will
overwrite a hand-punched field if `ING11` is regenerated, the same
exposure hand-edited Slater parameters already have.

	Only **f shells** are accepted, and only a **short list of values** for
	each. The cut keeps the first *n* terms of the CFP deck, and deck order
	is neither energy order nor consistent from shell to shell. Where the
	deck lists the terms by descending multiplicity, a cut that ends on a
	multiplicity-block boundary keeps whole blocks of the highest
	multiplicities -- by Hund's rule the parents lying lowest, which is the
	truncation the physics justifies. Hence:

	| Shell | Accepted values |
	|---|---|
	| f^6^, f^8^ | 1, 17, 73, 119 |
	| f^7^ | 1, 7, 47, 119 |
	| f^9^ | 3, 27, 73 |
	| f^10^ | 5, 27, 47 |
	| f^2^, f^3^, f^4^, f^5^, f^11^, f^12^ | the full term count only (7, 17, 47, 73, 17, 7), i.e. no cut |

	The last value of each row is the full count and means "no cut". For
	f^4^, f^5^ and f^11^ the deck runs the other way and for f^2^, f^3^ and
	f^12^ it alternates, so no leading part of them is a set of low terms:
	a prefix of f^4^ would discard the ^5^S..^5^I block, the ground terms.
	Anything else stops the run before any work is done, and the message
	names the parity set, the configuration and the values that shell
	allows.

	A cutoff and **true *jj* coupling** cannot both be had: the rotation
	needs the complete term list of the shell, and with terms missing the
	term map is not a bijection. When both are asked for, the cutoff wins
	-- true *jj* is switched off and a warning is written to `OUTG11` and
	to the screen. (The earlier line announcing that true *jj* will be
	added is printed before the configuration cards are read; the later
	message retracts it.) Cowan's *JJ* (*J*~1~*J*~2~) output is unaffected.

	The saving is what makes such cases possible at all: Dy VIII
	(f^7^p^4^, f^8^p^3^, f^6^p^5^ / f^7^p^3^d, f^6^p^4^d, f^7^p^3^s,
	f^6^p^4^s) on 8 threads spends about 2 hours in `CALCFC` alone with no
	cutoff; with cutoffs of 17 / 7 / 17 the whole run takes 3 min 35 s, of
	which `CALCFC` is 11.7 s.

9. The parameter `ICRIT2` is disabled in `RCE`: its value has no
effect. The program exits after completing the fit, whether the
parameters were linked or unlinked, and never starts another run with
the next value of `LSQM`.

10. *Suspended* execution is not supported in `RCE`: intermediate
results cannot be saved after the *n*-th iteration and resumed from
there. The new `RCEOUT`/`RCEINP` files (see below) achieve almost the
same thing; the only cost is one extra initial iteration when a new run
is started from an `RCEOUT` copied into `RCEINP`. The `I216` parameter is
disabled in `RCE`, and its position in `OUTGINE` is taken by `IW6`, which
directs some output to the screen when it is negative.

11. Earlier versions of `RCE` refused to fit when every energy level
was flagged "experimental": the run printed `Iteration No 0` and stopped,
and, worse, wrote `RCEOUT` with every level starred, so feeding it back as
`RCEINP` also discarded the experimental energies. The 2026 revision
removed that restriction. A case in which every level is known is the
best-determined one there is, and `RCE` now fits it like any other,
whatever iteration count `OUTGINE` carries. The only conditions that still
end a run after the zeroth cycle are the two that leave nothing to fit: no
varied parameter at all, or no known level anywhere.

12. `RCE` accepts an alternative input file format. On its first run it
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
   of the `LEVELS1` file in *LS* coupling: the energy; the status
   character -- an asterisk (`*`) for an "unknown" level, a blank for an
   ordinary "experimental" one, or a letter naming a relative-energy
   group (see the next section); the
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

5. If any level record carries a group letter, the card block giving the
   two standard deviations of each relative-energy group follows the
   parameter cards. See the next section.

The same structure (1-5) follows for the second parity, if there is one.


## Relative-energy groups in an RCE fit

An experimental level list is not one measurement. A compilation carries
systems of levels known to a hundredth of a wavenumber inside themselves
and tied to the ground level only through one resonance transition
measured to a few hundred -- the `+x` systems. Before this release `RCE`
had two things to do with such a level, and both are wrong. Fitted at
unit weight (a blank in column 14), its unknown offset is charged to the
parameters, and since some parameters are shared between configurations
-- ZETA of an open d shell is one linked group across every odd
configuration of a Ni-like ion -- the error does not stay where it was
made. Marked unknown (`*`), the accurately spaced levels inside the
system are thrown away with it.

A **relative-energy group** lets you say what is actually known: the
levels are precise relative to one another, and the whole system sits at
an uncertain absolute energy.

**Marking the levels.** Column 14 of a level record of `RCEINP` -- the
column that held the `*` -- now takes a letter:

	blank   an ordinary experimental level, fitted at unit weight
	*       not fitted
	A-Z     the base of that group -- the one level whose tie to the
	        ground level was measured
	a-z     another member of the same group

A level cannot carry both a letter and a `*`, the column being one
character wide, so every lettered level is a fitted one.

**The group cards.** The two standard deviations that say what a group is
worth have nowhere to live on a level record, so they come on a card
block at the end of the deck, after the parameter cards, one block per
parity section:

	 GROUP SIGMA_BASE  SIGMA_INT
	X         500.000     20.000
	Z                     15.000

The header is any line containing `GROUP`. One card per group follows:
the letter in column 1, `SIGMA_BASE` in columns 2-13, `SIGMA_INT` in
columns 14-25. A blank line, a `-1` line, or end of file ends the block.
**The sigmas are in the energy unit of the deck** -- `UENRGY` from
`ING11`, which is 1000 cm^-1^ in most decks, so 500 cm^-1^ is punched as
`0.500`. A blank `SIGMA_BASE` means the group has no measured tie to the
ground level at all, and its base then contributes no residual; a blank
`SIGMA_INT` means weight 1. Zero is not accepted for `SIGMA_BASE`,
because zero reads as "known exactly", which is the opposite of what a
blank field means.

`RCE` writes the letters and the card block back to `RCEOUT`, so the
usual `RCEOUT` -> `RCEINP` cycle carries the grouping forward with
everything else.

**What it does to the fit.** A member of a group contributes the
difference `(E_exp,i - E_exp,base) - (E_calc,i - E_calc,base)` instead of
its own residual, and the base contributes its absolute residual, so a
group of *n* levels still gives *n* rows. Each row is weighted by

	w = sigma_0^2^ / (sigma_0^2^ + sigma^2^)

with `sigma_int` for the differences and `sigma_base` for the base.
`sigma_0` is not supplied: it is the standard deviation the fit itself
achieves, taken from the previous cycle, so the weights are fixed
self-consistently rather than by a magic number. It never grows, and it
is frozen once it settles to within one per cent. `OUTE` prints the raw
`STDEV`, the adopted `sigma_0`, whether it is frozen, and the two weights
of every group, every cycle.

**Two things to know before using it.**

- `AVDEV` and `STDEV` of a grouped fit are weighted quantities over
differenced rows, and are **not** comparable with those of the same
deck without the letters. Judge a grouped fit by the O-C of the levels
you care about, not by the number at the foot of the page.
- `SIGMA_INT` is **not** the experimental uncertainty. It is weighed
against `sigma_0`, and `sigma_0` comes from the same fit those levels
are in, so a group bad enough to need downweighting also inflates the
quantity it is being weighed against. Set `SIGMA_INT` to ten times or
more the `AVDEV` you expect from the rest of the fit. `WORK/Mo15/`
is the worked example of exactly this, and its `README.txt` shows what
happens at both settings.

A deck with no letter in column 14 is read, fitted and printed exactly as
before, to the byte. The details are in `FOR/LSJJ/README.txt`.


## The calculated ground level at zero: `rce --ground`

A least-squares fit is easier to start with the calculated ground level
at zero, where the observed one is. On the first RCE run of a fit, type

        rce --ground

instead of `rce`. Before RCE starts, `add_eav --ground` takes the lowest
level listed in `OUTGINE` as the ground, adds minus its energy to every
Eav and every listed level, and keeps RCG's file as `OUTGINE.RCG`. It
reports the level it took, for example

        Ground level adopted: level 1 with J = 0 of configuration set 1, E_RCG = -14.4602

The *J* is taken from `genealogy.txt`; without it the level is named by
its position, "level 1 of J block 1". This replaces the old first step
of a fit -- `add_eav` on `ING11` and a second run of RCG -- because RCE's
first run takes its starting parameters and its calculated levels from
`OUTGINE` alone: `TAPE2E` holds neither, and a shift common to every Eav
moves no eigenvector.

- It acts only when `INVECT` in the first record of `OUTGINE` is 0. With
  any other value RCE reads `RCEINP` instead; `add_eav` then warns,
  shifts nothing, and RCE runs as usual.
- An `OUTGINE` whose levels are all zero -- what RCG writes when it made
  no diagonalization -- is an error, and RCE is not run.
- A ground already at zero is left alone, so a second `--ground` does
  nothing.
- An Eav at or below -4000 is Cowan's mark of a continuum configuration
  in an autoionization run, not an energy (see the autoionization
  section of [RCG_DOC.txt](Cowan/RCG_DOC.txt)); it is not moved, and
  neither are the continuum levels of its configuration set. (RCG's own
  treatment of such levels is described
  [below](#levels-at-or-below--4000-no-transitions).)
- `ING11` and `OUTG11` keep RCG's values, so `printout` run straight
  after this first RCE run shows RCG's Eav differing from those in
  `RCEOUT`, and `RCEINP.HF` holds the shifted Eav. The first `update11`
  carries the fitted values into `ING11`, and from then on they agree.

`add_eav --ground [folder]` can also be run by itself. The full account,
with the acceptance runs, is in `FOR/LSJJ/README.txt` under THE
CALCULATED GROUND LEVEL AT ZERO FOR RCE.

## Levels at or below -4000: no transitions

RCG takes any level whose energy is at or below -4000, in the deck's
energy unit (1000 cm-1 unless the `ING11` main card says otherwise, so
-4 000 000 cm-1), for a continuum state. That is Cowan's convention for
autoionization runs, in which the continuum configurations carry an Eav
of -9500 or -8500 (see the autoionization section of
[RCG_DOC.txt](Cowan/RCG_DOC.txt)). `SPECTR` computes no radiative
transition to or from such a level, and it applies the test to every
deck, whether or not it has a continuum configuration.

RCN2 puts the Eav of the first configuration at zero, so in practice
this touches only highly charged ions, whose ground level can lie
thousands of kK below that average. `WORK/W70` (W69+) is one: its ground
level, 2s2.2p 2P1/2, comes out at -8391.10 kK. It appears in the
eigenvalue print of `OUTG11` and in none of its line lists -- every
transition to the ground level is missing.

The remedy is to put the ground level near zero and run RCG again,
which is what a sensible output needs anyway: shift `ING11` after RCN2
with

    add_eav 8391.1026

(minus the ground energy RCG printed, in kK), then run `rcg`. In
`WORK/W70` the ground level then has its 36 line entries, and every
energy downstream is measured from it, to about 1 cm-1 (`ING11` holds an
Eav above 99999.99 only to 0.1). `rce --ground` shifts only
`OUTGINE`, so it does not bring the lines back by itself; the first
`update11` of a fit, which carries the shifted Eav into `ING11` and runs
RCG again, does.

