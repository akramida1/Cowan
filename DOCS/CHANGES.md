# Release notes

[Back to the README](../README.md)

## Version 26.2.3dev -- in development

Changes since version 26.2.2 of November 10, 2026, not yet released. The
next release is planned as 26.3.

- **Documentation corrected: `update_IDEN`.** Since 26.2.1 it has given
  each calculated level the root-mean-square of E_fixed - E_calc over
  the fixed levels of its own configuration, and applied `--ucalc` only
  to the configurations with none; `--nosd` switches that off. The
  documentation still described the older rule, under which `--ucalc`
  went to every calculated level, and did not mention `--nosd`. See
  [UTILITIES.md](UTILITIES.md#update_iden). The program is unchanged.
- **Windows executables rebuilt** to carry 26.2.3dev in their version
  resource, from sources unchanged since 26.2.1. On the Pr III
  regression case every result file is byte-identical to that of the
  26.2.2 executables.

## Version 26.2.2 -- November 10, 2026

Changes since version 26.2.1 of August 21, 2026:

- **New: `rce --ground`.** On the first RCE run of a fit it puts the
  calculated ground level at zero by rewriting `OUTGINE`, which used to
  take `add_eav` on `ING11` and a second run of RCG. The same is
  available as `add_eav --ground`. See
  [MANUAL_ADDITIONS.md](MANUAL_ADDITIONS.md#the-calculated-ground-level-at-zero-rce---ground).
- **Fixed:** the Python reader of `ING11` took the `-55555555.` card for
  a continuation of the last parameter card when nothing stood between
  them -- a deck of one configuration with no transitions -- so `add_eav`
  and the other utilities that read `ING11` failed on such a deck.
  Continuation cards are now read by count, as RCG reads them.
- **Fixed: the Windows `RCE20KD.exe` of 26.2.1 predated its own source.**
  It had been built before the relative-energy groups were added to
  `RCE20KD.F`, so on Windows RCE stopped with "Bad value during integer
  read" on any `RCEINP` carrying a group letter -- `WORK/Mo15` among
  them. Linux and macOS, which compile from the sources, were not
  affected.
- **All four Windows executables rebuilt** from the unchanged sources,
  with the same compiler (WinLibs gfortran 13.2.0 r8), and they carry
  26.2.2 in their version resource. On the Pr III regression case every
  result file -- `TAPE2N`, `ING11`, `TAPE2E`, `OUTGINE`, `genealogy.txt`
  and all of RCE's -- is byte-identical to that of the 26.2.1
  executables, and the listings differ only in their timing lines; on
  `WORK/Mo15` RCE now fits, to the AVDEV its README gives.
- **How to tell which version you have:** see the
  [README](../README.md). `install.sh` now installs `version.txt` on
  Linux and macOS, where it is the only record of the version.
- **Documentation:** `README.md` is now a front page; the documentation
  it held is in `DOCS/`, and R. D. Cowan's own manuals are in
  `DOCS/Cowan/`. `WORK/Mo15` now starts its fit with `rce --ground`.
- **Documented, not changed:** RCG computes no transitions for a level at
  or below -4000 in the deck's energy unit, Cowan's mark of a continuum
  state, even in a deck with no continuum configuration. In a highly
  charged ion that can be the ground level (`WORK/W70`); shifting `ING11`
  so the ground is near zero and running RCG again restores its lines.
  See
  [MANUAL_ADDITIONS.md](MANUAL_ADDITIONS.md#levels-at-or-below--4000-no-transitions).

The Fortran sources are unchanged from 26.2.1.

## Program versions

The date of the last change to each code, as of version 26.2.3dev.

| Code                  | Date     |
|-----------------------|----------|
| RCN  ver. 36k         | Aug 2026 |
| RCN2 ver. 2k          | Aug 2026 |
| RCG  ver. 11kd        | Aug 2026 |
| RCE  ver. 20kd        | Aug 2026 |
| **Python utility codes:** |        |
| printout.py                 | Aug 2026 |
| conv_out.py                 | Aug 2026 |
| add_eav.py                  | Nov 2026 |
| s11.py                      | Aug 2026 |
| subst_params_ing11.py       | Aug 2026 |
| subst_exp_levs.py           | Aug 2026 |
| insert_level_ids.py         | Aug 2026 |
| write_IDEN_level_ids.py     | Sep 2026 |
| update_IDEN.py              | Sep 2026 |
| reorder_ing11.py            | Aug 2026 |
| scale_param.py              | Aug 2026 |
| conf_list.py                | Aug 2026 |
| read_aa.py                  | Aug 2026 |
| calc_TP.py                  | Aug 2026 |
| extract_Icalc.py            | Aug 2026 |
| config_purity.py            | Aug 2026 |
| make_levels_list_prepub.py  | Aug 2026 |
| make_levels_list_pub.py     | Aug 2026 |
| toolbox/ (shared modules)   | Nov 2026 |
| **Command wrappers:**       |          |
| CODE/*.bat (Windows)        | Nov 2026 |
| CODE/unix/ (Linux, macOS)   | Nov 2026 |

Each of the Python codes above, except `config_purity`, the two
`make_levels_list` programs, the two IDEN2 programs and `toolbox/`,
replaces a Perl program of the same name, kept in `CODE/perl_legacy/`.
Two Perl files have no Python namesake: `conv_cowan.pl` was
`conv_out.pl`'s own library and moved with it, and `vacair.pl` is now
`toolbox/vacair.py`. The `config_purity` and two `make_levels_list`
programs, `write_IDEN_level_ids` and `update_IDEN`, and the `toolbox`
library itself, are new: they never had a Perl version.


## The fundamental constants

**Read this before comparing output from this release against output
from an earlier one.** The physical constants the programs use have been
brought up to CODATA 2022, and a small number of printed quantities move
as a result. Nothing here is a change of method: every formula in the
package was checked algebraically first and every one of them was
already right.

Until this release the constants were written as literals scattered
through the Fortran, and most of them were the values of the late 1960s.
Rinf appeared as `109737.31` in one routine, `109737.0` in another and
`109735.3` in a third; alpha as `137.036` in one place and `137.04` in
another; the Bohr radius as `0.529`; pi as `3.14159`. They are now
**named parameters in one include file, `FOR/fundconst.f`**, which every
routine that needs one includes. The rule that file states, and that the
test suite enforces with no exception list, is that *every* fundamental
constant used anywhere in the package lives there -- there is no
threshold of smallness below which one may stay written inline, because
a literal left in the Fortran is one that the next update of these values
cannot reach.

`fundconst.f` has two parts. The first holds CODATA values exactly as
NIST publishes them, each with its relative uncertainty in the comment
beside it. The second holds everything derived from them -- unit-scaled
variants, products, quotients -- *computed* from the first rather than
typed in again, so that two spellings of one constant cannot drift
apart. Each derived entry also carries, in a comment, the literal that
used to stand in the source, so the old numbers can be restored if you
ever need to reproduce a historical result. (One exception: `RINF` now
stands for three different legacy roundings and its comment can restore
only one of them.) The Python utility codes take their coefficients from
the same place, through `CODE/toolbox/fundconst.py`, which a test holds
to the Fortran file bit for bit.

**What moves.** Relative shifts, on the same input:

| Quantity | Shift | Note |
|---|---|---|
| Autoionization rates (`AA`, `SUMGAA`, the dielectronic-recombination table) | **+4.8x10^-4^** | a truncated constant, see below |
| *gf* (E1) in `OUTG11` | +1.9x10^-5^ | |
| *gA* (E1) in `OUTG11` | +8.6x10^-5^ | |
| M1/E1 ratio | +5.8x10^-5^ | |
| E2/E1 ratio | +6.7x10^-4^ | the old *a*~0~ was 0.529 and enters squared |
| Lande *g*-factors | -2.3x10^-7^ | below the printed field width almost everywhere |
| M1 matrix elements | -6.9x10^-7^ | |
| Maxwellian rate coefficient | -1.5x10^-5^ | |
| Photoionization cross-section coefficient | +8.2x10^-5^ | |

Nothing structural changed anywhere: no line count, no label, no column
and no field width moved, in any of the five regression cases the work
was verified against.

The largest of these deserves a sentence of its own, because it is the
one you will notice. **Autoionization rates printed by every earlier
version of `RCG` were about 0.05 % low.** The rate expression in
`SUBROUTINE ENERGY` is Fermi's golden rule, and its conversion factor
Ry/hbar was written as `2066` where the value is 2067.0686667586 -- a
truncation rather than a rounding. The rates rise by 4.81x10^-4^ rather
than the 5.2x10^-4^ that implies, because the Rydberg constant in the
same expression is squared in the denominator and its own correction
pulls the other way. This affects `WORK/AI` and `WORK/AI/DIEL` and any
calculation like them.

Two of the constants the utility codes carried were not stale but
wrong -- one E2 coefficient inconsistent with the E2 coefficient beside
it, and an eV conversion out by a factor of a thousand. Both are now
right and both are described under the summary of bugs fixed below,
together with three more the Fortran carried. The *gA*/log(*gf*)
scatter that is sometimes reported is **not** among them -- it is
quantization of the three-decimal log(*gf*) print field and there is
nothing in the arithmetic to correct. See item 3 of the
[additions to the manual](MANUAL_ADDITIONS.md#additions-to-the-manual), which explains why *gA*, not log(*gf*), is the quantity to
take numbers from.

The full account, with every measurement, every verification step and
the three sites that no test case in `WORK/` exercises, is in
`FOR/LSJJ/README.txt` under EVERY FUNDAMENTAL CONSTANT IN ONE FILE.


## Major bugs fixed
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
was fed back as `RCEINP`. See item 11 of the
[additions to the manual](MANUAL_ADDITIONS.md#additions-to-the-manual).
- Gaps in the numbering of linked parameter groups in `RCE` (a group
dropped while a fit is being developed leaves a hole) added a column to
the least-squares system that no parameter feeds. The fitted parameters
were unaffected, but the count of free parameters, and hence the
degrees of freedom and every standard deviation derived from them, were
too large. The groups are now renumbered consecutively before they are
counted.
- The conversion factor Ry/hbar in the autoionization rate computed by
`RCG` (`SUBROUTINE ENERGY`) was written as 2066 where the value is
2067.0686667586 -- a truncation, not a rounding. Together with the two
other stale constants on the same line this made every autoionization
rate `RCG` has ever printed about 0.05 % low. All three now come from
`FOR/fundconst.f`; see [The fundamental constants](#the-fundamental-constants) above.
- The line-selection threshold in the `TEXCIT` branch of `RCG`'s
`SUBROUTINE SPECTR` carried 0.6698 where the seven other sites of the
same constant carry 0.66698, a dropped digit worth +0.42 %. The value is
not printed, but it decides which lines survive the cutoff, so the
affected transition lists were missing -- or keeping -- lines they should
not have.
- `RCN` used two different values of the fine-structure constant: the
spin-orbit parameter took alpha^2^/2 from a pre-1969 alpha, 4x10^-5^ low,
while the mass-velocity and Darwin terms took (alpha^2^/2)Rinf evaluated
with the modern one. Half of it had been updated at some point and half
had not. Both now derive from the same constant.
- The E2 oscillator-strength coefficient used by `conv_out` was
internally inconsistent with the E2 transition-probability coefficient
beside it in the same table, by +3.6x10^-4^; the two are related exactly
and cannot both be right. The E2 rows of the absorption model were
therefore out by that much relative to the E1 and M1 rows. Both
coefficients now come from `CODE/toolbox/fundconst.py`, and so from the
same values the Fortran uses.
- The Auger *dE* column of `conv_out` was printing meV under an eV
heading: the eV-to-cm^-1^ conversion constant held the value per 1000
cm^-1^ while the code used it as the value per cm^-1^. The value was
wrong and the usage was right, in the Python and in the Perl it was
translated from. The column now reads in eV, as its heading always
claimed. This is a deliberate difference from `CODE/perl_legacy/`, which
still has the old behaviour.

