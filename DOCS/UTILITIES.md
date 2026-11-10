# Utility codes

[Back to the README](../README.md)

The utility codes prepare input for the four main programs and interpret
their output. They are written in Python 3; the Perl versions they
replaced are kept, unchanged and unused, in `CODE/perl_legacy/`, as the
reference the Python was accepted against. Each is run by its command
name from the working directory, on every system.

| Command | What it does |
|---|---|
| [`add_eav`](#add_eav) | shift every Eav of `ING11`; with `--ground`, put the calculated ground of `OUTGINE` at zero |
| [`update11`](#update11) | after an RCE fit: fitted parameters into `ING11`, RCG again, `printout` |
| [`s11`](#s11) | substitute the parameters of `RCEOUT` into `ING11` |
| [`scale_param`](#scale_param) | scale parameter values in `RCEINP` |
| [`subst_exp_levs`](#subst_exp_levs) | carry experimental levels into an RCE deck |
| [`subst_params_ing11`](#subst_params_ing11) | copy parameter cards from one `ING11` to another |
| [`reorder_ing11`](#reorder_ing11) | reorder the subshells of `ING11`, to change the order in which they are coupled |
| [`printout`](#printout) | the level list, with compositions and parameters |
| [`conv_out`](#conv_out) | level and line lists from the output of RCG and RCE |
| [`read_aa`](#read_aa) | autoionization output of RCG |
| [`conf_list`](#conf_list) | lists of configurations |
| [`insert_level_ids`](#insert_level_ids) | unique level identifiers in the `conv_out` tables |
| [`write_IDEN_level_ids`](#write_iden_level_ids), [`update_IDEN`](#update_iden) | the round trip through IDEN2 and LOPT |
| [`calc_TP`](#calc_tp) | Monte Carlo uncertainties of transition probabilities |
| [`extract_Icalc`](#extract_icalc) | the columns of a `calc_TP` line list a line-classification code needs, with synthetic intensities |
| [`config_purity`](#config_purity), [`make_levels_list_prepub`](#make_levels_list_prepub), [`make_levels_list_pub`](#make_levels_list_pub) | [from a fitted calculation to a published level list](#from-a-fitted-calculation-to-a-published-level-list) |

The docstring at the top of each `.py` file is the complete description
of that program. The supporting code -- the shared library `toolbox/`,
the `LSJJ.DAT` generator `lsjj/`, the test suite and `perl_legacy/` --
is described [at the end](#supporting-code).

## Running the utilities

Each utility is named below by its command name. That name is a `.bat`
file in `CODE/` on Windows and a shell script in `CODE/unix/` on Linux
and macOS; either way it does nothing but call the `.py` file of the
same name, which sits in `CODE/`, and hand back its exit code. So
`conf_list` is exactly a run of `conf_list.py`.

You can run the `.py` files directly by

    py <file_name.py> [parameters]          (Windows)
    python3 <file_name.py> [parameters]     (Linux, macOS)

from any directory -- no PYTHONPATH is needed, because every program
puts its own folder on the import path before it imports anything.

Every program prints its own usage when it is given arguments it cannot
make sense of, and the docstring at the top of each `.py` file is the
full description of what it does and why.

You can modify any of the `.py` files as needed. However, if you find
any bugs or introduce new useful options, please notify me, so that your
changes and bug fixes could be incorporated in future releases.

## Least-squares fitting with RCE

### `add_eav`

Adds a given energy shift to all Eav parameters in the `ING11` file.
Usage:

    add_eav <energy_shift_in_kK>

The same action can be achieved by including the rescale 'card' in
`ING11`, i.e. by adding at the top of `ING11` a line with zero in
position 5 and the value of the energy shift in positions 21 through 30.
(One difference: the rescale card leaves an Eav at or below -4000,
Cowan's mark of a continuum configuration, alone; `add_eav` shifts it.)

    add_eav --ground [folder]

works on `OUTGINE` instead: it puts the lowest calculated level at zero
for the first RCE run of a fit, moving every Eav and every listed level,
and keeps RCG's file as `OUTGINE.RCG`. `rce --ground` runs it before
RCE. See
[The calculated ground level at zero](MANUAL_ADDITIONS.md#the-calculated-ground-level-at-zero-rce---ground).

### `update11`

Run after the LSF fitting with RCE, to transfer the fitted parameters
into `ING11`, run RCG with those fitted parameters, and run
[`printout`](#printout). It takes the number of parity sets to
substitute from `IN36`, so it needs no argument.

### `s11`

Substitutes the parameters from the `RCEOUT` file (produced by the RCE
least-squares fitting program) back into the `ING11` file (input for the
RCG program). Usage:

    s11

with no parameters -- for calculations with one parity set, or for
substituting the fitted parameters into `ING11` for the first parity set
only;

    s11 2

for substituting the fitted parameters into `ING11` for two parity sets;

    s11 auto

for taking the number of parity sets from `IN36` -- what `update11`
uses.

### `scale_param`

Scales the parameter values in the `RCEINP` file by given factor(s).
Parameters:

1. `input_file_name` -- name of the input file (normally, a saved copy
   of the `RCEINP` file).
2. `output_file_name` -- name of the new `RCEINP` file to produce.
3. One or more of the following tokens separated by space:

       F=<F_factor>
       G=<G_factor>
       CI=<CI_factor>
       GROUP<number>=<factor>

   where the factors given are used to multiply all parameters from the
   corresponding group of the input file. The GROUP number is the
   absolute value of the numerical flag used to group parameters in
   `RCEINP`, e.g., 51, 52, etc.
4. `conf_set` -- 1 for the first parity set, 2 for the second parity
   set, or 12 for both.
5. `FIX` -- an optional parameter. If specified, all Slater parameters
   that were 'free' in the input file will be 'fixed' in the output file
   for the given parity set(s).

### `subst_exp_levs`

Substitutes experimental levels into `RCEINP` from another file. Usage:

    subst_exp_levs <source> <target> [params=<what>] [--matchby <how>]

The `<target>` is the `RCEINP` file the observed energies are put into.
It is rewritten in place, and what it held is kept as `<target>.bak`.

The `<source>` is either another RCE file or a list of experimental
levels; which one it is does not have to be stated, since an RCE file is
recognized by being readable as one.

An RCE source is the older use of this program: you made an LSF fitting
with a limited set of configurations and then discovered you need more
configurations to include. After running RCN/RCN2/RCG/RCE with an
extended configuration set, instead of manually editing your new
`RCEINP` file again, this program transfers the level values from the
previous calculation into it. Optionally it transfers the parameters
too:

- `params=flags` will set only the flags for the same parameters in the
  new `RCEINP` (i.e., make them fixed or free or grouped in the same way
  as they were in your old `RCEINP`).
- `params=noflags` will transfer only the values of the parameters but
  not the flags.
- `params=Eav` or `params=F` or `params=G` or `params=Z` or `params=CI`
  will transfer only the values and flags for the corresponding parameter
  group (Z stands for ZETA). Both values and flags will be transferred.
- `params=none` will skip all parameter values and flags (leave them
  unchanged in the new `RCEINP` file).

An experimental source is any tab-delimited list of levels -- a
spreadsheet saved as text, the output of a least-squares level
optimization, a published table typed out. It has no parameter table,
so `params=` does not apply to it. Its layout is not fixed: the columns
that hold the energy, the parity and the J value are named on the
command line, by column number or, when the file has a header line, by
column name:

    --ecol <col>  --jcol <col>  --pcol <col>

and everything else is optional. Energies are cm-1 unless `--units kK`
or `--units eV` says otherwise. A column named `exclude`, or one named
by `--excludecol`, drops the levels marked in it.

Up to three composition triads may be given --
`--perccol`/`--confcol`/`--termcol` and the numbered forms after them --
and where they are, a level is identified by its composition rather than
by its energy. The percentages are read as intensity per cent, which is
what a published table quotes.

The configuration labels are matched against the ones the target deck
itself uses, through `IN36`. RCE writes a configuration in six columns
and abbreviates it from the left to get there, so the label alone does
not say which subshells it means; `IN36` is where they are written out,
and the card whose subshells the source names is the card it means.
`IN36` is looked for beside the target deck and then in the directory
above it. Without one the program falls back on matching a source
configuration to the deck label its compacted form ends with, which is
right for the labels that really are abbreviations and wrong for a card
that was nicknamed something else.

Only LS term labels can be matched against an LS deck. But a level
designated in JK or pair coupling still writes its parent terms inside
the configuration -- `4f2.(3F<2>).6p` beside the pair label `(2,3/2)*`
-- and the last of those is exactly the term RCE prints in brackets in
its own label, `f26p (3F) 4D`. So a component is matched on the
configuration, that parent and the final term separately, to whatever
depth the source names them, and the target's components are added up
over the rest. Only a level that gives neither an LS parent nor an LS
term is matched on (parity, J, energy) alone. The program reports how
many levels went each way.

`--matchby` says how a source level is tied to a target level:
`heuristic`, the default, uses the composition where the source gives
one and the energy where it does not; `Ecalc` uses the energy alone,
against the target's fitted energies; `Eexp` uses the energy alone,
against the target's observed ones, which is how a deck fitted early in
an analysis has its observed energies replaced by optimized ones.

An energy written with a `+x` on it -- `12345.678+x` -- is a level of a
system measured precisely inside itself but tied to the ground level
only loosely. Such a system becomes a relative-energy group in the
target deck: the letter goes into column 14 of every one of its level
records, uppercase on the base and lowercase on the rest. Which level is
the base is a decision about the measurements, so it is asked for rather
than guessed: mark it in a column of the source and name that column
with `--basecol`. A system of one level is its own base. The two
standard deviations of the group go on a card block at the end of the
deck, which this program does not write -- add it by hand, as
[Relative-energy groups in an RCE fit](MANUAL_ADDITIONS.md#relative-energy-groups-in-an-rce-fit)
describes. RCE will stop and name the letter if it is missing.

Running `subst_exp_levs.py` with no arguments, or with `--help`, prints
the whole option list with the composition columns explained in full.

### `subst_params_ing11`

Substitutes parameter values from one `ING11` file into another one. The
substitution is done line by line. If the configuration name(s) on the
line read from the input file matches those in the output file, the
parameter line(s) for these configurations are substituted in place of
similar lines in the output file. The input `ING11` file is assumed to
have both parities, while the output `ING11` is supposed to have only
one parity set. The main purpose is to substitute fitted parameters to
set up `ING11` for autoionization calculations. Usage:

    subst_params_ing11 <input_ING11> <output_ING11> <conf_set>

where `<conf_set>` can be 1 or 2, for conf. set #1 or #2 in the
`<input_ING11>`.

### `reorder_ing11`

Reorders the subshells in the `ING11` file. This is needed if you want
to change the default summation order of the subshells. For example,
the configuration 1s.2s.2p of Li-like spectra is best described by first
combining the 2s and 2p subshells to form intermediate LS terms 3P\* or
1P\*, and then combining this intermediate term with the 2S term of the
1s subshell. The default order of shell summation in RCG is left to
right, i.e., 1s is combined with 2s to produce 3S of 1S, then this is
combined with the remaining 2p shell to produce the final 4P\* or 2P\*
term. Usage:

    reorder_ing11 <output_filename> <new_order>

where `<output_filename>` is the name of the new `ING11` file (the
original `ING11` must be present in the current directory; it is used as
input), and `<new_order>` is a sequence of digits, e.g., 231 to specify
that shells 2 and 3 will go first in the new `ING11`, followed by shell
3. You can recognize the shells in the top section of `ING11` and count
their numbers from left to right, starting with 1. The program
rearranges not only the shell designations, but also the corresponding
Slater and CI parameters.

## Level and line lists

### `printout`

Prints out the level list in a form similar to `LEVELS2`, but with
configuration and term labels including proper parent terms in the order
of summation of shells, in any of the coupling schemes RCG and RCE can
produce. The eigenvector compositions with parentage are read from the
`OUTG11` file, while the experimental energies are substituted from the
`RCEOUT` file. Thus, to produce correct results, you must first
substitute the fitted parameters from `RCEOUT` into `ING11` (e.g., by
running `s11 auto`) and re-run RCG. The program also produces a file
with fitted and HF parameter values and ratios LSF/HF. Usage:

    printout <levels_out_name> [params_out_name]

where `<levels_out_name>` is the name of the output levels file
(typically, `levels.txt`), and `[params_out_name]` is the optional name
of the parameters output file (if omitted, the default is
`params.txt`). The output files are ASCII tab-delimited files best
viewed with Excel or any other spreadsheet program.

This program takes into account the possible reordering of subshells
(see [`reorder_ing11`](#reorder_ing11)). For it to work, the files
`IN36`, `ING11`, `OUTG11`, `genealogy.txt`, `OUTE`, `RCEINP.HF` and
`LEVELS1` must be present in the current directory. `genealogy.txt` is
written by RCG only when IPGEN, columns 83-84 of the `ING11` main card,
is 2. The scaling factors for the HF parameters are read from `OUTG11`.
Therefore, for the ratios LSF/HF to be correct, the first run of RCE,
which produces `RCEINP.HF`, must be made with the sequence of runs
RCN/RCN2/RCG/RCE. The files `OUTE` and `RCEINP.HF` are skipped if you
specify an optional `--no-rce` parameter in the command line or RCE was
not run at all.

The Python version is no longer limited to LS and JJ. Its options are:

- `--couplings LS,JJ` -- which coupling schemes to print, in this order.
  Any of the eight RCG and RCE can produce may be named.
- `--from LSJK=JK` -- take a coupling from a sub-folder, for the common
  case where the run that has it is a second RCE run in a directory of
  its own. Merging the two level lists used to be hand work.
- `--from LSr=LS@reordered` -- the same, under a name of your own, so
  that one levels file can hold two blocks of the same coupling -- which
  is what a reordered run is.
- `--dominant-config` -- add a conf0 column: the level's leading
  configuration as RCE determines it, from the full-precision
  eigenvector. The weight of a configuration in a level is the same in
  every coupling scheme, so this is the one configuration name that does
  not depend on the scheme. The programs that
  [turn a fitted calculation into a published level list](#from-a-fitted-calculation-to-a-published-level-list)
  need it.
- `--drop-full-shells 5p,4d` -- leave the named subshells out of the
  printed labels where they are closed.
- `--components 5   --min-percent 1.5` -- how much of each eigenvector
  composition to print.
- `--match eigenvector|energy-order` -- how a level of a second folder
  is tied to the same level here.

Labels are built from `genealogy.txt`, so the middle of a coupling chain
is read rather than guessed, and the principal quantum numbers are taken
from `IN36` through the `ING11` slots, so a reordered deck is labeled
correctly. Every label carries the full electron count.

### `conv_out`

Converts both level and line lists from the output of RCG and RCE.
Input parameters:

1. **Output levels file name** (e.g., `levels.txt`). The format of this
   particular file is slightly different from a similar file created by
   the `printout` code. Namely, it contains several additional columns.
   The most important one is the unique level number, which is used to
   label the lower and upper levels in the simultaneously created
   transitions file.
2. **Output transitions file name** (e.g., `tr.txt`).
3. **`RCE` or `noRCE`** -- for including or omitting RCE data such as
   experimental energy levels.
4. **`LTE` or `DR` or `AUGER` or `CM` or `BF:<time_of_flight_ns>`.**
   This specifies how to scale the line intensities. The most popular
   use is LTE (local thermodynamic equilibrium, i.e. a Boltzmann
   population distribution). DR stands for 'dielectronic
   recombination'. This option, as well as AUGER, requires
   autoionization rates to be available in two separate files, one for
   the first parity set and one for the second. If this option is
   specified, the names of these files must be given at the end of the
   parameter list (parameters 7 and 8).

   CM is a special 'Cascade Matrix' option in calculations of Auger
   electron spectra. With it, the levels file carries five extra
   columns, n, l, occ, S and L, which classify each level for this
   model; they are written with this option only.

   BF stands for 'Beam-foil'. For this option, an additional parameter,
   time of flight in nanoseconds, must be given together with 'BF:'. It
   means the time required for the emitting ions to cross the field of
   view of the spectrometer. The intensity is estimated by integrating
   the line intensity from each of the decaying levels within this 'time
   of flight' period.
5. **Effective excitation temperature** in 1000 cm-1. This parameter is
   required with all options, but is used only in LTE.
6. **Max. scaled intensity** -- a number meaning the ratio of max and
   min intensities of the lines to be included in the transition list.
   The min. intensity is assigned the value of 1, and the rest are
   scaled accordingly. Lines with intensity smaller than 0.5 are
   omitted. All intensities are rounded to an integer number.
7. and 8. **Optional: autoionization rate files** for the first and
   second parity configuration sets. These files are generated by the
   utility code [`read_aa`](#read_aa): one run of it per parity, made in
   the directory holding that parity's `OUTG11`, with its output
   redirected to the file named here. That working arrangement -- two
   directories, even configurations in one and odd in the other -- is
   why there are two files and not one.
9. **Optional: `IDEN`.** It is used to generate input files `ENLEV.DAT`
   and `TRANS.DAT` for the visual line identification code IDEN2
   (V. I. Azarov, A. Kramida, and M. Ya. Vokhmentsev, Comput. Phys.
   Commun. 225, 149–153 (2018); DOI:10.1016/j.cpc.2017.12.012). This
   will work only if the observed line list `DLV.DAT` and an auxiliary
   file `numset.dat` are present in the working directory (for
   explanation of those files, see the above article on the IDEN2 code).

   If the option IDEN is used, it may be followed by additional optional
   parameters (added in the year 2024 release): `truncate` and/or
   `linear` or `convert2log`.

   - The option `truncate` can be specified if one wants to reduce the
     number of energy levels in the ENLEV file by deleting the levels
     that do not have any transitions from them.
   - The option `linear` should be specified if a linear intensity scale
     is desired for the calculated line intensities in the `TRANS.DAT`
     file (default is a log scale).
   - The option `convert2log` can be used to convert the observed
     intensities in the `DLV.DAT` file to a log scale matching the
     default scale of calculated intensities. A new `DLV.DAT` file will
     then be created. Options `linear` and `convert2log` are
     incompatible.

In the transitions output file, with the LTE, BF, and DR options, the
wavelengths are given in angstroms. Those between 2000 and 20000
angstroms (or, more precisely, those corresponding to wavenumbers
between 5000 and 50000 cm-1) are given in standard air, otherwise in
vacuum. With the AUGER options, the transition energy in eV is printed
instead.

For each line, the configuration and LS term label are printed for both
levels, as well as the Landé g factors, radiative lifetimes, gA and A
values, the value of the cancellation factor, and, if autoionization
data are included, the corresponding autoionization rates and branching
fractions.

The positional parameter list above is long and strictly ordered, and a
missing argument shifts everything after it. Every setting may therefore
be given as a long option instead -- `--model`, `--temperature`,
`--max-intensity`, `--ai`, `--couplings`, `--components`,
`--min-percent`, `--iden`, `--truncate`, `--linear`, `--convert2log` --
and the two forms may be mixed freely.

Better still, every setting may be written into a `conv_out.par` file in
the working directory, in the same layout as `calc_TP.par`: the value
first, then a semicolon, then the key, then prose saying what it is for.
Run `conv_out` with no par file and no arguments and it writes a fully
commented sample and stops. The par file supplies the defaults and the
command line overrides them, so the historical positional form above
keeps working exactly as it did.

This program began as a reimplementation of an older Borland Pascal code
of the same name, and a few of that program's options have still not
been carried over. What it gains in exchange is that it places no limit
on the number of levels or transitions, it is portable, and it can write
any of the coupling schemes RCG produces rather than LS and JJ alone.

### `read_aa`

Processes the output of RCG for autoionization calculations (see the
autoionization section of Cowan's [RCG_DOC.txt](Cowan/RCG_DOC.txt)).
Usage:

    read_aa <cores> <IP> [--set N] [--folder F]

`<cores>` is a comma-separated list of `name:energy` giving the core
states of the ionized atom with their energies above that ion's ground
state, in the run's own energy unit, and with no blanks anywhere in it.
A name is a regular expression matched against RCG's channel labels, so
one column of the table may sum several of them. `<IP>` is the
ionization energy of the atom being calculated, in the same unit, which
is what puts the core energies on this ion's scale. The table goes to
standard output; redirect it to the file [`conv_out`](#conv_out) is to
be given.

Running `read_aa.py` itself -- `py read_aa.py` from `CODE/`, or
`python3 read_aa.py` -- explains the rest and gives worked examples.

### `conf_list`

Prints out the lists of configurations included in the first and second
parity sets of `IN36`.

## Level identifiers and the round trip through IDEN2

### `insert_level_ids`

Inserts unique level identifiers into the level list and the transition
list that [`conv_out`](#conv_out) writes. The identifiers come from a
tab-delimited list of experimental levels -- the same kind of file
[`subst_exp_levs`](#subst_exp_levs) takes as its source. Usage:

    insert_level_ids <source> <levels file> [<lines file>]

The two target files may equally be the two a [`calc_TP`](#calc_tp) run
writes -- `out_levels` and `out_file` -- which carry the same `Lev_#`,
`lid1` and `lid2` columns for exactly this purpose.

Both target files are rewritten in place, and what they held is kept as
`<file>.bak`. The line list may be left out, and then only the level
list is done.

A level of the source is tied to a level of the list on (parity, J,
energy), the way `subst_exp_levs` ties one whose composition it cannot
use: the Ee column holds the observed energies of the RCE deck, which
came from a list like this one in the first place, so the energies are
not so much compared as recognized. Every pair within the tolerance is
considered, the closest first, and each level of either side is used
once, so two levels of one J cannot both take the level between them.

The columns of the source are named on the command line, by column
number or, when the file has a header line, by column name:

    --idcol <col>  --ecol <col>  --jcol <col>
    --pcol <col>

`--idcol` is the column holding the identifiers and defaults to one
named `level_id`. Energies are cm-1 unless `--units kK` or `--units eV`
says otherwise, and `--tol <value>`, in those same units, says how far
apart two energies may be and still be one level (1 by default). A
column named `exclude`, or one named by `--excludecol`, drops the levels
marked in it.

The identifier goes into the level list as a `level_id` column beside
Ee, and into the transition list as `id1` and `id2` -- the lower level's
and the upper level's -- beside the `lid1` and `lid2` columns they are
looked up through. Nothing is matched twice: `lid1` and `lid2` are the
`Lev_#` of the level list.

Columns of those names already in a file are updated rather than
duplicated. They keep their place, are emptied, and are filled with what
this run matched, so a level that has lost its identifier since the last
run loses it here too instead of keeping a stale one.

The program reports how many levels matched, which source levels found
nothing, how many levels of the list were left without an identifier,
and any match made against a calculated energy rather than an observed
one.

### `write_IDEN_level_ids`

Writes the table that translates level identifiers into the level
numbers IDEN2 uses. It is the first half of the round trip through
IDEN2; [`update_IDEN`](#update_iden) is the second.

A level of a published analysis is known by its identifier -- the key a
database holds it under, and the same name in every paper and every run.
IDEN2 knows nothing of it: it numbers its own levels, from the highest
energy down, over whatever `conv_out` gave it, and that number is the
only name a level has in `ENLEV.DAT` and `TRANS.DAT`. It means nothing
outside one pair of those files -- rebuild them and every number moves
-- so the correspondence has to be written down before it can be used.
Usage:

    write_IDEN_level_ids <source>
        [--idcol <col>] [--ecol <col>] [--jcol <col>]
        [--pcol <col>] [--excludecol <col>]
        [--units cm-1|kK|eV] [--skip <n>]
        [--tol <value>]

The `<source>` is the tab-delimited experimental level list
`subst_exp_levs` and `insert_level_ids` take, and its columns are named
the same way. It is the only argument: `ENLEV.DAT` is read from the
working directory under that name, because that is the name IDEN2 gives
it, and the table is written beside it as `IDEN_level_ids.txt` -- two
columns, `level_id` and `IDEN_id`, tab separated, in level order.
Neither input is touched.

Only the FIXED levels of `ENLEV.DAT` are matched -- the ones carrying
the star, which in IDEN2's terms are the experimental ones -- and the
match is the one `insert_level_ids` makes, on (parity, J, energy) with
every pair inside `--tol` offered and the closest taken first. The
energy compared is the fixed one and not the calculated one: it is the
observed value, and an experimental level list holds observed values.

Parity comes from the sixth character of the IDEN2 label, `_` for even
and `~` for odd. That reading is checked rather than assumed: the
matching is tried a second time with the marks reversed, and the run
says so if the wrong reading would tie more levels -- which is what an
`ENLEV.DAT` from a different calculation looks like.

A source level `ENLEV.DAT` has no counterpart for is reported and left
out of the table. That is the usual case rather than the exception:
`conv_out` writes at most 9999 levels and can truncate a large
calculation to the levels its printed transitions touch, so a
compilation reaches IDEN2 with part of itself missing.

### `update_IDEN`

Carries the result of a level optimization back into the two files
IDEN2 reads, so that the next look at a questionable level is taken
against the current answer rather than the one the files were built
with. Usage:

    update_IDEN <LOPT levels> [<LOPT lines>]
        [--ucalc <value>] [--ids <file>]
        [--wntol <value>]

`<LOPT levels>` is the tab-delimited level table LOPT writes, read by
its Designation, Energy, D1 and D2tot columns; `<LOPT lines>` is its
line table, read by L1, L2, wn_o and Weight. Leaving the line table out
updates the levels only. `ENLEV.DAT`, `TRANS.DAT` and `dlv.dat` are
taken from the working directory under the names IDEN2 gives them, and
what each rewritten file held is kept beside it as `<name>.bak`.

Three things are done in one pass:

1. Every FIXED level of `ENLEV.DAT` takes the energy LOPT gives it and,
   as its uncertainty, max(D1, D2tot) -- the larger of the purely
   statistical uncertainty and the total one relative to the ground
   level. That is what a level is worth when it is looked at on its own,
   which is what IDEN2 does with it.
2. With `--ucalc`, every level that is not fixed takes that uncertainty.
   The calculated levels come out of `conv_out` at 5000 cm-1, which is a
   placeholder; the uncertainty is the window IDEN2 searches around a
   level, and the point of the option is to narrow it once the fit is
   good enough to deserve it.
3. `TRANS.DAT` takes the same energies -- the `$` record's, the `+`
   record's lower level, and the transition energy that is their
   difference -- and, for a transition between two fixed levels, the
   assignment it should have.

The third is the one that matters. IDEN2 reads exactly one field of a
`+` record beyond the two level numbers: the assigned line's row number
in `dlv.dat`. The set of transitions carrying one has to be the set of
lines the optimization used, or the next session in IDEN2 is looking at
the previous round's answer.

A line counts as used when its Weight is above zero. That is the flag to
rely on and not the F column: P, M and virt all mean a line was
excluded, but the weight says so for an excluded line that carries no
letter as well. Where one transition appears on several rows -- the same
pair measured twice, a blend entered twice -- the row with the largest
weight speaks for it. A transition between two fixed levels that the
optimization did not use loses its assignment, and the observed
intensity, wavenumber and dwn beside it are cleared with it, so that the
record does not go on quoting an observation it no longer has. A
transition reaching only one fixed level, or none, was never in the fit
and keeps whatever IDEN2's own session gave it.

LOPT names a level by its identifier and IDEN2 by a number of its own,
so the two are tied through `IDEN_level_ids.txt`, which
[`write_IDEN_level_ids`](#write_iden_level_ids) writes and `--ids` can
rename. That table must have been made from this `ENLEV.DAT`: a table
naming levels the file has not got is refused by name, which is the one
mistake this pair of programs makes easy.

LOPT names a line by its observed wavenumber and IDEN2 by its row of
`dlv.dat`, so those two are tied by searching the observed list.
`--wntol` is how far apart they may be: 0.05 cm-1 by default, because
LOPT prints a wavenumber rounded to the decimals its uncertainty
deserves and the two then disagree in the last place and no further.
The `dlv.dat` in the working directory has to be the one `TRANS.DAT` was
assigned against -- the assignment is a row number of it and means
nothing against another list -- so every wn_o of the line table is
looked for in it before anything is written, and the count and the worst
offset are reported. Nothing is refused on the strength of that: a
tolerance set too tight shows as a few near misses and a mismatched file
as thousands, and the two read differently.

A line the optimization used whose transition is not in `TRANS.DAT`
cannot be assigned there, and the count is reported. `conv_out` writes
only the transitions above its intensity threshold and only between the
levels that survived the truncation, so some of the used lines will have
no record to go into. The remedy is a fuller `TRANS.DAT`, not a
fabricated record: `TRANS.DAT` holds calculated transitions, and a row
invented for an observation has no calculated intensity to put in it.

## Transition-probability uncertainties

### `calc_TP`

Runs Monte Carlo trials for estimation of uncertainties of E1, M1 and E2
transition rates (see A. Kramida, Atoms 2, 86-122 (2014),
https://doi.org/10.3390/atoms2020086). This was `calc_TP_M1E2`, which
covered the two forbidden multipoles only; the old name is retired and
is not accepted as an alias.

The settings that used to be hard-coded in the Perl source are now in
`calc_TP.par`, a commented parameter file. Run the program with no par
file present and it writes one into the working directory and stops, so
that you can read it through and fill it in; `--write-par` writes a
fresh one at any time. `multipoles` chooses which of E1, M1 and E2 to
do, and the default is all three. A deck with only one parity set has
no E1 transitions, so the run says so and drops E1 rather than failing.

`calc_TP` needs `numpy`. It pairs each trial's eigenvectors with the
unvaried run's by their overlap, block by block, which is the only way
the identification survives a trial that has moved the eigenvalues; that
is a matrix product per *J* block, and without `numpy` the program stops
and says to install it.

The RCG run it reads must print the lines it is asked about. For M1 and
E2 that needs two things on the `IN2` main card. IMAG and IQUAD, columns
49 and 50, switch the M1 and E2 spectra on. DMIN, columns 61-65, must be
**0.00**: a DMIN between 1 and 50 means "omit every transition weaker
than gA = 10\*\*DMIN", and forbidden lines are almost all weaker than the
5.00 or 10.0 that suits electric-dipole work, so the spectrum comes out
empty. See `WORK/Forbid/README.txt`.

For E1 the useful DMIN is whatever cut you mean to publish at, and it is
worth punching ISPECC (column 72) as **3** rather than 2. Only a level
that heads a printed block carries RCG's own total decay rate, SUMGA,
which is where the branching fractions come from and which -- unlike a
sum over the printed lines -- is not censored by DMIN. With 3 both level
sorts print, so every level has one; the repeated print of each line is
folded. `WORK/Pr3/README.txt` gives the deck.

The run writes two files. `out_file` is the line list; `out_levels`,
beside it, is the level list holding the `Lev_#` that the line list's
`lid1` and `lid2` columns name. The numbering is the one `conv_out`
uses, and the two files have the shape
[`insert_level_ids`](#insert_level_ids) reads, so the identifiers of a
published analysis go onto these lines the same way they go onto
`conv_out`'s -- which is what lets the two programs' tables be joined.
Join them on those identifiers, not on `Lev_#`: a level number is a
position in one run's energy ordering, and two runs whose eigenvalues
have moved do not order their near-degenerate levels alike.

Trial subdirectories are kept, so a trial can be read again without
being run again. A long run makes a great many of them: set
`keep_trial_output = 0` to delete each trial's `OUTG11` once it has been
read, and

    calc_TP --cleanup

to remove the trial directories afterwards. It lists what it will
delete and asks once; `--yes` skips the question. Warnings go to
`calc_TP.log` in the working directory, and the terminal gets one line
saying how many there were.

### `extract_Icalc`

Pulls out of a `calc_TP` line list the few columns a line classification
code needs, and, if a temperature is given, a synthetic intensity built
out of them. Reads `tp_E1.txt` with `id1` and `id2` already put on it by
[`insert_level_ids`](#insert_level_ids), and writes `Icalc.txt`. Usage:

    extract_Icalc [<tp file>] [<output file>]
                  [--T_eV <value>] [--d_off <value>]

The columns written are `id1`, `id2`, `rwn` (E_up_exp - E_low_exp), `gA`
((2\*J2+1) times the A(s-1) column), `u%gA` (std_A%, which is a relative
standard deviation and so is the same for gA as for A) and `Eup`
(E_up_exp).

`calc_TP` carries its energies in kK and this writes cm-1, so `rwn` and
`Eup` are its values times 1000.

Given BOTH `--T_eV`, a temperature in eV, and `--d_off`, a log-offset, a
seventh column follows:

    Icalc = exp(-Eup/T_cm + d_off) * gA / rwn

with T_cm the temperature converted to cm-1 -- a Boltzmann population of
the upper level times its rate, per unit wavenumber. `d_off` is the one
free normalisation, absorbing the partition function, the column density
and the choice of unit for `rwn` together, so it is given rather than
computed. One of the two options without the other is refused: the
column is then not written at all rather than written empty.

Only transitions with an identifier at both ends are taken, and of those
only the ones with an experimental energy at both ends -- a line to a
level the analysis does not hold has no identifier to be classified
under, and one to a level the fit did not observe has no experimental
energy to be classified by. The run reports how many it passed over, and
why.

## From a fitted calculation to a published level list

Three programs with no Perl ancestor, run in this order, each handing the
next a file. The middle file is the one a person edits:

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

They rest on the invariance noted under `--dominant-config` of
[`printout`](#printout): because recoupling is block-diagonal in the
configuration, the weight of a configuration in a level is the same in
every scheme. Purities of two schemes are therefore comparable, and a
composition may quote each configuration in its own scheme without the
percentages ceasing to add up.

Prepare the input with

    printout levels.txt --dominant-config
             --couplings LS,JJ,LSJK --from LSJK=JK
             --components 10 --min-percent 0.5

Parity is the first column of all three outputs: deleting a column for
an article table is easy, reconstructing one is not.

### `config_purity`

Which coupling scheme describes each configuration best. Writes
`coupling_choice.txt`: one row per configuration in `IN36` order, one
10-character column per coupling scheme found in `levels.txt`, and in
each cell the purity -- the unweighted mean, over the levels of that
configuration, of the leading component's percentage. A star marks the
chosen scheme, one per row, and the file is written with the star
already on the highest purity. That is a suggestion, not a verdict:
previous literature on the spectrum, and the schemes used for
neighboring isoelectronic ions in the database, routinely outweigh a few
percent gained by a rarely used scheme. Move the stars by hand.

The file is fixed-width because it is meant to be edited. It is never
silently discarded: writing over an existing one needs `--force`, and
`--merge` updates the numbers while leaving the stars where you put
them.

Options: `--levels`, `-o`, `--force`, `--merge`.

### `make_levels_list_prepub`

One coupling scheme per configuration, applied. For each level it takes
from each scheme's block only those components whose configuration was
assigned that scheme, concatenates them across schemes, sorts by
percentage and prints. A level mixing two configurations may therefore
read 44% in JJ for one and 35% in LS for the other. Writes
`levels_prepub.txt`.

A level belonging to a relative-energy group (see
[Relative-energy groups in an RCE fit](MANUAL_ADDITIONS.md#relative-energy-groups-in-an-rce-fit))
carries the group's letter in a suffix column after Ee, written as
` +x`, and `make_levels_list_pub` carries that column to the published
list. The level table `conv_out` itself writes gains a `grp` column when
a letter is present.

Options: `--levels`, `--choice`, `-o`, `--components`, `--min-percent`.

### `make_levels_list_pub`

Designations unique within each parity and J, which is what a published
or database level needs. The name comes from the leading component where
the level is reasonably pure; where it is not, from further down the
composition; failing that, from letter prefixes (a, b, ...) assigned per
configuration and term across all J, in order of increasing energy.
Where no meaningful name can be assigned at all, the name is the integer
part of the energy. Writes `levels_pub.txt`, in cm-1.

`--nlev` takes one number per configuration set, in `IN36` order, saying
how many levels to publish; the default is the position of the highest
observed level in each set, and a set with no observed level makes
`--nlev` mandatory.

Options: `--nlev`, `--n_comp`, `-o`, `--min-name-percent`,
`--letters-max`, `--leading-min`, `--max-outside`.

## Supporting code

### `toolbox/`

The shared library the programs are built on: the readers for `IN36`,
`ING11`, `OUTG11`, `genealogy.txt` and the RCE files, the ASD notation
formatter, and the pieces the individual programs assemble. Four of its
modules will also print what they have read, which is useful when a deck
is not doing what you expect:

    py -m toolbox.in36 [IN36]
    py -m toolbox.outg11 [OUTG11]
    py -m toolbox.genealogy
    py -m toolbox.jjkeys

(`python3 -m ...` on Linux and macOS.) They must be run with `-m` as
shown, from a directory that has `CODE/` on PYTHONPATH, or run the
launcher of the program you want instead. They are library modules, not
commands.

Two of its modules are worth naming:

- **`toolbox/vacair.py`** -- the module `conv_out` uses to convert
  wavelengths between vacuum and standard air. The functions are
  `Lvac()` and `Lair()`; their arguments are the air and the vacuum
  wavelength in angstroms, respectively. This was `vacair.pl`.
- **`toolbox/fundconst.py`** -- the Python copy of `FOR/fundconst.f`:
  the same constants, under the same names, derived in the same order.
  `toolbox/radiate.py` takes all of its conversion coefficients from
  here, so a change to the Fortran file reaches the utility codes too. A
  test compiles `FOR/fundconst.f` and checks every value against this
  module bit for bit.

### `lsjj/`

The generator `make_cfp` runs to build `LSJJ.DAT`: `mkdeck.py`, which
writes the deck, over `lsjj.py`, `cowan_cfp.py` and `recur.py`, which
are the algebra. It is here rather than in `FOR/` because `make_cfp` has
to find it beside itself in an installed tree. It is the one part of the
package that needs `numpy` at installation. The record of the derivation
and its validation stayed in `FOR/LSJJ/`.

### `tests/`

The test suite (`py -m pytest`, or `python3 -m pytest`, from `CODE/`).
It is a development file and `install.sh` does not install it; run it
from the unpacked tree. A large part of it runs the retired Perl
programs from `perl_legacy/` side by side with the Python and compares
the output, which is the evidence that the port is faithful. Those tests
read a finished run and skip themselves when the case in `WORK/` has not
been run -- the cases ship as inputs only. Run one as its `README.txt`
describes and its tests come back.

### `perl_legacy/`

The retired Perl utilities, unchanged. Nothing in the package calls them
and they are not on the PATH; they are kept as the reference the Python
was accepted against, and because they are the record of what the codes
used to do.
