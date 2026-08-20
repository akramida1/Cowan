WORK/Pr3 -- Pr III, a real least-squares fit in three coupling schemes
======================================================================

Thirty-nine configurations of doubly ionized praseodymium, both
parities, fitted to observed energies with RCE and then written out as
a publication-ready level list in which each configuration carries the
coupling scheme that describes it best.  It is the worked example of
everything downstream of a fit: the RCEOUT/RCEINP input format, the
update11 cycle, the cross-folder coupling merge, and the three programs
that turn levels.txt into levels_pub.txt.

INPUTS ONLY

The run writes about 412 MB, so nothing generated is kept here.  What
is in the folder:

    IN2, IN36     the deck, as for any other case in WORK/
    OUTGINE.fin   RCE's control file, converged
    RCEINP.fin    RCE's vector-input file, converged
    JK/ING11      the same fit, asked for in LSJK instead of JJ
    dlv.dat       the observed line list, for the IDEN2 output below
    numset.dat    the wavelength ranges dlv.dat is divided into

A real fit of a spectrum this size is dozens of RCE iterations with the
level flags and the parameter groupings adjusted by hand between them.
That history is not reproducible from IN36 and is not the point of the
example.  The two `.fin` files are the state it converged to, and the
run below reproduces its **last iteration**, which is what every number
in the published list comes from.

    * OUTGINE.fin differs from the OUTGINE that RCG writes in two of
      the seven fields of its first line: the third, the number of
      iterations, is 1 rather than 5, and the fifth, INVECT (columns
      21-25), is 1 rather than 0.  INVECT=1 tells RCE to take the
      levels and the parameters from RCEINP and to ignore everything in
      OUTGINE below the first control line of each parity.
    * RCEINP.fin is the RCEOUT of the iteration before, renamed.  It
      carries the observed energies, the fitted parameter values, and
      the fixed/free/group flag on every one of them.

THE RUN

Work in a copy of this directory, not in the distribution.

The commands below are written as they are typed on Windows.  On
Linux and macOS type them in lower case -- rcn, rcn2, rcg, rce, and
likewise for the utility codes -- and read `cp` wherever this file
says `copy`.  Nothing else differs.

    RCN
    RCN2
    RCG
    RCE                             <- do not skip this one

RCN2 writes ING11 straight from the G5INP card of IN2, and it needs no
editing.  RCN2 punches IPGEN=2 into columns 83-84 of the ING11 main
card itself (FOR/RCN2K.F:1593), so RCG writes genealogy.txt, which is
where the labels come from; nobody has to remember the column.

The first RCE run is the one that fits nothing.  There are no observed
energies in the OUTGINE that RCG wrote, so it stops after cycle 0 --
but on the way it writes RCEINP.HF, the copy of the vector-input file
holding the **Hartree-Fock** parameters.  That file is what printout
divides the fitted parameters by to get the LSF/HF ratios in
params.txt, and nothing later in this sequence recreates it.  Skip this
run and the ratios are wrong.  (Equally, if an RCEINP.HF left over from
some other deck is sitting in the directory, delete it first.)

Now put the converged fit in place and run the last iteration:

    copy OUTGINE.fin OUTGINE
    copy RCEINP.fin  RCEINP
    RCE

It reports

    Iteration No 1 finished with AVDEV=   0.133

-- a mean deviation of 133 cm-1 over the fitted levels.  Then

    update11 2

with the 2 because IN36 has two parity sets (RCN2 reports the split:
configurations 1-21 are the first, 22-39 the second).  update11 is a
batch file and does four things: it saves OUTGINE, runs s11 to
substitute the fitted Slater parameters from RCEOUT back into ING11,
re-runs RCG so that the eigenvectors and the transition rates belong to
the fitted parameters rather than to the Hartree-Fock ones, restores
OUTGINE, and runs `printout levels.txt`.

That last printout is the default one, LS only.  The three-scheme list
is made below and overwrites it.

THE SECOND COUPLING

RCG will write two coupling schemes into one OUTG11, and this deck asks
for LS and JJ.  A third scheme needs a second run, which is what JK/ is
for.  Its ING11 is the fitted ING11 that update11 just wrote, with two
fields of the main card changed:

    columns 31-39   KCPLD    0011111  ->  0111011
    columns 23-25   NLSMAX   (blank)  ->  999

KCPLD(I)=0 asks for coupling I, so 0011111 is LS + JJ and 0111011 is
LS + LSJK -- coupling 5 in RCG's list (`FOR/rcg11kd.f:1075`).  LS is
kept because it is what ties the two runs together: the folders are
matched level by level through their LS eigenvectors.

The 999 in columns 23-25 is NLSMAX, and since 2026 it is no longer
needed -- a blank field means "no limit" and permits couplings 3-7,
where it used to forbid them (`FOR/rcg11kd.f:735-756`).  It is punched
here because the file predates the change and because it is explicit.

The parameter cards of JK/ING11 hold the same fitted values that
update11 produces, to the last digit; only the exponent notation of the
zero and near-zero cards differs (`0.00e+01` where RCN2 wrote
`0.00E+001`), because the shipped file has been through an earlier s11
and the generated one has not.  Both read as the same number.  So it
can be used as it stands, or replaced by a copy of your own ING11 with
those two fields edited.

    cd JK
    RCG
    RCE
    cd ..

RCG in JK/ reads its own ING11 and needs nothing from the parent; RCE
after it runs another cycle 0, which is all that is wanted -- the
parameters are already the fitted ones, so its levels are the fitted
levels, and its job here is only to identify them.  IN36 is not copied:
printout reads it from the parent folder.

THE LEVEL LIST IN THREE SCHEMES

    printout levels.txt --dominant-config
             --couplings LS,JJ,LSJK --from LSJK=JK
             --components 10 --min-percent 0.5

(on one line).  `--from LSJK=JK` is the merge that used to be done by
hand: it reads the LSJK block out of JK/OUTG11 and pairs each of its
levels with the same level here through the LS eigenvectors.  It
reports

    LS, JJ: 25 J blocks, 5626 levels
    LSJK (LSJK) from JK, matched by eigenvector
    wrote levels.txt (16878 levels)
       LSJK J=12.5 of set 2 has one basis state, so RCG printed no
       LSJK section; it is named in LS

16878 is 3 x 5626: one block per scheme, every level in all three.  The
last message is not a fault.  A J block with a single basis state has
nothing to recouple, so RCG prints no alternate-coupling section for
it, and printout falls back to the LS name -- which is the same state
under any name.

`--components 10 --min-percent 0.5` are wider than printout's defaults
of 5 and 1.5% on purpose.  The next two programs take components from
this file and cannot recover one that was cut here.

FROM THE FIT TO A PUBLISHED LIST

Three programs, in order, the middle file edited by hand:

    config_purity

writes coupling_choice.txt: one row per configuration, one column per
scheme, and in each cell the purity -- the mean, over the levels of
that configuration, of the leading component's percentage.  On this
deck it reports

    wrote coupling_choice.txt (39 configurations)
       suggested: LS for 19, JJ for 13, LSJK for 7
       15 configurations are within 2% between two schemes and are the
       ones worth a decision: ...

A star is already on the highest purity in every row.  **It is a
suggestion.**  Highest purity is the main criterion but not the only
one: the literature on the spectrum, and the schemes used for the
neighboring isoelectronic ions, routinely outweigh the few percent
that the 15 near-ties above are decided by.  Move the stars, keeping
exactly one per row.  Nothing else in the file may be edited; it is
read back by column.  Re-running config_purity over a marked file needs
--force, or --merge to redo the numbers and leave the stars alone.

    make_levels_list_prepub

applies the marks.  For each level it takes from each scheme's block
only the components whose configuration was assigned that scheme,
concatenates them and sorts by percentage, so one composition may quote
44% of one configuration in JJ and 35% of the next in LS.  The
percentages still add up, because the weight of a configuration in a
level does not depend on the coupling scheme -- recoupling is
block-diagonal in the configuration.  That invariance is what makes the
whole of this section legitimate, and it is what Sugar's Pr III tables
already do.  Writes levels_prepub.txt (5626 levels).

The run prints a note about 2016 levels covering less of their weight
after the merge than one block alone did.  It is expected and the
message says why: the scheme a configuration was assigned may spread it
over more basis states, so more of it falls past printout's per-block
cut.  Raise --components if that weight is wanted.

    make_levels_list_pub

names them.  A published or database level needs a designation no other
level of the same parity and J carries; the leading component supplies
it where the level is reasonably pure, and where it is not the name is
taken from further down the composition, then from letter prefixes,
then from the energy.  It reports

    publishing 1036 (o 481, e 555), to the highest observed level of
    each set
    45 levels (4%) share a name with another of their parity and J and
    take a letter prefix
    wrote levels_pub.txt (1036 levels, 3 composition groups)

1036 of 5626, because levels above the highest observed one of their
set are extrapolations and are not published.  --nlev overrides that
count; on this deck it is not needed, since both sets have observed
levels.  levels_pub.txt is in cm-1 -- everything above it stays in
RCE's own 10^3 cm-1.

A LINE LIST WITH MODELED INTENSITIES

The same run also carries the transition array, and conv_out turns it
into a line list whose intensities follow a chosen excitation model.
Run it in the same directory, after update11, so that OUTG11 holds the
fitted eigenvectors and RCEOUT the fitted energies:

    conv_out levels_13kk.txt tr_13kk.txt RCE LTE 13 1e8 IDEN truncate

reading, in the order the arguments stand:

    levels_13kk.txt  the level list to write
    tr_13kk.txt      the line list to write
    RCE              take the observed energies from RCEOUT and compute
                     Ritz wavelengths with them; noRCE would use the
                     calculated energies alone
    LTE              the excitation model -- a Boltzmann population of
                     the upper level.  LTEABS, BF, DR, AUGER and CM are
                     the others
    13               the effective temperature, in the unit of the run,
                     so 13 kK
    1e8              the intensity given to the strongest line;
                     everything else is scaled to it.  Required in
                     every model, including the ones that ignore it
    IDEN             also write ENLEV.DAT and TRANS.DAT, the two input
                     files of the IDEN2 line-identification code
    truncate         leave out of ENLEV.DAT the levels that no
                     calculated line connects to any other

The two output names are the author's own and record the temperature;
nothing in the code reads them.  The same settings can be written into
a conv_out.par instead, and `conv_out --write-par` writes a commented
sample.

IDEN needs the observed lines.  dlv.dat is that list -- intensity,
wavenumber, wavelength, uncertainty -- in **descending wavenumber**,
and numset.dat divides it into wavelength ranges, each with a default
uncertainty for the lines in it that carry none of their own.  Every
calculated line between two observed levels is looked up there, and
the observed line nearest in wavenumber, within the uncertainty, is
the one it is assigned to.  The two files here are the author's, cut
to this spectrum; they are examples of the format, not a recommended
line list.

The run takes about seven seconds and reports

    LS, JJ: 25 J blocks, 5626 levels
    795 transitions
    wrote levels_13kk.txt (11252 levels)
    wrote tr_13kk.txt
    wrote ENLEV.DAT (555 levels) and TRANS.DAT (792 lines)

795 lines out of a spectrum this size is not an error, and it is worth
knowing where the number comes from: DMIN in columns 61-65 of the deck
is 10.0, and RCG keeps only lines with gA above 10**DMIN = 1E10 /s.
Its own report, in the update11 log, is

    2564398 LINES OMITTED, WITH MAX GA=   0.99836E+10
       5834 LINES OMITTED, WITH LAMBDA > 500000.0000 ANGSTR.

Lower DMIN and re-run RCG for a fuller list -- 0.00 keeps everything,
which is what a forbidden-transition calculation needs and what it
costs in size.

Energies in both output files are in the unit of the run, 1000 cm-1,
so a wavelength in angstroms is 1E5 divided by the difference.

WHAT IT COSTS

On a 2026 desktop the whole sequence is under two minutes:

    RCN 1 s, RCN2 <1 s, RCG 10 s, RCE 4 s, RCE again 5 s,
    update11 10 s, JK: RCG 23 s + RCE 3 s, the five utilities ~30 s,
    conv_out 7 s

and about 412 MB of disk, nearly all of it OUTG11 (85 MB here, 111 MB
in JK/) and OUTE (73 MB in each).  Neither is needed once levels.txt
and params.txt exist.

WHAT THIS DECK IS NOT

It is not the regression case.  Bit-exactness of the Fortran is checked
against a separate deck that carries E0MAX cards and a smaller DMIN;
this one is set up to be run, read and published from.  The two agree
where it matters: TAPE2E from the sequence above is byte-identical to
the regression run's.
