#!/usr/bin/env python3
"""Uncertainties on transition probabilities.

Grown out of `calc_TP_M1E2.pl`, the last of the Perl utilities, and now
covering the allowed lines as well as the forbidden ones.  It answers
the question Cowan's codes do not: *how well is this transition
probability known?*

There is no analytic answer, so it is obtained by resampling.  RCE has
fitted the Slater parameters to the observed levels and printed a
standard deviation for each of them in OUTE; RCN2 has computed the
radial transition integrals, whose accuracy is a fractional guess.
This program displaces all of them at random within those
uncertainties, runs RCG again in a trial subdirectory, pairs the
trial's eigenvectors with the unvaried run's, and collects the
line strength and the cancellation factor of every line.  Over many
trials the spread of ln S is the uncertainty, and it is reported as a
percentage, which is what a log-normal distribution's relative width
is.

The result is one tab-separated line list -- every transition between
the fitted levels, in whichever of E1, M1 and E2 the `multipoles`
parameter asks for, with A, S, the cancellation factor, the total
radiative decay rate of the upper level, the branching fraction, and
the standard deviations of the last two.

Beside it goes a level list, `out_levels`, holding the `Lev_#` that the
`lid1` and `lid2` columns of the line list name.  The numbering is
`levels.numbered`, the rule `conv_out` numbers by, and the two files
have the shape `insert_level_ids` reads, so the identifiers of a
published analysis can be put on these lines as they are put on
`conv_out`'s -- which is what joins the two programs' tables.  Joining
them on `Lev_#` alone would be wrong wherever the two came from
different calculations: a `Lev_#` is a position in one run's energy
ordering, and two runs whose eigenvalues have moved do not order their
near-degenerate levels alike.

WHY THE TRACKED QUANTITY IS S AND NOT A

RCG prints gA, and gA carries the trial's own eigenvalues inside it:
gA goes as sigma**3 for a dipole and sigma**5 for a quadrupole, so a
level the fit did not observe drags its lines with the cube or the
fifth power of the energy it happens to land on in that trial.  The
spread of A over the trials would then be the spread of the matrix
element and the spread of the energy together, and the two are not the
same uncertainty: the energy of an observed level is known to the fit,
and the energy of an unobserved one is a different question that this
program is not answering.

So each trial's printed gA is turned back into a line strength using
**that trial's own calculated energies** -- the wavelength RCG printed
on the same record -- and it is S that is accumulated.  A is then
quoted from the mean S at one fixed wavelength: the observed one where
both levels are observed on one scale, the unvaried run's fitted one
otherwise.

Where both levels are observed the wavelength does not move between
trials at all, and the relative spread of S is exactly the relative
spread of A; the two differ only when a level is unobserved and the
fit lets its energy move.  That is the case S was chosen for.

The S column of the line list is the unvaried run's, S_mean the
geometric mean over the trials, and std_A% the relative standard
deviation of S -- which is the relative standard deviation of A, in
the fixed-wavelength sense just described.

A line the trials printed **once, or not at all**, carries a blank in
every uncertainty column: one sample has no spread, and a zero there
would read as a line that is known exactly.  The n_trials column says
how many trials the line was seen in, and it is load-bearing -- a line
seen in fewer trials than were run was censored out of the rest by
DMIN or by the wavelength window, and its statistics are those of the
trials that kept it.

THE TRIALS ARE RUN IN A WIDER WINDOW THAN THE LIST IS WRITTEN IN

DMIN and the FLBDN/FLBDX wavelength window censor RCG's printed line
list, and a line just above either of them is printed in some trials
and not in others.  Its mean is then a mean over a censored sample --
biased high -- and its n_trials is smaller than the number of trials
run.  Measured on Pr III at DMIN=8.0, 9462 lines, 5 trials: 1407 lines
were seen in 1 to 4 trials and 1773 more appeared only in trials.

So the trial decks are run one decade lower (trial_dmin_offset) and in a
wider wavelength window.  The two ends of that window are not widened
by the same fraction.  A level energy moves by a per cent or two, and
trial_lambda_margin -- 5 per cent -- is the right allowance at the
short end; but a transition energy is a *difference* of two level
energies, and between two highly excited configurations it is small, so
the same absolute movement is a far larger fraction of it.  FLBDX cuts
on exactly those smallest transition energies.  The long end is
therefore widened by trial_lambda_margin scaled up by the ratio of the
two configurations' EXCITATION energies to the smallest gap between two
of them the spectrum actually connects, bounded by
trial_lambda_margin_max.  See `long_margin` -- and read its note on why
the energies have to be measured from the ground configuration.

DMIN itself is the deck's, and whose it is depends on the multipoles.
Where E1 is in the run the forbidden lines are part of the same
spectrum and share the DMIN chosen for it; where the forbidden lines
are all that is being analysed there is no reason for a threshold at
all, and the trial decks are punched at DMIN = 0.0.  The run says which
of the two it did.

Only the trial decks are touched; the unvaried run is read as it
stands.  E0MAX is left alone -- it is a scope choice, not a threshold
noise pushes a line across.  See `censor`.

A trial may also print a line with its two levels the OTHER WAY ROUND:
which end is upper is decided by energy in that run, and a level the
fit did not observe moves.  The reader tries the reversed key before it
calls a line new, so the base transition keeps that trial in its
statistics; the E_low/E_up columns stay the unvaried run's, which is
the calculation being published.  The count is in the log.

WHAT IS IN THE LIST: THE UNVARIED RUN'S TRANSITIONS, AND ONLY THOSE

The widened trial window brings lines in that the deck's own cuts keep
out, and those are not results: the list is the transitions the
unvaried run printed, every one of them and no others.  The wider
window is there to complete the statistics of the lines that ARE in
the list -- to stop a line near a cut from being averaged over the
trials that happened to keep it -- and not to add lines to it.  So
membership is decided by the unvaried run alone, which re-applies the
reporting cut exactly, that run having been made at it.  `min_bf`
still drops a line: it is a physics choice the user made, not
censoring the trials introduced.  See `write`.

THE TOTAL DECAY RATE COMES FROM RCG, NOT FROM ADDING THE LINES UP

The branching fraction is the line's A over the total radiative decay
rate of its upper level, and that total is RCG's own SUMGA rather than
a sum over the printed lines.  The reason is where SUMGA accumulates:
at label 580 of `rcg11kd.f`, which the DMIN test reaches by a `GO TO
578` that jumps past the *print* and into the *sum*.  SUMGA therefore
counts every line of the level, printed or not, and a branching
fraction taken from it is not censored by DMIN -- which matters most
exactly where DMIN is doing the most work.

Three things are worth knowing about it:

* It **is** censored by the FLBDN/FLBDX window, which drops a line at
  label 498 before it ever enters the arrays.  A branching fraction
  from a run with a window in force is a fraction of what is in the
  window, and the run says so once, in the log.
* It is printed as `E10.3` -- four significant figures.  Where nothing
  is censored the two agree to that: measured on WORK/Forbid, whose
  deck has DMIN = 0.0 and no window, the worst level differs from the
  sum over its lines by 4e-4.
* A field of the record can overflow into asterisks.  Where it is
  SUMGA itself, the whole level falls back to the sum over its printed
  lines and the log names it; where it is only a neighbouring SUMFI
  field, SUMGA is still good and is taken.  Nothing is substituted
  silently.  Only the Lund layout carries a SUMGA that can be tied to
  a level at all -- the wide layout names the level by its energy
  alone -- so a wide-layout run falls back throughout, with a warning.
* A level has a SUMGA only where it **heads** a printed block, and
  which levels do is ISPECC's business (column 72 of the main card):
  bit 1 prints the first-parity sort and bit 2 the second.  A deck
  punching 2 heads its blocks with the levels of one parity set only,
  and the other set falls back to the sum over its lines: measured on
  Pr III, 550 of 1046 decaying levels.  Punching **3** is the way out.
  Both sorts then print every line, and `read_transitions` folds the
  repeat -- the first print of a transition in one OUTG11 is the one
  that contributes, the rest are checked against it and dropped -- so
  one pass still contributes one value per transition and every level
  gets an uncensored total.  The run counts the levels that fell back
  and says so in the log.

SUMGA is not, however, the same quantity as the sum this program makes.
Where both levels of a line are observed this program quotes gA at the
OBSERVED wavelength (see the note at the end of this header), and RCG
summed SUMGA from its own calculated energies before this program saw
it.  Measured on Pr III at DMIN = 8.0, SUMGA over the sum of the
rescaled printed lines ran from 0.96 to 1.05.

So the total is neither one nor the other but

    total = sum_rescaled * SUMGA / sum_unrescaled

-- the sum of the printed lines on the reporting scale, times a factor
whose two terms are both on RCG's own scale and which is therefore
purely what DMIN removed.  `Lev` carries the two sums side by side for
it.  It is an approximation, and a deliberate one: the totals serve the
branching fractions, and the branching fractions serve as a cut-off for
deciding which lines are worth publishing, where a per cent either way
changes nothing (author, 2026-08-27).  Where there is no usable SUMGA
the rescaled sum is returned alone, which is all there is.

WARNINGS GO TO A FILE

A few hundred trials produce warnings by the hundred, and on the
terminal they bury the progress lines.  They go to `calc_TP.log` in the
working directory instead -- written and flushed as they happen, so a
run killed after six hours still says what it found -- and the terminal
gets one line at the end saying how many there were.  See `Log`.

WHY E1 IS NOT SIMPLY M1 WITH A DIFFERENT SECTION NAME

Two things change with it.  It crosses the two RCE parity sets, so both
have to be read and the two ends of a line looked up in different sets;
and it makes the line list large, so the printed list is censored -- by
DMIN, and by the FLBDN/FLBDX wavelength window -- in a way a forbidden
run never noticed.  A line that leaves the printed list in some trials
and not others has statistics built on the trials that kept it.

NUMPY IS REQUIRED

Alone in this package, `calc_TP` will not run without numpy: the
eigenvector matcher is a matrix product per J block and there is no
sensible pure-Python fallback for a case with five thousand of them.
A missing numpy stops the program with the `pip install numpy` line
rather than degrading.  `make_cfp`, where numpy is optional and only
true jj coupling is lost without it, is unaffected.  Measured on
Pr III, 5626 eigenvectors: `overlap.remap` over the whole run is 0.15 s
against 15.7 s for the trial around it.

Usage:

    calc_TP [folder] [--par <file>] [--write-par]
    calc_TP [folder] --cleanup [--yes]

The second form deletes the trial subdirectories of the working
directory and does nothing else.  They are kept by default: a trial can
be read again without being run again (`read_only`), and twenty
gigabytes of it is not this program's to throw away unasked.  For a run
long enough that the twenty gigabytes are the obstacle there is
`keep_trial_output = 0`, which deletes each trial's OUTG11 as soon as it
has been read and leaves the deck that made it; the disk then holds one
trial's output at a time rather than all of them.

There are no options for the physics.  Everything the Perl held in the
block of `my $...` at the top of the file -- the number of trials, the
variances, the cut-offs, the output file name -- is now read from

    calc_TP.par

in the working directory, or from the copy that ships beside this
script if the working directory has none.  `--write-par` drops a fresh,
fully commented copy of that file into the working directory; each
parameter carries its own explanation there, and the file is the
documentation.  See `toolbox/parfile.py` for the layout.

WHAT THE RUN NEEDS

  ING11, OUTG11   the unvaried RCG run, with the spectrum of every
                  multipole asked for, in the Lund layout (`ILUND`).
                  ING11 must carry the `//R1//` and `//R2//`
                  radial-integral cards if those integrals are to be
                  varied.
  RCEOUT, OUTE    the converged fit: the levels and their eigenvectors,
                  and the standard deviation of every parameter.
  LEVELS1         optional; the Lande g factors and the better energies.
  IN36            the configurations, for the level labels.

WHAT THE RUN WRITES

  out_file        the line list (`tp.txt` by default).
  out_levels      the level list its lid1/lid2 point into
                  (`tp_levels.txt` by default).
  calc_TP.log     the run's header, its warnings and its notes.

DIFFERENCES FROM THE PERL, ALL DELIBERATE

* **The varied parameters go straight into the deck.**  The Perl wrote
  a stand-in RCEOUT into each trial directory and shelled out to
  `s11.bat` to transfer them, then to `del *.bak`.  Here the deck is
  edited in memory and written out once, which removes the two `cmd`
  calls and the risk of a half-written RCEOUT being mistaken for the
  real one.
* **The trial count may be any number.**  The Perl drew its normal
  deviates in Box-Muller pairs and sent the two members of a pair to
  two consecutive trials, so the count had to be even.  See
  `toolbox/trials.py`.
* **`random_seed` makes a run repeatable**, which the Perl could not.
* **`fixed_param_variance` works.**  `calc_TP_M1E2.pl:582` tests
  `abs($flag == 100)`, which is `abs(($flag == 100))` -- 0 or 1 -- and a
  fixed parameter's flag is -100, so the branch never ran and the
  setting did nothing.  The test here is `abs(flag) == 100`, which is
  what the same file writes correctly four lines above.
* **The total decay rate of a trial is taken from that trial.**
  `calc_TP_M1E2.pl:838` averages `sum_A_hash[...][i]` over `i` from 0,
  which is the unvaried run's slot, and stops one short of the last
  trial.  The trials here are 1 to N throughout.
* **All the `//R2//` cards are varied**, not only the first
  uninterrupted run of them; see `ing11.read_r_integrals`.
* **Both parity sets can be done**, under `parity_sets`.  The Perl
  labelled every level it read as parity set 1.  E1 requires both.
* **A trial's eigenvectors are paired with the unvaried run's**, on
  overlap, rather than identified against RCE from scratch.  See
  `toolbox/overlap.py`: identification falls back to position above the
  highest observed level, and a trial moves the eigenvalues there.

ONE THING THAT IS KEPT, AND IS WORTH KNOWING

Where **both** levels of a line have observed energies, the line
strength S is held fixed and gA is recomputed at the observed
wavelength (`calc_TP_M1E2.pl:485-487`).  That is a rescaling, and
`conv_out.py` deliberately does not do it (author, 2026-08-15).  It is
kept here because it is what this program has always done and because
a forbidden line is where it matters most -- gA goes as the cube or the
fifth power of the wave number, and a small energy difference computed
from fitted levels can be wrong by a large fraction -- but it is a
question for the author, and it is written up in `FOR/LSJJ/README.txt`.
"""

import math
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import case as case_mod                               # noqa: E402
from toolbox import (files, ing11, levels, outg11, overlap,        # noqa: E402
                     parfile, radiate, rce, spectrum, trials, vacair)

PAR_FILE = 'calc_TP.par'
LOG_FILE = 'calc_TP.log'
CONFIG = 'COWAN.CFG'
RCG = 'rcg11kd.exe'

#: The same program built on Unix, where gfortran's -o leaves no
#: suffix.  Both names are looked for, .exe first.
RCG_UNIX = 'rcg11kd'

#: RCG's name for each multipole's spectrum section, and the order in
#: which the sections are asked for.  These are the three GAFROMSCOEF
#: writes; M2 and higher exist in `radiate` but RCG prints no spectrum
#: for them.
SPECTRUM_OF = {radiate.E1: 'ELEC DIP',
               radiate.M1: 'MAG DIP',
               radiate.E2: 'ELEC QUD'}

#: The multipoles that connect two levels of the SAME parity, and so of
#: one RCE parity set.  E1 is the odd one out: it crosses the two sets,
#: which is why it needs both of them read and why its levels are looked
#: up in different sets at the two ends of a line.
SAME_PARITY = (radiate.M1, radiate.E2)

#: The radial transition integral each multipole is built on -- the
#: //R1// cards for E1 and the //R2// cards for E2.  M1 has none: its
#: operator is L + 2S, which has no radial part to be uncertain about.
R_ORDER = {radiate.E1: 1, radiate.E2: 2}

#: A wavelength between these two wave numbers, in the run's own unit,
#: is quoted in standard air -- 2000 to 20000 angstroms.
AIR_FROM, AIR_TO = 5.0, 50.0

#: The lowest DMIN the trials are ever run at.  Still inside RCG's
#: gA range -- below 1.0 the field means something else entirely --
#: and 10**2 s-1 is no cut worth the name.
DMIN_FLOOR = 2.0

#: Two prints of one transition agree to better than this, relatively,
#: or the pairing is wrong and the log says so.  The printed gA is
#: `1PE10.3` -- four significant figures -- and the two level sorts run
#: the same arithmetic, so anything above a part in a thousand is not
#: rounding.
SAME_LINE = 1e-3

#: The largest FLBDX an `E5.0` field will hold: 99999 angstrom.  A trial
#: window widened past it is written at it instead, which is no cut at
#: all for any deck this package is used on.
FLBDX_FIELD_MAX = 99999.0

#: How many stored trial values in the written list before the run says
#: so.  Five million is a file of some tens of megabytes of trial
#: columns alone, which is where the list stops being a thing to open in
#: a spreadsheet.
TRIAL_CELL_WARN = 5000000

#: Bytes one stored trial value costs in the file: `%.5E` is eleven
#: characters and the tab makes twelve.
TRIAL_CELL_BYTES = 12

TAB = '\t'


class Usage(Exception):
    """The command line will not do what was asked."""


class RunError(Exception):
    """The run cannot be made, or cannot be trusted."""


class Log:
    """The run's warnings, gathered into one file.

    A run of a few hundred trials produces a warning per trial at least
    -- the eigenvectors that paired poorly, the lines a trial did not
    print -- and on the terminal they scroll past unread and bury the
    progress lines between them.  So they go to `calc_TP.log` in the
    working directory instead, written as they happen and flushed, so
    that a run killed after six hours still says what it had found; the
    terminal gets one line at the end saying how many there were.

    A note is written and not counted; a warning is written with its
    prefix and counted.  The file is overwritten by every run: it
    belongs to the run that wrote it and to no other.
    """

    def __init__(self, path):
        self.path = Path(path)
        self.count = 0
        self.fh = None

    def start(self, seed, folder):
        self.fh = open(self.path, 'w')
        self.fh.write('calc_TP %s -- %s\n'
                      % (folder, time.strftime('%Y-%m-%d %H:%M:%S')))
        self.fh.write('random seed: %r\n\n' % (seed,))
        self.fh.flush()

    def note(self, text):
        if self.fh is not None:
            self.fh.write(text + '\n')
            self.fh.flush()

    def warn(self, text):
        self.count += 1
        self.note('WARNING: ' + text)

    def finish(self):
        if self.fh is None:
            return
        self.fh.write('\n%d warning(s).\n' % self.count)
        self.fh.close()
        self.fh = None
        if self.count:
            print('%d warning(s); they are in %s' % (self.count, self.path))
        else:
            print('No warnings.  %s holds the run\'s header only.'
                  % self.path)


# --------------------------------------------------------------------------
# the parameter file
# --------------------------------------------------------------------------

S = parfile.Spec
SPECS = (
    S('num_trials', parfile.INT, 0, minimum=0, description="""
        Number of Monte-Carlo trials.  Each one is a complete RCG run in
        its own subdirectory trial1, trial2, ... and takes as long as the
        original did, so this is what sets the cost of the calculation.
        0 reads the unvaried run only and writes the line list with no
        statistics at all.  A few hundred is enough for a 1% answer.
        Unlike the Perl this need not be an even number."""),

    S('read_only', parfile.INT, 0, minimum=0, description="""
        Number of trial subdirectories that already hold a finished RCG
        run and are only to be read.  Use it to add trials to a run that
        was interrupted, or to re-do the statistics without re-running
        anything: read_only = num_trials reads every trial and runs RCG
        not at all.  Larger than num_trials is taken as num_trials."""),

    S('keep_trial_output', parfile.BOOL, True, description="""
        1 keeps each trial's OUTG11 after it has been read; 0 deletes
        it.  OUTG11 is by far the largest file a trial makes -- the
        Pr III list at the author's working DMIN = 3.0 is of the order
        of a gigabyte per trial, so a few hundred trials is tens of
        gigabytes -- and this program has no use for it once the trial
        has been folded into the statistics.  What is deleted is the
        output only: the trial's ING11 and its parameter record stay,
        so the trial can be re-run exactly as it was made.  It cannot
        be re-READ, though, so read_only can no longer reach it; keep
        the output if that is wanted.  Only the trials THIS run
        computed are deleted -- a trial that was there already and was
        read under read_only is left as it was found."""),

    S('store_trials', parfile.INT, 100, minimum=0, description="""
        How many individual trial values of A -- and of the cancellation
        factor -- are written out beside each transition.  The mean and
        the standard deviation always use every trial; this is only what
        is printed, so that the distribution can be looked at.  Larger
        than num_trials is taken as num_trials.  0 prints none."""),

    S('random_seed', parfile.TEXT, '', description="""
        Seed for the random number generator, so that a run can be
        repeated exactly.  Any text will do.  Blank draws a fresh seed
        from the system, which is what the Perl always did; the seed
        actually used is printed at the start of the run."""),

    S('multipoles', parfile.TEXT, 'E1,M1,E2', description="""
        Which multipoles to read from OUTG11, as a comma-separated list
        of E1, M1 and E2 in any order and any case.  A multipole RCG did
        not compute simply contributes no lines, so the default asks for
        all three.  E1 is what columns 41-42 of the ING11 main card
        switch on; M1 and E2 are IMAG and IQUAD, columns 49 and 50 of
        the G5INP card in IN2.  Note that E1 crosses the two parity sets
        and therefore needs both of them: asking for it sets
        parity_sets to `both` whatever the parameter says."""),

    S('r1_variance', parfile.FLOAT, 0.25, minimum=0.0, description="""
        Fractional standard deviation of the E1 radial transition
        integrals -- the //R1// cards at the end of ING11.  Same
        reasoning as r2_variance below, for the electric-dipole
        operator.  0 leaves them alone, and it is ignored when E1 is not
        among the multipoles."""),

    S('r2_variance', parfile.FLOAT, 0.25, minimum=0.0, description="""
        Fractional standard deviation of the E2 radial transition
        integrals -- the //R2// cards at the end of ING11.  These are
        Hartree-Fock values with no fitted uncertainty of their own, so
        this is an estimate of how far a single-configuration HF
        integral can be from the truth.  0 leaves them alone.

        0.25 is a scope-wide estimate of that error for a
        single-configuration Hartree-Fock integral, not a number this
        deck knows: it is the same for r1_variance and it is the first
        thing to revisit if the spread of the trials comes out plainly
        too wide or too narrow against a spectrum with measured
        transition probabilities."""),

    S('vary_fitted_params', parfile.BOOL, True, description="""
        1 to vary the RCE parameters that were free in the fit, each
        within the standard deviation OUTE gives it.  0 varies only the
        E2 integrals, which answers a different question: how much of
        the uncertainty comes from the radial integrals alone."""),

    S('var_grouped_params', parfile.BOOL, True, description="""
        1 to respect the fit's parameter grouping: parameters tied
        together in RCE by a shared group flag are displaced by one
        common ratio, because the fit determined their ratio and not
        their separate values.  0 displaces every parameter on its own,
        which overstates the freedom the fit actually had."""),

    S('fixed_param_variance', parfile.FLOAT, 0.0, minimum=0.0, description="""
        Fractional standard deviation applied to the parameters that
        were HELD FIXED in the fit (flag 100), which have no standard
        deviation of their own.  0 leaves them alone.  Note that this
        setting had no effect at all in the Perl -- see the header of
        this program."""),

    S('trial_dmin_offset', parfile.FLOAT, 1.0, minimum=0.0, description="""
        Decades below the deck's own DMIN at which the TRIALS are run.
        DMIN censors the printed list, and a line just above it drops
        out of some trials and not others -- which biases its mean high
        and gives it a trial count smaller than the number of trials.
        Running the trials one decade lower keeps those lines in every
        trial; the deck's DMIN is then applied again when the list is
        written, so the list is the one DMIN asked for and its
        statistics are complete.  0 runs the trials at the deck's own
        DMIN, which is what the Perl did.  Only whole decades have any
        effect: RCG takes the integer part of DMIN as the exponent."""),

    S('trial_lambda_margin', parfile.FLOAT, 0.05, minimum=0.0, description="""
        Fractional widening of the FLBDN/FLBDX wavelength window in the
        trial decks, for the same reason as trial_dmin_offset: a line
        near either end of the window leaves the printed list in some
        trials and not others.  0.05 runs the trials over a window 5 per
        cent wider at each end.  0 leaves the window alone, and it does
        nothing at all for a deck with no KCPL=3 card, which prints
        every wavelength anyway.

        It is the SHORT-wavelength end that is widened by this
        fraction.  The long end needs more, and how much more is worked
        out from the deck; see trial_lambda_margin_max."""),

    S('trial_lambda_margin_max', parfile.FLOAT, 10.0, minimum=0.0,
      description="""
        The largest fractional widening allowed at the LONG-wavelength
        end of the trial window.

        A level energy moves by a per cent or two between trials, but a
        transition energy is a difference of two of them, and between
        two highly excited configurations that difference is small -- so
        the same absolute movement is a much larger fraction of it.
        FLBDX cuts on exactly those smallest transition energies, so the
        symmetric trial_lambda_margin is far too small there.  The
        widening actually used is

            trial_lambda_margin * (Eav_up + Eav_low) / dEav_min

        with dEav_min the smallest gap between the average energies of
        two configurations the unvaried spectrum actually connects, and
        Eav_up and Eav_low the average energies of that same pair
        measured from the LOWEST configuration average energy in the run
        -- excitation energies, not Cowan's total electronic energies,
        which are dominated by a core binding no fitted parameter moves.
        That ratio has no bound of its own -- two configurations can sit
        on top of each other -- so this is the ceiling on it.  10.0
        means the long end may be widened elevenfold.  Written to the
        log, with the pair it came from and whether this ceiling bound
        it."""),

    S('merged_list', parfile.BOOL, True, description="""
        1 to write one line per transition, with the M1 and E2
        contributions added together and the column E2_frac saying what
        fraction of the total is E2.  0 writes the M1 list and then the
        E2 list separately, one after the other in the same file.
        (`$mearged_list` in the Perl.)"""),

    S('min_mixed_fraction', parfile.FLOAT, 0.02, minimum=0.0, maximum=1.0,
      description="""
        When merging, the smaller of the two multipoles is dropped from
        the line and the line keeps the name of the larger one if the
        smaller contributes less than this fraction of the total
        A-value, allowing two standard deviations.  Has no effect when
        merged_list = 0."""),

    S('min_bf', parfile.FLOAT, 0.0001, minimum=0.0, maximum=1.0,
      description="""
        Transitions with a radiative branching fraction smaller than
        this are left out of the list.  The branching fraction is the
        line's A over the total radiative decay rate of its upper level,
        counting only the channels this run holds -- that is, only the
        multipoles asked for, and only the lines RCG printed."""),

    S('out_file', parfile.TEXT, 'tp.txt', description="""
        Name of the tab-separated line list this program writes, in the
        working directory.  (`$out_tr_file` in the Perl, where it was a
        constant in the source and had to be edited for every ion.)"""),

    S('out_levels', parfile.TEXT, '', description="""
        Name of the tab-separated level list this program writes beside
        the line list, carrying the Lev_# the lid1 and lid2 columns of
        the line list refer to.  Blank derives it from out_file, which
        turns tp_E1.txt into tp_E1_levels.txt.  The two files are what
        insert_level_ids takes, so the identifiers of a published
        analysis can be put on these lines as they are put on
        conv_out's."""),

    S('print_dbf', parfile.BOOL, True, description="""
        1 to print the standard deviation of the branching fraction
        beside each line.  It needs the stored trials, so it is silent
        when store_trials = 0."""),

    S('print_cf_trials', parfile.BOOL, False, description="""
        1 to print the individual trial values of the cancellation
        factor as well as those of A.  This doubles the width of the
        file and is for looking at a suspect line, not for routine use."""),

    S('sort_trials', parfile.BOOL, True, description="""
        1 to print the stored trial values in ascending order, which is
        what a histogram or a percentile wants.  WARNING: sorting breaks
        the correspondence between the A values, the cancellation
        factors and the branching fractions of one and the same trial.
        0 keeps them in trial order."""),

    S('print_distribution', parfile.BOOL, False, description="""
        1 to print, on the terminal, histograms of how far the trials
        strayed from the unvaried run -- in ln A, in A, in the
        cancellation factor and in the input parameters -- together with
        the kurtosis of the input deviations.  This is a check that the
        sampling itself is sound: a normal distribution has kurtosis 0
        by Excel's definition, which is the one used."""),

    S('parity_sets', parfile.TEXT, '1', choices=('1', '2', 'both'),
      description="""
        Which RCE parity set to work on.  M1 and E2 connect levels of
        the same parity, so the two sets are independent and a run may
        do either or both.  E1 connects one set to the other and needs
        both, so asking for E1 among the multipoles overrides this to
        `both`.  The Perl could only do the first, and labelled
        everything it read as belonging to it."""),

    S('rcg', parfile.TEXT, '', description="""
        The RCG executable to run in each trial subdirectory.  Blank
        uses the one that sits beside this program, which is the
        arrangement the .bat files assume.  Not used when there are no
        trials to run (num_trials = read_only)."""),
)

TITLE = ('calc_TP.par -- input for calc_TP',
         '',
         'Uncertainties on E1, M1 and E2 transition probabilities, by',
         'repeating the RCG calculation with the fitted parameters and',
         'the radial transition integrals displaced within their',
         'errors.',
         '',
         'The quantity the trials are collected in is the LINE STRENGTH,',
         'not the transition probability: gA carries the eigenvalues of',
         'the trial in it, as the cube or the fifth power of the wave',
         'number, so its spread would be the matrix element and the',
         'energy together.  Where both levels are observed the wavelength',
         'does not move between trials and the two spreads are the same;',
         'they differ only when a level is unobserved and the fit lets',
         'its energy move.  See the header of calc_TP.py.')


# --------------------------------------------------------------------------
# the levels of one run
# --------------------------------------------------------------------------

class Lev:
    """One RCE level, with everything the line list asks of it."""

    __slots__ = ('one', 'kset', 'j', 'index', 'ar_tot', 'ar_raw', 'ga_tot')

    def __init__(self, one):
        self.one = one
        self.kset = one.kset
        self.j = one.j
        self.index = one.level.index       # its position in RCEOUT
        self.ar_tot = []                   # summed A per trial, 0 = unvaried
        self.ar_raw = []                   # the same sum, NOT rescaled to the
                                           # observed wavelength -- the scale
                                           # SUMGA is on.  See `total`.
        self.ga_tot = []                   # RCG's own SUMGA, None where none

    @property
    def key(self):
        return (self.kset, self.j, self.index)

    @property
    def observed(self):
        return self.one.level.observed

    @property
    def e_calc(self):
        return self.one.level.e_calc

    @property
    def e_obs(self):
        return self.one.level.e_obs

    @property
    def group(self):
        """The relative-energy group this level is in, or `''`."""
        return self.one.level.group

    @property
    def absolute(self):
        """Whether the observed energy is on the deck's absolute scale."""
        return self.one.level.absolute

    @property
    def g(self):
        return self.one.level.g

    def energy(self, other):
        """The energy to compare with another level's, in the same terms.

        Two energies are in the same terms where both are observed on
        one scale -- both absolute, or both in one relative-energy
        group.  Across a group boundary the offset is unknown and the
        fitted energies are the pair that can be compared.
        """
        return self.e_obs if rce.comparable(self, other) else self.e_calc

    def decay(self, trial, rate, raw):
        """Add one printed line to the sum over this level's decays.

        Twice over: `rate` is the line as this program quotes it --
        rescaled to the observed wavelength where both levels are
        observed -- and `raw` is the same line at RCG's own calculated
        energies, which is the scale RCG summed SUMGA on.  `total`
        needs both, one for the answer and one for the ratio that
        carries the censoring correction across.
        """
        while len(self.ar_tot) <= trial:
            self.ar_tot.append(0.0)
            self.ar_raw.append(0.0)
        self.ar_tot[trial] += rate
        self.ar_raw[trial] += raw

    def summed(self, trial, ga):
        """Add one section's SUMGA -- RCG's own total for this level.

        One per multipole: the sections are separate spectra and their
        totals add.  `None` marks the level as one whose SUMGA cannot
        be used at all, and it stays marked whatever else is added.
        """
        while len(self.ga_tot) <= trial:
            self.ga_tot.append(0.0)
        if ga is None:
            self.ga_tot[trial] = None
        elif self.ga_tot[trial] is not None:
            self.ga_tot[trial] += ga

    def total(self, trial):
        """The total radiative decay rate of this level in one trial.

        Not simply SUMGA.  SUMGA is RCG's own sum over the level's
        decays and it counts the lines DMIN kept out of the printed
        list -- that is the whole reason to read it -- but it is a sum
        at RCG's *calculated* energies, and the lines this program adds
        up are rescaled to the observed wavelength wherever both ends
        are observed.  The two are therefore not the same quantity, and
        on Pr III at DMIN = 8.0 they differed by up to 5 per cent.

        So the rescaling is kept and the censoring is taken from SUMGA
        as a ratio:

            total = sum_rescaled * SUMGA / sum_unrescaled

        where the two sums run over the same printed lines.  The factor
        `SUMGA / sum_unrescaled` is then purely what DMIN removed, both
        of its terms being on RCG's own scale, and it multiplies a sum
        that is on the reporting scale.

        This is an approximation, and it is one deliberately: the
        totals serve the branching fractions, and the branching
        fractions serve as a cut-off for deciding which lines are worth
        publishing, where a per cent either way changes nothing
        (author, 2026-08-27).

        Where there is no usable SUMGA -- no record, or a starred field
        -- or where the un-rescaled sum is zero, the rescaled sum is
        returned alone, which is all there is.  See `spectrum.totals`.
        """
        rescaled = self.ar_tot[trial] if trial < len(self.ar_tot) else 0.0
        if trial < len(self.ga_tot) and self.ga_tot[trial] is not None:
            sumga = self.ga_tot[trial] / (2 * self.j + 1)
            raw = self.ar_raw[trial] if trial < len(self.ar_raw) else 0.0
            if raw > 0:
                return rescaled * sumga / raw
            return sumga
        return rescaled


class Transition:
    """One pair of levels and one multipole, over every trial."""

    __slots__ = ('upper', 'lower', 'multipole', 'delta_e', 'wavelength',
                 'ga', 's', 'cf', 'log_s', 'cancel', 'a_values', 'cf_values',
                 'n_seen')

    def __init__(self, upper, lower, multipole):
        self.upper, self.lower = upper, lower
        self.multipole = multipole
        self.delta_e = 0.0            # from the fitted energies
        self.wavelength = 0.0         # vacuum, from the fitted energies
        self.ga = None                # the unvaried run's, None if it had none
        self.s = None
        self.cf = None
        self.log_s = trials.Accumulator()
        self.cancel = trials.Accumulator()
        self.a_values = []
        self.cf_values = []
        self.n_seen = 0               # trials this line was printed in

    @property
    def key(self):
        return (self.upper.key, self.lower.key, self.multipole)

    @property
    def a(self):
        """The unvaried run's transition probability, or None."""
        if self.ga is None:
            return None
        return self.ga / (2 * self.upper.j + 1)

    @property
    def lam(self):
        """The wavelength this line's S is turned back into gA at.

        The observed one where both levels are observed on one scale,
        the unvaried run's fitted one otherwise -- the same choice
        `row_text` prints in `wl_c`, and the one that makes the mean A
        the mean S expressed as a rate.
        """
        if rce.comparable(self.upper, self.lower):
            spread = abs(self.upper.e_obs - self.lower.e_obs)
            if spread:
                return 1e5 / spread
        return self.wavelength

    @property
    def s_mean(self):
        """The geometric mean of S over the trials -- exp(mean ln S)."""
        return math.exp(self.log_s.mean) if self.log_s.n else 0.0

    @property
    def a_mean(self):
        """The mean S written as a transition probability.

        S is what the trials are collected in; A follows from it at one
        fixed wavelength, so the two means say the same thing and the
        one printed is the one the reader wants.
        """
        lam = self.lam
        if not lam or not self.log_s.n:
            return 0.0
        return (radiate.ga_from_s(self.s_mean, lam, self.multipole)
                / (2 * self.upper.j + 1))

    @property
    def spread(self):
        """The relative standard deviation of S, in per cent."""
        return (math.exp(self.log_s.sigma) - 1.0) * 100.0


# --------------------------------------------------------------------------
# reading one OUTG11
# --------------------------------------------------------------------------

def _kset_of(deck):
    return {card.label: card.kset for card in deck}


def parse_multipoles(text):
    """The `multipoles` parameter as an ordered tuple of radiate names.

    The order is the one written, because it is the order the sections
    of the line list come out in when merged_list is 0.
    """
    out = []
    for word in text.replace(',', ' ').split():
        name = word.upper()
        if name not in SPECTRUM_OF:
            raise RunError(
                'multipoles names %r, which is not one of E1, M1, E2.' % word)
        if name not in out:
            out.append(name)
    if not out:
        raise RunError('multipoles is empty; name at least one of E1, M1, E2.')
    return tuple(out)


def read_transitions(folder, mapping, of_label, wanted, store, trial,
                     found, table, order, multipoles, pairs=None, log=None):
    """Fold one OUTG11's spectrum into the running statistics.

    `found` is `{(kset, J, RCE index): Lev}` and `table` the transitions
    keyed by `Transition.key`; both are built on the unvaried run and
    added to by every trial.  A transition a trial has and the unvaried
    run had not is added with no unvaried values, which is what the
    `gA is None` tests downstream mean.

    **A trial may put the two levels the other way round.**  Which end
    is upper is decided by energy *in that run*, and a level the fit did
    not observe moves between trials, so two levels of nearly equal
    energy can change places.  The same physical line then arrives under
    the reversed key.  It is looked up under the ordered key first and
    under the reversed one after, and a hit on the reversed key
    contributes to the transition as it stands; only a miss on both
    makes a new one.  The `E_low`/`E_up` columns keep the unvaried run's
    ordering -- they describe the calculation that is being published --
    and only the trials are allowed to disagree about which end is
    which.  `decay` is left alone: the decay rate belongs to whichever
    level is upper in that trial, which is what RCG computed.

    **A run may print the same transition more than once.**  ISPECC
    (column 72 of the main card) chooses the level sorts, and a deck
    asking for both prints every line in each of them.  That is worth
    having -- only a level that HEADS a printed block gets a SUMGA, so
    both sorts is what gives every level an uncensored total -- so the
    repeat is folded rather than refused: the first print of a
    transition in one OUTG11 is the one that contributes, and every
    later print of the same unordered pair is checked against it and
    otherwise ignored.  One OUTG11 still contributes exactly one value
    per transition.  A disagreement beyond `SAME_LINE` is a warning
    naming the transition, not a silent choice: the two sorts print the
    same computed number, and a difference means the pairing is wrong.

    `pairs`, where given, collects `((kset, config), (kset, config))` for
    every line read -- the configuration pairs the spectrum actually
    connects, which is what `censor` measures its long-wavelength
    widening on.

    Returns `(reversed_matches, repeated_prints)` for the log.
    """
    seen = {}
    reversed_hits = 0
    repeated = 0
    mismatched = 0
    sections = tuple(SPECTRUM_OF[name] for name in multipoles)
    for line in spectrum.lines(folder=folder, multipoles=sections):
        multipole = radiate.OF_SPECTRUM[line.multipole]
        kset = of_label.get(line.initial.label)
        if kset is None:
            raise RunError(
                'the OUTG11 spectrum of %s names the configuration %r, which '
                'is not on any ING11 card' % (folder, line.initial.label))
        if kset not in wanted:
            continue

        # RCG prints each E1 pair once, and the initial level of that
        # one line is always the one in the second parity set: the two
        # ends of the line are therefore looked up in different sets.
        # M1 and E2 stay inside one set, where `other` is `kset` and
        # this is the lookup the program has always made.
        other = kset if multipole in SAME_PARITY else 3 - kset
        if other not in wanted:
            continue
        first = mapping.get((kset, line.initial.j, line.initial.level))
        second = mapping.get((other, line.final.j, line.final.level))
        if first is None or second is None:
            continue

        a = found.get((first.kset, first.j, first.level.index))
        b = found.get((second.kset, second.j, second.level.index))
        if a is None:
            a = found.setdefault((first.kset, first.j, first.level.index),
                                 Lev(first))
        if b is None:
            b = found.setdefault((second.kset, second.j, second.level.index),
                                 Lev(second))
        if a is b:
            raise RunError(
                'the %s spectrum of %s has a line from level %d (set %d, '
                'J=%s) to itself' % (multipole, folder, a.index, a.kset, a.j))

        # The upper level is the one higher in energy, on the observed
        # energies where both are observed and on the fitted ones
        # otherwise (`calc_TP_M1E2.pl:453`).
        upper, lower = (a, b) if a.energy(b) >= b.energy(a) else (b, a)

        difference = abs(upper.e_calc - lower.e_calc)
        if not difference:
            continue
        fitted = 1e5 / difference
        # gA as RCG printed it, undone back to the line strength at the
        # wavelength RCG printed it AT -- this trial's own calculated
        # energies.  S is what the trial contributes; gA is then this
        # trial's S at the reporting wavelength, which is the observed
        # one where both levels are observed on one scale (the
        # rescaling `calc_TP_M1E2.pl:485-487` has always done) and the
        # trial's own otherwise.
        s = radiate.s_from_ga(line.ga, _printed(line), multipole)
        ga = line.ga
        if rce.comparable(upper, lower):
            spread = abs(upper.e_obs - lower.e_obs)
            if spread:
                ga = radiate.ga_from_s(s, 1e5 / spread, multipole)

        rate = ga / (2 * upper.j + 1)

        # The unordered key: the pair, whichever way round this run put
        # it.  Both the duplicate test and the trial-to-base lookup are
        # made on it, or a reversal would read as a duplicate print and
        # a duplicate print as a new line.
        first_key, second_key = upper.key, lower.key
        unordered = (min(first_key, second_key), max(first_key, second_key),
                     multipole)
        if unordered in seen:
            # A second level sort printing the same line again.  It
            # carries the same computed number, so it is checked and
            # dropped; adding it would count this run twice.
            repeated += 1
            before = seen[unordered]
            if before > 0 and s > 0 and abs(s - before) > SAME_LINE * max(s,
                                                                         before):
                mismatched += 1
                if log is not None and mismatched == 1:
                    log.warn(
                        '%s: the %s transition between level %d (set %d, '
                        'J=%s) and level %d (set %d, J=%s) is printed twice '
                        'with line strengths differing by %.2g per cent; the '
                        'first print is the one used.  The two level sorts '
                        'print the same computed number, so a difference '
                        'means the levels were paired wrongly.'
                        % (Path(folder).name, multipole, upper.index,
                           upper.kset, upper.j, lower.index, lower.kset,
                           lower.j, 100.0 * abs(s - before) / max(s, before)))
            continue
        seen[unordered] = s

        upper.decay(trial, rate, line.ga / (2 * upper.j + 1))
        if pairs is not None:
            pairs.add(((kset, line.initial.config), (other, line.final.config)))

        key = (first_key, second_key, multipole)
        one = table.get(key)
        if one is None:
            # This run may have swapped the two ends; the transition is
            # then in the table under the reversed key, and that is the
            # entry this line belongs to.
            one = table.get((second_key, first_key, multipole))
            if one is not None:
                reversed_hits += 1
        if one is None:
            one = Transition(upper, lower, multipole)
            one.delta_e = difference
            one.wavelength = fitted
            table[key] = one
            order.append(one)
            if not trial:
                one.ga, one.s, one.cf = ga, s, line.cf
        if trial:
            one.n_seen += 1
        if s > 0:
            # S, not A, is what the trials are collected in.  gA as RCG
            # prints it carries the trial's own eigenvalues in it -- a
            # level the fit did not observe moves, and its lines move
            # with the cube or the fifth power of the wave number -- so
            # the spread of A over the trials would be the spread of the
            # matrix element and the spread of the energy together.  S
            # is the matrix element alone.  A transition with no
            # strength at all has no logarithm and is left out of the
            # mean rather than stopping the run, which is what the
            # Perl's log() would do.
            one.log_s.add(math.log(s))
        one.cancel.add(line.cf)
        if trial and trial <= store:
            one.a_values.append((trial, rate))
            one.cf_values.append((trial, line.cf))
    return reversed_hits, repeated


def read_totals(folder, mapping, of_label, wanted, trial, found,
                multipoles, log):
    """Fold one OUTG11's SUMGA records into the levels of `found`.

    RCG prints its own summed decay rate at the end of every `* * *`
    block, and it is that, not the sum over the printed lines, that the
    branching fractions are taken from: SUMGA accumulates at label 580
    of `rcg11kd.f`, which the DMIN test jumps into rather than past, so
    it counts the lines DMIN kept out of the list.  What it does not
    count is a line the FLBDN/FLBDX window dropped -- that one never
    reaches the arrays at all (label 498) -- so a total from a run with
    a window in force is a total over what is in the window, and `run`
    says so once.

    A record whose SUMGA field overflowed into asterisks marks its
    level: the whole level then falls back to the sum over its printed
    lines, because a total missing one multipole is worse than a total
    that is merely censored.  A star in a neighbouring SUMFI field
    leaves SUMGA readable and is taken.

    Levels the run has no record for keep whatever `decay` summed.
    """
    sections = tuple(SPECTRUM_OF[name] for name in multipoles)
    seen = set()
    starred = []
    for one in spectrum.totals(folder=folder, multipoles=sections):
        kset = of_label.get(one.level.label)
        if kset is None or kset not in wanted:
            continue
        where = mapping.get((kset, one.level.j, one.level.level))
        if where is None:
            continue
        key = (where.kset, where.j, where.level.index)
        level = found.get(key)
        if level is None:
            continue
        # The two level sorts print the same level's total twice when
        # the multipole connects a parity set to itself.  The two
        # records carry the same sum; adding them would double it.
        if (one.multipole, key) in seen:
            continue
        seen.add((one.multipole, key))
        level.summed(trial, one.sumga)
        if one.sumga is None:
            starred.append((one.multipole, level))
    if starred:
        log.warn('%s: the SUMGA field overflowed for %d level(s); their '
                 'total decay rate is the sum over the printed lines '
                 'instead, which DMIN censors.  First: %s level %d '
                 '(set %d, J=%s).'
                 % (Path(folder).name, len(starred), starred[0][0],
                    starred[0][1].index, starred[0][1].kset,
                    starred[0][1].j))
    return len(seen)


def _printed(line):
    """The wavelength RCG printed, recomputed where the field overflowed."""
    if line.wavelength is not None and line.wavelength > 0:
        return line.wavelength
    return 1e5 / line.delta_e if line.delta_e else 0.0


# --------------------------------------------------------------------------
# building the trial decks
# --------------------------------------------------------------------------

def parameter_slots(deck, kset):
    """`[(card, position)]` in the order RCE prints its parameter table.

    Every single-configuration card first, in deck order, then every CI
    card -- the walk `s11.py` makes, and the one RCG asks for them in.
    """
    out = []
    for card in deck.parameters[kset - 1]:
        out.extend((card, i) for i in range(len(card.values)))
    for card in deck.ci[kset - 1]:
        out.extend((card, i) for i in range(len(card.values)))
    return out


def check_parameters(deck, csets, oute, kset):
    """That the deck, RCEOUT and OUTE describe one and the same fit.

    The deck's parameter fields carry no names, so they are matched to
    RCEOUT by position and count; OUTE is matched to RCEOUT by name,
    position for position, because OUTE is the only file with a standard
    deviation in it and RCEOUT the only one that ties a parameter to a
    configuration.  Everything downstream then indexes all three by the
    same integer.
    """
    cset = None
    for one in csets:
        if one.index == kset:
            cset = one
    if cset is None:
        raise RunError('RCEOUT has no parity set %d' % kset)
    fitted = oute.get(kset)
    if not fitted:
        raise RunError(
            'OUTE has no parameter table for parity set %d.  Without it '
            'there is no standard deviation to vary a parameter within; run '
            'RCE and keep its OUTE, or set vary_fitted_params to 0.' % kset)
    if len(fitted) != len(cset.parameters):
        raise RunError(
            'parity set %d: OUTE lists %d parameters and RCEOUT %d.  The two '
            'files are from different runs.'
            % (kset, len(fitted), len(cset.parameters)))
    for at, (a, b) in enumerate(zip(fitted, cset.parameters), 1):
        if a.name != b.name:
            raise RunError(
                'parity set %d, parameter %d: OUTE calls it %r and RCEOUT '
                '%r.  The two files are from different runs.'
                % (kset, at, a.name, b.name))
    slots = parameter_slots(deck, kset)
    if len(slots) != len(fitted):
        raise RunError(
            'parity set %d: ING11 has room for %d parameters and the fit has '
            '%d.  Re-run RCG from the deck the fit was made on.'
            % (kset, len(slots), len(fitted)))
    return cset, fitted, slots


def vary_parameters(fitted, options, rng, histogram):
    """New values for one parity set's parameters.

    A free parameter is displaced within its own standard deviation.
    Where `var_grouped_params` is on, parameters that share a group flag
    are displaced by one common ratio instead: the fit determined their
    ratio, not their separate values, so displacing them independently
    would grant the fit a freedom it did not have.
    """
    out = []
    ratios = {}
    for one in fitted:
        value, sigma, flag = one.value, one.sigma or 0.0, one.flag
        new = value
        if options['vary_fitted_params'] and sigma:
            if options['var_grouped_params'] and abs(flag) != 100 and value:
                ratio = ratios.get(flag)
                if ratio is None:
                    ratio = rng.normal(1.0, abs(sigma / value))
                    ratios[flag] = ratio
                    histogram.add((ratio - 1.0) / abs(sigma / value),
                                  keep=True)
                new = value * ratio
            else:
                new = rng.normal(value, abs(sigma))
                histogram.add((new - value) / sigma, keep=True)
        elif (options['fixed_param_variance'] and value
                and abs(flag) == 100):
            # The branch `calc_TP_M1E2.pl:582` could never reach.
            width = abs(value * options['fixed_param_variance'])
            new = rng.normal(value, width)
            histogram.add((new / value - 1.0)
                          / options['fixed_param_variance'], keep=True)
        out.append(new)
    return out


def vary_r_integrals(cards, fraction, rng, histogram):
    """New values for one order's radial transition integrals.

    `cards` are the //R1// or the //R2// cards of the deck and
    `fraction` the variance parameter that goes with that order.
    """
    out = []
    for card in cards:
        new = rng.normal(card.value, abs(fraction * card.value))
        if fraction and card.value:
            histogram.add((new / card.value - 1.0) / fraction, keep=True)
        out.append(new)
    return out


def write_trial(deck, folder, cards, values, per_set, ending):
    """Write one trial's ING11.  The deck is restored before returning.

    `cards` are every radial-integral card being varied -- the //R1//
    and the //R2// cards together -- and `values` their new values in
    the same order.
    """
    was_radial = [card.value for card in cards]
    was_par = {}
    for kset, (slots, fitted) in per_set.items():
        for (card, at), value in zip(slots, fitted):
            was_par.setdefault((id(card), at), (card, at, card.values[at]))
            card.values[at] = value
            card.edited = True
    for card, value in zip(cards, values):
        ing11.set_r_value(card, value)

    files.write_deck(Path(folder) / files.ING11, ing11.write_ing11(deck),
                     ending)

    for card, at, value in was_par.values():
        card.values[at] = value
        card.edited = False
    for card, value in zip(cards, was_radial):
        ing11.set_r_value(card, value)


def long_margin(pairs, configs, margin, ceiling, log=None):
    """The fractional widening of FLBDX the trial decks need.

    `trial_lambda_margin` is right at the short-wavelength end and badly
    wrong at the long one, and the reason is that a level energy and a
    transition energy do not move by the same fraction.  A level moves
    by a per cent or two; a transition energy is a *difference*, and
    between two highly excited configurations that difference is small,
    so the same absolute movement is a far larger fraction of it -- and
    it is the smallest transition energies in the run that FLBDX cuts
    on.

    Author's estimate, adopted (2026-08-27).  Take `dEav_min`, the
    smallest absolute difference of the two configuration average
    energies over the configuration PAIRS THAT ACTUALLY CARRY A LINE of
    the unvaried spectrum inside the window -- `pairs`, collected while
    that spectrum was read -- and take `Eav_up` and `Eav_low` from that
    same pair.  Then

        f_long = margin * (Eav_up + Eav_low) / dEav_min

    which is the fractional movement of the smallest transition energy
    the run has, if each end moves by `margin` of its own energy.  The
    author's phrasing is "even and odd", which is the E1 case; the rule
    here is the general one, the minimum running over whatever pairs the
    spectrum connects -- same-parity ones for M1 and E2.

    **The energies are EXCITATION energies**, measured from the lowest
    configuration average energy in the run.  `outg11`'s Eav is Cowan's
    own, the total electronic energy of the configuration: on Pr III it
    is about -2.03e6 kK, of which 350 kK is the whole spread of the
    thirty-nine configurations and the rest is the binding of 57
    electrons that no fitted parameter moves.  Taken literally the
    formula would then be
    `0.05 x 4.05e6 / 3.07 = 65000` and every deck would sit on the
    ceiling.  What "a level moves by a per cent or two" means is a per
    cent or two of its energy ABOVE THE GROUND CONFIGURATION, which is
    the number the fit is about, so that is what is used: the same
    Pr III pair gives 0.05 x 183.5 / 3.07 = 3.0.

    `dEav_min` can be tiny or zero, two configurations of nearly equal
    average energy being perfectly possible, so `f_long` is bounded by
    `ceiling` (`trial_lambda_margin_max`).  Returns `(fraction, pair)`,
    with `margin` itself and no pair where there is nothing to measure
    on -- an empty spectrum, or configurations the OUTG11 table printed
    no energy for.
    """
    # The zero: the lowest configuration average energy the run has.  A
    # configuration under a term cutoff prints neither name nor energy
    # (`outg11.ConfigCard`) and takes no part in it.
    energies = [card.eav for card in configs.values() if card.eav]
    if not energies:
        return margin, None
    ground = min(energies)

    best = None
    for (ka, ca), (kb, cb) in pairs:
        one, two = configs.get((ka, ca)), configs.get((kb, cb))
        if one is None or two is None or not one.eav or not two.eav:
            continue
        gap = abs(one.eav - two.eav)
        if best is None or gap < best[0]:
            up, low = (one, two) if one.eav >= two.eav else (two, one)
            best = (gap, up, low)
    if best is None:
        return margin, None
    gap, up, low = best
    total = (up.eav - ground) + (low.eav - ground)
    bound = gap <= 0 or margin * total > ceiling * gap
    fraction = ceiling if bound else margin * total / gap
    if log is not None:
        log.note('long-wavelength trial margin: %g (%g x (%g + %g) / %g, '
                 'excitation energies above the lowest configuration '
                 'average %g)%s.  The closest pair of configurations the '
                 'spectrum connects is set %d configuration %d and set %d '
                 'configuration %d.'
                 % (fraction, margin, up.eav - ground, low.eav - ground, gap,
                    ground, ' -- bound by trial_lambda_margin_max'
                    if bound else '',
                    up.kset, up.number, low.kset, low.number))
    return fraction, best


def censor(deck, options, multipoles, ntrials, pairs=(), configs=None,
           log=None):
    """Set the cuts the trial decks are run under; return the reporting cut.

    Two thresholds censor RCG's printed line list, and neither was ever
    a systematic the Perl noticed, because a forbidden run punched
    neither: DMIN, columns 61-65 of the main card, and the FLBDN/FLBDX
    wavelength window of the KCPL=3 card.  A line just above either of
    them is printed in some trials and not in others, and a mean taken
    over the trials that kept it is a mean over a censored sample --
    biased high, and with a trial count smaller than the number of
    trials.  Measured on Pr III at DMIN=8.0, 9462 lines, 5 trials: 1407
    lines were seen in 1 to 4 trials and 1773 more appeared only in
    trials.

    So the trials are run in a window wider than the one being reported
    -- `trial_dmin_offset` decades lower and `trial_lambda_margin`
    wider at each end.  Every line of the written list then has the
    statistics of every trial, and the lines the wider window brought in
    are dropped again at the end: the list holds the unvaried run's
    transitions and only those, which `write` decides by membership
    rather than by re-applying this threshold.

    Only the trial decks are touched.  The unvaried run is the author's
    own and is read as it stands; it is what the list is taken from.

    E0MAX is deliberately left alone.  It is a scope choice -- the
    highest level the run is about -- and not a threshold that noise
    pushes a line across: a level near it moves by its own uncertainty,
    which is a fraction of a per cent of the energy, so the set of
    levels E0MAX admits is the same in every trial.

    Returns `(kind, threshold)` for the cut to re-apply at write time,
    or None when the trials run under the deck's own cuts and there is
    nothing to re-apply.
    """
    if not ntrials:
        return None
    dmin = deck.main.dmin
    report = None

    if radiate.E1 not in multipoles:
        # Author's decision, 2026-08-27.  With E1 in the run the
        # forbidden multipoles share the spectrum, and the DMIN chosen
        # for E1 is the one that matters.  A forbidden-only run has few
        # enough lines to want all of them, so the trials are run with
        # no cut at all.
        if dmin:
            ing11.set_dmin(deck, 0.0)
            print('   E1 is not among the multipoles, so the trials are run '
                  'with DMIN = 0.0 in place of the punched %g.' % dmin)
            if log is not None:
                log.warn('the unvaried run was made at DMIN = %g and the '
                         'trials are run at 0.0, so the trials hold lines '
                         'the list cannot: the list is the unvaried run.'
                         % dmin)
    elif dmin and options['trial_dmin_offset']:
        kind, threshold = ing11.dmin_cut(dmin)
        if kind != 'gA':
            # Below 1.0 the field cuts on D0, the squared reduced matrix
            # element, which is not a quantity this program has; it
            # cannot be re-applied at write time, so it is not lowered.
            print('   DMIN = %g cuts on the matrix element rather than on '
                  'gA, so trial_dmin_offset is ignored.' % dmin)
        else:
            decades = round(options['trial_dmin_offset'])
            if abs(options['trial_dmin_offset'] - decades) > 1e-9:
                print('   trial_dmin_offset = %g rounded to %d: RCG takes the '
                      'integer part of DMIN as the exponent.'
                      % (options['trial_dmin_offset'], decades))
            if decades:
                lowered = max(abs(dmin) - decades, DMIN_FLOOR)
                # A negative DMIN also deletes upward transitions from
                # the level sorts (NEGTP); that is not a cut and is kept.
                ing11.set_dmin(deck, -lowered if dmin < 0 else lowered)
                print('   trials run at DMIN = %g (gA >= %.0E); the list is '
                      'written at the punched DMIN = %g (gA >= %.0E).'
                      % (lowered, ing11.dmin_cut(lowered)[1], dmin, threshold))
                report = (kind, threshold)

    margin = options['trial_lambda_margin']
    if margin and deck.control(3) is not None:
        flbdn, flbdx = deck.window
        if (flbdn, flbdx) != (ing11.FLBDN_DEFAULT, ing11.FLBDX_DEFAULT):
            # The two ends are not widened by the same fraction: see
            # `long_margin` for why FLBDX needs far more than FLBDN.
            far, _ = long_margin(pairs, configs or {}, margin,
                                 options['trial_lambda_margin_max'], log)
            widened = min(flbdx * (1.0 + far), FLBDX_FIELD_MAX)
            if widened < flbdx * (1.0 + far) and log is not None:
                log.warn('the trial FLBDX wanted %g A and the E5.0 field '
                         'holds at most %g; it is written at the maximum, '
                         'which cuts nothing an optical deck cares about.'
                         % (flbdx * (1.0 + far), FLBDX_FIELD_MAX))
            ing11.set_window(deck, flbdn / (1.0 + margin), widened)
            print('   trials run over %g to %g A (long end widened by %g); '
                  'the list is written over %g to %g A.'
                  % (deck.window + (far, flbdn, flbdx)))
    return report


def run_rcg(folder, executable):
    """Run RCG in one trial directory, the way rcg.bat does.

    RCG finds its permanent files -- SENIOR, the cfp tapes, LSJJ.DAT --
    through COWAN.CFG, which holds the directory they are in and which
    the batch file writes fresh in the working directory on every run
    (`FOR/rcg11kd.f:337-351`).  The same is done here, and the file is
    removed again afterwards, so a trial directory is left holding only
    what RCG produced.

    MSDOS is set to 1 in the source (`FOR/rcg11kd.f:313`), so the
    branch that reads COWAN.CFG is taken on every platform and the
    ELSE branch's compiled-in path is dead.  RCG cuts the line back
    to its last backslash if it has one and otherwise keeps it whole,
    repairing nothing: the separator this appends is what makes the
    path usable, on Unix as much as on Windows.
    """
    folder = Path(folder)
    config = folder / CONFIG
    with open(config, 'w') as fh:
        fh.write(str(Path(executable).resolve().parent) + os.sep + '\n')
    try:
        done = subprocess.run([str(executable)], cwd=str(folder))
    finally:
        if config.exists():
            config.unlink()
    if done.returncode:
        raise RunError('RCG failed in %s (exit code %d)'
                       % (folder, done.returncode))


def discard_output(folder, log):
    """Delete one trial's OUTG11; return the bytes it freed.

    `keep_trial_output = 0`.  The trial has been read by the time this
    is called and nothing else in the program looks at it again, so the
    file is dead weight -- and it is the one file in the trial that is
    large: at the author's working DMIN the Pr III spectrum is of the
    order of a gigabyte a trial, which is what makes a few hundred
    trials impractical to keep.

    Only OUTG11 goes.  ING11 and the parameter record stay behind, so
    the trial is a complete description of itself and can be run again;
    what it cannot be is re-READ, so `read_only` can no longer reach it.
    A failure to delete is a note, not an error: the statistics are
    already in hand and there is no reason to lose a run over a locked
    file.
    """
    path = Path(folder) / files.OUTG11
    try:
        size = path.stat().st_size
        path.unlink()
    except OSError as problem:
        log.warn('%s could not be deleted (%s); keep_trial_output = 0 has '
                 'not freed it.' % (path, problem))
        return 0
    return size


def find_rcg(named):
    if named:
        path = Path(named)
        if not path.is_file():
            raise RunError('the rcg parameter names %s, which is not there'
                           % path)
        return path
    for name in (RCG, RCG_UNIX):
        beside = files.code_dir() / name
        if beside.is_file():
            return beside
    found = shutil.which(RCG_UNIX) or shutil.which(RCG)
    if found:
        return Path(found)
    raise RunError(
        'RCG was not found beside %s nor on the PATH.  Name it with the rcg '
        'parameter of %s.' % (files.code_dir(), PAR_FILE))


# --------------------------------------------------------------------------
# merging M1 with E2
# --------------------------------------------------------------------------

class Row:
    """One printed line: one pair of levels, one or both multipoles."""

    __slots__ = ('upper', 'lower', 'kind', 'ga', 's', 'cf', 'ga_mean',
                 's_mean', 'spread', 'cf_mean', 'cf_sigma', 'e2_fraction',
                 'e2_sigma', 'a_values', 'cf_values', 'wavelength', 'n_seen',
                 'in_base')

    def __init__(self, one):
        self.upper, self.lower = one.upper, one.lower
        self.kind = one.multipole
        self.ga, self.s, self.cf = one.ga, one.s, one.cf
        # Whether the UNVARIED run printed this transition.  It is what
        # decides whether the row is written at all: see `write`.
        self.in_base = one.ga is not None
        self.ga_mean = one.a_mean * (2 * one.upper.j + 1)
        self.s_mean = one.s_mean if one.log_s.n else None
        self.spread = one.spread
        self.cf_mean = one.cancel.mean
        self.cf_sigma = one.cancel.sigma
        self.e2_fraction = 0.0 if one.multipole == radiate.M1 else 1.0
        self.e2_sigma = 0.0
        self.a_values = dict(one.a_values)
        self.cf_values = dict(one.cf_values)
        self.wavelength = one.lam
        self.n_seen = one.n_seen

    @property
    def strength(self):
        """The gA the branching fraction is taken from."""
        return self.ga if self.ga is not None else self.ga_mean


def merge(rows, ntrials, minimum):
    """Add the M1 and the E2 line of one level pair into one row.

    The smaller multipole is dropped, and the row keeps the name of the
    larger, when it contributes less than `minimum` of the total
    A-value allowing two standard deviations; otherwise the two are
    added and the row is called M1+E2.  This is
    `calc_TP_M1E2.pl:640-770`, with the level pair used as the key
    directly: the Perl keyed on the wavelength rounded to three decimals
    and then compared the level pairs anyway, appending a '1' to the key
    when they differed.
    """
    m1 = {(r.upper.key, r.lower.key): r for r in rows
          if r.kind == radiate.M1}
    out = []
    for row in rows:
        if row.kind == radiate.M1:
            continue
        other = m1.pop((row.upper.key, row.lower.key), None)
        if other is None:
            out.append(row)
            continue
        out.append(_combine(row, other, ntrials, minimum))
    out.extend(m1.values())
    return out


def _combine(e2, m1, ntrials, minimum):
    """One E2 row and the M1 row of the same level pair."""
    ga_e2, ga_m1 = e2.strength, m1.strength
    total = ga_e2 + ga_m1
    if total <= 0:
        return e2
    fraction = ga_e2 / total

    # The two spreads add in quadrature, weighted by the two strengths,
    # in the logarithm -- which is where the spread is normal.
    log_e2 = math.log(e2.spread / 100.0 + 1.0)
    log_m1 = math.log(m1.spread / 100.0 + 1.0)
    r = ga_m1 / ga_e2 if ga_e2 else 0.0
    spread = math.sqrt(log_e2 * log_e2 + r * r * log_m1 * log_m1) / (1.0 + r)

    both = Row.__new__(Row)
    both.upper, both.lower = e2.upper, e2.lower
    both.kind = 'M1+E2'
    # The level pair is in the unvaried run if either multipole of it
    # was, whichever name the merged row ends up under.
    both.in_base = e2.in_base or m1.in_base
    both.wavelength = e2.wavelength
    both.spread = (math.exp(spread) - 1.0) * 100.0
    both.ga = None if e2.ga is None else total
    both.ga_mean = e2.ga_mean + m1.ga_mean
    both.s = None
    # Two multipoles have no common line strength -- their S are in
    # different units of the operator -- so the merged row carries none,
    # and its trial count is the larger of the two: the row was printed
    # in at least that many trials.
    both.s_mean = None
    both.n_seen = max(e2.n_seen, m1.n_seen)
    both.cf = (None if e2.cf is None or m1.cf is None
               else (abs(e2.cf) * ga_e2 + abs(m1.cf) * ga_m1) / total)
    both.cf_sigma = math.sqrt((e2.cf_sigma * ga_e2 / total) ** 2
                              + (m1.cf_sigma * ga_m1 / total) ** 2)
    both.cf_mean = (abs(e2.cf_mean) * ga_e2 + abs(m1.cf_mean) * ga_m1) / total

    # Per trial: the two A values add, and the cancellation factors are
    # averaged over them.  The fraction of E2 is followed trial by trial
    # so that its own spread can be reported.
    a_values, cf_values = {}, {}
    fractions = trials.Accumulator()
    for at in sorted(set(e2.a_values) | set(m1.a_values)):
        a_e2 = e2.a_values.get(at, 0.0)
        a_m1 = m1.a_values.get(at, 0.0)
        if at in e2.a_values and at in m1.a_values and a_e2 + a_m1:
            fractions.add(a_e2 / (a_e2 + a_m1))
        a_values[at] = a_e2 + a_m1
        if a_values[at]:
            cf_values[at] = (abs(e2.cf_values.get(at, 0.0)) * a_e2
                             + abs(m1.cf_values.get(at, 0.0)) * a_m1) \
                / a_values[at]
        else:
            cf_values[at] = 0.0
    both.a_values, both.cf_values = a_values, cf_values

    mean = fractions.mean if fractions.n else fraction
    offset = fraction - mean
    variance = fractions.sigma ** 2 if fractions.n else 0.0
    both.e2_fraction = fraction
    both.e2_sigma = math.sqrt(variance + offset * offset)

    if fraction + 2 * both.e2_sigma < minimum:
        m1.e2_fraction, m1.e2_sigma = fraction, both.e2_sigma
        m1.n_seen = both.n_seen
        m1.in_base = both.in_base
        return m1
    if 1.0 - fraction + 2 * both.e2_sigma < minimum:
        # The line is E2 to within the cut-off.  The Perl keeps the
        # summed gA under the name E2 and leaves S at the pure E2 value;
        # that is kept, and is one of the questions in README.txt.
        e2.ga, e2.ga_mean = both.ga, both.ga_mean
        e2.cf, e2.cf_mean, e2.cf_sigma = both.cf, both.cf_mean, both.cf_sigma
        e2.spread = both.spread
        e2.a_values, e2.cf_values = a_values, cf_values
        e2.e2_fraction, e2.e2_sigma = fraction, both.e2_sigma
        e2.n_seen = both.n_seen
        e2.in_base = both.in_base
        return e2
    return both


# --------------------------------------------------------------------------
# the report
# --------------------------------------------------------------------------

def header(options):
    # E1 and E2 used to head the two level energies.  Both are now
    # multipole names as well, so the four energy columns say which end
    # of the line they belong to instead.
    # lid1 and lid2 sit between the calculated energies and the
    # wavelength, where conv_out puts them, and mean the same thing: the
    # Lev_# of the lower and of the upper level in the level list this
    # run writes beside the line list.
    out = ['type', 'conf1', 't1', 'J1', 'conf2', 't2', 'J2', 'lande1',
           'lande2', 'E_low', 'E_up', 'lid1', 'lid2', 'wl_c(A)',
           'E_low_exp', 'E_up_exp', 'wl_exp(A)', 'A(s-1)', 'S', 'cf',
           'sum A', 'BF']
    if options['merged_list']:
        out.append('E2_frac')
        if options['num_trials']:
            out.append('d_E2frac')
    if options['num_trials']:
        out.extend(['n_trials', 'd_BF', 'A_mean', 'S_mean', 'std_A%',
                    'cf_mean', 'std_cf'])
    if options['store_trials']:
        out.append('Trial A data')
        out.extend('%d' % n for n in range(2, options['store_trials'] + 1))
        if options['print_cf_trials']:
            out.append('Trial cf data')
            out.extend('%d' % n for n in range(2, options['store_trials'] + 1))
    return out


def _air(wavelength, difference):
    if AIR_FROM < difference < AIR_TO:
        return vacair.to_air(wavelength)
    return wavelength


def row_text(row, options, name, number):
    """One printed line, or None where it falls under the cut-offs.

    Columns 1 and 2 are the LOWER and the UPPER level in that order,
    which is what the Perl prints under the names conf1 and conf2.
    """
    upper, lower = row.upper, row.lower
    at = 0                              # the unvaried run's slot
    total = upper.total(at)
    if not total:
        # No decay at all in the unvaried run -- average the trials.
        kept = [x for x in (upper.total(n) for n in
                            range(1, max(len(upper.ar_tot),
                                         len(upper.ga_tot)))) if x > 0]
        total = sum(kept) / len(kept) if kept else 0.0
    branch = (row.strength / (2 * upper.j + 1) / total) if total else 1.0
    if branch < options['min_bf']:
        return None

    difference = abs(upper.e_calc - lower.e_calc)
    if not difference:
        return None
    wavelength = 1e5 / difference
    experimental = ''
    if rce.comparable(upper, lower):
        spread = abs(upper.e_obs - lower.e_obs)
        if spread:
            difference = spread
            experimental = '%.4f' % _air(1e5 / spread, spread)
    row.wavelength = 1e5 / difference     # vacuum, from the better energies
    wavelength = '%.4f' % _air(wavelength, difference)

    out = [row.kind, name(lower)[0], name(lower)[1], '%.1f' % lower.j,
           name(upper)[0], name(upper)[1], '%.1f' % upper.j,
           _g(lower), _g(upper),
           '%.4f' % lower.e_calc, '%.4f' % upper.e_calc,
           number(lower), number(upper), wavelength,
           '%.4f' % lower.e_obs if lower.observed else '',
           '%.4f' % upper.e_obs if upper.observed else '',
           experimental,
           '' if row.ga is None else '%.5E' % (row.ga / (2 * upper.j + 1)),
           '' if row.s is None else '%.4E' % row.s,
           '' if row.cf is None else '%.15g' % row.cf,
           '%.5E' % total, '%.15g' % branch]
    if options['merged_list']:
        out.append('%.15g' % row.e2_fraction)
        if options['num_trials']:
            out.append('' if row.n_seen < 2 else '%.15g' % row.e2_sigma)
    if options['num_trials']:
        out.extend(_statistics(row, upper, options, branch))
    if options['store_trials']:
        out.extend(_trial_columns(row.a_values, options, '%.5E'))
        if options['print_cf_trials']:
            out.extend(_trial_columns(row.cf_values, options, '%.7E'))
    return out


def _g(level):
    return '' if level.g is None else '%.3f' % level.g


def _statistics(row, upper, options, branch):
    """n_trials, d_BF, A_mean, S_mean, std_A%, cf_mean, std_cf.

    The spread quoted is the spread of the trials about the unvaried
    run, not about their own mean: the offset between the two is a real
    part of the uncertainty and is added in quadrature, in the
    logarithm for S and directly for the cancellation factor.

    A line the trials printed once, or not at all, gets a BLANK standard
    uncertainty in every column that carries one -- std_A%, std_cf and
    d_BF.  One sample has no spread, and a zero there would read as a
    line that is known exactly.  `n_trials` is printed beside them
    because it is what tells the two apart, and because a line seen in
    fewer trials than were run was censored out of the rest (DMIN, or
    the wavelength window) and its statistics are those of the trials
    that kept it.
    """
    spread, cf_sigma = row.spread, row.cf_sigma
    reference = row.strength
    if reference > 0 and row.ga_mean > 0:
        log = math.log(spread / 100.0 + 1.0)
        offset = math.log(row.ga_mean / reference)
        spread = (math.exp(math.sqrt(log * log + offset * offset)) - 1.0) * 100
        start = row.cf if row.cf is not None else row.cf_mean
        d_cf = start - row.cf_mean
        cf_sigma = math.sqrt(cf_sigma * cf_sigma + d_cf * d_cf)

    d_bf = ''
    if options['print_dbf'] and options['store_trials'] and row.n_seen > 1:
        fractions = []
        for at in range(1, options['store_trials'] + 1):
            a = row.a_values.get(at)
            if a is None:
                continue
            total = upper.total(at)
            if total:
                fractions.append(a / total)
        found = trials.stdev(branch, fractions)
        if found is not None:
            d_bf = '%.15g' % found

    a_mean = row.ga_mean / (2 * upper.j + 1)
    # S_mean is the mean of the tracked quantity itself.  A_mean beside
    # it is that same number at the VACUUM wavelength of the better
    # energies; the wl_c column may be in air
    # (`calc_TP_M1E2.pl:826-833`).
    s_mean = '' if row.s_mean is None else '%.4E' % row.s_mean
    enough = row.n_seen > 1
    return ['%d' % row.n_seen, d_bf, '%.5E' % a_mean, s_mean,
            '%.2f' % spread if enough else '',
            '%.7f' % row.cf_mean, '%.7f' % cf_sigma if enough else '']


def _trial_columns(values, options, shape):
    """The stored trials, in order or sorted, padded out to store_trials."""
    kept = [values[at] for at in sorted(values) if at]
    if options['sort_trials']:
        kept.sort()
    out = [shape % x for x in kept[:options['store_trials']]]
    return out + [''] * (options['store_trials'] - len(out))


# --------------------------------------------------------------------------
# putting the run together
# --------------------------------------------------------------------------

def run(folder, options, out, log):
    folder = Path(folder)
    deck = ing11.read_ing11(folder=folder)
    ending = files.line_ending(files.ING11, folder)
    of_label = _kset_of(deck)

    multipoles = parse_multipoles(options['multipoles'])
    asked = options['parity_sets'].lower()
    both = tuple(k for k in (1, 2) if deck.main.nconfigurations(k))

    if radiate.E1 in multipoles and len(both) < 2:
        # A one-parity deck has no E1 transition at all, so this is not
        # an error: the default multipoles ask for all three and a
        # forbidden-line deck must still run under them.
        print('   the deck has configurations in one parity set only, so '
              'there are no E1 transitions; E1 is dropped.')
        multipoles = tuple(n for n in multipoles if n != radiate.E1)
        if not multipoles:
            raise RunError(
                '%s has configurations in one parity set only and E1 was the '
                'only multipole asked for, so there is nothing to do.'
                % (folder / files.ING11))
    if radiate.E1 in multipoles and asked != 'both':
        # E1 has one level in each set, so a run restricted to one set
        # would identify one end of every line and drop it.
        print('   E1 was asked for, so both parity sets are read '
              '(parity_sets = %s ignored).' % options['parity_sets'])
        asked = 'both'

    wanted = {'1': (1,), '2': (2,), 'both': (1, 2)}[asked]
    wanted = tuple(k for k in wanted if k in both)
    if not wanted:
        raise RunError('%s has no configuration in parity set %s'
                       % (folder / files.ING11, options['parity_sets']))

    if options['merged_list'] and not (radiate.M1 in multipoles
                                       and radiate.E2 in multipoles):
        # Merging is the M1/E2 arrangement and nothing else: one row per
        # level pair with E2_frac saying how the two share it.  E1 does
        # not belong in that sum -- it connects the other parity set --
        # and a single forbidden multipole has nothing to merge with.
        # So the merge, and the E2_frac column with it, is off unless
        # the run is exactly M1 and E2.
        print('   the run is not M1 and E2 together, so the list is written '
              'one section per multipole (merged_list = 1 ignored).')
        options['merged_list'] = False

    ntrials = options['num_trials']
    already = min(options['read_only'], ntrials)
    store = min(options['store_trials'], ntrials)
    options['store_trials'] = store

    # The seed is drawn here, before anything is read, so that the log
    # can carry it in its header; no deviate is taken until the trials.
    rng = trials.Randomizer(options['random_seed'] or None)
    log.start(rng.seed, folder)
    if ntrials:
        print('Random seed: %r' % (rng.seed,))

    print('Reading the unvaried run...')
    base = case_mod.Case.read(folder, couplings=('LS',), use_emax=True)
    if base.unmatched:
        log.warn('%d fitted level(s) above E0MAX have no RCG eigenvector; '
                 'they carry no line.' % len(base.unmatched))
    mapping = {(one.kset, one.j, one.rcg): one for one in base.identified}

    radial = []                      # [(order, fraction, [cards])]
    for name in multipoles:
        order = R_ORDER.get(name)
        if order is None:                # M1 has no radial integral
            continue
        fraction = options['r%d_variance' % order]
        cards = ing11.read_r_integrals(deck, order=order)
        radial.append((order, fraction, cards))
        print('   %d %s radial integrals (//R%d// cards)'
              % (len(cards), name, order))
        if ntrials and fraction and not cards:
            log.warn('r%d_variance is set but the deck has no //R%d// card; '
                     'the %s line strengths will not move with it.'
                     % (order, order, name))
    r_cards = [card for _, _, cards in radial for card in cards]

    slots = {}
    if ntrials and (options['vary_fitted_params']
                    or options['fixed_param_variance']):
        oute = rce.read_oute(folder=folder)
        for kset in wanted:
            _, fitted, positions = check_parameters(deck, base.csets, oute,
                                                    kset)
            slots[kset] = (positions, fitted)
        print('   %d fitted parameters'
              % sum(len(v[0]) for v in slots.values()))

    found, table, order = {}, {}, []
    pairs = set()
    print('Reading transitions...')
    _, repeated = read_transitions(folder, mapping, of_label, wanted, store, 0,
                                   found, table, order, multipoles,
                                   pairs=pairs, log=log)
    print('   %d transitions' % len(order))
    # The stored trials are a rectangle: one column per trial for every
    # transition, twice that with print_cf_trials.  At the author's
    # working DMIN the list is of order 100000 lines, so a few hundred
    # stored trials is a file of gigabytes -- worth saying before the
    # trials are run rather than after.
    cells = store * len(order) * (2 if options['print_cf_trials'] else 1)
    if cells > TRIAL_CELL_WARN:
        log.warn('store_trials = %d over %d transition(s)%s is %.1f million '
                 'trial values in the list, about %.1f GB of tab-separated '
                 'text on top of the rest of the file.  Lower store_trials '
                 'if that is not wanted: the mean and the standard '
                 'deviation use every trial either way, and this is only '
                 'what is printed.'
                 % (store, len(order),
                    ' and the cancellation factors'
                    if options['print_cf_trials'] else '',
                    cells / 1e6, cells * TRIAL_CELL_BYTES / 1e9))
    if repeated:
        # A deck punching both level sorts (ISPECC=3) prints every line
        # in each of them.  That is worth having -- it is what gives
        # every level a SUMGA -- and the repeats are folded, one value
        # per transition either way.
        print('   %d repeated print(s) folded (both level sorts)' % repeated)
    # RCG's own totals, which the branching fractions are taken from.
    # They are read after the lines because a level enters `found` only
    # when a line of it was read, and a level with no line has no
    # branching fraction to want a total for.
    totals = read_totals(folder, mapping, of_label, wanted, 0, found,
                         multipoles, log)
    print('   %d total decay rates (SUMGA)' % totals)
    # A level gets a SUMGA record only where it HEADS a `* * *` block,
    # and which levels do that is ISPECC's business (column 72 of the
    # main card): bit 1 prints the first-parity sort and bit 2 the
    # second.  A deck punching 2 heads its blocks with the levels of one
    # set only, and the rest fall back to the sum over their printed
    # lines.  Punching 3 IS the way out: both sorts then print every
    # line, `read_transitions` folds the repeat, and every level has an
    # uncensored total.
    short = [one for one in found.values()
             if one.ar_tot and one.ar_tot[0] > 0
             and not (one.ga_tot and one.ga_tot[0] is not None)]
    if totals and short:
        log.warn('%d of %d decaying level(s) head no printed block of the '
                 'spectrum, so RCG printed no SUMGA for them and their '
                 'total decay rate is the sum over the printed lines, '
                 'which DMIN censors.  Which levels head a block is set by '
                 'ISPECC, column 72 of the ING11 main card: punch 3 for '
                 'both level sorts and every level gets a total.  The '
                 'repeated prints that come with it are folded.'
                 % (len(short), sum(1 for one in found.values()
                                    if one.ar_tot and one.ar_tot[0] > 0)))
    if not totals:
        log.warn('the unvaried run printed no SUMGA record this program '
                 'could read, so every total radiative decay rate is the '
                 'sum over the printed lines and DMIN censors it.  SUMGA '
                 'is printed by the Lund layout (ILUND), which is what '
                 'this program reads for its CF column anyway.')
    if deck.control(3) is not None and deck.window != (ing11.FLBDN_DEFAULT,
                                                       ing11.FLBDX_DEFAULT):
        log.warn('the deck punches a wavelength window of %g to %g A.  A '
                 'line outside it never enters RCG\'s arrays (label 498), '
                 'so SUMGA does not count it and every branching fraction '
                 'here is a fraction of what is inside the window.'
                 % deck.window)

    histograms = {name: trials.Histogram(name) for name in
                  ('parameters', 'lnA', 'A', 'cf')}

    # The configuration table, for the long-wavelength trial margin.
    # `read_outg11` stops at the end of that table, so this costs
    # nothing beside the spectrum just read.
    configs = outg11.read_outg11(folder=folder).configs
    censor(deck, options, multipoles, ntrials, pairs=pairs, configs=configs,
           log=log)

    executable = find_rcg(options['rcg']) if ntrials > already else None
    freed = 0
    for trial in range(1, ntrials + 1):
        where = folder / ('trial%d' % trial)
        where.mkdir(exist_ok=True)
        if trial > already:
            print('Trial %d of %d...' % (trial, ntrials))
            per_set = {kset: (positions,
                              vary_parameters(fitted, options, rng,
                                              histograms['parameters']))
                       for kset, (positions, fitted) in slots.items()}
            r_values = [value
                        for _, fraction, cards in radial
                        for value in vary_r_integrals(
                            cards, fraction, rng, histograms['parameters'])]
            write_trial(deck, where, r_cards, r_values, per_set, ending)
            for name in ('LEVELS1', 'LEVELS2'):
                source = files.resolve(name, folder, required=False)
                if source is not None:
                    shutil.copyfile(source, where / name)
            run_rcg(where, executable)
        else:
            print('Trial %d of %d (already computed)...' % (trial, ntrials))
        # The trial's eigenvectors are paired with the unvaried run's by
        # overlap rather than identified from scratch; `toolbox/overlap.py`
        # says why `levels.identify` cannot be used on a run whose
        # eigenvalues have moved.  A trial that moved nothing pairs one to
        # one and costs one matrix product per J block, which is nothing
        # beside the RCG run it follows, so this is not made conditional.
        blocks = levels.by_block(outg11.read_outg11(folder=where), ['LS'])
        trial_map, report = overlap.remap(mapping, base.blocks, blocks)
        log.note('trial %d: %s' % (trial, report))
        swapped, again = read_transitions(where, trial_map, of_label, wanted,
                                          store, trial, found, table, order,
                                          multipoles, log=log)
        if swapped or again:
            # A reversal means two levels of nearly equal energy changed
            # places when the fit moved them, which is worth knowing;
            # the line still counts towards the base transition.
            log.note('trial %d: %d line(s) matched with the two levels the '
                     'other way round, %d repeated print(s) folded.'
                     % (trial, swapped, again))
        read_totals(where, trial_map, of_label, wanted, trial, found,
                    multipoles, log)
        if not options['keep_trial_output'] and trial > already:
            # Everything this program wants of the trial is now folded
            # into the statistics, and OUTG11 is what makes a long run
            # cost tens of gigabytes.  Only a trial THIS run computed is
            # deleted: one that was there already was the user's file
            # before it was ours.  The ING11 and the parameter record
            # stay, so the trial can be re-run exactly as it was made.
            freed += discard_output(where, log)

    if freed:
        print('%.1f GB of trial output deleted (keep_trial_output = 0).'
              % (freed / 1e9))

    # `write` prints one section per multipole when the list is not
    # merged, and it must be the multipoles that survived above.
    options['multipoles'] = multipoles
    write(out, base, order, options, ntrials, log)
    if options['print_distribution']:
        distribution(order, options, histograms)


def write(out, base, order, options, ntrials, log):
    """The line list itself.

    **What is in it.**  The unvaried run's transitions, every one of
    them and no others (author, 2026-08-27).  The trials are run in a
    window wider than the deck punched, so a trial holds lines the
    unvaried run does not; those are not results.  The widened window
    is there to complete the statistics of the lines that ARE results,
    and membership is therefore decided by the unvaried run alone --
    which also re-applies the reporting cut exactly, the unvaried run
    having been made at it.  `min_bf` still drops a line: that is a
    physics choice the user made and not censoring the trials brought
    in.

    A level is named after the largest component of its LS eigenvector
    in the **unvaried** run, which is what `get_leading_LS_term`
    (`conv_cowan.pl:3242`) does.  A level that only a trial has a
    transition for carries that trial's eigenvector number, so the base
    run's own numbering is used wherever it has one: the identification
    is against the same RCE levels either way, but the name must not
    depend on which trial happened to see the line first.
    """
    names = {}
    numbers = {(one.kset, one.j, one.level.index): one.rcg
               for one in base.identified}

    def name(level):
        if level.key not in names:
            rcg = numbers.get(level.key, level.one.rcg)
            weights = base.amplitudes('LS', level.kset, level.j, rcg)
            found = base.designation('LS', level.kset, level.j,
                                     weights[0][1])
            names[level.key] = (found.configuration, found.term)
        return names[level.key]

    # Every level of the unvaried run is numbered, whether or not a line
    # survived to it, and the rule is `levels.numbered` -- the rule
    # conv_out numbers by.  So the two programs' Lev_# agree wherever
    # they are run on the same output and conv_out's upper_bound admits
    # every level, and the level list written below is what turns a lid
    # into a level.  A level only a trial saw is not in the base run's
    # ranking and its lid is left blank.
    ranking = levels.numbered(base)
    lids = {(one.kset, one.j, one.level.index): at
            for at, (one, _) in enumerate(ranking, 1)}

    def number(level):
        at = lids.get(level.key)
        return '' if at is None else '%d' % at

    groups = [''] if options['merged_list'] else list(options['multipoles'])
    trial_only = short = written = 0
    with open(out, 'w', newline='') as fh:
        fh.write(TAB.join(header(options)) + '\n')
        for at, kind in enumerate(groups):
            rows = [Row(one) for one in order
                    if not kind or one.multipole == kind]
            if options['merged_list']:
                rows = merge(rows, ntrials, options['min_mixed_fraction'])
            printed = []
            for row in rows:
                # A transition no trial would have shown had the trials
                # been run in the deck's own window.  See the docstring.
                if not row.in_base:
                    trial_only += 1
                    continue
                text = row_text(row, options, name, number)
                if text is not None:
                    printed.append((row.wavelength, row.upper.key,
                                    row.lower.key, text))
                    written += 1
                    short += 0 < row.n_seen < ntrials
            printed.sort(key=lambda x: x[:3])
            for _, _, _, text in printed:
                fh.write(TAB.join(text) + '\n')
    print('%s: %d level pairs written.' % (out, written))
    write_levels(levels_path(out, options), base, ranking)
    if trial_only:
        log.note('%d transition(s) appeared only in trials, the trials '
                 'having been run in a wider window than the deck punched; '
                 'they are not in the list.' % trial_only)
    if short:
        log.warn('%d written line(s) were printed in fewer than the %d '
                 'trials; their statistics are those of the trials that '
                 'kept them.  Raise trial_dmin_offset or '
                 'trial_lambda_margin.' % (short, ntrials))


def levels_path(out, options):
    """Where the level list goes: `out_levels`, or beside `out_file`."""
    named = options['out_levels'].strip()
    if named:
        return Path(out).parent / named
    out = Path(out)
    return out.with_name(out.stem + '_levels' + out.suffix)


def write_levels(out, base, ranking):
    """The level list the `lid1` and `lid2` columns point into.

    The columns conv_out writes, as far as the leading component:
    a `Lev_#` is only useful beside the level it names, and
    `insert_level_ids` reads `Lev_#`, `Ee`, `e_c`, `J` and `par` from a
    file of exactly this shape.  Running it on this pair of files puts
    the identifiers of a published analysis on these lines, which is the
    only way to join them to a conv_out list made from another
    calculation: a Lev_# is a position in one run's ordering and two
    runs whose eigenvalues have moved do not order their near-degenerate
    levels alike.

    Only the leading component is written.  The full purity is
    `printout`'s business and this file exists to be joined, not read.
    """
    parity = {}
    with open(out, 'w', newline='') as fh:
        fh.write(TAB.join(['Cpl', 'par', 'Lev_#', 'Ee', 'e_c', 'Ec', 'J',
                           'Lande_g', '%1', 'conf1', 'term1']) + '\n')
        for at, (one, _) in enumerate(ranking, 1):
            if one.kset not in parity:
                parity[one.kset] = ('e' if base.card(one.kset, 1).parity == 0
                                    else 'o')
            level = one.level
            weights = base.amplitudes('LS', one.kset, one.j, one.rcg)
            found = base.designation('LS', one.kset, one.j, weights[0][1])
            fh.write(TAB.join([
                'LS', parity[one.kset], '%d' % at,
                '%12.6f' % (level.e_obs if level.observed else level.e_calc),
                '' if level.observed else '*',
                '%12.4f' % level.e_calc, '%4.1f' % one.j,
                '' if level.g is None else '%.3f' % level.g,
                '%d' % int(weights[0][0] * weights[0][0] * 100 + 0.5),
                found.configuration, found.term]) + '\n')
    print('%s: %d levels written.' % (out, len(ranking)))


def distribution(order, options, histograms):
    """The histograms of how far the trials strayed, on the terminal."""
    for one in order:
        reference = one.a
        # Normalized by the spread of ln S, which is the spread the
        # sampling produced; where both levels are observed the two
        # logarithms differ by a constant and this is exactly ln A.
        if reference is None or reference <= 0 or not one.log_s.sigma:
            continue
        for _, a in one.a_values:
            if a > 0:
                histograms['lnA'].add(math.log(a / reference)
                                      / one.log_s.sigma, keep=True)
        if one.cancel.sigma:
            for _, cf in one.cf_values:
                histograms['cf'].add((cf - (one.cf or 0.0))
                                     / one.cancel.sigma)
        width = math.sqrt(sum((a - reference) ** 2 for _, a in one.a_values)
                          / len(one.a_values)) if one.a_values else 0.0
        if width:
            for _, a in one.a_values:
                histograms['A'].add((a - reference) / width)

    for name, title in (('lnA', 'logarithmic deviations'),
                        ('A', 'straight deviations'),
                        ('cf', 'CF deviations'),
                        ('parameters', 'input parameter deviations')):
        one = histograms[name]
        if not one:
            continue
        print('\nDistribution of %s:' % title)
        if name == 'parameters':
            kurtosis = trials.kurtosis(0.0, one.values)
            print('Kurtosis = %s' % ('' if kurtosis is None
                                     else '%.4f' % kurtosis))
        print('d/stdev\tN')
        for deviation, count in one.rows():
            print('%g\t%d' % (deviation, count))


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

def cleanup(folder, assume_yes):
    """Delete the trial subdirectories of one working directory.

    A run of a few hundred trials at the author's working DMIN leaves
    twenty gigabytes of trial output behind, and it is kept by default
    -- a trial can be re-read (`read_only`) or looked at, and deleting
    it is not this program's business to do on its own.  So the removal
    is a command-line flag that does nothing else: no par file is read,
    no RCG is run, no list is written.

    Only `trial*` directories directly inside the working directory go,
    a symbolic link among them is refused rather than followed, and
    nothing outside the working directory is touched.
    """
    folder = Path(folder)
    if not folder.is_dir():
        raise RunError('%s is not a directory' % folder)
    found = sorted(one for one in folder.glob('trial*')
                   if one.is_dir() and not one.is_symlink())
    linked = sorted(one.name for one in folder.glob('trial*')
                    if one.is_symlink())
    if linked:
        raise RunError(
            '%s holds symbolic link(s) named like trial directories (%s).  '
            'They are not followed and not removed; move them aside and run '
            '--cleanup again.' % (folder, ', '.join(linked)))
    if not found:
        raise RunError('%s holds no trial directory, so there is nothing '
                       'to clean up.' % folder)
    print('This will delete, from %s:' % folder)
    for one in found:
        print('   %s%s' % (one.name, os.sep))
    if not assume_yes:
        try:
            answer = input('%d directories.  Delete them? [y/N] '
                           % len(found))
        except EOFError:
            answer = ''
        if answer.strip().lower() not in ('y', 'yes'):
            print('Nothing was removed.')
            return 0
    for one in found:
        shutil.rmtree(one)
    print('%d trial directories removed.' % len(found))
    return 0


def main(argv):
    folder, named, fresh = '.', None, False
    clean = assume_yes = False
    rest = list(argv)
    while rest:
        arg = rest.pop(0)
        if arg in ('--write-par', '--write_par'):
            fresh = True
        elif arg == '--cleanup':
            clean = True
        elif arg in ('--yes', '-y'):
            assume_yes = True
        elif arg in ('--par', '--par-file'):
            if not rest:
                raise Usage('--par needs a file name')
            named = rest.pop(0)
        elif arg.startswith('--par='):
            named = arg.partition('=')[2]
        elif arg.startswith('-'):
            raise Usage('%r is not an option of this program' % arg)
        elif folder == '.':
            folder = arg
        else:
            raise Usage('%r is a second folder; only one may be given' % arg)

    if clean:
        if fresh or named:
            raise Usage('--cleanup does nothing but delete the trial '
                        'directories; it takes no other option')
        return cleanup(folder, assume_yes)
    if assume_yes:
        raise Usage('--yes answers the question --cleanup asks, and only '
                    'that; it does nothing on its own')

    if fresh:
        path = parfile.write(Path(folder) / PAR_FILE, SPECS, title=TITLE)
        print('%s written.  Read it through, then run calc_TP again.'
              % path)
        return 0

    # `parfile.find` is the shared rule -- the named file, the working
    # directory, then the copy beside the code, which this program does
    # ship.  Absent everywhere it is written rather than demanded: there
    # is nothing useful to run without it, and a run on the defaults
    # would spend hours computing the wrong thing.
    path = parfile.find(PAR_FILE, folder, named)
    if path is None:
        path = parfile.write(Path(folder) / PAR_FILE, SPECS, title=TITLE)
        print('There is no %s here and none beside %s, so a fresh one '
              'has been written:' % (PAR_FILE, files.code_dir()))
        print('  %s' % path)
        print('Every parameter is explained in it.  Set what this run '
              'needs, then run calc_TP again.')
        return 0
    print('Parameters from %s' % path)
    options = parfile.read(path, SPECS)
    log = Log(Path(folder) / LOG_FILE)
    try:
        run(folder, options, Path(folder) / options['out_file'], log)
    finally:
        log.finish()
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('calc_TP: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1].split(
                             'WHAT THE RUN NEEDS')[0]))
        sys.exit(2)
    except (RunError, parfile.ParFileError, case_mod.CaseError,
            ing11.Ing11Error, rce.RceError, levels.LevelsError,
            spectrum.SpectrumError, outg11.Outg11Error,
            OSError) as problem:
        sys.stderr.write('calc_TP: %s\n' % problem)
        sys.exit(1)
