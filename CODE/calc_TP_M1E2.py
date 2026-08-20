#!/usr/bin/env python3
"""Uncertainties on forbidden-line transition probabilities.

The Python replacement for `calc_TP_M1E2.pl`, and the last of the Perl
utilities.  It answers the question Cowan's codes do not: *how well is
this M1 or E2 transition probability known?*

There is no analytic answer, so it is obtained by resampling.  RCE has
fitted the Slater parameters to the observed levels and printed a
standard deviation for each of them in OUTE; RCN2 has computed the E2
radial integrals, whose accuracy is a fractional guess.  This program
displaces all of them at random within those uncertainties, runs RCG
again in a trial subdirectory, re-identifies the eigenvectors with the
fitted levels, and collects the transition probability and the
cancellation factor of every line.  Over many trials the spread of
ln A is the uncertainty, and it is reported as a percentage, which is
what a log-normal distribution's relative width is.

The result is one tab-separated line list -- every M1 and E2 transition
between the fitted levels, with A, S, the cancellation factor, the
total radiative decay rate of the upper level, the branching fraction,
and the standard deviations of the last two.

Usage:

    calc_TP_M1E2 [folder] [--par <file>] [--write-par]

There are no options for the physics.  Everything the Perl held in the
block of `my $...` at the top of the file -- the number of trials, the
variances, the cut-offs, the output file name -- is now read from

    calc_TP_M1E2.par

in the working directory, or from the copy that ships beside this
script if the working directory has none.  `--write-par` drops a fresh,
fully commented copy of that file into the working directory; each
parameter carries its own explanation there, and the file is the
documentation.  See `toolbox/parfile.py` for the layout.

WHAT THE RUN NEEDS

  ING11, OUTG11   the unvaried RCG run, with the M1 and/or E2 spectrum
                  in the Lund layout (`ILUND`).  ING11 must carry the
                  `//R2//` radial-integral cards if the E2 integrals
                  are to be varied.
  RCEOUT, OUTE    the converged fit: the levels and their eigenvectors,
                  and the standard deviation of every parameter.
  LEVELS1         optional; the Lande g factors and the better energies.
  IN36            the configurations, for the level labels.

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
  labelled every level it read as parity set 1.

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
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import case as case_mod                               # noqa: E402
from toolbox import (files, ing11, levels, outg11, parfile,        # noqa: E402
                     radiate, rce, spectrum, trials, vacair)

PAR_FILE = 'calc_TP_M1E2.par'
CONFIG = 'COWAN.CFG'
RCG = 'rcg11kd.exe'

#: The same program built on Unix, where gfortran's -o leaves no
#: suffix.  Both names are looked for, .exe first.
RCG_UNIX = 'rcg11kd'

#: The two multipoles this program is about.  E1 is left to `conv_out`.
MULTIPOLES = ('MAG DIP', 'ELEC QUD')

#: A wavelength between these two wave numbers, in the run's own unit,
#: is quoted in standard air -- 2000 to 20000 angstroms.
AIR_FROM, AIR_TO = 5.0, 50.0

TAB = '\t'


class Usage(Exception):
    """The command line will not do what was asked."""


class RunError(Exception):
    """The run cannot be made, or cannot be trusted."""


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

    S('r2_variance', parfile.FLOAT, 0.15, minimum=0.0, description="""
        Fractional standard deviation of the E2 radial transition
        integrals -- the //R2// cards at the end of ING11.  These are
        Hartree-Fock values with no fitted uncertainty of their own, so
        this is an estimate of how far a single-configuration HF
        integral can be from the truth.  0 leaves them alone."""),

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
        counting only the M1 and E2 channels this run holds."""),

    S('out_file', parfile.TEXT, 'tp_M1E2.txt', description="""
        Name of the tab-separated line list this program writes, in the
        working directory.  (`$out_tr_file` in the Perl, where it was a
        constant in the source and had to be edited for every ion.)"""),

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
        do either or both.  The Perl could only do the first, and
        labelled everything it read as belonging to it."""),

    S('rcg', parfile.TEXT, '', description="""
        The RCG executable to run in each trial subdirectory.  Blank
        uses the one that sits beside this program, which is the
        arrangement the .bat files assume.  Not used when there are no
        trials to run (num_trials = read_only)."""),
)

TITLE = ('calc_TP_M1E2.par -- input for calc_TP_M1E2',
         '',
         'Uncertainties on M1 and E2 transition probabilities, by',
         'repeating the RCG calculation with the fitted parameters and',
         'the E2 radial integrals displaced within their errors.')


# --------------------------------------------------------------------------
# the levels of one run
# --------------------------------------------------------------------------

class Lev:
    """One RCE level, with everything the line list asks of it."""

    __slots__ = ('one', 'kset', 'j', 'index', 'ar_tot')

    def __init__(self, one):
        self.one = one
        self.kset = one.kset
        self.j = one.j
        self.index = one.level.index       # its position in RCEOUT
        self.ar_tot = []                   # total A per trial, 0 = unvaried

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
    def g(self):
        return self.one.level.g

    def energy(self, other):
        """The energy to compare with another level's, in the same terms."""
        return (self.e_obs if self.observed and other.observed
                else self.e_calc)

    def decay(self, trial, rate):
        while len(self.ar_tot) <= trial:
            self.ar_tot.append(0.0)
        self.ar_tot[trial] += rate


class Transition:
    """One pair of levels and one multipole, over every trial."""

    __slots__ = ('upper', 'lower', 'multipole', 'delta_e', 'wavelength',
                 'ga', 's', 'cf', 'log_a', 'cancel', 'a_values', 'cf_values')

    def __init__(self, upper, lower, multipole):
        self.upper, self.lower = upper, lower
        self.multipole = multipole
        self.delta_e = 0.0            # from the fitted energies
        self.wavelength = 0.0
        self.ga = None                # the unvaried run's, None if it had none
        self.s = None
        self.cf = None
        self.log_a = trials.Accumulator()
        self.cancel = trials.Accumulator()
        self.a_values = []
        self.cf_values = []

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
    def a_mean(self):
        """The geometric mean of A over the trials -- exp(mean ln A)."""
        return math.exp(self.log_a.mean) if self.log_a.n else 0.0

    @property
    def spread(self):
        """The relative standard deviation of A, in per cent."""
        return (math.exp(self.log_a.sigma) - 1.0) * 100.0


# --------------------------------------------------------------------------
# reading one OUTG11
# --------------------------------------------------------------------------

def _kset_of(deck):
    return {card.label: card.kset for card in deck}


def _identify(csets, blocks, emax):
    """{(kset, J, eigenvector): Identified} for one RCG run."""
    out = {}
    for cset in csets:
        limit = emax[cset.index - 1] if cset.index <= 2 else 0.0
        for one in levels.identify(cset, blocks, emax=limit or None,
                                   skipped=[]):
            out[(one.kset, one.j, one.rcg)] = one
    return out


def read_transitions(folder, mapping, of_label, wanted, store, trial,
                     found, table, order):
    """Fold one OUTG11's M1 and E2 spectrum into the running statistics.

    `found` is `{(kset, J, RCE index): Lev}` and `table` the transitions
    keyed by `Transition.key`; both are built on the unvaried run and
    added to by every trial.  A transition a trial has and the unvaried
    run had not is added with no unvaried values, which is what the
    `gA is None` tests downstream mean.
    """
    seen = set()
    for line in spectrum.lines(folder=folder, multipoles=MULTIPOLES):
        multipole = radiate.OF_SPECTRUM[line.multipole]
        kset = of_label.get(line.initial.label)
        if kset is None:
            raise RunError(
                'the OUTG11 spectrum of %s names the configuration %r, which '
                'is not on any ING11 card' % (folder, line.initial.label))
        if kset not in wanted:
            continue
        first = mapping.get((kset, line.initial.j, line.initial.level))
        second = mapping.get((kset, line.final.j, line.final.level))
        if first is None or second is None:
            continue

        a = found.get((kset, first.j, first.level.index))
        b = found.get((kset, second.j, second.level.index))
        if a is None:
            a = found.setdefault((kset, first.j, first.level.index),
                                 Lev(first))
        if b is None:
            b = found.setdefault((kset, second.j, second.level.index),
                                 Lev(second))
        if a is b:
            raise RunError(
                'the %s spectrum of %s has a line from level %d (set %d, '
                'J=%s) to itself' % (multipole, folder, a.index, kset, a.j))

        # The upper level is the one higher in energy, on the observed
        # energies where both are observed and on the fitted ones
        # otherwise (`calc_TP_M1E2.pl:453`).
        upper, lower = (a, b) if a.energy(b) >= b.energy(a) else (b, a)

        difference = abs(upper.e_calc - lower.e_calc)
        if not difference:
            continue
        fitted = 1e5 / difference
        if upper.observed and lower.observed:
            spread = abs(upper.e_obs - lower.e_obs)
            # S is held fixed and gA moved to the observed wavelength;
            # see the header of this file.
            s = radiate.s_from_ga(line.ga, _printed(line), multipole)
            ga = (radiate.ga_from_s(s, 1e5 / spread, multipole)
                  if spread else line.ga)
        else:
            s = radiate.s_from_ga(line.ga, _printed(line), multipole)
            ga = line.ga

        rate = ga / (2 * upper.j + 1)
        upper.decay(trial, rate)

        key = (upper.key, lower.key, multipole)
        if key in seen:
            # `calc_TP_M1E2.pl:717` stops here too.  One pass over one
            # OUTG11 must contribute one value per transition, or the
            # statistics of that transition count a trial twice.
            raise RunError(
                'the %s spectrum of %s prints the transition between level '
                '%d (J=%s) and level %d (J=%s) of parity set %d twice'
                % (multipole, folder, upper.index, upper.j, lower.index,
                   lower.j, kset))
        seen.add(key)

        one = table.get(key)
        if one is None:
            one = Transition(upper, lower, multipole)
            one.delta_e = difference
            one.wavelength = fitted
            table[key] = one
            order.append(one)
            if not trial:
                one.ga, one.s, one.cf = ga, s, line.cf
        if rate > 0:
            # A transition with no rate at all has no logarithm; it is
            # left out of the mean rather than stopping the run, which
            # is what the Perl's log() would do.
            one.log_a.add(math.log(rate))
        one.cancel.add(line.cf)
        if trial and trial <= store:
            one.a_values.append((trial, rate))
            one.cf_values.append((trial, line.cf))


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


def vary_r_integrals(cards, options, rng, histogram):
    """New values for the //R2// radial transition integrals."""
    fraction = options['r2_variance']
    out = []
    for card in cards:
        new = rng.normal(card.value, abs(fraction * card.value))
        if fraction and card.value:
            histogram.add((new / card.value - 1.0) / fraction, keep=True)
        out.append(new)
    return out


def write_trial(deck, folder, r2, r2_values, per_set, ending):
    """Write one trial's ING11.  The deck is restored before returning."""
    was_r2 = [card.value for card in r2]
    was_par = {}
    for kset, (slots, values) in per_set.items():
        for (card, at), value in zip(slots, values):
            was_par.setdefault((id(card), at), (card, at, card.values[at]))
            card.values[at] = value
            card.edited = True
    for card, value in zip(r2, r2_values):
        ing11.set_r_value(card, value)

    files.write_deck(Path(folder) / files.ING11, ing11.write_ing11(deck),
                     ending)

    for card, at, value in was_par.values():
        card.values[at] = value
        card.edited = False
    for card, value in zip(r2, was_r2):
        ing11.set_r_value(card, value)


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
                 'spread', 'cf_mean', 'cf_sigma', 'e2_fraction',
                 'e2_sigma', 'a_values', 'cf_values', 'wavelength')

    def __init__(self, one):
        self.upper, self.lower = one.upper, one.lower
        self.kind = one.multipole
        self.ga, self.s, self.cf = one.ga, one.s, one.cf
        self.ga_mean = one.a_mean * (2 * one.upper.j + 1)
        self.spread = one.spread
        self.cf_mean = one.cancel.mean
        self.cf_sigma = one.cancel.sigma
        self.e2_fraction = 0.0 if one.multipole == radiate.M1 else 1.0
        self.e2_sigma = 0.0
        self.a_values = dict(one.a_values)
        self.cf_values = dict(one.cf_values)
        self.wavelength = 0.0

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
    both.wavelength = e2.wavelength
    both.spread = (math.exp(spread) - 1.0) * 100.0
    both.ga = None if e2.ga is None else total
    both.ga_mean = e2.ga_mean + m1.ga_mean
    both.s = None
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
        return e2
    return both


# --------------------------------------------------------------------------
# the report
# --------------------------------------------------------------------------

def header(options):
    out = ['type', 'conf1', 't1', 'J1', 'conf2', 't2', 'J2', 'lande1',
           'lande2', 'E1', 'E2', 'wl_c(A)', 'E1_exp', 'E2_exp', 'wl_exp(A)',
           'A(s-1)', 'S', 'cf', 'sum A', 'BF']
    if options['merged_list']:
        out.append('E2_frac')
        if options['num_trials']:
            out.append('d_E2frac')
    if options['num_trials']:
        out.extend(['d_BF', 'A_mean', 'S_mean', 'std_A%', 'cf_mean',
                    'std_cf'])
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


def row_text(row, options, name):
    """One printed line, or None where it falls under the cut-offs.

    Columns 1 and 2 are the LOWER and the UPPER level in that order,
    which is what the Perl prints under the names conf1 and conf2.
    """
    upper, lower = row.upper, row.lower
    at = 0                              # the unvaried run's slot
    total = upper.ar_tot[at] if at < len(upper.ar_tot) else 0.0
    if not total:
        # No M1 or E2 decay in the unvaried run -- average the trials.
        kept = [x for x in upper.ar_tot[1:] if x > 0]
        total = sum(kept) / len(kept) if kept else 0.0
    branch = (row.strength / (2 * upper.j + 1) / total) if total else 1.0
    if branch < options['min_bf']:
        return None

    difference = abs(upper.e_calc - lower.e_calc)
    if not difference:
        return None
    wavelength = 1e5 / difference
    experimental = ''
    if upper.observed and lower.observed:
        spread = abs(upper.e_obs - lower.e_obs)
        if spread:
            difference = spread
            experimental = '%.4f' % _air(1e5 / spread, spread)
    row.wavelength = 1e5 / difference     # vacuum, from the better energies
    wavelength = '%.4f' % _air(wavelength, difference)

    out = [row.kind, name(lower)[0], name(lower)[1], '%.1f' % lower.j,
           name(upper)[0], name(upper)[1], '%.1f' % upper.j,
           _g(lower), _g(upper),
           '%.4f' % lower.e_calc, '%.4f' % upper.e_calc, wavelength,
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
            out.append('%.15g' % row.e2_sigma)
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
    """d_BF, A_mean, S_mean, std_A%, cf_mean, std_cf.

    The spread quoted is the spread of the trials about the unvaried
    run, not about their own mean: the offset between the two is a real
    part of the uncertainty and is added in quadrature, in the
    logarithm for A and directly for the cancellation factor.
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
    if options['print_dbf'] and options['store_trials']:
        fractions = []
        for at in range(1, options['store_trials'] + 1):
            a = row.a_values.get(at)
            if a is None:
                continue
            total = (upper.ar_tot[at] if at < len(upper.ar_tot) else 0.0)
            if total:
                fractions.append(a / total)
        found = trials.stdev(branch, fractions)
        if found is not None:
            d_bf = '%.15g' % found

    a_mean = row.ga_mean / (2 * upper.j + 1)
    # S_mean goes with the VACUUM wavelength of the better energies,
    # which is the wavelength gA itself was moved to; the wl_c column
    # beside it may be in air (`calc_TP_M1E2.pl:826-833`).
    s_mean = ('' if row.kind == 'M1+E2' or not row.wavelength
              else '%.4E' % radiate.s_from_ga(row.ga_mean, row.wavelength,
                                              row.kind))
    return [d_bf, '%.5E' % a_mean, s_mean, '%.2f' % spread,
            '%.7f' % row.cf_mean, '%.7f' % cf_sigma]


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

def run(folder, options, out):
    folder = Path(folder)
    deck = ing11.read_ing11(folder=folder)
    ending = files.line_ending(files.ING11, folder)
    emax = deck.emax
    of_label = _kset_of(deck)

    wanted = {'1': (1,), '2': (2,), 'both': (1, 2)}[
        options['parity_sets'].lower()]
    wanted = tuple(k for k in wanted if deck.main.nconfigurations(k))
    if not wanted:
        raise RunError('%s has no configuration in parity set %s'
                       % (folder / files.ING11, options['parity_sets']))

    ntrials = options['num_trials']
    already = min(options['read_only'], ntrials)
    store = min(options['store_trials'], ntrials)
    options['store_trials'] = store

    print('Reading the unvaried run...')
    base = case_mod.Case.read(folder, couplings=('LS',), use_emax=True)
    if base.unmatched:
        print('   %d fitted levels above E0MAX have no RCG eigenvector'
              % len(base.unmatched))
    mapping = {(one.kset, one.j, one.rcg): one for one in base.identified}

    r2 = ing11.read_r_integrals(deck, order=2)
    print('   %d E2 radial integrals (//R2// cards)' % len(r2))
    if ntrials and options['r2_variance'] and not r2:
        print('   WARNING: r2_variance is set but the deck has no //R2// '
              'card; only the parameters will be varied.')

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
    print('Reading transitions...')
    read_transitions(folder, mapping, of_label, wanted, store, 0,
                     found, table, order)
    print('   %d transitions' % len(order))

    rng = trials.Randomizer(options['random_seed'] or None)
    if ntrials:
        print('Random seed: %r' % (rng.seed,))
    histograms = {name: trials.Histogram(name) for name in
                  ('parameters', 'lnA', 'A', 'cf')}

    executable = find_rcg(options['rcg']) if ntrials > already else None
    for trial in range(1, ntrials + 1):
        where = folder / ('trial%d' % trial)
        where.mkdir(exist_ok=True)
        if trial > already:
            print('Trial %d of %d...' % (trial, ntrials))
            per_set = {kset: (positions,
                              vary_parameters(fitted, options, rng,
                                              histograms['parameters']))
                       for kset, (positions, fitted) in slots.items()}
            write_trial(deck, where, r2,
                        vary_r_integrals(r2, options, rng,
                                         histograms['parameters']),
                        per_set, ending)
            for name in ('LEVELS1', 'LEVELS2'):
                source = files.resolve(name, folder, required=False)
                if source is not None:
                    shutil.copyfile(source, where / name)
            run_rcg(where, executable)
        else:
            print('Trial %d of %d (already computed)...' % (trial, ntrials))
        if slots:
            blocks = levels.by_block(outg11.read_outg11(folder=where), ['LS'])
            trial_map = _identify(base.csets, blocks, emax)
        else:
            trial_map = mapping
        read_transitions(where, trial_map, of_label, wanted, store, trial,
                         found, table, order)

    write(out, base, order, options, ntrials)
    if options['print_distribution']:
        distribution(order, options, histograms)


def write(out, base, order, options, ntrials):
    """The line list itself.

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

    groups = [''] if options['merged_list'] else [radiate.M1, radiate.E2]
    with open(out, 'w', newline='') as fh:
        fh.write(TAB.join(header(options)) + '\n')
        for at, kind in enumerate(groups):
            rows = [Row(one) for one in order
                    if not kind or one.multipole == kind]
            if options['merged_list']:
                rows = merge(rows, ntrials, options['min_mixed_fraction'])
            printed = []
            for row in rows:
                text = row_text(row, options, name)
                if text is not None:
                    printed.append((row.wavelength, row.upper.key,
                                    row.lower.key, text))
            printed.sort(key=lambda x: x[:3])
            for _, _, _, text in printed:
                fh.write(TAB.join(text) + '\n')
    print('%s: %d level pairs written.' % (out, len(order)))


def distribution(order, options, histograms):
    """The histograms of how far the trials strayed, on the terminal."""
    for one in order:
        reference = one.a
        if reference is None or reference <= 0 or not one.log_a.sigma:
            continue
        for _, a in one.a_values:
            if a > 0:
                histograms['lnA'].add(math.log(a / reference)
                                      / one.log_a.sigma, keep=True)
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

def main(argv):
    folder, named, fresh = '.', None, False
    rest = list(argv)
    while rest:
        arg = rest.pop(0)
        if arg in ('--write-par', '--write_par'):
            fresh = True
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

    if fresh:
        path = parfile.write(Path(folder) / PAR_FILE, SPECS, title=TITLE)
        print('%s written.  Read it through, then run calc_TP_M1E2 again.'
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
              'needs, then run calc_TP_M1E2 again.')
        return 0
    print('Parameters from %s' % path)
    options = parfile.read(path, SPECS)
    run(folder, options, Path(folder) / options['out_file'])
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('calc_TP_M1E2: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1].split(
                             'WHAT THE RUN NEEDS')[0]))
        sys.exit(2)
    except (RunError, parfile.ParFileError, case_mod.CaseError,
            ing11.Ing11Error, rce.RceError, levels.LevelsError,
            spectrum.SpectrumError, outg11.Outg11Error,
            OSError) as problem:
        sys.stderr.write('calc_TP_M1E2: %s\n' % problem)
        sys.exit(1)
