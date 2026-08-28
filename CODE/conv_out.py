#!/usr/bin/env python3
"""Write the level and line lists of a finished RCG/RCE run.

The Python replacement for `conv_out.pl`.  It produces the same two
tab-separated files -- a level list with one block per coupling and a
line list sorted by wavelength -- and, on request, the two input files
of the IDEN2 line-identification code.

What it does that the Perl could not is what `printout.py` does: labels
built from `genealogy.txt`, n taken from IN36 through the ING11 slots,
one electron count in every label, and any of RCG's couplings rather
than LS and Cowan's jj alone.  The transitions themselves are read from
the Lund layout of the OUTG11 spectrum, which prints each level's number
within its J block -- so a line is tied to an eigenvector by that number
instead of being matched back on energy and term label.

Usage, the historical order, which the .bat file passes straight through:

    conv_out <levels file> <lines file> [u<E>] RCE|noRCE
             LTE|LTEABS|BF:<ns>|DR|AUGER|CM <T in kK> <max intensity>
             [<AI file, first parity> <AI file, second parity>]
             [IDEN [truncate] [linear] [convert2log]]

and the same things as long options, which may be mixed in freely:

    --couplings LS,JJ      --drop-full-shells 5p,4d
    --components 5         --min-percent 1.5
    --upper-bound 500      --model LTE   --temperature 100
    --max-intensity 1000   --ai even.txt,odd.txt
    --no-rce --iden --truncate --linear --convert2log

Every one of those settings may instead be written into a `conv_out.par`
in the working directory, in the author's own layout -- value first,
then the key, then prose saying what it is for.  The par file supplies
defaults and the command line overrides them, so nothing here changes
for a run that does not have one.  `conv_out` with no par file and no
arguments writes a commented sample and stops; `--write-par` writes one
on request, and `--par <file>` reads one by name.

Energies are in the unit of the run -- 1000 cm-1 in every deck of this
package -- so a wavelength in angstroms is 1e5 / delta E, and the
effective temperature is given in the same unit.
"""

import math
import os
import sys
from pathlib import Path

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import case as case_mod                               # noqa: E402
from toolbox import (ai, asd, designate, fundconst, iden, in36,   # noqa: E402
                     levels, parfile, radiate, rce, spectrum, vacair)

#: The parameter file.  Unlike `calc_TP.par` this one is NOT
#: shipped in CODE/: it must be absent on a first run, or the sample
#: below would never be written.
PAR_FILE = 'conv_out.par'

#: The intensity models, and what each one means.
LTE = 'LTE'
LTEABS = 'LTEABS'
BF = 'BF'
DR = 'DR'
AUGER = 'AUGER'
CM = 'CM'
MODELS = (LTE, LTEABS, BF, DR, AUGER, CM)

#: A wavelength between these two wavenumbers is quoted in standard air,
#: as every published line list does between 2000 and 20000 angstroms.
AIR_FROM, AIR_TO = 5.0, 50.0

#: Cowan's own beam-foil default: the time an ion takes to cross the
#: spectrograph's entrance slit, in nanoseconds (`conv_out.pl:50`).
FLIGHT = 0.114


# --------------------------------------------------------------------------
# the parameter file
# --------------------------------------------------------------------------
#
# The command line is long, strictly ordered and tells the reader
# nothing: two file names, then keywords, then two bare numbers that
# are told apart only by the order they stand in.  The par file says
# what each setting is for, and is read before the command line so that
# the command line still has the last word.  Nothing here is required
# to use it: the historical positional form works exactly as it did.

S = parfile.Spec
SPECS = (
    S('levels', parfile.TEXT, None, description="""
        Where the level list goes -- the first of the two output files,
        one block per coupling scheme, tab separated.  Required: there
        is no default and no run without it."""),

    S('lines', parfile.TEXT, None, description="""
        Where the line list goes -- the transitions, sorted by
        wavelength, tab separated.  Required.  With model CM this file
        holds the cascade matrix instead, comma separated because it is
        read in a spreadsheet."""),

    S('model', parfile.TEXT, None, choices=MODELS, description="""
        What the printed intensities mean.  Required; one of:

        LTE -- relative emission intensities, from Boltzmann
        populations of the upper level of each transition at the
        effective temperature below.

        LTEABS -- relative reduced absorption intensities, from
        Boltzmann populations of the lower level.

        BF -- beam-foil: the population still radiating after the ions
        have crossed the viewing aperture of the spectrometer.  See
        flight.

        DR -- dielectronic recombination.  Built from autoionization
        rates, so it requires ai.

        AUGER -- the Auger electron spectrum: one row per
        autoionization channel rather than per radiative transition.
        Requires ai.

        CM -- cascade matrix.  Not a line list at all but a square
        table of A(i->j) over the levels, comma separated.  A matrix
        has one axis, so only the first entry of couplings is used,
        and the level list gains the five columns n, l, occ, S and L,
        which classify a level for this model."""),

    S('temperature', parfile.FLOAT, None, minimum=0.0, description="""
        Effective temperature, in the unit of the run -- 10^3 cm-1 in
        every deck of this package.  Required in every model, and that
        is the trap worth saying out loud: only LTE and LTEABS use it,
        but the run stops without it whatever the model was."""),

    S('max_intensity', parfile.FLOAT, None, minimum=1.0, description="""
        The intensity given to the strongest line; every other line is
        scaled to it.  Required, and at least 1."""),

    S('flight', parfile.FLOAT, FLIGHT, minimum=0.0, description="""
        Time of flight in nanoseconds -- how long an ion in the beam
        takes to cross the viewing aperture of the spectrometer.  Used
        by model BF and by nothing else; 0.114 is Cowan's own value.
        On the command line it is written into the model keyword
        instead, as BF:0.114."""),

    S('no_rce', parfile.BOOL, False, description="""
        1 to ignore the RCE output and work from the calculated
        energies alone.  0, the usual case, fetches the experimental
        energies from RCEOUT and computes Ritz wavelengths with them.
        On the command line this is the bare keyword RCE or noRCE."""),

    S('couplings', parfile.TEXT, 'LS,JJ', description="""
        Which coupling schemes the level list holds, comma separated
        and in the order the blocks are to be written: LS, JJ, TRUEjj,
        JJJK, LSLK, LSJK, LSJLKS, LSJLSJ.  Each becomes its own block.
        The line list and the cascade matrix have one label column
        apiece, so they use the first of these and no other."""),

    S('drop_full_shells', parfile.TEXT, '', description="""
        Subshells left out of every printed configuration label, comma
        separated, as 5p,4d.  A shell that is full in every
        configuration of the run distinguishes nothing and only makes
        the labels longer.  Blank prints them all."""),

    S('components', parfile.INT, 5, minimum=1, description="""
        How many eigenvector components are printed beside each level,
        in descending order of percentage."""),

    S('min_percent', parfile.FLOAT, 1.5, minimum=0.0, description="""
        Components weaker than this percentage are not printed, however
        many components asked for."""),

    S('upper_bound', parfile.FLOAT, None, description="""
        Levels above this energy, in the unit of the run, are left out
        of both files.  Blank keeps them all.  On the command line this
        is the bare keyword u500."""),

    S('ai', parfile.TEXT, '', description="""
        The two files of autoionization rates, comma separated, in the
        order the configuration sets stand in IN36 -- what the Perl
        called the even file and the odd file.  Required by DR and
        AUGER and unused by every other model.

        There are two of them because of how they are made: read_aa is
        run in two separate directories, one holding the even
        configurations and one the odd, and each writes the rates of
        the set it was run on.  Name them here with paths that reach
        them from this directory."""),

    S('iden', parfile.BOOL, False, description="""
        1 to write the two input files of the IDEN2 line-identification
        code, ENLEV.DAT and TRANS.DAT, beside the level and line lists.
        The files dlv.dat and numset.dat must be present."""),

    S('truncate', parfile.BOOL, False, description="""
        1 to drop from the IDEN2 level list those levels that no
        calculated transition connects to any other.  Only iden reads
        it."""),

    S('linear', parfile.BOOL, False, description="""
        1 to write the IDEN2 intensities in TRANS.DAT on a linear
        scale; 0, the default, writes them on a log scale.  Only iden
        reads it."""),

    S('convert2log', parfile.BOOL, False, description="""
        1 to write a new DLV.DAT with the observed intensities
        converted to a log scale, I_log = 10*ln(I).  This asks the
        opposite of linear and the two may not both be set.  Only iden
        reads it."""),
)

TITLE = ('conv_out.par -- input for conv_out',
         '',
         'The level list and the line list of a finished RCG/RCE run.',
         'Every setting here can also be given on the command line, and',
         'the command line has the last word.')


class Usage(Exception):
    """The command line will not do what was asked."""


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

class Options:

    #: What a setting is when neither the command line nor the par
    #: file has said anything about it.  Held here rather than assigned
    #: in `__init__` because `None` has to mean "nobody has said" until
    #: both sources have been consulted; see `_settle`.
    DEFAULT = {'couplings': ('LS', 'JJ'), 'drop': (), 'components': 5,
               'min_percent': 1.5, 'upper_bound': float('inf'),
               'flight': FLIGHT, 'ai': (), 'no_rce': False, 'iden': False,
               'truncate': False, 'linear': False, 'convert2log': False}

    #: Every setting, so that `__init__` can clear them all in one go.
    NAMES = ('levels', 'lines', 'couplings', 'drop', 'components',
             'min_percent', 'upper_bound', 'model', 'temperature',
             'max_intensity', 'flight', 'ai', 'no_rce', 'iden', 'truncate',
             'linear', 'convert2log')

    def __init__(self, argv, par=None):
        for name in self.NAMES:
            setattr(self, name, None)

        positional = []
        rest = list(argv)
        while rest:
            arg = rest.pop(0)
            if not arg.startswith('--'):
                if not self._keyword(arg):
                    positional.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.replace('-', '_')
            if name in ('no_rce', 'iden', 'truncate', 'linear', 'convert2log'):
                setattr(self, name, True)
                continue
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name.replace('_', '-'))
                value = rest.pop(0)
            self._option(name, value)

        self._positional(positional)
        self._settle(par or {})
        self._check()

    # -- the par file, and then the built-in defaults --------------------

    def _settle(self, par):
        """Fill from the par file what the command line did not set.

        This order is what keeps the par file purely additive: anything
        typed on the command line already holds a value by the time we
        get here and is never overwritten, and what neither source
        mentioned falls through to `DEFAULT`.  A blank field in the par
        file reads back as the spec's default, which for a setting with
        no default is `None` -- so blank there means "nobody has said"
        as well.
        """
        for name, value in _par_values(par).items():
            if getattr(self, name) is None:
                setattr(self, name, value)
        for name, value in self.DEFAULT.items():
            if getattr(self, name) is None:
                setattr(self, name, list(value) if isinstance(value, tuple)
                        else value)

    # -- the historical bare keywords, which may stand anywhere ----------

    def _keyword(self, arg):
        upper = arg.upper()
        if upper in ('IDEN', 'TRUNCATE', 'LINEAR', 'CONVERT2LOG'):
            setattr(self, upper.lower(), True)
            return True
        if upper in ('RCE', 'NORCE'):
            self.no_rce = upper == 'NORCE'
            return True
        if upper in MODELS or upper.startswith('BF:'):
            self._model(upper)
            return True
        if arg[:1] in 'uU' and _number(arg[1:]) is not None:
            self.upper_bound = float(arg[1:])
            return True
        return False

    def _option(self, name, value):
        if name == 'couplings':
            self.couplings = [_coupling(c) for c in value.split(',')
                              if c.strip()]
        elif name == 'drop_full_shells':
            self.drop = [_nl(s) for s in value.split(',') if s.strip()]
        elif name == 'components':
            self.components = int(value)
        elif name == 'min_percent':
            self.min_percent = float(value)
        elif name == 'upper_bound':
            self.upper_bound = float(value)
        elif name == 'model':
            self._model(value.upper())
        elif name == 'temperature':
            self.temperature = float(value)
        elif name == 'max_intensity':
            self.max_intensity = float(value)
        elif name == 'ai':
            self.ai = [s.strip() for s in value.split(',') if s.strip()]
        elif name in ('levels', 'lines'):
            setattr(self, name, value)
        else:
            raise Usage('--%s is not an option' % name.replace('_', '-'))

    def _model(self, text):
        if text.startswith('BF:'):
            flight = _number(text[3:])
            if flight is None or flight <= 0:
                raise Usage('BF wants a time of flight in nanoseconds, as '
                            'BF:0.114, not %r' % text)
            self.model, self.flight = BF, flight
            return
        if text not in MODELS:
            raise Usage('%r is not one of %s' % (text, ', '.join(MODELS)))
        if text == BF:
            raise Usage('BF wants a time of flight in nanoseconds, as BF:0.114')
        self.model = text

    def _positional(self, positional):
        """What is left over is the two file names, then the two numbers.

        The historical command line puts the effective temperature and
        the maximum scaled intensity after the model keyword, so the
        numbers are taken in that order from whatever the keywords did
        not claim; any AI file names follow them.
        """
        names = [one for one in positional if _number(one) is None]
        numbers = [float(one) for one in positional if _number(one) is not None]
        if self.levels is None and names:
            self.levels = names.pop(0)
        if self.lines is None and names:
            self.lines = names.pop(0)
        if self.temperature is None and numbers:
            self.temperature = numbers.pop(0)
        if self.max_intensity is None and numbers:
            self.max_intensity = numbers.pop(0)
        if self.ai is None and names:
            # `and names`, so that an empty command line leaves `ai`
            # unset rather than set to nothing -- what the par file has
            # to say about it still has to get its turn in `_settle`.
            self.ai = names[:2]
            names = names[2:]
        if names or numbers:
            raise Usage('I do not know what to do with %s'
                        % ', '.join(names + ['%g' % n for n in numbers]))

    def _check(self):
        if not self.levels or not self.lines:
            raise Usage('two output file names are needed, one for the '
                        'levels and one for the lines')
        if self.model is None:
            raise Usage('an intensity model is needed: %s' % ', '.join(MODELS))
        if self.temperature is None or self.temperature <= 0:
            raise Usage('an effective temperature in kK is needed; it is used '
                        'only by LTE and LTEABS, but it is always required')
        if self.max_intensity is None or self.max_intensity < 1:
            raise Usage('a maximum scaled intensity of 1 or more is needed')
        if len(self.ai) == 1:
            raise Usage('autoionization rates come as two files, one per '
                        'parity set; only %s was given' % self.ai[0])
        if self.model in (DR, AUGER) and not self.ai:
            raise Usage('the %s model is built from autoionization rates, so '
                        'it needs the two read_aa files' % self.model)
        if self.linear and self.convert2log:
            raise Usage('linear and convert2log ask for opposite things')
        if self.model == CM:
            self.couplings = [self.couplings[0]]

    @property
    def delimiter(self):
        """The cascade matrix is a spreadsheet; everything else is TSV."""
        return ',' if self.model == CM else '\t'


def _par_values(par):
    """A par file's settings, named and shaped as `Options` holds them.

    Only what the file actually says comes back: a key it left blank is
    `None` in `par` where the spec has no default, and is dropped here
    so that `Options.DEFAULT` gets its turn.  The comma-separated
    fields are split exactly as `Options._option` splits the same text
    on the command line, so a value means the same thing in both
    places.  `model` is not passed through `_model`: the par file keeps
    the time of flight in its own `flight` field rather than glued to
    the keyword, so a bare BF is right here and wrong there.
    """
    out = {}
    for key in ('levels', 'lines', 'model'):
        if par.get(key):
            out[key] = par[key].strip()
    if 'model' in out:
        out['model'] = out['model'].upper()
    for key in ('temperature', 'max_intensity', 'upper_bound', 'components',
                'min_percent', 'flight'):
        if par.get(key) is not None:
            out[key] = par[key]
    if par.get('couplings'):
        out['couplings'] = [_coupling(c) for c in par['couplings'].split(',')
                            if c.strip()]
    if par.get('drop_full_shells'):
        out['drop'] = [_nl(one) for one in par['drop_full_shells'].split(',')
                       if one.strip()]
    if par.get('ai'):
        out['ai'] = [one.strip() for one in par['ai'].split(',') if one.strip()]
    for key in ('no_rce', 'iden', 'truncate', 'linear', 'convert2log'):
        if par.get(key) is not None:
            out[key] = bool(par[key])
    return out


def _number(text):
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def _coupling(name):
    name = name.strip()
    for known in designate.TABLE_OF:
        if known.lower() == name.lower():
            return known
    raise Usage('%r is not one of RCG\'s couplings (%s)'
                % (name, ', '.join(sorted(designate.TABLE_OF))))


def _nl(text):
    text = text.strip()
    for i, ch in enumerate(text):
        if ch.isalpha():
            try:
                return (int(text[:i]), in36.LBCD.index(ch.lower()))
            except ValueError:
                break
    raise Usage('%r is not a subshell like 5p' % text)


# --------------------------------------------------------------------------
# the levels of a run, with everything a line list needs hung on them
# --------------------------------------------------------------------------

class Lev:
    """One level: where it came from, what it is called, what it does."""

    __slots__ = ('one', 'kset', 'j', 'rcg', 'key', 'energy', 'number',
                 'ar_tot', 'aa_tot', 'e_ai', 'dr', 'channels')

    def __init__(self, one, key):
        self.one = one                  # levels.Identified
        self.kset, self.j, self.rcg = one.kset, one.j, one.rcg
        self.key = key                  # levels.eigenvalue, tie-break and all
        self.energy = key[0]            # the RCG eigenvalue, for the physics
        self.number = 0                 # Lev_#, filled once all are known
        self.ar_tot = 0.0               # sum of A over the lines from here
        self.aa_tot = 0.0
        self.e_ai = 0.0
        self.dr = 0.0
        self.channels = {}              # channel name -> Auger rate

    @property
    def level(self):
        return self.one.level

    @property
    def e_calc(self):
        return self.one.level.e_calc

    @property
    def observed(self):
        return self.one.level.observed

    @property
    def group(self):
        """The relative-energy group this level is in, or `''`.

        What `rce.comparable` asks about: a difference of two observed
        energies is on one scale only inside one group.
        """
        return self.one.level.group

    @property
    def absolute(self):
        """Whether the observed energy is on the deck's absolute scale."""
        return self.one.level.absolute

    @property
    def e_obs(self):
        """The observed energy where there is one, else the fitted one."""
        level = self.one.level
        return level.e_obs if level.observed else level.e_calc

    @property
    def width(self):
        """The total decay rate, radiative and Auger together, in s-1."""
        return self.aa_tot + self.ar_tot

    @property
    def lifetime(self):
        """In nanoseconds.  A level with no way out lives for ever."""
        total = self.width
        return 1e9 / total if total > 1e-19 else 1e27

    @property
    def radiative_branch(self):
        total = self.width
        return self.ar_tot / total if total > 0 else 0.0


class Trans:
    """One transition, its two levels already resolved."""

    __slots__ = ('upper', 'lower', 'ga', 'cf', 'multipole', 'wavelength',
                 'rcg_wavelength', 'intensity')

    def __init__(self, upper, lower, ga, cf, multipole):
        self.upper, self.lower = upper, lower
        self.ga, self.cf, self.multipole = ga, cf, multipole
        self.wavelength = None          # printed: from the fitted energies
        self.rcg_wavelength = None      # modelled: RCG's own, see _intensities
        self.intensity = 0.0

    @property
    def a(self):
        """The transition probability of the upper level, A = gA / (2J+1)."""
        return self.ga / (2 * self.upper.j + 1)

    @property
    def branch(self):
        total = self.upper.width
        return self.a / total if total > 0 else 0.0


def _collect(case, options):
    """Every level of the run, keyed by (parity set, J, eigenvector).

    In Lev_# order.  `levels.numbered` is the one place that ordering is
    written down -- `calc_TP` numbers its levels by the same rule, so
    that its line list and this one can be joined -- and what is left
    for `_number_levels` is to drop the levels above `upper_bound` and
    count off the rest.
    """
    return {(one.kset, one.j, one.rcg): Lev(one, key)
            for one, key in levels.numbered(case, cm=options.model == CM)}


def _number_levels(found, options):
    """Give every level its Lev_#, and return them in that order.

    Within a parity set the order is the RCG eigenvalue's, as everywhere
    else in this package; the sets follow one another, first then
    second, which is what makes a level number tell the parity at a
    glance.  The cascade matrix is the exception: it is a single square
    matrix over all levels, so its rows run in energy order across both
    sets.

    Both branches carry `levels.eigenvalue`'s own tie-break, so that two
    levels of one set with the same printed eigenvalue -- common, RCG
    printing it to two decimals -- fall here in the order `printout`
    gives them, which is `levels.order_key`.  `found` arrives in that
    order from `_collect`; all that is left is the bound.
    """
    order = [one for one in found if one.e_calc <= options.upper_bound]
    for number, one in enumerate(order, 1):
        one.number = number
    return order


# --------------------------------------------------------------------------
# the transitions
# --------------------------------------------------------------------------

def _kset_of(case):
    """The parity set each ING11 configuration label belongs to."""
    out = {}
    for card in case.deck:
        out[card.label] = card.kset
    return out


def _read_transitions(case, found, options):
    """Read the OUTG11 spectrum and tie both ends of every line to a level.

    RCG prints each pair of levels once per multipole.  Where a run
    holds both M1 and E2 the same pair comes up twice, and the author's
    rule is to keep the stronger of the two and drop the other -- one
    line in a list may carry only one probability.  `conv_out.pl` added
    both to the upper level's total decay rate even when it dropped one
    of them; here the total takes only what survives.
    """
    of_label = _kset_of(case)
    ksets = sorted({card.kset for card in case.deck})
    by_pair = {}
    out = []
    missing = set()

    for line in spectrum.lines(folder='.'):
        kset = of_label.get(line.initial.label)
        if kset is None:
            raise case_mod.CaseError(
                'the OUTG11 spectrum names the configuration %r, which is '
                'not on any ING11 card' % line.initial.label)
        other = kset if line.same_parity else _other(ksets, kset)
        first = found.get((kset, line.initial.j, line.initial.level))
        second = found.get((other, line.final.j, line.final.level))
        if first is None or second is None:
            missing.add((kset, line.initial.j, line.initial.level)
                        if first is None
                        else (other, line.final.j, line.final.level))
            continue
        if (line.initial.energy > options.upper_bound
                or line.final.energy > options.upper_bound):
            continue

        upper, lower = ((first, second) if line.initial.energy
                        >= line.final.energy else (second, first))
        key = (upper.kset, upper.j, upper.rcg,
               lower.kset, lower.j, lower.rcg)
        already = by_pair.get(key)
        if already is None:
            one = Trans(upper, lower, line.ga, line.cf, line.multipole)
            by_pair[key] = one
            out.append(one)
            upper.ar_tot += one.a
        elif line.ga > already.ga:
            upper.ar_tot += (line.ga - already.ga) / (2 * upper.j + 1)
            already.ga, already.cf = line.ga, line.cf
            already.multipole = line.multipole
    if missing:
        print('   %d levels named in the spectrum are not in the level list'
              % len(missing))
    return out


def _other(ksets, kset):
    for one in ksets:
        if one != kset:
            return one
    return kset


def _wavelength(upper, lower):
    """The vacuum wavelength in angstroms, from the fitted energies.

    This is the wavelength the line list prints and sorts on: where RCE
    has fitted the levels its energies are the better ones.  It is *not*
    the wavelength the intensity is modelled at -- see `_intensities`.
    Two levels that come out equal have no line, and are dropped by the
    caller.
    """
    difference = abs(upper.e_calc - lower.e_calc)
    return 1e5 / difference if difference else None


def _rcg_wavelength(upper, lower):
    """The vacuum wavelength in angstroms, from RCG's own eigenvalues.

    The intensity model uses this one and not the fitted wavelength, on
    the author's instruction (2026-08-15).  RCG's A-value belongs to
    RCG's wavelength, and moving it onto a fitted or observed one is a
    rescaling: in a real case the observed levels are a small fraction
    of all, so it would be applied to some lines and not others, with no
    good answer for a line that has one known level or none.  Cowan's
    codes give unscaled A-values, that is what people expect, and the
    error this leaves in a relative intensity is far smaller than the
    error in Cowan's line strengths.  Scaling, where it is wanted at
    all, is a case-by-case job at publication time.
    """
    difference = abs(upper.energy - lower.energy)
    return 1e5 / difference if difference else None


def _intensities(transitions, options):
    """The modelled intensity of every line, scaled to the maximum asked for.

    Emission follows gA / lambda, weighted by the population of the
    upper level; absorption follows gf * lambda**2 weighted by the
    population of the lower one.  Everything else is a multiplier on
    that, and they are the author's models:

      LTE, LTEABS  a Boltzmann population at the effective temperature
      BF           a beam that has been in view for a flight time, so a
                   level contributes for as long as it lives
      DR           dielectronic recombination, so the population follows
                   the capture rate into the level
      AUGER        the electron spectrum, handled on the levels instead

    Everything here is RCG's: the wavelength is `_rcg_wavelength` and
    the populations follow RCG's eigenvalues, so gA is used exactly as
    RCG computed it and no line is rescaled to a fitted or observed
    wavelength.
    """
    t = options.temperature
    for one in transitions:
        if one.rcg_wavelength is None:
            continue
        if options.model == LTEABS:
            base = radiate.reduced_absorption(
                one.ga, one.rcg_wavelength, radiate.OF_SPECTRUM[one.multipole])
            base *= math.exp(-one.lower.energy / t)
        else:
            base = one.ga / one.rcg_wavelength
            if options.model == LTE:
                base *= math.exp(-one.upper.energy / t)
        if options.model == BF:
            life = one.upper.lifetime
            base *= life * (1 - math.exp(-options.flight / life))
        elif options.model == DR:
            base *= one.upper.lifetime * one.upper.dr
        else:
            base *= one.upper.radiative_branch
        one.intensity = base

    biggest = max([one.intensity for one in transitions], default=0.0)
    if biggest > 0:
        for one in transitions:
            one.intensity *= options.max_intensity / biggest


# --------------------------------------------------------------------------
# the level list
# --------------------------------------------------------------------------

def _level_header(options):
    out = ['Cpl', 'par', 'Lev_#', 'Ee', 'e_c', 'Ec', 'J']
    if not options.no_rce:
        out.append('Lande_g')
    if options.ai:
        out.extend(['Eai', 'AA_tot'])
    out.append('AR_tot')
    # n, l, occ, S and L classify a level for the cascade-matrix model and
    # are of no use to the others, so they are written only with CM.  The
    # Perl wrote them in every model.
    if options.model == CM:
        out.extend(['n', 'l', 'occ', 'S', 'L'])
    for at in range(1, options.components + 1):
        out.extend(['%%%d' % at, 'conf%d' % at, 'term%d' % at])
    return out


def _outer_shell(found):
    """The outermost open subshell of a designation: (n, l, occupation).

    `conv_cowan.pl` took the last subshell of the printed label, which
    was the outermost open one only because the Perl dropped the closed
    shells.  They are kept now, so the open one is asked for by name --
    it is what the column is for, being how a level is grouped in a
    published table.
    """
    open_shells = [tok for tok in found.level.tokens
                   if isinstance(tok, asd.Subshell) and tok.occ
                   and tok.occ != 2 * (2 * tok.lval + 1)]
    if not open_shells:
        open_shells = [tok for tok in found.level.tokens
                       if isinstance(tok, asd.Subshell)]
    if not open_shells:
        return '', '', ''
    last = open_shells[-1]
    return last.n, last.lval, _plain(last.occ)


def _term_quantum_numbers(found):
    """S and L of a level's leading LS term, as numbers.

    Blank in any other coupling: a J1J2 or K term has neither.
    """
    term = found.level.term
    if not isinstance(term, asd.LSTerm):
        return '', ''
    return '%4.1f' % ((term.mult - 1) / 2.0), '%d' % asd._l_value(term.l)


def _plain(value):
    return '%d' % value if float(value).is_integer() else '%g' % value


def _level_rows(case, order, coupling, options):
    out = []
    for one in order:
        parity = 'e' if case.card(one.kset, 1).parity == 0 else 'o'
        level = one.level
        row = [coupling, parity, '%d' % one.number,
               '%12.6f' % one.e_obs, '' if level.observed else '*',
               '%12.4f' % one.e_calc, '%4.1f' % one.j]
        if not options.no_rce:
            row.append('' if level.g is None else '%.3f' % level.g)
        if options.ai:
            row.extend(['%12.4f' % one.e_ai, '%9.2e' % one.aa_tot])
        row.append('%9.2e' % one.ar_tot)
        row.extend(_components(case, coupling, one, options))
        out.append(row)
    return out


def _components(case, coupling, one, options):
    """The eigenvector components, after the n, l, occ, S, L columns.

    Those five come first and only with model CM -- see `_level_header`.
    """
    weights = case.amplitudes(coupling, one.kset, one.j, one.rcg)
    out = []
    for rank, (amplitude, state) in enumerate(weights[:options.components], 1):
        percent = int(amplitude * amplitude * 100 + 0.5)
        if percent < options.min_percent and rank > 1:
            break
        found = case.designation(coupling, one.kset, one.j, state)
        if rank == 1 and options.model == CM:
            n, l, occ = _outer_shell(found)
            s, big_l = _term_quantum_numbers(found)
            out.extend(['%s' % n, '%s' % l, '%s' % occ, s, big_l])
        out.extend(['%d' % percent, found.configuration, found.term])
    if not out and options.model == CM:
        out = ['', '', '', '', '']
    return out


def _leading(case, coupling, one):
    """The designation of a level's largest component: its name."""
    weights = case.amplitudes(coupling, one.kset, one.j, one.rcg)
    found = case.designation(coupling, one.kset, one.j, weights[0][1])
    return found.configuration, found.term


# --------------------------------------------------------------------------
# the line list
# --------------------------------------------------------------------------

def _line_header(options):
    out = ['level1', '', 'J1', 'level2', '', 'J2', 'E1', 'E2', 'lid1', 'lid2',
           'wl(A)']
    if not options.no_rce:
        out.extend(['E1_exp', 'E2_exp', 'wl_exp(A)'])
    out.extend(['Icalc', 'g2A21', 'cf', 'lft1(ns)', 'lft2(ns)'])
    if options.ai:
        out.append('tot AI')
    out.extend(['sum A2i', 'A21', 'BRrad'])
    if options.model in (LTE, LTEABS):
        out.extend(['lande1', 'lande2'])
    return out


def _air(wavelength, difference):
    """The wavelength as a line list quotes it: air between 2000 and 20000 A."""
    if AIR_FROM < difference < AIR_TO:
        return vacair.to_air(wavelength)
    return wavelength


def _line_rows(case, transitions, coupling, options):
    """One row per line, in ascending wavelength.

    The lower level comes first, as `conv_out.pl` prints it: the columns
    called level1 and E1 are the lower ones, and lid1 is its number.
    """
    out = []
    for one in sorted(transitions, key=lambda t: (t.wavelength, t.upper.number,
                                                  t.lower.number)):
        if one.intensity < 0.5:
            continue
        upper, lower = one.upper, one.lower
        difference = abs(upper.e_calc - lower.e_calc)
        experimental = ''
        # The observed wavelength comes from a *difference* of two
        # observed energies, so it survives an unknown common offset --
        # but only where the two share it.  Across a group boundary the
        # offset does not cancel and the fitted energies are the better
        # answer, which is the branch this already had.
        if rce.comparable(upper, lower):
            spread = abs(upper.e_obs - lower.e_obs)
            difference = spread or difference
            if spread:
                experimental = '%.4f' % _air(1e5 / spread, spread)
        wavelength = _air(one.wavelength, difference)

        low = _leading(case, coupling, lower)
        high = _leading(case, coupling, upper)
        row = [low[0], low[1], '%.1f' % lower.j,
               high[0], high[1], '%.1f' % upper.j,
               '%.4f' % lower.e_calc, '%.4f' % upper.e_calc,
               '%d' % lower.number, '%d' % upper.number,
               '%.4f' % wavelength]
        if not options.no_rce:
            row.extend(['' if not lower.observed else '%.6f' % lower.e_obs,
                        '' if not upper.observed else '%.6f' % upper.e_obs,
                        experimental])
        # gA and cf are written back exactly as OUTG11 printed them
        # (1PE10.3 and F6.2), so the two files can be compared field by
        # field without a rounding of ours getting in the way.
        row.extend(['%.0f' % one.intensity, '%.3E' % one.ga, '%.2f' % one.cf,
                    '%10.2e' % lower.lifetime, '%10.2e' % upper.lifetime])
        if options.ai:
            row.append('%10.2e' % upper.aa_tot)
        row.extend(['%10.2e' % upper.ar_tot, '%10.2e' % one.a,
                    '%5.3f' % one.branch])
        if options.model in (LTE, LTEABS):
            row.extend(['' if lower.level.g is None else '%.3f' % lower.level.g,
                        '' if upper.level.g is None else '%.3f' % upper.level.g])
        out.append(row)
    return out


def _auger_rows(case, order, channels, coupling, options):
    """The Auger electron spectrum: one row per level and open channel."""
    out = []
    rows = []
    for one in order:
        if one.aa_tot <= 0:
            continue
        for channel in channels:
            rate = one.channels.get(channel.name, 0.0)
            if rate <= 0:
                continue
            intensity = (2 * one.j + 1) * rate / one.width
            rows.append((one, channel, rate, intensity))
    biggest = max([r[3] for r in rows], default=0.0)
    if biggest > 0:
        rows = [(a, b, c, d * options.max_intensity / biggest)
                for a, b, c, d in rows]
    for one, channel, rate, intensity in sorted(
            rows, key=lambda r: (r[0].e_calc - r[1].energy)):
        if intensity < 0.5:
            continue
        # e_calc and channel.energy are in the run's unit, 1000 cm-1,
        # so *1000 makes cm-1 and EV2CM (cm-1 per eV) makes eV.  Until
        # 2026-08-21 this divided by 8.0655 -- eV per 1000 cm-1 -- and
        # the column was a thousand times too large; the Perl still is.
        kinetic = (one.e_calc - channel.energy) * 1000 / fundconst.EV2CM
        configuration, term = _leading(case, coupling, one)
        row = [configuration, term, '%4.1f' % one.j, channel.text,
               '%.4f' % one.e_calc, '%.4f' % channel.energy,
               '%.3f' % kinetic, '%d' % one.number]
        if not options.no_rce:
            # The channel energy is an external absolute, so this
            # difference needs an absolute level energy and not merely
            # an observed one: a group level's offset has nothing to
            # cancel against here.
            row.extend(['' if not one.observed else '%.6f' % one.e_obs,
                        '' if not one.absolute else '%.3f'
                        % (abs(one.e_obs - channel.energy) * 1000
                           / fundconst.EV2CM)])
        row.extend(['%.0f' % intensity, '%9.2e' % rate,
                    '%10.2e' % one.lifetime])
        out.append(row)
    return out


def _auger_header(options):
    out = ['up lev', '', 'J_up', 'AI_channel', 'E_up', 'E_lo', 'dE(eV)',
           'id_up']
    if not options.no_rce:
        out.extend(['Eup_exp', 'dE_exp(eV)'])
    out.extend(['Icalc', 'AI_rate(s-1)', 'lft(ns)'])
    return out


def _cascade(order, transitions):
    """The cascade matrix: A(i -> j) for every pair, levels down the side."""
    numbers = [one.number for one in order]
    rows = [['CM'] + ['%d' % n for n in numbers]]
    by_upper = {}
    for one in transitions:
        by_upper.setdefault(one.upper.number, {})[one.lower.number] = one.a
    for one in order:
        found = by_upper.get(one.number, {})
        rows.append(['%d' % one.number]
                    + ['%10.3e' % found[n] if n in found else '0'
                       for n in numbers])
    return rows


# --------------------------------------------------------------------------
# putting it out
# --------------------------------------------------------------------------

def _write(path, header, rows, delimiter='\t', preamble=None):
    with open(path, 'w') as out:
        if preamble:
            out.write(delimiter.join(preamble) + '\n')
        out.write(delimiter.join(header) + '\n')
        for row in rows:
            out.write(delimiter.join(row).rstrip(delimiter) + '\n')


def _preamble(name, options):
    """The first line of a line list: what it is and how it was scaled."""
    if options.model in (LTE, LTEABS):
        what = ('LTE absorption' if options.model == LTEABS else 'LTE')
        what += ' eff. temperature %g kK' % options.temperature
    elif options.model == AUGER:
        what = 'AUGER electron spectrum'
    elif options.model == BF:
        what = 'BF:%g' % options.flight
    else:
        what = options.model
    return ['Spectrum', name, '', what, '', '', '', 'Max. intens.', '',
            '%g' % options.max_intensity]


def _par_arguments(argv):
    """Split `--par <file>` and `--write-par` off the command line."""
    rest, named, fresh = [], None, False
    argv = list(argv)
    while argv:
        arg = argv.pop(0)
        if arg in ('--write-par', '--write_par'):
            fresh = True
        elif arg in ('--par', '--par-file'):
            if not argv:
                raise Usage('--par needs a file name')
            named = argv.pop(0)
        elif arg.startswith('--par='):
            named = arg.partition('=')[2]
        else:
            rest.append(arg)
    return rest, named, fresh


def _sample(path):
    """Write a fresh, fully commented par file and say where it is."""
    path = parfile.write(path, SPECS, title=TITLE)
    print('%s written.  Every parameter is explained in it.  Fill in the '
          'five that have no' % path)
    print('default -- levels, lines, model, temperature and max_intensity '
          '-- then run conv_out again.')
    return path


def _parameters(argv):
    """The par file's settings, or nothing, or a written sample.

    Returns `(settings, wrote_a_sample)`.  The sample is written only
    when there is no par file AND no command line at all: that is a
    first run, with nothing to go on and nothing to lose.  A command
    line that is merely incomplete is a mistake to be reported, not an
    invitation to write a file into somebody's working directory --
    which is also what keeps the reference runs in WORK/ clean.
    """
    argv, named, fresh = _par_arguments(argv)
    if fresh:
        _sample(Path(PAR_FILE))
        return argv, {}, True
    path = parfile.find(PAR_FILE, '.', named, shipped=False)
    if path is None:
        if not argv:
            _sample(Path(PAR_FILE))
            return argv, {}, True
        return argv, {}, False
    print('Parameters from %s' % path)
    return argv, parfile.read(path, SPECS), False


def main(argv):
    argv, par, stop = _parameters(argv)
    if stop:
        return 0
    options = Options(argv, par)
    case = case_mod.Case.read('.', couplings=options.couplings,
                              with_rce=not options.no_rce, drop=options.drop)
    print('%s: %d J blocks, %d levels'
          % (', '.join(options.couplings), len(case.blocks),
             len(case.identified)))

    found = _collect(case, options)
    channels = _read_ai(case, found, options)
    transitions = _read_transitions(case, found, options)
    print('%d transitions' % len(transitions))

    for one in transitions:
        one.wavelength = _wavelength(one.upper, one.lower)
        one.rcg_wavelength = _rcg_wavelength(one.upper, one.lower)
    transitions = [one for one in transitions
                   if one.wavelength and one.rcg_wavelength]
    _intensities(transitions, options)

    order = _number_levels(list(found.values()), options)
    rows = []
    for coupling in options.couplings:
        rows.extend(_level_rows(case, order, coupling, options))
    _write(options.levels, _level_header(options), rows)
    print('wrote %s (%d levels)' % (options.levels, len(rows)))

    coupling = options.couplings[0]
    name = spectrum.name(folder='.')
    if options.model == CM:
        matrix = _cascade(order, transitions)
        _write(options.lines, matrix[0], matrix[1:], delimiter=',')
    elif options.model == AUGER:
        _write(options.lines, _auger_header(options),
               _auger_rows(case, order, channels, coupling, options),
               preamble=_preamble(name, options))
    else:
        _write(options.lines, _line_header(options),
               _line_rows(case, transitions, coupling, options),
               preamble=_preamble(name, options))
    print('wrote %s' % options.lines)

    if options.iden:
        iden.write(case, order, transitions, options)
    return 0


def _read_ai(case, found, options):
    """Hang the autoionization rates of the two files on the levels."""
    if not options.ai:
        return []
    channels = []
    for path, kset in zip(options.ai, sorted({k for k, _ in case.blocks})):
        theirs, rows = ai.read_ai(path)
        channels = channels or theirs
        here = [one for one in case.identified if one.kset == kset]
        if not here:
            here = levels.unfitted(case.blocks, kset)
        paired = ai.pair(rows, here, case.blocks, path)
        for key, row in paired.items():
            one = found.get(key)
            if one is None:
                continue
            one.aa_tot = row.total
            one.e_ai = row.energy
            one.dr = row.dr(theirs)
            one.channels = dict(row.rates)
    return channels


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('conv_out: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage,')[1]))
        sys.exit(2)
    except parfile.ParFileError as problem:
        sys.stderr.write('conv_out: %s\n' % problem)
        sys.exit(2)
