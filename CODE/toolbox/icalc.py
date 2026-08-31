"""The line list a classification program wants out of `calc_TP`.

`calc_TP` writes a wide table -- fifty columns and more once the trial
data are stored -- and a line classification code needs six of them:
which two published levels the transition joins, how far apart they
are, how strong it is, how well that strength is known, and where the
upper level sits.  This is the extractor for those, and, if a
temperature is given, for a synthetic intensity built out of them.

WHICH LINES ARE TAKEN

Only the ones between two levels of the published analysis: `id1` and
`id2`, which `insert_level_ids` writes on the lines it could match,
must both be there.  A line to a level the analysis does not hold has
no identifier to be classified under, and a line to a level the fit did
not observe has no experimental energy to be classified by.  In
practice those are the same lines -- `calc_TP` fills `E_low_exp` and
`E_up_exp` only for observed levels, and an identifier comes from the
same list the observed energies did -- but both are checked, because a
missing experimental energy would silently become a wavenumber of zero
or of the wrong thing, and the count of each is reported.

THE UNITS

`calc_TP` carries its energies in kK, the unit RCN2 and RCG print
parameters in; `rwn` and `Eup` are written here in cm-1, which is what
a classification list is held in and what the Boltzmann factor below
needs.  So both are the file's value times 1000.

WHAT gA IS

The line list's `A(s-1)` column is A, the transition probability of the
upper level; `gA` is that times the statistical weight 2*J2+1 of the
upper level, and `J2` is the column beside it.  `std_A%` is a
*relative* standard deviation, so it is the relative standard deviation
of gA as well and is carried across unchanged -- including the blank
`calc_TP` leaves on a line the trials printed once or not at all, which
must not be read as a line known exactly.

THE SYNTHETIC INTENSITY

Given a temperature in eV and a log-offset, the seventh column is

    Icalc = exp(-Eup/T + d_off) * gA / rwn

with Eup, T and rwn all in cm-1 -- a Boltzmann population of the upper
level times the transition rate, per unit wavenumber.  `d_off` is the
one free normalisation: it absorbs the partition function, the column
density and the choice of unit for rwn all at once, which is why it is
the caller's to set rather than something computed here.  Without both
of them the column is not written at all, rather than written empty.
"""

import math

from . import files
from . import fundconst
from . import levids


class IcalcError(Exception):
    """The extraction cannot be done as asked."""


#: The default names, in and out.
SOURCE = 'tp_E1.txt'
OUTPUT = 'Icalc.txt'

#: The columns written, in order.  The seventh is conditional.
COLUMNS = ('id1', 'id2', 'rwn', 'gA', 'u%gA', 'Eup')
INTENSITY = 'Icalc'

#: What has to be in the source's header.  `std_A%` is not here: a
#: `calc_TP` run with no trials has no uncertainties to print, and that
#: is a table this can still be used on -- `u%gA` comes out blank.
NEEDED = ('id1', 'id2', 'J2', 'E_low_exp', 'E_up_exp', 'A(s-1)')
SPREAD = 'std_A%'

#: kK -> cm-1.
KK = 1000.0


class Report:
    """What one run passed over, and why."""

    def __init__(self):
        self.read = 0           # data rows in the source
        self.written = 0        # rows in the output
        self.unidentified = 0   # no id1, or no id2
        self.unobserved = 0     # identified, but no experimental energy
        self.unstrengthed = 0   # identified, but no A
        self.degenerate = 0     # identified, but the two energies are equal

    def write(self, out):
        out.write('%d transition%s read, %d written.\n'
                  % (self.read, '' if self.read == 1 else 's', self.written))
        for count, why in ((self.unidentified,
                            'without an identifier at one end or both'),
                           (self.unobserved,
                            'without an experimental energy at one end'),
                           (self.unstrengthed, 'without a transition rate'),
                           (self.degenerate,
                            'whose two experimental energies are equal')):
            if count:
                out.write('  %d passed over, %s.\n' % (count, why))


def _fixed(value):
    """A number written plainly, without a tail of zeros.

    The energies come in with four decimals in kK, so a thousand times
    one of them is exact to the first decimal; `%.4f` keeps whatever the
    source had and the trimming takes the padding back off, so a value
    written as 138.3427 comes out as 138342.7 and not as 138342.7000.
    """
    text = '%.4f' % value
    if '.' in text:
        text = text.rstrip('0').rstrip('.')
    return text or '0'


def _number(text, name, row):
    try:
        return float(text)
    except ValueError:
        raise IcalcError('%s is %r on the transition %s, and that is not a '
                         'number' % (name, text, ' '.join(row[:7])))


def read(path):
    """The source table, checked for the columns this needs."""
    try:
        table = levids.Table.read(path, ('id1', 'id2'))
    except levids.LevelIdError:
        raise IcalcError(
            '%s has no header line naming id1 and id2 in its first few '
            'lines.  Is it a calc_TP line list that insert_level_ids has '
            'been run on?' % path)
    missing = [name for name in NEEDED if table.index(name) is None]
    if missing:
        raise IcalcError('%s has no %s column%s; is it a calc_TP line list?'
                         % (path, ', '.join(missing),
                            '' if len(missing) == 1 else 's'))
    return table


def extract(table, temperature=None, offset=None):
    """The output rows, and the account of what was left out.

    `temperature` is in cm-1 -- the caller converts, because it is the
    caller that was given eV -- and both it and `offset` have to be
    there for the intensity column to be made.
    """
    wanted = temperature is not None and offset is not None
    at = dict((name, table.index(name)) for name in NEEDED)
    at[SPREAD] = table.index(SPREAD)

    report = Report()
    rows = []
    for row in table.rows:
        if not any(field.strip() for field in row):
            continue
        report.read += 1

        id1 = table.field(row, at['id1'])
        id2 = table.field(row, at['id2'])
        if not id1 or not id2:
            report.unidentified += 1
            continue
        low = table.field(row, at['E_low_exp'])
        up = table.field(row, at['E_up_exp'])
        if not low or not up:
            report.unobserved += 1
            continue
        rate = table.field(row, at['A(s-1)'])
        if not rate:
            report.unstrengthed += 1
            continue

        e_low = _number(low, 'E_low_exp', row) * KK
        e_up = _number(up, 'E_up_exp', row) * KK
        rwn = abs(e_up - e_low)
        if not rwn:
            # Nothing to divide by, and no wavelength either.  calc_TP
            # drops such a line on its calculated energies; this is the
            # same guard on the experimental ones.
            report.degenerate += 1
            continue
        weight = 2.0 * _number(table.field(row, at['J2']), 'J2', row) + 1.0
        ga = weight * _number(rate, 'A(s-1)', row)

        out = [id1, id2, _fixed(rwn), '%.5E' % ga,
               table.field(row, at[SPREAD]), _fixed(e_up)]
        if wanted:
            out.append('%.5E' % (math.exp(-e_up / temperature + offset)
                                 * ga / rwn))
        rows.append(out)
        report.written += 1
    return rows, report


def write(path, rows, intensity, ending='\n'):
    """The tab-delimited output, header first."""
    header = list(COLUMNS) + ([INTENSITY] if intensity else [])
    lines = ['\t'.join(header)] + ['\t'.join(row) for row in rows]
    return files.write_deck(path, lines, ending)


def to_wavenumbers(ev):
    """A temperature in eV as one in cm-1."""
    return ev * fundconst.EV2CM
