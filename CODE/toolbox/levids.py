"""Level identifiers, from an experimental list into conv_out's tables.

A level of a published analysis has a name that outlives the
calculation it came from -- the identifier a database keeps it under,
`059003.000065` in the author's sheets.  `conv_out` cannot know it: its
level list is built from a finished RCG/RCE run, and the only name a
level has there is `Lev_#`, which is a position in this run's ordering
and means nothing outside it.  So the identifier has to be carried in
from the file that holds it, and the only thing the two files share is
the level itself.

They are tied together the way `subst_exp_levs` ties an experimental
list to a deck when the list carries no composition: on **(parity, J,
energy)**.  That is weaker evidence than an eigenvector, and it is
enough here for a reason it would not be there -- the observed energies
in `conv_out`'s `Ee` column came *from* a list like this one, through
the RCE deck, so the two sides are not being compared but recognised.
The tolerance is there for the case where the two files are different
vintages of the analysis and a level has moved by a hundredth of a
wavenumber since.

The tie is made once, against the level list, and the transition list
inherits it: `lid1` and `lid2` are `Lev_#`, so the identifiers of a
line's two levels are looked up rather than matched again.

Nothing here is specific to the author's sheet.  The source is any
tab-delimited file whose identifier, energy, J and parity columns are
named on the command line, by number or by header name, exactly as
`altlevsrc` takes them -- and this module borrows that reader's field
handling so that the two tools read one file the same way.
"""

from . import altlevsrc, files
from .altlevsrc import UNITS

#: The identifier column's name in the source, and its name in the
#: level list this writes.  The source's may be overridden; the
#: target's is fixed, because whatever reads the level list next looks
#: for it by that name.
SOURCE_ID = 'level_id'
LEVEL_ID = 'level_id'

#: What the identifiers are called in the transition list: the lower
#: level's and the upper level's, beside the `lid1`/`lid2` they were
#: found through.
LINE_IDS = ('id1', 'id2')

#: The columns this reads in the two files conv_out writes.
LEVEL_KEYS = ('Lev_#', 'Ee', 'J', 'par')
LINE_KEYS = ('lid1', 'lid2')

#: How far apart, in the source's own unit, two energies may be and
#: still be one level.  One wavenumber: far enough to survive a level
#: optimization being redone, far short of the spacing of two levels of
#: one J and parity in any spectrum this package is used on.
TOLERANCE = 1.0


class LevelIdError(Exception):
    """A file, or a column of it, that cannot be used."""


# --------------------------------------------------------------------------
# the source
# --------------------------------------------------------------------------

class Spec:
    """Which column of the source file is which.

    The same three ways of naming a column as `altlevsrc.Spec`: `None`
    for a column that is not there, a 1-based number, or a name looked
    up in the header line.  `idcol` is the one with a default, because
    the file that has an identifier column usually calls it `level_id`.
    """

    NAMES = ('idcol', 'ecol', 'jcol', 'pcol', 'excludecol')

    def __init__(self, **given):
        for name in self.NAMES:
            setattr(self, name, given.pop(name, None))
        if self.idcol is None:
            self.idcol = SOURCE_ID
        self.units = given.pop('units', 'cm-1')
        self.skip = given.pop('skip', None)         # None: work it out
        self.tolerance = given.pop('tolerance', TOLERANCE)
        if given:
            raise LevelIdError('%s is not a source option'
                               % ', '.join(sorted(given)))
        if str(self.units).lower() not in UNITS:
            raise LevelIdError('%s is not a unit; %s is what there is'
                               % (self.units, ', '.join(sorted(UNITS))))
        try:
            self.tolerance = float(self.tolerance)
        except (TypeError, ValueError):
            raise LevelIdError('--tol wants a number, not %r' % self.tolerance)
        if self.tolerance < 0:
            raise LevelIdError('--tol is not negative')

    @property
    def scale(self):
        return UNITS[str(self.units).lower()]


class SourceLevel:
    """One row of the source: an identifier and what it belongs to."""

    def __init__(self, ident, energy, j, parity, line):
        self.ident = ident
        self.energy = energy                # in the run's unit, 10^3 cm-1
        self.j = j                          # 2J, as an integer
        self.parity = parity                # 'e', 'o' or ''
        self.line = line
        self.matched = None


class Report:
    """What a run wants said about it."""

    def __init__(self):
        self.path = ''
        self.rows = 0
        self.skipped = 0
        self.levels = 0
        self.excluded = 0
        self.no_id = 0
        self.bad_j = []
        self.bad_parity = []
        self.repeated = {}          # {identifier: rows carrying it}
        self.parity_used = True
        self.matched = 0
        self.calculated = 0         # matched to a level with no observed E
        self.lost = []              # source levels no target level took
        self.blank = 0              # target levels left without an identifier
        self.lines = 0
        self.orphans = {}           # {Lev_# in the line list: rows}

    def write(self, out, tolerance, units):
        scale = UNITS[str(units).lower()]
        out.write('%s: %d rows read' % (self.path, self.rows))
        if self.skipped:
            out.write(', %d header line%s skipped'
                      % (self.skipped, '' if self.skipped == 1 else 's'))
        out.write('; %d levels' % self.levels)
        if self.excluded:
            out.write(', %d excluded by the exclude column' % self.excluded)
        out.write('.\n')
        if self.no_id:
            out.write('   %d row%s no identifier and %s left out.\n'
                      % (self.no_id, ' has' if self.no_id == 1 else 's have',
                         'was' if self.no_id == 1 else 'were'))
        if not self.parity_used:
            out.write('   No parity column was given, so the match is on '
                      '(J, energy) alone.\n')
        for ident, count in sorted(self.repeated.items()):
            out.write('   Warning: identifier %s is on %d rows; it is meant '
                      'to be unique.\n' % (ident, count))
        for number, text in self.bad_parity[:5]:
            out.write('   Warning: line %d: %r is not a parity.\n'
                      % (number, text))
        for number, text in self.bad_j[:5]:
            out.write('   Warning: line %d: %r is not a J.\n'
                      % (number, text))
        out.write('   %d matched within %g %s; %d source level%s unmatched, '
                  '%d target level%s left blank.\n'
                  % (self.matched, tolerance, units,
                     len(self.lost), '' if len(self.lost) == 1 else 's',
                     self.blank, '' if self.blank == 1 else 's'))
        if self.calculated:
            out.write('   Warning: %d of the matches are to a level whose Ee '
                      'is the calculated energy and not an observed one.\n'
                      % self.calculated)
        for level in self.lost[:20]:
            out.write('   unmatched: %s, parity %s, J = %s, E = %g %s.\n'
                      % (level.ident, level.parity or '?', j_text(level.j),
                         level.energy / scale, units))
        if len(self.lost) > 20:
            out.write('   ... and %d more.\n' % (len(self.lost) - 20))
        for number, count in sorted(self.orphans.items())[:5]:
            out.write('   Warning: the transitions name level %d, which the '
                      'level list does not have (%d row%s).\n'
                      % (number, count, '' if count == 1 else 's'))


def read_source(path, spec):
    """Read the identifiers out of a tab-delimited level list.

    Returns `(levels, report)`.  A row with no identifier names no
    level to carry across and is counted rather than matched.
    """
    rows = altlevsrc._rows(path)
    if not rows:
        raise LevelIdError('%s holds no rows' % path)

    report = Report()
    report.path = path
    header, first = altlevsrc._header(rows, spec)
    report.skipped = first

    columns = {}
    for name in Spec.NAMES:
        columns[name] = altlevsrc._column(getattr(spec, name), name,
                                          header, path)
    for name in ('idcol', 'ecol', 'jcol'):
        if columns[name] is None:
            raise LevelIdError(
                '--%s is needed to say which column of %s holds its %s'
                % (name, path,
                   {'idcol': 'level identifiers', 'ecol': 'energy',
                    'jcol': 'J'}[name]))
    report.parity_used = columns['pcol'] is not None

    out, seen = [], {}
    for number, row in enumerate(rows[first:], first + 1):
        if not any(field.strip() for field in row):
            continue
        report.rows += 1
        if altlevsrc._excluded(row, columns['excludecol']):
            report.excluded += 1
            continue
        ident = altlevsrc._field(row, columns['idcol'])
        if not ident:
            report.no_id += 1
            continue
        try:
            energy = float(altlevsrc._field(row, columns['ecol'])
                           .replace(',', '')) * spec.scale
        except ValueError:
            continue
        parity = ''
        if columns['pcol'] is not None:
            parity = altlevsrc._parity(altlevsrc._field(row, columns['pcol']))
            if parity is None:
                report.bad_parity.append(
                    (number, altlevsrc._field(row, columns['pcol'])))
                continue
        j = j_key(altlevsrc._field(row, columns['jcol']))
        if j is None:
            report.bad_j.append(
                (number, altlevsrc._field(row, columns['jcol'])))
            continue
        seen[ident] = seen.get(ident, 0) + 1
        out.append(SourceLevel(ident, energy, j, parity, number))
        report.levels += 1

    report.repeated = dict((ident, count) for ident, count in seen.items()
                           if count > 1)
    return out, report


def j_key(text):
    """2J as an integer, from any of the ways a file writes a J.

    `1/2`, `0.5` and `.5` are one level; the doubled value is what they
    are compared as, so that no equality test is ever made between two
    floats.  `None` when the field is not a J at all -- which includes
    the blank field of a row standing for a whole term, and such a row
    names no single level to give an identifier to.
    """
    value = str(text).strip().strip('"').strip().rstrip('*')
    if not value:
        return None
    try:
        if '/' in value:
            top, _, bottom = value.partition('/')
            number = float(top) / float(bottom)
        else:
            number = float(value)
    except (ValueError, ZeroDivisionError):
        return None
    doubled = int(round(2 * number))
    if abs(2 * number - doubled) > 1e-6 or doubled < 0:
        return None
    return doubled


def j_text(j):
    """2J back to the way a level list writes it: `3`, `3/2`."""
    return '%d' % (j // 2) if j % 2 == 0 else '%d/2' % j


# --------------------------------------------------------------------------
# a table conv_out wrote
# --------------------------------------------------------------------------

class Table:
    """A tab-delimited table with a header line near the top.

    The line list carries a preamble line above its header and the
    level list does not, so the header is found by what it names rather
    than by where it is.  Everything above it is kept and written back
    untouched.

    Rows are short: `conv_out` strips the trailing delimiters off every
    row it writes, so a row with nothing in its last columns has fewer
    fields than the header.  That is kept -- a row is padded only as far
    as an insertion needs it, and the trailing delimiters come off again
    on the way out -- so that a table this tool has not changed a field
    of comes back out as the bytes it went in as.
    """

    def __init__(self, path, delimiter='\t'):
        self.path = str(path)
        self.delimiter = delimiter
        self.ending = files.line_ending(self.path)
        with open(self.path, 'r', errors='replace') as handle:
            text = handle.read()
        lines = text.split('\n')
        if lines and lines[-1] == '':
            lines.pop()
        self.raw = [line.rstrip('\r') for line in lines]
        self.first = 0
        self.header = []
        self.rows = []

    @classmethod
    def read(cls, path, keys, delimiter='\t', within=5):
        """Read `path` and find the header line that names every key."""
        table = cls(path, delimiter)
        wanted = [key.lower() for key in keys]
        for i, line in enumerate(table.raw[:within]):
            named = [field.strip().lower() for field in line.split(delimiter)]
            if all(key in named for key in wanted):
                table.first = i
                table.header = line.split(delimiter)
                table.rows = [row.split(delimiter)
                              for row in table.raw[i + 1:]]
                return table
        raise LevelIdError(
            '%s has no header line naming %s in its first %d lines; is it '
            'the file conv_out wrote?'
            % (path, ', '.join(keys), within))

    def index(self, name):
        """Where a column is, by name, or `None`."""
        wanted = name.strip().lower()
        for i, field in enumerate(self.header):
            if field.strip().lower() == wanted:
                return i
        return None

    def field(self, row, index):
        if index is None or index >= len(row):
            return ''
        return row[index].strip()

    def put(self, row, index, value):
        while len(row) <= index:
            row.append('')
        row[index] = value

    def column(self, name, before=None):
        """Where a column is to be written, making it if it is not there.

        A column of that name already in the table keeps its place and
        is emptied, which is what makes a second run of the tool put
        the identifiers it has just matched where the first run's are
        rather than beside them.  A new one goes in front of the column
        named by `before`, or at the end if there is none.
        """
        at = self.index(name)
        if at is not None:
            for row in self.rows:
                if at < len(row):
                    row[at] = ''
            return at
        at = self.index(before) if before else None
        return self._insert(name, len(self.header) if at is None else at)

    def column_after(self, name, other):
        """The column `name`, made just after `other` if it is not there.

        By index and not by the name of whatever stands there now: a
        line list's header has blank fields in it -- the term columns
        beside `level1` and `level2` are unnamed -- and a blank is no
        name to insert in front of.
        """
        if self.index(name) is not None:
            return self.column(name)
        at = self.index(other)
        if at is None:
            raise LevelIdError('%s has no %s column' % (self.path, other))
        return self._insert(name, at + 1)

    def _insert(self, name, at):
        """Make a column at `at`, moving every row's fields along."""
        self.header.insert(at, name)
        for row in self.rows:
            if not any(field.strip() for field in row):
                continue
            while len(row) < at:
                row.append('')
            row.insert(at, '')
        return at

    def write(self, backup=None):
        out = list(self.raw[:self.first])
        for row in [self.header] + self.rows:
            out.append(self.delimiter.join(row).rstrip(self.delimiter))
        files.write_deck(self.path, out, self.ending, backup=backup)


# --------------------------------------------------------------------------
# the level list
# --------------------------------------------------------------------------

class TargetLevel:
    """One level of the level list, and every row that prints it.

    A level appears once per coupling block -- the same `Lev_#` and the
    same energy under a different `Cpl` -- so it is one level here and
    the identifier it is given goes on all of its rows.
    """

    def __init__(self, number, energy, j, parity, calculated):
        self.number = number
        self.energy = energy
        self.j = j
        self.parity = parity
        self.calculated = calculated
        self.rows = []
        self.ident = ''


def read_levels(table):
    """The levels of `table`: `(levels, {Lev_#: level})`.

    `e_c` holds a star on a level whose `Ee` is the calculated energy
    because no observed one was known, and that is worth carrying: such
    a level is a poorer thing to match an observed energy against, and
    the match to it is reported.
    """
    at = dict((name, table.index(name)) for name in LEVEL_KEYS)
    marker = table.index('e_c')
    out, by_number = [], {}
    for row in table.rows:
        if not any(field.strip() for field in row):
            continue
        try:
            number = int(table.field(row, at['Lev_#']))
            energy = float(table.field(row, at['Ee']))
        except ValueError:
            continue
        level = by_number.get(number)
        if level is None:
            level = TargetLevel(number, energy,
                                j_key(table.field(row, at['J'])),
                                altlevsrc._parity(table.field(row, at['par']))
                                or '',
                                table.field(row, marker) == '*')
            by_number[number] = level
            out.append(level)
        level.rows.append(row)
    return out, by_number


def match(source, targets, tolerance):
    """Tie source levels to target levels, the closest pair first.

    Every pair within the tolerance is offered, and they are taken in
    order of how far apart they are, each level of either side being
    used once.  Doing it that way rather than level by level is what
    keeps two neighbouring levels of one J from both taking the target
    that lies between them: the nearer one takes it and the other is
    reported unmatched, instead of the order of the file deciding.

    An `Ee` that is the calculated energy loses a tie to an observed
    one, because an observed energy is what a source level is another
    vintage of.

    Returns the source levels that found nothing.
    """
    groups = {}
    for level in targets:
        groups.setdefault(level.j, []).append(level)

    offers = []
    for i, one in enumerate(source):
        for level in groups.get(one.j, ()):
            if one.parity and level.parity and one.parity != level.parity:
                continue
            apart = abs(one.energy - level.energy)
            if apart <= tolerance:
                offers.append((apart, level.calculated, i, level.number,
                               one, level))
    offers.sort(key=lambda offer: offer[:4])

    for _, _, _, _, one, level in offers:
        if one.matched is not None or level.ident:
            continue
        one.matched = level
        level.ident = one.ident
    return [one for one in source if one.matched is None]


def insert_levels(table, targets, report):
    """Write the identifiers into the level list, beside `Ee`."""
    at = table.column(LEVEL_ID, before='Ee')
    for level in targets:
        if not level.ident:
            report.blank += 1
            continue
        report.matched += 1
        if level.calculated:
            report.calculated += 1
        for row in level.rows:
            table.put(row, at, level.ident)
    return at


# --------------------------------------------------------------------------
# the transition list
# --------------------------------------------------------------------------

def insert_lines(table, by_number, report):
    """Write the two levels' identifiers beside `lid1` and `lid2`.

    `id1` goes in after `lid2` and `id2` after `id1`, so the four
    columns end up in the order a reader expects and the second
    insertion moves neither of the `lid` columns nor the first.
    """
    table.column_after(LINE_IDS[0], LINE_KEYS[1])
    table.column_after(LINE_IDS[1], LINE_IDS[0])
    into = [table.index(name) for name in LINE_IDS]
    lid = [table.index(name) for name in LINE_KEYS]

    for row in table.rows:
        if not any(field.strip() for field in row):
            continue
        report.lines += 1
        for which in (0, 1):
            try:
                number = int(table.field(row, lid[which]))
            except ValueError:
                continue
            level = by_number.get(number)
            if level is None:
                report.orphans[number] = report.orphans.get(number, 0) + 1
                continue
            table.put(row, into[which], level.ident)


# --------------------------------------------------------------------------
# the IDEN2 translation table
# --------------------------------------------------------------------------
#
# IDEN2 numbers its levels itself, from the highest energy down, and that
# number is the only name a level has in `ENLEV.DAT` and `TRANS.DAT`.  It
# means nothing outside one pair of those files -- rebuild them and every
# number moves -- so anything that has to say the same thing to IDEN2 and
# to the rest of the analysis needs the two names side by side.  That is
# what this table is, and `write_IDEN_level_ids` is what writes it.

#: The file the table is kept in, and what its two columns are called.
IDEN_IDS = 'IDEN_level_ids.txt'
IDEN_ID = 'IDEN_id'


def write_iden_ids(pairs, path=IDEN_IDS, ending='\n'):
    """Write the `(identifier, IDEN2 number)` pairs as a two-column table.

    LF, and the tab, are what the other tab-delimited tables of the
    package are written with -- calc_TP's line list, conv_out's two --
    and every reader here takes either ending.
    """
    lines = ['%s\t%s' % (SOURCE_ID, IDEN_ID)]
    lines.extend('%s\t%d' % (ident, number) for ident, number in pairs)
    with open(path, 'w', newline='') as handle:
        handle.write(''.join(line + ending for line in lines))
    return path


def read_iden_ids(path=IDEN_IDS):
    """The table back as `{identifier: IDEN2 level number}`.

    A repeated identifier, or a repeated level number, is an error and
    not a last-one-wins: the table is a bijection, and a caller that
    translates a level list through a broken one would quietly put two
    levels' energies into one record.
    """
    try:
        rows = altlevsrc._rows(path)
    except OSError:
        raise LevelIdError(
            '%s is not here; run write_IDEN_level_ids to make it' % path)
    if not rows:
        raise LevelIdError('%s holds no rows' % path)

    header = [field.strip().lower() for field in rows[0]]
    try:
        left = header.index(SOURCE_ID.lower())
        right = header.index(IDEN_ID.lower())
    except ValueError:
        raise LevelIdError('%s has no header line naming %s and %s'
                           % (path, SOURCE_ID, IDEN_ID))

    out, taken = {}, {}
    for at, row in enumerate(rows[1:], 2):
        if not any(field.strip() for field in row):
            continue
        ident = altlevsrc._field(row, left)
        text = altlevsrc._field(row, right)
        if not ident and not text:
            continue
        try:
            number = int(text)
        except ValueError:
            raise LevelIdError('%s line %d: %r is not an IDEN2 level number'
                               % (path, at, text))
        if ident in out:
            raise LevelIdError('%s names %s twice' % (path, ident))
        if number in taken:
            raise LevelIdError('%s gives IDEN2 level %d to both %s and %s'
                               % (path, number, taken[number], ident))
        out[ident] = number
        taken[number] = ident
    if not out:
        raise LevelIdError('%s holds no levels' % path)
    return out
