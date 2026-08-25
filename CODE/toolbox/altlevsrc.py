"""An experimental level list that is not an RCE deck.

`subst_exp_levs` was written to carry observed energies from one RCE
deck into another, and that is still the case it is best at: both sides
speak the same eigenvectors, so a level is recognised by what it is made
of rather than by where it sits.  But the observed energies do not start
life in an RCE deck.  They come out of a level-optimization run, or off
a spreadsheet somebody keeps a spectrum's bookkeeping in, and those
files agree on almost nothing -- not the columns, not the units, not
whether a composition is given at all, and where it is given it is in
**intensity** per cent rather than the signed amplitude per cent RCE
prints.

So this module reads *any* tab-delimited level list, told by option
which column is which, and builds the same `Deck` the matcher already
takes.  Nothing here is specific to LOPT, which is what the branch it
replaces was written against and had hard-coded down to the field
number.

Three things are worth knowing before reading the code.

**The key space is the target deck's own configuration labels.**  RCE
writes a configuration in six columns and abbreviates from the left to
get there -- `4f.5d2` is `f5d2`, `5d.6s.6p` is `d6s6p` -- and no rule
recovers the missing digits in general.  It does not have to: the deck
names its own configurations, that list is short, and a source
configuration written out in full matches the deck label its compact
form *ends with*.  The longest such label wins, and a configuration no
label matches is reported rather than guessed at.

**The comparison is in intensity per cent, not amplitude.**  A source
sheet quotes what a published table quotes, so the target's amplitudes
are squared over 100 to meet it and the sign question does not arise.
`FLOOR` still cuts the target's components in amplitude space, so the
same components survive as in a deck-to-deck run.

**A row whose leading component is not in the target's coupling gets no
vector at all.**  Experimental designations are chosen for whichever
scheme makes the level look pure, so one file can hold LS, JK and pair
labels side by side; only LS can be matched against an LS deck.  Such a
row is matched on (parity, J, energy) instead, and counted in the report
so that the share of the file that got the weaker treatment is never a
surprise.
"""

import csv
import re

from . import purity
from .explevs import Deck, Level, _j_text
from .fundconst import EV2CM

#: Source energies, in the 10^3 cm-1 RCE works in.
UNITS = {'cm-1': 0.001, 'cm1': 0.001, 'kk': 1.0, 'ev': EV2CM / 1000.0}

#: How a parity column may say which parity a level is.
ODD = ('o', 'odd', 'u', '1', '*')
EVEN = ('e', 'even', 'g', '0')

_L = 'SPDFGHIKLMNOQ'
_LS = re.compile(r'^(\d+)([%s])(\d*)\*?$' % _L, re.I)
_SHELL = re.compile(r'(\d*)([spdfghiklmno])(\d*?)'
                    r'(?=\d[spdfghiklmno]|[spdfghiklmno]|$)')
_REPEAT = re.compile(r'^(\d+[%s])\d+$' % _L, re.I)
_ORBITALS = 'spdfghiklmno'


class AltSourceError(Exception):
    """A source level list the options do not describe."""


class Spec:
    """Which column of the source file is which.

    Every field is either `None` (the column is not there), a 1-based
    column number, or a column name to be looked up in the header.  The
    reader resolves them once it has read the file, because only then
    does it know whether there is a header at all.
    """

    NAMES = ('ecol', 'jcol', 'pcol', 'excludecol', 'altcol',
             'confcol', 'termcol', 'perccol', 'conf1col', 'term1col',
             'perc2col', 'conf2col', 'term2col',
             'perc3col', 'conf3col', 'term3col')

    def __init__(self, **given):
        for name in self.NAMES:
            setattr(self, name, given.pop(name, None))
        self.units = given.pop('units', 'cm-1')
        self.skip = given.pop('skip', None)         # None: work it out
        self.altvalue = given.pop('altvalue', 'or')
        self.parity_sets = given.pop('parity_sets', None)
        if given:
            raise AltSourceError('%s is not a source option'
                                 % ', '.join(sorted(given)))
        if str(self.units).lower() not in UNITS:
            raise AltSourceError('%s is not a unit; %s is what there is'
                                 % (self.units, ', '.join(sorted(UNITS))))

    @property
    def scale(self):
        return UNITS[str(self.units).lower()]


class Report:
    """What the reader wants said about the file it has just read."""

    def __init__(self):
        self.path = ''
        self.header = []
        self.skipped = 0
        self.rows = 0
        self.levels = 0
        self.excluded = 0
        self.by_vector = 0
        self.by_energy = 0
        self.partial = 0            # a vector of one component only
        self.couplings = {}         # {scheme: rows whose leading term is it}
        self.unknown_conf = {}      # {source configuration: rows}
        self.bad_parity = []
        self.bad_j = []
        self.parities = {}          # {parity: parity set}
        self.use_conf = False       # whether the keys name a configuration

    def write(self, out):
        out.write('%s: %d rows read' % (self.path, self.rows))
        if self.skipped:
            out.write(', %d header line%s skipped'
                      % (self.skipped, '' if self.skipped == 1 else 's'))
        out.write('; %d levels' % self.levels)
        if self.excluded:
            out.write(', %d excluded by the exclude column' % self.excluded)
        out.write('.\n')
        out.write('   parity sets: %s\n'
                  % ', '.join('%s -> %d' % (parity, kset) for parity, kset
                              in sorted(self.parities.items())))
        if self.by_vector or self.by_energy:
            out.write('   %d level%s carry a composition the target can be '
                      'matched on; %d %s matched on (parity, J, energy) '
                      'alone.\n'
                      % (self.by_vector, '' if self.by_vector == 1 else 's',
                         self.by_energy,
                         'is' if self.by_energy == 1 else 'are'))
        if self.couplings:
            out.write('   leading components by coupling scheme: %s\n'
                      % ', '.join('%s %d' % (name.upper(), count)
                                  for name, count
                                  in sorted(self.couplings.items())))
        if self.partial:
            out.write('   Warning: %d of them name one component only, which '
                      'is thin evidence; check those matches.\n'
                      % self.partial)
        for label, count in sorted(self.unknown_conf.items()):
            out.write('   Warning: configuration %s (%d level%s) is in no '
                      'component of the target deck.\n'
                      % (label, count, '' if count == 1 else 's'))
        for number, text in self.bad_parity[:5]:
            out.write('   Warning: line %d: %r is not a parity.\n'
                      % (number, text))
        for number, text in self.bad_j[:5]:
            out.write('   Warning: line %d: %r is not a J.\n'
                      % (number, text))


# --------------------------------------------------------------------------
# labels
# --------------------------------------------------------------------------

def coupling(term):
    """Which scheme a term label is written in, or None.

    The label says so itself: a bracket is the JK pair label `2[3/2]`,
    a parenthesis with a comma is the jj pair label `(2,3/2)`, and a
    multiplicity followed by an orbital letter is LS.  A blank, or the
    bare parity mark a sheet writes when it has no term to give, is
    none of them.
    """
    text = term.strip().strip('"').strip()
    if not text or text == '*':
        return None
    if '[' in text:
        return 'jk'
    if '(' in text and ',' in text:
        return 'jj'
    if _LS.match(text):
        return 'ls'
    return None


def term_key(term):
    """The LS term a component is matched on: `2P*` and `2P2` -> `2P`."""
    text = term.strip().strip('"')
    text = text.split()[-1] if text.split() else ''
    text = text.rstrip('*').upper()
    found = _REPEAT.match(text)
    return found.group(1) if found else text


def compact(configuration):
    """`4f2.(3F<2>).6p` -> `4f26p`: what a deck label abbreviates."""
    text = configuration.strip().strip('"')
    return purity.strip_intermediates(text).replace('.', '')


def deck_labels(target):
    """Every configuration label the target deck's components use."""
    out = set()
    for level in target.levels():
        for key in level.vector:
            label = key.partition(',')[0]
            if label:
                out.add(label)
    return out


def resolve_configuration(configuration, labels):
    """The deck label a written-out configuration abbreviates to.

    RCE cuts a configuration from the left to fit six columns, so the
    deck's label is a tail of the full form.  The longest label that is
    one wins; nothing else is a match.
    """
    text = compact(configuration)
    if not text:
        return ''
    best = ''
    for label in labels:
        if text.endswith(label) and len(label) > len(best):
            best = label
    return best


def parity_of_label(label):
    """0 even, 1 odd, from a deck's abbreviated configuration label.

    The digits are ambiguous -- `f25d` is 4f2.5d, so the `5` is the
    next subshell's principal quantum number and not part of the f
    occupancy -- and the last digit before a letter is always the
    letter's, which is what the lookahead in `_SHELL` says.

    A label is whatever column 3 of the IN36 card called it, so it need
    not be an abbreviation at all: `f4p5` is Pr III's 4f4.5p5, and read
    as a configuration it comes out even where the deck has it odd.
    That is why `parity_sets` takes a majority rather than trusting one
    label.
    """
    odd = 0
    for _, orbital, count in _SHELL.findall(label.lower()):
        if _ORBITALS.index(orbital) % 2:
            odd += int(count) if count else 1
    return odd % 2


def parity_sets(target):
    """{'e': parity set, 'o': parity set}, from the target's own labels.

    The old reader let the first row of the source settle which parity
    set was which, because a LOPT list does not say.  The target deck
    does say -- every component names a configuration, and a
    configuration has a parity -- so it is read off there instead, and
    the source is only asked which parity each of *its* levels is.
    """
    out = {}
    for kset in sorted(target.sets):
        votes = {0: 0, 1: 0}
        for j in target.sets[kset]:
            for level in target.group(kset, j):
                label = level.leading[1].partition(',')[0]
                if label:
                    votes[parity_of_label(label)] += 1
        if not (votes[0] or votes[1]):
            continue
        # A majority: a nickname that is not an abbreviation reads as
        # the wrong parity, and one of Pr III's twenty-one does.
        out['e' if votes[0] >= votes[1] else 'o'] = kset
    return out


# --------------------------------------------------------------------------
# the target, in the same key space
# --------------------------------------------------------------------------

def normalize_target(target, use_conf=True):
    """Give every target level the vector the source is compared to.

    In intensity per cent, keyed the way a source row is keyed, with
    the components that fall on one key added together -- an experiment
    that quotes `4f3 2P` is not distinguishing the parents RCE prints
    beside it.
    """
    for level in target.levels():
        out = {}
        for key, percent in level.vector.items():
            configuration, _, term = key.partition(',')
            mine = '%s,%s' % (configuration if use_conf else '',
                              term_key(term))
            out[mine] = out.get(mine, 0.0) + percent * percent / 100.0
        level.mvector = out
    return target


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def read_alt_lev_src(path, spec, target):
    """Read a tab-delimited level list into the Deck the matcher takes.

    `target` is the deck the levels are going into: it supplies the
    configuration labels the source's are resolved against and the
    parity of each parity set.  Returns `(deck, report)`.
    """
    rows = _rows(path)
    if not rows:
        raise AltSourceError('%s holds no rows' % path)

    report = Report()
    report.path = path
    header, first = _header(rows, spec)
    report.header, report.skipped = header, first

    columns = {}
    for name in Spec.NAMES:
        columns[name] = _column(getattr(spec, name), name, header, path)
    if columns['excludecol'] is None:
        named = [field.strip().strip('"').lower() for field in header]
        if 'exclude' in named:
            columns['excludecol'] = named.index('exclude')
    for name in ('ecol', 'jcol', 'pcol'):
        if columns[name] is None:
            raise AltSourceError(
                '%s does not parse as an RCE deck, so --%s is needed to say '
                'which column its %s is in'
                % (path, name.replace('_', '-'),
                   {'ecol': 'energy', 'jcol': 'J', 'pcol': 'parity'}[name]))

    wants_conf = any(columns[name] is not None
                     for name in ('confcol', 'conf1col', 'conf2col',
                                  'conf3col'))
    labels = deck_labels(target) if wants_conf else set()
    use_conf = bool(labels)
    report.use_conf = use_conf
    report.parities = _parity_sets(spec, target)

    deck = Deck(path)
    for number, row in enumerate(rows[first:], first + 1):
        if not any(field.strip() for field in row):
            continue
        report.rows += 1
        if _excluded(row, columns['excludecol']):
            report.excluded += 1
            continue
        try:
            value = float(_field(row, columns['ecol']).replace(',', '')) \
                * spec.scale
        except ValueError:
            continue
        parity = _parity(_field(row, columns['pcol']))
        kset = report.parities.get(parity) if parity else None
        if kset is None:
            report.bad_parity.append((number, _field(row, columns['pcol'])))
            continue

        vector, scheme = _vector(row, columns, spec, labels, use_conf, report)
        report.couplings[scheme] = report.couplings.get(scheme, 0) + 1
        if vector:
            report.by_vector += 1
            if len(vector) == 1:
                report.partial += 1
        else:
            report.by_energy += 1

        js = _js(_field(row, columns['jcol']),
                 _field(row, columns['termcol']))
        if not js:
            report.bad_j.append((number, _field(row, columns['jcol'])))
            continue
        for one in js:
            deck.add(Level(j=one, ee=value, ec=value, observed=True,
                           vector=dict(vector), kset=kset, line=number,
                           intensity=True))
            report.levels += 1
    return deck, report


def _rows(path):
    with open(path, newline='') as handle:
        return list(csv.reader(handle, delimiter='\t'))


def _header(rows, spec):
    """The header line, and the row the levels start on.

    With `--skip` given it is taken as read.  Without it, the leading
    rows that hold no number at all are the header, which covers a file
    with a header line and one with none alike -- and a file whose
    columns are named rather than numbered needs the header found
    before its columns can be resolved.
    """
    if spec.skip is not None:
        first = max(0, int(spec.skip))
    else:
        first = 0
        while first < len(rows) and not _numeric(rows[first]):
            first += 1
        if first >= len(rows):
            first = 0
    first = min(first, len(rows))
    return (rows[first - 1] if first else []), first


def _numeric(row):
    """Whether any field of a row is a number: a data row, not a header."""
    for field in row:
        try:
            float(field.strip().replace(',', ''))
            return True
        except ValueError:
            continue
    return False


def _column(given, name, header, path):
    """One column argument, as a 0-based index."""
    if given is None:
        return None
    text = str(given).strip()
    if not text:
        return None
    try:
        number = int(text)
    except ValueError:
        wanted = text.lower()
        for i, field in enumerate(header):
            if field.strip().strip('"').lower() == wanted:
                return i
        raise AltSourceError(
            '--%s: %s names no column of %s%s'
            % (name.replace('_', '-'), text, path,
               '' if header else ', which has no header line to name one in'))
    if number < 1:
        raise AltSourceError('--%s: %s is not a column number'
                             % (name.replace('_', '-'), text))
    return number - 1


def _field(row, index):
    if index is None or index >= len(row):
        return ''
    return row[index].strip().strip('"').strip()


def _excluded(row, index):
    """Anything non-blank and non-zero in the exclude column excludes."""
    text = _field(row, index)
    if not text:
        return False
    try:
        return float(text) != 0.0
    except ValueError:
        return True


def _parity(text):
    value = text.strip().strip('"').lower()
    if value in ODD:
        return 'o'
    if value in EVEN:
        return 'e'
    if value.endswith('*'):
        return 'o'
    return None


def _parity_sets(spec, target):
    if spec.parity_sets:
        out = {}
        for kset, piece in enumerate(str(spec.parity_sets).split(','), 1):
            parity = _parity(piece)
            if parity is None:
                raise AltSourceError('--parity-sets: %r is not a parity'
                                     % piece)
            out[parity] = kset
        return out
    out = parity_sets(target)
    if not out:
        raise AltSourceError(
            'the target deck names no configuration, so which parity set is '
            'which cannot be read off it; give --parity-sets')
    return out


def _percent(text):
    try:
        return abs(float(text.replace('%', '')))
    except ValueError:
        return 0.0


def _vector(row, columns, spec, labels, use_conf, report):
    """{key: intensity per cent} for one row, and its leading coupling.

    The triads, in the order the author's files write them:

    * the first is `(perc, conf, term)`, but `conf1`/`term1` override
      the labels wherever either is given -- a sheet that keeps `conf`
      and `term` for bookkeeping puts the component's real labels
      there;
    * the second and third are `(perc2, conf2, term2)` and the same
      again, a blank configuration or term meaning the one in the
      `conf` or `term` column;
    * unless the alternate-coupling column says otherwise, in which
      case the second triad is the *first* component written again in
      another scheme, and the third is the second component.

    Of the labels offered for one component, the first that is LS is
    taken.  If none is, the row has no vector at all: an experimental
    designation in another scheme cannot be compared with an LS deck,
    and the leading component is the one that decides.
    """
    conf = _field(row, columns['confcol'])
    term = _field(row, columns['termcol'])
    triads = [(_percent(_field(row, columns['perccol'])) or 100.0,
               [(_field(row, columns['conf1col']) or conf,
                 _field(row, columns['term1col']) or term),
                (conf, term)])]
    for n in ('2', '3'):
        one = _field(row, columns['conf%scol' % n])
        other = _field(row, columns['term%scol' % n])
        percent = _field(row, columns['perc%scol' % n])
        if not (one or other or percent):
            continue
        triads.append((_percent(percent), [(one or conf, other or term)]))

    alternate = (columns['altcol'] is not None
                 and _field(row, columns['altcol']).lower()
                 == str(spec.altvalue).strip().lower())
    if alternate and len(triads) > 1:
        # The second triad is the first component written again in
        # another coupling, not a second component: its labels join the
        # first component's, and its percentage -- the same level,
        # differently resolved -- is the one to believe.
        percent, offered = triads[1]
        triads[0] = (percent or triads[0][0], triads[0][1] + offered)
        del triads[1]

    out, scheme = {}, 'none'
    for i, (percent, offered) in enumerate(triads):
        label = None
        for one in offered:
            if one[1] and coupling(one[1]) == 'ls':
                label = one
                break
        if not i:
            for one in offered:
                if one[1]:
                    scheme = coupling(one[1]) or 'none'
                    break
            if label is None:
                return {}, scheme
        if label is None or percent <= 0:
            continue
        configuration = ''
        if use_conf:
            configuration = resolve_configuration(label[0], labels)
            if not configuration:
                report.unknown_conf[label[0]] = \
                    report.unknown_conf.get(label[0], 0) + 1
                if not i:
                    return {}, scheme
                continue
        key = '%s,%s' % (configuration, term_key(label[1]))
        out[key] = out.get(key, 0.0) + percent
    return out, scheme


def _js(j, term):
    """The J values one row stands for.

    A row that names a J is one level.  A row with the J column empty
    is a **term**, and stands for every J that term has, |L-S| to L+S;
    that is how a level list not yet resolved into levels still carries
    into a deck.
    """
    text = j.strip().strip('"').strip()
    if text:
        if '/' in text:
            return [text]
        try:
            return [_j_text(float(text))]
        except ValueError:
            return []
    found = _LS.match(term.strip().strip('"'))
    if not found:
        return []
    spin = (int(found.group(1)) - 1) / 2.0
    orbital = _L.index(found.group(2).upper())
    out, value = [], abs(orbital - spin)
    while value <= orbital + spin + 1e-9:
        out.append(_j_text(value))
        value += 1.0
    return out
