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
recovers the missing digits from the label alone.  IN36 is where they
are written out, so where there is one beside the deck a source
configuration is the card whose subshells it names, closed ones it
leaves out allowed to be left out (`ConfIndex`).  Without an IN36 the
weaker rule stands: the deck's label is a tail of the compacted full
form, and the longest such label wins.  Either way a configuration
nothing matches is reported rather than guessed at.

**The comparison is in intensity per cent, not amplitude.**  A source
sheet quotes what a published table quotes, so the target's amplitudes
are squared over 100 to meet it and the sign question does not arise.
`FLOOR` still cuts the target's components in amplitude space, so the
same components survive as in a deck-to-deck run.

**A component is matched to whatever depth it is named.**  Experimental
designations are chosen for whichever scheme makes the level look pure,
so one file can hold LS, JK and pair labels side by side, and only LS
can be matched against an LS deck.  But a label in another scheme still
writes its parent terms inside the configuration -- `4f2.(3F<2>).6p`
beside `(2,3/2)*` -- and the last of those is exactly what RCG has room
for in brackets in the eight characters RCE prints.  So a component is
keyed on the configuration, that parent and the final term separately,
with `explevs.ANY` for each one it does not give, and `explevs.project`
adds the target's components up over the fields it left open.  A row
whose leading component gives neither an LS parent nor an LS term is
matched on (parity, J, energy) instead, and counted in the report so
that the share of the file that got the weaker treatment is never a
surprise.
"""

import csv
import os
import re

from . import files, in36, purity
from .explevs import Deck, Level, _j_text, match_key
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

#: The `+x` an experimental compilation writes on an energy whose
#: system is tied to the ground level only loosely -- and the `-x` form,
#: for a compilation that writes it that way.  One letter, at the end of
#: the field, is the whole of it; a level whose energy is `12345.678+x`
#: is a member of system `x` and its energy is 12345.678.
#:
#: Until this was recognised such a value simply failed `float()` and the
#: row was dropped without a word, which is the defect `02_python_audit.md`
#: Sec. 4 found.  Both halves are now read: the number, and the letter.
_SUFFIX = re.compile(r'^(.*?)\s*[+-]\s*([A-Za-z])$')

#: A subshell of a configuration spelled out in full, `4f3` or `6p`.
_WRITTEN = re.compile(r'^(\d+)([%s])(\d*)$' % _ORBITALS)

#: What a parent term is written inside, and the J a J1J2 label adds.
_PARENS = re.compile(r'\(([^)]*)\)')
_JOFTERM = re.compile(r'<[^>]*>')

#: How RCG writes the eight-character label RCE prints: the last term
#: the chain writes, in brackets, and the final term after it.
_LABEL = re.compile(r'^\(([^)]*)\)\s*(.*)$')


class AltSourceError(Exception):
    """A source level list the options do not describe."""


class Spec:
    """Which column of the source file is which.

    Every field is either `None` (the column is not there), a 1-based
    column number, or a column name to be looked up in the header.  The
    reader resolves them once it has read the file, because only then
    does it know whether there is a header at all.
    """

    NAMES = ('ecol', 'jcol', 'pcol', 'excludecol', 'altcol', 'basecol',
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
        self.by_parent = 0          # components named to the parent only
        self.in36 = ''              # the IN36 the configurations came from
        self.couplings = {}         # {scheme: rows whose leading term is it}
        self.unknown_conf = {}      # {source configuration: rows}
        self.bad_parity = []
        self.bad_j = []
        self.bad_energy = []        # rows whose energy does not parse
        self.groups = {}            # {suffix letter: [(line, is base)]}
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
        if self.in36:
            out.write('   configurations read from %s\n' % self.in36)
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
        if self.by_parent:
            out.write('   %d component%s named only by the parent term '
                      'written in the configuration; those are matched '
                      'against the parent of the deck\'s own label.\n'
                      % (self.by_parent, '' if self.by_parent == 1 else 's'))
        if self.partial:
            out.write('   Warning: %d of them name one component only, which '
                      'is thin evidence; check those matches.\n'
                      % self.partial)
        for label, count in sorted(self.unknown_conf.items()):
            out.write('   Warning: configuration %s (%d level%s) is in no '
                      'component of the target deck.\n'
                      % (label, count, '' if count == 1 else 's'))
        if self.groups:
            out.write('   relative-energy systems: %s\n'
                      % ', '.join('+%s (%d level%s)'
                                  % (letter, len(rows),
                                     '' if len(rows) == 1 else 's')
                                  for letter, rows
                                  in sorted(self.groups.items())))
        for number, text in self.bad_energy[:5]:
            out.write('   Warning: line %d: %r is not an energy, so the '
                      'level was not read.\n' % (number, text))
        if len(self.bad_energy) > 5:
            out.write('   Warning: %d further rows have an energy that does '
                      'not parse.\n' % (len(self.bad_energy) - 5))
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
    """The LS term a component is matched on: `2P*` and `2P2` -> `2P`.

    A J1J2 parent carries its J as well (`3H<4>`), and that J is no
    part of the LS term RCG prints in its label, so it goes too.
    """
    text = _JOFTERM.sub('', term.strip().strip('"'))
    text = text.split()[-1] if text.split() else ''
    text = text.rstrip('*').upper()
    found = _REPEAT.match(text)
    return found.group(1) if found else text


def parent_term(configuration):
    """The LS parent a configuration label writes, or `''`.

    Spectroscopic notation writes inside the configuration exactly the
    terms of the coupling chain that the final term does not give back
    (`designate.py`), and the **last** of them is the one RCG has room
    for in the eight characters RCE prints -- `4f2.(3F<2>).6p` against
    `f26p   (3F) 4D`, `4f.(2F*).5d2.(3F)` against `f5d2   (3F) 2P`.

    So it is worth having even where the final term is not: a level
    designated `(2,3/2)*` is in pair coupling and cannot be compared
    with an LS deck at all, but the `3F` written beside it can.  Only
    an LS parent is taken; a J1J2 or pair parent is not one.
    """
    found = _PARENS.findall(configuration.strip().strip('"'))
    if not found:
        return ''
    text = _JOFTERM.sub('', found[-1]).strip()
    return term_key(text) if coupling(text) == 'ls' else ''


def split_label(label):
    """`(3F) 4D3` -> `('3F', '4D')`: what RCE prints for one component."""
    text = ' '.join(str(label).split())
    found = _LABEL.match(text)
    if found is None:
        return '', term_key(text)
    return term_key(found.group(1)), term_key(found.group(2))


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


def written_shells(configuration):
    """`{(n, l): occupation}` for a configuration spelled out in full.

    `None` where the label is not spelled out -- a token with no
    principal quantum number, or none this recognises at all -- which
    is what sends the caller back to matching on the written form.
    """
    text = purity.strip_intermediates(configuration.strip().strip('"'))
    if not text:
        return None
    out = {}
    for token in text.split('.'):
        found = _WRITTEN.match(token.strip().lower())
        if found is None:
            return None
        n, letter, count = found.groups()
        key = (int(n), _ORBITALS.index(letter))
        if key in out:
            return None
        out[key] = float(count) if count else 1.0
    return out or None


class ConfIndex:
    """The deck's configurations as IN36 spells them out.

    The deck itself cannot say what its labels mean.  RCE writes a
    configuration in six columns and abbreviates from the left to get
    there -- `4f.5d2` is `f5d2`, `5d.6s.6p` is `d6s6p` -- and nothing
    recovers the missing digits from the label alone.  The old reader
    therefore matched a written-out configuration against the deck
    label its compact form *ends with*, which works for the labels that
    really are abbreviations and quietly fails for the ones that are
    not: a card may be nicknamed anything, and Pr III has `f4p5` for
    4f4.5p5 among its twenty-one.

    IN36 is where the configurations are written out, and `in36`
    already reads it the way RCN does, so the mapping is made there
    instead: a source configuration is the card whose subshells it
    names, subshell for subshell, with the closed ones it leaves out
    allowed to be left out.  A source that fits two cards equally well
    -- `4f.6s2` and `4f.7s2` are told apart only by a closed shell --
    fits neither, and is reported rather than guessed at.
    """

    def __init__(self, deck):
        self.path = getattr(deck, 'path', '')
        self.cards = [card for card in deck.configurations if card.shells]
        self.parity = {}
        for card in self.cards:
            self.parity.setdefault(card.label, card.parity)

    def __bool__(self):
        return bool(self.cards)

    def resolve(self, configuration):
        """The label of the card a written-out configuration names."""
        wanted = written_shells(configuration)
        if wanted is None:
            return ''
        best, spare = '', None
        for card in self.cards:
            extra = self._extra(card, wanted)
            if extra is None:
                continue
            if spare is None or extra < spare:
                best, spare = card.label, extra
            elif extra == spare and card.label != best:
                best = ''                      # two cards fit alike
        return best

    @staticmethod
    def _extra(card, wanted):
        """How many closed subshells the card has and the source left out.

        `None` where the card is not the configuration at all: an
        occupation that differs, a subshell the source names and the
        card has not, or one the card has open and the source is silent
        about.
        """
        held = {(s.n, s.l): s.occ for s in card.shells if s.occ}
        for key, occ in wanted.items():
            if held.get(key) != occ:
                return None
        extra = 0
        for key, occ in held.items():
            if key in wanted:
                continue
            if occ != 4 * key[1] + 2:
                return None
            extra += 1
        return extra


def read_configurations(target):
    """The `ConfIndex` of the folder the target deck is in, or `None`.

    IN36 is RCN's input, so a folder holding only the files RCG and RCE
    rewrite has none of its own; the directory above is tried next,
    which is the sub-folder arrangement `toolbox.case` reads the same
    way.  Without one the reader falls back on the written form.
    """
    here = os.path.dirname(os.path.abspath(getattr(target, 'path', '') or '.'))
    for folder in (here, os.path.dirname(here)):
        if not folder or not os.path.isdir(folder):
            continue
        try:
            path = files.resolve(files.IN36, folder, required=False)
            if path is None:
                continue
            deck = in36.read_in36(path)
        except (in36.In36Error, OSError, UnicodeDecodeError):
            # A deck RCN itself would reject, or two files differing
            # only in case: the fallback still has something to say, so
            # this is not the place to stop the run.
            return None
        index = ConfIndex(deck)
        return index if index else None
    return None


def resolve_configuration(configuration, labels, index=None):
    """The deck label a written-out configuration stands for.

    With IN36 beside the deck (`index`) the answer is read off the
    cards; the label it gives is still required to be one the deck
    itself uses, so a card RCE never got is not put in a key that can
    match nothing.  Without IN36 the old rule stands: the deck's label
    is a tail of the compacted full form, and the longest one wins.
    """
    if index is not None:
        label = index.resolve(configuration)
        if label and (not labels or label in labels):
            return label
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


def parity_sets(target, index=None):
    """{'e': parity set, 'o': parity set}, from the target's own labels.

    The old reader let the first row of the source settle which parity
    set was which, because a LOPT list does not say.  The target deck
    does say -- every component names a configuration, and a
    configuration has a parity -- so it is read off there instead, and
    the source is only asked which parity each of *its* levels is.

    With IN36 beside the deck (`index`) a label's parity is the card's
    and there is nothing to weigh up.  Without it the parity has to be
    read out of the abbreviation itself, which a nickname that is not
    an abbreviation gets wrong, so a majority of the set decides.
    """
    out = {}
    for kset in sorted(target.sets):
        votes = {0: 0, 1: 0}
        for j in target.sets[kset]:
            for level in target.group(kset, j):
                label = level.leading[1].partition(',')[0]
                if not label:
                    continue
                if index is not None and label in index.parity:
                    votes[index.parity[label]] += 1
                elif index is None:
                    votes[parity_of_label(label)] += 1
        if not (votes[0] or votes[1]):
            continue
        out['e' if votes[0] >= votes[1] else 'o'] = kset
    return out


# --------------------------------------------------------------------------
# the target, in the same key space
# --------------------------------------------------------------------------

def normalize_target(target, use_conf=True):
    """Give every target level the vector the source is compared to.

    In intensity per cent, and keyed on all three of the fields a
    source row may name: the configuration, the parent term RCE writes
    in brackets, and the final term.  Nothing is added together here --
    `explevs.project` does that, for each source level, over exactly
    the fields that level left open, so a source that quotes `4f3 2P`
    is not distinguishing the parents and one that quotes
    `4f2.(3F).6p` is distinguishing nothing else.
    """
    for level in target.levels():
        out = {}
        for key, percent in level.vector.items():
            configuration, _, label = key.partition(',')
            parent, term = split_label(label)
            mine = match_key(configuration if use_conf else '', parent, term)
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

    index = read_configurations(target)
    if index is not None:
        report.in36 = index.path
    wants_conf = any(columns[name] is not None
                     for name in ('confcol', 'conf1col', 'conf2col',
                                  'conf3col'))
    labels = deck_labels(target) if wants_conf else set()
    use_conf = bool(labels)
    report.use_conf = use_conf
    report.parities = _parity_sets(spec, target, index)

    deck = Deck(path)
    for number, row in enumerate(rows[first:], first + 1):
        if not any(field.strip() for field in row):
            continue
        report.rows += 1
        if _excluded(row, columns['excludecol']):
            report.excluded += 1
            continue
        text = _field(row, columns['ecol'])
        value, suffix = energy_and_suffix(text)
        if value is None:
            # Reported, never dropped in silence: until this was here a
            # `12345.678+x` -- the very notation this feature is about
            # -- made the row vanish and the counts at the end of the
            # run simply failed to add up.
            report.bad_energy.append((number, text))
            continue
        value *= spec.scale
        parity = _parity(_field(row, columns['pcol']))
        kset = report.parities.get(parity) if parity else None
        if kset is None:
            report.bad_parity.append((number, _field(row, columns['pcol'])))
            continue

        vector, scheme = _vector(row, columns, spec, labels, use_conf,
                                 report, index)
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
        base = bool(_field(row, columns['basecol'])) if suffix else False
        status = (suffix.upper() if base else suffix) if suffix else ' '
        for one in js:
            deck.add(Level(j=one, ee=value, ec=value, observed=True,
                           status=status, vector=dict(vector), kset=kset,
                           line=number, intensity=True))
            report.levels += 1
            if suffix:
                report.groups.setdefault(suffix, []).append((number, base))
    _settle_bases(deck, report)
    return deck, report


def _settle_bases(deck, report):
    """Check that each suffix group names exactly one base level.

    Which level of a `+x` system carries the absolute tie to the ground
    level is a decision about the *measurements*, and this module has no
    way to make it: the theory in the target deck certainly must not,
    and the lowest level of the group is a convention, not a fact.  So
    it is asked for, with `--basecol`, and refused where it is missing
    or given twice -- the one exception being a group of a single level,
    where there is nothing to choose between.
    """
    for letter in sorted(report.groups):
        rows = report.groups[letter]
        bases = [number for number, base in rows if base]
        if len(bases) > 1:
            raise AltSourceError(
                'the +%s system has %d levels marked as its base, on lines '
                '%s; exactly one level of a group carries the tie to the '
                'ground level' % (letter, len(bases),
                                  ', '.join('%d' % n for n in bases)))
        if bases:
            continue
        if len(rows) > 1:
            raise AltSourceError(
                'the +%s system has %d levels and none of them is marked as '
                'its base.  Which one carries the measured tie to the ground '
                'level is a decision about the data, not one this can make: '
                'mark it with a column and name that column with --basecol.'
                % (letter, len(rows)))
        # A group of one: it is its own base, and nothing was chosen.
        for level in deck.levels():
            if level.status == letter:
                level.status = letter.upper()


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
    """Whether any field of a row is a number: a data row, not a header.

    A suffixed energy counts as one.  A file whose every energy carries
    a `+x` would otherwise look like header all the way down, which is
    the second place `02_python_audit.md` Sec. 4 found the same blind
    `float()`.
    """
    for field in row:
        if energy_and_suffix(field)[0] is not None:
            return True
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


def energy_and_suffix(text):
    """`'12345.678+x'` -> `(12345.678, 'x')`; `(None, '')` for a non-number.

    The suffix is the compilation's own notation for a level system
    known precisely inside itself and tied to the ground level only
    loosely; it is what becomes a group letter in the deck.  Anything
    that is not a number even after the suffix comes off gives `None`,
    which the caller reports -- it must never be dropped in silence.
    """
    body = text.strip().strip('"').strip()
    suffix = ''
    found = _SUFFIX.match(body)
    if found:
        body, suffix = found.group(1), found.group(2).lower()
    try:
        return float(body.replace(',', '')), suffix
    except ValueError:
        return None, ''


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


def _parity_sets(spec, target, index=None):
    if spec.parity_sets:
        out = {}
        for kset, piece in enumerate(str(spec.parity_sets).split(','), 1):
            parity = _parity(piece)
            if parity is None:
                raise AltSourceError('--parity-sets: %r is not a parity'
                                     % piece)
            out[parity] = kset
        return out
    out = parity_sets(target, index)
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


def _vector(row, columns, spec, labels, use_conf, report, index=None):
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
        # Of the labels offered for one component -- the same component
        # written in one scheme or another -- the first LS final term
        # is taken, and the first LS parent written inside any of the
        # configurations beside it.  Either alone is worth having: the
        # parent is what an LS deck's label carries in brackets, so a
        # component designated `(2,3/2)*` in pair coupling still says
        # something an LS deck can be asked about.
        term, written, parent = '', '', ''
        for one in offered:
            if one[1] and not term and coupling(one[1]) == 'ls':
                term = term_key(one[1])
            if one[0] and not written:
                written = one[0]
            parent = parent or parent_term(one[0])
        if not i:
            for one in offered:
                if one[1]:
                    scheme = coupling(one[1]) or 'none'
                    break
            if not (term or parent):
                return {}, scheme
        if not (term or parent) or percent <= 0:
            continue
        configuration = ''
        if use_conf:
            configuration = resolve_configuration(written, labels, index)
            if not configuration:
                report.unknown_conf[written] = \
                    report.unknown_conf.get(written, 0) + 1
                if not i:
                    return {}, scheme
                continue
        if not term:
            report.by_parent += 1
        key = match_key(configuration, parent, term)
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
