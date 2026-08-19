"""How well each coupling scheme describes each configuration.

A level's composition is written in some coupling scheme, and no one
scheme suits every configuration of a spectrum.  Which one to use is
decided per configuration, and the number that decides it is the
**purity**: the mean weight of the leading eigenvector component over
the levels of that configuration.  A scheme that concentrates a level
into one basis state is a scheme that names it well.

Comparing that number between schemes is legitimate because of an
invariance: a change of coupling recouples angular momenta *within* a
configuration, so the recoupling matrix is block-diagonal in the
configuration and the weight a configuration carries in a level is the
same in every scheme.  What changes is only how that weight is spread
over the configuration's own basis states -- which is exactly what
purity measures.  The same invariance is why a level has one dominant
configuration independent of the scheme, and RCE computes it once, from
full-precision coefficients, before any output pass (`NJMS`,
`FOR/RCE20KD.F:1157-1204`).  That is the `conf0` column, and levels are
grouped on it here, never on `conf1`: `conf1` is the leading
*component's* configuration, which moves with the scheme, and it carries
intermediate terms besides, so one configuration would appear under
several strings.

This module holds three things `config_purity.py` and
`make_levels_list_prepub.py` both need: the levels.txt reader, the
levels-to-configuration map, and the reader and writer of
`coupling_choice.txt`.
"""

import datetime
import os
import re

from dataclasses import dataclass, field

from . import asd, in36


class PurityError(Exception):
    """levels.txt or coupling_choice.txt will not support the request."""


#: The columns levels.txt always has, in order, before the components.
FIXED = ['cpl', 'par', 'Ee', 'e_c', 'Ec', 'Ee-Ec', 'J', 'Lande_g']

CHOICE = 'coupling_choice.txt'

#: Width of one scheme's field in coupling_choice.txt: the star slot,
#: three digits, and blanks enough to keep the columns apart by eye.
FIELD = 10


# --------------------------------------------------------------------------
# levels.txt
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Component:
    """One eigenvector component of a level, as levels.txt writes it."""

    percent: int
    configuration: str              # with its intermediate terms
    term: str

    @property
    def bare(self):
        """The configuration without its intermediate terms.

        `4f.5p6.5d2.(3F)` -> `4f.5p6.5d2`.  A parenthesised group is
        always an intermediate term of the subshell in front of it
        (`toolbox/asd.py`), so removing every one of them leaves the
        configuration and nothing else.
        """
        return canonical_configuration(
            strip_intermediates(self.configuration))


@dataclass
class Row:
    """One line of levels.txt."""

    coupling: str
    parity: str                     # 'e' or 'o'
    e_obs: float
    unfitted: bool                  # the e_c asterisk: Ee holds Ecalc
    e_calc: float
    j: float
    g: object = None                # float, or None where RCG printed none
    conf0: str = ''                 # empty unless --dominant-config was used
    components: list = field(default_factory=list)
    line: int = 0                   # 1-based line of the file

    @property
    def leading(self):
        return self.components[0] if self.components else None


@dataclass
class Table:
    """A whole levels.txt."""

    path: str
    header: list
    couplings: list                 # the cpl values, in the file's own order
    rows: list

    @property
    def has_conf0(self):
        return 'conf0' in self.header

    def by_coupling(self):
        out = {}
        for row in self.rows:
            out.setdefault(row.coupling, []).append(row)
        return out


def strip_intermediates(text):
    """`4f.5p6.5d2.(3F)` -> `4f.5p6.5d2`; also drops the separator."""
    out, depth = [], 0
    for ch in text:
        if ch == '(':
            depth += 1
            # The '.' that introduced the group goes with it.
            while out and out[-1] == asd.DOT:
                out.pop()
        elif ch == ')':
            depth = max(0, depth - 1)
        elif depth == 0:
            out.append(ch)
    return ''.join(out).strip(asd.DOT)


_SUBSHELL = re.compile(r'^(\d+)([a-z])([-+]?)(\d*)$')


def canonical_configuration(label):
    """The configuration a component belongs to, however it was written.

    Two coupling schemes and two decks can write one configuration three
    ways, and all three have to land on the deck card the choice file
    keys on:

    * true jj splits every open l>0 subshell into its j = l-1/2 and
      j = l+1/2 halves, so `2p-2.2p+` is the deck's `2p3` -- the split
      is a change of basis inside the configuration, which is the same
      invariance the rest of this module rests on;
    * a block merged from a reordered deck prints its subshells in that
      deck's summation order, so `2p2.1s` is `1s.2p2` -- `designate.py`
      keeps summation order wherever an intermediate term makes it
      visible, and that is right for a basis state and wrong for the
      name of a configuration.

    So the halves are added back together and the subshells put in
    (n, l) order, which is the order `deck_configurations` builds its
    labels in.  A label with a token this does not recognise is returned
    exactly as it came, rather than half-canonicalised.
    """
    parsed = []
    for token in label.split(asd.DOT):
        found = _SUBSHELL.match(token)
        if found is None:
            return label
        n, l, _, count = found.groups()
        parsed.append((int(n), asd.ORBITALS.index(l), n + l,
                       int(count) if count else 1))
    if not parsed:
        return label

    order, occupation = [], {}
    for n, l, key, count in sorted(parsed, key=lambda one: (one[0], one[1])):
        if key not in occupation:
            order.append(key)
            occupation[key] = 0
        occupation[key] += count
    return asd.DOT.join(
        key + ('' if occupation[key] == 1 else '%d' % occupation[key])
        for key in order)


def read_levels_txt(path='levels.txt'):
    """Read levels.txt into a `Table`.

    The header line is the specification: the fixed columns are found by
    name and the components by the `%n`/`confn`/`termn` triples, so a
    file written with `--dominant-config`, with more or fewer
    components, or with columns added later, all read the same way.
    """
    with open(path) as stream:
        lines = stream.read().splitlines()
    if not lines:
        raise PurityError('%s is empty' % path)
    header = lines[0].split('\t')
    at = {name: i for i, name in enumerate(header)}
    for name in FIXED:
        if name not in at:
            raise PurityError('%s has no %r column; it is not a levels.txt'
                              % (path, name))
    triples = []
    n = 1
    while ('%%%d' % n) in at:
        triples.append((at['%%%d' % n], at['conf%d' % n], at['term%d' % n]))
        n += 1

    rows, couplings = [], []
    for number, line in enumerate(lines[1:], 2):
        if not line.strip():
            continue
        cell = line.split('\t')

        def get(name, default=''):
            i = at.get(name)
            return cell[i].strip() if i is not None and i < len(cell) \
                else default

        row = Row(coupling=get('cpl'), parity=get('par'),
                  e_obs=_float(get('Ee'), path, number, 'Ee'),
                  unfitted=get('e_c') == '*',
                  e_calc=_float(get('Ec'), path, number, 'Ec'),
                  j=_float(get('J'), path, number, 'J'),
                  g=(float(get('Lande_g')) if get('Lande_g') else None),
                  conf0=get('conf0'), line=number)
        for percent, conf, term in triples:
            if percent >= len(cell) or not cell[percent].strip():
                break
            row.components.append(Component(
                percent=int(cell[percent]),
                configuration=cell[conf].strip() if conf < len(cell) else '',
                term=cell[term].strip() if term < len(cell) else ''))
        if row.coupling not in couplings:
            couplings.append(row.coupling)
        rows.append(row)
    if not rows:
        raise PurityError('%s carries no levels' % path)
    return Table(path=path, header=header, couplings=couplings, rows=rows)


def _float(text, path, number, name):
    try:
        return float(text)
    except ValueError:
        raise PurityError('%s line %d: %r is not a %s'
                          % (path, number, text, name))


# --------------------------------------------------------------------------
# the configurations of the deck
# --------------------------------------------------------------------------

@dataclass
class ConfigRow:
    """One configuration of the deck, as a row of the choice file."""

    index: int                      # 1-based position in IN36
    parity: str                     # 'e' or 'o'
    label: str                      # the bare ASD label, as conf0 writes it
    deck_label: str                 # the ING11/IN36 label, e.g. 'f5d2'
    nlev: int = 0
    purity: dict = field(default_factory=dict)   # coupling -> int
    chosen: str = ''                # the starred coupling


def deck_configurations(folder='.', drop=()):
    """The deck's configurations in IN36 order, labelled as `conf0` is.

    `printout.py` writes `conf0` by taking the configuration's occupied
    subshells less the closed core the whole run shares, in (n, l)
    order, through `asd.configuration`.  The same label is rebuilt here
    from IN36 alone rather than passed along, and every `conf0` found in
    levels.txt is then required to be one of these -- which makes the
    two independent constructions check each other instead of one of
    them being trusted.
    """
    deck = in36.read_in36(folder=folder)
    core = {(s.n, s.l) for s in in36.common_closed_shells(deck.configurations)}
    out = []
    for card in deck:
        tokens = [asd.Subshell(n=s.n, l=s.l, occ=s.occ)
                  for s in card.shells
                  if s.occ and (s.n, s.l) not in core
                  and not ((s.n, s.l) in drop and not s.is_open)]
        tokens.sort(key=lambda t: (t.n, t.lval))
        out.append(ConfigRow(index=card.index,
                             parity='o' if card.parity else 'e',
                             label=asd.configuration(tokens),
                             deck_label=card.label))
    return out


# --------------------------------------------------------------------------
# the purity table
# --------------------------------------------------------------------------

def purity_table(table, rows):
    """Fill `nlev` and `purity` of each `ConfigRow` from a `Table`.

    Returns `(kept, empty, nameless)`: the rows that have levels, the
    rows that have none, and the levels whose `conf0` matches no
    configuration of the deck.  Nothing is dropped silently; the caller
    decides what to say about each.
    """
    if not table.has_conf0:
        raise PurityError(
            '%s has no conf0 column, so its levels cannot be grouped by '
            'configuration.  Write it again with\n\n'
            '    printout levels.txt --couplings %s --dominant-config\n'
            % (table.path, ','.join(table.couplings)))
    by_label = {}
    for row in rows:
        if row.label in by_label:
            raise PurityError(
                'the configurations %r and %r of the deck have the same '
                'printed label %r, so their levels cannot be told apart'
                % (by_label[row.label].deck_label, row.deck_label, row.label))
        by_label[row.label] = row

    total = {}                      # (label, coupling) -> [sum, count]
    seen = {}                       # label -> set of (parity, Ec, J)
    nameless = []
    for level in table.rows:
        found = by_label.get(level.conf0)
        if found is None:
            nameless.append(level)
            continue
        seen.setdefault(level.conf0, set()).add(
            (level.parity, level.e_calc, level.j))
        if level.leading is None:
            # A block taken from a second folder can lose a level's
            # components; it has no purity to contribute, but it is
            # still one of the configuration's levels.
            continue
        cell = total.setdefault((level.conf0, level.coupling), [0, 0])
        cell[0] += level.leading.percent
        cell[1] += 1

    for row in rows:
        row.nlev = len(seen.get(row.label, ()))
        for coupling in table.couplings:
            cell = total.get((row.label, coupling))
            if cell and cell[1]:
                row.purity[coupling] = int(cell[0] / float(cell[1]) + 0.5)

    kept = [r for r in rows if r.nlev]
    empty = [r for r in rows if not r.nlev]
    return kept, empty, nameless


def suggest(row, couplings):
    """The coupling with the highest purity, ties going to file order."""
    best, chosen = None, ''
    for coupling in couplings:
        value = row.purity.get(coupling)
        if value is not None and (best is None or value > best):
            best, chosen = value, coupling
    return chosen


# --------------------------------------------------------------------------
# coupling_choice.txt
# --------------------------------------------------------------------------

_NOTE = """\
Purity of a configuration in a coupling scheme is the UNWEIGHTED mean,
over the levels whose dominant configuration (conf0) is that one, of the
leading eigenvector component's percentage.  Unweighted: no 1/(2J+1).

Mark the scheme each configuration is to be written in by putting a '*'
anywhere in that scheme's column -- exactly one star on every line.  A
star is already on the highest purity, and that is a suggestion, not the
answer: highest purity is the main criterion but not the only one.
Previous literature on the spectrum, and the schemes used for the
neighboring isoelectronic ions, routinely outweigh a few percent gained
by a rarely used scheme.  Nlev is beside the purities for the same
reason -- 88% over three levels is not evidence.

Everything left of the first scheme column is written by config_purity.py
and read back by position; edit the stars and nothing else.
"""


def write_choice(path, rows, couplings, folder='.', levels='levels.txt'):
    """Write coupling_choice.txt: one line per configuration."""
    width = max([len('configuration')] + [len(r.label) for r in rows])
    lines = []
    lines.append('# %s -- written by config_purity.py on %s'
                 % (os.path.basename(path),
                    datetime.date.today().isoformat()))
    lines.append('# folder: %s' % os.path.abspath(folder))
    lines.append('# levels: %s' % levels)
    lines.append('#')
    for line in _NOTE.splitlines():
        lines.append(('# ' + line).rstrip())
    lines.append('#')
    lines.append(_line('par', 'configuration', 'Nlev', width,
                       [c.ljust(FIELD) for c in couplings]))
    for row in rows:
        fields = []
        for coupling in couplings:
            value = row.purity.get(coupling)
            star = '*' if row.chosen == coupling else ' '
            fields.append((' ' * FIELD if value is None
                           else (star + '%3d' % value).ljust(FIELD)))
        lines.append(_line(row.parity, row.label, '%d' % row.nlev, width,
                           fields))
    with open(path, 'w') as out:
        out.write('\n'.join(lines) + '\n')


def _line(parity, label, nlev, width, fields):
    return ('%-4s%-*s%5s  %s'
            % (parity, width + 2, label, nlev, ''.join(fields))).rstrip()


def read_choice(path=CHOICE):
    """Read a marked coupling_choice.txt: {configuration label: coupling}.

    The header line fixes every column, so a scheme is selected by the
    field the star falls in and not by counting spaces.  A line with no
    star, or with two, stops the run and is named: a wrong scheme here
    silently rewrites a whole configuration's compositions, and guessing
    which of two stars was meant is not something to do on the author's
    behalf.
    """
    with open(path) as stream:
        lines = stream.read().splitlines()
    head = None
    for line in lines:
        if line.startswith('#') or not line.strip():
            continue
        head = line
        break
    if head is None:
        raise PurityError('%s carries no header line' % path)

    spans = _spans(head, path)
    couplings = [name for name, _, _ in spans[3:]]
    out = {}
    started = False
    for number, line in enumerate(lines, 1):
        if line.startswith('#') or not line.strip():
            continue
        if not started:
            started = True                      # the header itself
            continue
        label = line[spans[1][1]:spans[1][2]].strip()
        marked = [name for name, start, stop in spans[3:]
                  if '*' in line[start:stop]]
        if len(marked) != 1:
            raise PurityError(
                '%s line %d (%s): %s.  Every configuration takes exactly one '
                'coupling scheme, marked by one star.'
                % (path, number, label or 'no configuration',
                   'no scheme is starred' if not marked
                   else 'these schemes are all starred: %s'
                        % ', '.join(marked)))
        if label in out:
            raise PurityError('%s names the configuration %r twice'
                              % (path, label))
        out[label] = marked[0]
    if not out:
        raise PurityError('%s carries no configurations' % path)
    return out, couplings


def _spans(head, path):
    """[(name, start, stop)] from the header line, in column order."""
    names, at = [], 0
    while at < len(head):
        if head[at] == ' ':
            at += 1
            continue
        stop = at
        while stop < len(head) and head[stop] != ' ':
            stop += 1
        names.append((head[at:stop], at))
        at = stop
    if len(names) < 4 or [n for n, _ in names[:3]] != ['par', 'configuration',
                                                       'Nlev']:
        raise PurityError(
            '%s does not begin par / configuration / Nlev / <schemes>; its '
            'header line reads %r' % (path, head))
    out = []
    for i, (name, start) in enumerate(names):
        stop = names[i + 1][1] if i + 1 < len(names) else len(head) + FIELD
        out.append((name, start, stop))
    return out


def merge_choice(path, rows, couplings):
    """Carry the stars of an existing choice file onto fresh rows.

    Returns the labels that were starred but are no longer in the deck,
    so the caller can say what it dropped.  A configuration the file
    does not name keeps the suggestion it was given.
    """
    old, _ = read_choice(path)
    by_label = {r.label: r for r in rows}
    for label, coupling in old.items():
        row = by_label.get(label)
        if row is not None and coupling in couplings:
            row.chosen = coupling
    return sorted(label for label in old if label not in by_label)
