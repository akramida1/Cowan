"""Reading OUTG11: the configuration table, the J blocks, the vectors.

OUTG11 is the big file of a run - 92 MB for Pr III - so this module
streams it.  read_outg11() returns a Run whose blocks are read one at a
time as they are asked for; nothing but the configuration table and the
block being looked at is ever resident.

What the reader recovers, and the FORMAT that specifies each part:

    the configuration table    rcg11kd.f:3284 (heading), :3432 (a row)
    a J block heading          rcg11kd.f:11279-11290
    CONFIG. NO.                rcg11kd.f:11294
    G-VALUES                   rcg11kd.f:1259
    EIGENVECTORS (xxx)         rcg11kd.f:1273, :1785
    the column headings        rcg11kd.f:1318, :1320
    a row, LS                  rcg11kd.f:1350
    a row, JJ and TRUEjj       rcg11kd.f:1353
    a row, couplings 3 to 7    rcg11kd.f:1802, :1807, :1813, :1817
    PURITY                     rcg11kd.f:1379, :1842

Every one of those row FORMATs ends the same way: the state number
finishes in column 28 and the amplitudes follow as F9.5 from column 29.
Only the 24 columns in between differ by coupling, and this module keeps
them as written and leaves their reading to the coupling that owns them.

A section is printed eleven columns at a time.  The column headings are
written once, for the coupling the run was set up in (KCPL); the other
sections carry rows only, in the same column order, so the groups are
found by the blank line that separates them.
"""

import re
from dataclasses import dataclass, field as _field

from . import files

# COUPL, rcg11kd.f:1079-1080.  Index 8 is the true jj coupling added in
# 2026; 9 to 11 are unused and 12 is SL.
COUPLINGS = {1: 'LS', 2: 'JJ', 3: 'JJJK', 4: 'LSLK', 5: 'LSJK',
             6: 'LSJLKS', 7: 'LSJLSJ', 8: 'TRUEjj', 12: 'SL'}

# Where the amplitudes start, and how wide they are, in every row FORMAT.
_DATA = 28
_W = 9

_EIGENVECTORS = '   EIGENVECTORS   ('
_EIGENVALUES = '  EIGENVALUES '


class Outg11Error(Exception):
    """OUTG11 did not read as RCG writes it."""


@dataclass(frozen=True)
class ConfigCard:
    """One row of the configuration table at the head of the run."""
    kset: int            # configuration set, 1 or 2
    number: int          # its number within that set
    shells: tuple        # ((letter, occupation), ... eight slots, in the
                         # order of the ING11 card, which is the order
                         # RCG sums them in
    nelec: int
    parity: int          # as printed: -1 odd, +1 even
    terms: str           # 'ALL', or the term list where a cutoff applies
    name: str            # the 18-character 3A6 spectrum and label
    eav: float
    kinetic: float

    @property
    def label(self):
        """The 12 characters a row of an eigenvector section carries."""
        return self.name[6:18].rstrip()

    @property
    def spectrum(self):
        return self.name[:6].strip()

    @property
    def open_shells(self):
        return tuple((l, n) for l, n in self.shells
                     if n and n != 2 * (2 * _L.index(l) + 1))


_L = 'SPDFGHIKLMNOQRTUVWXYZABC'      # LSYMB, rcg11kd.f:1081-1083


@dataclass(frozen=True)
class Row:
    """One basis state of one eigenvector section."""
    cfg: int             # configuration number within the set
    desc: str            # the columns between the ':' or '=' and the
                         # state number, as written: the configuration
                         # label, then whatever the coupling puts after
                         # it.  Reading that tail belongs to the coupling
                         # that owns the FORMAT, not here.
    number: int          # its place in the basis, 1-based
    amplitudes: tuple

    @property
    def label(self):
        """The configuration label, 2A6 in most rows and A6A5 in some."""
        return self.desc[:12].rstrip()

    @property
    def bracketed(self):
        """The A7 or A8 state label of an LS, JJ or TRUEjj row."""
        at = self.desc.find('(')
        return self.desc[at + 1:].rstrip() if at >= 0 else ''


@dataclass
class Section:
    """One coupling of one J block."""
    coupling: str
    rows: list = _field(default_factory=list)
    column_labels: list = _field(default_factory=list)
    # the eight-character configuration label of each eigenvector, from
    # the heading; empty for a section that carries no heading
    column_names: list = _field(default_factory=list)
    # its designation, LDES, likewise
    purity: list = _field(default_factory=list)
    average_purity: float = None

    @property
    def nstates(self):
        return len(self.rows)


@dataclass
class JBlock:
    """One J block of one configuration set."""
    j: float
    kset: int = 0        # filled in from the configuration table
    eigenvalues: list = _field(default_factory=list)
    config_numbers: list = _field(default_factory=list)
    g_values: list = _field(default_factory=list)
    sections: dict = _field(default_factory=dict)

    @property
    def nlevels(self):
        return len(self.eigenvalues)

    def section(self, coupling):
        try:
            return self.sections[coupling]
        except KeyError:
            raise Outg11Error(
                'OUTG11 has no %s section for J=%.1f of configuration set %d.'
                '  That coupling is asked for in columns 31-39 of the ING11'
                ' main card, %s in column %d.'
                % (coupling, self.j, self.kset, coupling,
                   30 + _column_of(coupling)))


def _column_of(coupling):
    for index, name in COUPLINGS.items():
        if name == coupling:
            return index
    raise Outg11Error('%r is not one of RCG\'s couplings' % coupling)


class Run:
    """The configuration table of an OUTG11, and a way through its blocks."""

    def __init__(self, path, configs):
        self.path = path
        self.configs = configs          # {(kset, number): ConfigCard}

    def cards(self, kset):
        return [c for (k, _), c in sorted(self.configs.items()) if k == kset]

    def card(self, kset, number):
        try:
            return self.configs[(kset, number)]
        except KeyError:
            raise Outg11Error('OUTG11 has no configuration %d in set %d'
                              % (number, kset))

    def blocks(self, couplings=None):
        """Yield the J blocks in the order they are printed.

        `couplings` limits what is kept: with LS alone the eigenvector
        rows of the other sections are dropped as they are read, which
        is what keeps a Pr III pass small.  A requested coupling that
        never appears is the caller's error to report; JBlock.section
        names it.
        """
        for block in _blocks(self.path, self.configs, couplings):
            yield block


# --------------------------------------------------------------------------
# fixed-column helpers
# --------------------------------------------------------------------------

def _fields(line, start, width=_W):
    """Split a fixed-width numeric line, which may have no gaps in it.

    Eigenvalues are F9.2 and run together above 100000 cm-1, so they
    cannot be read by splitting on blanks; the column grid is the only
    thing that separates them.
    """
    out = []
    at = start
    while at < len(line):
        text = line[at:at + width].strip()
        if not text:
            break
        out.append(text)
        at += width
    return out


def _floats(line, start, width=_W):
    try:
        return [float(t) for t in _fields(line, start, width)]
    except ValueError:
        raise Outg11Error('a number will not read in %r' % line.rstrip())


def _ints(line, start, width=_W):
    try:
        return [int(t) for t in _fields(line, start, width)]
    except ValueError:
        raise Outg11Error('a number will not read in %r' % line.rstrip())


# Where the state number sits in a row, by coupling.  Every row FORMAT
# fills exactly 28 columns before the amplitudes, but divides them
# differently, and the number is the one field this module has to read
# for itself:
#   LS          I3,':',2A6,'(',A7,I3,1X          rcg11kd.f:1350
#   JJ, TRUEjj  I3,':',2A6,'(',A8,I3             rcg11kd.f:1353
#   JJJK, LSJK  I3,'=',A6,A5,F4.1,'K',F4.1,I3,1X rcg11kd.f:1802
#   LSLK        I3,'=',2A6,2X,A1,F4.1,I4,1X      rcg11kd.f:1807
#   LSJLKS      I3,'=',A6,A5,F4.1,A1,F3.1,I2,I3  rcg11kd.f:1813
#   LSJLSJ      I3,'=',A6,A5,F4.1,I2,A1,F3.1,I3  rcg11kd.f:1817
_NUMBER_AT = {'LS': (24, 27), 'JJ': (25, 28), 'TRUEjj': (25, 28),
              'JJJK': (24, 27), 'LSJK': (24, 27), 'LSLK': (23, 27),
              'LSJLKS': (25, 28), 'LSJLSJ': (25, 28), 'SL': (24, 27)}


def _all_floats(line):
    """The numbers of a line that holds nothing else, or None."""
    text = line.split()
    if not text:
        return None
    try:
        return [float(t) for t in text]
    except ValueError:
        return None


def _is_row(line, coupling):
    """True for a line of an eigenvector section of `coupling`.

    Every row FORMAT opens I3 then a single character - ':' for LS, JJ
    and TRUEjj, '=' for couplings 3 to 7.  That much is not enough: the
    Hamiltonian and dipole matrices are printed the same way, and the
    print of one block's matrix stands between the eigenvector sections
    of the block before it and the eigenvalues of its own.  What tells
    them apart is the bracket - the eigenvector FORMATs open the state
    label with 1H( in column 17 (rcg11kd.f:1350, :1353) and the matrix
    FORMATs with 2H ( in column 18 (:1466, :1485, :1488).
    """
    if len(line) <= _DATA or not line[:3].strip().isdigit():
        return False
    if coupling in ('LS', 'JJ', 'TRUEjj', 'SL'):
        return line[3] == ':' and line[16] == '('
    return line[3] == '='


def _row(line, coupling):
    at, to = _NUMBER_AT.get(coupling, (25, 28))
    try:
        return Row(cfg=int(line[:3]), desc=line[4:at],
                   number=int(line[at:to]),
                   amplitudes=tuple(_floats(line, _DATA)))
    except ValueError:
        raise Outg11Error('a row of the %s section will not read: %r'
                          % (coupling, line.rstrip()))


# --------------------------------------------------------------------------
# the configuration table
# --------------------------------------------------------------------------

def _config_row(line):
    """One row of the configuration table.

    Two FORMATs write it, and which one is used depends on the
    configuration:

      23  I4,I4,7X,8(A1,I2,3X),I5,I5,8X,'ALL',8X,3A6,F13.3,F9.4
          rcg11kd.f:3432 - the ordinary row, all terms of every
          subshell kept, with the 3A6 name and the average energy.
      21  I4,I4,7X,8(A1,I2,3X),I5,I5,I5,3X,13(I1,A1,A2)/81X,...
          rcg11kd.f:3420 - where the terms of a subshell are cut back,
          the subshell number and the surviving terms take the place of
          the name and the energy, which are then NOT PRINTED AT ALL.

    So a configuration under a term cutoff comes back with an empty name
    and a zero energy.  That is what the file holds; ing11.py reads the
    names from the deck, which always has them.

    Either way the subshell slots start in column 16 and run six columns
    each, in the order of the ING11 card - which is the order RCG sums
    them in, not a canonical one.
    """
    wide = line.ljust(133)
    shells = []
    for i in range(8):
        at = 15 + 6 * i
        letter = wide[at:at + 1]
        if not letter.strip():
            shells.append((' ', 0))
            continue
        shells.append((letter, int(wide[at + 1:at + 3])))
    common = dict(kset=int(wide[:4]), number=int(wide[4:8]),
                  shells=tuple(shells),
                  nelec=int(wide[63:68]), parity=int(wide[68:73]))
    if wide[81:84] == 'ALL':
        return ConfigCard(terms='ALL', name=wide[92:110],
                          eav=float(wide[110:123] or 0),
                          kinetic=float(wide[123:132] or 0), **common)
    return ConfigCard(terms=wide[73:78].strip() + ': ' + wide[81:].strip(),
                      name='', eav=0.0, kinetic=0.0, **common)


def read_outg11(path=files.OUTG11, folder=None):
    """Read the configuration table of an OUTG11 and return a Run."""
    if folder is not None:
        path = files.resolve(path, folder)
    configs = {}
    with open(path, 'r', errors='replace') as fh:
        for line in fh:
            if line.startswith('0  K   J') or _EIGENVALUES in line[:16]:
                if configs:
                    break               # the table is over
                continue
            if _looks_like_config_row(line):
                card = _config_row(line)
                configs[(card.kset, card.number)] = card
            elif configs and line.strip() and not line.startswith(' '):
                break
    if not configs:
        raise Outg11Error('%s carries no configuration table' % path)
    return Run(path, configs)


def _looks_like_config_row(line):
    if len(line) < 92:
        return False
    if not line[:4].strip().isdigit() or not line[4:8].strip().isdigit():
        return False
    return line[8:15].strip() == '' and line[15:16].strip() != ''


# --------------------------------------------------------------------------
# the autoionization sections
# --------------------------------------------------------------------------

# RCG prints an Auger rate per (level, channel) where a run has
# continuum configurations in it, with `FORMAT 38` at
# `FOR/rcg11kd.f:11373`:
#
#     (8H **AA** ,2A6,1H(,A8,1P,11E9.2/(8H **AA** ,21X,11E9.2))
#
# so the tag is columns 1-8, the target configuration 9-20, its term
# 22-29, and the rates 30 onward, eleven to a line, one per eigenvector
# of the J block.  A continuation blanks columns 9-29 and carries on.
#
# The eigenvalues the rates belong to come from `FORMAT 32` (`:11232`)
# or `33` (`:11235`), which put eleven F9.3 or F9.2 fields from column
# 29 and are told apart by nothing -- both are read as free-width
# 9-column fields.
#
# `FORMAT 39` (`:11381`) writes two ` TOT, E.LT.` summary rows after the
# channels; they are totals over channels below a cut, not channels, and
# are not read.  `**GAA**` is a different section again.
_AA_TAG = ' **AA** '
_AA_START = 29
_AA_CONFIG_NUMBERS = '   CONFIG. NO.'
_AA_HEADING = ' AA(UNITS'
_AA_BLOCK_END = 'AVEEIG'


@dataclass
class AugerBlock:
    """The Auger rates of one J block: one row per channel."""

    kset: int                       # parity set, 1-based, in file order
    j: float
    eigenvalues: list = _field(default_factory=list)
    rates: dict = _field(default_factory=dict)   # {channel: [rate per level]}

    @property
    def nlevels(self):
        return len(self.eigenvalues)


def _aa_channel(line):
    """The channel a `**AA**` row belongs to, or None on a continuation.

    The name is the target configuration and the term it is left in,
    with runs of blanks squeezed to one -- `18.58h (2S) 3H`.  That is the
    string `read_aa.pl:172-174` builds and matches the user's core
    patterns against, so it is built the same way here.
    """
    text = '%s(%s' % (line[8:20], line[21:29])
    return ' '.join(text.split()) or None


def read_auger(path=files.OUTG11, folder=None):
    """Every J block that carries Auger rates, in file order.

    A block is `0  EIGENVALUES (J=...)`, its eigenvalues, and -- where
    the run has a continuum in it -- an `AA(UNITS OF 1/SEC)` section
    whose rows are the channels.  The parity set is counted off the J
    values: RCG walks J upward within a set, so a J that does not
    increase begins the next one.

    Two channels may carry the same label, because the label is an LS
    term and two target states can share one; their rates are added,
    which is what `read_aa.pl:186-195` does.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    out, block, kset, previous = [], None, 0, None
    reading, channel, at = None, None, 0

    for line in files.read_lines(path, width=140):
        if _EIGENVALUES in line[:16] and '(J=' in line:
            j = float(line[23:27])
            if previous is None or j <= previous:
                kset += 1
            previous = j
            block = AugerBlock(kset=kset, j=j)
            out.append(block)
            reading, channel = 'eigenvalues', None
            continue
        if block is None:
            continue
        if line.startswith(_AA_CONFIG_NUMBERS):
            reading = None
        elif line.startswith(_AA_HEADING):
            reading, channel = 'rates', None
        elif _AA_BLOCK_END in line:
            reading, block, channel = None, None, None
        elif reading == 'eigenvalues' and line.strip():
            block.eigenvalues.extend(_floats(line, _AA_START))
        elif reading == 'rates' and line.startswith(_AA_TAG):
            named = _aa_channel(line)
            if named:
                channel, at = named, 0
            if channel is not None:
                at = _put(block.rates.setdefault(channel, []), at,
                          _floats(line, _AA_START))
    return [b for b in out if b.rates]


def _put(row, at, values):
    """Write one row's worth of rates, adding where a slot is taken.

    A continuation row carries on where the last one stopped, so the
    slots are empty and the values are appended.  A **second** row named
    for a channel already seen starts at slot zero and its rates are
    added to the ones there: two target states can carry the same LS
    term, and the rate to that term is the sum of the two.
    """
    for i, value in enumerate(values):
        k = at + i
        if k < len(row):
            row[k] += value
        else:
            row.append(value)
    return at + len(values)


# --------------------------------------------------------------------------
# the parameter sections
# --------------------------------------------------------------------------

_PARAMETERS = 'PARAMETER VALUES IN'


@dataclass(frozen=True)
class ParameterSection:
    """The names RCG gave one parameter card, and the values beside them."""

    name: str            # the card's own name, as ING11 punches it
    names: tuple = ()
    values: tuple = ()


def read_parameter_names(path=files.OUTG11, folder=None):
    """{card name: ParameterSection} for every parameter card of a run.

    An ING11 parameter card is a row of numbers with a one-digit type
    code each and no names at all: the name is *positional* - `F2(11)`,
    `ZETA 3`, `120D1112` - and RCG prints it here and nowhere else.
    Anything that has to know which parameter a value is, as
    `reorder_ing11.py` does, reads it back from OUTG11.

    A section is written by three FORMATs -- `FOR/rcg11kd.f:10876` for a
    single-configuration heading, `:10928` for a configuration-
    interaction one, then `:10880` for the names and `:10882` for the
    values:

        the heading         one line; two for a CI card, the second
                            at column 90 (`89X` in FORMAT 17)
        (blank)
        the names           A8 on a 15-column grid
        (blank)
        the values          F15.5 then F15.6, on the same grid

    The names are taken by splitting on runs of three or more blanks
    rather than by the grid, so that the grid may shift -- and it does,
    past 22 names, where the FORMAT runs out of line groups and reverts
    to its own beginning, which puts a blank line *inside* the names.  A
    blank line is therefore skipped rather than taken as the end; what
    ends the names is the first line that is nothing but numbers, which
    is the values.  No parameter name reads as a number.

    The key is columns 1-18 of the card, blanks off -- the same 18
    columns `ing11.ParameterCard.name` and `ing11.CiCard.name` hold.
    Both headings put them in columns 2-19 (`1X,3A6` and `1X,2A6,A7`),
    and taking the whole field instead would sometimes carry more: RCG
    echoes what it READ from the card with `2A6,A7`, and RCN2 punches
    the parameter count in columns 19-20, so a card with ten or more
    parameters has the first digit of that count inside its A7.  Pr III
    has one, `f5d6d    -f5d7d   1`.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    out = {}
    with open(path, 'r', errors='replace') as fh:
        lines = iter(fh)
        for raw in lines:
            if _PARAMETERS not in raw:
                continue
            name = raw[1:19].strip()
            section = _parameter_section(name, lines, path)
            if out.get(name, section) != section:
                raise Outg11Error(
                    '%s prints two different parameter sections for %s'
                    % (path, name))
            out[name] = section
    return out


def _parameter_section(name, lines, path):
    """One parameter section, the heading line already consumed."""
    names, values, first = [], [], True
    for raw in lines:
        line = raw.rstrip('\n').rstrip('\r')
        if first:
            first = False
            # FORMAT 17's second heading line, at column 90.  FORMAT 14
            # opens with a `/`, so a *blank* line here is the start of
            # the names and anything else is still the heading.
            if line.strip():
                continue
        if not line.strip():
            if values:
                break
            continue
        # A new heading, or the `1` of a Fortran page break, ends the
        # section: RCG follows the last value line with the energy
        # matrix and no blank line between them.
        if _PARAMETERS in line or line[:1] not in (' ', ''):
            break
        fields = [f for f in _SPACES.split(line.strip()) if f]
        numbers = _as_floats(fields)
        if numbers is None:
            if values:
                break
            names.extend(fields)
        else:
            values.extend(numbers)
    if len(names) != len(values):
        raise Outg11Error(
            '%s: the parameter section of %s has %d name(s) and %d value(s)'
            % (path, name, len(names), len(values)))
    return ParameterSection(name=name, names=tuple(names), values=tuple(values))


# Three blanks or more: the A8 name grid leaves seven between names, and
# no name has more than one blank inside it (`ZETA 3`, `T2(D 5)`).
_SPACES = re.compile(r'\s{3,}')


def _as_floats(fields):
    """The fields as numbers, or None if any of them is a name.

    The decimal point is the test, not `float()`: a CI name of exchange
    type reads as an exponential, so `float('672E2424')` succeeds and
    returns infinity.  FORMAT 15 is F15.5 and F15.6 and always writes
    the point; no parameter name carries one.
    """
    if not all('.' in f for f in fields):
        return None
    try:
        return [float(f) for f in fields]
    except ValueError:
        return None


# --------------------------------------------------------------------------
# the J blocks
# --------------------------------------------------------------------------

def _blocks(path, configs, couplings):
    """Stream the J blocks of an OUTG11."""
    block = None
    section = None
    mode = None          # what the lines now arriving belong to
    group_open = False   # inside a run of eigenvector rows
    row_at = 0

    with open(path, 'r', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n').rstrip('\r')

            if _EIGENVALUES in line[:16] and '(J=' in line:
                if block is not None:
                    yield _finish(block, configs)
                # format 30-33: '1  EIGENVALUES ',5X,'(J=',F4.1,')'
                at = line.index('(J=')
                block = JBlock(j=float(line[at + 3:at + 7]))
                section = None
                mode = 'eigenvalues'
                group_open = False
                continue

            if block is None:
                continue

            if line.startswith('   CONFIG. NO.'):
                mode = 'config'
                continue
            if line.startswith('   G-VALUES'):
                mode = 'g'
                continue
            if line.startswith(_EIGENVECTORS):
                name = line[19:25].strip()
                if not name:
                    raise Outg11Error('an EIGENVECTORS heading with no '
                                      'coupling: %r' % line)
                keep = couplings is None or name in couplings
                section = Section(coupling=name) if keep else None
                if keep:
                    block.sections[name] = section
                mode = 'vectors'
                group_open = False
                row_at = 0
                continue

            if mode == 'eigenvalues':
                if line.strip():
                    block.eigenvalues.extend(_floats(line, _DATA))
                continue
            if mode == 'config':
                if line.strip():
                    block.config_numbers.extend(_ints(line, 25))
                continue
            if mode == 'g':
                if line.strip():
                    block.g_values.extend(_floats(line, _DATA))
                continue
            if mode == 'purity':
                # The continuation of a PURITY line starts in column 29
                # for LS, JJ and TRUEjj (rcg11kd.f:1379) and in column 16
                # for couplings 3 to 7 (:1842).  A purity is between 0
                # and 1 in F9.5, so it never touches its neighbour and
                # blanks separate them whichever FORMAT wrote the line.
                # Nothing marks the end of the run of continuations - the
                # TIME line can follow it directly - so a line that is
                # not wholly numbers is the end of it.
                values = _all_floats(line)
                if values is None:
                    mode = 'vectors'
                else:
                    if section is not None:
                        section.purity.extend(values)
                    continue

            if mode != 'vectors' or section is None:
                continue

            if _is_row(line, section.coupling):
                row = _row(line, section.coupling)
                if not group_open:
                    group_open = True
                    row_at = 0
                if row_at < len(section.rows):
                    have = section.rows[row_at]
                    if have.number != row.number:
                        raise Outg11Error(
                            'the %s section of J=%.1f runs out of step at '
                            'state %d' % (section.coupling, block.j,
                                          row.number))
                    section.rows[row_at] = Row(
                        cfg=have.cfg, desc=have.desc, number=have.number,
                        amplitudes=have.amplitudes + row.amplitudes)
                else:
                    section.rows.append(row)
                row_at += 1
                continue

            group_open = False

            if 'PURITY=' in line:
                at = line.index('PURITY=') + 7
                section.average_purity = float(line[at:at + 5])
                section.purity.extend(_floats(line, _DATA))
                mode = 'purity'
                continue

            if _is_heading(line):
                # format 25: I21,9X,11(A8,1X) - the configuration label
                # of each eigenvector, then format 29: 30X,11('(',A8)
                section.column_labels.extend(
                    line[at:at + 8].rstrip()
                    for at in range(30, len(line) - 7, 9))
            elif line.startswith(' ' * 30) and line[30:31] == '(':
                section.column_names.extend(
                    line[at + 1:at + 9].rstrip()
                    for at in range(30, len(line) - 8, 9))

        if block is not None:
            yield _finish(block, configs)


def _is_heading(line):
    """Format 25: a group number right-justified in 21 columns, then names."""
    return (len(line) > 30 and line[:21].strip().isdigit()
            and line[21:30].strip() == '')


def _finish(block, configs):
    """Work out which configuration set a block belongs to.

    A block does not say.  Its rows name a configuration by a number
    that runs within one set, and the same number means a different
    configuration in the other set, so the number alone will not do -
    but the row also carries the 12-character label, and the pair
    settles it.  Where a run has one set only, that set is the answer.

    A configuration under a term cutoff has no name in the table (see
    _config_row), so it cannot be ruled out by its label and is taken to
    match whatever the row says.  In Pr III that is the whole of the
    J=11.5 block, whose one configuration is the truncated 4f3.
    """
    sets = sorted({k for k, _ in configs})
    if len(sets) == 1:
        block.kset = sets[0]
        return block
    any_rows = False
    for section in block.sections.values():
        for row in section.rows:
            any_rows = True
            hits = [k for k in sets
                    if (k, row.cfg) in configs
                    and configs[(k, row.cfg)].label in ('', row.label)]
            if len(hits) == 1:
                block.kset = hits[0]
                return block
    if any_rows:
        raise Outg11Error(
            'the configuration set of the J=%.1f block cannot be told from '
            'its rows: no (number, label) pair is unique to one set' % block.j)
    return block


# --------------------------------------------------------------------------

def _main(argv):
    """Summarize an OUTG11: python -m toolbox.outg11 [path] [coupling ...]."""
    path = argv[1] if len(argv) > 1 else files.OUTG11
    couplings = set(argv[2:]) or None
    run = read_outg11(path)
    print('%s: %d configurations in %d set(s)'
          % (run.path, len(run.configs),
             len({k for k, _ in run.configs})))
    for block in run.blocks(couplings):
        found = ', '.join('%s %d/%d' % (name, s.nstates,
                                        len(s.rows[0].amplitudes)
                                        if s.rows else 0)
                          for name, s in block.sections.items())
        print('  set %d  J=%4.1f  %4d levels  [%s]'
              % (block.kset, block.j, block.nlevels, found))
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(_main(sys.argv))
