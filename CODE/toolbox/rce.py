"""RCE's output: the fitted levels, and the parameters behind them.

Four files, each written by a `FORMAT` named below:

  * **RCEOUT** -- one section per parity set: the spectrum and its
    configuration names, then one line per level (observed energy,
    fitted energy, J, and up to ten leading components), then the
    parameter table with each parameter's flag and value.  This is the
    authoritative list.
  * **LEVELS1** (and `LEVELS1J`, the true-jj pass) -- the same levels
    again, with the Lande g factor and one more digit on the fitted
    energy.  `LEVELS2` and `LEVELS3` are sorted rearrangements of it.
  * **OUTE** -- the fitting listing.  Its last cycle is the only place a
    parameter's standard deviation is printed.
  * **RCEINP.HF** -- the starting (Hartree-Fock) parameter values, which
    is what the `LSF/HF` column of `params.txt` is divided by.

A fifth file is read here although RCE does not write it: **OUTGINE**,
RCG's input deck *for* RCE.  It is the only place a run that has never
been through RCE carries its parameters under their names, and
`read_outgine` is what lets `printout --no-rce` write params.txt at all.

Every one of these has *two or more* layouts, chosen inside RCE by how
large the energies are (`EMAX.LT.1000000.0` and friends).  The fields
after the energies never move, so the layout is settled here by the
column the first decimal point falls in, once per file, and a file that
mixes layouts says so rather than guessing line by line.
"""

from dataclasses import dataclass, field

from . import files
from .fortfmt import read_f, read_i, write_f_wide

OUTGINE = 'OUTGINE'
RCEOUT = 'RCEOUT'
RCEINP = 'RCEINP'
RCEINP_HF = 'RCEINP.HF'
OUTE = 'OUTE'
LEVELS = ('LEVELS1', 'LEVELS2', 'LEVELS3')
LEVELS_JJ = ('LEVELS1J', 'LEVELS2J', 'LEVELS3J')


class RceError(Exception):
    """An RCE output file that does not match the FORMAT that wrote it."""


# --------------------------------------------------------------------------
# column 14: the status character
# --------------------------------------------------------------------------

#: What RCE puts beside a level's observed energy, and what each value
#: means (`FOR/RCE20KD.F:626-660`):
#:
#:   blank   an ordinary experimental level, fitted at weight 1
#:   `*`     no observation; the observed column holds a copy of the
#:           fitted energy and the level is not fitted
#:   `A`-`Z` the **base** of a heteroscedastic group
#:   `a`-`z` another **member** of that group
#:
#: A group is a set of levels known precisely relative to one another
#: but tied to the ground level only loosely -- the `+x` systems of an
#: experimental compilation.  Every one of its levels carries an
#: experimental energy; what the letter says is that the *offset* of
#: that energy is uncertain, so a difference taken across a group
#: boundary is not the well-determined quantity a difference inside one
#: is.  `comparable` is where that distinction is asked about.
UNFITTED = '*'

_UPPER = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
_LOWER = 'abcdefghijklmnopqrstuvwxyz'


def group_letter(flag):
    """The group a status character names, lowercased, or `''`.

    Case says base from member and is asked about by `is_group_base`;
    the letter itself is the group's identity either way, so it is
    folded here and any two levels' groups compare as one spelling.
    """
    text = (flag or '')[:1]
    if text in _UPPER:
        return text.lower()
    return text if text in _LOWER else ''


def is_group_base(flag):
    """Whether a status character is a group's base: an uppercase letter."""
    return (flag or '')[:1] in _UPPER


def is_experimental(flag):
    """Whether the column-14 character says there is an observed energy.

    Blank and every group letter do; `*` does not.  Anything else is
    what RCE now rejects outright (`FOR/RCE20KD.F:649-657`), and it is
    read here the way RCE reads a level it will not fit -- not as an
    observation.
    """
    text = (flag or '').strip()
    return text == '' or bool(group_letter(text))


def comparable(one, other):
    """Whether a difference of two levels' observed energies is meaningful.

    Two absolute levels are on one scale, and so are two levels of one
    group -- that difference is exactly what the group's `sigma_int`
    makes well determined.  Across a group boundary the unknown offset
    enters, and the caller falls back to the fitted energies, which is
    the branch every difference site in this package already has.

    Takes anything carrying `observed` and `group`, so the wrappers
    `conv_out`, `calc_TP_M1E2` and `iden` put around a level answer it
    as well as a `Level` does.
    """
    if not getattr(one, 'observed', False):
        return False
    if not getattr(other, 'observed', False):
        return False
    return getattr(one, 'group', '') == getattr(other, 'group', '')


# --------------------------------------------------------------------------
# the pieces
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Component:
    """One eigenvector component of a fitted level.

    `percent` is what RCE prints: a signed percentage in RCEOUT (the
    amplitude's sign is kept, `IP1`) and a magnitude in LEVELS
    (`IP2`).  `label` is the 12-character basis-state name, which is
    what ties a level back to `genealogy.txt`.
    """

    percent: int
    configuration: str
    label: str

    def __str__(self):
        return '%4d%% %-6s %s' % (self.percent, self.configuration,
                                  self.label)


@dataclass
class Level:
    """One fitted level."""

    e_obs: float                    # observed energy, 0 where there is none
    e_calc: float                   # fitted energy
    j: float
    flag: str = ''                  # column 14: RCE's observed/computed mark
    configuration: str = ''         # the leading component's configuration
    conf0: str = ''                 # LEVELS only: the ING11 label of the
                                    # scheme-invariant dominant configuration
    g: float = None                 # Lande g, LEVELS only
    components: list = field(default_factory=list)
    kset: int = 0                   # parity set, 1 or 2
    index: int = 0                  # position within its file, 1-based

    @property
    def observed(self):
        """True where the level carries an observed energy.

        RCE marks a level it had no observation for by putting the
        computed value in the observed column and an asterisk beside
        it.  A group letter there is not that mark: the level has an
        observed energy and was fitted to it, only under a weight of
        its own.  See `is_experimental`.
        """
        return is_experimental(self.flag)

    @property
    def group(self):
        """The heteroscedastic group this level belongs to, or `''`."""
        return group_letter(self.flag)

    @property
    def is_base(self):
        """Whether this level is its group's base."""
        return is_group_base(self.flag)

    @property
    def absolute(self):
        """Whether the observed energy is on the deck's absolute scale.

        The second of the two questions `observed` used to answer on
        its own: an ordinary experimental level is absolute, a group
        level has an observed energy whose offset is not known.
        """
        return self.observed and not self.group

    @property
    def leading(self):
        return self.components[0] if self.components else None

    def __str__(self):
        return '%12.3f J=%s %s' % (self.e_calc, _j_text(self.j),
                                   self.configuration)


@dataclass(frozen=True)
class Parameter:
    """One fitted parameter, as RCE names it.

    The name is positional -- `F2(33)`, `ZETA 3`, `120D1112` -- and
    means nothing until a configuration's ING11 slots have been tied to
    their `nl`.  That is what `toolbox.parnames` does; this module only
    reports what the file says.
    """

    name: str
    flag: int                       # LF: 0 fixed, negative free, group code
    value: float
    maximum: float = 0.0
    configurations: tuple = ()      # the pair a CI parameter belongs to
    sigma: float = None             # standard deviation, OUTE only

    @property
    def is_ci(self):
        return len(self.configurations) == 2

    def __str__(self):
        return '%-10s %12.4f' % (self.name, self.value)


@dataclass(frozen=True)
class Group:
    """One heteroscedastic group, as its card in the deck gives it.

    The letter lives on the level records; the two standard deviations
    that say what the group is worth have nowhere to live there, so
    they come on a card block of their own at the end of the section
    (`FOR/RCE20KD.F:2492`, SUBROUTINE RELRD).

    `sigma_base` is `None` where the field was left blank, which is how
    the deck says the group has no measured tie to the ground level at
    all -- its offset is unknown rather than poorly known, and its base
    then contributes no residual.  That is why a blank and an explicit
    zero are different things here, the second being rejected by RCE:
    zero reads as "known exactly", which in the weight formula means
    precisely that.
    """

    letter: str                     # lowercase, the group's identity
    sigma_base: float = None        # None: no tie to ground
    sigma_int: float = 0.0

    @property
    def untied(self):
        return self.sigma_base is None

    def __str__(self):
        return '%s: sigma_base=%s sigma_int=%.3f' % (
            self.letter,
            'untied' if self.untied else '%.3f' % self.sigma_base,
            self.sigma_int)


@dataclass
class Cset:
    """One parity set of an RCE run."""

    index: int                      # 1 or 2, RCE's CSET
    spectrum: str = ''
    configurations: list = field(default_factory=list)
    levels: list = field(default_factory=list)
    parameters: list = field(default_factory=list)
    groups: dict = field(default_factory=dict)   # {letter: Group}

    def by_j(self):
        """{J: [levels]}, each list in the file's own order."""
        out = {}
        for level in self.levels:
            out.setdefault(level.j, []).append(level)
        return out

    def __len__(self):
        return len(self.levels)


# --------------------------------------------------------------------------
# shared field reading
# --------------------------------------------------------------------------

def _j_text(j):
    return '%d/2' % round(2 * j) if round(2 * j) % 2 else '%d' % round(j)


def _point_column(line, start, stop):
    """The 1-based column of the first '.' in a slice, or None.

    This is how a file's layout is settled: the energy fields differ in
    width between RCE's variants but the decimal point of the first one
    is at a different column in each, and everything after them is
    fixed.
    """
    at = line.find('.', start, stop)
    return None if at < 0 else at + 1


# The component group, in the two shapes RCE writes it: RCEOUT's
# `(I4,1H ,A6,1X,A12)` (`FOR/RCE20KD.F:2074`) and LEVELS's
# `(I4,2H% ,A6,1X,A12)` (`:2033`), which spends one more column on the
# per-cent sign it prints and RCEOUT does not.
RCEOUT_COMPONENT = (24, slice(5, 11), slice(12, 24))
LEVELS_COMPONENT = (25, slice(6, 12), slice(13, 25))


def _components(line, start, shape, nfields=10):
    """The components of one level, beginning at a column.

    A component whose field is blank ends the list: RCE writes only as
    many as the level has, and the record simply stops.
    """
    width, cfg, label = shape
    out = []
    for i in range(nfields):
        at = start + width * i
        chunk = line[at:at + width]
        if not chunk.strip():
            break
        out.append(Component(percent=read_i(chunk[:4]),
                             configuration=chunk[cfg].strip(),
                             label=chunk[label].rstrip()))
    return out


# --------------------------------------------------------------------------
# RCEOUT
# --------------------------------------------------------------------------

# `FORMAT (F13.6,A1,F10.3,I4,A2,10(I4,1H ,A6,1X,A12:))` at
# `FOR/RCE20KD.F:2074`, and the wide-energy variant `1601` at `:2075`,
# `F13.4,A1,F10.1`.  Both are 24 columns before the components, so only
# the split of columns 1-24 differs -- and the decimal point of the
# first field says which: column 7 for F13.6, column 9 for F13.4.
_RCEOUT_LAYOUTS = {7: (13, 6, 10, 3), 9: (13, 4, 10, 1)}

# `FORMAT (A10,I5,2F14.6,14X,A6,A8)` at `FOR/RCE20KD.F:3540`.
_PAR_NAME, _PAR_FLAG = slice(0, 10), slice(10, 15)
_PAR_VALUE, _PAR_MAX = slice(15, 29), slice(29, 43)
_PAR_C1, _PAR_C2 = slice(57, 63), slice(63, 71)


def read_rceout(path=RCEOUT, folder=None):
    """Read RCEOUT: a list of `Cset`, in the order RCE wrote them.

    A section is three header lines -- the run control card, the
    spectrum with its configuration names, and RCE's scaling card --
    then the levels, then `PARAMETER FLAG ...` and the parameter table.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    lines = files.read_lines(path, width=140)
    out, i = [], 0
    while i < len(lines):
        while i < len(lines) and not _is_control(lines[i]):
            i += 1
        if i >= len(lines):
            break
        cset = Cset(index=len(out) + 1)
        i += 1
        cset.spectrum, cset.configurations, i = _read_names(lines, i)
        i += 1                                  # RCE's scaling card
        cset.levels, i = _read_rceout_levels(lines, i, cset.index, path)
        cset.parameters, i = _read_parameters(lines, i)
        cset.groups, i = read_groups(lines, i)
        out.append(cset)
    if not out:
        raise RceError('%s carries no parity set' % path)
    return out


def _is_control(line):
    """RCE's run control card: eight small integers and two reals."""
    fields = line.split()
    if len(fields) < 8:
        return False
    try:
        return (int(fields[0]) in (1, 2) and int(fields[1]) == 0
                and float(fields[4]) is not None)
    except ValueError:
        return False


def _read_names(lines, i):
    """The spectrum and its configuration names.

    RCEOUT writes them with `FORMAT (7(A8,2X))` (`FOR/RCE20KD.F:479`)
    and LEVELS with `FORMAT (A8,10(2X,A8)/(8X,10(2X,A8)))` (`:1305`).
    The two agree on the first line -- spectrum in columns 1-8, names in
    11-18, 21-28 and so on -- and differ only in their continuation,
    which RCEOUT starts at column 1 and LEVELS indents by eight.  Taking
    a non-blank first field of a continuation as a name covers both.

    A configuration label never carries a decimal point, so the first
    line that does is the next record: the scaling card in RCEOUT, a
    level in LEVELS.
    """
    spectrum = lines[i][:8].strip()
    names = _names_of(lines[i], first=False)
    i += 1
    while i < len(lines) and lines[i].strip() and '.' not in lines[i]:
        names.extend(_names_of(lines[i], first=True))
        i += 1
    return spectrum, names, i


def _names_of(line, first):
    """The A8 fields of one names line, on the 10-column grid."""
    out = []
    for k in range(0 if first else 1, (len(line) + 9) // 10):
        name = line[10 * k:10 * k + 8].strip()
        if name:
            out.append(name)
    return out


def _read_rceout_levels(lines, i, kset, path):
    """Every level line of one section, up to the parameter table."""
    out, layout = [], None
    while i < len(lines) and 'PARAMETER FLAG' not in lines[i]:
        line = lines[i]
        if not line.strip():
            i += 1
            continue
        column = _point_column(line, 0, 13)
        if column is None:
            break
        if layout is None:
            layout = _RCEOUT_LAYOUTS.get(column)
            if layout is None:
                raise RceError(
                    '%s line %d: the first energy field has its decimal '
                    'point in column %d, which is neither of RCE\'s two '
                    'layouts (7 or 9)' % (path, i + 1, column))
        w1, d1, w2, d2 = layout
        # `I4,A2`: RCE halves the doubled J and blanks the A2 where it
        # comes out whole, so '/2' in columns 29-30 is what says the I4
        # is 2J rather than J (`FOR/RCE20KD.F:2055-2062`).
        two_j = read_i(line[24:28])
        if line[28:30] != '/2':
            two_j *= 2
        out.append(Level(e_obs=read_f(line[:w1], d1),
                         flag=line[w1:w1 + 1],
                         e_calc=read_f(line[w1 + 1:w1 + 1 + w2], d2),
                         j=two_j / 2.0,
                         components=_components(line, 30, RCEOUT_COMPONENT),
                         kset=kset, index=len(out) + 1))
        if out[-1].components:
            out[-1].configuration = out[-1].components[0].configuration
        i += 1
    return out, i


def level_lines(lines):
    """Yield `(line number, parity set, position)` for every level line.

    `read_rceout` builds objects; this walks the same file for a
    **rewriter**, which has to put a value back into the line it came
    from and leave every other column of it alone.  The two follow the
    same rules -- a section is a control card, two header lines, then
    levels until `PARAMETER FLAG` -- so the n-th line this yields for a
    parity set is the n-th level `read_rceout` returns for it.
    """
    i, kset = 0, 0
    while i < len(lines):
        while i < len(lines) and not _is_control(lines[i]):
            i += 1
        if i >= len(lines):
            return
        kset += 1
        i += 1
        _, _, i = _read_names(lines, i)
        i += 1                                  # RCE's scaling card
        position = 0
        while i < len(lines) and PARAMETER_HEADING not in lines[i]:
            if not lines[i].strip():
                i += 1
                continue
            if _point_column(lines[i], 0, 13) is None:
                break
            position += 1
            yield i, kset, position
            i += 1
        while i < len(lines) and PARAMETER_HEADING not in lines[i]:
            i += 1


def rceout_layout(line):
    """The energy-field widths of one RCEOUT level line, or None."""
    return _RCEOUT_LAYOUTS.get(_point_column(line, 0, 13))


def set_observed(line, energy, flag=' '):
    """A level line carrying an observed energy, and its status character.

    RCE puts the computed value in the observed column and an asterisk
    in the column behind it where it had no observation
    (`FOR/RCE20KD.F:2074`); giving the level an energy means writing the
    energy and putting the status the caller asks for in that column.
    Nothing else on the line moves -- the fitted energy, J and the whole
    eigenvector stay as they are, which is what
    `subst_exp_levs.pl:385-388` rebuilds field by field and can
    therefore only do for one of RCE's two layouts.

    `flag` defaults to a blank, which is the ordinary absolute level
    this has always written and is what clears an asterisk.  A group
    letter is passed when the energy being carried in belongs to a
    group: writing a blank there would silently promote a relative
    energy to an absolute one in the target deck.
    """
    layout = rceout_layout(line)
    if layout is None:
        raise RceError('%r is not a level line' % line[:30])
    if len(flag) != 1 or not (flag == ' ' or flag == UNFITTED
                              or group_letter(flag)):
        raise RceError('%r is not a status character for column 14' % flag)
    w1, d1, _, _ = layout
    return write_f_wide(energy, w1, d1) + flag + line[w1 + 1:]


def _read_parameters(lines, i):
    """The `PARAMETER FLAG ...` table that closes a section."""
    out = []
    while i < len(lines) and 'PARAMETER FLAG' not in lines[i]:
        i += 1
    if i >= len(lines):
        return out, i
    i += 1
    while i < len(lines) and lines[i][:10].strip() \
            and not _is_control(lines[i]) and not is_group_heading(lines[i]):
        line = lines[i]
        pair = tuple(x for x in (line[_PAR_C1].strip(),
                                 line[_PAR_C2].strip().lstrip('-')) if x)
        out.append(Parameter(name=line[_PAR_NAME].strip(),
                             flag=read_i(line[_PAR_FLAG]),
                             value=read_f(line[_PAR_VALUE], 6),
                             maximum=read_f(line[_PAR_MAX], 6),
                             configurations=pair if len(pair) == 2 else ()))
        i += 1
    return out, i


# --------------------------------------------------------------------------
# the group card block
# --------------------------------------------------------------------------

# RCE reads it in SUBROUTINE RELRD (`FOR/RCE20KD.F:2492`) and writes it
# in RELWR (`:2614`), immediately behind a section's parameter cards --
# the only place it can go, the parameter cards having been the last
# thing ever read from unit 20, so that a deck written before the
# feature existed hits end of file here and yields no groups.
#
#      GROUP SIGMA_BASE  SIGMA_INT
#     X         500.000     20.000
#     Z                     15.000
#
# The heading is any line carrying 'GROUP'; a card is the letter in
# column 1 and the two sigmas in columns 2-13 and 14-25; a blank line, a
# '-1', a line that is not a card, or end of file closes the block.
GROUP_HEADING = ' GROUP SIGMA_BASE  SIGMA_INT'

_GROUP_LETTER = slice(0, 1)
_GROUP_BASE, _GROUP_INT = slice(1, 13), slice(13, 25)


def is_group_heading(line):
    """Whether a line opens the group card block.

    The same test RELRD makes -- the word alone -- so that a heading
    written by hand in the author's own spacing is still recognised.
    """
    return 'GROUP' in line.upper()


def read_groups(lines, i=0):
    """The group card block beginning at `lines[i]`: ({letter: Group}, i).

    An empty dict and an unmoved `i` where there is no block, which is
    every deck written before the feature and every section of a
    lettered deck that has no group of its own.  The sigmas are read
    free-form out of their fixed columns, as RELRD reads them, so that
    a field punched `500` means 500 and not 0.5.
    """
    out = {}
    if i >= len(lines) or not is_group_heading(lines[i]):
        return out, i
    i += 1
    while i < len(lines):
        line = lines[i]
        if not line.strip() or line[:2] == '-1':
            return out, i + 1
        letter = group_letter(line[_GROUP_LETTER])
        if not letter:
            return out, i
        if letter in out:
            raise RceError('two group cards for group %s: %r'
                           % (letter, line[:40]))
        out[letter] = Group(letter=letter,
                            sigma_base=_group_sigma(line[_GROUP_BASE], letter,
                                                    'SIGMA_BASE', line),
                            sigma_int=_group_sigma(line[_GROUP_INT], letter,
                                                   'SIGMA_INT', line) or 0.0)
        i += 1
    return out, i


def _group_sigma(field, letter, name, line):
    """One sigma field: None where it is blank, as RELRD takes it."""
    text = field.strip()
    if not text:
        return None
    try:
        return float(text)
    except ValueError:
        raise RceError('the %s of group %s is not a number: %r'
                       % (name, letter, line[:40]))


def group_lines(groups):
    """The group card block for `groups`, as RELWR writes it.

    `(A1,2F12.3)`, or `(A1,12X,F12.3)` for an untied group, whose
    SIGMA_BASE field comes back blank -- exactly how it was spelt on
    input, which is what makes the block round-trip.  Nothing at all
    where there is no group.
    """
    if not groups:
        return []
    out = [GROUP_HEADING]
    for letter in sorted(groups):
        group = groups[letter]
        out.append('%s%s%12.3f'
                   % (letter.upper(),
                      ' ' * 12 if group.untied else '%12.3f' % group.sigma_base,
                      group.sigma_int))
    return out


# --------------------------------------------------------------------------
# editing the parameter table of an RCEINP in place
# --------------------------------------------------------------------------

# RCE writes a parameter with `FORMAT (A10,I5,2F14.6,14X,A6,A8)` at
# `FOR/RCE20KD.F:1470` and, where a denominator is in play, with
# `(A10,I5,3F14.6,A6,A8)` at `:1464`.  Both put the name in columns
# 1-10, the least-squares flag in 11-15 and the value in 16-29, and both
# put the configuration pair of a CI parameter in 58-71; only what lies
# between differs.  Nothing past column 29 is ever rewritten here, so
# the two layouts need not be told apart.
#
# RCE reads the table back with `FORMAT (A10,I5,2F14.4)` (`:642`), four
# decimals rather than six.  The field carries its decimal point, so
# reading it as F14.4 gives the same number either way, and the value is
# written as RCE writes it.
PARAMETER_HEADING = 'PARAMETER FLAG'
FLAG_FIXED = 100

_VALUE_WIDTH, _VALUE_DECIMALS = 14, 6


def is_control_card(line):
    """RCE's run control card, which is what ends a parameter table.

    `subst_exp_levs.pl:125` tests it as `^ {4}\\d {4}\\d {4}\\d`; the
    same card is what `_is_control` recognises for `read_rceout`, so
    there is one test here and not two.
    """
    return _is_control(line)


def is_parameter_line(line):
    """Whether a line **inside** a parameter table is a parameter.

    Only the name field is tested, because the caller is expected to
    have opened the table on its heading and to close it on the next
    control card.  Outside a table this says nothing useful: the
    configuration-names line that follows a control card also has
    something in columns 1-10.
    """
    return bool(line[_PAR_NAME].strip()) and not _is_control(line)


def parameter_flag(line):
    """The least-squares flag of a parameter line, LF."""
    return read_i(line[_PAR_FLAG])


def parameter_value(line):
    """The value of a parameter line, F14.6 in columns 16-29."""
    return read_f(line[_PAR_VALUE], _VALUE_DECIMALS)


def parameter_group(line):
    """The fitting group a parameter belongs to.

    RCE keeps the group in the magnitude of the flag and the free/fixed
    choice in its sign, and reserves 100 for "fixed on its own", so a
    magnitude of 100 or more is not a group (`scale_param.pl:61-62`).
    """
    group = abs(parameter_flag(line))
    return group if group < FLAG_FIXED else 0


def is_ci_parameter(line):
    """Whether the parameter belongs to a pair of configurations."""
    return bool(line[_PAR_C1].strip() and line[_PAR_C2].strip())


def edit_parameter_line(line, value=None, fixed=None):
    """A parameter line with a new value, a new flag, or both.

    Everything from column 30 on -- the maximum, the denominator where
    there is one, and the configuration pair -- is carried across
    untouched, so a line this does not change comes back byte for byte
    the same.

    `fixed=True` sets the flag to RCE's fixed code, keeping the sign it
    already had, which is what `scale_param.pl:65-67` does.
    """
    line = line.ljust(29)
    head = line[_PAR_NAME] + line[_PAR_FLAG]
    if fixed:
        flag = -FLAG_FIXED if parameter_flag(line) < 0 else FLAG_FIXED
        head = line[_PAR_NAME] + '%5d' % flag
    body = line[_PAR_VALUE]
    if value is not None:
        body = write_f_wide(value, _VALUE_WIDTH, _VALUE_DECIMALS)
    return head + body + line[29:]


# --------------------------------------------------------------------------
# LEVELS1 / LEVELS2 / LEVELS3, and their true-jj twins
# --------------------------------------------------------------------------

# `FORMAT (F13.6,A1,F11.4,F10.6,F7.3,1X,A6,A8,4X,3(I4,2H% ,A6,1X,A12))`
# at `FOR/RCE20KD.F:2033`; `1301` (`F13.4,A1,F10.1,F11.4,F7.3`) at
# `:2036` for a run whose energies pass 10**6, and `14`
# (`F13.3,A1,F10.0,F11.3,F7.3`) at `:2044` for a scaled run.  Columns 36
# onward are the same in all three, so only the point column separates
# them: 7 for F13.6, 9 for F13.4 and 10 for F13.3.
_LEVELS_LAYOUTS = {7: ((13, 6), (11, 4), (10, 6)),
                   9: ((13, 4), (10, 1), (11, 4)),
                   10: ((13, 3), (10, 0), (11, 3))}


def read_levels(path='LEVELS1', folder=None):
    """Read a LEVELS file: a list of `Cset`, one per section.

    LEVELS1 has one section per parity set per coupling scheme, so a
    two-parity run gives four; LEVELS1J, written by the true-jj pass,
    has one per parity set.  A section begins with the spectrum and its
    configuration names.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    lines = files.read_lines(path, width=140)
    out, i, layout = [], 0, None
    while i < len(lines):
        if not lines[i].strip():
            i += 1
            continue
        if _point_column(lines[i], 0, 13) is None:
            cset = Cset(index=len(out) + 1)
            cset.spectrum, cset.configurations, i = _read_names(lines, i)
            out.append(cset)
            continue
        if not out:
            raise RceError('%s begins with a level, not a section header'
                           % path)
        level, layout, i = _read_levels_line(lines, i, layout, path,
                                             out[-1].index)
        level.index = len(out[-1].levels) + 1
        out[-1].levels.append(level)
    return out


def _j_of(field, path, i):
    """The A8 J field of a LEVELS line: '  J= 0.5' -> 0.5.

    RCE writes it through `FORMAT (4X,F4.1)` on the way in
    (`FOR/RCE20KD.F:2056`), so the value is a plain decimal; only the
    'J=' in front of it has to come off.
    """
    text = field.strip()
    if text.startswith('J='):
        text = text[2:]
    try:
        return float(text)
    except ValueError:
        raise RceError('%s line %d: %r is not a J field' % (path, i + 1,
                                                            field))


def _read_levels_line(lines, i, layout, path, kset):
    line = lines[i]
    column = _point_column(line, 0, 13)
    if layout is None:
        layout = _LEVELS_LAYOUTS.get(column)
        if layout is None:
            raise RceError(
                '%s line %d: the first energy field has its decimal point '
                'in column %d, which is none of RCE\'s three layouts '
                '(7, 9 or 10)' % (path, i + 1, column))
    (w1, d1), (w2, d2), (w3, d3) = layout
    at = w1 + 1
    level = Level(e_obs=read_f(line[:w1], d1),
                  flag=line[w1:w1 + 1],
                  e_calc=read_f(line[at:at + w2], d2),
                  g=read_f(line[35:42], 3),
                  j=_j_of(line[49:57], path, i),
                  configuration=line[43:49].strip(),
                  # The same A6 field, kept under its own name: in a
                  # LEVELS file it is not the leading component's
                  # configuration but RCE's scheme-invariant dominant
                  # one (`NJMS`, `FOR/RCE20KD.F:1157-1204`), chosen from
                  # full-precision coefficients before any output pass.
                  conf0=line[43:49].strip(),
                  components=_components(line, 61, LEVELS_COMPONENT, 3),
                  kset=kset)
    return level, layout, i + 1


# --------------------------------------------------------------------------
# OUTE -- the only place a standard deviation is printed
# --------------------------------------------------------------------------

def read_oute(path=OUTE, folder=None):
    """The last fitting cycle of each parity set: {cset: [Parameter]}.

    A cycle is a `CYC NO=` line, then a `SLATER PRMTR ... CSET= n`
    heading that names the parity set, a blank, one line per parameter,
    and `Iteration No n finished`.  The parity is read from the heading
    rather than counted off the cycle numbers, which is what the Perl
    does and what makes it depend on the fit converging in a particular
    number of cycles.  Later cycles overwrite earlier ones, so what is
    left at the end is the converged fit.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    lines = files.read_lines(path, width=100)
    out = {}
    for i, line in enumerate(lines):
        if 'SLATER PRMTR' not in line:
            continue
        out[_cset_of(line)] = _read_oute_block(lines, i + 2)
    return out


def _cset_of(line):
    _, _, tail = line.partition('CSET=')
    return read_i(tail[:5]) if tail.strip() else 1


def _read_oute_block(lines, i):
    """`(.{15})(.{4})(.{13}).{39}(.{13})` -- name, LF, new x, ..., SDX."""
    out = []
    while i < len(lines) and 'finished' not in lines[i]:
        line = lines[i]
        if line[:15].strip():
            out.append(Parameter(name=line[:15].strip(),
                                 flag=read_i(line[15:19]),
                                 value=read_f(line[19:32], 4),
                                 sigma=read_f(line[71:84], 3)))
        i += 1
    return out


# `FORMAT (/' RELATIVE-ENERGY GROUPS - CYCLE',I3,':  STDEV=',F10.3,
# '   SIGMA_0=',F10.3,A10/' GROUP      W_BASE       W_INT')` and
# `(5X,A1,2F12.6)` at `FOR/RCE20KD.F:3148-3155`, SUBROUTINE RELSIG.
# Both the raw STDEV and the adopted SIGMA_0 are on the page because the
# ratchet that produces the second is irreversible, so a cycle that
# clamped it has to be visible.
_GROUP_CYCLE = 'RELATIVE-ENERGY GROUPS - CYCLE'


@dataclass
class GroupCycle:
    """One cycle's group summary, as OUTE prints it."""

    cycle: int
    stdev: float                    # the raw value the cycle achieved
    sigma_0: float                  # what the next cycle adopts, ratcheted
    frozen: bool                    # the weights are held from here on
    weights: dict = field(default_factory=dict)   # {letter: (w_base, w_int)}


def read_oute_groups(path=OUTE, folder=None):
    """The last group summary of OUTE, or None where there is no group.

    `sigma_0` and the weights are derived, not data: they are not
    written to RCEOUT, and a run resumed from that deck re-anneals to
    its own fixed point.  OUTE is where they are recorded, so it is
    where they are read.
    """
    if folder is not None:
        path = files.resolve(path, folder, required=False)
        if path is None:
            return None
    lines = files.read_lines(path, width=100)
    out = None
    for i, line in enumerate(lines):
        if _GROUP_CYCLE not in line:
            continue
        head = line.partition(_GROUP_CYCLE)[2]
        # The two labels are named rather than counted off: they are
        # what makes the line readable, and reading it by them means a
        # change of spacing in the FORMAT cannot silently misparse it.
        out = GroupCycle(cycle=read_i(head[:3]),
                         stdev=read_f(head.partition('STDEV=')[2][:10], 3),
                         sigma_0=read_f(head.partition('SIGMA_0=')[2][:10], 3),
                         frozen='FROZEN' in head)
        for text in lines[i + 2:]:
            letter = group_letter(text[5:6])
            if not letter:
                break
            out.weights[letter] = (read_f(text[6:18], 6),
                                   read_f(text[18:30], 6))
    return out


# --------------------------------------------------------------------------
# RCEINP.HF -- the starting values
# --------------------------------------------------------------------------

def read_hf(path=RCEINP_HF, folder=None):
    """The Hartree-Fock parameter values: {cset: [Parameter]}.

    RCEINP.HF is an RCEINP written before any fitting, so its parameter
    tables are in RCEOUT's own layout; only the values differ.
    """
    if folder is not None:
        path = files.resolve(path, folder, required=False)
        if path is None:
            return {}
    lines = files.read_lines(path, width=140)
    out, i, kset = {}, 0, 0
    while i < len(lines):
        if 'PARAMETER FLAG' in lines[i]:
            kset += 1
            out[kset], i = _read_parameters(lines, i)
        else:
            i += 1
    return out


# --------------------------------------------------------------------------
# OUTGINE -- the parameters of a run RCE has not seen
# --------------------------------------------------------------------------

# RCG writes one section per parity set in SUBROUTINE RCEINP
# (`FOR/rcg11kd.f:15704-15855`).  Only four of its records are wanted
# here, and each is on a 10-column grid:
#
#   `(4I5,F10.3,2I5,F10.4)`  the section head, which `_is_control`
#                            recognises exactly as it does RCEOUT's
#   `(A8,2X,6(A6,A4))` +     the spectrum and the configuration names
#     `(7(A6,A4))`
#   `(7A10)`   `:15739`      the parameter names, NPAR of them
#   ...                      the eigenvalues and their indices, per J
#   `(7I10)`   `:15833`      MULTQQ, the least-squares flag per parameter
#   `(7F10.4)` `:15835`      XMAX, all zero
#   `(7A10)`   `:15843`      the parameter values, as text
#   `(3F10.1,F5.2,I5,F10.3,F10.2,2F5.2,I2)`   the convergence card
#   `(4I5,...)` with -1      the section terminator
#
# The three parameter records are the last three blocks before the
# convergence card, so they are found by counting back from the
# terminator rather than by walking the eigenvalues, whose record count
# depends on the number of J blocks and on which of the three energy
# FORMATs RCG chose.  The all-zero block in the middle is checked, and
# is what says the count landed where it was meant to.
_OUTGINE_PER_LINE = 7
_OUTGINE_WIDTH = 10


def read_outgine(path=OUTGINE, folder=None):
    """The parameters RCG punched for RCE: {cset: [Parameter]}.

    The values are the ING11 ones -- already multiplied by the HR scale
    factors of the parameter cards, exactly as RCEINP.HF's are -- and
    the flag is MULTQQ, the group code RCE will fit them under.  There
    are no standard deviations and no CI configuration columns; the two
    configurations of a CI integral are decoded from its name instead
    (`parnames._ci`).
    """
    if folder is not None:
        path = files.resolve(path, folder, required=False)
        if path is None:
            return {}
    lines = files.read_lines(path, width=140)
    out, i = {}, 0
    while i < len(lines):
        if not _is_control(lines[i]):
            i += 1
            continue
        names, i = _outgine_names(lines, i + 1)
        end = i
        while end < len(lines) and lines[end].strip() != '-1':
            end += 1
        if end >= len(lines):
            raise RceError('%s: the section that begins with %s has no '
                           'terminator' % (path, names[0] if names else '?'))
        out[len(out) + 1] = _outgine_parameters(lines, end, names, path)
        i = end + 1
    return out


def _outgine_names(lines, i):
    """The parameter names of one section, and the line after them.

    They follow the configuration names on the same 10-column grid, so
    the two blocks are told apart by their content rather than by a
    count: the parameter block starts a record of its own and its first
    field is always the EAV of the first configuration, which no
    configuration label can be.  It ends at the first eigenvalue, which
    is the first line carrying a decimal point.
    """
    while i < len(lines) and not lines[i][:_OUTGINE_WIDTH].startswith('EAV'):
        if not lines[i].strip() or '.' in lines[i]:
            raise RceError('a parameter block with no EAV record')
        i += 1
    out = []
    while i < len(lines) and lines[i].strip() and '.' not in lines[i]:
        out.extend(_outgine_fields(lines[i]))
        i += 1
    return out, i


def _outgine_fields(line):
    """The non-blank 10-column fields of one record."""
    return [f for f in (line[k:k + _OUTGINE_WIDTH].strip()
                        for k in range(0, len(line), _OUTGINE_WIDTH)) if f]


def _outgine_parameters(lines, end, names, path):
    """The flags and the values, counted back from the terminator."""
    n = -(-len(names) // _OUTGINE_PER_LINE)         # records per block
    first = end - 1 - 3 * n                         # the convergence card is
    if first < 0:                                   # the line before `first`
        raise RceError('%s: a section too short to hold %d parameters'
                       % (path, len(names)))
    flags = _outgine_block(lines, first, n)
    zeros = _outgine_block(lines, first + n, n)
    values = _outgine_block(lines, first + 2 * n, n)
    for block in (flags, zeros, values):
        if len(block) != len(names):
            raise RceError('%s: %d parameter names but a block of %d'
                           % (path, len(names), len(block)))
    if any(float(z) for z in zeros):
        raise RceError('%s: the record before the parameter values is not '
                       'the block of zeros it should be' % path)
    return [Parameter(name=name, flag=int(flag), value=float(value))
            for name, flag, value in zip(names, flags, values)]


def _outgine_block(lines, first, n):
    out = []
    for line in lines[first:first + n]:
        out.extend(_outgine_fields(line))
    return out
