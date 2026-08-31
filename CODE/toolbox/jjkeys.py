"""Reading true-jj_keys.txt, and decoding a true jj label with it.

A label in the TRUEjj eigenvector sections of OUTG11 codes the basis
state one character per quantity, taking the occupied subshells in order
and leaving the empty ones out:

    Ai   the jj state of subshell i, keyed by configuration and subshell
    Ji   2 * the total J of subshell i
    Ki   2 * the running J after coupling it in

With n occupied subshells the label reads

    n=1   A1 J1              n=2   A1 J1 A2 J2
    n=3   A1 J1 A2 J2 K2 A3 J3

The first and last running J are never coded - the first is J1 and the
last is the J of the block - so within one configuration and one J block
equal labels mean the same column.  Eight characters do not stretch
further: with four occupied subshells the Ki are dropped, with five or
more the Ji as well, and such labels need not be unique.  decode() says
so rather than guessing.

The writers are SUBROUTINE LSJJKP (FOR/rcg11kd.f:16712-16918) for this
file and CHARACTER FUNCTION LSJJCC (:16682-16710) for the alphabet.
The entry FORMATs followed here are

    20  the heading of one subshell of one configuration   :16817
    23  a code that carries one partition                  :16896
    24/25 a code that carries two, one per J               :16899, :16902

A subshell with l=0, with one electron or with one hole is not rotated,
because its J1J2 state already is a jj state.  It then has one LS term
and so one code, which covers both values of J, and its entry gives
both - hence the two-line form.  The Ji character of the label says
which applies, and that is what decode() uses to pick.
"""

from dataclasses import dataclass

from . import files

# LSJJCC's alphabet: 1..9, then A..Z, then a..z, then punctuation.  '0'
# is left out so that a code is never taken for a blank or a zero, '*'
# means "no code left", and the brackets and the space are left out
# because the FORMATs around a label use them.  85 codes.
_PUN = '#%&@+-=?<>~^_|;:!/\\.,"${'
ALPHABET = ('123456789'
            + ''.join(chr(55 + n) for n in range(10, 36))
            + ''.join(chr(61 + n) for n in range(36, 62))
            + _PUN)

NO_CODE = '*'

# Where the key ran out, LSJJKP writes this and the labels around it need
# not be unique.
_VALUE = {c: i + 1 for i, c in enumerate(ALPHABET)}


class JJKeyError(Exception):
    """true-jj_keys.txt did not read as SUBROUTINE LSJJKP writes it."""


@dataclass(frozen=True)
class Part:
    """The j- or the j+ half of a subshell."""
    l: str               # subshell letter, lower case
    sign: str            # '-' for j = l-1/2, '+' for j = l+1/2
    nel: int
    seniority: int       # 0 where LSJJKP prints none, i.e. where J alone
                         # tells the states of that j**N apart
    j: float             # J of this half, NOT doubled: the key file
                         # prints J itself, the labels double it


@dataclass(frozen=True)
class Entry:
    """One code of one subshell of one configuration.

    `parts` holds the partition where the code carries just one.  Where
    the subshell was not rotated the code covers both values of J, and
    `by_j` maps each J to its own partition; `parts` is then empty.
    """
    code: str
    parts: tuple = ()
    by_j: tuple = ()     # of (J, (Part, ...)), in the order printed

    def partition(self, j=None):
        """The partition of this code, at subshell J `j` where it splits."""
        if not self.by_j:
            return self.parts
        if j is None:
            raise JJKeyError(
                'code %r covers two values of J; say which one' % self.code)
        for value, parts in self.by_j:
            if abs(value - j) < 1e-6:
                return parts
        raise JJKeyError('code %r carries no J=%.1f' % (self.code, j))


@dataclass(frozen=True)
class Subshell:
    """The key to one subshell of one configuration."""
    config: str          # the 18-character 3A6 name, spectrum and label
    index: int           # which subshell of the configuration, 1-based
    l: str               # its letter, lower case
    nel: int             # its occupation
    codes: dict          # {code character: Entry}


@dataclass(frozen=True)
class Position:
    """One decoded position of a label."""
    kind: str            # 'A', 'J' or 'K'
    char: str
    subshell: int        # which occupied subshell it belongs to, 1-based
    value: int = 0       # 2J for a J or K position, 0 for an A


class JJKeys:
    """The key file of one run."""

    def __init__(self, path, subshells):
        self.path = path
        self.subshells = subshells      # {(config, index): Subshell}

    def __len__(self):
        return len(self.subshells)

    def __iter__(self):
        return iter(self.subshells.values())

    def subshell(self, config, index):
        try:
            return self.subshells[(config.rstrip(), index)]
        except KeyError:
            raise JJKeyError('true-jj_keys.txt has no subshell %d of %r'
                             % (index, config))

    def configs(self):
        names = []
        for config, _ in self.subshells:
            if config not in names:
                names.append(config)
        return names


# --------------------------------------------------------------------------
# the alphabet
# --------------------------------------------------------------------------

def code_value(char):
    """The number a code character stands for, 1-based as LSJJCC counts.

    A J or K character reads 1 for 2J=0, 2 for 2J=1, and so on, so the
    doubled J is code_value(char) - 1.
    """
    if char == NO_CODE:
        raise JJKeyError('the label carries %r: the key ran out of codes '
                         'and the labels around it need not be unique'
                         % NO_CODE)
    try:
        return _VALUE[char]
    except KeyError:
        raise JJKeyError('%r is not one of the label codes' % char)


def layout(nopen):
    """The kinds of the eight label positions for `nopen` open subshells.

    Returns a list of (kind, subshell) pairs, one per used character, or
    raises where the label cannot name the state.  The forms are those
    LSJJKP documents at rcg11kd.f:16767-16797.
    """
    if nopen < 1:
        raise JJKeyError('a basis state with no open subshell has no label')
    if nopen <= 3:
        out = []
        for i in range(1, nopen + 1):
            out.append(('A', i))
            out.append(('J', i))
            if 1 < i < nopen:
                out.append(('K', i))
        return out
    if nopen == 4:
        # The Ki are dropped: A1 J1 A2 J2 A3 J3 A4 J4 fills the eight.
        return [(kind, i) for i in range(1, 5) for kind in ('A', 'J')]
    # Five or more: the Ji go as well, leaving the Ai alone, and only as
    # many of them as fit.
    return [('A', i) for i in range(1, min(nopen, 8) + 1)]


def decode(label, nopen):
    """Split a true jj label into its positions.

    `nopen` is the number of open subshells of the configuration, which
    the label itself does not carry: genealogy.txt or the configuration
    card supplies it.
    """
    plan = layout(nopen)
    text = label.rstrip()
    if len(text) < len(plan):
        raise JJKeyError('label %r is %d characters, but %d open subshells '
                         'need %d' % (label, len(text), nopen, len(plan)))
    out = []
    for (kind, sub), char in zip(plan, text):
        value = 0 if kind == 'A' else code_value(char) - 1
        out.append(Position(kind=kind, char=char, subshell=sub, value=value))
    return out


# --------------------------------------------------------------------------
# the file
# --------------------------------------------------------------------------

def _parts(text, where):
    """The partition inside the brackets of an entry.

    LSJJKP writes each half as the subshell letter, the sign, the
    occupation with no leading blank, then ' v' and one digit where the
    seniority is needed, then ' J=' and F4.1.  Halves are joined by
    ', ', and neither a comma nor a space can occur inside one.
    """
    out = []
    for piece in text.split(','):
        piece = piece.strip()
        if not piece:
            continue
        head, sep, jtext = piece.partition(' J=')
        if not sep:
            raise JJKeyError('no J in %r of %s' % (piece, where))
        seniority = 0
        if ' v' in head:
            head, _, vtext = head.partition(' v')
            seniority = int(vtext)
        try:
            out.append(Part(l=head[0], sign=head[1], nel=int(head[2:]),
                            seniority=seniority, j=float(jtext)))
        except ValueError:
            raise JJKeyError('%r of %s does not read as a subshell half'
                             % (piece, where))
    return tuple(out)


def read_jjkeys(path=files.JJKEYS, folder=None):
    """Read a true-jj_keys.txt into JJKeys."""
    if folder is not None:
        path = files.resolve(path, folder)
    lines = files.read_lines(path)

    subshells = {}
    current = None
    pending = None        # a code whose two-J entry is still being read

    for line in lines:
        if '   subshell' in line and line.startswith(' ') \
                and line.rstrip().endswith(')'):
            # format 20: /1X,3A6,'   subshell',I2,'   (',A1,A,')'
            wide = line.ljust(40)
            config = wide[1:19].rstrip()
            rest = wide[19:]
            head, _, tail = rest.partition('   subshell')
            index = int(tail[:2])
            inside = tail[2:].strip()
            if not (inside.startswith('(') and inside.endswith(')')):
                raise JJKeyError('the heading of subshell %d of %r does not '
                                 'end in its (l occupation)' % (index, config))
            inside = inside[1:-1]
            current = Subshell(config=config, index=index, l=inside[0],
                               nel=int(inside[1:]), codes={})
            subshells[(config, index)] = current
            pending = None
            continue

        if current is None:
            continue                        # the preamble

        body = line.rstrip()
        if not body.strip() or not body.startswith('     '):
            continue
        where = 'subshell %d of %r' % (current.index, current.config)

        if body[5:6] != ' ' and body[6:8] == ': ':
            # formats 23 and 24: a code, then either the partition or the
            # first of its two values of J
            code = body[5]
            rest = body[8:].strip()
            if rest.startswith('J='):
                jtext, _, tail = rest[2:].partition(': ')
                entry = Entry(code=code,
                              by_j=((float(jtext), _parts(_inside(tail, where),
                                                          where)),))
                pending = code
            else:
                entry = Entry(code=code, parts=_parts(_inside(rest, where),
                                                      where))
                pending = None
            current.codes[code] = entry
        elif body[8:10] == 'J=':
            # format 25: the second value of J of the pending code
            if pending is None:
                raise JJKeyError('a second J with no code before it in %s'
                                 % where)
            jtext, _, tail = body[10:].partition(': ')
            entry = current.codes[pending]
            current.codes[pending] = Entry(
                code=entry.code,
                by_j=entry.by_j + ((float(jtext),
                                    _parts(_inside(tail, where), where)),))
            pending = None
    return JJKeys(path, subshells)


def _inside(text, where):
    text = text.strip()
    if not (text.startswith('(') and text.endswith(')')):
        raise JJKeyError('%r of %s is not bracketed' % (text, where))
    return text[1:-1]


def _main(argv):
    """Summarize a true-jj_keys.txt: python -m toolbox.jjkeys [path]."""
    path = argv[1] if len(argv) > 1 else files.JJKEYS
    keys = read_jjkeys(path)
    print('%s: %d subshell keys over %d configurations'
          % (keys.path, len(keys), len(keys.configs())))
    for sub in keys:
        print('  %-18s subshell %d  (%s%d)  %d codes'
              % (sub.config, sub.index, sub.l, sub.nel, len(sub.codes)))
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(_main(sys.argv))
