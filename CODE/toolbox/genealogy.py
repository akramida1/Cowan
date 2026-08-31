"""Reading genealogy.txt, the coupling chain of every basis state.

RCG writes this file when columns 83-84 of the ING11 main card hold
IPGEN=2.  It is the only place the package states, in full, what a basis
state is: an eight-character label has no room for the subshell terms,
the subshell J's and the running quantum numbers that couple them, so an
eigenvector composition in OUTG11 cannot be read back to the states it
names without this table.

One table covers one J block of one configuration set, one line per
basis state, numbered as the eigenvector components are.  A field of a
line is one open subshell, taken in the order of the configuration;
empty and closed subshells are left out.  Three tables exist per block:

    LS       (term index, 2S+1, L, Cowan's alpha) then the running
             2S+1 and L after coupling that subshell in
    J1J2     (term index, 2S+1, L, J of the subshell) then the running
             J.  Cowan's jj, in which a subshell is still an LS term.
    true jj  (j- part, j+ part) then 2J of the subshell and 2K.  Every
             J in this table is doubled, as it is in the label key.

The writer is SUBROUTINE GENPRT, FOR/rcg11kd.f:16966-17275, and its
FORMAT statements are the specification followed here:

    11  block heading                        rcg11kd.f:17114
    13  a configuration of the block         rcg11kd.f:17130
    31  an LS field, 15 columns              rcg11kd.f:17187
    32  a J1J2 field, 20 columns             rcg11kd.f:17194
    42/43 a true jj part inside its field    rcg11kd.f:17246, :17250
    44  a true jj field, 30 columns          rcg11kd.f:17256
    45  a state line                         rcg11kd.f:17268
    46  its continuation                     rcg11kd.f:17269

This module reads; it does not interpret.  Turning a chain into
spectroscopic notation belongs to asd.py.
"""

from dataclasses import dataclass, field as _field

from . import files

# Table names as this module reports them, and the field width each one
# was written with (NW in GENPRT).
LS = 'LS'
J1J2 = 'J1J2'
TRUEJJ = 'true jj'

_WIDTH = {LS: 15, J1J2: 20, TRUEJJ: 30}

# The heading that opens each table (rcg11kd.f:17138, :17145, :17152).
_TABLE_OF = {
    ' LS basis': LS,
    ' J1J2 basis': J1J2,
    ' true jj basis': TRUEJJ,
}

# A state line carries its number in columns 1-7; a continuation line is
# blank there (format 46 is 22X).
_LEAD = 22


class GenealogyError(Exception):
    """genealogy.txt did not read as SUBROUTINE GENPRT writes it."""


@dataclass(frozen=True)
class LSField:
    """One open subshell of an LS basis state (format 31)."""
    term: int            # index into the term list of that subshell
    mult: int            # 2S+1 of the subshell term
    l: str               # its L, as a letter
    alpha: str           # Cowan's alpha, blank where he prints none
    run_mult: int        # 2S+1 after coupling this subshell in
    run_l: str           # and its L


@dataclass(frozen=True)
class J1J2Field:
    """One open subshell of a J1J2 basis state (format 32)."""
    term: int
    mult: int
    l: str
    j: float             # J of the subshell
    run_j: float         # J after coupling it in


@dataclass(frozen=True)
class JJPart:
    """The j- or the j+ half of a subshell in true jj (formats 42/43)."""
    l: str               # subshell letter, lower case
    sign: str            # '-' for j = l-1/2, '+' for j = l+1/2
    nel: int             # electrons in this half
    seniority: int       # 0 where GENPRT prints none, i.e. where J alone
    j2: int              # 2J of this half


@dataclass(frozen=True)
class JJField:
    """One open subshell of a true jj basis state (format 44).

    j2 and k2 are doubled, as everything in that table is.  They repeat
    the J1J2 line above: the rotation mixes only the subshell terms of
    columns that already share the whole chain of J.
    """
    parts: tuple         # of JJPart, j- first; a half with no electrons
                         # is left out, so this holds one or two
    j2: int              # 2J of the subshell
    k2: int              # 2K, K being the running J


@dataclass(frozen=True)
class State:
    """One line of one table: one basis state."""
    number: int          # its place in the eigenvector, 1-based
    cfg: int             # configuration number WITHIN ITS OWN SET
    label: str           # the eight characters OUTG11 prints
    fields: tuple        # one entry per open subshell, in card order


@dataclass
class Block:
    """One J block of one configuration set."""
    j: float
    kset: int            # configuration set, 1 or 2
    nstates: int         # basis states, from the heading
    configurations: dict = _field(default_factory=dict)
    # {cfg number: ((letter, occupation), ... one per subshell slot)}
    tables: dict = _field(default_factory=dict)
    # {LS|J1J2|TRUEJJ: [State, ...]}

    @property
    def key(self):
        return (self.kset, self.j)

    def state(self, table, number):
        """The basis state numbered `number` in the named table."""
        rows = self.tables.get(table)
        if rows is None:
            raise GenealogyError(
                'no %s table for J=%.1f of set %d' % (table, self.j, self.kset))
        if not 1 <= number <= len(rows):
            raise GenealogyError(
                'state %d is outside the %s table of J=%.1f of set %d, '
                'which holds %d' % (number, table, self.j, self.kset, len(rows)))
        row = rows[number - 1]
        if row.number != number:
            raise GenealogyError(
                'the %s table of J=%.1f of set %d is not in order at state %d'
                % (table, self.j, self.kset, number))
        return row


class Genealogy:
    """Every block of one genealogy.txt, keyed by (set, J)."""

    def __init__(self, path, blocks):
        self.path = path
        self.blocks = blocks

    def __iter__(self):
        return iter(self.blocks)

    def __len__(self):
        return len(self.blocks)

    def block(self, kset, j):
        for b in self.blocks:
            if b.kset == kset and abs(b.j - j) < 1e-6:
                return b
        raise GenealogyError(
            'genealogy.txt has no J=%.1f block for configuration set %d'
            % (j, kset))

    def tables(self):
        """The table names present, in the order GENPRT writes them."""
        seen = []
        for b in self.blocks:
            for name in (LS, J1J2, TRUEJJ):
                if name in b.tables and name not in seen:
                    seen.append(name)
        return seen


# --------------------------------------------------------------------------
# field decoding
# --------------------------------------------------------------------------

def _int(text, what, j, kset):
    text = text.strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        raise GenealogyError('%s is not a number in J=%.1f of set %d: %r'
                             % (what, j, kset, text))


def _float(text, what, j, kset):
    text = text.strip()
    if not text:
        return 0.0
    try:
        return float(text)
    except ValueError:
        raise GenealogyError('%s is not a number in J=%.1f of set %d: %r'
                             % (what, j, kset, text))


def _ls_field(fld, j, kset):
    # format 31: ('  ,I4,I3,A1,A1,')',I3,A1
    return LSField(term=_int(fld[1:5], 'a term index', j, kset),
                   mult=_int(fld[5:8], 'a multiplicity', j, kset),
                   l=fld[8],
                   alpha=fld[9].strip(),
                   run_mult=_int(fld[11:14], 'a running multiplicity', j, kset),
                   run_l=fld[14])


def _j1j2_field(fld, j, kset):
    # format 32: ('  ,I4,I3,A1,F5.1,')',F5.1
    return J1J2Field(term=_int(fld[1:5], 'a term index', j, kset),
                     mult=_int(fld[5:8], 'a multiplicity', j, kset),
                     l=fld[8],
                     j=_float(fld[9:14], 'a subshell J', j, kset),
                     run_j=_float(fld[15:20], 'a running J', j, kset))


def _jj_parts(text, j, kset):
    """The (j- part, j+ part) of format 44's A20, formats 42 and 43.

    A part is written 2A1,I2 for the subshell and its occupation, then
    'v' and one digit where the seniority is needed, then I3 for 2J.  A
    part with no electrons is left out, and the two are joined by a
    comma with no space, so the comma is the separator and nothing else
    in the text can be one.
    """
    parts = []
    for piece in text.split(','):
        piece = piece.rstrip()
        if not piece.strip():
            continue
        seniority = 0
        if 'v' in piece:
            head, _, tail = piece.partition('v')
            seniority = _int(tail[:1], 'a seniority', j, kset)
            j2 = _int(tail[1:], 'a subshell 2J', j, kset)
        else:
            head, j2 = piece[:4], _int(piece[4:], 'a subshell 2J', j, kset)
        parts.append(JJPart(l=head[0], sign=head[1],
                            nel=_int(head[2:4], 'an occupation', j, kset),
                            seniority=seniority, j2=j2))
    return tuple(parts)


def _jj_field(fld, j, kset):
    # format 44: ('  ,A20,')',I4,I4
    return JJField(parts=_jj_parts(fld[1:21], j, kset),
                   j2=_int(fld[22:26], 'a subshell 2J', j, kset),
                   k2=_int(fld[26:30], 'a running 2K', j, kset))


_DECODE = {LS: _ls_field, J1J2: _j1j2_field, TRUEJJ: _jj_field}


def _fields(line, table, j, kset):
    """The fields of one state line or continuation, left to right.

    GENPRT writes them as 8(2X,A) after a 22-column lead, and builds
    each one into a string first so that a line of fewer fields does not
    end on the bracket of a field that never came.  A blank chunk is
    therefore the end of the line and not a missing subshell.
    """
    width = _WIDTH[table]
    decode = _DECODE[table]
    out = []
    pos = _LEAD
    while pos + 2 + width <= len(line):
        chunk = line[pos + 2:pos + 2 + width]
        if not chunk.strip():
            break
        out.append(decode(chunk, j, kset))
        pos += 2 + width
    return out


# --------------------------------------------------------------------------
# the file
# --------------------------------------------------------------------------

def read_genealogy(path=files.GENEALOGY, folder=None):
    """Read a genealogy.txt into Blocks.

    A block's heading is written once per GENPRT call, so a J block with
    a true jj table appears twice - once for the LS and J1J2 tables and
    once, a few lines later, for the rotated one.  The two are merged
    here on (set, J), which is unique within a run.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    lines = files.read_lines(path)

    blocks = []
    index = {}
    block = None
    table = None
    reading_configurations = False

    for raw in lines:
        line = raw.rstrip()
        if line.startswith(' J ='):
            # format 11: ' J =',F5.1,'    configuration set',I2,
            #            '    basis states:',I7
            wide = line.ljust(56)
            j = _float(wide[4:9], 'a block J', 0.0, 0)
            kset = _int(wide[30:32], 'a configuration set', j, 0)
            key = (kset, j)
            block = index.get(key)
            if block is None:
                block = Block(j=j, kset=kset,
                              nstates=_int(wide[49:56], 'a state count', j, kset))
                index[key] = block
                blocks.append(block)
            table = None
            reading_configurations = False
            continue

        if block is None:
            continue                       # the preamble

        if line.startswith('    configurations:'):
            reading_configurations = True
            continue

        head = line[:14]
        for prefix, name in _TABLE_OF.items():
            if head.startswith(prefix):
                table = name
                block.tables.setdefault(table, [])
                reading_configurations = False
                break
        else:
            if line.startswith('  state '):
                continue                   # the column heading of a table
            if not line.strip():
                continue
            if reading_configurations:
                # format 13: I7,4X,8(A1,I2,2X)
                wide = line.ljust(51)
                cfg = _int(wide[:7], 'a configuration number', block.j,
                           block.kset)
                shells = []
                for i in range(8):
                    at = 11 + 5 * i
                    letter = wide[at:at + 1]
                    if not letter.strip():
                        break
                    shells.append(
                        (letter, _int(wide[at + 1:at + 3], 'an occupation',
                                      block.j, block.kset)))
                block.configurations[cfg] = tuple(shells)
                continue
            if table is None:
                continue
            rows = block.tables[table]
            wide = line.ljust(_LEAD + 2)
            if wide[:_LEAD].strip():
                # format 45: I7,I5,2X,A8,8(2X,A)
                number = _int(wide[:7], 'a state number', block.j, block.kset)
                cfg = _int(wide[7:12], 'a configuration number', block.j,
                           block.kset)
                rows.append(State(number=number, cfg=cfg,
                                  label=wide[14:22],
                                  fields=tuple(_fields(line, table, block.j,
                                                       block.kset))))
            else:
                # format 46: the same state, more subshells
                if not rows:
                    raise GenealogyError(
                        'a continuation line precedes any state in the %s '
                        'table of J=%.1f of set %d' % (table, block.j,
                                                       block.kset))
                last = rows[-1]
                rows[-1] = State(number=last.number, cfg=last.cfg,
                                 label=last.label,
                                 fields=last.fields + tuple(
                                     _fields(line, table, block.j, block.kset)))

    for b in blocks:
        for name, rows in b.tables.items():
            if len(rows) != b.nstates:
                raise GenealogyError(
                    'the %s table of J=%.1f of set %d holds %d states, but '
                    'its heading says %d' % (name, b.j, b.kset, len(rows),
                                             b.nstates))
    return Genealogy(path, blocks)


def _main(argv):
    """Summarize a genealogy.txt: python -m toolbox.genealogy [path]."""
    path = argv[1] if len(argv) > 1 else files.GENEALOGY
    gen = read_genealogy(path)
    print('%s: %d blocks, tables %s'
          % (gen.path, len(gen), ', '.join(gen.tables())))
    for b in gen:
        counts = ', '.join('%s %d' % (n, len(b.tables[n]))
                           for n in (LS, J1J2, TRUEJJ) if n in b.tables)
        print('  set %d  J=%4.1f  %5d states  %2d configurations  [%s]'
              % (b.kset, b.j, b.nstates, len(b.configurations), counts))
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(_main(sys.argv))
