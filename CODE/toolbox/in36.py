"""Reading IN36 the way RCN reads it.

RCN does not parse a configuration card by column.  Subroutine ANALYZ
(`FOR/RCN36K.F:1122-1238`) walks the card, splits it into blank-delimited
tokens, and classifies each token by where its orbital letter, decimal
point and minus sign fall - so `4f3`, `4f 3`, ` 4f3 ` and `4f 03` are all
the same subshell and the layout of a card is free.  Subroutine SETCFG
(`:1269-1405`) then builds the heaviest noble-gas core that fits the
electron count and applies the card's subshells to it, which is where the
principal quantum numbers come from.

Both are ported here as exactly as the two languages allow, including
ANALYZ's habit of resuming its scan from a cursor the classification
branch moved, because the alternative - a regular expression over the
same cards - is what the Perl did, and it is why WORK/U90 came out
labelled `1s` for the configuration `1s1 2s2`.

    confs = read_in36('IN36')
    confs[0].label        -> '2s'
    confs[0].shells       -> [Shell(1, 0, 2.0), Shell(2, 0, 1.0)]
    confs[0].parity       -> 0
    parities(confs)       -> {0} | {1} | {0, 1}
"""

from dataclasses import dataclass, field

from . import files
from .fortfmt import is_whole, read_f, read_i

# Orbital letters, RCN's own table (`FOR/RCN36K.F:1142-1143`).  j is not
# used, so the index of a letter is l.
LBCD = 'spdfghiklmnoqrtuvwxyz'

# Occupation strings recognised without a decimal point, in the order of
# the wwwnl table at `FOR/RCN36K.F:1264-1266`.  Entry j codes j-1 for
# j < 16 and j-16 from there on, which is what makes '7' and '07' equal.
_WWWNL = ['%-6s' % s for s in
          ('0', '1', '2', '3', '4', '5', '6', '7', '8', '9',
           '10', '11', '12', '13', '14',
           '00', '01', '02', '03', '04', '05', '06', '07', '08', '09')]

# Noble-gas electron counts and the number of subshells in each
# (`FOR/RCN36K.F:1319`, `:1263`).
_NE = (2, 10, 18, 36, 54, 86, 118)
_MORB = (1, 3, 5, 8, 11, 15, 19)

_CARD = 70          # width of inpcard, card columns 11-80
_MAX_ORB = 8        # nlbcd8/wwnl8 hold eight subshells per card


class In36Error(Exception):
    """A card RCN itself would have rejected."""


@dataclass(frozen=True)
class Shell:
    """One subshell of a configuration."""

    n: int
    l: int
    occ: float

    @property
    def full(self):
        """Occupation of the closed subshell, 2(2l+1)."""
        return 4 * self.l + 2

    @property
    def is_closed(self):
        return self.occ == self.full

    @property
    def is_open(self):
        return self.occ != 0 and self.occ != self.full

    @property
    def nl(self):
        """The subshell without its occupation, e.g. '4f'."""
        return '%d%s' % (self.n, LBCD[self.l])

    def __str__(self):
        if self.occ == 1:
            return self.nl
        if is_whole(self.occ):
            return '%s%d' % (self.nl, int(self.occ))
        return '%s%g' % (self.nl, self.occ)


@dataclass
class Configuration:
    """One configuration card, expanded as SETCFG expands it."""

    index: int                      # 1-based position among the cards
    spectrum: str                   # cols 11-16, e.g. 'Pr2+', 'B IV'
    label: str                      # cols 17-28, RCG's 3A6 name less the first
    iz: int
    nspect: int
    shells: list                    # full expansion, SETCFG's order
    card_shells: list               # only the subshells the card names
    kut: int = 0
    ee8: float = 0.0
    alfmin: float = 0.0
    alfmax: float = 0.0
    fractional: bool = False        # SETCFG's ifrac
    raw: str = ''                   # the card as read

    @property
    def noelec(self):
        """Number of electrons, `iz + 1 - nspect` (`FOR/RCN36K.F:476`)."""
        return self.iz + 1 - self.nspect

    @property
    def parity(self):
        """0 even, 1 odd."""
        return int(sum(s.l * s.occ for s in self.shells)) % 2

    @property
    def open_shells(self):
        return [s for s in self.shells if s.is_open]

    @property
    def closed_shells(self):
        return [s for s in self.shells if s.is_closed]

    def canonical(self):
        """The shells in display order: increasing n, then increasing l.

        SETCFG's order is the core followed by whatever the card added,
        which is a summation order, not a way to write a configuration.
        """
        return sorted(self.shells, key=lambda s: (s.n, s.l))

    def by_nl(self):
        """{'4f': Shell, ...} for matching against decks that lack n."""
        return {s.nl: s for s in self.shells}

    def __str__(self):
        return ' '.join(str(s) for s in self.canonical())


@dataclass
class In36:
    """A whole IN36 deck."""

    path: str = ''
    control: list = field(default_factory=list)     # the control card(s), raw
    configurations: list = field(default_factory=list)

    def __iter__(self):
        return iter(self.configurations)

    def __len__(self):
        return len(self.configurations)

    def by_label(self):
        """{label: Configuration}; a repeated label is an error."""
        out = {}
        for c in self.configurations:
            if c.label in out:
                raise In36Error(
                    'configuration label %r appears twice in %s (cards %d and %d)'
                    % (c.label, self.path, out[c.label].index, c.index))
            out[c.label] = c
        return out


# --------------------------------------------------------------------------
# ANALYZ
# --------------------------------------------------------------------------

def _put(buf, start, end, src):
    """Fortran substring assignment buf(start:end) = src, 1-based.

    The source is truncated to the range and the rest of the range is
    blank-filled, which is what makes ANALYZ's right-alignment work.
    """
    if end < start:
        return
    if start < 1 or end > len(buf):
        raise In36Error('card field %d:%d is outside the card' % (start, end))
    width = end - start + 1
    src = (src[:width]).ljust(width)
    buf[start - 1:end] = list(src)


def analyz(inpcard):
    """Port of `SUBROUTINE ANALYZ` (`FOR/RCN36K.F:1122-1238`).

    `inpcard` is card columns 11-80.  Returns
    (conf, nlbcd8, wwnl8, kut, ee8, alfmin, alfmax), where conf is the
    18-column 3A6 name, nlbcd8 holds up to eight 3-character subshell
    labels in RCN's `(i2,a1)` form (' 4f', '10s') and wwnl8 the matching
    6-character occupation fields.
    """
    card = list(inpcard.ljust(_CARD)[:_CARD])

    # The configuration name ends at three blanks or at column 18.
    i = 6
    nbl = 0
    while True:
        i += 1
        if i > _CARD:
            break
        if card[i - 1] != ' ':
            nbl = 0
            continue
        nbl += 1
        if nbl >= 3 or i >= 18:
            break
    i = min(i, 18)
    conf = ''.join(card[:i]).ljust(18)

    i3 = i
    nlbcd8 = ['   '] * _MAX_ORB
    wwnl8 = ['      '] * _MAX_ORB
    kut = 0
    ee8 = 0.0
    nalf = 0
    alfmin = alfmax = 0.0
    norb = 0

    while True:
        # 130: restart the scan from wherever the last branch left i3.
        i = i3
        i1 = i2 = i3 = idec = ineg = 0
        chtemp = [' '] * _CARD

        # 132: collect one blank-delimited token.
        done = False
        while True:
            i += 1
            if i > _CARD:
                done = True
                break
            c = card[i - 1]
            if c == ' ':
                if i3 != 0:
                    break
                continue
            if i1 == 0:
                i1 = i
            if i2 == 0:
                if 'A' <= c <= 'Z':             # ANALYZ folds in place
                    c = c.lower()
                    card[i - 1] = c
                if c in LBCD:
                    i2 = i
                if c == '.':
                    idec = i
                if c == '-':
                    ineg = i
            i3 = i
            if i >= _CARD and i1 != 0:
                break
        if done:
            break

        # 150: classify, in RCN's order (`:1198-1202`).
        if i2 != 0 and i1 >= i2 - 2:
            # 300: an orbital.  Right-align n before the letter in a
            # 9-column buffer, so ' 4f' and '3    ' come out of one read.
            i3 = min(i3 + 1, i2 + 6)
            _put(chtemp, 3 + i1 - i2, 9, ''.join(card[i1 - 1:i3]))
            if norb >= _MAX_ORB:
                raise In36Error(
                    'more than %d subshells on card: %r' % (_MAX_ORB, inpcard))
            nlbcd8[norb] = ''.join(chtemp[0:3])
            wwnl8[norb] = ''.join(chtemp[3:9])
            norb += 1
        elif norb == 0 and nalf == 0 and i3 > i1 and ineg == 0:
            # 200: alfmax, and alfmin if the field holds two numbers.
            if i2 > 0:
                i3 = i2 - 3
            _put(chtemp, i1 - i3 + 4, 4, ''.join(card[i1 - 1:i3]))
            alfmin = 0.01 * read_i(''.join(chtemp[0:2]))
            alfmax = 0.01 * read_i(''.join(chtemp[2:4]))
            nalf = 1
        elif idec > 0:
            ee8 = _read_ee8(card, i1, i3)
        elif ineg > 0 or i1 == i3:
            # 400: kut, the last two characters of the token.
            kut = read_i(''.join(card[i3 - 2:i3]))
        else:
            ee8 = _read_ee8(card, i1, i3)

    return conf, nlbcd8, wwnl8, kut, ee8, alfmin, alfmax


def _read_ee8(card, i1, i3):
    """500: the continuum energy, `format (f70.3)`."""
    return read_f(''.join(card[i1 - 1:i3]), 3)


# --------------------------------------------------------------------------
# SETCFG
# --------------------------------------------------------------------------

def _noble_core(noelec):
    """The core of `SETCFG` (`FOR/RCN36K.F:1326-1355`).

    Returns [(n, l, occ)] for the heaviest noble gas with no more than
    `noelec` electrons.  The m=m-1 corrections drop the f and g
    subshells that the naive n,l enumeration would insert ahead of a
    noble gas that does not have them.
    """
    necore = 0
    mcore = 0
    for k in range(7):
        if _NE[k] > noelec:
            break
        necore = _NE[k]
        mcore = _MORB[k]

    orb = []
    for n in range(1, 8):
        for lp1 in range(1, n + 1):
            if mcore == len(orb):
                return orb
            l = lp1 - 1
            orb.append([n, l, float(4 * lp1 - 2)])
            if ((necore == 54 and n == 4 and l == 3)
                    or (necore == 86 and n == 5 and l in (3, 4))
                    or (necore == 118 and n == 5 and l == 4)
                    or (necore == 118 and n == 6 and l in (3, 4, 5))):
                orb.pop()
    return orb


def _occupation(wwnl):
    """Decode one 6-character occupation field (`FOR/RCN36K.F:1362-1373`).

    Returns (occupation, fractional).  A field that matches the wwwnl
    table is an integer count; anything else is 1 electron unless it
    looks like a written-out number, in which case it is read as F6.3.
    """
    for j, entry in enumerate(_WWWNL, start=1):
        if wwnl == entry:
            return float(j - 1 if j < 16 else j - 16), False
    occ = 1.0
    if wwnl[2] != ' ' or wwnl[1] == '.' or wwnl[0] == '.':
        occ = read_f(wwnl, 3)
    return occ, not is_whole(occ)


def _split_label(label):
    """' 4f' -> (4, 3); the inverse of `write (nlbcd(m),22) n, lbcd(l+1)`."""
    n = read_i(label[0:2])
    letter = label[2]
    if letter not in LBCD:
        raise In36Error('orbital label %r is not in RCN\'s table' % label)
    l = LBCD.index(letter)
    if l + 1 > n:
        raise In36Error('orbital %r has l >= n' % label.strip())
    return n, l


def setcfg(noelec, nlbcd8, wwnl8, n1scf=1):
    """Port of `SUBROUTINE SETCFG` (`FOR/RCN36K.F:1269-1405`).

    Returns (shells, card_shells, fractional).  `shells` is in SETCFG's
    own order - the noble-gas core first, then whatever the card added.
    """
    orb = _noble_core(noelec)
    labels = ['%2d%s' % (n, LBCD[l]) for n, l, _ in orb]
    card_shells = []
    fractional = False

    for i in range(_MAX_ORB):
        if nlbcd8[i][0:3] == '   ':
            continue
        occ, frac = _occupation(wwnl8[i])
        fractional = fractional or frac
        n, l = _split_label(nlbcd8[i])
        card_shells.append(Shell(n, l, occ))

        for j, lab in enumerate(labels):
            if nlbcd8[i] != lab:
                continue
            if occ != 0.0 or j + 1 < n1scf:
                orb[j][2] = occ                 # 240
            else:
                del orb[j]                      # a zero deletes the subshell
                del labels[j]
            break
        else:
            # 237: a subshell outside the core.  Zero means absent.
            if occ == 0.0 and len(orb) + 1 >= n1scf:
                continue
            orb.append([n, l, occ])
            labels.append(nlbcd8[i])

    return [Shell(n, l, occ) for n, l, occ in orb], card_shells, fractional


# --------------------------------------------------------------------------
# the deck
# --------------------------------------------------------------------------

def read_in36(path='IN36', folder=None):
    """Read an IN36 deck; returns an In36.

    Card layout is `format (i1,i4,i3,i2,a70)` (`FOR/RCN36K.F:437-438`):
    nhftrm, iz, ibb1, nspect, then columns 11-80 for ANALYZ.  Reading
    stops at the first card with iz < 0; iz = 0 starts a further control
    card, and a card whose columns 11-16 are blank tells RCN to skip
    completed runs on TAPE2 rather than describing a configuration.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    lines = files.read_lines(path, width=80)

    deck = In36(path=str(path))
    i = 0
    expect_control = True
    index = 0
    while i < len(lines):
        line = lines[i]
        i += 1
        if expect_control:
            deck.control.append(line)
            expect_control = False
            continue
        if not line.strip():
            continue

        nhftrm = read_i(line[0:1])
        iz = read_i(line[1:5])
        ibb1 = read_i(line[5:8])
        nspect = read_i(line[8:10])
        inpcard = line[10:80]

        if iz < 0:
            break
        if iz == 0:
            expect_control = True
            continue
        if inpcard[0:6] == '      ':
            continue                            # 214: TAPE2 skip card

        # `n1scf` is the frozen-core guard set from ibb1 (`:439-440`); it
        # stops a zero occupation from deleting a frozen core subshell.
        # RCN also requires nconft > 0 there, which only arises when a run
        # is continued from TAPE2; no deck in WORK sets ibb1 at all.
        n1scf = max(1, ibb1 + 1) if 0 < ibb1 <= 20 else 1

        conf, nlbcd8, wwnl8, kut, ee8, alfmin, alfmax = analyz(inpcard)
        if nspect == 0:
            nspect = 1
        noelec = iz + 1 - nspect
        shells, card_shells, frac = setcfg(noelec, nlbcd8, wwnl8, n1scf)

        index += 1
        deck.configurations.append(Configuration(
            index=index,
            spectrum=conf[0:6].strip(),
            label=conf[6:18].strip(),
            iz=iz, nspect=nspect,
            shells=shells, card_shells=card_shells,
            kut=kut, ee8=ee8, alfmin=alfmin, alfmax=alfmax,
            fractional=frac, raw=line.rstrip()))

        # HF term cards follow their configuration card (`:487-494`).
        for _ in range(nhftrm):
            if i >= len(lines):
                break
            nfg = read_i(lines[i][10:15])
            i += 1
            if nfg > 3 and i < len(lines):
                i += 1

    return deck


def parities(configurations):
    """{0}, {1} or {0, 1} - one parity or two."""
    return {c.parity for c in configurations}


def common_closed_shells(configurations):
    """The closed subshells every configuration shares, in canonical order.

    These are the ones that may be dropped from a printed label without
    making the electron count differ from one level to the next.
    """
    if not configurations:
        return []
    shared = None
    for c in configurations:
        here = {(s.n, s.l) for s in c.shells if s.is_closed}
        shared = here if shared is None else (shared & here)
    return [Shell(n, l, 4 * l + 2) for n, l in sorted(shared)]


def _main(argv):
    """python -m toolbox.in36 [IN36] - show what RCN would make of a deck."""
    deck = read_in36(argv[1] if len(argv) > 1 else 'IN36')
    shared = common_closed_shells(deck.configurations)
    print('%s: %d configurations, %s'
          % (deck.path, len(deck),
             {frozenset({0}): 'even parity only',
              frozenset({1}): 'odd parity only'}.get(
                  frozenset(parities(deck.configurations)), 'both parities')))
    print('closed in every configuration: %s'
          % (' '.join(s.nl for s in shared) or '(none)'))
    for c in deck:
        print('%3d  %-8s %s  %s'
              % (c.index, c.label, '*' if c.parity else ' ', c))
    return 0


if __name__ == '__main__':
    import sys
    sys.exit(_main(sys.argv))
