"""The ING11 deck: RCG's input, and the only place `n` and the deck's
own summation order can be read together.

Three things here are needed by everything downstream:

  * **the main card**, because it says how many configurations of each
    parity there are, which couplings were asked for (`KCPLD`), and
    whether the genealogy was printed (`IPGEN`);
  * **the configuration cards**, because they carry the 18-character
    name that every eigenvector row and every parameter card refers to,
    and the *slot order*, which is the order RCG sums the subshells in;
  * **the parameter cards**, because their names are positional --
    `F2(11)`, `ZETA 3` -- and become readable only once a slot is tied
    to an `nl`.

A configuration card names a subshell by its L symbol and occupation
only (`F 3  P 6  D 0 ...`); `n` is in IN36 and nowhere else.  Tying the
two together is `match_shells()` below, and getting it right is what
fixes the labels the Perl gets wrong on a reordered deck.

Every FORMAT mirrored here is named by its `FOR/rcg11kd.f` line.  Where
this reader is stricter than the Fortran it says so; it never guesses.
"""

from dataclasses import dataclass, field

from . import files
from .fortfmt import read_f, read_i, write_f_wide, write_param
from .in36 import LBCD

# RCG's term-symbol table, `FOR/rcg11kd.f:1081-1083`.  A configuration
# card names a subshell with one of these, so the index is l.
LSYMB = 'SPDFGHIKLMNOQRTUVWXYZABC'

# The main card's DMIN field, columns 61-65, read as `E5.1` (`:899`),
# and the KCPL=3 card's wavelength window, columns 11-15 and 16-20,
# read as `2E5.0` since 2026 (`:727`).
_DMIN = slice(60, 65)
_FLBDN = slice(10, 15)
_FLBDX = slice(15, 20)

# What RCG supplies for a blank or impossible window, `:774-775`.
FLBDN_DEFAULT, FLBDX_DEFAULT = 0.001, 500000.0

# Subshell slots on a configuration card: 8 x (A1,I2,2X), `:3334`.
NOSUBC_MAX = 8

# The couplings, by their KCPLD digit; RCG's own table is COUPL at
# `FOR/rcg11kd.f:1079-1080`.  KCPLD(I)=0 asks for coupling I, for I in
# 1..7 -- zero means yes because a blank field is what Cowan's decks
# punch and LS is what they want (`:786-790`).  Column 8 is inverted:
# KCPLD(8)=8 asks for true jj, precisely so that a blank does not.
KCPLD_COUPLINGS = {1: 'LS', 2: 'JJ', 3: 'JJJK', 4: 'LSLK', 5: 'LSJK',
                   6: 'LSJLKS', 7: 'LSJLSJ'}
TRUE_JJ = 'TRUEjj'


class Ing11Error(Exception):
    """A deck that cannot be read as RCG would read it."""


# --------------------------------------------------------------------------
# the cards
# --------------------------------------------------------------------------

@dataclass
class MainCard:
    """The parameter line, `FOR/rcg11kd.f:619-698`.

    Only the fields any tool in this package has a use for are named;
    the rest of the card is kept in `raw` so that a writer can put back
    what it did not understand.
    """

    kcpl: int = 0
    nck: tuple = (0, 0)
    nocset: int = 0
    nsconf: tuple = ((0, 0, 0), (0, 0, 0))   # [set][0..2], 1-based in Fortran
    iabg: int = 0
    iv: int = 0
    nlsmax: int = 0                 # as RCG holds it: blank reads 999
    nlspch: int = 0                 # what was actually punched there
    nlt11: int = 0
    nevmax: int = 0
    kcpld: tuple = (0,) * 9
    ielpun: int = 0
    imag: int = 0
    iquad: int = 0
    uenrgy: float = 0.0
    ilncuv: int = 0
    iplev: int = 0
    icpc: int = 0
    icfc: int = 0
    idip: int = 0
    isenior: int = 0
    ipgen: int = 0
    dmin: float = 0.0
    raw: str = ''

    @property
    def nsets(self):
        """1 if the run has one parity, 2 if it has both."""
        return 2 if self.nsconf[1][1] else 1

    def nconfigurations(self, k):
        """Configuration cards in parity set k (k = 1 or 2), NSCONF(2,K)."""
        return self.nsconf[k - 1][1]

    @property
    def true_jj(self):
        """KCPLD(8)=8 -- the added true jj sections (`:786-812`)."""
        return self.kcpld[7] == 8

    @property
    def couplings(self):
        """The couplings this deck asks RCG to print, in RCG's order.

        Couplings 3 to 7 are gated by NLSMAX as well as by their KCPLD
        digit: a J block with more than NLSMAX basis states gets none of
        them (`:1739`, `:9619`, `:14993`).  Since 2026 columns 23-25 are
        a permission rather than a requirement -- blank reads as 999 and
        permits them, and only a punched 0 forbids them
        (`rcg11kd.f:897-910`, `:955-978`) -- which is what `nlsmax`
        already holds, blank having been turned into 999 on the way in.
        This list is therefore what the deck *asks* for; what a
        particular block actually carries is a question for `outg11.py`.
        """
        out = []
        for i in sorted(KCPLD_COUPLINGS):
            if self.kcpld[i - 1] == 0 and (i <= 2 or self.nlsmax):
                out.append(KCPLD_COUPLINGS[i])
        if self.true_jj:
            out.append(TRUE_JJ)
        return out

    @property
    def nlsmax_effective(self):
        """NLSMAX as RCG uses it: 999 means no limit (`:977-978`).

        A blank field is already 999 here, so it lifts the limit too;
        a punched 0 is a limit of 0, which is what forbids couplings
        3-7 outright.
        """
        if self.nlsmax == 999 or (self.nocset and self.nlsmax):
            return None                     # KMX, i.e. every block qualifies
        return self.nlsmax


@dataclass
class ConfigCard:
    """One configuration card, `FOR/rcg11kd.f:3316-3334`.

    `slots` is the card's own order -- the order RCG couples the
    subshells in, which a reordered deck deliberately permutes.  It is
    not display order and must never be sorted in place.
    """

    kset: int                       # parity set, 1 or 2
    index: int                      # position within the set, 1-based
    slots: list                     # [(L letter, occupation), ...], 8 long
    name: str                       # cols 41-58, the 3A6 RCG prints
    eav: float = 0.0                # cols 59-71, average energy
    ekin: float = 0.0               # cols 72-80, kinetic energy
    nlt: int = 0                    # cols 83-85, per-config term cutoff
    raw: str = ''

    @property
    def spectrum(self):
        """The first A6 of the name: 'Pr2+  ', 'U89+  '."""
        return self.name[:6].strip()

    @property
    def label(self):
        """The rest of the name: '4f3', '2s2p', 'f26p'."""
        return self.name[6:].strip()

    @property
    def occupied(self):
        """The slots that hold electrons, in card order: [(index, l, occ)]."""
        return [(i + 1, _l_of(sym), occ)
                for i, (sym, occ) in enumerate(self.slots) if occ]

    @property
    def nelectrons(self):
        """Electrons named on the card.

        This is *not* the electron count of the spectrum: RCN2 strips the
        closed core every configuration shares, so W70's `s22p` shows 3.
        The real count comes from `in36.py`.
        """
        return sum(occ for _, occ in self.slots)

    @property
    def parity(self):
        """0 even, 1 odd -- RCG's IPRITY, as a parity bit (`:3357`)."""
        return sum(_l_of(sym) * occ for sym, occ in self.slots) % 2

    def __str__(self):
        return '%s (set %d, #%d)' % (self.name.strip(), self.kset, self.index)


@dataclass
class ParameterCard:
    """A single-configuration parameter card, `FOR/rcg11kd.f:10787-10790`.

    `FORMAT (2A6,A7,1X,F10.5,4(F9.4,I1),A2,4I2)` -- an 18-character name
    matching a configuration card's, EAV, then four value/group pairs,
    continued 7 at a time by `FORMAT (7(F9.4,I1))` (`:10833`).
    """

    name: str                       # cols 1-18, matching a configuration card
    np: int = 0                     # cols 19-20, RCN2's parameter count
    values: list = field(default_factory=list)      # EAV first
    groups: list = field(default_factory=list)      # LSF group per value
    kind: str = ''                  # cols 71-72, 'hr' and friends
    factors: tuple = (0, 0, 0, 0)   # cols 73-80, the 4I2 scale percentages
    raw: list = field(default_factory=list)
    fields: list = field(default_factory=list)   # each value's ten columns
    edited: bool = False            # whether `raw` still says what the card is
    trailing: list = field(default_factory=list)    # raw lines that follow it

    @property
    def label(self):
        return self.name[6:].strip()


@dataclass
class CiCard:
    """A configuration-interaction parameter card, `FOR/rcg11kd.f:10900`.

    `FORMAT (2A6,A7,1X,5(F9.4,I1),A2,4I2)` -- the same header, but five
    values on the first line instead of EAV plus four, because a CI card
    has no average energy.
    """

    name: str                       # 'label1   -label2  ', as punched
    np: int = 0
    values: list = field(default_factory=list)
    groups: list = field(default_factory=list)
    kind: str = ''
    factors: tuple = (0, 0, 0, 0)
    raw: list = field(default_factory=list)
    fields: list = field(default_factory=list)   # each value's ten columns
    edited: bool = False            # whether `raw` still says what the card is
    trailing: list = field(default_factory=list)    # raw lines that follow it

    @property
    def labels(self):
        """The two configuration labels the card couples."""
        left, _, right = self.name.partition('-')
        return left.strip(), right.strip()


@dataclass
class Ing11:
    """A whole deck: the prelude, the main card, and one block per parity."""

    path: str = ''
    prelude: list = field(default_factory=list)     # KCPL=3/9/>10 cards, raw
    main: MainCard = None
    configurations: list = field(default_factory=list)   # every card, in order
    parameters: list = field(default_factory=list)       # per set: [cards]
    ci: list = field(default_factory=list)               # per set: [cards]
    trailer: list = field(default_factory=list)          # radial-integral and
                                                         # sentinel cards, raw

    def set(self, k):
        """The configuration cards of parity set k, in card order."""
        return [c for c in self.configurations if c.kset == k]

    @property
    def emax(self):
        """(E0MAX(1), E0MAX(2)) of the KCPL=0 rescale card, or (0, 0).

        RCG reads that card at `FOR/rcg11kd.f:673`:

            FORMAT (20X,F10.5,10I1,10X,F10.5,F5.2,2F10.5)

        so E0MAX(1) is columns 66-75 and E0MAX(2) columns 76-85.  They
        are the highest energy RCG will print a level at, and are what
        tells a level identification which levels it is allowed not to
        find (`levels.identify`, `emax`).

        `calc_TP_M1E2.pl:194` reads the same two numbers from columns
        61-70 and 71-80, five columns to the left.  The two readings
        agree for any card whose values are right-justified in ten-column
        fields, which is every deck in WORK/, because the fields the
        Perl reads then hold the same digits followed by blanks; they
        would not agree for a value punched hard against the left of its
        field.  The FORMAT is the specification, so it is what is used.
        """
        for line in self.prelude:
            if read_i(line[:5]) == 0 and line[:5].strip() in ('0', ''):
                padded = line.ljust(85)
                return read_f(padded[65:75], 5), read_f(padded[75:85], 5)
        return 0.0, 0.0

    @property
    def window(self):
        """(FLBDN, FLBDX), the printed wavelength window in angstrom.

        The optional first control card, KCPL=3, carries them in columns
        11-15 and 16-20; RCG reads it at `FOR/rcg11kd.f:725` and
        supplies 0.001 and 500000.0 for a blank or impossible field at
        label 95, which is what is returned for a deck with no such
        card.  Since 2026 the fields are `E5.0` and read as angstrom.
        """
        line = self.control(3)
        flbdn = flbdx = 0.0
        if line is not None:
            padded = line.ljust(20)
            flbdn = read_f(padded[_FLBDN])
            flbdx = read_f(padded[_FLBDX])
        if flbdn <= 0.0:
            flbdn = FLBDN_DEFAULT
        if flbdx <= flbdn:
            flbdx = FLBDX_DEFAULT
        return flbdn, flbdx

    def control(self, kcpl):
        """The prelude card with this KCPL, or None."""
        for line in self.prelude:
            if line[:5].strip() and read_i(line[:5]) == kcpl:
                return line
        return None

    def by_label(self):
        """{label: ConfigCard}; a repeated label is an error, not a pick."""
        out = {}
        for card in self.configurations:
            if card.label in out:
                raise Ing11Error('configuration label %r is on two cards of %s'
                                 % (card.label, self.path))
            out[card.label] = card
        return out

    def __iter__(self):
        return iter(self.configurations)

    def __len__(self):
        return len(self.configurations)


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def _l_of(symbol):
    """l behind a configuration card's L symbol; RCG's loop at `:3338`."""
    i = LSYMB.find(symbol.upper())
    if i < 0:
        raise Ing11Error('%r is not one of RCG\'s subshell symbols' % (symbol,))
    return i


def read_main(card):
    """Decode the main parameter line, `FORMAT` at `FOR/rcg11kd.f:696`.

    A configuration count of 100 or more does not fit its 2-digit field,
    so RCN2 punches `**` there and all four counts are read instead from
    columns 93-104 as 4I3 (`:701-710`).
    """
    c = card.ljust(104)
    chconf1, chconf2, chconfs = c[11:15], c[16:20], c[92:104]
    stars = any('**' in f for f in (chconf1[:2], chconf1[2:],
                                    chconf2[:2], chconf2[2:]))
    if stars:
        n = [read_i(chconfs[i:i + 3]) for i in range(0, 12, 3)]
        set1 = (read_i(c[10:11]), n[0], n[1])
        set2 = (read_i(c[15:16]), n[2], n[3])
    else:
        set1 = (read_i(c[10:11]), read_i(chconf1[:2]), read_i(chconf1[2:4]))
        set2 = (read_i(c[15:16]), read_i(chconf2[:2]), read_i(chconf2[2:4]))
    return MainCard(
        kcpl=read_i(c[0:5]),
        nck=(read_i(c[5:6]), read_i(c[6:8])),
        nocset=read_i(c[8:10]),
        nsconf=(set1, set2),
        iabg=read_i(c[20:21]),
        iv=read_i(c[21:22]),
        nlsmax=_nlsmax(c[22:25]),
        nlspch=read_i(c[22:25]),
        nlt11=read_i(c[25:27]),
        nevmax=read_i(c[27:30]),
        kcpld=tuple(read_i(c[30 + i:31 + i]) for i in range(9)),
        ielpun=read_i(c[39:40]),
        imag=read_i(c[48:49]),
        iquad=read_i(c[49:50]),
        uenrgy=read_f(c[50:60], 5),
        ilncuv=read_i(c[65:66]),
        iplev=read_i(c[66:67]),
        icpc=read_i(c[67:68]),
        icfc=read_i(c[68:69]),
        idip=read_i(c[69:70]),
        dmin=read_f(c[60:65], 1),
        isenior=read_i(c[80:82]),
        ipgen=read_i(c[82:84]),
        raw=card)


def _nlsmax(field):
    """NLSMAX as RCG holds it after reading columns 23-25.

    RCG reads the field as A3 and a blank one becomes 999, "no limit"
    (`rcg11kd.f:897-910`).  Only an explicitly punched 0 switches
    couplings 3-7 off, and `nlspch` is what tells the two apart.
    """
    return 999 if not field.strip() else read_i(field)


def read_config_card(card, kset, index):
    """Decode a configuration card, `FORMAT` at `FOR/rcg11kd.f:3334`.

    Cowan's card ends at column 80 and RCN2 still writes exactly that,
    so the record is normally short; the A3 cutoff field then reads as
    blanks, which is what `NLTCHK` treats as "use NLT11" (`:3326-3332`).
    """
    c = card.ljust(85)
    slots = []
    for i in range(NOSUBC_MAX):
        at = 5 * i
        symbol = c[at:at + 1]
        if symbol == ' ':
            symbol = LSYMB[0]       # a blank slot is S 0, as RCG reads it
        slots.append((symbol.upper(), read_i(c[at + 1:at + 3])))
    return ConfigCard(kset=kset, index=index, slots=slots,
                      name=c[40:58], eav=read_f(c[58:71], 3),
                      ekin=read_f(c[71:80], 4), nlt=read_i(c[82:85]),
                      raw=card)


def _values(card, start, count):
    """`count` consecutive (F9.4, I1) pairs beginning at column `start`."""
    values, groups = [], []
    for i in range(count):
        at = start + 10 * i
        values.append(read_f(card[at:at + 9], 4))
        groups.append(read_i(card[at + 9:at + 10]))
    return values, groups


def _tail(card):
    """The A2 kind and 4I2 scale factors that close a parameter card."""
    return card[70:72], tuple(read_i(card[72 + 2 * i:74 + 2 * i])
                              for i in range(4))


def _is_first_line(card):
    """True on the first line of a parameter or CI card.

    Columns 71-72 hold the A2 scaling tag ('hr'), which only the first
    line of a card carries; a continuation is seven value/group pairs
    and stops at column 70.  This is the test `s11.pl:270` uses.
    """
    return len(card) > 70 and card[70:72].strip() != ''


def read_ing11(path=files.ING11, folder=None):
    """Read a whole ING11 deck.

    The deck is a stream of records whose meaning depends on position,
    so it is read the way RCG reads it: the prelude cards first, then
    the main card, then NSCONF(2,K) configuration cards per parity set,
    then that set's parameter and CI cards.  Anything after the last
    set -- the radial-integral cards and Cowan's sentinels -- is kept
    verbatim in `trailer`; nothing in this package interprets it.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    # Trailing blanks are kept: `write_ing11` puts the cards it did not
    # rebuild back as the bytes they arrived as, and Cowan's sentinel
    # cards at the end of a deck are padded out to column 80.
    lines = files.read_lines(path, strip=False)
    deck = Ing11(path=str(path))
    i = 0

    # Prelude.  The main card is the one with KCPL = 1 (LS) or 2 (JJ);
    # everything before it is optional and is kept raw.  KCPL = 3 is
    # Cowan's first control card (`:485`), > 10 his second (`:549`), 9
    # the developer card, which consumes a second line (`:509-513`), and
    # 0 a rescale card, which RCG reads and then goes back for the main
    # card (`:726-745`).
    while i < len(lines):
        kcpl = read_i(lines[i][:5])
        if kcpl in (1, 2):
            break
        deck.prelude.append(lines[i])
        i += 1
        if kcpl == 9 and i < len(lines):    # LSJJSF consumes a second card
            deck.prelude.append(lines[i])
            i += 1
    if i >= len(lines):
        raise Ing11Error('%s ends before its main card' % path)

    deck.main = read_main(lines[i])
    i += 1

    # Both lists are indexed by parity set, so they are laid out first
    # and filled in place: appending as the sets are met puts set 2's
    # block at index 0 whenever set 1 is the empty one.
    deck.parameters = [[], []]
    deck.ci = [[], []]

    for k in (1, 2):
        n = deck.main.nconfigurations(k)
        if not n:
            continue
        for j in range(n):
            if i >= len(lines):
                raise Ing11Error(
                    '%s ends after %d of the %d configuration cards of '
                    'parity set %d' % (path, j, n, k))
            deck.configurations.append(read_config_card(lines[i], k, j + 1))
            i += 1

    for k in (1, 2):
        n = deck.main.nconfigurations(k)
        if not n:
            continue
        labels = set(c.label for c in deck.set(k))
        cards, ci, i = _read_parameter_block(lines, i, n, labels, path, k)
        deck.parameters[k - 1] = cards
        deck.ci[k - 1] = ci

    deck.trailer = lines[i:]
    return deck


def _read_parameter_block(lines, i, nconf, labels, path, k):
    """One parity set's parameter cards, then its CI cards.

    RCG reads exactly NSCONF(2,K) single-configuration cards
    (`:10782-10787`) and then one CI card per pair with parameters
    (`:10894-10900`), so the count is what separates the two groups --
    not a guess about the text.  A CI card is confirmed by its name,
    which is two of this set's own configuration labels with a '-'
    between them; anything else ends the block, which is how the
    radial-integral cards and Cowan's sentinels stay in the trailer.
    """
    cards, ci = [], []
    for j in range(nconf):
        if i >= len(lines) or not _is_first_line(lines[i]):
            raise Ing11Error(
                '%s: parity set %d has %d configuration cards but only %d '
                'parameter cards' % (path, k, nconf, j))
        rec, i = _read_card(lines, i, ParameterCard)
        i = _eat_r_integrals(lines, i, rec)
        cards.append(rec)
    while i < len(lines) and _is_first_line(lines[i]):
        if not _ci_labels(lines[i][:18], labels):
            break
        rec, i = _read_card(lines, i, CiCard)
        i = _eat_r_integrals(lines, i, rec)
        ci.append(rec)
    return cards, ci, i


def _eat_r_integrals(lines, i, card):
    """Take the radial-integral cards that follow a parameter card.

    RCN2 does not put every radial-integral card at the end of the deck.
    A same-parity block -- the //R2// cards a quadrupole run needs -- is
    written straight after the parameter card of the parity set it
    belongs to, so on an IQUAD deck one of them sits between parity set
    1's parameters and parity set 2's.  RCG reads them there (`SPECTR`,
    `FOR/rcg11kd.f:12824-12851`), so this reader has to as well; they
    are kept on the card they follow and `write_ing11` puts them back in
    the same place, which keeps the round trip byte for byte.
    """
    while i < len(lines) and is_r_integral(lines[i]):
        card.trailing.append(lines[i])
        i += 1
    return i


def _ci_labels(name, labels):
    """The two labels of a CI card's name, or None if it is not one."""
    left, sep, right = name.partition('-')
    if not sep:
        return None
    left, right = left.strip(), right.strip()
    if left in labels and right in labels:
        return left, right
    return None


def _read_card(lines, i, cls):
    """One parameter or CI card and its continuations.

    Both share `2A6,A7,1X` for the name -- 18 columns, with RCN2's own
    parameter count in 19-20, which RCG skips over as part of the A7 and
    the 1X.  A parameter card then has EAV as F10.5 followed by four
    (F9.4,I1) pairs (`:10790`); a CI card has no EAV and five pairs
    (`:10902`).  Continuations are seven pairs each (`:10834`).
    """
    card = lines[i].ljust(80)
    i += 1
    kind, factors = _tail(card)
    rec = cls(name=card[:18], np=read_i(card[18:20]), kind=kind,
              factors=factors, raw=[lines[i - 1]])
    if cls is ParameterCard:
        rec.values, rec.groups = _values(card, 30, 4)
        rec.values.insert(0, read_f(card[20:30], 5))
        rec.groups.insert(0, 0)     # EAV has no LSF group of its own
        rec.fields = [card[20 + 10 * j:30 + 10 * j] for j in range(5)]
    else:
        rec.values, rec.groups = _values(card, 20, 5)
        rec.fields = [card[20 + 10 * j:30 + 10 * j] for j in range(5)]
    while i < len(lines) and lines[i].strip() and not _is_first_line(lines[i]):
        cont = lines[i].ljust(70)
        values, groups = _values(cont, 0, 7)
        rec.values.extend(values)
        rec.groups.extend(groups)
        rec.fields.extend(cont[10 * j:10 * j + 10] for j in range(7))
        rec.raw.append(lines[i])
        i += 1
    del rec.values[max(rec.np, 1):]
    del rec.groups[max(rec.np, 1):]
    del rec.fields[max(rec.np, 1):]
    return rec, i


# --------------------------------------------------------------------------
# tying a slot to an nl
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class SlotShell:
    """One configuration-card slot, once IN36 has said what `n` it is."""

    slot: int                       # 1-based position on the card
    n: int
    l: int
    occ: int

    @property
    def nl(self):
        return '%d%s' % (self.n, LBCD[self.l])

    @property
    def is_open(self):
        return 0 < self.occ < 2 * (2 * self.l + 1)

    def __str__(self):
        return '%s%s' % (self.nl, self.occ if self.occ != 1 else '')


def _wanted_shells(configuration, core):
    """{(l, occupation): [n, ...]} for what the card should carry.

    The closed core the whole run shares is dropped, because RCN2 drops
    it when it writes the card.
    """
    dropped = {(s.n, s.l) for s in core}
    wanted = {}
    for shell in configuration.shells:
        if shell.occ and (shell.n, shell.l) not in dropped:
            wanted.setdefault((shell.l, int(shell.occ)), []).append(shell.n)
    for key in wanted:
        wanted[key].sort()
    return wanted


def slot_order(cards, by_label, core=()):
    """The order the deck's slots have to be read in to be ascending n,l.

    RCN2 writes a configuration's subshells into the card in ascending
    n,l, so reading the slots left to right normally gives them in that
    order and `match_shells` can hand out the `n` values in the same
    direction.  A **reordered deck permutes the slots on purpose**, and
    then left to right is the wrong direction: U90's reordered `2s2p`
    card is `S 1  P 1  S 1`, and taking the two S slots in card order
    calls the first one 1s when the deck means it to be 2s.  Nothing on
    that one card can tell them apart -- both are (l=0, occupation 1).

    The rest of the deck can.  The same permutation applies to every
    card, so a configuration where the occupations happen to differ
    fixes it for the ones where they do not: the reordered `2s` card is
    `S 1  P 0  S 2`, whose two S slots hold 1 and 2 electrons, so slot 3
    is 1s and slot 1 is 2s, and therefore slot 3 comes before slot 1
    everywhere.  Every such card contributes its pairwise constraints
    and the order is their topological sort, ties broken by slot number
    so that a deck which constrains nothing keeps Cowan's own order.

    Contradictory constraints mean the cards cannot all be one
    permutation of one canonical order, which is not a deck this reader
    can make sense of; it says so rather than picking one.
    """
    after = {}                          # slot -> slots that must follow it
    for card in cards:
        configuration = by_label.get(card.label)
        if configuration is None:
            continue
        wanted = _wanted_shells(configuration, core)
        # Only a slot whose (l, occupation) is unique in its
        # configuration says which subshell it is on its own; a repeated
        # one is exactly the ambiguity being resolved here.
        fixed = []
        for slot, l, occ in card.occupied:
            ns = wanted.get((l, int(occ)), ())
            if len(ns) == 1:
                fixed.append((slot, (ns[0], l)))
        for a, nl_a in fixed:
            for b, nl_b in fixed:
                if nl_a < nl_b:
                    after.setdefault(a, set()).add(b)

    slots = list(range(1, NOSUBC_MAX + 1))
    before = dict((s, set()) for s in slots)
    for a, bs in after.items():
        for b in bs:
            before[b].add(a)
    out = []
    remaining = set(slots)
    while remaining:
        ready = sorted(s for s in remaining if not before[s] & remaining)
        if not ready:
            raise Ing11Error(
                'the configuration cards do not agree on one subshell order: '
                'slots %s constrain each other in a circle'
                % ', '.join(str(s) for s in sorted(remaining)))
        out.append(ready[0])
        remaining.discard(ready[0])
    return tuple(out)


def match_shells(card, configuration, core=(), order=None):
    """Give every occupied slot of `card` the `n` IN36 gives it.

    RCN2 writes the ING11 card from the IN36 configuration, dropping the
    closed core the whole run shares -- `core`, which is what
    `in36.common_closed_shells()` returns -- so the card's occupied
    slots are a permutation of what is left, with the same
    (l, occupation) multiset.  Matching on that multiset and then
    assigning `n` in ascending order within each l is a permutation of
    the slots and never a positional assumption, which is the whole
    point: a reordered deck permutes them deliberately.

    `order` is the sequence the slots are handed their `n` values in --
    `slot_order()` above, which is the deck's own order and is what a
    reordered deck needs.  Left to right is the default and is right for
    any deck RCN2 wrote itself.  The list comes back in slot order
    either way, because that is the order genealogy.txt names subshells
    in.

    A mismatch prints both lists and raises; nothing is guessed.
    """
    wanted = _wanted_shells(configuration, core)

    occupied = card.occupied
    if order is not None:
        rank = {slot: i for i, slot in enumerate(order)}
        occupied = sorted(occupied, key=lambda s: rank.get(s[0], s[0]))

    used = dict((k, 0) for k in wanted)
    out = []
    for slot, l, occ in occupied:
        key = (l, int(occ))
        available = wanted.get(key, [])
        taken = used.get(key, 0)
        if taken >= len(available):
            raise Ing11Error(
                'configuration %s: the card has a %s%d subshell that %s does '
                'not.\n  card:  %s\n  IN36:  %s'
                % (card.name.strip(), LBCD[l], occ, configuration.label,
                   _card_text(card), _conf_text(configuration)))
        used[key] = taken + 1
        out.append(SlotShell(slot=slot, n=available[taken], l=l, occ=int(occ)))

    for key, ns in wanted.items():
        if used.get(key, 0) != len(ns):
            l, occ = key
            raise Ing11Error(
                'configuration %s: %s has a %s%d subshell the card does '
                'not.\n  card:  %s\n  IN36:  %s'
                % (card.name.strip(), configuration.label, LBCD[l], occ,
                   _card_text(card), _conf_text(configuration)))
    out.sort(key=lambda s: s.slot)
    return out


# --------------------------------------------------------------------------
# writing a parameter card back
# --------------------------------------------------------------------------

# A parameter or CI card always carries at least five value fields, even
# where RCN2's count in columns 19-20 is smaller: the first line's
# FORMAT has five of them and the A2 tag behind them is what marks a
# first line, so the field has to be filled.  RCN2 fills a spare field
# with a plain zero on a parameter card and with 0.0000 on a CI card,
# and gives the spares of a CI card group code 5 -- see `WORK/U90/ing11`
# lines 7 and 11, and `s11.pl:243-251`, which is where these came from.
FIELDS_PER_LINE = 5
FIELDS_PER_CONTINUATION = 7
_SPARE_PARAMETER = ('        0', 0)
_SPARE_CI = ('   0.0000', 5)


def format_value(value, index, ci):
    """One nine-column value field of a parameter or CI card.

    A CI card's fields are plain `F9.4` (`:10902`), a parameter card's
    are RCN2's scaled-integer encoding (`fortfmt.write_param`, the rule
    quoted in `s11.pl:202-221`), and the first field of a parameter card
    is EAV, which is never written in the exponential form because it is
    never small.

    `F9.4` is written through `write_f_wide`, which drops decimals
    rather than the number.  `s11.pl:228` writes a bare `%9.4f` there
    and so puts ten characters into a nine-column field for any CI
    parameter at or past 1000 in magnitude -- `WORK/U90`'s
    `131D1233 = -1369.505304` is one -- which shifts every field behind
    it and corrupts the card.  That is the defect the author reported.
    """
    if ci:
        return write_f_wide(value, 9, 4)
    return write_param(value, first=(index == 0))


def _unchanged(card, i, ci):
    """The ten columns field `i` arrived in, if it still says the same.

    A field whose value and group have not moved keeps the characters it
    was punched in.  That is not only tidiness: RCN2 and the Perl
    utilities have each written these fields in forms the other does not
    produce -- `  0.00e+0` against ` 0.00E+00` for a zero, both of which
    RCG reads as zero -- and re-encoding every field of every card it
    touches would make each of these programs rewrite parts of a deck
    that it has no business changing.  What is encoded here is what has
    actually been given a new number.
    """
    if i >= len(card.fields):
        return None
    text = card.fields[i]
    if i == 0 and not ci:
        # EAV is F10.5 across the whole field: RCG reads the tenth
        # column as part of the number (`:10790`), so it has no group.
        return text if read_f(text, 5) == card.values[i] else None
    if read_f(text[:9], 4) == card.values[i] \
            and read_i(text[9:10]) == card.groups[i] % 10:
        return text
    return None


def format_card(card, ci=False):
    """A parameter or CI card as lines, from its values and groups.

    Columns 1-20 (the 18-character name and RCN2's parameter count) and
    the tail past column 70 (the A2 scaling tag and its 4I2 factors) are
    carried across as the characters they already are: nothing in this
    package has a reason to re-encode them, and a card that comes back
    from `read_ing11` unchanged has to come out byte for byte the same.
    """
    header = card.raw[0][:20].ljust(20)
    tail = card.raw[0][70:]
    spare = _SPARE_CI if ci else _SPARE_PARAMETER
    count = max(len(card.values), FIELDS_PER_LINE)

    fields = []
    for i in range(count):
        if i < len(card.values):
            kept = _unchanged(card, i, ci)
            fields.append(kept if kept is not None else
                          '%s%d' % (format_value(card.values[i], i, ci),
                                    card.groups[i] % 10))
        else:
            fields.append('%s%d' % spare)

    lines = [header + ''.join(fields[:FIELDS_PER_LINE]) + tail]
    rest = fields[FIELDS_PER_LINE:]
    for at in range(0, len(rest), FIELDS_PER_CONTINUATION):
        lines.append(''.join(rest[at:at + FIELDS_PER_CONTINUATION]))
    return lines


def set_values(card, values, groups=None):
    """Give a card new values, and mark it for re-encoding.

    A card nobody has edited is written back out as the characters it
    arrived as, so a tool that touches three cards of a deck changes
    three cards of the file and nothing else.  Going through here is
    what says a card is no longer what its text says.
    """
    card.values[:] = list(values)
    if groups is not None:
        card.groups[:] = list(groups)
    card.edited = True
    return card


def take_card(target, source):
    """Replace one card with another deck's, text and all.

    Used where a whole card is being carried across between decks
    rather than recomputed -- `subst_params_ing11` -- and the point is
    that the values are moved as the characters RCN2 punched, never
    re-encoded.  That is the rule `toolbox/reorder.py` already works
    to, and the reason a value cannot lose a digit in transit.
    """
    target.raw = list(source.raw)
    target.name = source.name
    target.np = source.np
    target.values = list(source.values)
    target.groups = list(source.groups)
    target.kind, target.factors = source.kind, source.factors
    target.edited = False
    return target


def write_ing11(deck):
    """The whole deck as lines, ready to be written out.

    Everything the reader did not have to understand -- the prelude, the
    main card, the configuration cards and the trailer -- goes back as
    the text it came in as, and so does every parameter card that was
    not edited.  A deck nothing has touched therefore round-trips byte
    for byte; `test_ing11.py` holds that to all five reference decks.
    """
    lines = list(deck.prelude)
    lines.append(deck.main.raw)
    lines.extend(card.raw for card in deck.configurations)
    for k in (1, 2):
        for card in deck.parameters[k - 1]:
            lines.extend(format_card(card) if card.edited else card.raw)
            lines.extend(card.trailing)
        for card in deck.ci[k - 1]:
            lines.extend(format_card(card, ci=True) if card.edited
                         else card.raw)
            lines.extend(card.trailing)
    lines.extend(deck.trailer)
    return lines


def _card_text(card):
    return ' '.join('%s%d' % (sym, occ) for sym, occ in card.slots if occ)


def _conf_text(configuration):
    return ' '.join(str(s) for s in configuration.shells if s.occ)


# --------------------------------------------------------------------------
# the radial transition integrals of the trailer
# --------------------------------------------------------------------------

@dataclass
class RCard:
    """One reduced multipole matrix element of the deck's trailer.

    RCN2 punches these (`FOR/RCN2K.F:668`, FORMAT 43/44) and RCG reads
    them back at `FOR/rcg11kd.f:12851`:

        FORMAT (3A6,2X,3A6,F12.5, A6,A8, F6.4, A2,2E4.1)

    so the value is columns 39-50 and everything behind it -- the
    `(  //Rk//  )` tag, FRAC, the `hr` mark and the two overlap
    percentages -- is columns 51-80.  The multipole order k is a single
    digit in column 58; RCG reads it that way itself, with
    `FORMAT (50X,A1,6X,I1)` at `:12824`, and refuses the card if column
    51 is not a bracket.

    `where` and `at` are the list the card's line sits in and its index
    in it, which is where it is written back to.  There are two such
    lists: `Ing11.trailer`, and the `trailing` of a parameter card, for
    the same-parity blocks RCN2 writes in the middle of the deck.
    """

    where: list
    at: int
    order: int                      # k of //Rk//: 1 = dipole, 2 = quadrupole
    value: float
    raw: str

    @property
    def names(self):
        """The two configuration labels the integral connects."""
        return self.raw[0:18].strip(), self.raw[20:38].strip()


#: Columns 39-50, `F12.5`, and the two that identify the card.
_R_VALUE = slice(38, 50)
_R_BRACKET = 50
_R_ORDER = 57


def is_r_integral(line):
    """Whether a trailer line is a radial transition-integral card.

    The test is RCG's own: a bracket in column 51 and a digit in column
    58, inside the `(  //Rk//  )` tag.  Cowan's sentinel cards
    (`-55555555.`, `-99999999.`, `-1`) fail it, as does anything else a
    deck may carry after the parameters.
    """
    return (len(line) >= 64 and line[_R_BRACKET] == '('
            and line[54:57] == '//R' and line[_R_ORDER].isdigit()
            and line[58:60] == '//' and line[63] == ')')


def read_r_integrals(deck, order=None):
    """The trailer's radial-integral cards, or those of one multipole.

    `order=2` gives the E2 integrals, which is what a forbidden-line
    calculation varies.  Every card of that order is returned wherever
    it sits in the trailer; `calc_TP_M1E2.pl:203-216` instead took the
    first one and then consecutive lines only, so a deck whose R2 block
    is split by another multipole's would have lost the rest.
    """
    out = []
    for where in _r_lists(deck):
        for at, line in enumerate(where):
            if not is_r_integral(line):
                continue
            k = int(line[_R_ORDER])
            if order is not None and k != order:
                continue
            out.append(RCard(where=where, at=at, order=k,
                             value=read_f(line[_R_VALUE], 5), raw=line))
    return out


def _r_lists(deck):
    """Every list of raw lines a radial-integral card may sit in, in
    the order the deck writes them out."""
    for k in (1, 2):
        for card in deck.parameters[k - 1]:
            yield card.trailing
        for card in deck.ci[k - 1]:
            yield card.trailing
    yield deck.trailer


def set_r_value(card, value):
    """Put a new value in one radial-integral card, changing nothing else.

    Columns 1-38 and 51-80 come back as the bytes they arrived in, so a
    card whose value has not moved round-trips exactly -- the same rule
    the parameter cards are written under (`format_card`).  The value
    goes through `write_f_wide`, which gives up a decimal rather than
    the number when F12.5 will not hold it.
    """
    text = write_f_wide(value, 12, 5)
    card.where[card.at] = card.raw[:38] + text + card.raw[50:]
    card.value = value
    return card.where[card.at]
# --------------------------------------------------------------------------
# the cuts that censor the printed line list
# --------------------------------------------------------------------------

#: The range of DMIN in which the cut is on gA rather than on the
#: squared matrix element, `FOR/rcg11kd.f:13296`.
DMIN_GA_FROM, DMIN_GA_TO = 1.0, 50.0


def dmin_cut(dmin):
    """What a DMIN value means, as `(kind, threshold)`.

    RCG turns the field into a cut at `FOR/rcg11kd.f:13294-13302`.  A
    magnitude above 1 and not above 50 is an **exponent**: the cut is
    `gA < 10**int(|DMIN|)`, applied to the printed GA(SEC-1) column at
    `:13764`, `:13931`, `:14118` and `:14473`.  Anything else is a cut
    on `D0`, the squared reduced matrix element, at `:13656`.  Only the
    integer part matters in the first case, which is why an offset in
    decades is a whole number there.

    A negative DMIN additionally deletes upward transitions from the
    level sorts (NEGTP); the magnitude is what sets the cut, and that
    is what is returned.
    """
    size = abs(dmin)
    if DMIN_GA_FROM < size <= DMIN_GA_TO:
        return 'gA', 10.0 ** int(size)
    return 'D0', size


def set_dmin(deck, value):
    """Put a new DMIN in columns 61-65 of the main card.

    Everything else on the card comes back as the bytes it arrived in,
    the rule `set_r_value` and `format_card` are written under.  The
    field is `E5.1`, so a value is written with its decimal point in it
    and reads back as punched.
    """
    raw = deck.main.raw.ljust(65)
    deck.main.raw = (raw[:_DMIN.start] + write_f_wide(value, 5, 1)
                     + raw[_DMIN.stop:]).rstrip()
    deck.main.dmin = value
    return deck.main.raw


def set_window(deck, flbdn, flbdx):
    """Put a new FLBDN/FLBDX in the KCPL=3 card, columns 11-20.

    The fields are `E5.0` and hold angstrom, so a widened window has to
    fit five columns; `write_f_wide` gives up the decimal rather than
    the number, which costs at most half an angstrom in the last digit
    -- nothing beside the fractional margin a caller widens by.  A deck
    with no KCPL=3 card has no window to widen and is an error here:
    RCG's own defaults, 0.001 to 500000 A, print everything already.
    """
    for at, line in enumerate(deck.prelude):
        if line[:5].strip() and read_i(line[:5]) == 3:
            padded = line.ljust(20)
            deck.prelude[at] = (padded[:_FLBDN.start]
                                + write_f_wide(flbdn, 5, 1)
                                + write_f_wide(flbdx, 5, 1)
                                + padded[_FLBDX.stop:]).rstrip()
            return deck.prelude[at]
    raise Ing11Error('%s has no KCPL=3 card, so it has no wavelength window'
                     % (deck.path or files.ING11))
