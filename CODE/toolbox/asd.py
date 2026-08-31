"""NIST ASD notation, written inline.

This module is only about *how a designation is written*.  It knows
nothing about Cowan's files; the callers in phase 6 read `genealogy.txt`
and hand the pieces over already decided.  Keeping it separate is what
lets the notation be pinned to `ASD_syntax_examples.txt` and to the
Handbook (`Handbook/010.tex` sections 010-9-1 and 010-8-2 .. 010-8-5)
by tests that mention no Cowan output at all.

Everything is inline: no superscripts, no subscripts, no markup.  The
rules, with the source of each:

* **Tokens are joined by a period** -- every subshell and every
  parenthesised intermediate, without exception
  (`ASD_syntax_examples.txt`).
* **Odd parity is a `*` on the term**, never on the configuration:
  `2P*`, `2D2*`, `(1,1/2)*`, `2[3]*`.  A term inside the configuration
  carries the parity of the electrons summed into it so far, so
  `(2P*<1/2>)` is right and the configuration itself stays unmarked.
* **An intermediate LS term gives its J in angle brackets**:
  `4f2.(3H<4>).6s`, `2s.2p.(2P*<1/2>).(0).3s`.  Literature subscripts
  that value; angle brackets are the ASD inline stand-in.
* **An intermediate that *is* a J stays in parentheses** with nothing
  added: `1s.2s.(1).2p-`, `1s.2p-.(0).2p+`.  It is written only where
  it says something the final term does not: the last group's J is
  already the first half of `(J1,J2)`, so `1s.2p+2` carries no `(2)`.
* **A jj subshell** is n, l, sign, occupation: `2p-2`, `2p+2`, where
  `-` is j = l - 1/2 and `+` is j = l + 1/2.  An s subshell never takes
  a sign -- j = 1/2 is its only value (`010.tex:511-513`) -- so one
  given here is dropped rather than printed.
* **LK's bare total-L letter belongs to the term**, not to the
  configuration: `2s2.2p2.(3P).4f | D 2[3]* | 5/2`.  The Handbook
  writes it at the end of the configuration (`010.tex:550-558`), but
  it is a quantum number of the term, and leaving it in the other
  column collapses `D 2[3]*` and `F 2[3]*` into one entry in any list
  of distinct terms -- the same defect the author fixed by moving the
  parity mark off J.  Author's decision, 2026-08-14.
* Final terms by scheme (`010.tex:585-589`): LS `2S+1 L`, J1J2
  `(J1,J2)`, J1L2 and LS1 `2S2+1[K]`.
* A Nielson-Koster index follows the L letter directly, before the
  parity mark: `2D2*`.
* **Seniority rides on the subshell that owns it**, as a right
  subscript: `4f+4<2>`.  Angle brackets are the inline stand-in for a
  literature subscript here exactly as they are for a parent term's J;
  which side of the quantity they sit on says which subscript is meant,
  and a parent J is always on the right of its term, so the two never
  collide.  The alternative -- a left subscript on the group's J, as
  GRASP output is printed -- was rejected because every RCG coupling is
  a running sum: only the *last* subshell's own J is ever written, so a
  subshell in the middle of the configuration would have no anchor.
  Author's decision, 2026-08-14.

Every coupling RCG prints is a *running* sum -- each intermediate
couples everything to its left -- so a designation is a flat list of
tokens and the Handbook's explicitly bracketed form (`010.tex:492`,
`[(2s.2p-)0.2p+]3/2`) is not needed and is not implemented.  If a
nested tree ever has to be written, that is a question for the author,
not a rule to invent here.
"""

from dataclasses import dataclass
from fractions import Fraction

# Term letters, RCG's own table (`FOR/rcg11kd.f:1081-1083`); the index of
# a letter is L.  Orbital letters are lower case and stop at l = 20
# (`FOR/RCN36K.F:1142-1143`, mirrored in `in36.LBCD`).
L_LETTERS = 'SPDFGHIKLMNOQRTUVWXYZABC'
ORBITALS = 'spdfghiklmnoqrtuvwxyz'

DOT = '.'
SPACE = ' '


class AsdError(Exception):
    """A designation that cannot be written as ASD writes them."""


# --------------------------------------------------------------------------
# quantum numbers
# --------------------------------------------------------------------------

def format_j(j):
    """'0', '3', '1/2', '15/2' -- whole or half, never a decimal point."""
    try:
        f = Fraction(j)             # exact: 0.3 must not be rounded to 1/2
    except (TypeError, ValueError):
        raise AsdError('%r is not a J value' % (j,))
    if f.denominator not in (1, 2):
        raise AsdError('J = %r is neither whole nor half-integral' % (j,))
    if f < 0:
        raise AsdError('J = %r is negative' % (j,))
    if f.denominator == 1:
        return str(f.numerator)
    return '%d/%d' % (f.numerator, f.denominator)


def format_j2(j2):
    """The same, from a doubled J -- what the true jj table prints."""
    return format_j(Fraction(int(j2), 2))


def l_letter(l):
    """L as its capital letter; an int or a letter is accepted."""
    if isinstance(l, str):
        letter = l.strip().upper()
        if len(letter) != 1 or letter not in L_LETTERS:
            raise AsdError('%r is not an L symbol' % (l,))
        return letter
    if not 0 <= l < len(L_LETTERS):
        raise AsdError('L = %r has no letter in RCG\'s table' % (l,))
    return L_LETTERS[l]


def _orbital(l):
    """l as its lower-case orbital letter; an int or a letter is accepted."""
    if isinstance(l, str):
        letter = l.strip().lower()
        if len(letter) != 1 or letter not in ORBITALS:
            raise AsdError('%r is not an orbital letter' % (l,))
        return letter
    if not 0 <= l < len(ORBITALS):
        raise AsdError('l = %r has no orbital letter' % (l,))
    return ORBITALS[l]


def _l_value(l):
    """The integer l behind either spelling."""
    return ORBITALS.index(_orbital(l))


def _seniority(v):
    """A seniority as it is written: `<2>`, or nothing where there is none.

    Angle brackets are the inline stand-in for a literature subscript,
    exactly as they are for an intermediate term's J.  Which side of the
    quantity they sit on says which subscript is meant.
    """
    if v is None or v == '':
        return ''
    v = int(v)
    if v < 0:
        raise AsdError('seniority %r is negative' % (v,))
    return '<%d>' % v


def _occupation(occ):
    """Occupation as it is written after the subshell; 1 is not written."""
    if occ == 1:
        return ''
    if occ < 0:
        raise AsdError('occupation %r is negative' % (occ,))
    if float(occ).is_integer():
        return str(int(occ))
    return '%g' % occ


# --------------------------------------------------------------------------
# the pieces of a designation
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Subshell:
    """One subshell as it appears in a configuration: `4f3`, `2p-2`."""

    n: int
    l: object                       # int or letter, either spelling
    occ: float = 1
    sign: str = ''                  # '-' or '+' for a jj subshell
    v: object = None                # seniority, where J alone is not enough

    sep = DOT

    def __post_init__(self):
        if self.sign not in ('', '-', '+'):
            raise AsdError('%r is not a j sign' % (self.sign,))
        if self.n < 1:
            raise AsdError('n = %r' % (self.n,))

    @property
    def letter(self):
        return _orbital(self.l)

    @property
    def lval(self):
        return _l_value(self.l)

    @property
    def written_sign(self):
        """The sign as written -- an s subshell has only one j, so none."""
        return '' if self.lval == 0 else self.sign

    @property
    def text(self):
        return '%d%s%s%s%s' % (self.n, self.letter, self.written_sign,
                               _occupation(self.occ), _seniority(self.v))

    def __str__(self):
        return self.text


@dataclass(frozen=True)
class LSTerm:
    """`2S`, `2P*`, `2D2*`, and as an intermediate `3H<4>`, `2P*<1/2>`."""

    mult: int                       # 2S+1
    l: object
    index: str = ''                 # Nielson-Koster index, after the letter
    odd: bool = False
    j: object = None                # given only where the term is a parent

    sep = DOT

    @property
    def text(self):
        out = '%d%s%s' % (self.mult, l_letter(self.l), self.index)
        if self.odd:
            out += '*'
        if self.j is not None:
            out += '<%s>' % format_j(self.j)
        return out

    def __str__(self):
        return self.text


@dataclass(frozen=True)
class JTerm:
    """An intermediate that is itself a J: the `(1)` of `1s.2s.(1).2p-`.

    It carries no seniority: it is a running J, and a seniority belongs
    to one subshell's own state.  It is written on the subshell instead
    -- see `Subshell.v`.
    """

    j: object

    sep = DOT

    @property
    def text(self):
        return format_j(self.j)

    def __str__(self):
        return self.text


@dataclass(frozen=True)
class PairTerm:
    """The final term of J1J2 coupling: `(1,3/2)`, `(0,1/2)*`.

    J1 is the running J of everything but the last subshell and J2 is
    that subshell's own; neither takes a seniority, which is written on
    the subshell itself.
    """

    j1: object
    j2: object
    odd: bool = False

    sep = DOT

    @property
    def text(self):
        return '(%s,%s)%s' % (format_j(self.j1), format_j(self.j2),
                              '*' if self.odd else '')

    def __str__(self):
        return self.text


@dataclass(frozen=True)
class KTerm:
    """The final term of J1L2 (J1K) and LS1 (LK) coupling: `2[3]*`.

    LK's total-L letter, where there is one, stands in front of it and
    is part of the term: `D 2[3]*`.  See the module docstring for why it
    is not left at the end of the configuration.
    """

    mult: int                       # 2S2+1
    k: object
    odd: bool = False
    l: object = None                # LK's total L, written before the term

    sep = DOT

    @property
    def text(self):
        out = '%d[%s]%s' % (self.mult, format_j(self.k),
                            '*' if self.odd else '')
        if self.l is None:
            return out
        return l_letter(self.l) + SPACE + out

    def __str__(self):
        return self.text


@dataclass(frozen=True)
class Intermediate:
    """Any term standing inside a configuration, hence in parentheses.

    `space=True` sets it off with a space instead of a period, which is
    how the Handbook marks a term formed only after the ones to its left
    (`010.tex:393-395`).  Cowan's schemes never need it; it is here so
    that a caller which does can say so rather than fake it.
    """

    term: object
    space: bool = False

    @property
    def sep(self):
        return SPACE if self.space else DOT

    @property
    def text(self):
        return '(%s)' % self.term.text

    def __str__(self):
        return self.text


# --------------------------------------------------------------------------
# putting them together
# --------------------------------------------------------------------------

def parity(subshells):
    """0 even, 1 odd, from sum(l * occupation) over the subshells given."""
    return int(sum(s.lval * s.occ for s in subshells)) % 2


def configuration(tokens):
    """The configuration column: tokens joined by their own separators."""
    out = []
    for i, tok in enumerate(tokens):
        try:
            text = tok.text
        except AttributeError:
            raise AsdError('%r is not a designation token' % (tok,))
        out.append(text if i == 0 else tok.sep + text)
    return ''.join(out)


@dataclass
class Level:
    """One row of an ASD table: configuration, term, J."""

    tokens: list                    # Subshell / Intermediate / LSymbol
    term: object                    # LSTerm, PairTerm or KTerm
    j: object

    @property
    def configuration(self):
        return configuration(self.tokens)

    @property
    def term_text(self):
        return self.term.text

    @property
    def j_text(self):
        return format_j(self.j)

    def columns(self):
        return (self.configuration, self.term_text, self.j_text)

    def __str__(self):
        return ' | '.join(self.columns())
