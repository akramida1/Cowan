"""Turn one basis state into the three columns of an ASD table.

This is the module the Perl calls the designation builder
(`conv_cowan.pl:2670-3123`).  It is rebuilt here on `genealogy.txt`
rather than on the eight-character label, which is the whole reason the
Fortran writes that file: a label has room for the last subshell's term
and the running total, and for nothing in between, so the Perl had to
guess the middle of a chain from the configuration card.  GENPRT prints
the chain, one field per open subshell, so nothing has to be guessed.

WHAT GOES IN PARENTHESES
------------------------

Cowan couples subshell by subshell, left to right in the order of the
ING11 card.  Written out in full that chain is

    s1 (t1) s2 (t2) (r12) s3 (t3) (r123=T)

-- each subshell's own term, and the running term after each coupling.
Spectroscopic notation drops from that everything a reader can put back:

* a subshell whose occupation allows **one term only** (empty, closed,
  one electron or one hole) carries no term of its own;
* the **first** running term is the first subshell's own term, so it is
  never written twice;
* the **last** running term is the final term, which stands in the term
  column, so it is not repeated inside the configuration.

What survives is exactly the Handbook's notation.  Its example 6
(`010.tex:370-374`, `4f7(8S*) 5d(7D*) 6p  8F`) is the general case: 5d
is one electron and drops its own `2D`, and the `(7D*)` that remains is
the running term of `4f7.5d`.  The two-subshell cases the reference runs
are full of come out the same way -- `1s.2p2.(3P) | 4P` writes 2p2's own
term and drops both running terms, one being 1s's `2S` and the other the
final `4P`.

The same three rules serve every coupling; only what a term *is*
changes.  In J1J2 a subshell's own term is its LS term with its J
(`4f2.(3H<4>).6s | (4,1/2)`) and a running term is a bare J
(`1s.2s.(1).2p-`).  In the JK schemes the chain is LS or J1J2 up to the
second-to-last subshell and the last subshell joins it through K, which
is why those sections print only two numbers per row: the rest of the
chain is the LS or J1J2 chain that `genealogy.txt` already carries.

DISPLAY ORDER
-------------

Tokens come out in the order the deck sums them, because that is the
order the intermediates belong to.  Where no intermediate is written the
order is unobservable, and the subshells are then sorted by (n, l) so
that a reordered deck prints the same label as its parent -- which is
the `1s.(2S).2s2` defect in `WORK/U90/reordered`.  A deck that both
reorders its subshells and needs an intermediate would print in
summation order; `noncanonical()` reports those so they are visible.
"""

from dataclasses import dataclass

from . import asd
from .asd import (Intermediate, JTerm, KTerm, LSTerm, Level, PairTerm,
                  Subshell)
from . import genealogy as gen

# RCG's couplings, by the name the sections carry.
LS = 'LS'
JJ = 'JJ'
TRUEJJ = 'TRUEjj'
JJJK = 'JJJK'
LSLK = 'LSLK'
LSJK = 'LSJK'
LSJLKS = 'LSJLKS'
LSJLSJ = 'LSJLSJ'

#: Which `genealogy.txt` table carries the chain of each coupling.  The
#: JK schemes have no table of their own: their chain up to the
#: second-to-last subshell is the chain of the coupling they are built
#: from, and only the last step is new.
TABLE_OF = {LS: gen.LS, JJ: gen.J1J2, TRUEJJ: gen.TRUEJJ,
            JJJK: gen.J1J2, LSLK: gen.LS, LSJK: gen.LS,
            LSJLKS: gen.LS, LSJLSJ: gen.LS}

#: The couplings whose last subshell joins through K, and whether the
#: printed label carries J1 (as opposed to a bare L or nothing).
K_SCHEMES = (JJJK, LSLK, LSJK, LSJLKS, LSJLSJ)


class DesignateError(Exception):
    """A basis state will not read against the deck it came from."""


def one_term(l, occ):
    """True where l**occ has exactly one LS term, so none is written.

    Empty, closed, one electron or one hole.  It is the same rule that
    decides which subshells true jj leaves unrotated (`true-jj_keys.txt`
    states it there), for the same reason: with one state there is
    nothing to tell apart.
    """
    l = asd._l_value(l)
    return occ in (0, 1, 4 * l + 1, 4 * l + 2)


def one_state(l, occ):
    """True where l**occ has exactly one (LS term, J), so none is written.

    Stricter than `one_term`, and it is what J1J2 needs for a subshell
    that is not the last: one electron in a p, d or f subshell has one
    term but two j values, and in the middle of a chain neither the
    final term nor the running J need pin it down -- which is why
    `2s.2p.(2P*<1/2>).(0).3s` writes it and `1s2.2s` does not.  For the
    *last* subshell J comes back from the final term, so `one_term` is
    the test there, and `5d9.6d | (3/2,5/2)` stays as the author wants
    it.
    """
    l = asd._l_value(l)
    return l == 0 or occ in (0, 4 * l + 2)


@dataclass
class Designation:
    """One basis state, as columns and as the objects behind them."""

    level: object                   # asd.Level
    canonical: bool = True          # False where summation order was kept
                                    # and differs from (n, l) order

    @property
    def configuration(self):
        return self.level.configuration

    @property
    def term(self):
        return self.level.term_text

    def columns(self):
        return self.level.columns()

    def __str__(self):
        return str(self.level)


# --------------------------------------------------------------------------
# the pieces of a chain
# --------------------------------------------------------------------------

def _ls_term(field, odd=False, j=None):
    """A subshell's own LS term, or a running one, from an LS field."""
    return LSTerm(mult=field.mult, l=field.l, index=field.alpha.strip(),
                  odd=odd, j=j)


def _running_ls(field, odd=False, j=None):
    return LSTerm(mult=field.run_mult, l=field.run_l, odd=odd, j=j)


def _j1j2_own(field, odd=False):
    """A subshell's own term in J1J2: its LS term carrying its own J."""
    return LSTerm(mult=field.mult, l=field.l, odd=odd, j=field.j)


def _absorbed(field, ls):
    """True where a subshell's own term is already the running term.

    That happens where everything coupled in before it came to 1S (or to
    J=0), which is to say the running term after this subshell *is* its
    own term under another name.  Writing both would say the same thing
    twice, so the own term is dropped -- `4f3 | 4D*`, not
    `4f3.(4D) | 4D*`.  Cowan's index of a repeated term then rides on
    the term column, which is where the reference runs put it (`2D2*`).
    """
    if not ls:
        # Never in J1J2.  There the term column is a pair of J values and
        # carries no LS term at all, so a subshell's own term is the only
        # thing that says which parent a level came from -- `4f3` alone
        # has fifteen of them at J=3/2.
        return False
    return (field.mult, field.l) == (field.run_mult, field.run_l)


def _parities(open_shells):
    """Per open subshell: its own parity, and the running one after it.

    An intermediate carries the parity of what it stands for, and it is
    not a duplicate of the final term's -- `2s.2p.(3P*).3s | 4P*` has an
    odd parent and an odd total, and `4f2.(3F).5p6.6p | 4D*` an even
    parent and an odd total.  Only closed and empty subshells can be
    ignored here: their l*N is always even.
    """
    own, running, total = [], [], 0
    for s in open_shells:
        odd = int(asd._l_value(s.l) * s.occ) % 2
        total = (total + odd) % 2
        own.append(bool(odd))
        running.append(bool(total))
    return own, running


def _parity_of(shells):
    return int(sum(asd._l_value(s.l) * s.occ for s in shells)) % 2


# --------------------------------------------------------------------------
# the builder
# --------------------------------------------------------------------------

def designate(coupling, state, shells, outer=None, drop=()):
    """The ASD columns of one basis state.

    `state`   a `genealogy.State` from the table `TABLE_OF[coupling]`.
    `shells`  every slot of the configuration card, in card order, as
              `ing11.SlotShell` -- which is where n comes from.
    `outer`   the tail of the OUTG11 row, for the couplings whose last
              step is not in any genealogy table; see `read_outer()`.
    `drop`    (n, l) pairs to leave out of the label where they are
              closed.  The caller decides that; see `printout.py`.
    """
    open_shells = [s for s in shells if s.is_open]
    if len(open_shells) != len(state.fields):
        raise DesignateError(
            'basis state %d has %d fields but its configuration card has %d '
            'open subshells (%s)'
            % (state.number, len(state.fields), len(open_shells),
               ' '.join(str(s) for s in open_shells)))

    odd = bool(_parity_of(shells))
    if coupling == TRUEJJ:
        tokens, term = _true_jj(state, shells, open_shells, odd, drop)
    elif coupling in K_SCHEMES:
        tokens, term = _k_scheme(coupling, state, shells, open_shells,
                                 odd, drop, outer)
    else:
        tokens, term = _chain(coupling, state, shells, open_shells, odd, drop)

    tokens, canonical = _ordered(tokens)
    return Designation(Level(tokens=tokens, term=term, j=None),
                       canonical=canonical)


def _kept(shells, drop):
    """The slots that appear in the label, in card order."""
    out = []
    for s in shells:
        if not s.occ:
            continue                            # an empty slot is not a shell
        if (s.n, asd._l_value(s.l)) in drop and not s.is_open:
            continue
        out.append(s)
    return out


def _chain(coupling, state, shells, open_shells, odd, drop, upto=None,
           j1=None):
    """Tokens and final term for LS or J1J2, by the three rules above.

    `upto` stops the chain short, which is what the JK schemes need: they
    run the chain to the second-to-last subshell and take the last step
    from the printed label instead.  `j1` is that label's J1, written
    into the last intermediate where the scheme carries one.
    """
    ls = coupling in (LS, LSLK, LSJK, LSJLKS, LSJLSJ)
    fields = list(state.fields)
    last = len(fields) if upto is None else upto
    if last < 1:
        raise DesignateError('a basis state with no open subshell')

    own_odd, run_odd = _parities(open_shells)
    by_slot = {}
    for k, s in enumerate(open_shells[:last], start=1):
        field = fields[k - 1]
        parts = []
        unique = (one_term(s.l, s.occ) if ls or k == last
                  else one_state(s.l, s.occ))
        if not unique and not (k == last and _absorbed(field, ls)):
            parts.append(Intermediate(
                _j1j2_own(field, odd=own_odd[k - 1]) if not ls
                else _ls_term(field, odd=own_odd[k - 1])))
        if 2 <= k < last:
            # A running term of two or more subshells: the reader cannot
            # put it back, so it is written.  The first is the first
            # subshell's own term and the last is the final term.
            parts.append(Intermediate(
                JTerm(field.run_j) if not ls
                else _running_ls(field, odd=run_odd[k - 1])))
        by_slot[id(s)] = parts

    if upto is not None:
        # The chain stops short of the last subshell, so its running term
        # is no longer the final term and has to be written -- carrying
        # J1 where the scheme prints one (`2s2.2p2.(3P<0>).5f`).
        tail = fields[last - 1]
        parts = by_slot.setdefault(id(open_shells[last - 1]), [])
        if ls:
            running = _running_ls(tail, odd=run_odd[last - 1], j=j1)
        else:
            running = JTerm(tail.run_j if j1 is None else j1)
        if _absorbed(tail, ls) or one_term(open_shells[last - 1].l,
                                           open_shells[last - 1].occ):
            # The running term is this subshell's own term under another
            # name, so it replaces it rather than following it.
            parts[:] = [Intermediate(running)]
        else:
            parts.append(Intermediate(running))

    tokens = []
    for s in _kept(shells, drop):
        tokens.append(Subshell(n=s.n, l=s.l, occ=s.occ))
        tokens.extend(by_slot.get(id(s), ()))

    if upto is not None:
        return tokens, None
    tail = fields[-1]
    if ls:
        # Where the last subshell's own term was absorbed into the
        # running one, Cowan's index of a repeated term comes with it:
        # the term column is the only place left to write it.
        term = _running_ls(tail, odd=odd)
        if _absorbed(tail, ls) and tail.alpha.strip():
            term = LSTerm(mult=term.mult, l=term.l, index=tail.alpha.strip(),
                          odd=odd)
    else:
        j1 = fields[-2].run_j if len(fields) > 1 else 0
        term = PairTerm(j1=j1, j2=tail.j, odd=odd)
    return tokens, term


def _k_scheme(coupling, state, shells, open_shells, odd, drop, outer):
    """LSLK, LSJK, JJJK and the two LSJL schemes.

    The chain up to the second-to-last subshell is the LS or J1J2 chain;
    the last subshell joins it through K, and K, J1 and the total L are
    the two or three numbers the section's own FORMAT prints.
    """
    if outer is None:
        raise DesignateError(
            '%s carries its K in the OUTG11 row, so the row has to be '
            'passed in' % coupling)
    n = len(open_shells)
    if n > 1:
        tokens, _ = _chain(coupling, state, shells, open_shells, odd, drop,
                           upto=n - 1, j1=outer.j1)
    else:
        # One open subshell: there is no parent for K to be built on, so
        # the configuration is bare and the term carries the whole label.
        tokens = [Subshell(n=s.n, l=s.l, occ=s.occ)
                  for s in _kept(shells, drop)]
    mult = state.fields[-1].mult            # 2S+1 of the last subshell
    return tokens, KTerm(mult=mult, k=outer.k, odd=odd, l=outer.l)


def _true_jj(state, shells, open_shells, odd, drop):
    """The true jj chain: the groups are the j-halves, not the subshells.

    The last subshell is left whole where it has one term -- one
    electron or one hole -- because the final term's J2 already says
    which half it sits in, so writing `5d-4.(0).5d+5` in place of `5d9`
    says nothing (`true-jj_keys.txt` states the rule; the author settled
    the notation on 2026-08-14).  In the *middle* of a chain nothing
    fixes it, so `2s.2p.(1).3p` would be ambiguous and the half is
    written: `2s.2p+.(1).3p`.  An s subshell is never split, having one
    j.  GENPRT splits every subshell so that its table reads alike, and
    what does not need splitting is put back together here.
    """
    fields = list(state.fields)
    last = len(fields)
    by_slot = {}
    for k, s in enumerate(open_shells, start=1):
        field = fields[k - 1]
        parts = []
        whole = (one_term(s.l, s.occ) if k == last
                 else one_state(s.l, s.occ))
        if whole:
            shown = [Subshell(n=s.n, l=s.l, occ=s.occ)]
        else:
            shown = [Subshell(n=s.n, l=p.l, occ=p.nel, sign=p.sign,
                              v=p.seniority or None)
                     for p in field.parts if p.nel]
            if len(shown) > 1 and k < last:
                # Two halves coupling to the subshell's own J: without it
                # the halves' J values have nothing to couple through.
                parts.append(Intermediate(JTerm(field.j2 / 2.0)))
        if 2 <= k < last:
            parts.append(Intermediate(JTerm(field.k2 / 2.0)))
        by_slot[id(s)] = (shown, parts)

    tokens = []
    for s in _kept(shells, drop):
        shown, parts = by_slot.get(id(s), (None, ()))
        if shown is None:
            shown = [Subshell(n=s.n, l=s.l, occ=s.occ)]
        tokens.extend(shown)
        tokens.extend(parts)

    tail = fields[-1]
    k1 = fields[-2].k2 / 2.0 if len(fields) > 1 else 0
    return tokens, PairTerm(j1=k1, j2=tail.j2 / 2.0, odd=odd)


def _ordered(tokens):
    """Summation order, or (n, l) order where no intermediate is written."""
    if any(isinstance(t, Intermediate) for t in tokens):
        shells = [t for t in tokens if isinstance(t, Subshell)]
        keys = [(t.n, t.lval) for t in shells]
        return tokens, keys == sorted(keys)
    return sorted(tokens, key=lambda t: (t.n, t.lval)), True


# --------------------------------------------------------------------------
# the tail of a JK row
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Outer:
    """The numbers a JK section prints that no genealogy table carries."""
    j1: object = None               # running J of all but the last subshell
    k: object = None
    l: object = None                # LK's total L, a letter, else None
    mult: int = None                # where the FORMAT prints one


def read_outer(coupling, desc):
    """Decode the tail of an eigenvector row of a JK section.

    The four FORMATs are at `FOR/rcg11kd.f:1799-1817`; `desc` is what
    `outg11.Row` keeps, which starts at the column after the `=`.

      23  A6,A5,F4.1,1HK,F4.1      JJJK and LSJK
      24  2A6,2X,A1,F4.1           LSLK
      26  A6,A5,F4.1,A1,F3.1,I2    LSJLKS
      27  A6,A5,F4.1,I2,A1,F3.1    LSJLSJ
    """
    d = desc.ljust(28)
    try:
        if coupling in (JJJK, LSJK):
            return Outer(j1=float(d[11:15]), k=float(d[16:20]))
        if coupling == LSLK:
            return Outer(l=d[14:15], k=float(d[15:19]))
        if coupling == LSJLKS:
            return Outer(j1=float(d[11:15]), l=d[15:16],
                         k=float(d[16:19]), mult=int(d[19:21]))
        if coupling == LSJLSJ:
            return Outer(j1=float(d[11:15]), mult=int(d[15:17]),
                         l=d[17:18], k=float(d[18:21]))
    except ValueError:
        raise DesignateError('a %s row will not read: %r' % (coupling, desc))
    raise DesignateError('%s prints no K' % coupling)
