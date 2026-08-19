"""Permuting the subshell slots of an ING11 deck.

A configuration card names its subshells by slot: eight `(A1,I2,2X)`
fields in columns 1-40 (`FOR/rcg11kd.f:3334`).  The slot order is the
order RCG *sums* the subshells in, and it is the order every coupling
chain is built in, so summing them in another order gives the same
levels with different intermediate terms -- which is the whole reason
for reordering a deck.  `WORK/U90/reordered` is `WORK/U90` summed 2s,
2p, 1s instead of 1s, 2s, 2p.

Permuting the cards is the easy half.  The hard half is that a
parameter card carries no names: it is a row of `(F9.4,I1)` fields whose
*position* says what each one is, because RCG emits them in a fixed
order that depends on the slot numbers.  Move the slots and that order
changes, so the values have to move with it.

Two things make this exact rather than empirical:

**The emission order is in the Fortran.**  PFGD writes, per
configuration (`FOR/rcg11kd.f:4300-5000`):

    EAV
    for each slot i:      Fk(ii) for k ascending, then ALPHA, BETA,
                          GAMMA and the T family, which PFGD emits
                          inside the same `DO 199 I=1,NOSUBC` loop and
                          which therefore belong to the slot the last
                          Fk(ii) named
    for each slot i:      ZETA i                      (`DO 299 I=1,NOSUBC`)
    for i < j, k ascending:  Fk(ij)   (`DO 399 I=1,NOSCM1; DO 389 J=I+1`)
    for i < j, k ascending:  Gk(ij)   (`DO 497 I=1,NOSCM1; DO 489 J=I+1`)

and, for a pair of configurations, PRK writes the CI integrals in the
order of `DO IRHO=1,NOSUBC; DO ISIG=IRHO,NOSUBC; DO IRHOP=1,NOSUBC; DO
ISIGP=IRHOP,NOSUBC` (`:5333-5339`), each group being its D integrals by
ascending rank and then its E integrals (`:5781-5790`).  So a sort key
built from the *new* slot numbers reproduces the order RCG will use for
the reordered deck, and the sort is stable so that everything sharing a
key -- the ranks of one integral, ALPHA/BETA/GAMMA/T behind their
Fk(ii) -- keeps the order it already had.

**The D/E letter is positional too.**  `120D1323` is
`(configurations 1,2) rank 0, direct, slots 1,3 of the first and 2,3 of
the second` (`:5889`), and D means the pairing (n1 with n3, n2 with n4),
E the pairing (n1 with n4, n2 with n3).  RCG always writes n1 <= n2 and
n3 <= n4, so after a permutation the two pairs may have to be sorted
back -- and sorting exactly one of them exchanges the pairings, which
turns a D into an E and back.  Where n1 = n2 or n3 = n4 the two
pairings are the same integral and RCG calls it D.

None of this needs `n`: the permutation is on slot *positions*, and the
pairing it has to preserve is between positions.  So this module reads
neither IN36 nor the subshell identities -- which is also why it has no
counterpart of `reorder_ing11.pl`'s "last movable shell" restriction,
which existed because that script did have to re-attach `n`.

The names themselves come from OUTG11 (`outg11.read_parameter_names`),
which is the only file that carries them.

Values are moved as **raw ten-character field text**, never re-encoded.
RCN2 punches an ING11 field in whichever way it fits -- `822465.00` with
its point, `140570` without, `0.00e+01` in exponential form -- and
re-encoding is both lossy and unnecessary when the field is only
changing places.
"""

import re

from . import files, ing11, outg11

# Where a value field sits, by card.  A parameter card spends columns
# 21-30 on EAV as F10.5 and then has four `(F9.4,I1)` pairs (`:10790`);
# a CI card has no EAV and five pairs (`:10902`).  A continuation of
# either is seven pairs from column 1 (`:10834`).
_FIRST = {ing11.ParameterCard: 30, ing11.CiCard: 20}
_ON_FIRST = {ing11.ParameterCard: 4, ing11.CiCard: 5}
_ON_CONT = 7
_WIDTH = 10

# The names RCE and RCG write, one pattern each; see `parnames.py`,
# which reads the same names for a different purpose.
_EAV = re.compile(r'^EAV\b')
_FG = re.compile(r'^([FG])(\d)\((\d)(\d)\)$')
_ZETA = re.compile(r'^ZETA\s*(\d)$')
_ANCHORED = re.compile(r'^(ALPHA|BETA|GAMMA|T\d?\(.*\))$')
_CI = re.compile(r'^(..)(\d)([DE])(\d)(\d)(\d)(\d)$')


class ReorderError(Exception):
    """A deck or a permutation this module will not act on."""


# --------------------------------------------------------------------------
# the permutation
# --------------------------------------------------------------------------

def read_order(text):
    """The slot permutation behind a digit string such as `231`.

    The digits are the old slot numbers in their new order, which is
    how `reorder_ing11.pl` has always read them: `231` builds a deck
    whose first slot is the old second, whose second is the old third,
    and whose third is the old first.

    Returns `new[old] = new slot`, indexed 1 to 8.  Slots past the end
    of the string keep their places, so a partial permutation is
    allowed and leaves the tail of every card alone.
    """
    digits = [c for c in text if not c.isspace()]
    if len(digits) < 2:
        raise ReorderError(
            'a slot order needs at least two digits; %r has %d'
            % (text, len(digits)))
    try:
        old = [int(c) for c in digits]
    except ValueError:
        raise ReorderError('%r is not a string of slot numbers' % text)
    if sorted(old) != list(range(1, len(old) + 1)):
        raise ReorderError(
            '%r is not a permutation: it has to use each of the numbers 1 to '
            '%d exactly once' % (text, len(digits)))
    if len(old) > ing11.NOSUBC_MAX:
        raise ReorderError('a configuration card has only %d slots, and %r '
                           'names %d' % (ing11.NOSUBC_MAX, text, len(old)))
    new = {slot: slot for slot in range(1, ing11.NOSUBC_MAX + 1)}
    for position, slot in enumerate(old, start=1):
        new[slot] = position
    return new


# --------------------------------------------------------------------------
# the emission order
# --------------------------------------------------------------------------

def _single_key(name, new, anchor):
    """The sort key of one single-configuration parameter, and the anchor.

    The anchor is the slot of the most recent `Fk(ii)`, which is the
    subshell ALPHA, BETA, GAMMA and T belong to; it is carried along
    rather than looked up because that is exactly how PFGD emits them.
    """
    if _EAV.match(name):
        return (0,), anchor

    m = _FG.match(name)
    if m:
        letter, i, j = m.group(1), int(m.group(3)), int(m.group(4))
        if i == j:
            return (1, new[i]), new[i]
        a, b = sorted((new[i], new[j]))
        return (3 if letter == 'F' else 4, a, b), anchor

    m = _ZETA.match(name)
    if m:
        return (2, new[int(m.group(1))]), anchor

    if _ANCHORED.match(name):
        if anchor is None:
            raise ReorderError(
                '%s comes before any F integral, so there is no subshell to '
                'move it with' % name)
        return (1, anchor), anchor

    raise ReorderError('%r is not a parameter name this reader knows' % name)


def _ci_key(name, new):
    """The sort key of one configuration-interaction integral."""
    m = _CI.match(name)
    if not m:
        raise ReorderError(
            '%r is not a configuration-interaction parameter name' % name)
    letter = m.group(3)
    a1, a2 = new[int(m.group(4))], new[int(m.group(5))]
    b1, b2 = new[int(m.group(6))], new[int(m.group(7))]
    n1, n2 = sorted((a1, a2))
    n3, n4 = sorted((b1, b2))
    if n1 == n2 or n3 == n4:
        # Both pairings are the same integral, and RCG writes it as D.
        letter = 'D'
    elif (a1 > a2) != (b1 > b2):
        # Exactly one of the two pairs had to be sorted back, which
        # exchanges the pairings.
        letter = 'E' if letter == 'D' else 'D'
    return (n1, n2, n3, n4, 0 if letter == 'D' else 1)


def rename(name, new):
    """What RCG will call this parameter once the slots have moved.

    The writer does not need this -- an ING11 card carries no names --
    but a reordered deck that has been run through RCG can be checked
    against it, name for name, which is how the rules above are held to
    RCG's own output rather than to a reading of the Fortran.
    """
    if _EAV.match(name) or _ANCHORED.match(name):
        return name
    m = _FG.match(name)
    if m:
        a, b = sorted((new[int(m.group(3))], new[int(m.group(4))]))
        return '%s%s(%d%d)' % (m.group(1), m.group(2), a, b)
    m = _ZETA.match(name)
    if m:
        return 'ZETA %d' % new[int(m.group(1))]
    m = _CI.match(name)
    if m:
        letter = m.group(3)
        a1, a2 = new[int(m.group(4))], new[int(m.group(5))]
        b1, b2 = new[int(m.group(6))], new[int(m.group(7))]
        n1, n2 = sorted((a1, a2))
        n3, n4 = sorted((b1, b2))
        if n1 == n2 or n3 == n4:
            letter = 'D'
        elif (a1 > a2) != (b1 > b2):
            letter = 'E' if letter == 'D' else 'D'
        return '%s%s%s%d%d%d%d' % (m.group(1), m.group(2), letter,
                                   n1, n2, n3, n4)
    raise ReorderError('%r is not a parameter name this reader knows' % name)


def permutation(names, new, is_ci):
    """Where each of a card's value fields goes, as `[old index, ...]`.

    The result is in the *new* order and holds the old position of each
    field, which is what a caller needs to move the raw text.  EAV is
    not in it: it stays where it is, in the F10.5 field a parameter
    card keeps for it.
    """
    keys, anchor = [], None
    for at, name in enumerate(names):
        if is_ci:
            key = _ci_key(name, new)
        else:
            key, anchor = _single_key(name, new, anchor)
        keys.append((key, at))
    keys.sort(key=lambda pair: pair[0])         # stable: ties keep `at` order
    return [at for _, at in keys]


# --------------------------------------------------------------------------
# rewriting the deck
# --------------------------------------------------------------------------

def _slot_fields(card_text):
    """The eight `(A1,I2,2X)` slot fields of a configuration card."""
    text = card_text.ljust(40)
    return [text[5 * i:5 * i + 5] for i in range(ing11.NOSUBC_MAX)]


def reorder_config(card_text, new):
    """One configuration card with its slots permuted.

    Columns 41 onward -- the 3A6 name, the average energy and the
    kinetic energy -- are copied untouched: nothing there depends on
    which slot a subshell sits in.
    """
    fields = _slot_fields(card_text)
    old_of_new = {position: slot for slot, position in new.items()}
    moved = [fields[old_of_new[p] - 1] for p in range(1, ing11.NOSUBC_MAX + 1)]
    return ''.join(moved) + card_text.ljust(40)[40:]


def _field_places(card, count):
    """`[(line, column), ...]` for a card's first `count` value fields."""
    start = _FIRST[type(card)]
    on_first = _ON_FIRST[type(card)]
    places = []
    for i in range(min(count, on_first)):
        places.append((0, start + _WIDTH * i))
    line, left = 1, count - on_first
    while left > 0:
        for i in range(min(left, _ON_CONT)):
            places.append((line, _WIDTH * i))
        left -= _ON_CONT
        line += 1
    return places


def reorder_card(card, names, new, lines=None):
    """The lines of one parameter or CI card, with its values permuted.

    The card's own layout is not disturbed: the same number of values
    sits in the same fields on the same number of lines, and only their
    order changes, so every value is carried across as the ten
    characters RCN2 punched.

    `lines` is the card as the file holds it, trailing blanks and all;
    without it the reader's own copy is used, which has them off.
    """
    is_ci = isinstance(card, ing11.CiCard)
    if len(names) != len(card.values):
        raise ReorderError(
            'OUTG11 gives %d name(s) for %s but its ING11 card carries %d '
            'value(s)' % (len(names), card.name.strip(), len(card.values)))

    # A parameter card's first value is EAV, which is not one of the
    # `(F9.4,I1)` fields and never moves.
    movable = names if is_ci else names[1:]
    order = permutation(movable, new, is_ci)

    places = _field_places(card, len(movable))
    source = list(lines if lines is not None else card.raw)
    # Each line is padded only as far as its own last field, because
    # RCN2 stops a continuation line after its last value and a copy
    # that blank-fills it to column 70 differs from the original on
    # every one of those lines.
    needed = {}
    for line, at in places:
        needed[line] = max(needed.get(line, 0), at + _WIDTH)
    out = [text.ljust(needed.get(i, 0)) for i, text in enumerate(source)]
    was = [out[line][at:at + _WIDTH] for line, at in places]
    for (line, at), came_from in zip(places, order):
        out[line] = (out[line][:at] + was[came_from]
                     + out[line][at + _WIDTH:])
    return out


# --------------------------------------------------------------------------
# the whole deck
# --------------------------------------------------------------------------

def line_ending(path=files.ING11, folder=None):
    """The line terminator a deck uses; `files.line_ending` does it.

    Kept as a name here because `reorder_ing11.py` and `test_reorder.py`
    call it, and because every deck rewriter in the package needs it --
    which is why it now lives in `files.py` with the rest of the
    file-level work.
    """
    return files.line_ending(path, folder)


def reorder_deck(path, new, outg11_path=None, folder=None):
    """Read an ING11 and return the lines of its reordered copy.

    `folder`, when given, is where both ING11 and OUTG11 are looked up,
    case-insensitively.  Everything this module does not rewrite -- the
    prelude, the main card, the radial-integral cards and Cowan's
    sentinels -- is copied out of the file byte for byte, trailing
    blanks and all, which is why the file is read again here instead of
    being rebuilt from the deck object.
    """
    if folder is not None:
        path = files.resolve(path, folder)
        outg11_path = files.resolve(outg11_path or files.OUTG11, folder)
    deck = ing11.read_ing11(path)
    sections = outg11.read_parameter_names(outg11_path)

    with open(path, 'r', errors='replace') as fh:
        raw = fh.read().split('\n')
    if raw and raw[-1] == '':
        raw.pop()
    raw = [line.rstrip('\r') for line in raw]

    out = list(raw)
    at = len(deck.prelude) + 1                  # the prelude and the main card
    for card in deck.configurations:
        out[at] = reorder_config(raw[at], new)
        at += 1
    for k in (1, 2):
        if not deck.main.nconfigurations(k):
            continue
        for card in deck.parameters[k - 1] + deck.ci[k - 1]:
            names = sections.get(card.name.strip())
            if names is None:
                raise ReorderError(
                    'OUTG11 has no parameter section for %s, so there is no '
                    'way to tell what its values are.  Is it the OUTG11 of '
                    'this ING11?' % card.name.strip())
            n = len(card.raw)
            lines = reorder_card(card, names.names, new, raw[at:at + n])
            for offset, line in enumerate(lines):
                out[at + offset] = line
            at += n
    return out
