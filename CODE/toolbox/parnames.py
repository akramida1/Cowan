"""Turning RCE's positional parameter names into readable ones.

RCE names a parameter by the ING11 *slot* it belongs to: `F2(11)` is
F2 of the first subshell with itself, `ZETA 3` the spin-orbit parameter
of the third, `G1(13)` the exchange integral between the first and the
third, and `120D1112` a configuration-interaction integral between
configurations 1 and 2.  None of that can be read without the deck, and
the deck is where the slots are tied to their `nl` (`ing11.match_shells`).

Three names carry no slot at all -- `ALPHA`, `BETA`, `GAMMA` and the
`Tn(D 5)` family.  They are not ambiguous, though: PFGD emits them
inside the same `DO I=1,NOSUBC` loop that emits `Fk(II)` for the same
subshell (`FOR/rcg11kd.f:4383-4400`, the loop; `:4505` and `:4663`, the
two names), so the slot is the one named by the `Fk(ii)` most recently
seen.  That is exact, and it replaces the cursor-and-rescan heuristic
in `conv_cowan.pl:1676-1830`, which had to guess because it never
looked at the emission order.

This module says only what a parameter *is*.  Formatting it for a
report, and pairing it with an uncertainty, belongs to the caller.
"""

import re
from dataclasses import dataclass, replace

from .in36 import LBCD

# The names RCE writes, one pattern each.
_FG_SELF = re.compile(r'^([FG])(\d)\((\d)\3\)$')        # F2(11)
_FG_PAIR = re.compile(r'^([FG])(\d)\((\d)(\d)\)$')      # G1(13)
_ZETA = re.compile(r'^ZETA\s*(\d)$')                    # ZETA 3
_T_SHELL = re.compile(r'^T(\d?)\(([A-Z])\s*(\d)\)$')    # T2(D 5)
_ABG = re.compile(r'^(ALPHA|BETA|GAMMA)$')
_EAV = re.compile(r'^EAV\s*(\S*)$')
# A CI integral: two configuration codes, the rank, D or E, four slots.
_CI = re.compile(r'^(..)(\d)([DE])(\d)(\d)(\d)(\d)$')

_ABG_LETTER = {'ALPHA': 'A', 'BETA': 'B', 'GAMMA': 'G'}


class ParNameError(Exception):
    """A parameter name that does not belong to any subshell of the deck."""


@dataclass(frozen=True)
class Named:
    """One parameter, with its slots resolved to subshells.

    `kind` is what the parameter is -- 'Eav', 'F', 'G', 'ZETA', 'T',
    'ALPHA', 'BETA', 'GAMMA', 'R' -- and `shells` the subshells it
    involves, as `ing11.SlotShell`.  `text` is the readable name and
    `key` the order the report puts it in.
    """

    raw: str                        # exactly what RCE wrote
    kind: str
    rank: int = None
    shells: tuple = ()
    configurations: tuple = ()      # labels, for Eav and for a CI integral
    exchange: bool = False          # a CI integral's D/E
    key: str = ''

    @property
    def text(self):
        names = [s.nl for s in self.shells]
        if self.kind == 'Eav':
            return 'Eav(%s)' % (self.configurations[0] if self.configurations
                                else '?')
        if self.kind in ('F', 'G'):
            return '%s%d(%s,%s)' % (self.kind, self.rank, names[0], names[1])
        if self.kind == 'ZETA':
            # Lower case, as the author's reports have always spelled it:
            # the parameter is the Greek zeta, and RCE's upper-case name
            # is a punched-card limitation rather than a choice
            # (`conv_cowan.pl:2028`).
            return 'zeta(%s)' % names[0]
        if self.kind == 'T':
            return 'T%s(%s)' % (self.rank if self.rank else '', names[0])
        if self.kind in ('ALPHA', 'BETA', 'GAMMA'):
            return '%s(%s)' % (self.kind, names[0])
        if self.kind == 'R':
            return 'R%s%d(%s%s,%s%s)' % ('E' if self.exchange else 'D',
                                         self.rank, names[0], names[1],
                                         names[2], names[3])
        return self.raw

    def __str__(self):
        return self.text


def resolve(parameters, cards, shells):
    """Name every parameter of one parity set.

    `parameters` is that set's list in RCE's own order -- which is the
    order PFGD emitted them, configuration by configuration, and is what
    makes the slotless names resolvable.  `cards` is the parity set's
    `ing11.ConfigCard` list in card order, and `shells` maps a
    configuration label to its `ing11.SlotShell` list.

    Returns a list of `Named`, one per parameter, in the same order.
    """
    by_slot = {}
    for card in cards:
        by_slot[card.label] = {s.slot: s for s in shells[card.label]}
    order = [card.label for card in cards]

    out = []
    current = None                  # the configuration being described
    anchor = None                   # the slot of the last Fk(ii) seen
    for parameter in parameters:
        name = parameter.name.strip()
        eav = _EAV.match(name)
        if eav:
            current = _configuration_of(eav.group(1), order, parameter)
            anchor = None
            out.append(Named(raw=name, kind='Eav', configurations=(current,),
                             key=_key(order.index(current), 0, '0')))
            continue
        ci = _CI.match(name)
        if ci and (parameter.configurations or _ci_codes(ci, order)):
            out.append(_ci(name, ci, parameter, by_slot, order))
            continue
        if current is None:
            raise ParNameError(
                'parameter %r comes before any EAV, so there is no '
                'configuration to attach it to' % name)
        named, anchor = _single(name, current, by_slot[current], anchor,
                                order.index(current))
        # Every single-configuration parameter belongs to the
        # configuration of the EAV above it -- that is what makes the
        # slot numbers in its name mean anything -- so it is recorded
        # here rather than left for a caller to work out again.
        out.append(replace(named, configurations=(current,)))
    return out


def _configuration_of(label, order, parameter):
    """The configuration an `EAV <label>` names.

    RCE truncates the label to what fits its field, so a prefix match is
    what the Perl does too (`s11.pl:196`); an ambiguous prefix is an
    error rather than a pick.
    """
    if label in order:
        return label
    hits = [c for c in order if c.startswith(label)]
    if len(hits) == 1:
        return hits[0]
    raise ParNameError(
        'EAV %r matches %s of the deck\'s configurations (%s)'
        % (label, 'none' if not hits else '%d' % len(hits), ', '.join(order)))


def _single(name, label, slots, anchor, nconf):
    """One single-configuration parameter; returns it and the new anchor."""
    m = _FG_SELF.match(name)
    if m:
        letter, rank, slot = m.group(1), int(m.group(2)), int(m.group(3))
        shell = _slot(slots, slot, name, label)
        return (Named(raw=name, kind=letter, rank=rank, shells=(shell, shell),
                      key=_key(nconf, slot, '1%s%d' % (letter, rank))),
                slot)

    m = _FG_PAIR.match(name)
    if m:
        letter, rank = m.group(1), int(m.group(2))
        a, b = sorted((int(m.group(3)), int(m.group(4))))
        shells = _canonical((_slot(slots, a, name, label),
                             _slot(slots, b, name, label)))
        return (Named(raw=name, kind=letter, rank=rank, shells=shells,
                      key=_key(nconf, a, '3%d%s%d' % (b, letter, rank))),
                anchor)

    m = _ZETA.match(name)
    if m:
        slot = int(m.group(1))
        return (Named(raw=name, kind='ZETA', shells=(_slot(slots, slot, name,
                                                           label),),
                      key=_key(nconf, slot, '2')),
                anchor)

    m = _T_SHELL.match(name)
    if m:
        rank = int(m.group(1)) if m.group(1) else 0
        slot = _slot_of_shell(slots, m.group(2), int(m.group(3)), anchor,
                              name, label)
        return (Named(raw=name, kind='T', rank=rank,
                      shells=(slots[slot],),
                      key=_key(nconf, slot, '1T%d' % rank)),
                anchor)

    m = _ABG.match(name)
    if m:
        if anchor is None:
            raise ParNameError(
                '%s of configuration %s comes before any F integral, so '
                'there is no subshell to attach it to' % (name, label))
        kind = m.group(1)
        return (Named(raw=name, kind=kind, shells=(slots[anchor],),
                      key=_key(nconf, anchor, '1H%s' % _ABG_LETTER[kind])),
                anchor)

    raise ParNameError('%r (configuration %s) is not a parameter name this '
                       'reader knows' % (name, label))


def _ci(name, m, parameter, by_slot, order):
    """A configuration-interaction integral: `120D1112`.

    The two configurations come from RCEOUT's own columns where the
    file has them, rather than from the two code characters, which are
    alphanumeric once a run has more than nine configurations
    (`FOR/rcg11kd.f:6308-6322`).  OUTGINE has no such columns, so there
    the characters are decoded instead -- they are the configurations'
    serial numbers within the parity set, which is what `order` holds.
    """
    rank, letter = int(m.group(2)), m.group(3)
    n1, n2, n3, n4 = (int(m.group(i)) for i in range(4, 8))
    if parameter.configurations:
        left, right = parameter.configurations
        left = _configuration_of(left, order, parameter)
        right = _configuration_of(right, order, parameter)
    else:
        left, right = _ci_codes(m, order)
    a, b = sorted((n1, n2))
    c, d = sorted((n3, n4))
    shells = _canonical((_slot(by_slot[left], a, name, left),
                         _slot(by_slot[left], b, name, left)))         + _canonical((_slot(by_slot[right], c, name, right),
                      _slot(by_slot[right], d, name, right)))
    key = '%03d%03d%d%d%d%d%s%d' % (order.index(left), order.index(right),
                                    a, b, c, d, letter, rank)
    return Named(raw=name, kind='R', rank=rank, shells=shells,
                 configurations=(left, right), exchange=letter == 'E',
                 key=key)


def _ci_codes(m, order):
    """The two configurations a CI name encodes, or None if it cannot.

    RCG writes the serial number of each configuration as one printable
    character: `0`-`9` for 0-9, `A`-`Z` for 10-35, and one character per
    number from 36 up, with `*` standing for everything past 188
    (`FOR/rcg11kd.f:6308-6322`).  A `*`, or a number the parity set does
    not have, means the name cannot be decoded, and the caller then
    treats it as the single-configuration parameter it is not -- which
    raises rather than mislabels.
    """
    out = []
    for character in m.group(1):
        if character == '*':
            return None
        if character.isdigit():
            number = ord(character) - 48
        elif 'A' <= character <= 'Z':
            number = ord(character) - 55
        else:
            number = ord(character) - 61
        if not 1 <= number <= len(order):
            return None
        out.append(order[number - 1])
    return tuple(out)


def _canonical(shells):
    """A pair of subshells in n,l order, whatever order the deck sums them.

    The slot order is the deck's business and a reordered deck permutes
    it deliberately; a parameter's *name* should read the same either
    way, so it is written in the canonical order.
    """
    return tuple(sorted(shells, key=lambda s: (s.n, s.l)))


def _slot(slots, slot, name, label):
    shell = slots.get(slot)
    if shell is None:
        raise ParNameError(
            '%s names slot %d of configuration %s, which holds no electrons'
            % (name, slot, label))
    return shell


def _slot_of_shell(slots, symbol, occupation, anchor, name, label):
    """The slot a `Tn(D 5)` belongs to: the shell symbol and occupation.

    A T parameter names its subshell rather than its slot, so it is
    resolvable on its own; the anchor only breaks a tie where the same
    (symbol, occupation) appears twice in one configuration, and there
    the emission order says it is the one being described.
    """
    letter = symbol.lower()
    hits = [slot for slot, shell in sorted(slots.items())
            if LBCD[shell.l] == letter and shell.occ == occupation]
    if not hits:
        raise ParNameError(
            '%s names a %s%d subshell, which configuration %s does not have'
            % (name, letter, occupation, label))
    if anchor in hits:
        return anchor
    return hits[0]


def _key(nconf, slot, tail):
    """The order a report lists parameters in.

    Configuration first, then subshell, then Cowan's own grouping within
    a subshell: the F integrals, ALPHA/BETA/GAMMA and T, then ZETA, then
    the two-subshell integrals.  It is a string so that the caller can
    sort without knowing any of this.
    """
    return '%03d%02d%s' % (nconf, slot, tail)
