#!/usr/bin/env python3
"""Put RCE's fitted parameters back into ING11.

The Python replacement for `s11.pl`, and the step that closes the loop
of a least-squares fit: RCE fits the Slater parameters to the observed
levels and prints them in RCEOUT, and RCG has to be re-run with them
before the eigenvectors, the transition probabilities and everything
downstream describe the fitted calculation rather than the Hartree-Fock
one it started from.

Usage:

    s11 [2|auto] [folder]

With no argument only the first parity set is substituted, which is
what the Perl did; `2` does both, and `auto` -- what `update11` asks
for -- takes the number of parity sets from IN36, so that nothing has
to be typed to match the deck.  ING11 and RCEOUT are read from
`folder`, or from the working directory, and ING11 is rewritten in
place with the previous deck kept as `ing11.bak`.

A parameter card carries no names -- its values are read by position --
so the two files are matched by **position and count**: RCE prints its
parameters in the order RCG asked for them, and each ING11 card says in
columns 19-20 how many of them are its own.  The configuration name of
each card is checked against the one RCE printed beside the parameter
before anything is written, and a disagreement stops the run rather
than shifting the whole deck by one.

Two differences from `s11.pl`, both in how a value is punched:

* **a CI parameter no longer overflows its field.**  `s11.pl:228`
  writes `%9.4f` into nine columns, which takes ten characters for any
  value at or past 1000 in magnitude and pushes the rest of the card
  one column right.  `WORK/U90`'s `131D1233 = -1369.505304` is such a
  value, and the corrupted card is what the author reported.  The field
  is written through `fortfmt.write_f_wide`, which gives up a decimal
  rather than the number -- `-1369.505` -- which is what the author's
  own deck holds.
* **a parameter that is zero is punched as RCN2 punches it.**  RCN2
  writes `1PE9.2` (`FOR/RCN2K.F:1946`), giving ` 0.00E+00`; `s11.pl:231`
  writes `%10.2e` and then deletes a zero out of the exponent to make it
  fit, giving `  0.00e+0`.  RCG reads both as zero, so every deck in
  `WORK/` that has been through the Perl carries the second form and
  will change on the first run of this program.  The Fortran's own
  encoder is used because it is the specification, and because the
  Perl's squeeze fails outright for an exponent of -10 or beyond, where
  there is no zero in it to delete.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, in36, ing11, rce                     # noqa: E402

BACKUP = 'ing11.bak'


class Usage(Exception):
    """The command line will not do what was asked."""


class SubstitutionError(Exception):
    """RCEOUT and ING11 do not describe the same calculation."""


def _eav_label(parameter):
    """The configuration RCE printed an EAV for: 'EAV 2s' -> '2s'."""
    head, _, tail = parameter.name.partition(' ')
    return tail.strip() if head.upper() == 'EAV' else None


def _agrees(deck_label, rce_label):
    """Whether RCE's name and the deck's label are the same label.

    RCE holds a configuration name in an A6 and ING11 in the last
    twelve columns of its 3A6, so a long label is truncated in RCEOUT
    and only the beginning can be compared -- the test at
    `s11.pl:196` and `:183`.
    """
    return bool(rce_label) and deck_label.startswith(rce_label)


def _take(parameters, at, count, what):
    """The next `count` fitted parameters, or a complaint that there are not."""
    if at + count > len(parameters):
        raise SubstitutionError(
            'RCEOUT runs out of parameters at %s: it has %d where the deck '
            'asks for %d' % (what, len(parameters) - at, count))
    return parameters[at:at + count]


def substitute(deck, csets, kset):
    """Replace one parity set's values with the fitted ones; count them.

    The cards are walked in the order RCG emits them -- every
    single-configuration card, then every CI card -- because that is the
    order RCE prints its parameter table in.
    """
    if kset > len(csets):
        raise SubstitutionError(
            'ING11 has a parity set %d that RCEOUT does not' % kset)
    fitted = csets[kset - 1].parameters
    at = 0

    for card in deck.parameters[kset - 1]:
        wanted = _take(fitted, at, len(card.values), 'configuration %s'
                       % card.label)
        label = _eav_label(wanted[0])
        if not _agrees(card.label, label):
            raise SubstitutionError(
                'the deck\'s next parameter card is %s but RCEOUT\'s next '
                'parameter is %r; the two files are out of step'
                % (card.label, wanted[0].name))
        ing11.set_values(card, [p.value for p in wanted])
        at += len(wanted)

    for card in deck.ci[kset - 1]:
        wanted = _take(fitted, at, len(card.values), 'CI card %s'
                       % card.name.strip())
        left, right = card.labels
        pair = wanted[0].configurations
        if len(pair) != 2 or not (_agrees(left, pair[0])
                                  and _agrees(right, pair[1])):
            raise SubstitutionError(
                'the deck\'s next CI card couples %s and %s but RCEOUT\'s '
                'next parameter %r belongs to %s; the two files are out of '
                'step' % (left, right, wanted[0].name,
                          '-'.join(pair) if pair else 'no pair'))
        ing11.set_values(card, [p.value for p in wanted])
        at += len(wanted)

    if at != len(fitted):
        raise SubstitutionError(
            'parity set %d: RCEOUT has %d parameters and the deck has room '
            'for %d' % (kset, len(fitted), at))
    return at


def parity_sets(folder):
    """How many parities the run has, from IN36.

    RCN2 splits the configurations into a first and a second set by
    parity, so the `2` that had to be typed after `s11` in the Perl
    days is nothing but a count of the parities in IN36.  The Python
    parses IN36 the way RCN itself does (`toolbox.in36`), so the count
    can be read off instead of typed.

    Returns 0 when IN36 cannot answer: a sub-folder holding only ING11
    and RCEOUT, as `WORK/U90/reordered` does, or an IN36 that is the
    bare end card and names no configuration, as `WORK/AI/DIEL`'s is
    until `diel` writes the real one.  The caller then falls back to
    substituting every set: the deck's own second set is empty when
    there is only one parity, and `main` skips it.
    """
    path = files.resolve(files.IN36, folder, required=False)
    if path is None:
        return 0
    return len(in36.parities(in36.read_in36(path).configurations))


def main(argv):
    both, auto, folder = False, False, '.'
    for argument in argv:
        if argument == '2':
            both = True
        elif argument == 'auto':
            auto = True
        elif folder == '.':
            folder = argument
        else:
            raise Usage('%r is neither the option 2 nor auto nor a folder'
                        % argument)

    if auto and not both:
        found = parity_sets(folder)
        both = found != 1
        print('IN36: one parity; substituting the first parity set.'
              if found == 1 else
              'IN36: both parities; substituting both parity sets.'
              if found else
              'IN36 does not say how many parities there are; '
              'substituting every parity set the deck has.')

    path = files.resolve(files.ING11, folder)
    ending = files.line_ending(path)
    deck = ing11.read_ing11(path)
    csets = rce.read_rceout(folder=folder)

    done = []
    for kset in ((1, 2) if both else (1,)):
        if not deck.main.nconfigurations(kset):
            continue
        done.append((kset, substitute(deck, csets, kset)))
    if not done:
        raise SubstitutionError('%s has no configuration in parity set 1'
                                % path)

    files.write_deck(path, ing11.write_ing11(deck), ending, backup=BACKUP)
    print('%s: %s substituted; the old deck is %s.'
          % (path, ', '.join('%d parameters into parity set %d' % (n, k)
                             for k, n in done), BACKUP))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('s11: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (SubstitutionError, ing11.Ing11Error, rce.RceError,
            OSError) as problem:
        sys.stderr.write('s11: %s\n' % problem)
        sys.exit(1)
