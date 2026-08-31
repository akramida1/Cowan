#!/usr/bin/env python3
"""Add a fixed amount to every Eav parameter of ING11.

The Python replacement for `add_eav.pl`.  Eav is the average energy of a
configuration, so shifting all of them together moves the whole
calculation up or down the energy scale without changing a single level
separation -- which is what putting a calculation on an observed ground
level, or on an ionization limit, amounts to.

Usage:

    add_eav <shift> [folder]

`<shift>` is in the run's own energy unit, which is the UENRGY of the
ING11 main card and is 1000 cm-1 in every deck here.  It may be
negative.  ING11 is read from `folder`, or from the working directory,
and is rewritten in place with the previous deck kept as `ING11.BAK`.

Only the Eav on a **parameter card** is shifted, as in the Perl.  The
average energy in columns 59-71 of a configuration card is RCN2's
Hartree-Fock value, which RCG does not use once a parameter card carries
an Eav, and it is left alone.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, ing11                                # noqa: E402

BACKUP = 'ING11.BAK'


class Usage(Exception):
    """The command line will not do what was asked."""


def shift_eav(deck, amount):
    """Add `amount` to the Eav of every parameter card; count them.

    Eav is the first value of a parameter card and is the only one a CI
    card does not have, so the CI cards are untouched by construction
    rather than by a test on their text.
    """
    changed = 0
    for cards in deck.parameters:
        for card in cards:
            if not card.values:
                continue
            values = list(card.values)
            values[0] += amount
            ing11.set_values(card, values)
            changed += 1
    return changed


def main(argv):
    if not argv or len(argv) > 2:
        raise Usage('one or two arguments are needed, a shift and '
                    'optionally a folder')
    try:
        amount = float(argv[0])
    except ValueError:
        raise Usage('%r is not a number' % argv[0])
    if amount == 0.0:
        raise Usage('a shift of zero would not change anything')
    folder = argv[1] if len(argv) > 1 else '.'

    path = files.resolve(files.ING11, folder)
    ending = files.line_ending(path)
    deck = ing11.read_ing11(path)
    changed = shift_eav(deck, amount)
    files.write_deck(path, ing11.write_ing11(deck), ending, backup=BACKUP)
    print('%s: %g added to %d Eav parameters; the old deck is %s.'
          % (path, amount, changed, BACKUP))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('add_eav: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (ing11.Ing11Error, OSError) as problem:
        sys.stderr.write('add_eav: %s\n' % problem)
        sys.exit(1)
