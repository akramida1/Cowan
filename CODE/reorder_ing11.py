#!/usr/bin/env python3
"""Write a copy of ING11 whose subshell slots are in another order.

The Python replacement for `reorder_ing11.pl`.  The slot order of a
configuration card is the order RCG *sums* the subshells in, so a deck
summed 2s, 2p, 1s gives the same levels as one summed 1s, 2s, 2p but
builds them through different intermediate terms -- which is what a
reordered deck is for.  `WORK/U90/reordered` is `WORK/U90` reordered
this way.

Usage:

    reorder_ing11 <output file> <slot order>

The slot order is the old slot numbers written in their new order, so

    reorder_ing11 ing11.new 231

builds a deck whose first slot is the old second, whose second is the
old third and whose third is the old first.  Slots the string does not
name keep their places.

ING11 and OUTG11 are read from the working directory.  OUTG11 is needed
because a parameter card carries no names -- its values are identified
by position alone -- and RCG prints the names nowhere else; it must
therefore be the OUTG11 of this ING11.  The reasoning behind the move
is in `toolbox/reorder.py`.

What this does that the Perl did not:

* it needs no IN36 and no principal quantum numbers, because the
  permutation is on slot positions and so is the D/E pairing of a
  configuration-interaction integral.  The Perl's restriction against
  moving the last occupied subshell, which existed only because it had
  to re-attach `n`, is gone with it;
* values are moved as the ten characters RCN2 punched, never re-encoded,
  so no field can lose a digit in transit.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, reorder                              # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


def main(argv):
    if len(argv) != 2:
        raise Usage('two arguments are needed, an output file and a slot '
                    'order')
    out_file, order = argv
    new = reorder.read_order(order)
    lines = reorder.reorder_deck(files.ING11, new, folder='.')
    ending = reorder.line_ending(files.ING11, folder='.')
    with open(out_file, 'w', newline='') as fh:
        fh.write(''.join(line + ending for line in lines))
    moved = sum(1 for slot, to in new.items() if slot != to)
    print('%s written: %d of the %d slots moved.'
          % (out_file, moved, len(new)))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('reorder_ing11: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (reorder.ReorderError, OSError) as problem:
        sys.stderr.write('reorder_ing11: %s\n' % problem)
        sys.exit(1)
