#!/usr/bin/env python3
"""Add a fixed amount to every Eav parameter of ING11.

The Python replacement for `add_eav.pl`.  Eav is the average energy of a
configuration, so shifting all of them together moves the whole
calculation up or down the energy scale without changing a single level
separation -- which is what putting a calculation on an observed ground
level, or on an ionization limit, amounts to.

Usage:

    add_eav <shift> [folder]
    add_eav --ground [folder]

`<shift>` is in the run's own energy unit, which is the UENRGY of the
ING11 main card and is 1000 cm-1 in every deck here.  It may be
negative.  ING11 is read from `folder`, or from the working directory,
and is rewritten in place with the previous deck kept as `ING11.BAK`.

Only the Eav on a **parameter card** is shifted, as in the Perl.  The
average energy in columns 59-71 of a configuration card is RCN2's
Hartree-Fock value, which RCG does not use once a parameter card carries
an Eav, and it is left alone.

`--ground` works on **OUTGINE** instead, the file the first RCE run of
a fit (INVECT=0) takes its starting parameters and its calculated levels
from.  It takes the lowest level listed there for the ground and adds
minus its energy to every Eav and every level, so that RCE starts with
the calculated ground at zero; RCG's file is kept as `OUTGINE.RCG`.
Nothing else needs to change: TAPE2E holds neither parameter values nor
eigenvalues, and a common shift of the Eav moves no eigenvector.  `rce
--ground` runs this before RCE.

  * With INVECT not 0 in OUTGINE's first record RCE reads RCEINP and not
    OUTGINE, so nothing is shifted; that is a warning, not an error.
  * A ground already at zero is left alone, so a second `--ground` does
    nothing.
  * An OUTGINE whose levels are all zero (RCG made no diagonalization)
    is an error.
  * An Eav or level at or below -4000 is RCG's mark of a continuum
    configuration in an autoionization run, not an energy, and is left
    alone, as RCG's own DELEAV shift leaves it.

ING11 and OUTG11 keep RCG's values, so `printout` shows RCG's Eav
differing from RCEOUT's until the first `update11` carries the fitted,
shifted values into ING11.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, genealogy, ing11, rce                # noqa: E402

BACKUP = 'ING11.BAK'
OUTGINE_BACKUP = 'OUTGINE.RCG'


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


def ground(folder):
    """Put the calculated ground level of OUTGINE at zero."""
    path = files.resolve(rce.OUTGINE, folder)
    lines = files.read_lines(path, strip=False)
    invect = rce.outgine_invect(lines)
    if invect != 0:
        sys.stderr.write(
            'add_eav: warning: INVECT=%d in the first record of %s, so RCE '
            'takes its parameters\n  and levels from RCEINP, not from '
            'OUTGINE.  Nothing has been shifted.\n' % (invect, path))
        return 0
    low = rce.outgine_ground(lines, path)
    print('Ground level adopted: %s, E_RCG = %s'
          % (describe(low, folder), low.text))
    if low.energy == 0.0:
        print('%s: the ground level is already at zero; nothing to shift.'
              % path)
        return 0
    new, neav, nlev = rce.shift_outgine(lines, -low.energy, path)
    files.write_deck(path, new, files.line_ending(path),
                     backup=OUTGINE_BACKUP)
    shift = low.text[1:] if low.text.startswith('-') else '-' + low.text
    print('%s: %s added to %d Eav parameters and %d levels; RCG\'s file is '
          'kept as %s.' % (path, shift, neav, nlev, OUTGINE_BACKUP))
    return 0


def describe(low, folder):
    """Where the ground level is, with its J when genealogy.txt says.

    OUTGINE lists a block's levels by position and never names its J.
    genealogy.txt names every block in the same order; it is believed
    only if its blocks of this set agree with OUTGINE's in number and,
    for the block concerned, in size -- a file left over from another
    run does not.
    """
    where = 'level %d of J block %d' % (low.index, low.block)
    try:
        path = files.resolve(files.GENEALOGY, folder, required=False)
        heads = [h for h in genealogy.read_headings(path)
                 if h[0] == low.cset] if path else []
    except (genealogy.GenealogyError, OSError):
        heads = []
    if len(heads) == low.nblocks and heads[low.block - 1][2] == low.nlevels:
        where = 'level %d with J = %g' % (low.index, heads[low.block - 1][1])
    return '%s of configuration set %d' % (where, low.cset)


def main(argv):
    if argv and argv[0] == '--ground':
        if len(argv) > 2:
            raise Usage('--ground takes at most a folder')
        return ground(argv[1] if len(argv) > 1 else '.')
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
    except (ing11.Ing11Error, rce.RceError, OSError) as problem:
        sys.stderr.write('add_eav: %s\n' % problem)
        sys.exit(1)
