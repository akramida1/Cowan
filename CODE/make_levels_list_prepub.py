#!/usr/bin/env python3
"""Write levels_prepub.txt: each configuration in the scheme chosen for it.

The middle of the three tools that turn a finished fit into a published
level list:

    levels.txt --12--> coupling_choice.txt --13--> levels_prepub.txt
                                                   --14--> levels_pub.txt

`config_purity.py` measured the purity of every configuration in every
coupling scheme and wrote `coupling_choice.txt`; the author marked one
scheme per configuration with a star.  This applies the marks.  For each
level it takes from each scheme's block only those components whose
configuration was assigned that scheme, concatenates them, sorts
descending by percentage and prints.  One composition may therefore
quote one configuration in JJ and the next in LS, and the percentages
still add up, because the weight a configuration carries in a level does
not depend on the scheme -- see `toolbox/prepub.py` for why, and for the
check that says so about this particular run rather than in general.

Usage:

    make_levels_list_prepub [--levels levels.txt]
                            [--choice coupling_choice.txt]
                            [-o levels_prepub.txt]
                            [--components 5] [--min-percent 1.5]

Run `printout --couplings ... --dominant-config --components 10
--min-percent 0.5` for this pipeline: the components come from
printout's per-block cut, and a configuration cut out of the very block
it was assigned cannot be recovered here.  The trimming for publication
is this script's `--components` and `--min-percent`, applied after the
merge.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import prepub, purity                             # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


class Options:

    def __init__(self, argv):
        self.levels = 'levels.txt'
        self.choice = purity.CHOICE
        self.out = prepub.PREPUB
        self.components = 5
        self.min_percent = 1.5

        rest = list(argv)
        positional = []
        while rest:
            arg = rest.pop(0)
            if arg == '-o':
                if not rest:
                    raise Usage('-o needs a file name')
                self.out = rest.pop(0)
                continue
            if not arg.startswith('--'):
                positional.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.replace('-', '_')
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name.replace('_', '-'))
                value = rest.pop(0)
            if name == 'levels':
                self.levels = value
            elif name == 'choice':
                self.choice = value
            elif name in ('o', 'out'):
                self.out = value
            elif name == 'components':
                self.components = _number(int, value, '--components')
            elif name == 'min_percent':
                self.min_percent = _number(float, value, '--min-percent')
            else:
                raise Usage('--%s is not an option' % name.replace('_', '-'))

        if positional:
            self.levels = positional[0]
        if len(positional) > 1:
            raise Usage('too many file names: %s' % ', '.join(positional))
        if self.components < 1:
            raise Usage('--components must be at least 1')


def _number(kind, text, name):
    try:
        return kind(text)
    except ValueError:
        raise Usage('%s wants a number, not %r' % (name, text))


def main(argv):
    options = Options(argv)
    table = purity.read_levels_txt(options.levels)
    choice, columns = purity.read_choice(options.choice)

    # The choice file's columns and the levels file's blocks are two
    # readings of one run and must agree; a scheme marked in the choice
    # file that no block carries would silently take a configuration's
    # components out of the composition altogether.
    unknown = [name for name in set(choice.values())
               if name not in table.couplings]
    if unknown:
        raise purity.PurityError(
            '%s marks the scheme%s %s, which %s has no block for.  Its '
            'blocks are %s.'
            % (options.choice, '' if len(unknown) == 1 else 's',
               ', '.join(sorted(unknown)), options.levels,
               ', '.join(table.couplings)))
    for name in columns:
        if name not in table.couplings:
            print('   %s is a column of %s and not a block of %s; nothing '
                  'is marked in it' % (name, options.choice, options.levels))

    levels, missing = prepub.collect(table)
    print('%s: %d levels, %d couplings (%s); %s: %d configurations'
          % (options.levels, len(levels), len(table.couplings),
             ', '.join(table.couplings), options.choice, len(choice)))
    for coupling, absent in sorted(missing.items()):
        print('   %d levels have no %s block and take their components from '
              'the other schemes' % (len(absent), coupling))

    deck = _deck_labels(options)
    unassigned = prepub.apply_choice(levels, choice, table.couplings,
                                     deck=deck)
    if unassigned:
        print('   %d configurations appear in a composition and are not in '
              '%s, so they are written in %s: %s'
              % (len(unassigned), options.choice, table.couplings[0],
                 ', '.join('%s (%d components)' % (label, unassigned[label])
                           for label in sorted(unassigned))))

    _check(levels, table.couplings)

    prepub.trim(levels, components=options.components,
                min_percent=options.min_percent)
    ordered = prepub.order(levels)
    prepub.write_prepub(options.out, ordered, components=options.components)
    print('wrote %s (%d levels, up to %d components at %g%% and over)'
          % (options.out, len(ordered), options.components,
             options.min_percent))
    return 0


def _deck_labels(options):
    """The deck's configuration labels, where IN36 is beside us.

    They are only used to tell the two kinds of unnamed configuration
    apart -- one the choice file legitimately omits, one that says the
    files are from different runs -- so a folder without IN36 still
    runs, with the weaker of the two checks.
    """
    if not os.path.exists('IN36'):
        print('   no IN36 here, so a component configuration missing from '
              '%s cannot be checked against the deck' % options.choice)
        return None
    return [row.label for row in purity.deck_configurations('.')]


def _check(levels, couplings):
    """Say what the invariance and the truncation did to this run."""
    off = prepub.disagreements(levels, couplings)
    if off:
        worst = sorted(off, key=lambda d: -d.spread)[:5]
        print('   %d configuration weights differ between schemes by more '
              'than the printing accounts for:' % len(off))
        for one in worst:
            print('      %s Ec=%.4f J=%.1f  %s: %s (slack %g)'
                  % (one.level.parity, one.level.e_calc, one.level.j,
                     one.label,
                     ', '.join('%s %d%%' % (c, one.weights[c])
                               for c in couplings if c in one.weights),
                     one.slack))
    short = prepub.shortfalls(levels)
    if short:
        worst = max(short, key=lambda triple: triple[2] - triple[1])
        print('   %d levels come out of the merge covering less of their '
              'weight than one block alone did, the worst being %d%% against '
              '%d%% at Ec=%.4f J=%.1f.  Nothing is wrong: the scheme a '
              'configuration was assigned spreads it over more basis states '
              'than another scheme does, so more of it falls past the '
              'per-block cut of printout.  Raise printout --components if '
              'that weight is wanted.'
              % (len(short), worst[1], worst[2], worst[0].e_calc, worst[0].j))


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('make_levels_list_prepub: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except purity.PurityError as problem:
        sys.stderr.write('make_levels_list_prepub: %s\n' % problem)
        sys.exit(1)
