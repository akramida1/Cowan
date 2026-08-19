#!/usr/bin/env python3
"""Measure how well each coupling scheme describes each configuration.

Reads a levels.txt written with `printout --dominant-config` and writes
`coupling_choice.txt`: one line per configuration of the deck, one
column per coupling scheme the levels file holds, and in each cell the
**purity** of that configuration in that scheme -- the unweighted mean,
over the levels whose dominant configuration is that one, of the leading
eigenvector component's percentage.

The file is then **marked by hand**: a `*` anywhere in a scheme's column
says that configuration is to be written in that scheme, and
`make_levels_list_prepub.py` reads the marks back.  A star is written
already, on the highest purity, as a suggestion.  It is not the answer:
the highest purity is the main criterion but not the only one, and
previous literature on the spectrum and the schemes used for the
neighboring isoelectronic ions routinely outweigh a few percent gained
by a rarely used scheme.  That judgment is the author's and this script
does not make it.

Usage:

    config_purity [--levels levels.txt] [-o coupling_choice.txt]
                  [--drop-full-shells 5p,4d]   as printout was run
                  [--force]                    overwrite a marked file
                  [--merge]                    keep its stars, redo the sums

Why the file is fixed width and levels.txt is not: this one is edited by
a person, and the rest are read in a spreadsheet.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import in36, purity                               # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


class Options:

    def __init__(self, argv):
        self.levels = 'levels.txt'
        self.out = purity.CHOICE
        self.drop = []
        self.force = False
        self.merge = False

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
            if name in ('force', 'merge'):
                setattr(self, name, True)
                continue
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name.replace('_', '-'))
                value = rest.pop(0)
            if name == 'levels':
                self.levels = value
            elif name in ('o', 'out'):
                self.out = value
            elif name == 'drop_full_shells':
                self.drop = [_nl(s) for s in value.split(',') if s.strip()]
            else:
                raise Usage('--%s is not an option' % name.replace('_', '-'))

        if positional:
            self.levels = positional[0]
        if len(positional) > 1:
            raise Usage('too many file names: %s' % ', '.join(positional))
        if self.force and self.merge:
            raise Usage('--force throws the stars away and --merge keeps '
                        'them; ask for one of them')


def _nl(text):
    """`5p` -> (5, 1)."""
    text = text.strip()
    for i, ch in enumerate(text):
        if ch.isalpha():
            try:
                return (int(text[:i]), in36.LBCD.index(ch.lower()))
            except ValueError:
                break
    raise Usage('%r is not a subshell like 5p' % text)


def main(argv):
    options = Options(argv)
    table = purity.read_levels_txt(options.levels)
    rows = purity.deck_configurations('.', drop=options.drop)
    kept, empty, nameless = purity.purity_table(table, rows)

    print('%s: %d levels, %d couplings (%s)'
          % (options.levels, len(table.rows), len(table.couplings),
             ', '.join(table.couplings)))
    if nameless:
        raise purity.PurityError(
            '%d levels of %s have a conf0 that is not a configuration of '
            'IN36, the first being %r on line %d.  IN36 and levels.txt are '
            'not from the same run, or printout was given a different '
            '--drop-full-shells than this run was.'
            % (len(nameless), options.levels, nameless[0].conf0,
               nameless[0].line))
    if empty:
        print('   %d configurations dominate no level and are left out: %s'
              % (len(empty), ', '.join(r.label for r in empty)))

    # The suggestion, before anything a marked file says overrides it.
    for row in kept:
        row.chosen = purity.suggest(row, table.couplings)

    if os.path.exists(options.out):
        if options.merge:
            gone = purity.merge_choice(options.out, kept, table.couplings)
            print('   merged the marks of %s' % options.out)
            if gone:
                print('   these were marked and are no longer in the deck: %s'
                      % ', '.join(gone))
        elif not options.force:
            raise Usage(
                '%s exists.  It is the one file of this pipeline that is '
                'edited by hand, so it is not overwritten by accident: '
                '--merge keeps its stars and writes the new purities under '
                'them, --force throws them away.' % options.out)

    purity.write_choice(options.out, kept, table.couplings,
                        folder='.', levels=options.levels)
    print('wrote %s (%d configurations)' % (options.out, len(kept)))
    _report(kept, table.couplings)
    return 0


def _report(rows, couplings):
    """What the table says, in one line: where the schemes disagree."""
    count = {}
    for row in rows:
        if row.chosen:
            count[row.chosen] = count.get(row.chosen, 0) + 1
    print('   suggested: %s'
          % ', '.join('%s for %d' % (c, count[c])
                      for c in couplings if c in count))
    close = [r for r in rows
             if r.chosen and len([c for c in couplings
                                  if abs(r.purity.get(c, -99)
                                         - r.purity[r.chosen]) <= 2]) > 1]
    if close:
        print('   %d configurations are within 2%% between two schemes and '
              'are the ones worth a decision: %s'
              % (len(close), ', '.join(r.label for r in close)))


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('config_purity: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except purity.PurityError as problem:
        sys.stderr.write('config_purity: %s\n' % problem)
        sys.exit(1)
