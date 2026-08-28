#!/usr/bin/env python3
"""Put unique level identifiers into the tables `conv_out` writes.

A level of a published analysis is known by an identifier -- the key a
database holds it under -- and `conv_out` has no way of knowing it: the
only name a level has in its output is `Lev_#`, a position in this
run's ordering that means nothing outside it.  This carries the
identifiers in from the file that has them, the same tab-delimited
experimental level list `subst_exp_levs` takes.

Usage:

    insert_level_ids <source> <levels file> [<lines file>]

    the source's columns, by number or by header name:
        --idcol <col>      the identifier    (default: level_id)
        --ecol <col>       energy            --units cm-1|kK|eV
        --jcol <col>       J                 --skip <n>   header lines
        --pcol <col>       parity: o|odd|u|1|*, e|even|g|0
        --excludecol <col> exclude a level: anything non-blank, non-zero
        --tol <value>      how far apart two energies may be and still be
                           one level, in --units (default 1)

`<levels file>` and `<lines file>` are the two files of a `conv_out`
run, or the two a `calc_TP` run writes -- `out_levels` and `out_file`,
which carry the same `Lev_#`, `lid1` and `lid2` columns for exactly
this purpose.  Both are rewritten in place, with what they held kept as
`<file>.bak`; the line list may be left out, and then only the level
list is done.

A source level is tied to a level of the list on **(parity, J,
energy)**, as `subst_exp_levs` ties one whose composition it cannot
use.  The energies are not being compared so much as recognised: the
`Ee` column holds the observed energies of the RCE deck, which came
from a list like this one in the first place.  `--tol` is there for the
case where the two files are different vintages of the analysis.  Every
pair within it is considered, closest first, and each level of either
side is used once, so two levels of one J cannot both take the target
between them.

The identifier goes into the level list as a `level_id` column beside
`Ee`, and into the line list as `id1` and `id2` -- the lower level's
and the upper level's -- beside the `lid1` and `lid2` they are looked
up through.  Nothing is matched twice for the line list: `lid1` and
`lid2` are `Lev_#`, so the levels are already known.

Columns of those names already in a file are **updated**, not
duplicated: they keep their place, are emptied, and are filled with
what this run matched.  So a level that has lost its identifier since
the last run loses it here too, rather than keeping a stale one.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import levids                                        # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


#: The options that describe the source file, and nothing else.
SOURCE_OPTIONS = levids.Spec.NAMES + ('units', 'skip', 'tol')


class Options:
    """What the command line asked for."""

    def __init__(self, argv):
        self.files = []
        self.source = {}

        rest = list(argv)
        while rest:
            arg = rest.pop(0)
            if not arg.startswith('--'):
                self.files.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.strip().lower().replace('-', '_')
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name)
                value = rest.pop(0)
            if name not in SOURCE_OPTIONS:
                raise Usage('--%s is not an option' % name)
            self.source['tolerance' if name == 'tol' else name] = value

        if not 1 < len(self.files) < 4:
            raise Usage('a source file, a level list and, if it is wanted, a '
                        'line list are what this takes, and %d file name%s '
                        'given'
                        % (len(self.files),
                           ' was' if len(self.files) == 1 else 's were'))
        if 'skip' in self.source:
            self.source['skip'] = _number(self.source['skip'], '--skip')
        try:
            self.spec = levids.Spec(**self.source)
        except levids.LevelIdError as problem:
            # What the options say is part of the command line, so it
            # is checked here and reported as a usage error.
            raise Usage(str(problem))

        self.source_path = self.files[0]
        self.levels_path = self.files[1]
        self.lines_path = self.files[2] if len(self.files) > 2 else None


def _number(text, name):
    try:
        value = int(text)
    except ValueError:
        raise Usage('%s wants a number, not %r' % (name, text))
    if value < 0:
        raise Usage('%s is not negative' % name)
    return value


def main(argv):
    if not argv or '--help' in argv or '-h' in argv or '/?' in argv:
        sys.stdout.write(__doc__.split('Usage:')[1])
        return 0

    options = Options(argv)
    spec = options.spec
    source, report = levids.read_source(options.source_path, spec)

    levels = levids.Table.read(options.levels_path, levids.LEVEL_KEYS)
    targets, by_number = levids.read_levels(levels)
    # The tolerance is given in the source's unit and the energies of
    # both files are in the run's, so it is scaled with them.
    report.lost = levids.match(source, targets, spec.tolerance * spec.scale)
    levids.insert_levels(levels, targets, report)

    lines = None
    if options.lines_path:
        lines = levids.Table.read(options.lines_path, levids.LINE_KEYS)
        levids.insert_lines(lines, by_number, report)

    levels.write(backup=os.path.basename(options.levels_path) + '.bak')
    if lines is not None:
        lines.write(backup=os.path.basename(options.lines_path) + '.bak')

    report.write(sys.stdout, spec.tolerance, spec.units)
    print('\n%s: %d of %d levels given a %s.'
          % (options.levels_path, report.matched,
             report.matched + report.blank, levids.LEVEL_ID))
    if lines is not None:
        print('%s: %s written on %d transitions.'
              % (options.lines_path, ' and '.join(levids.LINE_IDS),
                 report.lines))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('insert_level_ids: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (levids.LevelIdError, OSError) as problem:
        sys.stderr.write('insert_level_ids: %s\n' % problem)
        sys.exit(1)
