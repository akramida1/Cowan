#!/usr/bin/env python3
"""Write the table that translates level identifiers into IDEN2's numbers.

A level of a published analysis is known by an identifier -- the key a
database holds it under -- and IDEN2 knows nothing of it.  IDEN2 numbers
its levels itself, from the highest energy down, and that number is the
only name a level has in `ENLEV.DAT` and `TRANS.DAT`.  It means nothing
outside one pair of those files: rebuild them with `conv_out` and every
number moves.

So anything that has to say the same thing to IDEN2 and to the rest of
the analysis -- a level optimization's output, above all -- needs the two
names side by side.  This writes them, as `IDEN_level_ids.txt`, and
`update_IDEN` is what reads it.

Usage:

    write_IDEN_level_ids <source>

    the source's columns, by number or by header name:
        --idcol <col>      the identifier    (default: level_id)
        --ecol <col>       energy            --units cm-1|kK|eV
        --jcol <col>       J                 --skip <n>   header lines
        --pcol <col>       parity: o|odd|u|1|*, e|even|g|0
        --excludecol <col> exclude a level: anything non-blank, non-zero
        --tol <value>      how far apart two energies may be and still be
                           one level, in --units (default 1)

`<source>` is the tab-delimited experimental level list `subst_exp_levs`
and `insert_level_ids` take, and it is the only argument: `ENLEV.DAT` is
read from the working directory under that name, because that is the
name IDEN2 gives it, and `IDEN_level_ids.txt` is written beside it.
Neither file is modified.

Only the **fixed** levels of `ENLEV.DAT` are matched -- the ones carrying
the star, which are the experimental ones -- and they are matched on
**(parity, J, energy)**, as `insert_level_ids` matches, against the
*fixed* energy `ENLEV.DAT` holds for them.  Every pair within `--tol` is
considered, closest first, and each level of either side is used once, so
two levels of one J cannot both take the level between them.

Parity comes from the sixth character of the IDEN2 label, `_` for even
and `~` for odd; `conv_out` puts it there from the parity of the
configuration set the level belongs to.  A source level `ENLEV.DAT` has
no counterpart for is reported and left out -- `conv_out` writes only as
many levels as IDEN2 can hold, so a large calculation reaches it
truncated, and that is the usual reason.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, iden, levids                           # noqa: E402


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

        if len(self.files) != 1:
            raise Usage('a source level list is what this takes, and %d file '
                        'name%s given'
                        % (len(self.files),
                           ' was' if len(self.files) == 1 else 's were'))
        if 'skip' in self.source:
            self.source['skip'] = _number(self.source['skip'], '--skip')
        try:
            self.spec = levids.Spec(**self.source)
        except levids.LevelIdError as problem:
            # What the options say is part of the command line, so it is
            # checked here and reported as a usage error.
            raise Usage(str(problem))

        self.source_path = self.files[0]


def _number(text, name):
    try:
        value = int(text)
    except ValueError:
        raise Usage('%s wants a number, not %r' % (name, text))
    if value < 0:
        raise Usage('%s is not negative' % name)
    return value


def _matched(source, levels, tolerance, flip=False):
    """Match once, and say how many levels the match tied.

    The source levels carry the result of a match, so they are cleared
    first: this is called more than once on the same list.
    """
    for one in source:
        one.matched = None
    targets = iden.fixed_targets(levels, flip=flip)
    lost = levids.match(source, targets, tolerance)
    return targets, lost


def main(argv):
    if not argv or '--help' in argv or '-h' in argv or '/?' in argv:
        sys.stdout.write(__doc__.split('Usage:')[1])
        return 0

    options = Options(argv)
    spec = options.spec
    source, report = levids.read_source(options.source_path, spec)

    path = files.resolve(iden.ENLEV, '.')
    levels = iden.read_enlev(path)
    # The tolerance is given in the source's unit and the energies of
    # both files are in the run's, so it is scaled with them.
    tolerance = spec.tolerance * spec.scale
    targets, report.lost = _matched(source, levels, tolerance)
    report.matched = sum(1 for one in targets if one.ident)
    report.blank = len(targets) - report.matched

    # The parity marks are the one thing here that is a convention rather
    # than a number, so the reading of them is measured instead of
    # trusted: if the file were marked the other way round, the wrong
    # reading would tie far fewer levels, and that is worth saying out
    # loud rather than leaving as a quiet half-empty table.
    _, other = _matched(source, levels, tolerance, flip=True)
    flipped = len(source) - len(other)
    targets, report.lost = _matched(source, levels, tolerance)

    pairs = sorted(((one.ident, one.number) for one in targets if one.ident),
                   key=lambda pair: pair[1])
    levids.write_iden_ids(pairs)

    report.write(sys.stdout, spec.tolerance, spec.units)
    if flipped > report.matched:
        sys.stdout.write(
            '   Warning: reading the parity marks the other way round -- `_` '
            'odd, `~` even -- would tie %d levels rather than %d.  Is %s the '
            'ENLEV.DAT of this calculation?\n'
            % (flipped, report.matched, path))
    print('\n%s: %d of %d fixed levels in %s given a %s.'
          % (levids.IDEN_IDS, report.matched, len(targets), path,
             levids.SOURCE_ID))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('write_IDEN_level_ids: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (levids.LevelIdError, iden.IdenError, OSError) as problem:
        sys.stderr.write('write_IDEN_level_ids: %s\n' % problem)
        sys.exit(1)
