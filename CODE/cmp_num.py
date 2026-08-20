#!/usr/bin/env python3
"""Compare two Fortran output files, or two run folders, numerically.

`diff` is the acceptance test of CLAUDE.md 1 and stays so for every
change that is meant to be bit-exact.  This is the tool for the changes
that are not: it requires the *text* around the numbers to match exactly,
and then reports how far the numbers moved -- max, rms, a decade
histogram and the worst offenders -- instead of answering only
"different".

Usage:

    cmp_num [options] OLD NEW
    cmp_num [options] OLDFOLDER NEWFOLDER

With two folders, every file present in both is compared, by the name it
has in OLD.  A Fortran unformatted file cannot be tokenized, so it is
handled by record layout where one is known and by md5 where none is:
TAPE2N is read with `toolbox/tape2n.py` and compared field by field like
a text file, and anything else unformatted gets an md5 answer instead --
identical, or differing by an amount this tool cannot quantify.  Neither
kind is ever silently skipped.

Options:

    --rtol R      relative tolerance, default 0 (exact)
    --atol A      absolute tolerance, default 0.  A pair passes if it is
                  within *either* tolerance, so --atol covers the fields
                  printed to a fixed number of decimals, where a value
                  near zero has no meaningful relative deviation.
    --top N       how many worst offenders to list, default 10
    --ignore RE   skip lines matching RE; repeatable
    --no-default-ignores
                  compare the timing lines too, which are wall clock and
                  never agree.  Only useful for checking this tool.
    --d-exponents read D as an exponent letter as well as E.  Off by
                  default because no output file here has one, while
                  configuration labels like `4d7` are everywhere.
    --quiet       print only the summary lines, no per-file detail
    --only NAME   in folder mode, compare just this one file; repeatable

Exit status:

    0   text matched and every number was within tolerance
    1   numbers outside tolerance, or an unformatted file differed;
        the text is intact either way
    2   structural mismatch, a file missing, or a bad command line

Two failing codes for two different kinds of answer: 1 says the change
did what it was supposed to do and here is how big it was, 2 says
something moved that was never allowed to.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import cmpnum                                     # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


def parse(argv):
    """Hand-rolled, to match the other scripts here and to keep the one
    dependency of this package a Python install and nothing else."""
    options = {'rtol': 0.0, 'atol': 0.0, 'top': 10, 'quiet': False,
               'ignores': list(cmpnum.DEFAULT_IGNORES), 'only': [],
               'd_exponents': False}
    rest = []
    i = 0
    while i < len(argv):
        word = argv[i]
        if word in ('--rtol', '--atol', '--top', '--ignore', '--only'):
            if i + 1 >= len(argv):
                raise Usage('%s needs a value' % word)
            value = argv[i + 1]
            i += 2
            try:
                if word == '--rtol':
                    options['rtol'] = float(value)
                elif word == '--atol':
                    options['atol'] = float(value)
                elif word == '--top':
                    options['top'] = int(value)
            except ValueError:
                raise Usage('%s wants a number, not %r' % (word, value))
            if word == '--ignore':
                options['ignores'].append(value)
            elif word == '--only':
                options['only'].append(value)
            continue
        if word == '--no-default-ignores':
            for pattern in cmpnum.DEFAULT_IGNORES:
                if pattern in options['ignores']:
                    options['ignores'].remove(pattern)
        elif word == '--d-exponents':
            options['d_exponents'] = True
        elif word == '--quiet':
            options['quiet'] = True
        elif word.startswith('-'):
            raise Usage('unknown option %r' % word)
        else:
            rest.append(word)
        i += 1
    if len(rest) != 2:
        raise Usage('two file or folder names are needed, got %d' % len(rest))
    return options, rest[0], rest[1]


def pairs(old, new, only):
    """The (old, new) file pairs to compare, and the names skipped.

    Folder mode walks OLD, not NEW, so a file the new run failed to write
    is reported as missing rather than quietly dropped from the
    comparison.
    """
    if os.path.isfile(old) and os.path.isfile(new):
        return [(old, new)], []
    if not (os.path.isdir(old) and os.path.isdir(new)):
        raise Usage('%s and %s must be two files or two folders' % (old, new))
    found, missing = [], []
    for name in sorted(os.listdir(old)):
        if only and name not in only:
            continue
        a = os.path.normpath(os.path.join(old, name))
        b = os.path.normpath(os.path.join(new, name))
        if not os.path.isfile(a):
            continue
        if not os.path.isfile(b):
            missing.append(name)
            continue
        found.append((a, b))
    return found, missing


def main(argv):
    options, old, new = parse(argv)
    todo, missing = pairs(old, new, options['only'])
    if not todo:
        raise Usage('nothing to compare in %s' % old)
    status = 0
    summary = []
    for a, b in todo:
        name = os.path.basename(a)
        report = cmpnum.compare_files(
            a, b, rtol=options['rtol'], atol=options['atol'],
            ignores=options['ignores'], d_exponents=options['d_exponents'])
        if not options['quiet']:
            print(cmpnum.format_report(report, top=options['top'],
                                       rtol=options['rtol'],
                                       atol=options['atol']))
            print('')
        if report.binary and not report.layout:
            if report.binary_same:
                summary.append('%-16s unformatted, md5 identical' % name)
            else:
                status = max(status, 1)
                summary.append('%-16s unformatted, md5 DIFFERS -- not '
                               'quantifiable here' % name)
        elif report.structural:
            status = 2
            summary.append('%-16s STRUCTURAL: %d line(s)'
                           % (name, len(report.structural)))
        elif report.offenders:
            status = max(status, 1)
            summary.append('%-16s %d of %d outside tolerance, max rel %.3e'
                           % (name, len(report.offenders), report.compared,
                              report.largest))
        else:
            summary.append('%-16s %d numbers, max rel %.3e -- within tolerance'
                           % (name, report.compared, report.largest))
    for name in missing:
        status = 2
        summary.append('%-16s MISSING from %s' % (name, new))
    print('cmp_num %s -> %s' % (old, new))
    for line in summary:
        print('  ' + line)
    return status


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('cmp_num: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except OSError as problem:
        sys.stderr.write('cmp_num: %s\n' % problem)
        sys.exit(2)
