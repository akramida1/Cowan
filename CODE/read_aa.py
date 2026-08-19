#!/usr/bin/env python3
"""Tabulate the autoionization rates RCG prints, one column per channel.

The Python replacement for `read_aa.pl`.  RCG prints an Auger rate for
every (level, target state) pair of a run that has a continuum in it;
this collects them into the table `conv_out` reads, summing the target
states that belong to one core and closing a channel the level has not
enough energy to reach:

    J   Energy   1s2:0   1s2s:120.1   ...   AA_tot

Usage:

    read_aa <cores> <IP> [--set N] [--folder F]

`<cores>` is a comma-separated list of `name:energy`, the core states of
the ionized atom with their energies **above that ion's ground state**,
in the run's own energy unit (1000 cm-1 in every deck here) -- and no
blanks anywhere in it.  `<IP>` is the ionization energy of the atom
being calculated, in the same unit, which is what puts the core
energies on this ion's scale.  A core name is a **regular expression**
matched against RCG's channel labels, so one column may sum several of
them:

    read_aa 1s2:0,1s2s:120.1,1s2p:125.7 400.0
    read_aa "p3.+4S:0,p3.+2D:41.25,p3.+2P:62.44" 500.0

The table goes to standard output.  `--folder` says where OUTG11 is;
`--set` says which parity set to tabulate, the first by default.  The
table has no parity column -- a row is a J and an energy -- so one set
is one file, which is what `toolbox/ai.py` expects to read back.  The
Perl had no such option: it stopped at the first J that did not
increase (`read_aa.pl:104`) and so could only ever reach the first set.

A rate is dropped where the level lies below the core it would leave
behind, because there is then no energy for the free electron.  RCG
prints that electron's kinetic energy on its own `KE FOR CORE` lines,
which would say the same thing more directly; they are not used, here
or in the Perl, where the test is left over as a comment
(`read_aa.pl:263`).  The two do not agree in every case and which is
right is a question for the author, so what the Perl does is what is
done.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, outg11                               # noqa: E402

# A level RCG has parked out of the way rather than calculated: the
# continuum placeholders come out at -9500 and the Perl drops anything
# below -5000 (`read_aa.pl:247`).
FLOOR = -5000.0


class Usage(Exception):
    """The command line will not do what was asked."""


class CoresError(Exception):
    """A core list that does not describe this run."""


def read_cores(text, ip):
    """[(name, energy)] from `name:energy,...`, put on this ion's scale.

    The energies are given above the ground state of the **ionized**
    atom, which is where a spectroscopist has them; adding the
    ionization energy puts them above the ground state of the atom being
    calculated, which is the scale RCG's eigenvalues are on.
    """
    out = []
    for item in text.split(','):
        name, sep, energy = item.rpartition(':')
        if not sep or not name:
            raise CoresError('%r is not a core: write it as name:energy'
                             % item)
        try:
            out.append((name, float(energy) + ip))
        except ValueError:
            raise CoresError('%r is not an energy' % energy)
    out.sort(key=lambda core: core[1])
    return out


def channels_of(cores, blocks):
    """{core name: [channel, ...]}, and a complaint if one matches none.

    A core name is a pattern, so one column may gather several of RCG's
    channels -- the terms of one core -- which is what the Perl's
    `$target =~ /$core/` does (`read_aa.pl:257`).
    """
    every = sorted(set(c for block in blocks for c in block.rates))
    out, missing = {}, []
    for name, _ in cores:
        try:
            pattern = re.compile(name)
        except re.error as problem:
            raise CoresError('%r is not a usable pattern: %s'
                             % (name, problem))
        found = [c for c in every if pattern.search(c)]
        out[name] = found
        if not found:
            missing.append(name)
    if missing:
        raise CoresError(
            'no autoionization channel of this run matches %s.\n'
            'The channels it has are:\n  %s'
            % (', '.join(repr(m) for m in missing), '\n  '.join(every)))
    return out


def rows(blocks, cores, channels):
    """One row per level: (J, energy, [rate per core], total)."""
    out = []
    for block in sorted(blocks, key=lambda b: b.j):
        for i, energy in enumerate(block.eigenvalues):
            if energy < FLOOR:
                continue
            rates = []
            for name, limit in cores:
                rate = 0.0
                if energy >= limit:
                    for channel in channels[name]:
                        row = block.rates.get(channel, ())
                        rate += row[i] if i < len(row) else 0.0
                rates.append(rate)
            out.append((block.j, energy, rates, sum(rates)))
    return out


def write(out, cores, table):
    """The table as `read_aa.pl` writes it: tab separated, %9.2e rates."""
    heading = ['J', 'Energy']
    heading += ['%s:%g' % (name, energy) for name, energy in cores]
    heading.append('AA_tot')
    out.write('\t'.join(heading) + '\n')
    for j, energy, rates, total in table:
        fields = [_j(j), '%g' % energy]
        fields += ['%9.2e' % rate for rate in rates]
        fields.append('%9.2e' % total)
        out.write('\t'.join(fields) + '\n')


def _j(j):
    """J as RCG's own heading spells it: 5.0, 4.5."""
    return '%.1f' % j


def main(argv):
    arguments, kset, folder = [], 1, '.'
    rest = list(argv)
    while rest:
        one = rest.pop(0)
        if one == '--set':
            if not rest:
                raise Usage('--set needs a parity set number')
            kset = int(rest.pop(0))
        elif one == '--folder':
            if not rest:
                raise Usage('--folder needs a directory')
            folder = rest.pop(0)
        else:
            arguments.append(one)
    if len(arguments) != 2:
        raise Usage('two arguments are needed, a list of cores and an '
                    'ionization energy')
    try:
        ip = float(arguments[1])
    except ValueError:
        raise Usage('%r is not an ionization energy' % arguments[1])

    cores = read_cores(arguments[0], ip)
    every = outg11.read_auger(files.OUTG11, folder=folder)
    blocks = [b for b in every if b.kset == kset]
    if not blocks:
        raise CoresError(
            'this run has no autoionization rates in parity set %d; it has '
            'them in %s' % (kset, ', '.join(str(k) for k in
                                            sorted(set(b.kset for b in every)))
                            or 'none'))
    write(sys.stdout, cores, rows(blocks, cores, channels_of(cores, blocks)))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('read_aa: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (CoresError, outg11.Outg11Error, OSError) as problem:
        sys.stderr.write('read_aa: %s\n' % problem)
        sys.exit(1)
