#!/usr/bin/env python3
"""Pull the classification columns out of a `calc_TP` line list.

A line classification code does not want the fifty columns `calc_TP`
writes.  It wants, for every transition between two levels of the
published analysis, the two level identifiers, the wavenumber, gA and
its uncertainty, and the upper level's energy -- and, where a plasma
temperature is assumed, a synthetic intensity built out of those.  This
writes that table.

Usage:

    extract_Icalc [<tp file>] [<output file>]

        <tp file>       the calc_TP line list, tp_E1.txt by default,
                        with id1 and id2 already put on it by
                        insert_level_ids
        <output file>   Icalc.txt by default

        --T_eV <value>    the temperature, in eV
        --d_off <value>   the log-offset of the intensity

The columns written are, in order:

    id1     the lower level's identifier
    id2     the upper level's identifier
    rwn     E_up_exp - E_low_exp, in cm-1
    gA      (2*J2+1) * A(s-1)
    u%gA    std_A%, the relative standard deviation of A and so of gA
    Eup     E_up_exp, in cm-1

`calc_TP` carries its energies in kK and this writes cm-1, so `rwn` and
`Eup` are its values times 1000.

Give BOTH `--T_eV` and `--d_off` and a seventh column follows:

    Icalc   exp(-Eup/T_cm + d_off) * gA / rwn,  T_cm = T_eV * (eV in cm-1)

a Boltzmann population of the upper level times its rate, per unit
wavenumber.  `d_off` is the free normalisation -- partition function,
column density and the choice of unit all at once -- so it is given
rather than computed.  Neither option alone will do: the column is then
not written at all, rather than written empty.

Only transitions with an identifier at BOTH ends are taken, and of
those only the ones with an experimental energy at both ends; the run
says how many of each it passed over.  See `CODE/toolbox/icalc.py`.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files                                         # noqa: E402
from toolbox import icalc                                         # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


#: The options, and the attribute each is kept under.
OPTIONS = {'t_ev': 'temperature', 'd_off': 'offset'}


class Options:
    """What the command line asked for."""

    def __init__(self, argv):
        files_named = []
        self.temperature = None
        self.offset = None

        rest = list(argv)
        while rest:
            arg = rest.pop(0)
            if not arg.startswith('--'):
                files_named.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.strip().lower().replace('-', '_')
            if name not in OPTIONS:
                raise Usage('--%s is not an option' % name)
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name)
                value = rest.pop(0)
            try:
                setattr(self, OPTIONS[name], float(value))
            except ValueError:
                raise Usage('--%s wants a number, not %r' % (name, value))

        if len(files_named) > 2:
            raise Usage('a line list and an output name are what this takes, '
                        'and %d file names were given' % len(files_named))
        self.source = files_named[0] if files_named else icalc.SOURCE
        self.output = (files_named[1] if len(files_named) > 1
                       else icalc.OUTPUT)

        if self.temperature is not None and self.temperature <= 0.0:
            raise Usage('--T_eV is a temperature; %g is not one'
                        % self.temperature)
        # Both or neither.  One of the two on its own is a command line
        # the user got wrong, and quietly writing the six-column table
        # would hide that until the classification code went looking for
        # the seventh.
        if (self.temperature is None) != (self.offset is None):
            raise Usage('--T_eV and --d_off go together; the Icalc column '
                        'needs both of them')

    @property
    def intensity(self):
        return self.temperature is not None


def main(argv):
    if '--help' in argv or '-h' in argv or '/?' in argv:
        sys.stdout.write(__doc__.split('Usage:')[1])
        return 0

    options = Options(argv)
    table = icalc.read(options.source)
    temperature = (icalc.to_wavenumbers(options.temperature)
                   if options.intensity else None)
    rows, report = icalc.extract(table, temperature, options.offset)

    # The output keeps the line ending the source has, as every rewrite
    # in this package does: a table read on the other platform should
    # not differ from its source on every line.
    icalc.write(options.output, rows, options.intensity,
                files.line_ending(options.source))

    report.write(sys.stdout)
    if options.intensity:
        print('Icalc at T = %g eV (%.1f cm-1), offset %g.'
              % (options.temperature, temperature, options.offset))
    print('%s written.' % options.output)
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('extract_Icalc: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (icalc.IcalcError, OSError) as problem:
        sys.stderr.write('extract_Icalc: %s\n' % problem)
        sys.exit(1)
