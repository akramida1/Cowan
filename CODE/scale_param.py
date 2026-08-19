#!/usr/bin/env python3
"""Scale groups of Slater parameters in a copy of RCEINP.

The Python replacement for `scale_param.pl`.  A least-squares fit is
started from Hartree-Fock parameters scaled down by a factor -- 0.85 for
the electrostatic ones is Cowan's own recommendation -- and this is how
that starting deck is made, and how a fit that will not converge is
restarted from a different scaling with some parameters held fixed.

Usage:

    scale_param <in file> <out file> <factor> [<factor> ...] <sets> [FIX|ALL]

A factor is one of

    F=<x>        every Fk parameter                  (direct integrals)
    G=<x>        every Gk parameter                  (exchange integrals)
    CI=<x>       every configuration-interaction parameter
    GROUP<n>=<x> every parameter whose fitting group is n

and any number of them may be given, in any order.  `<sets>` is 1 for
the first parity, 2 for the second, or 12 for both.

`FIX` additionally holds the parameters this run scaled, by setting
their least-squares flag to RCE's fixed code and keeping its sign.
`ALL` holds **every** parameter of the named sets, scaled or not.  The
Perl's usage message calls the second one FIX (`scale_param.pl:41-42`)
while its code tests for ALL (`:30`); both spellings are accepted here
and they do the two different things the code does.

A parameter matches at most one factor, in the order F, G, CI, group,
which is the Perl's order and matters where a group holds both kinds.

**This program does something the Perl no longer does.**
`scale_param.pl:58` recognises a parameter line with a regex that
requires the maximum-value field to read `0.0000`, and RCE has written
that field as F14.6 -- `0.000000` -- for as long as any deck here goes
back.  The regex therefore matches nothing, and the Perl copies its
input unchanged: verified on `WORK/W70`, `WORK/Pr3` and `WORK/U90`,
where it changes not one line for any combination of factors.  The
fields are taken here from the FORMAT that writes them
(`FOR/RCE20KD.F:1470`), so there is nothing to fall out of step with.
"""

import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, rce                                  # noqa: E402

_FACTOR = re.compile(r'^(F|G|CI|GROUP\s*(-?\d+))\s*=\s*(\S+)$', re.I)


class Usage(Exception):
    """The command line will not do what was asked."""


class Factors:
    """What to multiply each kind of parameter by."""

    def __init__(self):
        self.f = self.g = self.ci = 1.0
        self.groups = {}

    def add(self, argument):
        """Read one `F=`, `G=`, `CI=` or `GROUP<n>=` option."""
        found = _FACTOR.match(argument)
        if not found:
            return False
        kind, group, text = found.group(1).upper(), found.group(2), \
            found.group(3)
        try:
            value = float(text)
        except ValueError:
            raise Usage('%r is not a number' % text)
        if group is not None:
            self.groups[abs(int(group))] = value
        elif kind == 'F':
            self.f = value
        elif kind == 'G':
            self.g = value
        else:
            self.ci = value
        return True

    def of(self, line):
        """The factor that applies to one parameter line, or None.

        The name of a direct or exchange integral is `Fk(ij)` or
        `Gk(ij)`, so its first letter is the test -- and a CI
        parameter's name begins with a digit, which is why the letters
        can be tested first without catching one.
        """
        name = line[:10].strip().upper()
        if name.startswith('F') and self.f != 1.0:
            return self.f
        if name.startswith('G') and self.g != 1.0:
            return self.g
        if rce.is_ci_parameter(line) and self.ci != 1.0:
            return self.ci
        factor = self.groups.get(rce.parameter_group(line))
        return factor if factor not in (None, 1.0) else None

    @property
    def empty(self):
        return (self.f == 1.0 and self.g == 1.0 and self.ci == 1.0
                and not any(v != 1.0 for v in self.groups.values()))


def scale(lines, factors, sets, fix, fix_all):
    """Rewrite the parameter tables of the named parity sets.

    A table opens on RCE's `PARAMETER FLAG` heading and runs to its next
    run control card; everything else in the file, including the levels
    and the headings themselves, is passed through as it stands.
    """
    out, section, scaled, fixed = [], 0, 0, 0
    inside = False
    for line in lines:
        if rce.PARAMETER_HEADING in line:
            section += 1
            inside = True
            out.append(line)
            continue
        if inside and rce.is_control_card(line):
            inside = False
        if not inside or section not in sets \
                or not rce.is_parameter_line(line):
            out.append(line)
            continue
        factor = factors.of(line)
        value = None
        if factor is not None:
            value = rce.parameter_value(line) * factor
            scaled += 1
        hold = fix_all or (fix and factor is not None)
        if hold:
            fixed += 1
        out.append(rce.edit_parameter_line(line, value=value, fixed=hold))
    return out, scaled, fixed


def main(argv):
    if len(argv) < 4:
        raise Usage('an input file, an output file, at least one factor and '
                    'a parity set are needed')
    source, target, rest = argv[0], argv[1], list(argv[2:])

    factors = Factors()
    while rest and factors.add(rest[0]):
        rest.pop(0)
    if factors.empty:
        raise Usage('no factor was given, or every factor given is 1')
    if not rest:
        raise Usage('the parity set is missing')

    text = rest.pop(0)
    if text not in ('1', '2', '12'):
        raise Usage('%r is not a parity set; use 1, 2 or 12' % text)
    sets = set(int(c) for c in text)

    fix = fix_all = False
    if rest:
        word = rest.pop(0).upper()
        if word == 'ALL':
            fix = fix_all = True
        elif word == 'FIX':
            fix = True
        else:
            raise Usage('%r is neither FIX nor ALL' % word)
    if rest:
        raise Usage('%r was not expected' % rest[0])

    lines = files.read_lines(source, strip=False)
    out, scaled, fixed = scale(lines, factors, sets, fix, fix_all)
    files.write_deck(target, out, files.line_ending(source))
    print('%s written: %d parameters scaled, %d held fixed.'
          % (target, scaled, fixed))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('scale_param: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (rce.RceError, OSError) as problem:
        sys.stderr.write('scale_param: %s\n' % problem)
        sys.exit(1)
