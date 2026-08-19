#!/usr/bin/env python3
"""Carry observed energies, and optionally parameters, into an RCE deck.

The Python replacement for `subst_exp_levs.pl`.  A calculation is not
fitted once: a configuration is added, a parameter is freed, the scaling
is changed, and the whole deck is rebuilt -- and every observed level
has to be tied to its eigenvector again in the new deck.  The tie cannot
be made on energy, because the energies are what the new fit is going to
change; it is made on the **eigenvector**, which is what a level is.

Usage:

    subst_exp_levs <source> <target> [params=<what>]

`<source>` is the deck the observed energies are in and `<target>` is
the deck they are being put into; the target is rewritten in place, with
what it held kept as `<target>.bak`.  `<source>` may instead be a LOPT
level list -- a file of ten tab-separated columns, or one named `.LEV` --
in which case the parameter section is not touched and the argument in
third place is read as an uncertainty above which a level is left out.

`params=<what>` says how much of the source's parameter table to bring
across as well:

    (nothing)  every value and every least-squares flag
    flags      the flags only, leaving the values as they are
    noflags    the values only, leaving the flags as they are
    Eav        the Eav values only;  F, G, Z, CI likewise, by name
    none       nothing: the parameter section is copied unchanged

Three reports go to standard output when the run is done: the matches
worth checking by hand, any two levels of one J that came out with the
same observed energy, and the configurations experiment has reached.

**One difference from the Perl.**  `params=noflags` substitutes nothing
at all in `subst_exp_levs.pl`: the value is held by the test at `:448`,
which asks whether the parameter's *name* contains the word `noflags`
and is therefore always true, and the flag is held at `:451`.  Here it
does what its name says and what its pairing with `flags` means -- the
values, and not the flags.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import explevs, files, rce                         # noqa: E402

OPTIONS = ('flags', 'noflags', 'eav', 'f', 'g', 'z', 'ci', 'none')


class Usage(Exception):
    """The command line will not do what was asked."""


# --------------------------------------------------------------------------
# the parameter table
# --------------------------------------------------------------------------

class Parameters:
    """Which half of a parameter line to take from the source deck.

    `Eav` and `CI` name a **kind** of parameter and the rest name a
    letter its name begins with, which is why the Perl tests them in
    three separate branches (`subst_exp_levs.pl:402`, `:418`, `:443`)
    rather than one.  Matching `CI` against the name of a CI parameter
    would never succeed -- they are called `120D1112` -- and matching
    `F` against a name would catch `EAV 4f3` along with `F2(11)`.
    """

    def __init__(self, option):
        self.option = (option or '').lower()

    @property
    def skip(self):
        return self.option == 'none'

    def wants(self, kind, name):
        """(take the value, take the flag) for one parameter."""
        if self.option == '':
            return True, True
        if self.option == 'flags':
            return False, True
        if self.option == 'noflags':
            return True, False
        if kind in ('eav', 'ci'):
            return self.option == kind, False
        return self.option in name.lower(), False


def walk_parameters(lines):
    """Yield `(line number, key, kind, name)` for every parameter.

    A parameter's name is not unique in a file.  `F2(11)` belongs to
    every configuration that has an f shell, and `120D1112` to every
    pair of configurations it couples, so the name alone would make the
    parameters of one configuration overwrite another's.  Two things
    make the key unique, and they are the two `subst_exp_levs.pl:126-146`
    uses:

    * a single-configuration parameter is keyed by the configuration it
      belongs to, which is the one the preceding `EAV` line named;
    * a CI parameter is keyed by the pair of configurations printed
      beside it and by its **position** among that pair's parameters,
      which needs no assumption that a pair never repeats a name.
    """
    section, inside, configuration = 0, False, ''
    pair, index = '', 0
    for i, raw in enumerate(lines):
        # Padded for the column tests only; what is yielded is the line
        # itself, trailing blanks and all.
        line = raw.ljust(80)
        if rce.PARAMETER_HEADING in line:
            section += 1
            inside, configuration, pair, index = True, '', '', 0
            continue
        if inside and rce.is_control_card(line):
            inside = False
        if not inside or not rce.is_parameter_line(line):
            continue
        name = line[:10].strip()
        if rce.is_ci_parameter(line):
            here = line[57:71].strip()
            if here != pair:
                pair, index = here, 0
            index += 1
            yield i, (section, pair, index), 'ci', name
        elif name.upper().startswith('EAV'):
            configuration = name[3:].strip()
            yield i, (section, configuration, name), 'eav', name
        else:
            yield i, (section, configuration, name), 'plain', name


def substitute_parameters(lines, source_lines, wanted):
    """Put the source deck's parameter values or flags into `lines`."""
    if wanted.skip:
        return lines, 0
    source = {}
    for i, key, _, _ in walk_parameters(source_lines):
        source[key] = source_lines[i]

    out, done = list(lines), 0
    for i, key, kind, name in walk_parameters(lines):
        theirs = source.get(key)
        if theirs is None:
            continue
        take_value, take_flag = wanted.wants(kind, name)
        if not (take_value or take_flag):
            continue
        line = out[i]
        head = theirs[:15] if take_flag else line[:15]
        body = theirs[15:29] if take_value else line[15:29]
        out[i] = head + body + line[29:]
        done += 1
    return out, done


# --------------------------------------------------------------------------
# the reports
# --------------------------------------------------------------------------

def report(out, source, target, lost):
    for level in lost:
        out.write('Warning: no level of the target deck matches the source '
                  'level at parity set %d, J = %s, Ecalc = %g, Eobs = %g.\n'
                  % (level.kset, level.j, level.ec, level.ee))

    doubtful = explevs.suspicious(source)
    if doubtful:
        out.write('\nMatches worth checking:\n')
        for level, best in doubtful:
            out.write('J=%4s source: %13.6f%10.3f%s\n'
                      % (level.j, level.ee, level.ec,
                         explevs.vector_text(level)))
            out.write('J=%4s target: %13.6f%10.3f%s\n\n'
                      % (best.j, best.assigned, best.ec,
                         explevs.vector_text(best)))

    same = explevs.duplicates(target)
    if same:
        out.write('\nLevels of one J given the same observed energy:\n')
        for before, after in same:
            for one in (before, after):
                out.write('J=%4s%13.6f%10.3f%s\n'
                          % (one.j, one.assigned, one.ec,
                             explevs.vector_text(one)))
            out.write('\n')

    for kset, found in sorted(explevs.known_configurations(target).items()):
        out.write('\nConfigurations reached by observation, parity set %d:\n'
                  % kset)
        for configuration, _ in found:
            out.write('%s\n' % configuration)


# --------------------------------------------------------------------------

def _is_lopt(path):
    """Whether the source is a LOPT level list rather than an RCE deck."""
    if str(path).upper().endswith('.LEV'):
        return True
    for line in files.read_lines(path, strip=False)[:400]:
        if explevs._LOPT.match(line):
            return True
    return False


def main(argv):
    if not 2 <= len(argv) <= 3:
        raise Usage('two or three arguments are needed, a source deck, a '
                    'target deck and optionally params=<what>')
    source_path, target_path = argv[0], argv[1]
    third = argv[2] if len(argv) > 2 else ''

    lopt = _is_lopt(source_path)
    option, most = '', None
    if third:
        if third.lower().startswith('params'):
            option = third.partition('=')[2].strip().lower()
            if option not in OPTIONS:
                raise Usage('params=%s is none of %s'
                            % (option, ', '.join(OPTIONS)))
        elif lopt:
            try:
                most = float(third)
            except ValueError:
                raise Usage('%r is neither params=<what> nor an uncertainty'
                            % third)
        else:
            raise Usage('%r is not params=<what>' % third)

    source = explevs.read_lopt(source_path, most) if lopt \
        else explevs.read_deck(source_path)
    target = explevs.read_deck(target_path)
    lost = explevs.identify(source, target)

    lines = files.read_lines(target_path, strip=False)
    given = 0
    for level in target.levels():
        if level.assigned is None or level.line < 0:
            continue
        lines[level.line] = rce.set_observed(lines[level.line], level.assigned)
        given += 1

    parameters = 0
    if not lopt:
        lines, parameters = substitute_parameters(
            lines, files.read_lines(source_path, width=140),
            Parameters(option))

    files.write_deck(target_path, lines, files.line_ending(target_path),
                     backup=os.path.basename(target_path) + '.bak')
    report(sys.stdout, source, target, lost)
    print('\n%s: %d observed energies and %d parameters substituted; %d '
          'source levels went unmatched.'
          % (target_path, given, parameters, len(lost)))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('subst_exp_levs: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (explevs.ExpLevelsError, rce.RceError, OSError) as problem:
        sys.stderr.write('subst_exp_levs: %s\n' % problem)
        sys.exit(1)
