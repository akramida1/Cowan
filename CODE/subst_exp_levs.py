#!/usr/bin/env python3
"""Carry observed energies, and optionally parameters, into an RCE deck.

The Python replacement for `subst_exp_levs.pl`.  A calculation is not
fitted once: a configuration is added, a parameter is freed, the scaling
is changed, and the whole deck is rebuilt -- and every observed level
has to be tied to its eigenvector again in the new deck.  The tie cannot
be made on energy, because the energies are what the new fit is going to
change; it is made on the **eigenvector**, which is what a level is.

The source may be another RCE deck, in which case nothing below is
needed and the two sides speak the same eigenvectors.  Or it may be an
experimental level list -- a level-optimization run, a spreadsheet a
spectrum's bookkeeping is kept in -- and then it has to be described,
because no two such files are laid out alike.  It is recognised as one
by not parsing as an RCE deck, it must be **tab-delimited**, and its
columns are named by the options below, by number or by header name.

Usage:

    subst_exp_levs <source> <target> [params=<what>] [--matchby <how>]

    an experimental source needs, and takes:
        --ecol <col>       energy         --units cm-1|kK|eV
        --jcol <col>       J              --skip <n>   header lines
        --pcol <col>       parity: o|odd|u|1|*, e|even|g|0
        --excludecol <col> exclude a level: anything non-blank, non-zero
        --parity-sets o,e  which parity each of the deck's sets is
        --confcol  <col> --termcol  <col> --perccol  <col>
        --conf1col <col> --term1col <col>
        --conf2col <col> --term2col <col> --perc2col <col>
        --conf3col <col> --term3col <col> --perc3col <col>
        --altcol <col> --altvalue or

`<source>` is where the observed energies are and `<target>` is the
deck they are being put into; the target is rewritten in place, with
what it held kept as `<target>.bak`.  A column argument is a 1-based
column number, or a column's name when the file has a header line.

`params=<what>` says how much of an **RCE** source's parameter table to
bring across as well; an experimental source has none:

    (nothing)  every value and every least-squares flag
    flags      the flags only, leaving the values as they are
    noflags    the values only, leaving the flags as they are
    Eav        the Eav values only;  F, G, Z, CI likewise, by name
    none       nothing: the parameter section is copied unchanged

`--matchby <how>` says how a source level is tied to a target level:

    heuristic  the composition where the source gives one, the energy
               where it does not -- the default
    Ecalc      the energy alone, against the target's fitted energies
    Eexp       the energy alone, against the target's observed ones,
               which is how a deck fitted early in an analysis has its
               observed energies replaced by optimized values

**The composition columns.**  The first triad is `(perc, conf, term)`,
but `conf1`/`term1` override the labels wherever either is given, for
the sheets that keep `conf` and `term` for bookkeeping and put the
component's real labels beside them.  The second and third triads are
`(perc2, conf2, term2)` and the same again; a blank configuration or
term in either means the one in the `conf` or `term` column.  With
`--altcol` given, a row whose value there is `--altvalue` is one whose
second triad is the *first* component written again in another coupling
scheme rather than a second component.

The percentages are read as **intensity** per cent, which is what a
published table quotes, and the deck's signed amplitudes are squared to
meet them.  A component with labels and no percentage is taken as 100.

A configuration is resolved against the deck's own labels through
**IN36**, where the configurations are written out in full: RCE
abbreviates a configuration from the left to fit six columns and no
rule recovers the missing digits from the label alone, so the card
whose subshells the source names is the card it means.  IN36 is looked
for beside the target deck and then in the directory above it; without
one, a source configuration falls back on matching the deck label its
compacted form ends with.

Only LS labels can be matched against an LS deck.  But a component
designated in JK or pair coupling still writes its **parent terms**
inside the configuration -- `4f2.(3F<2>).6p` beside `(2,3/2)*` -- and
the last of those is exactly the term RCE prints in brackets in its own
label, `f26p   (3F) 4D`.  So a component is matched on the
configuration, that parent and the final term separately, to whatever
depth the source gives them, and the target's components are added up
over the rest.  A row whose leading component gives neither an LS
parent nor an LS term gets no composition at all and is matched on
(parity, J, energy); the report says how many went each way.

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

from toolbox import altlevsrc, explevs, files, rce                # noqa: E402

OPTIONS = ('flags', 'noflags', 'eav', 'f', 'g', 'z', 'ci', 'none')


class Usage(Exception):
    """The command line will not do what was asked."""


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

#: The options that only describe an experimental source file.
SOURCE_OPTIONS = altlevsrc.Spec.NAMES + ('units', 'skip', 'altvalue',
                                         'parity_sets')


class Options:
    """What the command line asked for."""

    def __init__(self, argv):
        self.files = []
        self.params = None
        self.matchby = 'heuristic'
        self.source = {}

        rest = list(argv)
        while rest:
            arg = rest.pop(0)
            if not arg.startswith('--'):
                if arg.lower().startswith('params'):
                    self.params = arg.partition('=')[2].strip().lower()
                    if self.params not in OPTIONS:
                        raise Usage('params=%s is none of %s'
                                    % (self.params, ', '.join(OPTIONS)))
                else:
                    self.files.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.strip().lower().replace('-', '_')
            if not value:
                if not rest:
                    raise Usage('--%s needs a value'
                                % name.replace('_', '-'))
                value = rest.pop(0)
            if name == 'matchby':
                self.matchby = value.strip().lower()
                if self.matchby not in explevs.MATCH_BY:
                    raise Usage('--matchby %s is none of %s'
                                % (value, ', '.join(explevs.MATCH_BY)))
            elif name in SOURCE_OPTIONS:
                self.source[name] = value
            else:
                raise Usage('--%s is not an option' % name.replace('_', '-'))

        if len(self.files) != 2:
            raise Usage('a source file and a target deck are needed, and '
                        '%d file name%s given'
                        % (len(self.files),
                           ' was' if len(self.files) == 1 else 's were'))
        if 'skip' in self.source:
            self.source['skip'] = _number(self.source['skip'], '--skip')
        try:
            self.spec = altlevsrc.Spec(**self.source)
        except altlevsrc.AltSourceError as problem:
            # What the options say is part of the command line, so it
            # is checked here and reported as a usage error.
            raise Usage(str(problem))


def _number(text, name):
    try:
        value = int(text)
    except ValueError:
        raise Usage('%s wants a number, not %r' % (name, text))
    if value < 0:
        raise Usage('%s is not negative' % name)
    return value


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

def report(out, source, target, lost, certain=explevs.CERTAIN):
    for level in lost:
        out.write('Warning: no level of the target deck matches the source '
                  'level at parity set %d, J = %s, Ecalc = %g, Eobs = %g.\n'
                  % (level.kset, level.j, level.ec, level.ee))

    doubtful = explevs.suspicious(source, certain)
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

def is_deck(path):
    """Whether the source is an RCE deck rather than a level list.

    The old test asked whether the file had exactly ten tab-separated
    fields, which is the layout of one version of one program's output
    and nothing else.  The question worth asking is the other way
    round: an RCE deck is a file `rce.read_rceout` can read, and
    everything else is a level list to be described by option.
    """
    try:
        rce.read_rceout(path)
        return True
    except (rce.RceError, UnicodeDecodeError):
        return False


def main(argv):
    if not argv or '--help' in argv or '-h' in argv or '/?' in argv:
        sys.stdout.write(__doc__.split('Usage:')[1])
        return 0

    options = Options(argv)
    source_path, target_path = options.files
    deck = is_deck(source_path)
    if deck and options.source:
        raise Usage('%s is an RCE deck, and %s describe%s a level list that '
                    'is not one'
                    % (source_path,
                       ', '.join('--' + name.replace('_', '-')
                                 for name in sorted(options.source)),
                       '' if len(options.source) > 1 else 's'))
    if not deck and options.params is not None:
        raise Usage('%s is not an RCE deck, so it has no parameter table for '
                    'params=%s to take' % (source_path, options.params))

    target = explevs.read_deck(target_path)
    if deck:
        source = explevs.read_deck(source_path)
        certain = explevs.CERTAIN
    else:
        spec = options.spec
        source, notes = altlevsrc.read_alt_lev_src(source_path, spec, target)
        altlevsrc.normalize_target(target, notes.use_conf)
        notes.write(sys.stdout)
        certain = explevs.CERTAIN_INTENSITY

    lost = explevs.identify(source, target, options.matchby, certain)

    lines = files.read_lines(target_path, strip=False)
    given = 0
    for level in target.levels():
        if level.assigned is None or level.line < 0:
            continue
        lines[level.line] = rce.set_observed(lines[level.line], level.assigned)
        given += 1

    parameters = 0
    if deck:
        lines, parameters = substitute_parameters(
            lines, files.read_lines(source_path, width=140),
            Parameters(options.params))

    files.write_deck(target_path, lines, files.line_ending(target_path),
                     backup=os.path.basename(target_path) + '.bak')
    report(sys.stdout, source, target, lost, certain)
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
    except (altlevsrc.AltSourceError, explevs.ExpLevelsError, rce.RceError,
            OSError) as problem:
        sys.stderr.write('subst_exp_levs: %s\n' % problem)
        sys.exit(1)
