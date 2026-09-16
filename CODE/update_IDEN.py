#!/usr/bin/env python3
"""Carry a level optimization's result back into the IDEN2 working files.

A round of the analysis goes out of IDEN2 as a set of line assignments,
through the classification suite and LOPT, and comes back as levels that
have moved and lines that were, or were not, used in the fit.  This puts
that result into the two files IDEN2 reads, so that the next look at a
questionable level is taken against the current answer rather than the
one the files were built with:

  1. every **fixed** level of `ENLEV.DAT` gets the energy LOPT gives it
     and, as its uncertainty, `max(D1, D2tot)` -- the larger of the
     purely statistical uncertainty and the total one relative to the
     ground level, which is what a level is worth when it is looked at
     on its own;
  2. every level that is *not* fixed -- a purely calculated one, with
     a blank flag -- gets an uncertainty that says how far the
     calculation can be wrong for a level of its kind, since that
     uncertainty is the window IDEN2 searches around it.  The measure is
     the root-mean-square of `E_fixed - E_calc`, the observed energy
     minus the calculated one, over the experimental levels that share
     the calculated level's **configuration**: the first field of its
     `ENLEV.DAT` label, the part before the `_` or `~` that opens the
     term.  A configuration with at least one experimental level is
     *known*, and every calculated level of it is given that
     root-mean-square.  A configuration with none is *unknown*: there is
     nothing to measure it by, so its calculated levels keep the
     uncertainty the file gives them unless `--ucalc` names one;
  3. `TRANS.DAT` gets the same energies, and the lines assigned to
     transitions between two fixed levels are made to be exactly the
     lines the optimization used.

Usage:

    update_IDEN <LOPT levels> [<LOPT lines>] [options]

        --ucalc <value>   uncertainty, in cm-1, for the calculated levels
                          of the unknown configurations -- the ones with
                          no experimental level to measure one from
        --nosd            leave the calculated levels of the known
                          configurations alone too, so that --ucalc
                          applies to every calculated level
        --ids <file>      the translation table
                          (default: IDEN_level_ids.txt)
        --wntol <value>   how far a LOPT wavenumber may lie from its row
                          of dlv.dat and still be that row, in cm-1
                          (default: 0.05)

`<LOPT levels>` is the tab-delimited level table LOPT writes, read by its
`Designation`, `Energy`, `D1` and `D2tot` columns.  `<LOPT lines>` is its
line table, read by `L1`, `L2`, `wn_o` and `Weight`; leaving it out does
steps 1 and 2 and leaves `TRANS.DAT` alone.

`ENLEV.DAT`, `TRANS.DAT` and `dlv.dat` are taken from the working
directory under the names IDEN2 gives them, and each rewritten file's
previous contents are kept beside it as `<name>.bak`.  `dlv.dat` has to
be the one `TRANS.DAT` was assigned against: the assignment is a row
number of it and means nothing against another list.  Every `wn_o` of the
line table is looked for in it before anything is written, and the count
and the worst offset are reported.

LOPT names a level by its identifier and IDEN2 by a number of its own, so
the two are tied through `IDEN_level_ids.txt`, which `write_IDEN_level_ids`
writes.  **That table has to have been made from this `ENLEV.DAT`**: the
numbering is a position in one file's ordering, and rebuilding the file
with `conv_out` moves every level in it.

A line is taken to have been used in the optimization when its `Weight`
is above zero.  That is the flag to rely on rather than the `F` column:
`P`, `M` and `virt` all mean a line was excluded, but the weight says so
for the excluded line that carries no letter as well, and a weight that
rounds to zero was worth nothing to the fit whatever it was before the
rounding.  A transition between two fixed levels that the optimization
did not use loses its assignment; one it used gets the row of `dlv.dat`
that its observed wavenumber falls on.

Transitions that reach only one fixed level, or none, are not the
optimization's to speak for and keep the assignment they have; their
energies are brought up to date like everything else.
"""

import csv
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from toolbox import files, iden, levids                           # noqa: E402


class Usage(Exception):
    """The command line will not do what was asked."""


#: The columns each of the two LOPT tables is read by.
LEVEL_COLUMNS = ('Designation', 'Energy', 'D1', 'D2tot')
LINE_COLUMNS = ('L1', 'L2', 'wn_o', 'Weight')

#: How far a wavenumber quoted by LOPT may lie from the row of `dlv.dat`
#: it stands for.  LOPT prints an observed wavenumber rounded to the
#: number of decimals its uncertainty deserves, so the two disagree in
#: the last place and no further; this is wide enough for that and far
#: short of the spacing of two resolved lines.
WNTOL = 0.05


class Options:
    """What the command line asked for."""

    def __init__(self, argv):
        self.files = []
        self.ucalc = None
        self.spread = True
        self.ids = levids.IDEN_IDS
        self.wntol = WNTOL

        rest = list(argv)
        while rest:
            arg = rest.pop(0)
            if not arg.startswith('--'):
                self.files.append(arg)
                continue
            name, _, value = arg[2:].partition('=')
            name = name.strip().lower().replace('-', '_')
            if name in ('nosd', 'no_sd'):
                if value:
                    raise Usage('--%s takes no value' % name)
                self.spread = False
                continue
            if not value:
                if not rest:
                    raise Usage('--%s needs a value' % name)
                value = rest.pop(0)
            if name == 'ids':
                self.ids = value
            elif name == 'ucalc':
                self.ucalc = _positive(value, '--ucalc')
            elif name == 'wntol':
                self.wntol = _positive(value, '--wntol')
            else:
                raise Usage('--%s is not an option' % name)

        if not 0 < len(self.files) < 3:
            raise Usage('a LOPT level table and, if the transitions are to be '
                        'done, a LOPT line table are what this takes, and %d '
                        'file name%s given'
                        % (len(self.files),
                           ' was' if len(self.files) == 1 else 's were'))
        self.levels_path = self.files[0]
        self.lines_path = self.files[1] if len(self.files) > 1 else None


def _positive(text, name):
    try:
        value = float(text)
    except ValueError:
        raise Usage('%s wants a number, not %r' % (name, text))
    if value <= 0:
        raise Usage('%s wants a positive number, and %g is not one'
                    % (name, value))
    return value


class TableError(Exception):
    """A LOPT table that cannot be read."""


def read_table(path, wanted):
    """The rows of a tab-delimited LOPT table, as `{column: field}`.

    The header is the first row, and the columns are found in it by name
    rather than by position: LOPT's tables carry columns this has no use
    for, and their number is not fixed.
    """
    with open(path, newline='', errors='replace') as handle:
        rows = list(csv.reader(handle, delimiter='\t'))
    if not rows:
        raise TableError('%s holds no rows' % path)
    header = [field.strip() for field in rows[0]]
    lowered = [field.lower() for field in header]
    at = {}
    for name in wanted:
        try:
            at[name] = lowered.index(name.lower())
        except ValueError:
            raise TableError('%s has no %s column; it names %s'
                             % (path, name, ', '.join(header) or 'nothing'))
    out = []
    for number, row in enumerate(rows[1:], 2):
        if not any(field.strip() for field in row):
            continue
        one = dict((name, row[i].strip() if i < len(row) else '')
                   for name, i in at.items())
        one['line'] = number
        out.append(one)
    if not out:
        raise TableError('%s holds no rows below its header' % path)
    return out


def _float(row, name, path):
    try:
        return float(row[name])
    except ValueError:
        raise TableError('%s line %d: %r is not a number in the %s column'
                         % (path, row['line'], row[name], name))


# --------------------------------------------------------------------------
# how far the calculation is wrong, configuration by configuration
# --------------------------------------------------------------------------

def split_label(label):
    """`("configuration", "term")` of an `ENLEV.DAT` label.

    A label is written as the configuration in IDEN2's own abbreviation,
    then a separator -- `_` for a level named in LS coupling, `~` for one
    named in jj -- then the term, as in `f25g _3H4G` or `p5f3d_4I6Ga`.
    Neither separator occurs inside either part, so the first one found
    divides them; the term keeps its separator, because that is how
    IDEN2 writes and sorts these labels.  A label with no separator at
    all is taken to be a configuration and nothing else.  The field is
    padded to a fixed width but not always by the same amount, so the
    split is made on the separator and never on a column.
    """
    for at, character in enumerate(label):
        if character in '_~':
            return label[:at].strip(), label[at:].strip()
    return label.strip(), ''


class Spread:
    """How far the calculation misses, over one configuration."""

    def __init__(self, residuals):
        self.count = len(residuals)
        self.rms = math.sqrt(sum(r * r for r in residuals) / len(residuals))
        self.worst = max(abs(r) for r in residuals)

    @property
    def usable(self):
        """Is this a window IDEN2 can search?

        A root-mean-square of zero is not.  It happens when a
        configuration's only experimental level was never moved off its
        calculated position, and a level given a zero uncertainty is one
        IDEN2 searches no window at all around.
        """
        return self.rms > 0


def configuration_spread(levels):
    """`{configuration: Spread}` over the experimental levels.

    The residual of one level is `E_fixed - E_calc`, its observed energy
    minus its calculated one, and only the fixed levels have both: for a
    level that has not been found, `E_fixed` is a copy of `E_calc` and
    the residual would be zero by construction, which would pull every
    root-mean-square towards nothing.  Being fixed is read from the star
    that IDEN2 sets, which is the same flag step 1 above works from.
    """
    residuals = {}
    for one in levels:
        if not one.fixed:
            continue
        configuration = split_label(one.label)[0]
        if configuration:
            residuals.setdefault(configuration, []).append(
                one.e_fixed - one.e_calc)
    return dict((configuration, Spread(found))
                for configuration, found in residuals.items())


# --------------------------------------------------------------------------
# the levels
# --------------------------------------------------------------------------

def update_levels(levels, rows, ids, path, out):
    """Step 1: the fixed levels take the energy the optimization gives them.

    The uncertainty that goes with the energy is `max(D1, D2tot)`, the
    larger of the purely statistical uncertainty and the total one
    relative to the ground level, which is what a level is worth when it
    is looked at on its own.
    """
    by_number = dict((one.number, one) for one in levels)
    named, unknown, absent = set(), [], []

    for row in rows:
        ident = row['Designation']
        if not ident:
            continue
        number = ids.get(ident)
        if number is None:
            unknown.append(ident)
            continue
        level = by_number.get(number)
        if level is None:
            raise TableError(
                'the translation table gives %s the IDEN2 level number %d, '
                'which %s has not got; the table was made from a different '
                'ENLEV.DAT' % (ident, number, iden.ENLEV))
        if not level.fixed:
            absent.append(ident)
            continue
        level.update(e_fixed=_float(row, 'Energy', path),
                     uncertainty=max(_float(row, 'D1', path),
                                     _float(row, 'D2tot', path)))
        named.add(number)

    fixed = [one for one in levels if one.fixed]
    missed = [one for one in fixed if one.number not in named]
    zero = [one for one in fixed if one.number in named and not one.uncertainty]

    out.write('%s: %d of the %d fixed levels of %s brought up to date.\n'
              % (path, len(named), len(fixed), iden.ENLEV))
    if unknown:
        out.write('   %d level%s of the table not in %s and left out'
                  % (len(unknown), '' if len(unknown) == 1 else 's',
                     levids.IDEN_IDS))
        out.write(': %s%s.\n' % (', '.join(unknown[:5]),
                                 ', ...' if len(unknown) > 5 else ''))
    if absent:
        out.write('   Warning: %d level%s of the table name%s a level of %s '
                  'that is not fixed: %s.\n'
                  % (len(absent), '' if len(absent) == 1 else 's',
                     's' if len(absent) == 1 else '', iden.ENLEV,
                     ', '.join(absent[:5])))
    for one in missed[:5]:
        out.write('   Warning: fixed level %d (%s) is in %s and not in the '
                  'LOPT table; it keeps the energy it had.\n'
                  % (one.number, one.label.strip(), iden.ENLEV))
    if len(missed) > 5:
        out.write('   ... and %d more fixed levels the LOPT table does not '
                  'name.\n' % (len(missed) - 5))
    if zero:
        out.write('   Warning: %d level%s came out with a zero uncertainty '
                  '(%s); IDEN2 searches no window at all around such a '
                  'level.\n'
                  % (len(zero), '' if len(zero) == 1 else 's',
                     ', '.join(str(one.number) for one in zero[:5])))

    return by_number


#: A root-mean-square resting on fewer experimental levels than this is
#: quoted with a warning: it is one or two numbers standing for a whole
#: configuration, and the next level found in that configuration can
#: change it several-fold.
THIN = 3


def update_calculated(levels, spread, ucalc, out):
    """Step 2: the window IDEN2 will search around each calculated level.

    A calculated level of a **known** configuration -- one with at least
    one experimental level to measure by -- is given the
    root-mean-square of `E_fixed - E_calc` over the experimental levels
    of that same configuration.  A calculated level of an **unknown**
    configuration has nothing to measure it by: it keeps the uncertainty
    the file already gives it, and takes `ucalc` only when `ucalc` was
    asked for.  `spread` of None is the case where the measurement was
    refused altogether (`--nosd`), and then every calculated level counts
    as unknown.
    """
    refused = spread is None
    spread = spread or {}
    measured, left = {}, {}
    for one in levels:
        if one.fixed:
            continue
        configuration = split_label(one.label)[0]
        found = spread.get(configuration)
        if found is not None and found.usable:
            one.update(uncertainty=found.rms)
            measured.setdefault(configuration, []).append(one)
        else:
            left.setdefault(configuration, []).append(one)

    if measured:
        total = sum(len(group) for group in measured.values())
        out.write('   %d calculated level%s of %d known configuration%s given '
                  'the rms of E_obs - E_calc of %s own configuration:\n'
                  % (total, '' if total == 1 else 's', len(measured),
                     '' if len(measured) == 1 else 's',
                     'its' if total == 1 else 'their'))
        out.write('      %-12s %6s %10s %10s %6s\n'
                  % ('config', 'exp.', 'rms', 'max|O-C|', 'calc.'))
        for configuration in sorted(measured, key=lambda c: -spread[c].rms):
            found = spread[configuration]
            out.write('      %-12s %6d %10.2f %10.2f %6d%s\n'
                      % (configuration, found.count, found.rms, found.worst,
                         len(measured[configuration]),
                         '   thin' if found.count < THIN else ''))
        thin = [c for c in measured if spread[c].count < THIN]
        if thin:
            out.write('      Warning: %d of these rms values rest on fewer '
                      'than %d experimental levels, and are that rough.\n'
                      % (len(thin), THIN))

    if not left:
        return
    total = sum(len(group) for group in left.values())
    if refused:
        why = ('%d calculated level%s of %d configuration%s, none of them '
               'measured because --nosd was given'
               % (total, '' if total == 1 else 's', len(left),
                  '' if len(left) == 1 else 's'))
    else:
        why = ('%d calculated level%s of %d configuration%s with no '
               'experimental level to measure by'
               % (total, '' if total == 1 else 's', len(left),
                  '' if len(left) == 1 else 's'))
    if ucalc is None:
        out.write('   %s, and %s the uncertainty %s had; --ucalc is what sets '
                  'one: %s.\n'
                  % (why, 'keeps' if total == 1 else 'keep',
                     'it' if total == 1 else 'they', _configurations(left)))
        return
    for group in left.values():
        for one in group:
            one.update(uncertainty=ucalc)
    out.write('   %s, given the --ucalc uncertainty of %g cm-1: %s.\n'
              % (why, ucalc, _configurations(left)))


def _configurations(groups):
    """The configuration names of a `{configuration: levels}`, listed."""
    names = sorted(name or '(no label)' for name in groups)
    if len(names) <= 8:
        return ', '.join(names)
    return '%s, ... and %d more' % (', '.join(names[:8]), len(names) - 8)


# --------------------------------------------------------------------------
# the transitions
# --------------------------------------------------------------------------

class Used:
    """One transition the optimization used, and the line it used."""

    def __init__(self, weight, wavenumber):
        self.weight = weight
        self.wavenumber = wavenumber


def check_wavenumbers(rows, nearest, tolerance, path, out):
    """Every `wn_o` of the line table against the rows of `dlv.dat`.

    The assignment `update_trans` writes is a row number of `dlv.dat`, so
    the two files have to be the pair that belong together: the `dlv.dat`
    beside `TRANS.DAT` must be the one the line list was measured from.
    That is checked here, over the whole table -- every row, used or not,
    and whether or not its levels are in the translation table -- rather
    than only over the lines that reach a transition `TRANS.DAT` lists,
    because a wrong `dlv.dat` would otherwise show up as nothing worse
    than a few unassigned lines.

    Nothing is refused on the strength of it.  The worst offset is
    reported so that a tolerance set too tight is visible as a small
    number of near misses, and a mismatched file as a large one.
    """
    missing, worst, where, unreadable = 0, 0.0, None, 0
    for row in rows:
        try:
            wavenumber = float(row['wn_o'])
        except ValueError:
            unreadable += 1
            continue
        at, apart = nearest(wavenumber)
        if apart is None or apart > tolerance:
            missing += 1
        if apart is not None and apart > worst:
            worst, where = apart, wavenumber

    out.write('dlv.dat: %d of the %d wavenumbers in %s fall on a row of it '
              'within %g cm-1.\n'
              % (len(rows) - missing - unreadable, len(rows), path,
                 tolerance))
    if where is not None:
        out.write('   the furthest is %.4f cm-1, %.4f cm-1 from its nearest '
                  'row.\n' % (where, worst))
    if unreadable:
        out.write('   Warning: %d row%s no wn_o that reads as a number.\n'
                  % (unreadable, ' has' if unreadable == 1 else 's have'))
    if missing:
        out.write('   Warning: %d wavenumber%s no row of dlv.dat within the '
                  'tolerance.  Is this the dlv.dat the line list was made '
                  'from?  --wntol widens the search.\n'
                  % (missing, ' has' if missing == 1 else 's have'))
    return missing


def read_used(rows, ids, path, out):
    """`{(lower, upper): Used}` over the IDEN2 numbering, and what was seen.

    A transition can appear more than once -- the same pair of levels
    measured in two spectra, or a blend entered twice -- so the row with
    the largest weight is the one that speaks for it.  `seen` holds every
    pair the table names at all, used or not, because a pair that is in
    the table with a zero weight and one that is not in it are the same
    thing here: neither was used.
    """
    used, seen, unknown, unweighted = {}, set(), 0, 0
    for row in rows:
        lower = ids.get(row['L1'])
        upper = ids.get(row['L2'])
        if lower is None or upper is None:
            unknown += 1
            continue
        seen.add((lower, upper))
        try:
            weight = float(row['Weight'])
        except ValueError:
            unweighted += 1
            continue
        if weight <= 0:
            continue
        wavenumber = _float(row, 'wn_o', path)
        one = used.get((lower, upper))
        if one is None or weight > one.weight:
            used[(lower, upper)] = Used(weight, wavenumber)

    out.write('%s: %d transitions between levels of %s, %d of them used in '
              'the optimization.\n'
              % (path, len(seen), levids.IDEN_IDS, len(used)))
    if unknown:
        out.write('   %d row%s a level the translation table does not have '
                  'and %s left out.\n'
                  % (unknown, ' names' if unknown == 1 else 's name',
                     'is' if unknown == 1 else 'are'))
    if unweighted:
        out.write('   Warning: %d row%s no weight and %s taken as unused.\n'
                  % (unweighted, ' has' if unweighted == 1 else 's have',
                     'is' if unweighted == 1 else 'are'))
    return used


def update_trans(records, by_number, used, observed, find, out):
    """Step 3: the energies, and the assignments of the fixed-fixed lines."""
    counted = dict(kept=0, added=0, moved=0, cleared=0, unresolved=0,
                   other=0, orphan=0)
    reached = set()

    for record in records:
        if record.kind == '$':
            level = by_number.get(record.upper)
            if level is None:
                counted['orphan'] += 1
                continue
            record.set_level_energy(level.e_fixed, level.fixed)
            continue
        if record.kind != '+':
            continue

        upper = by_number.get(record.upper)
        lower = by_number.get(record.lower)
        if upper is None or lower is None:
            counted['orphan'] += 1
            continue

        computed = abs(upper.e_fixed - lower.e_fixed)
        record.set_lower(lower.e_fixed, lower.fixed)
        record.set_transition(computed)

        if not (upper.fixed and lower.fixed):
            # Neither the optimization nor this tool has anything to say
            # about a transition that does not run between two observed
            # levels: it can carry no Ritz wavenumber, so it was never in
            # the fit, and whatever IDEN2's own session assigned to it
            # stands.
            counted['other'] += 1
            record.refresh(computed)
            continue

        one = used.get((record.lower, record.upper))
        reached.add((record.lower, record.upper))
        if one is None:
            if record.line:
                counted['cleared'] += 1
            record.assign(0)
            continue
        at = find(one.wavenumber)
        if not at:
            # The line the optimization used is in no row of dlv.dat, so
            # there is no number to write.  Saying nothing would leave
            # the record looking like one the optimization dropped, so
            # what it has is kept and the case is reported.
            counted['unresolved'] += 1
            record.refresh(computed)
            continue
        if at == record.line:
            counted['kept'] += 1
        elif record.line:
            counted['moved'] += 1
        else:
            counted['added'] += 1
        record.assign(at, observed[at - 1], computed)

    out.write('%s: %d assignments kept, %d added, %d moved, %d cleared.\n'
              % (iden.TRANS, counted['kept'], counted['added'],
                 counted['moved'], counted['cleared']))
    out.write('   %d transitions do not run between two fixed levels and '
              'keep what they had.\n' % counted['other'])
    if counted['unresolved']:
        out.write('   Warning: %d used line%s no row of dlv.dat within the '
                  'wavenumber tolerance; %s assignment %s left as it was.\n'
                  % (counted['unresolved'],
                     ' has' if counted['unresolved'] == 1 else 's have',
                     'its' if counted['unresolved'] == 1 else 'their',
                     'is' if counted['unresolved'] == 1 else 'was'))
    if counted['orphan']:
        out.write('   Warning: %d record%s a level %s has not got.\n'
                  % (counted['orphan'], ' names' if counted['orphan'] == 1
                     else 's name', iden.ENLEV))
    beyond = set(used) - reached
    if beyond:
        out.write('   %d of the lines the optimization used are on '
                  'transitions %s does not list, and cannot be assigned '
                  'there; %s writes only the transitions above its intensity '
                  'threshold.\n' % (len(beyond), iden.TRANS, 'conv_out'))
    return counted


# --------------------------------------------------------------------------

def main(argv):
    if not argv or '--help' in argv or '-h' in argv or '/?' in argv:
        sys.stdout.write(__doc__.split('Usage:')[1])
        return 0

    options = Options(argv)
    out = sys.stdout

    ids = levids.read_iden_ids(options.ids)
    enlev = files.resolve(iden.ENLEV, '.')
    levels = iden.read_enlev(enlev)
    rows = read_table(options.levels_path, LEVEL_COLUMNS)
    by_number = update_levels(levels, rows, ids, options.levels_path, out)
    # Step 2 measures the configurations on the energies step 1 has just
    # written, so that the window a calculated level is given reflects
    # the optimization this run is carrying back and not the one before
    # it.
    update_calculated(levels,
                      configuration_spread(levels) if options.spread else None,
                      options.ucalc, out)
    iden.write_enlev(levels, enlev,
                     backup=os.path.basename(str(enlev)) + '.bak')
    out.write('%s rewritten; what it held is %s.bak.\n'
              % (enlev, os.path.basename(str(enlev))))

    if options.lines_path:
        lines = read_table(options.lines_path, LINE_COLUMNS)
        observed = iden.read_dlv(files.resolve('dlv.dat', '.'))
        nearest = iden.dlv_nearest(observed)
        check_wavenumbers(lines, nearest, options.wntol,
                          options.lines_path, out)
        used = read_used(lines, ids, options.lines_path, out)
        find = iden.dlv_index(observed, options.wntol)
        trans = files.resolve(iden.TRANS, '.')
        records = iden.read_trans(trans)
        update_trans(records, by_number, used, observed, find, out)
        iden.write_trans(records, trans,
                         backup=os.path.basename(str(trans)) + '.bak')
        out.write('%s rewritten; what it held is %s.bak.\n'
                  % (trans, os.path.basename(str(trans))))
    return 0


if __name__ == '__main__':
    try:
        sys.exit(main(sys.argv[1:]))
    except Usage as problem:
        sys.stderr.write('update_IDEN: %s\n\n%s\n'
                         % (problem, __doc__.split('Usage:')[1]))
        sys.exit(2)
    except (TableError, levids.LevelIdError, iden.IdenError,
            OSError) as problem:
        sys.stderr.write('update_IDEN: %s\n' % problem)
        sys.exit(1)
