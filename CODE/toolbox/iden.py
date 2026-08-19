"""The two input files of the IDEN2 line-identification code.

IDEN2 takes a level list (`ENLEV.DAT`) and, under each level, the lines
that leave it (`TRANS.DAT`), with each calculated line tied to an
observed one where the two agree inside the observation's uncertainty.
Making that connection is what this module does:

  * `dlv.dat` is the observed line list -- intensity, wavenumber, an
    uncertainty, and the columns IDEN2 keeps for itself;
  * `numset.dat` divides it into wavelength ranges, each with a default
    uncertainty for the lines that carry none of their own;
  * every calculated line between two *observed* levels is looked up in
    that list, and the observed line nearest in wavenumber, inside the
    uncertainty, is the one it is assigned to.

Both files are fixed-column and are read by the columns
`conv_out.pl:ReadLines` and `:ReadNumset` count.  The observed list runs
in **descending wavenumber** -- ascending wavelength -- which is what
the search here relies on, as the Perl's did.

IDEN2 numbers its levels from the top down, so the numbering here is not
the `Lev_#` of the line list; the two are kept apart deliberately.
"""

import math
import shutil
from dataclasses import dataclass, field

from . import fortfmt, outg11

#: A line assigned to an observation keeps only the transitions within
#: this fraction of the strongest one assigned to it: a much weaker
#: candidate is not a competing identification (`conv_out.pl:$Athresh`).
THRESHOLD = 0.1

#: IDEN2's level number is I4 and its arrays are sized for this many.
MAX_LEVELS = 9999


class IdenError(Exception):
    """A file IDEN2 needs is missing or does not read."""


# --------------------------------------------------------------------------
# the observed line list
# --------------------------------------------------------------------------

@dataclass
class Observed:
    """One line of `dlv.dat`."""

    intensity: float
    wavenumber: float
    uncertainty: float          # of the wavelength, in angstroms
    middle: str                 # the 34 columns IDEN2 keeps for itself
    tail: str
    assigned: list = field(default_factory=list)

    def spread(self, default):
        """The wavenumber uncertainty this line is matched within.

        The wavelength uncertainty is the one that is quoted, so it is
        turned into a wavenumber one at this line's own wavenumber.
        """
        one = self.uncertainty or default
        if not one:
            raise IdenError(
                'the line at %.3f cm-1 has no uncertainty in dlv.dat and its '
                'range in numset.dat gives no default one' % self.wavenumber)
        return one * self.wavenumber / (1e8 / self.wavenumber)


def read_dlv(path='dlv.dat'):
    """The observed line list, in the order the file has it."""
    out = []
    try:
        handle = open(path, 'r', errors='replace')
    except OSError:
        raise IdenError('%s is needed for the IDEN2 files and is not here'
                        % path)
    with handle:
        for at, raw in enumerate(handle, 1):
            line = raw.rstrip('\n').rstrip('\r')
            if not line.strip():
                continue
            if len(line) < 20:
                raise IdenError('%s line %d is too short to be a line list row'
                                % (path, at))
            wide = line.ljust(60)
            try:
                intensity = float(wide[:5] or 0)
                wavenumber = float(wide[5:19] or 0)
                uncertainty = abs(float(wide[53:60] or 0))
            except ValueError:
                raise IdenError('%s line %d does not read in the fixed columns '
                                'of a dlv.dat row' % (path, at))
            if wavenumber <= 0:
                raise IdenError('%s line %d has no wavenumber' % (path, at))
            out.append(Observed(intensity=abs(intensity),
                                wavenumber=wavenumber,
                                uncertainty=uncertainty,
                                middle=wide[19:53], tail=line[60:]))
    if not out:
        raise IdenError('%s holds no lines' % path)
    return out


def write_dlv(lines, path='dlv.dat'):
    """Rewrite the observed list with its intensities on a log scale.

    IDEN2 compares intensities directly, so a list whose intensities run
    over several decades has to be put on the same footing as the
    calculated ones, which are logarithmic.  The original is kept beside
    it as `dlv.dat.bak`.
    """
    biggest = max(one.intensity for one in lines)
    if 999 - 10 * math.log(biggest) < 0:
        raise IdenError('the largest intensity in %s, %g, is too large to put '
                        'on a log scale' % (path, biggest))
    shutil.copyfile(path, path + '.bak')
    with open(path, 'w') as out:
        for one in lines:
            one.intensity = 10 * math.log(one.intensity)
            uncertainty = ('%7.4f' % one.uncertainty
                           if one.uncertainty else ' ')
            out.write('%5d%14.3f%34s%7s%s\n'
                      % (one.intensity, one.wavenumber, one.middle,
                         uncertainty, one.tail))


# --------------------------------------------------------------------------
# the wavelength ranges
# --------------------------------------------------------------------------

@dataclass
class Range:
    """One row of `numset.dat`: where a range starts and how good it is."""

    label: str
    wavelength: float           # where the range starts, in angstroms
    uncertainty: float          # the default for lines that carry none
    first: int = 0              # the dlv.dat line it starts at


def read_numset(path='numset.dat'):
    out = []
    try:
        handle = open(path, 'r', errors='replace')
    except OSError:
        raise IdenError('%s is needed for the IDEN2 files and is not here'
                        % path)
    with handle:
        for at, raw in enumerate(handle, 1):
            line = raw.rstrip('\n').rstrip('\r')
            if not line.strip() or line.startswith('$'):
                continue
            wide = line.ljust(22)
            try:
                wavelength = float(wide[5:15] or 0)
                uncertainty = float(wide[15:22] or 0)
            except ValueError:
                raise IdenError('%s line %d does not read in the fixed columns '
                                'of a numset.dat row' % (path, at))
            if wavelength <= 0:
                raise IdenError('%s line %d has no starting wavelength'
                                % (path, at))
            if uncertainty <= 0:
                raise IdenError('%s line %d has no default uncertainty'
                                % (path, at))
            out.append(Range(label=wide[:5].strip(), wavelength=wavelength,
                             uncertainty=uncertainty))
    if not out:
        raise IdenError('%s holds no ranges' % path)
    return out


def locate(ranges, lines, to_vacuum):
    """Find where each range starts in the observed list.

    A range is given by a wavelength, in air where a published list
    quotes air, so it is converted before it is compared; the line it
    starts at is the first one at or below that wavenumber, or the first
    one that agrees with it inside its own uncertainty.  A range whose
    first line carries no uncertainty of its own is given the range's.
    """
    default = None
    for one in ranges:
        wavelength = one.wavelength
        if 2000 < wavelength < 20000:
            wavelength = to_vacuum(wavelength)
        start = 1e8 / wavelength
        default = one.uncertainty
        first = 0
        for at, line in enumerate(lines, 1):
            if line.wavenumber <= start:
                first = at
                break
            spread = line.wavenumber * max(line.uncertainty, default) \
                / (1e8 / line.wavenumber)
            if not spread:
                raise IdenError(
                    'neither numset.dat range %s nor dlv.dat line %d gives an '
                    'uncertainty' % (one.label, at))
            if abs(line.wavenumber - start) <= spread:
                first = at
                break
        first = first or len(lines)
        if lines[first - 1].uncertainty:
            default = lines[first - 1].uncertainty
        else:
            lines[first - 1].uncertainty = default
        one.first = first
        one.uncertainty = default
    return ranges


def find(lines, wavenumber, first, default):
    """The observed line a calculated wavenumber belongs to, or 0.

    The list runs downward in wavenumber from `first`, so this is a
    bisection over a descending sequence, and the nearer of the two
    lines it closes on wins provided it is inside its uncertainty.
    """
    first = first or 1
    last = len(lines)
    start = lines[first - 1]
    if wavenumber > start.wavenumber:
        return first if abs(wavenumber - start.wavenumber) \
            <= start.spread(default) else 0

    while last > first:
        step = (last - first) // 2
        if not step:
            break
        middle = first + step
        if lines[middle - 1].wavenumber < wavenumber:
            last = middle
        else:
            first = middle
            if abs(wavenumber - lines[first - 1].wavenumber) \
                    <= lines[first - 1].spread(default):
                return first

    near, far = lines[first - 1], lines[last - 1]
    best, spread = abs(wavenumber - near.wavenumber), near.spread(default)
    other = abs(wavenumber - far.wavenumber)
    at = first
    if other <= far.spread(default) and (other < best or best > spread):
        at, best, spread = last, other, far.spread(default)
    return 0 if best > spread else at


# --------------------------------------------------------------------------
# writing the two files
# --------------------------------------------------------------------------

@dataclass
class IdenLevel:
    """One level as IDEN2 numbers and names it."""

    number: int
    lev: object                 # conv_out.Lev
    label: str


def _label(case, coupling, one):
    """`sp2  _3P4P`: the configuration, the parity, the term, as IDEN2 wants.

    IDEN2's label is eleven characters and is compared as text, so it is
    built from the ING11 label and the eigenvector row's own state
    label rather than from the spectroscopic designation -- the same
    thing `conv_out.pl:Create_IDEN_lev_map` does.
    """
    block = case.blocks[(one.kset, one.j)]
    try:
        rows = block.section(coupling).rows
    except outg11.Outg11Error:
        # A block of one basis state gets no section but LS's; the level
        # is that one state, so its label is the one row LS printed.
        rows = block.section('LS').rows
    weights = case.amplitudes(coupling, one.kset, one.j, one.rcg)
    row = rows[weights[0][1] - 1]
    configuration, term = row.label, row.bracketed
    for ch in ' ()>]}':
        term = term.replace(ch, '')
    parity = '_' if case.card(one.kset, 1).parity == 0 else '~'
    return '%-5s%1s%-5s' % (configuration[:5], parity, term[:5])


def _iden_levels(case, coupling, order):
    """Every level, numbered from the highest energy down, as IDEN2 wants."""
    down = sorted(order, key=lambda one: (-one.energy, one.kset, one.j))
    return [IdenLevel(number=at, lev=one, label=_label(case, coupling, one))
            for at, one in enumerate(down, 1)]


def write(case, order, transitions, options):
    """Write ENLEV.DAT and TRANS.DAT for one run."""
    from . import vacair

    coupling = options.couplings[0]
    lines = read_dlv()
    if options.convert2log:
        write_dlv(lines)
    ranges = locate(read_numset(), lines, vacair.to_vacuum)

    numbered = _iden_levels(case, coupling, order)
    by_level = {(one.lev.kset, one.lev.j, one.lev.rcg): one
                for one in numbered}
    printed = [one for one in transitions if one.intensity >= 0.5]
    assigned = _assign(printed, lines, ranges)

    truncate = options.truncate
    if len(numbered) > MAX_LEVELS and not truncate:
        print('   %d levels is more than the %d IDEN2 can hold; truncating'
              % (len(numbered), MAX_LEVELS))
        truncate = True

    rows = _trans_rows(printed, by_level, assigned, lines, options)
    if truncate:
        numbered, rows = _truncate(numbered, rows)
    _write_enlev(numbered)
    _write_trans(numbered, rows, truncate)
    print('wrote ENLEV.DAT (%d levels) and TRANS.DAT (%d lines)'
          % (len(numbered), len(rows)))


def _assign(transitions, lines, ranges):
    """Tie every calculated line between two observed levels to an observation."""
    out = {}
    for one in transitions:
        if not (one.upper.observed and one.lower.observed):
            continue
        spread = abs(one.upper.e_obs - one.lower.e_obs) * 1000
        if not spread:
            continue
        chosen = None
        for a_range in ranges:
            start = lines[a_range.first - 1]
            if spread <= start.wavenumber + start.spread(a_range.uncertainty):
                chosen = a_range
                break
        if chosen is None:
            continue
        at = find(lines, spread, chosen.first, chosen.uncertainty)
        if at:
            lines[at - 1].assigned.append(one)
            out[id(one)] = at

    for line in lines:
        line.assigned.sort(key=lambda t: -t.intensity)
        keep = []
        for one in line.assigned:
            if keep and one.intensity < THRESHOLD * keep[0].intensity:
                break
            keep.append(one)
        for one in line.assigned[len(keep):]:
            out.pop(id(one), None)
        line.assigned = keep
    return out


def _trans_rows(transitions, by_level, assigned, lines, options):
    """One row per line: the two IDEN2 level numbers and the numbers IDEN2 wants."""
    biggest = max([one.intensity for one in lines], default=0.0)
    if options.linear:
        shift = 0.0
    else:
        shift = ((10 * math.log(biggest) if options.convert2log else biggest)
                 - 10 * math.log(options.max_intensity))
    if ((biggest > 99999 and options.linear)
            or (shift > 99999 and not options.linear)):
        raise IdenError(
            'the largest observed intensity, %g, is too large for IDEN2'
            % biggest)

    out = []
    for one in transitions:
        upper = by_level.get((one.upper.kset, one.upper.j, one.upper.rcg))
        lower = by_level.get((one.lower.kset, one.lower.j, one.lower.rcg))
        if upper is None or lower is None:
            continue
        intensity = one.intensity
        if intensity < 1:
            continue
        if options.linear:
            intensity *= biggest / options.max_intensity
        else:
            intensity = 10 * math.log(intensity) + shift
        at = assigned.get(id(one), 0)
        observation = lines[at - 1] if at else None
        computed = abs(one.upper.e_obs - one.lower.e_obs) * 1000
        out.append((upper.number, lower.number, intensity, one.lower,
                    computed,
                    observation.intensity if observation else 0.0,
                    observation.wavenumber if observation else 0.0,
                    (observation.wavenumber - computed) if observation else 0.0,
                    at))
    out.sort(key=lambda row: (row[0], row[1]))
    return out


def _truncate(numbered, rows):
    """Drop every level no printed line touches, and renumber."""
    wanted = set()
    for row in rows:
        wanted.add(row[0])
        wanted.add(row[1])
    kept = [one for one in numbered if one.number in wanted]
    if len(kept) > MAX_LEVELS:
        raise IdenError('after truncation there are still %d levels, more '
                        'than the %d IDEN2 can hold' % (len(kept), MAX_LEVELS))
    renumbered = {}
    out = []
    for at, one in enumerate(kept, 1):
        renumbered[one.number] = at
        out.append(IdenLevel(number=at, lev=one.lev, label=one.label))
    rows = [(renumbered[row[0]], renumbered[row[1]]) + row[2:] for row in rows]
    rows.sort(key=lambda row: (row[0], row[1]))
    return out, rows


def _write_enlev(numbered, path='ENLEV.DAT'):
    with open(path, 'w') as out:
        for one in numbered:
            lev = one.lev
            observed = lev.observed
            # The three energy fields are F12.3 in `conv_out.pl`, which
            # asterisks out above 10^9 cm-1 -- U90 reaches 8x10^8 and a
            # heavier ion passes it.  `write_f_wide` drops decimals and
            # then goes exponential instead, as the Fortran does in the
            # same situation: IDEN2 parses the field as a number, and
            # any of these is a number.
            out.write('%4d%s%10.3f%s%2s%s%5.1f /%11s/\n'
                      % (one.number,
                         fortfmt.write_f_wide(lev.e_calc * 1000, 12, 3),
                         10 if observed else 5000,
                         fortfmt.write_f_wide(lev.e_obs * 1000, 12, 3),
                         '*' if observed else ' ',
                         fortfmt.write_f_wide((lev.e_obs - lev.e_calc) * 1000,
                                              12, 3),
                         lev.j, one.label))


def _write_trans(numbered, rows, truncate, path='TRANS.DAT'):
    by_number = {one.number: one for one in numbered}
    with open(path, 'w') as out:
        at = 0
        for row in rows:
            (upper, lower, intensity, low, computed, observed_intensity,
             wavenumber, difference, line_number) = row
            if upper != at:
                for number in range(at + 1, upper + 1):
                    if truncate and number != upper:
                        continue
                    one = by_number.get(number)
                    if one is None:
                        continue
                    # A star marks a level that *has* an observed energy;
                    # `conv_out.pl` carries the opposite flag, `exp_c`,
                    # which is the star RCE prints for one that has not.
                    out.write('$%4d    J=%4.1f  %11s%s%2s \n'
                              % (number, one.lev.j, one.label,
                                 fortfmt.write_f_wide(
                                     abs(one.lev.e_obs) * 1000, 13, 3),
                                 '*' if one.lev.observed else ' '))
                at = upper
            out.write('+%4d%5d%s%2s%s%6d%s%s%6d\n'
                      % (lower, intensity,
                         fortfmt.write_f_wide(low.e_obs * 1000, 12, 3),
                         '*' if low.observed else ' ',
                         fortfmt.write_f_wide(computed, 13, 3),
                         observed_intensity,
                         fortfmt.write_f_wide(wavenumber, 12, 3),
                         fortfmt.write_f_wide(difference, 11, 3),
                         line_number))
