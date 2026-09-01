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

import bisect
import math
import shutil
from dataclasses import dataclass, field

from . import files, fortfmt, outg11, rce

#: A line assigned to an observation keeps only the transitions within
#: this fraction of the strongest one assigned to it: a much weaker
#: candidate is not a competing identification (`conv_out.pl:$Athresh`).
THRESHOLD = 0.1

#: IDEN2's level number is I4 and its arrays are sized for this many.
MAX_LEVELS = 9999

#: A `dlv.dat` whose intensities are already logarithmic is now the
#: common case, and `convert2log` on one of those takes the log a second
#: time -- silently, since both scales are just numbers in an I5 field.
#: The shapes are far apart, though: line intensities are steeply
#: skewed, so a linear list of any size sits decades below its own
#: maximum at the median, and 10*ln compresses both the skew and the
#: range.  Neither number is proof on its own, so both have to say log,
#: and the result only warns -- reinterpreting the author's own data
#: would be the worse failure.
LOG_MEDIAN = 0.2                #: median/max above this looks logarithmic
LOG_RANGE = 50.0                #: max/min below this looks logarithmic

#: Below this many lines the distribution says nothing either way.
LOG_SAMPLE = 20


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


def looks_logarithmic(lines):
    """Whether the observed intensities look as though they are already logs.

    See `LOG_MEDIAN`.  A short list is left alone: the shape of a
    distribution of a dozen numbers says nothing.
    """
    values = sorted(one.intensity for one in lines if one.intensity > 0)
    if len(values) < LOG_SAMPLE:
        return False
    biggest, smallest = values[-1], values[0]
    median = values[len(values) // 2]
    return (median > LOG_MEDIAN * biggest
            and biggest < LOG_RANGE * smallest)


def write_dlv(lines, path='dlv.dat'):
    """Rewrite the observed list with its intensities on a log scale.

    IDEN2 compares intensities directly, so a list whose intensities run
    over several decades has to be put on the same footing as the
    calculated ones, which are logarithmic.  The original is kept beside
    it as `dlv.dat.bak`.
    """
    biggest = max(one.intensity for one in lines)
    if biggest <= 0:
        raise IdenError('every intensity in %s is zero; there is nothing to '
                        'put on a log scale' % path)
    if 999 - 10 * math.log(biggest) < 0:
        raise IdenError('the largest intensity in %s, %g, is too large to put '
                        'on a log scale' % (path, biggest))
    shutil.copyfile(path, path + '.bak')
    with open(path, 'w') as out:
        for one in lines:
            # An intensity below 1 has a negative log, and `read_dlv`
            # takes the field with abs(), which would read it back as a
            # strong line.  Cap at zero instead, which is IDEN2's "no
            # measured intensity" -- the truthful thing to say about a
            # line that weak.  Choosing a scale on which the observed
            # intensities clear 1 is the user's job.
            one.intensity = (max(0.0, 10 * math.log(one.intensity))
                             if one.intensity > 0 else 0.0)
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
        # An OUTG11 written before RCG printed the sections of a
        # one-state block has none for it but LS; the level is that
        # one state, so its label is the one row LS printed.
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
    # `write_dlv` replaces these with their logs, and the scale the
    # calculated intensities are put on is set by the observed maximum
    # *as the file gives it* -- `conv_out.pl` keeps it in $I_obs_max,
    # which WriteDLV does not touch.  So take it before the rewrite.
    biggest = max([one.intensity for one in lines], default=0.0)
    if options.convert2log:
        if looks_logarithmic(lines):
            print('   the intensities in dlv.dat already look logarithmic; '
                  'convert2log will take their log a second time')
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

    rows = _trans_rows(printed, by_level, assigned, lines, biggest, options)
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
        # A difference of two observed energies, matched against an
        # observed line list: safe inside one relative-energy group,
        # where the common offset cancels, and not across two.
        if not rce.comparable(one.upper, one.lower):
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


def _trans_rows(transitions, by_level, assigned, lines, biggest, options):
    """One row per line: the two IDEN2 level numbers and the numbers IDEN2 wants.

    `biggest` is the largest observed intensity **as dlv.dat gave it**,
    before `write_dlv` may have replaced the list with its logs.
    """
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
    # The ground level stays whatever the lines say.  It is the zero of
    # every energy in the analysis and the one level IDEN2 must be able
    # to reach; a run whose printed lines happen not to touch it -- a
    # wavelength window that stops short of the lowest opposite-parity
    # level, say -- would otherwise hand IDEN2 a spectrum with no
    # bottom.  `numbered` runs from the highest energy down, so this is
    # its last entry, and keeping it changes no other level's number.
    if numbered:
        wanted.add(numbered[-1].number)
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
            # IDEN2 fits **absolute** energies and the weight column is
            # the lever, so a level whose offset is unknown is not one
            # of the well-determined ones however well it is observed.
            observed = lev.absolute
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


# --------------------------------------------------------------------------
# reading the two files back
# --------------------------------------------------------------------------
#
# What `_write_enlev` and `_write_trans` put out is what these read in,
# and the column spans below are the single statement of where each field
# is: the writers above build their records with `%` formats, the readers
# here take them apart by name, and `WORK/Pr3/IDEN2_file_formats.md` is
# what both agree with.
#
# IDEN2 writes the files back after a session, so a reader has to survive
# a record the user's own run has been through.  Only version 2 of
# `ENLEV.DAT` -- the four-column level number -- is read: version 1 is the
# legacy layout, IDEN2 keeps whichever it was given, and nothing in this
# package writes version 1.

ENLEV = 'ENLEV.DAT'
TRANS = 'TRANS.DAT'

#: The parity the sixth character of a label stands for.  `_label` puts
#: it there from the configuration set's own parity, so it is the level's
#: parity and not a set number, whatever the set numbering happens to be.
PARITY_MARK = {'_': 'e', '~': 'o'}

#: `ENLEV.DAT` version 2, 0-based [start, end) column spans.
ENLEV_COLUMNS = {
    'number': (0, 4), 'e_calc': (4, 16), 'uncertainty': (16, 26),
    'e_fixed': (26, 38), 'fixed': (39, 40), 'difference': (40, 52),
    'j': (52, 57), 'label': (59, 70),
}

#: `TRANS.DAT`, the `$` record that opens an upper level's group.
TRANS_LEVEL_COLUMNS = {
    'number': (1, 5), 'j': (11, 15), 'label': (17, 28),
    'energy': (28, 41), 'fixed': (42, 43),
}

#: `TRANS.DAT`, the `+` record of one transition down from that level.
TRANS_LINE_COLUMNS = {
    'lower': (1, 5), 'intensity': (5, 10), 'e_lower': (10, 22),
    'fixed': (23, 24), 'energy': (24, 37), 'obs_intensity': (37, 43),
    'wavenumber': (43, 55), 'difference': (55, 66), 'line': (66, 72),
}

#: Both records are 72 columns wide; a version 1 ENLEV.DAT record is 71.
CARD = 72


def _span(card, columns, name):
    start, end = columns[name]
    return card[start:end]


def _put(card, columns, name, text):
    """Write `text` into one field, leaving every other column alone.

    A rewritten record must differ from the one it came from only in the
    fields that were changed: IDEN2 keeps columns of its own in these
    files -- the ones the format document marks as unused here -- and a
    record rebuilt from a `%` format would lose them.
    """
    start, end = columns[name]
    width = end - start
    text = str(text).rjust(width)
    if len(text) > width:
        raise IdenError('%r does not fit the %d columns of the %s field'
                        % (text, width, name))
    card = card.ljust(end)
    return card[:start] + text + card[end:]


def _number(card, columns, name, path, at, default=0.0):
    text = _span(card, columns, name).strip()
    if not text:
        return default
    try:
        return float(text)
    except ValueError:
        raise IdenError('%s line %d: %r is not a number in the %s field'
                        % (path, at, text, name))


def _integer(card, columns, name, path, at):
    text = _span(card, columns, name).strip()
    if not text:
        return 0
    try:
        return int(text)
    except ValueError:
        raise IdenError('%s line %d: %r is not a number in the %s field'
                        % (path, at, text, name))


def _j2(text, path, at):
    """2J as an integer, so that no two Js are ever compared as floats."""
    value = text.strip()
    try:
        doubled = int(round(2 * float(value)))
    except ValueError:
        raise IdenError('%s line %d: %r is not a J' % (path, at, value))
    if doubled < 0:
        raise IdenError('%s line %d: J is negative' % (path, at))
    return doubled


@dataclass
class EnlevLevel:
    """One record of `ENLEV.DAT`, and the columns it was read from."""

    number: int
    e_calc: float
    uncertainty: float
    e_fixed: float
    fixed: bool                 # the star: an experimental level
    j2: int                     # 2J
    label: str
    card: str                   # the record, with any new field written in
    at: int = 0                 # the line it came from
    width: int = 0              # the columns it filled in the file

    @property
    def parity(self):
        """`'e'`, `'o'`, or `''` if the label carries no parity mark."""
        return PARITY_MARK.get(self.label[5:6], '')

    def update(self, e_fixed=None, uncertainty=None):
        """Put a new fixed energy, a new uncertainty, or both, in place.

        `Ef - Ec` goes with the energy.  IDEN2 does not read that column
        and recomputes it, but a file left saying one thing in two places
        is a file that misleads whoever opens it in an editor.
        """
        if e_fixed is not None:
            self.e_fixed = e_fixed
            self.card = _put(self.card, ENLEV_COLUMNS, 'e_fixed',
                             fortfmt.write_f_wide(e_fixed, 12, 3))
            self.card = _put(self.card, ENLEV_COLUMNS, 'difference',
                             fortfmt.write_f_wide(e_fixed - self.e_calc,
                                                  12, 3))
        if uncertainty is not None:
            self.uncertainty = uncertainty
            self.card = _put(self.card, ENLEV_COLUMNS, 'uncertainty',
                             fortfmt.write_f_wide(uncertainty, 10, 3))


def enlev_version(cards, path=ENLEV):
    """1 or 2, from where the slashes around the label fall.

    The two layouts differ by one column throughout, so the level number
    alone cannot tell them apart -- but the label is fenced by slashes
    that no other field holds, and those are one column apart too.
    """
    for at, card in enumerate(cards, 1):
        if not card.strip():
            continue
        if card[58:59] == '/':
            return 2
        if card[57:58] == '/':
            return 1
        raise IdenError('%s line %d is not an ENLEV.DAT record: neither '
                        'column 58 nor column 59 holds the slash that opens '
                        'the label' % (path, at))
    raise IdenError('%s holds no levels' % path)


def read_enlev(path=ENLEV):
    """The levels of `ENLEV.DAT`, in the order the file has them."""
    try:
        cards = files.read_lines(path, strip=False)
    except OSError:
        raise IdenError('%s is not here' % path)
    version = enlev_version([card.ljust(CARD) for card in cards], path)
    if version != 2:
        raise IdenError(
            '%s is in the legacy version %d layout, which this does not read; '
            'open it in IDEN2 and save it, or rebuild it with conv_out, for '
            'the version 2 layout' % (path, version))

    out = []
    for at, card in enumerate(cards, 1):
        if not card.strip():
            continue
        width, card = len(card), card.ljust(CARD)
        text = _span(card, ENLEV_COLUMNS, 'number').strip()
        try:
            number = int(text)
        except ValueError:
            raise IdenError('%s line %d: %r is not a level number'
                            % (path, at, text))
        out.append(EnlevLevel(
            number=number,
            e_calc=_number(card, ENLEV_COLUMNS, 'e_calc', path, at),
            uncertainty=_number(card, ENLEV_COLUMNS, 'uncertainty', path, at),
            e_fixed=_number(card, ENLEV_COLUMNS, 'e_fixed', path, at),
            fixed=_span(card, ENLEV_COLUMNS, 'fixed') == '*',
            j2=_j2(_span(card, ENLEV_COLUMNS, 'j'), path, at),
            label=_span(card, ENLEV_COLUMNS, 'label'),
            card=card, at=at, width=width))
    if not out:
        raise IdenError('%s holds no levels' % path)
    return out


def _trim(one):
    """A record back as text, no narrower than it was in the file.

    Fortran pads a record to its FORMAT's full width, trailing blanks and
    all, and IDEN2 writes the files back the same way; a record this has
    not written into has to come back out as the bytes it went in as, so
    the blanks the file had are kept and only what a rewrite pushed past
    them is dropped.
    """
    text = one.card.rstrip()
    return text.ljust(one.width) if len(text) < one.width else text


def write_enlev(levels, path=ENLEV, backup=None):
    """Write the levels back, every record keeping its own columns."""
    return files.write_deck(path, [_trim(one) for one in levels],
                            files.line_ending(path), backup=backup)


def fixed_targets(levels, flip=False):
    """The experimental levels of `ENLEV.DAT` as `levids` match targets.

    Energies go over in the unit the matcher works in, 10^3 cm-1, and the
    *fixed* energy is the one to match on: it is the observed value,
    which is what an experimental level list holds.  `flip` reverses the
    parity marks and is there only so that a caller can measure how much
    worse the wrong reading of them would be.
    """
    from . import levids

    marks = PARITY_MARK
    if flip:
        marks = dict((mark, 'o' if parity == 'e' else 'e')
                     for mark, parity in PARITY_MARK.items())
    return [levids.TargetLevel(one.number, one.e_fixed * 0.001, one.j2,
                               marks.get(one.label[5:6], ''), False)
            for one in levels if one.fixed]


@dataclass
class TransRecord:
    """One record of `TRANS.DAT`: a `$` level or a `+` transition."""

    kind: str                   # '$', '+', or '' for anything else
    upper: int                  # the group's upper level, on both kinds
    lower: int                  # the '+' record's lower level, else 0
    line: int                   # the assigned observation, 0 if none
    card: str
    at: int = 0
    width: int = 0

    def set_level_energy(self, energy, fixed):
        """The `$` record's energy and star.  `_write_trans` writes |E|."""
        self.card = _put(self.card, TRANS_LEVEL_COLUMNS, 'energy',
                         fortfmt.write_f_wide(abs(energy), 13, 3))
        self.card = _put(self.card, TRANS_LEVEL_COLUMNS, 'fixed',
                         '*' if fixed else ' ')

    def set_lower(self, energy, fixed):
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'e_lower',
                         fortfmt.write_f_wide(energy, 12, 3))
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'fixed',
                         '*' if fixed else ' ')

    def set_transition(self, energy):
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'energy',
                         fortfmt.write_f_wide(energy, 13, 3))

    def refresh(self, computed):
        """Recompute `dwn` against a transition energy that has moved.

        The assignment itself is left as it is; this is for the record
        whose two levels have new energies but whose observed line has
        not changed, where the old difference would now be wrong.
        """
        wavenumber = 0.0
        if self.line:
            text = _span(self.card, TRANS_LINE_COLUMNS, 'wavenumber').strip()
            try:
                wavenumber = float(text or 0)
            except ValueError:
                return
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'difference',
                         fortfmt.write_f_wide(wavenumber - computed
                                              if self.line else 0.0, 11, 3))

    def assign(self, line, observation=None, computed=None):
        """Tie this transition to an observation, or to none at all.

        `line` is the row of `dlv.dat`, and it is the only one of these
        four fields IDEN2 reads; the intensity, the wavenumber and `dwn`
        beside it are written for the reader's sake, and are zeroed with
        it so that an unassigned record does not go on quoting the
        observation it has just lost.
        """
        self.line = line
        keep = observation if line else None
        intensity = keep.intensity if keep else 0.0
        wavenumber = keep.wavenumber if keep else 0.0
        difference = (wavenumber - computed
                      if keep and computed is not None else 0.0)
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'obs_intensity',
                         '%6d' % int(round(intensity)))
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'wavenumber',
                         fortfmt.write_f_wide(wavenumber, 12, 3))
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'difference',
                         fortfmt.write_f_wide(difference, 11, 3))
        self.card = _put(self.card, TRANS_LINE_COLUMNS, 'line', '%6d' % line)


def read_trans(path=TRANS):
    """The records of `TRANS.DAT`, in file order, `$` groups and all.

    Every record is kept, one this does not recognise included, so that
    the file can be written back with nothing dropped.
    """
    try:
        cards = files.read_lines(path, strip=False)
    except OSError:
        raise IdenError('%s is not here' % path)

    out, upper = [], 0
    for at, card in enumerate(cards, 1):
        if not card.strip():
            continue
        width, card = len(card), card.ljust(CARD)
        if card[0] == '$':
            upper = _integer(card, TRANS_LEVEL_COLUMNS, 'number', path, at)
            out.append(TransRecord('$', upper, 0, 0, card, at, width))
        elif card[0] == '+':
            if not upper:
                raise IdenError('%s line %d: a transition record before any '
                                'level record' % (path, at))
            out.append(TransRecord(
                '+', upper,
                _integer(card, TRANS_LINE_COLUMNS, 'lower', path, at),
                _integer(card, TRANS_LINE_COLUMNS, 'line', path, at),
                card, at, width))
        else:
            out.append(TransRecord('', upper, 0, 0, card, at, width))
    if not out:
        raise IdenError('%s holds no records' % path)
    return out


def write_trans(records, path=TRANS, backup=None):
    return files.write_deck(path, [_trim(one) for one in records],
                            files.line_ending(path), backup=backup)


def dlv_nearest(lines):
    """A lookup from a wavenumber to the nearest row of `dlv.dat`.

    Gives back `(row, apart)` -- the row number, counting from 1, and how
    far the row's wavenumber is from the one asked for -- or `(0, None)`
    for an empty list.  The list runs downward in wavenumber, so it is
    sorted once here and then bisected.
    """
    order = sorted(range(len(lines)), key=lambda i: lines[i].wavenumber)
    values = [lines[i].wavenumber for i in order]

    def nearest(wavenumber):
        at = bisect.bisect_left(values, wavenumber)
        best, chosen = None, 0
        for i in (at - 1, at):
            if 0 <= i < len(values):
                apart = abs(values[i] - wavenumber)
                if best is None or apart < best:
                    best, chosen = apart, order[i] + 1
        return chosen, best

    return nearest


def dlv_index(lines, tolerance):
    """A lookup from an observed wavenumber to its row of `dlv.dat`.

    The row number *is* the identifier IDEN2 uses, so this is how a line
    named by its wavenumber somewhere else -- in the output of a level
    optimization, say -- is found again here.  0 comes back for a
    wavenumber no row of the list is within `tolerance` of.
    """
    nearest = dlv_nearest(lines)

    def find(wavenumber):
        at, apart = nearest(wavenumber)
        return 0 if apart is None or apart > tolerance else at

    return find
