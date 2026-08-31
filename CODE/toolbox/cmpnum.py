"""Tolerance-based numeric comparison of two Fortran output files.

CLAUDE.md 1 makes byte-identical output the acceptance test, and `diff`
is the tool for that.  The constants work (`PLAN_fundconst.md`)
deliberately changes numbers, so from its Stage B onward `diff` says
only "different" and stops being an oracle.  This module is the
replacement: it separates a text file into the part that must not change
at all -- the text skeleton -- and the part that is allowed to move
within a tolerance.

  * The **skeleton** is everything that is not a number, and it is
    required to match exactly.  A changed heading, a lost line, a
    configuration label that moved, a field that overflowed to `*******`
    -- all of those are structural, and structural change is never
    acceptable in this package.

  * The **numeric vector** is then compared pairwise with a relative and
    an absolute tolerance, and reported as a distribution rather than as
    a pass/fail: max deviation with its line and column, rms, a decade
    histogram, and the worst offenders.  A constants change has a
    predictable signature -- everything at 1e-5, say -- and the histogram
    is what shows whether the diff really has that shape or whether
    something else came along with it.

Four things about the tokenizer are decisions rather than details:

  * **Digits glued to letters stay in the skeleton.**  `Sn7+`, `d65p`
    and the `d7` of `4d7` are labels, not data.  A number is recognized
    only where it does not follow a letter, so a configuration label
    cannot be silently reinterpreted as a number and compared as one.

  * **An exponent is at most three digits**, which is all the E and D
    edit descriptors ever write.  Without that limit the true-jj
    basis-state label `232E1518` reads as 232e1518 = Inf, and a pair of
    them makes a NaN.

  * **Only `E` is an exponent letter, not `D`.**  No output file this
    package writes contains a `D` exponent -- checked across OUTG11,
    OUT36, OUT2, OUTGINE and TAPE2E -- while `4d7` and its kind are
    everywhere.  Reading `D` as an exponent would turn that
    configuration label into 4e7.  Pass `d_exponents=True` if a future
    FORMAT ever needs it.

  * **A bare-sign exponent is two tokens, not one.**  Fortran writes
    `1.234-105` when a three-digit exponent will not fit the field, but
    two packed fields whose second is negative look exactly the same and
    are far more common here.  The ambiguity is resolved towards the
    common case; `-105` is read as a separate number.

  * **Runs of blanks collapse to one in the skeleton.**  These are
    fixed-width fields, so a number that gains a digit eats one of the
    blanks padding it.  Comparing whitespace literally would report every
    such number twice, once as a numeric change and once as a structural
    one.

Packed fields that run together with no separator at all
(`1000.0002000.000`) cannot be split without knowing the field widths,
and are tokenized as one long number.  Both files tokenize the same way,
so the comparison stays valid; the deviation reported for such a token
is just not physically meaningful.  `long_tokens` in the report counts
them so they are visible rather than silent.
"""

import hashlib
import math
import os
import re

#: Lines whose numbers are wall-clock time and never compare equal.
#: CLAUDE.md 1 already excludes exactly these three from the OUTG11
#: acceptance test, and they are excluded here for the same reason; the
#: two `t(...)=` forms are RCN's spelling of the same thing in OUT36.
DEFAULT_IGNORES = (
    r'TIME=',
    r'FINISHED .* AT',
    r'DIAGTIME',
    r't\(run\)=',
    r't\(job\)=',
)

# A label is matched *first* and stepped over, so that no part of one
# can be read as a number.  A lookbehind cannot do this job: rejecting a
# digit that follows a letter leaves the *next* digit of `d65p` looking
# perfectly innocent, and `5` comes back as a number.  And a lookbehind
# that also rejected a digit after a digit would quietly drop the sign of
# the second of two packed fields -- `418.430-238.989` would compare as
# +238.989 -- because that minus follows a digit.  Consuming the whole
# label in one bite avoids both.
_NUM = r'[+-]?(?:\d+(?:\.\d*)?|\.\d+)'
_LABEL = r'[A-Za-z_][A-Za-z0-9_]*'

# An exponent is at most three digits.  Fortran's E and D edit descriptors
# write E+dd or E+ddd and nothing wider -- when a three-digit exponent
# will not fit the field they drop the letter instead, which is the
# bare-sign form dealt with above.  So a letter followed by four or more
# digits is not an exponent, and the guard has to say so: the true-jj
# basis-state labels in `OUTG11` are things like `232E1518`, and read as
# a number that is 232e1518, i.e. +Inf.  Two Infs then give Inf-Inf = NaN
# in `deviation`, and a Pr III `OUTG11` comparison died on it.  With the
# width limited, `232E1518` tokenizes exactly as the `230D1518` beside it
# always has: the leading digits are a number and the rest is a label,
# which is the same convention as `4d7`.
_EXP_E = r'[Ee][+-]?\d{1,3}(?![A-Za-z0-9_])'
_EXP_ED = r'[EeDd][+-]?\d{1,3}(?![A-Za-z0-9_])'

_NUMBER_E = re.compile(r'(?P<num>' + _NUM + r'(?:' + _EXP_E + r')?)'
                       r'|(?P<label>' + _LABEL + r')')
_NUMBER_ED = re.compile(r'(?P<num>' + _NUM + r'(?:' + _EXP_ED + r')?)'
                        r'|(?P<label>' + _LABEL + r')')

_BLANKS = re.compile(r' +')

#: A token longer than this is almost certainly two packed fields read
#: as one.  The widest single field the FORMATs here write is around
#: fifteen characters.
LONG_TOKEN = 18


def looks_binary(path):
    """True when `path` is Fortran unformatted rather than text.

    A NUL byte in the first block is the test.  TAPE2N is written as one
    unformatted record per configuration, mixing integers, characters and
    doubles inside a single record; comparing it numerically needs a
    reader that knows that record's layout, which is not something a
    generic tokenizer can infer.  `toolbox/tape2n.py` is that reader, and
    `compare_tape2n` below uses it; every other unformatted file still
    gets the md5 answer.
    """
    with open(path, 'rb') as handle:
        block = handle.read(8192)
    return b'\0' in block


def md5(path):
    """The md5 of a file, read in blocks: OUTG11 is 122 MB for Pr III."""
    digest = hashlib.md5()
    with open(path, 'rb') as handle:
        for block in iter(lambda: handle.read(1 << 20), b''):
            digest.update(block)
    return digest.hexdigest()


def tokenize(line, d_exponents=False):
    """Split one line into (skeleton, values, columns, texts).

    `skeleton` is the tuple of text pieces between the numbers, with runs
    of blanks collapsed; `values` are the numbers in order; `columns` are
    their 1-based starting columns in the original line, which is what a
    report has to quote for the number to be findable in the file; and
    `texts` are the tokens as written.
    """
    pattern = _NUMBER_ED if d_exponents else _NUMBER_E
    pieces = []
    values = []
    columns = []
    texts = []
    at = 0
    for hit in pattern.finditer(line):
        if hit.lastgroup == 'label':
            continue            # a label: it stays in the skeleton piece
        text = hit.group(0)
        try:
            value = float(text.replace('D', 'E').replace('d', 'e'))
        except ValueError:                      # a lone sign, or a lone `.`
            continue
        pieces.append(_BLANKS.sub(' ', line[at:hit.start()]))
        values.append(value)
        columns.append(hit.start() + 1)
        texts.append(text)
        at = hit.end()
    pieces.append(_BLANKS.sub(' ', line[at:]))
    return tuple(pieces), values, columns, texts


def deviation(a, b):
    """Absolute and relative deviation of one pair.

    The relative measure is symmetric -- divided by the larger magnitude
    -- so that neither file is privileged as "the truth".  Stage B has no
    truth; it has two runs.
    """
    absolute = abs(a - b)
    scale = max(abs(a), abs(b))
    return absolute, (absolute / scale if scale else 0.0)


class Report(object):
    """What one file-to-file comparison found."""

    def __init__(self, name_a, name_b):
        self.name_a = name_a
        self.name_b = name_b
        self.lines = 0
        self.lines_ignored = 0
        self.compared = 0
        self.long_tokens = 0
        self.structural = []        # (line number, what changed)
        self.offenders = []         # (rel, abs, line, column, a, b, ta, tb)
        self.largest = 0.0
        self.sum_sq = 0.0
        self.histogram = {}
        self.binary = False         # Fortran unformatted
        self.binary_same = None     # md5 equality, always set for those
        self.layout = None          # name of the record layout, if one applied
        self.records = 0

    @property
    def ok(self):
        """True when nothing structural moved and nothing exceeded tolerance."""
        if self.binary and not self.layout:
            return bool(self.binary_same)
        return not self.structural and not self.offenders

    @property
    def rms(self):
        return math.sqrt(self.sum_sq / self.compared) if self.compared else 0.0

    def worst(self, count=10):
        return sorted(self.offenders, reverse=True)[:count]

    def record(self, relative):
        self.sum_sq += relative * relative
        if relative > self.largest:
            self.largest = relative
        key = None if relative == 0.0 else int(math.floor(math.log10(relative)))
        self.histogram[key] = self.histogram.get(key, 0) + 1


def compare_lines(lines_a, lines_b, report, rtol=0.0, atol=0.0,
                  ignores=DEFAULT_IGNORES, d_exponents=False):
    """Compare two already-read lists of lines into `report`.

    Every numeric pair goes into the histogram, whether or not it is
    within tolerance -- the distribution is the point of the exercise, and
    a histogram of only the failures would hide the shape of the change.
    """
    skip = [re.compile(p) for p in ignores]
    if len(lines_a) != len(lines_b):
        report.structural.append(
            (0, '%d lines in %s, %d in %s'
                % (len(lines_a), report.name_a, len(lines_b), report.name_b)))
        return report
    for number, (line_a, line_b) in enumerate(zip(lines_a, lines_b), 1):
        report.lines += 1
        if any(p.search(line_a) or p.search(line_b) for p in skip):
            report.lines_ignored += 1
            continue
        skel_a, values_a, columns, texts_a = tokenize(line_a, d_exponents)
        skel_b, values_b, _, texts_b = tokenize(line_b, d_exponents)
        if skel_a != skel_b:
            report.structural.append((number, _first_difference(skel_a, skel_b)))
            continue
        for i, (a, b) in enumerate(zip(values_a, values_b)):
            report.compared += 1
            if max(len(texts_a[i]), len(texts_b[i])) > LONG_TOKEN:
                report.long_tokens += 1
            absolute, relative = deviation(a, b)
            report.record(relative)
            if absolute > atol and relative > rtol:
                report.offenders.append((relative, absolute, number,
                                         columns[i], a, b,
                                         texts_a[i], texts_b[i]))
    return report


def _first_difference(skel_a, skel_b):
    """Say where two skeletons part company, in words a reader can act on."""
    if len(skel_a) != len(skel_b):
        return ('%d numbers on the line, %d on the other'
                % (len(skel_a) - 1, len(skel_b) - 1))
    for piece_a, piece_b in zip(skel_a, skel_b):
        if piece_a != piece_b:
            return 'text %r became %r' % (piece_a, piece_b)
    return 'skeletons differ'                                   # unreachable


def read_lines(path):
    """Read a text file for comparison.

    latin-1 because Fortran output is bytes and every byte must round-trip
    rather than raise; trailing blanks go because the FORMAT pads to the
    field width and that padding is not information.
    """
    with open(path, 'r', encoding='latin-1', newline=None) as handle:
        return [line.rstrip() for line in handle]


def compare_files(path_a, path_b, rtol=0.0, atol=0.0, ignores=DEFAULT_IGNORES,
                  d_exponents=False):
    """Compare two files by name.

    A Fortran unformatted file is not tokenized -- see `looks_binary` --
    but it is not skipped either.  It gets an md5 answer, and `binary` on
    the report says so, because a silent skip that still reports success
    is exactly the failure mode this whole exercise is guarding against.
    """
    report = Report(str(path_a), str(path_b))
    if looks_binary(path_a) or looks_binary(path_b):
        report.binary = True
        report.binary_same = md5(path_a) == md5(path_b)
        if _is_tape2n(path_a) and _is_tape2n(path_b):
            compare_tape2n(path_a, path_b, report, rtol=rtol, atol=atol)
        return report
    return compare_lines(read_lines(path_a), read_lines(path_b), report,
                         rtol=rtol, atol=atol, ignores=ignores,
                         d_exponents=d_exponents)


def _is_tape2n(path):
    """True when `path` is named like RCN's TAPE2N.

    The layout in `toolbox.tape2n` belongs to that one file and to no
    other, so it is applied by name rather than guessed at.  A file that
    is named TAPE2N and does not parse is reported as such -- see
    `compare_tape2n` -- and never quietly treated as something else.
    """
    return os.path.basename(str(path)).upper().startswith('TAPE2N')


def compare_tape2n(path_a, path_b, report, rtol=0.0, atol=0.0):
    """Compare two TAPE2N files record by record, into `report`.

    The division of labour is the one the text comparison already makes.
    The integers and the character fields of a record are structure: a
    difference in one of them means the run described a different atom,
    a different number of orbitals or a different mesh, and no tolerance
    applies to that.  Only the `real*8` fields are compared numerically.

    A file that will not parse is not a failure of the comparison but a
    failure of the reader to still match `RCN36K.F`, and it is reported
    in those words rather than falling back to md5 and calling the file
    "not quantifiable".
    """
    from toolbox import tape2n

    try:
        recs_a = tape2n.read_records(path_a)
        recs_b = tape2n.read_records(path_b)
    except tape2n.Tape2nError as problem:
        report.layout = 'TAPE2N'
        report.structural.append(
            (0, 'TAPE2N does not parse: %s.  toolbox/tape2n.py has drifted '
                'from the write statement in RCN36K.F' % problem))
        return report

    report.layout = 'TAPE2N'
    if len(recs_a) != len(recs_b):
        report.structural.append(
            (0, '%d records in %s, %d in %s'
                % (len(recs_a), report.name_a, len(recs_b), report.name_b)))
        return report

    for number, (rec_a, rec_b) in enumerate(zip(recs_a, recs_b), 1):
        report.records += 1
        struct_a = tape2n.structure(rec_a)
        struct_b = tape2n.structure(rec_b)
        if struct_a != struct_b:
            report.structural.append(
                (number, _first_field_difference(struct_a, struct_b)))
            continue
        for (label, a), (_, b) in zip(tape2n.flatten(rec_a),
                                      tape2n.flatten(rec_b)):
            report.compared += 1
            absolute, relative = deviation(a, b)
            report.record(relative)
            if absolute > atol and relative > rtol:
                report.offenders.append((relative, absolute, number, 0,
                                         a, b, label, ''))
    return report


def _first_field_difference(struct_a, struct_b):
    """Name the first integer or character field of a record that moved."""
    for (name, value_a), (_, value_b) in zip(struct_a, struct_b):
        if value_a != value_b:
            return 'field %s: %r became %r' % (name, value_a, value_b)
    return 'records differ in a field this reader does not name'


def format_report(report, top=10, rtol=0.0, atol=0.0):
    """The printed form: one block per comparison.

    No colour and no terminal-width assumptions -- this goes into a log
    at least as often as onto a screen.
    """
    out = ['%s' % report.name_a, '  vs %s' % report.name_b]
    if report.binary and not report.layout:
        out.append('  Fortran unformatted -- md5 %s'
                   % ('identical' if report.binary_same else
                      'DIFFERS.  How far the numbers moved cannot be said '
                      'without a reader for this file\'s record layout'))
        return '\n'.join(out)
    if report.layout:
        out.append('  Fortran unformatted, read with the %s record layout '
                   '-- md5 %s'
                   % (report.layout,
                      'identical' if report.binary_same else 'differs'))
    if report.structural:
        out.append('  STRUCTURAL MISMATCH -- text around the numbers changed, '
                   'which is never acceptable:')
        for number, reason in report.structural[:top]:
            where = 'line %d' % number if number else 'file'
            out.append('    %-12s %s' % (where, reason))
        if len(report.structural) > top:
            out.append('    ... and %d more' % (len(report.structural) - top))
    if report.layout:
        out.append('  records          %d compared' % report.records)
    else:
        out.append('  lines            %d compared, %d ignored as timing'
                   % (report.lines - report.lines_ignored,
                      report.lines_ignored))
    out.append('  numbers compared %d' % report.compared)
    if report.long_tokens:
        out.append('  long tokens      %d -- packed fields read as one number; '
                   'their deviations are not meaningful' % report.long_tokens)
    if report.compared:
        out.append('  max relative     %.3e' % report.largest)
        out.append('  rms relative     %.3e' % report.rms)
        out.append('  outside rtol=%g atol=%g: %d'
                   % (rtol, atol, len(report.offenders)))
        out.append('  relative deviation by decade:')
        for label, count in decades(report):
            out.append('    %-16s %d' % (label, count))
        worst = report.worst(top)
        if worst:
            out.append('  worst %d:' % len(worst))
            for rel, absolute, number, column, a, b, ta, tb in worst:
                if report.layout:
                    out.append('    rec %-5d %-20s %-17.9g %-17.9g '
                               'rel %.3e  abs %.3e'
                               % (number, ta, a, b, rel, absolute))
                else:
                    out.append('    line %-7d col %-4d %-16s %-16s '
                               'rel %.3e  abs %.3e'
                               % (number, column, ta, tb, rel, absolute))
            if len(report.offenders) > len(worst):
                out.append('    ... and %d more'
                           % (len(report.offenders) - len(worst)))
    return '\n'.join(out)


def decades(report):
    """Histogram rows: exact agreement first, then one row per decade."""
    rows = []
    if None in report.histogram:
        rows.append(('exact', report.histogram[None]))
    for key in sorted(k for k in report.histogram if k is not None):
        rows.append(('1e%d..1e%d' % (key, key + 1), report.histogram[key]))
    return rows
