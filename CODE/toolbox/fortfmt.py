"""Fortran fixed-field editing, the few descriptors these decks use.

Two things here are not what a naive str/float conversion would do, and
both change results:

  * On input, a blank field is zero, not an error (BLANK='NULL' is the
    default for formatted input), and an F w.d field with no decimal
    point in it has one implied d digits from the right.  IN36 relies on
    this: `0pf6.5` and `f70.3` fields are routinely punched without a
    point.

  * On output, RCE writes a Slater parameter as a scaled integer and
    falls back to a floating field only when that would overflow.  The
    thresholds are part of the file format, not a cosmetic choice.
"""

import math
import re

_EXP = re.compile(r'[EeDd]')


def read_i(field, default=0):
    """Fortran `I w` input: blanks are zero, embedded blanks ignored."""
    s = field.replace(' ', '')
    if not s or s in '+-':
        return default
    return int(s)


def read_f(field, d=0, default=0.0):
    """Fortran `F w.d` / `E w.d` input.

    Blanks are ignored, an all-blank field is zero, and if the mantissa
    carries no decimal point the last `d` of its digits are the
    fraction.  An exponent may be written E, D, or as a bare sign.
    """
    s = field.replace(' ', '')
    if not s or s in '+-':
        return default
    s = _EXP.sub('E', s)
    if 'E' not in s:
        # A sign after the first character starts a bare exponent: 1.5-8.
        for i in range(1, len(s)):
            if s[i] in '+-':
                s = s[:i] + 'E' + s[i:]
                break
    if 'E' in s:
        mant, _, expo = s.partition('E')
    else:
        mant, expo = s, ''
    if '.' not in mant and d:
        sign = ''
        if mant and mant[0] in '+-':
            sign, mant = mant[0], mant[1:]
        mant = mant.rjust(d + 1, '0')
        mant = sign + mant[:len(mant) - d] + '.' + mant[len(mant) - d:]
    if expo:
        if expo in '+-':
            expo = '0'
        return float(mant + 'E' + expo)
    return float(mant)


def write_i(value, w):
    """Fortran `I w` output: all asterisks when the field overflows."""
    s = '%*d' % (w, int(round(value)))
    return '*' * w if len(s) > w else s


def write_f(value, w, d):
    """Fortran `F w.d` output."""
    s = '%*.*f' % (w, d, value)
    return '*' * w if len(s) > w else s


def write_e(value, w, d):
    """Fortran `1P E w.d` output, e.g. 1PE9.2 -> ' 1.23E-04'."""
    s = '%*.*E' % (w, d, value)
    if len(s) > w:
        # gfortran drops the exponent's leading zero rather than overflow.
        s = re.sub(r'E([+-])0(\d\d)', r'E\1\2', '%.*E' % (d, value)).rjust(w)
    return '*' * w if len(s) > w else s


def write_f_wide(value, w, d):
    """`F w.d`, giving up decimals rather than the number.

    The rule the author uses in the Fortran where a value may outgrow
    its field: try `F w.d`, then `F w.d-1` and so on down to `F w.0`,
    and only when the integer part alone will not fit fall back to
    exponential with the widest mantissa the field allows.  A reader
    that parses the field as a number gets a number in every case,
    which a field of asterisks does not give it.
    """
    for less in range(d, -1, -1):
        s = '%*.*f' % (w, less, value)
        if len(s) <= w:
            return s
    for digits in range(w - 6, -1, -1):
        s = write_e(value, w, digits)
        if '*' not in s:
            return s
    return '*' * w


def write_param(value, first=False):
    """The 9-column energy-parameter field of ING11 / RCEOUT.

    Mirrors the encoder quoted in `CODE/s11.pl:202-221`: a parameter is
    normally punched as I9 holding 10000 x its value; a parameter that is
    zero to within 1e-3 (and is not the first, which is EAV) goes out as
    1PE9.2; and one too large for the scaled integer falls back to F9.1
    for EAV or F9.2 for the rest.
    """
    if -1e-3 <= value <= 1e-3 and not first:
        return write_e(value, 9, 2)
    if -9999.99 <= value <= 99999.99:
        return write_i(10000.0 * value, 9)
    return write_f(value, 9, 1 if first else 2)


def half_int(two_j):
    """Format a doubled angular momentum as ASD writes it: 0, 1/2, 11/2."""
    two_j = int(round(two_j))
    if two_j % 2 == 0:
        return str(two_j // 2)
    return '%d/2' % two_j


def parse_half_int(text):
    """Inverse of half_int(); returns the doubled value."""
    text = text.strip()
    if '/' in text:
        num, _, den = text.partition('/')
        return int(round(2.0 * int(num) / int(den)))
    return int(round(2.0 * float(text)))


def is_whole(x, tol=1e-9):
    """True if x is an integer to within tol (RCN's dmod(x,1.d0) test)."""
    return abs(x - math.floor(x + 0.5)) < tol
