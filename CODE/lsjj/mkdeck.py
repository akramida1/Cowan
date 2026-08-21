#!/usr/bin/env python3
"""Write LSJJ.DAT -- the data deck the Fortran side of the LS->jj recursion
reads.  This script IS the recovery procedure: LSJJ.DAT is a derived file
and is never edited or repaired, only rebuilt from here.

The supported way to run it is not to: make_cfp rebuilds the deck along with
the CFP decks, on the machine the package is installed on, and that is the
only route an installed user has.  By hand, from CODE/lsjj:

    python mkdeck.py                 # writes ../LSJJ.DAT, where RCG looks
    python mkdeck.py somewhere/LSJJ.DAT
    python mkdeck.py -t deck.txt     # the old text form, for inspection

It reads the CFP card deck and SENIOR through cowan_cfp.py, which takes them
from $LSJJ_CFP and $LSJJ_SENIOR when those are set -- which is how make_cfp
points it at the ING11 it has just written and at CODE/SENIOR.  Everything
else it needs is in this directory and is deterministic: the LS cfp come from
Cowan's own tables via cowan_cfp.py, the jj cfp from the recursion in
recur.py.  Neither depends on any earlier LSJJ.DAT, so a lost or corrupted
deck costs nothing but the minute this takes to run.  FOR/LSJJ/README.txt has
the derivations; FOR/LSJJ/validate.py and cmpfort.py check the result.

The default form is UNFORMATTED SEQUENTIAL, which is what RCG's LSJJIN reads:
native INTEGER*4 and REAL*8 in native byte order, each record wrapped in the
4-byte length markers gfortran (and every other common compiler on this
platform) uses.  That makes the file machine specific -- carried to a machine
of the other endianness it will be refused, not misread, and the answer is to
run this script there.  Nothing in it was ever meant to be read by a human:
it is a megabyte of twelve-digit cfp.

Record layout, identical in both forms so that the text version can be
diffed against the binary one field for field.  The first record is new in
version 1 and carries the magic and version LSJJIN checks:

    MAGIC  VERSION  NLREC  NJREC

    section 1, LS cfp, one record per (l, n):

        LQ  N  NT  NP
        (MULT(i), 2L(i), IA1(i), IA2(i),  i = 1..NT)
        (IPAR(k),                         k = 1..NP)
        ((CFP(i,k), k = 1..NP),           i = 1..NT)

    IA1, IA2 are the two characters of Cowan's alpha field, as ICHAR
    values, so that a term can be matched against RCG's MULT/LBCD/IALF
    exactly.  IPAR(k) is the index, in the term list of l^(n-1), of the
    k-th column of the cfp matrix.

    section 2, jj cfp, one record per (j, n):

        2J  N  NC  NPJ
        (2Jp(k), v_p(k),  k = 1..NPJ)
        (2J(i), v(i), (CFP(i,k), k = 1..NPJ),  i = 1..NC)

In the text form the four LS lines and the jj lines are wrapped over as many
physical lines as they need, which list-directed input does not care about;
in the binary form each is exactly one unformatted record.

The l values covered are those Cowan's cfp deck covers: p, d, f for every n,
and g, h for n <= 2.  The jj tables run to n = 2j+1 for j <= 7/2 and to
n = 2 for j = 9/2 and 11/2, which is what those l values need.
"""

import os
import struct
import sys

import cowan_cfp
from cowan_cfp import decks, _term_ls, LSYMB
import recur

_LS = 'SPDFGH'

# Keep in step with the PARAMETER line in SUBROUTINE LSJJIN.
MAGIC = 20260801
VERSION = 1


def ls_records():
    tab = decks()
    out = []
    for l in range(1, 6):
        lsym = _LS[l]
        n = 1
        prev = [(1, 'S', '')]
        while (lsym, n) in tab:
            d = tab[(lsym, n)]
            ipar = [prev.index(p) + 1 for p in d.parents]
            out.append((l, n, list(d.terms), ipar, d.cfp))
            prev = list(d.terms)
            n += 1
    return out


def jj_records():
    out = []
    for j2 in (1, 3, 5, 7, 9, 11):
        nmax = j2 + 1 if j2 <= 7 else 2
        for n in range(1, nmax + 1):
            s = recur.jj_cfp(j2, n)
            p = recur.jj_cfp(j2, n - 1)
            out.append((j2, n, s.levels, p.levels, s.cfp))
    return out


def _wrap(items, per, fh):
    for i in range(0, len(items), per):
        fh.write(' ' + ' '.join(items[i:i + per]) + '\n')


class _Text(object):
    """The original list-directed form.  Kept for inspection and diffing."""

    def __init__(self, fh):
        self.fh = fh

    def ints(self, vals):
        _wrap(['%d' % v for v in vals], 20, self.fh)

    def reals(self, vals):
        _wrap(['%18.12E' % v for v in vals], 4, self.fh)


class _Binary(object):
    """Unformatted sequential: each record is length, payload, length, with
    the length a native 4-byte integer.  One call here is one Fortran READ
    on the other side, so the two record streams must stay in step."""

    def __init__(self, fh):
        self.fh = fh

    def _rec(self, payload):
        n = struct.pack('=i', len(payload))
        self.fh.write(n)
        self.fh.write(payload)
        self.fh.write(n)

    def ints(self, vals):
        vals = list(vals)
        self._rec(struct.pack('=%di' % len(vals), *vals))

    def reals(self, vals):
        vals = list(vals)
        self._rec(struct.pack('=%dd' % len(vals), *vals))


def write(path, text=False):
    ls = ls_records()
    jj = jj_records()
    mode, kind = ('w', _Text) if text else ('wb', _Binary)
    with open(path, mode, newline='' if text else None) as fh:
        w = kind(fh)
        w.ints([MAGIC, VERSION, len(ls), len(jj)])
        for l, n, terms, ipar, cfp in ls:
            w.ints([l, n, len(terms), len(ipar)])
            fld = []
            for mult, lsy, alpha in terms:
                a = (alpha + '  ')[:2]
                fld += [mult, 2 * LSYMB.index(lsy), ord(a[0]), ord(a[1])]
            w.ints(fld)
            w.ints(ipar)
            # One record, NT*NP values, row-major -- LSJJIN reads the whole
            # cfp matrix of l^n in a single READ.  List-directed input did
            # not care how this was split over lines; an unformatted record
            # boundary is not so forgiving.
            w.reals([x for row in cfp for x in row])
        for j2, n, lev, plev, cfp in jj:
            w.ints([j2, n, len(lev), len(plev)])
            w.ints([x for J, v in plev for x in (J, v)])
            for i, (J, v) in enumerate(lev):
                w.ints([J, v])
                w.reals(cfp[i])
    print('wrote %s (%s): %d LS records, %d jj records, %d bytes'
          % (path, 'text' if text else 'binary',
             len(ls), len(jj), os.path.getsize(path)))


if __name__ == '__main__':
    argv = sys.argv[1:]
    text = False
    if argv and argv[0] in ('-t', '--text'):
        text = True
        argv = argv[1:]
    cowan_cfp.install()
    # The default target is CODE/, the parent of this package -- that is
    # where RCG opens LSJJ.DAT, and writing it into CODE/lsjj/ beside this
    # script would put it somewhere nothing looks.  make_cfp names the
    # target explicitly and does not rely on this.
    here = os.path.dirname(os.path.abspath(__file__))
    default = os.path.join(os.path.dirname(here), 'LSJJ.DAT')
    write(argv[0] if argv else default, text)
