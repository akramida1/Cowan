#!/usr/bin/env python3
"""Compare the Fortran LS->jj recursion against recur.py, element by element.

    python cmpfort.py FDUMP.TXT

FDUMP.TXT is what LSJJSF writes.  Rows are matched on (deck term index,
2J) and columns on (N1, v1, 2J1, N2, v2, 2J2, 2J), so a permutation inside
a J block would be caught, not hidden.
"""

import sys

import numpy as np

import cowan_cfp
import recur
from recur import Transform, jj_cfp

_LS = 'spdfgh'


def read_dump(path):
    """-> [(l, n, orth, cross, {2J: (rows, cols, matrix)})]"""
    out = []
    with open(path) as fh:
        lines = fh.read().splitlines()
    i, nl = 0, len(lines)
    while i < nl:
        ln = lines[i]
        if not ln.startswith(' SHELL'):
            i += 1
            continue
        l = _LS.index(ln[7])
        n = int(ln[9:11])
        orth = float(ln[ln.index('ORTH') + 4:ln.index('1V2')])
        cross = float(ln[ln.index('1V2') + 3:])
        i += 1
        blocks = {}
        while i < nl and lines[i].startswith(' 2J='):
            j2 = int(lines[i][4:7])
            nr = int(lines[i][7:12])
            i += 1
            rows = []
            while lines[i].startswith(' R'):
                rows.append(tuple(int(lines[i][2 + 5 * k:7 + 5 * k])
                                  for k in range(2)))
                i += 1
            cols = []
            while lines[i].startswith(' C'):
                cols.append(tuple(int(lines[i][2 + 5 * k:7 + 5 * k])
                                  for k in range(5)))
                i += 1
            vals = []
            need = nr * len(cols)
            while len(vals) < need:
                vals += [float(x) for x in lines[i].split()]
                i += 1
            m = np.array(vals).reshape(len(cols), nr).T
            blocks[j2] = (rows, cols, m)
        out.append((l, n, orth, cross, blocks))
    return out


def compare(rec):
    l, n, orth, cross, blocks = rec
    tr = Transform(l, n)
    jm2, jp2 = tr.jm2, tr.jp2
    worst = 0.0
    if sorted(blocks) != sorted(J for J in tr.T[n] if tr.T[n][J].size):
        return None, 'J blocks differ'
    for j2, (rows, cols, m) in blocks.items():
        pr = tr.ls[n][j2]
        pc = tr.jj[n][j2]
        if len(pr) != len(rows) or len(pc) != len(cols):
            return None, '2J=%d: %dx%d vs %dx%d' % (j2, len(rows), len(cols),
                                                    len(pr), len(pc))
        for i, (t, jj) in enumerate(rows):
            if (t - 1, jj) != pr[i]:
                return None, '2J=%d row %d: %s vs %s' % (j2, i, (t - 1, jj),
                                                         pr[i])
        for k, (n1, k1, n2, k2, jj) in enumerate(cols):
            j1, v1 = jj_cfp(jm2, n1).levels[k1 - 1]
            j2s, v2 = jj_cfp(jp2, n2).levels[k2 - 1]
            if (n1, v1, j1, n2, v2, j2s) != pc[k]:
                return None, '2J=%d col %d: %s vs %s' % (
                    j2, k, (n1, v1, j1, n2, v2, j2s), pc[k])
        worst = max(worst, float(np.abs(m - tr.T[n][j2]).max()))
    return worst, ''


def main(argv):
    cowan_cfp.install()
    recs = read_dump(argv[0] if argv else 'FDUMP.TXT')
    print('')
    print(' Fortran LS->jj recursion against recur.py')
    print('  shell  dim      orthogonal   1 vs 2    vs recur.py')
    print('-' * 60)
    npass = nfail = 0
    for rec in recs:
        l, n, orth, cross, blocks = rec
        dim = sum(len(b[0]) for b in blocks.values())
        dev, msg = compare(rec)
        ok = dev is not None and dev < 1e-9
        print('  %-5s %5d    %8.1e  %8.1e   %8.1e   %s'
              % ('%s^%d' % (_LS[l], n), dim, orth, cross,
                 dev if dev is not None else float('nan'),
                 'ok' if ok else 'FAIL ' + msg))
        if ok:
            npass += 1
        else:
            nfail += 1
    print('-' * 60)
    print(' %d passed, %d failed' % (npass, nfail))
    return 0 if not nfail else 1


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
