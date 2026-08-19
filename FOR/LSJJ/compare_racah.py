"""Cross-check the two LS constructions against each other.

The Racah side (lsjj.LSShell) diagonalizes the pair operator and, for f, the
Casimirs of SO(7) and G2.  The Cowan side (cowan_cfp.CowanLSShell) recurses on
Cowan's CFP decks.  They are completely independent: no CFPs in the first, no
diagonalization in the second.

Both must span the same space, so

    O = T_cowan @ T_racah^T

is orthogonal, and it must be block diagonal with respect to (L, S, v) and,
for f, also Racah's (W) and (U): a Cowan state may mix only with Racah states
carrying the same labels.  Where those labels are complete -- all of p and d,
and every f term outside the groups Racah's classification cannot separate --
each block is 1x1, so O is a signature matrix and the two constructions agree
state for state up to sign.

The (W) and (U) of a Cowan state are not assumed: they are measured, as the
expectation values of the two Casimir operators in the CFP-built state.  The
"Casimir spread" column is the largest r.m.s. deviation from an eigenstate
over the whole shell; it staying at the 1e-4 level of the 10-decimal CFP
cards is the statement that Cowan's alpha respects Racah's classification.

Usage:  python compare_racah.py [p|d|f|all]
"""

import sys
import time
import numpy as np

import lsjj
import cowan_cfp
from lsjj import LSShell, JJShell, jj_det_to_ls, dot, TOL, LSYM


def _split_suffix(suf):
    """The Racah suffix ' 210/20' -> W = '210', U = '20'."""
    s = suf.strip()
    return tuple(s.split('/')) if s else ('', '')


def ls_rows(shell, M2, seniority, racah=None):
    """(label, J2, key, state) for one LS basis at fixed M2.

    `key` is (L2, S2, v, W, U): the labels a Cowan state and a Racah state
    must share before they are allowed to overlap.
    """
    out = []
    for (L2, S2), built in sorted(shell.terms.items()):
        for J2 in range(abs(L2 - S2), L2 + S2 + 1, 2):
            if abs(M2) > J2:
                continue
            for alpha, suf, degen, col in built:
                st = {}
                for ML2 in range(-L2, L2 + 1, 2):
                    MS2 = M2 - ML2
                    if abs(MS2) > S2:
                        continue
                    c = lsjj.cg(L2, ML2, S2, MS2, J2, M2)
                    if abs(c) < TOL:
                        continue
                    for k, val in col[(ML2, MS2)].items():
                        st[k] = st.get(k, 0.0) + c * val
                st = {k: v for k, v in st.items() if abs(v) > TOL}
                if not st:
                    continue
                if seniority is None:            # Racah side
                    v = alpha
                    W, U = _split_suffix(suf)
                else:                            # Cowan side
                    key = (S2 + 1, LSYM[L2 // 2], alpha)
                    v = seniority[key]
                    W, U = racah.get(key, ('', '', 0.0))[:2]
                out.append((shell.format_label(L2, S2, alpha, suf, J2),
                            J2, (L2, S2, v, W, U), lsjj.normalize(st)))
    return out


def compatible(ck, rk):
    """May a Cowan state and a Racah state overlap?

    (L, S, v) must agree.  W and U are compared only where both sides carry
    them: the Racah builder attaches them only inside a space that seniority
    left degenerate, since elsewhere they are not needed to fix the basis.
    """
    if ck[:3] != rk[:3]:
        return False
    for a, b in ((ck[3], rk[3]), (ck[4], rk[4])):
        if a and b and a != b:
            return False
    return True


def compare(l, N, out=print):
    M2 = 0 if N % 2 == 0 else 1

    lsjj.set_ls_builder(LSShell)
    rac = lsjj.ls_shell(l, N)
    rac_rows = ls_rows(rac, M2, None)

    lsjj.set_ls_builder(cowan_cfp.CowanLSShell)
    cow = lsjj.ls_shell(l, N)
    rl = cow.racah_labels()
    cow_rows = ls_rows(cow, M2, cow.seniority(), rl)
    sigma = max([v[2] for v in rl.values()] or [0.0])

    ok = True
    worst_offblock = 0.0
    ndiag = nblock = 0

    for J2 in sorted({r[1] for r in rac_rows}):
        R = [r for r in rac_rows if r[1] == J2]
        C = [r for r in cow_rows if r[1] == J2]
        if len(R) != len(C):
            out('   *** J=%s: %d Racah rows vs %d Cowan rows'
                % (lsjj._jlab(J2), len(R), len(C)))
            ok = False
            continue
        O = np.array([[dot(c[3], r[3]) for r in R] for c in C])
        if np.abs(O @ O.T - np.eye(len(C))).max() > 1e-8:
            out('   *** J=%s: overlap matrix not orthogonal' % lsjj._jlab(J2))
            ok = False
        for i, c in enumerate(C):
            same = [j for j, r in enumerate(R)
                    if compatible(c[2], r[2])]
            for j in range(len(R)):
                if j not in same:
                    worst_offblock = max(worst_offblock, abs(O[i, j]))
            nblock += 1
            if len(same) == 1:
                ndiag += 1
    if worst_offblock > 1e-8:
        out('   *** mixing across (L,S,v,W,U): %.2e' % worst_offblock)
        ok = False
    return ok, worst_offblock, ndiag, nblock, sigma


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else 'all'
    shells = {'p': (1, 6), 'd': (2, 10), 'f': (3, 14)}
    todo = ([shells[c] for c in which] if which != 'all'
            else [shells['p'], shells['d'], shells['f']])

    lines = []
    out = lines.append
    out('')
    out('Cowan-CFP basis vs Racah-classified basis')
    out('  O = T_cowan T_racah^T must not mix different (L, S, v, W, U)')
    out('-' * 74)
    out('  shell   max |O| outside          pinned 1-to-1   Casimir      time')
    out('          (L,S,v,W,U)                              spread')
    npass = nfail = 0
    for l, full in todo:
        for N in range(1, full + 1):
            t0 = time.time()
            try:
                ok, dev, nd, nb, sig = compare(l, N, out)
            except Exception as exc:
                out('  %s^%-2d  *** %s: %s' % ('spdf'[l], N,
                                               type(exc).__name__, exc))
                nfail += 1
                continue
            out('  %s^%-2d    %8.1e              %5d / %-5d %8.1e %7.2f s'
                % ('spdf'[l], N, dev, nd, nb, sig, time.time() - t0))
            npass += ok
            nfail += not ok
    out('')
    out('=' * 74)
    out(' %d passed, %d failed' % (npass, nfail))
    text = '\n'.join(lines)
    print(text)
    with open('compare_%s.log' % which, 'w', newline='') as fh:
        fh.write(text + '\n')
    return 0 if nfail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
