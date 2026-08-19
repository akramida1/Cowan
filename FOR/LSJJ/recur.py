#!/usr/bin/env python3
"""
LS -> jj transformation by recursion on N, from coefficients of fractional
parentage alone.  No Slater determinants.

Why this exists
---------------
lsjj.py computes the transformation by building both bases in the m-scheme
and taking overlaps.  That is the right way to *establish* the numbers -- it
has no phase conventions to reconcile -- but it is not the way to compute
them inside RCG: f^7 alone has 3432 determinants, and the m-scheme route
needs a null-space solver and an eigensolver.

The recursion below needs neither.  It uses

  * Cowan's LS cfp deck, which RCG already reads (CUVFD, rcg11kd.f:2136),
  * a jj cfp table for j = 1/2, 3/2, 5/2, 7/2, generated once by jj_cfp()
    below and small enough to carry as data,
  * Clebsch-Gordan coefficients, i.e. what RCG's S9J/S6J already provide.

Everything else is small dense matrix products, blocked by J.

The recursion
-------------
Write  T_n( a L S J ; ( j1^N1 c1 J1, j2^N2 c2 J2 ) J )  for the overlap of the
two coupled states of n = N1 + N2 electrons.  Peel one creation operator off
the jj side, in turn from each subshell, and use

    a_{jm} | l^n a L S J M >
        = sqrt(n) sum_{ab,Jb} cfp(ab|a) Y(ab;j,Jb) < j m, Jb Mb | J M >
                              | l^(n-1) ab Jb Mb >

which is the LS cfp seen through the recoupling

    Y = < ( (l s) j , (Lb Sb) Jb ) J | ( (Lb l) L , (Sb s) S ) J >

(a 9j).  The result is

    T_n = sqrt(n) * sum_{k=1,2} sigma_k / sqrt(N_k)
              * sum_{ck-bar} cfp_jj( ck-bar | ck )
              * sum_{Jb} R_k(Jb)
              * sum_{ab} cfp_LS( ab | a ) Y( ab ; j_k, Jb ) T_(n-1)

with sigma_1 = 1, sigma_2 = (-1)^N1 (the creation operator of subshell 2 has
to be carried left past the N1 operators of subshell 1), and R_k the 6j-type
recoupling of the spectator subshell,

    R_1 = < ( j1, (J1b J2) Jb ) J | ( (j1 J1b) J1, J2 ) J >
    R_2 = < ( j2, (J1 J2b) Jb ) J | ( J1, (j2 J2b) J2 ) J >.

Conventions
-----------
Both cfp definitions put the creation operator on the left and the parent
first in the angular coupling:

    | l^n a L S; ML MS > = (1/sqrt n) sum_ab cfp
            sum_(m,mu) < Lb MLb, l m | L ML > < Sb MSb, 1/2 mu | S MS >
                       a+_(m,mu) | l^(n-1) ab; MLb MSb >

    | j^n v J; M >       = (1/sqrt n) sum_vb cfp
            sum_m < j m, Jb Mb | J M >  a+_(j,m) | j^(n-1) vb Jb; Mb >

The first is exactly cowan_cfp.py's, so the LS deck is used unchanged.  The
second fixes the jj phases; jj_cfp() extracts it from the m-scheme subshells
of lsjj.py, so the two routes are automatically in phase and can be compared
element by element.

The recoupling coefficients are evaluated by explicit Clebsch-Gordan sums
rather than by 6j/9j formulas.  They are small, and a summed-CG definition
cannot carry a phase-convention error.  A Fortran port would substitute
S6J/S9J and check against these.

Command line
------------
    python recur.py check l N     one shell against the m-scheme route
    python recur.py check all     p, d and f, every N
    python recur.py jjcfp         print the jj cfp tables
    python recur.py deck file     write the jj cfp deck
"""

import sys
from math import sqrt

import numpy as np

import lsjj
from lsjj import cg, dot, TOL, LSYM, JSubshell, _create, _jlab
import cowan_cfp
from cowan_cfp import decks, _term_ls, LSYMB


# --------------------------------------------------------------------------
#  recoupling coefficients, by explicit CG sums
# --------------------------------------------------------------------------

_YC = {}


def ycoef(l2, Lb2, Sb2, L2, S2, J2, Jb2, j2):
    """< ((l s) j, (Lb Sb) Jb) J | ((Lb l) L, (Sb s) S) J >."""
    key = (l2, Lb2, Sb2, L2, S2, J2, Jb2, j2)
    hit = _YC.get(key)
    if hit is not None:
        return hit
    M2 = J2
    tot = 0.0
    for MLb in range(-Lb2, Lb2 + 1, 2):
        for MSb in range(-Sb2, Sb2 + 1, 2):
            Mb = MLb + MSb
            if abs(Mb) > Jb2:
                continue
            c1 = cg(Lb2, MLb, Sb2, MSb, Jb2, Mb)
            if abs(c1) < TOL:
                continue
            for ml in range(-l2, l2 + 1, 2):
                for ms in (1, -1):
                    mj = ml + ms
                    if mj + Mb != M2 or abs(mj) > j2:
                        continue
                    a = (cg(l2, ml, 1, ms, j2, mj) * c1
                         * cg(j2, mj, Jb2, Mb, J2, M2))
                    if abs(a) < TOL:
                        continue
                    ML, MS = MLb + ml, MSb + ms
                    b = (cg(Lb2, MLb, l2, ml, L2, ML)
                         * cg(Sb2, MSb, 1, ms, S2, MS)
                         * cg(L2, ML, S2, MS, J2, M2))
                    tot += a * b
    _YC[key] = tot
    return tot


_R1C = {}


def r1coef(j2, Jb1, Jt2, J1, Jb, J):
    """< ( j1, (J1b J2) Jb ) J | ( (j1 J1b) J1, J2 ) J >."""
    key = (j2, Jb1, Jt2, J1, Jb, J)
    hit = _R1C.get(key)
    if hit is not None:
        return hit
    M, tot = J, 0.0
    for m in range(-j2, j2 + 1, 2):
        for M1b in range(-Jb1, Jb1 + 1, 2):
            M2s = M - m - M1b
            if abs(M2s) > Jt2:
                continue
            Mb = M1b + M2s
            if abs(Mb) > Jb:
                continue
            a = cg(Jb1, M1b, Jt2, M2s, Jb, Mb) * cg(j2, m, Jb, Mb, J, M)
            if abs(a) < TOL:
                continue
            M1 = m + M1b
            b = cg(j2, m, Jb1, M1b, J1, M1) * cg(J1, M1, Jt2, M2s, J, M)
            tot += a * b
    _R1C[key] = tot
    return tot


_R2C = {}


def r2coef(j2, J1, Jb2, Js2, Jb, J):
    """< ( j2, (J1 J2b) Jb ) J | ( J1, (j2 J2b) J2 ) J >."""
    key = (j2, J1, Jb2, Js2, Jb, J)
    hit = _R2C.get(key)
    if hit is not None:
        return hit
    M, tot = J, 0.0
    for m in range(-j2, j2 + 1, 2):
        for M1 in range(-J1, J1 + 1, 2):
            Mb2 = M - m - M1
            if abs(Mb2) > Jb2:
                continue
            Mb = M1 + Mb2
            if abs(Mb) > Jb:
                continue
            a = cg(J1, M1, Jb2, Mb2, Jb, Mb) * cg(j2, m, Jb, Mb, J, M)
            if abs(a) < TOL:
                continue
            M2s = m + Mb2
            b = cg(j2, m, Jb2, Mb2, Js2, M2s) * cg(J1, M1, Js2, M2s, J, M)
            tot += a * b
    _R2C[key] = tot
    return tot


# --------------------------------------------------------------------------
#  jj coefficients of fractional parentage
# --------------------------------------------------------------------------

class JSub(object):
    """Terms and cfp of one j^n subshell.

    levels   list of (J2, v), in the order lsjj.JSubshell produces them
    cfp      len(levels) x len(parent levels), cfp[i, k]
    """

    __slots__ = ('j2', 'n', 'levels', 'cfp')

    def __init__(self, j2, n, levels, cfp):
        self.j2, self.n, self.levels, self.cfp = j2, n, levels, cfp


_JSUB = {}


def jj_cfp(j2, n):
    """j^n terms and their cfp, extracted from the m-scheme subshell.

    Uses  < j^n vJ M | sum_m < j m, Jb Mb | J M > a+_(j,m) | j^(n-1) vb Jb Mb >
          = sqrt(n) * cfp,
    the one-subshell case of the identity quoted in the module docstring.
    Each row is checked to be normalized, which is the statement that the
    parent list is complete.
    """
    key = (j2, n)
    hit = _JSUB.get(key)
    if hit is not None:
        return hit
    if n == 0:
        out = JSub(j2, 0, [(0, 0)], np.zeros((1, 0)))
        _JSUB[key] = out
        return out

    sub = JSubshell(j2, n)
    par = JSubshell(j2, n - 1) if n > 1 else None
    plev = ([(J, v) for J, v, _d, _c in par.levels] if par else [(0, 0)])
    pcol = ([col for _J, _v, _d, col in par.levels] if par
            else [{0: {0: 1.0}}])
    mvals = list(range(j2, -j2 - 1, -2))

    levels, rows = [], []
    for J2, v, _d, col in sub.levels:
        state = col[J2]
        row = np.zeros(len(plev))
        for k, (Jb2, _vb) in enumerate(plev):
            W = {}
            for m in mvals:
                Mb = J2 - m
                if abs(Mb) > Jb2:
                    continue
                c = cg(j2, m, Jb2, Mb, J2, J2)
                if abs(c) < TOL:
                    continue
                p = mvals.index(m)
                for mask, amp in pcol[k][Mb].items():
                    s, m2 = _create(mask, p)
                    if s:
                        W[m2] = W.get(m2, 0.0) + c * amp * s
            if W:
                row[k] = dot(state, W) / sqrt(n)
        nrm = float(row @ row)
        if abs(nrm - 1.0) > 1e-8:
            raise AssertionError('j=%s^%d %s: cfp row norm %.10f'
                                 % (_jlab(j2), n, _jlab(J2), nrm))
        levels.append((J2, v))
        rows.append(row)
    out = JSub(j2, n, levels, np.array(rows) if rows else np.zeros((0, 0)))
    _JSUB[key] = out
    return out


# --------------------------------------------------------------------------
#  the LS side of the recursion
# --------------------------------------------------------------------------

class LSChain(object):
    """LS terms of l^0 ... l^Nmax and the deck that links them."""

    def __init__(self, l, nmax):
        self.l = l
        self.terms = {0: [(1, 'S', '')]}
        self.pmap = {}          # n -> deck column k -> row index in l^(n-1)
        self.cfp = {}           # n -> array
        lsym = LSYMB[l]
        tab = decks()
        for n in range(1, nmax + 1):
            d = tab.get((lsym, n))
            if d is None:
                raise KeyError('no CFP deck for %s^%d' % (lsym, n))
            self.terms[n] = list(d.terms)
            self.pmap[n] = [self.terms[n - 1].index(p) for p in d.parents]
            self.cfp[n] = d.cfp

    def rows(self, n):
        """[(term index, 2J)] for l^n, grouped later by J."""
        out = []
        for t, lab in enumerate(self.terms[n]):
            L2, S2 = _term_ls(lab)
            for J2 in range(abs(L2 - S2), L2 + S2 + 1, 2):
                out.append((t, J2))
        return out


# --------------------------------------------------------------------------
#  the transformation itself
# --------------------------------------------------------------------------

class Transform(object):
    """T_n for n = 0 ... N, blocked by total J.

    ls[n][J2]   list of (term index, 2J)
    jj[n][J2]   list of (N1, v1, 2J1, N2, v2, 2J2)
    T[n][J2]    array, rows over ls[n][J2], columns over jj[n][J2]
    """

    def __init__(self, l, N, verbose=False):
        self.l, self.N = l, N
        self.jm2 = 2 * l - 1 if l > 0 else -1
        self.jp2 = 2 * l + 1
        self.chain = LSChain(l, N)
        self.ls, self.jj, self.T = {}, {}, {}
        self.lsidx, self.jjidx = {}, {}
        self._build(verbose)

    # -- state lists -------------------------------------------------------
    def _ls_lists(self, n):
        by = {}
        for t, J2 in self.chain.rows(n):
            by.setdefault(J2, []).append((t, J2))
        self.ls[n] = by
        self.lsidx[n] = {J2: {r: i for i, r in enumerate(v)}
                         for J2, v in by.items()}

    def _jj_lists(self, n):
        by = {}
        nm_max = self.jm2 + 1 if self.jm2 > 0 else 0
        for N1 in range(0, min(nm_max, n) + 1):
            N2 = n - N1
            if N2 > self.jp2 + 1:
                continue
            s1 = jj_cfp(self.jm2, N1) if N1 else None
            l1 = s1.levels if s1 else [(0, 0)]
            s2 = jj_cfp(self.jp2, N2)
            for J1, v1 in l1:
                for J2s, v2 in s2.levels:
                    for J2 in range(abs(J1 - J2s), J1 + J2s + 1, 2):
                        by.setdefault(J2, []).append(
                            (N1, v1, J1, N2, v2, J2s))
        self.jj[n] = by
        self.jjidx[n] = {J2: {r: i for i, r in enumerate(v)}
                         for J2, v in by.items()}

    # -- the A matrix:  cfp_LS * Y ----------------------------------------
    def _amat(self, n, J2, Jb2, j2):
        key = (n, J2, Jb2, j2)
        hit = self._acache.get(key)
        if hit is not None:
            return hit
        rows = self.ls[n].get(J2, [])
        prow = self.ls[n - 1].get(Jb2, [])
        A = np.zeros((len(rows), len(prow)))
        if rows and prow:
            pidx = self.lsidx[n - 1][Jb2]
            cfp = self.chain.cfp[n]
            pmap = self.chain.pmap[n]
            l2 = 2 * self.l
            for i, (t, _) in enumerate(rows):
                L2, S2 = _term_ls(self.chain.terms[n][t])
                for k, pt in enumerate(pmap):
                    c = cfp[t, k]
                    if abs(c) < 1e-12:
                        continue
                    jpos = pidx.get((pt, Jb2))
                    if jpos is None:
                        continue
                    Lb2, Sb2 = _term_ls(self.chain.terms[n - 1][pt])
                    y = ycoef(l2, Lb2, Sb2, L2, S2, J2, Jb2, j2)
                    if abs(y) > TOL:
                        A[i, jpos] += c * y
        self._acache[key] = A
        return A

    # -- driver ------------------------------------------------------------
    def _build(self, verbose):
        self._acache = {}
        self.cross = 0.0
        self._ls_lists(0)
        self._jj_lists(0)
        self.T[0] = {0: np.ones((1, 1))}

        for n in range(1, self.N + 1):
            self._ls_lists(n)
            self._jj_lists(n)
            self.T[n] = {}
            for J2, cols in self.jj[n].items():
                rows = self.ls[n].get(J2, [])
                M = np.zeros((len(rows), len(cols)))
                for c, key in enumerate(cols):
                    self._column(n, J2, c, key, M)
                self.T[n][J2] = M
            if verbose:
                print('  n=%2d  %d J blocks' % (n, len(self.T[n])))
        self._acache = {}

    def _column(self, n, J2, c, key, M):
        """One column of T_n.

        Peeling a creation operator off subshell 1 and peeling one off
        subshell 2 are each, by themselves, a complete rewriting of the jj
        state -- the two subshells are separate factors, so there is no sum
        over "which subshell the last electron came from" as there is over
        "which electron" inside one shell.  Either route therefore gives the
        whole column.  Both are computed whenever both subshells are
        occupied and the difference is recorded in self.cross; it is a free
        and very sharp check on the recoupling coefficients and on the
        (-1)^N1 that carries a+ of subshell 2 past subshell 1.
        """
        N1, v1, J1, N2, v2, Js2 = key
        got = None
        for k in (1, 2):
            Nk = N1 if k == 1 else N2
            if Nk == 0:
                continue
            acc = np.zeros(M.shape[0])
            jk = self.jm2 if k == 1 else self.jp2
            sub = jj_cfp(jk, Nk)
            i = sub.levels.index((J1 if k == 1 else Js2,
                                  v1 if k == 1 else v2))
            par = jj_cfp(jk, Nk - 1)
            sig = 1.0 if k == 1 else (-1.0) ** N1
            pref = sqrt(float(n)) * sig / sqrt(float(Nk))
            for kb, (Jbk, vbk) in enumerate(par.levels):
                cj = sub.cfp[i, kb]
                if abs(cj) < 1e-12:
                    continue
                if k == 1:
                    other = Js2
                    bkey = (N1 - 1, vbk, Jbk, N2, v2, Js2)
                else:
                    other = J1
                    bkey = (N1, v1, J1, N2 - 1, vbk, Jbk)
                for Jb in range(abs(Jbk - other), Jbk + other + 1, 2):
                    if abs(Jb - jk) > J2 or Jb + jk < J2:
                        continue
                    pos = self.jjidx[n - 1].get(Jb, {}).get(bkey)
                    if pos is None:
                        continue
                    r = (r1coef(jk, Jbk, other, J1, Jb, J2) if k == 1
                         else r2coef(jk, other, Jbk, Js2, Jb, J2))
                    if abs(r) < TOL:
                        continue
                    A = self._amat(n, J2, Jb, jk)
                    if not A.size:
                        continue
                    acc += pref * cj * r * (A @ self.T[n - 1][Jb][:, pos])
            if got is None:
                got = acc
                M[:, c] = acc
            else:
                self.cross = max(self.cross, float(np.abs(acc - got).max()))

    # -- labels ------------------------------------------------------------
    def ls_label(self, n, row):
        t, J2 = row
        mult, lsy, alpha = self.chain.terms[n][t]
        return '%d%s%s)%s' % (mult, lsy, alpha, _jlab(J2))

    def jj_label(self, key, J2):
        N1, v1, J1, N2, v2, Js2 = key
        return '(%s^%d %d %s, %s^%d %d %s)%s' % (
            _jlab(self.jm2), N1, v1, _jlab(J1),
            _jlab(self.jp2), N2, v2, _jlab(Js2), _jlab(J2))


# --------------------------------------------------------------------------
#  checks
# --------------------------------------------------------------------------

def check_orthogonal(tr, out=print):
    """Every J block of T must be orthogonal."""
    dev = 0.0
    for J2, M in sorted(tr.T[tr.N].items()):
        if M.shape[0] != M.shape[1]:
            out(' *** J=%s block is %dx%d' % (_jlab(J2), M.shape[0],
                                              M.shape[1]))
            return None
        dev = max(dev, np.abs(M @ M.T - np.eye(M.shape[0])).max())
    return dev


def check_against_mscheme(l, N, out=print):
    """Compare with lsjj's determinant route, element by element."""
    cowan_cfp.install()
    tr = Transform(l, N)
    ls_states, jj_states = lsjj.transformation(l, N)

    worst = 0.0
    for J2 in sorted(tr.T[N]):
        rows = tr.ls[N].get(J2, [])
        cols = tr.jj[N][J2]
        rlab = [tr.ls_label(N, r) for r in rows]
        clab = [tr.jj_label(k, J2) for k in cols]
        L = [s for s in ls_states if s[1] == J2]
        Jt = [s for s in jj_states if s[1] == J2]
        if len(L) != len(rlab) or len(Jt) != len(clab):
            out(' *** J=%s: %d/%d rows, %d/%d columns'
                % (_jlab(J2), len(rlab), len(L), len(clab), len(Jt)))
            return None
        rpos = {s[0].strip(): i for i, s in enumerate(L)}
        cpos = {s[0].strip(): i for i, s in enumerate(Jt)}
        for i, lab in enumerate(rlab):
            if lab not in rpos:
                out(' *** J=%s: LS state %r not in the m-scheme list'
                    % (_jlab(J2), lab))
                return None
        for j, lab in enumerate(clab):
            if lab not in cpos:
                out(' *** J=%s: jj state %r not in the m-scheme list'
                    % (_jlab(J2), lab))
                return None
        ref = np.zeros((len(rlab), len(clab)))
        for i, la in enumerate(rlab):
            a = L[rpos[la]][3]
            for j, lb in enumerate(clab):
                ref[i, j] = dot(a, Jt[cpos[lb]][3])
        worst = max(worst, np.abs(tr.T[N][J2] - ref).max())
    return worst


def run(l, N, out=print):
    import time
    name = 'spdfgh'[l]
    t0 = time.time()
    cowan_cfp.install()
    tr = Transform(l, N)
    t1 = time.time()
    orth = check_orthogonal(tr, out)
    dev = check_against_mscheme(l, N, out)
    t2 = time.time()
    ok = (orth is not None and orth < 1e-9
          and dev is not None and dev < 1e-9
          and tr.cross < 1e-9)
    ndim = sum(M.shape[0] for M in tr.T[N].values())
    out('  %-5s %5d states   %8.1e  %8.1e  %8.1e   %6.2f %6.2f s   %s'
        % ('%s^%d' % (name, N), ndim,
           orth if orth is not None else float('nan'),
           tr.cross,
           dev if dev is not None else float('nan'),
           t1 - t0, t2 - t1, 'ok' if ok else 'FAIL'))
    return ok


def main(argv):
    argv = list(argv)
    if 'fromdeck' in argv:
        argv.remove('fromdeck')
        import os
        read_deck(os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               'JJCFP'))
        print('jj cfp taken from the deck, not from the m-scheme')
    if not argv or argv[0] == 'check':
        which = argv[1] if len(argv) > 1 else 'all'
        if which == 'all':
            cases = ([(1, n) for n in range(1, 7)]
                     + [(2, n) for n in range(1, 11)]
                     + [(3, n) for n in range(1, 15)])
        else:
            cases = [(int(argv[1]), int(argv[2]))]
        print('')
        print('Recursive LS->jj transformation, checked against the m-scheme')
        print('  shell  dim         orthogonal  1 vs 2  vs m-scheme'
              '     build  check')
        print('-' * 78)
        npass = nfail = 0
        for l, n in cases:
            if run(l, n):
                npass += 1
            else:
                nfail += 1
        print('-' * 78)
        print(' %d passed, %d failed' % (npass, nfail))
        return 0 if not nfail else 1
    if argv[0] == 'jjcfp':
        for j2 in (1, 3, 5, 7):
            print('')
            print(' j = %s' % _jlab(j2))
            for n in range(1, j2 + 2):
                s = jj_cfp(j2, n)
                p = jj_cfp(j2, n - 1)
                print('  n=%d  terms %s' % (
                    n, ' '.join('%sv%d' % (_jlab(J), v) for J, v in s.levels)))
                for i, (J, v) in enumerate(s.levels):
                    row = '   %sv%-2d :' % (_jlab(J), v)
                    for k, (Jb, vb) in enumerate(p.levels):
                        if abs(s.cfp[i, k]) > 1e-12:
                            row += '  %sv%d %10.7f' % (_jlab(Jb), vb,
                                                       s.cfp[i, k])
                    print(row)
        return 0
    if argv[0] == 'deck':
        write_deck(argv[1] if len(argv) > 1 else None)
        return 0
    print(__doc__)
    return 0


def read_deck(path):
    """Load a jj cfp deck written by write_deck() into the cache.

    After this the recursion runs from the two decks and Clebsch-Gordan
    coefficients alone: no determinants, no eigensolver.  That is the whole
    point -- it is what the Fortran side will do.
    """
    with open(path) as fh:
        lines = [ln for ln in fh.read().splitlines() if ln.strip()]
    i = 0
    while i < len(lines):
        j2 = int(lines[i][0:4])
        n = int(lines[i][4:8])
        nch = int(lines[i][8:12])
        npa = int(lines[i][12:16])
        i += 1
        plev = [(int(lines[i][6 * k:6 * k + 4]),
                 int(lines[i][6 * k + 4:6 * k + 6])) for k in range(npa)]
        i += 1
        levels, rows = [], []
        for _ in range(nch):
            rec = lines[i]
            i += 1
            levels.append((int(rec[0:4]), int(rec[4:6])))
            rows.append([float(rec[6 + 16 * k:22 + 16 * k])
                         for k in range(npa)])
        _JSUB[(j2, 0)] = JSub(j2, 0, [(0, 0)], np.zeros((1, 0)))
        if n - 1 > 0 and (j2, n - 1) not in _JSUB:
            raise ValueError('deck out of order at j=%s^%d' % (_jlab(j2), n))
        if n - 1 > 0 and _JSUB[(j2, n - 1)].levels != plev:
            raise ValueError('parent list of j=%s^%d does not match j=%s^%d'
                             % (_jlab(j2), n, _jlab(j2), n - 1))
        _JSUB[(j2, n)] = JSub(j2, n, levels, np.array(rows))
    return len(_JSUB)


def write_deck(path=None):
    """The jj cfp tables, in a fixed-format deck for the Fortran side."""
    lines = []
    for j2 in (1, 3, 5, 7):
        for n in range(1, j2 + 2):
            s = jj_cfp(j2, n)
            p = jj_cfp(j2, n - 1)
            lines.append('%4d%4d%4d%4d' % (j2, n, len(s.levels),
                                           len(p.levels)))
            lines.append(''.join('%4d%2d' % (J, v) for J, v in p.levels))
            for i, (J, v) in enumerate(s.levels):
                lines.append('%4d%2d' % (J, v)
                             + ''.join('%16.10f' % x for x in s.cfp[i]))
    text = '\n'.join(lines) + '\n'
    if path:
        with open(path, 'w', newline='') as fh:
            fh.write(text)
        print('wrote %s (%d lines)' % (path, len(lines)))
    else:
        print(text)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))
