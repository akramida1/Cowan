#!/usr/bin/env python3
"""
LS basis for l^N built from Cowan's own coefficients of fractional parentage.

Why
---
The LS side of lsjj.py is built by diagonalizing the pair operator and the
Casimir operators of SO(7) and G2.  That fixes the basis completely for p^N
and d^N, but for f^N Racah's labels are incomplete: 8 to 15 states per shell
share every label, and inside such a group the basis is an arbitrary
orthogonal rotation.  A transformation built that way is orthogonal and
M-independent but its rows do not correspond, state for state, to the ones
RCG prints.

This module removes the ambiguity by building the LS states the way Cowan
does, from the CFP decks that RCG itself reads (unit IR in CUVFD,
rcg11kd.f:2136-2171).  The recursion

    | l^n a L S ; ML MS >  =  sum_ab  ( l^(n-1) ab | } l^n a )
                             * sum_(m,mu) < Lb MLb, l m | L ML >
                                          < Sb MSb, 1/2 mu | S MS >
                                          a+_(m,mu) | l^(n-1) ab ; MLb MSb >

reproduces Cowan's states including his alpha labels and his signs, because
the cfp table is the definition of both.

Two facts make this self-checking, and both are asserted at build time:

  * the vector built above has norm exactly sqrt(n) -- the coupled creation
    operator satisfies < l^n aLS | a+ | l^(n-1) ab Lb Sb >^(LS) = sqrt(n) cfp,
    so cfp row orthonormality collapses the sum onto a single state;
  * within one (ML,MS) block the states are orthonormal and complete, i.e.
    their number equals the number of Slater determinants in the block.

If a deck were misread, either check would fail immediately.

Structure
---------
Everything is blocked by (ML, MS).  Determinants of l^n at fixed (ML,MS) are
few (at most ~120 even for f^7), so each block is a small dense matrix and
the whole recursion is a few hundred small matrix products.  Building all of
f^1 ... f^14 takes a few seconds.

The overall phase (-1)^(n-1) between "electron n in slot n" (Cowan) and
a+ acting on the left (here) is uniform over a shell, so it cancels out of
every relative sign and is not applied.
"""

import os
import re
import numpy as np
from itertools import combinations

import lsjj
from lsjj import LSShell, cg, TOL, LSYM

LSYMB = 'SPDFGHIKLMNOQRTUVWXYZABC'


# --------------------------------------------------------------------------
#  locating the data files
# --------------------------------------------------------------------------

_HERE = os.path.dirname(os.path.abspath(__file__))

_CFP_CANDIDATES = [
    os.environ.get('LSJJ_CFP'),
    os.path.join(_HERE, 'ING11'),
    os.path.join(_HERE, '..', 'CFP', 'ING11'),
    r'C:\cowan\FOR\g95\Claude\CFP\ING11',
]

_SENIOR_CANDIDATES = [
    os.environ.get('LSJJ_SENIOR'),
    os.path.join(_HERE, 'SENIOR'),
    os.path.join(_HERE, '..', 'CFP', 'SENIOR'),
    r'C:\cowan\FOR\g95\Claude\CFP\SENIOR',
]


def _find(cands, what):
    for c in cands:
        if c and os.path.isfile(c):
            return os.path.abspath(c)
    raise IOError('%s not found; tried:\n  %s'
                  % (what, '\n  '.join(str(c) for c in cands)))


# --------------------------------------------------------------------------
#  the CFP deck
# --------------------------------------------------------------------------

class Deck(object):
    """One (l, n) CFP table, exactly as it sits on the cards.

    terms    list of (mult, Lsym, alpha) for l^n, in deck order
    parents  list of (mult, Lsym, alpha) for l^(n-1), in column order
    cfp      NOT x NOP array, cfp[i, k] = ( l^(n-1) parent_k | } l^n term_i )
    """

    __slots__ = ('lsym', 'n', 'terms', 'parents', 'cfp')

    def __init__(self, lsym, n, terms, parents, cfp):
        self.lsym, self.n = lsym, n
        self.terms, self.parents, self.cfp = terms, parents, cfp

    def __repr__(self):
        return '<Deck %s^%d %dx%d>' % (self.lsym, self.n,
                                       len(self.terms), len(self.parents))


_HEADER = re.compile(r'^   ([SPDFGHIKL])\s*[-0-9]')


def _label(text):
    """4-character term field  ->  (mult, Lsym, alpha)."""
    return (int(text[0]), text[1], text[2:4].strip())


def read_deck_file(path=None):
    """Parse Cowan's CFP card deck.  Returns {(Lsym, n): Deck}."""
    path = path or _find(_CFP_CANDIDATES, 'CFP deck ING11')
    with open(path) as fh:
        # cards are fixed format; pad so that short lines still have all
        # their (blank) fields
        lines = [ln.ljust(300) for ln in fh.read().splitlines()]

    decks, i, nline = {}, 0, len(lines)
    while i < nline:
        line = lines[i]
        if not _HEADER.match(line):
            i += 1
            continue
        lsym = line[3]
        n = int(line[4:8])
        not_ = int(line[8:12])
        nop = int(line[12:16])
        i += 1

        # FORMAT (20(I1,A1,A2)): 20 parents per card, then continuation cards
        parents = []
        while len(parents) < nop:
            pline = lines[i]
            i += 1
            for k in range(min(20, nop - len(parents))):
                parents.append(_label(pline[4 * k:4 * k + 4]))

        terms, cfp = [], np.zeros((not_, nop))
        for m in range(not_):
            row = np.zeros(nop)
            first, prev = None, 0
            while True:
                rec = lines[i]
                i += 1
                lab = _label(rec[0:4])
                it = int(rec[4:8] or 0)
                if first is None:
                    first = lab
                elif lab != first:
                    raise ValueError('%s^%d: card out of order at %r'
                                     % (lsym, n, rec))
                if abs(it) != prev + 1:
                    raise ValueError('%s^%d: continuation out of order at %r'
                                     % (lsym, n, rec))
                prev = abs(it)
                for g in range(4):
                    fi = rec[8 + 16 * g:11 + 16 * g]
                    fv = rec[11 + 16 * g:24 + 16 * g]
                    # IPNT = 0 is the pad that ends a card (CUVFD line 2158),
                    # not parent number 0
                    k = int(fi) if fi.strip() else 0
                    if k == 0:
                        continue
                    row[k - 1] = float(fv) if fv.strip() else 0.0
                if it <= 0:
                    break
            terms.append(first)
            cfp[m] = row
        decks[(lsym, n)] = Deck(lsym, n, terms, parents, cfp)
    return decks


def read_seniority_file(path=None):
    """Parse Cowan's SENIOR file.  Returns {(Lsym, n, mult, L, alpha): v}."""
    path = path or _find(_SENIOR_CANDIDATES, 'seniority file SENIOR')
    out, key = {}, None
    with open(path) as fh:
        for line in fh.read().splitlines():
            if not line.strip():
                continue
            if line[0].isalpha() and len(line.strip()) <= 3:
                key = (line.strip()[0], int(line.strip()[1:]))
                continue
            for k in range(0, len(line.rstrip()), 6):
                rec = line[k:k + 6]
                if len(rec) < 6:
                    break
                mult, lsy, alf = _label(rec[0:4])
                out[(key[0], key[1], mult, lsy, alf)] = int(rec[4:6])
    return out


_DECK_CACHE = {}


def decks():
    if 'd' not in _DECK_CACHE:
        _DECK_CACHE['d'] = read_deck_file()
    return _DECK_CACHE['d']


# --------------------------------------------------------------------------
#  the recursion
# --------------------------------------------------------------------------

def _term_ls(lab):
    """(mult, Lsym, alpha) -> (2L, 2S)."""
    return 2 * LSYMB.index(lab[1]), lab[0] - 1


class _ShellStates(object):
    """All | l^n a L S ; ML MS > for one n, blocked by (ML, MS).

    block[(ML2, MS2)] = (rows, V) with rows a list of indices into `terms`
    and V a (len(rows) x ndet) orthonormal matrix over the determinants of
    that block.
    """

    __slots__ = ('n', 'terms', 'block', 'dets', 'pos')

    def __init__(self, n, terms):
        self.n, self.terms = n, terms
        self.block, self.dets, self.pos = {}, {}, {}


# the recursion is deterministic, so a shell built once can be reused as the
# parent of the next one: (l, n) -> _ShellStates
_CHAIN_CACHE = {}


def _det_blocks(norb, n, ml2, ms2):
    """Determinant masks of n electrons, grouped by (ML2, MS2)."""
    blocks, pos = {}, {}
    for combo in combinations(range(norb), n):
        mask = 0
        a = b = 0
        for i in combo:
            mask |= 1 << i
            a += ml2[i]
            b += ms2[i]
        lst = blocks.setdefault((a, b), [])
        pos[mask] = (a, b, len(lst))
        lst.append(mask)
    return blocks, pos


class CowanLSShell(LSShell):
    """LSShell whose terms come from Cowan's CFP decks instead of from a
    diagonalization.  Same interface; `terms` carries Cowan's alpha label in
    place of the Racah suffix, and the residual degeneracy is always 1."""

    def __init__(self, l, N, deck_table=None):
        self._decks = deck_table if deck_table is not None else decks()
        self._seniority = None
        LSShell.__init__(self, l, N)

    # -- the part that replaces the diagonalization -----------------------
    def _build_terms(self):
        l, N = self.l, self.N
        lsym = LSYMB[l]
        norb = len(self.orbs)
        ml2 = self.ml2
        ms2 = self.ms2

        # n = 0: the closed-shell vacuum, one state 1S
        cur = _ShellStates(0, [(1, 'S', '')])
        cur.dets = {(0, 0): [0]}
        cur.pos = {0: (0, 0, 0)}
        cur.block = {(0, 0): ([0], np.ones((1, 1)))}

        for n in range(1, N + 1):
            hit = _CHAIN_CACHE.get((l, n))
            if hit is not None:
                cur = hit
                continue
            deck = self._decks.get((lsym, n))
            if deck is None:
                raise KeyError('no CFP deck for %s^%d' % (lsym, n))
            cur = self._step(cur, deck, n, norb, ml2, ms2)
            _CHAIN_CACHE[(l, n)] = cur

        self._shell = cur
        return self._export(cur)

    def _step(self, prev, deck, n, norb, ml2, ms2):
        """One rung of the recursion: l^(n-1) -> l^n."""
        # column k of the deck is parent deck.parents[k]; map it onto the
        # row index of the previous shell (CUVFD does the same at line 2181)
        pmap = []
        for lab in deck.parents:
            try:
                pmap.append(prev.terms.index(lab))
            except ValueError:
                raise ValueError('parent %s of %s^%d not among terms of %s^%d'
                                 % (lab, deck.lsym, n, deck.lsym, n - 1))
        inv = {t: k for k, t in enumerate(pmap)}

        new = _ShellStates(n, list(deck.terms))
        new.dets, new.pos = _det_blocks(norb, n, ml2, ms2)

        tls = [_term_ls(t) for t in new.terms]
        pls = [_term_ls(t) for t in prev.terms]

        for (ML2, MS2), masks in new.dets.items():
            rows = [i for i, (L2, S2) in enumerate(tls)
                    if L2 >= abs(ML2) and S2 >= abs(MS2)]
            V = np.zeros((len(rows), len(masks)))
            rowpos = {t: i for i, t in enumerate(rows)}

            for p in range(norb):
                src = (ML2 - ml2[p], MS2 - ms2[p])
                if src not in prev.block:
                    continue
                prows, P = prev.block[src]
                # creation operator a+_p : source determinants -> this block
                si, di, sg = [], [], []
                bit = 1 << p
                low = bit - 1
                for j, m in enumerate(prev.dets[src]):
                    if m & bit:
                        continue
                    si.append(j)
                    di.append(new.pos[m | bit][2])
                    sg.append(-1.0 if bin(m & low).count('1') % 2 else 1.0)
                if not si:
                    continue
                Q = np.zeros((P.shape[0], len(masks)))
                Q[:, di] = P[:, si] * np.array(sg)

                # coefficient matrix for this orbital
                C = np.zeros((len(rows), P.shape[0]))
                any_c = False
                for ri, ti in enumerate(rows):
                    L2, S2 = tls[ti]
                    for pj, pi in enumerate(prows):
                        k = inv.get(pi)
                        if k is None:
                            continue
                        c = deck.cfp[ti, k]
                        if abs(c) < 1e-12:
                            continue
                        Lb2, Sb2 = pls[pi]
                        a = cg(Lb2, ML2 - ml2[p], 2 * self.l, ml2[p], L2, ML2)
                        if abs(a) < TOL:
                            continue
                        b = cg(Sb2, MS2 - ms2[p], 1, ms2[p], S2, MS2)
                        if abs(b) < TOL:
                            continue
                        C[ri, pj] = c * a * b
                        any_c = True
                if any_c:
                    V += C @ Q

            self._verify(V, rows, n, ML2, MS2, len(masks))
            new.block[(ML2, MS2)] = (rows, V / np.sqrt(n))
        return new

    @staticmethod
    def _verify(V, rows, n, ML2, MS2, ndet):
        if len(rows) != ndet:
            raise AssertionError(
                'n=%d block (2ML,2MS)=(%d,%d): %d terms but %d determinants'
                % (n, ML2, MS2, len(rows), ndet))
        if not rows:
            return
        G = V @ V.T
        dev = np.abs(G - n * np.eye(len(rows))).max()
        if dev > 1e-8:
            raise AssertionError(
                'n=%d block (2ML,2MS)=(%d,%d): Gram deviates from %d*I by %.2e'
                % (n, ML2, MS2, n, dev))

    # -- convert to the structure LSShell.states() expects ------------------
    def _export(self, shell):
        terms = {}
        cols = [dict() for _ in shell.terms]
        for (ML2, MS2), (rows, V) in shell.block.items():
            masks = shell.dets[(ML2, MS2)]
            for ri, ti in enumerate(rows):
                st = {masks[c]: V[ri, c] for c in range(len(masks))
                      if abs(V[ri, c]) > TOL}
                cols[ti][(ML2, MS2)] = st
        for ti, lab in enumerate(shell.terms):
            L2, S2 = _term_ls(lab)
            # `v` is replaced by Cowan's alpha; keep the slot so that
            # LSShell.states() can format it unchanged
            terms.setdefault((L2, S2), []).append(
                (lab[2] if lab[2] else '', '', 1, cols[ti]))
        return terms

    def format_label(self, L2, S2, alpha, suf, J2):
        """(2S+1) L alpha ) 2J -- the form RCG writes at rcg11kd.f:7550,
        FORMAT (I1,A1,1H),I2,2A1)."""
        return '%d%s%s)%s' % (S2 + 1, LSYM[L2 // 2], alpha, lsjj._jlab(J2))

    def term_label(self, L2, S2, alpha, suf):
        return '%d%s%s' % (S2 + 1, LSYM[L2 // 2], alpha)

    # -- Racah's group labels of the CFP-built states ----------------------
    def racah_labels(self):
        """{(mult, Lsym, alpha): (W, U, sigma)} for an f shell.

        W and U are Racah's SO(7) and G2 irreps, read off the expectation
        values of the two Casimir operators in the CFP-built state; sigma is
        the root-mean-square spread, so sigma ~ 0 means the state really is
        an eigenstate of both, i.e. Cowan's alpha respects Racah's
        classification.  Evaluated on the highest-weight component, which is
        equivalent because the Casimirs commute with L and S.
        """
        from lsjj import casimir_apply, dot, name_irrep
        if self.l != 3:
            return {}
        out = {}
        for (L2, S2), built in self.terms.items():
            for alpha, _suf, _d, col in built:
                st = col[(L2, S2)]
                labs, spread = [], 0.0
                for ks, which in (((1, 3, 5), 'W'), ((1, 5), 'U')):
                    img = casimir_apply(st, self.orbs, self.l, ks)
                    e = dot(st, img)
                    var = dot(img, img) - e * e
                    spread = max(spread, np.sqrt(max(0.0, var)))
                    labs.append(name_irrep(e, which))
                out[(S2 + 1, LSYM[L2 // 2], alpha)] = (labs[0], labs[1],
                                                       spread)
        return out

    # -- seniority, for cross-checking against Cowan's SENIOR file ---------
    def seniority(self):
        """{(mult, Lsym, alpha): v} from < A+A > of the built states."""
        if self._seniority is None:
            from lsjj import (pair_matrix, seniority_from_eigenvalue)
            out = {}
            for (L2, S2), built in self.terms.items():
                for alpha, _suf, _d, col in built:
                    hw = [col[(L2, S2)]]
                    e = pair_matrix(hw, self.pair_ops)[0, 0]
                    out[(S2 + 1, LSYM[L2 // 2], alpha)] = \
                        seniority_from_eigenvalue(e, self.N, 2 * self.l + 1)
            self._seniority = out
        return self._seniority


# --------------------------------------------------------------------------
#  installation into lsjj
# --------------------------------------------------------------------------

def install():
    """Make lsjj build its LS side from Cowan's CFPs."""
    lsjj.set_ls_builder(CowanLSShell)


def uninstall():
    lsjj.set_ls_builder(LSShell)


if __name__ == '__main__':
    import sys
    args = sys.argv[1:]
    install()
    if not args or args[0] == 'terms':
        l = int(args[1]) if len(args) > 1 else 3
        N = int(args[2]) if len(args) > 2 else 3
        sh = CowanLSShell(l, N)
        sen = sh.seniority()
        print('%s^%d  Cowan terms (alpha as printed by RCG), with the '
              'seniority recovered from <A+A>:' % (LSYMB[l], N))
        line = []
        for (L2, S2), built in sorted(sh.terms.items()):
            for alpha, _s, _d, _c in built:
                key = (S2 + 1, LSYM[L2 // 2], alpha)
                line.append('%d%s%s(v=%d)' % (S2 + 1, LSYM[L2 // 2], alpha,
                                              sen[key]))
        print('  ' + '  '.join(line))
    elif args[0] == 'show':
        lsjj.report(int(args[1]), int(args[2]), show_matrix=True)
    elif args[0] == 'dump':
        lsjj.dump(int(args[1]), int(args[2]),
                  args[3] if len(args) > 3 else None)
    else:
        print(__doc__)
