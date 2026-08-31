#!/usr/bin/env python3
"""
LS -> jj transformation matrices for a single shell of equivalent electrons.

Route 2 of the three discussed:  build both coupled bases explicitly in the
m-scheme (Slater determinants of one-electron states) and take overlaps.
The only angular-momentum algebra involved is the Clebsch-Gordan coefficient;
there are no coefficients of fractional parentage, no 9j symbols, and no
phase conventions to reconcile between the two schemes -- both sides are
constructed by the same highest-weight-plus-lowering procedure, so
Condon-Shortley phases follow automatically.

Bases
-----
LS side   | l^N  a L S J M >            a = seniority v (plus a residual
                                        label where v L S is not complete)
jj side   | ( j-^N1 v1 J1, j+^N2 v2 J2 ) J M >     j- = l-1/2, j+ = l+1/2

Both span the same antisymmetric space, so for each J the two bases have the
same dimension and the transformation is square and orthogonal.

Conventions
-----------
Half-integer quantum numbers are carried doubled (2j, 2m) as integers.
Determinants are bit masks over an ordered one-electron orbital list;
| D > = a+_{o1} a+_{o2} ... a+_{oN} | 0 >  with  o1 < o2 < ... < oN.

The sign of each basis vector is fixed by requiring that its largest
component (ties broken by lowest determinant index) be positive.  This is a
convention, not physics; matching another code's tables requires matching
its sign convention as well.
"""

from math import factorial, sqrt
from itertools import combinations
import numpy as np

TOL = 1.0e-9
SIGN_TOL = 1.0e-12          # tie-break window for fix_sign; see there


# --------------------------------------------------------------------------
#  Clebsch-Gordan coefficient, all arguments doubled
# --------------------------------------------------------------------------

def cg(j1, m1, j2, m2, J, M):
    """< j1 m1 j2 m2 | J M >, arguments doubled (2j, 2m)."""
    if m1 + m2 != M:
        return 0.0
    if J < abs(j1 - j2) or J > j1 + j2:
        return 0.0
    if (j1 + j2 + J) % 2:
        return 0.0
    if abs(m1) > j1 or abs(m2) > j2 or abs(M) > J:
        return 0.0
    for x in (j1 + m1, j1 - m1, j2 + m2, j2 - m2, J + M, J - M):
        if x % 2:
            return 0.0

    f = factorial
    pre = sqrt((J + 1)
               * f((J + j1 - j2) // 2) * f((J - j1 + j2) // 2)
               * f((j1 + j2 - J) // 2) / f((j1 + j2 + J) // 2 + 1))
    pre *= sqrt(f((J + M) // 2) * f((J - M) // 2)
                * f((j1 - m1) // 2) * f((j1 + m1) // 2)
                * f((j2 - m2) // 2) * f((j2 + m2) // 2))

    total = 0.0
    for k in range(0, 200):
        args = ((j1 + j2 - J) // 2 - k,
                (j1 - m1) // 2 - k,
                (j2 + m2) // 2 - k,
                (J - j2 + m1) // 2 + k,
                (J - j1 - m2) // 2 + k)
        if args[0] < 0 or args[1] < 0 or args[2] < 0:
            if args[0] < 0:
                break
            continue
        if args[3] < 0 or args[4] < 0:
            continue
        denom = f(k)
        for a in args:
            denom *= f(a)
        total += (-1.0) ** k / denom
    return pre * total


# --------------------------------------------------------------------------
#  second quantization on bit-mask determinants
# --------------------------------------------------------------------------

def _annihilate(mask, p):
    bit = 1 << p
    if not (mask & bit):
        return 0, 0
    sign = -1 if bin(mask & (bit - 1)).count('1') % 2 else 1
    return sign, mask ^ bit


def _create(mask, p):
    bit = 1 << p
    if mask & bit:
        return 0, 0
    sign = -1 if bin(mask & (bit - 1)).count('1') % 2 else 1
    return sign, mask | bit


def apply_one_body(ops, state):
    """ops: list of (p, q, value) meaning  value * a+_p a_q .  state: dict."""
    out = {}
    for mask, c in state.items():
        for p, q, v in ops:
            s1, m1 = _annihilate(mask, q)
            if s1 == 0:
                continue
            s2, m2 = _create(m1, p)
            if s2 == 0:
                continue
            out[m2] = out.get(m2, 0.0) + c * v * s1 * s2
    return {k: v for k, v in out.items() if abs(v) > TOL}


def apply_pair_annihilator(pairs, state):
    """pairs: list of (p, q, value) meaning  value * a_p a_q  (p != q)."""
    out = {}
    for mask, c in state.items():
        for p, q, v in pairs:
            s1, m1 = _annihilate(mask, q)
            if s1 == 0:
                continue
            s2, m2 = _annihilate(m1, p)
            if s2 == 0:
                continue
            out[m2] = out.get(m2, 0.0) + c * v * s1 * s2
    return {k: v for k, v in out.items() if abs(v) > TOL}


def dot(a, b):
    if len(a) > len(b):
        a, b = b, a
    return sum(c * b.get(k, 0.0) for k, c in a.items())


def normalize(state):
    n = sqrt(dot(state, state))
    return {k: v / n for k, v in state.items()}


def fix_sign(state):
    """Deterministic overall sign: largest component positive.

    The tie-break is deliberately tolerant.  Two components equal in exact
    arithmetic can differ in the last bit from one LAPACK build to the next;
    a strict max() then crowns a different component on a different machine
    and flips the sign of the whole vector -- and of every state derived
    from it down the lowering chain, which fixes no signs of its own.
    Comparing magnitudes only to SIGN_TOL, and settling the tie on the
    basis key, makes the convention depend on our own ordering rather than
    on the library's last bit.  With an exact tie this picks the same
    component the old code did, so nothing changes except at a near-tie.
    """
    if not state:
        return state
    top = max(abs(v) for v in state.values())
    key = min(k for k, v in state.items() if abs(v) >= top - SIGN_TOL)
    if state[key] < 0:
        return {k: -v for k, v in state.items()}
    return state


# --------------------------------------------------------------------------
#  generic machinery: build coupled states in a set of orbitals
# --------------------------------------------------------------------------

class OrbitalSet:
    """An ordered list of one-electron orbitals with raising/lowering ops.

    Each orbital carries a tuple of doubled projection quantum numbers.
    `ladders` is a list of (raise_ops, lower_ops, proj_index) describing each
    independent angular momentum whose ladder we need.
    """

    def __init__(self, orbitals, projections):
        self.orbitals = orbitals
        self.n = len(orbitals)
        self.index = {o: i for i, o in enumerate(orbitals)}
        self.projections = projections   # list of lists, doubled, per orbital

    def dets(self, N):
        out = []
        for combo in combinations(range(self.n), N):
            mask = 0
            for i in combo:
                mask |= 1 << i
            out.append(mask)
        return out

    def proj(self, mask, which):
        tot = 0
        i = 0
        m = mask
        while m:
            if m & 1:
                tot += self.projections[which][i]
            m >>= 1
            i += 1
        return tot


def build_ladder_ops(orbs, coord, jvals):
    """Raising/lowering one-body ops for one angular momentum.

    `orbs` is the ordered list of orbital label tuples, all entries doubled;
    `coord` selects which entry of the tuple is the projection being laddered.
    Every other entry must match, or the operator would connect orbitals that
    differ in an unrelated quantum number.
    """
    up, dn = [], []
    for q, oq in enumerate(orbs):
        j, m = jvals[q], oq[coord]
        for p, op in enumerate(orbs):
            if jvals[p] != j:
                continue
            if any(op[k] != oq[k] for k in range(len(oq)) if k != coord):
                continue
            if op[coord] == m + 2:
                up.append((p, q, sqrt((j * (j + 2) - m * (m + 2)) / 4.0)))
            if op[coord] == m - 2:
                dn.append((p, q, sqrt((j * (j + 2) - m * (m - 2)) / 4.0)))
    return up, dn


def highest_weight_space(dets, raise_op_list, proj_getter, targets):
    """Determinants with the given projections, annihilated by every raiser.

    Returns an orthonormal basis (list of dicts) of the highest-weight space.
    """
    sub = [d for d in dets
           if all(proj_getter(d, i) == t for i, t in enumerate(targets))]
    if not sub:
        return []
    idx = {d: i for i, d in enumerate(sub)}
    # matrix of each raising operator restricted to this block -> null space
    rows = []
    for ops in raise_op_list:
        for d in sub:
            res = apply_one_body(ops, {d: 1.0})
            rows.append((d, res))
    # collect the image space
    img_keys = sorted({k for _, res in rows for k in res})
    if img_keys:
        A = np.zeros((len(img_keys), len(sub)))
        kidx = {k: i for i, k in enumerate(img_keys)}
        for d, res in rows:
            for k, v in res.items():
                A[kidx[k], idx[d]] += v
        # null space of A
        u, s, vt = np.linalg.svd(A)
        rank = int((s > 1e-8).sum())
        null = vt[rank:]
    else:
        null = np.eye(len(sub))
    return [{sub[i]: v for i, v in enumerate(vec) if abs(v) > TOL}
            for vec in null]


def seniority_from_eigenvalue(e, N, Omega):
    """Seniority from the quasispin eigenvalue.

    With the pair operator normalized as built here (verified numerically
    against p^2 1S = 3, p^3 2P = 2, d^3 2D = 0 and 4),

        < A+ A >  =  Q(Q+1) - Qz(Qz-1),   Qz = (N-Omega)/2,  Q = (Omega-v)/2

    so v follows by inverting the quadratic.  Omega = 2l+1 for an LS shell
    l^N, and j+1/2 for a single j-shell.
    """
    qz = (N - Omega) / 2.0
    X = e + qz * (qz - 1.0)
    Q = (-1.0 + sqrt(max(0.0, 1.0 + 4.0 * X))) / 2.0
    v = Omega - 2.0 * Q
    vi = int(round(v))
    if abs(v - vi) > 1e-6 or vi < 0 or vi > N or (N - vi) % 2:
        raise ValueError('non-integer seniority %.6f (e=%.6f N=%d Omega=%d)'
                         % (v, e, N, Omega))
    return vi


def pair_matrix(space, pair_ops):
    """Matrix of A+A inside `space`."""
    k = len(space)
    images = [apply_pair_annihilator(pair_ops, s) for s in space]
    G = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            G[i, j] = dot(images[i], images[j])
    return G


def seniority_split(space, pair_ops, N, Omega, orbs=None, l=None):
    """Classify `space` by seniority, and for an f shell also by the
    Casimir operators of SO(7) and G2 (Racah's W and U labels).

    Returns a list of (v, extra_labels, degen, state) where `degen` is the
    size of the group that still shares every label -- the split inside such
    a group is an arbitrary rotation, not physics.
    """
    if not space:
        return []
    out = []
    for e, vecs in split_by_operator(space, pair_matrix(space, pair_ops)):
        v = seniority_from_eigenvalue(e, N, Omega)
        # Racah's group labels are only needed where seniority is not enough,
        # which for a single shell means l = 3.
        if len(vecs) > 1 and l == 3 and orbs is not None:
            level2 = []
            for w7, g7 in split_by_operator(
                    vecs, casimir_matrix(vecs, orbs, l, (1, 3, 5))):
                for w2, g2 in split_by_operator(
                        g7, casimir_matrix(g7, orbs, l, (1, 5))):
                    level2.append(((w7, w2), g2))
        else:
            level2 = [((), vecs)]
        for labels, grp in level2:
            for st in grp:
                out.append((v, labels, len(grp),
                            fix_sign(normalize(st))))
    return out


def threej(j1, m1, j2, m2, j3, m3):
    """3j symbol, arguments doubled."""
    if m1 + m2 + m3 != 0:
        return 0.0
    ph = (-1.0) ** ((j1 - j2 - m3) // 2)
    return ph * cg(j1, m1, j2, m2, j3, -m3) / sqrt(j3 + 1)


def unit_tensor_ops(orbs, l, k, q):
    """One-body ops for the unit tensor U^k_q acting on the orbital part.

    < l m' | u^k_q | l m > = (-1)^(l-m') * ( l  k  l ; -m'  q  m )
    with the reduced matrix element < l || u^k || l > defined as 1.
    Spin is a spectator.
    """
    ops = []
    for i, (ml2, ms2) in enumerate(orbs):
        for p, (ml2p, ms2p) in enumerate(orbs):
            if ms2p != ms2 or ml2p != ml2 + 2 * q:
                continue
            v = ((-1.0) ** (l - ml2p // 2)
                 * threej(2 * l, -ml2p, 2 * k, 2 * q, 2 * l, ml2))
            if abs(v) > TOL:
                ops.append((p, i, v))
    return ops


def casimir_apply(state, orbs, l, ks):
    """sum_k (2k+1) (U^k . U^k) applied to one state."""
    acc = {}
    for k in ks:
        tensors = {q: unit_tensor_ops(orbs, l, k, q)
                   for q in range(-k, k + 1)}
        for q in range(-k, k + 1):
            inner = apply_one_body(tensors[-q], state)
            if not inner:
                continue
            outer = apply_one_body(tensors[q], inner)
            ph = (2 * k + 1) * (-1.0) ** q
            for key, val in outer.items():
                acc[key] = acc.get(key, 0.0) + ph * val
    return {key: val for key, val in acc.items() if abs(val) > TOL}


def casimir_matrix(space, orbs, l, ks):
    """Matrix of  sum_k (2k+1) (U^k . U^k)  inside `space`.

    k = 1,3,5 gives the Casimir of SO(7); k = 1,5 that of G2 (up to the
    normalization and the one-body term, neither of which affects the
    eigenvectors or the ordering of the eigenvalues within a fixed N).
    """
    n = len(space)
    images = [casimir_apply(s, orbs, l, ks) for s in space]
    C = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            C[i, j] = dot(space[i], images[j])
    return 0.5 * (C + C.T)


#  The two Casimir matrices built above come out as exactly
#      4 * G(G2)   and   5 * G(R7)
#  in Racah's normalization, verified on f^3 where the v=3 doublets give
#  SO(7) = 9.0 for all (W = (210), G(R7) = 1.8) and G2 = 4.0, 4.6667, 7.0
#  for U = (11), (20), (21) (G(G2) = 1, 7/6, 7/4).
SO7_SCALE, G2_SCALE = 5.0, 4.0


def _so7_irreps(wmax=3):
    for w1 in range(wmax + 1):
        for w2 in range(w1 + 1):
            for w3 in range(w2 + 1):
                yield (w1, w2, w3), 0.1 * (w1 * (w1 + 5) + w2 * (w2 + 3)
                                           + w3 * (w3 + 1))


def _g2_irreps(umax=4):
    for u1 in range(umax + 1):
        for u2 in range(u1 + 1):
            yield (u1, u2), (u1 * u1 + u2 * u2 + u1 * u2
                             + 5 * u1 + 4 * u2) / 12.0


def name_irrep(eig, which):
    """Turn a Casimir eigenvalue into Racah's (W) or (U) label."""
    table = (list(_so7_irreps()) if which == 'W' else list(_g2_irreps()))
    scale = SO7_SCALE if which == 'W' else G2_SCALE
    target = eig / scale
    hits = [lab for lab, val in table if abs(val - target) < 1e-6]
    if not hits:
        return '?%.4f' % target
    if len(hits) > 1:
        return '|'.join(''.join(str(x) for x in h) for h in hits)
    return ''.join(str(x) for x in hits[0])


def split_by_operator(space, matrix):
    """Block a space by the eigenvalues of a commuting operator.

    Returns [(eigenvalue, [states])] with the states orthonormal.
    """
    w, V = np.linalg.eigh(matrix)
    order = np.argsort(w)
    w, V = w[order], V[:, order]
    groups, out = [], []
    for i, wi in enumerate(w):
        if groups and abs(wi - w[groups[-1][0]]) < 1e-6:
            groups[-1].append(i)
        else:
            groups.append([i])
    for grp in groups:
        vecs = []
        for col in grp:
            st = {}
            for i, ci in enumerate(V[:, col]):
                if abs(ci) < TOL:
                    continue
                for key, val in space[i].items():
                    st[key] = st.get(key, 0.0) + ci * val
            vecs.append({k: v for k, v in st.items() if abs(v) > TOL})
        out.append((w[grp[0]], vecs))
    return out


def lower_chain(state, lower_ops, j2, m2):
    """Apply the lowering operator once, with Condon-Shortley normalization."""
    res = apply_one_body(lower_ops, state)
    norm = sqrt((j2 * (j2 + 2) - m2 * (m2 - 2)) / 4.0)
    return {k: v / norm for k, v in res.items()}


# --------------------------------------------------------------------------
#  LS side:  l^N  a L S J M
# --------------------------------------------------------------------------

LSYM = 'SPDFGHIKLMNOQRTUVWXYZABC'   # LSYMB of rcg11kd.f, line 856


class LSShell:
    def __init__(self, l, N):
        self.l, self.N = l, N
        # orbital labels are (2*ml, 2*ms), both doubled
        orbs = [(ml2, ms2) for ml2 in range(2 * l, -2 * l - 1, -2)
                for ms2 in (1, -1)]
        self.orbs = orbs
        self.ml2 = [o[0] for o in orbs]
        self.ms2 = [o[1] for o in orbs]
        self.lv2 = [2 * l] * len(orbs)
        self.sv2 = [1] * len(orbs)
        self.set = OrbitalSet(orbs, [self.ml2, self.ms2])
        self.Lp, self.Lm = build_ladder_ops(orbs, 0, self.lv2)
        self.Sp, self.Sm = build_ladder_ops(orbs, 1, self.sv2)
        self.dets = self.set.dets(N)
        self.pair_ops = self._pair_ops()
        self.terms = self._build_terms()

    def _pair_ops(self):
        l = self.l
        ops = []
        for i, (ml2, ms2) in enumerate(self.orbs):
            jbar = self.set.index[(-ml2, -ms2)]
            ph = (-1) ** (l - ml2 // 2) * (1 if ms2 == 1 else -1)
            ops.append((jbar, i, 0.5 * ph))
        return ops

    def _build_terms(self):
        """dict (L2, S2) -> list of (v, degen, {(ML2,MS2): state})."""
        terms = {}
        maxL2 = 2 * self.l * self.N
        for L2 in range(0, maxL2 + 1):
            for S2 in range(self.N % 2, self.N + 1, 2):
                hw = highest_weight_space(
                    self.dets, [self.Lp, self.Sp],
                    lambda d, i: self.set.proj(d, i), (L2, S2))
                if not hw:
                    continue
                split = seniority_split(hw, self.pair_ops, self.N,
                                        2 * self.l + 1, self.orbs, self.l)
                built = []
                for v, labels, degen, st in split:
                    if labels:
                        # W and U separated by '/': a Casimir eigenvalue can
                        # be shared by two irreps (SO(7) (220) and (300), for
                        # one), in which case name_irrep returns '220|300'
                        # and a fixed-width split would be wrong
                        suf = ' %s/%s' % (name_irrep(labels[0], 'W'),
                                          name_irrep(labels[1], 'U'))
                    else:
                        suf = ''
                    built.append((v, suf, degen, self._ladder(st, L2, S2)))
                if built:
                    terms[(L2, S2)] = built
        return terms

    def _ladder(self, hw, L2, S2):
        """All (ML2, MS2) components of a term, from the highest weight."""
        col = {(L2, S2): hw}
        for ML2 in range(L2, -L2 - 1, -2):
            if ML2 != L2:
                col[(ML2, S2)] = lower_chain(col[(ML2 + 2, S2)],
                                             self.Lm, L2, ML2 + 2)
            for MS2 in range(S2 - 2, -S2 - 1, -2):
                col[(ML2, MS2)] = lower_chain(col[(ML2, MS2 + 2)],
                                              self.Sm, S2, MS2 + 2)
        return col

    def format_label(self, L2, S2, v, suf, J2):
        """(2S+1) L (seniority [+ Racah W U]) J."""
        return '%d%s(%s%s)%s' % (S2 + 1, LSYM[L2 // 2], v, suf, _jlab(J2))

    def term_label(self, L2, S2, v, suf):
        return '%d%s(%s%s)' % (S2 + 1, LSYM[L2 // 2], v, suf)

    def states(self, M2):
        """List of (label, J2, state) at fixed total M2."""
        out = []
        for (L2, S2), built in sorted(self.terms.items()):
            for J2 in range(abs(L2 - S2), L2 + S2 + 1, 2):
                if abs(M2) > J2:
                    continue
                for v, suf, degen, col in built:
                    st = {}
                    for ML2 in range(-L2, L2 + 1, 2):
                        MS2 = M2 - ML2
                        if abs(MS2) > S2:
                            continue
                        c = cg(L2, ML2, S2, MS2, J2, M2)
                        if abs(c) < TOL:
                            continue
                        for k, val in col[(ML2, MS2)].items():
                            st[k] = st.get(k, 0.0) + c * val
                    st = {k: v2 for k, v2 in st.items() if abs(v2) > TOL}
                    if not st:
                        continue
                    lab = '%-10s' % self.format_label(L2, S2, v, suf, J2)
                    # no sign fixing here: the phase is already determined by
                    # the highest-weight state and the lowering chain, and
                    # re-fixing it at each M would make it M dependent
                    out.append((lab, J2, degen, normalize(st)))
        return out


def _jlab(J2):
    return str(J2 // 2) if J2 % 2 == 0 else '%d/2' % J2


# Which class builds the LS side.  The default classifies the states by
# seniority and, for f, by the Casimir operators of SO(7) and G2.  Importing
# cowan_cfp and calling its install() swaps in a builder that reads Cowan's
# CFP decks instead, so that the alpha labels and the signs are RCG's.
LS_BUILDER = LSShell
_LS_CACHE = {}


def set_ls_builder(cls):
    global LS_BUILDER
    LS_BUILDER = cls
    _LS_CACHE.clear()


def ls_shell(l, N):
    key = (LS_BUILDER, l, N)
    if key not in _LS_CACHE:
        _LS_CACHE[key] = LS_BUILDER(l, N)
    return _LS_CACHE[key]


# --------------------------------------------------------------------------
#  jj side:  ( j-^N1 v1 J1 , j+^N2 v2 J2 ) J M
# --------------------------------------------------------------------------

class JSubshell:
    """A single j-shell, j^n, classified by seniority and J."""

    def __init__(self, j2, n):
        self.j2, self.n = j2, n
        self.mvals = list(range(j2, -j2 - 1, -2))
        self.jv2 = [j2] * len(self.mvals)
        self.set = OrbitalSet([(m,) for m in self.mvals], [self.mvals])
        self.Jp, self.Jm = build_ladder_ops(self.set.orbitals, 0, self.jv2)
        self.dets = self.set.dets(n)
        self.pair_ops = self._pair_ops()
        self.levels = self._build()

    def _pair_ops(self):
        ops = []
        for i, m in enumerate(self.mvals):
            ibar = self.mvals.index(-m)
            ph = (-1) ** ((self.j2 - m) // 2)
            ops.append((ibar, i, 0.5 * ph))
        return ops

    def _build(self):
        """list of (J2, v, degen, {M2: state})."""
        out = []
        maxJ2 = self.j2 * self.n
        for J2 in range(0, maxJ2 + 1):
            hw = highest_weight_space(self.dets, [self.Jp],
                                      lambda d, i: self.set.proj(d, i), (J2,))
            if not hw:
                continue
            for v, _lab, degen, st in seniority_split(
                    hw, self.pair_ops, self.n, (self.j2 + 1) // 2):
                col = {J2: st}
                for M2 in range(J2 - 2, -J2 - 1, -2):
                    col[M2] = lower_chain(col[M2 + 2], self.Jm, J2, M2 + 2)
                out.append((J2, v, degen, col))
        return out


class JJShell:
    def __init__(self, l, N):
        self.l, self.N = l, N
        self.jm2, self.jp2 = 2 * l - 1, 2 * l + 1
        self.nm_max = self.jm2 + 1
        self.np_max = self.jp2 + 1
        # global jj orbital order: all j- first (m descending), then all j+
        self.orbs = ([(self.jm2, m) for m in range(self.jm2, -self.jm2 - 1, -2)]
                     if l > 0 else [])
        self.orbs += [(self.jp2, m)
                      for m in range(self.jp2, -self.jp2 - 1, -2)]
        self.index = {o: i for i, o in enumerate(self.orbs)}
        self.nminus = len([o for o in self.orbs if o[0] == self.jm2]) \
            if l > 0 else 0

    def states(self, M2):
        """List of (label, J2, degen, state-in-jj-determinants)."""
        out = []
        for n1 in range(0, min(self.nm_max, self.N) + 1):
            n2 = self.N - n1
            if n2 < 0 or n2 > self.np_max:
                continue
            if self.l == 0 and n1:
                continue
            sub1 = JSubshell(self.jm2, n1) if self.jm2 > 0 else None
            sub2 = JSubshell(self.jp2, n2)
            lv1 = sub1.levels if sub1 else [(0, 0, 1, {0: {0: 1.0}})]
            for J1, v1, d1, col1 in lv1:
                for J2s, v2, d2, col2 in sub2.levels:
                    for J2 in range(abs(J1 - J2s), J1 + J2s + 1, 2):
                        if abs(M2) > J2:
                            continue
                        st = {}
                        for M1 in range(-J1, J1 + 1, 2):
                            M2s = M2 - M1
                            if abs(M2s) > J2s:
                                continue
                            c = cg(J1, M1, J2s, M2s, J2, M2)
                            if abs(c) < TOL:
                                continue
                            for k1, a in col1[M1].items():
                                for k2, b in col2[M2s].items():
                                    key = self._merge(k1, k2, n1)
                                    st[key] = st.get(key, 0.0) + c * a * b
                        st = {k: v for k, v in st.items() if abs(v) > TOL}
                        if not st:
                            continue
                        # j^n v J for each subshell, then the total J.
                        # The seniority is part of the label, not decoration:
                        # 7/2^4 J=2, for one, occurs with both v=2 and v=4.
                        lab = '(%s^%d %d %s, %s^%d %d %s)%s' % (
                            _jlab(self.jm2), n1, v1, _jlab(J1),
                            _jlab(self.jp2), n2, v2, _jlab(J2s), _jlab(J2))
                        out.append(('%-30s' % lab, J2, max(d1, d2),
                                    normalize(st)))
        return out

    def _merge(self, k1, k2, n1):
        """j- mask (local) and j+ mask (local) -> global jj determinant mask.

        The global order puts every j- orbital before every j+ orbital, and
        each block is already in ascending local order, so the concatenation
        is already sorted: the relative sign is +1.
        """
        return k1 | (k2 << self.nminus)


# --------------------------------------------------------------------------
#  jj determinants -> ml,ms determinants
# --------------------------------------------------------------------------

def spin_orbit_cg(l, j2, mj2, ms2):
    """< l, mj-ms ; 1/2, ms | j mj >, closed form."""
    ml2 = mj2 - ms2
    if abs(ml2) > 2 * l:
        return 0.0
    mj = mj2 / 2.0
    if j2 == 2 * l + 1:
        if ms2 == 1:
            return sqrt((l + mj + 0.5) / (2 * l + 1))
        return sqrt((l - mj + 0.5) / (2 * l + 1))
    else:
        if ms2 == 1:
            return -sqrt((l - mj + 0.5) / (2 * l + 1))
        return sqrt((l + mj + 0.5) / (2 * l + 1))


def jj_det_to_ls(mask, jj_orbs, ls_index, l):
    """Expand one jj determinant into ml,ms determinants."""
    occ = [jj_orbs[i] for i in range(len(jj_orbs)) if mask & (1 << i)]
    result = {(): 1.0}
    for (j2, mj2) in occ:
        nxt = {}
        for ms2 in (1, -1):
            c = spin_orbit_cg(l, j2, mj2, ms2)
            if abs(c) < TOL:
                continue
            p = ls_index[(mj2 - ms2, ms2)]
            for seq, coef in result.items():
                nxt[seq + (p,)] = nxt.get(seq + (p,), 0.0) + coef * c
        result = nxt
    out = {}
    for seq, coef in result.items():
        if len(set(seq)) != len(seq):
            continue
        perm = sorted(range(len(seq)), key=lambda i: seq[i])
        # parity of the permutation
        sign, seen = 1, [False] * len(seq)
        for i in range(len(seq)):
            if seen[i]:
                continue
            j, ln = i, 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                ln += 1
            if ln % 2 == 0:
                sign = -sign
        m = 0
        for p in seq:
            m |= 1 << p
        out[m] = out.get(m, 0.0) + coef * sign
    return {k: v for k, v in out.items() if abs(v) > TOL}


# --------------------------------------------------------------------------
#  driver
# --------------------------------------------------------------------------

def transformation(l, N, verbose=False):
    M2 = 0 if N % 2 == 0 else 1

    ls = ls_shell(l, N)
    ls_states = ls.states(M2)

    jj = JJShell(l, N)
    jj_states = jj.states(M2)

    # convert jj states into the ml,ms determinant space
    cache = {}
    jj_conv = []
    for lab, J2, deg, st in jj_states:
        acc = {}
        for mask, c in st.items():
            if mask not in cache:
                cache[mask] = jj_det_to_ls(mask, jj.orbs, ls.set.index, l)
            for k, v in cache[mask].items():
                acc[k] = acc.get(k, 0.0) + c * v
        jj_conv.append((lab, J2, deg, {k: v for k, v in acc.items()
                                       if abs(v) > TOL}))

    return ls_states, jj_conv


def report(l, N, show_matrix=False, out=print):
    name = 'spdfgh'[l]
    ls_states, jj_states = transformation(l, N)
    M2 = 0 if N % 2 == 0 else 1

    out('')
    out('=' * 74)
    out(' %s^%d      m-scheme dimension at 2M=%d : %d   (full shell: %d)'
        % (name, N, M2, len(ls_states),
           len(list(combinations(range(2 * (2 * l + 1)), N)))))
    out('=' * 74)

    ok = True
    # ---- per-J dimensions
    jls, jjj = {}, {}
    for _, J2, _, _ in ls_states:
        jls[J2] = jls.get(J2, 0) + 1
    for _, J2, _, _ in jj_states:
        jjj[J2] = jjj.get(J2, 0) + 1
    if jls != jjj:
        ok = False
        out(' *** DIMENSION MISMATCH')
        out('     LS: %s' % sorted(jls.items()))
        out('     jj: %s' % sorted(jjj.items()))
    else:
        out(' per-J dimensions agree: ' +
            ' '.join('J=%s:%d' % (_jlab(k), v) for k, v in sorted(jls.items())))

    # ---- orthonormality of each basis
    for tag, states in (('LS', ls_states), ('jj', jj_states)):
        n = len(states)
        G = np.zeros((n, n))
        for i in range(n):
            for j in range(i, n):
                G[i, j] = G[j, i] = dot(states[i][3], states[j][3])
        err = np.abs(G - np.eye(n)).max() if n else 0.0
        if err > 1e-8:
            ok = False
            out(' *** %s basis not orthonormal, max dev %.2e' % (tag, err))

    # ---- the transformation, block by J
    maxdev = 0.0
    for J2 in sorted(jls):
        Lst = [s for s in ls_states if s[1] == J2]
        Jst = [s for s in jj_states if s[1] == J2]
        T = np.zeros((len(Lst), len(Jst)))
        for i, a in enumerate(Lst):
            for j, b in enumerate(Jst):
                T[i, j] = dot(a[3], b[3])
        dev = np.abs(T @ T.T - np.eye(len(Lst))).max()
        maxdev = max(maxdev, dev)
        if show_matrix:
            out('')
            out(' J = %s' % _jlab(J2))
            for j, b in enumerate(Jst):
                out('   %2d = %s' % (j + 1, b[0].strip()))
            out(' ' * 12 + ''.join('%9d' % (j + 1) for j in range(len(Jst))))
            for i, a in enumerate(Lst):
                out(' %-11s' % a[0].strip()
                    + ''.join('%9.5f' % T[i, j] for j in range(len(Jst))))
    if maxdev > 1e-8:
        ok = False
        out(' *** transformation not orthogonal, max deviation %.2e' % maxdev)
    else:
        out(' transformation orthogonal to %.1e in every J block' % maxdev)

    # ---- M independence
    if N % 2 == 0 and N > 0:
        dev = check_m_independence(l, N)
        if dev is not None:
            if dev > 1e-8:
                ok = False
                out(' *** M dependence detected, %.2e' % dev)
            else:
                out(' independent of M to %.1e' % dev)

    # ---- unresolved labels
    amb = sorted({s[0].strip() for s in ls_states if s[2] > 1})
    if amb:
        out(' NOTE: %d LS state(s) share every label (v and, for f, Racah W '
            'and U);' % len(amb))
        out('       the split inside such a group is an arbitrary rotation:')
        out('       ' + ', '.join(amb))
    amb = sorted({s[0].strip() for s in jj_states if s[2] > 1})
    if amb:
        out(' NOTE: jj states with degenerate (v,J) in a subshell: '
            + ', '.join(amb))

    out(' %s' % ('PASS' if ok else 'FAIL'))
    return ok


def check_m_independence(l, N):
    """Recompute the J-block transformation at a different M and compare."""
    M2a, M2b = 0, 2
    ls = ls_shell(l, N)
    jj = JJShell(l, N)
    cache = {}

    def build(M2):
        lss = ls.states(M2)
        jjs = []
        for lab, J2, deg, st in jj.states(M2):
            acc = {}
            for mask, c in st.items():
                if mask not in cache:
                    cache[mask] = jj_det_to_ls(mask, jj.orbs, ls.set.index, l)
                for k, v in cache[mask].items():
                    acc[k] = acc.get(k, 0.0) + c * v
            jjs.append((lab, J2, deg, acc))
        return lss, jjs

    la, ja = build(M2a)
    lb, jb = build(M2b)
    if not lb:
        return None
    dev = 0.0
    for J2 in sorted({s[1] for s in lb}):
        La = [s for s in la if s[1] == J2]
        Ja = [s for s in ja if s[1] == J2]
        Lb = [s for s in lb if s[1] == J2]
        Jb = [s for s in jb if s[1] == J2]
        if not (La and Lb):
            continue
        Ta = np.array([[dot(a[3], b[3]) for b in Ja] for a in La])
        Tb = np.array([[dot(a[3], b[3]) for b in Jb] for a in Lb])
        dev = max(dev, np.abs(Ta - Tb).max())
    return dev


def term_list(l, N):
    """LS terms with seniority, for checking against published tables."""
    ls = ls_shell(l, N)
    out = []
    for (L2, S2), built in sorted(ls.terms.items()):
        for v, suf, degen, _ in built:
            # the Racah suffix is deliberately dropped: the published
            # tables this is checked against carry only the seniority
            out.append(ls.term_label(L2, S2, v, ''))
    return out


def dump(l, N, path=None):
    """Write the transformation as plain text, one J block at a time.

    Format, suitable for reading back from Fortran:
        l  N
        J2  nrows  ncols
        <row labels>
        <col labels>
        matrix, one row per line
    """
    ls_states, jj_states = transformation(l, N)
    lines = ['%3d %3d' % (l, N)]
    for J2 in sorted({s[1] for s in ls_states}):
        L = [s for s in ls_states if s[1] == J2]
        Jt = [s for s in jj_states if s[1] == J2]
        lines.append('%3d %4d %4d' % (J2, len(L), len(Jt)))
        for s in L:
            lines.append('  %s' % s[0].strip())
        for s in Jt:
            lines.append('  %s' % s[0].strip())
        for a in L:
            row = [dot(a[3], b[3]) for b in Jt]
            lines.append(''.join('%14.9f' % x for x in row))
    text = '\n'.join(lines) + '\n'
    if path:
        with open(path, 'w', newline='') as fh:
            fh.write(text)
    else:
        print(text)


if __name__ == '__main__':
    import sys
    # go through the imported module, not through this __main__ copy of it,
    # so that cowan_cfp.install() switches the builder these calls actually see
    import lsjj as L

    args = sys.argv[1:]
    if 'cowan' in args:
        args.remove('cowan')
        import cowan_cfp
        cowan_cfp.install()
    if not args or args[0] == 'demo':
        L.report(1, 2, show_matrix=True)
        L.report(2, 2, show_matrix=True)
        L.report(1, 3, show_matrix=True)
    elif args[0] == 'show':
        L.report(int(args[1]), int(args[2]), show_matrix=True)
    elif args[0] == 'dump':
        L.dump(int(args[1]), int(args[2]),
               args[3] if len(args) > 3 else None)
    elif args[0] == 'terms':
        print(' '.join(L.term_list(int(args[1]), int(args[2]))))
    else:
        print(__doc__)
