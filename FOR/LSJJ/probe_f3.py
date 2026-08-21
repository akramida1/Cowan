"""Check that the Casimir split of f^3 reproduces Racah's pattern:
every v=3 doublet has a single SO(7) eigenvalue (W=(210)) and is
separated only by G2 (U)."""
import os
import sys

# The modules this checks -- lsjj, cowan_cfp, recur -- moved to CODE/lsjj/,
# where make_cfp can find them in an installed tree.  This directory keeps
# the record of the derivation and its validation, so the scripts here reach
# across to the one copy rather than carry a second.
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                '..', '..', 'CODE', 'lsjj'))
import numpy as np
from lsjj import (LSShell, highest_weight_space, pair_matrix,
                  split_by_operator, casimir_matrix,
                  seniority_from_eigenvalue, LSYM)

l, N = 3, 3
sh = LSShell(l, N)
print('%-6s %-4s %-3s  %-22s %s' % ('term', 'dim', 'v', 'SO(7) Casimir',
                                    'G2 Casimir'))
for L2 in range(2 * l * N + 1):
    for S2 in range(N % 2, N + 1, 2):
        hw = highest_weight_space(sh.dets, [sh.Lp, sh.Sp],
                                  lambda d, i: sh.set.proj(d, i), (L2, S2))
        if len(hw) < 2:
            continue
        for e, vecs in split_by_operator(hw, pair_matrix(hw, sh.pair_ops)):
            v = seniority_from_eigenvalue(e, N, 2 * l + 1)
            if len(vecs) < 2:
                continue
            w7 = np.linalg.eigvalsh(casimir_matrix(vecs, sh.orbs, l, (1, 3, 5)))
            w2 = np.linalg.eigvalsh(casimir_matrix(vecs, sh.orbs, l, (1, 5)))
            print('%d%-5s %-4d %-3d  %-22s %s'
                  % (S2 + 1, LSYM[L2 // 2], len(vecs), v,
                     ' '.join('%.4f' % x for x in w7),
                     ' '.join('%.4f' % x for x in w2)))
