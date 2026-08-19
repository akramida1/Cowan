"""Validation sweep for lsjj.py:  p^N, d^N, f^N.

    python validate.py [p|d|f|all] [cowan]

With the `cowan` argument the LS side is built from Cowan's CFP decks
(cowan_cfp.py) instead of by diagonalization, so the labels are RCG's.  The
term-list check then compares the seniority recovered from <A+A> for each of
Cowan's alpha labels against the published tables.
"""
import sys
import time
from collections import Counter
import lsjj

# Published term lists with seniorities (Racah / Cowan tables) for spot checks.
KNOWN = {
    (1, 2): '1S(0) 3P(2) 1D(2)',
    (1, 3): '4S(3) 2P(1) 2D(3)',
    (2, 2): '1S(0) 3P(2) 1D(2) 3F(2) 1G(2)',
    (2, 3): '2D(1) 4P(3) 2P(3) 2D(3) 2F(3) 4F(3) 2G(3) 2H(3)',
    (2, 4): '1S(0) 1D(2) 3P(2) 3F(2) 1G(2) 5D(4) 3P(4) 1S(4) 3D(4) '
            '1D(4) 3F(4) 1F(4) 3G(4) 1G(4) 3H(4) 1I(4)',
    (3, 2): '1S(0) 3P(2) 1D(2) 3F(2) 1G(2) 3H(2) 1I(2)',
    (3, 3): '4S(3) 2P(3) 2D(3) 2D(3) 4D(3) 2F(1) 2F(3) 4F(3) 2G(3) 2G(3) '
            '4G(3) 2H(3) 2H(3) 4I(3) 2I(3) 2K(3) 2L(3)',
}


def norm(s):
    return Counter(s.split()) if isinstance(s, str) else Counter(s)


def main():
    args = [a for a in sys.argv[1:]]
    cowan = 'cowan' in args
    args = [a for a in args if a != 'cowan']
    which = args[0] if args else 'all'
    tag = which + ('_cowan' if cowan else '')
    if cowan:
        import cowan_cfp
        cowan_cfp.install()
    lines = []
    out = lines.append
    npass = nfail = 0

    shells = {'p': (1, 6), 'd': (2, 10), 'f': (3, 14)}
    todo = [shells[c] for c in which] if which != 'all' \
        else [shells['p'], shells['d'], shells['f']]

    # ---- term lists first
    out('')
    out('LS side built from %s'
        % ('Cowan CFP decks' if cowan else 'seniority + SO(7)/G2 Casimirs'))
    out('')
    out('LS term lists with seniority, checked against published tables')
    out('-' * 74)
    for (l, N), expect in sorted(KNOWN.items()):
        if not any(l == ll for ll, _ in todo):
            continue
        if cowan:
            sen = lsjj.ls_shell(l, N).seniority()
            got = ['%d%s(%d)' % (k[0], k[1], v) for k, v in sen.items()]
        else:
            got = lsjj.term_list(l, N)
        ok = norm(got) == norm(expect)
        out('  %s^%-2d  %2d terms  %s' % ('spdf'[l], N, len(got),
                                          'OK' if ok else 'MISMATCH'))
        if not ok:
            out('        got      : %s' % ' '.join(sorted(got)))
            out('        expected : %s' % ' '.join(sorted(expect.split())))
            nfail += 1
        else:
            npass += 1

    # ---- the transformations
    for l, full in todo:
        for N in range(1, full + 1):
            t0 = time.time()
            try:
                ok = lsjj.report(l, N, show_matrix=False, out=out)
            except Exception as exc:
                out(' *** %s^%d raised %s: %s'
                    % ('spdf'[l], N, type(exc).__name__, exc))
                ok = False
            out(' %.2f s' % (time.time() - t0))
            npass += ok
            nfail += not ok

    out('')
    out('=' * 74)
    out(' %d passed, %d failed' % (npass, nfail))
    text = '\n'.join(lines)
    print(text)
    with open('validate_%s.log' % tag, 'w', newline='') as fh:
        fh.write(text + '\n')
    return 0 if nfail == 0 else 1


if __name__ == '__main__':
    sys.exit(main())
