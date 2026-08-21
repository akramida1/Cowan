"""Every numeric literal in the Fortran is accounted for.

Sections 4, 5 and 6 of the constants work substituted the constants that
an inventory of *known spellings* had found.  An inventory built that way
can only find what someone thought to look for, and Section 10 found
three more sites by re-reading every literal in the package one by one.
One of them, the ``2066.D0`` in the autoionization rate, was wrong by
5.2e-04 -- the largest constant error in the package, and the only one
that reached printed output.

This test is what stops that from happening again.  It reads every
non-comment real literal with three or more significant digits out of
``FOR/`` and fails on any that is not in the inventory below.  A new
constant cannot be added to the Fortran without either going through
``fundconst.f`` or being written down here with a reason.

Scope, and why it is real literals only
---------------------------------------
An *integer* literal in this package is a statement label, an array
bound, a loop limit or a FORMAT width; there are hundreds and not one is
physics.  A physical constant is always written with a decimal point or
an exponent -- ``2066.D0`` and ``109735.3D0`` both are -- so restricting
the scan to real literals loses nothing and keeps the inventory small
enough to be read.  Comment lines, FORMAT statements, character strings
and everything past column 72 are skipped, which is what fixed-form
Fortran means by "code".

Pure mathematical constants are listed here but are deliberately *not*
candidates for ``fundconst.f``.  They cannot go stale.

There is no exception list
--------------------------
Section 13 folded in the last four inline constants -- the electron
g-factor in its three spellings and the Maxwellian rate-coefficient
prefactor -- and with them the holding-pen group that used to sit here.
What this test now asserts is the stronger statement, and it is the
author's standing rule: **``FOR/`` contains no physical constant outside
``fundconst.f``.**  Not "none above some threshold of smallness"; none.
A literal left in the Fortran is a literal that the next change of these
values cannot reach.  If you have come here to add an exception, that is
the sentence you have to argue with first.
"""

import os
import re

import pytest

from conftest import REPO

#: A numeric literal: an integer, a decimal, or either with a Fortran
#: exponent.  The look-behind keeps us out of the middle of a name and
#: off the second dot of a range like ``1.D0`` inside ``E12CS(1)``.
LITERAL = re.compile(r'(?<![A-Za-z0-9_.])(\d+\.\d*|\.\d+|\d+)'
                     r'([dDeEqQ][-+]?\d+)?')


def significant(token):
    """How many significant digits a literal carries."""
    mantissa = re.match(r'^[0-9.]+', token).group(0)
    mantissa = mantissa.replace('.', '').lstrip('0')
    return len(mantissa.rstrip('0'))


def sources():
    """Every Fortran source and include in the package."""
    found = []
    for folder in (REPO / 'FOR', REPO / 'FOR' / 'LSJJ'):
        if not folder.is_dir():
            continue
        for path in sorted(folder.iterdir()):
            if path.is_file() and path.suffix.lower() == '.f':
                found.append(path)
    return found


def literals(path):
    """(line number, token) for every real literal of >= 3 digits."""
    out = []
    with open(path, errors='replace') as handle:
        for number, raw in enumerate(handle, 1):
            line = raw.rstrip('\n')
            if not line.strip() or line[0] in 'CcDd*!':
                continue
            code = line[6:72] if len(line) > 6 else ''
            code = re.sub(r"'[^']*'", "''", code).split('!')[0]
            if re.match(r'\s*(FORMAT|format)', code):
                continue
            for found in LITERAL.finditer(code):
                token = found.group(0)
                if '.' not in token and found.group(2) is None:
                    continue                      # an integer: see above
                if significant(token) >= 3:
                    out.append((number, token))
    return out


# --------------------------------------------------------------------------
# the inventory
# --------------------------------------------------------------------------

def _group(reason, tokens):
    return {token: reason for token in tokens.split()}


#: CODATA 2022, transcribed in `FOR/fundconst.f`.  This is the one file
#: that is *supposed* to spell constants out; every other entry below is
#: a literal that survives because it is not a physical constant.
CODATA = _group('CODATA 2022, transcribed in fundconst.f', """
    299792458.D0 2.99792458D10 6.62607015D-34 1.602176634D-19
    1.380649D-23 6.02214076D23 9.1093837139D-31 8.8541878188D-12
    1.25663706127D-6 9.2740100657D-24 5.29177210544D-11 5.29177210544D-9
    4.3597447222060D-18 7.2973525643D-3 109737.31568157D0
    219474.63136314D0 3.14159265358979324D0
    2.00231930436092D0
""")

#: Pure mathematics.  Out of scope for `fundconst.f` on purpose: a
#: measured value can go stale and these cannot.
MATHEMATICS = _group('pure mathematics', """
    .333333333e+00 2.302585093D0 1.224744871D0 0.6875D0
""")

#: The Thomas-Fermi statistical model and Latter's tail correction.
#: Numbers of the model rather than of nature -- 0.88534138 is
#: (1/2)(3*pi/4)**(2/3) -- and they move only if the model does.
THOMAS_FERMI = _group('Thomas-Fermi / Latter tail, a model constant', """
    0.88534138 1.22177412 1.142
""")

#: Quadrature and finite-difference weights: Newton-Cotes 5/288 and
#: 2.5/288, the five-point Adams predictor, the 1/8 of a centred
#: derivative, and the 1/16 mesh step of the HF tail integration.
QUADRATURE = _group('quadrature or finite-difference weight', """
    0.00125D0 288.D0 288.0 1.09375 0.2734375 0.546875 0.21875
    0.0390625 0.125 0.0625
""")

#: Iteration control: damping factors, and the epsilons that make a
#: floating comparison or a truncation to integer come out right.
CONVERGENCE = _group('iteration control or rounding epsilon', """
    1.01D0 0.66667 1.25 0.999999 2.000001D0 1.000001D0 1.002D0
""")

#: Values chosen to fill a print field, or to mean "no value here".
#: They are read back by the parsers in CODE/, so they are as much a
#: part of the file format as the FORMAT statements are.
SENTINELS = _group('print-field limit or sentinel', """
    99999999.D0 9999.99D0 99999.99D0 9999.9999 999.D0 99.999D0
    9.9946D33 9999999.D0 9999.9D0
""")

#: Bookkeeping: a weight that turns a set of small integers into one
#: sortable key, and a percent-scaled overlap written into an integer.
BOOKKEEPING = _group('sort key or scaled integer code', """
    0.008403361D0 0.0084033613D0 100.4D0
""")

#: An empirical correction from the literature, fitted rather than
#: derived.  It has no CODATA value to be updated to.
EMPIRICAL = _group('empirical fit factor', """
    0.07702D0
""")

#: Three DATA tables.  Each is a block of tabulated numbers that is
#: meaningless one entry at a time: the Thomas-Fermi screening function
#: in RCN, and in RCG the average-energy coefficients and the Chebyshev
#: series for the exponential integral E1.
TABLES = _group('DATA table (Thomas-Fermi function; CAVEDN average-energy'
                ' coefficients; E12CS Chebyshev series for E1)', """
    0.05556 0.06311 0.07289 0.07866 0.08517 0.09255 0.10094 0.11053
    0.12152 0.13412 0.14852 0.16493 0.18354 0.19376 0.20468 0.21645
    0.22926 0.24348 0.25953 0.27789 0.29887 0.32262 0.34927 0.36377
    0.37916 0.39553 0.41304 0.43184 0.45214 0.47416 0.49813 0.52427
    0.55285 0.56814 0.58416 0.60097 0.61866 0.63732 0.65710 0.67816
    0.70066 0.72473 0.75044 0.76393 0.77785 0.79220 0.80700 0.82227
    0.83802 0.85429 0.87109 0.88844 0.90634 0.91549 0.92475 0.93413
    0.94359 0.95311 0.96266 0.97218 0.98164 0.99094
    .0000029856D0 .0000074299D0 .0000098278D0 .0000171733D0
    .0000535641D0 .0001073102D0 .0001099913D0 .0003091183D0
    .0008097559D0 .0013461707D0 .0019268451D0 .0072241015D0
    .0134322662D0 .0144191240D0 .0864811785D0 .1125352434D0
    .1303182079D0 .1892918000D0 .6057732466D0 0.0373902147D0
    0.0427239860D0 0.603022689D0 105.D0 108.88888887D0 116.66666667D0
    118.7939392D0 12.681818182D0 12.72792206D0 12.75D0 122.5D0 125.D0
    127.2792206D0 13.74772708D0 134.1640786D0 135.D0 136.D0
    14.96662955D0 141.D0 144.D0 148.492424D0 15.55555555D0 155.D0
    160.9968944D0 17.25D0 17.5D0 175.D0 178.1909089D0 19.09188309D0
    195.D0 2.25D0 2.863636364D0 225.D0 235.D0 248.33333333D0
    27.15015068D0 275.D0 3.75D0 315.D0 322.5D0 33.33333333D0 336.D0
    34.36931771D0 37.41657387D0 39.68626967D0 40.83333333D0 47.6235236D0
    5.031152949D0 52.5D0 54.99090834D0 55.27707984D0 6.136363636D0
    6.568181818D0 6.75D0 66.66666667D0 67.5D0 91.66666667D0
    92.86549413D0 97.22222222D0
""")

INVENTORY = {}
for part in (CODATA, MATHEMATICS, THOMAS_FERMI, QUADRATURE, CONVERGENCE,
             SENTINELS, BOOKKEEPING, EMPIRICAL, TABLES):
    INVENTORY.update(part)


# --------------------------------------------------------------------------
# the tests
# --------------------------------------------------------------------------

@pytest.fixture(scope='module')
def found():
    paths = sources()
    if not paths:
        pytest.skip('FOR/ is not present')
    return {path: literals(path) for path in paths}


def test_every_literal_is_accounted_for(found):
    """No real literal in FOR/ is missing from the inventory above."""
    strangers = []
    for path, hits in found.items():
        for number, token in hits:
            if token not in INVENTORY:
                strangers.append('%s:%d  %s' % (path.name, number, token))
    assert not strangers, (
        'numeric literals that no entry in this file accounts for:\n  '
        + '\n  '.join(strangers)
        + '\n\nIf one is a physical constant it belongs in '
          'FOR/fundconst.f, not in the Fortran.  If it is not, add it to '
          'the right group above with the reason, so the next one is '
          'still visible.')


def test_the_inventory_has_no_stale_entries(found):
    """Nothing is listed here that the Fortran no longer contains.

    An allow-list that keeps entries for literals that have gone stops
    describing the code, and a stale entry silently re-permits a value
    that was deliberately removed -- ``2066.D0`` would be waved through
    if it ever came back.
    """
    live = {token for hits in found.values() for _, token in hits}
    stale = sorted(set(INVENTORY) - live)
    assert not stale, (
        'listed here but no longer in FOR/: %s' % ', '.join(stale))


#: The five sites Section 10 substituted.  They are named rather than
#: described because the point is the exact spelling: if one of these
#: comes back, a constant has been un-substituted.
SECTION_10 = {'2066.D0': 'AURATE, the autoionization rate factor',
              '109735.3D0': 'RINF in the autoionization rate',
              '6.28318D0': 'TWOPIC',
              '8.06662D-18': 'PICSCN, the photoionization cross-section',
              '109.737D0': 'RINFKK in RCN2 G5INP'}


#: The four Section 13 substituted, the last physical constants that
#: were still written inline.  Same rule, same reason.
SECTION_13 = {'0.50116D0': 'GEM1H, (|g_e|-1)/2 in the Lande g-factor',
              '1.50116D0': 'GEP1H, 1+(|g_e|-1)/2 in the Lande g-factor',
              '1.00232D0': 'GEM1, |g_e|-1 in the M1 matrix element',
              '2.1717D-8': 'CMAXWL, the Maxwellian rate coefficient'}

SUBSTITUTED = dict(SECTION_10, **SECTION_13)


def test_the_substituted_literals_have_not_come_back(found):
    back = []
    for path, hits in found.items():
        for number, token in hits:
            if token in SUBSTITUTED:
                back.append('%s:%d  %s (should be %s)'
                            % (path.name, number, token,
                               SUBSTITUTED[token]))
    assert not back, '\n  '.join(back)
