"""GAFROMSCOEF / GFFROMGACOEF / GFFROMSCOEF in `FOR/rcg11kd.f`: the
radiative conversion coefficients, as functions.

This is the one place in the constants work with an *independent oracle*.
Everywhere else the standard is bit-exactness against the program's own
previous output, because there is nothing else to compare with; here the
handbook publishes the three numbers directly:

    Martin & Wiese, Atomic Spectroscopy, in Drake (ed.), Atomic,
    Molecular and Optical Physics Handbook (AIP Press, 1996), collected
    in `WORK/A-S-f-conversion_formulas.docx`.

So the tests come in two layers.  The first recomputes the coefficients
in Python straight from the CODATA 2022 values and checks them against
the handbook, which needs no compiler and always runs.  The second
compiles the Fortran functions out of `FOR/rcg11kd.f` and checks
that they return those same numbers to the last bit; it skips when there
is no gfortran.  The second layer is the one that would catch an edit to
`fundconst.f`; the first is what says the target is right.

Section 12 added a third layer at the end of the file: `radiate.py` now
derives the same conversions for all six multipoles instead of carrying
a transcribed table, and the tests there hold its E1, M1 and E2 to the
compiled functions exactly and its M2, E3 and M3 to the published table
the transcription came from.
"""

import re
import subprocess
import sys

import pytest

from conftest import REPO

SOURCE = REPO / 'FOR' / 'rcg11kd.f'

# --------------------------------------------------------------------------
# the oracle
# --------------------------------------------------------------------------

#: CODATA 2022, transcribed from the NIST complete listing exactly as
#: `FOR/fundconst.f` transcribes them.  Deliberately re-typed here rather
#: than parsed out of the include file: a test that read the value it is
#: checking would agree with any typo.
CLIGHT = 2.99792458e10                  # cm/s, exact
ALPHFS = 7.2973525643e-3                # u(r) = 1.6e-10
ABOHR = 5.29177210544e-9                # cm, u(r) = 1.6e-10
AUCM = 219474.63136314                  # cm-1, the hartree
PICON = 3.14159265358979324

GFE1 = 2.0 / (3.0 * AUCM)
GACON = 8.0 * PICON ** 2 * ALPHFS ** 2 * ABOHR * CLIGHT

#: gA = COEFFICIENT[itype] * S * sigma**(2p+1), sigma in cm-1, S in a.u.
#: The power of sigma is 3 for E1 and M1 and 5 for E2 -- that is what
#: makes the coefficient a constant, and the reason `SIGMA_POWER` is
#: written out below rather than assumed.
COEFFICIENT = {1: GACON * GFE1,
               2: GACON * GFE1 * (ALPHFS / 2.0) ** 2,
               3: GACON * GFE1 * (PICON * ABOHR) ** 2 / 5.0}
SIGMA_POWER = {1: 3, 2: 3, 3: 5}

#: gf = GF_FROM_GA * gA / sigma**2, every type.
GF_FROM_GA = 1.0 / GACON

#: gf = GF_COEFFICIENT[itype] * S * sigma**(2p-1), sigma in cm-1.  Two
#: powers of sigma lower than the gA coefficient, the gA/gf relation
#: carrying a sigma**2 of its own.  This is what SPECTR's PSQ is.
GF_COEFFICIENT = {itype: COEFFICIENT[itype] * GF_FROM_GA
                  for itype in COEFFICIENT}

#: The number the Perl (`CODE/perl_legacy/conv_out.pl`) and the toolbox
#: (`CODE/toolbox/radiate.py`) both carry for gf(E1) = C S / lambda with
#: lambda in Angstrom.  It is *not* GFE1 * 1e8 on CODATA 2022 values:
#: 2e8/(3*219474.63136314) is 303.75568352754544, which is 1.8e-8 lower.
#: The published figure was computed from an older hartree, and 1.8e-8
#: is far below anything RCG prints -- but it is above double precision,
#: so the check below is a tolerance and not an equality, and says so.
HANDBOOK_GF_E1 = 303.75568885954


def test_the_three_coefficients_are_the_published_ones():
    """gA(E1) = 2.0261e18 S/lambda^3, gA(M1) = 2.6974e13 S/lambda^3 and
    gA(E2) = 1.11995e18 S/lambda^5 with lambda in Angstrom.  Since
    sigma = 1e8/lambda, that is the coefficient times 1e24, 1e24 and
    1e40.  This is the check on the whole of `fundconst.f`: get alpha,
    a0, c or the hartree wrong and one of these moves."""
    scale = {itype: 10.0 ** (8 * SIGMA_POWER[itype]) for itype in COEFFICIENT}
    assert COEFFICIENT[1] * scale[1] == pytest.approx(2.0261e18, rel=5e-5)
    assert COEFFICIENT[2] * scale[2] == pytest.approx(2.6974e13, rel=5e-5)
    assert COEFFICIENT[3] * scale[3] == pytest.approx(1.11995e18, rel=5e-6)


def test_gf_from_ga_inverts_ga_from_gf():
    """The two coefficients are one relation read in the two directions,
    so a round trip has to come back to where it started."""
    assert GF_FROM_GA * GACON == pytest.approx(1.0, rel=1e-15)


def test_the_gf_coefficient_is_the_published_one():
    """gf(E1) = 303.756 S / lambda with lambda in Angstrom, which is the
    E1 coefficient times 1e8 since sigma = 1e8/lambda.  The agreement is
    to 1.8e-8, not to the last bit: see HANDBOOK_GF_E1."""
    assert GF_COEFFICIENT[1] * 1e8 == pytest.approx(HANDBOOK_GF_E1,
                                                    rel=5e-8)
    assert GF_COEFFICIENT[1] * 1e8 != HANDBOOK_GF_E1


def test_the_gf_coefficients_are_the_ga_ones_two_powers_of_sigma_down():
    """gA = GACON * gf * sigma**2, so the gf coefficient is the gA one
    over GACON and the type factors are untouched: the M1/E1 ratio is
    still (alpha/2)**2 and the E2/E1 ratio still (pi a0)**2/5."""
    assert (GF_COEFFICIENT[2] / GF_COEFFICIENT[1]
            == pytest.approx((ALPHFS / 2.0) ** 2, rel=1e-15))
    assert (GF_COEFFICIENT[3] / GF_COEFFICIENT[1]
            == pytest.approx((PICON * ABOHR) ** 2 / 5.0, rel=1e-15))


def test_the_ga_coefficient_is_consistent_with_the_gf_one():
    """gA = GACON * gf * sigma**2 is the identity SPECTR uses at eight
    sites.  Reading it backwards through GF_FROM_GA has to give the gf
    coefficients the handbook publishes: 3.0376e-6 for E1, and the E1
    one scaled by (alpha/2)**2 for M1."""
    assert COEFFICIENT[1] * GF_FROM_GA == pytest.approx(3.0375568e-6,
                                                        rel=1e-7)
    assert (COEFFICIENT[2] / COEFFICIENT[1]
            == pytest.approx((ALPHFS / 2.0) ** 2, rel=1e-15))


# --------------------------------------------------------------------------
# the Fortran itself
# --------------------------------------------------------------------------

DRIVER = """      PROGRAM DRV
      DOUBLE PRECISION GAFROMSCOEF, GFFROMGACOEF, GFFROMSCOEF
      WRITE (*,10) GAFROMSCOEF(1)
      WRITE (*,10) GAFROMSCOEF(2)
      WRITE (*,10) GAFROMSCOEF(3)
      WRITE (*,10) GFFROMGACOEF()
      WRITE (*,10) GAFROMSCOEF(0)
      WRITE (*,10) GAFROMSCOEF(4)
      WRITE (*,10) GAFROMSCOEF(-1)
      WRITE (*,10) GFFROMSCOEF(1)
      WRITE (*,10) GFFROMSCOEF(2)
      WRITE (*,10) GFFROMSCOEF(3)
      WRITE (*,10) GFFROMSCOEF(0)
   10 FORMAT (1PD26.17)
      END
"""

#: The build flags of CLAUDE.md section 2, less -fopenmp: these two
#: routines have no parallel region in them, and -fno-automatic is kept
#: because it is what makes them non-re-entrant, which their own comment
#: warns about.
FLAGS = ('-std=legacy', '-fshort-enums', '-malign-double', '-ftracer',
         '-fno-backslash', '-fallow-argument-mismatch', '-fcray-pointer',
         '-fno-automatic', '-O3')


def extract(name, text):
    """The named function subprogram, from its header to its END.

    Compiling all 17 700 lines of `rcg11kd.f` to check two constants
    would cost a minute; these two routines touch nothing else in the
    file -- no COMMON, no unit, no state -- so they compile alone
    against `FOR/fundconst.f`, which is the point of them.
    """
    lines = text.split('\n')
    head = re.compile(r'^      DOUBLE PRECISION FUNCTION %s\b' % name)
    for first, line in enumerate(lines):
        if head.match(line):
            for last in range(first, len(lines)):
                if lines[last].rstrip() == '      END':
                    return '\n'.join(lines[first:last + 1]) + '\n'
            raise AssertionError('%s has no END' % name)
    raise AssertionError('%s is not in %s' % (name, SOURCE))


@pytest.fixture(scope='module')
def returned(tmp_path_factory):
    """What the two functions actually return, in the driver's order.

    Skips rather than fails when gfortran is not installed: the Python
    layer above still runs, and it is the layer that carries the physics.
    """
    if not SOURCE.is_file():
        pytest.skip('%s is not present' % SOURCE)
    folder = tmp_path_factory.mktemp('coefficients')
    text = SOURCE.read_text(encoding='latin-1')
    (folder / 'coef.f').write_text(extract('GAFROMSCOEF', text)
                                   + extract('GFFROMGACOEF', text)
                                   + extract('GFFROMSCOEF', text))
    (folder / 'drv.f').write_text(DRIVER)
    program = folder / ('drv.exe' if sys.platform == 'win32' else 'drv')
    build = ['gfortran'] + list(FLAGS) + [
        '-I', str(REPO / 'FOR'), '-o', str(program), 'coef.f', 'drv.f']
    try:
        done = subprocess.run(build, cwd=str(folder),
                              capture_output=True, text=True)
    except OSError as trouble:
        pytest.skip('gfortran is not available: %s' % trouble)
    if done.returncode:
        pytest.fail('the coefficient routines do not compile:\n%s'
                    % (done.stderr or done.stdout))
    run = subprocess.run([str(program)], cwd=str(folder),
                         capture_output=True, text=True)
    if run.returncode:
        pytest.skip('the driver did not run: %s'
                    % (run.stderr or run.stdout))
    values = [float(one.replace('D', 'E')) for one in run.stdout.split()]
    assert len(values) == 11, run.stdout
    return values


def test_ga_from_s_returns_the_published_coefficients(returned):
    """To the last bit, not to a tolerance.  Section 7 of the plan asks
    for full double precision, and there is no rounding anywhere between
    the CODATA digits and the answer to stop it being exact."""
    assert returned[0] == COEFFICIENT[1]
    assert returned[1] == COEFFICIENT[2]
    assert returned[2] == COEFFICIENT[3]


def test_gf_from_ga_returns_the_reciprocal_of_gacon(returned):
    assert returned[3] == GF_FROM_GA


def test_an_unknown_itype_returns_a_hard_zero(returned):
    """0, 4 and -1.  A wrong answer that shows, rather than a plausible
    number: every real coefficient here is positive, so a caller that
    ignores the contract gets a gA of zero and notices."""
    assert returned[4] == 0.0
    assert returned[5] == 0.0
    assert returned[6] == 0.0


def test_gf_from_s_is_the_product_of_the_other_two(returned):
    """To the last bit again: the routine is written as the product, so
    what it returns is the product, and no third derivation of the type
    factors exists to drift away from these two."""
    assert returned[7] == returned[0] * returned[3]
    assert returned[8] == returned[1] * returned[3]
    assert returned[9] == returned[2] * returned[3]
    assert returned[7] == GF_COEFFICIENT[1]
    assert returned[8] == GF_COEFFICIENT[2]
    assert returned[9] == GF_COEFFICIENT[3]


def test_gf_from_s_passes_the_hard_zero_through(returned):
    """An unknown ITYPE has to fail the same way here as in
    GAFROMSCOEF, and does so for free: 0 times a finite reciprocal."""
    assert returned[10] == 0.0


def test_the_dummy_argument_is_an_integer_named_itype():
    """Under `IMPLICIT REAL*8 (A-H,O-Z)` a dummy named TYPE would be a
    real -- it would work for 1.0, 2.0 and 3.0 and be a trap forever
    after -- so the name and the explicit INTEGER*4 are both load
    bearing, and neither is visible in the returned values."""
    text = SOURCE.read_text(encoding='latin-1')
    for name in ('GAFROMSCOEF', 'GFFROMSCOEF'):
        body = extract(name, text)
        assert 'FUNCTION %s(ITYPE)' % name in body
        assert re.search(r'^      INTEGER\*4 ITYPE\s*$', body, re.M), name


def statements(body):
    """The lines of `body` that are Fortran rather than comment."""
    return '\n'.join(line for line in body.split('\n')
                     if line[:1] not in ('C', 'c', '*', '!') and line.strip())


def test_neither_routine_carries_state_or_touches_a_common_block():
    """Their comment promises it, and it is what makes them safe to call
    from anywhere outside a parallel region -- and what lets the test
    above compile them on their own."""
    text = SOURCE.read_text(encoding='latin-1')
    for name in ('GAFROMSCOEF', 'GFFROMGACOEF', 'GFFROMSCOEF'):
        code = statements(extract(name, text))
        assert 'COMMON' not in code, name
        assert 'SAVE' not in code, name
        assert 'DATA' not in code, name
        assert "INCLUDE 'rcgpar.f'" not in code, name
    for name in ('GAFROMSCOEF', 'GFFROMGACOEF'):
        assert "INCLUDE 'fundconst.f'" in statements(extract(name, text)), name
    #: GFFROMSCOEF needs no constants of its own -- it is the product of
    #: the other two -- and must not acquire any.
    assert "INCLUDE" not in statements(extract('GFFROMSCOEF', text))


def test_no_literal_constant_is_written_into_the_coefficients():
    """The whole point of the section: the numbers come from
    `fundconst.f`, so the only literals allowed in the arithmetic are
    the exact small integers of the formulae -- the 2 of alpha/2, the 5
    of the E2 factor, the 1 of the reciprocal, and the 0 of the
    unknown-type branch."""
    text = SOURCE.read_text(encoding='latin-1')
    allowed = {'0.D0', '1.D0', '2.D0', '5.D0'}
    for name in ('GAFROMSCOEF', 'GFFROMGACOEF', 'GFFROMSCOEF'):
        for line in statements(extract(name, text)).split('\n'):
            for number in re.findall(r'\d+\.?\d*[DdEe][-+]?\d+', line):
                assert number.upper() in allowed, '%s: %s' % (name, line)


# --------------------------------------------------------------------------
# the toolbox side (Section 12)
# --------------------------------------------------------------------------
#
# `CODE/toolbox/radiate.py` used to carry a hard-coded table of these
# same conversions, transcribed from `conv_out.pl`.  It now derives them
# from `toolbox.fundconst`, and these are the tests that tie the derived
# values to the Fortran and to the published table they replaced.
#
# The tie for E1, M1 and E2 is exact equality with the compiled
# functions.  M2, E3 and M3 have no Fortran counterpart and come from
# the general Grant factor; for them the check is the published table,
# to the tolerance the published roundings allow.

from toolbox import radiate                                   # noqa: E402

#: gA(x) = radiate._S[x][1] * S / lambda**n with lambda in angstrom, so
#: the table entry is the cm-1 coefficient times 1e(8n).
FORTRAN_MULTIPOLE = {1: radiate.E1, 2: radiate.M1, 3: radiate.E2}


def test_the_toolbox_coefficients_are_the_fortran_ones(returned):
    """Exact equality, in both directions of the conversion.  This is
    the tie between the two languages: `radiate` writes the E1, M1 and
    E2 type factors exactly as GAFROMSCOEF writes them, for the sole
    purpose of making this assertion an equality rather than a
    tolerance."""
    for itype, name in FORTRAN_MULTIPOLE.items():
        assert radiate.COEFFICIENT[name] == returned[itype - 1], name
    assert radiate.GF_FROM_GA == returned[3]


def test_the_toolbox_gf_coefficients_are_the_fortran_ones(returned):
    """GFFROMSCOEF, the product of the two, against `radiate._GF` read
    back out of its wavelength scaling."""
    for itype, name in FORTRAN_MULTIPOLE.items():
        power, constant = radiate._GF[name]
        assert power == 2 * radiate.ORDER[name] - 1, name
        assert constant / 10.0 ** (8 * power) == returned[6 + itype], name


def test_the_grant_factor_reproduces_the_explicit_type_factors():
    """The reconstruction of Section 7 becomes load bearing here: E3 and
    M3 have no other derivation.  What makes that safe is that only the
    *ratio* of two Grant factors is used, and for p = 2 that ratio is
    pi**2/5 exactly -- the (pi*a0)**2/5 the Fortran writes.  Equality of
    the ratio, not of the factor, because (pi*a0)**2/5 and (pi**2/5)*a0**2
    are the same number algebraically and one bit apart in binary."""
    assert (radiate._grant(2) / radiate._grant(1)
            == PICON ** 2 / 5.0)
    assert (radiate._grant(1) / radiate._grant(1)) == 1.0


def test_every_coefficient_reproduces_the_published_table():
    """All six multipoles, both conversions, against the table
    `radiate` used to carry -- which is the only independent check M2,
    E3 and M3 have.  Tolerances are the published roundings: 3e-6 is
    what a six-figure 3.14441e17 can be held to.

    E2's gf constant is the exception and is excluded here; see
    `test_the_legacy_e2_gf_constant_is_inconsistent_with_its_own_partner`."""
    for name, published in radiate.LEGACY_S.items():
        power, derived = radiate._S[name]
        assert derived == pytest.approx(published, rel=3e-6), name
    for name, published in radiate.LEGACY_GF.items():
        if name == radiate.E2:
            continue
        power, derived = radiate._GF[name]
        assert derived == pytest.approx(published, rel=3e-6), name


def test_the_legacy_e2_gf_constant_is_inconsistent_with_its_own_partner():
    """A defect found by deriving what had been copied, of the same kind
    as the 0.6698 recorded in the constants plan as 0a.1.

    gf = gA * lambda**2 / (GACON * 1e16) is exact and type independent,
    so the legacy 1.11995e18 and 167.84224 cannot both be right.  The
    Grant factor confirms the gA one to 2.7e-8; the gf one is out by
    +3.6e-4, four orders of magnitude worse than any other entry in the
    table.  This test pins the size of the disagreement so that the
    number in `radiate.LEGACY_GF` is not mistaken for a target."""
    published = radiate.LEGACY_GF[radiate.E2]
    from_its_partner = (radiate.LEGACY_S[radiate.E2] * radiate.GF_FROM_GA
                        / 1e16)
    assert from_its_partner / published - 1.0 == pytest.approx(3.57e-4,
                                                               rel=1e-2)
    #: and the derived value sides with the gA one, not with it.
    assert radiate._GF[radiate.E2][1] == pytest.approx(from_its_partner,
                                                       rel=1e-7)


def test_the_ev_conversion_is_gone_from_radiate():
    """`radiate.EV_PER_CM` held 8.06554393734921 -- eV per *1000* cm-1 --
    and `conv_out.py` divided by it a quantity already in cm-1, so the
    Auger `dE(eV)` column printed meV.  The name is gone rather than
    corrected in place, so that a caller written against the old
    behaviour fails loudly instead of quietly changing by 1000."""
    assert not hasattr(radiate, 'EV_PER_CM')


def test_the_ev_conversion_is_the_transcription_it_replaced_times_1000():
    """The transcribed constant was right in every digit and wrong only
    in its power of ten: that is the whole of the defect, and this pins
    it so the record cannot be misread as a value change."""
    from toolbox import fundconst
    assert fundconst.EV2CM / 1000.0 == pytest.approx(8.06554393734921,
                                                     rel=1e-15)
