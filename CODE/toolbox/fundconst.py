"""The Python side of `FOR/fundconst.f`: one source for the constants.

`FOR/fundconst.f` holds every fundamental physical constant the Fortran
uses.  The Python toolbox needs some of the same numbers -- the
radiative conversion coefficients of `radiate.py` above all -- and until
now carried its own transcribed copies of them, which a change to the
Fortran could not reach.  This module is the fix: the same CODATA 2022
values, entered with the same derivations, in the same order, under the
same names.

**This file is not the source of truth; `FOR/fundconst.f` is.**  What
ties them together is `CODE/tests/test_fundconst.py`, which compiles the
include file into a driver, prints every parameter it declares, and
asserts that each one equals its counterpart here *bit for bit*.  A name
added to the Fortran and not to this module fails that test, and so does
a digit changed in either.  The test skips where there is no gfortran,
so the toolbox stays importable on a machine that has no compiler --
`CODE/` is what goes on the path, `FOR/` is not distributed with it, and
a module that fails to import without the source tree would be a bad
trade for a value that changes once a decade.

The values are deliberately re-typed here from the NIST complete listing
(https://physics.nist.gov/cuu/Constants/Table/allascii.txt, the 2022
adjustment) rather than parsed out of the Fortran.  A test that read the
number it is checking would agree with any typo.

ORDER MATTERS, in the same way it does in the Fortran: a derived
quantity comes after what it is derived from, and the expressions are
written in the same association as the Fortran writes them, because
float arithmetic is not associative and the tie between the two
languages is exact equality rather than a tolerance.
"""

import math

# =========================================================================
# SECTION A -- CODATA 2022, exact as published.
# =========================================================================

#: Speed of light in vacuum, m/s and cm/s.  Exact by definition.
CLIGHM = 299792458.0
CLIGHT = 2.99792458e10

#: Planck constant, J s.  Exact by definition of the kilogram.
PLANCK = 6.62607015e-34

#: Elementary charge, C.  Exact by definition of the ampere.
ECHARG = 1.602176634e-19

#: Boltzmann constant, J K-1.  Exact by definition of the kelvin.
KBOLTZ = 1.380649e-23

#: Avogadro constant, mol-1.  Exact by definition of the mole.
AVOGAD = 6.02214076e23

#: Electron mass, kg.                       u(r) = 3.1e-10
EMASS = 9.1093837139e-31

#: Vacuum electric permittivity, F m-1.     u(r) = 1.6e-10
EPS0 = 8.8541878188e-12

#: Vacuum magnetic permeability, N A-2.     u(r) = 1.6e-10
AMU0 = 1.25663706127e-6

#: Bohr magneton, J T-1.                    u(r) = 3.1e-10
BOHRMG = 9.2740100657e-24

#: Electron g-factor, magnitude.            u(r) = 1.3e-13
#: CODATA 2022 publishes it as -2.00231930436092(36); the package uses
#: |g| everywhere, so the magnitude is what is stored.
GEFACT = 2.00231930436092

#: Bohr radius, m and cm.                   u(r) = 1.6e-10
ABOHRM = 5.29177210544e-11
ABOHR = 5.29177210544e-9

#: Hartree energy, J.                       u(r) = 1.1e-12
HARTJ = 4.3597447222060e-18

#: Fine-structure constant, and its square. u(r) = 1.6e-10
#: 1/ALPHFS = 137.035999177(21).
ALPHFS = 7.2973525643e-3
ALPHF2 = ALPHFS * ALPHFS

#: Rydberg constant, cm-1.                  u(r) = 1.1e-12
RINFCD = 109737.31568157

#: Hartree in cm-1, i.e. 2*Rinf: the atomic unit of energy in the units
#: RCN, RCN2 and RCG print.
AUCM = 219474.63136314


# =========================================================================
# SECTION B -- the constants the code actually uses.
#
# The comments here are shorter than the Fortran's on purpose: each entry
# there carries the legacy literal it replaced and the relative error
# that literal had, and duplicating that history in two places is exactly
# the kind of drift this file exists to stop.  Read `FOR/fundconst.f` for
# it.
# =========================================================================

#: Rinf in cm-1, and in 1000 cm-1 (the unit RCN2 carries parameters in).
RINF = RINFCD
RINFKK = RINFCD / 1.0e3

#: 1/alpha and 2/alpha.
ALPHIN = 1.0 / ALPHFS
TWOOAL = 2.0 / ALPHFS

#: alpha**2/2, the spin-orbit coefficient, and alpha**2/4.
ALPH2H = ALPHF2 / 2.0
ALPH2Q = ALPHF2 / 4.0

#: (alpha**2/2)*Rinf and alpha**2*Rinf, in cm-1.
ZETCON = ALPHF2 * RINFCD / 2.0
ZETCN2 = ALPHF2 * RINFCD

#: eV -> cm-1, and keV -> cm-1.  Exact: it follows from the exact e, h
#: and c, so it is written as that quotient rather than transcribed.
EV2CM = ECHARG / (PLANCK * CLIGHT)
EV2CMK = EV2CM * 1.0e3

#: Rydberg energy in eV and in keV.
RYEV = RINFCD / EV2CM
RYKEV = RYEV / 1.0e3

#: pi, and 2*pi.
PICON = 3.14159265358979324
TWOPIC = 2.0 * PICON

#: 4/3.
FOURO3 = 4.0 / 3.0

#: gf(E1) = GFE1 * sigma(cm-1) * S(a.u.), with GFE1 = 2/(3*2*Rinf).
GFE1 = 2.0 / (3.0 * AUCM)

#: gA = GACON * gf * sigma(cm-1)**2, with
#: GACON = 8*pi**2*alpha**2*a0(cm)*c(cm/s).
GACON = 8.0 * PICON ** 2 * ALPHF2 * ABOHR * CLIGHT

#: hbar, J s.  Exact from the exact h.
HBARJ = PLANCK / (2.0 * PICON)

#: Ry/hbar in units of 1e13 s-1: the factor that turns 2*pi*|V|**2, with
#: V in rydbergs, into an autoionization rate in s-1.
AURATE = HARTJ / (2.0 * HBARJ * 1.0e13)

#: 4*pi**2*alpha*a0**2, cm**2: the photoionization cross-section
#: coefficient CSIGMA.
PICSCN = 4.0 * PICON ** 2 * ALPHFS * ABOHR ** 2

#: The electron g-factor in the three arrangements the Fortran needs:
#: the two coefficients of the Lande g-factor, and the spin part of the
#: M1 matrix element.
GEM1H = (GEFACT - 1.0) / 2.0
GEP1H = 1.0 + GEM1H
GEM1 = GEFACT - 1.0

#: 2*sqrt(pi)*alpha*c*a0**2, cm3 s-1: the Maxwellian rate coefficient
#: per unit collision strength, with the temperature in rydbergs.
CMAXWL = 2.0 * math.sqrt(PICON) * ALPHFS * CLIGHT * ABOHR ** 2

