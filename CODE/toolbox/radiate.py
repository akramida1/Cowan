"""Line strength, oscillator strength and transition probability.

The three ways of quoting the same number, and the conversions between
them, for the six multipoles Cowan's codes can produce.

The constants are **derived**, not transcribed.  They used to be a table
copied out of `conv_out.pl:1650-1710`, which is where the author had
collected the NIST Atomic Spectroscopy conversions
(https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas);
that made them a second, independent statement of the same physics as
`FOR/rcg11kd.f`'s GAFROMSCOEF, unreachable by any change to
`FOR/fundconst.f`.  Section 12 of the constants plan replaced the table
with the formulae, built on `toolbox.fundconst`, which is itself held to
the Fortran bit for bit by `CODE/tests/test_fundconst.py`.

    gA = COEFFICIENT[multipole] * S * sigma**(2p+1)
    gf = COEFFICIENT[multipole] * S * sigma**(2p-1) / GACON

with sigma in cm-1, S in atomic units, gA = g(upper)*A in s-1 and
gf = g(lower)*f.  The tables below are the same coefficients rewritten
for lambda in angstroms, which is the unit this module's callers work
in: sigma = 1e8/lambda, so a coefficient carrying sigma**n becomes that
coefficient times 1e(8n).

E1, M1 and E2 are written exactly as `GAFROMSCOEF` writes them and
agree with it to the last bit -- `CODE/tests/test_coefficients.py` is
the tie.  E3 and M3 have no Fortran counterpart and come from the
general Grant factor `_grant` below; see its docstring for the check
that it reproduces the explicit forms.

Wavelengths are in angstroms, S in atomic units, gA in s-1.  A caller
working from Cowan's output has energies in the section's own unit --
1000 cm-1 in every deck here -- so a wavelength is `1e5 / delta E`.

`conv_out.py` and `calc_TP_M1E2.py` both use this module; it is the one
copy of these conversions on the Python side.

`CODE/perl_legacy/` still carries the old hard-coded table.  It is the
reference for *intent* (CLAUDE.md section 8) and is not maintained: its
constants are knowingly stale from here on, and one of them is worse
than stale -- see LEGACY_GF below.  Do not "restore" them.
"""

from . import fundconst as fc

E1, M1, E2, M2, E3, M3 = 'E1', 'M1', 'E2', 'M2', 'E3', 'M3'

#: The multipole order p of each name.  The power of sigma in the gA
#: relation is 2p+1 and in the gf relation 2p-1.
ORDER = {E1: 1, M1: 1, E2: 2, M2: 2, E3: 3, M3: 3}


def _double_factorial(n):
    out = 1
    while n > 1:
        out *= n
        n -= 2
    return out


def _grant(p):
    """The Grant factor Np of the general electric-multipole rate.

        Np = 2**(2p+2) * (p+1) * (2p+1) * pi**(2p+1) / (p * [(2p+1)!!]**2)

    and gA(Ep) = Np * S(cgs) * sigma**(2p+1) / hbar.  Reconstructed in
    Section 7 of the constants plan from `WORK/A-S-f-conversion_formulas.docx`,
    whose equation objects do not survive text extraction, and verified
    there against the document's explicit E1, E2 and E3 forms.

    Only the *ratio* of two Grant factors is used here, which is what
    makes the reconstruction safe: whatever prefactor the derivation
    carries divides out, and the E1 anchor comes from `GAFROMSCOEF`
    instead.  For p = 2 the ratio is exactly pi**2/5, which is the
    (pi*a0)**2/5 the Fortran writes; a test asserts it.
    """
    return (2.0 ** (2 * p + 2) * (p + 1) * (2 * p + 1)
            * fc.PICON ** (2 * p + 1) / (p * _double_factorial(2 * p + 1) ** 2))


#: gA = GA_E1 * S * sigma**3 for E1, the anchor for everything else.
#: This is GAFROMSCOEF(1) and equals it exactly.
GA_E1 = fc.GACON * fc.GFE1

#: The type factor: COEFFICIENT[x] = GA_E1 * FACTOR[x].
#:
#: E1, M1 and E2 are spelled as `GAFROMSCOEF` spells them rather than
#: through `_grant`, and that is deliberate.  The two are the same
#: number algebraically but not necessarily the same double -- (pi*a0)**2/5
#: and (pi**2/5)*a0**2 need not round alike -- and for these three the
#: requirement is exact equality with the Fortran, not agreement with
#: it.  E3 has no Fortran counterpart, so it takes the general form.
#: A magnetic multipole is its electric partner times (alpha/2)**2, for
#: every p; that is the same relation the Fortran uses for M1.
FACTOR = {E1: 1.0,
          E2: (fc.PICON * fc.ABOHR) ** 2 / 5.0,
          E3: _grant(3) * fc.ABOHR ** 4 / _grant(1)}
FACTOR[M1] = FACTOR[E1] * (fc.ALPHFS / 2.0) ** 2
FACTOR[M2] = FACTOR[E2] * (fc.ALPHFS / 2.0) ** 2
FACTOR[M3] = FACTOR[E3] * (fc.ALPHFS / 2.0) ** 2

#: gA = COEFFICIENT[x] * S * sigma**(2p+1), sigma in cm-1, S in a.u.
COEFFICIENT = dict((name, GA_E1 * FACTOR[name]) for name in ORDER)

#: gf = GF_FROM_GA * gA / sigma**2, every multipole and both parities.
#: The relation between A and f is not type dependent -- see
#: GFFROMGACOEF in `FOR/rcg11kd.f`, which is this same number.
GF_FROM_GA = 1.0 / fc.GACON

#: (power of the wavelength, constant) in  gf = S * constant / lambda**n
_GF = dict((name, (2 * ORDER[name] - 1,
                   COEFFICIENT[name] * GF_FROM_GA
                   * 10.0 ** (8 * (2 * ORDER[name] - 1))))
           for name in ORDER)

#: (power of the wavelength, divisor) in  S = gA * lambda**n / divisor
_S = dict((name, (2 * ORDER[name] + 1,
                  COEFFICIENT[name] * 10.0 ** (8 * (2 * ORDER[name] + 1))))
          for name in ORDER)

#: What the two tables above used to hold, transcribed from
#: `conv_out.pl` and kept here for the record.  `test_conv_out.py`
#: checks the new values against these, which is the only check on M2,
#: E3 and M3 there is -- and the reason the check has a tolerance is
#: that these are the published roundings, not the derivation.
#:
#: Five of the six gf constants agree with the derivation to 3e-7 or
#: better.  **The sixth, E2, is out by +3.6e-4**, and it is out against
#: its own gA partner: gf = gA*lambda**2/(GACON*1e16) is exact and type
#: independent, so 1.11995e18 and 167.84224 cannot both be right, and
#: the Grant factor confirms the gA one to 2.7e-8.  167.902 is the value
#: consistent with the rest of the table.  This is a defect in the
#: legacy table of the same kind as the 0.6698 recorded in the plan as
#: defect 0a.1, found the same way -- by deriving what had been copied.
LEGACY_GF = {E1: 303.75568885954, M1: 0.0040438504,
             E2: 167.84224,       M2: 0.002235255,
             E3: 47.140897,       M3: 0.000627579}
LEGACY_S = {E1: 2.0261269e18, M1: 269735e8,
            E2: 1.11995e18,   M2: 14909714e6,
            E3: 3.14441e17,   M3: 4.1861e12}

#: `EV_PER_CM` is gone.  It held 8.06554393734921, which is eV per
#: *1000* cm-1, and both of `conv_out.py`'s call sites divided by it a
#: quantity already converted to cm-1 -- so the Auger `dE(eV)` column
#: came out a thousand times too large, in meV under an eV heading.  The
#: Perl does the same (`conv_out.pl:733,965,968`) and its comment names
#: the constant "eV to cm-1", which is 8065.54: the value was the
#: mistake, not the usage.  Both sites now take `fundconst.EV2CM`
#: directly and the column reads in eV.  This is a deliberate,
#: documented divergence from `CODE/perl_legacy/`, decided by the author
#: on 2026-08-21; see Section 12 of PLAN_fundconst.md.

#: The multipole an OUTG11 spectrum section holds (`spectrum.MULTIPOLES`).
OF_SPECTRUM = {'ELEC DIP': E1, 'MAG DIP': M1, 'ELEC QUD': E2}


class RadiateError(Exception):
    """A multipole this module does not have the constants for."""


def _pair(table, multipole):
    try:
        return table[multipole or E1]
    except KeyError:
        raise RadiateError('%r is not one of %s'
                           % (multipole, ', '.join(sorted(table))))


def gf_from_s(s, wavelength, multipole=E1):
    power, constant = _pair(_GF, multipole)
    return s * constant / wavelength ** power


def s_from_gf(gf, wavelength, multipole=E1):
    power, constant = _pair(_GF, multipole)
    return gf * wavelength ** power / constant


def s_from_ga(ga, wavelength, multipole=E1):
    power, divisor = _pair(_S, multipole)
    return ga * wavelength ** power / divisor


def ga_from_s(s, wavelength, multipole=E1):
    power, divisor = _pair(_S, multipole)
    return s * divisor / wavelength ** power


def reduced_absorption(ga, wavelength, multipole=E1):
    """gf * lambda**2, the quantity an absorption spectrum is drawn in.

    See the NIST tables' note on relative intensities in absorption:
    what an absorption line's depth follows is the oscillator strength
    weighted by the wavelength, not the transition probability.
    """
    s = s_from_ga(ga, wavelength, multipole)
    return gf_from_s(s, wavelength, multipole) * wavelength * wavelength
