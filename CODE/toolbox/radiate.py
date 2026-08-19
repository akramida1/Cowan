"""Line strength, oscillator strength and transition probability.

The three ways of quoting the same number, and the conversions between
them, for the six multipoles Cowan's codes can produce.  The constants
are transcribed from `conv_out.pl:1650-1710`, which is where the author
collected them; they are the NIST Atomic Spectroscopy conversions
(https://www.nist.gov/pml/atomic-spectroscopy-compendium-basic-ideas).

Wavelengths are in angstroms, S in atomic units, gA in s-1.  A caller
working from Cowan's output has energies in the section's own unit --
1000 cm-1 in every deck here -- so a wavelength is `1e5 / delta E`.

`calc_TP_M1E2.py` (phase 10) needs the same constants; they live here so
that there is one copy of them.
"""

E1, M1, E2, M2, E3, M3 = 'E1', 'M1', 'E2', 'M2', 'E3', 'M3'

#: (power of the wavelength, constant) in  gf = S * constant / lambda**n
_GF = {E1: (1, 303.75568885954),
       M1: (1, 0.0040438504),
       E2: (3, 167.84224),
       M2: (3, 0.002235255),
       E3: (5, 47.140897),
       M3: (5, 0.000627579)}

#: (power of the wavelength, divisor) in  S = gA * lambda**n / divisor
_S = {E1: (3, 2.0261269e18),
      M1: (3, 269735e8),
      E2: (5, 1.11995e18),
      M2: (5, 14909714e6),
      E3: (7, 3.14441e17),
      M3: (7, 4.1861e12)}

#: CODATA 2018 and later.
EV_PER_CM = 8.06554393734921

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
