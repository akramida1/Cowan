"""Wavelength conversion between vacuum and standard air.

The index of refraction is E. R. Peck and K. Reeder, JOSA 62 (1972)
958-962.  A port of CODE/vacair.pl, including its convention that
wavelengths at or below 1500 A are left in vacuum.
"""

_TOL = 1e-12
_AIR_LIMIT = 1500.0


def n_air_minus_1(sigma):
    """(n - 1) of standard air at wavenumber `sigma` in cm-1."""
    s2 = sigma * sigma
    return (8060.51e-8
            + 2480990.0 / (132.274e8 - s2)
            + 17455.7 / (39.32957e8 - s2))


def to_air(lambda_vac):
    """Wavelength in air, in angstroms, from the vacuum wavelength."""
    return lambda_vac / (1.0 + n_air_minus_1(1e8 / lambda_vac))


def to_vacuum(lambda_air):
    """Wavelength in vacuum, in angstroms, from the air wavelength.

    Iterated because the index depends on the vacuum wavelength.  Below
    1500 A the wavelength is taken to be in vacuum already.
    """
    if lambda_air <= _AIR_LIMIT:
        return lambda_air
    previous = lambda_air
    while True:
        current = lambda_air * (1.0 + n_air_minus_1(1e8 / previous))
        if abs(current - previous) <= _TOL:
            return current
        previous = current
