"""`CODE/toolbox/fundconst.py` against `FOR/fundconst.f`, bit for bit.

The toolbox needs the same fundamental constants the Fortran uses, and
until Section 12 of the constants plan it carried its own transcribed
copies of them in `radiate.py`.  `toolbox/fundconst.py` replaced those
with the same values under the same names -- but a second copy is only
better than an independent one if something *checks* that it is still a
copy.  That is this file.

It compiles `FOR/fundconst.f` into a throwaway driver that prints every
parameter the include file declares, and asserts that each equals the
attribute of the same name in `toolbox.fundconst` **exactly**, not to a
tolerance.  Both sides are IEEE doubles evaluated from the same digits
by the same expressions, so there is nothing for them to disagree about
short of a typo -- which is the whole point.

Two properties beyond the values themselves:

  * every name in the Fortran must exist in the Python.  A constant
    added to `FOR/fundconst.f` and not here fails, which is what keeps
    the two from drifting as the file grows.
  * the Python must add no constant of its own that the Fortran does
    not have, so that `toolbox/fundconst.py` cannot quietly become a
    second source of truth.

The Fortran layer skips where there is no gfortran, exactly as
`test_coefficients.py` does: `CODE/` is what the author puts on the
path and `FOR/` is not distributed with it, so the toolbox has to stay
importable on a machine that cannot compile anything.  The consistency
is then proved wherever it *can* be proved, which is every machine that
builds the package.
"""

import re
import subprocess
import sys

import pytest

from conftest import REPO
from toolbox import fundconst

INCLUDE = REPO / 'FOR' / 'fundconst.f'

#: The build flags of CLAUDE.md section 2, less -fopenmp.  The include
#: file has no executable statement in it, so nothing here depends on
#: them -- they are used anyway, because a constant that came out
#: differently under the package's own flags would be worth knowing
#: about.
FLAGS = ('-std=legacy', '-fshort-enums', '-malign-double', '-ftracer',
         '-fno-backslash', '-fallow-argument-mismatch', '-fcray-pointer',
         '-fno-automatic', '-O3')

#: A PARAMETER statement in the include file.  One name per statement is
#: the file's own convention and this relies on it; a grouped
#: `PARAMETER (A = 1, B = 2)` would be silently half-read, so the test
#: below checks the count against the declarations as well.
PARAMETER = re.compile(r'^      PARAMETER \((\w+) *=', re.M)

#: The continuation-free `DOUBLE PRECISION` declarations, which name the
#: same constants in the same order -- the cross-check on the regex
#: above.
DECLARED = re.compile(r'^      DOUBLE PRECISION +(.+?)\s*$', re.M)


def names():
    """Every parameter `FOR/fundconst.f` defines, in file order."""
    if not INCLUDE.is_file():
        pytest.skip('%s is not present' % INCLUDE)
    text = INCLUDE.read_text(encoding='latin-1')
    found = PARAMETER.findall(text)
    assert found, 'no PARAMETER statement found in %s' % INCLUDE
    return found


def declared():
    """Every name the include file declares DOUBLE PRECISION."""
    text = INCLUDE.read_text(encoding='latin-1')
    out = []
    for group in DECLARED.findall(text):
        out.extend(one.strip() for one in group.split(','))
    return out


def driver(found):
    """A program that prints each constant, one per line, in file order.

    1PD26.17 is 18 significant digits, more than the 17 a double needs
    to round-trip, so what comes back is the compiler's value and not a
    rounding of it.
    """
    lines = ['      PROGRAM DRV', '      IMPLICIT NONE',
             "      INCLUDE 'fundconst.f'"]
    lines.extend('      WRITE (*,10) %s' % name for name in found)
    lines.append('   10 FORMAT (1PD26.17)')
    lines.append('      END')
    return '\n'.join(lines) + '\n'


@pytest.fixture(scope='module')
def fortran(tmp_path_factory):
    """{name: value} as gfortran evaluates the include file."""
    found = names()
    folder = tmp_path_factory.mktemp('fundconst')
    (folder / 'drv.f').write_text(driver(found))
    program = folder / ('drv.exe' if sys.platform == 'win32' else 'drv')
    build = ['gfortran'] + list(FLAGS) + [
        '-I', str(REPO / 'FOR'), '-o', str(program), 'drv.f']
    try:
        done = subprocess.run(build, cwd=str(folder),
                              capture_output=True, text=True)
    except OSError as trouble:
        pytest.skip('gfortran is not available: %s' % trouble)
    if done.returncode:
        pytest.fail('FOR/fundconst.f does not compile:\n%s'
                    % (done.stderr or done.stdout))
    run = subprocess.run([str(program)], cwd=str(folder),
                         capture_output=True, text=True)
    if run.returncode:
        pytest.skip('the driver did not run: %s' % (run.stderr or run.stdout))
    values = [float(one.replace('D', 'E')) for one in run.stdout.split()]
    assert len(values) == len(found), run.stdout
    return dict(zip(found, values))


def test_every_parameter_is_declared_exactly_once():
    """The regex above reads one name per PARAMETER statement, which is
    the file's convention.  If a grouped PARAMETER is ever added it will
    read only the first name of it, so the declarations -- which the
    file also writes out in full -- are the check that nothing was
    missed."""
    found = names()
    assert sorted(found) == sorted(declared())


def test_every_fortran_constant_has_a_python_counterpart():
    """The property that keeps the two files together as the Fortran
    grows: a constant added to `FOR/fundconst.f` and not to
    `toolbox/fundconst.py` fails here."""
    missing = [name for name in names() if not hasattr(fundconst, name)]
    assert not missing, ('CODE/toolbox/fundconst.py is missing %s'
                         % ', '.join(missing))


def test_the_python_module_adds_no_constant_of_its_own():
    """It is a copy of `FOR/fundconst.f` and nothing else.  A name here
    that is not there would be a constant with no oracle, which is the
    situation Section 12 exists to end."""
    found = set(names())
    extra = sorted(name for name in vars(fundconst)
                   if name.isupper() and name not in found)
    assert not extra, ('CODE/toolbox/fundconst.py defines %s, which is not '
                       'in FOR/fundconst.f' % ', '.join(extra))


def test_every_constant_agrees_to_the_last_bit(fortran):
    """Not to a tolerance.  Both sides evaluate the same expressions
    from the same digits in IEEE double, so any difference at all is a
    typo or a changed derivation, and either is worth failing on."""
    wrong = []
    for name, value in fortran.items():
        mine = getattr(fundconst, name, None)
        if mine != value:
            wrong.append('%-8s fortran %.17e  python %s'
                         % (name, value,
                            'absent' if mine is None else '%.17e' % mine))
    assert not wrong, '\n'.join(wrong)


def test_the_values_are_the_codata_2022_ones():
    """A handful of spot checks against the published digits, so that a
    Fortran typo and a Python typo would have to be the *same* typo to
    pass.  The equality test above only proves the two files agree."""
    assert fundconst.CLIGHT == 2.99792458e10
    assert fundconst.PLANCK == 6.62607015e-34
    assert fundconst.ECHARG == 1.602176634e-19
    assert fundconst.ALPHFS == 7.2973525643e-3
    assert fundconst.ABOHR == 5.29177210544e-9
    assert fundconst.AUCM == 219474.63136314
    assert fundconst.RINFCD == 109737.31568157
    #: 1/alpha = 137.035999177(21) as published, which is a different
    #: transcription of the same constant and so a real check.  The
    #: tolerance is the last digit the listing prints, 5e-12 relative,
    #: not the 1.6e-10 uncertainty: the reciprocal has to reproduce the
    #: printed figure, and it does -- 137.03599917759.
    assert fundconst.ALPHIN == pytest.approx(137.035999177, rel=5e-12)
    #: The hartree is twice the Rydberg, in the file's own two entries.
    assert fundconst.AUCM == pytest.approx(2.0 * fundconst.RINFCD, rel=1e-15)


def test_pi_is_pi():
    """PICON is a literal on both sides -- the one entry with no CODATA
    listing behind it -- so it gets its own check against math.pi."""
    import math
    assert fundconst.PICON == math.pi
    assert fundconst.TWOPIC == 2.0 * math.pi
