"""Make CODE/ importable and locate the reference runs in WORK/."""

import os
import sys
from pathlib import Path

import pytest

CODE = Path(__file__).resolve().parent.parent
REPO = CODE.parent
WORK = REPO / 'WORK'

#: Where the Perl utilities live since phase 11 retired them.  They are
#: not the shipped programs any more, but they are still the reference
#: every port from phase 6 onward was accepted against, so the
#: comparison tests keep running them from here.  `conv_out.pl` also
#: needs this on its `-I` path, because `conv_cowan.pl` moved with it.
PERL = CODE / 'perl_legacy'

if str(CODE) not in sys.path:
    sys.path.insert(0, str(CODE))


@pytest.fixture(scope='session')
def work():
    """The WORK directory of the repository."""
    if not WORK.is_dir():
        pytest.skip('WORK/ is not present')
    return WORK


def case(name):
    """Skip a test when its reference run is not there to be read.

    The reference cases in WORK/ ship as **inputs only** -- OUTG11 alone
    is 122 MB for Pr III -- so the distribution has the decks and their
    README.txt and nothing generated.  Every test below reads a finished
    run, so it is the presence of OUTG11, not of the folder, that says
    whether there is anything to compare against.  Run the case as its
    README.txt describes and the tests come back.
    """
    folder = WORK / name
    if not folder.is_dir():
        pytest.skip('WORK/%s is not present' % name)
    if not (folder / 'OUTG11').is_file():
        pytest.skip('WORK/%s has not been run: no OUTG11 in it '
                    '(see its README.txt)' % name)
    return folder


def perl_env():
    """The environment a Perl reference run needs.

    The decks and every output file in WORK/ have CRLF line endings, and
    the Perl utilities count columns from the start of the line and
    anchor with `$`.  Strawberry Perl -- the one the .bat files use --
    strips the CR on input, so those regexes match; a Perl built for
    Cygwin or MSYS does not, and every fixed-column regex in
    `conv_cowan.pl` then fails.  PERLIO=crlf asks for the layer by name,
    which makes the comparison run the same way whichever Perl is first
    on the PATH.
    """
    env = dict(os.environ)
    env['PERLIO'] = 'crlf'
    return env
