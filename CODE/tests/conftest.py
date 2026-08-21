"""Make CODE/ importable and locate the reference runs in WORK/."""

import os
import shutil
import subprocess
import sys
import tempfile
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




# --------------------------------------------------------------------------
# the permanent files a run of the Fortran needs
# --------------------------------------------------------------------------

#: The cfp binary decks RCG reads.  They are built on the machine the
#: package is installed on -- make_cfp.bat, or unix/make_cfp -- and are
#: not part of the distribution, so a source tree has none of them.
CFP_DECKS = ('FOR072', 'FOR073', 'FOR074')

#: The LS-to-jj deck, built by the same run of make_cfp since 2026 and
#: not distributed either.  It is kept apart from CFP_DECKS because it
#: is *optional*: its generator, CODE/lsjj/mkdeck.py, needs numpy, and
#: where numpy is absent make_cfp warns, builds the three above and
#: exits 0.  RCG then runs everything except true jj coupling.  So a
#: missing LSJJ.DAT is not a failed deck build and must not skip the
#: session fixture; it only makes the true-jj tests unrunnable.
LSJJ_DECK = 'LSJJ.DAT'

#: The four programs, spelled the way each platform's build leaves them:
#: gfortran is told `-o NAME.exe` on Windows and `-o NAME` on Unix.
PROGRAMS = ('RCN36K', 'RCN2K', 'rcg11kd', 'RCE20KD')
SUFFIX = '.exe' if os.name == 'nt' else ''
EXECUTABLES = tuple(name + SUFFIX for name in PROGRAMS)
RCG = 'rcg11kd' + SUFFIX


def scratch_root():
    """Where the permanent files a test run needs are kept.

    Outside the repository, and never CODE/: RCG opens FOR072, FOR073
    and FOR074 in the directory COWAN.CFG names, with STATUS='UNKNOWN',
    so pointing a run at a CODE/ that has no decks creates three empty
    ones there and leaves them behind when RCG then dies on them.
    Nothing may write into CODE/.

    The directory is kept between runs -- building the decks is a whole
    run of RCG -- and rebuilt whenever the executable or the cfp deck it
    was built from is newer than it, so a recompile is picked up on its
    own.  COWAN_TEST_CODE names it somewhere else.
    """
    named = os.environ.get('COWAN_TEST_CODE')
    if named:
        return Path(named)
    return Path(tempfile.gettempdir()) / 'cowan-tests' / 'CODE'


def _stale_deck(deck, executable, cfp):
    """True when the deck is missing, empty, or older than what built it.

    Empty counts as missing: three empty decks are exactly what an
    interrupted run leaves behind, RCG having opened them before finding
    it had nothing to read.
    """
    if not deck.is_file() or deck.stat().st_size == 0:
        return True
    when = deck.stat().st_mtime
    return when < executable.stat().st_mtime or when < cfp.stat().st_mtime


def _stale_lsjj(deck, sources):
    """The same question for LSJJ.DAT, which is allowed to be absent.

    Absent means "numpy is not installed here", not "the build failed",
    so it is *not* stale -- rebuilding would only run make_cfp again to
    the same warning.  A deck that is there and older than the generator
    that wrote it is stale, and that is what picks up an edit to
    CODE/lsjj/.
    """
    if not deck.is_file() or deck.stat().st_size == 0:
        return False
    when = deck.stat().st_mtime
    return any(when < one.stat().st_mtime
               for one in sources if one.is_file())


def lsjj_sources():
    """The generator LSJJ.DAT is built from, for the staleness check."""
    return sorted((CODE / 'lsjj').glob('*.py'))


def have_numpy():
    """True when LSJJ.DAT can be built here.

    numpy is the one thing in this package that has to be pip-installed,
    and it is optional: it buys true jj coupling and nothing else.  The
    tests that need the deck skip on this rather than fail, because an
    install without numpy is a supported install.
    """
    try:
        import numpy                                   # noqa: F401
    except ImportError:
        return False
    return True


def _make_cfp(into):
    """Build the cfp decks in `into` by running the shipped make_cfp.

    make_cfp takes the directory the decks are to go in.  It copies
    SENIOR there, and builds LSJJ.DAT there along with the three cfp
    decks -- it used to copy that one too, until it became a derived
    file like the others.  RCG is always the one in CODE/.
    It runs *in* a working directory, leaving ING11 and OUTG11 behind,
    so it is given one of its own beside the decks.

    The batch file calls `mark_time.bat` by a bare name, which the
    author's own setup finds on the PATH; CODE/ is put on the front of
    PATH here so that the call resolves without that assumption.
    """
    work = into.parent / 'cfp'
    shutil.rmtree(work, ignore_errors=True)
    work.mkdir(parents=True, exist_ok=True)
    if os.name == 'nt':
        command = ['cmd', '/c', str(CODE / 'make_cfp.bat'), str(into)]
    else:
        command = [str(CODE / 'unix' / 'make_cfp'), str(into)]
    env = dict(os.environ)
    env['PATH'] = str(CODE) + os.pathsep + env.get('PATH', '')
    env['COWAN_NO_CLEAR'] = '1'
    done = subprocess.run(command, cwd=str(work), env=env,
                          capture_output=True, text=True)
    missing = [one for one in CFP_DECKS if not (into / one).is_file()]
    if done.returncode or missing:
        pytest.skip('make_cfp could not build the cfp decks in %s: %s'
                    % (into, (done.stdout or done.stderr or '')[-2000:]))
    shutil.rmtree(work, ignore_errors=True)


@pytest.fixture(scope='session')
def cowan_code():
    """A directory of permanent files for RCG, built where a run may write.

    The programs themselves are run from CODE/ -- reading it is no
    trouble, it is writing into it that is not allowed -- and are told
    through COWAN.CFG that their permanent files are here instead.
    """
    for name in EXECUTABLES:
        if not (CODE / name).is_file():
            pytest.skip('%s is not built' % name)
    into = scratch_root()
    into.mkdir(parents=True, exist_ok=True)
    if (any(_stale_deck(into / one, CODE / RCG, CODE / 'ING11.CFP')
            for one in CFP_DECKS)
            or _stale_lsjj(into / LSJJ_DECK, lsjj_sources())):
        _make_cfp(into)
    return into


@pytest.fixture(scope='session')
def cowan_rcg(cowan_code):
    """A copy of RCG beside those files, for the tools that run it.

    calc_TP_M1E2 writes COWAN.CFG naming the directory the executable it
    was given sits in (`run_rcg`), so a tool pointed at CODE/rcg11kd
    would put the decks back in CODE/.  Pointing it at a copy beside the
    decks is the whole arrangement in one path.
    """
    copy = cowan_code / RCG
    source = CODE / RCG
    if (not copy.is_file() or copy.stat().st_size != source.stat().st_size
            or copy.stat().st_mtime < source.stat().st_mtime):
        try:
            shutil.copy2(source, copy)
        except OSError as trouble:
            pytest.skip('could not copy %s into %s: %s'
                        % (RCG, cowan_code, trouble))
    if not copy.is_file():
        pytest.skip('%s does not stay in %s -- an antivirus is removing it'
                    % (RCG, cowan_code))
    return copy
