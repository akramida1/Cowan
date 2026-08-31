"""LSJJ.DAT: that it is reproducible, and that make_cfp is what builds it.

The deck is a derived file since 2026 -- not committed, not distributed,
built on the installed machine by make_cfp along with FOR072, FOR073 and
FOR074.  These tests are what stands in for the copy that used to be in
the repository to diff against.

THERE IS DELIBERATELY NO HARD-CODED md5 HERE.  There was going to be
one: `CODE/lsjj/lsjj.py`'s tolerant `fix_sign` tie-break (SIGN_TOL) was
written to make the deck the same bytes on every machine, and it did
remove the whole-vector sign flips that used to differ between one
LAPACK build and the next.  What it cannot remove is the rest: the
coefficients are null spaces and degenerate eigenspaces returned by
`numpy.linalg.svd` and `eigh`, and two LAPACK builds disagree about a
few dozen of the ~58000 doubles in the last bit.  Measured on this deck,
Windows against WSL2: 28 doubles differ, none by more than 2.3e-16, and
not one of them is a sign flip or a genuinely different value.  Nothing
RCG prints moves with them -- on WORK/W70, the true jj case, TAPE2E is
byte-identical and OUTG11 differs only in a timing line.

So an md5 would be a tripwire that fires on every machine but the one it
was recorded on, which is worse than none.  What is checkable, and is
checked below, is that the build is deterministic *here*: the same
inputs give the same bytes, twice running, and the deck make_cfp left in
the scratch directory is the one mkdeck.py builds today.
"""

import filecmp
import os
import subprocess
import sys

import pytest

from conftest import CODE, LSJJ_DECK, have_numpy


needs_numpy = pytest.mark.skipif(
    not have_numpy(),
    reason='numpy is not installed, so LSJJ.DAT cannot be built here; '
           'that is a supported install -- it costs true jj coupling '
           'and nothing else')


def build_deck(into):
    """Run mkdeck.py the way make_cfp runs it, and return the path."""
    into.mkdir(parents=True, exist_ok=True)
    target = into / 'LSJJ.DAT'
    env = dict(os.environ)
    env['LSJJ_CFP'] = str(CODE / 'ING11.CFP')
    env['LSJJ_SENIOR'] = str(CODE / 'SENIOR')
    done = subprocess.run(
        [sys.executable, '-B', str(CODE / 'lsjj' / 'mkdeck.py'), str(target)],
        capture_output=True, text=True, env=env)
    assert done.returncode == 0, done.stdout + done.stderr
    assert target.is_file() and target.stat().st_size > 0
    return target


@needs_numpy
def test_the_deck_is_reproducible(tmp_path):
    """Twice from the same inputs, byte for byte.

    This is the half of reproducibility that does hold everywhere: it was
    demonstrated under WSL2 as well as on Windows.  It is what says the
    generator carries no dictionary ordering, no hash seed and no clock
    into the file.
    """
    one = build_deck(tmp_path / 'a')
    two = build_deck(tmp_path / 'b')
    assert filecmp.cmp(str(one), str(two), shallow=False)


@needs_numpy
def test_make_cfp_builds_the_deck_that_mkdeck_builds(cowan_code, tmp_path):
    """The deck in the scratch directory is the current generator's.

    `cowan_code` is the directory make_cfp filled for the whole test
    session, so this is the shipped route -- make_cfp.bat or
    unix/make_cfp -- checked against mkdeck.py run directly.  If they
    ever disagree it is the wiring in the wrapper that has drifted, not
    the algebra.
    """
    theirs = cowan_code / LSJJ_DECK
    if not theirs.is_file():
        pytest.skip('make_cfp did not build %s in %s' % (LSJJ_DECK,
                                                         cowan_code))
    mine = build_deck(tmp_path)
    assert filecmp.cmp(str(theirs), str(mine), shallow=False)


@needs_numpy
def test_the_deck_carries_the_magic_rcg_expects(tmp_path):
    """First record: magic, version, and the two record counts.

    LSJJIN refuses a deck whose first record does not say
    ILSMAG=20260801, ILSVER=1 (`FOR/rcg11kd.f`), which is what stops a
    stale deck, one from the other byte order, or one still in the old
    text form from being read as garbage.  Bump one of the two and this
    test says which file the other one lives in.
    """
    import struct

    deck = build_deck(tmp_path)
    with open(str(deck), 'rb') as fh:
        head = fh.read(24)
    length, magic, version, nlrec, njrec = struct.unpack('<5i', head[:20])
    assert length == 16
    assert magic == 20260801
    assert version == 1
    assert (nlrec, njrec) == (34, 24)
