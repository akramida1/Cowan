"""read_aa.py and the OUTG11 Auger sections behind it.

There is no autoionization run among the three reference cases, so the
data here is `WORK/AI/DIEL/OUTG11.reference`, the B IV example the
author ships with the package.  The plan rules that case out as
*validation* material and it is not used as such: nothing is diffed
against a stored expectation.  What is checked is that the rates read
back are the ones printed in the file, which the file itself settles.

`read_aa.pl` cannot be the reference here.  It stops at the first J that
does not increase (`read_aa.pl:104`), which is the start of the second
parity set -- and in this run the Auger rates are all in the second
set.  It therefore reads a set that has none, matches stray numbers as
eigenvalues, and prints two rows of zeros.  That is why `read_aa.py`
takes a `--set`.
"""

import os
import subprocess
import sys

import pytest

from conftest import WORK
from toolbox import outg11

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture(scope='module')
def diel(tmp_path_factory):
    """The B IV dielectronic example, as a working directory."""
    import shutil
    source = WORK / 'AI' / 'DIEL' / 'OUTG11.reference'
    if not source.is_file():
        pytest.skip('WORK/AI/DIEL is not present')
    folder = tmp_path_factory.mktemp('diel')
    shutil.copy(str(source), str(folder / 'OUTG11'))
    return folder


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'read_aa.py')] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


def table(out):
    return [line.split('\t') for line in out.strip().split('\n')]


# --------------------------------------------------------------------------
# the reader
# --------------------------------------------------------------------------

def test_the_blocks_that_carry_rates_are_the_second_parity_set(diel):
    blocks = outg11.read_auger(folder=str(diel))
    assert blocks
    assert set(b.kset for b in blocks) == {2}
    assert sorted(b.j for b in blocks) == [4.0, 5.0, 6.0, 7.0, 8.0]


def test_a_channel_has_one_rate_for_every_eigenvalue(diel):
    for block in outg11.read_auger(folder=str(diel)):
        for channel, rates in block.rates.items():
            assert len(rates) == block.nlevels, (block.j, channel)


def test_the_channel_label_is_the_target_and_its_term(diel):
    blocks = dict((b.j, b) for b in outg11.read_auger(folder=str(diel)))
    assert sorted(blocks[6.0].rates) == ['18.58h (2S) 3H', '18.58k (2S) 3K']


def test_the_rates_are_the_ones_printed(diel):
    """`FOR/rcg11kd.f:11373`, the J=6 block of the example."""
    blocks = dict((b.j, b) for b in outg11.read_auger(folder=str(diel)))
    assert blocks[6.0].rates['18.58k (2S) 3K'] == [
        0.0, 0.0, 5.51e+01, 1.02e+05, 4.89e+04, 2.03e-01]
    assert blocks[6.0].rates['18.58h (2S) 3H'] == [0.0] * 6


# --------------------------------------------------------------------------
# the table
# --------------------------------------------------------------------------

def test_the_table_has_one_column_per_core(diel):
    out = run(diel, '18.58h:0,18.58k:0', '2039.0', '--set', '2')
    assert out.returncode == 0, out.stderr
    rows = table(out.stdout)
    assert rows[0] == ['J', 'Energy', '18.58h:2039', '18.58k:2039', 'AA_tot']
    assert all(len(row) == 5 for row in rows)


def test_a_core_pattern_sums_the_terms_it_matches(diel):
    """The two K terms of the J=7 block add to what one column shows."""
    out = run(diel, '18.58k:0', '2039.0', '--set', '2')
    assert out.returncode == 0, out.stderr
    seven = [row for row in table(out.stdout) if row[0] == '7.0']
    assert [float(row[2]) for row in seven] == [1.02e+05, 4.89e+04, 1.51e+05]


def test_the_placeholder_levels_are_left_out(diel):
    """RCG parks the continuum states at -9500; they are not levels."""
    out = run(diel, '18.58k:0', '2039.0', '--set', '2')
    assert all(float(row[1]) > -5000 for row in table(out.stdout)[1:])


def test_a_channel_the_level_cannot_reach_is_closed(diel):
    """Above the core's energy the rate stands; below it there is none."""
    high = run(diel, '18.58k:0', '2039.0', '--set', '2')
    shut = run(diel, '18.58k:0', '9999.0', '--set', '2')
    assert any(float(row[2]) for row in table(high.stdout)[1:])
    assert not any(float(row[2]) for row in table(shut.stdout)[1:])


def test_the_total_is_the_sum_of_the_columns(diel):
    out = run(diel, '18.58h:0,18.58k:0', '2039.0', '--set', '2')
    for row in table(out.stdout)[1:]:
        assert float(row[4]) == pytest.approx(float(row[2]) + float(row[3]),
                                              rel=1e-2)


def test_a_set_with_no_rates_says_so_and_names_the_one_that_has_them(diel):
    out = run(diel, '18.58k:0', '2039.0', '--set', '1')
    assert out.returncode == 1
    assert 'parity set 1' in out.stderr and 'them in 2' in out.stderr


def test_a_core_that_matches_nothing_lists_what_there_is(diel):
    out = run(diel, '1s2:0', '2039.0', '--set', '2')
    assert out.returncode == 1
    assert "'1s2'" in out.stderr and '18.58k (2S) 3K' in out.stderr


def test_the_cores_need_an_energy_each(diel):
    assert run(diel, '18.58k', '2039.0', '--set', '2').returncode == 1
    assert run(diel, '18.58k:0', 'later').returncode == 2
