"""add_eav --ground: OUTGINE with the calculated ground level at zero.

The decks here are built record by record with the FORMATs RCG writes
OUTGINE with (`FOR/rcg11kd.f:15740-15893`), so they run without a
reference case.  The last test takes the OUTGINE of every reference run
that has one, and skips where the run has not been made.
"""

import os
import subprocess
import sys

import pytest

from conftest import case
from toolbox import files, rce

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CONVERGENCE = ('       1.0       0.0      -1.0 0.00 1000     1.000  10000.00'
               '-0.30 0.00-1')


def f4(x):
    return '%10.4f' % x


def section(kset, configs, names, values, blocks, decimals=4, invect=0):
    """One section of OUTGINE as RCG writes it; `values` are CVPAR text."""
    out = ['%5d%5d%5d%5d%5d%5d%5d' % (1, 30, 5, -6, invect, 0, 0),
           '%5d%5d%5d%5d%10.3f%5d%5d%10.4f' % (kset, 0, 1, 1, 0.0, 1, 3, 1.0),
           ('Xx0+      ' + ''.join('%-10s' % c for c in configs)).ljust(70)]
    for at in range(0, len(names), 7):
        out.append(''.join('%-10s' % n for n in names[at:at + 7]))
    for energies in blocks:
        for at in range(0, len(energies), 7):
            out.append(''.join('%10.*f' % (decimals, e)
                               for e in energies[at:at + 7]))
        numbers = list(range(1, len(energies) + 1))
        for at in range(0, len(numbers), 7):
            out.append(''.join('%10d' % n for n in numbers[at:at + 7]))
    for block in (['%10d' % 100] * len(names), [f4(0.0)] * len(names),
                  values):
        for at in range(0, len(block), 7):
            out.append(''.join(block[at:at + 7]))
    return out + [CONVERGENCE, '   -1']


# Pb I 6p2, one configuration, as RCG wrote it for F:\COWAN\Pb-like.
PB = section(1, ['6p2'], ['EAV 6p2', 'F2(11)', 'ALPHA', 'ZETA 1'],
             [f4(0.0), f4(23.7096), '  0.00E+00', f4(7.3892)],
             [[-14.4602, 15.6065], [-6.5398], [-3.4721, 7.1667]])


def write(folder, lines, name=rce.OUTGINE):
    files.write_deck(os.path.join(str(folder), name), lines, '\r\n')


def ground(folder, *extra):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'add_eav.py'), '--ground']
        + list(extra), cwd=str(folder), capture_output=True, text=True)


def lines_in(folder, name=rce.OUTGINE):
    return files.read_lines(os.path.join(str(folder), name), strip=False)


def heading(kset, j, n):
    """A genealogy.txt block heading, FORMAT 11 of GENPRT."""
    return ' J =%5.1f    configuration set%2d    basis states:%7d' % (j, kset,
                                                                      n)


# --------------------------------------------------------------------------
# the shift
# --------------------------------------------------------------------------

def test_the_ground_goes_to_zero_and_every_eav_and_level_with_it(tmp_path):
    write(tmp_path, PB)
    out = ground(tmp_path)
    assert out.returncode == 0, out.stderr
    assert 'E_RCG = -14.4602' in out.stdout
    assert '14.4602 added to 1 Eav parameters and 5 levels' in out.stdout
    new = lines_in(tmp_path)
    assert new[4:10] == ['    0.0000   30.0667', '         1         2',
                         '    7.9204', '         1',
                         '   10.9881   21.6269', '         1         2']
    assert new[12] == '   14.4602   23.7096  0.00E+00    7.3892'
    # Everything else is the characters it was.
    changed = [i for i, (a, b) in enumerate(zip(PB, new)) if a != b]
    assert changed == [4, 6, 8, 12]
    assert lines_in(tmp_path, 'OUTGINE.RCG') == PB


def test_a_second_ground_shift_is_zero_and_writes_nothing(tmp_path):
    write(tmp_path, PB)
    assert ground(tmp_path).returncode == 0
    once = lines_in(tmp_path)
    out = ground(tmp_path)
    assert out.returncode == 0
    assert 'already at zero' in out.stdout
    assert lines_in(tmp_path) == once
    assert lines_in(tmp_path, 'OUTGINE.RCG') == PB     # still RCG's file


def test_a_shift_of_zero_gives_back_the_file():
    assert rce.shift_outgine(PB, 0.0)[0] == PB


def test_the_ground_of_the_second_set_moves_the_first_as_well(tmp_path):
    first = section(1, ['a'], ['EAV a', 'ZETA 1'], [f4(5.0), f4(1.0)],
                    [[4.0, 6.0]])
    second = section(2, ['b'], ['EAV b', 'ZETA 1'], [f4(1.0), f4(1.0)],
                     [[-2.5, 3.0]])
    write(tmp_path, first + second)
    out = ground(tmp_path)
    assert 'configuration set 2, E_RCG = -2.5000' in out.stdout
    s1, s2 = rce.read_outgine(folder=tmp_path).values()
    assert (s1[0].value, s2[0].value) == (7.5, 3.5)
    low = rce.outgine_ground(lines_in(tmp_path))
    assert (low.cset, low.energy) == (2, 0.0)


# --------------------------------------------------------------------------
# when nothing is, or can be, shifted
# --------------------------------------------------------------------------

def test_invect_not_zero_is_a_warning_and_shifts_nothing(tmp_path):
    deck = list(PB)
    deck[0] = '%5d%5d%5d%5d%5d%5d%5d' % (1, 30, 5, -6, 1, 0, 0)
    write(tmp_path, deck)
    out = ground(tmp_path)
    assert out.returncode == 0
    assert 'warning' in out.stderr and 'INVECT=1' in out.stderr
    assert lines_in(tmp_path) == deck
    assert not (tmp_path / 'OUTGINE.RCG').exists()


def test_an_outgine_with_every_level_zero_is_an_error(tmp_path):
    # RCG writes zeros for both when it has made no diagonalization.
    write(tmp_path, section(1, ['6p2'], ['EAV 6p2', 'ZETA 1'],
                            [f4(0.0), f4(0.0)], [[0.0, 0.0], [0.0]]))
    out = ground(tmp_path)
    assert out.returncode == 1
    assert 'no diagonalization' in out.stderr


# --------------------------------------------------------------------------
# continuum configurations, and levels that only look like them
# --------------------------------------------------------------------------

def test_a_continuum_eav_and_its_levels_stay_where_they_are(tmp_path):
    # RCN2 punches -9500 for a continuum configuration of an
    # autoionization run, and its levels come out near that value.
    deck = section(1, ['d', 'c'], ['EAV d', 'ZETA 1', 'EAV c'],
                   [f4(2.0), f4(1.0), f4(-9500.0)],
                   [[-9512.25, -1.5, 3.0]])
    low = rce.outgine_ground(deck)
    assert (low.index, low.energy) == (2, -1.5)
    new, neav, nlev = rce.shift_outgine(deck, 1.5)
    assert (neav, nlev) == (1, 2)
    write(tmp_path, new)
    (params,) = rce.read_outgine(folder=tmp_path).values()
    assert [p.value for p in params] == [3.5, 1.0, -9500.0]
    # The continuum Eav counts in GRE as it does in RCG: |-9500| > 7000
    # asks for F10.3.
    assert new[4] == ' -9512.250     0.000     4.500'


def test_a_level_below_minus_4000_is_a_level_without_a_continuum_eav():
    # W70: W69+, every Eav >= 0, ground level -8391.10 from the spin-orbit
    # splitting of 2s2 2p.  Nothing in that section is a continuum.
    deck = section(1, ['s22p'], ['EAV s22p', 'ZETA 2'],
                   [f4(0.0), '  12586.65'], [[-8391.10], [4195.55]],
                   decimals=2)
    low = rce.outgine_ground(deck)
    assert low.energy == -8391.10 and low.text == '-8391.10'
    new, neav, nlev = rce.shift_outgine(deck, 8391.10)
    assert (neav, nlev) == (1, 2)
    # GRE is now the 12586.65 of ZETA 2, which asks for F10.3.
    assert new[4] == '     0.000' and new[6] == ' 12586.650'


# --------------------------------------------------------------------------
# encoding: the file stays one RCG could have written
# --------------------------------------------------------------------------

def test_the_level_format_follows_the_largest_parameter_after_the_shift():
    # GRE = 6990 asks for F10.4; shifted to 7010 it asks for F10.3.
    deck = section(1, ['a'], ['EAV a', 'ZETA 1'], [f4(6990.0), f4(1.0)],
                   [[-20.0, 1.23456]])
    new = rce.shift_outgine(deck, 20.0)[0]
    assert new[4] == '     0.000    21.235'


def test_an_eav_near_zero_is_written_1pe10_2_unless_it_is_the_first():
    deck = section(1, ['a', 'b'], ['EAV a', 'EAV b'],
                   [f4(-1.0), f4(-1.0004)], [[-1.0, 0.0]])
    new = rce.shift_outgine(deck, 1.0)[0]
    # `FOR/rcg11kd.f:15841-15848`: the exponent form for a small value,
    # but never for the first parameter of a section.
    assert new[-3] == '    0.0000 -4.00E-04'


# --------------------------------------------------------------------------
# naming the J of the ground level
# --------------------------------------------------------------------------

def test_the_j_of_the_ground_comes_from_genealogy(tmp_path):
    write(tmp_path, PB)
    files.write_deck(tmp_path / files.GENEALOGY,
                     [' Genealogy of the basis states', '',
                      heading(1, 0.0, 2), '', heading(1, 1.0, 1), '',
                      heading(1, 2.0, 2)], '\r\n')
    out = ground(tmp_path)
    assert ('Ground level adopted: level 1 with J = 0 of configuration '
            'set 1, E_RCG = -14.4602') in out.stdout


def test_a_genealogy_that_does_not_match_is_not_believed(tmp_path):
    write(tmp_path, PB)
    files.write_deck(tmp_path / files.GENEALOGY,
                     [heading(1, 0.0, 3), heading(1, 1.0, 1),
                      heading(1, 2.0, 2)], '\r\n')
    out = ground(tmp_path)
    assert 'level 1 of J block 1 of configuration set 1' in out.stdout


def test_without_genealogy_the_block_is_named_by_position(tmp_path):
    write(tmp_path, PB)
    out = ground(tmp_path)
    assert 'level 1 of J block 1 of configuration set 1' in out.stdout


# --------------------------------------------------------------------------
# the reference runs
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name', ['U90', 'W70', 'Pr3', 'Mo15'])
def test_a_reference_outgine_moves_only_its_eav_and_its_levels(name,
                                                               tmp_path):
    path = files.resolve(rce.OUTGINE, case(name), required=False)
    if path is None:
        pytest.skip('%s has no OUTGINE' % name)
    deck = files.read_lines(path, strip=False)
    assert rce.shift_outgine(deck, 0.0)[0] == deck
    write(tmp_path, deck)
    out = ground(tmp_path)
    assert out.returncode == 0, out.stderr
    if rce.outgine_invect(deck) != 0:
        assert 'warning' in out.stderr and lines_in(tmp_path) == deck
        return
    low = rce.outgine_ground(deck)
    before = rce.read_outgine(path=path)
    after = rce.read_outgine(folder=tmp_path)
    for k in before:
        for p, q in zip(before[k], after[k]):
            if p.name.startswith('EAV') and p.value > rce.CONTINUUM:
                # F10.2 past 99999.99 is RCG's own precision there.
                assert q.value == pytest.approx(p.value - low.energy,
                                                abs=0.005)
            else:
                assert q.value == p.value
    assert rce.outgine_ground(lines_in(tmp_path)).energy == 0.0
    assert 'already at zero' in ground(tmp_path).stdout
