"""Pairing two runs' eigenvectors by overlap, and the multipole table.

`toolbox/overlap.py` is what replaced a per-trial `levels.identify` in
`calc_TP`, and the thing it has to get right is a swap: two neighbouring
eigenvectors of one J block that change places when the parameters move.
The blocks here are built by hand, so the answer is known.
"""

import math

import pytest

import calc_TP
from toolbox import outg11, radiate

numpy = pytest.importorskip('numpy')

from toolbox import overlap                                   # noqa: E402


# --------------------------------------------------------------------------
# blocks built by hand
# --------------------------------------------------------------------------

class Row:
    def __init__(self, number, amplitudes):
        self.number = number
        self.amplitudes = list(amplitudes)
        self.label = 'c%d' % number
        self.bracketed = 's%d' % number


class Section:
    def __init__(self, rows):
        self.rows = rows


def block(columns, kset=1, j=2.0):
    """A JBlock whose LS eigenvectors are the given columns."""
    size = len(columns)
    rows = [Row(s + 1, [columns[k][s] for k in range(size)])
            for s in range(len(columns[0]))]
    return outg11.JBlock(j=j, kset=kset,
                         eigenvalues=[0.0] * size,
                         sections={'LS': Section(rows)})


def rotate(angle):
    """Two eigenvectors mixed by `angle`, which at pi/4 is a full swap."""
    c, s = math.cos(angle), math.sin(angle)
    return [[c, s], [-s, c]]


# --------------------------------------------------------------------------

def test_an_unmoved_block_pairs_one_to_one():
    one = block([[1.0, 0.0], [0.0, 1.0]])
    pairs, unpaired = overlap.pair_block(one, one)
    assert unpaired == 0
    assert sorted((a, b) for a, b, _ in pairs) == [(1, 1), (2, 2)]
    assert all(value == pytest.approx(1.0) for _, _, value in pairs)


def test_two_eigenvectors_that_changed_places_are_followed():
    base = block([[1.0, 0.0], [0.0, 1.0]])
    trial = block([[0.0, 1.0], [1.0, 0.0]])
    pairs, unpaired = overlap.pair_block(base, trial)
    assert unpaired == 0
    assert sorted((a, b) for a, b, _ in pairs) == [(1, 2), (2, 1)]


def test_a_small_rotation_does_not_swap_them():
    base = block([[1.0, 0.0], [0.0, 1.0]])
    trial = block(rotate(0.2))
    pairs, _ = overlap.pair_block(base, trial)
    assert sorted((a, b) for a, b, _ in pairs) == [(1, 1), (2, 2)]
    assert min(value for _, _, value in pairs) == pytest.approx(math.cos(0.2))


def test_the_sign_of_an_eigenvector_does_not_matter():
    # RCG fixes no phase, so a trial may print -v where the base run
    # printed v.  The overlap is taken in absolute value for that.
    base = block([[1.0, 0.0], [0.0, 1.0]])
    trial = block([[-1.0, 0.0], [0.0, -1.0]])
    pairs, _ = overlap.pair_block(base, trial)
    assert sorted((a, b) for a, b, _ in pairs) == [(1, 1), (2, 2)]


def test_remap_carries_the_identification_over():
    base = {(1, 2.0): block([[1.0, 0.0], [0.0, 1.0]])}
    trial = {(1, 2.0): block([[0.0, 1.0], [1.0, 0.0]])}
    mapping = {(1, 2.0, 1): 'first', (1, 2.0, 2): 'second'}
    out, report = overlap.remap(mapping, base, trial)
    assert out == {(1, 2.0, 2): 'first', (1, 2.0, 1): 'second'}
    assert report.paired == 2
    assert report.below(0.7) == 0


def test_a_trial_with_a_different_basis_stops_the_run():
    base = {(1, 2.0): block([[1.0, 0.0], [0.0, 1.0]])}
    trial = {(1, 2.0): block([[1.0]])}
    with pytest.raises(overlap.OverlapError) as problem:
        overlap.remap({}, base, trial)
    assert 'J=2' in str(problem.value)


def test_a_missing_block_stops_the_run():
    base = {(1, 2.0): block([[1.0, 0.0], [0.0, 1.0]])}
    with pytest.raises(overlap.OverlapError) as problem:
        overlap.remap({}, base, {})
    assert 'missing from the trial' in str(problem.value)


def test_the_report_reads_back():
    # Half-and-half is the worst a 2x2 rotation can do: both overlaps
    # come out at 1/sqrt(2), which is the point at which the pairing
    # stops meaning anything and the report has to say so.
    base = {(1, 2.0): block([[1.0, 0.0], [0.0, 1.0]])}
    trial = {(1, 2.0): block(rotate(math.pi / 4))}
    _, report = overlap.remap({}, base, trial)
    assert report.paired == 2
    assert report.median == pytest.approx(0.5 ** 0.5)
    assert report.below(0.8) == 2
    assert '0 below 0.3' in str(report)


# --------------------------------------------------------------------------
# the multipole table calc_TP dispatches on
# --------------------------------------------------------------------------

def test_every_multipole_has_a_spectrum_name_rcg_knows():
    for name, section in calc_TP.SPECTRUM_OF.items():
        assert radiate.OF_SPECTRUM[section] == name


def test_only_e1_crosses_the_parity_sets():
    crossing = [name for name in calc_TP.SPECTRUM_OF
                if name not in calc_TP.SAME_PARITY]
    assert crossing == [radiate.E1]


def test_the_multipole_list_is_parsed_in_the_order_written():
    assert calc_TP.parse_multipoles('E2, m1 ,E1') == ('E2', 'M1', 'E1')
    assert calc_TP.parse_multipoles('e1') == ('E1',)
    assert calc_TP.parse_multipoles('E1 E1 M1') == ('E1', 'M1')


def test_an_unknown_multipole_is_named_in_the_message():
    with pytest.raises(calc_TP.RunError) as problem:
        calc_TP.parse_multipoles('E1,M2')
    assert 'M2' in str(problem.value)
    assert 'E1, M1, E2' in str(problem.value)

    with pytest.raises(calc_TP.RunError):
        calc_TP.parse_multipoles('  ')
