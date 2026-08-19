"""Permuting the slots of an ING11 deck, and reorder_ing11.py end to end.

The acceptance criterion here is stronger than "it agrees with the
Perl".  `WORK/U90/reordered` was built by reordering `WORK/U90` and then
*run through RCG*, so its OUTG11 says what RCG itself calls every
parameter of the reordered deck and in what order it emits them.  The
rules in `toolbox/reorder.py` are held against that, name for name.
"""

import os
import subprocess
import sys

import pytest

from conftest import case
from toolbox import ing11, outg11, reorder
from toolbox.reorder import ReorderError

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# `WORK/U90/reordered` is `WORK/U90` with its slots in the order 2, 3, 1.
U90_ORDER = '231'


def run(folder, *arguments):
    """reorder_ing11.py in a working directory, as the .bat runs it."""
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'reorder_ing11.py')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True)


# --------------------------------------------------------------------------
# the command line's permutation
# --------------------------------------------------------------------------

def test_the_digits_are_the_old_slots_in_their_new_order():
    new = reorder.read_order('231')
    assert new[2] == 1 and new[3] == 2 and new[1] == 3


def test_slots_the_string_does_not_name_stay_where_they_are():
    new = reorder.read_order('21')
    assert [new[s] for s in range(1, 9)] == [2, 1, 3, 4, 5, 6, 7, 8]


@pytest.mark.parametrize('text', ['1', '', '1', '133', '124', '21x',
                                  '123456789'])
def test_a_string_that_is_not_a_permutation_is_refused(text):
    with pytest.raises(ReorderError):
        reorder.read_order(text)


# --------------------------------------------------------------------------
# the emission order
# --------------------------------------------------------------------------

def test_a_configuration_card_permutes_only_its_first_forty_columns():
    card = 'S 2  S 1  P 0  S 0  S 0  S 0  S 0  S 0  U89+  2s           -2380426.056   0.0000'
    moved = reorder.reorder_config(card, reorder.read_order('231'))
    assert moved[:40] == 'S 1  P 0  S 2  S 0  S 0  S 0  S 0  S 0  '
    assert moved[40:] == card[40:]


def test_alpha_beta_gamma_move_with_the_f_integral_above_them():
    # PFGD emits them inside the same subshell loop as Fk(ii), so they
    # follow whichever slot that Fk(ii) has become.
    names = ['F2(11)', 'ALPHA', 'BETA', 'ZETA 1', 'F2(22)']
    new = reorder.read_order('21')
    assert reorder.permutation(names, new, is_ci=False) == [4, 0, 1, 2, 3]


def test_a_two_subshell_integral_is_named_by_its_sorted_slots():
    new = reorder.read_order('231')
    assert reorder.rename('G1(13)', new) == 'G1(23)'
    assert reorder.rename('G0(12)', new) == 'G0(13)'
    assert reorder.rename('ZETA 3', new) == 'ZETA 2'
    assert reorder.rename('ALPHA', new) == 'ALPHA'


def test_sorting_one_pair_of_a_ci_integral_exchanges_d_and_e():
    new = reorder.read_order('231')
    # 1,3 of the first configuration becomes 3,2 and has to be sorted;
    # 2,3 of the second becomes 1,2 and does not.
    assert reorder.rename('120D1323', new) == '120E2312'
    assert reorder.rename('121E1323', new) == '121D2312'


def test_a_ci_integral_with_a_repeated_slot_stays_direct():
    # Both pairings are then the same integral, and RCG calls it D.
    new = reorder.read_order('231')
    assert reorder.rename('120D1112', new) == '120D3313'


# --------------------------------------------------------------------------
# against RCG's own output
# --------------------------------------------------------------------------

def test_the_names_are_the_ones_rcg_gave_the_reordered_run():
    parent = outg11.read_parameter_names(folder=case('U90'))
    child = outg11.read_parameter_names(folder=case('U90') / 'reordered')
    new = reorder.read_order(U90_ORDER)
    assert set(parent) == set(child)
    for key, section in parent.items():
        is_ci = '-' in key
        movable = section.names if is_ci else section.names[1:]
        order = reorder.permutation(movable, new, is_ci)
        got = ([] if is_ci else ['EAV']) \
            + [reorder.rename(movable[i], new) for i in order]
        assert got == list(child[key].names), key


def test_the_deck_is_the_one_the_author_reordered(tmp_path):
    """Every field lands where `WORK/U90/reordered/ing11` has it.

    The two decks were punched from different fits -- the reference one
    from the Hartree-Fock values scaled by 0.85, this one from the
    parent's fitted parameters -- so three values differ in their last
    digits.  Everything about *where* they sit agrees.
    """
    folder = case('U90')
    out = run(folder, str(tmp_path / 'ing11.new'), U90_ORDER)
    assert out.returncode == 0, out.stderr
    with open(tmp_path / 'ing11.new') as f:
        ours = f.read().split('\n')
    with open(folder / 'reordered' / 'ing11') as f:
        theirs = f.read().split('\n')
    assert len(ours) == len(theirs)
    differ = [i for i, (a, b) in enumerate(zip(ours, theirs)) if a != b]
    assert differ == [8, 9, 11]
    for i in differ:
        assert len(ours[i]) == len(theirs[i])
        # only inside the value fields, and only in the low digits
        assert ours[i][:20] == theirs[i][:20]
        assert ours[i][70:] == theirs[i][70:]


# --------------------------------------------------------------------------
# the deck stays a deck
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name, order', [('U90', '231'), ('U90', '321'),
                                         ('W70', '231'), ('W70', '312')])
def test_the_reordered_deck_still_reads(name, order, tmp_path):
    folder = case(name)
    out = run(folder, str(tmp_path / 'ing11.new'), order)
    assert out.returncode == 0, out.stderr
    before = ing11.read_ing11(folder / 'ing11')
    after = ing11.read_ing11(tmp_path / 'ing11.new')
    assert len(before) == len(after)
    for a, b in zip(before, after):
        assert a.name == b.name
        assert sorted(a.slots) == sorted(b.slots)
        assert a.parity == b.parity
    for k in (0, 1):
        for a, b in zip(before.parameters[k], after.parameters[k]):
            assert a.np == b.np and a.factors == b.factors
            assert sorted(a.values) == pytest.approx(sorted(b.values))
        for a, b in zip(before.ci[k], after.ci[k]):
            assert a.np == b.np
            assert sorted(a.values) == pytest.approx(sorted(b.values))


def test_reordering_twice_round_returns_the_original(tmp_path):
    """The permutation and its inverse compose to nothing.

    `312` undoes `231`, so the deck comes back byte for byte -- which is
    only true because the values are moved as text rather than
    re-encoded.
    """
    folder = case('U90')
    assert run(folder, str(tmp_path / 'a'), '231').returncode == 0
    # The second pass needs an OUTG11 for the deck it is reading, and
    # the reordered run has one.
    lines = reorder.reorder_deck(
        tmp_path / 'a', reorder.read_order('312'),
        outg11_path=folder / 'reordered' / 'OUTG11')
    with open(folder / 'ing11') as f:
        original = [line.rstrip('\n').rstrip('\r') for line in f]
    assert lines == original


def test_an_outg11_that_does_not_match_the_deck_is_refused(tmp_path):
    folder = case('U90')
    with pytest.raises(ReorderError) as problem:
        reorder.reorder_deck(folder / 'ing11', reorder.read_order('231'),
                             outg11_path=case('W70') / 'OUTG11')
    assert 'OUTG11' in str(problem.value)
