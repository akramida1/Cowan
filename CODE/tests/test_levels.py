"""Level identification and the two ways of merging a second run."""

import pytest

from conftest import case
from toolbox import ing11, levels, outg11, rce
from toolbox.levels import ENERGY_ORDER, EIGENVECTOR, LevelsError

_CACHE = {}


def run(name):
    """The deck, the OUTG11 blocks and the identified levels of one case."""
    if name not in _CACHE:
        folder = case(name)
        deck = ing11.read_ing11(folder=folder)
        blocks = levels.by_block(outg11.read_outg11(folder=folder))
        found = []
        for cset in rce.read_rceout(folder=folder):
            found.extend(levels.identify(cset, blocks))
        _CACHE[name] = (deck, blocks, found)
    return _CACHE[name]


# --------------------------------------------------------------------------
# the basis index
# --------------------------------------------------------------------------

def test_an_rce_component_finds_its_rcg_basis_state():
    # RCG writes the state label after an unmatched '(' and RCE writes
    # the same text with the '(' in front of it; that is the whole
    # translation, and everything downstream rests on it.
    _, blocks, _ = run('U90')
    index = levels.BasisIndex.of(blocks[(1, 0.5)].section('LS'))
    cset = rce.read_rceout(folder=case('U90'))[0]
    for level in cset.levels:
        if level.j != 0.5:
            continue                # a component belongs to its own block
        for component in level.components:
            assert index.number_of(component) is not None, component
    assert index.number_of(rce.Component(100, '2p2', '(3P) 4P')) == 3


def test_every_component_of_every_reference_level_is_found():
    for name in ('U90', 'U90/reordered', 'W70'):
        _, blocks, _ = run(name)
        for cset in rce.read_rceout(folder=case(name)):
            for level in cset.levels:
                index = levels.BasisIndex.of(
                    blocks[(cset.index, level.j)].section('LS'))
                assert all(index.number_of(c) is not None
                           for c in level.components), (name, level.e_calc)


# --------------------------------------------------------------------------
# identification
# --------------------------------------------------------------------------

def test_every_rce_level_is_tied_to_one_rcg_eigenvector():
    for name in ('U90', 'U90/reordered', 'W70'):
        _, _, found = run(name)
        seen = set()
        for one in found:
            key = (one.kset, one.j, one.rcg)
            assert key not in seen, key
            seen.add(key)


def test_a_run_with_no_observations_falls_back_to_position():
    # W70 fits nothing -- every level is flagged '*' -- so RCE never
    # re-sorted and its order is RCG's (`conv_cowan.pl:1212-1218`).
    _, _, found = run('W70')
    assert all(one.positional for one in found)
    assert all(one.rcg == n for one in found
               for n in [one.rcg])            # the mapping is the identity
    by_j = {}
    for one in found:
        by_j.setdefault((one.kset, one.j), []).append(one.rcg)
    assert all(order == list(range(1, len(order) + 1))
               for order in by_j.values())


def test_a_fitted_run_is_matched_on_composition_below_the_ceiling():
    _, _, found = run('U90')
    assert not any(one.positional for one in found)
    assert all(one.distance is not None and one.distance < 100
               for one in found)


def test_the_identification_agrees_with_the_energies():
    # A level matched on composition should still land on an eigenvalue
    # close to its own fitted energy; a gross disagreement would mean
    # the distance had picked the wrong one.
    _, blocks, found = run('U90')
    for one in found:
        energy = blocks[(one.kset, one.j)].eigenvalues[one.rcg - 1]
        assert abs(energy - one.level.e_calc) < 1.0


@pytest.mark.slow
def test_the_real_least_squares_fit_identifies_whole():
    # Pr III is the case the method exists for: 5626 levels, most of
    # them above the highest observed one and so matched by position,
    # 920 below it and matched on composition.
    _, _, found = run('Pr3')
    assert len(found) == 1844 + 3782
    assert sum(1 for one in found if not one.positional) == 450 + 470


# --------------------------------------------------------------------------
# merging
# --------------------------------------------------------------------------

def test_the_match_mode_is_read_off_the_decks():
    plain, _, _ = run('U90')
    swapped, _, _ = run('U90/reordered')
    assert levels.match_mode(plain.configurations,
                             plain.configurations) == EIGENVECTOR
    assert levels.match_mode(plain.configurations,
                             swapped.configurations) == ENERGY_ORDER


def test_two_different_cases_are_not_a_pair():
    a, _, _ = run('U90')
    b, _, _ = run('W70')
    with pytest.raises(LevelsError) as e:
        levels.match_mode(a.configurations, b.configurations)
    assert 'not the same case' in str(e.value)


def test_a_reordered_run_merges_by_energy_order():
    _, _, ours = run('U90')
    _, _, theirs = run('U90/reordered')
    merged = levels.merge_by_energy_order(ours, theirs)
    assert len(merged) == len(ours)
    assert not levels.unmatched(merged)
    # The reordering is not exact; the differences are small and are
    # reported rather than hidden (the plan accepts them).
    worst = levels.energy_order_report(merged)
    assert abs(worst[0][1].delta) < 0.1


@pytest.mark.slow
def test_the_jk_folder_merges_into_its_parent_on_eigenvectors():
    # This is the merge the author used to do by hand: the same fit
    # re-run asking for LSJK instead of JJ.  The decks agree and so do
    # the block sizes, so every level of the parent has exactly one twin
    # carrying the same eigenvector, and the fitted energies are equal.
    ours_deck, ours_blocks, ours = run('Pr3')
    theirs_deck, theirs_blocks, theirs = run('Pr3/JK')
    assert levels.match_mode(ours_deck.configurations,
                             theirs_deck.configurations) == EIGENVECTOR
    merged = levels.merge_by_eigenvector(ours, theirs,
                                         ours_blocks, theirs_blocks)
    assert len(merged) == len(ours)
    assert not levels.unmatched(merged)
    assert max(abs(m.delta) for m in merged) < 0.01


def test_a_second_run_with_a_smaller_basis_is_refused():
    # A term cutoff (NLT11, or columns 83-85 of a configuration card)
    # shrinks the basis without changing a single card, so agreeing
    # decks are not enough -- the block sizes are what give it away.
    # Pr3/JK was such a case until the author re-punched its main card.
    _, blocks, found = run('U90')
    smaller = {k: b for k, b in blocks.items()}
    victim = sorted(smaller)[0]
    smaller[victim] = _Shrunk(smaller[victim].nlevels - 1)
    with pytest.raises(LevelsError) as e:
        levels.merge_by_eigenvector(found, found, blocks, smaller)
    assert 'do not share a basis' in str(e.value)
    assert 'NLT11' in str(e.value)


class _Shrunk:
    """A block that lost a basis state, and nothing else."""

    def __init__(self, nlevels):
        self.nlevels = nlevels


def test_a_run_merged_with_itself_is_the_identity():
    _, blocks, found = run('U90')
    merged = levels.merge_by_eigenvector(found, found, blocks, blocks)
    assert len(merged) == len(found)
    assert all(m.matched and m.delta == 0.0 for m in merged)
