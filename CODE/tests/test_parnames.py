"""Resolving RCE's positional parameter names against the deck."""

import pytest

from conftest import case
from toolbox import in36, ing11, parnames, rce
from toolbox.parnames import ParNameError


def resolved(name):
    """{parity set: [Named]} for one reference run."""
    deck = ing11.read_ing11(folder=case(name))
    confs = in36.read_in36(folder=case(name.split('/')[0]))
    core = in36.common_closed_shells(confs.configurations)
    by = confs.by_label()
    shells = {c.label: ing11.match_shells(c, by[c.label], core) for c in deck}
    out = {}
    for cset in rce.read_rceout(folder=case(name)):
        out[cset.index] = parnames.resolve(cset.parameters,
                                           deck.set(cset.index), shells)
    return out


def texts(name, kset):
    return [n.text for n in resolved(name)[kset]]


# --------------------------------------------------------------------------
# the names one kind at a time
# --------------------------------------------------------------------------

def test_every_kind_of_name_is_resolved():
    names = texts('U90', 1)
    assert names[:7] == ['Eav(2s)', 'Eav(2s2)', 'Eav(2p2)', 'F2(2p,2p)',
                         'ALPHA(2p)', 'zeta(2p)', 'G1(1s,2p)']
    # A CI integral names two subshells from each configuration.
    assert names[7:] == ['RD0(1s1s,1s2s)', 'RD0(1s2s,2s2s)',
                         'RD1(1s2s,2p2p)', 'RD1(2s2s,2p2p)']


def test_alpha_beta_and_gamma_take_the_shell_of_the_f_before_them():
    # This is what makes them resolvable: PFGD emits them inside the
    # same per-subshell loop as Fk(ii) (`FOR/rcg11kd.f:4383-4400`).
    names = texts('Pr3', 1)
    assert names[:8] == ['Eav(4f3)', 'F2(4f,4f)', 'F4(4f,4f)', 'F6(4f,4f)',
                         'ALPHA(4f)', 'BETA(4f)', 'GAMMA(4f)', 'zeta(4f)']


def test_zeta_carries_its_own_slot():
    # W70's s24f has one open subshell besides the closed 2s2, and RCE
    # writes 'ZETA 3' for it whatever slot 3 happens to hold.
    names = texts('W70', 1)
    assert 'zeta(4f)' in names and 'zeta(3p)' in names


def test_every_parameter_of_every_reference_run_resolves():
    for name in ('U90', 'U90/reordered', 'W70', 'Pr3', 'Pr3/JK'):
        for kset, named in resolved(name).items():
            assert all(n.text for n in named), (name, kset)


# --------------------------------------------------------------------------
# reordering
# --------------------------------------------------------------------------

def test_a_reordered_deck_gives_the_same_single_configuration_names():
    # The slot numbers differ -- 'ZETA 3' becomes 'ZETA 2' -- but the
    # subshell behind them does not, so the readable name is the same.
    plain = {n.raw: n.text for n in resolved('U90')[2]}
    swapped = {n.raw: n.text for n in resolved('U90/reordered')[2]}
    assert plain['ZETA 3'] == swapped['ZETA 2'] == 'zeta(2p)'
    assert plain['G1(13)'] == swapped['G1(12)'] == 'G1(1s,2p)'
    assert plain['G0(12)'] == swapped['G0(13)'] == 'G0(1s,2s)'


def test_the_reordered_deck_names_different_ci_integrals():
    # Not a naming failure: reorder_ing11.pl does not permute the slot
    # digits of the CI cards, so the reordered deck really does carry
    # different integrals.  The plan accepts that -- a reordered run is
    # used only for percentage compositions -- and this pins it so that
    # it is a recorded fact rather than a surprise.
    plain = sorted(n.text for n in resolved('U90')[2] if n.kind == 'R')
    swapped = sorted(n.text for n in resolved('U90/reordered')[2]
                     if n.kind == 'R')
    assert plain == ['RD0(1s1s,1s2s)', 'RD0(1s2p,2s2p)', 'RE1(1s2p,2s2p)']
    assert swapped == ['RD0(1s1s,1s2s)', 'RD1(1s2p,1s2p)', 'RE0(1s2p,1s2p)']


# --------------------------------------------------------------------------
# ordering and failure
# --------------------------------------------------------------------------

def test_the_sort_key_groups_a_configuration_s_parameters_together():
    named = resolved('Pr3')[1]
    keys = [n.key for n in named[:8]]
    assert keys == sorted(keys)
    # Sorting by the key leaves the average energies in configuration
    # order, which is the order the report lists the blocks in.
    eavs = [n for n in sorted(named, key=lambda x: x.key) if n.kind == 'Eav']
    assert [n.configurations[0] for n in eavs[:3]] == ['4f3', 'f26p', 'f27p']


def test_a_name_that_belongs_to_no_subshell_is_named_not_guessed():
    deck = ing11.read_ing11(folder=case('U90'))
    confs = in36.read_in36(folder=case('U90'))
    core = in36.common_closed_shells(confs.configurations)
    by = confs.by_label()
    shells = {c.label: ing11.match_shells(c, by[c.label], core) for c in deck}
    bad = [rce.Parameter(name='EAV 2s', flag=0, value=0.0),
           rce.Parameter(name='ZETA 8', flag=0, value=0.0)]
    with pytest.raises(ParNameError) as e:
        parnames.resolve(bad, deck.set(1), shells)
    assert 'slot 8' in str(e.value)


def test_a_parameter_before_any_eav_is_an_error():
    deck = ing11.read_ing11(folder=case('U90'))
    with pytest.raises(ParNameError):
        parnames.resolve([rce.Parameter(name='ZETA 3', flag=0, value=0.0)],
                         deck.set(1), {c.label: [] for c in deck})


# --------------------------------------------------------------------------
# a CI name with no configuration columns beside it
# --------------------------------------------------------------------------

def test_a_ci_name_decodes_its_own_two_configurations():
    # OUTGINE has no columns 58-71, so the two code characters are what
    # is left.  RCG writes '0'-'9' for 0-9, 'A'-'Z' for 10-35 and one
    # character per number above that (`FOR/rcg11kd.f:6308-6322`).
    order = ['c%d' % i for i in range(1, 40)]
    match = parnames._CI.match('1K0D1112')
    assert parnames._ci_codes(match, order) == ('c1', 'c20')
    assert parnames._ci_codes(parnames._CI.match('Aa0D1112'), order) \
        == ('c10', 'c36')


def test_a_ci_name_the_parity_set_cannot_hold_is_not_decoded():
    # '*' stands for every serial number past 188, and a number the set
    # does not have is a name from somewhere else; either way the caller
    # is told nothing rather than the wrong thing.
    order = ['c1', 'c2']
    assert parnames._ci_codes(parnames._CI.match('*20D1112'), order) is None
    assert parnames._ci_codes(parnames._CI.match('190D1112'), order) is None


def test_the_ci_integrals_of_outgine_match_the_ones_rceout_names():
    # The same run read both ways: RCEOUT's own columns, and OUTGINE's
    # two characters.  They must name the same integrals.
    folder = case('Pr3/JK')
    deck = ing11.read_ing11(folder=folder)
    confs = in36.read_in36(folder=folder.parent)      # JK/ has no IN36
    core = in36.common_closed_shells(confs.configurations)
    by = confs.by_label()
    shells = {c.label: ing11.match_shells(c, by[c.label], core) for c in deck}
    starting = rce.read_outgine(folder=folder)
    named = parnames.resolve(starting[1], deck.set(1), shells)
    integrals = [n for n in named if n.kind == 'R']
    assert integrals and all(len(n.configurations) == 2 for n in integrals)
    assert integrals[0].text == 'RD2(4f4f,4f6p)'
    assert integrals[0].configurations == ('4f3', 'f26p')
