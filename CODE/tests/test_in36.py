"""IN36 reading, against hand-built cards and the reference decks."""

import pytest

from conftest import case
from toolbox import in36
from toolbox.in36 import Shell


def card(text):
    """One configuration card, columns 11-80, through ANALYZ and SETCFG."""
    return in36.analyz(text)


def occ(shells):
    return {s.nl: s.occ for s in shells}


# --------------------------------------------------------------------------
# ANALYZ
# --------------------------------------------------------------------------

def test_free_layout_of_a_subshell():
    # The same subshell written four ways, as RCN's tokeniser allows.
    for text in ('Pr2+  4f3              4f3',
                 'Pr2+  4f3              4f 3',
                 'Pr2+  4f3               4f3',
                 'Pr2+  4f3              4f 03'):
        conf, nlbcd8, wwnl8, _, _, _, _ = card(text)
        assert nlbcd8[0] == ' 4f', text
        _, _, frac = in36.setcfg(57, nlbcd8, wwnl8)
        assert not frac


def test_orbital_letter_is_case_folded():
    _, nlbcd8, wwnl8, _, _, _, _ = card('Pr2+  4f3              4F3')
    assert nlbcd8[0] == ' 4f'


def test_two_digit_principal_quantum_number():
    _, nlbcd8, wwnl8, _, _, _, _ = card('B IV  1s10s            1s   10s')
    assert nlbcd8[:2] == [' 1s', '10s']


def test_l_letters_skip_j():
    # 1s9l: l is the ninth letter of RCN's table, so l = 8.
    _, nlbcd8, _, _, _, _, _ = card('B IV  1s9l             1s   9l')
    assert nlbcd8[1] == ' 9l'
    assert in36.LBCD.index('l') == 8


def test_configuration_name_and_spectrum():
    conf, _, _, _, _, _, _ = card('Pr2+  f26p             4f2  5p6')
    assert conf[0:6].strip() == 'Pr2+'
    assert conf[6:18].strip() == 'f26p'
    # A spectrum name may hold a blank; the scan starts past it.
    conf, _, _, _, _, _, _ = card('B IV  1s2              1s2')
    assert conf[0:6].strip() == 'B IV'
    assert conf[6:18].strip() == '1s2'


def test_fractional_occupation():
    # The occupation has to touch the orbital letter: ANALYZ takes at most
    # six columns after it, and never across a blank.
    _, nlbcd8, wwnl8, _, _, _, _ = card('Xx    4f65             4f6.500')
    assert (nlbcd8[0], wwnl8[0]) == (' 4f', '6.500 ')
    shells, _, frac = in36.setcfg(57, nlbcd8, wwnl8)
    assert frac
    assert occ(shells)['4f'] == 6.5


def test_a_detached_number_is_not_an_occupation():
    # '4f 6.500' is two tokens: a 4f with one electron, and a continuum
    # energy.  This is why Pr3's '4f   5p6' means 4f1.
    _, nlbcd8, wwnl8, _, ee8, _, _ = card('Xx    4f               4f 6.500')
    shells, _, _ = in36.setcfg(57, nlbcd8, wwnl8)
    assert occ(shells)['4f'] == 1.0
    assert ee8 == pytest.approx(6.5)


def test_kut_and_continuum_energy():
    _, _, _, kut, ee8, _, _ = card('B IV  1s2              1s2  -3')
    assert kut == -3
    _, _, _, _, ee8, _, _ = card('B IV  1sep             1s   ep     0.250')
    assert ee8 == pytest.approx(0.25)


def test_alfmin_alfmax():
    _, _, _, _, _, alfmin, alfmax = card('B IV  1s2      2065     1s2')
    assert alfmin == pytest.approx(0.20)
    assert alfmax == pytest.approx(0.65)


# --------------------------------------------------------------------------
# SETCFG
# --------------------------------------------------------------------------

def test_noble_gas_core_sizes():
    assert len(in36._noble_core(2)) == 1                  # He
    assert len(in36._noble_core(57)) == 11                # Xe
    assert sum(s[2] for s in in36._noble_core(57)) == 54
    assert sum(s[2] for s in in36._noble_core(86)) == 86  # Rn
    assert in36._noble_core(3)[-1][:2] == [1, 0]          # He core only
    # The Xe core has no 4f; the Rn core has no 5f and no 5g.
    assert (4, 3) not in [(n, l) for n, l, _ in in36._noble_core(57)]
    assert (5, 3) not in [(n, l) for n, l, _ in in36._noble_core(87)]


def test_zero_occupation_deletes_a_core_subshell():
    # U89+ '2p2' is 1s1 2s0 2p2: the 2s0 is not there to be deleted, but
    # a 1s0 would delete the core 1s.
    _, nlbcd8, wwnl8, _, _, _, _ = card('U89+  2p2             1s1  2s0  2p2')
    shells, _, _ = in36.setcfg(3, nlbcd8, wwnl8)
    assert occ(shells) == {'1s': 1.0, '2p': 2.0}
    _, nlbcd8, wwnl8, _, _, _, _ = card('Xx    2s2             1s0  2s2')
    shells, _, _ = in36.setcfg(2, nlbcd8, wwnl8)
    assert occ(shells) == {'2s': 2.0}


def test_l_not_less_than_n_is_rejected():
    _, nlbcd8, wwnl8, _, _, _, _ = card('Xx    1p              1p')
    with pytest.raises(in36.In36Error):
        in36.setcfg(2, nlbcd8, wwnl8)


# --------------------------------------------------------------------------
# whole decks
# --------------------------------------------------------------------------

def test_u90_deck():
    deck = in36.read_in36('IN36', case('U90'))
    confs = deck.by_label()
    assert sorted(confs) == ['2p', '2p2', '2s', '2s2', '2s2p']
    for c in deck:
        assert c.noelec == 3
        assert sum(s.occ for s in c.shells) == 3

    # The label the Perl loses entirely: 1s1 2s2, not '1s'.
    assert str(confs['2s2']) == '1s 2s2'
    assert str(confs['2s']) == '1s2 2s'
    assert str(confs['2p2']) == '1s 2p2'
    assert str(confs['2s2p']) == '1s 2s 2p'

    assert confs['2s'].parity == 0
    assert confs['2p'].parity == 1
    assert in36.parities(deck.configurations) == {0, 1}


def test_w70_deck():
    deck = in36.read_in36('IN36', case('W70'))
    confs = deck.by_label()
    for c in deck:
        assert c.noelec == 5
        assert sum(s.occ for s in c.shells) == 5
    assert str(confs['s22p']) == '1s2 2s2 2p'
    assert str(confs['2p3']) == '1s2 2p3'
    assert str(confs['sp4f']) == '1s2 2s 2p 4f'
    assert confs['s22p'].parity == 1
    assert confs['sp4f'].parity == 0       # 2p and 4f together are even


def test_pr3_deck():
    deck = in36.read_in36('IN36', case('Pr3'))
    confs = deck.by_label()
    for c in deck:
        assert c.noelec == 57
        assert sum(s.occ for s in c.shells) == 57
    assert str(confs['4f3']).endswith('4d10 4f3 5s2 5p6')
    assert str(confs['f5d2']).endswith('4d10 4f 5s2 5p6 5d2')
    assert confs['4f3'].parity == 1
    assert confs['f5d2'].parity == 1       # 4f 5d2
    assert confs['f26p'].parity == 1        # 4f2 6p
    assert confs['f26s'].parity == 0        # 4f2 6s
    assert confs['f4p5'].parity == 1        # 4f4 5p5
    assert in36.parities(deck.configurations) == {0, 1}


def test_common_closed_shells_of_pr3():
    # 5p is closed in most of Pr III but open in the three 5p5
    # configurations carried for CI, so it is not droppable by default.
    deck = in36.read_in36('IN36', case('Pr3'))
    shared = in36.common_closed_shells(deck.configurations)
    assert [s.nl for s in shared] == \
        ['1s', '2s', '2p', '3s', '3p', '3d', '4s', '4p', '4d', '5s']
    assert '5p' not in [s.nl for s in shared]


def test_open_shells():
    deck = in36.read_in36('IN36', case('Pr3'))
    confs = deck.by_label()
    assert [s.nl for s in confs['4f3'].open_shells] == ['4f']
    assert [s.nl for s in confs['f5d2'].open_shells] == ['4f', '5d']
