"""Building a designation out of a coupling chain.

The rules are pinned twice over: against `ASD_syntax_examples.txt` and
the Handbook, which say what the notation is, and against the reference
runs, which say what the chains are.
"""

import pytest

from conftest import case
from toolbox import designate, genealogy, in36, ing11, outg11
from toolbox.designate import DesignateError


_CACHE = {}


def deck_of(name):
    """(shells by (kset, cfg), genealogy) of one reference run."""
    if name not in _CACHE:
        folder = case(name)
        deck = ing11.read_ing11(folder=folder)
        confs = in36.read_in36(folder=case(name.split('/')[0]))
        core = in36.common_closed_shells(confs.configurations)
        by = confs.by_label()
        shells = {(c.kset, c.index): ing11.match_shells(c, by[c.label], core)
                  for c in deck}
        _CACHE[name] = (shells, genealogy.read_genealogy(folder=folder))
    return _CACHE[name]


def label(name, coupling, kset, j, state, table=None, outer=None):
    shells, gen = deck_of(name)
    block = gen.block(kset, j)
    chain = block.state(table or designate.TABLE_OF[coupling], state)
    found = designate.designate(coupling, chain, shells[(kset, chain.cfg)],
                                outer=outer)
    return '%s | %s' % (found.configuration, found.term)


# --------------------------------------------------------------------------
# which subshells carry a term of their own
# --------------------------------------------------------------------------

def test_one_term_is_empty_closed_one_electron_or_one_hole():
    assert designate.one_term('p', 0) and designate.one_term('p', 6)
    assert designate.one_term('p', 1) and designate.one_term('p', 5)
    assert not designate.one_term('p', 2)
    assert designate.one_term('f', 13) and not designate.one_term('f', 12)
    # An s subshell has one term whatever it holds.
    assert all(designate.one_term('s', n) for n in (0, 1, 2))


def test_one_state_is_stricter_because_a_j_can_still_be_ambiguous():
    # 5d9 has one term, 2D, but two J values; in the middle of a J1J2
    # chain neither the running J nor the final term need pin it down.
    assert designate.one_term('d', 9) and not designate.one_state('d', 9)
    assert designate.one_state('s', 1) and designate.one_state('p', 6)


# --------------------------------------------------------------------------
# LS
# --------------------------------------------------------------------------

def test_a_single_open_subshell_writes_its_term_once():
    # 4f3's own term is the running term is the final term, so it is
    # written in the term column and nowhere else -- and Cowan's index
    # of a repeated term comes with it.
    assert label('Pr3', 'LS', 1, 1.5, 4) == '4f3.5p6 | 2D1*'
    assert label('Pr3', 'LS', 1, 1.5, 5) == '4f3.5p6 | 2D2*'


def test_the_last_subshell_keeps_its_own_term():
    # 1s.2p2: 1s has one term and drops it, 2p2's 3P is not the running
    # 4P and is written.  This is the reference runs' commonest label.
    assert label('U90', 'LS', 1, 1.5, 1) == '1s.2p2.(3P) | 4P'
    assert label('U90', 'LS', 1, 1.5, 2) == '1s.2p2.(1D) | 2D'


def test_a_running_term_of_two_subshells_is_written():
    # 2s.2p.3s: neither 2s nor 2p carries a term of its own, but the
    # term of the two of them together cannot be put back, so it is --
    # and it is odd where the total is odd, which is not a duplicate.
    assert label('W70', 'LS', 2, 0.5, 6) == '2s.2p.(3P*).3p | 4D'


def test_a_one_electron_subshell_in_the_middle_carries_no_term():
    # The Handbook writes 4f7(8S*) 5d(7D*) 6p: 5d's own 2D is dropped
    # and only the running term stays (`010.tex:370-374`).  The Perl
    # wrote both, which is where `4f.5p6.5d.(2D).(1D*).6s` came from.
    assert label('Pr3', 'LS', 1, 2.5, 92) == '4f.5p6.5d.(3G*).6s | 4G*'


# --------------------------------------------------------------------------
# J1J2
# --------------------------------------------------------------------------

def test_an_intermediate_j_is_written_in_angle_brackets():
    # The fix the author could not make consistently in Perl: an LS
    # term's J goes in angle brackets, a bare running J in parentheses.
    assert label('U90', 'JJ', 1, 1.5, 1) == '1s.2p2.(3P<2>) | (1/2,2)'


def test_a_subshell_term_is_kept_even_where_its_j_repeats_the_final_one():
    # `(0,3/2)` says nothing about which of 2p3's three terms this is.
    assert label('W70', 'JJ', 1, 1.5, 4) == '2p3.(2P*<3/2>) | (0,3/2)*'


def test_a_one_electron_subshell_in_the_middle_keeps_its_j():
    # 2s.2p.3s in J1J2: 2p's J is not the final term's and not the
    # running one, so `1s.2s.(1).2p-` is not the model here.
    assert (label('W70', 'JJ', 2, 0.5, 9)
            == '2s.2p.(2P*<1/2>).(1).3p | (1,3/2)')


# --------------------------------------------------------------------------
# true jj
# --------------------------------------------------------------------------

def test_a_split_subshell_is_written_as_its_halves():
    assert label('W70', 'TRUEjj', 1, 1.5, 2) == '2p+3 | (0,3/2)*'
    assert label('W70', 'TRUEjj', 1, 1.5, 4) == '2p-2.2p+ | (0,3/2)*'


def test_a_subshell_with_one_electron_is_not_split():
    # The author's rule, settled 2026-08-14: `5d9.6d | (3/2,5/2)` is
    # right as it stands, so a subshell whose j is fixed by the final
    # term is left whole.
    assert label('W70', 'TRUEjj', 1, 1.5, 1) == '2s2.2p | (0,3/2)*'


def test_the_running_k_of_two_groups_is_written():
    assert label('W70', 'TRUEjj', 1, 1.5, 8) == '2s.3s.(1).3p | (1,3/2)*'
    # A one-electron subshell in the middle IS split: nothing else says
    # which half it sits in.  Compare `2s2.2p` above, where it is last.
    assert label('W70', 'TRUEjj', 2, 0.5, 6) == '2s.2p+.(1).3p | (1,3/2)'


# --------------------------------------------------------------------------
# the JK schemes
# --------------------------------------------------------------------------

def test_lsjk_takes_its_parent_from_the_ls_chain_and_k_from_the_row():
    # `2s2.2p2.(3P<0>).5f | 2[3]*` is the shape; Pr III's is the same
    # with 4f2 as the parent.
    outer = designate.read_outer('LSJK', 'f26p        2.0K 1.0')
    text = label('Pr3/JK', 'LSJK', 1, 0.5, 3, table=genealogy.LS,
                 outer=outer)
    assert text == '4f2.(3F<2>).5p6.6p | 2[1]*'


def test_the_lsjk_row_format_is_read_as_j1_then_k():
    outer = designate.read_outer('LSJK', '4f3         0.5K 1.5')
    assert (outer.j1, outer.k) == (0.5, 1.5)


def test_lslk_prints_a_total_l_instead_of_j1():
    # FORMAT 24, `2A6,2X,A1,F4.1` (`rcg11kd.f:1805`).
    outer = designate.read_outer('LSLK', '%-14s%1s%4.1f'
                                 % ('2s2 2p2 4f', 'D', 3.0))
    assert outer.l == 'D' and outer.k == 3.0 and outer.j1 is None


def test_an_unknown_coupling_will_not_invent_a_k():
    with pytest.raises(DesignateError):
        designate.read_outer('LS', 'anything')


# --------------------------------------------------------------------------
# display order
# --------------------------------------------------------------------------

def test_a_reordered_deck_prints_in_n_l_order_where_nothing_is_observable():
    # `WORK/U90/reordered` sums 2s2 before 1s.  Both are closed or hold
    # one electron, so no intermediate is written and the order cannot
    # be seen; the label is the parent's.  The Perl printed
    # `1s.(2S).2s2`, which reads as though 1s were summed first.
    assert label('U90/reordered', 'LS', 1, 0.5, 2) == '1s.2s2 | 2S'
    assert label('U90', 'LS', 1, 0.5, 2) == '1s.2s2 | 2S'


def test_summation_order_is_kept_where_an_intermediate_makes_it_visible():
    shells, gen = deck_of('U90/reordered')
    block = gen.block(1, 1.5)
    chain = block.state(genealogy.J1J2, 1)
    found = designate.designate('JJ', chain, shells[(1, chain.cfg)])
    assert found.configuration.startswith('2p2.')
    assert not found.canonical


# --------------------------------------------------------------------------
# failure
# --------------------------------------------------------------------------

def test_a_chain_that_does_not_fit_its_card_is_an_error():
    shells, gen = deck_of('U90')
    chain = gen.block(1, 1.5).state(genealogy.LS, 1)
    with pytest.raises(DesignateError) as e:
        designate.designate('LS', chain, shells[(1, 1)])
    assert 'open subshell' in str(e.value)
