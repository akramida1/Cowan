"""The ING11 reader, against hand-built cards and the reference decks."""

import pytest

from conftest import case
from toolbox import files, in36, ing11
from toolbox.ing11 import ConfigCard, Ing11Error


def deck(name):
    return ing11.read_ing11(folder=case(name))


def raw(name):
    """The deck's own lines, trailing blanks and all."""
    return files.read_lines(files.resolve(files.ING11, case(name)),
                            strip=False)


def in36_of(name):
    """The IN36 of a run; a sub-folder deck shares its parent's."""
    return in36.read_in36(folder=case(name.split('/')[0]))


# --------------------------------------------------------------------------
# the main card
# --------------------------------------------------------------------------

W70_MAIN = ('    1   -1412124 8 80         001111181 00000000 0 1000.0000 '
            '10.0   1 02-6 0 0     2       2 12 12  8  8')


def test_the_main_card_is_read_field_by_field():
    m = ing11.read_main(W70_MAIN)
    assert m.kcpl == 1
    assert m.nocset == -1
    assert m.nsconf == ((4, 12, 12), (4, 8, 8))
    assert (m.nconfigurations(1), m.nconfigurations(2)) == (12, 8)
    assert m.nsets == 2
    assert m.iabg == 0            # W70 uses no ALPHA/BETA/GAMMA
    assert m.kcpld == (0, 0, 1, 1, 1, 1, 1, 8, 1)
    assert m.uenrgy == 1000.0
    assert m.ipgen == 2


def test_kcpld_says_yes_with_a_zero_and_true_jj_with_an_eight():
    # KCPLD(I)=0 asks for coupling I (`rcg11kd.f:786-790`); column 8 is
    # inverted so that a blank card does not switch true jj on.
    m = ing11.read_main(W70_MAIN)
    assert m.true_jj
    assert m.couplings == ['LS', 'JJ', 'TRUEjj']
    # Couplings 3-7 are gated by NLSMAX as well.  W70 leaves it blank,
    # which since 2026 reads as 999 and permits them; W70 asks for none
    # of them anyway, its KCPLD digits 3-7 all being 1.
    assert m.nlsmax == 999 and m.nlspch == 0


def test_a_hundred_configurations_are_read_from_the_overflow_field():
    # NSCONF(2,K) is a 2-digit field; RCN2 punches ** and puts all four
    # counts in columns 93-104 as 4I3 (`rcg11kd.f:701-710`).
    card = ('    1   -18**  8**  2').ljust(92) + '120120090090'
    m = ing11.read_main(card)
    assert m.nsconf == ((8, 120, 120), (8, 90, 90))


# --------------------------------------------------------------------------
# configuration cards
# --------------------------------------------------------------------------

PR3_CARD = ('F 3  P 6  D 0  S 0  P 0  D 0  S 0  P 0  Pr2+  4f3     '
            '     -2026070.397   0.0000')


def test_a_configuration_card_is_read_slot_by_slot():
    c = ing11.read_config_card(PR3_CARD, 1, 1)
    assert c.slots[:3] == [('F', 3), ('P', 6), ('D', 0)]
    assert (c.spectrum, c.label) == ('Pr2+', '4f3')
    assert c.eav == pytest.approx(-2026070.397)
    assert c.ekin == 0.0
    assert c.nlt == 0                   # blank field: use NLT11
    assert c.parity == 1                # f3 p6 -> 3*3 + 1*6 odd


def test_a_short_card_reads_as_no_term_cutoff():
    # Cowan's card ends at column 80 and RCN2 still writes exactly that,
    # so the A3 in 83-85 is off the end of the record (`:3326-3332`).
    assert ing11.read_config_card(PR3_CARD[:80], 1, 1).nlt == 0
    assert ing11.read_config_card(PR3_CARD[:80].ljust(84) + '9', 1, 1).nlt == 9


def test_the_electron_count_on_a_card_is_not_the_spectrum_s():
    # RCN2 strips the closed core the run shares, so W70's s22p shows 3.
    d = deck('W70')
    assert d.by_label()['s22p'].nelectrons == 3


# --------------------------------------------------------------------------
# whole decks
# --------------------------------------------------------------------------

def test_every_reference_deck_reads_whole():
    for name, nconf in (('W70', 20), ('Pr3', 39), ('U90', 5),
                        ('U90/reordered', 5), ('Pr3/JK', 39)):
        d = deck(name)
        assert len(d) == nconf
        assert len(d.configurations) == sum(d.main.nconfigurations(k)
                                            for k in (1, 2))
        for k in (1, 2):
            if d.main.nconfigurations(k):
                assert len(d.parameters[k - 1]) == d.main.nconfigurations(k)


def test_the_prelude_is_kept_but_not_interpreted():
    # Pr3 carries Cowan's KCPL=3 control card and a KCPL=0 rescale card
    # in front of the main card (`rcg11kd.f:485`, `:726`).
    d = deck('Pr3')
    assert len(d.prelude) == 2
    assert d.prelude[0].startswith('    3')
    assert d.main.kcpl == 1


def test_parameter_and_ci_cards_carry_their_own_count():
    d = deck('U90')
    par = {p.label: p for p in d.parameters[0]}
    assert par['2s'].np == 1 and len(par['2s'].values) == 1
    assert par['2p2'].np == 5 and len(par['2p2'].values) == 5
    assert par['2p2'].kind == 'hr'
    ci = {c.labels: c for c in d.ci[0]}
    assert set(ci) == {('2s', '2s2'), ('2s', '2p2'), ('2s2', '2p2')}
    assert ci[('2s', '2s2')].np == 2


def test_a_card_with_more_than_five_parameters_is_continued():
    # Seven (F9.4,I1) pairs per continuation line (`rcg11kd.f:10834`).
    d = deck('Pr3')
    par = {p.label: p for p in d.parameters[0]}
    assert par['f25f'].np == len(par['f25f'].values) == 22
    assert len(par['f25f'].raw) == 4


def test_the_trailer_is_kept_verbatim():
    # Cowan's sentinels are not this reader's business, but a writer has
    # to be able to put them back.
    d = deck('U90')
    assert d.trailer[-1].strip() == '-1'
    assert d.trailer[0].strip() == '-55555555.'


def test_the_radial_integral_cards_are_found_wherever_they_sit():
    # RCN2 writes a cross-parity block after the last parameter card of
    # the deck and a same-parity block in the middle of it, so they are
    # kept on the card they follow rather than in the trailer.
    d = deck('U90')
    cards = ing11.read_r_integrals(d)
    assert cards and all(one.order == 1 for one in cards)
    assert cards[0].raw.startswith('U89+  2s')
    assert '//R1//' in cards[0].raw
    assert not ing11.read_r_integrals(d, order=2)


def test_a_radial_integral_card_round_trips_when_its_value_is_put_back():
    d = deck('U90')
    for one in ing11.read_r_integrals(d):
        assert ing11.set_r_value(one, one.value) == one.raw
    assert ing11.write_ing11(d) == raw('U90')


def test_the_jk_deck_asks_for_lsjk_instead_of_jj():
    m = deck('Pr3/JK').main
    assert m.couplings == ['LS', 'LSJK']
    assert not m.true_jj


def test_a_blank_nlsmax_permits_couplings_3_to_7_and_a_punched_zero_does_not():
    # `rcg11kd.f:897-910` reads columns 23-25 as A3: blank is 999, "no
    # limit", and only a punched 0 switches couplings 3-7 off.  The card
    # here asks for LS and LSJK, as the JK deck does.
    card = ('    1   -182121818182 ' + '   '
            + '     0111011 1 00000000 0 1000.0000')
    blank = card
    assert ing11.read_main(blank).couplings == ['LS', 'LSJK']
    assert ing11.read_main(blank).nlsmax_effective is None
    zero = card[:22] + '  0' + card[25:]
    assert ing11.read_main(zero).couplings == ['LS']
    assert ing11.read_main(zero).nlsmax_effective == 0


def test_nlsmax_of_999_and_nocset_both_mean_no_limit():
    # `rcg11kd.f:756-762`.  Either 999 or a non-zero NLSMAX with NOCSET
    # set lifts the limit; the JK deck has both.
    m = deck('Pr3/JK').main
    assert m.nocset != 0 and m.nlsmax == 999
    assert m.nlsmax_effective is None
    assert deck('Pr3').main.nlsmax_effective == 0


def test_the_jk_deck_keeps_the_whole_term_list():
    # The JK deck once carried 999 in columns 24-26 instead of 23-25, so
    # RCG read NLSMAX=99 and NLT11=9 and cut 4f3 from 17 LS terms to 9.
    # The author re-punched the card and re-ran it; NLT11 back to 0 is
    # what says the two folders share a basis, which is what the
    # eigenvector merge in levels.py needs.
    assert deck('Pr3/JK').main.nlt11 == 0
    assert deck('Pr3').main.nlt11 == 0


# --------------------------------------------------------------------------
# tying a slot to an nl
# --------------------------------------------------------------------------

def matched(name):
    d = deck(name)
    confs = in36_of(name)
    core = in36.common_closed_shells(confs.configurations)
    by = confs.by_label()
    return {c.label: ing11.match_shells(c, by[c.label], core) for c in d}


def test_every_card_of_every_deck_finds_its_n():
    for name in ('W70', 'Pr3', 'U90', 'U90/reordered'):
        shells = matched(name)
        assert len(shells) == len(deck(name))


def test_the_reordered_deck_keeps_its_own_summation_order():
    # This is the bug: U90/reordered sums 2s before 1s, and the label has
    # to say so.  Matching is on the (l, occupation) multiset, so the
    # slots come back permuted rather than renumbered.
    plain = matched('U90')
    swapped = matched('U90/reordered')
    assert [str(s) for s in plain['2s2']] == ['1s', '2s2']
    assert [str(s) for s in swapped['2s2']] == ['2s2', '1s']
    assert [(s.slot, s.n) for s in swapped['2s2']] == [(1, 2), (3, 1)]
    # Same subshells either way, only the order differs.
    assert sorted(str(s) for s in plain['2s2']) == \
        sorted(str(s) for s in swapped['2s2'])


def test_the_shared_core_is_not_on_the_card():
    # Pr III shares [Xe]-less 1s2..5s2; the cards begin at 4f.
    shells = matched('Pr3')
    assert [str(s) for s in shells['4f3']] == ['4f3', '5p6']
    assert [str(s) for s in shells['f26p']] == ['4f2', '5p6', '6p']


def test_a_card_that_does_not_match_its_configuration_is_named():
    conf = in36_of('U90').by_label()['2s2']
    card = ConfigCard(kset=1, index=1, name='U89+  2s2         ',
                      slots=[('S', 1), ('P', 2)] + [('S', 0)] * 6)
    with pytest.raises(Ing11Error) as e:
        ing11.match_shells(card, conf)
    assert '2s2' in str(e.value) and 'p2' in str(e.value)
# --------------------------------------------------------------------------
# the cuts that censor the printed line list
# --------------------------------------------------------------------------

def test_dmin_is_read_off_the_main_card():
    # Columns 61-65, `E5.1`; W70 punches ` 10.0`, a gA cut at 1e10.
    assert ing11.read_main(W70_MAIN).dmin == 10.0
    assert ing11.read_main(W70_MAIN[:60] + ' 8.00' + W70_MAIN[65:]).dmin == 8.0
    assert ing11.read_main(W70_MAIN[:60] + '     ' + W70_MAIN[65:]).dmin == 0.0


def test_the_window_of_a_hand_built_kcpl3_card():
    # `2E5.0` in columns 11-20 since 2026: a punched integer is angstrom.
    d = ing11.Ing11(prelude=['    3 0  0  80016000'])
    assert d.window == (800.0, 16000.0)
    ing11.set_window(d, 761.9, 16800.0)
    assert d.window == (761.9, 16800.0)
    assert d.prelude[0][:10] == '    3 0  0'
    # A blank or impossible field takes RCG own defaults at label 95.
    assert ing11.Ing11(prelude=['    3']).window == (
        ing11.FLBDN_DEFAULT, ing11.FLBDX_DEFAULT)


def test_dmin_and_the_wavelength_window_are_read_off_the_deck():
    d = deck('Pr3/JK')
    assert d.main.dmin == 10.0                  # main card columns 61-65
    assert d.window == (800.0, 16000.0)         # KCPL=3 card, columns 11-20


def test_a_deck_with_no_kcpl3_card_reports_RCG_own_window():
    d = deck('U90')
    assert d.control(3) is None
    assert d.window == (ing11.FLBDN_DEFAULT, ing11.FLBDX_DEFAULT)
    with pytest.raises(Ing11Error):
        ing11.set_window(d, 100.0, 200.0)


def test_dmin_above_one_is_an_exponent_and_below_one_is_not():
    # `FOR/rcg11kd.f:13294-13302`: only the integer part is taken, and
    # the sign is the level-sort flag rather than part of the cut.
    assert ing11.dmin_cut(8.0) == ('gA', 1e8)
    assert ing11.dmin_cut(8.7) == ('gA', 1e8)
    assert ing11.dmin_cut(-8.0) == ('gA', 1e8)
    assert ing11.dmin_cut(0.5) == ('D0', 0.5)
    assert ing11.dmin_cut(1.0) == ('D0', 1.0)
    assert ing11.dmin_cut(60.0) == ('D0', 60.0)


def test_a_new_dmin_and_window_read_back_as_they_were_written(tmp_path):
    d = deck('Pr3/JK')
    ing11.set_dmin(d, 9.0)
    ing11.set_window(d, 800.0 / 1.05, 16000.0 * 1.05)
    files.write_deck(tmp_path / files.ING11, ing11.write_ing11(d), chr(10))
    back = ing11.read_ing11(folder=tmp_path)
    assert back.main.dmin == 9.0
    low, high = back.window
    assert 761.0 <= low <= 762.0 and high == 16800.0
    # Nothing else on either card moved.
    was = deck('Pr3/JK')
    assert back.main.raw[:60] == was.main.raw[:60]
    assert back.main.raw[65:] == was.main.raw[65:]
    assert back.prelude[0][:10] == was.prelude[0][:10]
    assert back.prelude[0][20:] == was.prelude[0][20:]


def test_writing_the_same_dmin_back_changes_nothing():
    d = deck('Pr3/JK')
    before = d.main.raw
    assert ing11.set_dmin(d, d.main.dmin) == before
