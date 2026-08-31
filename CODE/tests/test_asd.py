"""The ASD notation, pinned to the two documents that define it.

`ASD_syntax_examples.txt` is reproduced line for line: the table below
builds each designation out of `asd` objects, and the file is then read
and compared against it.  The Handbook cases that the sample file does
not reach -- J1K, LK, and the jj form of an equivalent-electron group --
follow.
"""

import pytest

from conftest import REPO
from toolbox import asd
from toolbox.asd import (Intermediate, JTerm, KTerm, LSTerm, Level, PairTerm,
                         Subshell)


def s(n, l, occ=1, sign=''):
    return Subshell(n, l, occ, sign)


# --------------------------------------------------------------------------
# ASD_syntax_examples.txt, in the order of the file
# --------------------------------------------------------------------------

EXAMPLES = [
    # Li-like
    Level([s(1, 's', 2), s(2, 's')], LSTerm(2, 'S'), 0.5),
    Level([s(1, 's', 2), s(2, 'p')], LSTerm(2, 'P', odd=True), 0.5),
    Level([s(1, 's'), s(2, 's', 2)], PairTerm(0.5, 0), 0.5),
    Level([s(1, 's'), s(2, 's'), Intermediate(JTerm(1)), s(2, 'p', 1, '-')],
          PairTerm(1, 0.5, odd=True), 1.5),
    Level([s(1, 's'), s(2, 's'), Intermediate(JTerm(0)), s(2, 'p', 1, '-')],
          PairTerm(0, 0.5, odd=True), 0.5),
    Level([s(1, 's'), s(2, 'p', 2, '-')], PairTerm(0.5, 0), 0.5),
    Level([s(1, 's'), s(2, 's'), Intermediate(JTerm(1)), s(2, 'p', 1, '+')],
          PairTerm(1, 1.5, odd=True), 2.5),
    Level([s(1, 's'), s(2, 'p', 1, '-'), Intermediate(JTerm(1)),
           s(2, 'p', 1, '+')], PairTerm(1, 1.5), 2.5),
    Level([s(1, 's'), s(2, 'p', 1, '-'), Intermediate(JTerm(1)),
           s(2, 'p', 1, '+')], PairTerm(1, 1.5), 1.5),
    Level([s(1, 's'), s(2, 's'), Intermediate(JTerm(0)), s(2, 'p', 1, '+')],
          PairTerm(0, 1.5, odd=True), 1.5),
    Level([s(1, 's'), s(2, 'p', 1, '-'), Intermediate(JTerm(0)),
           s(2, 'p', 1, '+')], PairTerm(0, 1.5), 1.5),
    # No trailing intermediate: the J of the last group is already the
    # first half of the final term, so writing it again says nothing.
    Level([s(1, 's'), s(2, 'p', 2, '+')], PairTerm(0.5, 2), 2.5),
    Level([s(1, 's'), s(2, 'p', 2, '+')], PairTerm(0.5, 0), 0.5),
    # Pr III
    Level([s(4, 'f', 3)], LSTerm(2, 'L', odd=True), 7.5),
    Level([s(4, 'f', 3)], LSTerm(2, 'D', index='2', odd=True), 1.5),
    Level([s(4, 'f', 2), Intermediate(LSTerm(3, 'H', j=4)), s(6, 's')],
          PairTerm(4, 0.5), 3.5),
    # LSLK: the total-L letter is part of the term, not of the
    # configuration (author's decision, 2026-08-14).
    Level([s(2, 's', 2), s(2, 'p', 2), Intermediate(LSTerm(3, 'P')),
           s(4, 'f')], KTerm(2, 3, odd=True, l='D'), 2.5),
    # LSJK
    Level([s(2, 's', 2), s(2, 'p', 2), Intermediate(LSTerm(3, 'P', j=0)),
           s(5, 'f')], KTerm(2, 3, odd=True), 2.5),
]

def file_rows():
    text = (REPO / 'ASD_syntax_examples.txt').read_text()
    rows = []
    for line in text.splitlines()[1:]:            # past the column heading
        if line.count('|') == 2:
            rows.append(tuple(part.strip() for part in line.split('|')))
    return rows


def test_every_line_of_the_sample_is_reproduced():
    rows = file_rows()
    assert len(rows) == len(EXAMPLES)
    for row, level in zip(rows, EXAMPLES):
        assert level.columns() == row


def test_the_sample_covers_what_it_is_meant_to():
    # A guard on the reader above: if the file grows a section, this
    # count changes and the table has to grow with it.
    assert len(file_rows()) == 18


# --------------------------------------------------------------------------
# the rules one at a time
# --------------------------------------------------------------------------

def test_an_s_subshell_never_takes_a_sign():
    # j = 1/2 is the only value, so no sign is written (010.tex:511-513).
    # The readers hand '+' over for every s subshell, so this is dropped
    # here rather than being the caller's business.
    assert s(1, 's', 1, '+').text == '1s'
    assert s(2, 's', 2, '+').text == '2s2'
    assert s(2, 'p', 1, '+').text == '2p+'
    assert s(2, 'p', 2, '-').text == '2p-2'


def test_occupation_one_is_not_written():
    assert s(6, 's').text == '6s'
    assert s(4, 'f', 13).text == '4f13'


def test_l_may_be_given_as_a_number_or_a_letter():
    assert s(4, 3, 3).text == s(4, 'f', 3).text == '4f3'
    assert LSTerm(3, 5).text == LSTerm(3, 'H').text == '3H'


def test_parity_is_marked_on_the_term_and_only_there():
    tokens = [s(2, 's', 2), s(2, 'p', 2), Intermediate(LSTerm(3, 'P')),
              s(4, 'f')]
    assert asd.parity(t for t in tokens if isinstance(t, Subshell)) == 1
    assert asd.configuration(tokens) == '2s2.2p2.(3P).4f'
    # The L letter of LK carries no mark of its own even where the level
    # is odd; the mark goes on the bracket, once.
    assert KTerm(2, 3, odd=True, l='D').text == 'D 2[3]*'
    assert KTerm(2, 3, odd=True).text == '2[3]*'


def test_an_intermediate_ls_term_gives_its_j_in_angle_brackets():
    # This is the JJ fix: the Perl wrote (2P*(1/2)).
    inter = Intermediate(LSTerm(2, 'P', odd=True, j=0.5))
    assert inter.text == '(2P*<1/2>)'
    assert Intermediate(LSTerm(3, 'H', j=4)).text == '(3H<4>)'
    # An intermediate that is itself a J keeps plain parentheses.
    assert Intermediate(JTerm(0.5)).text == '(1/2)'
    assert Intermediate(JTerm(2)).text == '(2)'


def test_the_nielson_koster_index_follows_the_letter():
    assert LSTerm(2, 'D', index='2', odd=True).text == '2D2*'
    # With a parent J as well, the order is letter, index, parity, J.
    # No line of the sample file shows this combination.
    assert LSTerm(2, 'H', index='2', odd=True, j=4.5).text == '2H2*<9/2>'


def test_seniority_rides_on_the_subshell_that_owns_it():
    # Nothing in the sample or the Handbook writes a seniority inline;
    # this is the convention the author settled on 2026-08-14.  The
    # anchor is the subshell because a middle subshell's own J is never
    # written by a running-sum coupling and so cannot carry one.
    assert s(4, 'f', 4, '+').text == '4f+4'
    assert Subshell(4, 'f', 4, '+', v=2).text == '4f+4<2>'
    assert Subshell(4, 'f', 4, '-', v=0).text == '4f-4<0>'


def test_a_seniority_of_none_writes_nothing():
    assert JTerm(2).text == '2'
    assert PairTerm(0.5, 2).text == '(1/2,2)'
    assert Subshell(4, 'f', 4, '+').text == '4f+4'
    with pytest.raises(asd.AsdError):
        Subshell(4, 'f', 4, '+', v=-1).text


def test_a_space_may_set_off_a_later_coupled_term():
    # 010.tex:393-395: a space marks a term formed only after those to
    # its left.  Cowan's schemes are all running sums and never need it.
    tokens = [s(4, 'f', 7), Intermediate(LSTerm(8, 'S', odd=True)),
              s(5, 'd'), Intermediate(LSTerm(7, 'D', odd=True), space=True),
              s(6, 'p')]
    assert asd.configuration(tokens) == '4f7.(8S*).5d (7D*).6p'


def test_j_is_written_whole_or_halved():
    assert asd.format_j(0) == '0'
    assert asd.format_j(4.0) == '4'
    assert asd.format_j(0.5) == '1/2'
    assert asd.format_j(7.5) == '15/2'
    assert asd.format_j2(0) == '0'
    assert asd.format_j2(15) == '15/2'
    assert asd.format_j2(8) == '4'


def test_a_j_that_is_neither_whole_nor_half_is_an_error():
    for bad in (0.3, -0.5, 'x'):
        with pytest.raises(asd.AsdError):
            asd.format_j(bad)


def test_a_bad_sign_or_letter_is_an_error():
    with pytest.raises(asd.AsdError):
        Subshell(2, 'p', 1, '*')
    with pytest.raises(asd.AsdError):
        Subshell(2, 'j').text
    with pytest.raises(asd.AsdError):
        LSTerm(2, 'J').text
    with pytest.raises(asd.AsdError):
        asd.configuration(['4f3'])


# --------------------------------------------------------------------------
# Handbook cases the sample file does not reach
# --------------------------------------------------------------------------

def test_jj_coupling_of_equivalent_electrons():
    # 010.tex:451, example 5: 4d3(5/2) 4d2(3/2) (9/2,2)_11/2, written
    # with the + and - signs of example 5 rather than the j values.
    level = Level([s(4, 'd', 3, '+'), s(4, 'd', 2, '-')],
                  PairTerm(4.5, 2), 5.5)
    assert level.columns() == ('4d+3.4d-2', '(9/2,2)', '11/2')


def test_j1j2_with_an_ls_parent():
    # 010.tex:481, example 1: 3d9(2D_5/2) 4p_3/2 (5/2,3/2)*_3.
    level = Level([s(3, 'd', 9), Intermediate(LSTerm(2, 'D', j=2.5)),
                   s(4, 'p', 1, '+')], PairTerm(2.5, 1.5, odd=True), 3)
    assert level.columns() == ('3d9.(2D<5/2>).4p+', '(5/2,3/2)*', '3')


def test_j1k_coupling():
    # 010.tex:523, example 1: 3p5(2P*_1/2) 5g  2[9/2]*_5.
    level = Level([s(3, 'p', 5), Intermediate(LSTerm(2, 'P', odd=True, j=0.5)),
                   s(5, 'g')], KTerm(2, 4.5, odd=True), 5)
    assert level.columns() == ('3p5.(2P*<1/2>).5g', '2[9/2]*', '5')


def test_j1l2_coupling_with_an_external_term():
    # 010.tex:525, example 3: 4f13(2F*_7/2) 5d2(1D) 1[7/2]*_7/2.  The
    # external term is an ordinary intermediate; only the final term
    # carries the bracket.
    level = Level([s(4, 'f', 13), Intermediate(LSTerm(2, 'F', odd=True, j=3.5)),
                   s(5, 'd', 2), Intermediate(LSTerm(1, 'D'))],
                  KTerm(1, 3.5, odd=True), 3.5)
    assert level.columns() == ('4f13.(2F*<7/2>).5d2.(1D)', '1[7/2]*', '7/2')


def test_lk_coupling():
    # 010.tex:550, example 1: 3s2 3p(2P*) 4f G  2[7/2]_3, with the G
    # moved into the term column.
    level = Level([s(3, 's', 2), s(3, 'p'),
                   Intermediate(LSTerm(2, 'P', odd=True)), s(4, 'f')],
                  KTerm(2, 3.5, l='G'), 3)
    assert level.columns() == ('3s2.3p.(2P*).4f', 'G 2[7/2]', '3')
