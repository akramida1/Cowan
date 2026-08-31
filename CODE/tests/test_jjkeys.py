"""true-jj_keys.txt: the alphabet, the label layouts, and the W70 key."""

import pytest

from conftest import case
from toolbox import jjkeys


# --------------------------------------------------------------------------
# the alphabet
# --------------------------------------------------------------------------

def test_alphabet_is_the_one_lsjjcc_builds():
    # 1..9, then A..Z, then a..z, then 24 punctuation marks: 85 codes.
    # '0' is left out so that a code is never taken for a blank or a
    # zero, and '*' is reserved for "no code left".
    assert len(jjkeys.ALPHABET) == 85
    assert jjkeys.ALPHABET[:9] == '123456789'
    assert jjkeys.ALPHABET[9] == 'A' and jjkeys.ALPHABET[34] == 'Z'
    assert jjkeys.ALPHABET[35] == 'a' and jjkeys.ALPHABET[60] == 'z'
    assert '0' not in jjkeys.ALPHABET
    assert '*' not in jjkeys.ALPHABET
    assert ' ' not in jjkeys.ALPHABET
    assert '(' not in jjkeys.ALPHABET and ')' not in jjkeys.ALPHABET


def test_code_values_count_from_one():
    # A J character reads 1 for 2J=0, so the doubled J is value - 1.
    assert jjkeys.code_value('1') == 1
    assert jjkeys.code_value('9') == 9
    assert jjkeys.code_value('A') == 10
    assert jjkeys.code_value('a') == 36
    assert jjkeys.code_value('z') == 61
    assert jjkeys.code_value('#') == 62
    assert jjkeys.code_value(jjkeys.ALPHABET[-1]) == 85


def test_the_overflow_mark_is_reported_not_decoded():
    with pytest.raises(jjkeys.JJKeyError):
        jjkeys.code_value('*')


# --------------------------------------------------------------------------
# the label layouts
# --------------------------------------------------------------------------

def test_layouts():
    kinds = lambda n: ''.join(k for k, _ in jjkeys.layout(n))
    assert kinds(1) == 'AJ'
    assert kinds(2) == 'AJAJ'
    assert kinds(3) == 'AJAJKAJ'
    # With four subshells the Ki no longer fit and are dropped.
    assert kinds(4) == 'AJAJAJAJ'
    # With five or more the Ji go as well, leaving the Ai alone.
    assert kinds(5) == 'AAAAA'
    assert len(jjkeys.layout(9)) == 8


def test_decoding_a_three_subshell_label():
    # W70's state 5 of J=1/2, set 1: s+1 s+1 p+1 with 2J = 1, 1, 3 and
    # 2K2 = 2, which genealogy.txt gives independently.
    got = jjkeys.decode('1212314', 3)
    assert [(p.kind, p.value) for p in got] == [
        ('A', 0), ('J', 1), ('A', 0), ('J', 1), ('K', 2), ('A', 0), ('J', 3)]
    assert [p.subshell for p in got] == [1, 1, 2, 2, 2, 3, 3]


def test_a_short_label_is_an_error_not_a_guess():
    with pytest.raises(jjkeys.JJKeyError):
        jjkeys.decode('12', 3)


# --------------------------------------------------------------------------
# the file
# --------------------------------------------------------------------------

def keys():
    return jjkeys.read_jjkeys(case('W70') / 'true-jj_keys.txt')


def test_w70_key_covers_every_rotated_subshell():
    k = keys()
    assert len(k) == 41
    assert len(k.configs()) == 20
    assert 'W69+  2p3' in k.configs()


def test_a_subshell_with_several_states():
    # 2p3 has four jj states, each with its own code and partition.
    sub = keys().subshell('W69+  2p3', 2)
    assert (sub.l, sub.nel) == ('p', 3)
    assert sorted(sub.codes) == ['1', '2', '3', '4']
    parts = sub.codes['1'].partition()
    assert [(p.sign, p.nel, p.j) for p in parts] == [('-', 1, 0.5),
                                                     ('+', 2, 0.0)]
    assert [(p.sign, p.nel, p.j)
            for p in sub.codes['4'].partition()] == [('-', 2, 0.0),
                                                     ('+', 1, 1.5)]
    # A code that carries one partition needs no J to be asked for.
    assert sub.codes['2'].by_j == ()


def test_an_unrotated_subshell_carries_both_values_of_j():
    # One electron: LSJJRT leaves the subshell alone, so it has one LS
    # term and one code, and that code covers both j- and j+.  The Ji
    # character of the label says which applies.
    sub = keys().subshell('W69+  s24f', 4)
    entry = sub.codes['1']
    assert entry.parts == ()
    assert [j for j, _ in entry.by_j] == [2.5, 3.5]
    assert entry.partition(2.5)[0].sign == '-'
    assert entry.partition(3.5)[0].sign == '+'
    with pytest.raises(jjkeys.JJKeyError):
        entry.partition()               # which J?  say so.
    with pytest.raises(jjkeys.JJKeyError):
        entry.partition(0.5)


def test_an_s_subshell_has_only_one_j():
    # l=0: there is no j-, so the entry has a single line even though
    # the subshell is one of the unrotated ones.
    sub = keys().subshell('W69+  ss3p', 1)
    assert (sub.l, sub.nel) == ('s', 1)
    parts = sub.codes['1'].partition()
    assert [(p.sign, p.j) for p in parts] == [('+', 0.5)]


def test_an_unknown_subshell_is_named():
    with pytest.raises(jjkeys.JJKeyError):
        keys().subshell('W69+  2p3', 7)
