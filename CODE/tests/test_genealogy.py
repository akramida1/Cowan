"""genealogy.txt, against hand-built lines and the reference runs."""

import pytest

from conftest import case
from toolbox import genealogy
from toolbox.genealogy import LS, J1J2, TRUEJJ


# --------------------------------------------------------------------------
# one field at a time
# --------------------------------------------------------------------------

def field(text, table):
    """One field, as if it stood in the first slot of a state line."""
    line = ' ' * 22 + '  ' + text
    return genealogy._fields(line, table, 0.5, 1)[0]


def test_ls_field():
    f = field('(   6  2S )  2S', LS)
    assert (f.term, f.mult, f.l, f.alpha) == (6, 2, 'S', '')
    assert (f.run_mult, f.run_l) == (2, 'S')


def test_ls_field_keeps_cowans_alpha():
    # Cowan's alpha is the character after the L of the subshell term.
    f = field('(  12  3H2)  3H', LS)
    assert (f.term, f.mult, f.l, f.alpha) == (12, 3, 'H', '2')
    assert (f.run_mult, f.run_l) == (3, 'H')


def test_j1j2_field():
    f = field('(   6  2P  1.5)  0.5', J1J2)
    assert (f.term, f.mult, f.l) == (6, 2, 'P')
    assert (f.j, f.run_j) == (1.5, 0.5)


def test_true_jj_field_one_part():
    f = field('(p+ 1  3             )   3   1', TRUEJJ)
    assert len(f.parts) == 1
    p = f.parts[0]
    assert (p.l, p.sign, p.nel, p.seniority, p.j2) == ('p', '+', 1, 0, 3)
    assert (f.j2, f.k2) == (3, 1)


def test_true_jj_field_two_parts():
    f = field('(p- 1  1,p+ 2  0     )   1   1', TRUEJJ)
    assert [(p.sign, p.nel, p.j2) for p in f.parts] == [('-', 1, 1),
                                                        ('+', 2, 0)]


def test_true_jj_field_with_a_seniority():
    # The seniority is printed only where J alone does not tell two
    # states of the same j**N apart, and then as 'v' and one digit.
    f = field('(f- 4v2  4,f+ 2  2   )   6   6', TRUEJJ)
    assert [(p.nel, p.seniority, p.j2) for p in f.parts] == [(4, 2, 4),
                                                             (2, 0, 2)]


def test_a_blank_chunk_ends_the_line():
    # GENPRT builds each field into a string first, so a line of fewer
    # fields simply stops; it does not end on an unclosed bracket.
    line = ' ' * 22 + '  ' + '(   1  2P )  2P' + ' ' * 40
    assert len(genealogy._fields(line, LS, 0.5, 1)) == 1


# --------------------------------------------------------------------------
# whole files
# --------------------------------------------------------------------------

def test_w70_has_all_three_tables():
    gen = genealogy.read_genealogy(case('W70') / 'genealogy.txt')
    assert gen.tables() == [LS, J1J2, TRUEJJ]
    assert len(gen) == 11
    block = gen.block(1, 0.5)
    assert block.nstates == 24
    assert all(len(block.tables[t]) == 24 for t in (LS, J1J2, TRUEJJ))
    # The configuration list of the block is printed once, with the
    # mode=1 call, and names a subshell by its slot.
    assert block.configurations[1] == (('S', 2), ('P', 1), ('S', 0), ('P', 0))


def test_w70_state_five_reads_as_the_file_has_it():
    gen = genealogy.read_genealogy(case('W70') / 'genealogy.txt')
    block = gen.block(1, 0.5)

    ls = block.state(LS, 5)
    # The label is A8 and keeps its trailing blanks; the eigenvector
    # sections of OUTG11 print the same eight columns.
    assert (ls.cfg, ls.label) == (6, '3S) 4P  ')
    assert [(f.mult, f.l) for f in ls.fields] == [(2, 'S'), (2, 'S'), (2, 'P')]
    assert [(f.run_mult, f.run_l) for f in ls.fields] == [(2, 'S'), (3, 'S'),
                                                          (4, 'P')]

    jj = block.state(J1J2, 5)
    assert jj.label == '2S 2) 3 '
    assert [f.j for f in jj.fields] == [0.5, 0.5, 1.5]
    assert [f.run_j for f in jj.fields] == [0.5, 1.0, 0.5]

    tj = block.state(TRUEJJ, 5)
    assert tj.label == '1212314 '
    assert [(f.j2, f.k2) for f in tj.fields] == [(1, 1), (1, 2), (3, 1)]
    # Every J of the true jj table is doubled, and its subshell J and
    # running J repeat the J1J2 line above it.
    assert [f.j2 for f in tj.fields] == [round(2 * f.j) for f in jj.fields]
    assert [f.k2 for f in tj.fields] == [round(2 * f.run_j) for f in jj.fields]


def test_the_one_by_one_block_still_has_its_chain():
    # W70's J=11/2 block is one state, so RCG prints no JJ or TRUEjj
    # eigenvector section for it; genealogy.txt carries the chain all
    # the same, which is what lets the parser fill the gap.
    gen = genealogy.read_genealogy(case('W70') / 'genealogy.txt')
    block = gen.block(2, 5.5)
    assert block.nstates == 1
    jj = block.state(J1J2, 1)
    assert [(f.j, f.run_j) for f in jj.fields] == [(0.5, 0.5), (1.5, 2.0),
                                                   (3.5, 5.5)]


def test_pr3_and_u90_read_whole():
    # No true jj in either: neither deck sets KCPLD(8)=8.
    for name, blocks in (('Pr3', 25), ('U90', 6)):
        gen = genealogy.read_genealogy(case(name) / 'genealogy.txt')
        assert len(gen) == blocks
        assert gen.tables() == [LS, J1J2]
        for b in gen:
            for table in (LS, J1J2):
                rows = b.tables[table]
                assert len(rows) == b.nstates
                assert [r.number for r in rows] == list(range(1, b.nstates + 1))


def test_states_carry_every_open_subshell():
    # Pr III runs to four open subshells, past the four fields the J1J2
    # table puts on a line, so the continuation lines must be joined on.
    gen = genealogy.read_genealogy(case('Pr3') / 'genealogy.txt')
    block = gen.block(1, 2.5)
    for table in (LS, J1J2):
        for state in block.tables[table]:
            shells = block.configurations[state.cfg]
            opens = sum(1 for letter, occ in shells
                        if occ and occ != 2 * (2 * 'SPDFGHIKLMNO'.index(letter)
                                               + 1))
            assert len(state.fields) == opens, (table, state.number)


def test_a_missing_block_is_named():
    gen = genealogy.read_genealogy(case('U90') / 'genealogy.txt')
    with pytest.raises(genealogy.GenealogyError):
        gen.block(1, 7.5)
