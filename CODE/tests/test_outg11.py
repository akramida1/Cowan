"""OUTG11: the configuration table, the J blocks, the eigenvectors."""

import pytest

from conftest import case
from toolbox import genealogy, outg11


def run_of(name, folder=None):
    return outg11.read_outg11((case(name) / (folder or '')) / 'OUTG11')


def blocks_of(name, folder=None, couplings=None):
    run = run_of(name, folder)
    return run, {(b.kset, b.j): b for b in run.blocks(couplings)}


# --------------------------------------------------------------------------
# the configuration table
# --------------------------------------------------------------------------

def test_w70_configuration_table():
    run = run_of('W70')
    assert len(run.configs) == 20
    assert len(run.cards(1)) == 12 and len(run.cards(2)) == 8
    card = run.card(1, 1)
    assert card.label == 's22p'
    assert card.spectrum == 'W69+'
    # Three, not the five of IN36: ING11 carries only the subshells RCG
    # works with, and RCN2 has already taken the 1s2 that every
    # configuration of this run shares out of the deck.  The electron
    # count of a label therefore cannot come from here.
    assert card.nelec == 3
    assert card.parity == -1
    assert card.terms == 'ALL'
    assert card.shells[:4] == (('S', 2), ('P', 1), ('S', 0), ('P', 0))
    assert card.open_shells == (('P', 1),)
    assert card.eav == pytest.approx(-1751028.878)


def test_a_truncated_configuration_has_no_name_in_the_table():
    # Where the terms of a subshell are cut back, FORMAT 21
    # (`rcg11kd.f:3420`) puts the subshell number and the surviving term
    # list where FORMAT 23 puts the 3A6 name and the average energy, so
    # neither is printed at all.
    #
    # No reference deck is truncated today -- Pr3/JK was, through a
    # mispunched main card, and the author has since corrected it -- so
    # the row is built here from the FORMAT rather than read from a run.
    # The rule still has to hold: NLT11 and the per-configuration cutoff
    # in columns 83-85 are both live features, and a tool that took this
    # table for the list of names would print blanks for them.
    line = ('%4d%4d%7s' % (1, 1, '')
            + ''.join('%1s%2d%3s' % (l, n, '')
                      for l, n in [('F', 3), ('P', 6)] + [(' ', 0)] * 6)
            + '%5d%5d%5d%3s' % (9, 1, 1, '')
            + '4I  2H  2G  2F  2D  4S  4F  4G')
    card = outg11._config_row(line)
    assert card.shells[0] == ('F', 3)
    assert card.name == '' and card.eav == 0.0
    assert card.terms.startswith('1: 4I')

    # An ordinary row of the same run keeps both.
    plain = run_of('Pr3').card(1, 1)
    assert plain.shells[0] == ('F', 3)
    assert plain.label == '4f3' and plain.eav != 0.0


# --------------------------------------------------------------------------
# the J blocks
# --------------------------------------------------------------------------

def test_w70_blocks_and_their_sets():
    run, blocks = blocks_of('W70')
    assert sorted(blocks) == [(1, 0.5), (1, 1.5), (1, 2.5), (1, 3.5), (1, 4.5),
                              (2, 0.5), (2, 1.5), (2, 2.5), (2, 3.5), (2, 4.5),
                              (2, 5.5)]
    b = blocks[(1, 0.5)]
    assert b.nlevels == 24
    assert len(b.config_numbers) == 24
    assert len(b.g_values) == 24
    assert b.eigenvalues[0] == 0.0
    assert b.eigenvalues[-1] == pytest.approx(206071.19)
    # Eigenvalues are F9.2 and run together above 100000 cm-1, so the
    # column grid, not the blanks, is what separates them.
    assert b.eigenvalues[6] == pytest.approx(100621.53)
    assert b.eigenvalues[7] == pytest.approx(104677.46)


def test_the_three_sections_of_a_w70_block():
    _, blocks = blocks_of('W70')
    b = blocks[(1, 0.5)]
    assert sorted(b.sections) == ['JJ', 'LS', 'TRUEjj']
    for name in ('LS', 'JJ', 'TRUEjj'):
        s = b.section(name)
        assert s.nstates == 24
        assert all(len(r.amplitudes) == 24 for r in s.rows)
        assert [r.number for r in s.rows] == list(range(1, 25))
    ls = b.section('LS')
    assert ls.rows[0].label == 's22p'
    assert ls.rows[0].bracketed == '1S) 2P'
    assert ls.rows[0].amplitudes[0] == pytest.approx(0.99976)
    assert ls.rows[0].amplitudes[1] == pytest.approx(-0.02135)
    assert b.section('JJ').rows[0].bracketed == '1S 0) 1'
    assert b.section('TRUEjj').rows[4].bracketed == '1212314'


def test_only_the_primary_coupling_carries_column_headings():
    # RCG writes the heading once, for the coupling the run was set up
    # in; the other sections are rows only, in the same column order.
    _, blocks = blocks_of('W70')
    b = blocks[(1, 0.5)]
    assert b.section('LS').column_labels[:3] == ['s22p', '2p3', 's23p']
    assert b.section('LS').column_names[:3] == ['1S) 2P', '2P) 2P', '1S) 2P']
    assert len(b.section('LS').column_labels) == 24
    assert b.section('JJ').column_labels == []
    assert b.section('TRUEjj').column_labels == []


def test_purity_is_read_past_its_first_line():
    _, blocks = blocks_of('W70')
    ls = blocks[(1, 0.5)].section('LS')
    assert ls.average_purity == pytest.approx(0.809)
    assert len(ls.purity) == 24
    assert ls.purity[0] == pytest.approx(0.99952)
    assert ls.purity[-1] == pytest.approx(0.89486)


def test_a_one_by_one_block_has_no_jj_section():
    # This is the whole of the "blank configuration" problem: RCG prints
    # no JJ and no TRUEjj section where the block is one state, so the
    # level must be built from genealogy.txt instead.  Every block that
    # lacks a section in the reference runs is of this kind.
    _, blocks = blocks_of('W70')
    b = blocks[(2, 5.5)]
    assert b.nlevels == 1
    assert sorted(b.sections) == ['LS']
    with pytest.raises(outg11.Outg11Error) as bad:
        b.section('JJ')
    assert 'column 32' in str(bad.value)     # KCPLD(2), ING11 col 32


def test_every_block_without_a_jj_section_is_one_by_one():
    for name, folder in (('W70', None), ('U90', None), ('U90', 'reordered')):
        _, blocks = blocks_of(name, folder)
        for key, b in blocks.items():
            if 'JJ' not in b.sections:
                assert b.nlevels == 1, (name, folder, key)


def test_couplings_may_be_limited():
    _, blocks = blocks_of('W70', couplings={'LS'})
    b = blocks[(1, 0.5)]
    assert sorted(b.sections) == ['LS']
    assert b.section('LS').nstates == 24


# --------------------------------------------------------------------------
# against the other readers
# --------------------------------------------------------------------------

def test_labels_agree_with_genealogy():
    _, blocks = blocks_of('W70')
    gen = genealogy.read_genealogy(case('W70') / 'genealogy.txt')
    pairs = ((genealogy.LS, 'LS'), (genealogy.J1J2, 'JJ'),
             (genealogy.TRUEJJ, 'TRUEjj'))
    checked = 0
    for gb in gen:
        b = blocks[(gb.kset, gb.j)]
        for table, coupling in pairs:
            if coupling not in b.sections:
                continue
            rows = b.section(coupling).rows
            states = gb.tables[table]
            assert len(rows) == len(states)
            for row, state in zip(rows, states):
                assert row.bracketed == state.label.rstrip()
                assert row.cfg == state.cfg
                checked += 1
    assert checked > 300


def test_pr3_reads_with_its_second_coupling():
    _, blocks = blocks_of('Pr3', 'JK', couplings={'LS', 'LSJK'})
    b = blocks[(1, 0.5)]
    assert sorted(b.sections) == ['LS', 'LSJK']
    assert b.section('LSJK').nstates == b.nlevels
    # Couplings 3 to 7 write '=' where LS writes ':', and put running
    # quantum numbers where LS puts a term label.
    row = b.section('LSJK').rows[0]
    assert row.label == '4f3'
    assert row.desc.rstrip().endswith('K 0.5')


def test_an_absent_coupling_names_its_ing11_column():
    _, blocks = blocks_of('U90')
    with pytest.raises(outg11.Outg11Error) as bad:
        blocks[(1, 0.5)].section('TRUEjj')
    assert 'TRUEjj' in str(bad.value)
    assert 'column 38' in str(bad.value)     # KCPLD(8), ING11 col 38
