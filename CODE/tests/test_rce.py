"""RCE's output files, against hand-built lines and the reference runs."""

import pytest

from conftest import case
from toolbox import rce


def out(name):
    return rce.read_rceout(folder=case(name))


# --------------------------------------------------------------------------
# one line at a time
# --------------------------------------------------------------------------

RCEOUT_LINE = ('     0.000000      0.000   1/2 -99 2s     (1S) 2S     '
               '   2 2p2    (3P) 4P     ')


def test_an_rceout_level_line_is_read_field_by_field():
    levels, _ = rce._read_rceout_levels([RCEOUT_LINE.ljust(140)], 0, 1, 'x')
    level, = levels
    assert level.e_obs == 0.0
    assert level.e_calc == 0.0
    assert level.j == 0.5
    assert level.observed                       # no flag in column 14
    assert [(c.percent, c.configuration, c.label) for c in level.components] \
        == [(-99, '2s', '(1S) 2S'), (2, '2p2', '(3P) 4P')]


def test_a_whole_j_is_written_without_the_halving():
    # `I4,A2`: RCE halves 2J and blanks the A2 where J comes out whole
    # (`FOR/RCE20KD.F:2055-2062`), so the A2 is what says which is which.
    line = ('%13.6f%1s%10.3f%4d%2s' % (0.0, ' ', 0.0, 2, '  ')
            + '%4d %-6s %-12s' % (100, '2p2', '(3P) 3P'))
    levels, _ = rce._read_rceout_levels([line.ljust(140)], 0, 1, 'x')
    assert levels[0].j == 2.0


def test_a_computed_level_is_flagged_in_column_14():
    line = ('     0.000000*     0.000   1/2 100 s22p  (1S) 2P     ')
    levels, _ = rce._read_rceout_levels([line.ljust(140)], 0, 1, 'x')
    assert levels[0].flag == '*'
    assert not levels[0].observed


def test_an_unknown_energy_layout_is_named_rather_than_guessed():
    # Neither of RCE's two layouts puts the point in column 5.
    line = '  0.12345678      0.000   1/2 100 2s     (1S) 2S     '
    with pytest.raises(rce.RceError) as e:
        rce._read_rceout_levels([line.ljust(140)], 0, 1, 'RCEOUT')
    assert 'column 4' in str(e.value) and '7 or 9' in str(e.value)


def test_a_levels_line_carries_the_g_factor():
    line = ('     0.000000      0.0000 -0.000019  2.002 2s      J= 0.5 '
            '    100% 2s     (1S) 2S')
    level, layout, _ = rce._read_levels_line([line.ljust(140)], 0, None,
                                             'LEVELS1', 1)
    assert level.g == pytest.approx(2.002)
    assert level.j == 0.5
    assert level.configuration == '2s'
    assert [(c.percent, c.configuration, c.label) for c in level.components] \
        == [(100, '2s', '(1S) 2S')]


# --------------------------------------------------------------------------
# whole files
# --------------------------------------------------------------------------

REFERENCE = [('U90', [(3, 10), (2, 9)]),
             ('U90/reordered', [(3, 10), (2, 9)]),
             ('W70', [(12, 95), (8, 74)]),
             ('Pr3', [(21, 1844), (18, 3782)]),
             # Same shape as its parent: the JK deck is the same case
             # re-run for a second coupling, nothing else.
             ('Pr3/JK', [(21, 1844), (18, 3782)])]


def test_every_reference_rceout_reads_whole():
    for name, shape in REFERENCE:
        csets = out(name)
        assert [(len(c.configurations), len(c)) for c in csets] == shape


def test_levels1_repeats_rceout_once_per_coupling_scheme():
    # LEVELS1 carries one section per parity set per coupling, so a
    # two-parity run gives four (`FOR/RCE20KD.F:1301-1305`).
    for name, shape in REFERENCE:
        levels = rce.read_levels(folder=case(name))
        assert len(levels) == 2 * len(shape)
        assert [len(s) for s in levels] == [n for _, n in shape
                                            for _ in (0, 1)]


def test_the_true_jj_pass_writes_its_own_level_list():
    # Only W70 sets KCPLD(8)=8, and its true-jj list goes to LEVELS1J so
    # that LEVELS1 comes out exactly as it always did (`:2017-2019`).
    levels = rce.read_levels('LEVELS1J', folder=case('W70'))
    assert [len(s) for s in levels] == [95, 74]


def test_a_level_agrees_between_rceout_and_levels1():
    # LEVELS1 prints one more digit of the fitted energy, which is the
    # only reason to read both.
    cset = out('U90')[0]
    section = rce.read_levels(folder=case('U90'))[0]
    for a, b in zip(cset.levels, section.levels):
        assert a.j == b.j
        assert a.e_obs == pytest.approx(b.e_obs, abs=5e-4)
        assert a.e_calc == pytest.approx(b.e_calc, abs=5e-3)


def test_the_parameter_table_carries_flags_and_ci_pairs():
    parameters = out('U90')[0].parameters
    by = {p.name: p for p in parameters}
    assert by['EAV 2s'].value == pytest.approx(1.405684)
    assert by['ZETA 3'].value == pytest.approx(23258.095202)
    assert not by['ZETA 3'].is_ci
    assert by['120D1112'].configurations == ('2s', '2s2')
    assert by['120D1112'].is_ci
    assert by['ALPHA'].flag == -100


def test_oute_gives_one_standard_deviation_per_parameter():
    # OUTE is the only file that prints SDX, and its heading names the
    # parity set, so nothing has to be counted (`FOR/RCE20KD.F` listing).
    oute = rce.read_oute(folder=case('U90'))
    assert set(oute) == {1, 2}
    assert len(oute[1]) == len(out('U90')[0].parameters)
    by = {p.name: p for p in oute[2]}
    assert by['EAV 2p'].sigma == pytest.approx(246.087)
    assert by['G1(23)'].value == pytest.approx(2010.7517)


def test_the_last_cycle_is_the_one_that_is_kept():
    # Every cycle writes the same table; what matters is the converged
    # one, so a later block overwrites an earlier one of the same set.
    oute = rce.read_oute(folder=case('U90'))
    assert oute[1][0].name == 'EAV 2s'
    assert oute[1][0].value == pytest.approx(1.4057)


def test_the_hartree_fock_values_read_in_the_same_shape():
    hf = rce.read_hf(folder=case('U90'))
    fitted = out('U90')
    for k, cset in zip((1, 2), fitted):
        assert [p.name for p in hf[k]] == [p.name for p in cset.parameters]
    # The starting value differs from the fitted one; that difference is
    # the LSF/HF column of params.txt.
    assert hf[1][0].value != fitted[0].parameters[0].value


def test_a_run_without_rceinp_hf_is_not_an_error():
    assert rce.read_hf(folder=case('W70') / 'nothing-here') == {}


# --------------------------------------------------------------------------
# OUTGINE, the deck's parameters under RCE's names
# --------------------------------------------------------------------------

def _outgine(names, flags, values, spectrum='X1+', configurations=('2s',)):
    """One OUTGINE section, in the records `read_outgine` counts on."""
    def block(items, form):
        return [''.join(form % x for x in items[k:k + 7])
                for k in range(0, len(items), 7)]
    head = ['    1   30    5   -6    0    0    0',
            '    1    0    1    1     0.000    1    3    1.0000',
            '%-8s  ' % spectrum + ''.join('%-10s' % c
                                          for c in configurations)]
    return (head + block(names, '%-10s')
            + ['    1.0000    2.0000', '         1         2']
            + block(flags, '%10d') + block([0.0] * len(names), '%10.4f')
            + block(values, '%10.4f')
            + ['       1.0       0.0      -1.0 0.00 1000',
               '   -1'])


def test_outgine_finds_the_parameters_by_counting_back(tmp_path):
    # The three parameter records are the last three blocks before the
    # convergence card, whatever the eigenvalues in between look like.
    names = ['EAV 2s', 'F2(11)', 'ZETA 1', '120D1112'] * 3
    values = [float(i) for i in range(12)]
    path = tmp_path / 'OUTGINE'
    path.write_text('\n'.join(_outgine(names, [-100] * 12, values)) + '\n')
    got = rce.read_outgine(folder=tmp_path)
    assert list(got) == [1]
    assert [p.name for p in got[1]] == names
    assert [p.value for p in got[1]] == values
    assert all(p.flag == -100 for p in got[1])


def test_outgine_refuses_a_section_whose_blocks_do_not_line_up(tmp_path):
    # The all-zero record is the check that the count landed right.
    lines = _outgine(['EAV 2s', 'F2(11)'], [-100, -100], [1.0, 2.0])
    lines.insert(-2, '    9.0000    9.0000')
    path = tmp_path / 'OUTGINE'
    path.write_text('\n'.join(lines) + '\n')
    with pytest.raises(rce.RceError):
        rce.read_outgine(folder=tmp_path)


def test_a_folder_with_no_outgine_is_not_an_error(tmp_path):
    # printout --no-rce reports that it cannot write params.txt; it does
    # not stop, because the level list needs none of this.
    assert rce.read_outgine(folder=tmp_path) == {}


def test_outgine_carries_the_ing11_values():
    from toolbox import ing11
    folder = case('Pr3/JK')
    got = rce.read_outgine(folder=folder)
    deck = ing11.read_ing11(folder=folder)
    # The first card of parity set 1: EAV and the parameters after it,
    # in the order RCG punched them.
    card = deck.parameters[0][0]
    values = [p.value for p in got[1][:len(card.values)]]
    assert values == [round(v, 4) for v in card.values]
