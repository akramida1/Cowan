"""insert_level_ids.py: unique level identifiers into conv_out's tables.

The two files this writes into are the ones `conv_out` writes, and the
shapes that matter are the awkward ones: a level list holds one row per
level **per coupling block**, so an identifier goes on several rows; a
line list has a preamble line above its header; and both have their
trailing tabs stripped off, so a row can be shorter than the header.  A
table built here reproduces all three.

The matching itself is `(parity, J, energy)`, and the two tests that
carry the weight are the ones that could go wrong silently: two source
levels of one J within the tolerance of one target level, where the
nearer must win and the other must be reported rather than displace it,
and two levels of one J and one energy in opposite parities, which must
not be swapped.
"""

import os
import subprocess
import sys

import pytest

from conftest import WORK
from toolbox import levids

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SOURCE = WORK / 'Pr3' / 'Pr3_exp_levs.txt'

#: The columns of `Pr3_exp_levs.txt`; its identifier column is called
#: `level_id`, which is the default, so it is not named here.
PR3 = ['--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity']

LEVEL_HEADER = ['Cpl', 'par', 'Lev_#', 'Ee', 'e_c', 'Ec', 'J', 'Lande_g',
                'AR_tot', '%1', 'conf1', 'term1']
LINE_HEADER = ['level1', '', 'J1', 'level2', '', 'J2', 'E1', 'E2',
               'lid1', 'lid2', 'wl(A)', 'Icalc']
PREAMBLE = ['Spectrum', 'Pr III', '', 'LTE', '', '', '', 'Max. intens.',
            '', '1000']


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'insert_level_ids.py')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True)


def write(path, rows, ending='\n'):
    """Write a table the way conv_out does: trailing tabs stripped."""
    with open(str(path), 'w', newline='') as handle:
        for row in rows:
            handle.write('\t'.join(row).rstrip('\t') + ending)
    return str(path)


def read(path):
    with open(str(path), newline='') as handle:
        return [line.rstrip('\r\n').split('\t')
                for line in handle if line.strip('\r\n')]


def source_file(path, levels, header=('level_id', 'energy', 'j_val',
                                      'parity')):
    """A source list: (identifier, energy in cm-1, J, parity) per level."""
    rows = [list(header)]
    for ident, energy, j, parity in levels:
        rows.append([ident, '%.6f' % energy, j, parity])
    return write(path, rows)


def level_file(path, levels, couplings=('LS', 'JJ')):
    """A level list: (Lev_#, energy in kK, J, parity, calculated) per level."""
    rows = [list(LEVEL_HEADER)]
    for coupling in couplings:
        for number, energy, j, parity, calculated in levels:
            rows.append([coupling, parity, '%d' % number, '%12.6f' % energy,
                         '*' if calculated else '', '%12.4f' % (energy + 0.5),
                         '%4.1f' % j, '1.234', '1.00e+08', '95', '4f3', '2P*'])
    return write(path, rows)


def line_file(path, pairs):
    """A line list, preamble and all: (lid1, lid2) per transition."""
    rows = [list(PREAMBLE), list(LINE_HEADER)]
    for low, high in pairs:
        rows.append(['4f3', '2P*', '0.5', '4f3', '4D*', '1.5', '1.0', '2.0',
                     '%d' % low, '%d' % high, '1234.5', '100'])
    return write(path, rows)


def column(rows, name):
    """One column of a table, by header name, blank where a row is short."""
    at = [field.strip() for field in rows[0]].index(name)
    return [row[at].strip() if at < len(row) else '' for row in rows[1:]]


# --------------------------------------------------------------------------
# the pieces
# --------------------------------------------------------------------------

@pytest.mark.parametrize('text, doubled', [
    ('1/2', 1), ('3/2', 3), ('0.5', 1), ('.5', 1), ('4', 8), (' 4.0', 8),
    ('4.0*', 8), ('0', 0), ('"5/2"', 5),
])
def test_j_key_reads_every_way_a_j_is_written(text, doubled):
    assert levids.j_key(text) == doubled


@pytest.mark.parametrize('text', ['', '  ', 'x', '1/0', '0.3', '2P*'])
def test_j_key_refuses_what_is_not_a_j(text):
    assert levids.j_key(text) is None


def test_j_text_is_the_way_a_level_list_writes_it():
    assert [levids.j_text(n) for n in (0, 1, 3, 8)] == ['0', '1/2', '3/2', '4']


def test_units_scale_the_source_energies(tmp_path):
    source = source_file(tmp_path / 'SRC', [('a', 18693.65, '1/2', 'o')])
    levels, _ = levids.read_source(source, levids.Spec(
        ecol='energy', jcol='j_val', pcol='parity'))
    assert levels[0].energy == pytest.approx(18.69365)
    levels, _ = levids.read_source(source, levids.Spec(
        ecol='energy', jcol='j_val', pcol='parity', units='kK'))
    assert levels[0].energy == pytest.approx(18693.65)


def test_a_column_may_be_named_by_number(tmp_path):
    source = source_file(tmp_path / 'SRC', [('a', 1000.0, '1/2', 'o')])
    levels, _ = levids.read_source(source, levids.Spec(
        idcol='1', ecol='2', jcol='3', pcol='4'))
    assert levels[0].ident == 'a'


def test_a_row_with_no_identifier_is_counted_and_left_out(tmp_path):
    source = source_file(tmp_path / 'SRC', [('a', 1000.0, '1/2', 'o'),
                                            ('', 2000.0, '1/2', 'o')])
    levels, report = levids.read_source(
        source, levids.Spec(ecol='energy', jcol='j_val', pcol='parity'))
    assert [one.ident for one in levels] == ['a']
    assert (report.rows, report.levels, report.no_id) == (2, 1, 1)


def test_a_repeated_identifier_is_reported(tmp_path):
    source = source_file(tmp_path / 'SRC', [('a', 1000.0, '1/2', 'o'),
                                            ('a', 2000.0, '3/2', 'o')])
    _, report = levids.read_source(
        source, levids.Spec(ecol='energy', jcol='j_val', pcol='parity'))
    assert report.repeated == {'a': 2}


def test_a_missing_energy_column_says_so(tmp_path):
    source = source_file(tmp_path / 'SRC', [('a', 1000.0, '1/2', 'o')])
    with pytest.raises(levids.LevelIdError) as problem:
        levids.read_source(source, levids.Spec(jcol='j_val', pcol='parity'))
    assert '--ecol' in str(problem.value)


def test_a_table_that_is_not_conv_outs_says_so(tmp_path):
    path = write(tmp_path / 'X', [['a', 'b'], ['1', '2']])
    with pytest.raises(levids.LevelIdError) as problem:
        levids.Table.read(path, levids.LEVEL_KEYS)
    assert 'Lev_#' in str(problem.value)


# --------------------------------------------------------------------------
# the level list
# --------------------------------------------------------------------------

def one_case(tmp_path, source_levels, target_levels, pairs=None):
    source = source_file(tmp_path / 'SRC', source_levels)
    levels = level_file(tmp_path / 'levels.txt', target_levels)
    lines = line_file(tmp_path / 'lines.txt', pairs) if pairs else None
    arguments = [source, 'levels.txt'] + (['lines.txt'] if lines else []) + PR3
    return run(tmp_path, *arguments), levels, lines


def test_the_identifier_goes_in_beside_ee_on_every_block(tmp_path):
    done, levels, _ = one_case(
        tmp_path,
        [('id-1', 18693.65, '1/2', 'o'), ('id-2', 23465.43, '3/2', 'e')],
        [(1, 18.69365, 0.5, 'o', False), (2, 23.46543, 1.5, 'e', False)])
    assert done.returncode == 0, done.stderr
    rows = read(levels)
    assert rows[0][2:5] == ['Lev_#', 'level_id', 'Ee']
    # Two coupling blocks, so each identifier is on two rows.
    assert column(rows, 'level_id') == ['id-1', 'id-2'] * 2


def test_a_level_no_source_level_matches_is_left_blank(tmp_path):
    done, levels, _ = one_case(
        tmp_path, [('id-1', 18693.65, '1/2', 'o')],
        [(1, 18.69365, 0.5, 'o', False), (2, 99.0, 0.5, 'o', False)])
    assert column(read(levels), 'level_id') == ['id-1', '', 'id-1', '']
    assert '1 target level left blank' in done.stdout


def test_a_source_level_no_level_matches_is_reported(tmp_path):
    done, _, _ = one_case(
        tmp_path, [('id-1', 18693.65, '1/2', 'o'), ('id-2', 50000.0, '1/2',
                                                    'o')],
        [(1, 18.69365, 0.5, 'o', False)])
    assert '1 source level unmatched' in done.stdout
    assert 'unmatched: id-2' in done.stdout


def test_parity_tells_two_levels_of_one_j_and_energy_apart(tmp_path):
    done, levels, _ = one_case(
        tmp_path,
        [('odd', 18693.65, '1/2', 'o'), ('even', 18693.65, '1/2', 'e')],
        [(1, 18.69365, 0.5, 'e', False), (2, 18.69365, 0.5, 'o', False)])
    assert done.returncode == 0, done.stderr
    assert column(read(levels), 'level_id')[:2] == ['even', 'odd']


def test_j_tells_two_levels_of_one_energy_apart(tmp_path):
    _, levels, _ = one_case(
        tmp_path, [('half', 18693.65, '1/2', 'o'),
                   ('three-halves', 18693.65, '3/2', 'o')],
        [(1, 18.69365, 1.5, 'o', False), (2, 18.69365, 0.5, 'o', False)])
    assert column(read(levels), 'level_id')[:2] == ['three-halves', 'half']


def test_the_nearer_of_two_source_levels_takes_the_target(tmp_path):
    """Both are within the tolerance; the closer one must win.

    Level by level, the first row of the file would take the target and
    the second would find it gone, whichever was nearer.  The offers
    are sorted instead, so the order of the source file cannot decide.
    """
    done, levels, _ = one_case(
        tmp_path, [('far', 18694.30, '1/2', 'o'), ('near', 18693.70, '1/2',
                                                   'o')],
        [(1, 18.69365, 0.5, 'o', False)])
    assert column(read(levels), 'level_id')[0] == 'near'
    assert 'unmatched: far' in done.stdout


def test_the_tolerance_is_in_the_sources_own_unit(tmp_path):
    """0.65 cm-1 apart: inside the default, outside --tol 0.5."""
    source = [('id-1', 18694.30, '1/2', 'o')]
    target = [(1, 18.69365, 0.5, 'o', False)]
    _, levels, _ = one_case(tmp_path, source, target)
    assert column(read(levels), 'level_id')[0] == 'id-1'

    tight = tmp_path / 'tight'
    tight.mkdir()
    written = source_file(tight / 'SRC', source)
    level_file(tight / 'levels.txt', target)
    done = run(tight, written, 'levels.txt', '--tol', '0.5', *PR3)
    assert column(read(tight / 'levels.txt'), 'level_id')[0] == ''
    assert 'within 0.5 cm-1' in done.stdout


def test_a_match_to_a_calculated_energy_is_flagged(tmp_path):
    done, levels, _ = one_case(
        tmp_path, [('id-1', 18693.65, '1/2', 'o')],
        [(1, 18.69365, 0.5, 'o', True)])
    assert column(read(levels), 'level_id')[0] == 'id-1'
    assert 'calculated energy' in done.stdout


def test_an_observed_energy_wins_a_tie_with_a_calculated_one(tmp_path):
    _, levels, _ = one_case(
        tmp_path, [('id-1', 18693.65, '1/2', 'o')],
        [(1, 18.69365, 0.5, 'o', True), (2, 18.69365, 0.5, 'o', False)])
    assert column(read(levels), 'level_id')[:2] == ['', 'id-1']


# --------------------------------------------------------------------------
# the line list
# --------------------------------------------------------------------------

def test_the_two_identifiers_go_in_beside_lid2(tmp_path):
    done, _, lines = one_case(
        tmp_path,
        [('id-1', 18693.65, '1/2', 'o'), ('id-2', 23465.43, '3/2', 'e')],
        [(1, 18.69365, 0.5, 'o', False), (2, 23.46543, 1.5, 'e', False)],
        pairs=[(1, 2), (2, 1)])
    assert done.returncode == 0, done.stderr
    rows = read(lines)
    assert rows[0] == PREAMBLE                     # kept as it was
    assert rows[1][8:12] == ['lid1', 'lid2', 'id1', 'id2']
    assert column(rows[1:], 'id1') == ['id-1', 'id-2']
    assert column(rows[1:], 'id2') == ['id-2', 'id-1']


def test_a_line_naming_a_level_the_list_has_not_got_is_reported(tmp_path):
    done, _, lines = one_case(
        tmp_path, [('id-1', 18693.65, '1/2', 'o')],
        [(1, 18.69365, 0.5, 'o', False)], pairs=[(1, 777)])
    assert 'level 777' in done.stdout
    assert column(read(lines)[1:], 'id2') == ['']


def test_a_level_with_no_identifier_leaves_the_lines_blank(tmp_path):
    _, _, lines = one_case(
        tmp_path, [('id-1', 18693.65, '1/2', 'o')],
        [(1, 18.69365, 0.5, 'o', False), (2, 99.0, 0.5, 'o', False)],
        pairs=[(1, 2)])
    rows = read(lines)[1:]
    assert (column(rows, 'id1'), column(rows, 'id2')) == (['id-1'], [''])


def test_a_blank_header_field_after_lid2_is_not_a_name_to_insert_before(
        tmp_path):
    """The term columns of a line list have no header of their own."""
    source = source_file(tmp_path / 'SRC', [('id-1', 18693.65, '1/2', 'o'),
                                            ('id-2', 23465.43, '3/2', 'e')])
    level_file(tmp_path / 'levels.txt', [(1, 18.69365, 0.5, 'o', False),
                                         (2, 23.46543, 1.5, 'e', False)])
    write(tmp_path / 'lines.txt',
          [list(PREAMBLE),
           ['lid1', 'lid2', '', 'wl(A)'],
           ['1', '2', 'x', '1234.5']])
    done = run(tmp_path, source, 'levels.txt', 'lines.txt', *PR3)
    assert done.returncode == 0, done.stderr
    rows = read(tmp_path / 'lines.txt')
    assert rows[1] == ['lid1', 'lid2', 'id1', 'id2', '', 'wl(A)']
    assert rows[2] == ['1', '2', 'id-1', 'id-2', 'x', '1234.5']


def test_a_line_list_without_the_lid_columns_says_so(tmp_path):
    source = source_file(tmp_path / 'SRC', [('id-1', 18693.65, '1/2', 'o')])
    level_file(tmp_path / 'levels.txt', [(1, 18.69365, 0.5, 'o', False)])
    write(tmp_path / 'lines.txt', [['a', 'b'], ['1', '2']])
    done = run(tmp_path, source, 'levels.txt', 'lines.txt', *PR3)
    assert done.returncode == 1
    assert 'lid1, lid2' in done.stderr


def test_the_line_list_may_be_left_out(tmp_path):
    done, levels, _ = one_case(tmp_path, [('id-1', 18693.65, '1/2', 'o')],
                               [(1, 18.69365, 0.5, 'o', False)])
    assert done.returncode == 0, done.stderr
    assert column(read(levels), 'level_id')[0] == 'id-1'
    assert not (tmp_path / 'lines.txt').exists()


# --------------------------------------------------------------------------
# running it twice
# --------------------------------------------------------------------------

def test_a_second_run_updates_the_columns_rather_than_repeating_them(tmp_path):
    arguments = ([source_file(tmp_path / 'SRC',
                              [('id-1', 18693.65, '1/2', 'o'),
                               ('id-2', 23465.43, '3/2', 'e')]),
                  'levels.txt', 'lines.txt'] + PR3)
    level_file(tmp_path / 'levels.txt', [(1, 18.69365, 0.5, 'o', False),
                                         (2, 23.46543, 1.5, 'e', False)])
    line_file(tmp_path / 'lines.txt', [(1, 2)])

    assert run(tmp_path, *arguments).returncode == 0
    first = ((tmp_path / 'levels.txt').read_bytes(),
             (tmp_path / 'lines.txt').read_bytes())
    assert run(tmp_path, *arguments).returncode == 0
    assert ((tmp_path / 'levels.txt').read_bytes(),
            (tmp_path / 'lines.txt').read_bytes()) == first
    assert read(tmp_path / 'levels.txt')[0].count('level_id') == 1
    assert read(tmp_path / 'lines.txt')[1].count('id1') == 1
    # The backup is what the second run was given.
    assert (tmp_path / 'levels.txt.bak').read_bytes() == first[0]


def test_an_identifier_that_no_longer_matches_is_cleared(tmp_path):
    """The column is emptied before it is filled, so nothing goes stale."""
    source = source_file(tmp_path / 'SRC', [('id-1', 18693.65, '1/2', 'o')])
    level_file(tmp_path / 'levels.txt', [(1, 18.69365, 0.5, 'o', False),
                                         (2, 99.0, 0.5, 'o', False)])
    assert run(tmp_path, source, 'levels.txt', *PR3).returncode == 0

    # The second level is given an identifier by hand, and the source
    # still says nothing about it.
    rows = read(tmp_path / 'levels.txt')
    at = rows[0].index('level_id')
    rows[2][at] = 'stale'
    write(tmp_path / 'levels.txt', rows)
    assert run(tmp_path, source, 'levels.txt', *PR3).returncode == 0
    assert column(read(tmp_path / 'levels.txt'), 'level_id') \
        == ['id-1', '', 'id-1', '']


def test_the_line_endings_of_the_two_files_are_kept(tmp_path):
    source = source_file(tmp_path / 'SRC', [('id-1', 18693.65, '1/2', 'o')])
    rows = [list(LEVEL_HEADER),
            ['LS', 'o', '1', '%12.6f' % 18.69365, '', '  19.1937', ' 0.5',
             '1.234', '1.00e+08', '95', '4f3', '2P*']]
    write(tmp_path / 'levels.txt', rows, ending='\r\n')
    assert run(tmp_path, source, 'levels.txt', *PR3).returncode == 0
    written = (tmp_path / 'levels.txt').read_bytes()
    assert b'\r\n' in written and b'\n\n' not in written


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

def test_no_arguments_prints_the_usage():
    done = run(CODE)
    assert done.returncode == 0
    assert 'insert_level_ids <source>' in done.stdout


@pytest.mark.parametrize('arguments, said', [
    (['SRC'], 'file name was given'),
    (['SRC', 'L', 'N', 'X'], 'file names were given'),
    (['SRC', 'L', '--nonesuch', '1'], '--nonesuch is not an option'),
    (['SRC', 'L', '--units', 'furlongs'], 'is not a unit'),
    (['SRC', 'L', '--tol', 'x'], '--tol wants a number'),
    (['SRC', 'L', '--skip'], '--skip needs a value'),
])
def test_a_command_line_that_will_not_do(tmp_path, arguments, said):
    done = run(tmp_path, *arguments)
    assert done.returncode == 2
    assert said in done.stderr


# --------------------------------------------------------------------------
# the author's own source file
# --------------------------------------------------------------------------

def test_the_pr3_sheet_reads_with_the_default_identifier_column():
    if not SOURCE.is_file():
        pytest.skip('WORK/Pr3 has no Pr3_exp_levs.txt')
    levels, report = levids.read_source(str(SOURCE), levids.Spec(
        ecol='energy', jcol='j_val', pcol='parity'))
    assert report.levels == len(levels) > 500
    assert not report.repeated                 # the identifiers are unique
    assert not report.bad_j and not report.bad_parity
    assert all(one.ident.startswith('059003.') for one in levels)
