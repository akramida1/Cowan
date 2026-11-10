"""write_IDEN_level_ids.py and update_IDEN.py: the round trip through IDEN2.

Neither `ENLEV.DAT` nor `TRANS.DAT` is committed -- they are generated,
and the Pr III pair is 1.6 MB -- so every fixture here is built by hand,
in the columns `WORK/Pr3/IDEN2_file_formats.md` documents and
`toolbox/iden.py` writes.

Three things could go wrong quietly and are what these weigh:

  * a record that nothing was changed in must come back out as the bytes
    it went in as, trailing blanks and line ending included, or the
    author's diff of a run is unreadable;
  * the parity marks `_` and `~` must be read as even and odd and not the
    other way round, which no count of matched levels would reveal on its
    own if every level of the file were the same parity;
  * a transition between two fixed levels must end up assigned exactly
    when the optimization used it -- an assignment left standing on a
    line LOPT threw out is the failure that would look like agreement.
"""

import math
import os
import subprocess
import sys

import pytest

from toolbox import iden, levids

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(tool, folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, tool)] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


# --------------------------------------------------------------------------
# the fixtures
# --------------------------------------------------------------------------

def enlev_card(number, e_calc, uncertainty, e_fixed, fixed, j, label):
    """One record in the version 2 layout, as `_write_enlev` writes it."""
    return ('%4d%12.3f%10.3f%12.3f%2s%12.3f%5.1f /%-11s/'
            % (number, e_calc, uncertainty, e_fixed, '*' if fixed else ' ',
               e_fixed - e_calc, j, label))


def enlev_file(path, levels, ending='\r\n'):
    """`levels`: (number, Ec, unc, Ef, fixed, J, label) per record."""
    with open(str(path), 'w', newline='') as handle:
        for one in levels:
            handle.write(enlev_card(*one) + ending)
    return str(path)


def level_card(number, j, label, energy, fixed):
    return ('$%4d    J=%4.1f  %-11s%13.3f%2s '
            % (number, j, label, abs(energy), '*' if fixed else ' '))


def line_card(lower, intensity, e_lower, fixed, energy,
              obs_intensity=0, wavenumber=0.0, difference=0.0, line=0):
    return ('+%4d%5d%12.3f%2s%13.3f%6d%12.3f%11.3f%6d'
            % (lower, intensity, e_lower, '*' if fixed else ' ', energy,
               obs_intensity, wavenumber, difference, line))


def trans_file(path, records, ending='\r\n'):
    with open(str(path), 'w', newline='') as handle:
        for card in records:
            handle.write(card + ending)
    return str(path)


def dlv_file(path, lines, ending='\r\n'):
    """`lines`: (intensity, wavenumber) per row, descending in wavenumber."""
    with open(str(path), 'w', newline='') as handle:
        for intensity, wavenumber in lines:
            handle.write('%5d%14.3f%34s%7.4f' % (intensity, wavenumber, ' ',
                                                 0.004) + ending)
    return str(path)


def table(path, header, rows, ending='\n'):
    with open(str(path), 'w', newline='') as handle:
        handle.write('\t'.join(header) + ending)
        for row in rows:
            handle.write('\t'.join(str(field) for field in row) + ending)
    return str(path)


def source_file(path, levels):
    """The experimental list: (identifier, energy, J, parity) per level."""
    return table(path, ('level_id', 'energy', 'j_val', 'parity'), levels)


#: A four-level spectrum: two odd levels low down and two even ones on
#: top, all four observed, with a fifth level that is calculated only.
#: IDEN2 numbers from the highest energy down, so level 1 is the top.
LEVELS = [
    (1, 30000.0, 5000.0, 30000.0, False, 2.0, 'd2   _3F2F '),
    (2, 25010.0, 10.0, 25000.0, True, 3.0, 'd2   _3F3F '),
    (3, 20020.0, 10.0, 20000.0, True, 2.0, 'd2   _3F3F '),
    (4, 10030.0, 10.0, 10000.0, True, 4.0, 'f2   ~3H3H '),
    (5, 40.0, 10.0, 0.0, True, 3.0, 'f2   ~3H3H '),
]

SOURCE = [
    ('X.000002', 25000.0, '3', 'e'),
    ('X.000003', 20000.0, '2', 'e'),
    ('X.000004', 10000.0, '4', 'o'),
    ('X.000005', 0.0, '3', 'o'),
]


@pytest.fixture
def folder(tmp_path):
    enlev_file(tmp_path / iden.ENLEV, LEVELS)
    source_file(tmp_path / 'levs.txt', SOURCE)
    return tmp_path


# --------------------------------------------------------------------------
# reading the files back
# --------------------------------------------------------------------------

def test_read_enlev_reads_the_documented_columns(folder):
    levels = iden.read_enlev(str(folder / iden.ENLEV))
    assert [one.number for one in levels] == [1, 2, 3, 4, 5]
    assert [one.fixed for one in levels] == [False, True, True, True, True]
    assert [one.j2 for one in levels] == [4, 6, 4, 8, 6]
    assert [one.parity for one in levels] == ['e', 'e', 'e', 'o', 'o']
    assert levels[1].e_calc == 25010.0
    assert levels[1].e_fixed == 25000.0
    assert levels[1].uncertainty == 10.0


def test_a_version_1_enlev_is_refused(tmp_path):
    """The legacy layout is one column narrower and is not read."""
    path = str(tmp_path / iden.ENLEV)
    with open(path, 'w', newline='') as handle:
        handle.write(' 36 1048359.985  1000.000 1048359.985'
                     '          .000  3.0 / 4 (2D) 3D /\r\n')
    with pytest.raises(iden.IdenError) as problem:
        iden.read_enlev(path)
    assert 'version 1' in str(problem.value)


@pytest.mark.parametrize('ending', ['\r\n', '\n'])
def test_a_file_nothing_changed_in_comes_back_unchanged(tmp_path, ending):
    """Trailing blanks and the line ending both survive a rewrite."""
    path = str(tmp_path / iden.ENLEV)
    enlev_file(path, LEVELS, ending)
    before = open(path, 'rb').read()
    iden.write_enlev(iden.read_enlev(path), path)
    assert open(path, 'rb').read() == before

    path = str(tmp_path / iden.TRANS)
    trans_file(path, [level_card(2, 3.0, 'd2   _3F3F ', 25000.0, True),
                      line_card(4, 12, 10000.0, True, 15000.0)], ending)
    before = open(path, 'rb').read()
    iden.write_trans(iden.read_trans(path), path)
    assert open(path, 'rb').read() == before


def test_only_the_field_written_to_moves(tmp_path):
    """A record keeps every column but the ones that were set."""
    path = str(tmp_path / iden.ENLEV)
    enlev_file(path, LEVELS)
    levels = iden.read_enlev(path)
    was = levels[1].card
    levels[1].update(e_fixed=25001.5)
    now = levels[1].card
    assert now[:16] == was[:16]                 # number and Ec
    assert now[16:26] == was[16:26]             # the uncertainty
    assert now[26:38].strip() == '25001.500'
    assert now[40:52].strip() == '-8.500'       # Ef - Ec, kept in step
    assert now[52:] == was[52:]                 # J and the label


# --------------------------------------------------------------------------
# write_IDEN_level_ids
# --------------------------------------------------------------------------

def test_the_table_ties_every_fixed_level(folder):
    done = run('write_IDEN_level_ids.py', folder, 'levs.txt',
               '--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity')
    assert done.returncode == 0, done.stderr
    assert levids.read_iden_ids(str(folder / levids.IDEN_IDS)) == {
        'X.000002': 2, 'X.000003': 3, 'X.000004': 4, 'X.000005': 5}


def test_the_table_is_two_columns_in_level_order(folder):
    run('write_IDEN_level_ids.py', folder, 'levs.txt',
        '--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity')
    rows = [line.rstrip('\r\n').split('\t')
            for line in open(str(folder / levids.IDEN_IDS))
            if line.strip()]
    assert rows[0] == [levids.SOURCE_ID, levids.IDEN_ID]
    assert [row[1] for row in rows[1:]] == ['2', '3', '4', '5']


def test_the_parity_marks_are_not_read_backwards(folder):
    """`~` is odd and `_` is even, and a level list proves which.

    Levels 3 and 4 sit 10000 cm-1 apart and could not be confused by
    energy; levels of one J in opposite parities are what would be, and
    the tie is made on the mark alone.
    """
    enlev_file(folder / iden.ENLEV, [
        (1, 20020.0, 10.0, 20000.0, True, 3.0, 'd2   _3F3F '),
        (2, 20025.0, 10.0, 20000.0, True, 3.0, 'f2   ~3H3H '),
    ])
    source_file(folder / 'levs.txt', [('EVEN', 20000.0, '3', 'e'),
                                      ('ODD', 20000.0, '3', 'o')])
    done = run('write_IDEN_level_ids.py', folder, 'levs.txt',
               '--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity')
    assert done.returncode == 0, done.stderr
    assert levids.read_iden_ids(str(folder / levids.IDEN_IDS)) == {
        'EVEN': 1, 'ODD': 2}


def test_a_source_level_the_file_has_not_got_is_reported(folder):
    source_file(folder / 'levs.txt', SOURCE + [('X.000009', 99999.0, '1', 'e')])
    done = run('write_IDEN_level_ids.py', folder, 'levs.txt',
               '--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity')
    assert done.returncode == 0, done.stderr
    assert 'unmatched: X.000009' in done.stdout
    assert 'X.000009' not in open(str(folder / levids.IDEN_IDS)).read()


def test_a_calculated_level_is_never_given_an_identifier(folder):
    """Level 1 has no star, so nothing may be tied to it."""
    source_file(folder / 'levs.txt', SOURCE + [('X.000001', 30000.0, '2', 'e')])
    run('write_IDEN_level_ids.py', folder, 'levs.txt',
        '--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity')
    assert 1 not in levids.read_iden_ids(
        str(folder / levids.IDEN_IDS)).values()


# --------------------------------------------------------------------------
# update_IDEN
# --------------------------------------------------------------------------

#: One line per transition of the little spectrum, in the descending
#: wavenumber order dlv.dat is in.
OBSERVED = [(100, 25001.501), (90, 15001.499), (80, 10000.0),
            (70, 5000.0)]

#: The LOPT level table: the two levels that moved and the two that did
#: not.  `D2tot` is the larger uncertainty on level 4, so it is the one
#: that has to come through.
LOPT_LEVELS = [
    ('X.000002', 25001.500, 0.020, 0.010, 0.000, 0.010),
    ('X.000003', 20000.400, 0.030, 0.000, 0.000, 0.000),
    ('X.000004', 10000.000, 0.010, 0.040, 0.000, 0.040),
    ('X.000005', 0.000, 0.000, 0.000, 0.000, 0.000),
]

LOPT_LEVEL_HEADER = ('Designation', 'Energy', 'D1', 'D2stat', 'D2sys',
                     'D2tot', 'N_lines', 'Comments')
LOPT_LINE_HEADER = ('Iobs', 'wn_o', 'L1', 'L2', 'E1', 'E2', 'F', 'Weight')


def lopt_files(folder, lines):
    table(folder / 'lev.txt', LOPT_LEVEL_HEADER,
          [row + ('3', '') for row in LOPT_LEVELS])
    table(folder / 'lin.txt', LOPT_LINE_HEADER, lines)
    dlv_file(folder / 'dlv.dat', OBSERVED)


@pytest.fixture
def spectrum(folder):
    """The whole working directory: the two decks and the table."""
    trans_file(folder / iden.TRANS, [
        level_card(2, 3.0, 'd2   _3F3F ', 25000.0, True),
        line_card(4, 12, 10000.0, True, 15000.0),
        line_card(5, 20, 0.0, True, 25000.0),
        level_card(3, 2.0, 'd2   _3F3F ', 20000.0, True),
        line_card(4, 15, 10000.0, True, 10000.0, 80, 10000.0, 0.0, 3),
        level_card(4, 4.0, 'f2   ~3H3H ', 10000.0, True),
        line_card(1, 8, 30000.0, False, 20000.0),
    ])
    run('write_IDEN_level_ids.py', folder, 'levs.txt',
        '--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity')
    return folder


def cards(path):
    return [line.rstrip('\r\n')
            for line in open(str(path)) if line.strip()]


def test_the_fixed_levels_take_their_lopt_energies(spectrum):
    lopt_files(spectrum, [(10, 15001.5, 'X.000004', 'X.000002',
                           10000.0, 25001.5, '', 1.0)])
    done = run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    assert done.returncode == 0, done.stderr

    levels = enlev_levels(spectrum)
    assert levels[2].e_fixed == 25001.5
    assert levels[2].uncertainty == 0.020        # D1 is the larger here
    assert levels[4].uncertainty == 0.040        # and D2tot here
    # The calculated level was not named by LOPT: its energy is left alone.
    assert levels[1].e_fixed == 30000.0


def enlev_levels(folder):
    return dict((one.number, one)
                for one in iden.read_enlev(str(folder / iden.ENLEV)))


#: The window a calculated d2 level is given: the rms of E_fixed - E_calc
#: over the two experimental d2 levels *after* LOPT has moved them,
#: 25001.5 - 25010 and 20000.4 - 20020, as ENLEV.DAT prints it.
D2_RMS = round(math.sqrt((8.5 ** 2 + 19.6 ** 2) / 2), 3)


def test_a_calculated_level_takes_the_rms_of_its_configuration(spectrum):
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt')
    assert done.returncode == 0, done.stderr
    assert enlev_levels(spectrum)[1].uncertainty == D2_RMS == 15.106
    # Two experimental levels are fewer than THIN, and the report says so.
    assert 'thin' in done.stdout


def test_ucalc_leaves_a_known_configuration_alone(spectrum):
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt', '--ucalc', '250')
    assert done.returncode == 0, done.stderr
    levels = enlev_levels(spectrum)
    assert levels[1].uncertainty == D2_RMS
    assert levels[2].uncertainty == 0.020


def test_with_nosd_ucalc_is_given_to_every_calculated_level(spectrum):
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt', '--ucalc', '250',
               '--nosd')
    assert done.returncode == 0, done.stderr
    levels = enlev_levels(spectrum)
    assert levels[1].uncertainty == 250.0
    assert levels[2].uncertainty == 0.020
    # Without the line table TRANS.DAT is not touched at all.
    assert not (spectrum / (iden.TRANS + '.bak')).exists()


def test_with_nosd_alone_a_calculated_level_keeps_its_uncertainty(spectrum):
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt', '--nosd')
    assert done.returncode == 0, done.stderr
    assert enlev_levels(spectrum)[1].uncertainty == 5000.0
    assert '--nosd' in done.stdout


def unknown_configuration(folder):
    """Level 1 moved into a configuration no experimental level is in."""
    enlev_file(folder / iden.ENLEV,
               [LEVELS[0][:6] + ('p2   _3P2P ',)] + LEVELS[1:])


def test_an_unknown_configuration_keeps_its_uncertainty(spectrum):
    unknown_configuration(spectrum)
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt')
    assert done.returncode == 0, done.stderr
    assert enlev_levels(spectrum)[1].uncertainty == 5000.0
    assert 'p2' in done.stdout


def test_an_unknown_configuration_takes_ucalc(spectrum):
    unknown_configuration(spectrum)
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt', '--ucalc', '250')
    assert done.returncode == 0, done.stderr
    assert enlev_levels(spectrum)[1].uncertainty == 250.0


def test_a_zero_rms_is_no_window(spectrum):
    # Both experimental d2 levels left exactly where RCE put them: the
    # rms is 0, which IDEN2 would read as no window at all, so the
    # configuration counts as unknown.
    lopt_files(spectrum, [])
    table(spectrum / 'lev.txt', LOPT_LEVEL_HEADER,
          [('X.000002', 25010.0) + LOPT_LEVELS[0][2:] + ('3', ''),
           ('X.000003', 20020.0) + LOPT_LEVELS[1][2:] + ('3', '')]
          + [row + ('3', '') for row in LOPT_LEVELS[2:]])
    done = run('update_IDEN.py', spectrum, 'lev.txt', '--ucalc', '250')
    assert done.returncode == 0, done.stderr
    assert enlev_levels(spectrum)[1].uncertainty == 250.0


def test_ucalc_wants_a_positive_number(spectrum):
    lopt_files(spectrum, [])
    done = run('update_IDEN.py', spectrum, 'lev.txt', '--ucalc', '0')
    assert done.returncode == 2
    assert '--ucalc' in done.stderr


def test_a_used_line_is_assigned_and_an_unused_one_is_not(spectrum):
    """The set of assignments becomes the set the optimization used.

    The 4 -> 2 transition was in the file with no line and is used, so it
    gains one; the 4 -> 3 transition was assigned line 3 and LOPT gave it
    a zero weight, so it loses it.
    """
    lopt_files(spectrum, [
        (10, 15001.5, 'X.000004', 'X.000002', 10000.0, 25001.5, '', 1.0),
        (20, 10000.0, 'X.000004', 'X.000003', 10000.0, 20000.4, 'P', 0.0),
        (30, 25001.5, 'X.000005', 'X.000002', 0.0, 25001.5, '', 0.5),
    ])
    done = run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    assert done.returncode == 0, done.stderr

    assigned = dict(((one.upper, one.lower), one.line)
                    for one in iden.read_trans(str(spectrum / iden.TRANS))
                    if one.kind == '+')
    assert assigned[(2, 4)] == 2                 # 15001.5 -> the second row
    assert assigned[(2, 5)] == 1                 # 25001.5 -> the first row
    assert assigned[(3, 4)] == 0                 # the weight was zero
    # The transition down to the calculated level is none of LOPT's
    # business and keeps what it had.
    assert assigned[(4, 1)] == 0


def test_every_wavenumber_is_looked_for_in_dlv(spectrum):
    """The whole line table, not only the rows that reach a transition.

    The assignment written into TRANS.DAT is a row number of dlv.dat, so
    the two have to be the pair that belong together.  A row LOPT did not
    use, and a row naming levels the table has not got, would never reach
    the assignment pass, so neither would say anything about a dlv.dat
    that is the wrong file -- hence a check over every row.
    """
    lopt_files(spectrum, [
        (10, 15001.5, 'X.000004', 'X.000002', 10000.0, 25001.5, '', 1.0),
        (20, 10000.0, 'X.000004', 'X.000003', 10000.0, 20000.4, 'P', 0.0),
        (30, 25001.5, 'X.000005', 'X.000002', 0.0, 25001.5, '', 0.5),
    ])
    done = run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    assert done.returncode == 0, done.stderr
    assert '3 of the 3 wavenumbers' in done.stdout
    # 15001.5 against the 15001.499 of dlv.dat is the worst of the three.
    assert '0.0010 cm-1 from its nearest row' in done.stdout
    assert 'no row of dlv.dat within the tolerance' not in done.stdout


def test_a_wavenumber_no_dlv_row_matches_is_reported(spectrum):
    lopt_files(spectrum, [
        (10, 15001.5, 'X.000004', 'X.000002', 10000.0, 25001.5, '', 1.0),
        # Unused, and between levels the table has not got: it reaches
        # the assignment pass by neither road, and is checked anyway.
        (20, 12345.6, 'Y.000001', 'Y.000002', 0.0, 12345.6, 'P', 0.0),
    ])
    done = run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    assert done.returncode == 0, done.stderr
    assert '1 of the 2 wavenumbers' in done.stdout
    assert ('1 wavenumber has no row of dlv.dat within the tolerance'
            in done.stdout)
    # It is a warning and not a refusal: TRANS.DAT is still written.
    assert (spectrum / (iden.TRANS + '.bak')).exists()


def test_wntol_widens_the_search(spectrum):
    lopt_files(spectrum, [
        (10, 15001.9, 'X.000004', 'X.000002', 10000.0, 25001.5, '', 1.0),
    ])
    done = run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    assert '0 of the 1 wavenumbers' in done.stdout
    done = run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt',
               '--wntol', '0.5')
    assert '1 of the 1 wavenumbers' in done.stdout
    assigned = dict(((one.upper, one.lower), one.line)
                    for one in iden.read_trans(str(spectrum / iden.TRANS))
                    if one.kind == '+')
    assert assigned[(2, 4)] == 2


def test_a_cleared_assignment_takes_its_observation_with_it(spectrum):
    lopt_files(spectrum, [(20, 10000.0, 'X.000004', 'X.000003',
                           10000.0, 20000.4, 'P', 0.0)])
    run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    card = [one for one in cards(spectrum / iden.TRANS)
            if one.startswith('+   4') and one[5:10].strip() == '15'][0]
    assert card[37:43].strip() == '0'            # the observed intensity
    assert float(card[43:55]) == 0.0             # the wavenumber
    assert float(card[55:66]) == 0.0             # dwn


def test_the_transition_energies_follow_the_new_levels(spectrum):
    lopt_files(spectrum, [(10, 15001.5, 'X.000004', 'X.000002',
                           10000.0, 25001.5, '', 1.0)])
    run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    records = iden.read_trans(str(spectrum / iden.TRANS))
    group = [one for one in records if one.kind == '$' and one.upper == 2][0]
    assert float(group.card[28:41]) == 25001.5
    row = [one for one in records
           if one.kind == '+' and (one.upper, one.lower) == (2, 4)][0]
    assert float(row.card[24:37]) == 15001.5
    assert float(row.card[55:66]) == pytest.approx(-0.001, abs=1e-6)


def test_a_second_run_changes_nothing(spectrum):
    lopt_files(spectrum, [
        (10, 15001.5, 'X.000004', 'X.000002', 10000.0, 25001.5, '', 1.0),
        (20, 10000.0, 'X.000004', 'X.000003', 10000.0, 20000.4, 'P', 0.0),
    ])
    run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt', '--ucalc', '250')
    first = (open(str(spectrum / iden.ENLEV), 'rb').read(),
             open(str(spectrum / iden.TRANS), 'rb').read())
    run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt', '--ucalc', '250')
    assert (open(str(spectrum / iden.ENLEV), 'rb').read(),
            open(str(spectrum / iden.TRANS), 'rb').read()) == first


def test_the_previous_files_are_kept(spectrum):
    lopt_files(spectrum, [(10, 15001.5, 'X.000004', 'X.000002',
                           10000.0, 25001.5, '', 1.0)])
    before = open(str(spectrum / iden.ENLEV), 'rb').read()
    run('update_IDEN.py', spectrum, 'lev.txt', 'lin.txt')
    assert open(str(spectrum / (iden.ENLEV + '.bak')), 'rb').read() == before
    assert (spectrum / (iden.TRANS + '.bak')).is_file()


def test_a_table_from_another_enlev_is_refused(spectrum):
    lopt_files(spectrum, [])
    levids.write_iden_ids([('X.000002', 99)],
                          str(spectrum / levids.IDEN_IDS))
    done = run('update_IDEN.py', spectrum, 'lev.txt')
    assert done.returncode == 1
    assert 'different ENLEV.DAT' in done.stderr


def test_a_table_that_is_not_a_bijection_is_refused(spectrum):
    lopt_files(spectrum, [])
    levids.write_iden_ids([('X.000002', 2), ('X.000003', 2)],
                          str(spectrum / levids.IDEN_IDS))
    done = run('update_IDEN.py', spectrum, 'lev.txt')
    assert done.returncode == 1
    assert 'both' in done.stderr


def test_a_missing_column_says_which_one(spectrum):
    lopt_files(spectrum, [])
    table(spectrum / 'lev.txt', ('Designation', 'Energy'),
          [('X.000002', 25001.5)])
    done = run('update_IDEN.py', spectrum, 'lev.txt')
    assert done.returncode == 1
    assert 'D1' in done.stderr


def test_both_tools_print_their_usage(tmp_path):
    for tool in ('write_IDEN_level_ids.py', 'update_IDEN.py'):
        done = run(tool, tmp_path, '--help')
        assert done.returncode == 0
        assert tool.replace('.py', ' ') in done.stdout
