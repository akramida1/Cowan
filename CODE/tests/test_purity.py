"""Purity per configuration, and the file the choice is marked in.

Phase 12 has no Perl ancestor, so its acceptance is not a diff.  What is
pinned here is the pipeline that carries the choice: the `conf0` column
printout now writes, the grouping that column makes possible, and the
one file of the whole chain that a person edits -- written, marked by
hand in every column of its field, read back, and merged.
"""

import os
import subprocess
import sys

import pytest

from conftest import case
from toolbox import purity
from toolbox.purity import PurityError

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(script, folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, script)] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


@pytest.fixture(scope='module')
def pr3(tmp_path_factory):
    """A Pr III levels.txt with conf0, and the deck beside it."""
    folder = case('Pr3')
    out = tmp_path_factory.mktemp('purity')
    for name in ('IN36',):
        with open(os.path.join(str(folder), name), 'rb') as f:
            data = f.read()
        with open(str(out / name), 'wb') as f:
            f.write(data)
    done = run('printout.py', folder,
               str(out / 'levels.txt'), str(out / 'params.txt'),
               '--couplings', 'LS,JJ', '--dominant-config')
    assert done.returncode == 0, done.stderr
    return out


# --------------------------------------------------------------------------
# the conf0 column
# --------------------------------------------------------------------------

def test_conf0_is_off_by_default(tmp_path):
    """levels.txt is unchanged unless the column is asked for."""
    done = run('printout.py', case('W70'), str(tmp_path / 'levels.txt'),
               str(tmp_path / 'p.txt'), '--couplings', 'LS')
    assert done.returncode == 0, done.stderr
    with open(str(tmp_path / 'levels.txt')) as f:
        header = f.readline().rstrip('\n').split('\t')
    assert 'conf0' not in header
    assert header[:8] == purity.FIXED


def test_conf0_column_is_the_invariant_configuration(pr3):
    """Every level of a coupling agrees with the same level of the other.

    RCE picks the dominant configuration once, before it prints
    anything, so the LS block and the JJ block of one level must name
    the same configuration -- which is the invariance the whole phase
    rests on, checked against the file rather than assumed.
    """
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    assert table.has_conf0
    assert table.couplings == ['LS', 'JJ']
    seen = {}
    for row in table.rows:
        key = (row.parity, row.e_calc, row.j)
        if key in seen:
            assert seen[key] == row.conf0
        seen[key] = row.conf0
    assert len(seen) * 2 == len(table.rows)


def test_conf0_is_a_configuration_of_the_deck(pr3):
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    rows = purity.deck_configurations(str(pr3))
    _, _, nameless = purity.purity_table(table, rows)
    assert nameless == []
    labels = {r.label for r in rows}
    assert '4f3.5p6' in labels
    # `p5f3d` in ING11, cut to six characters by RCE's A6 field and put
    # back together from the LEVELS header.
    assert '4f3.5p5.5d' in labels


def test_a_levels_file_without_conf0_says_what_to_run(tmp_path):
    done = run('printout.py', case('W70'), str(tmp_path / 'levels.txt'),
               str(tmp_path / 'p.txt'), '--couplings', 'LS')
    assert done.returncode == 0
    table = purity.read_levels_txt(str(tmp_path / 'levels.txt'))
    with pytest.raises(PurityError) as problem:
        purity.purity_table(table, purity.deck_configurations(str(case('W70'))))
    assert '--dominant-config' in str(problem.value)


# --------------------------------------------------------------------------
# block aliasing: two blocks of one coupling in one file
# --------------------------------------------------------------------------

def test_two_blocks_of_one_coupling_need_a_name_each(tmp_path):
    done = run('printout.py', case('U90'), str(tmp_path / 'levels.txt'),
               str(tmp_path / 'p.txt'), '--couplings', 'LS,LS')
    assert done.returncode == 2
    assert 'twice' in done.stderr


def test_an_aliased_block_keeps_its_own_name(tmp_path):
    done = run('printout.py', case('U90'), str(tmp_path / 'levels.txt'),
               str(tmp_path / 'p.txt'), '--couplings', 'LS,LSr',
               '--from', 'LSr=LS@reordered', '--dominant-config')
    assert done.returncode == 0, done.stderr
    table = purity.read_levels_txt(str(tmp_path / 'levels.txt'))
    assert table.couplings == ['LS', 'LSr']
    # The alias is a second reading of the same levels, so it groups on
    # the same configurations.
    assert ({r.conf0 for r in table.rows if r.coupling == 'LS'}
            == {r.conf0 for r in table.rows if r.coupling == 'LSr'})


# --------------------------------------------------------------------------
# the purity table
# --------------------------------------------------------------------------

def test_purity_is_the_unweighted_mean_of_the_leading_component(pr3):
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    kept, empty, _ = purity.purity_table(
        table, purity.deck_configurations(str(pr3)))
    assert kept and not [r for r in kept if not r.nlev]
    row = [r for r in kept if r.label == '4f3.5p6'][0]
    for coupling in table.couplings:
        want = [level.leading.percent for level in table.rows
                if level.coupling == coupling and level.conf0 == row.label]
        assert row.nlev == len(want)
        assert row.purity[coupling] == int(sum(want) / float(len(want)) + 0.5)
    assert empty == [] or all(r.nlev == 0 for r in empty)


def test_the_deck_order_is_kept(pr3):
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    rows = purity.deck_configurations(str(pr3))
    kept, _, _ = purity.purity_table(table, rows)
    assert [r.index for r in kept] == sorted(r.index for r in kept)


def test_jj_wins_where_the_open_shells_are_far_apart(pr3):
    """4f2.6p is a jj configuration and 4f.5d2 is an LS one.

    Not a formatting check: it is the one place the numbers are asserted
    to mean something, and both are configurations whose coupling is not
    in doubt.
    """
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    kept, _, _ = purity.purity_table(
        table, purity.deck_configurations(str(pr3)))
    by_label = {r.label: r for r in kept}
    assert (by_label['4f2.5p6.6p'].purity['JJ']
            > by_label['4f2.5p6.6p'].purity['LS'])
    assert (by_label['4f.5p6.5d2'].purity['LS']
            > by_label['4f.5p6.5d2'].purity['JJ'])


# --------------------------------------------------------------------------
# coupling_choice.txt
# --------------------------------------------------------------------------

def make_rows():
    return [purity.ConfigRow(index=1, parity='o', label='4f3', deck_label='f3',
                             nlev=41, purity={'LS': 78, 'JJ': 60}),
            purity.ConfigRow(index=2, parity='e', label='4f2.6p',
                             deck_label='f26p', nlev=67,
                             purity={'LS': 49, 'JJ': 64})]


def write(tmp_path, rows, couplings=('LS', 'JJ')):
    path = str(tmp_path / purity.CHOICE)
    for row in rows:
        row.chosen = purity.suggest(row, list(couplings))
    purity.write_choice(path, rows, list(couplings))
    return path


def test_the_star_goes_on_the_highest_purity(tmp_path):
    path = write(tmp_path, make_rows())
    found, couplings = purity.read_choice(path)
    assert couplings == ['LS', 'JJ']
    assert found == {'4f3': 'LS', '4f2.6p': 'JJ'}


@pytest.mark.parametrize('column', range(purity.FIELD))
def test_a_star_anywhere_in_the_field_selects_it(tmp_path, column):
    """A hand-placed star need not sit where the writer put one."""
    path = write(tmp_path, make_rows())
    lines = open(path).read().splitlines()
    head = [line for line in lines if line.startswith('par ')][0]
    start = head.index('JJ')
    at = [i for i, line in enumerate(lines) if line.startswith('o   4f3')][0]
    line = lines[at].ljust(start + purity.FIELD)
    line = line[:start].replace('*', ' ') + line[start:]
    lines[at] = line[:start + column] + '*' + line[start + column + 1:]
    open(path, 'w').write('\n'.join(lines) + '\n')
    found, _ = purity.read_choice(path)
    assert found['4f3'] == 'JJ'


def test_no_star_stops_and_names_the_line(tmp_path):
    path = write(tmp_path, make_rows())
    lines = open(path).read().splitlines()
    at = [i for i, line in enumerate(lines) if line.startswith('o   4f3')][0]
    lines[at] = lines[at].replace('*', ' ')
    open(path, 'w').write('\n'.join(lines) + '\n')
    with pytest.raises(PurityError) as problem:
        purity.read_choice(path)
    assert '4f3' in str(problem.value)
    assert 'no scheme is starred' in str(problem.value)


def test_two_stars_stop_and_name_both(tmp_path):
    path = write(tmp_path, make_rows())
    lines = open(path).read().splitlines()
    head = [line for line in lines if line.startswith('par ')][0]
    start = head.index('JJ')
    at = [i for i, line in enumerate(lines) if line.startswith('o   4f3')][0]
    line = lines[at].ljust(start + purity.FIELD)
    lines[at] = line[:start] + '*' + line[start + 1:]
    open(path, 'w').write('\n'.join(lines) + '\n')
    with pytest.raises(PurityError) as problem:
        purity.read_choice(path)
    assert 'LS' in str(problem.value) and 'JJ' in str(problem.value)


def test_merge_keeps_the_stars_while_the_numbers_change(tmp_path):
    rows = make_rows()
    path = write(tmp_path, rows)
    # The author overrules the suggestion on the first configuration.
    lines = open(path).read().splitlines()
    head = [line for line in lines if line.startswith('par ')][0]
    start = head.index('JJ')
    at = [i for i, line in enumerate(lines) if line.startswith('o   4f3')][0]
    line = lines[at].ljust(start + purity.FIELD).replace('*', ' ')
    lines[at] = line[:start] + '*' + line[start + 1:]
    open(path, 'w').write('\n'.join(lines) + '\n')

    fresh = make_rows()
    fresh[0].purity = {'LS': 71, 'JJ': 55}      # a new fit, new numbers
    for row in fresh:
        row.chosen = purity.suggest(row, ['LS', 'JJ'])
    gone = purity.merge_choice(path, fresh, ['LS', 'JJ'])
    assert gone == []
    assert fresh[0].chosen == 'JJ'              # the author's mark survives
    purity.write_choice(path, fresh, ['LS', 'JJ'])
    found, _ = purity.read_choice(path)
    assert found['4f3'] == 'JJ'
    assert ' 71' in open(path).read()


def test_merge_reports_what_left_the_deck(tmp_path):
    path = write(tmp_path, make_rows())
    fresh = make_rows()[:1]
    gone = purity.merge_choice(path, fresh, ['LS', 'JJ'])
    assert gone == ['4f2.6p']


# --------------------------------------------------------------------------
# the script
# --------------------------------------------------------------------------

def test_the_script_writes_and_will_not_overwrite(pr3):
    done = run('config_purity.py', pr3, '--levels', 'levels.txt')
    assert done.returncode == 0, done.stderr
    path = str(pr3 / purity.CHOICE)
    assert os.path.exists(path)
    first = open(path).read()

    again = run('config_purity.py', pr3, '--levels', 'levels.txt')
    assert again.returncode == 2
    assert '--merge' in again.stderr and '--force' in again.stderr
    assert open(path).read() == first

    forced = run('config_purity.py', pr3, '--levels', 'levels.txt', '--force')
    assert forced.returncode == 0, forced.stderr
    assert '--force' not in run('config_purity.py', pr3, '--levels',
                                'levels.txt', '--merge').stderr


def test_force_and_merge_are_not_both(pr3):
    done = run('config_purity.py', pr3, '--levels', 'levels.txt',
               '--force', '--merge')
    assert done.returncode == 2
    assert 'one of them' in done.stderr


def test_intermediates_come_off_a_component_configuration():
    assert purity.strip_intermediates('4f.5p6.5d2.(3F)') == '4f.5p6.5d2'
    assert purity.strip_intermediates('4f2.(3P*<0>).6p') == '4f2.6p'
    assert purity.strip_intermediates('4f3') == '4f3'
    assert purity.Component(44, '4f2.(3P*<0>).6p', '(0,3/2)*').bare \
        == '4f2.6p'


def test_a_true_jj_component_names_the_deck_configuration():
    """`2p-2.2p+` is the deck's `2p3`, and W70 is full of both.

    True jj splits the open subshells, so its component configurations
    are not the deck's until the halves are added back together --
    without which no scheme could be chosen for them and phase 13 would
    stop on every W70 run.
    """
    assert purity.canonical_configuration('2p-2.2p+') == '2p3'
    assert purity.canonical_configuration('2p+3') == '2p3'
    assert purity.canonical_configuration('2s.2p-.2p+') == '2s.2p2'
    assert purity.canonical_configuration('2s.2p+.3d') == '2s.2p.3d'
    assert purity.Component(44, '2p-2.(0).2p+', '(0,3/2)*').bare == '2p3'


def test_a_reordered_deck_names_the_same_configuration():
    """`2p2.1s` is `1s.2p2`: U90/reordered writes summation order.

    A basis state's tokens are printed in the deck's summation order
    wherever an intermediate makes that order visible (`designate.py`),
    which is right for the state and wrong for the name of the
    configuration it belongs to.
    """
    assert purity.canonical_configuration('2p2.1s') == '1s.2p2'
    assert purity.canonical_configuration('1s.2p2') == '1s.2p2'
    assert purity.Component(44, '2p2.(3P).1s', '4P').bare == '1s.2p2'


def test_an_unrecognised_label_is_left_alone():
    assert purity.canonical_configuration('') == ''
    assert purity.canonical_configuration('18.58h') == '18.58h'
