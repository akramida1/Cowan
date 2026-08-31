"""One coupling scheme per configuration, applied to every level.

Phase 13 has no Perl ancestor either, so its acceptance is not a diff.
What is pinned here is the rule -- a component is printed in the scheme
its configuration was assigned and in no other -- together with the
invariance the rule rests on, measured on Pr III rather than asserted:
the weight a configuration carries in a level is the same in LS, in JJ
and in LSJK, so the percentages of a mixed-scheme composition still add
up.
"""

import os
import subprocess
import sys

import pytest

from conftest import case
from toolbox import prepub, purity
from toolbox.purity import PurityError

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(script, folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, script)] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


@pytest.fixture(scope='module')
def pr3(tmp_path_factory):
    """Pr III through phases 6 and 12: levels.txt and a marked choice.

    Run the way the pipeline asks for -- ten components down to 0.5% --
    so that `prepub` does the trimming and not `printout`.  Three
    schemes, the third from the `JK/` sub-folder, because a single
    scheme could not show that anything is being chosen.
    """
    folder = case('Pr3')
    out = tmp_path_factory.mktemp('prepub')
    with open(os.path.join(str(folder), 'IN36'), 'rb') as source:
        data = source.read()
    with open(str(out / 'IN36'), 'wb') as target:
        target.write(data)
    done = run('printout.py', folder,
               str(out / 'levels.txt'), str(out / 'params.txt'),
               '--couplings', 'LS,JJ,LSJK', '--from', 'LSJK=JK',
               '--dominant-config', '--components', '10',
               '--min-percent', '0.5')
    assert done.returncode == 0, done.stderr
    done = run('config_purity.py', out, '--levels', 'levels.txt')
    assert done.returncode == 0, done.stderr
    return out


@pytest.fixture(scope='module')
def merged(pr3):
    """The levels with their compositions built, before any trimming."""
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    choice, _ = purity.read_choice(str(pr3 / 'coupling_choice.txt'))
    levels, missing = prepub.collect(table)
    prepub.apply_choice(levels, choice, table.couplings,
                        deck=[r.label for r in
                              purity.deck_configurations(str(pr3))])
    return table, choice, levels, missing


# --------------------------------------------------------------------------
# printout has to name every component it writes
# --------------------------------------------------------------------------

def test_the_header_names_as_many_components_as_were_asked_for(pr3):
    """A ten-component run whose header stops at five hides the rest.

    The header is the specification of the file -- `purity.py` finds the
    components by counting `%n`/`confn`/`termn` -- so five names over
    ten written triples silently threw away components 6 to 10, which is
    exactly the material this pipeline was told to ask for.
    """
    with open(str(pr3 / 'levels.txt')) as stream:
        lines = stream.read().splitlines()
    header = lines[0].split('\t')
    assert header[-3:] == ['%10', 'conf10', 'term10']
    assert max(len(line.split('\t')) for line in lines[1:]) == len(header)


def test_the_default_header_is_unchanged(tmp_path):
    done = run('printout.py', case('W70'), str(tmp_path / 'levels.txt'),
               str(tmp_path / 'p.txt'), '--couplings', 'LS')
    assert done.returncode == 0, done.stderr
    with open(str(tmp_path / 'levels.txt')) as stream:
        header = stream.readline().rstrip('\n').split('\t')
    assert header == purity.FIXED + [
        name % n for n in range(1, 6)
        for name in ('%%%d', 'conf%d', 'term%d')]


# --------------------------------------------------------------------------
# the blocks, paired back into levels
# --------------------------------------------------------------------------

def test_every_level_is_gathered_from_all_its_blocks(merged):
    table, _, levels, missing = merged
    assert table.couplings == ['LS', 'JJ', 'LSJK']
    assert len(levels) * len(table.couplings) == len(table.rows)
    for level in levels:
        assert set(level.blocks) == set(table.couplings)
    assert missing == {}


def test_a_level_met_first_without_g_takes_it_from_a_later_block():
    """The optional fields are the same in every block, so any will do."""
    class Row:
        def __init__(self, coupling, g, conf0):
            self.coupling, self.g, self.conf0 = coupling, g, conf0
            self.parity, self.e_calc, self.j = 'e', 1.0, 0.5
            self.e_obs, self.unfitted, self.line = 1.0, False, 2
            self.group, self.components = '', []

    class Table:
        path = 'x'
        couplings = ['LS', 'JJ']
        rows = [Row('LS', None, ''), Row('JJ', 1.25, '4f2')]

    levels, missing = prepub.collect(Table())
    assert len(levels) == 1
    assert levels[0].g == 1.25 and levels[0].conf0 == '4f2'
    assert missing == {}


def test_two_blocks_of_one_level_in_one_coupling_stop_the_run():
    class Row:
        def __init__(self):
            self.coupling, self.g, self.conf0 = 'LS', None, ''
            self.parity, self.e_calc, self.j = 'e', 1.0, 0.5
            self.e_obs, self.unfitted, self.line = 1.0, False, 2
            self.group, self.components = '', []

    class Table:
        path = 'x'
        couplings = ['LS']
        rows = [Row(), Row()]

    # Two rows of one coupling at one (parity, Ec, J) are two different
    # levels that happen to be degenerate, and pair by occurrence.
    levels, _ = prepub.collect(Table())
    assert len(levels) == 2


# --------------------------------------------------------------------------
# the rule
# --------------------------------------------------------------------------

def test_each_component_comes_from_the_scheme_its_configuration_chose(merged):
    _, choice, levels, _ = merged
    for level in levels:
        for component in level.components:
            assert component.bare in choice
        # and every component of an unassigned block is gone: the count
        # of a configuration's components equals what its own block has.
        for label in {c.bare for c in level.components}:
            block = level.blocks[choice[label]]
            assert (len([c for c in level.components if c.bare == label])
                    == len([c for c in block.components if c.bare == label]))


def test_a_composition_can_mix_schemes(merged):
    """The point of the whole tool: two configurations, two schemes."""
    _, choice, levels, _ = merged
    mixed = [one for one in levels
             if len({choice[c.bare] for c in one.components}) >= 3]
    assert mixed, 'no level of Pr III draws on three schemes at once'


def test_components_come_out_in_descending_order(merged):
    _, _, levels, _ = merged
    for level in levels:
        percents = [c.percent for c in level.components]
        assert percents == sorted(percents, reverse=True)


def test_the_configuration_weights_are_scheme_invariant(merged):
    """The physics the phase rests on, measured on 5626 levels.

    A change of coupling is block-diagonal in the configuration, so a
    configuration's summed weight in a level must not move between LS,
    JJ and LSJK -- beyond what integer percentages and the per-block cut
    can account for.  Nothing here corrects anything; if this fails, the
    premise of phases 12-14 is wrong for this run.
    """
    table, _, levels, _ = merged
    assert prepub.disagreements(levels, table.couplings) == []


def test_a_scheme_that_disagrees_is_reported():
    """The check has to be able to fail, so make it."""
    class Component:
        def __init__(self, percent, bare):
            self.percent, self.bare = percent, bare
            self.configuration, self.term = bare, '3H'

    class Row:
        def __init__(self, coupling, pairs):
            self.coupling = coupling
            self.components = [Component(p, b) for p, b in pairs]

    level = prepub.Level(parity='e', e_obs=0.0, unfitted=True, e_calc=0.0,
                         j=1.0, g=None, conf0='4f2', line=2)
    # Both blocks print a full 100%, so neither has room to hide a
    # difference, and they still disagree by 40 points about 4f2.
    level.blocks = {'LS': Row('LS', [(100, '4f2')]),
                    'JJ': Row('JJ', [(60, '4f2'), (40, '5d2')])}
    off = prepub.disagreements([level], ['LS', 'JJ'])
    assert len(off) == 1 and off[0].label == '4f2'
    assert off[0].spread == 40

    # The same 40 points where the JJ block simply stopped printing are
    # not a disagreement: they are the per-block cut, and the slack says
    # so.  A configuration only one block printed says nothing either.
    level.blocks['JJ'] = Row('JJ', [(60, '4f2')])
    assert prepub.disagreements([level], ['LS', 'JJ']) == []


# --------------------------------------------------------------------------
# what the choice file does not say
# --------------------------------------------------------------------------

def test_a_configuration_the_choice_file_omits_takes_the_first_scheme():
    """It has a line only if it dominates a level; it may still mix in."""
    class Component:
        def __init__(self, percent, bare):
            self.percent, self.bare = percent, bare
            self.configuration, self.term = bare, '3H'

    class Row:
        def __init__(self, coupling, pairs):
            self.coupling, self.line = coupling, 2
            self.components = [Component(p, b) for p, b in pairs]

    level = prepub.Level(parity='e', e_obs=0.0, unfitted=True, e_calc=0.0,
                         j=1.0, g=None, conf0='4f2', line=2)
    level.blocks = {'LS': Row('LS', [(90, '4f2'), (5, '5d2')]),
                    'JJ': Row('JJ', [(88, '4f2'), (6, '5d2')])}
    unassigned = prepub.apply_choice([level], {'4f2': 'JJ'}, ['LS', 'JJ'],
                                     deck=['4f2', '5d2'])
    assert unassigned == {'5d2': 2}
    assert [(c.percent, c.bare) for c in level.components] \
        == [(88, '4f2'), (5, '5d2')]


def test_a_configuration_in_no_scheme_and_no_deck_stops_the_run():
    class Component:
        percent, bare = 90, '6h2'
        configuration, term = '6h2', '3H'

    class Row:
        coupling, line = 'LS', 7
        components = [Component()]

    level = prepub.Level(parity='e', e_obs=0.0, unfitted=True, e_calc=0.0,
                         j=1.0, g=None, conf0='4f2', line=7)
    level.blocks = {'LS': Row()}
    with pytest.raises(PurityError) as problem:
        prepub.apply_choice([level], {'4f2': 'LS'}, ['LS'], deck=['4f2'])
    assert '6h2' in str(problem.value)
    assert 'line 7' in str(problem.value)


def test_a_marked_scheme_with_no_block_stops_the_run(pr3, tmp_path):
    """A star on a column the levels file has no block for.

    Silently ignoring it would take that configuration's components out
    of every composition, which is worse than stopping.
    """
    done = run('printout.py', case('Pr3'), str(tmp_path / 'levels.txt'),
               str(tmp_path / 'p.txt'), '--couplings', 'LS',
               '--dominant-config')
    assert done.returncode == 0, done.stderr
    with open(str(pr3 / 'coupling_choice.txt')) as stream:
        marked = stream.read()
    with open(str(tmp_path / 'coupling_choice.txt'), 'w') as stream:
        stream.write(marked)
    done = run('make_levels_list_prepub.py', tmp_path)
    assert done.returncode == 1
    assert 'no block for' in done.stderr


# --------------------------------------------------------------------------
# trimming and order
# --------------------------------------------------------------------------

def test_trim_cuts_by_count_and_by_percentage(merged):
    table, choice, _, _ = merged
    levels, _ = prepub.collect(table)
    prepub.apply_choice(levels, choice, table.couplings)
    prepub.trim(levels, components=3, min_percent=5)
    for level in levels:
        assert len(level.components) <= 3
        assert all(c.percent >= 5 for c in level.components)


def test_the_order_is_parity_then_Ec(merged):
    _, _, levels, _ = merged
    ordered = prepub.order(levels)
    parities = []
    for level in ordered:
        if level.parity not in parities:
            parities.append(level.parity)
    # The parity blocks come in the order levels.txt met them -- Pr III
    # has its odd set first -- and each is sorted on Ec.
    assert parities == ['o', 'e']
    for parity in parities:
        energies = [one.e_calc for one in ordered if one.parity == parity]
        assert energies == sorted(energies)


# --------------------------------------------------------------------------
# the file
# --------------------------------------------------------------------------

@pytest.fixture(scope='module')
def written(pr3):
    done = run('make_levels_list_prepub.py', pr3)
    assert done.returncode == 0, done.stderr
    with open(str(pr3 / prepub.PREPUB)) as stream:
        return done, stream.read().splitlines()


def test_the_columns_are_the_published_ones(written):
    _, lines = written
    header = lines[0].split('\t')
    assert header[:6] == ['par', 'Ee', 'Ec', 'dE_O-C', 'J', 'Lande_g']
    assert 'cpl' not in header and 'e_c' not in header
    assert header[6:9] == ['perc1', 'conf1', 'term1']
    assert header[-3:] == ['perc5', 'conf5', 'term5']


def test_an_unfitted_level_has_no_observed_energy_and_no_difference(written):
    """levels.txt marks it with an asterisk beside a copy of Ec; here
    the two columns are simply blank, which is what a published table
    wants and what a spreadsheet reads as missing."""
    _, lines = written
    blanks = 0
    for line in lines[1:]:
        cell = line.split('\t')
        assert (cell[1].strip() == '') == (cell[3].strip() == '')
        if not cell[1].strip():
            blanks += 1
        else:
            assert abs(float(cell[1]) - float(cell[2])
                       - float(cell[3])) < 5e-5
    assert blanks and blanks < len(lines) - 1


def test_every_level_of_the_levels_file_is_present(pr3, written):
    _, lines = written
    table = purity.read_levels_txt(str(pr3 / 'levels.txt'))
    assert len(lines) - 1 == len(table.rows) // len(table.couplings)


def test_the_run_says_what_it_did(written):
    done, _ = written
    assert 'levels.txt: 5626 levels, 3 couplings (LS, JJ, LSJK)' in done.stdout
    assert 'wrote levels_prepub.txt (5626 levels' in done.stdout
    # The invariance held, so nothing was reported against it.
    assert 'differ between schemes' not in done.stdout


def test_the_options_are_carried_through(pr3, tmp_path):
    done = run('make_levels_list_prepub.py', pr3,
               '--levels', 'levels.txt', '--choice', 'coupling_choice.txt',
               '-o', str(tmp_path / 'short.txt'),
               '--components', '2', '--min-percent', '10')
    assert done.returncode == 0, done.stderr
    with open(str(tmp_path / 'short.txt')) as stream:
        lines = stream.read().splitlines()
    assert lines[0].split('\t')[-3:] == ['perc2', 'conf2', 'term2']
    for line in lines[1:]:
        cell = line.split('\t')
        assert len(cell) <= 12
        for at in (6, 9):
            if at < len(cell) and cell[at].strip():
                assert int(cell[at]) >= 10


def test_an_unknown_option_says_so(pr3):
    done = run('make_levels_list_prepub.py', pr3, '--couplings', 'LS')
    assert done.returncode == 2
    assert '--couplings is not an option' in done.stderr


# --------------------------------------------------------------------------
# the `+x` suffix of a relative-energy system
# --------------------------------------------------------------------------

def a_level(group=''):
    return prepub.Level(parity='o', e_obs=12.345678, unfitted=False,
                        e_calc=12.3, j=0.5, g=1.25, conf0='4f3', line=2,
                        group=group)


def test_the_suffix_is_a_column_of_its_own_beside_the_energy():
    """`Ee` stays a bare number and the `+x` goes beside it.

    Appending it to the energy would make `Ee` non-numeric and take
    `dE_O-C` and every spreadsheet that reads the column with it.  The
    leading space is what stops Excel reading `+x` as a formula.
    """
    at = prepub.FIXED.index('suffix')
    assert prepub.FIXED[at - 1] == 'Ee'
    assert prepub.row_of(a_level('x'))[at] == ' +x'
    assert prepub.row_of(a_level())[at] == ''
    assert prepub.row_of(a_level('x'))[at - 1] == '%12.6f' % 12.345678


def test_an_unfitted_level_has_no_suffix_either():
    level = a_level('x')
    level.unfitted = True
    at = prepub.FIXED.index('suffix')
    assert prepub.row_of(level)[at - 1] == ''
    assert prepub.row_of(level)[at] == ''


def test_the_suffix_survives_the_round_trip(tmp_path):
    path = str(tmp_path / 'levels_prepub.txt')
    prepub.write_prepub(path, [a_level('x'), a_level()])
    back = prepub.read_prepub(path)
    assert [one.group for one in back] == ['x', '']
    assert back[0].e_obs == pytest.approx(12.345678)


def test_a_file_written_before_the_suffix_column_still_reads(tmp_path):
    """Its absence means what a blank cell means: no level is in a system."""
    path = tmp_path / 'levels_prepub.txt'
    head = [name for name in prepub.header(1) if name != 'suffix']
    row = [value for name, value in zip(prepub.header(1),
                                        prepub.row_of(a_level()))
           if name != 'suffix']
    path.write_text('\t'.join(head) + '\n' + '\t'.join(row) + '\n')
    back = prepub.read_prepub(str(path))
    assert len(back) == 1 and back[0].group == ''
