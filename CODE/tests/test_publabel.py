"""Designations unique within each (parity, J).

Phase 14 has no Perl ancestor either.  What is pinned here is the
progression -- leading component, then further down the composition,
then letters, then the energy -- and the one property the whole file
exists for: no two published levels of one parity and J carry the same
name.  The search is exercised on constructed groups, where the answer
can be worked out by hand, and then on all 1036 published levels of
Pr III.
"""

import os
import subprocess
import sys

import pytest

from conftest import case
from toolbox import prepub, publabel
from toolbox.publabel import PubError
from toolbox.purity import Component

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(script, folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, script)] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


@pytest.fixture(scope='module')
def pr3(tmp_path_factory):
    """Pr III through phases 6, 12 and 13: a levels_prepub.txt.

    The pipeline as it is meant to be run -- ten components down to
    0.5% out of `printout`, three schemes with the third from `JK/`,
    then `prepub` doing the trimming.
    """
    folder = case('Pr3')
    out = tmp_path_factory.mktemp('publabel')
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
    done = run('make_levels_list_prepub.py', out, '--levels', 'levels.txt')
    assert done.returncode == 0, done.stderr
    return out


@pytest.fixture(scope='module')
def named(pr3):
    """The published Pr III levels, named."""
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    published, counts = publabel.cut(levels)
    names, lettered, nameless = publabel.name_all(
        publabel.choose(published))
    return levels, published, counts, names, lettered, nameless


# --------------------------------------------------------------------------
# levels_prepub.txt reads back
# --------------------------------------------------------------------------

def test_prepub_reads_back_what_it_wrote(pr3, tmp_path):
    """Phase 13's file is phase 14's input, so it has to survive both."""
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    assert len(levels) == 5626
    again = str(tmp_path / 'again.txt')
    prepub.write_prepub(again, levels)
    with open(str(pr3 / 'levels_prepub.txt')) as first:
        with open(again) as second:
            assert [line.rstrip('\t') for line in first.read().splitlines()] \
                == [line.rstrip('\t') for line in second.read().splitlines()]


def test_an_unfitted_level_comes_back_unfitted(pr3):
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    unfitted = [one for one in levels if one.unfitted]
    assert unfitted
    # Blank Ee means Ee holds a copy of Ec, as levels.txt's asterisk does.
    assert all(one.e_obs == one.e_calc for one in unfitted)
    assert all(not one.unfitted for one in levels if not one.unfitted)


# --------------------------------------------------------------------------
# how many levels are published
# --------------------------------------------------------------------------

def test_the_cut_is_the_highest_observed_level_of_each_set(named):
    levels, published, counts, _, _, _ = named
    assert counts == {'o': 481, 'e': 555}
    assert len(published) == 1036
    for parity, count in counts.items():
        mine = [one for one in levels if one.parity == parity]
        assert not mine[count - 1].unfitted          # the highest observed
        assert all(one.unfitted for one in mine[count:])


def test_nlev_overrides_the_cut(pr3):
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    published, counts = publabel.cut(levels, [100, 120])
    assert counts == {'o': 100, 'e': 120}
    assert len(published) == 220


def test_a_set_with_no_observed_level_makes_nlev_mandatory():
    levels = [_level('e', 10.0, 0.5, [(90, 'a', 'X')], unfitted=True),
              _level('e', 20.0, 1.5, [(90, 'a', 'Y')], unfitted=True)]
    with pytest.raises(PubError) as problem:
        publabel.cut(levels)
    assert 'the e set has no observed level' in str(problem.value)
    assert '--nlev' in str(problem.value)


def test_nlev_wants_one_count_per_set(pr3):
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    with pytest.raises(PubError) as problem:
        publabel.cut(levels, [100])
    assert 'one count per configuration set' in str(problem.value)


# --------------------------------------------------------------------------
# the search
# --------------------------------------------------------------------------

def _level(parity, e_calc, j, components, unfitted=False, e_obs=None):
    return prepub.Level(
        parity=parity, e_obs=e_calc if e_obs is None else e_obs,
        unfitted=unfitted, e_calc=e_calc, j=j, g=None, conf0='', line=0,
        components=[Component(percent=p, configuration=c, term=t)
                    for p, c, t in components])


def test_a_pure_level_keeps_its_leading_component(named):
    *_, names, _, _ = named
    pure = [one for one in names
            if one.level.components and one.level.components[0].percent > 90]
    assert pure
    assert all(one.leading for one in pure)


def test_the_swap_greedy_cannot_solve():
    """Both levels lead on X; the one that loses least gives it up.

    Greedy in energy order would hand X to the first level and leave the
    second with Y at 30%.  The search takes the total: giving X to the
    second costs it 45 - 44 = 1 and gains the first 40 - 30 = 10.
    """
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X'), (40, 'c', 'Y')]),
             _level('e', 2.0, 2.5, [(45, 'c', 'X'), (30, 'c', 'Z')])]
    answer = publabel.choose(group)
    assert [one[1][1] for one in answer] == [('c', 'Y'), ('c', 'X')]


def test_the_leading_component_is_preferred_over_a_larger_total():
    """A level does not give up its lead to gain a point elsewhere."""
    group = [_level('e', 1.0, 2.5, [(40, 'c', 'X'), (39, 'c', 'Y')]),
             _level('e', 2.0, 2.5, [(38, 'c', 'Y')])]
    answer = publabel.choose(group)
    # Taking Y for the first would total 39 + nothing; keeping the lead
    # names both, and naming a level at all outweighs a better name.
    assert [one[1][1] for one in answer] == [('c', 'X'), ('c', 'Y')]


def test_a_component_below_min_name_percent_is_no_candidate():
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.0, 2.5, [(45, 'c', 'X'), (19, 'c', 'Y')])]
    answer = publabel.choose(group, min_percent=20.0)
    # Only one X to go round, and 45 keeps it; 19% is not offered to the
    # other, so it is left for the letters.
    assert [one[1] for one in answer] == [None, (245, ('c', 'X'), 0)]
    # Offered, and taking it names both levels, which comes first.
    lower = publabel.choose(group, min_percent=15.0)
    assert [one[1][1] for one in lower] == [('c', 'X'), ('c', 'Y')]


def test_levels_of_different_j_do_not_compete():
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.0, 3.5, [(45, 'c', 'X')])]
    answer = publabel.choose(group)
    assert [one[1][1] for one in answer] == [('c', 'X'), ('c', 'X')]


# --------------------------------------------------------------------------
# the letters
# --------------------------------------------------------------------------

def test_a_group_with_no_assignment_falls_to_letters():
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.0, 2.5, [(45, 'c', 'X')])]
    names, lettered, nameless = publabel.name_all(publabel.choose(group))
    assert not nameless
    assert [one.label for one in names] == [('c', 'a X'), ('c', 'b X')]
    assert len(lettered) == 2


def test_both_sides_of_a_contest_take_a_letter(named):
    """`4G*` beside `a 4G*` in one (parity, J) would read as two things."""
    *_, names, lettered, _ = named
    plain = set()
    marked = set()
    for one in names:
        key = (one.level.parity, one.level.j, one.conf, one.term)
        (marked if one.letter else plain).add(key)
    assert marked
    assert not (marked & plain)


def test_the_letters_run_across_j_and_not_within_one_j():
    """One family, two multiplets, and the letters are its own.

    Within a (parity, J) group the count would restart, so the second
    multiplet would be `b` where it is contested and `a` where it is
    not.  Counted over the family it is `b` throughout.
    """
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.0, 2.5, [(45, 'c', 'X')]),
             _level('e', 3.0, 3.5, [(44, 'c', 'X')]),
             _level('e', 4.0, 3.5, [(45, 'c', 'X')])]
    names, _, _ = publabel.name_all(publabel.choose(group))
    assert [one.label[1] for one in names] == ['a X', 'b X', 'c X', 'd X']
    assert publabel._letter(0) == 'a'
    assert publabel._letter(25) == 'z'
    assert publabel._letter(26) == 'aa'


def test_a_level_with_no_component_at_all_is_named_by_its_energy():
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 12.3456, 2.5, [], e_obs=12.3456)]
    names, _, nameless = publabel.name_all(publabel.choose(group))
    assert len(nameless) == 1
    assert names[1].label == ('12345', '')
    assert names[1].percent is None
    assert names[1].energy


# --------------------------------------------------------------------------
# the global fallback
# --------------------------------------------------------------------------

def test_too_many_letters_makes_every_name_an_energy():
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.0, 2.5, [(45, 'c', 'X')]),
             _level('e', 3.0, 3.5, [(90, 'c', 'Y')])]
    names, lettered, _ = publabel.name_all(publabel.choose(group))
    why, share = publabel.fallback(names, lettered)
    assert why == 'letters'
    assert round(share) == 67
    publabel.by_energy(names)
    assert [one.label for one in names] == [('1000', ''), ('2000', ''),
                                            ('3000', '')]
    assert all(one.percent is None for one in names)


def test_too_much_mixing_makes_every_name_an_energy():
    group = [_level('e', 1.0, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.0, 3.5, [(45, 'c', 'Y')]),
             _level('e', 3.0, 4.5, [(90, 'c', 'Z')])]
    names, lettered, _ = publabel.name_all(publabel.choose(group))
    assert not lettered
    why, share = publabel.fallback(names, lettered)
    assert why == 'mixing'
    assert round(share) == 67


def test_colliding_integer_energies_become_a_sequential_index():
    """Two levels of one parity inside the same whole cm-1."""
    group = [_level('e', 1.20001, 2.5, [(44, 'c', 'X')]),
             _level('e', 1.20009, 3.5, [(45, 'c', 'Y')]),
             _level('o', 1.9, 3.5, [(45, 'c', 'Z')])]
    names, _, _ = publabel.name_all(publabel.choose(group))
    assert publabel.by_energy(names) is True
    assert [one.label[0] for one in names] == ['1', '2', '1']


def test_distinct_integer_energies_stay_energies():
    group = [_level('e', 1.2, 2.5, [(44, 'c', 'X')]),
             _level('e', 2.8, 3.5, [(45, 'c', 'Y')])]
    names, _, _ = publabel.name_all(publabel.choose(group))
    assert publabel.by_energy(names) is False
    assert [one.label[0] for one in names] == ['1200', '2800']


# --------------------------------------------------------------------------
# the file
# --------------------------------------------------------------------------

def test_every_published_name_is_unique_in_its_parity_and_j(named):
    """The one property the whole file exists for."""
    *_, names, _, _ = named
    assert publabel.duplicates(names) == []


def test_the_columns(named):
    *_, names, _, _ = named
    assert publabel.header(3) == [
        'par', 'conf', 'term', 'perc',
        'perc1', 'conf1', 'term1', 'perc2', 'conf2', 'term2',
        'perc3', 'conf3', 'term3',
        'Ee', 'suffix', 'Ec', 'dE_O-C', 'J', 'Lande_g']
    assert publabel.header(1) == ['par', 'conf', 'term', 'perc',
                                  'perc1', 'conf1', 'term1',
                                  'Ee', 'suffix', 'Ec', 'dE_O-C', 'J',
                                  'Lande_g']
    assert all(len(publabel.row_of(one, 3)) == len(publabel.header(3))
               for one in names)


def test_group_one_is_blank_where_the_name_came_from_it(named):
    """The ASD rule: a level says once where its name came from."""
    *_, names, _, _ = named
    from_lead = [one for one in names if one.leading]
    from_below = [one for one in names if one.index > 0]
    assert from_lead and from_below
    for one in from_lead:
        assert publabel.row_of(one, 3)[4:7] == ['', '', '']
    for one in from_below:
        row = publabel.row_of(one, 3)
        assert row[5] == one.level.components[0].configuration
        assert row[6] == one.level.components[0].term


def test_ee_and_de_are_blank_together(pr3, tmp_path):
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    published, _ = publabel.cut(levels, [200, 200])
    names, _, _ = publabel.name_all(publabel.choose(published))
    out = str(tmp_path / 'levels_pub.txt')
    publabel.write_pub(out, names)
    with open(out) as stream:
        lines = stream.read().splitlines()
    at = {name: i for i, name in enumerate(lines[0].split('\t'))}
    blank = 0
    for line in lines[1:]:
        cell = line.split('\t')
        assert bool(cell[at['Ee']].strip()) == bool(cell[at['dE_O-C']].strip())
        blank += not cell[at['Ee']].strip()
    assert blank                    # there are unfitted levels among these


def test_the_energies_are_written_in_wavenumbers(pr3, tmp_path):
    """The published file is in cm-1, the pipeline above it in kK."""
    levels = prepub.read_prepub(str(pr3 / 'levels_prepub.txt'))
    published, _ = publabel.cut(levels, [200, 200])
    names, _, _ = publabel.name_all(publabel.choose(published))
    out = str(tmp_path / 'levels_pub.txt')
    publabel.write_pub(out, names)
    with open(out) as stream:
        lines = stream.read().splitlines()
    at = {name: i for i, name in enumerate(lines[0].split('	'))}
    seen = 0
    for line, one in zip(lines[1:], names):
        cell = line.split('	')
        assert float(cell[at['Ec']]) == round(one.level.e_calc * 1000, 3)
        if cell[at['Ee']].strip():
            seen += 1
            assert float(cell[at['Ee']]) == round(one.level.e_obs * 1000, 3)
            assert float(cell[at['dE_O-C']]) == round(
                (one.level.e_obs - one.level.e_calc) * 1000, 3)
    assert seen
    # Sugar's Pr III level at 82492.4 cm-1, not at 82.4924.
    assert any(abs(float(line.split('	')[at['Ee']] or 0) - 82492.400) < 5e-4
               for line in lines[1:])


def test_a_name_from_outside_the_printed_groups_is_found(named):
    *_, names, _, _ = named
    assert publabel.outside(names, 5) == []
    beyond = publabel.outside(names, 1)
    assert beyond
    assert all(one.index >= 1 for one in beyond)


# --------------------------------------------------------------------------
# the script
# --------------------------------------------------------------------------

def test_the_script_writes_the_file(pr3):
    done = run('make_levels_list_pub.py', pr3)
    assert done.returncode == 0, done.stderr
    assert 'publishing 1036' in done.stdout
    assert 'take a letter prefix' in done.stdout
    with open(str(pr3 / 'levels_pub.txt')) as stream:
        lines = stream.read().splitlines()
    assert len(lines) == 1037
    assert lines[0].split('\t') == publabel.header(3)


def test_the_script_takes_its_options(pr3):
    done = run('make_levels_list_pub.py', pr3, '--nlev', '100,120',
               '--n_comp', '1', '-o', 'small.txt')
    assert done.returncode == 0, done.stderr
    with open(str(pr3 / 'small.txt')) as stream:
        lines = stream.read().splitlines()
    assert len(lines) == 221
    assert lines[0].split('\t') == publabel.header(1)


def test_n_comp_is_between_one_and_five(pr3):
    done = run('make_levels_list_pub.py', pr3, '--n_comp', '6')
    assert done.returncode == 2
    assert '--n_comp is between 1 and 5' in done.stderr


def test_an_unknown_option_stops(pr3):
    done = run('make_levels_list_pub.py', pr3, '--nope', '1')
    assert done.returncode == 2
    assert '--nope is not an option' in done.stderr


# --------------------------------------------------------------------------
# the other two reference cases
# --------------------------------------------------------------------------

def _pipeline(folder, out, *couplings):
    """A whole run through phases 6, 12, 13 and 14."""
    with open(os.path.join(str(folder), 'IN36'), 'rb') as source:
        data = source.read()
    with open(str(out / 'IN36'), 'wb') as target:
        target.write(data)
    done = run('printout.py', folder,
               str(out / 'levels.txt'), str(out / 'params.txt'),
               *(list(couplings) + ['--dominant-config', '--components', '10',
                                    '--min-percent', '0.5']))
    assert done.returncode == 0, done.stderr
    for script in ('config_purity.py', 'make_levels_list_prepub.py'):
        done = run(script, out, '--levels', 'levels.txt')
        assert done.returncode == 0, done.stderr
    return out


def test_w70_publishes_with_a_true_jj_block(tmp_path):
    """No LSF, so `--nlev` is mandatory and the abort names the set."""
    out = _pipeline(case('W70'), tmp_path, '--couplings', 'LS,JJ,TRUEjj')
    done = run('make_levels_list_pub.py', out)
    assert done.returncode == 1
    assert 'the o set has no observed level' in done.stderr
    assert '--nlev' in done.stderr

    done = run('make_levels_list_pub.py', out, '--nlev', '80,89')
    assert done.returncode == 0, done.stderr
    levels = prepub.read_prepub(str(out / 'levels_prepub.txt'))
    published, _ = publabel.cut(levels, [80, 89])
    names, _, _ = publabel.name_all(publabel.choose(published))
    assert publabel.duplicates(names) == []
    # Near-pure jj: the relativistic labels do reach the published file.
    assert any('2p-' in one.conf or '2p+' in one.conf for one in names)


def test_u90_publishes_with_a_reordered_block(tmp_path):
    out = _pipeline(case('U90'), tmp_path,
                    '--couplings', 'LS,JJ,LSr', '--from', 'LSr=LS@reordered')
    done = run('make_levels_list_pub.py', out)
    assert done.returncode == 0, done.stderr
    levels = prepub.read_prepub(str(out / 'levels_prepub.txt'))
    published, counts = publabel.cut(levels)
    names, _, _ = publabel.name_all(publabel.choose(published))
    assert publabel.duplicates(names) == []
    assert len(names) == sum(counts.values())


def test_a_missing_prepub_file_says_what_writes_it(tmp_path):
    done = run('make_levels_list_pub.py', tmp_path)
    assert done.returncode == 1
    assert 'make_levels_list_prepub.py' in done.stderr
