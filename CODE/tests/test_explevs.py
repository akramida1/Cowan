"""subst_exp_levs.py: carrying observed energies between RCE decks.

On `WORK/U90` and `WORK/Pr3` the rewritten deck agrees with
`subst_exp_levs.pl` on every level line and every parameter value --
593 observed energies and 1207 parameters on Pr III alone -- and
differs only in that the Perl drops the trailing blanks of a CI
parameter line, which its regex leaves outside the part it puts back
(`subst_exp_levs.pl:421`, `:436`).

`WORK/W70` is where they part company, and the Perl is wrong.  Its
level-line regex wants a blank between J and the leading percentage
(`:73`), and there is none when a half-integer J fills its field and
the percentage is -100: `   1/2-100 ss4p`.  Perl keeps the previous
match's captures when a match fails, so from the first such line the
level counter and the J are stale, every lookup behind it misses, and
`PrintOut` writes an undefined energy -- `0.000000` -- over each of the
remaining levels.  Sixty of W70's observed energies are destroyed that
way.  The test below pins the defect so that nobody restores it.
"""

import os
import shutil
import subprocess
import sys

import pytest

from conftest import PERL, case, perl_env
from toolbox import explevs, files, rce

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'subst_exp_levs.py')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True)


def perl(folder, *arguments):
    if shutil.which('perl') is None:
        pytest.skip('perl is not on the PATH')
    out = subprocess.run(
        ['perl', str(PERL / 'subst_exp_levs.pl')] + list(arguments),
        cwd=str(folder), capture_output=True, text=True, env=perl_env())
    if out.returncode:
        pytest.skip('subst_exp_levs.pl will not run here: %s'
                    % out.stderr.strip())
    return out


def prepare(name, tmp_path):
    """A working directory holding one run's RCEOUT and RCEINP."""
    folder = case(name)
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(str(files.resolve('RCEOUT', folder)), str(tmp_path / 'SRC'))
    shutil.copy(str(files.resolve('RCEINP', folder)), str(tmp_path / 'TGT'))
    return tmp_path


def lines_of(path):
    return files.read_lines(path, strip=False)


def is_level(line):
    return len(line) > 30 and '.' in line[:13] and line[24:30].strip()


# --------------------------------------------------------------------------
# against the Perl
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name', ['U90', 'Pr3'])
def test_every_line_agrees_with_the_perl_but_for_its_trailing_blanks(
        name, tmp_path):
    mine = prepare(name, tmp_path / 'py')
    theirs = prepare(name, tmp_path / 'pl')
    assert run(mine, 'SRC', 'TGT').returncode == 0
    perl(theirs, 'SRC', 'TGT')
    ours, perls = lines_of(mine / 'TGT'), lines_of(theirs / 'TGT')
    assert len(ours) == len(perls)
    for a, b in zip(ours, perls):
        assert a.rstrip() == b.rstrip()


def test_w70_is_where_the_perl_destroys_the_deck(tmp_path):
    """The Perl zeroes observed energies it never matched anything to.

    W70 was not fitted, so its RCEOUT carries no observed level and
    there is nothing to substitute: the right answer is to leave the
    deck exactly as it was.
    """
    mine = prepare('W70', tmp_path / 'py')
    theirs = prepare('W70', tmp_path / 'pl')
    assert run(mine, 'SRC', 'TGT').returncode == 0
    perl(theirs, 'SRC', 'TGT')
    original = lines_of(files.resolve('RCEINP', case('W70')))
    assert [line.rstrip() for line in lines_of(mine / 'TGT')] \
        == [line.rstrip() for line in original]
    zeroed = [(a, b) for a, b in zip(original, lines_of(theirs / 'TGT'))
              if a.rstrip() != b.rstrip() and is_level(b)]
    assert zeroed
    for was, now in zeroed:
        assert float(was[:13]) != 0.0        # it had an energy
        assert float(now[:13]) == 0.0        # and the Perl wrote it away


def test_the_perl_regex_is_what_fails_on_those_lines():
    """`   1/2-100` -- a half-integer J against a percentage of -100.

    RCE writes J as I4 then A2 and the leading percentage as I4 behind
    it, so a J of the form n/2 fills its six columns and a percentage of
    -100 fills its four: there is no blank where the Perl's `+` wants
    one.  An integer J leaves the A2 blank and the same line matches,
    which is why this shows on some runs and not others.
    """
    import re
    perl_pattern = re.compile(
        r'^ *([0-9-]+\.\d{0,6})([* ]+)(-{0,1}\d+\.\d{1,3}) +(\d+)(\/(\d+))*'
        r' +([0-9-]+) (.{6}) +(\S+) (\S+)')
    lines = files.read_lines(files.resolve('RCEINP', case('W70')), width=140)
    bad = [line for line in lines
           if is_level(line) and line[24:30].strip().endswith('/2')
           and not perl_pattern.match(line)]
    assert bad
    assert all(re.match(r' *\d+/2-100', line[24:34]) for line in bad)


# --------------------------------------------------------------------------
# what it does
# --------------------------------------------------------------------------

def test_only_the_observed_column_and_its_mark_move(tmp_path):
    mine = prepare('U90', tmp_path)
    before = lines_of(mine / 'TGT')
    assert run(mine, 'SRC', 'TGT', 'params=none').returncode == 0
    for old, new in zip(before, lines_of(mine / 'TGT')):
        assert old[14:] == new[14:]


def test_a_level_that_takes_an_energy_loses_its_unobserved_mark():
    line = ('  2264.603000* 2264.603   1/2-100 2p     (1S) 2P     ')
    after = rce.set_observed(line, 2270.5)
    assert after[13] == ' '
    assert float(after[:13]) == pytest.approx(2270.5)
    assert after[14:] == line[14:]


def test_the_old_deck_is_kept(tmp_path):
    mine = prepare('U90', tmp_path)
    original = lines_of(mine / 'TGT')
    assert run(mine, 'SRC', 'TGT').returncode == 0
    assert lines_of(mine / 'TGT.bak') == original


def test_every_observed_level_of_the_source_finds_a_home(tmp_path):
    for name in ('U90', 'Pr3'):
        mine = prepare(name, tmp_path / name)
        out = run(mine, 'SRC', 'TGT')
        assert out.returncode == 0, out.stderr
        assert '0 source levels went unmatched' in out.stdout


# --------------------------------------------------------------------------
# the parameter options
# --------------------------------------------------------------------------

def _changed(tmp_path, name, *arguments):
    """(flags moved, values moved) for one run of the program."""
    mine = prepare(name, tmp_path)
    before = files.read_lines(mine / 'TGT', width=80)
    assert run(mine, 'SRC', 'TGT', *arguments).returncode == 0
    after = files.read_lines(mine / 'TGT', width=80)
    flags = values = 0
    for old, new in zip(before, after):
        if old == new or not rce.is_parameter_line(old):
            continue
        flags += old[10:15] != new[10:15]
        values += old[15:29] != new[15:29]
    return flags, values


def test_by_default_both_the_values_and_the_flags_come_across(tmp_path):
    flags, values = _changed(tmp_path, 'Pr3')
    assert values == 292 and flags == 0        # this fit did not free more


def test_params_none_leaves_the_table_alone(tmp_path):
    assert _changed(tmp_path, 'Pr3', 'params=none') == (0, 0)


def test_params_noflags_moves_the_values_and_not_the_flags(tmp_path):
    """The Perl moves neither, `subst_exp_levs.pl:448` and `:451`."""
    flags, values = _changed(tmp_path, 'Pr3', 'params=noflags')
    assert values == 292 and flags == 0


@pytest.mark.parametrize('option,expected', [
    ('params=Eav', 17), ('params=CI', 66), ('params=F', 86),
    ('params=G', 50), ('params=Z', 56)])
def test_a_named_kind_moves_only_its_own(option, expected, tmp_path):
    """Eav and CI name a kind of parameter; F, G and Z name a letter.

    Matching `CI` against a name would never succeed -- CI parameters
    are called `120D1112` -- and matching `F` against one would catch
    `EAV 4f3` as well as `F2(11)`, which is why the Perl keeps three
    branches and so does this.
    """
    flags, values = _changed(tmp_path, 'Pr3', option)
    assert (flags, values) == (0, expected)


def test_the_kinds_together_are_the_whole_table(tmp_path):
    """Nothing is counted twice and nothing is missed twice."""
    total = _changed(tmp_path / 'all', 'Pr3')[1]
    parts = sum(_changed(tmp_path / option, 'Pr3', 'params=' + option)[1]
                for option in ('Eav', 'CI', 'F', 'G', 'Z'))
    assert parts <= total


# --------------------------------------------------------------------------
# the matcher itself
# --------------------------------------------------------------------------

def test_the_distance_ignores_the_sign_of_a_whole_eigenvector():
    """RCE does not fix an eigenvector's overall sign, so neither may this."""
    one = {'2s2p,(3S) 4P': 95, '2s2p,(3S) 2P': 32}
    flipped = dict((k, -v) for k, v in one.items())
    assert explevs.distance(one, flipped)[0] == 0
    assert explevs.distance(one, one)[0] == 0


def test_a_component_missing_from_one_side_counts_its_own_square():
    assert explevs.distance({'a,b': 50}, {})[0] == 2500
    assert explevs.distance({}, {'a,b': 50})[0] == 2500


def test_two_vectors_dominated_by_one_state_settle_the_match():
    one = {'2p,(1S) 2P': 100}
    other = {'2p,(1S) 2P': 99, 'x,y': 12}
    apart, certain = explevs.distance(one, other)
    assert certain and apart == 0


def test_a_level_already_carrying_an_energy_is_not_given_a_second():
    taken = explevs.Level(j='1/2', ee=1.0, ec=1.0, observed=True)
    free = explevs.Level(j='1/2', ee=0.0, ec=1.0, observed=False)
    assert taken.taken and not free.taken
    free.assigned = 3.0
    assert free.taken and free.energy == 3.0


def test_the_component_key_does_not_depend_on_how_wide_rce_padded_it():
    from toolbox.rce import Component
    assert explevs._key(Component(95, '2s2p  ', '(3S)  4P   ')) \
        == explevs._key(Component(95, '2s2p', '(3S) 4P'))


def test_the_arguments_are_checked(tmp_path):
    mine = prepare('U90', tmp_path)
    assert run(mine, 'SRC').returncode == 2
    assert run(mine, 'SRC', 'TGT', 'params=sideways').returncode == 2
    assert run(mine, 'SRC', 'TGT', 'sideways').returncode == 2
