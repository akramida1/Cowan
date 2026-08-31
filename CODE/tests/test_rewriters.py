"""The deck rewriters: add_eav, s11, subst_params_ing11, scale_param.

The acceptance criterion for all four is the plan's: run the Perl and
the Python on the reference decks and account for every difference.
Three of them come out **byte-identical** to `s11.pl`, `add_eav.pl` and
`subst_params_ing11.pl`; the exceptions are listed one by one below and
written up in `FOR/LSJJ/README.txt`.

What holds the rest together is one property, tested first: a card that
nothing has edited comes back out as the bytes it came in as, and so
does a card that is re-encoded from values nothing has changed.  A
rewriter therefore changes the cards it means to change and no others,
whatever else is in the deck.
"""

import os
import shutil
import subprocess
import sys

import pytest

import s11 as s11_tool

from conftest import PERL, case, perl_env
from toolbox import files, ing11, rce

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Every deck of the reference runs, with the name it is spelled under.
DECKS = [('W70', 'ING11'), ('U90', 'ing11'), ('U90/reordered', 'ing11'),
         ('Pr3', 'ing11'), ('Pr3/JK', 'ing11')]


def run(script, folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, script)] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


def perl(script, folder, *arguments):
    """Run one of the Perl utilities as the reference, or skip.

    `add_eav.pl:7` opens its deck as `ing11.` -- with the trailing dot
    Windows drops and a POSIX filesystem does not -- so under a Cygwin
    or MSYS Perl it cannot open the file at all and there is no
    reference run to compare against.  That is reported as a skip with
    the Perl's own message rather than as a failure of the port.
    """
    if shutil.which('perl') is None:
        pytest.skip('perl is not on the PATH')
    out = subprocess.run(
        ['perl', str(PERL / script)] + list(arguments),
        cwd=str(folder), capture_output=True, text=True, env=perl_env())
    if out.returncode:
        pytest.skip('%s will not run under this perl: %s'
                    % (script, out.stderr.strip() or out.stdout.strip()))
    return out


def copy(folder, tmp_path, *names):
    """Put the named files of a reference run into a temporary directory.

    `conftest.case` gates on OUTG11, which says RCG has run; a folder
    where RCE has not run as well is half a reference run, and a test
    that wants RCEOUT out of it skips rather than fails.
    """
    tmp_path.mkdir(parents=True, exist_ok=True)
    for name in names:
        found = files.resolve(name, folder, required=False)
        if found is None:
            pytest.skip('%s has no %s (see its README.txt)' % (folder, name))
        shutil.copy(str(found), str(tmp_path / name))
    return tmp_path


def lines_of(path):
    return files.read_lines(path, strip=False)


def perl_deck(folder, deck):
    """The deck the Perl wrote, wherever it managed to leave it.

    `add_eav.pl` and `s11.pl` build the new deck in a temporary file and
    then move it into place by shelling out to `copy` and `ren`, which
    are cmd.exe built-ins; run under a Perl whose backticks go to a
    POSIX shell -- a Cygwin or MSYS one, which is what a `pytest` from a
    Git Bash prompt finds first -- the move fails silently and the deck
    is left as it was.  The temporary file is the output either way, so
    it is what gets compared, and the shell-out is not what is under
    test.  That both programs need one at all is the second reason the
    plan gives for replacing them.
    """
    for name in ('tmp.', 'tmp'):
        temporary = folder / name
        if temporary.is_file() and temporary.stat().st_size:
            return lines_of(temporary)
    return lines_of(folder / deck)


# --------------------------------------------------------------------------
# the property everything else rests on
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name,deck', DECKS)
def test_a_deck_nothing_edited_round_trips_byte_for_byte(name, deck):
    folder = case(name)
    path = files.resolve(deck, folder)
    assert ing11.write_ing11(ing11.read_ing11(path)) == lines_of(path)


@pytest.mark.parametrize('name,deck', DECKS)
def test_re_encoding_every_card_changes_nothing_either(name, deck):
    """A field whose value has not moved keeps the characters it had.

    RCN2 and the Perl have each written a zero parameter in a form the
    other does not produce -- ` 0.00E+00` against `  0.00e+0` -- and
    both read as zero.  Re-encoding a card must not choose between them
    for a value nobody has touched.
    """
    folder = case(name)
    path = files.resolve(deck, folder)
    parsed = ing11.read_ing11(path)
    for k in (0, 1):
        for card in parsed.parameters[k] + parsed.ci[k]:
            card.edited = True
    assert ing11.write_ing11(parsed) == lines_of(path)


def test_a_new_value_is_encoded_as_rcn2_encodes_it():
    """`fortfmt.write_param` is RCN2's own rule, `FOR/RCN2K.F:1941-1963`."""
    from toolbox.fortfmt import write_param
    assert write_param(1.405683, first=True) == '    14057'
    assert write_param(773862.017396, first=True) == ' 773862.0'
    assert write_param(0.0) == ' 0.00E+00'
    assert len(write_param(1e-10)) == 9


def test_a_ci_value_gives_up_a_decimal_rather_than_its_field():
    """The defect the author reported, `s11.pl:228`.

    `%9.4f` of -1369.505304 is ten characters, which pushes every field
    behind it one column right and corrupts the card.
    """
    assert ing11.format_value(-1369.505304, 1, ci=True) == '-1369.505'
    assert len(ing11.format_value(-1369.505304, 1, ci=True)) == 9


# --------------------------------------------------------------------------
# add_eav
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name,deck', DECKS)
def test_add_eav_agrees_with_the_perl_byte_for_byte(name, deck, tmp_path):
    folder = case(name)
    mine = copy(folder, tmp_path / 'py', deck)
    theirs = copy(folder, tmp_path / 'pl', deck)
    assert run('add_eav.py', mine, '10.0').returncode == 0
    perl('add_eav.pl', theirs, '10.0')
    assert lines_of(mine / deck) == perl_deck(theirs, deck)


def test_add_eav_moves_every_eav_and_nothing_else(tmp_path):
    folder = case('U90')
    mine = copy(folder, tmp_path, 'ing11')
    before = ing11.read_ing11(mine / 'ing11')
    assert run('add_eav.py', mine, '-2.5').returncode == 0
    after = ing11.read_ing11(mine / 'ing11')
    for k in (0, 1):
        for old, new in zip(before.parameters[k], after.parameters[k]):
            assert new.values[0] == pytest.approx(old.values[0] - 2.5,
                                                  abs=1e-4)
            assert new.values[1:] == old.values[1:]
        for old, new in zip(before.ci[k], after.ci[k]):
            assert new.values == old.values


def test_add_eav_keeps_the_old_deck(tmp_path):
    folder = case('U90')
    mine = copy(folder, tmp_path, 'ing11')
    assert run('add_eav.py', mine, '1.0').returncode == 0
    assert lines_of(mine / 'ING11.BAK') == lines_of(files.resolve('ing11',
                                                                  folder))


def test_add_eav_refuses_a_shift_that_is_not_a_number(tmp_path):
    mine = copy(case('U90'), tmp_path, 'ing11')
    assert run('add_eav.py', mine, 'sideways').returncode == 2


# --------------------------------------------------------------------------
# s11
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name,deck', DECKS)
def test_s11_agrees_with_the_perl_except_where_the_perl_overflows(
        name, deck, tmp_path):
    """W70, Pr3 and Pr3/JK are byte-identical.

    Both U90 decks differ on exactly one line, the CI card whose value
    passes 1000 in magnitude and which `s11.pl` writes ten characters
    wide.  That line is what the author reported and it is the only
    difference on any reference deck.
    """
    folder = case(name)
    mine = copy(folder, tmp_path / 'py', deck, 'RCEOUT')
    theirs = copy(folder, tmp_path / 'pl', deck, 'RCEOUT')
    assert run('s11.py', mine, '2').returncode == 0
    perl('s11.pl', theirs, '2')
    ours, perls = lines_of(mine / deck), perl_deck(theirs, deck)
    differ = [i for i, (a, b) in enumerate(zip(ours, perls)) if a != b]
    if name.startswith('U90'):
        assert len(differ) == 1
        assert perls[differ[0]].startswith('2s       -2p2')
        assert len(ours[differ[0]]) == len(perls[differ[0]]) - 1
    else:
        assert differ == []


@pytest.mark.parametrize('name', ['U90', 'U90/reordered'])
def test_s11_reproduces_the_authors_own_deck(name, tmp_path):
    """The stronger check: the deck the author keeps is what comes out.

    `WORK/U90/ing11` is the corrected deck -- the CI card the Perl
    corrupts, put right -- so reproducing it field for field says the
    encoding is right and not merely different from the Perl's.
    """
    folder = case(name)
    mine = copy(folder, tmp_path, 'ing11', 'RCEOUT')
    assert run('s11.py', mine, '2').returncode == 0
    theirs = ing11.read_ing11(files.resolve('ing11', folder))
    ours = ing11.read_ing11(mine / 'ing11')
    for k in (0, 1):
        for a, b in zip(theirs.parameters[k] + theirs.ci[k],
                        ours.parameters[k] + ours.ci[k]):
            assert a.name == b.name and a.groups == b.groups
            for x, y in zip(a.values, b.values):
                assert x == pytest.approx(y, abs=5e-4)


def test_s11_substitutes_one_parity_unless_asked_for_both(tmp_path):
    """Without the argument only the first parity set is substituted.

    What is checked is the run's own report and the text of the second
    set's cards, not their values: these decks already carry the fit
    being put back into them, so set 2 would hold the same numbers
    either way and only the report says which sets were walked.
    """
    folder = case('U90')
    one = copy(folder, tmp_path / 'one', 'ing11', 'RCEOUT')
    both = copy(folder, tmp_path / 'both', 'ing11', 'RCEOUT')
    first, second = run('s11.py', one), run('s11.py', both, '2')
    assert first.returncode == 0 and second.returncode == 0
    assert 'parity set 1' in first.stdout
    assert 'parity set 2' not in first.stdout
    assert 'parity set 1' in second.stdout and 'parity set 2' in second.stdout
    original = ing11.read_ing11(files.resolve('ing11', folder))
    after_one = ing11.read_ing11(one / 'ing11')
    assert [c.raw for c in after_one.parameters[1]] \
        == [c.raw for c in original.parameters[1]]


@pytest.mark.parametrize('name,parities', [
    ('U90', 2), ('Pr3', 2), ('W70', 2), ('AI', 1), ('Forbid', 1),
    ('AI/DIEL', 0),
])
def test_s11_auto_counts_the_parities_in_in36(work, name, parities):
    """`auto` reads off what the trailing 2 used to say.

    The decks that need both parity sets substituted are exactly the
    ones whose IN36 holds both parities.  AI/DIEL is the third answer:
    its shipped IN36 is the bare end card, `diel` writes the real one,
    and until then IN36 cannot say -- which is 0 here and the fallback
    in `main`.
    """
    folder = work / name
    if not folder.is_dir():
        pytest.skip('WORK/%s is not present' % name)
    assert s11_tool.parity_sets(folder) == parities


def test_s11_auto_substitutes_both_sets_when_in36_has_both(tmp_path):
    """`s11 auto` beside a two-parity IN36 is `s11 2`."""
    folder = case('U90')
    mine = copy(folder, tmp_path / 'auto', 'ing11', 'RCEOUT', 'IN36')
    theirs = copy(folder, tmp_path / 'two', 'ing11', 'RCEOUT')
    automatic, explicit = run('s11.py', mine, 'auto'), run('s11.py', theirs, '2')
    assert automatic.returncode == 0 and explicit.returncode == 0
    assert 'both parities' in automatic.stdout
    assert lines_of(mine / 'ing11') == lines_of(theirs / 'ing11')


def test_s11_auto_falls_back_when_there_is_no_in36(tmp_path):
    """A folder IN36 cannot answer for gets every set the deck has."""
    mine = copy(case('U90'), tmp_path, 'ing11', 'RCEOUT')
    out = run('s11.py', mine, 'auto')
    assert out.returncode == 0
    assert 'does not say how many parities' in out.stdout
    assert 'parity set 1' in out.stdout and 'parity set 2' in out.stdout


def test_s11_stops_when_rceout_does_not_match_the_deck(tmp_path):
    """A deck and an RCEOUT of two different runs must not be merged."""
    mine = copy(case('U90'), tmp_path, 'ing11')
    shutil.copy(str(files.resolve('RCEOUT', case('W70'))),
                str(tmp_path / 'RCEOUT'))
    out = run('s11.py', mine, '2')
    assert out.returncode == 1 and 'out of step' in out.stderr


# --------------------------------------------------------------------------
# subst_params_ing11
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name,source,target,kset', [
    ('U90', 'ing11', 'ing11.bak', '1'),
    ('U90', 'ing11', 'ing11.bak', '2'),
    ('Pr3', 'ing11', 'ing11.bak', '1'),
    ('Pr3', 'ing11', 'ing11.bak', '2'),
    ('W70', 'ING11', 'ING11', '1'),
])
def test_subst_params_agrees_with_the_perl_byte_for_byte(
        name, source, target, kset, tmp_path):
    folder = case(name)
    for who in ('py', 'pl'):
        where = tmp_path / who
        where.mkdir()
        shutil.copy(str(files.resolve(source, folder)), str(where / 'src'))
        shutil.copy(str(files.resolve(target, folder)), str(where / 'tgt'))
    assert run('subst_params_ing11.py', tmp_path / 'py',
               'src', 'tgt', kset).returncode == 0
    perl('subst_params_ing11.pl', tmp_path / 'pl', 'src', 'tgt', kset)
    assert lines_of(tmp_path / 'py' / 'tgt') \
        == lines_of(tmp_path / 'pl' / 'tgt')


def test_subst_params_carries_the_values_across_unre_encoded(tmp_path):
    folder = case('Pr3')
    shutil.copy(str(files.resolve('ing11', folder)), str(tmp_path / 'src'))
    shutil.copy(str(files.resolve('ing11.bak', folder)), str(tmp_path / 'tgt'))
    assert run('subst_params_ing11.py', tmp_path,
               'src', 'tgt', '1').returncode == 0
    source = ing11.read_ing11(tmp_path / 'src')
    target = ing11.read_ing11(tmp_path / 'tgt')
    by_name = dict((c.name, c) for c in target.parameters[0])
    for card in source.parameters[0]:
        assert by_name[card.name].raw == card.raw


def test_subst_params_stops_where_the_target_lacks_a_configuration(tmp_path):
    shutil.copy(str(files.resolve('ing11', case('Pr3'))),
                str(tmp_path / 'src'))
    shutil.copy(str(files.resolve('ing11', case('U90'))),
                str(tmp_path / 'tgt'))
    out = run('subst_params_ing11.py', tmp_path, 'src', 'tgt', '1')
    assert out.returncode == 1 and 'no configuration card' in out.stderr


# --------------------------------------------------------------------------
# scale_param
# --------------------------------------------------------------------------

def test_the_perl_scale_param_matches_nothing_any_more():
    """Why this one is not compared against the Perl.

    `scale_param.pl:58` needs the maximum-value field to read `0.0000`
    and RCE writes it as F14.6.  The regex matches no line of any
    reference deck, so the Perl copies its input unchanged and there is
    nothing to agree with.
    """
    import re
    perl_pattern = re.compile(
        r'^(.{10} )( *[0-9-]{1,4})( +[0-9.-]{6,13})( +0\.0000)'
        r'( +[^-]+-[^-]+)* *$')
    for name in ('W70', 'Pr3', 'U90'):
        lines = files.read_lines(files.resolve('RCEINP', case(name)))
        assert not any(perl_pattern.match(line) for line in lines)


def test_scale_param_scales_the_group_it_is_asked_for(tmp_path):
    folder = case('U90')
    shutil.copy(str(files.resolve('RCEINP', folder)), str(tmp_path / 'in'))
    assert run('scale_param.py', tmp_path,
               'in', 'out', 'F=0.85', '12').returncode == 0
    before = files.read_lines(tmp_path / 'in', width=80)
    after = files.read_lines(tmp_path / 'out', width=80)
    moved = 0
    for old, new in zip(before, after):
        if not rce.is_parameter_line(old):
            assert old == new
            continue
        if old[:10].strip().upper().startswith('F'):
            assert rce.parameter_value(new) == pytest.approx(
                0.85 * rce.parameter_value(old), abs=1e-4)
            moved += 1
        else:
            assert old == new
    assert moved


def test_scale_param_fixes_only_what_it_scaled(tmp_path):
    folder = case('U90')
    shutil.copy(str(files.resolve('RCEINP', folder)), str(tmp_path / 'in'))
    assert run('scale_param.py', tmp_path,
               'in', 'out', 'G=0.9', '12', 'FIX').returncode == 0
    for line in files.read_lines(tmp_path / 'out', width=80):
        if rce.is_parameter_line(line) \
                and line[:10].strip().upper().startswith('G'):
            assert abs(rce.parameter_flag(line)) == rce.FLAG_FIXED


def test_scale_param_all_fixes_every_parameter_of_the_named_sets(tmp_path):
    folder = case('U90')
    shutil.copy(str(files.resolve('RCEINP', folder)), str(tmp_path / 'in'))
    assert run('scale_param.py', tmp_path,
               'in', 'out', 'F=0.85', '1', 'ALL').returncode == 0
    section, inside = 0, False
    for line in files.read_lines(tmp_path / 'out', width=80):
        if rce.PARAMETER_HEADING in line:
            section, inside = section + 1, True
            continue
        if inside and rce.is_control_card(line):
            inside = False
        if inside and section == 1 and rce.is_parameter_line(line):
            assert abs(rce.parameter_flag(line)) == rce.FLAG_FIXED


def test_scale_param_keeps_the_sign_of_a_flag_it_fixes():
    """A negative flag means a free parameter; fixing it keeps the sign."""
    line = 'ZETA 3      -51  25302.130300      0.000000' + ' ' * 28
    assert rce.parameter_flag(rce.edit_parameter_line(line, fixed=True)) \
        == -rce.FLAG_FIXED
    plain = 'ZETA 3        0  25302.130300      0.000000' + ' ' * 28
    assert rce.parameter_flag(rce.edit_parameter_line(plain, fixed=True)) \
        == rce.FLAG_FIXED


def test_scale_param_needs_a_factor_and_a_parity_set(tmp_path):
    folder = case('U90')
    shutil.copy(str(files.resolve('RCEINP', folder)), str(tmp_path / 'in'))
    assert run('scale_param.py', tmp_path, 'in', 'out', '12').returncode == 2
    assert run('scale_param.py', tmp_path,
               'in', 'out', 'F=0.85').returncode == 2
    assert run('scale_param.py', tmp_path,
               'in', 'out', 'F=0.85', '3').returncode == 2
