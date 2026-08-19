"""The Case object and printout.py, end to end on the reference runs.

The acceptance criterion for this phase is the plan's: the Python output
must agree with `printout.pl` line for line except where the Perl is
known to be wrong.  Those four defects are asserted here as *changes*,
not tolerated as differences.
"""

import os
import re
import subprocess
import sys

import pytest

from conftest import case
from toolbox import case as case_mod
from toolbox.case import Case, CaseError

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
_CACHE = {}


def read(name, couplings=('LS',), **kw):
    key = (name, tuple(couplings), tuple(sorted(kw.items())))
    if key not in _CACHE:
        _CACHE[key] = Case.read(case(name), couplings=couplings, **kw)
    return _CACHE[key]


def run(folder, *arguments):
    """printout.py in a working directory, as the .bat wrapper runs it."""
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'printout.py')] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


def rows(path):
    with open(path) as f:
        return [line.rstrip('\n').split('\t') for line in f][1:]


# --------------------------------------------------------------------------
# the Case object
# --------------------------------------------------------------------------

def test_a_case_joins_its_files():
    one = read('U90', ('LS', 'JJ'))
    assert len(one.blocks) == 6
    assert len(one.identified) == 19
    assert one.shells[(1, 1)][0].nl == '1s'
    assert one.card(1, 1).label == '2s'


def test_a_coupling_the_deck_never_asked_for_is_a_hard_error():
    # The plan is explicit: name the coupling and the column, and stop.
    with pytest.raises(CaseError) as e:
        Case.read(case('U90'), couplings=('LSJK',))
    assert 'LSJK' in str(e.value) and 'columns 31-39' in str(e.value)


def test_true_jj_is_asked_for_only_where_it_is_wanted():
    # W70 is the one reference deck with KCPLD(8)=8.
    assert 'TRUEjj' in read('W70', ('LS', 'TRUEjj')).deck.main.couplings
    assert 'TRUEjj' not in read('U90', ('LS',)).deck.main.couplings


def test_the_energies_come_from_levels1_not_rceout():
    # RCEOUT prints F10.3 and no g factor; LEVELS1 prints F11.4 and one.
    one = read('U90', ('LS', 'JJ'))
    level = one.identified[0].level
    assert level.g is not None
    assert abs(one.cset(1).levels[0].e_calc
               - one.identified[0].level.e_calc) < 1e-9


def test_a_one_state_block_is_recovered_rather_than_dropped():
    # `WORK/W70` has a J block of one basis state, for which RCG prints
    # no JJ section at all.  Amplitude 1 by construction.
    one = read('W70', ('LS', 'JJ', 'TRUEjj'))
    assert one.amplitudes('JJ', 2, 5.5, 1) == [(1.0, 1)]
    found = one.designation('JJ', 2, 5.5, 1)
    assert found.configuration == '2s.2p.(2P*<3/2>).(2).4f'
    assert found.term == '(2,7/2)'


# --------------------------------------------------------------------------
# printout.py against the Perl's own output
# --------------------------------------------------------------------------

def test_u90_reproduces_the_perl_except_where_the_perl_is_wrong(tmp_path):
    folder = case('U90')
    out = run(folder, str(tmp_path / 'levels.txt'), str(tmp_path / 'p.txt'),
              '--couplings', 'LS,JJ')
    assert out.returncode == 0, out.stderr
    with open(folder / 'levels.txt') as f:
        theirs = [l.rstrip('\n\r') for l in f][1:]
    ours = [('\t'.join(r)) for r in rows(tmp_path / 'levels.txt')]
    assert len(ours) == len(theirs)

    changed = [(a, b) for a, b in zip(theirs, ours) if a != b]
    assert changed, 'nothing changed at all, which cannot be right'
    brackets = reattached = filled = 0
    for a, b in changed:
        fixed = re.sub(r'\(([0-9/]+)\)\)', r'<\1>)', a)
        if fixed != a:
            brackets += 1
        if '\t1s\t' in fixed:
            reattached += 1
            fixed = fixed.replace('\t1s\t', '\t1s.2s2\t')
        if fixed != b and b.startswith(fixed):
            filled += 1                     # the level the Perl left blank
            fixed = b
        assert fixed == b, (a, b)
    # Every one of the four defects the plan named, and nothing else.
    assert (brackets, reattached, filled) == (9, 2, 1)


def test_the_blank_level_the_perl_could_not_name_is_named(tmp_path):
    out = run(case('W70'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'JJ')
    assert out.returncode == 0, out.stderr
    # J=11/2 of the even set is one basis state, so RCG printed no JJ
    # section for it and the Perl left the whole row blank.
    named = [r for r in rows(tmp_path / 'levels.txt')
             if r[6].strip() == '5.5' and len(r) > 10]
    assert ([(r[8], r[9], r[10]) for r in named]
            == [('100', '2s.2p.(2P*<3/2>).(2).4f', '(2,7/2)')])


def test_every_label_of_a_run_names_the_same_number_of_electrons(tmp_path):
    # The author's rule: only --drop-full-shells may break this.
    out = run(case('W70'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'LS')
    assert out.returncode == 0, out.stderr
    counts = set()
    for row in rows(tmp_path / 'levels.txt'):
        for at in range(9, len(row), 3):
            counts.add(_electrons(row[at]))
    assert len(counts) == 1, sorted(counts)


def _electrons(configuration):
    total = 0
    for token in configuration.split('.'):
        if token.startswith('('):
            continue
        digits = ''
        for ch in reversed(token):
            if not ch.isdigit():
                break
            digits = ch + digits
        total += int(digits) if digits else 1
    return total


def test_drop_full_shells_drops_only_where_the_shell_is_closed(tmp_path):
    # Pr III has one configuration with an open 5p, so 5p6 is in every
    # label by default; --drop-full-shells 5p takes it out of all the
    # others and leaves 5p5 where it stands.
    out = run(case('Pr3'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'LS',
              '--drop-full-shells', '5p')
    assert out.returncode == 0, out.stderr
    labels = {r[9] for r in rows(tmp_path / 'levels.txt') if len(r) > 9}
    assert '4f3' in labels
    assert not any('5p6' in one for one in labels)
    assert any('5p5' in one for one in labels)


def test_the_parameter_table_names_every_parameter(tmp_path):
    out = run(case('U90'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'LS')
    assert out.returncode == 0, out.stderr
    table = rows(tmp_path / 'p.txt')
    assert len(table) == 21
    by_name = {r[3]: r for r in table}
    assert by_name['Eav(2s)'][1] == '1s2.2s'
    assert by_name['zeta(2p)'][6] == ''       # not linked to a group
    assert by_name['RD1(2s2s,2p2p)'][6] == '1'
    # A slot name that occurs twice must not collapse: the odd set has
    # two zeta(2p), one per configuration, with different uncertainties.
    odd = [r for r in table if r[0] == 'o' and r[3] == 'zeta(2p)']
    assert len(odd) == 2 and odd[0][4] != odd[1][4]


@pytest.mark.slow
def test_the_rows_are_ordered_by_the_rcg_eigenvalue(tmp_path):
    # Not by the fitted energy.  RCE lists its levels in the order of its
    # own Ecalc and RCG recomputes the same quantity from the rounded
    # parameters RCE punched back, so the two orders agree except where
    # a pair is close enough for the rounding to swap them -- and it is
    # RCG's order that printout.pl and RCE's own output follow.
    out = run(case('Pr3'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'LS')
    assert out.returncode == 0, out.stderr
    one = read('Pr3', ('LS',))
    key = {('%4.1f' % x.j, '%12.4f' % x.level.e_calc):
           one.blocks[(x.kset, x.j)].eigenvalues[x.rcg - 1]
           for x in one.identified}

    swapped = 0
    for parity in ('e', 'o'):
        block = [r for r in rows(tmp_path / 'levels.txt') if r[1] == parity]
        eigenvalues = [key[(r[6], r[4])] for r in block]
        assert eigenvalues == sorted(eigenvalues)
        fitted = [float(r[4]) for r in block]
        swapped += sum(1 for a, b in zip(fitted, fitted[1:]) if b < a)
    # ... and the two keys really are different keys: fourteen pairs of
    # this run print in descending fitted energy because that is the
    # order of the eigenvalues.
    assert swapped == 14


# --------------------------------------------------------------------------
# the cross-folder merge
# --------------------------------------------------------------------------

def test_a_reordered_folder_merges_by_energy_order(tmp_path):
    out = run(case('U90'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'LS,JJ',
              '--from', 'JJ=reordered')
    assert out.returncode == 0, out.stderr
    assert 'matched by energy-order' in out.stdout
    assert len(rows(tmp_path / 'levels.txt')) == 38


@pytest.mark.slow
def test_the_jk_folder_supplies_lsjk_for_the_parent_s_levels(tmp_path):
    # The merge that used to be hand work.
    out = run(case('Pr3'), str(tmp_path / 'levels.txt'),
              str(tmp_path / 'p.txt'), '--couplings', 'LS,LSJK',
              '--from', 'LSJK=JK')
    assert out.returncode == 0, out.stderr
    assert 'matched by eigenvector' in out.stdout
    table = rows(tmp_path / 'levels.txt')
    assert len(table) == 2 * 5626
    jk = [r for r in table if r[0] == 'LSJK' and len(r) > 10]
    # Every LSJK component carries a K term, bar the one J block that
    # holds a single basis state: RCG prints no LSJK section for it, so
    # it is named in LS and the run says so.
    plain = [r for r in jk if '[' not in r[10]]
    assert len(jk) > 5000 and len(plain) == 1
    assert 'has one basis state' in out.stdout


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

def test_an_unknown_coupling_is_refused_by_name(tmp_path):
    out = run(case('U90'), str(tmp_path / 'levels.txt'), '--couplings', 'SL')
    assert out.returncode == 2
    assert 'SL' in out.stderr


def test_from_must_name_a_requested_coupling(tmp_path):
    out = run(case('U90'), str(tmp_path / 'levels.txt'),
              '--couplings', 'LS', '--from', 'JJ=reordered')
    assert out.returncode == 2
    assert '--from' in out.stderr
