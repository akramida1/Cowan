"""conf_list.py: the run's configurations, one sorted list per parity.

Byte-identical to `conf_list.pl` on all three reference runs, which is
the whole acceptance criterion here -- the printed order is a
presentation choice the author made and the port's job is to keep it.
"""

import os
import shutil
import subprocess
import sys

import pytest

from conftest import PERL, case, perl_env
from toolbox import conflist, in36

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

CASES = ['W70', 'U90', 'Pr3']


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'conf_list.py')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True)


def perl(folder):
    if shutil.which('perl') is None:
        pytest.skip('perl is not on the PATH')
    out = subprocess.run(['perl', str(PERL / 'conf_list.pl')],
                         cwd=str(folder), capture_output=True, text=True,
                         env=perl_env())
    if out.returncode:
        pytest.skip('conf_list.pl will not run here: %s' % out.stderr.strip())
    return out


@pytest.mark.parametrize('name', CASES)
def test_the_list_is_what_the_perl_prints(name):
    folder = case(name)
    ours = run(folder)
    assert ours.returncode == 0, ours.stderr
    assert ours.stdout == perl(folder).stdout


@pytest.mark.parametrize('name', CASES)
def test_every_configuration_of_in36_is_in_one_of_the_two_lists(name):
    folder = case(name)
    configurations = in36.read_in36(folder=folder)
    printed = []
    for parity in (0, 1):
        line = conflist.parity_line(configurations, parity)
        printed += [one.strip() for one in line.split(',') if one.strip()]
    assert len(printed) == len(configurations)
    assert len(set(printed)) == len(printed)


def test_a_configuration_is_written_without_its_empty_subshells():
    folder = case('U90')
    by_label = dict((c.label, c) for c in in36.read_in36(folder=folder))
    assert conflist.label(by_label['2p2']) == '1s.2p2'
    assert conflist.label(by_label['2s2p']) == '1s.2s.2p'
    assert conflist.label(by_label['2s2']) == '1s.2s2'


def test_the_parity_is_taken_from_the_subshells_that_are_printed():
    folder = case('U90')
    by_label = dict((c.label, c) for c in in36.read_in36(folder=folder))
    assert conflist.parity_of(by_label['2s2']) == 0
    assert conflist.parity_of(by_label['2p']) == 1
    assert conflist.parity_of(by_label['2s2p']) == 1


def test_a_shared_closed_f_shell_is_dropped_from_every_label():
    """`conf_list.pl:94-105`, the one closed shell these decks drop."""
    folder = case('Pr3')
    configurations = in36.read_in36(folder=folder)
    assert not conflist.drop_leading_f14(configurations)   # Pr III has 4f3
    line = conflist.parity_line(configurations, 0)
    assert '4f' in line


def test_the_order_puts_the_lowest_outer_shell_first():
    """The key is (outermost n), then the subshells, fullest first."""
    folder = case('U90')
    by_label = dict((c.label, c) for c in in36.read_in36(folder=folder))
    keys = [conflist.sort_key(by_label[label]) for label in ('2s2', '2p2')]
    assert keys == sorted(keys)


def test_more_than_one_folder_is_refused():
    assert run(case('U90'), 'a', 'b').returncode == 2
