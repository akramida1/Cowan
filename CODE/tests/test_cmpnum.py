"""cmp_num.py: tolerance-based comparison of two Fortran outputs.

`diff` is the oracle for every change that is meant to be bit-exact
(CLAUDE.md 1).  This is the oracle for the ones that are not, so the
tests here are about the two ways it can lie: calling a structural change
numeric, and calling a numeric change structural.  The real files in
WORK/ are used for the second, because the tokenizer only has to survive
the FORMATs this package actually writes.
"""

import os
import subprocess
import sys

import pytest

from conftest import WORK
from toolbox import cmpnum

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: The text outputs of the Sn VIII run in WORK/.  Small, and between them
#: they carry every FORMAT the comparator has to survive: the 1PE fields
#: of the transition list, the fixed F fields of the level list, and the
#: labels that must not be mistaken for numbers.
OUTPUTS = ['OUTG11', 'OUT36', 'OUT2', 'OUTGINE', 'TAPE2E']


def run(*arguments):
    return subprocess.run([sys.executable, os.path.join(CODE, 'cmp_num.py')]
                          + list(arguments), capture_output=True, text=True)


def output(name):
    path = WORK / name
    if not path.is_file():
        pytest.skip('WORK/%s is not present: run the Sn VIII case first' % name)
    return path


# ---------------------------------------------------------------- tokenizer

def texts(line, **kw):
    return cmpnum.tokenize(line, **kw)[3]


def test_a_configuration_label_is_not_a_number():
    """`4d7` read as an E-format number is 4e7.  This is the whole reason
    the exponent letter is E and not D, and the reason a digit following
    a letter is not a number at all."""
    assert texts('Sn7+  4d7') == ['4']
    assert texts('Sn7+  d65p') == []
    assert texts('1 (4F) 4F') == ['1', '4', '4']


def test_d_is_an_exponent_only_when_asked_for():
    assert texts('4d7') == ['4']
    assert texts('4d7', d_exponents=True) == ['4d7']


def test_a_true_jj_basis_label_does_not_overflow_to_infinity():
    """`232E1518` is a basis-state label, not 232e1518.  Read as a number
    it is +Inf, and Inf-Inf in `deviation` is a NaN that took down a whole
    Pr III OUTG11 comparison.  An exponent is at most three digits, so the
    E spelling now tokenizes exactly as the D spelling beside it."""
    assert texts('  230D1518  232E1518') == ['230', '232']
    assert texts('  230D1518  232E1518', d_exponents=True) == ['230', '232']
    assert cmpnum.tokenize('232E1518')[1] == [232.0]
    assert texts('4.18E+03') == ['4.18E+03']
    assert texts('1.0E-105') == ['1.0E-105']


def test_the_1pe_fields_of_the_transition_list():
    line = ('    1    -29.947  4.5  1 (4F) 4F      418.430    238.989  '
            '-1.485 3.821E+09  0.09     1')
    assert texts(line) == ['1', '-29.947', '4.5', '1', '4', '4', '418.430',
                           '238.989', '-1.485', '3.821E+09', '0.09', '1']


def test_a_bare_sign_exponent_is_read_as_two_packed_fields():
    """Ambiguous by construction; resolved towards the common case.  Two
    fixed-width fields whose second is negative are everywhere in this
    output, a three-digit exponent is not."""
    assert texts('1.234-105') == ['1.234', '-105']


def test_packed_negative_fields_split_on_the_sign():
    assert texts('  418.430-238.989') == ['418.430', '-238.989']


def test_a_run_of_blanks_is_one_blank_in_the_skeleton():
    """A number that gains a digit eats a blank of its field padding.
    That must not read as a structural change as well as a numeric one."""
    a = cmpnum.tokenize('  1.0  2.0')[0]
    b = cmpnum.tokenize(' 11.0 22.0')[0]
    assert a == b


# ---------------------------------------------------------------- structure

def compare(a, b, **kw):
    report = cmpnum.Report('a', 'b')
    return cmpnum.compare_lines(a.splitlines(), b.splitlines(), report, **kw)


def test_changed_text_is_structural_not_numeric():
    report = compare('ELEC DIP SPECTRUM 1.5', 'ELEC QUAD SPECTRUM 1.5')
    assert report.structural
    assert not report.offenders


def test_a_lost_line_is_structural():
    report = compare('1.0\n2.0\n', '1.0\n')
    assert report.structural
    assert report.structural[0][0] == 0        # the file, not a line


def test_an_overflowed_field_is_structural():
    """A field that will not fit prints as asterisks.  That is a lost
    number, and it must never come back as a passing comparison."""
    report = compare(' 1.5  2.5', ' 1.5  ****')
    assert report.structural
    assert 'numbers on the line' in report.structural[0][1]


def test_a_line_that_is_structural_contributes_no_numbers():
    report = compare('A 1.0\nB 2.0', 'C 1.0\nB 2.0')
    assert report.compared == 1


# ---------------------------------------------------------------- tolerance

def test_within_tolerance_passes_and_still_enters_the_histogram():
    """The distribution is the deliverable, not the pass/fail.  A
    histogram of only the failures would hide the shape of the change."""
    report = compare(' 1.000000', ' 1.000005', rtol=1e-5)
    assert report.ok
    assert report.compared == 1
    assert report.histogram == {-6: 1}


def test_outside_tolerance_is_reported_with_its_line_and_column():
    report = compare('x\n  1.000  2.000', 'x\n  1.000  2.002', rtol=1e-4)
    assert len(report.offenders) == 1
    _, _, line, column, a, b, _, _ = report.offenders[0]
    assert (line, column, a, b) == (2, 10, 2.0, 2.002)


def test_either_tolerance_is_enough():
    """A field printed to three decimals has no meaningful relative
    deviation once it is near zero, which is what --atol is for."""
    report = compare(' 0.001', ' 0.002', rtol=1e-9, atol=1e-2)
    assert report.ok
    report = compare(' 0.001', ' 0.002', rtol=1e-9, atol=0.0)
    assert not report.ok


def test_two_zeros_are_exact_not_undefined():
    report = compare(' 0.000', ' 0.000')
    assert report.ok
    assert report.histogram == {None: 1}


def test_timing_lines_are_ignored_by_default():
    """CLAUDE.md 1 excludes exactly these from the OUTG11 acceptance
    test; they are wall clock and never agree."""
    a = '0FINISHED PLEV AT    0.000 MIN\n           TIME=    0.000 MIN'
    b = '0FINISHED PLEV AT    9.750 MIN\n           TIME=    3.125 MIN'
    assert compare(a, b).ok
    assert compare(a, b, ignores=()).offenders


def test_rms_is_over_every_pair_compared():
    report = compare(' 1.0  1.0', ' 1.0  2.0')
    assert report.compared == 2
    assert report.rms == pytest.approx((0.5 ** 2 / 2) ** 0.5)


# ---------------------------------------------------------------- real files

@pytest.mark.parametrize('name', OUTPUTS)
def test_a_file_compares_exactly_with_itself(name):
    """Stage A of the plan is byte-identical output, so this is the case
    the tool has to get right before any of the rest matters."""
    path = output(name)
    report = cmpnum.compare_files(path, path)
    assert report.ok
    assert not report.structural


@pytest.mark.parametrize('name', OUTPUTS)
def test_every_line_of_a_real_output_tokenizes_without_a_long_token(name):
    """Packed fields that run together are read as one number and their
    deviations are meaningless.  None of the FORMATs in use produce one;
    if that ever changes, this is where it shows up."""
    path = output(name)
    report = cmpnum.compare_files(path, path)
    assert report.long_tokens == 0


def test_an_unformatted_file_gets_an_md5_answer_not_a_silent_skip():
    """TAPE2N is written as one unformatted record per configuration
    (RCN36K.F:2968), mixing integers, characters and doubles.  It cannot
    be tokenized -- but a skip that still exits 0 is exactly the failure
    this whole exercise guards against."""
    path = output('TAPE2N')
    assert cmpnum.looks_binary(path)
    report = cmpnum.compare_files(path, path)
    assert report.binary and report.binary_same and report.ok


# ---------------------------------------------------------------- the script

def test_the_script_exits_zero_on_an_identical_pair():
    path = str(output('OUTG11'))
    done = run('--quiet', path, path)
    assert done.returncode == 0, done.stderr


def test_the_script_exits_two_on_a_structural_mismatch(tmp_path):
    a = tmp_path / 'a'
    b = tmp_path / 'b'
    a.write_text('ELEC DIP SPECTRUM  1.5\n')
    b.write_text('ELEC QUAD SPECTRUM  1.5\n')
    done = run(str(a), str(b))
    assert done.returncode == 2
    assert 'STRUCTURAL' in done.stdout


def test_the_script_exits_one_on_a_numeric_difference(tmp_path):
    a = tmp_path / 'a'
    b = tmp_path / 'b'
    a.write_text('  1.000\n')
    b.write_text('  1.002\n')
    assert run(str(a), str(b)).returncode == 1
    assert run('--rtol', '1e-2', str(a), str(b)).returncode == 0


def test_the_script_compares_a_folder_and_reports_what_is_missing(tmp_path):
    old = tmp_path / 'old'
    new = tmp_path / 'new'
    old.mkdir()
    new.mkdir()
    (old / 'OUT36').write_text('  1.000\n')
    (new / 'OUT36').write_text('  1.000\n')
    (old / 'OUT2').write_text('  2.000\n')
    done = run('--quiet', str(old), str(new))
    assert done.returncode == 2
    assert 'OUT2' in done.stdout and 'MISSING' in done.stdout


def test_the_script_rejects_a_bad_command_line():
    assert run('one-name-only').returncode == 2
    assert run('--rtol', 'wide', 'a', 'b').returncode == 2
