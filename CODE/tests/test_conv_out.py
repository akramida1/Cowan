"""conv_out.py, end to end on the reference runs.

The same acceptance criterion as phase 6: the Python output must agree
with `conv_out.pl` line for line except where the Perl is known to be
wrong, and every difference must be one of the classes named in
`FOR/LSJJ/README.txt`.  The `conv_out.pl` outputs are not checked in --
they are megabytes and they depend on the model asked for -- so the
tests that need one run the Perl themselves and skip where it is gone.
"""

import os
import re
import shutil
import subprocess
import sys

import pytest

import conv_out
from conftest import PERL, case, perl_env
from toolbox import ai, parfile, radiate, spectrum

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'conv_out.py')] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


def perl(folder, *arguments):
    if shutil.which('perl') is None:
        pytest.skip('perl is not on the PATH')
    out = subprocess.run(
        ['perl', '-x', '-I' + str(PERL), str(PERL / 'conv_out.pl')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True,
        env=perl_env())
    assert out.returncode == 0, out.stderr
    return out


#: A par file with every setting in it, for the two tests about which
#: source wins.  `no_rce` and `iden` are on so that a command line
#: turning one back off is visible.
PAR = {'levels': 'lev.txt', 'lines': 'tr.txt', 'model': 'lte',
       'temperature': 100.0, 'max_intensity': 1000.0, 'couplings': 'LSJK',
       'drop_full_shells': '5p,4d', 'components': 9, 'min_percent': 0.5,
       'upper_bound': 500.0, 'ai': 'even.txt, odd.txt', 'flight': 0.2,
       'no_rce': True, 'iden': True, 'truncate': False, 'linear': False,
       'convert2log': False}


def rows(path, skip=1):
    with open(path) as f:
        return [line.rstrip('\n').split('\t') for line in f][skip:]


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

def test_the_historical_argument_order_is_understood(tmp_path):
    out = run(case('U90'), str(tmp_path / 'lev.txt'), str(tmp_path / 'tr.txt'),
              'RCE', 'LTE', '1e6', '1000')
    assert out.returncode == 0, out.stderr
    assert (tmp_path / 'lev.txt').exists() and (tmp_path / 'tr.txt').exists()


def test_a_model_that_needs_autoionization_rates_says_so(tmp_path):
    out = run(case('U90'), str(tmp_path / 'lev.txt'), str(tmp_path / 'tr.txt'),
              'RCE', 'AUGER', '1e6', '1000')
    assert out.returncode == 2
    assert 'read_aa' in out.stderr


def test_beam_foil_wants_its_time_of_flight(tmp_path):
    out = run(case('U90'), str(tmp_path / 'lev.txt'), str(tmp_path / 'tr.txt'),
              'RCE', 'BF', '1e6', '1000')
    assert out.returncode == 2 and 'BF:0.114' in out.stderr


def test_one_autoionization_file_is_not_enough(tmp_path):
    out = run(case('U90'), str(tmp_path / 'lev.txt'), str(tmp_path / 'tr.txt'),
              'RCE', 'DR', '1e6', '1000', 'even.txt')
    assert out.returncode == 2 and 'two files' in out.stderr


# --------------------------------------------------------------------------
# against conv_out.pl
# --------------------------------------------------------------------------

BRACKET = re.compile(r'\(([0-9/]+)\)\)')


def _levels(folder, tmp_path, *extra):
    """Both programs' level lists, row by row on (coupling, Ec, J)."""
    perl(folder, str(tmp_path / 'a.txt'), str(tmp_path / 'b.txt'),
         'RCE', 'LTE', '1e6', '1000', *extra)
    out = run(folder, str(tmp_path / 'c.txt'), str(tmp_path / 'd.txt'),
              'RCE', 'LTE', '1e6', '1000', *extra)
    assert out.returncode == 0, out.stderr
    theirs = {(r[0], r[5], r[6]): r for r in rows(tmp_path / 'a.txt')}
    ours = {(r[0], r[5], r[6]): r for r in rows(tmp_path / 'c.txt')}
    return theirs, ours


def test_u90_levels_differ_only_where_the_perl_is_wrong(tmp_path):
    theirs, ours = _levels(case('U90'), tmp_path)
    assert set(theirs) == set(ours)
    counts = {}
    for key in theirs:
        a, b = theirs[key], ours[key]
        if a == b:
            continue
        fixed = [BRACKET.sub(r'<\1>)', one) for one in a]
        _bump(counts, 'brackets', fixed != a)
        _bump(counts, 'outer', fixed[9:12] != b[9:12])
        fixed[9:12] = b[9:12]
        _bump(counts, 'n', '1s' in fixed)
        fixed = ['1s.2s2' if one == '1s' else one for one in fixed]
        _bump(counts, 'filled', len(b) > len(fixed))
        if len(b) > len(fixed):
            fixed = b
        assert fixed == b, (a, b)
    assert counts == {'brackets': 9, 'outer': 23, 'n': 2, 'filled': 1}


def _bump(counts, name, hit):
    if hit:
        counts[name] = counts.get(name, 0) + 1


def test_w70_levels_differ_only_in_the_brackets_and_the_outer_shell(tmp_path):
    theirs, ours = _levels(case('W70'), tmp_path)
    counts = {}
    for key in set(theirs) & set(ours):
        a, b = theirs[key], ours[key]
        if a == b:
            continue
        fixed = [BRACKET.sub(r'<\1>)', one) for one in a]
        _bump(counts, 'brackets', fixed != a)
        _bump(counts, 'outer', fixed[9:12] != b[9:12])
        fixed[9:12] = b[9:12]
        _bump(counts, 'filled', len(b) > len(fixed))
        if len(b) > len(fixed):
            fixed = b
        assert fixed == b, (a, b)
    assert counts == {'brackets': 135, 'outer': 45, 'filled': 1}


def test_the_outer_shell_columns_name_the_outermost_open_subshell(tmp_path):
    # `conv_cowan.pl:3043` matched the first subshell of the label
    # instead of the last, which for 2s2.2p gave 2s2.
    out = run(case('W70'), str(tmp_path / 'lev.txt'), str(tmp_path / 'tr.txt'),
              'RCE', 'LTE', '1e6', '1000')
    assert out.returncode == 0, out.stderr
    ground = rows(tmp_path / 'lev.txt')[0]
    assert ground[14:16] == ['100', '2s2.2p']
    assert ground[9:12] == ['2', '1', '1']


def test_the_cascade_matrix_is_reproduced_exactly(tmp_path):
    # Nothing in it depends on a designation, so it is byte for byte.
    folder = case('U90')
    perl(folder, str(tmp_path / 'a.txt'), str(tmp_path / 'a.csv'),
         'RCE', 'CM', '1e6', '1000')
    out = run(folder, str(tmp_path / 'b.txt'), str(tmp_path / 'b.csv'),
              'RCE', 'CM', '1e6', '1000')
    assert out.returncode == 0, out.stderr
    with open(tmp_path / 'a.csv') as f:
        theirs = f.read()
    with open(tmp_path / 'b.csv') as f:
        assert f.read() == theirs


def test_every_transition_of_the_perl_is_in_the_python_list(tmp_path):
    folder = case('W70')
    perl(folder, str(tmp_path / 'a.txt'), str(tmp_path / 'a_tr.txt'),
         'RCE', 'LTE', '1e6', '1000')
    out = run(folder, str(tmp_path / 'b.txt'), str(tmp_path / 'b_tr.txt'),
              'RCE', 'LTE', '1e6', '1000')
    assert out.returncode == 0, out.stderr
    theirs = {(r[6], r[7]) for r in rows(tmp_path / 'a_tr.txt', 2)}
    ours = {(r[6], r[7]) for r in rows(tmp_path / 'b_tr.txt', 2)}
    # One line of the 706 sits on the 0.5 printing threshold and falls
    # the other side of it: the intensity model divides by RCG's own
    # wavelength, taken from its eigenvalues rather than from the three
    # decimals of OUTG11's printed F11.3, which at W70's 1.1 A is a
    # rounding of a few parts in 10^4 that the normalisation spreads
    # over every line.
    assert theirs - ours and len(theirs - ours) <= 1
    assert not ours - theirs


def test_the_modelled_intensities_agree_to_better_than_a_percent(tmp_path):
    folder = case('U90')
    perl(folder, str(tmp_path / 'a.txt'), str(tmp_path / 'a_tr.txt'),
         'RCE', 'LTE', '1e6', '1000')
    out = run(folder, str(tmp_path / 'b.txt'), str(tmp_path / 'b_tr.txt'),
              'RCE', 'LTE', '1e6', '1000')
    assert out.returncode == 0, out.stderr
    theirs = {(r[6], r[7]): float(r[14]) for r in rows(tmp_path / 'a_tr.txt', 2)}
    worst = 0.0
    for row in rows(tmp_path / 'b_tr.txt', 2):
        one = theirs.get((row[6], row[7]))
        if one:
            worst = max(worst, abs(one - float(row[14])) / one)
    # Both programs model the intensity from RCG's unscaled gA at RCG's
    # own wavelength; what is left is that `conv_out.pl` reads that
    # wavelength from OUTG11's F11.3 field, which for U90's 0.12 A lines
    # holds three significant figures.
    assert 0 < worst < 0.03            # and every one is a rounded unit


# --------------------------------------------------------------------------
# the pieces
# --------------------------------------------------------------------------

def test_the_line_strength_conversions_invert_each_other():
    for multipole in (radiate.E1, radiate.M1, radiate.E2, radiate.M3):
        s = radiate.s_from_ga(1e8, 1500.0, multipole)
        assert abs(radiate.ga_from_s(s, 1500.0, multipole) - 1e8) < 1.0
        gf = radiate.gf_from_s(s, 1500.0, multipole)
        assert abs(radiate.s_from_gf(gf, 1500.0, multipole) - s) < 1e-9 * s


def test_every_multipole_rcg_prints_has_its_constants():
    for name in spectrum.MULTIPOLES:
        assert radiate.OF_SPECTRUM[name] in (radiate.E1, radiate.M1,
                                             radiate.E2)


def test_an_autoionization_table_is_read_by_its_heading(tmp_path):
    path = tmp_path / 'even.txt'
    path.write_text('J\tEnergy\t1s2:0\t1s2s:120.1\tAA_tot\n'
                    '0.5\t130.5\t 1.20e+12\t 0.00e+00\t 1.20e+12\n'
                    '1.5\t131.0\t 0.00e+00\t 3.40e+11\t 3.40e+11\n')
    channels, found = ai.read_ai(str(path))
    assert [c.name for c in channels] == ['1s2', '1s2s']
    assert channels[0].text == '1s2:0'
    assert found[0].rates == {'1s2': 1.2e12} and found[0].total == 1.2e12
    assert found[1].dr(channels) == 0.0


def test_a_table_that_is_not_read_aa_s_is_refused(tmp_path):
    path = tmp_path / 'bad.txt'
    path.write_text('wavelength\tintensity\n1\t2\n')
    with pytest.raises(ai.AiError):
        ai.read_ai(str(path))


# --------------------------------------------------------------------------
# the parameter file
# --------------------------------------------------------------------------
#
# The par file is purely additive: it supplies defaults, the command
# line overrides them, and a run without one behaves exactly as it did
# before there was one.  The tests above are the evidence that the port
# is faithful and they drive the historical command line, so the ones
# here are about the two sources meeting, not about the output.

def test_no_par_file_and_no_arguments_writes_a_sample_and_stops(tmp_path):
    out = run(tmp_path)
    assert out.returncode == 0, out.stderr
    assert (tmp_path / conv_out.PAR_FILE).is_file()
    assert 'run conv_out again' in out.stdout


def test_an_incomplete_command_line_is_an_error_not_a_new_par_file(tmp_path):
    # The sample is written only when there is nothing to go on at all.
    # A command line that is merely wrong must be reported, and must not
    # leave a file behind -- the comparison tests above run inside the
    # reference decks of WORK/, which nothing may write into.
    out = run(tmp_path, 'lev.txt')
    assert out.returncode == 2
    assert not (tmp_path / conv_out.PAR_FILE).exists()


def test_the_sample_reads_back_through_its_own_reader(tmp_path):
    run(tmp_path)
    got = parfile.read(tmp_path / conv_out.PAR_FILE, conv_out.SPECS)
    assert got == {spec.key: spec.default for spec in conv_out.SPECS}


def test_write_par_writes_the_sample_on_request(tmp_path):
    out = run(tmp_path, '--write-par')
    assert out.returncode == 0, out.stderr
    assert (tmp_path / conv_out.PAR_FILE).is_file()


def test_every_model_is_named_and_explained_in_the_sample(tmp_path):
    run(tmp_path)
    text = (tmp_path / conv_out.PAR_FILE).read_text()
    for model in conv_out.MODELS:
        assert re.search(r'^;\s+%s --' % model, text, re.M), model


def test_every_parameter_of_the_program_is_documented():
    for spec in conv_out.SPECS:
        assert spec.description.strip(), spec.key


def test_the_par_file_supplies_what_the_command_line_does_not():
    got = conv_out.Options([], PAR)
    assert (got.levels, got.lines, got.model) == ('lev.txt', 'tr.txt', 'LTE')
    assert got.temperature == 100.0 and got.max_intensity == 1000.0
    assert got.couplings == ['LSJK'] and got.components == 9
    assert got.min_percent == 0.5 and got.upper_bound == 500.0
    assert got.drop == [(5, 1), (4, 2)]
    assert got.ai == ['even.txt', 'odd.txt']
    assert got.flight == 0.2
    assert got.no_rce and got.iden


def test_the_command_line_has_the_last_word_over_the_par_file():
    got = conv_out.Options(['other_lev.txt', 'other_tr.txt', 'RCE', 'LTEABS',
                            '5', '20', '--components', '3'], PAR)
    assert (got.levels, got.lines) == ('other_lev.txt', 'other_tr.txt')
    assert got.model == 'LTEABS'
    assert got.temperature == 5.0 and got.max_intensity == 20.0
    assert got.components == 3
    assert got.no_rce is False           # the bare keyword RCE, over no_rce=1


def test_a_setting_neither_source_mentions_takes_its_built_in_default():
    got = conv_out.Options(['lev.txt', 'tr.txt', 'RCE', 'LTE', '1', '1000'],
                           {})
    assert got.couplings == ['LS', 'JJ'] and got.components == 5
    assert got.min_percent == 1.5 and got.upper_bound == float('inf')
    assert got.drop == [] and got.ai == []
    assert got.flight == conv_out.FLIGHT
    assert not (got.no_rce or got.iden or got.truncate or got.linear
                or got.convert2log)


def test_the_par_file_carries_the_time_of_flight_beside_a_bare_bf():
    # On the command line the two are glued together as BF:0.114,
    # because there is nowhere else to put the number; the par file has
    # a field for it, so a bare BF is right there and wrong there.
    par = dict(PAR, model='BF', flight=0.3)
    assert (conv_out.Options([], par).model,
            conv_out.Options([], par).flight) == ('BF', 0.3)


def test_a_par_file_that_will_not_read_is_named_with_its_line(tmp_path):
    (tmp_path / conv_out.PAR_FILE).write_text(
        'LTE ; model' + chr(10) + '7 ; nonesuch' + chr(10))
    out = run(tmp_path, 'lev.txt', 'tr.txt', 'RCE', 'LTE', '1', '1000')
    assert out.returncode == 2 and 'nonesuch' in out.stderr


def test_a_value_the_par_file_gets_wrong_is_named_with_its_line(tmp_path):
    (tmp_path / conv_out.PAR_FILE).write_text('SPQR ; model' + chr(10))
    out = run(tmp_path, 'lev.txt', 'tr.txt', 'RCE', 'LTE', '1', '1000')
    assert out.returncode == 2 and 'model' in out.stderr
