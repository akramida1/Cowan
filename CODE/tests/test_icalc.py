"""extract_Icalc.py: the classification columns of a calc_TP line list.

Three things here can go wrong quietly and are what these tests are
for: the units, because `calc_TP` writes kK and this writes cm-1 and a
factor of 1000 in the Boltzmann exponent is not visible in the output;
the statistical weight, because gA and A differ by 2*J2+1 and both are
plausible numbers; and the filter, because a line without an identifier
or without an experimental energy has to be dropped rather than written
with a blank or a zero in it.

The source rows are built here rather than taken from a run: the shape
that matters is the awkward one `calc_TP` and `insert_level_ids` leave
behind -- trailing tabs stripped off, so a row can be shorter than the
header, and a blank in any of the columns this reads.
"""

import math
import os
import subprocess
import sys

import pytest

from toolbox import fundconst
from toolbox import icalc

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

#: A `calc_TP` header with `insert_level_ids`' identifiers on it, cut
#: down to the columns this tool reads and the ones it has to step over.
HEADER = ['type', 'conf1', 't1', 'J1', 'conf2', 't2', 'J2', 'E_low', 'E_up',
          'lid1', 'lid2', 'id1', 'id2', 'wl_c(A)', 'E_low_exp', 'E_up_exp',
          'A(s-1)', 'S', 'std_A%']


def row(id1='059003.000065', id2='059003.000094', j2='6.5',
        low='13.3521', up='138.3427', rate='1.29676E+06', spread='143.07'):
    return ['E1', '4f3.5p6', '4I*', '4.5', '4f2.5p6.6f', '4I', j2,
            '13.3567', '138.3132', 'Lev_1', 'Lev_9', id1, id2, '800.2785',
            low, up, rate, '4.5887E-03', spread]


def write(path, rows, ending='\n'):
    """Write a table the way calc_TP does: trailing tabs stripped."""
    with open(str(path), 'w', newline='') as handle:
        for one in rows:
            handle.write('\t'.join(one).rstrip('\t') + ending)
    return str(path)


def read(path):
    with open(str(path), newline='') as handle:
        return [line.rstrip('\r\n').split('\t')
                for line in handle if line.strip('\r\n')]


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'extract_Icalc.py')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True)


def extract(tmpdir, rows, *arguments):
    write(tmpdir.join('tp_E1.txt'), [HEADER] + rows)
    done = run(tmpdir, *arguments)
    assert done.returncode == 0, done.stderr
    return read(tmpdir.join('Icalc.txt')), done


# ----------------------------------------------------------------- shape

def test_default_names_and_header(tmpdir):
    """tp_E1.txt in, Icalc.txt out, and six columns without the options."""
    out, _ = extract(tmpdir, [row()])
    assert out[0] == list(icalc.COLUMNS)
    assert len(out) == 2


def test_the_six_columns(tmpdir):
    """gA carries the weight, the energies are in cm-1, u%gA is copied."""
    out, _ = extract(tmpdir, [row()])
    got = out[1]
    assert got[0] == '059003.000065' and got[1] == '059003.000094'
    assert got[2] == '124990.6'                  # 138.3427 - 13.3521, kK
    # The column is written %.5E, so it is compared to that precision.
    assert float(got[3]) == pytest.approx(14 * 1.29676e6, rel=1e-5)  # (2J+1)*A
    assert got[4] == '143.07'
    assert got[5] == '138342.7'


def test_energies_are_not_padded(tmpdir):
    """A whole number of cm-1 comes out without a tail of zeros."""
    out, _ = extract(tmpdir, [row(low='0.0000', up='100.0000')])
    assert out[1][2] == '100000' and out[1][5] == '100000'


def test_short_rows(tmpdir):
    """A row with its trailing tabs stripped is read, not tripped over."""
    short = row(spread='')
    while not short[-1]:
        short.pop()
    out, _ = extract(tmpdir, [short])
    assert len(out) == 2 and out[1][4] == ''


def test_crlf_is_kept(tmpdir):
    """The output keeps the line ending the source had."""
    write(tmpdir.join('tp_E1.txt'), [HEADER, row()], ending='\r\n')
    assert run(tmpdir).returncode == 0
    with open(str(tmpdir.join('Icalc.txt')), 'rb') as handle:
        assert b'\r\n' in handle.read()


def test_output_name(tmpdir):
    """The second file name is the output."""
    write(tmpdir.join('lines.txt'), [HEADER, row()])
    done = run(tmpdir, 'lines.txt', 'other.txt')
    assert done.returncode == 0, done.stderr
    assert os.path.exists(str(tmpdir.join('other.txt')))
    assert not os.path.exists(str(tmpdir.join('Icalc.txt')))


# ---------------------------------------------------------------- filter

@pytest.mark.parametrize('bad, why', [
    (dict(id1=''), 'without an identifier'),
    (dict(id2=''), 'without an identifier'),
    (dict(low=''), 'without an experimental energy'),
    (dict(up=''), 'without an experimental energy'),
    (dict(rate=''), 'without a transition rate'),
    (dict(low='138.3427'), 'whose two experimental energies are equal'),
])
def test_lines_left_out(tmpdir, bad, why):
    """Each of the five reasons drops the line and is reported."""
    out, done = extract(tmpdir, [row(), row(**bad)])
    assert len(out) == 2                        # the header and the good one
    assert '2 transitions read, 1 written.' in done.stdout
    assert why in done.stdout


def test_blank_uncertainty_stays_blank(tmpdir):
    """calc_TP's blank std_A% is a line with no spread, not a zero."""
    out, _ = extract(tmpdir, [row(spread='')])
    assert out[1][4] == ''


def test_no_std_column_at_all(tmpdir):
    """A calc_TP run with no trials has no std_A%; the table still works."""
    header = [name for name in HEADER if name != 'std_A%']
    write(tmpdir.join('tp_E1.txt'), [header, row()[:-1]])
    done = run(tmpdir)
    assert done.returncode == 0, done.stderr
    assert read(tmpdir.join('Icalc.txt'))[1][4] == ''


# ------------------------------------------------------------- intensity

def test_icalc(tmpdir):
    """exp(-Eup/T_cm + d_off)*gA/rwn, with everything in cm-1."""
    out, done = extract(tmpdir, [row()], '--T_eV', '1.5', '--d_off', '2.0')
    assert out[0][6] == icalc.INTENSITY
    t_cm = 1.5 * fundconst.EV2CM
    ga = 14 * 1.29676e6
    want = math.exp(-138342.7 / t_cm + 2.0) * ga / 124990.6
    assert float(out[1][6]) == pytest.approx(want, rel=1e-5)
    assert 'Icalc at T = 1.5 eV' in done.stdout


def test_offset_is_a_log_offset(tmpdir):
    """d_off enters the exponent: +1 multiplies the column by e."""
    plain, _ = extract(tmpdir, [row()], '--T_eV', '1', '--d_off', '0')
    raised, _ = extract(tmpdir, [row()], '--T_eV', '1', '--d_off', '1')
    assert (float(raised[1][6]) / float(plain[1][6])
            == pytest.approx(math.e, rel=1e-5))


def test_temperature_is_in_ev(tmpdir):
    """T_eV is converted with EV2CM and not used as a wavenumber."""
    out, _ = extract(tmpdir, [row()], '--T_eV', '1', '--d_off', '0')
    want = math.exp(-138342.7 / fundconst.EV2CM) * 14 * 1.29676e6 / 124990.6
    assert float(out[1][6]) == pytest.approx(want, rel=1e-5)


# ------------------------------------------------------------ complaints

@pytest.mark.parametrize('arguments', [
    ('--T_eV', '1.0'),                  # one of the pair without the other
    ('--d_off', '0.0'),
    ('--T_eV', '0'),                    # not a temperature
    ('--T_eV', 'warm', '--d_off', '0'),
    ('--nonsense', '1'),
    ('a.txt', 'b.txt', 'c.txt'),
])
def test_usage(tmpdir, arguments):
    write(tmpdir.join('tp_E1.txt'), [HEADER, row()])
    done = run(tmpdir, *arguments)
    assert done.returncode == 2
    assert 'extract_Icalc:' in done.stderr
    assert not os.path.exists(str(tmpdir.join('Icalc.txt')))


def test_no_identifiers_in_the_source(tmpdir):
    """A line list insert_level_ids has not been run on says so."""
    header = [name for name in HEADER if name not in ('id1', 'id2')]
    write(tmpdir.join('tp_E1.txt'), [header, row()[:11] + row()[13:]])
    done = run(tmpdir)
    assert done.returncode == 1
    assert 'insert_level_ids' in done.stderr


def test_missing_column(tmpdir):
    """A column this needs and cannot do without is named."""
    header = [name for name in HEADER if name != 'J2']
    write(tmpdir.join('tp_E1.txt'), [header, row()[:6] + row()[7:]])
    done = run(tmpdir)
    assert done.returncode == 1
    assert 'J2' in done.stderr


def test_help():
    done = subprocess.run(
        [sys.executable, os.path.join(CODE, 'extract_Icalc.py'), '--help'],
        capture_output=True, text=True)
    assert done.returncode == 0
    assert 'extract_Icalc' in done.stdout and 'Icalc' in done.stdout
