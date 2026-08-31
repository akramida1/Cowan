"""The two files IDEN2 reads, and the search that fills them.

There is no reference `dlv.dat` in this repository -- an observed line
list belongs to a measurement, not to a calculation -- so the tests
build one out of two lines the run itself predicts and check that they
come back assigned, with no residual.  The formats are the ones
`conv_out.pl:PrintENLEV` and `:WriteTRANS_DAT` write, which is what
IDEN2 parses.
"""

import math
import os
import shutil
import subprocess
import sys

import pytest

from conftest import case
from toolbox import iden
from toolbox.iden import IdenError

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


@pytest.fixture
def run_folder(tmp_path):
    """A copy of WORK/U90 -- the IDEN2 files are written beside the deck."""
    folder = tmp_path / 'U90'
    shutil.copytree(str(case('U90')), str(folder),
                    ignore=shutil.ignore_patterns('reordered'))
    return folder


def observed(folder, wavenumbers, uncertainty=1e-5):
    """A dlv.dat of the given wavenumbers, and a numset.dat over them."""
    with open(str(folder / 'dlv.dat'), 'w') as out:
        for at, one in enumerate(wavenumbers):
            out.write('%5d%14.3f%34s%7.5f\n'
                      % (100 - 10 * at, one, '', uncertainty))
    with open(str(folder / 'numset.dat'), 'w') as out:
        out.write('%-5s%10.4f%7.5f\n' % ('a', 0.1, uncertainty))


def conv_out(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'conv_out.py'), 'lev.txt',
         'tr.txt', 'RCE', 'LTE', '1e6', '1000'] + list(arguments),
        cwd=str(folder), capture_output=True, text=True)


#: Two lines of the U90 run, in cm-1, between levels that both carry an
#: observed energy.
LINES = (810407104.0, 774051580.0)


def test_a_predicted_line_finds_its_observation(run_folder):
    observed(run_folder, LINES)
    out = conv_out(run_folder, 'IDEN')
    assert out.returncode == 0, out.stderr
    with open(str(run_folder / 'TRANS.DAT')) as f:
        rows = [line.rstrip('\n') for line in f if line.startswith('+')]
    # +%4d%5d%12.3f%2s%13.3f%6d%12.3f%11.3f%6d.  U90's wavenumbers are
    # 8e8 cm-1, too wide for F12.3, so the writer drops a decimal rather
    # than the number and the columns stay where IDEN2 expects them.
    assert {len(line) for line in rows} == {72}
    assigned = [line for line in rows if int(line[66:72]) > 0]
    assert len(assigned) == 2
    for line in assigned:
        assert float(line[55:66]) == 0.0            # no residual
        assert float(line[43:55]) in LINES          # at the observed one
        assert any('%.2f' % one in line for one in LINES)


def test_the_level_file_is_written_from_the_top_down(run_folder):
    observed(run_folder, LINES)
    out = conv_out(run_folder, 'IDEN')
    assert out.returncode == 0, out.stderr
    with open(str(run_folder / 'ENLEV.DAT')) as f:
        rows = f.readlines()
    assert len(rows) == 19
    # %4d%12.3f%10.3f%12.3f%2s%12.3f%5.1f /%11s/.  U90's energies are
    # 8e8 cm-1, so the three F12.3 fields go out as F12.2; the record
    # keeps its width and every field stays in its own columns.
    rows = [line.rstrip('\n') for line in rows]
    assert {len(line) for line in rows} == {71}
    energies = [float(line[4:16]) for line in rows]
    assert energies == sorted(energies, reverse=True)
    first = rows[0]
    assert first[58] == '/' and first.endswith('/')
    assert first[38:40] == ' *'                     # this level is observed
    assert first[16:26] == '    10.000'             # so its weight is 10


def test_truncation_keeps_only_the_levels_a_line_touches(run_folder):
    observed(run_folder, LINES)
    out = conv_out(run_folder, 'IDEN', 'truncate')
    assert out.returncode == 0, out.stderr
    with open(str(run_folder / 'ENLEV.DAT')) as f:
        rows = f.readlines()
    assert 0 < len(rows) < 19
    assert [int(line[:4]) for line in rows] == list(range(1, len(rows) + 1))


def test_convert2log_rewrites_the_observed_list_and_keeps_the_original(
        run_folder):
    observed(run_folder, LINES)
    out = conv_out(run_folder, 'IDEN', 'convert2log')
    assert out.returncode == 0, out.stderr
    assert (run_folder / 'dlv.dat.bak').exists()
    with open(str(run_folder / 'dlv.dat')) as f:
        first = f.readline()
    assert int(first[:5]) == 46                     # 10 * ln(100)


def test_a_missing_observed_list_is_named(run_folder):
    out = conv_out(run_folder, 'IDEN')
    assert out.returncode != 0
    assert 'dlv.dat' in out.stderr


# --------------------------------------------------------------------------
# the readers on their own
# --------------------------------------------------------------------------

def test_a_dlv_row_reads_in_its_columns(tmp_path):
    path = tmp_path / 'dlv.dat'
    path.write_text('%5d%14.3f%34s%7.4f  tail\n' % (25, 486071.61, 'x' * 34,
                                                    0.007))
    one = iden.read_dlv(str(path))[0]
    assert (one.intensity, one.wavenumber, one.uncertainty) \
        == (25.0, 486071.61, 0.007)
    assert one.middle == 'x' * 34


def test_a_numset_row_needs_both_of_its_numbers(tmp_path):
    path = tmp_path / 'numset.dat'
    path.write_text('%-5s%10.4f%7.4f\n' % ('uv', 2000.0, 0.0))
    with pytest.raises(IdenError) as problem:
        iden.read_numset(str(path))
    assert 'default uncertainty' in str(problem.value)


def test_the_search_runs_downward_in_wavenumber():
    lines = [iden.Observed(intensity=1.0, wavenumber=w, uncertainty=0.1,
                           middle='', tail='')
             for w in (50000.0, 40000.0, 30000.0, 20000.0)]
    # 0.1 A at 30000 cm-1 is 0.9 cm-1, so only a near miss is accepted.
    assert iden.find(lines, 30000.5, 1, 0.1) == 3
    assert iden.find(lines, 35000.0, 1, 0.1) == 0
    assert iden.find(lines, 50000.0, 1, 0.1) == 1


def _skewed(count=300):
    """A linear observed list of the shape a real one has.

    Line intensities are steeply skewed -- a few strong lines and a long
    tail of weak ones -- and a measured list has a floor: the weakest
    line anyone reported, not an intensity of 1.  Both matter to
    `iden.looks_logarithmic`, and these numbers give the five decades
    and the floor near 20 that `WORK/Pr3/dlv.dat` has.
    """
    return [iden.Observed(intensity=2e6 / (at * at),
                          wavenumber=50000.0 - at, uncertainty=0.01,
                          middle='', tail='')
            for at in range(1, count + 1)]


def test_an_already_logarithmic_list_is_told_from_a_linear_one():
    linear = _skewed()
    assert not iden.looks_logarithmic(linear)
    logged = [iden.Observed(intensity=10 * math.log(one.intensity),
                            wavenumber=one.wavenumber,
                            uncertainty=one.uncertainty, middle='', tail='')
              for one in linear]
    assert iden.looks_logarithmic(logged)


def test_too_few_lines_to_judge_the_scale_are_left_alone():
    assert not iden.looks_logarithmic(_skewed(iden.LOG_SAMPLE - 1))


def test_convert2log_warns_when_the_list_is_already_logarithmic(run_folder):
    lines = _skewed()
    with open(str(run_folder / 'dlv.dat'), 'w') as out:
        for one in lines:
            out.write('%5d%14.3f%34s%7.5f\n'
                      % (10 * math.log(one.intensity), one.wavenumber, '',
                         1e-5))
    with open(str(run_folder / 'numset.dat'), 'w') as out:
        out.write('%-5s%10.4f%7.5f\n' % ('a', 0.1, 1e-5))
    out = conv_out(run_folder, 'IDEN', 'convert2log')
    assert 'already look logarithmic' in out.stdout


def test_a_log_intensity_below_one_is_capped_at_zero(tmp_path):
    path = tmp_path / 'dlv.dat'
    with open(str(path), 'w') as out:
        for intensity, wavenumber in ((100, 5000.0), (1, 4000.0),
                                      (0, 3000.0)):
            out.write('%5d%14.3f%34s%7.4f\n' % (intensity, wavenumber, '',
                                                0.01))
    lines = iden.read_dlv(str(path))
    lines[1].intensity = 0.4        # weaker than the I5 field can carry
    iden.write_dlv(lines, str(path))
    # 10*ln(0.4) is negative and read_dlv takes the field with abs(), so
    # an uncapped log would read back as a line of intensity 9.
    assert [one.intensity for one in iden.read_dlv(str(path))] \
        == [46.0, 0.0, 0.0]


def test_a_list_with_no_intensity_at_all_cannot_go_on_a_log_scale(tmp_path):
    path = tmp_path / 'dlv.dat'
    path.write_text('%5d%14.3f%34s%7.4f\n' % (0, 5000.0, '', 0.01))
    with pytest.raises(IdenError) as problem:
        iden.write_dlv(iden.read_dlv(str(path)), str(path))
    assert 'every intensity' in str(problem.value)
