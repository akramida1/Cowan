"""Reading the transition sections of OUTG11.

The columns are pinned against the FORMATs that write them -- the Lund
layout's `rcg11kd.f:1057` and `:1058`, and the wide layout's `:55` and
`:58` -- and against the reference runs, whose printed lines are counted
independently of the reader.

Every deck in WORK/ is ISPECC = 2, so the reference runs cover only the
Lund layout.  The wide one is exercised here on an excerpt copied
verbatim out of a WORK/Forbid run made with ISPECC = 4, and end to end
in `test_calc_tp.test_a_wavelength_sorted_run_gives_the_same_table`,
which builds that run and compares the two layouts line for line.
"""

import pytest

from conftest import case
from toolbox import spectrum
from toolbox.spectrum import SpectrumError


def read(name, **kw):
    return list(spectrum.lines(folder=case(name), **kw))


#: The head of the MAG DIP section of WORK/Forbid built with ISPECC = 4,
#: verbatim: the FORMAT 55 heading and its table of configurations, the
#: FORMAT 56 column header, and two FORMAT 58 lines.
WIDE = [
    '1          MAG DIP SPECTRUM     (ENERGIES IN UNITS OF  1000.0 CM-1'
    ' =    0.12 EV)',
    '',
    '    1      Sn7+  4d7            ---      Sn7+  4d7         ',
    '',
    '0           E      J   CONF              EP      JP  CONFP        '
    '   DELTA E   LAMBDA(A) S/PMAX**2   GF    LOG GF GA(SEC-1) CF,BRNCH',
    '',
    '    1   -10.4561  1.5  1 (4P) 4P       -8.6736  2.5  1 (4P) 4P    '
    '    1.7825  56102.1209   1.548   0.0000  -6.952 2.3643E-01  0.6970'
    '71     2     2',
    '    2   -19.8826  2.5  1 (4F) 4F      -18.0005  1.5  1 (4F) 4F    '
    '    1.8821  53132.5546   8.357   0.0000  -6.197 1.5027E+00  0.9883'
    '26     1     1',
]


def wide(tmp_path, body=WIDE):
    (tmp_path / 'OUTG11').write_text('\n'.join(body) + '\n', newline='')
    return list(spectrum.lines(folder=tmp_path))


# --------------------------------------------------------------------------
# the fixed columns
# --------------------------------------------------------------------------

def test_the_initial_level_line_reads_in_its_columns():
    first = read('W70')[0].initial
    # 5X,7H* * *  ,F10.3,F5.1,I4,2H (,A8,I6,1X,3A6
    assert (first.energy, first.j, first.config) == (2117.661, 0.5, 5)
    assert first.term == '3P) 4P'
    assert first.level == 1
    assert (first.spectrum, first.label) == ('W69+', 'sp2')


def test_the_final_level_line_reads_in_its_columns():
    one = read('W70')[0]
    # 1H ,I4,F11.3,F5.1,I3,1X,1H(,A8,F11.3,F11.3,F8.3,1PE10.3,0PF6.2,I6
    assert (one.final.energy, one.final.j, one.final.config) == (0.0, 0.5, 1)
    assert one.final.term == '1S) 2P' and one.final.level == 1
    assert (one.delta_e, one.wavelength) == (2117.661, 47.222)
    assert (one.log_gf, one.ga, one.cf) == (-1.744, 5.391e10, -0.94)


def test_the_multipole_says_whether_the_parity_changes():
    one = read('W70')[0]
    assert one.multipole == 'ELEC DIP' and not one.same_parity
    # SPID, rcg11kd.f:12812.  The other two connect a set to itself.
    assert set(spectrum.SAME_PARITY) == {'MAG DIP', 'ELEC QUD'}


# --------------------------------------------------------------------------
# what the reference runs hold
# --------------------------------------------------------------------------

@pytest.mark.parametrize('name,count', [('U90', 66), ('W70', 2205)])
def test_every_printed_line_is_read(name, count):
    # Counted independently: a line of the spectrum is the only thing in
    # OUTG11 with a level number in columns 81-86 and ' (' in 25-26.
    assert len(read(name)) == count


def test_a_pair_of_levels_is_printed_once():
    # Which is why the caller may treat a second appearance of a pair as
    # the other multipole rather than the reverse direction.
    seen = set()
    for one in read('W70'):
        key = frozenset([(one.initial.config, one.initial.j,
                          one.initial.level),
                         (one.final.config, one.final.j, one.final.level)])
        assert key not in seen
        seen.add(key)


@pytest.mark.slow
def test_a_serial_number_past_9999_does_not_lose_its_line():
    # Pr III prints 120274 lines, so the I4 counter of FORMAT 1058
    # overflows to asterisks for all but the first 9999 of them.
    assert len(read('Pr3')) == 120274


def test_a_coupling_filter_skips_the_sections_it_does_not_want():
    assert read('W70', multipoles=['MAG DIP']) == []


# --------------------------------------------------------------------------
# the wide layout, the one the wavelength sort prints
# --------------------------------------------------------------------------

def test_the_excerpt_is_the_layout_the_reader_expects():
    # If this fails the constants have drifted from the excerpt, and
    # every other wide-layout test below is passing for the wrong
    # reason.
    assert WIDE[4] == spectrum._WIDE_COLUMNS
    assert spectrum._configurations(WIDE[2]) == \
        (1, 'Sn7+  4d7         ', 'Sn7+  4d7         ')
    assert all(len(one) == 146 for one in WIDE[6:])


def test_a_wide_line_carries_both_levels_in_its_columns(tmp_path):
    one = wide(tmp_path)[0]
    # 1H ,I4,F11.4,F5.1,I3,2H (,A7,F13.4,F5.1,I3,2H (,A7,F13.4,F12.4,
    # F8.3,F9.4,F8.3,1PE11.4,0PF10.6,I6,I6
    assert (one.initial.energy, one.initial.j, one.initial.config) == \
        (-10.4561, 1.5, 1)
    assert one.initial.term == '4P) 4P' and one.initial.level == 2
    assert (one.final.energy, one.final.j, one.final.config) == \
        (-8.6736, 2.5, 1)
    assert one.final.term == '4P) 4P' and one.final.level == 2
    assert (one.delta_e, one.wavelength) == (1.7825, 56102.1209)
    assert (one.log_gf, one.ga, one.cf) == (-6.952, 2.3643e-01, 0.697071)
    assert one.multipole == 'MAG DIP' and one.same_parity


def test_a_wide_line_takes_both_names_from_the_heading(tmp_path):
    # It prints none beside either level, so they come from the FORMAT
    # 55 table, which numbers the configurations of each parity set the
    # way the CONF and CONFP columns index them.
    one = wide(tmp_path)[0]
    assert (one.initial.spectrum, one.initial.label) == ('Sn7+', '4d7')
    assert (one.final.spectrum, one.final.label) == ('Sn7+', '4d7')
    assert spectrum.name(folder=tmp_path) == 'Sn7+'


def test_every_wide_line_is_read(tmp_path):
    assert len(wide(tmp_path)) == 2


def test_a_layout_that_is_neither_is_refused_by_name(tmp_path):
    # FORMAT 57, the photoionization cross-section form: it prints
    # SIGMA(CM2) where FORMAT 56 prints LOG GF and GA, and no level
    # numbers at all.
    sigma = list(WIDE)
    sigma[4] = sigma[4].replace(' GA(SEC-1) CF,BRNCH', ' SIGMA(CM2)   CF')
    with pytest.raises(SpectrumError) as problem:
        wide(tmp_path, sigma)
    assert 'does not read' in str(problem.value)
    assert 'SIGMA(CM2)' in str(problem.value)


def test_a_file_with_no_spectrum_says_so():
    with pytest.raises(SpectrumError) as problem:
        list(spectrum.lines(path='genealogy.txt', folder=case('W70')))
    assert 'no transition section' in str(problem.value)
