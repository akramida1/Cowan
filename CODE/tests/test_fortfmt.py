"""Fortran field editing, at and past the points where it switches form."""

from toolbox import fortfmt as ff


def test_blank_fields_are_zero():
    assert ff.read_i('     ') == 0
    assert ff.read_f('      ', 3) == 0.0


def test_implied_decimal_point():
    # f6.3 with no point in the field: the last three digits are the
    # fraction.  IN36 punches occupations this way.
    assert ff.read_f('  1500', 3) == 1.5
    assert ff.read_f('     7', 3) == 0.007
    assert ff.read_f('  1.50', 3) == 1.5
    assert ff.read_f('  1.50', 0) == 1.5


def test_exponents():
    assert ff.read_f('5.E-08', 3) == 5e-08
    assert ff.read_f('1.E-11', 3) == 1e-11
    assert ff.read_f(' 1.5-8', 3) == 1.5e-08      # bare-sign exponent
    assert ff.read_f('2.5D+3', 3) == 2500.0


def test_integer_output_overflows_to_asterisks():
    assert ff.write_i(12345, 9) == '    12345'
    assert ff.write_i(-1234567890123, 9) == '*********'


def test_write_param_switches_form():
    # Normal: scaled integer, 10000 x value.
    assert ff.write_param(1.2345) == '    12345'
    # Zero to within 1e-3 and not EAV: floating point.
    assert ff.write_param(0.0).strip().upper().startswith('0.00E')
    # EAV is exempt from that rule.
    assert ff.write_param(0.0, first=True) == '        0'
    # Too large for the scaled integer.
    assert ff.write_param(123456.0) == '123456.00'
    assert ff.write_param(123456.0, first=True) == ' 123456.0'


def test_half_integers_round_trip():
    assert ff.half_int(0) == '0'
    assert ff.half_int(1) == '1/2'
    assert ff.half_int(11) == '11/2'
    assert ff.half_int(4) == '2'
    for text in ('0', '1/2', '11/2', '2'):
        assert ff.half_int(ff.parse_half_int(text)) == text
    assert ff.parse_half_int(' 5.5 ') == 11


def test_a_wide_value_gives_up_decimals_before_it_gives_up_digits():
    """The rule the Fortran uses where a value may outgrow its field."""
    # Fits as written.
    assert ff.write_f_wide(1234.5678, 12, 3) == '    1234.568'
    # U90's energies in cm-1: F12.3 would overflow, F12.2 does not.
    assert ff.write_f_wide(810407104.0, 12, 3) == '810407104.00'
    assert ff.write_f_wide(8104071040.0, 12, 3) == '8104071040.0'
    assert ff.write_f_wide(81040710400.0, 12, 3) == ' 81040710400'
    # Past that the integer part will not fit either, so exponential
    # with the widest mantissa the field allows -- still a number.
    wide = ff.write_f_wide(8.10407104e14, 12, 3)
    assert len(wide) == 12 and abs(float(wide) / 8.10407104e14 - 1) < 1e-6
    # A sign costs a digit, and nothing overflows to asterisks.
    assert ff.write_f_wide(-810407104.0, 12, 3) == '-810407104.0'
