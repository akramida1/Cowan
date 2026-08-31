"""The parameter-file reader, and the file calc_TP ships."""

import pytest

import calc_TP
from toolbox import parfile
from toolbox.parfile import ParFileError, Spec

SPECS = (
    Spec('count', parfile.INT, 7, minimum=0, maximum=99,
         description='How many.'),
    Spec('width', parfile.FLOAT, 0.5, minimum=0.0, description='How wide.'),
    Spec('loud', parfile.BOOL, False, description='Whether to shout.'),
    Spec('name', parfile.TEXT, 'out.txt', description='Where it goes.'),
    Spec('mode', parfile.TEXT, 'one', choices=('one', 'both'),
         description='Which.'),
)


def written(tmp_path, text, name='x.par'):
    path = tmp_path / name
    path.write_text(text)
    return path


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def test_a_value_is_everything_left_of_the_semicolon(tmp_path):
    got = parfile.read(written(tmp_path,
                               '  12   ; count  how many of them\n'), SPECS)
    assert got['count'] == 12


def test_a_key_left_out_takes_its_default(tmp_path):
    got = parfile.read(written(tmp_path, '3 ; count\n'), SPECS)
    assert got == {'count': 3, 'width': 0.5, 'loud': False,
                   'name': 'out.txt', 'mode': 'one'}


def test_a_blank_value_takes_the_default_too(tmp_path):
    # This is how a parameter is turned off without deleting the line
    # that documents it.
    got = parfile.read(written(tmp_path, '     ; count  nothing here\n'),
                       SPECS)
    assert got['count'] == 7


def test_comments_and_blank_lines_are_ignored(tmp_path):
    text = ('; a whole-line comment\n'
            '\n'
            '   ; an indented comment\n'
            '4 ; count\n')
    assert parfile.read(written(tmp_path, text), SPECS)['count'] == 4


def test_order_does_not_matter(tmp_path):
    a = parfile.read(written(tmp_path, '4 ; count\n1 ; loud\n'), SPECS)
    b = parfile.read(written(tmp_path, '1 ; loud\n4 ; count\n'), SPECS)
    assert a == b


def test_both_spellings_of_a_flag(tmp_path):
    for text in ('1', 'yes', 'TRUE', 'on'):
        assert parfile.read(written(tmp_path, '%s ; loud\n' % text),
                            SPECS)['loud'] is True
    for text in ('0', 'no', 'False', 'off'):
        assert parfile.read(written(tmp_path, '%s ; loud\n' % text),
                            SPECS)['loud'] is False


def test_text_keeps_its_spaces_inside(tmp_path):
    got = parfile.read(written(tmp_path, ' my file.txt ; name\n'), SPECS)
    assert got['name'] == 'my file.txt'


# --------------------------------------------------------------------------
# what it refuses
# --------------------------------------------------------------------------

@pytest.mark.parametrize('text, says', [
    ('4 ; counts  a typo\n', "'counts' is not a parameter"),
    ('4 ; count\n5 ; count\n', 'already given on line 1'),
    ('4\n', 'a value with no key'),
    ('x ; count\n', 'must be a whole number'),
    ('4.5 ; count\n', 'must be a whole number'),
    ('-1 ; count\n', 'must be at least 0'),
    ('100 ; count\n', 'must be at most 99'),
    ('maybe ; loud\n', 'must be 0 or 1'),
    ('neither ; mode\n', 'must be one of one, both'),
])
def test_a_bad_line_is_named_by_its_number(tmp_path, text, says):
    with pytest.raises(ParFileError) as problem:
        parfile.read(written(tmp_path, text), SPECS)
    assert says in str(problem.value)


# --------------------------------------------------------------------------
# writing
# --------------------------------------------------------------------------

def test_a_written_file_reads_back_as_the_values_it_was_given(tmp_path):
    values = {'count': 12, 'width': 1.25, 'loud': True, 'name': 'a.txt',
              'mode': 'both'}
    path = parfile.write(tmp_path / 'x.par', SPECS, values, title=('Hi',))
    assert parfile.read(path, SPECS) == values


def test_a_written_file_reads_back_as_the_defaults(tmp_path):
    path = parfile.write(tmp_path / 'x.par', SPECS)
    assert parfile.read(path, SPECS) == {s.key: s.default for s in SPECS}


def test_an_existing_file_is_never_replaced(tmp_path):
    parfile.write(tmp_path / 'x.par', SPECS)
    with pytest.raises(ParFileError):
        parfile.write(tmp_path / 'x.par', SPECS)


def test_every_description_line_of_a_written_file_is_a_comment(tmp_path):
    path = parfile.write(tmp_path / 'x.par', SPECS, title=('Hi', 'there'))
    for line in path.read_text().split('\n'):
        if line.strip() and not line.lstrip().startswith(';'):
            assert ';' in line, line          # a value line carries its key


# --------------------------------------------------------------------------
# the file the program ships
# --------------------------------------------------------------------------

def test_the_shipped_par_file_is_the_defaults():
    # CODE/calc_TP.par is what a run falls back on when the working
    # directory has none, so it has to say exactly what SPECS says.
    from toolbox import files
    path = files.code_dir() / calc_TP.PAR_FILE
    assert path.is_file()
    got = parfile.read(path, calc_TP.SPECS)
    assert got == {s.key: s.default for s in calc_TP.SPECS}


def test_the_shipped_par_file_is_what_write_par_would_write():
    from toolbox import files
    path = files.code_dir() / calc_TP.PAR_FILE
    fresh = parfile.template(calc_TP.SPECS, title=calc_TP.TITLE)
    assert path.read_text().split('\n')[:len(fresh)] == fresh


def test_every_parameter_of_the_program_is_documented():
    for spec in calc_TP.SPECS:
        assert spec.description.strip(), spec.key


# --------------------------------------------------------------------------
# finding the file, and a parameter with no default
# --------------------------------------------------------------------------
#
# `find` is one copy of a rule two programs need -- the named file, the
# working directory, then the copy beside the code -- so that phase 11
# does not leave two idioms behind.

def test_a_file_named_outright_is_the_one_used(tmp_path):
    named = written(tmp_path, '4 ; count\n', name='elsewhere.par')
    assert parfile.find('nothing.par', tmp_path, str(named)) == named


def test_a_file_named_outright_that_is_not_there_is_an_error(tmp_path):
    with pytest.raises(ParFileError):
        parfile.find('nothing.par', tmp_path, str(tmp_path / 'gone.par'))


def test_the_working_directory_comes_before_the_shipped_copy(tmp_path):
    here = written(tmp_path, '4 ; count\n', name='calc_TP.par')
    assert parfile.find('calc_TP.par', tmp_path) == here


def test_the_copy_beside_the_code_is_the_fallback(tmp_path):
    from toolbox import files
    found = parfile.find('calc_TP.par', tmp_path)
    assert found == files.code_dir() / 'calc_TP.par'


def test_a_program_that_ships_nothing_gets_nothing(tmp_path):
    # conv_out must not ship a par file, or the file would never be
    # absent and the sample would never be written.
    assert parfile.find('calc_TP.par', tmp_path, shipped=False) is None
    assert parfile.find('conv_out.par', tmp_path) is None


NEEDED = (Spec('who', parfile.TEXT, None, description='Who it is for.'),
          Spec('many', parfile.INT, None, description='How many.'),
          Spec('wide', parfile.FLOAT, None, minimum=0.0,
               description='How wide.'))


def test_a_parameter_with_no_default_is_written_blank(tmp_path):
    lines = parfile.template(NEEDED)
    for spec in NEEDED:
        assert ('%s; %s' % (' ' * 11, spec.key)) in lines


def test_a_parameter_with_no_default_reads_back_as_none(tmp_path):
    path = tmp_path / 'needed.par'
    parfile.write(path, NEEDED)
    assert parfile.read(path, NEEDED) == {'who': None, 'many': None,
                                          'wide': None}


def test_a_parameter_with_no_default_is_still_range_checked(tmp_path):
    with pytest.raises(ParFileError):
        parfile.read(written(tmp_path, '-1 ; wide\n'), NEEDED)
