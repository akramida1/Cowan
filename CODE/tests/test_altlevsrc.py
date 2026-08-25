"""subst_exp_levs.py: an experimental level list that is not an RCE deck.

`WORK/Pr3/Pr3_exp_levs.txt` is the reference case, and it is a hard one:
593 levels designated in a mixture of LS, JK and pair coupling, with the
`perc2/conf2/term2` triad meaning *the leading component written again
in another scheme* on the rows whose `remark` column says `or`.  It is
also the first real file this branch of the program has ever had -- the
LOPT reader it replaces was ported from the Perl line for line and never
run against data.

The test that carries the most weight is the last one:
`--matchby Eexp` against `RCEINP.fin` reproduces the deck byte for byte,
because that file's observed energies are the ones the sheet holds.  A
matcher that put one level's energy on another would not do that.
"""

import os
import shutil
import subprocess
import sys

import pytest

from conftest import WORK
from toolbox import altlevsrc, explevs

CODE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

SOURCE = WORK / 'Pr3' / 'Pr3_exp_levs.txt'
TARGET = WORK / 'Pr3' / 'RCEINP.fin'

#: The columns of `Pr3_exp_levs.txt`, as the author's sheet writes them.
PR3 = ['--ecol', 'energy', '--jcol', 'j_val', '--pcol', 'parity',
       '--confcol', 'conf', '--termcol', 'term', '--perccol', 'perc',
       '--term1col', 'term1',
       '--conf2col', 'conf2', '--term2col', 'term2', '--perc2col', 'perc2',
       '--altcol', 'remark', '--altvalue', 'or']


def run(folder, *arguments):
    return subprocess.run(
        [sys.executable, os.path.join(CODE, 'subst_exp_levs.py')]
        + list(arguments), cwd=str(folder), capture_output=True, text=True)


def pr3():
    if not (SOURCE.is_file() and TARGET.is_file()):
        pytest.skip('WORK/Pr3 has no Pr3_exp_levs.txt and RCEINP.fin')
    return SOURCE, TARGET


def prepare(tmp_path):
    source, target = pr3()
    tmp_path.mkdir(parents=True, exist_ok=True)
    shutil.copy(str(source), str(tmp_path / 'SRC'))
    shutil.copy(str(target), str(tmp_path / 'TGT'))
    return tmp_path


def write(path, rows):
    path.write_text('\n'.join('\t'.join(row) for row in rows) + '\n')
    return str(path)


# --------------------------------------------------------------------------
# labels
# --------------------------------------------------------------------------

def test_the_term_label_says_which_coupling_it_is_written_in():
    assert altlevsrc.coupling('2P*') == 'ls'
    assert altlevsrc.coupling('4D') == 'ls'
    assert altlevsrc.coupling('2[1]*') == 'jk'
    assert altlevsrc.coupling('(2,3/2)*') == 'jj'
    assert altlevsrc.coupling('*') is None
    assert altlevsrc.coupling('  ') is None


def test_the_term_key_drops_the_parity_mark_and_the_repeat_index():
    assert altlevsrc.term_key('2P*') == '2P'
    assert altlevsrc.term_key('(3F) 4D3') == '4D'
    assert altlevsrc.term_key('(2P) 2P') == '2P'
    assert altlevsrc.term_key('') == ''


def test_a_configuration_is_compacted_the_way_rce_abbreviates_it():
    assert altlevsrc.compact('4f2.(3F<2>).6p') == '4f26p'
    assert altlevsrc.compact('4f.(2F*).5d2.(3F)') == '4f5d2'
    assert altlevsrc.compact('4f3') == '4f3'


@pytest.mark.parametrize('written, label', [
    ('4f2.(3F<2>).6p', 'f26p'),
    ('4f.(2F*).5d2.(3F)', 'f5d2'),
    ('4f3', '4f3'),
    ('4f.5d.(3F*).6s', 'f5d6s'),
    ('5d.6s.6p', 'd6s6p'),
    ('4f2.(3H).5f', 'f25f'),
])
def test_the_deck_label_a_written_out_configuration_abbreviates_to(
        written, label):
    """The deck's label is a tail of the full form; the longest wins."""
    labels = {'4f3', 'f26p', 'f27p', 'f25f', 'f5d2', 'f5d6s', 'd6s6p',
              'f6p2', '6p3'}
    assert altlevsrc.resolve_configuration(written, labels) == label


def test_a_configuration_the_deck_does_not_have_is_not_guessed_at():
    assert altlevsrc.resolve_configuration('4f2.7d', {'4f3', 'f26p'}) == ''


def test_the_parity_of_an_abbreviated_label():
    """The digit after a subshell letter is the *next* subshell's n
    whenever another follows: `f25d` is 4f2.5d, not 4f25.d."""
    assert altlevsrc.parity_of_label('4f3') == 1
    assert altlevsrc.parity_of_label('f5d2') == 1
    assert altlevsrc.parity_of_label('d6s6p') == 1
    assert altlevsrc.parity_of_label('f26p') == 1
    assert altlevsrc.parity_of_label('f25d') == 0
    assert altlevsrc.parity_of_label('f6s2') == 1


def test_a_term_with_no_j_stands_for_every_j_it_has():
    """A row that names a term rather than a level."""
    assert altlevsrc._js('', '3P') == ['0', '1', '2']
    assert altlevsrc._js('', '4F') == ['3/2', '5/2', '7/2', '9/2']
    assert altlevsrc._js('5/2', '2D') == ['5/2']


# --------------------------------------------------------------------------
# reading, without a deck to match against
# --------------------------------------------------------------------------

class FakeDeck:
    """A target deck of two levels, enough to key and file a source row."""

    def __init__(self):
        odd = explevs.Level(j='1/2', ee=0.0, ec=1.0, observed=False,
                            vector={'4f3,(2P) 2P': 96}, kset=1)
        even = explevs.Level(j='1/2', ee=0.0, ec=1.0, observed=False,
                             vector={'f25d,(3H) 4I': 96}, kset=2)
        self.sets = {1: {'1/2': [odd]}, 2: {'1/2': [even]}}

    def levels(self):
        for kset in sorted(self.sets):
            for j in self.sets[kset]:
                for level in self.sets[kset][j]:
                    yield level

    def group(self, kset, j):
        return self.sets.get(kset, {}).get(j, [])


def read(path, **options):
    return altlevsrc.read_alt_lev_src(path, altlevsrc.Spec(**options),
                                      FakeDeck())


def test_the_parity_sets_are_read_off_the_target_not_off_the_source():
    assert altlevsrc.parity_sets(FakeDeck()) == {'o': 1, 'e': 2}


def test_a_header_line_is_found_without_being_told_about(tmp_path):
    path = write(tmp_path / 'src.txt',
                 [['energy', 'j_val', 'parity'],
                  ['18693.65', '1/2', 'o']])
    deck, notes = read(path, ecol='energy', jcol='j_val', pcol='parity')
    assert notes.skipped == 1 and notes.levels == 1
    assert list(deck.levels())[0].ee == pytest.approx(18.69365)


def test_a_file_with_no_header_needs_column_numbers_and_nothing_else(
        tmp_path):
    path = write(tmp_path / 'src.txt', [['18693.65', '1/2', 'o']])
    deck, notes = read(path, ecol=1, jcol=2, pcol=3)
    assert notes.skipped == 0 and notes.levels == 1


@pytest.mark.parametrize('units, energy, kk', [
    ('cm-1', '18693.65', 18.69365),
    ('kK', '18.69365', 18.69365),
    ('eV', '1.0', 8.06554393734921),
])
def test_the_energies_are_brought_into_the_deck_s_own_unit(
        units, energy, kk, tmp_path):
    path = write(tmp_path / 'src.txt', [[energy, '1/2', 'o']])
    deck, _ = read(path, ecol=1, jcol=2, pcol=3, units=units)
    assert list(deck.levels())[0].ee == pytest.approx(kk)


def test_the_exclude_column_is_found_by_its_name(tmp_path):
    path = write(tmp_path / 'src.txt',
                 [['energy', 'j_val', 'parity', 'exclude'],
                  ['18693.65', '1/2', 'o', ''],
                  ['23465.43', '1/2', 'o', '0'],
                  ['66324.84', '1/2', 'o', 'x'],
                  ['71501.00', '1/2', 'o', '1']])
    _, notes = read(path, ecol='energy', jcol='j_val', pcol='parity')
    assert (notes.rows, notes.levels, notes.excluded) == (4, 2, 2)


def test_a_column_name_that_is_in_no_header_is_an_error(tmp_path):
    path = write(tmp_path / 'src.txt', [['18693.65', '1/2', 'o']])
    with pytest.raises(altlevsrc.AltSourceError):
        read(path, ecol='energy', jcol=2, pcol=3)


def test_the_three_columns_that_have_no_default(tmp_path):
    path = write(tmp_path / 'src.txt', [['18693.65', '1/2', 'o']])
    with pytest.raises(altlevsrc.AltSourceError):
        read(path, ecol=1, jcol=2)


def test_labels_with_no_percentage_are_taken_as_the_whole_level(tmp_path):
    path = write(tmp_path / 'src.txt',
                 [['18693.65', '1/2', 'o', '4f3', '2P*']])
    deck, _ = read(path, ecol=1, jcol=2, pcol=3, confcol=4, termcol=5)
    assert list(deck.levels())[0].vector == {'4f3,2P': 100.0}


def test_conf1_and_term1_override_the_bookkeeping_labels(tmp_path):
    path = write(tmp_path / 'src.txt',
                 [['18693.65', '1/2', 'o', '4f3', '*', '95', '2P*']])
    deck, _ = read(path, ecol=1, jcol=2, pcol=3, confcol=4, termcol=5,
                   perccol=6, term1col=7)
    assert list(deck.levels())[0].vector == {'4f3,2P': 95.0}


def test_a_blank_configuration_in_a_later_triad_is_the_first_one(tmp_path):
    path = write(tmp_path / 'src.txt',
                 [['18693.65', '1/2', 'o', '4f3', '2P*', '95', '5', '',
                   '4D*']])
    deck, _ = read(path, ecol=1, jcol=2, pcol=3, confcol=4, termcol=5,
                   perccol=6, perc2col=7, conf2col=8, term2col=9)
    assert list(deck.levels())[0].vector == {'4f3,2P': 95.0, '4f3,4D': 5.0}


def test_the_alternate_coupling_flag_makes_one_component_of_two(tmp_path):
    """`remark = or`: the second triad is the first component again."""
    rows = [['111779.59', '1/2', 'o', '4f2.(3H<4>).5f', '2[1]*', '92', 'or',
             '92', '4f2.(3H).5f', '4D*']]
    path = write(tmp_path / 'src.txt', rows)
    deck, notes = read(path, ecol=1, jcol=2, pcol=3, confcol=4, termcol=5,
                       perccol=6, altcol=7, perc2col=8, conf2col=9,
                       term2col=10, altvalue='or')
    level = list(deck.levels())[0]
    assert notes.couplings == {'jk': 1}
    # One component, the LS label of it, and the deck knows no f25f.
    assert notes.unknown_conf and level.vector == {}


def test_a_leading_component_in_another_coupling_leaves_the_row_bare(
        tmp_path):
    path = write(tmp_path / 'src.txt',
                 [['66324.84', '1/2', 'o', '4f3', '(2,3/2)*', '81']])
    deck, notes = read(path, ecol=1, jcol=2, pcol=3, confcol=4, termcol=5,
                       perccol=6)
    assert list(deck.levels())[0].vector == {}
    assert (notes.by_vector, notes.by_energy) == (0, 1)
    assert notes.couplings == {'jj': 1}


# --------------------------------------------------------------------------
# matching
# --------------------------------------------------------------------------

def test_the_target_is_brought_into_intensity_per_cent():
    deck = FakeDeck()
    altlevsrc.normalize_target(deck)
    level = list(deck.levels())[0]
    assert level.vector == {'4f3,(2P) 2P': 96}
    assert level.mvector == {'4f3,2P': pytest.approx(92.16)}


def test_a_source_that_names_two_components_is_not_charged_for_ten():
    """`distance_named`, and why a partial vector needs it.

    The source says 24% of one basis state and 76% of something it does
    not name.  `distance` charges every component the target has and
    the source has not, so the target it likes best is whichever level
    is most thoroughly mixed -- here the one that shares nothing.
    """
    mine = {'f25f,2P': 24.0}
    right = {'f25f,2P': 30.0, 'f25f,4D': 20.0, 'f25f,2D': 15.0}
    wrong = {'f4p5,6F': 7.0, 'f4p5,2P2': 5.0, 'f4p5,4F': 5.0}
    # `distance` cannot tell the two apart: it charges the right level
    # for the components the source never claimed to list.
    assert (explevs.distance(mine, wrong)[0]
            / explevs.distance(mine, right)[0]) < 1.1
    assert (explevs.distance_named(mine, wrong)[0]
            / explevs.distance_named(mine, right)[0]) > 10


def test_two_vectors_agreeing_on_one_state_in_intensity_settle_the_match():
    assert explevs.distance_named({'a': 60.0}, {'a': 55.0})[1]
    assert not explevs.distance_named({'a': 40.0}, {'a': 45.0})[1]


# --------------------------------------------------------------------------
# the reference case
# --------------------------------------------------------------------------

def test_the_pr3_sheet_reads_as_the_author_describes_it(tmp_path):
    source, target = pr3()
    deck = explevs.read_deck(str(target))
    spec = altlevsrc.Spec(ecol='energy', jcol='j_val', pcol='parity',
                          confcol='conf', termcol='term', perccol='perc',
                          term1col='term1', conf2col='conf2',
                          term2col='term2', perc2col='perc2',
                          altcol='remark', altvalue='or')
    levels, notes = altlevsrc.read_alt_lev_src(str(source), spec, deck)
    assert notes.rows == 593 and notes.levels == 593
    assert notes.parities == {'o': 1, 'e': 2}
    # 343 rows expose an LS label for their leading component -- in
    # `term`, in `term1`, or in `term2` where `remark` says `or`; the
    # other 250 are designated in JK or pair coupling, or not at all.
    assert notes.by_vector == 343 and notes.by_energy == 250
    assert sum(notes.couplings.values()) == 593
    assert notes.couplings['ls'] == 226


def test_every_level_of_the_pr3_sheet_finds_a_home(tmp_path):
    mine = prepare(tmp_path)
    out = run(mine, 'SRC', 'TGT', *PR3)
    assert out.returncode == 0, out.stderr
    assert '593 observed energies' in out.stdout
    assert '0 source levels went unmatched' in out.stdout
    assert 'same observed energy' not in out.stdout


def test_matchby_eexp_puts_the_pr3_deck_back_exactly_as_it_was(tmp_path):
    """The one test with an oracle.

    `RCEINP.fin` already carries the observed energies the sheet holds,
    so replacing them with themselves must leave the file alone -- byte
    for byte, in the §1 sense.  Nothing but a matcher that tied every
    source level to the level it came from would manage it.
    """
    mine = prepare(tmp_path)
    out = run(mine, 'SRC', 'TGT', '--ecol', '1', '--jcol', '3',
              '--pcol', '4', '--matchby', 'Eexp')
    assert out.returncode == 0, out.stderr
    assert (mine / 'TGT').read_bytes() == TARGET.read_bytes()


# --------------------------------------------------------------------------
# the command line
# --------------------------------------------------------------------------

def test_no_arguments_and_help_print_the_syntax(tmp_path):
    for arguments in ([], ['--help']):
        out = run(tmp_path, *arguments)
        assert out.returncode == 0
        assert '--ecol' in out.stdout and '--matchby' in out.stdout


def test_the_uncertainty_options_are_gone(tmp_path):
    mine = prepare(tmp_path)
    assert run(mine, 'SRC', 'TGT', '--max-unc', '0.5').returncode == 2
    assert run(mine, 'SRC', 'TGT', '--unccol', '5').returncode == 2


def test_the_options_of_a_level_list_are_refused_for_an_rce_deck(tmp_path):
    mine = prepare(tmp_path)
    shutil.copy(str(TARGET), str(mine / 'SRC2'))
    out = run(mine, 'SRC2', 'TGT', '--ecol', '1')
    assert out.returncode == 2 and '--ecol' in out.stderr


def test_params_is_refused_for_a_level_list(tmp_path):
    mine = prepare(tmp_path)
    out = run(mine, 'SRC', 'TGT', 'params=none', '--ecol', '1',
              '--jcol', '3', '--pcol', '4')
    assert out.returncode == 2 and 'parameter table' in out.stderr


def test_a_source_that_is_not_a_deck_needs_its_three_columns(tmp_path):
    mine = prepare(tmp_path)
    out = run(mine, 'SRC', 'TGT')
    assert out.returncode == 1 and '--ecol' in out.stderr


@pytest.mark.parametrize('bad', [
    ['--matchby', 'sideways'],
    ['--units', 'furlongs'],
    ['--skip', 'two'],
    ['--sideways', '1'],
])
def test_the_arguments_are_checked(bad, tmp_path):
    mine = prepare(tmp_path)
    assert run(mine, 'SRC', 'TGT', '--ecol', '1', '--jcol', '3',
               '--pcol', '4', *bad).returncode == 2
