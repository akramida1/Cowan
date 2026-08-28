"""calc_TP: the merge, the parameter map, and one whole run.

The end-to-end test builds `WORK/Forbid` from its two input cards by
running RCN, RCN2, RCG and RCE into a temporary directory, then runs the
program over it.  It is the only reference case in this package with an
electric-quadrupole spectrum and a //R2// card in it; see that folder's
README for why it is inputs only and what it does not cover.
"""

import math
import os
import shutil
import subprocess

import pytest

import calc_TP
from calc_TP import Lev, Row, Transition
from conftest import CODE, EXECUTABLES, WORK, case
import insert_level_ids
from toolbox import ing11, levids, parfile, radiate, rce, vacair


# --------------------------------------------------------------------------
# building rows by hand
# --------------------------------------------------------------------------

def level(index, energy, j=2.0, kset=1, observed=False):
    one = rce.Level(e_obs=energy, e_calc=energy, j=j, kset=kset, index=index,
                    flag='' if observed else '*')
    from toolbox.levels import Identified
    return Lev(Identified(level=one, kset=kset, j=j, rcg=index))


def transition(multipole, ga, cf=0.5, upper=None, lower=None, values=()):
    """One transition with `values` as its per-trial A.

    The program tracks S, so the trial values go in as line strengths at
    the one wavelength this fixture works at; A comes back out of them
    unchanged, which is what the merge tests are written against.
    """
    one = Transition(upper, lower, multipole)
    one.ga, one.cf = ga, cf
    one.wavelength = 1000.0
    one.s = radiate.s_from_ga(ga, one.wavelength, multipole)
    weight = 2 * upper.j + 1
    for at, a in enumerate(values, 1):
        one.n_seen += 1
        one.log_s.add(math.log(radiate.s_from_ga(a * weight, one.wavelength,
                                                 multipole)))
        one.cancel.add(cf)
        one.a_values.append((at, a))
        one.cf_values.append((at, cf))
    return one


def rows(e2_ga, m1_ga):
    up, down = level(2, 100.0), level(1, 0.0)
    return [Row(transition(radiate.E2, e2_ga, upper=up, lower=down)),
            Row(transition(radiate.M1, m1_ga, upper=up, lower=down))]


# --------------------------------------------------------------------------
# the M1 + E2 merge
# --------------------------------------------------------------------------

def test_a_line_that_is_all_e2_keeps_the_name_e2():
    out = calc_TP.merge(rows(100.0, 1.0), 0, 0.02)
    assert len(out) == 1
    assert out[0].kind == radiate.E2
    assert out[0].ga == pytest.approx(101.0)
    assert out[0].e2_fraction == pytest.approx(100.0 / 101.0)


def test_a_line_that_is_all_m1_keeps_the_name_m1_and_its_own_strength():
    out = calc_TP.merge(rows(1.0, 100.0), 0, 0.02)
    assert len(out) == 1
    assert out[0].kind == radiate.M1
    assert out[0].ga == pytest.approx(100.0)
    assert out[0].e2_fraction == pytest.approx(1.0 / 101.0)


def test_a_mixed_line_is_called_m1_plus_e2_and_the_two_add():
    out = calc_TP.merge(rows(50.0, 50.0), 0, 0.02)
    assert len(out) == 1
    assert out[0].kind == 'M1+E2'
    assert out[0].ga == pytest.approx(100.0)
    assert out[0].e2_fraction == pytest.approx(0.5)
    # A merged line has no single line strength, so none is quoted.
    assert out[0].s is None


def test_a_multipole_with_no_partner_goes_through_untouched():
    up, down = level(2, 100.0), level(1, 0.0)
    only = [Row(transition(radiate.E2, 7.0, upper=up, lower=down))]
    out = calc_TP.merge(only, 0, 0.02)
    assert [r.kind for r in out] == [radiate.E2]
    other = [Row(transition(radiate.M1, 7.0, upper=up, lower=down))]
    assert [r.kind for r in calc_TP.merge(other, 0, 0.02)] == \
        [radiate.M1]


def test_two_different_level_pairs_are_not_merged_with_each_other():
    a, b, c = level(1, 0.0), level(2, 100.0), level(3, 200.0)
    given = [Row(transition(radiate.E2, 1.0, upper=b, lower=a)),
             Row(transition(radiate.M1, 1.0, upper=c, lower=a))]
    out = calc_TP.merge(given, 0, 0.02)
    assert len(out) == 2
    assert sorted(r.kind for r in out) == [radiate.E2, radiate.M1]


def test_the_trials_of_the_two_multipoles_add_trial_by_trial():
    up, down = level(2, 100.0), level(1, 0.0)
    given = [Row(transition(radiate.E2, 50.0, upper=up, lower=down,
                            values=(1.0, 2.0, 3.0))),
             Row(transition(radiate.M1, 50.0, upper=up, lower=down,
                            values=(10.0, 20.0, 30.0)))]
    out = calc_TP.merge(given, 3, 0.02)
    assert out[0].kind == 'M1+E2'
    assert out[0].a_values == {1: 11.0, 2: 22.0, 3: 33.0}


def test_a_wide_spread_keeps_a_small_multipole_in_the_merge():
    # The cut-off allows two standard deviations, so a fraction under
    # min_mixed_fraction survives when the trials disagree about it.
    up, down = level(2, 100.0), level(1, 0.0)
    given = [Row(transition(radiate.E2, 1.0, upper=up, lower=down,
                            values=(0.01, 4.0))),
             Row(transition(radiate.M1, 100.0, upper=up, lower=down,
                            values=(100.0, 100.0)))]
    out = calc_TP.merge(given, 2, 0.02)
    assert out[0].kind == 'M1+E2'


# --------------------------------------------------------------------------
# tying the deck, RCEOUT and OUTE together
# --------------------------------------------------------------------------

def test_the_deck_has_a_slot_for_every_fitted_parameter():
    folder = case('Pr3')
    deck = ing11.read_ing11(folder=folder)
    csets = rce.read_rceout(folder=folder)
    oute = rce.read_oute(folder=folder)
    for kset in (1, 2):
        cset, fitted, slots = calc_TP.check_parameters(
            deck, csets, oute, kset)
        assert len(slots) == len(fitted) == len(cset.parameters)


def test_oute_and_rceout_are_matched_name_by_name():
    folder = case('Pr3')
    deck = ing11.read_ing11(folder=folder)
    csets = rce.read_rceout(folder=folder)
    oute = rce.read_oute(folder=folder)
    oute[1] = oute[1][:-1]              # one parameter short
    with pytest.raises(calc_TP.RunError) as problem:
        calc_TP.check_parameters(deck, csets, oute, 1)
    assert 'different runs' in str(problem.value)


def test_a_missing_oute_is_a_named_error_not_a_crash():
    folder = case('Pr3')
    deck = ing11.read_ing11(folder=folder)
    csets = rce.read_rceout(folder=folder)
    with pytest.raises(calc_TP.RunError) as problem:
        calc_TP.check_parameters(deck, csets, {}, 1)
    assert 'OUTE' in str(problem.value)


# --------------------------------------------------------------------------
# varying the parameters
# --------------------------------------------------------------------------

def options(**changes):
    out = {s.key: s.default for s in calc_TP.SPECS}
    out.update(changes)
    return out


def parameter(name, flag, value, sigma):
    return rce.Parameter(name=name, flag=flag, value=value, sigma=sigma)


def test_a_parameter_with_no_standard_deviation_does_not_move():
    from toolbox import trials
    fitted = [parameter('EAV 4d7', -100, 21.7, 0.0),
              parameter('F2(11)', -100, 67.3, 0.0)]
    got = calc_TP.vary_parameters(fitted, options(), trials.Randomizer(1),
                                       trials.Histogram('x'))
    assert got == [21.7, 67.3]


def test_a_group_moves_by_one_common_ratio():
    from toolbox import trials
    fitted = [parameter('F2(11)', -3, 100.0, 5.0),
              parameter('F4(11)', -3, 50.0, 2.5),
              parameter('F6(11)', -4, 20.0, 1.0)]
    got = calc_TP.vary_parameters(fitted, options(),
                                       trials.Randomizer('g'),
                                       trials.Histogram('x'))
    assert got[0] / 100.0 == pytest.approx(got[1] / 50.0)
    assert got[2] / 20.0 != pytest.approx(got[0] / 100.0)


def test_without_grouping_every_parameter_moves_on_its_own():
    from toolbox import trials
    fitted = [parameter('F2(11)', -3, 100.0, 5.0),
              parameter('F4(11)', -3, 50.0, 2.5)]
    got = calc_TP.vary_parameters(fitted,
                                       options(var_grouped_params=False),
                                       trials.Randomizer('g'),
                                       trials.Histogram('x'))
    assert got[0] / 100.0 != pytest.approx(got[1] / 50.0)


def test_a_fixed_parameter_moves_only_when_it_is_asked_to():
    # The branch `calc_TP_M1E2.pl:582` could never reach: it tests
    # abs($flag == 100), which is 0 for a flag of -100.
    from toolbox import trials
    fitted = [parameter('ALPHA', -100, 10.0, 0.0)]
    assert calc_TP.vary_parameters(
        fitted, options(), trials.Randomizer(1), trials.Histogram('x')) \
        == [10.0]
    moved = calc_TP.vary_parameters(
        fitted, options(fixed_param_variance=0.1), trials.Randomizer(1),
        trials.Histogram('x'))
    assert moved != [10.0]


def test_the_r_integrals_move_by_the_fraction_asked_for():
    from toolbox import trials
    folder = case('U90')
    deck = ing11.read_ing11(folder=folder)
    cards = ing11.read_r_integrals(deck, order=1)
    rng = trials.Randomizer('r')
    got = calc_TP.vary_r_integrals(cards, options(r2_variance=0.15), rng,
                                        trials.Histogram('x'))
    moved = [(a, b.value) for a, b in zip(got, cards) if b.value]
    assert moved and all(a != b for a, b in moved)
    unmoved = calc_TP.vary_r_integrals(cards, options(r2_variance=0.0),
                                            rng, trials.Histogram('x'))
    assert unmoved == [one.value for one in cards]


# --------------------------------------------------------------------------
# one whole run
# --------------------------------------------------------------------------

def build(folder, code, ispecc=None, name='Forbid'):
    """RCN, RCN2, RCG, RCE over one WORK case's two input cards.

    `code` is the `cowan_code` fixture: the directory the permanent
    files are in, cfp decks and all.  The programs are run from CODE/
    but COWAN.CFG names that directory instead, because RCG creates
    FOR072, FOR073 and FOR074 wherever COWAN.CFG points and nothing may
    write into CODE/.

    `ispecc` overwrites column 72 of the IN2 main card, which chooses
    which of RCG's three line sorts are printed and so which layout the
    spectrum comes out in; the deck's own value is 2, the second-parity
    sort in the Lund layout.

    `name` is the case under WORK/ to build, and `''` is WORK itself --
    the Sn VIII installation test, which is the electric-dipole one.
    """
    source = WORK / name if name else WORK
    if not (source / 'IN2').is_file():
        pytest.skip('%s is not present' % source)
    folder.mkdir(parents=True, exist_ok=True)
    for name in ('IN2', 'IN36'):
        shutil.copyfile(source / name, folder / name)
    if ispecc is not None:
        deck = (folder / 'IN2').read_text().split('\n')
        card = deck[0].rstrip('\r')
        deck[0] = card[:71] + str(ispecc) + card[72:] + '\r'
        (folder / 'IN2').write_text('\n'.join(deck), newline='')
    (folder / 'COWAN.CFG').write_text(str(code) + os.sep + '\n')
    for name in EXECUTABLES:
        done = subprocess.run([str(CODE / name)], cwd=str(folder),
                              capture_output=True, text=True)
        assert done.returncode == 0, '%s: %s' % (name, done.stdout[-2000:])
    (folder / 'COWAN.CFG').unlink()
    return folder


def columns(path):
    return [line.split('\t') for line in
            path.read_text().rstrip('\n').split('\n')]


def test_a_whole_run_with_no_trials(tmp_path, cowan_code):
    folder = build(tmp_path, cowan_code)
    assert calc_TP.main([str(folder)]) == 0
    table = columns(folder / 'tp.txt')
    head, body = table[0], table[1:]
    assert head[:7] == ['type', 'conf1', 't1', 'J1', 'conf2', 't2', 'J2']
    assert body, 'the E2 spectrum produced no line'
    # The deck has DMIN = 0.00, so both multipoles print and the merge
    # is exercised on a real run: see WORK/Forbid/README.txt.
    assert {row[0] for row in body} == {radiate.M1, radiate.E2, 'M1+E2'}
    # Every level of this deck is 4d7, and none of them is observed.
    assert all(row[1] == '4d7' and row[4] == '4d7' for row in body)
    # The columns are found by name: lid1 and lid2 were added between
    # the calculated energies and the wavelength.
    at = head.index
    assert all(row[at('E_low_exp')] == '' and row[at('E_up_exp')] == ''
               for row in body)
    # Sorted by the vacuum wavelength the fitted energies give.  The
    # printed column is that wavelength in standard air between 2000 and
    # 20000 A, and in vacuum outside, as every line list quotes it.
    vacuum = [1e5 / abs(float(r[at('E_up')]) - float(r[at('E_low')]))
              for r in body]
    assert vacuum == sorted(vacuum)
    assert any(2000.0 < x < 20000.0 for x in vacuum), 'no air line in the case'
    for row, x in zip(body, vacuum):
        wanted = vacair.to_air(x) if 2000.0 < x < 20000.0 else x
        assert float(row[at('wl_c(A)')]) == pytest.approx(wanted, abs=5e-5)


def test_the_level_list_carries_the_lids(tmp_path, cowan_code):
    """`lid1` and `lid2` name a level of the file written beside them.

    The pair is what `insert_level_ids` reads, so it is read here the
    way that program reads it -- by column name, through `levids.Table`
    -- rather than by position.
    """
    folder = build(tmp_path, cowan_code)
    assert calc_TP.main([str(folder)]) == 0

    table = levids.Table.read(folder / 'tp_levels.txt', levids.LEVEL_KEYS)
    targets, by_number = levids.read_levels(table)
    assert targets, 'the level list is empty'
    # Every level of the run is numbered, whether or not a line reached
    # it, and the numbers run 1..N with none missing.
    assert sorted(by_number) == list(range(1, len(by_number) + 1))
    # This deck is 4d7 throughout: one parity set, even, and no
    # observed level.
    assert {one.parity for one in targets} == {'e'}
    assert all(one.calculated for one in targets)

    lines = levids.Table.read(folder / 'tp.txt', levids.LINE_KEYS)
    head = columns(folder / 'tp.txt')[0]
    at = head.index
    assert head[9:14] == ['E_low', 'E_up', 'lid1', 'lid2', 'wl_c(A)']
    seen = 0
    for row in lines.rows:
        for lid, jcol in (('lid1', 'J1'), ('lid2', 'J2')):
            number = int(lines.field(row, lines.index(lid)))
            level = by_number[number]
            assert level.j == levids.j_key(row[at(jcol)].strip())
            seen += 1
    assert seen, 'the run produced no line'


def test_insert_level_ids_reads_what_calc_TP_writes(tmp_path, cowan_code):
    """The whole point of the two new columns: the join can be made."""
    folder = build(tmp_path, cowan_code)
    assert calc_TP.main([str(folder)]) == 0

    # A source list of the kind `subst_exp_levs` takes, naming two of
    # the run's own levels by (parity, J, energy).
    table = columns(folder / 'tp_levels.txt')
    head, chosen = table[0], table[1:3]
    at = head.index
    source = folder / 'exp_levs.txt'
    source.write_text('\n'.join(
        ['\t'.join(['level_id', 'E', 'J', 'parity'])]
        + ['\t'.join(['%d' % (700 + n), row[at('Ee')].strip(),
                      row[at('J')].strip(), row[at('par')].strip()])
           for n, row in enumerate(chosen)]) + '\n')

    assert insert_level_ids.main(
        [str(source), str(folder / 'tp_levels.txt'), str(folder / 'tp.txt'),
         '--ecol', 'E', '--jcol', 'J', '--pcol', 'parity',
         '--units', 'kK', '--tol', '1']) == 0

    head = columns(folder / 'tp_levels.txt')[0]
    assert 'level_id' in head
    head = columns(folder / 'tp.txt')[0]
    assert head[head.index('lid2') + 1:head.index('lid2') + 3] == ['id1',
                                                                   'id2']
    wanted = {'%d' % (700 + n) for n in range(len(chosen))}
    got = set()
    for row in columns(folder / 'tp.txt')[1:]:
        got.update(x for x in row[head.index('id1'):head.index('id2') + 1]
                   if x.strip())
    assert got and got <= wanted


def test_a_wavelength_sorted_run_gives_the_same_table(tmp_path, cowan_code):
    """ISPECC = 4 prints the wide layout and nothing else.

    Column 72 of the main card chooses the sorts: 4 asks for the
    wavelength sort alone, and that one always prints the wide layout
    whatever ILUND says (`rcg11kd.f:14019`).  The same lines come out of
    it, to the extra digit the wider fields carry -- the wavelength sort
    prints F11.4 energies and an E11.4 gA where the Lund layout prints
    F11.3 and E10.3.
    """
    lund = build(tmp_path / 'lund', cowan_code, ispecc=2)
    wide = build(tmp_path / 'wide', cowan_code, ispecc=4)
    for folder in (lund, wide):
        assert calc_TP.main([str(folder)]) == 0

    def by_pair(folder):
        table = columns(folder / 'tp.txt')
        return {tuple(row[:7]): row for row in table[1:]}, table[0]

    one, head = by_pair(lund)
    two, head_again = by_pair(wide)
    assert head == head_again
    assert set(one) == set(two)
    assert one, 'the run produced no line'
    a = head.index('A(s-1)')
    for key, row in one.items():
        assert float(row[a]) == pytest.approx(float(two[key][a]), rel=2e-3)


def test_a_whole_run_with_trials_is_repeatable(tmp_path, cowan_code,
                                               cowan_rcg):
    folder = build(tmp_path, cowan_code)
    # The trials run RCG themselves, and it is the copy of it that has
    # to be run: run_rcg writes COWAN.CFG naming the directory the
    # executable sits in, so the one in CODE/ would put the cfp decks
    # back there.
    parfile.write(folder / calc_TP.PAR_FILE, calc_TP.SPECS,
                  {'num_trials': 2, 'store_trials': 2, 'random_seed': 'test',
                   'rcg': str(cowan_rcg)})
    assert calc_TP.main([str(folder)]) == 0
    first = (folder / 'tp.txt').read_text()
    head = columns(folder / 'tp.txt')[0]
    assert 'std_A%' in head and 'Trial A data' in head
    assert (folder / 'trial1' / 'OUTG11').is_file()
    assert not (folder / 'trial1' / 'COWAN.CFG').exists()

    # The same seed over the trials that are already there: the numbers
    # come back the same, and RCG is not run again.
    (folder / calc_TP.PAR_FILE).unlink()
    parfile.write(folder / calc_TP.PAR_FILE, calc_TP.SPECS,
                  {'num_trials': 2, 'store_trials': 2, 'random_seed': 'test',
                   'read_only': 2, 'rcg': str(cowan_rcg)})
    assert calc_TP.main([str(folder)]) == 0
    assert (folder / 'tp.txt').read_text() == first


def test_the_trial_decks_differ_from_the_original_only_in_the_r2_cards(
        tmp_path, cowan_code, cowan_rcg):
    folder = build(tmp_path, cowan_code)
    parfile.write(folder / calc_TP.PAR_FILE, calc_TP.SPECS,
                  {'num_trials': 1, 'store_trials': 1, 'random_seed': 'deck',
                   'rcg': str(cowan_rcg)})
    assert calc_TP.main([str(folder)]) == 0
    original = (folder / 'ING11').read_text().split('\n')
    trial = (folder / 'trial1' / 'ING11').read_text().split('\n')
    assert len(original) == len(trial)
    differ = [i for i, (a, b) in enumerate(zip(original, trial)) if a != b]
    assert differ, 'the trial deck is identical to the original'
    for i in differ:
        assert ing11.is_r_integral(original[i].rstrip('\r'))
        # Only the value changed; the tag and the scale factors did not.
        assert original[i][:38] == trial[i][:38]
        assert original[i][50:] == trial[i][50:]


def test_a_missing_parameter_file_is_written_rather_than_demanded(
        tmp_path, monkeypatch):
    # It used to refuse and tell the user to run --write-par.  There is
    # nothing useful to run without a par file, so the program writes
    # one and says so.  The name is moved aside for the test because the
    # copy CODE ships would otherwise be found first, which is the whole
    # point of shipping it.
    monkeypatch.setattr(calc_TP, 'PAR_FILE', 'absent_here.par')
    assert calc_TP.main([str(tmp_path)]) == 0
    got = parfile.read(tmp_path / 'absent_here.par', calc_TP.SPECS)
    assert got == {spec.key: spec.default for spec in calc_TP.SPECS}


def test_the_program_writes_a_parameter_file_on_request(tmp_path):
    assert calc_TP.main([str(tmp_path), '--write-par']) == 0
    got = parfile.read(tmp_path / calc_TP.PAR_FILE, calc_TP.SPECS)
    assert got == {s.key: s.default for s in calc_TP.SPECS}
# --------------------------------------------------------------------------
# censoring: the trial window and the reporting cut
# --------------------------------------------------------------------------

def censored(dmin, window=True, **changes):
    """A one-card deck with a DMIN and, by default, a KCPL=3 window."""
    prelude = ['    3 0  0  80016000'] if window else []
    card = ('    1   -1412124 8 80         001111181 00000000 0 1000.0000'
            '                 2       2 12 12  8  8')
    main = card[:60] + ('%5.1f' % dmin) + card[65:]
    deck = ing11.Ing11(prelude=prelude, main=ing11.read_main(main))
    settings = dict(trial_dmin_offset=1.0, trial_lambda_margin=0.05,
                    trial_lambda_margin_max=10.0)
    settings.update(changes)
    return deck, settings


def test_the_trials_are_run_a_decade_below_the_punched_dmin():
    deck, settings = censored(8.0)
    cut = calc_TP.censor(deck, settings, (radiate.E1,), 5)
    assert deck.main.dmin == 7.0
    assert cut == ('gA', 1e8)


def test_the_trial_window_is_wider_at_both_ends():
    deck, settings = censored(8.0)
    calc_TP.censor(deck, settings, (radiate.E1,), 5)
    low, high = deck.window
    assert low < 800.0 and high > 16000.0
    assert 800.0 / low == pytest.approx(1.05, abs=2e-3)
    assert high / 16000.0 == pytest.approx(1.05, abs=2e-3)


def test_a_forbidden_only_run_lifts_dmin_altogether():
    # Author decision, 2026-08-27: with E1 in the run the forbidden
    # multipoles share its DMIN; on their own they have no cut.
    deck, settings = censored(8.0)
    cut = calc_TP.censor(deck, settings, (radiate.M1, radiate.E2), 5)
    assert deck.main.dmin == 0.0
    assert cut is None


def test_a_zero_offset_leaves_the_deck_as_it_was():
    deck, settings = censored(8.0, trial_dmin_offset=0.0,
                              trial_lambda_margin=0.0)
    was = deck.main.raw, list(deck.prelude)
    assert calc_TP.censor(deck, settings, (radiate.E1,), 5) is None
    assert (deck.main.raw, deck.prelude) == was


def test_a_run_with_no_trials_touches_nothing():
    deck, settings = censored(8.0)
    was = deck.main.raw, list(deck.prelude)
    assert calc_TP.censor(deck, settings, (radiate.E1,), 0) is None
    assert (deck.main.raw, deck.prelude) == was


def test_a_dmin_that_cuts_on_the_matrix_element_is_not_lowered():
    # Below 1.0 the field cuts on D0, which this program cannot
    # re-apply at write time, so it is left alone.
    deck, settings = censored(0.5)
    assert calc_TP.censor(deck, settings, (radiate.E1,), 5) is None
    assert deck.main.dmin == 0.5


def test_the_trial_dmin_never_leaves_RCG_exponent_range():
    deck, settings = censored(8.0, trial_dmin_offset=20.0)
    cut = calc_TP.censor(deck, settings, (radiate.E1,), 5)
    assert deck.main.dmin == calc_TP.DMIN_FLOOR
    assert ing11.dmin_cut(deck.main.dmin)[0] == 'gA'
    assert cut == ('gA', 1e8)


def test_a_deck_with_no_window_card_keeps_its_dmin_treatment():
    deck, settings = censored(8.0, window=False)
    assert calc_TP.censor(deck, settings, (radiate.E1,), 5) == ('gA', 1e8)
    assert deck.main.dmin == 7.0


# --------------------------------------------------------------------------
# the total decay rate: RCG's own SUMGA, and the fall-back
# --------------------------------------------------------------------------

class Silent:
    """A log that keeps its warnings where a test can look at them."""

    def __init__(self):
        self.warnings = []

    def note(self, text):
        pass

    def warn(self, text):
        self.warnings.append(text)


def read_run(folder, multipoles=(radiate.M1, radiate.E2)):
    """The unvaried run of `folder`, as `calc_TP.run` reads it."""
    from toolbox import case as case_mod
    deck = ing11.read_ing11(folder=folder)
    of_label = calc_TP._kset_of(deck)
    base = case_mod.Case.read(folder, couplings=('LS',), use_emax=True)
    mapping = {(o.kset, o.j, o.rcg): o for o in base.identified}
    found, table, order = {}, {}, []
    calc_TP.read_transitions(folder, mapping, of_label, (1,), 0, 0,
                             found, table, order, multipoles)
    return mapping, of_label, found, multipoles


def test_a_level_prefers_rcgs_own_total_to_the_sum_of_its_lines():
    one = level(1, 0.0, j=3.5)
    one.decay(0, 1.0, 1.0)
    one.decay(0, 2.0, 2.0)
    assert one.total(0) == pytest.approx(3.0)
    one.summed(0, 80.0)                  # gA, over 2J+1 = 8
    assert one.total(0) == pytest.approx(10.0)


def test_the_totals_of_the_multipoles_add():
    one = level(1, 0.0, j=3.5)
    one.summed(0, 80.0)
    one.summed(0, 8.0)
    assert one.total(0) == pytest.approx(11.0)


def test_a_level_whose_sumga_overflowed_falls_back_to_its_lines():
    # And it stays fallen back: a total missing one multipole would be
    # worse than one that is merely censored by DMIN.
    one = level(1, 0.0, j=3.5)
    one.decay(0, 4.0, 4.0)
    one.summed(0, 80.0)
    one.summed(0, None)
    assert one.total(0) == pytest.approx(4.0)
    one.summed(0, 8.0)
    assert one.total(0) == pytest.approx(4.0)


def test_a_level_with_no_total_at_all_uses_its_lines():
    one = level(1, 0.0, j=3.5)
    one.decay(0, 4.0, 4.0)
    assert one.total(0) == pytest.approx(4.0)
    assert one.total(3) == 0.0


def test_the_totals_of_a_real_run_are_read_and_used(tmp_path, cowan_code):
    folder = build(tmp_path, cowan_code)
    mapping, of_label, found, multipoles = read_run(folder)
    log = Silent()
    assert calc_TP.read_totals(folder, mapping, of_label, (1,), 0, found,
                               multipoles, log)
    assert not log.warnings
    # WORK/Forbid punches DMIN = 0.0 and no wavelength window, so
    # nothing is censored and SUMGA is the sum over the printed lines,
    # to the four figures its E10.3 field carries.
    checked = 0
    for one in found.values():
        internal = one.ar_tot[0] if one.ar_tot else 0.0
        if not internal:
            continue
        assert one.total(0) == pytest.approx(internal, rel=2e-3)
        checked += 1
    assert checked, 'no level had a decay to compare'


def test_a_starred_sumga_is_reported_and_not_used(tmp_path, cowan_code):
    folder = build(tmp_path, cowan_code)
    text = (folder / 'OUTG11').read_text(errors='replace')
    at = text.index('SUMGA=')
    text = text[:at] + 'SUMGA=**********' + text[at + 16:]
    (folder / 'OUTG11').write_text(text, newline='')
    mapping, of_label, found, multipoles = read_run(folder)
    log = Silent()
    calc_TP.read_totals(folder, mapping, of_label, (1,), 0, found,
                        multipoles, log)
    assert len(log.warnings) == 1
    assert 'overflowed' in log.warnings[0]
    starred = [one for one in found.values()
               if one.ga_tot and one.ga_tot[0] is None]
    assert len(starred) == 1
    assert starred[0].total(0) == pytest.approx(starred[0].ar_tot[0])


# --------------------------------------------------------------------------
# what the list holds
# --------------------------------------------------------------------------

def test_a_row_knows_whether_the_unvaried_run_had_it():
    up, down = level(2, 100.0), level(1, 0.0)
    one = transition(radiate.E2, 7.0, upper=up, lower=down)
    assert Row(one).in_base
    one.ga = None                        # a transition only a trial had
    assert not Row(one).in_base


def test_a_merged_row_is_in_the_base_if_either_multipole_was():
    up, down = level(2, 100.0), level(1, 0.0)
    e2 = transition(radiate.E2, 50.0, upper=up, lower=down, values=(50.0,))
    m1 = transition(radiate.M1, 50.0, upper=up, lower=down, values=(50.0,))
    m1.ga = None
    out = calc_TP.merge([Row(e2), Row(m1)], 1, 0.02)
    assert len(out) == 1 and out[0].in_base
    e2.ga = None
    out = calc_TP.merge([Row(e2), Row(m1)], 1, 0.02)
    assert len(out) == 1 and not out[0].in_base


def test_the_list_holds_the_unvaried_runs_transitions(tmp_path, cowan_code,
                                                      cowan_rcg):
    """The unvaried run's transitions, all of them and only them.

    The trials are run in a wider window than the deck punched, so they
    hold lines the unvaried run does not; those are not results.  Here
    the same deck is read with trials and without them, and the two
    lists must hold the same level pairs.
    """
    folder = build(tmp_path, cowan_code)
    parfile.write(folder / calc_TP.PAR_FILE, calc_TP.SPECS,
                  {'num_trials': 1, 'store_trials': 1, 'random_seed': 'base',
                   'rcg': str(cowan_rcg), 'out_file': 'with_trials.txt'})
    assert calc_TP.main([str(folder)]) == 0
    # On the level pair, not on the name: a trial can move the E2
    # fraction of a line across min_mixed_fraction and rename it M1+E2,
    # which is the merge doing its job and not a change of membership.
    with_trials = {tuple(row[1:7])
                   for row in columns(folder / 'with_trials.txt')[1:]}

    (folder / calc_TP.PAR_FILE).unlink()
    parfile.write(folder / calc_TP.PAR_FILE, calc_TP.SPECS,
                  {'num_trials': 0, 'out_file': 'alone.txt'})
    assert calc_TP.main([str(folder)]) == 0
    alone = {tuple(row[1:7]) for row in columns(folder / 'alone.txt')[1:]}
    assert with_trials and with_trials == alone


# --------------------------------------------------------------------------
# the warning log
# --------------------------------------------------------------------------

def test_the_log_carries_the_seed_and_counts_only_warnings(tmp_path, capsys):
    log = calc_TP.Log(tmp_path / 'calc_TP.log')
    log.start('a seed', tmp_path)
    log.note('nothing to see')
    log.warn('something to see')
    log.finish()
    text = (tmp_path / 'calc_TP.log').read_text()
    assert "'a seed'" in text
    assert 'nothing to see' in text
    assert 'WARNING: something to see' in text
    assert '1 warning(s).' in text
    assert '1 warning(s)' in capsys.readouterr().out


def test_a_run_writes_its_log(tmp_path, cowan_code):
    folder = build(tmp_path, cowan_code)
    assert calc_TP.main([str(folder)]) == 0
    assert (folder / calc_TP.LOG_FILE).is_file()


# --------------------------------------------------------------------------
# --cleanup
# --------------------------------------------------------------------------

def trial_dirs(folder, how_many):
    for n in range(1, how_many + 1):
        one = folder / ('trial%d' % n)
        one.mkdir()
        (one / 'OUTG11').write_text('x', newline='')


def test_cleanup_removes_the_trial_directories_and_nothing_else(tmp_path):
    trial_dirs(tmp_path, 5)
    (tmp_path / 'tp.txt').write_text('kept', newline='')
    (tmp_path / 'ING11').write_text('kept', newline='')
    assert calc_TP.main([str(tmp_path), '--cleanup', '--yes']) == 0
    assert not list(tmp_path.glob('trial*'))
    assert (tmp_path / 'tp.txt').is_file()
    assert (tmp_path / 'ING11').is_file()


def test_cleanup_asks_before_it_deletes(tmp_path, monkeypatch):
    trial_dirs(tmp_path, 2)
    monkeypatch.setattr('builtins.input', lambda prompt: 'n')
    assert calc_TP.main([str(tmp_path), '--cleanup']) == 0
    assert len(list(tmp_path.glob('trial*'))) == 2
    monkeypatch.setattr('builtins.input', lambda prompt: 'y')
    assert calc_TP.main([str(tmp_path), '--cleanup']) == 0
    assert not list(tmp_path.glob('trial*'))


def test_cleanup_with_nothing_to_clean_says_so(tmp_path):
    with pytest.raises(calc_TP.RunError) as problem:
        calc_TP.main([str(tmp_path), '--cleanup', '--yes'])
    assert 'no trial directory' in str(problem.value)


def test_cleanup_runs_nothing_else(tmp_path):
    # No par file is read and no list is written: it is a flag that
    # deletes and does nothing more.
    trial_dirs(tmp_path, 1)
    assert calc_TP.main([str(tmp_path), '--cleanup', '--yes']) == 0
    assert not (tmp_path / 'tp.txt').exists()
    assert not (tmp_path / calc_TP.PAR_FILE).exists()


def test_cleanup_takes_no_other_option(tmp_path):
    trial_dirs(tmp_path, 1)
    with pytest.raises(calc_TP.Usage):
        calc_TP.main([str(tmp_path), '--cleanup', '--write-par'])
    with pytest.raises(calc_TP.Usage):
        calc_TP.main([str(tmp_path), '--yes'])


# --------------------------------------------------------------------------
# reading a spectrum that a trial has rearranged
# --------------------------------------------------------------------------
#
# Two things a real OUTG11 does that the fixtures above cannot: it puts
# the two levels of a line the other way round when the fit has moved an
# unobserved level past its neighbour, and it prints every line twice
# when ISPECC asks for both level sorts.  Both are fed in here as a
# synthetic spectrum, `spectrum.lines` being the one seam.

def synthetic(monkeypatch, printed):
    """Drive `read_transitions` off a made-up list of printed lines."""
    from toolbox import spectrum

    def ref(number, j, energy):
        return spectrum.LevelRef(energy=energy, j=j, config=1, term='x' * 8,
                                 level=number, name='Sn 8  4d10 5s   ')

    def lines(folder=None, multipoles=()):
        for a, b, ga in printed:
            up, low = ref(*a), ref(*b)
            delta = abs(up.energy - low.energy)
            yield spectrum.Line(multipole='ELEC QUD', initial=up, final=low,
                                delta_e=delta, wavelength=1e5 / delta,
                                log_gf=0.0, ga=ga, cf=0.5)

    monkeypatch.setattr(calc_TP.spectrum, 'lines', lines)


def pairing(monkeypatch, printed, energies, table=None, trial=0):
    """`read_transitions` over `printed`, with `energies` per RCE level."""
    synthetic(monkeypatch, printed)
    mapping, found = {}, {}
    for number, energy in energies.items():
        one = level(number, energy, j=2.0)
        mapping[(1, 2.0, number)] = one.one
        found[(1, 2.0, number)] = one
    order = []
    table = {} if table is None else table
    got = calc_TP.read_transitions('.', mapping, {'4d10 5s': 1}, (1,), 0,
                                   trial, found, table, order,
                                   (radiate.E2,), log=Silent())
    return got, table, order, found


def test_a_run_that_prints_a_transition_twice_counts_it_once(monkeypatch):
    # ISPECC=3 asks for both level sorts, and both print every line.
    # That is worth having -- it is what gives every level a SUMGA -- so
    # the repeat is folded, not refused.
    printed = [((1, 2.0, 10.0), (2, 2.0, 0.0), 1e8),
               ((2, 2.0, 0.0), (1, 2.0, 10.0), 1e8)]
    (swapped, again), table, order, found = pairing(
        monkeypatch, printed, {1: 10.0, 2: 0.0})
    assert (swapped, again) == (0, 1)
    assert len(order) == 1
    assert found[(1, 2.0, 1)].ar_tot[0] == pytest.approx(1e8 / 5)


def test_a_second_print_that_disagrees_is_warned_about(monkeypatch):
    printed = [((1, 2.0, 10.0), (2, 2.0, 0.0), 1e8),
               ((1, 2.0, 10.0), (2, 2.0, 0.0), 2e8)]
    from toolbox import spectrum
    synthetic(monkeypatch, printed)
    log = Silent()
    mapping, found = {}, {}
    for number, energy in ((1, 10.0), (2, 0.0)):
        one = level(number, energy, j=2.0)
        mapping[(1, 2.0, number)] = one.one
        found[(1, 2.0, number)] = one
    calc_TP.read_transitions('.', mapping, {'4d10 5s': 1}, (1,), 0, 0,
                             found, {}, [], (radiate.E2,), log=log)
    assert len(log.warnings) == 1
    assert 'printed twice' in log.warnings[0]


def test_a_trial_that_reverses_the_two_levels_still_matches_the_base(
        monkeypatch):
    # Level 1 is above level 2 in the base run and below it in the
    # trial: the same physical line, and it must not become a new one.
    printed = [((1, 2.0, 10.0), (2, 2.0, 9.0), 1e8)]
    _, table, order, _ = pairing(monkeypatch, printed, {1: 10.0, 2: 9.0})
    assert len(order) == 1

    reversed_print = [((2, 2.0, 10.0), (1, 2.0, 9.0), 1e8)]
    (swapped, again), table, order, _ = pairing(
        monkeypatch, reversed_print, {1: 9.0, 2: 10.0}, table=table, trial=1)
    assert (swapped, again) == (1, 0)
    assert order == []                        # no trial-only transition
    assert table[((1, 2.0, 1), (1, 2.0, 2), radiate.E2)].n_seen == 1


# --------------------------------------------------------------------------
# the blank standard uncertainty of a line the trials barely saw
# --------------------------------------------------------------------------

def defaults(**changed):
    """The shipped parameter defaults, with `changed` applied."""
    out = {spec.key: spec.default for spec in calc_TP.SPECS}
    out.update(changed)
    return out


def statistics(values, store=5):
    """The `_statistics` columns of a row whose trials are `values`."""
    up, down = level(2, 100.0), level(1, 0.0)
    one = transition(radiate.E2, 7.0, upper=up, lower=down, values=values)
    row = Row(one)
    up.decay(0, 7.0 / (2 * up.j + 1), 7.0 / (2 * up.j + 1))
    for at, a in enumerate(values, 1):
        up.decay(at, a, a)
    options = defaults(store_trials=store, print_dbf=1)
    got = calc_TP._statistics(row, up, options, 0.5)
    return dict(zip(['n_trials', 'd_BF', 'A_mean', 'S_mean', 'std_A%',
                     'cf_mean', 'std_cf'], got))


def test_a_line_no_trial_printed_quotes_no_uncertainty():
    got = statistics(())
    assert got['n_trials'] == '0'
    assert got['std_A%'] == '' and got['std_cf'] == '' and got['d_BF'] == ''


def test_a_line_one_trial_printed_quotes_no_uncertainty():
    # One sample has no spread, and a zero here would read as a line
    # that is known exactly.
    got = statistics((1.4,))
    assert got['n_trials'] == '1'
    assert got['std_A%'] == '' and got['std_cf'] == '' and got['d_BF'] == ''
    assert got['S_mean'] and got['A_mean']       # the mean itself stands


def test_two_trials_are_enough_for_a_spread():
    got = statistics((1.0, 2.0))
    assert got['n_trials'] == '2'
    assert float(got['std_A%']) > 0.0
    assert got['std_cf'] != '' and got['d_BF'] != ''


# --------------------------------------------------------------------------
# the long-wavelength trial margin
# --------------------------------------------------------------------------

def config(kset, number, eav):
    from toolbox.outg11 import ConfigCard
    return ConfigCard(kset=kset, number=number, shells=(), nelec=1,
                      parity=1 if kset == 1 else -1, terms='ALL',
                      name='X  1  4f2         ', eav=eav, kinetic=0.0)


#: Three configurations on a scale like Cowan's own: a large negative
#: total electronic energy, of which only the differences are physics.
GROUND = -2.0e6
CARDS = {(1, 1): config(1, 1, GROUND),
         (1, 2): config(1, 2, GROUND + 100.0),
         (2, 1): config(2, 1, GROUND + 110.0)}


def test_the_margin_is_the_movement_of_the_smallest_transition_energy():
    # The pair (1,2)-(2,1) is 10 kK apart and sits 100 and 110 kK above
    # the ground configuration, so 5 per cent of each end is 0.05*210/10.
    fraction, best = calc_TP.long_margin([((1, 2), (2, 1))], CARDS,
                                         0.05, 10.0)
    assert fraction == pytest.approx(0.05 * 210.0 / 10.0)
    assert (best[1].kset, best[1].number) == (2, 1)
    assert (best[2].kset, best[2].number) == (1, 2)


def test_the_energies_are_measured_from_the_lowest_configuration():
    # Taken on Cowan's own scale the same pair would give
    # 0.05*4e6/10 = 20000, and every deck would sit on the ceiling.
    fraction, _ = calc_TP.long_margin([((1, 2), (2, 1))], CARDS, 0.05, 1e9)
    assert fraction < 2.0


def test_the_closest_pair_of_the_spectrum_is_the_one_that_counts():
    pairs = [((1, 1), (2, 1)), ((1, 2), (2, 1))]
    fraction, best = calc_TP.long_margin(pairs, CARDS, 0.05, 1e9)
    assert (best[1].kset, best[1].number) == (2, 1)
    assert (best[2].kset, best[2].number) == (1, 2)
    assert fraction == pytest.approx(0.05 * 210.0 / 10.0)


def test_two_configurations_of_the_same_energy_are_held_by_the_ceiling():
    cards = dict(CARDS)
    cards[(2, 1)] = config(2, 1, GROUND + 100.0)
    log = Silent()
    fraction, _ = calc_TP.long_margin([((1, 2), (2, 1))], cards, 0.05, 10.0,
                                      log=log)
    assert fraction == 10.0


def test_a_spectrum_with_nothing_to_measure_on_keeps_the_plain_margin():
    assert calc_TP.long_margin([], CARDS, 0.05, 10.0) == (0.05, None)
    blank = {(1, 1): config(1, 1, 0.0)}
    assert calc_TP.long_margin([((1, 1), (1, 1))], blank, 0.05, 10.0) \
        == (0.05, None)


# --------------------------------------------------------------------------
# the total as the rescaled sum times the SUMGA ratio
# --------------------------------------------------------------------------

def test_the_censoring_correction_is_taken_as_a_ratio():
    # Two printed lines, rescaled to the observed wavelength: 3.0
    # against 2.0 at RCG's own energies.  SUMGA counts a third line
    # DMIN kept out of the list, so it is 3.0 where the printed lines
    # sum to 2.0 -- and the answer is the rescaled sum times 3/2.
    one = level(1, 0.0, j=3.5)
    one.decay(0, 2.0, 1.0)
    one.decay(0, 1.0, 1.0)
    one.summed(0, 3.0 * 8)                      # gA, over 2J+1 = 8
    assert one.total(0) == pytest.approx(3.0 * 3.0 / 2.0)


def test_nothing_rescaled_leaves_sumga_as_it_stands():
    one = level(1, 0.0, j=3.5)
    one.decay(0, 2.0, 2.0)
    one.summed(0, 3.0 * 8)
    assert one.total(0) == pytest.approx(3.0)


def test_a_level_whose_printed_lines_sum_to_nothing_takes_sumga_alone():
    one = level(1, 0.0, j=3.5)
    one.decay(0, 5.0, 0.0)                      # every end unobserved
    one.summed(0, 3.0 * 8)
    assert one.total(0) == pytest.approx(3.0)


# --------------------------------------------------------------------------
# RCG's own line strength, against the one radiate derives from gA
# --------------------------------------------------------------------------
#
# The Lund layout calc_TP normally reads prints no line strength at all,
# for any multipole, so the check is made in the wide layout the
# wavelength sort always prints (ISPECC >= 4: FORMAT 56 header, FORMAT 58
# record).  Its `S/PMAX**2` column is `S2(I)`, and the true line strength
# is
#
#     S = S2 * PMAX**2
#
# -- `PSQ = GFCF*(PMAX**2)` and `GF = PSQ*S2(I)*ERAS1` at
# `rcg11kd.f:13448,13896`, with PMAX the largest reduced matrix element of
# the section, computed once at label 150 and printed once at the end of
# it by FORMAT 96 as `FOR PMAX=` and an F8.5.  The E2 arrangement
# (`GFCF*1.D16` against the later `/FLAM**2`) cancels; the autoionization
# case IGEN >= 5 carries an extra factor 3 in PSQ and is not exercised
# here.  Every other quantity the check needs -- LAMBDA(A) and GA(SEC-1)
# -- is on the record itself, so it is self-contained in one OUTG11 and
# needs nothing from calc_TP's own flow.
#
# PMAX IS PRINTED FOR THE ELECTRIC DIPOLE SPECTRUM ONLY.  `IPTG` is what
# gates FORMAT 96, and `rcg11kd.f:13341` reads
#
#       IF (II.NE.1) IPTG=0
#
# with II the multipole flag of SPECTR -- 0 M1, 1 E1, 2 E2 -- so an M1 or
# an E2 section prints neither PMAX nor the PMUP table it is the maximum
# of.  The absolute check is therefore made on E1 and the two forbidden
# multipoles get the strongest form that survives without it: within one
# section PMAX is a single number, so `S(gA)/S2` must be CONSTANT over
# the section, and that constant is PMAX**2.  Constancy is not a weak
# test here -- the section spans a wide range of wavelengths and of
# 2J+1, so a wrong power of the wavelength or a missing statistical
# weight breaks it at once, and the wavelength power is exactly what
# differs between the multipoles.  What it cannot see is a wrong
# constant prefactor, and only printing PMAX for II /= 1 would fix that.
#
# Precision.  S2 is F8.3 -- three decimals, ABSOLUTE -- so a line is
# worth 5e-4/S2 in relative terms: 5e-4 at S2 = 1, half a per cent at
# S2 = 0.1, and nothing at all below S2 ~ 1e-3, which prints as 0.000.
# PMAX at F8.5 adds ~1e-5, doubled by the square, which is negligible
# beside that.  So the check is a part in 2000 on the strong lines and is
# asserted at one per cent, on lines with S2 >= 0.1 only.  It is not a
# bit-exactness check and is not to be sold as one: what it exists to
# catch -- a wrong multipole coefficient, a wrong power of the
# wavelength, a missing 2J+1, an angstrom/cm-1 slip -- is a factor of a
# few at least, and usually orders of magnitude.

_PMAX = 'FOR PMAX='
_S2_FLOOR = 0.1                  # below this the F8.3 field says too little


def strengths(path):
    """Every wide-layout section: (multipole, PMAX or None, [(S2, S)]).

    `S` is the line strength `radiate` makes of the printed gA and
    wavelength; PMAX arrives after the lines of its section, and is
    None for the sections RCG prints none for.
    """
    from toolbox import spectrum
    multipole, layout, held, pmax = None, None, [], None

    def section():
        return (multipole, pmax, held)

    with open(path, 'r', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n').rstrip('\r')
            found = spectrum._heading(line)
            if found is not None:
                if multipole is not None and held:
                    yield section()
                multipole, layout, held, pmax = found[0], None, [], None
                continue
            if multipole is None:
                continue
            if line == spectrum._WIDE_COLUMNS:
                layout = 'wide'
                continue
            at = line.find(_PMAX)
            if at >= 0:
                pmax = float(line[at + len(_PMAX):at + len(_PMAX) + 8])
                continue
            if layout is None or not spectrum._is_wide(line):
                continue
            # FORMAT 58: S2 is the F8.3 between LAMBDA(A) and the F9.4
            # GF, which `spectrum._wide` steps over.
            wavelength, s2 = line[76:88], line[88:96]
            if '*' in wavelength or '*' in s2:
                continue
            ga = float(line[113:124])
            if float(s2) < _S2_FLOOR or ga <= 0.0:
                continue
            lam = float(wavelength)
            held.append((float(s2),
                         radiate.s_from_ga(ga, lam,
                                           radiate.OF_SPECTRUM[multipole])))
    if multipole is not None and held:
        yield section()


@pytest.mark.parametrize('name, wanted',
                         [('', {radiate.E1}),
                          ('Forbid', {radiate.M1, radiate.E2})])
def test_rcgs_own_line_strength_is_the_one_radiate_derives(
        tmp_path, cowan_code, name, wanted):
    folder = build(tmp_path, cowan_code, ispecc=4, name=name)
    seen, checked = set(), 0
    for spectrum_name, pmax, lines in strengths(folder / 'OUTG11'):
        multipole = radiate.OF_SPECTRUM[spectrum_name]
        seen.add(multipole)
        ratios = [s / s2 for s2, s in lines]
        if pmax is None:
            # No printed PMAX: the section's own implied one, which must
            # be one number for all of it.
            assert min(ratios) > 0.0
            assert max(ratios) == pytest.approx(min(ratios), rel=0.01), \
                'S(gA)/S2 is not constant over the %s section' % spectrum_name
        else:
            for ratio in ratios:
                assert ratio == pytest.approx(pmax * pmax, rel=0.01)
        checked += len(ratios)
    assert wanted <= seen, 'the deck printed %s' % (seen or 'no spectrum')
    assert checked, 'no line was strong enough to check'
