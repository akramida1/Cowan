"""calc_TP_M1E2: the merge, the parameter map, and one whole run.

The end-to-end test builds `WORK/Forbid` from its two input cards by
running RCN, RCN2, RCG and RCE into a temporary directory, then runs the
program over it.  It is the only reference case in this package with an
electric-quadrupole spectrum and a //R2// card in it; see that folder's
README for why it is inputs only and what it does not cover.
"""

import os
import shutil
import subprocess

import pytest

import calc_TP_M1E2
from calc_TP_M1E2 import Lev, Row, Transition
from conftest import CODE, WORK, case
from toolbox import ing11, parfile, radiate, rce, vacair


# --------------------------------------------------------------------------
# building rows by hand
# --------------------------------------------------------------------------

def level(index, energy, j=2.0, kset=1, observed=False):
    one = rce.Level(e_obs=energy, e_calc=energy, j=j, kset=kset, index=index,
                    flag='' if observed else '*')
    from toolbox.levels import Identified
    return Lev(Identified(level=one, kset=kset, j=j, rcg=index))


def transition(multipole, ga, cf=0.5, upper=None, lower=None, values=()):
    one = Transition(upper, lower, multipole)
    one.ga, one.cf = ga, cf
    one.s = radiate.s_from_ga(ga, 1000.0, multipole)
    for at, a in enumerate(values, 1):
        one.log_a.add(__import__('math').log(a))
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
    out = calc_TP_M1E2.merge(rows(100.0, 1.0), 0, 0.02)
    assert len(out) == 1
    assert out[0].kind == radiate.E2
    assert out[0].ga == pytest.approx(101.0)
    assert out[0].e2_fraction == pytest.approx(100.0 / 101.0)


def test_a_line_that_is_all_m1_keeps_the_name_m1_and_its_own_strength():
    out = calc_TP_M1E2.merge(rows(1.0, 100.0), 0, 0.02)
    assert len(out) == 1
    assert out[0].kind == radiate.M1
    assert out[0].ga == pytest.approx(100.0)
    assert out[0].e2_fraction == pytest.approx(1.0 / 101.0)


def test_a_mixed_line_is_called_m1_plus_e2_and_the_two_add():
    out = calc_TP_M1E2.merge(rows(50.0, 50.0), 0, 0.02)
    assert len(out) == 1
    assert out[0].kind == 'M1+E2'
    assert out[0].ga == pytest.approx(100.0)
    assert out[0].e2_fraction == pytest.approx(0.5)
    # A merged line has no single line strength, so none is quoted.
    assert out[0].s is None


def test_a_multipole_with_no_partner_goes_through_untouched():
    up, down = level(2, 100.0), level(1, 0.0)
    only = [Row(transition(radiate.E2, 7.0, upper=up, lower=down))]
    out = calc_TP_M1E2.merge(only, 0, 0.02)
    assert [r.kind for r in out] == [radiate.E2]
    other = [Row(transition(radiate.M1, 7.0, upper=up, lower=down))]
    assert [r.kind for r in calc_TP_M1E2.merge(other, 0, 0.02)] == \
        [radiate.M1]


def test_two_different_level_pairs_are_not_merged_with_each_other():
    a, b, c = level(1, 0.0), level(2, 100.0), level(3, 200.0)
    given = [Row(transition(radiate.E2, 1.0, upper=b, lower=a)),
             Row(transition(radiate.M1, 1.0, upper=c, lower=a))]
    out = calc_TP_M1E2.merge(given, 0, 0.02)
    assert len(out) == 2
    assert sorted(r.kind for r in out) == [radiate.E2, radiate.M1]


def test_the_trials_of_the_two_multipoles_add_trial_by_trial():
    up, down = level(2, 100.0), level(1, 0.0)
    given = [Row(transition(radiate.E2, 50.0, upper=up, lower=down,
                            values=(1.0, 2.0, 3.0))),
             Row(transition(radiate.M1, 50.0, upper=up, lower=down,
                            values=(10.0, 20.0, 30.0)))]
    out = calc_TP_M1E2.merge(given, 3, 0.02)
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
    out = calc_TP_M1E2.merge(given, 2, 0.02)
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
        cset, fitted, slots = calc_TP_M1E2.check_parameters(
            deck, csets, oute, kset)
        assert len(slots) == len(fitted) == len(cset.parameters)


def test_oute_and_rceout_are_matched_name_by_name():
    folder = case('Pr3')
    deck = ing11.read_ing11(folder=folder)
    csets = rce.read_rceout(folder=folder)
    oute = rce.read_oute(folder=folder)
    oute[1] = oute[1][:-1]              # one parameter short
    with pytest.raises(calc_TP_M1E2.RunError) as problem:
        calc_TP_M1E2.check_parameters(deck, csets, oute, 1)
    assert 'different runs' in str(problem.value)


def test_a_missing_oute_is_a_named_error_not_a_crash():
    folder = case('Pr3')
    deck = ing11.read_ing11(folder=folder)
    csets = rce.read_rceout(folder=folder)
    with pytest.raises(calc_TP_M1E2.RunError) as problem:
        calc_TP_M1E2.check_parameters(deck, csets, {}, 1)
    assert 'OUTE' in str(problem.value)


# --------------------------------------------------------------------------
# varying the parameters
# --------------------------------------------------------------------------

def options(**changes):
    out = {s.key: s.default for s in calc_TP_M1E2.SPECS}
    out.update(changes)
    return out


def parameter(name, flag, value, sigma):
    return rce.Parameter(name=name, flag=flag, value=value, sigma=sigma)


def test_a_parameter_with_no_standard_deviation_does_not_move():
    from toolbox import trials
    fitted = [parameter('EAV 4d7', -100, 21.7, 0.0),
              parameter('F2(11)', -100, 67.3, 0.0)]
    got = calc_TP_M1E2.vary_parameters(fitted, options(), trials.Randomizer(1),
                                       trials.Histogram('x'))
    assert got == [21.7, 67.3]


def test_a_group_moves_by_one_common_ratio():
    from toolbox import trials
    fitted = [parameter('F2(11)', -3, 100.0, 5.0),
              parameter('F4(11)', -3, 50.0, 2.5),
              parameter('F6(11)', -4, 20.0, 1.0)]
    got = calc_TP_M1E2.vary_parameters(fitted, options(),
                                       trials.Randomizer('g'),
                                       trials.Histogram('x'))
    assert got[0] / 100.0 == pytest.approx(got[1] / 50.0)
    assert got[2] / 20.0 != pytest.approx(got[0] / 100.0)


def test_without_grouping_every_parameter_moves_on_its_own():
    from toolbox import trials
    fitted = [parameter('F2(11)', -3, 100.0, 5.0),
              parameter('F4(11)', -3, 50.0, 2.5)]
    got = calc_TP_M1E2.vary_parameters(fitted,
                                       options(var_grouped_params=False),
                                       trials.Randomizer('g'),
                                       trials.Histogram('x'))
    assert got[0] / 100.0 != pytest.approx(got[1] / 50.0)


def test_a_fixed_parameter_moves_only_when_it_is_asked_to():
    # The branch `calc_TP_M1E2.pl:582` could never reach: it tests
    # abs($flag == 100), which is 0 for a flag of -100.
    from toolbox import trials
    fitted = [parameter('ALPHA', -100, 10.0, 0.0)]
    assert calc_TP_M1E2.vary_parameters(
        fitted, options(), trials.Randomizer(1), trials.Histogram('x')) \
        == [10.0]
    moved = calc_TP_M1E2.vary_parameters(
        fitted, options(fixed_param_variance=0.1), trials.Randomizer(1),
        trials.Histogram('x'))
    assert moved != [10.0]


def test_the_r_integrals_move_by_the_fraction_asked_for():
    from toolbox import trials
    folder = case('U90')
    deck = ing11.read_ing11(folder=folder)
    cards = ing11.read_r_integrals(deck, order=1)
    rng = trials.Randomizer('r')
    got = calc_TP_M1E2.vary_r_integrals(cards, options(r2_variance=0.15), rng,
                                        trials.Histogram('x'))
    moved = [(a, b.value) for a, b in zip(got, cards) if b.value]
    assert moved and all(a != b for a, b in moved)
    unmoved = calc_TP_M1E2.vary_r_integrals(cards, options(r2_variance=0.0),
                                            rng, trials.Histogram('x'))
    assert unmoved == [one.value for one in cards]


# --------------------------------------------------------------------------
# one whole run
# --------------------------------------------------------------------------

EXECUTABLES = ('RCN36K.exe', 'RCN2K.exe', 'rcg11kd.exe', 'RCE20KD.exe')


def build(folder, ispecc=None):
    """RCN, RCN2, RCG, RCE over WORK/Forbid's two input cards.

    `ispecc` overwrites column 72 of the IN2 main card, which chooses
    which of RCG's three line sorts are printed and so which layout the
    spectrum comes out in; the deck's own value is 2, the second-parity
    sort in the Lund layout.
    """
    source = WORK / 'Forbid'
    if not source.is_dir():
        pytest.skip('WORK/Forbid is not present')
    for name in EXECUTABLES:
        if not (CODE / name).is_file():
            pytest.skip('%s is not built' % name)
    folder.mkdir(parents=True, exist_ok=True)
    for name in ('IN2', 'IN36'):
        shutil.copyfile(source / name, folder / name)
    if ispecc is not None:
        deck = (folder / 'IN2').read_text().split('\n')
        card = deck[0].rstrip('\r')
        deck[0] = card[:71] + str(ispecc) + card[72:] + '\r'
        (folder / 'IN2').write_text('\n'.join(deck), newline='')
    (folder / 'COWAN.CFG').write_text(str(CODE) + os.sep + '\n')
    for name in EXECUTABLES:
        done = subprocess.run([str(CODE / name)], cwd=str(folder),
                              capture_output=True, text=True)
        assert done.returncode == 0, '%s: %s' % (name, done.stdout[-2000:])
    (folder / 'COWAN.CFG').unlink()
    return folder


def columns(path):
    return [line.split('\t') for line in
            path.read_text().rstrip('\n').split('\n')]


def test_a_whole_run_with_no_trials(tmp_path):
    folder = build(tmp_path)
    assert calc_TP_M1E2.main([str(folder)]) == 0
    table = columns(folder / 'tp_M1E2.txt')
    head, body = table[0], table[1:]
    assert head[:7] == ['type', 'conf1', 't1', 'J1', 'conf2', 't2', 'J2']
    assert body, 'the E2 spectrum produced no line'
    # The deck has DMIN = 0.00, so both multipoles print and the merge
    # is exercised on a real run: see WORK/Forbid/README.txt.
    assert {row[0] for row in body} == {radiate.M1, radiate.E2, 'M1+E2'}
    # Every level of this deck is 4d7, and none of them is observed.
    assert all(row[1] == '4d7' and row[4] == '4d7' for row in body)
    assert all(row[12] == '' and row[13] == '' for row in body)
    # Sorted by the vacuum wavelength the fitted energies give.  The
    # printed column is that wavelength in standard air between 2000 and
    # 20000 A, and in vacuum outside, as every line list quotes it.
    vacuum = [1e5 / abs(float(r[10]) - float(r[9])) for r in body]
    assert vacuum == sorted(vacuum)
    assert any(2000.0 < x < 20000.0 for x in vacuum), 'no air line in the case'
    for row, x in zip(body, vacuum):
        wanted = vacair.to_air(x) if 2000.0 < x < 20000.0 else x
        assert float(row[11]) == pytest.approx(wanted, abs=5e-5)


def test_a_wavelength_sorted_run_gives_the_same_table(tmp_path):
    """ISPECC = 4 prints the wide layout and nothing else.

    Column 72 of the main card chooses the sorts: 4 asks for the
    wavelength sort alone, and that one always prints the wide layout
    whatever ILUND says (`rcg11kd.f:14019`).  The same lines come out of
    it, to the extra digit the wider fields carry -- the wavelength sort
    prints F11.4 energies and an E11.4 gA where the Lund layout prints
    F11.3 and E10.3.
    """
    lund = build(tmp_path / 'lund', ispecc=2)
    wide = build(tmp_path / 'wide', ispecc=4)
    for folder in (lund, wide):
        assert calc_TP_M1E2.main([str(folder)]) == 0

    def by_pair(folder):
        table = columns(folder / 'tp_M1E2.txt')
        return {tuple(row[:7]): row for row in table[1:]}, table[0]

    one, head = by_pair(lund)
    two, head_again = by_pair(wide)
    assert head == head_again
    assert set(one) == set(two)
    assert one, 'the run produced no line'
    for key, row in one.items():
        assert float(row[15]) == pytest.approx(float(two[key][15]), rel=2e-3)


def test_a_whole_run_with_trials_is_repeatable(tmp_path):
    folder = build(tmp_path)
    parfile.write(folder / calc_TP_M1E2.PAR_FILE, calc_TP_M1E2.SPECS,
                  {'num_trials': 2, 'store_trials': 2, 'random_seed': 'test'})
    assert calc_TP_M1E2.main([str(folder)]) == 0
    first = (folder / 'tp_M1E2.txt').read_text()
    head = columns(folder / 'tp_M1E2.txt')[0]
    assert 'std_A%' in head and 'Trial A data' in head
    assert (folder / 'trial1' / 'OUTG11').is_file()
    assert not (folder / 'trial1' / 'COWAN.CFG').exists()

    # The same seed over the trials that are already there: the numbers
    # come back the same, and RCG is not run again.
    (folder / calc_TP_M1E2.PAR_FILE).unlink()
    parfile.write(folder / calc_TP_M1E2.PAR_FILE, calc_TP_M1E2.SPECS,
                  {'num_trials': 2, 'store_trials': 2, 'random_seed': 'test',
                   'read_only': 2})
    assert calc_TP_M1E2.main([str(folder)]) == 0
    assert (folder / 'tp_M1E2.txt').read_text() == first


def test_the_trial_decks_differ_from_the_original_only_in_the_r2_cards(
        tmp_path):
    folder = build(tmp_path)
    parfile.write(folder / calc_TP_M1E2.PAR_FILE, calc_TP_M1E2.SPECS,
                  {'num_trials': 1, 'store_trials': 1, 'random_seed': 'deck'})
    assert calc_TP_M1E2.main([str(folder)]) == 0
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
    monkeypatch.setattr(calc_TP_M1E2, 'PAR_FILE', 'absent_here.par')
    assert calc_TP_M1E2.main([str(tmp_path)]) == 0
    got = parfile.read(tmp_path / 'absent_here.par', calc_TP_M1E2.SPECS)
    assert got == {spec.key: spec.default for spec in calc_TP_M1E2.SPECS}


def test_the_program_writes_a_parameter_file_on_request(tmp_path):
    assert calc_TP_M1E2.main([str(tmp_path), '--write-par']) == 0
    got = parfile.read(tmp_path / calc_TP_M1E2.PAR_FILE, calc_TP_M1E2.SPECS)
    assert got == {s.key: s.default for s in calc_TP_M1E2.SPECS}
