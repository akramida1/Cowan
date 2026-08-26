"""The relative-energy groups of `FOR/RCE20KD.F`, driven end to end.

`CODE/tests/test_rce.py` covers the Python side of the feature -- the
status model, the parsers, the group card block as `toolbox/rce.py`
reads and writes it.  This file covers the Fortran: it compiles
`RCE20KD.F` out of the repository, builds `TAPE2E` for `WORK/U90` by
running RCN, RCN2 and RCG once, and then drives RCE over a matrix of
decks that differ from `WORK/U90/RCEINP.fin` by a letter in column 14
and a group card block.

Why U90.  It is the smallest real fit in the package -- ten fitted
levels in configuration set 1, nine in set 2, five J blocks, eleven
parameters -- and the whole chain runs in about a second, so every case
below can afford a run of its own.  It also has the shape the feature
was designed for: one configuration (2p2) whose levels spread over
J=1/2, 3/2 and 5/2, so a group can span J blocks, and a second
configuration set that carries no letter at all and must therefore come
out untouched.  Its converged fit is `AVDEV= 191.931` for set 1 and
`AVDEV= 190.618` for set 2.

The numbers asserted below are not transcribed from a run of the code
being tested.  Two of them are independent:

  * the weighted residual sum is recomputed here, in Python, out of the
    `T-W` column RCE prints and the weights it reports, and `AVDEV` and
    `STDEV` are derived from it -- so `test_the_residual_sum_gives_avdev
    _and_stdev` would fail if the differencing, the `sqrt(w)` scaling or
    the degrees of freedom moved, without knowing what RCE was going to
    print;
  * `AVDEV= 190.618` for configuration set 2 is the value the legacy
    deck converges to and is fixed by `WORK/U90/README.txt`.

Everything else is a structural claim -- a stop happens, a warning does
not stop, a letter comes back on the record it went in on -- and is
checked as such.

The two skips.  There is no test here without gfortran, because the
point of the file is the Fortran; and there is none without the cfp
decks, which `conftest.cowan_code` builds.  Both are the same skips the
rest of the suite already takes.
"""

import hashlib
import math
import os
import platform
import re
import shutil
import subprocess
import sys

import pytest

from conftest import CODE, EXECUTABLES, REPO, WORK

SOURCE = REPO / 'FOR' / 'RCE20KD.F'
DECK = WORK / 'U90'

#: CLAUDE.md section 3's build line, as `CODE/unix/build.sh` spells it:
#: the same flags and the same -O3, with -malign-double added only on
#: x86 because gfortran rejects it on Apple Silicon and on ARM Linux.
#: -fopenmp is left off -- no parallel region here reassociates a sum, so
#: a serial build gives the same numbers, and this way the test does not
#: need libgomp beside it.  -fno-automatic stays regardless: it is what
#: the arithmetic of this package was accepted with.
FLAGS = ('-std=legacy', '-fshort-enums', '-ftracer', '-fno-backslash',
         '-fallow-argument-mismatch', '-fcray-pointer',
         '-fallow-invalid-boz', '-fno-automatic', '-O3')
if platform.machine().lower() not in ('arm64', 'aarch64'):
    FLAGS += ('-malign-double',)

#: The header of the group card block, spelled as `RELWR` writes it
#: (`FOR/RCE20KD.F`, FORMAT 10 of `RELWR`) and as `RELRD` recognises it.
HEADER = ' GROUP SIGMA_BASE  SIGMA_INT'

#: Line numbers in `WORK/U90/RCEINP.fin`, 1-based, of the ten fitted
#: levels of configuration set 1.  Written out rather than searched for
#: because every case below is "this deck with these lines lettered",
#: and a name is easier to read at the call site than a number.
GROUND = 4                      # 2s,  J=1/2, the level at 0
L_2S2 = 5                       # 2s2, J=1/2
HALF = (6, 7, 8)                # 2p2, J=1/2
THREE_HALVES = (9, 10, 11)      # 2p2, J=3/2
FIVE_HALVES = (12, 13)          # 2p2, J=5/2

#: Where set 1's parameter block ends: set 2's control card is line 26,
#: so a group card block for set 1 goes immediately before it.
SET2 = 26


# --------------------------------------------------------------------------
# the deck generator
# --------------------------------------------------------------------------

def deck(letters=(), cards=(), cards2=(), edits=()):
    """`WORK/U90/RCEINP.fin` with letters in column 14 and group cards.

    `letters` maps a 1-based line number of the file to the character
    that is to go in column 14 of that record.  `cards` are the group
    cards of configuration set 1, placed with their header immediately
    before set 2's control card, which is where `RELRD` looks for them
    and where `RELWR` puts them; `cards2` are set 2's, at end of file.
    `edits` is a sequence of (old, new) byte pairs applied to the whole
    deck afterwards, for the cases that need a parameter moved.

    Everything is done on bytes and the CRLF endings of the original are
    kept, because column 14 is a column and this file is the input whose
    columns are being tested.
    """
    lines = DECK.joinpath('RCEINP.fin').read_bytes().split(b'\r\n')
    if lines and lines[-1] == b'':
        lines = lines[:-1]
    for number, letter in dict(letters).items():
        line = lines[number - 1].ljust(14)
        assert line[13:14] == b' ', 'line %d already has a status' % number
        lines[number - 1] = line[:13] + letter.encode() + line[14:]
    out = []
    for number, line in enumerate(lines, 1):
        if number == SET2 and cards:
            out.append(HEADER.encode())
            out.extend(one.encode() for one in cards)
        out.append(line)
    if cards2:
        out.append(HEADER.encode())
        out.extend(one.encode() for one in cards2)
    text = b'\r\n'.join(out) + b'\r\n'
    for old, new in edits:
        assert old in text, old
        text = text.replace(old, new)
    return text


def card(letter, sigma_base, sigma_int):
    """One group card: letter in column 1, the sigmas in 2-13 and 14-25.

    `sigma_base=None` is the untied spelling -- a blank field, which is
    not the same as a zero and is rejected as one if written that way.
    """
    base = '' if sigma_base is None else '%12.3f' % sigma_base
    return '%s%-12s%12.3f' % (letter, base, sigma_int)


# --------------------------------------------------------------------------
# fixtures
# --------------------------------------------------------------------------

@pytest.fixture(scope='module')
def rce(tmp_path_factory):
    """RCE built from the repository, not the executable shipped in CODE/.

    `CODE/RCE20KD.exe` is a committed binary and is only replaced at a
    release, so through the development of a feature it is the *old*
    program.  Testing the source means compiling the source.
    """
    if not SOURCE.is_file():
        pytest.skip('%s is not present' % SOURCE)
    folder = tmp_path_factory.mktemp('rcebuild')
    program = folder / ('rce.exe' if sys.platform == 'win32' else 'rce')
    build = ['gfortran'] + list(FLAGS) + ['-o', str(program), str(SOURCE)]
    try:
        done = subprocess.run(build, cwd=str(REPO / 'FOR'),
                              capture_output=True, text=True)
    except OSError as trouble:
        pytest.skip('gfortran is not available: %s' % trouble)
    if done.returncode or not program.is_file():
        pytest.fail('RCE20KD.F does not compile:\n%s'
                    % (done.stderr or done.stdout))
    return program


@pytest.fixture(scope='module')
def tape2e(tmp_path_factory, cowan_code):
    """`TAPE2E` for the U90 deck, from one run of RCN, RCN2 and RCG.

    RCE reads the angular data from `TAPE2E` and its run control from
    the first line of `OUTGINE`; the deck ships the converged `OUTGINE`
    and `RCEINP` but not `TAPE2E`, which is generated (CLAUDE.md
    section 7).  This builds it once for the whole module.
    """
    if not DECK.is_dir():
        pytest.skip('WORK/U90 is not present')
    folder = tmp_path_factory.mktemp('u90')
    for name in ('IN2', 'IN36'):
        shutil.copyfile(DECK / name, folder / name)
    (folder / 'COWAN.CFG').write_text(str(cowan_code) + os.sep + '\n')
    for name in EXECUTABLES:
        done = subprocess.run([str(CODE / name)], cwd=str(folder),
                              capture_output=True, text=True)
        assert done.returncode == 0, '%s: %s' % (name, done.stdout[-2000:])
    (folder / 'COWAN.CFG').unlink()
    made = folder / 'TAPE2E'
    assert made.is_file() and made.stat().st_size, 'RCG wrote no TAPE2E'
    return made


class Run(object):
    """One run of RCE, and the files it left behind."""

    def __init__(self, folder, done):
        self.folder = folder
        self.stdout = done.stdout or ''
        self.stderr = done.stderr or ''
        self.oute = (folder / 'OUTE').read_text(encoding='latin-1')
        out = folder / 'RCEOUT'
        self.rceout = out.read_bytes() if out.is_file() else b''

    @property
    def message(self):
        """Everything RCE said, for a test that wants to match on it."""
        return self.stdout + self.stderr


def run(rce, tape2e, folder, text, outgine=None):
    """Run RCE on `text` in `folder`, with U90's converged `OUTGINE`.

    `OUTGINE.fin` carries `INVECT=1` in columns 21-25 of its first line,
    which is what makes RCE take the levels and the parameters from
    `RCEINP` instead of writing a fresh one (`WORK/U90/README.txt`).
    """
    folder.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(str(tape2e), str(folder / 'TAPE2E'))
    if outgine is None:
        outgine = DECK.joinpath('OUTGINE.fin').read_bytes()
    (folder / 'OUTGINE').write_bytes(outgine)
    (folder / 'RCEINP').write_bytes(text)
    done = subprocess.run([str(rce)], cwd=str(folder),
                          capture_output=True, text=True)
    return Run(folder, done)


# --------------------------------------------------------------------------
# readers for what RCE prints
# --------------------------------------------------------------------------

#: ` GROUP  LEVELS    BASE  J BLOCK    SIGMA_BASE     SIGMA_INT`
SUMMARY = re.compile(r'^ +([A-Z]) +(\d+) +(\d+) +J= *([-\d.]+) +'
                     r'(\(UNTIED\)|[\d.]+) +([\d.]+) *$')
#: ` RELATIVE-ENERGY GROUPS - CYCLE  3:  STDEV=  383.172  SIGMA_0=  383.172`
CYCLE = re.compile(r'RELATIVE-ENERGY GROUPS - CYCLE +(\d+): +'
                   r'STDEV= *([-\d.]+) +SIGMA_0= *([-\d.]+)(.*)$')
#: `     X    0.325657    0.996698`
WEIGHTS = re.compile(r'^ +([A-Z]) +([\d.]+) +([\d.]+) *$')
#: `        10. LEVELS FIT WITH   7. PARAMETERS`
COUNTS = re.compile(r'([\d.]+) LEVELS FIT WITH +([\d.]+) PARAMETERS')
AVDEV = re.compile(r'Iteration No +(\d+) finished with AVDEV=(.+)$')


def sections(text):
    """`OUTE` split into one string per configuration set.

    RCE prints ` Conf. set No n` when it starts one, and U90 has two:
    the set that carries the letters and the set that must not notice
    them.  Every reader below is happier given one of these than the
    whole file, because the last `CYC NO=` line and the last `T-W`
    table of the file belong to set 2.
    """
    parts = re.split(r'(?= Conf\. set No )', text)
    return [one for one in parts if ' Conf. set No ' in one]


def summary(text):
    """The group summary rows, as {letter: (levels, base, j, sb, si)}."""
    rows = {}
    for line in text.split('\n'):
        found = SUMMARY.match(line.rstrip('\r'))
        if found:
            base = (None if found.group(5) == '(UNTIED)'
                    else float(found.group(5)))
            rows[found.group(1)] = (int(found.group(2)), int(found.group(3)),
                                    float(found.group(4)), base,
                                    float(found.group(6)))
    return rows


def cycles(text):
    """One entry per printed cycle: (n, stdev, sigma_0, frozen, weights)."""
    out = []
    lines = text.split('\n')
    for number, line in enumerate(lines):
        found = CYCLE.search(line)
        if not found:
            continue
        weights = {}
        for after in lines[number + 1:number + 30]:
            row = WEIGHTS.match(after.rstrip('\r'))
            if row:
                weights[row.group(1)] = (float(row.group(2)),
                                         float(row.group(3)))
            elif weights:
                break
        out.append((int(found.group(1)), float(found.group(2)),
                    float(found.group(3)),
                    'FROZEN' in found.group(4), weights))
    return out


def avdevs(text):
    """The AVDEV of every iteration, split into one list per set.

    A new configuration set starts at its own `Iteration No 0`, and the
    U90 deck has two, so this is what separates "the set that carries
    the group" from "the set that must be untouched".
    """
    sets = []
    for line in text.split('\n'):
        found = AVDEV.search(line)
        if not found:
            continue
        if found.group(1) == '0':
            sets.append([])
        value = found.group(2).strip()
        sets[-1].append(None if '*' in value else float(value))
    return sets


def counts(text):
    """(levels fit, parameters) of every `LEVELS FIT WITH` line."""
    return [(float(one), float(two)) for one, two in COUNTS.findall(text)]


def residuals(text, which=0):
    """The last `T-W` table: [(T, letter, W, T-W)] in printed order.

    `FOR/RCE20KD.F` writes it as T in columns 1-11, the status character
    in 12, W in 13-22 and T-W in 23-31.  The letter is there because
    B3 gave every write site of the status character a source of its
    own; it is what lets this reader tell a base from a member without
    re-reading the deck.
    """
    blocks = [n for n, line in enumerate(text.split('\n'))
              if 'T-W    GKK' in line]
    if not blocks:
        return []
    out = []
    for line in text.split('\n')[blocks[-1] + 1:]:
        line = line.rstrip('\r')
        if 'Finished output' in line or 'T-W    GKK' in line:
            break
        try:
            energy = float(line[0:11])
            calculated = float(line[12:22])
            difference = float(line[22:31])
        except (ValueError, IndexError):
            continue
        out.append((energy, line[11], calculated, difference))
    return out


def digest(path):
    return hashlib.md5(path.read_bytes()).hexdigest()


def unfinished(got):
    """True when the run stopped before it wrote a fit.

    `RCEOUT` is not a witness on its own: RCE opens unit 21 with the
    other units before it has read a card, so a stop in `RELCHK` still
    leaves an empty file behind.  What a finished run always writes is
    the parameter block -- and, when the deck has a group, the group
    card block after it.
    """
    out = got.folder / 'RCEOUT'
    return not out.is_file() or b'PARAMETER FLAG' not in out.read_bytes()


# --------------------------------------------------------------------------
# the decks the cases below share
# --------------------------------------------------------------------------

#: One base and two members, all in J=1/2, sigma_base 500, sigma_int 20.
#: The smallest case with an internal difference in it, and the one the
#: residual sum is recomputed for.
SMALL = dict(letters={HALF[1]: 'X', HALF[0]: 'x', HALF[2]: 'x'},
             cards=[card('X', 500.0, 20.0)])

#: The same base with every other 2p2 level of the set as a member:
#: seven levels over three J blocks.  This is the shape the feature was
#: written for -- a quasi-isolated system tied to ground by one badly
#: measured interval.
TIED = dict(letters=dict([(HALF[1], 'X')]
                         + [(one, 'x') for one in
                            (HALF[0], HALF[2]) + THREE_HALVES + FIVE_HALVES]),
            cards=[card('X', 500.0, 20.0)])


# --------------------------------------------------------------------------
# a deck with no letters is the deck it always was
# --------------------------------------------------------------------------

def shipped_rce_has_the_feature():
    """True once `CODE/RCE20KD.exe` is a build that knows about groups.

    The committed binary is the pre-feature program until a release
    replaces it, which makes it an oracle for "nothing changed".  After
    that release it is not one any more, and the comparison below has to
    say so rather than quietly compare a binary with itself.
    """
    program = CODE / 'RCE20KD.exe'
    if not program.is_file():
        return None
    return b'SIGMA_BASE' in program.read_bytes()


def test_a_deck_with_no_letters_says_nothing_about_groups(rce, tape2e,
                                                          tmp_path):
    """The whole feature is silent unless a letter asks for it.  This is
    the additivity clause of CLAUDE.md section 1 in its weakest and most
    durable form: no group card block, no summary, no per-cycle report,
    and set 1 back at its documented `AVDEV`."""
    got = run(rce, tape2e, tmp_path / 'legacy', deck())
    assert 'RELATIVE-ENERGY' not in got.oute
    assert HEADER not in got.oute
    assert HEADER.encode() not in got.rceout
    assert summary(got.oute) == {}
    first, second = avdevs(got.oute)
    assert first[-1] == 191.931
    assert second[-1] == 190.618


def test_a_deck_with_no_letters_reproduces_the_shipped_binary(rce, tape2e,
                                                              tmp_path):
    """Byte for byte, over every file RCE writes but `OUTE`.

    `OUTE` carries timing lines and is excluded for that reason alone;
    everything else -- `RCEOUT`, `PARVALS`, `LEVELS1..3`, `RCEINP.HF` --
    has to be identical, because a legacy deck is every deck in the
    repository and none of them may move."""
    has = shipped_rce_has_the_feature()
    if has is None:
        pytest.skip('CODE/RCE20KD.exe is not there to compare against')
    if has:
        pytest.skip('CODE/RCE20KD.exe already carries the feature, so it is '
                    'no longer an independent oracle for the old behaviour')
    text = deck()
    new = run(rce, tape2e, tmp_path / 'new', text)
    old = run(CODE / 'RCE20KD.exe', tape2e, tmp_path / 'old', text)
    for name in ('RCEOUT', 'PARVALS', 'LEVELS1', 'LEVELS2', 'LEVELS3',
                 'RCEINP.HF'):
        one, two = new.folder / name, old.folder / name
        assert one.is_file() == two.is_file(), name
        if one.is_file():
            assert digest(one) == digest(two), name


# --------------------------------------------------------------------------
# the group summary and the group card block
# --------------------------------------------------------------------------

def test_the_summary_names_the_base_its_block_and_both_sigmas(rce, tape2e,
                                                              tmp_path):
    got = run(rce, tape2e, tmp_path / 'small', deck(**SMALL))
    rows = summary(got.oute)
    assert list(rows) == ['X']
    levels, base, j, sigma_base, sigma_int = rows['X']
    assert (levels, j, sigma_base, sigma_int) == (3, 0.5, 500.0, 20.0)
    # the base is the second 2p2 level of the J=1/2 block, which is the
    # fourth level of that block counting the 2s and 2s2 levels below it
    assert base == 4


def test_a_sigma_without_a_decimal_point_is_read_as_units(rce, tape2e,
                                                          tmp_path):
    """`RELRD` reads the two fields list-directed out of their columns
    rather than under F12.3, so `500` means 500 and not 0.5.  The test
    is that the two spellings give the same run, not merely the same
    printed sigma."""
    plain = run(rce, tape2e, tmp_path / 'plain', deck(**SMALL))
    bare = run(rce, tape2e, tmp_path / 'bare',
               deck(letters=SMALL['letters'],
                    cards=['X         500          20']))
    assert summary(bare.oute) == summary(plain.oute)
    assert cycles(bare.oute) == cycles(plain.oute)
    assert avdevs(bare.oute) == avdevs(plain.oute)


def test_the_group_block_of_one_set_does_not_reach_the_next(rce, tape2e,
                                                            tmp_path):
    """Set 2 of U90 carries no letter.  `RELRD` reads set 1's block, hits
    set 2's control card, backspaces over it and leaves it for the read
    that expects it -- so set 2 has to converge to its legacy value."""
    got = run(rce, tape2e, tmp_path / 'tied', deck(**TIED))
    assert avdevs(got.oute)[1][-1] == 190.618


def test_a_group_may_span_configurations_and_j_blocks(rce, tape2e, tmp_path):
    """Two groups at once: X ties a 2s2 level to a 2p2 one, so its
    internal difference mixes derivative rows belonging to two different
    `EAV`s; Y lives entirely in J=3/2.  Both are reported, both are
    weighted, and the run converges."""
    got = run(rce, tape2e, tmp_path / 'two',
              deck(letters={L_2S2: 'X', HALF[0]: 'x',
                            THREE_HALVES[1]: 'Y', THREE_HALVES[2]: 'y'},
                   cards=[card('X', 500.0, 20.0), card('Y', 2000.0, 100.0)]))
    rows = summary(got.oute)
    assert sorted(rows) == ['X', 'Y']
    assert rows['X'][0] == rows['Y'][0] == 2
    assert rows['X'][2] == 0.5 and rows['Y'][2] == 1.5
    last = cycles(got.oute)[-1][4]
    # the wider group is the more heavily downweighted one, in both
    # halves of its weight pair
    assert last['Y'][0] < last['X'][0]
    assert last['Y'][1] < last['X'][1]
    settled = avdevs(got.oute)[0]
    assert settled[-1] == pytest.approx(settled[-2], abs=0.05)


# --------------------------------------------------------------------------
# the arithmetic
# --------------------------------------------------------------------------

def test_the_residual_sum_gives_avdev_and_stdev(rce, tape2e, tmp_path):
    """The known-answer case, recomputed here rather than transcribed.

    `03_design.md` section 4 says the mechanism is: difference each
    member's row against its base's, scale every row of the group by
    `sqrt(w)`, and leave the blank levels alone.  So the residual that
    reaches `R` is `w_base * r_base**2` for a base, `w_int * (r_i -
    r_base)**2` for a member and `r**2` for everything else -- and
    `AVDEV = sqrt(R/TNO)`, `STDEV = sqrt(R/(TNO-PNO))`.

    All four inputs come off the page: the `T-W` column, the weights of
    the last cycle, and the `LEVELS FIT WITH` counts.  Nothing here
    knows what `AVDEV` was going to be."""
    got = run(rce, tape2e, tmp_path / 'small', deck(**SMALL))
    first = sections(got.oute)[0]
    weight = cycles(first)[-1][4]['X']
    rows = residuals(first)
    assert len(rows) == 10
    base = [one for one in rows if one[1] == 'X']
    assert len(base) == 1
    offset = base[0][3]
    total = 0.0
    for _, letter, _, difference in rows:
        if letter == 'X':
            total += weight[0] * difference ** 2
        elif letter == 'x':
            total += weight[1] * (difference - offset) ** 2
        else:
            total += difference ** 2
    levels, parameters = counts(first)[0]
    assert (levels, parameters) == (10.0, 7.0)
    # the T-W column is printed to three decimals, which is what sets
    # the tolerance: 1e-3 on a residual of 4e2 is 5e-6 relative
    assert (math.sqrt(total / levels)
            == pytest.approx(avdevs(got.oute)[0][-1], abs=0.01))
    printed = float(re.search(r'STDEV= *([\d.]+)',
                              first.split('CYC NO=')[-1]).group(1))
    assert (math.sqrt(total / (levels - parameters))
            == pytest.approx(printed, abs=0.02))


def test_sigma_0_falls_monotonically_and_then_freezes(rce, tape2e, tmp_path):
    """`sigma_0` is the deviation the model itself achieves, and it is
    re-estimated every cycle from `STDEV` under a ratchet -- it never
    grows -- and frozen once it stops moving by more than one per cent
    (`03_design.md` section 2.2).  Both halves are visible on the page,
    which is why they are printed there."""
    got = run(rce, tape2e, tmp_path / 'tied', deck(**TIED))
    seen = [one for one in cycles(got.oute) if one[0] is not None]
    adopted = [one[2] for one in seen]
    assert adopted == sorted(adopted, reverse=True)
    assert seen[0][3] is False
    assert seen[-1][3] is True
    # the ratchet, on the page: the raw STDEV of the last cycles is
    # below the sigma_0 that is still being used
    assert seen[-1][1] <= seen[-1][2]
    # and the weights settle with it
    assert seen[-1][4]['X'] == seen[-2][4]['X']


def test_a_tighter_internal_sigma_weighs_the_differences_more(rce, tape2e,
                                                              tmp_path):
    """`w = sigma_0**2/(sigma_0**2 + sigma**2)` is monotone decreasing in
    sigma, so halving `sigma_int` has to raise `W_INT`, and the fit has
    to follow the internal differences more closely."""
    loose = run(rce, tape2e, tmp_path / 'loose',
                deck(letters=TIED['letters'], cards=[card('X', 500.0, 200.0)]))
    tight = run(rce, tape2e, tmp_path / 'tight',
                deck(letters=TIED['letters'], cards=[card('X', 500.0, 20.0)]))
    assert (cycles(tight.oute)[-1][4]['X'][1]
            > cycles(loose.oute)[-1][4]['X'][1])


# --------------------------------------------------------------------------
# untied groups
# --------------------------------------------------------------------------

def test_a_blank_sigma_base_is_an_untied_group(rce, tape2e, tmp_path):
    """Blank is not zero.  Blank means the offset is unknown rather than
    poorly known, so the base contributes no residual at all -- `w = 0`
    exactly, in every cycle -- and the summary says `(UNTIED)` where the
    number would be."""
    got = run(rce, tape2e, tmp_path / 'untied',
              deck(letters=SMALL['letters'], cards=[card('X', None, 20.0)]))
    assert summary(got.oute)['X'][3] is None
    assert all(one[4]['X'][0] == 0.0 for one in cycles(got.oute))
    assert all(one[4]['X'][1] > 0.0 for one in cycles(got.oute))


def test_an_untied_group_gives_up_a_degree_of_freedom(rce, tape2e, tmp_path):
    """Its base's row is identically zero, so it is not a constraint and
    must not be counted as one: `TNO` drops by one per untied group."""
    tied = run(rce, tape2e, tmp_path / 'tied',
               deck(letters=SMALL['letters'], cards=[card('X', 500.0, 20.0)]))
    free = run(rce, tape2e, tmp_path / 'free',
               deck(letters=SMALL['letters'], cards=[card('X', None, 20.0)]))
    assert counts(tied.oute)[0] == (10.0, 7.0)
    assert counts(free.oute)[0] == (9.0, 7.0)


def test_an_untied_group_that_frees_a_parameter_stops(rce, tape2e, tmp_path):
    """The gauge test of `03_design.md` section 2.4.  With every fitted
    level of set 1 in an untied group, nothing holds `EAV 2s`: the
    parameter and the group's offset move together with no effect on the
    fit, and the normal equations are singular in that direction.  RCE
    stops, names the configuration and states the three remedies."""
    letters = dict([(GROUND, 'X')]
                   + [(one, 'x') for one in
                      (L_2S2,) + HALF + THREE_HALVES + FIVE_HALVES])
    got = run(rce, tape2e, tmp_path / 'gauge2',
              deck(letters=letters, cards=[card('X', None, 20.0)]))
    assert 'NOTHING FIXES THE POSITION OF CONFIGURATION 2s' in got.message
    assert 'UNTIED GROUP LEAVES A PARAMETER FREE' in got.message
    assert unfinished(got)


def test_an_untied_group_that_leaves_a_level_outside_it_runs(rce, tape2e,
                                                             tmp_path):
    """The negative half of the same rule, and the reason it is a
    counting test and not a blanket prohibition: put every 2p2 level in
    the untied group but leave the 2s2 level blank, and that level's
    2p2 components still tie the configuration to an absolute energy.
    Nothing is free, so nothing stops."""
    letters = dict([(HALF[0], 'X')]
                   + [(one, 'x') for one in
                      HALF[1:] + THREE_HALVES + FIVE_HALVES])
    got = run(rce, tape2e, tmp_path / 'gauge',
              deck(letters=letters, cards=[card('X', None, 20.0)]))
    assert 'LEAVES A PARAMETER FREE' not in got.message
    assert summary(got.oute)['X'][0] == 8
    assert avdevs(got.oute)[1][-1] == 190.618


# --------------------------------------------------------------------------
# the malformed-input set
# --------------------------------------------------------------------------

#: Each case is the deck with one thing wrong, the phrase RCE must print
#: about it, and the phrase it must stop with.  The stop text matters as
#: much as the message: `gfortran`'s `STOP 'string'` exits 0, so a script
#: driving RCE can only tell these apart by what was said.
BAD = [
    ('a letter with no group card',
     dict(letters=SMALL['letters']),
     'LEVELS CARRY THE LETTER X, BUT THERE IS NO GROUP CARD FOR IT',
     'GROUP WITHOUT A CARD IN RCEINP'),
    ('a group card no level carries',
     dict(letters=SMALL['letters'],
          cards=[card('X', 500.0, 20.0), card('Z', 100.0, 10.0)]),
     'A GROUP CARD WAS GIVEN FOR GROUP Z, BUT NO LEVEL RECORD CARRIES',
     'GROUP CARD WITHOUT LEVELS IN RCEINP'),
    ('members with no base',
     dict(letters={HALF[0]: 'x', HALF[2]: 'x'},
          cards=[card('X', 500.0, 20.0)]),
     'GROUP X HAS     2 MEMBERS AND NO BASE',
     'GROUP WITHOUT A BASE IN RCEINP'),
    ('two bases for one letter',
     dict(letters={HALF[0]: 'X', HALF[1]: 'X', HALF[2]: 'x'},
          cards=[card('X', 500.0, 20.0)]),
     'GROUP X HAS TWO BASES',
     'TWO BASES FOR ONE GROUP IN RCEINP'),
    ('an explicit zero sigma_base',
     dict(letters=SMALL['letters'], cards=[card('X', 0.0, 20.0)]),
     'IT MUST BE POSITIVE.  LEAVE THE FIELD BLANK',
     'BAD GROUP CARD IN RCEINP'),
    ('a negative sigma_int',
     dict(letters=SMALL['letters'], cards=[card('X', 500.0, -20.0)]),
     'SIGMA_INT OF GROUP X IS',
     'BAD GROUP CARD IN RCEINP'),
    ('two cards for one group',
     dict(letters=SMALL['letters'],
          cards=[card('X', 500.0, 20.0), card('X', 100.0, 10.0)]),
     'TWO GROUP CARDS FOR GROUP X',
     'DUPLICATE GROUP CARD IN RCEINP'),
    ('an unrecognised status character',
     dict(letters={HALF[1]: '?'}),
     'UNRECOGNIZED STATUS CHARACTER "?" IN COLUMN 14',
     'BAD STATUS CHARACTER IN RCEINP'),
]


@pytest.mark.parametrize('name,spec,said,stopped',
                         BAD, ids=[one[0] for one in BAD])
def test_a_malformed_group_stops_and_says_why(rce, tape2e, tmp_path,
                                              name, spec, said, stopped):
    got = run(rce, tape2e, tmp_path / 'bad', deck(**spec))
    assert said in got.message, got.message[-1500:]
    assert stopped in got.message
    assert unfinished(got)


def test_a_group_of_one_level_is_a_warning_and_not_a_stop(rce, tape2e,
                                                          tmp_path):
    """Legal but pointless: a base alone is one downweighted absolute
    level and no internal difference at all.  RCE says so and runs."""
    got = run(rce, tape2e, tmp_path / 'singleton',
              deck(letters={HALF[1]: 'X'}, cards=[card('X', 500.0, 20.0)]))
    assert 'WARNING: GROUP X IS ITS BASE ALONE' in got.message
    assert summary(got.oute)['X'][0] == 1
    assert avdevs(got.oute)[1][-1] == 190.618


# --------------------------------------------------------------------------
# the round trip
# --------------------------------------------------------------------------

def test_rceout_carries_the_letters_and_the_group_block(rce, tape2e,
                                                        tmp_path):
    """`RCEOUT` copied onto `RCEINP` is how this package is worked, so
    everything the fit needs has to survive the copy: the letters in
    column 14, which are regenerated from `IGRP` and not carried, and
    the group card block, which goes back where `RELRD` looks for it --
    after the parameter block, before the next set's control card."""
    got = run(rce, tape2e, tmp_path / 'small', deck(**SMALL))
    lines = got.rceout.split(b'\r\n')
    lettered = {line[:13].strip(): line[13:14]
                for line in lines if len(line) > 14 and line[13:14].strip()}
    assert lettered.get(b'810836.077000') == b'X'
    assert lettered.get(b'776509.696000') == b'x'
    assert lettered.get(b'846376.043000') == b'x'
    header = lines.index(HEADER.encode())
    assert lines[header + 1].startswith(b'X')
    assert lines[header + 2].startswith(b'    2')      # set 2's control card


def test_the_round_trip_deck_fits_to_the_same_answer(rce, tape2e, tmp_path):
    """The whole point of writing the block: feed `RCEOUT` back in and
    the same problem is fitted again, with the same groups and the same
    weights.  Only the parameter values have moved, so the fit starts at
    its own answer and stays there."""
    first = run(rce, tape2e, tmp_path / 'first', deck(**SMALL))
    again = run(rce, tape2e, tmp_path / 'again', first.rceout)
    assert summary(again.oute) == summary(first.oute)
    assert avdevs(again.oute)[0][-1] == avdevs(first.oute)[0][-1]
    # `sigma_0` is deliberately not written to the deck -- it is derived,
    # not data (`RELWR`) -- so the second run starts again at weight 1 and
    # re-anneals.  It arrives at the same fixed point, not at the same
    # bytes: the agreement is to the last few figures, not exact.
    was = cycles(sections(first.oute)[0])[-1][4]['X']
    now = cycles(sections(again.oute)[0])[-1][4]['X']
    assert now == pytest.approx(was, abs=1e-3)


# --------------------------------------------------------------------------
# identity under rematching
# --------------------------------------------------------------------------

def test_a_letter_follows_its_level_through_the_rematching(rce, tape2e,
                                                           tmp_path):
    """`03_design.md` section 6: the group lives on the `(K,M)` slot,
    which `NF` fixes when the deck is read and which the reordering
    never permutes -- not on an eigenvector, which strong CI can hand to
    a different level between one cycle and the next.

    The J=1/2 block of U90 is where that matters: 776509.696 is 63 per
    cent 2p2 (3P) 4P against 51 per cent (1S) 2S, and its separation
    from 810836.077 is set by `ZETA 3`.  Halving `ZETA 3` moves the
    calculated block far enough that cycle 0 cannot even print an
    `AVDEV`, so the recognition has real work to do on the way back.
    The letters must still come out on the levels they went in on, and
    the fit must still reach the answer the unperturbed deck reaches."""
    steady = run(rce, tape2e, tmp_path / 'steady', deck(**SMALL))
    shaken = run(rce, tape2e, tmp_path / 'shaken',
                 deck(edits=[(b'ZETA 3        0  25302.130300',
                              b'ZETA 3        0  12000.000000')], **SMALL))
    assert avdevs(shaken.oute)[0][0] is None      # AVDEV printed as ******
    assert summary(shaken.oute) == summary(steady.oute)
    lettered = {line[:13].strip(): line[13:14]
                for line in shaken.rceout.split(b'\r\n')
                if len(line) > 14 and line[13:14].strip()}
    assert lettered.get(b'810836.077000') == b'X'
    assert lettered.get(b'776509.696000') == b'x'
    assert lettered.get(b'846376.043000') == b'x'
    assert (avdevs(shaken.oute)[0][-1]
            == pytest.approx(avdevs(steady.oute)[0][-1], abs=0.05))
