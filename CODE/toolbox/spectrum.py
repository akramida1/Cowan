"""Reading the transition sections of OUTG11.

After the eigenvectors RCG prints a spectrum: for every initial level a
`* * *` header line, then one line per final level with the wavelength,
the oscillator strength and the transition probability.  A run may hold
three of these, one per multipole -- MAG DIP, ELEC DIP and ELEC QUD
(`SPID`, `rcg11kd.f:12812`) -- and the first two are printed for the
same pair of parity sets in both directions.

Two layouts are read, told apart by the column header each section
prints above its lines.

The **Lund** layout is what `ILUND` selects (`rcg11kd.f:13424-13440`)
and what the two level sorts use: a `* * *` line for the initial level,
then one line per final level.  The **wide** layout (FORMAT 58,
`rcg11kd.f:13473`) carries both levels on every line and is what the
wavelength sort always prints, whatever `ILUND` says -- so an ISPECC of
4 or more (column 72 of the ING11 main card) produces it, and an ISPECC
of exactly 4 produces nothing else.

Both carry the level's number within its J block -- the wide one since
A. Kramida added NLEVJ and NLEVJP to FORMAT 58 in 2018 -- and that
number is what ties a line to an eigenvector without guessing: the
alternative is to match on energy, configuration number and term label,
which is what `conv_out.pl:FindLev` did and what it needed its several
fall-backs for.  The remaining layouts do not carry it and are refused
by name rather than half-read: FORMAT 57, the photoionization
cross-section form (`CSIGMA`), FORMAT 51/53, the Boltzmann-weighted one
(`TEXCIT`), and FORMAT 88, the branching-ratio variant of the
wavelength sort (`IBRNCH`, `rcg11kd.f:13071`).

The wide layout has no `* * *` line to carry the initial level's
configuration name, so both names are taken from the section heading's
own table of configurations (FORMAT 55), which numbers them for each
parity set exactly as the `NCSER`/`NCSERP` columns index them.  It also
truncates the 8-character basis-state label to seven (`A7`); the eighth
character is blank on every label this package produces, but a label
read from a wide section is the shorter field all the same.

Energies are in the unit named in the section heading -- 1000 cm-1 in
every deck of this package -- and are left in it.  A wavelength in
angstroms is 1e5 / (delta E in kK), which is where the `1e5` in the
callers comes from.

Beside the transitions each `* * *` block closes with a SUMFI record
carrying RCG's own summed decay rate for the block's level (FORMAT
1059).  `totals()` reads those; see its docstring for why they are
worth reading rather than summing the printed lines.

This module streams.  A Pr III OUTG11 holds three million lines of
spectrum, so `lines()` and `totals()` are generators and nothing
accumulates here.
"""

from dataclasses import dataclass, replace

from . import files

#: `SPID`, `rcg11kd.f:12812`.  The two that connect a parity set to
#: itself are the ones whose operator is even.
MULTIPOLES = ('MAG DIP', 'ELEC DIP', 'ELEC QUD')
SAME_PARITY = ('MAG DIP', 'ELEC QUD')


class SpectrumError(Exception):
    """A transitions section that does not read as RCG writes it."""


@dataclass(frozen=True)
class LevelRef:
    """One end of a transition, as the spectrum prints it."""

    energy: float
    j: float
    config: int             # the configuration's number within its set
    term: str               # the 8-character label, as in an eigenvector row
    level: int              # its number within its (parity set, J) block
    name: str = ''          # the 18-character 3A6, on the initial line only

    @property
    def label(self):
        """The 12-character configuration label, where one was printed."""
        return self.name[6:18].strip()

    @property
    def spectrum(self):
        return self.name[:6].strip()


@dataclass(frozen=True)
class Line:
    """One printed transition.

    `initial` is the level of the `* * *` header and `final` the one on
    the line itself; which of the two is the upper level is the caller's
    business, and RCG prints both directions.
    """

    multipole: str
    initial: LevelRef
    final: LevelRef
    delta_e: float
    wavelength: float       # None where the F11.3 field overflowed
    log_gf: float
    ga: float
    cf: float

    @property
    def same_parity(self):
        return self.multipole in SAME_PARITY


@dataclass(frozen=True)
class Total:
    """RCG's own summed decay rate for one level, as it prints it.

    The record that closes each `* * *` block of the Lund layout:
    `SUMGA` is the total gA of that block's level -- the sum over its
    decays to lower levels -- and `lifetime` is (2J+1)/SUMGA, the
    radiative lifetime RCG quotes beside it.

    `starred` says the record held an asterisk somewhere: a fixed-point
    field on it overflowed.  `sumga` and `lifetime` are None when the
    overflow was in that field itself, and the caller has to decide
    what to do without it; `starred` with both numbers present means
    only a neighbouring SUMFI field overflowed and the two are good.
    """

    multipole: str
    section: int            # which printed sort of that multipole
    level: LevelRef         # the `* * *` level the record belongs to
    sumga: float
    lifetime: float
    starred: bool


# --------------------------------------------------------------------------
# the fixed columns
# --------------------------------------------------------------------------

#: FORMAT 55, `rcg11kd.f:13371`:
#:   1H1,9X,A8,9H SPECTRUM,5X,21H(ENERGIES IN UNITS OF,F8.1,7H CM-1 =,...
_HEADING = ' SPECTRUM     (ENERGIES IN UNITS OF'


def _heading(line):
    """(multipole, energy unit) of a section heading, or None."""
    if not line.startswith('1') or line[18:53] != _HEADING:
        return None
    multipole = line[10:18].strip()
    if multipole not in MULTIPOLES:
        raise SpectrumError('unknown spectrum type %r in a section heading'
                            % multipole)
    return multipole, float(line[53:61])


#: The tail of FORMAT 55, one record per configuration number:
#:   (I5,6X,3A6,12H   ---      ,3A6)
#: The first 3A6 is the configuration of the first parity set, the
#: second that of the second set -- `ELEM(L,JN,L1)` and `ELEM(L,JN,L2)`
#: -- and `NCSER`/`NCSERP` on a transition line index into them.
_DASHES = '   ---      '


def _configurations(line):
    """(number, first-set name, second-set name) of a heading row."""
    if len(line) < 41 or line[29:41] != _DASHES:
        return None
    number = line[:5].strip()
    if not number.isdigit():
        return None
    return int(number), line[11:29], line[41:59]


#: FORMAT 1056, `rcg11kd.f:13424`, and FORMAT 56, `rcg11kd.f:13407`:
#: the column header each layout prints above its lines, and the only
#: thing that tells them apart without guessing at field widths.
_LUND_COLUMNS = ('            E      J     LEVEL        DELTA E '
                 ' LAMBDA(A)  LOG GF GA(SEC-1) CF/BR')
_WIDE_COLUMNS = ('0           E      J   CONF              EP      JP'
                 '  CONFP           DELTA E   LAMBDA(A) S/PMAX**2   GF'
                 '    LOG GF GA(SEC-1) CF,BRNCH')


#: FORMAT 1057, `rcg11kd.f:13428`:
#:   5X,7H* * *  ,F10.3,F5.1,I4,2H (,A8,I6,1X,3A6,6H * * *
_STARS = '     * * *  '


def _initial(line):
    return LevelRef(energy=float(line[12:22]), j=float(line[22:27]),
                    config=int(line[27:31]), term=line[33:41].rstrip(),
                    level=int(line[41:47]), name=line[48:66])


#: FORMAT 1058, `rcg11kd.f:13437`:
#:   1H ,I4,F11.3,F5.1,I3,1X,1H(,A8,F11.3,F11.3,F8.3,1PE10.3,0PF6.2,I6
def _final(line):
    """The final level and the line's radiative data, or None.

    The wavelength field overflows its F11.3 for a transition between
    two nearly degenerate levels and prints as asterisks; the energies
    are still exact, so the caller recomputes it rather than losing the
    line.
    """
    wavelength = line[45:56]
    return (LevelRef(energy=float(line[5:16]), j=float(line[16:21]),
                     config=int(line[21:24]), term=line[26:34].rstrip(),
                     level=int(line[80:86])),
            float(line[34:45]),
            None if '*' in wavelength else float(wavelength),
            float(line[56:64]), float(line[64:74]), float(line[74:80]))


def _is_final(line):
    """A final-level line, told from the headings and sums around it.

    The leading I4 is the line's serial number over the whole spectrum
    and overflows to asterisks past 9999 -- Pr III reaches seven
    figures.  Nothing is read from it, so it is only tested for shape.
    """
    if len(line) < 86 or line[0] != ' ' or line[24:26] != ' (':
        return False
    counter = line[1:5].strip()
    if not (counter.isdigit() or set(counter) == {'*'}):
        return False
    return line[80:86].strip().isdigit()


#: FORMAT 58, `rcg11kd.f:13473`:
#:   1H ,I4,F11.4,F5.1,I3,2H (,A7,F13.4,F5.1,I3,2H (,A7,F13.4,F12.4,
#:   F8.3,F9.4,F8.3,1PE11.4,0PF10.6,I6,I6
#: S/PMAX**2 and GF sit between the wavelength and LOG GF here and are
#: not read; the Lund layout does not print them, and nothing in this
#: package wants them.
def _wide(line):
    """Both levels and the radiative data of a wide-layout line."""
    wavelength = line[76:88]
    return (LevelRef(energy=float(line[5:16]), j=float(line[16:21]),
                     config=int(line[21:24]), term=line[26:33].rstrip(),
                     level=int(line[134:140])),
            LevelRef(energy=float(line[33:46]), j=float(line[46:51]),
                     config=int(line[51:54]), term=line[56:63].rstrip(),
                     level=int(line[140:146])),
            float(line[63:76]),
            None if '*' in wavelength else float(wavelength),
            float(line[105:113]), float(line[113:124]), float(line[124:134]))


def _is_wide(line):
    """A wide-layout transition line, told from the headings around it.

    FORMAT 88, the branching-ratio variant, prints the first label as
    A8 and so puts its second `2H (` one column further along; it is
    rejected here rather than read with FORMAT 58's offsets, and
    `lines()` names it if a whole section turns out to be in it.
    """
    if len(line) < 146 or line[0] != ' ':
        return False
    if line[24:26] != ' (' or line[54:56] != ' (':
        return False
    counter = line[1:5].strip()
    if not (counter.isdigit() or set(counter) == {'*'}):
        return False
    return (line[134:140].strip().isdigit()
            and line[140:146].strip().isdigit())


#: FORMAT 1059, `rcg11kd.f:13986`, which closes every `* * *` block of
#: the Lund layout:
#:   /1X,6HSUMFI=,F7.4,F8.4,2H =,F8.4,2X,9HRAD.LIFE=,1PE10.3,4H SEC,
#:   2X,6HSUMGA=,E10.3/
#: The wide layout's FORMAT 59 and 62 carry the same two numbers but
#: name the level by its energy alone, with neither the J-block level
#: number nor the configuration; there is nothing to tie those to an
#: eigenvector, so only this one is read.
_SUMFI = ' SUMFI='
_RADLIFE = 'RAD.LIFE='
_SUMGA = 'SUMGA='


def _total(line):
    """(SUMGA, lifetime, starred) of a FORMAT 1059 record, or None.

    A field of the record overflows into asterisks -- the author has
    seen it -- and the two read here are the last two, so an overflow
    in the SUMFI fields ahead of them does not move their columns.
    Either is returned as None when it is the one that overflowed.
    """
    if not line.startswith(_SUMFI):
        return None
    if line[34:43] != _RADLIFE or line[59:65] != _SUMGA:
        return None
    return (_number(line[65:75]), _number(line[43:53]), '*' in line)


def _number(field):
    return None if '*' in field else float(field)


# --------------------------------------------------------------------------
# reading
# --------------------------------------------------------------------------

def lines(path=files.OUTG11, folder=None, multipoles=None):
    """Every printed transition of an OUTG11, in file order.

    `multipoles` filters the sections; a spectrum whose lines are not
    wanted is skipped without being parsed, which is most of the file.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    wanted = set(multipoles) if multipoles else None
    multipole = None
    layout = None
    initial = None
    first_set = {}
    second_set = {}
    scanned = read_any = False
    refused = set()

    with open(path, 'r', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n').rstrip('\r')
            found = _heading(line)
            if found is not None:
                multipole, _ = found
                layout = initial = None
                first_set, second_set = {}, {}
                continue
            if multipole is None:
                continue
            if wanted is not None and multipole not in wanted:
                continue

            names = _configurations(line)
            if names is not None:
                number, one, two = names
                first_set[number] = one
                second_set[number] = two
                continue
            if line == _LUND_COLUMNS:
                layout, initial = 'lund', None
                continue
            if line == _WIDE_COLUMNS:
                layout, initial = 'wide', None
                continue

            scanned = True
            if layout is None:
                # A column header of a layout this module does not read
                # -- kept so the error at the end can name it.
                if 'DELTA E' in line and 'LAMBDA(A)' in line:
                    refused.add(line.strip())
                continue

            if layout == 'lund':
                if line.startswith(_STARS):
                    read_any = True
                    initial = _initial(line)
                    continue
                if initial is None or not _is_final(line):
                    continue
                final, delta_e, wavelength, log_gf, ga, cf = _final(line)
            else:
                if not _is_wide(line):
                    continue
                one, two, delta_e, wavelength, log_gf, ga, cf = _wide(line)
                # The wide layout prints no name beside either level;
                # both come from the heading's table of configurations.
                first = replace(one, name=first_set.get(one.config, ''))
                final = replace(two, name=second_set.get(two.config, ''))

            read_any = True
            yield Line(multipole=multipole,
                       initial=initial if layout == 'lund' else first,
                       final=final,
                       delta_e=delta_e, wavelength=wavelength,
                       log_gf=log_gf, ga=ga, cf=cf)

    if multipole is None:
        raise SpectrumError(
            '%s carries no transition section.  RCG prints one only when it '
            'is asked to: see the ING11 main card.' % path)
    if scanned and not read_any:
        raise SpectrumError(
            '%s prints its transitions in a layout this package does not '
            'read -- neither the Lund one (FORMAT 1056) nor the wide one '
            'the wavelength sort uses (FORMAT 56/58).  The column header '
            'it does print is %s.  See rcg11kd.f:13407-13473 for the '
            'others: the photoionization cross-section form, the '
            'Boltzmann-weighted one, and the branching-ratio variant, none '
            'of which carry both level numbers (NLEVJ, NLEVJP).'
            % (path, ' and '.join(repr(one) for one in sorted(refused))
               or 'not one of them'))


def totals(path=files.OUTG11, folder=None, multipoles=None):
    """RCG's own total decay rate for every level of a Lund-layout run.

    One `Total` per FORMAT 1059 record, in file order, carrying the
    `* * *` level the record closes.  This is RCG's `SUMGA`, and it is
    worth reading rather than summing the printed lines because of
    where it accumulates: at label 580 (`rcg11kd.f:13956`), which the
    `GO TO 578` of the DMIN test jumps *past* the print but *into*.  So
    SUMGA counts every line of the block, printed or not, and a
    branching fraction taken from it is not censored by DMIN.

    It is censored by the FLBDN/FLBDX wavelength window, which drops a
    line at label 498 before it ever enters the arrays; a total from a
    run with a window in force is a total over what is in the window.

    A run may print the same level twice, once in each of the two level
    sorts, when the multipole connects a parity set to itself; the two
    records carry the same sum, and it is the caller that must not add
    them together.  `section` counts the sections of one multipole so
    that they can be told apart.
    """
    if folder is not None:
        path = files.resolve(path, folder)
    wanted = set(multipoles) if multipoles else None
    multipole = None
    sections = {}
    section = 0
    layout = None
    initial = None

    with open(path, 'r', errors='replace') as fh:
        for raw in fh:
            line = raw.rstrip('\n').rstrip('\r')
            found = _heading(line)
            if found is not None:
                multipole, _ = found
                sections[multipole] = section = sections.get(multipole, 0) + 1
                layout = initial = None
                continue
            if multipole is None:
                continue
            if wanted is not None and multipole not in wanted:
                continue
            if line == _LUND_COLUMNS:
                layout, initial = 'lund', None
                continue
            if line == _WIDE_COLUMNS:
                layout, initial = 'wide', None
                continue
            if layout != 'lund':
                continue
            if line.startswith(_STARS):
                initial = _initial(line)
                continue
            if initial is None:
                continue
            record = _total(line)
            if record is None:
                continue
            sumga, lifetime, starred = record
            yield Total(multipole=multipole, section=section, level=initial,
                        sumga=sumga, lifetime=lifetime, starred=starred)
            initial = None


def name(path=files.OUTG11, folder=None):
    """The ion's name, as the spectrum's own header lines carry it."""
    for line in lines(path, folder):
        return line.initial.spectrum
    return ''
