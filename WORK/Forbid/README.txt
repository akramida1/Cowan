WORK/Forbid -- a forbidden-line case for calc_TP
===============================================

Sn VIII 4d7, one configuration and one parity, with both forbidden
multipoles switched on.  It is the smallest deck that exercises
calc_TP end to end: an RCG run with MAG DIP and ELEC QUD spectra,
a //R2// radial-integral card in ING11 for the program to vary, and an
RCE parameter table for it to read standard deviations from.

Inputs only, as in WORK/ and WORK/AI.  The commands below are typed
the same way on every platform.  Build it with

    rcn
    rcn2
    rcg
    rce

in an empty directory holding these two files, and then

    calc_TP --write-par
    calc_TP

calc_TP.par is checked in here as a worked example -- M1 and E2, twenty
trials -- so the --write-par step can be skipped.  The program covers
E1 as well now, and this deck has one parity set and therefore no E1
transition at all: it says so and drops E1 from the run rather than
treating it as an error.

CODE/tests/test_calc_tp.py does exactly that, into a temporary
directory, and is the end-to-end test of the program.

WHAT MAKES THIS DECK DIFFERENT FROM THE OTHERS IN WORK/

Columns 49 and 50 of the G5INP card in IN2 are IMAG and IQUAD, RCG's
magnetic-dipole and electric-quadrupole switches (`FOR/rcg11kd.f:608-609`).
Both are 1 here, meaning "for parity set 1".  IQUAD also makes RCN2
compute the E2 radial integrals and punch them as //R2// cards
(`FOR/RCN2K.F:668`), which is where the E2 line strengths come from and
what calc_TP displaces.

RCN2 writes a same-parity //Rk// block **inside** the deck, straight
after the parameter card of the parity set it belongs to, not at the
end with the cross-parity block.  On a two-parity IQUAD deck one of them
therefore sits between parity set 1's parameters and parity set 2's.
`toolbox/ing11.py` keeps such a card on the card it follows so that the
deck still round-trips byte for byte.

DMIN MUST BE ZERO FOR A FORBIDDEN-LINE CALCULATION

Columns 61-65 of the same card are DMIN, and in this version of Cowan's
codes a DMIN between 1 and 50 means "omit every transition with
gA < 10**DMIN" (`FOR/rcg11kd.f:12859-12871`).  The other decks in WORK/
carry 5.00, 9.00 or 10.0, which is a sensible cut for electric-dipole
work: it throws away lines no experiment would see and keeps OUTG11 to
a manageable size.

It is fatal here.  Every M1 and E2 line of this spectrum is weaker than
gA = 1e5, so a deck with DMIN = 5.00 prints the section headings and
not one line under them -- RCG calls SPECTR for the M1 spectrum
(`:11622`) and nothing comes out.  This deck once carried 5.00, and the
missing M1 spectrum was mistaken for a fault in RCG.

**Use 0.00 for any forbidden-transition calculation.**  It is the
default (`:836`) and it is what this deck now has.  With it the run
prints 89 M1 lines and 131 E2 lines over 132 distinct level pairs, 88
of which carry both.  That is what makes the M1 + E2 merge -- the whole
point of calc_TP -- testable on a real run rather than on rows
built up by hand.

ISPECC AND THE TWO SPECTRUM LAYOUTS

Column 72 is ISPECC, which chooses RCG's line sorts: odd values add the
first-parity sort, bit 2 the second-parity sort, and 4 or more the
wavelength sort (`DOCS/RCG_DOC.txt:846`).  This deck carries 2, as every
deck in WORK/ does.

The two level sorts print the Lund layout, a `* * *` line for the
initial level and one line per final level.  The wavelength sort always
prints the wide layout instead, both levels on every line, whatever
ILUND says (`FORMAT 58`, `FOR/rcg11kd.f:13473`) -- so ISPECC = 4 gives
that layout and nothing else.  `toolbox/spectrum.py` reads both; the
wide one carries the level numbers NLEVJ and NLEVJP that tie a line to
an eigenvector because A. Kramida added them to FORMAT 58 in 2018.

Running this deck with a 4 in column 72 and comparing the two tables is
`test_calc_tp.test_a_wavelength_sorted_run_gives_the_same_table`.  The
same level pairs come out of both, with the transition probabilities
agreeing to the extra digit the wider fields carry: the wavelength sort
prints F11.4 energies and an E11.4 gA where the Lund layout prints
F11.3 and E10.3.

The wide layout is also the only one that prints a line strength, in
its S/PMAX**2 column, and
`test_calc_tp.test_rcgs_own_line_strength_is_the_one_radiate_derives`
uses this deck at ISPECC = 4 to check `toolbox/radiate` against it.
PMAX itself is printed for the electric-dipole spectrum only
(`FOR/rcg11kd.f:13341` sets IPTG = 0 for II /= 1), so on this deck the
check is that S(gA)/S2 is one constant over each section -- which is
what a single PMAX per section means, and which a wrong power of the
wavelength or a missing 2J+1 breaks by orders of magnitude.  The
absolute form of the check is made on the E1 spectrum of WORK/.

THE FIT IS NOT A FIT

Running RCE on this deck as it stands gives cycle 0: there are no
observed energies in it, so every level is marked computed and every
parameter comes back with a zero standard deviation.  That is enough to
exercise the whole machinery -- the E2 integrals still vary, RCG is
still re-run, the eigenvectors are still re-identified -- but the
parameter variation itself does nothing, because there is nothing to
vary within.  Adding observed energies would make it a real fit; they
would have to be real measurements, and there are none here.
