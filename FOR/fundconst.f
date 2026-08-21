C     ==================================================================
C     fundconst.f -- every fundamental physical constant used anywhere
C     in the RCN / RCN2 / RCG / RCE package, in one place.
C
C     A. Kramida 2026.
C
C     USAGE
C          INCLUDE 'fundconst.f'
C     anywhere in the specification part of a program unit, after any
C     IMPLICIT statement.  The file contains type declarations and
C     PARAMETER statements only -- no COMMON, no DATA, no executable
C     statement -- so it is safe to include in any routine, and two
C     routines that include it cannot disagree about a value.
C
C     Every name in Sections A and B is declared DOUBLE PRECISION
C     explicitly rather than relying on the host routine's IMPLICIT
C     REAL*8 (A-H,O-Z).  Not all of RCN36K.F carries that IMPLICIT card
C     (see the note there about commenting them out for 64-bit
C     machines), and a constant silently demoted to single precision is
C     exactly the class of bug this file exists to remove.
C
C     ORDER MATTERS.  A PARAMETER may only use names already defined
C     above it, so an entry that is derived from others is placed after
C     them: the pure numbers precede the radiative coefficients because
C     GACON is built from PICON.  Keep it that way when adding an entry.
C
C     None of these names occurs anywhere else in the package.  Checked
C     over FOR/*.F, FOR/*.f and FOR/LSJJ/ : zero hits for every one.
C     Re-run that check before adding a name.  A PARAMETER that shadows
C     a local variable which is only ever read is a silent wrong-answer
C     bug, not a compile error.
C
C     ------------------------------------------------------------------
C     THE RULE:  EVERY FUNDAMENTAL CONSTANT USED ANYWHERE IN THE
C     PACKAGE BELONGS IN THIS FILE.  There is no threshold of
C     smallness below which a constant may stay written inline: a
C     literal left in the Fortran is a literal that the next change
C     of these values cannot reach, and reaching all of them is the
C     whole purpose of the file.  CODE/tests/test_literals.py
C     enforces it over FOR/ and keeps no exception list.
C     ------------------------------------------------------------------
C
C     ------------------------------------------------------------------
C     EVERY VALUE IS CODATA 2022, AND EACH ENTRY KEEPS ITS OWN HISTORY.
C     ------------------------------------------------------------------
C     Section B now holds current values throughout, but the number the
C     Fortran used to carry is kept above each one as
C
C          C     LEGACY: PARAMETER (NAME = <the old literal>)
C
C     bit for bit as it was written, together with its relative error.
C     That line is not decoration.  It is what makes the substitution
C     step provable: a throwaway build with every LEGACY line activated
C     and the line below it commented out must reproduce the previous
C     output byte for byte, and only then is the value change read on
C     its own.  RCN, RCN2 and RCG were each done that way; keep the
C     convention if you add an entry, and keep the spelling exact,
C     because a script depends on it.
C
C     CODATA 2022 values are transcribed from the NIST complete listing,
C     https://physics.nist.gov/cuu/Constants/Table/allascii.txt , and
C     derived quantities are computed from them, not quoted from memory.
C
C     ------------------------------------------------------------------
C     WARNING -- DEFAULT-REAL LITERALS
C     ------------------------------------------------------------------
C     RCN36K.F writes its constants with no exponent suffix (137.036,
C     2.66246e-5, 2.921830).  Those are DEFAULT REAL literals: they are
C     rounded to single precision and only then promoted where they meet
C     a double.  Substituting a double-precision PARAMETER changes the
C     result even though the decimal digits are identical -- by about
C     3e-8 relative, which is larger than several of the constant
C     updates in this file.
C
C     Sites affected: all 22 in RCN36K.F, plus one in rcg11kd.f.
C     RCN2K.F and the rest of rcg11kd.f use D0 suffixes throughout and
C     are not affected.
C
C     Each of those sites must be decided explicitly when it is
C     substituted -- either accept the change and account for it, or
C     give the parameter the single-precision-rounded value.  Do not let
C     the decision hide inside the constants diff.
C
C     The 22 RCN36K.F sites are DONE.  All 22 were decided the same way,
C     and the reason is in that program's header under Vers. 36.4.2: the
C     substitution was made first with parameters carrying the single
C     rounding, so that it was provably a no-op; the promotion to double
C     precision was then a step of its own, worth about 5e-8 relative;
C     and the values were changed only after that.  There is no site
C     among the 22 where the single rounding was doing anything but
C     losing digits on the way into double arithmetic.
C
C     The one site outside RCN36K.F is done too:
C          rcg11kd.f, subroutine ENERGY   ERAS=ERAS*109737.0/UENRGY
C     and it cost nothing.  109737.0 is a whole number far inside
C     binary32's 24-bit significand, so it was exact as a single and
C     exact again on promotion: that site had no precision component at
C     all, only a value one, and it now takes RINF like the others.
C     ==================================================================
C
C
C     ==================================================================
C     SECTION A -- CODATA 2022, exact as published.
C
C     Nothing in the package uses these yet.  They are the reference set
C     the Section B values will be rewritten from, and the input to the
C     E1 / M1 / E2 conversion coefficients.  Being unused, they cannot
C     affect any existing result.
C     ==================================================================
C
C     Speed of light in vacuum, m/s and cm/s.
C     Exact by definition of the metre.
      DOUBLE PRECISION CLIGHM, CLIGHT
      PARAMETER (CLIGHM = 299792458.D0)
      PARAMETER (CLIGHT = 2.99792458D10)
C
C     Planck constant, J s.  Exact by definition of the kilogram.
      DOUBLE PRECISION PLANCK
      PARAMETER (PLANCK = 6.62607015D-34)
C
C     Elementary charge, C.  Exact by definition of the ampere.
      DOUBLE PRECISION ECHARG
      PARAMETER (ECHARG = 1.602176634D-19)
C
C     Boltzmann constant, J K-1.  Exact by definition of the kelvin.
      DOUBLE PRECISION KBOLTZ
      PARAMETER (KBOLTZ = 1.380649D-23)
C
C     Avogadro constant, mol-1.  Exact by definition of the mole.
      DOUBLE PRECISION AVOGAD
      PARAMETER (AVOGAD = 6.02214076D23)
C
C     Electron mass, kg.                       u(r) = 3.1e-10
      DOUBLE PRECISION EMASS
      PARAMETER (EMASS = 9.1093837139D-31)
C
C     Vacuum electric permittivity, F m-1.     u(r) = 1.6e-10
      DOUBLE PRECISION EPS0
      PARAMETER (EPS0 = 8.8541878188D-12)
C
C     Vacuum magnetic permeability, N A-2.     u(r) = 1.6e-10
      DOUBLE PRECISION AMU0
      PARAMETER (AMU0 = 1.25663706127D-6)
C
C     Bohr magneton, J T-1.                    u(r) = 3.1e-10
      DOUBLE PRECISION BOHRMG
      PARAMETER (BOHRMG = 9.2740100657D-24)
C
C     Electron g-factor, magnitude.            u(r) = 1.3e-13
C     CODATA 2022 publishes it as g_e = -2.00231930436092(36); the
C     sign is a convention about the direction of the moment relative
C     to the spin and the package uses |g| everywhere, so the
C     magnitude is what is stored.  A sign flipped here would be
C     silent, which is why it is said out loud.  This is the most
C     precisely known number in this file.
      DOUBLE PRECISION GEFACT
      PARAMETER (GEFACT = 2.00231930436092D0)
C
C     Bohr radius, m and cm.                   u(r) = 1.6e-10
      DOUBLE PRECISION ABOHRM, ABOHR
      PARAMETER (ABOHRM = 5.29177210544D-11)
      PARAMETER (ABOHR  = 5.29177210544D-9)
C
C     Hartree energy, J.                       u(r) = 1.1e-12
      DOUBLE PRECISION HARTJ
      PARAMETER (HARTJ = 4.3597447222060D-18)
C
C     Fine-structure constant, and its square. u(r) = 1.6e-10
C     1/ALPHFS = 137.035999177(21).
      DOUBLE PRECISION ALPHFS, ALPHF2
      PARAMETER (ALPHFS = 7.2973525643D-3)
      PARAMETER (ALPHF2 = ALPHFS*ALPHFS)
C
C     Rydberg constant, cm-1.                  u(r) = 1.1e-12
C     NIST gives 10 973 731.568 157(12) m-1.
      DOUBLE PRECISION RINFCD
      PARAMETER (RINFCD = 109737.31568157D0)
C
C     Hartree in cm-1, i.e. 2*Rinf: the atomic unit of energy in the
C     units RCN, RCN2 and RCG print.  NIST hartree-inverse metre
C     relationship 2.194 746 313 6314(24) e7 m-1.
      DOUBLE PRECISION AUCM
      PARAMETER (AUCM = 219474.63136314D0)
C
C
C     ==================================================================
C     SECTION B -- the constants the code actually uses.
C
C     Each entry gives its CODATA 2022 value, the legacy literal it
C     replaced, and the relative error that literal had.  Where several
C     sites rounded the same constant differently, each rounding was
C     given its own name so that the substitution step could stay
C     bit-exact; once the values were updated those names had become
C     equal, and they were folded into one -- RINFL0, ALPINL, PICONA
C     and GACONX are gone, and the sites that used them are listed
C     under RINF, ALPHIN, PICON and GACON.
C     ==================================================================
C
C     ---- Rydberg constant --------------------------------------------
C
C     Rinf in cm-1.  RCN36K.F 3892 3894 3895 3903 4122; RCN2K.F 457
C     1110; rcg11kd.f 11138 11649 11818 13801 14719 14862.
C     CODATA 2022: 109737.31568157      legacy error -5.18e-08
C     The first rcg11kd.f site is the exception: it carried 109737.0,
C     error -2.88e-06, and as a default-real literal at that -- which
C     cost nothing, 109737.0 being a whole number well inside
C     binary32's 24-bit significand.  It was folded into RINF once the
C     values were equal.
C     The autoionization site, rcg11kd.f 11608, is a second exception:
C     it carried 109735.3, error -1.84e-05, and was folded in the same
C     way (Section 10 of the plan).  RINF therefore now stands for
C     THREE distinct legacy roundings -- 109737.31, 109737.0 and
C     109735.3 -- and the single LEGACY line below can restore only
C     the first.  legacyize.py no longer reproduces the original text
C     of every site, and a Stage A build that must be bit-exact has to
C     give the odd sites temporary names of their own, as Sections 6
C     and 10 did.
      DOUBLE PRECISION RINF
C     LEGACY: PARAMETER (RINF = 109737.31D0)
      PARAMETER (RINF = RINFCD)
C
C     Rinf in 1000 cm-1, the unit RCN2 carries parameters in.
C     RCN36K.F 1103; RCN2K.F 1210 1219 1239 1438 1709 1718 1723 1735
C     1750 2205.  The 1718 site carried the coarser 109.737, error
C     -2.88e-06, and was folded in once the values were equal.
C     CODATA 2022: 109.73731568157      legacy error -5.18e-08
      DOUBLE PRECISION RINFKK
C     LEGACY: PARAMETER (RINFKK = 109.73731D0)
      PARAMETER (RINFKK = RINFCD/1.D3)
C
C     ---- Fine-structure constant -------------------------------------
C
C     1/alpha.  RCN36K.F 2129 2258 2738 3857.
C     CODATA 2022: 137.035999177        legacy error +6.01e-09
C     There used to be a site in rcg11kd.f as well, the M1/E1 ratio
C     (alpha/2)**2, which carried the coarser 137.04, error +2.92e-05,
C     and was folded in once the values were equal.  It went with the
C     inline PSQ derivation when SPECTR moved onto GFFROMSCOEF, so
C     ALPHIN now has no consumer in RCG at all.  It is kept because it
C     is a plain CODATA-derived name rather than a unit-scaled variant.
      DOUBLE PRECISION ALPHIN
C     LEGACY: PARAMETER (ALPHIN = 137.036D0)
      PARAMETER (ALPHIN = 1.D0/ALPHFS)
C
C     2/alpha.  RCN36K.F 2259.
C     CODATA 2022: 274.071998354        legacy error +6.01e-09
      DOUBLE PRECISION TWOOAL
C     LEGACY: PARAMETER (TWOOAL = 274.072D0)
      PARAMETER (TWOOAL = 2.D0/ALPHFS)
C
C     alpha**2/2, the spin-orbit coefficient.  RCN36K.F 3885 3886 3887;
C     RCN2K.F 456.
C     CODATA 2022: 2.6625677223848D-05  legacy error -4.05e-05
      DOUBLE PRECISION ALPH2H
C     LEGACY: PARAMETER (ALPH2H = 2.66246D-5)
      PARAMETER (ALPH2H = ALPHF2/2.D0)
C
C     alpha**2/4.  RCN36K.F 3889 3890 (mass-velocity, Darwin).
C     CODATA 2022: 1.3312838611924D-05  legacy error -4.05e-05
      DOUBLE PRECISION ALPH2Q
C     LEGACY: PARAMETER (ALPH2Q = 1.33123D-5)
      PARAMETER (ALPH2Q = ALPHF2/4.D0)
C
C     (alpha**2/2)*Rinf, i.e. ALPH2H*RINF, in cm-1.  RCN36K.F 5118 5143
C     5222 (sm, sn, vk).  This one was updated at some point and is
C     already right to all the digits it carries -- unlike ALPH2H above,
C     which is the same physical quantity from a pre-1969 alpha.  That
C     inconsistency is the defect recorded in the plan as 0a.2.
C     CODATA 2022: 2.9218303467490      legacy error -1.19e-07
      DOUBLE PRECISION ZETCON
C     LEGACY: PARAMETER (ZETCON = 2.921830D0)
      PARAMETER (ZETCON = ALPHF2*RINFCD/2.D0)
C
C     alpha**2*Rinf, twice the above.  RCN36K.F 5001.
C     CODATA 2022: 5.8436606934980      legacy error -1.19e-07
      DOUBLE PRECISION ZETCN2
C     LEGACY: PARAMETER (ZETCN2 = 5.843660D0)
      PARAMETER (ZETCN2 = ALPHF2*RINFCD)
C
C     ---- Energy conversions ------------------------------------------
C
C     eV -> cm-1.  rcg11kd.f 13157 14975.
C     CODATA 2022: 8065.543937349       legacy error -8.05e-06
C     The CODATA value is exact: it follows from the exact e, h and c,
C     so it is written here as that quotient rather than transcribed.
      DOUBLE PRECISION EV2CM
C     LEGACY: PARAMETER (EV2CM = 8065.479D0)
      PARAMETER (EV2CM = ECHARG/(PLANCK*CLIGHT))
C
C     keV -> cm-1, the same number times 1000.  rcg11kd.f 13155 14094
C     14150 14507.
C     CODATA 2022: 8065543.937349       legacy error -8.05e-06
      DOUBLE PRECISION EV2CMK
C     LEGACY: PARAMETER (EV2CMK = 8065479.D0)
      PARAMETER (EV2CMK = EV2CM*1.D3)
C
C     Rydberg energy in eV.  rcg11kd.f 14984.  Derived from RINFCD and
C     EV2CM rather than transcribed, so that the three cannot drift
C     apart; the quotient is 13.605693122991, which is the NIST
C     listing's 13.605 693 122 990 to within the 1.1e-12 uncertainty
C     on Rinf.
C     CODATA 2022: 13.605693122990      legacy error -5.09e-05
      DOUBLE PRECISION RYEV
C     LEGACY: PARAMETER (RYEV = 13.605D0)
      PARAMETER (RYEV = RINFCD/EV2CM)
C
C     Rydberg energy in keV.  rcg11kd.f 14078 14096.
C     CODATA 2022: 0.013605693122990    legacy error +7.86e-06
      DOUBLE PRECISION RYKEV
C     LEGACY: PARAMETER (RYKEV = 0.0136058D0)
      PARAMETER (RYKEV = RYEV/1.D3)
C
C     ---- Pure numbers ------------------------------------------------
C
C     pi, 13 digits.  RCN36K.F 2630 3158; rcg11kd.f 17848, the E2 type
C     factor inside GAFROMSCOEF.
C     exact pi: 3.14159265358979        legacy error +6.6e-14
C     The rcg11kd.f site carried the coarser 3.14159, error -8.45e-07,
C     and was folded in once the values were equal.
      DOUBLE PRECISION PICON
C     LEGACY: PARAMETER (PICON = 3.141592653590D0)
      PARAMETER (PICON = 3.14159265358979324D0)
C
C     2*pi.  rcg11kd.f 11649, the autoionization rate.
C     exact 2*pi: 6.28318530717959     legacy error -8.45e-07
C     Named TWOPIC and not TWOPI: RCN36K.F subroutine phshift has a
C     local variable twopi and that routine includes this file, so
C     TWOPI here would be an assignment to a PARAMETER.
      DOUBLE PRECISION TWOPIC
C     LEGACY: PARAMETER (TWOPIC = 6.28318D0)
      PARAMETER (TWOPIC = 2.D0*PICON)
C
C     4/3, written out to 10 digits.  RCN2K.F 649.
C     exact 4/3: 1.33333333333333       legacy error -2.50e-10
      DOUBLE PRECISION FOURO3
C     LEGACY: PARAMETER (FOURO3 = 1.333333333D0)
      PARAMETER (FOURO3 = 4.D0/3.D0)
C
C     ---- Radiative-transition coefficients ---------------------------
C
C     gf(E1) = GFE1 * sigma(cm-1) * S(a.u.), with GFE1 = 2/(3*2*Rinf).
C     rcg11kd.f 17841, inside GAFROMSCOEF.  SPECTR's PSQ used to carry
C     this inline; it now comes from GFFROMSCOEF.
C     CODATA 2022: 3.0375568352755D-06  legacy error -1.87e-05
      DOUBLE PRECISION GFE1
C     LEGACY: PARAMETER (GFE1 = 3.0375D-6)
      PARAMETER (GFE1 = 2.D0/(3.D0*AUCM))
C
C     gA = GACON * gf * sigma(cm-1)**2, with
C     GACON = 8*pi**2*alpha**2*a0(cm)*c(cm/s).
C     rcg11kd.f 13488 13593 13735 13838 13940 14014 14303 14758, and
C     17841 17879 in the conversion-coefficient functions.  The
C     first of the eight is the TEXCIT line-selection threshold in
C     SPECTR, which used to carry 0.6698 -- one digit dropped, +4.16e-03
C     against CODATA, sixty times the error of the 0.66698 beside it.
C     That is the defect recorded in the plan as 0a.1, and folding it
C     into GACON is the fix.  It is not a printed quantity but it
C     decides which lines survive the cut, so a run with TEXCIT>0 can
C     now select a slightly different set of lines.  None of the three
C     test cases uses TEXCIT.
C     CODATA 2022: 0.667025176604160    legacy error -6.77e-05
      DOUBLE PRECISION GACON
C     LEGACY: PARAMETER (GACON = 0.66698D0)
      PARAMETER (GACON = 8.D0*PICON**2*ALPHF2*ABOHR*CLIGHT)
C
C     ---- Autoionization and photoionization --------------------------
C
C     hbar, J s.  No site uses it directly; it is the input to AURATE.
C     It sits here and not in Section A, where a CODATA quantity would
C     otherwise belong, because it is built from PICON and Section A
C     stands above PICON.  The alternative was a second spelling of pi
C     in this file, which is the one thing the file exists to prevent.
C     CODATA 2022: 1.05457181764616D-34, exact from the exact h.
      DOUBLE PRECISION HBARJ
      PARAMETER (HBARJ = PLANCK/(2.D0*PICON))
C
C     Ry/hbar in units of 1e13 s-1: the factor that turns 2*pi*|V|**2,
C     with V in rydbergs, into an autoionization rate in s-1.  Fermi's
C     golden rule.  rcg11kd.f 11649.
C     CODATA 2022: 2067.0686667586      legacy error -5.17e-04
C     2066 is a truncation of 2067.07, not a rounding, and it is the
C     largest constant error found anywhere in this package.  Unlike
C     the rest it is exercised: both WORK/AI decks print AA tables
C     built from this line, and WORK/AI/DIEL's GAATOT, BRION and DR
C     rates all follow from them.  Autoionization rates printed by RCG
C     before 2026-08-20 are low by 5.2e-04 relative.  The printed AA
C     rises by less than that, 4.81e-04, because RINF appears squared
C     in the denominator of the same expression and its own correction
C     of -1.84e-05 pulls the other way.
      DOUBLE PRECISION AURATE
C     LEGACY: PARAMETER (AURATE = 2066.D0)
      PARAMETER (AURATE = HARTJ/(2.D0*HBARJ*1.D13))
C
C     4*pi**2*alpha*a0**2, cm**2: the photoionization cross-section
C     coefficient CSIGMA.  rcg11kd.f 13201.
C     CODATA 2022: 8.0672837256768D-18  legacy error -8.23e-05
C     It gates the SIGMA path, which none of the decks in WORK/
C     exercises; the value was verified algebraically instead.
      DOUBLE PRECISION PICSCN
C     LEGACY: PARAMETER (PICSCN = 8.06662D-18)
      PARAMETER (PICSCN = 4.D0*PICON**2*ALPHFS*ABOHR**2)
C
C     ---- Electron g-factor, as the code spells it ---------------------
C
C     Three names for one constant, in the three arrangements the
C     Fortran needs, exactly as RINF / RINFKK and ALPH2H / ALPH2Q are
C     three spellings of two.  All are derived from GEFACT above.
C
C     (|g|-1)/2, the coefficient of the S(S+1)-L(L+1) term of the
C     Lande g-factor.  rcg11kd.f 1367 (SPRIN), 9019 9142 (CALCFC),
C     15287 (RCEINP).
C     CODATA 2022: 0.50115965218046      legacy error +6.94e-07
      DOUBLE PRECISION GEM1H
C     LEGACY: PARAMETER (GEM1H = 0.50116D0)
      PARAMETER (GEM1H = (GEFACT-1.D0)/2.D0)
C
C     1 + (|g|-1)/2, the additive term of the same expression.
C     rcg11kd.f 1368 (SPRIN), 9020 9143 (CALCFC), 15287 (RCEINP).
C     CODATA 2022: 1.50115965218046      legacy error +2.32e-07
      DOUBLE PRECISION GEP1H
C     LEGACY: PARAMETER (GEP1H = 1.50116D0)
      PARAMETER (GEP1H = 1.D0+GEM1H)
C
C     |g|-1, the spin part of the magnetic-moment matrix element.
C     rcg11kd.f 10495 (MUPOLE).
C     CODATA 2022: 1.00231930436092      legacy error +6.94e-07
      DOUBLE PRECISION GEM1
C     LEGACY: PARAMETER (GEM1 = 1.00232D0)
      PARAMETER (GEM1 = GEFACT-1.D0)
C
C     ---- Collision rate coefficient -----------------------------------
C
C     2*sqrt(pi)*alpha*c*a0**2, cm3 s-1: the Maxwellian rate
C     coefficient per unit collision strength, with the temperature in
C     rydbergs.  rcg11kd.f 14986, C2=CMAXWL/(SWT0*DSQRT(TRYD)).
C     CODATA 2022: 2.1716669877207D-08  legacy error +1.52e-05
      DOUBLE PRECISION CMAXWL
C     LEGACY: PARAMETER (CMAXWL = 2.1717D-8)
      PARAMETER (CMAXWL = 2.D0*DSQRT(PICON)*ALPHFS*CLIGHT*ABOHR**2)
C
C
C
C
C     ==================================================================
C     end of fundconst.f
C     ==================================================================
