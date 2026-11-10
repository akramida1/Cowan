LS -> jj transformation for shells of equivalent electrons
==========================================================
Prototype of "route 2": build both coupled bases explicitly in the m-scheme
and take overlaps.  Standalone Python; nothing here touches RCG.

Why this exists
---------------
RCG's KCPL=2 is not jj coupling for a shell of equivalent electrons.  The
matrix built in CALCFC (rcg11kd.f, lines 7687-7710) vanishes unless every
subshell keeps the same aLS term index, and what remains is a product of 9j
symbols running down the coupling tree:

      IF (NALS(LS,I).NE.NALSJ(JJ,I)) GO TO 530
      ...
      B=S9J(SCRL(LS,I-1),FL(ML,I),SCRL(LS,I), ...)

That is a pure recoupling of the tree; it never changes the internal state of
a subshell.  Line 7697, IF (NOSUBC.LE.1) GO TO 530, makes the point: for a
single open shell the matrix is the identity and the "JJ" labels printed are
the LS labels.  So RCG reports d^2 as 1S 3P 1D 3F 1G with J attached, where
true jj coupling would report (3/2)^2, (3/2)(5/2), (5/2)^2 with seniorities.

The missing coefficient is

    < l^N a L S J | ( j-^N1 v1 J1 , j+^N2 v2 J2 ) J >

which is not a 9j: antisymmetrization mixes the two schemes.  This code
computes it by explicit overlap, so no coefficients of fractional parentage
and no phase conventions have to be reconciled between the schemes.

Method
------
Both bases are built by the same procedure, which is why their relative
phases are meaningful:

  1. Enumerate Slater determinants of the one-electron states as bit masks.
  2. Find the highest-weight space (annihilated by every raising operator) at
     the target projections.
  3. Classify it: seniority from the pair operator A+A, then for f shells the
     Casimir operators of SO(7) and G2, giving Racah's (W) and (U).
  4. Generate every M component by the lowering operators, with
     Condon-Shortley normalization.
  5. Couple to good J with Clebsch-Gordan coefficients.
  6. Expand the jj determinants over the ml,ms determinants (one spin-orbit
     Clebsch-Gordan per orbital) and take the overlaps.

The only angular-momentum algebra involved is the Clebsch-Gordan coefficient.

Two ways of fixing the LS basis
-------------------------------
Step 3 above is the only place where a choice has to be made, and there are
two builders for it.

  lsjj.LSShell            diagonalize the pair operator, and for f also the
                          Casimirs of SO(7) and G2.  Self-contained: needs no
                          data files.  Complete for p and d; for f it leaves
                          8 to 15 states per shell sharing every label.

  cowan_cfp.CowanLSShell  recurse on Cowan's CFP decks, so the alpha labels
                          and the signs are the ones RCG works with.
                          Complete for every shell, because the deck is what
                          defines alpha in the first place.

Add `cowan` to any command line, or call cowan_cfp.install() from Python, to
use the second.  See the "Cowan CFP mode" section below.

Files
-----
The four production modules live in CODE/lsjj/, not here.  They moved there
in 2026 so that make_cfp can find them from beside itself in an installed
tree: CODE/unix/install.sh installs CODE/ and nothing else, and on Windows
the author's habit is to copy CODE\ to C:\cowan\code without the rest of
the repository.  This directory keeps the record of the derivation and its
validation.  Every script here reaches across with a two-line sys.path
insert rather than carrying a second copy.

  CODE/lsjj/lsjj.py       the module and a small command line
  CODE/lsjj/cowan_cfp.py  the LS side built from Cowan's CFP decks
  CODE/lsjj/recur.py      the transformation by recursion on N, from cfp
                          alone -- the algorithm meant for RCG (see
                          "Machinery" below)
  CODE/lsjj/mkdeck.py     writes LSJJ.DAT, the Fortran-side deck

Here:

JJCFP             jj cfp deck for j = 1/2, 3/2, 5/2, 7/2, written by recur.py
validate.py       the validation sweep over p^N, d^N, f^N
compare_racah.py  checks the two LS builders against each other
probe_f3.py       prints the raw Casimir eigenvalues for f^3
recur_all.log     the recursion checked against the m-scheme, every shell
jjcfp.log         the jj cfp tables, printed
validate_all.log       full sweep, diagonalized LS side
validate_all_cowan.log full sweep, CFP-built LS side
compare_all.log        output of compare_racah.py
p2.dat d2.dat          dump format, Racah labels
d2_cowan.dat f3_cowan.dat  dump format, Cowan labels

A jj basis state is labelled  (j- ^n1 v1 J1, j+ ^n2 v2 J2) J .  The seniority
belongs there: 7/2^4 with J = 2, for instance, occurs with both v = 2 and
v = 4, so (j, n, J) alone would not name a state.

Command line
------------
    python lsjj.py demo            p^2, d^2, p^3 with matrices
    python lsjj.py show  l N       one case, all J blocks printed
    python lsjj.py dump  l N file  machine-readable, for reading from Fortran
    python lsjj.py terms l N       LS terms with seniority
    python validate.py all         the whole sweep (p, d and f)
    python compare_racah.py all    the two LS builders against each other
    python cmpfort.py OUTG11       the Fortran recursion against recur.py

and, for the modules that moved, from this directory:

    python ../../CODE/lsjj/recur.py check all      against the m-scheme
    python ../../CODE/lsjj/recur.py check 3 7      one shell
    python ../../CODE/lsjj/recur.py jjcfp          the jj cfp tables
    python ../../CODE/lsjj/recur.py deck JJCFP     rewrite the jj cfp deck
    python ../../CODE/lsjj/mkdeck.py -t deck.txt   the deck as text

LSJJ.DAT itself is not rebuilt by hand any more: make_cfp builds it with
the CFP decks, on the machine the package is installed on, which is the
supported route and the only one an installed user has.

Add `fromdeck` to a recur.py command line to read JJCFP instead of
rebuilding the jj cfp from the m-scheme.

Insert the word `cowan` anywhere in an lsjj.py or validate.py command line
to build the LS side from the CFP decks:

    python lsjj.py cowan show 3 3
    python lsjj.py cowan dump 3 7 f7.dat
    python validate.py all cowan
    python cowan_cfp.py terms 3 3     Cowan's alpha with v recovered from A+A

What is verified
----------------
For every l = 1, 2, 3 and every N from 1 to 4l+2:

  * the LS and jj bases have equal dimension in every J block
  * each basis is orthonormal to ~1e-15
  * the transformation is orthogonal in every J block to ~1e-15
  * the matrix is independent of M (recomputed at 2M = 0 and 2M = 2)

and the LS term lists with seniorities reproduce the published tables for
p^2, p^3, d^2, d^3, d^4, f^2 and f^3.

Spot checks against known analytic values:

    p^2  J=0    sqrt(2/3), 1/sqrt(3)          = 0.81650, 0.57735
    d^2  J=4    2/sqrt(5), 1/sqrt(5)          = 0.89443, 0.44721

The Casimir operators as built here come out as exactly 4*G(G2) and 5*G(R7)
in Racah's normalization.  The f^3 seniority-3 doublets give SO(7) = 9.0
throughout, i.e. W = (210), and G2 = 4.0, 4.6667, 7.0, i.e. U = (11), (20),
(21) -- reproducing Racah's published classification of f^3, in which those
doublets are separated by U alone.

Timing: the whole sweep runs in about a minute; f^7 (m-scheme dimension 3432,
327 states at 2M=1) takes a few seconds.

Where Racah's labels run out
----------------------------
Seniority alone is complete for p^N and d^N.  For f^N it is not, and adding
Racah's (W) and (U) removes most but not all of the remainder:

    shell   states sharing every label   the groups
    f^5      8                           2F, 2H (221)(31); 2I, 2K (221)(31)
    f^6     14                           1G, 1I, 1L (222|311)(40);
                                         3F, 3H, 3I, 3K (221)(31)
    f^7     14                           2F, 2H, 2I, 2K (221)(31);
                                         2G, 2I, 2L (222|311)(40)
    f^8     14                           same as f^6
    f^9      8                           same as f^5

(each group is a pair; "222|311" is a Casimir eigenvalue that two SO(7)
irreps share, so the label alone does not say which.)  This is the known
incompleteness of Racah's classification for f^N, not a defect of the method.
Inside such a group the diagonalizing builder produces an arbitrary
orthogonal rotation, and the transformation is still orthogonal and still
M-independent -- those checks do not validate the choice of alpha.  That is
what the CFP builder is for.

The overall sign of each basis vector is also a convention.  With the
diagonalizing builder it is fixed by making the largest component of the
highest-weight state positive.  Note that the sign must be fixed at the
highest weight, not after J coupling -- doing it afterwards makes the
transformation M-dependent, which is what the M-independence check is there
to catch.  With the CFP builder the LS-side signs are not a convention at
all: they are Cowan's, carried by the cfp values themselves.  The jj-side
signs remain ours in both cases; RCG has no true jj basis to match them to.


Cowan CFP mode (cowan_cfp.py)
-----------------------------
The LS states are built by the recursion that defines the cfp:

    | l^n a L S ; ML MS >  =  sum_ab  ( l^(n-1) ab | } l^n a )
                             * sum_(m,mu) < Lb MLb, l m | L ML >
                                          < Sb MSb, 1/2 mu | S MS >
                                          a+_(m,mu) | l^(n-1) ab ; MLb MSb >

starting from the closed-shell vacuum, using the same card deck RCG reads
(CUVFD, rcg11kd.f:2136-2171).  The deck is looked for as $LSJJ_CFP, then
./ING11, then ../CFP/ING11, then C:\cowan\FOR\g95\Claude\CFP\ING11.

Everything is blocked by (ML, MS), so each step is a handful of small dense
matrix products.  Building f^1 ... f^14 takes about 17 s in total.

Two identities are asserted at every step, so a misread deck cannot pass:

  * the vector built above has norm exactly sqrt(n), because
    < l^n aLS | a+ | l^(n-1) ab Lb Sb >^(LS) = sqrt(n) cfp and the cfp rows
    are orthonormal;
  * within one (ML,MS) block the states are orthonormal and their number
    equals the number of Slater determinants in the block.

What is checked, beyond the sweep above:

  * all 37 cases pass the full validation with the CFP-built LS side
    (validate_all_cowan.log), and the term lists with seniority still match
    the published tables for p^2, p^3, d^2, d^3, d^4, f^2, f^3;
  * the seniority recovered from < A+A > for each of Cowan's alpha labels
    reproduces all 31 entries of Cowan's own SENIOR file (d^3 ... d^7),
    which is an independent statement that alpha is attached to the right
    state;
  * compare_racah.py overlaps the CFP basis with the diagonalized one.  The
    overlap matrix never connects states of different (L, S, v) -- or, for f,
    different (W) and (U) -- to better than 8e-11, which is the precision of
    the 10-decimal cards.  Every p and d state, and every f state outside the
    groups listed above, is pinned one-to-one: the two constructions agree
    state for state up to sign.
  * the (W) and (U) of a CFP-built state are measured, not assumed, as the
    expectation values of the two Casimirs.  The r.m.s. deviation from an
    eigenstate stays below 2e-4 over every f shell, i.e. Cowan's alpha
    respects Racah's classification wherever Racah's classification decides
    anything.

The residual (-1)^(n-1) between "electron n in slot n", which is what the
cfp is defined against, and a+ acting on the left, which is what the code
does, is uniform over a shell and therefore drops out of every relative sign.
It is not applied.

The calculational machinery (recur.py)
--------------------------------------
The m-scheme route above is how the numbers were established: it has no
phase conventions to reconcile, so it can be trusted on its own.  It is not
how RCG should compute them.  f^7 alone has 3432 determinants, and the
m-scheme needs a null-space solver for the highest weight and an eigensolver
for the classification -- neither of which belongs in RCG.

recur.py computes the same matrix by recursion on N, from coefficients of
fractional parentage alone.  Its ingredients are exactly three:

  * Cowan's LS cfp deck, which RCG already reads (CUVFD, rcg11kd.f:2136);
  * a jj cfp deck for j = 1/2, 3/2, 5/2, 7/2 -- 90 cards, in the file JJCFP;
  * Clebsch-Gordan coefficients, i.e. what S6J and S9J already provide.

No determinants, no linear algebra beyond small dense matrix products.

Write T_n( a L S J ; ( j1^N1 c1 J1, j2^N2 c2 J2 ) J ) for the overlap.  Peel
one creation operator off one jj subshell and use

    a_(jm) | l^n a L S J M >
      = sqrt(n) sum_(ab,Jb) cfp(ab|a) Y(ab;j,Jb) < j m, Jb Mb | J M >
                            | l^(n-1) ab Jb Mb >,

    Y = < ( (l s) j , (Lb Sb) Jb ) J | ( (Lb l) L , (Sb s) S ) J >  (a 9j),

which gives

    T_n = sqrt(n/Nk) * sigma_k * sum_(ck-bar) cfp_jj( ck-bar | ck )
              * sum_Jb R_k(Jb) * sum_ab cfp_LS( ab | a ) Y( ab ; jk, Jb )
              * T_(n-1),

with R_k the 6j-type recoupling of the spectator subshell and sigma_1 = 1,
sigma_2 = (-1)^N1 -- the sign from carrying a+ of subshell 2 left past the
N1 operators of subshell 1.

The recoupling coefficients are evaluated here by explicit Clebsch-Gordan
sums rather than by 6j/9j formulas.  They are small, and a summed-CG
definition cannot carry a phase-convention error; a Fortran port substitutes
S6J/S9J and checks against these.

Three things check the result, and all three are reported by
`python recur.py check all`:

  * every J block of T is orthogonal;
  * peeling from subshell 1 and peeling from subshell 2 are each, by
    themselves, a complete rewriting of the jj state -- the two subshells are
    separate factors, so there is no sum over which subshell the last
    electron came from, as there is over which electron inside one shell.
    Both routes are computed whenever both subshells are occupied and
    compared.  This is a very sharp check: it is what caught the initial
    double counting, and it is what pins the (-1)^N1;
  * the result is compared element by element with the m-scheme matrix,
    after matching row and column labels.

All 30 shells agree with the m-scheme to ~1e-10, which is the precision of
the 10-decimal LS cards (recur_all.log).  f^7 builds in a few seconds.

The jj cfp deck is generated once, by projecting the m-scheme subshells:

    < j^n vJ M | sum_m < j m, Jb Mb | J M > a+_(j,m) | j^(n-1) vb Jb Mb >
        = sqrt(n) cfp,

each row checked to be normalized, which is the statement that the parent
list is complete.  `python recur.py deck JJCFP` rewrites it; adding the word
`fromdeck` to a check makes the run read the file instead of rebuilding it,
so the whole sweep can be run with no m-scheme input at all.

For l <= 3 the pair (v, J) names a jj subshell state uniquely, so no further
label is needed; the first case where it does not is j = 9/2 with n = 4,
which does not arise for the shells Cowan's decks cover beyond n = 2.

Porting to RCG
--------------
The transformation itself is the easy part.  The work in RCG is that a jj
basis state needs a different descriptor than the code currently has: RCG
identifies a state by one term index per subshell, NALS(L,I) / NALSJ(L,I),
and writes the label at line 7608 as (2S+1 L 2J) 2J.  A true jj state needs
(N1, N2, v1, v2, J1, J2) per subshell, so those arrays and the label writer
need a parallel jj-side description, and the basis generator has to enumerate
distributions rather than terms.

Also, CPL37 derives couplings 3-7 from this matrix (it reads it back off
tape 2 at line 8424) assuming the subshell terms are shared between the two
schemes.  Cleanest containment is to carry the true jj matrix as a separate
transformation under a KCPL value of its own and leave the existing one
feeding CPL37 untouched.

That new KCPL value is internal, and it does not displace Cowan's JJ.  The
first plan was for a plain KCPL=2 to select true jj for itself, but that
would have changed what an existing deck computes, and Cowan's LSJJ output
is parsed by other programs: OUTG11 here, LEVELS and OUTE out of RCE.  True
jj is therefore ADDITIONAL, asked for by KCPLD(8)=8 in the flag block of the
parameter line ING11 already uses -- see "Running it" below.  Nothing an
existing input deck says changes meaning, and nothing it produces moves.

Pure jj as computed here is sequential coupling of relativistic subshells.
It is worth having and it does not cover every real coupling; Gaigalas's
Coupling code does that, and can already be driven from the LS composition
tables in OUTG11.

Two output questions follow from a true jj basis, neither of which the
calculational side settles:

  * the labels.  Cowan's are LSJ.  A true jj state of p^4 reads as something
    like  p-2.0.p+2 , the 0 being the intermediate J of the p- subshell, and
    that does not fit the present label fields.
  * the leading configuration label, which RCE prints in LEVELS1 separately
    from the component percentages.  Cowan's configurations are all non-
    relativistic; in true jj each breaks into a set of relativistic sub-
    configurations, and the label ought to show that.  The names come from
    IN36, where they are arbitrary and not necessarily parseable, so they
    cannot be rewritten.  The only route left is to carry all the subshell
    notation in the TERM label and leave the nonrelativistic configuration
    label exactly as IN36 gave it -- presumably through some encoding, with
    the scheme printed at the end of the output so the user can parse it.
  * the genealogical data -- intermediate terms and quantum numbers -- needs
    a proper printed form in OUTG11 for jj and the related couplings, with
    minimal disturbance to the existing formats.

All of these are downstream of the matrix and can be taken up separately.


The Fortran side (rcg11kd.f)
----------------------------
The recursion is now in RCG itself, at the end of rcg11kd.f, as a block of
program units that touch nothing else in the file:

    LSJJIN(IU,IERR)          read LSJJ.DAT (once)
    LSJJU(l,N,IERR)          build U(l^N), all J, into /LSJJR/
    LSJJY, LSJJR1, LSJJR2    the recoupling coefficients
    LSJJKY, LSJJRC, LSJJFR   packed state key, route choice, deck lookup
    LSJJOR(DEV)              orthogonality residual
    LSJJCK(DEV,IERR)         the 1-vs-2 route residual
    LSJJSF(IU,IW,lmin,lmax)  self test and dump
    LSJJRT(TMX,LDT,NBLK,MODE,IERR)  rotate one block's matrix in place
    LSJJID(NLS)              name the rotated columns
    LSJJKP(IW)               print the key to those names
    LSJJCC(NC)               the code alphabet
    LSJJOP(I,JC,K)           is that subshell open?
    LSJJPAR.F                parameters and COMMON blocks (a new include,
                             alongside RCGPAR.F)
    LSJJCOM.F                the same without the IMPLICIT lines, for
                             routines that also include RCGPAR.F
    LSJJID.F                 the label arrays and the key table

Where recur.py evaluates the recoupling coefficients by explicit sums over
Clebsch-Gordan coefficients -- slow, but a summed-CG definition cannot carry
a phase-convention error -- the Fortran uses RCG's own S9J and S6J:

    Y  = (-1)^(j+Jb-J) sqrt((2Jb+1)(2j+1)(2L+1)(2S+1))
             * 9j{ Lb Sb Jb ; l s j ; L S J }
    R1 = (-1)^(j1+J1b+J2+J) sqrt((2J1+1)(2Jb+1)) * 6j{ j1 J1b J1 ; J2 J Jb }
    R2 = (-1)^(a+2b+3c+g-f) sqrt((2g+1)(2f+1)) * 6j{ a b g ; c J f },
             a=j2, b=J2b, c=J1, f=Jb, g=J2

The phases in front are exactly what makes this worth checking rather than
trusting, which is what cmpfort.py is for.  Angular momenta are carried
doubled, as integers, everywhere except inside S6J and S9J.

Data.  LSJJ.DAT, written by CODE/lsjj/mkdeck.py when make_cfp runs it,
holds Cowan's LS cfp (p, d, f for
every n; g, h for n <= 2) and the jj cfp (j = 1/2 .. 7/2 for every n; 9/2
and 11/2 for n <= 2).  It belongs in the system-file directory that
COWAN.CFG points at, next to SENIOR; rcg11kd.f opens it on unit 23 and says
nothing if it is not there.

It is an UNFORMATTED SEQUENTIAL file -- native INTEGER*4 and REAL*8 in
native byte order, 453 kB.  It began as list-directed text, which was handy
while the recursion was being written, but a megabyte of twelve-digit cfp is
not something anyone reads, and parsing it cost more than the rotation it
feeds.  Binary also carries the cfp exactly: the text form rounded to twelve
significant digits, which showed up in TAPE2E at the 1e-13 level, invisible
in OUTG11's five decimals.

Being binary makes it machine specific.  The first record carries a magic
number and a version, so a deck from the other byte order, a deck from an
older layout, or one still in text form is refused (LSJJIN IERR=3) rather
than read as garbage, and the message on OUTG11 says how to rebuild it.

RECOVERY.  LSJJ.DAT is derived, never edited and never repaired -- if it is
lost, corrupted, or moved to a machine of the other endianness, run

    make_cfp

which rebuilds it along with the CFP decks, in the system-file directory,
in a fifth of a second.  That is the supported route and the only one an
installed user has.  mkdeck.py depends on nothing but cowan_cfp.py (Cowan's
own LS cfp tables) and recur.py (the jj cfp recursion), both beside it in
CODE/lsjj/ and both deterministic; no earlier LSJJ.DAT is involved.
`python ../../CODE/lsjj/mkdeck.py -t deck.txt` writes the same records as
text when they need to be looked at, and `python cmpfort.py OUTG11` after a
KCPL=9 run checks the rebuilt deck end to end against recur.py.

It needs numpy, which is the one thing in this package that has to be
pip-installed.  Where numpy is absent make_cfp says so and carries on: the
CFP decks are built, the install succeeds, and everything except true jj
coupling works, RCG's own diagnostic covering the rest.

One trap worth recording, since it cost a debugging round: list-directed
input does not care how a list is split over lines, but an unformatted
record boundary is exactly where the reader says it is.  The LS cfp matrix
of l^n is ONE record of NT*NP values, not NT records of NP.  mkdeck.py and
LSJJIN have to agree record for record.

Running the self test.  A card with KCPL=9 at the head of ING11 -- a value
that used to be simply illegal, and 9 rather than 4 because 1 to 7 are the
representation numbers of RCG_DOC, where 4 is LSLK, so a card called KCPL=4
reads as a request for LSLK -- builds U for every shell the deck covers,
writes the matrices and the two residuals on OUTG11, consumes the card and
lets the run continue.  Then

    python cmpfort.py OUTG11

matches rows on (term index, 2J) and columns on (N1, v1, J1, N2, v2, J2, J)
and compares element by element with recur.py.

What it gives.  34 shells, 34 pass, from inside RCG: every J block
orthogonal to 2.6e-10, the two peeling routes agreeing to 1.5e-10 -- both at
the precision of the 10-decimal LS cards -- and agreement with recur.py to
1.0e-13, i.e. to double-precision round-off.  The whole sweep, p^1 through
h^2, takes about half a second.

The rotation (LSJJRT)
---------------------

CALCFC builds Cowan's J1J2 matrix as it always did.  LSJJRT then turns it
into the true jj matrix in place, when asked:

    TMX_jj(.,c) = sum_c' TMX(.,c') * prod_i U_i( t_i(c'), c_i(c) )

The rotation touches columns only.  It leaves every subshell J and the whole
intermediate J chain alone, so at fixed Ji the old labels (an LS term of
subshell i) and the new ones ((j-^N1 v1 J1, j+^N2 v2 J2) Ji) are two bases of
one space and the dimensions never change.

How it is applied.  The columns of one CAFCDOLP block are a Cartesian
product over subshells -- subshell 8 varying fastest -- so CAFCDOLP now
records each block's shape in /LSJJBK/ and LSJJRT walks one subshell at a
time.  A fibre (the columns differing only in the term of subshell i) is
multiplied by the block of U_i belonging to that subshell J.  Nothing bigger
than one fibre is ever held, which matters at KMX=20000.

Matching RCG's terms to the deck's.  The two term lists hold the same terms
in different orders -- RCG reads Cowan's binary cfp deck, LSJJ.DAT was built
from the source tables.  f^3, for instance, is 4I 2H 2H 2F 2F 4S 4F 4G 2K 2D
2G 2P 2G 4D 2L 2I 2D in the deck against 4I 4G 4F 4D 4S 2L 2K 2I 2H 2H 2G 2G
2F 2F 2D 2D 2P in RCG.  So LSJJRT builds the permutation explicitly, keyed on
(2S, 2L) and, where that repeats, on the one non-blank character of Cowan's
alpha field -- first non-blank rather than a fixed column, so it does not
care how either side pads.  The result must be a bijection.  It is not a
convenience check: the NLT11 term cutoff leaves RCG with fewer terms than the
deck, and then no rotation exists at all, which this catches (IERR=5).

Subshells left alone.  l=0, or N<=1, or N>=4l+1.  Those are exactly the cases
where RCG itself carries a single made-up term instead of a deck term list,
and every subshell J block is one-dimensional there, so the two bases
coincide.

Running it.  Column 38 of the main control card -- KCPLD(8) -- set to 8:

        1   -182121818182         0011111 1 ...
                                       ^ column 38

True jj is an ADDITIONAL coupling, not a replacement.  Cowan's LSJJ (J1J2)
output has to stay exactly as it is, because other programs parse it: OUTG11
here, LEVELS and OUTE out of RCE.  So a KCPLD(8)=8 run prints everything it
printed before, unchanged to the byte, and adds an eigenvector section and a
purity line to OUTG11 and a trailer to TAPE2E.

Why 8 and not 0.  KCPLD(I)=0 means "do coupling I" for I=1 to 7, and reading
an eighth column the same way would have been the obvious thing -- but 0 is
what a blank field gives, so every deck ever written would have switched
this on.  Couplings 3 to 7 got away with 0-means-yes only because NLSMAX
gated them as well -- no longer true; see "Blank NLSMAX is a permission,
not a prohibition" at the end of this file.  KCPLD(8) is not unused
either: ENERGY reads it as IKKK, which acts on the value 9 (and on
nothing else).  So 9 is out and the sense has to be inverted.  8 also
reads as "coupling 8", the number this coupling carries internally:
KCPL1=8 in SPRIN, COUPL(8), AVP(8,.).

There is no longer a KCPL=8 card; the only card in that range is the KCPL=9
self test.  (9, not 4: 1 to 7 are the representation numbers listed in
RCG_DOC, where 4 is LSLK, so a card called KCPL=4 reads as a request for
LSLK to anyone who knows that list.)

LSJJ.DAT must sit with SENIOR and the binary cfp decks in the directory
COWAN.CFG points at (c:\cowan\code in this installation).

What it gives, on Pr III.  With KCPLD(8) blank the run is byte-identical to
the reference -- OUTGINE, TAPE2E and OUTG11 all -- so nothing existing moved.
With it set to 8: normal exit, no LSJJRT error on any block, and RCG's own
orthonormality test passes on every rotated matrix in the run, which is a
real test rather than a formality (an orthogonal TMX times a block-diagonal
orthogonal U must stay orthogonal, and a mismatched fibre, a wrong term
permutation or a wrong J block would break it at 1e-7).  OUTGINE is
byte-identical to the unrotated run: the energies do not move, which is what
an orthogonal change of basis must do.

Correctness rests on four things that are each checked per block rather than
assumed: the term permutation is verified to be a bijection against RCG's own
(2S,2L,alpha) list for that subshell, which would fail if the wrong shell's U
were loaded; the U block dimension must equal the fibre size for that
subshell J; the block must be the full Cartesian product its recorded shape
implies, which catches KMX overflow; and the rotated matrix must stay
orthonormal.  The U matrices themselves were validated to 1e-13 against
recur.py beforehand.


Where the two matrices go
-------------------------

Nothing that existed before the rotation reads the rotated matrix.  Cowan's
J1J2 matrix goes everywhere it always went, and the true jj matrix travels
beside it on a scratch unit of its own, read only by the code that produces
the new output:

  - CALCFC builds the J1J2 matrix, runs the orthonormality check and the
    ICFC printout on it exactly as before, and only then saves it to unit 43,
    checks the rotation (LSJJRT MODE=0), rotates (MODE=+1), rechecks
    orthonormality, names the columns (LSJJID), writes the matrix and the
    labels PK1-packed to unit 42, and restores the saved copy from unit 43;
  - units IL (31/32) and IC (41) therefore carry Cowan's matrix, and CPL37,
    SPRN37, ENERGY, the diagonalization and the transition arrays are
    untouched;
  - SPRIN, after its LS and JJ passes, makes a third pass in true jj: it
    reads the block's matrix from unit 42 and re-expresses the eigenvectors
    it already has, C_jj = Ujj' C_LS.  That adds an EIGENVECTORS (TRUEjj
    COUPLING) section and accumulates AVP(8,.), printed on a line of its own
    so that the AVE PURITIES table keeps its shape;
  - RCEINP writes every TAPE2E record it always wrote, and then appends the
    true jj trailer described below.

Restoring from unit 43 rather than rotating back with MODE=-1.  MODE=-1 is
the exact inverse in exact arithmetic -- U is orthogonal, so it is U
transposed -- but the products are still floating point, and the round trip
came back changed in the last bits.  That was enough to move TAPE2E, which
would have made the ordinary results of a run depend on whether the extra
output had been asked for.  A saved copy makes it exact.  MODE=-1 is still
there and still tested by the round-trip check below; it is simply not used
in CALCFC any more.

Unit 42 rather than a second in-memory buffer because TMX is KMX**2 REAL*8
-- 3.2 GB at KMX=20000, on an allocation set that is already about 19 GB.
Unit 42 is written once per J block by CALCFC, in the order the blocks are
produced, and read back one record per block by whichever consumer runs --
SPRIN, or the trailer in RCEINP, never both, since RCEINP is the alternative
branch to the diagonalization.  It is rewound wherever IC is.  A block whose
rotation failed still writes its record, marked IOK=0: the readers are
sequential and have no way back, so a skipped record would put every later
block on the wrong matrix.

Is PK1 worth it on the rotated matrix?  Yes, and not marginally.  PK1 stores
12 bytes per nonzero -- a REAL*8 and two INTEGER*2 indices -- against 8
bytes per element dense, so it stops paying at 66.7 per cent density.
Measured on Pr III by logging every PK1 call, rotation off vs on:

    unit 41 (IC)       0.16%  ->  0.20%     14.0 MB packed, 4680 MB dense
    unit  3 (TMX)      2.62%  -> 13.11%      3.3 MB packed,   16.9 MB dense
    units 31/32 (IL)   2.4/2.7% -> 8.8/14.2%
    unit  2 (TAPE2E)   0.14%  ->  0.18%     12.6 MB packed, 4668 MB dense

The rotation costs a factor of about five in density -- as it must, since
the J1J2 matrix is nonzero only where every subshell term index agrees and a
fibre rotation destroys exactly that sparsity -- but the rotated matrix is
still 87 per cent zeros, five times below break-even.  (PK1's INTEGER*2
indices also cap any dimension at 32767, which KMX=20000 fits, but not by
much.)


The labels
----------

The columns of the rotated matrix are jj states, and printing them under LS
term names would have made the new section useless.  A jj basis state is,
for every subshell, its jj state and its own J, plus the chain of
intermediate J's that couples them -- and spelt out that does not come near
the eight characters the eigenvector formats allow, since one subshell of
p^4 alone reads p-2(0) p+2(2).  So it is coded, one character per quantity,
and the key goes into a file of its own, true-jj_keys.txt.  With n the number
of OPEN
subshells (empty and closed ones carry no quantum numbers at all and are
left out):

    n=1   A1 J1
    n=2   A1 J1 A2 J2
    n=3   A1 J1 A2 J2 K2 A3 J3
    n=4   Ai Ji only, the chain dropped
    n>=5  Ai only

Ai is the jj state of subshell i -- its (j-^N1 v1 J1, j+^N2 v2 J2) -- keyed
per configuration and subshell.  Ji is 2 * the total J of subshell i and Ki
is 2 * the running J after coupling it in, both coded directly: character 1
is 2J=0, 2 is 2J=1.  The first and last running J are never coded because
they are never free -- the first is J1, the last is the J of the block --
which is what makes n=3 fit.  For n<=3 the label names the basis state
completely: within one configuration and one J block, two columns with the
same label are the same column.  On Pr III every one of the 5625 labels in
the 24 jj sections is unique, and no configuration reaches n=4.

Leaving closed subshells out is what makes that work.  Pr III's
configurations declare eight subshells and open two or three, but several of
them close a p^6; counting it took f6p7p to n=4 and cost 246 duplicate
labels.  A closed subshell has one state, J=0, and leaves the running J
alone, so its field would have said nothing.

The alphabet is 1..9, then A..Z, then a..z, then the punctuation listed in
LSJJCC -- 85 codes.  '0' is left out so a code is never taken for a blank or
a zero, '*' means the key ran out, and the brackets and the space are left
out because the formats around the label use them.  The alphanumerics alone
would cover every shell in the deck except f^5 through f^9, whose (N1 v1 J1,
N2 v2 J2) count peaks at 70, so the punctuation is not decoration.  Keeping
Ji out of Ai, in a character of its own, is what keeps the key short: folded
in, one subshell's codes would run to its number of levels -- about 300 for
f^6 -- instead of its number of (N1 v1 J1, N2 v2 J2).

LSJJID builds the labels and the key table; LSJJKP prints the key; LSJJCC is
the alphabet and LSJJOP decides what counts as open.  LSJJID.F carries the
layout and the COMMON blocks.

The key file.  LSJJKP writes true-jj_keys.txt in the run directory (unit 24)
and leaves OUTG11 nothing but a line saying where the key is.  The key is
needed for as long as the level lists are kept, and OUTG11 is big and is
usually the first file deleted.  Nothing is written unless a block was
rotated, so a run without KCPLD(8)=8 produces no such file.

An entry reads

     1: (f-1 J= 2.5, f+2 J= 2.0)

-- the code, then the j- part and the j+ part, each naming the subshell in
lower case, how many electrons are in it, and its J.  A part with no
electrons is left out, so a state with all the electrons in j+ reads
(f+3 J= 1.5).  A seniority is shown, as ` v2`, only where J alone does not
tell the states of that j^N apart: on Pr III that is four entries out of
several hundred, the two pairs of f+^4 J=2 and J=4.  Seniority is the right
extra quantum number here -- the Nielson-Koster indices that replace it in LS
coupling are defined for LS terms of d^N and f^N and have no meaning in jj --
and for the subshells the jj deck carries (j<=7/2 for any N, 9/2 and 11/2 for
N<=2) the pair (v,J) is always unique.

A subshell the rotation left alone -- l=0, one electron, or one hole -- is
skipped by LSJJRT not because it does not split into j- and j+, but because
there is nothing to rotate: its J1J2 state already is a jj state.  An s
subshell has only one j.  One electron has j = J of the subshell.  One hole
sits in the part whose j equals J, the other part being full and carrying
J=0.  So the partition is not stored in NJJST, but it is determined, and
SUBROUTINE LSJJ1E works it out from (l, N, 2J).  Such a subshell has exactly
one LS term and so one code, which therefore covers both values of J, and its
entry gives both, each against the J it belongs to:

     1: J= 0.5: (p-1 J= 0.5, p+4 J= 0.0)
        J= 1.5: (p-2 J= 0.0, p+3 J= 1.5)

The Ji character of the label says which one applies.  (An earlier version of
this key printed `not split` here instead, which read as though the subshell
had no j- / j+ decomposition at all.  It has one; it is just implied by J
rather than stored.  GENPRT prints it in full, LSJJ1E being shared by both.)
These are the only cases LSJJRT skips, and each has exactly one LS term,
which is why the key no longer prints one: the old code took the term name
from /C2/, which by the end of the run holds the last configuration read, not
the one being keyed.


What RCE gets
-------------

A trailer on the end of TAPE2E, after everything RCE has ever read from unit
2, written by RCEINP at label 850.  A reader that stops when it has what it
came for never sees it, so a TAPE2E from a KCPLD(8)=8 run is still a valid
TAPE2E for an RCE that predates all of this.  The layout is

    IJJM2E, IJJV2E, NBLOCKS              magic 20260802, version 1

then, per J block, in the order CALCFC produced them:

    IOK, J, NLS, NOSUBC,
      ((NJJST(1..4,c,i), i=1..NOSUBC), FIDJ2(c), c=1..NLS)
    the transformation matrix, PK1-packed as everywhere else

IOK=0 marks a block whose rotation failed; its matrix is Cowan's J1J2 one
and its labels are blank.  NJJST(1..4,c,i) is N1, the level of j-^N1, N2,
the level of j+^N2 for subshell i of column c, or -1 where the column still
carries an LS term.  FIDJ2(c) is the eight-character label described above.
The version number changes if the layout ever does; RCE should check both it
and the magic before believing what follows.


The RCE side
------------

RCE20KD.F reads that trailer and prints the fitted levels in true jj as
well.  It is switched by the file and by nothing else: RDJJ2E, called
before anything is read from unit 2, looks for the magic record and sets
IJJON, so a TAPE2E written without KCPLD(8)=8 gives exactly the run it
always gave and no input card changes.  Finding the trailer costs one
pass over TAPE2E in which three integers are taken from each record --
reading fewer items than a record holds skips the rest -- so it is
twelve bytes per record, not the file.

The trailer covers every cset of the RCG run in the order CALCFC made
the blocks, which is the order they appear on TAPE2E, so one cset is a
run of consecutive trailer blocks.  MAPJJ2 finds which run by matching J
and matrix size for every block of the cset at once, starting where the
last cset ended.  A cset it cannot place goes without the true jj output
and says so.

What it then does is what RCE already does for coupling 2.  The
eigenvectors of the fit are re-expressed by the trailer's matrix, which
is <LS|TRUEjj> and therefore acts on vectors in the LS basis; ORDOUT2 is
called with it in place of the J1J2 matrix from TAPE2E and needs no
change at all.  The pass runs third, after the two couplings RCG wrote,
and takes the coupling-1 vectors again rather than the coupling-2 ones
the second pass leaves behind - IEVV2/NW2 hold them, which is free here
because the ITRANS=2 case those arrays belong to needs KCPL(2)=1 and
cannot arise alongside true jj.  If RCG's first coupling was not LS, the
matrices on TAPE2E are in some other basis already and the pass is
skipped with a message.

The level lists go to LEVELS1J, and SORT8TO7 - now a wrapper around
SRTLEV, which takes its three units and the number of sections per cset
as arguments - sorts them into LEVELS2J and LEVELS3J the same way it
sorts LEVELS1 into LEVELS2 and LEVELS3.  A third section per cset in
LEVELS1 itself would have interleaved new lines with old ones and given
SORT8TO7 six sections where it expects four, so the true jj list has
files of its own and LEVELS1, LEVELS2 and LEVELS3 come out untouched.
OUTE does get the new sections, between the coupling-2 section of a cset
and the start of the next one, as OUTG11 gets its TRUEjj section between
the JJ and JJJK ones.

The labels are the eight-character codes described above, whose key RCG
writes into true-jj_keys.txt.  RCE reads OUTGINE and TAPE2E, neither of
which carries the key, so it names that file rather than repeating it --
which is the point of having the key in a small file of its own.


Verified on Pr III
------------------

With KCPLD(8) blank: TAPE2E and OUTGINE byte-identical to the reference, and
OUTG11 identical apart from timing lines.

With KCPLD(8)=8: OUTGINE still byte-identical; TAPE2E identical over the
whole length of the reference file, with the trailer appended after it; and
OUTG11 identical apart from one line -- the echo of the input card, which now
shows the 8 -- with 208266 lines added and none removed or changed.  The AVE
PURITIES table is untouched (LS 0.381, JJ 0.423) and the true jj purity,
0.368, is on the line below it.  Every label is unique, no key overflowed,
and no orthonormality, LSJJRT or out-of-step message appears anywhere in the
run.

A weighted L1 checksum of TMX taken before the rotation and after MODE=-1
agrees to 1.3e-10 relative over all 25 J blocks, which is ordinary
double-precision accumulation noise on sums of 1e5 to 1e6 terms -- and is
the reason CALCFC restores a saved copy instead.

Couplings 3-7 are exercised too.  They used to be out of reach here --
NLSMAX is an I3 field and every Pr III block is larger than 999, so SPRN37
and CPL37 returned at once and those purity columns read 0.000 whatever was
asked for.  NLSMAX=999 now means KMX, the way Cowan elsewhere lets the
largest value a field can hold stand for the array bound, so the case runs
with all seven representations (105 s against 25 s) and:

  - with true jj off, OUTGINE and TAPE2E are byte-identical to the
    KCPLD=0011111 reference, and the only lines of OUTG11 that change are
    the echo of the card and the five AVE PURITIES rows, whose 3-7 columns
    now carry real numbers (0.395 0.385 0.405 0.393 0.415 on 4f3 - f4p5)
    while the LS and JJ columns stay at 0.381 and 0.423.  Nothing is
    deleted; the seven-representation sections are added.
  - with true jj on, it is still additive: OUTGINE byte-identical, TAPE2E
    identical over the whole reference with the trailer appended, and
    OUTG11 with 205972 lines added, none deleted, and nothing changed but
    the card echo and two elapsed times.  The TRUEjj section falls between
    the JJ section and the 3-7 sections, since SPRIN makes its third pass
    before SPRN37 runs.


RCE verified on Pr III
----------------------

The TAPE2E and OUTGINE of the KCPLD(8)=8 run above, read by RCE with
INV1=0 (build RCEINP from OUTGINE) and again with INV1=1 (read the
RCEINP so built).  Normal exit both times, 5 s, 25 blocks found.

  - against the same RCE built without any of this and run on the
    TAPE2E without a trailer: LEVELS1, LEVELS2, LEVELS3, RCEINP, RCEOUT
    and PARVALS byte-identical, OUTE identical apart from six elapsed
    times.  The new code is inert when there is no trailer.
  - with the trailer: LEVELS1, LEVELS2, LEVELS3, RCEINP, RCEOUT and
    PARVALS byte-identical to the run without it, and OUTE identical
    line for line once the two TRUEjj sections - 279779 lines - are
    taken out again.  Nothing is changed or deleted.
  - the compositions themselves are RCG's.  The fit ran zero iterations,
    so RCE's eigenvectors are the ones RCG diagonalized, and they agree
    with the TRUEjj sections of OUTG11 label for label: the lowest even
    J=1/2 level, 50 pct 4f3 "12" and 48 pct 4f3 "22" in LEVELS1J,
    against -0.70425 and 0.69000 in column 1 of the first TRUEjj block;
    the lowest odd J=1/2 level, 74/12/6 pct on f25d "3514", "1514",
    "2514", against 0.85844, -0.34130 and 0.23470.  The second cset
    matching as well is the test of the block mapping.
  - true jj is the less pure scheme here, as it is in RCG: averaged over
    all 5626 levels the largest component carries 0.296 in LS and 0.264
    in true jj.


Parsing of the labels -- the code that turns a label plus its key entry into
ordinary spectroscopic notation, e.g. 1s.2p-.(1).2p+ 5/2 -- is deliberately
left out of the Fortran package and will be written separately, after the
Fortran codes are finished.  The key file is meant to be its input, which is
why it is small, one entry per line, and printed in lower case.


The leading configuration label (made coupling-independent)
----------------------------------------------------------

RCE names a level by the configuration that carries the largest share of the
eigenvector: it sums BUFY(JOFS+I)**2 into VC(NC(I,M)) over the basis states I
of the J matrix, takes the largest, and prints CONFIG(NJM,K3).  NC(I,M) is the
configuration of basis state I in the FIRST coupling; the sum is correct in a
rotated basis as long as the rotation does not mix configurations and does not
move a state out of its index range.

Both rotations satisfy that.  A change of coupling recouples angular momenta
inside a configuration and cannot carry one configuration into another, so its
matrix is block-diagonal in the configuration.  That is not only the argument
but a measurement: instrumented builds swept both matrices RCE uses -- the
J1J2 one off TAPE2E and the true jj one off the trailer -- over all 25 J
blocks of the Pr III case, and the largest element outside the configuration
blocks is exactly zero in every one of them.  For LSJJRT the same fact is
visible in the code: it works one CAFCDOLP block at a time, each block belongs
to a single configuration (found at its head from the term index, as at label
349 of CALCFC), and the rotation multiplies fibres of columns inside that
block only.  No column ever leaves its block.

So the sum over a configuration is an invariant of the level, and the name
cannot depend on the coupling scheme.  It did.  On the Pr III case the LS pass
and the J1J2 pass disagreed on 137 of the 3781 odd-parity levels, and the true
jj pass followed J1J2 on all 137 (the even cset agreed throughout).  The cause
is not the rotation but what each pass had left to sum: the LS pass sums the
eigenvector as the diagonalizer left it, whereas every later pass sums the
vector ORDOUT2 rebuilds out of IEVV and NW -- the MAXC=50 largest components
only, each rounded to 1.D-4, IEVV being INTEGER*2 -- so the configuration sums
come out short, and short by different amounts.  Where two configurations run
close, that decides the winner.  All 137 are strongly mixed levels, e.g. the
odd J=1/2 level at 181.219 cm-1, whose largest single component is 16 pct.

This matters beyond the name.  SRTLEV groups LEVELS3 by the printed name, so a
level that changed name changed group, and the number of levels listed under a
configuration depended on the coupling scheme -- on the odd cset, p5f3d 2505
levels in LS against 2369 in J1J2, p5f3s 226 against 329, fd7p 86 against 118,
f26g 160 against 161.  Those lists are read to compare the average purity of a
configuration between coupling schemes, which cannot be done while the schemes
do not agree on which levels belong to it.

Fixed in RCE20KD.F by computing the name once, for every level of every block,
in a loop placed just after the ESORTW reordering and before the first output
pass -- that is, from the full-precision eigenvectors, all of their components,
and before anything has been transformed.  The values go into NJMS(J,M)
(INTEGER*2, IA by IX) and every pass reads its NJM from there instead of
summing its own.

Effect on the Pr III case: the LS sections of LEVELS1, LEVELS2 and LEVELS3 are
byte-identical to what they were, as are RCEOUT, RCEINP and PARVALS; the J1J2
and true jj sections change for exactly the 137 levels that were wrong; and
all three schemes now list the identical set of levels under every
configuration, in both parities.


The genealogy printout (IPGEN, columns 83-84 of the main card)
--------------------------------------------------------------

A basis-state label has eight characters, which is room for the term of one
subshell and the total.  The rest of what defines the state -- the term of
every other subshell, and the chain of running L, S and J that couples them --
is in /C3A/ and was never printed except by the ICFC/IPT dump at label 405 of
CALCFC, which covers some subshell counts and not others and overprints two
lines with 1H+.  Without it an eigenvector composition cannot be read back to
the states it names.  SUBROUTINE GENPRT, at the end of rcg11kd.f beside the
rest of the LSJJ block, prints it properly:

    IPGEN = 0 or blank   nothing.  A blank field is what every existing deck
                         has, so nothing new is printed unless asked for
            1            the tables go into OUTG11
            2            the tables go into genealogy.txt, and OUTG11 gets a
                         pointer to it, as it does for the label key

Columns 83-84 were read into IDUM and thrown away.  IPGEN sits at the end of
/C1/ in rcgpar.f, so nothing already in that block moves.

RCN2 punches IPGEN=2 on every control card it writes (RCN\RCN2K.F, the WRITE
at label 185 and FORMAT 18, where the blank 6X over columns 79-84 became
4X,I2).  So the genealogy comes out by default, in a file of its own, without
depending on ICFC or on anyone remembering to punch the column; the value is
set on one line just above the WRITE, if it ever has to change.  Nothing else
on the card moves: OUT2 gets the same field, from the same FORMAT, and the
card is otherwise identical to what RCN2 wrote before.

One table per coupling per J block, one line per basis state, numbered as the
eigenvector components are:

    LS       (term index, 2S+1, L, Cowan's alpha), then the running 2S+1 and
             L after coupling that subshell in
    J1J2     (term index, 2S+1, L, J of the subshell), then the running J.
             Cowan's jj, in which a subshell is still named by an LS term
    true jj  (j- part, j+ part), then 2J of the subshell and 2K, K being the
             running J.  A part reads as the subshell, its occupation, a
             seniority where J alone does not tell the states of that j**N
             apart, and 2J.  Every J in this table is doubled, as in the
             label key; the other two tables print J itself

A field is one open subshell.  Empty and closed subshells are left out, as
they are in the labels, and by the same test (LSJJOP), so the table and the
labels agree on which subshell is which.  The configurations of the block are
listed above it by their occupations -- the names of IN36 cannot be printed
there, ELEM not being filled until ENERGY, which runs after CALCFC -- and the
heading names the configuration set, because a configuration number runs
within its own set.

A subshell with l=0, with one electron, or with one hole is printed split all
the same, its two parts worked out from J by LSJJ1E (see the key section
above), so that every field of the table reads the same way.  Os48+ 5p^5 4d,
for instance, comes out as (p- 2  0,p+ 3  3) and (p- 1  1,p+ 4  0) against
(d- 1  3) and (d+ 1  5) -- the four relativistic families a jj code would
list for that configuration.

The subshell J and the running J of a true jj line repeat those of the J1J2
line above it.  That is not a mistake, and it is the same fact the leading
configuration label rests on: LSJJRT mixes only the subshell terms of columns
that already share the whole chain of J.  What the rotation replaces is the
description of the subshell itself.

The LS and J1J2 tables come from CALCFC once the labels are formed and
de-duplicated; the true jj table comes from the IJJROT block, after LSJJID
has named the columns and before any MODE=-1, which clears NJJST.  Fields go
into strings and not straight into the format because a field opens with a
bracket: with the bracket in the format, a line of fewer fields than the
format holds would end on the bracket of the field that never came.

Checked on the Pr III case, all 5626 basis states of both csets:

  - TAPE2E and OUTGINE are byte-identical to a run with IPGEN=0, and a run
    with KCPLD(8) blank still reproduces TAPE2E.ref and OUTGINE.ref exactly.
    OUTG11 differs by the two announcements and nothing else.
  - every table lists every state of its block once, numbered 1..NLS.
  - in the true jj table, every Ji character of the label equals the printed
    2J and every Ki equals the printed 2K, and within one (cset,
    configuration, subshell) the Ai character and the printed jj state are
    in one-to-one correspondence -- which is what the key file says they are.
  - in the LS table the running 2S+1 and L of the last open subshell equal
    the total term in Cowan's own label, on every line.
  - the running J of the last open subshell equals the J of the block, in
    both the J1J2 and the true jj tables, on every line.

genealogy.txt is 1.5 MB for that case, against 122 MB of OUTG11: the tables
are one line per state, not one per eigenvector component.

The unrotated subshells were checked separately, LSJJ1E being exercised over
every case it can meet (l = 0..3, one electron and one hole, both values of
J): the two parts always account for exactly N electrons and carry the J's
above.  End to end, RCN2 was run on the Os48+ 5p^5 4d deck and RCG on the
ING11 it wrote, with no column patched by hand: the genealogy came out
because of the punched IPGEN, and the p^5 fields read as quoted above.


===========================================================================
Speed and memory
===========================================================================

Everything above is about a new output.  What follows is about the cost of
producing it, and of producing everything else the two programs already
produced.  The two subjects met on one case: Gd III, whose largest J blocks
run to about 16000 basis states, took roughly twenty hours in RCG and spent
almost all of it in three loops that were written when a block of six
hundred was large.  None of the changes below alters a result.  Every one of
them is verified the way the LS -> jj work is -- the Pr III case run before
and after, with TAPE2E and OUTGINE compared byte for byte.

The general shape of the work is worth stating once, because it recurs:
where a loop is cubic and its operands are dense, it is threaded; where the
operands have structure the loop does not know about, the structure is used
instead, and the cost then falls by a factor no number of cores could reach.
CALCFC is entirely of the second kind.


Threading, and what it costs to be safe
---------------------------------------

RCG is built with -fopenmp added to the usual flags:

    gfortran -fshort-enums -malign-double -ftracer -fno-backslash
             -fallow-argument-mismatch -fcray-pointer -fno-automatic
             -O2 -fopenmp -o rcg11kd.exe rcg11kd.f

Without -fopenmp every !$OMP line is a comment again and the program builds
and runs serially, giving identical results.  That is not an accident of the
directives: no parallel region added anywhere in either program reassociates
a sum.  Each thread owns whole rows or whole columns, so the arithmetic a
result is made of is the same arithmetic in the same order at any thread
count, and the threaded executables are bit-identical to the serial ones.
This is why the threading could be applied to an eigensolver at all.

-fno-automatic must accompany -fopenmp.  -fopenmp implies -frecursive, which
moves locals from static storage to the stack, and parts of both programs
rely on locals keeping their values between calls; gfortran warns that
-fno-automatic overrides -frecursive, and that warning is expected.  A
consequence is that no routine called from inside a parallel region would be
thread safe, since its locals are shared.  None of the parallel regions
calls one.  That is also why the S6J build loops in CPL37 are left serial:
S6J is mathematically pure, but its locals are static.

-static is not possible with OpenMP on this toolchain, so an OpenMP build
must ship libgomp-1.dll, libgfortran-5.dll, libquadmath-0.dll and
libwinpthread-1.dll beside the executable.  OMP_NUM_THREADS sets the thread
count; unset, OpenMP takes every logical processor.  Do not add -ffast-math:
it changes the arithmetic of the eigensolver.  -march=native is not used for
the distributed executable and buys about one per cent, the eigensolver
being memory-bandwidth bound.

What was threaded in RCG: the Householder reduction and back-transformation
in TRED2 and the QL sweeps in TQL2 (batched through TQLROT), taken over from
RCE20KD where they were already in service; the six coupling-transform
products in SPRIN; the eigenvector transform at label 200 of SPRN37 (2.81 s
-> 0.63 s in the all-couplings test); both cubic products in CPL37; the
line-selection loop of SPECTR, which is O(NLS*NLSP**2) and was the slowest
single part of a large run -- its dot products are now computed one row at a
time, and CTA/CT4 are stored transposed so that both operands are walked
contiguously; and the eigenvector orthonormality check in ENERGY, which is
now a threaded scan that only counts bad pairs, the original printing loop
running only if the scan finds any.

Two memory faults turned up on the way.  In SPRIN and MUPOLE the POINTER
statement for TMXP was spelled TMPX, so TMXP was never mapped onto the C3LC4
block and survived as an 800 MB static array; fixing the typo cut .bss from
890 MB to 90 MB.  And SORT2's insertion sort, O(N**2) with N up to
KLAM=500000, is now the stable merge sort MSORTX, with ORDERI, ORDER2 and
ORDERCH no longer putting their workspace on the stack.

Deliberately not done, and why:

  - in SPECTR the E0MAX and DMINGA cuts are applied after the dot products,
    so about 2.4 M pairs per Pr III run are computed and thrown away.
    Hoisting them would save a great deal more but would change the printed
    SUMS2/SUMS2N sum rules;
  - the eigenvector print loops in SPRIN and SPRN37 are now the largest
    single item in an all-couplings run, about 8 s of 20.  The cost is
    gfortran's F9.5 conversion.  Buffering and threading were both tried and
    both lost; only a hand-written converter would help, and it would have
    to reproduce gfortran's rounding exactly.


CALCFC: the cubic loops
-----------------------

Three loops in CALCFC were cubic in NLS.  All three are gone, and none was
replaced by a faster approximation of itself: in each case the matrix or the
list being scanned has a structure that makes most of the work provably
empty, and the replacement does only the rest.

De-duplicating term names -- CFCDUP.  The loop at label 4715 looked for the
first pair L<LP with the same configuration and the same label, advanced one
disambiguating character, and then started the whole double scan again from
L=1.  On Pr III it restarted about 0.8*NLS times per J block, so it cost
about 0.4*NLS**3 comparisons: fine at NLS=600, hours at NLS=16000.  Only
levels sharing a configuration and a label ever meet in that test, so the
work splits into independent classes.  CFCDUP sorts once, which brings each
class together, and then runs Cowan's original loop -- body, restart and all
-- inside one class at a time.  The renaming is deliberately not short-cut
with a formula: it looks like plain lettering but is not, because a blank in
column 7 turns into 'a' on the left member of a pair and so stops matching
the blanks still waiting behind it, which leaves the last member of an
odd-sized class blank.  Running the original loop keeps every such detail by
construction.  The classes are independent only if no renamed label lands on
a label of another class, which cannot happen with the labels the format at
label 41 produces but is not worth assuming: the labels are saved first, and
if a second sort finds any duplicate left, everything is put back and the
original loop runs over the whole block.  That fallback has never fired.

Orthonormality of Cowan's J1J2 matrix -- CFCORT.  The triple loop at label
540 ran over every pair of rows and summed over every column, so it cost
NLS**3/2 multiply-adds, and it read TMX along rows of a column-major array,
one cache line per term.  Measured alone on a 7950X: 11 s at NLS=2000, 138 s
at NLS=4000, 1427 s at NLS=8000 -- hours per block at NLS=16000, which is
where most of the twenty hours went.  But TMX is block diagonal.  The loop at
label 505 leaves TMX(LS,JJ) at zero unless NALS(LS,I)=NALSJ(JJ,I) for every
subshell I, so a row has nonzeros only in the columns carrying its own
subshell-term signature, and two rows of different signature have disjoint
supports and a product that is exactly zero term by term.  Grouping the rows
and the columns by signature and summing only inside a group drops the cost
to the sum of N(g)**2*M(g), which is linear in NLS in practice.

The result is bit for bit the old one, and by construction rather than by
measurement alone: the skipped terms are exactly 0.D0 (verified by
instrumenting a run and counting the nonzeros outside the groups -- none),
adding 0.D0 changes no bits, the surviving terms are still summed in
increasing column order because the sort is stable, and the pairs are still
visited in increasing L and then increasing LP, so any message comes out
where it used to.  If the signature encoding overflows, CFCORT falls back to
the original full check.

Orthonormality of the rotated true jj matrix -- removed, and replaced by
CFCUOR.  The rotation used to be followed by a second full check, of the same
NLS**3/2 form, run on TMX*U.  It could never find anything but a bad U: the
matrix going in has just been checked at label 540, and

      (T*U)*(T*U)' = T*(U*U')*T' = T*T'      when U*U' = 1,

so orthonormality of the product is equivalent to orthogonality of the
factor.  CFCUOR tests the factor.  What that buys is the size of the object
tested.  U is a property of the LSJJ.DAT deck and not of the run; its J
blocks are at worst a few tens square (327 LS states in all for f**7, spread
over 2J = 1..31); and LSJJU caches on (l,N), so CFCUOR keeps a SAVEd marker
and does real work once per subshell l**N rather than once per J block of
every configuration set.  The whole test is microseconds.

Two details make it the same test and not merely a similar one.  The fibre
rotation in LSJJRT sums over the rows IRU(1..NI1) and not over I0+1..I0+NI1;
IRU is a permutation of that range, anything else having already been
rejected with IERR=7, and a permutation of the rows leaves every column dot
product alone, so summing in the natural order tests the same thing.  And
the tolerance is the 1.D-7 of the check it replaces: that check saw U's
errors through an orthonormal T, which neither magnifies nor hides them, so
the two trip at the same place.

Where the test now sits is an improvement in its own right.  CFCUOR is
called from LSJJRT immediately after LSJJU, inside the MODE=0 pass, so a bad
U comes back as IERR=8 before any column has been rotated.  CALCFC's
existing message then reports it -- "LSJJRT FAILED, IERR= 8 ... MATRIX LEFT
IN J1J2 COUPLING; NO TRUE jj OUTPUT FOR THIS J" -- instead of the old
behaviour of printing a warning about a matrix that had already been
spoiled.  The /CFCINT/ common block and the column-range machinery at the
tail of CFCORT existed only to feed the removed check, and are gone with it.

Measured on Pr III with timers around the two regions of CALCFC, summed over
all 50 blocks of the run: de-duplication 0.129 s -> 0.000 s, orthonormality
0.404 s -> 0.001 s.  Those are small numbers because Pr III is a small case;
the point is the exponent, and Pr III is what can be checked byte for byte.
The Pr III run with KCPLD(8)=8 gives TAPE2E and OUTGINE identical to the run
before these changes, with only timing lines differing in OUTG11.

CFCUOR was also checked to have real detection power, which a test that
never fires does not otherwise demonstrate.  A build with UMAT(1,1)
perturbed by 1.D-2 for l=3 prints

     *****THE LS -> jj MATRIX FOR l= 3** 3 IS NOT ORTHOGONAL.
     *****LSJJRT FAILED, IERR= 8  J=  0.5.  MATRIX LEFT IN J1J2 COUPLING;
     NO TRUE jj OUTPUT FOR THIS J.

once per J block, and the run still reaches a normal exit with its ordinary
output intact.


RCE: dimensions, and the arrays that follow them
------------------------------------------------

RCE20KD.F is built the same way, with -O3:

    gfortran -fshort-enums -malign-double -ftracer -fno-backslash
             -fallow-argument-mismatch -fcray-pointer -O3
             -fopenmp -fno-automatic -o RCE20KD.exe RCE20KD.F

IA, the largest J matrix, was raised from 10000 to 20000 to match KMX in
RCG: a case RCG can now diagonalize is a case RCE must be able to fit.  The
cost is about 9.2 GB of committed address space at start-up against 3.0 GB
before, and it is address space rather than memory -- the resident set stays
small until the blocks are actually that big.  The arrays that scale as
IA**2 are A and C, at 8*IA**2 each, and A1 at 8*IA*(IB+2).  BUFY, which
holds the eigenvectors of every block at once, grows as the sum of IM**2
over the blocks and is indexed by INTEGER*4, so that sum must stay under
2**31, i.e. 17.2 GB of BUFY.

A1 is where the doubling would have hurt, and does not.  D_PTR is a union of
DELX(IB), D(JB,1) and A1(IA,IA), and A1 was allocated at IA**2 -- 3.2 GB at
IA=20000.  But A1's second index is a radial-parameter index and not a level
index: in EIGEN it never exceeds 3+LMAX-I1 .LE. IB+2.  So IA*(IB+2) is
enough, and the allocation is 0.96 GB.  The zeroing loop and the I2 clamps
in EIGEN carry the matching bounds.

The large local arrays of the output routines are declared SAVE.  Under
-fopenmp gfortran would place them on the stack, which overflows at once --
60 MB for NJMS alone.

Two more cubic loops went the way CALCFC's did.  Reordering the eigenvectors
by calculated energy was a sequential extraction of minima that rescanned
the whole block for every level and, for every candidate, rescanned
everything picked so far: O(IM**3), and about three minutes per parity on
Gd III.  ESORTW does it as a heapsort, O(IM*log IM), and yields the same
permutation including ties -- equal energies come out in decreasing index
order, which is what the original produced, its test W(N1,M).LE.GRE having
kept the last of several equal energies.  And the MAXC largest amplitudes of
an eigenvector used to be picked by MAXC successive searches through all IM
amplitudes, IM*IM*MAXC operations per block; they are now collected in a
single pass with a bounded min-heap and heap-sorted into descending order,
IM*(IM+MAXC*log MAXC).  The ordering is the same total order as before --
decreasing absolute amplitude, ties to the lower basis index -- and zero
amplitudes are still never selected, so NW and IEVV are unchanged.

Finally, the timer.  ETIME returns CPU time summed over all threads, which
overstates elapsed time by the number of active threads once the OpenMP
regions are in use, so every printed time in a threaded run was wrong by a
factor of the thread count.  SECONDS now uses SYSTEM_CLOCK, guards against
counter wrap-around, and falls back to CPU time only if no clock is
available.


================================================================================
PER-CONFIGURATION TERM CUTOFF (columns 83-85 of a configuration card)
================================================================================

Restricting the number of parent LS terms of a large core shell is what made
lanthanide calculations feasible in the 1960s.  When the core has widely
separated terms, the levels built on the highest ones lie far above anything
observed and can be dropped.  RCG has always had a switch for this, NLT11 in
columns 26-27 of the main card, but it is one number for the whole run and it
exempts the last configuration of each parity, so it cannot serve two
configurations whose cores differ.  Dy VIII has f6, f7 and f8 cores; no single
NLT11 cuts all three sensibly, because the multiplicity-block boundaries of f6
and f8 are {1,17,73,119} and those of f7 are {1,7,47,119}, and they meet only
at 1 and 119.

Columns 83-85 of each configuration card now give that configuration its own
limit, applied with no exemption.  Blank means "use NLT11", which is what every
existing deck has, so nothing changes for them.  The card Cowan wrote ends at
column 80 and RCN2 still writes it that way, so old records simply blank-fill.
Note that RCN2 will overwrite a hand-punched field if ING11 is regenerated;
this is the same exposure that hand-edited Slater parameters already have.

WHY THE ALLOWED VALUES ARE A SHORT LIST

The cut keeps the first NLTV terms of the CFP deck, and deck order is not
energy order.  It is not even consistent in direction from shell to shell:

  f2    7 | 3:1-1  1:2-3  3:4-5  1:6-7
  f3   17 | 4:1-1  2:2-5  4:6-8  2:9-13  4:14-14  2:15-17
  f4   47 | 1:1-20  3:21-42  5:43-47
  f5   73 | 2:1-46  4:47-70  6:71-73
  f6  119 | 7:1-1  5:2-17  3:18-73  1:74-119
  f7  119 | 8:1-1  6:2-7  4:8-47  2:48-119
  f8  119 | 7:1-1  5:2-17  3:18-73  1:74-119
  f9   73 | 6:1-3  4:4-27  2:28-73
  f10  47 | 5:1-5  3:6-27  1:28-47
  f11  17 | 2:1-12  4:13-17
  f12   7 | 1:1-1  3:2-2  1:3-3  3:4-4  1:5-5  3:6-6  1:7-7

(multiplicity : first-last index, in deck order).  Within a block the order is
the CFP tabulation order -- not L, not seniority, not energy -- so a cut inside
a block keeps an arbitrary set.  A cut that ends ON a block boundary of a
descending deck keeps whole blocks of the highest multiplicities, which by
Hund's rule are the parents lying lowest, and that is the truncation the
physics justifies.  Hence the accepted values:

  f6, f8 : 1, 17, 73, 119        f9  : 3, 27, 73
  f7     : 1, 7, 47, 119         f10 : 5, 27, 47

For f4, f5 and f11 the deck runs the other way and for f2, f3 and f12 it
alternates, so no leading part of them is a set of low terms -- a prefix of f4
would discard the 5S..5I block, the ground terms.  Their only accepted value is
the full count, i.e. no cut.  d5 is the obvious trap outside the f shells: it
lists 6S last of its sixteen terms, so a prefix cut would throw the ground term
away.  That, and the fact that truncation buys nothing for a small shell, is
why only f shells are accepted at all.

Anything else stops the run, and the message names the configuration and lists
the values that shell allows.  The boundary table was read out of the CFP decks
themselves (cowan_cfp.py, Deck.terms, which is in deck order) rather than typed
from a reference, so it is what LOCDSK actually hands back.

TRUE jj

A cutoff and the true jj rotation cannot both be had.  LSJJRT rotates a whole
fibre of the LS basis at once and the deck in LSJJ.DAT carries the complete
term list of the shell; with terms missing, the term map is not a bijection and
the block shape does not match -- IERR=5 and IERR=6.  Rather than let the run
die inside PLEV, the cutoff wins: IJJROT is cleared in LNCUV and a warning goes
to OUTG11 and to stdout.  This also replaces the IERR=5/6 death that NLT11 plus
KCPLD(8)=8 used to produce.  Cowan's JJ (J1J2) output is unaffected either way.
The earlier "TRUE jj COUPLING WILL BE ADDED" line is still printed above it,
because it is emitted before the configuration cards are read; the later
message explicitly retracts it.

WHERE IT LIVES

  NLTCF(KC,2), ICUTNY   end of /C1/ in rcgpar.f
  NLTCHK                decode, validate, abort; the boundary table
  LNCUV                 reads columns 83-85, calls NLTCHK, trims the printed
                        term table, clears IJJROT
  PLEV                  the cut that matters, before the READ(ID2) that fills
                        the term arrays

MEASURED

Dy VIII (f7p4, f8p3, f6p5 / f7p3d, f6p4d, f7p3s, f6p4s), 8 threads:

  no cutoff        CALCFC alone about 20 hr before the CALCFC work,
                   about 2 hr 7 min after it
  cutoff 1         whole run 1.9 s
  17 / 7 / 17      whole run 3 min 35 s, CALCFC 11.7 s

Pr III with the field blank reproduces TAPE2E and OUTGINE byte for byte and
OUTG11 apart from the timing lines.


WALL-CLOCK TIMING IN RCG

SUBROUTINE SECONDS called ETIME, which returns CPU time summed over all
threads.  Serially that equals elapsed time, which is why it went unnoticed
for so long; under OpenMP it is elapsed time multiplied by roughly the number
of active threads.  ENERGY and the line sort are the first heavily threaded
stages, so every TIME= and ABS TIME= printed from there to the end of the run
was inflated, and the inflation grew as more of the run was spent in parallel
regions.  The earlier stages were unaffected, which made the listing internally
inconsistent as well as wrong.

SECONDS now calls SYSTEM_CLOCK and subtracts the count read at the first call.
That first call is CALL CLOCK(TIME0) at the top of the main program, so ABS
TIME keeps its old meaning of time since start-up.  INTEGER*8 counters, with a
wrap-around guard and a fall back to zero if the system reports no clock.

This is the same fix already made in SUBROUTINE CLOCK of RCE20KD.F.  In RCG the
change is in SECONDS rather than CLOCK because CLOCK here is a two-line wrapper
that divides by 60; SECONDS is its only caller-of-substance and is called from
nowhere else.  No CALL CLOCK lies inside a parallel region, so the SAVEd
variables are not a thread-safety problem.

Timing is printed for information only, so no computed result moves.  Verified
on Pr III: TAPE2E and OUTGINE byte-identical to the references, and the
reported total agrees with real elapsed time (0.18 min reported against 12 s
real, the balance being start-up before the first CLOCK call and the final
writes after the last).  Running at 1 and at 8 threads now reports the same
elapsed time instead of times differing by the thread count.


RCE REFUSED TO FIT WHEN EVERY LEVEL WAS KNOWN

Symptom, found on the Li-like U case in U90: RCE printed "Iteration No 0
finished with AVDEV=6791.912" and stopped.  No cycle 1.  Two configuration sets,
both behaved the same way.  The listing gave the tell: CYC NO= 0 with DET B=
0.000E+00 and STDEV overflowing its field, i.e. the least-squares matrix was
never formed at all - the run was ending, not failing.

The cause is in the block of guards just before label 1114 that decide whether
to exit after the 0-th cycle.  One of them read

    IF (IFIXL.EQ.0) IFIRST=1

and IFIXL counts the UNKNOWN levels - those written with an asterisk in column
14 of RCEINP, for which no experimental energy is given.  So the test read "no
unknown levels, therefore do not fit", i.e. it abandoned the fit precisely when
every level had a measured energy.  That is the best-determined case there is,
and the ordinary one for a few-electron ion; in U90 all 10 levels of set 1 and
all 9 of set 2 are known, so the fit could never start.

The comment above the guard records the intent as "condition that there are no
fixed levels", which is the opposite condition - nothing to fit TO - and that
one already exists two lines further down as IFIXL.EQ.IMSUM.  It was, however,
dead: IFIXL was cleared inside the DO 105 loop over J, so what survived the loop
was the count for the last J block alone, while IMSUM is a sum over all J.  With
MMAX>1 the two could not be equal.  IFIXL is now cleared once, before the loop,
so it is the total the comment says it is, and the IFIXL.EQ.IMSUM guard works.

A second consequence is worth recording because it corrupts input rather than
just output.  IFIRST also drives label 1342,

    IF (NF(J,M).LT.0.OR.IFIRST.NE.0) ERAS(1)='*'

so when the run gives up after cycle 0 it stars EVERY level in RCEOUT.  RCEOUT
is the file one copies back over RCEINP to continue, and a starred level is an
unknown one - so the aborted run also threw away the experimental energies of
the case.  In U90 the supplied RCEOUT had all 19 levels starred against an
RCEINP with none.  With the guard removed, RCEOUT comes back with no asterisks
and the fitted energies in the second column, ready to feed back in.

The other half of the 1/3/2025 change - dropping Cowan's original
IF (IVARP.GT.IMSUM-IFIXL) test, i.e. refusing to fit more parameters than there
are known levels - was deliberate and is left alone.

Verification.  U90 now converges: set 1 AVDEV 6791.912 -> 191.931 over four
cycles, set 2 4550.311 -> 190.618, both stable to the third cycle.  For runs
that did fit before, nothing moves.  Two controls, each run against a build of
the unmodified source with the same flags:

  - U90 with one level per set marked unknown (the case that DID fit before):
    OUTE, RCEOUT, LEVELS1, LEVELS2, LEVELS3 and PARVALS all byte-identical.
  - U90 with every level marked unknown, which exercises the IFIXL.EQ.IMSUM
    path that the fix brings to life: same six files byte-identical, both
    builds stopping after cycle 0 as they should.


THE PARAMETER COUNT: FREE PLUS GROUPS, AND GAPS IN THE GROUP NUMBERING

Two separate counts are involved and they were in different states.

IVARP, the count used by the go/no-go guards that decide whether a fit is
possible at all, was already correct: the loop after DO 116 counts each freely
varied parameter once and each linked group once, using LPG to recognise a group
already seen.  Nothing needed doing there.

PNO, the count that carries the degrees of freedom, was not.  PNO comes from
JLMAX, and JLMAX was set in DO 116 by

    JLMAX=MIN0(JLMAX,LFL)

i.e. the most negative group number, which is the LARGEST group number in use.
That equals the NUMBER of groups only if they run 1,2,3,... with no gaps.  In
practice gaps are common - a group dropped while a fit is being developed leaves
a hole in the numbering, and there is no reason for the user to renumber.

What a gap does: the V* condensation in CALCV sums each parameter into column
J=-LF(L), so a group number that no parameter carries leaves that column of V
identically zero.  B(J,J) then comes out zero and is forced to 1 by the guard in
CALCV, with D(J,1)=0, so DELX(J)=0.  The fitted parameters are therefore NOT
affected - which is why this could sit unnoticed.  What is affected is PNO,
which is too large by the number of gaps, and PNO sets the denominator in

    STDEV=DSQRT(R/(TNO-PNO))

and STDEV multiplies every SDX in the main loop.  So all the parameter
uncertainties came out too large, by a factor that depends on nothing more
physical than how the groups happen to be numbered.

Measured on U90.  Renumbering the single group from 1 to 2 - same fit, same
AVDEV 191.931, same DET B 1.573E+51 - moved the printed count from 7 parameters
to 8 and STDEV from 350.416 to 429.170, +22 pct on every uncertainty.  A sharper
version of the same deck, with two groups numbered 2 and 5 instead of 1 and 2:
the old code counted 10 parameters against 9 levels, TNO.LE.PNO, so it took the
STDEV=0 branch and every SDX in the listing came out exactly 0.000 - the
uncertainties were not merely wrong but silently gone.

Fix: the groups are renumbered 1,2,3,... in order of first appearance, right
after DO 116 and before anything counts them, and JLMAX is set from that count.
Then JLMAX, PNO and IVARP all count things that exist, and PNO equals IVARP in
the linked case.  Two details worth recording:

  - LPG is now built from LF, as DO 116 left it, not from LF0, the raw input
    flag.  The two can disagree: the old loop tested IABS(LF0).NE.100 and would
    read, say, LF0=105 as group 55, while DO 116 treats anything with
    IABS(LF).GE.100 as fixed.  Building from LF removes the second reading.
  - LPG is cleared over its whole length rather than 1..LMAX.  A group number
    may exceed LMAX, and under -fno-automatic the array is static, so it would
    otherwise still hold the previous configuration set's counts.

INOCHG is deliberately left on its original test, which is on LF0, so the
"nothing changed in OUTGINE" detection is untouched.

Verification.  A deck whose groups are already consecutive is bit-exact: OUTE,
RCEOUT, LEVELS1, LEVELS2, LEVELS3 and PARVALS all identical to the previous
build (OUTE up to the timing lines) on the U90 case, on the same case with some
levels unknown, and on the same case with all levels unknown.  A deck with gaps
now produces exactly the output of the equivalent gap-free deck - the same six
files byte-for-byte, the only difference being the group numbers echoed back in
RCEOUT, which are the user's own and are deliberately preserved.


================================================================================
RCG: THE COUPLING-3-7 TRANSFORM IN RCEINP
================================================================================

The case.  Pr3/JK is the Pr III deck with NLSMAX=999 and KCPLD=011101 1, i.e.
KCPRCE=5, LSJK.  With that card the run took 9 min 04 s, and RCEINP took 8 min
53 s of it - everything else in RCG together took 11 s.  The J=3.5 block alone
(NLS=617, NPARS=525) took 2 min.

Where it went.  When exactly one of KCPLD(3..7) is zero and NLS.LE.NLSMAX,
RCEINP sets ITRANS=1 and transforms EVERY coefficient matrix of the block into
that coupling, C <- T C T', with two full NLS**3 products each.  NPARS is 525 to
645 here, so the block cost is 2*NPARS*NLS**3 - about 2.5e11 multiply-adds for
J=3.5.  The transform is the whole of RCEINP; the packing and the I/O are not
in it (TAPE2E is 38 MB).

What was done.  T is block diagonal.  CPL37 leaves TMX(L,LP) at zero unless
NCFG(L)=NCFG(LP) and the subshell quantum numbers the coupling does not touch
agree - the same structure CFCORT exploits for the orthonormality check - and
PLEV numbers the basis states one configuration at a time, so the support of a
row is an interval of columns.  New SUBROUTINE RCEBND measures that interval,
once per J block, by scanning TMX for the first and last non-zero of each row.
It measures rather than assumes: a dense T would come back as 1..NLS on every
row and RCETRB would then do exactly the work the old triple loops did.

New SUBROUTINE RCETRB does the two products with the inner sum restricted to
that interval.  Two further points:

  - The intermediate is held transposed, CT4(LP,L) for the old CT4(L,LP).  The
    old second product read CT4 along a row of a column-major KMX by KMX array,
    one cache line per term.  Both products now walk columns.
  - Both loops are threaded over the index that owns a whole column of the
    output, so no sum is ever split across threads.  The second loop is
    triangular, hence SCHEDULE(DYNAMIC,8).

Bit-exactness.  Every term dropped has TMX=0.D0 exactly, by construction of
JLO/JHI, and adding 0.D0 to a sum that started at +0.D0 changes nothing (the
sum can never be -0.D0: it starts at +0.D0 and x+y=-0.D0 needs both operands
-0.D0).  The surviving terms are still visited in increasing LPP.  The
transpose is pure addressing.

Also fixed here: a debugging WRITE(*,*) left live inside the NLS**2 loop that
builds CT4 on the KCPL=2 path.  It printed one line per matrix element, so a
KCPL=2 run with couplings 3-7 spent essentially all of RCEINP in list-directed
output.  Commented out, not deleted.

Result.  Pr3/JK: 9 min 04 s -> 26.7 s, RCEINP 8 min 53 s -> 15 s.  TAPE2E,
OUTGINE and genealogy.txt byte-identical to the previous build; OUTG11 differs
in 264 lines, all of them timing lines.  The Pr3_2026_v0 regression
(KCPLD=0011111, ITRANS=0, so RCETRB never runs) reproduces TAPE2E.ref and
OUTGINE.ref byte for byte, and its OUTG11 differs from the previous build in
126 lines, all timing.


RCE: WHERE THE TIME ACTUALLY GOES, AND SUBROUTINE CALCV
-------------------------------------------------------

The question asked was whether the coupling transformation to schemes 3-7 can
be accelerated in RCE the way it was in RCG.  It cannot, because RCE does not
perform one.  Couplings 3-7 are applied by RCG, in RCEINP, before the matrices
are written; what reaches RCE on TAPE2E is already in the chosen scheme, and
KCPL(1)/KCPL(2) merely record which two schemes those are.  The only coupling
transformation RCE does at all is the change from the first scheme to the
second, and it is a ready-made matrix A read from TAPE2E and applied to the
eigenvectors in ORDOUT2 - the EJK loop, which an earlier session had already
made sparse in C and threaded over K.

Measured on the Pr3/JK TAPE2E with the author's fitted RCEINP (a 5-iteration
fit of each parity), by bracketing the phases with CLOCK:

  TAPE2E read, unpack, scratch write     0.6 s
  CALCV                                 20.9 s      68%
  CALCX                                  0.0 s
  EIGEN (TRED2/TQL2 + ORDER)             7.9 s      26%
  the whole output block, both passes    1.4 s       4%
      of which ORDOUT2, i.e. the           0.10 s
      coupling transformation itself
  ---------------------------------------------
  total wall                            31.5 s

So the coupling transformation is 0.3% of the run.  The bottleneck is CALCV.

CALCV builds the V matrix: for every fitted level K of every J block and every
free parameter L, the quadratic form V = Y(.,K)' C Y(.,K), where C is that
parameter's coefficient matrix.  It was swept densely - IM**2 per K, hence
(number of fitted levels)*IM**2 per (block, parameter), every least-squares
iteration.

C is very sparse.  It is the coefficient matrix of a single parameter, and WR1
already stores it as a list of non-zero elements for exactly that reason; on
the Pr III case a 141x141 C carries between 2 and 48 non-zeros out of 19881.
Its non-zeros are now collected once per (block, parameter), outside the K
loop, into compressed rows - IRSTA/IROWL/IRCOL/RVAL, indexed by the ordinal of
a non-empty row so that all-zero rows cost nothing - and only those are swept.
The cost becomes IM**2 for the scan, the same order as the RD1 call above it
already costs, plus O(non-zeros) per K instead of IM**2 per K.  MAXNZ4 is a
ceiling, not a prediction: a C with more non-zeros than that falls back to the
original dense loop, which is kept in place.

Bit-exactness.  Within a row the surviving products are visited in increasing
J, exactly as before, so nothing is reassociated.  What is dropped is the
addition of C(I,J)=0 to SUMJ, and SUMJ starts at +0.D0 and can never become
-0.D0, which makes SUMJ+(+/-0.D0) exactly SUMJ.  The same argument licenses
skipping an all-zero row: it contributes SUMJ*Y(I,K) with SUMJ exactly +0.D0.
The diagonal's extra factor FF is still applied after the multiplication by Y
rather than folded into C, because (C*Y)*FF and (C*FF)*Y are not the same
double.  No parallel region was added.

Result.

  Pr3_2026_v0 (LS) fit    17.9 s -> 7.4 s
  Pr3/JK (LSJK) fit       31.5 s -> 10.0 s
  Pr3/JK, INV1=0           3.7 s -> 3.6 s   (CALCV runs once there)

  CALCV itself            20.9 s -> 0.60 s

In all three runs LEVELS1, LEVELS2, LEVELS3, PARVALS, RCEOUT and RCEINP are
byte-identical to the previous build, and OUTE differs in no line that does
not contain TIME.

What is left.  EIGEN is now 74% of the remaining 10 s.  Its TRED2, TQL2 and
TQLROT are already threaded, but under IF(L.GT.1000), and the Pr III blocks
are at most 617, so they run serial here by design.  Lowering that threshold
would be the next thing to measure; it was not touched, because the threshold
is a tuning decision of an earlier session and the gain on blocks this size is
not obviously positive.

================================================================================
THE PYTHON PRINTOUT AGAINST THE PERL ONE                       (2026-08-15)
================================================================================

CODE/printout.py replaces printout.pl.  Both were run on the three reference
cases and the outputs diffed line by line.  Every difference is listed here,
with what causes it.  Nothing else differs: the level counts, the energies, the
g factors, the percentages and the order of the rows are the same throughout.

W70 (507 rows in LS + JJ; the true jj block is new and has nothing to diff
against):  136 rows differ.

  135  an intermediate LS term's J moves from parentheses into angle brackets,
       2p3.(2P*(3/2)) -> 2p3.(2P*<3/2>).  The NIST ASD rule; parentheses stay
       where the intermediate label IS a J, as in 2s.2p.(2P*<1/2>).(1).3p.
       conv_cowan.pl could not tell the two cases apart and converted both to
       parentheses.

    1  the level at 132313.7944, J=11/2, which the Perl printed with no
       components at all.  Its J block holds one basis state, so RCG forms no
       transformation matrix and prints no JJ section; genealogy.txt prints the
       chain regardless and a one-dimensional block has amplitude 1 by
       construction, so the row is recovered rather than dropped.  It reads
       2s.2p.(2P*<3/2>).(2).4f | (2,7/2), which is its neighbours' label at
       132307, 132308 and 132320.  No Fortran change was needed.

       (Energies in all of these files are in the unit the deck asked for,
       which in every reference run is 1000 cm-1, so a wavelength in
       angstroms is 1e5 / delta E.  They are quoted here as printed.)

U90 (38 rows in LS + JJ):  12 rows differ.  Nine are the bracket change, one is
the 1x1 level, and two are the n bug:

    2  1s -> 1s.2s2.  The configuration is 1s1 2s2; conv_cowan.pl:529-555
       reattached n from a fixed-offset parse of IN36 and lost the 2s2
       entirely.  n now comes from the IN36 expansion through the ING11 slots
       (toolbox/ing11.py match_shells), which is a permutation of the slots and
       never a positional assumption.

U90/reordered:  the label 1s.(2S).2s2 becomes 1s.2s2.  Both subshells are
closed or hold one electron, so no intermediate term is written and the
summation order cannot be seen; printing 1s first and then an intermediate
implied that 1s had been summed first, which in that deck it was not.  Where an
intermediate IS written the summation order is kept and the run reports how
many labels that affected.

Pr3 (16878 rows in LS + JJ + LSJK; the LSJK block is new):  9558 rows of the
LS + JJ part differ, in five classes.  A row is counted in every class that
applies to it, so the classes sum to more than the total.

  5211  the bracket change, as above.

  3870  a redundant parent term dropped:  4f.5p6.5d.(2D).(1D*).6s becomes
        4f.5p6.5d.(1D*).6s.  5d holds one electron, so 2D is its only term and
        saying so adds nothing; (1D*) is the running term of 4f.5d and is what
        the reader cannot put back.  This is the Handbook's own notation --
        010.tex:370-374 writes 4f7(8S*) 5d(7D*) 6p, with 5d bare.

        The same rule takes 5p5's own term off, one hole having one term:
        4f3.(2P*).5p5.(2P*).(3D).5d becomes 4f3.(2P*).5p5.(3D).5d.

   532  a closed subshell no longer dropped:  4f -> 4f.5p6.6s2, 4f3 -> 4f3.5p6.
        One configuration of this run has an open 5p, so 5p6 has to stay in
        every label or the labels name different numbers of electrons -- the
        author's rule.  --drop-full-shells 5p takes it out of all of them and
        leaves 5p5 standing where the subshell is open.

   254  a term written twice:  5p6.5d2.(3F).(3F).6p becomes 5p6.5d2.(3F).6p.
        The subshell's own term and the running term after it are the same
        thing here, 5p6 being closed.

   176  two levels of nearly equal energy exchanged.  Both programs order a
        parity block by the RCG eigenvalue, which OUTG11 prints to three
        decimals, and these are the rows where two eigenvalues are equal at
        that precision -- 180.926 for the pair at Ec = 180.9259 and
        180.9262.  printout.pl then leaves them in the order Perl's hash
        iteration
        happened to produce (printout.pl:88-99, the $i counter over an unsorted
        keys %{...}), which is not reproducible and not worth reproducing;
        printout.py breaks the tie on the fitted energy.  Perl's whole LS
        sequence was checked against the eigenvalues and is monotone in them,
        so nothing else separates the two orders.

        The ordering key itself is the RCG eigenvalue, not the fitted energy
        (printout.py _order), on the author's instruction of 2026-08-15: RCE
        prints its own level list in the order of its Ecalc, and RCG's
        eigenvalue is that same quantity recomputed from the rounded parameters
        RCE punched back.  Fourteen pairs of this run genuinely differ between
        the two keys and print in descending fitted energy.

params.txt differs in four ways, all of them deliberate:

  * the configuration column names the whole configuration rather than the
    ING11 label -- 1s2.2s where the Perl wrote 2s -- and it is filled in for
    every parameter, not only for EAV;
  * a parameter is named with its subshell:  Eav(2s), F2(2p,2p), zeta(2p),
    RD0(1s1s,1s2s).  ALPHA, BETA, GAMMA and T need no heuristic; PFGD emits
    them in the same per-subshell loop as Fk(ii) (rcg11kd.f:4383-4400);
  * values carry the digits RCEOUT prints rather than the five decimals of the
    rewritten ING11 deck, and the CI rows are given the same three decimals as
    the rest instead of one;
  * the standard deviation, the least-squares flag and the group are taken from
    OUTE by position rather than by name.  RCE names a parameter by its slot,
    so ZETA 3 occurs once per configuration that has a slot 3, and a lookup by
    name silently merges them -- which is what put the wrong uncertainty on the
    second zeta(2p) of the odd parity set of U90.

Cross-folder merging, which was hand work, now runs from the command line:

    printout levels.txt --couplings LS,JJ,LSJK --from LSJK=JK

Pr3 and Pr3/JK share a basis, so their levels are paired on the LS eigenvector
each was identified with; the fitted energies agree to better than 0.01 of the
run's energy unit, that is 10 cm-1.
U90 and U90/reordered do not share a basis -- the subshells are summed in a
different order -- so they are paired by energy order within each (parity, J)
group, and the run reports the largest difference per group; the worst is
0.025 in that same unit.  Which mode applies is decided by comparing the two
decks' configuration cards, and --match forces it.


================================================================================
THE PYTHON CONV_OUT AGAINST THE PERL ONE                       (2026-08-15)
================================================================================

CODE/conv_out.py replaces conv_out.pl: the level list with its outer-shell and
decay-rate columns, the line list with its intensity models, the cascade
matrix, and the two files IDEN2 reads.  Three new modules carry it:

  toolbox/spectrum.py   the Lund layout of the OUTG11 spectrum, FORMAT 1057
                        and 1058 (rcg11kd.f:13428, :13437).  Those are the
                        only transition FORMATs that print a level's number
                        within its J block, which is what ties a line to an
                        eigenvector; an OUTG11 written in the wide layout is
                        refused by name rather than half-read, because the
                        alternative is conv_out.pl's FindLev, which matched a
                        line back on energy, configuration number and term
                        label and needed several fall-backs to do it.
  toolbox/radiate.py    S, gf and gA for E1 to M3, the constants of
                        conv_out.pl:1650-1710.  calc_TP_M1E2.py (phase 10)
                        needs the same ones.
  toolbox/ai.py         the two read_aa tables, tied to the levels by
                        position and checked on J.
  toolbox/iden.py       dlv.dat, numset.dat, ENLEV.DAT, TRANS.DAT.

Both programs were run on all three reference cases with RCE LTE 1e6 1000 and
the outputs compared row by row, the level lists on (coupling, Ec, J) and the
line lists on the two levels' energies.  Every row of every file is present in
both.  The differences:

LEVEL LISTS

  the brackets, the closed shell kept, the redundant parent dropped, the term
  written once, the 1s -> 1s.2s2 reattachment and the blank 1x1 level -- all of
  them the same classes as printout.py, listed above, in the same numbers
  (U90: 9 brackets, 2 reattached, 1 filled;  W70: 135 and 1;  Pr3: 5399
  brackets, 3894 parents, 535 closed shells, 256 terms written once).

  the n, l and occ columns now name the outermost OPEN subshell.
  conv_cowan.pl:3043 matched /(\d+)([a-z]+)(\d*)/ against the whole label
  without /g, so it took the FIRST subshell: the ground level of W70,
  2s2.2p, came out as n=2, l=0, occ=2 -- 2s2, not 2p.  It is the outer open
  subshell that groups a level in a published table, so that is what is
  written, and closed shells are skipped whether or not they are printed.
  U90 23 rows, W70 45, Pr3 3135.

  212 rows of Pr3 carry a different Lev_#.  conv_out.pl numbers the levels in
  the order of RCE's fitted energy; the whole package now numbers them in the
  order of the RCG eigenvalue, which is the author's instruction of
  2026-08-15 and what printout.pl already did.  The two orders differ only
  where a pair is closer than the rounding between them.

  29 rows of Pr3 lose a doubled parity mark, 2H2** -> 2H2*, and 10 keep a
  Nielson-Koster index the Perl dropped, 2D -> 2D2.  conv_out.pl fixes the
  parity of a designation in two places with two different substitutions --
  conv_out.pl:271-274 in the level list and :994-999 in get_leading_term --
  so the same level could come out spelled two ways in the two files.  There
  is one designation here and it is built once.

  the header names %1 conf1 term1 through %5, which is what the rows hold.
  conv_out.pl's header stopped at %2 while its rows ran to five.

LINE LISTS

  Pr3 17594 lines, W70 707, U90 18 -- the same lines, in the same order, with
  the same gA, cf, lifetimes, branching fractions and wavelengths.  What
  differs is the modelled intensity, by at most 1 percent, and through it the
  handful of lines that sit on the 0.5 threshold at which a line is printed
  (one of W70's 708).

  Both programs model the intensity from RCG's gA at RCG's own wavelength.
  conv_out.py first used the fitted energies for the wavelength as well;
  the author stopped that on 2026-08-15, and the reasoning is worth keeping
  because it is easy to think the fitted wavelength is simply better.  It is
  a rescaling of A, and in a real case the observed levels are a small
  fraction of all, so it would be applied to some lines and not others, with
  no good answer for a line that has one known level or none.  Cowan's codes
  give unscaled A-values and that is what people expect; the error a
  wavelength leaves in a relative intensity is far smaller than the error in
  Cowan's line strengths, or in a measured intensity.  Where scaling is
  wanted it is a case-by-case job at publication time, and an easy one.

  What is left of the difference is that conv_out.pl reads RCG's wavelength
  from OUTG11's printed F11.3 field, which for the 0.12 A lines of U90 holds
  three significant figures, while conv_out.py computes it from the same
  eigenvalues the level list carries.  U90's largest disagreement is then
  0.5 percent, W70's 0.4 percent -- the normalisation to the maximum spreads
  one line's rounding over all of them -- and one of W70's 706 lines falls
  the other side of the 0.5 printing threshold.

  The cascade matrix is reproduced byte for byte: nothing in it depends on a
  designation.

WHAT IS NOT CHECKED AGAINST ANYTHING

  The DR and AUGER models and the IDEN2 files need data this repository does
  not hold -- a read_aa table and an observed line list.  They are written to
  the formats conv_out.pl writes and are exercised in the tests on a dlv.dat
  built out of two lines the U90 run itself predicts, which come back
  assigned with no residual.  They have not been run against real
  measurements.

  One thing that showed up there: U90's wavenumbers are 8e8 cm-1 and overrun
  the F12.3 fields of ENLEV.DAT and TRANS.DAT, which conv_out.pl fills with
  asterisks.  On the author's instruction (2026-08-15) these fields are now
  written by fortfmt.write_f_wide, which follows the rule used elsewhere in
  the Fortran: try F12.3, then F12.2, F12.1, F12.0, and only when the
  integer part alone will not fit fall back to exponential with the widest
  mantissa the 12 columns allow.  The record keeps its width, every field
  keeps its columns, and IDEN2 -- which is Java and parses the field as a
  number -- gets a number in every case.  Whether it can then work with
  numbers that large is its own business, not conv_out's.

  M1 and E2 in one run: where a pair of levels is printed under both
  multipoles, the stronger gA is kept and the weaker dropped, as
  conv_out.pl:554-566 intends.  conv_out.pl then added both to the upper
  level's total radiative rate; only the surviving one is added here.  No
  reference run has two multipoles, so this is reasoning, not measurement.

================================================================================
THE PYTHON REORDER_ING11 AGAINST THE PERL ONE                  (2026-08-17)
================================================================================

CODE/reorder_ing11.py replaces reorder_ing11.pl: it writes a copy of an ING11
whose subshell slots are in another order.  The work is in toolbox/reorder.py,
and CODE/reorder_ing11.py is the command line over it.

Permuting the configuration cards is the easy half -- eight (A1,I2,2X) fields
in columns 1-40, columns 41 onward untouched.  The hard half is that a
parameter card carries no names at all.  It is a row of (F9.4,I1) fields whose
POSITION says what each one is, because RCG emits them in an order that
depends on the slot numbers.  Move the slots and that order changes, and the
values have to move with it.

The Perl derived that order empirically, by comparing the parameter names and
their order in Cowan's output; its CI_class routine is a heuristic and the
plan recorded it as having gaps.  It has none that could be found.  It is
right on all three reference cases, and the two rules below -- taken from the
Fortran, not from the output -- reproduce it exactly.

WHERE THE ORDER COMES FROM

  PFGD writes, per configuration (rcg11kd.f:4300-5000):

    EAV
    for each slot i     Fk(ii) for k ascending, then ALPHA, BETA, GAMMA and
                        the T family, which PFGD emits inside the same
                        DO 199 I=1,NOSUBC loop and which therefore belong to
                        the slot the last Fk(ii) named
    for each slot i     ZETA i                       (DO 299 I=1,NOSUBC)
    for i < j, k asc.   Fk(ij)   (DO 399 I=1,NOSCM1; DO 389 J=I+1,NOSUBC)
    for i < j, k asc.   Gk(ij)   (DO 497 I=1,NOSCM1; DO 489 J=I+1,NOSUBC)

  and PRK writes the CI integrals of a pair of configurations in the order of
  DO IRHO=1,NOSUBC; DO ISIG=IRHO,NOSUBC; DO IRHOP=1,NOSUBC; DO
  ISIGP=IRHOP,NOSUBC (:5333-5339), each group being its D integrals by
  ascending rank and then its E integrals (:5781-5790).

  So a sort key built from the NEW slot numbers gives the order RCG will use
  for the reordered deck.  The sort is stable, which is what carries the rest:
  the ranks of one integral, and ALPHA/BETA/GAMMA/T behind their Fk(ii), share
  a key and keep the order they already had.  That also means the ALPHA family
  needs no rule of its own beyond the anchor -- the same finding as parnames.py
  (phase 5), used the other way round.

THE D/E LETTER IS POSITIONAL

  120D1323 is (configurations 1,2) rank 0, direct, slots 1,3 of the first
  configuration and 2,3 of the second (:5889).  D means the pairing (n1 with
  n3, n2 with n4) and E the pairing (n1 with n4, n2 with n3).  RCG always
  writes n1 <= n2 and n3 <= n4, so after a permutation one or both pairs may
  have to be sorted back -- and sorting exactly one of them exchanges the
  pairings, which turns a D into an E and back.  Where n1 = n2 or n3 = n4 the
  two pairings are the same integral and RCG calls it D; that case has to be
  tested first, because the sort test would otherwise flip a letter that
  cannot be flipped.  U90's 120D1112 is exactly that: reordered it is
  120D3313, and the naive rule makes it an E.

  None of this needs n.  The permutation is on slot POSITIONS and so is the
  pairing it has to preserve, so reorder_ing11.py reads no IN36 and resolves
  no subshell identities.  The Perl's refusal to move the last occupied
  subshell when it holds different nl in different configurations -- W70's
  slot 4 is 3p, 4p, 4f, 3d, 4s or 4d depending on the configuration -- existed
  only because the Perl did have to re-attach n.  It is gone.

THE NAMES COME FROM OUTG11

  outg11.read_parameter_names() is new: it reads the PARAMETER VALUES IN
  sections (FORMAT 13 at :10876, FORMAT 17 at :10928 for a CI card, then 14
  and 15 for the names and the values).  Three things it has to get right:

  - the names are taken by splitting on runs of three or more blanks, not on
    the 15-column grid, because past 22 names FORMAT 14 runs out of line
    groups and reverts, and the grid shifts by two columns (Pr3's f5d6d, 31
    parameters);
  - a field is a value if it holds a decimal point, not if float() accepts it.
    A CI name of exchange type reads as an exponential: float of 672E2424 is
    infinity, not an error.  W70's sp3p-sp4p has one;
  - the key is columns 1-18 of the card.  Both headings echo what RCG READ
    from the card with 2A6,A7, and RCN2 punches the parameter count in columns
    19-20, so a card with ten or more parameters has the first digit of that
    count inside its A7.  Pr3 has one: f5d6d    -f5d7d   1.

VALUES ARE MOVED AS TEXT

  A value is carried across as the ten characters RCN2 punched, never
  re-encoded.  RCN2 fills a field whichever way it fits -- 822465.00 with its
  point, 140570 without one, 0.00e+01 in exponential form -- and re-encoding
  is both lossy and pointless for a field that is only changing places.

  It is also where the one difference from the Perl is.  reorder_ing11.pl
  re-encodes every value at four decimals, and a value that then needs ten
  characters overflows its nine-column field and pushes the rest of the line
  one column right.  U90 has one, on the 2s-2p2 card:

    ours   2s       -2p2      1-1369.5055   0.00005 ...
    Perl's 2s       -2p2      1-1369.50505  0.00005 ...

  RCG reading 5(F9.4,I1) takes the value correctly and then reads 0 as the
  type code and 5.0000 as the next value, so the Perl's card is corrupt.  It
  happens for any value below -1000 on a parameter or CI card.  Everything
  else on that card is untouched by the reordering, so the deck the author
  keeps in WORK/U90/reordered does not have it.

WHAT WAS COMPARED

  Both programs on WORK/W70 (231, 312, 321, 213, 132), WORK/Pr3 (2134567,
  7654321, 3142567, 2345671, 1325476) and WORK/U90 (231, 312, 321, 213, 132).
  Byte identical everywhere except the one U90 line above, which the Perl
  corrupts.

  The stronger check is against RCG rather than against the Perl.
  WORK/U90/reordered was reordered and then RUN, so its OUTG11 says what RCG
  itself calls every parameter of the reordered deck and in what order it
  emits them.  reorder.rename() predicts all 9 sections, every name and every
  position, with nothing left over -- including the three CI cards where the
  letters flip and the one where the repeated slot holds a D.  That is
  test_reorder.py's acceptance criterion.

  reorder_ing11.py on WORK/U90 with 231 reproduces the author's
  WORK/U90/reordered/ing11 field for field.  Three values differ in their last
  digits, on lines 9, 10 and 12, and that is not the reordering: the reference
  deck was punched from the Hartree-Fock values scaled by 0.85 and ours from
  the parent's fitted parameters.  The line lengths, the columns and the order
  agree exactly.

A CORRECTION TO THE PLAN

  The plan recorded, under the phase 5 findings, that reorder_ing11.pl does
  not permute the slot digits of the CI cards -- WORK/U90/reordered appearing
  to ask for R1(1s2p,1s2p) where the parent has R1(1s2p,2s2p) -- and set it
  down to be fixed in phase 8.  That was misattributed.  The Perl's deck is
  correct; every one of its CI integrals was checked by hand against the names
  RCG generated for the reordered run.  The wrong names came from our own
  ing11.match_shells, which handed the slots their n values left to right and
  so called the reordered 2s2p card's first slot 1s when the deck means 2s.
  ing11.slot_order() fixes it, by taking the pairwise constraints that the
  configurations with unequal occupations impose and topologically sorting
  them; see the docstring there.  The reordered case now names the same
  integrals as its parent, and its 2s2p levels read 2s.2p.(3P*).1s rather than
  1s.2p.(3P*).2s.

A BUG IN ING11.PY FOUND ON THE WAY

  read_ing11 built deck.parameters and deck.ci by appending as it met each
  parity set, so a run with an empty set 1 put set 2's block at index 0 while
  every caller indexes them by kset-1.  No reference case has one parity, which
  is why it had not shown.  Both lists are laid out first now and filled in
  place.


================================================================================
THE SEVEN SIMPLE REWRITERS AGAINST THE PERL ONES               (2026-08-17)
================================================================================

Phase 9 of the plan: s11, add_eav, scale_param, subst_exp_levs,
subst_params_ing11, conf_list and read_aa.  Four of them edit a deck in place,
two read one and print, and one is a level-by-level merge of two RCE decks.

Three new modules carry the shared work:

  toolbox/conflist.py   the printed order of a configuration list, which is a
                        presentation choice of the author's and is kept exactly
  toolbox/explevs.py    matching the levels of two RCE decks on their
                        eigenvectors, which is what subst_exp_levs is

and, inside the modules that already owned the files,

  ing11.format_card / write_ing11 / set_values / take_card
  rce.edit_parameter_line / level_lines / set_observed
  outg11.read_auger

WHAT AGREES WITH THE PERL

  conf_list             byte identical on W70, Pr3 and U90
  add_eav               byte identical on all five decks (W70, U90,
                        U90/reordered, Pr3, Pr3/JK)
  subst_params_ing11    byte identical on five runs: U90 sets 1 and 2, Pr3
                        sets 1 and 2, W70 set 1
  s11                   byte identical on W70, Pr3 and Pr3/JK; both U90 decks
                        differ on one line, which is the defect below
  subst_exp_levs        every level line and every parameter value agrees on
                        U90 and Pr3 -- 593 observed energies and 1207
                        parameters on Pr III -- and only the trailing blanks
                        of a CI parameter line differ.  W70 is where the Perl
                        destroys the deck; see below.
  scale_param           nothing to agree with: the Perl matches no line
  read_aa               nothing to agree with: the Perl reads the wrong
                        parity set

The whole deck-rewriting family rests on one property, which is what makes
"byte identical" reachable at all:

  A CARD NOBODY EDITED COMES BACK AS THE BYTES IT ARRIVED IN, AND SO DOES A
  FIELD WHOSE VALUE HAS NOT MOVED.

read_ing11 keeps each value's ten columns beside the value it read from them
(ParameterCard.fields), and format_card puts the original characters back
wherever the value and the group code are unchanged.  Both halves are needed.
Without the first, a rewriter reformats every card it walks past; without the
second, it reformats every field of every card it does touch.  RCN2 and the
Perl have each punched a zero parameter in a form the other does not produce
-- RCN2's 1PE9.2 gives " 0.00E+00" and s11.pl's %10.2e-with-a-zero-deleted
gives "  0.00e+0" -- and RCG reads both as zero, so there is no way to choose
between them for a value nobody has changed.  Not choosing is the answer.

Tested on all five decks two ways: read and write with nothing marked, and
read and write with every card marked for re-encoding.  Both are byte
identical to the deck they started from.  Pr III's deck has 181 fields in the
Perl's exponential form and they survive untouched.

THE DEFECT THE AUTHOR REPORTED  (s11.pl:228)

  s11.pl punches a CI parameter with sprintf("%9.4f"), which needs ten
  characters for any value at or past 1000 in magnitude.  The field is nine
  columns, so the tenth character pushes the group code and all four
  remaining fields one column right and the card is corrupt.  WORK/U90 has
  one such parameter, 131D1233 = -1369.505304, and the card the Perl writes is

    2s       -2p2      1-1369.50535   0.00005   0.00005   0.00005   0.00005hr

  which RCG reads as the value -1369.5053, then 5 as the group code, then
  5.0000 as the next parameter.  The deck the author keeps in WORK/U90 has the
  card put right by hand.  s11.py writes the field through
  fortfmt.write_f_wide, which gives up a decimal rather than the number:

    2s       -2p2      1-1369.5055   0.00005   0.00005   0.00005   0.00005hr

  -- which is the author's own card, field for field.  The same holds for
  WORK/U90/reordered.  It is the only line of any reference deck on which
  s11.py and s11.pl disagree, and our whole U90 output is byte identical to
  the author's deck but for the one zero parameter noted above.

  The Perl's card is not merely different: it cannot be read back.
  ing11.read_ing11 refuses it, because the "hr" tag it looks for in columns
  71-72 has been pushed to 72-73 and the scale factors behind it no longer
  read as integers.

A SECOND DEFECT OF THE SAME KIND  (subst_exp_levs.pl:73, :364)

  The Perl recognises a level line with a regex that wants a blank between J
  and the leading eigenvector percentage.  RCE writes J as I4 then A2 and the
  percentage as I4 behind it, so a half-integer J fills its six columns and a
  percentage of -100 fills its four, and there is no blank:

    204712.520000*204712.524   1/2-100 ss4p   (1S) 2P        1 ss4p   (3S) 4P

  An integer J leaves the A2 blank and the same line matches, which is why
  this shows on some runs and not others.  Perl keeps the previous match's
  captures when a match fails, so from the first such line the level counter
  and the J are stale; every lookup behind it misses, and PrintOut writes an
  undefined energy -- 0.000000 -- over each remaining level.

  WORK/W70 has eleven such lines and the Perl destroys fifty observed
  energies below them.  The right answer for W70 is to change nothing at all:
  W70 was never fitted, so its RCEOUT carries no observed level and there is
  nothing to substitute.  That is what subst_exp_levs.py does.

  Nothing here reads a level line by regex.  rce.level_lines walks the file
  the way read_rceout does and rce.set_observed puts the energy back into the
  columns it came out of, leaving every other column of the line alone.

TWO PROGRAMS THAT NO LONGER DID ANYTHING

  scale_param.pl:58 recognises a parameter line with a regex that requires the
  maximum-value field to read 0.0000.  RCE writes that field as F14.6 --
  0.000000 -- and has for as long as any deck here goes back, so the regex
  matches nothing and the Perl copies its input unchanged.  Verified on all
  three RCEINPs, for every combination of F=, G=, CI= and GROUP=.
  scale_param.py takes the fields from the FORMAT that writes them
  (FOR/RCE20KD.F:1470, and :1464 for the variant with a denominator) and so
  has nothing to fall out of step with.

  read_aa.pl stops at the first J that does not increase (:104), which is the
  start of the second parity set, and so can only ever reach the first.  In
  WORK/AI/DIEL -- the only run in the package that has autoionization output
  at all -- every AA section is in the second set.  The Perl therefore reads a
  set that has none, matches stray numbers off the G-VALUES and PURITY lines
  as eigenvalues, and prints two rows of zeros.  read_aa.py takes a --set,
  the first by default, and refuses a set that has no rates by name, saying
  which set does.

  Its output on that case is right against the file: the J=6 block's
  18.58k (2S) 3K channel reads 0, 0, 5.51e+01, 1.02e+05, 4.89e+04, 2.03e-01,
  which is FORMAT 38's row at rcg11kd.f:11373, and the J=7 rows of the 3K and
  1K terms sum into one 18.58k column as the core pattern asks.

THREE MORE THINGS THE PERL GETS WRONG, ALL IN subst_exp_levs.pl

  * params=noflags substitutes nothing.  The value is held by the test at
    :448, which asks whether the parameter's NAME contains the word "noflags"
    and is therefore always true, and the flag is held at :451.  Here it does
    what its name says and what its pairing with "flags" means: the values,
    and not the flags.

  * a parameter's name is not unique in a file.  F2(11) belongs to every
    configuration with an f shell and 120D1112 to every pair it couples.  The
    Perl scopes a single-configuration parameter by the configuration the
    preceding EAV line named and a CI parameter by its pair and its position
    within it, and so does walk_parameters here.  Keying by name alone -- the
    first thing tried -- put one configuration's parameters on another and
    showed up at once on Pr III, where 292 values came out wrong.

  * Eav and CI name a KIND of parameter and F, G and Z name a letter its name
    begins with, which is why the Perl keeps three branches (:402, :418, :443).
    Matching "CI" against a CI parameter's name would never succeed -- they
    are called 120D1112 -- and matching "F" against a name would catch
    EAV 4f3 along with F2(11).  On Pr III the five options move 17, 66, 86,
    50 and 56 values and none of them the same value twice.

WHAT WAS DROPPED, AND WHY

  s11.pl reads an "Eav shift" off the deck's prelude (:21-23) and subtracts it
  from every fitted Eav.  The regex captures ONE character -- the repetition
  is outside the group -- so the shift is a single character of whatever card
  it lands on.  On Pr III it lands on a blank, which is true in Perl and zero
  in arithmetic; on the other decks nothing matches at all.  It has therefore
  never done anything on any deck here, and it is not ported.  Shifting every
  Eav together is what add_eav is for, and it does it properly.

  read_aa.pl reads RCG's KE FOR CORE lines into a hash (:136-158) and then
  looks them up under a key built a different way, so the lookup always
  misses; the only use of the result is commented out (:263).  The free
  electron's kinetic energy is the direct statement of whether a channel is
  open, and the test actually used -- the level lying above the core's energy
  -- is the indirect one.  The two need not agree in every case.  What the
  Perl does is what is done, and the KE lines are not read; which of the two
  is wanted is a question for the author.

  read_aa.pl's block that reads the configuration table (:61-96) is
  unreachable: $confs_read is set to 1 before the loop and never cleared.

  conf_list.pl:94-105 drops the first subshell of every configuration when all
  of them begin with f14.  That rule is kept exactly.  It is tempting to
  generalise it to the whole shared closed core, and it must not be: the plan
  forbids labels of unequal electron count except at the user's explicit
  request, and this rule is already the one exception the author made.

ONE THING NOT VERIFIED

  subst_exp_levs also reads a LOPT level list -- ten tab-separated columns, or
  a file named .LEV -- which is how a freshly optimized set of observed
  energies is brought into a fit.  There is no LOPT list among the reference
  cases, so toolbox/explevs.read_lopt follows subst_exp_levs.pl line for line
  and HAS NOT BEEN RUN AGAINST REAL DATA.  Its two pieces of real logic are
  unit-tested against the Perl's arithmetic -- the term-to-J expansion at
  :206-221, and the parity taken from the term's own asterisk with the first
  level read fixing which set is which -- but the column layout is taken on
  trust.  A sample file is the one thing this phase needs and does not have.

WHAT THE PERL DID THAT NOTHING HERE DOES

  Four of these programs move their output into place by shelling out to
  cmd.exe's copy, del and ren.  add_eav.pl also opens its deck as "ing11.",
  with a trailing dot that Windows drops and a POSIX filesystem does not, so
  under a Cygwin or MSYS Perl it cannot open the file at all.  Everything here
  goes through files.write_deck: the new deck is written beside the old one
  and moved into place, so a run that fails partway leaves the deck it was
  given intact, and the whole family runs on Linux and macOS unchanged.

  files.read_lines grew a strip= option for this.  A reader wants the trailing
  blanks gone so it can index a column; a rewriter needs them, because Cowan's
  sentinel cards at the end of a deck are padded out to column 80 and a card
  nobody touched has to come back as it was.

TESTS

  tests/test_rewriters.py   42, the four deck editors and the round-trip
                            property they rest on
  tests/test_explevs.py     24, subst_exp_levs and the eigenvector matcher
  tests/test_auger.py       12, read_aa and outg11.read_auger
  tests/test_conflist.py    11, conf_list

  309 pass in all, up from 220.

  A note for whoever runs them next.  The Perl comparisons need a NATIVE
  Windows Perl.  Run under the Cygwin Perl that a pytest started from a Git
  Bash prompt finds first, add_eav.pl cannot open its deck (the trailing dot),
  s11.pl's ren fails silently, and conv_cowan.pl fails on every OUTG11,
  because Cygwin's Perl does not strip the CR of a CRLF file in text mode and
  its " *$" anchors then never match.  The tests skip rather than fail where
  the Perl will not run, and the five test_conv_out.py tests that appear to
  fail under Cygwin pass under Strawberry Perl.  Nothing in phase 9 changed
  them.


================================================================================
calc_TP_M1E2 IN PYTHON, AND A PARAMETER FILE FOR IT           (2026-08-18)
================================================================================

Phase 10 of the plan, and the last of the Perl programs.  calc_TP_M1E2 puts an
uncertainty on a forbidden-line transition probability, which is not something
Cowan's codes compute: it displaces the fitted Slater parameters within the
standard deviations OUTE gives them and the E2 radial integrals within a
fractional guess, runs RCG again in a trial subdirectory, re-identifies the
eigenvectors with the fitted levels, and reports the spread of ln A over the
trials.

Two new modules carry the shared work:

  toolbox/parfile.py   a parameter file in the author's own layout
  toolbox/trials.py    the deviates, the running sums, and the two Excel
                       shape statistics the author checks the sampling with

and, in the modules that already owned the files,

  ing11.read_r_integrals / set_r_value / is_r_integral / Ing11.emax
  levels.identify(emax=, skipped=)
  case.Case.read(use_emax=)

--------------------------------------------------------------------------
THE THING THAT HAS TO BE SAID FIRST: THE PERL NO LONGER READS ITS OWN INPUT
--------------------------------------------------------------------------

calc_TP_M1E2.pl reads the transitions section of OUTG11 with its own reader
(ReadTransitions_M1E2, :338-461), not with conv_cowan.pl's.  That reader wants

  the header  "0           E      J   CONF              EP      JP  CONFP ..."

and rows of the shape

  ^.{5}(.{11})(.{5})(.{3})(.{9})(.{13})(.{5})(.{3})(.{9})(.{13})(.{12})...

which is FORMAT 56 and FORMAT 58 of SPECTR -- the WIDE layout, where both
levels of a line are on one row.  RCG only writes that layout when ILUND is 0,
and ILUND is set to 1 unconditionally at rcg11kd.f:12770.  Every OUTG11 this
package can produce is therefore in the Lund layout, and calc_TP_M1E2.pl
matches nothing in it: it runs to completion, prints "Done.", and writes an
output file with a header and no lines.  Verified on a purpose-built Sn VIII
E2 deck (WORK/Forbid) under Strawberry Perl -- one header line, nothing else.

conv_out.pl was updated for the Lund layout at some point and calc_TP_M1E2.pl
was not.  So:

  THERE IS NO PERL OUTPUT TO DIFF AGAINST FOR THIS PROGRAM, ON ANY DECK THE
  CURRENT FORTRAN CAN MAKE.  The port is verified against the Fortran FORMAT
  statements and against a real end-to-end run, not against the Perl.

This also corrects a statement in the plan (phase 7, note 1d): the wide layout
is NOT missing the level numbers.  FORMAT 53 and FORMAT 58 both end in 2I6,
NLEVJ and NLEVJP, which the author added, and calc_TP_M1E2.pl's last two
capture groups read them.  What is true is that the wide layout is unreachable
while ILUND=1, so spectrum.py reading only the Lund layout costs nothing.

--------------------------------------------------------------------------
WORK/Forbid -- A REFERENCE CASE THAT DID NOT EXIST
--------------------------------------------------------------------------

None of W70, Pr3 or U90 has a magnetic-dipole or electric-quadrupole spectrum,
and no ING11 in WORK/ carries a //R2// card, so there was nothing to run this
program on at all.  WORK/Forbid is the smallest deck that exercises it: Sn VIII
4d7, one configuration, one parity, IMAG and IQUAD set in columns 49 and 50 of
the G5INP card of IN2.  Inputs only, as WORK/AI is; its README says how to
build it and what it does not cover.

tests/test_calc_tp.py builds it -- RCN, RCN2, RCG, RCE into a temporary
directory -- and then runs the program over it three ways: with no trials, with
trials and a seed, and once to check what the trial decks contain.  That is a
real end-to-end test, with RCG actually re-run in the trial directories.

WHAT THE CASE DOES NOT COVER

  * The parameter variation itself.  Running RCE on a deck with no observed
    energies gives cycle 0 and a zero standard deviation for every parameter,
    so vary_parameters is exercised but has nothing to vary within.  Real
    observed energies would be needed, and inventing them is not an option.
  * The M1 + E2 merge.  No MAG DIP spectrum comes out of this deck even with
    IMAG set; RCG calls SPECTR for it (rcg11kd.f:11539) and prints nothing, and
    why is not known.  The merge is covered by unit tests on rows built by hand
    instead.

A deck of the author's that produces both multipoles -- the Fe V one the Perl's
header points at, C:\work\cowan\ti-like\fe5\Forbid -- would close both gaps.

--------------------------------------------------------------------------
THE PARAMETER FILE
--------------------------------------------------------------------------

The Perl held twenty settings as "my $..." at the top of the source
(calc_TP_M1E2.pl:20-40), including the name of the output file, so every ion
meant editing the program.  They are now calc_TP_M1E2.par, in the author's own
layout:

    <value>    ; <key>  <description of the parameter and its options>

The value is everything left of the first semicolon, the key is the first word
after it, and the rest is prose.  A line starting with a semicolon is a
comment, so a description can run on for as many lines as it needs; a blank
value means "use the default"; order does not matter.

WHY NOT configparser, JSON OR TOML.  The stock formats put the name first and
the description nowhere -- "num_trials = 100", with the explanation in a
comment nothing checks, or in a manual.  This layout puts the values in one
column where the eye finds them and makes the file its own documentation,
which is the same choice Cowan's fixed-form decks make.  The cost is
toolbox/parfile.py, which is 200 lines; the gain is that parfile.template()
writes a fresh, fully commented file from the same table parfile.read()
validates against, so the file and the program cannot drift apart.  There is a
test that asserts exactly that of the shipped CODE/calc_TP_M1E2.par.

What the module adds over a hand-rolled reader is the checking.  An unknown
key, a repeated key, a value that is not a number, a number outside its
documented range, and a value with no key at all are each named with the line
they are on and stop the run.  A parameter file is read once at the start of a
calculation that then takes hours; a typo must not be found by the results
being quietly wrong.

  calc_TP_M1E2 --write-par     puts a fresh commented copy in the working
                               directory (an existing one is never replaced)
  calc_TP_M1E2                 uses ./calc_TP_M1E2.par, or the copy beside
                               the script, and says which
  calc_TP_M1E2 --par <file>    uses a named one

The whole command line is now [folder] [--par <file>] [--write-par].  Nothing
about the physics is an option.

--------------------------------------------------------------------------
WHAT THE PYTHON DOES DIFFERENTLY, AND WHY
--------------------------------------------------------------------------

Every one of these is deliberate; the first four are defects of the Perl.

1. fixed_param_variance did nothing.  calc_TP_M1E2.pl:582 tests
   abs($flag == 100), which Perl reads as abs(($flag == 100)) -- 0 or 1 -- and
   a fixed parameter's flag is -100, so the branch never ran.  The correct
   idiom, abs($flag) != 100, is four lines above it in the same routine.  The
   Python tests abs(flag) == 100, and there is a test for it.

2. The trial total decay rate was taken from the wrong slots.
   calc_TP_M1E2.pl:838 falls back to averaging sum_A_hash[...][$i] over $i from
   0 to num_trials-1 when the unvaried run has no decay for a level: [0] is the
   unvaried run's own slot, which is the zero it is standing in for, and the
   last trial is never reached.  The Python averages trials 1..N.

3. Only the first uninterrupted run of //R2// cards was varied.  read_R2
   (:203-216) skips to the first //R2// line and then takes consecutive lines
   only; everything after the first non-R2 line is footer.  RCN2 does not write
   them all together -- see below -- so a deck can have R2 cards the Perl would
   never touch.  ing11.read_r_integrals finds them wherever they sit.

4. Only parity set 1 was ever done.  Every level read was labelled parity 1
   whatever it was.  parity_sets = 1 | 2 | both.

5. The varied parameters go straight into the deck.  The Perl wrote a stand-in
   RCEOUT into each trial directory and shelled out to s11.bat to transfer
   them, then to "del *.bak".  That stand-in RCEOUT has no configuration names
   in its CI parameter lines, which s11.py checks, so the route could not have
   been reused anyway.  The Python edits the deck in memory and writes it once.
   The correspondence it needs is checked instead, once, at the start: OUTE
   against RCEOUT name by name, and RCEOUT against the deck's parameter slots
   by count and by configuration name.  A mismatch stops the run and says which
   parameter.

6. The trial count may be any number.  The Perl drew normal deviates in
   Box-Muller pairs and sent the two members of a pair to two consecutive trial
   directories, so the count had to be even and read_only with it.
   random.gauss gives the same guarantee one at a time.

7. random_seed makes a run repeatable, which the Perl could not be.  A
   Monte-Carlo run that cannot be repeated cannot be checked, and these runs
   take hours.  Tested: the same seed gives the same file.

8. The merge is keyed on the level pair.  The Perl keyed its line hash on the
   wavelength rounded to three decimals, compared the level pairs anyway when
   two collided, and appended a '1' to the key when they differed.  merge()
   pairs the M1 row with the E2 row of the same two levels directly and sorts
   at print time.

9. A transition printed twice in one OUTG11 stops the run.  The Perl dies on
   this too (:717), but only after the whole file has been read; here it is
   caught in the pass that would have counted the trial twice.

WHAT IS KEPT THAT IS WORTH A LOOK

  Where BOTH levels of a line have observed energies, S is held fixed and gA is
  recomputed at the observed wavelength (calc_TP_M1E2.pl:485-487).  That is a
  rescaling, and conv_out.py deliberately does not do it (author, 2026-08-15).
  It is kept because it is what this program has always done and because a
  forbidden line is where it matters most -- gA goes as the cube or the fifth
  power of the wave number.  It is a question for the author.

  When the merge decides a line is E2 to within min_mixed_fraction, the row
  keeps the SUMMED gA of both multipoles under the name E2, while S stays the
  pure E2 value computed in the unvaried run.  The two then disagree by up to
  min_mixed_fraction.  This is what the Perl does -- the merge arithmetic runs
  before the cut-off test -- and it is left alone; it is the second question.

--------------------------------------------------------------------------
TWO THINGS THE FORTRAN DOES THAT THE READERS DID NOT KNOW
--------------------------------------------------------------------------

RCN2 WRITES A SAME-PARITY //Rk// BLOCK IN THE MIDDLE OF THE DECK.

  The R1 cards of W70, Pr3 and U90 all sit after the last parameter card,
  because they are cross-parity blocks.  An IQUAD deck also has a same-parity
  block, and RCN2 writes that one straight after the parameter card of the
  parity set it belongs to -- so on a two-parity quadrupole deck a //R2// card
  sits between parity set 1's parameters and parity set 2's.  RCG reads it
  there (SPECTR, rcg11kd.f:12824-12851); read_ing11 could not, and stopped with
  a float conversion error on "0.00(4E//R2/".

  _eat_r_integrals now takes such cards onto the parameter card they follow,
  and write_ing11 puts them back in the same place, so the deck still round
  trips byte for byte.  is_r_integral is RCG's own test: a bracket in column
  51, "//R" in 55-57, a digit in 58, "//" in 59-60 and a bracket in 64.

  set_r_value writes columns 39-50 through fortfmt.write_f_wide and puts
  columns 1-38 and 51-80 back as the bytes they came in as.  The Perl rebuilt
  the card field by field, which uppercased the 'hr' tag in columns 71-72 and
  stripped the trailing blanks of the 2I4 on every card it wrote.  There is a
  test that a trial deck differs from its original ONLY in the value fields of
  the //R2// cards.

E0MAX IS IN COLUMNS 66-85, NOT 61-80.

  The KCPL=0 rescale card is read at rcg11kd.f:673 as

    FORMAT (20X,F10.5,10I1,10X,F10.5,F5.2,2F10.5)

  so E0MAX(1) is columns 66-75 and E0MAX(2) columns 76-85; TEXCIT, F5.2, has
  columns 61-65.  calc_TP_M1E2.pl:194 reads the two from columns 61-70 and
  71-80.  The two readings agree for any card whose values are right-justified
  in ten-column fields, which is every deck in WORK/ -- the digits fall inside
  both windows and the rest is blanks -- but not for one punched to the left.
  Ing11.emax follows the FORMAT.

  E0MAX is used for one thing: levels.identify(emax=) leaves out an RCE level
  above it that has no matching RCG eigenvector, instead of stopping.  RCG
  prints nothing above E0MAX, so such a level is the deck saying so.  That is
  the condition at conv_cowan.pl:1299.  Every other caller passes no emax and
  is unchanged.

--------------------------------------------------------------------------
THE PERL COMPARISON TESTS RUN UNDER EITHER PERL NOW
--------------------------------------------------------------------------

The phase 9 note said the five test_conv_out.py comparisons fail under a Cygwin
or MSYS Perl, because that Perl does not strip the CR of a CRLF file in text
mode and conv_cowan.pl's fixed-column regexes then never match.  That is right,
and it is now fixed rather than described: conftest.perl_env() puts PERLIO=crlf
in the environment of every Perl subprocess, which asks for the layer by name.
The five tests pass under both Perls, and three that used to skip now run.

--------------------------------------------------------------------------
TESTS
--------------------------------------------------------------------------

  tests/test_parfile.py     23, the reader, what it refuses, the writer, and
                            the shipped CODE/calc_TP_M1E2.par
  tests/test_trials.py      21, the deviates, the accumulator, stdev and
                            kurtosis against hand-computed values
  tests/test_calc_tp.py     19, the merge, the parameter map, the variation,
                            and three whole runs over WORK/Forbid
  tests/test_ing11.py       +2, radial-integral cards wherever they sit

  369 pass in all, up from 309, with 5 skipped where a reference run is not
  unpacked.  Nothing that passed before changed.

==========================================================================
PHASE 12 -- PURITY PER CONFIGURATION, AND THE COUPLING CHOICE FILE
==========================================================================

The first of three tools that turn a finished fit into a published level
list.  It has no Perl ancestor: the author has been making this decision by
hand and by eye, and this only puts a number beside it.

  levels.txt --12--> coupling_choice.txt --13--> levels_prepub.txt --14--> levels_pub.txt

WHAT PURITY IS

  Purity of a configuration in a coupling scheme is the UNWEIGHTED mean, over
  the levels whose dominant configuration is that one, of the leading
  eigenvector component's percentage.  Unweighted: no 1/(2J+1).  A scheme
  that concentrates a level onto one basis state is a scheme that names it.

  Comparing that number between two schemes is legitimate because of an
  invariance.  A change of coupling recouples angular momenta WITHIN a
  configuration, so the recoupling matrix is block-diagonal in the
  configuration and the weight a configuration carries in a level is the same
  in every scheme.  Only the spread of that weight over the configuration's
  own basis states moves, and that is what purity measures.

  The same invariance is why a level has one dominant configuration whatever
  the scheme.  RCE computes it once, from full-precision coefficients, before
  any output pass (NJMS, FOR/RCE20KD.F:1157-1204).

--------------------------------------------------------------------------
TWO ADDITIONS TO printout.py THAT PHASE 12 NEEDED
--------------------------------------------------------------------------

Both are off by default, and with them off levels.txt is byte for byte what
it was.  Nothing but the phase 12-14 pipeline wants either.

  --dominant-config     adds a conf0 column after Lande_g: RCE's
                        scheme-invariant dominant configuration.

    It is read from LEVELS1, not from RCEOUT.  RCEOUT's configuration field
    is the LEADING COMPONENT's configuration, which moves with the scheme;
    LEVELS prints NJMS.  rce.Level now carries both under separate names, and
    case._refine(), which was already taking the energies and the g factors
    from LEVELS1, carries conf0 across with them.

    RCE writes that field as A6 (FOR/RCE20KD.F:1988, columns 44-49) while the
    label it names is the A8 of ING11, so a name longer than six characters
    arrives cut.  case._full_label() puts it back by matching against the
    section's own header names, which list the full names of exactly the
    configurations that can appear in that section (FOR/RCE20KD.F:1305).  Two
    names of one parity set agreeing in six characters stop the run by name.
    Nothing wider is searched and nothing is guessed: a wrong configuration
    here silently regroups a whole purity table.

  --from LSr=LS@reordered   names a block

    A levels.txt has one cpl value per block, so it could not hold two blocks
    of the same coupling -- which is exactly what U90 + U90/reordered is.  The
    options now build a list of Blocks (name, coupling, folder); the name goes
    in the cpl column and is the author's to choose.  The old --from JJ=folder
    form is unchanged and means Block('JJ', 'JJ', 'folder').  Asking for one
    coupling twice without naming the second is an error that says how.

--------------------------------------------------------------------------
toolbox/purity.py AND config_purity.py
--------------------------------------------------------------------------

  purity.py   the levels.txt reader (there was none: nothing had needed to
              read the file back), the levels-to-configuration map, the
              purity table, and the reader and writer of coupling_choice.txt.
              Phase 13 needs every one of them.

  config_purity.py   a thin front end.  --levels, -o, --drop-full-shells,
              --force, --merge.

  Grouping is on conf0 and never on conf1.  conf1 is the leading COMPONENT's
  configuration, so it moves with the scheme, and it carries intermediate
  terms besides: Pr III prints 4f.5p6.5d2.(3F) and 4f.5p6.5d2.(1S) as two
  strings for one configuration.

  The rows are the deck's configurations in IN36 order.  Their labels are
  rebuilt from IN36 alone, the way printout.py builds conf0 -- occupied
  subshells less the closed core the whole run shares, in (n,l) order -- and
  then every conf0 found in levels.txt is REQUIRED to be one of them.  The two
  constructions check each other; neither is trusted.  A conf0 that matches
  nothing stops the run and names it, because that means IN36 and levels.txt
  are not from the same run, or printout was given a different
  --drop-full-shells.

  Configurations that dominate no level are left out and reported.

THE FILE

  Fixed width, because it is the one file of the pipeline a person edits; the
  bulk files stay tab-delimited because they are read in a spreadsheet.  A
  provenance and definition block in '#' comments, then

    par  configuration  Nlev  <one 10-character field per scheme>

  in each field a star slot and the purity as %3d.  A star ANYWHERE inside a
  scheme's field selects it; exactly one star per line, and no star or two
  stars stop the run naming the line and the schemes.  The field boundaries
  come from the header line, so nothing counts spaces.

  Nlev is beside the purities on purpose: 88% over three levels is not
  evidence.

  A star is written already, on the highest purity, and the header says it is
  a suggestion.  The choice is not purely numerical -- previous literature on
  the spectrum, and the schemes used for the neighbouring isoelectronic ions
  in the database, routinely outweigh a few percent gained by a rarely used
  scheme.  The script prints which configurations are within 2% between two
  schemes, since those are the ones that want a decision.

  A marked file is never discarded silently: writing over one needs --force,
  and --merge keeps the stars where they are and writes the new purities under
  them, reporting anything that was marked and has left the deck.

WHAT IT SAYS ABOUT Pr III

  Two couplings, 39 configurations with levels.  The numbers come out where
  physics says they should: 4f2.5p6.6p is 49 in LS and 64 in JJ, 4f2.5p6.nf
  the same way, while 4f.5p6.5d2 is 46 in LS against 41 in JJ and every
  4f.5p6.5d.nl is LS.  4f3.5p5.5d, 2490 levels of a hole spectrum, is 11 in
  both -- which is the honest answer, not a defect.

TESTS

  tests/test_purity.py      27.  conf0 off by default; the two blocks of one
                            level agreeing on conf0, which is the invariance
                            itself checked against the file; the A6 cut put
                            back (p5f3d -> 4f3.5p5.5d); the missing-column
                            message naming the command to re-run; block
                            aliasing, and asking for one coupling twice;
                            purity against the mean recomputed from the file;
                            the star in each of the ten columns of its field;
                            no star and two stars; --merge keeping the mark
                            while the numbers change; --force and its refusal.

  396 pass in all, up from 369, with 5 skipped where a reference run is not
  unpacked.  Nothing that passed before changed.

==========================================================================
PHASE 13 -- ONE COUPLING SCHEME PER CONFIGURATION, APPLIED
==========================================================================

  levels.txt --12--> coupling_choice.txt --13--> levels_prepub.txt --14--> levels_pub.txt

The second of the three tools with no Perl ancestor.  Phase 12 measured the
purity of every configuration in every scheme and the author starred one
scheme per configuration; this applies the stars.  The whole tool is one
rule:

  for each level, take from each scheme's block only those components whose
  configuration was assigned that scheme; concatenate across schemes; sort
  descending by percentage; print.

So one composition may quote one configuration in LSJK and the next in LS
and the next in JJ, and the percentages still add up.  From WORK/Pr3, a
level whose dominant configuration is 4f.5p6.5d.7d:

    26%  4f.5p6.5d.(3G*<4>).7d          2[2]*      LSJK
    17%  4f.5p6.5d.(3F*<4>).7d          2[2]*      LSJK
    15%  4f.5p6.6s.(3F*).7s             4F*        LS
     6%  4f.(2F*<5/2>).5p6.6p2.(1D<2>)  (5/2,2)*   JJ
     6%  4f.5p6.5d.(1D*<2>).7d          2[2]*      LSJK

That is legitimate, not sloppy, and for the reason phase 12 rests on: a
change of coupling recouples angular momenta WITHIN a configuration, so the
recoupling matrix is block-diagonal in the configuration and the weight a
configuration carries in a level is the same in every scheme.  Only the
spread of that weight over the configuration's own basis states moves.
Sugar's Pr III tables are written this way.

--------------------------------------------------------------------------
THE INVARIANCE IS MEASURED, NOT ASSUMED
--------------------------------------------------------------------------

prepub.disagreements() sums each configuration's weight in each block of a
level and compares the blocks.  It cannot demand equality: a block loses
whatever fell below printout's cut or past its component count, and every
percentage is an integer rounded on its own, so a configuration spread over
ten components legitimately prints as 102%.  Two blocks may differ by

  sum over the two blocks of (what the block did not print
                              + half its count of components in this
                                configuration)

and no further.  On WORK/Pr3 with three schemes -- LS and JJ from the main
folder, LSJK from Pr3/JK -- 5626 levels and 39 configurations produce ZERO
exceedances.  That is the physics of phases 12-14 confirmed against a real
run rather than argued.

--------------------------------------------------------------------------
A DEFECT IN printout.py THAT THIS PHASE FOUND
--------------------------------------------------------------------------

printout.header() wrote five %n/confn/termn triples whatever --components
said, while _rows() wrote as many as were asked for.  With the
--components 10 the phase 12-14 pipeline is told to use, WORK/Pr3 came out
with a 24-field header over 39-field rows, and purity.read_levels_txt --
which finds the components by counting the header's triples, as it should --
silently dropped components 6 to 10.  That is exactly the material this
pipeline asks for.

The header now follows --components.  With the default of five it is the
file it has always been, and the byte-identity of the Perl comparisons is
unaffected.

--------------------------------------------------------------------------
WHAT THE TOOL DOES ABOUT WHAT IT IS NOT TOLD
--------------------------------------------------------------------------

  A component's configuration is its conf column stripped of the
  parenthesised intermediate terms (purity.strip_intermediates), matched
  against the labels of coupling_choice.txt.

  A configuration the choice file does not name.  config_purity.py writes a
  line only for a configuration that dominates at least one level, so a
  configuration that is never more than a small admixture legitimately has
  none.  Those take the FIRST scheme of the file and are reported by name
  and count.  Dropping them would take their weight out of the level and the
  percentages would stop adding up, which is worse than a scheme chosen by
  default.  A label in no scheme AND in no deck is the other case entirely
  -- levels.txt and IN36 are not from one run -- and that stops the run,
  naming the configuration and the line.  Where the folder has no IN36 the
  run continues with the weaker check and says so.

  A scheme starred in the choice file that levels.txt has no block for stops
  the run: ignoring the star would silently empty that configuration out of
  every composition.  A choice-file COLUMN with no block is only reported,
  since nothing is marked in it.

  Levels a block does not carry -- which the cross-folder merge can produce
  -- take their components from the other schemes, reported by count.

  Truncation is reported the honest way round.  Not "this level is under
  100%": a level of the 2490-level hole configuration 4f3.5p5.5d has 11%
  purity in every scheme and most of its weight is below the cut wherever
  you look, which the merge did not cause.  What the merge can cost is a
  configuration cut out of the very block it was assigned while it survived
  in another, and the symptom is a merged composition covering LESS than one
  block alone did -- the merge should normally cover more, taking each
  configuration from where it is most concentrated.  2016 of Pr III's 5626
  levels are in that position, which is the choice working as intended (the
  starred scheme spreads that configuration wider) and not an error; the
  remedy, if the weight is wanted, is a larger printout --components.

--------------------------------------------------------------------------
THE FILE
--------------------------------------------------------------------------

  par  Ee  Ec  dE_O-C  J  Lande_g  perc1 conf1 term1  perc2 ...

  Tab-delimited like levels.txt: these are read in a spreadsheet, and only
  coupling_choice.txt is edited by hand.  Parity first, as it is in all
  three phase 12-14 outputs -- it is wanted in the database load, and
  deleting a column for an article table is easier than reconstructing one.

  No cpl column: one line now carries every scheme, and which scheme a
  component is written in is implied by its label.  No e_c column either:
  levels.txt marks an unfitted level with an asterisk beside a copy of Ec,
  and here Ee and dE_O-C are simply blank together, which is what a
  published table wants and what a spreadsheet reads as missing.

  Ordered by parity then Ec, the parity blocks in the order levels.txt met
  them -- that is the order of IN36's configuration sets, and it is the
  author's.

  Options: --levels, --choice, -o, --components (5), --min-percent (1.5).
  The trimming is done here, after the merge; run printout with
  --components 10 --min-percent 0.5 and let prepub cut.

--------------------------------------------------------------------------
CODE AND TESTS
--------------------------------------------------------------------------

  toolbox/prepub.py               collect() pairs the blocks back into
                                  levels, apply_choice() is the rule,
                                  disagreements() is the self-check,
                                  shortfalls() the truncation report, and
                                  order/header/row_of/write_prepub the file.
  make_levels_list_prepub.py      a thin front end, as config_purity.py is.
  make_levels_list_prepub.bat     no Perl ancestor, so no _py suffix.

  tests/test_prepub.py      21.  The header defect and the unchanged
                            default header; the blocks paired back with
                            nothing missing; every component coming from
                            its own configuration's scheme and the counts
                            matching that block's; a Pr III level drawing on
                            all three schemes at once; descending order; the
                            invariance over all 5626 levels; the self-check
                            made to fire, and shown not to fire on a plain
                            per-block cut; the omitted configuration taking
                            the first scheme; a configuration in no deck
                            stopping the run; a star on a scheme with no
                            block stopping the run; trimming by count and by
                            percentage; the order; the columns; Ee and
                            dE_O-C blank together and consistent where they
                            are not; every level present; the options.

  417 pass in all, up from 396, with 5 skipped where a reference run is not
  unpacked.  Nothing that passed before changed.

==========================================================================
PHASE 14 -- DESIGNATIONS UNIQUE WITHIN EACH (PARITY, J)
==========================================================================

  levels.txt --12--> coupling_choice.txt --13--> levels_prepub.txt --14--> levels_pub.txt

The last of the three tools with no Perl ancestor.  Phase 13 wrote every
configuration's components in the scheme chosen for it; this names the levels.
A published or database level needs a name no other level of the same parity
and J carries, because that triple is what a citation, a line list and an ASD
load all key on.

The name comes from the first of these that works:

  1. the leading component;
  2. a component further down the composition;
  3. the leading component with a letter prefix, a 4F, b 4F;
  4. the energy -- for one level that has no component left, or for the whole
     file when naming has stopped meaning anything.

That progression is not a loss.  The label of a strongly mixed level has
little physical meaning in any coupling scheme, and the composition printed
beside it is what carries the physics.

--------------------------------------------------------------------------
WHAT IS PUBLISHED, AND WHY IT IS DECIDED FIRST
--------------------------------------------------------------------------

Uniqueness is required only among the levels that will actually be printed.
Demanding it of a calculated tail nobody will cite would spend letters and
fallbacks on levels nobody asked for, so the truncation happens before the
search.

--nlev takes one count per configuration set, in IN36 order.  RCE is
hard-limited to two sets (CONFIG(IC,2), RCE20KD.F:170), so in practice that is
one per parity, and levels_prepub.txt carries the parity rather than a set
number -- the sets are the parities in the order the file meets them, which is
IN36's order.  The default for a set is the position of its highest observed
level, which is where the calculation stops being anchored to measurement.

A set with no observed level has no such position.  Rather than publish all of
it or none of it, the run stops and names the set: only the author knows how
far an unanchored calculation is worth printing.  W70 is exactly that case --
no LSF, so nothing is observed -- and it is where the abort was tested.

--------------------------------------------------------------------------
THE SEARCH
--------------------------------------------------------------------------

Within one (parity, J) group each level offers its components at or above
--min-name-percent (20) as candidates, and an assignment is wanted that gives
every level a different name.  That is a bipartite matching and it is solved
as one.  The group is first split into clusters that share no candidate name,
which in a real spectrum leaves nearly every level a cluster of its own;
each cluster is then searched with a node budget, most-constrained-first,
bounded by the best the remaining levels could still be worth, and seeded with
a greedy pass so that a budget running out leaves a usable answer rather than
none.  Pr III's 1036 published levels take 0.15 s.

The objective is the total chosen percentage, with a strong preference for the
leading component (LEADING) under a much larger preference for naming a level
at all (NAMED).  There is no >50% pre-pass: a level with a large leading
component keeps it because the search prefers it, which needs no threshold and
cannot miss the case where two levels of one group lead on the same basis
state.  Every level is offered the same leading bonus, so where two of them
want one name the bonus cancels and the percentages decide which keeps it --
the one that loses less by giving it up.

--------------------------------------------------------------------------
THE LETTERS
--------------------------------------------------------------------------

A level the search leaves unnamed takes its leading component anyway, with a
letter prefix on the term.

BOTH sides of a contest take a letter, not only the one that lost it.  The
first run of this on Pr III put `4G*` beside `a 4G*` in one (parity, J) -- 18
of 27 lettered levels had an unlettered twin -- which reads as two different
things when it is two levels of one name that had to be told apart.  So the
levels that need letters are decided per (parity, J): if one level of a group
falls through to the name (conf, term), every level of that group with that
name is lettered.

The letters themselves count per (parity, configuration, term) across ALL J,
in increasing energy, and not within one (parity, J) group.  Within a group
the letter a level got would depend on how many same-named levels happened to
lie below it in that J, which changes from J to J as levels cross, so one
multiplet would come out `a 4F` at J=3/2 and `b 4F` at J=5/2.  Counted over
the family the letter is a property of the level.  Uniqueness within a
(parity, J) holds either way, since two levels of one group that want one name
are two members of one family and take two letters.

Pr III: 45 of 1036 published levels (4%) take a letter, in a/b pairs.

--------------------------------------------------------------------------
THE GLOBAL FALLBACK, AND WHAT THE UNITS DO TO IT
--------------------------------------------------------------------------

If more than --letters-max (50%) of the published levels need letters, or more
than --leading-min (50%) of them have a leading component under 50%, then no
naming scheme is going to mean anything and the tool stops pretending: every
level is named by the integer part of its energy -- Ee where known, Ec
otherwise -- with a blank term and a blank percentage, a sequence number
having no percentage to quote.  Both options are shares of the published
levels; both default to half of them.

Note what the integer part is of, which is the unit question of this phase.

RCE writes its energies in 10^3 cm-1 and every file up to levels_prepub.txt
carries them in it -- spectrum.py, radiate.py and explevs.KK all say so, and
the two reference cases agree once they are read properly: Pr3's Ee=82.492400
is Sugar's 82492.4 cm-1, and W70's Ec=132251.7250 is 1.32e8 cm-1, which is the
16 keV an n=2 to n=4 transition in Li-like W69+ ought to be.  (An earlier
draft of this section read W70 as already being in cm-1.  It is not; there is
one unit in this package and it is kK.)

levels_pub.txt is the file that gets published and the only one a database
loads, so it is the one converted: the energies are multiplied by publabel.CM1
on the way into it and are in cm-1 there, to three decimals, which is exactly
what the kK figures hold (6 decimals on Ee, 4 on Ec) -- nothing lost, nothing
invented.  Everything above it is left in kK, so the intermediates still diff
against what they always did.  The Ec quoted in this tool's own diagnostics is
converted too, so that a level named in a message can be found in the file.

The integer part of the energy is therefore the whole cm-1 a spectroscopist
would quote, and it separates levels at the resolution the list has.  In kK it
would not have: a thousand levels spread over 200 kK collide almost
everywhere, and the sequential index would have been the usual outcome of this
fallback rather than the rare one.

--------------------------------------------------------------------------
TWO DEFECTS IN PHASE 13 THAT THIS PHASE FOUND
--------------------------------------------------------------------------

Both are in what a component's configuration is, and both stopped a whole run.
Phase 13 was developed on Pr3, where neither can occur; W70 and U90 are where
they live.

  True jj splits every open l>0 subshell into its j = l-1/2 and j = l+1/2
  halves, so a TRUEjj component names `2p-2.2p+` where the deck's card, the
  choice file and conf0 all say `2p3`.  Phase 13 matched the label as written,
  found it in neither the choice file nor the deck, and stopped every W70 run
  with "levels.txt and IN36 are not from the same run" -- which was not true.

  A block merged from a reordered deck prints its subshells in that deck's
  summation order, so U90/reordered gives `2p2.1s` where the main deck says
  `1s.2p2`.  designate.py keeps summation order wherever an intermediate term
  makes it visible, which is right for a basis state and wrong for the name of
  the configuration it belongs to.  Same abort, on every U90 run with an LSr
  block.

purity.canonical_configuration() now does both: the halves are added back
together and the subshells put in (n, l) order, which is the order
deck_configurations builds its labels in.  A label with a token it does not
recognise comes back exactly as it went in, rather than half-canonicalised.
Component.bare goes through it, so phase 12's purity table and phase 13's
merge and self-check all see one configuration however its components happen
to be written -- which is the invariance those phases rest on, now applied to
the label and not only to the weight.  Pr III's levels_prepub.txt is unchanged
to the byte.

--------------------------------------------------------------------------
THE FILE
--------------------------------------------------------------------------

  par  conf  term  perc  perc1 conf1 term1 ... percN confN termN
       Ee  Ec  dE_O-C  J  Lande_g

  conf/term/perc are the chosen name and the weight it carries; the groups are
  the composition in descending order, up to --n_comp (3, between 1 and 5).
  The energies are in cm-1 here and in kK everywhere above -- see the units
  paragraph above.
  Tab-delimited, ordered as levels_prepub.txt was, parity first as in all
  three phase 12-14 outputs.  The column order is the plan's sentence read
  literally -- name, composition, then the energies, J and g that prepub
  carries.

  THE GROUP 1 COLUMNS ARE BLANK on every row named after its own leading
  component.  That is the ASD rule -- a level says once where its name came
  from -- and it makes the filled group 1 exactly the rows whose name came
  from further down.  Reading rule: the leading component is group 1 where it
  is filled and conf/term/perc where it is not.  Only group 1 is ever blanked,
  so a name taken from component 2 does appear twice, once as the name and
  once in its own group; the composition is printed whole.

  Two outcomes are accepted rather than refused, both the author's rulings.  A
  name may come from a component past the --n_comp printed groups: reported,
  with the levels listed, once it passes --max-outside (5%) of the published
  levels.  And a published conf may differ from RCE's conf0 -- silently, and
  necessarily so here, levels_prepub.txt not carrying conf0 at all.

  Uniqueness is not trusted, it is checked: duplicates() runs on the way out
  of every run and a collision stops it as a defect in the search.

  Options: --prepub, -o, --nlev, --n_comp, --min-name-percent, --letters-max,
  --leading-min, --max-outside.

--------------------------------------------------------------------------
CODE AND TESTS
--------------------------------------------------------------------------

  toolbox/publabel.py         cut() truncates, choose() searches, name_all()
                              adds the letters, by_energy()/fallback() the
                              last resort, and header/row_of/write_pub the
                              file.  duplicates() is the acceptance test.
  toolbox/prepub.py           gained read_prepub(): phase 13's file is phase
                              14's input, so it has to be read as well as
                              written.
  toolbox/purity.py           gained canonical_configuration(), above.
  make_levels_list_pub.py     a thin front end, as the other two are.
  make_levels_list_pub.bat    no Perl ancestor, so no _py suffix.

  tests/test_publabel.py   31.  levels_prepub.txt read back byte for byte and
                           an unfitted level still unfitted; the cut at the
                           highest observed level and --nlev overriding it; a
                           set with no observed level stopping the run; the
                           swap greedy cannot solve; the leading component
                           preferred; --min-name-percent withholding a
                           candidate; different J not competing; a group with
                           no assignment falling to letters, both sides of a
                           contest lettered, and the letters running across J;
                           a level with no component named by its energy; both
                           global fallbacks and both energy-label forms; every
                           one of Pr III's 1036 published names unique; the
                           columns; group 1 blank where the name came from it;
                           Ee and dE_O-C blank together; a name from outside
                           the printed groups; W70 end to end with a TRUEjj
                           block, U90 end to end with a reordered block, and
                           the script's options.
  tests/test_publabel.py   +1.  the energies written in wavenumbers, Sugar's
                           82492.400 among them.
  tests/test_purity.py     +3.  canonical_configuration on the true jj and
                           reordered labels, and a label it does not
                           recognise left alone.

  452 pass in all, up from 417, with 5 skipped where a reference run is not
  unpacked.  Nothing that passed before changed.

==========================================================================
PHASE 15 -- A PARAMETER FILE FOR conv_out
==========================================================================

conv_out has the worst command line in the package:

  conv_out lev.txt tr.txt u500 RCE LTE 100 1000 even.txt odd.txt IDEN linear

Two file names, then keywords, then two bare numbers that are told apart only
by the order they stand in -- so a missing argument shifts everything after
it, and nothing on the line says what 100 and 1000 are.  Long options were
added with the port (--temperature, --max-intensity, ...), but the settings
still had to be retyped in full on every run.

They can now be written once into a conv_out.par in the working directory, in
the layout calc_TP_M1E2.par already uses: the value first, then the key, then
prose saying what the parameter is for.

  100        ; temperature
  ;             Effective temperature, in the unit of the run -- 10^3 cm-1 in
  ;             every deck of this package.  Required in every model, and that
  ;             is the trap worth saying out loud: only LTE and LTEABS use it,
  ;             but the run stops without it whatever the model was.

--------------------------------------------------------------------------
PURELY ADDITIVE, AND HOW THAT IS ENFORCED
--------------------------------------------------------------------------

The par file supplies defaults and the command line overrides them.  The
historical positional form works exactly as it did, which matters more here
than anywhere else in the port: the conv_out tests drive that command line and
compare against conv_out.pl, and they are the only evidence the port is
faithful.  They pass untouched.

What makes it additive is an order, not a rule anyone has to remember.  Every
setting of Options starts as None -- "nobody has said" -- rather than at its
default.  The command line is parsed first and assigns what it mentions; then
_settle() fills from the par file only what is still None; then
Options.DEFAULT fills what neither source mentioned.  A setting typed on the
command line therefore already holds a value by the time the par file is
consulted, and can never be overwritten by it.

That change had one consequence worth naming.  _positional() used to end with

    if not self.ai:
        self.ai = names[:2]

which, with self.ai starting as [], ran unconditionally and left ai == [] on
every run.  Under the new order that is a value, and it would have shut the
par file's ai out.  It is now `if self.ai is None and names:`.

A parameter with no default gained a meaning in parfile as well.  Spec.show()
writes None as a blank field, and a blank field reads back as the default, so
a required parameter is written blank, reads back as None, and the program can
say what is missing.  levels, lines, model, temperature and max_intensity are
the five.

--------------------------------------------------------------------------
WRITE THE SAMPLE, DO NOT DEMAND IT
--------------------------------------------------------------------------

conv_out with no par file and no arguments writes a fully commented sample and
stops.  Not write-and-run: levels, lines, model, temperature and max_intensity
have no defaults, so there is nothing to run.

The trigger is narrow on purpose -- no par file AND no arguments at all.  A
command line that is merely incomplete is a mistake to be reported, exactly as
before.  Writing a file on any failed command line would have been wrong twice
over: it would turn a returncode 2 into a 0, and the comparison tests run
inside the reference decks of WORK/, which nothing in this package may write
into.  test_an_incomplete_command_line_is_an_error_not_a_new_par_file is the
guard on both.

calc_TP_M1E2 changed the same way.  It used to refuse and tell the user to run
--write-par; it now writes the file itself and says so.  In practice its
branch rarely fires, because CODE/calc_TP_M1E2.par is shipped and is found as
the fallback -- which is the whole point of shipping it.

conv_out must NOT ship a conv_out.par in CODE/, or the file would never be
absent and the sample would never be written.  That is what
parfile.find(..., shipped=False) says, at the one call site that needs it.

--------------------------------------------------------------------------
ONE COPY OF THE FINDING RULE
--------------------------------------------------------------------------

calc_TP_M1E2.find_par was the pattern -- the named file, then the working
directory, then the copy beside the code -- and it is now parfile.find(), used
by both programs, so phase 11 does not leave two idioms behind.  It returns
None rather than raising when there is nothing to find, because the two
programs do the same thing with that answer and it is not an error in either.
A file named outright that is not there IS an error: the user asked for that
file and no other.

--------------------------------------------------------------------------
THE PROSE IS THE DELIVERABLE
--------------------------------------------------------------------------

Three parts of it are ten years unused and had to be reconstructed from the
sources rather than remembered.

  RCE/noRCE, LTE, LTEABS, BF, IDEN, truncate, linear and convert2log come
  from the author's own usage text, conv_out.pl:72-85.

  CM, the cascade matrix.  The Perl says only "CM is 'Cascade Matrix'
  option".  What it is, is in conv_out.py _cascade(): a square table of
  A(i->j) over the levels, levels down the side and across the top.  Two
  behaviours in _check() follow from that and are in the prose: the file is
  comma separated rather than tab separated, because it is read in a
  spreadsheet, and only the first coupling is used, a matrix having one axis.
  A third belongs with them: the level list's n, l, occ, S and L columns
  classify a level for this model and nothing else reads them, so they are
  written with CM and with no other model.

  DR and AUGER, and why ai is two files.  Both models are built from
  autoionization rates and both require ai.  The rates come as two files
  because of how they are made: read_aa is run in two separate directories,
  one holding the even configurations and one the odd, and each writes the
  rates of the set it was run on.  No source states that working arrangement;
  it is recorded here and in the sample's prose.

flight is the one setting spelled differently in the two places.  On the
command line the time of flight has nowhere to go but the model keyword,
BF:0.114; the par file has a field for it, so `model = BF` is right there and
would be rejected on the command line.  _par_values() therefore assigns the
par file's model directly instead of putting it through _model().

--------------------------------------------------------------------------
CODE AND TESTS
--------------------------------------------------------------------------

  toolbox/parfile.py       gained find(), the shared finding rule, and a Spec
                           whose default is None -- written blank, read back
                           as None.  No new value kinds: model is a TEXT with
                           choices, and couplings, drop_full_shells and ai are
                           comma-separated TEXT split exactly as
                           Options._option splits the same text.
  conv_out.py              SPECS and TITLE, _par_values(), Options._settle(),
                           _par_arguments()/_sample()/_parameters(), and
                           --par / --write-par.  No par file is shipped.
  calc_TP_M1E2.py          find_par() deleted in favour of parfile.find();
                           a missing par file is written rather than demanded.

  tests/test_conv_out.py   +13.  the sample written on a bare run and read
                           back through its own reader; an incomplete command
                           line still an error and no file left behind;
                           --write-par; every one of the six models named and
                           explained in the sample; every parameter
                           documented; the par file supplying what the command
                           line does not; the command line having the last
                           word, including the bare keyword RCE over
                           no_rce=1; the built-in defaults where neither
                           source speaks; a bare BF with its flight beside it;
                           and an unknown key and a bad value each named with
                           their line.
  tests/test_parfile.py    +8.  find() on all four outcomes and shipped=False;
                           a parameter with no default written blank, read
                           back as None, and still range-checked.
  tests/test_calc_tp.py    +1.  a missing par file written rather than
                           demanded.

  473 pass in all, up from 452, with the same 5 skipped.  Nothing that passed
  before changed.

==========================================================================
PHASE 11 -- RETIREMENT
==========================================================================

The Perl is out of the package.  Thirteen .pl files moved to
CODE/perl_legacy/, the eleven X_py.bat launchers became X.bat, and README.md
in the repository root now describes the Python.  DOCS/ is untouched, as ever.

Nothing in the port changed in this phase.  The suite was 473 passed and 5
skipped before the move and is 473 passed and 5 skipped after it, the same 5
being add_eav.pl refusing to open `ing11.` under this perl -- now reported
against its new path, which is the evidence the comparison tests really are
reaching perl_legacy/ and not quietly skipping.

--------------------------------------------------------------------------
THE MOVE
--------------------------------------------------------------------------

  Copied, verified byte-identical, and only then removed, because nothing in
  this repository is committed and a .pl deleted by accident is gone:

    add_eav      calc_TP_M1E2   conf_list    conv_out     conv_cowan
    printout     read_aa        reorder_ing11             s11
    scale_param  subst_exp_levs subst_params_ing11        vacair

  Thirteen in, thirteen out, and zero .pl left in CODE/.  conv_cowan.pl was
  conv_out.pl's own library and had to move with it; vacair.pl never had a
  bat and is now toolbox/vacair.py.

  tests/conftest.py gained one constant, PERL = CODE / 'perl_legacy', and the
  six call sites go through it -- four of them, because three of the six are
  the one perl() helper in test_rewriters.py.  test_conv_out.py's `-I` moves
  with the rest, or conv_out.pl cannot find conv_cowan.pl.  These comparison
  tests are the evidence behind every phase from 6 onward and had to keep
  running after the move, not be retired with the programs they compare
  against.

--------------------------------------------------------------------------
THE BATS
--------------------------------------------------------------------------

  Eleven X_py.bat became X.bat.  The Perl form's -x, the %COMSPEC% test and
  the goto script_failed_so_exit_with_non_zero_val label did not carry over.
  `python` became `py`, the launcher that is reliably present on Windows, in
  all fourteen -- the eleven above and the three phase 12-14 programs, which
  had never been through a swap and still said `python`.

  11a was already settled and nothing was done to it: every front end opens
  with sys.path.insert(0, dirname(abspath(__file__))), so a bat needs nothing
  but the interpreter and %~dps0, and PYTHONPATH need not be set at all.

  11b was declined by the author: in36, outg11, genealogy and jjkeys stay
  library modules with an inspection entry point, run as `py -m toolbox.in36`.
  Four more bats in CODE/ would be four more names on the system PATH, since
  CODE/ is copied to C:\cowan\code.  The -m form is documented in README.md
  instead.

--------------------------------------------------------------------------
A DEFECT IN THE SETTLED BAT FORM, FOUND AND FIXED HERE
--------------------------------------------------------------------------

  The two-line form the plan called settled -- the program, then

      if %errorlevel% == 9009 echo You do not have Python in your PATH.

  -- always exits 0.  An `if` whose condition is false clears ERRORLEVEL in
  cmd, so the trailing test threw away the exit code of every command in the
  package.  Measured, not reasoned: make_levels_list_pub returns 1 on a
  missing levels_prepub.txt, and returned 0 through its bat.  The same bat
  with the `if` removed returned 1.

  This mattered.  update11.bat chains s11 -> rcg -> printout and could not
  see a failed step, and neither could any batch file a user writes.  The
  Perl bats had explicit scaffolding for exactly this, which the plan
  correctly dropped as Perl-era -- but nothing replaced the thing it did.

  Note that `exit /b %errorlevel%` on a line after the `if` does not fix it:
  cmd expands the whole line at parse time, after the `if` has already
  cleared the variable.  The exit code has to be saved before the test:

      setlocal
      py "%~dps0X.py" %*
      set _rc=%errorlevel%
      if %_rc% == 9009 echo You do not have Python in your PATH.
      endlocal & exit /b %_rc%

  `endlocal & exit /b %_rc%` is one line, so %_rc% is expanded before endlocal
  discards it; setlocal keeps _rc out of the caller's environment, and it
  does not leak.  Verified: 1 on failure, 0 on success, and no _rc left
  behind.  All fourteen bats carry it.

--------------------------------------------------------------------------
ACCEPTANCE
--------------------------------------------------------------------------

  Each of the fourteen commands was run from a working directory that is
  neither CODE/ nor on PYTHONPATH, with PYTHONPATH empty, on a copy of
  WORK/Pr3 -- a real run, not a --help:

    conf_list                the 21 configurations of both parities
    printout                 16878 levels in LS, JJ and LSJK, the last from
                             Pr3/JK, with conf0
    config_purity            39 configurations; LS for 19, JJ for 13, LSJK
                             for 7
    make_levels_list_prepub  5626 levels
    make_levels_list_pub     1036 published, 45 (4%) taking a letter -- the
                             numbers phase 14 recorded
    s11                      665 + 542 parameters substituted
    add_eav                  39 Eav parameters shifted
    scale_param              404 parameters scaled in RCEINP
    reorder_ing11            3 of 8 slots moved
    subst_params_ing11       166 cards carried across
    subst_exp_levs           RCEINP into RCEINP.HF, with its match report
    conv_out                 11252 levels, 120274 transitions; and, with no
                             par file and no arguments, the sample written
                             and the run stopped
    calc_TP_M1E2             --write-par, 6975 bytes
    read_aa                  reached the OUTG11 reader and reported, rightly,
                             that this run has no autoionization rates --
                             Pr3 has none, and the AI reference case has no
                             OUTG11 to read

  Every .pl is present in perl_legacy/ and none in CODE/.  473 pass, 5 skip.

--------------------------------------------------------------------------
README.md
--------------------------------------------------------------------------

  Targeted replacements, not a rewrite: everything the author wrote about the
  physics, the file formats and the calculation chain is still true and was
  left as it stands.  What changed:

  - Strawberry Perl is now stated to be needed only for the gfortran
    toolchain and its runtime libraries; the utilities need Python 3 and the
    standard library, and nothing to pip install.  Installation step 0 asks
    for Python first and Perl only if the Fortran is to be recompiled.
  - The CURRENT VERSIONS table lists the fourteen Python programs and
    toolbox/, and says where the Perl went.
  - RCN.BAT/RCN2.BAT/RCG.BAT/RCE.BAT are written in their real lowercase, and
    mark_time.bat is named.
  - Section f.1 describes running a .py directly, the docstrings, the four
    `py -m toolbox.X` inspection entry points, tests/ and perl_legacy/.
  - printout's entry drops "works more-or-less well with LS and JJ coupling
    only" -- it is no longer true -- and gains --couplings, --from (both
    forms), --dominant-config, --drop-full-shells, --components,
    --min-percent and --match.
  - conv_out's entry gains the long options and conv_out.par, and its
    description of the two autoionization files now says what read_aa is and
    why there are two of them, in place of "it is rather difficult to give a
    proper description of this code ... email for advice".
  - calc_TP_M1E2's entry says its constants are in calc_TP_M1E2.par.
  - vacair.pl's entry becomes toolbox/vacair.py.
  - A new section f.2 documents config_purity, make_levels_list_prepub and
    make_levels_list_pub, the file each writes, the invariance they rest on,
    and the printout command that prepares their input.

  The date at the foot is now August 18, 2026.

--------------------------------------------------------------------------
CODE AND TESTS
--------------------------------------------------------------------------

  CODE/perl_legacy/        13 .pl files, unchanged, unused, not on the PATH.
  CODE/*.bat               14 launchers rewritten: py, the exit code, and a
                           one-line description of what the program does.
  tests/conftest.py        PERL = CODE / 'perl_legacy'.
  tests/test_conflist.py   \
  tests/test_conv_out.py    | one line each, routing the Perl reference run
  tests/test_explevs.py     | through PERL.  test_conv_out.py's -I too.
  tests/test_rewriters.py  /

  No test was added, removed or rewritten.  473 pass and 5 skip, exactly as
  before the move.


==========================================================================
PHASE 16 -- THE WAVELENGTH-SORTED SPECTRUM, AND DMIN
==========================================================================

Two corrections to WORK/Forbid, and one new layout in the OUTG11 reader.


DMIN: WHY WORK/Forbid PRINTED NO M1 SPECTRUM

  WORK/Forbid/README.txt used to say that no M1 spectrum comes out of the
  deck even with IMAG set, and that it was not known why.  It is known.

  Columns 61-65 of the IN2 main card are DMIN, and rcg11kd.f:12859-12871
  reads a DMIN between 1 and 50 as "omit every transition with
  gA < 10**DMIN".  The deck carried 5.00, which is a reasonable cut for
  electric-dipole work and is what the other decks in WORK/ carry -- 5.00,
  9.00, 10.0.  Every M1 line of Sn VIII 4d7 is weaker than gA = 1e5, so
  SPECTR was called (rcg11kd.f:11622), found nothing above the cut, and
  printed the heading with no lines under it.  Nothing was wrong with RCG.

  0.00 is the right value for any forbidden-transition calculation, and is
  RCG's own default (rcg11kd.f:836).  The deck now carries it, and prints
  89 M1 lines and 131 E2 lines over 132 distinct level pairs, 88 of which
  carry both.  That is what makes the M1 + E2 merge -- the thing
  calc_TP_M1E2 exists for -- testable on a real run; before, it was
  covered only by unit tests on rows built up by hand.

  test_calc_tp.test_a_whole_run_with_no_trials asserted that every row of
  tp_M1E2.txt was E2, which was a true statement about a deck that could
  not print M1.  It now asserts that all three kinds -- M1, E2 and M1+E2
  -- come out of it.


ISPECC: THE SECOND SPECTRUM LAYOUT

  Column 72 is ISPECC, which chooses RCG's three line sorts.  The two level
  sorts print the Lund layout, a `* * *` line for the initial level and one
  line per final level, and that is the only layout toolbox/spectrum.py
  read.  The wavelength sort prints the wide layout instead -- both levels
  on every line -- whatever ILUND says, at rcg11kd.f:14019.  So an ISPECC
  of 4 or more produces it, and an ISPECC of exactly 4 produces nothing
  else.  Every deck in WORK/ is ISPECC = 2, which is why it had never come
  up.

  The old refusal message said the wide layout does not carry the level
  numbers NLEVJ and NLEVJP that tie a line to an eigenvector.  That was
  true of Cowan's FORMAT 58 and stopped being true when A. Kramida added
  I6,I6 to it in 2018 (rcg11kd.f:13473, and the ORDERI calls at :13977).
  The wide layout carries both, and carries them in more digits than the
  Lund one does: F11.4 energies and an E11.4 gA against F11.3 and E10.3.

  toolbox/spectrum.py now reads both.  What it took:

  - The layout is decided per section, by the column header the section
    prints -- FORMAT 1056 for Lund, FORMAT 56 for wide -- and not by
    guessing at field widths.  A file may hold both: ISPECC = 7 prints
    three sorts of the same multipole, each under its own FORMAT 55
    heading.

  - The wide layout prints no configuration name beside either level.  Both
    come from the heading's own table, the tail of FORMAT 55, which numbers
    the configurations of each parity set exactly as the CONF and CONFP
    columns index them (ELEM(L,JN,L1) and ELEM(L,JN,L2), rcg11kd.f:14004).
    So Line.initial.label and Line.final.label are both filled, where the
    Lund layout fills only the first.

  - It truncates the 8-character basis-state label to seven (A7).  The
    eighth character is blank on every label in every reference run --
    checked over Pr3's 120274 printed lines, W70's 2205 and U90's 66 --
    so nothing is lost, but the field is the shorter one.

  - Three layouts remain unread and are now refused by name rather than by
    a message about level numbers: FORMAT 57, the photoionization
    cross-section form (CSIGMA); FORMAT 51/53, the Boltzmann-weighted one
    (TEXCIT); and FORMAT 88, the branching-ratio variant of the wavelength
    sort (IBRNCH, rcg11kd.f:13071), which puts its first label in A8 and so
    shifts every column after it by one.  _is_wide rejects that shift
    rather than reading FORMAT 88 with FORMAT 58's offsets.

  Checked on WORK/Forbid built twice, once with ISPECC = 2 and once with 4:
  220 transitions each, the same 132 level pairs, and a tp_M1E2.txt with
  the same rows and transition probabilities agreeing to 4e-4 -- the extra
  digit.  That is
  test_calc_tp.test_a_wavelength_sorted_run_gives_the_same_table.

  Checked on a two-parity cross-parity run as well, WORK/'s Sn VIII deck
  built the same two ways: conv_out reads both, 199 levels and 1695
  transitions from each, and lines.txt differs in 76 of 509 rows, every one
  of them in the last digit of A.

  test_spectrum.py pins the wide columns against an excerpt copied verbatim
  out of the ISPECC = 4 run, with a test that the excerpt still matches the
  reader's own constants -- otherwise the column tests could pass for the
  wrong reason.


DOCUMENTATION

  WORK/Forbid/README.txt gains the DMIN explanation, as a mistake to be
  avoided rather than a property of this deck, and a section on the two
  layouts.  Its opening no longer says the deck is E2-only.

  README.md's WORK addendum now lists all seven examples the distribution
  carries -- WORK/IN36 and IN2, AI, AI/DIEL, Forbid, Pr3, W70, U90 -- with
  a line on what each is for, in place of the two lines that described only
  the Sn VIII installation test.  calc_TP_M1E2's entry gains the IMAG,
  IQUAD and DMIN card columns its input depends on.  The date at the foot
  is now August 19, 2026.

  479 pass and 5 skip: six more than phase 11, five in test_spectrum.py and
  one in test_calc_tp.py.


CFP DECK REPRODUCIBILITY (stage 8a of the Unix port)

  The question was whether FOR072/FOR073/FOR074 come out of make_cfp with
  non-reproducible bytes in them, because CODE/FOR074 as committed differed
  from a fresh regeneration in eight bytes and nobody had explained why.
  If the decks were genuinely unstable there would be no way to tell a real
  regeneration difference from garbage, and the Unix make_cfp could not be
  verified against them.

  They are not unstable.  make_cfp was run twice from a copy of CODE with
  the committed rcg11kd.exe: FOR072, FOR073 and FOR074 are byte-identical
  between the two runs, and FOR072 and FOR073 were already byte-identical
  to the committed copies.  Nothing in LNCUV needed zeroing -- the
  Cray-pointer heap arrays were the suspects and they are not the cause.

  The eight bytes are at offset 307776 of FOR074, in the last block of the
  file, h2 (l = 5).  Read as a record they are the AA field -- the
  CHARACTER*8 parameter label -- of the fourth and fifth records written by
  the DO 492 loop of CUVFD:

        DO 492 I=1,IMAX          IMAX = FLL(1)
          AA=ABGNAM(J)
          ...
          WRITE (ID4) KPAR,K,AA,CAVEI(J),NOPC, (PC(N1), N1=1,NOPC)

  IMAX is FLL(1), so the loop runs l times, but ABGNAM holds three names,
  'ALPHA', 'BETA', 'GAMMA', and CAVEI three values.  For l <= 3 that is
  exact.  For g (l = 4) the loop makes one pass past the end of both
  arrays, and for h (l = 5) two.  The label written is then whatever static
  storage follows ABGNAM.  That address is fixed for a given executable, so
  the bytes are the same on every run of the same build -- which is why the
  regenerations agree -- but it moves when the compiler relays out statics,
  which is why the committed deck, made by an older build, carried
  'ZETA 1  ' where the current build writes 07 00 00 00 00 00 00 00.

  The label is not the only thing that overruns.  ABG is dimensioned
  ABG(KLSI,3) in COMMON/C3C/, and the same J indexes it at

          PC(NOPC)=ABG(M,J)

  so the extra passes write ABG(M,4) and ABG(M,5) as well -- that is, the
  bytes that follow ABG in the common block: CAVEI(1..3) and then NALSJI.
  Dumping the deck shows exactly that.  The g2 block carries a fourth
  parameter record whose PC values are -37.6471, -1.23077, -1.07563 --
  which are CAVEI(1), CAVEI(2), CAVEI(3), the ALPHA, BETA and GAMMA
  averages of the three records just above it -- followed by NALSJI's
  integers read as REAL*8 denormals of order 1e-313.  The h2 block carries
  the same fourth record and a fifth whose PC values are all zero.  None of
  it is physics.

  The two overruns behave differently across builds, and that is the whole
  of the eight bytes.  ABG's overrun lands inside a COMMON block, whose
  layout the standard fixes, so it produces the same bytes in every build.
  ABGNAM is a local static, and what follows it is the compiler's choice,
  so the label moved: 'ZETA 1  ' from the old build, 07 00 00 00 00 00 00
  00 from this one.

  The records are read back.  NPARI = NORCD-2 in PFGD, and NORCD was
  written as 1+NPAR with NPAR = 1 + 2*FLL(1), so the count already budgets
  l parameters of the alpha/beta/gamma family for an l(n) shell.  That is
  the real shape of the defect: RCG reserves l slots but only three such
  operators are defined.  For d the fit is exact once line

          IF (LL(1).EQ.LSYMB(3).AND.I.EQ.2) J=3

  sends the second slot to GAMMA's column; for f it is exact as written;
  for g there is one slot too many and for h two.  A g2 or h2
  configuration run with IABG > 0 therefore acquires one or two extra
  parameters with garbage names -- reaching PARNAM at label 188 of PFGD --
  and values of 1e-313 or zero.  They perturb nothing numerically, but they
  are there.  Nothing in the distribution's examples fits alpha/beta/gamma
  on a g or h shell, and Pr III has no such shell at all, which is why two
  Pr III runs had already shown the difference to be harmless without
  explaining it.  (FOR074 has blocks only through g2 and h2, so no
  higher-occupancy case exists to worry about.)

  Fixing this is not a matter of clamping IMAX to 3: NPAR would have to
  change with it, and that changes the deck format for g and h.  The
  question underneath it is which effective configuration-interaction
  operators the fourth slot of a g shell and the fourth and fifth of an h
  shell ought to be -- Cowan defined alpha, beta and gamma and stopped.
  That is physics, not porting, and there is no settled answer in the
  literature to look it up in; open g and h shells do not occur in
  practice, which is why nobody has needed one.

  WHAT WAS DONE INSTEAD (A. Kramida's decision, 19 August 2026)

  Write the missing slots blank, and refuse any run that would actually
  use them.  Two changes, neither of which moves a byte for any other
  shell:

    * CUVFD.  The DO 492 loop sets J=0 for I>3 and writes AA blank,
      CAVE zero and every PC zero, instead of reading past the end of
      ABGNAM, CAVEI and ABG.  The record COUNT is untouched, so the deck
      format and everything that counts records (NPARI = NORCD-2 in
      PFGD, NPARI = NORCD-1 on the zeta path) are unaffected.  For
      l <= 3 the loop takes the same branch as before, expression for
      expression, so d and f are bit-exact.

    * SUBROUTINE ABGCHK, called from LNCUV beside NLTCHK, once per
      configuration card while the card can still be named.  With
      IABG >= 1 a subshell of l >= 4 and 2 <= n <= 4l stops the run:

          *****CONFIGURATION ERROR, PARITY SET 1, CONFIGURATION 2 (...):
                THE SUBSHELL G** 2 CANNOT BE HANDLED WITH IABG.GE.1.
                ...
                PUNCH IABG=0 IN COLUMN 21 OF THE MAIN CARD TO RUN
                WITHOUT THEM.

      IABG=0 is deliberately left working: PFGD sets NPARI=FLL(I) there
      and reads the F**k records only, so g**2 and h**2 without
      alpha/beta/gamma are an ordinary calculation and are not refused.
      Verified both ways on a Pr III deck with one card's first subshell
      changed to G 2 -- IABG=2 stops in ABGCHK naming the card, IABG=0
      passes it and stops later for an unrelated reason of the doctored
      deck.

  FOR074 was regenerated with the new CUVFD and is committed:
  d1f3a5a724663106ad6d9323461e566b, 308372 bytes as before.  It differs
  from the previous deck in 106 bytes, all of them inside the four
  parameter records at the ends of the g2 and h2 blocks; FOR072 and
  FOR073 are unchanged.  CODE/rcg11kd.exe was rebuilt from the new
  source with gfortran_O3_OMP.bat's flags -- one warning, the expected
  -fno-automatic one.

  Regression: WORK/Pr3/JK/ING11 run twice, once with the previous
  executable and the previous decks and once with the new ones.  TAPE2E
  7fc8dc27eb58bcc30467c3073982af7a, OUTGINE b4c556904b9a3a92408fd4a0998d8634
  and genealogy.txt 14d00e55a4c025da09cfdf9e46fc68e5 are identical, and
  OUTG11 is identical once the timing lines are stripped.  (Pr III has no
  g or h shell, which is the point: nothing that any real deck reads has
  moved.)

  CODE/FOR074 has been refreshed from the clean regeneration, so the three
  decks in the repository now match what the committed executable produces.
  No Fortran was touched.  Pr III was run against the refreshed decks:
  TAPE2E md5 062a5445d38519c6f5be44bfd55c0965 and OUTGINE md5
  2d904d54ca00b80ddadc508210c6e0db, both unchanged, and OUTG11 identical
  once the timing lines are stripped.


END-TO-END VERIFICATION ON LINUX (stage 11 of the Unix port)

  Nothing in stages 1-10 proved the package runs on Linux, because until
  this point no Linux machine had been involved.  It was run under WSL 2,
  Ubuntu 24.04.2 LTS, gfortran 13.3.0, python3 3.12.3, x86-64, 32 cores.
  The tree was copied to ~/cowan-src first: /mnt/f is DrvFs and reports
  every file 777, which would have made the permissions test meaningless.

  Build.  build.sh produced all four programs in 14 seconds, with the
  expected -fno-automatic/-frecursive warning and nothing else.  The
  -fno-automatic .bss is 110 MB for rcg11kd and 323 MB for RCE20KD, far
  under the 2 GB at which the small code model runs out and
  -mcmodel=medium becomes necessary.  libgfortran and libgomp come from
  the system; nothing is shipped beside the executables, as it is on
  Windows.

  Install.  install.sh laid down four executables, fourteen tools,
  twenty-seven toolbox modules, seven data files and twenty-five
  wrappers, and no .bat, .exe, .dll, tests/ or perl_legacy/.  On a real
  Linux filesystem the modes are what the script intends, including the
  644 on the .py tools that MSYS could not show.  make_cfp ran through
  its permission dance and regenerated FOR072, FOR073 and FOR074
  BYTE-IDENTICAL to the committed decks -- which is both
  the confirmation of stage 8a's determinism finding and the proof of
  install.sh's claim that the decks are portable as bytes.  RCG then read
  them at mode 444 without complaint, confirming that libgfortran retries
  a refused read-write open (STATUS='UNKNOWN') read-only; a reinstall put
  them back to 644.

  Pr III, Linux against Windows, the same deck run through the whole
  README sequence on both:

      TAPE2E    byte-identical, 062a5445d38519c6f5be44bfd55c0965
      OUTGINE   byte-identical
      OUTG11    84 MB, 170 differing lines, every one a timing line
      OUT36     differences are timing lines only
      OUT2      differences are timing lines only

  All fifteen final products are identical once CR is stripped: ING11,
  RCEOUT, RCEINP.HF, genealogy.txt, levels.txt, params.txt,
  levels_prepub.txt, levels_pub.txt, levels_13kk.txt, tr_13kk.txt,
  ENLEV.DAT, TRANS.DAT.  The one exception is coupling_choice.txt, which
  differs in its "# folder:" provenance comment and nowhere else.  Every
  number WORK/Pr3/README.txt documents came back: AVDEV = 0.133; 25 J
  blocks and 5626 levels; 16878 levels written with the LSJK J=12.5 note;
  39 configurations split LS 19, JJ 13, LSJK 7, with 15 within 2%; 2016
  levels merged; 1036 published (o 481, e 555) and 45 levels (4%); 795
  transitions, 11252 levels, ENLEV.DAT 555 levels and TRANS.DAT 792
  lines.  Thread invariance holds on Linux at 1, 4 and 32 threads for
  both RCG and RCE.

  TAPE2N does differ between the platforms, and only TAPE2N.  With the
  record layout read out of RCN2K.F -- 4-byte markers, a 1306-byte
  header, then R(1801), RU(1801) and PNL(1801,NPNL) as doubles, then IW6
  -- the relative difference is 1.2e-13 over values above 1 and 1.3e-9
  over all of them.  That is ordinary last-bit rounding in RCN's radial
  wavefunctions between two gfortran versions, and it does not propagate:
  RCN2 writes the parameters to six figures, and TAPE2E downstream is
  byte-identical.  Do not chase it.

TWO DEFECTS STAGE 11 FOUND

  RCE indexed past the end of COMMON/C1/.  RCE died with SIGSEGV on the
  first J block of every Linux run.  Threading and stack size were ruled
  out (OMP_NUM_THREADS=1 and ulimit -s 524288 both still failed); a debug
  build put it at RCE20KD.F:511 in MAIN__.  gfortran's own backtrace
  handler pre-empts ASan, so the loop was instrumented with a
  write(0,...) instead, which showed N reaching 201 against NPAR(IC,IC)
  with IC=200.

  In the DO 3145 loop, N is a configuration index and was incremented
  with nothing bounding it.  Once the real configurations are past,
  NPAR(N,N) reads zero, the test fires on every pass, and N advances
  until it indexes off the end of the block.  Windows walks out of bounds
  too -- an instrumented MinGW build reaches N=201 as well -- but the
  word beyond happens to be large there, so the runaway stops and the run
  finishes.  Under Linux it is small, so N runs on until the read leaves
  the mapped region.

  The guard is IF (N.LT.IC) around the existing test, nested rather than
  written with .AND. because Fortran does not guarantee short-circuit
  evaluation.  It cannot change a number: NMAX is 21 here, the test at
  the head of the loop takes the UNPK2 branch whatever N is once N
  exceeds NMAX, the ELSE branch is reachable only while N is at most
  NMAX, and N is used nowhere else in the loop body.  That reasoning was
  checked rather than asserted -- the guarded source was built on Windows
  with the production flags and run against the committed executable on
  Pr III: RCEOUT, TAPE2E and OUTGINE byte-identical, OUTE (73 MB)
  differing in 22 lines, all of them timing.  CLAUDE.md section 1's
  standard, met.  CODE/RCE20KD.exe was rebuilt from the fixed source and
  reproduces the same result.  The reasoning is at the point of change
  and in the header change list of RCE20KD.F.

  The Python shims wrote __pycache__ into the installed package.  All
  fourteen tool wrappers in CODE/unix exec'd the interpreter without -B,
  so the first run of any tool dropped a .pyc beside it -- the one breach
  of the guarantee the whole install layout exists to make, that the code
  directory is never written to during a run.  All fourteen now pass -B,
  with the reason in a comment.  Recompiling these few modules on each
  run costs a few milliseconds.

  Not done: the pytest suite was not run on Linux.  pytest is not
  installed there, and on Ubuntu 24.04 neither pip nor python3 -m venv
  works without sudo.  The tools were exercised end to end on real Pr III
  data with byte-matching output, which is a stronger check of the same
  code, but it is not the unit suite and should not be reported as it.

  Not done: macOS.  Nothing in the package is known to depend on
  anything macOS lacks, and build.sh, install.sh and _common.sh handle it
  deliberately -- BSD readlink, the zsh default shell, Homebrew gfortran
  -- but no run on macOS, Intel or Apple Silicon, has been made.
  README.md says so in as many words.

DOCUMENTATION FOR THE PORT (stage 10)

  README.md stops describing itself as a Windows package.  The title, the
  system and hardware and software requirements, the "running the
  programs" section and the addendum are now written for all three
  systems; INSTALLATION became INSTALLATION ON WINDOWS and a new
  INSTALLING ON LINUX AND macOS sits beside it, covering gfortran and
  python3 from apt, dnf or brew, build.sh, install.sh, the PATH line for
  ~/.bashrc or ~/.zshrc, and the two Unix-only settings (OMP_STACKSIZE
  and .bss size) that install.sh prints.  RECOMPILING UNDER WINDOWS
  became RECOMPILING THE FORTRAN with a subsection per system and a third
  on the three flags that are load-bearing.  Package-relative paths are
  written with forward slashes throughout; the backslashes that remain
  are all inside genuinely Windows-only instructions.  The utility codes
  are listed by command name rather than by .bat name, since the name is
  now a .bat on one system and a shell script on another.

  The COMPILING FOR OTHER OPERATING SYSTEMS section was the one the
  author asked to have rewritten, and it says two things it did not say
  before.  First, that the MSDOS flag is not a Windows detail to be
  removed but the mechanism of the code/workspace split, that it is
  implemented and verified for Windows, Linux and macOS with gfortran and
  for no other combination, and what to do about it elsewhere (MSDOS=0
  and the hard-coded FIL=).  Second, a warning that building with another
  compiler drops the parallelization -- the threading is OpenMP as
  gfortran implements it -- and that on a large case this is worth more
  than a factor of ten in wall-clock time on a many-core machine.

  The WORK walkthroughs gained a note that their command lines are
  written as they are typed on Windows and are lower case on Linux and
  macOS, and that "copy" reads "cp" there: WORK/Pr3, WORK/U90, WORK/W70
  and WORK/Forbid.  WORK/AI/DIEL/README.txt had a command block naming
  executables that no longer exist under those names (rcn36.exe,
  rcg11.exe) and hand-writing COWAN.CFG; it now names the wrappers, says
  why the wrapper and not the executable, and gives the gcard
  concatenation for both shells.


EVERY FUNDAMENTAL CONSTANT IN ONE FILE (FOR/fundconst.f)

  The package carried its physical constants as literals scattered through
  the Fortran: 109737.31 in one routine and 109737.0 in another, alpha as
  137.036 here and 137.04 there, a0 as 0.529, pi as 3.14159.  Most of them
  were the values of the late 1960s.  This work replaced all of them with
  named parameters from a single include file, FOR/fundconst.f, and brought
  the values to CODATA 2022.

  The rule the author set, and it is written at the head of fundconst.f:
  every fundamental constant used anywhere in the package belongs in that
  file.  There is no threshold of smallness below which one may stay
  inline.  A literal left in the Fortran is a literal the next change of
  these values cannot reach, and reaching all of them is the whole point.
  CODE/tests/test_literals.py enforces it over FOR/ and keeps no exception
  list.

  The full plan, with every measurement, is PLAN_fundconst.md in the
  repository root.  What follows is the account.

THE THING THAT STARTED IT, WHICH TURNED OUT NOT TO BE A CONSTANT

  WORK/gA_gf_check.txt shows the gA printed by RCG disagreeing with the gA
  recomputed from the printed LOG GF by about 1e-3, with random sign.  This
  has been looked into more than once and the answer keeps being forgotten,
  so it is written down here in full.

  It is quantization of the three-decimal LOG GF print field, and nothing
  else.  Half a unit in the last decimal of a base-10 logarithm is
  ln(10) * 0.0005 = 1.151e-3 in relative gf.  Recomputing the exact log gf
  implied by the printed gA and dE, with the modern
  8*pi^2*alpha^2*a0*c = 0.66702518:

    sigma (cm-1)   log gf printed   log gf implied     delta   half-ulp
      418430          -1.485          -1.48521       +0.00021   inside
      411695          -2.049          -2.04940       +0.00040   inside
      408365          -2.745          -2.74521       +0.00021   inside
      397156          -2.285          -2.28511       +0.00011   inside
      396287          -3.276          -3.27563       -0.00037   inside
      384932          -4.247          -4.24688       -0.00012   inside
      368063          -4.965          -4.96458       -0.00042   inside

  All seven are inside +-0.0005.  The printed LOG GF is the correctly
  rounded value of the gf that is exactly consistent with the printed gA.
  Every observed feature follows: random sign, about 0.07% rms, mean near
  zero, the second decimal affected.

  So gA is the primary quantity in an RCG transition list and log gf is
  not.  The naive expectation is the opposite, because gf comes from S with
  a lower power of the transition energy (sigma rather than sigma^3), so
  energy rounding propagates into it less strongly.  That is true and it is
  beside the point: the dominant error is not in the energy, it is
  contributed by the log10 taken at print time, and that step is the one
  routinely overlooked.  gA is printed in floating form and carries its own
  digits; log gf does not.

  The author's decision was that nothing is to be corrected numerically and
  the LOG GF field is not to be widened -- widening it would move OUTG11
  columns and break the Python and Perl parsers to fix something that is not
  broken.  Instead the argument above is now permanent in three places: the
  header of rcg11kd.f, the top of SUBROUTINE SPECTR which is what writes the
  list, and README.md.  The wording is the same argument in each, so that
  whichever is found first is enough.

TWO DEFECTS THE SURVEY FOUND ON THE WAY

  1. rcg11kd.f, SPECTR: GA=0.6698D0*PSQ*ERAS*(TEMP**3).  A dropped digit --
     the seven other sites all write 0.66698D0.  It is +0.42%, and it drives
     the TEXCIT-weighted line-selection threshold, so it silently changed
     which lines survived the cutoff.  Not itself a printed quantity.
     Fixed, with the reasoning at the line.

  2. RCN36K.F used two different fine-structure constants.  zeta1 took
     2.66246e-5 = alpha^2/2 from a pre-1969 alpha, 4.0e-5 low, while
     sm, sn and vk took 2.921830 = (alpha^2/2)*Rinf evaluated with the
     modern alpha and exact to all eight digits.  Half of it had been
     updated at some point and half had not.  Both now derive from ALPHFS.

  Neither is a formula error.  Every conversion in the package was checked
  algebraically against the Grant and Martin-Wiese forms in
  WORK/A-S-f-conversion_formulas.docx first, and every one of them has the
  right shape.  This was a pure constants job from the beginning.

HOW fundconst.f IS BUILT

  Two sections.  Section A is CODATA transcriptions and nothing else -- the
  numbers as NIST publishes them, each with its relative uncertainty in the
  trailing comment.  Section B is everything derived: unit-scaled variants,
  products, quotients.  A Section B entry is computed from Section A rather
  than transcribed, so that two spellings of one constant cannot drift
  apart, and each carries a "C     LEGACY:" line above it holding the
  literal that used to be in the source, bit for bit.

  Three properties of the file are load-bearing:

  - A PARAMETER may only use names defined above it.  That is why HBARJ is
    in Section B and not Section A: it is PLANCK/(2*PICON) and PICON lives
    in Section B, and putting it in A would have required a second spelling
    of pi in the one file that exists to prevent exactly that.

  - The include is textual, so a PARAMETER whose name collides with a local
    variable in a host routine is an assignment to a parameter.  If the
    local is only ever read, that is a silent wrong answer and not a compile
    error.  Every candidate name is therefore collision-checked against
    FOR/*.F, FOR/*.f and FOR/LSJJ/ before it is written down.  The rule is
    at the head of the file and it earned its place: it forced ALPHFS and
    PICON (a local ALPHA, a local PI), then TWOPIC (subroutine phshift in
    RCN36K.F has a local twopi, and phshift includes the file).

  - The LEGACY lines plus CODE/legacyize.py can put the old numbers back,
    which is what made the two-stage protocol below cheap.  That is no
    longer completely true; see the RINF caveat further down.

THE TWO-STAGE PROTOCOL, AND WHY IT IS NOT NEGOTIABLE

  CLAUDE.md's acceptance standard is bit-exactness, and this work changes
  numbers on purpose, so diff stops being an oracle the moment the values
  move.  Every section was therefore done in two stages, and the ones that
  folded a coarse rounding into a shared name in three:

    stage A   literals replaced by parameters still holding the LEGACY
              values.  Must be BYTE-IDENTICAL to the committed build.
              This is the plumbing check, and it is the whole reason the
              LEGACY lines exist.
    stage B   values updated in fundconst.f only.  Everything that moves
              here is attributable by construction, because nothing else
              changed.
    stage B2  temporary names that existed only to preserve a coarser
              rounding folded into the shared name, now equal by
              definition.  Byte-identical again, which is what proves the
              fold.

  Stage A is worth more than it looks.  Adding an include to a routine that
  never had one is the exact circumstance in which a shadowed local goes
  unnoticed, and a stage A that is byte-identical rules it out.

  Where diff could not be the judge, CODE/cmp_num.py is.  It tokenizes two
  files into a text skeleton plus a numeric vector, requires the skeletons
  to match exactly -- structural change is never acceptable -- and then
  compares the numbers with relative and absolute tolerances.  Its exit
  codes are three-valued and the distinction is the point: 0 within
  tolerance, 1 the numbers moved and the text is intact, 2 something moved
  that was never allowed to.  Four tokenizer decisions in it are worth
  knowing before reading a report, and they are documented in the plan; the
  shortest of them is that E is the only exponent letter and D is not,
  because no output file this package writes has a D exponent while 4d7
  appears everywhere.

  Five regression cases in the end, not the usual one: WORK (Sn VIII),
  WORK/Forbid (the only M1/E2 deck), WORK/AI and WORK/AI/DIEL (which had
  never been used as regression cases before this work and had to be
  baselined first), and WORK/Pr3.  Without the two AI decks the largest
  change made here would have been unexercised.

WHAT MOVED, PROGRAM BY PROGRAM

  RCN (RCN36K.F).  The single-precision hazards are here and were the
  reason this program took longest: several of its constants feed REAL*4
  paths where a longer literal does not survive the assignment.  Its
  spin-orbit, mass-velocity and Darwin coefficients now all come from
  ALPHFS, which closes defect 2 above.  TAPE2N moves, and the author's
  decision was to promote the new TAPE2N as the benchmark from that point
  on rather than carry two.

  RCN2 (RCN2K.F).  Four routines -- ZETA2, DIP, SLI2 and G5INP.  Small and
  uneventful.  The one intended visible difference anywhere in RCN2 is
  OUT2's timing layout.

  RCG (rcg11kd.f).  36 sites in nine routines: SPRIN, CALCFC, MUPOLE,
  ENERGY, SPECTR, WNDIST, BORN, RCOEFF and RCEINP.  Each carries
  INCLUDE 'fundconst.f' after its INCLUDE 'rcgpar.f' -- that include, not an
  IMPLICIT card, is where RCG's IMPLICIT REAL*8 comes from, so the order
  matters.  The value change moves the transition list and nothing else:
  OUTGINE, TAPE2E and genealogy.txt are byte-identical, and in OUTG11 the
  shifts are 1.9e-5 in gf(E1), 8.6e-5 in gA(E1), 5.8e-5 in the M1/E1 ratio
  and 6.7e-4 in the E2/E1 one -- the last because its a0 was 0.529, 3.4e-4
  low, and enters squared.  No line count and no sum-rule structure moved,
  in any case.

  RCE (RCE20KD.F).  Nothing.  This was verified mechanically rather than by
  eye: the literal scanner was run over all six files of the RCE family and
  found not one real literal with three or more significant digits anywhere
  in them.  The most precise reals in the whole of RCE are 0.99D0 (a
  step-damping clamp), 0.001E0 (a divide-by-zero guard), 0.5D0, 2.D0 and
  60.D0 (seconds to minutes).  RCE therefore gets no INCLUDE 'fundconst.f':
  adding one would buy a compile-time dependency for nothing.

  One consequence of that is easy to mistake for a regression, so it is
  said plainly.  RCE's numbers still move even though RCE did not change,
  because it reads TAPE2E from RCG and RCG's constants moved.  Anything that
  shifts in RCE output after this work is inherited from upstream.  That is
  the expected behaviour of a pure least-squares fitter.

THE CONVERSION COEFFICIENTS AS ROUTINES

  Three functions at the end of rcg11kd.f:

    GAFROMSCOEF(ITYPE)   gA = C * S * sigma**(2p+1)
    GFFROMGACOEF()       gf = C * gA / sigma**2
    GFFROMSCOEF(ITYPE)   gf = C * S * sigma**(2p-1)

  ITYPE = 1, 2, 3 for E1, M1, E2; sigma in cm-1, S in atomic units.  They
  are built entirely from fundconst.f and hold no state, so they compile and
  test on their own.  The third is written as the product of the other two,
  so the type factors are derived once in the file rather than twice.
  -fno-automatic makes them non-re-entrant: do not call them from inside a
  parallel region.

  They are also the one independent check on fundconst.f in the whole
  package.  What they return -- 2.02612688449478D-06, 2.69735002205593D-11
  and 1.11995002964537D-22 -- are the Martin and Wiese handbook
  coefficients 2.0261e18, 2.6974e13 and 1.11995e18 (S and lambda in
  angstrom) to every digit that handbook prints.  There is no other oracle
  for these numbers anywhere in this project, which is why
  CODE/tests/test_coefficients.py compiles the functions out of the source
  and checks them against the handbook rather than against themselves.

  SPECTR now takes its PSQ from GFFROMSCOEF instead of deriving the same
  three type factors inline, which was two independent derivations of one
  piece of physics in one file.  Its E2 branch keeps the angstrom scaling it
  has always used -- the coefficient is multiplied by 1D16 at the one point
  PSQ is formed -- because PSQ also feeds the BRNCH, TEXCIT and unit 20
  paths, which do not divide by lambda^2 and whose arithmetic was not to be
  disturbed.  The reassociation costs 0 ulp on E1, 2 on M1 and 4 on E2,
  measured rather than estimated, and all five cases stay byte-identical.
  ABOHR8 went out of fundconst.f with the inline derivation, that having
  been its only site.

THE SWEEP, AND THE ONE LARGE ERROR

  The inventory that drove the first six sections was built by searching for
  the known spellings of the constants, which left anything spelled
  unexpectedly invisible to it.  So afterwards every non-comment real
  literal of three or more significant digits in all four programs was read
  one by one.  That turned up three more sites carrying physics, and the
  first of them was worth more than every value change made up to that
  point.

  The autoionization rate in ENERGY was

    AA(LP,I)=((SUM*UENRGY/109735.3D0)**2)*6.28318D0*2066.D0

  which is Fermi's golden rule, 2*pi*|V|^2 in rydbergs converted to units of
  1e13 /s.  Ry/hbar is 2067.0686667586.  The literal was 2066 -- a
  truncation, not a rounding, and 5.2e-04 low.  It is the largest constant
  error that was left anywhere in the package, and unlike the 0.6698 defect
  it is exercised, by WORK/AI and WORK/AI/DIEL.  The line now reads

    AA(LP,I)=((SUM*UENRGY/RINF)**2)*TWOPIC*AURATE

  and the printed rates rise by 4.81e-04, not the 5.2e-04 the literal alone
  implies, because RINF sits squared in the same denominator and its own
  1.8e-05 pulls the other way:

    (2067.0686667586/2066) * (109735.3/109737.31568157)**2
                           * (6.28318530717959/6.28318) = 1.000481354

  which is what the runs show, off SUMGAA: 1.96195E+06 -> 1.96290E+06,
  6.48954E+06 -> 6.49266E+06, 9.05517E+06 -> 9.05952E+06 in AI/DIEL, and
  1.58222E+15 -> 1.58298E+15 in AI.  AI/DIEL's GM*FRBAR and its whole
  T(KEV) dielectronic-recombination table move with them.

  This is a visible change to printed output and the largest this work
  makes.  Autoionization rates printed by every earlier version of RCG were
  0.05 percent low.  It is said here in plain words rather than left in a
  table of constants, because anyone comparing against old results will see
  it.

  The second site was CSIGMA in SPECTR, 8.06662D-18 for 4*pi^2*alpha*a0^2 in
  cm^2, 8.2e-05 low, now PICSCN.  It gates the photoionization
  cross-section path, which no deck in WORK/ turns on, so it was checked
  algebraically and by nothing else.  That is stated rather than glossed:
  an unexercised site must not be reported as a verified one.  The third
  site was 109.737 in RCN2K.F's G5INP, 2.9e-06 low.

  base -> stageB2 over the five cases differs in two files only: AI/OUTG11
  and AI/DIEL/OUTG11.  TAPE2E, OUTGINE, ING11, TAPE2N and genealogy.txt are
  byte-identical in every comparison, including that one.

  THE RINF CAVEAT.  RINF now stands for three distinct legacy roundings --
  109737.31, 109737.0 and 109735.3 -- and its single LEGACY line can restore
  only the first.  legacyize.py alone can therefore no longer reproduce the
  package's original text.  A bit-exact reversion has to give the odd sites
  temporary names of their own, which is what the staged work did (RINFL3 =
  109735.3D0, RINFK7 = 109.737D0).  The caveat is recorded under RINF in
  fundconst.f itself, not only here.

THE LAST FOUR, WHICH IS WHERE THE RULE GOT WRITTEN DOWN

  The sweep also found four literals that are physical constants and were
  still inline: the electron g-factor spelled 0.50116D0, 1.50116D0 and
  1.00232D0 in the Lande g-factor and the magnetic-moment matrix element,
  and 2.1717D-8, the Maxwellian rate coefficient in RCOEFF.  All four are
  at least an order of magnitude smaller than the 2066 error.  They were
  held in a group of their own in test_literals.py, PHYSICS_STILL_INLINE,
  which was a holding pen with an expiry date rather than a category.

  The author's instruction emptied it: every constant, no threshold.  They
  are now GEFACT in Section A, with the sign convention said out loud
  because the code uses |g| throughout and a sign flipped there would be
  silent, and GEM1H, GEP1H, GEM1 and CMAXWL in Section B, each derived from
  it -- CMAXWL from PICON, ALPHFS, CLIGHT and ABOHR as
  2*sqrt(pi)*alpha*c*a0^2.  DSQRT inside a PARAMETER is a gfortran
  extension; it was verified to compile and to fold before CMAXWL was
  written that way.  SPRIN, CALCFC, MUPOLE and RCEINP gained the include;
  each was checked against the whole name list of the file, not only the
  four new names, and there were no shadowed locals.  PHYSICS_STILL_INLINE
  is deleted -- the group, its comment block, the docstring paragraph and
  its entry in the inventory loop -- and test_literals.py now asserts that
  FOR/ contains no physical constant outside fundconst.f, with no exception
  list at all.

  RCN and RCN2 were not re-run for this and did not need to be.  Both
  include fundconst.f, so the change could in principle have reached them.
  It does not: RCN36K.F and RCN2K.F were compiled against the old and the
  new fundconst.f with the production flags and the object files are
  bit-identical.  That is a stronger statement than a matching run.

  The g-factor move is the second of the two places in this whole revision
  where a number the author reads directly changes, and it is small.  In
  five regression cases exactly one printed Lande g-factor moves: Sn VIII,
  J=1.5, the fourth level, 2.100 -> 2.099.  That was verified rather than
  assumed -- a throwaway build with F20.14 in place of F9.3 reads it as
  2.09950010910943 -> 2.09949934602454, which crosses 2.0995.  Forbid's
  MUPOLE sum rules SUMS21 and SUMS20 fall by 2.3e-07 and 1.4e-06, the
  latter bounded by twice the former because the M1 element enters squared.
  No gf, no gA, no M1 or E2 matrix element and no line count moves anywhere.

  TAPE2E differs in all five cases, in one thing only: the g-matrix RCEINP
  packs onto unit 15.  An instrumented build dumping that matrix at E26.18
  measures the whole of it -- max |delta| 1.62e-06 (Sn VIII, Pr III),
  1.16e-06 (Forbid), 6.96e-07 (AI), 1.39e-07 (AI/DIEL), on elements of
  magnitude at most 3.34.  OUTGINE is the same information printed and is
  byte-identical everywhere, so nothing RCE reads at printed precision
  moves.

  Three sites are not exercised by any deck in WORK/, and each was
  established positively with an instrumented build that produced no output
  at all, rather than inferred from the absence of a diff: CMAXWL, because
  RCOEFF is never reached; TAPE19, because CALCFC's WRITE (19,91) sits on a
  path none of the five decks takes; and CALCFC's P5 g-factors, which are
  under IF (KCPL.LT.2) while every deck here runs KCPL=2.  CMAXWL was
  checked algebraically instead: 2*sqrt(pi)*alpha*c*a0^2 is
  2.1716669877207e-08 against the 2.1717e-8 the source carried.

THE PYTHON SIDE

  CODE/toolbox/radiate.py carried its own table of the same coefficients,
  transcribed from conv_out.pl, consumed by conv_out.py and calc_TP_M1E2.py,
  and unreachable by a change to fundconst.f.  It is now a derived copy.

  The mechanism is a generated CODE/toolbox/fundconst.py, not a parse of the
  Fortran at import time.  Parsing would be truer to "one source" and worse
  engineering: CODE/ is what gets copied to the directory on the path,
  FOR/ is not distributed with it, and a module that fails to import when
  the source tree is absent is a bad trade for a value that changes once a
  decade.  What holds the copy honest is CODE/tests/test_fundconst.py, which
  compiles FOR/fundconst.f into a driver printing every parameter it
  declares at 1PD26.17 and asserts each equals its Python counterpart with
  ==, not a tolerance.  All 42 agree bit for bit.  Two structural tests go
  with the values: a Fortran name with no Python counterpart fails, and a
  Python name that is not in the Fortran fails, so the copy cannot drift as
  the file grows and cannot quietly become a second source of truth.  It
  skips without gfortran, so the toolbox stays importable anywhere while the
  tree stays provably consistent wherever it can be checked.

  radiate.py's E1, M1 and E2 are now spelled exactly as GAFROMSCOEF spells
  them and equal the compiled functions to the last bit, asserted in both
  directions.  The spelling matters: the general Grant factor gives the E2
  type factor as 5.52754142998615672D-17 against the Fortran's ...795D-17 --
  the same number algebraically and one bit apart.  M2, E3 and M3 have no
  Fortran counterpart and take the general Grant factor

    Np = 2^(2p+2) * (p+1) * (2p+1) * pi^(2p+1) / (p * [(2p+1)!!]^2)

  but only through the ratio Np(p)/Np(1), so the reconstruction's prefactor
  divides out and the absolute scale still comes from GAFROMSCOEF.  For p=2
  that ratio is pi^2/5 exactly in floating point, and a test asserts it.

  Two defects came out of doing this.

  1. THE LEGACY gf(E2) IS WRONG, by +3.6e-04, and not because it is stale.
     It is an internal inconsistency in the legacy table: gf is
     gA*lambda^2/(GACON*1e16), which is exact and type-independent, so
     1.11995e18 and 167.84224 cannot both be right, and the Grant factor
     confirms the gA one to 2.65e-08.  167.90221 is the value consistent
     with the rest of the table.  Found the same way as the 0.6698 defect --
     by deriving what had been copied.  It reaches the E2 rows of
     conv_out.py's absorption model, which move by 3.6e-04 against its E1
     and M1 rows.  Pinned by a test so the legacy number cannot later be
     mistaken for a target.

  2. EV_PER_CM WAS OUT BY A FACTOR OF A THOUSAND.  It held 8.06554393734921
     -- eV per 1000 cm-1 -- and both call sites in conv_out.py divided by it
     a quantity the *1000 had already converted to cm-1.  The Auger dE(eV)
     column has been printing meV under an eV heading since it was written.
     The Perl does the same (conv_out.pl:733, 965, 968) and its comment
     names the constant "eV to cm-1", which is 8065.54: the value was the
     mistake and the usage was right, in both languages.  WORK/AI/DIEL's own
     README confirms it independently -- its free electron is 18.58 Ry =
     2038.92 kK = 252.8 eV, and the old expression gave 252 800.  Put to the
     author 2026-08-21, who chose to fix it.  Both sites now take
     fundconst.EV2CM directly and the column reads in eV.  radiate.EV_PER_CM
     is deleted rather than corrected in place, so a caller written against
     the old behaviour fails loudly instead of quietly changing by a factor
     of 1000.  This is a deliberate, documented divergence from
     CODE/perl_legacy/.

  No test expectation needed changing, and that was checked rather than
  assumed: every expectation in test_conv_out.py and test_calc_tp.py that
  goes through radiate is a round trip, a ratio of two quantities converted
  the same way, or a comparison of two rows of one table -- all insensitive
  to the constant's value.  Not one absolute S, gf or gA figure is written
  down in either file.  A test expectation edited to match new output would
  not have been evidence of anything.

  CODE/perl_legacy/ is untouched.  It is the reference for intent and is not
  maintained; radiate.py says in a comment that its constants are now
  knowingly stale, so the next reader does not "restore" them.

WHAT AN OLD RESULT WILL LOOK LIKE AFTER THIS

  Everything below is a value change and none of it is a bug being fixed,
  except where it says so.

    autoionization rates (WORK/AI, AI/DIEL)     +4.81e-04   the 2066 fix
    gf(E1) in OUTG11                            +1.9e-05
    gA(E1) in OUTG11                            +8.6e-05
    M1/E1 ratio                                 +5.8e-05
    E2/E1 ratio                                 +6.7e-04    a0 was 0.529
    conv_out.py E2 absorption rows              +3.6e-04    a real defect
    conv_out.py Auger dE column                 x1000       a real defect
    Lande g-factors                             -2.3e-07    one field moves
    M1 matrix elements                          -6.9e-07
    Maxwellian rate coefficient                 -1.5e-05    unexercised
    photoionization cross-section coefficient   +8.2e-05    unexercised

  Structurally nothing moved anywhere: no line count, no label, no field
  width, in any of the five cases, at any stage.  Full suite at the end of
  the work: 309 passed, 252 skipped.


LSJJ.DAT IS BUILT BY make_cfp, AND THE SIGN CONVENTION IS FIXED
===============================================================

Two changes, done together because the second one moves the bytes of the
deck and the first one is what makes that harmless.

WHY THE DECK STOPPED BEING DISTRIBUTED

  LSJJ.DAT was committed in two places -- CODE/LSJJ.DAT and
  FOR/LSJJ/LSJJ.DAT -- and it never should have been.  It is a derived
  file of exactly the same kind as FOR072, FOR073 and FOR074: unformatted
  sequential, native INTEGER*4 and REAL*8 in native byte order, wrapped in
  the record markers of the compiler that wrote it.  The CFP decks have
  been built on the installed machine since the beginning, by make_cfp,
  and shipping the fourth deck while building the other three was an
  inconsistency, not a decision.

  So make_cfp builds all four now, on both platforms.  It costs 0.2 s
  against the RCG run make_cfp already does, and it needs nothing that is
  not already in the working directory at that moment: the ING11 that
  make_cfp itself has just written from CODE/ING11.CFP, and SENIOR.
  cowan_cfp.py takes both from $LSJJ_CFP and $LSJJ_SENIOR when those are
  set, which is why no file has to be moved to make it work.

  The generator moved to CODE/lsjj/ for this.  make_cfp has to find it
  from beside itself in an *installed* tree, and FOR/LSJJ/ is not
  installed: CODE/unix/install.sh copies CODE/ and nothing else, and on
  Windows the author copies CODE\ to C:\cowan\code by hand.  Four modules
  moved -- lsjj.py, cowan_cfp.py, recur.py, mkdeck.py -- and the scripts
  left here (validate.py, cmpfort.py, compare_racah.py, probe_f3.py) each
  gained a two-line sys.path insert pointing across at the one copy.  One
  source of truth, no duplicated module.

numpy

  This is the one place where the package's "nothing to pip install"
  promise now has an exception, and it is worth saying why it could not be
  avoided.  lsjj.py uses np.linalg.svd for the null space of the coupling
  operator and np.linalg.eigh to block a degenerate space by a commuting
  operator.  Both return a *basis*, not a number.  Any other
  implementation returns a different, equally valid basis -- and the deck
  would then no longer be the deck that was validated against Cowan's own
  coefficients by validate.py and cmpfort.py.  There is no way to check
  that, which makes rewriting the algebra a section 1 violation with
  nothing to accept it against.

  So numpy stays, and instead the absence of it degrades gracefully.
  make_cfp attempts the rebuild, and where numpy (or Python) is missing it
  prints a specific warning naming `pip install numpy`, builds the three
  CFP decks, and exits 0.  install.sh repeats the warning in its closing
  summary, since that is the text a user actually reads.  RCG already
  tolerated a missing LSJJ.DAT -- unit 23 is opened with ERR= and read
  only by a run that asks for true jj coupling -- so the cost of not
  having numpy is exactly one feature.

THE SIGN CONVENTION

  The deck was reproducible on this machine and not across machines.
  Built under WSL2 (numpy 2.5.2, its own LAPACK) it came out the same
  463628 bytes with 55 of the ~58000 doubles different: 24 a pure overall
  sign flip of a basis vector, 31 last-bit rounding, none a genuinely
  different value.

  The 24 were a defect, and a small one to state.  fix_sign() fixed the
  overall phase of each highest-weight state by making its largest
  component positive:

      key = max(state, key=lambda k: (abs(state[k]), -k))

  These components are square roots of small rationals, so two of them are
  often equal in exact arithmetic.  A strict max() then settles the tie on
  whichever one LAPACK returned a last bit larger -- and a different LAPACK
  crowns the other, flipping the whole vector, and with it every state
  derived from it down the lowering chain, which deliberately fixes no
  signs of its own.  A convention hinging on 1e-16.

  The fix is a tolerant tie-break: compare magnitudes only to
  SIGN_TOL = 1e-12 and settle the tie on the basis key, which is ours and
  is fixed.  1e-12 sits four orders above the last-bit noise on normalised
  components and many orders below the smallest genuine gap between two
  distinct component magnitudes anywhere in the deck.  With an exact tie it
  picks the component the old code did, so nothing changes except at a
  near-tie.

  Result: no sign flips left between the two machines.  Fourteen vectors
  changed sign against the deck that used to be committed -- those were the
  ties this machine had been losing.  JJCFP is written from the same
  convention and was rewritten with `recur.py deck JJCFP`; eight of its 90
  cards changed sign, and jjcfp.log with them.

  WHAT IT DOES NOT FIX, AND WHY THERE IS NO md5.  The other 31 doubles are
  raw LAPACK arithmetic -- svd and eigh returning values that differ in the
  last bit -- and no tie-break reaches them.  Measured after the fix,
  Windows against WSL2: 28 doubles differ across 21 records, |delta| <=
  2.3e-16, no sign flips and no genuinely different values.  So the deck is
  still not byte-portable, and a hard-coded md5 would be a tripwire that
  fires on every machine but the one that recorded it.
  CODE/tests/test_lsjjdeck.py says so at length and checks what is
  checkable instead: that two builds from the same inputs are byte
  identical, and that the deck make_cfp leaves behind is the one mkdeck.py
  builds today.

ACCEPTANCE

  The deck changed, so this is a section 1 change and was run as one.

    validate.py all, validate.py all cowan, compare_racah.py all
        identical to the committed logs, line for line.

    recur.py fromdeck check all, after rewriting JJCFP
        identical to recur_all.log except in the two timing columns --
        every residual, to the printed digit.  Before JJCFP was rewritten
        it failed f^4 through f^10 with a deviation of order 1, which is
        the phase mismatch showing itself and is how the stale deck was
        caught.

    KCPL=9 self test on WORK/W70's ING11, then cmpfort.py
        34 shells, 34 pass.  Agreement with recur.py *improved* by an
        order of magnitude, from 4e-14 to 3e-15; the orthogonality and
        1-vs-2 residuals did not move.  cmpfort_all.log rewritten.

    WORK/W70, the true jj case, old deck against new
        TAPE2E byte identical -- and that is where the true jj trailer
        lives, so the transformation RCE reads is unchanged.  OUTG11
        differs in one line, the printed timing.  true-jj_keys.txt,
        genealogy.txt, OUTGINE identical.

    WORK/Pr3, the regression case
        TAPE2E, OUTGINE and genealogy.txt byte identical; OUTG11 differing
        only in timing lines.  This also covers the rcg11kd.f change
        below, the executables having been rebuilt for it.

  rcg11kd.f changed in five places, all of them text: the runtime
  diagnostic at label 65 and the comment blocks that used to say
  `cd LSJJ ; python mkdeck.py ../LSJJ.DAT` now name make_cfp, which is the
  route an installed user has, and the paths that said LSJJ/ now say
  CODE/lsjj/ or FOR/LSJJ/ as the case may be.

  Full suite from CODE/: 312 passed, 252 skipped.


A VERSION RESOURCE ON THE WINDOWS EXECUTABLES
=============================================

Not part of the LS->jj work at all; recorded here because this file is
where the account of everything else is kept.

THE PROBLEM

  Norton blocks the download of the package: the source zip from
  codeload.github.com, a Release asset from objects.githubusercontent.com,
  and RCN2K.exe on its own, all with Win64:Evo-gen [Trj].  Since the
  Release asset is blocked as well, the verdict is not attached to the
  host -- Norton unpacks the archive, scans the PE files inside, and reacts
  to their content.  Win64:Evo-gen is a heuristic family, not a signature
  match.

  The author ruled out false-positive submissions and code signing, both
  free and paid, for the same reason: each is a chore at every release.
  What is left has to be automatic.

WHICH FEATURE, MEASURED

  All four executables were built within a minute of each other, so hash
  novelty and prevalence are excluded: three download and one does not.
  Variants of RCN2K.exe were built with the shipped flags and copied to
  Downloads to be scanned:

    plain                            6 trials, 0 survived
    + VERSIONINFO (0x398 bytes)      4 trials, 4 survived
    + a string-table resource only   2 trials, 2 survived
    + 1 KB of data, no resource      2 trials, 0 survived
    + -s (strip)                     1 trial,  1 survived

  The two controls separate the cause cleanly.  A resource carrying no
  version information at all survives, so the *content* of the resource is
  not being read; an equal-sized padding of ordinary data with no .rsrc is
  deleted, so it is not mere perturbation.  The discriminating feature is
  the presence of a resource directory.  A PE without one is unusual for
  legitimate Windows software and common in malware.

  Why RCN2K alone, when RCN36K has no resource either, is not explained by
  this and does not need to be: a heuristic score combines features, and
  RCN2K is the smallest of the four at 176 KB.

WHAT WAS BUILT

  CODE/version.txt     the identity strings, key=value.  Nothing bumps it
                       automatically and nothing needs it bumped -- the
                       resource's purpose is to exist, not to be accurate
                       to the day.  It is also the first and only place
                       the package version, 26.2.1, is written down.

  CODE/mkverres.bat    reads version.txt, maps the source stem to a
                       FileDescription and an OriginalFilename, writes a
                       .rc into %TEMP% and compiles it with windres, which
                       ships with the same MinGW as gfortran.  Any failure
                       -- no windres, no version.txt, an unknown stem --
                       prints one line and exits 0 with no .o.

  gfortran_O3_OMP.bat and gfortran_O3.bat call it and link %1_ver.o if it
  is there.  gfortran_dbg.bat does not: a debug build is never
  distributed.  All four programs get the resource, not just RCN2K -- it
  is the same one line in the build, the information is true of each
  program, and it removes the feature before a future edit pushes another
  binary over the same threshold.

  ELF has no VERSIONINFO, so CODE/unix/build.sh is untouched.  This is one
  of the few places where CLAUDE.md section 8's "Windows and Unix move
  together" does not apply, because there is no other half to move.

ACCEPTANCE

  Each executable grew by exactly 1042 bytes and now shows its strings in
  Properties -> Details.  WORK/Pr3 with the rebuilt executables: TAPE2E,
  OUTGINE and genealogy.txt byte identical, OUTG11 differing only in
  timing lines.  Full suite from CODE/: 312 passed, 252 skipped.

  WHAT IS NOT YET KNOWN.  Surviving an on-access scan is a strong proxy
  but it is not the download verdict, which is a different code path with
  a cloud lookup and can disagree.  The conclusive test is cutting a
  release and downloading the zip, and only the author can run it.  If it
  still fails, -s (strip) is the one cheap thing left to try -- it cleared
  the local detection once, and plausibly acts through a different feature
  -- and it was deliberately held back so that a failed download test
  stays attributable to one change.

  A note that does not change with any of this: the file is still unsigned
  and still has near-zero prevalence.  What the measurement says is that
  those are not the operative terms.


THE ONE-STATE J BLOCK AND THE TRUE jj SCRATCH FILE
==================================================

Reported by the author on 2026-08-21, from the Sn VIII installation test
in WORK/ with true jj switched on -- the same IN2 switches W70 uses,
KCPLD(8)=8 in column 38.  `printout` stopped with

    toolbox.levels.LevelsError: OUTG11 has no J=9/2 block in
    configuration set 2, but RCEOUT lists 26 levels for it

which is `identify` saying that RCE has levels RCG never printed.  The
same sequence in WORK/W70 had always worked.

WHAT WAS ACTUALLY WRONG

  Not the Python.  RCG had died, and the message is what a truncated
  OUTG11 looks like from downstream.  The run ends with

    DONE ENERGY FOR  J= 2.5  AT TIME=  0.00 MIN,  MATRIX SIZE=   39
    UNPK1 error: wrong record in file          42
    BB1(5),BB1(6),IM1,IM2,... =  37  8  36  36  20000  20000  1296

  and above that, three of

    *****TRUE jj RECORD OUT OF STEP: GOT J=  5.5, NLS=     1;
     EXPECTED J=  0.5, NLS=    21.  NO TRUE jj OUTPUT FOR THIS J.

  Unit 42 is the true jj scratch file: CALCFC writes one record pair per
  J block -- the labels, then the PK1-packed transformation matrix --
  and SPRIN reads one pair per block, in the same order.  The stream had
  slipped by exactly one record, and the messages track it: the second
  configuration set read the first set's last matrix, then its own
  first, then its own second, and at the fourth block a packed row index
  of 37 arrived in a 36 x 36 matrix and UNPK1 stopped the run.

THE MECHANISM

  SPRIN prints a block twice, once in the coupling the run asked for and
  once in the other, and the true jj pass is entered from the second of
  the two.  A block of ONE basis state has no second pass: there is
  nothing to transform, so at label 289

    289 IF ((NLS-1).LE.0) GO TO 400

  took it straight out of the routine -- past the read.  Its record
  stayed on the file and every later block took the record of the block
  before it.

  Sn VIII has such a block: J=11/2 of the first set, the single
  4d7 (2H) 2H state.  It is the sixth of fourteen, so eight blocks
  followed it.

WHY W70 NEVER SHOWED IT

  W70 has one too -- J=11/2 of its second set, the single sp4f state --
  and it is the LAST block of the run, so nothing followed it to read
  the record it left behind.  The only symptom there was the
  absence of its own true jj section, which nobody had reason to look
  for: `printout` fills that level from genealogy.txt (the fallback in
  toolbox/case.py, written for Cowan's second coupling, which really
  does print nothing for a one-state block) and its level counts come
  out right either way.  The reference case for the feature could not
  have found this.

THE FIX

  Three changes in rcg11kd.f, all inside the KPAR=4/IJJROT path:

  - Label 289 no longer leaves.  A one-state block enters the true jj
    pass, which consumes the record and prints the section.  VECT(1,1)
    is set from C(1,1) first, because the product at 287 reads the LS
    vectors from VECT when KCPL=1 and from C otherwise, and with no
    second pass to have moved them they are in C either way.

  - The other early exit, the term-keep card (NOLSKP), consumes its
    record and discards it.  There is no true jj output to make from a
    truncated LS basis -- LSJJRT will have refused the block and IOK=0
    says so -- but the record still has to come off the file.

  - The read is now SUBROUTINE LSJJRD, called from both places.  It
    unpacks the matrix at the dimension the RECORD carries rather than
    the caller's, so a reader out of step gets the existing OUT OF STEP
    message instead of stopping the run inside UNPK1 -- which is what
    turned a missing section into a dead run.

  RCEINP, which copies the same records onto the tail of TAPE2E, was
  never affected: it rewinds and reads them all, so RCE always had the
  one-state block's matrix.  That is the second reason for printing the
  section rather than discarding the record: with it, OUTG11 says what
  TAPE2E carries.

  Printing it is exact.  A one-dimensional block has one eigenvector,
  the rotation is a 1 x 1 orthogonal matrix, and all it can do is change
  a sign; the content of the section is the label, which is what a jj
  section is for.

  On the Python side only comments moved, plus one test:
  CODE/tests/test_outg11.py::test_a_one_by_one_block_has_no_jj_section
  now expects LS and TRUEjj and no JJ.  The genealogy fallback in
  toolbox/case.py and toolbox/iden.py stays as it is -- it is still
  needed for Cowan's second coupling, and for any OUTG11 written before
  this fix.

ACCEPTANCE

  WORK/Pr3, both stages of its README -- the Hartree-Fock run and the
  fitted run after update11 -- TAPE2E and OUTGINE byte identical against
  the previous executable, OUTG11 differing only in timing lines.  The
  publication chain that follows reproduces every number the README
  quotes: 5626 levels, 16878 in three schemes, 39 configurations
  suggested 19/13/7, 1036 published, 45 letter prefixes.

  WORK/W70: TAPE2E and OUTGINE byte identical, and OUTG11 differs by
  exactly the seven added lines of the true jj section of J=11/2 --
  nothing deleted, nothing changed.  levels.txt from `printout` is
  byte identical.

  WORK/ (Sn VIII) without true jj: TAPE2E and OUTGINE byte identical,
  OUTG11 identical modulo timing.  With KCPLD(8)=8 it now runs to
  NORMAL EXIT, reports 14 J blocks appended to TAPE2E, and the chain
  that failed -- RCE, then `printout ... --couplings LS,JJ,TRUEjj` --
  gives 14 J blocks, 199 levels, 597 rows.

  Full suite from CODE/: 450 passed, 107 skipped, with the reference
  runs of Sn VIII, W70 and Pr III (JK folder included) present.  The
  seven that fail do so with the previous executable as well and are
  not about this change: five want a W70 run whose energies the
  committed IN2/IN36 do not reproduce, and two want the author's own
  Pr3 ING11 -- the one with Cowan's two prelude cards, already once
  through s11.pl, whose zero parameters therefore read `  0.00e+01`
  where a freshly generated deck has RCN2's `0.00E+001`.  Both are
  worth a look, and neither belongs to the LSJJRD fix.
UPDATE11 NO LONGER NEEDS THE PARITY COUNT TYPED AFTER IT
========================================================

The trailing `2` of `update11 2` is Perl-era: s11.pl had no way of knowing
how many parity sets the run has, so the user told it.  The Python does
know -- toolbox/in36.py parses IN36 the way RCN's own ANALYZ and SETCFG do
-- so the count can be read off the deck instead.

WHAT WAS ADDED

  s11.py grows one function and one command word.  `s11 auto` calls
  parity_sets(folder), which is len(in36.parities(...)) over IN36's
  configurations, and substitutes both parity sets when there are two
  parities and the first only when there is one.  It says which it chose
  before it starts, so the choice is in the log:

      IN36: both parities; substituting both parity sets.

  IN36 cannot always answer.  A sub-folder holding only ING11 and RCEOUT,
  as WORK/U90/reordered does, has no IN36 at all, and WORK/AI/DIEL's
  shipped IN36 is the bare end card until `diel` writes the real one.
  Both give 0, and the fallback is to walk every set: main() already skips
  a parity set the deck has no configurations in, so a one-parity deck
  comes out the same either way.

  `s11` with no argument and `s11 2` are untouched -- the first still does
  the first set only, the second still does both -- so nothing that was
  typed before changes meaning.

  update11.bat and CODE/unix/update11 pass `auto` when the user typed
  nothing, and pass through whatever the user did type.  `update11 2`
  therefore still works and still means both.

ACCEPTANCE

  WORK/Bi56, a real fitted run of both parities: `s11 auto` beside its
  IN36 and `s11 2` without it produce byte-identical ING11 decks, both
  reporting 1 parameter into set 1 and 2416 into set 2.  Counting alone,
  over the shipped decks: WORK/, U90, Pr3 and W70 give two parities, AI
  and Forbid one, AI/DIEL none until diel has run.

  CODE/tests/test_rewriters.py gains three tests -- the count over every
  shipped deck, `auto` against explicit `2` on U90, and the no-IN36
  fallback.  Full suite from CODE/ unchanged otherwise.

  READMEs: WORK/Pr3, WORK/U90 and README.md now say `update11`, and
  README.md's s11 entry lists `auto`.


PRINTOUT IN A FOLDER RCE HAS NOT RUN IN
=======================================

The JK/ arrangement of WORK/Pr3 was written to be *read from the folder
above*: `printout --from LSJK=JK` merges the second coupling into the
parent's level list, and the sub-folder itself was never a place one ran
printout.  Running it there is the obvious thing to try -- the run is in
that folder -- and three separate things stopped it.

  IN36.  RCG and RCE do not read it; it is RCN's input, so JK/ has none,
  and Case.read stopped with `IN36 not found in .`.  The `parent`
  argument already existed for exactly this, but only the caller that
  passes `--from` was setting it.  Case now falls back to the directory
  above when the run's own folder has no IN36 (`case._in36_folder`); an
  explicit `parent` still wins, and if neither has one the error names
  the run's folder, which is where a reader looks first.

  NLSMAX.  Columns 23-25 became a *permission* in 2026 -- blank reads as
  999, "no limit", and only a punched 0 forbids couplings 3-7
  (rcg11kd.f:897-910 and :955-978).  ing11.py had not followed: it kept
  Cowan's rule that a blank field forbids them, so a deck that asked for
  LSJK with columns 23-25 left alone was reported as asking for LS
  alone, and printout refused the coupling RCG had just written.
  MainCard now holds NLSMAX as RCG holds it, blank turned into 999 on
  the way in, with `nlspch` keeping what was actually punched for the
  one place that has to tell the two apart.

  RCEOUT.  `--no-rce` was there for a run with no fit, but it wrote no
  params.txt at all, and without it the missing RCEOUT came out as a
  traceback.  It became advice first and then nothing at all: a folder
  RCG has run in and RCE has not is a whole case of its own -- the
  second coupling of a fit, most often -- so printout assumes --no-rce
  where there is no RCEOUT to read, and says on the first line that it
  has.  Asking for the option explicitly is now only a way of skipping
  OUTE and RCEINP.HF where RCE *has* run.

THE PARAMETERS OF AN UNFITTED RUN

  params.txt needs three things: the parameter values, their names, and
  the Hartree-Fock values the LSF/HF column divides by.  ING11 has the
  values but not the names -- its cards carry numbers in slot order and
  nothing else -- so the names have to come from somewhere.  They come
  from OUTGINE, RCG's input deck *for* RCE, which lists every parameter
  under the positional name RCE will use (`(7A10)` at :15739) and its
  value (`(7A10)` at :15843) and its least-squares flag (`(7I10)` at
  :15833).  Those three records are the last three blocks before the
  convergence card, so `rce.read_outgine` counts back from the section
  terminator rather than walking the eigenvalues in between, whose
  record count depends on the J blocks and on which of RCG's three
  energy FORMATs the run chose.  The all-zero record between the flags
  and the values is checked, and is what says the count landed right.

  One thing OUTGINE does not have is RCEOUT's columns 58-71, the two
  configurations of a CI integral.  They are recoverable: RCG encodes
  both serial numbers in the first two characters of the name --
  '0'-'9', then 'A'-'Z', then one character per number up to 188, with
  '*' past that (:6308-6322) -- and `parnames._ci_codes` decodes them
  against the parity set's card order.  Where RCEOUT's columns are
  there they are still what is used; the decode is the fallback, and a
  '*' or a number the set does not have decodes to nothing rather than
  to the wrong configuration.

  The HF column is then the deck's value divided by the HR scale factor
  of its card -- the same `_factor` division the fitted path makes, so
  the LSF/HF column comes out as the scale factor itself.  Two things
  about it are not the real Hartree-Fock value, and printout says both
  at the end of the run rather than leaving them to be discovered.  Eav
  carries no scale factor, so its HF column is the deck's own number
  repeated.  And where the deck has already been through s11 -- which
  is exactly the case in WORK/Pr3/JK -- the values are a fit's while
  the factors are still the original ones, so the division recovers
  nothing.  Only RCEINP.HF, written by an RCE run over the unfitted
  deck, has those numbers, and that is what the message names.

ACCEPTANCE

  In WORK/Pr3/JK, with RCG run and RCE not:

      printout levels.txt params.txt --couplings LS,LSJK --no-rce

      LS, LSJK: 25 J blocks, 5626 levels
      wrote levels.txt (11252 levels)
         LSJK J=12.5 of set 2 has one basis state, so RCG printed no
         LSJK section; it is named in LS
      wrote params.txt (1207 parameters)

  11252 is 2 x 5626, both couplings over every level; the one-basis-state
  message is the same one the parent folder's three-scheme run gives and
  is not a fault.  The 1207 parameters reproduce ING11's values exactly,
  card for card, and their LSF/HF column is 0.85 or 1.00 throughout,
  which is what `hr85998585` on those cards says it should be.

  Without the option, in the same folder, the run is the same one with a
  line in front of it:

      printout levels.txt params.txt --couplings LS,LSJK

      no RCEOUT here, so reporting the calculated levels alone, as
      --no-rce does
      ...
      wrote params.txt (1207 parameters)
         --no-rce, so the values are the deck's own and the HF column is
         each one divided by the scale factor of its card; Eav has no
         scale factor and is repeated as it stands.  RCEINP.HF, from an
         RCE run over the unfitted deck, is what makes every HF value
         right

  CODE/tests/ gains eleven tests: one in test_ing11.py over both the
  blank and the punched-0 NLSMAX (and one assertion there changed, the
  W70 card's blank field now reading 999), four in test_rce.py for
  read_outgine (the count-back, the misaligned section it refuses, a folder with no
  OUTGINE, and the ING11 values), three in test_parnames.py for the CI
  decode, and three in test_printout.py for the IN36 fallback, the
  --no-rce run and the same run with no option at all.  One helper
  changed: test_rewriters.copy skips instead of failing when a folder
  that has been through RCG has not been through RCE, which is what
  WORK/Pr3/JK now is.  Full suite from CODE/: 308 passed, 277 skipped.


AN EXPERIMENTAL LEVEL LIST THAT IS NOT AN RCE DECK
==================================================

subst_exp_levs has always had two kinds of source.  One is another RCE
deck, and that is the case it is good at: both sides speak the same
eigenvectors, so a level is recognised by what it is made of and the
match survives the fit being rebuilt.  The other was "a LOPT level
list", and it was a fiction.  It was ported from subst_exp_levs.pl line
for line, it had never been run against data, and every assumption in it
was a hard-coded field number of one output of one version of one
program: ten tab-separated columns, the designation in field 0, the
second uncertainty in field 4, the configuration and term in fields 8
and 9, and the parity taken from an asterisk on the designation with the
first row deciding which parity set was which.  The LOPT the author runs
now writes different columns; and the observed energies do not
necessarily come from LOPT at all.  They come off a spreadsheet, in
whatever shape the bookkeeping of that spectrum happens to have.

So the branch is gone and CODE/toolbox/altlevsrc.py is in its place: a
reader for *any* tab-delimited level list, told by option which column
is which.  Nothing in it is specific to a program.

WHAT IS RECOGNISED, AND HOW

  The old test asked whether the file had exactly ten tab fields.  The
  question worth asking is the other way round: an RCE deck is a file
  rce.read_rceout can read, and everything else is a level list to be
  described.  subst_exp_levs.is_deck is that one call.

  Tab-delimited is required, and quoted fields are unquoted, because
  these files come out of Excel and Excel quotes a field that holds a
  comma.  A header line is found without being told about -- the leading
  rows that hold no number at all are the header -- and once it has been
  found, a column may be named rather than numbered.  --skip forces it.

  --ecol, --jcol and --pcol are required for a source that is not a
  deck; the rest are optional and are how much of the composition the
  file happens to carry:

      --confcol  --termcol  --perccol   the first component
      --conf1col --term1col             its labels, where conf/term are
                                        bookkeeping and these are real
      --conf2col --term2col --perc2col  the second
      --conf3col --term3col --perc3col  the third
      --altcol --altvalue               the second triad is the first
                                        component in another coupling
      --excludecol                      leave this level starred
      --units cm-1|kK|eV  --skip  --parity-sets

  Energies are cm-1 by default, and kK and eV are converted with
  fundconst.EV2CM -- section 4 of CLAUDE.md, no constant written at the
  point of use.

  There is no uncertainty option any more.  Reading LOPT's second
  uncertainty existed so that a poorly determined level could be left
  out of the deck, and a column of uncertainties is one way a file might
  say that and not the general one.  --excludecol says it directly:
  anything non-blank and non-zero in that column excludes the level, a
  column called `exclude` is found without being named, and if there is
  no such column nothing is excluded and nothing is said about it.

THE THREE THINGS THAT MADE IT HARD

  1. Percentages.  A source sheet quotes intensity per cent, unsigned,
  the way a published table does; RCE prints a signed amplitude per
  cent.  So the target's amplitudes are squared over 100 to meet the
  source rather than the other way round, the sign question does not
  arise, and CERTAIN -- the share of one basis state that settles a
  match -- becomes 49 rather than 70, which is the same statement in the
  other space.  FLOOR still cuts the target's components in amplitude
  space, so exactly the components a deck-to-deck run would keep are the
  ones kept here.

  2. Partial vectors.  A table quotes two or three components and stops.
  explevs.distance charges every component one side has and the other
  does not, which is right when both vectors are complete and wrong
  here: it makes the nearest level whichever one is most thoroughly
  mixed, because that is the level whose components are all small.
  Several of Pr III's matches went to f4p5 levels 100 kK away for
  exactly that reason.  explevs.distance_named asks the question the
  source actually answered -- these components, and this much of
  everything else -- by comparing the named components one by one and
  lumping the rest into a single remainder on both sides.  It is used
  only for a source that is not a deck; the deck-to-deck distance is
  untouched, and so is every number the Perl comparison tests pin.

  3. Labels.  RCE writes a configuration in six columns and the source
  writes it out in full: `f5d2` against `4f.(2F*).5d2.(3F)`, `f26p`
  against `4f2.(3F<2>).6p`.  There is no rule that recovers the missing
  digits, and there does not need to be one -- the deck names its own
  configurations, the list is short, and the source's compact form ends
  with the deck's label.  The longest label that is a tail of it wins;
  a configuration no label matches is reported, not guessed at.  Terms
  are compared on the total LS term alone, RCE's `(3F) 4D3` and the
  sheet's `4D*` both reducing to `4D`, and target components that fall
  on one key are added together -- an experiment quoting `4f3 2P` is not
  distinguishing the parents RCE prints beside it.

  A label is whatever column 3 of the IN36 card called it, so it need
  not be an abbreviation at all: Pr III's `f4p5` is 4f4.5p5.  That one
  matters twice.  It is why the parity of each parity set is taken as a
  majority over the deck's labels rather than from one of them (`f4p5`
  reads even where the deck has it odd, and it is one label in
  twenty-one), and it is why an unresolvable configuration is a warning
  rather than an error.

MIXED COUPLING, WHICH IS WHAT THE REFERENCE FILE IS

  WORK/Pr3/Pr3_exp_levs.txt is the author's own sheet for Pr III and the
  first real file this part of the program has ever had.  593 levels, 14
  columns, and the designations are in whatever scheme made each level
  look pure: LS, JK (`2[1]*`), pair (`(2,3/2)*`), or nothing at all.
  The term label says which -- a bracket is JK, a parenthesis with a
  comma is pair, a multiplicity and an orbital letter is LS -- and only
  LS can be compared with an LS deck.

  Two conventions of that sheet no column list expresses, and both are
  options rather than special cases:

    term1 is the leading component's real label, used where term holds
    only the parity mark.  Hence --conf1col/--term1col: the first triad
    is (perc, conf, term) unless either of those is given, and a blank
    configuration or term in any triad means the one in the conf or term
    column.

    remark = "or" makes the perc2/conf2/term2 triad the *same* component
    written again in another scheme rather than a second component.
    Hence --altcol/--altvalue: the second triad's labels join the
    first's, its percentage is the one to believe, and the third triad
    becomes the second component.

  Of the labels offered for the leading component the first that is LS
  is taken.  If none is, the row gets no vector at all and is matched on
  (parity, J, energy), which is what the author asked for and is the
  honest answer: a JK designation says nothing an LS deck can be
  compared against, and using a second component of 25% instead would be
  guessing.  Over the whole sheet that is 343 rows matched on
  composition and 250 on energy, and the run says so:

      Pr3_exp_levs.txt: 593 rows read, 1 header line skipped; 593 levels.
         parity sets: e -> 2, o -> 1
         343 levels carry a composition the target can be matched on;
         250 are matched on (parity, J, energy) alone.
         leading components by coupling scheme: JJ 151, JK 116, LS 226,
         NONE 100
         Warning: 218 of them name one component only, which is thin
         evidence; check those matches.

MATCHING BY ENERGY ON PURPOSE

  --matchby has three settings.  `heuristic`, the default, is the above:
  composition where the source gives one, energy where it does not.
  `Ecalc` is energy alone against the deck's fitted energies.  `Eexp` is
  energy alone against the deck's *observed* energies, and only levels
  that already carry one are candidates -- which is how a deck fitted
  early in an analysis has its observed energies replaced by optimized
  values.  The adjustments there are small, so nearest wins; two source
  levels of one (parity, J) landing on one target level is what the
  duplicate report has always caught, and it still does.

ACCEPTANCE

  WORK/Pr3/RCEINP.fin already carries the observed energies the sheet
  holds, so `--matchby Eexp` over it must put the deck back exactly as
  it was.  It does -- byte for byte, all 593 energies -- and nothing but
  a matcher that tied every source level to the level it came from would
  manage that.  It is the one oracle this branch has and it is pinned.

  The full run, with the composition columns, substitutes 593 energies,
  leaves nothing unmatched, and reports no two levels of one J given the
  same energy.

  CODE/tests/test_altlevsrc.py, 31 tests and 47 cases: the coupling
  classifier, the term and configuration keys, the parity of an
  abbreviated label, the header and column resolution, the units, the
  exclude column, each of the two sheet conventions above,
  distance_named against the case that needed it, and the two runs over
  the reference file.
  CODE/tests/test_explevs.py keeps every deck-to-deck test it had; the
  one test of the LOPT reader moved across, being the only part of it
  that was not LOPT-specific.  Full suite from CODE/: 340 passed, 286
  skipped.


LEVEL IDENTIFIERS IN THE conv_out TABLES
========================================

A level of a published analysis is known by an identifier -- the key the
database holds it under, `059003.000065` in the author's sheets -- and
that name outlives every calculation the level appears in.  conv_out
cannot know it.  The only name a level has in its output is `Lev_#`, and
that is a position in this run's ordering: add a configuration, refit,
run again, and every number has moved.  Anything downstream that has to
say "this level, the one the last paper called ..." has had to be done
by hand.

CODE/insert_level_ids.py carries the identifiers in from the file that
has them.  The source is the same tab-delimited experimental level list
subst_exp_levs takes -- WORK/Pr3/Pr3_exp_levs.txt is the reference one,
and its `level_id` column is where the name of the option's default came
from -- and the targets are the two files of a conv_out run.

  insert_level_ids <source> <levels file> [<lines file>]
      --idcol <col>  --ecol <col>  --jcol <col>  --pcol <col>
      --excludecol <col>  --units cm-1|kK|eV  --skip <n>  --tol <value>

WHAT IT IS MATCHED ON

  (parity, J, energy), which is what subst_exp_levs falls back to for a
  source level whose composition it cannot use.  It is weaker evidence
  there than an eigenvector and it is the right evidence here, for a
  reason that is worth writing down: the `Ee` column of the level list
  holds the observed energies of the RCE deck, and those energies came
  *from* a list like this one, through subst_exp_levs.  The two sides
  are not being compared.  They are being recognised, and the tolerance
  exists only for the case where the two files are different vintages of
  the analysis and a level has moved by a hundredth of a wavenumber
  since.  It is one wavenumber by default, in the source's own unit, and
  --tol changes it.

  The match is made over the offers rather than level by level.  Every
  (source level, list level) pair within the tolerance is collected, they
  are taken in order of how far apart they are, and each level of either
  side is used once.  Level by level, two source levels of one J that
  straddle one list level would be settled by which of them came first in
  the file; this way the nearer one takes it and the other is reported
  unmatched.  An `Ee` that is the calculated energy -- the star in `e_c`
  -- loses a tie to an observed one, and a match made to one is counted
  and reported, because it is the one kind of match that is more likely
  to be an accident than a recognition.

WHERE THE COLUMNS GO

  `level_id` goes into the level list beside `Ee`.  A level is printed
  once per coupling block, the same `Lev_#` under a different `Cpl`, so
  the identifier goes on every row of that level and the blocks stay in
  step.

  `id1` and `id2` -- the lower level's and the upper level's -- go into
  the transition list beside `lid1` and `lid2`.  Nothing is matched a
  second time to put them there: `lid1` and `lid2` *are* the `Lev_#` of
  the level list, so the two levels of a line are already known and the
  identifiers are looked up.  A line naming a level the list has not got
  is reported rather than guessed at, which is what a truncated level
  list would produce.

  Columns of those names already in a file are updated, not duplicated.
  They keep their place, are emptied, and are filled with what this run
  matched -- so a level that has lost its identifier since the last run
  loses it here too, instead of keeping a stale one that nothing would
  ever contradict.  A second run over unchanged files reproduces them
  byte for byte, and that is a test.

  Both files are rewritten through toolbox/files.write_deck, so the line
  endings survive and what was there is kept as `<file>.bak`.  Rows are
  written the way conv_out writes them, trailing delimiters stripped, and
  a row is padded only as far as an insertion needs: a table this has not
  changed a field of comes back out as the bytes it went in as.

  The header is found by the columns it names rather than by where it is,
  because the line list carries a preamble line above its header and the
  level list does not.  A file that is neither -- the cascade matrix of
  model CM, the AUGER spectrum, which has an `id_up` and no `lid1` --
  says so and stops.

TESTS

  CODE/tests/test_levids.py, 48 cases.  The J reader in every way a J is
  written; the units and the column resolution; the identifier column by
  name and by number; a repeated identifier, which is reported because
  these are meant to be unique.  Then the matching: parity telling two
  levels of one J and one energy apart, J telling two levels of one
  energy apart, the nearer of two source levels taking the target and the
  further one being reported, the tolerance in the source's unit, and an
  observed energy winning a tie with a calculated one.  Then the tables:
  the identifier on every coupling block, the two columns landing after
  `lid2`, the blank header field they must not be inserted in front of,
  the preamble kept, an orphaned `lid`, a second run reproducing the
  files byte for byte, a stale identifier cleared, and CRLF surviving.
  Full suite from CODE/: 388 passed, 286 skipped.


RELATIVE-ENERGY GROUPS IN RCE
=============================

An experimental level list is not one measurement.  A compilation carries
systems of levels that are known to a hundredth of a wavenumber inside
themselves and tied to the ground level only through one resonance
transition measured to a few hundred -- the `+x` systems.  Until now RCE
had one way of treating a level: either it was fitted at unit weight
against its calculated energy, or it was not fitted at all (`*` in column
14).  Both answers are wrong for a `+x` system.  Fitted at unit weight,
its unknown offset is charged to the parameters, and since some of those
parameters are shared -- ZETA of an open d shell is one linked group
across every odd configuration of a Ni-like ion -- the error does not stay
where it was made.  Left out, the tens of accurately spaced levels inside
the system are thrown away with it.

The feature fits the system as what it is: precise internally, loose
absolutely, with the two accuracies stated separately.

THE MODEL

  For a group g with base b and members i,

      E_exp,k = delta_g + eps_k        for every k in g

  where delta_g is the common offset -- the error of the tie to ground --
  with standard deviation sigma_base, and eps_k is the internal error,
  standard deviation sigma_int.  Whitening that covariance gives n
  residuals for a group of n levels, not n-1.  The base stays in the fit;
  it is simply the weakest constraint present:

      member i      (E_exp,i - E_exp,b) - (E_calc,i - E_calc,b)   w_int
      base b        E_exp,b - E_calc,b                            w_base
      base of an untied group   (no residual)                     0
      blank level k E_exp,k - E_calc,k                            1

  with

      w = sigma_0^2 / (sigma_0^2 + sigma^2)

  for each of the two sigmas.  The form saturates at 1 as sigma -> 0 and
  vanishes as sigma grows, so a perfectly known quantity enters the fit as
  it always did and an unknown one does not enter at all.  Both limits are
  special cases of one number, which is the sign that the parametrization
  is the right one.

  Two approximations are made rather than hidden.  The exact whitening of
  sigma_base^2 * 11' + sigma_int^2 * I is not diagonal -- the base and each
  internal difference are correlated -- and the off-diagonal terms are
  neglected; they are O(sigma_int^2 / sigma_base^2), about 0.2% at 20
  against 500 cm-1.  And a blank level is implicitly sigma = 0, its whole
  error taken as model error.  This package has never carried per-level
  experimental uncertainties and this does not introduce them.

sigma_0 IS TAKEN FROM THE FIT, NOT FROM A CARD

  sigma_0 is the deviation the model itself achieves, and it is read off
  the previous cycle's STDEV.  That is an ordinary variance-component
  argument -- the weights are fixed self-consistently by the scatter of
  the fit -- where a constant punched on a card would be a magic number.
  It also anneals in the safe direction: poor starting parameters give a
  large STDEV, hence weights near 1 and the bases fully weighted, and a
  base is released only as the fit earns the right to release it.  Before
  the first cycle STDEV does not exist, sigma_0 is zero, and every weight
  is 1, which is the sigma_0 -> infinity limit of the same formula and
  needs no starting decision.

  The feedback is positive in both directions, so it is guarded rather
  than trusted.  sigma_0 never grows (the ratchet), which removes the
  upward branch outright and makes the sequence monotone, bounded below,
  hence convergent; and once it stops moving by more than one per cent the
  weights are frozen for the rest of the run, which removes the slight
  optimism that weights estimated from the same data put into SDX.  The
  ratchet is irreversible, so an anomalously low STDEV in an early cycle
  would lock in weights that are too small.  Nothing detects that
  automatically, which is why both the raw STDEV and the adopted sigma_0
  are printed every cycle -- a clamp is then on the page, and the remedy is
  to restart from better parameters.

THE INPUT

  The letter goes in column 14 of a level record of `RCEINP`, the column
  that already held `*`:

      blank   an ordinary experimental level, fitted at unit weight
      *       not fitted
      A-Z     the base of that group
      a-z     another member of that group

  A level cannot carry both a letter and `*` -- the column is one
  character -- so every lettered level is a fitted one and a base can never
  also be a member, without anything having to check it.

  The two sigmas have nowhere to live on a level record, so they come on a
  card block at the end of the deck, after the parameter cards, one block
  per configuration set:

       GROUP SIGMA_BASE  SIGMA_INT
      X         500.000     20.000
      Z                     15.000

  The header is any line containing `GROUP`; one card per group follows,
  the letter in column 1 and the two sigmas in columns 2-13 and 14-25; a
  blank line, a `-1` line, end of file, or anything that is not a card ends
  the block.  The sigmas are in the deck's own energy unit, the UENRGY of
  `ING11` -- kK in most decks, so 500 cm-1 is punched as 0.500.

  They are read list-directed from those fixed columns rather than under
  `F12.3`, so that a card punched `500` means 500 and not 0.5; a value with
  an explicit decimal point reads the same either way.

  A blank SIGMA_BASE means the group has no measured tie to the ground
  level at all -- its offset is unknown, not merely poorly known -- and its
  base then contributes no residual.  That is why a SIGMA_BASE of zero is
  rejected instead of being taken for the same thing: zero reads as "known
  exactly", which in the weight formula means precisely that.  A blank
  SIGMA_INT is weight 1.

  The block is at the end because the parameter cards were the last thing
  ever read from unit 20, so no old deck can contain anything there.  An
  old deck hits end of file and yields no groups, which is why the read is
  on `IOSTAT`.  A deck with more than one configuration set has the next
  set's control card in that position; that is not a header, so it is
  backspaced over and left for the read that expects it.

WHAT IS REJECTED, AND WHAT IS ONLY WARNED ABOUT

  RELRD stops on a duplicate group card, on a sigma that is not a number,
  on a SIGMA_BASE that is not positive, and on a negative SIGMA_INT.
  RELCHK, which runs once per configuration set after the levels are in,
  stops on a member letter with no base anywhere, on two bases for one
  letter, on a letter with no group card, on a group card for a letter that
  is on no level record, and on an untied group that leaves a parameter
  free.  Every one of them names the group and the card or the level.

  Two things are legal and are said out loud instead.  A group that is its
  base alone downweights one absolute level and gains no internal
  difference -- legal, rarely what was meant.  And if every fitted level
  carries a letter, no residual enters at unit weight and nothing holds the
  scale of the weights up from below; the fit still runs, but sigma_0 is
  then estimated from data that is itself all downweighted.

  THE GAUGE TEST is the one that needs explaining.  A group with a
  SIGMA_BASE imposes an absolute constraint, weak but real, and needs no
  analysis.  A group without one has a free offset, and that becomes a free
  direction in parameter space if some parameter can shift that group alone
  -- in practice the EAV of a configuration whose fitted levels all lie in
  untied groups.  The fit would then be singular in that direction and
  nothing downstream would say so: LSS catches only an exact zero pivot,
  and the B(L,L)=0 -> 1 guard in CALCV repairs a null column silently.  So
  it is decided by counting, before the first cycle, and RCE stops and
  names the configuration rather than choosing a gauge on the user's
  behalf.  A level counts as touching a configuration if any eigenvector
  component listed for it in `RCEINP` belongs to that configuration; only
  the largest MAXCR components are listed, so a very small admixture can be
  missed, which errs towards not stopping -- the right direction for a test
  that aborts a run.  An EAV that is fixed now but carries a MAX.VALUE, and
  so may be released later in the run, is a warning rather than a stop.

  The test sees an exactly free direction and not a nearly free one.  A
  deck can therefore pass it and still oscillate, which is what the untied
  `WORK/U90` trial did.

WHAT IT PRINTS

  `OUTE` gains a group summary per configuration set --

       RELATIVE-ENERGY GROUPS - HETEROSCEDASTIC WEIGHTING OF THE ...
          2 GROUP(S):
       GROUP  LEVELS    BASE  J BLOCK    SIGMA_BASE     SIGMA_INT

  -- and, every cycle, the STDEV it just achieved, the sigma_0 it adopted,
  whether that is frozen, and the two weights of every group.  Those five
  numbers are the whole state of the mechanism, and they are on the page so
  that a fit that went somewhere unexpected can be read rather than guessed
  at.

  Nothing else moves.  `RCEOUT` is what it was plus the letters in column
  14 and the group card block; copied onto `RCEINP` it reports the same
  groups and converges to the same AVDEV.  A deck with no letter anywhere
  gives byte-identical output, which is the acceptance standard of
  CLAUDE.md section 1 and was verified on all six committed decks at every
  stage of the work.

TWO THINGS THAT ARE NOT TRUE

  w = 1 does NOT reproduce the unweighted fit.  It reads as though it must
  -- the weights are all 1 and the residuals are the same numbers
  rearranged -- and it was measured not to.  Differencing is a linear map
  on the residual vector whose A'A is not the identity: it reparametrizes
  the residual space but not the least-squares problem.  Two internal
  differences of one group are correlated by one half, at every weight, and
  that correlation is the first of the two approximations above, which does
  not vanish as sigma -> 0.  On `WORK/U90`, a group at exactly unit weight
  converges to AVDEV = 210.160 against the letterless deck's 191.931, and
  each solution is the minimum of its own objective and worse under the
  other.

  AVDEV and STDEV of a lettered deck are NOT comparable with those of the
  same deck without letters.  They are weighted quantities over differenced
  rows.  Compare two fits by the O-C of the levels you care about --
  `WORK/Mo15/README.txt` does exactly that -- and not by the number at the
  foot of the page.

THE PYTHON SIDE

  `CODE/toolbox/rce.py` keeps two questions apart that column 14 used to
  answer at once.  `observed` means "has an experimental energy", and blank
  and every group letter answer yes while `*` answers no.  `absolute`
  answers "on what scale", and only a blank level is absolute.
  `comparable(a, b)` is the one primitive the seven difference sites call:
  two levels of one group may be differenced, because that difference is
  what sigma_int makes well determined, while across a group boundary the
  unknown offset makes it meaningless.

  `subst_exp_levs` reads a `+x` suffix on a source energy -- `12345.678+x`
  -- and writes the letters into column 14 of the target deck.  Which level
  is the base is a decision about the measurements, so it is asked for with
  `--basecol` and never guessed; a system of one level is its own base.
  Before this, such an energy simply failed to parse and the row was
  dropped without a word.

  In the other direction, `prepub` writes the letter of every level of a
  group into a `suffix` column after `Ee`, as ` +x` with its leading space,
  and `pub` carries that column to the published list.  `conv_out`'s
  `levels.txt` gains a `grp` column when a letter exists.

  No Python writes the group card block.  A deck that acquires its first
  letters from a source therefore needs one added by hand before RCE will
  run it -- RCE stops and says which letter has no card.

TWO NOTES FOR ANYONE SCRIPTING RCE

  gfortran's `STOP 'string'` exits 0.  Every one of the validation stops
  above therefore leaves the exit status at success, and a script cannot
  tell a rejected deck from a finished one by it; it has to read what was
  said.  This is not new with this feature -- `STOP 'NORMAL EXIT'` and
  every other STOP in the package have always behaved this way -- but it is
  worth writing down.

  `RCEOUT` exists even after a stop.  RCE opens unit 21 with the other
  units before it has read a card, so a stop in RELCHK still leaves an
  empty `RCEOUT` behind.  "Did the run finish?" is therefore asked as "does
  `RCEOUT` contain a parameter block?".

THE WORKED EXAMPLE

  `WORK/Mo15/` is Ni-like Mo XV: 3d10 plus 3d9.nl, with 22 known levels, of
  which the five 3d9.4f and 3d9.5f levels are known to only 400-600 cm-1
  while the 4s and 4p levels are known to under 0.1 cm-1.  ZETA(3d) is one
  linked group across all the odd configurations, so fitting the f levels
  at unit weight pulls a parameter the 4p levels depend on.  Its README.txt
  gives the run step by step and the results; the deck itself,
  `RCEINP.example`, is the format above in use.

  One piece of guidance out of it, because it is the thing a first user
  gets wrong.  SIGMA_INT is not the experimental uncertainty.  It is
  weighed against sigma_0, and sigma_0 is estimated from the same fit those
  levels are in, so a group bad enough to need downweighting also inflates
  the sigma_0 it is being weighed against.  In Mo XV the f levels are
  quoted at 500 cm-1; punched as 0.500 the anneal stalls, sigma_0 freezing
  at 0.666 kK and w_int at 0.640, barely a downweight at all, and the 12
  fitted 4p levels come out at a mean |O-C| of 148 cm-1 against 87 for the
  unweighted fit -- worse than doing nothing.  Punched as 5.000 the anneal
  breaks out: sigma_0 falls to 0.093 kK, w_int to 3.5e-4, and the 4p levels
  come out at 58 cm-1, better than either.  The rule to take away is that
  SIGMA_INT is a knob calibrated against the residual scale of the levels
  you trust, and belongs at ten times or more the AVDEV you expect from the
  rest of the fit.  It is not a licence to inflate every group: the
  well-measured 3d9.4s group in the same deck is quoted honestly at
  0.080/0.015 and behaves.

TESTS

  `CODE/tests/test_relgroups.py`, 25 cases, is the executable form of the
  verification matrix and the place to look for what was checked, rather
  than a restatement of it here.  It compiles `RCE20KD.F` itself with
  `build.sh`'s flags and builds `TAPE2E` for `WORK/U90` from one run of
  RCN, RCN2 and RCG, so nothing in it depends on a generated file being
  present.  The committed pre-feature `CODE/RCE20KD.exe` serves as an
  independent oracle for the legacy bit-exactness case, and that test skips
  with a reason once a release replaces it.  The residual-sum case
  recomputes R, AVDEV and STDEV from RCE's own `T-W` column and the weights
  it reported, without being told the answer.  Full suite from CODE/: 440
  passed, 286 skipped.


ONE calc_TP FOR E1, M1 AND E2
=============================

`calc_TP_M1E2` is retired.  What replaces it is `calc_TP`, the same
program covering the allowed lines as well as the forbidden ones, and no
alias is kept: no `.par` file for the old name exists anywhere, so there
is nothing to be compatible with.  `CODE/calc_TP.py`, `calc_TP.bat`,
`CODE/unix/calc_TP`, `CODE/calc_TP.par`, and the explicit lists in
`CODE/unix/install.sh` all carry the new name.  `CODE/perl_legacy/` is
reference and was left alone, so the Perl there still reads
`calc_TP_M1E2.pl` and the Python still cites it by that name at the
points where the two differ.

The program answers the question Cowan's codes do not: how well is this
transition probability known?  There is no analytic answer, so it is
obtained by resampling -- displace the fitted Slater parameters within
the standard deviations RCE printed, displace the radial transition
integrals by a fractional guess, re-run RCG, and collect the spread over
many trials.  The Perl did that for M1 and E2.  Doing it for E1 is not a
matter of adding a third section name, and the six things below are why.

WHAT THE MULTIPOLE COSTS

  E1 crosses the two RCE parity sets, so both have to be read and the two
  ends of a line looked up in different ones; M1 and E2 never leave one.
  `SPECTRUM_OF` names RCG's three sections, `SAME_PARITY` says which stay
  inside a set, and the `multipoles` parameter -- TEXT, default
  `E1,M1,E2` -- chooses.  A one-parity deck asked for all three says so
  and drops E1 rather than failing, which is what lets `WORK/Forbid` run
  under the default.

  `merged_list` and the `E2_frac` column are gated off whenever E1 is in
  the run.  The merge adds the M1 and the E2 line of one level pair into
  one row, and there is no third partner for it to take; generalising it
  was possible and gating it was less risk.

THE EIGENVECTORS OF A TRIAL ARE PAIRED WITH THE UNVARIED RUN'S

  The Perl re-identified each trial's levels against RCE from scratch.
  `levels.identify` matches positionally above the highest observed
  fitted energy, and a varied trial moves the eigenvalues there, so the
  positions are not the same list twice.  Measured: Pr III trial 2 of 5
  died with "RCE level 62 ... falls on RCG eigenvector 62, which is
  already taken".

  `CODE/toolbox/overlap.py` pairs a trial's eigenvectors with the
  unvaried run's instead, on |<base|trial>| within each (parity set, J)
  block, largest overlaps assigned first, and then composes with the base
  run's identification.  The LS basis is the deck's and is identical in
  both runs, so a block is one matrix product.  Per trial the run reports
  how many were paired, the median overlap, and the counts below 0.7 and
  0.3; on Pr III, 5626 eigenvectors, median 0.995-1.000, 0-33 below 0.7,
  none below 0.3.

  numpy is a hard requirement for `calc_TP` because of this, and a
  missing one stops the program with the pip line rather than degrading.
  It earns it: `overlap.remap` over all 5626 eigenvectors is 0.15 s
  against 15.7 s for the trial around it.  `make_cfp`, where numpy is
  optional and only true jj coupling is lost without it, is unchanged.

THE TRACKED QUANTITY IS S, NOT A

  gA carries the trial's own eigenvalues inside it -- sigma**3 for a
  dipole, sigma**5 for a quadrupole -- so a level the fit did not observe
  drags its lines with the cube or the fifth power of wherever it lands
  in that trial.  The spread of A over the trials would then be the
  spread of the matrix element and the spread of the energy together, and
  the second is a different question.  Each trial's printed gA is
  therefore turned back into a line strength at THAT TRIAL'S OWN
  calculated energies, and it is S that is accumulated; A is quoted from
  the mean S at one fixed wavelength.

  Where both levels are observed the wavelength does not move between
  trials at all and the relative spread of S is exactly that of A.  The
  two differ only where a level is unobserved, which is the case S was
  chosen for.

  The two level-energy columns were `E1` and `E2`; with E1 a multipole
  name as well they are now `E_low` and `E_up`.

  A line the trials printed once, or not at all, carries a BLANK standard
  uncertainty in every column that has one.  One sample has no spread and
  a zero there would read as a line known exactly.  `n_trials` is printed
  beside them, and is load-bearing: a line seen in fewer trials than were
  run was censored out of the rest.

CENSORING, WHICH A FORBIDDEN RUN NEVER NOTICED

  DMIN and the FLBDN/FLBDX window censor RCG's printed list, and a line
  just above either is printed in some trials and not others.  Its mean
  is then a mean over a censored sample -- the geometric mean of a
  censored sample is biased high -- and 15 per cent of a real list is
  affected.  Measured on Pr III at DMIN 8.0, 9462 lines, 5 trials: 8017
  seen in all five, 1407 in one to four, 38 in none, and 1773 more that
  appeared only in trials.

  Three answers, and they are separate.

  `trial_dmin_offset` (1.0 decade) runs the trials below the reporting
  cut, and the reporting cut is applied at write time, so a line that
  dips below DMIN in one trial still has that trial's value.  Whose DMIN
  it is depends on the multipoles: with E1 in the run the forbidden lines
  are part of the same spectrum and share it; with the forbidden lines
  alone there is no reason for a threshold and the trial decks are
  punched at 0.0.

  The window is widened, and NOT SYMMETRICALLY.  A level energy moves by
  a per cent or two and `trial_lambda_margin` (5 per cent) is right at
  the short end, but a transition energy is a DIFFERENCE, and between two
  highly excited configurations it is small, so the same absolute
  movement is a far larger fraction of it -- and FLBDX cuts on exactly
  those smallest transition energies.  The long end is widened by

      f_long = trial_lambda_margin * (Eav_up + Eav_low) / dEav_min

  over the configuration pairs the unvaried spectrum actually connects,
  bounded by `trial_lambda_margin_max` (10.0).  The energies in it are
  EXCITATION energies, measured from the lowest configuration average
  energy of the run, and that is not a detail: `outg11`'s Eav is Cowan's
  own, the total electronic energy of the configuration, about -2.026e6
  kK on Pr III, of which 350 kK is the entire spread of the 39
  configurations and the rest is the binding of 57 electrons no fitted
  parameter moves.  Taken literally the formula gave 65000 and every deck
  sat on the ceiling.  Measured from the ground configuration the same
  Pr III pair gives 2.985 -- FLBDX 16000 -> 63760 A.  `FLBDX_FIELD_MAX`
  (99999 A) is the other guard: the E5.0 field cannot hold more, and
  `write_f_wide` would have punched five asterisks.

  E0MAX is left alone.  It is a scope choice, not a threshold noise
  pushes a line across.  At DMIN 3.0 that is worth qualifying: the list
  then reaches the top of a 150 kK window, and 403 of the 1838 short
  lines had E_up above 149 kK against 4615 of all 70428.  Whether the
  trial decks should carry a raised E0MAX is the author's decision and
  is not taken here.

WHAT THE LIST HOLDS

  The unvaried run's transitions, all of them and only those.  The wider
  trial window brings in lines the deck's own cuts keep out, and those
  are not results -- the reporting cut IS the unvaried run's deck, and
  the widened window exists to complete the statistics of the lines that
  are in it, not to add lines to it.  Membership is therefore decided by
  the base run alone, and the gA threshold is subsumed by it.  `min_bf`
  still drops a line: that is a physics choice the user made, not
  censoring the trials introduced.

  A trial may print a line with its two levels the OTHER WAY ROUND --
  which end is upper is decided by energy in that run, and an unobserved
  level moves.  The reader tries the reversed key before it calls a line
  new, so the base transition keeps that trial in its statistics instead
  of losing it to the membership rule.  `E_low`/`E_up` stay the unvaried
  run's, which is the calculation being published.  Pr III shows 0
  reversals at either DMIN, so the path is exercised by a synthetic
  spectrum in the tests.

THE LEVEL LIST, AND WHY THE JOIN IS NOT ON Lev_#

  The line list names its two levels by their leading LS designation and
  their J, which is what a reader wants and not enough to join on.  What
  the author wanted was the uncertainties of a calc_TP run carried into
  a conv_out line list of the same ion -- S, the energy-adjusted gA,
  std_A% -- and the two tables had no column in common.

  calc_TP now writes a level list beside its line list (`out_levels`,
  `tp_levels.txt` by default) and puts `lid1` and `lid2` in the line
  list, between the calculated energies and the wavelength, where
  conv_out puts them.  The numbering is `levels.numbered`, which is the
  rule conv_out numbers by and is now the one place it is written down:
  conv_out's `_collect` takes its levels from it in that order and
  `_number_levels` does nothing but drop what is above `upper_bound` and
  count off the rest.  The sort key and the level set are unchanged, so
  conv_out's tables are what they were.

  The two files then have the shape `insert_level_ids` reads -- `Lev_#`,
  `Ee`, `e_c`, `J`, `par` in one and `lid1`, `lid2` in the other -- so
  the identifiers of a published analysis go onto calc_TP's lines as
  they go onto conv_out's, and THAT is the column the two tables are
  joined on.

  Not `Lev_#`.  A level number is a position in one run's energy
  ordering, and the measurement is on the author's Pr III pair -- a
  calc_TP run made in a scratch folder of its own, which was not kept.
  Applying conv_out's own numbering rule to the calc_TP run gives a
  different number for 121 of the 5626 levels.
  Every one is an adjacent transposition of two levels of DIFFERENT J
  about a thousandth of a kK apart --

      position 395 -> conv_out Lev_# 396   o  J=7.5  Ec=133.0495
      position 396 -> conv_out Lev_# 395   o  J=3.5  Ec=133.0506

  -- so a join on the number would silently mispair exactly the levels
  most easily mispaired.  The identifiers do not move, being a property
  of the analysis rather than of the run.

  What moved the eigenvalues in that pair is worth recording, because it
  looked at first like a bug in one of the two programs.  The observed
  energies are bit-identical in both files, all 593 of them, and both
  files agree on which levels are observed.  What differs is 3752 of the
  5033 UNFITTED levels, whose `Ee` column holds the calculated energy
  and not an observed one -- their `Ee` differences are the same numbers
  as their `Ec` differences, to the digit -- and 376 of the fitted
  levels' `Ec`, by up to 0.0044 kK.  The two runs are different
  calculations: the conv_out tables predate the change of the
  fundamental constants (section 4 of CLAUDE.md), and the constants move
  every calculated energy a little.  Two runs of one ion in one folder
  are not the same run, and the join has to survive that.

THE TOTAL DECAY RATE COMES FROM RCG

  The branching fraction wants the total radiative decay rate of the
  upper level, and summing the printed lines gets a DMIN-censored total
  exactly where DMIN is doing the most work.  RCG's own SUMGA
  accumulates at label 580, which the DMIN test jumps past the print and
  into, so it counts every line of the level.  `spectrum.totals` reads it
  from FORMAT 1059, the record that closes each `* * *` block.

  Three facts about it belong with the number.  It IS censored by the
  wavelength window, which drops a line at label 498 before it enters the
  arrays -- the run warns once when a window is in force.  It is `E10.3`,
  four figures, and where nothing is censored the two agree to that: on
  `WORK/Forbid`, DMIN 0.0 and no window, all 36 records match the sum
  over the lines they close, worst case 3.8e-4.  And a field of the
  record can overflow into asterisks: where SUMGA itself is starred the
  level falls back to the sum over its printed lines and the log names
  it; where only a neighbouring SUMFI field is, SUMGA is still good and
  is taken.  Nothing is substituted silently.

  A level has a SUMGA only where it HEADS a printed block, and ISPECC
  (column 72) chooses which sorts print.  A deck punching 2 leaves half
  its levels without one -- 550 of Pr III's 1046.  Punching 3 gives all
  of them, at the cost of printing every transition twice, so
  `read_transitions` no longer refuses a repeat: the first print of a
  transition in one OUTG11 contributes, every later print of the same
  unordered pair and multipole is checked against it and dropped, and a
  disagreement beyond `SAME_LINE` is a warning naming the transition.
  One OUTG11 still contributes one value per transition.  Measured on
  Pr III at ISPECC=3: 9462 transitions from 18924 prints, every repeat
  consistent, 1213 SUMGA records, no level falling back.

  SUMGA is not, however, the quantity this program sums.  Where both
  levels are observed the printed line is rescaled to the observed
  wavelength -- which `calc_TP` has always done and `conv_out`
  deliberately does not -- and RCG summed SUMGA from its own energies.
  On Pr III the ratio ran 0.96 to 1.05.  So the total is neither:

      total = sum_rescaled * SUMGA / sum_unrescaled

  the two terms of the factor both being on RCG's own scale, which makes
  it purely what DMIN removed, applied to a sum on the reporting scale.
  `Lev` carries the two sums side by side for it.  After this, SUMGA over
  the un-rescaled sum runs 0.9995 to 1.0007 on Pr III, median 1.0000,
  none outside a per cent -- the residue being the four figures of the
  field.  It is an approximation and is documented as one: the totals
  serve the branching fractions, the branching fractions serve as a
  publication cut-off, and a per cent in a threshold of that kind changes
  nothing (author).

WARNINGS, AND --cleanup

  A few hundred trials produce warnings by the hundred and on the
  terminal they bury the progress lines.  They go to `calc_TP.log` in the
  working directory, written and flushed as they happen so that a run
  killed after six hours still says what it found, and the terminal gets
  one line at the end.  The censoring counts, the reversed-match counts,
  the low-overlap counts and the computed `f_long` all go there.

  `calc_TP --cleanup` deletes the `trial*/` subdirectories of the working
  directory and does nothing else -- no reading, no RCG, no output file.
  It is a command-line flag and not a par spec.  Trials are kept by
  default: a trial can be read again without being run again.  Because it
  is destructive it lists what it will remove and asks once; `--yes`
  skips the question, and it refuses to run where there is no `trial*` at
  all.  For a run where the disk rather than the re-reading is the
  problem there is `keep_trial_output = 0`, which deletes each trial's
  OUTG11 as soon as it has been read and leaves the deck that made it.

SCALE, MEASURED

  Pr III, the author's working DMIN of 3.0, window 800-16000 A,
  ISPECC=3, 25 trials, `keep_trial_output = 0`:

    6 min 32 s total, 15.7 s a trial, of which RCG is ~11 s and Python
    ~4 s.  Peak resident 1784 MB.

    `overlap.remap` over 5626 eigenvectors is 0.15 s.  Reading the
    eigenvectors it works on is 1.4 s and `Case.read` of the unvaried run
    3.1 s, so the matrix products are one per cent of a trial.

    2.8 GB of trial OUTG11 written and deleted as it went; the directory
    ended at 596 MB, of which 270 MB is the unvaried run.  Peak occupancy
    is one trial's OUTG11 above that, so the footprint is bounded by the
    base run and does not grow with the trial count.  At 200 trials the
    same arrangement writes and deletes ~22 GB and ends at ~2.9 GB.

    120274 base transitions read from 240548 prints, 70428 written after
    `min_bf` = 1e-4, 25172 trial-only rows dropped, 1254 SUMGA records
    with no level falling back, 0 reversals.

  1838 of the 70428 written lines (2.6 per cent) were seen in fewer than
  25 trials, against 0.9 per cent at DMIN 8.0.  Only 64 are above 8000 A,
  so the long-wavelength window is not the cause; 524 are within one
  decade of the reporting cut, which `trial_dmin_offset` = 1 is exactly
  the answer to and cannot improve on, the trials having already run at
  `DMIN_FLOOR`.  The remaining 1275 are neither, and their median std_A%
  is 324 per cent: these are cancellation-dominated lines whose gA swings
  by more than a decade between trials, and no window setting reaches
  them.  The wide spread the file quotes IS the answer for such a line.

  The run itself was made in a scratch folder of its own, not in
  `WORK/Pr3/`, and was deleted afterwards rather than committed.  That
  is the rule and not an accident of this run: the reference cases ship
  as inputs only and the suite looks them up by name, so a run left in
  `WORK/Pr3` makes its tests execute against a deck that is no longer
  the documented one -- measured, 4 failed and 32 errors.

RCG'S OWN LINE STRENGTH, AS A CHECK ON radiate

  Nothing in the flow above compares `toolbox/radiate` with RCG, because
  the Lund layout `calc_TP` reads prints no line strength at all, for any
  multipole.  The wide layout that the wavelength sort always prints does
  -- `S/PMAX**2` is `S2(I)`, and the line strength is `S2*PMAX**2` from
  `PSQ = GFCF*(PMAX**2)` and `GF = PSQ*S2(I)*ERAS1` -- and the record
  carries LAMBDA(A) and GA(SEC-1) beside it, so the check is
  self-contained in one OUTG11:

      radiate.s_from_ga(GA, LAMBDA, multipole) == S2 * PMAX**2

  It is a test, not part of the production path.  Precision: S2 is F8.3,
  three decimals ABSOLUTE, so a line is worth 5e-4/S2 -- 5e-4 at S2 = 1,
  half a per cent at 0.1, nothing below 1e-3, which prints as 0.000.
  PMAX at F8.5 adds ~1e-5, doubled by the square.  So it is a part in
  2000 on the strong lines, asserted at one per cent on S2 >= 0.1, and it
  is not a bit-exactness check and is not to be sold as one.  Every error
  it exists to catch -- a wrong multipole coefficient, a wrong power of
  the wavelength, a missing 2J+1, an angstrom/cm-1 slip -- is a factor of
  a few at least.

  ONE LIMITATION, FOUND IN THE WRITING.  PMAX is printed by FORMAT 96
  under `IPTG`, and `rcg11kd.f:13341` is

        IF (II.NE.1) IPTG=0

  with II the multipole flag of SPECTR -- 0 M1, 1 E1, 2 E2 -- so an M1 or
  an E2 section prints neither PMAX nor the PMUP table it is the maximum
  of.  The absolute check is therefore made on the E1 spectrum of
  `WORK/`, and the two forbidden multipoles get the strongest form that
  survives without PMAX: within one section PMAX is a single number, so
  `S(gA)/S2` must be CONSTANT over the section, and that constant is
  PMAX**2.  That is not weak.  The section spans a wide range of
  wavelength and of 2J+1, and the wavelength power is exactly what
  differs between the multipoles: reading `WORK/Forbid`'s E2 section with
  the E1 coefficient spreads the ratio by a factor of 1860, and its M1
  section with the E2 one by 475, against a tolerance of 1 per cent.
  What it cannot see is a wrong constant prefactor.  Printing PMAX for
  II /= 1 would close that, and is an additive output change nobody has
  asked for; it is written down here rather than made.

TESTS

  `CODE/tests/test_calc_tp.py` and `CODE/tests/test_overlap.py`.  Between
  them: the multipole table and the parsing of the `multipoles`
  parameter, the parity-crossing map, the overlap matcher on a permuted
  synthetic block and on a sign flip, the DMIN policy switch and the
  trial window, `long_margin` including its ceiling and its measurement
  from the ground configuration, the SUMGA parser and its starred
  fall-back on a deliberately corrupted OUTG11, `total()` as the rescaled
  sum times the SUMGA ratio and each of its fall-backs, the blank
  uncertainty at n_trials <= 1, a doubly-printed transition folded and a
  mismatched one warned about, a reversed-order line matching its base,
  the membership rule, the log, `--cleanup`, and the S2 cross-check
  above.  The end-to-end cases build `WORK/Forbid` from its two input
  cards by running RCN, RCN2, RCG and RCE into a temporary directory.
  The level list, the lid columns and the insert_level_ids round trip
  over calc_TP's own two files are covered there too.
  Full suite from CODE/: 503 passed, 291 skipped.


THE ROUND TRIP THROUGH IDEN2
============================

A questionable level is looked at in IDEN2.  Getting it there is the
problem this pair of tools solves.

The analysis runs in a circle.  conv_out writes ENLEV.DAT and TRANS.DAT
from a finished RCG/RCE run; IDEN2 is where the observed lines are
looked at against the calculated ones; what comes out of that goes
through the author's line-classification suite and then through LOPT,
which optimises the levels by least squares and says which lines it
used.  Then the circle closes -- and it closed by hand, because the two
ends do not speak the same language about a level.

  * LOPT knows a level by its identifier, `059003.000065`, the key the
    database holds it under.  It is the same name in every paper and
    every run.
  * IDEN2 knows a level by a number, and it numbers its own levels, from
    the highest energy down, over whatever conv_out gave it.  The number
    is a position in one pair of files.  Rebuild them -- add a
    configuration, refit, run conv_out again -- and every number moves.

So the first tool writes down the correspondence and the second uses it.

  write_IDEN_level_ids <source>
      --idcol <col>  --ecol <col>  --jcol <col>  --pcol <col>
      --excludecol <col>  --units cm-1|kK|eV  --skip <n>  --tol <value>

  update_IDEN <LOPT levels> [<LOPT lines>]
      --ucalc <value>  --ids <file>  --wntol <value>

WRITING THE TABLE

  The source is the same tab-delimited experimental level list
  subst_exp_levs and insert_level_ids take, and the columns are named
  the same way -- WORK/Pr3/Pr3_exp_levs.txt is the reference one.  It is
  the only argument.  ENLEV.DAT is read from the working directory under
  that name, because that is the name IDEN2 gives it, and the table goes
  beside it as IDEN_level_ids.txt: two columns, `level_id` and
  `IDEN_id`, tab separated, in level order.  Neither input is touched.

  Only the FIXED levels of ENLEV.DAT are matched -- the ones carrying
  the star in column 40, which in IDEN2's terms are the experimental
  ones -- and the match is `insert_level_ids`' own, on (parity, J,
  energy) with every pair inside --tol offered and the closest taken
  first.  The energy compared is the FIXED energy, not the calculated
  one: it is the observed value, and an experimental level list holds
  observed values.

  PARITY.  IDEN2's label is eleven characters, `f25g _3H4I `, and its
  sixth is `_` or `~`.  conv_out puts it there in `iden._label` from
  `case.card(kset, 1).parity`, the parity of the configuration set the
  level's block belongs to: `_` for parity 0, even, and `~` for parity
  1, odd.  It is the level's parity and not a set number, whatever the
  IN36 set numbering happens to be, and reading it that way is checked
  rather than assumed -- the tool matches a second time with the marks
  reversed and says so if the wrong reading would tie more levels, which
  is what would happen if it were handed the ENLEV.DAT of a different
  calculation.

  A source level ENLEV.DAT has no counterpart for is reported and left
  out of the table.  That is the usual case rather than the exception:
  conv_out writes at most 9999 levels and truncates a large calculation
  to the ones its printed transitions touch, so a compilation reaches
  IDEN2 with part of itself missing.  Pr III is 593 experimental levels
  and 492 of them are in its ENLEV.DAT.

WHAT update_IDEN DOES WITH IT

  Three things, in one pass, with ENLEV.DAT.bak and TRANS.DAT.bak left
  behind:

  1. Every fixed level of ENLEV.DAT takes the energy LOPT gives it and,
     as its uncertainty, max(D1, D2tot) -- the larger of the purely
     statistical uncertainty and the total one relative to the ground
     level.  That is what a level is worth when it is looked at on its
     own, which is what IDEN2 does with it.  `Ef - Ec` is carried along
     although IDEN2 recomputes it and never reads it: a file that says
     one thing in two places misleads whoever opens it in an editor.

  2. With --ucalc, every level that is NOT fixed takes that uncertainty.
     The calculated levels come out of conv_out at 5000 cm-1, which is a
     placeholder; the uncertainty is the window IDEN2 searches around a
     level, and the point of the option is to narrow it once the fit is
     good enough to deserve it.

  3. TRANS.DAT takes the same energies -- the `$` record's, the `+`
     record's lower level, and the transition energy that is their
     difference -- and, for a transition between two FIXED levels, the
     assignment it should have.

  THE THIRD IS THE ONE THAT MATTERS.  IDEN2 reads exactly one field of a
  `+` record beyond the two level numbers: the assigned line's index in
  dlv.dat, in columns 67-72.  The set of transitions that carry one has
  to be the set of lines the optimisation used, or the next session in
  IDEN2 is looking at the previous round's answer.

  A line counts as used when its `Weight` in the LOPT line table is
  above zero.  That is the flag to rely on and not the `F` column: `P`,
  `M` and `virt` all mean a line was excluded, but the weight says so
  for an excluded line that carries no letter as well, and a weight that
  rounds to zero was worth nothing to the fit whatever it was before the
  rounding.  Where one transition appears on several rows -- the same
  pair measured twice, a blend entered twice -- the row with the largest
  weight speaks for it.

  LOPT names the line by its observed wavenumber and IDEN2 by its row of
  dlv.dat, so the two are tied by bisecting the observed list, which
  runs downward in wavenumber.  --wntol is how far apart they may be:
  0.05 cm-1 by default, because LOPT prints a wavenumber rounded to the
  decimals its uncertainty deserves and the two then disagree in the
  last place and no further.  On Pr III the worst of 7612 rows is
  0.03 cm-1 out, and none is unresolved.

  That is checked over the whole line table before anything is
  written, and the count and the worst offset are reported.  The
  check the assignment pass makes on its own would not do: it sees
  only the used lines whose transition TRANS.DAT lists between two
  fixed levels, which on Pr III is 2823 of 7612 rows, so a dlv.dat
  that is simply the wrong file -- a later, renumbered list -- would
  show up as a handful of unassigned lines rather than as what it is.
  dlv.dat has to be beside TRANS.DAT in any case, the assignment
  being a row number of it, and nothing is refused on the strength
  of the check: a tolerance set too tight is a few near misses and a
  mismatched file is thousands, and the two read differently.

  A transition between two fixed levels that the optimisation did not
  use loses its assignment, and the observed intensity, wavenumber and
  dwn beside it are zeroed with it, so that the record does not go on
  quoting an observation it no longer has.  A transition that reaches
  only one fixed level, or none, can carry no Ritz wavenumber and was
  never in the fit; it is none of LOPT's business and keeps whatever
  IDEN2's own session gave it.

  WHAT IT CANNOT DO.  A line the optimisation used whose transition is
  not in TRANS.DAT cannot be assigned there, and the count is reported.
  conv_out writes only the transitions above its intensity threshold and
  only between levels that survived the truncation, so on the Pr III
  files of August 2026 that is 2823 of the 2954 used lines -- the tools
  put 84 new assignments in and cleared 1, on top of the 47 conv_out's
  own automatic matching had made.  The remedy is a fuller TRANS.DAT,
  not a fabricated record: TRANS.DAT holds calculated transitions, and a
  row invented for an observation has no calculated intensity to put in
  it.

THE FILES, AND WHY THE READERS ARE WHERE THEY ARE

  toolbox/iden.py wrote ENLEV.DAT and TRANS.DAT and now reads them back.
  The column spans are stated once, as ENLEV_COLUMNS, TRANS_LEVEL_COLUMNS
  and TRANS_LINE_COLUMNS, and the `%` formats of `_write_enlev` and
  `_write_trans` a few hundred lines above them are what they have to
  agree with; WORK/Pr3/IDEN2_file_formats.md is the third statement of
  the same thing and the one a user reads.  Splitting the readers into a
  module of their own would have made three places to change a column
  instead of two.

  ENLEV.DAT has two layouts, differing by one column throughout.  Only
  version 2 -- the four-column level number -- is read; version 1 is the
  legacy one, IDEN2 keeps whichever it was given, and nothing in this
  package writes it.  They are told apart by the slash that opens the
  label, column 58 in version 1 and 59 in version 2, because no other
  field can hold one.

  A record is rewritten field by field into the text it was read as,
  never rebuilt from a format.  Two reasons.  The format document marks
  several fields "not used by the program", and IDEN2 writes into them;
  and a record nothing was changed in has to come back out as the bytes
  it went in as, trailing blanks included, or the author's diff of a run
  is unreadable.  ENLEV.DAT records are 71 columns and TRANS.DAT's are
  72 and 44, and each keeps the width it had.

  IDEN_level_ids.txt is required to be a bijection: a repeated
  identifier or a repeated level number is an error and not a
  last-one-wins, because translating through a broken table would put
  two levels' energies quietly into one record.  A table whose numbers
  ENLEV.DAT has not got is refused by name -- it was made from a
  different ENLEV.DAT, which is the one mistake this pair of tools makes
  easy.

TESTS

  `CODE/tests/test_iden2.py`.  Neither ENLEV.DAT nor TRANS.DAT is
  committed -- they are generated, and the Pr III pair is 1.6 MB -- so
  every fixture is built in the documented columns.  Between them: the
  column reader against a hand-built record, the version 1 refusal, a
  file nothing was changed in coming back byte for byte under both line
  endings, an update touching one field and no other, the table over a
  spectrum whose two levels differ only in parity, a source level that
  is not in ENLEV.DAT, a calculated level that must never be given an
  identifier, max(D1, D2tot) resolving to each of the two in turn,
  --ucalc reaching only the calculated levels, an assignment gained and
  an assignment cleared in one run with its observation cleared with it,
  the transition energies and dwn following the new levels, a second run
  changing nothing, the backups, and each of the three ways a table or a
  LOPT file can be wrong.


THE CALCULATED GROUND LEVEL AT ZERO FOR RCE
===========================================

A fit with RCE is easier to start with the calculated ground level at
zero, where the observed one is.  The first step used to be: run RCG,
read the lowest level, `add_eav -E0` on ING11, run RCG again, then RCE.
The second RCG run had one purpose, to write OUTGINE again with the
shifted values, and it is no longer made:

  add_eav --ground [folder]
  rce --ground

WHY OUTGINE IS ENOUGH

  The first RCE run of a fit, INVECT=0, takes its starting parameters
  and its calculated levels from OUTGINE (RCE20KD.F:390-470).  TAPE2E
  holds the matrices, the transformation matrices and the labels, but no
  parameter value and no eigenvalue (rcg11kd.f:15390, :15633), and a
  shift common to every Eav moves every level by that amount and no
  eigenvector at all.  So OUTGINE is rewritten -- every Eav and every
  listed level moved by minus the lowest level -- and RCG's file is kept
  as OUTGINE.RCG.  (Not OUTGINE.BAK: update11 uses that name for its own
  save and restore.)  After the first RCE run RCE reads RCEINP
  (INVECT=1), and update11 restores the OUTGINE it saved, so the shift
  matters on the first run only, which is where it is wanted.

  Doing it inside RCG instead was considered and rejected.  RCG does not
  know the ground energy until both parities are diagonalized -- ENERGY
  prints parity 1 before it reads parity 2 -- so moving OUTG11's levels
  costs a second ENERGY pass; and RCG cannot tell a first run from a
  later one, since the files a folder holds prove nothing about how its
  ING11 came to be.

  ING11 and OUTG11 keep RCG's values.  printout run straight after the
  first RCE run therefore shows RCG's Eav differing from RCEOUT's, and
  RCEINP.HF holds the shifted Eav; the first update11 carries the fitted
  values into ING11 and from then on they agree.  An LSF/HF ratio means
  nothing for an Eav in any case.

THE RULES

  * The ground is the lowest level listed, over all configuration sets;
    of equal lowest energies, the first.  It is reported on the screen,
    e.g. "Ground level adopted: level 1 with J = 0 of configuration set
    1, E_RCG = -14.4602".  OUTGINE lists the levels of a J block by
    position and never names the J.  genealogy.txt names every block, in
    the same order, and is believed when its blocks of the set agree
    with OUTGINE's in number and, for the block concerned, in size;
    otherwise the block is named by its position.  Naming the level by
    its configuration and term would mean reading OUTG11, which is not
    worth it here.
  * INVECT not 0 in OUTGINE's first record -- the only one RCE reads it
    from -- means RCE reads RCEINP and not OUTGINE: a warning, nothing
    shifted, exit 0, and rce --ground goes on to run RCE.
  * Levels all zero is what RCG writes with no diagonalization (IEWR=0):
    an error, exit 1, and rce --ground does not run RCE.
  * A ground already at zero is left alone, so a second --ground writes
    nothing: the shift is idempotent by construction.
  * The rewritten file is one RCG could have written.  An Eav is
    re-encoded by the rule RCG writes CVPAR with (rcg11kd.f:15836-15866):
    1PE10.2 near zero except for the first parameter of a set, F10.4,
    and past 99999.99 F10.1 for that first one and F10.2 for the rest.
    The level lines of a set take the FORMAT RCG picks from the largest
    |parameter|, recomputed from the shifted values (:15798-15802).
    Every other field keeps its characters, and a shift of zero gives
    the file back unchanged.

THE EAV <= -4000 EXCEPTION, AND W70

  RCG's own Eav shift, DELEAV of the rescale card, leaves an Eav at or
  below -4000 alone (rcg11kd.f:11279) and does not rescale it.  Cowan's
  documentation says what it is (DOCS/Cowan/RCG_DOC.txt, the autoionization
  section): not an energy but the mark of a continuum configuration.
  With 9 in column 75 of G5INP, RCN2 punches -9500 or -8500 as the Eav
  of the continuum configurations; RCG splits the basis into discrete
  and continuum states on Eav > -4000 (:11481), zeroes the CI between
  them, computes autoionization rates from the saved CI elements, and
  sorts the rates into the < -8000 and < -4000 sums used for dielectronic
  recombination.  A continuum Eav enters the free-electron kinetic
  energy (:11302-11306) by its value, so a shift that moved it would move
  those energies.

  add_eav --ground follows RCG: an Eav at or below -4000 is not moved.
  A level is a continuum level only in a configuration set that has
  such an Eav, and then only at or below -4000 -- because a real level
  can lie there.  W70 (W69+) has every Eav >= 0 and its ground level,
  2s2 2p 2P1/2, at -8391.10: the 2p spin-orbit splitting puts it 8391 kK
  below the average.  A test on the level energy alone took it for a
  continuum level; the first version of this tool did exactly that, and
  W70 is why the rule is on the Eav.

  Two things found on the way, neither changed in the code:

  * RCG's SPECTR skips every level at or below -4000 (rcg11kd.f:13593,
    :13600, :13613) whether or not the deck has a continuum
    configuration.  In W70 the ground level gets no transitions at all:
    -8391.10 appears in OUTG11's eigenvalue print and in no line list.
    Any level at or below -4000 on the deck's energy scale loses its
    lines this way -- in practice a highly charged ion, whose first Eav
    is 0 and whose ground level lies far below it.  Shifting ING11 so
    the ground is near zero and running RCG again brings the lines
    back: in W70, add_eav 8391.1026 then rcg gives the ground level 36
    line entries where it had none (2131 -> 2205 rows in all), which is
    also what the first update11 of a --ground fit does.  The author's
    decision: RCG stays as it is, since its output needs that shift to
    be sensible anyway, and the behaviour is documented instead -- in
    DOCS/MANUAL_ADDITIONS.md, "Levels at or below -4000: no
    transitions", and in WORK/W70/README.txt.
  * add_eav <shift> on ING11 moves every Eav, a continuum one included,
    where DELEAV does not.  And ING11's Eav field holds a value past
    99999.99 only to 0.1 (F9.1 and a group digit, RCN2's own rule), so
    the old route lost up to 0.05 there; OUTGINE carries the digits.

ACCEPTANCE

  RCN, RCN2 and RCG were run in a scratch directory on U90, Sn VIII
  (WORK/), W70 and Mo15, and add_eav --ground on each OUTGINE:

  * Outside the Eav fields and the level lines no line changed, and
    every moved value moved by the shift exactly -- in W70 to 0.0044,
    where an Eav pushed past 99999.99 goes out as F10.2, RCG's own
    precision there.  A second --ground changed nothing.
  * RCE was run on the shifted and on the unshifted OUTGINE of each,
    with the same TAPE2E.  The RCEINPs it wrote agree in every J, flag,
    composition and non-Eav parameter; the calculated levels and the Eav
    differ by the shift to 0.0009 (0.005 in W70, from those F10.2
    values).
  * Against the old route -- add_eav <shift> on ING11, RCG again -- the
    OUTGINE is identical in Sn VIII; Mo15 differs only in the sign of a
    zero ground (RCG's rediagonalized one is -0.0000); U90 and W70
    differ by 0.03-0.05 in levels whose Eav the old route rounded to
    0.1 in ING11.
  * rce.bat --ground was run on U90: the shift, then RCE, whose own
    message confirms the first block is J = 0.5.  The error path (exit
    1, RCE not run) and the warning path (RCE run) were checked through
    rce.bat and, under Git Bash, through unix/rce.

TESTS

  `CODE/tests/test_outgine_ground.py`.  The decks are built record by
  record in RCG's FORMATs: the shift and the characters it leaves alone,
  the second shift that is zero, a shift of zero giving the file back,
  a ground in the second set moving the first, the INVECT warning, the
  all-zero error, a continuum Eav and its levels left where they are, a
  W70-like level below -4000 with no continuum Eav, the level FORMAT
  following the largest parameter after the shift, 1PE10.2 for a small
  Eav but not the first, and the J named from genealogy.txt, refused
  from one that does not match, and replaced by the block's position
  without one.  The reference runs' OUTGINEs are taken where they exist
  and skipped where they do not.

THE ING11 READER, FIXED THE SAME DAY

  toolbox/ing11.py took every line blank in columns 71-72 after a
  parameter card for a continuation.  The -55555555. card is blank there
  too, and follows the last parameter card directly in a deck with no
  radial cards -- one configuration, no transitions -- so add_eav failed
  on it.  Continuation cards are now read by count, as RCG reads them
  (rcg11kd.f:11267, :11344): none unless the card has more than five
  parameters, then seven to a line.  The old and the new reader agree
  card for card on W70, U90, Pr3, Pr3/JK, Mo15 and Sn VIII, and each
  still writes back byte for byte.


VERSION 26.2.2: THE WINDOWS EXECUTABLES REBUILT
===============================================

26.2.1 is commit 1e46e71 (tagged v26.2.1), the last before 2026-11-10.
26.2.2 carries the --ground work and the ING11 reader fix above, the
documentation moved out of README.md into DOCS/ (Cowan's manuals into
DOCS/Cowan/), and the four Windows executables rebuilt so that their
version resource says 26.2.2.  The repository is public now and users
run released versions, so CODE/version.txt no longer says the number
need not be kept current; each release bumps it, gets an entry in
DOCS/CHANGES.md and a tag, and rebuilds the executables.

WHAT THE REBUILD FOUND

  The committed RCE20KD.exe was older than its source.  It was last
  committed in e63a6c1 (2026-08-21), and the relative-energy groups went
  into RCE20KD.F in c4e5104 (2026-08-26) without a new executable.  The
  binary holds a message string of the pre-c4e5104 source and none of
  the strings c4e5104 added.  Run on WORK/Mo15's RCEINP.example it stops
  at RCE20KD.F line 361, "Bad value during integer read" -- the group
  letter in column 14.  Windows users of 26.2.1 could not use groups at
  all; Linux and macOS build from the sources and were not affected.
  Keeping the pre-feature binary was deliberate while the feature was
  developed: CODE/tests/test_relgroups.py compares a deck with no
  letters against it, byte for byte, as the oracle for "nothing
  changed", until a release replaces it.  26.2.1 went out without the
  replacement.  With the rebuilt binary that test now skips, as it was
  written to (one more skip in the suite); the Pr III comparison below
  is the same check, made against the 26.2.1 binary itself.
  rcg11kd.exe was committed with its last source change (a8183bc), and
  RCN36K.exe and RCN2K.exe after theirs.

  The rebuilt RCN36K.exe and RCN2K.exe are 17 KB and 5 KB smaller than
  the ones they replace, from sources unchanged since those were built;
  so the old two were not built exactly as gfortran_O3.bat builds now.
  Their output is identical (below), so this was not pursued.

THE BUILD

  From FOR/, gfortran_O3.bat RCN36K and RCN2K, gfortran_O3_OMP.bat
  rcg11kd and RCE20KD -- serial and OpenMP as before, checked against
  the DLL imports of the old binaries -- with GNU Fortran 13.2.0 (WinLibs,
  Brecht Sanders r8), the compiler DOCS/BUILDING.md pins.  The one
  warning is the expected -fno-automatic/-frecursive one.  All four carry
  FileVersion and ProductVersion 26.2.2 and their program's description.

ACCEPTANCE

  Pr III (WORK/Pr3, the full sequence of its README, both RCE runs) with
  the 26.2.1 executables and with the rebuilt ones, same machine:
  TAPE2N, ING11, TAPE2E, OUTGINE, genealogy.txt, RCEINP.HF, RCEINP,
  RCEOUT, LEVELS1-3 and PARVALS byte-identical.  OUT36, OUT2, OUTE and
  OUTG11 differ only in timing lines: every differing line, its digits
  masked, is one of the time formats, and the AVEEIG and LEVELS FIT
  lines that carry a time are identical with the time removed.  So the
  new RCE source, run without groups, is bit-exact with the old binary,
  as a deck with no letter in column 14 must be.  WORK/Mo15, fitting
  pass: the old RCE stops as above, the new one fits, final odd-set
  AVDEV 0.068 kK -- the 68 cm-1 Mo15's README gives for that case.
