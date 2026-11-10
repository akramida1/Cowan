WORK/Mo15 -- Mo XV, the relative-energy-group example
======================================================================

Nickel-like molybdenum, Mo14+.  Twenty-two observed levels of 3d10,
3d9.4s, 3d9.4p, 3d9.4f and 3d9.5f are fitted against a deck that also
carries the unobserved members of the 3d9.ns, 3d9.np, 3d9.nd, 3d9.nf
Rydberg series and the 3p5.nl and 3d8.nl.n'l' configurations.  It is a
small, fast, real least-squares fit, and it is the worked example of
the **relative-energy groups** -- telling RCE that a set of levels is
known precisely with respect to one of its own members but only
loosely with respect to the ground level.

INPUTS ONLY

    IN2, IN36           the deck, as for any other case in WORK/
    Mo15_explev.txt     the observed levels, tab-delimited
    RCEINP.example      RCE's vector-input file, ready to fit

Nothing generated is kept here; the sequence below rebuilds all of it.


WHY THIS CASE NEEDS GROUPS
----------------------------------------------------------------------

Column B1 of Mo15_explev.txt is the uncertainty of each level *with
respect to the level marked B1 in that column* -- 3d9.(2D<5/2>).4s
(5/2,1/2) J=3 at 1694848 cm-1, the base of the optimization that
produced the list.  Read down that column and the structure of the
measurement is plain:

    3d10 1S0            80 cm-1     the ground level, tied only loosely
    3d9.4s          10 - 16         four levels, tight against the base
    3d9.4p           8 - 20         twelve levels, tight against the base
    3d9.4f, 3d9.5f 400 - 600        five J=1 levels, from resonance lines

The absolute uncertainty of every one of these is about 80 cm-1,
because that is what the tie to the ground level is worth.  Their
internal spacings are worth several times more than that -- except for
the five f levels, which were determined directly from observed
resonance transitions and are 400-600 cm-1 either way.

The f levels are the problem.  ZETA(nd) is linked in a single group
across every odd configuration in this deck, 3d9.4p and 3d9.4f and
3d9.5f together, so five levels known to half a thousand wavenumbers
sit in the same normal equations as twelve known to fifteen, and they
carry a long lever arm: they lie some 800 000 cm-1 above the 4p levels.
Left unweighted they pull the fitted ZETA(nd) scaling and spoil the
3d9.4p fit.  That is what the group markup is for.


THE MARKUP
----------------------------------------------------------------------

Two things go into the deck.

**Column 14 of a level record** carries the status letter: a capital
marks the level a group is measured against, the matching lowercase
marks its members, a blank means the level is fitted on its absolute
energy as always, and `*` means not fitted.  In RCEINP.example:

    X   3d9.(2D<5/2>).4s (5/2,1/2) J=3  at 1694.848 kK   the LOPT base
    x   the other three 3d9.4s levels
        the twelve 3d9.4p levels -- blank, i.e. absolute
    Y   3d9.4f 3P*1  at 2749.100 kK
    y   the other two 3d9.4f and both 3d9.5f levels

**A group-card block** at the end of each parameter set says what the
two ties are worth.  A header line, then one card per group: the letter
in column 1, sigma_base in columns 2-13, sigma_int in columns 14-25,
read list-directed, in the deck's own energy unit (1000 cm-1 here).

     GROUP SIGMA_BASE  SIGMA_INT
    X       0.080       0.015
     GROUP SIGMA_BASE  SIGMA_INT
    Y       0.500       5.000

RCE then minimizes

    R = sum over blanks of r^2
      + w_base * r_base^2
      + sum over members of w_int * (r_i - r_base)^2

with w = sigma_0^2 / (sigma_0^2 + sigma^2), sigma_0 re-estimated from
the fit's own STDEV each cycle under a ratchet that never lets it grow
and freezes it at 1 %.  The mechanism is in FOR/RCEREL.F.

Note that the AVDEV and STDEV printed by a run with groups are weighted
quantities computed over differenced rows.  They are **not** comparable
with the AVDEV of a run without groups.  Compare residuals, not AVDEV.


THE RUN
----------------------------------------------------------------------

Work in a copy of this directory, not in the distribution.  The
commands are typed the same way on every platform; only `copy` is
Windows' own -- read `cp` wherever this file says it.

    rcn
    rcn2
    rcg
    rce --ground

`--ground` puts the calculated ground level at zero before RCE starts.
It takes the lowest level RCG listed in OUTGINE -- here the calculated
3d10 1S0, at -5.8714 kK -- and shifts every Eav and every listed level
of OUTGINE by the same amount, which moves the whole calculation along
the energy scale without changing a single separation.  RCE then starts
on the observed scale.  It reports

    Ground level adopted: level 1 with J = 0 of configuration set 1, E_RCG = -5.8714

and keeps RCG's file as OUTGINE.RCG.  ING11 and OUTG11 are not touched:
they keep the ground at -5.8714 until a fitted RCEOUT is carried back
into ING11 with update11.  (Earlier editions of this file shifted ING11
with `add_eav 5.8714` and ran RCG again; the fit below writes the same
RCEOUT, PARVALS and LEVELS files either way, byte for byte.)

This first pass fits nothing -- RCG's OUTGINE has no observed energies
in it -- but it is the pass that writes the RCEINP the observed
energies were substituted into, and it leaves TAPE2E and OUTGINE in the
state the fitting pass needs.  Now install the prepared deck and turn
the fit on:

    copy RCEINP.example RCEINP

and edit the **first line of OUTGINE**: the fifth field, columns 21-25,
is INVECT.  RCG always writes 0 there; change it to 1, which is what
tells RCE to take the levels and the parameters from RCEINP instead of
starting over from OUTGINE.  Columns 11-15 of the same line hold the
number of least-squares iterations, 5, which is enough for this case.

    rce

About a minute.  The results are in OUTE, the fitted deck in RCEOUT,
and the level tables in LEVELS1 (even) and LEVELS2/LEVELS3 (odd).


HOW RCEINP.example WAS BUILT
----------------------------------------------------------------------

From the RCEINP that the first `rce` above writes, in three steps.

1.  The observed energies, with **subst_exp_levs**, matching on the
    calculated energies and on the composition.  Typed on one line:

        subst_exp_levs Mo15_explev.txt RCEINP --units cm-1
            --ecol energy --jcol j_val --pcol parity --parity-sets e,o
            --confcol conf --termcol term --perccol perc
            --conf2col conf2 --term2col term2 --perc2col perc2
            --matchby Ecalc

    All 22 levels match and none is left over.

2.  The least-squares flags, by hand: Eav of 3d10 and 3d9.4s free in
    the even set, Eav of 3d9.4p and 3d9.4f free in the odd one, the
    Eav of 3d9.5f through 3d9.9f linked as one, and the ZETA, F2, G1,
    G2 and G3 scalings linked across each Rydberg series through the
    negative link codes -51 to -56.  Everything else fixed at its
    Hartree-Fock value.

3.  The nine status letters and the two group cards, by hand.

Step 1 is reproducible; steps 2 and 3 are the analysis.  Note that
subst_exp_levs can write the letters itself when the source quotes an
energy as `12345.678+x`, which is how a compilation writes a level
system tied loosely to the ground level; this source does not, so they
were placed by hand.  The group cards it never writes.


WHAT THE WEIGHTING DOES
----------------------------------------------------------------------

Three fits, one executable, the same TAPE2E, differing only in the
deck.  "none" is the same deck with the letters and the group cards
taken out; "500/500" is the Y card set to the experimental uncertainty;
"500/5000" is RCEINP.example as shipped.

Mean |O-C| over the fitted levels, cm-1:

    configuration    n      none    500/500   500/5000
    3d9.4s           4        24         25         25
    3d9.4p          12        87        148         58
    3d9.4f           3       757       1062       1180
    3d9.5f           2       456        520        410

Within each group, rms of the absolute residual and rms of the
residual measured from the group's own base, cm-1:

                          none    500/500   500/5000
    X  3d9.4s  absolute     24         29         29
               internal     42         34         34
    Y  4f, 5f  absolute    693       1063       1155
               internal   1484       1189       1291

Iterations, and final AVDEV in cm-1:

    even set    4 / 22      3 / 25      3 / 25
    odd set     4 / 385     6 / 485     5 / 68

The linked odd-parity ZETA(3d), in kK:

    10.7807     10.8760     10.7158

Read the four tables together and the picture is this.

Group X does what it is asked to do, and does it cleanly.  The four
3d9.4s levels give up 5 cm-1 of absolute accuracy and buy 8 cm-1 of
internal accuracy, which is the right trade when the internal spacings
are known to 15 cm-1 and the tie to the ground level only to 80.

Group Y is not there to improve the description of 3d9.4f and 3d9.5f.
It is there to stop those five levels from distorting the parameters
they share with 3d9.4p.  At 500/5000 the 3d9.4p fit improves from 87
to 58 cm-1 -- better than it was with no groups at all -- while the f
residuals grow, which is exactly what should happen to levels that are
no longer being fitted on their internal structure.


WHY SIGMA_INT IS 5000 AND NOT 500
----------------------------------------------------------------------

The experimental uncertainty of the f levels is 400-600 cm-1, so the
honest-looking group card is `Y 0.500 0.500`.  It does not work, and
the reason is worth understanding before setting up a case of your own.

The weight is not sigma against the experimental error.  It is sigma
against **sigma_0, the fit's own standard deviation**, re-estimated
from STDEV every cycle.  And STDEV is inflated by the very levels the
weighting is meant to suppress.  With sigma_int = 0.500 the anneal
settles into a self-consistent stalemate:

    cycle 0   sigma_0 = 2.887 kK   w = 0.971
    cycle 1             0.803          0.720
    cycle 3             0.677          0.647
    cycle 5             0.666 FROZEN   0.640

A weight of 0.64 is barely a downweight.  The f levels still dominate
ZETA(nd), which moves the wrong way, +0.9 %, and the 3d9.4p fit comes
out worse than it was with no weighting at all: 87 -> 148 cm-1.

Raise sigma_int to 5.000 kK and the anneal breaks out of that state.
Once the f levels stop constraining the shared parameters, the odd fit
is set by the 3d9.4p levels alone, STDEV collapses, and sigma_0 with
it:

    cycle 0   sigma_0 = 2.887 kK   w_int = 0.250
    cycle 1             0.432               0.0074
    cycle 2             0.118               0.00055
    cycle 4             0.093 FROZEN        0.00035

w_int = 0.093^2 / (0.093^2 + 5^2) = 3.5e-4, which is what "let these
levels fix Eav and nothing else" actually costs.

The guidance, then:

  * sigma_base and sigma_int are **weight-setting knobs, not
    measurements**.  Quoting the experimental uncertainty is the right
    thing to do only when the group's levels are of a quality
    comparable with the rest of the fit.

  * What matters is the ratio of sigma_int to the residual scale of
    the levels you trust.  Here the odd fit stands at about 90 cm-1
    once it is healthy; sigma_int has to be well above that, not
    merely above the 15 cm-1 the 4p levels are measured to.

  * A useful rule: for a group meant to have essentially no say in the
    shared parameters, set sigma_int to ten times or more the AVDEV
    you expect from the rest of the fit.  Ten times gives w = 0.01,
    fifty times gives 4e-4.

  * Watch the SIGMA_0 column in OUTE.  If it freezes at a value close
    to your sigma_int, the weighting has stalled and the group is
    still steering the fit.  If it settles well below sigma_int, the
    group has been let go.

  * Do not read this as a licence to inflate every group.  Group X is
    quoted honestly at 0.080 / 0.015 and behaves, because its levels
    are as good as the rest of the even set.  It is the mismatched
    group that needs the larger number.


ONE MORE THING WORTH KNOWING
----------------------------------------------------------------------

The Y base's absolute residual comes out at -1.3 cm-1 with 500/500 and
+0.3 cm-1 with 500/5000, even though sigma_base = 500 cm-1 should leave
it loosely tied.  That is not the weighting.  3d9.4f has a free Eav and
only one fitted level, so the fit can put the base exactly on its
observed energy at no cost to anything else: the internal-difference
term does not depend on that Eav.  The same happens to the X base.

So wherever a group's configuration carries a free Eav and no competing
observed level, the base's absolute residual is not evidence about
sigma_base.  Look at the members' internal residuals instead.
