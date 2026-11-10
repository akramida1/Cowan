WORK/W70 -- W LXX, the true jj example
======================================================================

Twenty configurations of boron-like tungsten, W69+, both parities, in a
spectrum where LS coupling has stopped meaning anything.  At Z=74 with
69 charges removed the spin-orbit interaction is larger than the
electrostatic one, and the deck asks RCG for the coupling that says so:
**true jj**, KCPLD(8)=8.  It is also the example of a run with no
least-squares fit at all, which changes what the publication programs
can be told and what they have to be told.

INPUTS ONLY

    IN2, IN36     the deck, and nothing else

There is no OUTGINE.fin and no RCEINP.fin here, because there is nothing
observed to fit: no level of W69+ has been measured.  Everything below
is produced from these two files.

WHAT ASKS FOR TRUE jj

One field of the G5INP card of IN2, carried through to the main card of
ING11 by RCN2 at the same columns:

    columns 31-38   KCPLD   00111118

KCPLD(I)=0 asks for coupling I, so the two zeroes are LS and JJ.  The
eighth position is different: true jj is asked for by punching **8**
there, not 0, in column 38 of the ING11 main card.  Blank or 0 leaves
it off, and then nothing below that names TRUEjj will work.

JJ and true jj are not the same coupling, and this deck is the clearest
place to see why.  JJ is Cowan's J1J2: each subshell is coupled to its
own J and those are coupled together, but the subshell's own state is
still built in LS.  True jj builds it in jj -- the electrons of a
subshell are distributed over j = l-1/2 and l+1/2 first.  Where no
subshell holds more than one electron the two are the same coupling and
give the same numbers to the digit.  Where one does, they are not.

THE RUN

Work in a copy of this directory, not in the distribution.

The commands below are typed the same way on every platform.

    rcn
    rcn2
    rcg
    rce

RCN2 writes ING11 with no editing needed.  RCG reports

    TRUE jj DATA FOR   11 J BLOCKS APPENDED TO TAPE2E

and writes two files a run without true jj does not have:

    true-jj_keys.txt   the key to the 8-character basis-state labels
    genealogy.txt      one line per basis state, the coupling chain
                       the label has no room for

RCE then writes LEVELS1J, LEVELS2J and LEVELS3J beside its usual
LEVELS1/2/3 -- the true jj level tables.

The RCE run fits nothing.  It stops after cycle 0, which is all that is
wanted, and its job here is to identify the levels and to write OUTE
and RCEINP.HF.  **Ignore the AVDEV it prints.**  It reads 0.003, and
with no observed energy in the deck there is nothing for that number to
be the deviation from.

THE LEVEL LIST IN THREE SCHEMES

    printout levels.txt --dominant-config
             --couplings LS,JJ,TRUEjj --components 10 --min-percent 0.5

(on one line).  No --from and no second folder: RCG wrote all three
blocks into one OUTG11, because true jj is appended rather than
counted against the two coupling sections RCG prints.  It reports

    LS, JJ, TRUEjj: 11 J blocks, 169 levels
    wrote levels.txt (507 levels)

507 is 3 x 169, every level in all three schemes.

WHICH SCHEME EACH CONFIGURATION WANTS

    config_purity

writes coupling_choice.txt, one row per configuration and one column
per scheme, each cell the mean leading-component percentage over the
levels of that configuration.  On this deck:

    par configuration   Nlev  LS        JJ        TRUEjj
    o   2s2.2p             2  *100       100       100
    o   2p3                5    71        71      *100
    o   2s.2p.3d          23    54      * 97        97
    e   2s.2p2             8    65        73      * 97
    e   2s.2p.4f          24    53      *100       100

(five of the twenty rows).  Two things are worth reading off it.

First, LS is finished.  Wherever a level has two open subshells it is
around 50-70% pure in LS and 97-100% pure in jj.  LS keeps 100% only
where there is nothing to couple -- 2s2.2p, 2s2.3p and the rest of the
one-electron-outside-a-closed-shell configurations, where every scheme
agrees because there is only one basis state.

Second, JJ and TRUEjj agree to the printed digit in eighteen of the
twenty rows, and differ in exactly the two that hold equivalent
electrons: 2p3 (71 against 100) and 2s.2p2 (73 against 97).  That is
the whole content of the distinction.  J1J2 has to build 2p3 out of LS
terms of the 2p shell and cannot express 2p-2.2p+ as one state; true jj
can, and does, at 100%.

The stars config_purity punches are a suggestion.  Move them if the
literature on the isoelectronic sequence says otherwise, keeping
exactly one star per row and editing nothing else in the file.

FROM THE RUN TO A PUBLISHED LIST

    make_levels_list_prepub
    make_levels_list_pub --nlev 95,74

The first applies the stars and writes levels_prepub.txt (169 levels).
The second names them and writes levels_pub.txt.

**--nlev is not optional here.**  Normally each configuration set is
published down to its highest observed level, and with nothing observed
there is no such level; the program says so and stops:

    make_levels_list_pub: the o set has no observed level, so there is
    no highest observed level to cut at.  Say how many of its levels
    are to be published with --nlev (one count per set, in the order
    o/e).

95 and 74 are all of them, in IN36 set order -- configurations 1-12 are
the odd set and 13-20 the even one, which RCN2 prints when it runs.
The run then reports

    levels_prepub.txt: 169 levels; publishing 169 (o 95, e 74)
    wrote levels_pub.txt (169 levels, 3 composition groups)

and no level takes a letter prefix: at this ionization every level is
pure enough that its leading component names it uniquely.  In the
published file the true jj configurations appear split, 2p3 coming out
as 2p-2.2p+, 2p-.2p+2 and 2p+3 -- three configurations in jj where LS
sees one.

One caution about the energies.  levels_pub.txt is in cm-1, as always,
but with nothing observed RCE has no ground level to measure from and
the numbers are RCG's own, on the scale its center of gravity sets.
The lowest level of the odd set comes out at -8391102.6 cm-1 and the
highest at 198726751.8.  Subtract the ground level yourself before
publishing.  (Everything above levels_pub.txt is in RCE's 1000 cm-1.
A level at Ec=132251.7250 in levels.txt is 1.32e8 cm-1, about 16 keV,
which is the right order for an n=2 hole in W69+.)

The same scale costs the ground level its lines.  RCG takes any level
at or below -4000 in the deck's unit -- here -4 000 000 cm-1 -- for a
continuum state, Cowan's convention for autoionization runs, and
computes no transition to or from it.  So 2s2.2p 2P1/2, at -8391.10 kK,
appears in OUTG11's eigenvalue print and in none of its line lists:
every transition to the ground level is missing.  To have them, put the
ground level near zero and run RCG on that -- after rcn2,

    add_eav 8391.1026
    rcg

and the rest as above.  The ground level then has its 36 entries in the
line lists, and every energy downstream is measured from it -- to about
1 cm-1 here, since ING11 holds an Eav above 99999.99 only to 0.1 -- which
also takes care of the subtraction above.  See "Levels at or below -4000" in
DOCS/MANUAL_ADDITIONS.md.

WHAT IT COSTS

Seconds, and about 3 MB:

    RCN 0.1 s, RCN2 0.1 s, RCG 1 s, RCE 1 s, the four utilities ~5 s

TAPE2N is 1.6 MB of it and OUTG11 0.5 MB.  Neither is needed once
levels.txt exists.

WHAT THIS DECK IS NOT

It is not a fit and it is not a recommendation.  The deck carries no
observed energies, so nothing in it has been tested against a
measurement; it exists to exercise true jj and to show what a
calculation with no experimental input can and cannot be asked for.
