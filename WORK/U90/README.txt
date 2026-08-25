WORK/U90 -- U XC, the reordered-subshell example
======================================================================

Five configurations of lithium-like uranium, U89+, both parities: the
ground 1s2.2s and 1s2.2p, and the K-shell-excited 1s.2s2, 1s.2p2 and
1s.2s.2p that sit about 96 keV above them.  It is a small, fast, real
least-squares fit, and it is the worked example of **reorder_ing11** --
reading a finished calculation in a different intermediate coupling by
permuting the subshell slots of its deck.

INPUTS ONLY

    IN2, IN36     the deck, as for any other case in WORK/
    OUTGINE.fin   RCE's control file, converged
    RCEINP.fin    RCE's vector-input file, converged

As in WORK/Pr3, the two .fin files are the state the fit converged to,
and the sequence below reproduces its last pass rather than the dozens
of hand-adjusted iterations that got there.  OUTGINE.fin asks for five
RCE cycles and sets INVECT=1 in columns 21-25 of its first line, which
is what tells RCE to take the levels and the parameters from RCEINP.

WHAT REORDERING IS

RCG builds the LS basis by coupling the subshells of a configuration in
the order the ING11 configuration card lists them, and the intermediate
terms printed inside a label -- the (3P) in 1s.2p2.(3P) 4P -- are the
running terms of that build.  Permute the slots and the basis changes.
The eigenstates do not, and neither do the energies: what changes is
the set of states they are expanded in, and therefore the percentages
and the names.

reorder_ing11 does that permutation on a finished, fitted ING11.  It is
useful because the deck's order is rarely a physical choice.  Here the
deck runs 1s, 2s, 2p, so 1s.2s.2p is built by coupling 1s and 2s first
and its label names the (1s2s) parent term.  Slot order 231 puts 2s and
2p first and names the (2s2p) term instead.  Both are exact; they are
different questions.

THE RUN

Work in a copy of this directory, not in the distribution.

The commands below are typed the same way on every platform.  Only
`copy` is Windows' own: read `cp` wherever this file says it.  Nothing
else differs.

    rcn
    rcn2
    rcg
    rce                             <- do not skip this one

The first RCE run fits nothing -- there are no observed energies in the
OUTGINE that RCG wrote -- but it writes RCEINP.HF, the Hartree-Fock
parameter file that printout later divides the fitted values by to get
the LSF/HF ratios.  Nothing else recreates it.

    copy OUTGINE.fin OUTGINE
    copy RCEINP.fin  RCEINP
    rce

Now it fits, five cycles per parity set, and converges at the third:

    Iteration No 0 finished with AVDEV=6791.912
    Iteration No 1 finished with AVDEV= 221.783
    Iteration No 2 finished with AVDEV= 192.032
    Iteration No 3 finished with AVDEV= 191.931

and 190.618 for the second set.  Those numbers are in the unit of the
run, 1000 cm-1, so the fit is good to about 1.9E5 cm-1 -- 24 eV out of
96 keV, 0.025%.  Read them as cm-1 and they look catastrophic; they are
not.  There is one unit in this package and it is 1000 cm-1 everywhere
above levels_pub.txt.

    update11

Nothing has to be typed after it: update11 reads IN36 and substitutes
both parity sets because both parities are there (RCN2 reports the
split: configurations 1-3 are the first, 4-5 the second).  The trailing
2 of the Perl days is still accepted and still means both.  update11
saves OUTGINE, runs s11 to substitute the fitted Slater parameters from
RCEOUT back into ING11, re-runs RCG so that the eigenvectors belong to
the fitted parameters, restores OUTGINE, and runs printout.

THE REORDERED RUN

Build the permuted deck in the parent directory, where the fitted ING11
and its OUTG11 are.  Make the sub-folder first -- reorder_ing11 writes
the file it is given and does not create the directory above it, and an
empty folder does not survive a git checkout or a zip:

    mkdir reordered
    reorder_ing11 reordered/ING11 231

    reordered/ING11 written: 3 of the 8 slots moved.

231 is the old slot numbers written in their new order: the new first
slot is the old second (2s), the new second is the old third (2p), and
the new third is the old first (1s).  Slots the string does not name
keep their places, which is why only 3 of the 8 moved.  OUTG11 has to
be there because a parameter card carries no names -- its values are
identified by position alone -- and RCG prints the names nowhere else.

Then run RCG and RCE in the sub-folder:

    cd reordered
    rcg
    rce
    cd ..

RCG reads its own ING11 and needs nothing from the parent.  RCE runs
another cycle 0: the parameters are already the fitted ones, so its
levels are the fitted levels and its job here is only to identify them.

THE LEVEL LIST

    printout levels.txt --dominant-config
             --couplings LS,JJ,LSr --from LSr=LS@reordered
             --components 10 --min-percent 0.5

(on one line).  LSr is a name of our own: --from LSr=LS@reordered says
"take coupling LS from folder reordered and call the block LSr".  The
short form --from LSJK=JK, used in WORK/Pr3, names the coupling and the
folder with one word and works only when they are the same.  A
reordered run is the case where they are not: both folders hold LS.

    LS, JJ: 6 J blocks, 19 levels
    LSr (LS) from reordered, matched by energy-order
       set 2 J=0.5  worst dE = -0.0252 cm-1
       set 2 J=1.5  worst dE = -0.0251 cm-1
       set 2 J=2.5  worst dE = -0.0250 cm-1
       set 1 J=0.5  worst dE = -0.0138 cm-1
       set 1 J=1.5  worst dE = +0.0096 cm-1
    wrote levels.txt (57 levels)
       15 basis states are printed in the deck's summation order, not
       in n,l order: an intermediate term makes that order visible, so
       reordering the subshells would change what the label says

Two of those messages carry the whole method.

**matched by energy-order.**  The Pr3/JK merge pairs levels by their LS
eigenvectors, because both folders expand in the same basis.  A
reordered folder does not, so that is impossible, and printout falls
back to the n-th level of a given (parity, J) being the n-th level of
that (parity, J) in the parent.  It chooses the mode by comparing the
two decks' configuration cards and reports which it used; --match
forces it.

**worst dE.**  The check that the fallback is safe.  If the reordering
were exact the two runs would give identical energies; the largest
disagreement anywhere here is 0.03 cm-1 out of 8E8, which is roundoff.
A reordering that had gone wrong would show up in this line instead of
silently mislabeling levels.

WHAT THE REORDERING SHOWS

The odd-parity 1s.2s.2p levels, in the deck's order and in 231:

    Ec           J    LS                        LSr
    773687.92  0.5    90 1s.2s.(3S).2p 4P*      90 2s.2p.(3P*).1s 4P*
                      10 1s.2s.(3S).2p 2P*       7 2s.2p.(1P*).1s 2P*
                                                 4 2s.2p.(3P*).1s 2P*
    776078.40  0.5   100 1s.2s.(1S).2p 2P*      74 2s.2p.(3P*).1s 2P*
                                                25 2s.2p.(1P*).1s 2P*
    809143.29  1.5    57 1s.2s.(3S).2p 2P*      42 2s.2p.(3P*).1s 4P*
                      42 1s.2s.(3S).2p 4P*      34 2s.2p.(1P*).1s 2P*
                       2 1s.2s.(1S).2p 2P*      24 2s.2p.(3P*).1s 2P*

Same energies to the digit, same eigenstates, different expansions.  A
level that is a single pure state in one basis is a two-state mixture
in the other.

Which is right is a question with an answer, and config_purity gives
it:

    config_purity

    par configuration   Nlev  LS        JJ        LSr
    e   1s2.2s             1  *100       100       100
    e   1s.2s2             1  * 87        87        87
    e   1s.2p2             8    64      * 72        64
    o   1s2.2p             2  *100       100       100
    o   1s.2s.2p           7    85      * 99        71

For this ion the deck's own order wins over 231 -- 85% against 71% for
the only configuration where the two differ -- and jj beats both at
99%, which at Z=92 is no surprise.  That is a result, not a failure of
the tool: reordering is how the question gets asked instead of assumed,
and here it confirms that the (1s2s) parent term is the meaningful one
for a K-hole spectrum.

FROM THE FIT TO A PUBLISHED LIST

    make_levels_list_prepub
    make_levels_list_pub

reporting

    levels_prepub.txt: 19 levels; publishing 19 (e 10, o 9), to the
    highest observed level of each set
    wrote levels_pub.txt (19 levels, 3 composition groups)

No --nlev is needed, unlike WORK/W70: both sets have observed levels,
so the default cut -- the highest observed level of each set -- exists.
levels_pub.txt is in cm-1; everything above it stays in 1000 cm-1.

WHAT IT COSTS

Seconds, and under 1 MB including the reordered folder.  The whole case
is small enough to run while reading this file.

WHAT THIS DECK IS NOT

It is not a complete K- and L-shell calculation of U XC and it is not a
recommended data set.  Five configurations and nineteen levels are as
much as is needed to make a reordered run legible; a real calculation
of this spectrum would be much larger and the reordering would matter
correspondingly more.
