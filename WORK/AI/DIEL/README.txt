Working DIEL example: autoionization / dielectronic recombination, B IV
=======================================================================
Verified 2026-08-11 with RCN36K.F v36.4.1, RCN2K.F and rcg11kd.f.

What DIEL does
--------------
DIEL removes the handwork of setting up an autoionization run: it finds
the kinetic energy of the captured electron by differencing average
energies, works out which orbital angular momenta the free electron may
have, and writes the "99l" continuum configuration cards.  It lives
entirely inside RCN; RCN2, RCG and RCE see nothing unusual.

    INDIEL  --+
    IN36(-1)--+--> RCN36 --> IN36 rebuilt, with the 99l continuum cards
                        \--> TAPE2N, bound and continuum wavefunctions
                                |
                        RCN2 ---+---> ING11
                                |
      [optional KCPL=3 card] ---+---> RCG --> OUTG11

Internally it is two passes.  Pass 1 writes the control card, the ion
cards and the autoionizing cards to tape10 and calls RCN just to get the
average energies; its exit card carries iz=-999, which makes RCN return
without writing the TAPE2N end-of-file so that DIEL can read the results
back.  ("-999" in a partly written IN36 is therefore not an error.)
Pass 2 rewinds INDIEL, forms eps = Eav(atom) - Eav(ion), applies the
triangle rules, writes the decay, autoionizing and continuum cards to
tape11, copies tape11 over tape10, and lets RCN run it normally.

RCE plays no part in this: it fits parameters to observed levels, and an
autoionizing configuration has none.

The example
-----------
Recombining ion   B V  1s          (the recombined ion is B IV)
Capture           B V 1s + e-  ->  B IV 2p 9i        (autoionizing)
Radiative decay   B IV 2p 9i   ->  B IV 1s 9i        (stabilizing, E1)
Free-electron kinetic energy found by DIEL:  18.58 Ry, l = 5 and 7

Files
-----
INDIEL              input to DIEL (unit 12)
IN36                contains nothing but the "-1" exit card
IN2                 RCN2 input; its G5INP card also carries the RCG
                    control card
gcard               optional RCG "KCPL=3" card, prepended to ING11 to
                    obtain the DR rate as a function of temperature
IN36.generated      the IN36 that DIEL builds - for comparison
OUT36.reference     RCN output
ING11.reference     RCN2 output
OUT2.reference      RCN2 printed output
OUTG11.reference    RCG output

Commands - run in an empty scratch directory holding INDIEL, IN36, IN2
----------------------------------------------------------------------
    RCN                     # DIEL + preliminary + final RCN
    RCN2                    # -> ING11, OUT2
    RCG                     # -> OUTG11, TAPE2E, OUTGINE

For the optional DR rates, put the gcard on the front of ING11 after
RCN2 and before RCG.  On Windows:

    copy /b gcard+ING11 tmp && move /y tmp ING11

and on Linux or macOS:

    cat gcard ING11 > tmp && mv tmp ING11

These are the wrapper names -- rcn.bat, rcn2.bat, rcg.bat on Windows,
and the lower-case rcn, rcn2, rcg of CODE/unix on Linux and macOS.
Call the wrappers rather than the executables: RCG has to be told
where the package's data files are, and writing the COWAN.CFG that
tells it is the wrapper's job.

What to look for in OUTG11
--------------------------
    GAATOT=      total autoionization (Auger) rate of the level, s-1
    BRION=       branching ratio to autoionization
    GM*FRBAR= / GM*FABAR=      summed DR strength
    T(KEV), BRNCH= ...         DR rate vs electron temperature; printed
                               only when the gcard is present

The gcard layout is rcg11kd.f statement 7:
    FORMAT (I5,I2,I3,2E5.1,2F10.5,I5,5F5.3,F10.5)
    KCPL=3, ILNCUV, ID2T, lambda bounds, DELEKEV, EIONRY, NPTKEV,
    then up to five electron temperatures in keV.
EIONRY is the ionization energy of the recombining ion in rydbergs -
25.03 here, taken straight from the RCN output for B V 1s.  Note that
the echoed header line labels NPTKEV as "NTKEV=", so NTKEV=0 there does
not mean the temperature table is missing.

Rules DIEL enforces, and what happens if they are broken
--------------------------------------------------------
1. The autoionizing configuration(s) must have parity OPPOSITE to the
   radiative-decay configuration(s), and all decay cards must come
   first.  DIEL classifies the atom cards by parity alone: the parity of
   the first one defines the decay group, and the first card of the
   other parity begins the autoionizing group.  A set with only one
   parity used to produce no continuum cards and then abort with a bare
   end-of-file on TAPE2N; as of v36.4.1 it stops with an explanatory
   message (exit status 16).
2. NSPECT on the ion card(s) must differ from NSPECT on the atom cards;
   the ion block ends at the first card whose NSPECT changes.
3. Every orbital that appears as a non-final orbital on any card must
   appear explicitly, even with occupation 0, on every card.
4. At most 5 ion configurations and 10 autoionizing configurations.
5. EMX on the control card should bracket the kinetic energy, or the
   continuum normalization is poor.  If DIEL's eps comes out far from
   the guessed EMX, edit the regenerated IN36 and rerun from RCN2.
6. The label prefix a, b, c ... on the generated continuum cards is the
   index of the ION configuration, not of the input set.
