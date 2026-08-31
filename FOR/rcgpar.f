      IMPLICIT REAL*8 (A-H,O-Z)
      IMPLICIT INTEGER*4 (I-N)
      PARAMETER (
C     A. Kramida 2026: upper bound on KMX.  KMX**2 must fit INTEGER*4
C     (it is passed to ALLOC, whose element count NMLOCS is INTEGER*4),
C     so KMX <= 46340.  gfortran enforces this at compile time via the
C     KLC42=KMX**2 parameter below.  Static memory is NOT a constraint:
C     .bss grows by only 560 bytes per unit of KMX (all the KMX**2
C     arrays are heap-allocated Cray pointees), against a ~2 GB
C     -mcmodel=small limit.  Heap is the practical constraint: peak
C     commit is roughly 40*KMX**2 + 1.0e5*KMX + 1.6 GB bytes, i.e.
C     ~6 GB at KMX=10000 and ~19 GB at KMX=20000.
     +  KS=8,
     +  KC=200,
     +  KMX=20000,
     +  KLSI=550,
     +  KPC=4000000,
C     A. Kramida 2026: KLC1 is only the INITIAL size of NOPCCC in PRK;
C     that array is grown by RALLOC on demand, so KLC1 need not be
C     exact.  The old form KMX**2*2-(KLSI**2)*9 silently wrapped around
C     in INTEGER*4 for KMX >= 32789 (gfortran diagnoses overflow in **
C     but NOT in *), yielding a negative KLC1 and an "ALLOC: Invalid
C     memory request" abort at run time with no hint of the cause.
C     Capping KMX at 32000 inside this one expression keeps the value
C     identical for KMX <= 32000 and safely bounded above it.
     +  KLC1=MIN(KMX,32000)**2*2-(KLSI**2)*9,
     +  KLC12=KMX**2-(KLSI**2)*9,
     +  KLC2=KMX**2-(KLSI**2)*7-KPC,
     +  KLS1=KLSI+1,
     +  KLS2=KLSI+51,
     +  KLSC=50000,
     +  KJK=300,
     +  KJP=2000,
     +  KISER=4000000,
     +  KPR=10000,
     +  KJJT=50,
     +  KLSU=17,
     +  KNT=800,
     +  KLAM=500000)
      PARAMETER (
     +  KBGKS=41,
     +  KTRAN=1000,
     +  KENRGS=21,
     +  KCASE=50)
      PARAMETER (
     +  KEXC=12800,
     +  KCC=KISER/4,
     +  KLC4=2*KISER+1,
     +  KLC42=KMX**2,
     +  KMX2=KMX/2,
     +  MAXBC0=100000,                                           !800000,
     +  KW=104,
     +  KP=2*KW+2,
     +  KMULT=100,
     +  KEE11=KLAM*9/2+KC*6+KPR*2+23-KLSC*(KS+1)*5/2-KJK*KS,
     +  KEE12=KMX*5/2+KC*13/2+51+2*KPR-KLSC*(KS+1)*5/2-KJK*KS,
     +  KEE1=KEE11,                                             ! A. Kramida 7/31/2012: Set KEE1 to the greatest of KEE11 and KEE12
     +  KEE2=KLAM*9/2+KC*6+KPR*2+23-KMX*5/2-KC*13/2-51-2*KPR,
     +  KEE3=KMX*5/2+KC*13/2+51+2*KPR-KLAM*9/2-KC*6-KPR*2-23,
     +  KNSJK=1000,KCCDUM=(KPR-KMX)*KMX)
      COMMON/IDSKNO/IR,IW,ID2,ID3,ID4,ILA,ILB,IL,IC,IE
      COMMON/C1/KCPL,NOCSET,NSCONF(3,2),IABG,IV,NLSMAX,NLT11,NEVMAX,
     1  KCPLD(9),SJN,SJX,IMAG,IQUAD,UENRGY,DMIN,IELPUN,ILNCUV,IPLEV,
     2  ICPC,ICFC,IDIP,IENGYD,ISPECC,IPCT,ICTC,ICTBCD,KCPL1,COUPL(12),
     3  LSYMB(24),NOSUBC,KK,TIME0,TIME,DELT,NDIFFT,NOTOTT,
     4  SJ8MN,SJ8MX,IGDLEV,NSCRJ8,SCRJ8,SCRJ8P,NLS,NLSP,NJJ,
     5  NNN,NDIFSJ,KPAR,NPAR,IX,CAVE,SUMS20,NOPC,J1A,J1B,
     6  NIJK(KS,KC,2),NOTSJ1(KS,KNT),NOT,NLT,NORCD,A,LCI,LCIP,J1X,
     7  ICNVRT(2),NLASTT(KS),NDIFFJ(KS),NOTOTJ(KS),NTOTTJ(KNT),
     8  NTOTJJ(KS,KNT),LLIJK(KS,KC,2),IIPRTY,FLLIJK(KS,KC,2),
     9  SSOPI2,SOPI2(KC,KC),NCK(2),TEXCIT,EMEAN,ZERO,HALF,ONE,THRHAL,
     *  TWO,IPLOTC,IW6,NOICWR,NOIC(2),NOICMUP(2),IQS,INWOL,ICPTR,
     *  ISENIOR,MAXBC,E0MAX(2),
C     A. Kramida 2026: IPGEN is the genealogy-print switch, columns
C     83-84 of the main card.  See label 100 of the main program and
C     SUBROUTINE GENPRT.  Appended at the end of /C1/ so that nothing
C     already in the block moves.
     *  IPGEN,
C     A. Kramida 2026: NLTCF(J,K) is the PER-CONFIGURATION limit on the
C     number of LS terms kept for the FIRST subshell of configuration J
C     of parity set K, punched in columns 83-85 of that configuration's
C     card in ING11.  0 means the field was blank, in which case the
C     global NLT11 (columns 26-27 of the main card) applies exactly as
C     it always has.  Read and validated in SUBROUTINE NLTCHK, called
C     from label 200 of SUBROUTINE LNCUV; applied in LNCUV (the printed
C     term table) and in SUBROUTINE PLEV (the basis itself).
C     ICUTNY is set non-zero by LNCUV when any configuration actually
C     loses terms, which is what turns the true-jj rotation off.
C     Appended at the end of /C1/ so that nothing already in the block
C     moves.
     *  NLTCF(KC,2),ICUTNY,
C     A. Kramida 2026: NLSPCH is the value the user actually PUNCHED in
C     columns 23-25 of the main card, before blank was turned into 999
C     and 999 into KMX.  0 means the field was blank.  Only SUBROUTINE
C     ENERGY needs it, to keep the KCPLD(8)=9 Born second pass gated on
C     what was punched rather than on the derived NLSMAX, which is now
C     non-zero in every deck.  See label 105 of the main program.
C     Appended at the end of /C1/ so that nothing already in the block
C     moves.
     *  NLSPCH
