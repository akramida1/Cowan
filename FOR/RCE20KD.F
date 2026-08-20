C      PROGRAM RCE20K   ! Version designed to work with Gerrit van het Hof's
C                      ! TROUT program (TROUTK version).
C
C     ---------------------------------------------------------------
C     A. Kramida 2026.  WHAT IS NEW IN THIS VERSION.  EACH ITEM IS
C     EXPLAINED IN FULL AT ITS POINT OF CHANGE, AND LSJJ/README.txt IS
C     THE LONG ACCOUNT; ONLY THE ANNOUNCEMENT IS HERE.
C
C       - TRUE jj COUPLING.  WHEN RCG IS RUN WITH KCPLD(8)=8 IT APPENDS
C         A TRAILER TO TAPE2E CARRYING, PER J BLOCK, THE TRANSFORMATION
C         TO THE TRUE jj BASIS - SEQUENTIAL COUPLING OF THE
C         RELATIVISTIC SUBSHELLS - AND EIGHT-CHARACTER LABELS FOR ITS
C         STATES.  IF THAT TRAILER IS THERE, RCE PRINTS THE FITTED
C         LEVELS IN TRUE jj AS WELL, AS A THIRD COUPLING AFTER THE TWO
C         RCG WROTE.  THERE IS NOTHING TO SWITCH ON, AND NOTHING
C         PRINTED BEFORE IS CHANGED: THE TRUE jj LEVEL LISTS GO TO
C         LEVELS1J, LEVELS2J AND LEVELS3J, OUTE GAINS THE NEW
C         EIGENVECTOR SECTIONS AND LOSES NOTHING, AND RCEINP, RCEOUT
C         AND PARVALS ARE UNTOUCHED.  A TAPE2E WITHOUT THE TRAILER
C         GIVES THE RUN IT ALWAYS GAVE, TO THE BYTE.  SEE RDJJ2E,
C         MAPJJ2 AND RCEJJ2E.F.
C
C       - SPEED.  THE EIGENSOLVER (TRED2, TQL2/TQLROT) IS THREADED WITH
C         OpenMP, AND THE THREE CUBIC LOOPS OF THE OUTPUT PHASE THAT
C         THEN BECAME THE BOTTLENECK ARE BIT-EXACT REWRITES: ESORTW FOR
C         THE ENERGY REORDERING, AND IN ORDOUT2 THE COMPRESSED-COLUMN
C         g-FACTOR SWEEP AND THE BOUNDED MIN-HEAP SELTOP.  THE
C         TRANSFORMATION TO THE SECOND COUPLING SCHEME IS THREADED OVER
C         ITS K LOOP.  SUBROUTINE CLOCK NOW REPORTS WALL-CLOCK TIME
C         (SYSTEM_CLOCK) INSTEAD OF CPU TIME SUMMED OVER THREADS.
C
C       - THREE CORRECTIONS, EACH BIT-EXACT ON THE CASES IT DID NOT
C         AFFECT: RCE NO LONGER REFUSES TO FIT WHEN EVERY LEVEL IS
C         KNOWN (SEE THE NOTES AT THE DO 105 LOOP AND AT IFIRST); GAPS
C         IN THE LINKED-GROUP NUMBERING NO LONGER INFLATE THE PARAMETER
C         COUNT AND HENCE THE UNCERTAINTIES (SEE BELOW THE DO 116
C         LOOP); AND NEGATIVE EIGENVECTOR AMPLITUDES ARE ROUNDED BY
C         MAGNITUDE IN ORDER, ORDOUT2 AND OUT8 (SEE THE NOTES THERE -
C         THE PRINTED COMPOSITIONS IN RCEOUT, RCEINP, LEVELS1..3 AND
C         OUTE MOVE BY A UNIT IN THE LAST DIGIT, NOTHING ELSE DOES).
C
C     THE THREADED CODE IS BIT-IDENTICAL TO THE SERIAL CODE AT ANY
C     THREAD COUNT: IT ONLY REDISTRIBUTES WORK, NEVER REASSOCIATES A
C     SUM.
C
C     BUILD (A.Kramida 2026)
C
C       gfortran -fshort-enums -malign-double -ftracer -fno-backslash
C                -fallow-argument-mismatch -fcray-pointer -O3
C                -fopenmp -fno-automatic  -o RCE20KD.exe RCE20KD.F
C
C       -fno-automatic MUST accompany -fopenmp, which otherwise implies
C                      -frecursive and moves locals from static storage
C                      onto the stack; parts of this program rely on
C                      locals keeping their values between calls, and
C                      the results change without it.  (gfortran warns
C                      that -fno-automatic overrides -frecursive - that
C                      warning is expected.)
C       Do NOT add -ffast-math: it alters the arithmetic of the
C                      eigensolver.  -march=native is likewise not used
C                      for the distributed executable: it buys about
C                      1 pct and enables fused multiply-add, which
C                      makes the results depend on the build machine.
C       Without -fopenmp the OpenMP directives are ordinary comments
C                      and everything runs serially, giving identical
C                      results.  A -static OpenMP build is not possible
C                      with this toolchain, so an OpenMP executable
C                      must be shipped with libgomp-1.dll,
C                      libgfortran-5.dll, libquadmath-0.dll and
C                      libwinpthread-1.dll alongside it.
C
C     Set the number of threads with the OMP_NUM_THREADS environment
C     variable.  If it is not set, OpenMP uses every logical processor
C     on the machine.  Only blocks larger than 1000 levels are
C     threaded; smaller blocks run serially because the fork/join cost
C     exceeds the work.
C     ---------------------------------------------------------------
C     PROGRAM RCE20K(INPUT=OUTGINE=TAPE10=INPUT,OUTPUT=TAPE9=OUTE,
C    1  TAPE2=TAPE2E,TAPE11=PARVALS,TAPE8=LEVELS1,
C    2  TAPE7=LEVELS2,TAPE12=LEVELS3)
C
C     ---------------------------------------------------------------
C          THIS PROGRAM COPYRIGHT BY LOS ALAMOS NATIONAL LABORATORY
C               AND THE U.S. GOVERNMENT.  THE PROGRAM MAY BE FREELY
C               DISTRIBUTED GRATIS ANYWHERE IN THE WORLD, BUT
C               MAY NOT BE SOLD.
C     ---------------------------------------------------------------
C
C     RCE MOD 19  NOV, 1971  (C-SET INPUT VIA TAPE FROM RCG MODS 5 TO 9)
C               MOD 19 IS THE SAME AS MOD 18 EXCEPT WITH LSQM=+1 REMOVED
C               AND WITH RECODING TO REDUCE ARRAY SIZES
C               MOD 20 IS DATED MARCH 1990
C
C
C     PROGRAM TO VARY AS MANY AS IB PARAMETERS TO MAKE
C          LEAST-SQUARES FIT OF ATOMIC ENERGY LEVELS FOR UP TO
C          IC CONFIGURATIONS INVOLVING NO MORE THAN IX J-MATRICES,
C          EACH UP TO IA-BY-IA IN SIZE; SUM OF ALL MATRIX SIZES
C          (TOTAL NUMBER OF ENERGY LEVELS) NO GREATER THAN IV.
C
C          NOTE---REQUIRED STORAGE FOR D(IB,IB+1) SPILLS OVER
C               INTO V(IV,IB).   DO NOT REORDER D AND V IN MEMORY.
C               IV MUST BE .GE. IB.
C
C          THE VERSION FOR A VAX 8200 MADE AT LUND UNIVERSITY,
C               MARCH 1988
C
C               THE VERSION FOR IBM-PC 386 OR 486 MADE AT THE
C          INSTITUTE FOR SPECTROSCOPY RUSSIAN ACADEMY OF SCIENSES
C                     BY A.E.KRAMIDA, NOVEMBER 1993
C          MODIFICATIONS MADE:
C               - BINARY TAPE 2 TAKEN FROM RCG IS IN COMPRESSED
C                 FORM;
C               - FOR ITERATIONS, MATRIXES ARE NOT SWAPPED TO DISK,
C                 BUT TO MEMORY BUFFERS IN COMPRESSED FORM (SEE
C                 SUBROUTINES WR1 AND RD1). BUFFERS SIZES ARE MAXB1
C                 AND MAXB2.
C               - FOR EIGENVECTOR REORDERING AFTER ITERATIONS, A NEW
C                 ALGORITHM IS IMPLEMENTED (SEE SUBROUTINE ORDER)
C               - ALTERNATE INPUT FILE IN SO-CALLED VECTOR FORMAT
C                 CAN BE USED INSTEAD OF THE OLD STANDARD 'OUTGINE' FILE.
C                 NAME OF THE FILE IS 'RCEINP'; IT IS CREATED BY THIS
C                 RCE PROGRAM IN THE FIRST RUN, AND THEN MAY BE EDITED.
C
C          THIS VERSION 20K IS MADE BY A.KRAMIDA AT NIST, USA,
C              SEPTEMBER 1994 - JANUARY 1995
c          Change by A. Kramida on 31 Jan. 2011:
c              temporary fix for calculations in alternate coupling schemes.
c              When the first coupling scheme is not LS,
c              this version will work correctly only if the energy ordering
c              in RCEINP is not disturbed (the order of experimental levels
c              must be the same as the order of theoretical levels.
c          Change on 19 Aug. 2026:
c              Bounded N against NPAR's declared dimension IC in the
c              DO 3145 loop over parameters.  The unbounded index was
c              harmless under MinGW by luck and fatal under Linux; it
c              never affected a computed number.  See the comment at the
c              point of change.
c          Change by A. Kramida on 22 Jan. 2026:
c              Dynamic memory allocation routines borrowed from GRASP2K.
c              Packing/unpacking routines UNPK2, WR1, LD1 now should work
c              with any sizes of input arrays that fit in memory.
c              Recompiling may be needed only if dimensions of fixed-size 
c              arrays exceed the defined dimensions.
c              Include file renamed to RCE20KD_par.f.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG,LPG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC_PTR,BUFA_PTR,BUFY_PTR,
c     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
     1  IBC_PTR,IBA_PTR,ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
      INTEGER*8 BUFY_PTR,BUFA_PTR,BUFC_PTR,IBA_PTR,IBC_PTR
      POINTER (BUFY_PTR,BUFY(1)),(BUFA_PTR,BUFA(1)),(BUFC_PTR,BUFC(1))
      POINTER (IBA_PTR,IBA),(IBC_PTR,IBC)
      DIMENSION IBC(IX,IB),IBA(IX,1)
c
      INTEGER*2 NW,NW0,NW2
      INTEGER*2 IEVV,IEVV0,ievv2
      INTEGER*2 II1,II2,IDECCH,NMAX2
      INTEGER*1 III1,III2
      COMMON/VECT/ IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      integer*8 IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
C     A.Kramida 2026: component index moved to the FIRST position so that
C     the MAXC components of one eigenvector are contiguous in memory.
      DIMENSION IEVV(MAXC,IA,IX),IEVV0(MAXC,IA,IX),NW(MAXC,IA,IX),
     1 NW0(MAXC,IA,IX)
      POINTER(IEVV_PTR,IEVV),(IEVV0_PTR,IEVV0),(NW_PTR,NW),(NW0_PTR,NW0)
c
      CHARACTER C1*8,Y1*12
      DIMENSION C1(MAXC),Y1(MAXC)
      COMMON/TIMING/TINIT,ITINIT
C       A.Kramida 2026: KCPL(3)=8 IS THE TRUE jj COUPLING - SEE RDJJ2E.
C       COUPLINGS 1 AND 2 ARE STILL THE TWO RCG WROTE ON TAPE2E.
      INCLUDE 'RCEJJ2E.F'
      DIMENSION KCPL(3)
      dimension IEVV2(MAXC,IA,IX),NW2(MAXC,IA,IX)
C       A.Kramida 2026: NJMS(J,M) IS THE LEADING CONFIGURATION OF LEVEL
C       J OF BLOCK M, WORKED OUT ONCE FOR ALL THE OUTPUT PASSES.  SEE
C       THE LOOP THAT FILLS IT, JUST BEFORE THE PASSES.
      INTEGER*2 NJMS(IA,IX)
      SAVE NJMS
C       A.Kramida 2026: SAVE keeps these large arrays in static storage.
C       Under -fopenmp gfortran would otherwise place them on the stack,
C       which overflows at once (60 MB here).  The same applies to the
C       big local arrays of ORDER, ORDOUT1 and ORDOUT2 below.
      SAVE IEVV2,NW2
C
      ITINIT=0                      ! INITIALIZE IBM-PC TIMING
C
C     CALL ERRSET(208,256,-1,1,0,0)
C     CALL ERRSET(212,256 -1,1,0,0)
C
C
C          READ PRINT, REORDER, AND END-OF-RUN CONTROL NUMBERS
C
c      write(6,*) "IA=",IA,", IV=",IV
c      stop(0)
c
c     Allocate initial buffers BUFY,BUFA,BUFC,IBA,IBC
      IBY_SIZE=MAXBY
      CALL ALLOC(BUFY_PTR,IBY_SIZE,8)
      MAXB2U=MAXB2
      CALL ALLOC(BUFA_PTR,MAXB2U,12)
      MAXB1U=MAXB1
      CALL ALLOC(BUFC_PTR,MAXB1U,12)
      CALL ALLOC(IBC_PTR,IX*IB,4)
      CALL ALLOC(IBA_PTR,IX,4)
C       A.Kramida 2026: D_PTR IS A UNION OF DELX(IB), D(JB,1) AND
C       A1(IA,IA).  A1 IS THE ONLY LARGE MEMBER, AND ITS SECOND INDEX IS
C       A RADIAL-PARAMETER INDEX, NOT A LEVEL INDEX: IN EIGEN IT NEVER
C       EXCEEDS 3+LMAX-I1 .LE. IB+2.  SO IA*(IB+2) IS ENOUGH, AND THE
C       OLD IA**2 WAS 3.2 GB AT IA=20000 AGAINST 0.96 GB NOW.  SEE THE
C       MATCHING BOUNDS ON THE ZEROING LOOP AND ON THE I2 CLAMPS IN
C       SUBROUTINE EIGEN.
      CALL ALLOC(D_PTR,max(IA*(IB+2),JB,IB),8)
      CALL ALLOC(V_PTR,IV*JB,8)
      i_C2_dim=IA*IX
      CALL ALLOC(TW_PTR,i_C2_dim,8)
      CALL ALLOC(T_PTR,i_C2_dim,8)
      CALL ALLOC(EVAL_PTR,i_C2_dim,8)
      CALL ALLOC(W_PTR,i_C2_dim,8)
      CALL ALLOC(KW_PTR,i_C2_dim,4)
      CALL ALLOC(NF_PTR,i_C2_dim,4)
      CALL ALLOC(NC_PTR,i_C2_dim,4)
      i_CMLC_dim=IA**2
      CALL ALLOC(C_PTR,i_CMLC_dim,8)
      CALL ALLOC(A_PTR,i_CMLC_dim,8)
      i_VECT_dim=IA*MAXC*IX
      CALL ALLOC(IEVV_PTR,i_VECT_dim/2,4)    ! /2 because all VECT arrays are INTEGER*2
      CALL ALLOC(IEVV0_PTR,i_VECT_dim/2,4)   ! /2 because all VECT arrays are INTEGER*2
      CALL ALLOC(NW_PTR,i_VECT_dim/2,4)      ! /2 because all VECT arrays are INTEGER*2
      CALL ALLOC(NW0_PTR,i_VECT_dim/2,4)     ! /2 because all VECT arrays are INTEGER*2
c
      OPEN(10,FILE='OUTGINE',STATUS='OLD')
      OPEN(9,FILE='OUTE',STATUS='UNKNOWN')
      OPEN(2,FILE='TAPE2E',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(11,FILE='PARVALS',STATUS='UNKNOWN')
      OPEN(8,FILE='LEVELS1',STATUS='UNKNOWN')
      OPEN(7,FILE='LEVELS2',STATUS='UNKNOWN')
      OPEN(12,FILE='LEVELS3',STATUS='UNKNOWN')
      REWIND 2
      REWIND 8
      OPEN(21,FILE='RCEOUT',STATUS='UNKNOWN')
      IDENOM=0
      NOCSET=-999
      K3=0
      IW=9
C       A.Kramida 2026: LOOK FOR THE TRUE jj TRAILER ON TAPE2E BEFORE
C       ANYTHING ELSE IS READ FROM IT.  RDJJ2E REWINDS UNIT 2 AGAIN
C       WHATEVER IT FINDS, SO THE READS THAT FOLLOW ARE UNAFFECTED.
      KCPL(3)=8
      CALL RDJJ2E
      GO TO 94
C
   95 CALL SORT8TO7
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME
      STOP '(NORMAL EXIT)'
C
   94 READ (10,8,END=95) NOCYPR,NOCYCR,NOCYCE,IW6,INV1,IORD,ISAMPL
    8 FORMAT (8I5)
      IF (IW6.GT.0) IW=6           ! DIRECT OUTPUT TO THE SCREEN
C ******************************************************* A.KRAMIDA v
C       INV1 (FIFTH FIELD OF THE FIRST RECORD OF 'OUTGINE') SELECTS THE
C       INPUT FILE AND, SINCE A.Kramida 2026, ALSO THE LEVEL-ORDERING
C       METHOD USED AFTER EACH LEAST-SQUARES ITERATION:
C         0 - READ 'OUTGINE', BUILD 'RCEINP' FOR SUBSEQUENT FITTING;
C             FULL EIGENVECTOR RECOGNITION
C         1 - READ 'RCEINP'; FULL EIGENVECTOR RECOGNITION
C         2 - READ 'RCEINP'; RELAXED RECOGNITION - EIGENVECTORS ARE
C             MATCHED ONLY FOR LEVELS TIED TO AN EXPERIMENTAL ENERGY,
C             ALL OTHERS ARE ORDERED BY INCREASING ENERGY.  BLOCKS WITH
C             NO EXPERIMENTAL LEVEL AT ALL ARE ORDERED BY ENERGY AT
C             ONCE.  INTENDED FOR CASES TOO LARGE FOR FULL RECOGNITION.
C         3 - READ 'RCEINP'; EIGENVECTORS IGNORED, EVERY LEVEL ORDERED
C             BY INCREASING ENERGY (COWAN'S ORIGINAL SORTING), WHILE
C             STILL USING THE RICHER 'RCEINP' INPUT.
C       VALUES 2 AND 3 TRADE ACCURACY OF THE LEVEL MATCHING FOR SPEED;
C       THE FIT ITSELF USES ONLY THE LEVELS THAT CARRY AN EXPERIMENTAL
C       ENERGY, WHICH MODE 2 STILL MATCHES THROUGH THEIR EIGENVECTORS.
      IF (NOCSET.EQ.-999) THEN
        INVECT=INV1
        IF (INVECT.EQ.0) THEN
C
C         PREPARE NEW VECTOR-INPUT FILE FOR SUBSEQUENT FITTING
C
          OPEN(20,FILE='RCEINP',STATUS='UNKNOWN')
          REWIND(20)
        ELSE
C
C         OPEN THE EXISTING VECTOR-INPUT FILE FOR THE FITTING
C
          OPEN(20,FILE='RCEINP',STATUS='OLD')
        ENDIF
        IF (INVECT.EQ.2) WRITE (6,9)
    9   FORMAT (' LEVEL ORDERING: relaxed - eigenvectors matched only',
     1    ' for levels with an experimental energy,'/
     2    '                 the rest ordered by increasing energy',
     3    ' (INVECT=2)')
        IF (INVECT.GE.3) WRITE (6,11)
   11   FORMAT (' LEVEL ORDERING: by increasing energy only,',
     1    ' eigenvectors ignored (INVECT=3)')
      ENDIF
C ******************************************************* A.KRAMIDA ^
      IF (NOCYPR.EQ.-1.AND.NOCYCR.EQ.0.AND.NOCYCE.EQ.0) GO TO 94
      IF (I216.NE.2) I216=16
      K3=K3+1
C
C     LOCATE THE DESIRED CSET ON UNIT 2
C
      IDN=IDENOM
      IF (INVECT.EQ.0) THEN  ! A.KRAMIDA
        READ (10,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S ! R.COWAN
        WRITE (20,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S ! A.KRAMIDA
      ELSE                   ! A.KRAMIDA
        READ (20,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S !A.KRAMIDA
      ENDIF
      WRITE (21,10) NCDES,IDENOM,IPRNV,IPRNA,WMAX,IPRNX,IPRNSQ,S ! A.KRAMIDA
   10 FORMAT(4I5,F10.4,2I5,F10.4)
      IF (WMAX.EQ.0.D0) WMAX=1.0D10
      IF (S.EQ.0.D0) S=1.D0
      IDN1=0
      IF (NCDES.LT.0) THEN
        GO TO 95
      ELSE
        IF (NCDES.GT.0) THEN
  949     IF (NOCSET.GT.NCDES) THEN
            REWIND 2
          ELSE
            IF (NOCSET.EQ.NCDES) NCDES=0
          ENDIF
C
          IF (NCDES.NE.0) THEN
   96       READ (2) NOCSET,MMAX,LMAX,NMAX,
     *        ((NPAR(L,M),M=L,NMAX),L=1,NMAX)
            IF (NOCSET.LT.NCDES) THEN
              LMAX = LMAX-NMAX+1
              DO M=1,MMAX
                READ (2) N
                READ (2) N       ! ADDITIONAL READ BECAUSE OF USE OF PK1
                DO L=1,LMAX
                  READ (2) F
                ENDDO
              ENDDO
              GO TO 96
            ELSE
              IF (NOCSET.GT.NCDES) GO TO 949
            ENDIF
          ENDIF
        ELSE
          IF (IDENOM.LE.0) THEN
            IF (IDENOM.EQ.0) IDENOM=IDN
            IDN1=1
          ENDIF
        ENDIF
      ENDIF
      IF (IDENOM.NE.0.AND.IDN1.EQ.0.AND.INVECT.EQ.0)
     1  READ (10,14) (DENOM(L), L=1,LMAX)
      N3(K3)=NMAX
C
      DO M=1,MMAX
        DO K=1,IA
          KW(K,M)=0
        ENDDO
      ENDDO
      CALL CLOCK(TIME0)
      TIME=0.D0
      WRITE (IW,7) TIME,TIME0
    7 FORMAT (/6H0TIME=,F9.4,19H MIN,   TOTAL TIME=,F9.4,
     1  F14.0,16H LEVELS FIT WITH,F5.0,11H PARAMETERS)
      IF (IW6.LT.0) WRITE (6,7) TIME,TIME0
C
C     REWIND 4
C     REWIND 15
      DO I=1,IC
        CONFIG(I,K3)=' '
      ENDDO
      IF (INVECT.EQ.0) THEN                        ! A.KRAMIDA
        READ (10,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! R.COWAN
        WRITE (20,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! A.KRAMIDA
      ELSE
        READ (20,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! A.KRAMIDA
      ENDIF
      WRITE (21,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)   ! A.KRAMIDA
   12 FORMAT (7(A8,2X))
      WRITE (11,12) ELIDEN,(CONFIG(I,K3),I=1,NMAX)
C     I=2*LMAX
      WRITE (IW,9710) K3
 9710 FORMAT (14H Conf. set No ,I1)            ! A.KRAMIDA
      IF (INVECT.EQ.0) THEN                    ! A.KRAMIDA
        READ (10,13) (XIDEN(L), L=1,LMAX)      ! R.COWAN
      ENDIF
   13 FORMAT (7A10)
C9714 CONTINUE
  971 TNO=0.D0
      IMSUM=0
      IMMAX=0
      JLMAX=0
      KKLM=0
C ****************************************** A.KRAMIDA v
C               INITIALIZE BUFFER INDEX OFFSETS
      IAOFFS=0
      ICOFFS=0
      ICOFS1=0
      IAOFS1=0
      IA2=1
C ****************************************** A.KRAMIDA ^
C
      INOCHG = 1  ! Added by A. Kramida 8/16/04
      EMAX = 0.D0                   ! Added by A. Kramida 8/20/04
C     A. Kramida 2026: IFIXL COUNTS THE UNKNOWN (NON-EXPERIMENTAL)
C       LEVELS.  IT USED TO BE CLEARED INSIDE THE DO 105 LOOP OVER J,
C       SO WHAT SURVIVED THE LOOP WAS THE COUNT FOR THE LAST J BLOCK
C       ALONE, NOT THE TOTAL.  THE TEST "IFIXL.EQ.IMSUM" BELOW - NO
C       KNOWN LEVEL ANYWHERE, HENCE NOTHING TO FIT - COMPARES IT
C       AGAINST A SUM OVER ALL J, SO IT COULD NEVER FIRE WHEN MMAX>1.
C       CLEARED HERE INSTEAD, SO IFIXL IS THE TOTAL IT IS DOCUMENTED
C       TO BE.
      IFIXL=0
      DO 105 M=1,MMAX
        IF (NCDES.NE.0) THEN
C
C      TRANSFER DATA FROM TAPE 2 TO DISC 4.  USE A FOR TRANSFORMATION
C      MATRICES HERE AND IN 'ORDER'.
C
C         READ (2) FJ(M),IM,((A(L,LP),L=1,IM),LP=1,IM),     ! R.COWAN
C    1      (NC(L,M),L=1,IM), ((YIDEN(L,N,M),L=1,IM),N=1,2) ! R.COWAN
C         IMSQ=IM*IM                                        ! R.COWAN
C  *********************************************************  A.KRAMIDA v
          READ (2) KCPL(1),KCPL(2),FJ(M),IM,
     1      (NC(L,M),L=1,IM), ((YIDEN(L,N,M),L=1,IM),N=1,2)
          IMSQ=IM*IM
          CALL UNPK2(2,IM,A,IA,BUFA_PTR,IMSQ,IAOFS1,MAXB2U)
*********************************************************  A.KRAMIDA ^
          IMAX(M) = IM
          IMSUM=IMSUM+IM
          IMMAX=MAX0(IMMAX,IM)
          IMSQ1=(IMSQ+IM)/2
          N=1
          K=0
C         IFS=4
          LMAXP1=LMAX+1
          DO 3145 L1=1,LMAXP1
            K=K+1
c            write (*,*) 'M,L1,MAXB1U,MAXB2U=',M,L1,MAXB1U,MAXB2U
            if (M.eq.2 .AND. L1.eq.64) then
              idbg=1
            endif
            IF ((K.GT.1).OR.(N.GT.NMAX)) THEN
C             READ (2) ((C(I,J),I=J,IM), J=1,IM)           ! R.COWAN
C  **************************************************
              CALL UNPK2(2,IM,C,IA,BUFC_PTR,IMSQ,ICOFS1,MAXB1U) ! A.KRAMIDA
C  **************************************************
              IF (L1.EQ.LMAXP1) THEN
                DO L=1,IM
                  GKK(L,M)=C(L,L)
                ENDDO
C               IFS=15
              ENDIF
            ELSE
              DO I=1,IM
                DO J=1,I
                  C(I,J)=0.D0
                  IF (J.EQ.I.AND.NC(I,M).EQ.N) C(I,J)=1.D0
                ENDDO
              ENDDO
            ENDIF
C           CALL WR(IFS,IM,C,IA,Y,IMSQ1)
            CALL WR1(M,L1,IM,C,IA,BUFC_PTR,IBC,IMSQ1,
     *               ICOFFS,MAXB1U,MMAX,LMAXP1)
            ICOFS1=ICOFFS
C           A.KRAMIDA 2026: N IS A CONFIGURATION INDEX AND NPAR IS
C           DIMENSIONED NPAR(IC,IC).  NOTHING BOUNDED N HERE, SO ONCE THE
C           REAL CONFIGURATIONS WERE PAST, NPAR(N,N) READ ZERO, THE TEST
C           FIRED ON EVERY PASS AND N ADVANCED UNTIL IT INDEXED OFF THE
C           END OF COMMON/C1/.  UNDER MINGW THE WORD BEYOND HAPPENED TO
C           BE LARGE, THE RUNAWAY STOPPED AT N=201 AND THE RUN FINISHED;
C           UNDER LINUX IT IS SMALL, SO N RAN ON UNTIL THE READ LEFT THE
C           MAPPED REGION AND RCE DIED WITH SIGSEGV ON THE FIRST J BLOCK.
C           THE GUARD CANNOT CHANGE A RESULT.  FOR N.GT.NMAX THE TEST AT
C           THE HEAD OF THIS LOOP TAKES THE UNPK2 BRANCH WHATEVER N IS,
C           THE ELSE BRANCH BELOW IS REACHED ONLY WHILE N.LE.NMAX, AND N
C           IS USED NOWHERE ELSE IN THE LOOP BODY.  VERIFIED TO CLAUDE.MD
C           SS1 ON PR III: RCEOUT, TAPE2E AND OUTGINE BYTE-IDENTICAL,
C           OUTE DIFFERING IN 22 LINES, ALL OF THEM TIMING.
            IF (N.LT.IC) THEN
              IF (K.GE.NPAR(N,N)) THEN
                N=N+1
                K=0
              ENDIF
            ENDIF
 3145     CONTINUE
C         CALL WR(5,IM,A,IA,Y,IMSQ)
          CALL WR1(M,IA2,IM,A,IA,BUFA_PTR,IBA,IMSQ,
     *             IAOFFS,MAXB2U,MMAX,IA2)
          IAOFS1=IAOFFS
        ENDIF
C
 1000   IM=IMAX(M)
        KKLM=MAX0(KKLM,IM)
        IF (INVECT.EQ.0) THEN       ! A.KRAMIDA
          KX=0
 3190     KN=KX+1
          KX=KX+7
          IF (KX.GT.IM) KX=IM
          READ (10,1014) (W(K,M), K=KN,KX), TEMP
 1014     FORMAT (8F10.5)
          IF (TEMP.LE.0.D0.AND.KX.LT.IM) GO TO 3190
          READ (10,15) (NF(K,M), K=1,KX)
   15     FORMAT (7I10)
          IF (KX.LT.IM) THEN
            KN=KX+1
            DO K=KN,IM
              W(K,M)=TEMP
              NF(K,M)=-K
            ENDDO
          ENDIF
C ****************************************** A.KRAMIDA v
        ELSE
          IF (M.EQ.1)
     1      READ(20,19)F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2
C
C           READ THE EIGENVECTORS FROM THE VECTOR-INPUT FILE
C
          MAXI=MIN0(MAXCR,IM)
          DO 3300 K=1,IM
            DO I=1,MAXI
              NW(I,K,M)=0 ! -128
              IEVV(I,K,M)=0
            ENDDO
            READ (20,3193) T(K,M),ERAS(1),W(K,M),J1,J2,
     1        (IEVV(I,K,M),C1(I),Y1(I), I=1,MAXI)
 3193       FORMAT (F13.6,A1,F10.3,I4,1X,I1,10(I4,1X,A6,1X,A12))
C  16   FORMAT (F13.6,A1,F10.3,I4,A2,10(I4,1H ,A6,1X,A12))
            DO I=1,MAXI
              IEVV(I,K,M)=IEVV(I,K,M)*100
            ENDDO
            NF(K,M)=K
            IF (ERAS(1)(1:1).EQ.'*') THEN
              NF(K,M)=-K
            ELSE
              W(K,M)=T(K,M)
              INOCHG = 0  ! Added by A. Kramida 8/16/04
            ENDIF
            DO 3250 I1=1,MAXI
              IF (IEVV(I1,K,M).EQ.0) GO TO 3250
              I2=0
              DO 3195 I=1,IM
                IF (YIDEN(I,1,M).NE.Y1(I1)(1:8)) GO TO 3195
                IF (CONFIG(NC(I,M),K3)(1:6).NE.C1(I1)(1:6)) GO TO 3195
                I2=I ! -128
                GO TO 3200
 3195         CONTINUE
 3200         NW(I1,K,M)=I2
 3250       CONTINUE
 3300     CONTINUE
        ENDIF
C ****************************************** A.KRAMIDA ^
C
c        write (6,*) " read M=",M," KX=",KX," W=",W(KX,M)," NF=",NF(KX,M)
C       A. Kramida 2026: IFIXL IS NOW CLEARED ONCE, BEFORE THE DO 105
C         LOOP, SO THAT IT ACCUMULATES OVER ALL J.
        DO K=1,IM
          NFKM=NF(K,M)
          IF (NFKM.LT.0) IFIXL=IFIXL+1
          I=IABS(NFKM)
          T(K,M)=W(I,M)
          if (DABS(T(K,M)).GT.EMAX) EMAX = DABS(T(K,M))   ! Added by A. Kramida 8/20/04
          TW(K,M)=0.D0
          IF (NFKM.GE.0) THEN
            IF (NFKM.GT.0) NF(K,M)=0
            TW(K,M)=T(K,M)
            TNO=TNO+1.D0
          ENDIF
        ENDDO
  105 CONTINUE
      WRITE (IW,16) NMAX,IC,LMAX,IB,MMAX,IX,IMMAX,IA,IMSUM,IV,ICOFFS,
     1   MAXB1,IAOFFS,MAXB2
   16 FORMAT (1H ,22X,'SIZE',2X,'DIMENSION'/
     1  ' NO. OF CONFIGURATIONS=',I3,I9/
     2  ' NO. OF PARAMETERS=',I7,I9/
     3  ' NO. OF J-VALUES=',I9,I9/
     4  ' LARGEST J-MATRIX=',I8,I9/
     5  ' TOTAL NO. OF LEVELS=',I5,I9/
     7  ' C-BUFFER SIZE=',I11,I9/
     8  ' A-BUFFER SIZE=',I11,I9)
 1104 IF (NMAX.GT.IC.OR.LMAX.GT.IB.OR.MMAX.GT.IX)
     1  STOP '(DIMENSION ERROR: IC, IB OR IX TOO SMALL)'
      IF (IMMAX.GT.IA.OR.IMSUM.GT.IV)
     1  STOP '(DIMENSION ERROR: IA OR IV TOO SMALL)'
C       A.Kramida 2026: FIND THE TRAILER BLOCKS THAT BELONG TO THE CSET
C       JUST READ AND TAKE THEIR TRUE jj LABELS INTO YIDEN(.,3,.).  WITH
C       NCDES=0 THE CSET IN MEMORY IS THE ONE MAPPED LAST, SO THE MAP
C       STANDS.
      IF (IJJON.NE.0.AND.NCDES.NE.0) CALL MAPJJ2
C
C       NF(K,M) - ARRAY OF LEVEL FLAGS FOR MOMENTUM J=M:
C         NF(K,M)=0 FOR KNOWN LEVELS
C         NF(K,M)=-K FOR UNKNOWN LEVELS
C       T(K,M)  - ARRAY OF INPUT LEVEL ENERGIES
C       W(K,M) - ARRAY OF CALCULATED LEVELS
C       TW(K,M) - DIFFERENCES T-W FOR KNOWN LEVELS, 0 FOR UNKNOWN
C       LF(L) - ARRAY OF PARAMETER FLAGS
C       LPG(L) - Array of counts for parameter groups
C
      IF (INVECT.EQ.0) THEN             ! A.KRAMIDA
        READ (10,15) (LF(L), L=1,LMAX)
        READ (10,14) (XMAX(L), L=1,LMAX)
   14   FORMAT (7F10.5)
        READ (10,14) (X(L), L=1,LMAX)
      ELSE
C ****************************************** A.KRAMIDA v
        READ (20,3411) L
 3411   FORMAT (A1)
        IF (IDENOM.EQ.0) THEN
          DO L=1,LMAX
            READ (20,3410) XIDEN(L),LF(L),X(L),XMAX(L)
          ENDDO
 3410     FORMAT (A10,I5,2F14.4)
        ELSE
          DO L=1,LMAX
            READ (20,3430) XIDEN(L),LF(L),X(L),XMAX(L),IDENOM
          ENDDO
 3430     FORMAT (A10,3F14.4)
        ENDIF
      ENDIF
      I=0
      DO L=1,LMAX
        IF (XIDEN(L)(1:3).EQ.'EAV') I=I+1
        IF (XIDEN(L)(1:1).EQ.'1') I=0
        IF (I.NE.0) THEN
          WRITE (NX(L),3444) I
 3444     FORMAT (I2)
        ELSE
          WRITE (NX(L),3445) '  '
 3445     FORMAT (A2)
        ENDIF
      ENDDO
      IF (LMAX.EQ.1.AND.X(1).EQ.0.D0) X(1)=1.D-6  ! TO AVOID 1/0 IN EIGEN
C ****************************************** A.KRAMIDA ^
      NOCYC=0
      LSQM=0
      JLMAX=LMAX
      SUMXSQ=0.D0
      IVARP=0
      DO 116 L=1,LMAX
        LF0(L)=LF(L)
        XMAX0(L)=XMAX(L)
        SDX(L) = 0.D0
        X0(L)=0.D0
        DELXC(L)=0.D0
        DELX(L)=0.D0
        IF (XMAX(L).LE.0) THEN
          XMAX(L)=1.0D10
        ELSE
          XM=XMAX(L)            ! CONVERT REAL*4 TO REAL*8
          XMAX(L)=DMAX1(XM,X(L))
        ENDIF
        ISUM=0
        IF (LF(L).LT.-50) THEN
          XMIN(L)=-XMAX(L)
          IF (LF(L).GT.-99) THEN
            LF(L)=LF(L)+50
            LFL=LF(L)           ! CONVERT INTEREG*2 TO INTEGER*4
            JLMAX=MIN0(JLMAX,LFL)
c            IVARP=IVARP+1
            ISUM=1
          ELSE
            IF (LF(L).EQ.-99) THEN
              LF(L)=0
c              IVARP=IVARP+1
              ISUM=1
            ELSE
              GO TO 116
            ENDIF
          ENDIF
        ENDIF
        IF (ISUM.EQ.0) THEN
          XMIN(L)=DMAX1(1.0D-6, X(L)-XMAX(L)+X(L) )
          IF (LF(L).GE.100) GO TO 116
          IF (LF(L).LT.0) THEN
            LFL=LF(L)           ! CONVERT INTEREG*2 TO INTEGER*4
            JLMAX=MIN0(JLMAX,LFL)
          ELSE
            LF(L)=0
          ENDIF
          ISUM=1
c          IVARP=IVARP+1
        ENDIF
        IF (ISUM.NE.0) SUMXSQ=SUMXSQ+X(L)*X(L)
  116 CONTINUE
C**********************************************    A.KRAMIDA v
c     COUNT VARIABLE PARAMETERS
C
C     A. Kramida 2026: THE GROUPS ARE RENUMBERED HERE BEFORE THEY ARE
C       COUNTED.  REASON:
C
C       ON ENTRY THE DO 116 LOOP ABOVE HAS LEFT LF(L)=0 FOR A FREELY
C       VARIED PARAMETER, LF(L)=-G FOR A MEMBER OF LINKED GROUP G, AND
C       IABS(LF(L)).GE.100 FOR A FIXED ONE.  G IS THE GROUP NUMBER THE
C       USER PUNCHED, EITHER AS -G (G=1..50) OR AS -(G+50), AND NOTHING
C       OBLIGES THE NUMBERING TO BE CONSECUTIVE - GAPS ARE COMMON, SINCE
C       A GROUP DROPPED WHILE A FIT IS BEING DEVELOPED LEAVES A HOLE.
C       DO 116 NEVERTHELESS TOOK THE LARGEST GROUP NUMBER AS THE NUMBER
C       OF GROUPS (JLMAX=MIN0(JLMAX,LFL)), SO EACH HOLE ADDED A COLUMN
C       TO THE LEAST-SQUARES SYSTEM THAT NO PARAMETER FEEDS.  THE V*
C       CONDENSATION IN CALCV LEAVES THAT COLUMN ZERO, B(L,L)=0 IS THEN
C       FORCED TO 1 THERE AND D(L,1) IS ZERO, SO DELX COMES OUT ZERO AND
C       THE FITTED PARAMETERS WERE NOT AFFECTED.  PNO WAS: IT IS SET
C       FROM JLMAX AND IS THUS TOO LARGE BY THE NUMBER OF HOLES, AND PNO
C       CARRIES THE DEGREES OF FREEDOM IN STDEV=DSQRT(R/(TNO-PNO)),
C       WHICH SCALES EVERY SDX.  MEASURED ON U90: RENUMBERING THE ONE
C       GROUP FROM 1 TO 2 - THE SAME FIT, TO THE SAME AVDEV 191.931 AND
C       THE SAME DET B 1.573E+51 - RAISED THE PRINTED PARAMETER COUNT
C       FROM 7 TO 8 AND STDEV FROM 350.416 TO 429.170.
C
C       RENUMBERING THE GROUPS 1,2,3,... IN ORDER OF FIRST APPEARANCE
C       REMOVES THE HOLES AT SOURCE, SO JLMAX, PNO AND IVARP ARE ALL
C       COUNTS OF THINGS THAT EXIST.  LPG IS BUILT FROM LF RATHER THAN
C       FROM LF0, SO IT CANNOT DISAGREE WITH WHAT DO 116 ACTUALLY DID,
C       AND IT IS CLEARED OVER ITS WHOLE LENGTH BECAUSE A GROUP NUMBER
C       MAY EXCEED LMAX AND, WITH -fno-automatic, IT WOULD OTHERWISE
C       STILL HOLD THE PREVIOUS CONFIGURATION SET'S COUNTS.
C
C       IVARP IS THEN THE FREELY VARIED PARAMETERS PLUS THE NUMBER OF
C       LINKED GROUPS - ONE DEGREE OF FREEDOM PER GROUP, NOT ONE PER
C       MEMBER - WHICH IS WHAT THE GUARDS BELOW NEED.
      DO L=1,IB
        LPG(L) = 0
      ENDDO
      IPG = 0
      DO L=1,LMAX
        LFL=LF(L)
        IF (LFL.LT.0.AND.LFL.GT.-100) THEN
          IF (LPG(-LFL).EQ.0) THEN
            IPG = IPG+1
            LPG(-LFL) = IPG        ! OLD GROUP -LFL -> NEW GROUP IPG
          ENDIF
        ENDIF
      ENDDO
      IVARP=IPG
      DO L=1,LMAX
        LFL=LF(L)
        IF (LFL.LT.0.AND.LFL.GT.-100) THEN
          LF(L) = -LPG(-LFL)
        ELSE
          IF (IABS(LFL).LT.100) IVARP=IVARP+1   ! FREELY VARIED
        ENDIF
      ENDDO
      IF (IPG.GT.0) JLMAX=-IPG
C     INOCHG IS LEFT ON ITS ORIGINAL TEST, WHICH IS ON LF0 AND NOT ON
C     LF, SO THAT THE "NOTHING CHANGED IN OUTGINE" DETECTION IS
C     UNAFFECTED BY ANY OF THE ABOVE.
      DO L=1,LMAX
        IF (ABS(LF0(L)).NE.100) INOCHG=0
      ENDDO
C           IF NO PARAMETERS ARE VARIED,
C           OR NO LEVELS ARE TO FIT, OR
C           LEVELS ARE TO BE FIT WITH TOO MANY FREE PARAMETERS,
C           OR IT SEEMS THAT NOTHING WAS CHANGED IN OUTGINE FILE
C           AFTER RUNNING RCG,
C           THEN EXIT AFTER THE 0-TH CYCLE
C        IVARP = NUMBER OF VARIED PARAMETERS
C        IFIXL  = NUMBER OF UNKNOWN LEVELS (NOT "EXPERIMENTAL")
      IFIRST=0
c     A. Kramida 1/3/2025 change: 
c     Remove the condition that levels are to be fit with too many parameters
c     Replace it with condition that there are no fixed levels
c     IF (IVARP.GT.IMSUM-IFIXL) IFIRST=1
C     A. Kramida 2026: THE 1/3/2025 REPLACEMENT READ
C
C       IF (IFIXL.EQ.0) IFIRST=1
C
C       AND IS REMOVED HERE.  IFIXL COUNTS THE UNKNOWN LEVELS, SO THAT
C       TEST ABANDONED THE FIT WHENEVER EVERY LEVEL CARRIED AN
C       EXPERIMENTAL ENERGY - THE BEST-DETERMINED CASE THERE IS, AND
C       THE USUAL ONE FOR A SMALL ION.  THE INTENT RECORDED IN THE
C       COMMENT ABOVE, "NO FIXED LEVELS", IS THE OPPOSITE CONDITION -
C       NO KNOWN LEVEL TO FIT TO - AND THAT IS ALREADY TESTED BELOW AS
C       IFIXL.EQ.IMSUM, WHICH NOW WORKS BECAUSE IFIXL IS A TOTAL OVER
C       ALL J (SEE THE NOTE AT THE DO 105 LOOP).  SYMPTOM: RCE PRINTED
C       ITERATION 0 AND STOPPED, WITH EVERY LEVEL STARRED IN 'RCEOUT'
C       (LABEL 1342 STARS EVERYTHING WHEN IFIRST IS SET), SO FEEDING
C       'RCEOUT' BACK AS 'RCEINP' ALSO LOST THE EXPERIMENTAL ENERGIES.
C       THE DELIBERATE PART OF THE 2025 CHANGE - DROPPING COWAN'S
C       TOO-MANY-PARAMETERS TEST ON THE LINE ABOVE - IS LEFT AS IT IS.
      IF (IVARP.GE.LMAX.AND.IFIXL.EQ.0.AND.
     1  ((NOCYCE .EQ. 5 .AND. INOCHG .EQ. 1) .OR.
     2    NOCYCE .EQ. 0) ) IFIRST=1
      IF (IVARP.EQ.0.OR.IFIXL.EQ.IMSUM) IFIRST=1
      IF (IFIRST.NE.0) NOCYCE=0
C**********************************************    A.KRAMIDA ^
 1114 MIK=0
      IF (SUMXSQ.EQ.0.D0) THEN
C       REWIND I216
C       DO 1115 I=1,KKLM
C         DO 1114 K=1,KKLM
C1114     Y(I,K)=0.D0
C1115   Y(I,I)=1.D0
C       DO 1116 M=1,MMAX
C         IM=IMAX(M)
C         IMSQ=IM*IM
C         CALL WR(I216,IM,Y,IA,C,IMSQ)
C1116   CONTINUE
        DO M=1,MMAX
          IM=IMAX(M)
          DO K=1,IM
            DO I=1,IM
              MIK=MIK+1
c             IF (MIK.LE.MAXBY) THEN
c               IF (I.NE.K) THEN
c                 BUFY(MIK)=0.D0
c               ELSE
c                 BUFY(MIK)=1.D0
c               ENDIF
c             ENDIF
              if (MIK.GT.IBY_SIZE) then
*               Double the allocation of storage for the arrays
*
                NEWSIZ = 2*IBY_SIZE
                CALL RALLOC(BUFY_PTR,IBY_SIZE,NEWSIZ,8)
                IBY_SIZE = NEWSIZ
c                GOTO 1114
              endif
              IF (I.NE.K) THEN
                BUFY(MIK)=0.D0
              ELSE
                BUFY(MIK)=1.D0
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ELSE
        DO M=1,MMAX
          IM=IMAX(M)
          DO K=1,IM
            DO I=1,IM
              MIK=MIK+1
            ENDDO
          ENDDO
        ENDDO
      ENDIF
      WRITE (IW,1117) MIK,MAXBY
 1117 FORMAT (' Y-BUFFER SIZE=',I11,I9)
c     IF (MIK.GT.MAXBY) THEN
c       WRITE (IW,1119) MIK,MAXBY
c1119   FORMAT (15H Y buffer size=,I10, 12H, dimension=,I10)
c       STOP '(DIMENSION OF Y-BUFFER EXCEEDED.)'
c     ENDIF
C
C
C          READ ITERATION CONTROL CARD AND FIND STARTING POINT
C
  120 IF (INVECT.EQ.0) THEN    ! A.KRAMIDA
        READ (10,19) F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2 !R.COWAN
   19   FORMAT (3F10.1,F5.2,I5,F10.3,F10.2,2F5.2,I2)            !R.COWAN
C *************************************************** A.KRAMIDA v
        WRITE (20,19) F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2
      ENDIF
      WRITE (21,19) F0,FM,F1,DELF,NOFMAX,GCOE,DXMAX,CRIT1,CRIT,ICRIT2
C *************************************************** A.KRAMIDA ^
 2120 F=F0
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME
C
C        SET UP FOR CALCN WITH LINKED PARAMS
C
      IF (JLMAX.LT.0.AND.SUMXSQ.GT.0) THEN
        F=1.D0
        LSQM=-1
        JLMAX=-JLMAX
        DO L=1,LMAX
          LFL=LF(L)
          IABSLF=IABS(LFL)-100
          IF ((IABSLF.LT.0.AND.LF(L).GE.0).OR.
     *        (IABSLF.GE.0.AND.XMAX(L).LE.0.9D10)) THEN
            JLMAX=JLMAX+1
            LF(L)=-JLMAX
          ENDIF
        ENDDO
        PNO=JLMAX
      ENDIF
C
      FF=-1.D0
      G=1.D0
      IF (F.GE.F1) THEN
        FF=0.D0
        IF (F.GE.0.9999D0) THEN
          NOFMAX=1000
          FF=1.D0
          IF (LSQM.GE.0) THEN
            DO L=1,LMAX
              LFL=LF(L)
              IF (IABS(LFL).GE.100) THEN
                IF (XMAX(L).LE.0.9D10) LF(L)=0
              ELSE
                XMAX(L)=1.0D10
                IF (XMIN(L).GE.0) THEN
                  XMIN(L)=1.0D-6
                ELSE
                  XMIN(L)=-1.0D10
                ENDIF
              ENDIF
            ENDDO
          ENDIF
        ENDIF
      ENDIF
C
  400 IF (LSQM.GE.0) THEN
        JLMAX=0
        PNO=0.D0
        DO L=1,LMAX
          LFL=LF(L)
          IF (IABS(LFL).LT.100) THEN
            JLMAX=JLMAX+1
            LF(L)=JLMAX
            PNO=PNO+1.D0
          ENDIF
        ENDDO
      ENDIF
C
C     REWIND 3
C
      NOF=0
C     N5=-1
C ****************************************v  A.KRAMIDA
      DO M=1,MMAX
        DO K=1,IMMAX
          DO I=1,MAXC
            IEVV0(I,K,M)=IEVV(I,K,M)
            NW0(I,K,M)=NW(I,K,M)
          ENDDO
        ENDDO
      ENDDO
C ****************************************^  A.KRAMIDA
      IF (SUMXSQ.NE.0.OR.NOCYCE.EQ.0) THEN
        SUMXSQ=0.D0
        DETB=0.D0
C       N5=-2
C       IF (NOCYCE.LE.2) N5=1-NOCYCE
C       IF (NOCYCE.EQ.0) N5=10
        ABDXMX=0.D0
        CALL CLOCK(TIME)
        DELT=TIME-TIME0
        WRITE (IW,7) DELT,TIME
        GO TO 180
      ENDIF
C
C
C                    ITERATION CONTROL
C
C
  140 NOCYC=NOCYC+1
      NOF=NOF+1
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME
      IF (NOF.LE.NOFMAX) GO TO 170
      NOF=1
      IF (FF.GT.0) GO TO 170
      IF (FF.LT.0) THEN
        F=F*FM
        IF (F.LT.F1) GO TO 170
        F=F1-DELF
        FF=0.D0
      ENDIF
      IF (F.LT.0.9999D0) THEN
        F=DMIN1((F+DELF),1.0D0)
        IF (F.GT.0.9999D0) THEN
          F=1.D0
          NOFMAX=NOFMAX+4
        ENDIF
        GO TO 170
      ENDIF
      NOFMAX=1000
      FF=1.D0
      DO L=1,LMAX
        LFL=LF(L)
        IF (IABS(LFL).GE.100) THEN
          IF (XMAX(L).LE.0.9D10) THEN
            JLMAX=JLMAX+1
            IF (LSQM.LT.0) THEN
              LF(L)=-JLMAX
            ELSE
              LF(L)=JLMAX
            ENDIF
  158       PNO=PNO+1.D0
          ENDIF
        ELSE
          XMAX(L)=1.0D10
          IF (XMIN(L).LT.0) THEN
            XMIN(L)=-1.0D10
          ELSE
            XMIN(L)=1.0D-6
          ENDIF
        ENDIF
      ENDDO
C
C
C                    CALCULATE NEW PARAMETERS AND EIGENVALUES
C
  170 CALL CALCV                ! READS Y ARRAY FROM Y-BUFFER
                                ! AND C ARRAY FROM C-BUFFER
      DO M=1,MMAX
        DO K=1,IMMAX
          DO I=1,MAXC
            IEVV0(I,K,M)=IEVV(I,K,M)
            NW0(I,K,M)=NW(I,K,M)
          ENDDO
        ENDDO
      ENDDO
      CALL CALCX
C
  180 NSQ=0
      CALL EIGEN(kcpl,IBY_SIZE)        ! WRITES Y ARRAY TO Y-BUFFER
                                       ! READS C ARRAY FROM C-BUFFER
C
      DO L=1,LMAX
        SDX(L)=STDEV*SDX(L)
      ENDDO
C
C
C
C                    OUTPUT
C
      CALL CLOCK(TIME)
      DELT=TIME-TIME0
      WRITE (IW,7) DELT,TIME,TNO,PNO
C
      LEAVE=0                             ! SET END-OF-CYCLE CONDITION
      IF (NOCYC.GE.NOCYCE) THEN
        LEAVE=1
      ELSE
        IF (NOCYC.NE.0.AND.F.GT.0.9999D0) THEN
          ABDXMX=0.D0
          DO L=1,LMAX
            Z=DABS(DELX(L))
            IF (Z.GT.ABDXMX) ABDXMX=Z
          ENDDO
          IF (ABDXMX.LE.3.D-3) LEAVE=1
        ENDIF
      ENDIF
C
C     IF (IPRNX.EQ.0.OR.LEAVE) THEN
C            WRITE NEW AND OLD PARAMETERS TO 'OUTE'
      WRITE (IW,71) ELIDEN, (CONFIG(L,K3),L=1,18),NOCYC,F,GEFF,CRIT,
     1  DETB,AVDEV,STDEV  !, SD1,SD2
   71 FORMAT (1H , 16A8/9X,3A8/ 8H CYC NO=,I2, 4H  F=,F5.3,
     1  4H  G=,F5.3, 7H  CRIT=,F4.2, 8H  DET B=,1PE10.3,
     2  8H  AVDEV=,0PF8.3, 8H  STDEV=,F8.3/)
C     WRITE (IW,18)
C  18 FORMAT (1H0)
C     WRITE (IW,20)
   20 FORMAT (1H )
      WRITE (IW,72) NOCSET
   72 FORMAT (44H SLATER PRMTR    LF      NEW X         DEL X,
     1  39H        DEL XC      OLD X           SDX,
     2  34H         RSS W      DIFF     CSET=,I3)
C     DO L=1,LMAX
C
C
C       WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
C    1    DELXC(L), X0(L),   SDX(L)
C     ENDDO
      WRITE (IW,20)
C
      DO L=1,LMAX
C
C
        IF (L.LE.MMAX) THEN
          DIFF=ANORM(L)-SEVSQ(L)
          WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1      DELXC(L), X0(L),   SDX(L), SEVSQ(L), DIFF
   73     FORMAT (1H ,A10,2X,I6, 4F13.4, F13.3, F15.3, F9.3)
        ELSE
          WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1      DELXC(L), X0(L),   SDX(L)
        ENDIF
      ENDDO
      IF (IDENOM.NE.0.AND.LEAVE.NE.0) THEN
C
C       CONVERT SUPERSCRIPTED PARAMETERS TO SUBSCRIPTED
C
        WRITE (IW,71) ELIDEN, (CONFIG(L,K3),L=1,15),NOCYC,F,GEFF,
     1    CRIT,DETB,AVDEV,STDEV  !, SD1,SD2
        WRITE (IW,80)
   80   FORMAT (13H  SUBSCRIPTED )
        WRITE (IW,72) NOCSET
        DO L=1,LMAX
          IF (DENOM(L).NE.0.D0) THEN
            X(L) = X(L)/DENOM(L)
            DELX(L) = DELX(L)/DENOM(L)
            DELXC(L) = DELXC(L)/DENOM(L)
            X0(L) = X0(L)/DENOM(L)
            SDX(L) = SDX(L)/DENOM(L)
          ENDIF
C
C
          IF (L.LE.MMAX) THEN
            DIFF=ANORM(L)-SEVSQ(L)
            WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1        DELXC(L), X0(L),   SDX(L), SEVSQ(L), DIFF
          ELSE
            WRITE (IW,73) XIDEN(L),LF(L),X(L),DELX(L),
     1        DELXC(L), X0(L),   SDX(L)
          ENDIF
        ENDDO
      ENDIF
      IF (IPRNV.EQ.0.AND.NOCYCE.EQ.0) CALL CALCV  ! FOR V-PRINTING
C
      WRITE (11,39) NOCYC,AVDEV    ! WRITE TO 'PARVALS'
   39 FORMAT (32H PARAMETER VALUES,  CYCLE NUMBER,I3,9H,  AVDEV=,F8.3)
      WRITE (11,40) (X(L), L=1,LMAX)
C
C *************************************************** A.KRAMIDA v
      WRITE (IW,9711) NOCYC,AVDEV   ! WRITE TO OUTE
      IF (IW6.LT.0) WRITE (6,9711) NOCYC,AVDEV   ! WRITE TO THE SCREEN
9711  FORMAT (14H Iteration No ,I1,21H finished with AVDEV=,0PF8.3)
   40 FORMAT (7F10.4)
C
      IF (LEAVE.EQ.0) GO TO 140
C
C     IF (LEAVE.NE.0.OR.IPRNSQ.NE.2) THEN
C
C       REORDER EIGENVECTORS BY CALCULATED ENERGY
C
C       A.Kramida 2026: SEE ESORTW - THE SEQUENTIAL EXTRACTION OF MINIMA
C         THAT USED TO STAND HERE COST O(IM**3) AND DOMINATED THE
C         POST-ITERATION OUTPUT PHASE FOR LARGE BLOCKS.
        DO M=1,MMAX
          CALL ESORTW(IMAX(M),KW(1,M),W(1,M))
        ENDDO
C ***************************************** A.Kramida 2026 v
C       THE LEADING CONFIGURATION OF EVERY LEVEL, WORKED OUT HERE ONCE,
C       FROM THE EIGENVECTORS AS THE DIAGONALIZATION LEFT THEM.
C         EACH OUTPUT PASS USED TO WORK IT OUT FOR ITSELF, FROM BUFY AS
C       THAT PASS FOUND IT.  IN THE FIRST PASS BUFY IS THE EIGENVECTOR
C       ITSELF, BUT EVERY LATER PASS SEES THE VECTOR ORDOUT2 REBUILT
C       FROM IEVV AND NW - THE MAXC LARGEST COMPONENTS ONLY, EACH
C       ROUNDED TO 1.D-4 - SO THE SUMS OVER THE CONFIGURATIONS CAME OUT
C       SHORT, AND SHORT BY DIFFERENT AMOUNTS.  A LEVEL WHOSE TWO
C       LEADING CONFIGURATIONS RUN CLOSE COULD THEREFORE BE NAMED AFTER
C       ONE OF THEM IN ONE COUPLING SCHEME AND AFTER THE OTHER IN THE
C       NEXT, AND SINCE 'LEVELS3' IS GROUPED BY THAT NAME, EVEN THE
C       NUMBER OF LEVELS PER CONFIGURATION DEPENDED ON THE SCHEME.  THAT
C       LEAVES THE CONFIGURATION PURITIES OF TWO SCHEMES NOT COMPARABLE,
C       WHICH IS WHAT THOSE LISTS ARE READ FOR.
C         THE NAME CANNOT DEPEND ON THE SCHEME.  A CHANGE OF COUPLING
C       RECOUPLES ANGULAR MOMENTA WITHIN A CONFIGURATION AND CANNOT
C       CARRY ONE CONFIGURATION INTO ANOTHER, SO ITS MATRIX IS
C       BLOCK-DIAGONAL IN THE CONFIGURATION - VERIFIED ON THIS CODE'S
C       OWN MATRICES, WHERE THE LARGEST ELEMENT OUTSIDE THE BLOCKS IS
C       EXACTLY ZERO - AND THE SUM OF THE SQUARES OVER A CONFIGURATION
C       IS AN INVARIANT OF THE LEVEL.  SO IT IS FOUND ONCE, BEFORE ANY
C       PASS TRANSFORMS ANYTHING, FROM THE FULL-PRECISION COEFFICIENTS
C       AND ALL OF THEM, AND THE SAME VALUE SERVES EVERY PASS.
        NYOFS=0
        DO MZ=1,MMAX
          IMZ=IMAX(MZ)
          DO JZ=1,IMZ
            DO KK=1,IC
              VC(KK)=0.D0
            ENDDO
            JOFS=NYOFS+(JZ-1)*IMZ
            DO I=1,IMZ
              NCC=NC(I,MZ)
              VC(NCC)=VC(NCC)+BUFY(JOFS+I)**2
            ENDDO
            GRE=0.D0
            NJM=1
            DO KK=1,IC
              IF (VC(KK).GT.GRE) THEN
                GRE=VC(KK)
                NJM=KK
              ENDIF
            ENDDO
            NJMS(JZ,MZ)=NJM
          ENDDO
          NYOFS=NYOFS+IMZ*IMZ
        ENDDO
C ***************************************** A.Kramida 2026 ^
C
C               WRITE TO 'OUTE' AND 'LEVELS1' AND 'RCEOUT'/'RCEINP'
C
C ***************************************** A.Kramida 2026 v
C       A THIRD PASS IN TRUE jj FOLLOWS THE TWO COUPLINGS RCG WROTE ON
C       TAPE2E, WHENEVER THE FILE CARRIES A TRUE jj TRAILER AND THE
C       BLOCKS OF THIS CSET WERE FOUND IN IT.  IT ADDS TO THE OUTPUT AND
C       CHANGES NONE OF IT: THE TWO PASSES BEFORE IT RUN EXACTLY AS THEY
C       DID, AND THE LEVEL LISTS OF THE THIRD GO TO 'LEVELS1J' RATHER
C       THAN 'LEVELS1', SO THAT SORT8TO7 STILL SEES THE FOUR SECTIONS IT
C       EXPECTS AND 'LEVELS2' AND 'LEVELS3' COME OUT UNTOUCHED.
C         THE EIGENVECTORS ARE THE ONES THE FIT PRODUCED, RE-EXPRESSED
C       IN THE TRUE jj BASIS BY THE TRAILER'S TRANSFORMATION MATRIX, THE
C       SAME WAY COUPLING 2 IS REACHED THROUGH THE MATRIX ON TAPE2E.
C       THAT MATRIX IS <LS|TRUEjj>, SO IT ACTS ON VECTORS IN THE LS
C       BASIS; IF RCG'S FIRST COUPLING WAS NOT LS THE MATRICES ON TAPE2E
C       ARE IN SOME OTHER BASIS ALREADY AND THE PASS IS SKIPPED.
        NOUTMX=2
        IF (IJJON.NE.0.AND.IJJCS.NE.0) THEN
          IF (KCPL(1).EQ.1) THEN
            NOUTMX=3
          ELSE
            WRITE (IW,7901) KCPL(1)
            IF (IW6.LT.0) WRITE (6,7901) KCPL(1)
 7901       FORMAT (' TRUE jj OUTPUT SKIPPED: RCG WROTE COUPLING',I2,
     1        ' FIRST ON TAPE2E, NOT LS')
          ENDIF
        ENDIF
        DO NOUT=1,NOUTMX
          IU8=8
          IF (NOUT.EQ.3) THEN
            IU8=13
            WRITE (IW,7101) NBKJ2
 7101       FORMAT (/' TRUE jj COUPLING, FROM THE TAPE2E TRAILER',
     1        ' (',I4,' J BLOCKS).  THE KEY TO THE LABELS IS IN THE'/
     2        ' FILE true-jj_keys.txt WRITTEN BY THE RCG RUN THAT MADE',
     3        ' THIS TAPE2E.'/
     4        ' LEVEL LISTS GO TO LEVELS1J, LEVELS2J, LEVELS3J.')
          ENDIF
C ***************************************** A.Kramida 2026 ^
          DO NSQ1=0,1
            IF (NSQ1.EQ.1) THEN
C                   WRITE TO 'LEVELS1', FOR BOTH COUPLING SCHEMES
              WRITE (IU8,77) ELIDEN, (CONFIG(I,K3),I=1,NMAX)
   77         FORMAT (A8,10(2X,A8)/(8X,10(2X,A8)))
            ENDIF
            IF (IPRNSQ.NE.2) THEN
C
C               PRINT EIGENVECTOR COMPONENTS FOR BOTH COUPLING SCHEMES
C
              IF (NSQ1.EQ.0) THEN
                WRITE (IW,74) NORUN,NOCYC
   74           FORMAT (/1H ,A2,I3,
     1            '  T         W        T-W    GKK    G    CONF',
     2            '     EIGENVECTOR COMPONENTS')
              ELSE
                IF (IPRNSQ.EQ.0) THEN
                  WRITE (IW,75) NORUN,NOCYC
   75             FORMAT (/1H ,A2,I3,
     1              '  T         W        T-W    GKK    G    CONF',
     2              '     EIGENVECTOR COMPONENTS SQUARED')
                ENDIF
              ENDIF
            ENDIF
C           MIK=0
            NYOFS=0
            DO M=1,MMAX
              IM=IMAX(M)
              IMSQ=IM*IM
              ITRANS=0
              IF (NOUT.EQ.2.AND.KCPL(1).EQ.1) ITRANS=1
              IF (NOUT.EQ.1.AND.KCPL(1).EQ.2) ITRANS=1
              IF (KCPL(2).EQ.1.AND.KCPL(1).NE.2) THEN
                IF (NOUT.EQ.1) THEN
                  ITRANS=2
                ELSE
                  ITRANS=1         ! A.Kramida 1/31/2011 (before if was ITRANs=0)
                ENDIF
              ENDIF
C ***************************************** A.Kramida 2026 v
C               THE THIRD PASS TAKES THE COUPLING-1 (LS) VECTORS AGAIN,
C               NOT THE COUPLING-2 ONES THE SECOND PASS LEFT BEHIND, SO
C               THEY ARE PUT ASIDE BEFORE THE SECOND PASS TRANSFORMS
C               THEM.  IEVV2/NW2 ARE FREE HERE: THEY BELONG TO THE
C               ITRANS=2 CASE ABOVE, WHICH NEEDS KCPL(2)=1 AND SO CANNOT
C               ARISE WHEN THE TRUE jj PASS RUNS.
              IF (NOUTMX.EQ.3.AND.NOUT.EQ.2.AND.NSQ1.EQ.0) THEN
                DO KKKK=1,IMMAX
                  DO IIII=1,MAXC
                    IEVV2(IIII,KKKK,M)=IEVV(IIII,KKKK,M)
                    NW2(IIII,KKKK,M)=NW(IIII,KKKK,M)
                  ENDDO
                ENDDO
              ENDIF
              IF (NOUT.EQ.3) THEN
                ITRANS=0
                IF (NSQ1.EQ.0) THEN
                  IMSQ1=(IMSQ+IM)/2
                  DO KKKK=1,IMMAX
                    DO IIII=1,MAXC
                      IEVV(IIII,KKKK,M)=IEVV2(IIII,KKKK,M)
                      NW(IIII,KKKK,M)=NW2(IIII,KKKK,M)
                    ENDDO
                  ENDDO
                  CALL RD1(M,LMAXP1,IM,C,IA,BUFC,IBC,IMSQ1,
     *              ICOFFS,MAXB1,MMAX,LMAXP1)
C                   THE TRUE jj MATRIX IN PLACE OF THE J1J2 ONE.  ORDOUT2
C                   THEN DOES EXACTLY WHAT IT DOES FOR COUPLING 2.
                  CALL RD1(MJJ(M),1,IM,A,IA,BUFJ,IBAJ,IMSQ,
     *              JOFFS,MAXBJU,NBJJX,1)
                  CALL ORDOUT2
                ENDIF
              ENDIF
C ***************************************** A.Kramida 2026 ^
              IF (ITRANS.NE.0.AND.NSQ1.EQ.0) THEN
c             IF (NOUT.EQ.2.AND.NSQ1.EQ.0) THEN
                IMSQ1=(IMSQ+IM)/2
C               CALL RD(I216,IM,Y,IA,EVEC,IMSQ)  ! R.COWAN
C               CALL RD(15,IM,C,IA,EVEC,IMSQ1)   ! R.COWAN
C               CALL RD(15,IM,A,IA,EVEC,IMSQ)    ! R.COWAN
C               DO K=1,IM
C                 DO I=1,IM
C                   MIK=MIK+1
C                   Y(I,K)=BUFY(MIK)
C                 ENDDO
C               ENDDO
                CALL RD1(M,LMAXP1,IM,C,IA,BUFC,IBC,IMSQ1,   ! A.KRAMIDA
     *            ICOFFS,MAXB1,MMAX,LMAXP1)
                CALL RD1(M,IA2,IM,A,IA,BUFA,IBA,IMSQ,       ! A.KRAMIDA
     *            IAOFFS,MAXB2,MMAX,IA2)
                DO K=1,ITRANS
c                                                           ! A.Kramida 1/31/2011 V
                  if (itrans.eq.2.AND.k.eq.1) then
c                     Save IEVV to IEVV2
                    DO KKKK=1,IMMAX
                      DO IIII=1,MAXC
                        IEVV2(IIII,KKKK,M)=IEVV(IIII,KKKK,M)
                        NW2(IIII,KKKK,M)=NW(IIII,KKKK,M)
                      enddo
                    enddo
                  else
                    if (itrans.eq.1.and.kcpl(1).ne.1) then
c                      Restore IEVV from IEVV2
                      DO KKKK=1,IMMAX
                        DO IIII=1,MAXC
                          IEVV(IIII,KKKK,M)=IEVV2(IIII,KKKK,M)
                          NW(IIII,KKKK,M)=NW2(IIII,KKKK,M)
                        enddo
                      enddo
                    endif
                  endif
c                                                           ! End A.Kramida 1/31/2011 ^
                  CALL ORDOUT2
                ENDDO
              ENDIF
C
C               Find pure-coupling-term labels for levels
              IF (CRIT1.GT.0.AND.NSQ1.EQ.1) CALL ORDOUT1
C
              DO K=1,IM
                KWKM=KW(K,M)
                J=IABS(KWKM)
                IF (W(J,M).LE.WMAX.OR.TW(J,M).NE.0.D0) THEN
C                   Configuration name for each level (by max.
C                   contribution).  A.Kramida 2026: taken from NJMS,
C                   which was filled from the untransformed
C                   eigenvectors before the first pass, so that every
C                   coupling scheme names the level the same way.  It
C                   used to be summed here, from a BUFY that the later
C                   passes had rebuilt from rounded and truncated
C                   components.
                  NJM=NJMS(J,M)
C
C                   Print the level to all output files needed
                  ERAS(1)=' '
                  IF (NF(J,M).LT.0.OR.IFIRST.NE.0) ERAS(1)='*'
                  IF (NF(J,M).LT.0) ERAS(1)='*'
C                 IF (IPRNSQ.EQ.0.OR.(NSQ1.EQ.0.AND.IPRNSQ.NE.2))
C    1              WRITE (IW,76) T(J,M),ERAS(1), W(J,M), TW(J,M),
C    2                GKK(J,M),GFAC(J,M), CONFIG(NJM,K3),
C    3                (EVEC(I,J), I=1,IM)
C  76               FORMAT (1H ,F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,
C    1                A6,8F10.3:/(49X,8F10.3:))
C                 IF (LEAVE.NE.0)
                  CALL OUT8(J,NJM,NOUT,NSQ1,LEAVE,IPRNSQ,KCPL(NOUT))
                ENDIF
              ENDDO
              NYOFS=NYOFS+IMSQ
            ENDDO
          ENDDO
          WRITE (IW,78) NOUT
   78     FORMAT (29H Finished output for coupling,I2)
        ENDDO
C     ENDIF
      IF (LEAVE.EQ.0) GO TO 140
C
C          WRITE PARAMETER VALUES FOR RCG INPUT
C
  900 WRITE (11,91)
   91 FORMAT (31H PARAMETER VALUES FOR RCG INPUT)
      L2=0
      DO NAZ=1,NMAX
        L1=L2+1
        L2=L1+NPAR(NAZ,NAZ)-1
        WRITE (11,92) ELIDEN,CONFIG(NAZ,K3), (X(L),L=L1,L2)
   92   FORMAT (A6,A8,6X,F10.4,4(F9.4,1X),10HLS99999999:/
     1    (7(F9.4,1X),10HLS99999999))
      ENDDO
      NXM1=NMAX-1
      IF (NXM1.EQ.0) GO TO 969
      DO N1=1,NXM1
        N1P1=N1+1
        DO N2=N1P1,NMAX
          NP12=NPAR(N1,N2)
          IF (NP12.NE.0) THEN
            L1=L2+1
            L2=L1+NP12-1
            WRITE (11,93) CONFIG(N1,K3),CONFIG(N2,K3), (X(L),L=L1,L2)
   93       FORMAT (A8,1H-,A8,3X,5(F9.4,1X),10HLS99999999:/
     1        (7(F9.4,1X),10HLS99999999))
          ENDIF
        ENDDO
      ENDDO
C *************************************************** A.KRAMIDA v
C
C          WRITE PARAMETER VALUES FOR RCE INPUT
C
  969 IF (INVECT.EQ.0) WRITE (20,3490)
      WRITE (21,3490)
 3490 FORMAT (51H PARAMETER FLAG      VALUE     MAX.VALUE      DENOM)
      LCI=0
      IF (IDENOM.NE.0) THEN
        DO L=1,LMAX
          X(L)=X(L)*DENOM(L)
        ENDDO
      ENDIF
      DO L=1,LMAX
        IF (XIDEN(L)(1:1).EQ.'1') LCI=1
        C1(1)=' '
        C1(2)=' '
        IF (LCI.NE.0) THEN
          READ (XIDEN(L),3505) III1,III2
c          write (6, 3506) XIDEN(L),II1,II2
 3505     FORMAT(2A1)
 3506     FORMAT(1H ,A8,2I3)
c           A. Kramida: decode alphanumeric codes of config. numbers into integer numbers
          NMAX2=NMAX
          II1=IDECCH(III1,NMAX2)
          II2=IDECCH(III2,NMAX2)
          C1(1)=CONFIG(II1,K3)
          C1(2)(1:1)='-'
          C1(2)(2:7)=CONFIG(II2,K3)
        ENDIF
        IF (IDENOM.NE.0) THEN
          IF (INVECT.EQ.0) WRITE (20, 3520) XIDEN(L),
     1      LF0(L),X(L),XMAX0(L),DENOM(L),C1(1),C1(2)
          WRITE (21, 3520) XIDEN(L),LF0(L),X(L),
     1      XMAX0(L),DENOM(L),C1(1),C1(2)
 3520     FORMAT (A10,I5,3F14.6,A6,A8)
        ELSE
          IF (INVECT.EQ.0) WRITE (20, 3540) XIDEN(L),
     1      LF0(L),X(L),XMAX0(L),C1(1),C1(2)
          WRITE (21, 3540) XIDEN(L),LF0(L),X(L),
     1      XMAX0(L),C1(1),C1(2)
 3540     FORMAT (A10,I5,2F14.6,14X,A6,A8)
        ENDIF
      ENDDO
C *************************************************** A.KRAMIDA ^
C
      WRITE (11,93)
C
C          SET FOR NEXT VALUE OF LSQM OR GO TO NEXT RUN
C
C 970 IF (LSQM.LT.0) THEN             ! A.KRAMIDA: ICRIT2 IS DISABLED
C       LSQM=LSQM+1                   ! BECAUSE THE EIGENVECTORS SAMPLE
C       IF (LSQM.LE.ICRIT2) GO TO 400 ! IS IN 1ST COUPLING SCHEME BUT
C     ENDIF                           ! CURRENT EIGENVECTORS HAVE BEEN
C                                     ! REARRANGED TO 2ND COUPLING SCHEME
C *************************************************** A.KRAMIDA v
      IF (INVECT.NE.0) THEN
 9709   READ (10,2110) ERAS(1)
 2110   FORMAT (A8)
        IF (ERAS(1).NE.'   -1   ') GO TO 9709
      ENDIF
C *************************************************** A.KRAMIDA ^
      GO TO 94
C
c     Deallocate memory
      CALL DALLOC(BUFY_PTR)
      CALL DALLOC(BUFC_PTR)
      CALL DALLOC(BUFA_PTR)
      CALL DALLOC(IBC_PTR)
      CALL DALLOC(IBA_PTR)
      CALL DALLOC(D_PTR)
      CALL DALLOC(V_PTR)
      CALL DALLOC(TW_PTR)
      CALL DALLOC(T_PTR)
      CALL DALLOC(EVAL_PTR)
      CALL DALLOC(W_PTR)
      CALL DALLOC(KW_PTR)
      CALL DALLOC(NF_PTR)
      CALL DALLOC(NC_PTR)
      CALL DALLOC(C_PTR)
      CALL DALLOC(A_PTR)
      CALL DALLOC(IEVV_PTR)
      CALL DALLOC(IEVV0_PTR)
      CALL DALLOC(NW_PTR)
      CALL DALLOC(NW0_PTR)
c
      END
C
************************************************************************
*                                                                      *
      SUBROUTINE ALLOC (PTR,NMLOCS,LENGTH)
*                                                                      *
*   Borrowed from GRASP2K [AK]                                         *
*   This  routine allocates  NMLOCS memory locations, each of length   *
*   LENGTH bytes to the pointee of the pointer PTR.   LENGTH is  not   *
*   required for the  CFT77  version --- the  8-byte word  length is   *
*   assumed common to all data types; it has been retained to ensure   *
*   compatibility with the code that calls this subprogram.            *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This version by Farid A. Parpia.      Last revision: 17 Sep 1992   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      *
************************************************************************
*
c
cbieron DEC pointer
c
c      INTEGER PTR                 ! DEC and CRAY require ?
      INTEGER*8 PTR
c      COMMON /mpi/ myid, nprocs, ierr
      pointer (PTR,ptrdummy)
      INTEGER NMLOCS, LENGTH
      INTEGER*8 NBYTES
      INTEGER*8 MALLOC
c      INTEGER*8 MMM
c
*
*   Compute the size of the memory region in bytes; this is retained for
*   calling compatibility and debugging purposes in the Cray version
*   Also, all allocations are multiples of 8 bytes for alignment
*   purposes.
*
cgd   NBYTES = NMLOCS*LENGTH
      NBYTES = NMLOCS
      NBYTES = NBYTES*LENGTH
      IF (MOD(NBYTES,8) .NE. 0) NBYTES = 8*(NBYTES/8 + 1)
*
      IF (NBYTES. LE. 0) THEN
         PRINT *, 'ALLOC: Invalid memory request:'
         PRINT *, ', PTR = ',PTR,', NMLOCS = ',NMLOCS,
     :            ', LENGTH = ',LENGTH,'.'
         STOP
*
      ELSE

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!         PTR = MMM (NBYTES)                ! SUN & DEC
!         PTR = MALLOC(%VAL(NBYTES))        ! IBM
!         IF (PTR .EQ. 0) THEN
!
!         IABORT = 0
!         CALL HPALLOC(PTR,NMLOCS,IERR,IABORT) ! CRAY
!         IF (ierr.ne.0) THEN
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

! Choose one from above

         PTR = MALLOC(NBYTES)                ! SUN & DEC
         IF (PTR .EQ. 0) THEN

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

            PRINT *, 'ALLOC: Unable to allocate memory:'
            PRINT *, ', PTR = ',PTR,', NMLOCS = ',NMLOCS,
     :               ', LENGTH = ',LENGTH,'.'
            STOP
         ENDIF
      ENDIF
      RETURN
      END
c------------------------------------
c
************************************************************************
*                                                                      *
      SUBROUTINE DALLOC (PTR)
*                                                                      *
*   Borrowed from GRASP2K [AK]                                         *
*   This  subprogram deallocates the  memory locations that are asso-  *
*   ciated  with the pointer PTR.   This version  is specific to Sun   *
*   FORTRAN 77  (Sun OS)  or  Cray CFT77 (UNICOS)  depending on  the   *
*   action of the preprocessor PREPRO.                                 *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This revision by Farid A. Parpia.     Last revision: 17 Sep 1992   *
*   Updated by Anders Ynnerman            Last revision: 26 Jan 1994   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      *
************************************************************************
*
c
cbieron DEC pointer
c
c      INTEGER PTR                        ! DEC and CRAY require ?
      integer*8 PTR
      POINTER (PTR,ptrdummy)

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!      CALL FREE (PTR)                    ! SUN & DEC
!      CALL FREE (%val(PTR))              ! IBM
!      CALL DISCLAIM (%val(PTR))          ! IBM
!      IERR = MZDALLOC(PTR,NSIZE)         ! TSS, Needs nsize 
!      IABORT = 0
!      CALL HPDEALLC(PTR,IERR,IABORT)     ! CRAY 
!      IF (IERR.NE.0) THEN
!         PRINT *, 'DEALLOCATION ERROR ',IERR
!      END IF
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      CALL FREE (PTR)                    ! SUN & DEC

!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

      RETURN
      END
c------------------------------------
c
************************************************************************
*                                                                      *
      SUBROUTINE RALLOC (PTR,OLDSIZ,NEWSIZ,LENGTH)
*                                                                      *
*   Borrowed from GRASP2K [AK]                                         *
*   This  subprogram reallocates  the memory for the pointee of PTR.   *
*   OLDSIZ is the number of storage locations of length LENGTH (4 or   *
*   8) bytes that are associated with  PTR upon entry. NEWSIZ is the   *
*   number of storage locations associated with  PTR  on exit.  This   *
*   version is specific to  SUN FORTRAN 77  (Sun  OS) or  Cray CFT77   *
*   (UNICOS)  depending on the action of the preprocessor PREPRO.      *
*                                                                      *
*   Call(s) to: [LIB92]: ALLOC, DALLOC.                                *
*                                                                      *
*   Original code by Charlotte F. Fischer.                             *
*                                                                      *
*   This revision by Farid A. Parpia.     Last revision: 17 Sep 1992   *
*   Updated for IBM and DEC by Bieron     Last revision: 05 Apr 1995   *
*                                                                      *
************************************************************************
*
      IMPLICIT REAL*8           (A-H,O-Z)
      INTEGER OLDSIZ
c      COMMON/iounit/istdi,istdo,istde
*
      INTEGER*8 PTR
      POINTER (PTR,ptrdummy)
*
      INTEGER*8 NBYTE
      INTEGER*4 IA,ITA
      CHARACTER*12 C,CT
      POINTER (IPTR,IA(1)),(ITPTR,ITA(1))
      POINTER (RPTR,RA(1)),(RTPTR,RTA(1))
      POINTER(CPTR,C(1)),(CTPTR,CT(1))
*
c      IF (NEWSIZ*LENGTH .LE. 0) THEN
      NBYTE=NEWSIZ
      NBYTE=NBYTE*LENGTH
      IF (NBYTE .LE. 0) THEN
         WRITE(*,*) 'RALLOC: Invalid memory request:'
         WRITE(*,*) ' NEWSIZ = ',NEWSIZ,', LENGTH = ',LENGTH,'.'
      STOP
      ENDIF
*
*   M is the number of array elements that are to be copied into the
*   new array
*
      M = MIN (OLDSIZ,NEWSIZ)
*
      IF (LENGTH .EQ. 8) THEN
*
*   Arrays of length 8 bytes are assumed to be REAL*8          ;
*   these are treated in the first branch
*
         RTPTR = PTR
         CALL ALLOC (RPTR,NEWSIZ,LENGTH)
*
         DO I = 1,M
            RA(I) = RTA(I)
         END DO
*
*   Reset PTR so that it is associated with the new array
*
         PTR = RPTR
*
*   Deallocate the array that used to be associated with PTR
*
         CALL DALLOC (RTPTR)
*
      ELSEIF (LENGTH .EQ. 12) THEN
*
*   Character*12 arrays are treated in the third branch; the procedure is
*   essentially the same as in the first branch
*
*   Associate CT with PTR
*
         CTPTR = PTR
*
*   Associate C with CPTR
*
         CALL ALLOC (CPTR,NEWSIZ,LENGTH)
*
*   Copy the contents of CT into IA
*
         DO I = 1,M
            C(I) = CT(I)
         END DO
*
*   Reset PTR so that it is associated with the new array
*
         PTR = CPTR
*
*   Release CT
*
         CALL DALLOC (CTPTR)
*
      ELSEIF (LENGTH .EQ. 4) THEN
*
*   Integer arrays are treated in the second branch; the procedure is
*   essentially the same as in the first branch
*
*   Associate ITA with PTR
*
         ITPTR = PTR
*
*   Associate IA with IPTR
*
         CALL ALLOC (IPTR,NEWSIZ,LENGTH)
*
*   Copy the contents of ITA into IA
*
         DO I = 1,M
            IA(I) = ITA(I)
         end do
*
*   Reset PTR so that it is associated with the new array
*
         PTR = IPTR
*
*   Release ITA
*
         CALL DALLOC (ITPTR)
*
      ELSE
*
         WRITE(*,*) 'RALLOC: Invalid argument LENGTH = ',LENGTH,'.'
         STOP
*
      ENDIF
*
      RETURN
      END
c------------------------------------
c
      SUBROUTINE UNPK2(LD,IM,A,IA,B_PTR,IMSQ,IOFFS,MAXLEN)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 B,BB
      INTEGER*8 B_PTR
      POINTER (B_PTR,B(1))
      INTEGER*2 BB1,I,J
      REAL*8 R8
      DIMENSION A(IA,IA),BB1(6)
      EQUIVALENCE (R8,BB1(1)),(BB,BB1(1)),(I,BB1(5)),(J,BB1(6))
C
      READ (LD) IOFFS2,(B(K),K=IOFFS+1,min(IOFFS+IOFFS2,maxlen))
      if (IOFFS+IOFFS2.gt.maxlen) then
c        write (6,*) "Attempted use of array B in UNPK2 beyond ",MAXLEN
c        stop '(Execution aborted.)'
        NEWSIZ=max(2*maxlen,ioffs+ioffs2)
        call RALLOC(B_PTR,maxlen,NEWSIZ,12)
        maxlen=NEWSIZ
        BACKSPACE LD
        READ (LD) IOFFS2,(B(K),K=IOFFS+1,IOFFS+IOFFS2)
      endif
C
      DO J=1,IM
        DO I=1,IM
          A(I,J)=0.D0
        ENDDO
      ENDDO
C
      IM2=IM*IM
      IF (IMSQ.LT.IM2) THEN
        DO K=IOFFS+1,IOFFS+IOFFS2
          BB=B(K)
          A(I,J)=R8
          A(J,I)=R8
        ENDDO
      ELSE
        DO K=IOFFS+1,IOFFS+IOFFS2
          BB=B(K)
          A(I,J)=R8
        ENDDO
      ENDIF
      IOFFS=IOFFS+IOFFS2
  300 RETURN
      END
C
C
      SUBROUTINE OUT8(J,NJM,NOUT,NSQ1,LEAVE,IPRNSQ,KCPL)
C        OUTPUT LEVELS LIST
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C  *******************************************************v A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC_PTR,BUFA_PTR,BUFY_PTR,
c     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
     1  IBC_PTR,IBA_PTR,ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
      INTEGER*8 BUFY_PTR,BUFA_PTR,BUFC_PTR,IBA_PTR,IBC_PTR
      POINTER (BUFY_PTR,BUFY(1)),(BUFA_PTR,BUFA(1)),(BUFC_PTR,BUFC(1))
      POINTER (IBA_PTR,IBA),(IBC_PTR,IBC)
      DIMENSION IBC(IX,IB),IBA(IX,1)
c
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0
      INTEGER*2 II1,II2
      COMMON/VECT/ IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      integer*8 IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      DIMENSION IEVV(MAXC,IA,IX),IEVV0(MAXC,IA,IX),NW(MAXC,IA,IX),
     1 NW0(MAXC,IA,IX)
      POINTER(IEVV_PTR,IEVV),(IEVV0_PTR,IEVV0),(NW_PTR,NW),(NW0_PTR,NW0)
c
C  *******************************************************^ A.KRAMIDA
C
      CHARACTER C1*8,Y1*12,Y11*3,AJ1*2,AJ2*2
      DIMENSION IP1(MAXC),IP2(MAXC),C1(MAXC),Y1(MAXC)
C *********************************************** R.COWAN v
C     DO 120 L=1,3
C       GRE=0.D0
C       DO 110 I=1,IM
C         IF (L.GT.1.AND.I.EQ.I1(1)) GO TO 110
C         IF (L.GT.2.AND.I.EQ.I1(2)) GO TO 110
C         IF (EVEC(I,J).LT.GRE) GO TO 110
C         GRE=EVEC(I,J)
C         I1(L)=I
C 110   CONTINUE
C 120 CONTINUE
C *********************************************** R.COWAN ^
C *********************************************** A.KRAMIDA v
C     IF (LEAVE.EQ.0.AND.NOUT.EQ.2) RETURN
C     IF (LEAVE.EQ.0.AND.IPRNSQ.NE.0.AND.NSQ1.EQ.1) RETURN
C     IF (LEAVE.EQ.0.AND.IPRNSQ.GT.1) RETURN
C
      IF (LEAVE.EQ.0) RETURN
C        IPRNSQ: 0 - PRINT AMPLITUDES AND SQUARES TO 'OUTE' : NSQ1=0,1
C                1 - PRINT ONLY AMPLITUDES TO 'OUTE'        : NSQ1=0
C                2 - NO EIGENVECTORS PRINT TO 'OUTE'
C                3-  PRINT TO 'OUTE' IN FORMAT READABLE BY TROUT
C     IF (NOUT.EQ.2.AND.NSQ1.EQ.0) RETURN
      L=MIN0(IM,MAXC,10)
      LM=0
      LM1=0
      LM2=0
      DO 128 I2=1,L
C       I=I1(I2)
C       IP1(I2)=100.D0*EVEC(I,J)
        I=NW(I2,J,M)
C       I=I+128
        IF (I.EQ.0) GO TO 129
        IP1(I2)=IEVV(I2,J,M)
C       P1=EVEC(I,J)
C       IP2(I2)=1.D2*P1+0.5D0
        IP2(I2)=DFLOAT(IP1(I2))*(IP1(I2)*1.D-6)+0.5D0
C       IP1(I2)=IP1(I2)*1.D-2+0.5D0
C ***************************************** A.Kramida 2026 v
C         IP2 IS A SQUARE AND SO NEVER NEGATIVE - ITS +0.5D0 IS RIGHT.
C         IP1 IS THE SIGNED PERCENTAGE AMPLITUDE, AND HERE THE FIXED
C         +0.5D0 TRUNCATED TOWARDS ZERO FOR EVERY NEGATIVE COMPONENT,
C         NOT JUST AT A ROUNDING BOUNDARY: IEVV=-9110 GAVE -90 WHERE
C         -91 IS CORRECT.  DSIGN GIVES THE OFFSET THE SIGN OF THE
C         AMPLITUDE.  POSITIVE AMPLITUDES KEEP THE IDENTICAL EXPRESSION.
        PAMP=IP1(I2)*1.D-2
        IP1(I2)=PAMP+DSIGN(0.5D0,PAMP)
C ***************************************** A.Kramida 2026 ^
        IF (IP1(I2).EQ.0) GO TO 129
        IF (I2.LE.3.AND.IP2(I2).GT.0) LM=I2
        IF (IP2(I2).NE.0) LM2=I2
        LM1=I2
        LC=NC(I,M)
        C1(I2)=CONFIG(LC,K3)
        Y1(I2)=YIDEN(I,NOUT,M)
        IF (KCPL.EQ.2.AND.NSQ1.EQ.1) THEN

C             UNPACK JJ-NOTATION ONLY FOR PRINTING TO 'LEVELS1..3'
C
          READ (Y1(I2),1101) Y11,J1,J2
 1101     FORMAT (A3,I2,1X,I2)
C          YJ1=DFLOAT(J1)/2
C          YJ2=DFLOAT(J2)/2
          IF (MOD(J1,2).EQ.0) THEN
            J1=J1/2
            AJ1='  '
          ELSE
            AJ1='/2'
          ENDIF
          IF (MOD(J2,2).EQ.0) THEN
            J2=J2/2
            AJ2='  '
          ELSE
            AJ2='/2'
          ENDIF
          WRITE (Y1(I2),1102) Y11,J1,AJ1,J2,AJ2
 1102     FORMAT (A3,I2,A2,1H),I2,A2)
C         WRITE (Y1(I2),1102) Y11,YJ1,YJ2
C1102     FORMAT (A3,F4.1,1H),F4.1)
        ENDIF
  128 CONTINUE
  129 TS=T(J,M)
      WS=W(J,M)
      TWS=TW(J,M)
C       A.Kramida 2026: THE TRUE jj PASS WRITES ITS LEVEL LIST TO
C       'LEVELS1J', LEAVING 'LEVELS1' AS IT ALWAYS WAS.
      IU8=8
      IF (NOUT.EQ.3) IU8=13
C     IF (LEAVE.NE.0.AND.NSQ1.EQ.1) THEN
      IF (NSQ1.EQ.1) THEN
C           WRITE TO 'LEVELS1'
        IF (S.LE.1.D0) THEN
          if (EMAX.LT.1000000.0) then  ! Condition added by A. Kramida 08/20/04
            WRITE (IU8,13) TS,ERAS(1),WS,TWS,GFAC(J,M),CONFIG(NJM,K3),
     1        FJ(M), (IP2(I),C1(I),Y1(I), I=1,LM)
          else
            WRITE (IU8,1301) TS,ERAS(1),WS,TWS,GFAC(J,M),CONFIG(NJM,K3),
     1        FJ(M), (IP2(I),C1(I),Y1(I), I=1,LM)
          endif
C             Mod A. Kramida 12/09/2009
C  13     FORMAT (F13.6,A1,F10.3,F11.6,F7.3,1X,A6,A8,4X,
   13     FORMAT (F13.6,A1,F11.4,F10.6,F7.3,1X,A6,A8,4X,
C             End mod A. Kramida 12/09/2009
     1      3(I4,2H% ,A6,1X,A12))
 1301     FORMAT (F13.4,A1,F10.1,F11.4,F7.3,1X,A6,A8,4X,
     1      3(I4,2H% ,A6,1X,A12))
        ELSE
            ! Print scaled energies
          TS1=S*TS
          WS1=S*WS
          TWS1=S*TWS
          WRITE (IU8,14) TS1,ERAS(1),WS1,TWS1,GFAC(J,M),CONFIG(NJM,K3),
     1      FJ(M), (IP2(I),C1(I),Y1(I), I=1,LM)
   14     FORMAT (F13.3,A1,F10.0,F11.3,F7.3,1X,A6,A8,4X,
     1      3(I4,2H% ,A6,1X,A12))
        ENDIF
      ENDIF
C
      IF (NOUT.EQ.1.AND.NSQ1.EQ.0) THEN
C       IF (LEAVE.NE.0) THEN
C            WRITE TO 'RCEOUT';
C            WRITE TO 'RCEINP' IF INVECT.EQ.0
          READ (FJ(M),15) FJ1
   15     FORMAT (4X,F4.1)
          J1=2.D0*FJ1
          IF (MOD(J1,2).EQ.0) THEN
            J1=J1/2
            AJ1='  '
          ELSE
            AJ1='/2'
          ENDIF
          IM1=MIN0(LM1,10)
          IF (INVECT.EQ.0) THEN
C            WRITE TO 'RCEINP'
            if (EMAX.LT.1000000.0) then  ! Condition added by A. Kramida 08/20/04
              WRITE (20,16)
     1          TS,ERAS(1),WS,J1,AJ1,(IP1(I),C1(I),Y1(I), I=1,IM1)
            else
              WRITE (20,1601)
     1          TS,ERAS(1),WS,J1,AJ1,(IP1(I),C1(I),Y1(I), I=1,IM1)
            endif
          ENDIF
   16     FORMAT (F13.6,A1,F10.3,I4,A2,10(I4,1H ,A6,1X,A12:))
 1601     FORMAT (F13.4,A1,F10.1,I4,A2,10(I4,1H ,A6,1X,A12:))
C            WRITE TO 'RCEOUT';
          if (EMAX.LT.1000000.0) then  ! Condition added by A. Kramida 08/20/04
            WRITE (21,16)TS,ERAS(1),WS,J1,AJ1,(IP1(I),C1(I),Y1(I),
     1        I=1,IM1)
          else
            WRITE (21,1601)TS,ERAS(1),WS,J1,AJ1,(IP1(I),C1(I),Y1(I),
     1        I=1,IM1)
          endif
C       ENDIF
      ENDIF
      IF (NSQ1.EQ.0.AND.IPRNSQ.NE.2) THEN
        IF (IPRNSQ.EQ.3) THEN
          IF (J.EQ.1) THEN
            WRITE (IW,91)
   91       FORMAT (1H )
            WRITE (IW,90) FJ(M), (YIDEN(I,NOUT,M), I=1,IM)
   90       FORMAT(1H ,A10,41X,8A10:/(50X,8A10:))
          ENDIF
          JOFS=NYOFS+(J-1)*IM
          if (EMAX.LT.10000.0) then  ! Condition added by A. Kramida 08/04/08
            WRITE (IW,7603) TS,ERAS(1), WS, TWS, GKK(J,M),
     1        GFAC(J,M), CONFIG(NJM,K3), (BUFY(I+JOFS), I=1,IM)
          else
            if (EMAX.LT.100000.0) then  ! Condition added by A. Kramida 08/04/08
              WRITE (IW,7602) TS,ERAS(1), WS, TWS, GKK(J,M),
     1          GFAC(J,M), CONFIG(NJM,K3), (BUFY(I+JOFS), I=1,IM)
            else
              if (EMAX.LT.1000000.0) then  ! Condition added by A. Kramida 08/20/04
                WRITE (IW,76) TS,ERAS(1), WS, TWS, GKK(J,M),
     1            GFAC(J,M), CONFIG(NJM,K3), (BUFY(I+JOFS), I=1,IM)
              else
                WRITE (IW,7601) TS,ERAS(1), WS, TWS, GKK(J,M),
     1            GFAC(J,M), CONFIG(NJM,K3), (BUFY(I+JOFS), I=1,IM)
              endif
            endif
          endif
   76     FORMAT (1H ,F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,8F10.3:/
     1      (49X,8F10.3:))
 7601     FORMAT (1H ,F10.1,A1,F10.1,F9.1,F7.3,F6.3,1X,A6,8F10.3:/
     1      (49X,8F10.3:))
 7602     FORMAT (1H ,F10.4,A1,F10.3,F9.1,F7.3,F6.3,1X,A6,8F10.3:/
     1      (49X,8F10.3:))
 7603     FORMAT (1H ,F10.5,A1,F10.3,F9.1,F7.3,F6.3,1X,A6,8F10.3:/
     1      (49X,8F10.3:))
        ELSE
          IF (J.EQ.1) WRITE (IW,20) FJ(M)
   20     FORMAT (A8)
          if (EMAX.LT.10000.0) then  ! Condition added by A. Kramida 08/04/08
            WRITE (IW,1703) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1        CONFIG(NJM,K3),(IP1(I),C1(I),Y1(I), I=1,LM1)
          else
            if (EMAX.LT.100000.0) then  ! Condition added by A. Kramida 08/04/08
              WRITE (IW,1702) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1          CONFIG(NJM,K3),(IP1(I),C1(I),Y1(I), I=1,LM1)
            else
              if (EMAX.LT.1000000.0) then  ! Condition added by A. Kramida 08/20/04
                WRITE (IW,17) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1            CONFIG(NJM,K3),(IP1(I),C1(I),Y1(I), I=1,LM1)
              else
                WRITE (IW,1701) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1            CONFIG(NJM,K3),(IP1(I),C1(I),Y1(I), I=1,LM1)
              endif
            endif
          endif
   17     FORMAT (F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,
     1      5(I4,1H ,A6,1X,A12:),
     2      (50X,5(I4,1H ,A6,1X,A12:)))
 1701     FORMAT (F10.1,A1,F10.1,F9.1,F7.3,F6.3,1X,A6,
     1      5(I4,1H ,A6,1X,A12:),
     2      (50X,5(I4,1H ,A6,1X,A12:)))
 1702     FORMAT (F10.4,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,
     1      5(I4,1H ,A6,1X,A12:),
     2      (50X,5(I4,1H ,A6,1X,A12:)))
 1703     FORMAT (F10.5,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,
     1      5(I4,1H ,A6,1X,A12:),
     2      (50X,5(I4,1H ,A6,1X,A12:)))
        ENDIF
      ENDIF
C
C     IF (NSQ1.EQ.1.AND.IPRNSQ.EQ.0.AND.NOUT.EQ.1) THEN
      IF (NSQ1.EQ.1.AND.IPRNSQ.EQ.0) THEN
        IF (J.EQ.1) WRITE (IW,20) FJ(M)
        if (EMAX.LT.1000000.0) then  ! Condition added by A. Kramida 08/20/04
          WRITE (IW,19) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1      CONFIG(NJM,K3),(IP2(I),C1(I),Y1(I), I=1,LM2)
        else
          WRITE (IW,1901) TS,ERAS(1),WS,TWS,GKK(J,M),GFAC(J,M),
     1      CONFIG(NJM,K3),(IP2(I),C1(I),Y1(I), I=1,LM2)
        endif
   19   FORMAT (F10.3,A1,F10.3,F9.3,F7.3,F6.3,1X,A6,
     1   5(I4,1H%,A6,1X,A12:),
     2   (/50X,5(I4,1H%,A6,1X,A12:)))
 1901   FORMAT (F10.1,A1,F10.1,F9.1,F7.3,F6.3,1X,A6,
     1   5(I4,1H%,A6,1X,A12:),
     2   (/50X,5(I4,1H%,A6,1X,A12:)))
      ENDIF
C *********************************************** A.KRAMIDA ^
  900 RETURN
      END
C
      SUBROUTINE RDJJ2E
C
C     A.Kramida 2026: READ THE TRUE jj TRAILER OF TAPE2E, IF THERE IS
C     ONE.  RCG APPENDS IT AFTER EVERYTHING RCE HAS EVER READ FROM UNIT
C     2, SO IT IS FOUND BY READING RECORD AFTER RECORD UNTIL THE MAGIC
C     ONE TURNS UP.  THREE INTEGERS ARE TAKEN FROM EACH RECORD ON THE
C     WAY: EVERY RECORD OF A TAPE2E IS LONGER THAN THAT, AND READING
C     FEWER ITEMS THAN A RECORD HOLDS SIMPLY SKIPS THE REST, SO THE SCAN
C     COSTS TWELVE BYTES PER RECORD RATHER THAN THE WHOLE FILE.  A
C     SHORT RECORD, IF ONE IS EVER MET, RAISES IOSTAT AND IS PASSED
C     OVER; ONLY END-OF-FILE STOPS THE SCAN.
C
C     UNIT 2 IS REWOUND ON THE WAY OUT WHATEVER HAPPENS, SO A CALLER
C     THAT KNOWS NOTHING OF THIS ROUTINE IS UNAFFECTED BY IT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      INCLUDE 'RCEJJ2E.F'
C       THE MAGIC AND VERSION NUMBERS OF THE TRAILER.  THEY MUST AGREE
C       WITH THE PARAMETER STATEMENT IN rcg11kd.f (SUBROUTINE ENERGY).
      PARAMETER (MAGJJ=20260802, MVERJJ=1)
C
      COMMON/CMLC/C_PTR,A_PTR
      INTEGER*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
      INTEGER*8 NJJS_PTR
      DIMENSION NJJS(1)
      POINTER (NJJS_PTR,NJJS)
C
      IJJON=0
      IJJCS=0
      NBKJ2=0
      REWIND 2
   10 READ (2,IOSTAT=IOS) I1,I2,I3
      IF (IOS.LT.0) GO TO 90
      IF (IOS.GT.0) GO TO 10
      IF (I1.NE.MAGJJ) GO TO 10
      IF (I2.NE.MVERJJ) THEN
        WRITE (6,11) I2,MVERJJ
   11   FORMAT (' TRUE jj TRAILER ON TAPE2E IS VERSION',I4,
     1    ', THIS RCE KNOWS VERSION',I4,' - IGNORED')
        GO TO 90
      ENDIF
      NB=I3
      IF (NB.LE.0) GO TO 90
      IF (NB.GT.NBJJX) THEN
        WRITE (6,12) NB,NBJJX
   12   FORMAT (' TRUE jj TRAILER HOLDS',I5,' J BLOCKS, MORE THAN THE',
     1    I5,' NBJJX ALLOWS - IGNORED')
        GO TO 90
      ENDIF
C
C       THE MATRICES GO INTO A BUFFER OF THEIR OWN, PACKED AS BUFA IS,
C       AND THE LABELS INTO ANOTHER.  BOTH GROW IF THEY HAVE TO.  THE
C       ROTATED MATRICES ARE ABOUT FIVE TIMES DENSER THAN THE J1J2 ONES
C       - A FIBRE ROTATION DESTROYS EXACTLY THE SPARSITY THAT MAKES J1J2
C       SPARSE - SO THIS BUFFER ENDS UP LARGER THAN BUFA.
      MAXBJU=MAXB2
      CALL ALLOC(BUFJ_PTR,MAXBJU,12)
      MAXYJ=IV
      CALL ALLOC(YJJ_PTR,MAXYJ,8)
      MXNJJS=4*IA
      CALL ALLOC(NJJS_PTR,MXNJJS,4)
      JOFFS=0
      JOFS1=0
      LYOFS=0
C
      DO 50 IBK=1,NB
        READ (2,IOSTAT=IOS) IOK,SJ8J,NLSB,NSBJ
        IF (IOS.NE.0) GO TO 80
        BACKSPACE 2
        IF (NLSB.LE.0.OR.NLSB.GT.IA.OR.NSBJ.LE.0) THEN
          WRITE (6,13) IBK,NLSB,NSBJ
   13     FORMAT (' TRUE jj TRAILER, BLOCK',I4,': MATRIX SIZE',I7,
     1      ' AND',I4,' SUBSHELLS MAKE NO SENSE - IGNORED')
          GO TO 90
        ENDIF
        NEED=4*NLSB*NSBJ
        IF (NEED.GT.MXNJJS) THEN
          CALL RALLOC(NJJS_PTR,MXNJJS,NEED,4)
          MXNJJS=NEED
        ENDIF
        IF (LYOFS+NLSB.GT.MAXYJ) THEN
          NEWSIZ=MAX0(2*MAXYJ,LYOFS+NLSB)
          CALL RALLOC(YJJ_PTR,MAXYJ,NEWSIZ,8)
          MAXYJ=NEWSIZ
        ENDIF
        READ (2,IOSTAT=IOS) IOK,SJ8J,NLSB,NSBJ,
     1    (((NJJS(K1+4*(I-1)+4*NSBJ*(L-1)),K1=1,4),I=1,NSBJ),
     2      YJJ(LYOFS+L), L=1,NLSB)
        IF (IOS.NE.0) GO TO 80
        NLSJ2(IBK)=NLSB
        NJJJ2(IBK)=NINT(2.D0*SJ8J)
        IOKJ2(IBK)=IOK
        LOFJ2(IBK)=LYOFS
        LYOFS=LYOFS+NLSB
        IMSQ=NLSB*NLSB
        CALL UNPK2(2,NLSB,A,IA,BUFJ_PTR,IMSQ,JOFS1,MAXBJU)
        CALL WR1(IBK,1,NLSB,A,IA,BUFJ_PTR,IBAJ,IMSQ,JOFFS,MAXBJU,
     *           NBJJX,1)
        JOFS1=JOFFS
   50 CONTINUE
      NBKJ2=NB
      IJJON=1
      OPEN(13,FILE='LEVELS1J',STATUS='UNKNOWN')
      OPEN(14,FILE='LEVELS2J',STATUS='UNKNOWN')
      OPEN(15,FILE='LEVELS3J',STATUS='UNKNOWN')
      REWIND 13
      WRITE (6,14) NB
   14 FORMAT (' TRUE jj DATA FOR',I5,' J BLOCKS READ FROM TAPE2E')
      GO TO 90
C
   80 WRITE (6,81) IBK,NB
   81 FORMAT (' TRUE jj TRAILER BREAKS OFF AT BLOCK',I5,' OF',I5,
     1  ' - IGNORED')
C
   90 REWIND 2
      RETURN
      END
C
      SUBROUTINE MAPJJ2
C
C     A.Kramida 2026: MATCH THE J MATRICES OF THE CSET NOW IN MEMORY
C     WITH THE TRAILER BLOCKS THAT BELONG TO THEM, AND TAKE THEIR TRUE
C     jj LABELS INTO YIDEN(.,3,.).
C
C     THE TRAILER COVERS EVERY CSET OF THE RCG RUN, IN THE ORDER CALCFC
C     PRODUCED THE BLOCKS, WHICH IS THE ORDER THEY APPEAR ON TAPE2E; SO
C     A CSET OCCUPIES MMAX CONSECUTIVE BLOCKS OF IT.  WHICH ONES IS
C     SETTLED BY LOOKING FOR A RUN WHOSE J AND MATRIX SIZE AGREE WITH
C     THE CSET'S THROUGHOUT, STARTING WHERE THE LAST CSET LEFT OFF.
C     GETTING BOTH RIGHT FOR EVERY BLOCK OF A RUN LEAVES NO ROOM FOR AN
C     ACCIDENT WORTH WORRYING ABOUT, AND A CSET THAT CANNOT BE PLACED
C     AT ALL SIMPLY GOES WITHOUT THE TRUE jj OUTPUT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      INCLUDE 'RCEJJ2E.F'
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
      DIMENSION NJ2(IX)
      SAVE IOFSAV
      DATA IOFSAV/0/
C
      IJJCS=0
      IF (IJJON.EQ.0) RETURN
      IF (MMAX.GT.NBKJ2) GO TO 80
      DO M1=1,MMAX
        READ (FJ(M1),10) FJ1
   10   FORMAT (4X,F4.1)
        NJ2(M1)=NINT(2.D0*FJ1)
      ENDDO
      NTRY=NBKJ2-MMAX+1
      DO 30 ITRY=1,NTRY
        IOFF=MOD(IOFSAV+ITRY-1,NTRY)
        DO M1=1,MMAX
          IF (NLSJ2(IOFF+M1).NE.IMAX(M1)) GO TO 30
          IF (NJJJ2(IOFF+M1).NE.NJ2(M1)) GO TO 30
        ENDDO
        GO TO 40
   30 CONTINUE
      GO TO 80
C
   40 IJJCS=1
      IOFSAV=IOFF+MMAX
      NBAD=0
      DO M1=1,MMAX
        MJJ(M1)=IOFF+M1
        L0=LOFJ2(MJJ(M1))
        DO L=1,IMAX(M1)
          YIDEN(L,3,M1)=YJJ(L0+L)
        ENDDO
        IF (IOKJ2(MJJ(M1)).EQ.0) NBAD=NBAD+1
      ENDDO
      IF (NBAD.NE.0) THEN
        WRITE (IW,41) NBAD
        IF (IW6.LT.0) WRITE (6,41) NBAD
   41   FORMAT (' TRUE jj: THE ROTATION FAILED IN RCG FOR',I4,
     1    ' J BLOCKS OF THIS CSET.  THEIR LABELS ARE BLANK AND THEIR'/
     2    ' EIGENVECTORS ARE THOSE OF COUPLING 2.')
      ENDIF
      RETURN
C
   80 WRITE (IW,81) NOCSET
      IF (IW6.LT.0) WRITE (6,81) NOCSET
   81 FORMAT (' TRUE jj: NO RUN OF TRAILER BLOCKS MATCHES CSET',I4,
     1  ' - NO TRUE jj OUTPUT FOR IT')
      RETURN
      END
C
      SUBROUTINE SORT8TO7
C
C     A.Kramida 2026: THE BODY OF THIS ROUTINE IS NOW SRTLEV, SO THAT
C     THE TRUE jj LEVEL LIST CAN BE SORTED THE SAME WAY INTO ITS OWN
C     FILES.  'LEVELS1' STILL CARRIES TWO SECTIONS PER CSET AND GIVES
C     THE SAME 'LEVELS2' AND 'LEVELS3' AS BEFORE; 'LEVELS1J' CARRIES ONE
C     PER CSET AND GIVES 'LEVELS2J' AND 'LEVELS3J'.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'RCE20KD_PAR.F'
      INCLUDE 'RCEJJ2E.F'
C
      CALL SRTLEV(8,7,12,2)
      IF (IJJON.NE.0) CALL SRTLEV(13,14,15,1)
      STOP 'NORMAL EXIT'
      END
C
      SUBROUTINE SRTLEV(IU1,IU2,IU3,NSEC)
C
C     SORT THE LEVEL LIST ON UNIT IU1 BY ENERGY INTO IU2, AND GROUP IT
C     BY CONFIGURATION INTO IU3.  NSEC IS THE NUMBER OF SECTIONS THE
C     FILE CARRIES PER CONFIGURATION SET (ONE PER COUPLING SCHEME).
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F1,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
      DIMENSION F(IV),JNK(4),JXK(4)
      CHARACTER A(IV)*136,FCH(IV)*6
C
      JNK(3)=0
      JXK(3)=0
      K=1
      REWIND IU1
      REWIND IU2
C     WRITE (IW,9)
C   9 FORMAT (1H1)
      I=0
      JN=1
      JX=0
      IEND=0
  100 I=I+1
      READ (IU1,10,END=500) A(I)
   10 FORMAT (A136)
      IF (A(I)(52:53).EQ.'J=') THEN
        JX=I
        GO TO 100
      ENDIF
      IF (JX.LE.0) THEN
        JN=I+1
        GO TO 100
      ENDIF
      GO TO 120
  500 IEND=1
  120 IF (JX.NE.JN) THEN
        DO J=JN,JX
          READ (A(J),13) F(J)
        ENDDO
   13   FORMAT (14X,F10.3)
        N=JX-JN+1
        CALL SORT(N,F(JN),A(JN))
      ENDIF
      DO J=1,JX
        WRITE (IU2,10) A(J)
      ENDDO
      JNM1=JN-1
C     DO 160 J=1,JNM1
C 160 WRITE (IW,16) A(J)
C  16 FORMAT (1X,A136)
C     DO 170 J=JN,JX
C 170 WRITE (IW,10) A(J)
      JNK(K)=JN
      JXK(K)=JX
      IF (IEND.LE.0) THEN
        K=K+1
        A(1)=A(I)
        I=1
        JN=2
        JX=0
        GO TO 100
      ENDIF
C
C
  600 REWIND IU2
      REWIND IU3
      DO K=1,4
        IF (JXK(K).EQ.0) GO TO 900
        K33=1+(K-1)/NSEC
        NMAX=N3(K33)
        JN=JNK(K)
        JX=JXK(K)
        JNM1=JN-1
        DO J=1,JX
          READ (IU2,10) A(J)
        ENDDO
        DO J=1,JNM1
          WRITE (IU3,10) A(J)
        ENDDO
C       DO 625 J=JN,JX
C 625   READ (A(J),62) FCH(J)
C  62   FORMAT (43X,A6)
        DO N=1,NMAX
          ERAS1=CONFIG(N,K33)
          DO J=JN,JX
            IF (A(J)(44:51).EQ.ERAS1) WRITE (IU3,10) A(J)
          ENDDO
          WRITE (IU3,10)
        ENDDO
      ENDDO
  900 RETURN
      END
C
      SUBROUTINE SORT(N,F,A)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      DIMENSION F(IV)
      CHARACTER*136 A(IV),A1
C
      J=1
   60 IF (F(J).GT.F(J+1)) THEN
        M=J
        ERAS=F(M+1)
        A1=A(M+1)
   70   F(M+1)=F(M)
        A(M+1)=A(M)
        M=M-1
        IF (M.EQ.0) GO TO 90
        IF (F(M).GT.ERAS) GO TO 70
   90   F(M+1)=ERAS
        A(M+1)=A1
      ENDIF
   80 J=J+1
      IF (J.LT.N) GO TO 60
      RETURN
      END
C
C ***************************************** A.Kramida 2026 v
      SUBROUTINE ESORTW(IM,KWM,WM)
C
C          FILL KWM(1..IM) WITH THE BLOCK'S LEVEL INDICES ORDERED BY
C          INCREASING CALCULATED ENERGY WM.  EQUAL ENERGIES COME OUT IN
C          DECREASING INDEX ORDER, WHICH IS WHAT THE ORIGINAL SEQUENTIAL
C          EXTRACTION OF MINIMA PRODUCED: ITS TEST 'W(N1,M).LE.GRE' KEPT
C          THE LAST OF SEVERAL EQUAL ENERGIES.
C
C          THAT EXTRACTION RESCANNED THE WHOLE BLOCK FOR EVERY LEVEL AND,
C          FOR EVERY CANDIDATE, RESCANNED EVERYTHING PICKED SO FAR -
C          O(IM**3) ALTOGETHER, WHICH DOMINATED THE POST-ITERATION OUTPUT
C          PHASE FOR LARGE BLOCKS (ABOUT 3 MIN PER PARITY IN GD III).
C          THE HEAPSORT BELOW COSTS O(IM*LOG(IM)) AND YIELDS THE SAME
C          PERMUTATION, TIES INCLUDED: IDXSRT ORDERS BY DECREASING KEY
C          AND, WITH ITIE=0, BREAKS TIES BY INCREASING INDEX, SO
C          REVERSING ITS RESULT GIVES EXACTLY THE ORDER REQUIRED.
C
C          THE WORKING INDEX ARRAY IS INTEGER*4 BECAUSE THAT IS WHAT
C          IDXSRT EXPECTS; KW ITSELF IS INTEGER*2, WHICH IS WIDE ENOUGH
C          FOR IA=20000.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER*2 KWM(*)
      DIMENSION WM(*)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      INTEGER*4 IDX(IA)
      SAVE IDX
C
      DO K=1,IM
        IDX(K)=K
      ENDDO
      CALL IDXSRT(IM,IDX,WM,0)
      DO K=1,IM
        KWM(K)=IDX(IM+1-K)
      ENDDO
      RETURN
      END
C
      SUBROUTINE IDXSRT(N,IDX,KEY,ITIE)
C
C       ORDER IDX(1..N) SO THAT KEY(IDX(1)) >= KEY(IDX(2)) >= ...
C       EQUAL KEYS ARE ORDERED BY DECREASING INDEX WHEN ITIE=1 AND BY
C       INCREASING INDEX WHEN ITIE=0.  THAT REPRODUCES THE ORDER IN
C       WHICH THE ORIGINAL SEQUENTIAL SEARCHES PICKED THEIR MAXIMA:
C       'IF (T(N1,M).GE.GRE)' KEPT THE LAST OF SEVERAL EQUAL VALUES
C       (ITIE=1) WHILE 'IF (EVAL(N1,M).GT.GRE)' KEPT THE FIRST (ITIE=0).
C       HEAPSORT IN PLACE, O(N*LOG(N)), REPLACING REPEATED SEARCHES
C       THAT COST O(N**3) ALTOGETHER.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDX(*)
      REAL*8 KEY(*)
C
      IF (N.LE.1) RETURN
      DO L=N/2,1,-1
        CALL IDXSFT(IDX,KEY,ITIE,L,N)
      ENDDO
      DO L=N,2,-1
        IT=IDX(1)
        IDX(1)=IDX(L)
        IDX(L)=IT
        CALL IDXSFT(IDX,KEY,ITIE,1,L-1)
      ENDDO
      RETURN
      END
C
      SUBROUTINE IDXSFT(IDX,KEY,ITIE,LROOT,LEND)
C
C       SIFT IDX(LROOT) DOWN THROUGH IDX(1..LEND).  USED BY IDXSRT.
C       THE HEAP IS ORDERED SO THAT ITS ROOT IS THE ELEMENT THAT HAS TO
C       BE PLACED LAST, WHICH LEAVES IDX IN DECREASING ORDER OF KEY.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION IDX(*)
      REAL*8 KEY(*)
C
      LP=LROOT
   10 LC=2*LP
      IF (LC.GT.LEND) RETURN
      IF (LC.LT.LEND) THEN
        IF (ILATER(IDX(LC+1),IDX(LC),KEY,ITIE).NE.0) LC=LC+1
      ENDIF
      IF (ILATER(IDX(LC),IDX(LP),KEY,ITIE).EQ.0) RETURN
      IT=IDX(LC)
      IDX(LC)=IDX(LP)
      IDX(LP)=IT
      LP=LC
      GO TO 10
      END
C
      FUNCTION ILATER(IA1,IB1,KEY,ITIE)
C
C       1 IF IA1 HAS TO BE PLACED AFTER IB1, 0 OTHERWISE.  USED BY
C       IDXSFT.  SEE IDXSRT FOR THE MEANING OF ITIE.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 KEY(*)
C
      ILATER=0
      IF (KEY(IA1).LT.KEY(IB1)) THEN
        ILATER=1
      ELSE IF (KEY(IA1).EQ.KEY(IB1)) THEN
        IF (ITIE.NE.0) THEN
          IF (IA1.LT.IB1) ILATER=1     ! EQUAL KEYS: LARGER INDEX FIRST
        ELSE
          IF (IA1.GT.IB1) ILATER=1     ! EQUAL KEYS: SMALLER INDEX FIRST
        ENDIF
      ENDIF
      RETURN
      END
C
      SUBROUTINE EVDIST(IVN,JVS,M,N1F,IDD)
C
C       SQUARED DISTANCE BETWEEN NEW EIGENVECTOR IVN AND SAMPLE
C       EIGENVECTOR JVS IN THE SPACE OF THEIR N1F LARGEST COMPONENTS,
C       MINIMIZED OVER THE TWO POSSIBLE OVERALL SIGNS.  THIS IS THE
C       SAME MEASURE THE FULL MATCHER IN ORDER USES; IT IS EVALUATED
C       DIRECTLY HERE BECAUSE THE RELAXED MATCHING (INVECT=2) NEEDS IT
C       FOR ONLY A FEW SLOTS AND CANNOT AMORTIZE AN INVERTED INDEX.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      INCLUDE 'RCE20KD_PAR.F'
      INTEGER*2 NW,NW0,IEVV,IEVV0
      COMMON/VECT/ IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      INTEGER*8 IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      DIMENSION IEVV(MAXC,IA,IX),IEVV0(MAXC,IA,IX),NW(MAXC,IA,IX),
     1 NW0(MAXC,IA,IX)
      POINTER(IEVV_PTR,IEVV),(IEVV0_PTR,IEVV0),(NW_PTR,NW),(NW0_PTR,NW0)
C
      IDD=0
      IDDN=0
      NFOUND=0
      DO K=1,N1F
        N1=NW(K,IVN,M)
        IA0=IEVV(K,IVN,M)
        N2=0
        DO K1=1,N1F
          NK1=NW0(K1,JVS,M)
          IF (NK1.EQ.0) GO TO 46
          IF (NK1.EQ.N1) THEN
            N2=K1
            NFOUND=NFOUND+1
            GO TO 46
          ENDIF
        ENDDO
   46   IA1=0
        IF (N2.NE.0) IA1=IEVV0(N2,JVS,M)
        IA10=IA1-IA0
        IDD=IDD+IA10*IA10
        IA10=IA1+IA0                   ! THE INVERTED-SIGN VARIANT
        IDDN=IDDN+IA10*IA10
      ENDDO
      IF (NFOUND.LT.N1F) THEN
C           ADD THE COMPONENTS OF NW0(JVS) MISSING IN NW(IVN)
        DO K=1,N1F
          N1=NW0(K,JVS,M)
          IA0=IEVV0(K,JVS,M)
          N2=0
          DO K1=1,N1F
            NK1=NW(K1,IVN,M)
            IF (NK1.EQ.0) GO TO 49
            IF (NK1.EQ.N1) GO TO 50
          ENDDO
   49     IDD=IDD+IA0*IA0
          IDDN=IDDN+IA0*IA0
   50     CONTINUE
        ENDDO
      ENDIF
      IF (IDDN.LT.IDD) IDD=IDDN
      RETURN
      END
C ***************************************** A.Kramida 2026 ^
C
      SUBROUTINE CLOCK(T)
C
C          TIMING ROUTINE (RESULTS FOR INFORMATION ONLY)
C               FOR OTHER MACHINES, SET T=0.0 OR ADD
C               APPROPRIATE CODE.
C
      REAL*8 T
c      REAL ETIME,T1,T2
c      DIMENSION T2(2)
C
C------------------------------------
C          USE THIS CODE FOR A CRAY
C     T=SECOND(T)
C------------------------------------
C          USE THIS CODE FOR A VAX OR MACINTOSH
C     T1=SECNDS(0.0)
C     T=DBLE(T1)
C------------------------------------
C          USE THIS CODE FOR A SUN
C     T=ETIME(T2)
C-------------------------------------
C          USE THIS CODE FOR IBM RISC
C     IT=MCLOCK()
C     T=IT
C     T=T/100.0
C-------------------------------------
C          USE THIS CODE FOR IBM-PC LAHEY FORTRAN COMPILER
c      REAL*8 TINIT,T3
c      COMMON/TIMING/TINIT,ITINIT
c      INTEGER HH,MM,SS,HS
c      CHARACTER TIMESTR*11
c      CALL TIME(TIMESTR)
c      READ (TIMESTR,10) HH,MM,SS,HS
c   10 FORMAT (3(I2,1X),I2)
c      T3=DFLOAT(3600*HH+60*MM+SS)+HS*1.D-2
c      IF (ITINIT.EQ.0) THEN
c        TINIT=T3
c        ITINIT=1
c      ENDIF
c      T=T3-TINIT
c------------------------------------
C          USE THIS CODE FOR the G95 compiler on Windows systems
C      REAL*8 T
c      real*4 tarray(2)
c      T=dble(etime(tarray))
c     T=dble(tarray(1))
C-------------------------------------
C          USE THIS CODE FOR LINUX OR macOS
C     NOTHING HAS TO BE UNCOMMENTED: THE ACTIVE CODE AT THE END OF
C     THIS ROUTINE IS SYSTEM_CLOCK, WHICH IS FORTRAN STANDARD AND
C     NEEDS NO LIBRARY, AND IS THE RIGHT CALL ON LINUX, macOS AND
C     WINDOWS ALIKE.  THIS BLOCK RECORDS THE ALTERNATIVES AND THE ONE
C     TRAP, BECAUSE THE BLOCKS ABOVE ARE ORGANIZED BY MACHINE AND A
C     READER LOOKING FOR "LINUX" WOULD OTHERWISE FIND NOTHING HERE.
C       THE TRAP: DECLARE SYSTEM_CLOCK'S ARGUMENTS INTEGER*8, NOT THE
C       DEFAULT INTEGER.  gfortran GIVES THE INTRINSIC MILLISECOND
C       RESOLUTION FOR INTEGER*4 ARGUMENTS AND NANOSECOND RESOLUTION
C       (clock_gettime) FOR INTEGER*8 ONES, ON BOTH LINUX AND macOS.
C       THAT IS WHY THE ACTIVE CODE BELOW DECLARES THEM INTEGER*8.
C     WALL CLOCK THE OTHER WAY, WHEN -fopenmp IS IN USE ANYWAY:
C     REAL*8 OMP_GET_WTIME
C     T=OMP_GET_WTIME()
C     CPU TIME, IF CPU TIME IS WHAT IS WANTED:
C     CALL CPU_TIME(T)
C       CPU_TIME IS STANDARD AND WORKS EVERYWHERE TOO, BUT IT SUMS
C       OVER THE THREADS -- WHICH IS EXACTLY THE BUG DESCRIBED IN THE
C       2026 NOTE BELOW.  DO NOT USE IT IN A THREADED BUILD.
C     THE GNU EXTENSIONS IN THE BLOCKS ABOVE ALSO COMPILE UNDER
C     gfortran ON LINUX AND macOS -- ETIME AND MCLOCK (CPU TIME, SO
C     THEY CARRY THAT SAME BUG) AND SECNDS (TIME OF DAY, SO IT WRAPS
C     AT MIDNIGHT).  NONE OF THEM IS PORTABLE TO ANOTHER COMPILER.
C     PREFER THE STANDARD INTRINSICS.
C-------------------------------------
C     A.Kramida 2026: ETIME RETURNS CPU TIME SUMMED OVER ALL THREADS,
C       WHICH OVERSTATES ELAPSED TIME BY THE NUMBER OF ACTIVE THREADS
C       ONCE THE OPENMP REGIONS ARE IN USE.  USE THE WALL CLOCK INSTEAD.
C       ICR=0 MEANS NO CLOCK IS AVAILABLE; FALL BACK TO CPU TIME THEN.
      INTEGER*8 ICNT,ICR,ICMAX
      INTEGER*8 ICNT0
      LOGICAL LFIRST
      SAVE ICNT0,LFIRST
      DATA LFIRST/.TRUE./
      CALL SYSTEM_CLOCK(ICNT,ICR,ICMAX)
      IF (ICR.EQ.0) THEN
        T=0.D0
      ELSE
        IF (LFIRST) THEN
          ICNT0=ICNT
          LFIRST=.FALSE.
        ENDIF
C       GUARD AGAINST COUNTER WRAP-AROUND
        IF (ICNT.LT.ICNT0) ICNT=ICNT+ICMAX
        T=DBLE(ICNT-ICNT0)/DBLE(ICR)
      ENDIF
      T=T/60.D0
      RETURN
      END
C
      SUBROUTINE CALCV
C
C          CALCULATES V MATRIX
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC_PTR,BUFA_PTR,BUFY_PTR,
c     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
     1  IBC_PTR,IBA_PTR,ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
      INTEGER*8 BUFY_PTR,BUFA_PTR,BUFC_PTR,IBA_PTR,IBC_PTR
      POINTER (BUFY_PTR,BUFY(1)),(BUFA_PTR,BUFA(1)),(BUFC_PTR,BUFC(1))
      POINTER (IBA_PTR,IBA),(IBC_PTR,IBC)
      DIMENSION IBC(IX,IB),IBA(IX,1)
c
      INTEGER*2 LF1
      DIMENSION LF1(IB)
C ***************************************** A.Kramida 2026 v
C       COMPRESSED-ROW FORM OF THE NON-ZERO ELEMENTS OF C, HELD FOR THE
C       WHOLE K LOOP BELOW.  IRSTA/IROWL ARE INDEXED BY THE ORDINAL OF A
C       NON-EMPTY ROW, NOT BY THE ROW ITSELF, SO ROWS THAT ARE ENTIRELY
C       ZERO COST NOTHING AT ALL.  MAXNZ4 IS A CEILING, NOT A
C       PREDICTION: IF A C EVER HAS MORE NON-ZEROS THAN THIS, THE
C       ORIGINAL DENSE LOOP IS USED FOR IT.  SAVEd - THESE ARE FAR TOO
C       BIG FOR THE STACK, AND -fno-automatic MAKES THEM STATIC ANYWAY.
      PARAMETER (MAXNZ4=2000000)
      INTEGER*4 IRSTA(IA+1),IROWL(IA),IRCOL(MAXNZ4)
      REAL*8 RVAL(MAXNZ4)
      SAVE IRSTA,IROWL,IRCOL,RVAL
C ***************************************** A.Kramida 2026 ^
C
C
C          CALCULATE V MATRIX
C
      IM=0
      FF = 1.D0/F
C     REWIND I216                    ! R.COWAN
C     REWIND 4                      ! R.COWAN
      KMX=0
      LMAXP1=LMAX+1
      MIK=0
      NYOFS=0
      DO M=1,MMAX
        IM=IMAX(M)
        IMSQ=IM*IM
        IMSQ1=(IMSQ+IM)/2
C       CALL RD(I216,IM,Y,IA,EVEC,IMSQ)            ! R.COWAN
C       DO 100 K=1,IM                             ! A.KRAMIDA
C       DO 100 I=1,IM
C         MIK=MIK+1
C 100   Y(I,K)=BUFY(MIK)
        JL=0
        DO L=1,LMAX
          LFL=LF(L)
          IF (IABS(LFL).LT.100.OR.NOCYCE.EQ.0) THEN
            JL=JL+1
            IF (JL.GT.JB) STOP 'JB DIMENSION TOO SMALL'
C           CALL RD(4,IM,C,IA,EVEC,IMSQ1)           ! R.COWAN
            CALL RD1(M,L,IM,C,IA,BUFC,IBC,IMSQ1,    ! A.KRAMIDA
     *               ICOFFS,MAXB1,MMAX,LMAXP1)
            KM=KMX
C ***************************************** A.Kramida 2026 v
C           V(KM,JL) IS THE QUADRATIC FORM Y(.,K)' * C * Y(.,K) FOR EACH
C           FITTED LEVEL K OF THIS J BLOCK.  IT WAS SWEPT DENSELY - IM**2
C           PER K, HENCE (NUMBER OF FITTED LEVELS)*IM**2 PER (BLOCK,
C           PARAMETER), REPEATED EVERY LEAST-SQUARES ITERATION.  THAT ONE
C           LOOP WAS 68% OF A Pr III FIT.
C             C IS THE COEFFICIENT MATRIX OF A SINGLE PARAMETER AND IS
C           VERY SPARSE - A FEW TENS OF NON-ZEROS IN IM**2 IS TYPICAL,
C           WHICH IS WHY WR1 STORES IT AS A LIST OF NON-ZERO ELEMENTS IN
C           THE FIRST PLACE.  ITS NON-ZEROS ARE THEREFORE COLLECTED ONCE
C           PER (BLOCK,PARAMETER), OUTSIDE THE K LOOP, AND ONLY THOSE ARE
C           SWEPT.  THE COST BECOMES IM**2 FOR THE SCAN - THE SAME ORDER
C           AS THE RD1 CALL JUST ABOVE ALREADY COSTS - PLUS O(NON-ZEROS)
C           PER K INSTEAD OF IM**2 PER K.
C             THIS IS BIT-EXACT.  WITHIN A ROW THE SURVIVING PRODUCTS ARE
C           VISITED IN INCREASING J, EXACTLY AS BEFORE, SO NOTHING IS
C           REASSOCIATED; WHAT IS DROPPED IS THE ADDITION OF C(I,J)=0 TO
C           SUMJ, AND SUMJ STARTS AT +0.D0 AND SO CAN NEVER BE -0.D0
C           (X+Y=-0.D0 NEEDS BOTH X AND Y TO BE -0.D0), WHICH MAKES
C           SUMJ+(+/-0.D0) EXACTLY SUMJ.  THE SAME ARGUMENT APPLIES TO
C           SUMIJ, AND IS WHY A ROW OF C THAT IS ENTIRELY ZERO CAN BE
C           SKIPPED ALTOGETHER: IT CONTRIBUTES SUMJ*Y(I,K) WITH SUMJ
C           EXACTLY +0.D0, I.E. +/-0.D0, TO A SUMIJ THAT IS NEVER -0.D0.
C           THE DIAGONAL'S EXTRA FACTOR FF IS STILL APPLIED AFTER THE
C           MULTIPLICATION BY Y, NOT FOLDED INTO C, BECAUSE (C*Y)*FF AND
C           (C*FF)*Y ARE NOT THE SAME DOUBLE.
            NNZC=0
            NROWC=0
            DO I=1,IM
              IST=NNZC+1
              DO J=1,IM
                TEMP=C(I,J)
                IF (TEMP.NE.0.D0) THEN
                  IF (NNZC.GE.MAXNZ4) THEN
                    NNZC=-1
                    GO TO 101
                  ENDIF
                  NNZC=NNZC+1
                  IRCOL(NNZC)=J
                  RVAL(NNZC)=TEMP
                ENDIF
              ENDDO
              IF (NNZC.GE.IST) THEN
                NROWC=NROWC+1
                IROWL(NROWC)=I
                IRSTA(NROWC)=IST
                IRSTA(NROWC+1)=NNZC+1
              ENDIF
            ENDDO
  101       CONTINUE
C ***************************************** A.Kramida 2026 ^
            DO K=1,IM
              IF (NF(K,M).GE.0) THEN
                KM=KM+1
                V(KM,JL)=0.D0
                SUMIJ=0.D0
                KOFS=NYOFS+(K-1)*IM
                IF (NNZC.GE.0) THEN
                  DO IR=1,NROWC
                    I=IROWL(IR)
                    SUMJ=0.D0
                    DO IP=IRSTA(IR),IRSTA(IR+1)-1
                      JC=IRCOL(IP)
                      TEMP=RVAL(IP)*BUFY(KOFS+JC)  !Y(J,K)
                      IF (JC.EQ.I) TEMP=TEMP*FF
                      SUMJ=SUMJ+TEMP
                    ENDDO
                    SUMIJ=SUMIJ+SUMJ*BUFY(KOFS+I)  !Y(I,K)
                  ENDDO
                ELSE
                  DO I=1,IM
                    SUMJ=0.D0
                    DO J=1,IM
                      TEMP=C(I,J)
                      IF (TEMP.NE.0.D0) THEN
                        TEMP=TEMP*BUFY(KOFS+J)     !Y(J,K)
                        IF (I.EQ.J) TEMP=TEMP*FF
                      ENDIF
                      SUMJ=SUMJ+TEMP
                    ENDDO
                    SUMIJ=SUMIJ+SUMJ*BUFY(KOFS+I)  !Y(I,K)
                  ENDDO
                ENDIF
                V(KM,JL)=SUMIJ*F
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        KMX=KM
        NYOFS=NYOFS+IMSQ
      ENDDO
C
C     IF (IPRNV.EQ.0.AND.N5.NE.0) THEN
      IF (IPRNV.EQ.0) THEN
        WRITE (IW,20)
   20   FORMAT (/34H      T        T-W         V(KM,L))
        KM=0
        DO M=1,MMAX
          IM=IMAX(M)
          WRITE (IW,21) FJ(M)
   21     FORMAT (A10)
          DO K=1,IM
            IF (NF(K,M).GE.0) THEN
              KM=KM+1
              LKM=0
              JL=0
              DO L=1,LMAX
                LFL=LF(L)
                IF (IABS(LFL).LT.100) THEN
                  JL=JL+1
                  LF1(L)=JL
                  IF (DABS(V(KM,JL)).GE.1.D-4) THEN
                    LKM=LKM+1
                    LDEG(LKM)=L
                  ENDIF
                ENDIF
              ENDDO
              WRITE (IW,22) T(K,M),TW(K,M), ( V( KM,LF1( LDEG(L) ) ),
     1          NX( LDEG(L) ), XIDEN( LDEG(L) ) , L=1,LKM)
   22         FORMAT (2F10.3,3(2X,F9.4,1H ,A2,A8)/
     1          1(20X,3(2X,F9.4,1H ,A2,A8)) )
            ENDIF
C
          ENDDO
        ENDDO
      ENDIF
C
C          IF LSQM NON-ZERO, CALCULATE V*
C
  230 IF (LSQM.NE.0) THEN
C
  260   DO KM=1,KMX
          DO J=1,JLMAX
            VC(J)=0.D0
          ENDDO
          JL=0
          DO L=1,LMAX
            IF (LF(L).GT.-100.AND.LF(L).LT.0) THEN
              JL=JL+1
              J=-LF(L)
              VC(J)=VC(J)+V(KM,JL)*X(L)
            ENDIF
          ENDDO
          DO J=1,JLMAX
            V(KM,J)=VC(J)
          ENDDO
        ENDDO
        IF (IPRNV.EQ.0) THEN
          WRITE (IW,27)
   27     FORMAT (/35H      T        T-W         V*(KM,L))
          KM=0
          DO M=1,MMAX
            IM=IMAX(M)
            WRITE (IW,21) FJ(M)
            DO K=1,IM
              IF (NF(K,M).GE.0) THEN
                KM=KM+1
                LKM=0
                DO L=1,LMAX
                  LFL=LF(L)
                  J=IABS(LFL)
                  IF (J.LT.100.AND.DABS(V(KM,J)).GE.1.D-4) THEN
                    LKM=LKM+1
                    LDEG(LKM)=L
                  ENDIF
                ENDDO
                LFL=LF(LDEG(L))
                WRITE (IW,22) T(K,M),TW(K,M), (V(KM,IABS(LFL)),
     1             NX(LDEG(L)), XIDEN(LDEG(L)), L=1,LKM)
              ENDIF
C
            ENDDO
          ENDDO
        ENDIF
      ENDIF
C
C
C          CALCULATE B AND D MATRICES
C
      DO L=1,JLMAX
        DL1=0.D0
        DO J=L,JLMAX
          BLJ=0.D0
          KM=0
          DO M=1,MMAX
            IM=IMAX(M)
            DO K=1,IM
              IF (NF(K,M).GE.0) THEN
                KM=KM+1
                IF (J.EQ.L) THEN
                  DL1=DL1+V(KM,L)*TW(K,M)
                ENDIF
                BLJ=BLJ+V(KM,L)*V(KM,J)
              ENDIF
            ENDDO
          ENDDO
          B(L,J)=BLJ
          B(J,L)=BLJ
        ENDDO
        D(L,1)=DL1
        IF (B(L,L).EQ.0.D0) B(L,L)=1.D0
      ENDDO
c      if (JLMAX.EQ.43) write (6,300) B(32,32)
c 300  format (' B(32,32)=',F20.3)
C
C     IF (IPRNV.EQ.0.AND.N5.NE.0) THEN
C       WRITE (IW,32)
C  32   FORMAT (///30H           D            B(L,J))
C       WRITE (IW,21)
C       DO 325 L=1,JLMAX
C         IF (LSQM.GE.0) THEN
C           WRITE (IW,33) D(L,1), (B(L,J), J=1,JLMAX)
C         ELSE
C           WRITE (IW,34) D(L,1), (B(L,J), J=1,JLMAX)
C  34       FORMAT (1H ,F15.4,1PE14.3,8E11.3/(E32.3,8E11.3))
C       ENDIF
C 325   CONTINUE
C  33   FORMAT (1H , F15.4, F12.3, 10F9.3/(F30.3, 10F9.3))
C     ENDIF
      RETURN
      END
C
      SUBROUTINE CALCX
C
C                 LEAST-SQUARE CALCN OF NEW PARAMETERS (MOD 19)
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C
C          CALCULATE DELX AND NEW X
C
  340 DO L=1,JLMAX
        DO J=L,JLMAX
          IF (L.EQ.J) THEN
            D(L,L+1)=1.D0
          ELSE
            D(L,J+1)=0.D0
            D(J,L+1)=0.D0
          ENDIF
        ENDDO
      ENDDO
      JLMAX1=JLMAX+1
c  350 CALL LSS(JLMAX,JLMAX1,JB,B,D,ERAS,DETB)
  350 CALL LSS(JLMAX,JLMAX1,JB,B,D,DETB)
      SUMDX0=SUMDX
      SUMDX=0.D0
      DO L=1,LMAX
        LFL=LF(L)
        IF (IABS(LFL).GE.100) THEN
          DELXC(L)=0.D0
          SDX(L)=0.D0
        ELSE
          IF (LSQM.LT.0) THEN
            J=-LF(L)
            IF (DELX(J).LE.-1.D0.AND.XMIN(L).GT.0.D0) DELX(J)=-0.99D0
            DELXC(L)=X(L)*DELX(J)
            DLL=D(J,J+1)
            IF (DLL.LT.0) DLL=0
            SDX(L)=DABS(X(L))*DSQRT(DLL)
          ELSE
            J=LF(L)
            DELXC(L)=DELX(J)
            DLL=D(J,J+1)
            IF (DLL.LT.0) DLL=0
            SDX(L)=DSQRT(DLL)
          ENDIF
        ENDIF
        SUMDX=SUMDX+DELXC(L)
      ENDDO
      GEFF=1.D0
      IF (NOF.GT.3) THEN
c        write (6,3551) NOF,SUMDX0,SUMDX
c 3551   format (' NOF=',I5,'    SUMDX0=',F20.10,'    SUMDX=',F20.10)
        ERS1=0.5D0*GCOE*SUMDX/SUMDX0
        GEFF=1.D0+SIGN(DMIN1(0.5D0,DABS(ERS1)),ERS1)
      ENDIF
  357 DXMAX1=DXMAX
      IF (F.LT.0.9999D0) THEN
        DO L=1,LMAX
          GEFF=DMIN1(GEFF, DXMAX1/(ABS(DELXC(L))+0.001E0))
        ENDDO
      ELSE
        IF (F.GE.0.7D0) THEN
          DXMAX1=2.D0*DXMAX
          DO L=1,LMAX
            GEFF=DMIN1(GEFF, DXMAX1/(ABS(DELXC(L))+0.001E0))
          ENDDO
        ENDIF
      ENDIF
  365 DO L=1,LMAX
        X0(L)=X(L)
        DELX(L)=DELXC(L)*GEFF
        X(L)=X0(L)+DELX(L)
      ENDDO
  370 DO L=1,LMAX
        IF (X(L).GT.XMAX(L)) THEN
          X(L)=XMAX(L)
          DELX(L)=X(L)-X0(L)
        ELSE
          IF (X(L).LE.XMIN(L)) THEN
            X(L)=XMIN(L)
            DELX(L)=X(L)-X0(L)
          ENDIF
        ENDIF
      ENDDO
C     DELX15=0.D0
C     DO 380 L=1,LMAX
C       IF (LF(L).EQ.-15) DELX15=DELX(L)
C 380 CONTINUE
C     DO 390 L=1,LMAX
C       IF (LF(L).EQ.115) THEN
C         DELX(L)=DELX15
C         X(L)=X0(L)+DELX15
C       ENDIF
C 390 CONTINUE
C
      RETURN
      END
C
      SUBROUTINE EIGEN(kcpl,IBY_SIZE)               ! A.Kramida 01/31/2011: before that there were no parameters
C
C                CALCN OF NEW EIGENVALUES AND VECTORS
C
C
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
c
      DIMENSION KCPL(2)                    ! A.Kramida 01/31/2011: before that there were no parameters
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC_PTR,BUFA_PTR,BUFY_PTR,
c     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
     1  IBC_PTR,IBA_PTR,ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
      INTEGER*8 BUFY_PTR,BUFA_PTR,BUFC_PTR,IBA_PTR,IBC_PTR
      POINTER (BUFY_PTR,BUFY),(BUFA_PTR,BUFA(1)),(BUFC_PTR,BUFC(1))
      POINTER (IBA_PTR,IBA),(IBC_PTR,IBC)
      DIMENSION IBC(IX,IB),IBA(IX,1),BUFY(IBY_SIZE)
C
      R=0.D0
      TNO=0.D0
C ======================================== R.COWAN
C     REWIND I216
C     REWIND 4
C     REWIND 15
C ======================================== R.COWAN
      LMAXP1=LMAX+1
C     MIK=0
      NYOFS=0                  ! NYOFS IS OFFSET IN THE Y-BUFFER
      DO 520 M=1,MMAX
        IM=IMAX(M)
        IMSQ=IM*IM
        IMSQ1=(IMSQ+IM)/2
C
C            CALCULATE A MATRIX
C
        IF (IPRNA.EQ.0) THEN
C           A.Kramida 2026: IB+2 COLUMNS, NOT IA - THAT IS ALL A1 HOLDS
C           AND ALL THAT IS ALLOCATED FOR IT.  ALSO SPARES THE ZEROING
C           OF IA*IM WORDS ONCE PER BLOCK PER CYCLE.
          DO L=1,IB+2
            DO I=1,IM
              A1(I,L)=0.D0
            ENDDO
          ENDDO
        ENDIF
        DO I=1,IM
          DO J=1,I
            A(I,J)=0.D0
          ENDDO
        ENDDO
        DO 402 L=1,LMAX
C ======================================== R.COWAN
C         CALL RD(4,IM,C,IA,EVEC,IMSQ1)
C ======================================== R.COWAN
C *********************************************** A.KRAMIDA v
          CALL RD1(M,L,IM,C,IA,BUFC,IBC,IMSQ1,
     *           ICOFFS,MAXB1,MMAX,LMAXP1)
C *********************************************** A.KRAMIDA ^
C *********************************************** A.KRAMIDA  v
          IF (IPRNA.EQ.0) THEN
            DO I=1,IM
C                SET UP A1(I,J) TO HOLD NON-ZERO ANGULAR COEFFICIENTS
C                FOR SLATER INTEGRALS (IN THE DIAGONAL OF THE MATRIX
C                TO BE DIAGONALIZED)
              IF (DABS(C(I,I)).GT.1.D-5) THEN
                IF (A1(I,1).EQ.0.D0) A1(I,1)=DFLOAT(L)
                I1=A1(I,1)
                A1(I,2)=DFLOAT(L)
                I2=3+L-I1
C                 A.Kramida 2026: CLAMP TO THE A1 COLUMN COUNT, IB+2.
                IF (I2.GT.IB+2) I2=IB+2
                A1(I,I2)=C(I,I)
              ENDIF
            ENDDO
          ENDIF
C *********************************************** A.KRAMIDA  ^
          EVKM=X(L)
          DO I=1,IM
            DO J=1,I
              TEMP=C(I,J)
              IF (TEMP.NE.0.D0) A(I,J)=A(I,J)+TEMP*EVKM
            ENDDO
          ENDDO
  402   CONTINUE
        EVKM=0.D0
        DO I=1,IM
          DO J=1,I-1
            TEMP=A(I,J)
            IF (F.NE.1.D0) THEN
              TEMP=TEMP*F
              A(I,J)=TEMP
              A(J,I)=TEMP
            ENDIF
            TEMP=TEMP*TEMP
            EVKM=EVKM+TEMP+TEMP
          ENDDO
          TEMP=A(I,I)
          EVKM=EVKM+TEMP*TEMP
        ENDDO
        ANORM(M)=DSQRT(EVKM)
C
C            IF IPRNA IS ZERO, WRITE TAPE 9  ('OUTE')
C
C       IF (IPRNA.EQ.0.AND.N5.NE.0) THEN
C         IF (M.EQ.1) THEN
C           WRITE (IW,40) NORUN,NOCYC
C  40       FORMAT (1H1,'RUN NO=',A6,'     CYC NO=',I3,'        A(I,J)')
C         ENDIF
C         WRITE (IW,41)
C  41     FORMAT (1H0)
C         DO 408 I=1,IM
C 408     WRITE (IW,42) (A(I,J), J=1,IM)
C  42     FORMAT (1H ,F10.3,11F11.3/(4X,11F11.3))
C       ENDIF
C *********************************************** A.KRAMIDA  v
        IF (IPRNA.EQ.0) THEN
          IF (M.EQ.1) THEN
            WRITE (IW,40)
   40       FORMAT (/46H      T         ANG.COEFS FOR RADIAL INTEGRALS)
          ENDIF
          WRITE (IW,41) FJ(M)
   41     FORMAT (A8)
          DO I=1,IM
            IF (NF(I,M).LT.0) GO TO 420
            I1=A1(I,1)
            I2=A1(I,2)
            IF (I1.LE.0) GO TO 420
            WRITE (IW,4003)
 4003       FORMAT (1H )
            IF (I2.GT.IB+2) I2=IB+2
            WRITE (IW,4004) T(I,M),(A1(I,J),XIDEN(I1+J-3), J=3,3+I2-I1)
 4004       FORMAT (F10.3,3(2X,F9.4,1H ,A8)/
     1        1(10X,3(2X,F9.4,1H ,A8)))
  420       CONTINUE
          ENDDO
C
        ENDIF
C *********************************************** A.KRAMIDA  ^
C
C            CALCULATE EIGENVALUES AND VECTORS
C
C 500   CALL TRED2(IA,IM,A,EVAL(1,M),EV,EVEC)
C       CALL TQL2(IA,IM,EVAL(1,M),EV,EVEC,MM)
  500   if (IBY_SIZE.LT.NYOFS+IM**2) then
          new_size=max(2*IBY_SIZE,NYOFS+IM**2)
          CALL RALLOC(BUFY_PTR,IBY_SIZE,new_size,8)
          IBY_SIZE=new_size
        endif
        CALL TRED2(IA,IM,A,EVAL(1,M),EV,BUFY(NYOFS+1))
        CALL TQL2(IA,IM,EVAL(1,M),EV,BUFY(NYOFS+1),MM)
        SEVSQ(M)=0.D0
        DO  K=1,IM
          EVKM=EVAL(K,M)
          SEVSQ(M)=SEVSQ(M)+EVKM*EVKM
        ENDDO
        SEVSQ(M)=DSQRT(SEVSQ(M))
        IF (SEVSQ(M).EQ.0.D0.OR.ANORM(M).NE.0.D0.AND.
     *  DABS(SEVSQ(M)/ANORM(M)-1.D0).GT.0.000002D0) THEN
          WRITE (IW,50) M
   50     FORMAT ('0SEVSQ AND ANORM DISAGREE FOR M=',I2)
        ENDIF
C
C            READ G-VALUE MATRIX INTO C ARRAY
C
C ======================================== R.COWAN
C       CALL RD(15,IM,C,IA,A,IMSQ1)
C       READ (15)
C ======================================== R.COWAN
C
C *********************************************** A.KRAMIDA
        CALL RD1(M,LMAXP1,IM,C,IA,BUFC,IBC,IMSQ1,
     *           ICOFFS,MAXB1,MMAX,LMAXP1)
C *********************************************** A.KRAMIDA
C
        N=-10
        CALL ORDER(kcpl)                  ! A.Kramida 01/31/2011: before that there were no parameters to ORDER
        NYOFS=NYOFS+IMSQ         ! INCREASE THE Y-BUFFER OFFSET
C
C       CALL WR(I216,IM,Y,IA,C,IMSQ)       ! R.COWAN
C       DO 515 K=1,IM
C       DO 515 I=1,IM
C         MIK=MIK+1
C 515   BUFY(MIK)=Y(I,K)
        CALL CLOCK(TIME)
        DELT=TIME-TIME0
        IF (IW6.LT.0) WRITE (6,52) NOCYC,FJ(M),DELT
   52   FORMAT (' CYCLE NO.=',I2, '    FINISHED DIAGONALIZATION FOR',
     1    A8,'  AT TIME=',F8.4,' MIN')
  520 CONTINUE
      IPRNA=1                   ! A.KRAMIDA: A-MATRIX DOES NOT CHANGE.
C                               ! IT IS ENOUGH TO HAVE IT PRINTED ONCE.
C
      IF (TNO.NE.0.D0) THEN
        AVDEV=DSQRT(R/TNO)
      ELSE
        AVDEV=0.D0
        NOCYCE=0
      ENDIF
      IF (TNO.GT.PNO) THEN
        STDEV=DSQRT(R/(TNO-PNO))
      ELSE
        STDEV=0.D0
      ENDIF
      RETURN
C
      END
C
      SUBROUTINE ORDER(kcpl)                         ! A.Kramida 01/31/2011: before that there were no parameters
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
      DIMENSION KCPL(2)                              ! A.Kramida 01/31/2011: before that there were no parameters
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C ********************************************************v  A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC_PTR,BUFA_PTR,BUFY_PTR,
c     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
     1  IBC_PTR,IBA_PTR,ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
      INTEGER*8 BUFY_PTR,BUFA_PTR,BUFC_PTR,IBA_PTR,IBC_PTR
      POINTER (BUFY_PTR,BUFY(1)),(BUFA_PTR,BUFA(1)),(BUFC_PTR,BUFC(1))
      POINTER (IBA_PTR,IBA),(IBC_PTR,IBC)
      DIMENSION IBC(IX,IB),IBA(IX,1)
c
C       DETERMINE EIGENVECTORS ORDERING FROM FIRST MAXCR COMPONENTS
C       WITH LARGEST AMPLITUDES
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0,IEVV1,NW1
      INTEGER*2 II1,II2
      COMMON/VECT/ IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      integer*8 IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      DIMENSION IEVV(MAXC,IA,IX),IEVV0(MAXC,IA,IX),NW(MAXC,IA,IX),
     1 NW0(MAXC,IA,IX)
      POINTER(IEVV_PTR,IEVV),(IEVV0_PTR,IEVV0),(NW_PTR,NW),(NW0_PTR,NW0)
c
C ********************************************************^  A.KRAMIDA
C
      REAL*4 GFACT
      DIMENSION GFACT(IA),IDM(IA),IEVV1(MAXC,IA),NW1(MAXC,IA)
      SAVE GFACT,IDM,IEVV1,NW1
C ***************************************** A.Kramida 2026 v
C       WORK SPACE FOR THE FAST EIGENVECTOR MATCHER (SEE BELOW).
C       THE FORMER ICOMP(IA,IA) INTEGER*1 MATRIX (100 MB, ACCESSED WITH
C       A 10000-BYTE STRIDE ON THE INNERMOST TEST) IS REPLACED BY THE
C       BIT MAP IXCL: BIT (I-1)*IM+(J-1) IS SET WHEN J IS EXCLUDED FROM
C       THE COMPARISON LIST OF ROW I.  FOR A FIXED I THE BITS OF ALL J
C       ARE CONTIGUOUS, AND ONLY (IM*IM+31)/32 WORDS ARE TOUCHED.
      INTEGER*4 IXCL(IA*(IA/32+1))
C       INVERTED INDEX OVER THE OLD (SAMPLE) VECTORS:
C         IPSTA(N)   - START IN IPOST OF THE LIST OF OLD VECTORS HAVING
C                      BASIS STATE N AMONG THEIR N1F LARGEST COMPONENTS
C         IPOST(P)   - THE OLD-VECTOR INDICES THEMSELVES
C       IPOS(N)      - SCATTER MAP: POSITION OF BASIS STATE N AMONG THE
C                      COMPONENTS OF THE CURRENT NEW VECTOR I, ELSE 0
C       ISQ(J)       - SUM OF SQUARED COMPONENTS OF OLD VECTOR J
C       IDOT(J)      - SCALAR PRODUCT OF NEW VECTOR I WITH OLD VECTOR J
C       ICAND(1..NCAND) - THE OLD VECTORS WITH IDOT(J).NE.0
      INTEGER*4 IPSTA(IA+1),IPOST(IA*MAXCR),IPOSV(IA*MAXCR),IPOS(IA)
      INTEGER*4 IPOSF(IA*MAXCR)
      INTEGER*4 ISQ(IA),ICAND(IA),IDOT(IA),IEXT(IA)
C       BOUNDED MIN-HEAP USED WHEN SELECTING THE MAXC LARGEST COMPONENTS
      INTEGER*4 IHEAP(MAXC+1)
      REAL*8 AHEAP(MAXC+1)
C       COMPRESSED-COLUMN FORM OF THE NON-ZERO ELEMENTS OF THE G-MATRIX
      PARAMETER (MAXNZ=1000000)
      INTEGER*4 ICSTA(IA+1),ICROW(MAXNZ)
      REAL*8 CVAL(MAXNZ)
C       JFIX   - SLOTS CARRYING AN EXPERIMENTAL LEVEL (NF(J,M)=0)
C       KVUSED - INDEX OF THE FIXED SLOT THAT CLAIMED EIGENVECTOR I
C       KSLOT,KVEC - SLOTS AND EIGENVECTORS ORDERED BY ENERGY
      INTEGER*4 JFIX(IA),KVUSED(IA),KSLOT(IA),KVEC(IA)
      SAVE JFIX,KVUSED,KSLOT,KVEC
      SAVE IXCL,IPSTA,IPOST,IPOSV,IPOSF,IPOS,ISQ,ICAND,IDOT,IEXT
      SAVE ICSTA,ICROW,CVAL
C ***************************************** A.Kramida 2026 ^
C
C ***************************************** A.Kramida 2026 v
C       THE G-FACTOR OF EIGENVECTOR K IS THE QUADRATIC FORM
C            GFK = SUM(J) SUM(I) EVEC(I,K)*C(I,J)*EVEC(J,K)
C       WHICH USED TO BE EVALUATED BY A FULL IM*IM SWEEP FOR EACH OF THE
C       IM EIGENVECTORS - IM**3 OPERATIONS PER (PARITY,J) BLOCK - EVEN
C       THOUGH THE TEST IN THE INNER LOOP ALREADY SKIPPED THE ZEROS OF
C       C.  C IS IN FACT VERY SPARSE (STRICTLY DIAGONAL IN EVERY CASE
C       SEEN SO FAR), SO ITS NON-ZERO ELEMENTS ARE NOW COLLECTED ONCE
C       PER BLOCK IN COMPRESSED-COLUMN FORM AND ONLY THOSE ARE SWEPT,
C       WHICH COSTS IM*(IM+NNZ) OPERATIONS INSTEAD.
C       THE PRODUCTS ARE FORMED AND ADDED IN EXACTLY THE ORDER AND THE
C       ASSOCIATION USED BEFORE, SO THE G-FACTORS COME OUT BIT FOR BIT
C       THE SAME.  SHOULD C EVER BE TOO DENSE TO INDEX WITHIN MAXNZ,
C       THE ORIGINAL LOOPS ARE USED INSTEAD.
      NNZ=0
      DO 612 J=1,IM
        ICSTA(J)=NNZ+1
        DO 611 I=1,IM
          GRE=C(I,J)
          IF (GRE.EQ.0.D0) GO TO 611
          IF (NNZ.GE.MAXNZ) THEN
            NNZ=-1
            GO TO 613
          ENDIF
          NNZ=NNZ+1
          ICROW(NNZ)=I
          CVAL(NNZ)=GRE
  611   CONTINUE
  612 CONTINUE
      ICSTA(IM+1)=NNZ+1
  613 CONTINUE
      IF (NNZ.GE.0) THEN
        DO K=1,IM
          GFK=0.D0
          KOFS=NYOFS+(K-1)*IM
          DO J=1,IM
C           EJK=EVEC(J,K)
            EJK=BUFY(J+KOFS)
            DO IP=ICSTA(J),ICSTA(J+1)-1
C             GFK=GFK+EVEC(ICROW(IP),K)*CVAL(IP)*EJK
              GFK=GFK+BUFY(ICROW(IP)+KOFS)*CVAL(IP)*EJK
            ENDDO
          ENDDO
          GFAC(K,M)=GFK
          GFACT(K)=GFK
        ENDDO
      ELSE
        DO K=1,IM
          GFK=0.D0
          KOFS=NYOFS+(K-1)*IM
          DO J=1,IM
C           EJK=EVEC(J,K)
            EJK=BUFY(J+KOFS)
            DO I=1,IM
              GRE=C(I,J)
C             IF (GRE.NE.0.D0) GFK=GFK+EVEC(I,K)*GRE*EJK
              IF (GRE.NE.0.D0) GFK=GFK+BUFY(I+KOFS)*GRE*EJK
            ENDDO
          ENDDO
          GFAC(K,M)=GFK
          GFACT(K)=GFK
        ENDDO
      ENDIF
C ***************************************** A.Kramida 2026 ^
C
C ========================================= R.COWAN v
C     IF (N.GT.0) GO TO 670   ! OUTPUT PHASE IS MOVED TO ORDOUT PROCEDURES
C
C          SET TO ORDER EIGENVECTORS BY LARGEST COMPONENT
C
C     DO 603 I=1,IM
C 603 KW(I,M)=0
C 604 DO 614 J=1,IM
C       GRE=0.D0
C       DO 606 I=1,IM
C         EV(I,1)=EVEC(I,J)**2
C         IF (EV(I,1).GT.GRE) THEN
C           GRE=EV(I,1)
C           K=I
C         ENDIF
C 606   CONTINUE
C       DO 608 I=1,IM
C         IF (I.NE.K.AND.EV(I,1)/GRE.GE.CRIT) GO TO 614
C 608   CONTINUE
C       IF (GRE.GT.0.5D0) THEN
C         KW(K,M)=J
C       ELSE
C         KWKM=KW(K,M)
C         IF (KWKM.LE.0) THEN
C           IF (KWKM.LT.0) THEN
C             KWKM=-100
C           ELSE
C             KWKM=-J
C           ENDIF
C           KW(K,M)=KWKM
C         ENDIF
C       ENDIF
C 614 CONTINUE
C     DO 622 K=1,IM
C       IF (KW(K,M).LE.-100) THEN
C         KW(K,M)=0
C       ELSE
C         KW(K,M)=IABS(KW(K,M))
C       ENDIF
C 622 CONTINUE
C
C          FOR EIGENVECTORS WITH NO DOMINANT COMPONENT,
C              SET TO ORDER BY MAGNITUDE OF EVAL (N NEG) OR GFAC (N POS)
C
C     N1I=1
C 631 DO 633 N1=N1I,IM
C       IF (NC(N1,M).EQ.NC(N1I,M)) N1F=IM
C 633 CONTINUE
C 631 N1F=IM
C 634 ICNTR=0
C     GRE=-1.0D30
C     DO 641 N1=N1I,N1F
C       IF (KW(N1,M).EQ.0) THEN
C         ICNTR=ICNTR+1
C         IF (N.GT.0) THEN
C           IF (GKK(N1,M).GE.GRE) THEN
C             GRE=GKK(N1,M)
C             K=N1
C           ENDIF
C         ELSE
C           IF (T(N1,M).GE.GRE) THEN
C             GRE=T(N1,M)
C             K=N1
C           ENDIF
C         ENDIF
C       ENDIF
C 641 CONTINUE
C     IF (ICNTR.NE.0) THEN
C       GRE=-1.0D30
C       DO 649 N1=N1I,N1F
C         DO 643 I=N1I,N1F
C           IF (IABS(KW(I,M)).EQ.N1) GO TO 649
C 643     CONTINUE
C         IF (N.LE.0) THEN
C           IF (EVAL(N1,M).GT.GRE) THEN
C             GRE=EVAL(N1,M)
C             J=N1
C           ENDIF
C         ELSE
C           IF (GFAC(N1,M).GT.GRE) THEN
C             GRE=GFAC(N1,M)
C             J=N1
C           ENDIF
C         ENDIF
C 649   CONTINUE
C       KW(K,M)=-J
C       GO TO 634
C     ENDIF
C 650 N1I=N1F+1
C     IF (N1I.LE.IM) GO TO 631
C ========================================= R.COWAN  ^
C
C ***************************************** A.KRAMIDA v
C       MAP EIGENVECTORS INTO IEVV(I,J,M) AND NW(I,J,M), J=1..MAXC.
C       K=NW(I,J,M) KEEPS THE REFERENCE TO EVEC(I,K) FOR CURRENT M:
C          IEVV(I,J,M)=EVEC(I,K)
C       IEVV(I,J,M) IS SORTED BY AMPLITUDE: IEVV(I,J,M)>=IEVV(I,J+1,M)
C
C       A.Kramida 2026: THE N1F LARGEST AMPLITUDES USED TO BE PICKED BY
C       N1F SUCCESSIVE SEARCHES THROUGH ALL IM AMPLITUDES (IM*IM*MAXC
C       OPERATIONS PER (PARITY,J) BLOCK).  THEY ARE NOW COLLECTED IN A
C       SINGLE PASS WITH A BOUNDED MIN-HEAP OF SIZE N1F AND HEAP-SORTED
C       INTO DESCENDING ORDER (IM*(IM+MAXC*LOG(MAXC)) OPERATIONS).
C       THE ORDERING IS THE SAME TOTAL ORDER AS BEFORE - DECREASING
C       ABSOLUTE AMPLITUDE, TIES BROKEN IN FAVOUR OF THE LOWER BASIS
C       INDEX - AND ZERO AMPLITUDES ARE STILL NEVER SELECTED, SO THE
C       CONTENTS OF NW AND IEVV ARE UNCHANGED.
C       AHEAP(P) HOLDS |AMPLITUDE|, IHEAP(P) THE BASIS INDEX.  THE HEAP
C       ROOT IS THE SMALLEST ELEMENT UNDER THAT ORDER.
      N1F=MIN0(IM,MAXC)
      DO N2=1,IM
        N2OFS=NYOFS+(N2-1)*IM
        NH=0
        DO N1=1,IM
C         EV1=DABS(EVEC(N1,N2))
          EV1=DABS(BUFY(N1+N2OFS))
          IF (EV1.LE.0.D0) GO TO 626
          IF (NH.GE.N1F) THEN
C             KEEP THE COMPONENT ONLY IF IT BEATS THE WEAKEST ONE HELD
            IF (EV1.LT.AHEAP(1)) GO TO 626
            IF (EV1.EQ.AHEAP(1).AND.N1.GT.IHEAP(1)) GO TO 626
            AHEAP(1)=EV1
            IHEAP(1)=N1
            NHP=1
            GO TO 624                 ! SIFT THE NEW ROOT DOWN
          ENDIF
C           HEAP NOT FULL YET - APPEND AND SIFT UP
          NH=NH+1
          AHEAP(NH)=EV1
          IHEAP(NH)=N1
          NHC=NH
  623     IF (NHC.GT.1) THEN
            NHP=NHC/2
            IF (AHEAP(NHC).LT.AHEAP(NHP).OR.(AHEAP(NHC).EQ.AHEAP(NHP)
     1          .AND.IHEAP(NHC).GT.IHEAP(NHP))) THEN
              EVT=AHEAP(NHC)
              AHEAP(NHC)=AHEAP(NHP)
              AHEAP(NHP)=EVT
              IHT=IHEAP(NHC)
              IHEAP(NHC)=IHEAP(NHP)
              IHEAP(NHP)=IHT
              NHC=NHP
              GO TO 623
            ENDIF
          ENDIF
          GO TO 626
C           SIFT AHEAP(NHP) DOWN TO ITS PLACE WITHIN AHEAP(1..NH)
  624     NHC=2*NHP
          IF (NHC.GT.NH) GO TO 626
          IF (NHC.LT.NH) THEN
            IF (AHEAP(NHC+1).LT.AHEAP(NHC).OR.(AHEAP(NHC+1).EQ.
     1          AHEAP(NHC).AND.IHEAP(NHC+1).GT.IHEAP(NHC))) NHC=NHC+1
          ENDIF
          IF (AHEAP(NHC).LT.AHEAP(NHP).OR.(AHEAP(NHC).EQ.AHEAP(NHP)
     1        .AND.IHEAP(NHC).GT.IHEAP(NHP))) THEN
            EVT=AHEAP(NHC)
            AHEAP(NHC)=AHEAP(NHP)
            AHEAP(NHP)=EVT
            IHT=IHEAP(NHC)
            IHEAP(NHC)=IHEAP(NHP)
            IHEAP(NHP)=IHT
            NHP=NHC
            GO TO 624
          ENDIF
  626     CONTINUE
        ENDDO
C           HEAP-SORT: MOVE THE ROOT (SMALLEST) TO THE BACK REPEATEDLY,
C           WHICH LEAVES AHEAP(1..NH) IN DECREASING ORDER
        NHS=NH
  627   IF (NHS.GT.1) THEN
          EVT=AHEAP(1)
          AHEAP(1)=AHEAP(NHS)
          AHEAP(NHS)=EVT
          IHT=IHEAP(1)
          IHEAP(1)=IHEAP(NHS)
          IHEAP(NHS)=IHT
          NHS=NHS-1
          NHP=1
  628     NHC=2*NHP
          IF (NHC.LE.NHS) THEN
            IF (NHC.LT.NHS) THEN
              IF (AHEAP(NHC+1).LT.AHEAP(NHC).OR.(AHEAP(NHC+1).EQ.
     1            AHEAP(NHC).AND.IHEAP(NHC+1).GT.IHEAP(NHC))) NHC=NHC+1
            ENDIF
            IF (AHEAP(NHC).LT.AHEAP(NHP).OR.(AHEAP(NHC).EQ.AHEAP(NHP)
     1          .AND.IHEAP(NHC).GT.IHEAP(NHP))) THEN
              EVT=AHEAP(NHC)
              AHEAP(NHC)=AHEAP(NHP)
              AHEAP(NHP)=EVT
              IHT=IHEAP(NHC)
              IHEAP(NHC)=IHEAP(NHP)
              IHEAP(NHP)=IHT
              NHP=NHC
              GO TO 628
            ENDIF
          ENDIF
          GO TO 627
        ENDIF
        DO K=1,NH
          K1=IHEAP(K)
          NW(K,N2,M)=K1 ! -128 FOR 1-BYTE INTEGER STORAGE
C         IEVV(K,N2,M)=100.D0*EVEC(K1,N2)+0.5
C ***************************************** A.Kramida 2026 v
C           ROUND TO NEAREST, NOT TOWARDS +INFINITY.  THE CONVERSION TO
C           INTEGER TRUNCATES TOWARDS ZERO, SO A FIXED +0.5D0 ROUNDED
C           POSITIVE AMPLITUDES CORRECTLY BUT MOVED NEGATIVE ONES ONE
C           UNIT TOWARDS ZERO (-0.005450 CAME OUT AS -54 INSTEAD OF
C           -55).  DSIGN GIVES THE OFFSET THE SIGN OF THE AMPLITUDE, SO
C           BOTH SIGNS ROUND BY MAGNITUDE.  POSITIVE AMPLITUDES KEEP THE
C           IDENTICAL EXPRESSION AND HENCE THE IDENTICAL RESULT.
          EVAMP=1.D4*BUFY(K1+N2OFS)
          IEVV(K,N2,M)=EVAMP+DSIGN(0.5D0,EVAMP)
C ***************************************** A.Kramida 2026 ^
        ENDDO
        DO I=NH+1,MAXC
          IEVV(I,N2,M)=0
          NW(I,N2,M)=0 ! -128 FOR 1-BYTE INTEGER STORAGE
        ENDDO
      ENDDO
      N1F=MIN0(IM,MAXCR)
C ***************************************** A.Kramida 2026 v
C       COLLECT THE SLOTS THAT CARRY AN EXPERIMENTAL LEVEL.  NF IS SET
C       WHEN THE INPUT IS READ AND IS NEVER PERMUTED BY THE REORDERING,
C       SO SLOT J ALWAYS MEANS THE SAME OBSERVED LEVEL T(J,M).
      NFIXL=0
      DO J=1,IM
        IF (NF(J,M).EQ.0) THEN
          NFIXL=NFIXL+1
          JFIX(NFIXL)=J
        ENDIF
      ENDDO
C
C       SELECT THE ORDERING METHOD.  INVECT COMES FROM THE FIFTH FIELD
C       OF THE FIRST RECORD OF 'OUTGINE' (SEE THE MAIN PROGRAM):
C         0,1 - FULL EIGENVECTOR RECOGNITION FOR EVERY LEVEL
C         2   - RELAXED: EIGENVECTOR RECOGNITION ONLY FOR THE LEVELS
C               THAT ARE MATCHED TO AN EXPERIMENTAL ENERGY, THE REST
C               ORDERED BY INCREASING ENERGY
C         3   - EIGENVECTORS IGNORED, EVERYTHING ORDERED BY ENERGY
C       A BLOCK WITHOUT ANY EXPERIMENTAL LEVEL HAS NOTHING TO RECOGNIZE
C       AGAINST, SO UNDER INVECT=2 IT IS ORDERED BY ENERGY AT ONCE.
      IESORT=0
c               ! Changed by A. Kramida 01/31/2011
c     IF ((NOCYC.EQ.0.AND.INVECT.EQ.0).OR.NOCYC.LT.IORD) THEN
      IF ((NOCYC.EQ.0.AND.(INVECT.EQ.0.or.kcpl(1).ne.1)).
     +    OR.NOCYC.LT.IORD) IESORT=1
c               ! End changes by A. Kramida 01/31/2011
      IF (INVECT.GE.3) IESORT=1
      IF (INVECT.EQ.2.AND.NFIXL.EQ.0) IESORT=1
      IF (IESORT.NE.0) THEN
C           NO REORDERING!
C           THE SLOT HOLDING THE R-TH LARGEST INPUT ENERGY TAKES THE
C           R-TH LARGEST EIGENVALUE.  THE ORIGINAL CODE OBTAINED THIS
C           BY IM SUCCESSIVE SEARCHES, EACH CONTAINING A NESTED SEARCH
C           FOR THE ALREADY USED EIGENVALUES - O(IM**3) ALTOGETHER,
C           WHICH IS PROHIBITIVE FOR LARGE BLOCKS.  TWO INDEX SORTS
C           GIVE THE SAME ASSIGNMENT, TIES INCLUDED.   A.Kramida 2026
        DO K=1,IM
          KSLOT(K)=K
          KVEC(K)=K
        ENDDO
        CALL IDXSRT(IM,KSLOT,T(1,M),1)
        CALL IDXSRT(IM,KVEC,EVAL(1,M),0)
        DO K=1,IM
          KW(KSLOT(K),M)=KVEC(K)
        ENDDO
      ELSE IF (INVECT.EQ.2) THEN
C
C           RELAXED MATCHING.  ONLY THE NFIXL SLOTS THAT CARRY AN
C           EXPERIMENTAL LEVEL ARE MATCHED THROUGH THEIR EIGENVECTORS,
C           WHICH ARE THE ONLY ONES WHOSE IDENTITY THE FIT DEPENDS ON.
C           THE COST DROPS FROM O(IM**2) DISTANCES TO O(NFIXL*IM), AND
C           THE BACKTRACKING SHRINKS IN THE SAME PROPORTION.  THE LEVELS
C           WITHOUT AN OBSERVED ENERGY ARE THEN ORDERED BY INCREASING
C           ENERGY; THEIR EIGENVECTOR LABELS NO LONGER TRACK A LEVEL
C           FROM CYCLE TO CYCLE, WHICH IS OF NO CONSEQUENCE TO THE FIT.
        DO K=1,IM
          KW(K,M)=0
          KVUSED(K)=0
        ENDDO
        NXWRD=(NFIXL*IM+31)/32        ! CLEAR THE EXCLUSION BIT MAP
        DO K=1,NXWRD
          IXCL(K)=0
        ENDDO
        IFX=0
  671   IF (IFX.EQ.NFIXL) GO TO 680
        IFX=IFX+1
  672   IDM(IFX)=1000000000
        IMIN=0
        J=JFIX(IFX)
        IXROW=(IFX-1)*IM              ! BIT OFFSET OF ROW IFX
        DO I=1,IM
          IXB=IXROW+I-1
          IF (BTEST(IXCL(IXB/32+1),IAND(IXB,31))) GO TO 674
          CALL EVDIST(I,J,M,N1F,IDD)
          IF (IDD.LT.IDM(IFX)) THEN
            IMIN=I
            IDM(IFX)=IDD
C             A DISTANCE OF 0.10 ALREADY IDENTIFIES THE EIGENVECTOR
            IF (IDM(IFX).LT.1000000) GO TO 675
          ENDIF
  674     CONTINUE
        ENDDO
  675   IF (IMIN.EQ.0) THEN
          WRITE (IW,*) ' Failed to identify eigenvectors in RCEINP.'
          STOP '(Fix some parameters or change MAXC and recompile.)'
        ENDIF
        IFK=KVUSED(IMIN)
        IF (IFK.NE.0) THEN
C             EIGENVECTOR IMIN IS ALREADY CLAIMED - THE CLOSER SLOT
C             KEEPS IT AND THE OTHER ONE LOOKS AGAIN
          IF (IDM(IFX).GT.IDM(IFK)) THEN
            IXB=(IFX-1)*IM+IMIN-1
          ELSE
            IXB=(IFK-1)*IM+IMIN-1
            DO K=IFK,NFIXL
              IF (KW(JFIX(K),M).NE.0) THEN
                KVUSED(KW(JFIX(K),M))=0
                KW(JFIX(K),M)=0
              ENDIF
            ENDDO
            IFX=IFK
          ENDIF
          IXCL(IXB/32+1)=IBSET(IXCL(IXB/32+1),IAND(IXB,31))
          GO TO 672
        ENDIF
        KVUSED(IMIN)=IFX
        KW(J,M)=IMIN
        GO TO 671
  680   CONTINUE
C           ORDER WHAT IS LEFT BY INCREASING ENERGY
        NS=0
        NV=0
        DO K=1,IM
          IF (KW(K,M).EQ.0) THEN
            NS=NS+1
            KSLOT(NS)=K
          ENDIF
          IF (KVUSED(K).EQ.0) THEN
            NV=NV+1
            KVEC(NV)=K
          ENDIF
        ENDDO
        CALL IDXSRT(NS,KSLOT,T(1,M),1)
        CALL IDXSRT(NV,KVEC,EVAL(1,M),0)
        DO K=1,NS
          KW(KSLOT(K),M)=KVEC(K)
        ENDDO
C ***************************************** A.Kramida 2026 ^
      ELSE
C       IDENTIFY THE NEW EIGENVECTORS WITH THE OLD ONES
        DO J=1,IM
          IDM(J)=0
          KW(J,M)=0
        ENDDO
        NXWRD=(IM*IM+31)/32           ! CLEAR THE EXCLUSION BIT MAP
        DO J=1,NXWRD
          IXCL(J)=0
        ENDDO
C ***************************************** A.Kramida 2026 v
C       PREPARATION FOR THE FAST DISTANCE EVALUATION.
C
C       WRITING SI = SUM(K) IEVV(K,I,M)**2, SJ = SUM(K) IEVV0(K,J,M)**2
C       AND DOT = SUM OVER THE COMPONENTS SHARED BY I AND J OF THE
C       PRODUCT IEVV(.,I,M)*IEVV0(.,J,M), THE DISTANCE COMPUTED BY THE
C       ORIGINAL DOUBLE LOOP IS EXACTLY
C            IDD = SI + SJ - 2*DOT   (JS=1)
C            IDD = SI + SJ + 2*DOT   (JS=2)
C       SO BOTH SIGN VARIANTS FOLLOW FROM ONE SCALAR PRODUCT, AND DOT
C       VANISHES FOR EVERY J THAT SHARES NO COMPONENT WITH I.  THE FEW
C       J THAT DO SHARE ONE ARE FOUND THROUGH AN INVERTED INDEX BUILT
C       ONCE PER (PARITY,J) BLOCK, SO THE INNER K,K1 LOOPS DISAPPEAR
C       AND EVERY J IS STILL VISITED IN THE ORIGINAL ORDER.  THE RESULT,
C       INCLUDING TIES AND THE EARLY EXIT BELOW, IS BIT-IDENTICAL.
C
C       A ZERO IN NW0(1..N1F,J,M) - WHICH OCCURS WHEN A COMPONENT READ
C       FROM 'RCEINP' CARRIES A TERM NAME THAT IS NOT RECOGNIZED -
C       TERMINATES THE SEARCH THROUGH THE COMPONENTS OF J IN THE FIRST
C       OF THE TWO ORIGINAL LOOPS BUT NOT IN THE SECOND.  COMPONENTS OF
C       J BEYOND THAT ZERO THEREFORE CANNOT PAIR UP WITH A COMPONENT OF
C       I, YET THEY ARE STILL DROPPED FROM THE SECOND SUM WHEN THEY
C       HAPPEN TO OCCUR IN I.  THAT IS ACCOUNTED FOR EXACTLY BY THE
C       CORRECTION IEXT(J), SO THE DISTANCE BECOMES
C            IDD = SI + SJ -+ 2*IDOT(J) - IEXT(J)
C       WITH IDOT RESTRICTED TO THE COMPONENTS OF J AHEAD OF ITS FIRST
C       ZERO, AND IEXT SUMMING THE SQUARED COMPONENTS BEYOND IT THAT
C       APPEAR AMONG THE COMPONENTS OF I.
C
C       WHAT THE IDENTITY STILL REQUIRES IS THAT NO BASIS INDEX BE
C       REPEATED (THE ORIGINAL LOOPS WOULD PAIR UP ONLY ITS FIRST
C       OCCURRENCE) AND THAT NW(1..N1F,.,M), WHICH IS BUILT ABOVE AND
C       IS SEARCHED WITHOUT TRUNCATION, HOLD NO ZERO.  IN THE UNLIKELY
C       CASE THAT EITHER FAILS, THE ORIGINAL LOOPS ARE USED INSTEAD.
        IFAST=1
        DO J=1,IM
          DO K=1,N1F
            N1=NW(K,J,M)
            IF (N1.EQ.0) IFAST=0
            DO K1=1,K-1
              IF (NW(K1,J,M).EQ.N1) IFAST=0
            ENDDO
          ENDDO
          DO K=1,N1F
            N2=NW0(K,J,M)
            IF (N2.EQ.0) GO TO 619       ! REST OF THE LIST IS UNPAIRABLE
            DO K1=1,K-1
              IF (NW0(K1,J,M).EQ.N2) IFAST=0
            ENDDO
          ENDDO
  619     CONTINUE
        ENDDO
        IF (IFAST.NE.0) THEN
C           SUMS OF SQUARED COMPONENTS OF THE OLD (SAMPLE) VECTORS
          DO J=1,IM
            IS=0
            DO K=1,N1F
              IAV=IEVV0(K,J,M)
              IS=IS+IAV*IAV
            ENDDO
            ISQ(J)=IS
          ENDDO
C           INVERTED INDEX: FOR EACH BASIS STATE, THE OLD VECTORS THAT
C           CARRY IT AMONG THEIR N1F LARGEST COMPONENTS
          DO N1=1,IM
            IPSTA(N1)=0
          ENDDO
          DO J=1,IM
            DO K=1,N1F
              N1=NW0(K,J,M)
              IF (N1.NE.0) IPSTA(N1)=IPSTA(N1)+1
            ENDDO
          ENDDO
          NP=1
          DO N1=1,IM
            NPC=IPSTA(N1)
            IPSTA(N1)=NP
            IPOS(N1)=NP
            NP=NP+NPC
          ENDDO
          IPSTA(IM+1)=NP
C           IPOSF(P)=0 MARKS A COMPONENT AHEAD OF THE FIRST ZERO OF ITS
C           LIST (IPOSV = THE AMPLITUDE, USED FOR IDOT), IPOSF(P)=1 ONE
C           BEYOND IT (IPOSV = THE SQUARED AMPLITUDE, USED FOR IEXT)
          DO J=1,IM
            IPRE=0
            DO K=1,N1F
              N1=NW0(K,J,M)
              IF (N1.EQ.0) THEN
                IPRE=1
                GO TO 620
              ENDIF
              IP=IPOS(N1)
              IPOST(IP)=J
              IAV=IEVV0(K,J,M)
              IPOSF(IP)=IPRE
              IF (IPRE.EQ.0) THEN
                IPOSV(IP)=IAV
              ELSE
                IPOSV(IP)=IAV*IAV
              ENDIF
              IPOS(N1)=IP+1
  620         CONTINUE
            ENDDO
          ENDDO
          DO N1=1,IM
            IPOS(N1)=0
            IDOT(N1)=0
            IEXT(N1)=0
          ENDDO
        ENDIF
        ITAG=0
        NCAND=0
C ***************************************** A.Kramida 2026 ^
        I=0
c        DO 660 I=1,IM
  641   IF (I.EQ.IM) GO TO 660
        I=I+1
  642   IDM(I)=1000000000
        JMIN=0
        IXROW=(I-1)*IM                ! BIT OFFSET OF ROW I
C ***************************************** A.Kramida 2026 v
C         ACCUMULATE THE SCALAR PRODUCTS OF THE NEW VECTOR I WITH THE
C         OLD VECTORS SHARING AT LEAST ONE OF ITS N1F COMPONENTS
        IF (IFAST.NE.0) THEN
          DO IP=1,NCAND
            IDOT(ICAND(IP))=0
            IEXT(ICAND(IP))=0
          ENDDO
          NCAND=0
          ITAG=ITAG+1
          ISQI=0
          DO K=1,N1F
            IAV=IEVV(K,I,M)
            ISQI=ISQI+IAV*IAV
            N1=NW(K,I,M)
            DO IP=IPSTA(N1),IPSTA(N1+1)-1
              J=IPOST(IP)
              IF (IPOS(J).NE.ITAG) THEN
                IPOS(J)=ITAG
                NCAND=NCAND+1
                ICAND(NCAND)=J
              ENDIF
              IF (IPOSF(IP).EQ.0) THEN
                IDOT(J)=IDOT(J)+IAV*IPOSV(IP)
              ELSE
                IEXT(J)=IEXT(J)+IPOSV(IP)
              ENDIF
            ENDDO
          ENDDO
        ENDIF
C ***************************************** A.Kramida 2026 ^
C
C         FIND IEVV0(J,...,M) CLOSEST TO IEVV(I,...,M)
        IFAIL=0
        DO JS=1,2
          ISGN=2*(2*JS-3)             ! -2 FOR JS=1, +2 FOR JS=2
          DO J=1,IM
            IXB=IXROW+J-1
            IF (BTEST(IXCL(IXB/32+1),IAND(IXB,31))) GO TO 652 ! SKIP COMPARISOM
C             CALC DISTANCE OF VECTOR IEVV0(I,...,M) FROM IEVV(J,...,M)
C             IN IM-DIMENSIONAL SPACE
C
            IF (IFAST.NE.0) THEN
              IDD=ISQI+ISQ(J)+ISGN*IDOT(J)-IEXT(J)
              GO TO 651
            ENDIF
            IDD=0
            NFOUND=0
            DO K=1,N1F
              N1=NW(K,I,M)
              IA0=IEVV(K,I,M)
              N2=0
              DO K1=1,N1F
                NK1=NW0(K1,J,M)
                IF (NK1.EQ.0) GOTO 646   ! -128 FOR INTEGER*1 NW,NW0
                IF (NK1.EQ.N1) THEN
                  N2=K1
                  NFOUND=NFOUND+1
                  GOTO 646
                ENDIF
              ENDDO
  646         IA1=0
              IF (N2.NE.0) IA1=IEVV0(N2,J,M)
              IF (JS.EQ.2) IA1=-IA1     ! CHECK FOR INVERTED SIGN
              IA10=IA1-IA0
              IDD=IDD+IA10*IA10
            ENDDO
            IF (NFOUND.LT.N1F) THEN
C                  ADD COMPONENTS OF NW0(J) MISSING IN NW(I)
              DO K=1,N1F
                N1=NW0(K,J,M)
                IA0=IEVV0(K,J,M)
                N2=0
                DO K1=1,N1F
                  NK1=NW(K1,I,M)
                  IF (NK1.EQ.0) GOTO 649   ! -128 FOR INTEGER*1 NW,NW0
                  IF (NK1.EQ.N1) GOTO 650  ! SKIP COMPONENT PROCESSED EARLIER
                ENDDO
  649           IDD=IDD+IA0*IA0
  650           CONTINUE
              ENDDO
            ENDIF
  651       IF (IDD.LT.IDM(I)) THEN
              JMIN=J
              IDM(I)=IDD
C          DISTANCE OF 0.10 IS WELL ENOUGH TO IDENTIFY THE EIGENVECTORS.
C          IT CORRESPONDS TO VALUE OF PROJECTION (SCALAR MULT.) OF 0.95
              IF (IDM(I).LT.1000000) GO TO 653
            ENDIF
  652       CONTINUE
          ENDDO
        ENDDO
        IF (JMIN.EQ.0) IFAIL=1    ! COMPARISON LIST IS EMPTY FOR I
  653   IF (IFAIL.EQ.0) THEN
          KWJM=KW(JMIN,M)
          JK=IABS(KWJM)
          IF (JK.NE.0) THEN
C           IF (ISORT.EQ.0) THEN
C            COMPARE DISTANCES IDM(I) AND IDM(JK)
              ID1=IDM(JK)
              IF (IDM(I).GT.ID1) THEN
C                 EXCLUDE JMIN FROM COMPARE LIST FOR CURRENT I
C                 AND REDO THE CYCLE FOR CURRENT I
                IEXR=I
              ELSE
C                 EXCLUDE JMIN FROM COMPARE LIST FOR I=JK
C                 AND RETURN BACK TO I=JK.
C                 NOTE THAT JK CAN NOT BE GREATER THAN I
                IEXR=JK
                DO J=1,IM
                  KWJM=KW(J,M)
                  IF (KWJM.NE.0.AND.KWJM.GE.JK) KW(J,M)=0
                ENDDO
                I=JK
              ENDIF
              IXB=(IEXR-1)*IM+JMIN-1
              IXCL(IXB/32+1)=IBSET(IXCL(IXB/32+1),IAND(IXB,31))
              GO TO 642
C           ELSE
C             KW(JMIN,M)=-JK
C           ENDIF
          ENDIF
        ENDIF
        IF (IFAIL.NE.0) THEN
C      
C             IDENTIFICATION FAILED
C      
          WRITE (IW,*) ' Failed to identify eigenvectors in RCEINP.'
          STOP '(Fix some parameters or change MAXC and recompile.)'
        ENDIF
        IF (KW(JMIN,M).GE.0) KW(JMIN,M)=I
        GO TO 641
  660   CONTINUE
C       IF (ISORT.NE.0) THEN
C            FOR EIGENVECTORS NOT MAPPED AT THIS POINT,
C            SET TO ORDER BY MAGNITUDE OF EVAL (N NEG) OR GFAC (N POS)
C
C         DO I=1,IM
C           IF (KW(I,M).LT.0) KW(I,M)=0
C         ENDDO
C 662     ICNTR=0
C         GRE=-1.0D30
C         DO N1=1,IM
C             FIND THE GREATEST EIGENVALUE (NUMBER K) NOT MAPPED YET
C           IF (KW(N1,M).EQ.0) THEN
C             ICNTR=ICNTR+1
C             IF (N.GT.0) THEN
C               IF (GKK(N1,M).GE.GRE) THEN
C                 GRE=GKK(N1,M)
C                 K=N1
C               ENDIF
C             ELSE
C               IF (T(N1,M).GE.GRE) THEN
C                 GRE=T(N1,M)
C                 K=N1
C               ENDIF
C             ENDIF
C           ENDIF
C         ENDDO
C         IF (ICNTR.NE.0) THEN
C           GRE=-1.0D30
C           DO 665 N1=1,IM
C             DO I=1,IM
C               IF (KW(I,M).EQ.N1) GO TO 665
C             ENDDO
C             IF (N.LE.0) THEN
C               IF (EVAL(N1,M).GT.GRE) THEN
C                 GRE=EVAL(N1,M)
C                 J=N1
C               ENDIF
C             ELSE
C               IF (GFAC(N1,M).GT.GRE) THEN
C                 GRE=GFAC(N1,M)
C                 J=N1
C               ENDIF
C             ENDIF
C 665       CONTINUE
C           KW(K,M)=J
C           GO TO 662
C         ENDIF
C       ENDIF
      ENDIF
C ***************************************** A.KRAMIDA ^
C
C          REARRANGE EIGENVALUES, VECTORS, AND GFACTORS
C
      N1F=MIN0(IM,MAXC)               ! A.KRAMIDA
      DO K=1,IM
        KWKM=KW(K,M)
        J=IABS(KWKM)
        W(K,M)=EVAL(J,M)
        GFAC(K,M)=GFACT(J)
C ***************************************** A.KRAMIDA v
        DO N1=1,N1F
          IEVV1(N1,K)=IEVV(N1,J,M)
          NW1(N1,K)=NW(N1,J,M)
        ENDDO
C ***************************************** A.KRAMIDA ^
        IF (NF(K,M).EQ.0) THEN
          TW(K,M)=T(K,M)-W(K,M)
          GRE=TW(K,M)
          R=R+GRE*GRE
          TNO=TNO+1.D0
        ENDIF
      ENDDO
C           REORDER EIGENVECTORS DIRECTLY IN THE Y-BUFFER
      DO K=1,IM
        KWKM=KW(K,M)
        J=IABS(KWKM)
        IF (J.NE.K) THEN
          KOFS=NYOFS+(K-1)*IM
          JOFS=NYOFS+(J-1)*IM
          DO I=1,IM
            GRE=BUFY(I+KOFS)
            BUFY(I+KOFS)=BUFY(I+JOFS)
            BUFY(I+JOFS)=GRE
          ENDDO
          ICTRL=0
          DO I=K,IM
            IF (KW(I,M).EQ.K) THEN
              ICTRL=1
              KW(I,M)=J
              GO TO 655
            ENDIF
          ENDDO
655       CONTINUE
          IF (ICTRL.EQ.0) THEN
            STOP '**** EGENVECTORS REORDERING FAILED. ****'
          ENDIF
        ENDIF
      ENDDO
      DO K=1,IM
        DO N1=1,N1F
          IEVV(N1,K,M)=IEVV1(N1,K)
          NW(N1,K,M)=NW1(N1,K)
        ENDDO
      ENDDO
      IF (NOCYC.EQ.IORD.OR.(NOCYC.GT.IORD.AND.ISAMPL.EQ.1)) THEN
        DO K=1,IM
          DO N1=1,N1F
            IEVV0(N1,K,M)=IEVV(N1,K,M)
            NW0(N1,K,M)=NW(N1,K,M)
          ENDDO
        ENDDO
      ENDIF
      RETURN
      END
C
      SUBROUTINE ORDOUT1
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C ********************************************************v  A.KRAMIDA
C       DETERMINE EIGENVECTORS ORDERING FROM FIRST MAXC COMPONENTS
C       WITH LARGEST AMPLITUDES
      INTEGER*2 NW,NW0,KW1
      INTEGER*2 IEVV,IEVV0
      INTEGER*2 II1,II2
      COMMON/VECT/ IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      integer*8 IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      DIMENSION IEVV(MAXC,IA,IX),IEVV0(MAXC,IA,IX),NW(MAXC,IA,IX),
     1 NW0(MAXC,IA,IX)
      POINTER(IEVV_PTR,IEVV),(IEVV0_PTR,IEVV0),(NW_PTR,NW),(NW0_PTR,NW0)
c
C ********************************************************^  A.KRAMIDA
      DIMENSION IDM(IA),KW1(IA)
      SAVE IDM,KW1
C       A.Kramida 2026: ICOMP(IA,IA) REPLACED BY THE BIT MAP IXCL.
      INTEGER*4 IXCL(IA*(IA/32+1))
      SAVE IXCL
C
      N1F=MIN0(IM,MAXC)
      DO J=1,IM
        NW0(1,J,M)=J
        IEVV0(1,J,M)=10000
        IDM(J)=0
        KW1(J)=0
        DO NM1=2,N1F
          NW0(NM1,J,M)=0
          IEVV0(NM1,J,M)=0
        ENDDO
      ENDDO
C       CLEAR THE EXCLUSION BIT MAP (REPLACES ICOMP).   A.Kramida 2026
      NXWRD=(IM*IM+31)/32
      DO J=1,NXWRD
        IXCL(J)=0
      ENDDO
C       IDENTIFY THE EIGENVECTORS WITH THE SAMPLE ONES
      I=0
  641 IF (I.EQ.IM) GO TO 660
      I=I+1
  642 IDM(I)=1000000000
      JMIN=0
      IXROW=(I-1)*IM                  ! BIT OFFSET OF ROW I
C
C       FIND IEVV0(J,...,M) CLOSEST TO IEVV(I,...,M)
      IFAIL=0
      DO JS=1,2
        DO J=1,IM
          IXB=IXROW+J-1
          IF (BTEST(IXCL(IXB/32+1),IAND(IXB,31))) GO TO 652 ! SKIP COMPARISOM
C         CALC DISTANCE OF VECTOR IEVV0(I,...,M) FROM IEVV(J,...,M)
C         IN IM-DIMENSIONAL SPACE
C
          IDD=0
          NFOUND=0
C          NK1=NW0(1,J,M)
          IA11=IEVV0(1,J,M)
          IAM=IEVV(1,I,M)
          DO K=1,N1F
            N1=NW(K,I,M)
            IA0=IEVV(K,I,M)
            IA1=0
            IF (N1.EQ.J) THEN
              NFOUND=1
              IF (DFLOAT(IA0*IA0)/DFLOAT(IAM*IAM).LT.CRIT1) GO TO 652
              IF (JS.EQ.2) THEN
                IA1=-IA11     ! CHECK FOR INVERTED SIGN
              ELSE
                IA1=IA11
              ENDIF
            ENDIF
            IA10=IA1-IA0
            IDD=IDD+IA10*IA10
          ENDDO
          IF (NFOUND.EQ.0) GO TO 652   !  IDD=IDD+IA11*IA11
          IF (IDD.LT.IDM(I)) THEN
            JMIN=J
            IDM(I)=IDD
C
C           IDD is a square of (distance*10000).
C           Distance of 0.76 (IDD=58000000) is well enough to identify
C           the eigenvectors. It corresponds to value of projection
C           (scalar mult.) of 0.71, and to 50% composition of the component.
C           Projection = 1.0 - (distance**2)/2
C
            IF (IDM(I).LT.58000000) GO TO 653
          ENDIF
  652     CONTINUE
        ENDDO
      ENDDO
      IF (JMIN.EQ.0) IFAIL=1    ! COMPARISON LIST IS EMPTY FOR I
  653 IF (IFAIL.EQ.0) THEN
        KW1JM=KW1(JMIN)
        JK=IABS(KW1JM)
        IF (JK.NE.0) THEN
C           COMPARE DISTANCES IDM(I) AND IDM(JK)
          ID1=IDM(JK)
          IF (IDM(I).GT.ID1) THEN
C               EXCLUDE JMIN FROM COMPARE LIST FOR CURRENT I
C               AND REDO THE CYCLE FOR CURRENT I
            IEXR=I
          ELSE
C               EXCLUDE JMIN FROM COMPARE LIST FOR I=JK
C               AND RETURN BACK TO I=JK.
C               NOTE THAT JK CAN NOT BE GREATER THAN I
            DO J=1,IM
              IF (KW1(J).GE.JK) KW1(J)=0
            ENDDO
            I=JK
            IEXR=JK
          ENDIF
          IXB=(IEXR-1)*IM+JMIN-1
          IXCL(IXB/32+1)=IBSET(IXCL(IXB/32+1),IAND(IXB,31))
          GO TO 642
        ENDIF
      ENDIF
 654  IF (IFAIL.NE.0) THEN
C
C             IDENTIFICATION FAILED
C
          WRITE (IW,670) FJ(M)
 670      FORMAT
     1      (34H Failed to label eigenvectors for ,A8)
          RETURN
      ENDIF
C       IF WE GET HERE, IDENTIFICATION I -> JMIN IS FOUND
      IF (KW1(JMIN).GE.0) KW1(JMIN)=I
      GO TO 641
 660  CONTINUE
      DO J=1,IM            ! LOOP THRU TERMS
        K=KW1(J)           ! K IS INDEX OF EIGENVECTOR LABELED BY TERM J
        DO N1=1,N1F    ! FIND THE JTH TERM AMONG ALL COMPONENTS
          NN1=NW(N1,K,M)
          IF (NN1.EQ.0) GO TO 680
          IF (NN1.EQ.J.AND.N1.NE.1) THEN
            IA1=IEVV(N1,K,M)
            NK1=NW(N1,K,M)
            DO N2=N1-1,1,-1
              IEVV(N2+1,K,M)=IEVV(N2,K,M)
              NW(N2+1,K,M)=NW(N2,K,M)
            ENDDO
            IEVV(1,K,M)=IA1
            NW(1,K,M)=NK1
            GO TO 680
          ENDIF
        ENDDO
 680    CONTINUE
      ENDDO
      RETURN
      END

      SUBROUTINE ORDOUT2
C
      IMPLICIT REAL*8 (A-H,O-Z)
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      REAL*4 XMIN,XMAX,XMAX0,DELXC,SDX,DENOM
      INTEGER*2 LF,LF0,LDEG
      COMMON/LBLOC/X(IB),X0(IB),D_PTR,V_PTR,B(JB,JB),
     1  VC(JB),DELXC(IB),XMIN(IB),XMAX(IB),SDX(IB),DENOM(IB),
     2  LF(IB),LF0(IB),XMAX0(IB),LDEG(IB)
      DIMENSION DELX(IB),A1(IA,IA),LPG(IB),V(IV,JB),D(JB,1)
      INTEGER*8 D_PTR,V_PTR
      POINTER (D_PTR,D),(D_PTR,DELX),(D_PTR,A1),(V_PTR,V)
C
      CHARACTER NORUN*4,ERAS*8
      COMMON/C1/IMAX(IX),ANORM(IX),SEVSQ(IX),NPAR(IC,IC),ERAS(40),GEFF,
     1  GCOE,CRIT,DETB,R,AVDEV,STDEV,PNO,TNO,F,G,DXMAX,TIME0,CRIT1,
     2  ICRIT2,IM,M,K3,NOCYC,NOF,NOCYCE,N,NOCSET,NCDES,NORUN,NMAX,MMAX,
     3  LMAX,LSQM,IMSQ,JLMAX,NYOFS,NSQ,ISORT,IPRNV,IPRNA,IW,IW6,INVECT,
     4  IORD,ISAMPL,idummy,SUMDX
      CHARACTER FJ*8
      INTEGER*2 KW,NF
      COMMON/C2/ TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR,
     1  FJ(IX),EV(IA,3),N3(2),S,EMAX
      INTEGER*8 TW_PTR,T_PTR,EVAL_PTR,W_PTR,KW_PTR,NF_PTR,NC_PTR
      DIMENSION TW(IA,IX),T(IA,IX),EVAL(IA,IX),W(IA,IX),KW(IA,IX),
     1 NF(IA,IX),NC(IA,IX)
      POINTER (TW_PTR,TW),(T_PTR,T),(EVAL_PTR,EVAL),(W_PTR,W),
     1  (KW_PTR,KW),(NF_PTR,NF),(NC_PTR,NC)
c
      CHARACTER ELIDEN*8,CONFIG(IC,2)*8,ERAS1*6,XIDEN(IB)*10,NX(IB)*2,
     1  YIDEN(IA,3,IX)*10
      COMMON/CHAR/ELIDEN,CONFIG,ERAS1,XIDEN,NX,YIDEN
C
c     COMMON/CMLC1/C(IA,IA)
c     COMMON/CMLC2/A(IA,IA)
      COMMON/CMLC/C_PTR,A_PTR
      integer*8 C_PTR,A_PTR
      DIMENSION C(IA,IA),A(IA,IA)
      POINTER (C_PTR,C),(A_PTR,A)
C
C
C
      REAL*4 GFAC,GKK
      COMMON/CMLC5/GKK(IA,IX),GFAC(IA,IX)
C ********************************************************v  A.KRAMIDA
      CHARACTER*12 BUFA,BUFC
      COMMON/BUFS/BUFC_PTR,BUFA_PTR,BUFY_PTR,
c     1  IBC(IX,IB),IBA(IX,1),ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
     1  IBC_PTR,IBA_PTR,ICOFFS,IAOFFS,IA2,MAXB1U,MAXB2U
      INTEGER*8 BUFY_PTR,BUFA_PTR,BUFC_PTR,IBA_PTR,IBC_PTR
      POINTER (BUFY_PTR,BUFY(1)),(BUFA_PTR,BUFA(1)),(BUFC_PTR,BUFC(1))
      POINTER (IBA_PTR,IBA),(IBC_PTR,IBC)
      DIMENSION IBC(IX,IB),IBA(IX,1)
c
C       DETERMINE EIGENVECTORS ORDERING FROM FIRST MAXC COMPONENTS
C       WITH LARGEST AMPLITUDES
      INTEGER*2 NW,NW0
      INTEGER*2 IEVV,IEVV0
      INTEGER*2 II1,II2
      COMMON/VECT/ IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      integer*8 IEVV_PTR,IEVV0_PTR,NW_PTR,NW0_PTR
      DIMENSION IEVV(MAXC,IA,IX),IEVV0(MAXC,IA,IX),NW(MAXC,IA,IX),
     1 NW0(MAXC,IA,IX)
      POINTER(IEVV_PTR,IEVV),(IEVV0_PTR,IEVV0),(NW_PTR,NW),(NW0_PTR,NW0)
c
C ********************************************************^  A.KRAMIDA
C
      REAL*4 GFACT
      INTEGER*2 KE
      DIMENSION GFACT(IA),KE(IA)
      SAVE GFACT,KE
C ***************************************** A.Kramida 2026 v
C       COMPRESSED-COLUMN FORM OF THE NON-ZERO ELEMENTS OF THE G-MATRIX
      PARAMETER (MAXNZ2=1000000)
      INTEGER*4 ICSTA(IA+1),ICROW(MAXNZ2)
      REAL*8 CVAL(MAXNZ2)
C       BOUNDED MIN-HEAP USED WHEN SELECTING THE MAXC LARGEST COMPONENTS
      INTEGER*4 IHEAP(MAXC+1)
      REAL*8 AHEAP(MAXC+1)
      SAVE ICSTA,ICROW,CVAL
!$    INTEGER OMP_GET_MAX_THREADS
C ***************************************** A.Kramida 2026 ^
C
C       REARRANGE EVEC(I,J),IEVV(I,K,J),NW(I,K,J) FOR THE SECOND
C       COUPLING SCHEME, USING TRANSFORMATION MATRIX A
C
C ***************************************** A.Kramida 2026 v
C       THE G-FACTOR IN THE SECOND COUPLING SCHEME IS THE QUADRATIC FORM
C            GKK = SUM(J) SUM(I) A(I,K)*C(I,J)*A(J,K)
C       WHICH WAS SWEPT DENSELY - IM**3 PRODUCTS PER (PARITY,J) BLOCK,
C       WITHOUT EVEN THE ZERO TEST THAT ORDER HAD.  AS IN ORDER, THE
C       NON-ZERO ELEMENTS OF C (WHICH IS STRICTLY DIAGONAL IN EVERY CASE
C       SEEN SO FAR) ARE COLLECTED ONCE PER BLOCK IN COMPRESSED-COLUMN
C       FORM AND ONLY THOSE ARE SWEPT, COSTING IM*(IM+NNZ) INSTEAD.
C       SKIPPING A ZERO ELEMENT REMOVES AN ADDITION OF EXACTLY ZERO, AND
C       THE SURVIVING PRODUCTS KEEP THEIR ORIGINAL ORDER AND ASSOCIATION,
C       SO GKK COMES OUT BIT FOR BIT THE SAME.  SHOULD C EVER BE TOO
C       DENSE TO INDEX WITHIN MAXNZ2, THE ORIGINAL LOOP IS USED INSTEAD.
C       THE GKK SUM AND THE EJK/BUFY SUM ARE INDEPENDENT OF EACH OTHER,
C       SO SPLITTING THE J LOOP THAT USED TO CARRY BOTH CHANGES NOTHING.
      NNZ=0
      DO 4612 J=1,IM
        ICSTA(J)=NNZ+1
        DO 4611 I=1,IM
          GRE=C(I,J)
          IF (GRE.EQ.0.D0) GO TO 4611
          IF (NNZ.GE.MAXNZ2) THEN
            NNZ=-1
            GO TO 4613
          ENDIF
          NNZ=NNZ+1
          ICROW(NNZ)=I
          CVAL(NNZ)=GRE
 4611   CONTINUE
 4612 CONTINUE
      ICSTA(IM+1)=NNZ+1
 4613 CONTINUE
C       THE TRANSFORMATION TO THE SECOND COUPLING SCHEME ITSELF - THE
C       EJK LOOP BELOW - IS WHAT COSTS MOST HERE ONCE THE G-FACTORS ARE
C       SPARSE.  IT IS IM*IM*MAXC PRODUCTS PER (PARITY,J) BLOCK AND, MORE
C       TO THE POINT, IT TOUCHES MAXC SCATTERED ELEMENTS OF EVERY COLUMN
C       OF A FOR EVERY K, SO IT IS BOUND BY MEMORY BANDWIDTH RATHER THAN
C       BY ARITHMETIC.  EACH K IS INDEPENDENT OF EVERY OTHER K - IT
C       WRITES ONLY ITS OWN SLICE BUFY(KOFS+1..KOFS+IM) AND ITS OWN
C       GKK(K,M) - SO THE K LOOP IS SPLIT OVER THREADS.  NOTHING IS
C       REORDERED OR REASSOCIATED WITHIN A K, HENCE THE RESULT IS
C       BIT-IDENTICAL AT ANY THREAD COUNT.
C       MEASURED AT IM=4000: 3.9 S SERIAL -> 0.69 S ON 8 THREADS ->
C       0.59 S ON 32.  IT FLATTENS PAST ~8 THREADS, AS A
C       BANDWIDTH-BOUND KERNEL MUST.  HOISTING THE INVARIANT NW/IEVV
C       LOOKUPS OUT OF THE J LOOP WAS ALSO TRIED AND GAINED NOTHING
C       (3.92 -> 3.82 S), CONFIRMING THAT A, NOT THE INDEX ARRAYS, IS
C       WHAT THE LOOP IS WAITING FOR.
      NF1=MIN0(MAXC,IM)
      LPAR = 0
!$    IF (IM.GT.200 .AND. OMP_GET_MAX_THREADS().GT.1) LPAR = 1
!$OMP PARALLEL DO IF(LPAR.NE.0) DEFAULT(SHARED)
!$OMP&  PRIVATE(K,GKKM,KOFS,J,AJK,IP,I,EJK,I1) SCHEDULE(STATIC)
      DO K=1,IM
        GKKM=0.D0
        KOFS=NYOFS+(K-1)*IM
        IF (NNZ.GE.0) THEN
          DO J=1,IM
            AJK=A(J,K)
            DO IP=ICSTA(J),ICSTA(J+1)-1
              GKKM=GKKM+A(ICROW(IP),K)*CVAL(IP)*AJK
            ENDDO
          ENDDO
        ELSE
          DO J=1,IM
            AJK=A(J,K)
            DO I=1,IM
              GKKM=GKKM+A(I,K)*C(I,J)*AJK
            ENDDO
          ENDDO
        ENDIF
        DO J=1,IM
          EJK=0.D0
C         EJK=EJK+A(I,J)*Y(I,K)
          DO I=1,NF1
            I1=NW(I,K,M)
            IF (I1.GT.0) EJK=EJK+A(I1,J)*DFLOAT(IEVV(I,K,M))*1.D-4
          ENDDO
C         EVEC(J,K)=EJK
          BUFY(J+KOFS)=EJK
        ENDDO
        GKK(K,M)=GKKM
      ENDDO
!$OMP END PARALLEL DO
C ***************************************** A.Kramida 2026 ^
C ***************************************** A.KRAMIDA v
C       MAP EIGENVECTORS INTO IEVV(I,J,M) AND NW(I,J,M), J=1..MAXC.
C       K=NW(I,J,M) KEEPS THE REFERENCE TO EVEC(I,K) FOR CURRENT M:
C          IEVV(I,J,M)=EVEC(I,K);
C       IEVV(I,J,M) IS SORTED BY AMPLITUDE: IEVV(I,J,M)>=IEVV(I,J+1,M);
C       EVEC(I,K) is located inside linear array BUFY using offset pointer.
C
C       A.Kramida 2026: AS IN ORDER, THE N1F LARGEST AMPLITUDES USED TO
C       BE PICKED BY N1F SUCCESSIVE SEARCHES THROUGH ALL IM AMPLITUDES
C       (IM*IM*MAXC OPERATIONS PER (PARITY,J) BLOCK).  SELTOP COLLECTS
C       THEM IN ONE PASS WITH A BOUNDED MIN-HEAP UNDER THE SAME TOTAL
C       ORDER - DECREASING ABSOLUTE AMPLITUDE, TIES IN FAVOUR OF THE
C       LOWER BASIS INDEX, ZERO AMPLITUDES NEVER SELECTED - SO NW AND
C       IEVV ARE UNCHANGED.
      N1F=MIN0(IM,MAXC)
      DO N2=1,IM
        N2OFS=NYOFS+(N2-1)*IM
        CALL SELTOP(BUFY(N2OFS+1),IM,N1F,NH,IHEAP,AHEAP)
        DO K=1,NH
          K1=IHEAP(K)
          NW(K,N2,M)=K1       ! -128 for 1-byte type
C         IEVV(K,N2,M)=100.D0*EVEC(K1,N2)+0.5D0
C ***************************************** A.Kramida 2026 v
C           ROUND TO NEAREST BY MAGNITUDE - SEE THE SAME PLACE IN ORDER.
          EVAMP=1.D4*BUFY(K1+N2OFS)
          IEVV(K,N2,M)=EVAMP+DSIGN(0.5D0,EVAMP)
C ***************************************** A.Kramida 2026 ^
        ENDDO
        DO I=NH+1,MAXC
          IEVV(I,N2,M)=0
          NW(I,N2,M)=0 ! -128 FOR INTEGER*1 NW
        ENDDO
      ENDDO
  680 RETURN
      END
C
C ***************************************** A.Kramida 2026 v
      SUBROUTINE SELTOP(V,NV,NTOP,NH,IHEAP,AHEAP)
C
C          SELECT THE NTOP LARGEST OF ABS(V(1..NV)) AND RETURN THEIR
C          INDICES IN IHEAP(1..NH), ORDERED BY DECREASING ABS(V) WITH
C          TIES BROKEN IN FAVOUR OF THE LOWER INDEX.  ELEMENTS EQUAL TO
C          ZERO ARE NEVER SELECTED, SO NH CAN COME BACK LESS THAN NTOP.
C          IHEAP AND AHEAP ARE SCRATCH OF LENGTH NTOP+1 AT LEAST.
C
C          THIS IS THE SAME TOTAL ORDER THAT NTOP SUCCESSIVE 'FIND THE
C          LARGEST NOT YET TAKEN' SEARCHES PRODUCE, BUT IT COSTS
C          NV+NTOP*LOG(NTOP) OPERATIONS INSTEAD OF NV*NTOP.  A BOUNDED
C          MIN-HEAP HOLDS THE BEST NTOP SEEN SO FAR; ITS ROOT IS THE
C          WEAKEST OF THEM, HENCE THE FIRST TO BE DISPLACED.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(*),IHEAP(*),AHEAP(*)
C
      NH=0
      DO N1=1,NV
        EV1=DABS(V(N1))
        IF (EV1.LE.0.D0) GO TO 30
        IF (NH.GE.NTOP) THEN
C            KEEP THE COMPONENT ONLY IF IT BEATS THE WEAKEST ONE HELD
          IF (EV1.LT.AHEAP(1)) GO TO 30
          IF (EV1.EQ.AHEAP(1).AND.N1.GT.IHEAP(1)) GO TO 30
          AHEAP(1)=EV1
          IHEAP(1)=N1
          CALL SELSFT(AHEAP,IHEAP,1,NH)
          GO TO 30
        ENDIF
C          HEAP NOT FULL YET - APPEND AND SIFT UP
        NH=NH+1
        AHEAP(NH)=EV1
        IHEAP(NH)=N1
        NHC=NH
   10   IF (NHC.GT.1) THEN
          NHP=NHC/2
          IF (ISWKR(AHEAP,IHEAP,NHC,NHP).NE.0) THEN
            CALL SELSWP(AHEAP,IHEAP,NHC,NHP)
            NHC=NHP
            GO TO 10
          ENDIF
        ENDIF
   30   CONTINUE
      ENDDO
C          HEAP-SORT: MOVE THE ROOT (THE WEAKEST ENTRY) TO THE BACK
C          REPEATEDLY, WHICH LEAVES IHEAP(1..NH) IN DECREASING ORDER
      DO NHS=NH,2,-1
        CALL SELSWP(AHEAP,IHEAP,1,NHS)
        CALL SELSFT(AHEAP,IHEAP,1,NHS-1)
      ENDDO
      RETURN
      END
C
      FUNCTION ISWKR(AHEAP,IHEAP,I1,I2)
C
C          1 IF HEAP ENTRY I1 IS WEAKER THAN ENTRY I2, AND SO BELONGS
C          CLOSER TO THE ROOT: SMALLER AMPLITUDE, OR AN EQUAL AMPLITUDE
C          WITH A LARGER BASIS INDEX.  USED BY SELTOP AND SELSFT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AHEAP(*),IHEAP(*)
C
      ISWKR=0
      IF (AHEAP(I1).LT.AHEAP(I2)) THEN
        ISWKR=1
      ELSE IF (AHEAP(I1).EQ.AHEAP(I2).AND.IHEAP(I1).GT.IHEAP(I2)) THEN
        ISWKR=1
      ENDIF
      RETURN
      END
C
      SUBROUTINE SELSWP(AHEAP,IHEAP,I1,I2)
C
C          EXCHANGE HEAP ENTRIES I1 AND I2.  USED BY SELTOP AND SELSFT.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AHEAP(*),IHEAP(*)
C
      EVT=AHEAP(I1)
      AHEAP(I1)=AHEAP(I2)
      AHEAP(I2)=EVT
      IHT=IHEAP(I1)
      IHEAP(I1)=IHEAP(I2)
      IHEAP(I2)=IHT
      RETURN
      END
C
      SUBROUTINE SELSFT(AHEAP,IHEAP,LROOT,LEND)
C
C          SIFT ENTRY LROOT DOWN THROUGH THE HEAP 1..LEND.  USED BY
C          SELTOP.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION AHEAP(*),IHEAP(*)
C
      NHP=LROOT
   10 NHC=2*NHP
      IF (NHC.GT.LEND) RETURN
      IF (NHC.LT.LEND) THEN
        IF (ISWKR(AHEAP,IHEAP,NHC+1,NHC).NE.0) NHC=NHC+1
      ENDIF
      IF (ISWKR(AHEAP,IHEAP,NHC,NHP).EQ.0) RETURN
      CALL SELSWP(AHEAP,IHEAP,NHC,NHP)
      NHP=NHC
      GO TO 10
      END
C ***************************************** A.Kramida 2026 ^
C
C
C     SUBROUTINE WR(IFS,IM,A,IA,B,IMSQ)
C
C     IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSION A(IA,IA),B(IMSQ)
C
C     IM2=IM*IM
C     K=0
C     DO 100 J=1,IM
C       I1=J
C       IF (IMSQ.GE.IM2) I1=1
C       DO 100 I=I1,IM
C         K=K+1
C 100 B(K)=A(I,J)
C     WRITE (IFS) (B(I),I=1,IMSQ)
C     RETURN
C     END
C
C     SUBROUTINE RD(IFS,IM,A,IA,B,IMSQ)
C
C     IMPLICIT REAL*8 (A-H,O-Z)
C     DIMENSION A(IA,IA),B(IMSQ)
C
C     READ (IFS) (B(I),I=1,IMSQ)
C     IM2=IM*IM
C     K=0
C     DO 100 J=1,IM
C       I1=J
C       IF (IMSQ.GE.IM2) I1=1
C       DO 100 I=I1,IM
C         K=K+1
C 100 A(I,J)=B(K)
C     IF (IMSQ.LT.IM2) THEN
C       DO 200 J=1,IM
C       DO 200 I=J,IM
C 200     A(J,I)=A(I,J)
C     ENDIF
C 300 RETURN
C     END
C
      SUBROUTINE WR1(IND1,IND2,IM,A,IA,B_PTR,IBIND,IMSQ,
     *               IOFFS,MAXLEN,MAXI1,MAXI2)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 B,BB
      INTEGER*8 B_PTR
      POINTER (B_PTR,B(1))
      INTEGER*2 BB1,I,J,I1
      REAL*8 R8
      DIMENSION A(IA,IA),BB1(6),IBIND(MAXI1,MAXI2)
      EQUIVALENCE (R8,BB1(1)),(BB,BB1(1)),(I,BB1(5)),(J,BB1(6))
C
      IM2=IM*IM
      IBIND(IND1,IND2)=IOFFS+1
      IFILL=0
      DO J=1,IM
        I1=J
        IF (IMSQ.GE.IM2) I1=1
        DO I=I1,IM
          R8=A(I,J)
          IF (DABS(R8).GT.1.D-30) THEN
            IF (IOFFS.GE.MAXLEN) THEN
c              WRITE (IW,10) MAXLEN
c   10         FORMAT(11H More than ,I9,28H non-zero elements in matrix)
c              STOP '(Dimension error)'
              NEWSIZ=MAXLEN*2
              CALL RALLOC(B_PTR,MAXLEN,NEWSIZ,12)
              MAXLEN=NEWSIZ
            ENDIF
            IOFFS=IOFFS+1
            B(IOFFS)=BB
            IFILL=1
          ENDIF
        ENDDO
      ENDDO
      IF (IFILL.EQ.0) THEN
        IOFFS=IOFFS+1
        IF (IOFFS.GT.MAXLEN) THEN
          NEWSIZ=MAXLEN*2
          CALL RALLOC(B_PTR,MAXLEN,NEWSIZ,12)
          MAXLEN=NEWSIZ
        ENDIF
        I=1
        J=1
        R8=0.D0
        B(IOFFS)=BB
      ENDIF
      RETURN
      END
C
      SUBROUTINE RD1(IND1,IND2,IM,A,IA,B,IBIND,IMSQ,
     *               MOFFS,MAXLEN,MAXI1,MAXI2)
C
C     PARAMETER (IB=200)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*12 B,BB
      INTEGER*2 BB1,I,J
      REAL*8 R8
      DIMENSION A(IA,IA),B(MAXLEN),BB1(6),IBIND(MAXI1,MAXI2)
      EQUIVALENCE (R8,BB1(1)),(BB,BB1(1)),(I,BB1(5)),(J,BB1(6))
C
      IOFFS1=IBIND(IND1,IND2)
      IF (IND2.LT.MAXI2) THEN
        IOFFS2=IBIND(IND1,IND2+1)
      ELSE
        IF (IND1.LT.MAXI1) THEN
          IOFFS2=IBIND(IND1+1,1)
        ELSE
          IOFFS2=MOFFS+1
        ENDIF
      ENDIF
      DO J=1,IM
        DO I=1,IM
          A(I,J)=0.D0
        ENDDO
      ENDDO
C
      IM2=IM*IM
      IF (IMSQ.LT.IM2) THEN
        DO K=IOFFS1,IOFFS2-1
          BB=B(K)
          IF (I.EQ.0) RETURN
          A(I,J)=R8
          A(J,I)=R8
        ENDDO
      ELSE
        DO K=IOFFS1,IOFFS2-1
          BB=B(K)
          IF (I.EQ.0) RETURN
          A(I,J)=R8
        ENDDO
      ENDIF
  300 RETURN
      END
C
C
c      SUBROUTINE LSS(N,M,I,A,B,ERAS,DET)
      SUBROUTINE LSS(N,M,I,A,B,DET)
C
      IMPLICIT REAL*8 (A-H,O-Z)
c      CHARACTER*8 ERAS
      DIMENSION A(I,N), B(I,M)    !,D(N)
C
      NN=N
      Nsing = 0
      Jbreak = NN
c     if (NN.EQ.43) write (6,991) A(32,32)
c 991 format (' A(32,32)=',F20.3)
      IF (NN.LE.0.OR.NN.GT.I) THEN
        WRITE (6,99) N
   99   FORMAT ('0LSS CALLED WITH N.LE.0 OR N.GT.I.  N= ',I4)
        DET=0.0D0
      ELSE
        MM=M
        SN=1.0D0
        DO J=1,NN
          L=J-1
          IF (J.NE.NN) THEN
            T=DABS(A(J,J))
            M1=J
            M2=J+1
            DO K=M2,NN
              X=DABS(A(K,J))
              IF (X.GT.T) THEN
                T=X
                M1=K
              ENDIF
            ENDDO
c            if (NN.EQ.43.and.J.eq.32) write (6,992) T,X,M1
c992   format (' OLD MAX=',F20.3,'; NEW MAX=',F20.3,'; M1=',I4)
            IF (M1.NE.J) THEN
              DO K=1,NN
                T=A(J,K)
                A(J,K)=A(M1,K)
                A(M1,K)=T
              ENDDO
              SN=-SN
              IF (MM.GT.0) THEN
                DO K=1,MM
                  T=B(J,K)
                  B(J,K)=B(M1,K)
                  B(M1,K)=T
                ENDDO
              ENDIF
            ENDIF
4           IF (A(J,J).EQ.0.D0) then
              Nsing = 4
              Jbreak = J
              GO TO 13
            endif
            AJJ1=1.D0/A(J,J)
            DO K=M2,NN
              S1=0.0D0
              S2=0.0D0
              IF (L.NE.0) S1=DOTPRO(L,A(J,1),I,A(1,K),1)
              A(J,K)=(A(J,K)-S1)*AJJ1
              S2=DOTPRO(J,A(K,1),I,A(1,M2),1)
              A(K,M2)=A(K,M2)-S2
            ENDDO
          ENDIF
          IF (MM.GT.0) THEN
c           IF (A(J,J).EQ.0.D0) GO TO 13
7           IF (A(J,J).EQ.0.D0) then
              Nsing = 7
              Jbreak = J
              GO TO 13
            endif
            AJJ1=1.D0/A(J,J)
            DO K=1,MM
              S1=0.0D0
              IF (L.NE.0) S1=DOTPRO(L,A(J,1),I,B(1,K),1)
              B(J,K)=(B(J,K)-S1)*AJJ1
            ENDDO
          ENDIF
        ENDDO
        DET=A(1,1)*SN
c       IF (DET.EQ.0.D0) GO TO 13
91      IF (DET.EQ.0.D0) then
          Nsing = 91
          GO TO 13
        endif
        IF (N.EQ.1) GO TO 15
        DO J=2,NN
          DET=DET*A(J,J)
        ENDDO
c       IF (DET.EQ.0.D0) GO TO 13
101     IF (DET.EQ.0.D0) then
          Nsing = 101
          GO TO 13
        endif
        IF (MM.EQ.0) GO TO 15
        M3=NN-1
        DO J=1,MM
          DO L=1,M3
            M1=NN-L
            S1=0.0D0
            M2=M1+1
            K=NN-M2+1
            S1=DOTPRO(K,A(M1,M2),I,B(M2,J),1)
            B(M1,J)=B(M1,J)-S1
          ENDDO
        ENDDO
        GO TO 15
   13   WRITE (6,98) N,Nsing,Jbreak
   98   FORMAT ('0LSS SINGULAR SYSTEM. INPUT DESTROYED.  N=',I4
     1    ,', AT LABEL',I4,', J=',I4 )
        DET=0.0D0
      ENDIF
15    RETURN
      END
C
      FUNCTION DOTPRO(N,X,IX,Y,IY)
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION X(1),Y(1)
C
      DOTPRO=0.0D0
      IF (N.LE.0) RETURN
      SUM=0.0D0
      J=1-IX
      K=1-IY
      DO I=1,N
        J=J+IX
        K=K+IY
        SUM=SUM+X(J)*Y(K)
      ENDDO
      DOTPRO=SUM
      RETURN
      END
C
      SUBROUTINE TRED2(NM,N,A,D,E,Z)                                    TRED2  3
C     DIMENSION A(NM,N),D(N),E(N),Z(NM,N)                               TRED2 56
C***BEGIN PROLOGUE  TRED2                                               TRED2  4
C***DATE WRITTEN   760101   (YYMMDD)                                    TRED2  5
C***REVISION DATE  830518   (YYMMDD)                                    TRED2  6
C***CATEGORY NO.  D4C1B1                                                TRED2  7
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK                          TRED2  8
C***AUTHOR  SMITH, B. T., ET AL.                                        TRED2  9
C***PURPOSE  REDUCE REAL SYMMETRIC MATRIX TO SYMMETRIC TRIDIAGONAL      TRED2 10
C            MATRIX USING AND ACCUMULATING ORTHOGONAL TRANSFORMATION    TRED2 11
C***DESCRIPTION                                                         TRED2 12
C                                                                       TRED2 13
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TRED2,    TRED2 14
C     NUM. MATH. 11, 181-195(1968) BY MARTIN, REINSCH, AND WILKINSON.   TRED2 15
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).   TRED2 16
C                                                                       TRED2 17
C     THIS SUBROUTINE REDUCES A REAL SYMMETRIC MATRIX TO A              TRED2 18
C     SYMMETRIC TRIDIAGONAL MATRIX USING AND ACCUMULATING               TRED2 19
C     ORTHOGONAL SIMILARITY TRANSFORMATIONS.                            TRED2 20
C                                                                       TRED2 21
C     ON INPUT                                                          TRED2 22
C                                                                       TRED2 23
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TRED2 24
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TRED2 25
C          DIMENSION STATEMENT.                                         TRED2 26
C                                                                       TRED2 27
C        N IS THE ORDER OF THE MATRIX.                                  TRED2 28
C                                                                       TRED2 29
C        A CONTAINS THE REAL SYMMETRIC INPUT MATRIX.  ONLY THE          TRED2 30
C          LOWER TRIANGLE OF THE MATRIX NEED BE SUPPLIED.               TRED2 31
C                                                                       TRED2 32
C     ON OUTPUT                                                         TRED2 33
C                                                                       TRED2 34
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE TRIDIAGONAL MATRIX.    TRED2 35
C                                                                       TRED2 36
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE TRIDIAGONAL         TRED2 37
C          MATRIX IN ITS LAST N-1 POSITIONS.  E(1) IS SET TO ZERO.      TRED2 38
C                                                                       TRED2 39
C        Z CONTAINS THE ORTHOGONAL TRANSFORMATION MATRIX                TRED2 40
C          PRODUCED IN THE REDUCTION.                                   TRED2 41
C                                                                       TRED2 42
C        A AND Z MAY COINCIDE.  IF DISTINCT, A IS UNALTERED.            TRED2 43
C                                                                       TRED2 44
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        TRED2 45
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         TRED2 46
C     ------------------------------------------------------------------TRED2 47
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, TRED2 48
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    TRED2 49
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    TRED2 50
C                 1976.                                                 TRED2 51
C***ROUTINES CALLED  (NONE)                                             TRED2 52
C***END PROLOGUE  TRED2                                                 TRED2 53
C                                                                       TRED2 54
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,N,II,NM,JP1                                       TRED2 55
      DIMENSION A(NM,N),D(N),E(N),Z(N,N)                                TRED2 56
C     REAL A(NM,N),D(N),E(N),Z(NM,N)                                    TRED2 56
C     REAL F,G,H,HH,SCALE                                               TRED2 57
C                                                                       TRED2 58
C***FIRST EXECUTABLE STATEMENT  TRED2                                   TRED2 59
      DO I = 1, N                                                       TRED2 60
C                                                                       TRED2 61
         DO J = 1, I                                                    TRED2 62
            Z(I,J) = A(I,J)                                             TRED2 63
         ENDDO
      ENDDO                                                             TRED2 64
C                                                                       TRED2 65
      IF (N .EQ. 1) GO TO 320                                           TRED2 66
C     .......... FOR I=N STEP -1 UNTIL 2 DO -- ..........               TRED2 67
      DO II = 2, N                                                      TRED2 68
         I = N + 2 - II                                                 TRED2 69
         L = I - 1                                                      TRED2 70
         H = 0.D0                                                       TRED2 71
         SCALE = 0.D0                                                   TRED2 72
         IF (L .LT. 2) GO TO 130                                        TRED2 73
C     .......... SCALE ROW (ALGOL TOL THEN NOT NEEDED) ..........       TRED2 74
         DO K = 1, L                                                    TRED2 75
           SCALE = SCALE + DABS(Z(I,K))                                 TRED2 76
         ENDDO
C                                                                       TRED2 77
         IF (SCALE .NE. 0.D0) GO TO 140                                 TRED2 78
  130    E(I) = Z(I,L)                                                  TRED2 79
         GO TO 290                                                      TRED2 80
C                                                                       TRED2 81
  140    DO K = 1, L                                                    TRED2 82
            Z(I,K) = Z(I,K) / SCALE                                     TRED2 83
            H = H + Z(I,K) * Z(I,K)                                     TRED2 84
         ENDDO                                                          TRED2 85
C                                                                       TRED2 86
         F = Z(I,L)                                                     TRED2 87
         G = -SIGN(DSQRT(H),F)                                          TRED2 88
         E(I) = SCALE * G                                               TRED2 89
         H = H - F * G                                                  TRED2 90
         Z(I,L) = F - G                                                 TRED2 91
         F = 0.D0                                                       TRED2 92
C                                                                       TRED2 93
C ***************************************** A.Kramida 2026 v
C       THE TWO O(L**2) LOOPS BELOW CARRY THE BULK OF THE WORK OF THE
C       HOUSEHOLDER REDUCTION AND ARE NOW RUN IN PARALLEL OVER J.  EACH
C       J STILL ACCUMULATES ITS OWN DOT PRODUCT SERIALLY, AND EVERY
C       ELEMENT OF Z RECEIVES THE SAME OPERATIONS IN THE SAME ORDER AS
C       BEFORE, SO THE RESULT IS BIT FOR BIT UNCHANGED.  TWO SERIAL
C       DEPENDENCIES HAD TO BE LIFTED OUT FIRST:
C         - THE SUM F, NOW ACCUMULATED IN A LOOP OF ITS OWN IN THE
C           ORIGINAL ORDER OF J;
C         - THE UPDATE OF E(J), WHICH THE SECOND LOOP READ FOR K.LE.J
C           ONLY AFTER IT HAD BEEN OVERWRITTEN, AND WHICH IS NOW
C           APPLIED TO THE WHOLE OF E BEFOREHAND.
C       WITHOUT -fopenmp THE DIRECTIVES ARE COMMENTS AND THIS RUNS
C       SERIALLY.  BLOCKS BELOW L=1000 STAY SERIAL: MEASURED CROSSOVER
C       IS NEAR L=1100 (BELOW IT THE THREAD OVERHEAD PER HOUSEHOLDER
C       STEP EXCEEDS THE GAIN), ABOVE IT THE SPEEDUP REACHES 5.4x ON 32
C       THREADS AT L=4000.
!$OMP PARALLEL DO IF(L.GT.1000) DEFAULT(SHARED) PRIVATE(J,K,G,JP1)
!$OMP&  SCHEDULE(STATIC,1)
         DO J = 1, L                                                    TRED2 94
            Z(J,I) = Z(I,J) / H                                         TRED2 95
            G = 0.D0                                                    TRED2 96
C     .......... FORM ELEMENT OF A*U ..........                         TRED2 97
            DO K = 1, J                                                 TRED2 98
              G = G + Z(J,K) * Z(I,K)                                   TRED2 99
            ENDDO
C                                                                       TRED2100
            JP1 = J + 1                                                 TRED2101
            IF (L .GE. JP1) THEN                                        TRED2102
C                                                                       TRED2103
              DO K = JP1, L                                             TRED2104
                G = G + Z(K,J) * Z(I,K)                                 TRED2105
              ENDDO
            ENDIF
C     .......... FORM ELEMENT OF P ..........                           TRED2106
            E(J) = G / H                                                TRED2107
         ENDDO                                                          TRED2109
!$OMP END PARALLEL DO
         DO J = 1, L
            F = F + E(J) * Z(I,J)                                       TRED2108
         ENDDO
C                                                                       TRED2110
         HH = F / (H + H)                                               TRED2111
         DO J = 1, L
            E(J) = E(J) - HH * Z(I,J)                                   TRED2115
         ENDDO
C     .......... FORM REDUCED A ..........                              TRED2112
!$OMP PARALLEL DO IF(L.GT.1000) DEFAULT(SHARED) PRIVATE(J,K,F,G)
!$OMP&  SCHEDULE(STATIC,1)
         DO J = 1, L                                                    TRED2113
            F = Z(I,J)                                                  TRED2114
            G = E(J)                                                    TRED2116
C                                                                       TRED2117
            DO K = 1, J                                                 TRED2118
              Z(J,K) = Z(J,K) - F * E(K) - G * Z(I,K)                   TRED2119
            ENDDO
         ENDDO                                                          TRED2120
C                                                                       TRED2121
  290    D(I) = H                                                       TRED2122
      ENDDO                                                             TRED2123
C                                                                       TRED2124
  320 D(1) = 0.D0                                                       TRED2125
      E(1) = 0.D0                                                       TRED2126
C     .......... ACCUMULATION OF TRANSFORMATION MATRICES ..........     TRED2127
      DO I = 1, N                                                       TRED2128
         L = I - 1                                                      TRED2129
         IF (D(I) .EQ. 0.D0) GO TO 380                                  TRED2130
C                                                                       TRED2131
C       BACK-TRANSFORMATION: EACH J WORKS ON COLUMN J ALONE AND ONLY
C       READS COLUMN I AND ROW I, NEITHER OF WHICH IT WRITES, SO THE
C       ITERATIONS ARE INDEPENDENT AND THE RESULT IS UNCHANGED.
!$OMP PARALLEL DO IF(L.GT.1000) DEFAULT(SHARED) PRIVATE(J,K,G)
!$OMP&  SCHEDULE(STATIC)
         DO J = 1, L                                                    TRED2132
            G = 0.D0                                                    TRED2133
C                                                                       TRED2134
            DO K = 1, L                                                 TRED2135
              G = G + Z(I,K) * Z(K,J)                                   TRED2136
            ENDDO
C                                                                       TRED2137
            DO K = 1, L                                                 TRED2138
               Z(K,J) = Z(K,J) - G * Z(K,I)                             TRED2139
            ENDDO
         ENDDO                                                          TRED2140
!$OMP END PARALLEL DO
C                                                                       TRED2141
  380    D(I) = Z(I,I)                                                  TRED2142
         Z(I,I) = 1.D0                                                  TRED2143
         IF (L .GE. 1) THEN                                             TRED2144
C                                                                       TRED2145
           DO J = 1, L                                                  TRED2146
             Z(I,J) = 0.D0                                              TRED2147
             Z(J,I) = 0.D0                                              TRED2148
           ENDDO                                                        TRED2149
         ENDIF
C                                                                       TRED2150
      ENDDO                                                             TRED2151
C                                                                       TRED2152
      RETURN                                                            TRED2153
      END                                                               TRED2154
C
      SUBROUTINE TQL2(NM,N,D,E,Z,IERR)                                  TQL2   3
C***BEGIN PROLOGUE  TQL2                                                TQL2   4
C***DATE WRITTEN   760101   (YYMMDD)                                    TQL2   5
C***REVISION DATE  830518   (YYMMDD)                                    TQL2   6
C***CATEGORY NO.  D4A5,D4C2A                                            TQL2   7
C***KEYWORDS  EIGENVALUES,EIGENVECTORS,EISPACK                          TQL2   8
C***AUTHOR  SMITH, B. T., ET AL.                                        TQL2   9
C***PURPOSE  COMPUTE EIGENVALUES AND EIGENVECTORS OF SYMMETRIC          TQL2  10
C            TRIDIAGONAL MATRIX.                                        TQL2  11
C***DESCRIPTION                                                         TQL2  12
C                                                                       TQL2  13
C     THIS SUBROUTINE IS A TRANSLATION OF THE ALGOL PROCEDURE TQL2,     TQL2  14
C     NUM. MATH. 11, 293-306(1968) BY BOWDLER, MARTIN, REINSCH, AND     TQL2  15
C     WILKINSON.                                                        TQL2  16
C     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).   TQL2  17
C                                                                       TQL2  18
C     THIS SUBROUTINE FINDS THE EIGENVALUES AND EIGENVECTORS            TQL2  19
C     OF A SYMMETRIC TRIDIAGONAL MATRIX BY THE QL METHOD.               TQL2  20
C     THE EIGENVECTORS OF A FULL SYMMETRIC MATRIX CAN ALSO              TQL2  21
C     BE FOUND IF  TRED2  HAS BEEN USED TO REDUCE THIS                  TQL2  22
C     FULL MATRIX TO TRIDIAGONAL FORM.                                  TQL2  23
C                                                                       TQL2  24
C     ON INPUT                                                          TQL2  25
C                                                                       TQL2  26
C        NM MUST BE SET TO THE ROW DIMENSION OF TWO-DIMENSIONAL         TQL2  27
C          ARRAY PARAMETERS AS DECLARED IN THE CALLING PROGRAM          TQL2  28
C          DIMENSION STATEMENT.                                         TQL2  29
C                                                                       TQL2  30
C        N IS THE ORDER OF THE MATRIX.                                  TQL2  31
C                                                                       TQL2  32
C        D CONTAINS THE DIAGONAL ELEMENTS OF THE INPUT MATRIX.          TQL2  33
C                                                                       TQL2  34
C        E CONTAINS THE SUBDIAGONAL ELEMENTS OF THE INPUT MATRIX        TQL2  35
C          IN ITS LAST N-1 POSITIONS.  E(1) IS ARBITRARY.               TQL2  36
C                                                                       TQL2  37
C        Z CONTAINS THE TRANSFORMATION MATRIX PRODUCED IN THE           TQL2  38
C          REDUCTION BY  TRED2, IF PERFORMED.  IF THE EIGENVECTORS      TQL2  39
C          OF THE TRIDIAGONAL MATRIX ARE DESIRED, Z MUST CONTAIN        TQL2  40
C          THE IDENTITY MATRIX.                                         TQL2  41
C                                                                       TQL2  42
C      ON OUTPUT                                                        TQL2  43
C                                                                       TQL2  44
C        D CONTAINS THE EIGENVALUES IN ASCENDING ORDER.  IF AN          TQL2  45
C          ERROR EXIT IS MADE, THE EIGENVALUES ARE CORRECT BUT          TQL2  46
C          UNORDERED FOR INDICES 1,2,...,IERR-1.                        TQL2  47
C                                                                       TQL2  48
C        E HAS BEEN DESTROYED.                                          TQL2  49
C                                                                       TQL2  50
C        Z CONTAINS ORTHONORMAL EIGENVECTORS OF THE SYMMETRIC           TQL2  51
C          TRIDIAGONAL (OR FULL) MATRIX.  IF AN ERROR EXIT IS MADE,     TQL2  52
C          Z CONTAINS THE EIGENVECTORS ASSOCIATED WITH THE STORED       TQL2  53
C          EIGENVALUES.                                                 TQL2  54
C                                                                       TQL2  55
C        IERR IS SET TO                                                 TQL2  56
C          ZERO       FOR NORMAL RETURN,                                TQL2  57
C          J          IF THE J-TH EIGENVALUE HAS NOT BEEN               TQL2  58
C                     DETERMINED AFTER 30 ITERATIONS.                   TQL2  59
C                                                                       TQL2  60
C     CALLS PYTHAG(A,B) FOR SQRT(A**2 + B**2).                          TQL2  61
C                                                                       TQL2  62
C     QUESTIONS AND COMMENTS SHOULD BE DIRECTED TO B. S. GARBOW,        TQL2  63
C     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY         TQL2  64
C     ------------------------------------------------------------------TQL2  65
C***REFERENCES  B. T. SMITH, J. M. BOYLE, J. J. DONGARRA, B. S. GARBOW, TQL2  66
C                 Y. IKEBE, V. C. KLEMA, C. B. MOLER, *MATRIX EIGEN-    TQL2  67
C                 SYSTEM ROUTINES - EISPACK GUIDE*, SPRINGER-VERLAG,    TQL2  68
C                 1976.                                                 TQL2  69
C***ROUTINES CALLED  PYTHAG                                             TQL2  70
C***END PROLOGUE  TQL2                                                  TQL2  71
C                                                                       TQL2  72
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR                          TQL2  73
      DIMENSION D(N),E(N),Z(N,N)                                        TQL2  74
C     A.Kramida 2026: BUFFERS HOLDING THE GIVENS ROTATIONS OF ONE QL
C       SWEEP SO THAT THEY CAN BE APPLIED TO Z IN ONE BLOCKED, THREADED
C       PASS (TQLROT) INSTEAD OF ONE THIN LOOP PER ROTATION.
C        PARAMETERS IA,IB,JB,IC,IX,IV,MAXB1,MAXB2,MAXBY
C        ARE TAKEN FROM INCLUDE FILE
      INCLUDE 'RCE20KD_PAR.F'
      DIMENSION CSAV(IA),SSAV(IA)
      SAVE CSAV,SSAV
C     REAL D(N),E(N),Z(NM,N)                                            TQL2  74
C     REAL B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2                             TQL2  75
C     REAL PYTHAG                                                       TQL2  76
C                                                                       TQL2  77
C***FIRST EXECUTABLE STATEMENT  TQL2                                    TQL2  78
      IERR = 0                                                          TQL2  79
      IF (N .EQ. 1) GO TO 1001                                          TQL2  80
C     A.Kramida 2026: BATCH THE ROTATIONS UNLESS THE BUFFERS ARE TOO
C       SMALL (CANNOT HAPPEN FOR N SET BY THIS PROGRAM, WHERE N.LE.IA).
      IBAT = 1
      IF (N .GT. IA) IBAT = 0
C                                                                       TQL2  81
      ONE=1.D0
      DO I = 2, N                                                       TQL2  82
        E(I-1) = E(I)                                                   TQL2  83
      ENDDO
C                                                                       TQL2  84
      F = 0.D0                                                          TQL2  85
      B = 0.D0                                                          TQL2  86
      E(N) = 0.D0                                                       TQL2  87
C                                                                       TQL2  88
      DO L = 1, N                                                       TQL2  89
         J = 0                                                          TQL2  90
         H = DABS(D(L)) + DABS(E(L))                                    TQL2  91
         IF (B .LT. H) B = H                                            TQL2  92
C     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........         TQL2  93
         DO M = L, N                                                    TQL2  94
            IF (B + DABS(E(M)) .EQ. B) GO TO 120                        TQL2  95
C     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT               TQL2  96
C                THROUGH THE BOTTOM OF THE LOOP ..........              TQL2  97
         ENDDO                                                          TQL2  98
C                                                                       TQL2  99
  120    IF (M .EQ. L) GO TO 220                                        TQL2 100
  130    IF (J .EQ. 30) GO TO 1000                                      TQL2 101
         J = J + 1                                                      TQL2 102
C     .......... FORM SHIFT ..........                                  TQL2 103
         L1 = L + 1                                                     TQL2 104
         L2 = L1 + 1                                                    TQL2 105
         G = D(L)                                                       TQL2 106
         P = (D(L1) - G) / (2.D0 * E(L))                                TQL2 107
C        R = PYTHAG(P,1.D0)                                             TQL2 108
           PP=DMAX1(DABS(P),ONE)                                        PYTHAG14
           QQ=DMIN1(DABS(P),ONE)                                        PYTHAG15
           IF (QQ.EQ.0.D0) GO TO 20                                     PYTHAG16
   10      CONTINUE                                                     PYTHAG17
           RR=(QQ/PP)**2                                                PYTHAG18
           TT=4.D0+RR                                                   PYTHAG19
           IF (TT.EQ.4.D0) GO TO 20                                     PYTHAG20
           SS=RR/TT                                                     PYTHAG21
           PP=PP+2.D0*PP*SS                                             PYTHAG22
           QQ=QQ*SS                                                     PYTHAG23
           GO TO 10                                                     PYTHAG24
   20    R=PP                                                           PYTHAG25
         D(L) = E(L) / (P + SIGN(R,P))                                  TQL2 109
         D(L1) = E(L) * (P + SIGN(R,P))                                 TQL2 110
         DL1 = D(L1)                                                    TQL2 111
         H = G - D(L)                                                   TQL2 112
         IF (L2 .GT. N) GO TO 145                                       TQL2 113
C                                                                       TQL2 114
         DO I = L2, N                                                   TQL2 115
           D(I) = D(I) - H                                              TQL2 116
         ENDDO
C                                                                       TQL2 117
  145    F = F + H                                                      TQL2 118
C     .......... QL TRANSFORMATION ..........                           TQL2 119
         P = D(M)                                                       TQL2 120
         C = 1.D0                                                       TQL2 121
         C2 = C                                                         TQL2 122
         EL1 = E(L1)                                                    TQL2 123
         S = 0.D0                                                       TQL2 124
         MML = M - L                                                    TQL2 125
C     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........             TQL2 126
         DO II = 1, MML                                                 TQL2 127
            C3 = C2                                                     TQL2 128
            C2 = C                                                      TQL2 129
            S2 = S                                                      TQL2 130
            I = M - II                                                  TQL2 131
            G = C * E(I)                                                TQL2 132
            H = C * P                                                   TQL2 133
            IF (DABS(P) .LT. DABS(E(I))) GO TO 150                      TQL2 134
            C = E(I) / P                                                TQL2 135
            R = DSQRT(C*C+1.D0)                                         TQL2 136
            E(I+1) = S * P * R                                          TQL2 137
            S = C / R                                                   TQL2 138
            C = 1.D0 / R                                                TQL2 139
            GO TO 160                                                   TQL2 140
  150       C = P / E(I)                                                TQL2 141
            R = DSQRT(C*C+1.D0)                                         TQL2 142
            E(I+1) = S * E(I) * R                                       TQL2 143
            S = 1.D0 / R                                                TQL2 144
            C = C * S                                                   TQL2 145
  160       P = C * D(I) - S * G                                        TQL2 146
            D(I+1) = H + S * (C * G + S * D(I))                         TQL2 147
C     .......... FORM VECTOR ..........                                 TQL2 148
C           A.Kramida 2026: RECORD THE ROTATION ONLY; THE WHOLE SWEEP IS
C           APPLIED TO Z AFTER THE LOOP.  Z IS NOT READ BY THE SCALAR
C           RECURRENCE ABOVE, SO DEFERRING ITS UPDATE IS EXACT.
            IF (IBAT .EQ. 0) THEN
              DO K = 1, N                                               TQL2 149
                 H = Z(K,I+1)                                           TQL2 150
                 Z(K,I+1) = S * Z(K,I) + C * H                          TQL2 151
                 Z(K,I) = C * Z(K,I) - S * H                            TQL2 152
              ENDDO                                                     TQL2 153
            ELSE
              CSAV(I) = C
              SSAV(I) = S
            ENDIF
C                                                                       TQL2 154
         ENDDO                                                          TQL2 155
         IF (IBAT .NE. 0) CALL TQLROT(N,Z,CSAV,SSAV,L,M-1)
C                                                                       TQL2 156
         P = -S * S2 * C3 * EL1 * E(L) / DL1                            TQL2 157
         E(L) = S * P                                                   TQL2 158
         D(L) = C * P                                                   TQL2 159
         IF (B + DABS(E(L)) .GT. B) GO TO 130                           TQL2 160
  220    D(L) = D(L) + F                                                TQL2 161
      ENDDO                                                             TQL2 162
C     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........          TQL2 163
      DO II = 2, N                                                      TQL2 164
         I = II - 1                                                     TQL2 165
         K = I                                                          TQL2 166
         P = D(I)                                                       TQL2 167
C                                                                       TQL2 168
         DO J = II, N                                                   TQL2 169
            IF (D(J) .LT. P) THEN                                       TQL2 170
              K = J                                                     TQL2 171
              P = D(J)                                                  TQL2 172
            ENDIF
         ENDDO                                                          TQL2 173
C                                                                       TQL2 174
         IF (K .NE. I) THEN                                             TQL2 175
           D(K) = D(I)                                                  TQL2 176
           D(I) = P                                                     TQL2 177
C                                                                       TQL2 178
           DO J = 1, N                                                  TQL2 179
              P = Z(J,I)                                                TQL2 180
              Z(J,I) = Z(J,K)                                           TQL2 181
              Z(J,K) = P                                                TQL2 182
           ENDDO                                                        TQL2 183
         ENDIF
C                                                                       TQL2 184
      ENDDO                                                             TQL2 185
C                                                                       TQL2 186
      GO TO 1001                                                        TQL2 187
C     .......... SET ERROR -- NO CONVERGENCE TO AN                      TQL2 188
C                EIGENVALUE AFTER 30 ITERATIONS ..........              TQL2 189
 1000 IERR = L                                                          TQL2 190
 1001 RETURN                                                            TQL2 191
      END                                                               TQL2 192

C ***************************************** A.Kramida 2026 v
      SUBROUTINE TQLROT(N,Z,CS,SN,ILO,IHI)
C
C          APPLIES ONE QL SWEEP'S WORTH OF GIVENS ROTATIONS, ACCUMULATED
C          BY TQL2 IN CS/SN, TO THE EIGENVECTOR MATRIX Z.  ROTATION I
C          MIXES COLUMNS I AND I+1 AND THE ROTATIONS MUST BE APPLIED FOR
C          I = IHI DOWN TO ILO, THE ORDER IN WHICH TQL2 PRODUCES THEM.
C
C          EVERY ROW OF Z IS TRANSFORMED INDEPENDENTLY OF EVERY OTHER
C          ROW, SO THE ROW RANGE IS SPLIT OVER THREADS.  WITHIN A ROW
C          THE OPERATIONS AND THEIR ORDER ARE EXACTLY THOSE OF THE
C          ORIGINAL IN-LINE LOOP, HENCE THE RESULT IS BIT-IDENTICAL AT
C          ANY THREAD COUNT.
C
C          THE ROWS ARE PROCESSED IN BLOCKS OF NB: A BLOCK'S SLICE OF
C          THE TWO ACTIVE COLUMNS THEN STAYS IN CACHE WHILE THE WHOLE
C          SWEEP MARCHES ACROSS IT.  BLOCKING ALSO SPEEDS UP THE SERIAL
C          PATH, WHICH OTHERWISE STREAMS ALL OF Z ONCE PER ROTATION.
C
C          THE BEST BLOCK SIZE DEPENDS ON HOW MANY THREADS ARE RUNNING
C          (MEASURED AT N=4000): SERIALLY, KBLKS=64 IS ~10 PCT FASTER
C          THAN NO BLOCKING, WHILE THE NARROW KBLKP=8 IS ALMOST TWICE AS
C          SLOW; WITH THREADS THE ORDER REVERSES, KBLKP=8 BEING NEARLY
C          TWICE AS FAST AS 64 BECAUSE EACH THREAD THEN KEEPS ITS SLICE
C          IN L1 ACROSS THE WHOLE SWEEP.  SO PICK NB AT RUN TIME.
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION Z(N,N),CS(*),SN(*)
      PARAMETER (KBLKS=64,KBLKP=8)
!$    INTEGER OMP_GET_MAX_THREADS
C          THREADING PAYS ONLY FOR BIG BLOCKS - SEE THE HEADER NOTES ON
C          THE TRED2 GUARD; BELOW THIS THE FORK/JOIN COSTS MORE THAN THE
C          SWEEP ITSELF.
      NTH = 1
!$    NTH = OMP_GET_MAX_THREADS()
      LPAR = 0
      IF (NTH.GT.1 .AND. N.GT.1000 .AND. (IHI-ILO+1).GT.64) LPAR = 1
      NB = KBLKS
      IF (LPAR .NE. 0) NB = KBLKP
!$OMP PARALLEL DO IF(LPAR.NE.0)
!$OMP&  DEFAULT(SHARED) PRIVATE(KB,KE,I,K,C,S,H) SCHEDULE(STATIC)
      DO KB = 1, N, NB
        KE = MIN0(KB+NB-1,N)
        DO I = IHI, ILO, -1
          C = CS(I)
          S = SN(I)
          DO K = KB, KE
            H = Z(K,I+1)
            Z(K,I+1) = S * Z(K,I) + C * H
            Z(K,I) = C * Z(K,I) - S * H
          ENDDO
        ENDDO
      ENDDO
!$OMP END PARALLEL DO
      RETURN
      END
C ***************************************** A.Kramida 2026 ^

      FUNCTION IDECCH(I1,NMAX)
      INTEGER*1 I1
      INTEGER*2 IDECCH,I,NMAX
c           A. Kramida: decode alphanumeric codes of config. numbers into integer numbers
      if (I1.LT.0) then
        I = 256+I1
      else
        I = I1
      endif
      IF (I.EQ.42) THEN
        ! Star signifies too many configs:
        ! designation of lower and upper conf will be invalid
        ! in CI param. section RCEINP/RCEOUT
        if (NMAX.LT.188) then
          I = NMAX
        else
          I=188
        endif
      ELSE
        IF (I.LE.254.AND.I.GE.97) THEN
          I=I-61     ! lower-case roman letters and extended printable characters
        ELSE
          IF (I.LE.90.AND.I.GE.65) THEN
            I=I-55   ! upper-case roman letters
          ELSE
            IF (I.LE.57.AND.I.GE.49) I=I-48  ! digits
          ENDIF
        ENDIF
      ENDIF
      IDECCH = I
      RETURN
      END
