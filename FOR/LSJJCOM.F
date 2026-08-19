C     Parameters and COMMON blocks of the LS -> jj recursion.  Split out
C     of LSJJPAR.F so that LSJJRT, which also needs RCGPAR.F, can pull
C     these in without a second IMPLICIT statement.  Include LSJJPAR.F
C     (which adds the IMPLICIT lines) unless RCGPAR.F is already there.
      PARAMETER (
C          LS cfp deck: p,d,f for every n, g,h for n<=2 -> 34 records,
C          119 terms at most, 55318 cfp in all
     +  MXLR=40,
     +  MXTM=120,
     +  MXCF=60000,
C          jj cfp deck: j=1/2..7/2 for every n, 9/2 and 11/2 for n<=2
     +  MXJR=30,
     +  MXJJ=6,
     +  MXJN=12,
     +  MXJL=10,
C          working basis: 327 LS states is the worst case (f^7), where
C          2J <= 31
     +  MXST=340,
     +  MXJI=64,
     +  MXKEY=1200000)
C
C     THE TWO DECKS
C
C       LRL,LRN        l and n of LS record IR
C       LRNT,LRNP      number of terms of l^n, number of cfp columns
C       LRTRM(1..4)    multiplicity, 2L, and the two characters of
C                      Cowan's alpha field as ICHAR values
C       LRPAR(K)       index in the term list of l^(n-1) of column K
C       LRCFP          cfp, row-major, at offset LROF(IR)
C       JNLV(N,IJ)     number of states of j^N, IJ=(2j+1)/2
C       JLJ,JLV        their 2J and seniority
C       JCFP(M,K,N,IJ) jj cfp
C
      REAL*8 LRCFP,JCFP
      COMMON /LSJJD/ LRCFP(MXCF),JCFP(MXJL,MXJL,0:MXJN,MXJJ),
     1  NLREC,NJREC,ILSJRD,
     2  LRL(MXLR),LRN(MXLR),LRNT(MXLR),LRNP(MXLR),LROF(MXLR),
     3  LRTRM(4,MXTM,MXLR),LRPAR(MXTM,MXLR),
     4  JNLV(0:MXJN,MXJJ),JLJ(MXJL,0:MXJN,MXJJ),JLV(MXJL,0:MXJN,MXJJ)
C
C     THE ANSWER  (block /LSJJR/; /LSJJU/ would clash with the routine)
C
C       IUROW(1,I)     term index in LS cfp deck record IUREC
C       IUROW(2,I)     2J
C       IUCOL(1..5,K)  N1, level of j-^N1, N2, level of j+^N2, 2J
C       UMAT(I,K)      < l^n a L S J | ( j-^N1 v1 J1, j+^N2 v2 J2 ) J >
C
C     Rows and columns are grouped by total J: block 2J starts at
C     IUOFR(2J)/IUOFC(2J) and holds NUJR(2J)/NUJC(2J) of them, and the
C     two counts are equal because both are bases of the same space.
C
      COMMON /LSJJR/ UMAT(MXST,MXST),
     1  IULQ,IUNEL,IUJM2,IUJP2,IUREC,IUALT,NUROW,NUCOL,
     2  IUROW(2,MXST),IUCOL(5,MXST),
     3  NUJR(0:MXJI),IUOFR(0:MXJI),NUJC(0:MXJI),IUOFC(0:MXJI)
C
C     SCRATCH OF THE RECURSION.  SLOTS 1 AND 2 HOLD n AND n-1.
C
      COMMON /LSJJS/ TT(MXST,MXST,2),AMAT(MXST,MXST),
     1  BMAT(MXST,MXST),UALT(MXST,MXST),NROW(2),NCOL(2),
     2  IRT(MXST,2),IRJ(MXST,2),ICL(5,MXST,2),
     3  NRJ(0:MXJI,2),IOFR(0:MXJI,2),NCJ(0:MXJI,2),IOFC(0:MXJI,2),
     4  IPOSR(MXTM),KMAP(MXKEY)
