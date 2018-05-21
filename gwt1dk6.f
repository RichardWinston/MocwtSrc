C     Last change: RBW, May 25, 2016
C  DK6DF  READ DK OPTION 
C     ******************************************************************
C
      SUBROUTINE DK6DF(INDK,IDKTIM,IDKRF,
     1              IDKZO,IDKFO,IDKZS,IDKFS,IOUT,IOUTS,DECAY)
C
C     ******************************************************************
C     READ IN DK OPTIONS
C     ******************************************************************
      DOUBLE PRECISION DECAY
      WRITE(IOUTS,1000) INDK
      READ(INDK,*) IDKRF,IDKTIM,IDKFO,IDKFS,IDKZO,IDKZS
      IF(IDKRF.NE.1) THEN
       IDKRF=0
       WRITE(IOUTS,1060)
      ELSE
       WRITE(IOUTS,1061)
      END IF
      IF(IDKTIM.NE.1) THEN
       IDKTIM=0
       WRITE(IOUTS,1010)
      ELSE
       WRITE(IOUTS,1011)
      END IF
      IF(IDKFO.NE.1) THEN
       IDKFO=0
       WRITE(IOUTS,1030)
      ELSE
       WRITE(IOUTS,1031)
       IF(DECAY.NE.0.D0) THEN
        DECAY=0.D0
        WRITE(IOUTS,1032)
       END IF
      END IF
      IF(IDKFS.NE.1) THEN
       IDKFS=0
       WRITE(IOUTS,1050)
      ELSE
       WRITE(IOUTS,1051)
       IF(DECAY.NE.0.D0) THEN
         DECAY=0.D0
         WRITE(IOUTS,1032)
       END IF
      END IF
      IF(IDKZO.NE.1) THEN
       IDKZO=0
       WRITE(IOUTS,1020)
      ELSE
       WRITE(IOUTS,1021)
      END IF
      IF(IDKZS.NE.1) THEN
       IDKZS=0
       WRITE(IOUTS,1040)
      ELSE
       WRITE(IOUTS,1041)
      END IF
      WRITE(IOUTS,*)
      RETURN
1000  FORMAT(/' INPUT FOR DK OPTIONS READ FROM UNIT',I10)
1010  FORMAT(' DECAY AND GROWTH RATES DO NOT CHANGE IN TIME')
1011  FORMAT(' DECAY AND GROWTH RATES CHANGE FOR EACH STRESS PERIOD')
1020  FORMAT(' NO SPATIALLY-VARIABLE ZERO-ORDER GROWTH')
1021  FORMAT(' SPATIALLY-VARIABLE ZERO-ORDER GROWTH')
1030  FORMAT(' NO SPATIALLY-VARIABLE FIRST-ORDER DECAY')
1031  FORMAT(' SPATIALLY-VARIABLE FIRST-ORDER DECAY')
1032  FORMAT('  UNIFORM DECAY FROM MOC INPUT FILE RESET TO ZERO')
1040  FORMAT(' NO DISTINCT SPATIALLY-VARIABLE ZERO-ORDER',
     1' GROWTH FOR SORBED SOLUTE'/
     2'  SORBED MASS GROWS AT SAME RATE AS DISSOLVED')
1041  FORMAT(' DISTINCT SPATIALLY-VARIABLE ZERO-ORDER',
     1' GROWTH FOR SORBED SOLUTE'/
     2'  SORBED MASS GROWS AT DIFFERENT RATE THAN DISSOLVED')
1050  FORMAT(' NO DISTINCT SPATIALLY-VARIABLE FIRST-ORDER',
     1' DECAY FOR SORBED SOLUTE'/
     2'  SORBED MASS DECAYS AT SAME RATE AS DISSOLVED')
1051  FORMAT(' DISTINCT SPATIALLY-VARIABLE FIRST-ORDER',
     1' DECAY FOR SORBED SOLUTE'/
     2'  SORBED MASS DECAYS AT DIFFERENT RATE THAN DISSOLVED')
1060  FORMAT(' NO SPATIALLY-VARIABLE RETARDATION FACTOR')
1061  FORMAT(' SPATIALLY-VARIABLE RETARDATION FACTOR'/
     1'  LAYER-CONSTANT RETARDATION FACTORS FROM MOC INPUT FILE',
     2' WILL NOT BE USED')
      END
C
C DK6AL ALLOCATE SPACE FOR DK OPTIONS
C     ******************************************************************
C
      SUBROUTINE DK6AL(ISUM,LENX,NSCOL,NSROW,NSLAY,
     1              IDKTIM,
     2   LSDKZO,IDKZO,LSDKFO,IDKFO,
     3   LSDKZS,IDKZS,LSDKFS,IDKFS,
     4   INDK,IOUT,IOUTS)   
C
C     ******************************************************************
C
C  ARRAYS ARE ENTIRE TRANSPORT SUBGRID, IF NEEDED
      ISIZ=NSCOL*NSROW*NSLAY
      ISUMIN=ISUM
      LSDKZO=ISUM
      IF(IDKZO.EQ.1) ISUM=ISUM+ISIZ
      LSDKFO=ISUM
      IF(IDKFO.EQ.1) ISUM=ISUM+ISIZ
      LSDKZS=ISUM
      IF(IDKZS.EQ.1) ISUM=ISUM+ISIZ
      LSDKFS=ISUM
      IF(IDKFS.EQ.1) ISUM=ISUM+ISIZ
C
      ISP=ISUM-ISUMIN
      WRITE(IOUTS,101) ISP
  101 FORMAT(1X,I8,' ELEMENTS IN X ARRAY ARE USED BY DK')
      ISUM1=ISUM-1
      WRITE(IOUT,102) ISUM1,LENX
  102 FORMAT(1X,I8,' ELEMENTS OF X ARRAY USED OUT OF',I8)
      IF(ISUM1.GT.LENX) WRITE(IOUT,103)
  103 FORMAT(1X,'   ***X ARRAY MUST BE DIMENSIONED LARGER***')
C
      RETURN
      END
C
C  DK6RP READ DK INPUT FILE 
C     ******************************************************************
C
      SUBROUTINE DK6RP(RF,DKZO,DKFO,DKZS,DKFS,
     1                 IDKRF,IDKZO,IDKFO,IDKZS,IDKFS,
     2                  KPER,IOUTS,INDK,
     3     NSCOL,NSROW,NSLAY,NODESS,MOCTYPE)
C
C     ******************************************************************
C     READ DK PARAMETERS FOR A STRESS PERIOD
C     ******************************************************************
C
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*24 ANAME(5)
C
      DIMENSION RF(NODESS),DKZO(NODESS),
     *  DKFO(NODESS)
      DIMENSION DKZS(NODESS),DKFS(NODESS)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
      DATA ANAME(1) /'  SPACE VAR. RETARD FCTR'/
      DATA ANAME(2) /'  ZERO-ORDER GROWTH RATE'/
      DATA ANAME(3) /' FIRST-ORDER DECAY COEF.'/
      DATA ANAME(4) /'ZERO-ORDER SORBED GROWTH'/
      DATA ANAME(5) /'FIRST-ORDER SORBED DECAY'/
C     ------------------------------------------------------------------
      NSCR=NSCOL*NSROW
C
      WRITE(IOUTS,*)
      WRITE(IOUTS,*) ' SIMPLE REACTION (DK) PACKAGE INPUT'
      IF(KPER.GT.1) WRITE(IOUTS,'(2X,17HFOR STRESS PERIOD,I6)') KPER
C
C4------READ RETARDATION FACTOR
      IF(IDKRF.EQ.1.AND.KPER.EQ.0) THEN
        DO 100 KS=1,NSLAY
C  KK IS NUMBER OF FLOW LAYER
          KK=KS+ISLAY1-1
          LOC=1+(KS-1)*NSCR
          CALL U2DREL(RF(LOC),ANAME(1),NSROW,NSCOL,KK,INDK,IOUTS)
  100   CONTINUE
      END IF
C
C-------READ RATE COEFFICIENTS LAYER BY LAYER
      DO 200 KS=1,NSLAY
C  KK IS NUMBER OF FLOW LAYER
      KK=KS+ISLAY1-1
      LOC=1+(KS-1)*NSCR
      IF(IDKFO.EQ.1)
     *CALL U2DREL(DKFO(LOC),ANAME(3),NSROW,NSCOL,KK,INDK,IOUTS)
      IF(IDKFS.EQ.1)
     *CALL U2DREL(DKFS(LOC),ANAME(5),NSROW,NSCOL,KK,INDK,IOUTS)
      IF(IDKZO.EQ.1) THEN
        CALL U2DREL(DKZO(LOC),ANAME(2),NSROW,NSCOL,KK,INDK,IOUTS)
C  CHECK FOR ZERO-ORDER LOSS; PRINT WARNING IF MOCTYPE NE 1
        IF(MOCTYPE.NE.1) THEN
          ILOSS=0
          LOC2=LOC-1
          DO 150 IS=1,NSROW
          DO 150 JS=1,NSCOL
            LOC2=LOC2+1
            IF(DKZO(LOC2).LT.0.0) ILOSS=1
150       CONTINUE
          IF(ILOSS.EQ.1) WRITE(IOUTS,4000)
4000  FORMAT(' **** WARNING **** Zero-order loss (DKZO<0) ',
     *'is nonlinear and may cause oscillation and negative ',
     *'concentrations'/
     *'   Zero-order loss will be set to zero where Conc < 0')
        END IF
      END IF
      IF(IDKZS.EQ.1) THEN
        CALL U2DREL(DKZS(LOC),ANAME(4),NSROW,NSCOL,KK,INDK,IOUTS)
C  CHECK FOR ZERO-ORDER LOSS; PRINT WARNING IF MOCTYPE NE 1
        IF(MOCTYPE.NE.1) THEN
          ILOSS=0
          LOC2=LOC-1
          DO 160 IS=1,NSROW
          DO 160 JS=1,NSCOL
            LOC2=LOC2+1
            IF(DKZS(LOC2).LT.0.0) ILOSS=1
160       CONTINUE
          IF(ILOSS.EQ.1) WRITE(IOUTS,4010)
4010  FORMAT(' **** WARNING **** Zero-order loss (DKZS<0) ',
     *'is nonlinear and may cause oscillation and negative ',
     *'concentrations'/
     *'   Zero-order loss will be set to 0.0 where Conc < 0')
        END IF
      END IF
  200 CONTINUE
C
      RETURN
      END       
C
C  DK6DK FIRST-ORDER DECAY ON PARTICLES, COMPUTE DECAY FACTOR
C  IN AN INDIVIDUAL CELL
C
C     ******************************************************************
C
      SUBROUTINE DK6DK(DCYT,DKFO,DKFS,RF,TIMV,
     *       IDKFO,IDKFS,
     *       JS,IS,KS,NSCOL,NSROW,NSLAY)
C
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DOUBLE PRECISION DCYFCT,DCYT
C
      DIMENSION DKFO(NSCOL,NSROW,NSLAY),
     *  DKFS(NSCOL,NSROW,NSLAY)
      DIMENSION RF(NSCOL,NSROW,NSLAY)
C     ------------------------------------------------------------------
C
C  COMPUTE DECAY TERMS
      RFP=RF(JS,IS,KS)
      DCYT=1.D0
      IF(IDKFO.EQ.1) THEN
        IF(IDKFS.EQ.1) THEN
          DECAYP=(DKFO(JS,IS,KS)+(RFP-1.0)*DKFS(JS,IS,KS))/RFP
        ELSE
          DECAYP=DKFO(JS,IS,KS)
        END IF
      ELSE
        DECAYP=(RFP-1.0)*DKFS(JS,IS,KS)/RFP
      END IF
      DCYFCT=TIMV*DECAYP
      DCYT=DEXP(-DCYFCT)
C
      RETURN
      END       
C
C  DK6AP  CHANGE IN CONC DUE TO ZERO-ORDER GROWTH  EXPLICIT DISP
C********************************************************
C
      SUBROUTINE DK6AP(IBOUND,CNCNC,
     *         CAVG,DKZO,DKZS,RF,THCK,POR,
     *         DKZOIN,DKZOOUT,DKZSIN,DKZSOUT,
     *         NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,KSTP,KPER,
     *         IMOV,TIMV,IDKZO,IDKZS)
C
C********************************************************
C
      DIMENSION IBOUND(NCOL,NROW,NLAY),RF(NSCOL,NSROW,NSLAY),
     *  CNCNC(NSCOL,NSROW,NSLAY),CAVG(NSCOL,NSROW,NSLAY)
      DIMENSION DKZO(NSCOL,NSROW,NSLAY),DKZS(NSCOL,NSROW,NSLAY)
      DIMENSION THCK(NSCOL,NSROW,NSLAY),POR(NSCOL,NSROW,NSLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C3D  CHANGES CONCENTRATIONS DUE TO ZERO-ORDER SINK/SOURCES FOR ONE STEP ONLY
C
C  INITIALIZE GLOBAL SINK/SOURCE BUDGET
      DKZOIN=0.0
      DKZOOUT=0.0
      DKZSIN=0.0
      DKZSOUT=0.0
      AREA=CDEL*RDEL
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           ZERO-ORDER SINK/SOURCES IN WATER AND SORBED PHASES
C           ...WITH DECAY OF RECHARGE DURING TIME INCREMENT
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
      VOLUME=THCK(JS,IS,KS)*POR(JS,IS,KS)*AREA
C
      IF(IDKZO.EQ.1) THEN
        IF(IDKZS.EQ.1) THEN
          ZOM=TIMV*DKZO(JS,IS,KS)
          ZOCNC=ZOM/RF(JS,IS,KS)
          ZSM=TIMV*(RF(JS,IS,KS)-1.0)*DKZS(JS,IS,KS)
          ZSCNC=ZSM/RF(JS,IS,KS)
C  BOTH ZO AND ZS
          ZOSCNC=ZOCNC+ZSCNC
          IF(ZOSCNC.LT.0.0) THEN
            CTEST1=CAVG(JS,IS,KS)+CNCNC(JS,IS,KS)
C  SKIP REST IF CONCENTRATION LESS THAN ZERO
            IF(CTEST1.LE.0.0) GO TO 65
            CTEST2=CTEST1+ZOSCNC
            IF(CTEST2.LT.0.0) THEN
C  ADJUST FOR NEGATIVE CONCS
              RATIO=-CTEST1/ZOSCNC
              ZOM=ZOM*RATIO
              ZSM=ZSM*RATIO
              ZOCNC=ZOCNC*RATIO
              ZSCNC=ZSCNC*RATIO
            END IF
          END IF
          CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+ZOCNC+ZSCNC
          IF(ZOCNC.GT.0.0) THEN
            DKZOIN=DKZOIN+ZOM*VOLUME
          ELSE
            DKZOOUT=DKZOOUT+ZOM*VOLUME
          END IF
          IF(ZSCNC.GT.0.0) THEN
            DKZSIN=DKZSIN+ZSM*VOLUME
          ELSE
            DKZSOUT=DKZSOUT+ZSM*VOLUME
          END IF
        ELSE
C  ONLY ZO, USE SAME RATE FOR SORBED
          ZOCNC=TIMV*DKZO(JS,IS,KS)
          ZOM=ZOCNC*RF(JS,IS,KS)
          IF(ZOCNC.GT.0.0) THEN
            DKZOIN=DKZOIN+ZOM*VOLUME
              CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+ZOCNC
          ELSE
            CTEST1=CAVG(JS,IS,KS)+CNCNC(JS,IS,KS)
            IF(CTEST1.LE.0.0) GO TO 65
            CTEST2=CTEST1+ZOCNC
            IF(CTEST2.GE.0.0) THEN
              DKZOOUT=DKZOOUT+ZOM*VOLUME
              CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+ZOCNC
            ELSE
              DKZOOUT=DKZOOUT-CTEST1*VOLUME*RF(JS,IS,KS)
              CNCNC(JS,IS,KS)=-CAVG(JS,IS,KS)
            END IF
          END IF
        END IF
      ELSE
C  ONLY ZS, NO GROWTH IN DISSOLVED PHASE
        ZSM=TIMV*(RF(JS,KS,IS)-1.0)*DKZS(JS,IS,KS)
        ZSCNC=ZSM/RF(JS,IS,KS)
        IF(ZSCNC.GT.0.0) THEN
          DKZSIN=DKZSIN+ZSM*VOLUME
          CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+ZSCNC
        ELSE
          CTEST1=CAVG(JS,IS,KS)+CNCNC(JS,IS,KS)
          IF(CTEST1.LE.0.0) GO TO 65
          CTEST2=CTEST1+ZSCNC
          IF(CTEST2.GE.0.0) THEN
            DKZSOUT=DKZSOUT+ZSM*VOLUME
            CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+ZSCNC
          ELSE
            DKZSOUT=DKZSOUT-CTEST1*VOLUME*(RF(JS,IS,KS)-1.0)
            CNCNC(JS,IS,KS)=-CAVG(JS,IS,KS)
          END IF
        END IF
      END IF
C
   65 CONTINUE
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C  DK6IAP  CHANGE IN CONC DUE TO ZERO-ORDER GROWTH, IMPLICIT-DISP
C********************************************************
C
      SUBROUTINE DK6IAP(IBOUND,CAVG,DKZO,DKZS,RF,
     *         THCK,POR,DKZOIN,DKZOOUT,DKZSIN,DKZSOUT,
     *         NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,KSTP,KPER,
     *         TIMV,IDKZO,IDKZS,
     *  RS,FDTMTH,VAD,MRNO,NODESS,IIMOV)
C
C********************************************************
C
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *  THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *  CAVG(NSCOL,NSROW,NSLAY)
      DIMENSION DKZO(NSCOL,NSROW,NSLAY),DKZS(NSCOL,NSROW,NSLAY)
      DIMENSION RS(NODESS),VAD(NODESS),MRNO(*)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C3D  CHANGES CONCENTRATIONS DUE TO ZERO-ORDER GROWTH/LOSS FOR ONE STEP ONLY
C
C  INITIALIZE GLOBAL SINK/SOURCE BUDGET
      DKZOIN=0.0
      DKZOOUT=0.0
      DKZSIN=0.0
      DKZSOUT=0.0
      AREA=CDEL*RDEL
C
C     ***************************************************************
C
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
      VOLUME=THCK(JS,IS,KS)*POR(JS,IS,KS)*AREA
C.....Calculate natural node number
      M=(KS-1)*NSROW*NSCOL+(IS-1)*NSCOL+JS
C
      IF(IDKZO.EQ.1) THEN
        IF(IDKZS.EQ.1) THEN
          ZOM=TIMV*DKZO(JS,IS,KS)
          ZOCNC=ZOM/RF(JS,IS,KS)
          ZSM=TIMV*(RF(JS,KS,IS)-1.0)*DKZS(JS,IS,KS)
          ZSCNC=ZSM/RF(JS,IS,KS)
C  BOTH ZO (IN SOLUTION) AND ZS (SORBED)
          ZOSCNC=ZOCNC+ZSCNC
          IF(ZOSCNC.LT.0.0) THEN
            CTEST1=CAVG(JS,IS,KS)
C  SKIP REST IF CONCENTRATION LESS THAN ZERO
            IF(CTEST1.LE.0.0) GO TO 65
          END IF
          IF(ZOCNC.GT.0.0) THEN
            DKZOIN=DKZOIN+ZOM*VOLUME
          ELSE
            DKZOOUT=DKZOOUT+ZOM*VOLUME
          END IF
          IF(ZSCNC.GT.0.0) THEN
            DKZSIN=DKZSIN+ZSM*VOLUME
          ELSE
            DKZSOUT=DKZSOUT+ZSM*VOLUME
          END IF
        ELSE
C  ONLY ZO, USE SAME RATE FOR SORBED
          ZOSCNC=TIMV*DKZO(JS,IS,KS)
          ZOM=ZOSCNC*RF(JS,IS,KS)
          IF(ZOSCNC.GT.0.0) THEN
            DKZOIN=DKZOIN+ZOM*VOLUME
          ELSE
            CTEST1=CAVG(JS,IS,KS)
            IF(CTEST1.LE.0.0) GO TO 65
              DKZOOUT=DKZOOUT+ZOM*VOLUME
          END IF
        END IF
      ELSE
C  ONLY ZS, NO GROWTH IN DISSOLVED PHASE
        ZSM=TIMV*(RF(JS,KS,IS)-1.0)*DKZS(JS,IS,KS)
        IF(ZSM.GT.0.0) THEN
          DKZSIN=DKZSIN+ZSM*VOLUME
        ELSE
          CTEST1=CAVG(JS,IS,KS)
          IF(CTEST1.LE.0.0) GO TO 65
            DKZSOUT=DKZSOUT+ZSM*VOLUME
        END IF
        ZOSCNC=ZSM/RF(JS,IS,KS)
      END IF
C
   40 RS(M)=RS(M)+ZOSCNC
C
   65 CONTINUE
C
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C DK6MB   CALCULATE SOLUTE MASS BALANCE
C
C     **************************************************************** 
C
      SUBROUTINE DK6MB(SBVL,DKZOIN,DKZOOUT,DKZSIN,DKZSOUT,IDKZO,IDKZS,
     * NIUNIT)
C
C     **************************************************************** 
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION SBVL(6,NIUNIT)
C
      IF(IDKZO.EQ.1) THEN
        SBVL(1,13)=SBVL(1,13)+DKZOIN
        SBVL(2,13)=SBVL(2,13)+DKZOOUT
      END IF
      IF(IDKZS.EQ.1) THEN
        SBVL(1,15)=SBVL(1,15)+DKZSIN
        SBVL(2,15)=SBVL(2,15)+DKZSOUT
      END IF
      RETURN
      END
