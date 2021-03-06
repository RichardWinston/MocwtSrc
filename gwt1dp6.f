C     Last change:  RBW, Dec. 8, 2014
C  DP6DF  READ DP OPTION 
C
C     ******************************************************************
C
      SUBROUTINE DP6DF(INDP,IDPZO,IDPFO,IDPTIM,IDPPS,IOUT,IOUTS)
C
C     ******************************************************************
C
C     READ IN DP (DOUBLE POROSITY) OPTIONS
C     ******************************************************************
C
      WRITE(IOUTS,1000) INDP
      READ(INDP,*) IDPFO,IDPZO,IDPTIM,IDPPS
      IF(IDPFO.NE.1) THEN
       IDPFO=0
       WRITE(IOUTS,1020)
      ELSE
       WRITE(IOUTS,1021)
      END IF
      IF(IDPZO.NE.1) THEN
       IDPZO=0
       WRITE(IOUTS,1010)
      ELSE
       WRITE(IOUTS,1011)
      END IF
        IF(IDPTIM.NE.1) THEN
          IDPTIM=0
          WRITE(IOUTS,1040)
        ELSE
          WRITE(IOUTS,1041)
        END IF
      IF(IDPPS.EQ.1) THEN
        WRITE(IOUTS,1030)
      ELSE IF(IDPPS.EQ.2) THEN
        WRITE(IOUTS,1031)
      ELSE IF(IDPPS.EQ.3) THEN
        WRITE(IOUTS,1032)
      ELSE
        IDPPS=0
        WRITE(IOUTS,1033)
      END IF
      WRITE(IOUTS,*)
      RETURN
1000  FORMAT(/' INPUT FOR DP OPTIONS READ FROM UNIT',I10)
1010  FORMAT(' NO ZERO-ORDER REACTION IN DOUBLE POROSITY')
1011  FORMAT(' ZERO-ORDER REACTION IN DOUBLE POROSITY')
1020  FORMAT(' NO FIRST-ORDER REACTION IN DOUBLE POROSITY')
1021  FORMAT(' FIRST-ORDER REACTION IN DOUBLE POROSITY')
1030  FORMAT(' DOUBLE POROSITY CONCENTRATIONS WILL BE PRINTED ',
     *'USING CONCENTRATION FORMATS AND FREQUENCY') 
1031  FORMAT(' DOUBLE POROSITY CONCENTRATIONS WILL BE SAVED USING ',
     *'CONCENTRATION SAVE FREQUENCY') 
1032  FORMAT(' DOUBLE POROSITY CONCENTRATIONS WILL BE PRINTED AND ',
     *'SAVED USING CONCENTRATION FORMATS AND FREQUENCY') 
1033  FORMAT(' DOUBLE POROSITY CONCENTRATIONS WILL NOT BE PRINTED ',
     *'OR SAVED') 
1040  FORMAT(' DOUBLE POROSITY COEFFICIENTS WILL NOT CHANGE IN TIME.')
1041  FORMAT(' DOUBLE POROSITY COEFFICIENTS WILL CHANGE IN TIME.')
      END
C
C  MOC1AL  ALLOCATE SPACE IN X ARRAY FOR GWT ARRAYS
C***************************************************************
C
      SUBROUTINE DP6AL(ISUM,LENX,
     *   nscol,nsrow,nslay,
     *   LSDPCON,LSDPXRAT,LSDPPOR,
     *   LSDPZO,IDPZO,LSDPFO,IDPFO,
     *   iout,iouts)
C
C***************************************************************
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C-----VERSION
C     ******************************************************************
C     ALLOCATE ARRAY STORAGE FOR DOUBLE POROSITY PACKAGE
C     ******************************************************************
C
C4------COMPUTE DIMENSIONS FOR ARRAYS.
      WRITE(IOUTS,*) ' DP6AL'
      WRITE(IOUT ,*) ' DP6AL'
      NSRC=NSROW*NSCOL
      ISSIZ=NSRC*NSLAY
C
      ISOLD=ISUM
C
      LSDPCON=isum
      isum=isum+issiz
      LSDPXRAT=isum
      isum=isum+issiz
      LSDPPOR=isum
      isum=isum+issiz
      LSDPZO=ISUM
      IF(IDPZO.EQ.1) ISUM=ISUM+ISSIZ
      LSDPFO=ISUM
      IF(IDPFO.EQ.1) ISUM=ISUM+ISSIZ
C
C6------PRINT THE AMOUNT OF SPACE USED BY THE MOC PACKAGE.
      ISP=ISUM-ISOLD
      WRITE(IOUT,101) ISP
      WRITE(IOUTS,101) ISP
  101 FORMAT(1X,I8,' ELEMENTS IN X ARRAY ARE USED BY DP')
      ISUM1=ISUM-1
      WRITE(IOUT,102) ISUM1,LENX
  102 FORMAT(1X,I8,' ELEMENTS OF X ARRAY USED OUT OF',I8)
      IF(ISUM1.GT.LENX) WRITE(IOUT,103)
  103 FORMAT(1X,'   ***X ARRAY MUST BE DIMENSIONED LARGER***')
C
C  ZERO BUDGET TERMS
      BDDPI=0.0
      BDDPO=0.0
C
C7------RETURN
      RETURN
      END
C
C     ******************************************************************
C
      SUBROUTINE DP6RP(IBOUND,
     *   DPCON,DPXRAT,DPPOR,
     *   DPZO,DPFO,IDPZO,IDPFO,
     *   NSCOL,NSROW,NSLAY,nodess,
     *   NCOL,NROW,NLAY,inDP,KPER,IOUT,IOUTS,MOCTYPE)
C
C     ******************************************************************
C
C
C-----VERSION
C     ******************************************************************
C     READ AND INITIALIZE GWT SOLUTE TRANSPORT ARRAYS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*24 ANAME(5)
C
      DIMENSION IBOUND(NCOL,NROW,NLAY)
      DIMENSION DPCON(NODESS),DPXRAT(NODESS),
     *  DPPOR(NODESS)
      DIMENSION DPZO(NODESS),DPFO(NODESS)
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
      DATA ANAME(1) /'DOUBLE POROSITY T=0 CONC'/
      DATA ANAME(2) /'DOUBLE POROSITY POROSITY'/
      DATA ANAME(3) /'DBL PORO XCHNG RATE COEF'/
      DATA ANAME(4) /'DOUBLE PORSTY DECAY COEF'/
      DATA ANAME(5) /'DP ZERO-ORDR GROWTH RATE'/
C     ------------------------------------------------------------------
      NSCR=NSCOL*NSROW
      NCR=NCOL*NROW
C
      WRITE(IOUTS,*)
      WRITE(IOUTS,*) ' DOUBLE POROSITY PACKAGE INPUT'
C
C4------READ STARTING CONCENTRATIONS AND POROSITY
      IF(KPER.EQ.0) THEN
        DO 300 KS=1,NSLAY
C  KK IS NUMBER OF FLOW LAYER
          KK=KS+ISLAY1-1
          LOC=1+(KS-1)*NSCR
          CALL U2DREL(DPCON(LOC),ANAME(1),NSROW,NSCOL,KK,INDP,IOUTS)
          CALL U2DREL(DPPOR(LOC),ANAME(2),NSROW,NSCOL,KK,INDP,IOUTS)
  300   CONTINUE
        END IF
        DO 310 KS=1,NSLAY
        KK=KS+ISLAY1-1
        LOC=1+(KS-1)*NSCR
        CALL U2DREL(DPXRAT(LOC),ANAME(3),NSROW,NSCOL,KK,INDP,IOUTS)
        IF(IDPFO.EQ.1)
     *    CALL U2DREL(DPFO(LOC),ANAME(4),NSROW,NSCOL,KK,INDP,IOUTS)
        IF(IDPZO.EQ.1) THEN
          CALL U2DREL(DPZO(LOC),ANAME(5),NSROW,NSCOL,KK,INDP,IOUTS)
C  CHECK FOR ZERO-ORDER LOSS; PRINT WARNING IF MOCTYPE NE 1
          IF(MOCTYPE.NE.1) THEN
              ILOSS=0
              LOC2=LOC-1
              DO 150 IS=1,NSROW
              DO 150 JS=1,NSCOL
                LOC2=LOC2+1
                IF(DPZO(LOC2).LT.0.0) ILOSS=1
150         CONTINUE
            IF(ILOSS.EQ.1) WRITE(IOUTS,4010)
4010  FORMAT(' **** WARNING **** Zero-order loss (DPZO<0) ',
     *'is nonlinear and may cause oscillation and negative ',
     *'concentrations'/
     *'   Zero-order loss will be set to 0.0 where Conc < 0')
            END IF
          END IF
  310 CONTINUE
C
C5------SET CONC IN NO FLOW CELLS TO CNOFLO
      LOCS=0
      DO 403 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 402 IS=1,NSROW
      I=IS+ISROW1-1
      DO 401 JS=1,NSCOL
      J=JS+ISCOL1-1
      LOCS=LOCS+1
      IF(IBOUND(J,I,K).EQ.0) THEN
        DPCON(LOCS)=CNOFLO
        DPXRAT(LOCS)=0.0
        DPPOR(LOCS)=0.0
        IF(IDPFO.EQ.1) DPFO(LOCS)=0.0
        IF(IDPZO.EQ.1) DPZO(LOCS)=0.0
      END IF
  401 CONTINUE
  402 CONTINUE
  403 CONTINUE
C
C9------RETURN
 1000 RETURN
      END
C
C  DP6AP  CHANGE IN CONC DUE TO DOUBLE POROSITY
C********************************************************
C
      SUBROUTINE DP6AP(IBOUND,
     *  DPCON,DPXRAT,DPPOR,DPZO,DPFO,
     *  RF,CAVG,THCK,POR,CNCNC,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *  KSTP,KPER,IOUTS,IMOV,TIMV,
     *  RMASS,ZODPIN,ZODPOUT,FODPIN,FODPOUT,
     *  IDPZO,IDPFO)
C
C********************************************************
C
C  CHANGES CONCENTRATIONS DUE TO DOUBLE POROSITY FOR ONE STEP ONLY
C
      DIMENSION IBOUND(NCOL,NROW,NLAY)
      DIMENSION DPCON(NSCOL,NSROW,NSLAY),
     *  DPXRAT(NSCOL,NSROW,NSLAY),DPPOR(NSCOL,NSROW,NSLAY),
     *  DPZO(NSCOL,NSROW,NSLAY),DPFO(NSCOL,NSROW,NSLAY)
      DIMENSION RF(NSCOL,NSROW,NSLAY),
     *  CAVG(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),CNCNC(NSCOL,NSROW,NSLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           DOUBLE POROSITY
C  NO ADJUSTMENTS FOR DECAY
C
      RMASS=0.0
      ZODPIN=0.0
      ZODPOUT=0.0
      FODPIN=0.0
      FODPOUT=0.0
C
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
C  SKIP IF ZERO POROSITY FOR IMMOBILE ZONE
      IF(DPPOR(JS,IS,KS).EQ.0.0) GO TO 65
      C1=CAVG(JS,IS,KS)
C  LIMIT CONCENTRATION GRADIENT BY C=0
      IF(C1.LT.0.0) C1=0.0
C  CRANK-NICOLSON (SP?) FD
      RATERM=DPXRAT(JS,IS,KS)/DPPOR(JS,IS,KS)
      ANEW=1.0/TIMV+RATERM/2.0
      AOLD=1.0/TIMV-RATERM/2.0
      IF(IDPFO.EQ.1) THEN
        ANEW=ANEW+DPFO(JS,IS,KS)/2.0
        AOLD=AOLD-DPFO(JS,IS,KS)/2.0
      END IF
      IF(IDPZO.EQ.1) THEN
        ZORAT=DPZO(JS,IS,KS)
        CMNEW=(DPCON(JS,IS,KS)*AOLD+RATERM*C1+ZORAT)/ANEW
        IF(ZORAT.LT.0.0) THEN
          IF(CMNEW.LT.0.0) THEN
C  RESET CONCENTRATION TO ZERO IF NEGATIVE DUE TO ZERO-ORDER SINK
            CMNEW=0.0
C  ADJUST ZORAT FOR MASS BALANCE
            ZORAT=-AOLD*DPCON(JS,IS,KS)-RATERM*C1
          END IF
C  MASS BALANCE
          ZODPOUT=ZODPOUT+ZORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        ELSE
          ZODPIN=ZODPIN+ZORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        END IF
      ELSE
        CMNEW=(DPCON(JS,IS,KS)*AOLD+RATERM*C1)/ANEW
      END IF
      CMAVG=(CMNEW+DPCON(JS,IS,KS))/2.0
C  FIRST-ORDER MASS BALANCE
      IF(IDPFO.EQ.1) THEN
        FORAT=-DPFO(JS,IS,KS)*CMAVG
        IF(FORAT.GT.0.0) THEN
          FODPIN=FODPIN+FORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        ELSE
          FODPOUT=FODPOUT+FORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        END IF
      END IF
C  CURRENT MASS IN DOUBLE POROSITY
      RMASS=RMASS+CMNEW*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
C
C    xchng is the mass per unit vol aquifer exchanged
C     between the IMMOBILE POROSITY and AQUIFER.  Postive for 
C     mass increase to AQUIFER and decrease in IMMOBILE POROSITY
C
C  COMPUTER TERMS FOR MAIN CONCENTRATION SOLUTION
      XCHNG=TIMV*DPXRAT(JS,IS,KS)*( CMAVG - C1 )
c     compute change in concentration of aquifer, explicit
      CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+XCHNG/
     *       (POR(JS,IS,KS)*RF(JS,IS,KS))
C
      DPCON(JS,IS,KS)=CMNEW
C
   65 CONTINUE
C
C  DOUBLE POROSITY MASS BUDGET
      RMASS=RMASS*CDEL*RDEL
      IF(IDPZO.EQ.1) THEN
        ZODPIN=ZODPIN*CDEL*RDEL*TIMV
        ZODPOUT=ZODPOUT*CDEL*RDEL*TIMV
      END IF
      IF(IDPFO.EQ.1) THEN
        FODPIN=FODPIN*CDEL*RDEL*TIMV
        FODPOUT=FODPOUT*CDEL*RDEL*TIMV
      END IF
C
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C  DP6IAP  CHANGE IN CONC DUE TO DOUBLE POROSITY
C********************************************************
C
      SUBROUTINE DP6IAP(IBOUND,
     *  DPCON,DPXRAT,DPPOR,DPZO,DPFO,
     *  RF,CAVG,THCK,POR,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *  KSTP,KPER,IOUTS,IMOV,TIMV,
     *  ZODPIN,ZODPOUT,IDPZO,IDPFO,
     *  RS,FDTMTH,VAD,MRNO,NODESS)
C
C********************************************************
C
C  CHANGES CONCENTRATIONS DUE TO DOUBLE POROSITY FOR ONE STEP ONLY
C
      DIMENSION IBOUND(NCOL,NROW,NLAY)
      DIMENSION DPCON(NSCOL,NSROW,NSLAY),
     *  DPXRAT(NSCOL,NSROW,NSLAY),DPPOR(NSCOL,NSROW,NSLAY),
     *  DPZO(NSCOL,NSROW,NSLAY),DPFO(NSCOL,NSROW,NSLAY)
      DIMENSION RF(NSCOL,NSROW,NSLAY),
     *  CAVG(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY)
      DIMENSION RS(NODESS),VAD(NODESS),MRNO(*)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           DOUBLE POROSITY
C  NO ADJUSTMENTS FOR DECAY
C
      ZODPIN=0.0
      ZODPOUT=0.0
C
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
C  SKIP IF ZERO POROSITY FOR IMMOBILE ZONE
      IF(DPPOR(JS,IS,KS).EQ.0.0) GO TO 65
      C1=CAVG(JS,IS,KS)
C  CRANK-NICOLSON FD
      RATERM=DPXRAT(JS,IS,KS)/DPPOR(JS,IS,KS)
      ANEW=1.0/TIMV+RATERM/2.0
      AOLD=1.0/TIMV-RATERM/2.0
      IF(IDPFO.EQ.1) THEN
        ANEW=ANEW+DPFO(JS,IS,KS)/2.0
        AOLD=AOLD-DPFO(JS,IS,KS)/2.0
      END IF
      BTRFP=DPXRAT(JS,IS,KS)*TIMV/(RF(JS,IS,KS)*POR(JS,IS,KS))
        BTRFP2=BTRFP/2.D0
        FTERM=RATERM/(2.D0*ANEW)
C.....Calculate natural node number
      M=(KS-1)*NSROW*NSCOL+(IS-1)*NSCOL+JS
C.....Add double porosity term for right hand side vector
        RS(M)=RS(M)+BTRFP*CAVG(JS,IS,KS)*(FTERM-1.D0) +
     *           BTRFP2*DPCON(JS,IS,KS)*(1.D0+AOLD/ANEW)
      IF(IDPZO.EQ.1) THEN
        ZORAT=DPZO(JS,IS,KS)
        IF(ZORAT.LT.0.0) THEN
          IF(DPCON(JS,IS,KS).LE.0.0) THEN
C  TURN OFF LOSS IF DP CONC LE ZERO
            ZORAT=0.D0
          END IF
C  EXPLICIT MASS BALANCE FOR ZERO-ORDER GROWTH/LOSS
          ZODPOUT=ZODPOUT+ZORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        ELSE
          ZODPIN=ZODPIN+ZORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        END IF
        RS(M)=RS(M)+BTRFP2*ZORAT/ANEW
      END IF
C
      IF(IMOV.EQ.1) THEN
        MA=MRNO(M)
          VAD(MA)=VAD(MA)+BTRFP*FDTMTH*(1.D0-FTERM)
      END IF
   65 CONTINUE
C
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C  DP6IUP  UPDATE CONC DUE TO DOUBLE POROSITY
C********************************************************
C
      SUBROUTINE DP6IUP(IBOUND,
     *  DPCON,DPXRAT,DPPOR,DPZO,DPFO,
     *  RF,CAVG,CNCNC,THCK,POR,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *  KSTP,KPER,IOUTS,IMOV,TIMV,
     *  RMASS,FODPIN,FODPOUT,ZODPIN,ZODPOUT,
     *  IDPZO,IDPFO,
     *  FDTMTH)
C
C********************************************************
C
C  UPDATE DOUBLE POROSITY CONCENTRATIONS (IMPLICIT-DISPERSION)
C
      DIMENSION IBOUND(NCOL,NROW,NLAY)
      DIMENSION DPCON(NSCOL,NSROW,NSLAY),
     *  DPXRAT(NSCOL,NSROW,NSLAY),DPPOR(NSCOL,NSROW,NSLAY),
     *  DPZO(NSCOL,NSROW,NSLAY),DPFO(NSCOL,NSROW,NSLAY)
      DIMENSION RF(NSCOL,NSROW,NSLAY),
     *  CAVG(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),CNCNC(NSCOL,NSROW,NSLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           DOUBLE POROSITY
C  NO ADJUSTMENTS FOR DECAY
C
      RMASS=0.0
      FODPIN=0.0
      FODPOUT=0.0
C
C
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
C  SKIP IF ZERO POROSITY FOR IMMOBILE ZONE
      IF(DPPOR(JS,IS,KS).EQ.0.0) GO TO 65
      C1=CAVG(JS,IS,KS)+FDTMTH*CNCNC(JS,IS,KS)
C  CRANK-NICOLSON FD
      RATERM=DPXRAT(JS,IS,KS)/DPPOR(JS,IS,KS)
      ANEW=1.0/TIMV+RATERM/2.0
      AOLD=1.0/TIMV-RATERM/2.0
      IF(IDPFO.EQ.1) THEN
        ANEW=ANEW+DPFO(JS,IS,KS)/2.0
        AOLD=AOLD-DPFO(JS,IS,KS)/2.0
      END IF
      IF(IDPZO.EQ.1) THEN
        ZORAT=DPZO(JS,IS,KS)
        IF(ZORAT.LT.0.0) THEN
          IF(DPCON(JS,IS,KS).LE.0.0) THEN
C  TURN OFF LOSS IF DP CONC LE ZERO
            ZORAT=0.D0
          END IF
C  EXPLICIT MASS BALANCE FOR ZERO-ORDER GROWTH/LOSS
          ZODPOUT=ZODPOUT+ZORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        ELSE
          ZODPIN=ZODPIN+ZORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        END IF
        ELSE
        ZORAT=0.D0
      END IF
      CMNEW=(RATERM*(FDTMTH*CNCNC(JS,IS,KS)+CAVG(JS,IS,KS)) +
     *  ZORAT+DPCON(JS,IS,KS)*AOLD)/ANEW
C
C  MASS BALANCE
      CMAVG=(CMNEW+DPCON(JS,IS,KS))/2.0
C  FIRST-ORDER MASS BALANCE
      IF(IDPFO.EQ.1) THEN
        FORAT=-DPFO(JS,IS,KS)*CMAVG
        IF(FORAT.GT.0.0) THEN
          FODPIN=FODPIN+FORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        ELSE
          FODPOUT=FODPOUT+FORAT*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
        END IF
      END IF
C  CURRENT MASS IN DOUBLE POROSITY
      RMASS=RMASS+CMNEW*DPPOR(JS,IS,KS)*THCK(JS,IS,KS)
C
      DPCON(JS,IS,KS)=CMNEW
   65 CONTINUE
C
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C***************************************************************
C
      SUBROUTINE SDP6C(DPCON,
     *    KSTP,NSTP,KPER,NPER,IMOV,NMOV,
     *    NSCOL,NSROW,NSLAY,IOUTS,JUNIT,PERTIM,TOTIM,SUMTCH,
     *    NPNTCL,ICONFM,IDPPS,NIUNIT)
C
C-----from VERSION 1653 15MAY1987 SBAS1H
C     *******************************************************
C     PRINT AND RECORD CONCS
C     *******************************************************
C
C        SPECIFICATIONS
C     -------------------------------------------------------
      DIMENSION DPCON(NSCOL,NSROW,NSLAY),
     *  JUNIT(NIUNIT)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     -------------------------------------------------------
C CHECK FLAGS FOR OUTPUT
      IPRNT=0
      IF (NPNTCL.EQ.0.AND.KPER.EQ.NPER.AND.KSTP.EQ.NSTP.
     *AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-2.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-1.AND.IMOV.EQ.NMOV) IPRNT=1
      IF(NPNTCL.GE.1) THEN
         IF(MOD(IMOV,NPNTCL).EQ.0) IPRNT=1
      ENDIF
      IF (KPER.EQ.NPER.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
C INITIAL CONDITION
      IF(IMOV.EQ.0) IPRNT=1
C SKIP IF NO OUTPUT
      IF (IPRNT.EQ.0) RETURN
C FOR EACH LAYER: PRINT DPCON IF REQUESTED.
      DO 80 KS=1,NSLAY
      KKS=KS
C
      IF(IDPPS.EQ.1.OR.IDPPS.EQ.3) THEN
       IF(JUNIT(2).LE.0) THEN
         IF(ICONFM.LT.0) CALL ULAPRS(DPCON(1,1,KKS),'DOUBLE POROSTY C',
     *       KSTP,KPER,NSCOL,NSROW,KKS,-ICONFM,IOUTS)
         IF(ICONFM.GE.0) CALL ULAPRW(DPCON(1,1,KKS),'DOUBLE POROSTY C',
     *       KSTP,KPER,NSCOL,NSROW,KKS,ICONFM,IOUTS)
       ELSE
C PRINT TO SEPARATE TEXT FILE
         IF(KS.EQ.1.AND.IMOV.EQ.0) WRITE(IOUTS,51) JUNIT(2)
         WRITE(JUNIT(2),55) IMOV, KSTP, KPER, SUMTCH
  51  FORMAT('DOUBLE POROSITY CONCENTRATION DATA WILL BE SAVED ON ',
     *  'UNIT ',I3,' IN ASCII FORMAT')
  55  FORMAT('DOUBLE POROSITY CONCENTRATIONS AT NODES IN SUBGRID.  ',
     * 'IMOV=',I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
         WRITE(JUNIT(2),*) 'SUBGRID LAYER ', KS
         DO 60 IS=1,NSROW
         DO 59 JS=1,NSCOL
  59  CONTINUE
           WRITE(JUNIT(2),'(1P10E12.4)') 
     *          (DPCON(JS,IS,KS),JS=1,NSCOL)
  60  CONTINUE
       ENDIF
      END IF
       IF(JUNIT(3).GT.0.AND.(IDPPS.EQ.2.OR.IDPPS.EQ.3)) THEN
         IF(KS.EQ.1.AND.IMOV.EQ.0) WRITE(IOUTS,65) JUNIT(3)
C PRINT BINARY FILE
  65  FORMAT('DOUBLE POROSITY CONCENTRATION DATA WILL BE SAVED ON ',
     * 'UNIT ',I3,' IN BINARY FORMAT')
         CALL ULASAV(DPCON(1,1,KKS),'DOUBLE POROSTY C',
     *     KSTP,KPER,SUMTCH,TOTIM,NSCOL,NSROW,KKS,JUNIT(3))
       ENDIF
  80  CONTINUE
C
   50 RETURN
      END
C
C DP6BD   CALCULATE SOLUTE MASS BALANCE
C
C     **************************************************************** 
C
      SUBROUTINE DP6MB(SBVL,
     *    RMASS,ZODPIN,ZODPOUT,FODPIN,FODPOUT,
     *    IDPZO,IDPFO,NIUNIT)
C
C     **************************************************************** 
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION SBVL(6,NIUNIT)
C
      SBVL(2,18)=RMASS
      IF(IDPZO.EQ.1) THEN
        SBVL(1,19)=SBVL(1,19)+ZODPIN
        SBVL(2,19)=SBVL(2,19)+ZODPOUT
      END IF
      IF(IDPFO.EQ.1) THEN
        SBVL(1,20)=SBVL(1,20)+FODPIN
        SBVL(2,20)=SBVL(2,20)+FODPOUT
      END IF
      RETURN
      END
C
C  DP6ST  STABILITY LIMIT FOR EXPLICIT
C********************************************************
C
      SUBROUTINE GWT1DP6ST(IBOUND,RF,POR,
     *    TIMDP,JMAXDP,IMAXDP,KMAXDP,
     *  DPXRAT,DPPOR,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS)
C
C********************************************************
C
C  SET STABILITY TIME LIMIT FOR DOUBLE POROSITY
C
      DIMENSION IBOUND(NCOL,NROW,NLAY)
      DIMENSION 
     *  DPXRAT(NSCOL,NSROW,NSLAY),DPPOR(NSCOL,NSROW,NSLAY)
      DIMENSION RF(NSCOL,NSROW,NSLAY),POR(NSCOL,NSROW,NSLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C     ***************************************************************
C     COMPUTE STABILITY TIME LIMIT, TIME FOR COMPLETE DEPLETION
C       OF SOLUTE FROM AQUIFER INTO DOUBLE POROSITY
C
      JMAXDP=0
      IMAXDP=0
      KMAXDP=0
      TIMDP=1.E30
C
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
C  SKIP IF ZERO POROSITY FOR IMMOBILE ZONE
C    OR ZERO RATE COEFFICIENT
      IF(DPPOR(JS,IS,KS).EQ.0.0.OR.
     *   DPXRAT(JS,IS,KS).EQ.0.0) GO TO 65
      TIM=POR(JS,IS,KS)*RF(JS,IS,KS)/DPXRAT(JS,IS,KS)
        TIM2=2.*DPPOR(JS,IS,KS)/DPXRAT(JS,IS,KS)
      IF(TIM.LT.TIMDP) THEN
        TIMDP=TIM
        JMAXDP=J
        IMAXDP=I
        KMAXDP=K
      END IF
C
   65 CONTINUE
C
C
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C
C     ******************************************************************
C
      SUBROUTINE GWT1DP6INIT(SBVL,IBOUND,THCK,
     *   DPCON,DPXRAT,DPPOR,
     *   DPZO,DPFO,IDPZO,IDPFO,
     *   NSCOL,NSROW,NSLAY,
     *   NCOL,NROW,NLAY,INDP,KPER,IOUT,IOUTS,NIUNIT)
C
C     ******************************************************************
C
C-----VERSION
C     ******************************************************************
C     INITIALIZE MASS IN DOUBLE POROSITY
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION SBVL(6,NIUNIT)
      DIMENSION IBOUND(NCOL,NROW,NLAY),THCK(NSCOL,NSROW,NSLAY)
      DIMENSION DPCON(NSCOL,NSROW,NSLAY),DPXRAT(NSCOL,NSROW,NSLAY),
     *  DPPOR(NSCOL,NSROW,NSLAY),
     *  DPZO(NSCOL,NSROW,NSLAY),DPFO(NSCOL,NSROW,NSLAY)
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ------------------------------------------------------------------
C
C5------SET CONC IN NO FLOW CELLS TO CNOFLO
C       SET INITIAL MASS IN DOUBLE POROSITY
      DO 403 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 402 IS=1,NSROW
      I=IS+ISROW1-1
      DO 401 JS=1,NSCOL
      J=JS+ISCOL1-1
      LOCS=LOCS+1
      IF(IBOUND(J,I,K).EQ.0) THEN
        DPCON(JS,IS,KS)=CNOFLO
        DPXRAT(JS,IS,KS)=0.0
        DPPOR(JS,IS,KS)=0.0
        IF(IDPFO.EQ.1) DPFO(JS,IS,KS)=0.0
        IF(IDPZO.EQ.1) DPZO(JS,IS,KS)=0.0
	ELSE
	  SBVL(1,18)=SBVL(1,18)+THCK(JS,IS,KS)*
     *    DPCON(JS,IS,KS)*DPPOR(JS,IS,KS)
      END IF
  401 CONTINUE
  402 CONTINUE
  403 CONTINUE
      SBVL(1,18)=SBVL(1,18)*CDEL*RDEL
C
C9------RETURN
 1000 RETURN
      END
