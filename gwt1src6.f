C  SRC6  TRANSPORT BC'S AT FLUID SINK/SOURCES
C       FIXED HEAD NODES
C
C  GWT1SRC6FM  INIALIZE FIXED HEADS, REST MODULAR FORM
C     ***************************************************************
C
      SUBROUTINE GWT1SRC6FM(
     *  EVTFLO,CHDFLO,
     *  BUFF,IBOUND,
     *  SRCFLO,SRCSOL,SNKFLO,CFXHBC,SRCDCY,
     *  NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,
     *  IOUTS,IFXHED,ICONLY,MOCTYPE)
C
C     ***************************************************************
C
      DIMENSION BUFF(NCOL,NROW,NLAY),IBOUND(NCOL,NROW,NLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),CFXHBC(NSCOL,NSROW,NSLAY),
     *  CHDFLO(NSCOL,NSROW,NSLAY),EVTFLO(NSCOL,NSROW,2)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH

C
C     ***************************************************************
cea
C  INITIALIZE MASS AND WATER VOLUME FOR THIS FLOW STEP
      IF(MOCTYPE.EQ.3) THEN
        WATVOL=0.0
        STINIT=0.0
        ADINIT=0.0
      ENDIF
C  FIXED HEAD CONC USE SUBGRID INDEX 
C  USE BUFF FOR FIXED HEAD FLUXES   
C  RETURN IF ICONLY = 1, NO FLUID SINK/SOURCES WITHIN SUBGRID
      IF(ICONLY.EQ.1) RETURN
C
C  INITIALIZE SINK/SOURCE ARRAYS
      DO 30 KS=1,NSLAY
      DO 20 IS=1,NSROW
      DO 12 JS=1,NSCOL
      SRCFLO(JS,IS,KS)=0.0D0
      SRCSOL(JS,IS,KS)=0.0D0
      SNKFLO(JS,IS,KS)=0.0D0
      CHDFLO(JS,IS,KS)=0.0D0
   12 CONTINUE
   20 CONTINUE
   30 CONTINUE
      DO 81 IS=1,NSROW
      DO 81 JS=1,NSCOL
      EVTFLO(JS,IS,1)=0.0
      EVTFLO(JS,IS,2)=1.0
   81 CONTINUE
cgzh allow srcdcy to accumulate
c      SRCDCY=0.0
C
C  RETURN IF NO FIXED HEAD NODES WITHIN TRANSPORT SUBGRID
      IF(IFXHED.LE.0) RETURN
      DO 31 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 21 IS=1,NSROW
      I=IS+ISROW1-1
      DO 11 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF CELL IS NOT FIXED HEAD NODE
      IF(IBOUND(J,I,K).GE.0) GO TO 11
C  SAVE FLUXES TO CHDFLO FOR MASS BALANCE   
      CHDFLO(JS,IS,KS)=BUFF(J,I,K)
C  FOR A SINK, ONLY NEED FLUID FLUX
      IF(BUFF(J,I,K).LT.0.0) THEN
         SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+BUFF(J,I,K)
C  FOR A SOURCE, NEED FLUID FLUX AND SOLUTE FLUX
      ELSE IF(BUFF(J,I,K).GT.0.0) THEN
         RATE=BUFF(J,I,K)
         SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+RATE
         RATEC=RATE*CFXHBC(JS,IS,KS)
         SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+RATEC
      END IF
   11 CONTINUE
   21 CONTINUE
   31 CONTINUE
C
      RETURN
      END
C
C
C  SRC6AP  CHANGE IN CONC DUE TO SINK/SOURCES
C********************************************************
C
      SUBROUTINE SRC6AP(IBOUND,
     *  SRCFLO,SRCSOL,SNKFLO,RF,
     *  CONC,THCK,POR,CNCNC,CNOLD,
     *  CAVG,EVTFLO,IEVT,SRCDCY,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NEVTOP,
     *  KSTP,KPER,IOUTS,IMOV,TIMV,DECAY,INAGE)
C
C********************************************************
C
      DOUBLE PRECISION DECAY
      DOUBLE PRECISION DCYT2,DCYFCT
C
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *  CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),CNCNC(NSCOL,NSROW,NSLAY),
     *  CNOLD(NSCOL,NSROW,NSLAY),CAVG(NSCOL,NSROW,NSLAY),
     *  EVTFLO(NSCOL,NSROW,2),IEVT(NROW,NCOL)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C3D  CHANGES CONCENTRATIONS DUE TO DISPERSION AND SOURCES FOR ONE STEP ONLY
C
      DCYT2=1.D0
      IF(DECAY.NE.0.D0) THEN
         DCYFCT=TIMV*DECAY
         DCYT2=EXP(-DCYFCT*0.5D0)
      END IF
      TOVERA=TIMV/(CDEL*RDEL)
C  INITIALIZE GLOBAL MASS FLUX RATE TO SINKS
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           MIXING AT SOURCE CELLS...
C           ...WITH DECAY OF RECHARGE DURING TIME INCREMENT
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
      PRBINV=1.D0/(POR(JS,IS,KS)*THCK(JS,IS,KS)*rf(js,is,ks))
      C1=CAVG(JS,IS,KS)
      IF(ABS(C1).LT.1.0E-20) C1=0.D0
      C2=C1
      KEV=1
      ETFLO=0.0
      ETCONC=0.0
C  IF AGE PACKAGE IS ACTIVE, DO NOT ALLOW ET TO CONCENTRATE AGE MASS 
      IF(INAGE.GT.0) ETCONC=C2
      KEV=EVTFLO(JS,IS,2)
      IF(KEV.EQ.K) ETFLO=EVTFLO(JS,IS,1)  
      DIV=SRCFLO(JS,IS,KS)+SNKFLO(JS,IS,KS)+ETFLO
      IF(DECAY.EQ.0.D0) GO TO 37
      IF(DCYT2.GT.0.9) GO TO 37
      IF(DIV.LE.0.0) GO TO 37
      IF(CNOLD(JS,IS,KS).LE.0.0.OR.CONC(JS,IS,KS).LE.0.0) GO TO 37
C     NEXT CALC IS EQUIVALENT TO C1=EXP((ALOG(CNOLD)+ALOG(CONC))*0.5)
      C1=SQRT(CNOLD(JS,IS,KS)*CONC(JS,IS,KS))
   37 CONTINUE
C  CONSISTENT WITH MOC-2D:
C  C1 IS ADJUSTED IF DECAY TIMES CRITICAL
C  C2 IS NOT ADJUSTED FOR DECAY
C  SRCSOL COULD BE SCALED AFTER VELO
C  CHANGE IN SIGN CONVENTION FOR FLUID SOURCES FROM MOC-2D
C  HERE SRCFLO IS POSITIVE FOR INFLOW TO AQUIFER
C
      SRC1=SRCSOL(JS,IS,KS)*DCYT2
cea   SRCDCY=SRCDCY+(SRC1-SRCSOL(JS,IS,KS))
c srcdcy now stores mass decayed instead of mass rate
      SRCDCY=SRCDCY+(SRC1-SRCSOL(JS,IS,KS))*TIMV
   40 CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+
     *      tovera*(-C1*DIV+SRC1+ETFLO*ETCONC
     *                  +SNKFLO(JS,IS,KS)*C2)*PRBINV
C
   65 CONTINUE
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C  SRC6IAP  CHANGE IN CONC DUE TO SINK/SOURCES for first step of implicit assembly
C********************************************************
C
      SUBROUTINE SRC6IAP(IBOUND,
     *  SRCFLO,SRCSOL,SNKFLO,RF,
     *  CONC,THCK,POR,CNOLD,
     *  CAVG,EVTFLO,IEVT,SRCDCY,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NEVTOP,
     *  TIMV,DECAY,
     X  RS,FDTMTH,VAD,MRNO,NODESS,IIMOV,INAGE)
C
C********************************************************
C
      DOUBLE PRECISION DECAY
      DOUBLE PRECISION DCYT,DCYT2,DCYFCT
c      DOUBLE PRECISION RS
C
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *  CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),
     *  CNOLD(NSCOL,NSROW,NSLAY),CAVG(NSCOL,NSROW,NSLAY),
     *  EVTFLO(NSCOL,NSROW,2),IEVT(NROW,NCOL)
      DIMENSION RS(NODESS),VAD(NODESS),MRNO(*)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C3D  CHANGES CONCENTRATIONS DUE TO DISPERSION AND SOURCES FOR ONE STEP ONLY
C
      DCYT=1.D0
      DCYT2=1.D0
      IF(DECAY.NE.0.D0) THEN
         DCYFCT=TIMV*DECAY
         DCYT2=EXP(-DCYFCT*0.5D0)
         DCYT=EXP(-DCYFCT)
      END IF
      TOVERA=TIMV/(CDEL*RDEL)
C  INITIALIZE GLOBAL MASS FLUX RATE TO SINKS
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           MIXING AT SOURCE CELLS...
C           ...WITH DECAY OF RECHARGE DURING TIME INCREMENT
      DO 65 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 65 IS=1,NSROW
      I=IS+ISROW1-1
      DO 65 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 65
C.....Calculate natural node number
      M=(KS-1)*NSROW*NSCOL+(IS-1)*NSCOL+JS
      PRBINV=1.D0/(POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
      C1=CAVG(JS,IS,KS)
      IF(ABS(C1).LT.1.0E-20) C1=0.D0
      C2=C1
      KEV=1
      ETFLO=0.0
      ETCONC=0.0
C  IF AGE PACKAGE IS ACTIVE, DO NOT ALLOW ET TO CONCENTRATE AGE MASS 
      IF(INAGE.GT.0) ETCONC=C2
      KEV=EVTFLO(JS,IS,2)
      IF(KEV.EQ.K) ETFLO=EVTFLO(JS,IS,1)  
      DIV=SRCFLO(JS,IS,KS)+SNKFLO(JS,IS,KS)+ETFLO
      IF(DECAY.EQ.0.D0) GO TO 37
      IF(DCYT2.GT.0.9) GO TO 37
      IF(DIV.LE.0.0) GO TO 37
      IF(CNOLD(JS,IS,KS).LE.0.0.OR.CONC(JS,IS,KS).LE.0.0) GO TO 37
C     NEXT CALC IS EQUIVALENT TO C1=EXP((ALOG(CNOLD)+ALOG(CONC))*0.5)
C     USED FOR >>> DECAY RATE
      C1=SQRT(CNOLD(JS,IS,KS)*CONC(JS,IS,KS))
   37 CONTINUE
C  CONSISTENT WITH MOC-2D:
C  C1 IS ADJUSTED IF DECAY TIMES CRITICAL
C  C2 IS NOT ADJUSTED FOR DECAY
C  SRCSOL COULD BE SCALED AFTER VELO
C  CHANGE IN SIGN CONVENTION FOR FLUID SOURCES FROM MOC-2D
C  HERE SRCFLO IS POSITIVE FOR INFLOW TO AQUIFER
C
      SRC1=SRCSOL(JS,IS,KS)*DCYT2
cea   SRCDCY=SRCDCY+(SRC1-SRCSOL(JS,IS,KS))
c srcdcy now stores mass decayed instead of mass rate
      SRCDCY=SRCDCY+(SRC1-SRCSOL(JS,IS,KS))*TIMV
C.....Apply the decayed change in source concentration
c**** tovera should be just delt, prbinv should be (por*Rf)^-1
c *** prbinv no longer needed, just delt and thickness multiply sources
   40 RS(M)=TOVERA*(SRC1+ETFLO*ETCONC+SNKFLO(JS,IS,KS)*C2
     *  -DIV*C1)
c..   40 RS(M)=TOVERA*PRBINV*(SRC1+ETFLO*ETCONC+SNKFLO(JS,IS,KS)*C2
c..     *  -DIV*C1)
cold      *  -FDTMTH*DIV*CONC(JS,IS,KS)-(1.0-FDTMTH)*DIV*C1)
C
      IF(IIMOV.EQ.1) THEN
        MA=MRNO(M)
c orig line
c...        VAD(MA)=VAD(MA)+FDTMTH*TOVERA*PRBINV*DIV
C THIS WAS CHANGED SO SINKS ARE NOT INCLUDED IN VAD (SEE THETA-W
C TERM IN EQ. 12 OF WRIR 98-4234)  ALSO, FOR CELLS WITH BOTH SRC
C AND SNK, SRCFLO SHOULD GIVES BETTER RESULTS
c...       VAD(MA)=VAD(MA)+FDTMTH*TOVERA*PRBINV*SRCFLO(JS,IS,KS)
       VAD(MA)=VAD(MA)+FDTMTH*TOVERA*SRCFLO(JS,IS,KS)
c       vad(ma)=vad(ma)+fdtmth*timv*thck(js,is,ks)*srcflo(js,is,ks)
c**** same change as above
      END IF
   65 CONTINUE
C
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
C
C  ELLSRC  CHANGE IN CONC DUE TO SINK/SOURCES
C********************************************************
C
      SUBROUTINE ELLSRC(IBOUND,SRCFLO,SRCSOL,SNKFLO,RF,
     *  THCK,POR,DIAGDS,RHS,CNOLD,
     *  EVTFLO,IEVT,SRCDCY,NODESS,IMOV,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NEVTOP,
     *  IOUTS,TIMV,DECAY,
     *  VC,VR,VL,VCMAX,VRMAX,VLMAX,
     *  LBNDY,NTFACE,RHSO,
     *  DELCOL,DELROW,XFOR,XBAC,YFOR,YBAC,AGER8,SAGE,
     *  INAGE)
C
C********************************************************
C
      DOUBLE PRECISION DECAY
C
      DIMENSION SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *  THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),IBOUND(NCOL,NROW,NLAY),
     *  EVTFLO(NSCOL,NSROW,2),IEVT(NROW,NCOL),DIAGDS(27,NODESS),
     *  RHS(NODESS),CNOLD(NODESS)
C
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
C     ***************************************************************
C
C3D  DECAY AND SOURCES FOR ONE STEP ONLY
C
        NSRC=NSROW*NSCOL
        DO 131 KS=1,NSLAY
           K=ISLAY1+KS-1
           DO 131 IS=1,NSROW
              I=ISROW1+IS-1
              DO 131 JS=1,NSCOL
                 J=ISCOL1+JS-1
                 NODE=(KS-1)*NSRC+(IS-1)*NSCOL+JS
C
      IF (IBOUND(J,I,K).EQ.0) GOTO 131 
C
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
C           MIXING AT SOURCE CELLS...
C           ...WITH DECAY OF RECHARGE DURING TIME INCREMENT
C
C   SINK CELL PROCESSING
      KEV=1
      ETFLO=0.0
      KEV=EVTFLO(JS,IS,2)
      IF(KEV.EQ.K) ETFLO=EVTFLO(JS,IS,1)
      IF ((SNKFLO(JS,IS,KS).GE.0D0).AND.(ETFLO.EQ.0D0)) GOTO 100
C
      TORF=TIMV/RF(JS,IS,KS)
      PRINV=1.D0/POR(JS,IS,KS)
      ETCONC=0.0
cea
      IF(INAGE.GT.0) ETCONC=CNOLD(NODE)
C
      IF (IMOV.EQ.1) DIAGDS(14,NODE)=DIAGDS(14,NODE)-
     *          TORF*0.5D0*SNKFLO(JS,IS,KS)*PRINV
C
      RHS(NODE)=RHS(NODE)+TORF*
     *         (0.5D0*SNKFLO(JS,IS,KS)*CNOLD(NODE)+ETFLO*ETCONC)*PRINV
C
100   CONTINUE
C
C  SOURCE CELL PROCESSING
C  SRCSOL COULD BE SCALED AFTER VELO
C  CHANGE IN SIGN CONVENTION FOR FLUID SOURCES FROM MOC-2D
C  HERE SRCFLO IS POSITIVE FOR INFLOW TO AQUIFER
C
      IF (SRCSOL(JS,IS,KS).LE.0D0) GOTO 131
      SRC1=SRCSOL(JS,IS,KS)/RF(JS,IS,KS)
C
      AREA=CINV*RINV*BINV
      YHAT=-.5
      ZHAT=-.5-HBINV
      DO 800 NYL=1,NSL
      ZHAT=ZHAT+BINV
C
      YHAT=-.5-HRINV
      DO 800 NXR=1,NSR
      YHAT =YHAT+RINV
C
      XHAT=-.5-HCINV
      DO 800 NYC=1,NSC
      XHAT=XHAT+CINV
C
cea
      AGEVOL=SRCFLO(JS,IS,KS)
C
      CALL BDYRHS(VC,VR,VL,RF,
     *           THCK,POR,IBOUND,
     *           NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *           IOUTS,TIMV,VCMAX,VRMAX,VLMAX,
     *           NODESS,LBNDY,NTFACE,SRC1,AREA,
     *           J,I,K,XHAT,YHAT,ZHAT,RHS,RHSO,
     *           DELCOL,DELROW,XFOR,XBAC,YFOR,YBAC,DECAY,SRCDCY,
     *           AGER8,AGEVOL,SAGE)
C
800   CONTINUE
C
131   CONTINUE
C     ****************************************************************
      RETURN
C     ****************************************************************
C
      END
