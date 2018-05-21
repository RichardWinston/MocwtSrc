C
C*************************************************************************
C
C SMOC6MBRP_WT    PRINT SEPARATE MASS BALANCE REPORT: MOCWT
      SUBROUTINE SMOC6MBRP_WT(SBVL,SRCDCY,IPERGWT,
     *  KSTP,KPER,IMOV,NSTP,NMOV,SUMTCH,NIUNIT,INDP,MULTSS,
     *  IOUTMBRPT,IOUTMBITEM,NSOL)
C
C     *******************************************************
C     PRINT "SPREADSHEET" OF MASS BALANCE FIGURES: WEIGHT PARTICLES OPTION
C     *******************************************************
C
      DOUBLE PRECISION SBVL
      DIMENSION SBVL(6,NIUNIT)
      DOUBLE PRECISION CSTM2,TOTMIN,TOTMOT,RESID,TIN,TOUT,
     *  DENOM,ERR,SMIN,SMPR,ERR2,TOTCUMIN,TOTCUMOT,ERRMIN
C
C*************************************************************************
C
C  LOOP OVER ALL SOLUTES 
C when this is implemented, the MBRP unit numbers will be in an array over NSOL
C   (as will be SBVL)
      DO 1000 ISOL=1,NSOL
C  ZERO ACCUMULATORS
      TOTMIN=0.0
      TOTMOT=0.0
      TOTCUMIN=0.0
      TOTCUMOT=0.0
C  ADD DECAY TERMS
      TOTMIN=TOTMIN+SBVL(3,3)
      TOTMOT=TOTMOT+SBVL(4,3)
C  LOOP THROUGH PACKAGES, ADDING APPROPRIATE MASS FLUXES
      DO 50 L=4,12
         TOTMIN=TOTMIN+SBVL(1,L)
         TOTMOT=TOTMOT+SBVL(4,L)
 50      CONTINUE
C  SIMPLE REACTIONS SINK/SOURCES AND AGE MASS INCREASE
      DO 55 L=13,17
        TOTMIN=TOTMIN+SBVL(1,L)
        TOTMOT=TOTMOT+SBVL(2,L)
 55   CONTINUE
C  DOUBLE POROSITY DIFFUSION SINK/SOURCES
      DO 56 L=19,20
        TOTMIN=TOTMIN+SBVL(1,L)
        TOTMOT=TOTMOT+SBVL(2,L)
 56   CONTINUE
C  Lake and Stream Packages
      DO 58 L=21,22
        TOTMIN=TOTMIN+SBVL(1,L)
        TOTMOT=TOTMOT+SBVL(2,L)
 58   CONTINUE
C  DRT Package
      DO 60 L=23,23
        TOTMIN=TOTMIN+SBVL(1,L)
        TOTMOT=TOTMOT+SBVL(2,L)
 60   CONTINUE
C  MNW Package
      DO 61 L=24,24
        TOTMIN=TOTMIN+SBVL(1,L)
        TOTMOT=TOTMOT+SBVL(4,L)
 61   CONTINUE
C  CCBD Package
      DO 62 L=26,26
        TOTMIN=TOTMIN+SBVL(1,L)
        TOTMOT=TOTMOT+SBVL(2,L)
 62   CONTINUE
C  CALCULATE CHANGE IN MASS STORED
      CSTM2=SBVL(4,1)-SBVL(1,1)+SBVL(2,2)-SBVL(1,2)
      IF(INDP.GT.0) CSTM2=CSTM2-SBVL(1,18)+SBVL(2,18)
C
C     CUMULATIVE BUDGET NUMBERS FOR CASE OF MULTIPLE STEADY STRESS PERIODS
      IF(MULTSS.EQ.1) THEN
C     CALC TEMPORARY CUMULATIVE TOTALS
       IF(KSTP.LT.NSTP.OR.IMOV.LT.NMOV) THEN
        TOTCUMIN=TOTCUMIN+SBVL(5,3)+SBVL(3,3)
        TOTCUMOT=TOTCUMOT+SBVL(6,3)+SBVL(4,3)
        DO 8000 I=4,12
          TOTCUMIN=TOTCUMIN+SBVL(1,I)+SBVL(5,I)
          TOTCUMOT=TOTCUMOT+SBVL(4,I)+SBVL(6,I)
 8000  CONTINUE
        DO 8001 I=13,23
          TOTCUMIN=TOTCUMIN+SBVL(1,I)+SBVL(5,I)
          TOTCUMOT=TOTCUMOT+SBVL(2,I)+SBVL(6,I)
 8001  CONTINUE
        DO 8002 I=24,24
          TOTCUMIN=TOTCUMIN+SBVL(1,I)+SBVL(5,I)
          TOTCUMOT=TOTCUMOT+SBVL(4,I)+SBVL(6,I)
 8002  CONTINUE
        DO 8003 I=26,26
          TOTCUMIN=TOTCUMIN+SBVL(1,I)+SBVL(5,I)
          TOTCUMOT=TOTCUMOT+SBVL(2,I)+SBVL(6,I)
 8003  CONTINUE
       ELSE
        TOTCUMIN=TOTCUMIN+SBVL(5,3)
        TOTCUMOT=TOTCUMOT+SBVL(6,3)
        DO 8004 I=3,24
          TOTCUMIN=TOTCUMIN+SBVL(5,I)
          TOTCUMOT=TOTCUMOT+SBVL(6,I)
 8004  CONTINUE
        DO 8005 I=26,26
          TOTCUMIN=TOTCUMIN+SBVL(5,I)
          TOTCUMOT=TOTCUMOT+SBVL(6,I)
 8005  CONTINUE
       END IF
      END IF
C  CALCULATE RESIDUAL
c initialize ERR in case of no solute at all in system
      ERR=0.0
      RESID=TOTMIN+TOTMOT+SRCDCY-CSTM2
c for error calc, take decayed source out of total mass in
      TIN=ABS(TOTMIN)+SRCDCY
      TOUT=ABS(TOTMOT)
C  ZERO MASS BALANCE FORMULA FLAGS
      IMBF=0
      IMBF2=0
      IPREC=0
C  CHECK FLAGS FOR MASS FLUXES, CALCULATE WITH ACCORDING
C  FORMULA
      IF( (SBVL(1,1).NE.0.0).OR.
     *   (TIN.EQ.0.0.AND.TOUT.EQ.0.0) ) GOTO 97
cgzh mbprec
c check CSTM2 vs present mass, if result is less than precision,
c don't use the flux in/out error as it will not be close enough       
      PRECSM=SBVL(4,1)*1E-7
	IF(ABS(CSTM2).LT.PRECSM) THEN
        IPREC=1
        GO TO 97
      END IF
      IF(TIN.GT.TOUT) THEN
         DENOM=TOTMIN
         IMBF=1
      ELSE
         DENOM=TOTMOT
         IMBF=2
      ENDIF
      ERR=(RESID*100.0)/DENOM
  97  IF(SBVL(1,1).EQ.0.) GOTO 95
C  SMIN=STORED MASS INITIALLY
      SMIN=SBVL(1,1)+SBVL(1,2)
C  SMPR=STORED MASS PRESENT
      SMPR=SBVL(4,1)+SBVL(2,2)
      IF(INDP.GT.0) THEN
        SMIN=SMIN+SBVL(1,17)
        SMPR=SMPR+SBVL(2,17)
      END IF
      IF(IMBF.EQ.1) THEN
C  RATIO MAY BE ADJUSTED
         IF((TIN/SMIN).LT.0.50) THEN
            DENOM=SMPR
            IMBF2=1
         ENDIF
      ENDIF
      IF(IMBF.EQ.2) THEN
C  RATIO MAY BE ADJUSTED
         IF((TOUT/SMIN).LT.0.50) THEN
            DENOM=SMPR
            IMBF2=1
         ENDIF
      ENDIF    
      IF(TIN.EQ.0.0.AND.TOUT.EQ.0.0) THEN
         DENOM=SMPR
         IMBF2=1
      ENDIF
c    mbprec use stored mass error calculation when precision concern with 
c    flux in/out
      IF(IPREC.EQ.1) THEN
         DENOM=SMPR
	   IMBF2=1
      END IF
  95  CONTINUE
C
	IF(IMBF2.GT.0) THEN
         ERR2=(RESID*100.0)/DENOM
	ENDIF 
C PRINT OUT LOWEST ERROR
      ERRMIN=MIN(ERR,ERR2)

C
C  WRITE HEADERS
      IF(KPER.EQ.IPERGWT.AND.KSTP.EQ.1.AND.IMOV.EQ.1) THEN
       IF(IOUTMBRPT.GT.0) THEN
	  WRITE (IOUTMBRPT,100) NSOL
        WRITE (IOUTMBRPT,105)
       END IF
       IF(IOUTMBITEM.GT.0) THEN
	  WRITE (IOUTMBITEM,100) NSOL
        WRITE (IOUTMBITEM,108)
       END IF
      END IF
C
C  WRITE DATA
      IF(IOUTMBRPT.GT.0) THEN
       WRITE (IOUTMBRPT,110) KPER,KSTP,IMOV,SUMTCH,
c     initdiss, presdiss, initsorb, pressorb
     & SBVL(1,1),SBVL(4,1),SBVL(1,2),SBVL(2,2),
c     dpinit,dppres
     & SBVL(1,18),SBVL(2,18),
c     chnge stored mass
     & CSTM2,
c     total mass in, out, srcdecay
     & TOTMIN,TOTMOT,SRCDCY,
c     total cumulative mass in, out, 
     & TOTCUMIN,TOTCUMOT,
c     resid, error
     & RESID,ERRMIN
      END IF
      IF(IOUTMBITEM.GT.0) THEN
       WRITE (IOUTMBITEM,112) KPER,KSTP,IMOV,SUMTCH,
c     initdiss, presdiss, initsorb, pressorb
     & SBVL(1,1),SBVL(4,1),SBVL(1,2),SBVL(2,2),
c     dpinit,dppres
     & SBVL(1,18),SBVL(2,18),
c     chnge stored mass
     & CSTM2,
c     total mass in, out, srcdecay
     & TOTMIN,TOTMOT,SRCDCY,
c     total cumulative mass in, out, 
     & TOTCUMIN,TOTCUMOT,
c     resid, error
     & RESID,ERRMIN,
c     decay, ch_in, ch_out,subgridin, subgridout
     &SBVL(3,4),SBVL(1,4),SBVL(4,4),SBVL(1,5),SBVL(4,5),
c     rch in, rch out
     &SBVL(1,6),SBVL(4,6),
c     wel in, wel out
     &SBVL(1,7),SBVL(4,7),
c     riv in, riv out
     &SBVL(1,8),SBVL(4,8),
c     drn in, drn out
     &SBVL(1,9),SBVL(4,9),
c     ghb in, ghb out
     &SBVL(1,10),SBVL(4,10),
c     evt 
     &SBVL(2,11),
c     fhb in, fhb out
     &SBVL(1,12),SBVL(4,12),
c     zero order growth/decay
     &SBVL(1,13),SBVL(2,13),
c     first order growth/decay
     &SBVL(1,14),SBVL(2,14),
c     zero order growth/decay, sorbed phase
     &SBVL(1,15),SBVL(2,15),
c     first order growth/decay, sorbed phase
     &SBVL(1,16),SBVL(2,16),
c     age mass increase, dp 0-order growth, dp 1st-order decay
     &SBVL(1,17),SBVL(1,19),SBVL(2,20),
c     sfr in, sfr out
     &SBVL(1,21),SBVL(2,21),
c     lak in, lak out
     &SBVL(1,22),SBVL(2,22),
c     drt in, drt out
     &SBVL(1,23),SBVL(4,23),
c     mnw in, mnw out, mnw borehole
     &SBVL(1,24),SBVL(4,24),SBVL(3,24),
c     mnw inQNET, mnw outQNET
     &SBVL(1,25),SBVL(2,25),
c     ccbd in, ccbd out
     &SBVL(1,26),SBVL(2,26)

      END IF
 1000 CONTINUE
C
 100  FORMAT ('  MASS BALANCE REPORT (MBRP) FILE FOR SOLUTE #',I3,'.  
     & SEE DOCUMENTATION FOR EXPLANATION OF ABBREVIATIONS.')
 105  FORMAT ('  StressPed FlowTS TransportTS Tot.TIME DISS_INIT
     & DISS_PRES
     & SORB_INIT SORB_PRES DP_INIT DP_PRES DELTA_MASS MASS_IN MASS_OUT
     & SRC_DECAY CUMU_MASS_IN CUMU_MASS_OUT RESIDUAL %_ERROR')
 108  FORMAT ('  StressPed FlowTS TransportTS Tot.TIME DISS_INIT
     & DISS_PRES
     & SORB_INIT SORB_PRES DP_INIT DP_PRES DELTA_MASS MASS_IN MASS_OUT
     & SRC_DECAY CUMU_MASS_IN CUMU_MASS_OUT RESIDUAL %_ERROR DECAY
     & CH_IN CH_OUT SUBGRID_IN SUBGRID_OUT RECH_IN RECH_OUT WEL_IN
     & WEL_OUT RIV_IN RIV_OUT DRN_IN DRN_OUT GHB_IN GHB_OUT EVT
     & FHB_IN FHB_OUT ZERO_GROWTH ZERO_DECAY FIRST_GROWTH FIRST_DECAY
     & ZERO_GROW_SORB ZERO_DECAY_SORB FIRST_GROW_SORB FIRST_DECAY_SORB
     & AGE_MASS DP_ZERO_GROWTH DP_FIRST_DECAY SFR_IN SFR_OUT LAK_IN
     & LAK_OUT DRT_IN DRT_OUT MNW_IN MNW_OUT MNW_BH MNW_IN_NET
     & MNW_OUT_NET CCBD_IN CCBD_OUT')
 110  FORMAT (3I6,1P16E18.5)
 112  FORMAT (3I6,1P58E18.5)
C
      RETURN
      END
C