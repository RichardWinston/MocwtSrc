C     Last change: RBW Dec. 8, 2014
C     Support for weighted particles (MOCWT) added.
C  MOCMB  GWT MASS BALANCE
C
C  SMOC6IM  CALCULATE INITIAL MASS (MOC AND MOCIMP)
C*************************************************************************
C
      SUBROUTINE SMOC6IM(IBOUND,CONC,THCK,RF,POR,SBVL,
     *  NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,
     *  IOUTS,JRF,NIUNIT)
C*************************************************************************
C
C
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *  CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  RF(NSCOL,NSROW,NSLAY),POR(NSCOL,NSROW,NSLAY),
     *  SBVL(6,NIUNIT)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C*************************************************************************
C
C  CALCULATE INITIAL SOLUTE MASS STORED FOR CHEMICAL MASS BALANCE
C
C********************************************************************
      AREA=CDEL*RDEL
C  CALCULATE INITIAL DISSOLVED MASS IN STORAGE
      DO 973 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 973 IS=1,NSROW
      I=IS+ISROW1-1
      DO 973 JS=1,NSCOL
      J=JS+ISCOL1-1
      IF(IBOUND(J,I,K).EQ.0) GO TO 973
      SBVL(1,1)=SBVL(1,1)+
     1  CONC(JS,IS,KS)*THCK(JS,IS,KS)*POR(JS,IS,KS)*AREA
  973 CONTINUE
C
C  CALCULATE INITIAL SOLUTE MASS SORBED ON SOLIDS
      IF (JRF.EQ.1) THEN
      DO 975 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 974 IS=1,NSROW
      I=IS+ISROW1-1
      DO 974 JS=1,NSCOL
      J=JS+ISCOL1-1
      IF(IBOUND(J,I,K).EQ.0) GO TO 974
      SBVL(1,2)=SBVL(1,2)+
     1  CONC(JS,IS,KS)*THCK(JS,IS,KS)*POR(JS,IS,KS)*AREA*
     *    (RF(JS,IS,KS)-1.0)
 974  CONTINUE
 975  CONTINUE
      END IF
C
C*************************************************************************
      RETURN
      END
C
C
C SMOC5I CALCULATE INITIAL STORED AND ADSORBED MASS (ELLAM)
C***********************************************************************
C
      SUBROUTINE SMOC5I(JRF,NODES,NODESS,NSCOL,NSROW,NSLAY,
     *          NCOL,NROW,NLAY,NZIN,TEMP,JAS,IAS,AS,NELTS,
     *          RW,RW1,LENW,leniw,KKSTP,KKPER,
     *          CNOLD,RF,DIAGST,POR,IBOUND,SBVL,NRFLG,THCK,NIUNIT,
     *          INAGE,AGER8,DELT,IPERGWT)
C***********************************************************************
C  SMOC5I 
C  - CONVERTS STORAGE MATRIX, DIAGST, TO SLAP COLUMN FORMAT ACCEPTED BY
C  SOLVER
C  - CALLS SMOC5A TO CALCULATE INITIAL MASS IN TRANSPORT DOMAIN USING
C  INITIAL NODAL CONCENTRATIONS STORED IN CNOLD AND INITIAL BOUNDARY
C  MASS STORED IN STINIT, ADINIT.  TOTAL STORED AND ADSORBED MASS IS
C  RETURNED IN STINIT AND ADINIT.
C  - UPDATES SBVL WITH STORED AND ADSORBED MASS, AND OLMASS WITH TOTAL 
C  MASS, FOR USE BY MOCBD.  OLMASS IS USED TO CALCULATE DECAY.
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
      DIMENSION SBVL(6,NIUNIT),CNOLD(NODESS),DIAGST(27,NODESS),RW(LENW)
      DIMENSION NZIN(NODESS),TEMP(LENW),JAS(1+NODESS),IAS(27*NODESS),
     *          AS(27*NODESS),RF(NSCOL,NSROW,NSLAY),POR(NODESS),
     *          IBOUND(NCOL,NROW,NLAY),IDISP(-1:1,-4:4)
cea
      DIMENSION RW1(NODESS)
C
C
      CALL GETIDISP(NSCOL,NSROW,NSLAY,  IDISP)
      IL=1
      IF (NSLAY.EQ.1) IL=0
      IR=1
      IF (NSROW.EQ.1) IR=0
      IC=1
      IF (NSCOL.EQ.1) IC=0
      NE=1
      DO NN=1,NODESS
        JAS(NN)=NE
        IAS(NE)=NN
        RW(NE)=DIAGST(14,NN)
        IF (RW(NE).EQ.0D0) RW(NE)=1D0
        NE=NE+1
        DO NL=-IL,IL
          DO NR=-IR,IR
            DO NC=-IC,IC
              IF (NL.NE.0 .OR. NR.NE.0 .OR. NC.NE.0) THEN
                NNROW=NN+IDISP(NC,3*NL+NR)
                IF (NNROW.GT.0 .AND. NNROW.LE.NODESS) THEN
                  NP=14-3*(3*NL+NR)-NC
                  T= DIAGST(NP,NNROW)
                  IF (T.NE.0) THEN
                    IAS(NE)=NNROW
                    RW(NE)=T
                    NE=NE+1
                  ENDIF
                ENDIF
              ENDIF
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      NELTS=NE-1
      JAS(NODESS+1)=NE
      DO NN=1,NELTS
        AS(NN)=RW(NN)
      ENDDO
C
C  IF FIRST STRESS PERIOD AND FIRST FLOW TIME STEP,
C  THEN SAVE INITIAL MASS
C
      IF(NRFLG.EQ.1) RETURN
      CALL SMOC5A(JRF,NODES,NODESS,NSCOL,NSROW,NSLAY,
     *          NCOL,NROW,NLAY,NZIN,TEMP,JAS,IAS,AS,NELTS,
     *          CNOLD,RF,POR,IBOUND,STINIT,ADINIT,LENW,THCK)
C
      IF (KKPER.EQ.IPERGWT.AND.KKSTP.EQ.1) THEN
         SBVL(1,1)=STINIT
         SBVL(1,2)=ADINIT
      ENDIF
      OLMASS   =STINIT   +ADINIT
C
cea
C  FIND WATER VOLUME THIS FLOW TIME PERIOD FOR AGE
C  USE SOLVER WORK AREA FOR VECTOR OF ONES
C  SMOC5A ACCUMULATES BOUNDARY + INTERIOR
C  INTO WATVOL
         IF(INAGE.GT.0) THEN
            DO NN=1,NODESS
             RW1(NN)=1.0
            ENDDO
            CALL SMOC5A(JRF,NODES,NODESS,NSCOL,NSROW,NSLAY,
     *          NCOL,NROW,NLAY,NZIN,TEMP,JAS,IAS,AS,NELTS,
     *          RW1,RF,POR,IBOUND,WATVOL,TMP,LENW,THCK)
         ENDIF
      RETURN
      END
C 
C
C SMOC5A  ACCUMULATE STORED AND ADSORBED MASS (ELLAM)
C**********************************************************************
C
      SUBROUTINE SMOC5A(JRF,NODES,NODESS,NSCOL,NSROW,NSLAY,
     *          NCOL,NROW,NLAY,NZIN,TEMP,JAS,IAS,AS,NELTS,
     *          C,RF,POR,IBOUND,STORED,SORBED,LENW,THCK)
C
C  STORED,SORBED HOLD MASS VALUES FROM BOUNDARY NODES UPON ENTRY
C
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      DIMENSION C(NODESS),RF(NSCOL,NSROW,NSLAY),
     *          POR(NODESS),IBOUND(NCOL,NROW,NLAY),
     *          TEMP(LENW),NZIN(NODESS),THCK(NODESS)
C
C   MULTIPLY AS*C 
C      WHERE C IS INITIAL CONCENTRATIONS STORED IN CNOLD
C      IF SMOC5A IS CALLED BY SMOC5I
C      AND C IS CURRENT CONCENTRATION STORED IN CONC
C      IF SMOC5A IS CALLED BY MOCBD
C   RESULT IS MASS/POROSITY
C
C   USE SOLVER WORK AREA FOR PRODUCT
      ISYM=0
      CALL SSMV(NODESS, C, TEMP, NELTS, IAS, JAS, AS, ISYM )
C
C   CONVERT RESULT TO MASS IN CELL BY DIVISION BY CELL POSOSITY
C   STORE CNOFLO IN INACTIVE CELLS FOR REPORTING
C   NOTE CELLS WITH ZERO TOTAL MASS IN NZIN
      TOT=0D0
      N=0
      DO KS=1,NSLAY
        K=ISLAY1+KS-1
        DO IS=1,NSROW
          I=ISROW1+IS-1
          DO JS=1,NSCOL
            J=ISCOL1+JS-1
            N=N+1
            IF (IBOUND(J,I,K).NE.0) THEN
C
C  IF CELL IS ACTIVE, FIND MASS IN CELL, AND
C  ACCUMULATE TOTAL FOR MASS BALANCE
C
              TEMP(N)=TEMP(N)*POR(N)
              TOT=TOT+TEMP(N)
C
C  DETERMINE IF CELL HOLDS MASS, FOR USE BY MASS TRACKING
C  RHS INTEGRATION WON'T TREAT CELLS WITH ZERO MASS
C  AZER0=NUMERICAL ZERO, DETERMINED BY MAX NORM OF RHS AND
C  SET IN ELSOLV (MAY BE SET TO 0D0)
C  IF SMOC5A IS CALLED BY SMOC5I, NZIN INFORMATION IS NEVER 
C  USED, SINCE MASS HAS ALREADY BEEN TRACKED FOR FIRST TIME
C  STEP.
C
cea              IF (ABS(TEMP(N)).GT.AZERO) THEN
cea                NZIN(N)=1
cea              ELSE
cea                NZIN(N)=0
cea              ENDIF
            ELSE
              C(N)=CNOFLO
              TOT=TOT+THCK(N)
	    ENDIF
          ENDDO
        ENDDO
      ENDDO
C
      STORED=STORED+TOT
C
C  IF MASS IS ADSORBED, CALCULATE HOW MUCH      
      IF (JRF.EQ.1) THEN
         TOT=0D0
	 NSCR=NSCOL*NSROW
         N=0
         DO KS=1,NSLAY
           K=ISLAY1+KS-1
           DO IS=1,NSROW
             I=ISROW1+IS-1
             DO JS=1,NSCOL
               J=ISCOL1+JS-1
               N=N+1
               IF (IBOUND(J,I,K).NE.0) THEN
                  TOT=TOT+TEMP(N)*(RF(JS,IS,KS) - 1)
               ENDIF
             ENDDO
           ENDDO
         ENDDO
	 SORBED=SORBED+TOT
      ENDIF
C
cea
C  FIND WATER VOLUME THIS TIME PERIOD 
C  USE SOLVER WORK AREA FOR VECTOR OF ONES
C  CALCULATE AVE CONC FOR OUTPUT
         DO NN=1,NODESS
           TEMP(NN+NODESS)=1.0
         ENDDO
cea multiply conc*ones for 1/por * integral of 1*dvol = cell vol/por
            CALL SSMV(NODESS, TEMP(1+NODESS), TEMP(2*NODESS +1), 
     *       NELTS, IAS, JAS, AS, ISYM )
      DO NN=1,NODESS
cea cell vol:
           TEMP(NN+2*NODESS)=TEMP(NN+2*NODESS)*POR(NN)
      ENDDO

C
C
      RETURN
      END
C
C SMOC6BD   CALCULATE SOLUTE MASS BALANCE
C***********************************************************************
C
      SUBROUTINE SMOC6BD(IBOUND,RF,
     *   CONC,CNOLD,CFXHBC,CINFL,
     *   THCK,POR,CHDFLO,CTCFLO,LOCBDY,
     *   RECH,CRECH,IRCH,NRCHOP,
     *   WELL,MXWELL,NWELLS,NWELVL,IWELLC,
     *   RIVR,MXRIVR,NRIVER,NRIVVL,IRIVRC,
     *   DRAI,MXDRN,NDRAIN,NDRNVL,
     *   BNDS,MXBND,NBOUND,NGHBVL,IBNDSC,
     *   DRTF,MXDRT,NDRTCL,NDRTVL,IDRTFL,
     *   EVTFLO,IEVT,NEVTOP,
     *   IFLLOC,BDFV,NFLW,IFHBD4,IFHBFC,
     *   IHDLOC,BDHV,NHED,NFHBX2,IFHBHC,
     *   DECAY,TIMV,CSTIN,CSTOUT,ALKIN,ALKOUT,NSOL,
     *   NCOL,NROW,NLAY,
     *   NSCOL,NSROW,NSLAY,
     *   IUNIT,IOUTS,NCINFL,JRF,SBVL,
     *   NFACES,IABOVE,IBELOW,MOCTYPE,
     *   DKFO,DKFS,INDK,IDKFO,IDKFS,NIUNIT,INAGE,
CMOCWT
     *   SGMSOUT,SNKMSOUT,SNKFLO,SUMSGMS,SUMBFMS,SUMMASS,
cgzh cbdy
     *   CINFLA,CINFLB,CINXY,
cgzh mnw
     *   nwell2,mxwel2,well2,MNWid)

C
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL,SUMMASS
      DOUBLE PRECISION well2
cgzh debug double sgmsout + snkmsout
      DOUBLE PRECISION SGMSOUT,SNKMSOUT,SUMSGMS,SUMBFMS
      DIMENSION SBVL(6,NIUNIT),IUNIT(NIUNIT),
     *  WELL(NWELVL,MXWELL),RIVR(NRIVVL,MXRIVR),
     *  DRAI(NDRNVL,MXDRN),BNDS(NGHBVL,MXBND),
     *  RECH(NCOL,NROW),CRECH(NSCOL,NSROW),IRCH(NCOL,NROW),
     *  EVTFLO(NSCOL,NSROW,2),IEVT(NCOL,NROW),
     *  DRTF(NDRTVL,MXDRT)
      DIMENSION BDFV(IFHBD4,NFLW),IFLLOC(4,NFLW)
      DIMENSION BDHV(NFHBX2,NHED),IHDLOC(4,NHED)
      DIMENSION well2(18,mxwel2),MNWid(mxwel2)
C
      DIMENSION IBOUND(NCOL,NROW,NLAY),RF(NSCOL,NSROW,NSLAY),
     *  CHDFLO(NSCOL,NSROW,NSLAY),CTCFLO(NFACES),LOCBDY(3,NFACES),
     *  CONC(NSCOL,NSROW,NSLAY),CNOLD(NSCOL,NSROW,NSLAY),
     *  CINFL(NCINFL),CFXHBC(NSCOL,NSROW,NSLAY),
     *  THCK(NSCOL,NSROW,NSLAY),POR(NSCOL,NSROW,NSLAY),
CMOCWT
     *  SGMSOUT(NSCOL,NSROW,NSLAY),SNKMSOUT(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),SUMSGMS(NSCOL,NSROW,NSLAY),
     *  SUMBFMS(NSCOL,NSROW,NSLAY),SUMMASS(NSCOL,NSROW,NSLAY)
C
cgzh cbdy
      DIMENSION CINFLA(NSCOL,NSROW),CINFLB(NSCOL,NSROW),
     * CINXY(NSCOL,NSROW,NSLAY)
      DIMENSION DKFO(NSCOL,NSROW,NSLAY),DKFS(NSCOL,NSROW,NSLAY)
C
      DIMENSION ALKIN(NSOL),ALKOUT(NSOL)
      DIMENSION CSTIN(NSOL),CSTOUT(NSOL)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      DOUBLE PRECISION DECAY,DCYFCT,DCYT,DCYT2
cgzh debug double mass balance numbers
      DOUBLE PRECISION TOTOUT,TOTIN,TOTOUTWT,STORM,ADSORB,ADSORBMASS,
     *  DMASS,CBE,DMASS2,DMASS1
CMOCWT
      INCLUDE 'ptwt.inc'
C
C     ****************************************************************   
C
C  ---COMPUTE MASS BALANCE FOR SOLUTE---                              
C
C     ****************************************************************   
      TOTOUT=0.0 
      TOTIN=0.0
CMOCWT
      TOTOUTWT=0.0
      STORM=0.0                            
      ADSORB=0.0                                                    
      ADSORBMASS=0.0
      DMASS=0.0 
C  COMPUTE DECAY TERM 
      DCYFCT=TIMV*DECAY
      DCYT=1.D0
      DCYT2=1.D0
      IF (DECAY.NE.0.D0) THEN
         DCYT=DEXP(-DCYFCT)
         DCYT2=DEXP(-DCYFCT*0.5D0)
      END IF
C
      AREA=CDEL*RDEL
C
C  LOOP OVER CELLS FOR MASS IN STORAGE, MASS ADSORBED, MASS
C  DECAYED, AND CONSTANT HEAD NODES
      DO 270 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 270 IS=1,NSROW
      I=IS+ISROW1-1
      DO 268 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF NO FLOW NODE
         ICH=IBOUND(J,I,K)
         IF (ICH.EQ.0) GO TO 268
        CBE=CONC(JS,IS,KS)*THCK(JS,IS,KS)*POR(JS,IS,KS)
C
C  MASS OF SOLUTE IN STORAGE
C
        STORM=STORM+CBE    
C
C  COMPUTE MASS ADSORBED
C
C  IF RF.NE.1
        IF (JRF.EQ.1) THEN
          IF(PTWTON.EQ.0) THEN
            ADSORB=ADSORB+CBE*(RF(JS,IS,KS)-1.0)
CMOCWT  FOR WEIGHTED OPTION, USE MASS ON PTS FOR THIS CALC
          ELSE
            ADSORBMASS=ADSORBMASS+SUMMASS(JS,IS,KS)*(RF(JS,IS,KS)-1.0)
          END IF
        END IF
C
C  MASS DECAYED
C
        IF (CNOLD(JS,IS,KS).GT.0.0) THEN
          IF (DECAY.NE.0.D0) THEN
            DELDCY=CNOLD(JS,IS,KS)*(1.0-DCYT)
            DMASS=DMASS-DELDCY*THCK(JS,IS,KS)*AREA*
     *            POR(JS,IS,KS)*RF(JS,IS,KS)
          ELSE
            IF(INDK.GT.0) THEN
              IF(IDKFO.EQ.1) THEN
                DCYFCT=DKFO(JS,IS,KS)*TIMV
                DCYT=DEXP(-DCYFCT)
                DELDCY=CNOLD(JS,IS,KS)*(1.0-DCYT)
                DMASS1=-DELDCY*THCK(JS,IS,KS)*AREA*POR(JS,IS,KS)
                IF(IDKFS.EQ.1) THEN
                  DCYFCT=DKFS(JS,IS,KS)*TIMV
                  DCYT=DEXP(-DCYFCT)
                  DELDCY=CNOLD(JS,IS,KS)*(1.0-DCYT)
                  DMASS2=-DELDCY*THCK(JS,IS,KS)*AREA*
     *              POR(JS,IS,KS)*(RF(JS,IS,KS)-1.0)
                ELSE
                  DMASS2=-DELDCY*THCK(JS,IS,KS)*AREA*
     *              POR(JS,IS,KS)*(RF(JS,IS,KS)-1.0)
                END IF
              ELSE IF(IDKFS.EQ.1) THEN
                  DCYFCT=DKFS(JS,IS,KS)*TIMV
                  DCYT=DEXP(-DCYFCT)
                  DELDCY=CNOLD(JS,IS,KS)*(1.0-DCYT)
                  DMASS2=-DELDCY*THCK(JS,IS,KS)*AREA*
     *              POR(JS,IS,KS)*(RF(JS,IS,KS)-1.0)
              END IF
              IF (DMASS1.LT.0.0) SBVL(2,14)=SBVL(2,14)+DMASS1
              IF (DMASS1.GT.0.0) SBVL(1,14)=SBVL(1,14)+DMASS1
              IF (DMASS2.LT.0.0) SBVL(2,16)=SBVL(2,16)+DMASS2
              IF (DMASS2.GT.0.0) SBVL(1,16)=SBVL(1,16)+DMASS2
            END IF
          END IF
        END IF
C
C  CONSTANT-HEAD NODES
C
        IF (CHDFLO(JS,IS,KS).LT.0.) THEN 
          IF(MOCTYPE.EQ.1) THEN
            TOTOUT=TOTOUT+CHDFLO(JS,IS,KS)*CNOLD(JS,IS,KS)
          ELSE
            TOTOUT=TOTOUT+(CHDFLO(JS,IS,KS)*(CNOLD(JS,IS,KS)+
     *                      CONC(JS,IS,KS))*0.5)
          END IF
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT CONSTANT HEAD = RATIO OF CHDFLO TO TOTAL SNKFLO * TOTAL MASS OUT OF
C   SINKS IN CELL (SNKMSOUT) (CALCULATED ON PARTICLES IN MOVEWT)
          IF(PTWTON.EQ.1) THEN
            TOTOUTWT=TOTOUTWT+(CHDFLO(JS,IS,KS)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
	    END IF
        ELSE     
          TOTIN=TOTIN+CHDFLO(JS,IS,KS)*CFXHBC(JS,IS,KS)
        ENDIF
  268 CONTINUE        
  270 CONTINUE        
C
C  END OF LOOP OVER CELLS
C
C  TOTAL CONSTANT HEAD NODES
      SBVL(1,4)=SBVL(1,4)+TOTIN*TIMV
	SBVL(2,4)=SBVL(2,4)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
      IF (PTWTON.EQ.1) THEN 
	  SBVL(4,4)=SBVL(4,4)+TOTOUTWT
        TOTOUTWT=0.0
      END IF
      TOTIN=0.0
      TOTOUT=0.0
C
C  TOTAL STORED MASS AND PRESENT MASS SORBED
C
      SBVL(2,1)=STORM*AREA
C  IF RF.NE.1
      IF (JRF.EQ.1) THEN
cgzh bug fix
        IF(PTWTON.EQ.0) THEN
         SBVL(2,2)=ADSORB*AREA  
CMOCWT  FOR WEIGHTED OPTION, USE MASS ON PTS FOR THIS CALC
        ELSE
         SBVL(2,2)=ADSORBMASS  
        END IF
      END IF
C
C  MASS DECAYED AND MASS "CREATED"
C   
C MOCWT FOR WEIGHTED OPTION, MASS DECAY IS TRACKED ON PARTICLES IN MOVEWT ROUTINE
      IF (DMASS.LT.0.0) SBVL(2,3)=SBVL(2,3)+DMASS
      IF (DMASS.GT.0.0) SBVL(1,3)=SBVL(1,3)+DMASS
C
C  CALCULATE SOLUTE FLUX ACROSS SUBGRID BOUNDARIES
C
      NSROWF=0
      NSCOLF=0
      IF(ISROW1.GT.1) NSROWF=1
      IF(ISROW2.LT.NROW) NSROWF=NSROWF+1
      IF(ISCOL1.GT.1) NSCOLF=1
      IF(ISCOL2.LT.NCOL) NSCOLF=NSCOLF+1
      IDTOPM=((NSROWF*NSCOL+NSCOLF*NSROW)*NSLAY)+1
C  GET LOCATION OF EACH SUBGRID BOUNDARY FACE
      DO 287 IE=1,NFACES
         J=LOCBDY(1,IE)
         I=LOCBDY(2,IE)
         K=LOCBDY(3,IE)
C 284-LESS THAN;  285-EQUAL TO;  286-GREATER THAN 0
C LESS THAN 0, CTCFLO POINTS OUT OF AQUIFER 
C GREATER THAN 0, CTCFLO POINTS INTO AQUIFER 
         IF (CTCFLO(IE)) 284,285,286
 284     CONTINUE
           JS=J-ISCOL1+1
           IS=I-ISROW1+1
           KS=K-ISLAY1+1 
           IF(MOCTYPE.EQ.1) THEN
             TOTOUT=TOTOUT+CTCFLO(IE)*CNOLD(JS,IS,KS)
           ELSE
             TOTOUT=TOTOUT+(CTCFLO(IE)*(CNOLD(JS,IS,KS)+
     *                    CONC(JS,IS,KS))*0.5)
           END IF
           GO TO 287
 286     JS=J-ISCOL1+1
         IS=I-ISROW1+1
         KS=K-ISLAY1+1 
         IF (IE.LT.IDTOPM) THEN
cgzh cinxy
corig           CPRIME=CINFL(KS)
            CPRIME=CINXY(JS,IS,KS)
         ELSE
cgzh cbdy
            CPRIME=0.0
            IF (MOD(IE,2).EQ.0) THEN
              CPRIME=CINFLB(JS,IS)  
            ELSE
              CPRIME=CINFLA(JS,IS)  
            END IF
         END IF
         TOTIN=TOTIN+CTCFLO(IE)*CPRIME
         GO TO 287
 285     CONTINUE
 287     CONTINUE
      SBVL(1,5)=SBVL(1,5)+TOTIN*TIMV
      SBVL(2,5)=SBVL(2,5)+TOTOUT*TIMV

CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT SUBGRID BOUNDARY = MASS FROM DEPARTED PARTICLES (SUMSGMS) + 
C   MASS REMOVED ON PARTICLES IN CELL (SGMSOUT) (CALCULATED IN MOVEWT)
C SUMSGMS IS NEGATIVE IF LEAVING SUBGRID
C SGMSOUT IS NEGATIVE IF LEAVING SUBGRID
      IF(PTWTON.EQ.1) THEN
cgzh debug output
c       if(js.eq.5.and.is.eq.5) then
c      write(iouts,*) 'Mass on pts that left 5,5=',SUMSGMS(JS,IS,KS)
c      write(iouts,*)'SG Mass leaving,applied to exisiting pts at 5,5=',
c     *    SGMSOUT(JS,IS,KS)
c	 end if
c       temp1=0.0
c	temp2=0.0
cgzh debug
C  LOOP OVER CELLS FOR MASS OUT SUBGRID, WEIGHTED PARTICLES OPTION
      DO 222 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 222 IS=1,NSROW
      I=IS+ISROW1-1
      DO 222 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF NO FLOW NODE
        IF (IBOUND(J,I,K).EQ.0) GO TO 222 
	  TOTOUTWT=TOTOUTWT+SUMSGMS(JS,IS,KS)+SGMSOUT(JS,IS,KS)
cgzh debug
c        temp1=temp1+SUMSGMS(JS,IS,KS)
c        temp2=temp2+SGMSOUT(JS,IS,KS)
  222 CONTINUE
      END IF
cgzh debug output
c        write(iouts,*) 'TOTOUTWT',TOTOUTWT
c        write(iouts,*) 'out subgrid: leaving points',temp1
c        write(iouts,*) 'out subgrid: extracted mass',temp2
cgzh
CMOCWT   SGMSOUT HAS TIMV IN IT, SUMSGMS IS COMPLETE TOO
      IF (PTWTON.EQ.1) THEN 
        SBVL(4,5)=SBVL(4,5)+TOTOUTWT
        TOTOUTWT=0.0
      END IF
      TOTIN=0.0
      TOTOUT=0.0
C
C  RECHARGE
C
      IF (IUNIT(8).GT.0) THEN
        KR=1
        DO 20 IS=1,NSROW
          I=IS+ISROW1-1
        DO 20 JS=1,NSCOL
          J=JS+ISCOL1-1
          IF(NRCHOP.NE.1) KR=IRCH(J,I)
          IF(IBOUND(J,I,KR).LE.0) GO TO 20
          KS=KR-ISLAY1+1
          RCHRAT=RECH(J,I)
          IF (KS.LT.1.OR.KS.GT.NSLAY) GO TO 20
C IF < 0, OUT OF SUBGRID
          IF(RCHRAT.LT.0.0) THEN
            TOTOUT=TOTOUT+RCHRAT*(CNOLD(JS,IS,KS)+CONC(JS,IS,KS))*0.5
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT RECH = RATIO OF RECH TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
              IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(RCHRAT/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
cgzh debug
c      write(iouts,*) 'Mass out rech=',SNKMSOUT(JS,IS,KS)
              ENDIF
C IF > 0, INTO SUBGRID
          ELSE IF(RCHRAT.GT.0.0) THEN
            TOTIN=TOTIN+RCHRAT*CRECH(JS,IS)
          END IF
   20 CONTINUE
        SBVL(1,6)=SBVL(1,6)+TOTIN*TIMV
        SBVL(2,6)=SBVL(2,6)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,6)=SBVL(4,6)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
        TOTIN=0.0
        TOTOUT=0.0
      END IF
C        
C  WELLS
C
      IF(IUNIT(2).GT.0)THEN
         DO 310 L=1,NWELLS
            K=WELL(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 310
            I=WELL(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 310
            J=WELL(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 310
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 310
C 304-LESS THAN;  310-EQUAL TO;  306-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
            IF (WELL(NWELVL,L)) 304,310,306
 304          IF(MOCTYPE.EQ.1) THEN
                TOTOUT=TOTOUT+(WELL(NWELVL,L)*CNOLD(JS,IS,KS))
              ELSE
                TOTOUT=TOTOUT+(WELL(NWELVL,L)*(CNOLD(JS,IS,KS)+
     *                    CONC(JS,IS,KS))*0.5)
              ENDIF
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT WELLS = RATIO OF WELLS TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
              IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(WELL(NWELVL,L)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
cgzh debug
c      write(iouts,*) 'Mass out well=',SNKMSOUT(JS,IS,KS)
              ENDIF
            GO TO 310
C IF > 0, INTO SUBGRID
 306        TOTIN=TOTIN+(WELL(NWELVL,L)*WELL(IWELLC,L))
 310     CONTINUE
         SBVL(1,7)=SBVL(1,7)+TOTIN*TIMV
         SBVL(2,7)=SBVL(2,7)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,7)=SBVL(4,7)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C
C  RIVERS
C
      IF(IUNIT(4).GT.0)THEN
         DO 350 L=1,NRIVER
            K=RIVR(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 350
            I=RIVR(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 350
            J=RIVR(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 350
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 350
            IF (RIVR(NRIVVL,L)) 354,350,356
 354        TOTOUT=TOTOUT+(RIVR(NRIVVL,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
cgzh bug fix!
corig            GO TO 350
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT RIVERS = RATIO OF RIVERS TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
              IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(RIVR(NRIVVL,L)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
cgzh debug
c      write(iouts,*) 'Mass out river=',SNKMSOUT(JS,IS,KS)
              ENDIF
            GO TO 350
 356        TOTIN=TOTIN+(RIVR(NRIVVL,L)*RIVR(IRIVRC,L))
 350     CONTINUE
         SBVL(1,8)=SBVL(1,8)+TOTIN*TIMV
         SBVL(2,8)=SBVL(2,8)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,8)=SBVL(4,8)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C
C  DRAINS
C
      IF(IUNIT(3).GT.0)THEN
         DO 360 L=1,NDRAIN
            K=DRAI(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 360
            I=DRAI(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 360
            J=DRAI(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 360
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 360
            IF (DRAI(NDRNVL,L)) 364,360,366
 364        TOTOUT=TOTOUT+(DRAI(NDRNVL,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT DRAINS = RATIO OF DRAINS TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
            IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(DRAI(NDRNVL,L)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
            ENDIF
            GO TO 360
 366        WRITE (IOUTS,368) L
 368        FORMAT (/1H ,' *** WARNING -- DRAIN FLUX IS INTO AQUIFER AT 
     * DRAIN NUMBER ',I4,'  ***'/)      
 360     CONTINUE
         SBVL(1,9)=0.0
         SBVL(2,9)=SBVL(2,9)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,9)=SBVL(4,9)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
         TOTOUT=0.0
      END IF
C
C  GENERAL HEAD BOUNDARIES
C
      IF(IUNIT(7).GT.0) THEN
         DO 370 L=1,NBOUND
            K=BNDS(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 370
            I=BNDS(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 370
            J=BNDS(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 370
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 370
            IF (BNDS(NGHBVL,L)) 374,370,376
 374        TOTOUT=TOTOUT+(BNDS(NGHBVL,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT GHB = RATIO OF GHB TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
            IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(BNDS(NGHBVL,L)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
            ENDIF
            GO TO 370
 376        TOTIN=TOTIN+(BNDS(NGHBVL,L)*BNDS(IBNDSC,L))
 370     CONTINUE
         SBVL(1,10)=SBVL(1,10)+TOTIN*TIMV
         SBVL(2,10)=SBVL(2,10)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,10)=SBVL(4,10)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C
C              EVAPOTRANSPIRATION
C
      IF (IUNIT(5).GT.0) THEN
        KE=1
        DO 40 IS=1,NSROW
          I=IS+ISROW1-1
        DO 40 JS=1,NSCOL
          J=JS+ISCOL1-1
          KE=EVTFLO(JS,IS,2)
          IF(IBOUND(J,I,KE).LE.0) GO TO 40
          KS=KE-ISLAY1+1
          EVTRAT=EVTFLO(JS,IS,1)
C  **CONCENTRATION ASSOCIATED WITH E-T IS SET TO ZERO**
          ETCONC=0.0
C  IF AGE PACKAGE IS ACTIVE, DO NOT ALLOW ET TO CONCENTRATE AGE MASS 
          IF(INAGE.GT.0) ETCONC=(CNOLD(JS,IS,KS)+CONC(JS,IS,KS))*0.5
          IF (KS.LT.1.OR.KS.GT.NSLAY) GO TO 40
          IF(EVTRAT.LT.0.0) THEN
            TOTOUT=TOTOUT+EVTRAT*ETCONC
          ELSE IF(EVTRAT.GT.0.0) THEN
          WRITE(IOUTS,*) '***WARNING***  EVT RATE > 0 AT SUBGRID NODE '
     *                ,JS,IS,KS    
          END IF
   40   CONTINUE
        SBVL(1,11)=0.0
        SBVL(2,11)=SBVL(2,11)+TOTOUT*TIMV
        TOTOUT=0.0
      END IF
C
C  FLOW BOUNDARIES (FHB PACKAGE)
C
      IF(IUNIT(16).GT.0) THEN
C  FLOW
         DO 380 L=1,NFLW
            K=IFLLOC(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 380
            I=IFLLOC(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 380
            J=IFLLOC(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 380
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 380
            IF (BDFV(1,L)) 384,380,386
 384        TOTOUT=TOTOUT+(BDFV(1,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT FHB = RATIO OF FHB TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
            IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(BDFV(1,L)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
            ENDIF
            GO TO 380
 386        TOTIN=TOTIN+(BDFV(1,L)*BDFV(IFHBFC+2,L))
 380     CONTINUE
         SBVL(1,12)=SBVL(1,12)+TOTIN*TIMV
         SBVL(2,12)=SBVL(2,12)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,12)=SBVL(4,12)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
         TOTIN=0.0
         TOTOUT=0.0
      ENDIF
C
C---(Note--following code must be revised when multiple solutes are allowed)
C     STREAM PACKAGE
C
      IF (IUNIT(44).GT.0) THEN
C         SBVL(1,21)=SBVL(1,21)+CSTIN(NSOL)*TIMV
C         SBVL(2,21)=SBVL(2,21)+CSTOUT(NSOL)*TIMV
         SBVL(1,21)=SBVL(1,21)+CSTIN(1)*TIMV
         SBVL(2,21)=SBVL(2,21)+CSTOUT(1)*TIMV
      END IF
C
C     LAKE PACKAGE
C
      IF (IUNIT(22).GT.0) THEN
C         SBVL(1,22)=SBVL(1,22)+ALKIN(NSOL)*TIMV
C         SBVL(2,22)=SBVL(2,22)+ALKOUT(NSOL)*TIMV
         SBVL(1,22)=SBVL(1,22)+ALKIN(1)*TIMV
CMOCWT  ALKOUT ARRAY BUILT WITH SNKMSOUT WHICH HAS TIMV IN IT
        IF(PTWTON.EQ.1) THEN
         SBVL(2,22)=SBVL(2,22)+ALKOUT(1)
        ELSE 
         SBVL(2,22)=SBVL(2,22)+ALKOUT(1)*TIMV
	  END IF
      END IF
C
C     DRT PACKAGE
C
      IF (IUNIT(40).GT.0) THEN
	   DO 400 ID=1,NDRTCL
	     K=DRTF(1,ID)
           KS=K-ISLAY1+1
           IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 400
           I=DRTF(2,ID)
           IS=I-ISROW1+1
           IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 400
           J=DRTF(3,ID)
           JS=J-ISCOL1+1
           IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 400
           IF(IBOUND(J,I,K).LE.0) GO TO 400
C 404-LESS THAN;  410-EQUAL TO;  406-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
           IF (DRTF(NDRTVL,ID)) 404,410,406
 404         TOTOUT=TOTOUT+(DRTF(NDRTVL,ID)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
CMOCWT
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT DRT = RATIO OF DRT TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
            IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(DRTF(NDRTVL,ID)/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
            ENDIF
             GO TO 410
C IF > 0, INTO SUBGRID
 406         WRITE (IOUTS,468) ID
 468         FORMAT (/1H ,' *** WARNING -- DRT FLUX IS INTO AQUIFER AT 
     * DRAIN NUMBER ',I4,'  ***'/)      
 410        CONTINUE
c drain returns (sources)
           IF (IDRTFL.GT.0) THEN
             KK = DRTF(6,ID)
             KKS=KK-ISLAY1+1
             IF(KKS.LT.1.OR.KKS.GT.NSLAY) GO TO 400
             IF (KK.NE.0) THEN
              II = DRTF(7,ID)
              IIS=II-ISROW1+1
              IF(IIS.LT.1.OR.IIS.GT.NSROW) GO TO 400
              JJ = DRTF(8,ID)
              JJS=JJ-ISCOL1+1
              IF(JJS.LT.1.OR.JJS.GT.NSCOL) GO TO 400
	       ENDIF
	     ENDIF
C 414-LESS THAN;  420-EQUAL TO;  416-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
           IF (DRTF(NDRTVL-1,ID)) 414,400,416
 414         WRITE (IOUTS,468) ID
 478         FORMAT (/1H ,' *** WARNING -- RETURN FLUX IS OUT AQUIFER AT 
     * DRAIN NUMBER ',I4,'  ***'/)      
             GO TO 400
 416         TOTIN=TOTIN+(DRTF(NDRTVL-1,ID)*CONC(JS,IS,KS))
 400  CONTINUE
         SBVL(1,23)=SBVL(1,23)+TOTIN*TIMV
         SBVL(2,23)=SBVL(2,23)+TOTOUT*TIMV
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,23)=SBVL(4,23)+TOTOUTWT
           TOTOUTWT=0.0
         END IF
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C
C
C     MNW PACKAGE
C
      IF (IUNIT(50).GT.0) THEN
        m = 0
c Loop over all MNW locations
        do while( m .lt. nwell2 )
          m = m + 1
C Only simple MNWs here
         if(MNWid(m).eq.0) then
c get node location 
          n = INT(well2(1,m))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
          IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 512
          IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 512
          IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 512
c   A very large # in WL reference array (8,m) triggers multi-node calculation
c
         if( well2(8,m) .gt. 1.0E30 ) then
c   Set ne = last node in this MNW
          ne  = ifrl(well2(7,m))
c   Save m for other loops
          msave=m
c   LGet C' value from 1st spot of MNWs
          Cinput=well2(4,m)
c   End Loop 1
c   Reset m for next loop
          m = msave  
c   Set ne = last node in this MNW
          ne  = ifrl(well2(7,m))
c   Loop 2: Have C' now, so apply to mass balance counters
c   Loop over nodes in this MNW
          do iin = m, ne
           qnode = well2(17,iin)
c get node location 
           n = INT(well2(1,iin))
           k = (n-1)/(ncol*nrow) + 1
           j = mod((n-1),ncol*nrow)/ncol + 1
           i = mod((n-1),ncol) + 1
c
           JS=I-ISCOL1+1
           IS=J-ISROW1+1
           KS=K-ISLAY1+1
           IF (qnode) 604,610,606
C Only do non-weighted sink calculation for simple nodes (complex nodes
C are handled in GWT1MNW1cx
 604         IF(PTWTON.EQ.0) 
     *          TOTOUT=TOTOUT+(qnode*(CNOLD(JS,IS,KS)+
     *          CONC(JS,IS,KS))*0.5)
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT simple MNWs = RATIO OF MNW TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
            IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(qnode/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
            ENDIF
            GO TO 610
C IF > 0, INTO SUBGRID
 606   TOTIN=TOTIN+(qnode*Cinput)
 610        CONTINUE
          enddo
c   End Loop 2
c   set counter to last node of MNW
          m = ne  
 612        CONTINUE
c
         else
c   if not a part of a MNW, just grab C' and process
		 Cinput=well2(4,m)
C 504-LESS THAN;  510-EQUAL TO;  506-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
           qnode = well2(17,m)
           IF (qnode) 504,510,506
C Only do non-weighted sink calculation for simple nodes (complex nodes
C are handled in GWT1MNW1cx
 504         IF(PTWTON.EQ.0.AND.MNWid(m).EQ.0) 
     *          TOTOUT=TOTOUT+(qnode*(CNOLD(JS,IS,KS)+
     *          CONC(JS,IS,KS))*0.5)
C WITH PARTICLE WEIGHTING OPTION,          
C TOTAL OUT simple MNWs = RATIO OF MNW TO TOTAL SNKFLO * TOTAL MASS OUT
C   SINKS IN CELL (SNKMSOUT) (CALCULATED IN MOVEWT)
            IF(PTWTON.EQ.1) THEN
                TOTOUTWT=TOTOUTWT+(qnode/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
            ENDIF
            GO TO 510
C IF > 0, INTO SUBGRID
 506   TOTIN=TOTIN+(qnode*Cinput)
 510        CONTINUE

         end if
 512        CONTINUE
       end if
      end do
C
         SBVL(1,24)=SBVL(1,24)+TOTIN*TIMV
         SBVL(2,24)=SBVL(2,24)+TOTOUT*TIMV
c   Separate mass balance item: Net mass flux into aquifer by all MNWs 
c   (for complex nodes the sbvl 1 and 2 slots may be different and are 
c    defined in the cx routine)
         SBVL(1,25)=SBVL(1,25)+TOTIN*TIMV
C
CMOCWT   SNKMSOUT HAS TIMV IN IT
         IF (PTWTON.EQ.1) THEN 
           SBVL(4,24)=SBVL(4,24)+TOTOUTWT
c   Separate mass balance item: Net mass flux out of aquifer by all MNWs 
           SBVL(2,25)=SBVL(2,25)+TOTOUTWT
           TOTOUTWT=0.0
         ELSE
c   Separate mass balance item: Net mass flux out of aquifer by all MNWs 
           SBVL(2,25)=SBVL(2,25)+TOTOUT*TIMV
         END IF
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C     ****************************************************************
      RETURN   
C     **************************************************************** 
      END
C
C
C SMOC6BE   CALCULATE SOLUTE MASS BALANCE (ELLAM)
C***********************************************************************
C
      SUBROUTINE SMOC6BE(IBOUND,RF,
     *   CONC,CBNDY,IDBNDY,CNOLD,CFXHBC,CINFL,
     *   THCK,POR,CHDFLO,CTCFLO,LOCBDY,
     *   RECH,CRECH,IRCH,NRCHOP,
     *   WELL,MXWELL,NWELLS,NWELVL,IWELLC,
     *   RIVR,MXRIVR,NRIVER,NRIVVL,IRIVRC,
     *   DRAI,MXDRAI,NDRAIN,NDRNVL,
     *   BNDS,MXBNDS,NBOUND,NGHBVL,IBNDSC,
     *   DRTF,MXDRT,NDRTCL,NDRTVL,IDRTFL,
     *   EVTFLO,IEVT,NEVTOP,
     *   NZIN,TEMP,JAS,IAS,AS,NELTS,
     *   LENW,DECAY,TIMV,CSTIN,CSTOUT,ALKIN,ALKOUT,NSOL,
     *   XFOR,XBAC,YFOR,YBAC,
     *   NCOL,NROW,NLAY,NODESS,NODES,
     *   NSCOL,NSROW,NSLAY,
     *   IUNIT,IOUTS,NCINFL,JRF,SBVL,NTFACE,
     *   NFACES,IABOVE,IBELOW,DELCOL,DELROW,NIUNIT,INAGE,
     *   DAGE,SRCAGE,SAGE,
     *   CINFLA,CINFLB,CINXY,
cgzh mnw
     *   nwell2,mxwel2,well2,MNWid)
C
C
C  STMASS,ADMASS HOLD MASS FROM BOUNDARY NODES UPON ENTRY
C
C
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION SBVL(6,NIUNIT),IUNIT(NIUNIT),
     *  WELL(NWELVL,MXWELL),RIVR(NRIVVL,MXRIVR),
     *  DRAI(NDRNVL,MXDRAI),BNDS(NGHBVL,MXBNDS),
     *  RECH(NCOL,NROW),CRECH(NSCOL,NSROW),IRCH(NCOL,NROW),
     *  EVTFLO(NSCOL,NSROW,2),IEVT(NCOL,NROW),
     *  DRTF(NDRTVL,MXDRT)
      DIMENSION well2(18,mxwel2),MNWid(mxwel2)
C
      DIMENSION IBOUND(NCOL,NROW,NLAY),RF(NSCOL,NSROW,NSLAY),
     *  CHDFLO(NSCOL,NSROW,NSLAY),CTCFLO(NFACES),LOCBDY(3,NFACES),
     *  CONC(NSCOL,NSROW,NSLAY),CNOLD(NSCOL,NSROW,NSLAY),
     *  CINFL(NCINFL),CFXHBC(NSCOL,NSROW,NSLAY),
     *  THCK(NSCOL,NSROW,NSLAY),POR(NSCOL,NSROW,NSLAY)
C
      DIMENSION CBNDY(NTFACE),TEMP(LENW),NZIN(NODESS),
     *          IDBNDY(NFACES),DELCOL(NCOL),DELROW(NROW)
C
cgzh cbdy
      DIMENSION CINFLA(NSCOL,NSROW),CINFLB(NSCOL,NSROW),
     * CINXY(NSCOL,NSROW,NSLAY)
      DIMENSION ALKIN(NSOL),ALKOUT(NSOL)
      DIMENSION CSTIN(NSOL),CSTOUT(NSOL)
      DOUBLE PRECISION DECAY
CMOCWT
      INCLUDE 'ptwt.inc'

      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
C      DOUBLE PRECISION DCYFCT,DCYT,DCYT2
C
C     ****************************************************************   
C
C  ---COMPUTE MASS BALANCE FOR SOLUTE---                              
C
C     ****************************************************************   
      TOTOUT=0.D0 
      TOTIN=0.D0
      DMASS=0.D0
C
C   COMPUTE DECAY OF OLD TIME LEVEL MASS
C
      IF (DECAY.NE.0.D0) DMASS= EXP(-DECAY*TIMV)*OLMASS - OLMASS
C
C  MASS DECAYED AND MASS "CREATED"
C  
      IF (DMASS.LT.0D0) SBVL(2,3)=SBVL(2,3)+DMASS
      IF (DMASS.GT.0D0) SBVL(1,3)=SBVL(1,3)+DMASS
C
C  DECAYED OR "CREATED" SOURCE MASS
C
cgzh not used, SRCDCY accumulates      SRCCUM=SRCCUM+SRCDCY
C
C   ACCUMULATE SOURCE AGE MASS
C
cea
      SRCAGE=SRCAGE+SAGE  
C
C
C   CALCULATE MASS IN STORAGE AND ADSORBED
C
      CALL SMOC5A(JRF,NODES,NODESS,NSCOL,NSROW,NSLAY,
     *          NCOL,NROW,NLAY,NZIN,TEMP,JAS,IAS,AS,NELTS,
     *          CONC,RF,POR,IBOUND,STMASS,ADMASS,LENW,THCK)
C
C  TOTAL STORED MASS AND PRESENT MASS SORBED
C
      SBVL(2,1)=STMASS
C  IF RF.NE.1
      IF (JRF.EQ.1) THEN
         SBVL(2,2)=ADMASS  
      END IF
      OLMASS=STMASS+ADMASS
C
cea
C  ACCUMULATE AGE MASS INCREASE THIS TRANSPORT TIME STEP
C
      SBVL(1,17)=SBVL(1,17)+WATVOL*DAGE
C
C  LOOP OVER CELLS FOR  CONSTANT HEAD NODES
      NODE=0
      DO 270 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 270 IS=1,NSROW
      I=IS+ISROW1-1
      DO 268 JS=1,NSCOL
      J=JS+ISCOL1-1
      IF (JRF.EQ.1) RFVAL=RF(JS,IS,KS)
      NODE=NODE+1
C  SKIP IF NO FLOW NODE
         ICH=IBOUND(J,I,K)
         IF (ICH.EQ.0) GO TO 268
C
C  CONSTANT-HEAD NODES
C
        IF (CHDFLO(JS,IS,KS).LT.0.D0) THEN 
c          TOTOUT=TOTOUT+CHDFLO(JS,IS,KS)*CONC(JS,IS,KS)
          TOTOUT=TOTOUT+CHDFLO(JS,IS,KS)*(CNOLD(JS,IS,KS)
     *       +CONC(JS,IS,KS))*0.5
        ELSE     
          TOTIN=TOTIN+CHDFLO(JS,IS,KS)*CFXHBC(JS,IS,KS)
        ENDIF
  268 CONTINUE        
  270 CONTINUE        
C
C  END OF LOOP OVER CELLS
C
C  TOTAL CONSTANT HEAD NODES
      SBVL(1,4)=SBVL(1,4)+TOTIN*TIMV
      SBVL(2,4)=SBVL(2,4)+TOTOUT*TIMV
      TOTIN=0.D0
      TOTOUT=0.D0
C
C  CALCULATE SOLUTE FLUX ACROSS SUBGRID BOUNDARIES
C
C  GET LOCATION OF EACH SUBGRID BOUNDARY FACE
c  idmax in ellam.inc by smoc5by
      DO 287 IE=1,IDMAX
C 284-LESS THAN;  285-EQUAL TO;  286-GREATER THAN 0
C LESS THAN 0, CTCFLO POINTS OUT OF AQUIFER 
C GREATER THAN 0, CTCFLO POINTS INTO OF AQUIFER
	IF (CTCFLO(IE)) 284,285,286
 284  CONTINUE
	TOTOUT=TOTOUT+CTCFLO(IE)*CBNDY(IDBNDY(IE))
      GO TO 287
 286    K=LOCBDY(3,IE)
	  KS=K-ISLAY1+1
          J=LOCBDY(1,IE)
          I=LOCBDY(2,IE)
          JS=J-ISCOL1+1
          IS=I-ISROW1+1
      IF (IE.LT.IDTOP) THEN
corig	  CPRIME=CINFL(KS)
            CPRIME=CINXY(JS,IS,KS)
	ELSE
cgzh cbdy
            CPRIME=0.0
            IF (MOD(IE,2).EQ.0) THEN
              CPRIME=CINFLB(JS,IS)  
            ELSE
              CPRIME=CINFLA(JS,IS)  
            END IF
	ENDIF
	TOTIN=TOTIN+CTCFLO(IE)*CPRIME
	GO TO 287
 285  CONTINUE
 287  CONTINUE
      SBVL(1,5)=SBVL(1,5)+TOTIN*TIMV
      SBVL(2,5)=SBVL(2,5)+TOTOUT*TIMV
      TOTIN=0.D0
      TOTOUT=0.D0
C
C  RECHARGE
C
      IF (IUNIT(8).GT.0) THEN
        KR=1
        DO 20 IS=1,NSROW
          I=IS+ISROW1-1
        DO 20 JS=1,NSCOL
          J=JS+ISCOL1-1
          IF(NRCHOP.NE.1) KR=IRCH(J,I)
          IF(IBOUND(J,I,KR).LE.0) GO TO 20
          KS=KR-ISLAY1+1
          RCHRAT=RECH(J,I)
          IF (KS.LT.1.OR.KS.GT.NSLAY) GO TO 20
C IF < 0, OUT OF SUBGRID
          IF(RCHRAT.LT.0D0) THEN
            TOTOUT=TOTOUT+RCHRAT*(CNOLD(JS,IS,KS)+CONC(JS,IS,KS))*0.5
C IF > 0, INTO SUBGRID
          ELSE IF(RCHRAT.GT.0D0) THEN
            TOTIN=TOTIN+RCHRAT*CRECH(JS,IS)
          END IF
   20 CONTINUE
        SBVL(1,6)=SBVL(1,6)+TOTIN*TIMV
        SBVL(2,6)=SBVL(2,6)+TOTOUT*TIMV
        TOTIN=0.D0
        TOTOUT=0.D0
      END IF
C        
C  WELLS
C
      IF(IUNIT(2).GT.0)THEN
         DO 310 L=1,NWELLS
            K=WELL(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 310
            I=WELL(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 310
            J=WELL(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 310
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 310
C 304-LESS THAN;  310-EQUAL TO;  306-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
            IF (WELL(NWELVL,L)) 304,310,306
 304        TOTOUT=TOTOUT+WELL(NWELVL,L)*(CNOLD(JS,IS,KS)+
     *    CONC(JS,IS,KS))*0.5
            GO TO 310
C IF > 0, INTO SUBGRID
 306        TOTIN=TOTIN+(WELL(NWELVL,L)*WELL(IWELLC,L))
 310     CONTINUE
         SBVL(1,7)=SBVL(1,7)+TOTIN*TIMV
         SBVL(2,7)=SBVL(2,7)+TOTOUT*TIMV
         TOTIN=0.D0
         TOTOUT=0.D0
      END IF
C
C  RIVERS
C
      IF(IUNIT(4).GT.0)THEN
         DO 350 L=1,NRIVER
            K=RIVR(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 350
            I=RIVR(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 350
            J=RIVR(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 350
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 350
            IF (RIVR(NRIVVL,L)) 354,350,356
 354        TOTOUT=TOTOUT+(RIVR(NRIVVL,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
            GO TO 350
 356        TOTIN=TOTIN+(RIVR(NRIVVL,L)*RIVR(IRIVRC,L))
 350     CONTINUE
         SBVL(1,8)=SBVL(1,8)+TOTIN*TIMV
         SBVL(2,8)=SBVL(2,8)+TOTOUT*TIMV
         TOTIN=0.D0
         TOTOUT=0.D0
      END IF
C
C  DRAINS
C
      IF(IUNIT(3).GT.0)THEN
         DO 360 L=1,NDRAIN
            K=DRAI(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 360
            I=DRAI(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 360
            J=DRAI(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 360
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 360
            IF (DRAI(NDRNVL,L)) 364,360,366
 364        TOTOUT=TOTOUT+(DRAI(NDRNVL,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
            GO TO 360
 366        WRITE (IOUTS,368) NDRAIN
 368        FORMAT (/1H ,' *** WARNING -- DRAIN FLUX IS INTO AQUIFER AT 
     * DRAIN NUMBER ',I4,'  ***'/)      
 360     CONTINUE
         SBVL(1,9)=0.D0
         SBVL(2,9)=SBVL(2,9)+TOTOUT*TIMV
         TOTOUT=0.D0
      END IF
C
C  GENERAL HEAD BOUNDARIES
C
      IF(IUNIT(7).GT.0) THEN
         DO 370 L=1,NBOUND
            K=BNDS(1,L)
            KS=K-ISLAY1+1
            IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 370
            I=BNDS(2,L)
            IS=I-ISROW1+1
            IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 370
            J=BNDS(3,L)
            IF(IBOUND(J,I,K).LE.0) GO TO 370
            JS=J-ISCOL1+1
            IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 370
            IF (BNDS(NGHBVL,L)) 374,370,376
 374        TOTOUT=TOTOUT+(BNDS(NGHBVL,L)*(CNOLD(JS,IS,KS)+
     *            CONC(JS,IS,KS))*0.5)
            GO TO 370
 376        TOTIN=TOTIN+(BNDS(NGHBVL,L)*BNDS(IBNDSC,L))
 370     CONTINUE
         SBVL(1,10)=SBVL(1,10)+TOTIN*TIMV
         SBVL(2,10)=SBVL(2,10)+TOTOUT*TIMV
         TOTIN=0.D0
         TOTOUT=0.D0
      END IF
C
C              EVAPOTRANSPIRATION
C
      IF (IUNIT(5).GT.0) THEN
        KE=1
        DO 40 IS=1,NSROW
          I=IS+ISROW1-1
        DO 40 JS=1,NSCOL
          J=JS+ISCOL1-1
          KE=EVTFLO(JS,IS,2)
          IF(IBOUND(J,I,KE).LE.0) GO TO 40
          KS=KE-ISLAY1+1
          EVTRAT=EVTFLO(JS,IS,1)
C  **CONCENTRATION ASSOCIATED WITH E-T IS SET TO ZERO**
          ETCONC=0.D0
C  IF AGE PACKAGE IS ACTIVE, DO NOT ALLOW ET TO CONCENTRATE AGE MASS 
          IF(INAGE.GT.0) ETCONC=CNOLD(JS,IS,KS)
          IF (KS.LT.1.OR.KS.GT.NSLAY) GO TO 40
          IF(EVTRAT.LT.0D0) THEN
            TOTOUT=TOTOUT+EVTRAT*ETCONC
          ELSE IF(EVTRAT.GT.0D0) THEN
          WRITE(IOUTS,*) '***WARNING***  EVT RATE > 0 AT SUBGRID NODE '
     *                ,JS,IS,KS    
          END IF
   40   CONTINUE
        SBVL(1,11)=0.D0
        SBVL(2,11)=SBVL(2,11)+TOTOUT*TIMV
        TOTOUT=0.D0
      END IF
C
C---(Note--following code must be revised when multiple solutes are allowed)
C     STREAM PACKAGE
C
      IF (IUNIT(44).GT.0) THEN
C         SBVL(1,21)=SBVL(1,21)+CSTIN(NSOL)*TIMV
C         SBVL(2,21)=SBVL(2,21)+CSTOUT(NSOL)*TIMV
         SBVL(1,21)=SBVL(1,21)+CSTIN(1)*TIMV
         SBVL(2,21)=SBVL(2,21)+CSTOUT(1)*TIMV
      END IF
C
C     LAKE PACKAGE
C
      IF (IUNIT(22).GT.0) THEN
C         SBVL(1,22)=SBVL(1,22)+ALKIN(NSOL)*TIMV
C         SBVL(2,22)=SBVL(2,22)+ALKOUT(NSOL)*TIMV
         SBVL(1,22)=SBVL(1,22)+ALKIN(1)*TIMV
         SBVL(2,22)=SBVL(2,22)+ALKOUT(1)*TIMV
      END IF
C
C     DRT PACKAGE
C
      IF (IUNIT(40).GT.0) THEN
	   DO 400 ID=1,NDRTCL
	     K=DRTF(1,ID)
           KS=K-ISLAY1+1
           IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 400
           I=DRTF(2,ID)
           IS=I-ISROW1+1
           IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 400
           J=DRTF(3,ID)
           JS=J-ISCOL1+1
           IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 400
           IF(IBOUND(J,I,K).LE.0) GO TO 400
C 404-LESS THAN;  410-EQUAL TO;  406-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
           IF (DRTF(NDRTVL,ID)) 404,410,406
 404         TOTOUT=TOTOUT+(DRTF(NDRTVL,ID)*CONC(JS,IS,KS))
             GO TO 410
C IF > 0, INTO SUBGRID
 406         WRITE (IOUTS,468) ID
 468         FORMAT (/1H ,' *** WARNING -- DRT FLUX IS INTO AQUIFER AT 
     * DRAIN NUMBER ',I4,'  ***'/)      
 410        CONTINUE
c drain returns (sources)
           IF (IDRTFL.GT.0) THEN
             KK = DRTF(6,ID)
             KKS=KK-ISLAY1+1
             IF(KKS.LT.1.OR.KKS.GT.NSLAY) GO TO 400
             IF (KK.NE.0) THEN
              II = DRTF(7,ID)
              IIS=II-ISROW1+1
              IF(IIS.LT.1.OR.IIS.GT.NSROW) GO TO 400
              JJ = DRTF(8,ID)
              JJS=JJ-ISCOL1+1
              IF(JJS.LT.1.OR.JJS.GT.NSCOL) GO TO 400
	       ENDIF
	     ENDIF
C 414-LESS THAN;  420-EQUAL TO;  416-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
           IF (DRTF(NDRTVL-1,ID)) 414,400,416
 414         WRITE (IOUTS,468) ID
 478         FORMAT (/1H ,' *** WARNING -- RETURN FLUX IS OUT AQUIFER AT 
     * DRAIN NUMBER ',I4,'  ***'/)      
             GO TO 400
 416         TOTIN=TOTIN+(DRTF(NDRTVL-1,ID)*CONC(JS,IS,KS))
 400  CONTINUE
         SBVL(1,23)=SBVL(1,23)+TOTIN*TIMV
         SBVL(2,23)=SBVL(2,23)+TOTOUT*TIMV
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C
C
C     MNW PACKAGE
C
      IF (IUNIT(50).GT.0) THEN
        m = 0
c Loop over all MNW locations
        do while( m .lt. nwell2 )
          m = m + 1
C Only simple MNWs here
         if(MNWid(m).eq.0) then
c get node location 
          n = INT(well2(1,m))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
          IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 512
          IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 512
          IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 512
c   A very large # in WL reference array (8,m) triggers multi-node calculation
c
         if( well2(8,m) .gt. 1.0E30 ) then
c   Set ne = last node in this MNW
          ne  = ifrl(well2(7,m))
c   Save m for other loops
          msave=m
c   LGet C' value from 1st spot of MNWs
          Cinput=well2(4,m)
c   End Loop 1
c   Reset m for next loop
          m = msave  
c   Set ne = last node in this MNW
          ne  = ifrl(well2(7,m))
c   Loop 2: Have C' now, so apply to mass balance counters
c   Loop over nodes in this MNW
          do iin = m, ne
           qnode = well2(17,iin)
c get node location 
           n = INT(well2(1,iin))
           k = (n-1)/(ncol*nrow) + 1
           j = mod((n-1),ncol*nrow)/ncol + 1
           i = mod((n-1),ncol) + 1
c
           JS=I-ISCOL1+1
           IS=J-ISROW1+1
           KS=K-ISLAY1+1
           IF (qnode) 604,610,606
 604         TOTOUT=TOTOUT+(qnode*(CNOLD(JS,IS,KS)+
     *          CONC(JS,IS,KS))*0.5)
            GO TO 610
C IF > 0, INTO SUBGRID
 606   TOTIN=TOTIN+(qnode*Cinput)
 610        CONTINUE
          enddo
c   End Loop 2
c   set counter to last node of MNW
          m = ne  
 612        CONTINUE
c
         else
c   if not a part of a MNW, just grab C' and process
		 Cinput=well2(4,m)
C 504-LESS THAN;  510-EQUAL TO;  506-GREATER THAN 0
C IF < 0, OUT OF SUBGRID
           qnode = well2(17,m)
           IF (qnode) 504,510,506
C Only do non-weighted sink calculation for simple nodes (complex nodes
C are handled in GWT1MNW1cx
 504         IF(MNWid(m).EQ.0) 
     *          TOTOUT=TOTOUT+(qnode*(CNOLD(JS,IS,KS)+
     *          CONC(JS,IS,KS))*0.5)
            GO TO 510
C IF > 0, INTO SUBGRID
 506   TOTIN=TOTIN+(qnode*Cinput)
 510        CONTINUE

         end if
 512        CONTINUE
       end if
      end do
C
         SBVL(1,24)=SBVL(1,24)+TOTIN*TIMV
         SBVL(2,24)=SBVL(2,24)+TOTOUT*TIMV
c   Separate mass balance item: Net mass flux into aquifer by all MNWs 
c   (for complex nodes the sbvl 1 and 2 slots may be different and are 
c    defined in the cx routine)
         SBVL(1,25)=SBVL(1,25)+TOTIN*TIMV
c   Separate mass balance item: Net mass flux out of aquifer by all MNWs 
         SBVL(2,25)=SBVL(2,25)+TOTOUT*TIMV
         TOTIN=0.0
         TOTOUT=0.0
      END IF
C     ****************************************************************
      RETURN   
C     **************************************************************** 
      END
C
C SMOC6M    PRINT OUTPUT FOR MASS BALANCE
C*************************************************************************
C
      SUBROUTINE SMOC6M(SBVL,SRCDCY,TIMV,
     * IOUTS,KPER,NPER,KSTP,NSTP,IMOV,NMOV,SUMTCH,ICONLY,
     * INDK,IDKZO,IDKFO,IDKZS,IDKFS,
     * INAGE,
     * INDP,IDPZO,IDPFO,INSFRUNIT,INLAKUNIT,INDRT,NIUNIT,
     * SRCAGE,INMNW)
C
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
cgzh debug double mass balance numbers
      DOUBLE PRECISION CSTM2,TOTMIN,TOTMOT,RESID,TIN,TOUT,
     *  DENOM,ERR,SMIN,SMPR,ERR2
      DIMENSION SBVL(6,NIUNIT)
C
C*************************************************************************
C
C  ZERO ACCUMULATORS
      TOTMIN=0.0
      TOTMOT=0.0
C  ADD MASS FLUX FROM DECAY FROM SOURCES TO THE TOTAL MASS OUT
cgzh decayed src never in aquifer, don't include in amount out
c      TOTMOT=TOTMOT+SRCCUM
C  accumulated in SRCDCY, no need for srccum any more
c      SRCCUM=SRCDCY*TIMV
      TOTMIN=TOTMIN+SRCAGE
C  LOOP THROUGH PACKAGES, ADDING APPROPRIATE MASS FLUXES
      DO 50 L=3,12
         TOTMIN=TOTMIN+SBVL(1,L)
         TOTMOT=TOTMOT+SBVL(2,L)
 50      CONTINUE
cgzh decayed src never in aquifer, don't include in amount out
c      TOTMOT=TOTMOT+SRCDCY
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
        TOTMOT=TOTMOT+SBVL(2,L)
 61   CONTINUE
C  CALCULATE CHANGE IN MASS STORED
cgzh orig       CSTM2=SBVL(1,1)-SBVL(2,1)+SBVL(1,2)-SBVL(2,2)
      CSTM2=SBVL(2,1)-SBVL(1,1)+SBVL(2,2)-SBVL(1,2)
cgzh orig      IF(INDP.GT.0) CSTM2=CSTM2+SBVL(1,18)-SBVL(2,18)
      IF(INDP.GT.0) CSTM2=CSTM2-SBVL(1,18)+SBVL(2,18)
      WRITE (IOUTS,100)
      WRITE (IOUTS,105) KPER,NPER,KSTP,NSTP,IMOV,NMOV
      WRITE (IOUTS,107) SUMTCH
cgzh mnw
      IF(INMNW.NE.0) THEN
       WRITE(IOUTS,*) 
       WRITE(IOUTS,*) '     External MNW mass flux (based on wellhead
     & flow) (L**3)(M/VOL)'
       WRITE(IOUTS,*) '     -----------------------------------------
     &--------------------'
       WRITE(IOUTS,'(A,1PE11.4)')
     &'       Cumulative mass flux into aquifer by MNWs   = ',SBVL(1,25)
       WRITE(IOUTS,'(A,1PE11.4)')
     &'       Cumulative mass flux out of aquifer by MNWs = ',SBVL(2,25)
       WRITE(IOUTS,*) 
      END IF
c
      WRITE (IOUTS,110) SBVL(1,1),SBVL(1,2),SBVL(2,1),SBVL(2,2)
      IF(INDP.GT.0) WRITE (IOUTS,112) SBVL(1,18),SBVL(2,18)
      WRITE (IOUTS,111) CSTM2
      WRITE (IOUTS,115) (SBVL(1,L),L=3,12)
      IF(INSFRUNIT.GT.0) WRITE (IOUTS,133) SBVL(1,21)
      IF(INLAKUNIT.GT.0) WRITE (IOUTS,134) SBVL(1,22)
      IF(INDK.GT.0) THEN
        IF(IDKZO.EQ.1) WRITE (IOUTS,116) SBVL(1,13)
        IF(IDKFO.EQ.1) WRITE (IOUTS,117) SBVL(1,14)
        IF(IDKZS.EQ.1) WRITE (IOUTS,118) SBVL(1,15)
        IF(IDKFS.EQ.1) WRITE (IOUTS,119) SBVL(1,16)
      END IF
      IF(INAGE.GT.0) WRITE (IOUTS,120) SBVL(1,17)
      IF(INDP.GT.0) THEN
        IF(IDPZO.EQ.1) WRITE (IOUTS,121) SBVL(1,19)
        IF(IDPFO.EQ.1) WRITE (IOUTS,122) SBVL(1,20)
      END IF
      IF(INDRT.GT.0) WRITE (IOUTS,137) SBVL(1,23)
      IF(INMNW.GT.0) WRITE (IOUTS,140) SBVL(1,24)
      WRITE (IOUTS,123) TOTMIN
      WRITE (IOUTS,125) (SBVL(2,L),L=3,12)
      IF(INSFRUNIT.GT.0) WRITE (IOUTS,135) SBVL(2,21)
      IF(INLAKUNIT.GT.0) WRITE (IOUTS,136) SBVL(2,22)
      IF(INDK.GT.0) THEN
        IF(IDKZO.EQ.1) WRITE (IOUTS,116) SBVL(2,13)
        IF(IDKFO.EQ.1) WRITE (IOUTS,117) SBVL(2,14)
        IF(IDKZS.EQ.1) WRITE (IOUTS,118) SBVL(2,15)
        IF(IDKFS.EQ.1) WRITE (IOUTS,119) SBVL(2,16)
      END IF
      IF(INDP.GT.0) THEN
        IF(IDPZO.EQ.1) WRITE (IOUTS,121) SBVL(2,19)
        IF(IDPFO.EQ.1) WRITE (IOUTS,122) SBVL(2,20)
      END IF
      IF(INDRT.GT.0) WRITE (IOUTS,138) SBVL(2,23)
      IF(INMNW.GT.0) WRITE (IOUTS,140) SBVL(2,24)
      WRITE (IOUTS,126) TOTMOT
      WRITE (IOUTS,139) SRCDCY
      IF(INAGE.GT.0) THEN
        WRITE (IOUTS,141) SRCAGE
      ENDIF
C  CALCULATE RESIDUAL
cgzh initialize ERR in case of no solute at all in system
      ERR=0.0
cgzh orig      RESID=TOTMIN+TOTMOT+CSTM2
cgzh debug keep source decay as separate itemization, but must be in resid calc
      RESID=TOTMIN+TOTMOT+SRCDCY-CSTM2
cgzh for error calc, take decayed source out of total mass in
      TIN=ABS(TOTMIN)+SRCDCY
      TOUT=ABS(TOTMOT)
C  ZERO MASS BALANCE FORMULA FLAGS
      IMBF=0
      IMBF2=0
      IPREC=0
C  CHECK FLAGS FOR MASS FLUXES, CALCULATE WITH ACCORDING
C  FORMULA
      IF(((ICONLY.EQ.1).AND.(SBVL(1,1).NE.0.0)).OR.
     *   ((TIN.EQ.0.0).AND.(TOUT.EQ.0.0))) GOTO 97
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
      SMPR=SBVL(2,1)+SBVL(2,2)
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
      IF((ICONLY.EQ.1).OR.((TIN.EQ.0.0).AND.(TOUT.EQ.0.0))) THEN
         DENOM=SMPR
         IMBF2=1
      ENDIF
cgzh mbprec use stored mass error calculation when precision concern with 
c    flux in/out
      IF(IPREC.EQ.1) THEN
         DENOM=SMPR
	   IMBF2=1
      END IF
  95  CONTINUE
      WRITE(IOUTS,150) RESID
C  PRINT OUT MASS BALANCE 
      IF(IMBF.EQ.1) WRITE (IOUTS,130) ERR
      IF(IMBF.EQ.2) WRITE (IOUTS,131) ERR
	IF(IMBF2.GT.0) THEN
         ERR2=(RESID*100.0)/DENOM
         WRITE(IOUTS,132) ERR2
	ENDIF 
cgzh mbprec
      IF(IPREC.EQ.1) THEN
        WRITE(IOUTS,*) '***NOTE*** THE CHANGE IN MASS STORED IS LESS
     * THAN 7 ORDERS OF MAGNITUDE SMALLER'
        WRITE(IOUTS,*) 'THAN THE MASS DISSOLVED.  THEREFORE, WITH
     * SINGLE-PRECISION, IT MAY NOT BE '
        WRITE(IOUTS,*) 'POSSIBLE TO ACCURATELY COMPUTE SOLUTE BUDGET
     * AND MASS BALANCE RESIDUAL NUMBERS.'
      END IF
      WRITE(IOUTS,'(//)')
 100  FORMAT (//1H ,10X,'SOLUTE BUDGET AND MASS BALANCE FOR TRANSPORT SU
     *BGRID'/)
 105  FORMAT (1H ,5X,'VALUES CALCULATED AT END OF: '/
     *15X,'STRESS PERIOD  ',I6,'  OUT OF ',I6/
     *14X,'FLOW TIME STEP  ',I6,'  OUT OF ',I6/
     *04X,'TRANSPORT TIME INCREMENT  ',I6,'  OUT OF ',I6/)
 107  FORMAT (1H ,5X,'ELAPSED TIME = ', 1PE11.4)
 110  FORMAT (/1H ,05X,'CHEMICAL MASS IN STORAGE: '/
     *10X,'INITIAL:   MASS DISSOLVED = ',1PE11.4,5X,
     *'MASS SORBED = ',1PE11.4/
     *10X,'PRESENT:   MASS DISSOLVED = ',1PE11.4,5X,
     *'MASS SORBED = ',1PE11.4)
 111  FORMAT (/
     *15X,'CHANGE IN MASS STORED = ',1PE11.4//)
 112  FORMAT (
     *10X,'INITIAL:MASS DOUBLE POROS = ',1PE11.4/
     *10X,'PRESENT:MASS DOUBLE POROS = ',1PE11.4)
 115  FORMAT (//05X,'CUMULATIVE SOLUTE MASS  (L**3)(M/VOL)'/
     *          05X,'----------------------'//10X,'IN:'/10X,'---'/
     *05X,'                DECAY = ',1PE11.4/
     *05X,'        CONSTANT HEAD = ',1PE11.4/
     *05X,'     SUBGRID BOUNDARY = ',1PE11.4/
     *05X,'             RECHARGE = ',1PE11.4/
     *05X,'                WELLS = ',1PE11.4/
     *05X,'               RIVERS = ',1PE11.4/
     *05X,'               DRAINS = ',1PE11.4/
     *05X,'GENL. HEAD-DEP. BDYS. = ',1PE11.4/
     *05X,'   EVAPOTRANSPIRATION = ',1PE11.4/
     *05X,' SPECIFIED FLOW (FHB) = ',1PE11.4)
 116  FORMAT (
     *05X,'    ZERO-ORDER GROWTH = ',1PE11.4)
 117  FORMAT (
     *05X,'    FIRST-ORDER DECAY = ',1PE11.4)
 118  FORMAT (
     *05X,'0-ORDER GROWTH SORBED = ',1PE11.4)
 119  FORMAT (
     *05X,' 1-ORDER DECAY SORBED = ',1PE11.4)
 120  FORMAT (
     *05X,'    AGE-MASS INCREASE = ',1PE11.4)
 121  FORMAT (
     *05X,' DP ZERO-ORDER GROWTH = ',1PE11.4)
 122  FORMAT (
     *05X,' DP FIRST-ORDER DECAY = ',1PE11.4)
 123  FORMAT (/
     *05X,'             TOTAL IN = ',1PE11.4/)
 125  FORMAT (//09X,'OUT:'/09X,'----'/
     *05X,'                DECAY = ',1PE11.4/
     *05X,'        CONSTANT HEAD = ',1PE11.4/
     *05X,'     SUBGRID BOUNDARY = ',1PE11.4/
     *05X,'             RECHARGE = ',1PE11.4/
     *05X,'                WELLS = ',1PE11.4/
     *05X,'               RIVERS = ',1PE11.4/
     *05X,'               DRAINS = ',1PE11.4/
     *05X,'GENL. HEAD-DEP. BDYS. = ',1PE11.4/
     *05X,'   EVAPOTRANSPIRATION = ',1PE11.4/
     *05X,' SPECIFIED FLOW (FHB) = ',1PE11.4)
 126  FORMAT (/
     *05X,'            TOTAL OUT = ',1PE11.4/)
 139  FORMAT (/05X,'    SOURCE-TERM DECAY = ',1PE11.4/)
 130  FORMAT (05X,'  PERCENT DISCREPANCY = ',1PE11.4,
     *' RELATIVE TO MASS FLUX IN'/)
 131  FORMAT (05X,'  PERCENT DISCREPANCY = ',1PE11.4,
     *' RELATIVE TO MASS FLUX OUT'/)
 132  FORMAT (05X,'  PERCENT DISCREPANCY = ',1PE11.4,
     *' RELATIVE TO INITIAL MASS STORED'/)
 133  FORMAT (
     *05X,'LOSING STREAM REACHES = ',1PE11.4)
 134  FORMAT (
     *05X,'    LOSING LAKE CELLS = ',1PE11.4)
 135  FORMAT (
     *05X,' GAINING STREAM RCHS. = ',1PE11.4)
 136  FORMAT (
     *05X,'   GAINING LAKE CELLS = ',1PE11.4)
 137  FORMAT (
     *05X,' DRAIN RETURNS(DRT)   = ',1PE11.4)
 138  FORMAT (
     *05X,' DRAINS (DRT)         = ',1PE11.4)
 140  FORMAT (
     *05X,'     MULTI-NODE WELLS = ',1PE11.4)
 141  FORMAT (/05X,'    SOURCE-TERM AGE   = ',1PE11.4/)
 150  FORMAT (//05X,'             RESIDUAL = ',1PE11.4/)
      RETURN
      END
C
C
C  SMOC5BY   CALCULATE FLOWS ACROSS BOUNDARIES OF SUBGRID
C*************************************************************************
C
      SUBROUTINE SMOC5BY(CTCFLO,BUFF,LOCBDY,
     *  IDIR,ID,NCOL,NROW,NLAY,IC1,IC2,IR1,IR2,IL1,IL2,NFACES)
C
      DIMENSION BUFF(NCOL,NROW,NLAY),CTCFLO(NFACES),LOCBDY(3,NFACES)
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
C*************************************************************************
C
C  DEPENDING ON IDIR (DIRECTION OF FLOW),  SAVE LOCATION OF CELL FACES
C  ON THE SUBGRID BOUNDARY, IF THEY EXIST
C
C  COLUMN DIRECTION
      IF (IDIR.EQ.1) THEN
      ID=0
         K1=IL1
         K2=IL2
         I1=IR1
         I2=IR2
         J1=IC1-1
         J2=IC2
C RETURN IF ONLY ONE COLUMN
      IF (NCOL.EQ.1) THEN
       RETURN
      ELSE
         DO 31 K=K1,K2
         DO 32 I=I1,I2
C IF FIRST SUBGRID COLUMN > 1, ADD LOCATION AND FLUX TO LIST
            IF(IC1.GT.1) THEN
               ID=ID+1
               LOCBDY(1,ID)=J1+1
               LOCBDY(2,ID)=I
               LOCBDY(3,ID)=K
C BUFF HAS SAME SIGN CONVENTION, <0 IS OUT, >0 IN IN
               CTCFLO(ID)=BUFF(J1,I,K)
            ENDIF
C IF LAST COLUMN < TOTAL GRID COLUMNS, ADD LOCATION AND FLUX TO LIST
            IF(IC2.LT.NCOL) THEN
               ID=ID+1
               LOCBDY(1,ID)=J2
               LOCBDY(2,ID)=I
               LOCBDY(3,ID)=K
C MUST FLIP BUFF'S SIGN CONVENTION ON THE OTHER SIDE
               CTCFLO(ID)=-BUFF(J2,I,K)
            ENDIF
 32         CONTINUE
 31         CONTINUE
        END IF
      END IF
C  NOW ROW DIRECTION
      IF (IDIR.EQ.2) THEN
         K1=IL1
         K2=IL2
         I1=IR1-1
         I2=IR2
         J1=IC1
         J2=IC2
      IF(NROW.EQ.1) THEN
         RETURN
      ELSE
         DO 51 K=K1,K2
         DO 51 J=J1,J2
            IF(IR1.GT.1) THEN
               ID=ID+1
               LOCBDY(1,ID)=J
               LOCBDY(2,ID)=I1+1
               LOCBDY(3,ID)=K
               CTCFLO(ID)=BUFF(J,I1,K)
            ENDIF
            IF(IR2.LT.NROW) THEN
               ID=ID+1
               LOCBDY(1,ID)=J
               LOCBDY(2,ID)=I2
               LOCBDY(3,ID)=K
               CTCFLO(ID)=-BUFF(J,I2,K)
            ENDIF
 51       CONTINUE
        END IF
      END IF
C AND LAYER DIRECTION
      IF (IDIR.EQ.3) THEN
         K1=IL1-1
         K2=IL2
         I1=IR1
         I2=IR2
         J1=IC1
         J2=IC2
cellam
         IDTOP=ID+1
cellam
         IF (NLAY.EQ.1) THEN
            RETURN
         ELSE
         DO 64 I=I1,I2
         DO 64 J=J1,J2
            IF(IL1.GT.1) THEN
               ID=ID+1
               LOCBDY(1,ID)=J
               LOCBDY(2,ID)=I  
               LOCBDY(3,ID)=K1+1
               CTCFLO(ID)=BUFF(J,I,K1)
            ENDIF
            IF(IL2.LT.NLAY) THEN
               ID=ID+1
               LOCBDY(1,ID)=J
               LOCBDY(2,ID)=I
               LOCBDY(3,ID)=K2
               CTCFLO(ID)=-BUFF(J,I,K2)
            ENDIF
 64      CONTINUE
        END IF
      END IF
C
cellam
      IDMAX=ID
cellam
      RETURN
      END
