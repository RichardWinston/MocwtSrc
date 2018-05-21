C
C  MOVEWT: MOVE PARTICLES USING WEIGHTED OPTION
C     ***************************************************************
C
      SUBROUTINE MOVEWT(PC,PR,PL,PCONC,IPTID,
     *  VC,VR,VL,
     *  RF,THCK,POR,CINFL,
! RBW CINFL isn't used
     *  SRCFLO,SNKFLO,CONC,
     *  CAVG,IBOUND,
     *  CNOLD,SUMMASS,
     *  NPCELL,NPOLD,LIMBO,
     *  PNEWC,PNEWR,PNEWL,IGENPT,
     *  NEWPTS,NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,NPMAX,NLIMBO,
     *  IOUTS,DECAY,IMOV,TIMV,NP,
     *  NCINFL,
! RBW NCINFL is only used to give the size of CINFL (which isn't used).     
     *  IABOVE,IBELOW,
! RBW IABOVE and IBELOW aren't used     
     *  VCMAX,VRMAX,TLMIN,ICONLY,
     *  DKFO,DKFS,INDK,IDKFO,IDKFS,
CMOCWT
     *  PTWT,SUMWT,SRCAVC,BDYSRC,BDYSNK,BDYSOL,SUMSGMS,SUMSGWT,
! RBW BDYSOL isn't used.     
     *  SGMSOUT,SNKMSOUT,
     *  SUMBFMS,SUMBFWT,
! RBW SUMBFMS and SUMBFWT are set to zero but aren't used otherwise.     
     *  RESIDWT,RESIDC,CELVOL,NZERO,SBVL,NIUNIT,NSKIP,
! RBW RESIDC isn't used     
cgzh varpt
     *  INIPDL,INIPDA,PCORIG,PRORIG,PLORIG,
     *  NPTCOLA,NPTROWA,NPTLAYA,IDIM,NMOV,INTRPL,
     *  REMCRIT,IRAND,
cgzh randpt
     *  ISEEDPT,
cgzh randbd
     *  ISEEDBD,
cgzh mnw
     *  SRCMNW,SNKMNW,INMNW,
c rbw mnw2     
     *  INMNW2,MNWMAX,MNW2, MNWNOD, NODTOT,NMNWVL,
cgzh srcfix
     *  ISRCID,MAXSRC,NPCELLMAX,SRCFAC,NPORIG,WTFAC,HNEW,INBFLX,
cgzh srcfix2
     *  SS_SRC,SS_SNK,ISRCFIX,
     &  TOTFLUXIN,TOTFLUXOT,ICONTRIBIN,ICONTRIBOT,    
! RBW ICONTRIBOT isn't used.       
cgzh ptob
     *  INUNITPTOB,SUMTCH,IPTOBLST,NUMPTOB,IPTOBMNW,NUMPTOB_MNW,
cgzh ccbd
     *  INCCBD,CCBDY,
cgzh et
cgzh debug kkper is debug
     *  INEVT,NEVTOP,IEVT,EVTFLO,INAGE,kkper)
! RBW NEVTOP and IEVT aren't used.     
C
C     ***************************************************************
CMOCWT
C VARIABLES USED JUST IN THIS SUBROUTINE
c      REAL SOURCEQ,FRACWT,FRACMS,RESIDMS,RESIDCL,NEWC,NEWR,NEWL
c      REAL TEMPMS,TEMPWT,SINKQ,RESIDWT,REMVWT,REMVMS
      REAL RESIDCL,NEWC,NEWR,NEWL
cgzh crit
      REAL CRITWT,CRITMS,NPSUM
cgzh mnw
      REAL SRCMNW,SNKMNW
      DOUBLE PRECISION CELVOL,wtdiff,dpcrit
      DOUBLE PRECISION HNEW,HMIN
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
cgzh debug double sgmsout + snkmsout
      DOUBLE PRECISION SGMSOUT,SNKMSOUT
cgzh debug double ptwt, sumwt
      DOUBLE PRECISION PTWT,SUMWT,SUMSGWT,SUMSGMS,SUMBFWT,SUMBFMS
      DOUBLE PRECISION FRACWT,FRACMS,TEMPMS,TEMPWT,
     * SOURCEQ,SINKQ,REMVWT,REMVMS,RESIDWT,SUMMASS,FRAC1,FRAC2,
cgzh srcfix2
     * SS_WT,SS_MS,DEFICIT,DEFMS,WTOUT,WTOFF,WTIN,
     * DEFICIT2,TEWT,TEMS,PTWTOFF,
cgzh et
     * FRACET,FRACSNK,TEMPPTWT,ETFLO
cgzh debug temp variables
      DOUBLE PRECISION SUMVOL,SUMVOL2,SMCELL,CUMMASS,CUMCELL,cumcell2,
     * sumtemp,sumrem,sumremms,addwt,addms,sumPT,sumvol3
cgzh srcnew
      DOUBLE PRECISION SRCMS,SRCWT,QMNWSINK,QSS_SINK
cgzh drycell
      DOUBLE PRECISION DRYWT,DRYMS,TEMPDRYWT,TEMPDRYMS
cgzh ccbd
      DOUBLE PRECISION TOTIN,TOTOUT,DIFFMASS,CHNGMASS
      INTEGER NPNEW 
C TEMPORARY ARRAYS
      ALLOCATABLE TEMPMS(:,:,:),TEMPWT(:,:,:),SINKQ(:,:,:),
     * REMVWT(:,:,:),REMVMS(:,:,:),NPNEW(:,:,:),SOURCEQ(:,:,:),
     * SGPTC(:,:,:),FRAC1(:,:,:),FRAC2(:,:,:),
cgzh srcfix2
     * SS_WT(:,:,:),SS_MS(:,:,:)
cgzh crit
      ALLOCATABLE CRITWT(:,:,:),CRITMS(:,:,:),NPSUM(:,:,:)
cgzh debug
      ALLOCATABLE SMCELL(:,:,:)
cgzh debug
      ALLOCATABLE QMNWSINK(:,:,:)
cgzh debug
      ALLOCATABLE QSS_SINK(:,:,:)
cgzh srcnew
      ALLOCATABLE SRCMS(:,:,:),SRCWT(:,:,:)
cgzh drycell
      ALLOCATABLE DRYWT(:,:,:),DRYMS(:,:,:),
     &  TEMPDRYWT(:,:,:),TEMPDRYMS(:,:,:)
cgzh ccbd
      ALLOCATABLE CHNGMASS(:,:,:)
cgzh zerocell
      ALLOCATABLE ZEROCELL(:,:,:)
cgzh srcfix
c     WTOUT: VOLUME OF PTS THAT LEFT A STRONG SRC and landed in a 
c       non-strong src
c     LIST: IP #'s of pts that left a strong src
c     DEFICIT: difference between WTOUT and SOURCEQ + SINKQ.  
c       IF DEFICIT<0, not enough left the src cell 
c       IF DEFICIT>0, too much left the src cell 
c     DEFMS: mass of culled wt from src cell pts
c
      ALLOCATABLE WTOUT(:,:,:),LIST(:,:),DEFICIT(:,:,:),
     * WTIN(:,:,:),LIST2(:,:),
     * DEFMS(:,:,:),AVC(:,:,:),AVR(:,:,:),AVL(:,:,:),LUMP(:,:,:),
     * DEFICIT2(:,:,:)
cgzh debug median
      ALLOCATABLE IVCOUNT(:,:,:)
C
      DOUBLE PRECISION DECAY
      DOUBLE PRECISION TINY,SMALL,PRECNO
      DOUBLE PRECISION OLDR,OLDC,OLDL,VCP,VRP,VLP,DISTC,DISTR,DISTL
cgzh debug oldc2
      DOUBLE PRECISION OLDC2
      DOUBLE PRECISION DBLTMP,TESTCK,DVDC,DVDL,DVDR
      DOUBLE PRECISION DCYFCT,DCYT,DCYT2
      DOUBLE PRECISION TOTFLUXIN,TOTFLUXOT
      PARAMETER (PERCNT=0.001)
      PARAMETER (TINY=1.D-20)
      PARAMETER (SMALL=1.D-4)
      PARAMETER (HUGE=1.E20)
C  CONSTANT FOR CHECK OF SINGLE PRECISION INTERVAL
c      PARAMETER (PRECNO=5.D-7)
cgzh debug  increasing for lejeune problem
      PARAMETER (PRECNO=5.D-5)
      DIMENSION
     *  PC(NPMAX),PR(NPMAX),PL(NPMAX),PCONC(NPMAX),IPTID(NPMAX),
     *  VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *  VL(NSCOL,NSROW,NSLAY+1),
     *  RF(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),CINFL(NCINFL),
     *  SRCFLO(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),CONC(NSCOL,NSROW,NSLAY),
     *  CAVG(NSCOL,NSROW,NSLAY),IBOUND(NCOL,NROW,NLAY),
     *  CNOLD(NSCOL,NSROW,NSLAY),
     *  SUMMASS(NSCOL,NSROW,NSLAY),
     *  NPCELL(NSCOL,NSROW,NSLAY),NPOLD(NSCOL,NSROW,NSLAY),
     *  LIMBO(NLIMBO),
     *  PNEWC(NEWPTS),PNEWR(NEWPTS),PNEWL(NEWPTS),
     *  IGENPT(NSCOL,NSROW,NSLAY)
cgzh varpt
      DIMENSION PCORIG(NPMAX),PRORIG(NPMAX),PLORIG(NPMAX),
     *  NPTCOLA(NSCOL,NSROW,NSLAY),NPTROWA(NSCOL,NSROW,NSLAY),
     *  NPTLAYA(NSCOL,NSROW,NSLAY)
      DIMENSION DKFO(NSCOL,NSROW,NSLAY),
     *    DKFS(NSCOL,NSROW,NSLAY)
CMOCWT
      DIMENSION PTWT(NPMAX),SBVL(6,NIUNIT)
      DIMENSION SUMWT(NSCOL,NSROW,NSLAY),SRCAVC(NSCOL,NSROW,NSLAY),
     *  BDYSRC(NSCOL,NSROW,NSLAY),BDYSNK(NSCOL,NSROW,NSLAY),
     *  BDYSOL(NSCOL,NSROW,NSLAY),SUMSGMS(NSCOL,NSROW,NSLAY),
     *  SUMSGWT(NSCOL,NSROW,NSLAY),RESIDWT(NSCOL,NSROW,NSLAY),
     *  SUMBFMS(NSCOL,NSROW,NSLAY),SUMBFWT(NSCOL,NSROW,NSLAY),
     *  RESIDC(NSCOL,NSROW,NSLAY),
     *  SGMSOUT(NSCOL,NSROW,NSLAY),SNKMSOUT(NSCOL,NSROW,NSLAY),
     *  CELVOL(NSCOL,NSROW,NSLAY),WTFAC(NSCOL,NSROW,NSLAY)
cgzh boundary face
      DIMENSION BDFACE(6,5)
cgzh mnw
      DIMENSION SRCMNW(NSCOL,NSROW,NSLAY),SNKMNW(NSCOL,NSROW,NSLAY),
     * IPTOBMNW(NSCOL,NSROW,NSLAY),IPTOBLST(4,NUMPTOB)
cgzh srcfix
      DIMENSION ISRCID(NSCOL,NSROW,NSLAY),SRCFAC(MAXSRC,6,3),
     * NPORIG(NSCOL,NSROW,NSLAY),HNEW(NCOL,NROW,NLAY)
cgzh srcfix2
      DIMENSION SS_SRC(NSCOL,NSROW,NSLAY),SS_SNK(NSCOL,NSROW,NSLAY)
      DIMENSION TOTFLUXIN(NSCOL,NSROW,NSLAY),
     & TOTFLUXOT(NSCOL,NSROW,NSLAY),ICONTRIBIN(NSCOL,NSROW,NSLAY,6),
     & ICONTRIBOT(NSCOL,NSROW,NSLAY,6)
cgzh ccbd
      DIMENSION CCBDY(NSCOL,NSROW,NSLAY)               
cgzh et
      DIMENSION EVTFLO(NSCOL,NSROW,2),IEVT(NROW,NCOL)             
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C1 ALLOCATE TEMPORARY ARRAYS
      ALLOCATE(TEMPMS(NSCOL,NSROW,NSLAY),TEMPWT(NSCOL,NSROW,NSLAY),
     * SINKQ(NSCOL,NSROW,NSLAY),REMVMS(NSCOL,NSROW,NSLAY),
     * REMVWT(NSCOL,NSROW,NSLAY),SOURCEQ(NSCOL,NSROW,NSLAY),
     * SS_WT(NSCOL,NSROW,NSLAY),SS_MS(NSCOL,NSROW,NSLAY),
     * SGPTC(NSCOL,NSROW,NSLAY),NPNEW(NSCOL,NSROW,NSLAY),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,1700)ISTAT
 1700   FORMAT(1X,'ALLOCATION OF GROUP 1 ARRAYS FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
C
      ALLOCATE (NPSUM(NSCOL,NSROW,NSLAY),CRITWT(NSCOL,NSROW,NSLAY),
     * CRITMS(NSCOL,NSROW,NSLAY),SRCMS(NSCOL,NSROW,NSLAY),
     * SRCWT(NSCOL,NSROW,NSLAY),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,1701)ISTAT
 1701   FORMAT(1X,'ALLOCATION OF GROUP 2 ARRAYS FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
C
      ALLOCATE(WTOUT(NSCOL,NSROW,NSLAY),LIST(MAXSRC,NPCELLMAX),
     * WTIN(NSCOL,NSROW,NSLAY),LIST2(MAXSRC,NPCELLMAX),
     * DEFICIT(NSCOL,NSROW,NSLAY),DEFMS(NSCOL,NSROW,NSLAY),
     * FRAC1(NSCOL,NSROW,NSLAY),FRAC2(NSCOL,NSROW,NSLAY),
     * AVC(NSCOL,NSROW,NSLAY),AVR(NSCOL,NSROW,NSLAY),
     * AVL(NSCOL,NSROW,NSLAY),LUMP(NSCOL,NSROW,NSLAY),
     * DEFICIT2(NSCOL,NSROW,NSLAY),
     * STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,1702)ISTAT
 1702   FORMAT(1X,'ALLOCATION OF GROUP 3 ARRAYS FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
cgzh debug
      write(iouts,*) 'debug: maxsrc,npcellmax',maxsrc,npcellmax
        CALL USTOP(' ')
      ENDIF
cgzh drycell
      ALLOCATE (DRYWT(NSCOL,NSROW,NSLAY),DRYMS(NSCOL,NSROW,NSLAY),
     &  TEMPDRYWT(NSCOL,NSROW,NSLAY),TEMPDRYMS(NSCOL,NSROW,NSLAY))
cgzh ccbd
      ALLOCATE (CHNGMASS(NSCOL,NSROW,NSLAY))
cgzh zerocell
      ALLOCATE (ZEROCELL(NSCOL,NSROW,NSLAY))
cgzh debug
      ALLOCATE(SMCELL(NSCOL,NSROW,NSLAY))
cgzh debug
      ALLOCATE(QMNWSINK(NSCOL,NSROW,NSLAY))
cgzh debug
      ALLOCATE(QSS_SINK(NSCOL,NSROW,NSLAY))
cgzh debug median
      ALLOCATE(IVCOUNT(NSCOL,NSROW,NSLAY))
c RBW debug
! It isn't clear that sumvol3 or wtof are used for anything.
      sumvol3 = 0
      wtof = 0
      wton=0.0
C RBW end debug            
C
cgzh debug
      sum2=0.
      sum3=0.
C2  COMPUTE DECAY TERMS
      DCYT=1.D0
      DCYT2=1.D0
      IF(DECAY.NE.0.D0) THEN
        DCYFCT=DBLE(TIMV)*DECAY
        DCYT=EXP(-DCYFCT)
        DCYT2=EXP(-DCYFCT*0.5D0)
      END IF
C
C  VCRIT IS MAXIMUM RELATIVE VELOCITY MULTIPLIED
C    BY CRITERION FOR IGNORING V CHANGE
      ARINV=1.0/(CDEL*RDEL)
      VCMAXR=VCMAX/CDEL
      VRMAXR=VRMAX/RDEL
      VLMAXR=1.0/TLMIN
C  USER CAN CHANGE 0.01 (1 PERCENT) CRITERION BELOW
      VCRIT=MAX(VCMAXR,VRMAXR,VLMAXR)*PERCNT
C  FLAG TO CHECK MULTIPLE CALLS TO SMOC5GP
      INPXFL=0
C  COMMENT FOLLOWING LINE FOR LESS DETAILED OUTPUT OF MOVE LOOP
      WRITE(IOUTS,670) NP,IMOV
C
C3D  ONE MOVE ONLY, RETURN HERE IF GENERATE PARTICLES
C
   10 NPTM=NP
C
C3  CLEAR SUMMASS & SUMWT, SET SOURCEQ
! SUMWT is not set to zero here. 
      DO 11 KS=1,NSLAY
      DO 11 IS=1,NSROW
      DO 11 JS=1,NSCOL
        SUMMASS(JS,IS,KS)=0.0
        NPNEW(JS,IS,KS)=0
        SUMSGMS(JS,IS,KS)=0.0
        SUMSGWT(JS,IS,KS)=0.0
        SUMBFMS(JS,IS,KS)=0.0
        SUMBFWT(JS,IS,KS)=0.0
        REMVWT(JS,IS,KS)=0.0
        REMVMS(JS,IS,KS)=0.0
        TEMPWT(JS,IS,KS)=0.0
        TEMPMS(JS,IS,KS)=0.0
        SINKQ(JS,IS,KS)=0.0
        SGMSOUT(JS,IS,KS)=0.0
        SNKMSOUT(JS,IS,KS)=0.0
        SRCMS(JS,IS,KS)=0.0
        SRCWT(JS,IS,KS)=0.0
cgzh srcfix
        WTOUT(JS,IS,KS)=0.0
        WTIN(JS,IS,KS)=0.0
        DEFICIT(JS,IS,KS)=0.0
        DEFMS(JS,IS,KS)=0.0
        AVC(JS,IS,KS)=0.0
        AVR(JS,IS,KS)=0.0
        AVL(JS,IS,KS)=0.0
        LUMP(JS,IS,KS)=0.0
cgzh debug median
        IVCOUNT(JS,IS,KS)=0.0
        DEFICIT2(JS,IS,KS)=0.0
cgzh srcfix2
        SS_WT(JS,IS,KS)=0.0
        SS_MS(JS,IS,KS)=0.0
   11 CONTINUE
cgzh srcfix
      LIST=0
      LIST2=0
cgzh drycell
      DRYMS=0.0   
      DRYWT=0.0     
      TEMPDRYMS=0.0   
      TEMPDRYWT=0.0     
      NDRYCELL=0
cgzh ccbd
      CHNGMASS=0.D0     
C4 CALCULATE SOURCEQ AND SINKQ
      DO 12 KS=1,NSLAY
        K = KS + ISLAY1 -1
      DO 12 IS=1,NSROW
      DO 12 JS=1,NSCOL
cgzh debug retard sourceq
c          SOURCEQ(JS,IS,KS)=TIMV*(SRCFLO(JS,IS,KS)+BDYSRC(JS,IS,KS))
c     &                            /RF(JS,IS,KS)
cgzh srcfix2
C  INCLUDE FLUX BETWEEN STRONG SOURCES (SS_SRC)
C  INMNW is for MNW1. INMNW2 is for MNW2
        SOURCEQ(JS,IS,KS)=TIMV*(SRCFLO(JS,IS,KS)+BDYSRC(JS,IS,KS)
     &                            +SS_SRC(JS,IS,KS))/RF(JS,IS,KS)
        IF((INMNW.GT.0).or.(INMNW2.GT.0))
     &    SOURCEQ(JS,IS,KS)=SOURCEQ(JS,IS,KS)+
     &                         (TIMV*SRCMNW(JS,IS,KS))/RF(JS,IS,KS)
cgzh debug
cgzh srcfix2 output
c      if(ss_src(js,is,ks).ne.0.0) then
c      write(iouts,*) 'SRC/SNK at cell js,is,ks',js,is,ks
c      write(iouts,*) 'SS_SRC*TIMV =',SS_SRC(JS,IS,KS)*TIMV
c      write(iouts,*) 'SS_SNK*TIMV =',SS_SNK(JS,IS,KS)*TIMV
c      if(js.eq.32.or.js.eq.62.or.js.eq.92) then
c      write(iouts,*) 'SOURCEQ =',SOURCEQ(JS,IS,KS)
c      write(iouts,*) 'SRCAVC =',SRCAVC(JS,IS,KS)
c      tm=SOURCEQ(JS,IS,KS)*SRCAVC(JS,IS,KS)
c      write(iouts,*) 'srcmass =',tm
c      end if
C
C  SINKQ IS VOLUME LEAVING CELL AS CALCULATED BY FD METHOD
C    PLUS RESIDUAL VOLUME (RESIDWT) DUE TO A SINK AT A CELL WITH
C    NO PARTICLES (CARRIED OVER FROM ANY PREVIOUS MOVES)
C
C    IF SINKQ = 0, DO NOT ADJUST PARTICLE WEIGHTS
C    IF SINKQ < 0, THEN VOLUME WILL BE REMOVED FROM PARTICLES
C    IF SINKQ > 0, THEN VOLUME WILL BE ADDED TO PARTICLES
C    IF NO PARTICLES, ADD TO RESIDUAL AND ATTEMPT TO APPLY IT NEXT MOVE
C
C  IF ET IS ACTIVE, GET ETFLO
cgzh debug
        if(imov.eq.7) then
          continue
        endif
        ETFLO=0.0
        IF(INEVT.GT.0) THEN
          KEV=EVTFLO(JS,IS,2)
          IF(KEV.EQ.K) ETFLO=EVTFLO(JS,IS,1)
        END IF
        SINKQ(JS,IS,KS)=TIMV*
     *    (SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS)+ETFLO)
     *                 +RESIDWT(JS,IS,KS)
   12 CONTINUE
C
C5    Handle MNW sources.
C Only do next few loops if MNW1 is on
      IF ((INMNW.GT.0).or.(INMNW2.GT.0)) THEN
C  FIRST MNW CELL LOOP
C  LOOP OVER CELLS TO:
C  CHECK FOR MNW SINKS WITH NOT ENOUGH PARTICLE WEIGHT AND CARRY OVER A RESIDUAL
        DO 160 KS=1,NSLAY
        DO 160 IS=1,NSROW
        DO 160 JS=1,NSCOL
C
C  DO CELLWISE CALCULATIONS FOR MNW SINKS
C
          QMNWSINK(JS,IS,KS)=SNKMNW(JS,IS,KS)*TIMV
          IF (QMNWSINK(JS,IS,KS).LT.0.0) THEN
C  IF NO PARTICLES IN CELL, CAN'T APPLY SINK VOLUME, SO STORE IN 
C    RESIDUAL WEIGHT AND TRY AGAIN NEXT MOVE
C  SET QMNWSINK=0.0 SINCE NO MASS REMOVED
C    (Use NPOLD as NPCELL is zeroed out before particle move)
            IF (NPOLD(JS,IS,KS).EQ.0) THEN
cgzh debug
c         write(iouts,*) '*WARNING* Residual due to NPCELL=0 at MNW sink'
c         write(iouts,*) 'Sink cell (js,is,ks):',js,is,ks
c         write(iouts,*) 'Residual= ',QMNWSINK(JS,IS,KS)
cgzh mnw residuals added in to residwt, it may have volume carried over
              RESIDWT(JS,IS,KS)=RESIDWT(JS,IS,KS)+QMNWSINK(JS,IS,KS)
              QMNWSINK(JS,IS,KS)=0.0
            ELSE
C  IF THERE ARE PARTICLES, BUT VOLUME ON ALL PARTICLES 
C    IS LESS THAN VOLUME TO BE REMOVED, REMOVE ALL VOLUME AND STORE RESIDUAL
              RESIDCL=SUMWT(JS,IS,KS)+QMNWSINK(JS,IS,KS)
              IF(RESIDCL.LT.0.0) THEN
cgzh debug
c         write(iouts,*) '***WARNING*** Residual due to not enough wt',
c     * ' on pts at MNW sink'
c         write(iouts,*) 'Sink cell (js,is,ks):',js,is,ks
c         write(iouts,*) 'Residual= ',RESIDCL
c         write(iouts,*) 'SUMWT=  ',SUMWT(js,is,ks)
c         write(iouts,*) 'QMNWSINK=  ',QMNWSINK(JS,IS,KS)
c         write(iouts,*) 'CONC=  ',CONC(JS,IS,KS)
                 RESIDWT(JS,IS,KS)=RESIDWT(JS,IS,KS)+RESIDCL
                     QMNWSINK(JS,IS,KS)=-SUMWT(JS,IS,KS)
              END IF
            END IF
          END IF
 160    CONTINUE
C
C
CGWT----PRINT PARTICLE OBSERVATION DATA
C       SEND IN QMNWSINK, WILL ONLY PRINT PTS IN MNWS HERE (CALL BELOW HANDLES
C       OTHER SINKS)
! RBW has NPCELL been updated at this point? Earlier NPOLD was used instead of NPCELL
        IF (INUNITPTOB.GT.0)
     *    CALL SPTOB5O(SUMTCH,IPTOBLST,IMOV,PC,PR,PL,NPMAX,NP,TIMV,
     *        NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,NUMPTOB,NUMPTOB_MNW,
     *        IPTOBMNW,NPCELL,PTWT,PCONC,SUMWT,SNKFLO,SNKMNW,1,MNWMAX, 
     *        MNW2, MNWNOD, NODTOT, NMNWVL)
C

C        ---MNW LOOP OVER PARTICLES: ADJUST WEIGHTS DUE TO SINKS---
        DO 390 IP=1,NPTM
C
          IF(PC(IP).EQ.0.0) GO TO 390
          J=PC(IP)+0.5
          JS=J-ISCOL1+1
          I=ABS(PR(IP))+0.5
          IS=I-ISROW1+1
          K=PL(IP)+0.5
          KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
          IF(IBOUND(J,I,K).EQ.0)       THEN
            WRITE(IOUTS,*) '**WARNING*** PARTICLE IN NO-FLOW CELL'
            WRITE(IOUTS,*) 'MNW LOOP 2'
            WRITE(IOUTS,*) 'J,I,K,IP,PTWT',J,I,K,IP,PTWT(IP)
C        STOP '**ERROR*** PARTICLE IN NO-FLOW CELL'        
            GO TO 390
          END IF
C
C  UPDATE WEIGHTS OF PTS. IN SINK CELLS
          IF (QMNWSINK(JS,IS,KS).LT.0.0) THEN
cgzh debug output
c      write(iouts,*) 'QMNWSINK js is ks=',QMNWSINK(JS,IS,KS),js,is,ks
C  COMPUTE FRACTION OF VOLUME TO BE REMOVED  
            IF(SUMWT(JS,IS,KS).GT.0.0) THEN
              FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*QMNWSINK(JS,IS,KS)
            ELSE
              FRACWT=0.0
            END IF
C  ADJUST VOLUME OF PARTICLE BY FRACTION OF VOLUME
C    NEW WEIGHT IS (OLD WEIGHT + FRACTION OF OLD WEIGHT)
            PTWT(IP)=PTWT(IP)+FRACWT
C   FIXED FOR ROUNDOFF ERROR CAUSING SMALL NEGATIVE WEIGHTS
            IF(PTWT(IP).LT.(-0.001)) THEN
              WRITE(IOUTS,*) 
     *          '***ERROR*** PTWT << 0, PTWT,IP',PTWT(ip),ip
            END IF
            if(ptwt(ip).lt.0.0) then
              write(iouts,*) 
     *           '**WARNING** PTWT < 0, ip,ptwt', ip,ptwt(ip)
              PTWT(IP)=0.0
            end if
C   TRACK VOLUME REMOVED FROM CELL
cgzh  don't need to track as this is calculated after pts are moved?
c        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+FRACWT
C  REMOVING WATER, ADJUST MASS BASED ON PT CONC, NO NEED TO ADJUST PCONC
c        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+FRACWT*PCONC(IP)
            SUM2=SUM2+FRACWT*PCONC(IP)*TIMV
            SUM3=SUM3+FRACWT*TIMV
          END IF
 390    CONTINUE
C END INMNW>0
      END IF
cgzh debug output
c      write(iouts,*) 'mass out mnws=', sum2
c      write(iouts,*) 'wt out mnws=', sum3
C     END MNW SINK LOOPS
C
C6 Handle strong source to strong source terms.
C  LOOP OVER CELLS TO APPLY SS_SNK, EXPLICIT "STRONG SOURCE-TO-STRONG SOURCE" TERM
C  (SS_SRC IS HANDLED IN SOURCEQ BELOW)
C  DOING THIS BEFORE ADVECTION WILL ENSURE THAT COLD IS USED FOR MASS TRANSFER
      IF(ISRCFIX.GT.0) THEN
C  LOOP OVER CELLS TO:
C  CHECK FOR SS-TO-SS SINKS WITH NOT ENOUGH PARTICLE WEIGHT AND CARRY OVER A RESIDUAL
        DO 170 KS=1,NSLAY
        DO 170 IS=1,NSROW
        DO 170 JS=1,NSCOL
C
C  DO CELLWISE CALCULATIONS FOR SS-TO-SS SINKS
C
          QSS_SINK(JS,IS,KS)=SS_SNK(JS,IS,KS)*TIMV
          IF (QSS_SINK(JS,IS,KS).LT.0.0) THEN
cgzh srcfix2 output
c      if(ss_snk(js,is,ks).ne.0.0) then
c      write(iouts,*) 'SS_SNK =',SS_SNK(JS,IS,KS)
c      end if
cgzh debug output
C  IF NO PARTICLES IN CELL, CAN'T APPLY SINK VOLUME, SO STORE IN 
C    RESIDUAL WEIGHT AND TRY AGAIN NEXT MOVE
C  SET QSS_SINK=0.0 SINCE NO MASS REMOVED
C    (Use NPOLD as NPCELL is zeroed out before particle move)
            IF (NPOLD(JS,IS,KS).EQ.0) THEN
cgzh debug
c         write(iouts,*) '*WARNING* Residual due to NPCELL=0 at SS_sink'
c         write(iouts,*) 'Sink cell (js,is,ks):',js,is,ks
c         write(iouts,*) 'Residual= ',QSS_SINK(JS,IS,KS)
cgzh SS-TO-SS residuals added in to residwt, it may have volume carried over
              RESIDWT(JS,IS,KS)=RESIDWT(JS,IS,KS)+QSS_SINK(JS,IS,KS)
              QSS_SINK(JS,IS,KS)=0.0
            ELSE
C  IF THERE ARE PARTICLES, BUT VOLUME ON ALL PARTICLES 
C    IS LESS THAN VOLUME TO BE REMOVED, REMOVE ALL VOLUME AND STORE RESIDUAL
              RESIDCL=SUMWT(JS,IS,KS)+QSS_SINK(JS,IS,KS)
              IF(RESIDCL.LT.0.0) THEN
cgzh debug
c         write(iouts,*) '***WARNING*** Residual due to not enough wt',
c     * ' on pts at SS-TO-SS sink'
c         write(iouts,*) 'Sink cell (js,is,ks):',js,is,ks
c         write(iouts,*) 'Residual= ',RESIDCL
c         write(iouts,*) 'SUMWT=  ',SUMWT(js,is,ks)
c         write(iouts,*) 'QSS_SINK=  ',QSS_SINK(JS,IS,KS)
c         write(iouts,*) 'CONC=  ',CONC(JS,IS,KS)
                 RESIDWT(JS,IS,KS)=RESIDWT(JS,IS,KS)+RESIDCL
                     QSS_SINK(JS,IS,KS)=-SUMWT(JS,IS,KS)
              END IF
            END IF
          END IF
 170    CONTINUE
C
C        ---SS-TO-SS LOOP OVER PARTICLES: ADJUST WEIGHTS DUE TO SINKS---
        DO 400 IP=1,NPTM
C
          NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
          IF(NEWC.LE.0.0D0) GO TO 400
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
          J=INT(NEWC+0.5D0)
          JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
          NEWR=PR(IP)
          IF(NEWR.LT.0.0D0) THEN
             NEWR=-NEWR
          END IF
          I=INT(NEWR+0.5D0)
          IS=I-ISROW1+1
          NEWL=PL(IP)
          K=INT(NEWL+0.5D0)
          KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
          IF(IBOUND(J,I,K).EQ.0)       THEN
            WRITE(IOUTS,*) '**WARNING*** PARTICLE IN NO-FLOW CELL'
            WRITE(IOUTS,*) 'SS-TO-SS LOOP 2'
            WRITE(IOUTS,*) 'J,I,K,IP,PTWT',J,I,K,IP,PTWT(IP)
C        STOP '**ERROR*** PARTICLE IN NO-FLOW CELL'        
            GO TO 400
          END IF
C
C  UPDATE WEIGHTS OF PTS. IN SINK CELLS
          IF (QSS_SINK(JS,IS,KS).LT.0.0) THEN
cgzh debug
c      if(js.eq.61) then
c      write(iouts,*) 'QSS_SINK<0 61',QSS_SINK(js,is,ks)
c      write(iouts,*) 'ptwt before',ptwt(ip)
c      endif
cgzh debug output
c      write(iouts,*) 'QSS_SINK js is ks=',QSS_SINK(JS,IS,KS),js,is,ks
C  COMPUTE FRACTION OF VOLUME TO BE REMOVED  
            IF(SUMWT(JS,IS,KS).GT.0.0) THEN
             FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*QSS_SINK(JS,IS,KS)
            ELSE
             FRACWT=0.0
            END IF
C  ADJUST VOLUME OF PARTICLE BY FRACTION OF VOLUME
C    NEW WEIGHT IS (OLD WEIGHT + FRACTION OF OLD WEIGHT)
            PTWT(IP)=PTWT(IP)+FRACWT
C   FIXED FOR ROUNDOFF ERROR CAUSING SMALL NEGATIVE WEIGHTS
            IF(PTWT(IP).LT.(-0.001)) THEN
              WRITE(IOUTS,*) 
     *          '***ERROR*** PTWT << 0, PTWT,IP',PTWT(ip),ip
            END IF
            if(ptwt(ip).lt.0.0) then
              write(iouts,*) 
     *          '**WARNING** PTWT < 0, ip,ptwt', ip,ptwt(ip)
              PTWT(IP)=0.0
            end if
cgzh debug
c      if(js.eq.61) then
c       write(iouts,*) 'ip,fracwt 61',ip,fracwt
c       write(iouts,*) 'ptwt after',ptwt(ip)
c      end if
C   TRACK VOLUME REMOVED FROM CELL
cgzh  don't need to track as this is calculated before pts are moved, i.e. before sumwt defined?
c        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+FRACWT
C  REMOVING WATER, ADJUST MASS BASED ON PT CONC, NO NEED TO ADJUST PCONC
c        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+FRACWT*PCONC(IP)
            SUM2=SUM2+FRACWT*PCONC(IP)*TIMV
            SUM3=SUM3+FRACWT*TIMV
          END IF
 400    CONTINUE
C END ISRCFIX>0
      END IF
C
C  DUE TO MNW AND SS_SNK LOOPS, MOVED SUMWT INIT HERE
      SUMWT=0.0
C
C  UNCOMMENT FOLLOWING LINES FOR MORE DETAILED OUTPUT OF MOVE LOOP
C     WRITE(IOUTS,*) ' START LOOP OVER PARTICLES'
C     WRITE(IOUTS,*) ' NP',NP
C        ---MOVE EACH PARTICLE---
cgzh debug output
c      write(iouts,*) 'np at start=',np
cgzh debug output
c        write(iouts,*) 'loop 0 pcorig(121)',pcorig(121)
c        write(iouts,*) 'loop 0 pc(121)',pc(121)
      CUMMASS3=0.0
cgzh debug zinnpt
c      if(imov.eq.68) then
c       write(iouts,*) 'Start SUMWT row 4=',SUMWT(4,1,1)
c      end if
c
C
C7     **** ADVECT PARTICLES ****
C
cgzh debug output
      DO 590 IP=1,NP
cgzh debug output
c      if(ip.eq.183) write(iouts,*) 'A: ptwt 183=',ptwt(183)
        OLDC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
        IF(OLDC.LE.0.0D0) GO TO 590
cgzh debug

        CUMMASS3=CUMMASS3+PTWT(IP)*PCONC(IP)
C     ***************************************************************
C           ---COMPUTE OLD LOCATION---
        J=INT(OLDC+0.5D0)
        JS=J-ISCOL1+1
C  IORIG SET TO 1 FOR PARTICLES ORIGINATING (CREATED) IN THIS CELL
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
        IORIG=0
        OLDR=PR(IP)
        IF(OLDR.LT.0.0D0) THEN
           OLDR=-OLDR
           IORIG=1
        END IF
        I=INT(OLDR+0.5D0)
        IS=I-ISROW1+1
        OLDL=PL(IP)
        K=INT(OLDL+0.5D0)
        KS=K-ISLAY1+1
cgzh boundary flux
C     ONLY CHECK PARTICLE LOCATION OUTSIDE GRID IF VELO IS NOT POINTING
C     OUT OF THAT FACE
        IPERROR=0
        IF(JS.LT.1.AND.VC(JS,IS,KS).GE.0.0) IPERROR=1
        IF(JS.GT.NSCOL.AND.VC(JS+1,IS,KS).LE.0.0) IPERROR=1
        IF(IS.LT.1.AND.VR(JS,IS,KS).GE.0.0) IPERROR=1
        IF(IS.GT.NSROW.AND.VR(JS,IS+1,KS).LE.0.0) IPERROR=1
        IF(KS.LT.1.AND.VL(JS,IS,KS).GE.0.0) IPERROR=1
        IF(KS.GT.NSLAY.AND.VL(JS,IS,KS+1).LE.0.0) IPERROR=1
        IF(IPERROR.EQ.1) THEN
           WRITE(IOUTS,*) ' IP,JS,IS,KS=',IP,JS,IS,KS
           WRITE(IOUTS,*) ' OLDC,OLDR,OLDL=',OLDC,OLDR,OLDL
           WRITE(IOUTS,*) ' PARTICLE ERROR MOVE, STOPPING'
           STOP ' PARTICLE ERROR IN MOVE'
        END IF
C1  USE ANALYTIC EXPRESSION WITHIN BLOCK FOR LINEAR V
        TSTEP=TIMV
C1  SAVE INITIAL LOCATION INDICES
        J1=J
        JS1=JS
        I1=I
        IS1=IS
        K1=K
        KS1=KS
cgzh debug
        if(ip.eq.5.and.imov.eq.181) then
          continue
        endif
C  INITIALIZE COUNT OF CELL BOUNDARY CROSSINGS
        NCROSS=0
C7-1 Handle particles in inactive and rewetted cells.
C  REMOVE PARTICLE IF LOCATED IN INACTIVE CELL
        IF(IBOUND(J,I,K).EQ.0) THEN
c        WRITE(IOUTS,*) 'PARTICLE REMOVED FROM NO-FLOW CELL: JS,IS,KS=',
c     * JS,IS,KS
cgzh drycell
           DRYMS(JS,IS,KS)=DRYMS(JS,IS,KS)+PTWT(IP)*PCONC(IP)        
           DRYWT(JS,IS,KS)=DRYWT(JS,IS,KS)+PTWT(IP)      
           NDRYCELL=NDRYCELL+1
           GO TO 545
        END IF
cgzh rewet
C  SKIP ADVECTION FOR PARTICLES ADDED AT REWETTED CELLS
C  CONC IS FLAGGED TO -888.0 FOR THESE PARTICLES 
C RBW, In a previous version, the concentration was set to -888.
C However, it was later changed so that the concentration was set to zero.
C IREWET(JS,IS,KS)=1 is the test
        IF(PCONC(IP).EQ.-888.0) THEN
          PCONC(IP)=0.0
          GO TO 590
        END IF
C7-2 Loop over partial time steps
C1  THIS IS BEGINNING OF LOOP FOR PARTIAL TIME STEPS
C1  CONST CONVERTS Q IN REAL UNITS TO RETARDED V IN RELATIVE UNITS
   32   CONST=ARINV/(RF(JS,IS,KS)*THCK(JS,IS,KS)*POR(JS,IS,KS))
C1  FIRST FOR C
        JSV1=JS
        JSV2=JS+1
cgzh intrpl 
C7-2a Calculate velocity in X direction (perpendicular to columns) and time to reach cell boundary. 
C     LINEAR INTERPOLATION
        ICSTVC = -1 
        IF(INTRPL.EQ.1) THEN
C  MOVTIM COMPUTES TIME TO MOVE TO CELL BOUNDARY IN GIVEN DIRECTION
          FAC=OLDC-J+0.5D0
          CALL MOVTIM(VC(JSV1,IS,KS),VC(JSV2,IS,KS),VCP,VCRIT,
     *            FAC,CONST,TIMC,DVDC,DISTC,ICSTVC,
     *            TINY,SMALL,HUGE)
C     BILINEAR INTERPOLATION
        ELSEIF (INTRPL.EQ.2) THEN
CBIL
CBIL  INTERPOLATE FLUX IN BOTH HORIZONTAL DIRECTIONS
CBIL  THIS IS STRAIGHT BILINEAR INTERPOLATION
CBIL  WILL NOT INTERPOLATE ADJACENT TO A NO-FLOW CELL
CBIL
C  FACR IS ZERO AT THE BOTTOM ROW BOUNDARY AND 1 AT THE TOP
          FACR=OLDR-I+0.5D0
C  IF EXACTLY AT THE MIDDLE IN THE ROW DIRECTION, USE THE VX FOR THIS CELL
          IF(FACR.EQ.0.5) THEN
             VC1=VC(JSV1,IS,KS)
             VC2=VC(JSV2,IS,KS)
C  IF IN THE UPPER PART OF THE CELL, TOWARDS INCREASING ROW NUMBER
C    THEN INTERPOLATE IN Y
          ELSE IF(FACR.GT.0.5) THEN
            IS2=IS+1
C  IF THE UPPER ROW IS OUTSIDE THE SUBGRID, THEN DON'T INTERPOLATE
            IF(IS2.GT.NSROW) THEN
              VC1=VC(JSV1,IS,KS)
              VC2=VC(JSV2,IS,KS)
            ELSE
C  IF THE ADJACENT CELL ABOVE IS NO-FLOW, DON'T INTERPOLATE
              IF(THCK(JS,IS2,KS).EQ.0.0) THEN
                 VC1=VC(JSV1,IS,KS)
                 VC2=VC(JSV2,IS,KS)
              ELSE
C  IF THIS IS THE FIRST COLUMN, INTERPOLATE FOR THE LEFT VALUE
                IF(JS.EQ.1) THEN
                  VC1=VC(JSV1,IS2,KS)+(1.0-(FACR-0.50))*
     *               (VC(JSV1,IS,KS)-VC(JSV1,IS2,KS))
                ELSE
C  CHECK TO SEE IF THE CELL ABOVE AND TO THE LEFT IS NO-FLOW
C   OR IF THE CELL TO THE LEFT IS NO-FLOW, IF SO, DON'T INTERPOLATE
                  IF(THCK(JS-1,IS2,KS).EQ.0.0.OR.
     *               THCK(JS-1,IS,KS).EQ.0.0) THEN
                     VC1=VC(JSV1,IS,KS)
                  ELSE
                     VC1=VC(JSV1,IS2,KS)+(1.0-(FACR-0.5))*
     *                  (VC(JSV1,IS,KS)-VC(JSV1,IS2,KS))
                  END IF
                END IF
C  IF THIS IS THE LAST COLUMN, INTERPOLATE FOR THE RIGHT VALUE
                IF(JS.EQ.NSCOL) THEN
                   VC2=VC(JSV2,IS2,KS)+(1.0-(FACR-0.5))*
     *              (VC(JSV2,IS,KS)-VC(JSV2,IS2,KS))
                ELSE
C  CHECK TO SEE IF THE CELL ABOVE AND TO THE RIGHT IS NO-FLOW
C    OR IF THE CELL TO THE RIGHT IS NO-FLOW, DON'T INTERPOLATE
                  IF(THCK(JS+1,IS2,KS).EQ.0.0.OR.
     *              THCK(JS+1,IS,KS).EQ.0.0) THEN
                    VC2=VC(JSV2,IS,KS)
                  ELSE
                    VC2=VC(JSV2,IS2,KS)+(1.0-(FACR-0.5))*
     *                 (VC(JSV2,IS,KS)-VC(JSV2,IS2,KS))
                  END IF
                END IF
              END IF
            END IF
C  IF THE PARTICLE IS IN THE LOWER HALF OF THE CELL, THEN INTERPOLATE
C   WITH THE VX FROM THE LOWER ROW
          ELSE
            IS2=IS-1
            IF(IS2.LT.1) THEN
              VC1=VC(JSV1,IS,KS)
              VC2=VC(JSV2,IS,KS)
            ELSE
              IF(THCK(JS,IS2,KS).EQ.0.0) THEN
                 VC1=VC(JSV1,IS,KS)
                 VC2=VC(JSV2,IS,KS)
              ELSE
                 IF(JS.EQ.1) THEN
                    VC1=VC(JSV1,IS2,KS)+(FACR+0.5)*
     *               (VC(JSV1,IS,KS)-VC(JSV1,IS2,KS))
                 ELSE
                    IF(THCK(JS-1,IS2,KS).EQ.0.0.OR.
     *               THCK(JS-1,IS,KS).EQ.0.0) THEN
                       VC1=VC(JSV1,IS,KS)
                    ELSE
                       VC1=VC(JSV1,IS2,KS)+(FACR+0.5)*
     *                  (VC(JSV1,IS,KS)-VC(JSV1,IS2,KS))
                    END IF
                 END IF
                 IF(JS.EQ.NSCOL) THEN
                    VC2=VC(JSV2,IS2,KS)+(FACR+0.5)*
     *               (VC(JSV2,IS,KS)-VC(JSV2,IS2,KS))
                 ELSE
                    IF(THCK(JS+1,IS2,KS).EQ.0.0.OR.
     *               THCK(JS+1,IS,KS).EQ.0.0) THEN
                       VC2=VC(JSV2,IS,KS)
                    ELSE
                       VC2=VC(JSV2,IS2,KS)+(FACR+0.5)*
     *                  (VC(JSV2,IS,KS)-VC(JSV2,IS2,KS))
                    END IF
                 END IF
              END IF
            END IF
          END IF
C
          FACC=OLDC-J+0.5D0
C  INTERPOLATE IN X DIRECTION
          VCP=(VC1+FACC*(VC2-VC1))*CONST
C  TIME TO REACH X BOUNDARY
          IF(VCP.EQ.0.0D0) THEN
             TIMC=HUGE
          ELSE IF(VCP.GT.0.0D0) THEN
             DISTC=J+0.5D0-OLDC
             IF(DISTC.LT.SMALL.AND.VC2.EQ.0.0) THEN
                VCP=0.0D0
                TIMC=HUGE
             ELSE      
                TIMC=DISTC/VCP
             ENDIF
          ELSE
             DISTC=J-0.5D0-OLDC
             IF(ABS(DISTC).LT.SMALL.AND.VC1.EQ.0.0) THEN
                VCP=0.0D0 
                TIMC=HUGE
             ELSE
                TIMC=DISTC/VCP
             ENDIF
          END IF
        END IF
C  NOW FOR R
cgzh intrpl 
C7-2b Calculate velocity in Y direction (perpendicular to rows) and time to reach cell boundary. 
C     LINEAR INTERPOLATION
        ICSTVR = -1
        IF(INTRPL.EQ.1) THEN
          ISV1=IS
          ISV2=IS+1
          FAC=OLDR-I+0.5D0
          CALL MOVTIM(VR(JS,ISV1,KS),VR(JS,ISV2,KS),VRP,VCRIT,
     *            FAC,CONST,TIMR,DVDR,DISTR,ICSTVR,
     *            TINY,SMALL,HUGE)
C     BILINEAR INTERPOLATION
        ELSEIF(INTRPL.EQ.2) THEN
          ISV1=IS
          ISV2=IS+1
          IF(FACC.EQ.0.5) THEN
             VR1=VR(JS,ISV1,KS)
             VR2=VR(JS,ISV2,KS)
          ELSE IF(FACC.GT.0.5) THEN
             JS2=JS+1
             IF(JS2.GT.NSCOL) THEN
                VR1=VR(JS,ISV1,KS)
                VR2=VR(JS,ISV2,KS)
             ELSE
                IF(THCK(JS2,IS,KS).EQ.0.0) THEN
                   VR1=VR(JS,ISV1,KS)
                   VR2=VR(JS,ISV2,KS)
                ELSE
                   IF(IS.EQ.1) THEN
                      VR1=VR(JS2,ISV1,KS)+(1.0-(FACC-0.5))*
     *                  (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
                   ELSE
                      IF(THCK(JS2,IS-1,KS).EQ.0.0.OR.
     *                  THCK(JS,IS-1,KS).EQ.0.0) THEN
                          VR1=VR(JS,ISV1,KS)
                      ELSE
                         VR1=VR(JS2,ISV1,KS)+(1.0-(FACC-0.5))*
     *                     (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
                      END IF
                   END IF
                   IF(IS.EQ.NSROW) THEN
                      VR2=VR(JS2,ISV2,KS)+(1.0-(FACC-0.5))*
     *               (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
                   ELSE
                      IF(THCK(JS2,IS+1,KS).EQ.0.0.OR.
     *                   THCK(JS,IS+1,KS).EQ.0.0) THEN
                            VR2=VR(JS,ISV2,KS)
                      ELSE
                         VR2=VR(JS2,ISV2,KS)+(1.0-(FACC-0.5))*
     *                     (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
                      END IF
                   END IF
                END IF
             END IF
          ELSE
             JS2=JS-1
             IF(JS2.LT.1) THEN
                VR1=VR(JS,ISV1,KS)
                VR2=VR(JS,ISV2,KS)
             ELSE
                IF(THCK(JS2,IS,KS).EQ.0.0) THEN
                   VR1=VR(JS,ISV1,KS)
                   VR2=VR(JS,ISV2,KS)
                ELSE
                   IF(IS.EQ.1) THEN
                      VR1=VR(JS2,ISV1,KS)+(FACC+0.5)*
     *                  (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
                   ELSE
                     IF(THCK(JS2,IS-1,KS).EQ.0.0.OR.
     *                 THCK(JS,IS-1,KS).EQ.0.0) THEN
                          VR1=VR(JS,ISV1,KS)
                     ELSE
                         VR1=VR(JS2,ISV1,KS)+(FACC+0.5)*
     *                     (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
                     END IF
                   END IF
                   IF(IS.EQ.NSROW) THEN
                      VR2=VR(JS2,ISV2,KS)+(FACC+0.5)*
     *                  (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
                   ELSE
                      IF(THCK(JS2,IS+1,KS).EQ.0.0.OR.
     *                  THCK(JS,IS+1,KS).EQ.0.0) THEN
                          VR2=VR(JS,ISV2,KS)
                      ELSE
                         VR2=VR(JS2,ISV2,KS)+(FACC+0.5)*
     *                     (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
                      END IF
                   END IF
                END IF
             END IF
          END IF
          VRP=(VR1+FACR*(VR2-VR1))*CONST
          IF(VRP.EQ.0.0D0) THEN
             TIMR=HUGE
          ELSE IF(VRP.GT.0.0D0) THEN
             DISTR=I+0.5D0-OLDR
             IF(DISTR.LT.SMALL.AND.VR2.EQ.0.0) THEN
                VRP=0.0D0
                TIMR=HUGE
             ELSE
                TIMR=DISTR/VRP
             ENDIF
          ELSE
             DISTR=I-0.5D0-OLDR
             IF(ABS(DISTR).LT.SMALL.AND.VR1.EQ.0.0) THEN
                VRP=0.0D0
                TIMR=HUGE
             ELSE
                TIMR=DISTR/VRP
             ENDIF
          END IF
        END IF
C  AND FOR L
C7-2c Calculate velocity in Z direction (perpendicular to layers) and time to reach cell boundary. 
        KSV1=KS
        KSV2=KS+1
        FAC=OLDL-K+0.5D0
        CALL MOVTIM(VL(JS,IS,KSV1),VL(JS,IS,KSV2),VLP,VCRIT,
     *            FAC,CONST,TIML,DVDL,DISTL,ICSTVL,
     *    TINY,SMALL,HUGE)
C
C7-2d Calculate length of partial time step.
C  TSTEP2 IS MINIMUM OF TOTAL TSTEP OR TIMES TO CELL BOUNDARIES
        if(ip.eq.208) then
          continue
        endif
        TSTEP2=MIN(TSTEP,TIMC,TIMR,TIML)
C
C7-2e Move particle in X direction
C  CHECK TO SEE IF REACHES C BOUNDARY
        TDIFF=TIMC-TSTEP2
        IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
C  IF TDIFF SMALL, ASSUME REACHES C BOUNDARY
cgzh debug output
c      if(ip.eq.12857) then
c       write(iouts,*) 'tdiff,small,',tdiff,small
c       write(iouts,*) 'oldc before',oldc
c      end if
cgzh debug end
        IF(TDIFF.LT.SMALL) THEN
C  MOVE TO BOUNDARY IF STEP TIME IS GREATER OR EQUAL TO
C    TIME TO C BOUNDARY
          IF(DISTC.GT.0.0D0) THEN
             J=J+1
             JS=JS+1
             OLDC=OLDC+DISTC+SMALL
          ELSE
             J=J-1
             JS=JS-1
             OLDC=OLDC+DISTC-SMALL
          END IF
        ELSE
          IF(INTRPL.EQ.1.AND.ICSTVC.EQ.0) THEN
C1  LINEAR VARIATION IN VC
             OLDC=VCP*(EXP(DVDC*DBLE(TSTEP2))-1.0D0)/DVDC+OLDC
          ELSE IF(ICSTVC.EQ.1.OR.INTRPL.EQ.2) THEN
C1  CONSTANT VC
             OLDC=VCP*DBLE(TSTEP2)+OLDC
          END IF
cgzh debug output
c      if(ip.eq.12857) then
c       write(iouts,*) 'oldc after',oldc
c      end if
cgzh debug end
C   CHECK FOR ERROR WHEN CONVERTING BACK TO SINGLE PRECISION
C   ONLY IN FORWARD DIRECTION
c orig line           IF (JS.EQ.NSCOL) THEN
C   CHECK TO SEE IF PARTICLE IS VERY CLOSE (VCLOSE) TO CELL BOUNDARY
          VCLOSE=(DBLE(J)+0.5-OLDC)
cgzh debug output
c      if(ip.eq.12857) then
c          write(iouts,*) 'vclose:',vclose
c      end if
C   VCLOSE>0 CHECK MEANS ONLY WHEN PARTICLE IS TO LEFT OF COLUMN BOUNDARY
          IF (VCLOSE.LT.SMALL.AND.VCLOSE.GT.0.0) THEN
C   OLDC,DBLTMP ARE DOUBLES; SNGTMP IS SINGLE PRECISION 
C   SNGTMP=OLDC SAVED IN SINGLE PRECISION
            SNGTMP=OLDC
C   DBLTMP=SNGTMP SAVED AS A DOUBLE
            DBLTMP=SNGTMP
C   TEST DIFFERENCE BETWEEN OLDC SAVED AS SINGLE VS DOUBLE
C   MULT BY APPROPRIATE PERCENTAGE IF ROUNDING ERROR  
            TESTCK=DBLTMP-OLDC
corig              IF (TESTCK.LT.SMALL.AND.TESTCK.GT.0.0)
corig     *            OLDC=DBLTMP*(1.D0-PRECNO)

cgzh  in some bigger proglems, testck was not computing correctly
cgzh  so we changed this to always move the particle back if gets
cgzh  "very close" to the cell boundary
            OLDC=DBLTMP-SMALL
          ENDIF
        END IF
C
C7-2f Move particle in Y direction
C  NOW FOR R
        TDIFF=TIMR-TSTEP2
        IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
cgzh debug output
c      if(ip.eq.112626) then
c       write(iouts,*) 'tdiff,small,',tdiff,small
c       write(iouts,*) 'oldr before',oldr
c      end if
cgzh debug end
        IF(TDIFF.LT.SMALL) THEN
           IF(DISTR.GT.0.0D0) THEN
              I=I+1
              IS=IS+1
              OLDR=OLDR+DISTR+SMALL
           ELSE
              I=I-1
              IS=IS-1
              OLDR=OLDR+DISTR-SMALL
           END IF
        ELSE
           IF(ICSTVR.EQ.0.AND.INTRPL.EQ.1) THEN
              OLDR=VRP*(EXP(DVDR*DBLE(TSTEP2))-1.D0)/DVDR+OLDR
           ELSE IF(ICSTVR.EQ.1.OR.INTRPL.EQ.2) THEN
              OLDR=VRP*DBLE(TSTEP2)+OLDR
           END IF
cgzh debug output
c      if(ip.eq.112626) then
c       write(iouts,*) 'oldr after',oldr
c      end if
cgzh debug end
C   CHECK FOR ERROR WHEN CONVERTING BACK TO SINGLE PRECISION
c orig line           IF (IS.EQ.NSROW) THEN
C   CHECK TO SEE IF PARTICLE IS VERY CLOSE (VCLOSE) TO CELL BOUNDARY
           VCLOSE=(DBLE(I)+0.5-OLDR)
cgzh debug output
c      if(ip.eq.112626) then
c          write(iouts,*) 'i,is,vclose:',i,is,vclose
c      end if
C   VCLOSE>0 CHECK MEANS ONLY WHEN PARTICLE IS ABOVE ROW BOUNDARY
           IF (VCLOSE.LT.SMALL.AND.VCLOSE.GT.0.0) THEN
              SNGTMP=OLDR
              DBLTMP=SNGTMP
              TESTCK=DBLTMP-OLDR
corig              IF (TESTCK.LT.SMALL.AND.TESTCK.GT.0.0)
corig     *         OLDR=DBLTMP*(1.D0-PRECNO)
cgzh  in some bigger proglems, testck was not computing correctly
cgzh  so we changed this to always move the particle back if gets
cgzh  "very close" to the cell boundary
              OLDR=DBLTMP-SMALL

cgzh debug output
c      if(ip.eq.112626) then
c           write(iouts,*) 'oldr after prec fix', oldr 
c      end if
           ENDIF
        END IF
C
C7-2g Move particle in Z direction
C  NOW L
C
        TDIFF=TIML-TSTEP2
        IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
        IF(TDIFF.LT.SMALL) THEN
           IF(DISTL.GT.0.0D0) THEN
              K=K+1
              KS=KS+1
              OLDL=OLDL+DISTL+SMALL
           ELSE
              K=K-1
              KS=KS-1
              OLDL=OLDL+DISTL-SMALL
           END IF
        ELSE
           IF(ICSTVL.EQ.0) THEN
              OLDL=VLP*(EXP(DVDL*DBLE(TSTEP2))-1.0D0)/DVDL+OLDL
           ELSE IF(ICSTVL.EQ.1) THEN
              OLDL=VLP*DBLE(TSTEP2)+OLDL
           END IF
C   CHECK FOR ERROR WHEN CONVERTING BACK TO SINGLE PRECISION
c orig           IF (KS.EQ.NSLAY) THEN
C   CHECK TO SEE IF PARTICLE IS VERY CLOSE (VCLOSE) TO CELL BOUNDARY
           VCLOSE=(DBLE(K)+0.5-OLDL)
C   VCLOSE>0 CHECK MEANS ONLY WHEN PARTICLE IS ABOVE LAYER BOUNDARY
           IF (VCLOSE.LT.SMALL.AND.VCLOSE.GT.0.0) THEN
              SNGTMP=OLDL
              DBLTMP=SNGTMP
              TESTCK=DBLTMP-OLDL
corig              IF (TESTCK.LT.SMALL.AND.TESTCK.GT.0.0)
corig     *         OLDL=DBLTMP*(1.D0-PRECNO)
cgzh  in some bigger proglems, testck was not computing correctly
cgzh  so we changed this to always move the particle back if gets
cgzh  "very close" to the cell boundary
              OLDL=DBLTMP-SMALL
           ENDIF
        END IF
C7-2h  INCREASE COUNT OF CELL BOUNDARY CROSSINGS
cgzh old check      IF(TSTEP2.LT.TSTEP) NCROSS=NCROSS+1
cgzh instead of checking vs time step, just check change in cell location
        IF(J.NE.J1.OR.I.NE.I1.OR.K.NE.K1) NCROSS=NCROSS+1
C
cgzh srcfix2
C
C7-2i Handle particles that leave grid or subgrid or that go between 
C  strong sources when volume balancing is used.
C  ISRCFIX = 1 means to use volume balancing in source cells so that 
C  the combined volumes of the particles in a cell equals 
C  the volume of the cell.
C  IF SRCFIX IS ON, REMOVE PARTICLES THAT LEAVE GRID, SUBGRID AND PTS FROM SS to SS
C  THIS FLUX IS ACCOUNTED FOR USING SS_SNK AND SS_SRC
cgzh debug output
cgzh debug output
        if(imov.eq.1.and.js1.eq.120.and.is1.eq.1) then
          continue
        end if
c      if(ip.eq.183) write(iouts,*) 'B: ptwt 183=',ptwt(183)
        IF(ISRCFIX.GT.0) THEN
C  1 of 3) CHECK TO SEE IF LEFT FLOW DOMAIN OR ENTERED NO-FLOW CELL
C  THIS SHOULD ONLY HAPPEN IF BFLX PACKAGE IS ON
          ILEFT=0
          IF(J.LT.1.OR.J.GT.NCOL.OR.
     *       I.LT.1.OR.I.GT.NROW.OR.
     *       K.LT.1.OR.K.GT.NLAY) THEN
             ILEFT=1
          ELSE IF(IBOUND(J,I,K).EQ.0) THEN
             ILEFT=1
          END IF
          IF(ILEFT.EQ.1) THEN
C Track this weight in SS_WT for cell, to be added back into strong source cell later
            SS_WT(JS1,IS1,KS1)=SS_WT(JS1,IS1,KS1)+PTWT(IP)
            SS_MS(JS1,IS1,KS1)=SS_MS(JS1,IS1,KS1)+PTWT(IP)*PCONC(IP)
cgzh debug
c      if(js1.eq.61) then
c       write(iouts,*) 'A:ss_wt updated in 61: ss_wt=',SS_WT(JS1,IS1,KS1)
c      end if
C ZERO OUT PARTICLE
            GO TO 545
          END IF
C  2 of 3) CHECK TO SEE IF MOVED FROM STRONG SOURCE TO OUTSIDE SUBGRID
          IF(IGENPT(JS1,IS1,KS1).EQ.1.AND.
     *       (JS.LT.1.OR.JS.GT.NSCOL.OR.
     *        IS.LT.1.OR.IS.GT.NSROW.OR.
     *        KS.LT.1.OR.KS.GT.NSLAY)) THEN
C  NLOC=1 IF PARTICLE LEFT SUBGRID
C            TRACK SOLUTE MASS & VOLUME THAT LEAVES SUBGRID ON PTS.
                NLOC=1
C Track this weight in SS_WT for cell, to be added back into strong source cell later
                SS_WT(JS1,IS1,KS1)=SS_WT(JS1,IS1,KS1)+PTWT(IP)
                SS_MS(JS1,IS1,KS1)=SS_MS(JS1,IS1,KS1)+PTWT(IP)*PCONC(IP)
cgzh debug output
c         write(iouts,*) 'Particle left subgrid SS: js,is,weight,conc= ',
c     * js1,is1,PTWT(IP),PCONC(IP)
C ZERO OUT PARTICLE
                GO TO 545
          END IF
C  3 of 3) OR, IF PARTICLE MOVED FROM ONE STRONG SOURCE TO ANOTHER, REMOVE IT
C  THIS FLUX IS ALSO ACCOUNTED FOR USING SS_SNK AND SS_SRC
cgzh debug
          IF (NCROSS.GT.0) THEN
cgzh only if particle is still in subgrid
             IF(JS.GE.1.AND.JS.LE.NSCOL.AND.IS.GE.1.AND.IS.LE.NSROW.
     &           AND.KS.GE.1.AND.KS.LE.NSLAY) THEN  
               IF(IGENPT(JS1,IS1,KS1).EQ.1.
     &             AND.IGENPT(JS,IS,KS).EQ.1) THEN
cgzh debug output
c      if(imov.eq.257.or.imov.eq.96) then
c      write(iouts,*) 'Particle moved between strong sources'
c      write(iouts,*) 'ip,ptwt,pconc',ip,ptwt(ip),pconc(ip)
c      write(iouts,*) 'orig: js1,is1,ks1',js1,is1,ks1
c      write(iouts,*) 'new:  js,is,ks',js,is,ks
c      write(iouts,*) 'mass=',PTWT(IP)*PCONC(IP)
c      end if
C
C Track this weight in SS_WT for cell, to be added back into strong source cell later
                 SS_WT(JS1,IS1,KS1)=SS_WT(JS1,IS1,KS1)+PTWT(IP)
                 SS_MS(JS1,IS1,KS1)=
     *             SS_MS(JS1,IS1,KS1)+PTWT(IP)*PCONC(IP)
cgzh debug
c      if(js1.eq.61) then
c       write(iouts,*) 'B:ss_wt updated in 61: ss_wt=',SS_WT(JS1,IS1,KS1)
c       write(iouts,*) 'B:ip,ptwt=',ip,ptwt(IP)
c       write(iouts,*) 'B:js,js1=',js,js1
c      end if
C ZERO OUT PARTICLE
                 GO TO 545
               END IF
             END IF
          END IF
        END IF
C
C7-2j  Handle particles that leave grid or subgrid.
C  CHECK TO SEE IF LEFT FLOW DOMAIN OR ENTERED NO-FLOW CELL
C  THIS SHOULD ONLY HAPPEN IF BFLX PACKAGE IS ON
        ILEFT=0
        IF(J.LT.1.OR.J.GT.NCOL.OR.
     *     I.LT.1.OR.I.GT.NROW.OR.
     *     K.LT.1.OR.K.GT.NLAY) THEN
           ILEFT=1
        ELSE IF(IBOUND(J,I,K).EQ.0) THEN
           ILEFT=1
        END IF
        IF(ILEFT.EQ.1) THEN
cgzh check BFLX and if not on, write warning
          IF(INBFLX.EQ.0) THEN
            write(iouts,*) 
     *        '***WARNING*** Particle crossed a no-flow boundary'
            write(iouts,*) 
     *        'This should not happen with BFLX package off'
            write(iouts,*) 'Originating cell (js,is,ks):',js1,is1,ks1
            write(iouts,*) 'Current cell (js,is,ks):',js,is,ks
          ENDIF

cgzh debug output
c      write(IOUTS,*) 'PARTICLE IN CELL WITH IBOUND=0; NCROSS=',NCROSS
c      write(IOUTS,*) 'CELL JS,IS,KS',JS,IS,KS
c      write(IOUTS,*) 'ptwt,pconc',ptwt(ip),pconc(ip)
        
c if particle crossed two or more cell boundaries, 
c   then attribute weight leaving to intermediate cell 
c
cgzh debug put in separate itemization: assume BFLX package is on
cgzh 8/23/06 scratch that, now put into SNKMSOUT, will be divvyed up
cgzh  correctly in MB routine still, adjust SINKQ accordingly for
cgzh  this removed weight
          IF (NCROSS.GT.1) THEN
c            SUMSGMS(JS2,IS2,KS2)=SUMSGMS(JS2,IS2,KS2)-PCONC(IP)*PTWT(IP)
c              SUMSGWT(JS2,IS2,KS2)=SUMSGWT(JS2,IS2,KS2)-PTWT(IP)
            SNKMSOUT(JS2,IS2,KS2)=
     *        SNKMSOUT(JS2,IS2,KS2)-PCONC(IP)*PTWT(IP)
            SINKQ(JS2,IS2,KS2)=SINKQ(JS2,IS2,KS2)-PTWT(IP)
cgzh debug
c      write(iouts,*) 'bfms in j,k',SNKMSOUT(JS2,IS2,KS2),js2,ks2
          ELSE
c otherwise, attribute weight leaving to originating cell
            SNKMSOUT(JS1,IS1,KS1)=
     *        SNKMSOUT(JS1,IS1,KS1)-PCONC(IP)*PTWT(IP)
            SINKQ(JS1,IS1,KS1)=SINKQ(JS1,IS1,KS1)-PTWT(IP)
cgzh debug
c      write(iouts,*) 'bfms in j,k',SNKMSOUT(JS1,IS1,KS1),js1,ks1
c            SUMSGMS(JS1,IS1,KS1)=SUMSGMS(JS1,IS1,KS1)-PCONC(IP)*PTWT(IP)
c            SUMSGWT(JS1,IS1,KS1)=SUMSGWT(JS1,IS1,KS1)-PTWT(IP)
          END IF   
C ZERO OUT PARTICLE
          GO TO 545
        
C  CHECK TO SEE IF STILL IN ACTIVE CELL OF TRANSPORT SUBGRID
        ELSE IF(JS.LT.1.OR.JS.GT.NSCOL.OR.
     *    IS.LT.1.OR.IS.GT.NSROW.OR.
     *    KS.LT.1.OR.KS.GT.NSLAY) THEN
C  NLOC=1 IF PARTICLE LEFT SUBGRID
C            TRACK SOLUTE MASS & VOLUME THAT LEAVES SUBGRID ON PTS.
            NLOC=1
cgzh debug from angle2k bug
c if particle crossed two or more cell boundaries, and
c   if originating cell is not a sink, 
c   then attribute weight leaving to intermediate cell 
cgzh (we assume *it* is a sink)
            IF (NCROSS.GT.1.AND.
     *         (SNKFLO(JS1,IS1,KS1).EQ.0.0.AND.
     *          BDYSNK(JS1,IS1,KS1).EQ.0.0)) THEN
                  SUMSGMS(JS2,IS2,KS2)=
     *              SUMSGMS(JS2,IS2,KS2)-PCONC(IP)*PTWT(IP)
                  SUMSGWT(JS2,IS2,KS2)=SUMSGWT(JS2,IS2,KS2)-PTWT(IP)
cgzh debug output
c           write(iouts,*) 'Particle left subgrid, attributed to ',
c     * 'intermediate cell: js2,is2,weight= ',
c     * js2,is2,PTWT(IP)
            ELSE
c otherwise, attribute weight leaving to originating cell
              SUMSGMS(JS1,IS1,KS1)=
     *          SUMSGMS(JS1,IS1,KS1)-PCONC(IP)*PTWT(IP)
              SUMSGWT(JS1,IS1,KS1)=SUMSGWT(JS1,IS1,KS1)-PTWT(IP)
cgzh debug output
c          if(js1.eq.5.and.is1.eq.1) then
c           write(iouts,*) 'Particle left subgrid: js,is,weight,conc= ',
c     * js1,is1,PTWT(IP),PCONC(IP)
c          end if
            END IF
C ZERO OUT PARTICLE
            GO TO 545
        END IF
C
C7-2k  NOW CHECK TO SEE IF MOVE COMPLETED
C
        IF(TSTEP2.LT.TSTEP) THEN
C
C  REDUCE STEP SIZE AND TAKE A NEW STEP
          TSTEP=TSTEP-TSTEP2
C  SAVE INTERMEDIATE LOCATION IN CASE NEEDED FOR SUMSGWT
          JS2=JS
          IS2=IS
          KS2=KS
C1  RETURN TO TOP OF LOOP FOR NEXT STEP IF STILL IN TRANSPORT SUBGRID
          GO TO 32
        ENDIF
C1  END
C7-3 Store particle location
C  END OF MOVE FOR THIS PARTICLE
        PC(IP)=OLDC
        PR(IP)=OLDR
        PL(IP)=OLDL
cgzh debug
c      if(ip.eq.110593) then
c       write(iouts,'(a,i4,f16.9)') 'move:js,pr of 110593',js,pr(ip)
c      end if
cgzh debug end
C  NLOC=0 SIGNIFIES PARTICLE IS WITHIN SUBGRID
        NLOC=0
C     ***************************************************************
C      ---SUM WEIGHTED CONCENTRATIONS & WEIGHTS, & COUNT PARTICLES---
C           ---DECAY PARTICLES---
C7-4 handle decay of particle concentrations.
        IF(DECAY.NE.0.D0) THEN
C ONLY DECAY PARTICLES WITH POSITIVE CONCENTRATION TO AVOID "CREATING"
C MASS ON PTS 
c        IF(PCONC(IP).GT.0.0) THEN
            PCOLD=PCONC(IP)
            PCONC(IP)=DBLE(PCONC(IP))*DCYT
            PCNCNC=PCOLD-PCONC(IP)
cgzh debug output roland
c      if(imov.eq.15.and.ip.eq.761) then 
c        write(*,*) 'after decay 2,pcncnc',pcncnc
c        write(*,*) 'after decay 2,ptwt',ptwt(ip)
c        write(*,*) 'after decay 2,pcold,pconc',pcold,pconc(ip)
c      end if
cgzh debug output roland
            IF(PCNCNC.GT.0.0) THEN
cgzh debug add RF to decayed mass
              SBVL(4,3)=SBVL(4,3)-(PCNCNC)*PTWT(IP)*RF(JS,IS,KS)
            ELSE
              SBVL(3,3)=SBVL(3,3)-(PCNCNC)*PTWT(IP)*RF(JS,IS,KS)
            END IF
c        END IF
        ELSE IF(INDK.GT.0) THEN
          IF(IDKFO.EQ.1.OR.IDKFS.EQ.1) THEN
            CALL DK6DK(DCYT,DKFO,DKFS,RF,TIMV,
     *          IDKFO,IDKFS,
     *          JS,IS,KS,NSCOL,NSROW,NSLAY)
            PCONC(IP)=PCONC(IP)*DCYT
          END IF
        END IF
C7-5 Update number of particles, volume, and concentration in cell, 
        NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)+1
        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+PTWT(IP)
cgzh debug output
c      if(js.eq.61) then
c      if(ip.eq.180) then
c       write(iouts,*) 'pt in 61 in 590, ip, ptwt',ip,ptwt(ip)
c      end if
c      SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)

        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+PCONC(IP)*PTWT(IP)
cgzh debug output
c      if(js.eq.32.or.js.eq.62.or.js.eq.92) then
c       write(iouts,*) '0: js, summass',js,summass(js,is,ks)
c      end if
cgzh debug output
cgzh debug zinnpt
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.17) then
c       write(iouts,*) 'pt in 17: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.18) then
c       write(iouts,*) 'pt in 18: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
c      if(summass(js,is,ks).lt.0) then
c       write(iouts,*) 'js,is,ks',js,is,ks
c       write(iouts,*) 'ip,ptwt,pconc'
c       write(iouts,*) ip,ptwt(IP),pconc(ip)
c      stop
c      endif
c      if(js.eq.140.and.ks.eq.2.and.imov.gt.5670) then
c      if(js.eq.26.and.ks.eq.3) then
c       if(pconc(ip).gt.1e-5.or.pconc(ip).lt.-1e-5) then
c       write(iouts,*) 'move a: ip,ptwt,pconc'
c       write(iouts,*) ip,ptwt(IP),pconc(ip)
c       write(iouts,*) 'pc,pr,pl',pc(ip),pr(ip),pl(ip)
c       end if
c      end if
c      if(imov.eq.3.and.is.eq.2.and.is1.eq.1) then
c      continue
c      endif
cgzh debug end
C     ****************************************************************
C7-6        ---CHECK FOR CHANGE IN CELL LOCATION---
        IF(J.EQ.J1.AND.I.EQ.I1.AND.K.EQ.K1) THEN
C  SET R COORDINATE NEGATIVE IF STILL IN ORIGINAL CELL
           IF(IORIG.EQ.1) PR(IP)=-PR(IP)
C  DONE WITH THIS PARTICLE
           GO TO 590
        END IF
C
C7-7  SUM WEIGHTS AND TRACK IPs OF PTS ENTERING SOURCE CELL
C  ONLY TRACK PARTICLES IF THEY ORIGINATE IN A NON-STRONG SOURCE
C  JS,IS,KS here is new location, JS1,IS1,KS1 is originating cell
        IF(IGENPT(JS,IS,KS).EQ.1.AND.IGENPT(JS1,IS1,KS1).EQ.0) THEN
cgzh if there is no flow into the cell, this will not work, so 
cgzh put particle back into originating cell
cgzh (cf 3ddisp problem at move 181, flow at 45 degrees moves pt
cgzh  into cell across corner)
          IF(ISRCFIX.GT.0) THEN
            IF(TOTFLUXIN(JS,IS,KS).EQ.0.0) THEN
cgzh debug output
c      write(iouts,*) 'TOTFLUXIN = 0, WTIN >0, kicking out particle'
              SS_WT(JS1,IS1,KS1)=SS_WT(JS1,IS1,KS1)+PTWT(IP)
              SS_MS(JS1,IS1,KS1)=SS_MS(JS1,IS1,KS1)+PTWT(IP)*PCONC(IP)
c          SUMWT(JS1,IS1,KS1)=SUMWT(JS1,IS1,KS1)+PTWT(IP)
c          SUMMASS(JS1,IS1,KS1)=SUMMASS(JS1,IS1,KS1)+PCONC(IP)*PTWT(IP)
c          NPCELL(JS1,IS1,KS1)=NPCELL(JS1,IS1,KS1)+1
              SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)-PTWT(IP)
              SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)-PCONC(IP)*PTWT(IP)
              NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)-1
C ZERO OUT PARTICLE
              GO TO 545
            END IF
C
C7-8 Update weight in for new cell
            WTIN(JS,IS,KS)=WTIN(JS,IS,KS)+PTWT(IP)
C7-9 Store particle number in LIST2
C  FOR EACH STRONG SOURCE, ISRCID HAS SEQUENTIAL NONZERO INTEGERS 
            DO 6444 LLST=1,NPCELLMAX
              IF(LIST2(ISRCID(JS,IS,KS),LLST).EQ.0) THEN
                LIST2(ISRCID(JS,IS,KS),LLST)=IP
                GO TO 6445
              END IF
cgzh debug output if no more room, quit
              if(LLST.EQ.NPCELLMAX) then
                write(iouts,*) 
     *            '***ERROR*** NPCELLMAX EXCEEDED in 6445 loop'
                write(iouts,*) 'js,is,ks,npcellmax',js,is,ks,npcellmax
                STOP 'ERROR IN MOVEWT: NPCELLMAX'
              end if
 6444       CONTINUE        
          END IF
 6445     CONTINUE        
        END IF
C           ---CHECK FOR STRONG FLUID SOURCE AT OLD LOCATION---
C  IGENPT=1: STRONG SOURCE, ADD A PARTICLE
C  NO SINK/SOURCES IF ICONLY=1
cgzh debug
c      if(imov.eq.1460.and.ip.eq.121) then
c       write(iouts,*) 'stop the debugger'
c      endif
C7-10
        IF(ICONLY.NE.1) THEN
CMOCWT  IGENLK NOT NEEDED FOR WEIGHTED VERSION; ALL CELLS AUTOMATICALLY
C       DETERMINED USING GENCRIT 
C  JS,IS,KS here is new location, JS1,IS1,KS1 is originating cell
           IF(IGENPT(JS1,IS1,KS1).EQ.1) THEN
c     *     OR.IGENLK(J1,I1,K1).NE.0) THEN
               GO TO 350
           END IF
        END IF
C  IGENPT=0: WEAK SOURCE, DO NOT ADD A PARTICLE, SOURCE WILL
C    BE ACCOUNTED FOR IN NEXT LOOP OVER PARTICLES (1590 loop)
cgzh nloc should always be = 0 here, change to 590
corig      GO TO 540
        GO TO 590
  350   CONTINUE
cgzh srcfix
C
C7-11 if volume balancing is used and the new location is not a strong
C     source, Update WTOUT and store the particle number in LIST
        IF(ISRCFIX.GT.0) THEN
C  SUM WEIGHTS AND TRACK IPs OF PTS LEAVING SOURCE CELL
C  ONLY TRACK PARTICLES IF THEY LAND IN A NON-STRONG SOURCE
C  JS,IS,KS here is new location, JS1,IS1,KS1 is originating cell
          IF(IGENPT(JS,IS,KS).EQ.0) THEN
            WTOUT(JS1,IS1,KS1)=WTOUT(JS1,IS1,KS1)+PTWT(IP)
C  FOR EACH STRONG SOURCE, ISRCID HAS SEQUENTIAL NONZERO INTEGERS 
            DO 4444 LLST=1,NPCELLMAX
              IF(LIST(ISRCID(JS1,IS1,KS1),LLST).EQ.0) THEN
                LIST(ISRCID(JS1,IS1,KS1),LLST)=IP
cgzh debug output
c          write(iouts,*) 'IP= ',IP,' added to list'
                GO TO 4445
              END IF
cgzh debug output if no more room, quit
              if(LLST.EQ.NPCELLMAX) then
                write(iouts,*) 
     *            '***ERROR*** NPCELLMAX EXCEEDED in 4445 loop'
                write(iouts,*) 'js,is,ks,npcellmax',js,is,ks,npcellmax
                STOP 'ERROR IN MOVEWT: NPCELLMAX'
              end if
 4444       CONTINUE        
 4445       CONTINUE        
          END IF
        END IF
cgzh srcfix end
C7-12  ALSO DO NOT CREATE A PARTICLE IF THE PARTICLE THAT LEFT STRONG SOURCE CELL
C   WAS NOT CREATED (I.E. DID NOT ORIGINATE) IN THE SOURCE CELL
corig  350 IF(IORIG.EQ.0) GO TO 540
        IF(IORIG.EQ.0) GO TO 590
C     ****************************************************************
C7-13           ---CREATE NEW PARTICLES---
C     ****************************************************************
        IF(NPTM.EQ.NPMAX) THEN
C7-14        ---RESTART MOVE IF PT. LIMIT EXCEEDED---
          WRITE(IOUTS,700) IMOV,IP
C
C DON'T CALL SMOC5GP TWICE in ONE MOVE
          IF(INPXFL.EQ.1) THEN
            WRITE(IOUTS,710) IMOV, IP
            STOP
          ENDIF
cgzh debug
c this is not working yet with weighted particles (regen), so stop
          write(*,*)
          STOP 'NPMAX EXCEEDED IN MOVEWT'
C
cgzh varpt  this not set up for variable point distributions
          CALL SMOC5GP(PC,PR,PL,PCONC,
     *                   CONC,IPTID,NPCELL,
     *     IBOUND,PNEWC,PNEWR,PNEWL,LIMBO,
     *     NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *     NEWPTS,NPMAX,NLIMBO,
     *     IOUTS,NP,WTFAC)
C
CMOCWT
cgzh debug  need to go over this??  not conserving mass
CGWT----DETERMINE INITIAL PARTICLE WEIGHTS AFTER REGEN
           CALL PTWT1INITWT(IBOUND,PC,PR,PL,NPCELL,CELVOL,PTWT,SUMWT,
     *       NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NPMAX,NP,IOUTS)
C
C FLAG FIRST CALL OF SMOC5GP
          INPXFL=1
C
          DO 610 KS=1,NSLAY
          DO 610 IS=1,NSROW
          DO 610 JS=1,NSCOL
            NPCELL(JS,IS,KS)=0
  610     CONTINUE
C RBW     Does NPOLD need to be updated here? It is used in C5 and C6
C RBW     consider calling MOC5AP here with TIMV = 0.
          GO TO 10
        ENDIF
C     ****************************************************************
C7-15           ---GENERATE NEW PARTICLE---
C ONLY AT STRONG FLUID SOURCE CELLS
C  INCREASE NUMBER OF PARTICLES
        NPTM=NPTM+1
C  SET NEW PARTICLE ID
        IPN=NPTM
C  GET RELATIVE LOCATION OF OLD PARTICLE IN CELL
cgzh debug output
c      write(iouts,*) 'New pt ipn at js is ks:',ipn,js1,is1,ks1
cgzh old particle ID is IP
C
c above will be user input eventually
cgzh debug  bug fix 6/11/09 (this did say EQ.1)
        IF(IRAND.GT.0) THEN
          DO 78 IDIR=1,3
   68       CONTINUE
            RANDLOC=RAN(ISEEDPT)
C DON'T USE RANDOM NUMBERS THAT MAY RESULT IN ROUNDOFF ERRORS
            IF(RANDLOC.LT.1.E-5.OR.RANDLOC.GT.0.9999) GO TO 68
cgzh debug check dimensionality here?  only do for dim-width >1
            IF(IDIR.EQ.1) PC(IPN)=J1+RANDLOC-0.5
            IF(IDIR.EQ.2) PR(IPN)=-I1-RANDLOC+0.5
            IF(IDIR.EQ.3) PL(IPN)=K1+RANDLOC-0.5
   78     CONTINUE
        ELSE
          IF(INIPDL.EQ.0.AND.INIPDA.EQ.0) THEN
            ITEM=IPTID(IP)
            PC(IPN)=J1+PNEWC(ITEM)
            PR(IPN)=-I1-PNEWR(ITEM)
            PL(IPN)=K1+PNEWL(ITEM)
            IPTID(IPN)=ITEM
          ELSE
C  PARTICLE DEFINITION WITH IPDL OR IPDA PACKAGE
cgzh varpt  
            PC(IPN)=PCORIG(IP)
            PR(IPN)=PRORIG(IP)
            PL(IPN)=PLORIG(IP)
cgzh store new origin
            PCORIG(IPN)=PCORIG(IP)
            PRORIG(IPN)=PRORIG(IP)
            PLORIG(IPN)=PLORIG(IP)
          END IF
cgzh randpt end if
        END IF
cgzh debug output
c      write(iouts,*) 'New pt ipn at pc,pr,pl:',ipn,
c     * PC(IPN),PR(IPN),PL(IPN)
cgzh debug skipping boundary face   
c      goto 681
cgzh boundary face
C FOR NEW PARTICLES AT SUBGRID BOUNDARY OR BFLX FACE PLACE PT ACCORDING
C TO MAGNITUDE OF VELOCITY OF SUBGRID FACE.  IF MORE THAN ONE
C CONTRIBUTING FACE, USE PROBABILITY FUNCTION TO PLACE PT.
C
C  IBDFC: SIGNIFIES NUMBER OF INFLOW BOUNDARY FACES
C  SUMVEL: SUM OF VELO'S IN FROM BOUNDARIES, TO BE USED BY
C          PROBABILITY FUNCTION 
C
C  BDFACE ARRAY SLOTS:
C    1: 0=NOT A FACE; 1=INFLOW FACE
C    2: NORMALIZED DISTANCE LIMITER (velocity*time/(cell-length-dimension))
C        This is the distance the pt would travel for the full time step,
C        normalized with the cell dimension (will be between 0 and 1, 
C        representing fraction of cell length pt would travel)
cgzh normalizing only works if CELDIS <= 1.0
C    3: MAGNITUDE OF VEL ON FACE 
C    4: LOWER BOUND FOR PROBABILITY FUNCTION
C    5: UPPER BOUND FOR PROBABILITY FUNCTION
C
C  JS1,IS1,KS1 is original cell
C  Velocity components are checked on the outside (subgrid) face
C    (remember that velo is stored on faces)
C
C  These six checks look for a positive inflow from the boundary face
C  BDFACE(I,1) is set to 1 (TRUE) if there is flow in 
C  BDFACE(I,2) is set to normalized distance limiter 
        IBDFC=0
C
C7-16 Check for inflow from boundary face.
C Check left column face
        IF(VC(JS1,IS1,KS1).GT.0.0) THEN
C Check subgrid boundary
          IF(JS1.EQ.1) THEN
            BDFACE(1,1)=1
C Check BFLX boundary
          ELSE IF (IBOUND(J1-1,I1,K1).EQ.0) THEN
            BDFACE(1,1)=1
          END IF
C If boundary flux found, set other parameters 
          IF(BDFACE(1,1).EQ.1) THEN
            IBDFC=IBDFC+1
C Convert VC to linear velocity
            VCLIN=VC(JS1,IS1,KS1)/
     *       (THCK(JS1,IS1,KS1)*POR(JS1,IS1,KS1)*rf(js1,is1,ks1)*RDEL)
            BDFACE(1,2)=(TIMV*VCLIN)/CDEL
            BDFACE(1,3)=VC(JS1,IS1,KS1)
          END IF
        END IF      
C Check right column face
        IF(VC(JS1+1,IS1,KS1).LT.0.0) THEN
          IF(JS1.EQ.NSCOL) THEN
            BDFACE(2,1)=1
          ELSE IF (IBOUND(J1+1,I1,K1).EQ.0) THEN
            BDFACE(2,1)=1
          END IF
          IF(BDFACE(2,1).EQ.1) THEN
            IBDFC=IBDFC+1
            VCLIN=VC(JS1+1,IS1,KS1)/
     *       (THCK(JS1,IS1,KS1)*POR(JS1,IS1,KS1)*rf(js1,is1,ks1)*RDEL)
            BDFACE(2,2)=ABS((TIMV*VCLIN)/CDEL)
            BDFACE(2,3)=ABS(VC(JS1+1,IS1,KS1))
          END IF
        END IF
C Check upper row face
        IF(VR(JS1,IS1,KS1).GT.0.0) THEN       
          IF(IS1.EQ.1) THEN       
            BDFACE(3,1)=1
          ELSE IF (IBOUND(J1,I1-1,K1).EQ.0) THEN
            BDFACE(3,1)=1
          END IF
          IF(BDFACE(3,1).EQ.1) THEN
            IBDFC=IBDFC+1
            VRLIN=VR(JS1,IS1,KS1)/
     *       (THCK(JS1,IS1,KS1)*POR(JS1,IS1,KS1)*rf(js1,is1,ks1)*CDEL)
            BDFACE(3,2)=(TIMV*VRLIN)/RDEL
            BDFACE(3,3)=VR(JS1,IS1,KS1)
          END IF
        END IF

C Check lower row face
        IF(VR(JS1,IS1+1,KS1).LT.0.0) THEN
          IF(IS1.EQ.NSROW) THEN
            BDFACE(4,1)=1
          ELSE IF (IBOUND(J1,I1+1,K1).EQ.0) THEN
            BDFACE(4,1)=1
          END IF
          IF(BDFACE(4,1).EQ.1) THEN
            IBDFC=IBDFC+1
            VRLIN=VR(JS1,IS1+1,KS1)/
     *       (THCK(JS1,IS1,KS1)*POR(JS1,IS1,KS1)*rf(js1,is1,ks1)*CDEL)
            BDFACE(4,2)=ABS((TIMV*VRLIN)/RDEL)
            BDFACE(4,3)=ABS(VR(JS1,IS1+1,KS1))
          END IF
        END IF
C Check upper layer face
        IF(VL(JS1,IS1,KS1).GT.0.0) THEN
          IF(KS1.EQ.1) THEN
            BDFACE(5,1)=1
          ELSE IF (IBOUND(J1,I1,K1-1).EQ.0) THEN
            BDFACE(5,1)=1
          END IF
          IF(BDFACE(5,1).EQ.1) THEN
              IBDFC=IBDFC+1
            VLLIN=VL(JS1,IS1,KS1)/
     *       (POR(JS1,IS1,KS1)*rf(js1,is1,ks1)*CDEL*RDEL)
            BDFACE(5,2)=(TIMV*VLLIN)/THCK(JS1,IS1,KS1)
            BDFACE(5,3)=VL(JS1,IS1,KS1)
          END IF
        END IF
C Check lower layer face
        IF(VL(JS1,IS1,KS1+1).LT.0.0) THEN
          IF(KS1.EQ.NSLAY) THEN
            BDFACE(6,1)=1
          ELSE IF (IBOUND(J1,I1,K1+1).EQ.0) THEN
            BDFACE(6,1)=1
          END IF
          IF(BDFACE(6,1).EQ.1) THEN
            IBDFC=IBDFC+1
            VLLIN=VL(JS1,IS1,KS1+1)/
     *       (POR(JS1,IS1,KS1)*rf(js1,is1,ks1)*CDEL*RDEL)
            BDFACE(6,2)=ABS((TIMV*VLLIN)/THCK(JS1,IS1,KS1))
            BDFACE(6,3)=ABS(VL(JS1,IS1,KS1+1))
          END IF
        END IF
C
C7-17 If there is a boundary face, restrict particle placement due to boundary face flux
C If no boundary faces, skip boundary face adjustment
        IF(IBDFC.NE.0) THEN
C
C Sum inflow boundary velos for use with probability function
C
          SUMVEL=0.0
          DO 650 I=1,6
 650        SUMVEL=SUMVEL+BDFACE(I,3)
C
C Set boundaries for probability function
C
          BOUNDLAST=0.0
          DO 660 I=1,6
            IF(BDFACE(I,1).GT.0) THEN
C   Set lower bound to last bound set (for first one, this is always 0.0)
              BDFACE(I,4)=BOUNDLAST
C   Set upper bound to fraction of flow this face accounts for 
C   BDFACE (I,3) is VEL
C   (If this is the only face, the fraction should be 1.0)
              BDFACE(I,5)=BOUNDLAST+(BDFACE(I,3)/SUMVEL)
C   Mark this as the last bound for next interval
              BOUNDLAST=BDFACE(I,5)
            END IF
 660      CONTINUE
C
C Restrict particle placement due to boundary face flux
cgzh randbd
C   Get random number; then check against intervals of probability function
C     defined above
          RANDVAL=RAN(ISEEDBD)
C   Check to see if flow in across this face
          DO 680 I=1,6
            IF(BDFACE(I,1).GT.0) THEN
C   Check probability function
              IF(BDFACE(I,4).LE.RANDVAL
     *          .AND.BDFACE(I,5).GT.RANDVAL) THEN
C   Recalculate position: remember particle positions are relative to
C     flow grid so use J instead of JS e.g.
C   Start at cell face (J1-0.5);
C     add adjustment: (PC-(J1-0.5)) is distance from face point is currently
C                     BDFACE(1,2) is adjustment factor  
C   Reverse all signs for forward faces
C     
cgzh debug output
c      write(iouts,*) 'before: pc,pr,pl',pc(ipn),pr(ipn),pl(ipn)
c      write(iouts,*) 'limiter=TIMV*VC(JS1,IS1,KS1))/CDEL',BDFACE(I,2)
c      write(iouts,*) 'TIMV,VCLIN,CDEL', TIMV,VCLIN,CDEL
c      write(iouts,*) 'before: pc',pc(ipn)
c           IF(I.EQ.1) PC(IPN)=J1-0.5+((PC(IPN)-(J1-0.5))*BDFACE(I,2))
c           IF(I.EQ.2) PC(IPN)=J1+0.5-((J1+0.5-PC(IPN))*BDFACE(I,2))
c           IF(I.EQ.3) PR(IPN)=(I1-0.5+((-PR(IPN)-(I1-0.5))*BDFACE(I,2)))
c     & *(-1)
c           IF(I.EQ.4) PR(IPN)=(I1+0.5-((I1+0.5+PR(IPN))*BDFACE(I,2)))
c     & *(-1)
c           IF(I.EQ.5) PL(IPN)=K1-0.5+((PL(IPN)-(K1-0.5))*BDFACE(I,2))
cgzh 7/12/2004 fixed typo here
c           IF(I.EQ.6) PL(IPN)=K1+0.5-((K1+0.5-PL(IPN))*BDFACE(I,2))
C REDUCE EXTENT TO AVOID ROUNDOFF ERROR
                IF(I.EQ.1) PC(IPN)=J1-0.4999+
     &            ((PC(IPN)-(J1-0.4999))*BDFACE(I,2))
                IF(I.EQ.2) PC(IPN)=J1+0.4999-
     &            ((J1+0.4999-PC(IPN))*BDFACE(I,2))
                IF(I.EQ.3) PR(IPN)=(I1-0.4999+
     &            ((-PR(IPN)-(I1-0.4999))*BDFACE(I,2)))*(-1)
                IF(I.EQ.4) PR(IPN)=(I1+0.4999-
     &            ((I1+0.4999+PR(IPN))*BDFACE(I,2)))*(-1)
                IF(I.EQ.5) PL(IPN)=K1-0.4999+
     &            ((PL(IPN)-(K1-0.4999))*BDFACE(I,2))
                IF(I.EQ.6) PL(IPN)=K1+0.4999-
     &            ((K1+0.4999-PL(IPN))*BDFACE(I,2))
C         
              END IF
            END IF
 680      CONTINUE
C ENDIF IBDFC.NE.0
        END IF
cgzh debug skipping boundary face
 681    CONTINUE
cgzh debug
        if(ipn.eq.361) then
          continue
        endif
C7-18 assign concentration to new particle and update other variables.		
        PCONC(IPN)=SRCAVC(JS1,IS1,KS1)
cgzh debug output
c      write(iouts,*) 'B: New pt ipn at pc,pr,pl:',ipn,
c     * PC(IPN),PR(IPN),PL(IPN),PCONC(IPN)
cgzh end randpt
C  FLAG THE WEIGHT OF THE PARTICLE
C  THIS WILL BE UPDATED IN SECOND LOOP OVER PARTICLES WHEN
C    THE NUMBER OF NEW PARTICLES IN THE SOURCE CELL IS KNOWN
        PTWT(IPN)=-99.0
        NPNEW(JS1,IS1,KS1)=NPNEW(JS1,IS1,KS1)+1
        NPCELL(JS1,IS1,KS1)=NPCELL(JS1,IS1,KS1)+1
C     ****************************************************************
C7-19           ---CHECK IF NEW LOCATION IS OUTSIDE OF SUBGRID---
  540   IF(NLOC.EQ.0) GO TO 590
C7-20  IF PARTICLE IS OUTSIDE SUBGRID OR IN INACTIVE CELL, 
C    ZERO ITS COMPONENTS AND STORE
C    ITS ID IN LIMBO SO AN "ACTIVE" PARTICLE CAN LATER USE THAT ID
  545   PC(IP)=0.0
        PR(IP)=0.0
        PL(IP)=0.0
        PCONC(IP)=0.0
        PTWT(IP)=0.0
C  
cgzh debug  only create limbo locations if not a new particle
        IF(IP.LE.NP) THEN
          DO 570 ID=1,NLIMBO
            IF(LIMBO(ID).GT.0) GO TO 570
            LIMBO(ID)=IP
            GO TO 590
cgzh debug output
c       write(iouts,*) 'debugging - PT PUT IN LIMBO, SUBGRID, ip=',ip
  570     CONTINUE
        END IF
C
  590 CONTINUE
cgzh debug output
cgzh debug zinnpt
c      if(imov.eq.68) then
c       write(iouts,*) 'After adv SUMWT row 4=',SUMWT(4,1,1)
c      end if
cgzh debug output
c      write(iouts,*) 'nptm after 1st loop=',nptm
cgzh debug pt loop
      CUMMASS=0.0
      CUMCELL=0.0
      CUMCELL2=0.0
      SUMVOL2=0.0
      SMCELL=0.0
C8     ---UPDATE CONC. OF PTS. IN ACTIVE CELLS---
      DO 1810 IP=1,NPTM
        
        IF(PC(IP).EQ.0.0) GO TO 1810
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
        IF(IBOUND(J,I,K).NE.0) THEN
cgzh debug output
c      if(js.eq.70.or.js.eq.71) then
c      write(iouts,*) 'js,ip,ptwt: ',js,ip,ptwt(ip)
c      end if
cgzh debug end
C       
CMOCWT  SUM WEIGHTS AND MASS 
          IF(PTWT(IP).GT.0.0) THEN
            CUMMASS=CUMMASS+PTWT(IP)*PCONC(IP)
            SMCELL(JS,IS,KS)=SMCELL(JS,IS,KS)+PTWT(IP)*PCONC(IP)
          END IF
cgzh debug
C RBW This includes new particles in source cells whose weights are less than zero
C RBW this probably should be fixed. It may not be needed.
            SUMVOL2=SUMVOL2+PTWT(IP)
cgzh debug
        ELSE
          WRITE(IOUTS,*) '**ERROR*** PARTICLE IN NO-FLOW CELL'
          WRITE(IOUTS,*) 'LOOP 1'
          WRITE(IOUTS,*) 'J,I,K,IP,PTWT',J,I,K,IP,PTWT(IP)
          STOP '**ERROR*** PARTICLE IN NO-FLOW CELL'
        END IF
 1810 CONTINUE
C9 calculate fractions to be added or deleted from particle volumes.
      SUMVOL=0.0
C   
      DO 1102 KS=1,NSLAY
      DO 1102 IS=1,NSROW
      DO 1102 JS=1,NSCOL
cgzh debug
        if(is.eq.7.and.ks.eq.1) then
          continue
        endif
cgzh debug output
c      if(js.eq.3.and.is.eq.1) then
c       write(iouts,*) 'a: js, summass',js,summass(js,is,ks)
c      end if
cgzh debug
c      if(imov.eq.8) then
c      write(iouts,*) 'sinkq components'
c      write(iouts,*) 'js,is,ks',js,is,ks
c      write(iouts,*) 'SNKFLO(JS,IS,KS)=',SNKFLO(JS,IS,KS)
c      write(iouts,*) 'BDYSNK(JS,IS,KS)=',BDYSNK(JS,IS,KS)
c      write(iouts,*) 'RESIDWT(JS,IS,KS)=',RESIDWT(JS,IS,KS)
c      end if
cgzh srcfix flag
        IF(ISRCFIX.GT.0) THEN
C
C  Calculate deficit between pts leaving strong source and source flux;
C    this difference will be corrected for
C
          IF(IGENPT(JS,IS,KS).GT.0.AND.SOURCEQ(JS,IS,KS).GT.0) THEN
corig        DEFICIT(JS,IS,KS)=WTOUT(JS,IS,KS)-SOURCEQ(JS,IS,KS)
corig     &                   -WTIN(JS,IS,KS)
cgzh srcfix2  include sinks and any weight that will be added back in from flow
cgzh srcfix2  to other strong sources (SS_WT)
c
C Begin 4 cases for new srcfix implementation
C This is case 1: too much left
            IF(WTOUT(JS,IS,KS).GT.TOTFLUXOT(JS,IS,KS)) THEN
c      write(iouts,*) 'Case 1 (too much left), cell js,is,ks=',js,is,ks
C  If too much left, redistribute the deficit on the pts in
C    the cell (this includes any new pts...) by removing weight 
C    proportionally from the pts that left this cell
              DEFICIT(JS,IS,KS)=WTOUT(JS,IS,KS)-TOTFLUXOT(JS,IS,KS)
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) 
c     * write(iouts,*) 'deficit111=',deficit(js,is,ks)
c
cgzh debug srcfix2   is this condition (wtout<deficit) ever met?  should not be with srcfix2
c
              IF(WTOUT(JS,IS,KS).LT.DEFICIT(JS,IS,KS)) THEN
cgzh debug output
c            write(iouts,*) 'Not enough to cover deficit in cell:',
c     &   js,is,ks
c            write(iouts,*) 'Deficit =',DEFICIT(js,is,ks)
c            DEFICIT(JS,IS,KS)=WTOUT(JS,IS,KS)
c            write(iouts,*) 'Deficit changed to wtout:',wtout(js,is,ks)
c            write(iouts,*) 'CELVOL=',CELVOL(js,is,ks)
c            write(iouts,*) 'SOURCEQ=',SOURCEQ(js,is,ks)
c            write(iouts,*) 'SINKQ=',SINKQ(js,is,ks)
c            write(iouts,*) 'SUMWT=',SUMWT(js,is,ks)
c            write(iouts,*) 'SS_WT=',SS_WT(js,is,ks)
c            write(iouts,*) 'vl+',vl(js,is,ks+1)*timv
c            write(iouts,*) 'vc-',vc(js,is,ks)*timv
c            write(iouts,*) 'vc+',vc(js+1,is,ks)*timv
c            write(iouts,*) 'vr-',vr(js,is,ks)*timv
c            write(iouts,*) 'vr+',vr(js,is+1,ks)*timv

c            write(*,*) 'Not enough to cover deficit in cell:',
c     &   js,is,ks
                stop 'unanticipated condition in srcfix: deficit3'
              END IF
C     
              LLST=1
              DO WHILE (LIST(ISRCID(JS,IS,KS),LLST).GT.0)
                IPLST=LIST(ISRCID(JS,IS,KS),LLST)
                WTOFF=(PTWT(IPLST)/WTOUT(JS,IS,KS))*DEFICIT(JS,IS,KS)
                PTWT(IPLST)=PTWT(IPLST)-WTOFF
cgzh debug output
c      if(imov.eq.68) 
c                 write(iouts,*) 'Shaved wt off leaving ip:',wtoff,iplst
                DEFMS(JS,IS,KS)=DEFMS(JS,IS,KS)+WTOFF*PCONC(IPLST)
cgzh debug output
c      if(imov.eq.95.or.imov.eq.96) 
c     &       write(iouts,*) 'Shaved ms off leaving ip:',
c     &  wtoff*pconc(iplst),iplst
                J=PC(IPLST)+0.5
                JJS=J-ISCOL1+1
                I=ABS(PR(IPLST))+0.5
                IIS=I-ISROW1+1
                K=PL(IPLST)+0.5
                KKS=K-ISLAY1+1
                TEMPWT(JJS,IIS,KKS)=TEMPWT(JJS,IIS,KKS)-WTOFF
                TEMPMS(JJS,IIS,KKS)=
     &            TEMPMS(JJS,IIS,KKS)-WTOFF*PCONC(IPLST)
c       
                LLST=LLST+1
C
C  IF BEYOND NPCELLMAX, STOP AND PRINT MESSAGE
                IF(LLST.GT.NPCELLMAX) THEN
                  WRITE(IOUTS,*) 
     *              '**ERROR** IN MOVEWT: NPCELLMAX in 1102 loop'
                  write(iouts,*) 'js,is,ks,npcellmax',
     *               js,is,ks,npcellmax
                  STOP 'ERROR IN MOVEWT: NPCELLMAX'
                END IF
c       
              END DO
C  Use SS_WT to distribute weight back onto pts in SRC cell 
cgzh debug (unfortunately this puts some weight on any new particles too)
              SS_WT(JS,IS,KS)=SS_WT(JS,IS,KS)+DEFICIT(JS,IS,KS)
cgzh debug output
c      if(imov.eq.68) 
cgzh debug output
c      write(iouts,*) 'Deficit added to SS_ms=',defms(js,is,ks)
              SS_MS(JS,IS,KS)=SS_MS(JS,IS,KS)+DEFMS(JS,IS,KS)
cgzh reset DEFMS (for case 3)
              DEFMS(JS,IS,KS)=0.0
cgzh debug
c      if(imov.eq.95.or.imov.eq.96) then
c      write(iouts,*) 'Deficit added to SS_wt=',deficit(js,is,ks)
c      write(iouts,*) 'Defms added to SS_ms=',DEFMS(js,is,ks)
c      write(iouts,*) 'js,is,ks=',js,is,ks
c      end if
cgzh debug
c move inclusion in sumwt to 890 or 1550
c          TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+DEFICIT(JS,IS,KS)
c          TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+DEFMS(JS,IS,KS)
cgzh debug IMPORTANT: if the "new" cell is a source, must update SRCWT and SRCMS too!
c
cgzh debug srcwt counted in 1590, 890, or 1550 only
c          SRCWT(JS,IS,KS)=SRCWT(JS,IS,KS)+DEFICIT(JS,IS,KS)
c          SRCMS(JS,IS,KS)=SRCMS(JS,IS,KS)+DEFMS(JS,IS,KS)
c      if(js.eq.2) then
c      write(iouts,*) 'srcwt update2'
c      end if
c
cgzh aq
c cycle through faces, look for inflow and outflow from non-ss cell and save that flux
c          IF(DEFICIT2(JS,IS,KS).GT.0) THEN
c moved totfluxin calc to ptwt file
            ENDIF
C This is case 2: not enough left
            IF(WTOUT(JS,IS,KS).LT.TOTFLUXOT(JS,IS,KS)) THEN
              DEFICIT(JS,IS,KS)=WTOUT(JS,IS,KS)-TOTFLUXOT(JS,IS,KS)
cgzh debug output
c      write(iouts,*) 'Case 2 (not enough left), cell js,is,ks=',js,is,ks
c      write(iouts,*) 'Deficit =',DEFICIT(js,is,ks)
c      write(iouts,*) 'wtout,totfluxout,celvol,sumwt,sourceq ',
c     & WTOUT(JS,IS,KS),TOTFLUXOT(JS,IS,KS),celvol(js,is,ks),
c     & sumwt(js,is,ks),sourceq(js,is,ks)*timv
            ENDIF
C This is case 3: too much entered
            IF(WTIN(JS,IS,KS).GT.TOTFLUXIN(JS,IS,KS)) THEN
c      write(iouts,*) 'Case 3 (too much entered),cell js,is,ks=',js,is,ks
              DEFICIT2(JS,IS,KS)=WTIN(JS,IS,KS)-TOTFLUXIN(JS,IS,KS)
c
cgzh debug srcfix2   is this condition (wtin<deficit2) ever met?  should not be with srcfix2
c
              IF(WTIN(JS,IS,KS).LT.DEFICIT2(JS,IS,KS)) THEN
cgzh debug output
c            write(iouts,*) 'Not enough to cover deficit in cell:',
c     &   js,is,ks
c            write(iouts,*) 'Deficit =',DEFICIT(js,is,ks)
                DEFICIT(JS,IS,KS)=WTOUT(JS,IS,KS)
c            write(iouts,*) 'Deficit changed to wtout:',wtout(js,is,ks)
c            write(iouts,*) 'CELVOL=',CELVOL(js,is,ks)
c            write(iouts,*) 'SOURCEQ=',SOURCEQ(js,is,ks)
c            write(iouts,*) 'SINKQ=',SINKQ(js,is,ks)
c            write(iouts,*) 'SUMWT=',SUMWT(js,is,ks)
c            write(iouts,*) 'SS_WT=',SS_WT(js,is,ks)
c            write(iouts,*) 'vl+',vl(js,is,ks+1)*timv
c            write(iouts,*) 'vc-',vc(js,is,ks)*timv
c            write(iouts,*) 'vc+',vc(js+1,is,ks)*timv
c            write(iouts,*) 'vr-',vr(js,is,ks)*timv
c            write(iouts,*) 'vr+',vr(js,is+1,ks)*timv

c            write(*,*) 'Not enough to cover deficit in cell:',
c     &   js,is,ks
                stop 'unanticipated condition in srcfix: deficit4'
              END IF
C     
              LLST=1
              DO WHILE (LIST2(ISRCID(JS,IS,KS),LLST).GT.0)
                IPLST=LIST2(ISRCID(JS,IS,KS),LLST)
                WTOFF=(PTWT(IPLST)/WTIN(JS,IS,KS))*DEFICIT2(JS,IS,KS)
                PTWT(IPLST)=PTWT(IPLST)-WTOFF
cgzh debug output
c      if(imov.eq.12.and.js.eq.11.and.is.eq.13) then 
c                write(iouts,*) 'Shaved wt off entering ip:',wtoff,iplst
c                write(iouts,*) 'pconc(iplst)',pconc(iplst)
c      end if
                DEFMS(JS,IS,KS)=DEFMS(JS,IS,KS)+WTOFF*PCONC(IPLST)
                J=PC(IPLST)+0.5
                JJS=J-ISCOL1+1
                I=ABS(PR(IPLST))+0.5
                IIS=I-ISROW1+1
                K=PL(IPLST)+0.5
                KKS=K-ISLAY1+1
                TEMPWT(JJS,IIS,KKS)=TEMPWT(JJS,IIS,KKS)-WTOFF
                TEMPMS(JJS,IIS,KKS)=
     *            TEMPMS(JJS,IIS,KKS)-WTOFF*PCONC(IPLST)
c        
                LLST=LLST+1
C
C  IF BEYOND NPCELLMAX, STOP AND PRINT MESSAGE
                IF(LLST.GT.NPCELLMAX) THEN
                  WRITE(IOUTS,*) 
     *              '**ERROR** IN MOVEWT: NPCELLMAX in 1102 loop'
                  write(iouts,*) 
     *              'js,is,ks,npcellmax',js,is,ks,npcellmax
                  STOP 'ERROR IN MOVEWT: NPCELLMAX'
                END IF
c         
              END DO
c        
C   Get random number; then check against intervals of probability function
              RANDVAL=RAN(ISEEDBD)
              ISRC=ISRCID(JS,IS,KS)
              DO 888 IFAC=1,6
C   If random value falls below calc'd prob, choose this face by stopping loop
                IF(RANDVAL.LT.SRCFAC(ISRC,IFAC,3)) GO TO 889
  888         CONTINUE
cgzh debug If it gets here, it never found and interval in srcfac, this is an error
              write(iouts,*) 
     &          ' ***ERROR*** SRCFAC NOT DEFINED FOR CELL (JS,IS,KS)',
     &          js,is,ks
              STOP 'ERROR: SRCFAC 2'
  889         CONTINUE
C       
              SELECT CASE (IFAC)
C IFAC=1  LEFT COLUMN FACE
                CASE (1) 
                  JJS=JS-1
                  IIS=IS
                  KKS=KS
C IFAC=2  RIGHT COLUMN FACE
                CASE (2) 
                  JJS=JS+1
                  IIS=IS
                  KKS=KS
C IFAC=3  UPPER ROW FACE
                CASE (3) 
                  JJS=JS
                  IIS=IS-1
                  KKS=KS
C IFAC=4  LOWER ROW FACE
                CASE (4) 
                  JJS=JS
                  IIS=IS+1
                  KKS=KS
C IFAC=5  UPPER LAYER FACE
                CASE (5) 
                  JJS=JS
                  IIS=IS
                  KKS=KS-1
C IFAC=6  LOWER LAYER FACE
                CASE (6) 
                  JJS=JS
                  IIS=IS
                  KKS=KS+1
              END SELECT
c
c

C  Use SS_WT to distribute weight back onto pts in randomly chosen contributing cell 

              SS_WT(JJS,IIS,KKS)=SS_WT(JJS,IIS,KKS)+DEFICIT2(JS,IS,KS)
              SS_MS(JJS,IIS,KKS)=SS_MS(JJS,IIS,KKS)+DEFMS(JS,IS,KS)
cgzh reset DEFMS
              DEFMS(JS,IS,KS)=0.0
cgzh debug
c      if(imov.eq.12.or.imov.eq.96) then
c      write(iouts,*) 'c3:Deficit2 added to SS_wt=',deficit2(js,is,ks)
c      write(iouts,*) 'js,is,ks=',js,is,ks
c      end if
            ENDIF
C This is case 4: not enough entered
            IF(WTIN(JS,IS,KS).LT.TOTFLUXIN(JS,IS,KS)) THEN
              DEFICIT2(JS,IS,KS)=WTIN(JS,IS,KS)-TOTFLUXIN(JS,IS,KS)
c      write(iouts,*) 'Case 4 (not enough enter),cell js,is,ks=',js,is,ks
c      write(iouts,*) 'deficit2=',DEFICIT2(JS,IS,KS)
            ENDIF
cgzh old definition of deficit
c        DEFICIT(JS,IS,KS)=CELVOL(JS,IS,KS)-
c     &                    (SOURCEQ(JS,IS,KS)+SINKQ(JS,IS,KS))
c     &                    -QSS_SINK(JS,IS,KS))
c     &                    -SUMWT(JS,IS,KS)-SS_WT(JS,IS,KS)
cgzh precision
            DEFCHECK=5E-8
cgzh debug end

cgzh IGENPT>0 end
          END IF
cgzh srcfix end
        END IF
C
C  CALCULATE FRACTIONAL WEIGHT TO BE ADDED TO EACH PARTICLE IN SOURCE CELLS
cgzh srcfix2 only do this with srcfix off, timing issue for adding pt in 1250 loop
        IF(ISRCFIX.EQ.0.AND.NPCELL(JS,IS,KS).GT.0
cgzh debug
cgzh 5/2/07 avoid divide by zero--but does this cause an imbalance if frac's are not defined?
cgzh 6/11/09 added non-strong source to this check
     *     .AND.SOURCEQ(JS,IS,KS).GT.0.AND.IGENPT(JS,IS,KS).EQ.0) THEN
           FRAC1(JS,IS,KS)=SOURCEQ(JS,IS,KS)/NPCELL(JS,IS,KS)
           FRAC2(JS,IS,KS)=SRCAVC(JS,IS,KS)*FRAC1(JS,IS,KS)
        END IF
C
C  SET PARTICLE REMOVAL ("REDISTRIBUTION") CRITERIA
cgzh if pt is less that critwt, its weight will be redistributed on other pts
cgzh this changed to CELVOL and not the current WT
cgzh debug
        if(imov.eq.256.and.js.eq.1.and.ks.eq.2.and.kkper.eq.1) then
c    if(js.eq.1.and.ks.eq.2.and.kkper.eq.2) then
          continue
        endif
c only do this once: celvol constant during time step
        CRITWT(JS,IS,KS)=ABS(CELVOL(JS,IS,KS)*REMCRIT)
c update critms every move as concentration may change
        CRITMS(JS,IS,KS)=ABS(SUMMASS(JS,IS,KS)*REMCRIT)
C  CALCULATE FRACTIONAL WEIGHT TO BE ADDED TO EACH PARTICLE IN SOURCE CELLS
        IF(NPCELL(JS,IS,KS).GT.0) THEN
          FRAC1(JS,IS,KS)=SOURCEQ(JS,IS,KS)/NPCELL(JS,IS,KS)
          FRAC2(JS,IS,KS)=SRCAVC(JS,IS,KS)*FRAC1(JS,IS,KS)
        ELSE
          CRITWT(JS,IS,KS)=0
          CRITMS(JS,IS,KS)=0
        END IF

c
        SUMVOL=SUMVOL+SUMWT(JS,IS,KS)
        CUMCELL=CUMCELL+SMCELL(JS,IS,KS)
        CUMCELL2=CUMCELL2+SUMMASS(JS,IS,KS)
 1102 CONTINUE
C10
cgzh update sumwt and summass again here
C  ADD WEIGHTS ALTERED BY 1102 LOOP TO SUM ARRAYS
      DO 295 KS=1,NSLAY
      DO 295 IS=1,NSROW
      DO 295 JS=1,NSCOL
C       
        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEMPWT(JS,IS,KS)
        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
C  REINITIALIZE 
        TEMPWT(JS,IS,KS)=0.0
        TEMPMS(JS,IS,KS)=0.0
 295  CONTINUE
C
cgzh debug output
c      WRITE(IOUTS,*) 'After 1st pt loop, SUMVOL=   ',SUMVOL
c      WRITE(IOUTS,*) 'After 1st pt loop, SUMVOL2=  ',SUMVOL2
          sumvol2=0.0
c      WRITE(IOUTS,*) 'After 1st pt loop, CUMMASS3=  ',CUMMASS3
c      WRITE(IOUTS,*) 'After 1st pt loop, CUMCELL=  ',CUMCELL
c      WRITE(IOUTS,*) 'After 1st pt loop,CUMCELL2=  ',CUMCELL2
      cumcell=0.
      cumcell2=0.
c      WRITE(IOUTS,*) 'After 1st pt loop, CUMMASS=  ',CUMMASS
cgzh debug output
c      WRITE(IOUTS,*) 'After 1st pt loop, SUMWT 17 14=  ',SUMWT(17,14,2)
C     ---END OF FIRST LOOP ON OLD PARTICLES---
cgzh debug SKIPPING DEBUG OUTPUT
      goto 999
      WRITE(IOUTS,330)
      DO 102 KS=1,NSLAY
      DO 102 IS=1,NSROW
  102   WRITE(IOUTS,331) (NPCELL(JS,IS,KS),JS=1,NSCOL)
  330 FORMAT (/,2X,'NPCELL AT END OF FIRST LOOP',/)
c      write(iouts,*) 'ptwt at end of first loop'
c      write(iouts,111) (PTWT(IP),IP=1,NPTM)
  111 format(20F8.3)
      write(iouts,*)
cgzh debug
 999  continue
C
C
C
C11  LOOP OVER PARTICLES TO FIX VOLUME AT STRONG SOURCES -- DEFICIT2
C
      IF(ISRCFIX.GT.0) THEN
        DO 490 IP=1,NPTM
          NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
          IF(NEWC.LE.0.0D0) GO TO 490
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
          J=INT(NEWC+0.5D0)
          JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
          NEWR=PR(IP)
          IF(NEWR.LT.0.0D0) THEN
             NEWR=-NEWR
          END IF
          I=INT(NEWR+0.5D0)
          IS=I-ISROW1+1
          NEWL=PL(IP)
          K=INT(NEWL+0.5D0)
          KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
          IF(IBOUND(J,I,K).EQ.0) THEN
            GO TO 490
          END IF
C  INITIALIZE TEMPORARY PT REDUCTION HOLDER
          PTWTOFF=0.D0
C
          IF(JS.LT.NSCOL) THEN
cgzh Case 4: not enough entered SS cell
            IF(DEFICIT2(JS+1,IS,KS).LT.0.0) THEN
cgzh aqu if this cell contributes to the right col face, remove weight based on
c    ratio to total flux into that cell
              IF(ICONTRIBIN(JS,IS,KS,1).GT.0) THEN
cgzh debug
                if(imov.eq.95) then
                  continue
                end if
c wtoff is percentage of face influx to total influx of deficit2
                WTOFF=DABS(((VC(JS+1,IS,KS)*TIMV)/
     &            TOTFLUXIN(JS+1,IS,KS))*DEFICIT2(JS+1,IS,KS))
c fracwt is weight to be removed from this particle, weighted for total sumwt of cell
c if not enough exists, do not continue
                IF(WTOFF.GT.SUMWT(JS,IS,KS)) THEN
                  WTOFF=SUMWT(JS,IS,KS)
cgzh debug output?
c          write(iouts,*) 'Case 4 deficit residual: not enough entered'
c          write(iouts,*) 'and contributing cell lacks enough weight'
c          write(iouts,*) 'js,is,ks,wtoff,sumwt',
c     *                    js,is,ks,wtoff,sumwt(js,is,ks)
c          stop 'case 4 residual'
                END IF
                IF(SUMWT(JS,IS,KS).GT.0.0) THEN
                  FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*WTOFF
                  FRACMS=FRACWT*PCONC(IP)
                ELSE
                  FRACWT=0.0
                  FRACMS=0.0
                END IF
c          PTWT(IP)=PTWT(IP)-FRACWT
cgzh instead of altering ptwt, keep track of this particle on all faces, then do adjustment
                PTWTOFF=PTWTOFF-FRACWT
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)-FRACWT
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)-FRACMS
c this weight and mass will be added back on via SS_WT in ss cell cell (here, js+1)
                SS_WT(JS+1,IS,KS)=SS_WT(JS+1,IS,KS)+FRACWT
                SS_MS(JS+1,IS,KS)=SS_MS(JS+1,IS,KS)+FRACMS
cgzh debug output
c      if((js+1).eq.61) then
c       write(iouts,*) 'C1:ss_wt update ',ss_wt(js+1,is,ks),js+1,is,ks
c      endif
              END IF
            END IF
          END IF
C
          IF(JS.GT.1) THEN
            IF(DEFICIT2(JS-1,IS,KS).LT.0.0) THEN
cgzh aqu if this cell contributes to the left col face
              IF(ICONTRIBIN(JS,IS,KS,2).GT.0) THEN
                WTOFF=DABS(((VC(JS,IS,KS)*TIMV)/
     &            TOTFLUXIN(JS-1,IS,KS))*DEFICIT2(JS-1,IS,KS))
                IF(WTOFF.GT.SUMWT(JS,IS,KS)) THEN
                  WTOFF=SUMWT(JS,IS,KS)
cgzh debug output?
c          write(iouts,*) 'Case 4 deficit residual: not enough entered'
c          write(iouts,*) 'and contributing cell lacks enough weight'
c          write(iouts,*) 'js,is,ks,wtoff,sumwt',
c     *                    js,is,ks,wtoff,sumwt(js,is,ks)
c          stop 'case 4 residual'
                END IF
                IF(SUMWT(JS,IS,KS).GT.0.0) THEN
                  FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*WTOFF
                  FRACMS=FRACWT*PCONC(IP)
                ELSE
                  FRACWT=0.0
                  FRACMS=0.0
                END IF
c         PTWT(IP)=PTWT(IP)-FRACWT
                PTWTOFF=PTWTOFF-FRACWT
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)-FRACWT
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)-FRACMS
                SS_WT(JS-1,IS,KS)=SS_WT(JS-1,IS,KS)+FRACWT
                SS_MS(JS-1,IS,KS)=SS_MS(JS-1,IS,KS)+FRACMS
              END IF
            END IF
          END IF
cgzh aqu if this cell contributes to the top row face
          IF(IS.LT.NSROW) THEN
            IF(DEFICIT2(JS,IS+1,KS).LT.0.0) THEN
              IF(ICONTRIBIN(JS,IS,KS,3).GT.0) THEN
                WTOFF=DABS(((VR(JS,IS+1,KS)*TIMV)/
     &            TOTFLUXIN(JS,IS+1,KS))*DEFICIT2(JS,IS+1,KS))
                IF(WTOFF.GT.SUMWT(JS,IS,KS)) THEN
                  WTOFF=SUMWT(JS,IS,KS)
cgzh debug output?
c          write(iouts,*) 'Case 4 deficit residual: not enough entered'
c          write(iouts,*) 'and contributing cell lacks enough weight'
c          write(iouts,*) 'js,is,ks,wtoff,sumwt',
c     *                    js,is,ks,wtoff,sumwt(js,is,ks)
c          stop 'case 4 residual'
                END IF
                IF(SUMWT(JS,IS,KS).GT.0.0) THEN
                  FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*WTOFF
                  FRACMS=FRACWT*PCONC(IP)
                ELSE
                  FRACWT=0.0
                  FRACMS=0.0
                END IF
c         PTWT(IP)=PTWT(IP)-FRACWT
                PTWTOFF=PTWTOFF-FRACWT
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)-FRACWT
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)-FRACMS
                 SS_WT(JS,IS+1,KS)=SS_WT(JS,IS+1,KS)+FRACWT
                 SS_MS(JS,IS+1,KS)=SS_MS(JS,IS+1,KS)+FRACMS
c       write(iouts,*) 'C3:ss_wt update ',ss_wt(js,is+1,ks),js,is+1,ks
              END IF
            END IF
          END IF
cgzh aqu if this cell contributes to the bottom row face
          IF(IS.GT.1) THEN
            IF(DEFICIT2(JS,IS-1,KS).LT.0.0) THEN
              IF(ICONTRIBIN(JS,IS,KS,4).GT.0) THEN
                WTOFF=DABS(((VR(JS,IS,KS)*TIMV)/TOTFLUXIN(JS,IS-1,KS))
     &            *DEFICIT2(JS,IS-1,KS))
                IF(WTOFF.GT.SUMWT(JS,IS,KS)) THEN
                  WTOFF=SUMWT(JS,IS,KS)
cgzh debug output?
c          write(iouts,*) 'Case 4 deficit residual: not enough entered'
c          write(iouts,*) 'and contributing cell lacks enough weight'
c          write(iouts,*) 'js,is,ks,wtoff,sumwt',
c    *                    js,is,ks,wtoff,sumwt(js,is,ks)
c          stop 'case 4 residual'
                END IF
                IF(SUMWT(JS,IS,KS).GT.0.0) THEN
                  FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*WTOFF
                  FRACMS=FRACWT*PCONC(IP)
                ELSE
                  FRACWT=0.0
                  FRACMS=0.0
                END IF
c         PTWT(IP)=PTWT(IP)-FRACWT
                PTWTOFF=PTWTOFF-FRACWT
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)-FRACWT
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)-FRACMS
                SS_WT(JS,IS-1,KS)=SS_WT(JS,IS-1,KS)+FRACWT
                SS_MS(JS,IS-1,KS)=SS_MS(JS,IS-1,KS)+FRACMS
c       write(iouts,*) 'C4:ss_wt update ',ss_wt(js,is-1,ks),js,is-1,ks
              END IF
            END IF
          END IF
cgzh aqu if this cell contributes to the upper layer face
          IF(KS.LT.NSLAY) THEN
            IF(DEFICIT2(JS,IS,KS+1).LT.0.0) THEN
              IF(ICONTRIBIN(JS,IS,KS,5).GT.0) THEN
                WTOFF=DABS(((VL(JS,IS,KS+1)*TIMV)/TOTFLUXIN(JS,IS,KS+1))
     &            *DEFICIT2(JS,IS,KS+1))
                IF(WTOFF.GT.SUMWT(JS,IS,KS)) THEN
                  WTOFF=SUMWT(JS,IS,KS)
cgzh debug output?
c          write(iouts,*) 'Case 4 deficit residual: not enough entered'
c          write(iouts,*) 'and contributing cell lacks enough weight'
c          write(iouts,*) 'js,is,ks,wtoff,sumwt',
c     *                    js,is,ks,wtoff,sumwt(js,is,ks)
c          stop 'case 4 residual'
                END IF
                IF(SUMWT(JS,IS,KS).GT.0.0) THEN
                  FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*WTOFF
                  FRACMS=FRACWT*PCONC(IP)
                ELSE
                  FRACWT=0.0
                  FRACMS=0.0
                END IF
c         PTWT(IP)=PTWT(IP)-FRACWT
                PTWTOFF=PTWTOFF-FRACWT
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)-FRACWT
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)-FRACMS
                 SS_WT(JS,IS,KS+1)=SS_WT(JS,IS,KS+1)+FRACWT
                 SS_MS(JS,IS,KS+1)=SS_MS(JS,IS,KS+1)+FRACMS
c        write(iouts,*) 'C5:ss_wt update ',ss_wt(js,is,ks+1),js,is,ks+1
             END IF
            END IF
          END IF
cgzh aqu if this cell contributes to the lower layer face
          IF(KS.GT.1) THEN
            IF(DEFICIT2(JS,IS,KS-1).LT.0.0) THEN
              IF(ICONTRIBIN(JS,IS,KS,6).GT.0) THEN
                WTOFF=DABS(((VL(JS,IS,KS)*TIMV)/TOTFLUXIN(JS,IS,KS-1))
     &            *DEFICIT2(JS,IS,KS-1))
                IF(WTOFF.GT.SUMWT(JS,IS,KS)) THEN
                  WTOFF=SUMWT(JS,IS,KS)
cgzh debug output?
c          write(iouts,*) 'Case 4 deficit residual: not enough entered'
c          write(iouts,*) 'and contributing cell lacks enough weight'
c          write(iouts,*) 'js,is,ks,wtoff,sumwt',
c     *                    js,is,ks,wtoff,sumwt(js,is,ks)
c          stop 'case 4 residual'
                END IF
                IF(SUMWT(JS,IS,KS).GT.0.0) THEN
                  FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*WTOFF
                  FRACMS=FRACWT*PCONC(IP)
                ELSE
                  FRACWT=0.0
                  FRACMS=0.0
                END IF
c                PTWT(IP)=PTWT(IP)-FRACWT
                PTWTOFF=PTWTOFF-FRACWT
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)-FRACWT
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)-FRACMS
                SS_WT(JS,IS,KS-1)=SS_WT(JS,IS,KS-1)+FRACWT
                SS_MS(JS,IS,KS-1)=SS_MS(JS,IS,KS-1)+FRACMS
c               write(iouts,*) 'C6:ss_wt update ',ss_wt(js,is,ks-1),js,is,ks-1
              END IF
            END IF
          END IF
C  REDUCE PARTICLE WEIGHT BASED ON REDUCTION TO ALL FACES (PTWTOFF is negative)
          PTWT(IP)=PTWT(IP)+PTWTOFF
C
  490   CONTINUE
      END IF
C
cgzh update sumwt and summass again here
C12  ADD WEIGHTS ALTERED BY 490 LOOP TO SUM ARRAYS
      DO 95 KS=1,NSLAY
      DO 95 IS=1,NSROW
      DO 95 JS=1,NSCOL
C
        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEMPWT(JS,IS,KS)
C  ADD MASS ALTERED BY SOURCES (OR "EXCESS" SINKS) TO SUM ARRAYS
        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
C  REINITIALIZE 
        TEMPWT(JS,IS,KS)=0.0
        TEMPMS(JS,IS,KS)=0.0
   95 CONTINUE
C
cgzh debug output
c        write(iouts,*) 'loop 1 pcorig(121)',pcorig(121)
c        write(iouts,*) 'loop 1 pc(121)',pc(121)
C     ****************************************************************
C
C13     LOOP 1 OVER CELLS TO FIX VOLUME AT STRONG SOURCES -- CREATE NEW PTS
cgzh srcfix2
cgzh debug
      IF(ISRCFIX.GT.0) THEN
        DO 1250 KS=1,NSLAY
        DO 1250 IS=1,NSROW
        DO 1250 JS=1,NSCOL
cgzh don't think the isrc check is needed here; is ss_wt ne 0, apply
c       ISRC=0
c        IF(INMNW.GT.0) THEN
c          IF (SRCMNW(JS,IS,KS).GT.0) ISRC=1
c        END IF
c        IF(SRCFLO(JS,IS,KS).GT.0.OR.BDYSRC(JS,IS,KS).GT.0) ISRC=1
cgzh srcfix2
c        IF(SS_SRC(JS,IS,KS).GT.0) ISRC=1
c        IF(ISRC.EQ.1) THEN
cgzh took out ss check: for new routine contributing cells could have ss_wt>0
c          IF(IGENPT(JS,IS,KS).GT.0.AND.SS_WT(JS,IS,KS).GT.0) THEN
            IF(SS_WT(JS,IS,KS).GT.0) THEN
C  ADD PARTICLE WITH WEIGHT EQUAL TO SS_WT
              IF (NPCELL(JS,IS,KS).EQ.0) THEN
                TEMPCP=SS_MS(JS,IS,KS)/SS_WT(JS,IS,KS)
cgzh debug
c            WRITE(IOUTS,*) 
c            WRITE(IOUTS,*) 'Calling ADDPTS:SS_WT>0, NPCELL=0'
c            WRITE(IOUTS,*) 'Creating 1 new particle at (js,is,ks):'
c            WRITE(IOUTS,*) js,is,ks
cgzh debug
cgzh varpt 
cgzh debug
                IF(INIPDL.EQ.0.AND.INIPDA.EQ.0) THEN
                  CALL ADDPTS(PC,PR,PL,PCONC,TEMPCP,IPTID,NPCELL,
     *             IBOUND,PNEWC,PNEWR,PNEWL,
     *             PTWT,SS_WT(JS,IS,KS),
     *             NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
C NEWPTS=2 means NPTPND=1
     *             2,NPMAX,
     *             JS,IS,KS,
     *             IOUTS,NPTM,WTFAC)
                ELSE
cgzh varpt 
cgzh set number of particles in each direction equal to 1 
                  NPTCOL=1
                  NPTROW=1
                  NPTLAY=1
                  CALL ADDPTSVAR(PC,PR,PL,PCONC,TEMPCP,
     *              NPCELL,PCORIG,PRORIG,PLORIG,
     *              PTWT,SS_WT(JS,IS,KS),
     *              NSCOL,NSROW,NSLAY,
     *              NPMAX,NPTCOL,NPTROW,NPTLAY,
     *              JS,IS,KS,
     *              IOUTS,NPTM,IDIM,WTFAC)
                END IF
c
C           TRACK AMOUNT ADDED BY THIS PROCESS
                SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+SS_WT(JS,IS,KS)
                SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)
     &                       +SS_WT(JS,IS,KS)*TEMPCP
cgzh srcfix2 debug
c so it is not double-counted, remove from srcwt, will be added back in in 1590 loop
c          SRCWT(JS,IS,KS)=SRCWT(JS,IS,KS)-SS_WT(JS,IS,KS)
c          SRCMS(JS,IS,KS)=SRCMS(JS,IS,KS)-SS_WT(JS,IS,KS)*TEMPCP
c      if(js.eq.2) then
c      write(iouts,*) 'srcwt update3'
c      end if
 
C  ZERO THESE OUT SO THEY DON"T TRIGGER IN NEXT PARTICLE LOOP (890)
                SS_WT(JS,IS,KS)=0.0
                SS_MS(JS,IS,KS)=0.0
              END IF
c isrc>0 uneeded now          END IF
            END IF
 1250   CONTINUE
      END IF
C     ****************************************************************
C
C14        ---SECOND LOOP OVER PARTICLES: ADJUST WEIGHTS DUE TO SOURCES---
C                   REMOVE LOW-WEIGHT PARTICLES AT SINKS
cgzh debug
      sumdum=0.0
cgzh debug output
      if(imov.eq.12) then
        sumvol3=0.0
        wtof=0.0
        wton=0.0
      endif
cgzh debug end
cgzh debug
c      if(imov.gt.184) then
c        write(55,*) 'imov=',imov
c        write(55,*) 'frac1,frac2',frac1(1,7,1),frac2(1,7,1)
c      end if
C  Loop over NPTM because it includes added particles
      DO 1590 IP=1,NPTM
cgzh debug
c      if(imov.eq.2.and.js.eq.1.and.is.eq.2) then
        if(imov.eq.2.and.(ip.eq.103)) then
          continue
        end if
        NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
        IF(NEWC.LE.0.0D0) GO TO 1590
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
        J=INT(NEWC+0.5D0)
        JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
        NEWR=PR(IP)
        IF(NEWR.LT.0.0D0) THEN
           NEWR=-NEWR
        END IF
        I=INT(NEWR+0.5D0)
        IS=I-ISROW1+1
        NEWL=PL(IP)
        K=INT(NEWL+0.5D0)
        KS=K-ISLAY1+1
        IF(JS.LT.1.OR.JS.GT.NSCOL.OR.IS.LT.1.OR.IS.GT.NSROW.OR.
     *    KS.LT.1.OR.KS.GT.NSLAY) THEN
           WRITE(IOUTS,*) ' IP,JS,IS,KS=',IP,JS,IS,KS
           WRITE(IOUTS,*) ' NEWC,NEWR,NEWL=',NEWC,NEWR,NEWL
           WRITE(IOUTS,*) ' SECOND LOOP PARTICLE ERROR MOVE, STOPPING'
           STOP ' PARTICLE ERROR IN MOVE'
        END IF
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
        IF(IBOUND(J,I,K).EQ.0) THEN
cgzh debug output?
            WRITE (*,*) ' PT. ',IP,' LOCATED IN INACTIVE CELL ',J,I,K
            GO TO 1590
        END IF
C      UPDATE WEIGHTS AND CONCENTRATIONS OF PTS. IN SOURCE CELLS
C     SOURCEQ IS TOTAL VOLUME INTO CELL

        IF (SOURCEQ(JS,IS,KS).GT.0.0) THEN
cgzh debug
c      if(i.eq.8.and.imov.ge.184) then
c        write(55,*) 'ip, ptwt, pconc',ip,ptwt(ip),pconc(ip)
c      end if
c      if(ip.eq.481) then
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) 
c       write(iouts,*) 'sourceq,npnew,ptwt,ip',
c     *  sourceq(js,is,ks),npnew(js,is,ks),ptwt(ip),ip
c      end if
cgzh debug output
c      if(ip.eq.13190.and.imov.eq.1460) then
c        write(iouts,*) 'stop the debugger'
c      end if
          if(js.eq.1.and.is.eq.1.and.ks.eq.1)
     *      sumdum=sumdum+ptwt(ip)
c     *  imov.lt.313) then
c        write(iouts,*) 'SOURCEQ,NPNEW at 17,14,2', SOURCEQ(JS,IS,KS),
c     *NPNEW(JS,IS,KS)
c        write(iouts,*) 'ip,ptwt,sumwt', ip,ptwt(ip),sumwt(js,is,ks)
c      end if
C
C            IF NEW PARTICLES IN THIS CELL, ADD ALL WEIGHT TO THEM (STRONG SOURCE)
          IF (NPNEW(JS,IS,KS).GT.0) THEN
C         PTWT IS SET TO -99.0 TO FLAG THIS PARTICULAR PARTICLE AS NEW
            IF (PTWT(IP).LT.-98.) THEN
C         NEW WEIGHT IS VOLUME IN DIVIDED AMONG NEW PARTICLES
              PTWT(IP)=SOURCEQ(JS,IS,KS)/REAL(NPNEW(JS,IS,KS))
C         TRACK VOLUME AND MASS ADDED TO CELL
              TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+PTWT(IP)
C         (PCONC for new particles is updated in first loop using SRCAVC)
              TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+PTWT(IP)*PCONC(IP)
cgzh debug output
c      if(js.eq.4.and.is.eq.1.and.ks.eq.1) then
c      write(iouts,*) 'debugging'
c      write(iouts,*) 'ADDING PARTICLE AT STRONG SOURCE'
c      write(iouts,*) 'Cell:',js,is,ks,' ip=',ip
c      write(iouts,*) 'Weight:',ptwt(ip)
c      write(iouts,*) 'NPNEW:',NPNEW(js,is,ks)
c      write(iouts,*) 'Pconc: ',pconc(ip)
c      sumdum=sumdum+ptwt(ip)
c      end if
cgzh debug output
            END IF      
          ELSE
C      IF NO NEW PARTICLES, ADD WEIGHT TO EXISTING PARTICLES (WEAK SOURCE)
C         NEW WEIGHT IS TOTAL WEIGHT DIVIDED AMONG NUMBER OF PARTICLES IN CELL
C
cgzh srcfix  only do the weak add if not fixing a deficit in a strong source
cgzh         (when we do that fix, we prefer to add a new pt, see next cell loop)
c
cgzh debug update: we made this inactive due to poor efficiency with this option
c            IF(DEFICIT(JS,IS,KS).GE.0.0) THEN
cgzh debug tried this originally (proportionally), now do it equally
Cold      NEW WEIGHT IS PROPORTIONAL TO THE RELATIVE WEIGHT OF THE PARTICLE
Cold        IN THE CELL, I.E. IF A PARTICLE ACCOUNTS FOR MORE OF THE EXISTING VOLUME
Cold        IT WILL RECEIVE MORE VOLUME FROM THE SOURCE
cgzh debug
C
C      IF NO PARTICLES AT ALL, WE WILL ACCOUNT FOR THIS SOURCE BY ADDING PARTICLES
C        IN NEXT LOOP OVER CELLS
cold (proportional--see above)                      FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*SOURCEQ(JS,IS,KS)
C      CALCULATE NEW PARTICLE CONCENTRATION: SUM EXISTING AND NEW MASS AND
C        DIVIDE BY SUM OF EXISTING AND NEW VOLUME
C
cgzh srcfix2 do this here if srcfix on, timing issue for adding pt in 1250 loop
cgzh added strong source flag to this, do not want to do it for weak sources
            IF(ISRCFIX.EQ.1.AND.IGENPT(JS,IS,KS).GT.0) THEN
              FRAC1(JS,IS,KS)=SOURCEQ(JS,IS,KS)/NPCELL(JS,IS,KS)
              FRAC2(JS,IS,KS)=SRCAVC(JS,IS,KS)*FRAC1(JS,IS,KS)
            END IF
            PCONC(IP)=(PCONC(IP)*PTWT(IP)+FRAC2(JS,IS,KS))/
     *                       (PTWT(IP)+FRAC1(JS,IS,KS))
C          
            PTWT(IP)=PTWT(IP)+FRAC1(JS,IS,KS)
cgzh debug output
c      if(imov.eq.1) then
c       if(js.eq.2.and.is.eq.2.and.ks.eq.2) then
c        write(iouts,*) '222 pconc,ip',pconc(ip),ip
c        write(iouts,*) 'frac1,frac2',FRAC1(JS,IS,KS),FRAC2(JS,IS,KS)
c      write(iouts,*) 'sourceq,srcavc',SOURCEQ(JS,IS,KS),SRCAVC(JS,IS,KS)
c      end if
c      end if
C      TRACK VOLUME AND MASS ADDED TO CELL
C      DO NOT PUT DIRECTLY INTO SUMWT AND SUMMASS BECAUSE SUMWT IS USED FOR THE 
C        WEAK SOURCE CALCULATION FOR ANY OTHER SUBSEQUENT PARTICLES IN THIS CELL
            TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+FRAC1(JS,IS,KS)
            TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+FRAC2(JS,IS,KS)
cgzh debug output
c      if(js.eq.4.and.is.eq.1.and.ks.eq.1) then
c     *  imov.lt.313) then
c      write(iouts,*) 'ADDING WEIGHT TO PTS AT SOURCE ("WEAKLY")'
c      write(iouts,*) 'Imov:',imov
c      write(iouts,*) 'Cell:',js,is,ks,' ip=',ip
c      write(iouts,*) 'Weight:',FRACWT
c      write(iouts,*) 'NPCELL:',NPCELL(js,is,ks)
c      write(iouts,*) 'sourceq:',SOURCEQ(JS,IS,KS)
c      write(iouts,*) 'srcavc:',srcavc(JS,IS,KS)
c      write(iouts,*) 'srcmnw:',srcmnw(JS,IS,KS)
c      end if
cgzh debug output
c           no deficit<0
c            END IF
c         npnew check
          END IF
C     sourceq gt 0
        END IF
cgzh debug
        SUMVOL3=SUMVOL3+PTWT(IP)
C  REMOVE PARTICLE IF ITS WEIGHT IS LESS THAN CRITERIA
C  CRITWT IS FRACTION (REMCRIT) OF SUM OF PARTICLE WEIGHTS IN CELL  
C  CRITMS IS FRACTION (REMCRIT) OF SUM OF PARTICLE MASSES IN CELL  
C  (REMCRIT) --> SET BY USER
cgzh debug only remove low weight pts if sink in cell
cgzh debug changed to: only if net sink
        SNKNET=SRCFLO(JS,IS,KS)+BDYSRC(JS,IS,KS)
     *         +SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS)
C       INMNW is for MNW1 
        IF((INMNW.GT.0).or.(INMNW2.GT.0))
     *    SNKNET=SNKNET+SRCMNW(JS,IS,KS)+SNKMNW(JS,IS,KS)
cgzh debug clark
        if(kkper.eq.60.and.js.eq.128.and.is.eq.64.and.ks.eq.12) then
          write(iouts,*) 'srcflo,snkflo',
     *      SRCFLO(JS,IS,KS),SNKFLO(JS,IS,KS)
          write(iouts,*) 'BDYSRC,BDYSNK',
     *      BDYSRC(JS,IS,KS),BDYSNK(JS,IS,KS)
          write(iouts,*) 'SRCMNW,SNKMNW',
     *      SRCMNW(JS,IS,KS),SNKMNW(JS,IS,KS)
          write(iouts,*) 'vl,vl+1',VL(JS,IS,KS),VL(JS,IS,KS+1)
          write(iouts,*) 'vc,vc+1',Vc(JS,IS,KS),Vc(JS+1,IS,KS)
          write(iouts,*) 'vr,vr+1',Vr(JS,IS,KS),Vr(JS,IS+1,KS)
          write(iouts,*) 'critwt,critms',
     *      critwt(JS,IS,KS),critms(JS,IS,KS)
          write(iouts,*) 'ptwt,ptmass',PTWT(ip),PTWT(IP)*PCONC(IP)
          write(iouts,*)'snknet,npcel,nptpnd',snknet,npcell(JS,IS,KS),
     &      nporig(js,is,ks)
          stop
        end if

        IF(SNKNET.LT.0.0) THEN
cdebug below line is old, replace with above 3 lines
c      IF(SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS).LT.0.0) THEN
C 
          PTMASS=PTWT(IP)*PCONC(IP)
          IF(PTWT(IP).LT.CRITWT(JS,IS,KS)
cgzh debug changed this to LE for no mass in cell case
c     *    .AND.PTMASS.LT.CRITMS(JS,IS,KS)
     *      .AND.PTMASS.LE.CRITMS(JS,IS,KS)
cgzh debug avoid removing unless NPTPND particles there
c     *    .AND.NPCELL(JS,IS,KS).GT.1) THEN
     *      .AND.NPCELL(JS,IS,KS).GT.(NPORIG(JS,IS,KS)+8)) THEN
C  ONLY REMOVE IF THERE IS AT LEAST ONE OTHER PARTICLE IN CELL
C  TRACK VOLUME AND MASS REMOVED, WILL BE ADDED BACK TO REMAINING 
C    PARTICLES IN NEXT PARTICLE LOOP; THIS CONSERVES MASS SO NOT IN MASS BAL
cgzh debug output
c      if(js.eq.2.and.is.eq.2.and.ks.eq.3) then
c         write(IOUTS,*) 'Particle',IP,' removed from',js,is,ks
c         write(IOUTS,*) 'SNKFLO=',SNKFLO(JS,IS,KS)
c         write(IOUTS,*) 'BDYSNK=',BDYSNK(JS,IS,KS)
c         write(IOUTS,*) 'PTWT=',PTWT(IP),' CRITWT=',CRITWT(JS,IS,KS)
c         write(IOUTS,*) 'PTMASS=',PTMASS,' CRITMS=',CRITMS(JS,IS,KS)
c      end if
              REMVWT(JS,IS,KS)=REMVWT(JS,IS,KS)+PTWT(IP)
              REMVMS(JS,IS,KS)=REMVMS(JS,IS,KS)+(PTWT(IP)*PCONC(IP))
              NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)-1
              PC(IP)=0.0
              PR(IP)=0.0
              PL(IP)=0.0
              PCONC(IP)=0.0
              PTWT(IP)=0.0
C  
cgzh debug  only create limbo locations if not a new particle
              IF(IP.LE.NP) THEN
                DO 1570 ID=1,NLIMBO
                  IF(LIMBO(ID).GT.0) GO TO 1570
                  LIMBO(ID)=IP
cgzh debug output
c      write(iouts,*) 'debugging - PT PUT IN LIMBO, CRITWT, ip=',ip
                  GO TO 1580
 1570           CONTINUE
              END IF
 1580       CONTINUE
          END IF
        END IF
cgzh srcnew
cgzh sum wt and mass for recalc of pconc in next loop
C   TRACK WT AND MASS AT SOURCE CELL FOR MIXING 
C   ONLY MIX STRONG SOURCES
        IF(SOURCEQ(JS,IS,KS).GT.0.0.AND.IGENPT(JS,IS,KS).EQ.1) THEN
C         INMNW is for MNW1 
          IF((INMNW.GT.0).or.(INMNW2.GT.0)) THEN
            SRC2=SRCFLO(JS,IS,KS)+SRCMNW(JS,IS,KS)
          ELSE
            SRC2=SRCFLO(JS,IS,KS)
          END IF
cgzh srcfix2 include ss_src in the check for mixing near a boundary?
          IF(ISRCFIX.GT.0) THEN
          SRC2=SRC2+SS_SRC(JS,IS,KS)
        END IF
C   ONLY MIX IF BOUNDARY SOURCE IS NOT STRONGEST SOURCE
C otherwise, concentration gradient coming in from bflx, e.g., would
c not be maintained.  check this carefully after srcfix2 revisions
cgzh this will then not mix for finite, angle2k etc
          IF(BDYSRC(JS,IS,KS).LT.SRC2) THEN
             SRCWT(JS,IS,KS)=SRCWT(JS,IS,KS)+PTWT(IP)
             SRCMS(JS,IS,KS)=SRCMS(JS,IS,KS)+PTWT(IP)*PCONC(IP)
c      if(js.eq.1.and.is.eq.1) then
c      write(iouts,*) 'srcwt update1',srcwt(1,1,1)
c      end if
          END IF
        END IF
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) then
c       write(iouts,*) '2nd loop ip,ptwt',ip,ptwt(ip)
c      end if
cgzh debug zinnpt
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c    &   and.is.eq.17) then
c       write(iouts,*) 'pt in 17: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.18) then
c       write(iouts,*) 'pt in 18: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
c         SUMVOL2=SUMVOL2+PTWT(IP)*PCONC(IP)
cgzh debug end
 1590 CONTINUE
cgzh debug
C
C
cgzh srcfix2
C15  ADD WEIGHTS ALTERED BY SOURCES TO SUM ARRAYS, MOVED THIS UP FOR 
C    SRCFIX2
      DO 1595 KS=1,NSLAY
      DO 1595 IS=1,NSROW
      DO 1595 JS=1,NSCOL
C
        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEMPWT(JS,IS,KS)
C  ADD MASS ALTERED BY SOURCES (OR "EXCESS" SINKS) TO SUM ARRAYS
        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
C  REINITIALIZE 
        TEMPWT(JS,IS,KS)=0.0
        TEMPMS(JS,IS,KS)=0.0
 1595 CONTINUE
C
C
C
C16  LOOP OVER PARTICLES SRCFIX: SHAVE WEIGHT OFF PTS IF NOT ENOUGH LEFT
C
      IF(ISRCFIX.GT.0) THEN
        DO 690 IP=1,NPTM
          NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
          IF(NEWC.LE.0.0D0) GO TO 690
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
          J=INT(NEWC+0.5D0)
          JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
          NEWR=PR(IP)
          IF(NEWR.LT.0.0D0) THEN
             NEWR=-NEWR
          END IF
          I=INT(NEWR+0.5D0)
          IS=I-ISROW1+1
          NEWL=PL(IP)
          K=INT(NEWL+0.5D0)
          KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
          IF(IBOUND(J,I,K).EQ.0) THEN
                GO TO 690
          END IF
C
          IF (SOURCEQ(JS,IS,KS).GT.0.0) THEN
cgzh srcfix
C If not enough left a strong source, shave that much off existing particles,
C in next loop over cells we will make a new particle in cell outside source with
C this mass
cgzh precision
            DEFCHECK=-5E-8
            IF(IGENPT(JS,IS,KS).GT.0.AND.
     &         DEFICIT(JS,IS,KS).LT.DEFCHECK) THEN
              IF(PTWT(IP).GT.0.AND.SUMWT(JS,IS,KS).GT.0.0) THEN
                WTOFF=(PTWT(IP)/SUMWT(JS,IS,KS))*DEFICIT(JS,IS,KS)
                PTWT(IP)=PTWT(IP)+WTOFF
c  DEFMS is always positive
                DEFMS(JS,IS,KS)=DEFMS(JS,IS,KS)-WTOFF*PCONC(IP)
cgzh debug
c      if(imov.eq.96.and.js.eq.3.and.ks.eq.2) then
c       write(iouts,*) 'mass removed from pt=',WTOFF*PCONC(IP)
c      end if
                TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+WTOFF
                TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+WTOFF*PCONC(IP)
cgzh remove this weight from SRCWT (sum of weights in src cells)
                SRCWT(JS,IS,KS)=SRCWT(JS,IS,KS)+WTOFF
                SRCMS(JS,IS,KS)=SRCMS(JS,IS,KS)+WTOFF*PCONC(IP)
              END IF
            END IF
          END IF
          SUMVOL2=SUMVOL2+PTWT(IP)*PCONC(IP)
  690   CONTINUE
cgzh srcfix end
cgzh srcfix2
C17  ADD WEIGHTS ALTERED BY SRCFIX TO SUM ARRAYS
        DO 1795 KS=1,NSLAY
        DO 1795 IS=1,NSROW
        DO 1795 JS=1,NSCOL
cgzh debug output
c      if(js.eq.32.or.js.eq.62.or.js.eq.92) then
c       write(iouts,*) 'c: js, summass',js,summass(js,is,ks)
c      end if
C
          SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEMPWT(JS,IS,KS)
C  ADD MASS ALTERED BY SOURCES (OR "EXCESS" SINKS) TO SUM ARRAYS
          SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
C  REINITIALIZE 
          TEMPWT(JS,IS,KS)=0.0
          TEMPMS(JS,IS,KS)=0.0
 1795   CONTINUE
c end if srcfix on
      END IF
C
C
cgzh debug output
c        write(iouts,*) 'loop 2 pcorig(121)',pcorig(121)
c        write(iouts,*) 'loop 2 pc(121)',pc(121)
c        write(iouts,*) 'loop 2 remvwt 223',remvwt(2,2,3)
cgzh debug
      SUMVOL=0.0
      sumtempwt=0.0
      sumsrc=0.0
      sumsrcq=0.0
      sumtemp=0.0
      sumrem=0.0
      sumremms=0.0
      DO 2102 KS=1,NSLAY
      DO 2102 IS=1,NSROW
      DO 2102 JS=1,NSCOL
        SUMVOL=SUMVOL+SUMWT(JS,IS,KS)
        sumtempwt=sumtempwt+tempwt(js,is,ks)
        sumsrc=sumsrc+(srcflo(js,is,ks)+bdysrc(js,is,ks))*timv
        if(js.eq.21.and.is.eq.20.and.ks.eq.13) then
          sumsrcq=sumsrcq+sourceq(js,is,ks)
          sumtemp=sumtemp+tempms(js,is,ks)
        end if
        sumrem=sumrem+remvwt(js,is,ks)
        sumremms=sumremms+remvms(js,is,ks)
         CUMCELL2=CUMCELL2+SUMMASS(JS,IS,KS)
 2102 CONTINUE
cgzh debug output
cgzh debug zinnpt
c      if(imov.eq.68) then
c       write(iouts,*) 'After 2nd ploop SUMWT row 4=',SUMWT(4,1,1)
c      end if
c      write(iouts,*) 'ptwt at end of second loop'
c      write(iouts,111) (PTWT(IP),IP=1,NPTM)
c      WRITE(IOUTS,*) 'After 2nd pt loop, SUMVOL=   ',SUMVOL
c      WRITE(IOUTS,*) 'After 2nd pt loop, SUMVOL2=   ',SUMVOL2
c      WRITE(IOUTS,*) 'After 2nd pt loop, SUMVOL3=   ',SUMVOL3
c      WRITE(IOUTS,*) 'After 2nd pt loop, CUMCELL2=   ',CUMCELL2
      cumcell2=0.
c      WRITE(IOUTS,*) 'After 2nd pt loop, TEMPWT=   ',SUMtempwt
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumsrc=   ',sumsrc
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumsrcq=   ',sumsrcq
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumtempms=',sumtemp
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumremwt= ',sumrem
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumremms= ',sumremms
cgzh debug output
cgzh debug output
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumwt 11 18 = ',sumwt(11,18,1)
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumwt 1 1 = ',sumwt(1,1,1)
c      WRITE(IOUTS,*) 'After 2nd pt loop, sumdum = ',sumdum
C
C
C     ****************************************************************
C
C18     LOOP 2 OVER CELLS TO FIX VOLUME AT STRONG SOURCES
cgzh srcfix
      IF(ISRCFIX.GT.0) THEN
        DO 1550 KS=1,NSLAY
        DO 1550 IS=1,NSROW
        DO 1550 JS=1,NSCOL
          ISRC=0
C         INMNW is for MNW1 
          IF((INMNW.GT.0).or.(INMNW2.GT.0)) THEN
            IF (SRCMNW(JS,IS,KS).GT.0) ISRC=1
          END IF
          IF(SRCFLO(JS,IS,KS).GT.0.OR.BDYSRC(JS,IS,KS).GT.0) ISRC=1
cgzh srcfix2
          IF(SS_SRC(JS,IS,KS).GT.0) ISRC=1
          IF(ISRC.EQ.1) THEN
C  deactivated this due to poor breakthrough with this method
c          IF (IGENPT(JS,IS,KS).GT.0.AND.DEFICIT(JS,IS,KS).LT.0.0) THEN
C  Use REMVWT to distribute weight onto pts in forward cell 
cgzh fixed for 1-d fint_age problem
c            REMVWT(JS+1,IS,KS)=REMVWT(JS+1,IS,KS)-DEFICIT(JS,IS,KS)
c            REMVMS(JS+1,IS,KS)=REMVMS(JS+1,IS,KS)+DEFMS(JS,IS,KS)
c            TEMPWT(JS+1,IS,KS)=TEMPWT(JS+1,IS,KS)-DEFICIT(JS,IS,KS)
c            TEMPMS(JS+1,IS,KS)=TEMPMS(JS+1,IS,KS)+DEFMS(JS,IS,KS)
c            REMVWT(js,IS,KS)=REMVWT(js,IS,KS)-DEFICIT(JS,IS,KS)
c            REMVMS(js,IS,KS)=REMVMS(js,IS,KS)+DEFMS(JS,IS,KS)
c            TEMPWT(js,IS,KS)=TEMPWT(js,IS,KS)-DEFICIT(JS,IS,KS)
c            TEMPMS(js,IS,KS)=TEMPMS(js,IS,KS)+DEFMS(JS,IS,KS)
cgzh debug output
c            write(iouts,*) 'Deficit added to remvwt=',deficit(js,is,ks)
c            write(iouts,*) 'Deficit added to remvms=',defms(js,is,ks)
c          END IF
c
c
cdebug skip this:  this skip is *OFF*
c      go to 789
c
c  If not enough left, create new particle in neighboring cell with deficit weight
c
cgzh precision
            DEFCHECK=-5E-8
            IF(IGENPT(JS,IS,KS).GT.0.AND.
     &         DEFICIT(JS,IS,KS).LT.DEFCHECK) THEN
            NPTM=NPTM+1
cgzh debug also need "restart" code here once restart is available
            IF(NPTM.EQ.NPMAX) THEN
              write(*,*)
              STOP 'NPMAX EXCEEDED IN MOVEWT'
            END IF
C
            IPN=NPTM
            J=JS+ISCOL1-1
            I=IS+ISROW1-1
            K=KS+ISLAY1-1
C SET WEIGHT AND CONCENTRATION OF NEW PT
            PTWT(IPN)=DABS(DEFICIT(JS,IS,KS))
cgzh debug
            wton=wton+ABS(DEFICIT(JS,IS,KS))
            PCONC(IPN)=DEFMS(JS,IS,KS)/PTWT(IPN)
C   Get random number; then check against intervals of probability function
            RANDVAL=RAN(ISEEDBD)
            ISRC=ISRCID(JS,IS,KS)
            DO 488 IFAC=1,6
C   If random value falls below calc'd prob, choose this face by stopping loop
              IF(RANDVAL.LT.SRCFAC(ISRC,IFAC,1)) GO TO 489
  488       CONTINUE
cgzh debug If it gets here, it never found and interval in srcfac, this is an error
            write(iouts,*) ' ***ERROR*** SRCFAC NOT DEFINED FOR CELL
     & (JS,IS,KS)',js,is,ks
            STOP 'ERROR: SRCFAC 1'
  489       CONTINUE
C
            SELECT CASE (IFAC)
C IFAC=1  LEFT COLUMN FACE
              CASE (1) 
                PC(IPN)=J-0.5+SRCFAC(ISRC,IFAC,2)
                PR(IPN)=I-0.4999+RAN(ISEEDBD)*0.99
                PL(IPN)=K-0.4999+RAN(ISEEDBD)*0.99
                NPCELL(JS-1,IS,KS)=NPCELL(JS-1,IS,KS)+1
                SUMWT(JS-1,IS,KS)=SUMWT(JS-1,IS,KS)+PTWT(IPN)
                SUMMASS(JS-1,IS,KS)=
     *            SUMMASS(JS-1,IS,KS)+PTWT(IPN)*PCONC(IPN)
C Particle added here, so make lumping a possibility here
                LUMP(JS-1,IS,KS)=1
C IFAC=2  RIGHT COLUMN FACE
              CASE (2) 
                PC(IPN)=J+0.5+SRCFAC(ISRC,IFAC,2)
                PR(IPN)=I-0.4999+RAN(ISEEDBD)*0.99
                PL(IPN)=K-0.4999+RAN(ISEEDBD)*0.99
                NPCELL(JS+1,IS,KS)=NPCELL(JS+1,IS,KS)+1
                SUMWT(JS+1,IS,KS)=SUMWT(JS+1,IS,KS)+PTWT(IPN)
                SUMMASS(JS+1,IS,KS)=
     *            SUMMASS(JS+1,IS,KS)+PTWT(IPN)*PCONC(IPN)
                LUMP(JS+1,IS,KS)=1
C IFAC=3  UPPER ROW FACE
              CASE (3) 
                PC(IPN)=J-0.4999+RAN(ISEEDBD)*0.99
                PR(IPN)=I-0.5+SRCFAC(ISRC,IFAC,2)
                PL(IPN)=K-0.4999+RAN(ISEEDBD)*0.99
                NPCELL(JS,IS-1,KS)=NPCELL(JS,IS-1,KS)+1
                SUMWT(JS,IS-1,KS)=SUMWT(JS,IS-1,KS)+PTWT(IPN)
                SUMMASS(JS,IS-1,KS)=
     *            SUMMASS(JS,IS-1,KS)+PTWT(IPN)*PCONC(IPN)
                LUMP(JS,IS-1,KS)=1
C IFAC=4  LOWER ROW FACE
              CASE (4) 
                PC(IPN)=J-0.4999+RAN(ISEEDBD)*0.99
                PR(IPN)=I+0.5+SRCFAC(ISRC,IFAC,2)
                PL(IPN)=K-0.4999+RAN(ISEEDBD)*0.99
                NPCELL(JS,IS+1,KS)=NPCELL(JS,IS+1,KS)+1
                SUMWT(JS,IS+1,KS)=SUMWT(JS,IS+1,KS)+PTWT(IPN)
                SUMMASS(JS,IS+1,KS)=
     *             SUMMASS(JS,IS+1,KS)+PTWT(IPN)*PCONC(IPN)
                LUMP(JS,IS+1,KS)=1
C IFAC=5  UPPER LAYER FACE
              CASE (5) 
                PC(IPN)=J-0.4999+RAN(ISEEDBD)*0.99
                PR(IPN)=I-0.4999+RAN(ISEEDBD)*0.99
                PL(IPN)=K-0.5+SRCFAC(ISRC,IFAC,2)
                NPCELL(JS,IS,KS-1)=NPCELL(JS,IS,KS-1)+1
                SUMWT(JS,IS,KS-1)=SUMWT(JS,IS,KS-1)+PTWT(IPN)
                SUMMASS(JS,IS,KS-1)=
     *            SUMMASS(JS,IS,KS-1)+PTWT(IPN)*PCONC(IPN)
                LUMP(JS,IS,KS-1)=1
C IFAC=6  LOWER LAYER FACE
              CASE (6) 
                PC(IPN)=J-0.4999+RAN(ISEEDBD)*0.99
                PR(IPN)=I-0.4999+RAN(ISEEDBD)*0.99
                PL(IPN)=K+0.5+SRCFAC(ISRC,IFAC,2)
                NPCELL(JS,IS,KS+1)=NPCELL(JS,IS,KS+1)+1
                SUMWT(JS,IS,KS+1)=SUMWT(JS,IS,KS+1)+PTWT(IPN)
                SUMMASS(JS,IS,KS+1)=
     *            SUMMASS(JS,IS,KS+1)+PTWT(IPN)*PCONC(IPN)
                LUMP(JS,IS,KS+1)=1
            END SELECT
cgzh store new origin 
cgzh debug this only works for IPDA and IDPL!!!!
c
c don't need these?  new pt is not originating in source 
c
c            PCORIG(IPN)=PC(IPN)
c            PRORIG(IPN)=PR(IPN)
c            PLORIG(IPN)=PL(IPN)
cgzh debug output
            if(imov.eq.1.and.js.eq.3.and.is.eq.1) then
            continue
            endif
c        write(iouts,*) 'Added new pt at forward cell, ip=',ipn
c        write(iouts,*) 'originating cell js,is,ks',js,is,ks
c        write(iouts,*) 'new pt at pc,pr,pl',pc(ipn),pr(ipn),pl(ipn)
c        write(iouts,*) 'ptwt,pconc=',PTWT(IPN),PCONC(IPN)
c        write(iouts,*) 'mass=',PTWT(IPN)*PCONC(IPN)

c      end if
cgzh debug output end
c  Also for this case, if there were no new particles, create one in src cell
c
c skip this add! (npnew never lt 0): this skip is *ON*
c
c            IF(NPNEW(JS,IS,KS).LT.0) THEN
c            IF(NPNEW(JS,IS,KS).EQ.0) THEN
cgzh debug check
c      diffwt=sourceq(js,is,ks)+deficit(js,is,ks)
c      if(diffwt.gt.1e-4) then
c        write(iouts,*) 'No new pts, sourceq.ne.deficit: diffwt=',diffwt
c        write(iouts,*) 'js,is,ks,sourceq,deficit',js,is,ks,
c    *    sourceq(js,is,ks),deficit(js,is,ks)
c      end if        
cgzh debug check end
c              NPTM=NPTM+1
c              IPN2=NPTM
c  left-center of cell
c              PC(IPN2)=JS-0.25
c              PC(IPN2)=JS-0.5+(0.5*TIMV*VCLIN)
c              PR(IPN2)=-IS
c              PL(IPN2)=KS
c              PTWT(IPN2)=SOURCEQ(JS,IS,KS)
c              PCONC(IPN2)=SRCAVC(JS,IS,KS)
c              NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)+1
c              TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+PTWT(IPN2)
c              TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+PTWT(IPN2)*PCONC(IPN2)
c              SRCWT(JS,IS,KS)=SRCWT(JS,IS,KS)+PTWT(IPN2)
c              SRCMS(JS,IS,KS)=SRCMS(JS,IS,KS)+PTWT(IPN2)*PCONC(IPN2)
C  PARTICLE DEFINITION WITH IPDL OR IPDA PACKAGE
cgzh store new origin 
cgzh debug this only works for IPDA and IDPL!!!!
c              PCORIG(IPN2)=PC(IPN2)
c              PRORIG(IPN2)=PR(IPN2)
c              PLORIG(IPN2)=PL(IPN2)
cgzh debug output
c        write(iouts,*) 'Added new pt at src cell, ip=', ipn2
c        write(iouts,*) 'js,is,ks',js,is,ks
c        write(iouts,*) 'pc,pr,pl',pc(ipn2),pr(ipn2),pl(ipn2)
c      write(iouts,*) 'ptwt,pconc=',PTWT(IPN2),PCONC(IPN2)
c      write(iouts,*) 'sourceq,srcavc=', 
c     * SOURCEQ(JS,IS,KS),SRCAVC(JS,IS,KS)
cgzh debug output end
c
          END IF
          END IF
  789   continue
 1550   CONTINUE
cgzh srcfix end
      END IF
C
C19  LOOP OVER PARTICLES TO FIX VOLUME AT STRONG SOURCES -- DISTRIBUTED
C
      IF(ISRCFIX.GT.0) THEN
        DO 890 IP=1,NPTM
          NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
          IF(NEWC.LE.0.0D0) GO TO 890
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
          J=INT(NEWC+0.5D0)
          JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
          NEWR=PR(IP)
          IF(NEWR.LT.0.0D0) THEN
             NEWR=-NEWR
          END IF
          I=INT(NEWR+0.5D0)
          IS=I-ISROW1+1
          NEWL=PL(IP)
          K=INT(NEWL+0.5D0)
          KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
          IF(IBOUND(J,I,K).EQ.0) THEN
                GO TO 890
          END IF
cgzh srcfix2
          IF(SS_WT(JS,IS,KS).GT.0) THEN
cgzh srcfix2
C  CALCULATE NEW PARTICLE CONCENTRATION: SUM EXISTING AND NEW MASS AND
C    DIVIDE BY SUM OF EXISTING AND NEW VOLUME
C
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1.and.imov.eq.1) then
c      write(iouts,*) '111 Distributing TO PTS VIA SS_WT at',js,is,ks
c      write(iouts,*) 'ss_wt,ss_ms',SS_WT(JS,IS,KS),SS_MS(JS,IS,KS)
c      write(iouts,*) 'ip,ptwt,pconc before',ip,ptwt(ip),pconc(ip)
c      end if
            TEWT=SS_WT(JS,IS,KS)/NPCELL(JS,IS,KS)
            TEMS=SS_MS(JS,IS,KS)/NPCELL(JS,IS,KS)
            PCONC(IP)=(PCONC(IP)*PTWT(IP)+TEMS)
     *                 /(PTWT(IP)+TEWT)
            PTWT(IP)=PTWT(IP)+SS_WT(JS,IS,KS)/NPCELL(JS,IS,KS)
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1.and.imov.eq.1) then
c      write(iouts,*) 'ip,ptwt,pconc after ',ip,ptwt(ip),pconc(ip)
c      end if
C  ADD WEIGHT AND MASS TO SUM ARRAYS
            SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEWT
            SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMS
C  COUNT THIS MASS AND WEIGHT FOR MIXING AT SOURCE
C   TRACK WT AND MASS AT SOURCE CELL FOR MIXING 
C   ONLY MIX STRONG SOURCES
C           INMNW is for MNW1 
            IF((INMNW.GT.0).or.(INMNW2.GT.0)) THEN
             SRC2=SRCFLO(JS,IS,KS)+SRCMNW(JS,IS,KS)
            ELSE
             SRC2=SRCFLO(JS,IS,KS)
            END IF
cgzh srcfix2 include ss_src in the check for mixing near a boundary?
            IF(ISRCFIX.GT.0) THEN
              SRC2=SRC2+SS_SRC(JS,IS,KS)
            END IF
            IF(BDYSRC(JS,IS,KS).LT.SRC2) THEN
              SRCWT(JS,IS,KS)=SRCWT(JS,IS,KS)+TEWT
              SRCMS(JS,IS,KS)=SRCMS(JS,IS,KS)+TEMS
c      if(js.eq.1.and.is.eq.1) then
c      write(iouts,*) 'srcwt update1',srcwt(1,1,1)
c      end if
            END IF 
          END IF 
C
  890   CONTINUE
C END LOOP OVER PARTICLES TO FIX VOLUME AT STRONG SOURCES
      END IF
C
C20  LOOP OVER CELLS TO:
C  CHECK FOR SOURCES WITH NO PARTICLES AND CREATE THEM THERE
C  UPDATE SUM OF WEIGHTS AND MASS IN CELLS
C  CALCULATE VOLUME TO BE REMOVED AT SINKS, CORRECTING FOR VOLUME THAT 
C    IS REMOVED BY PARTICLES LEAVING THE SUBGRID
C  CHECK FOR SINKS WITH NO PARTICLES AND CARRY OVER A RESIDUAL
cgzh drycell
C  PUT DRY WEIGHT AND MASS INTO CELLWISE ARRAY FOR DISTRIBUTION TO PTS
cgzh debug
      DO 1600 KS=1,NSLAY
      DO 1600 IS=1,NSROW
      DO 1600 JS=1,NSCOL
cgzh debug output
c      if(js.eq.32.or.js.eq.62.or.js.eq.92) then
c       write(iouts,*) 'e: js, summass',js,summass(js,is,ks)
c      end if
C
cgzh drycell
        IF(IMOV.EQ.1.AND.NDRYCELL.GT.0) THEN
          J=JS+ISCOL1-1
          I=IS+ISROW1-1
          K=KS+ISLAY1-1
          IF(DRYWT(JS,IS,KS).GT.0.0) THEN
            IDONE=0
C CHECK FOR CELL BELOW DRY CELL
            IF(KS.LT.NSLAY) THEN
              IF(IBOUND(J,I,K+1).NE.0) THEN
                TEMPDRYWT(JS,IS,KS+1)=TEMPDRYWT(JS,IS,KS+1)+
     &            DRYWT(JS,IS,KS)
                TEMPDRYMS(JS,IS,KS+1)=TEMPDRYMS(JS,IS,KS+1)+
     &            DRYMS(JS,IS,KS)
                IDONE=1
              END IF
            END IF
C IF FIRST CONDITION FOUND, SKIP OTHER CHECKS
            IF(IDONE.EQ.0) THEN
              IDONE=1
C FIND LATERAL NEIGHBOR WITH LOWEST HEAD, PUT WT AND MASS THERE
              HMIN=1E20
              ICASE=0
              IF(JS.GT.1) THEN
                IF(IBOUND(J-1,I,K).NE.0) THEN
                  IF(HNEW(J-1,I,K).LT.HMIN) THEN
                    HMIN=HNEW(J-1,I,K)
                    ICASE=1
                  ENDIF
                END IF
              END IF
              IF(IS.GT.1) THEN
                IF(IBOUND(J,I-1,K).NE.0) THEN
                  IF(HNEW(J,I-1,K).LT.HMIN) THEN
                    HMIN=HNEW(J,I-1,K)
                    ICASE=2
                  ENDIF
                END IF
              END IF
              IF(JS.LT.NSCOL) THEN
                IF(IBOUND(J+1,I,K).NE.0) THEN
                  IF(HNEW(J+1,I,K).LT.HMIN) THEN
                    HMIN=HNEW(J+1,I,K)
                    ICASE=3
                  ENDIF
                END IF
              END IF
              IF(IS.LT.NSROW) THEN
                IF(IBOUND(J,I+1,K).NE.0) THEN
                  IF(HNEW(J,I+1,K).LT.HMIN) THEN
                    HMIN=HNEW(J,I+1,K)
                    ICASE=4
                  ENDIF
                END IF
              END IF
              IF(JS.GT.1.AND.IS.GT.1) THEN
                IF(IBOUND(J-1,I-1,K).NE.0) THEN
                  IF(HNEW(J-1,I-1,K).LT.HMIN) THEN
                    HMIN=HNEW(J-1,I-1,K)
                    ICASE=5
                  ENDIF
                END IF
              END IF
              IF(JS.LT.NSCOL.AND.IS.GT.1) THEN
                IF(IBOUND(J+1,I-1,K).NE.0) THEN
                  IF(HNEW(J+1,I-1,K).LT.HMIN) THEN
                    HMIN=HNEW(J+1,I-1,K)
                    ICASE=6
                  ENDIF
                END IF
              END IF
              IF(JS.LT.NSCOL.AND.IS.LT.NSROW) THEN
                IF(IBOUND(J+1,I+1,K).NE.0) THEN
                  IF(HNEW(J+1,I+1,K).LT.HMIN) THEN
                    HMIN=HNEW(J+1,I+1,K)
                    ICASE=7
                  ENDIF
                END IF
              END IF
              IF(JS.GT.1.AND.IS.LT.NSROW) THEN
                IF(IBOUND(J-1,I+1,K).NE.0) THEN
                  IF(HNEW(J-1,I+1,K).LT.HMIN) THEN
                    HMIN=HNEW(J-1,I+1,K)
                    ICASE=8
                  ENDIF
                END IF
              END IF
C
              SELECT CASE(ICASE)
                CASE(0)
                  IDONE=0
                CASE(1)
                  TEMPDRYWT(JS-1,IS,KS)=
     &              TEMPDRYWT(JS-1,IS,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS-1,IS,KS)=
     &              TEMPDRYMS(JS-1,IS,KS)+DRYMS(JS,IS,KS)
                CASE(2)
                  TEMPDRYWT(JS,IS-1,KS)=
     &              TEMPDRYWT(JS,IS-1,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS,IS-1,KS)=
     &              TEMPDRYMS(JS,IS-1,KS)+DRYMS(JS,IS,KS)
                CASE(3)
                  TEMPDRYWT(JS+1,IS,KS)=
     &              TEMPDRYWT(JS+1,IS,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS+1,IS,KS)=
     &              TEMPDRYMS(JS+1,IS,KS)+DRYMS(JS,IS,KS)
                CASE(4)
                  TEMPDRYWT(JS,IS+1,KS)=
     &              TEMPDRYWT(JS,IS+1,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS,IS+1,KS)= 
     &              TEMPDRYMS(JS,IS+1,KS)+DRYMS(JS,IS,KS)
                CASE(5)
                  TEMPDRYWT(JS-1,IS-1,KS)=
     &              TEMPDRYWT(JS-1,IS-1,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS-1,IS-1,KS)=
     &              TEMPDRYMS(JS-1,IS-1,KS)+DRYMS(JS,IS,KS)
                CASE(6)
                  TEMPDRYWT(JS+1,IS-1,KS)=
     &              TEMPDRYWT(JS+1,IS-1,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS+1,IS-1,KS)=
     &              TEMPDRYMS(JS+1,IS-1,KS)+DRYMS(JS,IS,KS)
                CASE(7)
                  TEMPDRYWT(JS+1,IS+1,KS)=
     &              TEMPDRYWT(JS+1,IS+1,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS+1,IS+1,KS)=
     &              TEMPDRYMS(JS+1,IS+1,KS)+DRYMS(JS,IS,KS)
                CASE(8)
                  TEMPDRYWT(JS-1,IS+1,KS)=
     &              TEMPDRYWT(JS-1,IS+1,KS)+DRYWT(JS,IS,KS)
                  TEMPDRYMS(JS-1,IS+1,KS)=
     &              TEMPDRYMS(JS-1,IS+1,KS)+DRYMS(JS,IS,KS)
              END SELECT
            END IF                      
            IF(IDONE.EQ.0) THEN
                WRITE(IOUTS,*) '***WARNING***'
                WRITE(IOUTS,*) 'NO NEIGHTBOR CELL FOUND FOR PARTICLES
     & TRAPPED IN DRY CELL AT JS,IS,KS=', JS,IS,KS
            ENDIF
          END IF
        END IF
cgzh drycell end
        ISRC=0
C       INMNW is for MNW1 
        IF((INMNW.GT.0).or.(INMNW2.GT.0)) THEN
          IF (SRCMNW(JS,IS,KS).GT.0) ISRC=1
        END IF
        IF(SRCFLO(JS,IS,KS).GT.0.OR.BDYSRC(JS,IS,KS).GT.0) ISRC=1
        IF(ISRC.EQ.1) THEN
C
          IF(NPCELL(JS,IS,KS).EQ.0) THEN
C         GENERATE INITIAL DISTRIBUTION OF PARTICLES IN THIS CELL
C           AND EQUALLY DISTRIBUTE SOURCE FLUID (SOURCEQ) TO WEIGHTS
C         SET PCONC WITH SRCAVC, NPTM UPDATED
C
cgzh debug
c      WRITE(IOUTS,*) 
c      WRITE(IOUTS,*)'Calling ADDPTS: NPCELL=0 at weak source.'
c      WRITE(IOUTS,*) 'Creating',NEWPTS-1,' new particles at (js,is,ks):'
c      WRITE(IOUTS,*) js,is,ks
c      WRITE(IOUTS,*) 'SOURCEQ=', SOURCEQ(JS,IS,KS)
c      WRITE(IOUTS,*) 'SRCAVC= ', srcavc(js,is,ks)
cgzh varpt 
            IF(INIPDL.EQ.0.AND.INIPDA.EQ.0) THEN
              CALL ADDPTS(PC,PR,PL,PCONC,SRCAVC(JS,IS,KS),IPTID,
     *          NPCELL,IBOUND,PNEWC,PNEWR,PNEWL,
     *          PTWT,SOURCEQ(JS,IS,KS),
     *          NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *          NEWPTS,NPMAX,
     *          JS,IS,KS,
     *          IOUTS,NPTM,WTFAC)
            ELSE
cgzh varpt 
cgzh get particle distribution for this cell
              NPTCOL=NPTCOLA(JS,IS,KS)
              NPTROW=NPTROWA(JS,IS,KS)
              NPTLAY=NPTLAYA(JS,IS,KS)
              CALL ADDPTSVAR(PC,PR,PL,PCONC,SRCAVC(JS,IS,KS),
     *          NPCELL,PCORIG,PRORIG,PLORIG,
     *          PTWT,SOURCEQ(JS,IS,KS),
     *          NSCOL,NSROW,NSLAY,
     *          NPMAX,NPTCOL,NPTROW,NPTLAY,
     *          JS,IS,KS,
     *          IOUTS,NPTM,IDIM,WTFAC)
            END IF
C          TRACK AMOUNT ADDED BY THIS PROCESS
cgzh stet
cgzh srcfix2 changed these from tempwt due to moving tempwt addition above
c            TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+SOURCEQ(JS,IS,KS)
c            TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+
c     *                         SOURCEQ(JS,IS,KS)*SRCAVC(JS,IS,KS)
            SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+SOURCEQ(JS,IS,KS)
            SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+
     *                         SOURCEQ(JS,IS,KS)*SRCAVC(JS,IS,KS)

          END IF      
        END IF
cgzh
C  ADD WEIGHTS ALTERED BY SOURCES TO SUM ARRAYS, NEED THIS HERE FOR 
C    USE WITH SINKQ CHECK
c        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEMPWT(JS,IS,KS)
C
C  DO CELLWISE CALCULATIONS FOR SINKS
C
cgzh srcfix2
C  SINKQ CALCULATED ABOVE
C
C  SINKQ IS VOLUME LEAVING CELL AS CALCULATED BY FD METHOD
C    PLUS RESIDUAL VOLUME (RESIDWT) DUE TO A SINK AT A CELL WITH
C    NO PARTICLES (CARRIED OVER FROM ANY PREVIOUS MOVES)
C
C    IF SINKQ = 0, DO NOT ADJUST PARTICLE WEIGHTS
C    IF SINKQ < 0, THEN VOLUME WILL BE REMOVED FROM PARTICLES
C    IF SINKQ > 0, THEN VOLUME WILL BE ADDED TO PARTICLES
C    IF NO PARTICLES, ADD TO RESIDUAL AND ATTEMPT TO APPLY IT NEXT MOVE
C
c        SINKQ(JS,IS,KS)=TIMV*(SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS))
c     *                 +RESIDWT(JS,IS,KS)
cgzh srcfix2 this has been moved to top
cgzh debug output
c      if(imov.eq.8) then
c      write(iouts,*) 'sinkq components'
c      write(iouts,*) 'js,is,ks',js,is,ks
c      write(iouts,*) 'SNKFLO(JS,IS,KS)=',SNKFLO(JS,IS,KS)
c      write(iouts,*) 'BDYSNK(JS,IS,KS)=',BDYSNK(JS,IS,KS)
c      write(iouts,*) 'RESIDWT(JS,IS,KS)=',RESIDWT(JS,IS,KS)
c      end if
C  RESIDUAL FROM LAST MOVE IS NOW ACCOUNTED FOR IN SINKQ, SO REINITIALIZE IT
        RESIDWT(JS,IS,KS)=0.0
C
C  IF PARTICLES LEFT THE SUBGRID AT THIS CELL DURING THIS MOVE, 
C    SUBTRACT THAT WEIGHT (SUMSGWT) FROM THE TOTAL VOLUME LEAVING SINK
        IF(SUMSGWT(JS,IS,KS).LT.0) THEN
          SINKQ(JS,IS,KS)=SINKQ(JS,IS,KS)-SUMSGWT(JS,IS,KS)
        END IF
C
C  IF NET AFTER APPLYING SUMSGWT IS TO REMOVE WATER...
        IF (SINKQ(JS,IS,KS).LT.0.0) THEN
C  IF NO PARTICLES IN CELL, CAN'T APPLY SINK VOLUME, SO STORE IN 
C    RESIDUAL WEIGHT AND TRY AGAIN NEXT MOVE
C  SET SINKQ=0.0 SINCE NO MASS REMOVED
          IF (NPCELL(JS,IS,KS).EQ.0) THEN
cgzh debug
c         write(iouts,*) '***WARNING*** Residual due to NPCELL=0 at sink'
c         write(iouts,*) 'Sink cell (js,is,ks):',js,is,ks
c         write(iouts,*) 'Residual= ',SINKQ(js,is,ks)
            RESIDWT(JS,IS,KS)=SINKQ(JS,IS,KS)
            SINKQ(JS,IS,KS)=0.0
          ELSE
C  IF THERE ARE PARTICLES, BUT VOLUME ON ALL PARTICLES 
C    IS LESS THAN VOLUME TO BE REMOVED, REMOVE ALL VOLUME AND STORE RESIDUAL
            RESIDCL=SUMWT(JS,IS,KS)+SINKQ(JS,IS,KS)
            IF(RESIDCL.LT.0.0) THEN
cgzh debug
c         write(iouts,*) '***WARNING*** Residual due to not enough wt',
c     * ' on pts'
c         write(iouts,*) 'Sink cell (js,is,ks):',js,is,ks
c         write(iouts,*) 'Residual= ',RESIDCL
c         write(iouts,*) 'Tempwt= ',tempwt(js,is,ks)
c         write(iouts,*) 'Sumwt=  ',sumwt(js,is,ks)
c         write(iouts,*) 'Sinkq=  ',sinkq(js,is,ks)
               RESIDWT(JS,IS,KS)=RESIDCL
                   SINKQ(JS,IS,KS)=-SUMWT(JS,IS,KS)
            END IF
          END IF
C  IF ADDING WATER BACK IN (DUE TO LARGE SUMSGWT)
        ELSE IF (SINKQ(JS,IS,KS).GT.0.0) THEN
C     ALTER SUMSG ARRAYS TO SHOW THAT MASS DIDN'T LEAVE
          SGPTC(JS,IS,KS)=SUMSGMS(JS,IS,KS)/SUMSGWT(JS,IS,KS)
          SUMSGWT(JS,IS,KS)=SUMSGWT(JS,IS,KS)+SINKQ(JS,IS,KS)
          SUMSGMS(JS,IS,KS)=SUMSGMS(JS,IS,KS)+
     &                      SINKQ(JS,IS,KS)*SGPTC(JS,IS,KS)

C  IF NO PARTICLES TO ADD WEIGHT BACK ONTO,
C  ADD PARTICLE WITH WEIGHT EQUAL TO LEFTOVER AND CONC=PTS THAT LEFT
C  SET SINKQ=0 (ALL ACCOUNTED FOR VIA SUMSGWT)
          IF (NPCELL(JS,IS,KS).EQ.0) THEN
cgzh debug
c            WRITE(IOUTS,*) 
c            WRITE(IOUTS,*) 'Too much left on particle crossing subgrid'
c            WRITE(IOUTS,*) '...and NPCELL=0 at that sink!'
c            WRITE(IOUTS,*) 'Calling ADDPTS: SUMSGWT>SINKQ'
c            WRITE(IOUTS,*) 'Creating 1 new particle at (js,is,ks):'
c            WRITE(IOUTS,*) js,is,ks
c            WRITE(IOUTS,*) 'SUMSGWT= ',SUMSGWT(JS,IS,KS) 
c            WRITE(IOUTS,*) 'SINKQ (includes SUMSGWT) = ',SINKQ(JS,IS,KS)
c            WRITE(IOUTS,*) 'Setting SINKQ=0.0'
cgzh debug
cgzh varpt 
            IF(INIPDL.EQ.0.AND.INIPDA.EQ.0) THEN
              CALL ADDPTS(PC,PR,PL,PCONC,SGPTC(JS,IS,KS),IPTID,NPCELL,
     *          IBOUND,PNEWC,PNEWR,PNEWL,
     *          PTWT,SINKQ(JS,IS,KS),
     *          NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
C NEWPTS=2 means NPTPND=1
     *          2,NPMAX,
     *          JS,IS,KS,
     *          IOUTS,NPTM,WTFAC)
            ELSE
cgzh varpt 
cgzh set number of particles in each direction equal to 1 
              NPTCOL=1
              NPTROW=1
              NPTLAY=1
              CALL ADDPTSVAR(PC,PR,PL,PCONC,SGPTC(JS,IS,KS),
     *          NPCELL,PCORIG,PRORIG,PLORIG,
     *          PTWT,SINKQ(JS,IS,KS),
     *          NSCOL,NSROW,NSLAY,
     *          NPMAX,NPTCOL,NPTROW,NPTLAY,
     *          JS,IS,KS,
     *          IOUTS,NPTM,IDIM,WTFAC)
            END IF
C           TRACK AMOUNT ADDED BY THIS PROCESS
C           TEMPWT ALREADY USED, SO PUT DIRECTLY IN SUMWT
cgzh srcfix2 changed these from tempwt due to moving tempwt addition above
            SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+SINKQ(JS,IS,KS)
            SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)
     &                       +SINKQ(JS,IS,KS)*SGPTC(JS,IS,KS)
            SINKQ(JS,IS,KS)=0.0
          END IF
        END IF
C  ADD MASS ALTERED BY SOURCES (OR "EXCESS" SINKS) TO SUM ARRAYS
cgzh srcfix2 moved above
        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
C  REINITIALIZE FOR USE WITH SINK LOOP
cgzh stet
        TEMPWT(JS,IS,KS)=0.0
        TEMPMS(JS,IS,KS)=0.0
C  SET NPCELL COUNTER FOR SINK LOOP
        NPSUM(JS,IS,KS)=NPCELL(JS,IS,KS)
C  COMPUTE FRACTION OF VOLUME TO BE ADDED BACK IN FROM PTS 'REMOVED' DUE
C    TO RELATIVELY LOW WTS
cgzh debug
        if(imov.eq.11.and.js.eq.1) then
          continue
        end if
        IF(NPCELL(JS,IS,KS).GT.0.AND.REMVWT(JS,IS,KS).GT.0.0) THEN
          FRAC1(JS,IS,KS)=REMVWT(JS,IS,KS)/NPCELL(JS,IS,KS)
          FRAC2(JS,IS,KS)=REMVMS(JS,IS,KS)/NPCELL(JS,IS,KS)
        END IF
cgzh debug output wtcheck
c      if(js.eq.14.and.is.eq.1.and.ks.eq.1) then
c       write(87,*) SRCWT(JS,IS,KS),SRCMS(JS,IS,KS),NPNEW(JS,IS,KS)
c      end if
cgzh
 1600 CONTINUE
C
cgzh debug
C  UPDATE ACCUMULATORS AFTER THIS LOOP: srcfix CODE IN 1600 LOOP ALTERS 
C  JS+1, IS+1 etc SO UNABLE TO DO UPDATE CORRECTLY IN SAME LOOP
c      DO 1650 KS=1,NSLAY
c      DO 1650 IS=1,NSROW
c      DO 1650 JS=1,NSCOL
C  ADD MASS ALTERED BY SOURCES (OR "EXCESS" SINKS) TO SUM ARRAYS
c        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
C  REINITIALIZE FOR USE WITH SINK LOOP
c        TEMPWT(JS,IS,KS)=0.0
c        TEMPMS(JS,IS,KS)=0.0
c 1650      CONTINUE

cgzh debug output
c        write(iouts,*) 'loop 1600 sumwt,remvwt 111',sumwt(1,1,1),
c     *  remvwt(1,1,1)
cgzh debug skipping output
c      go to 3103
cgzh debug pt loop
C21 loop over particles to compute cumulative mass, mass per cell and cumulative volume.
      CUMMASS=0.0
      CUMCELL=0.0
      CUMCELL2=0.0
      SMCELL=0.0
      sumvol2=0.0
      DO 180 IP=1,NPTM
        IF(PC(IP).EQ.0.0) GO TO 180
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
C       ---SUM MASS---
        IF(IBOUND(J,I,K).NE.0) THEN
          CUMMASS=CUMMASS+PTWT(IP)*PCONC(IP)
          SMCELL(JS,IS,KS)=SMCELL(JS,IS,KS)+PTWT(IP)*PCONC(IP)
          sumvol2=sumvol2+ptwt(ip)
        END IF
  180 CONTINUE
C22 loop over cells to compute cumulative mass, and cumulative volume.  
      SUMVOL=0.0
      DO 3102 KS=1,NSLAY
      DO 3102 IS=1,NSROW
      DO 3102 JS=1,NSCOL
       SUMVOL=SUMVOL+SUMWT(JS,IS,KS)
       CUMCELL=CUMCELL+SMCELL(JS,IS,KS)
       CUMCELL2=CUMCELL2+SUMMASS(JS,IS,KS)
 3102 CONTINUE
 3103 CONTINUE
cgzh debug output
c      WRITE(IOUTS,*) 'After 1st cell loop, wton=  ',wton
c      WRITE(IOUTS,*) 'After 1st cell loop, wtof= ',wtof
c      WRITE(IOUTS,*) 'After 1st cell loop, SUMVOL=  ',SUMVOL
c      WRITE(IOUTS,*) 'After 1st cell loop, SUMVOL2= ',SUMVOL2
      sumvol=0.0
      sumvol2=0.0
c      WRITE(IOUTS,*) 'After 1st cell loop,CUMCELL2= ',CUMCELL2
c      WRITE(IOUTS,*) 'After 1st cell loop, CUMMASS= ',CUMMASS
cgzh debug output
c      WRITE(IOUTS,*) 'After 1st cl loop, sumwt 411 = ',sumwt(4,1,1)
cgzh debug output
C
C
C23 GWT----PRINT PARTICLE OBSERVATION DATA
C       SEND IN SINKQ, WILL PRINT PTS NOT IN MNWS HERE (CALL ABOVE HANDLES
C       MNW SINKS)
      IF (INUNITPTOB.GT.0)
     *  CALL SPTOB5O(SUMTCH,IPTOBLST,IMOV,PC,PR,PL,NPMAX,NP,TIMV,
     *        NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,NUMPTOB,NUMPTOB_MNW,
     *        IPTOBMNW,NPCELL,PTWT,PCONC,SUMWT,SNKFLO,SNKMNW,0,0,0.0,
     *        0.0,0,0)
C
C     ****************************************************************
C
C24        ---THIRD LOOP OVER PARTICLES: ADJUST WEIGHTS DUE TO SINKS---
cgzh srcnew
C                 MIX MASS AT SOURCES WITH NEW PARTICLES
cgzh debug
cgzh debug
      addwt=0.0
      addms=0.0
      sumdum=0.0
      SMCELL=0.0
cgzh debug output
c      write(iouts,*) 'nptm in 3rd loop=',nptm
      DO 2590 IP=1,NPTM
cgzh debug output
c      if(ip.eq.363) then
c      write(iouts,*) 'pc(363)=',pc(ip)
c      write(iouts,*) 'pcorig(363)=',pcorig(ip)
c      end if
        NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
        IF(NEWC.LE.0.0D0) GO TO 2590
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
        J=INT(NEWC+0.5D0)
        JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
        NEWR=PR(IP)
        IF(NEWR.LT.0.0D0) THEN
           NEWR=-NEWR
        END IF
        I=INT(NEWR+0.5D0)
        IS=I-ISROW1+1
        NEWL=PL(IP)
        K=INT(NEWL+0.5D0)
        KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
        IF(IBOUND(J,I,K).EQ.0)       THEN
          WRITE(IOUTS,*) '**ERROR*** PARTICLE IN NO-FLOW CELL'
          WRITE(IOUTS,*) 'LOOP 3'
          WRITE(IOUTS,*) 'J,I,K,IP,PTWT',J,I,K,IP,PTWT(IP)
          STOP '**ERROR*** PARTICLE IN NO-FLOW CELL'        
c        GO TO 2590
        END IF
cgzh srcnew
C FULLY MIX MASS AT SOURCE BY EQUALIZING PCONC
C     ON ALL PARTICLES IN CELL 
C ONLY MIX STRONG SOURCES
cgzh mixing means summing total mass and weight and using that to recalculate
cgzh equal particle concentrations
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) then
c       write(iouts,*) 'A 3rd loop ip,ptwt',ip,ptwt(ip)
c      end if
        if(js.eq.2) then
          continue
        end if
        IF(IGENPT(JS,IS,KS).EQ.1) THEN
C ONLY IF SRCWT>0, THIS TAKES CARE OF SOURCE WITH ZERO PTS WHEN SRCWT CALC'D
          IF(SRCWT(JS,IS,KS).GT.0.0.AND.NPCELL(JS,IS,KS).GT.1) THEN 
            PCONC(IP)=SRCMS(JS,IS,KS)/SRCWT(JS,IS,KS)
cgzh debug also distribute weight?
            PTWT(IP)=SRCWT(JS,IS,KS)/NPSUM(JS,IS,KS)
cgzh debug
c       write(iouts,*) 'SOURCE CELL PARTICLE MIXED'
c       write(iouts,*) 'js,is,ks,npcell',js,is,ks,npcell(js,is,ks)
c       write(iouts,*) 'srcwt,srcms',SRCWT(JS,IS,KS),SRCMS(JS,IS,KS)
c       write(iouts,*) 'sumwt-sinkq',sumwt(JS,IS,KS)-sinkq(js,is,ks)
          ENDIF
        END IF
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) then
c       write(iouts,*) 'after srcwt 3rd loop ip,ptwt',ip,ptwt(ip)
c      end if
C  IF PARTICLES WERE REMOVED FROM THIS CELL IN 1590 LOOP,
C    REDISTRIBUTE THE WEIGHT AND SOLUTE MASS ON REMAINING PARTICLES IN CELL
        IF (REMVWT(JS,IS,KS).GT.0.0) THEN
C  CALCULATE NEW PARTICLE CONCENTRATION: SUM EXISTING AND NEW MASS AND
C    DIVIDE BY SUM OF EXISTING AND NEW VOLUME
C
          PCONC(IP)=(PCONC(IP)*PTWT(IP)+FRAC2(JS,IS,KS))/
     *             (PTWT(IP)+FRAC1(JS,IS,KS))
C
C  ADJUST VOLUME OF PARTICLE BY FRACTION OF VOLUME
          PTWT(IP)=PTWT(IP)+FRAC1(JS,IS,KS)
cgzh debug output
c         write(IOUTS,*) 'Particle',IP,' gained wt in',js,is,ks
c         write(IOUTS,*) 'NPSUM=',NPSUM(JS,IS,KS),'
c     *         REMVWT=',REMVWT(JS,IS,KS)
c         write(IOUTS,*) 'FRACWT=',FRACWT
c         write(IOUTS,*) 'PTWT after =',ptwt(ip)
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) then
c       write(iouts,*) 'B 3rd loop ip,ptwt',ip,ptwt(ip)
c      end if
        END IF
C
C INITIALIZE FRACWT FOR THIS SECTION; USED IN LOGIC CHECK
        FRACWT=0.0
C  UPDATE WEIGHTS OF PTS. IN SINK CELLS
cgzh debug
        if(js.eq.3.and.is.eq.1) then
          continue
        end if
        IF (SINKQ(JS,IS,KS).NE.0.0) THEN
cgzh debug
          if(imov.eq.7) then
            continue
          endif
C  COMPUTE FRACTION OF VOLUME TO BE REMOVED (OR ADDED IF SINKQ>0) 
cgzh debug
C ONLY REMOVE WEIGHT IF THE CELL HAS A SINK (CHECK SINKRT)
C SINKQ MAY ONLY HAVE RESIDWT IN IT, DON'T REMOVE IN THIS CASE
c        IF(SINKQ(JS,IS,KS).LT.0.0) THEN
C ONLY IF SUMWT>0
C  IF ET IS ACTIVE, GET ETFLO
          ETFLO=0.0
          IF(INEVT.GT.0) THEN
            KEV=EVTFLO(JS,IS,2)
            IF(KEV.EQ.K) ETFLO=EVTFLO(JS,IS,1)*TIMV
          END IF
          SINKRT=SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS)+ETFLO
          IF(SINKQ(JS,IS,KS).LT.0.0.AND.SINKRT.LT.0.0) THEN
            IF(SUMWT(JS,IS,KS).GT.0.0) THEN
              FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*SINKQ(JS,IS,KS)
            END IF
          ELSEIF(SINKQ(JS,IS,KS).GT.0.0) THEN
C ADDING VOLUME
cgzh debug SINK FLAG
            SNKFLAG=1
            IF(SNKFLAG.EQ.0) THEN
cgzh this option seems bad as it would not add water if none was there
cgzh (in fact if snkflag=0 is implemented we need a check on sumwt>0 to avoid divby0)
              FRACWT=(PTWT(IP)/SUMWT(JS,IS,KS))*SINKQ(JS,IS,KS)
            ELSE
              FRACWT=SINKQ(JS,IS,KS)/REAL(NPSUM(JS,IS,KS))
            ENDIF
          END IF        
C  IF ADDING WATER, ADJUST MASS BASED ON PTS-LEAVING CONC
          IF (FRACWT.GT.0.0) THEN
            FRACMS=FRACWT*SGPTC(JS,IS,KS)
cgzh debug
cgzh debug IMPORTANT if mass being added on, should not be included in SNKMSOUT! 
c          write(IOUTS,*)'Mass added back onto pt at js,is,ks: ',js,is,ks
c          write(IOUTS,*) 'fracms=',fracms
cgzh debug
C  CALCULATE NEW PARTICLE CONCENTRATION: SUM EXISTING AND NEW MASS AND
C    DIVIDE BY SUM OF EXISTING AND NEW VOLUME   
C
cgzh debug div by zero check
            WT1=PTWT(IP)+FRACWT
            IF(WT1.GT.0) THEN
              PCONC(IP)=(PCONC(IP)*PTWT(IP)+FRACMS)/(WT1)
            ELSE
              PCONC(IP)=0.0
            END IF
            PTWT(IP)=PTWT(IP)+FRACWT
C
C  TRACK MASS ADDED TO CELL
            TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+FRACMS
C  IF REMOVING WATER, ADJUST MASS BASED ON PT CONC, NO NEED TO ADJUST PCONC 
C  UNLESS ET ACTIVE
          ELSE
C  IF NO ET OR AGE PACKAGE ON, CALC MASS AND ADJUST PTWT NORMALLY
C  IF AGE PACKAGE IS ACTIVE, DO NOT ALLOW ET TO CONCENTRATE AGE MASS 
            IF(ETFLO.EQ.0.OR.INAGE.GT.0) THEN
              TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+FRACWT*PCONC(IP)
              PTWT(IP)=PTWT(IP)+FRACWT
            ELSEIF(ETFLO.EQ.SINKQ(JS,IS,KS)) THEN
C  IF ONLY ET IN CELL, MASS OUT BASED ON ET SPECIFIER (CURRENTLY HARDWIRED TO 0)
C  RECALC PCONC BASED ON NEW WEIGHT IN CELL
              ETCONC=0.0
              IF(INAGE.GT.0) THEN 
                ETCONC=PCONC(IP)
              ENDIF
              TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+FRACWT*ETCONC
C  NEW PCONC IS OLD MASS DIVIDED BY NEW WEIGHT (THAT IS, REMOVE THE WATER, LEAVE THE SOLUTE)
              TEMPPTWT=PTWT(IP)+FRACWT
              PCONC(IP)=(PTWT(IP)*PCONC(IP))/TEMPPTWT
                  PTWT(IP)=TEMPPTWT
            ELSE          
C  IF ET AND OTHER SINKS IN CELL, APPLY 1/2 OF ET TO GET INTERMEDIATE PCONC, USE THAT
C  PCONC TO DETERMINE MASS LEAVING OTHER SINKS, THEN APPLY REMAINING ET
C
C  CALCULATE FRACTION OF ET IN TOTAL SINKQ
              FRACET=(ETFLO/SINKQ(JS,IS,KS))*FRACWT
C  CALCULATE FRACTION OF NON-ET SINKS IN TOTAL SINKQ
              FRACSNK=FRACWT-FRACET
C  INTERMEDIATE WT: REMOVE 1/2 ET
              TEMPPTWT=PTWT(IP)+(0.D5*FRACET)
C  INTERMEDIATE PCONC IS OLD MASS DIVIDED BY NEW WEIGHT (THAT IS, REMOVE THE WATER, LEAVE THE SOLUTE)
              if (TEMPPTWT.ne.0) then
                PCONC(IP)=(PTWT(IP)*PCONC(IP))/TEMPPTWT
              endif
C  BASE MASS LEAVING OTHER SINKS ON INTERMEDIATE PCONC
              TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+FRACSNK*PCONC(IP)
C  REMOVE NON-ET SINK WT
              PTWT(IP)=TEMPPTWT+FRACSNK
C  INTERMEDIATE WT: REMOVE OTHER 1/2 ET
              TEMPPTWT=PTWT(IP)+(0.D5*FRACET)
C  FINAL PCONC IS OLD MASS DIVIDED BY NEW WEIGHT (THAT IS, REMOVE THE WATER, LEAVE THE SOLUTE)
              if (TEMPPTWT.ne.0) then
                PCONC(IP)=(PTWT(IP)*PCONC(IP))/TEMPPTWT
              endif
C  SET FINAL WT
                PTWT(IP)=TEMPPTWT
            END IF
          END IF
C  ADJUST VOLUME OF PARTICLE BY FRACTION OF VOLUME
C    NEW WEIGHT IS (OLD WEIGHT + FRACTION OF OLD WEIGHT)
cgzh debug output
c      if(ip.eq.8857.and.imov.eq.36) then 
c        write(iouts,*) 'C: ptwt 8857=',ptwt(8857)
c      end if
C   FIXED FOR ROUNDOFF ERROR CAUSING SMALL NEGATIVE WEIGHTS
          IF(PTWT(IP).LT.(-0.001)) THEN
            WRITE(IOUTS,*) '***ERROR*** PTWT << 0, PTWT,IP',PTWT(ip),ip
          END IF
          if(ptwt(ip).lt.0.0) then
cgzh debug output?
c          write(iouts,*) '**WARNING** PTWT < 0, ip,ptwt', ip,ptwt(ip)
            PTWT(IP)=0.0
          end if
cgzh debug remove pts less than "small" fraction of critwt
          IF(PTWT(IP).LE.(CRITWT(JS,IS,KS)*0.0000001)) THEN
cgzh debug output
c            WRITE(IOUTS,*) 'Zero-weight particle removed at:',js,is,ks
              wtof=wtof+ptwt(ip)
c      write(iouts,*) 'wtof=',wtof
c      write(iouts,*) 'srcflo',srcflo(js,is,ks)
      
cgzh debug output
c      if (ip.eq.585902.or.ip.eq.585884) then
c       write(iouts,*) 'debugging - PT PUT IN LIMBO, low weight, ip=',ip
c       write(iouts,*) 'ptwt,sumwt',ptwt(ip),sumwt(js,is,ks)
c       write(iouts,*) 'pc,pr,pl',pc(ip),pr(ip),pl(ip)
c       write(iouts,*) 'js,is,ks',js,is,ks
c       write(iouts,*) 'npcell,srcflo,snkflo',npcell(js,is,ks)
c     * ,srcflo(js,is,ks),snkflo(js,is,ks)
c      end if
              NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)-1
              PC(IP)=0.0
              PR(IP)=0.0
              PL(IP)=0.0
              PCONC(IP)=0.0
              PTWT(IP)=0.0
C  
cgzh debug  only create limbo locations if not a new particle
              IF(IP.LE.NP) THEN
                DO 2570 ID=1,NLIMBO
                  IF(LIMBO(ID).GT.0) GO TO 2570
                  LIMBO(ID)=IP
                  GO TO 2580
 2570           CONTINUE
              END IF
 2580         CONTINUE
          END IF

C   TRACK VOLUME REMOVED FROM CELL
          TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+FRACWT
C end sinkq<0        
        END IF
cgzh debug output
        SMCELL(JS,IS,KS)=SMCELL(JS,IS,KS)+PTWT(IP)*PCONC(IP)
cgzh debug
        if(js.eq.1.and.is.eq.1.and.ks.eq.1)then
          sumdum=sumdum+ptwt(ip)
        end if
cgzh debug output
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1) then
c       write(iouts,*) 'B 3rd loop ip,ptwt',ip,ptwt(ip)
c      end if
cgzh debug zinnpt
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.17) then
c       write(iouts,*) 'pt in 17: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.18) then
c       write(iouts,*) 'pt in 18: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
        sumvol2=sumvol2+ptwt(ip)
 2590 CONTINUE
cgzh debug output
c        write(iouts,*) 'loop 3 pcorig(121)',pcorig(121)
c        write(iouts,*) 'loop 3 pc(121)',pc(121)
cgzh debug output
c         WRITE(IOUTS,8330)
c         WRITE(IOUTS,*) (PTWT(IP),IP=1,NPTM)
c 8330 FORMAT (/,2X,'PTWT AT END OF THIRD PT LOOP',/)
c      write(iouts,*)
C25 loop over cells to compute cumulative mass, and cumulative volume.
      CUMCELL=0.0
      SUMVOL=0.0
      DO 4112 KS=1,NSLAY
      DO 4112 IS=1,NSROW
      DO 4112 JS=1,NSCOL
        CUMCELL=CUMCELL+SMCELL(JS,IS,KS)
        sumvol=sumvol+sumwt(js,is,ks)
c        write(iouts,*) 'tempwt,sinkq js is ks',
c     & tempwt(JS,IS,KS),sinkq(JS,IS,KS)
 4112 CONTINUE
cgzh debug output
cgzh debug zinnpt
c      if(imov.eq.68) then
c       write(iouts,*) 'After 3rd ploop SUMWT row 4=',SUMWT(4,1,1)
c      end if
c      WRITE(IOUTS,*) 'After 3rd pt loop, SUMVOL=   ',SUMVOL
c      WRITE(IOUTS,*) 'After 3rd pt loop, SUMVOL2=   ',SUMVOL2
c      WRITE(IOUTS,*) 'After 3rd pt loop, CUMCELL= ',CUMCELL
cgzh debug output
cgzh debug output fort.72
cgzh debug output
c      WRITE(IOUTS,*) 'After 3rd pt loop, sumwt 11 18 = ',sumwt(11,18,1)
c      write(72,*) 'Sum of weights at 1,1=',SRCWT(1,1,1)
c      WRITE(IOUTS,*) 'After 3rd pt loop, SUMWT 17 7=  ',SUMWT(17,7,1)
c     WRITE(IOUTS,*) 'After 3rd pt loop, sumdum = ',sumdum
c      WRITE(IOUTS,*) 'After 3rd pt loop, npcel 17 7 = ',NPCELL(17,7,1)
cgzh debug output
c      WRITE(IOUTS,*) 'After 3rd pt loop, sumwt 1 1 = ',sumwt(1,1,1)
C     ****************************************************************
C
cgzh debug skipping 4th loop
      GO TO 8888
C26        ---FOURTH LOOP OVER PARTICLES: REMOVE ZERO-WEIGHT PARTICLES---
C   ACTUALLY, IF WEIGHT IS VERY SMALL RELATIVE TO INITIAL PT WTS
C
cgzh debug
      sumdum=0.0
      DO 3590 IP=1,NPTM
        NEWC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
        IF(NEWC.LE.0.0D0) GO TO 3590
C     ***************************************************************
C           ---COMPUTE NEW LOCATION---
        J=INT(NEWC+0.5D0)
        JS=J-ISCOL1+1
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
        NEWR=PR(IP)
        IF(NEWR.LT.0.0D0) THEN
           NEWR=-NEWR
        END IF
        I=INT(NEWR+0.5D0)
        IS=I-ISROW1+1
        NEWL=PL(IP)
        K=INT(NEWL+0.5D0)
        KS=K-ISLAY1+1
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
        IF(IBOUND(J,I,K).EQ.0)       THEN
          WRITE(IOUTS,*) '**ERROR*** PARTICLE IN NO-FLOW CELL'
          WRITE(IOUTS,*) 'LOOP 3'
          WRITE(IOUTS,*) 'J,I,K,IP,PTWT',J,I,K,IP,PTWT(IP)
          STOP '**ERROR*** PARTICLE IN NO-FLOW CELL'      
c        GO TO 3590
        END IF
cgzh orig          IF(PTWT(IP).LE.0.0) THEN
cgzh debug remove pts less than "small" fraction of critwt
        IF(PTWT(IP).LE.(CRITWT(JS,IS,KS)*0.0000001)) THEN
cgzh debug output
c            WRITE(IOUTS,*) 'Zero-weight particle removed at:',js,is,ks
c            write(iouts,*) 'ip,ptwt,npcell,sumwt'
c            write(iouts,*) ip,ptwt(ip),npcell(js,is,ks),sumwt(js,is,ks)
cgzh debug output
c      if(ip.eq.8857.and.imov.eq.36) then 
c        write(iouts,*) 'D: ptwt 8857=',ptwt(8857)
c      end if
c      if(ptwt(ip).lt.0.0) then
c        write(iouts,*) '***ERROR*** PTWT < 0!!'
c      end if
            NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)-1
            PC(IP)=0.0
            PR(IP)=0.0
            PL(IP)=0.0
            PCONC(IP)=0.0
            PTWT(IP)=0.0
C  
c            DO 2570 ID=1,NLIMBO
c              IF(LIMBO(ID).GT.0) GO TO 2570
c              LIMBO(ID)=IP
cgzh debug output
c      if(ip.eq.121) then
c      write(iouts,*) 'debugging - PT PUT IN LIMBO, zero-weight, ip=',ip
c      end if
c              GO TO 2580
c 2570       CONTINUE
c 2580 CONTINUE
        END IF
cgzh debug output
c       SMCELL(JS,IS,KS)=SMCELL(JS,IS,KS)+PTWT(IP)*PCONC(IP)
cgzh debug output
c      if(js.eq.140.and.ks.eq.2.and.imov.gt.5670) then
c       if(pconc(ip).gt.1e-3.or.pconc(ip).lt.-1e-3) 
c     * write(iouts,*) 'MOVE ptwt*pconc=',
c     * ptwt(ip)*pconc(ip)
c       if(pconc(ip).gt.1e-5.or.pconc(ip).lt.-1e-5) then
c       write(iouts,*) 'move b: ip,ptwt,pconc'
c       write(iouts,*) ip,ptwt(IP),pconc(ip)
c       endif
c      end if
c      if(js.eq.17.and.is.eq.7.and.ks.eq.1.and.imov.eq.10) then
c      sumdum=sumdum+ptwt(ip)
c      end if
cgzh debug end
cgzh debug zinnpt
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.17) then
c       write(iouts,*) 'pt in 17: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
c      if(imov.eq.68.and.js.eq.30.and.ks.eq.18.
c     &   and.is.eq.18) then
c       write(iouts,*) 'pt in 18: PC,PR,PL',pc(ip),pr(ip),pl(ip)      
c      end if
 3590 CONTINUE
cgzh debug output
c        write(iouts,*) 'loop 4 pcorig(121)',pcorig(121)
c        write(iouts,*) 'loop 4 pc(121)',pc(121)
 8888 continue
C27 loop over cells to compute cumulative mass, and cumulative volume.
      SUMVOL=0.0
      sumtemp=0.0
      DO 4102 KS=1,NSLAY
      DO 4102 IS=1,NSROW
      DO 4102 JS=1,NSCOL
       CUMCELL=CUMCELL+SUMMASS(JS,IS,KS)
       SUMVOL=SUMVOL+SUMWT(JS,IS,KS)
       sumtemp=sumtemp+tempwt(js,is,ks)
 4102 CONTINUE
cgzh debug output
cgzh debug zinnpt
c      if(imov.eq.68) then
c       write(iouts,*) 'After 4th ploop SUMWT row 4=',SUMWT(4,1,1)
c      end if
c      WRITE(IOUTS,*) 'After 4th pt loop, CUMCELL= ',CUMCELL
c      WRITE(IOUTS,*) 'After 4th pt loop, SUMVOL=   ',SUMVOL
c      WRITE(IOUTS,*) 'After 4th pt loop, sumtempwt=',sumtemp
c      WRITE(IOUTS,*) 'After 4th pt loop, addms=',addms
c      WRITE(IOUTS,*) 'addwt - remwt=',addwt-sumrem
c      WRITE(IOUTS,*) 'addms - remms=',addms-sumremms
cgzh debug output
cgzh debug output
c      WRITE(IOUTS,*) 'After 4th pt loop, SUMWT 17 7=  ',SUMWT(17,7,1)
c      WRITE(IOUTS,*) 'After 4th pt loop, sumdum = ',sumdum
C
C  LOOP OVER CELLS TO:
C  UPDATE SUM OF WEIGHTS AND MASS IN CELLS
C  TRACK MASS LEAVING SINKS FOR MASS BALANCE
cgzh drycell
C28  ADD NEW PARTICLES DUE TO DRY CELLS
cgzh debug 
      remsinks=0.0
      DO 2600 KS=1,NSLAY           
      DO 2600 IS=1,NSROW
      DO 2600 JS=1,NSCOL
C
cgzh drycell
        IF(TEMPDRYWT(JS,IS,KS).GT.0.0) THEN
c            WRITE(IOUTS,*) 
c            WRITE(IOUTS,*) 'Calling ADDPTS: DRYWT>0'
c            WRITE(IOUTS,*) 'Creating 1 new particle at (js,is,ks):'
c            WRITE(IOUTS,*) js,is,ks
          C1=TEMPDRYMS(JS,IS,KS)/TEMPDRYWT(JS,IS,KS)
C RBW Need to check if NPCELL is updated correctly if there are already
C one or more particles in the cell to which a particle is being added.          
          CALL ADDPTS(PC,PR,PL,PCONC,C1,IPTID,NPCELL,
     *          IBOUND,PNEWC,PNEWR,PNEWL,
     *          PTWT,TEMPDRYWT(JS,IS,KS),
     *          NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
C NEWPTS=2 means NPTPND=1
     *          2,NPMAX,
     *          JS,IS,KS,
     *          IOUTS,NPTM,WTFAC)
          TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+TEMPDRYWT(JS,IS,KS)
          TEMPMS(JS,IS,KS)=TEMPMS(JS,IS,KS)+TEMPDRYMS(JS,IS,KS)
          NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)+1
        END IF
cgzh drycell


        SUMWT(JS,IS,KS)=SUMWT(JS,IS,KS)+TEMPWT(JS,IS,KS)
C  ADD WEIGHTS AND MASS ALTERED BY SINKS TO "SUM" ARRAYS
        SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+TEMPMS(JS,IS,KS)
cgzh debug output
        if(tempms(js,is,ks).ne.0.0) then
         remsinks=remsinks+TEMPMS(JS,IS,KS)
        end if
C TRACK MASS LEAVING SUBGRID OR THROUGH SINKS
C   IF MASS WAS ADDED BACK IN, THIS IS SKIPPED (SINKQ GE 0)
C FRACSG = RATIO OF SINKQ ATTRIBUTED TO SUBGRID BOUNDARY
C FRACSNK = RATIO OF SINKQ ATTRIBUTED TO FLUID SINKS
        SINKRT=SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS)
C 9/4/04 ADDED CHECK FOR SINK IN CELL (IN LATER SP'S, SINKQ NE 0 POSSIBLE
C   DUE TO RESIDUALS)
        IF(SINKRT.LT.0.0.AND.SINKQ(JS,IS,KS).LT.0.0) THEN
C SUBTRACT SUMSGWT FROM BDYSNK FOR ITEMIZATION
C IF MORE LEFT SUMSGWT THAN BDYSNK (the MIN below), SET TEMP TERM TO ZERO 
C   SO THE REST WILL GO OUT SNKFLO IN MASS BALANCE 
          BDYTEMP=MIN(((BDYSNK(JS,IS,KS)*TIMV)-SUMSGWT(JS,IS,KS)),0.0)
corig          FRACSG=BDYTEMP/((TIMV*(SNKFLO(JS,IS,KS)+BDYSNK(JS,IS,KS)))
          FRACSG=BDYTEMP/((TIMV*(SINKRT))
     *          -SUMSGWT(JS,IS,KS))
          FRACSNK=1.0-FRACSG
C SGMSOUT  = Mass removed from particles in this cell due to flux out via subgrid
C SNKMSOUT = Mass removed this cell by fluid sinks.
C   NEGATIVE VALUE SIGNIFIES LEAVING AQUIFER
c          SGMSOUT(JS,IS,KS)=FRACSG*TEMPMS(JS,IS,KS)
c          SNKMSOUT(JS,IS,KS)=FRACSNK*TEMPMS(JS,IS,KS)
cgzh debug apply Rf to sink mass budget numbers
          SGMSOUT(JS,IS,KS)=FRACSG*TEMPMS(JS,IS,KS)*RF(JS,IS,KS)
          SNKMSOUT(JS,IS,KS)=SNKMSOUT(JS,IS,KS)+
     &      FRACSNK*TEMPMS(JS,IS,KS)*RF(JS,IS,KS)
cgzh debug output--this can trigger when conc<0.0
c      if(tempms(js,is,ks).gt.0.0) then
c       write(iouts,*) 'tempms>0 at sink??'
c      end if
c
c
c
c      if(FRACSNK.eq.1.0) then
c       write(iouts,*) 'fracsnk=1 cell js,is,ks=',
c     * FRACSNK,JS,IS,KS
c       write(iouts,*) 'bdysnk,snkflo',
c     * bdysnk(js,is,ks),snkflo(js,is,ks)
c      end if
c      end if
cgzh debug end
        END IF
C  
 2600 CONTINUE
cgzh debug 
c      write(iouts,*) 'remsinks=',remsinks
cgzh debug skipping output
c      go to 5103
cgzh debug pt loop
C29 loop over particles to compute cumulative mass, mass per cell and cumulative volume.
      CUMMASS=0.0
      CUMCELL=0.0
      CUMCELL2=0.0
      SMCELL=0.0
      sumvol2=0.0
      summass2=0.0
      kpt=0
      DO 1800 IP=1,NPTM
        IF(PC(IP).EQ.0.0) GO TO 1800
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
        IF(IBOUND(J,I,K).NE.0) THEN
cgzh debug
         CUMMASS=CUMMASS+PTWT(IP)*PCONC(IP)
         SMCELL(JS,IS,KS)=SMCELL(JS,IS,KS)+PTWT(IP)*PCONC(IP)
         sumvol2=sumvol2+ptwt(ip)
        END IF
 1800 CONTINUE
C30 loop over cells to compute cumulative mass, and cumulative volume.
      SUMVOL=0.0
      DO 5102 KS=1,NSLAY
      DO 5102 IS=1,NSROW
      DO 5102 JS=1,NSCOL
        CUMCELL=CUMCELL+SMCELL(JS,IS,KS)
        CUMCELL2=CUMCELL2+SUMMASS(JS,IS,KS)
        SUMVOL=SUMVOL+SUMWT(JS,IS,KS)
 5102 CONTINUE
 5103 CONTINUE
c      WRITE(IOUTS,*) 'After 2nd cell loop, SUMVOL= ',SUMVOL
c      WRITE(IOUTS,*) 'After 2nd cell loop, SUMVOL2=',SUMVOL2
c      WRITE(IOUTS,*) 'After 2nd cell loop, CUMCELL= ',CUMCELL
c      WRITE(IOUTS,*) 'After 2nd cell loop,CUMCELL2= ',CUMCELL2
c      WRITE(IOUTS,*) 'After 2nd cell loop, CUMMASS= ',CUMMASS
cgzh debug output
c      WRITE(IOUTS,*) 'After 2nd cl loop, SUMWT 17 7=  ',SUMWT(17,7,1)
C
cgzh debug loop start
c      sumdum=0.0
c      DO 23590 IP=1,NPTM
c      NEWC=PC(IP)
c      J=INT(NEWC+0.5D0)
c      JS=J-ISCOL1+1
c      NEWR=PR(IP)
c      IF(NEWR.LT.0.0D0) THEN
c         NEWR=-NEWR
c      END IF
c      I=INT(NEWR+0.5D0)
c      IS=I-ISROW1+1
c      NEWL=PL(IP)
c      K=INT(NEWL+0.5D0)
c      KS=K-ISLAY1+1
cgzh debug output
c      if(js.eq.17.and.is.eq.7.and.ks.eq.1.and.imov.eq.10) then
c        write(iouts,*) 'pt in 17,7: ip, ptwt', ip, ptwt(ip)
c      end if
c23590      CONTINUE
cgzh debug loop end
c  begin lumping loops
C31  LUMP PARTICLES WITH LOW WEIGHTS IN CELLS THE HAD PTS CREATED BY SRCFIX
C  CRITWT IS FRACTION (REMCRIT) OF SUM OF PARTICLE WEIGHTS IN CELL  
C  CRITMS IS FRACTION (REMCRIT) OF SUM OF PARTICLE MASSES IN CELL  
C  (REMCRIT) --> SET BY USER
cgzh debug ilump=0 off, ilump=1 on
      ilump=1
cgzh debug temp iii
       iii=0

      if(ilump.gt.0) then
        IF(ISRCFIX.GT.0) THEN
C  Reset NPNEW (used to count pts removed from a cell) and others
          NPNEW=0
          TEMPWT=0
          TEMPMS=0
C  Set NPSUM = NPCELL now
          NPSUM=NPCELL
          DO 9580 IP=1,NPTM
            IF(PC(IP).EQ.0.0) GO TO 9580
            J=PC(IP)+0.5
            JS=J-ISCOL1+1
            I=ABS(PR(IP))+0.5
            IS=I-ISROW1+1
            K=PL(IP)+0.5
            KS=K-ISLAY1+1
            PTMASS=PTWT(IP)*PCONC(IP)
C  Lump is set when the srcfix code added pts to a cell
cgzh debug
c      if(is.eq.2.and.ks.eq.3) then
c       if (iii.lt.1) write(55,*) 'imov,critwt=',imov,CRITWT(JS,IS,KS)
c       iii=iii+1
c       write(55,*) 'ip,ptwt(ip)',iii,ip,ptwt(ip)
c      end if
cgzh debug
            if(imov.eq.256.and.js.eq.1.and.is.eq.1.and.ks.eq.2) then
              continue
            endif
            IF(LUMP(JS,IS,KS).GT.0) THEN
cgzh debug output
c         if(is.eq.2.and.ks.eq.6.and.imov.gt.28) then
c         write(55,*) 'imov,ptwt,v2',imov,ptwt(ip),V2(JS,IS,KS)
c         endif
              IF(PTWT(IP).LT.CRITWT(JS,IS,KS)
c     *     .AND.PTMASS.LT.CRITMS(JS,IS,KS)
c  allow particles to build up a little bit before lumping
     *         .AND.NPSUM(JS,IS,KS).GT.(NPORIG(JS,IS,KS)+8)) THEN
cgzh median only remove half of the particles
                  IVCOUNT(JS,IS,KS)=IVCOUNT(JS,IS,KS)+1
                  IF(MOD(IVCOUNT(JS,IS,KS),2).EQ.0) THEN
cgzh debug output
                    write(55,*) 'pt lumped in cell is,ks,imov',
     *                is,ks,imov
C  TRACK VOLUME AND MASS REMOVED
cgzh debug output
c           if(is.eq.2.and.ks.eq.3)
c     & write(55,*), 'pt removed iii=',iii
                    TEMPWT(JS,IS,KS)=TEMPWT(JS,IS,KS)+PTWT(IP)
                    TEMPMS(JS,IS,KS)=
     *                TEMPMS(JS,IS,KS)+(PTWT(IP)*PCONC(IP))
C  Sum coordinates (will use to get average coordinate for new pt)
                    AVC(JS,IS,KS)=AVC(JS,IS,KS)+PC(IP)
                    AVR(JS,IS,KS)=AVR(JS,IS,KS)+ABS(PR(IP))
                    AVL(JS,IS,KS)=AVL(JS,IS,KS)+PL(IP)
                    NPNEW(JS,IS,KS)=NPNEW(JS,IS,KS)+1
C  Remove particle
                    NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)-1
                    PC(IP)=0.0
                    PR(IP)=0.0
                    PL(IP)=0.0
                    PCONC(IP)=0.0
                    PTWT(IP)=0.0
C  
cgzh only create limbo locations if not a new particle
                    IF(IP.LE.NP) THEN
                      DO 9570 ID=1,NLIMBO
                        IF(LIMBO(ID).GT.0) GO TO 9570
                        LIMBO(ID)=IP
                        GO TO 9580
 9570                 CONTINUE
                    END IF
c  once lumped, undo flag
                    LUMP(JS,IS,KS)=0
                  END IF
              END IF
            END IF
 9580     CONTINUE
C          
          DO 1111 KS=1,NSLAY
          DO 1111 IS=1,NSROW
          DO 1111 JS=1,NSCOL
C NPNEW here is number of particles removed from a lumping cell
C If > 0, create a new particle
            IF(NPNEW(JS,IS,KS).GT.0) THEN
cgzh debug output
c      write(iouts,*) 'Lumping',NPNEW(JS,IS,KS),' particles at js,is,ks='
c     * ,js,is,ks
              NPTM=NPTM+1
cgzh debug also need "restart" code here once restart is available
              IF(NPTM.EQ.NPMAX) THEN
                write(*,*)
                STOP 'NPMAX EXCEEDED IN MOVEWT'
              END IF
C
              IPN=NPTM
              PC(IPN)=AVC(JS,IS,KS)/NPNEW(JS,IS,KS)
              PR(IPN)=AVR(JS,IS,KS)/NPNEW(JS,IS,KS)
              PL(IPN)=AVL(JS,IS,KS)/NPNEW(JS,IS,KS)
              PTWT(IPN)=TEMPWT(JS,IS,KS)
              PCONC(IPN)=TEMPMS(JS,IS,KS)/TEMPWT(JS,IS,KS)
              NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)+1
C     A LUMPED PARTICLE SHOULD NOT BE ADJUSTED FOR WATER TABLE, ALL LUMPED PTS
C       HAVE BEEN
cgzh debug output
c         write(IOUTS,*) 'New lumped particle',IPN,' at ',js,is,ks
c         write(IOUTS,*) 'PTWT=',PTWT(IPN),'PTMASS=',PTMASS
c         write(IOUTS,*) 'PC=',PC(IPN),'PR=',PR(IPN),'PL=',PL(IPN)
             END IF
 1111     CONTINUE
c end srcfix>0, lumping loops
        END IF
cgzh debug end if for ilump
      END IF
C 
C
cgzh ccbd
C
C32        ---APPLY CONSTANT-CONCENTRATION BOUNDARY PACKAGE---
C  THIS SECTION WILL CHANGE CONCENTRATION OF PARTICLES IN CONSTANT
C  CONCENTRATION BOUNDARY CELLS TO THE VALUE PRESCRIBED BY THE USER.
C  THE NET CHANGE IN MASS FOR EACH SUCH CELL IS COMPUTED AND REPORTED
C  IN THE MASS BALANCE AS EITHER ENTERING OR LEAVING THE SYSTEM.  
C
C  THIS SECTION ACCOUNTS FOR ADVECTION, DECAY, SOURCES, AND SINKS 
C
      IF(INCCBD.GT.0) THEN 
C  LOOP OVER PARTICLES
        DO 9880 IP=1,NPTM
          IF(PC(IP).EQ.0.0) GO TO 9880
          J=PC(IP)+0.5
          JS=J-ISCOL1+1
          I=ABS(PR(IP))+0.5
          IS=I-ISROW1+1
          K=PL(IP)+0.5
          KS=K-ISLAY1+1
C  IF CCBD IN THIS CELL, ADJUST CONC AND TRACK MASS
          CCCONC=CCBDY(JS,IS,KS)
          IF(CCCONC.GE.0.0) THEN
C  CALCULATE DIFFERENCE BETWEEN NEW AND OLD PARTICLE CONCENTRATION
            DIFFCONC=CCCONC-PCONC(IP)
C  CALCULATE CHANGE IN MASS DUE TO REASSIGNING PCONC
c
c    positive chngmass is mass added
c    negative chngmass is mass removed
            IF(DIFFCONC.NE.0.0) THEN
              CHNGMASS(JS,IS,KS)=CHNGMASS(JS,IS,KS)+DIFFCONC*PTWT(IP)
C  REASSIGN NEW PARTICLE CONCENTRATION
              PCONC(IP)=CCCONC
C  UPDATE SUMMASS
              SUMMASS(JS,IS,KS)=SUMMASS(JS,IS,KS)+DIFFCONC*PTWT(IP)
C
            END IF        
          END IF        
 9880   CONTINUE
C
C33  LOOP OVER CELLS TO SUM MASS ADDED OR REMOVED VIA CCBD PACKAGE
        TOTIN=0.D0
        TOTOUT=0.D0
        DO 2111 KS=1,NSLAY
        DO 2111 IS=1,NSROW
        DO 2111 JS=1,NSCOL
          DIFFMASS=CHNGMASS(JS,IS,KS)
          IF(DIFFMASS.GT.0.D0) THEN
            TOTIN=TOTIN+DIFFMASS
          ELSEIF(DIFFMASS.LT.0.D0) THEN
            TOTOUT=TOTOUT+DIFFMASS
          END IF
 2111   CONTINUE
C
C  UPDATE SBVL ARRAY
       SBVL(1,26)=SBVL(1,26)+TOTIN
       SBVL(2,26)=SBVL(2,26)+TOTOUT
C
      END IF
cgzh ccbd
C END IF for CCBD package
C
C34        ---INSERT PARTICLES INTO LIMBO LOCATIONS---
C  NPTM = NP WHEN NO NEW PARTICLES WERE CREATED
      IF(NPTM.EQ.NP) GO TO 620
C  START WITH LAST PARTICLE ADDED, WHICH HAS AN ID = NPTM
      IP=NPTM
cgzh debug output
c      if(imov.eq.25.and.ip.eq.17522) then
c      write(iouts,*) 'c: 17522 pc=',pc(17522),pconc(17522)
c      end if
      DO 595 IL=1,NLIMBO
        IPL=LIMBO(IL)
C  LOOK FOR NON-ZERO ID NUMBERS THAT WERE STORED IN LIMBO
        IF(IPL.EQ.0) GO TO 595
cgzh debug
c      if(ipl.eq.4500) then
c       if(imov.eq.10) then
c       write(iouts,*) 'pt being put into a limbo ipl,ip:', ipl,ip
c       write(iouts,*) 'pconc(ipl),pconc(ip):', pconc(ipl),pconc(ip)
c       write(iouts,*) 'at this ip: pc,pr,pl,wt', 
c     *  pc(ip),pr(ip),pl(ip),ptwt(ip)
c        end if  
C  IF THERE IS AN ID IN LIMBO, USE IT BY SAVING THE INFO
C    ABOUT THE PARTICLE THAT WAS CREATED IN THAT ID SPOT
c      IF(PC(IP).LT.0.4) THEN
c        write(iouts,*) 'PC<0.4, IP, IPL:',IP,IPL
c        GO TO 594
c      end if
        PR(IPL)=PR(IP)
        PR(IP)=0.0
        PC(IPL)=PC(IP)
        PC(IP)=0.0
        PL(IPL)=PL(IP)
        PL(IP)=0.0
        PCONC(IPL)=PCONC(IP)
        PCONC(IP)=0.0
        PTWT(IPL)=PTWT(IP)
        PTWT(IP)=0.0
cgzh varpt
        IF(INIPDL.EQ.0.AND.INIPDA.EQ.0) THEN
          IPTID(IPL)=IPTID(IP)
          IPTID(IP)=0
        ELSE
C  PARTICLE DEFINITION WITH IPDL OR IPDA PACKAGE
cgzh varpt  
cgzh store origin
          PCORIG(IPL)=PCORIG(IP)
          PRORIG(IPL)=PRORIG(IP)
          PLORIG(IPL)=PLORIG(IP)
          PCORIG(IP)=0.0
          PRORIG(IP)=0.0
          PLORIG(IP)=0.0
        END IF
c      if(ipl.eq.121) then
c        write(iouts,*) 'debugging - using 121 in limbo'
c      end if
        LIMBO(IL)=0
C  GO BACK TO THE NEXT PARTICLE THAT WAS CREATED
  594   IP=IP-1
C  IF IP=NP, NOT A PARTICLE THAT WAS CREATED THIS MOVE
        IF(IP.EQ.NP) GOTO 596
  595 CONTINUE
C  OUT OF LOOP, SET THE NUMBER OF PARTICLE SPOTS USED 
C    THIS WILL BE SMALLER IF SOME NEW PARTICLES WERE INSERTED INTO
C    LIMBO LOCATIONS
  596 NPTM=IP
  620 CONTINUE
cgzh debug output
c      if(imov.eq.25) then
c      write(iouts,*) 'd: 17522 pc=',pc(17522),pconc(17522)c
c      write(iouts,*) 'npcell',npcell(30,10,1)
c      end if
cgzh debug output
c        write(iouts,*) 'loop 5 pcorig(121)',pcorig(121)
c        write(iouts,*) 'loop 5 pc(121)',pc(121)
cgzh debug output
c      WRITE(IOUTS,*) 'After 5th loop, SUMWT 17 14=  ',SUMWT(17,14,2)
cgzh debug loop start
      sumdum=0.0
      DO 13590 IP=1,NPTM
        NEWC=PC(IP)
        J=INT(NEWC+0.5D0)
        JS=J-ISCOL1+1
        NEWR=PR(IP)
        IF(NEWR.LT.0.0D0) THEN
           NEWR=-NEWR
        END IF
        I=INT(NEWR+0.5D0)
        IS=I-ISROW1+1
        NEWL=PL(IP)
        K=INT(NEWL+0.5D0)
        KS=K-ISLAY1+1
        if(js.eq.1.and.is.eq.1.and.ks.eq.1) then
          sumdum=sumdum+ptwt(ip)
c          write(iouts,*) 'pt in 17,7: ip, ptwt', ip, ptwt(ip)
        end if
13590 CONTINUE
cgzh debug loop end
c      WRITE(IOUTS,*) 'After debug loop, sumdum = ',sumdum
c      WRITE(IOUTS,*) 'After debug loop, npcel 17 7 = ',NPCELL(17,7,1)

C        ---ADJUST NUMBER OF PARTICLES---
      NP=NPTM
cgzh debug
c      write(iouts,*) 'NP at end of move = ',NP
C
C35  COMPUTE NODE CONCENTRATIONS
C     ****************************************************************
C     ---CONC. CHANGE AT NODES DUE TO ADVECTION---
C                                  ... AND DECAY---
      NZERO=0
      NSKIP=0
      DO 90 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 90 IS=1,NSROW
      I=IS+ISROW1-1
      DO 90 JS=1,NSCOL
      J=JS+ISCOL1-1
        IF(IBOUND(J,I,K).EQ.0) GO TO 90
cgzh zerocell array is false if there are particles in cell, true in none
        ZEROCELL(JS,IS,KS)=0
        IF(NPCELL(JS,IS,KS).EQ.0) THEN
           NZERO=NZERO+1
           ZEROCELL(JS,IS,KS)=1
cgzh debug output
c       if(abs(cnold(js,is,ks)).gt.1e-6) then
c         write(iouts,*) 'zero pts in js is ks',js,is,ks
c         write(iouts,*) 'cnold=',cnold(js,is,ks)
c       end if
cgzh debug
C        IF NO PARTICLES IN CELL, THEN NO MASS IN CELL
           CONC(JS,IS,KS)=0.0
        ELSE
          IF (SUMWT(JS,IS,KS).GT.0.0) THEN
              CONC(JS,IS,KS)=SUMMASS(JS,IS,KS)/SUMWT(JS,IS,KS)
C IF SUMWT IS MUCH LESS THAN CELVOL, DO NOT USE FOR DISPERSION
C CALCULATION (ACCOMPLISHED BY USING NPCELL AS "FLAG")
              VOLRATIO=
     &            SUMWT(JS,IS,KS)/CELVOL(JS,IS,KS)
C  SET RATIO LIMIT FOR EXCLUDING CELLS (user input??)  
corig              RATIOLIMIT=0.10
              RATIOLIMIT=0.25
              IF(VOLRATIO.LT.RATIOLIMIT) THEN
! RBW Should there be a check to see if NPCELL(JS,IS,KS) is greater than zero before doing this?              
! RBW It seems like there is a danger that a negative value could be changed to a positive one.
                NPCELL(JS,IS,KS)=-NPCELL(JS,IS,KS)
                NSKIP=NSKIP+1
cgzh debug output
c                    WRITE (IOUTS,*) 'SUMWT/CELVOL= ',VOLRATIO,' (too low)'
c                    WRITE (IOUTS,*) 'SKIPPING DISP FOR CELL ',JS,IS,KS
              END IF
          ELSE
                  WRITE (IOUTS,*) 'SUMWT.LE.0.0 AT: ',JS,IS,KS,
     *       '; ASSUME NO CHANGE IN CONC. BY ADVECTION'
                  WRITE (IOUTS,*) 'SUMWT = ',SUMWT(JS,IS,KS)
                  WRITE (IOUTS,*) 'NPCELL = ',NPCELL(JS,IS,KS)
cgzh debug  this should not happen (see NPCELL check above)
c           STOP?
         END IF
        END IF
cgzh debug output   
c      if(js.eq.1.and.is.eq.1) then
c       if(summass(js,is,ks).gt.0.0) write(iouts,*) 'MOVE summass=   ',
c     * summass(js,is,ks),js,is,ks
c      end if
cgzh debug output
cgzh debug add clause for no particles (and hence no mass) in cell during 
cgzh last move
        IF(NPOLD(JS,IS,KS).GT.0) THEN
          CAVG(JS,IS,KS)=0.5*(CONC(JS,IS,KS)+CNOLD(JS,IS,KS))
        ELSE
c if no mass last move, use current advective conc
          CAVG(JS,IS,KS)=CONC(JS,IS,KS) 
        END IF
   90 CONTINUE
cgzh debug ouput
c      write(43,*)
c      IF(NSKIP.GT.0) THEN
cgzh debug output?  or keep this?
c      WRITE(IOUTS,*)'SUMWT/CELVOL RATIO TOO LOW, SKIPPING',NSKIP,'CELLS'
c      END IF
C  ADD NZERO TO NSKIP FOR CHECK ON RECALC OF CELL CONNECTIONS
c this is just below output on line above; can be moved if output is done away with
      NSKIP=NSKIP+NZERO
c      write(43,*) 'Cavg around 24,21, move=',imov
c      write(43,*) (CAVG(JS,20,1),JS=23,25)
c      write(43,*) (CAVG(JS,21,1),JS=23,25)
c      write(43,*) (CAVG(JS,22,1),JS=23,25)
c      write(43,*)
c      write(43,*)
c      write(43,*) 'npold around 24,21, move=',imov
c      write(43,*) (npold(JS,20,1),JS=23,25)
c      write(43,*) (npold(JS,21,1),JS=23,25)
c      write(43,*) (npold(JS,22,1),JS=23,25)
c      write(43,*)

C36        ---CHECK NUMBER OF CELLS VOID OF PTS.---
      IF(NZERO.GT.0) WRITE(IOUTS,290) NZERO
      
      ! Calling RefreshEmptyCells here leads to mass balance errors
      ! in the RMA model with 16 particles per cell.
C      IF(NZERO.GT.0) then
!          call RefreshEmptyCellsV4(PC,PR,PL,IPTID,PCONC,PTWT,
!     *      NPCELL, CONC, CNOLD, CELVOL, SUMWT, 
!     *      IBOUND,
!     *      PNEWC,PNEWR,PNEWL,
!     *      NEWPTS, ! number of initial particles per cell plus 1. 
!     *      IOUTS,  ! unit number of listing file for GWT
!     *      THCK, VC, VR, VL, POR,    
!     *      NPMAX, NP, NSCOL,NSROW,NSLAY,ISCOL1,ISROW1,ISLAY1, 
!     *      NCOL,NROW,NLAY,
!     *      NZERO, ! number of empty cells.
!     *      TIMV) ! length of transport time step.
!          NPTM = NP
C          call RefreshEmptyCellsV2(PC,PR,PL,IPTID,PCONC,PTWT,
C     *      NPCELL, CONC, CNOLD, CELVOL, SUMWT, 
C     *      IBOUND,
C     *      PNEWC,PNEWR,PNEWL,
C     *      NEWPTS, ! number of initial particles per cell plus 1. 
C     *      IOUTS,  ! unit number of listing file for GWT
C     *      NPMAX, NP, NSCOL,NSROW,NSLAY,ISCOL1,ISROW1,ISLAY1, 
C     *      NCOL,NROW,NLAY)
C          NPTM = NP
C      endif
      
cgzh debug output
c      if(imov.eq.nmov)then
c         WRITE(80,300)
c         WRITE(80,320)
c         DO 1100 KS=1,NSLAY
c         WRITE(80,322) KS
c         DO 1100 IS=1,NSROW
c 1100    WRITE(80,331) (NPCELL(JS,IS,KS),JS=1,NSCOL)
c      end if
cgzh debug end 
      IF(NZERO.GT.NZCRIT) THEN
        WRITE(IOUTS,300)
cgzh use zerocell array to write modpath locations, can be copy/pasted and used as 
c    modpath input file
        WRITE(IOUTS,322) 
        write(*,*)
        DO 100 KS=1,NSLAY
        DO 100 IS=1,NSROW
        DO 100 JS=1,NSCOL
          K=KS+ISLAY1-1
          I=IS+ISROW1-1
          J=JS+ISCOL1-1
  100     IF(ZEROCELL(JS,IS,KS).GT.0) 
     &      WRITE(IOUTS,*) J,I,K,' 0.5   0.5   0.5'
c no regeneration with weighted particles, so stop
        write(*,*)
!  RBW the next line is commented out because now weighted particle regeneration is working.
!        STOP 'FZERO EXCEEDED IN MOVEWT'

cgzh debug
c this is not working yet with weighted particles (regen), so stop
!          write(*,*)
!  The previous write statment can be removed too but is left in for now as a way of 
!  easily seeing that RefreshEmptyCells is being called.          
!          STOP 'NPMAX EXCEEDED IN MOVEWT'
C   
C RBW begin
!          call RefreshEmptyCellsV3(PC,PR,PL,IPTID,PCONC,PTWT,
!     *      NPCELL, CONC, CNOLD, CELVOL, SUMWT, 
!     *      IBOUND,
!     *      PNEWC,PNEWR,PNEWL,
!     *      NEWPTS, ! number of initial particles per cell plus 1. 
!     *      IOUTS,  ! unit number of listing file for GWT
!     *      NPMAX, NP, NSCOL,NSROW,NSLAY,ISCOL1,ISROW1,ISLAY1, 
!     *      NCOL,NROW,NLAY)
     
!          call RefreshEmptyCellsV4(PC,PR,PL,IPTID,PCONC,PTWT,
!     *      NPCELL, CONC, CNOLD, CELVOL, SUMWT, 
!     *      IBOUND,
!     *      PNEWC,PNEWR,PNEWL,
!     *      NEWPTS, ! number of initial particles per cell plus 1. 
!     *      IOUTS,  ! unit number of listing file for GWT
!     *      THCK, VC, VR, VL, POR,    
!     *      NPMAX, NP, NSCOL,NSROW,NSLAY,ISCOL1,ISROW1,ISLAY1, 
!     *      NCOL,NROW,NLAY,
!     *      NZERO, ! number of empty cells.
!     *      TIMV) ! length of transport time step.
!          NPTM = NP
!          IF(NZERO.GT.NZCRIT) THEN
!            STOP 'FZERO EXCEEDED IN MOVEWT'
!          ENDIF
          
          ! Is it necessary to update CAVG, NZERO?
          ! maybe need to update SUMMASS, LIMBO, SUMSGMS, SUMSGWT, SGMSOUT, 
          ! SNKMSOUT, RESIDWT, NSKIP, PCORIG,PRORIG,PLORIG
C RBW end



      END IF
cgzh angle2k bug
C ADD CELLS SKIPPED DUE TO LOW WEIGHT TO NZERO; THIS IS USED AS A FLAG
C TO CALL ROUTINES TO UPDATE CELL CONNECTIONS IN MVOT      
cgzh debug this doesn't work, fix (if NSKIP>0 and NZERO=0 we will call BOOM erroneously)
c      NZERO=NZERO+NSKIP
C
cgzh debug  SKIPPING DEBUG OUTPUT
c      goto 9999
cgzh debug
c         WRITE(IOUTS,3320)
c         DO 101 KS=1,NSLAY
c         DO 101 IS=1,NSROW
c  101    WRITE(IOUTS,331) (NPCELL(JS,IS,KS),JS=1,NSCOL)
c 3320 FORMAT (/,2X,'NUMBER OF PARTICLES PER CELL AT END OF MOVE',/)
c      write(iouts,*)
c         WRITE(IOUTS,3321)
      sumwtte=0.0
      sumptte=0.0
      SUMVOL=0.0
c      WRITE(IOUTS,*) 'sumwt at 411', SUMWT(4,1,1)
c      if(imov.eq.46.or.imov.eq.47.or.imov.eq.48) then
C37 loop over cells to compute cumulative mass, and cumulative volume.
       DO 1201 KS=1,NSLAY
       DO 1201 IS=1,NSROW
       DO 1201 JS=1,NSCOL
c      if(igenpt(js,is,ks).gt.0) 
c      WRITE(IOUTS,333) js,is,ks,SUMWT(JS,IS,KS)
c     &WRITE(IOUTS,333) js,is,ks,SUMWT(JS,IS,KS)
         sumwtte=sumwtte+SUMWT(JS,IS,KS)
         SUMVOL=SUMVOL+SUMWT(JS,IS,KS)
c          WRITE(72,*) 'sumwt move=',imov
 1201  continue
c      end if
  333 FORMAT (5X,3I3,2X,1E12.6)
 3321 FORMAT (/,2X,'SUM OF WEIGHTS OF PARTICLES AT END OF MOVE',/)
c      write(iouts,*)
C38 compute discrepancy in weights? Is this just for debugging?
cgzh debug pt loop
cgzh debug output 
c      WRITE(IOUTS,*) 'After move, SUMVOL= ',SUMVOL
c      if(imov.eq.27) then 
c        write(iouts,*) 'sumwtte at 27=',sumwtte
c      end if
      SMCELL=0.0
      DO 9180 IP=1,NPTM
        IF(PC(IP).EQ.0.0) GO TO 9180
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
C       ---SUM WTS---
        IF(IBOUND(J,I,K).NE.0) THEN
          sumptte=sumptte+PTWT(IP)
C RBW SMCELL is being used for the particle weights instead of the particle mass.
          SMCELL(JS,IS,KS)=SMCELL(JS,IS,KS)+PTWT(IP)
        END IF
 9180 CONTINUE
        write(iouts,*)
c         WRITE(IOUTS,3329)
        DO 1209 KS=1,NSLAY
        DO 1209 IS=1,NSROW
c      WRITE(IOUTS,333) (SMCELL(JS,IS,KS),JS=1,NSCOL)
        DO 1209 JS=1,NSCOL
          badbalance=SUMWT(JS,IS,KS)-SMCELL(JS,IS,KS)
          if(sumwt(js,is,ks).gt.0.0.and.smcell(js,is,ks).gt.0.0) then
c      if(abs(badbalance/SUMWT(JS,IS,KS)).gt.0.00001) then 
c     write(iouts,*) 'Bad weight balance at ',js,is,ks
c      write(iouts,*) 'Sum from cells (SUMWT)=',SUMWT(JS,IS,KS)
c      write(iouts,*) 'Sum from pts (PTWTs)  =',SMCELL(JS,IS,KS)
c      write(iouts,*) 'NPCELL  =',NPCELL(JS,IS,KS)
cgzh debug this stop should probably be active...took down for zinn run
c       stop 'Weight imbalance error'
c      end if
cgzh debug srcfix
            IF(ISRCFIX.GT.0) THEN
              if(igenpt(js,is,ks).eq.1) then
                wtdiff=celvol(js,is,ks)-sumwt(js,is,ks)
                dpcrit=celvol(js,is,ks)*1D-6
cgzh debug output
c      if(js.eq.3.and.is.eq.1.and.ks.eq.1) then
c       write(777,*) imov,wtdiff
c      end if
c       if(abs(wtdiff).gt.dpcrit) then 
c        write(iouts,*) 'Weight not constant in strong source,
c     * diff=',wtdiff,' at js,is,ks=',js,is,ks
c        write(iouts,*) 'sumwt,celvol,dpcrit',
c     &  sumwt(js,is,ks),celvol(js,is,ks),
c     &  dpcrit
c        stop 'Weight not constant in strong source'
c      end if
              END IF
            end if
          end if
 1209   continue
 3329 FORMAT (/,2X,'SUM OF WEIGHTS OF PARTICLES AT END OF MOVE (PTs)',/)
c      write(iouts,*) 'sumwtte,sumptte',sumwtte,sumptte
cgzh debug  skip to 9999 for no debug output here
 9999 CONTINUE
cgzh debug output
c      write(iouts,*) 'move: imov,SUMWT,SUMMASS,npcell 111',
c     * imov,SUMWT(1,1,1),SUMMASS(17,1,1),npcell(js,is,ks)
C39  RELEASE MEMORY
      DEALLOCATE(TEMPMS,TEMPWT,SINKQ,REMVMS,REMVWT,NPNEW,SGPTC,SOURCEQ)
cgzh crit
      DEALLOCATE(CRITWT,CRITMS,NPSUM)
cgzh debug
      DEALLOCATE(SMCELL)
cgzh debug
      DEALLOCATE(QMNWSINK)
cgzh debug
      DEALLOCATE(QSS_SINK)
cgzh srcnew
      DEALLOCATE(SRCMS,SRCWT)
cgzh frac1
      DEALLOCATE(FRAC1,FRAC2)
cgzh srcfix
      DEALLOCATE(WTOUT,LIST,WTIN,LIST2,DEFICIT,DEFMS,AVC,AVR,AVL,LUMP)
      DEALLOCATE(DEFICIT2)
cgzh srcfix2
      DEALLOCATE(SS_WT,SS_MS)
cgzh drycell
      DEALLOCATE (DRYWT,DRYMS,TEMPDRYWT,TEMPDRYMS)
cgzh ccbd
      DEALLOCATE (CHNGMASS)
cgzh zerocell
      DEALLOCATE (ZEROCELL)
cgzh debug median
      DEALLOCATE(IVCOUNT)
C
      RETURN
C     ****************************************************************
C
  290 FORMAT (/,5X,'NUMBER OF CELLS WITH ZERO PARTICLES  =',I6)
  300 FORMAT (/,5X,'*** FZERO EXCEEDED: SIMULATION TERMINATED  ***',/)
  322 FORMAT(/,' FOLLOWING IS A LIST OF CELLS (J,I,K) THAT HAVE ZERO',/,
     & ' PARTICLES.  COPY AND PASTE THIS LIST INTO THE MODPATH',/,
     & ' STARTING LOCATIONS FILE AND BACKTRACK TO DETERMINE SOURCE .',/,
     & ' CELLS.  INCREASE NUMBER OF PARTICLES IN THOSE CELLS.'/)
cgzh debug (orig line)  331 FORMAT (5X,30I4)
  331 FORMAT (5X,30I8)
  670 FORMAT (2X,'NP',7X,'=',I9,' AT START OF MOVE',
     1 10X,'IMOV     =',I13)
  700 FORMAT(/,5X,' ***  WARNING ***',10X,'# Particles = NPMAX -- IMOV='
     1,I4,2X,'PT. NO.=',I8,5X,'REGENERATING PARTICLES',/)
  710 FORMAT (/,5X,' *** WARNING ***',10X,'NPMAX EXCEEDED -- IMOV='
     1,I4,2X,'PT. NO.=',I8,5X,'MUST INCREASE NPMAX')
      END
C     
C     Create new particles in cells that have none.
C     Used with weighted particles.
      Subroutine RefreshEmptyCellsV4(PC,PR,PL,IPTID,PCONC,PTWT,
     *  NPCELL, CONC, CELVOL, SUMWT,
     *  IBOUND,
     *  PNEWC,PNEWR,PNEWL,
     *  NEWPTS, ! number of initial particles per cell plus 1. 
     *  IOUTS,  ! unit number of listing file for GWT
     *  THCK, VC, VR, VL, POR,    
     *  NPMAX, NP, NSCOL,NSROW,NSLAY,ISCOL1,ISROW1,ISLAY1, 
     *  NCOL,NROW,NLAY,
     *  NZERO, ! number of empty cells.
     *  TIMV) ! length of transport time step.
      ! 1. Make a copy of PTWT 
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      DIMENSION
     *  PC(NPMAX),PR(NPMAX),PL(NPMAX), ! PC, PR, and PL store the positions of the particles
     *  IPTID(NPMAX),              ! IPTID indicates the starting location of the particle within the cell 
     *  PCONC(NPMAX),              ! PCONC stores the concentration of each particle. The concentration can be negative.
     *  PTWT(NPMAX),               ! PTWT stores the weight of each particle
     *  NPCELL(NSCOL,NSROW,NSLAY), ! NPCELL stores the number of particles in a cell.
     *  CONC(NSCOL,NSROW,NSLAY),   ! CONC stores the average concentration in a cell
     *  CELVOL(NSCOL,NSROW,NSLAY), ! CELVOL stores the volume of pore space in a cell. It might or might not account for the water table being below the top of the cell.
     *  SUMWT(NSCOL,NSROW,NSLAY),  ! SUMWT stores the sum of the weights of particles in a cell.
     *  IBOUND(NCOL,NROW,NLAY),    ! IBOUND indicates whether a cell is an active cell (>0), inactive cell (=0), or specified head cell (<0)
     *  PNEWC(NEWPTS),PNEWR(NEWPTS),PNEWL(NEWPTS) ! These store the positions of new particles
      DIMENSION THCK(NSCOL,NSROW,NSLAY), ! cell thickness (possibly saturated thickness)
     *  VC(NSCOL+1,NSROW,NSLAY),         ! velocity in the column direction
     *  VR(NSCOL,NSROW+1,NSLAY),         ! velocity in the row direction
     *  VL(NSCOL,NSROW,NSLAY+1),         ! velocity in the layer direction
     *  POR(NSCOL,NSROW,NSLAY)           ! porosity
      DOUBLE PRECISION PTWT, CELVOL   
      double precision SUMWT
      intent (in out) PC,PR,PL,PCONC,PTWT,IPTID,NP,NPCELL,CONC,SUMWT,
     *  NZERO
      intent (in) CELVOL,NSCOL,NSROW,NSLAY,ISCOL1,ISROW1,ISLAY1, 
     *  NCOL,NROW,NLAY,NPMAX,IBOUND,PNEWC,PNEWR,PNEWL,NEWPTS,IOUTS,
     *  THCK,VC,VR,VL,POR
C     Local variables        
      integer J,I,K      ! Column, Row and Layer of full grid respectively
      integer JS, IS, KS ! Column, Row and Layer of transport subgrid respectively
      DOUBLE PRECISION dt, dc, dr, dp 
      DOUBLE PRECISION TopFacePoreArea, FrontFacePoreArea, 
     *  SideFacePoreArea
      DOUBLE PRECISION New_PTWT(NPMAX)
      integer IP ! particle iterator
      integer ICell ! cell iterator
      double precision Flux,Velocity,ParticleWeight
      real FlowFraction
      integer NewNZERO
      type TEmptyCell
        integer JS, IS, KS            ! Column, Row and Layer of the empty cell in subgrid.
        Integer N1_JS, N1_IS, N1_KS   ! Column, Row and Layer of the neighbor 1 cell in subgrid. Neighbor 1 is the cell with the highest flux into the empty cell.
        Integer N2_JS, N2_IS, N2_KS   ! Column, Row and Layer of the neighbor 2 cell in subgrid. Neighbor 2 is the cell with the second highest flux into the empty cell.
        double precision Flux1, Flux2 ! Flux into empty cell from neighboring cells 1 and 2 respectively.
        double precision NewWeight    ! Desired weight to be assigned to new particle in the cell.
        double precision ActualWeight ! Actual weight available to be assigned to a new particle.
        double precision NewConc      ! Ultimately, this is concentration to be assigned to new particle in the cell. In an intermediate step it is the sum of the weights*concentrations of particles from which mass was removed.
      end type TEmptyCell
      type(TEmptyCell) EmptyCells(NZERO) 
      double precision SumOfMass(NSCOL,NSROW,NSLAY)
      parameter HalfCell = 0.499
      ! Quit if too many particles would need to be created.
      if ((NZERO+NP).gt.NPMAX) then
        WRITE(IOUTS,300)
  300   Format(/,5X, 
     *    '*** TOO MANY PARTICLES WOULD NEED TO BE CREATED: ***',/)
        WRITE(IOUTS,301)
  301   Format(/,5X, 
     *    '*** SIMULATION TERMINATED ***',/)
        STOP 'TOO MANY PARTICLES WOULD NEED TO BE CREATED'
      endif
      ! 1. Make a copy of PTWT 
      ! New particle weights will be assigned in New_PTWT and then
      ! assigned back to PTWT when the calculation is completed.
      ! This is needed because a single particle might have its weight adjusted
      ! more than once in this subroutine and the amount by which it should be
      ! adjusted depends on the original particle weight.
      New_PTWT = PTWT

      ! 2. initialize EmptyCells
      do ICell=1,NZERO
        EmptyCells(ICell)%JS = 0
        EmptyCells(ICell)%IS = 0
        EmptyCells(ICell)%KS = 0
        EmptyCells(ICell)%N1_JS = 0
        EmptyCells(ICell)%N1_IS = 0
        EmptyCells(ICell)%N1_KS = 0
        EmptyCells(ICell)%N2_JS = 0
        EmptyCells(ICell)%N2_IS = 0
        EmptyCells(ICell)%N2_KS = 0
        EmptyCells(ICell)%Flux1 = 0
        EmptyCells(ICell)%Flux2 = 0
        EmptyCells(ICell)%ActualWeight = 0;
        EmptyCells(ICell)%NewWeight = 0
        EmptyCells(ICell)%NewConc = 0;
      end do  
      ! 3. Fill EmptyCells with locations and flow rates of cells that lack particles.
      ICell = 1
      LayerLoop1: Do KS=1,NSLAY 
        K = KS + ISLAY1 -1
        RowLoop1: Do IS=1, NSROW
          I=IS+ISROW1-1
          ColumnLoop1: Do JS=1, NSCOL
            J=JS+ISCOL1-1
            IF (IBOUND(J,I,K).NE.0) then
              if (NPCELL(JS,IS,KS).eq.0) then
                EmptyCells(ICell)%JS = JS
                EmptyCells(ICell)%IS = IS
                EmptyCells(ICell)%KS = KS
                dt=dble(THCK(JS,IS,KS))
                dc=dble(CDEL)
                dr=dble(RDEL)
                dp=dble(POR(JS,IS,KS))
                TopFacePoreArea = dc*dr*dp
                FrontFacePoreArea = dc*dt*dp
                SideFacePoreArea = dr*dt*dp
                ! Compute flow through each of the faces with neighboring cells.
                ! Store values for the cells with the top 2 flows in EmptyCells.
                
                ! left face 
                Velocity = dble(VC(JS,IS,KS))
                Flux = SideFacePoreArea*Velocity
                if (Flux.gt.0) then
                  if (Flux.gt.EmptyCells(ICell)%Flux1) then
                    EmptyCells(ICell)%N2_JS = EmptyCells(ICell)%N1_JS
                    EmptyCells(ICell)%N2_IS = EmptyCells(ICell)%N1_IS
                    EmptyCells(ICell)%N2_KS = EmptyCells(ICell)%N1_KS
                    EmptyCells(ICell)%Flux2 = EmptyCells(ICell)%Flux1

                    EmptyCells(ICell)%N1_JS = JS-1
                    EmptyCells(ICell)%N1_IS = IS
                    EmptyCells(ICell)%N1_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  else if (Flux.gt.EmptyCells(ICell)%Flux2) then
                    EmptyCells(ICell)%N2_JS = JS-1
                    EmptyCells(ICell)%N2_IS = IS
                    EmptyCells(ICell)%N2_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  endif
                endif
                ! right face 
                Velocity = dble(VC(JS+1,IS,KS))
                Flux = -SideFacePoreArea*Velocity
                if (Flux.gt.0) then
                  if (Flux.gt.EmptyCells(ICell)%Flux1) then
                    EmptyCells(ICell)%N2_JS = EmptyCells(ICell)%N1_JS
                    EmptyCells(ICell)%N2_IS = EmptyCells(ICell)%N1_IS
                    EmptyCells(ICell)%N2_KS = EmptyCells(ICell)%N1_KS
                    EmptyCells(ICell)%Flux2 = EmptyCells(ICell)%Flux1

                    EmptyCells(ICell)%N1_JS = JS+1
                    EmptyCells(ICell)%N1_IS = IS
                    EmptyCells(ICell)%N1_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  else if (Flux.gt.EmptyCells(ICell)%Flux2) then
                    EmptyCells(ICell)%N2_JS = JS+1
                    EmptyCells(ICell)%N2_IS = IS
                    EmptyCells(ICell)%N2_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  endif
                endif
                ! back face 
                Velocity = dble(VR(JS,IS,KS))
                Flux = FrontFacePoreArea*Velocity
                if (Flux.gt.0) then
                  if (Flux.gt.EmptyCells(ICell)%Flux1) then
                    EmptyCells(ICell)%N2_JS = EmptyCells(ICell)%N1_JS
                    EmptyCells(ICell)%N2_IS = EmptyCells(ICell)%N1_IS
                    EmptyCells(ICell)%N2_KS = EmptyCells(ICell)%N1_KS
                    EmptyCells(ICell)%Flux2 = EmptyCells(ICell)%Flux1

                    EmptyCells(ICell)%N1_JS = JS
                    EmptyCells(ICell)%N1_IS = IS-1
                    EmptyCells(ICell)%N1_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  else if (Flux.gt.EmptyCells(ICell)%Flux2) then
                    EmptyCells(ICell)%N2_JS = JS
                    EmptyCells(ICell)%N2_IS = IS-1
                    EmptyCells(ICell)%N2_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  endif
                endif
                ! front face 
                Velocity = dble(VR(JS,IS+1,KS))
                Flux = -FrontFacePoreArea*Velocity
                if (Flux.gt.0) then
                  if (Flux.gt.EmptyCells(ICell)%Flux1) then
                    EmptyCells(ICell)%N2_JS = EmptyCells(ICell)%N1_JS
                    EmptyCells(ICell)%N2_IS = EmptyCells(ICell)%N1_IS
                    EmptyCells(ICell)%N2_KS = EmptyCells(ICell)%N1_KS
                    EmptyCells(ICell)%Flux2 = EmptyCells(ICell)%Flux1

                    EmptyCells(ICell)%N1_JS = JS
                    EmptyCells(ICell)%N1_IS = IS+1
                    EmptyCells(ICell)%N1_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  else if (Flux.gt.EmptyCells(ICell)%Flux2) then
                    EmptyCells(ICell)%N2_JS = JS
                    EmptyCells(ICell)%N2_IS = IS+1
                    EmptyCells(ICell)%N2_KS = KS
                    EmptyCells(ICell)%Flux1 = Flux
                  endif
                endif
                ! Top face 
                Velocity = dble(VL(JS,IS,KS))
                Flux = TopFacePoreArea*Velocity
                if (Flux.gt.0) then
                  if (Flux.gt.EmptyCells(ICell)%Flux1) then
                    EmptyCells(ICell)%N2_JS = EmptyCells(ICell)%N1_JS
                    EmptyCells(ICell)%N2_IS = EmptyCells(ICell)%N1_IS
                    EmptyCells(ICell)%N2_KS = EmptyCells(ICell)%N1_KS
                    EmptyCells(ICell)%Flux2 = EmptyCells(ICell)%Flux1

                    EmptyCells(ICell)%N1_JS = JS
                    EmptyCells(ICell)%N1_IS = IS
                    EmptyCells(ICell)%N1_KS = KS-1
                    EmptyCells(ICell)%Flux1 = Flux
                  else if (Flux.gt.EmptyCells(ICell)%Flux2) then
                    EmptyCells(ICell)%N2_JS = JS
                    EmptyCells(ICell)%N2_IS = IS
                    EmptyCells(ICell)%N2_KS = KS-1
                    EmptyCells(ICell)%Flux1 = Flux
                  endif
                endif
                ! Bottom face 
                Velocity = dble(VL(JS,IS,KS+1))
                Flux = -TopFacePoreArea*Velocity
                if (Flux.gt.0) then
                  if (Flux.gt.EmptyCells(ICell)%Flux1) then
                    EmptyCells(ICell)%N2_JS = EmptyCells(ICell)%N1_JS
                    EmptyCells(ICell)%N2_IS = EmptyCells(ICell)%N1_IS
                    EmptyCells(ICell)%N2_KS = EmptyCells(ICell)%N1_KS
                    EmptyCells(ICell)%Flux2 = EmptyCells(ICell)%Flux1

                    EmptyCells(ICell)%N1_JS = JS
                    EmptyCells(ICell)%N1_IS = IS
                    EmptyCells(ICell)%N1_KS = KS+1
                    EmptyCells(ICell)%Flux1 = Flux
                  else if (Flux.gt.EmptyCells(ICell)%Flux2) then
                    EmptyCells(ICell)%N2_JS = JS
                    EmptyCells(ICell)%N2_IS = IS
                    EmptyCells(ICell)%N2_KS = KS+1
                    EmptyCells(ICell)%Flux1 = Flux
                  endif
                endif
                
                ICell = ICell + 1
              endif
            endif
          end do ColumnLoop1
        end do RowLoop1
      end do LayerLoop1
      
!     4. Calculate SUMWT, SumOfMass, TotalParticleWeight, TotalParticleMass    
!        TotalParticleWeight and TotalParticleMass aren't used and could be deleted.
!      MAXID = 0       
      SUMWT = 0
      SumOfMass = 0
      TotalParticleWeight = 0
      TotalParticleMass = 0
      ParticleLoop0:  DO IP=1,NP
        IF(PC(IP).EQ.0.0) Cycle ParticleLoop0
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
        IF ((IBOUND(J,I,K).NE.0).and.(IPTID(IP).ne.0)) then
          SUMWT(JS,IS,KS)= SUMWT(JS,IS,KS) + PTWT(IP)
          SumOfMass(JS,IS,KS)= SumOfMass(JS,IS,KS) + PTWT(IP)*PCONC(IP)
          TotalParticleWeight = TotalParticleWeight + PTWT(IP)
          TotalParticleMass = TotalParticleWeight + PTWT(IP)*PCONC(IP)
        endif
      enddo ParticleLoop0
      Continue
      
      ! 5. Compute the weight to be assigned to the new particles in the empty cells.
      ! Limit the weight to no more than half the pore volume of the cell.
      ! Limit the weight to no more than half the weight in the upstream cell.
      ! Would it be better to limit it to no more than half the pore volume of the source cell?
      do ICell=1,NZERO
        if (EmptyCells(ICell)%Flux1.gt.0) then
          EmptyCells(ICell)%NewWeight = TIMV*EmptyCells(ICell)%Flux1
          ! don't let the new weight be greater than half the cell volume.
          JS= EmptyCells(ICell)%JS
          IS= EmptyCells(ICell)%IS
          KS= EmptyCells(ICell)%KS
!          WeightNew = CELVOL(JS,IS,KS)
          ! Limit the weight to no more than half the pore volume of the cell.
          if (EmptyCells(ICell)%NewWeight.gt.CELVOL(JS,IS,KS)/2) then
            EmptyCells(ICell)%NewWeight = CELVOL(JS,IS,KS)/2
          endif
          ! don't let the new weight be greater than half the weight of the particles in neighbor 1.
          JS= EmptyCells(ICell)%N1_JS
          IS= EmptyCells(ICell)%N1_IS
          KS= EmptyCells(ICell)%N1_KS
          if ((JS.NE.0).and.
     *      (EmptyCells(ICell)%NewWeight.gt.SUMWT(JS,IS,KS)/2)) then
               EmptyCells(ICell)%NewWeight = SUMWT(JS,IS,KS)/2
          endif
!          if (WeightNew .gt. EmptyCells(ICell)%NewWeight) then
!            EmptyCells(ICell)%NewWeight = WeightNew
!          endif
          ! don't let the new weight be greater than half the weight of the particles in neighbor 2.
!          JS= EmptyCells(ICell)%N2_JS
!          IS= EmptyCells(ICell)%N2_IS
!          KS= EmptyCells(ICell)%N2_KS
!          if ((JS.NE.0).and.
!     *      (EmptyCells(ICell)%NewWeight.gt.SUMWT(JS,IS,KS)/2)) then
!               EmptyCells(ICell)%NewWeight = SUMWT(JS,IS,KS)/2
!          endif
          ! don't let the new weight be less than zero.
          if (EmptyCells(ICell)%NewWeight.lt.0) then
            EmptyCells(ICell)%NewWeight = 0
          endif
        else
          EmptyCells(ICell)%NewWeight = 0
        endif
      end do  
      
      
!      WeightRemoved = 0
!      TotalMassRemoved = 0
      ! 6. remove some weight from existing particles in cells upstream of empty cells.
      ParticleLoop1:  DO IP=1,NP
        IF(PC(IP).EQ.0.0) Cycle ParticleLoop1
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
        IF ((IBOUND(J,I,K).NE.0).and.(IPTID(IP).ne.0)) then
          !if (PCONC(IP).gt.0) then
            ! ignore negative particle concentrations here?
            do ICell=1,NZERO
              if ((JS.eq.EmptyCells(ICell)%N1_JS)
     *          .and.(IS.eq.EmptyCells(ICell)%N1_IS)
     *          .and.(KS.eq.EmptyCells(ICell)%N1_KS)
     *          .and.(EmptyCells(ICell)%NewWeight.gt.0)) then
                  ! The weight removed from this particle is 
                  ! proportional to the fraction of the weight
                  ! in this cell represented by the particle. 
                  ParticleWeight = (PTWT(IP)/SUMWT(JS,IS,KS)
     *              * EmptyCells(ICell)%NewWeight)
                  EmptyCells(ICell)%ActualWeight = 
     *               EmptyCells(ICell)%ActualWeight + ParticleWeight
!                  WeightRemoved = WeightRemoved + ParticleWeight
!                  TotalMassRemoved = TotalMassRemoved 
!     *               + ParticleWeight*PCONC(IP)
                  New_PTWT(IP) = New_PTWT(IP) - ParticleWeight
                  EmptyCells(ICell)%NewConc = EmptyCells(ICell)%NewConc
     *              + ParticleWeight*PCONC(IP) 
              endif
            end do
          !endif
        end if
      end do ParticleLoop1
      
!      WeightAdded = 0
!      TotalMassAdded = 0
      ! 7. Create new particles just inside the edge of the face with the greatest flux
      ! The exact position will be determined by the relative fluxes on the faces with 
      ! the two greatest fluxes.
      NewNZERO = 0;
!      MAXID = MAXID+1
      do ICell=1,NZERO
        if (EmptyCells(ICell)%ActualWeight.ne.0) then
          JS = EmptyCells(ICell)%JS
          J = JS+ISCOL1-1
          IS = EmptyCells(ICell)%IS
          I = IS+ISROW1-1
          KS = EmptyCells(ICell)%KS
          K = KS+ISLAY1-1
          EmptyCells(ICell)%NewConc = 
     *      EmptyCells(ICell)%NewConc/EmptyCells(ICell)%ActualWeight
          NPCELL(JS,IS,KS) = 1
          NP = NP+1
          IPTID(NP) = 1
!          MAXID = MAXID+1
          PC(NP) = J
          PR(NP) = I
          PL(NP) = K
          New_PTWT(NP) = EmptyCells(ICell)%ActualWeight
          PCONC(NP) = EmptyCells(ICell)%NewConc
!          WeightAdded = WeightAdded + New_PTWT(NP);
!          TotalMassAdded = TotalMassAdded + New_PTWT(NP)*PCONC(NP)
          
          
          if (EmptyCells(ICell)%N1_JS.gt.JS) then
            ! right face
            PC(NP) = PC(NP) + HalfCell
          else if (EmptyCells(ICell)%N1_JS.lt.JS) then
            ! left face
            PC(NP) = PC(NP) - HalfCell
          else if (EmptyCells(ICell)%N1_IS.gt.IS) then
            ! front face
            PR(NP) = PR(NP) + HalfCell
          else if (EmptyCells(ICell)%N1_IS.lt.IS) then
            ! back face
            PR(NP) = PR(NP) - HalfCell
          else if (EmptyCells(ICell)%N1_KS.gt.KS) then
            ! lower face
            PL(NP) = PL(NP) + HalfCell
          else if (EmptyCells(ICell)%N1_KS.lt.KS) then
            ! upper face
            PL(NP) = PL(NP) - HalfCell
          else
            WRITE(IOUTS,400)
  400       Format(/,5X, 
     *         '*** Error when creating particles ***',/)
            STOP  '*** Error when creating particles ***'
          endif
          FlowFraction = EmptyCells(ICell)%Flux2/
     *      (EmptyCells(ICell)%Flux1+EmptyCells(ICell)%Flux2)
          if (FlowFraction.gt.HalfCell) then
            FlowFraction = HalfCell
          endif
          if (EmptyCells(ICell)%N2_JS.ne.0) then
            if (EmptyCells(ICell)%N2_JS.gt.JS) then
              ! right face
              if (EmptyCells(ICell)%N1_JS.eq.JS) then
                 PC(NP) = PC(NP) + FlowFraction
               endif
            else if (EmptyCells(ICell)%N2_JS.lt.JS) then
              ! left face
              if (EmptyCells(ICell)%N1_JS.eq.JS) then
                PC(NP) = PC(NP) - FlowFraction
              endif
            else if (EmptyCells(ICell)%N2_IS.gt.IS) then
              ! front face
              if (EmptyCells(ICell)%N1_IS.eq.IS) then
                PR(NP) = PR(NP) + FlowFraction
              endif
            else if (EmptyCells(ICell)%N2_IS.lt.IS) then
              ! back face
              if (EmptyCells(ICell)%N1_IS.eq.IS) then
                PR(NP) = PR(NP) - FlowFraction
              endif
            else if (EmptyCells(ICell)%N2_KS.gt.KS) then
              ! lower face
              if (EmptyCells(ICell)%N1_KS.eq.KS) then
                PL(NP) = PL(NP) + FlowFraction
              endif
            else if (EmptyCells(ICell)%N2_KS.lt.KS) then 
              ! upper face
              if (EmptyCells(ICell)%N1_KS.eq.KS) then
                PL(NP) = PL(NP) - FlowFraction
              endif
            else
              WRITE(IOUTS,400)
              STOP  '*** Error when creating particles ***'
            endif
          endif
        else
          NewNZERO = NewNZERO + 1
        endif 
      enddo
      NZERO = NewNZERO;
      ! 8. Copy New_PTWT back into PTWT
      PTWT = New_PTWT
      
!     9. Update SUMWT and CONC            
      SUMWT = 0
      SumOfMass = 0
      ParticleLoop3:  DO IP=1,NP
        IF(PC(IP).EQ.0.0) Cycle ParticleLoop3
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
        IF ((IBOUND(J,I,K).NE.0).and.(IPTID(IP).ne.0)) then
          SUMWT(JS,IS,KS)= SUMWT(JS,IS,KS) + PTWT(IP)
          SumOfMass(JS,IS,KS)= SumOfMass(JS,IS,KS) + PTWT(IP)*PCONC(IP)
        endif
      enddo ParticleLoop3
!      RevisedParticleWeight = 0
!      RevisedParticleMass = 0
      ParticleLoop4:  DO IP=1,NP
        IF(PC(IP).EQ.0.0) Cycle ParticleLoop4
        J=PC(IP)+0.5
        JS=J-ISCOL1+1
        I=ABS(PR(IP))+0.5
        IS=I-ISROW1+1
        K=PL(IP)+0.5
        KS=K-ISLAY1+1
        IF ((IBOUND(J,I,K).NE.0).and.(IPTID(IP).ne.0)) then
!          RevisedParticleWeight = RevisedParticleWeight + PTWT(IP)
!          RevisedParticleMass = RevisedParticleMass + PTWT(IP)*PCONC(IP)
          if (SUMWT(JS,IS,KS).gt.0) then
            CONC(JS,IS,KS)= SumOfMass(JS,IS,KS)/SUMWT(JS,IS,KS)
          endif
        endif
      enddo ParticleLoop4
      continue
      
      return
      end;     
