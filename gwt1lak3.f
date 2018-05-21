C     Last change: RBW Dec. 8, 2014
C     Support for weighted particles (MOCWT) added.
      SUBROUTINE GWT1LAK3BD(NLAKES,LKNODE,ILAKE,MXLKND,INLAKUNIT,
     1  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NSOL,
     2  IOUTS,SRCFLO,SRCSOL,SNKFLO,CLAKE,RATECO,FLOB,IGENLK,BUFF,
     *  INBFLX,IBOUND,VC,VR,VL,NLAKESAR)
C                                                                      C
C     *****************************************************************C
C       CALCULATE SOLUTE SINK/SOURCE FLUX TERMS FOR ...                C
C                AQUIFER-LAKE INTERACTION                              C
C     *****************************************************************C
C                                                                      C
C     SPECIFICATIONS:                                                  C
C     -----------------------------------------------------------------C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
C
      DIMENSION ILAKE(5,MXLKND),BUFF(NCOL,NROW,NLAY)
      DIMENSION CLAKE(NLAKESAR,NSOL),RATECO(MXLKND,NSOL)
      DIMENSION SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),IGENLK(NSCOL,NSROW,NSLAY)
      DIMENSION FLOB(MXLKND)
      DIMENSION VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1),IBOUND(NCOL,NROW,NLAY)
      CHARACTER*8 PACKAGE
C     -----------------------------------------------------------------C
C  INITIALIZE LAKE REGENERATION FLAG + BUFF
      IGENLK=0
C  RETURN IF LAK PACKAGE IS OFF                                                                    C
      IF(INLAKUNIT.EQ.0) RETURN
C                                                                      C
      IF(NLAKES.EQ.0) GO TO 600
C                                                                      
C2      THERE ARE LAKES.  
C                            
C---Initialize
      DO 50 L=1,LKNODE
               DO 50 ISOL=1,NSOL
 50               RATECO(L,ISOL)=0.0
C
C---LOOP OVER LAKE NODES
      DO 500 L=1,LKNODE
C                                                                      C
C5---DETERMINE LAKE LOCATION.                                          C
      IL=ILAKE(1,L)
      IR=ILAKE(2,L)
      IC=ILAKE(3,L)
C---IGRID=flag for transport subgrid (0=outside; 1=inside)
      IGRID=1
C
C---Check Transport Subgrid                                            C
      KS=IL-ISLAY1+1
      IS=IR-ISROW1+1
      JS=IC-ISCOL1+1
C6---- CHECK FOR CELLS OUTSIDE of transport subgrid
         IF(KS.LT.1.OR.KS.GT.NSLAY) IGRID=0
         IF(IS.LT.1.OR.IS.GT.NSROW) IGRID=0
         IF(JS.LT.1.OR.JS.GT.NSCOL) IGRID=0
         IF (IGRID.EQ.0) GO TO 500
C                                                                      C
C       SET IGENLK
         IF(BUFF(IC,IR,IL).LT.0.0) IGENLK(JS,IS,KS)=1
c 8/5/03 include this line (orig commented out) 
         IF(BUFF(IC,IR,IL).GT.0.0) IGENLK(JS,IS,KS)=-1
C                                                                      C
C9------FORMULATE SINK AND SOURCE TERMS FOR TRANSPORT
C
C---(SRCSOL array must be re-dimensioned for multiple constituents)
       	  RATE=FLOB(L)
            IF(INBFLX.GT.0) THEN
              IFACE=ILAKE(5,L)
              PACKAGE='    LAKE'
              IPCK=L
              CALL GWT1BFLX5LK(IFACE,RATE,PACKAGE,IPCK,
     *          VC,VR,VL,
     *          NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *          IOUTS,
     *          IC,IR,IL,JS,IS,KS)
            END IF
C
            IF(RATE.LT.0.0) THEN
               SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+RATE
            ELSE
C
C7------DETERMINE LAKE NUMBER.                                         C
               LAKE=ILAKE(4,L)
               DO 115 ISOL=1,NSOL
               SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+RATE
               RATEC=RATE*CLAKE(LAKE,ISOL)
               RATECO(L,ISOL)=RATEC
 115           SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+RATEC
            END IF
C
  500 CONTINUE
  600 CONTINUE
C
C31-----RETURN.                                                        C
      RETURN
      END
C
      SUBROUTINE CILK3BD(LKNODE,MXLKND,ILAKE,KLK,CLKSUM,
     1      NCOL,NROW,NLAY,NLAKES,NSSAR,NTRB,NDV,ISS,
     2      ITRB,IDIV,STROUT,NSLAY,NSROW,NSCOL,NSOL,TIMV,KCNT,
     3      MSUB,MSUB1,LKDONE,SOLPPT,CSWIN,CSWOUT,CLAKE,CGWL,
     4      WTHDRW,CWDRAW,CAUG,CONC,CNOLD,IOUTS,IIMOV,
     5      PRECIP,VOL,VOLOLD,CPPT,CGWIN,CGWOUT,CLKOLD,
     6      CSRUN,CRUNF,CSLAKE,FLOB,RNF,
CMOCWT
     7      SNKMSOUT,SNKFLO,SUMBFMS,SUMBFWT)
C
C     ******************************************************************
C     CALCULATE BUDGET TERMS (EXCEPT STREAM INFLOW) FOR ALL LAKES;
C     UPDATE SOLUTE CONC. IN HEADWATER & ISOLATED LAKES AND LAKE SYSTEMS
C     ******************************************************************
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      DOUBLE PRECISION PRECIP
      DIMENSION ILAKE(5,MXLKND),RNF(NLAKES),MSUB(NLAKES,NLAKES),
     1 MSUB1(NLAKES),STROUT(NSSAR),ITRB(NLAKES,NTRB),
     2 IDIV(NLAKES,NDV),FLOB(MXLKND)
       DIMENSION CLAKE(NLAKES,NSOL),CPPT(NLAKES,NSOL),
     7 VOLOLD(NLAKES),CAUG(NLAKES,NSOL),CLKOLD(NLAKES,NSOL),
     8 WTHDRW(NLAKES),CWDRAW(NLAKES,NSOL)
      DIMENSION PRECIP(NLAKES),CRUNF(NLAKES,NSOL),
     1 SOLPPT(NLAKES,NSOL),
     2 VOL(NLAKES),CSLAKE(NLAKES,NSOL),CSWIN(NLAKES,NSOL),
     3 CSWOUT(NLAKES,NSOL),CSRUN(NLAKES,NSOL),CGWIN(NLAKES,NSOL),
     4 CGWOUT(NLAKES,NSOL),CONC(NSCOL,NSROW,NSLAY),
     5 CNOLD(NSCOL,NSROW,NSLAY),CGWL(NSOL)
      DIMENSION LKDONE(NLAKES),KLK(NLAKES),CLKSUM(NSOL)
CMOCWT
      DIMENSION SNKMSOUT(NSCOL,NSROW,NSLAY),SNKFLO(NSCOL,NSROW,NSLAY),
     1          SUMBFMS(NSCOL,NSROW,NSLAY),SUMBFWT(NSCOL,NSROW,NSLAY)
      INCLUDE 'ptwt.inc'
cgzh debug double snkmsout
      DOUBLE PRECISION SNKMSOUT,SUMBFMS,SUMBFWT
C     ------------------------------------------------------------------
C
      RATIN = 0.
      RATOUT = 0.
C
C-------- TRACK UPDATING OF LAKE CONCENTRATIONS FOR THIS TIME INCREMENT
C---          (KLK > 0 indicates lake has tributary stream inflow)
C---          (LKDONE > 0 indicates updating completed for that lake)
C---          Update VOLOLD(LK) if IMOV>1 (only for transient flow cases)
C------
C2------PROCESS EACH CELL IN THE ILAKE LIST.
C----     Initialize variables and save previous lake concentration
C
      DO 30 LK=1,NLAKES
         KLK(LK)=0
         LKDONE(LK)=0
         IF (IIMOV.GT.1) VOLOLD(LK)=VOL(LK)
         DO 100 ISOL=1,NSOL
            CLKOLD(LK,ISOL)=CLAKE(LK,ISOL)
            CGWIN(LK,ISOL)=0.0
            CGWOUT(LK,ISOL)=0.0
            SOLPPT(LK,ISOL)=0.0
            CSWIN(LK,ISOL)=0.0
 100        CSWOUT(LK,ISOL)=0.0
 30      CONTINUE
C
C1A------IF NO LAKE NODES, KEEP ZERO IN ACCUMULATORS.
      IF (LKNODE.EQ.0) GO TO 1200
C
C2A --- FLAG LAKES TO SKIP BECAUSE THEY HAVE INFLOWING STREAM REACHES. 
      DO 200 LK=1,NLAKES
         DO 190 ITRIB=1,NTRB
            INODE=ITRB(LK,ITRIB)
            IF (INODE.LE.0) GO TO 190
            KLK(LK)=1
            GO TO 200
 190     CONTINUE
 200  CONTINUE
C   --- Check if flagged lake is part of a connected lake set;
C            if so, flag all lakes in that lake set.
      DO 250 IIC=1,KCNT
         JIC=MSUB1(IIC)
         DO 240 LIC=1,JIC
            LK=MSUB(LIC,IIC)
            IF (KLK(LK).GT.0) THEN
               DO 230 LIC2=1,JIC
                  LK2=MSUB(LIC2,IIC)
 230              KLK(LK2)=1
               GO TO 250
            END IF
 240     CONTINUE
 250  CONTINUE
C
C2B --- SUM UP SOLUTE OUTFLOWS INTO OUTFLOWING STREAM REACHES. 
         DO 300 LK=1,NLAKES
            DO 295 IDV=1,NDV
               INODE=IDIV(LK,IDV)
               IF (INODE.LE.0) GO TO 295
               DO 290 ISOL=1,NSOL
 290              CSWOUT(LK,ISOL)=CSWOUT(LK,ISOL)+
     *                 STROUT(INODE)*CLAKE(LK,ISOL)*TIMV
 295        CONTINUE
 300     CONTINUE
C
C   MASTER NODE LOOP -- COMPUTE LAKEBED SEEPAGE TERMS AND BUDGET TERMS
C
        DO 900 L=1,LKNODE
         IL=ILAKE(1,L)
         IR=ILAKE(2,L)
         IC=ILAKE(3,L)
C
C---IGRID=flag for transport subgrid (0=outside; 1=inside)
         IGRID=1
C
C---Check Transport Subgrid                                            C
         KS=IL-ISLAY1+1
         IS=IR-ISROW1+1
         JS=IC-ISCOL1+1
C
C6---- CHECK FOR CELLS OUTSIDE of transport subgrid
         IF(KS.LT.1.OR.KS.GT.NSLAY) IGRID=0
         IF(IS.LT.1.OR.IS.GT.NSROW) IGRID=0
         IF(JS.LT.1.OR.JS.GT.NSCOL) IGRID=0
C                                                                      C
C4------DETERMINE LAKE AND NODAL LAYER,ROW,COLUMN NUMBER.
         LAKE=ILAKE(4,L) 
C
         RATE=FLOB(L)
C
C------ SEE IF FLOW IS INTO AQUIFER OR INTO LAKE.
         IF (RATE) 885,900,890
C
C------ AQUIFER IS DISCHARGING TO LAKE.
 885     RATOUT=RATOUT-RATE
         DO 888 ISOL=1,NSOL
            IF (IGRID.EQ.0) THEN
               CGWL(ISOL)=CLAKE(LAKE,ISOL)
            ELSE
c   *** following line must be re-dimensioned for multiple solutes ***
CMOCWT
             IF(PTWTON.EQ.0) THEN
               CGWL(ISOL)=(CONC(JS,IS,KS)+CNOLD(JS,IS,KS))*0.5
             ELSE
               CGWL(ISOL)=SNKMSOUT(JS,IS,KS)/SNKFLO(JS,IS,KS)
cgzh bflx
c skip this?? not using cgwl below
c               CGWL(ISOL)=(SNKMSOUT(JS,IS,KS)+SUMBFMS(JS,IS,KS))
c     &                   /((SNKFLO(JS,IS,KS)*TIMV)+SUMBFWT(JS,IS,KS))
cgzh debug
c      if(snkmsout(js,is,ks).ne.0.and.js.eq.16.and.is.eq.7) then
c	 write(iouts,*) 'snkmsout,js,is,ks',snkmsout(js,is,ks),js,is,ks
c	endif
c       write(*,*) 'snkms,sumbf',SNKMSOUT(JS,IS,KS),SUMBFMS(JS,IS,KS)
c       if(sumbfms(js,is,ks).gt.0) 
c     * write(*,*) 'snkms,sumbf',SNKMSOUT(JS,IS,KS),SUMBFMS(JS,IS,KS)
             END IF
            END IF
CMOCWT
            IF(PTWTON.EQ.0) THEN
              CGWIN(LAKE,ISOL)=CGWIN(LAKE,ISOL)-RATE*CGWL(ISOL)*TIMV
            ELSE
CMOCWT SNKMSOUT HAS TIMV
              CGWIN(LAKE,ISOL)=CGWIN(LAKE,ISOL)-SNKMSOUT(JS,IS,KS)
            END IF
 888      CONTINUE     
         GO TO 900
C
C------ AQUIFER IS RECHARGED FROM LAKE.
 890    RATIN=RATIN+RATE
        DO 898 ISOL=1,NSOL
           CGWL(ISOL)=CLAKE(LAKE,ISOL)
 898       CGWOUT(LAKE,ISOL)=CGWOUT(LAKE,ISOL)+RATE*CGWL(ISOL)*TIMV
 900  CONTINUE
C
C11------CONSIDER OTHER LAKE SOURCES & SINKS.
C---         COMPUTE PREVIOUS SOLUTE MASS IN LAKES; IN RAINFALL;
C---              IN WITHDRAWALS.
C---         COMPUTE NEW CONCENTRATION IN HEADWATER & ISOLATED LAKES 
C
      DO 1000 LAKE=1,NLAKES
C
C8a------COMPUTE SOLUTE MASS IN RAINFALL RECHARGE FOR EACH LAKE 
         DO 580 ISOL=1,NSOL
 580        SOLPPT(LAKE,ISOL)=PRECIP(LAKE)*CPPT(LAKE,ISOL)*TIMV
C
         IF (WTHDRW(LAKE).GE.0.0) THEN
            DO 1008 ISOL=1,NSOL
 1008          CWDRAW(LAKE,ISOL)=WTHDRW(LAKE)*CLAKE(LAKE,ISOL)*TIMV
         ELSE
            DO 1018 ISOL=1,NSOL
 1018          CWDRAW(LAKE,ISOL)=WTHDRW(LAKE)*CAUG(LAKE,ISOL)*TIMV
         END IF
C
         IF(RNF(LAKE).LT.0.0) RUNF = -RNF(LAKE)
         IF(RNF(LAKE).GE.0.0) RUNF = RNF(LAKE)*PRECIP(LAKE)
         VOLUME=VOLOLD(LAKE)
         IF (ISS.NE.0) VOLUME=VOL(LAKE)
         DO 1028 ISOL=1,NSOL
            CSRUN(LAKE,ISOL) = RUNF*CRUNF(LAKE,ISOL)*TIMV
C---   MASS OF SOLUTE IN LAKE DURING TRANSPORT TIME INCREMENT
 1028       CSLAKE(LAKE,ISOL)=VOLUME*CLAKE(LAKE,ISOL)
C
C---   CHECK FOR TRIBUTARY STREAM INFLOW TO LAKE
         IF (KLK(LAKE).GT.0) GO TO 1000
C---   PRELIMINARY EST. OF NEW CONCENTRATION IN LAKE
         DO 1038 ISOL=1,NSOL
           IF (VOL(LAKE).GT.0.0) THEN
            CLAKE(LAKE,ISOL)=(CSLAKE(LAKE,ISOL)+CSWIN(LAKE,ISOL)
     *             -CSWOUT(LAKE,ISOL)+SOLPPT(LAKE,ISOL)
     *             +CGWIN(LAKE,ISOL)-CGWOUT(LAKE,ISOL)
     *             -CWDRAW(LAKE,ISOL)+CSRUN(LAKE,ISOL))
     *             /VOL(LAKE)
           ELSE
            CLAKE(LAKE,ISOL)=0.0
           END IF
 1038    CONTINUE
         LKDONE(LAKE)=1
C
 1000   CONTINUE
C
C-----ADJUST CONCENTRATIONS OF COALESCENT MULTIPLE-LAKE SYSTEMS
C   --- Check if updated lake is part of a connected lake set;
C            if so, account for mixing of all lakes in that lake set.
C            (No lakes in updated lake set can have tributary inflow.)
      DO 1250 IIC=1,KCNT
         VOLSUM=0.0
         DO 1210 ISOL=1,NSOL
 1210       CLKSUM(ISOL)=0.0
         JIC=MSUB1(IIC)
         LK=MSUB(1,IIC)
         IF (KLK(LK).EQ.0) THEN
            DO 1230 LIC=1,JIC
               LK2=MSUB(LIC,IIC)
               VOLSUM=VOLSUM+VOL(LK2)
               DO 1220 ISOL=1,NSOL
 1220             CLKSUM(ISOL)=CLKSUM(ISOL)+CLAKE(LK2,ISOL)*VOL(LK2)
 1230       CONTINUE
            DO 1240 LIC=1,JIC
               LK2=MSUB(LIC,IIC)
               DO 1235 ISOL=1,NSOL
 1235             CLAKE(LK2,ISOL)=CLKSUM(ISOL)/VOLSUM
               IF (KLK(LK2).GT.0.OR.LKDONE(LK2).EQ.0) 
     1             WRITE(IOUTS,1237) LK2,IIC
 1237 FORMAT (/1X,'*** WARNING: LAKE NUMBER ',I3,' IS PART OF LAKE SET',
     1I3,' BUT FAILS TEST1 ON TRIBUTARY INFLOW')
 1240        CONTINUE
         ELSE
            DO 1245 LIC=1,JIC
               LK2=MSUB(LIC,IIC)
               IF (KLK(LK2).EQ.0.OR.LKDONE(LK2).GT.0)
     1             WRITE(IOUTS,1238) LK2,IIC
 1238 FORMAT (/1X,'*** WARNING: LAKE NUMBER ',I3,' IS PART OF LAKE SET',
     1I3,' BUT FAILS TEST2 ON TRIBUTARY INFLOW')
 1245       CONTINUE
         END IF
 1250 CONTINUE
C
 1200 CONTINUE
C12-----RETURN.
      RETURN
      END
      SUBROUTINE CLAK3BD(KLK,CLKSUM,NLAKES,NSS,NTRB,
     1      ITRB,STRIN,NSOL,TIMV,KCNT,MSUB,MSUB1,LKDONE,
     2      SOLPPT,CSWIN,CSWOUT,CLAKE,CWDRAW,MDFLG,NSTRM,
     3      VOL,CGWIN,CGWOUT,COUT,CSRUN,CSLAKE,IOUTS,NSTRMAR,CNTRIB)
C
C     ******************************************************************
C     CALCULATE STREAM INFLOW BUDGET TERMS FOR ALL LAKES;
C     UPDATE SOLUTE CONC. IN LAKES AND LAKE SYSTEMS WITH TRIB. INFLOW
C     ******************************************************************
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C     ------------------------------------------------------------------
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION MSUB(NLAKES,NLAKES),MSUB1(NLAKES),
     1 ITRB(NLAKES,NTRB),STRIN(NSS)
      DIMENSION CLAKE(NLAKES,NSOL),CWDRAW(NLAKES,NSOL)
      DIMENSION SOLPPT(NLAKES,NSOL),CGWOUT(NLAKES,NSOL),
     1 VOL(NLAKES),CSLAKE(NLAKES,NSOL),CSWIN(NLAKES,NSOL),
     2 CSWOUT(NLAKES,NSOL),CSRUN(NLAKES,NSOL),CGWIN(NLAKES,NSOL)
      DIMENSION LKDONE(NLAKES),KLK(NLAKES),CLKSUM(NSOL)
      DIMENSION COUT(NSTRMAR,NSOL),CNTRIB(NSS,NSOL)
C     ------------------------------------------------------------------
C
C-------- TRACK UPDATING OF LAKE CONCENTRATIONS FOR THIS TIME INCREMENT
C---          (KLK > 0 indicates lake has tributary stream inflow)
C---          (LKDONE > 0 indicates updating completed for that lake)
      IF (LKDONE(MDFLG).GT.0) RETURN
      IF (KLK(MDFLG).EQ.0) WRITE (IOUTS,35) MDFLG
 35   FORMAT (1X,'*** WARNING ***  UNANTICIPATED KLK=0 FOR LAKE',I3)
C
C   --- Check if lake is part of a connected lake set;
C            if so, process all lakes in that lake set.
      NSET=0
      DO 250 IIC=1,KCNT
         IIC2=IIC
         JIC=MSUB1(IIC)
         DO 240 LIC=1,JIC
            LK=MSUB(LIC,IIC)
            IF (LK.EQ.MDFLG) THEN
               NSET=JIC
               GO TO 260
            END IF
 240     CONTINUE
 250  CONTINUE
      NSET=1
 260  CONTINUE
C
C2A --- UPDATE SOLUTE INFLOWS FROM INFLOWING STREAM REACHES. 
C---COMPUTE NEW CONCENTRATION IN LAKES ACCOUNTING FOR TRIBUTARY INFLOW
      DO 400 NLK=1,NSET
         IF (NSET.EQ.1) THEN
            LK=MDFLG
            DO 300 ITRIB=1,NTRB
               INODE=ITRB(LK,ITRIB)
               IF (INODE.LE.0) GO TO 300
               DO 290 ISOL=1,NSOL
 290             CSWIN(LK,ISOL)=CSWIN(LK,ISOL)+
cgzh 5/03/07 bug fix: replace cout (by reach) with cntrib (by segment)
corig     1                          STRIN(INODE)*COUT(INODE,ISOL)*TIMV
     1                          STRIN(INODE)*CNTRIB(INODE,ISOL)*TIMV
 300        CONTINUE
         DO 305 ISOL=1,NSOL
           IF (VOL(LK).GT.0.0) THEN
              CLAKE(LK,ISOL)=(CSLAKE(LK,ISOL)+CSWIN(LK,ISOL)
     *                        -CSWOUT(LK,ISOL)+SOLPPT(LK,ISOL)
     *                        +CGWIN(LK,ISOL)-CGWOUT(LK,ISOL)
     *                        -CWDRAW(LK,ISOL)+CSRUN(LK,ISOL))
     *                        /VOL(LK)
           ELSE
              CLAKE(LK,ISOL)=0.0
           END IF
 305    CONTINUE
        LKDONE(LK)=1
C
        ELSE
           DO 350 LIC=1,NSET
               LK=MSUB(LIC,IIC2)
            DO 340 ITRIB=1,NTRB
               INODE=ITRB(LK,ITRIB)
               IF (INODE.LE.0) GO TO 340
               DO 330 ISOL=1,NSOL
 330             CSWIN(LK,ISOL)=CSWIN(LK,ISOL)+
cgzh 5/03/07 bug fix: replace cout (by reach) with cntrib (by segment)
corig     1                          STRIN(INODE)*COUT(INODE,ISOL)*TIMV
     1                          STRIN(INODE)*CNTRIB(INODE,ISOL)*TIMV
 340        CONTINUE
C
            DO 345 ISOL=1,NSOL
               IF (VOL(LK).GT.0.0) THEN
                  CLAKE(LK,ISOL)=(CSLAKE(LK,ISOL)+CSWIN(LK,ISOL)
     *                    -CSWOUT(LK,ISOL)+SOLPPT(LK,ISOL)
     *                    +CGWIN(LK,ISOL)-CGWOUT(LK,ISOL)
     *                    -CWDRAW(LK,ISOL)+CSRUN(LK,ISOL))
     *                    /VOL(LK)
               ELSE
                  CLAKE(LK,ISOL)=0.0
               END IF
 345        CONTINUE
            LKDONE(LK)=1
C
 350       CONTINUE
        END IF
 400  CONTINUE
C
C-----ADJUST CONCENTRATIONS OF COALESCENT MULTIPLE-LAKE SYSTEMS
C   --- Check if updated lake is part of a connected lake set;
C            if so, account for mixing of all lakes in that lake set.
      IF (NSET.EQ.1) GO TO 1260
      IIC=IIC2
      VOLSUM=0.0
      DO 1210 ISOL=1,NSOL
 1210    CLKSUM(ISOL)=0.0
      JIC=MSUB1(IIC)
      IF (JIC.NE.NSET) WRITE (IOUTS,1215) NSET
 1215 FORMAT (/1X,'*** WARNING ***  NSET.NE.JIC, NSET =',I3/)
      DO 1230 LIC=1,JIC
         LK2=MSUB(LIC,IIC)
         VOLSUM=VOLSUM+VOL(LK2)
         DO 1220 ISOL=1,NSOL
 1220       CLKSUM(ISOL)=CLKSUM(ISOL)+CLAKE(LK2,ISOL)*VOL(LK2)
 1230 CONTINUE
      DO 1240 LIC=1,JIC
         LK2=MSUB(LIC,IIC)
         DO 1235 ISOL=1,NSOL
 1235       CLAKE(LK2,ISOL)=CLKSUM(ISOL)/VOLSUM
 1240 CONTINUE
 1260 CONTINUE
C
C12-----RETURN.
      RETURN
      END
      SUBROUTINE CLAK3AD(NLAKES,LKNODE,ILAKE,MXLKND,
     1  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NSOL,LKDONE,
     2  IOUTS,SRCFLO,SRCSOL,SNKFLO,CLAKE,RATECO,FLOB,CLKOLD,
     3  CNOLD,CONC,LWRT,KKPER,NPER,KKSTP,NSTP,IMOV,NMOV,NPNTCL,
     4  VOL,SOLPPT,CSWIN,CSWOUT,CWDRAW,CSRUN,CGWIN,CGWOUT,CSLAKE,
     5  ALKIN,ALKOUT,
CMOCWT
     6  SNKMSOUT,SUMBFMS,INBFLX)
C                                                                      C
C     *****************************************************************C
C       CALCULATE SOLUTE SINK/SOURCE FLUX TERMS FOR ...                C
C                AQUIFER-LAKE INTERACTION                              C
C       Update mass-balance (budget) terms for lake solute             C
C       Write output for lakes                                         C
C     *****************************************************************C
C                                                                      C
C     SPECIFICATIONS:                                                  C
C     -----------------------------------------------------------------C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
      DIMENSION ILAKE(5,MXLKND),CWDRAW(NLAKES,NSOL)
      DIMENSION SOLPPT(NLAKES,NSOL),CNOLD(NSCOL,NSROW,NSLAY),
     1 VOL(NLAKES),CSLAKE(NLAKES,NSOL),CSWIN(NLAKES,NSOL),
     2 CSWOUT(NLAKES,NSOL),CSRUN(NLAKES,NSOL),CGWIN(NLAKES,NSOL),
     3 CGWOUT(NLAKES,NSOL),CONC(NSCOL,NSROW,NSLAY)
      DIMENSION CLAKE(NLAKES,NSOL),RATECO(MXLKND,NSOL)
      DIMENSION SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY)
      DIMENSION FLOB(MXLKND),CLKOLD(NLAKES,NSOL),LKDONE(NLAKES),
     1  ALKIN(NSOL),ALKOUT(NSOL)
CMOCWT
      DIMENSION SNKMSOUT(NSCOL,NSROW,NSLAY),SUMBFMS(NSCOL,NSROW,NSLAY)
      INCLUDE 'ptwt.inc'
cgzh debug double snkmsout
      DOUBLE PRECISION SNKMSOUT,SUMBFMS
C     -----------------------------------------------------------------C
C                                                                      C
      IF(NLAKES.EQ.0) GO TO 600
C                                                                      
C2      THERE ARE LAKES.  
C                            
C-------- CHECK UPDATING OF LAKE CONCENTRATIONS FOR THIS TIME INCREMENT
C---          (LKDONE > 0 indicates updating completed for that lake)
      DO 30 LK=1,NLAKES
         IF (LKDONE(LK).LE.0) WRITE(IOUTS,25) LK
 25    FORMAT(/1X,'*** WARNING *** LAKE NUMBER',I3,' WAS NOT UPDATED')
 30   CONTINUE
C
      DO 50 ISOL=1,NSOL
         ALKIN(ISOL)=0.0
 50      ALKOUT(ISOL)=0.0
C
C---LOOP OVER LAKE NODES
      DO 500 L=1,LKNODE
C                                                                      C
C5---DETERMINE LAKE LOCATION.                                          C
      IL=ILAKE(1,L)
      IR=ILAKE(2,L)
      IC=ILAKE(3,L)
C---IGRID=flag for transport subgrid (0=outside; 1=inside)
      IGRID=1
C
C---Check Transport Subgrid                                            C
      KS=IL-ISLAY1+1
      IS=IR-ISROW1+1
      JS=IC-ISCOL1+1
C6---- CHECK FOR CELLS OUTSIDE of transport subgrid
         IF(KS.LT.1.OR.KS.GT.NSLAY) IGRID=0
         IF(IS.LT.1.OR.IS.GT.NSROW) IGRID=0
         IF(JS.LT.1.OR.JS.GT.NSCOL) IGRID=0
         IF (IGRID.EQ.0) GO TO 500
C                                                                      C
C7------DETERMINE LAKE NUMBER.                                         C
      LAKE=ILAKE(4,L)
C
C9------FORMULATE SOLUTE SOURCE TERMS FOR NEXT TRANSPORT TIME INCREMENT
C    ---         and SUM BUDGET TERMS
c   *** following code must be re-dimensioned for multiple solutes ***
       	    RATE=FLOB(L)
         IF(IGRID.GT.0) THEN
            IF(RATE.LT.0.0) THEN
C---           for gaining lake cell (aquifer discharge):
               DO 112 ISOL=1,NSOL
CMOCWT
                 IF(PTWTON.EQ.0) THEN
                  SQ=RATE*(CNOLD(JS,IS,KS)+CONC(JS,IS,KS))*0.5
                 ELSE
CMOCWT
                  SQ=(RATE/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
cgzh bflx
C IF MASS ENTERED THE NO FLOW CELL NEXT TO THE LAKE, ADD THAT FOR MB
c                  IF(INBFLX.GT.0) THEN
c                    IF(SUMBFMS(JS,IS,KS).LT.0.0) THEN
c                      SQ=SQ+SUMBFMS(JS,IS,KS)
c	              END IF
c                  END IF
                 END IF
 112              ALKOUT(ISOL)=ALKOUT(ISOL)+SQ
            ELSE
C---           for losing lake cell (aquifer recharge):
               DO 115 ISOL=1,NSOL
                  SQ=RATECO(L,ISOL)
                  ALKIN(ISOL)=ALKIN(ISOL)+SQ
                  SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)-SQ
                  RATEC=RATE*CLAKE(LAKE,ISOL)
                  RATECO(L,ISOL)=RATEC
 115              SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+RATEC
            END IF
         END IF
C
  500 CONTINUE
  600 CONTINUE
C
C30-----PRINT LAKE FLUID AND SOLUTE BUDGET DATA FOR EACH LAKE
C         CHECK FLAGS FOR OUTPUT
      IPRNT=0
      IF(LWRT.GT.0) THEN
        IF (KKPER.EQ.NPER.AND.KKSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) GO TO 625
        GO TO 800
      END IF
 625  IF (KKPER.EQ.NPER.AND.KKSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-2.AND.KKSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-1.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.GE.1) THEN
         IF(MOD(IMOV,NPNTCL).EQ.0) IPRNT=1
      ENDIF
C SKIP IF NO OUTPUT
      IF (IPRNT.EQ.0) GO TO 800
      WRITE(IOUTS,650)
 650  FORMAT(//1X,'SOLUTE BUDGETS FOR LAKES FOR THIS TIME INCREMENT:',/'
     1 Lake   Lake   Solute  Concen-     Ppt.    Stream    Stream  Withd
     2rawal  Runoff      GW        GW    Solute Mass'/'  No.  Volume    
     3No.   tration   Mass In   Mass In  Mass Out  Net Mass   Mass In   
     4Mass In  Mass Out   in Lake'/)
      DO 690 LK=1,NLAKES
         DO 680 ISOL=1,NSOL
            IF (ISOL.EQ.1) WRITE(IOUTS,675) LK,VOL(LK),ISOL,
     1         CLAKE(LK,ISOL),SOLPPT(LK,ISOL),CSWIN(LK,ISOL),
     2         CSWOUT(LK,ISOL),CWDRAW(LK,ISOL),CSRUN(LK,ISOL),
     3         CGWIN(LK,ISOL),CGWOUT(LK,ISOL),CSLAKE(LK,ISOL)
            IF (ISOL.GT.1) WRITE(IOUTS,676) ISOL,
     1         CLAKE(LK,ISOL),SOLPPT(LK,ISOL),CSWIN(LK,ISOL),
     2         CSWOUT(LK,ISOL),CWDRAW(LK,ISOL),CSRUN(LK,ISOL),
     3         CGWIN(LK,ISOL),CGWOUT(LK,ISOL),CSLAKE(LK,ISOL)
 675        FORMAT(1X,I3,1PE10.2,3X,I2,1X,9E10.2)
 676        FORMAT(16X,I3,1X,1P9E10.2)
 680     CONTINUE
 690  CONTINUE
      WRITE(IOUTS,677)
 677  FORMAT(//)
 800  CONTINUE
C
C31-----RETURN.                                                        C
      RETURN
      END
