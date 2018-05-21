      SUBROUTINE GWT1SFR2BD(NSTRM,STRM,ISTRM,INLAKUNIT,SEG,
     1  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NSOL,ISTCB1,ICBCFL,
     2  IOUTS,NSS,NSEGDIM,SGOTFLW,IDIVAR,IPTFLG,IOTSG,
     3  SRCFLO,SRCSOL,SNKFLO,CONC,CONCPPT,
     4  ISEG,COUT,CONCQ,CONCRUN,CNTRIB,SOLIN,CQIN,CGW,CLAKE,NLAKES,
     5  NLAKESAR)
C                                                                      C
C     *****************************************************************C
C     ROUTE SOLUTE DOWNSTREAM AND ...                                  C
C               ACCOUNT FOR SOLUTE SINK/SOURCES AT STREAMS             C
C     *****************************************************************C
C                                                                      C
C     SPECIFICATIONS:                                                  C
C     -----------------------------------------------------------------C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
      DIMENSION STRM(24,NSTRM),ISTRM(5,NSTRM),ISEG(4,NSEGDIM),
     2          SGOTFLW(NSS),IDIVAR(2,NSS),IOTSG(NSS),SEG(26,NSS)
      DIMENSION COUT(NSTRM,NSOL),CONCQ(NSEGDIM,NSOL),SOLIN(NSOL),
     1          CONCRUN(NSS,NSOL),CNTRIB(NSS,NSOL),CQIN(NSOL),
     2          CGW(NSOL),CLAKE(NLAKESAR,NSOL),CONCPPT(NSEGDIM,NSOL)
      DIMENSION SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),CONC(NSCOL,NSROW,NSLAY)
C     -----------------------------------------------------------------C
C                                                                      C
      IF(NSTRM.EQ.0) GO TO 600
cgzh debug
      write(*,*) concq
C                                                                      
C2      THERE ARE REACHES.  
C                                                                      C
C4------IF THERE ARE STREAMS, SET INPUTS OF FLOW AND SOLUTE.           C
   10 DO 500 L=1,NSTRM
      LL=L-1
C                                                                      C
C5---DETERMINE REACH LOCATION.                                         C
      IL=ISTRM(1,L)
      IR=ISTRM(2,L)
      IC=ISTRM(3,L)
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
C                                                                      C
C7------DETERMINE SEGMENT AND REACH NUMBER.                            C
      ISTSG=ISTRM(4,L)
      NREACH=ISTRM(5,L)
C
C-----Calculate streamflow conc. into each reach
      IF (NREACH.EQ.1) THEN
C-- STORE CONC. IN OUTFLOW FROM PREVIOUS SEGMENT
         IF (ISTSG.GT.1) THEN
             IFLG=ISTRM(4,LL)
             DO 15 ISOL=1,NSOL
   15        CNTRIB(IFLG,ISOL)=COUT(LL,ISOL)
         END IF
C-- DETERMINE INFLOW CONC. FOR 3 SEGMENT TYPES      
         IF (ISEG(3,ISTSG).EQ.5) THEN
             DO 25 ISOL=1,NSOL
   25        CQIN(ISOL)=CONCQ(ISTSG,ISOL)
C
C--          for type-6 (diversionary) stream segment:
         ELSE IF (ISEG(3,ISTSG).EQ.6) THEN
             NDFLG=IDIVAR(1,ISTSG)
             IF (NDFLG.GT.ISTSG) THEN
               WRITE(IOUTS,*) '***ERROR***  STREAM SEGMENT NUMBER ',
     1 NDFLG,' OUT OF ORDER; PROCESSING STOPPED'
               STOP
             END IF
             IF (NDFLG.LT.0) THEN
               IF (INLAKUNIT.LE.0) THEN
                  WRITE (IOUTS,130) ISTSG
  130 FORMAT(1X,'*** ERROR *** NEG. FLAG INDICATES ',
     *       'LAKE DIVERSION FOR STREAM SEGMENT ',I3,
     *       ' BUT LAKE PACKAGE INACTIVE; PROCESSING STOPPED')
                  STOP
               END IF
               MDFLG=IABS(NDFLG)
               DO 32 ISOL=1,NSOL
   32          CQIN(ISOL)=CLAKE(MDFLG,ISOL)
             ELSE
               DO 35 ISOL=1,NSOL
   35          CQIN(ISOL)=CNTRIB(NDFLG,ISOL)
             END IF
C
C--          for type-7 (tributary) stream segment:
         ELSE IF (ISEG(3,ISTSG).EQ.7) THEN
             ITRIB=1
             FLOWIN=0.0
             DO 45 ISOL=1,NSOL
   45        SOLIN(ISOL)=0.0
            DO WHILE (ITRIB.LE.NSS)
                IF(ISTSG.EQ.IOTSG(ITRIB)) THEN
                  TRBFLW=SGOTFLW(ITRIB)
                  FLOWIN=FLOWIN+TRBFLW
                  DO 55 ISOL=1,NSOL
 55                  SOLIN(ISOL)=SOLIN(ISOL)+
     1                    SGOTFLW(ITRIB)*CNTRIB(ITRIB,ISOL)
                END IF
                ITRIB=ITRIB+1      
            END DO
            FLOWIN=FLOWIN+SEG(2,ISTSG)
            IF (SEG(2,ISTSG).GT.0.0) THEN
               DO 57 ISOL=1,NSOL
 57               SOLIN(ISOL)=SOLIN(ISOL)+
     1                 SEG(2,ISTSG)*CONCQ(ISTSG,ISOL)
            END IF
               IF (FLOWIN.LT.0) THEN
                  FLOWIN=0.0
                  WRITE (IOUTS,4) ISTSG
    4           FORMAT (//2X,'*** WARNING *** FLOW INTO TRIBUTARY ',
     1              'STREAM SEGMENT No. ',I3,' WAS NEGATIVE; ',
     2              'FLOWIN RE-SET = 0.0'/)
               END IF
cgzh 4/20/07 avoid divide by zero if no flow into reach
             DO 65 ISOL=1,NSOL
              IF(FLOWIN.GT.0.0) THEN
                CQIN(ISOL)=SOLIN(ISOL)/FLOWIN
              ELSE
                CQIN(ISOL)=0.0
              ENDIF
   65 CONTINUE
         ELSE
C---     if iseg(3).ne.5 or 6 or 7:
             WRITE(IOUTS,*) '*** WARNING *** ISEG(3)<5 or ISEG(3)>7' 
             write (iouts,*)  ISEG(3,ISTSG),L,ISTSG,NREACH
             stop
             DO 75 ISOL=1,NSOL
   75        CQIN(ISOL)=COUT(LL,ISOL)
         END IF
      ELSE 
C---     if nreach.ne.1:
         DO 85 ISOL=1,NSOL
   85    CQIN(ISOL)=COUT(LL,ISOL)
      END IF
C
C8------Calculate conc. out of stream cell; route solute downstream
      IF (STRM(11,L).GE.0.0) THEN
C---        for losing stream:
         QFCT=STRM(10,L)+STRM(12,L)+STRM(13,L)+STRM(14,L)
         IF (QFCT.GT.0.0) THEN
           DO 95 ISOL=1,NSOL
   95       COUT(L,ISOL)=(STRM(10,L)*CQIN(ISOL)+STRM(12,L)*
     1           CONCRUN(ISTSG,ISOL)+STRM(14,L)*CONCPPT(ISTSG,ISOL))/
     2           QFCT
         ELSE
           DO 96 ISOL=1,NSOL
   96       COUT(L,ISOL)=0.0
         END IF
      ELSE 
C---        for gaining stream:
         QOUT=STRM(10,L)-STRM(11,L)+STRM(12,L)-STRM(13,L)+STRM(14,L)
         IF (QOUT.GT.0.0) THEN
           DO 105 ISOL=1,NSOL
           IF (IGRID.EQ.0) THEN
             CGW(ISOL)=CQIN(ISOL)
           ELSE
c   *** following line must be re-dimensioned for multiple solutes ***
             CGW(ISOL)=CONC(JS,IS,KS)
           END IF
           COUT(L,ISOL)=(STRM(10,L)*CQIN(ISOL)-STRM(11,L)*CGW(ISOL)+
     1        STRM(12,L)*CONCRUN(ISTSG,ISOL)+
     2        STRM(14,L)*CONCPPT(ISTSG,ISOL))/QOUT
  105      CONTINUE
         ELSE
           DO 106 ISOL=1,NSOL
  106       COUT(L,ISOL)=0.0
         END IF
      END IF
C
C9------FORMULATE SINK AND SOURCE TERMS FOR TRANSPORT
c   *** following code must be re-dimensioned for multiple solutes ***
         IF(IGRID.GT.0) THEN
            RATE=STRM(11,L)
            IF(RATE.LT.0.0) THEN
               SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+RATE
            ELSE
               DO 115 ISOL=1,NSOL
               SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+RATE
               RATEC=RATE*COUT(L,ISOL)
 115           SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+RATEC
            END IF
         END IF
  500 CONTINUE
  600 CONTINUE
C     
C31-----RETURN.                                                        C
      RETURN
      END
      SUBROUTINE CSTR3BD(NSTRM,STRM,ISTRM,SEG,
     1  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NSOL,ISTCB1,ICBCFL,
     2  IOUTS,NSS,SGOTFLW,IDIVAR,IPTFLG,IOTSG,
     3  SRCFLO,SRCSOL,SNKFLO,CSTIN,CSTOUT,COUTOLD,
     4  ISEG,COUT,CONCQ,CONCRUN,CNTRIB,SOLIN,CQIN,CGW,
     5  KLK,CLKSUM,NLAKES,NTRB,INLAKUNIT,CONC,CNOLD,CONCPPT,
     1  ITRB,STRIN,TIMV,KCNT,MSUB,MSUB1,LKDONE,
     2  SOLPPT,CSWIN,CSWOUT,CLAKE,CWDRAW,
     3  VOL,CGWIN,CGWOUT,CSRUN,CSLAKE,
     4  NPNTCL,KKPER,NPER,KKSTP,NSTP,IIMOV,NMOV,NSEGDIM,
     5  NLAKESAR,NSTRMAR)
C                                                                      C
C     *****************************************************************C
C     ROUTE SOLUTE DOWNSTREAM AND UPDATE SOLUTE SINK/SOURCES TERMS     C
C            AT STREAMS FOR NEXT TIME INCREMENT...                     C
C         AND CALCULATE SOLUTE BUDGETS FOR STREAM-AQUIFER INTERACTION  C
C            BASED ON S-W CONCENTRATIONS AT START OF TIME INCREMENT.   C
C         Update lake concentrations if lake package active.           C
C     *****************************************************************C
C                                                                      C
C     SPECIFICATIONS:                                                  C
C     -----------------------------------------------------------------C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
      DIMENSION STRM(24,NSTRM),ISTRM(5,NSTRM),ISEG(4,NSEGDIM),
     2          SGOTFLW(NSS),IDIVAR(2,NSS),IOTSG(NSS),SEG(26,NSS)
      DIMENSION COUT(NSTRM,NSOL),CONCQ(NSEGDIM,NSOL),SOLIN(NSOL),
     1          CONCRUN(NSEGDIM,NSOL),CNTRIB(NSS,NSOL),CQIN(NSOL),
     2          CGW(NSOL),CONCPPT(NSEGDIM,NSOL)
      DIMENSION SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),
     *  CONC(NSCOL,NSROW,NSLAY),CNOLD(NSCOL,NSROW,NSLAY)
      DIMENSION COUTOLD(NSOL),CSTIN(NSOL),CSTOUT(NSOL)
      DIMENSION MSUB(NLAKESAR,NLAKESAR),MSUB1(NLAKESAR),
     1 ITRB(NLAKESAR,NTRB),STRIN(NSS)
      DIMENSION CLAKE(NLAKESAR,NSOL),CWDRAW(NLAKESAR,NSOL)
      DIMENSION SOLPPT(NLAKESAR,NSOL),CGWOUT(NLAKESAR,NSOL),
     1 VOL(NLAKESAR),CSLAKE(NLAKESAR,NSOL),CSWIN(NLAKESAR,NSOL),
     2 CSWOUT(NLAKESAR,NSOL),CSRUN(NLAKESAR,NSOL),CGWIN(NLAKESAR,NSOL)
      DIMENSION LKDONE(NLAKESAR),KLK(NLAKESAR),CLKSUM(NSOL)
C     -----------------------------------------------------------------C
C                                                                      C
      IF(NSTRM.EQ.0) GO TO 600
      DO 5 ISOL=1,NSOL
      COUTOLD(ISOL)=0.0
      CSTIN(ISOL)=0.0
    5 CSTOUT(ISOL)=0.0
C                                                                      
C2      THERE ARE REACHES.  
C                                                                      C
C4------IF THERE ARE STREAMS, SET INPUTS OF FLOW AND SOLUTE.           C
   10 DO 500 L=1,NSTRM
      LL=L-1
C                                                                      C
C5---DETERMINE REACH LOCATION.                                         C
      IL=ISTRM(1,L)
      IR=ISTRM(2,L)
      IC=ISTRM(3,L)
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
C                                                                      C
C7------DETERMINE SEGMENT AND REACH NUMBER.                            C
      ISTSG=ISTRM(4,L)
      NREACH=ISTRM(5,L)
C
C-----Calculate streamflow conc. into each reach
      IF (NREACH.EQ.1) THEN
C-- STORE CONC. IN OUTFLOW FROM PREVIOUS SEGMENT
         IF (ISTSG.GT.1) THEN
             IFLG=ISTRM(4,LL)
             DO 15 ISOL=1,NSOL
   15        CNTRIB(IFLG,ISOL)=COUT(LL,ISOL)
         END IF
C-- DETERMINE INFLOW CONC. FOR 3 SEGMENT TYPES      
         IF (ISEG(3,ISTSG).EQ.5) THEN
             DO 25 ISOL=1,NSOL
   25        CQIN(ISOL)=CONCQ(ISTSG,ISOL)
C
C--          for type-6 (diversionary) stream segment:
         ELSE IF (ISEG(3,ISTSG).EQ.6) THEN
             NDFLG=IDIVAR(1,ISTSG)
             IF (NDFLG.LT.0) THEN
               MDFLG=IABS(NDFLG)
C----    Update lake concentration
               IF (INLAKUNIT.EQ.0) THEN
                  WRITE (IOUTS,630) ISTSG,NDFLG
 630             FORMAT(//1X,'*** ERROR ***  SOURCE FOR STREAM SEGMENT',
     1            I4,' IS LAKE',I4/15X,'LAKE PACKAGE INACTIVE; ',
     2            'EXECUTION STOPPING'/)
                  STOP
               END IF
               CALL CLAK3BD(KLK,CLKSUM,NLAKES,NSS,NTRB,
     1      ITRB,STRIN,NSOL,TIMV,KCNT,MSUB,MSUB1,LKDONE,
     2      SOLPPT,CSWIN,CSWOUT,CLAKE,CWDRAW,MDFLG,NSTRM,
     3      VOL,CGWIN,CGWOUT,COUT,CSRUN,CSLAKE,IOUTS,NSTRMAR,CNTRIB)
               DO 32 ISOL=1,NSOL
 32               CQIN(ISOL)=CLAKE(MDFLG,ISOL)
             ELSE
               DO 35 ISOL=1,NSOL
   35          CQIN(ISOL)=CNTRIB(NDFLG,ISOL)
             END IF
C
C--          for type-7 (tributary) stream segment:
         ELSE IF (ISEG(3,ISTSG).EQ.7) THEN
             ITRIB=1
             FLOWIN=0.0
             DO 45 ISOL=1,NSOL
   45        SOLIN(ISOL)=0.0
            DO WHILE (ITRIB.LE.NSS)
                IF(ISTSG.EQ.IOTSG(ITRIB)) THEN
                  TRBFLW=SGOTFLW(ITRIB)
                  FLOWIN=FLOWIN+TRBFLW
                  DO 55 ISOL=1,NSOL
 55                  SOLIN(ISOL)=SOLIN(ISOL)+
     1                    SGOTFLW(ITRIB)*CNTRIB(ITRIB,ISOL)
                END IF
                ITRIB=ITRIB+1      
            END DO
            FLOWIN=FLOWIN+SEG(2,ISTSG)
            IF (SEG(2,ISTSG).GT.0.0) THEN
               DO 57 ISOL=1,NSOL
 57               SOLIN(ISOL)=SOLIN(ISOL)+
     1                 SEG(2,ISTSG)*CONCQ(ISTSG,ISOL)
            END IF
               IF (FLOWIN.LT.0) THEN
                  FLOWIN=0.0
               END IF
cgzh 4/20/07 avoid divide by zero if no flow into reach
             DO 65 ISOL=1,NSOL
              IF(FLOWIN.GT.0.0) THEN
                CQIN(ISOL)=SOLIN(ISOL)/FLOWIN
              ELSE
                CQIN(ISOL)=0.0
              ENDIF
   65 CONTINUE
         ELSE
             DO 75 ISOL=1,NSOL
   75        CQIN(ISOL)=COUT(LL,ISOL)
         END IF
      ELSE 
         DO 85 ISOL=1,NSOL
   85    CQIN(ISOL)=COUT(LL,ISOL)
      END IF
C
C8------Calculate conc. out of stream cell; route solute downstream
      IF (STRM(11,L).GE.0.0) THEN
C---        for losing stream:
         QFCT=STRM(10,L)+STRM(12,L)-STRM(13,L)+STRM(14,L)
         IF (QFCT.GT.0.0) THEN
          DO 95 ISOL=1,NSOL
            COUTOLD(ISOL)=COUT(L,ISOL)
   95       COUT(L,ISOL)=(STRM(10,L)*CQIN(ISOL)+STRM(12,L)*
     1           CONCRUN(ISTSG,ISOL)+STRM(14,L)*CONCPPT(ISTSG,ISOL))/
     2           QFCT
         ELSE
           DO 96 ISOL=1,NSOL
            COUTOLD(ISOL)=COUT(L,ISOL)
   96       COUT(L,ISOL)=0.0
         END IF
	      ELSE 
C---        for gaining stream:
         QOUT=STRM(10,L)-STRM(11,L)+STRM(12,L)-STRM(13,L)+STRM(14,L)
         IF (QOUT.GT.0.0) THEN
          DO 105 ISOL=1,NSOL
           IF (IGRID.EQ.0) THEN
              CGW(ISOL)=CQIN(ISOL)
           ELSE
c   *** following line must be re-dimensioned for multiple solutes ***
              CGW(ISOL)=CONC(JS,IS,KS)
           END IF
           COUTOLD(ISOL)=COUT(L,ISOL)
           COUT(L,ISOL)=(STRM(10,L)*CQIN(ISOL)-STRM(11,L)*CGW(ISOL)+
     1        STRM(12,L)*CONCRUN(ISTSG,ISOL)+
     2        STRM(14,L)*CONCPPT(ISTSG,ISOL))/QOUT
  105     CONTINUE
         ELSE
           DO 106 ISOL=1,NSOL
            COUTOLD(ISOL)=COUT(L,ISOL)
  106       COUT(L,ISOL)=0.0
         END IF
      END IF
C
C9------FORMULATE SINK AND SOURCE TERMS FOR NEXT TRANSPORT TIME INCREMENT
C    ---         and SUM BUDGET TERMS
c   *** following code must be re-dimensioned for multiple solutes ***
         IF(IGRID.GT.0) THEN
            RATE=STRM(11,L)
            IF(RATE.LT.0.0) THEN
C---           for gaining stream:
               DO 112 ISOL=1,NSOL
  112          CSTOUT(ISOL)=CSTOUT(ISOL)+
     *                      RATE*(CNOLD(JS,IS,KS)+CONC(JS,IS,KS))*0.5
            ELSE
C---           for losing stream:
               DO 115 ISOL=1,NSOL
               RATECOLD=RATE*COUTOLD(ISOL)
               CSTIN(ISOL)=CSTIN(ISOL)+RATECOLD
               SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)-RATECOLD
               RATEC=RATE*COUT(L,ISOL)
  115          SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+RATEC
            END IF
         END IF
  500 CONTINUE
C
cgzh 4/26/07 tributary flow into lake but no outflow into stream?  process here
C loop over lakes
      DO 900 LK=1,NLAKES
C if lake concentration has not been calculated,
        IF (LKDONE(LK).EQ.0) THEN
C check to see if there are inflowing tributaries (if so set INTRIB=1)
         INTRIB=0
          DO 910 ITRIB=1,NTRB
            INODE=ITRB(LK,ITRIB)
            IF (INODE.LE.0) GO TO 910
            INTRIB=1
 910      CONTINUE
C if there are inflowing tributaries, update lake concentration
          IF(INTRIB.EQ.1) THEN       
               CALL CLAK3BD(KLK,CLKSUM,NLAKES,NSS,NTRB,
     1      ITRB,STRIN,NSOL,TIMV,KCNT,MSUB,MSUB1,LKDONE,
cgzh debug  send in LK for MGFLG
C   2      SOLPPT,CSWIN,CSWOUT,CLAKE,CWDRAW,MDFLG,NSTRM,
     2      SOLPPT,CSWIN,CSWOUT,CLAKE,CWDRAW,LK,NSTRM,
     3      VOL,CGWIN,CGWOUT,COUT,CSRUN,CSLAKE,IOUTS,NSTRMAR,CNTRIB)
          ELSE
C otherwise, print warning
            WRITE(IOUTS,925) LK
          END IF
	  END IF
C
 900   CONTINUE
 925   FORMAT(/1X,'*** WARNING *** LAKE NUMBER',I3,' WAS NOT UPDATED:',
     & ' TRIBUTARY CHECK')
C
  600 CONTINUE
C
C-------CHECK FLAGS FOR OUTPUT
      IPRNT=0
  625 IF(ISTCB1.EQ.0) GO TO 800
      IF (NPNTCL.EQ.0.AND.KKPER.EQ.NPER.AND.KKSTP.EQ.NSTP.
     *AND.IIMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-2.AND.KKSTP.EQ.NSTP.AND.IIMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-1.AND.IIMOV.EQ.NMOV) IPRNT=1
      IF(NPNTCL.GE.1) THEN
         IF(MOD(IIMOV,NPNTCL).EQ.0) IPRNT=1
      ENDIF
      IF (KKPER.EQ.NPER.AND.KKSTP.EQ.NSTP.AND.IIMOV.EQ.NMOV) IPRNT=1
      IF (IPRNT.EQ.0) GO TO 800
      IF(IPTFLG.GT.0) GO TO 800
C30-----PRINT STREAMFLOW RATES AND CONCENTRATIONS FOR EACH REACH.      C
      WRITE(IOUTS,650)
  650 FORMAT(///6X,'STREAM PACKAGE INFORMATION: '//
     * 11X,'LAYER',5X,'ROW',4X,'COLUMN',3X,'STREAM',3X,
     1'REACH',5X,'FLOW INTO',4X,'FLOW OUT OF',4X,'CONC. OUT OF'/37X,
     2      'NUMBER',3X,'NUMBER',3X,'STREAM REACH',2X,'STREAM REACH',
     3      3X,'STREAM REACH'/)
      DO 690 L=1,NSTRM
      IL=ISTRM(1,L)
      IR=ISTRM(2,L)
      IC=ISTRM(3,L)
      WRITE(IOUTS,675)IL,IR,IC,ISTRM(4,L),ISTRM(5,L),
     1     STRM(10,L),STRM(9,L),COUT(L,1)
  675 FORMAT(1X,5X,5I9,5X,1PE10.3,4X,E10.3,5X,E10.3)
      IF (NSOL.GT.1) THEN
         DO 680 ISOL=2,NSOL
  680       WRITE (IOUTS,685) COUT(L,ISOL)
  685 FORMAT(72X,'SOLUTE No. ',I3,2X,'COUT = ',G9.3)
      ENDIF
  690 CONTINUE
      WRITE (IOUTS,88)
   88 FORMAT (/)
  800 CONTINUE
C                                                                      C
C31-----RETURN.                                                        C
      RETURN
      END
