C  SOLUTE SINK/SOURCES AT WELLS 
C  
C     ***************************************************************
C
      SUBROUTINE GWT1WEL5BD(NWELLS,MXWELL,
     *  WELL,IBOUND,NCOL,NROW,NLAY,KSTP,KPER,IPERGWT,
     *  SRCFLO,SRCSOL,SNKFLO,
     *  NSCOL,NSROW,NSLAY,IOUTS,NWELVL,IWELAL,IWELLC,ICONLY,NOPRWL,
     *  VC,VR,VL,IFACEC)
C
C     ***************************************************************
C
C  CWEL5FM  FIND OUT WHICH OF THE WELL AUXILIARY PARAMETERS IS CONC
C  IF NOT FOUND, AND WELLS ARE WITHIN SUBGRID, THEN STOP  
C
      DIMENSION VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1)  
      CHARACTER*8 PACKAGE
C
      DIMENSION WELL(NWELVL,MXWELL)
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY)
      COMMON /WELCOM/welaux(20)
      CHARACTER*16 WELAUX
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
C
C  CHECK EXISTENCE OF CBCALLOCATE OPTION
      IF (IWELAL.NE.1) THEN
         WRITE(IOUTS,*) ' ERROR, CBCALLOCATE OPTION MUST BE ',
     *        'SPECIFIED IN WELL INPUT DATA'
         STOP
      ENDIF
C
C RETURN IF NOT SOLVING TRANSPORT YET
      IF(KPER.LT.IPERGWT) RETURN
C  SET POINTER TO CONCENTRATION DATA FOR WELLS
C  SET POINTER TO IFACE VALUE FOR WELLS
      IF(KPER.GE.IPERGWT.AND.KSTP.EQ.1) THEN
C COUNT FROM THE TOTAL NUMBER OF STANDARD DATA INPUTS FOR WELLS (4)
         NAUX=NWELVL-4-IWELAL
         IWELLC=0
         IFACEC=0
         IF(NAUX.GT.0) THEN
            DO 8 IAUX=1,NAUX
               IF(WELAUX(IAUX).EQ.'CONCENTRATION'.OR.
     *            WELAUX(IAUX).EQ.'CONC') THEN
                  IWELLC=IAUX+4
               END IF
               IF(WELAUX(IAUX).EQ.'IFACE') THEN
                  IFACEC=IAUX+4
               END IF
               IF(WELAUX(IAUX).EQ.'TEST1') THEN
                  Itest1=IAUX+4
               END IF
               IF(WELAUX(IAUX).EQ.'TEST2') THEN
                  Itest2=IAUX+4
               END IF
               IF(WELAUX(IAUX).EQ.'TEST3') THEN
                  Itest3=IAUX+4
               END IF
               IF(WELAUX(IAUX).EQ.'TEST4') THEN
                  Itest4=IAUX+4
               END IF
 8          CONTINUE
         END IF
         IF(IWELLC.EQ.0) WRITE(IOUTS,*) ' WARNING, ',
     *         'NO CONCENTRATION DATA READ FOR WELL NODES'
         IF(IFACEC.EQ.0) WRITE(IOUTS,*) ' IFACE NOT DEFINED: ALL WELLS', 
     *         ' APPLIED AS DISTRIBUTED SOURCES/SINKS'
      END IF
C
C  RETURN IF NO ACTIVE WELL NODES THIS PERIOD
      IF(NWELLS.LE.0) RETURN
C
C5------PRINT WELL CONCENTRATIONS AND FORMULATE
C        SINK AND SOURCE TERMS FOR TRANSPORT
      IWELSG=0
      IERSRC=0
      ITHKSG=0
      IPRNT=0
      DO 10 IW=1,NWELLS
         K=WELL(1,IW)
         KS=K-ISLAY1+1
         IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 10
         I=WELL(2,IW)
         IS=I-ISROW1+1
         IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 10
         J=WELL(3,IW)
C  SKIP WELLS IN FIXED HEAD CELLS
         IF(IBOUND(J,I,K).LE.0) THEN
           ITHKSG=ITHKSG+1
           GO TO 10
         ENDIF
         JS=J-ISCOL1+1
         IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 10
         IF(KSTP.EQ.1) THEN
           IF(IFACEC.LE.0) THEN
             IF(IWELLC.LE.0) THEN
               IF(IPRNT.EQ.0.AND.NOPRWL.NE.1) WRITE(IOUTS,29)
   29          FORMAT(/1H ,' WELL CONCENTRATIONS FOR THIS STRESS PERIOD'
     1         /'  (ONLY FOR WELL NODES WITHIN TRANSPORT SUBGRID)'//
     2         1H ,1X,'WELL NO.  LAYER    ROW   COLUMN    VOL/T   ',
     3         /1X,68('-'))
               IPRNT=1
               IF (NOPRWL.NE.1) WRITE(IOUTS,'(1X,I6,3I8,2X,1PE12.4)') 
     1                         IW,K,I,J,WELL(4,IW)
             ELSE
               IF(IPRNT.EQ.0.AND.NOPRWL.NE.1) WRITE(IOUTS,129)
  129          FORMAT(/1H ,' WELL CONCENTRATIONS FOR THIS STRESS PERIOD'
     1         /'  (ONLY FOR WELL NODES WITHIN TRANSPORT SUBGRID)'//
     2         1H ,1X,'WELL NO.  LAYER    ROW   COLUMN    VOL/T   ',
     3         ' SOURCE CONC    MASS/T'/1X,68('-'))
               IPRNT=1
                RATEC=WELL(4,IW)*WELL(IWELLC,IW)
                IF (NOPRWL.NE.1) WRITE(IOUTS,'(1X,I6,3I8,2X,1P3E12.4)') 
     1		                 IW,K,I,J,WELL(4,IW),WELL(IWELLC,IW),RATEC
             END IF
           ELSE
             IFACE=WELL(IFACEC,IW)
             IF(IWELLC.LE.0) THEN
               IF(IPRNT.EQ.0.AND.NOPRWL.NE.1) WRITE(IOUTS,229)
  229          FORMAT(/1H ,' WELL CONCENTRATIONS FOR THIS STRESS PERIOD'
     1         /'  (ONLY FOR WELL NODES WITHIN TRANSPORT SUBGRID)'//
     2         1H ,1X,'WELL NO.  LAYER    ROW   COLUMN    VOL/T   ',
     3         '  IFACE    '/1X,68('-'))
               IPRNT=1
               IF (NOPRWL.NE.1) WRITE(IOUTS,'(1X,I6,3I8,2X,1PE12.4,I8)') 
     1                         IW,K,I,J,WELL(4,IW),IFACE
             ELSE
               IF(IPRNT.EQ.0.AND.NOPRWL.NE.1) WRITE(IOUTS,329)
  329          FORMAT(/1H ,' WELL CONCENTRATIONS FOR THIS STRESS PERIOD'
     1         /'  (ONLY FOR WELL NODES WITHIN TRANSPORT SUBGRID)'//
     2         1H ,1X,'WELL NO.  LAYER    ROW   COLUMN    VOL/T   ',
     3         ' SOURCE CONC    MASS/T     IFACE'/1X,76('-'))
               IPRNT=1
                RATEC=WELL(4,IW)*WELL(IWELLC,IW)
              IF (NOPRWL.NE.1) WRITE(IOUTS,'(1X,I6,3I8,2X,1P3E12.4,I8)') 
     1		                 IW,K,I,J,WELL(4,IW),WELL(IWELLC,IW),RATEC,
     2                         IFACE
             END IF
           END IF
         END IF
C
         IF(ICONLY.EQ.1) THEN
            IWELSG=IWELSG+1
         ELSE
            RATE=WELL(NWELVL,IW)
            IFACE=0
            IF(IFACEC.GT.0) IFACE=WELL(IFACEC,IW)
            IF(IFACE.LT.0) IFACE=-1
            IF(IFACE.NE.0) THEN
              PACKAGE='    WELL'
              IPCK=IW
              CALL GWT1BFLX5PCK(IFACE,IBOUND,RATE,PACKAGE,IPCK,
     *         VC,VR,VL,
     *         NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *         IOUTS,
     *         J,I,K,JS,IS,KS)
            END IF
C
            IF(RATE.LT.0.0) THEN
               SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+RATE
            ELSE IF(RATE.GT.0.0) THEN
               IF(IWELLC.LE.0) THEN
                 IERSRC=IERSRC+1
               ELSE
                 SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+RATE
                 RATEC=RATE*WELL(IWELLC,IW)
                 SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+RATEC
               END IF
            END IF
         END IF
cgzh debug
c      write(*,*) 'test1',well(itest1,iw)
c      write(*,*) 'test2',well(itest2,iw)
c      write(*,*) 'test3',well(itest3,iw)
c      write(*,*) 'test4',well(itest4,iw)
   10 CONTINUE
C
C  STOP IF WELL WITHIN TRANSPORT SUBGRID AND ICONLY=1
      IF(IWELSG.GT.0) THEN
         WRITE(IOUTS,*) ' ERROR,',IWELSG,' WELL NODES WITHIN',
     *      ' SUBGRID, BUT ICONLY=1, STOPPING'
         STOP
      END IF
C  STOP IF SOURCE WELL NODE WITHIN SUBGRID, BUT NO CONC DATA
      IF(IERSRC.GT.0) THEN
         WRITE(IOUTS,*) ' ERROR:',IERSRC,' SOURCE WELL NODES',
     *     ' WITHIN SUBGRID, BUT NO CONCENTRATION DATA, STOPPING'
         STOP
      END IF
C  WARN IF WELL NODE WITHIN SUBGRID, BUT NO FLOW ASSOC WITH IT
      IF(ITHKSG.GT.0) THEN
         WRITE(IOUTS,*) ' WARNING:',ITHKSG,' WELL NODES',
     *     ' SKIPPED WITHIN SUBGRID (IBOUND <= 0)'
      END IF
C
      RETURN
      END
