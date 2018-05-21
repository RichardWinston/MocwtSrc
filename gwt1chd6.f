      SUBROUTINE GWT1CHD6RP(CHDS,NCHDS,MXCHD,IBOUND,NCOL,NROW,NLAY,
     1          NSCOL,NSROW,NSLAY,IOUTS,NCHDVL,NNPCHD,CFXHBC,KKPER,
     *          IPERGWT,ICONLY,NOPRCH)
C
C     ******************************************************************
C     READ CONCENTRATION DATA FOR CHD
C     ******************************************************************
C
C     SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION CHDS(NCHDVL,MXCHD),IBOUND(NCOL,NROW,NLAY)
      DIMENSION CFXHBC(NSCOL,NSROW,NSLAY)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      COMMON /CHDCOM/CHDAUX(5)
      CHARACTER*16 CHDAUX
C     ------------------------------------------------------------------
C
C RETURN IF NOT SOLVING TRANSPORT YET
      IF(KKPER.LT.IPERGWT) RETURN
C  SET POINTER TO CONCENTRATION DATA FOR CHD
      IF(KKPER.EQ.IPERGWT) THEN
C  COUNT FROM THE TOTAL NUMBER OF STANDARD DATA INPUTS (5)
         NAUX=NCHDVL-5
         ICHDC=0
         IF(NAUX.GT.0) THEN
            DO 8 IAUX=1,NAUX
               IF(CHDAUX(IAUX).EQ.'CONCENTRATION'.OR.
     *            CHDAUX(IAUX).EQ.'CONC') THEN
                  ICHDC=IAUX+5
                  GOTO 9
               END IF
 8          CONTINUE
            WRITE(IOUTS,*) ' WARNING, NO CONCENTRATION DATA READ',
     *         ' FOR CHD NODES; ZONCON VALUES WILL BE USED'
 9          CONTINUE
         END IF
      END IF
C
C  RETURN IF NO ACTIVE CHD NODES THIS PERIOD
      IF(NNPCHD.LE.0) RETURN
C
C  SAVE CHD CONCENTRATIONS IN CFXHBC
C        
      ICHDSG=0
      IPRNT=0
      DO 10 IC=1,NNPCHD
         K=CHDS(1,IC)
         KS=K-ISLAY1+1
         IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 10
         I=CHDS(2,IC)
         IS=I-ISROW1+1
         IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 10
         J=CHDS(3,IC)
         JS=J-ISCOL1+1
         IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 10
           IF(IPRNT.EQ.0.AND.NOPRCH.NE.1) WRITE(IOUTS,29)
   29      FORMAT(/1H ,' CHD CONCENTRATIONS FOR THIS STRESS PERIOD'
     1     /'  (ONLY FOR CHD NODES WITHIN TRANSPORT SUBGRID)'//
     2     1H ,1X,'  LAYER    ROW   COLUMN  ',
     3     ' SOURCE CONC '/1X,68('-'))
           IPRNT=1
           IF(ICHDC.GT.0) THEN
	       CFXHBC(JS,IS,KS)=CHDS(ICHDC,IC)
             IF(NOPRCH.NE.1) WRITE(IOUTS,'(1X,3I8,2X,1P3E12.4)') K,I,J,
     *           CHDS(ICHDC,IC)
           END IF
         IF(ICONLY.EQ.1) THEN
            IWELSG=IWELSG+1
         END IF
   10 CONTINUE
C
C  STOP IF WELL WITHIN TRANSPORT SUBGRID AND ICONLY=1
      IF(ICHDSG.GT.0) THEN
         WRITE(IOUTS,*) ' ERROR,',ICHDSG,' CHD NODES WITHIN',
     *      ' SUBGRID, BUT ICONLY=1, STOPPING'
         STOP
      END IF
C
      RETURN
      END
