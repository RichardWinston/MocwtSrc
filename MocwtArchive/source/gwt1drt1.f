C  CDRT1  SOLUTE SINKS AT DRAIN RETURNS 
C     ******************************************************************
C
      SUBROUTINE GWT1DRT1BD(NDRTCL,MXDRT,DRTF,NCOL,NROW,NLAY,
     &                      IBOUND,IOUTS,NDRTVL,IDRTAL,IDRTFL,
     &                      SRCFLO,SRCSOL,SNKFLO,CONC,
     &                      NSCOL,NSROW,NSLAY,ICONLY,KSTP,KPER,NOPRDT,
     &                      VC,VR,VL,IPERGWT,IFACEC)
C
C     ******************************************************************
C
      DIMENSION VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1)  
      CHARACTER*8 PACKAGE
C
      DIMENSION DRTF(NDRTVL,MXDRT),IBOUND(NCOL,NROW,NLAY)
      DIMENSION SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),CONC(NSCOL,NSROW,NSLAY)
C
      COMMON /DRTCOM/DRTAUX(5)
      CHARACTER*16 DRTAUX
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ******************************************************************
C RETURN IF NOT SOLVING TRANSPORT YET
      IF(KPER.LT.IPERGWT) RETURN
C  SET POINTER TO IFACE VALUE FOR DRAIN RETURNS
      IF(KPER.EQ.IPERGWT.AND.KSTP.EQ.1) THEN
C COUNT FROM THE TOTAL NUMBER OF STANDARD DATA INPUTS FOR DRAINS (5)
         NAUX=NDRTVL-5-IDRTAL
         IFACEC=0
         IF(NAUX.GT.0) THEN
           DO 8 IAUX=1,NAUX
             IF(DRTAUX(IAUX).EQ.'IFACE') THEN
               IFACEC=IAUX+5
             END IF
 8         CONTINUE
         END IF
         IF(IFACEC.EQ.0) WRITE(IOUTS,*) ' IFACE NOT DEFINED: ', 
     *         'ALL DRAIN RETURNS APPLIED AS DISTRIBUTED SINKS'
      END IF
C
C
C  RETURN IF NO ACTIVE DRAIN RETURNS THIS PERIOD
      IF(NDRTCL.LE.0) RETURN
C
C  ADD DRAINS AND RETURNS
      IDRTSG=0
      ITHKSG=0
      IPRNT=0
C  CHECK EXISTENCE OF CBCALLOCATE OPTION
      IF (IDRTAL.EQ.0) THEN
         WRITE(IOUTS,*) ' ERROR, CBCALLOCATE OPTION MUST BE ',
     *        'SPECIFIED IN DRT INPUT DATA'
         STOP
      ENDIF
C
      DO 10 ID=1,NDRTCL
         K=DRTF(1,ID)
         KS=K-ISLAY1+1
         IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 10
         I=DRTF(2,ID)
         IS=I-ISROW1+1
         IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 10
         J=DRTF(3,ID)
         JS=J-ISCOL1+1
         IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 10
C  SKIP DRAINS IN FIXED HEAD CELLS
         IF(IBOUND(J,I,K).LE.0) THEN
           ITHKSG=ITHKSG+1
           GO TO 10
         ENDIF
c drain returns (sources)
         IF (IDRTFL.GT.0) THEN
           KK = DRTF(6,ID)
           KKS=KK-ISLAY1+1
           IF(KKS.LT.1.OR.KKS.GT.NSLAY) GO TO 10
           IF (KK.NE.0) THEN
            II = DRTF(7,ID)
            IIS=II-ISROW1+1
            IF(IIS.LT.1.OR.IIS.GT.NSROW) GO TO 10
            JJ = DRTF(8,ID)
            JJS=JJ-ISCOL1+1
            IF(JJS.LT.1.OR.JJS.GT.NSCOL) GO TO 10
	     ENDIF
	   ENDIF
         IF(KSTP.EQ.1) THEN
           IF(IFACEC.LE.0) THEN
            IF(IPRNT.EQ.0.AND.NOPRDT.NE.1) WRITE(IOUTS,1000)
 1000      FORMAT(//' DRAIN RETURNS WITHIN SUBGRID THIS STRESS PERIOD'//
     1      '  DRAIN NO.   LAYER      ROW     COLUMN'/
     2      1X,40('-'))
            IPRNT=1
            IF(NOPRDT.NE.1) WRITE(IOUTS,'(1X,I6,3I10)') ID,K,I,J
           ELSE
            IFACE=DRTF(IFACEC,ID)
            IF(IPRNT.EQ.0.AND.NOPRDT.NE.1) WRITE(IOUTS,1001)
 1001      FORMAT(//' DRAIN RETURNS WITHIN SUBGRID THIS STRESS PERIOD'//
     1      '  DRAIN NO.   LAYER      ROW     COLUMN    IFACE'/
     2      1X,50('-'))
            IPRNT=1
            IF(NOPRDT.NE.1) WRITE(IOUTS,'(1X,I6,4I10)') ID,K,I,J,IFACE
           END IF
         END IF
         IF(ICONLY.EQ.1) THEN
            IDRTSG=IDRTSG+1
         ELSE
c drain return
            RATE=DRTF(NDRTVL,ID)
            IFACE=0
            IF(IFACEC.GT.0) IFACE=DRTF(IFACEC,ID)
            IF(IFACE.LT.0) IFACE=-1
            IF(IFACE.NE.0) THEN
              PACKAGE=' DRT-SNK'
              IPCK=ID
              CALL GWT1BFLX5PCK(IFACE,IBOUND,RATE,PACKAGE,IPCK,
     *         VC,VR,VL,
     *         NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *         IOUTS,
     *         J,I,K,JS,IS,KS)
            END IF
C
            SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+RATE
c drain return
c use CONC at donor cell
	      IF(IDRTFL.GT.0) THEN
	        RATE=DRTF(NDRTVL-1,ID)
              IF(IFACE.NE.0) THEN
                PACKAGE=' DRT-SRC'
                CALL GWT1BFLX5PCK(IFACE,IBOUND,RATE,PACKAGE,IPCK,
     *           VC,VR,VL,
     *           NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *           IOUTS,
     *           J,I,K,JS,IS,KS)
              END IF
              SRCFLO(JJS,IIS,KKS)=SRCFLO(JJS,IIS,KKS)+RATE
              RATEC=RATE*CONC(JS,IS,KS)
              SRCSOL(JJS,IIS,KKS)=SRCSOL(JJS,IIS,KKS)+RATEC
	      END IF
	   END IF
   10 CONTINUE
C
      IF(IDRTSG.GT.0) THEN
         WRITE(IOUTS,*) ' ERROR,',IDRTSG,' DRAINS WITHIN SUBGRID, ',
     *        'BUT ICONLY=1, STOPPING'
         STOP
      END IF
C  WARN IF DRAIN NODE WITHIN SUBGRID, BUT NO FLOW ASSOC WITH IT
      IF(ITHKSG.GT.0) THEN
         WRITE(IOUTS,*) ' WARNING:',ITHKSG,' DRAIN NODES',
     *     ' SKIPPED WITHIN SUBGRID (IBOUND <= 0)'
      ENDIF
C
      RETURN
      END
