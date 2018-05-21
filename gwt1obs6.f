C     Last change: RBW Dec. 8, 2014
C     Support for weighted particles (MOCWT) added.
C  OBS5 OBSERVATION WELLS  
C  8/95
C
C  GWT1OBS5DF READ NUMBER OF OBSERVATION WELLS 
C
C     ******************************************************************
C
      SUBROUTINE GWT1OBS5DF(NUMOBS,
     *                  IOUTS,INOBS,IOBSFL)
C
C     ******************************************************************
C
C     READ IN NUMOBS
C     ******************************************************************
C
      READ(INOBS,*) NUMOBS, IOBSFL
      RETURN
      END
C
C
C GWT1OBS5AL ALLOCATE SPACE FOR OBSERVATION WELLS
C
C     ******************************************************************
C
      SUBROUTINE GWT1OBS5AL(ISUM,LSOBSW,NUMOBS,
     *                  IOUT,IOUTS)    
C
C     ******************************************************************
C
C     ARRAY IS LAY,ROW,COL,UNIT#
      LSOBSW=ISUM
      ISUM=ISUM+NUMOBS*4
      ISP=NUMOBS*4
      WRITE(IOUTS,101) ISP
  101 FORMAT(1X,I8,' ELEMENTS IN X ARRAY ARE USED BY OBS')
C memory checks, lenx not needed for dynamic memory allocation
c      ISUM1=ISUM-1
c      WRITE(IOUT,102) ISUM1,LENX
c  102 FORMAT(1X,I8,' ELEMENTS OF X ARRAY USED OUT OF',I8)
c      IF(ISUM1.GT.LENX) WRITE(IOUT,103)
c  103 FORMAT(1X,'   ***X ARRAY MUST BE DIMENSIONED LARGER***')
C
      RETURN
      END
C
C
C  OBS5RP READ OBSERVATION WELL INPUT FILE 
C
C     ******************************************************************
C
      SUBROUTINE GWT1OBS5RP(IOBLST,NUMOBS,
     *                  IOUTS,INOBS,IOBSFL)
C
C     ******************************************************************
C
C     READ OBSERVATION WELL LOCATIONS
C     ******************************************************************
C
      DIMENSION IOBLST(4,NUMOBS)               
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ******************************************************************
C
      WRITE (IOUTS,140) NUMOBS
      WRITE (IOUTS,150)
C READ THE FIRST RECORD
      IOB=1
         READ(INOBS,*) IOBLST(1,IOB),IOBLST(2,IOB),IOBLST(3,IOB),
     *                 IOBLST(4,IOB)
         K=IOBLST(1,IOB)
         I=IOBLST(2,IOB)
         J=IOBLST(3,IOB)
         IOBUN=IOBLST(4,IOB)
C
         WRITE(IOUTS,'(5I8,5X,A40)') IOB,K,I,J,IOBUN            
         IF(K.LT.ISLAY1.OR.K.GT.ISLAY2.OR. 
     *      I.LT.ISROW1.OR.I.GT.ISROW2.OR. 
     *      J.LT.ISCOL1.OR.J.GT.ISCOL2) THEN
            WRITE(IOUTS,*) '***ERROR***   OBSERVATION WELL OUTSIDE',
     *                      ' SUBGRID'
            STOP
         ENDIF
C CYCLE THROUGH THE REMAINING RECORDS, READ ACCORDING TO IOBSFL
      DO 139 IOB=2,NUMOBS
        IF (IOBSFL.LE.0) THEN 
         READ(INOBS,*) IOBLST(1,IOB),IOBLST(2,IOB),IOBLST(3,IOB),
     *                 IOBLST(4,IOB)
         IOBUN=IOBLST(4,IOB)
        ELSE
         READ(INOBS,*) IOBLST(1,IOB),IOBLST(2,IOB),IOBLST(3,IOB)
        ENDIF
         K=IOBLST(1,IOB)
         I=IOBLST(2,IOB)
         J=IOBLST(3,IOB)
C
         WRITE(IOUTS,'(5I8,5X,A40)') IOB,K,I,J,IOBUN            
         IF(K.LT.ISLAY1.OR.K.GT.ISLAY2.OR. 
     *      I.LT.ISROW1.OR.I.GT.ISROW2.OR. 
     *      J.LT.ISCOL1.OR.J.GT.ISCOL2) THEN
            WRITE(IOUTS,*) '***ERROR***   OBSERVATION WELL OUTSIDE',
     *                      ' SUBGRID'
            STOP
         ENDIF
C
C      FLAG = 0 : SEPARATE OUTPUT FILES
C      FLAG = 1 : ONE OUPUT FILE (HEADER WRITTEN LATER)
  139 CONTINUE
      IF (IOBSFL.GE.1) THEN
            IOBUN=IOBLST(4,1)
            WRITE(IOUTS,160) IOBUN
            WRITE(IOUTS,'(/)')
      ELSE
            WRITE(IOUTS,'(/)')
            WRITE(IOUTS,*) 'OBSERVATION WELL DATA WILL BE'
     *                     ,' WRITTEN ON UNIT NUMBERS LISTED ABOVE' 
      ENDIF
  140 FORMAT(///'COORDINATES FOR',I4,' OBSERVATION WELLS:')
  150 FORMAT(/'  WELL #   LAYER     ROW  COLUMN    ',
     *                           'UNIT')
  160 FORMAT('ALL OBSERVATION WELL DATA WILL BE WRITTEN ON UNIT ',
     *        I3)
      RETURN
      END       
C
C
C*************************************************************************
C SOBS5O   PRINT OBSERVATION WELL DATA
C     *************************************************************
      SUBROUTINE SOBS5O(HNEW,CONC,SUMTCH,IOBLST,IMOV,
     *                  NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,NUMOBS,IOBSFL,
     *                  SUMWT,NPCELL)
C     ***************************************************************
cgzh debug zinn
      DOUBLE PRECISION HNEW
      DOUBLE PRECISION SUMWT
      CHARACTER*50  LFRMAT
C
      DIMENSION HNEW(NCOL,NROW,NLAY),CONC(NSCOL,NSROW,NSLAY),
     *          IOBLST(4,NUMOBS)
      DIMENSION NPCELL(NSCOL,NSROW,NSLAY),SUMWT(NSCOL,NSROW,NSLAY)
C     ***************************************************************
C
      ALLOCATABLE OBTEMP(:,:)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c make this 2 after zinn
	ALLOCATE (OBTEMP(4,NUMOBS))
C
C  LOOP OVER OBS WELLS
      DO 10 IOB=1,NUMOBS
         K=IOBLST(1,IOB)
         KS=K-ISLAY1+1
         I=IOBLST(2,IOB)
         IS=I-ISROW1+1
         J=IOBLST(3,IOB)
         JS=J-ISCOL1+1
         IF(IOBSFL.EQ.0.OR.IOB.EQ.1) IOBUN=IOBLST(4,IOB)
C  WRITE HEADER
         IF(IMOV.EQ.0) THEN
            IF(IOBSFL.EQ.0) THEN
              WRITE(IOBUN,*) '"OBSERVATION WELL DATA"'
              WRITE(IOBUN,*) '"NODE  (K,I,J)  ',K,I,J,'"'
              WRITE(IOBUN,*) '"   TIME               HEAD        CONC."'
            ELSE IF(IOB.EQ.1) THEN
              WRITE(IOBUN,*) '"OBSERVATION WELL DATA"'
              WRITE(IOBUN,*) '"TIME, THEN HEAD AND CONC. FOR ',
     *                              'EACH OBS. WELL AT NODE (K,I,J)"'
            ENDIF
         ENDIF
         OBTEMP(1,IOB)=HNEW(J,I,K)
         OBTEMP(2,IOB)=CONC(JS,IS,KS)
cgzh debug zinn
C         OBTEMP(3,IOB)=NPCELL(JS,IS,KS)
C         OBTEMP(4,IOB)=SUMWT(JS,IS,KS)
         IF(IOBSFL.EQ.0) 
     *      WRITE(IOBUN,45) SUMTCH,OBTEMP(1,IOB),OBTEMP(2,IOB)
  10  CONTINUE
C
C  ALL TO ONE FILE
      IF(IOBSFL.GE.1) THEN
         IUNIT=IOBLST(4,1)
         IF(IMOV.EQ.0) THEN
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,15) NUMOBS
  15  FORMAT('(A20,',I4,'(A12,2(I3,'',''),I3,A1), ''"'')')
            WRITE(IUNIT,LFRMAT) ' "   TIME:          ',('   H & C AT ',
     *         IOBLST(1,IOB),IOBLST(2,IOB),IOBLST(3,IOB),' ',
     *         IOB=1,NUMOBS)
         ENDIF
C  WRITE TIME: HEADS + CONCENTRATIONS (NUMOBS*2 IS NUMBER OF DATA, H AND C)
         WRITE(LFRMAT,25) NUMOBS*2
cgzh debug zinn (4 spots)
corig         WRITE(LFRMAT,25) NUMOBS*4
cgzh debug more prec for wellbore
corig  25  FORMAT('(1PE11.4,1P',I4,'(1X,E11.3))')
  25  FORMAT('(1PE19.8,1P',I4,'(1X,E11.3))')
cgzh debug temp OBS format change--more precision needed for testing
c  25  FORMAT('(1PE11.4,1P',I4,'(1X,E12.6))')
         WRITE(IUNIT,LFRMAT) SUMTCH,
     *            (OBTEMP(1,IOB),OBTEMP(2,IOB),IOB=1,NUMOBS)
cgzh debug zinn (4 spots)
c     *            (OBTEMP(1,IOB),OBTEMP(2,IOB),
c     *     OBTEMP(3,IOB),OBTEMP(4,IOB),IOB=1,NUMOBS)
      ENDIF
  45  FORMAT(1PE19.8,2E12.4)
  50  FORMAT(3X,'H & C AT ',3I3,1X)
      DEALLOCATE(OBTEMP)
      RETURN
      END
