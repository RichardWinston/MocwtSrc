C
C  SMOC5Z    SET ARRAY TO CONSTANT VALUE
C*************************************************************************
C
      SUBROUTINE SMOC5Z(ARRAY,NSCOL,NSROW,NSLAY,VALUE)
C
C*************************************************************************
C
      DIMENSION ARRAY(NSCOL,NSROW,NSLAY)
C
      DO 10 KS=1,NSLAY
      DO 10 IS=1,NSROW
      DO 10 JS=1,NSCOL
      ARRAY(JS,IS,KS)=VALUE
10    CONTINUE
C
      RETURN
      END
C
CC
C  SMOC5ZD    SET DP ARRAY TO CONSTANT VALUE
C*************************************************************************
C
      SUBROUTINE SMOC5ZD(ARRAY,NSCOL,NSROW,NSLAY,VALUE)
C
C*************************************************************************
C
      DOUBLE PRECISION ARRAY,VALUE
      DIMENSION ARRAY(NSCOL,NSROW,NSLAY)
C
      DO 10 KS=1,NSLAY
      DO 10 IS=1,NSROW
      DO 10 JS=1,NSCOL
      ARRAY(JS,IS,KS)=VALUE
10    CONTINUE
C
      RETURN
      END
C
C***************************************************************
C
      SUBROUTINE SMOC6O(INUNIT,INMOC,IOUTS,JUNIT,DUNIT,MOCTYPE,NIUNIT)
C
C-----from VERSION 0818 15JULY1993 SBAS5O
C     ******************************************************************
C     OPEN FILES.
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      DIMENSION JUNIT(NIUNIT)
      CHARACTER*4 DUNIT(NIUNIT)
      CHARACTER*80 LINE
      CHARACTER*7 FILSTAT
      CHARACTER*20 FILACT, FMTARG, ACCARG
	LOGICAL LOP
      INCLUDE 'openspec.inc'
CMOCWT
      INCLUDE 'ptwt.inc'
C     ---------------------------------------------------------------
C
C1------INITIALIZE CONSTANTS.  ILIST IS SET TO 1 ONCE THE LISTING
C1------FILE HAS BEEN OPENED; UNTIL THEN ERROR MESSAGES ARE WRITTEN
C1------ TO "*" UNIT.
      ILIST=0
      NFILE=0
      IOUTS=0
      INMOC=0
      MOCTYPE=0
      DO 5 I=1,40
      JUNIT(I)=0
5     CONTINUE
C
C2------READ A LINE; IGNORE BLANK LINES AND PRINT COMMENT LINES.
10    READ(INUNIT,'(A)',END=1000) LINE
      IF(LINE.EQ.' ') GO TO 10
      IF(LINE(1:1).EQ.'#') THEN
        IF(NFILE.NE.0.and.IOUTS.NE.0) WRITE(IOUTS,'(A)') LINE
        GO TO 10
      END IF
C
C3------DECODE THE FILE TYPE AND UNIT NUMBER.
      LLOC=1
      CALL URWORD(LINE,LLOC,ITYP1,ITYP2,1,N,R,IOUTS,INUNIT)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IU,R,IOUTS,INUNIT)
C
C4------CHECK FOR A VALID FILE TYPE.
      FMTARG='FORMATTED'
      ACCARG='SEQUENTIAL'
      FILSTAT='UNKNOWN'
      FILACT=' '
C
C4A-----FIRST ENTRY MUST BE FILE-TYPE "CLST".
      IF(NFILE.EQ.0) THEN
         IF(LINE(ITYP1:ITYP2).NE.'CLST') THEN
            WRITE(*,*) ' FIRST ENTRY IN GWT NAME FILE MUST BE "CLST".'
            STOP
	   ELSE
            FILSTAT='REPLACE'
	      ILIST=1
            IOUTS=IU
         END IF
C
C4B-----CHECK FOR "GWT" FILE TYPE.
      ELSEIF(LINE(ITYP1:ITYP2).EQ.'MOC'.AND.MOCTYPE.EQ.0) THEN
         INMOC=IU
         MOCTYPE=1
CMOCWT
         PTWTON=0
         FILSTAT='OLD    '
         FILACT=ACTION(1)
      ELSEIF(LINE(ITYP1:ITYP2).EQ.'MOCIMP'.AND.MOCTYPE.EQ.0) THEN
         INMOC=IU
         MOCTYPE=2
CMOCWT
         PTWTON=0
         FILSTAT='OLD    '
         FILACT=ACTION(1)
      ELSEIF(LINE(ITYP1:ITYP2).EQ.'ELLAM'.AND.MOCTYPE.EQ.0) THEN
         INMOC=IU
         MOCTYPE=3
CMOCWT
         PTWTON=0
         FILSTAT='OLD    '
         FILACT=ACTION(1)
      ELSEIF(LINE(ITYP1:ITYP2).EQ.'MOCWT'.AND.MOCTYPE.EQ.0) THEN
         INMOC=IU
         MOCTYPE=1
CMOCWT
         PTWTON=1
         FILSTAT='OLD    '
         FILACT=ACTION(1)
      ELSEIF(LINE(ITYP1:ITYP2).EQ.'MOCWTI'.AND.MOCTYPE.EQ.0) THEN
         INMOC=IU
         MOCTYPE=2
CMOCWT
         PTWTON=1
         FILSTAT='OLD    '
         FILACT=ACTION(1)
C
C4C-----CHECK FOR "UNFORMATTED" FILE TYPE.
      ELSE IF(LINE(ITYP1:ITYP2).EQ.'DATA(BINARY)') THEN
         FMTARG=FORM
         ACCARG=ACCESS
C
C4D-----CHECK FOR "FORMATTED" FILE TYPE.
      ELSE IF(LINE(ITYP1:ITYP2).EQ.'DATA') THEN
         FMTARG='FORMATTED'
C
C4E-----CHECK FOR MAJOR OPTIONS.
      ELSE
         DO 20 I=1,40
            IF(LINE(ITYP1:ITYP2).EQ.DUNIT(I)) THEN
               JUNIT(I)=IU
               FILSTAT='OLD   '
C       CHECK FOR BINARY FILE
               IF(I.EQ.3.OR.I.EQ.5.OR.I.EQ.7) THEN
                 FILSTAT='UNKNOWN'
                 FMTARG=FORM
                 ACCARG=ACCESS
	         END IF
C       OUTPUT FILES MAY NOT BE 'OLD'
               IF(I.EQ.2.OR.I.EQ.4.OR.I.EQ.6
cgzh mbrpt
     &         .OR.I.EQ.26.OR.I.EQ.27
cgzh zinn%
     &         .OR.I.EQ.19.OR.I.EQ.20.OR.I.EQ.21.OR.I.EQ.22) THEN
                 FILSTAT='UNKNOWN'
                 FMTARG='FORMATTED'
	         END IF
C       CHECK FOR PRTP FILE, NO FILE READ NEEDED SO SKIP OUT
               IF(I.EQ.24) THEN
                 GO TO 10
	         END IF
               GO TO 30
            END IF
20       CONTINUE
         WRITE(IOUTS,21) LINE(ITYP1:ITYP2)
21       FORMAT(1X,'ILLEGAL FILE TYPE IN NAME FILE: ',A)
         STOP
30       CONTINUE
      END IF
C
C5------DETERMINE FILE NAME AND WRITE THE FILE NAME IF THE FILE IS NOT 
C5------THE LISTING FILE.  THEN OPEN THE FILE.
      CALL URWORD(LINE,LLOC,INAM1,INAM2,0,N,R,IOUTS,INUNIT)
      INQUIRE(UNIT=IU,OPENED=LOP)
      IF (LOP) CLOSE (UNIT=IU)
C-----IF FILE STATUS IS AMBIGUOUS CHECK FOR "REPLACE" OR "OLD" OPTION
      IF (FILSTAT.EQ.'UNKNOWN') THEN
        CALL URWORD(LINE,LLOC,IOPT1,IOPT2,1,N,R,IOUT,INUNIT)
        IF (LINE(IOPT1:IOPT2).EQ.'REPLACE' .OR.
     &      LINE(IOPT1:IOPT2).EQ.'OLD')
     &      FILSTAT = LINE(IOPT1:IOPT2)
      ENDIF
      IF (FILACT.EQ.' ') FILACT=ACTION(2)

      IF(NFILE.NE.0) WRITE(IOUTS,36) LINE(INAM1:INAM2),
     1     LINE(ITYP1:ITYP2),IU,FILSTAT,FMTARG,ACCARG
36    FORMAT(1X,/1X,'OPENING ',A,/
     &  1X,'FILE TYPE:',A,'   UNIT',I4,3X,'STATUS:',A,/
     &  1X,'FORMAT:',A,3X,'ACCESS:',A)
      OPEN(UNIT=IU,FILE=LINE(INAM1:INAM2),FORM=FMTARG,
     1         ACCESS=ACCARG,STATUS=FILSTAT,ACTION=FILACT)
C
C6------IF THE OPENED FILE IS THE LISTING FILE, WRITE ITS NAME.
C6------GO BACK AND READ NEXT RECORD.
C1------IDENTIFY GWT PACKAGE.
      NFILE=NFILE+1
      IF(NFILE.EQ.1) WRITE(IOUTS,1)
    1 FORMAT(1H ,//
     *'              U.S. GEOLOGICAL SURVEY'/
     *'       Ground-Water Transport Process (GWT)'/
     *'           GWT (Version 3.0)  5/13/2016'///)
      IF(NFILE.EQ.1) WRITE(IOUTS,37) LINE(INAM1:INAM2),IU
37    FORMAT(1X,'LISTING FILE: ',A,3X,'UNIT',I4)
      GO TO 10
C
C7------END OF NAME FILE.  RETURN PROVIDED THAT LISTING FILE AND GWT
C7------FILES HAVE BEEN OPENED.
1000  IF(ILIST.EQ.0) THEN
         WRITE(*,*) ' GWT NAME FILE IS EMPTY.'
         STOP
      ELSE IF(INMOC.EQ.0) THEN
         WRITE(*,*) ' MAIN GWT INPUT FILE HAS NOT BEEN OPENED.'
         STOP
      END IF
C8------RETURN PROVIDED THAT ELLAM NOT BEING USED WITH DP,DK,AGE
      IF(MOCTYPE.EQ.3) THEN
cea*******************************************************************
C        IF(JUNIT(9).GT.0) THEN
C          WRITE(IOUTS,*) ' ELLAM AND AGE PACKAGE ARE NOT COMPATIBLE.'
C          STOP
C        ENDIF
cea*******************************************************************
        IF(JUNIT(10).GT.0) THEN
          WRITE(IOUTS,*) ' ELLAM AND DP PACKAGE ARE NOT COMPATIBLE.'
          STOP
        ENDIF
        IF(JUNIT(11).GT.0) THEN
          WRITE(IOUTS,*) ' ELLAM AND DK PACKAGE ARE NOT COMPATIBLE.'
          STOP
        ENDIF
      END IF
CMOCWT
C CHECK MOCWT OPTIONS FOR COMPATIBILITY
      IF ((PTWTON.EQ.0).AND.(JUNIT(28).GT.0)) THEN
        WRITE(IOUTS,*) '***WARNING*** VBAL PACKAGE ONLY COMPATIBLE ',
     *  'WITH MOCWT OR MOCWTI OPTIONS.'
        WRITE(IOUTS,*) 'VBAL PACKAGE WILL BE IGNORED.'
      ENDIF
C IPDL and IPDA only with MOCWT
      IF((PTWTON.EQ.0).AND.(JUNIT(13).GT.0.OR.JUNIT(14).GT.0)) THEN
        WRITE(IOUTS,*) ' IPDL AND IPDA PACKAGES ONLY COMPATIBLE WITH',
     *  ' MOCWT OR MOCWTI OPTIONS.'
        STOP
      ENDIF
      IF(JUNIT(13).GT.0.AND.JUNIT(14).GT.0) THEN
        WRITE(IOUTS,*) ' IPDL AND IPDA PACKAGES CANNOT BE USED',
     *    ' SIMULTANEOUSLY.'
        STOP
      ENDIF
      IF(JUNIT(9).GT.0.AND.JUNIT(25).GT.0) THEN
        WRITE(IOUTS,*) ' AGE AND CCBD PACKAGES ARE NOT COMPATIBLE.'
        STOP
      ENDIF
cgzh debug mocwt:more restrictions?  DP?  AGE?
cgzh debug cbdy only with mocwt?
      CLOSE(UNIT=INUNIT)
      RETURN
C
      END
C
C***************************************************************
C
      SUBROUTINE SMOC6INTERP(CONC,NPCELL,IBOUND,
     *    NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS)
C
C     *******************************************************
C     INTERPOLATE CONCS FOR CELLS WITH NPCELL=0
C     *******************************************************
C
      DIMENSION CONC(NSCOL,NSROW,NSLAY),NPCELL(NSCOL,NSROW,NSLAY),
     * IBOUND(NCOL,NROW,NLAY)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     -------------------------------------------------------
      DO 10 KS=1,NSLAY  
	K=KS+ISLAY1-1
	DO 10 IS=1,NSROW
	I=IS+ISROW1-1
	DO 10 JS=1,NSCOL
	J=JS+ISCOL1-1
       IF(IBOUND(J,I,K).GT.0) THEN
        IF(NPCELL(JS,IS,KS).EQ.0) THEN
          NCELLS=0
          CONCSUM=0.0
C ONLY USE CELL IF IF HAS PARTICLES
C CHECK COLUMN NEIGHBORS
          IF(JS.GT.1) THEN
		  IF(NPCELL(JS-1,IS,KS).GT.0) THEN
              NCELLS=NCELLS+1		 
              CONCSUM=CONCSUM+CONC(JS-1,IS,KS)
            END IF		 
          END IF		 
          IF(JS.LT.NSCOL) THEN
		  IF(NPCELL(JS+1,IS,KS).GT.0) THEN
              NCELLS=NCELLS+1		 
              CONCSUM=CONCSUM+CONC(JS+1,IS,KS)
            END IF		 
          END IF		 
C CHECK ROW NEIGHBORS
          IF(IS.GT.1) THEN
		  IF(NPCELL(JS,IS-1,KS).GT.0) THEN
              NCELLS=NCELLS+1		 
              CONCSUM=CONCSUM+CONC(JS,IS-1,KS)
            END IF		 
          END IF		 
          IF(IS.LT.NSROW) THEN
		  IF(NPCELL(JS,IS+1,KS).GT.0) THEN
              NCELLS=NCELLS+1		 
              CONCSUM=CONCSUM+CONC(JS,IS+1,KS)
            END IF		 
          END IF		 
C CHECK LAYER NEIGHBORS
          IF(KS.GT.1) THEN
		  IF(NPCELL(JS,IS,KS-1).GT.0) THEN
              NCELLS=NCELLS+1		 
              CONCSUM=CONCSUM+CONC(JS,IS,KS-1)
            END IF		 
          END IF		 
          IF(KS.LT.NSLAY) THEN
		  IF(NPCELL(JS,IS,KS+1).GT.0) THEN
              NCELLS=NCELLS+1		 
              CONCSUM=CONCSUM+CONC(JS,IS,KS+1)
            END IF		 
          END IF		 
C SET CONCENTRATION= SUM OF CONCENTRATIONS / NUMBER OF CELLS USED
          IF(NCELLS.GT.0) CONC(JS,IS,KS)=CONCSUM/NCELLS
        END IF
       END IF
  10  CONTINUE
C
      RETURN
      END
C
C
C***************************************************************
C
      SUBROUTINE SMOC6C(CONC,THCK,SBVL,SRCDCY,
     *    KSTP,NSTP,KPER,NPER,IMOV,NMOV,TIMV,
     *    DELR,DELC,NCOL,NROW,NLAY,
     *    NSCOL,NSROW,NSLAY,IOUTS,JUNIT,PERTIM,TOTIM,SUMTCH,
     *    NPNTCL,ICONFM,ICONLY,
     2 IDKZO,IDKFO,IDKZS,IDKFS,
     3 IDPZO,IDPFO,INSFRUNIT,INLAKUNIT,INDRT,NIUNIT,
C MOCWT
     4 DECAY,IUNIT,MULTSS,
     * SRCAGE,INMNW)
C
C-----from VERSION 1653 15MAY1987 SBAS1H
C     *******************************************************
C     PRINT AND RECORD CONCS
C     *******************************************************
C
C        SPECIFICATIONS
C     -------------------------------------------------------
      DOUBLE PRECISION DECAY
cgzh debug double sbvl
      DOUBLE PRECISION SBVL
      DIMENSION CONC(NSCOL,NSROW,NSLAY),SBVL(6,NIUNIT),
     *  IUNIT(NIUNIT),JUNIT(NIUNIT),THCK(NSCOL,NSROW,NSLAY),
     *  DELR(NCOL),DELC(NROW)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
cea
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
CMOCWT
      INCLUDE 'ptwt.inc'
C
C     -------------------------------------------------------
C CHECK FLAGS FOR OUTPUT
      IPRNT=0
      COLPRT=0
      IF (NPNTCL.EQ.0.AND.KPER.EQ.NPER.AND.KSTP.EQ.NSTP.
     *AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-2.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-1.AND.IMOV.EQ.NMOV) IPRNT=1
      IF(NPNTCL.GE.1) THEN
         IF(MOD(IMOV,NPNTCL).EQ.0) IPRNT=1
      ENDIF
      IF (KPER.EQ.NPER.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
C INITIAL CONDITION
      IF(IMOV.EQ.0) IPRNT=1
C SKIP IF NO OUTPUT
      IF (IPRNT.EQ.0) RETURN
C FOR EACH LAYER: PRINT CONC IF REQUESTED.
      DO 80 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
C
       IF(JUNIT(2).LE.0.AND.JUNIT(3).LE.0) THEN
C PRINT TO MAIN OUTPUT FILE
         WRITE(IOUTS,49) KS,KK
  49  FORMAT(//,'  CONCENTRATIONS: SUBGRID LAYER ',I3,
     *                       ' = MODFLOW LAYER ',I3)
         IF(ICONFM.LT.0) CALL CULAPRS(CONC(1,1,KKS),
     *     'TRANSPORT TIME INCREMENT ',IMOV,KSTP,KPER,
     *     NSCOL,NSROW,KK,-ICONFM,IOUTS)
         IF(ICONFM.GE.0) CALL CULAPRW(CONC(1,1,KKS),
     *     'TRANSPORT TIME INCREMENT ',IMOV,KSTP,KPER,
     *     NSCOL,NSROW,KK,ICONFM,IOUTS)
       ELSE 
C PRINT TO SEPARATE TEXT FILE
         IF(KS.EQ.1.AND.IMOV.EQ.0) THEN
            IF(ICONFM.GE.0) THEN
		    IF(JUNIT(2).GT.0) WRITE(IOUTS,51) JUNIT(2)
            ELSE
              IF(JUNIT(2).GT.0) WRITE(IOUTS,52) JUNIT(2)
            ENDIF
         ENDIF
  51  FORMAT('CONCENTRATION DATA WILL BE SAVED ON UNIT ',I3,
     *' IN MATRIX FORMAT')
  52  FORMAT('CONCENTRATION DATA WILL BE SAVED ON UNIT ',I3,
     *' IN X,Y,Z,CONC FORMAT')
C
         IF(JUNIT(2).GT.0) THEN
	     IF(ICONFM.GE.0) THEN
C PRINT IN MATRIX FORM
C HEADER LINE
           WRITE(JUNIT(2),55) KS, IMOV, KSTP, KPER, SUMTCH
  55  FORMAT('CONCENTRATIONS AT NODES IN SUBGRID LAYER ',I4,', IMOV='
     * ,I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
C
             DO 60 IS=1,NSROW
               WRITE(JUNIT(2),'(1P10E12.4)') (CONC(JS,IS,KS),JS=1,NSCOL)
  60         CONTINUE
C PRINT AS COL (X), ROW (Y), LAY (Z), CONC (SEPARATE LOOP BELOW)
           ELSEIF(ICONFM.LT.0) THEN
             COLPRT=1
           ENDIF
	   ENDIF
       ENDIF
       IF(JUNIT(3).GT.0) THEN
         IF(KS.EQ.1.AND.IMOV.EQ.0) WRITE(IOUTS,65) JUNIT(3)
C PRINT BINARY FILE
  65  FORMAT('CONCENTRATION DATA WILL BE SAVED ON UNIT ',I3,
     *' IN BINARY FORMAT')
         CALL ULASAV(CONC(1,1,KKS),'CONCENTRATIONS  ',
     *     KSTP,KPER,SUMTCH,TOTIM,NSCOL,NSROW,KK,JUNIT(3))
       ENDIF
  80  CONTINUE
      IF(COLPRT.GT.0) THEN
C PRINT IN X,Y,Z,C FORMAT
C HEADER LINE
         WRITE(JUNIT(2),56) IMOV, KSTP, KPER, SUMTCH
  56  FORMAT('CONCENTRATIONS AT NODES (X,Y,Z,CONC): IMOV='
     * ,I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
        TOTROW=0.0
        DO 61 I=1,NROW
        IS=I-ISROW1+1
        TOTROW=TOTROW+DELC(I)
         TOTCOL=0.0
         DO 61 J=1,NCOL
         JS=J-ISCOL1+1
         TOTCOL=TOTCOL+DELR(J)
          TOTTHK=0.0
          DO 61 K=1,NLAY
          KS=K-ISLAY1+1
            IF(K.GE.ISLAY1.AND.K.LE.ISLAY2.AND.
     *         I.GE.ISROW1.AND.I.LE.ISROW2.AND.
     *         J.GE.ISCOL1.AND.J.LE.ISCOL2) THEN
              TOTTHK=TOTTHK+THCK(JS,IS,KS)
              X=TOTCOL-(0.5*DELR(J))
              Y=TOTROW-(0.5*DELC(I))
              Z=TOTTHK-(0.5*THCK(JS,IS,KS))
              WRITE(JUNIT(2),'(1P4E12.4)') X,Y,Z,CONC(JS,IS,KS)
            END IF
  61    CONTINUE
      ENDIF
C PRINT MASS BALANCE
      IF(IMOV.GT.0) THEN
	  IF(PTWTON.NE.1) THEN
         CALL SMOC6M(SBVL,SRCDCY,TIMV,IOUTS,KPER,NPER,KSTP,NSTP,
     *   IMOV,NMOV,SUMTCH,ICONLY,
     *   JUNIT(11),IDKZO,IDKFO,IDKZS,IDKFS,
     *   JUNIT(9),
     *   JUNIT(10),IDPZO,IDPFO,INSFRUNIT,INLAKUNIT,INDRT,NIUNIT,
     *   SRCAGE,INMNW)      
C
C PRINT MASS BALANCE, WEIGHTED PARTICLES
        ELSE
         CALL SMOC6MWT(SBVL,IUNIT,SRCDCY,TIMV,IOUTS,KPER,NPER,KSTP,NSTP,
     *   IMOV,NMOV,SUMTCH,ICONLY,
     *   JUNIT(11),IDKZO,IDKFO,IDKZS,IDKFS,
     *   JUNIT(9),
     *   JUNIT(10),IDPZO,IDPFO,INSFRUNIT,INLAKUNIT,INDRT,NIUNIT,DECAY,
     *   MULTSS,INMNW,JUNIT(25))
        END IF
      END IF
C
   50 RETURN
      END
C
C***************************************************************
C
      SUBROUTINE SMOC6V(THCK,BUFF,POR,RF,
     *                VC,VR,VL,IBOUND,DELT,KSTP,KPER,NSTP,NPER,
     *                NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,NTIMV,
     *                TIMV,VCMAX,VRMAX,VLMAX,TLMIN,IDIR,TOTIM,
     *                MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *                MAXVLJ,MAXVLI,MAXVLK,ITCD,
     *                NPNTVL,IVELFM,JUNIT,SUMTCH,
     *                DKZO,DKFO,DKZS,
     *                DKFS,IDKZO,IDKFO,IDKZS,IDKFS,
cellam
     *                MOCTYPE,TCMIN,TRMIN,DELCOL,DELROW,NIUNIT)
cellam
C
C-----from VERSION 1653 15MAY1987 SBAS1H
C     *******************************************************
C     PRINT VELOCITIES
C     *******************************************************
C
C        SPECIFICATIONS
C     -------------------------------------------------------
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
cellam
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
      DIMENSION DELCOL(NCOL),DELROW(NROW)
cellam
C
      DIMENSION THCK(NSCOL,NSROW,NSLAY),BUFF(NSCOL,NSROW,NSLAY),
     *   POR(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *   VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1),
     *   IBOUND(NCOL,NROW,NLAY),
     *   JUNIT(NIUNIT)
C
C     -------------------------------------------------------
C     ---CALCULATE MAXIMUM VELOCITIES
C
            CALL VELOMAX(THCK,POR,RF,
     *         VC,VR,VL,IBOUND,
     *         NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,
     *         VCMAX,VRMAX,VLMAX,TLMIN,IDIR,
     *         MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *         MAXVLJ,MAXVLI,MAXVLK,
     *         MOCTYPE,TCMIN,TRMIN,DELCOL,DELROW)
C
C     ---COMPUTE NEXT TIME STEP---
cellam
  115  IF (CELDIS.NE.0.0) THEN
cellam
        IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
          TDELCB=CELDIS*CDEL/VCMAX
          TDELRB=CELDIS*RDEL/VRMAX
        ELSEIF (MOCTYPE.EQ.3) THEN
          TDELCB=CELDIS*TCMIN
          TDELRB=CELDIS*TRMIN
        ENDIF
        TDELLB=CELDIS*TLMIN
        ITCD=0
        IF(TDELRB.LT.TDELCB.AND.TDELRB.LE.TDELLB) ITCD=1
        IF(TDELLB.LT.TDELCB.AND.TDELLB.LT.TDELRB) ITCD=2
        TIMV=AMIN1(TDELCB,TDELRB,TDELLB)
cellam
        IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
          IF(AMAX1(VCMAX,VRMAX).LE.1.0E-20.AND.TLMIN.GE.1.0E+20)
     *      WRITE(IOUTS,570)
        ELSEIF(MOCTYPE.EQ.3) THEN
          IF(TCMIN.GE.1.0E+20.AND.TRMIN.GE.1.0E+20.AND.TLMIN.GE.1.0E+20)
     *      WRITE(IOUTS,570)
        ENDIF
cellam
        NTIMV=DELT/TIMV
        NTIMV=NTIMV+1
C
cellam
       ELSEIF(MOCTYPE.EQ.3) THEN
        TIMV=DELT
        NTIMV=1
       ENDIF
cellam
  570 FORMAT(/5X,47H*** WARNING ***  DECREASE CRITERIA IN VELO SUBR)
C
C
C  CHECK FLAGS FOR PRINTING VELOCITIES
      IPRNT=0
      IF(NPNTVL.EQ.0.AND.KPER.EQ.NPER.AND.KSTP.EQ.NSTP) IPRNT=1
      IF(NPNTVL.EQ.-1.AND.KSTP.EQ.NSTP) IPRNT=1
      IF(NPNTVL.GE.1) THEN
        IF(MOD(KSTP,NPNTVL).EQ.0.OR.
     *   (KPER.EQ.NPER.AND.KSTP.EQ.NSTP)) IPRNT=1
      ENDIF
C SKIP IF NO OUTPUT
      IF (IPRNT.EQ.0) GOTO 950
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        DCINV=1.D0/CDEL
        DRINV=1.D0/RDEL
        ARINV=DCINV*DRINV
      ENDIF
C LOOP OVER ALL THREE DIRECTIONS
      DO 888 IDIR=1,3
C DON'T PRINT TO BOTH MAIN AND SEPARATE OUTPUT FILES
cljk      IF(JUNIT(6).GT.0.OR.JUNIT(7).GT.0) GOTO 24
      IF(JUNIT(6).GT.0) GOTO 24
      IF(IDIR.EQ.1) THEN
         IF((JUNIT(6).LE.0).and.(junit(7).le.0)) WRITE(IOUTS,320)
         IF((JUNIT(6).LE.0).and.(junit(7).le.0)) WRITE(IOUTS,330)
         DO 50 KS=1,NSLAY
         K=KS+ISLAY1-1
CV  CONVERT BLOCK FACE VOLUMETRIC FLUXES TO RETARDED NODE VELOCITIES
cellam
         IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) DRIN2=DRINV/2.D0
cellam
         DO 52 IS=1,NSROW
            I=IS+ISROW1-1
            DO 51 JS=1,NSCOL
               J=JS+ISCOL1-1
               IF(IBOUND(J,I,K).EQ.0) THEN
                  BUFF(JS,IS,KS)=0.D0
               ELSE
cellam
                 IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
                   BUFF(JS,IS,KS)=(VC(JS,IS,KS)+VC(JS+1,IS,KS))*DRIN2/
     *                      (POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
                 ELSEIF(MOCTYPE.EQ.3) THEN
                   RIN2=1/(2.D0*RF(JS,IS,KS))
                   BUFF(JS,IS,KS)=(VC(JS,IS,KS)+VC(JS+1,IS,KS))*RIN2/
cgzh bug fix     *            (RF(JS,IS,KS)*DELROW(I)*POR(JS,IS,KS)*THCK(JS,IS,KS))
     *            (DELROW(I)*POR(JS,IS,KS)*THCK(JS,IS,KS))
                 END IF
cellam
               END IF
 51         CONTINUE
 52      CONTINUE
         KKS=KS
         KK=KS+ISLAY1-1
C OUTPUT TO TEXT FILE
Cljk prevent from printing if printing binary
cljk         IF (JUNIT(6).LE.0) THEN
         IF ((JUNIT(6).LE.0).and.(junit(7).le.0)) THEN
         IF(IVELFM.LT.0) CALL ULAPRS(BUFF(1,1,KKS),'VELOCITY (COL)  ',
     *     KSTP,KPER,NSCOL,NSROW,KK,-IVELFM,IOUTS)
         IF(IVELFM.GE.0) CALL ULAPRW(BUFF(1,1,KKS),'VELOCITY (COL)  ',
     *     KSTP,KPER,NSCOL,NSROW,KK,IVELFM,IOUTS)
         ENDIF
C OUTPUT TO BINARY FILE   
      IF (JUNIT(7).GT.0) CALL ULASAV(BUFF(1,1,KKS),'VELOCITY (COL)  ',
     *  KSTP,KPER,SUMTCH,TOTIM,NSCOL,NSROW,KK,JUNIT(7))
   50    CONTINUE
C
      ELSEIF(IDIR.EQ.2) THEN
         IF((JUNIT(6).LE.0).and.(junit(7).le.0)) WRITE(IOUTS,360)
         IF((JUNIT(6).LE.0).and.(junit(7).le.0)) WRITE(IOUTS,330)
         DO 70 KS=1,NSLAY
C KK IS NUMBER OF FLOW LAYER
         K=KS+ISLAY1-1
         IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) DCIN2=DCINV/2.D0
         DO 72 IS=1,NSROW
            I=IS+ISROW1-1
            DO 71 JS=1,NSCOL
               J=JS+ISCOL1-1
               IF(IBOUND(J,I,K).EQ.0) THEN
                  BUFF(JS,IS,KS)=0.D0
               ELSE
cellam
                 IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
                  BUFF(JS,IS,KS)=(VR(JS,IS,KS)+VR(JS,IS+1,KS))*DCIN2/
     *                      (POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
                 ELSEIF(MOCTYPE.EQ.3) THEN
                  RIN2=1/(2.D0*RF(JS,IS,KS))
                  BUFF(JS,IS,KS)=(VR(JS,IS,KS)+VR(JS,IS+1,KS))*RIN2/
cgzh bug fix     *            (DELCOL(J)*POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
     *            (DELCOL(J)*POR(JS,IS,KS)*THCK(JS,IS,KS))
                 END IF
cellam
               END IF
 71         CONTINUE
 72      CONTINUE
         KKS=KS
         KK=KS+ISLAY1-1
cljk         IF (JUNIT(6).LE.0) THEN
         IF ((JUNIT(6).LE.0).and.(junit(7).le.0)) THEN
         IF(IVELFM.LT.0) CALL ULAPRS(BUFF(1,1,KKS),'VELOCITY (ROW)  ',
     *     KSTP,KPER,NSCOL,NSROW,KK,-IVELFM,IOUTS)
         IF(IVELFM.GE.0) CALL ULAPRW(BUFF(1,1,KKS),'VELOCITY (ROW)  ',
     *     KSTP,KPER,NSCOL,NSROW,KK,IVELFM,IOUTS)
         ENDIF
         IF (JUNIT(7).GT.0) CALL ULASAV(BUFF(1,1,KKS),
     *    'VELOCITY (ROW)  ',KSTP,KPER,SUMTCH,TOTIM,NSCOL,NSROW,KK,
     *    JUNIT(7))
   70    CONTINUE
C
      ELSE
         IF(JUNIT(6).LE.0.and.junit(7).le.0) WRITE(IOUTS,365)
         IF(JUNIT(6).LE.0.and.junit(7).le.0) WRITE(IOUTS,330)
         DO 80 KS=1,NSLAY
         K=KS+ISLAY1-1
         IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) ARIN2=ARINV/2.D0
CV
         DO 82 IS=1,NSROW
            I=IS+ISROW1-1
            DO 81 JS=1,NSCOL
               J=JS+ISCOL1-1
               IF(IBOUND(J,I,K).EQ.0) THEN
                  BUFF(JS,IS,KS)=0.D0
               ELSE
cellam
                 IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
                  BUFF(JS,IS,KS)=(VL(JS,IS,KS)+VL(JS,IS,KS+1))*ARIN2/
     *                           (POR(JS,IS,KS)*RF(JS,IS,KS))
                 ELSEIF(MOCTYPE.EQ.3) THEN
                  RIN2=1/(2.D0*RF(JS,IS,KS))
                  BUFF(JS,IS,KS)=(VL(JS,IS,KS)+VL(JS,IS,KS+1))*RIN2/
cgzh bug fix     *               (DELCOL(J)*DELROW(I)*POR(JS,IS,KS)*RF(JS,IS,KS))
     *               (DELCOL(J)*DELROW(I)*POR(JS,IS,KS))
                 END IF
cellam
               END IF
 81         CONTINUE
 82      CONTINUE
         KKS=KS
         KK=KS+ISLAY1-1
cljk         IF (JUNIT(6).LE.0) THEN
         IF ((JUNIT(6).LE.0).and.(junit(7).le.0)) THEN
         IF(IVELFM.LT.0) CALL ULAPRS(BUFF(1,1,KKS),'VELOCITY (LAYER)',
     *     KSTP,KPER,NSCOL,NSROW,KK,-IVELFM,IOUTS)
         IF(IVELFM.GE.0) CALL ULAPRW(BUFF(1,1,KKS),'VELOCITY (LAYER)',
     *     KSTP,KPER,NSCOL,NSROW,KK,IVELFM,IOUTS)
         ENDIF
         IF (JUNIT(7).GT.0) CALL ULASAV(BUFF(1,1,KKS),
     *     'VELOCITY (LAYER)',
     *     KSTP,KPER,SUMTCH,TOTIM,NSCOL,NSROW,KK,JUNIT(7))
 80      CONTINUE
      ENDIF
C
C TEXT (ASCII) OUTPUT FILE NEEDS RECONFIGURED DO LOOP STRUCTURE
C
 24   CONTINUE
      IF (JUNIT(6).GT.0.AND.IDIR.EQ.3) THEN
         WRITE(IOUTS,1133) JUNIT(6)
 1133 FORMAT (//'VELOCITY DATA BEING WRITTEN IN ASCII FORMAT ON',
     * ' UNIT', I3)
         DO 900 KS=1,NSLAY
C WRITE VC HEADER
          WRITE(JUNIT(6),85) KS, KSTP, KPER, TOTIM
  85  FORMAT('VEL. IN COLUMN DIRECTION (in 10E12.4), SUBGRID LAYER '
     * ,I4, 
     *', KSTP=',I5,', KPER=',I5,', TOTIM=',G15.6)
C
          K=KS+ISLAY1-1
          IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) DRIN2=DRINV/2.D0
CV
          DO 92 IS=1,NSROW
             I=IS+ISROW1-1
           DO 91 JS=1,NSCOL
              J=JS+ISCOL1-1
              IF(IBOUND(J,I,K).EQ.0) THEN
                 BUFF(JS,IS,KS)=0.D0
              ELSE
cellam
                 IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
                   BUFF(JS,IS,KS)=(VC(JS,IS,KS)+VC(JS+1,IS,KS))*DRIN2/
     *                      (POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
                 ELSEIF(MOCTYPE.EQ.3) THEN
                   RIN2=1/(2.D0*RF(JS,IS,KS))
                   BUFF(JS,IS,KS)=(VC(JS,IS,KS)+VC(JS+1,IS,KS))*RIN2/
cgzh bug fix     *            (RF(JS,IS,KS)*DELROW(I)*POR(JS,IS,KS)*THCK(JS,IS,KS))
     *            (DELROW(I)*POR(JS,IS,KS)*THCK(JS,IS,KS))
                 END IF
cellam
              END IF
 91        CONTINUE
           WRITE(JUNIT(6),'(1P10E12.4)') (BUFF(JS,IS,KS),JS=1,NSCOL)
 92       CONTINUE
C WRITE VR HEADER
          WRITE(JUNIT(6),95) KS, KSTP, KPER, TOTIM
  95  FORMAT('VEL. IN    ROW DIRECTION (in 10E12.4), SUBGRID LAYER '
     * ,I4,
     *', KSTP=',I5,', KPER=',I5,', TOTIM=',G15.6)
C
          IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) DCIN2=DCINV/2.D0
CV
          DO 102 IS=1,NSROW
             I=IS+ISROW1-1
           DO 101 JS=1,NSCOL
              J=JS+ISCOL1-1
              IF(IBOUND(J,I,K).EQ.0) THEN
                 BUFF(JS,IS,KS)=0.D0
              ELSE
cellam
                 IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
                  BUFF(JS,IS,KS)=(VR(JS,IS,KS)+VR(JS,IS+1,KS))*DCIN2/
     *                      (POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
                 ELSEIF(MOCTYPE.EQ.3) THEN
                  CIN2=1/(2.D0*RF(JS,IS,KS))
                  BUFF(JS,IS,KS)=(VR(JS,IS,KS)+VR(JS,IS+1,KS))*CIN2/
cgzh bug fix     *            (DELCOL(J)*POR(JS,IS,KS)*THCK(JS,IS,KS)*RF(JS,IS,KS))
     *            (DELCOL(J)*POR(JS,IS,KS)*THCK(JS,IS,KS))
                 END IF
cellam
              END IF
 101       CONTINUE
            WRITE(JUNIT(6),'(1P10E12.4)') (BUFF(JS,IS,KS),JS=1,NSCOL)
 102      CONTINUE
C WRITE VL HEADER
          WRITE(JUNIT(6),105) KS, KSTP, KPER, TOTIM
 105  FORMAT('VEL. IN  LAYER DIRECTION (in 10E12.4), SUBGRID LAYER '
     * ,I4,
     *', KSTP=',I5,', KPER=',I5,', TOTIM=',G15.6)
C
          IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) ARIN2=ARINV/2.D0
CV
          DO 112 IS=1,NSROW
             I=IS+ISROW1-1
           DO 111 JS=1,NSCOL
              J=JS+ISCOL1-1
              IF(IBOUND(J,I,K).EQ.0) THEN
                 BUFF(JS,IS,KS)=0.D0
              ELSE
cellam
                 IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
                  BUFF(JS,IS,KS)=(VL(JS,IS,KS)+VL(JS,IS,KS+1))*ARIN2/
     *                           (POR(JS,IS,KS)*RF(JS,IS,KS))
                 ELSEIF(MOCTYPE.EQ.3) THEN
                  RIN2=1/(2.D0*RF(JS,IS,KS))
                  BUFF(JS,IS,KS)=(VL(JS,IS,KS)+VL(JS,IS,KS+1))*RIN2/
cgzh bug fix     *               (DELCOL(J)*DELROW(I)*POR(JS,IS,KS)*RF(JS,IS,KS))
     *               (DELCOL(J)*DELROW(I)*POR(JS,IS,KS))
                 END IF
cellam
              END IF
 111       CONTINUE
           WRITE(JUNIT(6),'(1P10E12.4)') (BUFF(JS,IS,KS),JS=1,NSCOL)
 112      CONTINUE
 900     CONTINUE
      ENDIF
C ENDIF FOR ASCII OUTPUT
 888  CONTINUE
C
C SKIP TO HERE IF NO OUTPUT
  950 CONTINUE
C     ****************************************************************
C
  320 FORMAT(1H ///,
     * 'EFFECTIVE MEAN SOLUTE VELOCITIES IN COLUMN DIRECTION')
  330 FORMAT(1H ,25X,8HAT NODES/)
  360 FORMAT(1H ///,
     * 'EFFECTIVE MEAN SOLUTE VELOCITIES IN ROW DIRECTION')
  365 FORMAT(1H ///,
     * 'EFFECTIVE MEAN SOLUTE VELOCITIES IN LAYER DIRECTION')
      RETURN
      END
C
C
C***************************************************************
C
      SUBROUTINE SMOC5D(DISPCC,
     *  DISPCR,DISPCL,DISPRR,DISPRC,DISPRL,DISPLL,
     *  DISPLC,DISPLR,
     *  KSTP,NSTP,KPER,NPER,NSCOL,NSROW,NSLAY,
     *  IOUTS,NODISP,IDSPFM,NPNTDL)
C
C-----from VERSION 1653 15MAY1987 SBAS1H
C     *******************************************************
C     PRINT DISPERSION COEFFICIENTS
C     *******************************************************
C
C        SPECIFICATIONS
C     -------------------------------------------------------
C
      DIMENSION
     *  DISPCC(NSCOL,NSROW,NSLAY),DISPCR(NSCOL,NSROW,NSLAY),
     *  DISPCL(NSCOL,NSROW,NSLAY),DISPRR(NSCOL,NSROW,NSLAY),
     *  DISPRC(NSCOL,NSROW,NSLAY),DISPRL(NSCOL,NSROW,NSLAY),
     *  DISPLL(NSCOL,NSROW,NSLAY),DISPLC(NSCOL,NSROW,NSLAY),
     *  DISPLR(NSCOL,NSROW,NSLAY)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ---PRINT DISPERSION EQUATION COEFFICIENTS---
C     -------------------------------------------------------
C  DETERMINE IF AND WHEN DISPERSION COEFFICIENTS SHOULD BE PRINTED
      IF(NPNTDL.EQ.0) RETURN
      IPRNT=0
      IF(NPNTDL.GT.0) THEN
         IF(MOD(KSTP,NPNTDL).EQ.0) IPRNT=1
      END IF
      IF(NPNTDL.EQ.-1.AND.KSTP.EQ.NSTP.AND.KPER.EQ.NPER) IPRNT=1
      IF(NPNTDL.EQ.-2.AND.KSTP.EQ.NSTP) IPRNT=1
      IF(IPRNT.EQ.0) RETURN
C
      WRITE(IOUTS,450)
      WRITE(IOUTS,460)
      DO 80 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      WRITE(IOUTS,7763) KK
 7763 FORMAT(/' FOR FLOW LAYER K=',I3/)
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPCC(1,1,KKS),'DISPCC IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPCC(1,1,KKS),'DISPCC IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
   80 CONTINUE
      WRITE(IOUTS,470)
      DO 90 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      WRITE(IOUTS,7763) KK
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPRR(1,1,KKS),'DISPRR IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPRR(1,1,KKS),'DISPRR IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
   90 CONTINUE
      WRITE(IOUTS,491)
      DO 100 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      KKP1=KK+1
      WRITE(IOUTS,7764) KK,KKP1
 7764 FORMAT(/' BETWEEN FLOW LAYERS K=',I3,'   AND K+1=',I3/)
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPLL(1,1,KKS),'DISPLL IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPLL(1,1,KKS),'DISPLL IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  100 CONTINUE
C  SKIP REST OF PRINTING IF NO DISPERSION
      IF(NODISP.EQ.1) GO TO 170
      WRITE(IOUTS,480)
      DO 110 KS=1,NSLAY
      KK=KS
      KK=KS+ISLAY1-1
      WRITE(IOUTS,7763) KK
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPCR(1,1,KKS),'DISPCR IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPCR(1,1,KKS),'DISPCR IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  110 CONTINUE
      WRITE(IOUTS,490)
      DO 120 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      WRITE(IOUTS,7763) KK
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPRC(1,1,KKS),'DISPRC IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPRC(1,1,KKS),'DISPRC IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  120 CONTINUE
      WRITE(IOUTS,492)
      DO 130 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      WRITE(IOUTS,7763) KK
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPCL(1,1,KKS),'DISPCL IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPCL(1,1,KKS),'DISPCL IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  130 CONTINUE
      WRITE(IOUTS,493)
      DO 140 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      WRITE(IOUTS,7763) KK
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPRL(1,1,KKS),'DISPRL IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPRL(1,1,KKS),'DISPRL IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  140 CONTINUE
      WRITE(IOUTS,494)
      DO 150 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      KKP1=KK+1
      WRITE(IOUTS,7764) KK,KKP1
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPLC(1,1,KKS),'DISPLC IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPLC(1,1,KKS),'DISPLC IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  150 CONTINUE
      WRITE(IOUTS,495)
      DO 160 KS=1,NSLAY
      KKS=KS
      KK=KS+ISLAY1-1
      KKP1=KK+1
      WRITE(IOUTS,7764) KK,KKP1
      IF(IDSPFM.LT.0)
     *CALL ULAPRS(DISPLR(1,1,KKS),'DISPLR IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,-IDSPFM,IOUTS)
      IF(IDSPFM.GE.0)
     *CALL ULAPRW(DISPLR(1,1,KKS),'DISPLR IN DSP1FM',
     *            KSTP,KPER,NSCOL,NSROW,KK,IDSPFM,IOUTS)
  160 CONTINUE
C     ****************************************************************
  170 CONTINUE
      RETURN
C     ****************************************************************
C
C
C
  330 FORMAT(1H ,25X,8HAT NODES/)
C      ACTIVATE NEXT LINE IF BOUNDARY VELOCITIES PRINTED
  350 FORMAT(1H ,1P10E12.3)
  450 FORMAT(1H ///,32HDISPERSION EQUATION COEFFICIENTS,10X,33H=(D-IJ)*(
     1B)*(POROS)/(GRID FACTOR))
  460 FORMAT(/35X,14HCC COEFFICIENT/)
  470 FORMAT(/35X,14HRR COEFFICIENT/)
  480 FORMAT(/35X,14HCR COEFFICIENT/)
  490 FORMAT(/35X,14HRC COEFFICIENT/)
  491 FORMAT(/35X,14HLL COEFFICIENT/)
  492 FORMAT(/35X,14HCL COEFFICIENT/)
  493 FORMAT(/35X,14HRL COEFFICIENT/)
  494 FORMAT(/35X,14HLC COEFFICIENT/)
  495 FORMAT(/35X,14HLR COEFFICIENT/)
      END
C
C*************************************************************************
C SMOC5P   DUMP PARTICLE LOCATIONS FOR PLOTTING
C     *************************************************************
      SUBROUTINE SMOC5P(PC,PR,PL,PCONC,
     *                  IOUTPD,KPER,NPER,KSTP,NSTP,IMOV,NMOV,TIMV,NP,
cgzh debug ptwt to pt files
     *                  SUMTCH,NPNTPL,IFLAG,PTWT,
cgzh wtfac
     *                  NSCOL,NSROW,NSLAY,WTFAC,IPRTFM,IOUTS)
C     ***************************************************************
      ALLOCATABLE PLtemp(:)
      DIMENSION PC(NP),PR(NP),PL(NP),PCONC(NP)
cgzh debug ptwt to pt files
      DOUBLE PRECISION PTWT,dist
      DIMENSION PTWT(NP)
      DIMENSION WTFAC(NSCOL,NSROW,NSLAY)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      ALLOCATE(PLtemp(NP))
C     ***************************************************************
C CHECK FLAGS FOR OUTPUT    
      IPRNT=0
      IF (NPNTPL.EQ.0.AND.KPER.EQ.NPER.AND.KSTP.EQ.NSTP.
     *AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTPL.EQ.-2.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTPL.EQ.-1.AND.IMOV.EQ.NMOV) IPRNT=1
      IF(NPNTPL.GE.1) THEN
        IF(MOD(IMOV,NPNTPL).EQ.0.OR.
     *   (KPER.EQ.NPER.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV)) IPRNT=1
      ENDIF
C INITIAL CONDITION        
      IF(IMOV.EQ.0) IPRNT=1
C SKIP IF NO OUTPUT
      IF (IPRNT.EQ.0) THEN
	  DEALLOCATE(PLtemp)
	  RETURN
      END IF
C
C  set pltemp(ip)=pl(ip)
      pltemp=pl
cgzh wtfac
C ALTER PL COORDINATE IF WTFAC OPTION ON
      IF(IPRTFM.GT.0) THEN
   20 CONTINUE
       DO 534 ip=1,np
       OLDC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
       IF(OLDC.LE.0.0D0) GO TO 534
C     ***************************************************************
C           ---COMPUTE OLD LOCATION---
       J=INT(OLDC+0.5D0)
       JS=J-ISCOL1+1
C  IORIG SET TO 1 FOR PARTICLES ORIGINATING IN THIS CELL
C   ORIGINAL PARTICLES HAVE NEGATIVE R COORDINATE
       OLDR=PR(IP)
       IF(OLDR.LT.0.0D0) THEN
         OLDR=-OLDR
       END IF
       I=INT(OLDR+0.5D0)
       IS=I-ISROW1+1
       OLDL=PL(IP)
       K=INT(OLDL+0.5D0)
       KS=K-ISLAY1+1
       IF(JS.LT.1.OR.JS.GT.NSCOL.OR.IS.LT.1.OR.IS.GT.NSROW.OR.
     *  KS.LT.1.OR.KS.GT.NSLAY) THEN
         WRITE(IOUTS,*) ' IP,JS,IS,KS=',IP,JS,IS,KS
         WRITE(IOUTS,*) ' OLDC,OLDR,OLDL=',OLDC,OLDR,OLDL
         WRITE(IOUTS,*) ' PARTICLE ERROR MOVE, STOPPING'
         STOP ' PARTICLE ERROR IN SMOC5P'
       END IF
C  WT ADJUSTMENT
C  THIS REDUCES K LOCATION OF PARTICLES BASED ON WTFAC, THE PERCENTAGE
C    SATURATED THICKNESS TO CELL DIMENSION
C  CALCULATION PARTS:
C    cellface : forward cell face location
C    ptht: particle's "height" relative to lower cell face   
C    *WTFAC   : CALCULATE FRACTION OF DISTANCE PARTICLE SHOULD MOVE
      IF(WTFAC(JS,IS,KS).GT.0.001) THEN
cgzh debug
c      write(iouts,*) 'adjustment made j,k,wtfac',j,k,wtfac(js,is,ks)
c	write(iouts,*) 'old pl',pl(ip)
c
        cellface=K+0.5
        ptht=cellface-pl(ip)
        dist=ptht*(1.0-WTFAC(JS,IS,KS))
c        dist=(1.0-(PL(IP)-0.5-(K-1)))*(1.0-WTFAC(JS,IS,KS))
c	write(iouts,*) 'cellface,ptht,dist',dist
        PLtemp(IP)=PL(IP)+dist
c
c	write(iouts,*) 'new pl',pltemp(ip)
	END IF
  534  CONTINUE 
	END IF
      IF(IFLAG.EQ.0) THEN
cgzh debug
c      npi=0
c      do 9876 ip=1,np
c      if (pl(ip).ge.10.5.and.pl(ip).le.25.5.
c     *and.pc(ip).ge.0.5.and.pc(ip).le.150.5.and.
c     *   abs(pr(ip)).ge.48.5.and.abs(pr(ip)).le.50.5) then
c        npi=npi+1
c      end if
 9876  continue
c       WRITE(IOUTPD,11) KPER,KSTP,IMOV,NPI,TIMV,SUMTCH 
       WRITE(IOUTPD,11) KPER,KSTP,IMOV,NP,TIMV,SUMTCH 

cgzh orig line, uncomment me! (not "NPI")      WRITE(IOUTPD,11) KPER,KSTP,IMOV,NP,TIMV,SUMTCH 
cgzh output ptwt       WRITE(IOUTPD,12) (PC(IP),ABS(PR(IP)),PL(IP),PTWT(IP),IP=1,NP)
c       WRITE(IOUTPD,12) (PC(IP),ABS(PR(IP)),PL(IP),PTWT(IP),IP=1,NP)
cgzh before wtfac change       WRITE(IOUTPD,12) (PC(IP),ABS(PR(IP)),PL(IP),PCONC(IP),IP=1,NP)
      WRITE(IOUTPD,12) (PC(IP),ABS(PR(IP)),PLtemp(IP),PCONC(IP),IP=1,NP)
c uncomment me too! orig line        WRITE(IOUTPD,12) (PC(IP),ABS(PR(IP)),PL(IP),PCONC(IP),IP=1,NP)
cgzh debug 
c      do 987 ip=1,np
c      if (pl(ip).ge.10.5.and.pl(ip).le.25.5.
c     *and.pc(ip).ge.0.5.and.pc(ip).le.150.5.and.
c     *   abs(pr(ip)).ge.48.5.and.abs(pr(ip)).le.50.5) then
c       WRITE(IOUTPD,12) PC(IP),ABS(PR(IP)),PL(IP),PCONC(IP)
c      end if
c 987  continue
cgzh debug end if
      ELSE
C PRINT UNFORMATTED BINARY OUTPUT
         WRITE(IOUTPD) IMOV,NP,TIMV,SUMTCH
c       IF(PC(IP).GT.0.4)
         WRITE(IOUTPD) (PC(IP),ABS(PR(IP)),PLtemp(IP),PCONC(IP),IP=1,NP)
cgzh debug output ptwt
c         WRITE(IOUTPD) (PC(IP),ABS(PR(IP)),PL(IP),PTWT(IP),IP=1,NP)
      END IF
  11  FORMAT(4I10,1PE12.4,1PE12.4)
  12  FORMAT(1P4E12.4)
      DEALLOCATE(PLtemp)
      RETURN
      END
C
*DECK R1MACH
      REAL FUNCTION R1MACH (I)
C***BEGIN PROLOGUE  R1MACH
C***PURPOSE  Return floating point machine dependent constants.
C***LIBRARY   SLATEC
C***CATEGORY  R1
C***TYPE      SINGLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Fox, P. A., (Bell Labs)
C           Hall, A. D., (Bell Labs)
C           Schryer, N. L., (Bell Labs)
C***DESCRIPTION
C
C   R1MACH can be used to obtain machine-dependent parameters for the
C   local machine environment.  It is a function subprogram with one
C   (input) argument, and can be referenced as follows:
C
C        A = R1MACH(I)
C
C   where I=1,...,5.  The (output) value of A above is determined by
C   the (input) value of I.  The results for various values of I are
C   discussed below.
C
C   R1MACH(1) = B**(EMIN-1), the smallest positive magnitude.
C   R1MACH(2) = B**EMAX*(1 - B**(-T)), the largest magnitude.
C   R1MACH(3) = B**(-T), the smallest relative spacing.
C   R1MACH(4) = B**(1-T), the largest relative spacing.
C   R1MACH(5) = LOG10(B)
C
C   Assume single precision numbers are represented in the T-digit,
C   base-B form
C
C              sign (B**E)*( (X(1)/B) + ... + (X(T)/B**T) )
C
C   where 0 .LE. X(I) .LT. B for I=1,...,T, 0 .LT. X(1), and
C   EMIN .LE. E .LE. EMAX.
C
C   The values of B, T, EMIN and EMAX are provided in I1MACH as
C   follows:
C   I1MACH(10) = B, the base.
C   I1MACH(11) = T, the number of base-B digits.
C   I1MACH(12) = EMIN, the smallest exponent E.
C   I1MACH(13) = EMAX, the largest exponent E.
C
C   To alter this function for a particular environment, the desired
C   set of DATA statements should be activated by removing the C from
C   column 1.  Also, the values of R1MACH(1) - R1MACH(4) should be
C   checked for consistency with the local operating system.
C
C***REFERENCES  P. A. Fox, A. D. Hall and N. L. Schryer, Framework for
C                 a portable library, ACM Transactions on Mathematical
C                 Software 4, 2 (June 1978), pp. 177-188.
C***ROUTINES CALLED  XERMSG
C***REVISION HISTORY  (YYMMDD)
C   790101  DATE WRITTEN
C   890213  REVISION DATE from Version 3.2
C   891214  Prologue converted to Version 4.0 format.  (BAB)
C   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
C   900618  Added DEC RISC constants.  (WRB)
C   900723  Added IBM RS 6000 constants.  (WRB)
C   910710  Added HP 730 constants.  (SMR)
C   911114  Added Convex IEEE constants.  (WRB)
C   920121  Added SUN -r8 compiler option constants.  (WRB)
C   920229  Added Touchstone Delta i860 constants.  (WRB)
C   920501  Reformatted the REFERENCES section.  (WRB)
C   920625  Added CONVEX -p8 and -pd8 compiler option constants.
C           (BKS, WRB)
C   930201  Added DEC Alpha and SGI constants.  (RWC and WRB)
C***END PROLOGUE  R1MACH
C
      INTEGER SMALL(2)
      INTEGER LARGE(2)
      INTEGER RIGHT(2)
      INTEGER DIVER(2)
      INTEGER LOG10(2)
C
      REAL RMACH(5)
      SAVE RMACH
C
      EQUIVALENCE (RMACH(1),SMALL(1))
      EQUIVALENCE (RMACH(2),LARGE(1))
      EQUIVALENCE (RMACH(3),RIGHT(1))
      EQUIVALENCE (RMACH(4),DIVER(1))
      EQUIVALENCE (RMACH(5),LOG10(1))
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING THE 68020/68881 COMPILER OPTION
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE AMIGA
C     ABSOFT FORTRAN COMPILER USING SOFTWARE FLOATING POINT
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7EFFFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE APOLLO
C
C     DATA SMALL(1) / 16#00800000 /
C     DATA LARGE(1) / 16#7FFFFFFF /
C     DATA RIGHT(1) / 16#33800000 /
C     DATA DIVER(1) / 16#34000000 /
C     DATA LOG10(1) / 16#3E9A209B /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM
C
C     DATA RMACH(1) / Z400800000 /
C     DATA RMACH(2) / Z5FFFFFFFF /
C     DATA RMACH(3) / Z4E9800000 /
C     DATA RMACH(4) / Z4EA800000 /
C     DATA RMACH(5) / Z500E730E8 /
C
C     MACHINE CONSTANTS FOR THE BURROUGHS 5700/6700/7700 SYSTEMS
C
C     DATA RMACH(1) / O1771000000000000 /
C     DATA RMACH(2) / O0777777777777777 /
C     DATA RMACH(3) / O1311000000000000 /
C     DATA RMACH(4) / O1301000000000000 /
C     DATA RMACH(5) / O1157163034761675 /
C
C     MACHINE CONSTANTS FOR THE CDC 170/180 SERIES USING NOS/VE
C
C     DATA RMACH(1) / Z"3001800000000000" /
C     DATA RMACH(2) / Z"4FFEFFFFFFFFFFFE" /
C     DATA RMACH(3) / Z"3FD2800000000000" /
C     DATA RMACH(4) / Z"3FD3800000000000" /
C     DATA RMACH(5) / Z"3FFF9A209A84FBCF" /
C
C     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
C
C     DATA RMACH(1) / 00564000000000000000B /
C     DATA RMACH(2) / 37767777777777777776B /
C     DATA RMACH(3) / 16414000000000000000B /
C     DATA RMACH(4) / 16424000000000000000B /
C     DATA RMACH(5) / 17164642023241175720B /
C
C     MACHINE CONSTANTS FOR THE CELERITY C1260
C
C     DATA SMALL(1) / Z'00800000' /
C     DATA LARGE(1) / Z'7F7FFFFF' /
C     DATA RIGHT(1) / Z'33800000' /
C     DATA DIVER(1) / Z'34000000' /
C     DATA LOG10(1) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fn COMPILER OPTION
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7FFFFFFF' /
C     DATA RMACH(3) / Z'34800000' /
C     DATA RMACH(4) / Z'35000000' /
C     DATA RMACH(5) / Z'3F9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -fi COMPILER OPTION
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE CONVEX
C     USING THE -p8 OR -pd8 COMPILER OPTION
C
C     DATA RMACH(1) / Z'0010000000000000' /
C     DATA RMACH(2) / Z'7FFFFFFFFFFFFFFF' /
C     DATA RMACH(3) / Z'3CC0000000000000' /
C     DATA RMACH(4) / Z'3CD0000000000000' /
C     DATA RMACH(5) / Z'3FF34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE CRAY
C
C     DATA RMACH(1) / 200034000000000000000B /
C     DATA RMACH(2) / 577767777777777777776B /
C     DATA RMACH(3) / 377224000000000000000B /
C     DATA RMACH(4) / 377234000000000000000B /
C     DATA RMACH(5) / 377774642023241175720B /
C
C     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200
C     NOTE - IT MAY BE APPROPRIATE TO INCLUDE THE FOLLOWING CARD -
C     STATIC RMACH(5)
C
C     DATA SMALL /    20K,       0 /
C     DATA LARGE / 77777K, 177777K /
C     DATA RIGHT / 35420K,       0 /
C     DATA DIVER / 36020K,       0 /
C     DATA LOG10 / 40423K,  42023K /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING G_FLOAT
C
C     DATA RMACH(1) / '00000080'X /
C     DATA RMACH(2) / 'FFFF7FFF'X /
C     DATA RMACH(3) / '00003480'X /
C     DATA RMACH(4) / '00003500'X /
C     DATA RMACH(5) / '209B3F9A'X /
C
C     MACHINE CONSTANTS FOR THE DEC ALPHA
C     USING IEEE_FLOAT
C
C     DATA RMACH(1) / '00800000'X /
C     DATA RMACH(2) / '7F7FFFFF'X /
C     DATA RMACH(3) / '33800000'X /
C     DATA RMACH(4) / '34000000'X /
C     DATA RMACH(5) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE DEC RISC
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE DEC VAX
C     (EXPRESSED IN INTEGER AND HEXADECIMAL)
C     THE HEX FORMAT BELOW MAY NOT BE SUITABLE FOR UNIX SYSTEMS
C     THE INTEGER FORMAT SHOULD BE OK FOR UNIX SYSTEMS
C
C     DATA SMALL(1) /       128 /
C     DATA LARGE(1) /    -32769 /
C     DATA RIGHT(1) /     13440 /
C     DATA DIVER(1) /     13568 /
C     DATA LOG10(1) / 547045274 /
C
C     DATA SMALL(1) / Z00000080 /
C     DATA LARGE(1) / ZFFFF7FFF /
C     DATA RIGHT(1) / Z00003480 /
C     DATA DIVER(1) / Z00003500 /
C     DATA LOG10(1) / Z209B3F9A /
C
C     MACHINE CONSTANTS FOR THE ELXSI 6400
C     (ASSUMING REAL*4 IS THE DEFAULT REAL)
C
C     DATA SMALL(1) / '00800000'X /
C     DATA LARGE(1) / '7F7FFFFF'X /
C     DATA RIGHT(1) / '33800000'X /
C     DATA DIVER(1) / '34000000'X /
C     DATA LOG10(1) / '3E9A209B'X /
C
C     MACHINE CONSTANTS FOR THE HARRIS 220
C
C     DATA SMALL(1), SMALL(2) / '20000000, '00000201 /
C     DATA LARGE(1), LARGE(2) / '37777777, '00000177 /
C     DATA RIGHT(1), RIGHT(2) / '20000000, '00000352 /
C     DATA DIVER(1), DIVER(2) / '20000000, '00000353 /
C     DATA LOG10(1), LOG10(2) / '23210115, '00000377 /
C
C     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000 SERIES
C
C     DATA RMACH(1) / O402400000000 /
C     DATA RMACH(2) / O376777777777 /
C     DATA RMACH(3) / O714400000000 /
C     DATA RMACH(4) / O716400000000 /
C     DATA RMACH(5) / O776464202324 /
C
C     MACHINE CONSTANTS FOR THE HP 730
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     3 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 2100
C     4 WORD DOUBLE PRECISION WITH FTN4
C
C     DATA SMALL(1), SMALL(2) / 40000B,       1 /
C     DATA LARGE(1), LARGE(2) / 77777B, 177776B /
C     DATA RIGHT(1), RIGHT(2) / 40000B,    325B /
C     DATA DIVER(1), DIVER(2) / 40000B,    327B /
C     DATA LOG10(1), LOG10(2) / 46420B,  46777B /
C
C     MACHINE CONSTANTS FOR THE HP 9000
C
C     DATA SMALL(1) / 00004000000B /
C     DATA LARGE(1) / 17677777777B /
C     DATA RIGHT(1) / 06340000000B /
C     DATA DIVER(1) / 06400000000B /
C     DATA LOG10(1) / 07646420233B /
C
C     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
C     THE XEROX SIGMA 5/7/9, THE SEL SYSTEMS 85/86  AND
C     THE PERKIN ELMER (INTERDATA) 7/32.
C
C     DATA RMACH(1) / Z00100000 /
C     DATA RMACH(2) / Z7FFFFFFF /
C     DATA RMACH(3) / Z3B100000 /
C     DATA RMACH(4) / Z3C100000 /
C     DATA RMACH(5) / Z41134413 /
C
C     MACHINE CONSTANTS FOR THE IBM PC
C
C     DATA SMALL(1) / 1.18E-38      /
C     DATA LARGE(1) / 3.40E+38      /
C     DATA RIGHT(1) / 0.595E-07     /
C     DATA DIVER(1) / 1.19E-07      /
C     DATA LOG10(1) / 0.30102999566 /
C
C     MACHINE CONSTANTS FOR THE IBM RS 6000
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE INTEL i860
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE PDP-10 (KA OR KI PROCESSOR)
C
C     DATA RMACH(1) / "000400000000 /
C     DATA RMACH(2) / "377777777777 /
C     DATA RMACH(3) / "146400000000 /
C     DATA RMACH(4) / "147400000000 /
C     DATA RMACH(5) / "177464202324 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     32-BIT INTEGERS (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1) /    8388608 /
C     DATA LARGE(1) / 2147483647 /
C     DATA RIGHT(1) /  880803840 /
C     DATA DIVER(1) /  889192448 /
C     DATA LOG10(1) / 1067065499 /
C
C     DATA RMACH(1) / O00040000000 /
C     DATA RMACH(2) / O17777777777 /
C     DATA RMACH(3) / O06440000000 /
C     DATA RMACH(4) / O06500000000 /
C     DATA RMACH(5) / O07746420233 /
C
C     MACHINE CONSTANTS FOR PDP-11 FORTRAN SUPPORTING
C     16-BIT INTEGERS  (EXPRESSED IN INTEGER AND OCTAL).
C
C     DATA SMALL(1), SMALL(2) /   128,     0 /
C     DATA LARGE(1), LARGE(2) / 32767,    -1 /
C     DATA RIGHT(1), RIGHT(2) / 13440,     0 /
C     DATA DIVER(1), DIVER(2) / 13568,     0 /
C     DATA LOG10(1), LOG10(2) / 16282,  8347 /
C
C     DATA SMALL(1), SMALL(2) / O000200, O000000 /
C     DATA LARGE(1), LARGE(2) / O077777, O177777 /
C     DATA RIGHT(1), RIGHT(2) / O032200, O000000 /
C     DATA DIVER(1), DIVER(2) / O032400, O000000 /
C     DATA LOG10(1), LOG10(2) / O037632, O020233 /
C
C     MACHINE CONSTANTS FOR THE SILICON GRAPHICS
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C
C     DATA RMACH(1) / Z'00800000' /
C     DATA RMACH(2) / Z'7F7FFFFF' /
C     DATA RMACH(3) / Z'33800000' /
C     DATA RMACH(4) / Z'34000000' /
C     DATA RMACH(5) / Z'3E9A209B' /
C
C     MACHINE CONSTANTS FOR THE SUN
C     USING THE -r8 COMPILER OPTION
C
C     DATA RMACH(1) / Z'0010000000000000' /
C     DATA RMACH(2) / Z'7FEFFFFFFFFFFFFF' /
C     DATA RMACH(3) / Z'3CA0000000000000' /
C     DATA RMACH(4) / Z'3CB0000000000000' /
C     DATA RMACH(5) / Z'3FD34413509F79FF' /
C
C     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES
C
C     DATA RMACH(1) / O000400000000 /
C     DATA RMACH(2) / O377777777777 /
C     DATA RMACH(3) / O146400000000 /
C     DATA RMACH(4) / O147400000000 /
C     DATA RMACH(5) / O177464202324 /
C
C     MACHINE CONSTANTS FOR THE Z80 MICROPROCESSOR
C
C     DATA SMALL(1), SMALL(2) /     0,    256/
C     DATA LARGE(1), LARGE(2) /    -1,   -129/
C     DATA RIGHT(1), RIGHT(2) /     0,  26880/
C     DATA DIVER(1), DIVER(2) /     0,  27136/
C     DATA LOG10(1), LOG10(2) /  8347,  32538/
C
C***FIRST EXECUTABLE STATEMENT  R1MACH
      IF (I .LT. 1 .OR. I .GT. 5) CALL XERMSG ('SLATEC', 'R1MACH',
     +   'I OUT OF BOUNDS', 1, 2)
C
      R1MACH = RMACH(I)
      RETURN
C
      END
C
        SUBROUTINE EVALTF(PC,PR,PL,TFVAL,IBOUND,
     *           NUMTF,JS,IS,KS,J,I,K,NCOL,NROW,NLAY,
     *           NSCOL,NSROW,NSLAY,NNZER,NSCOORD,THCK,
     *           XFOR,XBAC,YFOR,YBAC)
C***********************************************************************
C  EVALTF EVALUATES ANY NONZERO TEST FUNCTIONS
C  AT (J+XHAT,I+YHAT,K+ZHAT)=(PC,PR,PL), WHERE 
C  -0.5 <= XHAT <= 0.5,
C  -0.5 <= YHAT <= 0.5,
C  -0.5 <= ZHAT <= 0.5.
C
C  TEST FUNCTION IS DETERMINED BY NSC, NSR AND NSL.  FOR INTERIOR CELLS
C  SURROUNDED BY ACTIVE NEIGHBORS ON UNIFORM GRID, TEST FUNCTION IS 
C  PRODUCT  
C  F(XHAT)*G(YHAT)*H(ZHAT) WHERE
C  F(XHAT)= 0                                      XHAT<= -0.5-1/NSC
C           NSC*XHAT + 0.5*(NSC+1)    -0.5-1/NSC<  XHAT<  -0.5+1/NSC
C           1                         -0.5+1/NSC<= XHAT<=  0.5-1/NSC
C          -NSC*XHAT + 0.5*(NSC+1)     0.5-1/NSC<  XHAT<   0.5+1/NSC
C           0                          0.5+1/NSC<= XHAT 
C  AND SIMILARLY FOR G AND H.
C
C  IF GRID IS NONUNIFORM, TEST FUNCTIONS ARE SCALED BY GRID RATIOS:
C  F(XHAT)=
C           0                          XHAT<= -0.5-(1/NSC)*DELCOL(JS-1)
C
C           1-F   FOR CELL AT JS-1   
C           -0.5-(1/NSC)*DELCOL(JS-1)< XHAT<  -0.5
C
C           2*NSC*(1-XBAC(JS))*XHAT + XBAC(JS) + NSC*(1-XBAC(JS))   
C                                -0.5<=XHAT< -0.5+(1/NSC)*DELCOL(JS)
C           1  0.5+(1/NSC)*DELCOL(JS)<=XHAT<= 0.5-(1/NSC)*DELCOL(JS)
C
C           2*NSC*(XFOR(JS)-1)*XHAT + XFOR(JS) + NSC*(1-XFOR(JS))
C              0.5-(1/NSC)*DELCOL(JS)< XHAT< 0.5
C
C           1-F   FOR CELL AT JS+1
C                                  0.5<=XHAT< 0.5+(1/NSC)*DELCOL(JS+1)
C           0  .5+(1/NSC)*DELCOL(JS+1)<=XHAT     
C
C  AND SIMILARLY FOR G AND H.
C
C  THE SUM OF THE VALUES OF THE TEST FUNCTIONS AT ANY POINT IN THE
C  TRANSPORT DOMAIN IS 1.
C
C  TEST FUNCTIONS EXTEND FROM CENTER OF BOUNDARY CELLS TO BOUNDARY AT
C  VALUE OF ONE.
C  ALL TEST FUNCTIONS ARE ZERO IN INACTIVE CELLS.
C
C  THERE IS NO TEST FUNCTION ASSOCIATED WITH AN INACTIVE CELL.  IN 
C  A CELL ADJACENT TO AN INACTIVE CELL, WHERE THE TEST FUNCTION, IF THE
C  CELL WERE ACTIVE, WOULD BE NOZERO, THAT VALUE IS DISTRIBUTED 
C  PROPORTIONALLY TO OTHER TEST FUNCTIONS WHICH ARE NONZERO AT THAT 
C  POINT.
C
C  EVALTF RETURNS:
C  NNZER  THE NUMBER OF NONZERO TEST FUNCTIONS
C  FOR EACH NONZERO TEST FUNCTION, N=1 TO NNZER
C  NUMTF(N)   NODE NUMBER ASSOCIATED WITH EACH FUNCTION
C  TFVAL(N)   VALUE FOR EACH TEST FUNCTION
C  NSCOORD(3,N)  TRANSPORT GRID COORDINATES OF EACH NODE, IF NOT TF1=1
C**********************************************************************
C
      DOUBLE PRECISION PRECNO
      parameter (precno=5.d-7)
      DIMENSION IBOUND(NCOL,NROW,NLAY),NUMTF(8),TFVAL(8),NSCOORD(3,8)
      DIMENSION NTEMP(8),TEMP(8),NACT(8),NSTEMP(3,8),
     *    THCK(NSCOL,NSROW,NSLAY),XFOR(NSCOL),XBAC(NSCOL),YFOR(NSROW),
     *    YBAC(NSROW)
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
      XHAT=PC-J
      YHAT=PR-I
      ZHAT=PL-K
c
      if(pr.gt.1d10) print *,pc,pr,pl,js,is,ks,j,i,k
      if(abs(xhat).gt..5)  print *,'bad pc',pc,xhat
      if(abs(yhat).gt..5)  print *,'bad pr',pr,yhat
      if(abs(zhat).gt..5)  print *,'eval bad pl',pl,zhat
      if ((js.gt.nscol).or.(is.gt.nsrow).or.(ks.gt.nslay))        
     *   print *, 'evaltf error - point not insubgrid',
     *             js,xhat,is,yhat,ks,zhat
      if (iscol1+js-1.ne.j .or. isrow1+is-1.ne.i .or. islay1+ks-1.ne.k)
     *   print *,'evaltf bad point',iscol1,isrow1,islay1,js,is,ks,j,i,k,
     *           iscol1+js-1,isrow1+is-1,islay1+ks-1
C
      NODE=(KS-1)*NSROW*NSCOL+(IS-1)*NSCOL+JS
C
C
C  FIRST TEST FUNCTION ALWAYS FOR CURRENT NODE
      NSCOORD(1,1)=JS
      NSCOORD(2,1)=IS
      NSCOORD(3,1)=KS
      NUMTF(1)=NODE
C
C********************************************************************
C
C   POINT IN CELL CENTER
C   POINT IN CELL CENTER
      IF ((ABS(XHAT).LE..5-HCINV).AND.(ABS(YHAT).LE..5-HRINV).AND.
     *            (ABS(ZHAT).LE..5-HBINV)) GOTO 5000
c      IF ((ABS(XHAT).LE..5-CINV).AND.(ABS(YHAT).LE..5-RINV).AND.
c     *            (ABS(ZHAT).LE..5-BINV)) GOTO 5000 
C
C********************************************************************
C
C   HANDLE FACES,  EDGES AND CORNERS
C
C  ASSUMING NODE IS INTERIOR,
C  FIND NEIGHBOR WITH DIFFERENT X VALUE
C                     DIFFERENT X,Y VALUES
C                     DIFFERENT Y VALUE
C  EVALUATE PORTION OF THE PRODUCT TEST FUNCTION IN THE 
C                                                 X, Y OR Z DIRECTION
C  THEN ADJUST LAYER AND DO THE SAME
C
C***********************************************************************
      DO 15 N=2,8
      NTEMP(N)=NODE
C
15    CONTINUE
      IF (XHAT.LE.0) THEN
            NX=-1
            INCRE=-1
            XBNS=NCTF-NCTF*XBAC(JS)
            TCX=MIN(1.D0, XBNS*XHAT+XBAC(JS)+0.5D0*XBNS)
      ELSE
                  NX=1
            INCRE=1
            XFNS=NCTF*XFOR(JS)-NCTF
                  TCX=MIN(1.D0, XFNS*XHAT+XFOR(JS)-0.5D0*XFNS)
      ENDIF
C
      JB=J+NX
      JSB=JS+NX
      NTEMP(2)=NTEMP(2)+INCRE
      NTEMP(3)=NTEMP(3)+INCRE
      NTEMP(6)=NTEMP(6)+INCRE
      NTEMP(7)=NTEMP(7)+INCRE
C
      IF (YHAT.LE.0) THEN
            NY=-1
            INCRE=-NSCOL
            YBNS=NRTF-NRTF*YBAC(IS)
            TRY=MIN(1.D0, YBNS*YHAT+YBAC(IS)+0.5D0*YBNS)
      ELSE            
                  NY=1
            INCRE=NSCOL 
            YFNS=NRTF*YFOR(IS)-NRTF
            TRY=MIN(1.D0, YFNS*YHAT+YFOR(IS)-0.5D0*YFNS)
      ENDIF
C
      IB=I+NY
      ISB=IS+NY
      NTEMP(3)=NTEMP(3)+INCRE
      NTEMP(4)=NTEMP(4)+INCRE
      NTEMP(7)=NTEMP(7)+INCRE
      NTEMP(8)=NTEMP(8)+INCRE
C
      IF (ZHAT.LE.0) THEN
            NZ=-1
            INCRE=-NSCOL*NSROW
            ZBAC=TH(0,0,-1,JS,IS,KS,THCK,NSCOL,NSROW,NSLAY)
            ZBNS=NLTF-NLTF*ZBAC
            TLZ=MIN(1.D0, ZBNS*ZHAT+ZBAC+0.5D0*ZBNS)
      ELSE
            NZ=1
            INCRE=NSCOL*NSROW
            ZFOR=TH(0,0,1,JS,IS,KS,THCK,NSCOL,NSROW,NSLAY)
            ZFNS=NLTF*ZFOR-NLTF
            TLZ=MIN(1.D0, ZFNS*ZHAT+ZFOR-0.5D0*ZFNS)
      ENDIF
C
      KB=K+NZ
      KSB=KS+NZ
      NTEMP(5)=NTEMP(5)+INCRE
      NTEMP(6)=NTEMP(6)+INCRE
      NTEMP(7)=NTEMP(7)+INCRE
      NTEMP(8)=NTEMP(8)+INCRE
C
C   DETERMINE TRANSPORT GRID COORDINATES OF NEIGHBORING
C   NODES (IN CASE EVALTF IS CALLED BY BNDYTF)
C
      DO 85 N=2,3
      N1=N+2
      NSTEMP(1,N)=JSB
      NSTEMP(1,N1)=JS
85    CONTINUE
      NSTEMP(1,6)=JSB
      NSTEMP(1,7)=JSB
      NSTEMP(1,8)=JS
      NSTEMP(2,2)=IS
      DO 90 N=3,4
      N1=N+2
      NSTEMP(2,N)=ISB
      NSTEMP(2,N1)=IS
90    CONTINUE
      NSTEMP(2,7)=ISB
      NSTEMP(2,8)=ISB
      DO 95 N=1,4
      N1=N+4
      NSTEMP(3,N)=KS
      NSTEMP(3,N1)=KSB
95    CONTINUE
C
C   CALCULATE TEST FUNCTION VALUES
C
C   XHAT>=0
C
C      TLZ    = MIN(1.0,(.5*( NLTF ))*(-ABS( ZHAT )+.5)+.5) 
C      TRY    = MIN(1.0,(.5*( NRTF ))*(-ABS( YHAT )+.5)+.5) 
C      TCX    = MIN(1.0,(.5*( NCTF ))*(-ABS( XHAT )+.5)+.5) 
      TCXX   = 1.D0-TCX
      TRYY   = 1.D0-TRY 
      TLZZ   = 1.D0-TLZ 
C
      TFVAL(1)=TCX*          TRY*       TLZ
      TEMP(2)= TCXX*        TRY*       TLZ
      TEMP(3)= TCXX*         TRYY*           TLZ
      TEMP(4)= TCX*            TRYY*           TLZ
      TEMP(5)= TCX*           TRY*       TLZZ
      TEMP(6)= TCXX*        TRY*       TLZZ
      TEMP(7)= TCXX*         TRYY*           TLZZ
      TEMP(8)= TCX*            TRYY*           TLZZ
C
C   DETERMINE IF NEIGHBORS ARE ACTIVE, 
C   ASSUMING NODE IS INTERIOR TO FLOW GRID
C   THIS SHOULD BE A LEGITIMATE PROCEDURE, SINCE STORAGE IS NOT
C   ACCESSED OUTSIDE OF X ARRAY, BUT IF CODE HAS BEEN COMPILED WITH
C   SUBSCRIPT CHECKING, IT WILL BLOW UP HERE
C
c      if (jb.gt.0 .and. jb.le.ncol) NACT(2)=IBOUND(JB,I,K)
c      if (jb.gt.0 .and. jb.le.ncol .and. ib.gt.0
c     *    .and. ib.le.nrow) NACT(3)=IBOUND(JB,IB,K)
c      if (ib.gt.0 .and. ib.le.nrow) NACT(4)=IBOUND(J,IB,K)
c      if (kb.gt.0 .and. kb.le.nlay) NACT(5)=IBOUND(J,I,KB)
c      if (jb.gt.0 .and. jb.le.ncol .and. kb.gt.0 .and. kb.le.nlay)
c     *        NACT(6)=IBOUND(JB,I,KB)
c      if (jb.gt.0 .and. jb.le.ncol .and. ib.gt.0
c     *    .and. ib.le.nrow .and. kb.gt.0 .and. kb.le.nlay)
c     *         NACT(7)=IBOUND(JB,IB,KB)
c      if (ib.gt.0 .and. ib.le.nrow.and. kb.gt.0 .and. kb.le.nlay)
c     *       NACT(8)=IBOUND(J,IB,KB)
C
C
C   DETERMINE IF ON SUBGRID BOUNDARY
C   CORNER, EDGE
C 
      IF ((JB.LT.ISCOL1).OR.(JB.GT.ISCOL2).OR.
     *           (ABS(XHAT).LE..5-CINV)) THEN
         NACT(2)=0
         NACT(3)=0
         NACT(6)=0
         NACT(7)=0
	ELSE
         NACT(2)=IBOUND(JB,I,K)
         if(ib.gt.0.and.ib.le.nrow) NACT(3)=IBOUND(JB,IB,K)
         if(kb.gt.0.and.kb.le.nlay) NACT(6)=IBOUND(JB,I,KB)
         if(ib.gt.0.and.ib.le.nrow.and.
     *      kb.gt.0.and.kb.le.nlay) NACT(7)=IBOUND(JB,IB,KB)
      ENDIF
C 
      IF ((IB.LT.ISROW1).OR.(IB.GT.ISROW2).OR.
     *          (ABS(YHAT).LE..5-RINV)) THEN
         NACT(3)=0
         NACT(4)=0
         NACT(7)=0
         NACT(8)=0
	ELSE
         NACT(4)=IBOUND(J,IB,K)
         if(kb.gt.0.and.kb.le.nlay) NACT(8)=IBOUND(J,IB,KB)
      ENDIF
C
      IF ((KB.LT.ISLAY1).OR.(KB.GT.ISLAY2).OR.
     *           (ABS(ZHAT).LE..5-BINV)) THEN
         NACT(5)=0
         NACT(6)=0
         NACT(7)=0
         NACT(8)=0
	ELSE
         NACT(5)=IBOUND(J,I,KB)
      ENDIF
C
C****************************************************************
C   DETERMINE ACTIVE NEIGHBORING CELLS
C
61000 NNZER=1
      NNAC=8
        sum=tfval(1)
        if ((sum.gt.1+precno).or.(sum.lt.0.0-precno)) then
         print *,'evaltf bad tf'
        endif
      DO 65 N=2,8
          sum=sum+temp(n)
          if ((temp(n).gt.1+precno).or.(temp(n).lt.0.0-precno)) then
          print *,'evaltf bad tf'
        endif
        IF (NACT(N).NE.0) THEN
            NNZER=NNZER+1
            NUMTF(NNZER)=NTEMP(N)
            DO 60 M=1,3
            NSCOORD(M,NNZER)=NSTEMP(M,N)
60          CONTINUE
            TFVAL(NNZER)=TEMP(N)
        ELSE
            NUMTF(NNAC)=NTEMP(N)
            TFVAL(NNAC)=TEMP(N)
            NNAC=NNAC-1
        ENDIF
65    CONTINUE
      if ((sum.gt.1+precno).or.(sum.lt.0.0-precno)) then
         print *,'evaltf tfs ne 1'
        endif
C
C   DISTRIBUTE TFVALUES NNZER+1 TO 8 INTO TF S 1 TO NNZER
C
      IF (NNZER.EQ.8) RETURN
C
      SUMA=0.D0
      qs=0.D0
C
      DO 30 N=1,NNZER
        SUMA=SUMA+TFVAL(N)
30    CONTINUE
C
C   FIND FRACTION OF TOTAL ACTIVE VALUE IS TFVAL(N)
      SINV=1/SUMA
      DO 40 N=1,NNZER-1
         TEMP(N)=TFVAL(N)*SINV
         qs=qs+temp(n)
40    CONTINUE
C
c***********************************************************************
c  debug - delete above refs to qs also
c
      temp(nnzer)=1-qs
C   FIND TOTAL VALUE OF TFUNCTIONS FOR INACTIVE CELLS
      SUMI=0
      DO 45 N=NNZER+1,8
        SUMI=SUMI+TFVAL(N)
45    CONTINUE
      if (abs(suma+sumi-1).gt.precno) then
          print *,'sum tfs ne 1',
     *                           sumi-1,pc,pr,pl,xhat,yhat,zhat
      endif
        val=0.D0
      do 4 m=1,nnzer
4      val=val+sumi*temp(m)
      if (abs(val-sumi).gt.precno) 
     *          print *,'distr-actual inact',val-sumi
c
c**********************************************************************
C   ADD FRACTION OF INACTIVE TOTAL TO ACTIVE TVALUES
      VAL=0.D0
      SUMI=1. - SUMA
      DO 50 N=1,NNZER-1
        TFVAL(N)=TFVAL(N)+SUMI*TEMP(N)
        VAL=VAL+TFVAL(N)
50    CONTINUE
      TFVAL(NNZER)=1-VAL
C
      RETURN
C
C***********************************************************************
C
C  ALL SURROUNDING CELLS INACTIVE 
C  OR (XHAT,YHAT,ZHAT) IN CENTER AREA, TF=1
5000  NNZER=1
      TFVAL(1)=1
      RETURN
C
C
      END
C*************************************************************************
      SUBROUTINE BNDYTF(PC,PR,PL,TFVAL,IBOUND,
     *           NUMTF,JS,IS,KS,J,I,K,NCOL,NROW,NLAY,
     *           NSCOL,NSROW,NSLAY,NODESS,NNZER,LBNDY,NTFACE,
     *           THCK,XFOR,XBAC,YFOR,YBAC)
C*************************************************************************
C
C     BNDYTF EVALUATES ANY NONZERO TEST FUNCTIONS AT
C     THE SPACE-TIME BOUNDARY
C
C     INPUT COORDINATES PC,PR,PL ARE FOR POINT BARELY
C     (SMALL=10-4) OUTSIDE OF SUBGRID BOUNDARY IN ONE
C     OR MORE COORDINATE DIRECTION
C     OR, UNDER UPDATED VERSION OF MOVE, EXACTLY ON SUBGRID BOUNDARY
C     WITH JS,IS,KS,J,I,K OUSIDE OF SUBGRID
C
C     NNZER   THE NUMBER OF NONZERO TEST FUNCTIONS
C     FOR N=1,NNZER
C     NUMTF(N)   BOUNDARY NODE NUMBER ASSOCIATED WITH EACH FUNCTION
C     TFVAL (N)  VALUE OF EACH NONZERO TEST FUNCTION
C
C     ARE RETURNED
C*************************************************************************
C
      PARAMETER (SMALL=1.D-4)
      PARAMETER (PRECNO=5.D-7)
      DIMENSION IBOUND(NCOL,NROW,NLAY),NUMTF(8),TFVAL(8),
     *      NSCOORD(3,8),NUMTFB(8),TFVALB(8),LBNDY(NTFACE),
     *    THCK(NSCOL,NSROW,NSLAY),XFOR(NSCOL),XBAC(NSCOL),YFOR(NSROW),
     *    YBAC(NSROW)
      DIMENSION NF(3),NEAR(3)
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
C
      XHAT=PC-J
      YHAT=PR-I
      ZHAT=PL-K
      NF(1)=0
      NF(2)=0
      NF(3)=0
C
C   CHECK FOR BOUNDARY IN X DIRECTION
C   IF POINT IS ACROSS BOUNDARY, NOTE THAT: NF(1) MEANS ACROSS IN X DIR
C                                            NEAR(1) HOLDS WHICH END OF GRID
C
      IF ((JS.EQ.0).OR.(JS.EQ.NSCOL+1)) THEN
         IF (JS.EQ.0) THEN
            JS=1
            J=ISCOL1
            NEAR(1)=1
         ELSE
            JS=NSCOL
            J=ISCOL2
            NEAR(1)=0
         ENDIF
         XHAT=-XHAT
         PC=J+XHAT
         NF(1)=1
      ENDIF
C
C   NOW CHECK IF ACROSS Y BOUNDARY
C
      IF ((IS.EQ.0).OR.(IS.EQ.NSROW+1)) THEN
         IF (IS.EQ.0) THEN
            IS=1
            I=ISROW1
            NEAR(2)=1
         ELSE
            IS=NSROW
            I=ISROW2
            NEAR(2)=0
         ENDIF
         YHAT=-YHAT
         PR=I+YHAT
         NF(2)=1
      ENDIF
C
C   CHECK IF ACROSS Z BOUNDARY
C
      IF ((KS.EQ.0).OR.(KS.EQ.NSLAY+1)) THEN
         IF (KS.EQ.0) THEN
            KS=1
            K=ISLAY1
            NEAR(3)=1
         ELSE
            KS=NSLAY
            K=ISLAY2
            NEAR(3)=0
         ENDIF
         ZHAT=-ZHAT
         PL=K+ZHAT
         NF(3)=1
      ENDIF
C
C   POINT IS NOW IN TRANSPORT SUBGRID
C   IF IT'S PLACED IN AN ACTIVE CELL, GO FIND ASSOCIATED FACES
C
      IF (IBOUND(J,I,K).NE.0) GOTO 300
C
C   IF CELL IS INACTIVE, POINT IS JUST ACROSS AN INTERIOR CELL FACE
C   MOVE IT TO THE OTHER SIDE
C
      IF ((NF(1).EQ.0).AND.(ABS(XHAT).GE..5-SMALL-PRECNO)) THEN
         IF (XHAT.GT.0) THEN
            JS=JS+1
            J=J+1
         ELSE
            JS=JS-1
            J=J-1
         ENDIF
         XHAT=-XHAT
         PC=J+XHAT
      ENDIF
C
      IF (IBOUND(J,I,K).NE.0) GOTO 300
C
C   IF CELL IS INACTIVE, POINT IS JUST ACROSS AN INTERIOR CELL FACE
C   MOVE IT TO THE OTHER SIDE
C
      IF ((NF(2).EQ.0).AND.(ABS(YHAT).GE..5-SMALL-PRECNO)) THEN
         IF (YHAT.GT.0) THEN
            IS=IS+1
            I=I+1
         ELSE
            IS=IS-1
            I=I-1
         ENDIF
         YHAT=-YHAT
         PR=J+YHAT
      ENDIF
C
      IF (IBOUND(J,I,K).NE.0) GOTO 300
C
C   IF CELL IS INACTIVE, POINT IS JUST ACROSS AN INTERIOR CELL FACE
C   MOVE IT TO THE OTHER SIDE
C
      if ((nf(3).ne.0).or.(abs(zhat).lt..5-small-precno))
     *             print *,'ERROR BNDYTF'
      IF (ZHAT.GT.0) THEN
         KS=KS+1
         K=K+1
      ELSE
         KS=KS-1
         K=K-1
      ENDIF
      ZHAT=-ZHAT
      PL=K+ZHAT
      if (ibound(j,i,k).eq.0) print *,'error bndytf- point moved
     *         to inactive cell'
C
C   DISTRIBUTE WEIGHT TO SPACE TIME BOUNDARY NODES
C   USE EVALTF TO FIND SPACIAL BOUNDARY FACES AND DISTRIBUTE 
C   WEIGHT SPACIALLY
C
300      CALL EVALTF(PC,PR,PL,TFVALB,IBOUND,
     *           NUMTFB,JS,IS,KS,J,I,K,NCOL,NROW,NLAY,
     *           NSCOL,NSROW,NSLAY,NNZER,NSCOORD,THCK,
     *           XFOR,XBAC,YFOR,YBAC)
C
      D=NF(1)+NF(2)+NF(3)
      FRAC=1/D
      NPT2=0
C
      DO 50 NPT1=1,NNZER
      IF (NF(1).NE.0) THEN
         NPT2=NPT2+1
         NUMTF(NPT2)=NODESS+NBJ(NEAR(1),NSCOORD(2,NPT1),
     *                  NSCOORD(3,NPT1),NSROW)
         TFVAL(NPT2)=FRAC*TFVALB(NPT1)
      ENDIF
      IF (NF(2).NE.0) THEN
         NPT2=NPT2+1
         NUMTF(NPT2)=NODESS+NBI(NEAR(2),NSCOORD(1,NPT1),
     *                  NSCOORD(3,NPT1),NSCOL,NSROW,NSLAY)
         TFVAL(NPT2)=FRAC*TFVALB(NPT1)
      ENDIF
      IF (NF(3).NE.0) THEN
         NPT2=NPT2+1
         NUMTF(NPT2)=NODESS+NBK(NEAR(3),NSCOORD(1,NPT1),
     *                  NSCOORD(2,NPT1),NSCOL,NSROW,NSLAY)
         TFVAL(NPT2)=FRAC*TFVALB(NPT1)
      ENDIF
50      CONTINUE
      IF (D.EQ.3) TFVAL(3)=1-TFVAL(1)-TFVAL(2)
      NNZER=NPT2
C
      DO 55 N=1,NNZER
      IF (LBNDY(NUMTF(N)-NODESS).NE.0) THEN
c        print *,
c     *    'warning-bndytf outflow distributed to no- or inflow bndy',
c     *    ' bndy nbr=',numtf(n)-nodess,' lbndy=',LBNDY(numtf(n)-nodess)
        TFVAL(1)=TFVAL(1)+TFVAL(N)
      ENDIF
55    CONTINUE  
C
      RETURN
      END
C
C
C  SMOC6I   COPY IBOUND ARRAY
C*************************************************************************
C
      SUBROUTINE GWT1IBOU6(IBOUND,IBOUNC,NODES)
C
C*************************************************************************
C
      DIMENSION IBOUND(NODES),IBOUNC(NODES) 
C
      DO 10 N=1,NODES
          IBOUNC(N)=IBOUND(N)
10    CONTINUE
C
      RETURN
      END
C
C
C
      SUBROUTINE CULAPRS(BUF,TEXT,IMOV,KSTP,KPER,
     *                   NCOL,NROW,ILAY,IPRN,IOUT)
C
C     MODIFIED FROM
C-----VERSION 0755 01NOV1995 ULAPRS
C     BY GZHORNBE ON 07FEB2002
C     ******************************************************************
C     PRINT A 1 LAYER ARRAY IN STRIPS
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*25 TEXT
      DIMENSION BUF(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------MAKE SURE THE FORMAT CODE (IP OR IPRN) IS BETWEEN 1
C1------AND 21.
      IP=IPRN
      IF(IP.LT.1 .OR. IP.GT.21) IP=12
C
C2------DETERMINE THE NUMBER OF VALUES (NCAP) PRINTED ON ONE LINE.
      NCAP=10
      IF(IP.EQ.1) NCAP=11
      IF(IP.EQ.2) NCAP=9
      IF(IP.GT.2 .AND. IP.LT.7) NCAP=15
      IF(IP.GT.6 .AND. IP.LT.12) NCAP=20
      IF(IP.EQ.19) NCAP=5
      IF(IP.EQ.20) NCAP=6
      IF(IP.EQ.21) NCAP=7
C
C3------CALCULATE THE NUMBER OF STRIPS (NSTRIP).
      NCPF=129/NCAP
      IF(IP.GE.13 .AND. IP.LE.18) NCPF=7
      IF(IP.EQ.19) NCPF=13
      IF(IP.EQ.20) NCPF=12
      IF(IP.EQ.21) NCPF=10
      ISP=0
      IF(NCAP.GT.12 .OR. IP.GE.13) ISP=3
      NSTRIP=(NCOL-1)/NCAP + 1
      J1=1-NCAP
      J2=0
C
C4------LOOP THROUGH THE STRIPS.
      DO 2000 N=1,NSTRIP
C
C5------CALCULATE THE FIRST(J1) & THE LAST(J2) COLUMNS FOR THIS STRIP
      J1=J1+NCAP
      J2=J2+NCAP
      IF(J2.GT.NCOL) J2=NCOL
C
C6-------PRINT TITLE ON EACH STRIP DEPENDING ON ILAY
      IF(ILAY.GT.0) THEN
         WRITE(IOUT,1) TEXT,IMOV,KSTP,KPER
    1    FORMAT(/2X,'AT END OF ',A,I6,' IN TIME STEP',I4,
     1     ' IN STRESS PERIOD',I3/2X,78('-'))
      ELSE IF(ILAY.LT.0) THEN
         WRITE(IOUT,2) TEXT,KSTP,KPER
    2    FORMAT('1',/2X,A,' FOR CROSS SECTION AT END OF TIME STEP',I3,
     1     ' IN STRESS PERIOD',I3/2X,77('-'))
      END IF
C
C7------PRINT COLUMN NUMBERS ABOVE THE STRIP
      CALL UCOLNO(J1,J2,ISP,NCAP,NCPF,IOUT)
C
C8------LOOP THROUGH THE ROWS PRINTING COLS J1 THRU J2 WITH FORMAT IP
      DO 1000 I=1,NROW
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     1      180,190,200,210), IP
C
C------------FORMAT 11G10.3
   10 WRITE(IOUT,11) I,(BUF(J,I),J=J1,J2)
   11 FORMAT(1X,I3,2X,1PG10.3,10(1X,G10.3))
      GO TO 1000
C
C------------FORMAT 9G13.6
   20 WRITE(IOUT,21) I,(BUF(J,I),J=J1,J2)
   21 FORMAT(1X,I3,2X,1PG13.6,8(1X,G13.6))
      GO TO 1000
C
C------------FORMAT 15F7.1
   30 WRITE(IOUT,31) I,(BUF(J,I),J=J1,J2)
   31 FORMAT(1X,I3,1X,15(1X,F7.1))
      GO TO 1000
C
C------------FORMAT 15F7.2
   40 WRITE(IOUT,41) I,(BUF(J,I),J=J1,J2)
   41 FORMAT(1X,I3,1X,15(1X,F7.2))
      GO TO 1000
C
C------------FORMAT 15F7.3
   50 WRITE(IOUT,51) I,(BUF(J,I),J=J1,J2)
   51 FORMAT(1X,I3,1X,15(1X,F7.3))
      GO TO 1000
C
C------------FORMAT 15F7.4
   60 WRITE(IOUT,61) I,(BUF(J,I),J=J1,J2)
   61 FORMAT(1X,I3,1X,15(1X,F7.4))
      GO TO 1000
C
C------------FORMAT 20F5.0
   70 WRITE(IOUT,71) I,(BUF(J,I),J=J1,J2)
   71 FORMAT(1X,I3,1X,20(1X,F5.0))
      GO TO 1000
C
C------------FORMAT 20F5.1
   80 WRITE(IOUT,81) I,(BUF(J,I),J=J1,J2)
   81 FORMAT(1X,I3,1X,20(1X,F5.1))
      GO TO 1000
C
C------------FORMAT 20F5.2
   90 WRITE(IOUT,91) I,(BUF(J,I),J=J1,J2)
   91 FORMAT(1X,I3,1X,20(1X,F5.2))
      GO TO 1000
C
C------------FORMAT 20F5.3
  100 WRITE(IOUT,101) I,(BUF(J,I),J=J1,J2)
  101 FORMAT(1X,I3,1X,20(1X,F5.3))
      GO TO 1000
C
C------------FORMAT 20F5.4
  110 WRITE(IOUT,111) I,(BUF(J,I),J=J1,J2)
  111 FORMAT(1X,I3,1X,20(1X,F5.4))
      GO TO 1000
C
C------------FORMAT 10G11.4
  120 WRITE(IOUT,121) I,(BUF(J,I),J=J1,J2)
  121 FORMAT(1X,I3,2X,1PG11.4,9(1X,G11.4))
      GO TO 1000
C
C------------FORMAT 10F6.0
  130 WRITE(IOUT,131) I,(BUF(J,I),J=J1,J2)
  131 FORMAT(1X,I3,1X,10(1X,F6.0))
      GO TO 1000
C
C------------FORMAT 10F6.1
  140 WRITE(IOUT,141) I,(BUF(J,I),J=J1,J2)
  141 FORMAT(1X,I3,1X,10(1X,F6.1))
      GO TO 1000
C
C------------FORMAT 10F6.2
  150 WRITE(IOUT,151) I,(BUF(J,I),J=J1,J2)
  151 FORMAT(1X,I3,1X,10(1X,F6.2))
      GO TO 1000
C
C------------FORMAT 10F6.3
  160 WRITE(IOUT,161) I,(BUF(J,I),J=J1,J2)
  161 FORMAT(1X,I3,1X,10(1X,F6.3))
      GO TO 1000
C
C------------FORMAT 10F6.4
  170 WRITE(IOUT,171) I,(BUF(J,I),J=J1,J2)
  171 FORMAT(1X,I3,1X,10(1X,F6.4))
      GO TO 1000
C
C------------FORMAT 10F6.5
  180 WRITE(IOUT,181) I,(BUF(J,I),J=J1,J2)
  181 FORMAT(1X,I3,1X,10(1X,F6.5))
C
C------------FORMAT 5G12.5
  190 WRITE(IOUT,191) I,(BUF(J,I),J=J1,J2)
  191 FORMAT(1X,I3,1X,1PG12.5,4(1X,G12.5))
C
C------------FORMAT 6G11.4
  200 WRITE(IOUT,201) I,(BUF(J,I),J=J1,J2)
  201 FORMAT(1X,I3,1X,1PG11.4,5(1X,G11.4))
C
C------------FORMAT 7G9.2
  210 WRITE(IOUT,211) I,(BUF(J,I),J=J1,J2)
  211 FORMAT(1X,I3,1X,1PG9.2,6(1X,G9.2))
C
 1000 CONTINUE
 2000 CONTINUE
C
C9------RETURN
      RETURN
      END
C
C
      SUBROUTINE CULAPRW(BUF,TEXT,IMOV,KSTP,KPER,
     *                   NCOL,NROW,ILAY,IPRN,IOUT)
C
C     MODIFIED FROM
C-----VERSION 0758 01NOV1995 ULAPRW
C     BY GZHORNBE ON 07FEB2002
C     ******************************************************************
C     PRINT 1 LAYER ARRAY
C     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      CHARACTER*25 TEXT
      DIMENSION BUF(NCOL,NROW)
C     ------------------------------------------------------------------
C
C1------PRINT A HEADER DEPENDING ON ILAY
      IF(ILAY.GT.0) THEN
         WRITE(IOUT,1) TEXT,IMOV,KSTP,KPER
    1    FORMAT(/2X,'AT END OF ',A,I6,' IN TIME STEP',I4,
     1     ' IN STRESS PERIOD',I3/2X,78('-'))
      ELSE IF(ILAY.LT.0) THEN
         WRITE(IOUT,2) TEXT,KSTP,KPER
    2    FORMAT('1',/2X,A,' FOR CROSS SECTION AT END OF TIME STEP',I3,
     1     ' IN STRESS PERIOD',I3/2X,77('-'))
      END IF
C
C2------MAKE SURE THE FORMAT CODE (IP OR IPRN) IS
C2------BETWEEN 1 AND 21.
    5 IP=IPRN
      IF(IP.LT.1 .OR. IP.GT.21) IP=12
C
C3------CALL THE UTILITY MODULE UCOLNO TO PRINT COLUMN NUMBERS.
      IF(IP.EQ.1) CALL UCOLNO(1,NCOL,0,11,11,IOUT)
      IF(IP.EQ.2) CALL UCOLNO(1,NCOL,0,9,14,IOUT)
      IF(IP.GE.3 .AND. IP.LE.6) CALL UCOLNO(1,NCOL,3,15,8,IOUT)
      IF(IP.GE.7 .AND. IP.LE.11) CALL UCOLNO(1,NCOL,3,20,6,IOUT)
      IF(IP.EQ.12) CALL UCOLNO(1,NCOL,0,10,12,IOUT)
      IF(IP.GE.13 .AND. IP.LE.18) CALL UCOLNO(1,NCOL,3,10,7,IOUT)
      IF(IP.EQ.19) CALL UCOLNO(1,NCOL,0,5,13,IOUT)
      IF(IP.EQ.20) CALL UCOLNO(1,NCOL,0,6,12,IOUT)
      IF(IP.EQ.21) CALL UCOLNO(1,NCOL,0,7,10,IOUT)
C
C4------LOOP THROUGH THE ROWS PRINTING EACH ONE IN ITS ENTIRETY.
      DO 1000 I=1,NROW
      GO TO(10,20,30,40,50,60,70,80,90,100,110,120,130,140,150,160,170,
     1      180,190,200,210), IP
C
C------------ FORMAT 11G10.3
   10 WRITE(IOUT,11) I,(BUF(J,I),J=1,NCOL)
   11 FORMAT(1X,I3,2X,1PG10.3,10(1X,G10.3):/(5X,11(1X,G10.3)))
      GO TO 1000
C
C------------ FORMAT 9G13.6
   20 WRITE(IOUT,21) I,(BUF(J,I),J=1,NCOL)
   21 FORMAT(1X,I3,2X,1PG13.6,8(1X,G13.6):/(5X,9(1X,G13.6)))
      GO TO 1000
C
C------------ FORMAT 15F7.1
   30 WRITE(IOUT,31) I,(BUF(J,I),J=1,NCOL)
   31 FORMAT(1X,I3,1X,15(1X,F7.1):/(5X,15(1X,F7.1)))
      GO TO 1000
C
C------------ FORMAT 15F7.2
   40 WRITE(IOUT,41) I,(BUF(J,I),J=1,NCOL)
   41 FORMAT(1X,I3,1X,15(1X,F7.2):/(5X,15(1X,F7.2)))
      GO TO 1000
C
C------------ FORMAT 15F7.3
   50 WRITE(IOUT,51) I,(BUF(J,I),J=1,NCOL)
   51 FORMAT(1X,I3,1X,15(1X,F7.3):/(5X,15(1X,F7.3)))
      GO TO 1000
C
C------------ FORMAT 15F7.4
   60 WRITE(IOUT,61) I,(BUF(J,I),J=1,NCOL)
   61 FORMAT(1X,I3,1X,15(1X,F7.4):/(5X,15(1X,F7.4)))
      GO TO 1000
C
C------------ FORMAT 20F5.0
   70 WRITE(IOUT,71) I,(BUF(J,I),J=1,NCOL)
   71 FORMAT(1X,I3,1X,20(1X,F5.0):/(5X,20(1X,F5.0)))
      GO TO 1000
C
C------------ FORMAT 20F5.1
   80 WRITE(IOUT,81) I,(BUF(J,I),J=1,NCOL)
   81 FORMAT(1X,I3,1X,20(1X,F5.1):/(5X,20(1X,F5.1)))
      GO TO 1000
C
C------------ FORMAT 20F5.2
   90 WRITE(IOUT,91) I,(BUF(J,I),J=1,NCOL)
   91 FORMAT(1X,I3,1X,20(1X,F5.2):/(5X,20(1X,F5.2)))
      GO TO 1000
C
C------------ FORMAT 20F5.3
  100 WRITE(IOUT,101) I,(BUF(J,I),J=1,NCOL)
  101 FORMAT(1X,I3,1X,20(1X,F5.3):/(5X,20(1X,F5.3)))
      GO TO 1000
C
C------------ FORMAT 20F5.4
  110 WRITE(IOUT,111) I,(BUF(J,I),J=1,NCOL)
  111 FORMAT(1X,I3,1X,20(1X,F5.4):/(5X,20(1X,F5.4)))
      GO TO 1000
C
C------------ FORMAT 10G11.4
  120 WRITE(IOUT,121) I,(BUF(J,I),J=1,NCOL)
  121 FORMAT(1X,I3,2X,1PG11.4,9(1X,G11.4):/(5X,10(1X,G11.4)))
      GO TO 1000
C
C------------ FORMAT 10F6.0
  130 WRITE(IOUT,131) I,(BUF(J,I),J=1,NCOL)
  131 FORMAT(1X,I3,1X,10(1X,F6.0):/(5X,10(1X,F6.0)))
      GO TO 1000
C
C------------ FORMAT 10F6.1
  140 WRITE(IOUT,141) I,(BUF(J,I),J=1,NCOL)
  141 FORMAT(1X,I3,1X,10(1X,F6.1):/(5X,10(1X,F6.1)))
      GO TO 1000
C
C------------ FORMAT 10F6.2
  150 WRITE(IOUT,151) I,(BUF(J,I),J=1,NCOL)
  151 FORMAT(1X,I3,1X,10(1X,F6.2):/(5X,10(1X,F6.2)))
      GO TO 1000
C
C------------ FORMAT 10F6.3
  160 WRITE(IOUT,161) I,(BUF(J,I),J=1,NCOL)
  161 FORMAT(1X,I3,1X,10(1X,F6.3):/(5X,10(1X,F6.3)))
      GO TO 1000
C
C------------ FORMAT 10F6.4
  170 WRITE(IOUT,171) I,(BUF(J,I),J=1,NCOL)
  171 FORMAT(1X,I3,1X,10(1X,F6.4):/(5X,10(1X,F6.4)))
      GO TO 1000
C
C------------ FORMAT 10F6.5
  180 WRITE(IOUT,181) I,(BUF(J,I),J=1,NCOL)
  181 FORMAT(1X,I3,1X,10(1X,F6.5):/(5X,10(1X,F6.5)))
      GO TO 1000
C
C------------FORMAT 5G12.5
  190 WRITE(IOUT,191) I,(BUF(J,I),J=1,NCOL)
  191 FORMAT(1X,I3,2X,1PG12.5,4(1X,G12.5):/(5X,5(1X,G12.5)))
      GO TO 1000
C
C------------FORMAT 6G11.4
  200 WRITE(IOUT,201) I,(BUF(J,I),J=1,NCOL)
  201 FORMAT(1X,I3,2X,1PG11.4,5(1X,G11.4):/(5X,6(1X,G11.4)))
      GO TO 1000
C
C------------FORMAT 7G9.2
  210 WRITE(IOUT,211) I,(BUF(J,I),J=1,NCOL)
  211 FORMAT(1X,I3,2X,1PG9.2,6(1X,G9.2):/(5X,7(1X,G9.2)))
C
 1000 CONTINUE
C
C5------RETURN
      RETURN
      END
