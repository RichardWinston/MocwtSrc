C     Last change: RBW, July 29, 2015
C     Changed to allow for array bounds checking.
      SUBROUTINE CELCONEC(IBOUND,CIN,NPCELL,
     X  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NODESS)
C.....Static initialization for solver
CMOCWT
      INCLUDE 'ptwt.inc'
      INTEGER CIN,IBOUND,NPCELL
      DIMENSION CIN(6,NODESS),IBOUND(NCOL,NROW,NLAY),
     X  NPCELL(NSCOL,NSROW,NSLAY)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C...
C.....Create cell connection list
               NX=NSCOL
               NY=NSROW
               NZ=NSLAY
               NXY=NX*NY
               NXYZ=NSCOL*NSROW*NSLAY
      DO 70 M=1,NXYZ
               DO 63 IC=1,6
   63             CIN(IC,M)=0
C.....Skip excluded cells
               CALL MTOIJK(M,IS,JS,KS,NSCOL,NSROW)
               K=KS+ISLAY1-1
               I=JS+ISROW1-1
               J=IS+ISCOL1-1
               IF(IBOUND(J,I,K).EQ.0) GO TO 70
CMOCWT Skip NPCELL<1 cells for weighted code
               IF(PTWTON.EQ.1.AND.NPCELL(IS,JS,KS).LT.1) GO TO 70
C    
C.....Left
CMOCWT
               IF(IS.GT.1) THEN
CMOCWT Out of the 4 possible combinations of PTWTON and NPCELL,
C      only NPCELL=0 and PTWTON=1 should be excluded.  So, include
C      all times it is NOT this
			   IF(.NOT.(PTWTON.EQ.1.AND.NPCELL(IS-1,JS,KS).LT.1).AND.
     X             IBOUND(J-1,I,K).NE.0) CIN(3,M)=M-1
               END IF
C.....Right
               IF(IS.LT.NX) THEN
			   IF(.NOT.(PTWTON.EQ.1.AND.NPCELL(IS+1,JS,KS).LT.1).AND.
     X             IBOUND(J+1,I,K).NE.0) CIN(4,M)=M+1
               END IF
C.....Front
               IF(JS.GT.1) THEN
			   IF(.NOT.(PTWTON.EQ.1.AND.NPCELL(IS,JS-1,KS).LT.1).AND.
     X             IBOUND(J,I-1,K).NE.0) CIN(2,M)=M-NX
               END IF
C.....Back
               IF(JS.LT.NY) THEN
			   IF(.NOT.(PTWTON.EQ.1.AND.NPCELL(IS,JS+1,KS).LT.1).AND.
     X             IBOUND(J,I+1,K).NE.0) CIN(5,M)=M+NX
               END IF
C.....Below
               IF(KS.GT.1) THEN
			   IF(.NOT.(PTWTON.EQ.1.AND.NPCELL(IS,JS,KS-1).LT.1).AND.
     X             IBOUND(J,I,K-1).NE.0) CIN(1,M)=M-NXY
               END IF
C.....Above
               IF(KS.LT.NZ) THEN
			   IF(.NOT.(PTWTON.EQ.1.AND.NPCELL(IS,JS,KS+1).LT.1).AND.
     X             IBOUND(J,I,K+1).NE.0) CIN(6,M)=M+NXY
               END IF
   70          CONTINUE
      RETURN
	END
C
      SUBROUTINE INIT(IBOUND,CI,CIN,CIR,CIRH,CIRL,IND,MRNO,MRNO0,
     X     EPSSLV,MAXIT,IDIREC,NPCELL,
     X     NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NBN,NRN,NODESS,ISSIZH)
C.....Static initialization for solver
      INTEGER CIN,IBOUND,NPCELL
      DIMENSION CIN(6,NODESS),IBOUND(NCOL,NROW,NLAY),
     *  NPCELL(NSCOL,NSROW,NSLAY)
      INTEGER CI,CIR,CIRH,CIRL,MRNO,MRNO0,IND
c      DOUBLE PRECISION RS,VA,VAD,RHS,AP,PP,RA,RR,SS,XX,WW,WORK
      DIMENSION CI(6,NODESS),CIR(19,ISSIZH),CIRH(10,ISSIZH),
     X     CIRL(10,ISSIZH),MRNO(NODESS),MRNO0(NODESS),IND(NODESS)
      COMMON/CAIS/MAR(6,6),MAR1(19,19)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C...
               NX=NSCOL
               NY=NSROW
               NZ=NSLAY
               NXY=NX*NY
               NXYZ=NSCOL*NSROW*NSLAY
C.....Create cell connection list
      CALL CELCONEC(IBOUND,CIN,NPCELL,
     X  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NXYZ)
C
C.....Establish red-black ordering
        CALL REORDR(IDIREC,CI,CIR,CIRH,CIRL,IND,MRNO,MRNO0,
     X     NBN,NRN,NX,NY,NZ,NXYZ,ISSIZH)
         RETURN
         END
      SUBROUTINE LDCI(CI,MRNO,NX,NY,NZ,NODESS)
C.....Loads the CI array for connected nodes based on the renumbered
C.....     mesh
      INTEGER CI,I,MRNO,NX,NY,NZ,NXY,NXYZ,KX,KY,KZ,M
      DIMENSION CI(6,NODESS),I(6),MRNO(NODESS)
C...
      NXY = NX*NY
      NXYZ = NXY*NZ
      DO 1000 M=1,NXYZ
         CALL MTOIJK(M,KX,KY,KZ,NX,NY)
         MA = MRNO(M)
C.....MA is 0 for excluded cells
         IF(MA.GT.0) THEN
            I(1) = M-NXY
            I(2) = M-NX
            I(3) = M-1
            I(4) = M+1
            I(5) = M+NX
            I(6) = M+NXY
C.....I is natural node number
            DO 10 J=1,6
               IF ( (I(J) .GE. 1) .AND. (I(J) .LE. NXYZ)) 
     X              CI(J,MA) = MRNO(I(J))
   10       CONTINUE
C.....Modify for ends of node rows
            IF (KZ .EQ. 1) CI(1,MA) = 0
            IF (KY .EQ. 1) CI(2,MA) = 0
            IF (KX .EQ. 1) CI(3,MA) = 0
            IF (KX .EQ. NX) CI(4,MA) = 0
            IF (KY .EQ. NY) CI(5,MA) = 0
            IF (KZ .EQ. NZ) CI(6,MA) = 0
         ENDIF
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE LDCIR(CIR,CIRL,CIRH,CI,NBN,NRN,NXYZ)
C.....Loads the CIR, CIRL, CIRH index arrays for the 
C.....     conjugate gradient solver
      INTEGER CI,CIR,CIRL,CIRH,NBN,NRN,NXYZ
      DIMENSION CIR(19,NRN),CI(6,NXYZ),CIRL(10,NRN),CIRH(10,NRN)
c...not sure if this loop 99 is needed
      DO 99 K=1,NBN
         CIRL(1,K)=0
         CIRH(1,K)=0
   99 CONTINUE
      DO 100 K=1,NBN
         KK = NRN+K
C.....Lower level
         IX = CI(1,KK)
         IF (IX .EQ. 0) THEN
            DO 5 J=1,5
               CIR(J,K) = 0
    5       CONTINUE
         ELSE
            DO 10 J=1,5
               IY = CI(J,IX)
               IF (IY .EQ. 0) THEN
                  CIR(J,K) = 0
               ELSE
                  IND = IY-NRN
                  CIR(J,K) = IND
                  IF (IND .LT. K) THEN
                     CIRL(1,K) = CIRL(1,K)+1
                     CIRL(CIRL(1,K)+1,K) = J
                  ENDIF
                  IF (IND .GT. K) THEN
                     CIRH(1,K) = CIRH(1,K)+1
                     CIRH(CIRH(1,K)+1,K) = J
                  ENDIF
               ENDIF
   10       CONTINUE
         ENDIF
C.....Same level
         IX = CI(2,KK)
         IF (IX .EQ. 0) THEN
            DO 15 J=6,8
               CIR(J,K) = 0
   15       CONTINUE
         ELSE
            DO 20 J=2,4
               IY = CI(J,IX)
               IF (IY .EQ. 0) THEN
                  CIR(J+4,K) = 0
               ELSE
                  IND = IY-NRN
                  CIR(J+4,K) = IND
                  IF (IND .LT. K) THEN
                     CIRL(1,K) = CIRL(1,K)+1
                     CIRL(CIRL(1,K)+1,K) = J+4
                  ENDIF
                  IF (IND .GT. K) THEN
                     CIRH(1,K) = CIRH(1,K)+1
                     CIRH(CIRH(1,K)+1,K) = J+4
                  ENDIF
               ENDIF
   20       CONTINUE
         ENDIF
         IX = CI(3,KK)
         IF (IX .EQ. 0) THEN
            CIR(9,K) = 0
         ELSE
            IY = CI(3,IX)
            IF (IY .EQ. 0) THEN
               CIR(9,K) = 0
            ELSE
               IND = IY-NRN
               CIR(9,K) = IND
               IF (IND .LT. K) THEN
                  CIRL(1,K) = CIRL(1,K)+1
                  CIRL(CIRL(1,K)+1,K) = 9
               ENDIF
               IF (IND .GT. K) THEN
                  CIRH(1,K) = CIRH(1,K)+1
                  CIRH(CIRH(1,K)+1,K) = 9
               ENDIF
            ENDIF
         ENDIF
C...
         CIR(10,K) = K
C...
         IX = CI(4,KK)
         IF (IX .EQ. 0) THEN
            CIR(11,K) = 0
         ELSE
            IY = CI(4,IX)
            IF (IY .EQ. 0) THEN
               CIR(11,K) = 0
            ELSE
               IND = IY-NRN
               CIR(11,K) = IND
               IF (IND .LT. K) THEN
                  CIRL(1,K) = CIRL(1,K)+1
                  CIRL(CIRL(1,K)+1,K) = 11
               ENDIF
               IF (IND .GT. K) THEN
                  CIRH(1,K) = CIRH(1,K)+1
                  CIRH(CIRH(1,K)+1,K) = 11
               ENDIF
            ENDIF
         ENDIF
         IX = CI(5,KK)
         IF (IX .EQ. 0) THEN
            DO 25 J=12,14
               CIR(J,K) = 0
   25       CONTINUE
         ELSE
            DO 30 J=3,5
               IY = CI(J,IX)
               IF (IY .EQ. 0) THEN
                  CIR(J+9,K) = 0
               ELSE
                  IND = IY-NRN
                  CIR(J+9,K) = IND
                  IF (IND .LT. K) THEN
                     CIRL(1,K) = CIRL(1,K)+1
                     CIRL(CIRL(1,K)+1,K) = J+9
                  ENDIF
                  IF (IND .GT. K) THEN
                     CIRH(1,K) = CIRH(1,K)+1
                     CIRH(CIRH(1,K)+1,K) = J+9
                  ENDIF
               ENDIF
   30       CONTINUE
         ENDIF
C.....Upper level
         IX = CI(6,KK)
         IF (IX .EQ. 0) THEN
            DO 35 J=15,19
               CIR(J,K) = 0
   35       CONTINUE
         ELSE
            DO 40 J=2,6
               IY = CI(J,IX)
               IF (IY .EQ. 0) THEN
                  CIR(J+13,K) = 0
               ELSE
                  IND = IY-NRN
                  CIR(J+13,K) = IND
                  IF (IND .LT. K) THEN
                     CIRL(1,K) = CIRL(1,K)+1
                     CIRL(CIRL(1,K)+1,K) = J+13
                  ENDIF
                  IF (IND .GT. K) THEN
                     CIRH(1,K) = CIRH(1,K)+1
                     CIRH(CIRH(1,K)+1,K) = J+13
                  ENDIF
               ENDIF
   40       CONTINUE
         ENDIF
  100 CONTINUE
      RETURN
      END
      SUBROUTINE LDIND(IDIREC,IND,NX,NY,NZ,NODESS)
C.....Loads the IND array with a natural numbering in the permutation
C.....     of x, y, and z specified by IDIREC
      INTEGER IND,IDIREC,NX,NY,NZ
      DIMENSION IND(NODESS)
C...
      IF (IDIREC.EQ.1) THEN
         IX = 1
         IY = NX
         IZ = NX*NY
      ELSEIF (IDIREC.EQ.2) THEN
         IX = 1
         IY = NX*NZ
         IZ = NX
      ELSEIF (IDIREC.EQ.3) THEN
         IX = NY
         IY = 1
         IZ = NX*NY
      ELSEIF (IDIREC.EQ.4) THEN
         IX = NZ*NY
         IY = 1
         IZ = NY
      ELSEIF (IDIREC.EQ.5) THEN
         IX = NZ
         IY = NX*NZ
         IZ = 1
      ELSEIF (IDIREC.EQ.6) THEN
         IX = NY*NZ
         IY = NZ
         IZ = 1
      ENDIF
      L = 0
      LK = 1
      DO 30 K=1,NZ
         LJ = LK
         DO 20 J=1,NY
            LI = LJ
            DO 10 I=1,NX
               L = L+1
               IND(LI) = L
               LI = LI+IX
   10       CONTINUE
            LJ = LJ+IY
   20    CONTINUE
         LK = LK+IZ
   30 CONTINUE
      RETURN
      END
C **********************************************************************  B  20
C                                                                         B  30
      BLOCK DATA                                                          B  40
      COMMON /CAIS/ MAR(6,6),MAR1(19,19)                                  B 130
      DATA MAR /10,05,04,03,02,01, 15,10,08,07,06,02, 16,12,10,09,07,03,  B 210
     1          17,13,11,10,08,04, 18,14,13,12,10,05, 19,18,17,16,15,10/  B 220
      DATA MAR1/10, 5, 4, 3, 2, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,0   B 230
     1         ,15,10, 8, 7, 6, 5, 4, 3, 0, 2, 0, 0, 0, 0, 1, 0, 0, 0,0   B 240
     1         ,16,12,10, 9, 7, 0, 5, 0, 4, 3, 0, 2, 0, 0, 0, 1, 0, 0,0   B 250
     1         ,17,13,11,10, 8, 0, 0, 5, 0, 4, 3, 0, 2, 0, 0, 0, 1, 0,0   B 260
     1         ,18,14,13,12,10, 0, 0, 0, 0, 5, 0, 4, 3, 2, 0, 0, 0, 1,0   B 270
     1         , 0,15, 0, 0, 0,10, 8, 7, 0, 6, 0, 0, 0, 0, 2, 0, 0, 0,0   B 280
     1         , 0,16,15, 0, 0,12,10, 9, 8, 7, 0, 6, 0, 0, 3, 2, 0, 0,0   B 290
     1         , 0,17, 0,15, 0,13,11,10, 0, 8, 7, 0, 6, 0, 4, 0, 2, 0,0   B 300
     1         , 0, 0,16, 0, 0, 0,12, 0,10, 9, 0, 7, 0, 0, 0, 3, 0, 0,0   B 310
     1         ,19,18,17,16,15,14,13,12,11,10, 9, 8, 7, 6, 5, 4, 3, 2,1   B 320
     1         , 0, 0, 0,17, 0, 0, 0,13, 0,11,10, 0, 8, 0, 0, 0, 4, 0,0   B 330
     1         , 0, 0,18, 0,16, 0,14, 0,13,12, 0,10, 9, 7, 0, 5, 0, 3,0   B 340
     1         , 0, 0, 0,18,17, 0, 0,14, 0,13,12,11,10, 8, 0, 0, 5, 4,0   B 350
     1         , 0, 0, 0, 0,18, 0, 0, 0, 0,14, 0,13,12,10, 0, 0, 0, 5,0   B 360
     1         , 0,19, 0, 0, 0,18,17,16, 0,15, 0, 0, 0, 0,10, 8, 7, 6,2   B 370
     1         , 0, 0,19, 0, 0, 0,18, 0,17,16, 0,15, 0, 0,12,10, 9, 7,3   B 380
     1         , 0, 0, 0,19, 0, 0, 0,18, 0,17,16, 0,15, 0,13,11,10, 8,4   B 390
     1         , 0, 0, 0, 0,19, 0, 0, 0, 0,18, 0,17,16,15,14,13,12,10,5   B 400
     1        , 0, 0, 0, 0, 0, 0, 0, 0, 0,19, 0, 0, 0, 0,18,17,16,15,10/  B 410
      END                                                                 B 420
C
      SUBROUTINE RBORD(MRNO0,NX,NY,NZ,NODESS)
C.....Maps the natural node number, M, into the red-black node
C.....     number, MRNO0, based on renumbering sweeps in the x, then y,
C.....     then z-directions.
      INTEGER MRNO0,NXYZ,NX,NY,NZ
      DIMENSION MRNO0(NODESS)
C...
      NXY = NX*NY
      NXYZ = NXY*NZ
      NRN = (NXYZ+MOD(NXYZ,2))/2
      DO 10 M=1,NXYZ
C.....Calculate the red-black number for a given natural
C.....     node number,M
      MMOD = MOD(M,2)
      M1 = (M+MMOD)/2
      CALL MTOIJK(M,IX,IY,IZ,NX,NY)
C.....Is sum of indices even or odd?
      IF (MOD(IX+IY+IZ,2) .EQ. 1) THEN
C.....It's a red node
        MRNO0(M) = M1 
      ELSE
C.....It's a black node
        MRNO0(M) = NRN+M1
      ENDIF
   10 CONTINUE 
      RETURN
      END
      SUBROUTINE REORDR(IDIREC,CI,CIR,CIRH,CIRL,IND,MRNO,MRNO0,
     X     NBN,NRN,NX,NY,NZ,NXYZ,ISSIZH)
C.....Performs the red-black reordering of the nodes based on the
C.....     selected permutation of the coordinate directions
      INTEGER NX,NY,NZ,NXYZ,IDIREC,IND,NBN,NRN,MRNO
      INTEGER CI,CIR,CIRL,CIRH
      DIMENSION CIR(19,ISSIZH),CI(6,NXYZ),CIRL(10,ISSIZH),
     *  CIRH(10,ISSIZH)
      DIMENSION MRNO(NXYZ),MRNO0(NXYZ),IND(NXYZ)
C.....Poninters for RCG solver
      COMMON/CAIS/MAR(6,6),MAR1(19,19)
C...
      IF(IDIREC.EQ.1) THEN
         N1 = NX
         N2 = NY
         N3 = NZ
      ELSEIF(IDIREC.EQ.2) THEN
         N1 = NX
         N2 = NZ
         N3 = NY
      ELSEIF(IDIREC.EQ.3) THEN
         N1 = NY
         N2 = NX
         N3 = NZ
      ELSEIF(IDIREC.EQ.4) THEN
         N1 = NY
         N2 = NZ
         N3 = NX
      ELSEIF(IDIREC.EQ.5) THEN
         N1 = NZ
         N2 = NX
         N3 = NY
      ELSEIF(IDIREC.EQ.6) THEN
         N1 = NZ
         N2 = NY
         N3 = NX
      ENDIF
C.....Load MRNO0 using RBORD and N1,N2,N3     
      CALL RBORD(MRNO0,N1,N2,N3,NXYZ)
C.....Now swap indices from MRNO0 to MRNO
      CALL LDIND(IDIREC,IND,NX,NY,NZ,NXYZ)
      DO 50 I=1,NXYZ
         MRNO(IND(I)) = MRNO0(I)
   50 CONTINUE
C.....Load indexing arrays used by the solver, CGRIES 
      NMOD = MOD(NXYZ,2)
      NRN = (NXYZ+NMOD)/2
      NBN = NXYZ-NRN
      CALL LDCI(CI,MRNO,NX,NY,NZ,NXYZ)
      CALL LDCIR(CIR,CIRL,CIRH,CI,NBN,NRN,NXYZ)
      RETURN
      END
      SUBROUTINE ABMULT(VA,VAD,Y,X,CI,NRN,NXYZ,NBN)
C.....Multiplies A_rb*y for black nodes
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION VA,VAD,X,Y,S
      INTEGER CI,NRN
      DIMENSION CI(6,NXYZ)
      DIMENSION X(NRN),Y(NRN),VA(6,NRN),VAD(NXYZ)
C...
      DO 10 IRN=1,NRN
         S = 0.0
         DO 5 JC=1,6
            JCOL = CI(JC,IRN)
            IF (JCOL.GT.0) S = S+VA(JC,IRN)*Y(JCOL-NRN)
    5    CONTINUE
         X(IRN) = S/VAD(IRN)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE ARMULT(VA,Y,X,CI,NBN,NRN,NXYZ)
C.....Multiplies A_br^T*y for red nodes
C.....      A_br^T(=A_rb)
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION VA,X,Y,S
      INTEGER NBN,NRN,CI
      DIMENSION X(NRN),Y(NRN),VA(6,NRN),CI(6,NXYZ)
C...
      DO 10 I=1,NBN
         IB = I+NRN
         S = 0.0
         DO 5 J=1,6
            JR = CI(J,IB)
            IF (JR.GT.0) S = S+VA(7-J,JR)*Y(JR)
    5    CONTINUE
         X(I) = S
   10 CONTINUE
      RETURN
      END
      SUBROUTINE CGRIES(VA,VAD,RHS,AP,PP,RA,RR,SS,XX,WW,WORK,
     X     EPSSLV,MAXIT,CI,CIR,CIRL,CIRH,NXYZ,NRN,NBN,
     X     IOUTS,ITRN,NCXIT)
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION VA,VAD,RHS
c      DOUBLE PRECISION AP,PP,RA,RR,SS,XX,WW,WORK
c      DOUBLE PRECISION ALPHA,DELTA,H,R1,R00,RAT
c      DOUBLE PRECISION XIP
      INTEGER MAXIT,NBN,NRN,NXYZ
      INTEGER CI,CIR,CIRH,CIRL,MAR,MAR1
c      INCLUDE 'param0.inc'
      COMMON/CAIS/MAR(6,6),MAR1(19,19)
      DIMENSION VA(6,NRN),VAD(NXYZ),RHS(0:NXYZ)
      DIMENSION AP(NRN),PP(NRN),RA(10,NRN),RR(NRN),SS(NRN),
     X     XX(NRN),WW(NRN),
     X     WORK(NRN)
      DIMENSION CI(6,NXYZ),CIR(19,NRN),CIRH(10,NRN),CIRL(10,NRN)
C...
      R00 = SQRT(XIP(RHS,RHS,NXYZ))
      DO 10 I=1,NBN
         XX(I) = 0.0
   10 CONTINUE
      IF (R00 .LE. EPSSLV ) THEN
        IF(ITRN.EQ.NCXIT) WRITE(IOUTS,*)
     *   '  Initial residual < EPSSLV;',
     *   ' delta-C = 0 accepted; iterative solver skipped'
        RETURN
      ENDIF
      DO 20 IB=1,NBN
         WORK(IB) = 0.0
         DO 15 J=1,10
            RA(J,IB) = 0.0
   15    CONTINUE
   20 CONTINUE
C.....Multiply the 2 off-diagonal blocks,
C.....     subract from the black diagonal and factor
C.....     the result, RA
      CALL FORMR(VA,VAD,RA,CI,NBN,NRN,NXYZ)
      CALL RFACTM(RA,CIR,CIRH,NBN,WORK,NRN)
C.....Divide the red right hand side, RHS_r, by diagonal, D_r
      DO 700 I=1,NRN
         RHS(I) = RHS(I)/VAD(I)
  700 CONTINUE
C.....Form the RHS of the reduced system
      NRNP1 = NRN+1
      CALL ARMULT(VA,RHS(1),WW,CI,NBN,NRN,NXYZ)
      CALL VPSV(RHS(NRNP1),RHS(NRNP1),WW,-1.0,NBN)
C.....The right hand side is now stored in the bottom half of RHS
C.....Form the initial residual
      CALL ABMULT(VA,VAD,XX,RR,CI,NRN,NXYZ,NBN)
      CALL ARMULT(VA,RR,WW,CI,NBN,NRN,NXYZ)
      CALL DBMULT(VAD,XX,RR,NBN,NRN,NXYZ)
      CALL VPSV(WW,RR,WW,-1.0,NBN)
C.....WW now has R*X_black
      CALL VPSV(RR,RHS(NRNP1),WW,-1.0,NBN)
C.....Form the initial search direction
      ICOUNT = 0
  799 CONTINUE
      CALL LSOLV(RA,RR,WW,CIR,CIRL,NBN,NRN)
      CALL USOLV(RA,WW,PP,CIR,CIRH,NBN,NRN)
      ICOUNT = ICOUNT+1
      IF(ICOUNT.GT.MAXIT) GO TO 999
      CALL ABMULT(VA,VAD,PP,AP,CI,NRN,NXYZ,NBN)
      CALL ARMULT(VA,AP,WW,CI,NBN,NRN,NXYZ)
      CALL DBMULT(VAD,PP,AP,NBN,NRN,NXYZ)
      CALL VPSV(AP,AP,WW,-1.0,NBN)
      DELTA = XIP(PP,AP,NBN)
cgzh debug output
      if(delta.eq.0.0) then
	write(iouts,*) '***ERROR*** Delta=', delta,' in MOCIMP'
      stop 'Error in MOCIMP'
      end if
      ALPHA = XIP(PP,RR,NBN)/DELTA
C.....  Update the solution
      CALL VPSV(XX,XX,PP,ALPHA,NBN)
C.....  Update the residual
      CALL VPSV(RR,RR,AP,-ALPHA,NBN)
      R1 = SQRT(XIP(RR,RR,NBN))
      RAT=R1/R00
C..        WRITE(FUCLOG,*) 'i, current relative residual: ', ICOUNT, RAT
C..        WRITE(FUCLOG,*) 'current residual: ', R1
      IF (RAT.LE.EPSSLV) GOTO 1010
      CALL LSOLV(RA,RR,WW,CIR,CIRL,NBN,NRN)
      CALL USOLV(RA,WW,SS,CIR,CIRH,NBN,NRN)
      H = XIP(SS,AP,NBN)/DELTA
      CALL VPSV(PP,SS,PP,-H,NBN)
 1000 CONTINUE
      GO TO 799
  999 CONTINUE
      WRITE(IOUTS,9002) '****ERROR****   Conjugate-gradient solver ',
     X     'reached maximum iterations: ',MAXIT
      WRITE(IOUTS,9003) ICOUNT, RAT
      WRITE(IOUTS,*)
      WRITE(IOUTS,*) 'If residual is acceptable, increase EPSSLV'
      WRITE(IOUTS,*) 'If residual is unacceptable, increase MAXIT'
      STOP
 9002 FORMAT(/TR10,2A,I5)
 1010 CONTINUE
      IF(ITRN.EQ.NCXIT) WRITE(IOUTS,9003) ICOUNT, RAT
 9003 FORMAT('  No. of solver iterations = ',I4,' Relative residual = ',
     *1pE12.4)
c      IF(ICOUNT.EQ.1)
c     *  WRITE(IOUTS,*) '***WARNING***   Number of iterations = 1;',
c     *  '    Suggest decreasing EPSSLV'
C.....Form the red solution from the black half
      CALL ABMULT(VA,VAD,XX,WW,CI,NRN,NXYZ,NBN)
      CALL VPSV(RHS(1),RHS(1),WW,-1.0,NRN)
C.....Set rhs_black equal to x_black
      CALL VPSV(RHS(NRNP1),XX,WW,0.0,NBN)
      RETURN
      END
      SUBROUTINE DBMULT(VAD,Y,X,NBN,NRN,NXYZ)
C.....Multiplies Diag_A*y for the black nodes
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION VAD,X,Y
      INTEGER NBN,NRN
      DIMENSION X(NRN),Y(NRN),VAD(NXYZ)
C...
      DO 10 I=1,NBN
         J = I+NRN
         X(I) = Y(I)*VAD(J)
   10 CONTINUE
      RETURN
      END
      SUBROUTINE FORMR(VA,VAD,RA,CI,NBN,NRN,NXYZ)
C.....Form the product of the off-diagonal blocks,
C.....Scale and subtract from the black diagonal to get
C.....     the reduced matrix, RA.
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION VA,VAD,RA,DD
      INTEGER NBN,NRN,CI
      DIMENSION VA(6,NRN),VAD(NXYZ),RA(10,NRN),CI(6,NXYZ)
      COMMON/CAIS/MAR(6,6),MAR1(19,19)
C...    
      DO 1 K=1,NBN
         RA(10,K) = VAD(K+NRN)
    1 CONTINUE
      DO 10 K=1,NRN
         DO 5 IC=1,6
            IROW = CI(IC,K)
            IF (IROW.GT.0) THEN
               DD = VA(IC,K)
               NROW = IROW-NRN
               DO 20 J=1,IC
                  IF (CI(J,K).GT.0) 
     X     RA(MAR(IC,J),NROW) = RA(MAR(IC,J),NROW)-DD*VA(J,K)/VAD(K)
   20          CONTINUE
            ENDIF
    5    CONTINUE
   10 CONTINUE
      RETURN
      END
      SUBROUTINE LSOLV(RA,RR,WW,CIR,CIRL,NBN,NRN)
C.....Solves the lower triangular matrix equation, L*ww=rr
C.....     for the black-node equations
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION RA,RR,WW,S
      INTEGER CIRL,CIR
      DIMENSION RR(NRN),WW(NRN),RA(10,NRN),CIRL(10,NRN),CIR(19,NRN)
C...
      WW(1) = RR(1)
      DO 10 K=2,NBN
         S = RR(K)
         DO 5 JJ=1,CIRL(1,K)
            J = CIRL(JJ+1,K)
            JCOL = CIR(J,K)
! RBW Code changed to avoid Fortran trickery            
!            S = S-RA(J,K)*WW(JCOL)
! RBW new code
            if (J.gt.10) then
              if (K.lt.NRN) then
                S = S-RA(J-10,K+1)*WW(JCOL)
              endif
            else
              S = S-RA(J,K)*WW(JCOL)
            endif 
! RBW end new code            
    5    CONTINUE
         WW(K) = S
   10 CONTINUE
      RETURN
      END
      SUBROUTINE RFACTM(RA,CIR,CIRH,NBN,WORK,NRN)
C.....Factors the reduced matrix, RA, to modified incomplete LDL^T factors
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION RA,WORK,DD
      INTEGER CIR,CIRH,NBN
      DIMENSION RA(10,NRN),CIR(19,NRN),CIRH(10,NRN),WORK(NRN)
      COMMON/CAIS/MAR(6,6),MAR1(19,19)
C...
      DO 30 K=1,NBN-1
         RA(10,K) = RA(10,K)+WORK(K)
         DO 20 II=1,CIRH(1,K)
            I = CIRH(II+1,K)
            IROW = CIR(I,K)
! RBW Changed to avoid Fortran trickery   
! original code         
!            DD = RA(MAR1(I,10),IROW)/RA(10,K)
!            RA(MAR1(I,10),IROW) = DD
! RBW new code
            if (MAR1(I,10).GT.10) then
              if (IROW.lt.NRN) then
                DD = RA(MAR1(I,10)-10,IROW+1)/RA(10,K)
                RA(MAR1(I,10)-10,IROW+1) = DD
              endif
            else
              DD = RA(MAR1(I,10),IROW)/RA(10,K)
              RA(MAR1(I,10),IROW) = DD
            endif
! RBW end new code         
            DO 10 JJ=1,II
               J = CIRH(JJ+1,K)
               JROW = CIR(J,K)
               L = MAR1(I,J)
               IF (L .EQ. 0) THEN
! RBW Changed to avoid Fortran trickery   
! original code         
!                  WORK(IROW) = WORK(IROW)-DD*RA(10,K)*RA(20-J,JROW)
!                  RA(10,JROW) = RA(10,JROW)-DD*RA(10,K)*RA(20-J,JROW)
! RBW new code
                  if ((20-J).gt.10) then
                    if (JROW.lt.NRN) then
                      WORK(IROW) = 
     *                  WORK(IROW)-DD*RA(10,K)*RA(20-J-10,JROW+1)
                      RA(10,JROW) = 
     *                  RA(10,JROW)-DD*RA(10,K)*RA(20-J-10,JROW+1)
                    endif 
                  else
                    WORK(IROW) = WORK(IROW)-DD*RA(10,K)*RA(20-J,JROW)
                    RA(10,JROW) = RA(10,JROW)-DD*RA(10,K)*RA(20-J,JROW)
                  endif
! RBW end new code         
                  GO TO 10
               END IF
! RBW Changed to avoid Fortran trickery   
! original code         
!               RA(L,IROW) = RA(L,IROW)-DD*RA(10,K)*RA(20-J,JROW)
! RBW new code
               if (L.gt.10) then
                 if ((20-J).gt.10) then
                   if ((IROW.lt.NRN).and.(JROW.lt.NRN)) then
                     RA(L-10,IROW+1) = 
     *                 RA(L-10,IROW+1)-DD*RA(10,K)*RA(20-J-10,JROW+1)
                   endif 
                 else
                   if (IROW.lt.NRN) then
                     RA(L-10,IROW+1) = 
     *                 RA(L-10,IROW+1)-DD*RA(10,K)*RA(20-J,JROW)
                   endif 
                 endif
               else
                 if ((20-J).gt.10) then
                   if (JROW.lt.NRN) then
                     RA(L,IROW) = 
     *                 RA(L,IROW)-DD*RA(10,K)*RA(20-J-10,JROW+1)
                   endif 
                 else
                   RA(L,IROW) = RA(L,IROW)-DD*RA(10,K)*RA(20-J,JROW)
                 endif
               endif
! RBW end new code         
   10       CONTINUE
   20    CONTINUE
   30 CONTINUE
      RA(10,NBN) = RA(10,NBN)+WORK(NBN)
      RETURN
      END
      SUBROUTINE USOLV(RA,WW,YY,CIR,CIRH,NBN,NRN)
C.....Solves the upper triangular matrix equation, U*yy=ww
C.....     for the black-node equations
c      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
c      DOUBLE PRECISION RA,WW,YY,S
      INTEGER CIRH,CIR
      DIMENSION YY(NRN),WW(NRN),RA(10,NRN),CIRH(10,NRN),CIR(19,NRN)
C...
      YY(NBN) = WW(NBN)/RA(10,NBN)
      DO 1000 K=NBN-1,1,-1
         S = WW(K)/RA(10,K)
         DO 500 JJ=1,CIRH(1,K)
            J = CIRH(JJ+1,K)
            JCOL = CIR(J,K)
C RBW code changed to avoid fortran trickery
C original code            
!            S = S-RA(20-J,JCOL)*YY(JCOL)
! RBW new code
            if ((20-j).gt.10) then
              if (JCOL.lt.NRN) then
                S = S-RA(20-J-10,JCOL+1)*YY(JCOL)
              endif
            else
              S = S-RA(20-J,JCOL)*YY(JCOL)
            endif
! RBW end new code            
  500    CONTINUE 
         YY(K) = S
 1000 CONTINUE
      RETURN
      END
      SUBROUTINE VPSV(X,Y,Z,S,N)
C.....Calculates a vector plus a scalar times a vector
C.....     x=y+s*z
c      DOUBLE PRECISION X,Y,Z,S
      INTEGER N
      DIMENSION X(N),Y(N),Z(N)
C...
      DO 10 I=1,N
        X(I) = Y(I)+S*Z(I)
10    CONTINUE
      RETURN
      END
      FUNCTION XIP(X,Y,N)
C.....Forms the inner or dot product of two vectors
c      DOUBLE PRECISION XIP,X,Y
      INTEGER N
      DIMENSION X(*),Y(*)
C...
      XIP = 0.0
      DO 10 I=1,N
        XIP = XIP+X(I)*Y(I)
10    CONTINUE
      RETURN
      END
