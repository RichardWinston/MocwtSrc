C  IF CHANGES MADE TO MOVE, CHECK IF MOVEBI SHOULD ALSO BE CHANGED
C     ***************************************************************
C
      SUBROUTINE MOVE(PC,PR,PL,PCONC,IPTID,
     *  VC,VR,VL,
     *  RF,THCK,POR,CINFL,
     *  SRCFLO,SNKFLO,CONC,
     *  CAVG,IBOUND,
     *  CNOLD,SUMC,
     *  NPCELL,NPOLD,LIMBO,
     *  PNEWC,PNEWR,PNEWL,IGENPT,IGENLK,
     *  NEWPTS,NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,NPMAX,NLIMBO,
     *  IOUTS,DECAY,IMOV,TIMV,NP,
     *  NCINFL,
     *  IABOVE,IBELOW,
     *  VCMAX,VRMAX,TLMIN,ICONLY,INTRPL,
     *  DKFO,DKFS,INDK,IDKFO,IDKFS,INBFLX,CINFLA,CINFLB,CINXY,
     *  WTFAC,
cgzh ptob
     *  INUNITPTOB,SUMTCH,IPTOBLST,
     *  NUMPTOB,NUMPTOB_MNW,
     *  IPTOBMNW,PTWT,SUMWT)
C
C     ***************************************************************
cgzh ptob
      DIMENSION IPTOBLST(4,NUMPTOB),
     * IPTOBMNW(NSCOL,NSROW,NSLAY),SUMWT(NSCOL,NSROW,NSLAY),PTWT(NPMAX)
      DOUBLE PRECISION SINKTEMP
      ALLOCATABLE SINKTEMP(:,:,:)
C
      DOUBLE PRECISION DECAY
      DOUBLE PRECISION TINY,SMALL,PRECNO
      DOUBLE PRECISION OLDR,OLDC,OLDL,VCP,VRP,VLP,DISTC,DISTR,DISTL
      DOUBLE PRECISION DBLTMP,TESTCK,DVDC,DVDL,DVDR
      DOUBLE PRECISION DCYFCT,DCYT,DCYT2,SUMC
      PARAMETER (PERCNT=0.001)
      PARAMETER (TINY=1.D-20)
      PARAMETER (SMALL=1.D-4)
      PARAMETER (HUGE=1.E20)
C  CONSTANT FOR CHECK OF SINGLE PRECISION INTERVAL
      PARAMETER (PRECNO=5.D-7)
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
     *  SUMC(NSCOL,NSROW,NSLAY),
     *  NPCELL(NSCOL,NSROW,NSLAY),NPOLD(NSCOL,NSROW,NSLAY),
     *  LIMBO(NLIMBO),
     *  PNEWC(NEWPTS),PNEWR(NEWPTS),PNEWL(NEWPTS),
     *  IGENPT(NSCOL,NSROW,NSLAY),IGENLK(NSCOL,NSROW,NSLAY)
      DIMENSION DKFO(NSCOL,NSROW,NSLAY),
     *    DKFS(NSCOL,NSROW,NSLAY)
cgzh cbdy
      DIMENSION CINFLA(NSCOL,NSROW),CINFLB(NSCOL,NSROW),
     * CINXY(NSCOL,NSROW,NSLAY),WTFAC(NSCOL,NSROW,NSLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ***************************************************************
cgzh ptob
      ALLOCATE (SINKTEMP(NSCOL,NSROW,NSLAY))
C3D  ONE MOVE ONLY, RETURN HERE IF GENERATE PARTICLES
C
C  COMPUTE DECAY TERMS
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
   10 NPTM=NP
C
C  CLEAR SUMC
      DO 11 KS=1,NSLAY
      DO 11 IS=1,NSROW
      DO 11 JS=1,NSCOL
      SUMC(JS,IS,KS)=0.0
   11 CONTINUE
C
C  UNCOMMENT FOLLOWING LINES FOR MORE DETAILED OUTPUT OF MOVE LOOP
C     WRITE(IOUTS,*) ' START LOOP OVER PARTICLES'
C     WRITE(IOUTS,*) ' NP',NP
C
C        ---MOVE EACH PARTICLE---
      DO 590 IP=1,NP
      OLDC=PC(IP)
C  SKIP PARTICLE IF C COORDINATE IS ZERO
      IF(OLDC.LE.0.0D0) GO TO 590
C     ***************************************************************
C           ---COMPUTE OLD LOCATION---
      J=INT(OLDC+0.5D0)
      JS=J-ISCOL1+1
C  IORIG SET TO 1 FOR PARTICLES ORIGINATING IN THIS CELL
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
      IF(JS.LT.1.OR.JS.GT.NSCOL.OR.IS.LT.1.OR.IS.GT.NSROW.OR.
     *  KS.LT.1.OR.KS.GT.NSLAY) THEN
         WRITE(IOUTS,*) ' IP,JS,IS,KS=',IP,JS,IS,KS
         WRITE(IOUTS,*) ' OLDC,OLDR,OLDL=',OLDC,OLDR,OLDL
         WRITE(IOUTS,*) ' PARTICLE ERROR MOVE, STOPPING'
         STOP ' PARTICLE ERROR IN MOVE'
      END IF
C  SKIP PARTICLE IF LOCATED IN INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 590
C  ISOURC=0 SIGNIFIES CELL IS NOT A SOURCE (USED TO BE CALLED KFLAG)
      ISOURC=0
C1  USE ANALYTIC EXPRESSION WITHIN BLOCK FOR LINEAR V
      TSTEP=TIMV
C1  SAVE INITIAL LOCATION INDICES
      J1=J
      JS1=JS
      I1=I
      IS1=IS
      K1=K
      KS1=KS
C1  THIS IS BEGINNING OF LOOP FOR PARTIAL TIME STEPS
C1  CONST CONVERTS Q IN REAL UNITS TO RETARDED V IN RELATIVE UNITS
   32 CONST=ARINV/(RF(JS,IS,KS)*THCK(JS,IS,KS)*POR(JS,IS,KS))
C1  FIRST FOR C
      JSV1=JS
      JSV2=JS+1
cgzh intrpl 
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
     *               (VC(JSV2,IS,KS)-VC(JSV2,IS2,KS))
               ELSE
C  CHECK TO SEE IF THE CELL ABOVE AND TO THE RIGHT IS NO-FLOW
C    OR IF THE CELL TO THE RIGHT IS NO-FLOW, DON'T INTERPOLATE
                  IF(THCK(JS+1,IS2,KS).EQ.0.0.OR.
     *               THCK(JS+1,IS,KS).EQ.0.0) THEN
                     VC2=VC(JSV2,IS,KS)
                  ELSE
                     VC2=VC(JSV2,IS2,KS)+(1.0-(FACR-0.5))*
     *                  (VC(JSV2,IS,KS)-VC(JSV2,IS2,KS))
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
     *               (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
               ELSE
                  IF(THCK(JS2,IS-1,KS).EQ.0.0.OR.
     *               THCK(JS,IS-1,KS).EQ.0.0) THEN
                     VR1=VR(JS,ISV1,KS)
                  ELSE
                     VR1=VR(JS2,ISV1,KS)+(1.0-(FACC-0.5))*
     *                  (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
                  END IF
               END IF
               IF(IS.EQ.NSROW) THEN
                  VR2=VR(JS2,ISV2,KS)+(1.0-(FACC-0.5))*
     *               (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
               ELSE
                  IF(THCK(JS2,IS+1,KS).EQ.0.0.OR.
     *               THCK(JS,IS+1,KS).EQ.0.0) THEN
                     VR2=VR(JS,ISV2,KS)
                  ELSE
                     VR2=VR(JS2,ISV2,KS)+(1.0-(FACC-0.5))*
     *                  (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
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
     *               (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
               ELSE
                  IF(THCK(JS2,IS-1,KS).EQ.0.0.OR.
     *            THCK(JS,IS-1,KS).EQ.0.0) THEN
                     VR1=VR(JS,ISV1,KS)
                  ELSE
                     VR1=VR(JS2,ISV1,KS)+(FACC+0.5)*
     *                  (VR(JS,ISV1,KS)-VR(JS2,ISV1,KS))
                  END IF
               END IF
               IF(IS.EQ.NSROW) THEN
                  VR2=VR(JS2,ISV2,KS)+(FACC+0.5)*
     *               (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
               ELSE
                  IF(THCK(JS2,IS+1,KS).EQ.0.0.OR.
     *               THCK(JS,IS+1,KS).EQ.0.0) THEN
                     VR2=VR(JS,ISV2,KS)
                  ELSE
                     VR2=VR(JS2,ISV2,KS)+(FACC+0.5)*
     *                  (VR(JS,ISV2,KS)-VR(JS2,ISV2,KS))
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
      KSV1=KS
      KSV2=KS+1
      FAC=OLDL-K+0.5D0
      CALL MOVTIM(VL(JS,IS,KSV1),VL(JS,IS,KSV2),VLP,VCRIT,
     *            FAC,CONST,TIML,DVDL,DISTL,ICSTVL,
     *    TINY,SMALL,HUGE)
C
C
C  TSTEP2 IS MINIMUM OF TOTAL TSTEP OR TIMES TO CELL BOUNDARIES
      TSTEP2=MIN(TSTEP,TIMC,TIMR,TIML)
C
C  CHECK TO SEE IF REACHES C BOUNDARY
      TDIFF=TIMC-TSTEP2
      IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
C  IF TDIFF SMALL, ASSUME REACHES C BOUNDARY
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
C   CHECK FOR ERROR WHEN CONVERTING BACK TO SINGLE PRECISION
C   ONLY IN FORWARD DIRECTION
c orig           IF (JS.EQ.NSCOL) THEN
C   CHECK TO SEE IF PARTICLE IS VERY CLOSE (VCLOSE) TO CELL BOUNDARY
         VCLOSE=(DBLE(J)+0.5-OLDC)
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
              IF (TESTCK.LT.SMALL.AND.TESTCK.GT.0.0)
     *            OLDC=DBLTMP*(1.D0-PRECNO)
           ENDIF
      END IF
C
C  NOW FOR R
      TDIFF=TIMR-TSTEP2
      IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
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
C   CHECK FOR ERROR WHEN CONVERTING BACK TO SINGLE PRECISION
c orig line           IF (IS.EQ.NSROW) THEN
C   CHECK TO SEE IF PARTICLE IS VERY CLOSE (VCLOSE) TO CELL BOUNDARY
         VCLOSE=(DBLE(I)+0.5-OLDR)
C   VCLOSE>0 CHECK MEANS ONLY WHEN PARTICLE IS ABOVE ROW BOUNDARY
           IF (VCLOSE.LT.SMALL.AND.VCLOSE.GT.0.0) THEN
              SNGTMP=OLDR
              DBLTMP=SNGTMP
              TESTCK=DBLTMP-OLDR
              IF (TESTCK.LT.SMALL.AND.TESTCK.GT.0.0)
     *         OLDR=DBLTMP*(1.D0-PRECNO)
           ENDIF
      END IF
C
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
              IF (TESTCK.LT.SMALL.AND.TESTCK.GT.0.0)
     *         OLDL=DBLTMP*(1.D0-PRECNO)
           ENDIF
      END IF
C
C  CHECK TO SEE IF STILL IN TRANSPORT SUBGRID
      IF(JS.LT.1.OR.JS.GT.NSCOL.OR.
     *   IS.LT.1.OR.IS.GT.NSROW.OR.
     *   KS.LT.1.OR.KS.GT.NSLAY) THEN
C  NLOC=1 IF PARTICLE IS NOT IN SUBGRID
         NLOC=1
         GO TO 348
      END IF
C  CHECK TO SEE IF NEW LOCATION IS INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) THEN
C  ONLY PRINT WARNING IF BFLX PACKAGE IS OFF 
        IF(INBFLX.EQ.0) 
     *   WRITE(IOUTS,*) '***WARNING*** PARTICLE ENTERED INACTIVE CELL;'
     * ,' WILL BE REMOVED'
	  GOTO 565
	END IF
C
C  NOW CHECK TO SEE IF MOVE COMPLETED
C
      IF(TSTEP2.LT.TSTEP) THEN
C
C  REDUCE STEP SIZE AND TAKE A NEW STEP
        TSTEP=TSTEP-TSTEP2
C1  RETURN TO TOP OF LOOP FOR NEXT STEP IF STILL IN TRANSPORT SUBGRID
        GO TO 32
      ENDIF
C1  END
C
C  END OF MOVE FOR THIS PARTICLE
      PC(IP)=OLDC
      PR(IP)=OLDR
      PL(IP)=OLDL
C  NLOC=0 SIGNIFIES PARTICLE IS WITHIN SUBGRID
      NLOC=0
C     ***************************************************************
C        ---SUM CONCENTRATIONS AND COUNT PARTICLES---
C           ---DECAY PARTICLES---
      IF(DECAY.NE.0.D0) THEN
        PCONC(IP)=DBLE(PCONC(IP))*DCYT
      ELSE IF(INDK.GT.0) THEN
        IF(IDKFO.EQ.1.OR.IDKFS.EQ.1) THEN
          CALL DK6DK(DCYT,DKFO,DKFS,RF,TIMV,
     *        IDKFO,IDKFS,
     *        JS,IS,KS,NSCOL,NSROW,NSLAY)
          PCONC(IP)=PCONC(IP)*DCYT
       END IF
      END IF
      SUMC(JS,IS,KS)=SUMC(JS,IS,KS)+PCONC(IP)
      NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)+1
C     ****************************************************************
C        ---CHECK FOR CHANGE IN CELL LOCATION---
      IF(J.EQ.J1.AND.I.EQ.I1.AND.K.EQ.K1) THEN
C  SET R COORDINATE NEGATIVE IF STILL IN ORIGINAL CELL
         IF(IORIG.EQ.1) PR(IP)=-PR(IP)
C  DONE WITH THIS PARTICLE
         GO TO 590
      END IF
C           ---CHECK FOR CONST.-HEAD BDY. OR SOURCE AT OLD LOCATION---
C  CHECK FOR FIXED HEAD BC'S
C  FOR CELLS THAT ARE SINKS AND SOURCES, SOURCES ARE FIRST
      IBD=0
C  NO SINK/SOURCES IF ICONLY=1
      IF(ICONLY.NE.1) THEN
         IF(IBOUND(J1,I1,K1).LT.0.OR.IGENPT(JS1,IS1,KS1).EQ.1.
     *   OR.IGENLK(JS1,IS1,KS1).NE.0) THEN
            IF(SRCFLO(JS1,IS1,KS1).GT.0.0) GO TO 350
            IF(SNKFLO(JS1,IS1,KS1).LT.0.0) GO TO 360
         END IF
      END IF
C  CHECK IF PARTICLE SHOULD ENTER AT SUBGRID BOUNDARY
  348 IF(ISUBGD.EQ.0) GO TO 540
C  IBD=0 SIGNIFIES PARTICLE DID NOT LEAVE INFLOW BOUNDARY CELL
C  Check boundary from which particle moves in from, and if both
C  use the one with the larger velocity to fix position, smaller to mirror
C
C  JS1,IS1,KS1 is original cell
C  JS,IS,KS is new cell
C  Velocity components are checked on the outside (subgrid) face
C    (remember that velo is stored on faces)
C
      IBD=0
C
C  These six checks look for a positive inflow from the rearward face
C  in the direction in which the particle crossed a boundary.
C  If found, that coordinate will be mirrored (see IBD logic at 525)
C  For layer direction, if particle crossed boundary in x or y, then
C  replace it at original location (this avoids corner-checks for layer
C  direction in next section)
C
c orig lines
c      IF(KS1.EQ.1.AND.VL(JS1,IS1,KS1).GT.0.0.AND.KS.GT.KS1) IBD=3
c      IF(KS1.EQ.NSLAY.AND.VL(JS1,IS1,KS1+1).LT.0.0.AND.KS.LT.KS1) IBD=3
      IF(KS1.EQ.1.AND.VL(JS1,IS1,KS1).GT.0.0) THEN
        IF(KS.GT.KS1) THEN
          IBD=3
        ELSE IF(JS.NE.JS1.OR.IS.NE.IS1) THEN
          IBD=4
        ENDIF
      ENDIF
      IF(KS1.EQ.NSLAY.AND.VL(JS1,IS1,KS1+1).LT.0.0) THEN
        IF(KS.LT.KS1) THEN
          IBD=3
        ELSE IF(JS.NE.JS1.OR.IS.NE.IS1) THEN
          IBD=4
        ENDIF
      ENDIF
      IF(JS1.EQ.1.AND.VC(JS1,IS1,KS1).GT.0.0.AND.JS.GT.JS1) IBD=1
      IF(JS1.EQ.NSCOL.AND.VC(JS1+1,IS1,KS1).LT.0.0.AND.JS.LT.JS1) IBD=1
      IF(IS1.EQ.1.AND.VR(JS1,IS1,KS1).GT.0.0.AND.IS.GT.IS1) IBD=2
      IF(IS1.EQ.NSROW.AND.VR(JS1,IS1+1,KS1).LT.0.0.AND.IS.LT.IS1) IBD=2
C
C  These four checks look for particles that originated in a corner cell.
C  Currently this only works for the XY-plane.
C  If the particle crossed a boundary in a direction (e.g. COL) but the 
C  rearward velocity in that direction (e.g. VC) is zero, we check the 
C  velocity in the other direction (e.g. VR), and if that is positive
C  we replace the particle at its original location (IBD=4)
C  If the particle crossed both boundaries we check to see if VC>VR,
C  and if so set IBD=1 (it would have been set to 2 during the 6 checks 
C  above)
C
c orig lines
c      IF(JS1.EQ.1.AND.IS1.EQ.1.AND.JS.GT.JS1.AND.IS.GT.IS1) THEN
c        IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
c      ENDIF
c      IF(JS1.EQ.NSCOL.AND.IS1.EQ.NSROW.AND.JS.LT.JS1.AND.IS.LT.IS1) THEN
c        IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
c      ENDIF
c      IF(JS1.EQ.NSCOL.AND.IS1.EQ.1.AND.JS.LT.JS1.AND.IS.GT.IS1) THEN
c        IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
c      ENDIF
c      IF(JS1.EQ.1.AND.IS1.EQ.NSROW.AND.JS.GT.JS1.AND.IS.LT.IS1) THEN
c        IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
c      ENDIF
C
C  TOP LEFT 
      IF(JS1.EQ.1.AND.IS1.EQ.1) THEN
	  IF(JS.GT.JS1.AND.IS.GT.IS1) THEN
          IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
        ELSE
	    IF(JS.GT.JS1.AND.VR(JS1,IS1,KS1).GT.0.0) IBD=4
	    IF(IS.GT.IS1.AND.VC(JS1,IS1,KS1).GT.0.0) IBD=4
        ENDIF
      ENDIF
C  BOTTOM RIGHT 
      IF(JS1.EQ.NSCOL.AND.IS1.EQ.NSROW) THEN
	  IF(JS.LT.JS1.AND.IS.LT.IS1) THEN
          IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
        ELSE
	    IF(JS.LT.JS1.AND.VR(JS1,IS1+1,KS1).GT.0.0) IBD=4
	    IF(IS.LT.IS1.AND.VC(JS1+1,IS1,KS1).GT.0.0) IBD=4
        ENDIF
      ENDIF
C  TOP RIGHT 
      IF(JS1.EQ.NSCOL.AND.IS1.EQ.1) THEN
	  IF(JS.LT.JS1.AND.IS.GT.IS1) THEN
          IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
        ELSE
	    IF(JS.LT.JS1.AND.VR(JS1,IS1,KS1).GT.0.0) IBD=4
	    IF(IS.GT.IS1.AND.VC(JS1+1,IS1,KS1).GT.0.0) IBD=4
        ENDIF
      ENDIF
C  BOTTOM LEFT 
      IF(JS1.EQ.1.AND.IS1.EQ.NSROW) THEN
	  IF(JS.GT.JS1.AND.IS.LT.IS1) THEN
          IF(ABS(VC(JS1,IS1,KS1)).GT.ABS(VR(JS1,IS1,KS1))) IBD=1
	  ELSE
	    IF(JS.GT.JS1.AND.VR(JS1,IS1+1,KS1).GT.0.0) IBD=4
	    IF(IS.LT.IS1.AND.VC(JS1,IS1,KS1).GT.0.0) IBD=4
        ENDIF
      ENDIF
C  CHECK IGENPT FOR SUBGRID BOUNDARY
C  TREAT SUBGRID IGENPT LIKE SOURCE, USE PTID
      IF(IGENPT(JS1,IS1,KS1).EQ.1.AND.
     *     (JS1.EQ.1.OR.JS1.EQ.NSCOL.OR.
     *      IS1.EQ.1.OR.IS1.EQ.NSROW.OR.
     *      KS1.EQ.1.OR.KS1.EQ.NSLAY)) IBD=4
      IF(IBD.EQ.0) GO TO 540
C     ****************************************************************
C           ---CREATE NEW PARTICLES---
  350 IF(IORIG.EQ.0) GO TO 540
C  ISOURC=1 SIGNIFIES OLD CELL IS A SOURCE
      ISOURC=1
  360 CONTINUE
C     ****************************************************************
      IF(NPTM.EQ.NPMAX) THEN
C        ---RESTART MOVE IF PT. LIMIT EXCEEDED---
        WRITE(IOUTS,700) IMOV,IP
C
C DON'T CALL SMOC5GP TWICE in ONE MOVE
        IF(INPXFL.EQ.1) THEN
          WRITE(IOUTS,710) IMOV, IP
          STOP
        ENDIF
C
        CALL SMOC5GP(PC,PR,PL,PCONC,
     *                 CONC,IPTID,NPCELL,
     *   IBOUND,PNEWC,PNEWR,PNEWL,LIMBO,
     *   NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *   NEWPTS,NPMAX,NLIMBO,
     *   IOUTS,NP,WTFAC)
C
C FLAG FIRST CALL OF SMOC5GP
        INPXFL=1
        DO 610 KS=1,NSLAY
        DO 610 IS=1,NSROW
        DO 610 JS=1,NSCOL
        SUMC(JS,IS,KS)=0.0
        NPCELL(JS,IS,KS)=0
  610   CONTINUE
        GO TO 10
      ENDIF
C     ****************************************************************
C           ---GENERATE NEW TEMPORARY PARTICLE---
      NPTM=NPTM+1
      IPN=NPTM
C
  390 ITEM=IPTID(IP)
      IF(IBD.GT.0) GO TO 525
      IF(ISOURC.EQ.1) GO TO 399
C  FOR NEW PARTICLES AT SINKS USE NODE CONCENTRATION
CHALF  NEW PARTS AT SINK IN CELL FOR ENTIRE TIME STEP, USE FULL DECAY
  398 CONTINUE
      IF(DECAY.EQ.0.D0) THEN
        IF(INDK.GT.0) THEN
          IF(IDKFO.EQ.1.OR.IDKFS.EQ.1) THEN
            CALL DK6DK(DCYT,DKFO,DKFS,RF,TIMV,
     *        IDKFO,IDKFS,
     *        JS1,IS1,KS1,NSCOL,NSROW,NSLAY)
            SUMC(JS1,IS1,KS1)=DBLE(SUMC(JS1,IS1,KS1))+
     *          DBLE(CONC(JS1,IS1,KS1))*DCYT
          ELSE
            SUMC(JS1,IS1,KS1)=SUMC(JS1,IS1,KS1)+CONC(JS1,IS1,KS1)
          END IF
        ELSE
          SUMC(JS1,IS1,KS1)=SUMC(JS1,IS1,KS1)+CONC(JS1,IS1,KS1)
        END IF
      ELSE
        SUMC(JS1,IS1,KS1)=DBLE(SUMC(JS1,IS1,KS1))+
     *                  DBLE(CONC(JS1,IS1,KS1))*DCYT
      END IF
      NPCELL(JS1,IS1,KS1)=NPCELL(JS1,IS1,KS1)+1
      IF(NPOLD(JS1,IS1,KS1).GT.0)
     *  NPOLD(JS1,IS1,KS1)=NPOLD(JS1,IS1,KS1)-1
C  INSERT NEW PARTICLE AT NODE
      PC(IPN)=J1
      PR(IPN)=I1
      PL(IPN)=K1
      IPTID(IPN)=NEWPTS
C  DON'T DECAY THIS
C  PCONC WILL BE SET IN CNCON TO NODE CONC
      PCONC(IPN)=CONC(JS1,IS1,KS1)
      GO TO 540
  399 CONTINUE
      PC(IPN)=J1+PNEWC(ITEM)
      PR(IPN)=-I1-PNEWR(ITEM)
      PL(IPN)=K1+PNEWL(ITEM)
      IPTID(IPN)=ITEM
      PCONC(IPN)=CONC(JS1,IS1,KS1)
      GO TO 540
C
C  PARTICLES AT SUBGRID BOUNDARIES
cgzh debug
      if(imov.eq.15.and.ipn.eq.766) then
	write(iouts,*)
	continue
	end if
  525 IF(IBD.EQ.1) THEN
C C-dominant flow
C  ALONG C, MIRROR C BUT FIX R AND L
      PC(IPN)=SNGL(OLDC)+JS1-JS
      PR(IPN)=-I1-PNEWR(ITEM)
      IF(KS1.EQ.1.AND.KS.GT.1.AND.ISLAY1.GT.1) THEN
         KB=ISLAY1-1
      ELSE IF(KS1.EQ.NSLAY.AND.KS.LT.NSLAY.AND.ISLAY2.LT.NLAY) THEN
         KB=ISLAY2+1
      ELSE
cgzh bug fix 8/15/06
         KB=KS1-ISLAY1+1
      END IF
      PL(IPN)=K1+PNEWL(ITEM)
C ELSE FOR IBD=2
C R-dominant flow
      ELSE IF(IBD.EQ.2) THEN
      PC(IPN)=J1+PNEWC(ITEM)
      PR(IPN)=-SNGL(OLDR)-IS1+IS
      IF(KS1.EQ.1.AND.KS.GT.1.AND.ISLAY1.GT.1) THEN
         KB=ISLAY1-1
      ELSE IF(KS1.EQ.NSLAY.AND.KS.LT.NSLAY.AND.ISLAY2.LT.NLAY) THEN
         KB=ISLAY2+1
      ELSE
cgzh bug fix 8/15/06
         KB=KS1-ISLAY1+1
      END IF
      PL(IPN)=K1+PNEWL(ITEM)
C ELSE FOR IBD=3
C L-dominant flow
      ELSE IF(IBD.EQ.3) THEN
      PC(IPN)=J1+PNEWC(ITEM)
      PR(IPN)=-I1-PNEWR(ITEM)
      IF(KS1.EQ.1.AND.KS.GT.1.AND.ISLAY1.GT.1) THEN
         KB=ISLAY1-1
      ELSE IF(KS1.EQ.NSLAY.AND.KS.LT.NSLAY.AND.ISLAY2.LT.NLAY) THEN
         KB=ISLAY2+1
      ELSE
cgzh bug fix 8/15/06
         KB=KS1-ISLAY1+1
      END IF
      PL(IPN)=SNGL(OLDL)+KS1-KS
C ELSE FOR IBD=4
C IGENPT=1 
      ELSE IF(IBD.EQ.4) THEN
         PC(IPN)=J1+PNEWC(ITEM)
         PR(IPN)=-I1-PNEWR(ITEM)
         PL(IPN)=K1+PNEWL(ITEM)
      IF(KS1.EQ.1.AND.KS.GT.1.AND.ISLAY1.GT.1) THEN
         KB=ISLAY1-1
      ELSE IF(KS1.EQ.NSLAY.AND.KS.LT.NSLAY.AND.ISLAY2.LT.NLAY) THEN
         KB=ISLAY2+1
	ELSE IF(KS1.EQ.KS) THEN
         KB=KS
      ELSE
cgzh bug fix 8/15/06
         KB=KS1-ISLAY1+1
      END IF
      ELSE
         WRITE(IOUTS,*) '***ERROR***  IBD>4 Not expected'
      END IF
      IPTID(IPN)=ITEM
C  FOR NEW PARTICLES AT SUBGRID BOUNDARY, USE EXTERNAL CONCENTRATION (INPUT)
CHALF  ENTERING PARTICLES ONLY DECAY IN SUBGRID, USE HALF DECAY
C
      IF(INDK.GT.0) THEN
        IF(IDKFO.EQ.1.OR.IDKFS.EQ.1) THEN
          TIMV2=TIMV*0.5
          CALL DK6DK(DCYT2,DKFO,DKFS,RF,TIMV2,
     *        IDKFO,IDKFS,
     *        JS1,IS1,KS1,NSCOL,NSROW,NSLAY)
        END IF
      END IF
      IF ((KB.EQ.ISLAY1-1).AND.(IABOVE.EQ.1)) THEN
         PCONC(IPN)=DBLE(CINFLA(JS1,IS1))*DCYT2
      ELSEIF ((KB.EQ.ISLAY2+1).AND.(IBELOW.EQ.1)) THEN       
         PCONC(IPN)=DBLE(CINFLB(JS1,IS1))*DCYT2
      ELSE
cgzh cinxy         PCONC(IPN)=DBLE(CINFL(KB))*DCYT2 
         PCONC(IPN)=DBLE(CINXY(JS1,IS1,KB))*DCYT2 
      ENDIF
      SUMC(JS1,IS1,KS1)=SUMC(JS1,IS1,KS1)+PCONC(IPN)
      NPCELL(JS1,IS1,KS1)=NPCELL(JS1,IS1,KS1)+1
C     ****************************************************************
C           ---CHECK FOR DISCHARGE BOUNDARY AT NEW LOCATION---
  540 IF(NLOC.GT.0) GO TO 565
C  CHECK FOR FIXED HEAD BC'S (SINKS ONLY)
      IF(IBOUND(J,I,K).GE.0.AND.IGENPT(JS,IS,KS).EQ.0.
     1AND.IGENLK(JS,IS,KS).EQ.0) GO TO 590
      IF(ICONLY.EQ.1) GO TO 590
      IF(SNKFLO(JS,IS,KS).GE.0.0) GO TO 590
C     ****************************************************************
C           ---PUT PT. IN LIMBO IF PT. DENSITY NOT INCREASED---
      IF(NPOLD(JS,IS,KS).LT.1) GO TO 590
      SUMC(JS,IS,KS)=SUMC(JS,IS,KS)-CONC(JS,IS,KS)
      NPCELL(JS,IS,KS)=NPCELL(JS,IS,KS)-1
      NPOLD(JS,IS,KS)=NPOLD(JS,IS,KS)-1
  565 CONTINUE
      PC(IP)=0.0
      PR(IP)=0.0
      PL(IP)=0.0
      PCONC(IP)=0.0
C  ZERO IPTID
      DO 570 ID=1,NLIMBO
      IF(LIMBO(ID).GT.0) GO TO 570
      LIMBO(ID)=IP
      GO TO 590
  570 CONTINUE
C
  590 CONTINUE
C     ---END OF LOOP ON OLD PARTICLES---
C     ****************************************************************
C
C        ---INSERT TEMPORARY PARTICLES INTO LIMBO LOCATIONS---
      IF(NPTM.EQ.NP) GO TO 620
      IP=NPTM
      DO 595 IL=1,NLIMBO
      IPL=LIMBO(IL)
      IF(IPL.EQ.0) GO TO 595
      PR(IPL)=PR(IP)
      PR(IP)=0.0
      PC(IPL)=PC(IP)
      PC(IP)=0.0
      PL(IPL)=PL(IP)
      PL(IP)=0.0
      PCONC(IPL)=PCONC(IP)
      PCONC(IP)=0.0
      IPTID(IPL)=IPTID(IP)
      IPTID(IP)=0
      LIMBO(IL)=0
      IP=IP-1
      IF(IP.EQ.NP) GOTO 596
  595 CONTINUE
  596 NPTM=IP
  620 CONTINUE
C        ---ADJUST NUMBER OF PARTICLES---
      NP=NPTM
C
CGWT----PRINT PARTICLE OBSERVATION DATA
C       SEND IN QMNWSINK, WILL ONLY PRINT PTS IN MNWS HERE (CALL BELOW HANDLES
C       OTHER SINKS)
      IF (INUNITPTOB.GT.0) THEN
       CALL SPTOB5O(SUMTCH,IPTOBLST,IMOV,PC,PR,PL,NPMAX,NP,TIMV,
     *        NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,NUMPTOB,NUMPTOB_MNW,
cgzh dummies sent in for arrays
     *        IPTOBMNW,NPCELL,0.0,PCONC,0.0,SNKFLO,0.0,0)
      END IF
C
C  COMMENT FOLLOWING LINE FOR LESS DETAILED OUTPUT OF MOVE LOOP
      WRITE(IOUTS,670) NP,IMOV
C
C  COMPUTE NODE CONCENTRATIONS
C     ****************************************************************
C     ---CONC. CHANGE AT NODES DUE TO ADVECTION---
C                                  ... AND DECAY---
      NZERO=0
      DO 90 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 90 IS=1,NSROW
      I=IS+ISROW1-1
      DO 90 JS=1,NSCOL
      J=JS+ISCOL1-1
      IF(IBOUND(J,I,K).EQ.0) GO TO 90
      IF(NPCELL(JS,IS,KS).EQ.0) THEN
         IF(IBOUND(J,I,K).GT.0) NZERO=NZERO+1
      ELSE
         CONC(JS,IS,KS)=SUMC(JS,IS,KS)/REAL(NPCELL(JS,IS,KS))
      END IF
      CAVG(JS,IS,KS)=0.5*(CONC(JS,IS,KS)+CNOLD(JS,IS,KS))
   90 CONTINUE
C        ---CHECK NUMBER OF CELLS VOID OF PTS.---
      IF(NZERO.GT.0) WRITE(IOUTS,290) NZERO
      IF(NZERO.GT.NZCRIT) THEN 
         WRITE(IOUTS,300)
         WRITE(IOUTS,320)
         DO 100 KS=1,NSLAY
         WRITE(IOUTS,322) KS
         DO 100 IS=1,NSROW
  100    WRITE(IOUTS,331) (NPCELL(JS,IS,KS),JS=1,NSCOL)
C  CALL SMOC5GP HERE FOR ZERO PARTICLE CELLS
C   HOWEVER, DO NOT MOVE PARTICLES AGAIN
         CALL SMOC5GP(PC,PR,PL,PCONC,
     *                 CONC,IPTID,NPCELL,
     *   IBOUND,PNEWC,PNEWR,PNEWL,LIMBO,
     *   NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *   NEWPTS,NPMAX,NLIMBO,
     *   IOUTS,NP,WTFAC)
      END IF
      DEALLOCATE(SINKTEMP)
      RETURN
C     ****************************************************************
C
  290 FORMAT (/,5X,'NUMBER OF CELLS WITH ZERO PARTICLES  =',I6)
  300 FORMAT (/,5X,'*** FZERO EXCEEDED  ---  REGENERATING ',
     * 'PARTICLES ***',/)
  320 FORMAT (/,2X,'NUMBER OF PARTICLES PER CELL',/)
  322 FORMAT(/,' FOR SUBGRID LAYER (KS) =',I6,/)
  331 FORMAT (5X,30I4)
  670 FORMAT (2X,'NP',7X,'=',I9,' AT START OF MOVE',
     1 10X,'IMOV     =',I13)
  700 FORMAT(/,5X,' ***   NOTE   ***',10X,'# Particles = NPMAX -- IMOV='
     1,I4,2X,'PT. NO.=',I8,5X,'REGENERATING PARTICLES',/)
  710 FORMAT (/,5X,' *** ERROR ***',10X,'NPMAX EXCEEDED -- IMOV='
     1,I4,2X,'PT. NO.=',I8,5X,'MUST INCREASE NPMAX')
      END
C
C
C     ***************************************************************
C
      SUBROUTINE ELMOVE(PC,PR,PL,
     *  VC,VR,VL,
     *  RF,THCK,POR,
     *  IBOUND,
     *  NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,
     *  IOUTS,TIMV,TSTEP2,
     *  VCMAX,VRMAX,VLMAX,JS,IS,KS,J,I,K,NLOC,
     *  DELCOL,DELROW)
C
C     ***************************************************************
      DOUBLE PRECISION TINY,SMALL,PRECNO
      DOUBLE PRECISION OLDR,OLDC,OLDL,VCP,VRP,VLP,DISTC,DISTR,DISTL
      DOUBLE PRECISION DBLTMP,TESTCK,DVDC,DVDL,DVDR
      PARAMETER (PERCNT=0.01)
      PARAMETER (TINY=1.D-20)
      PARAMETER (SMALL=1.D-4)
      PARAMETER (HUGE=1.E20)
C  CONSTANT FOR CHECK OF SINGLE PRECISION INTERVAL
      PARAMETER (PRECNO=5.D-7)
      DIMENSION
     *  VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *  VL(NSCOL,NSROW,NSLAY+1),
     *  RF(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),
     *  IBOUND(NCOL,NROW,NLAY),DELCOL(NCOL),DELROW(NROW)
C
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
C     ***************************************************************
C3D  ONE MOVE ONLY
C
C  VCRIT IS MAXIMUM RELATIVE VELOCITY MULTIPLIED
C    BY CRITERION FOR IGNORING V CHANGE
C  USER CAN CHANGE 0.01 (1 PERCENT) CRITERION BELOW
      VCRIT=MAX(VCMAX,VRMAX,VLMAX)*PERCNT
C
C     ***************************************************************
      IF ((PC.LT.ISCOL1-.5).
     * OR.(PC.GT.ISCOL2+.5).
     * OR.(PR.LT.ISROW1-.5).
     * OR.(PR.GT.ISROW2+.5).
     * OR.(PL.LT.ISLAY1-.5).
     * OR.(PL.GT.ISLAY2+.5)) THEN
             WRITE(IOUTS,*) ' PC,PR,PL=',PC,PR,PL
             WRITE(IOUTS,*) ' MOVE ERROR, STOPPING'
             STOP ' ERROR IN MOVE'
c      ELSE                               
c              IF (ISCOL1-PC.EQ..5) PC=PC+SMALL
c        IF (ISROW1-PR.EQ..5) PR=PR+SMALL
c              IF (ISLAY1-PL.EQ..5) PL=PL+SMALL
c              IF (PC-ISCOL2.EQ..5) PC=PC-SMALL
c              IF (PR-ISROW2.EQ..5) PR=PR-SMALL
c              IF (PL-ISLAY2.EQ..5) PL=PL-SMALL
      ENDIF
c***************************************************************
C
      OLDC=PC
      J=INT(OLDC+0.5D0)
      IF (J.GT.ISCOL2) J=ISCOL2
      JS=J-ISCOL1+1
      OLDR=PR
      I=INT(OLDR+0.5D0)
      IF (I.GT.ISROW2) I=ISROW2
      IS=I-ISROW1+1
      OLDL=PL
      K=INT(OLDL+0.5D0)
      IF (K.GT.ISLAY2) K=ISLAY2
      KS=K-ISLAY1+1
C  
C  ****************************************   
C     KEEP FOR DEBUGGING
C  POINT LOCATED IN INACTIVE CELL SHOULD NOT OCCUR
      IF(IBOUND(J,I,K).EQ.0) THEN
         WRITE(IOUTS,*) ' J,I,K=',J,I,K
         WRITE(IOUTS,*) ' OLDC,OLDR,OLDL=',OLDC,OLDR,OLDL
         WRITE(IOUTS,*) 'MOVE POINT FROM INACTIVE CELL' 
      ENDIF
C  ****************************************
C
C1  USE ANALYTIC EXPRESSION WITHIN BLOCK FOR LINEAR V
      TSTEP=TIMV
C1  THIS IS BEGINNING OF LOOP FOR PARTIAL TIME STEPS
C1  CONST CONVERTS Q IN REAL UNITS TO RETARDED V IN RELATIVE UNITS
   32 ARINV=1/(DELCOL(J)*DELROW(I))
cellam
      CONST=ARINV/(RF(JS,IS,KS)*THCK(JS,IS,KS)*POR(JS,IS,KS))
C1  FIRST FOR C
      JSV1=JS
      JSV2=JS+1
      V1=VC(JSV1,IS,KS)
      V2=VC(JSV2,IS,KS)
C  EMOVTIM COMPUTES TIME TO MOVE TO CELL BOUNDARY IN GIVEN DIRECTION
      FACC=OLDC-J+0.5D0
      CALL EMOVTIM(V1,V2,VCP,VCRIT,
     *            FACC,CONST,TIMC,DVDC,DISTC,ICSTVC,
     *   TINY,SMALL,HUGE)
C  NOW FOR R
      ISV1=IS
      ISV2=IS+1
      V1=VR(JS,ISV1,KS)
      V2=VR(JS,ISV2,KS)
      FACR=OLDR-I+0.5D0
      CALL EMOVTIM(V1,V2,VRP,VCRIT,
     *            FACR,CONST,TIMR,DVDR,DISTR,ICSTVR,
     *   TINY,SMALL,HUGE)
C  AND FOR L
      KSV1=KS
      KSV2=KS+1
      V1=VL(JS,IS,KSV1)
      V2=VL(JS,IS,KSV2)
      FACL=OLDL-K+0.5D0
      CALL EMOVTIM(V1,V2,VLP,VCRIT,
     *            FACL,CONST,TIML,DVDL,DISTL,ICSTVL,
     *    TINY,SMALL,HUGE)
C
C
C  TSTEP2 IS MINIMUM OF TOTAL TSTEP OR TIMES TO CELL BOUNDARIES
      TSTEP2=MIN(TSTEP,TIMC,TIMR,TIML)
C
C  CHECK TO SEE IF REACHES C BOUNDARY
      TDIFF=TIMC-TSTEP2
      IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
C  IF TDIFF SMALL, ASSUME REACHES C BOUNDARY
      IF(TDIFF.LT.SMALL) THEN
C  MOVE TO BOUNDARY IF STEP TIME IS GREATER OR EQUAL TO
C    TIME TO C BOUNDARY
         IF(DISTC.GT.0.0D0) THEN
            J=J+1
            JS=JS+1
c            OLDC=OLDC+DISTC+SMALL
c            OLDC=OLDC+DISTC
          if(abs(oldc+distc-j+0.5D0).gt.small) print *,'cdiff gt small'
            OLDC=J-0.5D0
         ELSE
            IF(DISTC.LT.0.0D0) THEN
              J=J-1
              JS=JS-1
c             OLDC=OLDC+DISTC-SMALL
c             OLDC=OLDC+DISTC
          if(abs(oldc+distc-j-0.5D0).gt.small) print *,'cdiff gt small'
              OLDC=J+0.5D0
            ELSE
              IF (FACC.EQ.1D0 .AND. VCP.GT.0D0) THEN
                J=J+1
                JS=JS+1
              ELSE
                IF (FACC.EQ.0D0 .AND. VCP.LT.0D0) THEN
                  J=J-1
                  JS=JS-1
                else
        print *,'error in assumptions for move: dist = 0d0 and
     *              fac/velocity unexpected',distc,facc,vcp
                endif
              ENDIF
            ENDIF
         END IF
      ELSE
         IF(ICSTVC.EQ.0) THEN
C1  LINEAR VARIATION IN VC
            OLDC=VCP*(EXP(DVDC*DBLE(TSTEP2))-1.0D0)/DVDC+OLDC
         ELSE IF(ICSTVC.EQ.1) THEN
C1  CONSTANT VC
            OLDC=VCP*TSTEP2+OLDC
         END IF
      END IF
C
C  NOW FOR R
      TDIFF=TIMR-TSTEP2
      IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
      IF(TDIFF.LT.SMALL) THEN
         IF(DISTR.GT.0.0D0) THEN
            I=I+1
            IS=IS+1
c            OLDR=OLDR+DISTR+SMALL
C            OLDR=OLDR+DISTR
          if(abs(oldr+distr-i+0.5D0).gt.small) print *,'rdiff gt small'
            OLDR=I-0.5D0
         ELSE
            IF(DISTR.LT.0.0D0) THEN
              I=I-1
              IS=IS-1
c              OLDR=OLDR+DISTR-SMALL
C              OLDR=OLDR+DISTR
          if(abs(oldr+distr-i-0.5D0).gt.small) print *,'rdiff gt small'
              OLDR=I+0.5D0
            ELSE
              IF (FACR.EQ.1D0 .AND. VRP.GT.0D0) THEN
                I=I+1
                IS=IS+1
              ELSE
                IF (FACR.EQ.0D0 .AND. VRP.LT.0D0) THEN
                  I=I-1
                  IS=IS-1
                else
        print *,'error in assumptions for move: dist = 0d0 and
     *              fac/velocity unexpected',distr,facr,vrp
                endif
              ENDIF
            ENDIF
         END IF
      ELSE
         IF(ICSTVR.EQ.0) THEN
            OLDR=VRP*(EXP(DVDR*DBLE(TSTEP2))-1.D0)/DVDR+OLDR
         ELSE IF(ICSTVR.EQ.1) THEN
            OLDR=VRP*DBLE(TSTEP2)+OLDR
         END IF
      END IF
        if (abs(oldr-i).gt..5) print *,'oldr-i gt .5',oldr,i,tdiff,
     *  ICSTVR,
     *  distr
C
C  NOW L
C
      TDIFF=TIML-TSTEP2
      IF(TDIFF.LT.HUGE) TDIFF=TDIFF/TIMV
      IF(TDIFF.LT.SMALL) THEN
         IF(DISTL.GT.0.0D0) THEN
            K=K+1
            KS=KS+1
c            OLDL=OLDL+DISTL+SMALL
C            OLDL=OLDL+DISTL
          if(abs(oldl+distl-k+0.5D0).gt.small) print *,'ldiff gt small'
            OLDL=K-0.5D0
         ELSE
            IF(DISTL.LT.0.0D0) THEN
              K=K-1
              KS=KS-1
c              OLDL=OLDL+DISTL-SMALL
C              OLDL=OLDL+DISTL
          if(abs(oldl+distl-k-0.5D0).gt.small) print *,'ldiff gt small'     
              OLDL=K+0.5D0
            ELSE
              IF (FACL.EQ.1D0 .AND. VLP.GT.0D0) THEN
                K=K+1
                KS=KS+1
              ELSE
                IF (FACL.EQ.0D0 .AND. VLP.LT.0D0) THEN
                  K=K-1
                  KS=KS-1
                else
       print *,'error in assumptions for move: dist = 0d0 and
     *              fac/velocity unexpected',distl,facl,vlp
                endif
              ENDIF
            ENDIF
         END IF
      ELSE
         IF(ICSTVL.EQ.0) THEN
            OLDL=VLP*(EXP(DVDL*DBLE(TSTEP2))-1.0D0)/DVDL+OLDL
         ELSE IF(ICSTVL.EQ.1) THEN
            OLDL=VLP*DBLE(TSTEP2)+OLDL
         END IF
      END IF
C
C  CHECK TO SEE IF STILL IN TRANSPORT SUBGRID
      NLOC=0
      IF(JS.LT.1.OR.JS.GT.NSCOL.OR.
     *   IS.LT.1.OR.IS.GT.NSROW.OR.
     *   KS.LT.1.OR.KS.GT.NSLAY) THEN
         NLOC=1
         GO TO 348
      END IF
C
C  NOW CHECK TO SEE IF MOVE COMPLETED
C
      IF(TSTEP2.LT.TSTEP) THEN
C
C  REDUCE STEP SIZE AND TAKE A NEW STEP
        TSTEP=TSTEP-TSTEP2
C1  RETURN TO TOP OF LOOP FOR NEXT STEP IF STILL IN TRANSPORT SUBGRID
        GO TO 32
      ENDIF
C1  END
C
C  END OF MOVE FOR THIS PARTICLE
348   PC=OLDC
      PR=OLDR
      PL=OLDL
C
      if (iscol1+js-1.ne.j .or. isrow1+is-1.ne.i .or. islay1+ks-1.ne.k)
     *   print *,'move error',iscol1,isrow1,islay1,js,is,ks,j,i,k,
     *           iscol1+js-1,isrow1+is-1,islay1+ks-1
      RETURN
      END
C
C  MOVTIM  COMPUTES TIME TO BOUNDARY OF CELL
C     ***************************************************************
C
      SUBROUTINE MOVTIM(V1,V2,VP,VCRIT,
     *              FAC,CONST,TIME,DVDL,DIST,ICNST,
     *   TINY,SMALL,HUGE)
C
C     ***************************************************************
      DOUBLE PRECISION TINY,SMALL
      DOUBLE PRECISION VP,DVDL,DIST,RATIO,VB
C
C  CHECK TO SEE IF BOUNDARY VELOCITIES ARE EQUAL
      V12DIF=V2-V1
      IF(ABS(V12DIF).LT.TINY) THEN
C  CHECK TO SEE IF BOTH V'S ARE ZERO
         IF(ABS(V1).LT.TINY) THEN
            TIME=HUGE
            ICNST=-1
            RETURN
         END IF
C  FOR NOZERO EQUAL V'S, USE CONSTANT V
         VP=DBLE(V1)*DBLE(CONST)
         IF(VP.GT.0.0D0) THEN
            DIST=1.0-FAC
         ELSE
            DIST=-FAC
         END IF
         TIME=DIST/VP
         ICNST=1
         RETURN
      END IF
C  UNEQUAL V'S
C  INTERPOLATE TO GET V OF PARTICLE
      VP=(V1+FAC*V12DIF)*CONST
      IF(VP.GT.0.0D0) THEN
         IF(V2.LE.0.0) THEN
C  IF BOUNDARY V2 IN OPPOSITE DIRECTION
C   FIND ZERO V POSITION AND SET TIME TO INFINITY
            DIST=-V1/V12DIF-FAC
            TIME=HUGE
            ICNST=0
            DVDL=-VP/DIST
            RETURN
         END IF
C  CHECK TO SEE IF CONSTANT V APPROX OK
         DIST=1.0-FAC
         IF(VP.LT.VCRIT) THEN
            ICNST=1
            TIME=DIST/VP
            RETURN
         END IF
C  SET BOUNDARY V TO V2
         VB=V2*CONST
      ELSE IF(VP.LT.0.0D0) THEN
         IF(V1.GE.0.0) THEN
C  IF POSITIVE V1, FIND ZERO V POSITION AND SET TIME TO INFINITY
            DIST=-V1/V12DIF-FAC
            TIME=HUGE
            ICNST=0
            DVDL=-VP/DIST
            RETURN
         END IF
         DIST=-FAC
C  CHECK IF CONSTANT V APPROX OK
         IF(-VP.LT.VCRIT) THEN
            ICNST=1
            TIME=DIST/VP
            RETURN
         END IF
C  SET BOUNDARY V TO V1
         VB=V1*CONST
      ELSE
         TIME=HUGE
         ICNST=-1
         RETURN
      END IF
C
C  UNEQUAL V'S, BOTH IN SAME DIRECTION, RATIO > 0
      RATIO=VB/VP
C  CRITERION FOR USE OF LOG FORM CAN BE SET HERE BY USER
      IF(ABS(RATIO-1.0D0).LT.SMALL.OR.
     *  ABS(VP).LT.VCRIT) THEN
         TIME=DIST/VP
         ICNST=1
      ELSE
         DVDL=(VB-VP)/DIST
         TIME=LOG(RATIO)/DVDL
         ICNST=0
      END IF
      RETURN
      END
C
C
C  SMOC5GP  
C     ***************************************************************
C
      SUBROUTINE SMOC5GP(PC,PR,PL,PCONC,
     *                 CONC,IPTID,NPCELL,
     *   IBOUND,PNEWC,PNEWR,PNEWL,LIMBO,
     *   NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *   NEWPTS,NPMAX,NLIMBO,
     *   IOUTS,NP,WTFAC)
C
C     ***************************************************************
C
C     GENERATE PARTICLE POSITIONS
C
cgzh vlwt
      DIMENSION PC(NPMAX),PR(NPMAX),PL(NPMAX),PCONC(NPMAX),
     *     CONC(NSCOL,NSROW,NSLAY),IPTID(NPMAX),
     *   NPCELL(NSCOL,NSROW,NSLAY),WTFAC(NSCOL,NSROW,NSLAY),
     *   PNEWC(NEWPTS),PNEWR(NEWPTS),PNEWL(NEWPTS),LIMBO(NLIMBO)
      DIMENSION IBOUND(NCOL,NROW,NLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
CMOCWT
      INCLUDE 'ptwt.inc'
C     ***************************************************************
C
C  INITIALIZE PARTICLE ARRAYS
   10 DO 20 IP=1,NPMAX
      IPTID(IP)=0
      PC(IP)=0.0
      PR(IP)=0.0
      PL(IP)=0.0
      PCONC(IP)=0.0
   20 CONTINUE
C
C     ---SET UP LIMBO ARRAY---
      DO 40 IN=1,NLIMBO
   40 LIMBO(IN)=0
      IND=0
C
C  NPTPND IS NEWPTS MINUS 1
      NPTPND=NEWPTS-1
C     ***************************************************************
C     ---INSERT PARTICLES---
C        ---TRACK PARTICLE LOCATIONS IN COORDINATES OF PRIMARY GRID---
      NACTIV=0
      DO 410 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 410 IS=1,NSROW
      I=IS+ISROW1-1
      DO 410 JS=1,NSCOL
      J=JS+ISCOL1-1
      NPCELL(JS,IS,KS)=0
      IF(IBOUND(J,I,K).EQ.0) THEN
         IF(CONC(JS,IS,KS).NE.CNOFLO) THEN
cgzh condition this: only for unweighted 
          IF(PTWTON.EQ.0) THEN
            CONC(JS,IS,KS)=CNOFLO
            WRITE(IOUTS,4498) J,I,K
 4498 FORMAT(' *** WARNING *** TRANSPORT CELL WENT DRY (J,I,K)=',3I4)
          END IF
         END IF
         GO TO 410
      END IF
      NACTIV=NACTIV+1
      NPCELL(JS,IS,KS)=NPTPND
      C1=CONC(JS,IS,KS)
C
C  USE NEW LOCATIONS
      DO 280 IP=1,NPTPND
      IND=IND+1
      IF(IND.GT.NPMAX) THEN
         WRITE(IOUTS,500) IND
         WRITE(IOUTS,*) ' ** ERROR ** -- NPMAX EXCEEDED, STOPPING'
         WRITE(IOUTS,*) '       INCREASE NPMAX AND RERUN'
         STOP
      END IF
      PC(IND)=J+PNEWC(IP)
      PR(IND)=-I-PNEWR(IP)
      PL(IND)=K+PNEWL(IP)
      PCONC(IND)=C1
      IPTID(IND)=IP
  280 CONTINUE
C
  410 CONTINUE
      NP=IND
      WRITE(IOUTS,500) NP
  500 FORMAT(/,' TOTAL NUMBER OF PARTICLES GENERATED =',I10)
      NZCRIT=INT(FZERO*NACTIV)
      WRITE(IOUTS,501) NACTIV
  501 FORMAT(' TOTAL NUMBER OF ACTIVE NODES (NACTIV) =',I10)
      WRITE(IOUTS,502) NZCRIT
  502 FORMAT(' MAX. NUMBER OF CELLS THAT CAN BE VOID OF '
     1       ,'PARTICLES (NZCRIT) = ',I6/5X,
     2       '(IF NZCRIT EXCEEDED, PARTICLES ARE REGENERATED)')
C
C     ****************************************************************
      RETURN
C     ****************************************************************
      END
C
C
C  EMOVTIM  COMPUTES TIME TO BOUNDARY OF CELL (ELLAM)
C     ***************************************************************
C
      SUBROUTINE EMOVTIM(V1,V2,VP,VCRIT,
     *              FAC,CONST,TIME,DVDL,DIST,ICNST,
     *   TINY,SMALL,HUGE)
C
C     ***************************************************************
      DOUBLE PRECISION TINY,SMALL
      DOUBLE PRECISION VP,DVDL,DIST,RATIO,VB
C
C  CHECK TO SEE IF BOUNDARY VELOCITIES ARE EQUAL
      V12DIF=V2-V1
C  CHECK TO SEE IF BOTH V'S ARE ZERO
        if(abs(V12DIF)-tiny) 100,100,200
c  velocities are equal
100     continue
         IF(ABS(V1).LT.TINY) THEN
            TIME=HUGE
            ICNST=-1
         else
C  FOR NOZERO EQUAL V'S, USE CONSTANT V
            VP=DBLE(V1)*DBLE(CONST)
            IF(VP.GT.0.0D0) THEN
               DIST=1D0-FAC
            ELSE
               DIST=-FAC
            END IF
            TIME=DIST/VP
            ICNST=1
         endif
c
         RETURN
c
200     continue
C  UNEQUAL V'S
C  INTERPOLATE TO GET V OF PARTICLE
      VP=(V1+FAC*V12DIF)*CONST
      IF(VP) 400,500,300
c*** case VP > 0
300     continue
         IF(V2.LE.0.D0) THEN
C  IF BOUNDARY V2 IN OPPOSITE DIRECTION
C   FIND ZERO V POSITION AND SET TIME TO INFINITY
            DIST=-V1/V12DIF-FAC
            TIME=HUGE
            ICNST=0
            DVDL=-VP/DIST
            RETURN
         END IF
C  CHECK TO SEE IF CONSTANT V APPROX OK
         DIST=1D0-FAC
         IF(VP.LT.VCRIT) THEN
            ICNST=1
            TIME=DIST/VP
            RETURN
         END IF
C  SET BOUNDARY V TO V2
         VB=V2*CONST
         goto 600
c
c*** case VP < 0
400     continue
         IF(V1.GE.0.D0) THEN
C  IF POSITIVE V1, FIND ZERO V POSITION AND SET TIME TO INFINITY
            DIST=-V1/V12DIF-FAC
            TIME=HUGE
            ICNST=0
            DVDL=-VP/DIST
            RETURN
         END IF
         DIST=-FAC
C  CHECK IF CONSTANT V APPROX OK
         IF(-VP.LT.VCRIT) THEN
            ICNST=1
            TIME=DIST/VP
            RETURN
         END IF
C  SET BOUNDARY V TO V1
         VB=V1*CONST
        goto 600
c*** case VP=0
500     continue
         TIME=HUGE
         ICNST=-1
        return
c
600     continue
C
C  UNEQUAL V'S, BOTH IN SAME DIRECTION, RATIO > 0
      RATIO=VB/VP
C  CRITERION FOR USE OF LOG FORM CAN BE SET HERE BY USER
      IF(ABS(RATIO-1.0D0).LT.SMALL.OR.
     *  ABS(VP).LT.VCRIT) THEN
         TIME=DIST/VP
         ICNST=1
      ELSE
         DVDL=(VB-VP)/DIST
         TIME=LOG(RATIO)/DVDL
         ICNST=0
      END IF
      RETURN
      END
