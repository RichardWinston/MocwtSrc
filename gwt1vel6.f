C     Last change: RBW Dec. 8, 2014
C     Support for weighted particles (MOCWT) added.
C VELO   CALCULATE VELOCITIES
C**************************************************************
C
      SUBROUTINE VELO(BUFF,VC,VR,VL,IBOUND,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,
     *  IDIR,kkstp,
cgzh vlwt
     *  THCK,POR,
     *  LAYHDT,HNEW,HOLD,BOTM,NBOTM,DELT,ISS,DELCOL,DELROW)
C
C**************************************************************
C
      DIMENSION BUFF(NCOL,NROW,NLAY),
     *   VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1),
     *   IBOUND(NCOL,NROW,NLAY),
     *   THCK(NSCOL,NSROW,NSLAY),POR(NSCOL,NSROW,NSLAY)
C
cgzh vlwt
      DOUBLE PRECISION HNEW
      DIMENSION DELCOL(NCOL),DELROW(NROW)
C
      DIMENSION HNEW(NCOL,NROW,NLAY),
     *  HOLD(NCOL,NROW,NLAY),LAYHDT(NLAY),
     *  BOTM(NCOL,NROW,0:NBOTM)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
C
C  USE MODFLOW COMPUTATION OF FLUX TERMS   
C  MODFLOW T TERMS YIELD VOLUMETRIC FLUX THROUGH FACE  L3/T
C  INDEXING FOR VELOCITY VX(IXV, ,) IS FLOW FROM IXV-1 TO IXV
C     ***************************************************************
      IF(IDIR.EQ.1) THEN
C
C
C        ---VELOCITIES AT CELL BOUNDARIES---
C           TERM SAVED IS VELOCITY*THICKNESS*POROSITY*RDEL = VOL FLUX
C
         JSV2=NSCOL+1
         DO 10 KS=1,NSLAY
         K=KS+ISLAY1-1
         DO 10 IS=1,NSROW
         I=IS+ISROW1-1
         DO 10 JSV=1,JSV2
            VC(JSV,IS,KS)=0.0
            J=JSV+ISCOL1-1
            IF(J.EQ.1.OR.J.GT.NCOL) GO TO 10
            IF(IBOUND(J,I,K).EQ.0) GO TO 10
C  SKIP REST IF VELOCITY IS ZERO
            IF(BUFF(J-1,I,K).EQ.0.0) GO TO 10
            VC(JSV,IS,KS)=BUFF(J-1,I,K)
 10      CONTINUE
C  ROW DIRECTION
      ELSE IF(IDIR.EQ.2) THEN
         ISV2=NSROW+1
         DO 11 KS=1,NSLAY
         K=KS+ISLAY1-1
         DO 11 ISV=1,ISV2
         IS=ISV
         I=IS+ISROW1-1
         DO 11 JS=1,NSCOL
            VR(JS,ISV,KS)=0.0
            IF(I.EQ.1.OR.I.GT.NROW) GO TO 11
            J=JS+ISCOL1-1
            IF(IBOUND(J,I,K).EQ.0) GO TO 11
            IF(BUFF(J,I-1,K).EQ.0.0) GO TO 11
            VR(JS,ISV,KS)=BUFF(J,I-1,K)
 11      CONTINUE
      ELSE
C
C  Z TERM = VELOCITY*POROSITY*CDEL*RDEL
         KSV2=NSLAY+1
         DO 12 KSV=1,KSV2
         KS=KSV
         K=KS+ISLAY1-1
         DO 12 IS=1,NSROW
         I=IS+ISROW1-1
         DO 12 JS=1,NSCOL
         J=JS+ISCOL1-1
           VL(JS,IS,KSV)=0.0
           IF(K.EQ.1.OR.K.GT.NLAY) GO TO 12
           IF(IBOUND(J,I,K).EQ.0) GO TO 12
cgzh alter velocities due to water table drop/rise
c SKIPPING THIS!!!!!!!!!!
       go to 444
C BEGIN WATER TABLE VELOCITY ADJUSTMENT
C ONLY DO FOR CONVERTIBLE LAYERS WITH CHANGES IN HEAD IN TRANSIENT MODE
            IF(LAYHDT(K).EQ.1.AND.(HNEW(J,I,K).NE.HOLD(J,I,K))
     &        .AND.ISS.EQ.0) THEN
C IF IN FIRST LAYER, MUST BE WT
cgzh debug
      if(js.eq.4.and.ks.eq.3.and.(kkstp.eq.8.or.kkstp.eq.9)) then
	continue
	endif
              IUPDATE=0
              IF(K.EQ.1) THEN
			  IUPDATE=1
C IF NOT IN FIRST LAYER, ONLY WT IF CELL ABOVE IS DRY
			ELSE IF (IBOUND(J,I,K-1).EQ.0) THEN
			  IUPDATE=1
			END IF 
              IF(IUPDATE.EQ.1) THEN
  	          TTOP=BOTM(J,I,LBOTM(K)-1)
C DO NOT PROCEED IF BOTH HEADS ARE ABOVE TOP
                IF(.NOT.(HNEW(J,I,K).GT.TTOP.AND
     &		          .HOLD(J,I,K).GT.TTOP)) THEN
C SPECIAL CASE: HNEW IS > top & HOLD < top, INCREASE VEL USING TOP, NOT HNEW
cgzh debug (if the cell above is rewettable, however, can use HNEW
	            IF(HNEW(J,I,K).GT.TTOP.AND.HOLD(J,I,K).LT.TTOP) THEN
                    H1=TTOP
                  ELSE
                    H1=HNEW(J,I,K)
                  ENDIF
                  VLWT=(HOLD(J,I,K)-H1)/DELT
cgzh debug output
c      write(iouts,*) 'VLWT AT CELL JS,IS,KS',VLWT,JS,IS,KS      
C VL IS NOT REALLY A VOL FLUX, SO HERE MULTIPLY BY CONST SO THE VLWT WILL COME OUT
C   RIGHT IN INTERPOLATION IN MOVE
c Rf here too
                  CONST=DELCOL(J)*DELROW(I)*THCK(JS,IS,KS)*POR(JS,IS,KS)
                  FLXWT=VLWT*CONST
cgzh debug output
c      write(iouts,*) 'WT VOL FLUX AT CELL JS,IS,KS',FLXWT,JS,IS,KS
                  VL(JS,IS,KSV)=VL(JS,IS,KSV)+FLXWT
                END IF  
              END IF
            END IF
C END WATER TABLE VELOCITY ADJUSTMENT
 444  CONTINUE
c continue for skip
            IF(K.EQ.1.OR.K.GT.NLAY) GO TO 12
            IF(BUFF(J,I,K-1).EQ.0.0) GO TO 12
            VL(JS,IS,KSV)=BUFF(J,I,K-1)
 12      CONTINUE
      END IF
C
C***************************************************
C
C  END LOOPS TO COMPUTE VELOCITIES
C
C***************************************************
      RETURN
C***************************************************
C
      END
C
C
C VELO   CALCULATE MAXIMUM VELOCITIES
C**************************************************************
C
      SUBROUTINE VELOMAX(THCK,POR,RF,
     *                VC,VR,VL,IBOUND,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,IOUTS,
     *  VCMAX,VRMAX,VLMAX,TLMIN,IDIR,
     *  MAXVCJ,MAXVCI,MAXVCK,MAXVRJ,MAXVRI,MAXVRK,
     *  MAXVLJ,MAXVLI,MAXVLK,
cellam
     *  MOCTYPE,TCMIN,TRMIN,DELCOL,DELROW)
cellam
C**************************************************************
C
      DIMENSION THCK(NSCOL,NSROW,NSLAY),
     *   POR(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *   VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1),
     *   IBOUND(NCOL,NROW,NLAY)
cellam
      DIMENSION DELCOL(NCOL),DELROW(NROW)
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
cellam
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C  USE MODFLOW COMPUTATION OF FLUX TERMS   
C  MODFLOW T TERMS YIELD VOLUMETRIC FLUX THROUGH FACE  L3/T
C  INDEXING FOR VELOCITY VX(IXV, ,) IS FLOW FROM IXV-1 TO IXV
C     ***************************************************************
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        DCINV=1.D0/CDEL
        DRINV=1.D0/RDEL
      ENDIF
C
C  INITIALIZE MAX VELOCITY VALUES AND LOCATIONS FIRST TIME THROUGH
         VCMAX=1.0E-20
         MAXVCJ=0
         MAXVCI=0
         MAXVCK=0
         VRMAX=1.0E-20
         MAXVRJ=0
         MAXVRI=0
         MAXVRK=0
         VLMAX=1.0E-20
cellam
         IF(MOCTYPE.EQ.3) THEN
           TCMIN=1.0D+20
           TRMIN=1.0D+20
         ENDIF
cellam
C  USE TIME FOR VERTICAL BECAUSE THICKNESS MAY BE VARIABLE
         TLMIN=1.0E+20
         MAXVLJ=0
         MAXVLI=0
         MAXVLK=0
C
C
C        ---VELOCITIES AT CELL BOUNDARIES---
C           TERM SAVED IS VELOCITY*THICKNESS*POROSITY*RDEL = VOL FLUX
C
         JSV2=NSCOL+1
         DO 10 KS=1,NSLAY
         K=KS+ISLAY1-1
         DO 10 IS=1,NSROW
         I=IS+ISROW1-1
         DO 10 JSV=1,JSV2
            J=JSV+ISCOL1-1
            IF(J.EQ.1.OR.J.GT.NCOL) GO TO 10
            IF(IBOUND(J,I,K).EQ.0) GO TO 10
C  SKIP REST IF VELOCITY IS ZERO
            IF(VC(JSV,IS,KS).EQ.0.0) GO TO 10
C  IF THIS CELL WITHIN SUBGRID, COMPUTE ABSOLUTE MAXIMUM V
C   USE MINIMUM OF TWO ADJACENT THICK*POR*RF
            JS=JSV
            IF(JS.GT.1.AND.JS.LE.NSCOL) THEN
               PORBrf=THCK(JS,IS,KS)*POR(JS,IS,KS)*rf(js,is,ks)
c boundary flux  thickness = 0 but vc > 0 can occur; check ibound for backward cell 
               IF(IBOUND(J-1,I,K).EQ.0) THEN
                 ABVC=ABS(VC(JSV,IS,KS))/PORBrf
               ELSE
                 PRBrf2=THCK(JS-1,IS,KS)*POR(JS-1,IS,KS)*rf(js-1,is,ks)
                 ABVC=ABS(VC(JSV,IS,KS))/AMIN1(PORBrf,PRBrf2)  
               END IF
C   IF CELL OUTSIDE SUBGRID, USE THICK*POR*RF FOR ADJACENT CELL WITHIN SUBGRID
            ELSE IF(JS.GT.NSCOL) THEN
                ABVC=ABS(VC(JSV,IS,KS))/
     *            (THCK(JS-1,IS,KS)*POR(JS-1,IS,KS)*rf(js-1,is,ks))
            ELSE
C  THIS REACHED ONLY IF JS=1, USE THIS THICK*POR*RF
               ABVC=ABS(VC(JSV,IS,KS))/
     *            (THCK(JS,IS,KS)*POR(JS,IS,KS)*rf(js,is,ks))
            END IF
C  CONVERT TO LINEAR VELOCITY
cellam
            IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN 
              ABVC=ABVC*DRINV
            ELSEIF(MOCTYPE.EQ.3) THEN
              RDRINV=1.0/DELROW(I)
              ABVC=ABVC*RDRINV
              TCOL=DELCOL(J)/ABVC
              IF(TCOL.LT.TCMIN) TCMIN=TCOL
            ENDIF
cellam
            IF(ABVC.GT.VCMAX) THEN
               VCMAX=ABVC
               MAXVCI=I
               MAXVCJ=J
               MAXVCK=K
            END IF
 10      CONTINUE
C
         ISV2=NSROW+1
         DO 11 KS=1,NSLAY
         K=KS+ISLAY1-1
         DO 11 ISV=1,ISV2
         IS=ISV
         I=IS+ISROW1-1
         DO 11 JS=1,NSCOL
            IF(I.EQ.1.OR.I.GT.NROW) GO TO 11
            J=JS+ISCOL1-1
            IF(IBOUND(J,I,K).EQ.0) GO TO 11
            IF(VR(JS,ISV,KS).EQ.0.0) GO TO 11
            IF(IS.GT.1.AND.IS.LE.NSROW) THEN
               PORBrf=THCK(JS,IS,KS)*POR(JS,IS,KS)*rf(js,is,ks)
               IF(IBOUND(J,I-1,K).EQ.0) THEN
                 ABVR=ABS(VR(JS,ISV,KS))/PORBrf
               ELSE
                 PRBrf2=THCK(JS,IS-1,KS)*POR(JS,IS-1,KS)*rf(js,is-1,ks)
                 ABVR=ABS(VR(JS,ISV,KS))/AMIN1(PORBrf,PRBrf2)
               END IF
            ELSE IF(IS.GT.NSROW) THEN
                ABVR=ABS(VR(JS,ISV,KS))/
     *           (THCK(JS,IS-1,KS)*POR(JS,IS-1,KS)*rf(js,is-1,ks))
            ELSE
               ABVR=ABS(VR(JS,ISV,KS))/
     *          (THCK(JS,IS,KS)*POR(JS,IS,KS)*rf(js,is,ks))
            END IF
C  CONVERT
cellam
            IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN 
              ABVR=ABVR*DCINV
            ELSEIF(MOCTYPE.EQ.3) THEN
              RDCINV=1.0/DELCOL(J)
              ABVR=ABVR*RDCINV
              TROW=DELROW(I)/ABVR
              IF(TROW.LT.TRMIN) TRMIN=TROW
            ENDIF
cellam
            IF(ABVR.GT.VRMAX) THEN
               VRMAX=ABVR
               MAXVRI=I
               MAXVRJ=J
               MAXVRK=K
            END IF
 11      CONTINUE
C
C  Z TERM = VELOCITY*POROSITY*CDEL*RDEL
         KSV2=NSLAY+1
         DO 12 KSV=1,KSV2
         KS=KSV
         K=KS+ISLAY1-1
         DO 12 IS=1,NSROW
         I=IS+ISROW1-1

         DO 12 JS=1,NSCOL
            IF(K.EQ.1.OR.K.GT.NLAY) GO TO 12
            J=JS+ISCOL1-1
            IF(IBOUND(J,I,K).EQ.0) GO TO 12
            IF(VL(JS,IS,KSV).EQ.0.0) GO TO 12
            IF(KS.GT.1.AND.KS.LE.NSLAY) THEN          
                PORBRF=THCK(JS,IS,KS)*POR(JS,IS,KS)*RF(js,is,KS)
                rf1=rf(js,is,ks)
               IF(IBOUND(J,I,K-1).EQ.0) THEN
                ABVL=ABS(VL(JS,IS,KSV))/(POR(JS,IS,KS)*RF1)
                TLAYB=THCK(JS,IS,KS)/ABVL
               ELSE 
			  PRBRF2=THCK(JS,IS,KS-1)*POR(JS,IS,KS-1)*RF(js,is,KS-1)
                rf2=rf(js,is,ks-1)
                IF(PRBRF2.LT.PORBRF) THEN
                  ABVL=ABS(VL(JS,IS,KSV))/(POR(JS,IS,KS-1)*RF2)
                  TLAYB=THCK(JS,IS,KS-1)/ABVL
                ELSE
                  ABVL=ABS(VL(JS,IS,KSV))/(POR(JS,IS,KS)*RF1)
                  TLAYB=THCK(JS,IS,KS)/ABVL
                END IF
               END IF
		  ELSE IF(KS.GT.NSLAY) THEN
               ABVL=ABS(VL(JS,IS,KSV))/(POR(JS,IS,KS-1)*RF(js,is,KS-1))
               TLAYB=THCK(JS,IS,KS-1)/ABVL
            ELSE
               ABVL=ABS(VL(JS,IS,KSV))/(POR(JS,IS,KS)*RF(js,is,KS))
               TLAYB=THCK(JS,IS,KS)/ABVL
            END IF
cellam
            IF(MOCTYPE.EQ.3) THEN
              AREA=DELCOL(J)*DELROW(I)
              TLAYB=TLAYB*AREA
            END IF
cellam
            IF(TLAYB.LT.TLMIN) THEN
               TLMIN=TLAYB
cellam
            IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
               VLMAX=ABVL
            ELSEIF(MOCTYPE.EQ.3) THEN
               VLMAX=ABVL/AREA
            END IF
cellam
               MAXVLI=I
               MAXVLJ=J
               MAXVLK=K
            END IF
 12      CONTINUE
C
C  SCALE VERTICAL TERMS BY AREA
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        TLMIN=TLMIN*CDEL*RDEL
        VLMAX=VLMAX/(CDEL*RDEL)
      END IF
C
C***************************************************
C
C  END LOOPS TO COMPUTE MAXIMUM VELOCITIES
C
C***************************************************
      RETURN
C***************************************************
C
      END
