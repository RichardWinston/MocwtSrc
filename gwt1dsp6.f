C     Last change:  RBW, July 24, 2015
C     Support for weighted particles (MOCWT) added.
C  DSP6   DISPERSION
C
C  DSP6FM  DISPERSION COEFFICIENTS
C***********************************************
C
      SUBROUTINE DSP6FM(DISPCC,
     *  DISPCR,DISPCL,DISPRR,DISPRC,DISPRL,DISPLL,
     *  DISPLC,DISPLR,THCK,ALONG,ATRANH,ATRANV,
     *  POR,RF,VC,VR,VL,IBOUND,MOCTYPE,
     *  DELT,NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *  IOUTS,NTIMD,TIMDC,NODISP,DIFFUS,
cellam
     *  DELCOL,DELROW)
cellam
C
C FOR ELLAM, MODIFIED TO
C  - YIELD CORRECT CELL DIMENSION IN DENOMINATOR
C  - HANDLE NONUNIFORM GRID.
C  THIS VERSION RETAINS GWT APPROXIMATION IGNORING VARIABLE
C  CELL THICKNESS IN DENOMINATOR
C***********************************************
C
      DIMENSION
     *  DISPCC(NSCOL,NSROW,NSLAY),DISPCR(NSCOL,NSROW,NSLAY),
     *  DISPCL(NSCOL,NSROW,NSLAY),DISPRR(NSCOL,NSROW,NSLAY),
     *  DISPRC(NSCOL,NSROW,NSLAY),DISPRL(NSCOL,NSROW,NSLAY),
     *  DISPLL(NSCOL,NSROW,NSLAY),DISPLC(NSCOL,NSROW,NSLAY),
     *  DISPLR(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  ALONG(NSLAY),ATRANH(NSLAY),ATRANV(NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY),
     *  VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *  VL(NSCOL,NSROW,NSLAY+1),IBOUND(NCOL,NROW,NLAY)
C
cellam
      DIMENSION DELCOL(NCOL),DELROW(NROW)
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
cellam
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C  BUFF HOLDS AVG CONC; DSP BY FACE, NOT BLOCK
C  INDEXING FOR VELOCITY VX(IXV, ,) IS FLOW FROM IXV-1 TO IXV
C
C     ***************************************************************
C
C     ---COMPUTE DISPERSION COEFFICIENTS---
C  SPATIAL TERMS IN HORIZONTAL DIRECTIONS
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
       DCINV=1.D0/CDEL
       CC2INV=DCINV*DCINV
       DRINV=1.D0/RDEL
       RR2INV=DRINV*DRINV
       CR2INV=0.25D0*DCINV*DRINV
       ARINV=DCINV*DRINV
      ENDIF
cellam
C
C  VERTICAL DIMENSION DEPENDS ON CELLS
C
C*****************************************
C
C  BEGIN LOOP FOR DISPERSION COEFFICIENTS
C
C*****************************************
C
      DO 60 KS=1,NSLAY
      KSV1=KS
      KSV2=KS+1
      K=KS+ISLAY1-1
C  DISPERSIVITIES FOR THIS LAYER
      IF(NODISP.NE.1) THEN
         ALONGL=ALONG(KS)
         ATRNH=ATRANH(KS)
         ATRNV=ATRANV(KS)
C  TAKE HARMONIC MEAN OF DISPERSIVITIES FOR BETWEEN LAYERS
         IF(KS.LT.NSLAY) THEN
            IF(ALONGL.EQ.0.0.OR.ALONG(KS+1).EQ.0.0) THEN
               ALONGV=0.0D0
            ELSE
               ALONGV=2.D0*ALONGL*ALONG(KS+1)/(ALONGL+ALONG(KS+1))
            END IF
            IF(ATRNV.EQ.0.0.OR.ATRANV(KS+1).EQ.0.0) THEN
               ATRNVV=0.0D0
            ELSE
               ATRNVV=2.D0*ATRNV*ATRANV(KS+1)/(ATRNV+ATRANV(KS+1))
            END IF
         END IF
      END IF
      DO 60 IS=1,NSROW
      ISV1=IS
      ISV2=IS+1
      I=IS+ISROW1-1
cellam
      IF(MOCTYPE.EQ.3) THEN
        DIINV=1.D0/DELROW(I)
        R4INV=0.25D0*DIINV
      END IF
cellam
      DO 60 JS=1,NSCOL
      JSV1=JS
      JSV2=JS+1
      J=JS+ISCOL1-1
C  INITIALIZE DISPERSION COEFFICIENTS TO ZERO
      DISPCC(JS,IS,KS)=0.0D0
      DISPRR(JS,IS,KS)=0.0D0
      DISPLL(JS,IS,KS)=0.0D0
      IF(NODISP.NE.1) THEN
         DISPCR(JS,IS,KS)=0.0D0
         DISPRC(JS,IS,KS)=0.0D0
         DISPCL(JS,IS,KS)=0.0D0
         DISPRL(JS,IS,KS)=0.0D0
         DISPLC(JS,IS,KS)=0.0D0
         DISPLR(JS,IS,KS)=0.0D0
      END IF
C
C  SKIP REST IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 60
C
      THCKNS=THCK(JS,IS,KS)
      PORB1=THCKNS*POR(JS,IS,KS)
C  COMPUTE DELTA Z (L) AT KS+1/2 FROM CELL THICKNESSES
C   AND SET SPATIAL TERMS FOR DISP COEFFICIENTS
      IF(KS.LT.NSLAY) THEN
         DELINV=2.D0/(THCKNS+THCK(JS,IS,KS+1))
cellam
         IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
cgzh orig           DL2INV=DELINV*DELINV
           CL2INV=0.25D0*DCINV*DELINV
           C2INV=0.5D0*DCINV
           c4inv = 0.25d0*dcinv
           RL2INV=0.25D0*DRINV*DELINV
           R2INV=0.5D0*DRINV
           r4inv = 0.25d0*drinv
         ELSEIF(MOCTYPE.EQ.3) THEN
           L4INV=0.25D0*DELINV
         END IF
cellam
      END IF
C        ---FORWARD COEFFICIENTS: COL-DIRECTION---
C  SKIP IF THIS IS LAST SOLUTE TRANSPORT COLUMN, NO DISPERSION
C   ACROSS SUBGRID BOUNDARY, OR IF CELL IN NEXT COLUMN IS INACTIVE
cellam
      IF(MOCTYPE.EQ.3) THEN
        DJINV=1.D0/DELCOL(J)
        C4INV=0.25D0*DJINV
      END IF
cellam
      IF(JS.EQ.NSCOL) GO TO 20
cellam
      IF(MOCTYPE.EQ.3) THEN
        DCINV=2.D0/(DELCOL(J)+DELCOL(J+1))
        DJ1INV=1.D0/DELCOL(J+1)
      END IF
cellam
      IF(IBOUND(J+1,I,K).EQ.0) GO TO 20
C  IF DIFFUSION ONLY, SET TERM AND SKIP TO ROW TERMS
C
C  MOLECULAR DIFFUSION COEFFICIENT INCLUDES TORTUOSITY
C   HAS SAME UNITS AS (ALPHA)V, MULTIPLY BY POROSITY LIKE SUTRA
C   ALSO MULTIPLY BY THICKNESS FOR HORIZONTAL TERMS USING
C   HARMONIC MEAN INTERBLOCK VALUE
      IF(DIFFUS.NE.0.0) THEN
cellam
         IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
           PORB2=POR(JS+1,IS,KS)*THCK(JS+1,IS,KS)
           DISPCC(JS,IS,KS)=DIFFUS*2.D0*PORB1*PORB2/
     *                    (PORB1+PORB2)*CC2INV
         ELSEIF(MOCTYPE.EQ.3) THEN
           PORB2=POR(JS+1,IS,KS)*THCK(JS+1,IS,KS)
           DISPCC(JS,IS,KS)=DIFFUS*2.D0*PORB1*PORB2/
     *                    (PORB1+PORB2)*DCINV
         ENDIF
cellam
      END IF
      IF(NODISP.EQ.1) GO TO 20
C  COMPUTE VELOCITY COMPONENTS AT (J+1/2,I,K)
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        VC2=VC(JSV2,IS,KS)*DRINV
        VC22=VC2*VC2
C  BINLINEAR INTERPOLATION FOR VR AND VL
        VRC2=(VR(JS,ISV1,KS)+VR(JS+1,ISV1,KS)+
     *        VR(JS,ISV2,KS)+VR(JS+1,ISV2,KS))*0.25D0*DCINV
        VRC22=VRC2*VRC2
        VLC2=(VL(JS,IS,KSV1)+VL(JS+1,IS,KSV1)+
     *        VL(JS,IS,KSV2)+VL(JS+1,IS,KSV2))*0.25D0*THCKNS*ARINV
        VLC22=VLC2*VLC2
      ELSEIF(MOCTYPE.EQ.3) THEN
        VC2=VC(JSV2,IS,KS)*DIINV
        VC22=VC2*VC2
C  BINLINEAR INTERPOLATION FOR VR AND VL
        VRC2=(VR(JS,ISV1,KS)+VR(JS,ISV2,KS))*0.25D0*DJINV+
     *     (VR(JS+1,ISV1,KS)+VR(JS+1,ISV2,KS))*0.25D0*DJ1INV
        VRC22=VRC2*VRC2
        VLC2=((VL(JS,IS,KSV1)+VL(JS,IS,KSV2))*DJINV+
     *  (VL(JS+1,IS,KSV1)+VL(JS+1,IS,KSV2))*DJ1INV)*0.25D0*THCKNS*DIINV
        VLC22=VLC2*VLC2
      ENDIF
cellam
C  MAGNITUDE OF VELOCITY AT (J+1/2,I,K)
      VMGC22=VC22+VRC22+VLC22
      VMGC2=SQRT(VMGC22)
C  SKIP IF VELOCITY IS MACHINE ZERO
      IF(VMGC2.LT.1.0E-20) GO TO 20
C
C           ---CC COEFFICIENT---
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        DISPCC(JS,IS,KS)=DISPCC(JS,IS,KS)+
     *         (ALONGL*VC22+ATRNH*VRC22+ATRNV*VLC22)/VMGC2*CC2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPCC(JS,IS,KS)=DISPCC(JS,IS,KS)+
     *         (ALONGL*VC22+ATRNH*VRC22+ATRNV*VLC22)/VMGC2*DCINV
      ENDIF
cellam
C  SKIP IF FIRST OR LAST ROW IN SUBGRID
      IF(IS.EQ.1.OR.IS.EQ.NSROW) GO TO 10
C           ---CR COEFFICIENT---
C  SKIP IF ANY OF CELLS USED FOR D(CONC)/DR ARE INACTIVE
C   SKIP TO ROW IF DISPERSIVITIES ARE EQUAL (NO CROSS TERMS)
      IF(IBOUND(J,I-1,K).EQ.0.OR.IBOUND(J+1,I-1,K).EQ.0.OR.
     *   IBOUND(J,I+1,K).EQ.0.OR.IBOUND(J+1,I+1,K).EQ.0.OR.
     *   ALONGL.EQ.ATRNH)
     *    GO TO 10
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        DISPCR(JS,IS,KS)=(ALONGL-ATRNH)*VC2*VRC2/VMGC2*CR2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPCR(JS,IS,KS)=(ALONGL-ATRNH)*VC2*VRC2/VMGC2*R4INV
      ENDIF
cellam
C           ---CL COEFFICIENT---
C  SKIP IF FIRST OR LAST LAYER IN SUBGRID
   10 IF(KS.EQ.1.OR.KS.EQ.NSLAY) GO TO 20
C  SKIP IF ANY OF CELLS USED FOR D(CONC)/DL ARE INACTIVE
      IF(IBOUND(J,I,K-1).EQ.0.OR.IBOUND(J+1,I,K-1).EQ.0.OR.
     *   IBOUND(J,I,K+1).EQ.0.OR.IBOUND(J+1,I,K+1).EQ.0.OR.
     *   IBOUND(J+1,I,K).EQ.0.OR.
     *   ALONGL.EQ.ATRNV)
     *    GO TO 20
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        DISPCL(JS,IS,KS)=(ALONGL-ATRNV)*VC2*VLC2/VMGC2*C2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPCL(JS,IS,KS)=(ALONGL-ATRNV)*VC2*VLC2/VMGC2*L4INV
      ENDIF
cellam
C        ---FORWARD COEFFICIENTS: ROW-DIRECTION---
   20 CONTINUE
cellam
      IF(MOCTYPE.EQ.3) ARINV=DJINV*DIINV
cellam
      IF(IS.EQ.NSROW) GO TO 40
cellam
      IF(MOCTYPE.EQ.3) THEN
        DRINV=2.D0/(DELROW(I)+DELROW(I+1))
        DI1INV=1.D0/DELROW(I+1)
      END IF
cellam
      IF(IBOUND(J,I+1,K).EQ.0) GO TO 40
      IF(DIFFUS.NE.0.0) THEN
cellam
        IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
          PORB2=POR(JS,IS+1,KS)*THCK(JS,IS+1,KS)
          DISPRR(JS,IS,KS)=DIFFUS*2.D0*PORB1*PORB2/
     *                    (PORB1+PORB2)*RR2INV
        ELSEIF(MOCTYPE.EQ.3) THEN
          PORB2=POR(JS,IS+1,KS)*THCK(JS,IS+1,KS)
          DISPRR(JS,IS,KS)=DIFFUS*2.D0*PORB1*PORB2/
     *                    (PORB1+PORB2)*DRINV
        END IF
cellam
      END IF
      IF(NODISP.EQ.1) GO TO 40
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        VR2=VR(JS,ISV2,KS)*DCINV
        VR22=VR2*VR2
        VCR2=(VC(JSV1,IS,KS)+VC(JSV1,IS+1,KS)+
     *      VC(JSV2,IS,KS)+VC(JSV2,IS+1,KS))*0.25*DRINV
        VCR22=VCR2*VCR2
        VLR2=(VL(JS,IS,KSV1)+VL(JS,IS+1,KSV1)+
     *      VL(JS,IS,KSV2)+VL(JS,IS+1,KSV2))*0.25*THCKNS*ARINV
        VLR22=VLR2*VLR2
      ELSEIF(MOCTYPE.EQ.3) THEN
        VR2=VR(JS,ISV2,KS)*DJINV
        VR22=VR2*VR2
        VCR2=(VC(JSV1,IS,KS)+VC(JSV2,IS,KS))*0.25D0*DIINV+
     *     (VC(JSV1,IS+1,KS)+VC(JSV2,IS+1,KS))*0.25D0*DI1INV
        VCR22=VCR2*VCR2
        DOVC=THCKNS*DJINV
        VLR2=(VL(JS,IS,KSV1)+VL(JS,IS,KSV2))*0.25D0*DOVC*DIINV+
     *     (VL(JS,IS+1,KSV1)+VL(JS,IS+1,KSV2))*0.25D0*DOVC*DI1INV
        VLR22=VLR2*VLR2
      END IF
cellam
      VMGR22=VCR22+VR22+VLR22
      VMGR2=SQRT(VMGR22)
      IF(VMGR2.LT.1.0E-20) GO TO 40
C           ---RR COEFFICIENT---
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN      
        DISPRR(JS,IS,KS)=DISPRR(JS,IS,KS)+
     *         (ALONGL*VR22+ATRNH*VCR22+ATRNV*VLR22)/VMGR2*RR2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPRR(JS,IS,KS)=DISPRR(JS,IS,KS)+
     *         (ALONGL*VR22+ATRNH*VCR22+ATRNV*VLR22)/VMGR2*DRINV
      ENDIF
cellam
C   SKIP TO LAYER IF DISPERSIVITIES ARE EQUAL (NO CROSS TERMS)
      IF(JS.EQ.1.OR.JS.EQ.NSCOL) GO TO 30
C           ---RC COEFFICIENT---
      IF(IBOUND(J-1,I,K).EQ.0.OR.IBOUND(J-1,I+1,K).EQ.0.OR.
     *   IBOUND(J+1,I,K).EQ.0.OR.IBOUND(J+1,I+1,K).EQ.0.OR.
     *   ALONGL.EQ.ATRNH)
     *    GO TO 30
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN      
        DISPRC(JS,IS,KS)=(ALONGL-ATRNH)*VCR2*VR2/VMGR2*CR2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPRC(JS,IS,KS)=(ALONGL-ATRNH)*VCR2*VR2/VMGR2*C4INV
      ENDIF
cellam
C           ---RL COEFFICIENT---
   30 IF(KS.EQ.1.OR.KS.EQ.NSLAY) GO TO 40
      IF(IBOUND(J,I,K-1).EQ.0.OR.IBOUND(J,I+1,K-1).EQ.0.OR.
     *   IBOUND(J,I,K+1).EQ.0.OR.IBOUND(J,I+1,K+1).EQ.0.OR.
     *   IBOUND(J,I+1,K).EQ.0.OR.
     *   ALONGL.EQ.ATRNV)
     *    GO TO 40
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN      
        DISPRL(JS,IS,KS)=(ALONGL-ATRNV)*VR2*VLR2/VMGR2*R2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPRL(JS,IS,KS)=(ALONGL-ATRNV)*VR2*VLR2/VMGR2*L4INV
      ENDIF
cellam
   40 CONTINUE
C  Z TERMS
C  HERE WE DON'T WANT THICKNESS INCLUDED, BECAUSE MAY VARY BY LAYER
      IF(KS.EQ.NSLAY) GO TO 60
      IF(IBOUND(J,I,K+1).EQ.0) GO TO 60
      IF(DIFFUS.NE.0.0) THEN
cellam
        IF(MOCTYPE.EQ.1) THEN
cgzh orig         DISPLL(JS,IS,KS)=DIFFUS*2.D0*POR(JS,IS,KS)*POR(JS,IS,KS+1)/
cgzh orig     *                    (POR(JS,IS,KS)+POR(JS,IS,KS+1))*DL2INV
         DISPLL(JS,IS,KS)=DIFFUS*2.D0*POR(JS,IS,KS)*POR(JS,IS,KS+1)/
     *                    (POR(JS,IS,KS)+POR(JS,IS,KS+1))*DELINV
        ELSEIF(MOCTYPE.EQ.2.or.MOCTYPE.EQ.3) THEN
         DISPLL(JS,IS,KS)=DIFFUS*2.D0*POR(JS,IS,KS)*POR(JS,IS,KS+1)/
     *                    (POR(JS,IS,KS)+POR(JS,IS,KS+1))*DELINV
        END IF
cellam
      END IF
      IF(NODISP.EQ.1) GO TO 60
C
      VL2=VL(JS,IS,KSV2)*ARINV
      VL22=VL2*VL2
cellam
      IF(MOCTYPE.EQ.1.OR.MOCTYPE.EQ.2) THEN
        VCL2=((VC(JSV1,IS,KS)+VC(JSV2,IS,KS))/THCKNS+
     *      (VC(JSV1,IS,KS+1)+VC(JSV2,IS,KS+1))/THCK(JS,IS,KS+1))
     *      *0.25*DRINV
        VCL22=VCL2*VCL2
        VRL2=((VR(JS,ISV1,KS)+VR(JS,ISV2,KS))/THCKNS+
     *      (VR(JS,ISV1,KS+1)+VR(JS,ISV2,KS+1))/THCK(JS,IS,KS+1))
     *      *0.25*DCINV
      ELSEIF(MOCTYPE.EQ.3) THEN
        VCL2=((VC(JSV1,IS,KS)+VC(JSV2,IS,KS))/THCKNS+
     *      (VC(JSV1,IS,KS+1)+VC(JSV2,IS,KS+1))/THCK(JS,IS,KS+1))
     *      *0.25*DIINV
        VCL22=VCL2*VCL2
        VRL2=((VR(JS,ISV1,KS)+VR(JS,ISV2,KS))/THCKNS+
     *      (VR(JS,ISV1,KS+1)+VR(JS,ISV2,KS+1))/THCK(JS,IS,KS+1))
     *      *0.25*DJINV
      END IF
cellam
      VRL22=VRL2*VRL2
      VMGL22=VCL22+VRL22+VL22
      VMGL2=SQRT(VMGL22)
      IF(VMGL2.LT.1.0E-20) GO TO 60
C           ---LL COEFFICIENT---
cellam
      IF(MOCTYPE.EQ.1) THEN
cgzh orig        DISPLL(JS,IS,KS)=DISPLL(JS,IS,KS)+
cgzh orig     *           (ALONGV*VL22+ATRNVV*(VCL22+VRL22))/VMGL2*DL2INV
        DISPLL(JS,IS,KS)=DISPLL(JS,IS,KS)+
     *           (ALONGV*VL22+ATRNVV*(VCL22+VRL22))/VMGL2*DELINV
      ELSEIF(MOCTYPE.EQ.2.or.MOCTYPE.EQ.3) THEN
        DISPLL(JS,IS,KS)=DISPLL(JS,IS,KS)+
     *           (ALONGV*VL22+ATRNVV*(VCL22+VRL22))/VMGL2*DELINV
      END IF
cellam
C   SKIP TO END IF DISPERSIVITIES ARE EQUAL (NO CROSS TERMS)
      IF(ALONGV.EQ.ATRNVV) GO TO 60
      DSPDIF=ALONGV-ATRNVV
C           ---LC COEFFICIENT---
      IF(JS.EQ.1.OR.JS.EQ.NSCOL) GO TO 50
      IF(IBOUND(J-1,I,K).EQ.0.OR.IBOUND(J-1,I,K+1).EQ.0.OR.
     *   IBOUND(J+1,I,K).EQ.0.OR.IBOUND(J+1,I,K+1).EQ.0)
     *    GO TO 50
cellam
      IF(MOCTYPE.EQ.1) THEN
        DISPLC(JS,IS,KS)=DSPDIF*VCL2*VL2/VMGL2*C2INV
      ELSEIF(MOCTYPE.EQ.2) THEN
        DISPLC(JS,IS,KS)=DSPDIF*VCL2*VL2/VMGL2*CL2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPLC(JS,IS,KS)=DSPDIF*VCL2*VL2/VMGL2*C4INV
      END IF
cellam
C           ---LR COEFFICIENT---
   50 IF(IS.EQ.1.OR.IS.EQ.NSROW) GO TO 60
      IF(IBOUND(J,I-1,K).EQ.0.OR.IBOUND(J,I-1,K+1).EQ.0.OR.
     *   IBOUND(J,I+1,K).EQ.0.OR.IBOUND(J,I+1,K+1).EQ.0)
     *    GO TO 60
cellam
      IF(MOCTYPE.EQ.1) THEN
        DISPLR(JS,IS,KS)=DSPDIF*VRL2*VL2/VMGL2*R2INV
      ELSEIF(MOCTYPE.EQ.2) THEN
        DISPLR(JS,IS,KS)=DSPDIF*VRL2*VL2/VMGL2*RL2INV
      ELSEIF(MOCTYPE.EQ.3) THEN
        DISPLR(JS,IS,KS)=DSPDIF*VRL2*VL2/VMGL2*R4INV
      END IF
cellam
C
   60 CONTINUE
C
C*****************************************
C
C  END LOOP FOR DISPERSION COEFFICIENTS
C
C*****************************************
C
      if (MOCTYPE.ne.1) return
C
C     ****************************************************************
C     ---CHECK FOR STABILITY OF EXPLICIT METHOD---
C     ****************************************************************
C
      TIMDIS=1.E-20
      DO 70 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 70 IS=1,NSROW
      I=IS+ISROW1-1
      DO 70 JS=1,NSCOL
      J=JS+ISCOL1-1
      IF(IBOUND(J,I,K).EQ.0) GO TO 70
      IF(JS.EQ.1) THEN
         DISPC=DISPCC(JS,IS,KS)
      ELSE
         DISPC=MAX(DISPCC(JS,IS,KS),DISPCC(JS-1,IS,KS))
      END IF
      IF(IS.EQ.1) THEN
       DISPR=DISPRR(JS,IS,KS)
      ELSE
       DISPR=MAX(DISPRR(JS,IS,KS),DISPRR(JS,IS-1,KS))
      END IF
      IF(KS.EQ.1) THEN
       DISPL=DISPLL(JS,IS,KS)
      ELSE
       DISPL=MAX(DISPLL(JS,IS,KS),DISPLL(JS,IS,KS-1))
      END IF
cgzh orig      TDCO=((DISPC+DISPR)/THCK(JS,IS,KS)+DISPL)/
cgzh orig     *    (POR(JS,IS,KS)*RF(js,is,KS))
c the 1/b term for displ has been moved to dsp6ap, so not in displ here
      TDCO=((DISPC+DISPR+DISPL)/THCK(JS,IS,KS))/
     *    (POR(JS,IS,KS)*RF(js,is,KS))
      IF(TDCO.GT.TIMDIS) TIMDIS=TDCO
   70 CONTINUE
      TIMDC=0.5/TIMDIS
C  CHECK STABILITY 
      NTIMD=DELT/TIMDC
      NTIMD=NTIMD+1
      RETURN
      END
C********************************************************
C
      SUBROUTINE DSP6AP(
     *  DISPCC,DISPCR,DISPCL,
     *  DISPRR,DISPRC,DISPRL,
     *  DISPLL,DISPLC,DISPLR,
     *  THCK,POR,RF,CONC,CNCNC,
CMOCWT  npcell
     *  IBOUND,CAVG,NPCELL,SUMWT,CELVOL,
     *  NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,
     *  IOUTS,KSTP,KPER,IMOV,TIMV,
     *  NODISP,
c chfb
     *  LOCHFB,ONCHFB)
C
C********************************************************
C
C  CHANGES CONCENTRATIONS DUE TO DISPERSION AND SOURCES FOR ONE STEP ONLY
C  COMPUTE EXPLICIT CHANGE IN NODE CONC FOR FULL TIME STEP
CLMT  LIMIT CROSS PRODUCT TERMS
CLMT  CHECK FOR NEGATIVE CONCS IN ADJACENT CELLS
CFWD  COMPUTE BY FACE INSTEAD OF BY BLOCK
C
CLMT
      DOUBLE PRECISION CELVOL
      DOUBLE PRECISION ADJFCT,DXXCUR,DXXFOR,DYYCUR,DYYFOR,DZZCUR,DZZFOR
      DOUBLE PRECISION TEMFCT
CMOCWT
      DOUBLE PRECISION SUMWT
      ALLOCATABLE ADJFCT(:,:,:),
     *  DXXCUR(:,:,:),DXXFOR(:,:,:),
     *  DYYCUR(:,:,:),DYYFOR(:,:,:),
     *  DZZCUR(:,:,:),DZZFOR(:,:,:)
      DIMENSION
     *  DISPCC(NSCOL,NSROW,NSLAY),DISPCR(NSCOL,NSROW,NSLAY),
     *  DISPCL(NSCOL,NSROW,NSLAY),DISPRR(NSCOL,NSROW,NSLAY),
     *  DISPRC(NSCOL,NSROW,NSLAY),DISPRL(NSCOL,NSROW,NSLAY),
     *  DISPLL(NSCOL,NSROW,NSLAY),DISPLC(NSCOL,NSROW,NSLAY),
     *  DISPLR(NSCOL,NSROW,NSLAY),  THCK(NSCOL,NSROW,NSLAY),
     *     POR(NSCOL,NSROW,NSLAY),    RF(NSCOL,NSROW,NSLAY),
     *    CONC(NSCOL,NSROW,NSLAY), CNCNC(NSCOL,NSROW,NSLAY),
     *  IBOUND(NCOL,NROW,NLAY),CAVG(NSCOL,NSROW,NSLAY),
c chfb
     *  LOCHFB(NSCOL,NSROW,NSLAY,2)
CMOCWT
      DIMENSION  NPCELL(NSCOL,NSROW,NSLAY),SUMWT(NSCOL,NSROW,NSLAY),
     *  CELVOL(NSCOL,NSROW,NSLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
CMOCWT
      INCLUDE 'ptwt.inc'
CLMT
      ALLOCATE(ADJFCT(NSCOL,NSROW,NSLAY),
     * DXXCUR(NSCOL,NSROW,NSLAY),DXXFOR(NSCOL,NSROW,NSLAY),
     * DYYCUR(NSCOL,NSROW,NSLAY),DYYFOR(NSCOL,NSROW,NSLAY),
     * DZZCUR(NSCOL,NSROW,NSLAY),DZZFOR(NSCOL,NSROW,NSLAY))
C
C INITIALIZE TEMP ARRAYS
      DXXCUR=0.0
      DYYCUR=0.0
      DZZCUR=0.0
      DXXFOR=0.0
      DYYFOR=0.0
      DZZFOR=0.0
C     ***************************************************************
C  SKIP TO NODE CONCS IF NO DISP,DIFFUS,OR SINK/SOURCE OR DECAY
C     ***************************************************************
C     ---CONC. CHANGE DUE TO:
      DO 20 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 20 IS=1,NSROW
      I=IS+ISROW1-1
      DO 20 JS=1,NSCOL
      J=JS+ISCOL1-1
CLMT  SET ADJUSTMENT FACTOR TO -1.0 TO FLAG INACTIVE CELL FOR ADJUSTMENT LOOP
      ADJFCT(JS,IS,KS)=-1.D0
C     SKIP IF INACTIVE CELL
      IF(IBOUND(J,I,K).EQ.0) GO TO 20
CMOCWT  SKIP IF NO PARTICLES IN CELL
      IF(PTWTON.EQ.1.AND.NPCELL(JS,IS,KS).LT.1) GO TO 20
C
CFWD  DELTA T / (RF POR THCK) FOR THIS CELL
      TOPRF=timv/(POR(JS,IS,KS)*rf(js,is,ks))
      TRFPRB=TOPRF/THCK(JS,IS,KS)
      C2=CAVG(JS,IS,KS)
C         ---CONC. CHANGE DUE TO DISPERSION FOR TIMV---
CLMT  CHECK STABILITY OF TERMS
c orig      CELLMT=-MAX(C2,0.0)/TRFPRB
c ptlmt
c      PTLMT=SUMMASS(JS,IS,KS)/SUMWT(JS,IS,KS)
c      CELLMT=-MAX(PTLMT,C2,0.0)/TRFPRB
cgzh debug output
c      if (ptlmt.gt.c2.and.ptlmt.gt.0.0) 
c     *      write(iouts,*) 'ptlmt for js,is,ks',js,is,ks
C        ---DISPERSION WITH TENSOR COEFFICIENTS---
CFM  DISP FLUX IN FORWARD COLUMN DIRECTION
      IF(JS.LT.NSCOL) THEN
      IF(IBOUND(J+1,I,K).NE.0) THEN
CMOCWT  SKIP TO ROW DIRECTION IF NO PARTICLES IN FORWARD CELL
       IF(PTWTON.EQ.1.AND.NPCELL(JS+1,IS,KS).LT.1) GO TO 15
         tPRBIN2=timv/(POR(JS+1,IS,KS)*THCK(JS+1,IS,KS)*rf(js+1,is,ks))
CFM  FLUX DUE TO DXX
         DSPXX2=DISPCC(JS,IS,KS)*(CAVG(JS+1,IS,KS)-C2)
cgzh debug output
c      if(js.eq.2.and.is.eq.1.and.imov.eq.5) then
c       write(iouts,*) 'DSPXX2 2,1',DSPXX2
c       write(iouts,*) 'DISPCC,CAVG,C2 2,1',DISPCC(JS,IS,KS),
c     *CAVG(JS+1,IS,KS),C2
c	end if
cgzh debug output
c      if(dspxx2.ne.0.0) then
c        write(iouts,*) 'DSPXX2 js is ks',DSPXX2, js, is, ks
c        write(iouts,*) 'CAVG(JS+1),C2',CAVG(JS+1,IS,KS),C2
c        write(iouts,*) 'DSPCUR',DSPXX2*TRFPRB
c      end if
         IF(NODISP.NE.1) THEN
CFWD  FLUX DUE TO DXY
c orig
c            IF(IS.GT.1.AND.IS.LT.NSROW)
c     *        DSPXX2=DSPXX2+DISPCR(JS,IS,KS)*
c     *         (CAVG(JS,IS+1,KS)+CAVG(JS+1,IS+1,KS)
c     *         -CAVG(JS,IS-1,KS)-CAVG(JS+1,IS-1,KS))
c
c define C for cells around node in question
c for current layer, relative to C2 which is at JS,IS,KS:
c
c         C11  C12
c    C20 (C2)  C22
c    C30  C31  C32
c
c for layer above (k-1)
c
c         C41  C42
c    C50  C51  C52
c    C60  C61  C62
c
c for layer below (k+1)
c
c         C71  C72
c    C80  C81  C82
c    C90  C91  C92
c
            IF(IS.GT.1.AND.IS.LT.NSROW) THEN
c if disprc=0, skip redefinition of cavg terms
		    IF(DISPCR(JS,IS,KS).NE.0.0) THEN
                CMULT=1.0
	          C11=CAVG(JS,IS-1,KS)
	          C12=CAVG(JS+1,IS-1,KS)
	          C31=CAVG(JS,IS+1,KS)
	          C32=CAVG(JS+1,IS+1,KS)
c chfb
               IF(ONCHFB.EQ.1) THEN
c  check two cells above, if either hfb then use 
c   current row + alter length
	          C22=CAVG(JS+1,IS,KS)
	          IF(LOCHFB(JS,IS-1,KS,2).EQ.1.OR.
     *		     LOCHFB(JS+1,IS-1,KS,2).EQ.1) THEN
                      C11=C2
                      C12=C22
c                     account for new delta R here
			        CMULT=2.0
     	          END IF
c  then check two cells below, if either hfb then use 
c   current row + alter length
	          IF(LOCHFB(JS,IS,KS,2).EQ.1.OR.
     *		     LOCHFB(JS+1,IS,KS,2).EQ.1) THEN
                      C31=C2
                      C32=C22
c                     account for new delta R here
			        CMULT=2.0
	          END IF
c  if both above and below, grad should=0 when we do difference below
c chfb end
               END IF
C MOCWT
 	         IF(PTWTON.EQ.1) THEN
	          IF(NPCELL(JS,IS+1,KS).LT.1) C31=C2
	          IF(NPCELL(JS+1,IS+1,KS).LT.1) C32=CAVG(JS+1,IS,KS)
	          IF(NPCELL(JS,IS-1,KS).LT.1) C11=C2
	          IF(NPCELL(JS+1,IS-1,KS).LT.1) C12=CAVG(JS+1,IS,KS)
	         END IF
               DSPXX2=DSPXX2+DISPCR(JS,IS,KS)*(C31+C32-C11-C12)*CMULT
c disp ne 0 end
	        END IF
	       END IF
CFWD  FLUX DUE TO DXZ
c            IF(KS.GT.1.AND.KS.LT.NSLAY)
c     *        DSPXX2=DSPXX2+DISPCL(JS,IS,KS)*
c     *         (CAVG(JS,IS,KS+1)+CAVG(JS+1,IS,KS+1)
c     *         -CAVG(JS,IS,KS-1)-CAVG(JS+1,IS,KS-1))
            IF(KS.GT.1.AND.KS.LT.NSLAY) THEN
c disp fix  the check below should eradicate divide by zero thickness 
c disp fix    due to check in DSP6FM routine for IBOUND of cells used here
             IF(DISPCL(JS,IS,KS).NE.0.0) THEN
	        C81=CAVG(JS,IS,KS+1)
	        C82=CAVG(JS+1,IS,KS+1)
	        C51=CAVG(JS,IS,KS-1)
	        C52=CAVG(JS+1,IS,KS-1)
	        C22=CAVG(JS+1,IS,KS)
	        B2=THCK(JS,IS,KS)
	        B81=THCK(JS,IS,KS+1)
	        B82=THCK(JS+1,IS,KS+1)
	        B51=THCK(JS,IS,KS-1)
	        B52=THCK(JS+1,IS,KS-1)
	        B22=THCK(JS+1,IS,KS)
C MOCWT
 	        IF(PTWTON.EQ.1) THEN
	         IF(NPCELL(JS,IS,KS+1).LT.1) C81=C2
	         IF(NPCELL(JS+1,IS,KS+1).LT.1) C82=CAVG(JS+1,IS,KS)
	         IF(NPCELL(JS,IS,KS-1).LT.1) C51=C2
	         IF(NPCELL(JS+1,IS,KS-1).LT.1) C52=CAVG(JS+1,IS,KS)
	        END IF
c orig              DSPXX2=DSPXX2+DISPCL(JS,IS,KS)*(C81+C82-C51-C52)
c disp fix
              DSPXX2=DSPXX2+DISPCL(JS,IS,KS)*
     *        (((C82-C22)/(B82+B22))+((C22-C52)/(B22+B52))+
     *         ((C81-C2)/(B81+B2))+((C2-C51)/(B2+B51)))
	       END IF
            END IF
         END IF
C
CLMT SAVE FLUX FOR USE IN ADJUSTMENT LOOP
C
C    DXXCUR IS THE DXX FLUX ACROSS THE FACE INTO/OUT OF CURRENT CELL
C    DXXFOR IS THE DXX FLUX ACROSS THE FACE INTO/OUT OF FORWARD CELL
C    NOTE THE SIGN IS THE SAME FOR BOTH FLUXES, SO WHEN APPLIED TO CNCNC
C    OPPOSITE SIGNS MUST BE USED
C
         DXXCUR(JS,IS,KS)=DSPXX2*TRFPRB
         DXXFOR(JS,IS,KS)=DSPXX2*tPRBIN2
cgzh debug output
c      if(js.eq.1.and.is.eq.1) then
c        write(iouts,*) 'first loop, cell=', js,is,ks
c        write(iouts,*) 'dxxcur,dxxfor',DXXCUR(JS,IS,KS),DXXFOR(JS,IS,KS)
c        write(71,*) 'dxxcur,dxxfor',DXXCUR(JS,IS,KS),DXXFOR(JS,IS,KS)
c        write(iouts,*) 'conc,cavg,sumwt',conc(JS,IS,KS),cavg(JS,IS,KS),
c     * sumwt(JS,IS,KS)
c      end if
c      if(js.eq.4.and.is.eq.1) then
c        write(iouts,*) 'first loop, cell=', js,is,ks
c        write(iouts,*) 'dxxcur,dxxfor',DXXCUR(JS,IS,KS),DXXFOR(JS,IS,KS)
c        write(iouts,*) 'conc,cavg,sumwt',conc(JS,IS,KS),cavg(JS,IS,KS),
c     * sumwt(JS,IS,KS)
c      end if
c      if(js.eq.24.and.ks.eq.1) then
c        write(iouts,*) 'first loop, cell=', js,is,ks
c        write(iouts,*) 'dxxcur,dxxfor',DXXCUR(JS,IS,KS),DXXFOR(JS,IS,KS)
c      end if
CFM  ADD FLUX TO FORWARD ADJACENT BLOCK
cgzh debug we do this later...
c         CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+DXXCUR(JS,IS,KS)
c         CNCNC(JS+1,IS,KS)=CNCNC(JS+1,IS,KS)-DXXFOR(JS,IS,KS)
cgzh debug output
c      if(js.eq.138.and.ks.eq.2.and.imov.gt.2670) then
c        write(iouts,*) 'dspxx2 at js   = ',DSPXX2*TRFPRB
c      end if
c      if(js+1.eq.138.and.ks.eq.2.and.imov.gt.2670) then
c        write(iouts,*) 'dspxx2 at js+1 = ',(-DSPXX2*timv*tPRBIN2/
c     *     rf(js+1,is,ks))
c      end if
      END IF
      END IF
C
CFM  DISP FLUX IN FORWARD ROW DIRECTION
 15   IF(IS.LT.NSROW) THEN
      IF(IBOUND(J,I+1,K).NE.0) THEN
CMOCWT  SKIP TO LAYER DIRECTION IF NO PARTICLES IN FORWARD CELL
       IF(PTWTON.EQ.1.AND.NPCELL(JS,IS+1,KS).LT.1) GO TO 17
        tPRBIN2=timv/(POR(JS,IS+1,KS)*THCK(JS,IS+1,KS)*rf(js,is+1,ks))
         DSPYY2=DISPRR(JS,IS,KS)*(CAVG(JS,IS+1,KS)-C2)
cgzh debug output
c      if(js.eq.1.and.is.eq.2.and.imov.eq.5) then
c       write(iouts,*) 'DSPYY2 1,2',DSPYY2
c       write(iouts,*) 'DISPRR,CAVG,C2 2,1',DISPRR(JS,IS,KS),
c     *CAVG(JS,IS+1,KS),C2
c	end if
         IF(NODISP.NE.1) THEN
c orig
c            IF(JS.GT.1.AND.JS.LT.NSCOL)
c     *        DSPYY2=DSPYY2+DISPRC(JS,IS,KS)*
c     *         (CAVG(JS+1,IS,KS)+CAVG(JS+1,IS+1,KS)
c     *         -CAVG(JS-1,IS,KS)-CAVG(JS-1,IS+1,KS))
C
            IF(JS.GT.1.AND.JS.LT.NSCOL) THEN
		    IF(DISPRC(JS,IS,KS).NE.0.0) THEN
 			 CMULT=1.0
 	         C22=CAVG(JS+1,IS,KS)
	         C32=CAVG(JS+1,IS+1,KS)
	         C20=CAVG(JS-1,IS,KS)
	         C30=CAVG(JS-1,IS+1,KS)
c chfb
               IF (ONCHFB.EQ.1) THEN
	          C31=CAVG(JS,IS+1,KS)
	          IF(LOCHFB(JS-1,IS,KS,1).EQ.1.OR.
     *		     LOCHFB(JS-1,IS+1,KS,1).EQ.1) THEN
                      C20=C2
                      C30=C31
			        CMULT=2.0
	          END IF
	          IF(LOCHFB(JS,IS,KS,1).EQ.1.OR.
     *		     LOCHFB(JS,IS+1,KS,1).EQ.1) THEN
                      C22=C2
                      C32=C31
			        CMULT=2.0
	          END IF
c chfb end
               END IF
C
	         IF(PTWTON.EQ.1) THEN
	          IF(NPCELL(JS+1,IS,KS).LT.1) C22=C2
	          IF(NPCELL(JS+1,IS+1,KS).LT.1) C32=CAVG(JS,IS+1,KS)
	          IF(NPCELL(JS-1,IS,KS).LT.1) C20=C2
	          IF(NPCELL(JS-1,IS+1,KS).LT.1) C30=CAVG(JS,IS+1,KS)
	         END IF
C
               DSPYY2=DSPYY2+DISPRC(JS,IS,KS)*(C22+C32-C20-C30)*CMULT
c disp ne 0 end
              END IF
	      END IF
C
c orig
c            IF(KS.GT.1.AND.KS.LT.NSLAY)
c     *        DSPYY2=DSPYY2+DISPRL(JS,IS,KS)*
c     *         (CAVG(JS,IS,KS+1)+CAVG(JS,IS+1,KS+1)
c     *         -CAVG(JS,IS,KS-1)-CAVG(JS,IS+1,KS-1))
            IF(KS.GT.1.AND.KS.LT.NSLAY) THEN
c disp fix  the check below should eradicate divide by zero thickness 
c disp fix    due to check in DSP6FM routine for IBOUND of cells used here
             IF(DISPRL(JS,IS,KS).NE.0.0) THEN
	        C81=CAVG(JS,IS,KS+1)
	        C91=CAVG(JS,IS+1,KS+1)
	        C51=CAVG(JS,IS,KS-1)
	        C61=CAVG(JS,IS+1,KS-1)
	        C31=CAVG(JS,IS+1,KS)
	        B2=THCK(JS,IS,KS)
	        B81=THCK(JS,IS,KS+1)
	        B91=THCK(JS,IS+1,KS+1)
	        B51=THCK(JS,IS,KS-1)
	        B61=THCK(JS,IS+1,KS-1)
	        B31=THCK(JS,IS+1,KS)
	        IF(PTWTON.EQ.1) THEN
	          IF(NPCELL(JS,IS,KS+1).LT.1) C81=C2
	          IF(NPCELL(JS,IS+1,KS+1).LT.1) C91=CAVG(JS,IS+1,KS)
	          IF(NPCELL(JS,IS,KS-1).LT.1) C51=C2
	          IF(NPCELL(JS,IS+1,KS-1).LT.1) C61=CAVG(JS,IS+1,KS)
	        END IF
c orig              DSPYY2=DSPYY2+DISPRL(JS,IS,KS)*(C3+C4-C5-C6)
c disp fix
              DSPYY2=DSPYY2+DISPRL(JS,IS,KS)*
     *        (((C91-C31)/(B91+B31))+((C31-C61)/(B31+B61))+
     *         ((C81-C2)/(B81+B2))+((C2-C51)/(B2+B51)))
	       END IF
	      END IF
         END IF
C
CLMT SAVE FLUX FOR USE IN ADJUSTMENT LOOP
C
C    DYYCUR IS THE DYY FLUX ACROSS THE FACE INTO/OUT OF CURRENT CELL
C    DYYFOR IS THE DYY FLUX ACROSS THE FACE INTO/OUT OF FORWARD CELL
C
         DYYCUR(JS,IS,KS)=DSPYY2*TRFPRB
         DYYFOR(JS,IS,KS)=DSPYY2*tPRBIN2
CFM  ADD FLUX TO FORWARD ADJACENT BLOCK
c	   CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+DYYCUR(JS,IS,KS)
c         CNCNC(JS,IS+1,KS)=CNCNC(JS,IS+1,KS)-DYYFOR(JS,IS,KS) 
cgzh debug output
cgzh debug output
c      if(js.eq.1.and.is.eq.1) then
c        write(iouts,*) 'first loop, cell=', js,is,ks
c        write(iouts,*) 'dyycur,dyyfor',DyyCUR(JS,IS,KS),DyyFOR(JS,IS,KS)
c      end if
c      if(js.eq.138.and.ks.eq.2.and.imov.gt.2670) then
c        write(iouts,*) 'dspyy2 at js   = ',DSPYY2*TRFPRB
c      end if
c      if(js+1.eq.138.and.ks.eq.2.and.imov.gt.2670) then
c        write(iouts,*) 'dspyy2 at is+1 = ',(-DSPYY2*tPRBIN2)
c      end if
      END IF
      END IF
CFM  Z TERM
CFM  THIS VERSION INCLUDES CORRECTION FOR THICKNESS VARIABILITY 
CFM   WHEN INTERPOLATING FOR GRADZ
  17  IF(KS.LT.NSLAY) THEN
      IF(IBOUND(J,I,K+1).NE.0) THEN
CMOCWT  SKIP TO ADJUSTMENT FACTOR CALCULATION IF NO PARTICLES IN FORWARD CELL
       IF(PTWTON.EQ.1.AND.NPCELL(JS,IS,KS+1).LT.1) GO TO 19
         TRFPRB2=TIMV/(RF(js,is,KS+1)*POR(JS,IS,KS+1)*THCK(JS,IS,KS+1))
         DSPZZ2=DISPLL(JS,IS,KS)*(CAVG(JS,IS,KS+1)-C2)
         IF(NODISP.NE.1) THEN
c disp fix
            SUMTHCK=THCK(JS,IS,KS)+THCK(JS,IS,KS+1)
            THCKWT1=THCK(JS,IS,KS)/SUMTHCK
            THCKWT2=THCK(JS,IS,KS+1)/SUMTHCK
c orig
c            IF(JS.GT.1.AND.JS.LT.NSCOL)
c     *        DSPZZ2=DSPZZ2+DISPLC(JS,IS,KS)*
c     *         (CAVG(JS+1,IS,KS)+CAVG(JS+1,IS,KS+1)
c     *         -CAVG(JS-1,IS,KS)-CAVG(JS-1,IS,KS+1))
            IF(JS.GT.1.AND.JS.LT.NSCOL) THEN
c if displc=0, skip redefinition of cavg terms
		    IF(DISPLC(JS,IS,KS).NE.0.0) THEN
 			 CMULT=1.0
	         C22=CAVG(JS+1,IS,KS)
	         C82=CAVG(JS+1,IS,KS+1)
	         C20=CAVG(JS-1,IS,KS)
	         C80=CAVG(JS-1,IS,KS+1)
c chfb
               IF(ONCHFB.EQ.1) THEN
	            C81=CAVG(JS,IS,KS+1)
	          IF(LOCHFB(JS-1,IS,KS,1).EQ.1.OR.
     *		     LOCHFB(JS-1,IS,KS+1,1).EQ.1) THEN
                      C20=C2
                      C80=C81
			        CMULT=2.0
	          END IF
	          IF(LOCHFB(JS,IS,KS,1).EQ.1.OR.
     *		     LOCHFB(JS,IS,KS+1,1).EQ.1) THEN
                      C22=C2
                      C82=C81
			        CMULT=2.0
	          END IF
c chfb end
	         END IF
	         IF(PTWTON.EQ.1) THEN
	          IF(NPCELL(JS+1,IS,KS).LT.1) C22=C2
	          IF(NPCELL(JS+1,IS,KS+1).LT.1) C82=CAVG(JS,IS,KS+1)
	          IF(NPCELL(JS-1,IS,KS).LT.1) C20=C2
	          IF(NPCELL(JS-1,IS,KS+1).LT.1) C80=CAVG(JS,IS,KS+1)
	         END IF
c orig              DSPZZ2=DSPZZ2+DISPLC(JS,IS,KS)*(C3+C4-C5-C6)
c disp fix
               DSPZZ2=DSPZZ2+DISPLC(JS,IS,KS)*(THCKWT1*(C82-C80)+
     *           THCKWT2*(C22-C20))*CMULT
c disp ne 0 end
              END IF
	      END IF
c            IF(IS.GT.1.AND.IS.LT.NSROW)
c     *        DSPZZ2=DSPZZ2+DISPLR(JS,IS,KS)*
c     *         (CAVG(JS,IS+1,KS)+CAVG(JS,IS+1,KS+1)
c     *         -CAVG(JS,IS-1,KS)-CAVG(JS,IS-1,KS+1))
            IF(IS.GT.1.AND.IS.LT.NSROW) THEN
c if displr=0, skip redefinition of cavg terms
		    IF(DISPLR(JS,IS,KS).NE.0.0) THEN
			 CMULT=1.0
	         C31=CAVG(JS,IS+1,KS)
	         C91=CAVG(JS,IS+1,KS+1)
	         C11=CAVG(JS,IS-1,KS)
	         C71=CAVG(JS,IS-1,KS+1)
c chfb
               IF(ONCHFB.EQ.1) THEN
	          C81=CAVG(JS,IS,KS+1)
cgzh 4/25/07 bug fix
c	          IF(LOCHFB(JS-1,IS,KS,2).EQ.1.OR.
c     *		     LOCHFB(JS-1,IS,KS+1,2).EQ.1) THEN
	          IF(LOCHFB(JS,IS-1,KS,2).EQ.1.OR.
     *		     LOCHFB(JS,IS-1,KS+1,2).EQ.1) THEN
                      C31=C2
                      C91=C81
			        CMULT=2.0
	          END IF
	          IF(LOCHFB(JS,IS,KS,2).EQ.1.OR.
     *		     LOCHFB(JS,IS,KS+1,2).EQ.1) THEN
                      C11=C2
                      C81=C81
			        CMULT=2.0
	          END IF
c chfb end
               END IF 
	         IF(PTWTON.EQ.1) THEN
	          IF(NPCELL(JS,IS+1,KS).LT.1) C31=C2
	          IF(NPCELL(JS,IS+1,KS+1).LT.1) C91=CAVG(JS,IS,KS+1)
	          IF(NPCELL(JS,IS-1,KS).LT.1) C11=C2
	          IF(NPCELL(JS,IS-1,KS+1).LT.1) C71=CAVG(JS,IS,KS+1)
	         END IF
c orig              DSPZZ2=DSPZZ2+DISPLR(JS,IS,KS)*(C3+C4-C5-C6)
c disp fix
               DSPZZ2=DSPZZ2+DISPLR(JS,IS,KS)*(THCKWT1*(C91-C71)+
     *           THCKWT2*(C31-C11))*CMULT
c disp ne 0 end
	        END IF
	      END IF
         END IF
C
CLMT SAVE FLUX FOR USE IN ADJUSTMENT LOOP
C
C    DZZCUR IS THE DZZ FLUX ACROSS THE FACE INTO/OUT OF CURRENT CELL
C    DZZFOR IS THE DZZ FLUX ACROSS THE FACE INTO/OUT OF FORWARD CELL
C
c orig         DZZCUR(JS,IS,KS)=DSPZZ2*TOPRF
c disp fix
         DZZCUR(JS,IS,KS)=DSPZZ2*TRFPRB
cgzh dzz
cgzh orig         DZZFOR(JS,IS,KS)=DSPZZ2*TIMV*RFPOR2*
C multiply by ratio of thicknes to balance mass in forward cell
cgzh don't do this mult any more, thinkness ratios already handled for dz terms
cgzh (note that now, all dzz terms * trfprb or ___b2, this accounts for diff)
c         DZZFOR(JS,IS,KS)=DSPZZ2*TRFPRB2*
c     *  (THCK(JS,IS,KS)/THCK(JS,IS,KS+1))
cgzh dzz
         DZZFOR(JS,IS,KS)=DSPZZ2*TRFPRB2
cgzh debug output
c      if(js.eq.1.and.ks.eq.2) then
c        write(iouts,*) 'first loop, cell=', js,is,ks
c        write(iouts,*) 'dzzcur,dzzfor',DzzCUR(JS,IS,KS),DzzFOR(JS,IS,KS)
c      end if
CFM  ADD FLUX TO FORWARD ADJACENT BLOCK
c	   CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+DZZCUR(JS,IS,KS)
c         CNCNC(JS,IS,KS+1)=CNCNC(JS,IS,KS+1)-DZZFOR(JS,IS,KS) 
cgzh debug output
c      if(js.eq.138.and.ks.eq.2.and.imov.gt.2670) then
c        write(iouts,*) 'dspzz2 at js   = ',(DSPZZ2*timv/
c     *     (POR(JS,IS,KS)*rf(js,is,ks)))
c      end if
c      if(js+1.eq.138.and.ks.eq.2.and.imov.gt.2670) then
c        write(iouts,*) 'dspzz2 at is+1 = ',(-DSPZZ2*TIMV*RFPOR2*
c     *  (THCK(JS,IS,KS)/THCK(JS,IS,KS+1)))
c      end if
      END IF
      END IF
C
CLMT  CALCULATE ADJUSTMENT FACTOR FOR CURRENT CELL
C
C     INITIALIZE ADJUSTMENT FACTOR
 19   ADJFCT(JS,IS,KS)=1.0
cgzh debug skipping adjustment
c      if(js.eq.1.and.is.eq.1.and.ks.eq.1)write(*,*) 'skipping adj'
c      go to 79
C     IF REDUCTION IN MASS
cgzh debug
c this was taken out since adjusted face flows may not be reflected in CNCNC
c for the CNEW check below
c      IF(CNCNC(JS,IS,KS).LT.0.0) THEN
c
C     CALCULATE ADJUSTED BACKWARD FACE FLOWS
C     THIS RESULTS IN A MORE PRECISE CALCULATION OF THE ADJUSTMENT FACTOR
      BACKFX=0.0
      BACKFY=0.0
      BACKFZ=0.0
      IF(JS.GT.1) THEN
        IF(ADJFCT(JS-1,IS,KS).NE.-1.0) 
cgzh debug taking out * by 1.0 
c     *   BACKFX=DXXFOR(JS-1,IS,KS)*ADJFCT(JS-1,IS,KS)	       
     *   BACKFX=DXXFOR(JS-1,IS,KS)       
      END IF
      IF(IS.GT.1) THEN
        IF(ADJFCT(JS,IS-1,KS).NE.-1.0) 
cgzh debug changed to YY bug fix 
     *   BACKFY=DYYFOR(JS,IS-1,KS)	       
      END IF
      IF(KS.GT.1) THEN
        IF(ADJFCT(JS,IS,KS-1).NE.-1.0) 
cgzh debug changed to ZZ bug fix 
     *   BACKFZ=DZZFOR(JS,IS,KS-1)	       
      END IF
C
C     CALCULATE APPROXIMATE CHANGE IN CONC AT CELL USING ADJUSTED BACKWARD
C       FACE FLUXES
      TCNCNC=-BACKFX-BACKFY-BACKFZ+DXXCUR(JS,IS,KS)+DYYCUR(JS,IS,KS)
     *    +DZZCUR(JS,IS,KS)
C       CALCULATE (POTENTIAL) NEW CONCENTRATION
C       HERE, CONC ARRAY HAS NEW CONCENTRATION AT CELLS DUE TO 
C         ADVECTION OF PARTICLES 
cgzh debug
C       CHECK FOR NEGATIVE CONCENTRATION (UNWEIGHTED)
        IF(PTWTON.EQ.0) THEN
         CNEW=CONC(JS,IS,KS)+TCNCNC
        ELSE
C       CHECK FOR NEGATIVE MASS (WEIGHTED)
         CNEW=(CONC(JS,IS,KS)*SUMWT(JS,IS,KS))+(TCNCNC*CELVOL(JS,IS,KS))
        END IF
C       IF THIS RESULTS IN A NEGATIVE CONC OR MASS, SET ADJUSTMENT FACTOR 
        IF(CNEW.LT.0.0) THEN
cgzh debug output
c      write(iouts,*) 'cnew<0:'
c	write(iouts,*) 'SUMWT(JS,IS,KS),TCNCNC,CELVOL(JS,IS,KS)'
c	write(iouts,*) SUMWT(JS,IS,KS),TCNCNC,CELVOL(JS,IS,KS)
C         IF CONCENTRATION IS ALREADY ZERO, SET FACTOR TO ZERO TO
C           AVOID REMOVING EVEN MORE MASS
          IF(CONC(JS,IS,KS).LE.0.0) THEN
            ADJFCT(JS,IS,KS)=0.0
cgzh debug output
c            write(iouts,*) 'Adjustment set to zero'
          ELSE
C         FACTOR IS RATIO OF CONC AVAILABLE TO PROPOSED CONC CHANGE;
C           WILL BE MULTIPLIED BY CNCNC IN NEXT LOOP
C         E.G. CONC=1, CNCNC=-3, ADJFCT = 1/3
C         NEW CNCNC=-3*(1/3)=-1  
cgzh debug
c changed this to disallow negative factors, which could "swing" the
c   calculations too much ("unrealistically")
c
c orig          ADJFCT(JS,IS,KS)=CONC(JS,IS,KS)/ABS(CNCNC(JS,IS,KS))
c            ADJFCT(JS,IS,KS)=ABS(CONC(JS,IS,KS)/CNCNC(JS,IS,KS))
C RBW begin change
!            ADJFCT(JS,IS,KS)=ABS(CONC(JS,IS,KS)/TCNCNC)
            if (TCNCNC.ne.0) then
              ADJFCT(JS,IS,KS)=ABS(CONC(JS,IS,KS)/TCNCNC)
             else
               ADJFCT(JS,IS,KS)=0
             endif
C RBW end change
c            write(iouts,*) 'Adjustment set to ',ADJFCT(JS,IS,KS)
cgzh debug
C           IF WEIGHT ON PARTICLES IS LESS THAN CELL VOLUME, REDUCE
C           ADJUSTMENT FACTOR TO AVOID REMOVING MORE MASS THAN EXISTS
C           ON PARTICLES
cgzh  check volume here...should we be checking mass instead?
            IF(PTWTON.GT.0) THEN
              IF(SUMWT(JS,IS,KS).LT.CELVOL(JS,IS,KS)) THEN
c      write(iouts,*) 'sumwt<celvol: adj precorrection=',ADJFCT(JS,IS,KS) 
                ADJFCT(JS,IS,KS)=
     *          ADJFCT(JS,IS,KS)*(SUMWT(JS,IS,KS)/CELVOL(JS,IS,KS))
c      write(iouts,*) 'sumwt<celvol: adj corrected=',ADJFCT(JS,IS,KS) 
              END IF
            END IF
          END IF
cgzh debug output
c          write(iouts,*) 'adjustment needed at js,is,ks',JS,IS,KS 
c          write(iouts,*) 'conc,cncnc',CONC(JS,IS,KS),CNCNC(JS,IS,KS)
c          write(iouts,*) 'cnew,tcncnc',cnew,tcncnc
c          write(iouts,*) 'adjustment factor=',ADJFCT(JS,IS,KS) 
c          write(iouts,*) BACKFX,BACKFY,BACKFZ,DXXCUR(JS,IS,KS),
c     *    DYYCUR(JS,IS,KS),DZZCUR(JS,IS,KS)
c          write(iouts,*) 
        END IF
cgzh debug CNCNC lt 0 check
c      END IF
cgzh debug adj off
c      ADJFCT(JS,IS,KS)=1.0
   20 CONTINUE
cgzh debug skipping adjustment
c      return
C     ****************************************************************
C
CLMT  LIMIT DISP BY AVAILABLE MASS 
C     LOOP OVER CELLS TO UPDATE ESTIMATE OF ADJUSTMENT FACTOR
C
      DO 40 KS=1,NSLAY
      DO 40 IS=1,NSROW
      DO 40 JS=1,NSCOL
C  SKIP IF INACTIVE CELL
      IF(ADJFCT(JS,IS,KS).EQ.-1.D0) GO TO 40
C  SKIP IF NO PARTICLES IN CELL
      IF(PTWTON.EQ.1.AND.NPCELL(JS,IS,KS).LT.1) GO TO 40
      F1F=0.0
      F1B=0.0
      F2F=0.0
      F2B=0.0
      F3F=0.0
      F3B=0.0
C
C  *************************
C  START DXX TERM ADJUSTMENT
C  *************************
C     FORWARD DIRECTION: SKIP CELLS IN LAST COLUMN
      IF(JS.LT.NSCOL) THEN
C      SKIP CELL IF FORWARD CELL IS NO FLOW
       IF(ADJFCT(JS+1,IS,KS).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN FORWARD CELL (MOCWT)
CMOCWT
cgzh debug -> fix for unweighted routine
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS+1,IS,KS).GT.0)) THEN
C
C         THE FLUX CALCULATED AT CURRENT CELL WAS A CONTRIBUTOR TO BOTH CELLS          
C          SO IF ADJUSTMENT IS REQUIRED AT EITHER CURRENT OR FORWARD CELL,
C          AN ADJUSTMENT IS REQUIRED
C
          F1F=DXXCUR(JS,IS,KS)
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS+1,IS,KS).LT.1.0) THEN
C
C          IN CASE BOTH CURRENT AND FORWARD CELL NEED ADJUSTMENT,
C           NEED TO CHECK LARGEST ADJUSTMENT (THAT IS, THE SMALLEST FACTOR)
C           REMEMBER WE WILL MULTIPLY THE DISPERSIVE FLUX (CNCNC) BY THE FACTOR
C           MIN WORKS BECAUSE DEFAULT VALUE OF ADJFCT IS 1.0, AND NO FLOW
C           CELLS (WHERE ADJFCT=-1.0) HAVE BEEN SKIPPED
C
              TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS+1,IS,KS))
              F1F=DXXCUR(JS,IS,KS)*TEMFCT
cgzh debug
c             write(*,*) 'AAAAA js is ks ',js,is,ks
          END IF
C
         END IF
       END IF
      END IF
C
C     BACKWARD DIRECTION: SKIP CELLS IN FIRST COLUMN
      IF(JS.GT.1) THEN
C      SKIP CELL IF BACKWARD CELL IS NO FLOW
       IF(ADJFCT(JS-1,IS,KS).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN BACKWARD CELL (MOCWT)
CMOCWT
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS-1,IS,KS).GT.0)) THEN
C
          F1B=DXXFOR(JS-1,IS,KS)
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS-1,IS,KS).LT.1.0) THEN
              TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS-1,IS,KS))
              F1B=DXXFOR(JS-1,IS,KS)*TEMFCT
cgzh debug
c             write(*,*) 'BBBBB js is ks ',js,is,ks
          END IF
C
         END IF
       END IF
      END IF
C
C  *************************
C  START DYY TERM ADJUSTMENT
C  *************************
C
C     FORWARD DIRECTION: SKIP CELLS IN LAST COLUMN
      IF(IS.LT.NSROW) THEN
C      SKIP CELL IF FORWARD CELL IS NO FLOW
       IF(ADJFCT(JS,IS+1,KS).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN FORWARD CELL (MOCWT)
CMOCWT
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS,IS+1,KS).GT.0)) THEN
C
          F2F=DYYCUR(JS,IS,KS)
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS,IS+1,KS).LT.1.0) THEN
              TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS,IS+1,KS))
              F2F=DYYCUR(JS,IS,KS)*TEMFCT
          END IF
C
         END IF
       END IF
      END IF
C
C     BACKWARD DIRECTION: SKIP CELLS IN FIRST COLUMN
      IF(IS.GT.1) THEN
C      SKIP CELL IF BACKWARD CELL IS NO FLOW
       IF(ADJFCT(JS,IS-1,KS).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN BACKWARD CELL (MOCWT)
CMOCWT
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS,IS-1,KS).GT.0)) THEN
C
          F2B=DYYFOR(JS,IS-1,KS)
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS,IS-1,KS).LT.1.0) THEN
              TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS,IS-1,KS))
              F2B=DYYFOR(JS,IS-1,KS)*TEMFCT
          END IF
C
         END IF
       END IF
      END IF
C
C
C  *************************
C  START DZZ TERM ADJUSTMENT
C  *************************
C
C     FORWARD DIRECTION: SKIP CELLS IN LAST COLUMN
      IF(KS.LT.NSLAY) THEN
C      SKIP CELL IF FORWARD CELL IS NO FLOW
       IF(ADJFCT(JS,IS,KS+1).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN FORWARD CELL (MOCWT)
CMOCWT
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS,IS,KS+1).GT.0)) THEN
C
          F3F=DZZCUR(JS,IS,KS)
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS,IS,KS+1).LT.1.0) THEN
              TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS,IS,KS+1))
              F3F=DZZCUR(JS,IS,KS)*TEMFCT
          END IF
C
         END IF
       END IF
      END IF
C
C     BACKWARD DIRECTION: SKIP CELLS IN FIRST COLUMN
      IF(KS.GT.1) THEN
C      SKIP CELL IF BACKWARD CELL IS NO FLOW
       IF(ADJFCT(JS,IS,KS-1).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN BACKWARD CELL (MOCWT)
CMOCWT
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS,IS,KS-1).GT.0)) THEN
C
          F3B=DZZFOR(JS,IS,KS-1)
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS,IS,KS-1).LT.1.0) THEN
              TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS,IS,KS-1))
              F3B=DZZFOR(JS,IS,KS-1)*TEMFCT
          END IF
C
         END IF
       END IF
      END IF
C
        TCNCNC=F1F-F1B+F2F-F2B+F3F-F3B
c        CNEW=(CONC(JS,IS,KS)*SUMWT(JS,IS,KS))+(TCNCNC*CELVOL(JS,IS,KS))
cgzh ?? did not include check by mass here, because tcncnc already has adjustment?
        CNEW=CONC(JS,IS,KS)+TCNCNC
C       IF THIS WOULD RESULT IN A NEGATIVE CONCENTRATION, SET ADJUSTMENT FACTOR 
        IF(CNEW.LT.0.0) THEN
C         IF CONCENTRATION IS ALREADY ZERO, SET FACTOR TO ZERO TO
C           AVOID REMOVING EVEN MORE MASS
          IF(CONC(JS,IS,KS).LE.0.0) THEN
            ADJFCT(JS,IS,KS)=0.0
          ELSE
C         FACTOR IS RATIO OF CONC AVAILABLE TO PROPOSED CONC CHANGE;
C           WILL BE MULTIPLIED BY CNCNC IN NEXT LOOP
C         E.G. CONC=1, CNCNC=-3, ADJFCT = 1/3
C         NEW CNCNC=-3*(1/3)=-1  
            ADJFCT(JS,IS,KS)=ABS(CONC(JS,IS,KS)/TCNCNC)
C
C           IF WEIGHT ON PARTICLES IS LESS THAN CELL VOLUME, REDUCE
C           ADJUSTMENT FACTOR TO AVOID REMOVING MORE MASS THAN EXISTS
C           ON PARTICLES
            IF(PTWTON.GT.0) THEN
              IF(SUMWT(JS,IS,KS).LT.CELVOL(JS,IS,KS)) THEN
		      ADJFCT(JS,IS,KS)=
     *             ADJFCT(JS,IS,KS)*(SUMWT(JS,IS,KS)/CELVOL(JS,IS,KS))
              END IF
            END IF
          END IF
        END IF
C  END LOOP FOR ADJUSTMENT OF ADJFCT
cgzh debug adj off
c      ADJFCT(JS,IS,KS)=1.0
   40 CONTINUE
C
C     ****************************************************************
C
CLMT  LIMIT DISP BY AVAILABLE MASS 
C
C     LOOP OVER CELLS TO ADJUST CHANGE IN CONC BY ADJUSTMENT FACTOR
cgzh debug continue
 79   continue
C
      DO 80 KS=1,NSLAY
      DO 80 IS=1,NSROW
      DO 80 JS=1,NSCOL
C  SKIP IF INACTIVE CELL
cgzh debug      IF(ADJFCT(JS,IS,KS).LT.0.0) GO TO 80
      IF(ADJFCT(JS,IS,KS).EQ.-1.D0) GO TO 80
C  SKIP IF NO PARTICLES IN CELL
      IF(PTWTON.EQ.1.AND.NPCELL(JS,IS,KS).LT.1) GO TO 80
C
C  *************************
C  START DXX TERM ADJUSTMENT
C  *************************
C     SKIP CELLS IN LAST COLUMN
      IF(JS.LT.NSCOL) THEN
C      SKIP CELL IF FORWARD CELL IS NO FLOW
       IF(ADJFCT(JS+1,IS,KS).GT.-1.0) THEN
C        SKIP IF NO PARTICLES IN FORWARD CELL (MOCWT)
CMOCWT
cgzh debug -> fix for unweighted routine
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS+1,IS,KS).GT.0)) THEN
C
C         THE FLUX CALCULATED AT CURRENT CELL WAS A CONTRIBUTOR TO BOTH CELLS          
C          SO IF ADJUSTMENT IS REQUIRED AT EITHER CURRENT OR FORWARD CELL,
C          AN ADJUSTMENT IS REQUIRED
C
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS+1,IS,KS).LT.1.0) THEN
C
C          IN CASE BOTH CURRENT AND FORWARD CELL NEED ADJUSTMENT,
C           NEED TO CHECK LARGEST ADJUSTMENT (THAT IS, THE SMALLEST FACTOR)
C           REMEMBER WE WILL MULTIPLY THE DISPERSIVE FLUX (CNCNC) BY THE FACTOR
C           MIN WORKS BECAUSE DEFAULT VALUE OF ADJFCT IS 1.0, AND NO FLOW
C           CELLS (WHERE ADJFCT=0.0) HAVE BEEN SKIPPED
C
           TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS+1,IS,KS))
C
C          AT CURRENT CELL, SUBTRACT THE DXX FLUX FROM CNCNC
C            (1) CNCNC=CNCNC-DXXFLX
C          THEN ADD BACK IN THE DXX FLUX TIMES THE ADJUSTMENT FACTOR
C            (2) CNCNC=CNCNC+(DXXFLX*FACTOR)
C          THIS SIMPLIFIES TO
C            (3) CNCNC=CNCNC+DXXFLX(FACTOR-1.0)       
C          DXXCUR HAS THE DXX FLUX ACROSS THE FACE INTO/OUT OF CURRENT CELL
C
cgzh debug output
c          write(iouts,*) 'face flux adjustment at js,is,ks',JS,IS,KS 
c          write(iouts,*) 'dxx cncnc BEFORE adjustment at current cell' 
c          write(iouts,*) 'DXXCUR= ',DXXCUR(JS,IS,KS)
c
cgzh debug  new way: just adjust the face flux
           DXXCUR(JS,IS,KS)=DXXCUR(JS,IS,KS)*TEMFCT
c           CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+
c     *                       DXXCUR(JS,IS,KS)*(TEMFCT-1.D0)
cgzh debug output
c          write(iouts,*) 'dxx cncnc AFTER adjustment at current cell' 
c          write(iouts,*) 'DXXCUR= ',DXXCUR(JS,IS,KS)
c          write(iouts,*) 'TEMFCT= ',TEMFCT
c
cgzh debug output
c      if(js.eq.14.and.ks.eq.2) then
c        write(iouts,*)
c        write(iouts,*) 'second loop, cell=', js,is,ks
c        write(iouts,*) 'temfct,dxxcur*(temfct-1)',temfct,
c     *  DxxCUR(JS,IS,KS)*(temfct-1.0)
c        write(iouts,*)
c      end if
c      if(js.eq.22.and.ks.eq.2) then
c        write(iouts,*)
c        write(iouts,*) 'second loop, cell=', js,is,ks
c        write(iouts,*) 'temfct,dxxcur*(temfct-1)',temfct,
c     *  DxxCUR(JS,IS,KS)*(temfct-1.0)
c        write(iouts,*)
c      end if
C
C          AT FORWARD CELL, ADD THE DXX FLUX TO CNCNC
C          SINCE DXXFOR HAS SAME SIGN AS DXXCUR, JUST NEED TO CHANGE
C            SIGN WHEN APPLIED TO CNCNC
C
cgzh debug output
c          write(iouts,*) 'face flux adjustment at js,is,ks',JS,IS,KS 
c          write(iouts,*) 'dxx cncnc BEFORE adjustment at forward cell' 
c          write(iouts,*) 'DXXFOR= ',DXXFOR(JS,IS,KS)
c
cgzh debug  new way: just adjust the face flux
           DXXFOR(JS,IS,KS)=DXXFOR(JS,IS,KS)*TEMFCT
c           CNCNC(JS+1,IS,KS)=CNCNC(JS+1,IS,KS)-
c     *                       DXXFOR(JS,IS,KS)*(TEMFCT-1.D0)
cgzh debug output
c          write(iouts,*) 'dxx cncnc AFTER adjustment at forward cell' 
c          write(iouts,*) 'DXXFOR= ',DXXFOR(JS,IS,KS)
c
cgzh debug output
c      if(js.eq.14.and.ks.eq.2) then
c        write(iouts,*)
c        write(iouts,*) 'second loop, cell=', js,is,ks
c        write(iouts,*) 'temfct,dxxfor*(1-temfct)',temfct,
c     *  Dxxfor(JS,IS,KS)*(1.0-temfct)
c        write(iouts,*)
c      end if
c      if(js.eq.22.and.ks.eq.2) then
c        write(iouts,*)
c        write(iouts,*) 'second loop, cell=', js+1,is,ks
c        write(iouts,*) 'temfct,dxxfor*(1-temfct)',temfct,
c     *  Dxxfor(JS,IS,KS)*(1.0-temfct)
c        write(iouts,*)
c      end if
C
          END IF
         END IF
       END IF
      END IF
C
C  *************************
C  START DYY TERM ADJUSTMENT
C  *************************
C
C  SEE COMMENTS IN DXX SECTION
C
      IF(IS.LT.NSROW) THEN
       IF(ADJFCT(JS,IS+1,KS).GT.-1.0) THEN
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS,IS+1,KS).GT.0)) THEN
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS,IS+1,KS).LT.1.0) THEN
           TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS,IS+1,KS))
C          AT CURRENT CELL, SUBTRACT THE DYY FLUX FROM CNCNC
cgzh debug output
c          write(iouts,*) 'face flux adjustment at js,is,ks',JS,IS,KS 
c          write(iouts,*) 'dyy cncnc BEFORE adjustment at current cell' 
c          write(iouts,*) 'cncnc= ',CNCNC(JS,IS,KS)
c
cgzh debug  new way: just adjust the face flux
           DYYCUR(JS,IS,KS)=DYYCUR(JS,IS,KS)*TEMFCT
c           CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+
c     *                       DYYCUR(JS,IS,KS)*(TEMFCT-1.D0)
cgzh debug output
c          write(iouts,*) 'dyy cncnc AFTER adjustment at current cell' 
c          write(iouts,*) 'cncnc= ',CNCNC(JS,IS,KS)
c
C          AT FORWARD CELL, ADD THE DYY FLUX TO CNCNC
cgzh debug output
c          write(iouts,*) 'face flux adjustment at js,is,ks',JS,IS+1,KS 
c          write(iouts,*) 'dyy cncnc BEFORE adjustment at forward cell' 
c          write(iouts,*) 'cncnc= ',CNCNC(JS,IS+1,KS)
c
c bug fix           CNCNC(JS,IS+1,KS)=CNCNC(JS,IS+1,KS)-
c     *                       DYYFOR(JS,IS,KS)*(1.0-TEMFCT)
cgzh debug  new way: just adjust the face flux
           DYYFOR(JS,IS,KS)=DYYFOR(JS,IS,KS)*TEMFCT
c           CNCNC(JS,IS+1,KS)=CNCNC(JS,IS+1,KS)-
c     *                       DYYFOR(JS,IS,KS)*(TEMFCT-1.D0)
cgzh debug output
c          write(iouts,*) 'dyy cncnc AFTER adjustment at forward cell' 
c          write(iouts,*) 'cncnc= ',CNCNC(JS,IS+1,KS)
c
          END IF
         END IF
       END IF
      END IF
C
C  *************************
C  START DZZ TERM ADJUSTMENT
C  *************************
C
C  SEE COMMENTS IN DXX SECTION
C
c
      IF(KS.LT.NSLAY) THEN
       IF(ADJFCT(JS,IS,KS+1).GT.-1.0) THEN
         IF(PTWTON.EQ.0.OR.
     *     (PTWTON.EQ.1.AND.NPCELL(JS,IS,KS+1).GT.0)) THEN
          IF(ADJFCT(JS,IS,KS).LT.1.0.OR.ADJFCT(JS,IS,KS+1).LT.1.0) THEN
           TEMFCT=MIN(ADJFCT(JS,IS,KS),ADJFCT(JS,IS,KS+1))
C          AT CURRENT CELL, SUBTRACT THE DZZ FLUX FROM CNCNC
cgzh debug output
c          write(iouts,*) 'face flux adjustment at js,is,ks',JS,IS,KS 
c          write(iouts,*) 'DZZCUR BEFORE adjustment at current cell' 
c          write(iouts,*) 'DZZCUR= ',DZZCUR(JS,IS,KS)
c
cgzh debug  new way: just adjust the face flux
           DZZCUR(JS,IS,KS)=DZZCUR(JS,IS,KS)*TEMFCT
c           CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+
c     *                       DZZCUR(JS,IS,KS)*(TEMFCT-1.D0)
cgzh debug output
c          write(iouts,*) 'DZZCUR AFTER adjustment at current cell' 
c          write(iouts,*) 'DZZCUR= ',DZZCUR(JS,IS,KS)
c
cgzh debug output
c      if(js.eq.24.and.ks.eq.1) then
c        write(iouts,*)
c        write(iouts,*) 'second loop, cell=', js,is,ks
c        write(iouts,*) 'temfct,dzzcur*(temfct-1)',temfct,
c     *  DzzCUR(JS,IS,KS)*(temfct-1.0)
c        write(iouts,*)
c      end if
c      if(js.eq.22.and.ks.eq.2) then
c        write(iouts,*)
c        write(iouts,*) 'second loop, cell=', js,is,ks
c        write(iouts,*) 'temfct,dzzcur*(temfct-1)',temfct,
c     *  DzzCUR(JS,IS,KS)*(temfct-1.0)
c        write(iouts,*)
c      end if
C          AT FORWARD CELL, ADD THE DZZ FLUX TO CNCNC
cgzh debug output
c         write(iouts,*) 'face flux adjustment at js,is,ks+1',JS,IS,KS 
c         write(iouts,*) 'DZZFOR BEFORE adjustment at forward cell' 
c         write(iouts,*) 'DZZFOR= ',DZZFOR(JS,IS,KS)
c
c bug fix           CNCNC(JS,IS,KS+1)=CNCNC(JS,IS,KS+1)-
c     *                       DZZFOR(JS,IS,KS)*(1.0-TEMFCT)
cgzh debug  new way: just adjust the face flux
           DZZFOR(JS,IS,KS)=DZZFOR(JS,IS,KS)*TEMFCT
c           CNCNC(JS,IS,KS+1)=CNCNC(JS,IS,KS+1)-
c     *                       DZZFOR(JS,IS,KS)*(TEMFCT-1.D0)
cgzh debug output
c          write(iouts,*) 'DZZFOR AFTER adjustment at forward cell' 
c          write(iouts,*) 'DZZFOR= ',DZZFOR(JS,IS,KS)
cgzh debug output
c      if(js.eq.22.and.ks.eq.1) then
c        write(iouts,*) 
c        write(iouts,*) 'second loop, cell=', js,is,ks+1
c        write(iouts,*) 'temfct,dzzfor*(1-temfct)',temfct,
c     *  Dzzfor(JS,IS,KS)*(1.0-temfct)
c        write(iouts,*) 
c      end if
c      if(js.eq.22.and.ks.eq.2) then
c        write(iouts,*) 
c        write(iouts,*) 'second loop, cell=', js,is,ks+1
c        write(iouts,*) 'temfct,dzzfor*(1-temfct)',temfct,
c     *  Dzzfor(JS,IS,KS)*(1.0-temfct)
c        write(iouts,*) 
c      end if
c
          END IF
         END IF
       END IF
      END IF
cgzh debug new way: update cncnc at end
C     FORWARD FACE FLUXES
cgzh debug
c        write(72,*) 'A cncnc,dxxcur',CNCNC(JS,IS,KS),DXXCUR(JS,IS,KS)
      CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+DXXCUR(JS,IS,KS)
cgzh debug
c        write(72,*) 'B cncnc,dxxcur',CNCNC(JS,IS,KS),DXXCUR(JS,IS,KS)
      CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+DYYCUR(JS,IS,KS)
      CNCNC(JS,IS,KS)=CNCNC(JS,IS,KS)+DZZCUR(JS,IS,KS)
cgzh debug
cgzh debug output
c      if(DXXCUR(js,is,ks).ne.0.0) then
c        write(iouts,*) 'cncnc,dxxcur at js is ks'
c        write(iouts,*) CNCNC(JS,IS,KS),DXXCUR(JS,IS,KS),js,is,ks
c      end if
c      if(DyyCUR(js,is,ks).ne.0.0) then
c        write(iouts,*) 'cncnc,dyycur at js is ks'
c        write(iouts,*) CNCNC(JS,IS,KS),DyyCUR(JS,IS,KS),js,is,ks
c      end if
c      if(DzzCUR(js,is,ks).ne.0.0) then
c        write(iouts,*) 'cncnc,dzzcur at js is ks'
c        write(iouts,*) CNCNC(JS,IS,KS),DzzCUR(JS,IS,KS),js,is,ks
c      end if
c        write(72,*) 'C cncnc,dxxcur',CNCNC(JS,IS,KS),DXXCUR(JS,IS,KS)
C     BACKWARD FACE FLUXES
      IF(JS.LT.NSCOL) CNCNC(JS+1,IS,KS)=
     &  CNCNC(JS+1,IS,KS)-DXXFOR(JS,IS,KS)
      IF(IS.LT.NSROW) CNCNC(JS,IS+1,KS)=
     &  CNCNC(JS,IS+1,KS)-DYYFOR(JS,IS,KS) 
      IF(KS.LT.NSLAY) CNCNC(JS,IS,KS+1)=
     &  CNCNC(JS,IS,KS+1)-DZZFOR(JS,IS,KS) 
C  END LOOP FOR CNCNC ADJUSTMENT
   80 CONTINUE
C     ****************************************************************
C  RELEASE MEMORY
      DEALLOCATE(ADJFCT,DXXCUR,DXXFOR,DYYCUR,DYYFOR,DZZCUR,DZZFOR)
C
      RETURN
C     ****************************************************************
C
      END
C
C********************************************************
C
      SUBROUTINE CR5DSP(C,DC,
     X     TSXY,TSXZ,TSYX,TSYZ,TSZX,TSZY,
     X     RS,
     X     RHS,CIN,IBC,MRNO,
     X     NX,NY,NZ,NXY,NXYZ,NRN,
c chfb
     X     LOCHFB,ONCHFB,
     X     THCK,xd_mask)
C.....Assembles the matrix coefficients and right hand side vector
C.....     for the solute equation
      INTEGER CIN,IBC,MRNO,NX,NY,NZ,NXYZ
      INTEGER CELLNO
      INTEGER C10,C11,C12,C20,C21,C22,C30,C31,C32
      INTEGER C40,C41,C42,C50,C51,C52,C60,C61,C62
      INTEGER C70,C71,C72,C80,C81,C82,C90,C91,C92
      DIMENSION C(*),DC(*),
     X     TSXY(*),TSXZ(*),TSYX(*),TSYZ(*),
     X     TSZX(*),TSZY(*),
     X     RS(*),
     X     RHS(0:*),CIN(6,*),IBC(*),MRNO(*),
c chfb
     X     LOCHFB(NX,NY,NZ,2)
      DIMENSION THCK(NX,NY,NZ)
      LOGICAL xd_mask
      dimension xd_mask(0:nx+1,0:ny+1,0:nz+1)
C.....Function for cell number
      CELLNO(I,J,K)=(K-1)*NXY+(J-1)*NX+I
C.....Function for updated concentration
      CNP(MM)=C(MM)+DC(MM)
C.....Function for cross derivative dispersive flux term for xy and yx
      crdf1(xtc,cpp,cpm,c0p,c0m)=xtc*(cpp-cpm+c0p-c0m)
C.....Function for cross derivative dispersive flux term for xz and xy
      crdf2(xtc,cpp,cp0,cpm,c0p,c00,c0m,dxpp,dxpm,dx0p,dx0m) = 
     x     xtc*((cpp-cp0)/dxpp + (cp0-cpm)/dxpm + (c0p-c00)/dx0p
     x     + (c00-c0m)/dx0m)
C.....Function for cross derivative dispersive flux term for zx and zy
      crdf3(xtc,cpp,cmp,cp0,cm0,dx0,dxp) = 
     x     xtc*(dx0*(cpp-cmp) + dxp*(cp0-cm0))
C...
      DO 340 M=1,NXYZ
         MA=MRNO(M)
C.....Decode M into I,J,K
         CALL MTOIJK(M,I,J,K,NX,NY)
         js = i
         is = j
         ks = k
         MIMJK=CIN(3,M)
         mimjpk = 0
         IF(xd_mask(i-1,j+1,k)) MIMJPK=CELLNO(I-1,J+1,K)
         mimjmk = 0
         IF(xd_mask(i-1,j-1,k)) MIMJMK=CELLNO(I-1,J-1,K)
         mimjkp = 0
         IF(xd_mask(i-1,j,k+1)) MIMJKP=CELLNO(I-1,J,K+1)
         mimjkm = 0
         IF(xd_mask(i-1,j,k-1)) MIMJKM=CELLNO(I-1,J,K-1)
         MIPJK=CIN(4,M)
         mipjpk = 0
         IF(xd_mask(i+1,j+1,k)) MIPJPK=CELLNO(I+1,J+1,K)
         mipjmk = 0
         IF(xd_mask(i+1,j-1,k)) MIPJMK=CELLNO(I+1,J-1,K)
         mipjkp = 0
         IF(xd_mask(i+1,j,k+1)) MIPJKP=CELLNO(I+1,J,K+1)
         mipjkm = 0
         IF(xd_mask(i+1,j,k-1)) MIPJKM=CELLNO(I+1,J,K-1)
         MIJMK=CIN(2,M)
         MIJPK=CIN(5,M)
         MIJKM=CIN(1,M)
         mijmkm = 0
         IF(xd_mask(i,j-1,k-1)) MIJMKM=CELLNO(I,J-1,K-1)
         mijpkm = 0
         IF(xd_mask(i,j+1,k-1)) MIJPKM=CELLNO(I,J+1,K-1)
         MIJKP=CIN(6,M)
         mijmkp =0
         IF(xd_mask(i,j-1,k+1)) MIJMKP=CELLNO(I,J-1,K+1)
         mijpkp = 0
         IF(xd_mask(i,j+1,k+1)) MIJPKP=CELLNO(I,J+1,K+1)
cgzh debug
         mijk = 0
         IF(xd_mask(i,j,k)) mijk=CELLNO(I,J,K)
C.....Zero coefficients are used to suppress equation terms that are
C.....     not present due to geometry of boundaries and/or equations
C.....     that are not being solved
C
C.....Calculate explicit cross-derivative dispersive flux terms
         UCROSC=0.
c ... If any of the 6 nodes needed for a gradient are excluded, 
c ...      that flux term is zero.
c
c  create local cell scheme for C gradient calc's
c  this will enable cells across HFB's to be excluded
c  for current layer, relative to C2 which is at I,J,K:
c
c    C10  C11  C12
c    C20  C2   C22
c    C30  C31  C32
c
c for layer above (k-1)
c
c    C40  C41  C42
c    C50  C51  C52
c    C60  C61  C62
c
c for layer below (k+1)
c
c    C70  C71  C72
c    C80  C81  C82
c    C90  C91  C92
c
C.....X-direction
         FTXYDP=0.
         FTXYDM=0.
         IF(MIJMK.GT.0 .AND. MIJPK.GT.0) THEN
c BEGIN for xy, cell m (i,j,k)
            IF(MIPJK.GT.0 .and. mipjpk.gt.0 .and. mipjmk.gt.0) THEN
              C11=MIJMK
              C12=MIPJMK
		    C31=MIJPK
	        C32=MIPJPK
cgzh orig...looks like a bug	        CMULT=2.0
	        CMULT=1.0
c chfb 
              IF(ONCHFB.EQ.1) THEN 
               IF(LOCHFB(I,J-1,K,2).EQ.1.OR.
     X           LOCHFB(I+1,J-1,K,2).EQ.1) THEN
	             C11=CELLNO(I,J,K)
	             C12=MIPJK
	             CMULT=2.0
	         END IF
               IF(LOCHFB(I,J,K,2).EQ.1.OR.
     X           LOCHFB(I+1,J,K,2).EQ.1) THEN
	             C31=CELLNO(I,J,K)
	             C32=MIPJK
	             CMULT=2.0
	         END IF
c chfb end
	        END IF
cgzh
              if(mijk.gt.0)
     x        ftxydp=crdf1((TSXY(M)*CMULT),CNP(C32),CNP(C12),
     x           CNP(C31),CNP(C11))
c END for xy, cell m (i,j,k)
            END IF

c BEGIN for xy, cell m-1 (i-1,j,k)
            if(MIMJK.GT.0 .and. mimjpk.gt.0 .and. mimjmk.gt.0) THEN
			C10=MIMJMK
			C11=MIJMK
			C30=MIMJPK
			C31=MIJPK
			CMULT=1.0
c chfb 
              IF(ONCHFB.EQ.1) THEN 
                IF(LOCHFB(I-1,J-1,K,2).EQ.1.OR.
     X            LOCHFB(I,J-1,K,2).EQ.1) THEN
	              C10=MIMJK
	              C11=CELLNO(I,J,K)
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I-1,J,K,2).EQ.1.OR.
     X            LOCHFB(I,J,K,2).EQ.1) THEN
	              C30=MIMJK
	              C31=CELLNO(I,J,K)
	              CMULT=2.0
	          END IF
c chfb end
	        END IF
cgzh
              if(mijk.gt.0)
     x        ftxydm=crdf1((TSXY(M-1)*CMULT),cnp(C31),cnp(C11),
     x           cnp(C30),cnp(C10))
            END IF
         ENDIF
         FTXZDP=0.
         FTXZDM=0.
c*** some changes
         IF(MIJKM.GT.0 .AND. MIJKP.GT.0) THEN
            dthp = thck(js,is,ks+1)+thck(js,is,ks)
            dthm = thck(js,is,ks)+thck(js,is,ks-1)
            IF(MIPJK.GT.0 .and. mipjkp.gt.0 .and. mipjkm.gt.0) then
               dthpp = thck(js+1,is,ks+1)+thck(js+1,is,ks)
               dthpm = thck(js+1,is,ks)+thck(js+1,is,ks-1)
cgzh orig line...potential bug               FTXZDP=CRDF2(TSXY(M),CNP(MIPJKP),CNP(mipjk),CNP(MIPJKM),
cgzh
               if(mijk.gt.0)
     x          FTXZDP=CRDF2(TSXZ(M),CNP(MIPJKP),CNP(mipjk),CNP(MIPJKM),
     x		 	CNP(mijkp),cnp(m),cnp(mijkm),dthpp,dthpm,dthp,dthm)
            end if
            IF(MIMJK.GT.0 .and. mimjkp.gt.0 .and. mimjkm.gt.0) then
               dthmp = thck(js-1,is,ks+1)+thck(js-1,is,ks)
               dthmm = thck(js-1,is,ks)+thck(js-1,is,ks-1)
cgzh
              if(mijk.gt.0)
     x        ftxzdm=CRDF2(TSXZ(M-1),CNP(MIJKP),CNP(M),CNP(MIJKM),
     x         cnp(mimjkp),cnp(mimjk),cnp(mimjkm),dthp,dthm,dthmp,dthmm)
            end if
         ENDIF
C.....Y-direction
         FTYXDP=0.
         FTYXDM=0.
         IF(mimjk.gt.0 .and. mipjk.gt.0) THEN
c BEGIN for yx, cell m (i,j,k)
            IF(MIJPK.GT.0 .and. mipjpk.gt.0 .and. mimjpk.gt.0) THEN
			C20=MIMJK
			C30=MIMJPK
			C22=MIPJK
			C32=MIPJPK
			CMULT=1.0
c chfb
              IF(ONCHFB.EQ.1) THEN
                IF(LOCHFB(I-1,J,K,1).EQ.1.OR.
     X            LOCHFB(I-1,J+1,K,1).EQ.1) THEN
	              C20=CELLNO(I,J,K)
	              C30=MIJPK
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I,J,K,1).EQ.1.OR.
     X            LOCHFB(I,J+1,K,1).EQ.1) THEN
	              C22=CELLNO(I,J,K)
	              C32=MIJPK
	              CMULT=2.0
	          END IF
c chfb end
              END IF
cgzh
               if(mijk.gt.0)
     x         FTYXDP=CRDF1((TSYX(M)*CMULT),CNP(C32),CNP(C30),
     x           CNP(C22),CNP(C20))
c END for yx, cell m (i,j,k)
            END IF
c BEGIN for yx, cell m-1 (i,j-1,k)
cgzh bug fix for 1.8
            IF(MIJMK.GT.0 .AND. mipjmk.gt.0 .and. mimjmk.gt.0) THEN
			C10=MIMJMK
			C20=MIMJK
			C12=MIPJMK
			C22=MIPJK
			CMULT=1.0
c chfb
              IF(ONCHFB.EQ.1) THEN
                IF(LOCHFB(I-1,J-1,K,1).EQ.1.OR.
     X            LOCHFB(I-1,J,K,1).EQ.1) THEN
	              C10=MIJMK
	              C20=CELLNO(I,J,K)
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I,J-1,K,1).EQ.1.OR.
     X            LOCHFB(I,J,K,1).EQ.1) THEN
	              C12=MIJMK
	              C22=CELLNO(I,J,K)
	              CMULT=2.0
	          END IF
c chfb end
              END IF
cgzh
              if(mijk.gt.0)
     x        FTYXDM=CRDF1((TSYX(MIJMK)*CMULT),CNP(C22),CNP(C20),
     X           CNP(C12),CNP(C10))
c END for yx, cell m-1 (i,j-1,k)
            END IF
         ENDIF
         FTYZDP=0.
         FTYZDM=0.
c*** start changes again
         IF(MIJKM.GT.0 .AND. MIJKP.GT.0) THEN
            dthp = thck(js,is,ks+1)+thck(js,is,ks)
            dthm = thck(js,is,ks)+thck(js,is,ks-1)
            IF(MIJPK.GT.0 .and. mijpkp.gt.0 .and. mijpkm.gt.0) then
               dthpp = thck(js,is+1,ks+1)+thck(js,is+1,ks)
               dthpm = thck(js,is+1,ks)+thck(js,is+1,ks-1)
cgzh
               if(mijk.gt.0)
     x         ftyzdp=crdf2(tsyz(m),cnp(mijpkp),cnp(mijpk),
     x              cnp(mijpkm),cnp(mijkp),cnp(m),cnp(mijkm),
     x              dthpp,dthpm,dthp,dthm)
            end if
            IF(MIJMK.GT.0 .and. mijmkp.gt.0 .and. mijmkm.gt.0) then
               dthmp = thck(js,is-1,ks+1)+thck(js,is-1,ks)
               dthmm = thck(js,is-1,ks)+thck(js,is-1,ks-1)
cgzh
              if(mijk.gt.0)
     x        ftyzdm=crdf2(tsyz(mijmk),cnp(mijkp),cnp(m),
     x              cnp(mijkm),cnp(mijmkp),cnp(mijmk),cnp(mijmkm),
     x              dthp,dthm,dthmp,dthmm)
            end if
         ENDIF
C.....Z-direction
         FTZXDP=0.
         FTZXDM=0.
         IF(MIMJK.GT.0 .AND. MIPJK.GT.0) THEN
c BEGIN for zx, cell m (i,j,k)
            IF(MIJKP.GT.0 .and. mipjkp.gt.0 .and. mimjkp.gt.0) THEN
			C82=MIPJKP
			C22=MIPJK
			C80=MIMJKP
			C20=MIMJK
			CMULT=1.0
c chfb
              IF(ONCHFB.EQ.1) THEN
                IF(LOCHFB(I-1,J,K,1).EQ.1.OR.
     X            LOCHFB(I-1,J,K+1,1).EQ.1) THEN
	              C20=CELLNO(I,J,K)
	              C80=MIJKP
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I,J,K,1).EQ.1.OR.
     X            LOCHFB(I,J,K+1,1).EQ.1) THEN
	              C22=CELLNO(I,J,K)
	              C82=MIJKP
	              CMULT=2.0
	          END IF
c chfb end
              END IF
cgzh
               if(mijk.gt.0)
     x         ftzxdp=crdf3((tszx(m)*CMULT),cnp(C82),cnp(C80),
     x           cnp(C22),cnp(C20),thck(js,is,ks),thck(js,is,ks+1))
c END for zx, cell m (i,j,k)
            END IF
c BEGIN for zx, cell m-1 (i,j,k-1)
cgzh bug fix for 1.8
            IF(mijkm.gt.0 .and. mipjkm.gt.0 .and. mimjkm.gt.0) THEN
			C22=MIPJK
			C52=MIPJKM
			C50=MIMJKM
			C20=MIMJK
			CMULT=1.0
c chfb
              IF(ONCHFB.EQ.1) THEN
                IF(LOCHFB(I-1,J,K,1).EQ.1.OR.
     X            LOCHFB(I-1,J,K-1,1).EQ.1) THEN
	              C20=CELLNO(I,J,K)
	              C50=MIJKM
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I,J,K,1).EQ.1.OR.
     X            LOCHFB(I,J,K-1,1).EQ.1) THEN
	              C22=CELLNO(I,J,K)
	              C52=MIJKM
	              CMULT=2.0
	          END IF
c chfb end
              END IF
cgzh
              if(mijk.gt.0)
     x        ftzxdm=crdf3((tszx(mijkm)*CMULT),cnp(C22),
     x           cnp(C20),cnp(C52),cnp(C50),
     x           thck(js,is,ks-1),thck(js,is,ks))
c END for zx, cell m-1 (i,j,k)
            END IF
        ENDIF
         FTZYDP=0.
         FTZYDM=0.
         IF(MIJMK.GT.0 .AND. MIJPK.GT.0) THEN
c BEGIN for zy, cell m (i,j,k)
            IF(MIJKP.GT.0 .and. mijpkp.gt.0 .and. mijmkp.gt.0) THEN
			C91=MIJPKP
			C31=MIJPK
			C71=MIJMKP
			C11=MIJMK
			CMULT=1.0
c chfb
              IF(ONCHFB.EQ.1) THEN
                IF(LOCHFB(I-1,J,K,2).EQ.1.OR.
     X            LOCHFB(I-1,J,K+1,2).EQ.1) THEN
	              C11=CELLNO(I,J,K)
	              C71=MIJKP
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I,J,K,2).EQ.1.OR.
     X            LOCHFB(I,J,K+1,2).EQ.1) THEN
	              C31=CELLNO(I,J,K)
	              C91=MIJKP
	              CMULT=2.0
	          END IF
c chfb end
              END IF
cgzh
               if(mijk.gt.0)
     x         ftzydp=crdf3((tszy(m)*CMULT),cnp(C91),cnp(C71),
     x           cnp(C31),cnp(C11),thck(js,is,ks),thck(js,is,ks+1))
c END for zy, cell m (i,j,k)
            END IF
c BEGIN for zy, cell m-1 (i,j,k-1)
corig            IF(mijpkm.gt.0 .and. mijmkm.gt.0) THEN
            IF(mijkm.gt.0 .and. mijpkm.gt.0 .and. mijmkm.gt.0) THEN
 	        C31=MIJPK
			C61=MIJPKM
			C11=MIJMK
			C41=MIJMKM
			CMULT=1.0
c chfb
              IF(ONCHFB.EQ.1) THEN
                IF(LOCHFB(I-1,J,K,2).EQ.1.OR.
     X            LOCHFB(I-1,J,K-1,2).EQ.1) THEN
	              C11=CELLNO(I,J,K)
	              C41=MIJKM
	              CMULT=2.0
	          END IF
                IF(LOCHFB(I,J,K,2).EQ.1.OR.
     X            LOCHFB(I,J,K-1,2).EQ.1) THEN
	              C31=CELLNO(I,J,K)
	              C61=MIJKM
	              CMULT=2.0
	          END IF
c chfb end
              END IF
cgzh
              if(mijk.gt.0)
     x        ftzydm=crdf3((tszy(mijkm)*CMULT),cnp(C31),
     x           cnp(C11),cnp(C61),cnp(C41),
     x           thck(js,is,ks-1),thck(js,is,ks))
c END for zy, cell m-1 (i,j,k)
            END IF
         ENDIF
         UCROSC=FTXYDP+FTXZDP+FTYXDP+FTYZDP+FTZXDP+FTZYDP-
     x        FTXYDM-FTXZDM-FTYXDM-FTYZDM-FTZXDM-FTZYDM
C.....Build right-hand-side
         RHS(MA)=RS(M)+UCROSC
  340 CONTINUE
      RETURN
      END
C
C
      SUBROUTINE RHSN(C,TSX,TSY,TSZ,RS,
     X     NX,NY,NZ,NXY,NXYZ,
c dzz
     X     NPCELL)
CMOCWT
      INCLUDE 'ptwt.inc'
C.....Calculates right hand side terms at time level n,
      INTEGER NX,NY,NZ,NXYZ
      DIMENSION C(NXYZ),TSX(NXYZ),TSY(NXYZ),TSZ(NXYZ),RS(NXYZ)
c dzz
      DIMENSION NPCELL(NX,NY,NZ)
      INTEGER CELLNO
C...
C.....Cell number function
      CELLNO(I,J,K)=(K-1)*NXY+(J-1)*NX+I
C...
C.....Calculate right hand side dispersive flux terms
C.....      (not cross-dispersive flux terms)
C.....Inactive cells are excluded by zero flow rate and transmissivity
c ... Assembly is by cell face with facial area adjustments in horizontal
      DO 252 K=1,NZ
         DO 251 J=1,NY
            DO 250 I=1,NX
               M=CELLNO(I,J,K)
               TS1=TSX(M)
               TS2=TSY(M)
               TS3=TSZ(M)
cgzh account for npcell=0
cgzh debug only for mocwt (npcell checks)
               IF(PTWTON.EQ.1.AND.NPCELL(I,J,K).LT.1) THEN
cgzh debug output
c      if(i.eq.20.and.j.eq.7.and.k.eq.1) then
c	 write(21,*) 'ts zeroed for 20,7'
c	end if
			   TS1=0.0 
			   TS2=0.0 
			   TS3=0.0 
               END IF
C
               IF(I.LT.NX) THEN
C.....X-direction
corig                  RS(M+1)=RS(M+1)-TSX(M)*(C(M+1)-C(M))
corig                  RS(M)=RS(M)+TSX(M)*(C(M+1)-C(M))
cgzh debug only for mocwt (npcell checks)
                  IF(PTWTON.EQ.1.AND.NPCELL(I+1,J,K).LT.1) TS1=0.0
cgzh debug output
c      if(i.eq.19.and.j.eq.7.and.k.eq.1.and.(NPCELL(I+1,J,K).LT.1)) then
c	 write(21,*) 'tsjp zeroed for 19,7'
c	end if
				RS(M+1)=RS(M+1)-TS1*(C(M+1)-C(M))
                  RS(M)=RS(M)+TS1*(C(M+1)-C(M))
               ENDIF
               IF(J.LT.NY) THEN
C.....Y-direction
                  MIJPK=CELLNO(I,J+1,K)
corig                  RS(MIJPK)=RS(MIJPK)-TSY(M)*(C(MIJPK)-C(M))
corig                  RS(M)=RS(M)+TSY(M)*(C(MIJPK)-C(M))
cgzh debug only for mocwt (npcell checks)
                  IF(PTWTON.EQ.1.AND.NPCELL(I,J+1,K).LT.1) TS2=0.0
cgzh debug output
c      if(i.eq.20.and.j.eq.6.and.k.eq.1.and.(NPCELL(I,J+1,K).LT.1)) then
c	 write(21,*) 'tsip zeroed for 20,6'
c	end if
                  RS(MIJPK)=RS(MIJPK)-TS2*(C(MIJPK)-C(M))
                  RS(M)=RS(M)+TS2*(C(MIJPK)-C(M))
               ENDIF
               IF(K.LT.NZ) THEN
C.....Z-direction
                  MIJKP=CELLNO(I,J,K+1)
corig                  RS(MIJKP)=RS(MIJKP)-TSZ(M)*(C(MIJKP)-C(M))
corig                  RS(M)=RS(M)+TSZ(M)*(C(MIJKP)-C(M))
cgzh debug only for mocwt (npcell checks)
                  IF(PTWTON.EQ.1.AND.NPCELL(I,J,K+1).LT.1) TS3=0.0
                  RS(MIJKP)=RS(MIJKP)-TS3*(C(MIJKP)-C(M))
                  RS(M)=RS(M)+TS3*(C(MIJKP)-C(M))
               ENDIF
c	      end if
  250       CONTINUE
  251    CONTINUE
  252 CONTINUE
      RETURN
      END
C
C
CMOCWT ADDED NPCELL<1 TO EXCLUDED CELLS FOR MOCWT
      SUBROUTINE ASEMBL(
     X     TSX,TSY,TSZ,
     X     VA,VAD,RHS,CIN,IBOUND,MRNO,FDTMTH,
     X     NX,NY,NZ,NXY,NXYZ,NRN,
     X     NCOL,NROW,NLAY,
c dzz
     X     RF,POR,THCK,
     X     NPCELL)
C.....Assembles the matrix coefficients 
C.....     for the solute equation
CMOCWT
      INCLUDE 'ptwt.inc'
      INTEGER CIN,MRNO,NX,NY,NZ,NXYZ
      DIMENSION 
     X     TSX(NXYZ),TSY(NXYZ),TSZ(NXYZ),
     X     VA(6,*),
     X     VAD(NXYZ),RHS(0:NXYZ),CIN(6,*),MRNO(*)
      DIMENSION IBOUND(NCOL,NROW,NLAY),NPCELL(NX,NY,NZ)
c dzz
      DIMENSION RF(NX,NY,NZ),POR(NX,NY,NZ),THCK(NX,NY,NZ)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      PARAMETER(ICXM=3,ICXP=4,ICYM=2,ICYP=5,ICZM=1,ICZP=6)
c initialize RHS(0)
      RHS(0)=0.0
C...
C.....Compute and assemble coefficients in difference equations
C.....     cell-by-cell
      DO 340 M=1,NXYZ
         MA=MRNO(M)
C.....Asemble VA for red nodes only and VAD and RHS for all nodes
         IF(MA.LE.NRN) THEN
         DO 10 IC=1,6
            VA(IC,MA)=0.
   10    CONTINUE
         ENDIF
C.....Decode M into I,J,K
         CALL MTOIJK(M,I,J,K,NX,NY)
C.....Solve trivial equation for excluded cells, iterative solver
         II=I+ISCOL1-1
         JJ=J+ISROW1-1
         KK=K+ISLAY1-1
         IF(IBOUND(II,JJ,KK).EQ.0) THEN
            VAD(MA)=1.
            RHS(MA)=0.
            GO TO 340
         ENDIF
CMOCWT  EXCLUDE NPCELL<1 CELLS FROM RHS
         IF(PTWTON.EQ.1.AND.NPCELL(I,J,K).LT.1) THEN
            VAD(MA)=1.
            RHS(MA)=0.
cgzh debug output
c      write(21,*) 'j,i,k skipped in ASEMBL',i,j,k
            GO TO 340
         ENDIF
C.....Conductances treated explicitly
CMOCWT  NPCELL<1 for below branches is done by redifining CIN in
C       CELCONEC, which should be called just before this routine
         TSXM=0.
         MIMJK=CIN(3,M)
         IF(MIMJK.GT.0) THEN
            TSXM=TSX(MIMJK)
         ENDIF
         TSXP=0.
         MIPJK=CIN(4,M)
         IF(MIPJK.GT.0) THEN
            TSXP=TSX(M)
         ENDIF
         TSYM=0.
         MIJMK=CIN(2,M)
         IF(MIJMK.GT.0) TSYM=TSY(MIJMK)
         TSYP=0.
         MIJPK=CIN(5,M)
         IF(MIJPK.GT.0) TSYP=TSY(M)
         TSZM=0.
         MIJKM=CIN(1,M)
         IF(MIJKM.GT.0) TSZM = TSZ(MIJKM)
         TSZP=0.
         MIJKP=CIN(6,M)
         IF(MIJKP.GT.0) TSZP=TSZ(M)
C.....Zero coefficients are used to suppress equation terms that are
C.....     not present due to geometry of boundaries and/or equations
C.....     that are not being solved
c ... Capacitance term appears in diagonal
         vad(ma) = rf(i,j,k)*por(i,j,k)*thck(i,j,k)
C.....Install conductance terms for red nodes only
         IF(MA.LE.NRN) THEN
C.....X-direction
         VA(ICXM,MA)=-FDTMTH*TSXM
         VA(ICXP,MA)=-FDTMTH*TSXP
C.....Y-direction
         VA(ICYM,MA)=-FDTMTH*TSYM
         VA(ICYP,MA)=-FDTMTH*TSYP
C.....Z-direction
         VA(ICZM,MA)=-FDTMTH*TSZM
         VA(ICZP,MA)=-FDTMTH*TSZP
         ENDIF
cgzh debug output
c      if(i.eq.20.and.j.eq.6.and.k.eq.1) write(21,*) 
c     *  'VA(ICYP,MA)',VA(ICYP,MA)
c      if(i.eq.19.and.j.eq.7.and.k.eq.1) write(21,*) 
c     *  'VA(ICXP,MA)',VA(ICXP,MA)
c      if(i.eq.21.and.j.eq.7.and.k.eq.1) write(21,*) 
c     *  'VA(ICXM,MA)',VA(ICXM,MA)
c      if(i.eq.20.and.j.eq.8.and.k.eq.1) write(21,*) 
c     *  'VA(ICYM,MA)',VA(ICYM,MA)
c      if(i.eq.20.and.j.eq.7.and.k.eq.2) write(21,*) 
c     *  'VA(ICZM,MA)',VA(ICZM,MA)
c      if(i.eq.20.and.j.eq.7.and.k.eq.0) write(21,*) 
c     *  'VA(ICZP,MA)',VA(ICZP,MA)
C.....Diagonal term
         VAD(MA)=VAD(MA)+FDTMTH*(TSXM+TSXP+TSYM+TSYP+TSZM+TSZP)
cgzh debug output
c      if(iimov.eq.80) then
c	write(78,*)
c	do 888 m=1,nodess
c	 do 888 p=1,6
c         MA=MRNO(M)
c         IF(MA.LE.NRN) THEN
c         CALL MTOIJK(M,I,J,K,NSCOL,NSROW)
C      PARAMETER(ICXM=3,ICXP=4,ICYM=2,ICYP=5,ICZM=1,ICZP=6)
c 	   if (va(p,ma).ne.0.0) then
c	     write(78,*) 'ma, va',ma,va(p,ma)
c 	     write(78,*) 'j,i,k',i,j,k
c           if(p.eq.1) write(78,*) 'k-minus, tszm=',tszm
c           if(p.eq.2) write(78,*) 'i-minus, tsym=',tsym
c           if(p.eq.3) write(78,*) 'j-minus, tsxm=',tsxm
c           if(p.eq.4) write(78,*) 'j-plus, tsxm=',tsxp
c           if(p.eq.5) write(78,*) 'i-plus, tsym=',tsyp
c           if(p.eq.6) write(78,*) 'k-plus, tszm=',tszp
c           write(78,*) '--------------------------------'
c         end if
c	   end if
c888   continue
c      end if
cgzh debug output end
  340 CONTINUE
c           write(78,*) '*********END OF ASEMBL CALL**********'
      RETURN
      END
C
C     ****************************************************************
C
      SUBROUTINE LOADT(
     *     IBOUND,NCOL,NROW,NLAY,
     *     DISPCC,DISPCR,DISPCL,DISPRR,DISPRC,DISPRL,
     *     DISPLL,DISPLC,DISPLR,
     *     NSCOL,NSROW,NSLAY,
     X     THCK,POR,RF,TIMV)
C.....Loads the conductance factors into one-dimensional arrays from
C.....     the 3-dimensional arrays
C.....     in GWT with extra factor
      DIMENSION
     *  IBOUND(NCOL,NROW,NLAY),
     *  THCK(NSCOL,NSROW,NSLAY),
     *  POR(NSCOL,NSROW,NSLAY),RF(NSCOL,NSROW,NSLAY)
      DIMENSION
     *  DISPCC(NSCOL,NSROW,NSLAY),DISPCR(NSCOL,NSROW,NSLAY),
     *  DISPCL(NSCOL,NSROW,NSLAY),DISPRR(NSCOL,NSROW,NSLAY),
     *  DISPRC(NSCOL,NSROW,NSLAY),DISPRL(NSCOL,NSROW,NSLAY),
     *  DISPLL(NSCOL,NSROW,NSLAY),DISPLC(NSCOL,NSROW,NSLAY),
     *  DISPLR(NSCOL,NSROW,NSLAY)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C...
      DO 10 IZ=1,NSLAY
         KK=IZ+ISLAY1-1
         DO 20 IY=1,NSROW
            II=IY+ISROW1-1
            DO 30 IX=1,NSCOL
               JJ=IX+ISCOL1-1
               IF(IBOUND(JJ,II,KK).NE.0) THEN
c orig
c               TORF2=TIMV/(POR(IX,IY,IZ)*RF(IX,IY,IZ))
c               TRFPRB=TORF2/THCK(IX,IY,IZ)
c               DISPCC(IX,IY,IZ)=DISPCC(IX,IY,IZ)*TRFPRB
c               DISPCR(IX,IY,IZ)=DISPCR(IX,IY,IZ)*TRFPRB
c               DISPCL(IX,IY,IZ)=DISPCL(IX,IY,IZ)*TRFPRB
c               DISPRR(IX,IY,IZ)=DISPRR(IX,IY,IZ)*TRFPRB
c               DISPRC(IX,IY,IZ)=DISPRC(IX,IY,IZ)*TRFPRB
c               DISPRL(IX,IY,IZ)=DISPRL(IX,IY,IZ)*TRFPRB
c               DISPLL(IX,IY,IZ)=DISPLL(IX,IY,IZ)*TORF2
c               DISPLC(IX,IY,IZ)=DISPLC(IX,IY,IZ)*TORF2
c               DISPLR(IX,IY,IZ)=DISPLR(IX,IY,IZ)*TORF2
                  dispcc(ix,iy,iz)=dispcc(ix,iy,iz)*timv
                  dispcr(ix,iy,iz)=dispcr(ix,iy,iz)*timv
                  dispcl(ix,iy,iz)=dispcl(ix,iy,iz)*timv
                  disprr(ix,iy,iz)=disprr(ix,iy,iz)*timv
                  disprc(ix,iy,iz)=disprc(ix,iy,iz)*timv
                  disprl(ix,iy,iz)=disprl(ix,iy,iz)*timv
                  displl(ix,iy,iz)=displl(ix,iy,iz)*timv
                  displc(ix,iy,iz)=displc(ix,iy,iz)*timv
                  displr(ix,iy,iz)=displr(ix,iy,iz)*timv
               ENDIF
   30       CONTINUE
   20    CONTINUE
   10 CONTINUE
      RETURN
      END
C
      SUBROUTINE MTOIJK(M,I,J,K,NX,NY)
C.....Returns the index (I,J,K) of the point with
C.....     natural index M.
C...
      NXY = NX*NY
      IMOD = MOD(M,NXY)
      K = (M-IMOD)/NXY + MIN(1,IMOD)
      KR = IMOD
      IF (KR .EQ. 0) KR = NXY
      IMOD = MOD(KR,NX)
      J = (KR-IMOD)/NX + MIN(1,IMOD)
      I = IMOD
      IF (I .EQ. 0) I = NX
      RETURN
      END
C
C
C
C     ELLDIS BUILD DISPERSION INTEGRAL 
C***********************************************************************
C
      SUBROUTINE ELLDIS(IBOUND,DISPCC,DISPCR,DISPCL,DISPRC,DISPRR,
     *               DISPRL,DISPLC,DISPLR,DISPLL,
     *      POR,THCK,RF,DIAGDS,TIMV,
     *      NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NODESS,
     *      DELCOL,DELROW)
C
C***********************************************************************
C   ELLDIS IS CALLED AT THE FIRST TRANSPORT TIME STEP TO USE DISPERSION
C   COEFFICIENTS CALCULATED BY DSP6 TO BUILD THE LHS DISPERSION MATRIX,
C   DIAGDS.  DIAGDS HAS THE SAME 27 DIAGONAL FORM AS THE STORAGE MATRIX,
C   AND LIKE DIAGST IS STORED AS A 27 X NODESS ARRAY.
C   AT THE FIRST SOLVE (IMOV=1) DIAGST AND DIAGDS ARE SUMMED TO PRODUCE
C   THE SYSTEM MATRIX FOR SOLUTION, AND CONVERTED TO THE SLAP COLUMN
C   FORMAT USED BY THE SOLVER.  THE CONVERTED FORM OF DIAGDS+DIAGST IS
C   THEN STORED IN DIAGDS FOR USE IN FUTURE TRANSPORT SOLVES DURING 
C   DURING SINGLE FLOW TIME STEP.
C   
      DIMENSION  IBOUND(NCOL,NROW,NLAY),DISPCC(NSCOL,NSROW,NSLAY),
     *    DISPCR(NSCOL,NSROW,NSLAY),
     *    DISPCL(NSCOL,NSROW,NSLAY),DISPRC(NSCOL,NSROW,NSLAY),
     *    DISPRR(NSCOL,NSROW,NSLAY),DISPRL(NSCOL,NSROW,NSLAY),
     *    DISPLC(NSCOL,NSROW,NSLAY),DISPLR(NSCOL,NSROW,NSLAY),
     *    DISPLL(NSCOL,NSROW,NSLAY)
      DIMENSION POR(NODESS),THCK(NSCOL,NSROW,NSLAY),
     *    RF(NSCOL,NSROW,NSLAY),
     *    DIAGDS(27,NODESS),DELCOL(NCOL),DELROW(NROW)
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
C   DISPCC CONTAINS POR*THCK*DCC/CDEL
C   DISPCR CONTAINS POR*THCK*DCR/(4*RDEL)
C   DISPCL CONTAINS POR*THCK*DCL/(4*B)
C   DISPRC CONTAINS POR*THCK*DRC/(4*CDEL)
C   DISPRR CONTAINS POR*THCK*DRR/RDEL
C   DISPRL CONTAINS POR*THCK*DRL/(4*B)
C   DISPLC CONTAINS POR*DLC/(4*CDEL)
C   DISPLR CONTAINS POR*DLR/(4*RDEL)
C   DISPLL CONTAINS POR*DLL/B
C
          NSRC=NSROW*NSCOL
        DO 131 KS=1,NSLAY
          K=ISLAY1+KS-1
        DO 131 IS=1,NSROW
          I=ISROW1+IS-1
        DO 131 JS=1,NSCOL
          J=ISCOL1+JS-1
          NODE=(KS-1)*NSRC+(IS-1)*NSCOL+JS
C
      IF (IBOUND(J,I,K).EQ.0) THEN
        DIAGDS(14,NODE)=0D0
        GOTO 131
      ENDIF
C
      RFPINV=1.D0/(RF(JS,IS,KS)*POR(NODE))
      TRFPR=TIMV*RFPINV
C
      JSP=JS+1
      JSM=JS-1
      ISP=IS+1
      ISM=IS-1
      KSP=KS+1
      KSM=KS-1
C
C  ADD DISPCX COEFFS TO LHS
C
      RTRFPR=DELROW(I)*TRFPR
C
C  FORWARD COEFFS
      IF (JS.NE.NSCOL) THEN
        IF (KS.NE.1) THEN
          DIAGDS(5,NODE)=DIAGDS(5,NODE)+RTRFPR*DISPCL(JS,IS,KS)
          DIAGDS(6,NODE)=DIAGDS(6,NODE)+RTRFPR*DISPCL(JS,IS,KS)
        ENDIF
        IF (IS.NE.1) THEN
          DIAGDS(11,NODE)=DIAGDS(11,NODE)+RTRFPR*DISPCR(JS,IS,KS)
          DIAGDS(12,NODE)=DIAGDS(12,NODE)+RTRFPR*DISPCR(JS,IS,KS)
        ENDIF
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+RTRFPR*DISPCC(JS,IS,KS)
        DIAGDS(15,NODE)=DIAGDS(15,NODE)-RTRFPR*DISPCC(JS,IS,KS)
        IF (IS.NE.NSROW) THEN
          DIAGDS(17,NODE)=DIAGDS(17,NODE)-RTRFPR*DISPCR(JS,IS,KS)
          DIAGDS(18,NODE)=DIAGDS(18,NODE)-RTRFPR*DISPCR(JS,IS,KS)
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(23,NODE)=DIAGDS(23,NODE)-RTRFPR*DISPCL(JS,IS,KS)
          DIAGDS(24,NODE)=DIAGDS(24,NODE)-RTRFPR*DISPCL(JS,IS,KS)
        ENDIF
      ENDIF
C  BACKWARD COEFFS
      IF (JS.NE.1) THEN
        IF (KS.NE.1) THEN
          DIAGDS(4,NODE)=DIAGDS(4,NODE)-RTRFPR*DISPCL(JSM,IS,KS)
          DIAGDS(5,NODE)=DIAGDS(5,NODE)-RTRFPR*DISPCL(JSM,IS,KS)
        ENDIF
        IF (IS.NE.1) THEN
          DIAGDS(10,NODE)=DIAGDS(10,NODE)-RTRFPR*DISPCR(JSM,IS,KS)
          DIAGDS(11,NODE)=DIAGDS(11,NODE)-RTRFPR*DISPCR(JSM,IS,KS)
        ENDIF
        DIAGDS(13,NODE)=DIAGDS(13,NODE)-RTRFPR*DISPCC(JSM,IS,KS)
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+RTRFPR*DISPCC(JSM,IS,KS)
        IF (IS.NE.NSROW) THEN
          DIAGDS(16,NODE)=DIAGDS(16,NODE)+RTRFPR*DISPCR(JSM,IS,KS)
          DIAGDS(17,NODE)=DIAGDS(17,NODE)+RTRFPR*DISPCR(JSM,IS,KS)
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(22,NODE)=DIAGDS(22,NODE)+RTRFPR*DISPCL(JSM,IS,KS)
          DIAGDS(23,NODE)=DIAGDS(23,NODE)+RTRFPR*DISPCL(JSM,IS,KS)
        ENDIF
      ENDIF
C
C  ADD DISPRX COEFFS TO LHS
C
      CTRFPR=DELCOL(J)*TRFPR
C
C  FORWARD COEFFS
      IF (IS.NE.NSROW) THEN
        IF (KS.NE.1) THEN
          DIAGDS(5,NODE)=DIAGDS(5,NODE)+CTRFPR*DISPRL(JS,IS,KS)
          DIAGDS(8,NODE)=DIAGDS(8,NODE)+CTRFPR*DISPRL(JS,IS,KS)
        ENDIF
        IF (JS.NE.1) THEN
          DIAGDS(13,NODE)=DIAGDS(13,NODE)+CTRFPR*DISPRC(JS,IS,KS)
          DIAGDS(16,NODE)=DIAGDS(16,NODE)+CTRFPR*DISPRC(JS,IS,KS)
        ENDIF
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CTRFPR*DISPRR(JS,IS,KS)
        DIAGDS(17,NODE)=DIAGDS(17,NODE)-CTRFPR*DISPRR(JS,IS,KS)
        IF (JS.NE.NSCOL) THEN
          DIAGDS(15,NODE)=DIAGDS(15,NODE)-CTRFPR*DISPRC(JS,IS,KS)
          DIAGDS(18,NODE)=DIAGDS(18,NODE)-CTRFPR*DISPRC(JS,IS,KS)
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(23,NODE)=DIAGDS(23,NODE)-CTRFPR*DISPRL(JS,IS,KS)
          DIAGDS(26,NODE)=DIAGDS(26,NODE)-CTRFPR*DISPRL(JS,IS,KS)
        ENDIF
      ENDIF
C  BACKWARD COEFFS
      IF (IS.NE.1) THEN
        IF (KS.NE.1) THEN
          DIAGDS(2,NODE)=DIAGDS(2,NODE)-CTRFPR*DISPRL(JS,ISM,KS)
          DIAGDS(5,NODE)=DIAGDS(5,NODE)-CTRFPR*DISPRL(JS,ISM,KS)
        ENDIF
        IF (JS.NE.1) THEN
          DIAGDS(10,NODE)=DIAGDS(10,NODE)-CTRFPR*DISPRC(JS,ISM,KS)
          DIAGDS(13,NODE)=DIAGDS(13,NODE)-CTRFPR*DISPRC(JS,ISM,KS)
        ENDIF
        DIAGDS(11,NODE)=DIAGDS(11,NODE)-CTRFPR*DISPRR(JS,ISM,KS)
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CTRFPR*DISPRR(JS,ISM,KS)
        IF (JS.NE.NSCOL) THEN
          DIAGDS(12,NODE)=DIAGDS(12,NODE)+CTRFPR*DISPRC(JS,ISM,KS)
          DIAGDS(15,NODE)=DIAGDS(15,NODE)+CTRFPR*DISPRC(JS,ISM,KS)
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(20,NODE)=DIAGDS(20,NODE)+CTRFPR*DISPRL(JS,ISM,KS)
          DIAGDS(23,NODE)=DIAGDS(23,NODE)+CTRFPR*DISPRL(JS,ISM,KS)
        ENDIF
      ENDIF
C
C  ADD DISPLX COEFFS TO LHS
C
      CRTRFPR=DELCOL(J)*DELROW(I)*TRFPR
C
C  FORWARD COEFFS
      IF (KS.NE.NSLAY) THEN
        IF (JS.NE.1) THEN
          DIAGDS(13,NODE)=DIAGDS(13,NODE)+CRTRFPR*DISPLC(JS,IS,KS)
          DIAGDS(22,NODE)=DIAGDS(22,NODE)+CRTRFPR*DISPLC(JS,IS,KS)
        ENDIF
        IF (IS.NE.1) THEN
          DIAGDS(11,NODE)=DIAGDS(11,NODE)+CRTRFPR*DISPLR(JS,IS,KS)
          DIAGDS(20,NODE)=DIAGDS(20,NODE)+CRTRFPR*DISPLR(JS,IS,KS)
        ENDIF
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLL(JS,IS,KS)
        DIAGDS(23,NODE)=DIAGDS(23,NODE)-CRTRFPR*DISPLL(JS,IS,KS)
        IF (JS.NE.NSCOL) THEN
          DIAGDS(15,NODE)=DIAGDS(15,NODE)-CRTRFPR*DISPLC(JS,IS,KS)
          DIAGDS(24,NODE)=DIAGDS(24,NODE)-CRTRFPR*DISPLC(JS,IS,KS)
        ENDIF
        IF (IS.NE.NSROW) THEN
          DIAGDS(17,NODE)=DIAGDS(17,NODE)-CRTRFPR*DISPLR(JS,IS,KS)
          DIAGDS(26,NODE)=DIAGDS(26,NODE)-CRTRFPR*DISPLR(JS,IS,KS)
        ENDIF
      ENDIF
C  BACKWARD COEFFS
      IF (KS.NE.1) THEN
        IF (JS.NE.1) THEN
          DIAGDS(4,NODE)=DIAGDS(4,NODE)-CRTRFPR*DISPLC(JS,IS,KSM)
          DIAGDS(13,NODE)=DIAGDS(13,NODE)-CRTRFPR*DISPLC(JS,IS,KSM)
        ENDIF
        IF (IS.NE.1) THEN
          DIAGDS(2,NODE)=DIAGDS(2,NODE)-CRTRFPR*DISPLR(JS,IS,KSM)
          DIAGDS(11,NODE)=DIAGDS(11,NODE)-CRTRFPR*DISPLR(JS,IS,KSM)
        ENDIF
c orig lines
c        DIAGDS(5,NODE)=DIAGDS(5,NODE)-CRTRFPR*DISPLL(JS,IS,KSM)
c        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLL(JS,IS,KSM)
        IF(THCK(JS,IS,KSM).GT.0.0) THEN
          DIAGDS(5,NODE)=DIAGDS(5,NODE)-CRTRFPR*DISPLL(JS,IS,KSM)
c dzz fix
     X    *(THCK(JS,IS,KS)/THCK(JS,IS,KSM))
          DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLL(JS,IS,KSM)
c dzz fix
     X    *(THCK(JS,IS,KS)/THCK(JS,IS,KSM))
        ELSE
          DIAGDS(5,NODE)=DIAGDS(5,NODE)-CRTRFPR*DISPLL(JS,IS,KSM)
          DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLL(JS,IS,KSM)
        ENDIF
C
        IF (JS.NE.NSCOL) THEN
          DIAGDS(6,NODE)=DIAGDS(6,NODE)+CRTRFPR*DISPLC(JS,IS,KSM)
          DIAGDS(15,NODE)=DIAGDS(15,NODE)+CRTRFPR*DISPLC(JS,IS,KSM)
        ENDIF
        IF (IS.NE.NSROW) THEN
          DIAGDS(8,NODE)=DIAGDS(8,NODE)+CRTRFPR*DISPLR(JS,IS,KSM)
          DIAGDS(17,NODE)=DIAGDS(17,NODE)+CRTRFPR*DISPLR(JS,IS,KSM)
        ENDIF
      ENDIF
C
131   CONTINUE
C
      RETURN
      END
C
C
C
C
C     ELLDIS BUILD DISPERSION INTEGRAL 
C***********************************************************************
C
      SUBROUTINE ELLDISH(IBOUND,DISPCC,DISPCR,DISPCL,DISPRC,DISPRR,
     *               DISPRL,DISPLC,DISPLR,DISPLL,
     *      POR,THCK,RF,DIAGDS,TIMV,
     *      NCOL,NROW,NLAY,NSCOL,NSROW,NSLAY,NODESS,
     *      DELCOL,DELROW,LOCHFB)
C
C***********************************************************************
C   ELLDIS IS CALLED AT THE FIRST TRANSPORT TIME STEP TO USE DISPERSION
C   COEFFICIENTS CALCULATED BY DSP6 TO BUILD THE LHS DISPERSION MATRIX,
C   DIAGDS.  DIAGDS HAS THE SAME 27 DIAGONAL FORM AS THE STORAGE MATRIX,
C   AND LIKE DIAGST IS STORED AS A 27 X NODESS ARRAY.
C   AT THE FIRST SOLVE (IMOV=1) DIAGST AND DIAGDS ARE SUMMED TO PRODUCE
C   THE SYSTEM MATRIX FOR SOLUTION, AND CONVERTED TO THE SLAP COLUMN
C   FORMAT USED BY THE SOLVER.  THE CONVERTED FORM OF DIAGDS+DIAGST IS
C   THEN STORED IN DIAGDS FOR USE IN FUTURE TRANSPORT SOLVES DURING 
C   DURING SINGLE FLOW TIME STEP.
C   
      DIMENSION  IBOUND(NCOL,NROW,NLAY),DISPCC(NSCOL,NSROW,NSLAY),
     *    DISPCR(NSCOL,NSROW,NSLAY),
     *    DISPCL(NSCOL,NSROW,NSLAY),DISPRC(NSCOL,NSROW,NSLAY),
     *    DISPRR(NSCOL,NSROW,NSLAY),DISPRL(NSCOL,NSROW,NSLAY),
     *    DISPLC(NSCOL,NSROW,NSLAY),DISPLR(NSCOL,NSROW,NSLAY),
     *    DISPLL(NSCOL,NSROW,NSLAY),
     *    LOCHFB(NSCOL,NSROW,NSLAY,2)
      DIMENSION POR(NODESS),THCK(NODESS),RF(NSCOL,NSROW,NSLAY),
     *    DIAGDS(27,NODESS),DELCOL(NCOL),DELROW(NROW)
      COMMON /SUBGRD/ ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,
     *  ISLAY2,ISUBGD
      COMMON /ELLAM/ CINV,RINV,BINV,HCINV,HRINV,HBINV,WATVOL,
     *               STINIT,ADINIT,STMASS,ADMASS,OLMASS,
     *               AZERO,
     *               NSC,NSR,NSL,NT,NCTF,NRTF,NLTF,
     *               NEIGHB(8,2),NSLOPE(3),
     *               IDTOP,IDMAX,NCOEF,NSCH,NSRH,NSLH
C
C   DISPCC CONTAINS POR*THCK*DCC/CDEL
C   DISPCR CONTAINS POR*THCK*DCR/(4*RDEL)
C   DISPCL CONTAINS POR*THCK*DCL/(4*B)
C   DISPRC CONTAINS POR*THCK*DRC/(4*CDEL)
C   DISPRR CONTAINS POR*THCK*DRR/RDEL
C   DISPRL CONTAINS POR*THCK*DRL/(4*B)
C   DISPLC CONTAINS POR*DLC/(4*CDEL)
C   DISPLR CONTAINS POR*DLR/(4*RDEL)
C   DISPLL CONTAINS POR*DLL/B
C
          NSRC=NSROW*NSCOL
        DO 131 KS=1,NSLAY
          K=ISLAY1+KS-1
        DO 131 IS=1,NSROW
          I=ISROW1+IS-1
        DO 131 JS=1,NSCOL
          J=ISCOL1+JS-1
          NODE=(KS-1)*NSRC+(IS-1)*NSCOL+JS
C
      IF (IBOUND(J,I,K).EQ.0) THEN
        DIAGDS(14,NODE)=0D0
        GOTO 131
      ENDIF
C
      RFPINV=1.D0/(RF(JS,IS,KS)*POR(NODE))
      TRFPR=TIMV*RFPINV
C
      JSP=JS+1
      JSM=JS-1
      ISP=IS+1
      ISM=IS-1
      KSP=KS+1
      KSM=KS-1
C
C  ADD DISPCX COEFFS TO LHS
C
      RTRFPR=DELROW(I)*TRFPR
C
C  FORWARD COEFFS
c
c in current layer (js, is, ks is 14):
c 10 11 12
c 13 14 15
c 16 17 18
c
c check for HFB location: if exists, use appropriate node and 
c decrease grid factor by 2 (so a total of 2Y instead of 4Y in
c denominator)
c if doesn't exist, set as normal
      IF (JS.NE.NSCOL) THEN
        IF (KS.NE.1) THEN
          DIAGDS(5,NODE)=DIAGDS(5,NODE)+RTRFPR*DISPCL(JS,IS,KS)
          DIAGDS(6,NODE)=DIAGDS(6,NODE)+RTRFPR*DISPCL(JS,IS,KS)
        ENDIF
        IF (IS.NE.1) THEN
c check for hfb below 11 or 12
	    IF(LOCHFB(JS,IS-1,KS,2).EQ.1.OR.
     *	   LOCHFB(JS+1,IS-1,KS,2).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)+RTRFPR*DISPCR(JS,IS,KS)*2.0
            DIAGDS(15,NODE)=DIAGDS(15,NODE)+RTRFPR*DISPCR(JS,IS,KS)*2.0
          ELSE
            DIAGDS(11,NODE)=DIAGDS(11,NODE)+RTRFPR*DISPCR(JS,IS,KS)
            DIAGDS(12,NODE)=DIAGDS(12,NODE)+RTRFPR*DISPCR(JS,IS,KS)
          END IF
        ENDIF
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+RTRFPR*DISPCC(JS,IS,KS)
        DIAGDS(15,NODE)=DIAGDS(15,NODE)-RTRFPR*DISPCC(JS,IS,KS)
        IF (IS.NE.NSROW) THEN
c check for hfb below 14 or 15
	    IF(LOCHFB(JS,IS,KS,2).EQ.1.OR.
     *	   LOCHFB(JS+1,IS,KS,2).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)-RTRFPR*DISPCR(JS,IS,KS)*2.0
            DIAGDS(15,NODE)=DIAGDS(15,NODE)-RTRFPR*DISPCR(JS,IS,KS)*2.0
          ELSE
            DIAGDS(17,NODE)=DIAGDS(17,NODE)-RTRFPR*DISPCR(JS,IS,KS)
            DIAGDS(18,NODE)=DIAGDS(18,NODE)-RTRFPR*DISPCR(JS,IS,KS)
          END IF
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(23,NODE)=DIAGDS(23,NODE)-RTRFPR*DISPCL(JS,IS,KS)
          DIAGDS(24,NODE)=DIAGDS(24,NODE)-RTRFPR*DISPCL(JS,IS,KS)
        ENDIF
      ENDIF
C  BACKWARD COEFFS
      IF (JS.NE.1) THEN
        IF (KS.NE.1) THEN
          DIAGDS(4,NODE)=DIAGDS(4,NODE)-RTRFPR*DISPCL(JSM,IS,KS)
          DIAGDS(5,NODE)=DIAGDS(5,NODE)-RTRFPR*DISPCL(JSM,IS,KS)
        ENDIF
        IF (IS.NE.1) THEN
c check for hfb below 10 or 11
	    IF(LOCHFB(JS-1,IS-1,KS,2).EQ.1.OR.
     *	   LOCHFB(JS,IS-1,KS,2).EQ.1) THEN
            DIAGDS(13,NODE)=DIAGDS(13,NODE)-RTRFPR*DISPCR(JSM,IS,KS)*2.0
            DIAGDS(14,NODE)=DIAGDS(14,NODE)-RTRFPR*DISPCR(JSM,IS,KS)*2.0
          ELSE
            DIAGDS(10,NODE)=DIAGDS(10,NODE)-RTRFPR*DISPCR(JSM,IS,KS)
            DIAGDS(11,NODE)=DIAGDS(11,NODE)-RTRFPR*DISPCR(JSM,IS,KS)
          END IF
        ENDIF
        DIAGDS(13,NODE)=DIAGDS(13,NODE)-RTRFPR*DISPCC(JSM,IS,KS)
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+RTRFPR*DISPCC(JSM,IS,KS)
        IF (IS.NE.NSROW) THEN
c check for hfb below 13 or 14
	    IF(LOCHFB(JS-1,IS,KS,2).EQ.1.OR.
     *	   LOCHFB(JS,IS,KS,2).EQ.1) THEN
            DIAGDS(13,NODE)=DIAGDS(13,NODE)+RTRFPR*DISPCR(JSM,IS,KS)*2.0
            DIAGDS(14,NODE)=DIAGDS(14,NODE)+RTRFPR*DISPCR(JSM,IS,KS)*2.0
          ELSE
            DIAGDS(16,NODE)=DIAGDS(16,NODE)+RTRFPR*DISPCR(JSM,IS,KS)
            DIAGDS(17,NODE)=DIAGDS(17,NODE)+RTRFPR*DISPCR(JSM,IS,KS)
          END IF
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(22,NODE)=DIAGDS(22,NODE)+RTRFPR*DISPCL(JSM,IS,KS)
          DIAGDS(23,NODE)=DIAGDS(23,NODE)+RTRFPR*DISPCL(JSM,IS,KS)
        ENDIF
      ENDIF
C
C  ADD DISPRX COEFFS TO LHS
C
      CTRFPR=DELCOL(J)*TRFPR
C
C  FORWARD COEFFS
      IF (IS.NE.NSROW) THEN
        IF (KS.NE.1) THEN
          DIAGDS(5,NODE)=DIAGDS(5,NODE)+CTRFPR*DISPRL(JS,IS,KS)
          DIAGDS(8,NODE)=DIAGDS(8,NODE)+CTRFPR*DISPRL(JS,IS,KS)
        ENDIF
        IF (JS.NE.1) THEN
c check for hfb to right of 13 or 16
	    IF(LOCHFB(JS-1,IS,KS,1).EQ.1.OR.
     *	   LOCHFB(JS-1,IS+1,KS,1).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)+CTRFPR*DISPRC(JS,IS,KS)*2.0
            DIAGDS(17,NODE)=DIAGDS(17,NODE)+CTRFPR*DISPRC(JS,IS,KS)*2.0
          ELSE
            DIAGDS(13,NODE)=DIAGDS(13,NODE)+CTRFPR*DISPRC(JS,IS,KS)
            DIAGDS(16,NODE)=DIAGDS(16,NODE)+CTRFPR*DISPRC(JS,IS,KS)
          END IF
        ENDIF
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CTRFPR*DISPRR(JS,IS,KS)
        DIAGDS(17,NODE)=DIAGDS(17,NODE)-CTRFPR*DISPRR(JS,IS,KS)
        IF (JS.NE.NSCOL) THEN
c check for hfb to right of 14 or 17
	    IF(LOCHFB(JS,IS,KS,1).EQ.1.OR.
     *	   LOCHFB(JS,IS+1,KS,1).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)-CTRFPR*DISPRC(JS,IS,KS)*2.0
            DIAGDS(17,NODE)=DIAGDS(17,NODE)-CTRFPR*DISPRC(JS,IS,KS)*2.0
          ELSE 
		  DIAGDS(15,NODE)=DIAGDS(15,NODE)-CTRFPR*DISPRC(JS,IS,KS)
            DIAGDS(18,NODE)=DIAGDS(18,NODE)-CTRFPR*DISPRC(JS,IS,KS)
          END IF
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(23,NODE)=DIAGDS(23,NODE)-CTRFPR*DISPRL(JS,IS,KS)
          DIAGDS(26,NODE)=DIAGDS(26,NODE)-CTRFPR*DISPRL(JS,IS,KS)
        ENDIF
      ENDIF
C  BACKWARD COEFFS
      IF (IS.NE.1) THEN
        IF (KS.NE.1) THEN
          DIAGDS(2,NODE)=DIAGDS(2,NODE)-CTRFPR*DISPRL(JS,ISM,KS)
          DIAGDS(5,NODE)=DIAGDS(5,NODE)-CTRFPR*DISPRL(JS,ISM,KS)
        ENDIF
        IF (JS.NE.1) THEN
c check for hfb to right of 10 or 13
	    IF(LOCHFB(JS-1,IS-1,KS,1).EQ.1.OR.
     *	   LOCHFB(JS-1,IS,KS,1).EQ.1) THEN
            DIAGDS(11,NODE)=DIAGDS(11,NODE)-CTRFPR*DISPRC(JS,ISM,KS)*2.0
            DIAGDS(14,NODE)=DIAGDS(14,NODE)-CTRFPR*DISPRC(JS,ISM,KS)*2.0
          ELSE 
            DIAGDS(10,NODE)=DIAGDS(10,NODE)-CTRFPR*DISPRC(JS,ISM,KS)
            DIAGDS(13,NODE)=DIAGDS(13,NODE)-CTRFPR*DISPRC(JS,ISM,KS)
          END IF
        ENDIF
        DIAGDS(11,NODE)=DIAGDS(11,NODE)-CTRFPR*DISPRR(JS,ISM,KS)
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CTRFPR*DISPRR(JS,ISM,KS)
        IF (JS.NE.NSCOL) THEN
c check for hfb to right of 11 or 14
	    IF(LOCHFB(JS,IS-1,KS,1).EQ.1.OR.
     *	   LOCHFB(JS,IS,KS,1).EQ.1) THEN
            DIAGDS(11,NODE)=DIAGDS(11,NODE)+CTRFPR*DISPRC(JS,ISM,KS)*2.0
            DIAGDS(14,NODE)=DIAGDS(14,NODE)+CTRFPR*DISPRC(JS,ISM,KS)*2.0
          ELSE 
            DIAGDS(12,NODE)=DIAGDS(12,NODE)+CTRFPR*DISPRC(JS,ISM,KS)
            DIAGDS(15,NODE)=DIAGDS(15,NODE)+CTRFPR*DISPRC(JS,ISM,KS)
          END IF
        ENDIF
        IF (KS.NE.NSLAY) THEN
          DIAGDS(20,NODE)=DIAGDS(20,NODE)+CTRFPR*DISPRL(JS,ISM,KS)
          DIAGDS(23,NODE)=DIAGDS(23,NODE)+CTRFPR*DISPRL(JS,ISM,KS)
        ENDIF
      ENDIF
C
C  ADD DISPLX COEFFS TO LHS
C
      CRTRFPR=DELCOL(J)*DELROW(I)*TRFPR
C
C  FORWARD COEFFS
      IF (KS.NE.NSLAY) THEN
        IF (JS.NE.1) THEN
c check for hfb to right of 22 or 13
	    IF(LOCHFB(JS-1,IS,KS+1,1).EQ.1.OR.
     *	   LOCHFB(JS-1,IS,KS,1).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLC(JS,IS,KS)*2.0
            DIAGDS(23,NODE)=DIAGDS(23,NODE)+CRTRFPR*DISPLC(JS,IS,KS)*2.0
          ELSE 
            DIAGDS(13,NODE)=DIAGDS(13,NODE)+CRTRFPR*DISPLC(JS,IS,KS)
            DIAGDS(22,NODE)=DIAGDS(22,NODE)+CRTRFPR*DISPLC(JS,IS,KS)
          END IF
        ENDIF
        IF (IS.NE.1) THEN
c check for hfb below 20 or 11
	    IF(LOCHFB(JS,IS-1,KS+1,2).EQ.1.OR.
     *	   LOCHFB(JS,IS-1,KS,2).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLR(JS,IS,KS)*2.0
            DIAGDS(23,NODE)=DIAGDS(23,NODE)+CRTRFPR*DISPLR(JS,IS,KS)*2.0
          ELSE 
            DIAGDS(11,NODE)=DIAGDS(11,NODE)+CRTRFPR*DISPLR(JS,IS,KS)
            DIAGDS(20,NODE)=DIAGDS(20,NODE)+CRTRFPR*DISPLR(JS,IS,KS)
          END IF
        ENDIF
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLL(JS,IS,KS)
        DIAGDS(23,NODE)=DIAGDS(23,NODE)-CRTRFPR*DISPLL(JS,IS,KS)
        IF (JS.NE.NSCOL) THEN
c check for hfb to right of 23 or 14
	    IF(LOCHFB(JS,IS,KS+1,1).EQ.1.OR.
     *	   LOCHFB(JS,IS,KS,1).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)-CRTRFPR*DISPLC(JS,IS,KS)*2.0
            DIAGDS(23,NODE)=DIAGDS(23,NODE)-CRTRFPR*DISPLC(JS,IS,KS)*2.0
          ELSE 
            DIAGDS(15,NODE)=DIAGDS(15,NODE)-CRTRFPR*DISPLC(JS,IS,KS)
            DIAGDS(24,NODE)=DIAGDS(24,NODE)-CRTRFPR*DISPLC(JS,IS,KS)
          END IF
        ENDIF
        IF (IS.NE.NSROW) THEN
c check for hfb below 23 or 14
	    IF(LOCHFB(JS,IS,KS+1,2).EQ.1.OR.
     *	   LOCHFB(JS,IS,KS,2).EQ.1) THEN
            DIAGDS(14,NODE)=DIAGDS(14,NODE)-CRTRFPR*DISPLR(JS,IS,KS)*2.0
            DIAGDS(23,NODE)=DIAGDS(23,NODE)-CRTRFPR*DISPLR(JS,IS,KS)*2.0
          ELSE 
            DIAGDS(17,NODE)=DIAGDS(17,NODE)-CRTRFPR*DISPLR(JS,IS,KS)
            DIAGDS(26,NODE)=DIAGDS(26,NODE)-CRTRFPR*DISPLR(JS,IS,KS)
          END IF
        ENDIF
      ENDIF
C  BACKWARD COEFFS
      IF (KS.NE.1) THEN
        IF (JS.NE.1) THEN
c check for hfb to right of 4 or 13
	    IF(LOCHFB(JS-1,IS,KS-1,1).EQ.1.OR.
     *	   LOCHFB(JS-1,IS,KS,1).EQ.1) THEN
           DIAGDS(14,NODE)=DIAGDS(14,NODE)-CRTRFPR*DISPLC(JS,IS,KSM)*2.0
           DIAGDS(23,NODE)=DIAGDS(23,NODE)-CRTRFPR*DISPLC(JS,IS,KSM)*2.0
          ELSE 
            DIAGDS(4,NODE)=DIAGDS(4,NODE)-CRTRFPR*DISPLC(JS,IS,KSM)
            DIAGDS(13,NODE)=DIAGDS(13,NODE)-CRTRFPR*DISPLC(JS,IS,KSM)
          END IF
        ENDIF
        IF (IS.NE.1) THEN
c check for hfb below 2 or 11
	    IF(LOCHFB(JS,IS-1,KS-1,2).EQ.1.OR.
     *    LOCHFB(JS,IS-1,KS,2).EQ.1) THEN
           DIAGDS(14,NODE)=DIAGDS(14,NODE)-CRTRFPR*DISPLR(JS,IS,KSM)*2.0
           DIAGDS(23,NODE)=DIAGDS(23,NODE)-CRTRFPR*DISPLR(JS,IS,KSM)*2.0
          ELSE 
            DIAGDS(2,NODE)=DIAGDS(2,NODE)-CRTRFPR*DISPLR(JS,IS,KSM)
            DIAGDS(11,NODE)=DIAGDS(11,NODE)-CRTRFPR*DISPLR(JS,IS,KSM)
          END IF
        ENDIF
        DIAGDS(5,NODE)=DIAGDS(5,NODE)-CRTRFPR*DISPLL(JS,IS,KSM)
        DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLL(JS,IS,KSM)
        IF (JS.NE.NSCOL) THEN
c check for hfb to right of 5 or 14
	    IF(LOCHFB(JS,IS,KS-1,1).EQ.1.OR.
     *	   LOCHFB(JS,IS,KS,1).EQ.1) THEN
           DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLC(JS,IS,KSM)*2.0
           DIAGDS(23,NODE)=DIAGDS(23,NODE)+CRTRFPR*DISPLC(JS,IS,KSM)*2.0
          ELSE 
            DIAGDS(6,NODE)=DIAGDS(6,NODE)+CRTRFPR*DISPLC(JS,IS,KSM)
            DIAGDS(15,NODE)=DIAGDS(15,NODE)+CRTRFPR*DISPLC(JS,IS,KSM)
          END IF
        ENDIF
        IF (IS.NE.NSROW) THEN
c check for hfb below 5 or 14
	    IF(LOCHFB(JS,IS,KS-1,2).EQ.1.OR.
     *	   LOCHFB(JS,IS,KS,2).EQ.1) THEN
           DIAGDS(14,NODE)=DIAGDS(14,NODE)+CRTRFPR*DISPLR(JS,IS,KSM)*2.0
           DIAGDS(23,NODE)=DIAGDS(23,NODE)+CRTRFPR*DISPLR(JS,IS,KSM)*2.0
          ELSE 
            DIAGDS(8,NODE)=DIAGDS(8,NODE)+CRTRFPR*DISPLR(JS,IS,KSM)
            DIAGDS(17,NODE)=DIAGDS(17,NODE)+CRTRFPR*DISPLR(JS,IS,KSM)
          END IF
        ENDIF
      ENDIF
C
131   CONTINUE
C
      RETURN
      END
C
C
      SUBROUTINE MASKXE(MRNO,xd_mask,xe_mask,NPCELL,
     X     NX,NY,NZ,NXY,NXYZ,NRN)
C.....Redefine xd_mask temp array (xe_mask) to include NPCELL=0
      INTEGER CIN,MRNO,NX,NY,NZ,NXYZ
      logical xd_mask,xe_mask
      DIMENSION MRNO(*)
      DIMENSION NPCELL(NX,NY,NZ)
      dimension xd_mask(0:NX+1,0:NY+1,0:NZ+1)
      dimension xe_mask(0:NX+1,0:NY+1,0:NZ+1)
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C...
      DO 340 M=1,NXYZ
         MA=MRNO(M)
C.....Decode M into I,J,K
          CALL MTOIJK(M,I,J,K,NX,NY)
C
          IF(NPCELL(I,J,K).LT.1) THEN
C Set xe_mask (received in CR5DSP as xd_mask) to false 
           xe_mask(I,J,K)=.FALSE.
          ENDIF
  340 CONTINUE
C
      RETURN
      END
C