C     Last change:  Dec. 8, 2014
C
C   GWT1BFLX5DF -- READ BOUNDARY FLUX FLAGS
C     ****************************************************************
C
      SUBROUTINE GWT1BFLX5DF(INBFLX,INRECH,INEVT,INETS,
     *  IRCHTP,IEVTTP,NCHNDS,IOUTS,LSBFCH)
C
C READ FLAG FOR BOUNDARY OR DISTRIBUTED SOURCE FOR RCH AND EVT
      WRITE(IOUTS,*)
      IF(INBFLX.GT.0) THEN
        WRITE(IOUTS,*) 'BOUNDARY FLUX PACKAGE ACTIVATED '
        IF(INRECH.GT.0) THEN
          READ(INBFLX,*) IRCHTP
        ELSE
          WRITE(IOUTS,*) 'NO RECHARGE IN SIMULATION: IRCHTP NOT READ'
        END IF
        IF(INEVT.GT.0.OR.INETS.GT.0) THEN
          READ(INBFLX,*) IEVTTP
        ELSE
          WRITE(IOUTS,*) 'NO ET IN SIMULATION: IEVTTP NOT READ'
        END IF
        READ(INBFLX,*) NCHNDS
        IF(NCHNDS.LE.0) THEN
          SELECT CASE(NCHNDS)
          CASE(0)
            LSBFCH=1
          CASE(-1)
            WRITE(IOUTS,91) NCHNDS
          CASE(-2)
            WRITE(IOUTS,92) NCHNDS
          CASE(-3)
            WRITE(IOUTS,93) NCHNDS
          CASE(-4)
            WRITE(IOUTS,94) NCHNDS
	    END SELECT 
        ELSE
          WRITE(IOUTS,'(A,I8,A)') ' NCHNDS= ',NCHNDS,
     *  '; IFACE VALUES READ BELOW'
   	  END IF  
   91 FORMAT(' NCHNDS=',I4,'; SOURCE/SINK TERM FOR ALL CONSTANT-HEAD ',
     * /,'   NODES DISTRIBUTED UNIFORMLY ACROSS ALL NO-FLOW FACES 1-4')
   92 FORMAT(' NCHNDS=',I4,'; SOURCE/SINK TERM FOR ALL CONSTANT-HEAD ',
     * /,'   NODES DISTRIBUTED UNIFORMLY ACROSS ALL NO-FLOW FACES 1-6')
   93 FORMAT(' NCHNDS=',I4,'; SOURCE/SINK TERM FOR ALL CONSTANT-HEAD ',
     * /,'   NODES ASSIGNED TO TOP FACE OF UPPERMOST ACTIVE LAYER')
   94 FORMAT(' NCHNDS=',I4,'; SOURCE/SINK TERM FOR ALL CONSTANT-HEAD ',
     * /,'   NODES ASSIGNED TO BOTTOM FACE OF LOWERMOST ACTIVE LAYER')
      ELSE
        IRCHTP=0
        IEVTTP=0
        LSBFCH=1
        WRITE(IOUTS,*) 'BOUNDARY FLUX PACKAGE INACTIVE; ALL SOURCE/SINK 
     *TERMS FOR RECHARGE, ET,'
        WRITE(IOUTS,*) ' AND CONSTANT-HEAD CELLS ARE DISTRIBUTED '
      END IF
C
      IF(INRECH.GT.0) THEN
        WRITE(IOUTS,*) 'IRCHTP=',IRCHTP
        IF(IRCHTP.EQ.0) THEN
          WRITE(IOUTS,*) 'RECHARGE APPLIED AS DISTRIBUTED STRESS'
        ELSE
          WRITE(IOUTS,*) 'RECHARGE APPLIED AS BOUNDARY FLUX ON TOP FACE'
        END IF
      END IF
C
      IF(INEVT.GT.0) THEN
        WRITE(IOUTS,*) 'IEVTTP=',IEVTTP
        IF(IEVTTP.EQ.0) THEN
          WRITE(IOUTS,*) ' EVAPOTRANSPIRATION APPLIED AS DISTRIBUTED 
     *    STRESS'
        ELSE
          WRITE(IOUTS,*) ' EVAPOTRANSPIRATION APPLIED AS BOUNDARY FLUX 
     *    ON TOP FACE'
        END IF
      END IF
C
      WRITE(IOUTS,*)
      RETURN
      END
C
C
C GWT1BFLX5AL ALLOCATE SPACE FOR BOUNDARY FLUX CONSTANT-HEAD CELLS
C
C     ******************************************************************
C
      SUBROUTINE GWT1BFLX5AL(ISUMI,LSBFCH,INBFLX,NCHNDS,IOUTS)    
C
C     ******************************************************************
C
C     ARRAY IS LAY,ROW,COL,IFACE
      IF(INBFLX.GT.0.AND.NCHNDS.GT.0) THEN
        LSBFCH=ISUMI
        ISUMI=ISUMI+NCHNDS*4
        ISP=NCHNDS*4
        WRITE(IOUTS,101) ISP
      ELSE
        LSBFCH=1
      END IF
  101 FORMAT(1X,I8,' ELEMENTS IN IX ARRAY ARE USED BY BFLX')
      RETURN
      END
C
C
C  GWT1BFLX5RP READ BOUNDARY FLUX INFO FOR CONSTANT-HEAD CELLS
C
C     ******************************************************************
C
      SUBROUTINE GWT1BFLX5RP(BFCH,NCHNDS,IOUTS,INBFLX,IBOUND,
     *  NCOL,NROW,NLAY)
C
C     ******************************************************************
C
C     READ BOUNDARY FLUX INFO FOR CONSTANT-HEAD CELLS
C     ******************************************************************
C
      INTEGER BFCH
      DIMENSION BFCH(4,NCHNDS),IBOUND(NCOL,NROW,NLAY)               
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ******************************************************************
C
      WRITE(IOUTS,89) NCHNDS
   89 FORMAT(/,'NCHNDS=',I4,'; CONSTANT-HEAD CELLS WHERE BOUNDARY FLUX',
     *  /,' ASSIGNED ACCORDING TO IFACE VALUES')
      WRITE (IOUTS,150)
  150 FORMAT('   LAYER     ROW  COLUMN   IFACE')
C READ THE FIRST RECORD
      IERR=0
      DO 139 ICH=1,NCHNDS
         READ(INBFLX,*) BFCH(1,ICH),BFCH(2,ICH),
     *                 BFCH(3,ICH),BFCH(4,ICH)
         K=BFCH(1,ICH)
         I=BFCH(2,ICH)
         J=BFCH(3,ICH)
         IF(IBOUND(J,I,K).GE.0) THEN
           WRITE(IOUTS,'(A,3I4,A)') 'BFLX: CELL',J,I,K,
     *       ' IS NOT A CONSTANT-HEAD CELL'
         END IF
         IFACE=BFCH(4,ICH)
C
         WRITE(IOUTS,'(5I8)') K,I,J,IFACE            
         IF(K.LT.ISLAY1.OR.K.GT.ISLAY2.OR. 
     *      I.LT.ISROW1.OR.I.GT.ISROW2.OR. 
     *      J.LT.ISCOL1.OR.J.GT.ISCOL2) THEN
            WRITE(IOUTS,*) '***WARNING***  IFACE DEFINED FOR CONSTANT-',
c            WRITE(IOUTS,*) '***ERROR***  IFACE DEFINED FOR CONSTANT-',
     *                      'HEAD CELL OUTSIDE SUBGRID'
            IERR=1
         ENDIF
  139 CONTINUE
      WRITE(IOUTS,*) 
c      IF(IERR.EQ.1) THEN
c        WRITE(IOUTS,*) 'STOPPING BECAUSE IFACE LOCATION ERROR
c     & IN BFLX PACKAGE'
c        STOP 'ERROR IN BLFX PACKAGE'
c	END IF
      RETURN
      END       
C
C  GWT1BFLX5FM  APPLY CONSTANT-HEAD BOUNDARY FLUXES
C     ***************************************************************
C
      SUBROUTINE GWT1BFLX5FM(
     *  CHDFLO,IBOUND,
     *  NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,
     *  IOUTS,IFXHED,ICONLY,
     *  BFCH,NCHNDS,VC,VR,VL,THCK)
C
C     ***************************************************************
C
      INTEGER BFCH
C TEMPORARY ARRAY
      ALLOCATABLE IIFACE(:,:,:)
      DIMENSION IBOUND(NCOL,NROW,NLAY),CHDFLO(NSCOL,NSROW,NSLAY)
      DIMENSION VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1),BFCH(4,NCHNDS),THCK(NSCOL,NSROW,NSLAY) 
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C ALLOCATE TEMPORARY ARRAY
      ALLOCATE(IIFACE(NSCOL,NSROW,NSLAY))
C     ***************************************************************
C  RETURN IF ICONLY = 1, NO FLUID SINK/SOURCES WITHIN SUBGRID
      IF(ICONLY.EQ.1) THEN
	  DEALLOCATE(IIFACE)
	  RETURN
      END IF
C 
C  INITIALIZE TEMP ARRAY
      IIFACE=0
C
C  RETURN IF NO FIXED HEAD NODES WITHIN TRANSPORT SUBGRID
      IF(IFXHED.LE.0) RETURN
C  FILL TEMPORARY ARRAY WITH IFACE VALUES (BFLX PACKAGE)
      IF(NCHNDS.GT.0) THEN
        DO 82 ICH=1,NCHNDS
          K=BFCH(1,ICH)
          KS=K-ISLAY1+1
          I=BFCH(2,ICH)
          IS=I-ISROW1+1
          J=BFCH(3,ICH)
          JS=J-ISCOL1+1
          IFACE=BFCH(4,ICH)
          IF(IFACE.LT.0) IFACE=-1
          IF(K.LT.ISLAY1.OR.K.GT.ISLAY2.OR. 
     *       I.LT.ISROW1.OR.I.GT.ISROW2.OR. 
     *       J.LT.ISCOL1.OR.J.GT.ISCOL2) GO TO 82
          IIFACE(JS,IS,KS)=IFACE
   82   CONTINUE
      END IF        
C
      DO 31 KS=1,NSLAY
      K=KS+ISLAY1-1
      DO 21 IS=1,NSROW
      I=IS+ISROW1-1
      DO 11 JS=1,NSCOL
      J=JS+ISCOL1-1
C  SKIP IF CELL IS NOT FIXED HEAD NODE
      IF(IBOUND(J,I,K).GE.0) GO TO 11
      IF(NCHNDS.GT.0) THEN
	  IFACE=IIFACE(JS,IS,KS)
      ELSE
        IFACE=NCHNDS
      END IF
	SELECT CASE(IFACE)
C  ASSIGN FLUX TO BOTTOM FACE OF LOWERMOST ACTIVE LAYER
      CASE(-4)
        IBOT=1
C       CHECK FOR ALL INACTIVE CELLS BELOW
        IF(K.LT.NLAY) THEN
          DO 101 KK=K+1,NLAY
            IF(IBOUND(J,I,KK).NE.0) IBOT=0
  101     CONTINUE
        END IF
        IF(IBOT.EQ.1) THEN
          VL(JS,IS,KS+1)=VL(JS,IS,KS+1)-CHDFLO(JS,IS,KS)
          WRITE(IOUTS,'(2A,3I4)') ' BOUNDARY FLUX',
     * ' APPLIED TO BOTTOM FACE OF CONSTANT-HEAD CELL',J,I,K
c uncomment these lines for extra BFLX info
c        ELSE
c          WRITE(IOUTS,*)
c          WRITE(IOUTS,'(2A,3I4)') ' BOUNDARY FLUX NOT ',
c     * 'APPLIED TO CELL',J,I,K
        END IF
C  ASSIGN FLUX TO TOP FACE OF UPPERMOST LAYER
      CASE(-3)
        ITOP=1
C       CHECK FOR ALL INACTIVE CELLS ABOVE
        IF(K.GT.1) THEN
          DO 102 KK=1,K-1
            IF(IBOUND(J,I,KK).NE.0) ITOP=0
  102     CONTINUE
        END IF
        IF(ITOP.EQ.1) THEN
		VL(JS,IS,KS)=VL(JS,IS,KS)+CHDFLO(JS,IS,KS)
          WRITE(IOUTS,'(2A,3I4)') ' BOUNDARY FLUX',
     * ' APPLIED TO TOP FACE OF CONSTANT-HEAD CELL',J,I,K
c uncomment these lines for extra BFLX info
c        ELSE
c          WRITE(IOUTS,'(2A,3I4)') ' BOUNDARY FLUX NOT ',
c     * 'APPLIED TO CELL',J,I,K
        END IF
C  ASSIGN FLUX TO ALL FACES BORDERING INACTIVE CELLS
      CASE(-2)	
        I1=0
        I2=0
        I3=0
        I4=0
        I5=0
        I6=0
C INITIALIZE TOTAL AREA OF BOUNDARY FOR DENOMINATOR
        BOUNDA=0.0
        NFACE=0
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE OR 1D
        IF(NCOL.GT.1) THEN
          IF(J.GT.1) THEN
            IF(IBOUND(J-1,I,K).EQ.0) THEN
              I1=1
              BOUNDA=BOUNDA+CDEL*THCK(JS,IS,KS)
              NFACE=NFACE+1
            END IF 
          ELSE
            I1=1
            BOUNDA=BOUNDA+CDEL*THCK(JS,IS,KS)
            NFACE=NFACE+1
          END IF 
          IF(J.LT.NCOL) THEN
            IF(IBOUND(J+1,I,K).EQ.0) THEN
              I2=1
              BOUNDA=BOUNDA+CDEL*THCK(JS,IS,KS)
              NFACE=NFACE+1
            END IF 
          ELSE
            I2=1
            BOUNDA=BOUNDA+CDEL*THCK(JS,IS,KS)
            NFACE=NFACE+1
          END IF 
        END IF
        IF(NROW.GT.1) THEN
          IF(I.GT.1) THEN
            IF(IBOUND(J,I-1,K).EQ.0) THEN
              I3=1
              BOUNDA=BOUNDA+RDEL*THCK(JS,IS,KS)
              NFACE=NFACE+1
            END IF 
          ELSE
            I3=1
            BOUNDA=BOUNDA+RDEL*THCK(JS,IS,KS)
            NFACE=NFACE+1
          END IF 
          IF(I.LT.NROW) THEN
            IF(IBOUND(J,I+1,K).EQ.0) THEN
              I4=1
              BOUNDA=BOUNDA+RDEL*THCK(JS,IS,KS)
              NFACE=NFACE+1
            END IF 
          ELSE
            I4=1
            BOUNDA=BOUNDA+RDEL*THCK(JS,IS,KS)
            NFACE=NFACE+1
          END IF 
        END IF
        IF(NLAY.GT.1) THEN
          IF(K.GT.1) THEN
            IF(IBOUND(J,I,K-1).EQ.0) THEN
              I6=1
              BOUNDA=BOUNDA+RDEL*CDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I6=1
            BOUNDA=BOUNDA+RDEL*CDEL
            NFACE=NFACE+1
          END IF 
          IF(K.LT.NLAY) THEN
            IF(IBOUND(J,I,K+1).EQ.0) THEN
              I5=1
              NFACE=NFACE+1
              BOUNDA=BOUNDA+RDEL*CDEL
            END IF 
          ELSE
            I5=1
            NFACE=NFACE+1
            BOUNDA=BOUNDA+RDEL*CDEL
          END IF 
        END IF
C  
        IF(NFACE.EQ.0) THEN
          WRITE(IOUTS,9) J,I,K
          WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
   9      FORMAT('***WARNING*** IFACE VALUE FOR CONSTANT-HEAD CELL',3I4, 
     *       ' HAS NO NO-FLOW BOUNDARIES')
        ELSE
          RATE=CHDFLO(JS,IS,KS)
          IF(I1.EQ.1) THEN
	      FRACA=(CDEL*THCK(JS,IS,KS))/BOUNDA
		  VC(JS,IS,KS)=VC(JS,IS,KS)+RATE*FRACA                 
          END IF
		IF(I2.EQ.1) THEN
	      FRACA=(CDEL*THCK(JS,IS,KS))/BOUNDA
		  VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-RATE*FRACA                 
          END IF
		IF(I3.EQ.1) THEN
	      FRACA=(RDEL*THCK(JS,IS,KS))/BOUNDA
		  VR(JS,IS,KS)=VR(JS,IS,KS)+RATE*FRACA                 
          END IF
          IF(I4.EQ.1) THEN
	      FRACA=(RDEL*THCK(JS,IS,KS))/BOUNDA
		  VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-RATE*FRACA                 
          END IF
          IF(I6.EQ.1) THEN
	      FRACA=(CDEL*RDEL)/BOUNDA
		  VL(JS,IS,KS)=VL(JS,IS,KS)+RATE*FRACA                 
          END IF
          IF(I5.EQ.1) THEN
	      FRACA=(CDEL*RDEL)/BOUNDA
		  VL(JS,IS,KS+1)=VL(JS,IS,KS+1)-RATE*FRACA                 
          END IF
c          WRITE(IOUTS,'(2A,I2,A,3I4)') 'BOUNDARY FLUX',
c     * ' APPLIED TO',NFACE,' FACE(S) OF CONSTANT-HEAD CELL',J,I,K
        END IF
C  ASSIGN FLUX TO LATERAL FACES (1-4)
      CASE(-1)
        I1=0
        I2=0
        I3=0
        I4=0
C INITIALIZE TOTAL LENGTH OF BOUNDARY FOR DENOMINATOR
C (THICKNESS IS CONSTANT SO CAN BE LEFT OUT)
        BOUNDL=0.0
        NFACE=0
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE OR 1D
        IF(NCOL.GT.1) THEN
          IF(J.GT.1) THEN
            IF(IBOUND(J-1,I,K).EQ.0) THEN
              I1=1
              BOUNDL=BOUNDL+CDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I1=1
            BOUNDL=BOUNDL+CDEL
            NFACE=NFACE+1
          END IF 
          IF(J.LT.NCOL) THEN
            IF(IBOUND(J+1,I,K).EQ.0) THEN
              I2=1
              BOUNDL=BOUNDL+CDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I2=1
            BOUNDL=BOUNDL+CDEL
            NFACE=NFACE+1
          END IF 
        END IF
        IF(NROW.GT.1) THEN
          IF(I.GT.1) THEN
            IF(IBOUND(J,I-1,K).EQ.0) THEN
              I3=1
              BOUNDL=BOUNDL+RDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I3=1
            BOUNDL=BOUNDL+RDEL
            NFACE=NFACE+1
          END IF 
          IF(I.LT.NROW) THEN
            IF(IBOUND(J,I+1,K).EQ.0) THEN
              I4=1
              BOUNDL=BOUNDL+RDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I4=1
            BOUNDL=BOUNDL+RDEL
            NFACE=NFACE+1
          END IF 
        END IF
C  
        IF(NFACE.EQ.0) THEN
          WRITE(IOUTS,9) J,I,K
          WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
        ELSE
          RATE=CHDFLO(JS,IS,KS)
          IF(I1.EQ.1) THEN
	      FRACL=CDEL/BOUNDL
		  VC(JS,IS,KS)=VC(JS,IS,KS)+RATE*FRACL                 
          END IF
		IF(I2.EQ.1) THEN
	      FRACL=CDEL/BOUNDL
		  VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-RATE*FRACL               
          END IF
		IF(I3.EQ.1) THEN
	      FRACL=RDEL/BOUNDL
		  VR(JS,IS,KS)=VR(JS,IS,KS)+RATE*FRACL                
          END IF
		IF(I4.EQ.1) THEN
	      FRACL=RDEL/BOUNDL
		  VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-RATE*FRACL  
          END IF              
c          WRITE(IOUTS,'(2A,I2,A,3I4)') 'BOUNDARY FLUX',
c     * ' APPLIED TO',NFACE,' FACE(S) OF CONSTANT-HEAD CELL',J,I,K
        END IF
C IFACE=1
      CASE (1)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
        IF(J.GT.1) THEN
          IF(IBOUND(J-1,I,K).NE.0) THEN
            WRITE(IOUTS,111) J,I,K
            WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
          ELSE
            VC(JS,IS,KS)=VC(JS,IS,KS)+CHDFLO(JS,IS,KS)                          
          END IF 
        ELSE
          VC(JS,IS,KS)=VC(JS,IS,KS)+CHDFLO(JS,IS,KS)                          
        END IF 
C IFACE=2
      CASE (2)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
        IF(J.LT.NCOL) THEN
          IF(IBOUND(J+1,I,K).NE.0) THEN
            WRITE(IOUTS,111) J,I,K
            WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
          ELSE
            VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-CHDFLO(JS,IS,KS)                          
          END IF 
        ELSE
          VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-CHDFLO(JS,IS,KS)                          
        END IF 
C IFACE=3
      CASE (3)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
        IF(I.LT.NROW) THEN
          IF(IBOUND(J,I+1,K).NE.0) THEN
            WRITE(IOUTS,111) J,I,K
            WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
          ELSE
            VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-CHDFLO(JS,IS,KS)                          
          END IF 
        ELSE
          VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-CHDFLO(JS,IS,KS)                          
        END IF 
C IFACE=4
      CASE (4)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
        IF(I.GT.1) THEN
          IF(IBOUND(J,I-1,K).NE.0) THEN
            WRITE(IOUTS,111) J,I,K
            WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
          ELSE
            VR(JS,IS,KS)=VR(JS,IS,KS)+CHDFLO(JS,IS,KS)                          
          END IF 
        ELSE
          VR(JS,IS,KS)=VR(JS,IS,KS)+CHDFLO(JS,IS,KS)                          
        END IF 
C IFACE=5
      CASE (5)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
        IF(K.LT.NLAY) THEN
          IF(IBOUND(J,I,K+1).NE.0) THEN
            WRITE(IOUTS,111) J,I,K
            WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
          ELSE
            VL(JS,IS,KS+1)=VL(JS,IS,KS+1)-CHDFLO(JS,IS,KS)                          
          END IF 
        ELSE
          VL(JS,IS,KS+1)=VL(JS,IS,KS+1)-CHDFLO(JS,IS,KS)                          
        END IF 
C IFACE=6
      CASE (6)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
        IF(K.GT.1) THEN
          IF(IBOUND(J,I,K-1).NE.0) THEN
            WRITE(IOUTS,111) J,I,K
            WRITE(IOUTS,*) 'CONSTANT-HEAD CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
          ELSE
            VL(JS,IS,KS)=VL(JS,IS,KS)+CHDFLO(JS,IS,KS)                          
          END IF 
        ELSE
          VL(JS,IS,KS)=VL(JS,IS,KS)+CHDFLO(JS,IS,KS)                          
        END IF 
C
      END SELECT
  111   FORMAT('***WARNING*** IFACE VALUE FOR CONSTANT-HEAD CELL',3I4, 
     *       ' NOT ON NO-FLOW BOUNDARY')
   11 CONTINUE
   21 CONTINUE
   31 CONTINUE
C
      DEALLOCATE(IIFACE)
C
      RETURN
      END
C
C  GWT1BFLX5PCK  APPLY STRESS PACKAGE BOUNDARY FLUXES
C     ***************************************************************
C
      SUBROUTINE GWT1BFLX5PCK(
     *  IFACE,IBOUND,RATE,PACKAGE,IPCK,
     *  VC,VR,VL,
     *  NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,
     *  IOUTS,
     *  J,I,K,JS,IS,KS)
C
C     ***************************************************************
C
      CHARACTER*8 PACKAGE
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *   VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1)
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
C
C     ***************************************************************
            SELECT CASE (IFACE)
C IFACE<0
            CASE (-1)
C ASSIGN FLUX TO LATERAL FACES (1-4)
              I1=0
              I2=0
              I3=0
              I4=0
              NFACE=0
C INITIALIZE TOTAL LENGTH OF BOUNDARY FOR DENOMINATOR
C (THICKNESS IS CONSTANT SO CAN BE LEFT OUT)
        BOUNDL=0.0
        NFACE=0
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE OR 1D
        IF(NCOL.GT.1) THEN
          IF(J.GT.1) THEN
            IF(IBOUND(J-1,I,K).EQ.0) THEN
              I1=1
              BOUNDL=BOUNDL+CDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I1=1
            BOUNDL=BOUNDL+CDEL
            NFACE=NFACE+1
          END IF 
          IF(J.LT.NCOL) THEN
            IF(IBOUND(J+1,I,K).EQ.0) THEN
              I2=1
              BOUNDL=BOUNDL+CDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I2=1
            BOUNDL=BOUNDL+CDEL
            NFACE=NFACE+1
          END IF 
        END IF
        IF(NROW.GT.1) THEN
          IF(I.GT.1) THEN
            IF(IBOUND(J,I-1,K).EQ.0) THEN
              I3=1
              BOUNDL=BOUNDL+RDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I3=1
            BOUNDL=BOUNDL+RDEL
            NFACE=NFACE+1
          END IF 
          IF(I.LT.NROW) THEN
            IF(IBOUND(J,I+1,K).EQ.0) THEN
              I4=1
              BOUNDL=BOUNDL+RDEL
              NFACE=NFACE+1
            END IF 
          ELSE
            I4=1
            BOUNDL=BOUNDL+RDEL
            NFACE=NFACE+1
          END IF 
        END IF
C  
        IF(NFACE.EQ.0) THEN
          WRITE(IOUTS,9) PACKAGE,IPCK
          WRITE(IOUTS,'(A,A,A)') PACKAGE,' CELL TREATED AS DISTRIBUTED',
     *             ' SOURCE/SINK TERM'
   9      FORMAT('***WARNING*** IFACE VALUE FOR ',A,I4, 
     *       ' HAS NO NO-FLOW BOUNDARIES')
        ELSE
          IWARN=0
          IF(I1.EQ.1) THEN
	      FRACL=CDEL/BOUNDL
		  VC(JS,IS,KS)=VC(JS,IS,KS)+RATE*FRACL                 
          END IF
		IF(I2.EQ.1) THEN
	      FRACL=CDEL/BOUNDL
		  VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-RATE*FRACL               
          END IF
		IF(I3.EQ.1) THEN
	      FRACL=RDEL/BOUNDL
		  VR(JS,IS,KS)=VR(JS,IS,KS)+RATE*FRACL                
          END IF
		IF(I4.EQ.1) THEN
	      FRACL=RDEL/BOUNDL
		  VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-RATE*FRACL  
          END IF              
c          WRITE(IOUTS,'(2A,I2,A,A8,A,3I4)') ' BOUNDARY FLUX',
c     * ' APPLIED TO',NFACE,' FACE(S) OF',PACKAGE,' CELL',J,I,K
        END IF
C
C IFACE=1
            CASE (1)
C DO NOT ASSIGN FLUX TO BOUNDARY IF NOT A NO-FLOW FACE
              IF(J.GT.1) THEN
                IF(IBOUND(J-1,I,K).NE.0) THEN
                   IWARN=1
                ELSE
                   VC(JS,IS,KS)=VC(JS,IS,KS)+RATE                          
                END IF 
              ELSE
                VC(JS,IS,KS)=VC(JS,IS,KS)+RATE                          
              END IF 
C IFACE=2
            CASE (2)
              IF(J.LT.NCOL) THEN
                IF(IBOUND(J+1,I,K).NE.0) THEN
                   IWARN=1
                ELSE
                   VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-RATE                          
                END IF 
              ELSE
                VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-RATE                          
              END IF 
C IFACE=3
            CASE (3)
              IF(I.LT.NROW) THEN
                IF(IBOUND(J,I+1,K).NE.0) THEN
                   IWARN=1
                ELSE
                   VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-RATE                          
                END IF 
              ELSE
                VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-RATE                          
              END IF 
C IFACE=4
            CASE (4)
              IF(I.GT.1) THEN
                IF(IBOUND(J,I-1,K).NE.0) THEN
                   IWARN=1
                ELSE
                   VR(JS,IS,KS)=VR(JS,IS,KS)+RATE                          
                END IF 
              ELSE
                VR(JS,IS,KS)=VR(JS,IS,KS)+RATE                          
              END IF 
C IFACE=5
            CASE (5)
              IF(K.LT.NLAY) THEN
                IF(IBOUND(J,I,K+1).NE.0) THEN
                   IWARN=1
                ELSE
                   VL(JS,IS,KS+1)=VL(JS,IS,KS+1)-RATE            
                END IF 
              ELSE
                VL(JS,IS,KS+1)=VL(JS,IS,KS+1)-RATE            
              END IF 
C IFACE=6
            CASE (6)
              IF(K.GT.1) THEN
                IF(IBOUND(J,I,K-1).NE.0) THEN
                   IWARN=1
                ELSE
                   VL(JS,IS,KS)=VL(JS,IS,KS)+RATE                          
                END IF 
              ELSE
                VL(JS,IS,KS)=VL(JS,IS,KS)+RATE                          
              END IF 
            END SELECT
C
            IF(IWARN.EQ.1) THEN
              WRITE(IOUTS,11) PACKAGE,IPCK
              WRITE(IOUTS,'(A,A8,A)') 'THIS ',PACKAGE,' CELL TREATED AS 
     *DISTRIBUTED SOURCE/SINK TERM'         
            END IF
C
   11       FORMAT('***WARNING*** IFACE VALUE FOR ',A,I4, 
     *       ' NOT ON NO-FLOW BOUNDARY')
C
      RETURN
      END
C
C
      SUBROUTINE GWT1BFLX5LK(
     *  IFACE,RATE,PACKAGE,IPCK,
     *  VC,VR,VL,
     *  NCOL,NROW,NLAY,
     *  NSCOL,NSROW,NSLAY,
     *  IOUTS,
     *  J,I,K,JS,IS,KS)
C
C     ***************************************************************
C
      CHARACTER*8 PACKAGE
      DIMENSION IBOUND(NCOL,NROW,NLAY),
     *   VC(NSCOL+1,NSROW,NSLAY),VR(NSCOL,NSROW+1,NSLAY),
     *   VL(NSCOL,NSROW,NSLAY+1)
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
C
C     ***************************************************************
            IWARN=1
            SELECT CASE (IFACE)
C IFACE=1
            CASE (1)
              VC(JS,IS,KS)=VC(JS,IS,KS)+RATE
			IWARN=0                          
C IFACE=2
            CASE (2)
              VC(JS+1,IS,KS)=VC(JS+1,IS,KS)-RATE
			IWARN=0                          
C IFACE=3
            CASE (3)
              VR(JS,IS+1,KS)=VR(JS,IS+1,KS)-RATE                          
			IWARN=0                          
C IFACE=4
            CASE (4)
              VR(JS,IS,KS)=VR(JS,IS,KS)+RATE                          
			IWARN=0                          
C IFACE=6
            CASE (6)
              VL(JS,IS,KS)=VL(JS,IS,KS)+RATE                          
			IWARN=0                          
            END SELECT
C
            IF(IWARN.EQ.1) THEN
              WRITE(IOUTS,11) PACKAGE,IPCK
              WRITE(IOUTS,'(A,A8,A)') 'THIS ',PACKAGE,' CELL TREATED AS 
     *DISTRIBUTED SOURCE/SINK TERM'         
            END IF
C
   11       FORMAT('***WARNING*** IFACE VALUE FOR ',A,I4, 
     *       ' NOT VALID')
C
      RETURN
      END
