CC
C***************************************************************
C
      SUBROUTINE SMOC6PCT(CONC,CINIT,CTNMO,
     *    KSTP,NSTP,KPER,NPER,IMOV,NMOV,TIMV,
     *    NSCOL,NSROW,NSLAY,IOUTS,JUNIT,SUMTCH,
     *    NPNTCL,ICONFM,NIUNIT,NODESS)
C
C     *******************************************************
C     PRINT AND RECORD GAMMA AND PERCENT CHANGE DATA
C     *******************************************************
C
C        SPECIFICATIONS
C     -------------------------------------------------------
      real pchngtot,gzchngto,gnchngto,mpchng,mgzchng,mgnchng,
     &  pvarsum,gzvarsum,gnvarsum,pchngvr,gzchngvr,gnchngvr
      real pchng,gammnzero,gammnnmo
      allocatable pchng(:,:,:),gamnzero(:,:,:),gamnnmo(:,:,:)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),CINIT(NSCOL,NSROW,NSLAY),
     *  CTNMO(NSCOL,NSROW,NSLAY),JUNIT(NIUNIT)
C
      allocate(pchng(nscol,nsrow,nslay),gamnzero(nscol,nsrow,nslay),
     &  gamnnmo(nscol,nsrow,nslay))
C     -------------------------------------------------------
C RETURN IF NO OUTPUT SELECTED
      IF(JUNIT(19).LE.0.AND.JUNIT(20).LE.0.AND.
     &   JUNIT(21).LE.0.AND.JUNIT(22).LE.0) RETURN
C CHECK FLAGS FOR OUTPUT
      IPRNT=0
      IF (NPNTCL.EQ.0.AND.KPER.EQ.NPER.AND.KSTP.EQ.NSTP.
     *AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-2.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
      IF (NPNTCL.EQ.-1.AND.IMOV.EQ.NMOV) IPRNT=1
      IF(NPNTCL.GE.1) THEN
         IF(MOD(IMOV,NPNTCL).EQ.0) IPRNT=1
      ENDIF
      IF (KPER.EQ.NPER.AND.KSTP.EQ.NSTP.AND.IMOV.EQ.NMOV) IPRNT=1
C SKIP IF NO OUTPUT
      IF (IPRNT.EQ.0) RETURN
C Initialize
      pchng=0
      gamnzero=0
      gamnnmo=0
	pchngtot=0.
      gzchngto=0.
      gnchngto=0.
      mpchng=0.
	mgzchng=0.
	mgnchng=0.
      pvarsum=0.
      gzvarsum=0.
      gnvarsum=0.
      pchngvr=0.
      gzchngvr=0.
      gnchngvr=0.
C
C FOR EACH LAYER: CALC AND PRINT GAMMA AND PERCENT CHANGE IF REQUESTED.
      DO 80 KS=1,NSLAY
C PRINT TO SEPARATE TEXT FILES
         IF(KS.EQ.1.AND.IMOV.EQ.1) THEN
	     IF(JUNIT(19).GT.0) THEN
		   WRITE(IOUTS,51) JUNIT(19)
C HEADER LINE
             WRITE(JUNIT(19),55) KS, IMOV, KSTP, KPER, SUMTCH
           ELSE
             WRITE(iouts,*) 'Percent change output not enabled'
           END IF
		 IF(JUNIT(20).GT.0) THEN
             WRITE(IOUTS,52) JUNIT(20)
             WRITE(JUNIT(20),56) KS, IMOV, KSTP, KPER, SUMTCH
           ELSE
             WRITE(iouts,*) 'Gamma0 output not enabled'
           ENDIF
		 IF(JUNIT(21).GT.0) THEN
             WRITE(IOUTS,53) JUNIT(21)
             WRITE(JUNIT(21),57) KS, IMOV, KSTP, KPER, SUMTCH
           ELSE
             WRITE(iouts,*) 'Gamma1 output not enabled'
           ENDIF
		 IF(JUNIT(22).GT.0) THEN
             WRITE(IOUTS,54) JUNIT(21)
           ELSE
             WRITE(iouts,*) 'Means and variances output not enabled'
           ENDIF
         ENDIF
  51  FORMAT('PERCENT CHANGE DATA WILL BE SAVED ON UNIT ',I3)
  52  FORMAT('GAMMA0 DATA WILL BE SAVED ON UNIT ',I3)
  53  FORMAT('GAMMA1 DATA WILL BE SAVED ON UNIT ',I3)
  54  FORMAT('MEANS AND VARIANCES DATA WILL BE SAVED ON UNIT ',I3)
C
  55  FORMAT('PERCENT CHANGE AT NODES IN SUBGRID LAYER ',I4,', IMOV='
     * ,I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
  56  FORMAT('GAMMA0         AT NODES IN SUBGRID LAYER ',I4,', IMOV='
     * ,I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
  57  FORMAT('GAMMA1         AT NODES IN SUBGRID LAYER ',I4,', IMOV='
     * ,I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
C
      DO 80 IS=1,NSROW               
       DO 60 JS=1,NSCOL
C  Calculate percent change
         IF(CTNMO(JS,IS,KS).GT.0) 
     &     pchng(JS,IS,KS)=(CONC(JS,IS,KS)-CTNMO(JS,IS,KS))/
     &     CTNMO(JS,IS,KS)*100
C  Sum pct change
         pchngtot=pchngtot+pchng(JS,IS,KS)
C  Calculate gamma w/ respect to initial concentration
         IF(CONC(JS,IS,KS).GE.CINIT(JS,IS,KS).AND.
     &     CINIT(JS,IS,KS).GT.0) THEN
	     gamnzero(JS,IS,KS)=(CONC(JS,IS,KS)-CINIT(JS,IS,KS))/
     &     CINIT(JS,IS,KS)
         ELSE IF (CONC(JS,IS,KS).GT.0) THEN
	     gamnzero(JS,IS,KS)=(CONC(JS,IS,KS)-CINIT(JS,IS,KS))/
     &     CONC(JS,IS,KS)
         END IF
C  Sum gamma w/ respect to initial concentration
         gzchngto=gzchngto+gamnzero(JS,IS,KS)
C  Calculate gamma w/ respect to previous time's concentration
         IF(CONC(JS,IS,KS).GE.CTNMO(JS,IS,KS).AND.
     &     CTNMO(JS,IS,KS).GT.0) THEN
	     gamnnmo(JS,IS,KS)=(CONC(JS,IS,KS)-CTNMO(JS,IS,KS))/
     &     CTNMO(JS,IS,KS)
         ELSE IF (CONC(JS,IS,KS).GT.0) THEN
	     gamnnmo(JS,IS,KS)=(CONC(JS,IS,KS)-CTNMO(JS,IS,KS))/
     &     CONC(JS,IS,KS)
         END IF
C  Sum gamma w/ respect to previous time's concentration
         gnchngto=gnchngto+gamnnmo(JS,IS,KS)
C  Save this time's concentration for next time we do this loop
         CTNMO(JS,IS,KS)=CONC(JS,IS,KS)
C  End column loop for calculation
  60   CONTINUE
C  Write output out in rows to output files
        IF(JUNIT(19).GT.0)
     *   WRITE(JUNIT(19),'(1P10E12.4)') (pchng(JS,IS,KS),JS=1,NSCOL)
        IF(JUNIT(20).GT.0)
     *   WRITE(JUNIT(20),'(1P10E12.4)') (gamnzero(JS,IS,KS),JS=1,NSCOL)
        IF(JUNIT(21).GT.0)
     *   WRITE(JUNIT(21),'(1P10E12.4)') (gamnnmo(JS,IS,KS),JS=1,NSCOL)
  80  CONTINUE
C
C  Calculate means and variances
      mpchng=pchngtot/nodess
	mgzchng=gzchngto/nodess
	mgnchng=gnchngto/nodess
C
      DO 100 KS=1,NSLAY
      DO 100 IS=1,NSROW
      DO 100 JS=1,NSCOL
         pvarsum=pvarsum+(pchng(js,is,ks)-mpchng)**2
         gzvarsum=gzvarsum+(gamnzero(js,is,ks)-mgzchng)**2
         gnvarsum=gnvarsum+(gamnnmo(js,is,ks)-mgnchng)**2
 100  CONTINUE
      pchngvr=pvarsum/(nodess-1)
      gzchngvr=gzvarsum/(nodess-1)
      gnchngvr=gnvarsum/(nodess-1)
C  Output means and variances
C 
      WRITE(JUNIT(22),58) IMOV, KSTP, KPER, SUMTCH
  58  FORMAT('MEANS AND VARS AT IMOV='
     * ,I5,', KSTP=',I5,', KPER=',I5,', SUMTCH=',1PE10.4)
      WRITE(JUNIT(22),*) 'Mean percent change =',mpchng
      WRITE(JUNIT(22),*) 'Mean gamma zero     =',mgzchng
      WRITE(JUNIT(22),*) 'Mean gamma n        =',mgnchng
      WRITE(JUNIT(22),*) 'Var percent change  =',pchngvr
      WRITE(JUNIT(22),*) 'Var gamma zero      =',gzchngvr
      WRITE(JUNIT(22),*) 'Var gamma n         =',gnchngvr
C
      deallocate(pchng,gamnzero,gamnnmo)
   50 RETURN
      END
C
C
C
C GWT1PCT5AL ALLOCATE SPACE FOR PCT CHANGE AND GAMMA OUTPUT
C
C     ******************************************************************
C
      SUBROUTINE GWT1PCT5AL(ISUM,LSCINI,LSCTNM,JUNIT,NIUNIT,IOUTS,
     *  NODESS)    
C
C     ******************************************************************
C
      DIMENSION JUNIT(NIUNIT)
C  CHECK FOR ANY OF THE PCT CHANGE OR GAMMA OUTPUT
      IF(JUNIT(19).GT.0.OR.JUNIT(20).GT.0.OR.
     &   JUNIT(21).GT.0.OR.JUNIT(22).GT.0) THEN
        LSCINI=ISUM
        ISUM=ISUM+NODESS
        LSCTNM=ISUM
        ISUM=ISUM+NODESS
        WRITE(IOUTS,101) NODESS*2
      ELSE
        LSCINI=1
        LSCTNM=1
      END IF
  101 FORMAT(1X,I8,' ELEMENTS IN IX ARRAY ARE USED BY PCT AND GAMMA')
      RETURN
      END
C
C