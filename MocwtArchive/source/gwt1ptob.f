C  PTOB  OBSERVATION OF PARTICLE DISTRIBUTIONS AT SPECIFIED CELLS
C  8/95
C
C  GWT1PTOB5DF READ NUMBER OF PARTICLE OBSERVATIONS 
C
C     ******************************************************************
C
      SUBROUTINE GWT1PTOB5DF(NUMPTOB,NUMPTOB_MNW,
     *                  IOUTS,INPTOB)
C
C     ******************************************************************
C
C     READ IN NUMPTOB
C     ******************************************************************
C
      READ(INPTOB,*) NUMPTOB, NUMPTOB_MNW
      WRITE(IOUTS,*) 'READING ',NUMPTOB,' NON-MNW PARTICLE OBSERVATION
     * LOCATIONS'
      WRITE(IOUTS,*) 'READING ',NUMPTOB_MNW,' MNW PARTICLE OBSERVATION
     * LOCATIONS'
      RETURN
      END
C
C
C GWT1PTOB5AL ALLOCATE SPACE FOR PARTICLE OBSERVATIONS
C
C     ******************************************************************
C
      SUBROUTINE GWT1PTOB5AL(ISUM,LSPTOB,LSPTMN,LSPTUN,NUMPTOB,
     *                  NUMPTOB_MNW,NSCOL,NSROW,NSLAY,
     *                  IOUT,IOUTS,INUNITPTOB)    
C
C     ******************************************************************
C
CMOCWT
      INCLUDE 'ptwt.inc'
C     ARRAY IS LAY,ROW,COL,UNIT#
      IF(INUNITPTOB.EQ.0) THEN
        LSPTOB=1
	  LSPTMN=1
	  LSPTUN=1
        NUMPTOB=0
        NUMPTOB_MNW=0
        RETURN
      END IF
      IF(NUMPTOB.GT.0) THEN
        LSPTOB=ISUM
        ISUM=ISUM+NUMPTOB*4
        ISTE=NUMPTOB*4
      ELSE
        LSPTOB=1
        ISTE=0
	END IF
      IF(NUMPTOB_MNW.GT.0) THEN
        LSPTMN=ISUM
        ISUM=ISUM+NUMPTOB_MNW*2
        ISTE=ISTE+NUMPTOB_MNW*2
        LSPTUN=ISUM
        ISUM=ISUM+NSCOL*NSROW*NSLAY
        ISTE=ISTE+NSCOL*NSROW*NSLAY
      ELSE
        LSPTMN=1
        LSPTUN=1
	END IF
      WRITE(IOUTS,101) ISTE
  101 FORMAT(1X,I8,' ELEMENTS IN X ARRAY ARE USED BY PTOB')
C
      RETURN
      END
C
C
C  PTOB5RP READ PARTICLE OBSERVATION INPUT FILE 
C
C     ******************************************************************
C
      SUBROUTINE GWT1PTOB5RP(IPTOBLST,IPTMNLST,NUMPTOB,NUMPTOB_MNW,
     *                  MNWsite,mxwel2,nwell2,IPTOBMNW,well2,
     *                  PTOBLST,ncol,nrow,nlay,nscol,nsrow,nslay,
     *                  IOUTS,INPTOB)
C
C     ******************************************************************
C
C     READ PARTICLE OBSERVATION LOCATIONS
C     ******************************************************************
C
      double precision well2
      CHARACTER*32 MNWsite,PTOBLST,Msite
      DIMENSION IPTOBLST(4,NUMPTOB),PTOBLST(NUMPTOB_MNW)               
      DIMENSION IPTMNLST(NUMPTOB_MNW,2),well2(18,mxwel2)             
      dimension MNWsite(mxwel2),IPTOBMNW(NSCOL,NSROW,NSLAY)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ******************************************************************
C
C INITIALIZE
      IPTOBLST=0
      IPTMNLST=0
      IPTOBMNW=0
C
      IF(NUMPTOB.gt.0) THEN
      WRITE (IOUTS,140) NUMPTOB
      WRITE (IOUTS,150)
      end if
C READ THE NUMPTOB RECORDS
      DO 139 IPTOB=1,NUMPTOB
         READ(INPTOB,*) IPTOBLST(1,IPTOB),IPTOBLST(2,IPTOB),
     *                 IPTOBLST(3,IPTOB),IPTOBLST(4,IPTOB)
         IPTOBUN=IPTOBLST(4,IPTOB)
         K=IPTOBLST(1,IPTOB)
         I=IPTOBLST(2,IPTOB)
         J=IPTOBLST(3,IPTOB)
C
         WRITE(IOUTS,'(5I8,5X,A40)') IPTOB,K,I,J,IPTOBUN            
         IF(K.LT.ISLAY1.OR.K.GT.ISLAY2.OR. 
     *      I.LT.ISROW1.OR.I.GT.ISROW2.OR. 
     *      J.LT.ISCOL1.OR.J.GT.ISCOL2) THEN
            WRITE(IOUTS,*) '***ERROR***   PARTICLE OBSERVATION',
     *                      ' LOCATION OUTSIDE SUBGRID'
            STOP
         ENDIF
C
  139 CONTINUE
C  WRITE NUMPTOB HEADERS
      DO 149 IPTOB=1,NUMPTOB
         IHEADER=0
         IPTOBUN=IPTOBLST(4,IPTOB)
C  WRITE HEADER IF NOT ALREADY IN FILE
         DO 11 IPTE=1,NUMPTOB
C  CYCLE THROUGH UNIT NUMBERS, LOOK FOR MATCH IN PRECEEDING OBS
           IF(IPTE.LT.IPTOB) THEN
C  IF MATCH FOUND, SET HEADER FLAG TO TRUE
             IF(IPTOBUN.EQ.IPTOBLST(4,IPTE)) IHEADER=1
           END IF
   11    CONTINUE
         IF(IHEADER.EQ.0) THEN
           WRITE(IPTOBUN,'(256A)') '"PARTICLE OBSERVATION DATA"'
           WRITE(IPTOBUN,'(256A)') '"  J,  I,  K,         TIME          
     &PTWT        PCONC.      QSINK       QMNWSINK    
     &WT-REMV     MASS_OUT"'
         END IF
  149 CONTINUE
C
C READ THE NUMPTOB_MNW RECORDS
c and store them
      DO 239 IOB=1,NUMPTOB_MNW
         IS_SITE=0
         READ(INPTOB,*) PTOBLST(IOB),IPTMNLST(IOB,1)
         call UPCASE(PTOBLST(IOB))
c initialize header flag
         IPTMNLST(IOB,2)=0
C
  239 CONTINUE
      IF(NUMPTOB.gt.0) THEN
            WRITE(IOUTS,'(/)')
            WRITE(IOUTS,*) 'NON-MNW PARTICLE OBSERVATION DATA WILL BE'
     *                     ,' WRITTEN ON UNIT NUMBERS LISTED ABOVE' 
      END IF
  140 FORMAT(///'COORDINATES FOR',I4,' NON-MNW PARTICLE OBSERVATIONS:')
  150 FORMAT(/'  OBS #   LAYER     ROW  COLUMN    ',
     *                           'UNIT')
      RETURN
      END       
C
C
C  PTOB5RPS DEFINE PTOB MNW LOCATIONS FOR THIS STRESS PERIOD
C
C     ******************************************************************
C
      SUBROUTINE GWT1PTOB5RPS(IPTMNLST,NUMPTOB_MNW,
     *                  MNWsite,mxwel2,nwell2,IPTOBMNW,well2,
     *                  PTOBLST,ncol,nrow,nlay,nscol,nsrow,nslay,
     *                  IOUTS,INPTOB)
C
C     ******************************************************************
C
C     DEFINE PTOB MNW LOCATIONS FOR THIS STRESS PERIOD
C     ******************************************************************
C
      double precision well2
      CHARACTER*32 MNWsite,PTOBLST,Msite,TEMPSITE
      DIMENSION PTOBLST(NUMPTOB_MNW)               
      DIMENSION IPTMNLST(NUMPTOB_MNW,2),well2(18,mxwel2)             
      dimension MNWsite(mxwel2),IPTOBMNW(NSCOL,NSROW,NSLAY)
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ******************************************************************
c initialize array that carries active mnws
      IPTOBMNW=0
c write info line to output file
C
      IF(NUMPTOB_MNW.gt.0) THEN
      WRITE (IOUTS,140) NUMPTOB_MNW
      WRITE (IOUTS,150)
      end if
  140 FORMAT(/'SITE-IDS FOR',I4,' MNW PARTICLE OBSERVATIONS:')
  150 FORMAT(/'  SITE-ID          UNIT')

c each stress period,
c check site vs list of site names in MNWSITE
      DO 239 IOB=1,NUMPTOB_MNW
         IS_SITE=0
         TEMPSITE=PTOBLST(IOB)
c Loop over all MNW locations
         m = 0
         do while( m .lt. nwell2 )
          m = m + 1
          MSITE=MNWSITE(m)
          call UPCASE(MSITE)
          IF(TEMPSITE.EQ.MSITE) THEN
		  IS_SITE=1
C Define IPTOBMNW, which stores the MNW locations as non-zero unit #s in each cell
c get node location 
            n = INT(well2(1,m))
            k = (n-1)/(ncol*nrow) + 1
            j = mod((n-1),ncol*nrow)/ncol + 1
            i = mod((n-1),ncol) + 1
            JS=I-ISCOL1+1
            IS=J-ISROW1+1
            KS=K-ISLAY1+1
            IPTOBMNW(JS,IS,KS)=IPTMNLST(IOB,1)
          END IF
         end do      
C
         WRITE(IOUTS,'(I8,3X,A12,I8)') IOB,TEMPSITE,IPTMNLST(IOB,1)        
         IF(IS_SITE.EQ.0) THEN
            WRITE(IOUTS,*) '  MNW PARTICLE OBSERVATION SITE ',TEMPSITE,
     *                      '  NOT ACTIVE THIS STRESS PERIOD'
            WRITE(IOUTS,*) 
         ENDIF
  239 CONTINUE
C  WRITE NUMPTOB_MNW HEADERS
      DO 249 IPTOB=1,NUMPTOB_MNW
c only proceed if header flag is false
       IF(IPTMNLST(IPTOB,2).EQ.0) THEN
         IHEADER=0
         IPTOBUN=IPTMNLST(IPTOB,1)
C  WRITE HEADER IF NOT ALREADY IN FILE
         DO 111 IPTE=1,NUMPTOB_MNW
C  CYCLE THROUGH UNIT NUMBERS, LOOK FOR MATCH IN PRECEEDING OBS
           IF(IPTE.LT.IPTOB) THEN
C  IF MATCH FOUND, SET HEADER FLAG TO TRUE
             IF(IPTOBUN.EQ.IPTMNLST(IPTE,1)) IHEADER=1
           END IF
  111    CONTINUE
         IF(IHEADER.EQ.0) THEN
          WRITE(IPTOBUN,'(256A)') '"PARTICLE OBSERVATION DATA, SITE: ',
     & PTOBLST(IPTOB),'"'
         WRITE(IPTOBUN,'(256A)') '"  J,  I,  K,         TIME          
     &PTWT        PCONC.      QSINK       QMNWSINK    
     &WT-REMV     MASS_OUT"'
c set header flag to true
         IPTMNLST(IPTOB,2)=1
         END IF
        END IF
  249 CONTINUE
      RETURN
	END
C
C
C*************************************************************************
C SPTOB5O   PRINT PARTICLE OBSERVATION DATA
C     *************************************************************
      SUBROUTINE SPTOB5O(SUMTCH,IPTOBLST,IMOV,PC,PR,PL,NPMAX,NP,TIMV,
     *        NSCOL,NSROW,NSLAY,NCOL,NROW,NLAY,NUMPTOB,NUMPTOB_MNW,
     *        IPTOBMNW,NPCELL,PTWT,PCONC,SUMWT,SNKFLO,SNKMNW,MNWFLAG)
C     ***************************************************************
      DOUBLE PRECISION SUMWT,PTWT,WTREM,MSOUT
      CHARACTER*50  LFRMAT
	REAL L
C
      DIMENSION IPTOBLST(4,NUMPTOB)
      DIMENSION PC(NPMAX),PR(NPMAX),PL(NPMAX)
      DIMENSION NPCELL(NSCOL,NSROW,NSLAY),SUMWT(NSCOL,NSROW,NSLAY)
	DIMENSION PTWT(NPMAX),PCONC(NPMAX),
     & SNKFLO(NSCOL,NSROW,NSLAY),IPTOBMNW(NSCOL,NSROW,NSLAY),
     & SNKMNW(NSCOL,NSROW,NSLAY)
C     ***************************************************************
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
CMOCWT
      INCLUDE 'ptwt.inc'
C
C  LOOP OVER NUMPTOB LOCATIONS
c does not include IPTOBMNW
      if(mnwflag.eq.0) then
      DO 10 IPTOB=1,NUMPTOB
         KO=IPTOBLST(1,IPTOB)
         KSO=KO-ISLAY1+1 
         IO=IPTOBLST(2,IPTOB)
         ISO=IO-ISROW1+1
         JO=IPTOBLST(3,IPTOB)
         JSO=JO-ISCOL1+1
         IPTOBUN=IPTOBLST(4,IPTOB)
           DO 9 IP=1,NP
      IF(PC(IP).EQ.0.0) GO TO 9
      J=PC(IP)+0.5
      JS=J-ISCOL1+1
      I=ABS(PR(IP))+0.5
      IS=I-ISROW1+1
      K=PL(IP)+0.5
      KS=K-ISLAY1+1
             IF(JSO.EQ.JS.AND.ISO.EQ.IS.AND.KSO.EQ.KS) THEN
C If a sink and using weighted pts, calculate weight removed on this particle
Cljk avoid divide by zero if no particles
cljk              IF(SNKFLO(JS,IS,KS).LT.0.0.AND.PTWTON.GT.0) THEN
              IF((SNKFLO(JS,IS,KS).LT.0.0.AND.PTWTON.GT.0).and.
     +         (sumwt(js,is,ks).gt.0)) THEN
                WTREM=-((PTWT(IP)/SUMWT(JS,IS,KS))*SNKFLO(JS,IS,KS))
cljk multiply by time
     +          *TIMV
                MSOUT=WTREM*PCONC(IP)
	        ELSE
                WTREM=0.0
                MSOUT=0.0
			END IF              
C
      IF(PTWTON.EQ.1) THEN
                WRITE(IPTOBUN,45) JO,IO,KO,SUMTCH,PTWT(IP),PCONC(IP),
     & SNKFLO(JS,IS,KS),0.0,WTREM,MSOUT
      ELSE
                WRITE(IPTOBUN,55) JO,IO,KO,SUMTCH,'     NA    ',
     & PCONC(IP),SNKFLO(JS,IS,KS),0.0,'     NA    ','     NA    '
	END IF
	        END IF
   9       CONTINUE
  10  CONTINUE
      end if
  45  FORMAT(3I4,2X,1PE18.10,6E12.4)
  55  FORMAT(3I4,2X,1PE18.10,A11,3E12.4,2A11)
C
C
C
C  DO NUMPTOB_MNW LOCATIONS WITHIN PARTICLE LOOP
      IF(NUMPTOB_MNW.GT.0) THEN
           DO 19 IP=1,NP
      IF(PC(IP).EQ.0.0) GO TO 19
      J=PC(IP)+0.5
      JS=J-ISCOL1+1
      I=ABS(PR(IP))+0.5
      IS=I-ISROW1+1
      K=PL(IP)+0.5
      KS=K-ISLAY1+1
C If printing out this well, IPTOBMNW will have a non-zero unit #
             IF(IPTOBMNW(JS,IS,KS).GT.0) THEN
              IPTOBUN=IPTOBMNW(JS,IS,KS)
C If a sink and using weighted pts, calculate weight removed on this particle
             IF(MNWFLAG.EQ.0) THEN
Cljk avoid divide by zero if no particles
cljk              IF(SNKFLO(JS,IS,KS).LT.0.0.AND.PTWTON.GT.0) THEN
              IF((SNKFLO(JS,IS,KS).LT.0.0.AND.PTWTON.GT.0).and.
     +         (sumwt(js,is,ks).gt.0)) THEN
			  WTREM=-((PTWT(IP)/SUMWT(JS,IS,KS))*SNKFLO(JS,IS,KS))
     +                *TIMV
                MSOUT=WTREM*PCONC(IP)
	        ELSE
                WTREM=0.0
                MSOUT=0.0
			END IF              
             ELSE
Cljk avoid divide by zero if no particles
cljk              IF(SNKFLO(JS,IS,KS).LT.0.0.AND.PTWTON.GT.0) THEN
              IF((SNKFLO(JS,IS,KS).LT.0.0.AND.PTWTON.GT.0).and.
     +         (sumwt(js,is,ks).gt.0)) THEN
			  WTREM=-((PTWT(IP)/SUMWT(JS,IS,KS))*SNKMNW(JS,IS,KS))
     +                *TIMV
                MSOUT=WTREM*PCONC(IP)
	        ELSE
                WTREM=0.0
                MSOUT=0.0
			END IF              
             END IF
C
             IF(MNWFLAG.EQ.0) THEN
              IF(SNKFLO(JS,IS,KS).LT.0.0) THEN
               IF(PTWTON.EQ.1) THEN
                WRITE(IPTOBUN,145) J,I,K,SUMTCH,PTWT(IP),PCONC(IP),
     & SNKFLO(JS,IS,KS),0.0,WTREM,MSOUT,' after_advection'
               ELSE
                WRITE(IPTOBUN,155) J,I,K,SUMTCH,'     NA    ',PCONC(IP),
     & SNKFLO(JS,IS,KS),0.0,'     NA    ','     NA    ',
     & ' after_advection'
               END IF
 	        END IF
Cgzh MNWFLAG is passed in >0 (=1) when called before advection for complex MNWs
             ELSE
              IF(SNKMNW(JS,IS,KS).LT.0.0) THEN
               IF(PTWTON.EQ.1) THEN
                WRITE(IPTOBUN,145) J,I,K,SUMTCH,PTWT(IP),PCONC(IP),
     & 0.0,SNKMNW(JS,IS,KS),WTREM,MSOUT,'before_advection'
               ELSE
                WRITE(IPTOBUN,155) J,I,K,SUMTCH,'     NA    ',PCONC(IP),
     & 0.0,SNKMNW(JS,IS,KS),'     NA    ','     NA    ',
     & 'before_advection'
               END IF
              END IF
	       END IF
             ENDIF
C
 145  FORMAT(3I4,2X,1PE18.10,6E12.4,X,16A)
 155  FORMAT(3I4,2X,1PE18.10,A11,3E12.4,2A11,X,16A)
   19      CONTINUE
      END IF
      RETURN
      END
