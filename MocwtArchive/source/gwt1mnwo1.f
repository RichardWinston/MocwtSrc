c
c
C  MNW WELLS DESIGNATED FOR OBSERVATION 
C  1/05
C
C
C GWT1MNWO5AL READ INIT DATA AND ALLOCATE SPACE FOR MNW WELLS DESIGNATED FOR OBSERVATION
C
C     ******************************************************************
C
      SUBROUTINE GWT1MNWO5AL(ISUM,ISUMI,LSMNWU,LSMNWO,LSIQZE,
     *                  MNWOBS,mxwel2,
     *                  IOUT,IOUTS,INMNWO,INMNW)    
C
C     ******************************************************************
C
      IF(INMNWO.GT.0.AND.INMNW.LE.0) THEN
        WRITE(IOUTS,*) '***ERROR*** : MNWO PACKAGE CAN ONLY BE 
     *USED IF MNW PACKAGE IS ACTIVE'
        STOP 'MNWO ERROR'
      END IF  
C
      IF(INMNWO.EQ.0) THEN
        LSMNWU=1
        LSMNWO=1
        LSIQZE=1
      ELSE
        READ(INMNWO,*) MNWOBS
        IF(MNWOBS.LE.0) THEN
          WRITE(IOUTS,*) 'MNWOBS MUST BE > 0'
          STOP
        END IF
        LSMNWU=ISUMI
        ISUMI=ISUMI+MNWOBS*2
        LSIQZE=ISUMI
        ISUMI=ISUMI+mxwel2
        LSMNWO=ISUM
        ISUM=ISUM+mxwel2*7
        ISPI=MNWOBS*2+mxwel2
        ISPX=mxwel2*6
      WRITE(IOUTS,102) ISPI
  102 FORMAT(1X,I8,' ELEMENTS IN IR ARRAY ARE USED BY MNWO')
      WRITE(IOUTS,103) ISPX
  103 FORMAT(1X,I8,' ELEMENTS IN RX ARRAY ARE USED BY MNWO')
	END IF
C
      RETURN
      END
C
C  MNWO5RP READ INPUT FILE FOR MNW WELLS DESIGNATED FOR OBSERVATION 
C
C     ******************************************************************
C
      SUBROUTINE GWT1MNWO5RP(MNWOLST,MNWUNIT,MNWOBS,
     *                  IOUTS,INMNWO,
     *                  well2,MNWSITE,mxwel2,nwell2,MNWO)
C
C     ******************************************************************
C
C     READ LOCATIONS OF MNW WELLS DESIGNATED FOR OBSERVATION 
C     ******************************************************************
C
      real MNWO
      double precision well2
      DIMENSION MNWOLST(MNWOBS)               
      DIMENSION MNWUNIT(MNWOBS,2)               
      DIMENSION MNWO(mxwel2,7)               
      dimension well2(18,mxwel2)
      dimension MNWsite(mxwel2)
      character*32 MNWsite,MNWOLST
      CHARACTER*32 SITE,MSITE
C
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
C
C     ******************************************************************
C
      IF(MNWOBS.EQ.1) THEN
       WRITE (IOUTS,120) MNWOBS
	ELSEIF(MNWOBS.GT.1) THEN
       WRITE (IOUTS,140) MNWOBS
	ELSEIF(MNWOBS.LT.1) THEN
       RETURN
	END IF
      WRITE (IOUTS,150)
C
C  Initialize data array
      MNWO=0.0
C READ THE FIRST RECORD
      IOB=1
      IS_SITE=0
c MNWUNIT(IOB,2) is MNWOflag
        READ(INMNWO,*) MNWOLST(IOB),MNWUNIT(IOB,1),MNWUNIT(IOB,2)
        SITE=MNWOLST(IOB)
        call UPCASE(SITE)
c check site vs list of site names in MNWSITE
c Loop over all MNW locations
         m = 0
         do while( m .lt. nwell2 )
          m = m + 1
          MSITE=MNWSITE(m)
          call UPCASE(MSITE)
          IF(SITE.EQ.MSITE) IS_SITE=1
         end do      
C
         WRITE(IOUTS,'(I8,3X,A12,2I8)') IOB,SITE,MNWUNIT(IOB,1),
     &  MNWUNIT(IOB,2)        
c         IF(IS_SITE.EQ.0) THEN
c            WRITE(IOUTS,*) '***ERROR***   SITE FOR MNW ',
c     *'WELL DESIGNATED FOR OBSERVATION NOT FOUND'
cgzh debug            STOP 'MNWO ERROR'
c         ENDIF
C CYCLE THROUGH THE REMAINING RECORDS
      DO 139 IOB=2,MNWOBS
        IS_SITE=0
         READ(INMNWO,*) MNWOLST(IOB),MNWUNIT(IOB,1),MNWUNIT(IOB,2)
         SITE=MNWOLST(IOB)
c check site vs list of site names in MNWSITE
c Loop over all MNW locations
         call UPCASE(SITE)
         m = 0
         do while( m .lt. nwell2 )
          m = m + 1
          MSITE=MNWSITE(m)
          call UPCASE(MSITE)
          IF(SITE.EQ.MSITE) IS_SITE=1
         end do      
C
         WRITE(IOUTS,'(I8,3X,A12,I8,I10)') IOB,SITE,MNWUNIT(IOB,1),
     * MNWUNIT(IOB,2)        
c         IF(IS_SITE.EQ.0) THEN
c            WRITE(IOUTS,*) '***ERROR***   SITE FOR MNW ',
c     *'WELL DESIGNATED FOR OBSERVATION NOT FOUND'
cgzh debug            STOP
c         ENDIF
C
  139 CONTINUE
            WRITE(IOUTS,'(/)')
            WRITE(IOUTS,'(140A)') 'DATA FOR MNW WELLS DESIGNATED FOR
     * OBSERVATION WILL BE WRITTEN ON UNIT NUMBERS LISTED ABOVE' 
            WRITE(IOUTS,'(/)')
            WRITE(IOUTS,*) 'MNWO OUTPUT:'
      WRITE(iouts,*) 'MNWOflag = 0: TIME, CONCENTRATION AT WELL'
      WRITE(iouts,*) 'MNWOflag = 1: TIME, CONCENTRATION AT WELL,
     * PLUS EXTENDED MNW MASS FLUX INFO'
      WRITE(iouts,*) 'MNWOflag = 2: TIME, CONCENTRATION AT WELL,
     * PLUS CONCENTRATIONS AT WELL NODES'
      WRITE(iouts,*) 'MNWOflag = 3: TIME, CONCENTRATION AT WELL,
     * MASS FLUX INFO, CONC. AT WELL NODES'


  120 FORMAT(///'SITE ID FOR',I4,
     * ' MNW WELL DESIGNATED FOR OBSERVATION:')
  140 FORMAT(///'SITE IDS FOR',I4,
     * ' MNW WELLS DESIGNATED FOR OBSERVATION:')
  150 FORMAT(/'  WELL #   SITE ID         UNIT  MNWOflag')
      RETURN
      END       
C
C
C*************************************************************************
C SMNWO5O   PRINT DATA FOR MNW WELLS DESIGNATED FOR OBSERVATION
C     *************************************************************
      SUBROUTINE SMNWO5O(SUMTCH,MNWOLST,MNWUNIT,MNWOBS,
     *                IMOV,well2,MNWSITE,mxwel2,nwell2,CWELL,MNWO,
     *                IQzero,KSTP,KPER)
C     ***************************************************************
      CHARACTER*50  LFRMAT
C
C MNWO array (for each MNW well)
C
C If Qnet not equal to 0 as defined by user in 1st SP
C   1: Qnet * Cwell
C   2: Cumulative Qnet * Cwell
C If Qnet equal to 0 as defined by user in 1st SP
C   1: Intra-aquifer mass-transfer
C   2: Cumulative intra-aquifer mass-transfer
C
C   3: Qsnk * Caq
C   4: Cumulative Qsnk * Caq
C   5: Qsrc * C'
C   6: Cumulative Qsrc * C'
C   7: been active before? 0=false 1=true (for header print)
      real MNWO
      double precision well2
      DIMENSION MNWOLST(MNWOBS)               
      DIMENSION MNWUNIT(MNWOBS,2)               
      dimension MNWsite(mxwel2)
      dimension IQzero(mxwel2)
      dimension CWELL(mxwel2)
      dimension MNWO(mxwel2,7)
      dimension well2(18,mxwel2)
      character*32 MNWsite,MNWOLST
      CHARACTER*32 SITE,MSITE
C     ***************************************************************
C
      ALLOCATABLE CONCMNW(:,:),mmap(:,:),mnwsingle(:)
C
      COMMON /SUBGRD/ 
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD 
cljk change mmap to dim (x,3) 
cljk        ALLOCATE (CONCMNW(mxwel2,7),mmap(mxwel2,2),mnwsingle(mxwel2)) 
      ALLOCATE (CONCMNW(mxwel2,7),mmap(mxwel2,3),mnwsingle(mxwel2)) 
      mnwsingle=0 
C FIRST LOOP OVER MNW OBS SITES: 
C build CONCMNW array with CWELL and MNWO 
      do 100 iob=1,MNWOBS 
c initialize 
      mmap(iob,1)=0 
      mmap(iob,2)=0 
      mmap(iob,3)=0
c Loop over all MNW locations
      m = 0 
        do while( m .lt. nwell2 )
         m = m + 1
c   A very large # in WL reference array (8,m) flags multi-node wells
c
        if( well2(8,m) .gt. 1.0E30 ) then
c   Set ne = number of last node in this MNW
          ne  = ifrl(well2(7,m))
c   If this is the site for this observation (iob), save conc
          SITE=MNWSITE(m)
          MSITE=MNWOLST(iob)
          MNWOflag=MNWUNIT(IOB,2)
          call UPCASE(SITE)
          call UPCASE(MSITE)
          if(SITE.EQ.MSITE) then
		  CONCMNW(iob,1)=well2(12,m)
c   print additional data if requested
            if(MNWOflag.eq.1.or.MNWOflag.eq.3) then
		    CONCMNW(iob,2)=MNWO(m,1)
		    CONCMNW(iob,3)=MNWO(m,2)
		    CONCMNW(iob,4)=MNWO(m,3)
		    CONCMNW(iob,5)=MNWO(m,4)
		    CONCMNW(iob,6)=MNWO(m,5)
		    CONCMNW(iob,7)=MNWO(m,6)
            end if
            mmap(iob,1)=m 
            mmap(iob,2)=ne 
cljk add flag for printing 
            mmap(iob,3)=1 
c   found a match--so skip out to next iob
            go to 100
          end if
c   set counter to last node of MNW
          m = ne  
        else
C   Single-node wells
c   If this is the site for this observation (iob), save conc
          SITE=MNWSITE(m)
          MSITE=MNWOLST(iob)
          call UPCASE(SITE)
          call UPCASE(MSITE)
          if(SITE.EQ.MSITE) then
		  CONCMNW(iob,1)=well2(12,m)
c   print additional data if requested
cgzh debug turned this off as backing out this info not simple
c            if(MNWOflag.eq.1) then
c		    CONCMNW(iob,2)=MNWO(m,1)
c		    CONCMNW(iob,3)=MNWO(m,2)
c		    CONCMNW(iob,4)=MNWO(m,3)
c		    CONCMNW(iob,5)=MNWO(m,4)
c		    CONCMNW(iob,6)=MNWO(m,5)
c		    CONCMNW(iob,7)=MNWO(m,6)
c           end if
            mmap(iob,1)=m 
            mmap(iob,2)=m 
cljk add flag for printing 
            mmap(iob,3)=1 
            mnwsingle(m)=1
c   found a match--so skip out
            go to 100
          end if
        end if
c       end if multi-node
        enddo
c   End Loop over all wells      
 100  end do
c   End loop over MNW obs sites
c
c
c  additional columns for Qnet ne 0 wells: QnetCwell, sumQnetCwell
c  additional columns for Qnet eq 0 wells: intra-aq mass, sum intra-aq mass
c  additional columns for all wells QsnkCaq, sumQsnkCaq, QsrcCwell, sumQsrcCwell
c
c
C  SECOND LOOP OVER printable MNW OBS SITES 
      DO 10 IOB=1,MNWOBS 
c check print flag
        if (mmap(iob,3).eq.1) then 
         IOBUN=MNWUNIT(IOB,1) 
         MNWOflag=MNWUNIT(IOB,2) 
         miob=mmap(iob,1) 
         neiob=mmap(iob,2) 
         numob=neiob-miob+1
C  WRITE HEADER
cgzh         IF(IMOV.EQ.1.AND.KSTP.EQ.1.AND.KPER.EQ.1) THEN
c print header if header flag is 0
         IF(MNWO(miob,7).EQ.0) THEN
           WRITE(IOBUN,*) 'DATA FOR MNW WELL DESIGNATED FOR OBSERVATION'
               WRITE(IOBUN,*) 'SITE  ',MNWOLST(IOB)
           IF(MNWOflag.eq.0.or.mnwsingle(m).eq.1) then
               WRITE(IOBUN,*) '   TIME     WELL-CONC'
           ELSEIF(MNWOflag.eq.1) then
             if(IQzero(miob).eq.0) then
               WRITE(IOBUN,'(120A)') '   TIME      WELL-CONC   NET MASS
     & FLUX IN WELL    MASS INTO BOREHOLE      MASS OUT OF BOREHOLE'
               WRITE(IOBUN,'(120A)') '                         PRESENT 
     & t      CUM.     PRESENT t      CUM.     PRESENT t      CUM. '
             else
               WRITE(IOBUN,'(120A)') '   TIME      WELL-CONC   MASS FLUX
     & IN BOREHOLE     MASS INTO BOREHOLE      MASS OUT OF BOREHOLE'
               WRITE(IOBUN,'(120A)') '                         PRESENT 
     & t      CUM.     PRESENT t      CUM.     PRESENT t      CUM. '
             end if
           ELSEIF(MNWOflag.eq.2) then
             WRITE(IOBUN,'(120A)') '   TIME      WELL-CONC',
     & '   CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,14) numob-1
  14  FORMAT('(A35,',I4,'I12)')
               WRITE(IOBUN,LFRMAT) '                         NODE NUM=1'
     & ,(IOB2,IOB2=2,numob)
           ELSEIF(MNWOflag.eq.3) then
             if(IQzero(miob).eq.0) then
               WRITE(IOBUN,'(120A)') '   TIME      WELL-CONC   NET MASS
     & FLUX IN WELL    MASS INTO BOREHOLE      MASS OUT OF BOREHOLE',
     &'   CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,15) numob
  15  FORMAT('(A107,',I4,'I12)')
               WRITE(IOBUN,LFRMAT) '                         PRESENT 
     & t      CUM.     PRESENT t      CUM.     PRESENT t      CUM.    
     &NODE NUM=1',(IOB2,IOB2=2,numob)
             else
               WRITE(IOBUN,'(120A)') '   TIME      WELL-CONC   MASS FLUX
     & IN BOREHOLE     MASS INTO BOREHOLE      MASS OUT OF BOREHOLE',
     &'   CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,15) numob
               WRITE(IOBUN,LFRMAT) '                         PRESENT 
cgzh add NODE #s here
     & t      CUM.     PRESENT t      CUM.     PRESENT t      CUM.    
     &NODE NUM=1',(IOB2,IOB2=2,numob)
             end if
           END IF
c set header print flag to true
           MNWO(miob,7)=1
         ENDIF
      TIME=ABS(SUMTCH)
          IF(MNWOflag.eq.0.or.mnwsingle(m).eq.1) THEN
            WRITE(IOBUN,45) TIME,CONCMNW(IOB,1)
          ELSEIF(MNWOflag.eq.1) THEN
            WRITE(IOBUN,145) TIME,(CONCMNW(IOB,I),I=1,7)
  145  FORMAT(1PE11.4,8E12.4)
         ELSEIF(MNWOflag.eq.2) THEN
C Create format for MNWOflag=2 write
            WRITE (LFRMAT,155) numob
  155  FORMAT ('(1PE11.4,E12.4,',I3,'(E12.4))')
            WRITE(IOBUN,LFRMAT) TIME,CONCMNW(IOB,1),
     &        (Cwell(iin),iin=miob,neiob)
         ELSEIF(MNWOflag.eq.3) THEN
C Create format for MNWOflag=3 write
            WRITE (LFRMAT,165) numob
  165  FORMAT ('(1PE11.4,7E12.4,',I3,'(E12.4))')
           WRITE(IOBUN,LFRMAT) TIME,(CONCMNW(IOB,I),I=1,7),
     &        (Cwell(iin),iin=miob,neiob)
          ENDIF
        end if 
  10  CONTINUE
  45  FORMAT(1PE11.4,2E12.4)
  50  FORMAT(3X,'H & C AT ',3I3,1X)
      DEALLOCATE(CONCMNW,mmap,mnwsingle)
      RETURN
      END
