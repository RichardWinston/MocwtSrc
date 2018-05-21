c   GZH  20080208 
C
C
C GWF1MNWIAL READ INIT DATA AND ALLOCATE SPACE FOR MNW WELLS DESIGNATED FOR OBSERVATION
C
C     ******************************************************************
C
      SUBROUTINE GWF1MNWIAL(INMNWI,INMNW2,IOUT,LCMNIO,
     & Wel1flag,QSUMflag,BYNDflag,ISUMZ,MNWOBS,LSMNWO,MNWMAX)    
C
C     ******************************************************************
      IMPLICIT NONE
      INTEGER INMNWI,INMNW2,IOUT,LCMNIO,Wel1flag,QSUMflag,BYNDflag,
     & ISUMZ,MNWOBS,LSMNWO,MNWMAX
C
      IF(INMNWI.GT.0.AND.INMNW2.LE.0) THEN
        WRITE(IOUT,*) '***ERROR*** : MNWI PACKAGE CAN ONLY BE 
     *USED IF MNW2 PACKAGE IS ACTIVE'
        STOP 'MNWI ERROR'
      END IF  
C
      IF(INMNWI.EQ.0) THEN
        LCMNIO=1
      ELSE
        READ(INMNWI,*) Wel1flag,QSUMflag,BYNDflag
        WRITE(IOUT,*) 'MNWI Package input:'
        WRITE(IOUT,*) 'Wel1flag = ',Wel1flag
        WRITE(IOUT,*) 'QSUMflag = ',QSUMflag
        WRITE(IOUT,*) 'BYNDflag = ',BYNDflag
        WRITE(IOUT,*) 
C
        READ(INMNWI,*) MNWOBS
        IF(MNWOBS.LT.0) THEN
          WRITE(IOUT,*) 'MNWOBS MUST BE > 0'
          STOP
        END IF
C
        LCMNIO=ISUMZ
        ISUMZ=ISUMZ+MNWOBS*6+1
        LSMNWO=ISUMZ
        ISUMZ=ISUMZ+MNWMAX*7
      WRITE(IOUT,103) MNWOBS*6+MNWMAX*7
  103 FORMAT(1X,I8,' ELEMENTS IN RZ ARRAY ARE USED BY MNWI')
	END IF
C
      RETURN
      END
c
c_________________________________________________________________________________
c
C
C  GWF1MNWIRP READ INPUT FILE FOR MNW2 WELLS DESIGNATED FOR OBSERVATION 
C
C     ******************************************************************
C
      SUBROUTINE GWF1MNWIRP(MNWOBS,IOUT,MNWILST,INMNWI,MNWIID,MNWMAX,
     & MNW2,WELLID,GWTUNIT,MNWO, NMNWVL)
C
C     ******************************************************************
C
C     READ LOCATIONS OF MNW2 WELLS DESIGNATED FOR OBSERVATION 
C     ******************************************************************
C
c        specifications:
c     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER MNWOBS,IOUT,IOB,IS_SITE,INMNWI,iw,MNWMAX,GWTUNIT,NMNWVL
      DOUBLE PRECISION MNWILST,MNW2,MNWO
      CHARACTER*20 WELLID,SITE,MNWIID,MSITE
      DIMENSION mnw2(NMNWVL,mNWMAX),WELLID(MNWMAX),
     & MNWIID(MNWMAX),MNWILST(6,MNWOBS),MNWO(7,MNWMAX)
C
C     ******************************************************************
C
      IF(MNWOBS.EQ.0) THEN
       RETURN
	ENDIF
      IF(MNWOBS.EQ.1) THEN
       WRITE (IOUT,120) MNWOBS
	ELSEIF(MNWOBS.GT.1) THEN
       WRITE (IOUT,140) MNWOBS
	ELSEIF(MNWOBS.LT.1) THEN
       RETURN
	END IF
      WRITE (IOUT,150)
      IF(MNWOBS.GT.MNWMAX) then
        write(iout,*) '***ERROR*** MNWOBS > MNWMAX'
        STOP 'MNWI ERROR'
      end if
C
C  Initialize data array
      MNWILST=0.0
	MNWO=0.d0
C READ THE FIRST RECORD
      IOB=1
      IS_SITE=0
c
c MNWILST(1,IOB) is Well # in MNW list
c MNWILST(2,IOB) is net volume in/out well
c MNWILST(3,IOB) is unit number for output
c MNWILST(4,IOB4) is QNDflag
c MNWILST(5,IOB) is QBHflag
c MNWILST(6,IOB) is CONCflag
c
      if(GWTUNIT.GT.0) then
       READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),
     & MNWILST(5,IOB),MNWILST(6,IOB)
      else
       READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),
     & MNWILST(5,IOB)
       MNWILST(6,IOB)=0
      end if
      SITE=MNWIID(IOB)
      call UPCASE(SITE)
c check site vs list of site names in MNWSITE
c Loop over all MNW locations
c   Loop over all wells
      do iw=1,MNWMAX
        MSITE=WELLID(iw)
        call UPCASE(MSITE)
        IF(SITE.EQ.MSITE) THEN
          IS_SITE=1
          MNWILST(1,IOB)=iw
        END IF
      end do      
C
      WRITE(IOUT,'(I8,3X,A12,3I8)') IOB,SITE,INT(MNWILST(3,IOB)),
     &  INT(MNWILST(4,IOB)),INT(MNWILST(5,IOB))       
      IF(IS_SITE.EQ.0) THEN
         WRITE(IOUT,*) '***ERROR***   SITE FOR MNWI ',
     *'WELL DESIGNATED FOR OBSERVATION NOT FOUND'
         STOP 'MNWI ERROR'
      ENDIF
C CYCLE THROUGH THE REMAINING RECORDS
      DO IOB=2,MNWOBS
        IS_SITE=0
      if(GWTUNIT.GT.0) then
       READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),
     & MNWILST(5,IOB),MNWILST(6,IOB)
      else
       READ(INMNWI,*) MNWIID(IOB),MNWILST(3,IOB),MNWILST(4,IOB),
     & MNWILST(5,IOB)
       MNWILST(6,IOB)=0
      end if
c check site vs list of site names in WELLID
        SITE=MNWIID(IOB)
        call UPCASE(SITE)
c   Loop over all wells
        do iw=1,MNWMAX
          MSITE=WELLID(iw)
          call UPCASE(MSITE)
          IF(SITE.EQ.MSITE) THEN
            IS_SITE=1
            MNWILST(1,IOB)=iw
          END IF
        end do      
C
        WRITE(IOUT,'(I8,3X,A12,3I8)') IOB,SITE,INT(MNWILST(3,IOB)),
     &  INT(MNWILST(4,IOB)),INT(MNWILST(5,IOB))       
        IF(IS_SITE.EQ.0) THEN
         WRITE(IOUT,*) '***ERROR***   SITE FOR MNWI ',
     *'WELL DESIGNATED FOR OBSERVATION NOT FOUND'
          STOP 'MNWI ERROR'
        ENDIF
C
      END DO
            WRITE(IOUT,'(140A)') 'DATA FOR MNW WELLS DESIGNATED FOR
     * OBSERVATION WILL BE WRITTEN ON UNIT NUMBERS LISTED ABOVE' 
            WRITE(IOUT,'(/)')
  120 FORMAT(///'SITE ID FOR',I4,
     * ' MNW2 WELL DESIGNATED FOR OBSERVATION:')
  140 FORMAT(///'SITE IDS FOR',I4,
     * ' MNW2 WELLS DESIGNATED FOR OBSERVATION:')
  150 FORMAT(/'  WELL #   SITE ID         UNIT  QNDflag QBHflag')
      RETURN
      END       
C
c_________________________________________________________________________________
c
      subroutine GWF1MNWIOT(Wel1flag,QSUMflag,BYNDflag,nstp,kstp,
     & nmnw2,MNWMAX,MNW2,MNWNOD,NODTOT,WELLID,TOTIM,hnew,ncol,nrow,nlay,
     & MNWOBS,DELT,MNWIID,MNWILST,kper,ntotnod,iout,hdry,gwtunit,
     & gwtflag,imov,MNWO, NMNWVL)
C     VERSION 20070923 GZH
c
c     ******************************************************************
c    Sort well output into useful tables
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      IMPLICIT NONE
      ALLOCATABLE QBH(:)
      INTEGER Wel1flag,QSUMflag,BYNDflag,nstp,kstp,nmnw2,
     & MNWMAX,NODTOT,iw,INODE,firstnode,lastnode,il,ir,ic,ioutqsum,
     & ncol,nrow,nlay,MNWOBS,iwobs,QNDflag,QBHflag,CONCflag,i,istat,
     & iout,kper,ibh,ind,NNODES,numvals,ntotnod,nd,gwtunit,nnod,Inod,
     & imov,gwtflag,NMNWVL
      REAL TOTIM,DELT,hdry
      DOUBLE PRECISION mnw2,MNWNOD,q,hwell,qin,qout,qnet,hnew,hcell,
     & MNWILST,QBH,small,cnode,concw,verysmall,MNWO
      CHARACTER*20 WELLID,obssite,MNWIID,blank
      CHARACTER*150 LFRMAT
      DIMENSION mnw2(NMNWVL,mNWMAX),mnwnod(34,NODTOT),WELLID(MNWMAX),
     & hnew(ncol,nrow,nlay),MNWIID(MNWMAX),MNWILST(6,MNWOBS),
     & MNWO(MNWMAX,7) 
c
c------------------------------------------------------------------
      small=1D-5
     
      ALLOCATE(QBH(NODTOT),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,1700)ISTAT
 1700   FORMAT(1X,'ALLOCATION OF MNW ROUTING ARRAY FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
c gwtflag = 1 if called from MVOT
c
C
C MNWO array (for each MNW well)
C
C If Qnet not equal to 0 
C   1: Qnet * Cwell
C   2: Cumulative Qnet * Cwell
C
C If Qnet equal to 0 
C   1: Intra-aquifer mass-transfer
C   2: Cumulative intra-aquifer mass-transfer
C
C   3: Qsnk * Caq
C   4: Cumulative Qsnk * Caq
C   5: Qsrc * C'
C   6: Cumulative Qsrc * C'
C   7: been active before? 0=false 1=true (for header print)
C
c Print WEL1 file
      if(Wel1flag.gt.0.and.gwtflag.eq.0) then
c print max number of wells, set IWELCB=0
        if(kper.eq.1.and.kstp.eq.1)
     &  write(Wel1flag,'(2i10)') NODTOT,0
c only write at end of stress period (nstp=kstp)
        if(nstp.eq.kstp) then
c write number of wells
          write(Wel1flag,'(i10)') ntotnod 
c   Loop over all wells
          do iw=1,MNWMAX
c   Loop over nodes in well
              firstnode=MNW2(4,iw)
              lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
              do INODE=firstnode,lastnode
                il=MNWNOD(1,INODE)              
                ir=MNWNOD(2,INODE)              
                ic=MNWNOD(3,INODE)              
                q=MNWNOD(4,INODE)
c   Write L,R,C,q
                write(Wel1flag,'(i9,2i10,1x,g15.6)') il,ir,ic,q       
              end do
          end do
        end if
      end if
c  Print QSUM file
      if(QSUMflag.gt.0.and.gwtflag.eq.0) then
c   Write header
        if(kstp.eq.1.and.kper.eq.1) then
         if(gwtunit.le.0) then
           write(QSUMflag,'(200A)') 'WELLID                    Totim    
     &        Qin           Qout           Qnet          hwell'
         else
           write(QSUMflag,'(200A)') 'WELLID                    Totim    
     &        Qin           Qout           Qnet          hwell
     &       Conc'
         end if
        end if
c   Loop over all wells
        do iw=1,MNWMAX
          qin=0.0D0
          qout=0.0D0
          qnet=0.0D0
c   Only operate on active wells (MNW2(1,iw)=1)
          if (MNW2(1,iw).EQ.1) then
c   Loop over nodes in well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
            hwell=MNW2(17,iw)
            do INODE=firstnode,lastnode
              il=MNWNOD(1,INODE)              
              ir=MNWNOD(2,INODE)              
              ic=MNWNOD(3,INODE)              
              q=MNWNOD(4,INODE)
              if(q.lt.0.0D0) then
                qin=qin+q
              else
                qout=qout+q
              end if
              qnet=qnet+q
            end do
            if(qnet.lt.(small*qin)) qnet=0.d0
            if(gwtunit.le.0) then
             write(QSUMflag,'(A20,5(1x,1Pg14.6))')
     &              WELLID(iw),totim,qin,qout,qnet,hwell
            else
             concw=MNW2(31,iw)
             if(MNW2(32,iw).eq.1) then 
              write(QSUMflag,'(A20,6(1x,1Pg14.6))')
     &              WELLID(iw),totim,qin,qout,qnet,hwell,concw
             else
              write(QSUMflag,'(A20,5(1x,1Pg14.6),A11)')
     &              WELLID(iw),totim,qin,qout,qnet,hwell,' N/A'
             end if
            end if
          end if
        end do     
      end if
c   Print BYND (ByNode) file
      if(BYNDflag.gt.0.and.gwtflag.eq.0) then
c   Write header
        if(kstp.eq.1.and.kper.eq.1) then
         if(gwtunit.le.0) then
           write(BYNDflag,'(101A)') 'WELLID                NODE   Lay 
     &  Row   Col        Totim        Q-node         hwell         hcell
     &   Seepage elev.'
         else
           write(BYNDflag,'(117A)') 'WELLID                NODE   Lay 
     &  Row   Col        Totim        Q-node         hwell         hcell
     &   Seepage elev.  Concentration'
         endif
        end if

c   Loop over all wells
        do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
          if (MNW2(1,iw).EQ.1) then
c   Loop over nodes in well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
            hwell=MNW2(17,iw)
            do INODE=firstnode,lastnode
              il=MNWNOD(1,INODE)              
              ir=MNWNOD(2,INODE)              
              ic=MNWNOD(3,INODE)              
              q=MNWNOD(4,INODE)
              hcell=hnew(ic,ir,il)
              nd=INODE-firstnode+1
              cnode=MNWNOD(32,INODE)
              blank='              '
c   If no seepage face in cell, don't print seepage elev.
              if(MNWNOD(15,INODE).EQ.hwell.or.
     &           MNWNOD(15,INODE).eq.Hdry) then
               if(gwtunit.le.0) then
                write(BYNDflag,'(A20,4i6,1x,1P4e14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell
               else
                write(BYNDflag,'(A20,4i6,1x,1P4e14.6,A14,1Pe14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,blank,
     &              cnode
               end if
              else
c   If seepage face in cell, MNWNOD(15) will hold the bottom elev of the cell,
c   which is used with hcell to get the gradient used to calculate the Q for the
c   seepage face.  
               if(gwtunit.le.0) then
                write(BYNDflag,'(A20,4i6,1x,1P5e14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,
     &              MNWNOD(15,INODE)
               else
                write(BYNDflag,'(A20,4i6,1x,1P6e14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,
     &              MNWNOD(15,INODE),cnode
               end if
              end if
            end do
         end if
        end do           
      end if
c
c  Print MNWOBS files
      if((gwtunit.gt.0.and.gwtflag.gt.0).or.
     &   (gwtunit.le.0)) then
      do iwobs=1,MNWOBS    
        qnet=0.d0
        obssite=MNWIID(iwobs)
        iw=MNWILST(1,iwobs)  
c   Loop over nodes in well
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
          nnod=lastnode-firstnode+1
          hwell=MNW2(17,iw)
          qin=0.D0
          qout=0.D0
          do INODE=firstnode,lastnode
            q=MNWNOD(4,INODE)
            if(q.lt.0.0D0) then
              qin=qin+q
            else
              qout=qout+q
            end if
            qnet=qnet+q
          end do
          if (abs(qnet).lt.verysmall) qnet=0.d0
c   Cumulative volume for this well
          MNWILST(2,iwobs)=MNWILST(2,iwobs)+qnet*DELT
            
c get NNODES
          NNODES=INT(ABS(MNW2(2,iw)))
c
c  Print according to flags
          QNDflag=MNWILST(4,iwobs)
          QBHflag=MNWILST(5,iwobs)
          CONCflag=MNWILST(6,iwobs)
c
          if (gwtunit.eq.0) imov=1
c
C  Create format for header--single node
C  Ignore ND and BH flags
C  Create format for header--single node
          if(NNODES.eq.1) then
            if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
             if(gwtunit.eq.0) then
          Write (MNWILST(3,iwobs),'(A)') 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell'
             else
          Write (MNWILST(3,iwobs),'(A)') 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     &    Well-Conc'
             end if
            end if
            if(gwtunit.eq.0) then
            write(MNWILST(3,iwobs),'(A20,1x,1P6e14.6)')
     &             WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell
            else
            write(MNWILST(3,iwobs),'(A20,1x,1P7e14.6)')
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              MNW2(31,iw)
            end if
          else
c
c  all multi-node well output below
          if(QNDflag.eq.0) then
            if(QBHflag.eq.0) then
              if(gwtunit.eq.0) then            
c  QNDflag=0, QBHflag=0, gwtunit=0
c write header
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
          Write (MNWILST(3,iwobs),'(120A)') 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet      Cum.Vol.         hwell
     &   '
                end if
c   Only operate on active wells (MNW2(1,iw)=1)
                if (MNW2(1,iw).EQ.1) then

                write(MNWILST(3,iwobs),'(A20,1x,1P6e14.6)')
     &             WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell
                end if
              else
c  QNDflag=0, QBHflag=0, gwtunit.gt.0
c write header based on concflag
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
                 SELECT CASE (concflag) 
c
                 CASE (0)
          Write (MNWILST(3,iwobs),'(320A)') 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     &   Well-Conc'
c
                 CASE (1)
          Write (MNWILST(3,iwobs),'(320A)') 'WELLID                
     &     TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     &      Well-Conc    NET MASS FLUX TO/FROM WELL   MASS INTO 
     &BOREHOLE
     &         MASS OUT OF BOREHOLE'
               WRITE(MNWILST(3,iwobs),'(320A)') '                       
     &                          
     &                                                                  
     &       PRESENT 
     & t      CUM.        PRESENT t      CUM.         PRESENT t      
     &CUM. '
c
                 CASE (2)
          Write (MNWILST(3,iwobs),'(320A)') 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     &   Well-Conc    CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,114) nnod-1
 114  FORMAT('(A130,',I4,'I12)')
               WRITE(MNWILST(3,iwobs),LFRMAT) '                         
     &                              
     &                                                       NODE No.=1'
     & ,(Inod,Inod=2,nnod)
c
                 CASE (3)
          Write (MNWILST(3,iwobs),'(320A)') 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     &   Well-Conc    NET MASS FLUX TO/FROM WELL     MASS INTO BOREHOLE 
     &         MASS OUT OF BOREHOLE     
     &CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,115) nnod-1
 115  FORMAT('(A214,',I4,'I12)')
               WRITE(MNWILST(3,iwobs),LFRMAT) '                         
     &                                                                  
     &                               PRESENT 
     & t      CUM.         PRESENT t      CUM.         PRESENT t      
     &CUM.
     &      NODE No.=1'
     & ,(Inod,Inod=2,nnod)
                 END SELECT
                end if
c
                if (MNW2(1,iw).EQ.1) then
                 SELECT CASE (concflag) 
c
                 CASE (0)
                 write(MNWILST(3,iwobs),'(A20,1x,1P7e14.6)')
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &  MNW2(31,iw)
c
                 CASE (1)
                 write(MNWILST(3,iwobs),'(A20,1x,1P13e14.6)')
     &           WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &  MNW2(31,iw), (MNWO(iw,ind),ind=1,6)
c
                 CASE (2)
C Create format for write
                WRITE (LFRMAT,356) NNODES
  356           FORMAT ('(A20,1p7e14.6,1P',I4,'E12.4)')
                 write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &  MNW2(31,iw),(MNWNOD(32,inod),inod=firstnode,lastnode)
c
c
                 CASE (3)
C Create format for write
                WRITE (LFRMAT,357) NNODES
  357           FORMAT ('(A20,1p13e14.6,1P',I4,'E12.4)')
                 write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &  MNW2(31,iw),(MNWO(iw,ind),ind=1,6),
     & (MNWNOD(32,inod),inod=firstnode,lastnode)
c
                 END SELECT
c
                end if
              end if
            else
              if(gwtunit.eq.0) then            
c  QNDflag=0, QBHflag=1, gwtunit=0
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
C  Create format for header  
            WRITE(LFRMAT,14) NNODES-1
  14  FORMAT('(A,',I4,'I12)')
C  WRITE FORMAT FOR HEADER LINE
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & QBH_seg-->1',(ibh,ibh=2,NNODES)
                end if
                if (MNW2(1,iw).EQ.1) then
c
C Create format for write
                WRITE (LFRMAT,355) NNODES
  355           FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))')
             write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(27,i),i=firstnode,lastnode)
                end if
              else
c  QNDflag=0, QBHflag=1, gwtunit=1
c write header based on concflag
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,15) NNODES-1
  15  FORMAT('(A,',I4,'I12,A)')
                 SELECT CASE (concflag) 
c
                 CASE (0)
C  Create format for header  
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & QBH_seg-->1',(ibh,ibh=2,NNODES),'     Well-Conc'
c
                 CASE (1)
C  Create format for header  
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & QBH_seg-->1',(ibh,ibh=2,NNODES),'     Well-Conc      NET MASS 
     & FLUX IN WELL    MASS INTO BOREHOLE     
     &   MASS OUT OF BOREHOLE'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,111) (NNODES-1)*12
 111  FORMAT('(A137,',I4,'x,A72)')
            WRITE(MNWILST(3,iwobs),LFRMAT) ' ',
     &'PRESENT 
     & t      CUM.     PRESENT t      CUM.        PRESENT t      CUM. '
c
                 CASE (2)
C  Create format for header  
            WRITE(LFRMAT,135) NNODES-1
 135  FORMAT('(A,',I4,'I12,A)')
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & QBH_seg-->1',(ibh,ibh=2,NNODES),'   Well-Conc 
     &   CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,112) (NNODES-1)*12,NNODES-1
 112  FORMAT('(A132,',I4,'x,A10,'I4,'I12)')
               WRITE(MNWILST(3,iwobs),LFRMAT) ' ','NODE No.=1'
     & ,(Inod,Inod=2,nnod)
c
                 CASE (3)
C  Create format for header  
            WRITE(LFRMAT,615) NNODES-1
 615  FORMAT('(A,',I4,'I12,A)')
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & QBH_seg-->1',(ibh,ibh=2,NNODES),'     Well-Conc  
     & NET MASS FLUX TO/FROM WELL   MASS INTO BOREHOLE     
     &     MASS OUT OF BOREHOLE      CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,715) (NNODES-1)*12,NNODES-1
 715  FORMAT('(A136,',I4,'x,A92,',I4,'I12)')
               WRITE(MNWILST(3,iwobs),LFRMAT) ' ','PRESENT 
     & t      CUM.      PRESENT t      CUM.         PRESENT t      CUM. 
     &        NODE No.=1'
     & ,(Inod,Inod=2,nnod)
c
                 END SELECT
                end if
c                
                if (MNW2(1,iw).EQ.1) then

                 SELECT CASE (concflag) 
c
                 CASE (0)
C Create format for write
                WRITE (LFRMAT,255) NNODES
  255           FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4),1pe14.6)')
                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw)
                 CASE (1)
C Create format for write
                WRITE (LFRMAT,256) NNODES
  256           FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4),1p7e14.6)')
                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWO(iw,ind),ind=1,6)

                 CASE (2)
C Create format for write
                WRITE (LFRMAT,257) NNODES,NNODES
  257          FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4),1pe14.6,
     &',I4,'(1pE12.4))')
               write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWNOD(32,inod),inod=firstnode,lastnode)
                 CASE (3)
C Create format for write
                WRITE (LFRMAT,258) NNODES,NNODES
  258          FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4),1p8e14.6,
     &',I4,'(1pE12.4))')
               write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWO(iw,ind),ind=1,6),
     &              (MNWNOD(32,inod),inod=firstnode,lastnode)
                 END SELECT
                end if
              end if
            end if  
          else
c
            if(QBHflag.eq.0) then
              if(gwtunit.eq.0) then            
c  QNDflag=1, QBHflag=0, gwtunit=0
c write "smart" header
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
C  Create format for header  
            WRITE(LFRMAT,14) NNODES-1
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES)
                end if
                if (MNW2(1,iw).EQ.1) then
C Create format for write
                WRITE (LFRMAT,155) NNODES
  155           FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))')

                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode)
                end if
              else
c  QNDflag=1, QBHflag=0, gwtunit=1
c
c write "smart" header
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
                 SELECT CASE (concflag) 
c
                 CASE (0)
C  Create format for header  
            WRITE(LFRMAT,15) NNODES-1
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),'   Well-Conc'
                 CASE (1)
C  Create format for header  
            WRITE(LFRMAT,15) NNODES-1
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),'     Well-Conc     NET MASS 
     & FLUX IN WELL     MASS INTO BOREHOLE     
     &    MASS OUT OF BOREHOLE'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,111) (NNODES-1)*12
            WRITE(MNWILST(3,iwobs),LFRMAT) ' ',
     &'PRESENT 
     & t      CUM.     PRESENT t      CUM.         PRESENT t      CUM. '
                 CASE (2)
C  Create format for header  
            WRITE(LFRMAT,15) NNODES-1
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),'   Well-Conc 
     &   CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,112) (NNODES-1)*12,NNODES-1
               WRITE(MNWILST(3,iwobs),LFRMAT) ' ','NODE No.=1'
     & ,(Inod,Inod=2,nnod)
                 CASE (3)
C  Create format for header  
            WRITE(LFRMAT,15) NNODES-1
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),'     Well-Conc  
     & NET MASS FLUX TO/FROM WELL   MASS INTO BOREHOLE     
     &     MASS OUT OF BOREHOLE      CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,715) (NNODES-1)*12,NNODES-1
               WRITE(MNWILST(3,iwobs),LFRMAT) ' ','PRESENT 
     & t      CUM.      PRESENT t      CUM.         PRESENT t      CUM. 
     &        NODE No.=1'
     & ,(Inod,Inod=2,nnod)
                 END SELECT
                end if
c
                if (MNW2(1,iw).EQ.1) then
                 SELECT CASE (concflag) 
c
                 CASE (0)
C Create format for write
                WRITE (LFRMAT,255) NNODES
                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),MNW2(31,iw)
c
                 CASE (1)
C Create format for write
                WRITE (LFRMAT,256) NNODES
                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWO(iw,ind),ind=1,6)
c
                 CASE (2)
C Create format for write
                WRITE (LFRMAT,257) NNODES,NNODES
                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWNOD(32,inod),inod=firstnode,lastnode)
c
                 CASE (3)
C Create format for write
                WRITE (LFRMAT,258) NNODES,NNODES
                write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWO(iw,ind),ind=1,6),
     &              (MNWNOD(32,inod),inod=firstnode,lastnode)
c
                 END SELECT
                end if
              end if
            else
              if(gwtunit.eq.0) then            
c  QNDflag=1, QBHflag=1, gwtunit=0
C Create format for write

c write "smart" header
                if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
C  Create format for header  
          WRITE(LFRMAT,16) NNODES-1,NNODES-1
  16  FORMAT('(A,',I4,'I12,A,',I4,'I12)')
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),
     &' QBH_seg-->1',(ibh,ibh=2,NNODES)
                end if
                if (MNW2(1,iw).EQ.1) then
       numvals=NNODES+NNODES
               WRITE (LFRMAT,156) numvals
  156          FORMAT ('(A20,1p6e14.6,',I4,'(1pE12.4))')
                  write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),
     &              (MNWNOD(27,i),i=firstnode,lastnode)

                end if
              else
c  QNDflag=1, QBHflag=1, gwtunit=1
       if(kstp.eq.1.and.kper.eq.1.and.imov.eq.1) then
                 SELECT CASE (concflag) 
c
                 CASE (0)
          WRITE(LFRMAT,17) NNODES-1,NNODES-1
  17  FORMAT('(A,',I4,'I12,A,',I4,'I12,A)')
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),
     &' QBH_seg-->1',(ibh,ibh=2,NNODES),'     Well-Conc'
c
                 CASE (1)
          WRITE(LFRMAT,117) NNODES-1,NNODES-1
 117  FORMAT('(A,',I4,'I12,A,',I4,'I12,A)')
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),
     &' QBH_seg-->1',(ibh,ibh=2,NNODES),'     Well-Conc    NET MASS 
     & FLUX IN WELL     MASS INTO BOREHOLE     
     &      MASS OUT OF BOREHOLE'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,417) (NNODES-1)*12,(NNODES-1)*12
 417  FORMAT('(A136,',I4,'x,A12,',I4,'x,A)')
            WRITE(MNWILST(3,iwobs),LFRMAT) ' ',' ',
     &'PRESENT 
     & t      CUM.     PRESENT t      CUM.          PRESENT t      CUM.'
c
                 CASE (2)
          WRITE(LFRMAT,117) NNODES-1,NNODES-1
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),
     &' QBH_seg-->1',(ibh,ibh=2,NNODES),'   Well-Conc 
     &   CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,517) (NNODES-1)*12,(NNODES-1)*12,NNODES-1
 517  FORMAT('(A132,',I4,'x,A12,',I4,'x,A,',I4,'I12)')
               WRITE(MNWILST(3,iwobs),LFRMAT) ' ',' ','NODE No.=1'
     & ,(Inod,Inod=2,nnod)
c
                 CASE (3)
          WRITE(LFRMAT,317) NNODES-1,NNODES-1
 317  FORMAT('(A,',I4,'I12,A,',I4,'I12,A)')
          Write (MNWILST(3,iwobs),LFRMAT) 'WELLID                
     &       TOTIM           Qin
     &          Qout          Qnet         QCumu         hwell
     & Flow@Nd-->1',(ind,ind=2,NNODES),
     &' QBH_seg-->1',(ibh,ibh=2,NNODES),'     Well-Conc  
     & NET MASS FLUX TO/FROM WELL   MASS INTO BOREHOLE     
     &     MASS OUT OF BOREHOLE      CONCENTRATION AT NODES IN WELL...'
C  WRITE FORMAT FOR HEADER LINE
            WRITE(LFRMAT,815) (NNODES-1)*12,(NNODES-1)*12,NNODES-1
 815  FORMAT('(A135,',I4,'x,A12,',I4,'x,A,',I4,'I12)')
               WRITE(MNWILST(3,iwobs),LFRMAT) ' ',' ','PRESENT 
     & t      CUM.       PRESENT t      CUM.         PRESENT t      CUM.
     &       NODE No.=1'
     & ,(Inod,Inod=2,nnod)
     
                 END SELECT
        end if
                if (MNW2(1,iw).EQ.1) then
                 SELECT CASE (concflag) 
c
                 CASE (0)
                WRITE (LFRMAT,173) NNODES, 
     & NNODES
  173         FORMAT 
     & ('(A20,1p6e14.6,',I4,'(1pE12.4),',I4,'(1pE12.4),1pe14.6)')
                  write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw)
c
                 CASE (1)
                WRITE (LFRMAT,273) NNODES, 
     & NNODES
  273         FORMAT 
     & ('(A20,1p6e14.6,',I4,'(1pE12.4),',I4,'(1pE12.4),1p7e14.6)')
                  write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWO(iw,ind),ind=1,6)
c
                 CASE (2)
                WRITE (LFRMAT,373) NNODES,NNODES,NNODES
  373         FORMAT 
     & ('(A20,1p6e14.6,',I4,'(1pE12.4),',I4,'(1pE12.4),1pe14.6,
     &',I4,'(1pE12.4))')
                  write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWNOD(32,inod),inod=firstnode,lastnode)
c
                 CASE (3)
                WRITE (LFRMAT,473) NNODES,NNODES,NNODES
  473         FORMAT 
     & ('(A20,1p6e14.6,',I4,'(1pE12.4),',I4,'(1pE12.4),1p7e14.6,
     &',I4,'(1pE12.4))')
                  write(MNWILST(3,iwobs),LFRMAT)
     &            WELLID(iw),totim,qin,qout,qnet,MNWILST(2,IWOBS),hwell,
     &              (MNWNOD(4,i),i=firstnode,lastnode),
     &              (MNWNOD(27,i),i=firstnode,lastnode),MNW2(31,iw),
     &              (MNWO(iw,ind),ind=1,6),
     &              (MNWNOD(32,inod),inod=firstnode,lastnode)
c
                  END SELECT
                end if
              end if
            end if  
          end if  
          end if   
      end do
      end if
c
      return
      end
