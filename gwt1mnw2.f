C
C
C GWT1MNW2AL ALLOCATE SPACE FOR MNWid
C
C     ******************************************************************
C
      SUBROUTINE GWT1MNW2AL(ISUMI,LSMNWI,mnwmax,IOUTS)    
C
C     ******************************************************************
C
C     ARRAY IS LAY,ROW,COL,IFACE
      IF(mnwmax.GT.0) THEN
        LSMNWI=ISUMI
        ISUMI=ISUMI+mnwmax
        WRITE(IOUTS,101) mnwmax
      ELSE
	  LSMNWI=1
      END IF
  101 FORMAT(1X,I8,' ELEMENTS IN IR ARRAY ARE USED BY GWT-MNW2')
      RETURN
      END
C
c
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW2bd(nmnw2,MNWMAX,mnw2,NODTOT,MNWNOD,ibound,
     +             ncol,nrow,nlay,WELLID,MNWPRNT,
     +             iouts,conc,srcflo,srcsol,snkflo,nscol,nsrow,nslay,
     +             SRCMNW,SNKMNW,SOLMNW,KSTP,
     +             kkper,IPERGWT,NMNWVL)
C
c     ******************************************************************
c     calculate solute mass terms for mnw wells
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      CHARACTER*20 WELLID
      INTEGER IBOUND,NMNW2,IW,MNWMAX,firstnode,lastnode,INODE,il,ir,ic,
     & NODTOT,ncol,nrow,nlay,kstp,nodes,
     & kper,MNWPRNT
      DOUBLE PRECISION mnw2,qdes,
     & small,hwell,MNWNOD,QCUT,qnet
      dimension mnw2(NMNWVL,MNWMAX),mnwnod(34,NODTOT)
      dimension ibound(NCOL,NROW,NLAY)
      dimension WELLID(mnwmax)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY)
C
      DIMENSION SRCMNW(NSCOL,NSROW,NSLAY),SNKMNW(NSCOL,NSROW,NSLAY),
     *  SOLMNW(NSCOL,NSROW,NSLAY)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      character*16 text
      text = '            MNW2'
c     ------------------------------------------------------------------
c
C Initialize SRCMNW,SNKMNW,AND SOLMNW
      SRCMNW=0.
      SNKMNW=0.
      SOLMNW=0.
c-----if there are no wells do not process
      if(nmnw2.gt.0) then
c Initialize counter used to label MNWs
        id=0
c Initialize flag for printing conc only for 1st node of MNW
        idcheck=0
c Initialize print flag for header line to main file
        IPRNT=0
c
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
c        process multi-node wells
         if(ABS(MNW2(2,iw)).gt.0) then
c   Loop 1: check for MNW's that go across subgrid boundary
c   Flag active wells for transport
c
c   Loop over nodes in well
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
          do INODE=firstnode,lastnode
            k=MNWNOD(1,INODE)              
            i=MNWNOD(2,INODE)              
            j=MNWNOD(3,INODE)              
c
            JS=J-ISCOL1+1
            IS=I-ISROW1+1
            KS=K-ISLAY1+1
c set active transport well flag
           MNW2(32,iw)=0
           IF(KS.GE.1.AND.KS.LE.NSLAY
     &       .AND.IS.GE.1.AND.IS.LE.NSROW
     &       .AND.JS.GE.1.AND.JS.LE.NSCOL) then 
               MNW2(32,iw)=1
           END IF            

            if (INODE.eq.firstnode) then
c   Record where first node is (in or out)
              IF(KS.LT.1.OR.KS.GT.NSLAY
     &       .OR.IS.LT.1.OR.IS.GT.NSROW
     &       .OR.JS.LT.1.OR.JS.GT.NSCOL) then 
                Ifirst=1
              ELSE
                Ifirst=0
              END IF
            else
c   Check other nodes vs the first
              IF(KS.LT.1.OR.KS.GT.NSLAY
     &       .OR.IS.LT.1.OR.IS.GT.NSROW
     &       .OR.JS.LT.1.OR.JS.GT.NSCOL) then 
                Icurr=1
              ELSE
                Icurr=0
              END IF
              
c   If some nodes in and some nodes out of subgrid, stop
              if(Icurr.ne.Ifirst) then
                write(iouts,*)                 
                write(iouts,*)                 
                write(iouts,*) '***ERROR*** Multi-node well has nodes
     & inside and outside the GWT transport subgrid'    
                write(iouts,*) 'Error occurred on well at col, row,
     & lay=', i,j,k    
                STOP 'ERROR IN GWT (MNW) '          
              end if
            end if
c            
          enddo
c   End Loop 1
c
C   print header for output
      IF(KSTP.EQ.1.AND.IPRNT.EQ.0) WRITE(IOUTS,129)
      IPRNT=1
  129  FORMAT(/1H ,'  MNW CONCENTRATIONS FOR THIS STRESS PERIOD'
     1 /'  (ONLY FOR MNW NODES WITHIN TRANSPORT SUBGRID)'//
     2 1H ,47X,'EXTERNAL',/
     3 1H ,1X,'WELL NO.  LAYER    ROW   COLUMN    VOL/T   ',
     4 '  SOURCE CONC   SITEID'/1X,67('-'))
c
c   Loop 3: process each source node in simple MNW
c   Loop over nodes in this MNW
          do INODE=firstnode,lastnode
c get node location 
            k=MNWNOD(1,INODE)              
            i=MNWNOD(2,INODE)              
            j=MNWNOD(3,INODE)              
c
            JS=J-ISCOL1+1
            IS=I-ISROW1+1
            KS=K-ISLAY1+1
            IF(MNW2(32,iw).eq.1) then 
            if( ibound(j,i,k) .ne. 0 ) then
             qnode = MNWNOD(4,INODE)
			  IF(qnode.LT.0.0) THEN
                 SNKMNW(JS,IS,KS)=SNKMNW(JS,IS,KS)+qnode
! RBW Is SNKFLO is needed by MNWI
!                 IF(PTWTON.EQ.0) SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+qnode
                 SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+qnode
C   print output (only print Cinput for first well in list)
                 IF(KSTP.EQ.1) then
                   WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 iw,K,I,J,qnode,MNW2(12,iw),WELLID(iw)
                 END IF
                ELSE IF(qnode.GT.0.0) THEN
                 SRCMNW(JS,IS,KS)=SRCMNW(JS,IS,KS)+qnode
! RBW Is SRCFLO is needed by MNWI
!                 IF(PTWTON.EQ.0) SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+qnode
                 SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+qnode
C   print output (only print Cinput for first well in list)
                 IF(KSTP.EQ.1) then
                   WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 iw,K,I,J,qnode,MNW2(12,iw),WELLID(iw)
                 end if
c   SOLMNW (based on Cwell, not Cinput) is updated each move in GWT1MNW2cx
                END IF
              end if
            end if
c   End Loop 3
          enddo
C
        else       !  End of multi-node conditioning IF statement
c     Single-node well
         MNW2(32,iw)=1
         INODE=MNW2(4,iw)
         k=MNWNOD(1,INODE)              
         i=MNWNOD(2,INODE)              
         j=MNWNOD(3,INODE)              
c
         JS=J-ISCOL1+1
         IS=I-ISROW1+1
         KS=K-ISLAY1+1
         IF(KS.LT.1.OR.KS.GT.NSLAY) MNW2(32,iw)=0
         IF(IS.LT.1.OR.IS.GT.NSROW) MNW2(32,iw)=0
         IF(JS.LT.1.OR.JS.GT.NSCOL) MNW2(32,iw)=0
         IF(MNW2(32,iw).eq.0) GO TO 20
         if( ibound(j,i,k) .eq. 0 ) GO TO 20
C
          qnode = MNWNOD(4,INODE)
          IF(qnode.LT.0.0) THEN
             SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+qnode
C   print output
                 IF(KSTP.EQ.1) 
     *             WRITE(IOUTS,'(1X,I6,3I8,2X,1P1E12.4,16X,A)') 
     1		                 iw,K,I,J,qnode,WELLID(iw)
          ELSE IF(qnode.GT.0.0) THEN
             SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+qnode
             qnodeC=qnode*MNW2(12,iw)
             if(MNW2(12,iw).lt.0.0) then
               write(iouts,*) '***ERROR*** MNW concentration < 0 at
     & col,row,lay=', i,j,k
               write(iouts,*) 'Execution stopping.'
               STOP 'ERROR in GWT (MNW)'
             end if
              SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+qnodeC
C   print output
                 IF(KSTP.EQ.1) 
     *             WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 iw,K,I,J,qnode,MNW2(12,iw),WELLID(iw)
          END IF
  20  CONTINUE
        endif
c
        end if
c       end if active well=true       
      enddo
c     end loop over all wells
      endif
c     end if nmnw2>0 
c-----return
      return
      end
c
c_______
c
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW2cx(nmnw2,MNWMAX,mnw2,NODTOT,MNWNOD,ibound,
     +          ncol,nrow,nlay,WELLID,MNWPRNT,
     +          iouts,conc,srcflo,srcsol,snkflo,nscol,nsrow,nslay,
     +          MNWid,SRCMNW,SNKMNW,SOLMNW,SRCVOL,
     *          RF,DECAY,TIMV,SRCDCY,
     *          DKFO,DKFS,INDK,IDKFO,IDKFS,
     *          SBVL,NIUNIT,Cwell,THCK,imov,
     *          BOTM,NBOTM,MNWIin,MNWO, NMNWVL)
C
c     ******************************************************************
c     calculate solute mass terms for all mnw wells
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      CHARACTER*20 WELLID
      DOUBLE PRECISION mnw2,qdes,
     & small,hwell,MNWNOD,QCUT,qnet
      double precision QSnk,Msnk,qnode,qwt
      dimension WELLID(MNWMAX)
      dimension Cwell(NODTOT)
      REAL MNWD
      DOUBLE PRECISION SBVL
      DOUBLE PRECISION DECAY,MNWO
      DOUBLE PRECISION DCYFCT,DCYT,DCYT2
      dimension mnw2(NMNWVL,MNWMAX),mnwnod(34,NODTOT)
      dimension ibound(NCOL,NROW,NLAY)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY),BOTM(NCOL,NROW,0:NBOTM),
     *  SRCVOL(NSCOL,NSROW,NSLAY)
      DIMENSION SRCMNW(NSCOL,NSROW,NSLAY),
     *  SNKMNW(NSCOL,NSROW,NSLAY),SOLMNW(NSCOL,NSROW,NSLAY)
C
      dimension MNWid(MNWMAX),MNWO(MNWMAX,7)
      DIMENSION SBVL(6,NIUNIT)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      character*16 text
      text = '            MNW2'
c     -----------------------------------------------------------------c
c-----if there are no wells do not process
      if(nmnw2.gt.0) then
C     if not using weighted method, reset SRCSOL to value that didn't include last 
c       move's MNW rates (so when added back it, we get the right number)
        IF(PTWTON.EQ.0) THEN
C LOOP OVER CELLS
          DO 20 KS=1,NSLAY
          DO 20 IS=1,NSROW
          DO 20 JS=1,NSCOL
           SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)-SOLMNW(JS,IS,KS)
  20      CONTINUE
        END IF
c     initialize SOLMNW and redefine each move, 
c       as it changes according to Cwell
      SOLMNW=0.0
c
c
c   Initialize sum terms
          Qnet = 0.000
          QSnk=0.0
          MSnk=0.0
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
c  only if in subgrid
        IF(MNW2(32,iw).eq.1) then 
c mnwo reset non-cumulative Qsnk and Qsrc counters
          if(MNWIin.gt.0) then
            MNWO(iw,3)=0.0
            MNWO(iw,5)=0.0
          end if
c   Loop over nodes in well
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
          do INODE=firstnode,lastnode
            k=MNWNOD(1,INODE)              
            i=MNWNOD(2,INODE)              
            j=MNWNOD(3,INODE)              
c
            JS=J-ISCOL1+1
            IS=I-ISROW1+1
            KS=K-ISLAY1+1
c
c   A very large # in WL reference array (8,m) flags multi-node wells
c
c        process multi-node wells
              if(abs(MNW2(2,iw)).gt.0) then
c
c   Loop 1: compute Qnet, get C' values
c
                if( ibound(j,i,k) .ne. 0 ) then
c   Qnet is net flow in or out of MNW well
                  qnode = MNWNOD(4,INODE)
                  Qnet  = Qnet  + qnode
c   For losing node (sink), sum mass coming out.  
              if(qnode.lt.0.0) then
cgzh check SUMWT and restrict qnode if necessary
c                qwt=(-1)*qnode*TIMV
c                IF(qwt.GT.SUMWT(JS,IS,KS)) qnode=-SUMWT(JS,IS,KS)/TIMV
c   These are pieces of eq. 1 and eq. 3 of GWT-MNW report (Qin)
                    QSnk=QSnk+qnode
c   here, at beginning of move, CONC is C from last time step
                    MSnk=MSnk+qnode*CONC(JS,IS,KS)
c mnwo
                    if (MNWIin.gt.0) then 
c m,3 is zeroed out above, so it is only cumulative over this time period
c m,5 is npt zeroed out, so it is cumulative over the whole simulation
                      MNWO(iw,3)=MNWO(iw,3)+qnode*CONC(JS,IS,KS)*TIMV
                      MNWO(iw,4)=MNWO(iw,4)+qnode*CONC(JS,IS,KS)*TIMV
                    end if
c  add all nodes to MNW item budget
                    SBVL(2,24)=SBVL(2,24)+qnode*CONC(JS,IS,KS)*TIMV
                  end if
                end if
              end if
          enddo
c   End Loop 1
c
         Qnet=mnw2(18,iw)
c
c     Call routing routine to determine concentrations in well
          CALL GWT1MNW2rt(nmnw2,MNWMAX,mnw2,NODTOT,MNWNOD,
     +             MNWid,CONC,WELLID,
     +             ncol,nrow,nlay,iouts,nscol,nsrow,nslay,
     +             Cwell,THCK,BOTM,NBOTM,iw,Qnet,IBOUND,NMNWVL)
c
c   Budget and MNWO stuff for differing Qnets
c
          if(Qnet.gt.0.0) then
            QCT=Qnet*MNW2(12,iw)*TIMV
            if (MNWIin.gt.0) then 
              MNWO(iw,1)=QCT    
              MNWO(iw,2)=MNWO(iw,2)+QCT    
            end if
c   Separate mass balance item: Net mass flux into aquifer by all MNWs 
            SBVL(1,25)=SBVL(1,25)+QCT
c   Also track flow through borehole
c    For net sources, that is just the sum of the sinks 
            SBVL(3,24)=SBVL(3,24)-MSnk*TIMV
c
          else if (Qnet.lt.0.0) then
c
            QCT=Qnet*MNW2(31,iw)*TIMV
            if (MNWIin.gt.0) then 
              MNWO(iw,1)=QCT    
              MNWO(iw,2)=MNWO(iw,2)+QCT    
            end if
c   Separate mass balance item: Net mass flux out of aquifer by all MNWs 
            SBVL(2,25)=SBVL(2,25)+QCT
c   Also track flow through borehole
c    For net sinks, this is just the "remaining" sources in the well 
            SBVL(3,24)=SBVL(3,24)+QCT-MSnk*TIMV
c
          else if (Qnet.eq.0) then
c   Also track flow through borehole
c    For Qnet=0 this is either the sum source mass or sink mass (here we use sinks) 
            SBVL(3,24)=SBVL(3,24)-MSnk*TIMV
c   Mass balance for unweighted versions
	                SBVL(2,24)=SBVL(2,24)+MSnk*TIMV
            if (MNWIin.gt.0) then 
			    MNWO(iw,1)=0.0
            end if
		  end if            
c     Save CWELL in MNWNOD(32,INODE) for output via MNW package:QSUM and lst file
c
c   Loop 2: now that we have Cwell, use it for each source node
c   Accumulate SOLMNW
c   Loop over nodes in this MNW
          do INODE= firstnode, lastnode
c     Save CWELL in MNWNOD(32,INODE) for output via MNW package:BYNODE (individual rates)
		   MNWNOD(32,INODE)=Cwell(INODE)             
c get node location 
           k=MNWNOD(1,INODE)              
           i=MNWNOD(2,INODE)              
           j=MNWNOD(3,INODE)              
c
           JS=J-ISCOL1+1
           IS=I-ISROW1+1
           KS=K-ISLAY1+1
           if( ibound(j,i,k) .ne. 0 ) then
             qnode = MNWNOD(4,INODE)
                IF(qnode.GT.0.0) THEN
                   qnodeC=qnode*Cwell(INODE)
                 if (MNWIin.gt.0) then 
c m,5 is zeroed out above, so it is only cumulative over this time period
c m,6 is npt zeroed out, so it is cumulative over the whole simulation
                   MNWO(iw,5)=MNWO(iw,5)+qnodeC*TIMV
                   MNWO(iw,6)=MNWO(iw,6)+qnodeC*TIMV
                 end if
c  add all nodes to MNW item budget
                 SBVL(1,24)=SBVL(1,24)+qnodeC*TIMV
                 SOLMNW(JS,IS,KS)=SOLMNW(JS,IS,KS)+qnodeC
c this doesn't accumlate SRCSOL: SRCMNW is taken out at the beginning 
c of this routine
c debug output
                 IF(PTWTON.EQ.0) SRCSOL(JS,IS,KS)=
     *                           SRCSOL(JS,IS,KS)+qnodeC
                END IF
            end if
c   End Loop 2
          enddo
c
        else
C   For single-node well, set sinks to CONC and srcs to Cinput
          qnode = MNWNOD(4,INODE)
		if(qnode.LT.0.0) then
		  Cwell(INODE)=CONC(JS,IS,KS)
          else
            Cwell(INODE)=MNW2(12,iw)
          end if
c     Save CWELL in MNWNOD(32,INODE) for output via MNW package:QSUM and lst file (summed rates)
          MNW2(31,iw)=Cwell(INODE)
c end if inside subgrid
         end if
c end if active well
       end if
c   End Loop over all wells
      enddo
C
C RE-COMPUTE SRCAVC FOR MOCWT
      IF(PTWTON.EQ.1) THEN
      DCYT2=1.D0
C ONLY DECAY OVER 1/2 TIMV
      IF(DECAY.NE.0.D0) THEN
        DCYFCT=DBLE(TIMV)*DECAY
        DCYT2=DEXP(-DCYFCT*0.5D0)
      END IF
C LOOP OVER CELLS
      DO 30 KS=1,NSLAY
      DO 30 IS=1,NSROW
      DO 30 JS=1,NSCOL
C ONLY IF MNW SOURCE AT THIS CELL (ONLY > 0 for a complex well)
        IF(SOLMNW(JS,IS,KS).GT.0.) THEN
C GET SPATIALLY-VARYING DECAY RATE
          IF(INDK.GT.0) THEN
            IF(IDKFO.EQ.1.OR.IDKFS.EQ.1) THEN
              TIMV2=TIMV*0.5
              CALL DK6DK(DCYT2,DKFO,DKFS,RF,TIMV2,
     *          IDKFO,IDKFS,
     *          JS,IS,KS,NSCOL,NSROW,NSLAY)
            END IF
          END IF 
C DECAY SOURCE (DCYT2=1.0 if no decay)
          MNWD=SOLMNW(JS,IS,KS)*DCYT2
c orig          SRCDCY=SRCDCY+(SRCD-SRCSOL(JS,IS,KS)+(BDYD-BDYSOL(JS,IS,KS)))
cgzh convert SRCDCY to mass by mult by time
          SRCDCY=SRCDCY+TIMV*(MNWD-SOLMNW(JS,IS,KS))
          SRCVOL(JS,IS,KS)=SRCVOL(JS,IS,KS)+SRCMNW(JS,IS,KS)
          SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+MNWD
cgzh store decayed mass for SRCAVC calculation
          SOLMNW(JS,IS,KS)=MNWD
        END IF
   30 CONTINUE   
      END IF
c  end if for nwell2>0
      endif
c-----return
      return
      end
c
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW2rt(nmnw2,MNWMAX,mnw2,NODTOT,MNWNOD,
     +             MNWid,CONC,WELLID,
     +             ncol,nrow,nlay,iouts,nscol,nsrow,nslay,
     +             Cwell,THCK,BOTM,NBOTM,iw,Qnet,IBOUND,NMNWVL)
C
c     ******************************************************************
c     route solute in mnw wells 
c     ******************************************************************
c
      CHARACTER*20 WELLID
      double precision MNW2,Qnet,MNWNOD,verysmall,Caq
      double precision q4,qqnet,Cprime,q27top,C27t,q27bot,C27b
      integer CC,firstnode,lastnode,inode
      real LBH
C TEMPORARY ARRAYS
      ALLOCATABLE CC(:),LBH(:),betweennodes(:)
      dimension mnw2(NMNWVL,MNWMAX),mnwnod(34,NODTOT),WELLID(MNWMAX)
      dimension Cwell(NODTOT)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY)
      dimension MNWid(MNWMAX),BOTM(NCOL,NROW,0:NBOTM)
      dimension IBOUND(NCOL,NROW,NLAY)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
c
      ALLOCATE(CC(NODTOT),LBH(NODTOT),
     *   betweennodes(NODTOT),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,1700)ISTAT
 1700   FORMAT(1X,'ALLOCATION OF MNW ROUTING ARRAYS FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
C Initialize
	C=0.0
	CC=0
cgzh 3/31/09 recoded routing with PUMPLOC capability
c          
c   Get node and cell info
      firstnode=MNW2(4,iw)
      lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
      INODE=firstnode
 8    k=MNWNOD(1,INODE)              
      i=MNWNOD(2,INODE)              
      j=MNWNOD(3,INODE)              
c   First check for no-flow cells at top of well, set these to CNOFLO and
c     continue from the first non-dry cell
      if(IBOUND(J,I,K).EQ.0) then
        Cwell(INODE)=CNOFLO
        CC(INODE)=6
c  If all nodes are no-flow, report this and return
        if(INODE.eq.LASTNODE) then
          write(iouts,*) 'For MNW well ',WELLID(iw)
          write(iouts,*) 'All nodes in this MNW well are no-flow'
          GOTO 111
	  end if
        INODE=INODE+1
        GO TO 8
	end if
c
c   Check first active node (not IBOUND=0) for flow down
c     if so, define Cwell here 
      if(MNWNOD(27,INODE+1).LE.0.0) then
         knd=INODE
      if(cc(knd).eq.0.0) then
         q4=0.0
         Caq=0.0
         qqnet=0.0
         Cprime=0.0
         q27top=0.0
         C27t=0.0
         if (MNWNOD(4,knd).lt.0.0) then
            q4=-MNWNOD(4,knd)
            CALL GETNODE(MNWNOD,NODTOT,knd,j,i,k,js,is,ks)
            if (CONC(JS,IS,KS).gt.0.0) Caq=CONC(JS,IS,KS)
         end if
         if (MNWNOD(27,knd).lt.0.0) then
            q27top=-MNWNOD(27,knd)
            C27t=MNW2(12,iw)
         end if
         if (MNWNOD(33,knd).ne.0.0.and.Qnet.gt.0.0) then
            qqnet=qnet
            Cprime=MNW2(12,iw)
         end if
c rbw begin change         
!         Cwell(knd)=(q4*Caq+qqnet*Cprime+q27top*C27t)/
!     &               (q4+qqnet+q27top)
         if (q4+qqnet+q27top.ne.0) then
           Cwell(knd)=(q4*Caq+qqnet*Cprime+q27top*C27t)/
     &               (q4+qqnet+q27top)
         else
           Cwell(knd)=Caq
         endif
c rbw end change     
         cc(knd)=1
      end if
      else
c     If no flow downward, skip sweep down, go to sweep up
        GO TO 12          
      end if
c   If only one node, skip sweep down
      if(firstnode.eq.lastnode) GO TO 12
c   If first node defined, proceed down well
      INODE=firstnode+1
c
c   Sweep down
c   Check for downward flow (<0) or no-flow (=0) on lower borehole face
c   For both of those conditions, depending on pump location, can define Cwell here
c
      do while (MNWNOD(27,INODE).LE.0.0)
c   Check node for flow out to aquifer
        if(MNWNOD(4,INODE).GE.0) then
c   if so, and no injection well here, set C to C at last node
		  if(MNWNOD(33,inode).eq.0.or.
     &      (MNWNOD(33,inode).eq.1.and.Qnet.LE.0.0)) then
		    Cwell(INODE)=Cwell(INODE-1)
            CC(INODE)=1
          else
c   if injection well, mix that with QBH from above
            Cwell(INODE)=(Qnet*MNW2(12,iw)-
     &        MNWNOD(27,INODE)*CWELL(INODE-1))/(Qnet-MNWNOD(27,INODE))
          end if
        else
c   Otherwise (there is flow in from aquifer here) 
c   if so, and no injection well here, mix C last node and Caq
          if(MNWNOD(33,inode).eq.0.or.
     &      (MNWNOD(33,inode).eq.1.and.Qnet.LE.0.0)) then
            CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
            Caq = CONC(JS,IS,KS)
            if (Caq.lt.0.0) Caq=0.0
c   Flux in from above = MNWNOD(27,INODE) (negative)
c   Flux in from aq = MNWNOD(4,INODE) (negative)
c   all terms are negative, Cwell should be positive
              Cwell(INODE)=(MNWNOD(27,INODE)*Cwell(INODE-1)+ 
     &         MNWNOD(4,INODE)*Caq)/(MNWNOD(27,INODE)+MNWNOD(4,INODE)) 
              CC(INODE)=1
          else
c   if flow in from aq, and injection well here, mix C last node, Caq, and Cinj
            CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
            Caq = CONC(JS,IS,KS)
            if (Caq.lt.0.0) Caq=0.0
            Cwell(INODE)=(Qnet*MNW2(12,iw)-MNWNOD(4,INODE)*
     &        Caq-MNWNOD(27,INODE)*CWELL(INODE-1))/
     &        (Qnet-MNWNOD(4,INODE)-MNWNOD(27,INODE))
          end if
        end if
c   If the lower borehole face is no flow (=0), leave sweep down
        if(MNWNOD(27,INODE+1).eq.0.0) GO TO 10
        INODE=INODE+1
c   Leave if past last node
        if(INODE.gt.lastnode) GO TO 10
      end do          
  10  CONTINUE 
c
c  Sweep up
  12  INODE=lastnode
c
c  Check last node for for flow in from aq and flow up
c  if so, try to define Cwell here 
      if(MNWNOD(4,INODE).LE.0.0.AND.MNWNOD(27,INODE).GE.0.0) then
c   If pump is not located here, or pump is here and withdrawal/zero, use CONC
        if(MNWNOD(33,inode).eq.0.or.
     &    (MNWNOD(33,inode).eq.1.and.Qnet.LE.0.0)) then
          CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
          caq = CONC(JS,IS,KS)
          if (CONC(JS,IS,KS).lt.0.0) Caq=0.0
		  Cwell(INODE)=Caq  
C   CC flags nodes where C is defined
          CC(INODE)=1 
C   If pump is here, and it is injection, mix C' and CONC              
        elseif (MNWNOD(33,inode).eq.1.and.Qnet.gt.0.0) then
          CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
          caq = CONC(JS,IS,KS)
          if (CONC(JS,IS,KS).lt.0.0) Caq=0.0
          Cwell(INODE)=(Qnet*MNW2(12,iw)-MNWNOD(4,INODE)*
     &      Caq)/(Qnet-MNWNOD(4,INODE))
C   CC flags nodes where C is defined
          CC(INODE)=1
        end if
      else
c     If no flow out to aq or no flow downward, skip sweep up
        GO TO 20          
      end if
c   If only one node, skip sweep up
      if(firstnode.eq.lastnode) GO TO 20
c   If last node defined, proceed up well
      INODE=lastnode-1
c
c   Begin sweep up
c   Check for upward flow (>0) or no-flow (=0) on upper borehole face
c   For both of those conditions, depending on pump location, can define Cwell here
c
      do while (MNWNOD(27,INODE).GE.0.0)
c   Check node for flow out to aquifer
        if(MNWNOD(4,INODE).GE.0) then
c   if so, and no injection well here, set C to C at last node
		  if(MNWNOD(33,inode).eq.0.or.
     &      (MNWNOD(33,inode).eq.1.and.Qnet.LE.0.0)) then
		    Cwell(INODE)=Cwell(INODE+1)
            CC(INODE)=1
          else
c   if injection well, mix that with QBH from below
            Cwell(INODE)=(Qnet*MNW2(12,iw)-
     &      MNWNOD(27,INODE+1)*CWELL(INODE+1))/(Qnet-MNWNOD(27,INODE+1))
          end if
        else
c   Otherwise (there is flow in from aquifer here) 
c   if so, and no injection well here, mix C last node and Caq
          if(MNWNOD(33,inode).eq.0.or.
     &      (MNWNOD(33,inode).eq.1.and.Qnet.LE.0.0)) then
            CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
            Caq = CONC(JS,IS,KS)
            if (Caq.lt.0.0) Caq=0.0
c   Flux in from below = MNWNOD(27,INODE+1) (positive)
c   Flux in from aq = MNWNOD(4,INODE) (negative)
              Cwell(INODE)=(MNWNOD(27,INODE+1)*Cwell(INODE+1)- 
     &         MNWNOD(4,INODE)*Caq)/(MNWNOD(27,INODE+1)-MNWNOD(4,INODE))
              CC(INODE)=1
          else
c   if flow in from aq, and injection well here, mix C last node, Caq, and Cinj
            CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
            Caq = CONC(JS,IS,KS)
            if (Caq.lt.0.0) Caq=0.0
            Cwell(INODE)=(Qnet*MNW2(12,iw)-MNWNOD(4,INODE)*
     &        Caq+MNWNOD(27,INODE+1)*CWELL(INODE+1))/
     &        (Qnet-MNWNOD(4,INODE)+MNWNOD(27,INODE+1))
          end if
        end if
c   Proceed to next node up
        INODE=INODE-1
c   Leave if past first node
        if(INODE.lt.firstnode) GO TO 20
      end do          
C
c   End of initial sweep up and down
C

  20  CONTINUE 
c     Check to see if all node C's defined using CC flag array
      IOK=1
      do 30 INODE=firstnode,lastnode
        if(CC(INODE).EQ.0) IOK=0
  30  continue
c     If all C's defined, skip to end
      if(IOK.EQ.1) GO TO 140
C
c     Flag strong sources to borehole with a 2 in CC
      do INODE = firstnode+1, lastnode-1
        if(CC(INODE).EQ.0) then
c     If flow down from above and up from below
          if(MNWNOD(27,INODE).GE.0.
     &      AND.MNWNOD(27,INODE+1).LE.0.0) then
            CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
            Caq = CONC(JS,IS,KS)
            if (Caq.lt.0.0) Caq=0.0
c     Check pump location 
            if (MNWNOD(33,INODE).eq.1.and.Qnet.gt.0.0) then
              Cwell(INODE)=(Qnet*MNW2(12,iw)-
     &        MNWNOD(4,INODE)*Caq)/(Qnet-MNWNOD(4,INODE))
            else
              Cwell(INODE)=Caq
            end if
c     For a node that is stagnation point, still set Cwell here but 
c     do not flag as strong source
            if(MNWNOD(27,INODE).EQ.0.0
     &        .AND.MNWNOD(27,INODE+1).EQ.0.0) THEN
		      CC(INODE)=4
c
            else			  
		      CC(INODE)=2
		    end if
		  end if
        end if
      enddo
c
c
c                 
c     From strong sources, sweep up for undefined nodes 
          do INODE = firstnode+1, lastnode-1 
c     CC=2 flags strong source 
            if(CC(INODE).EQ.2) then 
              knd=INODE-1 
              knd1=knd+1 
c     Only proceed if the strong source was not "one-sided" in the down direction 
c       So, check to make sure flow is going up from strong source 
             if(MNWNOD(27,INODE).gt.0) then     
              do while (knd.ge.firstnode) 
c     If undefined, try to define 
                if(CC(knd).eq.0) then 
c       If flow from above is downward, strong sink.  Cwell(knd) may not 
c         be defined, so abort and handle later. 
                  if(MNWNOD(27,knd).lt.0.0) then 
                    goto 2220 
c       If flow is upwards or no flow 
                  else if(MNWNOD(27,knd).ge.0.0) then 
c         If flow is in from aquifer, mix 
                   if(MNWNOD(4,knd).lt.0) then 
                      CALL GETNODE(MNWNOD,NODTOT,KND,j,i,k,js,is,ks) 
                           if(CONC(JS,IS,KS).ge.0.0) then 
                            Caq=CONC(JS,IS,KS) 
                           else 
                            Caq=0.0 
                           end if 
c                Check pump location 
                    if (MNWNOD(33,knd).ne.0.and.Qnet.gt.0.0) then 
                     Cwell(knd)=(Qnet*MNW2(12,iw)-MNWNOD(4,knd)*Caq 
     &                 +MNWNOD(27,knd1)*Cwell(knd1))/ 
     &                 (Qnet-MNWNOD(4,knd)+MNWNOD(27,knd1)) 
                     CC(knd)=3 
                    else 
                      Cwell(knd)=(MNWNOD(27,knd1)*Cwell(knd1)- 
     &                      MNWNOD(4,knd)*Caq)/ 
     &                     (MNWNOD(27,knd1)-MNWNOD(4,knd)) 
                      CC(knd)=3 
                    end if 
c        If no flow in from aquifer, just set = Cwell from below 
c        [else if(MNWNOD(4,k).ge.0) then] 
                   else 
                    if (MNWNOD(33,knd).ne.0.and.Qnet.gt.0.0) then 
                     Cwell(knd)=(Qnet*MNW2(12,iw) 
     &                 +MNWNOD(27,knd1)*Cwell(knd1))/ 
     &                 (Qnet+MNWNOD(27,knd1)) 
                    else 
                      Cwell(knd)=Cwell(knd1) 
                      CC(knd)=3 
                      end if 
                   end if 
                    end if 
c        If no flow into the cell above, leave sweep up 
                  if(MNWNOD(27,knd).eq.0.0) GO TO 2220           
               ELSE 
                    go to 2220 
                 end if 
                knd=knd-1 
              end do 
 2220         continue 
             end if 
            end if 
c         end sweep up for strong sources 
          enddo 
c                 
C     now check first node and calculate Cwell(1) if not already defined
         knd=firstnode
      if(cc(knd).eq.0.0) then
         q4=0.0
         Caq=0.0
         qqnet=0.0
         Cprime=0.0
         q27top=0.0
         C27t=0.0
         q27bot=0.0
         C27b=0.0
         if (MNWNOD(4,knd).lt.0.0) then
            q4=-MNWNOD(4,knd)
            CALL GETNODE(MNWNOD,NODTOT,knd,j,i,k,js,is,ks)
            if (CONC(JS,IS,KS).gt.0.0) Caq=CONC(JS,IS,KS)
         end if
         if (MNWNOD(27,knd).lt.0.0) then
            q27top=-MNWNOD(27,knd)
            C27t=MNW2(12,iw)
         end if
         if (MNWNOD(27,knd+1).gt.0.0) then
            q27bot=MNWNOD(27,knd+1)
            C27b=Cwell(knd+1)
         end if
         if (MNWNOD(33,knd).ne.0.0.and.Qnet.gt.0.0) then
            qqnet=qnet
            Cprime=MNW2(12,iw)
         end if
         Cwell(knd)=(q4*Caq+qqnet*Cprime+q27top*C27t+q27bot*C27b)/
     &               (q4+qqnet+q27top+q27bot)
         cc(knd)=3
      end if
c                 
c     From strong sources, sweep down for undefined nodes 
          do INODE = firstnode+1, lastnode-1 
c     CC=2 flags strong source 
            if(CC(INODE).EQ.2) then 
              knd=INODE+1     
              knd1=knd+1 
              knd1m=knd-1 
c     Only proceed if the strong source was not "one-sided" in the up direction 
c       So, check to make sure flow is going down from strong source 
             if(MNWNOD(27,KND).le.0) then 
              do while (knd.lt.lastnode) 
c     If undefined, try to define 
                if(CC(knd).eq.0) then 
c       If flow from below is upward, strong sink.  Cwell(knd) may not 
c         be defined, so abort and handle later. 
c 
                  if(MNWNOD(27,knd1).gt.0.0) then 
                    go to 2221 
c       If flow is upwards or no flow 
c                     [else if(MNWNOD(27,knd+1).le.0.0) then] 
                  else 
c         If flow is in from aquifer, mix 
                   if(MNWNOD(4,knd).lt.0) then 
                      CALL GETNODE(MNWNOD,NODTOT,KND,j,i,k,js,is,ks) 
                           if(CONC(JS,IS,KS).ge.0.0) then 
                            Caq=CONC(JS,IS,KS) 
                           else 
                            Caq=0.0 
                           end if 
c                Check pump location 
                    if (MNWNOD(33,knd).ne.0.and.Qnet.gt.0.0) then 
                     Cwell(knd)=(Qnet*MNW2(12,iw)-MNWNOD(4,knd)*Caq 
     &                 +MNWNOD(27,knd)*Cwell(knd1m))/ 
     &                 (Qnet-MNWNOD(4,knd)+MNWNOD(27,knd)) 
                     CC(knd)=3 
                    else 
                      Cwell(knd)=(MNWNOD(27,knd)*Cwell(knd1m)- 
     &                      MNWNOD(4,knd)*Caq)/ 
     &                     (MNWNOD(27,knd)-MNWNOD(4,knd)) 
                      CC(knd)=3 
                    end if 
c        If no flow in from aquifer, just set = Cwell from above 
c          [else if(MNWNOD(4,k).ge.0) then] 
                   else 
                    if (MNWNOD(33,knd).ne.0.and.Qnet.gt.0.0) then 
                     Cwell(knd)=(Qnet*MNW2(12,iw) 
     &                 +MNWNOD(27,knd)*Cwell(knd1m))/ 
     &                 (Qnet+MNWNOD(27,knd)) 
                    else 
                      Cwell(knd)=Cwell(knd1m) 
                      CC(knd)=3 
                      end if 
                   end if 
                    end if 
c        If no flow into the cell below, leave sweep down 
                  if(MNWNOD(27,knd1).eq.0.0) GO TO 2221           
               ELSE 
                    go to 2221 
                 end if 
                knd=knd+1 
              end do 
 2221         continue 
c        If at last node, just set equal to Cwell above (if not pumploc) 
              if(knd.eq.lastnode.and.CC(knd).eq.0) then 
               if (MNWNOD(33,knd).ne.0.and.Qnet.gt.0.0) then 
                          Cwell(knd)=(Qnet*MNW2(12,iw) 
     &                 +MNWNOD(27,knd)*Cwell(knd-1))/ 
     &                 (Qnet+MNWNOD(27,knd)) 
                CC(knd)=3 
               else 
                          Cwell(knd)=Cwell(knd-1) 
                CC(knd)=3 
               end if 
              end if 
             end if 
            end if 
c         end sweep down for strong sources 
          enddo 
C Strong sinks
c
      do INODE = firstnode+1, lastnode-1
c If upper and lower nodes are sources to node, 
c node is strong sink in borehole
        if((MNWNOD(27,INODE).lt.0).and.(MNWNOD(27,INODE+1).gt.0)) then
c   If pump is not located here, or pump is here and withdrawal/zero, use 
c      conc in well above and below
          if(MNWNOD(33,inode).eq.0.or.
     &    (MNWNOD(33,inode).eq.1.and.Qnet.LE.0.0)) then
c   If also flow in from aquifer, use that too
            if(MNWNOD(4,INODE).LT.0.0) then
              CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
              Caq = CONC(JS,IS,KS)
              if (Caq.lt.0.0) Caq=0.0
              Cwell(INODE)=(-MNWNOD(27,INODE)*Cwell(INODE-1)+
     &                      MNWNOD(27,INODE+1)*Cwell(INODE+1)-
     &                      MNWNOD(4,INODE)*Caq)/
     &                     (-MNWNOD(27,INODE)+MNWNOD(27,INODE+1)
     &                     -MNWNOD(4,INODE))
              CC(INODE)=5
            else
c   If flow out to aq, in from above & below
              Cwell(INODE)=(-MNWNOD(27,INODE)*Cwell(INODE-1)+
     &                      MNWNOD(27,INODE+1)*Cwell(INODE+1))/
     &                     (-MNWNOD(27,INODE)+MNWNOD(27,INODE+1))
              CC(INODE)=5
            endif
          else if (MNWNOD(33,inode).eq.1.and.Qnet.gt.0.0) then
c   If pump here and injection, mix with C'
            Cwell(INODE)=(-MNWNOD(27,INODE)*Cwell(INODE-1)+
     &                    MNWNOD(27,INODE+1)*Cwell(INODE+1)+
     &                    Qnet*MNW2(12,iw))/
     &                   (-MNWNOD(27,INODE)+MNWNOD(27,INODE+1)+Qnet)
            CC(INODE)=5
          end if
        end if
      end do
C
C
 100  CONTINUE
c     Check to see if all node C's defined using CC flag array
      do 40 INODE=firstnode,lastnode-1
        if(CC(INODE).EQ.0) then
cljk just set Cwell to zero for small flows instead of stopping
           if (abs(MNWNOD(4,INODE)).lt.0.00005) then
             Cwell(INODE) = 0.0
           else
             write(iouts,*)  
             write(iouts,*) '***ERROR*** C at node in MNW not
     & defined by routing algorithm'
             write(iouts,*) 'firstnode,inode',firstnode,inode
             write(iouts,*) 'MNWNOD(4,INODE)',MNWNOD(4,INODE)
             stop 'MNW routing error'
           end if
         end if
  40  continue
 140  continue
C
c     for a withdrawal well--use 1st node for MNWO output conc
c
         if(Qnet.lt.0) then
          do inode=firstnode, lastnode
            if(MNWNOD(33,inode).eq.1) MNW2(31,iw)=Cwell(inode)
          end do
         else 
c     otherwise, compute avg c for mnwo based on borehole length, this output goes to MNWO
c     save in first node of MNW (m)
c
c     compute length associated with each section (LBH)
c
c     determine distance between nodes (betweennodes)
c
          do INODE = firstnode,lastnode-1
            CALL GETNODE(MNWNOD,NODTOT,INODE,j1,i1,k1,js1,is1,ks1)
            CALL GETNODE(MNWNOD,NODTOT,INODE+1,j2,i2,k2,js2,is2,ks2)
C     convert to real coodinates
            x1=js1*CDEL+(0.5*CDEL)
            y1=is1*RDEL+(0.5*RDEL)
            z1=(0.5)*THCK(JS1,IS1,ks1)
		    if(ks1.gt.1) then
		      do 45 kk=1,ks1-1
               z1=z1+THCK(js1,is1,kk)
   45         continue
            end if
		    x2=js2*CDEL+(0.5*CDEL)
            y2=is2*RDEL+(0.5*RDEL)
            z2=(0.5)*THCK(JS2,IS2,ks2)
		    if(ks2.gt.1) then
		      do 46 kk=1,ks2-1
               z2=z2+THCK(js2,is2,kk)
   46         continue
            end if
c     caculate distance between nodes
        betweennodes(INODE)=SQRT(((x1-x2)**2)+((y2-y2)**2)+((z1-z2)**2))
	    end do
c
c     estimate length of borehole
c
c first node segment
          LBH(firstnode)=(0.5*betweennodes(firstnode))
c extend to the water table (head in well), or land surface, whichever is lower
c set top equal to water level in well node 1
          top=MNW2(17,iw)          
c get top of land surface
          CALL GETNODE(MNWNOD,NODTOT,firstnode,j,i,k,js,is,ks)
          j=js+ISCOL1-1
          i=is+ISROW1-1
          surftop=BOTM(J,I,LBOTM(1)-1)
          if(surftop.lt.top) top=surftop
          heightnode=BOTM(J,I,LBOTM(ks)-1)-(0.5*THCK(JS,IS,KS))
          LBH(firstnode)=LBH(firstnode)+top-heightnode
c last node segment
c
c check last two nodes, compute length based on difference in coords in lay, col, row
          CALL GETNODE(MNWNOD,NODTOT,lastnode,j,i,k,js1,is1,ks1)
          CALL GETNODE(MNWNOD,NODTOT,lastnode-1,j,i,k,js2,is2,ks2)
c js1 etc is last node
c js2 etc is next-to-last node
          dirj=0.
          diri=0.
          dirk=0.
          if(js1.ne.js2) dirj=1.
          if(is1.ne.is2) diri=1.
          if(ks1.ne.ks2) dirk=1.
C     convert to real coodinates
          x1=js1*CDEL+(0.5*CDEL)
          y1=is1*RDEL+(0.5*RDEL)
          z1=(0.5)*THCK(JS,IS,ks1)
            if(ks1.gt.1) then
		    do 47 kk=1,ks1-1
               z1=z1+THCK(js1,is1,kk)
   47         continue
            end if
          cellfacex=x1-(0.5*CDEL)
          cellfacey=y1-(0.5*RDEL)
          cellfacez=z1-((0.5)*THCK(JS,IS,ks1))
          finalseg=SQRT( (((x1-cellfacex)**2)*dirj)+
     &                   (((y1-cellfacey)**2)*diri)+
     &                   (((z1-cellfacez)**2)*dirk) )
          LBH(lastnode)=(betweennodes(lastnode-1)*0.5)+finalseg

c
          do INODE = firstnode+1,lastnode-1
            LBH(INODE)=0.5*(betweennodes(INODE-1)+betweennodes(INODE))
          end do
c     sum mass and length
          suml=0.0
          summ=0.0
          do INODE = firstnode,lastnode
            suml=suml+LBH(INODE)
            summ=summ+LBH(INODE)*Cwell(INODE)
          end do
c     calulate average, save in well2(12), which MNWO and MNW uses to output
          MNW2(31,iw)=summ/suml   
          
         end if
 111  CONTINUE
c
      DEALLOCATE(CC,LBH,betweennodes)
C
      RETURN
C
      END
C
C
C
C  GETNODE      TRANSLATE MNW NODE TO SUBGRID COORDINATES
C
C     ******************************************************************
C
      SUBROUTINE GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
C
C     ******************************************************************
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c
      double precision MNWNOD
      dimension mnwnod(34,NODTOT)
c
c get node location 
          k=MNWNOD(1,INODE)              
          i=MNWNOD(2,INODE)              
          j=MNWNOD(3,INODE)              
c
          IS=I-ISROW1+1
          JS=J-ISCOL1+1
          KS=K-ISLAY1+1
C
      RETURN
      END
C
C
C
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW2sm(nmnw2,MNWMAX,mnw2,mnwnod,nodtot,
     +             ibound,ncol,nrow,nlay,
     +             iouts,CAVG,nscol,nsrow,nslay,
     +             MNWid,TIMV,SNKFLO,SNKMSOUT,
     *             SBVL,NIUNIT,Cwell,imov,
     *             MNWO,MNWIin,NMNWVL)
C
c     ******************************************************************
c     calculate solute mass leaving MNW sinks
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      double precision MNW2,MNWNOD
      dimension mnw2(NMNWVL,MNWMAX),mnwnod(34,NODTOT)
      double precision well2,QSnk,Msnk,qnode,SNKMSOUT
      dimension Cwell(NODTOT)
      REAL MNWD
      DOUBLE PRECISION SBVL,MNWO
      dimension well2(18,MNWMAX),ibound(ncol,nrow,nlay)
      DIMENSION CAVG(NSCOL,NSROW,NSLAY),SNKMSOUT(NSCOL,NSROW,NSLAY)
      DIMENSION SNKFLO(NSCOL,NSROW,NSLAY)
C
      dimension MNWid(MNWMAX),MNWO(MNWMAX,7)
CMOCWT
      DIMENSION SBVL(6,NIUNIT)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c     -----------------------------------------------------------------c
c-----if there are no wells do not process
      if(nmnw2.gt.0) then
c Loop over all MNW locations
C Sum sink fluxes
c   Loop over all wells
       do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
        if (MNW2(1,iw).EQ.1) then
c   Loop over nodes in well
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
CLJK REVERSE ORDER TO get node by node conc
c   Loop 1: compute Qnet, get C' values
c   Initialize sum terms
          Qnet = 0.000
          QSnk=0.0
          MSnk=0.0
          do INODE=lastnode,firstnode,-1
c get node location 
           CALL GETNODE(MNWNOD,NODTOT,INODE,j,i,k,js,is,ks)
           IF(MNW2(32,iw).eq.1) then 
c
c   Process multi-node wells
c
!            if(MNW2(2,iw).gt.0) then
            if(Abs(MNW2(2,iw)).gt.1) then
c
c  if sink
             if (mnw2(18,iw).lt.0) then
c
             if( ibound(j,i,k) .ne. 0 ) then
              qnode = MNWNOD(4,INODE)
c   For losing node (sink), sum mass coming out.  
              if(qnode.lt.0.0) then
c
               QSnk=QSnk+qnode
c
c   here, define mass that left well
                if(ptwton.eq.0) then
                  MSnk=MSnk+qnode*CAVG(JS,IS,KS)
                else
                  Msnk=Msnk+(qnode/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
                end if
cljk                set node-by-node conc for well with all negative q (next 3 lines) 
                if (Qsnk.lt.0) then 
                        if(ptwton.eq.1) then
                          Cwell(INODE)=(Msnk/TIMV)/Qsnk
                        else
                          Cwell(INODE)=Msnk/Qsnk
                        end if
cljk              MNWNOD(32,INODE) = cavg(JS,IS,KS)
                  MNWNOD(32,INODE) = Cwell(INODE)
                end if  
                if (MNWIin.gt.0) then 
c m,3 is zeroed out above, so it is only cumulative over this time period
c m,5 is npt zeroed out, so it is cumulative over the whole simulation
                        if(ptwton.eq.0) then
                    MNWO(iw,3)=MNWO(iw,3)+qnode*CAVG(JS,IS,KS)*TIMV
                    MNWO(iw,4)=MNWO(iw,4)+qnode*CAVG(JS,IS,KS)*TIMV
                        else
                          MNWO(iw,3)=MNWO(iw,3)+(qnode/SNKFLO(JS,IS,KS))
     *                         *SNKMSOUT(JS,IS,KS)
                          MNWO(iw,4)=MNWO(iw,4)+(qnode/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
                        end if
                end if
               end if
              end if
             end if
c
c   For weighted method that builds Msnk with SNKMSOUT, TIMV is included
             if(Qsnk.ne.0.0) then
        if(ptwton.eq.1) then
          Cwell(INODE)=(Msnk/TIMV)/Qsnk
        else
              Cwell(INODE)=Msnk/Qsnk
        end if
cljk       Save CWELL in MNWNOD(32,INODE) for output via MNW package:QSUM and lst file (summed rates) 
c
               if(PUMPLOC.ne.0) then
                 if(MNWNOD(33,INODE).eq.1) MNW2(31,iw)=Cwell(INODE)
               else
                 MNW2(31,iw)=Cwell(INODE) 
               end if
             end if
            else
C   For single-node well, set sinks to CAVG 
             qnode = MNWNOD(4,INODE)
		     if(qnode.LT.0.0) then
              if(ptwton.eq.0) then
                Cwell(firstnode)=CAVG(JS,IS,KS)
              else
                if (SNKFLO(JS,IS,KS).eq.0) then
                  Cwell(firstnode) = 0
                else
                  Cwell(firstnode)=
     1              SNKMSOUT(JS,IS,KS)/(SNKFLO(JS,IS,KS)*TIMV)
                endif
              end if
             end if
c   Save well concentration 
             MNW2(31,iw)=Cwell(INODE)
c   endif multi/single
            end if
c   endif in subgrid           
           end if
c   End Loop over nodes
          enddo
c
c   Endif active wells
        endif
c   End Loop over all wells
       enddo
C

c  end if for nmnw2>0
      endif
c-----return
      return
      end
