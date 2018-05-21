C
C
C GWT1MNW1AL ALLOCATE SPACE FOR MNWid
C
C     ******************************************************************
C
      SUBROUTINE GWT1MNW1AL(ISUMI,LSMNWI,mxwel2,IOUTS)    
C
C     ******************************************************************
C
C     ARRAY IS LAY,ROW,COL,IFACE
      IF(mxwel2.GT.0) THEN
        LSMNWI=ISUMI
        ISUMI=ISUMI+mxwel2
        WRITE(IOUTS,101) mxwel2
      ELSE
	  LSMNWI=1
      END IF
  101 FORMAT(1X,I8,' ELEMENTS IN IR ARRAY ARE USED BY GWT-MNW')
      RETURN
      END
C
c
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW1bd(nwell2,mxwel2,well2,ibound,ncol,nrow,nodes,
     +             iouts,conc,srcflo,srcsol,snkflo,nscol,nsrow,nslay,
     +             MNWid,SRCMNW,SNKMNW,SOLMNW,icomplex,KSTP,MNWsite,
     +             kkper,IQzero,MNWOon,IPERGWT)
C
c     ******************************************************************
c     calculate solute mass terms for mnw wells
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      double precision well2
      dimension MNWsite(mxwel2),IQzero(mxwel2)
      character*32 MNWsite
      dimension well2(18,mxwel2),ibound(nodes)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),
     *  SRCFLO(NSCOL,NSROW,NSLAY),SRCSOL(NSCOL,NSROW,NSLAY),
     *  SNKFLO(NSCOL,NSROW,NSLAY)
C
      dimension MNWid(mxwel2)
      DIMENSION SRCMNW(NSCOL,NSROW,NSLAY),SNKMNW(NSCOL,NSROW,NSLAY),
     *  SOLMNW(NSCOL,NSROW,NSLAY)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      character*16 text
      text = '             MNW'
c     ------------------------------------------------------------------
c
C Initialize SRCMNW,SNKMNW,AND SOLMNW
      SRCMNW=0.
      SNKMNW=0.
      SOLMNW=0.
      if(MNWOon.gt.0) IQzero=0
c-----if there are no wells do not process
      if(nwell2.gt.0) then
c Initialize complex MNW flags
        MNWid=0
c Initialize counter used to label MNWs
        id=0
c Initialize flag for printing conc only for 1st node of MNW
        idcheck=0
c Initialize flag for whether or not there are any complex MNWs
        icomplex=0
c Initialize print flag for header line to main file
        IPRNT=0
c
        m = 0
c Loop over all MNW locations
        do while( m .lt. nwell2 )
          m = m + 1
c get node location 
          n = INT(well2(1,m))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
c   A very large # in WL reference array (8,m) flags multi-node wells
c
        if( well2(8,m) .gt. 1.0E30 ) then
c Toggle flag for whether or not there are any complex MNWs
          icomplex=1
c   Set ne = number of last node in this MNW
          ne  = ifrl(well2(7,m))
c   Save m for other loops
          msave=m
c   Set flag to check for complex MNW 
          IQ=0
c   Set flag to check for Qnet=0 
          IQnet=0
c   Loop 1: check for MNW's that go across subgrid boundary
c   Flag MNW's with flow in AND out (complex)
c   Check for MNW's with all Q's=0 (Qnet=0)
c
c   Loop over nodes in this MNW
          do iin = m, ne
c get node location 
          n = INT(well2(1,iin))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
            if (iin.eq.m) then
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
c   Don't need to do this check if already flagged as complex
cgzh debug
c      write(iouts,*) 'MNWid,iin check 1:', MNWid(iin),iin

            IF(MNWid(iin).eq.0) then
             IF(well2(17,iin).NE.0) THEN
c   Record first flow that is positive or negative
              IF(IQ.EQ.0) THEN
                if(well2(17,iin).GT.0.0) then
                  IQ=1
                ELSE
                  IQ=-1              
                END IF
              ELSE
c   Record next flow in list as positive or negative
                if(well2(17,iin).GT.0.0) then
                  IQnext=1
                ELSE
                  IQnext=-1              
                END IF
c   If some nodes inflow and some nodes outflow, flag this well as complex MNW
                if(IQnext.ne.IQ) then
                  m2=msave
                  do iin2 = m2, ne
                     if(m2.eq.iin2) id=id+1                   
                     MNWid(iin2)=id
c get node location 
                     n = INT(well2(1,iin2))
                     k = (n-1)/(ncol*nrow) + 1
                     j = mod((n-1),ncol*nrow)/ncol + 1
                     i = mod((n-1),ncol) + 1
cgzh this will print all nodes of a complex mnw...
cgzh debug this is covered by print statements below
c                  if(m2.eq.iin2) write(iouts,*) 'The following nodes
c     & are part of a single complex MNW:'
c				write(iouts,'(A,3I6)') 'Layer, Row, Column = ', k,j,i                  
                  enddo
                end if            
              END IF
             END IF
            END IF
C  Check desired flow rate, if all zeros then IQnet will remain zero, flag
c  in MNWid will be set below
             IF(well2(2,iin).NE.0.0) IQnet=1
cgzh debug
c      write(iouts,*) 'IQnet check 1:', IQnet

c            
          enddo
c   End Loop 1
cgzh debug moving this out as an initialization above
c   Set Qnet flag for MNWO; default is 0 (not a Qnet=0), will be reset appropriately below
c          if(MNWOon.gt.0.and.kkper.eq.1) IQzero(m)=0
c
c   Set MNWid = -1 for Qnet = 0.  This will be when IQnet never gets
c   set to 1 in loop 1
cgzh debug
          iz=0
          if(IQnet.eq.0) then
            m2=msave
            do iin2 = m2, ne
c   Set Qnet flag for MNWO for 1st SP
               if(MNWOon.gt.0.and.kkper.eq.IPERGWT) IQzero(iin2)=1
c   Set MNWid flag for Qnet=0, this can change with time
               MNWid(iin2)=-1
cgzh debug ouput
c      write(iouts,*) 'Node ', iin2,' set to Qnet=0'
            iz=iz+1
            enddo
          end if            
c      if(iz.gt.0) write(iouts,*) iz,' nodes set to Qnet=0 for this well'
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
c   Reset m
          m=msave
c
c   Loop 3: process each source node in simple MNW
c   Loop over nodes in this MNW
          do iin = m, ne
c get node location 
          n = INT(well2(1,iin))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
            IF(KS.GE.1.AND.KS.LE.NSLAY
     &       .AND.IS.GE.1.AND.IS.LE.NSROW
     &       .AND.JS.GE.1.AND.JS.LE.NSCOL) then 
            if( ibound(n) .ne. 0 ) then
c   For non-complex MNWs
              qnode = well2(17,iin)
              if(MNWid(iin).eq.0) then      
                IF(qnode.LT.0.0) THEN
                 SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+qnode
c     Save Conc in well2(11,m) for output via MNW package:BYNODE (individual rates)
		       well2(11,iin)=CONC(JS,IS,KS)           
C   print output
                 IF(KSTP.EQ.1) 
     *		     WRITE(IOUTS,'(1X,I6,3I8,2X,1P1E12.4,16X,A)') 
     1		                 iin,K,I,J,qnode,MNWsite(iin)
                ELSE IF(qnode.GT.0.0) THEN
                 SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+qnode
                 qnodeC=qnode*well2(4,m)
                 SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+qnodeC
c     Save Conc in well2(11,m) for output via MNW package:BYNODE (individual rates)
		       well2(11,iin)=well2(4,m)  
C   print output
                 IF(KSTP.EQ.1) 
     *                 WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 iin,K,I,J,qnode,well2(4,m),MNWsite(iin)
                END IF
c   For complex MNWs
              else                 
			  IF(qnode.LT.0.0) THEN
                 SNKMNW(JS,IS,KS)=SNKMNW(JS,IS,KS)+qnode
                 IF(PTWTON.EQ.0) SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+qnode
C   print output (only print Cinput for first well in list)
                 IF(KSTP.EQ.1) then
                 IF(MNWid(iin).ne.idcheck) THEN
                   WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 iin,K,I,J,qnode,well2(4,m),MNWsite(iin)
                 ELSE
                   WRITE(IOUTS,'(1X,I6,3I8,2X,1P1E12.4,16X,A)') 
     1		                 iin,K,I,J,qnode,MNWsite(iin)
                 END IF
                 end if
			   idcheck=MNWid(iin)
                ELSE IF(qnode.GT.0.0) THEN
                 SRCMNW(JS,IS,KS)=SRCMNW(JS,IS,KS)+qnode
                 IF(PTWTON.EQ.0) SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+qnode
C   print output (only print Cinput for first well in list)
                 IF(KSTP.EQ.1) then
                 IF(MNWid(iin).ne.idcheck) THEN
                   WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 iin,K,I,J,qnode,well2(4,m),MNWsite(iin)
                 ELSE
                   WRITE(IOUTS,'(1X,I6,3I8,2X,1P1E12.4,16X,A)') 
     1		                 iin,K,I,J,qnode,MNWsite(iin)
                 END IF
                 end if
                 idcheck=MNWid(iin)
c   SOLMNW (based on Cwell, not Cinput) is updated each move in GWT1MNW1cx
                END IF
              end if
            end if
            end if
c   End Loop 3
          enddo
c   set counter to last node of MNW
          m = ne  
C
        else       !  End of multi-node conditioning IF statement
c     Single-node well
         IF(KS.LT.1.OR.KS.GT.NSLAY) GO TO 20
         IF(IS.LT.1.OR.IS.GT.NSROW) GO TO 20
         IF(JS.LT.1.OR.JS.GT.NSCOL) GO TO 20
         if( ibound(n) .eq. 0 ) GO TO 20
C
          qnode = well2(17,m)
          IF(qnode.LT.0.0) THEN
             SNKFLO(JS,IS,KS)=SNKFLO(JS,IS,KS)+qnode
C   print output
                 IF(KSTP.EQ.1) 
     *             WRITE(IOUTS,'(1X,I6,3I8,2X,1P1E12.4,16X,A)') 
     1		                 m,K,I,J,qnode,MNWsite(m)
          ELSE IF(qnode.GT.0.0) THEN
             SRCFLO(JS,IS,KS)=SRCFLO(JS,IS,KS)+qnode
             qnodeC=qnode*well2(4,m)
             if(well2(4,m).lt.0.0) then
               write(iouts,*) '***ERROR*** MNW concentration < 0 at
     & col,row,lay=', i,j,k
               write(iouts,*) 'Execution stopping.'
               STOP 'ERROR in GWT (MNW)'
             end if
              SRCSOL(JS,IS,KS)=SRCSOL(JS,IS,KS)+qnodeC
C   print output
                 IF(KSTP.EQ.1) 
     *             WRITE(IOUTS,'(1X,I6,3I8,2X,1P2E12.4,4X,A)') 
     1		                 m,K,I,J,qnode,well2(4,m),MNWsite(m)
          END IF
  20  CONTINUE
        endif
c
        enddo
c   End Loop over all wells
        endif
c-----return
      return
      end
c
c_______
c
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW1cx(nwell2,mxwel2,well2,ibound,ncol,nrow,nlay,
     +          nodes,iouts,conc,srcflo,srcsol,snkflo,nscol,nsrow,nslay,
     +             MNWid,SRCMNW,SNKMNW,SOLMNW,SRCVOL,
     *             RF,DECAY,TIMV,SRCDCY,
     *             DKFO,DKFS,INDK,IDKFO,IDKFS,
     *             SBVL,NIUNIT,MNWSITE,Cwell,SUMWT,THCK,imov,
     *             MNWO,MNWOon,IQzero,BOTM,NBOTM)
C
c     ******************************************************************
c     calculate solute mass terms for complex mnw wells
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      double precision well2,Qnet,QSnk,Msnk,qnode,qwt
      dimension MNWsite(mxwel2)
      dimension Cwell(mxwel2)
      character*32 MNWsite
      REAL MNWD,MNWO
      DOUBLE PRECISION SBVL,SUMWT
      DOUBLE PRECISION DECAY
      DOUBLE PRECISION DCYFCT,DCYT,DCYT2
      DIMENSION SUMWT(NSCOL,NSROW,NSLAY)
      dimension well2(18,mxwel2),ibound(nodes)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY),
     *  SRCSOL(NSCOL,NSROW,NSLAY),SRCVOL(NSCOL,NSROW,NSLAY),
     *  BOTM(NCOL,NROW,0:NBOTM)
      DIMENSION SRCMNW(NSCOL,NSROW,NSLAY),
     *  SNKMNW(NSCOL,NSROW,NSLAY),SOLMNW(NSCOL,NSROW,NSLAY)
C
      dimension MNWid(mxwel2),MNWO(mxwel2,7),IQzero(mxwel2)
CMOCWT
      DIMENSION SBVL(6,NIUNIT)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      character*16 text
      text = '             MNW'
c     -----------------------------------------------------------------c
c-----if there are no wells do not process
      if(nwell2.gt.0) then
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
        m = 0
c Loop over all MNW locations
C Compute Qnet and sum sink fluxes
        do while( m .lt. nwell2 )
         m = m + 1
c
c get node location 
          n = INT(well2(1,m))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
            IF(KS.GE.1.AND.KS.LE.NSLAY
     &       .AND.IS.GE.1.AND.IS.LE.NSROW
     &       .AND.JS.GE.1.AND.JS.LE.NSCOL) then 
c
cgzh mnwo reset non-cumulative Qsnk and Qsrc counters
          if(MNWOon.gt.0) then
		  MNWO(m,3)=0.0
		  MNWO(m,5)=0.0
          end if
c   A very large # in WL reference array (8,m) flags multi-node wells
c
        if( well2(8,m) .gt. 1.0E30 ) then
c   Set ne = number of last node in this MNW
          ne  = ifrl(well2(7,m))
c   Save m for other loops
          msave=m
c
c only process Cwell for complex wells
cgzh this moved 3/2/2006, simple MNW sinks now handled after advection
        if (MNWid(m).ne.0) then
c   Loop 1: compute Qnet, get C' values
c   Initialize sum terms
          Qnet = 0.000
          QSnk=0.0
          MSnk=0.0
C   Cinput is now well2(4,m)
          Cinput=0.0
c   Loop over nodes in this MNW
          do iin = m, ne
c get node location 
          n = INT(well2(1,iin))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
c
            if( ibound(n) .ne. 0 ) then
c   Qnet is net flow in or out of MNW well
              qnode = well2(17,iin)
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
cljk                set node-by-node conc for well with all negative q (next 3 lines) 
                    if ((mnwid(m).eq.0).and.(qnet.lt.0)) then 
                      well2(11,iin) = conc(js,is,ks) 
                    end if  
cgzh mnwo
                if (MNWOon.gt.0) then 
c m,3 is zeroed out above, so it is only cumulative over this time period
c m,5 is npt zeroed out, so it is cumulative over the whole simulation
                  MNWO(m,3)=MNWO(m,3)+qnode*CONC(JS,IS,KS)*TIMV
                  MNWO(m,4)=MNWO(m,4)+qnode*CONC(JS,IS,KS)*TIMV
                end if
cgzh  add all nodes to MNW item budget
                SBVL(4,24)=SBVL(4,24)+qnode*CONC(JS,IS,KS)*TIMV
              end if
            end if
          enddo
c   End Loop 1
c
c   Apply Qnet=0 flag to Qnet
      IF(MNWid(m).lt.0) Qnet=0.
c
c only process Cwell for complex wells
cgzh debug with MNWid including -1 for qnet=0 changed this to ne 0
c         if (MNWid(m).gt.0) then
cgzh this moved 3/2/2006, simple MNW sinks now handled after advection
c        if (MNWid(m).ne.0) then
cgzh debug CUT large section dealing with complex Qnet >< 0
cgzh now, route all complex wells
c     Call routing routine to determine concentrations in well
            CALL GWT1MNW1rt(nwell2,mxwel2,well2,MNWid,CONC,
     +             ncol,nrow,nlay,iouts,nscol,nsrow,nslay,
     +             Cwell,THCK,BOTM,NBOTM,m,Qnet,IBOUND)
c
c   Budget and MNWO stuff for differing Qnets
c
          if(Qnet.gt.0.0) then
            QCT=Qnet*well2(4,m)*TIMV
cgzh mnwo
            if (MNWOon.gt.0.and.IQzero(m).eq.0) then 
              MNWO(m,1)=QCT    
              MNWO(m,2)=MNWO(m,2)+QCT    
            end if
c   Separate mass balance item: Net mass flux into aquifer by all MNWs 
            SBVL(1,25)=SBVL(1,25)+QCT
c   Also track flow through borehole
c    For net sources, that is just the sum of the sinks 
            SBVL(3,24)=SBVL(3,24)-MSnk*TIMV
cgzh mnwo
            if (MNWOon.gt.0.and.IQzero(m).eq.1) then 
              MNWO(m,1)=-MSnk*TIMV    
              MNWO(m,2)=MNWO(m,2)-MSnk*TIMV
            end if
c
          else if (Qnet.lt.0.0) then
c
            QCT=Qnet*Cwell(m)*TIMV
cgzh mnwo
            if (MNWOon.gt.0) then 
              MNWO(m,1)=QCT    
              MNWO(m,2)=MNWO(m,2)+QCT    
            end if
c   Separate mass balance item: Net mass flux out of aquifer by all MNWs 
            SBVL(2,25)=SBVL(2,25)+QCT
c   Also track flow through borehole
c    For net sinks, this is just the "remaining" sources in the well 
            SBVL(3,24)=SBVL(3,24)+QCT-MSnk*TIMV
cgzh mnwo
            if (MNWOon.gt.0.and.IQzero(m).eq.1) then 
              MNWO(m,1)=-MSnk*TIMV    
              MNWO(m,2)=MNWO(m,2)-MSnk*TIMV
            end if
c
          else if (Qnet.eq.0.0) then
c   Also track flow through borehole
c    For Qnet=0 this is either the sum source mass or sink mass (here we use sinks) 
            SBVL(3,24)=SBVL(3,24)-MSnk*TIMV
c   Mass balance for unweighted versions
            SBVL(2,24)=SBVL(2,24)+MSnk*TIMV
cgzh mnwo
            if (MNWOon.gt.0) then 
c     only include if this well was Qnet=0 from the start 
                if(IQzero(m).eq.1) then
			    MNWO(m,1)=-MSnk*TIMV
                  MNWO(m,2)=MNWO(m,2)-MSnk*TIMV
c     if this is a well that became Qnet 0 after the beginning, set non-cumulative spot to 0
                else
			    MNWO(m,1)=0
                end if
		  end if
		end if            
c     Save CWELL in well2(12,m) for output via MNW package:QSUM and lst file
cgzh removing due to inclusion in routing routine          well2(12,m)=Cwell(m)
c   Reset m
          m=msave
c
c   Loop 2: now that we have Cwell, use it for each source node
c   Accumulate SOLMNW
c   Loop over nodes in this MNW
          do iin = m, ne
c     Save CWELL in well2(11,m) for output via MNW package:BYNODE (individual rates)
		   well2(11,iin)=Cwell(iin)             
c get node location 
          n = INT(well2(1,iin))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
            if( ibound(n) .ne. 0 ) then
                qnode = well2(17,iin)
                IF(qnode.GT.0.0) THEN
                   qnodeC=qnode*Cwell(iin)
cgzh mnwo
                 if (MNWOon.gt.0) then 
c m,5 is zeroed out above, so it is only cumulative over this time period
c m,6 is npt zeroed out, so it is cumulative over the whole simulation
                   MNWO(m,5)=MNWO(m,5)+qnodeC*TIMV
                   MNWO(m,6)=MNWO(m,6)+qnodeC*TIMV
                 end if
cgzh debug add all nodes to MNW item budget
                 SBVL(1,24)=SBVL(1,24)+qnodeC*TIMV
                 SOLMNW(JS,IS,KS)=SOLMNW(JS,IS,KS)+qnodeC
cgzh this doesn't accumlate SRCSOL: SRCMNW is taken out at the beginning 
cgzh of this routine
                 IF(PTWTON.EQ.0) SRCSOL(JS,IS,KS)=
     *                           SRCSOL(JS,IS,KS)+qnodeC
                END IF
            end if
c   End Loop 2
          enddo
c
         else
C   For non-complex well, mix water in well and compute Conc to report
cljk            **** This should be lt rather than gt  ************************** 
cljk            **** if Qsnk is negative (pumping) use msnk/qsnk **************** 
cljk            **** if Qsnk is postive (injection) use C from MNW file  ******** 
cljk            If(Qsnk.GT.0.0) THEN 
            If(Qsnk.lT.0.0) THEN 
		    Cwell(m)=MSnk/QSnk
            else
              Cwell(m)=well2(4,m)
            end if
cljk          **** add this line so printed output will be right  ***** 
cljk          set well2(11,m) above (individual rates) 
cljk          Save CWELL in well2(12,m) for output via MNW package:QSUM and lst file (summed rates) 
            well2(12,m)=Cwell(m) 
         endif
c   set counter to last node of MNW
          m = ne  
        else
C   For single-node well, set sinks to CONC and srcs to Cinput
          qnode = well2(17,m)
		if(qnode.LT.0.0) then
		  Cwell(m)=CONC(JS,IS,KS)
          else
            Cwell(m)=well2(4,m)
          end if
c     Save CWELL in well2(11,m) for output via MNW package:BYNODE file (inidivual rates)
          well2(11,m)=Cwell(m)
c     Save CWELL in well2(12,m) for output via MNW package:QSUM and lst file (summed rates)
          well2(12,m)=Cwell(m)
c end if inside subgrid
        end if
c       end if multi-node
        end if
        enddo
c   End Loop over all wells
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
      SUBROUTINE GWT1MNW1rt(nwell2,mxwel2,well2,MNWid,CONC,
     +             ncol,nrow,nlay,iouts,nscol,nsrow,nslay,
     +             Cwell,THCK,BOTM,NBOTM,m,Qnet,IBOUND)
C
c     ******************************************************************
c     route solute in complex mnw well 
c     ******************************************************************
c
      double precision well2,Qnet
      integer CC
      real LBH
C TEMPORARY ARRAYS
      ALLOCATABLE QBH(:),CC(:),LBH(:),betweennodes(:)
      dimension Cwell(mxwel2)
      dimension well2(18,mxwel2)
      DIMENSION CONC(NSCOL,NSROW,NSLAY),THCK(NSCOL,NSROW,NSLAY)
      dimension MNWid(mxwel2),BOTM(NCOL,NROW,0:NBOTM)
      dimension IBOUND(NCOL,NROW,NLAY)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
c
      ALLOCATE(QBH(mxwel2),CC(mxwel2),LBH(mxwel2),
     *   betweennodes(mxwel2),STAT=ISTAT)
      IF (ISTAT.NE.0) THEN
        WRITE(*,1700)ISTAT
 1700   FORMAT(1X,'ALLOCATION OF MNW ROUTING ARRAYS FAILED,',
     &  ' RETURNED ERROR MESSAGE NUMBER: ',I6)
        CALL USTOP(' ')
      ENDIF
C Initialize
      QBH=0.0
	C=0.0
	CC=0
c   Set ne = last node in this MNW
          ne  = ifrl(well2(7,m))
c   Save m for other loops
          msave=m
c   QBH is flow in between nodes in borehole
c   Set flux from first node equal to flux at that node and Qnet, which is applied 
c    at first node only
          QBH(m)=well2(17,m)-Qnet
c   Loop over other nodes in this MNW to set QBH
          do iin = m+1, ne-1
c   QBH between successive nodes is Q at previous node - Q at node
            QBH(iin)=QBH(iin-1)+well2(17,iin)
          enddo
c   Set QBH at last node =0, this will be used in sweep down nodes
c
          QBH(ne)=0.0
c
c   Calculate distribution of conc in well bore
c
c
cgzh 10/5/06
c   First check for no-flow cells at top of well, set these to CNOFLO and
c     continue from the first non-dry cell
c          
c
   8      CALL GETNODE1(well2,mxwel2,m,js,is,ks,ncol,nrow,nlay)
          J=JS+ISCOL1-1
          I=IS+ISROW1-1
          K=KS+ISLAY1-1
          if(IBOUND(J,I,K).EQ.0) then
            Cwell(m)=CNOFLO
            CC(m)=6
c  If all nodes are no-flow, report this and return
            if(m.eq.ne) then
             write(iouts,*) 'For MNW well index m=',m
             write(iouts,*) 'All nodes in this MNW well are no-flow'
             GOTO 111
	      end if
            m=m+1
            GO TO 8
	    end if
c
c   Check first active node for no flow in borehole.  
		if(QBH(m).eq.0.0) then
c     If flow in from aquifer, or no flow, use CONC
            if(well2(17,m).LE.0.0) then
              CALL GETNODE1(well2,mxwel2,m,js,is,ks,ncol,nrow,nlay)
		    Cwell(m)=CONC(JS,IS,KS)              
c     Otherwise flow must be from Qnet, so use C'input
            else
              Cwell(m)=well2(4,m)
            end if
C   CC flags nodes where C is defined
            CC(m)=1
c     If no flow downward, skip over sweep down
            GO TO 12
          end if
c
c   Check first node for downward flow, if so then define Cwell there
c
          if(QBH(m).LT.0) then
c   If outflow to aquifer or no flow, Cwell is C'input here
           if(well2(17,m).GE.0) then
             Cwell(m)=well2(4,m)
           else
c   Otherwise, inflow from aquifer
            CALL GETNODE1(well2,mxwel2,m,js,is,ks,ncol,nrow,nlay)
c   Mix with Qnet*C' if a net source (note aq flux sign is negative)
            if(Qnet.gt.0.) then
              Cwell(m)=(Qnet*well2(4,m)-well2(17,m)*CONC(JS,IS,KS))/
     &                   (Qnet-well2(17,m))
cgzh debug output
              write(iouts,*) 'Case 1: Qnet>0 in routing routine'
c   If Qnet is 0 or outward (<0), use aquifer conc
            else
		    Cwell(m)=CONC(JS,IS,KS)
            end if
C   CC flags nodes where C is defined
            CC(m)=1
cgzh debug
      if(cwell(m).gt.0.0) then
	 continue
	end if
c   If first node defined, proceed down well
            iin=m+1
c   If only one node, skip sweep down
            if(iin.gt.ne) GO TO 12
c
c   Begin sweep down
c   Check for downward flow (<0) or no-flow (=0) on lower borehole face
c   For both of those conditions, can define Cwell here
c
            do while (QBH(iin).LE.0)
c   Check node for inflow from aquifer, if none then set C to C at last node
              if(well2(17,iin).GT.0) then
			  Cwell(iin)=Cwell(iin-1)
                CC(iin)=1
              else
c   Otherwise (there is flow in from aquifer here) mix C at last node with Caq 
                CALL GETNODE1(well2,mxwel2,iin,js,is,ks,
     &               ncol,nrow,nlay)
c   Flux in from above = QBH(iin-1) (negative)
c   Flux in from aq = well2(17,iin) (negative)
c   all terms are negative, Cwell should be positive
                Cwell(iin)=(QBH(iin-1)*Cwell(iin-1)+ 
     &                  well2(17,iin)*CONC(JS,IS,KS))/ 
     &                 (QBH(iin-1)+well2(17,iin)) 
                CC(iin)=1
              end if
c   If the lower borehole face was no flow (=0), leave sweep down
c   Last node is caught here as QBH(ne)=0.0
              if(QBH(iin).eq.0.0) GO TO 10
              iin=iin+1
            end do          
  10        CONTINUE 
           end if
          end if
c
c   Check last node for inflow or no flow, if so then set C to Caq
c
  12      if(well2(17,ne).LE.0) then
            CALL GETNODE1(well2,mxwel2,ne,js,is,ks,ncol,nrow,nlay)
		  Cwell(ne)=CONC(JS,IS,KS)
            CC(ne)=1
c   If last node defined, proceed up well
            iin=ne-1
c   If only one node, skip sweep up
            if(iin.lt.m) GO TO 20
c   If borehole flow on upward face of last node is zero, skip sweep
            if(QBH(iin).eq.0.0) GO TO 20
c   If at first node, go to special section to handle that case
            if(iin.eq.m) GO TO 18
c
c   Begin sweep up
c   Check upper borehole face for upward or no flow
c
            do while (QBH(iin-1).GE.0)
c   Check for no inflow, if so just use last well conc
              if(well2(17,iin).GE.0) then
 			  Cwell(iin)=Cwell(iin+1)
                CC(iin)=1
              else
c   Otherwise (there is flow in from aquifer here) mix C at last node with Caq
c   Note sign difference in two flux terms requires "-" not "+" 
                CALL GETNODE1(well2,mxwel2,iin,js,is,ks,
     &               ncol,nrow,nlay)
                Cwell(iin)=(QBH(iin)*Cwell(iin+1)-
     &                  well2(17,iin)*CONC(JS,IS,KS))/
     &                 (QBH(iin)-well2(17,iin)) 
                
                CC(iin)=1
              end if
              iin=iin-1
c   if reach first node, go to special section for that case
              if(iin.eq.m) GO TO 18
            end do
c   if never got to first node, skip that section
            GO TO 20
c
c   Begin section for special case: swept up all the way to first node, iin=m
c
  18        CONTINUE
c   Check for no aquifer inflow
            if(well2(17,m).GT.0) then
c   Mix water coming up with Qnet>0 water 
              if(Qnet.gt.0.0) then
                Cwell(m)=(Qnet*well2(4,m)+QBH(m)*Cwell(m+1))/
     &                   (Qnet+QBH(m))
cgzh debug output
c                write(iouts,*) 'Case 2: Qnet>0 in routing routine'
              else
c   If Qnet<=0, just use Cwell from below
 			  Cwell(m)=Cwell(m+1)
              end if
                CC(m)=1
            else
c   Otherwise well2(17,m)<0 (there is flow in from aquifer here) mix C at last node with Caq 
              CALL GETNODE1(well2,mxwel2,m,js,is,ks,
     &               ncol,nrow,nlay)
              if(Qnet.gt.0.0) then
cgzh debug this would be 3 inputs with no outputs, error!  
                STOP 'Qnet ERROR in MNW routing routine'
              else
cgzh debug if no flow at all, don't divide by zero, set to no-flow conc
                if((QBH(m)-well2(17,m)).EQ.0) then
cgzh debug or should this be cwell=conc?
                  Cwell(m)=CNOFLO
                else
                  Cwell(m)=(QBH(m)*Cwell(m+1)-
     &                  well2(17,m)*CONC(JS,IS,KS))/
     &                 (QBH(m)-well2(17,m))                 
			  end if
              end if
              CC(m)=1
            end if
c
c     end section for special case: swept up all the way to first node
c
  20        CONTINUE 
          end if
C
c     Check to see if all node C's defined using CC flag array
          IOK=1
          do 30 iin=m,ne
           if(CC(iin).EQ.0) IOK=0
  30      continue
c     If all C's defined, skip 
          if(IOK.EQ.1) GO TO 140
C
c     Flag strong sources to borehole with a 2 in CC
          do iin = m+1, ne-1
            if(CC(iin).EQ.0) then
              if(QBH(iin-1).GE.0.AND.QBH(iin).LE.0.0) then
                CALL GETNODE1(well2,mxwel2,iin,js,is,ks,
     &               ncol,nrow,nlay)
			  Cwell(iin)=CONC(JS,IS,KS)
c       For a node that is stagnation point, still set Cwell here but 
c         do not flag as strong source
                if(QBH(iin-1).EQ.0.0.AND.QBH(iin).EQ.0.0) THEN
			    CC(iin)=4
                else			  
			    CC(iin)=2
			  end if
			end if
            end if
          enddo
c                
c     From strong sources, sweep up for undefined nodes
          do iin = m+1, ne-1
c     CC=2 flags strong source
            if(CC(iin).EQ.2) then
              k=iin-1
c     Only procede if the strong source was not "one-sided" in the down direction
c       So, check to make sure flow is going up from strong source
             if(QBH(iin-1).gt.0) then     
              do while (k.gt.m) 
c     If undefined, try to define
                if(CC(k).eq.0) then
c       If flow from above is downward, strong sink.  Cwell(k-1) may not
c         be defined, so abort and handle later.
                  if(QBH(k-1).lt.0.0) then
                    goto 2220
c       If flow is upwards or no flow
                  else if(QBH(k-1).ge.0.0) then
c         If flow is in from aquifer, mix
                    if(well2(17,k).lt.0) then 
                      CALL GETNODE1(well2,mxwel2,k,js,is,ks,
     &                  ncol,nrow,nlay)
                      Cwell(k)=(QBH(k)*Cwell(k+1)-
     &                      well2(17,k)*CONC(JS,IS,KS))/
     &                     (QBH(k)-well2(17,k))
                      CC(k)=3
c        If no flow in from aquifer, just set = Cwell from below
                    else if(well2(17,k).ge.0) then
                      Cwell(k)=Cwell(k+1)
                      CC(k)=3
	              end if
	            end if
c        If the cell above was no flow, leave sweep up
                  if(QBH(k-1).eq.0.0) GO TO 2220          
			  end if              
                k=k-1
              end do
 2220         continue
c can get to k=m from do loop (swept all the way up) or if iin=m+1 (str src in node 2)
              if(k.eq.m.and.CC(k).eq.0) then
c sweeping up from a strong source and reaching the 1st node
c   has several possible outcomes depending on Qnet and flow from well node
c
c   if qnet=0, well flow must go out node into aquifer, so use Cwell below
               if(Qnet.eq.0.0) then
			  Cwell(k)=Cwell(k+1)
cgzh debug output
                write(iouts,*) 'Qnet=0 in node 1 above a strong source'
c   if Qnet is out and well flow is out, use Cwell below
               elseif(Qnet.lt.0.0.and.well2(17,m).ge.0.0) then
			  Cwell(k)=Cwell(k+1)
cgzh debug output
                write(iouts,*) 'Qnet<0 and well flux ge 0 
     & in node 1 above a strong source'
c   if Qnet is out and well flow is in, mix QBH and well flow
               elseif(Qnet.lt.0.0.and.well2(17,m).lt.0.0) then
                CALL GETNODE1(well2,mxwel2,k,js,is,ks,ncol,nrow,nlay)
                Cwell(k)=(Cwell(k+1)*QBH(k)-well2(17,m)*CONC(JS,IS,KS))/
     &                    (QBH(k)-well2(17,m))
cgzh debug output
                write(iouts,*) 'Qnet<0 and well flux lt 0 
     & in node 1 above a strong source'
c   if Qnet is in, well flow must be out, so mix QBH and Qnet
               elseif(Qnet.gt.0.0) then
                Cwell(k)=(Cwell(k+1)*QBH(k)+Qnet*well2(4,m))/
     &                    (QBH(k)+Qnet)                
cgzh debug output
                write(iouts,*) 'Qnet>0 (well flux must be out) 
     & in node 1 above a strong source'
               else
                 write(iouts,*) 'Error in MNW routing routine, stopping'
                 STOP
               end if
                CC(k)=3
              end if
             end if
            end if
c         end sweep up for strong sources
          enddo
c                
c                
c     From strong sources, sweep down for undefined nodes
          do iin = m+1, ne-1
            if(CC(iin).EQ.2) then
              k=iin+1     
c     Only procede if the strong source was not "one-sided" in the up direction
c       So, check to make sure flow is going down from strong source
             if(QBH(iin).lt.0) then
              do while (k.lt.ne) 
                if(CC(k).eq.0) then
c       If flow from below is upward, strong sink.  Cwell(k+1) may not
c         be defined, so abort and handle later.
c
cgzh the indices below were fixed, k+1 was wrong in last version
c                  if(QBH(k+1).gt.0.0) then
                  if(QBH(k).gt.0.0) then
                    goto 2221
c       If flow is downwards or no flow
c                  else if(QBH(k+1).le.0.0) then
                  else if(QBH(k).le.0.0) then
c         If flow is in from aquifer, mix
                    if(well2(17,k).lt.0) then 
                      CALL GETNODE1(well2,mxwel2,k,js,is,ks,
     &                  ncol,nrow,nlay)
                      Cwell(k)=(QBH(k-1)*Cwell(k-1)- 
     &                      well2(17,k)*CONC(JS,IS,KS))/ 
     &                     (QBH(k-1)-well2(17,k)) 
                      CC(k)=3
c        If no flow in from aquifer, just set = Cwell from above
                    else if(well2(17,k).ge.0) then
                      Cwell(k)=Cwell(k-1)
                      CC(k)=3
	              end if
	            end if
c        If the cell below was no flow, leave sweep down
                  if(QBH(k).eq.0.0) GO TO 2221          
                end if              
                k=k+1
              end do
 2221         continue
c        If at last node, just set equal to Cwell above
              if(k.eq.ne.and.CC(k).eq.0) then
			  Cwell(k)=Cwell(k-1)
                CC(k)=3
              end if
             end if
            end if
c         end sweep down for strong sources
          enddo
c                
C Strong sinks
c Handle case where upper and lower nodes are sources to node, 
c node is strong sink in borehole
c
          do iin = m+1,ne-1
cgzh changed le to lt, get rid of =0 case, then not a strong sink, handled above
            if((QBH(iin-1).lt.0).and.(QBH(iin).gt.0)) then
                    Cwell(iin)=(QBH(iin)*Cwell(iin+1)-
     &                    QBH(iin-1)*Cwell(iin-1))/
     &                   (QBH(iin)-QBH(iin-1))
              CC(iin)=5
            end if
          end do
C
C
 100  CONTINUE
c     Check to see if all node C's defined using CC flag array
          do 40 iin=m,ne
           if(CC(iin).EQ.0) then
cljk just set Cwell to zero for small flows instead of stopping
cgzh should not need this, notify us if you get the error below! (commented out cwell=0 setting)
cgzh 10/2/06 brian clark problem did get an error here...uncommenting
             if (abs(well2(17,iin)).lt.0.00005) then
               Cwell(iin) = 0.0
             else
               write(iouts,*)  
               write(iouts,*) '***ERROR*** C at node in MNW not
     & defined by routing algorithm'
               write(iouts,*) 'm,iin',m,iin
               write(iouts,*) 'well2(17,iin)',well2(17,iin)
               stop 'MNW routing error'
             end if
           end if
  40      continue
 140      continue
C
c     for a withdrawal well--use 1st node for MNWO output conc
c
         if(Qnet.lt.0) then
          well2(12,m)=Cwell(m)
         else 
c     otherwise, compute avg c for mnwo based on borehole length, this output goes to MNWO
c     save in first node of MNW (m)
c
c     compute length associated with each section (LBH)
c
c     determine distance between nodes (betweennodes)
          do iin = m,ne-1
            CALL GETNODE1(well2,mxwel2,iin,js1,is1,ks1,ncol,nrow,nlay)
            CALL GETNODE1(well2,mxwel2,iin+1,js2,is2,ks2,ncol,nrow,nlay)
C     convert to real coodinates
            x1=js1*CDEL+(0.5*CDEL)
            y1=is1*RDEL+(0.5*RDEL)
            z1=(0.5)*THCK(JS,IS,ks1)
		  if(ks1.gt.1) then
		    do 45 kk=1,ks1-1
               z1=z1+THCK(js1,is1,kk)
   45         continue
            end if
		  x2=js2*CDEL+(0.5*CDEL)
            y2=is2*RDEL+(0.5*RDEL)
            z2=(0.5)*THCK(JS,IS,ks2)
		  if(ks2.gt.1) then
		    do 46 kk=1,ks2-1
               z2=z2+THCK(js1,is1,kk)
   46         continue
            end if
c     caculate distance between nodes
          betweennodes(iin)=SQRT(((x1-x2)**2)+((y2-y2)**2)+((z1-z2)**2))
	    end do
c
c     estimate length of borehole
c
c first node segment
          LBH(m)=(0.5*betweennodes(m))
c extend to the water table (head in well), or land surface, whichever is lower
c set top equal to water level in well node 1
          top=well2(10,m)          
c get top of land surface
          CALL GETNODE1(well2,mxwel2,m,js,is,ks,ncol,nrow,nlay)
          j=js+ISCOL1-1
          i=is+ISROW1-1
          surftop=BOTM(J,I,LBOTM(1)-1)
          if(surftop.lt.top) top=surftop
          heightnode=BOTM(J,I,LBOTM(ks)-1)-(0.5*THCK(JS,IS,KS))
          LBH(m)=LBH(m)+top-heightnode
c last node segment
c
c check last two nodes, compute length based on difference in coords in lay, col, row
          CALL GETNODE1(well2,mxwel2,ne,js1,is1,ks1,ncol,nrow,nlay)
          CALL GETNODE1(well2,mxwel2,ne-1,js2,is2,ks2,ncol,nrow,nlay)
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
          LBH(ne)=betweennodes(ne-1)+finalseg

cgzh old last segment stuff below
c          LBH(ne)=betweennodes(ne-1)
c check that this length does not exceed 1/2 of the "appropriate" cell dim
c          CALL GETNODE1(well2,mxwel2,ne-1,js1,is1,ks1,ncol,nrow,nlay)
c          CALL GETNODE1(well2,mxwel2,ne,js2,is2,ks2,ncol,nrow,nlay)
c          if(js1.eq.js2.and.is1.eq.is2.and.abs(ks1-ks2).eq.1) then
c            zlimit=THCK(js2,is2,ks2)
c            if(LBH(ne).gt.zlimit) LBH(ne)=zlimit
c          else
c            zlimit=CDEL
c            if(RDEL.gt.zlimit) zlimit=RDEL
c            if(THCK(js2,is2,ks2).gt.zlimit) zlimit=THCK(js2,is2,ks2)
c            if(LBH(ne).gt.zlimit) LBH(ne)=zlimit
c          end if
cgzh end old last segment stuff
c
          do iin = m+1,ne-1
            LBH(iin)=0.5*(betweennodes(iin-1)+betweennodes(iin))
          end do
c     sum mass and length
          suml=0.0
          summ=0.0
          do iin = m,ne
            suml=suml+LBH(iin)
            summ=summ+LBH(iin)*Cwell(iin)
          end do
c     calulate average, save in well2(12), which MNWO and MNW uses to output
          well2(12,m)=summ/suml   
         end if
 111  CONTINUE
c
      DEALLOCATE(QBH,CC,LBH,betweennodes)
C
      RETURN
C
      END
C
C
C
C  GETNODE1      TRANSLATE MNW NODE TO SUBGRID COORDINATES
C
C     ******************************************************************
C
      SUBROUTINE GETNODE1(well2,mxwel2,m,js,is,ks,ncol,nrow,nlay)
C
C     ******************************************************************
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c
      double precision well2
      dimension well2(18,mxwel2)
c
c get node location 
          n = INT(well2(1,m))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
C
      RETURN
      END
C
C
c_________________________________________________________________________________
c
      SUBROUTINE GWT1MNW1sm(nwell2,mxwel2,well2,ibound,ncol,nrow,nodes,
     +             iouts,CAVG,nscol,nsrow,nslay,
     +             MNWid,TIMV,SNKFLO,SNKMSOUT,
     *             SBVL,NIUNIT,Cwell,imov,
     *             MNWO,MNWOon)
C
c     ******************************************************************
c     calculate solute mass leaving simple MNW sinks
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      double precision well2,QSnk,Msnk,qnode,SNKMSOUT
      dimension Cwell(mxwel2)
      REAL MNWD,MNWO
      DOUBLE PRECISION SBVL
      dimension well2(18,mxwel2),ibound(nodes)
      DIMENSION CAVG(NSCOL,NSROW,NSLAY),SNKMSOUT(NSCOL,NSROW,NSLAY)
      DIMENSION SNKFLO(NSCOL,NSROW,NSLAY)
C
      dimension MNWid(mxwel2),MNWO(mxwel2,7)
CMOCWT
      DIMENSION SBVL(6,NIUNIT)
CMOCWT
      INCLUDE 'ptwt.inc'
c
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c     -----------------------------------------------------------------c
c-----if there are no wells do not process
      if(nwell2.gt.0) then
        m = 0
c Loop over all MNW locations
C Sum simple sink fluxes
        do while( m .lt. nwell2 )
         m = m + 1
c get node location 
          n = INT(well2(1,m))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
            IF(KS.GE.1.AND.KS.LE.NSLAY
     &       .AND.IS.GE.1.AND.IS.LE.NSROW
     &       .AND.JS.GE.1.AND.JS.LE.NSCOL) then 
c
c   A very large # in WL reference array (8,m) flags multi-node wells
c
        if( well2(8,m) .gt. 1.0E30 ) then
c   Set ne = number of last node in this MNW
          ne  = ifrl(well2(7,m))
c
c only process simple sinks here
        if (MNWid(m).eq.0) then
c   Loop 1: compute Qnet, get C' values
c   Initialize sum terms
          Qnet = 0.000
          QSnk=0.0
          MSnk=0.0
CLJK                do iin = m, ne
CLJK REVERSE ORDER TO get node by node conc
                do iin =ne,m,-1
c get node location 
          n = INT(well2(1,iin))
          k = (n-1)/(ncol*nrow) + 1
          j = mod((n-1),ncol*nrow)/ncol + 1
          i = mod((n-1),ncol) + 1
c
          JS=I-ISCOL1+1
          IS=J-ISROW1+1
          KS=K-ISLAY1+1
c
            if( ibound(n) .ne. 0 ) then
              qnode = well2(17,iin)
c   For losing node (sink), sum mass coming out.  
              if(qnode.lt.0.0) then
c
                QSnk=QSnk+qnode
c
c   here, define mass that left well
                if(ptwton.eq.0) then
                  MSnk=MSnk+qnode*CAVG(JS,IS,KS)
cgzh  add all nodes to MNW item budget
c                  SBVL(4,24)=SBVL(4,24)+qnode*CAVG(JS,IS,KS)*TIMV
                else
                  Msnk=Msnk+(qnode/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
c                  SBVL(4,24)=SBVL(4,24)+(qnode/SNKFLO(JS,IS,KS))
c     *                        *SNKMSOUT(JS,IS,KS)
                end if
cljk                  set node-by-node conc for well with all negative q (next 3 lines)
                      if (Qsnk.lt.0) then
cljk set cwell (for mnwo) and well2(11,iin) (for bynode)
cljk as long as pump is assumed above well and looping order is reversed, this should
cljk give the correct node-by-node concentrations  (used highlighted logic (orange) below                      
                        if(ptwton.eq.1) then
                          Cwell(iin)=(Msnk/TIMV)/Qsnk
                        else
                          Cwell(iin)=Msnk/Qsnk
                        end if
cljk                        well2(11,iin) = cavg(JS,IS,KS)
                        well2(11,iin) = cwell(iin)
                      end if
cgzh                  mnwo
                      if (MNWOon.gt.0) then
c                       m,3 is zeroed out above, so it is only cumulative over this time period
c                       m,5 is npt zeroed out, so it is cumulative over the whole simulation
                        if(ptwton.eq.0) then
                          MNWO(m,3)=MNWO(m,3)+qnode*CAVG(JS,IS,KS)*TIMV
                          MNWO(m,4)=MNWO(m,4)+qnode*CAVG(JS,IS,KS)*TIMV
                        else
                          MNWO(m,3)=MNWO(m,3)+(qnode/SNKFLO(JS,IS,KS))
     *                         *SNKMSOUT(JS,IS,KS)
                          MNWO(m,4)=MNWO(m,4)+(qnode/SNKFLO(JS,IS,KS))
     *                        *SNKMSOUT(JS,IS,KS)
                        end if
                      end if
                    end if
                  end if
                enddo
c
c   Save well concentration in first spot (m) of Cwell
c   For weighted method that build Msnk with SNKMSOUT, TIMV is included
      if(Qsnk.ne.0.0) then
        if(ptwton.eq.1) then
          Cwell(m)=(Msnk/TIMV)/Qsnk
        else
          Cwell(m)=Msnk/Qsnk
        end if
cljk          Save CWELL in well2(12,m) for output via MNW package:QSUM and lst file (summed rates) 
        well2(12,m)=Cwell(m) 
      end if
c end if simple well
        end if
c
      m=ne
        else
C   For single-node well, set sinks to CAVG 
          qnode = well2(17,m)
		if(qnode.LT.0.0) then
            if(ptwton.eq.0) then
              Cwell(m)=CAVG(JS,IS,KS)
            else
              Cwell(m)=SNKMSOUT(JS,IS,KS)/(SNKFLO(JS,IS,KS)*TIMV)
            end if
          end if
c     Save CWELL in well2(11,m) for output via MNW package:BYNODE file (inidivual rates)
          well2(11,m)=Cwell(m)
c     Save CWELL in well2(12,m) for output via MNW package:QSUM and lst file (summed rates)
          well2(12,m)=Cwell(m)
c end if inside subgrid
        end if
c       end if multi-node
        end if
        enddo
c   End Loop over all wells
C

c  end if for nwell2>0
      endif
c-----return
      return
      end
