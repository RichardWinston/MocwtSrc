c   GZH  20090405  -- Conversion from MNW1 to MNW2
c
      SUBROUTINE GWF1MNW2AL(isumz,lcmnw2,lcmnwn,lcmnwi,lcmnwc,mnwmax,
     +       in,iout,IMNWCB,MNWPRNT,NODTOT,IFREFM,NLAY,NMNWVL)
C     VERSION 20090405 GZH
c
c----- MNW2 
c     ******************************************************************
c     allocate array storage for MNW2 package
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER isumz,lcmnw2,lcmnwn,lcmnwi,mnwmax, 
     &        in,iout,IMNWCB,MNWPRNT,NODTOT,IFREFM,NLAY,
     &        lloc,istart,istop,isp,ispmnwn,ispmnwi,
     &        lcmnwc,ispmnwc,n,NAUX,NMNWVL,NMNW2,ntotnod
      REAL R
      CHARACTER*200 LINE
      COMMON /MNW2COM/MNWAUX(5)
      CHARACTER*16 MNWAUX
c
c1------identify package and initialize mnwmax
      write(iout,1) in
    1 format(/,1x,'MNW2 -- MULTI-NODE WELL PACKAGE, VERSION 2,',
     +' 04/05/2009.',/,4X,'INPUT READ FROM UNIT ',i3)
      NMNW2=0
      ntotnod=0
c
c2------read max number of mnw wells,
c2------unit or flag for cell-by-cell flow terms,
c2------output print flag, and aux variables
      CALL URDCOM(IN,IOUT,LINE)
      LLOC=1
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MNWMAX,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,IMNWCB,R,IOUT,IN)
      CALL URWORD(LINE,LLOC,ISTART,ISTOP,2,MNWPRNT,R,IOUT,IN)
c      READ(in,*) MNWMAX,IMNWCB,MNWPRNT
c
      write(iout,3) MNWMAX
    3 format(1h ,'MAXIMUM OF',i5,' ACTIVE MULTI-NODE WELLS AT ONE TIME')
      write(iout,*) 
      if(IMNWCB.gt.0) write(iout,9) IMNWCB
    9 format(1x, 'CELL-BY-CELL FLOWS WILL BE RECORDED ON UNIT', i3)
      if(IMNWCB.lt.0) write(iout,*) 'IMNWCB = ',IMNWCB
      if(IMNWCB.lt.0) write(iout,8)
    8 format(1x,'CELL-BY-CELL FLOWS WILL BE PRINTED WHEN ICBCFL NOT 0')
      write(iout,*) 'MNWPRNT = ',MNWPRNT
cdebug NoMoIter can be set here and used to force the solution to stop after
cdebug a certain amount of flow soplution iterations (used for debugging)
c      NoMoIter=9999
c      write(iout,7) NoMoIter
c    7 format(1x,'Flow rates will not be estimated after the',i4,'th',
c     +          ' iteration')
c
C3------READ AUXILIARY VARIABLES 
      NAUX=0
   10 CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
c      IF(LINE(ISTART:ISTOP).EQ.'CBCALLOCATE' .OR.
c     1   LINE(ISTART:ISTOP).EQ.'CBC') THEN
c         IMNWAL=1
c         WRITE(IOUT,11)
c   11    FORMAT(1X,'MEMORY IS ALLOCATED FOR CELL-BY-CELL BUDGET TERMS')
c         GO TO 10
      IF(LINE(ISTART:ISTOP).EQ.'AUXILIARY' .OR.
     1        LINE(ISTART:ISTOP).EQ.'AUX') THEN
         CALL URWORD(LINE,LLOC,ISTART,ISTOP,1,N,R,IOUT,IN)
         IF(NAUX.LT.5) THEN
            NAUX=NAUX+1
            MNWAUX(NAUX)=LINE(ISTART:ISTOP)
            WRITE(IOUT,12) MNWAUX(NAUX)
   12       FORMAT(1X,'AUXILIARY MNW2 VARIABLE: ',A)
         END IF
         GO TO 10
      END IF
      NMNWVL=32+NAUX
c3------set lcwel2 equal to location of well list in z array.
c
c4------add amount of space used by well list to isumz.
c     # of fields * max wells * max layers +50
      lcmnw2 = isumz    
c     NMNWVL is the number of records in MNW2
      isp=MNWMAX*(NMNWVL)
      isumz=isumz+isp
c
      lcmnwn=isumz
c approximate number of nodes= max mnw wells * number of layers, this works well
c if all are mostly vertical wells.  add 10*nlay+25 for extra room.  ispmnwn is 
c passed out to RP routine to check allocation while reading actual # nodes used
      NODTOT=(MNWMAX*NLAY)+(10*NLAY)+25
c  33 is the number of records in MNWNOD
      ispmnwn  = 33 * NODTOT   
      isumz = isumz+ispmnwn
c  11 is the number of records in MNWINT
      lcmnwi=isumz
      ispmnwi  = 11 * NODTOT   
      isumz = isumz+ispmnwi
c  27 is the hard-wired number of entries allowed in the Capacity Table (CapTable)
c  2 is the number of fields in the Capacity Table (CapTable): Lift (1) and Q (2)
      lcmnwc=isumz
      ispmnwc  = MNWMAX*27*2  
      isumz = isumz+ispmnwc
c
c5------print number of spaces in z array used by well package.
      write(iout,4) isp+ispmnwn+ispmnwi+ispmnwc
    4 format(1x,i6,' ELEMENTS IN Z ARRAY ARE USED FOR MNW2')
      write(iout,*) 
c7------return
      return
      end
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2RP(WELLID,MNW2,MNWNOD,MNWINT,
     +  BOTM,NBOTM,ibound,mnwmax,NODTOT,nrow,ncol,nlay,INTTOT,in,iout,
     +  igwtunit,nmnw2,small,hclose,kkper,MNWPRNT,ntotnod,CapTable,
     +  NMNWVL)
C     VERSION 20090405 GZH
c
c     ******************************************************************
c     read mnw2 locations, stress rates, conc, well char., and limits
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ntotnod,INTTOT,mnwid,mnwmax,in,nnodes,id,iout,
     &  nodnum,inode,il,ir,ic,intnum,nintvl,intnodes,iint,nodecount,
     &  IRlast,IClast,k,PUMPLAY,PUMPROW,PUMPCOL,firstnode,lastnode,
     &  ITMP,iread,igwtunit,ifound,iw,nodtot,nrow,ncol,nlay,nbotm,
     &  IBOUND,LBOTM,LAYCBD,nmnw2,PUMPLOC,Qlimit,i,j,kkper,MNWPRNT,
     &  PPFLAG,QCUT,nod,PUMPCAP,index,NMNWVL,NAUX,IAUX
      REAL BOTM,hclose
      DOUBLE PRECISION MNW2,MNWNOD,MNWINT,CapTable,CapMult
      DOUBLE PRECISION Rw,Rskin,Kskin,B,C,P,CWC,RwNode,RskinNode,
     & KskinNode,BNode,CNode,PNode,CWCNode,Ztop,Zbotm,Zbotmlast,
     & Zpump,Hlim,Qfrcmn,Qfrcmx,Qdes,Cprime,PP,
     & small,Qtemp,Hlift,LIFTq0,LIFTqdes,LIFTn,Qn,HWtol
      CHARACTER*20 WELLID,WELLNAME,LOSSTYPE

      dimension WELLID(mnwmax)
      dimension mnw2(NMNWVL,mnwmax),MNWNOD(33,NODTOT),MNWINT(11,NODTOT)
      dimension CapTable(mnwmax,27,2)
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),ibound(ncol,nrow,nlay)
     &          
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
      COMMON /MNW2COM/MNWAUX(5)
      CHARACTER*16 MNWAUX
C
c
c
c------------------------------------------------------------------
c     The 11 rows of the MNW2 array store:
c      Row #  = Description
c------------------------------------------------------------------
c         1   = IACTIV (0: inactive; 1: active)
c         2   = NNODES (number of nodes in this well)
c         3   = LOSSTYPE (0: none; 1:THIEM; 2: SKIN; 3: GENERAL; 4:SPEC. COND)
c         4   = NODNUM (number, in node list (MNWNOD), of first node of well)
c         5   = QDES (desired flow rate for this stress period)
c         6   = QLIMIT (pumpage constraint flag, QLIMIT>0 turns on constraint)
c         7   = HLIM (limiting water level for pumpage constraint)
c         8   = QCUT (pump cutoff flag: QCUT>0 limit by rate, QCUT<0 limit by 
c                     fraction of QDES)
c         9   = Qfrcmn (minimum rate for well to remain active)
c        10   = Qfrcmx (maximum rate to reactivate)
c        11   = Cprime (concentration of inflow)
c        12   = 
c------------------------------------------------------------------
c
c  If past first stress period, skip reading data sets1+2
      IF(KKPER.GT.1) GOTO 888
c     set defaults
      ntotnod=0
      INTTOT=0
      small=DBLE(hclose)
c     initialize
      WELLID=' '
      MNW2=0.0D0
      MNWNOD=0.0D0
      MNWINT=0.0D0
      write(iout,*) 
      write(iout,*) 'MNW2 Input:'
c0----read MNW info: wellid, nnodes, losstype, pumploc, Qlimit
      DO 1 MNWID=1,MNWMAX
c  initialize CapFlag2 at beginning of simulation
      mnw2(27,MNWID)=0
c
      write(iout,*) 
      write(iout,*) 
      write(iout,*) 'WELLID             NNODES    LOSSTYPE',
     &' PUMPLOC  Qlimit  PPFLAG PUMPCAP'
c     read Data Set 2a
        read(in,*) WELLID(MNWID),NNODES
c     read Data Set 2b
        read(in,*) LOSSTYPE,PUMPLOC,Qlimit,PPFLAG,PUMPCAP
c     convert wellid and losstype to uppercase
      call UPCASE(WELLID(MNWID))
      call UPCASE(LOSSTYPE)
c        write(iout,'(1x,A20,1x,I4,5x,A10,2x,I3,2(5x,I3))')
c     & WELLID(MNWID),INT(NNODES),LOSSTYPE,INT(PUMPLOC),Qlimit,PPFLAG
c     write output with Format depending on LOSSTYPE
        if(LOSSTYPE.EQ.'NONE') then
         write(iout,'(1x,A20,1x,I4,6x,A10,1x,I3,3(5x,I3))')
     & WELLID(MNWID),INT(NNODES),LOSSTYPE,INT(PUMPLOC),Qlimit,PPFLAG,
     & PUMPCAP
        elseif(LOSSTYPE.EQ.'THIEM') then
         write(iout,'(1x,A20,1x,I4,6x,A10,1x,I3,3(5x,I3))')
     & WELLID(MNWID),INT(NNODES),LOSSTYPE,INT(PUMPLOC),Qlimit,PPFLAG,
     & PUMPCAP
        elseif(LOSSTYPE.EQ.'SKIN') then
         write(iout,'(1x,A20,1x,I4,6x,A10,1x,I3,3(5x,I3))')
     & WELLID(MNWID),INT(NNODES),LOSSTYPE,INT(PUMPLOC),Qlimit,PPFLAG,
     & PUMPCAP
        elseif(LOSSTYPE.EQ.'GENERAL') then
         write(iout,'(1x,A20,1x,I4,5x,A10,2x,I3,3(5x,I3))')
     & WELLID(MNWID),INT(NNODES),LOSSTYPE,INT(PUMPLOC),Qlimit,PPFLAG,
     & PUMPCAP
        elseif(LOSSTYPE.EQ.'SPECIFYCWC') then
         write(iout,'(1x,A20,1x,I4,2x,A10,5x,I3,3(5x,I3))')
     & WELLID(MNWID),INT(NNODES),LOSSTYPE,INT(PUMPLOC),Qlimit,PPFLAG,
     & PUMPCAP
        else
	   write(iout,*) '***ERROR*** LOSSTYPE Not Recognized; stopping.'
	   STOP 'MNW2 ERROR - LOSSTYPE'
	  end if
c     check WELLID vs existing names
        do 2 ID=1,MNWID-1
	    if(WELLID(MNWID).EQ.WELLID(ID)) then
            write(iout,*) '***ERROR*** (MNW2) non-unique MNWID:,',
     *        WELLID(MNWID)
            STOP 'MNW2 ERROR - WELLID'
          end if                
    2   continue
c     set PUMPLOC and PUMPCAP
        MNW2(11,MNWID)=PUMPLOC
        MNW2(22,MNWID)=PUMPCAP
c     CapTable has max 27 entires, so PUMPCAP must not be > 25
        if(PUMPCAP.GT.25) then
          write(iout,*) '***ERROR*** PUMPCAP cannot be greater than 25'
          STOP 'MNW2 ERROR - PUMPCAP'
        end if
c     set IACTIV=0, defaulting well to inactive
        MNW2(1,MNWID)=0
c     set number of nodes
        MNW2(2,MNWID)=NNODES
c     define LOSSTYPE using integers in MNW2 array
        if(LOSSTYPE.EQ.'NONE') then
c     for none, NNODES must be 1
          if(NNODES.NE.1) then
            STOP 'MNW2 ERROR -  OPTION: NONE REQUIRES NNODES=1'
	    end if
          MNW2(3,MNWID)=0
        elseif(LOSSTYPE.EQ.'THIEM') then
          MNW2(3,MNWID)=1 
        elseif(LOSSTYPE.EQ.'SKIN') then
          MNW2(3,MNWID)=2
        elseif(LOSSTYPE.EQ.'GENERAL') then
          MNW2(3,MNWID)=3
        elseif(LOSSTYPE.EQ.'SPECIFYCWC') then
          MNW2(3,MNWID)=4	
	  end if
c     initialize QDES=0
        MNW2(5,MNWID)=0.0
c     set Qlimit flag.  Qlimit.ne.0 means limit or constraint is on.  Qlimit<0 means
c          read it every stress period
        MNW2(6,MNWID)=Qlimit
c     set PPFLAG flag.  PPFLAG>0 means calculate partial penetration effect.
        MNW2(19,MNWID)=PPFLAG
c     end read Data Set 2a and 2b
c
c     warning if LOSSTYPE=SPECIFYCWC and PPFLAG>0 (no PP correction done for this LOSSTYPE)
        if(LOSSTYPE.EQ.'SPECIFYCWC'.and.PPFLAG.GT.0) then
          write(iout,*)
          write(iout,*) '***WARNING*** Partial penetration not',
     & ' calculated for LOSSTYPE = SPECIFYCWC'
          write(iout,*)
        end if

c
c     read Data Set 2c, depending on LOSSTYPE (MNW2(3,MNWID)
        SELECT CASE (INT(MNW2(3,MNWID)))
          CASE (1)
            READ(in,*) Rw
c     don't allow Rw = 0
            if(Rw.eq.0.0) then
              write(iout,*) '***ERROR*** Rw=0.0; Rw read=',Rw
              STOP 'MNW2 ERROR - Rw'
            endif
          CASE (2)
            READ(in,*) Rw,Rskin,Kskin
            if(Rw.eq.0.0) then
              write(iout,*) '***ERROR*** Rw=0.0; Rw read=',Rw
              STOP 'MNW2 ERROR - Rw'
            endif
            if(Rskin.eq.0.0) then
              write(iout,*) '***ERROR*** Rskin=0.0; Rskin read=',Rskin
              STOP 'MNW2 ERROR'
            endif
            if(Kskin.eq.0.0) then
              write(iout,*) '***ERROR*** Kskin=0.0; Kskin read=',Kskin
              STOP 'MNW2 ERROR - Kskin'
            endif
          CASE (3)
            READ(in,*) Rw,B,C,P
            if(Rw.eq.0.0) then
              write(iout,*) '***ERROR*** Rw=0.0; Rw read=',Rw
              STOP 'MNW2 ERROR - Rw'
            endif
            if(P.gt.0.0.and.(P.lt.1.0.or.P.gt.3.5)) then
              write(iout,*) '***ERROR*** P=',P,' exceeds 1 <= P <=3.5'
              STOP 'MNW2 ERROR - P'
            endif
          CASE (4)
            READ(in,*) CWC            
        END SELECT
c     end read Data Set 2c
c
c     set HWflag (horizontal well flag) to 0 as default
c     can be set to 1 (true) if NNODES>0 and any R,C is different than first
            MNW2(21,MNWID)=0
c     read Data Set 2d, the list of nodes
c
c     if nnodes>0, read IL,IR,IC and LOSSTYPE variables set < 0
        IF(NNODES.GT.0) THEN
c     calculate first node number in nodal list= (total nodes so far) + 1
c     this will be used to access the information in the nodal array (MNWNOD)
          NODNUM=ntotnod+1
          MNW2(4,MNWID)=NODNUM
c     count nodes (will uses this to check vs. allocation)
          ntotnod=ntotnod+NNODES
          if(ntotnod.gt.nodtot) then
            write(iout,*)
            write(iout,*) 'MNW2 NODE ARRAY ALLOCATION INSUFFICIENT'
            STOP 'MNW2 ALLOCATION ERROR'
          end if
c     loop over nodes 
          DO 3 INODE=1,NNODES
c     If PPFLAG=0, don't read PP variable
           IF(MNW2(19,mnwid).eq.0) then           
c     access the LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, read IL,IR,IC only
              CASE (0)
                READ(in,*) IL,IR,IC
c
c     LOSSTYPE=THIEM, read IL,IR,IC,{Rw}
              CASE (1)
                IF(Rw.GT.0.0) THEN            
                  READ(in,*) IL,IR,IC
c     Rw at each node is the same if Rw>0.0  
                  RwNode=Rw         
                ELSE
c     If Rw<0, read in separate Rw for each node
                  READ(in,*) IL,IR,IC,RwNode
                END IF
c     Set Rw in node list, spot is 1st node (NODNUM) + current step (INODE) - 1
                MNWNOD(5,NODNUM+INODE-1)=RwNode             
c
c     LOSSTYPE=SKIN, read IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
                IF(Rw.GT.0.0) THEN            
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC
                      RwNode=Rw             
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,KskinNode
                      RwNode=Rw             
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,RskinNode
                      RwNode=Rw             
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,RskinNode,KskinNode
                      RwNode=Rw             
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,RwNode
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,RwNode,KskinNode
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,RwNode,RskinNode
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,RwNode,RskinNode,KskinNode
                    ENDIF
                  ENDIF
                END IF
c     Set vars in node list, spot is 1st node (NODNUM) + current step (INODE) - 1
                MNWNOD(5,NODNUM+INODE-1)=RwNode             
                MNWNOD(6,NODNUM+INODE-1)=RskinNode             
                MNWNOD(7,NODNUM+INODE-1)=KskinNode             
c
c     LOSSTYPE=GENERAL, read IL,IR,IC,{Rw B C P}
              CASE (3)
                IF(Rw.GT.0.0) THEN            
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,PNode
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,CNode
                        RwNode=Rw             
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,CNode,PNode
                        RwNode=Rw             
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,BNode
                        RwNode=Rw             
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,BNode,PNode
                        RwNode=Rw             
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,BNode,CNode
                        RwNode=Rw             
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,BNode,CNode,PNode
                        RwNode=Rw             
                      END IF
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,RwNode
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,PNode
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,Rwnode,CNode
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,CNode,PNode
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,Rwnode,BNode
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,BNode,PNode
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,Rwnode,BNode,CNode
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,BNode,CNode,PNode
                      END IF
                    ENDIF
                  ENDIF
                END IF
c     Set vars in node list, spot is 1st node (NODNUM) + current step (INODE) - 1
                MNWNOD(5,NODNUM+INODE-1)=RwNode             
                MNWNOD(8,NODNUM+INODE-1)=BNode             
                MNWNOD(9,NODNUM+INODE-1)=CNode            
                MNWNOD(10,NODNUM+INODE-1)=PNode             
c
c     LOSSTYPE=SPECIFYcwc, read IL,IR,IC,{CWC}
              CASE (4)
                IF(CWC.GT.0.0) THEN            
                  READ(in,*) IL,IR,IC
                  CWCNode=CWC         
                ELSE
                  READ(in,*) IL,IR,IC,CWCNode
                END IF
                MNWNOD(11,NODNUM+INODE-1)=CWCNode             
            END SELECT
c    ELSE if PPFLAG NE 0, read PP flag
           ELSE
c     access the LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, read IL,IR,IC only
              CASE (0)
                READ(in,*) IL,IR,IC,PP
c
c     LOSSTYPE=THIEM, read IL,IR,IC,{Rw}
              CASE (1)
                IF(Rw.GT.0.0) THEN           
                  READ(in,*) IL,IR,IC,PP
c     Rw at each node is the same if Rw>0.0  
                  RwNode=Rw         
                ELSE
c     If Rw<0, read in separate Rw for each node
                  READ(in,*) IL,IR,IC,RwNode,PP
                END IF
c     Set Rw in node list, spot is 1st node (NODNUM) + current step (INODE) - 1
                MNWNOD(5,NODNUM+INODE-1)=RwNode             
c
c     LOSSTYPE=SKIN, read IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
                IF(Rw.GT.0.0) THEN            
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,PP
                      RwNode=Rw             
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,KskinNode,PP
                      RwNode=Rw             
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,RskinNode,PP
                      RwNode=Rw             
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,RskinNode,KskinNode,PP
                      RwNode=Rw             
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,RwNode,PP
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,RwNode,KskinNode,PP
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) IL,IR,IC,RwNode,RskinNode,PP
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) IL,IR,IC,RwNode,RskinNode,KskinNode,PP
                    ENDIF
                  ENDIF
                END IF
c     Set vars in node list, spot is 1st node (NODNUM) + current step (INODE) - 1
                MNWNOD(5,NODNUM+INODE-1)=RwNode             
                MNWNOD(6,NODNUM+INODE-1)=RskinNode             
                MNWNOD(7,NODNUM+INODE-1)=KskinNode             
c
c     LOSSTYPE=GENERAL, read IL,IR,IC,{Rw B C P}
              CASE (3)
                IF(Rw.GT.0.0) THEN            
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,PP
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,PNode,PP
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,CNode,PP
                        RwNode=Rw             
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,CNode,PNode,PP
                        RwNode=Rw             
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,BNode,PP
                        RwNode=Rw             
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,BNode,PNode,PP
                        RwNode=Rw             
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,BNode,CNode,PP
                        RwNode=Rw             
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,BNode,CNode,PNode,PP
                        RwNode=Rw             
                      END IF
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,RwNode,PP
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,PNode,PP
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,Rwnode,CNode,PP
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,CNode,PNode,PP
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,Rwnode,BNode,PP
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,BNode,PNode,PP
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) IL,IR,IC,Rwnode,BNode,CNode,PP
                        PNode=P
                      ELSE
                        READ(in,*) IL,IR,IC,Rwnode,BNode,CNode,PNode,PP
                      END IF
                    ENDIF
                  ENDIF
                END IF
c     Set vars in node list, spot is 1st node (NODNUM) + current step (INODE) - 1
                MNWNOD(5,NODNUM+INODE-1)=RwNode             
                MNWNOD(8,NODNUM+INODE-1)=BNode             
                MNWNOD(9,NODNUM+INODE-1)=CNode            
                MNWNOD(10,NODNUM+INODE-1)=PNode             
c
c     LOSSTYPE=SPECIFYcwc, read IL,IR,IC,{CWC}
              CASE (4)
                IF(CWC.GT.0.0) THEN            
                  READ(in,*) IL,IR,IC,PP
                  CWCNode=CWC         
                ELSE
                  READ(in,*) IL,IR,IC,CWCNode,PP
                END IF
                MNWNOD(11,NODNUM+INODE-1)=CWCNode             
            END SELECT            
           END IF
c     save node location, set Qdes=0.0, set PP, flag ZPD
            MNWNOD(1,NODNUM+INODE-1)=IL             
            MNWNOD(2,NODNUM+INODE-1)=IR            
            MNWNOD(3,NODNUM+INODE-1)=IC           
            MNWNOD(4,NODNUM+INODE-1)=0.0           
            MNWNOD(19,NODNUM+INODE-1)=PP          
c     save IR and IC to check vs the subsequent nodes for vert/horiz
            IF(INODE.EQ.1) THEN
              IRlast=IR
              IClast=IC
            ELSE
c     if any node is a different R,C, this is a horizontal well
              IF((IR.NE.IRlast).OR.(IC.NE.IClast)) THEN
c       set HWflag to true
                MNW2(21,MNWID)=1
              END IF
            END IF

c     if partial penetration ne 0, set ZPD to 1d30.  this will act as a flag
c     until ZDP (and ZPL) are set for the well)
            if(pp.ne.0.d0) MNWNOD(20,NODNUM+INODE-1)=1d30         
    3     CONTINUE
c
c       end nnodes>0 read statements
c
c     if nnodes<0, read in Ztop and Zbot which define intervals
        ELSE
c
c     calculate first interval number in interval list= (total ints so far) + 1
c     this will be used to access the information in the interval array (MNWINT)
          INTNUM=INTTOT+1
          MNW2(13,MNWID)=INTNUM
c     the abs value of NNODES represents the number of intervals to de defined
          NINTVL=ABS(NNODES)
c     count intervals to check vs. allocation
          INTTOT=INTTOT+NINTVL
c     initialize interval node counter
          intnodes=0
c     loop over the intervals in this well
          DO 4 IINT=1,NINTVL
c
c     access the LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, read Ztop,Zbotm,IR,IC only
              CASE (0)
                READ(in,*) Ztop,Zbotm,IR,IC
c
c     LOSSTYPE=THIEM, read Ztop,Zbotm,IR,IC,{Rw}
              CASE (1)
                IF(Rw.GT.0.0) THEN            
                  READ(in,*) Ztop,Zbotm,IR,IC
c     Rw at each node is the same if Rw>0.0  
                  RwNode=Rw         
                ELSE
c     If Rw<0, read in separate Rw for each node
                  READ(in,*) Ztop,Zbotm,IR,IC,RwNode
                END IF
c     Set Rw in interval list, spot is 1st int (INTNUM) + current step (IINT) - 1
                MNWINT(5,INTNUM+IINT-1)=RwNode             
c
c     LOSSTYPE=SKIN, read Ztop,Zbotm,IR,IC,{Rw Rskin Kskin}
              CASE (2)
                IF(Rw.GT.0.0) THEN            
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) Ztop,Zbotm,IR,IC
                      RwNode=Rw             
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) Ztop,Zbotm,IR,IC,KskinNode
                      RwNode=Rw             
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) Ztop,Zbotm,IR,IC,RskinNode
                      RwNode=Rw             
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) Ztop,Zbotm,IR,IC,RskinNode,KskinNode
                      RwNode=Rw             
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(Rskin.GT.0.0) THEN
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) Ztop,Zbotm,IR,IC,RwNode
                      RskinNode=Rskin            
                      KskinNode=Kskin
                    ELSE             
                      READ(in,*) Ztop,Zbotm,IR,IC,RwNode,KskinNode
                      RskinNode=Rskin            
                    ENDIF
c                 else Rskin<0
                  ELSE
                    IF(Kskin.GT.0.0) THEN
                      READ(in,*) Ztop,Zbotm,IR,IC,RwNode,RskinNode
                        KskinNode=Kskin
                    ELSE             
                      READ(in,*) Ztop,Zbotm,IR,IC,RwNode,RskinNode,
     &                  KskinNode
                    ENDIF
                  ENDIF
                END IF
c     Set vars for interval
                MNWINT(5,INTNUM+IINT-1)=RwNode             
                MNWINT(6,INTNUM+IINT-1)=RskinNode             
                MNWINT(7,INTNUM+IINT-1)=KskinNode             
c
c     LOSSTYPE=GENERAL, read Ztop,Zbotm,IR,IC,{Rw B C P}
              CASE (3)
                IF(Rw.GT.0.0) THEN            
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,PNode
                        RwNode=Rw             
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,CNode
                        RwNode=Rw             
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,CNode,PNode
                        RwNode=Rw             
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,BNode
                        RwNode=Rw             
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,BNode,PNode
                        RwNode=Rw             
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,BNode,CNode
                        RwNode=Rw             
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,BNode,CNode,PNode
                        RwNode=Rw             
                      END IF
                    ENDIF
                  ENDIF
c               else Rw<0
                ELSE
                  IF(B.GE.0.0) THEN
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,RwNode
                        BNode=B            
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,PNode
                        BNode=B            
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,CNode
                        BNode=B            
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,CNode,PNode
                        BNode=B            
                      END IF
                    ENDIF
c                 else B<0
                  ELSE
                    IF(C.GE.0.0) THEN
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,BNode
                        CNode=C
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,BNode,PNode
                        CNode=C
                      END IF
c                   else C<0
                    ELSE             
                      IF(P.GE.0.0) THEN
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,BNode,CNode
                        PNode=P
                      ELSE
                        READ(in,*) Ztop,Zbotm,IR,IC,Rwnode,BNode,CNode,
     &                    PNode
                      END IF
                    ENDIF
                  ENDIF
                END IF
c     Set vars for interval
                MNWINT(5,INTNUM+IINT-1)=RwNode             
                MNWINT(8,INTNUM+IINT-1)=BNode             
                MNWINT(9,INTNUM+IINT-1)=CNode            
                MNWINT(10,INTNUM+IINT-1)=PNode             
c
c     LOSSTYPE=SPECIFYcwc, read Ztop,Zbotm,IR,IC,{CWC}
              CASE (4)
                IF(CWC.GT.0.0) THEN            
                  READ(in,*) Ztop,Zbotm,IR,IC
                  CWCNode=CWC         
                ELSE
                  READ(in,*) Ztop,Zbotm,IR,IC,CWCNode
                END IF
c     Set var for interval
                MNWINT(11,INTNUM+IINT-1)=CWCNode             
            END SELECT
c     Set vars for interval
            MNWINT(1,INTNUM+IINT-1)=Ztop            
            MNWINT(2,INTNUM+IINT-1)=Zbotm            
            MNWINT(3,INTNUM+IINT-1)=IR            
            MNWINT(4,INTNUM+IINT-1)=IC            
c     Do error-checking on interval elevations
c     If beyond the first interval, do check on elevation
            IF(IINT.GT.1) THEN
              IF(Ztop.GT.Zbotmlast) THEN
                WRITE(iout,*) '***ERROR*** Top of interval is ',
     &           'above bottom of last interval'                
                WRITE(iout,*) 'Well: ',WELLID(MNWID)
                STOP 'MNW2 ERROR - Intervals'
	        END IF
	      END IF
c     Save bottom of last interval for above check

            Zbotmlast=Zbotm
c
c     create nodes from interval information
c     MNWNOD will access MNWINT by pointing to the first and last interval
c       that intersects the node
c
c     set node counter for this interval
            nodecount=0
c     save IR and IC to check vs the subsequent nodes
            IF(IINT.EQ.1) THEN
              IRlast=IR
              IClast=IC
            ELSE
              IF((IR.NE.IRlast).OR.(IC.NE.IClast)) THEN
                write(iout,*) '***ERROR*** Row,Col must be constant in',
     &            'a vertical well (NNODES < 0)'
                STOP 'MNW2 ERROR - Vertical'
              END IF
            END IF
c
c     find first layer that the top of this interval penetrates, set = IL
            K=1
c     botm(...k) points to the bottom of layer K
            DO WHILE (Ztop.le.BOTM(IC,IR,LBOTM(K)))
              K=K+1
            END DO
            IF(K.LE.NLAY) then
              IL=K
            ELSE
              write(iout,*) '***ERROR*** MNW: ',
     &                      'Ztop below bottom of model'
              STOP 'MNW2 ERROR - Ztop'
            END IF
c
c     now that we have coordinates, start creating cells
c
c     if we haven't create any cells, create the first
            IF(intnodes.eq.0) THEN
              IF(IBOUND(IC,IR,IL).NE.0) THEN
c     calculate first node number in nodal list= (total nodes so far) +1
c     this will be used to access the information in the nodal array (MNWNOD)
                NODNUM=ntotnod+1
c     increase total node count 
                ntotnod=ntotnod+1          
          if(ntotnod.gt.nodtot) then
            write(iout,*)
            write(iout,*) 'MNW2 NODE ARRAY ALLOCATION INSUFFICIENT'
            STOP 'MNW2 ALLOCATION ERROR'
          end if
c     increase count of nodes in this interval
                nodecount=nodecount+1   
c     increase count of nodes in this well
                intnodes=intnodes+1   
c     mark this as the first node in the node list for this well
                if(intnodes.eq.1) MNW2(4,MNWID)=NODNUM
c     create node in node list from interval coordinates
                MNWNOD(1,NODNUM)=IL             
                MNWNOD(2,NODNUM)=IR            
                MNWNOD(3,NODNUM)=IC           
                MNWNOD(4,NODNUM)=0.0           
c     set first and last intervals in this node= interval #1
c     INTNUM is location in int list of 1st interval for this well
c     IINT is loop counter for intervals
                MNWNOD(12,NODNUM)=INTNUM+IINT-1
                MNWNOD(13,NODNUM)=INTNUM+IINT-1
c     if partial penetration ne 0, set ZPD to 1d30.  this will act as a flag
c     until ZDP (and ZPL) are set for the well)
                MNWNOD(20,NODNUM)=1d30         
              ELSE
                WRITE(iout,*) '***ERROR*** MNW2 screen in no-flow node'
                STOP 'MNW2 - screen in no-flow node' 
              END IF
c     if a node has been created (nodecount>0), then check to see if this interval
c     is still in that node (still NODNUM)
            ELSE
              IF(MNWNOD(1,NODNUM).EQ.IL) THEN
c     if interval is still in previous node, re-set "last int" for that node
c     do not increase nodecount, as it is in same node
                MNWNOD(13,NODNUM)=INTNUM+IINT-1
c     if top of this interval is in a new node, create a node in that layer
              ELSE
                IF(IBOUND(IC,IR,IL).NE.0) THEN
                  NODNUM=NODNUM+1   
c     increase total node count 
                  ntotnod=ntotnod+1          
          if(ntotnod.gt.nodtot) then
            write(iout,*)
            write(iout,*) 'MNW2 NODE ARRAY ALLOCATION INSUFFICIENT'
            STOP 'MNW2 ALLOCATION ERROR'
          end if
c     increase count of nodes in this interval
                  nodecount=nodecount+1   
c     increase count of nodes in this well
                  intnodes=intnodes+1   
                  MNWNOD(1,NODNUM)=IL             
                  MNWNOD(2,NODNUM)=IR            
                  MNWNOD(3,NODNUM)=IC           
                  MNWNOD(4,NODNUM)=0.0           
c     set first and last intervals in this node= interval #1
c     INTNUM is location in int list of 1st interval for this well
c     IINT is loop counter for intervals
                  MNWNOD(12,NODNUM)=INTNUM+IINT-1
                  MNWNOD(13,NODNUM)=INTNUM+IINT-1
c     if partial penetration ne 0, set ZPD to 1d30.  this will act as a flag
c     until ZDP (and ZPL) are set for the well)
                  MNWNOD(20,NODNUM)=1d30         
                ELSE
                 WRITE(iout,*) '***ERROR*** MNW2 screen in no-flow node'
                 STOP 'MNW2 - screen in no-flow node' 
                END IF
	        END IF
            END IF
c
c     check nodecount here before looping over lower layers, possibility that the
c     interval defines nothing is then caught
            IF(nodecount.gt.0) THEN
c     add more nodes if the bottom of the interval penetrates below this layer,
c     as long as it isn't the last layer
              K=IL
              DO WHILE(Zbotm.lt.BOTM(IC,IR,LBOTM(K)).and.
     &               ((K+1).LE.NLAY))                  
                K=K+1
                IL=K
			    NODNUM=NODNUM+1
c     increase total node count 
                ntotnod=ntotnod+1          
          if(ntotnod.gt.nodtot) then
            write(iout,*)
            write(iout,*) 'MNW2 NODE ARRAY ALLOCATION INSUFFICIENT'
            STOP 'MNW2 ALLOCATION ERROR'
          end if
c     increase count of nodes in this interval
                nodecount=nodecount+1   
c     increase count of nodes in this well
                intnodes=intnodes+1   
                MNWNOD(1,NODNUM)=IL             
                MNWNOD(2,NODNUM)=IR            
                MNWNOD(3,NODNUM)=IC           
                MNWNOD(4,NODNUM)=0.0           
c     set first and last intervals in this node= interval number
                MNWNOD(12,NODNUM)=INTNUM+IINT-1
                MNWNOD(13,NODNUM)=INTNUM+IINT-1
c     if partial penetration ne 0, set ZPD (MNWNOD(20) to 1d30.  this will act as a flag
c     until ZDP (and ZPL) are set for the well
                MNWNOD(20,NODNUM)=1d30         
              END DO
	      END IF
c

    4     CONTINUE
c         end loop over intervals
c     reset NUMNODES to -(number of nodes) instead of -(# of intervals)
            MNW2(2,MNWID)=-1*intnodes
c
c     Print interval information
          if(MNWPRNT.gt.0) then
            write(iout,*)
c     write info line
            write(iout,'(100A)') ' NNODES < 0: well defined using open',
     &' intervals as described below'
c     write header depending on LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, read Ztop,Zbotm,IR,IC only
              CASE (0)
            write(iout,'(100A)') ' Interval      Ztop       Zbotm    ',
     &'  Row  Col'
c     LOSSTYPE=THIEM, read Ztop,Zbotm,IR,IC,{Rw}
              CASE (1)
            write(iout,'(100A)') ' Interval      Ztop       Zbotm    ',
     &'  Row  Col      Rw     '

c     LOSSTYPE=SKIN, read Ztop,Zbotm,IR,IC,{Rw Rskin Kskin}
              CASE (2)
            write(iout,'(100A)') ' Interval      Ztop       Zbotm    ',
     &'  Row  Col      Rw     Rskin    ',
     &' Kskin '

c     LOSSTYPE=GENERAL, read Ztop,Zbotm,IR,IC,{Rw B C P}
              CASE (3)
            write(iout,'(100A)') ' Interval      Ztop       Zbotm    ',
     &'  Row  Col      Rw     B         C         P  '

c     LOSSTYPE=SPECIFYcwc, read Ztop,Zbotm,IR,IC,{CWC}
              CASE (4)
            write(iout,'(100A)') ' Interval      Ztop       Zbotm    ',
     &'  Row  Col      spec.CWC'
            END SELECT
c
c     get first interval in this well
            INTNUM=MNW2(13,MNWID)
c     write data depending on LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, write Ztop,Zbotm,IR,IC only
              CASE (0)
            do IINT=INTNUM,(INTNUM+NINTVL-1)
            write(iout,'(1x,I4,6x,1P2G12.5,1x,I4,1x,I4)')
     &IINT-intnum+1,(MNWINT(j,iint),j=1,2),(INT(MNWINT(j,iint)),j=3,4)
            end do
c     LOSSTYPE=THIEM, write Ztop,Zbotm,IR,IC,{Rw}
              CASE (1)
            do IINT=INTNUM,(INTNUM+NINTVL-1)
            write(iout,'(1x,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P1G10.4)')
     &IINT-intnum+1,(MNWINT(j,iint),j=1,2),(INT(MNWINT(j,iint)),j=3,4),
     &(MNWINT(5,iint))
            end do
c     LOSSTYPE=SKIN, write Ztop,Zbotm,IR,IC,{Rw Rskin Kskin}
              CASE (2)
            do IINT=INTNUM,(INTNUM+NINTVL-1)
            write(iout,'(1x,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P3G10.4)')
     &IINT-intnum+1,(MNWINT(j,iint),j=1,2),(INT(MNWINT(j,iint)),j=3,4),
     &(MNWINT(j,iint),j=5,7)
            end do
c     LOSSTYPE=GENERAL, write Ztop,Zbotm,IR,IC,{Rw B C P}
              CASE (3)
            do IINT=INTNUM,(INTNUM+NINTVL-1)
            write(iout,'(1x,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P4G10.4)')
     &IINT-intnum+1,(MNWINT(j,iint),j=1,2),(INT(MNWINT(j,iint)),j=3,4),
     &(MNWINT(5,iint)),(MNWINT(j,iint),j=8,10)
            end do
c     LOSSTYPE=SPECIFYcwc, write Ztop,Zbotm,IR,IC,{CWC}
              CASE (4)
            do IINT=INTNUM,(INTNUM+NINTVL-1)
            write(iout,'(1x,I4,6x,1P2G12.5,1x,I4,1x,I4,2x,1P1G10.4)')
     &IINT-intnum+1,(MNWINT(j,iint),j=1,2),(INT(MNWINT(j,iint)),j=3,4),
     &(MNWINT(11,iint))
            end do
            END SELECT
          end if
c
        END IF
c       end read Data Set 2d
c
c     loop over wells to print out Node Info
c
      write(iout,*)
c
      if(MNWPRNT.gt.0) then
        firstnode=MNW2(4,MNWID)
        lastnode=MNW2(4,MNWID)+ABS(MNW2(2,MNWID))-1
c     write info line
        if(nnodes.lt.0) then
c     for MNWs defined by intervals
         write(iout,'(100A)') ' The following',
     &' nodes were assigned to this well based on above open interval',
     &' information'
         write(iout,'(100A)') ' Node  Lay  Row  Col '
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
c  if more than one interval made up this node, write composite 
          if(MNWNOD(12,INODE).ne.MNWNOD(13,INODE)) then

           write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,9A)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),
     & 'COMPOSITE'

          else
           write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3)       
          end if
         end do
c
        else
c     for MNWs defined by nodes
c
      if(PPFLAG.EQ.0) then
c     write header depending on LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, read IL,IR,IC only
              CASE (0)
      write(iout,'(100A)') ' Node  Lay  Row  Col'
c     LOSSTYPE=THIEM, read IL,IR,IC,{Rw}
              CASE (1)
      write(iout,'(100A)') ' Node  Lay  Row  Col      Rw     '
c     LOSSTYPE=SKIN, read IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
      write(iout,'(100A)') ' Node  Lay  Row  Col      Rw     Rskin    ',
     &' Kskin'
c     LOSSTYPE=GENERAL, read IL,IR,IC,{Rw B C P}
              CASE (3)
      write(iout,'(100A)') ' Node  Lay  Row  Col      Rw     B
     &         C          P  '
c     LOSSTYPE=SPECIFYcwc, read IL,IR,IC,{CWC}
              CASE (4)
      write(iout,'(100A)') ' Node  Lay  Row  Col    spec.CWC'
            END SELECT
c
c     write data depending on LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, write IL,IR,IC only
              CASE (0)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3)       
         end do
c     LOSSTYPE=THIEM, write IL,IR,IC,{Rw}
              CASE (1)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P1G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),(MNWNOD(5,INODE))       
         end do
c     LOSSTYPE=SKIN, write IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P3G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),(MNWNOD(j,INODE),j=5,7)       
         end do
c     LOSSTYPE=GENERAL, write IL,IR,IC,{Rw B C P}
              CASE (3)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P4G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),(MNWNOD(5,INODE)),
     &(MNWNOD(j,INODE),j=8,10)       
         end do
c     LOSSTYPE=SPECIFYcwc, write IL,IR,IC,{CWC}
              CASE (4)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P1G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),(MNWNOD(11,INODE))       
         end do
            END SELECT
c     If PPFLAG>0 print PP input
      else
c     write header depending on LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, read IL,IR,IC only
              CASE (0)
      write(iout,'(100A)') ' Node  Lay  Row  Col        PP'
c     LOSSTYPE=THIEM, read IL,IR,IC,{Rw}
              CASE (1)
      write(iout,'(100A)') ' Node  Lay  Row  Col      Rw        PP'
c     LOSSTYPE=SKIN, read IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
      write(iout,'(100A)') ' Node  Lay  Row  Col      Rw     Rskin    ',
     &' Kskin        PP'
c     LOSSTYPE=GENERAL, read IL,IR,IC,{Rw B C P}
              CASE (3)
      write(iout,'(100A)') ' Node  Lay  Row  Col      Rw     B
     &         C          P        PP'
c     LOSSTYPE=SPECIFYcwc, read IL,IR,IC,{CWC}
              CASE (4)
      write(iout,'(100A)') ' Node  Lay  Row  Col    spec.CWC        PP'
            END SELECT
c
c     write data depending on LOSSTYPE
            SELECT CASE (INT(MNW2(3,MNWID)))
c     LOSSTYPE=NONE, write IL,IR,IC only
              CASE (0)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,G10.3)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),MNWNOD(19,INODE)       
         end do
c     LOSSTYPE=THIEM, write IL,IR,IC,{Rw}
              CASE (1)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P2G10.3)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),MNWNOD(5,INODE),
     & MNWNOD(19,INODE)       
         end do
c     LOSSTYPE=SKIN, write IL,IR,IC,{Rw Rskin Kskin}
              CASE (2)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P4G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),(MNWNOD(j,INODE),j=5,7),
     & MNWNOD(19,INODE)      
         end do
c     LOSSTYPE=GENERAL, write IL,IR,IC,{Rw B C P}
              CASE (3)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P5G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),(MNWNOD(5,INODE)),
     &(MNWNOD(j,INODE),j=8,10),MNWNOD(19,INODE)       
         end do
c     LOSSTYPE=SPECIFYcwc, write IL,IR,IC,{CWC}
              CASE (4)
         do INODE=firstnode,lastnode
          nod=INODE-firstnode+1
          write(iout,'(1x,I4,1x,I4,1x,I4,1x,I4,2x,1P2G10.4)')
     &nod,(INT(MNWNOD(i,INODE)),i=1,3),MNWNOD(11,INODE),
     &MNWNOD(19,INODE)      
         end do
            END SELECT
      end if      
        end if
      end if
c
c   check well nodes in grid
c   Loop over nodes in well
      firstnode=MNW2(4,MNWID)
      lastnode=MNW2(4,MNWID)+ABS(MNW2(2,MNWID))-1
      do INODE=firstnode,lastnode
        il=MNWNOD(1,INODE)              
        ir=MNWNOD(2,INODE)              
        ic=MNWNOD(3,INODE) 

        if(il.lt.1.or.il.gt.nlay.or.
     &     ir.lt.1.or.ir.gt.nrow.or.
     &     ic.lt.1.or.ic.gt.ncol) then
          write(iout,*) 
          write(iout,*) 'MNW2 Node not in grid; Layer, Row, Col='
          write(iout,*) il,ir,ic
          STOP 'MNW2 ERROR - Grid'
        end if
      end do 
c
c     read Data Set 2e, PUMPLOC
c
        IF(PUMPLOC.GT.0.0) THEN
c     if PUMPLOC>0, read PUMPLAY,PUMPROW,PUMPCOL 
          READ(in,*) PUMPLAY,PUMPROW,PUMPCOL    
          MNW2(14,MNWID)=PUMPLAY
          MNW2(15,MNWID)=PUMPROW
          MNW2(16,MNWID)=PUMPCOL

        ELSEIF(PUMPLOC.LT.0) THEN
c     if PUMPLOC<0, read Zpump and calulate PUMPLAY,PUMPROW,PUMPCOL
          READ(in,*) Zpump
c         loop over nodes in this well
          firstnode=MNW2(4,MNWID)
          lastnode=MNW2(4,MNWID)+ABS(MNW2(2,MNWID))-1
          ifound=0
          DO INODE=firstnode,lastnode
            IL=MNWNOD(1,INODE)
            IF(Zpump.LT.BOTM(IC,IR,LBOTM(IL)-1).AND.
     &         Zpump.GT.BOTM(IC,IR,LBOTM(IL))) THEN
              IR=MNWNOD(2,INODE)
              IC=MNWNOD(3,INODE)
              MNW2(14,MNWID)=IL
              MNW2(15,MNWID)=IR
              MNW2(16,MNWID)=IC
              ifound=1
            END IF
          END DO
c     if PUMPLOC not in a node, assume it is at top and print warning
          if(ifound.eq.0) then
            write(iout,*) '***WARNING*** Pump location specified but 
     & not found within a MNW2 node'
            write(iout,*) ' Pump assumed to be at top node'
            MNW2(11,MNWID)=0
          end if
	  END IF
c
c     end read Data Set 2f
c
        IF(Qlimit.GT.0.0) THEN
          READ(in,*) Hlim,QCUT
          BACKSPACE in
          IF(QCUT.NE.0) THEN
            READ(in,*) Hlim,QCUT,Qfrcmn,Qfrcmx
          ELSE
            READ(in,*) Hlim,QCUT
          END IF
c     write info line
            write(iout,*) 
            write(iout,'(100A)') ' Qlimit > 0 : this well will
     & be constrained'
            write(iout,'(A,1PG13.5)') '     Hlim = ',Hlim 
            write(iout,1111) QCUT
 1111 FORMAT('     QCUT = ',I4)
           if(QCUT.lt.0) then 
           write(iout,'(A,1PG13.5,A)') '   Qfrcmn = ',Qfrcmn
           write(iout,'(A,1PG13.5,A)') '   Qfrcmx = ',Qfrcmx
           elseif(QCUT.gt.0) then
            write(iout,'(A,1PG13.5,A)') '   Qfrcmn = ',Qfrcmn, ' L**3/T'
            write(iout,'(A,1PG13.5,A)') '   Qfrcmx = ',Qfrcmx, ' L**3/T'
           end if
c     process min and max Q's based on QCUT
          MNW2(7,MNWID)=Hlim
          MNW2(8,MNWID)=QCUT
          MNW2(9,MNWID)=Qfrcmn
          MNW2(10,MNWID)=Qfrcmx
c     error check Qfrcmn
          if(QCUT.GT.0) then
           if(Qfrcmn.gt.1.d0) then
            write(iout,*) '***ERROR*** Qfrcmn > 1 is out of range'
            STOP 'MNW2 ERROR - Qfrcmn'
           end if
c     error check Qfrcmx
           if(Qfrcmx.gt.1.d0) then
            write(iout,*) '***ERROR*** Qfrcmx > 1 is out of range'
            STOP 'MNW2 ERROR - Qfrcmx'
           end if
	    end if
	  END IF
c
c     end read Data Set 2f
c
c     read Data Set 2g, PUMPCAP data
c
        IF(PUMPCAP.GT.0.0) THEN
c     if PUMPCAP>0, read Hlift, LIFTq0, LIFTqdes
          READ(in,*) Hlift,LIFTq0,LIFTqdes,HWtol   
          mnw2(23,MNWID)=Hlift
          mnw2(28,MNWID)=HWtol
c    CapTable(WELLID,index,type) where
c      index counts the number of values (PUMPCAP+2)
c      type 1 = Lift values
c      type 2 = Q values
          CapTable(MNWID,1,1)=LIFTq0
          CapTable(MNWID,1,2)=0.d0
          CapTable(MNWID,PUMPCAP+2,1)=abs(LIFTqdes)
c  when we know qdes, set this
c          CapTable(MNWID,PUMPCAP+1,2)=qdes
          DO 382 index=2,PUMPCAP+1
            READ(in,*) Liftn,Qn              
            CapTable(MNWID,index,1)=LIFTn
            CapTable(MNWID,index,2)=abs(Qn)
 382      CONTINUE
	  END IF
c     end read Data Set 2g
c
c     check consistency of CapTable
        DO 383 index=2,PUMPCAP
          if(CapTable(MNWID,index,1).GT.CapTable(MNWID,index-1,1)) then
            write(iout,*) '***ERROR*** Lift values in capacity table 
     &must be in descending order'
            STOP 'MNW2 ERROR - CapTable'
          end if
          if(CapTable(MNWID,index,2).LT.CapTable(MNWID,index-1,2)) then
            write(iout,*) '***ERROR*** Q values in capacity table 
     &must be in ascending order'
            STOP 'MNW2 ERROR - CapTable'
          end if
 383    CONTINUE

c
c     end loop over MNWMAX
    1 CONTINUE
c
c     read Data Set 3, ITMP (number of wells or flag saying reuse well data)
c
c  Skip to here unless in first SP
  888 CONTINUE
      READ(in,*) ITMP
      write(iout,*) 
c     return if there are no wells to read ........
      if( itmp .ge. 0 ) then
c     reset # of wells and active well flag
         nmnw2=0
         do iw=1,MNWMAX
           MNW2(1,iw)=0
	   end do
         if(itmp.eq.0) return
      end if    
      if( itmp.lt.0 ) then
c        if itmp less than zero, reuse data. print message and return.
        write(iout,6)
    6   format(1h0,'REUSING MNW2 INFORMATION FROM LAST STRESS PERIOD')
        return
      else
c  If itmp > 0, read ITMP wells
c
c
       if(itmp.gt.1) then
        write(iout,*) 'MNW2: ',ITMP, ' active wells in stress period ',
     & kkper
       else
        write(iout,*) 'MNW2: ',ITMP, ' active well in stress period ',
     & kkper
       end if
        write(iout,*) 
        do iread=1,ITMP
c  read data set 4a
c  read WELLNAME and then backspace to check for PUMPCAP
         read(in,*) WELLNAME
         backspace(in)
         call UPCASE(WELLNAME)
c  check for existence of well
         ifound=0
         MNWID=0
         do iw=1,MNWMAX
 	      if(WELLID(iw).EQ.WELLNAME) then
              ifound=1
c             set IACTIV=1 to turn well on for this SP
              MNW2(1,iw)=1
              nmnw2=nmnw2+1
              MNWID=iw
	      end if
	   end do
         if (ifound.eq.0) then
            write(iout,*) '***ERROR*** Well name not found in list'
            STOP 'MNW2 ERROR - WELLID'
         end if
c   continue reading data set 4a        
         PUMPCAP=MNW2(22,MNWID)
         NAUX=NMNWVL-32
         if(PUMPCAP.EQ.0) then
          if(igwtunit.le.0) then
            read(in,*) WELLNAME,Qdes,(MNW2(32+IAUX,MNWID),IAUX=1,NAUX)
          else
            read(in,*) WELLNAME,Qdes,Cprime,
     &                 (MNW2(32+IAUX,MNWID),IAUX=1,NAUX)
	    endif
	   else
          if(igwtunit.le.0) then
            read(in,*) WELLNAME,Qdes,CapMult,
     &                 (MNW2(32+IAUX,MNWID),IAUX=1,NAUX)
         else
            read(in,*) WELLNAME,Qdes,CapMult,Cprime,
     &                 (MNW2(32+IAUX,MNWID),IAUX=1,NAUX)
	    endif
         end if
c     read Data Set 4b, Qlimit
c
        Qlimit=MNW2(6,MNWID)
        IF(Qlimit.LT.0.0) THEN
          READ(in,*) Hlim,QCUT
          BACKSPACE in
          IF(QCUT.NE.0) THEN
            READ(in,*) Hlim,QCUT,Qfrcmn,Qfrcmx
          ELSE
            READ(in,*) Hlim,QCUT
            Qfrcmn=0.0
            Qfrcmx=0.0
          END IF
c     write info line
            write(iout,*) 
            write(iout,'(100A)') ' Qlimit < 0 : this well will
     & be constrained'
            write(iout,'(A,1PG13.5)') '     Hlim = ',Hlim 
            write(iout,2111) QCUT
 2111 FORMAT('     QCUT = ',I4)
           if(QCUT.lt.0) then 
           write(iout,'(A,1PG13.5,A)') '   Qfrcmn = ',Qfrcmn
           write(iout,'(A,1PG13.5,A)') '   Qfrcmx = ',Qfrcmx
           elseif(QCUT.gt.0) then
            write(iout,'(A,1PG13.5,A)') '   Qfrcmn = ',Qfrcmn, ' L**3/T'
            write(iout,'(A,1PG13.5,A)') '   Qfrcmx = ',Qfrcmx, ' L**3/T'
           end if
c     process min and max Q's based on QCUT
          MNW2(7,MNWID)=Hlim
          MNW2(8,MNWID)=QCUT
          MNW2(9,MNWID)=Qfrcmn
          MNW2(10,MNWID)=Qfrcmx
c     error check Qfrcmn
          if(QCUT.GT.0) then
           if(Qfrcmn.gt.1.d0) then
            write(iout,*) '***ERROR*** Qfrcmn > 1 is out of range'
            STOP 'MNW2 ERROR - Qfrcmn'
           end if
c     error check Qfrcmx
           if(Qfrcmx.gt.1.d0) then
            write(iout,*) '***ERROR*** Qfrcmx > 1 is out of range'
            STOP 'MNW2 ERROR - Qfrcmx'
           end if
	    end if
	  END IF
c
c     end read Data Set 4b
c
c        set desired flow rate for well
         MNW2(5,MNWID)=Qdes
c        if transport, set cprime for well
         if(igwtunit.gt.0) MNW2(12,MNWID)=Cprime
c        set all Qact (actual flow rate at node) = Qdes/NNODES as launching point
c        Loop over nodes in well
         firstnode=MNW2(4,MNWID)
c          put total Qdes in first node, allow it to go from there
         MNWNOD(4,firstnode)=Qdes
         write(iout,2112) WELLNAME,Qdes,kkper
 2112 FORMAT('MNW2 Well ', A20,' active, desired Q =',
     & 1pe12.4,' for stress period ',I4)
         if(PUMPCAP.GT.0) then
c  if Qdes is not negative for this stress period, do not apply pump capacity restraints
           if(Qdes.GE.0.d0) then
              MNW2(25,MNWID)=0
c   Initialize CapFlag2
              mnw2(27,MNWID)=0
           else
              MNW2(25,MNWID)=1
c  only set Qdes at upper end of CapTable if not set already
             if (CapTable(MNWID,PUMPCAP+2,2).LE.0.d0) then
               CapTable(MNWID,PUMPCAP+2,2)=abs(qdes)
             end if
             MNW2(24,MNWID)=CapMult
             write(iout,*) 
             write(iout,1114) mnw2(23,MNWID)
 1114 FORMAT('Reference head for calculating lift = ', 1pE12.4)
             write(iout,1113) CapMult
 1113 FORMAT('Pump Capacity Multiplier= ', 1pE12.4)
c   zero capacity use flag if CapMult=0
             if( CapMult.eq.0.d0) then
              MNW2(25,MNWID)=0
c   Initialize CapFlag2
              mnw2(27,MNWID)=0
             end if
c  now that we have Qdes, write Capacity table
             write(iout,*) 
             write(iout,*) 'Well Capacity Table'
             write(iout,*) ' Lift     Discharge'
             do index=1,PUMPCAP+2
                write(iout,'(1p2G10.3)') CapTable(MNWID,index,1),
     &                      CapTable(MNWID,index,2)
             end do
	     end if
	   end if
   	  end do
      end if
      return
c
      end
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2ad(nmnw2,MNWMAX,mnw2,NODTOT,MNWNOD,ibound,sc1,
     +                      delr,delc,BOTM,NBOTM,
     +                      hnew,hold,hstrt,ncol,nrow,nlay,small,kkstp,
     +                      CV,MNWINT,INTTOT,IOUT,WELLID,INBCF,ssflag,
     +                      LAYHDT,HK,VKA,kkper,MNWPRNT,NMNWVL)
C
C     VERSION 20090405 GZH modified from:
C     VERSION 20020819 KJH
c     
c----- MNW1 by K.J. Halford
c----- MNW2 by G.Z. Hornberger
c
c     ******************************************************************
c     Update Qact for wells that were constrained
c     ******************************************************************
c
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*20 WELLID
      INTEGER IBOUND,NMNW2,IW,MNWMAX,firstnode,lastnode,INODE,il,ir,ic,
     & NODTOT,ncol,nrow,nlay,NBOTM,LBOTM,LAYCBD,ipole,INTTOT,IOUT,kkstp,
     & itflag,INBCF,ssflag,LAYHDT,NNODES,kkper,MNWPRNT,NMNWVL
      REAL DELR,DELC,BOTM,CV,hold,hstrt,sc1,HK,VKA
      DOUBLE PRECISION mnw2,hnew,qoff,qon,qdes,csum,chsum,qact,Qsmall,
     & small,hwell,MNWNOD,verysmall,hlim,qpot,cond,ratio,hmax,hsim,
     & MNWINT,QCUT,qnet
      dimension mnw2(NMNWVL,MNWMAX),MNWNOD(33,NODTOT),MNWINT(11,INTTOT)
      dimension ibound(NCOL,NROW,NLAY),hnew(NCOL,NROW,NLAY),
     &  hstrt(ncol,nrow,nlay),SC1(ncol,nrow,nlay),
     &  CV(NCOL,NROW,NLAY),delr(ncol), delc(nrow),hold(NCOL,NROW,NLAY),
     & LAYHDT(NLAY),HK(NCOL,NROW,NLAY),VKA(NCOL,NROW,NLAY)
      dimension WELLID(mnwmax)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
C
      verysmall = 1.0D-8
c1------if number of wells <= 0 then return.
      if(nmnw2.le.0) return
c
c   Compute cell-to-well conductance for each well node
c
c   Send in ITFLAG=0, this means don't calculate partial penetration effects
c   because we are not currently iterating on a solution
      ITFLAG=0
      call SMNW2COND(delr,delc,hnew,hold,hstrt,sc1,
     +                 ncol,nrow,nlay,kkstp,ssflag,kkper,
     &                 BOTM,NBOTM,ibound,MNWINT,INTTOT,
     &                 MNWNOD,NODTOT,nmnw2,MNWMAX,MNW2,CV,IOUT,ITFLAG,
c sending in kiter=0; this is before iter loop
     &                 INBCF,LAYHDT,WELLID,0,HK,VKA,MNWPRNT,NMNWVL)

c
c
c   Allow constrained wells a new chance with the next time step
c
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
        if (MNW2(1,iw).EQ.1) then
          qdes = mnw2(5,iw)
c   Define qnet for well, which may be updated if well was restricted in previous time step
          qnet = mnw2(5,iw)
C  SET Q
c  if qlimit was hit in a previous FM routine, set qnet to updated Q, then reset flags
          if(mnw2(20,iw).ne.0.or.MNW2(27,iw).gt.0) then
            qnet=mnw2(18,iw)
            mnw2(20,iw)=0
            mnw2(27,iw)=0
c  default Q used = Qdes (this may be changed if limits kick in below)
          else
            mnw2(18,iw)=mnw2(5,iw)
          end if
c
c   Retrieve QCUT, Qfrcmn, Qfrcmx
          if(MNW2(6,iw).gt.0) then
            QCUT = MNW2(8,iw)
            if(QCUT.ne.0) then
             qoff = mnw2(9,iw)
             qon  = mnw2(10,iw)
             if(QCUT.GT.0) then
c     convert rate into fraction of Qdes (ratio is used to compare)
              if(Qdes.NE.0) then
                qoff=abs(qoff/Qdes)
                qon=abs(qon/Qdes)
              else
                qoff=0.0
                qon=0.0
              end if
             elseif (QCUT.LT.0) then
c     convert percentage into fraction (ratio is used to compare)
c              qoff=qoff*0.01
c              qon=qon*0.01
             end if
            end if
          end if
c
c   Compute hwell / Qpot for multi-node well (not single-cell wells)
c   
c         if NNODES>1, that is, this is a multi-node well
c      write(*,*) 'MNW2(2,iw)',MNW2(2,iw)
          NNODES=abs(MNW2(2,iw))
          if(NNODES.gt.1) then
c      write(*,*) 'MNW2(6,iw)',MNW2(6,iw)
            csum = 0.000D0
            chsum = 0.000D0
            qact = 0.0000D0
            Qsmall = small*abs(qdes)
c   Loop over nodes in well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+NNODES-1
            hwell = mnw2(17,iw)
            do INODE=firstnode,lastnode
              il=MNWNOD(1,INODE)              
              ir=MNWNOD(2,INODE)              
              ic=MNWNOD(3,INODE)              
              if(IBOUND(ic,ir,il).ne.0) then
                csum  = csum  + MNWNOD(14,INODE)
                chsum = chsum + MNWNOD(14,INODE)*hnew(ic,ir,il)
                qact  = qact  + MNWNOD(4,INODE)
              else
                qact  = 0.0000D0
              end if            
            end do
c---div0 ---  CSUM could go to zero if the entire well is dry
            if( csum .gt. verysmall ) then
c           for limit procedure, use qnet here which may have been
c           restricted in previous time step
c              hwell = ( qdes + chsum ) / csum
              hwell = ( qnet + chsum ) / csum
            else
              hwell = hnew(ic,ir,il)
            endif
c      Test Hlim constraint if QLIMIT flag is set
            if(MNW2(6,iw).GT.0) then
              hlim = mnw2(7,iw)
              ipole = 0
              if( abs(qdes).gt.verysmall ) ipole = qdes / abs(qdes)
              hmax = ipole*( hlim )
              hsim = ipole*( hwell )
c      Potential Q is...
              if( hsim .gt. hmax ) then
                hwell = hlim
              endif
              qpot = hlim*csum - chsum
            end if
            cond = csum
c   check cond<0, reset to 0 and print warning
              if(cond.lt.0.d0) then
                write(iout,*) '***WARNING*** CWC<0 reset to CWC=0'
                write(iout,*) 'In Well ',WELLID(iw),' Node ',INODE
                cond=0.d0
              end if
c   Else, dealing with a single-cell well
          else
            qact = mnw2(5,iw)
            Qsmall = small
c     Compute hwell / Qpot for single-node well
            firstnode=MNW2(4,iw)
            il=MNWNOD(1,firstnode)              
            ir=MNWNOD(2,firstnode)              
            ic=MNWNOD(3,firstnode)              
            cond = MNWNOD(14,firstnode)
c   check cond<0, reset to 0 and print warning
              if(cond.lt.0.d0) then
                write(iout,*) '***WARNING*** CWC<0 reset to CWC=0'
                write(iout,*) 'In Well ',WELLID(iw),' Node ',INODE
                cond=0.d0
              end if
            if(MNW2(6,iw).GT.0) then
              hlim = mnw2(7,iw)
              qpot = ( hlim - hnew(ic,ir,il) )*cond
            end if
          end if
c
c  Compute ratio of potential/desired flow rates
          if(MNW2(6,iw).GT.0) then
            ratio = 1.00D0
            if( abs(qdes) .gt. small ) ratio =  qpot / qdes
            if( ratio .gt. 0.9999D0 ) then
              ratio =  1.000D0
              Qpot = Qdes
            endif
c  Check if potential flow rate is below cutoff
c  If so, set Qact = 0.0
            firstnode=MNW2(4,iw)
            if( ratio .lt. Qoff ) then
c              mnw2(18,iw) = 0.0
              mnw2(30,iw) = 0.0
c  set qact in node too
              MNWNOD(4,firstnode)=0.0
c  Check if potential flow rate is above restart threshold
c  If so, set Qact = Qpot
            elseif( ratio.gt.Qon .and. (abs(qact).lt.Qsmall)) then
c                mnw2(18,iw) = Qpot
                mnw2(30,iw) = Qpot
                MNWNOD(4,firstnode)=Qpot
            else
c  Otherwise leave the flow rate alone
            endif
c  End if, QLimit>0
c  with QCUT=0 now, set q=qpot if neither situation holds?
            if(QCUT.EQ.0.and.ratio.gt.0.D0) then
c              mnw2(18,iw) = Qpot
              mnw2(30,iw) = Qpot
              MNWNOD(4,firstnode)=Qpot              
            end if
          endif
c  End if, active wells
        end if
c  End do, loop over all wells
      end do
c
      return
      end
c
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2fm(nmnw2,ncol,nrow,nlay,delr,delc,sc1,
     & BOTM,NBOTM,ibound,MNWNOD,NODTOT,MNWMAX,MNW2,kiter,hnew,hold,
     $ hstrt,hcof,rhs,small,CV,MNWINT,INTTOT,iout,kkstp,INBCF,ssflag,
     & LAYHDT,WELLID,HK,VKA,kkper,hdry,CapTable,MNWPRNT,
     & hclose,NMNWVL)
C     VERSION 20020819 KJH
C     VERSION 20090405 GZH
c
c----- MNW1 by K.J. Halford
c----- MNW2 by G.Z. Hornberger
c
c     ******************************************************************
c     add well flow to source term
c     ******************************************************************
c
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER nmnw2,ncol,nrow,nlay,nbotm,ibound,NODTOT,MNWMAX,iw,LBOTM,
     & kiter,firstnode,lastnode,inode,il,ir,ic,ipole,iqslv,LAYCBD,
     & INTTOT,iout,ITFLAG,kkstp,INBCF,ssflag,LAYHDT,kkper,NNODES,
     & MNWPRNT,PUMPCAP,NMNWVL
      REAL delr,delc,botm,hcof,rhs,CV,hold,hstrt,sc1,HK,VKA,hdry,hclose
      DOUBLE PRECISION verysmall,MNWNOD,MNW2,hnew,qdes,hlim,hwell,qact,
     & cond,dhc2w,hmax,hsim,ratio,small,MNWINT,CapTable,qactCap,Hlift,
     & CapMult,lastQ,qtemp,temppct,lastH,htemp,HWtol
      dimension mnw2(NMNWVL,MNWMAX),MNWNOD(33,NODTOT),MNWINT(11,INTTOT)
      DIMENSION delr(ncol), delc(nrow), CV(NCOL,NROW,NLAY)
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),IBOUND(NCOL,NROW,NLAY),
     & hnew(NCOL,NROW,NLAY),hcof(NCOL,NROW,NLAY),rhs(NCOL,NROW,NLAY),
     & hold(NCOL,NROW,NLAY),hstrt(ncol,nrow,nlay),SC1(ncol,nrow,nlay),
     & LAYHDT(nlay),HK(NCOL,NROW,NLAY),VKA(NCOL,NROW,NLAY)
      dimension CapTable(mnwmax,27,2)
      CHARACTER*20 WELLID
      dimension WELLID(mnwmax)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
C
      verysmall = 1.0D-20
c
c                 CR( i, j, k)    ------>   CR  i + 1/2
c                 CC( i, j, k)    ------>   CC  j + 1/2
c                 CV( i, j, k)    ------>   CV  k + 1/2
c
c1------if number of wells <= 0 then return.
      if(nmnw2.le.0) return
c

c   Compute cell-to-well conductance for each well node
c
c   Send in ITFLAG=1, this means calculate partial penetration effects
c   because we are currently iterating on a solution
      ITFLAG=1
      call SMNW2COND(delr,delc,hnew,hold,hstrt,sc1,
     +                 ncol,nrow,nlay,kkstp,ssflag,kkper,
     &                 BOTM,NBOTM,ibound,MNWINT,INTTOT,
     &                 MNWNOD,NODTOT,nmnw2,MNWMAX,MNW2,CV,IOUT,ITFLAG,
     &                 INBCF,LAYHDT,WELLID,kiter,HK,VKA,MNWPRNT,NMNWVL)
c
c   Prepare components and limits of a multi-node well
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
        if (MNW2(1,iw).EQ.1) then
c   Check capacity calculation flag
         PUMPCAP=MNW2(22,iw)
         CapMult=MNW2(24,iw)
         qdes = mnw2(5,iw)
c   If CapFlag was turned off (e.g. with % Q check), turn back on at beginning of new time step (kiter=1)
         if(kiter.eq.1.and.PUMPCAP.GT.0.and.CapMult.ne.0.d0.
     &      and.qdes.lt.0.d0) then
            MNW2(25,iw)=1
         end if
c  compute well capacity restraint
c         if capacity should be calculated, wait until after first iteration so
c         the lift can be computed with a new hwell (kiter>1)
          if(MNW2(25,iw).GT.0.and.kiter.gt.1) then
            hwell = mnw2(17,iw)
c           save hwell before update in SEEP routine, to check vs HWtol later
            lastH=hwell
            Hlift=MNW2(23,iw)
c           call MNW2CAPACITY to look up q based on capacity table (returns qactCap)
            call MNW2CAPACITY(qactCap,WELLID,mnwmax,Hlift,
     & hwell,iw,PUMPCAP,CapTable,iout)
c           get last Q that was looked up in table; this will be used to see if
c           new Q is less than 1% or more than 25% from new Q
            lastQ=mnw2(26,iw)
c           save last Q that was looked up in table for next iteration
            mnw2(26,iw)=qactCap
c           only do following if capacity Q is less than Qdes
            if (abs(qactCap).lt.abs(mnw2(5,iw))) then
c            if capacity Q from table is zero at beginning first time step, 
c            (kiter=2 is first time here due to kiter>1 check above)
c            set qpot to 0; "allow" to go to zero
             if(qactCap.eq.0.d0.and.kiter.eq.2.and.kkstp.gt.1) then
              mnw2(29,iw)=0.d0
c             set CapFlag2 to true, this is used in AD routine to say that the
c             Q has been constrained
              mnw2(27,iw)=1
             else
c             only perform % checks after second iteration (lets soln get going)
              if(kiter.gt.2) then
c              % check divides by lastQ, so avoid this if it is zero
               if(lastQ.ne.0.0) then
c               this is the percentage change, below which the Q will be locked in
c               for the time step
                temppct=0.01d0
c               calculated % change in Q looked up in capacity tables
                qtemp=abs(qactCap-lastQ)/abs(lastQ)
c
c               if Q is changing more than 1%, continue
                if (qtemp.gt.temppct) then
c                 if Q is changing more than 25% after 1st TS, constrain change
c                 in Q to 25% increase/decrease
                  if(qtemp.gt.0.25d0.and.kkstp.gt.1) then
c                   check to see which way Q is moving; adjust accordingly
                    if(lastQ.lt.qactCap) then
                      mnw2(29,iw)=lastQ*0.75d0
                    else
                      mnw2(29,iw)=lastQ*1.25d0
                    end if
c                   set CapFlag2 to true, this is used in AD routine to say that the
c                   Q has been constrained
                    mnw2(27,iw)=1
c                   save last Q
                    mnw2(26,iw)=mnw2(29,iw)
c  write message that Q was constrained
                    if(mnwprnt.gt.1) then
      write(iout,*) 'Capacity Q was changing more than 25%, constrained
     & to 25% increase/decrease for well ',WELLID(iw)
                    end if
                  else
c               if Q is changing less than 1%, set equal to Q just looked up
                    mnw2(29,iw)=qactCap
c                   set CapFlag2 to true, this is used in AD routine to say that the
c                   Q has been constrained
                    mnw2(27,iw)=1
                  end if
                else
c  if Qcap not changing more than 1%, stop doing capacity constraint
                  mnw2(29,iw)=qactCap
c  set CapFlag to false to discontinue doing capacity checks this TS
                  mnw2(25,iw)=0
c                 set CapFlag2 to true, this is used in AD routine to say that the
c                 Q has been constrained
                  mnw2(27,iw)=1
c  write message that Q was constrained
                  if(mnwprnt.gt.1) then
      write(iout,*) 'Capacity Q was changing less than 1%, capacity 
     & calculations stopped this TS   for well ',WELLID(iw)
                  end if
                end if
c              else lastQ=0               
               else
c                if lastQ is zero, but new Q is not, and at first time this
c                happens in TS, only let it turn on at half the Q looked up
                 if(qactCap.lt.0.d0.and.kiter.eq.3) then
                   mnw2(29,iw)=qactCap*0.5d0
c                  set CapFlag2 to true, this is used in AD routine to say that the
c                  Q has been constrained
                   mnw2(27,iw)=1
c  write message that Q was constrained
                   if(mnwprnt.gt.1) then
      write(iout,*) 'Capacity Q went from zero to nonzero, new Q  
     & constrained to half of value looked up for well ',WELLID(iw)
                   end if
                 end if
               end if
              end if
             end if
            end if
          end if
c         if NNODES>1, that is, this is a multi-node well
          NNODES=abs(MNW2(2,iw))
          if(NNODES.gt.1) then
            CALL SMNW2SEEP(MNW2,MNWNOD,iw,MNWMAX,NODTOT,IBOUND,
     &   ncol,nrow,nlay,hnew,kiter,BOTM,NBOTM,kkstp,
     &   kkper,small,hdry,NMNWVL)
c  for PUMPCAP option, check hwell vs last H
            if(MNW2(25,iw).GT.0.and.kiter.gt.2) then
              hwell = mnw2(17,iw)
              htemp=abs(lastH-hwell)
              HWtol=mnw2(28,iw)
c  if htemp (head change) is less than closure criterion stop doing capacity constraint
              if (htemp.le.HWtol) then
                mnw2(29,iw)=qactCap
                mnw2(25,iw)=0
                mnw2(27,iw)=1
c  write message that Q was constrained
                if(mnwprnt.gt.1) then
      write(iout,*) 'Capacity Q calculations stopped due to hwell
     & change less than HWtol for well ',WELLID(iw)
                end if
              end if
            end if
          end if
        end if
      end do
c
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1 and IBOUND>0
        if (MNW2(1,iw).EQ.1) then
          qdes = mnw2(5,iw)
c   If Capacity restrictions are set, Qdes here is actually the retricted Qpot
          if(mnw2(27,iw).ne.0) qdes = mnw2(29,iw)
          hwell = mnw2(17,iw)
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
c   Loop over nodes in well
          do INODE=firstnode,lastnode
            il=MNWNOD(1,INODE)              
            ir=MNWNOD(2,INODE)              
            ic=MNWNOD(3,INODE)              
            if(IBOUND(ic,ir,il).ne.0) then
              qact = MNWNOD(4,INODE)
              cond = MNWNOD(14,INODE)
c   Process single-node wells (MNW2(2,iw).eq.1, i.e. NUMNODES)
              if( abs(MNW2(2,iw)).eq.1.and.cond.gt.verysmall ) then
                dhc2w = Qact / cond
                hwell = hnew(ic,ir,il) + dhc2w
                MNW2(17,iw) = hwell
                iqslv=0
c   Test DD constraints, Hlim is assumed to be a Max/Min for Injection/Production wells
c   Set hwell in mnw2 array
                if(MNW2(6,iw).GT.0) then
                  hlim = mnw2(7,iw)                
                  ipole = 0
                  if( abs(qdes).gt.verysmall ) ipole = qdes / abs(qdes)
                  hmax = ipole*( hlim )
                  hsim = ipole*( hwell )
c   Calculate ratio of actual flow to desired flow
                  ratio = 1.00D0
                  if( abs(qdes) .gt. verysmall ) ratio =  qact / qdes
                  if( abs(ratio).gt. 1.00D0 ) qact = qdes
                  if( ratio .lt. verysmall ) qact = 0.0D0
c   If iqslv=0, Well will be simulated as a specified rate 
c   If iqslv=1, Well will be simulated as a GHB 
                  iqslv = 0
                  if( hsim.gt.hmax .and. hmax.gt.verysmall ) iqslv = 1
                  if((qdes-qact)**2 .gt. small           ) iqslv = 1
                  if(abs(qact).lt.verysmall .and.hsim.gt.hmax) iqslv = 0
                  if(abs(qact).lt.verysmall .and.hsim.lt.hmax) iqslv = 1
                  if(abs(qdes).lt.verysmall .or. 
     &             ratio.gt.1.0D0-verysmall ) iqslv =0
                else
                  hlim=mnw2(17,iw)
                end if
c 
              elseif( cond.lt.verysmall ) then
                qact = 0.0D0
                iqslv = 0
              else
c Process multi-node wells, Constraints were already tested when allocating flow
                if( mod(kiter,2).eq.0 .and. abs(qact).gt.small ) then
c                  hlim = mnw2(17,iw)    !!   Water level in wellbore 
                  hlim = MNWNOD(15,INODE)    !!   Water level in wellbore, or bottom of cell if seepage face cell
                  iqslv = 1
                else
                  qact = MNWNOD(4,INODE)
                  iqslv = 0
                endif
              end if
c
c   Modify HCOF and RHS arrays
cdebug replace below line with this if debugging with No Mo Iter
cdebug              if(iqslv.ne.0.and.kiter.gt.1.and.kiter.lt.NoMoIter) then
              if(iqslv.ne.0.and.kiter.gt.1) then
                qact = ( hlim - hnew(ic,ir,il) ) * cond
                hcof(ic,ir,il) = hcof(ic,ir,il) - cond
                rhs(ic,ir,il)  = rhs(ic,ir,il)  - cond * hlim
              else
c  Specify Q and solve for head;  add Q to RHS accumulator.
                rhs(ic,ir,il) = rhs(ic,ir,il) - qact
              endif
              MNWNOD(4,INODE) = qact
            end if
          end do
        end if
      end do
c
      return
      end
c
c
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2bd(IMNWCB,ICBCFL,NAUX,KSTP,KPER,
     & NCOL,NROW,NLAY,nmnw2,mnw2,iout,DELT,PERTIM,TOTIM,IBOUND,
     & MNWMAX,NODTOT,msum,HDRY,VBNM,VBVL,BUFF,MNWNOD,hnew,WELLID,
     & MNWPRNT,IGWTON,NMNWVL,MNWOBS,gwtunit)

C     VERSION 20030710 KJH
C     VERSION 20090405 GZH
c
c----- MNW1 by K.J. Halford        1/31/98
c----- MNW2 by G.Z. Hornberger     9/13/07
c     ******************************************************************
c     calculate volumetric budget for multi-node wells
c     ******************************************************************
c
c        specifications:
c     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER ibd,IMNWCB,ICBCFL,NAUX,NCOL,NROW,NLAY,KSTP,KPER,nmnw2,
     & iout,ibound,MNWMAX,imult,iw,firstnode,lastnode,inode,il,ir,ic,
     & NODTOT,ioch,ipole,ioc,NMNWVL,msum,nd,iweldry,MNWPRNT,MNWOBS,
     & gwtunit,igwton

      REAL HDRY, VBVL, BUFF, PERTIM, TOTIM, DELT
      DOUBLE PRECISION ratin,ratout,mnw2,MNWNOD,hnew,DryTest,
     & q,qdes,hlim,hwell,s,sNL,sL,qnet,qin,qout,verysmall,hcell,small
      CHARACTER*20 WELLID
      CHARACTER*16 text,vbnm(msum)
      DIMENSION ibound(ncol,nrow,nlay),mnw2(NMNWVL,MNWMAX),
     & buff(NCOL,NROW,NLAY),MNWNOD(33,NODTOT),hnew(NCOL,NROW,NLAY),
     & vbvl(4,msum)
      dimension WELLID(mnwmax)
      COMMON /MNW2COM/MNWAUX(5)
      CHARACTER*16 MNWAUX
c             ----+----1----+-
      text = '            MNW2'
      small = 1.D-10
      verysmall = 1.D-25
c     ------------------------------------------------------------------
c
c  clear ratin and ratout accumulators.
      ratin=0.D0
      ratout=0.D0
      ibd=0
      IF(IMNWCB.GT.0) IBD=ICBCFL
C
C2-----IF CELL-BY-CELL FLOWS WILL BE SAVED AS A LIST, WRITE HEADER.
      IF(IBD.EQ.2) THEN
         NAUX = 0   !!   Set to zero -- Change order to dump
c         IF(IAUXSV.EQ.0) NAUX=0
         CALL UBDSV4(KSTP,KPER,TEXT,NAUX,MNWAUX,IMNWCB,NCOL,NROW,NLAY,
     1          nmnw2,IOUT,DELT,PERTIM,TOTIM,IBOUND)
      END IF
c  clear the buffer.
      buff=0.000000000
c
c2------if there are no wells do not accumulate flow
      if(nmnw2.gt.0) then
c  test for dry wells
        imult = 0
c  Loop over all wells
        do iw=1,MNWMAX
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
c   Loop over nodes in well
          do INODE=firstnode,lastnode
            il=MNWNOD(1,INODE)              
            ir=MNWNOD(2,INODE)              
            ic=MNWNOD(3,INODE)              
            DryTest = Hnew(ic,ir,il) - Hdry
            if(ABS(DryTest).lt.verysmall) then
              MNWNOD(4,INODE)= 0.0D0
              MNW2(17,iw)=Hdry
              iweldry=1
              if(MNWPRNT.gt.1) then
                write(iout,*) 'MNW2 node in dry cell, Q set to 0.0'
                write(iout,*) 'Well: ',WELLID(iw),' Node: ',INODE
              end if
            else
              iweldry=0
            endif
              q = MNWNOD(4,INODE)
c
              buff(ic,ir,il) = buff(ic,ir,il) + q
              if( q.ge.0.0D0 ) then
c -----pumping rate is positive(recharge). add it to ratin.
                ratin = ratin + q
              else
c -----pumping rate is negative(discharge). add it to ratout.
                ratout = ratout - q
              endif
          enddo
c warn about dry well
            if(MNWPRNT.gt.1.and.iweldry.eq.1) then
              write(iout,*) 'Note-- the following MNW2 well went dry:',
     &  WELLID(iw)
            end if
c if MNWOBS or GWT turned on, calculate borehole flow
          if(MNWOBS.gt.0.or.gwtunit.gt.0)  
     &      call GWF1MNW2BH(MNWMAX,mnw2,NODTOT,MNWNOD,iout,NMNWVL,iw)

        enddo
c6------if cell-by-cell flows will be saved call ubudsv to record them
        if( abs(IMNWCB).gt.0 .and. icbcfl.ne.0 ) then           
          ioc = abs(IMNWCB)
          if( ibd.eq.2 ) then   !!  Write COMPACT budget
c   Loop over all wells
            do iw=1,MNWMAX
c   active well check
             if(MNW2(1,iw).gt.0) then
c   Loop over nodes in well
              firstnode=MNW2(4,iw)
              lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
              do INODE=firstnode,lastnode
                ic=MNWNOD(3,INODE)
                ir=MNWNOD(2,INODE)
                il=MNWNOD(1,INODE)
                call UBDSVB(ioc,ncol,nrow,IC,IR,IL,real(Q),real(IL),
     +                    NMNWVL,NAUX,32,IBOUND,NLAY)
              end do
             endif
            enddo
          else                  !!  Write full 3D array
            call ubudsv(kstp,kper,text,ioc, buff,ncol,nrow,nlay,iout)
          endif
        endif
      endif
c
c7------move rates into vbvl for printing by module bas1ot.
      vbvl(3,msum)=ratin
      vbvl(4,msum)=ratout
c
c8------move rates times time step length into vbvl accumulators.
      vbvl(1,msum) = vbvl(1,msum) + ratin*delt
      vbvl(2,msum) = vbvl(2,msum) + ratout*delt
c
c9------move budget term labels into vbnm for printing.
      vbnm(msum) = text
c
c10-----increment budget term counter(msum).
      msum = msum + 1
c
c11-----return
      return
      end
c
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2BCF(delr,delc,cr,cc,hy,hnew,
     +                 ncol,nrow,nlay,Hdry,small,
     &                 LAYHDT,BOTM,NBOTM,TRPY,
     &   MNWNOD,NODTOT,nmnw2,MNWMAX,MNW2,WELLID,MNWPRNT,iout,NMNWVL)
C     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
C     VERSION 20090405 GZH        -- MNW2
c
c----- MNW1 by K.J. Halford
c
c     ******************************************************************
c     Compute transmissivities used to calculate cell-to-well conductance
c     ******************************************************************
C     Note: BCF, when LAYCON=0 or 2, does not save cell-by-cell
C     Transmissivity (T) values.  Instead, it converts the cell-by-cell
C     T values to branch conductances CR and CC, using harmonic
C     averaging.  When BCF is used, the method used in this routine to
C     generate cell-specific values of Tx and Ty is an approximation
C     based on CR and CC.  When LPF or HUF is used, cell-by-cell
C     hydraulic-conductivity values are stored, this approximation is
C     not needed, and the values generated for Tx and Ty are exact --
C     ERB 1/29/01.
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*20 WELLID
      dimension WELLID(mnwmax)
      INTEGER ix,iy,iz,LBOTM,LAYHDT,ncol,nrow,nlay,LAYCON,NODTOT,NBOTM,
     & nmnw2,MNWMAX,iw,firstnode,lastnode,INODE,LAYCBD,iout,MNWPRNT,
     & NMNWVL
      REAL DELR,DELC,BOTM,TRPY,CR,CC,hy,Hdry
      DOUBLE PRECISION verysmall,dx,dy,top,bot,ah,dxp,Txp,dxm,Txm,dyp,
     & Typ,dym,Tym,small,Txx,div,Tyy,upper,hnew,TempKX,thick,MNWNOD,MNW2
      dimension delr(ncol), delc(nrow),
     & cr(NCOL,NROW,NLAY),cc(NCOL,NROW,NLAY),hy(NCOL,NROW,NLAY),
     & hnew(NCOL,NROW,NLAY),MNWNOD(33,NODTOT),mnw2(NMNWVL,MNWMAX)
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),
     &          LAYHDT(NLAY), TRPY(NLAY)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
      COMMON /BCFCOM/LAYCON(999)
C     ------------------------------------------------------------------
      verysmall = 1.D-25
c
c1------if number of wells <= 0 then return.
      if(nmnw2.le.0) return
C
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
c   Loop over nodes in well
        firstnode=MNW2(4,iw)
        lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
        do INODE=firstnode,lastnode
         ix=MNWNOD(3,INODE)
         iy=MNWNOD(2,INODE)
         iz=MNWNOD(1,INODE)
         dx   = delr(ix)
         dy   = delc(iy)
         top = BOTM(IX,IY,LBOTM(IZ)-1)
         bot = BOTM(IX,IY,LBOTM(IZ))
C
C     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
         AH = TRPY(IZ)
C
         if (LAYHDT(IZ).EQ.0) then
C       THICKNESS IS NOT HEAD-DEPENDENT
          dxp  = dx
          Txp  = 0.00000D0
          if( ix .lt. ncol ) then
            dxp = delr(ix+1)
            Txp  = cr(ix,iy,iz) * (dx+dxp) / 2.D0
          endif
          dxm = dx
          Txm  = Txp
          if( ix .gt. 1  ) then
            dxm = delr(ix-1)
            Txm  = cr(ix-1,iy,iz) * (dx+dxm) / 2.D0
          endif
          if( Txp.lt.small ) Txp = Txm
          if( Txm.lt.small ) Txm = Txp
c
          dyp  = dy
          Typ  = 0.00000D0
          if( iy .lt. nrow ) then
            dyp = delc(iy+1)
            Typ  = cc(ix,iy,iz) * (dy+dyp) / 2.D0
          endif
          dym = dy
          Tym  = Typ
          if( iy .gt. 1 ) then
            dym = delc(iy-1)
            Tym  = cc(ix,iy-1,iz) * (dy+dym) / 2.D0
          endif
          if( Typ.lt.small ) Typ = Tym
          if( Tym.lt.small ) Tym = Typ
          Txp = Txp / dy
          Txm = Txm / dy
          Typ = Typ / dx
          Tym = Tym / dx
c
c  Eliminate zero values .....
c
          if( Typ.lt.small .or. nrow.lt.2 )  then
            Typ = Txp
            Tym = Txm
          endif
c
          if( Txp.lt.small .or. ncol.lt.2 )  then
            Txp = Typ
            Txm = Tym
          endif
c
c   Assuming expansion of grid is slight, if present, & that Txx and Tyy of the adjacent
c   cells are about the same value.
          Txx = 0.00000000D0
          div  = Txp + Txm
          if( div.gt.small ) Txx  = 2*Txp*Txm / div
          Tyy = 0.00000000D0
          div  = Typ + Tym
          if( div.gt.small ) Tyy  = 2*Typ*Tym / div
          if( Txx.gt.small .and. Tyy.lt.small ) Tyy = Txx
          if( Tyy.gt.small .and. Txx.lt.small ) Txx = Tyy
         else
C       THICKNESS IS HEAD-DEPENDENT
c  Estimate T to well in an unconfined system
c
          upper = hnew(ix,iy,iz)
          if (LAYCON(IZ).EQ.3) then
           if( upper.gt.top ) upper = top
          endif
          TempKX = hy(ix,iy,iz)      !! BCF Hydraulic Conductivity array
          thick = upper - bot
c   set thickness / conductance to 0 if cell is dry
          if( (hnew(ix,iy,iz)-Hdry )**2 .lt. verysmall ) 
     &       thick = 0.0000000000D0
          Txx = TempKX * thick
          if( Txx .lt.verysmall ) Txx = 0.000000000000000D0
          Tyy = Txx * AH
         endif
         MNWNOD(16,INODE)=Txx
         MNWNOD(17,INODE)=Tyy
        end do 
       end if 
      end do 
      return
      end
c
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2LPF(delr,delc,hnew,
     +                 ncol,nrow,nlay,Hdry,
     &                 LAYHDT,BOTM,NBOTM,HK,HANI,
     &     MNWNOD,NODTOT,nmnw2,MNWMAX,MNW2,WELLID,MNWPRNT,iout,NMNWVL)
C
C     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
C     VERSION 20090405 GZH        -- MNW2
c
c----- MNW1 by K.J. Halford
c
c     ******************************************************************
c     Compute transmissivities used to calculate cell-to-well conductance
c     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
C     ------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*20 WELLID
      dimension WELLID(mnwmax)
      INTEGER ix,iy,iz,LBOTM,LAYHDT,ncol,nrow,nlay,LAYCON,LAYWET,LAYVKA,
     &        LAYAVG,LAYTYP,NODTOT,NBOTM,LAYCBD,
     & nmnw2,MNWMAX,iw,firstnode,lastnode,INODE,iout,MNWPRNT,NMNWVL
      REAL DELR,DELC,BOTM,CHANI,HANI,Hdry,HK
      DOUBLE PRECISION verysmall,dx,dy,top,bot,ah,
     & Txx,Tyy,upper,hnew,TempKX,thick,MNWNOD,MNW2
      dimension delr(ncol), delc(nrow),
     & hnew(NCOL,NROW,NLAY),HK(NCOL,NROW,NLAY),
     & MNWNOD(33,NODTOT),mnw2(NMNWVL,MNWMAX)
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),
     &          LAYHDT(NLAY),HANI(NCOL,NROW,NLAY)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
      COMMON /LPFCOM/LAYTYP(999),LAYAVG(999),CHANI(999),LAYVKA(999),
     1               LAYWET(999)
C     ------------------------------------------------------------------
      verysmall = 1.D-25
c
c1------if number of wells <= 0 then return.
      if(nmnw2.le.0) return
C
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
c   Loop over nodes in well
        firstnode=MNW2(4,iw)
        lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
        do INODE=firstnode,lastnode
         ix=MNWNOD(3,INODE)
         iy=MNWNOD(2,INODE)
         iz=MNWNOD(1,INODE)
         dx   = delr(ix)
         dy   = delc(iy)
         top = BOTM(IX,IY,LBOTM(IZ)-1)
         bot = BOTM(IX,IY,LBOTM(IZ))
C
C     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
         AH = 1.0D0
         IF (CHANI(IZ).GT.0.0) THEN
          AH = CHANI(IZ)
         ELSE
          AH = HANI(IX,IY,IZ)
         ENDIF
C
         if (LAYHDT(IZ).EQ.0) then
C       THICKNESS IS NOT HEAD-DEPENDENT
          THICK = TOP-BOT
          TXX = HK(ix,iy,iz)*THICK
          TYY = TXX*AH
         else
C       THICKNESS IS HEAD-DEPENDENT
c  Estimate T to well in an unconfined system
c
          upper = hnew(ix,iy,iz)
          if( upper.gt.top ) upper = top
          TempKX = hk(ix,iy,iz)       !!LPF Hydraulic Conductivity array
          thick = upper - bot
c   set thickness / conductance to 0 if cell is dry
          if( (hnew(ix,iy,iz)-Hdry )**2 .lt. verysmall ) 
     &          thick = 0.0000000000D0
          Txx = TempKX * thick
          if( Txx .lt.verysmall ) Txx = 0.000000000000000D0
          Tyy = Txx * AH
         endif
         MNWNOD(16,INODE)=Txx
         MNWNOD(17,INODE)=Tyy
        end do
       end if
      end do 
      return
      end
c
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2HUF(delr,delc,hnew,
     +                 ncol,nrow,nlay,Hdry,
     &                 LAYHDT,BOTM,NBOTM,HK,HKCC,
     &     MNWNOD,NODTOT,nmnw2,MNWMAX,MNW2,WELLID,MNWPRNT,iout,NMNWVL)
C     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
C     VERSION 20090405 GZH        -- MNW2
c
c----- MNW1 by K.J. Halford
c
c     ******************************************************************
c     Compute transmissivities used to calculate cell-to-well conductance
c     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*20 WELLID
      dimension WELLID(mnwmax)
      INTEGER ix,iy,iz,LBOTM,LAYHDT,ncol,nrow,nlay,LAYCON,NODTOT,
     &nmnw2,MNWMAX,iw,firstnode,lastnode,INODE,LAYCBD,NBOTM,iout,MNWPRNT
     &,NMNWVL
      REAL DELR,DELC,BOTM,TRPY,HK,HKCC,Hdry,KY
      DOUBLE PRECISION verysmall,dx,dy,top,bot,ah,
     & Txx,Tyy,upper,hnew,TempKX,thick,MNWNOD,MNW2
      dimension delr(ncol), delc(nrow),
     & hnew(NCOL,NROW,NLAY),MNWNOD(33,NODTOT),mnw2(NMNWVL,MNWMAX)
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),
     &          LAYHDT(NLAY),HK(NCOL,NROW,NLAY),
     &          HKCC(NCOL,NROW,NLAY)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)

C     ------------------------------------------------------------------
C
      verysmall = 1.D-25
c
c1------if number of wells <= 0 then return.
      if(nmnw2.le.0) return
C
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
c   Loop over nodes in well
        firstnode=MNW2(4,iw)
        lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
        do INODE=firstnode,lastnode
         ix=MNWNOD(3,INODE)
         iy=MNWNOD(2,INODE)
         iz=MNWNOD(1,INODE)
         dx   = delr(ix)
         dy   = delc(iy)
         top = BOTM(IX,IY,LBOTM(IZ)-1)
         bot = BOTM(IX,IY,LBOTM(IZ))
C
C     FIND HORIZONTAL ANISOTROPY, THE RATIO Ky/Kx
         TempKX = HK(ix,iy,iz)
         KY = HKCC(IX,IY,IZ)
         AH = KY/TempKX
C
         if (LAYHDT(IZ).EQ.0) then
C       THICKNESS IS NOT HEAD-DEPENDENT
          THICK = TOP-BOT
          TXX = HK(ix,iy,iz)*THICK
          TYY = TXX*AH
         else
C       THICKNESS IS HEAD-DEPENDENT
c  Estimate T to well in an unconfined system
c
          upper = hnew(ix,iy,iz)
          if( upper.gt.top ) upper = top
          TempKX = hk(ix,iy,iz)      !!HUF Hydraulic Conductivity array
          thick = upper - bot
c   set thickness / conductance to 0 if cell is dry
          if( (hnew(ix,iy,iz)-Hdry )**2 .lt. verysmall ) 
     &           thick = 0.0000000000D0
          Txx = TempKX * thick
          if( Txx .lt.verysmall ) Txx = 0.000000000000000D0
          Tyy = Txx * AH
         endif
         MNWNOD(16,INODE)=Txx
         MNWNOD(17,INODE)=Tyy
        end do 
       end if 
      end do 
      return
      end
c
c
c_________________________________________________________________________________
c
      SUBROUTINE SMNW2COND(delr,delc,hnew,hold,hstrt,SC1,
     +                 ncol,nrow,nlay,kkstp,ssflag,kkper,
     &                 BOTM,NBOTM,ibound,MNWINT,INTTOT,
     &                 MNWNOD,NODTOT,nmnw2,MNWMAX,MNW2,CV,IOUT,ITFLAG,
     &                 INBCF,LAYHDT,WELLID,kiter,HK,VKA,MNWPRNT,NMNWVL)
C     VERSION 20060704 KJH
c
c----- MNW1 by K.J. Halford
c----- MNW2 by G.Z. Hornberger
c
c     ******************************************************************
c     Calculate all Cell-to-well conductance terms
c     ******************************************************************
c
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*20 WELLID
      CHARACTER*9 ctext
      INTEGER nmnw2,iw,MNWMAX,firstnode,lastnode,inode,ix,iy,iz,ibound,
     & NODTOT,LBOTM,ncol,nrow,nlay,NBOTM,LAYCBD,firstint,lastint,iint,
     & INTTOT,LOSSTYPE,IOUT,ITFLAG,kkstp,INBCF,ssflag,layhdt,NNODES,
     & ISOLNFLAG,kiter,irecalc,LAYVKA,LAYTYP,LAYAVG,LAYWET,ipr,nod,
     & kkper,MNWPRNT,NMNWVL,nd
      REAL delr,delc,botm, CV,hold,hstrt,SC1,HK,VKA,CHANI
      DOUBLE PRECISION verysmall,MNW2,cond,dx,dy,top,bot,thck,
     & Txx,Tyy,MNWNOD,rw,Qact,Rskin,Kskin,B,C,CF,PLoss,cel2wel2,alpha,
     & Kz,totlength,lengthint,MNWINT,ratio,CWC,hnew,ztop,zbotm,dhp,SS,
     & ZPD,ZPL,ABC,ABCD,lengthratio,T,Kh,QQ,dpp,topscreen,bottomscreen,
     & Skin
      dimension mnw2(NMNWVL,MNWMAX),MNWNOD(33,NODTOT),MNWINT(11,INTTOT)
      dimension delr(ncol), delc(nrow), CV(NCOL,NROW,NLAY),
     & SC1(ncol,nrow,nlay)
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),IBOUND(ncol,nrow,nlay),
     & hnew(ncol,nrow,nlay),hold(ncol,nrow,nlay),hstrt(ncol,nrow,nlay),
     & LAYHDT(NLAY),HK(NCOL,NROW,NLAY),VKA(NCOL,NROW,NLAY)
      dimension WELLID(mnwmax)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
      COMMON /LPFCOM/LAYTYP(999),LAYAVG(999),CHANI(999),LAYVKA(999),
     1               LAYWET(999)
C
      verysmall = 1.0D-20
c
c1------if number of wells <= 0 then return.
      if(nmnw2.le.0) return
c
c   set print flag for well output 
c   if transient, print every TS; if steady, every SP
      ipr=0
      if(ssflag.eq.0) then
        if(kiter.eq.1) ipr=1
      else
        if(kkstp.eq.1.and.kiter.eq.1) ipr=1
      end if
c   now check mnwprnt and SPs
      if(mnwprnt.eq.0) ipr=0
      if(kkper.gt.1.and.mnwprnt.lt.2) ipr=0
c   print header for well output
c   if transient, by kiter=1 , if not, by tstep=1 

      if(ipr.eq.1) then
        write(iout,*) 
        write(iout,'(120A)') 'MNW2 Well Conductance and Screen (Open 
     &Interval) Data'
        write(iout,'(120A)') '                              M O D E L
     &  L A Y E R     W E L L  S C R E E N   Penetration    SKIN     
     &  CALCULATED'
        write(iout,'(120A)') 'WELLID        Node    CWC*    top elev   
     &bot  elev     top elev    bot elev    fraction     COEFF.
     &          B'
      end if
c   Compute cell-to-well conductance for each well node
C
c   Loop over all wells
      do iw=1,MNWMAX
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
        LOSSTYPE=INT(MNW2(3,iw))
        NNODES=INT(MNW2(2,iw))
        firstnode=MNW2(4,iw)
        lastnode=MNW2(4,iw)+ABS(NNODES)-1
        alpha=1.0
c   determine well characteristics for horizontal wells
        if(MNW2(21,iw).GT.0) then
          CALL MNW2HORIZ(LOSSTYPE,NNODES,firstnode,lastnode,MNWNOD,
     & NODTOT,LAYHDT,DELR,DELC,BOTM,NBOTM,NCOL,NROW,NLAY,HNEW,WELLID,
     & MNWMAX,IW,kkstp,kkper,ipr,alpha,INBCF,VKA,HK,iout)
        else
c   for all other wells, define CWC in node loop   
c   Loop over nodes in well
        do INODE=firstnode,lastnode
         nod=INODE-firstnode+1
         ix=MNWNOD(3,INODE)
         iy=MNWNOD(2,INODE)
         iz=MNWNOD(1,INODE)
c set flag for deciding whether to recalculate CWC (1=true)
         irecalc=1
c
c-----if the cell is inactive or specified then bypass processing.
         if(ibound(ix,iy,iz).lt.1 ) irecalc=0
c if confined (THICKNESS IS NOT HEAD-DEPENDENT), don't recalculate CWC
         if(LAYHDT(IZ).EQ.0.and.kiter.gt.1) irecalc=0        
c
c if GENERAL, always recalculate
         if(LOSSTYPE.eq.3.and.MNWNOD(9,INODE).GT.0.d0) irecalc=1
c
         if(irecalc.eq.1) then
c-----if the cell is inactive or specified then bypass processing.
c         if(ibound(ix,iy,iz).ne.0 ) then
            if(LAYHDT(IZ).EQ.0) then
c if confined (THICKNESS IS NOT HEAD-DEPENDENT), don't use hnew=top
              top=BOTM(IX,IY,LBOTM(IZ)-1)
            else
              top = hnew(ix,iy,iz)
              if(top.gt.(BOTM(IX,IY,LBOTM(IZ)-1)))
     &          top=BOTM(IX,IY,LBOTM(IZ)-1)
            end if
            bot = BOTM(IX,IY,LBOTM(IZ))
            thck = top-bot
c     Check for SPECIFIED CONDUCTANCE option (LOSSTYPE=4) for node-defined well
          if(LOSSTYPE.EQ.4.and.NNODES.GT.0) then
            cond = MNWNOD(11,INODE)
          else
            dx   = delr(ix)
            dy   = delc(iy)
            Txx = MNWNOD(16,INODE)
            Tyy = MNWNOD(17,INODE)
            Qact = MNWNOD(4,INODE)
c
c           If this is not a vertical well (i.e. NNODES>0)

            if(NNODES.GT.0) then
c
              rw = MNWNOD(5,INODE)
              Rskin = MNWNOD(6,INODE)
              Kskin = MNWNOD(7,INODE)
              B = MNWNOD(8,INODE)
              Cf = MNWNOD(9,INODE)
              PLoss = MNWNOD(10,INODE)           
c   compute conductance term for node
              cond = cel2wel2(LOSSTYPE,Txx,Tyy,dx,dy,
     &                 rw,Rskin,Kskin,B,Cf,PLoss,thck,Qact,
     &                 IOUT,WELLID(iw),Skin)
c   check cond<0, reset to 0 and print warning
              if(cond.lt.0.d0) then
                write(iout,*) '***WARNING*** CWC<0 reset to CWC=0'
                write(iout,*) 'In Well ',WELLID(iw),' Node ',INODE
                cond=0.d0
              end if
c   check cond<0, reset to 0 and print warning
              if(cond.lt.0.d0) then
                write(iout,*) '***WARNING*** CWC<0 reset to CWC=0'
                write(iout,*) 'In Well ',WELLID(iw),' Node ',INODE
                cond=0.d0
              end if
c   check PPFLAG, if on, alpha defined for each node
              if(MNW2(19,iw).GT.0) then
                alpha=MNWNOD(19,INODE)
              else
                alpha=1.0D0
              end if
            else
c   else this is a vertical well: process it
c   get first and last interval intersecting this node
              firstint=MNWNOD(12,INODE)
              lastint=MNWNOD(13,INODE)
c   initialize total length of borehole within cell
              totlength=0.0D0
c   initialize conductance; will be summed for multiple intervals
              cond=0.D0
c   initialize specified conductance; will be summed for multiple intervals
              CWC=0.D0
              do iint=firstint,lastint    
c   length of interval is ztop-zbotm     
                ztop=MNWINT(1,iint)
                zbotm=MNWINT(2,iint)
c   check boundaries/saturated thickness
                if(ztop.ge.top) ztop=top
                if(zbotm.le.bot) zbotm=bot
                if(ztop.gt.zbotm) then
                  lengthint=ztop-zbotm
                else
                  lengthint=0.D0
                endif
c   calculate total length of borehole within cell
                totlength=totlength+lengthint
                if(LOSSTYPE.EQ.4) then
                  if(totlength.gt.0D0) then
                    lengthratio=lengthint/totlength
                    CWC = CWC + lengthratio*(MNWINT(11,iint))
                  endif
                else
c   calculate weighting ratio based on full thickness of node
                  ratio=lengthint/thck
c   maximum ratio is 1.0
                  if(ratio.gt.1.d0) ratio=1.d0 
c   use length-weighted ratios for each interval to determine CWC of that interval
                  if(ratio.gt.0.D0) then
                    rw = MNWINT(5,iint)
                    Rskin = MNWINT(6,iint)
                    Kskin = MNWINT(7,iint)
                    B = MNWINT(8,iint)
                    Cf = MNWINT(9,iint)
                    Ploss = MNWINT(10,iint)
c  calculate cond, weight it by length in cell (*ratio) and sum to get effective CWC
                    cond = cond + ratio*(cel2wel2(LOSSTYPE,Txx,Tyy,dx,
     &                 dy,rw,Rskin,Kskin,B,Cf,PLoss,thck,Qact,
     &                 IOUT,WELLID(iw),Skin))
c   check cond<0, reset to 0 and print warning
              if(cond.lt.0.d0) then
                write(iout,*) '***WARNING*** CWC<0 reset to CWC=0'
                write(iout,*) 'In Well ',WELLID(iw),' Node ',INODE
                cond=0.d0
              end if
                  end if
                end if
              end do
c   calculate alpha for partial penetration effect if PPFLAG is on
              if(MNW2(19,iw).GT.0) then
                alpha=totlength/(thck)
                if(alpha.gt.0.99.and.alpha.lt.1.0) then
                  if (MNWPRNT.gt.1.and.kiter.eq.1) then
                  nd=INODE-firstnode+1
                  write(iout,*) 'Penetration fraction > 0.99 for node ',
     & nd,' of well ',wellid(iw)
                  write(iout,*) 'Value reset to 1.0 for this well'
                  end if
                  alpha=1.0
                end if
              else
                alpha=1.0
              endif
            end if  
c
c
c     Correct conductance calculation for partial penetration effect
c
c     prepare variables for partial penetration calculation
c           only do partial penetration effect if PP>0 and alpha <1.0
            IF(MNW2(19,iw).GT.0.and.alpha.lt.1.D0) then
c
c  use saved partial penetration effect if steady state and past 1st iter
              if(ssflag.gt.0.and.kiter.gt.1) then
                dhp=MNWNOD(18,INODE)
              else
c      if transient, update dhp
                T = (Txx*Tyy)**0.5D0
                Kh = T/thck
                QQ=Qact*(-1.D0)
c  for 1-layer problems, Kz is not defined; assume=Kh     
C  also, FOR BCF, must assume Kh=Kz as only input is VCONT
                if(NLAY.EQ.1.or.INBCF.gt.0) then
                  Kz=Kh
                else
C  FOR LPF==check LAYVKA and VKA
                  IF(LAYVKA(iz).EQ.0) THEN
                    Kz=VKA(ix,iy,iz)
                  ELSE
                    Kz=HK(ix,iy,iz)/VKA(ix,iy,iz)
                  END IF
                end if
c
c   calculate specific storage; if steady state, send in a mock SS value
                if(ssflag.ne.0) then
                  SS=1d-5
c
                else
                  SS=SC1(IX,IY,IZ)/(thck*dx*dy)
                end if
c
c    determine location of well screen in cell
c
c    only calculate this once for each well, then save topscreen and bottomscreen
c    topscreen (MNWNOD(20) is flagged as 1d30 until it is set
                if(MNWNOD(20,INODE).eq.1d30) then
c             if a vertical well
                 if(NNODES.LT.0) then
c    if firstint=lastint for this node, it is the only interval, so use exact
c    location of borehole for analytical calculation
                  if(firstint.eq.lastint ) then            
                    topscreen=ztop
                    bottomscreen=ztop-totlength
c    if multiple screens in a confined (constant thck) cell, assume in middle
c    (calculation: from the top, go down 1/2 the amount of "unscreened" aquifer
                  else  
                    IF(LAYHDT(IZ).EQ.0) then
                      topscreen=top-((thck-totlength)/2)
                      bottomscreen=topscreen-totlength
c    if multiple screens in an unconfined (WT) cell, assume at bottom of last screen
                    else
c    (zbotm works here as it is the last thing set in the interval loop above)
                      topscreen=zbotm+totlength 
                      bottomscreen=zbotm
                    end if
                  end if
c             save top and bottom of screen
                  MNWNOD(20,INODE)=topscreen
                  MNWNOD(21,INODE)=bottomscreen
c             else if not a vertical well
c          
                 else
c             alpha specified; calculate length of screen
                  totlength=thck*alpha
c                 if confined (constant thck), assume borehole in middle             
                  IF(LAYHDT(IZ).EQ.0) then
                    topscreen=top-((thck-totlength)/2)
                    bottomscreen=topscreen-totlength
c             if unconfined, assume borehole at bottom of cell             
                  else
                    topscreen=bot+totlength
                    bottomscreen=bot                    
                  end if                                          
c             save top and bottom of screen
                  MNWNOD(20,INODE)=topscreen
                  MNWNOD(21,INODE)=bottomscreen

                 end if
c             if topscreen and bottomscreen have been calculated, retrieve them
                else
                  topscreen=MNWNOD(20,INODE)
                  bottomscreen=MNWNOD(21,INODE)
                end if
c
c from top and bottom of screen info, calculate ZPD and ZPL for PPC routine  
                ZPD=top-topscreen
                ZPL=top-bottomscreen
c if ZPD is less that zero, the screen is at the "top", so set ZPD=0
                if(ZPD.lt.0.D0) ZPD=0.D0
c
C
c calculate dhp (Delta-H due to Penetration) using analytical solution          
c
                CALL PPC(dhp,ISOLNFLAG,thck,Kh,Kz,SS,QQ,rw,ZPD,ZPL)
c          
c  if analyitcal solution failed, report no partial penetration and set dhp=0.0
                if(ISOLNFLAG.EQ.-1.AND.ITFLAG.GT.0.and.QQ.ne.0.D0) then
c  if alpha <= 0.2, shut well off if PPC did not converge
                  if(alpha.lt.0.2) then
                   if (MNWPRNT.gt.1) then
                    nd=INODE-firstnode+1
                    write(iout,*) 'Partial penetration solution did not
     & converge; penetration fraction < 0.2,      resetting CWC= 0.0 for
     & node '
     & ,nd,' of well ',wellid(iw)
                   end if
                   cond=0.0
                  else
c  if alpha > 0.2, set PPC effect = 0 if did not converge
                   if (MNWPRNT.gt.1) then
                    nd=INODE-firstnode+1
                    write(iout,*) 'Partial penetration solution did not
     & converge; penetration fraction > 0.2,      assume full 
     & penetration for
     & node ',nd,' of well ',wellid(iw)
                   end if
                   dhp=0.0 
                  endif
                end if
                if(ISOLNFLAG.EQ.0) then
                    write(iout,*) 'Partial penetration effect negligible
     &  for node '
     & ,nd,' of well ',wellid(iw),'.  No PP effect for this node'
                end if
c  store partial penetration effect (dhp)         
                MNWNOD(18,INODE)=dhp
              end if              
c             end if recalc dhp
c
c  correct partially penetrating node-defined cells by ratio of screenlength/satthck
              if(NNODES.GT.0) then
                ratio=(topscreen-bottomscreen)/thck
                cond=cond*ratio
              end if
c             
c  re-calculate conductance to include partial penetration
c    calculate dpp (partial penetration effect with specific Q) if Q and dhp "align" correctly 
c    eg if removing water (Q<0), dhp should be positive
c    (dhp>0 signifies drawdown).  Q is either <> 0 so no div 0 problem

              if(ITFLAG.EQ.1.
     &          .AND.(Qact.lt.0.D0.AND.dhp.gt.0.D0)
     &          .OR.(Qact.gt.0.D0.AND.dhp.lt.0.D0)) then
                  dpp=dhp/(Qact*(-1.D0))
                ABC=1/cond
                ABCD=ABC+dpp
                cond=1/ABCD
cgzh add dhp.ne.0 for PPC that didn't converge
              else if (ITFLAG.EQ.1.and.Qact.ne.0.d0.and.dhp.ne.0.0) then
                dpp=0.d0
                write(iout,*) '***WARNING***  Partial penetration term
     & (dpp) set to 0.0 due to misalignment of dhp= ',dhp,' and Q=',Qact
              end if     

            END IF              
c           end if PP effect

          endif
c         endif LOSSTYP EQ 4 and NNODES GT 0
c        Save conductance of each node
          MNWNOD(14,INODE) = cond
c        endif irecalc=1
         else
c        if irecalc=0, use saved cond
          cond= MNWNOD(14,INODE)
         endif
c        output node info  
c  if more than one interval made up this node, write composite 
          if(MNWNOD(12,INODE).ne.MNWNOD(13,INODE)) then
            ctext='COMPOSITE'
          else
            ctext='         '
          end if        
c only write screen info for cells that have partial penetration
         if(ipr.eq.1) then
           if(MNW2(19,iw).GT.0.and.alpha.lt.1.0D0) then
            if(LOSSTYPE.eq.2) then
             write(iout,'(A15,I3,1P7G12.5,1PG12.4,9A)') 
     & WELLID(iw),nod,cond,
     & top,bot,topscreen,bottomscreen,alpha,Skin,B,ctext
            else 
             write(iout,'(A15,I3,1P6G12.5,12A,12A,9A)') 
     & WELLID(iw),nod,cond,
     & top,bot,topscreen,bottomscreen,alpha,'     N/A    ',
     & '     N/A    ',ctext
            end if
           else
c for no partial penetration, just repeat top and bot of layer
            if(LOSSTYPE.eq.2) then
             write(iout,'(A15,I3,1P7G12.5,1PG12.4,9A)') 
     & WELLID(iw),nod,cond,
     & top,bot,top,bot,alpha,Skin,B,ctext
            else
             write(iout,'(A15,I3,1P6G12.5,12A,12A,9A)') 
     & WELLID(iw),nod,cond,
     & top,bot,top,bot,alpha,'     N/A    ',
     & '     N/A    ',ctext
            end if
           end if
         end if
        enddo
c       enddo loop over nodes
       endif
c      endif horizontal well check
       endif
c      endif active node
      enddo
c     enddo loop over wells
c
c     write note about CWC values
c     if transient, by kiter=1 , if not, by tstep=1 
      if(ipr.eq.1) then
      write(iout,'(120A)') '* Cell-to-well conductance values (CWC) may 
     &change during the course of a stress period'
      write(iout,*) 
      end if
c
      return
      end
c
c
c
c
c_________________________________________________________________________________
c
      DOUBLE PRECISION function cel2wel2(LOSSTYPE,Txx,Tyy,dx,dy,
     &                 rw,Rskin,Kskin,B,Cf,PLoss,thck,Q,IOUT,WELLNAME,
     &                 Skin)
C
C     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
C     VERSION 20090405 GZH        -- MNW2
c
c----- MNW1 by K.J. Halford
c
c     ******************************************************************
c     Compute conductance term to define head loss from cell to wellbore
c      Methodology is described in full by Peaceman (1983)
c     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*20 WELLNAME
      INTEGER LOSSTYPE,IOUT
      DOUBLE PRECISION pi,verysmall,rw,Txx,Tyy,yx4,xy4,ro,dx,dy,Tpi2,A,
     & Ploss,B,Rskin,Kskin,C,Cf,Q,thck,T,Tskin,Skin
C     ------------------------------------------------------------------
C
      pi = 3.1415926535897932D0
      verysmall = 1.D-25

      if( rw.lt.verysmall .or. Txx.lt.verysmall .or. Tyy.lt.verysmall ) 
     &  then
        cel2wel2 = ( Txx * Tyy )** 0.5D0
      else
        yx4 = (Tyy/Txx)**0.25D0
        xy4 = (Txx/Tyy)**0.25D0
        ro = 0.28D0 *((yx4*dx)**2 +(xy4*dy)**2)**0.5D0 / (yx4+xy4)
c
        Tpi2 = 2.D0*pi * (Txx*Tyy)**0.5D0
c       if ro/rw is <1, 'A' term will be negative.  Warn user and cut off flow from this node
        if (ro/rw.lt.1.D0) then
          write(iout,*) 
     &      '     Ro/Rw =  ',Ro/Rw, 
     &      '***WARNING*** Ro/Rw < 1, CWC set = 0.0 for well ',WELLNAME
          cel2wel2 = 0.D0
          GOTO 888
        end if
        A = log(ro/rw) / Tpi2
c       For the "NONE" option, multiply the Kh by 1000 to equivalate Hnew and hwell
        if(LOSSTYPE.EQ.0) then
          cel2wel2=1.0D3*((Txx*Tyy)**0.5D0)/thck    
c
c       THEIM option (LOSSTYPE.EQ.1) only needs A, so no need to calculate  B or C
c
c       SKIN (LINEAR) option, calculate B, C=0
        elseif(LOSSTYPE.EQ.2) then
c         average T in aquifer assumed to be sqrt of Txx*Tyy
          T  = (Txx*Tyy)**0.5D0
          Tskin = Kskin*thck
          if(Tskin.gt.0.D0.and.rw.gt.0.D0) then
c         this is from eqs 3 and 5 in orig MNW report
            Skin = ((T/Tskin)-1)*(DLOG(Rskin/rw))
            B = Skin / Tpi2
          else
            B = 0.D0
          end if
          C = 0.D0
c       GENERAL option, calculate B and C
       else if (LOSSTYPE.EQ.3) then
          if(Cf.NE.0.0) then
            C = Cf * abs(Q)**(PLoss-1)
          else
            C = 0.D0
          end if
       else
          B = 0.D0
          C = 0.D0
       end if
        cel2wel2 = A + B + C 
        cel2wel2 = 1.000000D0 / cel2wel2
      endif
c
 888  return
      end
c
c
c
c_________________________________________________________________________________
c
c
c_________________________________________________________________________________
c
      subroutine SMNW2SEEP(mnw2,MNWNOD,iw,MNWMAX,NODTOT,IBOUND,
     & ncol,nrow,nlay,hnew,kiter,BOTM,NBOTM,kkstp,
     & kkper,small,hdry,NMNWVL)
c     
      IMPLICIT NONE
      INTEGER kSeep,firstnode,lastnode,iw,mnwmax,INODE,il,ir,ic,NODTOT,
     & IBOUND,ncol,nrow,nlay,NBOTM,LBOTM,ipole,kiter,kkper,nnodes,
     & LAYCBD,NMNWVL,kkstp
      REAL BOTM,hdry
      DOUBLE PRECISION mnw2,qdes,qact,csum,chsum,Qseep,MNWNOD,hnew,
     & hwell,hlim,verysmall,hmax,hsim,bottom,qcut,qoff,qon,qsmall,small,
     & qpot,cond,ratio
      DIMENSION BOTM(NCOL,NROW,0:NBOTM),mnw2(NMNWVL,MNWMAX),
     & MNWNOD(33,NODTOT),IBOUND(NCOL,NROW,NLAY),hnew(NCOL,NROW,NLAY)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
C
      verysmall = 1.D-25
      Qsmall = verysmall
c
c QFM  rechecking Q cutoff within FM routine
c skip this whole thing it 1st or 2nd iter
      if(kiter.gt.2) then
c   Retrieve Qfrcmn, Qfrcmx, Qdes
         qdes = mnw2(5,iw)
         if(mnw2(6,iw).gt.0) then
          QCUT = MNW2(8,iw)
          if(QCUT.ne.0) then
           qoff = mnw2(9,iw)
           qon  = mnw2(10,iw)
           if(QCUT.GT.0) then
c     convert rate into fraction of Qdes (ratio is used to compare)
            if(Qdes.NE.0) then
              qoff=abs(qoff/Qdes)
              qon=abs(qon/Qdes)
            else
              qoff=0.0
              qon=0.0
            end if
           end if
          end if
c
c   Compute hwell / Qpot for multi-node well (not single-cell wells)
c   
          NNODES=abs(MNW2(2,iw))
            csum = 0.000D0
            chsum = 0.000D0
            qact = 0.0000D0
            Qsmall = small*abs(qdes)
c   Loop over nodes in well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+NNODES-1
            hwell = mnw2(17,iw)
            do INODE=firstnode,lastnode
              il=MNWNOD(1,INODE)              
              ir=MNWNOD(2,INODE)              
              ic=MNWNOD(3,INODE)              
              if(IBOUND(ic,ir,il).ne.0) then
                csum  = csum  + MNWNOD(14,INODE)
                chsum = chsum + MNWNOD(14,INODE)*hnew(ic,ir,il)
                qact  = qact  + MNWNOD(4,INODE)
              else
                qact  = 0.0000D0
              end if            
            end do
c---div0 ---  CSUM could go to zero if the entire well is dry
            if( csum .gt. verysmall ) then
              hwell = ( qdes + chsum ) / csum
            else
              hwell = hnew(ic,ir,il)
            endif
c      Test Hlim constraint if QLIMIT flag is set
            if(MNW2(6,iw).GT.0) then
              hlim = mnw2(7,iw)
              ipole = 0
              if( abs(qdes).gt.verysmall ) ipole = qdes / abs(qdes)
              hmax = ipole*( hlim )
              hsim = ipole*( hwell )
c      Potential Q is...
              if( hsim .gt. hmax ) then
                hwell = hlim
              endif
              qpot = hlim*csum - chsum
            end if
            cond = csum
c    
c  Compute ratio of potential/desired flow rates
            if(QCUT.ne.0) then
             ratio = 1.00D0
             if( abs(qdes) .gt. small ) ratio =  qpot / qdes
             if( ratio .gt. 0.9999D0 ) then
              ratio =  1.000D0
              Qpot = Qdes
             endif
c  Check if potential flow rate is below cutoff
c  If so, set Qact = 0.0
c  If q-limit condition is met, do not perform this check on subsequent iterations to
c    avoid oscillations
c  Flag to perform this check is mnw2(20,iw).  If > 0, limit has been met
c  Do not enforce q-limit oscillation-stopper until the 3rd iteration; hardwire=0
c    for first two iterations
             if(kiter.le.2) mnw2(20,iw) = 0
c
             if( ratio .lt. Qoff . and. mnw2(20,iw) .ge. 0) then
c              mnw2(18,iw) = 0.D0
              mnw2(30,iw) = 0.D0
c  Set flag to avoid rechecking q-limit
              mnw2(20,iw) = -1
              
c  Check if potential flow rate is above restart threshold
c  If so, set Qact = Qpot
c  If q-limit condition is met, do not perform this check on subsequent iterations to
c    avoid oscillations
c  Flag to perform this check is mnw2(20,iw).  If > 0, limit has been met
             elseif( ratio.gt.Qon .and. abs(qact).lt.Qsmall .and. 
     &              mnw2(20,iw) .le. 0 ) then
c              mnw2(18,iw) = Qpot
              mnw2(30,iw) = Qpot
c  Set flag to avoid rechecking q-limit
              mnw2(20,iw) = 1
             elseif( ratio.lt.Qon.and.ratio.gt.Qoff ) then
              
              mnw2(20,iw) = 2
             else
c  Otherwise leave the flow rate alone
             endif
            endif
c  End if, QLimit>0
         endif
c  End if, kiter < 3
      end if
c
c
      qdes = mnw2(5,iw)
      qact = qdes
c   If Q in last TS was restricted, update qact
      if((abs(mnw2(18,iw)-qdes).gt.Qsmall))  then
        qact = mnw2(18,iw)
      end if
c   If restrictions were set this time step, update qact
      if((mnw2(20,iw).ne.0.or.mnw2(27,iw).ne.0)) then
        if((abs(mnw2(29,iw)).lt.abs(mnw2(18,iw)))) then
           if(mnw2(27,iw).ne.0) then
             mnw2(18,iw)=mnw2(29,iw)
           else
             mnw2(18,iw)=mnw2(30,iw)
           end if
        else
           if(mnw2(27,iw).eq.0) then
             mnw2(18,iw)=mnw2(30,iw)
           else
             mnw2(18,iw)=mnw2(29,iw)
           end if
        end if
        qact = mnw2(18,iw)
      end if
c
c  Make 2 passes to find seepage faces
      do kSeep = 1, 2          
        csum = 0.000D0
        chsum = 0.000D0
        Qseep = 0.000D0
        firstnode=MNW2(4,iw)
        lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
c  Loop over nodes in well
        do INODE=firstnode,lastnode
            il=MNWNOD(1,INODE)              
            ir=MNWNOD(2,INODE)              
            ic=MNWNOD(3,INODE)              
c  First time through, set the hwell_by_node to big; this will be used as a flag for 
c  whether or not there is a seepage face in the cell
            if( kSeep.eq.1 ) then
              MNWNOD(15,INODE) = 1.0D31
            end if
            if(IBOUND(ic,ir,il).ne.0) then
              Bottom = BOTM(ic,ir,LBOTM(il))
c
c  only allow seepage or flow into a cell if the head in the cell is greater than the bottom of the cell
c
             if(hnew(ic,ir,il).gt.bottom) then
c  Second time through loop, now that we have a guess for hwell, check to see if there is a seepage face
c  If so, use the bottom (instead of the low hwell) as the gradient driver in the qact calc
c  Set the hwell_by_node to the bottom to flag this as a seepage face cell
c  Sum the seepage (there may be more than one node's worth
              if( kSeep.gt.1 .and. hwell.lt.Bottom )then
                MNWNOD(4,INODE) = (bottom - hnew(ic,ir,il))
     & * MNWNOD(14,INODE)
                MNWNOD(15,INODE) = bottom
                Qseep = Qseep + MNWNOD(4,INODE)
              else 
                csum  = csum  + MNWNOD(14,INODE)
                chsum = chsum + MNWNOD(14,INODE)*hnew(ic,ir,il)
              endif
             else
              MNWNOD(4,INODE) = 0.0D0
              MNWNOD(14,INODE) = 0.0D0
              MNWNOD(15,INODE) = hdry
             end if
            else
             MNWNOD(4,INODE) = 0.0D0
             MNWNOD(15,INODE) = hdry
            end if
c          if(kkper.eq.2.and.inode.eq.1) write(iout,*) 
c     &'seep1: hnew, bottom, mnwnod4',hnew(ic,ir,il),
c     &BOTM(ic,ir,LBOTM(il)),MNWNOD(4,INODE)
        end do           
c  End loop over nodes in well
c---div0 ---  CSUM could go to verysmall if the entire well is dry
        if( csum .gt. verysmall ) then
          hwell = ( qact - Qseep + chsum ) / csum
        else
          hwell = hnew(ic,ir,il)
        endif
c   Because the q/hwell may now be different due to seepage flow, need to re-check constraints
c   Test Hlim constraints, Hlim is assumed to be a Max/Min for Injection/Production wells
        if(MNW2(6,iw).GT.0) then
          hlim = mnw2(7,iw)
          ipole = 0
          if( abs(qdes).gt.verysmall ) ipole = qdes / abs(qdes)
          hmax = ipole*( hlim )
          hsim = ipole*( hwell )
c
          if( hsim .gt. hmax ) then
            hwell = hlim
            qact = hwell*csum - chsum + Qseep
c      Hlim constraints that stop production are not tested until after the 2nd iteration
            if( kiter.gt.2 ) then
              ratio = 1.00D0
              if( abs(qdes) .gt. small ) ratio =  qact / qdes
              if( ratio .lt. 0.00001D0 ) then
                qact  = 0.000D0
                if (csum .gt. 0.0D0) then
                  hwell = ( qact - Qseep + chsum ) / csum
                else
                  hwell = hnew(ic,ir,il)
                endif
              endif  
            endif   !! potentially stop all production after Kiter = 2 
          endif     !! End Hlim exceedence loop
        endif  !! Qlimit>0
      enddo  !!  kSeep -- End of seepage face corrector here
c
c  Loop over nodes in well, assign flow rates and water levels
c  use qdes to sum q's at end
      qdes=0
      do INODE=firstnode,lastnode
        il=MNWNOD(1,INODE)              
        ir=MNWNOD(2,INODE)              
        ic=MNWNOD(3,INODE)              
c  Qseep flag (MNWNOD(15): if still 1E31 here, there is no seepage face so use actual hwell to calc Q
c (if not set here, hwell < bottom and MNWNOD(15) stores value of bottom to calculate seepage, see above)
        if(MNWNOD(15,INODE) .gt. 1.0E30 )then
          qact = ( hwell - hnew(ic,ir,il) ) * MNWNOD(14,INODE)
          MNWNOD(4,INODE) = qact
          MNWNOD(15,INODE) = hwell
          qdes=qdes+qact
        endif
      end do           
c  Set hwell 
      MNW2(17,iw) = hwell
c                           
c          
      return      
      end
c
c_________________________________________________________________________________
c
c
c_________________________________________________________________________________
c
      SUBROUTINE GWF1MNW2BH(MNWMAX,mnw2,NODTOT,MNWNOD,iout,NMNWVL,iw)
C
c     ******************************************************************
c     compute borehole flow in mnw well 
c     ******************************************************************
c
      IMPLICIT NONE
      integer iw,MNWMAX,NODTOT,firstnode,lastnode,INODE,PUMPLOC,
     & nodepump,ifound,ir,irp,il,ilp,ic,icp,iout,NMNWVL
      double precision mnw2,Qnet,MNWNOD,diff,q,q27
      dimension MNWNOD(33,NODTOT),mnw2(NMNWVL,MNWMAX)
c   QBH (MNWNOD(27,m) is flow between nodes (saved at upper face) in borehole
c   Only operate on active wells (MNW2(1,iw)=1)
       if (MNW2(1,iw).EQ.1) then
        firstnode=MNW2(4,iw)
        lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
        PUMPLOC=MNW2(11,iw)
c
        if(PUMPLOC.eq.0) then
         nodepump=firstnode
        else
          do inode=firstnode,lastnode
C Initialize
          ifound=0
	    MNWNOD(27,INODE)=0.D0
            MNWNOD(33,INODE)=0
	    
c get node coordinates
            il=MNWNOD(1,INODE)              
            ir=MNWNOD(2,INODE)              
            ic=MNWNOD(3,INODE)              
c get pump coordinates
            ilp=MNW2(14,iw)              
            irp=MNW2(15,iw)              
            icp=MNW2(16,iw)              
            if(il.eq.ilp.and.ir.eq.irp.and.ic.eq.icp) then
              ifound=1
              nodepump=inode 
c set pump node flag
              MNWNOD(33,inode)=1
            end if   
          end do
          if(ifound.eq.0) then
            write(iout,*) '***ERROR*** Pump location specified but 
     & not found, MNW2'
            STOP 'MNW2 ERROR - PUMPLOC2'
          end if
        end if
c get qnet
c        Qnet=0.D0
c        do INODE=firstnode,lastnode
c          q=MNWNOD(4,INODE)
c          Qnet=Qnet+q
c        end do
       Qnet=mnw2(18,iw)          
c   Set flux at top of first node equal to Qnet (if pump location not specified)
c     or set to zero if pump is specified 
c
        q27=0.0
        if(PUMPLOC.eq.0) then
          MNWNOD(27,firstnode)=-Qnet
          q27=-qnet
        else
          MNWNOD(27,firstnode)=0.d0
          if(nodepump.eq.firstnode) q27=-qnet
        end if
c
c   Loop over nodes in well
        do INODE=firstnode+1,lastnode
c   Loop over other nodes in this MNW to set QBH
c   QBH between successive nodes is QBH at previous node - Qaq at node
c   If pump location is specified, apply Qnet at that node too
          if(PUMPLOC.NE.0.and.nodepump.eq.inode) then
            MNWNOD(27,INODE)=q27+MNWNOD(4,INODE-1)-Qnet
          else          
            MNWNOD(27,INODE)=q27+MNWNOD(4,INODE-1)
          end if
          q27=MNWNOD(27,INODE)
        enddo
       end if
C
      RETURN
C
      END
C
C
C MNW2HORIZ
C
C Process CWC for horizontal MNWs
C
C     ******************************************************************
C
      SUBROUTINE MNW2HORIZ(LOSSTYPE,NNODES,firstnode,lastnode,MNWNOD,
     & NODTOT,LAYHDT,DELR,DELC,BOTM,NBOTM,NCOL,NROW,NLAY,HNEW,WELLID,
     & MNWMAX,IW,kkstp,kkper,ipr,alpha,INBCF,VKA,HK,iout)
C
C     ******************************************************************
      IMPLICIT NONE
      ALLOCATABLE ivert1(:),ivert2(:)
      INTEGER Wel1flag,QSUMflag,BYNDflag,nstp,kstp,nmnw2,INBCF
      CHARACTER*20 WELLID
      INTEGER L1,R1,C1,L2,R2,C2,L,R,C,LAYHDT,LAYCBD,LOSSTYPE,NNODES,
     & firstnode,lastnode,NODTOT,is_intersection,INODE,NCOL,NROW,NLAY,
     & NBOTM,LBOTM,kkstp,kkper,ivert1,ivert2,idone,iout,mnwmax,IW,nod,
     & ipr
      INTEGER LAYTYP,LAYAVG,LAYVKA,LAYWET
      REAL DELR,DELC,BOTM,CHANI,
     & x1face1,x1face2,y1face1,y1face2,z1face1,z1face2,
     & x1face,y1face,z1face,
     & x2face1,x2face2,y2face1,y2face2,z2face1,z2face2,
     & x2face,y2face,z2face,
     & m,lxf,lyf,lzf,lbf,
     & zwt,ywt,xwt,
     & zi,yi,xi,zi2,yi2,xi2,HK,VKA
      DOUBLE PRECISION z1,y1,x1,z2,y2,x2,top1,bot1,top2,bot2,
     & betweennodes,omega_opp,omega,theta_opp,theta_hyp,
     & theta,thck1,thck2,lw,cel2wel2SEG,dx1,dx2,dy1,dy2,
     & cel2wel2,alpha,T,Kh,Kz,Txx1,Tyy1,Skin
      DOUBLE PRECISION MNWNOD,hnew(ncol,nrow,nlay),
     & Txx,Tyy,rw,Rskin,Kskin,B,Cf,PLoss,Qact,cond1,cond2,cond
      DOUBLE PRECISION dgr_to_rad,pi
      dimension WELLID(mnwmax)      
      DIMENSION mnwnod(33,NODTOT)
      DIMENSION LAYHDT(NLAY),delr(ncol),delc(nrow),
     & BOTM(NCOL,NROW,0:NBOTM),HK(NCOL,NROW,NLAY),VKA(NCOL,NROW,NLAY)
      COMMON /DISCOM/LBOTM(999),LAYCBD(999)
      COMMON /LPFCOM/LAYTYP(999),LAYAVG(999),CHANI(999),LAYVKA(999),
     1               LAYWET(999)
      ALLOCATE(ivert1(NODTOT),ivert2(NODTOT))
c convert degree trig func modified from http://techpubs.sgi.com
      pi = 3.1415926535897932D0
      dgr_to_rad = (pi/180.D0)
c     compute borehole length and screen orientation
c
c     compute length associated with each section 
c
c     compute CWC for each node
c
c   Initialize flags
      ivert1=0
      ivert2=0
c   Loop over "segments"
        do INODE=firstnode,lastnode-1
         nod=INODE-firstnode+1
c   Initialize flags
         is_intersection=0
c   Define node and next node
         L1=MNWNOD(1,INODE)              
         R1=MNWNOD(2,INODE)              
         C1=MNWNOD(3,INODE)              
         L2=MNWNOD(1,INODE+1)              
         R2=MNWNOD(2,INODE+1)              
         C2=MNWNOD(3,INODE+1)
         dx1=DELR(C1)             
         dx2=DELR(C2)             
         dy1=DELC(R1)             
         dy2=DELC(R2)             
C     convert to real coodinates
         x1=0 
         do C=1,C1-1
           x1=x1+DELR(C)
         end do
         x1=x1+0.5D0*DELR(C1)
         x2=0 
         do C=1,C2-1
           x2=x2+DELR(C)
         end do
         x2=x2+0.5D0*DELR(C2)
         y1=0 
         do R=1,R1-1
           y1=y1+DELC(R)
         end do
         y1=y1+0.5D0*DELC(R1)
         y2=0 
         do R=1,R2-1
           y2=y2+DELC(R)
         end do
         y2=y2+0.5D0*DELC(R2)
          if(LAYHDT(L1).EQ.0) then
c if confined (THICKNESS IS NOT HEAD-DEPENDENT), don't use hnew=top
           top1=BOTM(C1,R1,LBOTM(L1)-1)
          else
           top1 = hnew(C1,R1,L1)
           if(top1.gt.(BOTM(C1,R1,LBOTM(L1)-1)))
     &       top1=BOTM(C1,R1,LBOTM(L1)-1)
          end if
          bot1 = BOTM(C1,R1,LBOTM(L1))
          thck1 = (top1-bot1)/2.d0
          z1 = 0.5D0*(top1+bot1)
c
          if(LAYHDT(L2).EQ.0) then
c if confined (THICKNESS IS NOT HEAD-DEPENDENT), don't use hnew=top
           top2=BOTM(C2,R2,LBOTM(L2)-1)
          else
           top2 = hnew(C2,R2,L2)
           if(top2.gt.(BOTM(C2,R2,LBOTM(L2)-1)))
     &      top2=BOTM(C2,R2,LBOTM(L2)-1)
          end if
          bot2 = BOTM(C2,R2,LBOTM(L2))
          thck2 = (top2-bot2)/2.d0
          z2 = 0.5D0*(top2+bot2)
c   save z coords as we don't want z screen elevations to change for WT cases
c
         if(kkstp.eq.1.and.kkper.eq.1) then
           MNWNOD(26,INODE)=z1
           MNWNOD(26,INODE+1)=z2
         else
           z1=MNWNOD(26,INODE)
           z2=MNWNOD(26,INODE+1)
         end if
c     caculate distance between nodes
      betweennodes=SQRT(((x1-x2)**2)+((y1-y2)**2)+((z1-z2)**2))
c
c
c   estimate length of borehole segments
c
c   in first node, use vertical section up to top or WT for segment 1 
      if(INODE.eq.1) then
        MNWNOD(23,INODE)=0.D0
        if(z1.lt.top1) MNWNOD(23,INODE)=top1-z1   
        ivert1(INODE)=1
      end if
c   if this is a vertical segment, define lengths with elevations, skip other calc
      if(x1.eq.x2.and.y1.eq.y2) then   
        if(top1.le.bot1) MNWNOD(24,INODE)=0.D0         
        if(z1.gt.top1) then
          MNWNOD(24,INODE)=top1-bot1  
        else
          MNWNOD(24,INODE)=z1-bot1  
        endif       
        MNWNOD(23,INODE+1)=top2-z2         
c  if blank spaces inbetween, save that length
        if(bot1.ne.top2) MNWNOD(25,INODE)=bot1-top2
        ivert2(INODE)=1
        ivert1(INODE+1)=1
c
c     if not vertical, calculate theta and omega for segment
c
      else        
c
      if(z1.eq.z2) then
        omega=0.d0
      else  
        omega_opp=SQRT(((x1-x2)**2)+((y1-y2)**2))
        omega=DASIN((dgr_to_rad * omega_opp)/betweennodes)
      end if
      MNWNOD(28,INODE)=omega
c
      theta_opp=abs(y2-y1)
      theta_hyp=SQRT(((x1-x2)**2)+((y1-y2)**2))
      theta=DASIN((dgr_to_rad * theta_opp)/(dgr_to_rad * theta_hyp))
c     correct for right quadrant
      if(y2.ge.y1) then
        if(x2.ge.x1) then
          theta=360.D0-theta
        else if (x2.le.x1) then
          theta=270.D0-theta
        end if
      else if (y2.le.y1) then
        if (x2.le.x1) then
          theta=180.D0-theta
        end if        
      end if      
      MNWNOD(29,INODE)=theta
c   define first cell's limits to test for first intersection
c   only for non-vertical sections
          x1face1=x1-0.5D0*DELR(C1)
          x1face2=x1+0.5D0*DELR(C1)
          y1face1=y1-0.5D0*DELC(R1)
          y1face2=y1+0.5D0*DELC(R1)
        z1face1=z1-0.5D0*(BOTM(C1,R1,LBOTM(L1)-1)-BOTM(C1,R1,LBOTM(L1)))
        z1face2=z1+0.5D0*(BOTM(C1,R1,LBOTM(L1)-1)-BOTM(C1,R1,LBOTM(L1)))
c   define possible face of intersection in x direction, first cell
          if(x2.gt.x1) then
            x1face=x1face2
          else if(x2.lt.x1) then
            x1face=x1face1
          else
            x1face=0
          end if           
c   define possible face of intersection in y direction, first cell
          if(y2.gt.y1) then
            y1face=y1face2
          else if(y2.lt.y1) then
            y1face=y1face1
          else
            y1face=0
          end if           
c   define possible face of intersection in z direction, first cell
          if(z2.gt.z1) then
            z1face=z1face2
          else if(z2.lt.z1) then
            z1face=z1face1
          else
            z1face=0
          end if           
c   define second cell's limits to test for last intersection
          x2face1=x2-0.5D0*DELR(C2)
          x2face2=x2+0.5D0*DELR(C2)
          y2face1=y2-0.5D0*DELC(R2)
          y2face2=y2+0.5D0*DELC(R2)
        z2face1=z2-0.5D0*(BOTM(C2,R2,LBOTM(L2)-1)-BOTM(C2,R2,LBOTM(L2)))
        z2face2=z2+0.5D0*(BOTM(C2,R2,LBOTM(L2)-1)-BOTM(C2,R2,LBOTM(L2)))
c   define possible face of intersection in x direction, second cell
          if(x2.gt.x1) then
            x2face=x2face1
          else if(x2.lt.x1) then
            x2face=x2face2
          else
            x2face=0
          end if           
c   define possible face of intersection in y direction, second cell
          if(y2.gt.y1) then
            y2face=y2face1
          else if(y2.lt.y1) then
            y2face=y2face2
          else
            y2face=0
          end if           
c   define possible face of intersection in z direction, second cell
          if(z2.gt.z1) then
            z2face=z2face1
          else if(z2.lt.z1) then
            z2face=z2face2
          else        
            z2face=0
          end if           
c
c   if 1st z-coord is greater than the WT, start from intersection with WT
          if(z1.gt.HNEW(C1,R1,L1)) then
            zwt=HNEW(C1,R1,L1)                
c   at wt face, determine intersection with line segment
            m=(zwt-z1)/(z2-z1)
            xwt=x1+m*(x2-x1)             
            ywt=y1+m*(y2-y1)
c   redefine 1st point
            x1=xwt
            y1=ywt
            z1=zwt
          end if
c   at x face, determine intersection with line segment
c   xi=intersection point for x face
c   m is "slope" in parameterization of 3d line segment,
c     define m for known x (at the face) and then use that m to solve for
c     other coordinates to give point of intersection
c   xi=x1 + (x2-x1)*m
c   xi-x1/(x2-x1)=m  
c
          if(x1face.ne.0) then
          is_intersection=0
          idone=0
          m=(x1face-x1)/(x2-x1)
          yi=y1+m*(y2-y1)             
          zi=z1+m*(z2-z1)
          if(yi.ge.y1face1.and.yi.le.y1face2.and.
     &       zi.ge.z1face1.and.zi.le.z1face2) then
c       if x1face intersection point lies within cell, this is exit point
             xi=x1face
             lxf=SQRT(((x1-xi)**2)+((y1-yi)**2)+((z1-zi)**2))     
             MNWNOD(24,INODE)=lxf
             is_intersection=1
          end if
c       if exit point is on boundary with second cell, done with both segments
          if(is_intersection.eq.1) then
            if(x2face.eq.xi) then
              lxf=SQRT(((x2-xi)**2)+((y2-yi)**2)+((z2-zi)**2))         
              MNWNOD(23,INODE+1)=lxf
              idone=1
            end if
          end if
          else
            lxf=0.d0
          end if
c
c   at y face, determine intersection with line segment
          if(y1face.ne.0) then
          m=(y1face-y1)/(y2-y1)
          xi=x1+m*(x2-x1)             
          zi=z1+m*(z2-z1)
          if(xi.ge.x1face1.and.xi.le.x1face2.and.
     &       zi.ge.z1face1.and.zi.le.z1face2) then
c       if yface intersection point lies within cell, this is exit point
             yi=y1face
             lyf=SQRT(((x1-xi)**2)+((y1-yi)**2)+((z1-zi)**2))     
             MNWNOD(24,INODE)=lyf
             is_intersection=1
          end if
c       if exit point is on boundary with second cell, done with both segments
          if(is_intersection.eq.1) then
            if(y2face.eq.yi) then
              lyf=SQRT(((x2-xi)**2)+((y2-yi)**2)+((z2-zi)**2))         
              MNWNOD(23,INODE+1)=lyf
              idone=1
            end if
          end if
          else
            lyf=0.d0
          end if
c
c   at z face, determine intersection with line segment
          if(z1face.ne.0) then
          m=(z1face-z1)/(z2-z1)
          xi=x1+m*(x2-x1)             
          yi=y1+m*(y2-y1)
          if(xi.ge.x1face1.and.xi.le.x1face2.and.
     &       yi.ge.y1face1.and.yi.le.y1face2) then
c       if zface intersection point lies within cell, this is exit point
             zi=z1face
             lzf=SQRT(((x1-xi)**2)+((y1-yi)**2)+((z1-zi)**2))     
             MNWNOD(24,INODE)=lzf             
             is_intersection=1
          end if
c       if exit point is on boundary with second cell, done with both segments
          if(is_intersection.eq.1) then
            if(z2face.eq.zi) then
              lzf=SQRT(((x2-xi)**2)+((y2-yi)**2)+((z2-zi)**2))         
              MNWNOD(23,INODE+1)=lzf
              idone=1
            end if
          end if
          else
            lzf=0.d0
          end if
c   if idone still=0, then there are blank spaces inbetween nodes.  Calculate
c   length of that segment by getting intersection out of last node
          if(idone.eq.0) then
            is_intersection=0
c   at x face, determine intersection with line segment
            if(x2face.ne.0) then
            m=(x2face-x2)/(x2-x1)
            yi2=y2+m*(y2-y1)             
            zi2=z2+m*(z2-z1)
            if(yi2.ge.y2face1.and.yi2.le.y2face2.and.
     &       zi2.ge.z2face1.and.zi2.le.z2face2) then
c       if x2face intersection point lies within cell, this is exit point
             xi2=x2face 
             lxf=SQRT(((x2-xi2)**2)+((y2-yi2)**2)+((z2-zi2)**2))     
             MNWNOD(23,INODE+1)=lxf
             is_intersection=1
            end if
            else
             lxf=0.d0
            end if
c   at y face, determine intersection with line segment
            if(y2face.ne.0) then
            m=(y2face-y2)/(y2-y1)
            xi2=x2+m*(x2-x1)             
            zi2=z2+m*(z2-z1)
            if(xi2.ge.x2face1.and.xi2.le.x2face2.and.
     &       zi2.ge.z2face1.and.zi2.le.z2face2) then
c       if y2face intersection point lies within cell, this is exit point
             yi2=y2face 
             lyf=SQRT(((x2-xi2)**2)+((y2-yi2)**2)+((z2-zi2)**2))     
             MNWNOD(23,INODE+1)=lyf
             is_intersection=1
            end if
            else
             lyf=0.d0
            end if
c   at z face, determine intersection with line segment
            if(z2face.ne.0) then
            m=(z2face-z2)/(z2-z1)
            xi2=x2+m*(x2-x1)             
            yi2=y2+m*(y2-y1)
            if(xi2.ge.x2face1.and.xi2.le.x2face2.and.
     &       yi2.ge.y2face1.and.yi2.le.y2face2) then
c       if z2face intersection point lies within cell, this is exit point
             zi2=z2face 
             lzf=SQRT(((x2-xi2)**2)+((y2-yi2)**2)+((z2-zi2)**2))     
             MNWNOD(23,INODE+1)=lzf
             is_intersection=1
            end if
            else
             lzf=0.d0
            end if
c  now that we have both node exit intersection points, blank distance is betweem
c  them.  Save in MNWNOD(25) of the first node between them            
            lbf=SQRT(((xi-xi2)**2)+((yi-yi2)**2)+((zi-zi2)**2))     
            MNWNOD(25,INODE)=lbf 
          end if
c
c    for last segment, continue the line to the exit intersection of the last cell
c   define possible face of intersection in x direction, second cell
         if(INODE.eq.(lastnode-1)) then
          if(x2.gt.x1) then
            x2face=x2face2
          else if(x2.lt.x1) then
            x2face=x2face1
          else
            x2face=0
          end if           
c   define possible face of intersection in y direction, second cell
          if(y2.gt.y1) then
            y2face=y2face2
          else if(y2.lt.y1) then
            y2face=y2face1
          else
            y2face=0
          end if           
c   define possible face of intersection in z direction, second cell
          if(z2.gt.z1) then
            z2face=z2face2
          else if(z2.lt.z1) then
            z2face=z2face1
          else        
            z2face=0
          end if  
          if(x2face.ne.0) then
          m=(x2face-x2)/(x2-x1)
          yi=y2+m*(y2-y1)             
          zi=z2+m*(z2-z1)
          if(yi.ge.y2face1.and.yi.le.y2face2.and.
     &       zi.ge.z2face1.and.zi.le.z2face2) then
c       if x2face intersection point lies within cell, this is exit point
             xi=x2face
             lxf=SQRT(((x2-xi)**2)+((y2-yi)**2)+((z2-zi)**2))     
             MNWNOD(24,INODE+1)=lxf
          end if
          else
            lxf=0.d0
          end if
c
c   at y face, determine intersection with line segment
          if(y2face.ne.0) then
          m=(y2face-y2)/(y2-y1)
          xi=x2+m*(x2-x1)             
          zi=z2+m*(z2-z1)
          if(xi.ge.x2face1.and.xi.le.x2face2.and.
     &       zi.ge.z2face1.and.zi.le.z2face2) then
c       if yface intersection point lies within cell, this is exit point
             yi=y2face
             lyf=SQRT(((x2-xi)**2)+((y2-yi)**2)+((z2-zi)**2))     
             MNWNOD(24,INODE+1)=lyf
          end if
          else
            lyf=0.d0
          end if
c
c   at z face, determine intersection with line segment
          if(z2face.ne.0) then
          m=(z2face-z2)/(z2-z1)
          xi=x2+m*(x2-x1)             
          yi=y2+m*(y2-y1)
          if(xi.ge.x2face1.and.xi.le.x2face2.and.
     &       yi.ge.y2face1.and.yi.le.y2face2) then
c       if zface intersection point lies within cell, this is exit point
             zi=z2face
             lzf=SQRT(((x2-xi)**2)+((y2-yi)**2)+((z2-zi)**2))     
             MNWNOD(24,INODE+1)=lyf             
          end if
          else
            lzf=0.d0
          end if
        end if
      end if
c     
      lw=MNWNOD(23,INODE)

      if (lw.gt.0.D0) then
          Txx = MNWNOD(16,INODE)
          Tyy = MNWNOD(17,INODE)
          Txx1 = Txx*0.5d0
          Tyy1 = Tyy*0.5d0
          rw  = MNWNOD(5,INODE)
          Rskin = MNWNOD(6,INODE)  
          Kskin = MNWNOD(7,INODE)
          B = MNWNOD(8,INODE)
          Cf = MNWNOD(9,INODE)
          PLoss = MNWNOD(10,INODE)
          Qact = MNWNOD(4,INODE)
c   compute conductance term for segment
          if(ivert1(INODE).eq.0) then
c  for 1-layer problems, Kz is not defined; assume=Kh     
C  also, FOR BCF, must assume Kh=Kz as only input is VCONT
                if(NLAY.EQ.1.or.INBCF.gt.0) then
                 T = (Txx*Tyy)**0.5D0
                 Kh = T/thck1
                 Kz=Kh
                else
C  FOR LPF==check LAYVKA and VKA
                  IF(LAYVKA(L1).EQ.0) THEN
                    Kz=VKA(C1,R1,L1)
                  ELSE
                    Kz=HK(C1,R1,L1)/VKA(C1,R1,L1)
                  END IF
                end if
             cond1 = cel2wel2SEG(lw,theta,omega,LOSSTYPE,
     &         Txx,Tyy,dx1,dy1,rw,Rskin,Kskin,B,Cf,PLoss,thck1,Qact,
     &         IOUT,WELLID(iw),Kz)
c   if a vertical segment, use original function
          else
              cond1 = cel2wel2(LOSSTYPE,Txx1,Tyy1,dx1,dy1,
     &                 rw,Rskin,Kskin,B,Cf,PLoss,thck1,Qact,
     &                 IOUT,WELLID(iw),Skin)
         end if
      else
         cond1=0.D0 
      end if
      MNWNOD(30,INODE)=cond1
c calculate CWC of second segment of node
      lw=MNWNOD(24,INODE)
      if(lw.gt.0D0) then                      
          Txx = MNWNOD(16,INODE)
          Tyy = MNWNOD(17,INODE)
          Txx1 = Txx*0.5d0
          Tyy1 = Tyy*0.5d0
          rw  = MNWNOD(5,INODE)
          Rskin = MNWNOD(6,INODE)  
          Kskin = MNWNOD(7,INODE)
          B = MNWNOD(8,INODE)
          Cf = MNWNOD(9,INODE)
          PLoss = MNWNOD(10,INODE)
          Qact = MNWNOD(4,INODE)
c   compute conductance term for segment
          if(ivert2(INODE).eq.0) then
c  for 1-layer problems, Kz is not defined; assume=Kh     
C  also, FOR BCF, must assume Kh=Kz as only input is VCONT
                if(NLAY.EQ.1.or.INBCF.gt.0) then
                 T = (Txx*Tyy)**0.5D0
                 Kh = T/thck1
                 Kz=Kh
                else
C  FOR LPF==check LAYVKA and VKA
                  IF(LAYVKA(L1).EQ.0) THEN
                    Kz=VKA(C1,R1,L1)
                  ELSE
                    Kz=HK(C1,R1,L1)/VKA(C1,R1,L1)
                  END IF
                end if
             cond2 = cel2wel2SEG(lw,theta,omega,LOSSTYPE,
     &         Txx,Tyy,dx1,dy1,rw,Rskin,Kskin,B,Cf,PLoss,thck1,Qact,
     &         IOUT,WELLID(iw),Kz)
c   if a vertical segment, use original function
          else
             cond2 = cel2wel2(LOSSTYPE,Txx1,Tyy1,dx1,dy1,
     &                 rw,Rskin,Kskin,B,Cf,PLoss,thck1,Qact,
     &                 IOUT,WELLID(iw),Skin)
          end if
      else
         cond2=0.D0 
      end if
      MNWNOD(31,INODE)=cond2
c sum cond for cell to get resultant CWC for node
      cond=cond1+cond2
c     Save conductance of each node
      MNWNOD(14,INODE) = cond
      if(ipr.eq.1) then 
       write(iout,'(A15,I3,1P6G12.5,9A)') WELLID(iw),nod,cond,
     & top1,bot1,top1,bot1,alpha,'         '
      end if
c
c process last node separately
c
      if(INODE.EQ.lastnode-1) then
c calculate CWC of first segment in node
      lw=MNWNOD(23,INODE+1)
      if (lw.gt.0.D0) then
          Txx = MNWNOD(16,INODE+1)
          Tyy = MNWNOD(17,INODE+1)
          Txx1 = Txx*0.5d0
          Tyy1 = Tyy*0.5d0
          rw  = MNWNOD(5,INODE+1)
          Rskin = MNWNOD(6,INODE+1)  
          Kskin = MNWNOD(7,INODE+1)
          B = MNWNOD(8,INODE+1)
          Cf = MNWNOD(9,INODE+1)
          PLoss = MNWNOD(10,INODE+1)
          Qact = MNWNOD(4,INODE+1)
c   compute conductance term for segment
          if(ivert1(INODE+1).eq.0) then
c  for 1-layer problems, Kz is not defined; assume=Kh     
C  also, FOR BCF, must assume Kh=Kz as only input is VCONT
                if(NLAY.EQ.1.or.INBCF.gt.0) then
                 T = (Txx*Tyy)**0.5D0
                 Kh = T/thck2
                 Kz=Kh
                else
C  FOR LPF==check LAYVKA and VKA
                  IF(LAYVKA(L2).EQ.0) THEN
                    Kz=VKA(C2,R2,L2)
                  ELSE
                    Kz=HK(C2,R2,L2)/VKA(C2,R2,L2)
                  END IF
                end if
             cond1 = cel2wel2SEG(lw,theta,omega,LOSSTYPE,
     &         Txx,Tyy,dx2,dy2,rw,Rskin,Kskin,B,Cf,PLoss,thck2,Qact,
     &         IOUT,WELLID(iw),Kz)
c   if a vertical segment, use original function
          else
              cond1 = cel2wel2(LOSSTYPE,Txx1,Tyy1,dx2,dy2,
     &                 rw,Rskin,Kskin,B,Cf,PLoss,thck2,Qact,
     &                 IOUT,WELLID(iw),Skin)
         end if
      else
         cond1=0.D0 
      end if
      MNWNOD(30,INODE+1)=cond1
c calculate CWC of second segment of node
c it is the same as the other segment in this node
      MNWNOD(24,INODE+1)=MNWNOD(23,INODE+1)
      cond2=cond1
      MNWNOD(31,INODE+1)=cond2
c sum cond for cell to get resultant CWC for node
      cond=cond1+cond2
c     Save conductance of each node
      MNWNOD(14,INODE+1) = cond      
      if(ipr.eq.1) then 
       write(iout,'(A15,I3,1P6G12.5,9A)') WELLID(iw),nod+1,cond,
     & top2,bot2,top2,bot2,alpha,'         '
      end if
      end if
c     end loop over "segments"
      end do
c     print segment info
      if(ipr.eq.1) then
       write(iout,*)
       write(iout,*) 'MNW2 Non-vertical Well:   Segment Information'
       write(iout,'(A)')  'Node   L   R   C   Segment   Length   
     &   DEG.TILT   MAP-ANGLE   CWC-segment'            
       do INODE=firstnode,lastnode
         L=MNWNOD(1,INODE)              
         R=MNWNOD(2,INODE)              
         C=MNWNOD(3,INODE)              
c segment 1
         lw=MNWNOD(23,INODE)
         if(inode.gt.1) then
           omega=MNWNOD(28,INODE-1)
           theta=MNWNOD(29,INODE-1)
         else
           omega=0.d0
           theta=0.d0
         end if
         cond1=MNWNOD(30,INODE)
         write(iout,'(4I4,I8,1pG16.6,1p3G12.5)')
     & INODE,L,R,C,1,lw,omega,theta,cond1
c segment 2
         lw=MNWNOD(24,INODE)
         if(inode.lt.lastnode) then
           omega=MNWNOD(28,INODE)
           theta=MNWNOD(29,INODE)
         else
           omega=MNWNOD(28,INODE-1)
           theta=MNWNOD(29,INODE-1)        
         end if
         cond2=MNWNOD(31,INODE)
         write(iout,'(4I4,I8,1pG16.6,1p3G12.5)') 
     & INODE,L,R,C,2,lw,omega,theta,cond2
c closed casings
         if(MNWNOD(25,INODE).GT.0) then
           write(iout,'(A,1pG16.6)') '   Closed casing length ',
     & MNWNOD(25,INODE)
         end if

       end do
      end if
c     
      DEALLOCATE(ivert1,ivert2)
C
      RETURN
C
      END
C
c
c_________________________________________________________________________________
c
      DOUBLE PRECISION function cel2wel2SEG(lw,theta,omega,LOSSTYPE,Txx,
     &  Tyy,dx,dy,rw,Rskin,Kskin,B,Cf,PLoss,thck,Q,IOUT,WELLNAME,Kz)
c
C
C     VERSION 20030327 KJH        -- Patched Hyd.K term in LPF solution
C     VERSION 20090405 GZH        -- MNW2
c
c----- MNW1 by K.J. Halford
c
c     ******************************************************************
c     Compute conductance term to define head loss from cell to wellbore
c      Methodology is described in full by Peaceman (1983)
c     ******************************************************************
C
C        SPECIFICATIONS:
C     ------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER LOSSTYPE,IOUT,i
      CHARACTER*20 WELLNAME
      DOUBLE PRECISION pi,verysmall,rw,Txx,Tyy,yx4,xy4,ro,dx,dy,Tpi2,A,
     & Ploss,B,Rskin,Kskin,C,Cf,Q,thck,T,Tskin,x1,x2,x3,x4,
     & roz,roy,rox,zx4,xz4,zy4,yz4,Az,Ay,Ax,Tpi2z,Tpi2y,Tpi2x,
     & theta,omega,kz,ky,kx,lw,bz,by,bx,clz,cly,clx,CLi,
     & numerator,denom1,denom2,lwz,lwy,lwx,
     & ct,st,cw,sw
c convert degree trig func modified from http://techpubs.sgi.com
      DOUBLE PRECISION dgr_to_rad
c      dsind = sin(dgr_to_rad * dgr_argument)
c      dcosd = cos(dgr_to_rad * dgr_argument)
C     ------------------------------------------------------------------
c define parameters
c
C
      pi = 3.1415926535897932D0
      dgr_to_rad = (pi/180.D0)
      verysmall = 1.D-25
c

      Kx=Txx/thck
      Ky=Tyy/thck
c    this makes conductance very small
      if( rw.lt.verysmall .or. Txx.lt.verysmall .or. Tyy.lt.verysmall ) 
     &  then
        cel2wel2SEG = ( Txx * Tyy )** 0.5D0
c       For the "NONE" option, multiply the Kh by 1000 to equivalate Hnew and hwell
      else if(LOSSTYPE.EQ.0) then
        cel2wel2SEG=1.0D3*((Kx*Ky)**0.5D0)   
      else 
c
c    define ro (effective radius) for each direction 
        yx4 = (Ky/Kx)**0.25D0
        xy4 = (Kx/Ky)**0.25D0
        roz = 0.28D0 *((yx4*dx)**2 +(xy4*dy)**2)**0.5D0 / (yx4+xy4)
        Tpi2z = 2.D0*pi * thck *(Kx*Ky)**0.5D0
c
        zx4 = (Kz/Kx)**0.25D0
        xz4 = (Kx/Kz)**0.25D0
        roy = 0.28D0 *((zx4*dx)**2 +(xz4*thck)**2)**0.5D0 / (zx4+xz4)
        Tpi2y = 2.D0*pi * dy * (Kx*Kz)**0.5D0
c
        yz4 = (Kz/Ky)**0.25D0
        zy4 = (Ky/Kz)**0.25D0
        rox = 0.28D0 *((yz4*dy)**2 +(zy4*thck)**2)**0.5D0 / (yz4+zy4)
        Tpi2x = 2.D0*pi * dx * (Kz*Ky)**0.5D0
c
c       if ro/rw is <1, 'A' term will be negative.  Warn user and cut off flow from this node
        if (rox/rw.lt.1.D0.or.roy/rw.lt.1.or.roz/rw.lt.1) then
          write(iout,*) 
     &      '     Ro_x/Rw =  ',Rox/Rw, 
     &      '     Ro_y/Rw =  ',Roy/Rw, 
     &      '     Ro_z/Rw =  ',Roz/Rw, 
     &      '***WARNING*** At least one value of Ro/Rw < 1,
     & CWC set = 0.0 for well '
          cel2wel2SEG = 0.D0
          GOTO 888
        end if
        Az = log(roz/rw) / Tpi2z
        Ay = log(roy/rw) / Tpi2y
        Ax = log(rox/rw) / Tpi2x
c
c       THEIM option (LOSSTYPE.EQ.1) only needs A, so no need to calculate  B or C
c
c       SKIN (LINEAR) option, calculate B, C=0
        if(LOSSTYPE.EQ.2) then
c         average T in aquifer assumed to be sqrt of Txx*Tyy
          if(Kskin.gt.0.D0.and.rw.gt.0.D0) then
c         this is from eqs 3 and 5 in orig MNW report
            lwz=thck
            Bz=(thck*(Kx*Ky)**0.5D0/(Kskin*lw)-1)*(DLOG(Rskin/rw))/Tpi2z
            lwy=dy
            By = (dy*(Kx*Kz)**0.5D0/(Kskin*lw)-1)*(DLOG(Rskin/rw))/Tpi2y
            lwx=dx
            Bx = (dx*(Ky*Kz)**0.5D0/(Kskin*lw)-1)*(DLOG(Rskin/rw))/Tpi2x
          else
            Bx = 0.D0
            By = 0.D0
            Bz = 0.D0
          end if
          C = 0.D0
c       NONLINEAR option, calculate B and C
       else if (LOSSTYPE.EQ.3) then
          B = B / Tpi2z 
          if(Cf.NE.0.0) then
            C = Cf * abs(Q)**(PLoss-1)
          else
            C = 0.D0
          end if
       else
          Bx = 0.D0
          By = 0.D0
          Bz = 0.D0
          C = 0.D0
       end if
c       these are per length
       CLz = Az + Bz + C 
       CLz = 1.000000D0 / CLz / thck
       CLy = Ay + By + C 
       CLy = 1.000000D0 / CLy / dy
       CLx = Ax + Bx + C 
       CLx = 1.000000D0 / CLx / dx
c calculate CWC for slanted well (from 2.45b in SUTRA doc)
       numerator=(CLz*CLy*CLx)
c      dsind = sin(dgr_to_rad * dgr_argument)
c      dcosd = cos(dgr_to_rad * dgr_argument)
       x1=dcos(dgr_to_rad * theta)
       x2=dsin(dgr_to_rad * theta)
       x3=dcos(dgr_to_rad * omega)
       x4=dsin(dgr_to_rad * omega)
       denom1=CLz*((CLy*(dcos(dgr_to_rad * theta)**2))
     &              +CLx*(dsin(dgr_to_rad * theta)**2))
     &              *dsin(dgr_to_rad * omega)**2
       denom2=CLx*Cly*(dcos(dgr_to_rad * omega)**2)
c 
       if((denom1+denom2).eq.0) then
         write(iout,*) '***ERROR*** MNW2 slanted well error'
         STOP 'MNW2 -- slanted well'
       end if
	 CLi=numerator/(denom1+denom2)
       cel2wel2SEG=lw*(numerator/(denom1+denom2))
      endif
c
 888  continue
      end
c
C
C
C MNW2CAPACITY
C
C Compute Qact restrained by pumping capacity
C
C     ******************************************************************
C
      SUBROUTINE MNW2CAPACITY(qactCap,WELLID,mnwmax,Hlift,hwell,iw,
     & PUMPCAP,CapTable,iout)
C
C     ******************************************************************
      IMPLICIT NONE
      CHARACTER*20 WELLID
      INTEGER mnwmax,iw,PUMPCAP,idone,index,iout,ifirstL,isecondL
      DOUBLE PRECISION qactCap,LIFTact,Hlift,hwell,CapTable,m,b,
     & L1,L2,Q1,Q2
      DIMENSION WELLID(mnwmax),CapTable(mnwmax,27,2)    
C 
C
C     Compute lift
      LIFTact=Hlift-hwell
C     
      qactCap=0.d0
c     progress flag: idone=1 mean have interp points; idone=2 means have Q value
      idone=0
c     if actual lift is greater than first value in table, use Q for first value
      if(LIFTact.gt.CapTable(iw,1,1)) then
        qactCap=CapTable(iw,1,2)
        idone=2
      end if
c     if actual lift is less than final value in table, use Q for final value
      if(LIFTact.lt.CapTable(iw,PUMPCAP+2,1)) then
        qactCap=CapTable(iw,PUMPCAP+2,2)
        idone=2
      end if
C     Loop over CapTable to check for table entry matches or to find encompassing Lift values
      if(idone.eq.0) then
        do index=1,PUMPCAP+2
c     if actual lift equals one of the table entries, set Q and done
          if(LIFTact.eq.CapTable(iw,index,1)) then
            qactCap=CapTable(iw,index,2)
            idone=2
          end if
c     if LIFTact is an intermediate value, find first entry it is less than; this
c     will define which two value to use in interpolation
          if(idone.eq.0) then
           if(index.lt.(PUMPCAP+2)) then
            if(LIFTact.gt.CapTable(iw,index+1,1)) then
              ifirstL=index
              isecondL=index+1
              idone=1
            end if
           else
c     if table is constructed properly, this should never be executed (index=PUMPCAP+2)
            write(iout,*) '***ERROR*** MNW2 Capacity table read error'
            STOP 'MNW2 ERROR - CapTable'
           end if  
          end if
        end do
      end if
c     error check; idone should be set by now
      if(idone.eq.0) then
        write(iout,*) '***ERROR*** MNW2 Capacity table read error'
        STOP 'MNW2 ERROR - CapTable'
      end if      
C     Interpolate Q value from table
      if(idone.eq.1) then
c     define points
        L1=CapTable(iw,ifirstL,1)
        L2=CapTable(iw,isecondL,1)
        Q1=CapTable(iw,ifirstL,2)
        Q2=CapTable(iw,isecondL,2)
c     calculate slope and intercept of line between the two points
        m=(Q2-Q1)/(L2-L1)
        b=Q1-(m*L1)
c     interpolate by finding Q on the line segment for actual lift
        qactCap=m*LIFTact+b
      end if
c     convert discharge to MODFLOW sign convention
      qactCap=qactCap*(-1.d0)
c      
      RETURN
      END
C
C
C
C GWF1MNW2bo
C
C Separate output for GWT version so conc can be printed
C
C     ******************************************************************
C
      SUBROUTINE GWF1MNW2bo(IMNWCB,ICBCFL,NAUX,KSTP,KPER,
     & NCOL,NROW,NLAY,nmnw2,mnw2,iout,DELT,PERTIM,TOTIM,IBOUND,
     & MNWMAX,NODTOT,msum,HDRY,VBNM,VBVL,BUFF,MNWNOD,hnew,WELLID,
     & MNWPRNT,IGWTON,MNWOBS,gwtunit,NMNWVL)
C
C     ******************************************************************
C
      IMPLICIT NONE
      INTEGER ibd,IMNWCB,ICBCFL,NAUX,NCOL,NROW,NLAY,KSTP,KPER,nmnw2,
     & iout,ibound,MNWMAX,imult,iw,firstnode,lastnode,inode,il,ir,ic,
     & NODTOT,ioch,ipole,ioc,nmnwvl,msum,nd,iweldry,MNWPRNT,IGWTON,
     & MNWOBS,gwtunit
C
      REAL HDRY, VBVL, BUFF, PERTIM, TOTIM, DELT
      DOUBLE PRECISION ratin,ratout,mnw2,MNWNOD,hnew,DryTest,
     & q,qdes,hlim,hwell,s,sNL,sL,qnet,qin,qout,verysmall,hcell,small,
     & cnode
cgzh debug auxtext will be replaced with common block mnwaux
      CHARACTER*20 WELLID,blank
      CHARACTER*16 text,vbnm(msum),AUXTXT(5)
      DIMENSION ibound(ncol,nrow,nlay),mnw2(32,MNWMAX),
     & buff(NCOL,NROW,NLAY),mnwnod(33,NODTOT),hnew(NCOL,NROW,NLAY),
     & vbvl(4,msum)
      dimension WELLID(mnwmax)     
c -----print the header for individual rates if requested(IMNWCB<0).
      if( IMNWCB.lt.0 .and. icbcfl.ne.0 ) then
       if(gwtunit.le.0) then
        write(iout,'(/,1x,a,9h PERIOD =,i5,8h  STEP =,i5)')
     +           'Nodal information for all MNW2 wells,',  kper,kstp
           write(iout,'(101A)') 'WELLID                NODE   Lay 
     &  Row   Col        Totim        Q-node         hwell         hcell
     &   Seepage elev.'
       else
        write(iout,'(/,1x,a,9h PERIOD =,i5,8h  STEP =,i5)')
     +           'Nodal information for all MNW2 wells,',  kper,kstp
           write(iout,'(117A)') 'WELLID                NODE   Lay 
     &  Row   Col        Totim        Q-node         hwell         hcell
     &   Seepage elev.  Concentration'
       end if
      endif
c
c2------if there are no wells do not accumulate flow
      if(nmnw2.gt.0) then
c  test for dry wells
        imult = 0
c  Loop over all wells
        do iw=1,MNWMAX
          firstnode=MNW2(4,iw)
          lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
c   Loop over nodes in well
          do INODE=firstnode,lastnode
            il=MNWNOD(1,INODE)              
            ir=MNWNOD(2,INODE)              
            ic=MNWNOD(3,INODE)              
            DryTest = Hnew(ic,ir,il) - Hdry
            if(ABS(DryTest).lt.verysmall) then
              MNWNOD(4,INODE)= 0.0D0
              MNW2(17,iw)=Hdry
              iweldry=1
              if(MNWPRNT.gt.1) then
                write(iout,*) 'MNW2 node in dry cell, Q set to 0.0'
                write(iout,*) 'Well: ',WELLID(iw),' Node: ',INODE
              end if
            else
              iweldry=0
            endif
            q = MNWNOD(4,INODE)
c
c    Report all wells with production less than the desired rate......
            if(ibound(ic,ir,il).ne.0 .or.ABS(DryTest).lt.verysmall) then
              if(MNW2(6,iw).GT.0) then
                hlim = mnw2(7,iw)
              else
                hlim = 0.0D0
              end if
c
              hwell = mnw2(17,iw)
c
              ioch = 0
              if( IMNWCB.lt.0 .and. icbcfl.ne.0 ) ioch = 1
c -----print the individual rates if requested(IMNWCB<0).
              if( ioch.eq.1 ) then
              q=MNWNOD(4,INODE)
              hcell=hnew(ic,ir,il)
              nd=INODE-firstnode+1
              blank='              '
              if(gwtunit.gt.0) cnode=MNWNOD(32,INODE)
c   If no seepage face in cell (true for single-node wells), don't print seepage elev.
              if(MNWNOD(15,INODE).EQ.hwell.or.
     &           MNWNOD(15,INODE).eq.Hdry.or.
     &           firstnode.eq.lastnode) then
               if(gwtunit.le.0) then
                write(iout,'(A20,4i6,1x,1P4e14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell
               else
                write(iout,'(A20,4i6,1x,1P4e14.6,A14,1Pe14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,
     &              blank,cnode
               endif   
              else
c   If seepage face in cell, MNWNOD(15) will hold the bottom elev of the cell,
c   which is used with hcell to get the gradient used to calculate the Q for the
c   seepage face.  
               if(gwtunit.le.0) then
                write(iout,'(A20,4i6,1x,1P5e14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,
     &              MNWNOD(15,INODE)
               else
                write(iout,'(A20,4i6,1x,1P6e14.6)')
     &              WELLID(iw),nd,il,ir,ic,totim,q,hwell,hcell,
     &              MNWNOD(15,INODE),cnode
               endif   
              end if
              endif
          endif
         enddo
       enddo
c
c   Sum components of  multi-node wells
c
c -----print the header for multi-node rates if requested(IMNWCB<0).
        if(  ioch.eq.1  ) then
          write(iout,'(/,5x,a )') 'Summary information for MNW2 wells'
           write(iout,'(200A)') 'WELLID                    Totim    
     &        Qin           Qout           Qnet          hwell'
        endif
c
c  Loop over all wells
        do iw=1,MNWMAX
c        if NNODES>1, that is, this is a multi-node well
            firstnode=MNW2(4,iw)
            lastnode=MNW2(4,iw)+ABS(MNW2(2,iw))-1
c   Loop over nodes in well
            qnet = 0.000D0
            qin  = 0.000D0
            qout = 0.000D0
            do INODE=firstnode,lastnode
              il=MNWNOD(1,INODE)              
              ir=MNWNOD(2,INODE)              
              ic=MNWNOD(3,INODE)              
              if( ibound(ic,ir,il).eq.0 ) MNWNOD(4,INODE) = 0.0D0
              if( MNWNOD(4,INODE).le.0.0D0 ) then
                qin = qin  + MNWNOD(4,INODE)
              else
                qout = qout  + MNWNOD(4,INODE)
              endif
              qnet  = qnet  + MNWNOD(4,INODE)
            enddo
            mnw2(18,iw) = qnet
c
c  if Q is constrained, print message
            if(mnwprnt.gt.0) then
            qdes = mnw2(5,iw)
             if(MNW2(6,iw).ne.0.or.MNW2(22,iw).ne.0) then
              if(abs(qnet-qdes).gt.small) then
              write(iout,*)
               if(abs((mnw2(29,iw)-qnet)).lt.small) then
              write(iout,*) WELLID(iw),' Qdes has been updated to',
     & qnet,' because of Well Capacity restraint'
cgzh debug output
c       write(*,*) 'mnw2(29,iw),qnet',mnw2(29,iw),qnet
               else
cgzh debug output
c       write(*,*) 'mnw2(29,iw),qnet',mnw2(29,iw),qnet
              write(iout,*) WELLID(iw),' Qdes has been updated to',
     & qnet,' because of Hlim constraint'
               end if
              end if
             end if
            end if
cgzh debug output
c           write(*,*) 'mnw2(18,iw) in BD=',mnw2(18,iw)
c -----print the summed rates if requested(IMNWCB<0).
            hwell = mnw2(17,iw)
            if(  ioch.eq.1  ) then
            if(qnet.lt.(small*qin)) qnet=0.d0
            write(iout,'(A20,5(1x,1Pg14.6))')
     &              WELLID(iw),totim,qin,qout,qnet,hwell
            endif
        enddo
            if(  ioch.eq.1  ) write(iout,*)
      end if
c
c  ----- END  MULTI-NODE reporting section -------------
      return
      end