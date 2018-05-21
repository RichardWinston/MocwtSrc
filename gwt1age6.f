C     Last change:  Dec. 8, 2014
C     Weighted particle subroutine added.
C  AGE6DF  READ AGING RATE
C*******************************************
C
      SUBROUTINE AGE6DF(INAGE,AGER8,IOUT,IOUTS)
C
C*******************************************
      READ(INAGE,*) AGER8
      WRITE(IOUTS,111) AGER8
      RETURN
C
111     format(//' **************'/
     1'   AGE OPTION   '/' **************'//
     2'  CONCENTRATION OUTPUT IS GROUND-WATER AGE;',
     3' SEE:'/
     4'  Goode, D.J., 1999, Age, double porosity, and simple reaction'/
     5'       modifications to the MOC3D ground-water transport '/
     6'       model: U.S. Geological Survey WRI 99-4041, 34 p.'/
     8'  SOURCE CONCENTRATIONS WILL BE USED AS SOURCE AGES AND',
     9' SHOULD BE CHECKED'/
     A'   GENERALLY SOURCE AGES SHOULD BE SET TO ZERO'/
     8'  INITIAL CONCENTRATION WILL BE USED AS INITIAL AGE AND',
     9' SHOULD BE CHECKED'/
     B'  RETARDATION FACTOR (RF) WILL BE USED AND SHOULD',
     C' BE CHECKED'/
     D'   GENERALLY, RF=1.0 FOR AGE SIMULATION'/
     E' AGING RATE =',1PG11.4/
     F'   (RATIO OF OUTPUT AGE UNITS TO MODEL TIME UNITS, ',
     G'USUALLY 1.0)'/)
C
      END
C
C  AGE6AD
c***********************************************************************
c
      subroutine AGE6AP(pc,pr,pl,pconc,conc,ibound,AMASIN,
     1  THCK,POR,DAGE,ncol,nrow,nlay,np,nscol,nsrow,nslay)
c
c***********************************************************************
c
c     increase age (concentration) of particles and nodes by length of
c       transport time step; generally ager8 = timv
c
      DOUBLE PRECISION AMASIN,DAMAS
      dimension pc(np),pr(np),pl(np)
      dimension pconc(np),conc(nscol,nsrow,nslay)
      dimension ibound(ncol,nrow,nlay)
      DIMENSION THCK(NsCOL,NsROW,NsLAY),POR(NsCOL,NsROW,NsLAY)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c
      do 10 ip=1,np
        if(pc(ip).gt.0.0) pconc(ip)=pconc(ip)+DAGE
10      continue
c
C  ZERO MASS TERM
      DAMAS=0.0
      do 20 ks=1,nslay
      k=ks+islay1-1
      do 20 is=1,nsrow
      i=is+isrow1-1
      do 20 js=1,nscol
      j=js+iscol1-1
        if(ibound(j,i,k).ne.0) THEN
          conc(js,is,ks)=conc(js,is,ks)+DAGE
          DAMAS=DAMAS+POR(Js,Is,Ks)*THCK(Js,Is,Ks)*DAGE
        END IF
20      continue
C
      AMASIN=AMASIN+DAMAS*CDEL*RDEL
c
      return
      end

c***********************************************************************
C
c***********************************************************************
c
      subroutine AGE6APWT(pc,pr,pl,pconc,AMASIN,
     1  ptwt,DAGE,np,conc,NSCOL,NSROW,NSLAY,IBOUND,ncol,nrow,nlay)
c
c***********************************************************************
c
c     increase age (concentration) of particles and nodes by length of
c       transport time step; generally ager8 = timv
c
      DOUBLE PRECISION AMASIN,ptwt,DAMAS
      dimension pc(np),pr(np),pl(np)
      dimension pconc(np),ptwt(np)
      dimension conc(nscol,nsrow,nslay)
      dimension IBOUND(ncol,nrow,nlay)
C
      COMMON /GWT/ CDEL,RDEL,CNOFLO,CELDIS,FZERO,NZCRIT
      COMMON /SUBGRD/
     *  ISCOL1,ISCOL2,ISROW1,ISROW2,ISLAY1,ISLAY2,ISUBGD
c
C  ZERO MASS TERM
      DAMAS=0.0
      do 10 ip=1,np
        if(pc(ip).gt.0.0) then
	    pconc(ip)=pconc(ip)+DAGE
          DAMAS=DAMAS+ptwt(ip)*DAGE
        end if
10      continue
C
      AMASIN=AMASIN+DAMAS
C
      do 20 ks=1,nslay
      k=ks+islay1-1
      do 20 is=1,nsrow
      i=is+isrow1-1
      do 20 js=1,nscol
      j=js+iscol1-1
        if(ibound(j,i,k).ne.0) THEN
          conc(js,is,ks)=conc(js,is,ks)+DAGE
        END IF
20      continue
c
      return
      end
c***********************************************************************
