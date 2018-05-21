C
C      ********************************************************
C      *                                                      *
C      *                   **** POINT3 ****                   *
C      *                                                      *
C      *   THREE-DIMENSIONAL GROUND-WATER SOLUTE TRANSPORT    *
C      *                                                      *
C      *    MODEL FOR AN AQUIFER OF INFINITE EXTENT WITH A    *
C      *                                                      *
C      *   CONTINUOUS POINT SOURCE AT X=XC, Y=YC, AND Z=ZC    *
C      *                                                      *
C      *        GROUND-WATER FLOW IN X-DIRECTION ONLY         *
C      *                                                      *
C      *            VERSION CURRENT AS OF 04/01/90            *
C      *              CODE CLEANUP 04/02/96 - RSR             *
C      *                                                      *
C      ********************************************************
C
C
C      Although this program has been used by the U.S. Geological
C      Survey, no warranty, expressed or implied, is made by the USGS
C      as to the accuracy and functioning of the program and related
C      program material, nor shall the fact of distribution constitute
C      any such warranty, and no responsibility is assumed by the USGS
C      in connection therewith.
C
      INCLUDE 'dimens.inc'
      INCLUDE 'pltdat.inc'
      INCLUDE 'units.inc'
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
      DOUBLE PRECISION cxy(MAXX,MAXY), x(MAXX), y(MAXY), z(MAXZ), 
     &                 t(MAXT), c0, cn, dk, dx, dy, dz, por, qm, vx, xc, 
     &                 xx, yc, yy, zc, zz
      INTEGER it, ix, iy, iz, j, np, np1, np2, npage, nt, nx, ny, nz,
     &        ifx
C
C---EXTERNALS
      EXTERNAL OFILE, TITLE, CNRML3, PLTNC3, RDXYZT, RDPLOT, WCONH3
C
C     PROGRAM VARIABLES
C
C     C0       SOLUTE CONCENTRATION IN INJECTED FLUID [M/L**3]
C     DX       LONGITUDINAL DISPERSION COEFFICIENT [L**2/T]
C     DY       TRANSVERSE (Y-DIRECTION) DISPERSION COEFFICIENT [L**2/T]
C     DZ       TRANSVERSE (Z-DIRECTION) DISPERSION COEFFICIENT [L**2/T]
C     VX       GROUND-WATER VELOCITY IN X-DIRECTION [L/T]
C     DK       FIRST-ORDER SOLUTE DECAY CONSTANT [1/T]
C     X        X-POSITION AT WHICH CONCENTRATION IS EVALUATED [L]
C     Y        Y-POSITION AT WHICH CONCENTRATION IS EVALUATED [L]
C     Z        Z-POSITION AT WHICH CONCENTRATION IS EVALUATED [L]
C     T        TIME AT WHICH CONCENTRATION IS EVALUATED [T]
C     CN       NORMALIZED CONCENTRATION C/C0 [DIMENSIONLESS]
C     CXY      SOLUTE CONCENTRATION C(X,Y,Z,T) [M/L**3]
C     XC       X-COORDINATE OF CONTINUOUS POINT SOURCE [L]
C     YC       Y-COORDINATE OF CONTINUOUS POINT SOURCE [L]
C     ZC       Z-COORDINATE OF CONTINUOUS POINT SOURCE [L]
C     QM       FLUID INJECTION RATE [L**3/T]
C     POR      AQUIFER POROSITY [DIMENSIONLESS]
C
C     NX       NUMBER OF X-POSITIONS AT WHICH SOLUTION IS EVALUATED
C     NY       NUMBER OF Y-POSITIONS AT WHICH SOLUTION IS EVALUATED
C     NZ       NUMBER OF Z-POSITIONS AT WHICH SOLUTION IS EVALUATED
C     NT       NUMBER OF TIME VALUES AT WHICH SOLUTION IS EVALUATED
C
C     IPLT     PLOT CONTROL. IF IPLT>0, CONTOUR MAPS ARE PLOTTED
C     XSCLP    SCALING FACTOR TO CONVERT X TO PLOTTER INCHES
C     YSCLP    SCALING FACTOR TO CONVERT Y TO PLOTTER INCHES
C     DELTA    CONTOUR INCREMENT FOR PLOT. (VALUE BETWEEN 0 AND 1.0)
C
      Versn = 
     & '$Id: point3.f,v 1.1 1996/04/03 19:47:51 rsregan Exp rsregan $'
      Versn = 
     &    '@(#)POINT3 - Analytical solution for 3-D GW solute transport'
      Versn = 
     &      '@(#)POINT3 - in an infinite system with uniform flow and a'
      Versn = '@(#)POINT3 - continuous point solute source'
      Versn = '@(#)POINT3 - USGS TWRI, book 3, chap. B7, E.J. Wexler'
      Versn = '@(#)POINT3 - Contact: h2osoft@usgs.gov'
      Versn = '@(#)POINT3 - Version: 1.1x 1996/04/03'
C
C---DEFINE INPUT/OUTPUT FILES AND PRINT TITLE PAGE
      CALL OFILE
      CALL TITLE
      WRITE (Io,9020)
C
C---READ IN MODEL PARAMETERS
      READ (In,9005) nx, ny, nz, nt, Iplt
      WRITE (Io,9025) nx, ny, nz, nt
      READ (In,9010) Cunits, Vunits, Dunits, Kunits, Lunits, Qunits, 
     &               Tunits
      READ (In,9015) c0, vx, dx, dy, dz, dk
      WRITE (Io,9030) c0, Cunits, vx, Vunits, dx, Dunits, dy, Dunits, 
     &                dz, Dunits, dk, Kunits
      READ (In,9015) xc, yc, zc, qm, por
      WRITE (Io,9035) xc, Lunits, yc, Lunits, zc, Lunits, qm, Qunits, 
     &                por
      CALL RDXYZT(nx,ny,nz,nt,x,y,z,t)
      IF (Iplt.GT.0) CALL RDPLOT(Cunits,Xsclp,Ysclp,Delta)
C
C---BEGIN TIME LOOP
      DO 60 it = 1, nt
C
C---BEGIN Z LOOP
        DO 50 iz = 1, nz
          zz = z(iz) - zc
C
C---BEGIN X LOOP
          DO 20 ix = 1, nx
            xx = x(ix) - xc
C
C---CALCULATE NORMALIZED CONCENTRATION FOR ALL Y AT X=X(IX) AND Z=Z(IZ)
            DO 10 iy = 1, ny
              yy = y(iy) - yc
              CALL CNRML3(qm,por,dk,t(it),xx,yy,zz,dx,dy,dz,vx,cn)
              cxy(ix,iy) = c0*cn
   10       CONTINUE
   20     CONTINUE
C
C---PRINT OUT TABLES OF CONCENTRATION VALUES
          npage = 1 + (ny-1)/9
          WRITE (Io,9040) t(it), Tunits, z(iz), Lunits, Lunits
          ifx = 0
          DO 40 np = 1, npage
            np1 = (np-1)*9
            np2 = 9
            IF ((np1+np2).GT.ny) np2 = ny - np1
            CALL WCONH3(ifx,np1,np2,t(it),z(iz),y)
            ifx = 45
            DO 30 ix = 1, nx
              WRITE (Io,9045) x(ix), (cxy(ix,np1+j),j=1,np2)
              CALL WCONH3(ix,np1,np2,t(it),z(iz),y)
   30       CONTINUE
   40     CONTINUE
C
C---PLOT CONTOUR MAPS OF NORMALIZED CONCENTRATION FOR EACH TIME VALUE.
          IF (Iplt.GT.0) CALL PLTNC3(nx,ny,nz,iz,x,y,t(it),z(iz),cxy,c0)
   50   CONTINUE
   60 CONTINUE
      CLOSE (In)
      CLOSE (Io)
      STOP
C
C---FORMAT STATEMENTS
 9005 FORMAT (20I4)
 9010 FORMAT (8A10)
 9015 FORMAT (8F10.0)
 9020 FORMAT (/////,30X,'ANALYTICAL SOLUTION TO THE THREE-DIMENSIONAL',
     &        /,29X,'ADVECTIVE-DISPERSIVE SOLUTE TRANSPORT EQUATION',/,
     &        35X,'FOR AN AQUIFER OF INFINITE EXTENT',/,31X,
     &        'WITH A CONTINUOUS POINT SOURCE AT XC,YC,ZC',////,41X,
     &        'INPUT DATA',/,41X,10('-'))
 9025 FORMAT (/,26X,'NUMBER OF X-COORDINATES (NX) =',I5,/,26X,
     &        'NUMBER OF Y-COORDINATES (NY) =',I5,/,26X,
     &        'NUMBER OF Z-COORDINATES (NZ) =',I5,/,26X,
     &        'NUMBER OF TIME VALUES (NT) =',I5)
 9030 FORMAT (/,26X,'SOLUTE CONCENTRATION IN INJECTED FLUID (C0) =',
     &        1P1E13.6,1X,A10,/,26X,
     &        'GROUND-WATER VELOCITY IN X-DIRECTION (VX) =',1P1E13.6,1X,
     &        A10,/,26X,'DISPERSION IN THE X-DIRECTION (DX) =',1P1E13.6,
     &        1X,A10,/,26X,'DISPERSION IN THE Y-DIRECTION (DY) =',
     &        1P1E13.6,1X,A10,/,26X,
     &        'DISPERSION IN THE Z-DIRECTION (DZ) =',1P1E13.6,1X,A10,/,
     &        26X,'FIRST-ORDER SOLUTE DECAY RATE (DK) =',1P1E13.6,1X,
     &        A10)
 9035 FORMAT (/,26X,'AQUIFER IS OF INFINITE EXTENT',//,26X,
     &        'CONTINUOUS POINT SOURCE IS AT X =',1P1E13.6,1X,A10,/,56X,
     &        'Y =',1P1E13.6,1X,A10,/,56X,'Z =',1P1E13.6,1X,A10,/,26X,
     &        'FLUID INJECTION RATE (QM) =',1P1E13.6,1X,A10,/,26X,
     &        'AQUIFER POROSITY (POR) =',1P1E13.6)
 9040 FORMAT ('1',//,16X,'SOLUTE CONCENTRATION AT TIME =',F12.4,1X,A10,
     &        /,36X,'AND AT Z =',F12.4,1X,A10,//,26X,
     &        'Y-COORDINATE, IN ',A10)
 9045 FORMAT (6X,F12.4,'  !',9F12.6)
      END
C
C
C
      SUBROUTINE CNRML3(QM,POR,DK,T,X,Y,Z,DX,DY,DZ,VX,CN)
      DOUBLE PRECISION QM, POR, DK, T, X, Y, Z, DX, DY, DZ, VX, CN
C
C       THIS ROUTINE CALCULATES SOLUTE CONCENTRATION AT X,Y,Z BASED ON
C       THE ANALYTIC SOLUTION TO THE THREE-DIMENSIONAL ADVECTIVE-
C       DISPERSIVE SOLUTE TRANSPORT EQUATION FOR AN AQUIFER OF
C       INFINITE EXTENT WITH A CONTINUOUS POINT SOURCE LOCATED AT
C       X=XC, Y=YC, AND Z=ZC.  A CLOSED FORM SOLUTION WAS OBTAINED.
C
      DOUBLE PRECISION beta, gamma, pi, rtdxt, x1, x2, y1, y2, z1, z2, 
     &                 z3, z4
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
C---EXTERNALS
      EXTERNAL EXERFC
C
C---INTRINSICS
      DOUBLE PRECISION DSQRT
      INTRINSIC DSQRT
C
C---FOR T=0, ALL CONCENTRATIONS EQUAL 0.0
      IF (T.LE.0.0D0) THEN
        CN = 0.0D0
C
C---CHECK FOR X=Y=Z=0
      ELSEIF (X.EQ.0.0D0 .AND. Y.EQ.0.0D0 .AND. Z.EQ.0.0D0) THEN
        WRITE (Io,9005)
      ELSE
C
        beta = DSQRT(VX*VX+4.0D0*DX*DK)
        gamma = DSQRT(X*X+Y*Y*DX/DY+Z*Z*DX/DZ)
        rtdxt = 2.0D0*DSQRT(DX*T)
C
C---TERM 1
        x1 = (VX*X-gamma*beta)/(2.0D0*DX)
        y1 = (gamma-beta*T)/rtdxt
        CALL EXERFC(x1,y1,z1)
C
C---TERM 2
        x2 = (VX*X+gamma*beta)/(2.0D0*DX)
        y2 = (gamma+beta*T)/rtdxt
        CALL EXERFC(x2,y2,z2)
C
C---TERM 3
        z3 = z1 + z2
        z4 = DSQRT(DY*DZ)
        pi = 3.14159265358979D0
        CN = QM*z3/(8.0D0*POR*pi*gamma*z4)
      ENDIF
C
      RETURN
C
C---FORMAT STATEMENTS
 9005 FORMAT (/,6X,
     &'**** WARNING ****  A SOLUTION CAN NOT BE COMPUTED FOR X=XC,Y=YC,Z
     &=ZC',/)
      END
