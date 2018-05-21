C
C      ********************************************************
C      *                                                      *
C      *                   **** FINITE ****                   *
C      *                                                      *
C      * ONE-DIMENSIONAL GROUND-WATER SOLUTE-TRANSPORT MODEL  *
C      *                                                      *
C      *   FOR A FINITE SYSTEM WITH A FIRST- OR THIRD-TYPE    *
C      *                                                      *
C      *              BOUNDARY CONDITION AT X=0               *
C      *                                                      *
C      *            VERSION CURRENT AS OF 04/01/90            *
C      *              CODE CLEANUP 03/27/96 - RSR             *
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
      INCLUDE 'units.inc'
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
      DOUBLE PRECISION cxt(MAXX,MAXT), x(MAXX), t(MAXT), root(MAXRT), a, 
     &                 c, c0, cn, dk, dx, vx, xl
      REAL xsclp
      INTEGER iplt, it, ix, j, nbc, np, np1, np2, npage, nroot, nt, nx, 
     &        ifx
      CHARACTER*1 ierr(MAXX,MAXT)
C
C---EXTERNALS
      EXTERNAL OFILE, TITLE, CNORM1, CNORM3, ROOT1, ROOT3, READXT, 
     &         WCONHD
C
C     PROGRAM VARIABLES
C
C     C0       SOLUTE CONCENTRATION AT THE INFLOW BOUNDARY [M/L**3]
C     DX       LONGITUDINAL DISPERSION COEFFICIENT [L**2/T]
C     VX       GROUND-WATER VELOCITY IN X-DIRECTION [L/T]
C     DK       FIRST-ORDER SOLUTE DECAY CONSTANT [1/T]
C     X        X-POSITION AT WHICH CONCENTRATION IS EVALUATED [L]
C     T        TIME AT WHICH CONCENTRATION IS EVALUATED [T]
C     CN       NORMALIZED CONCENTRATION C/C0 [DIMENSIONLESS]
C     CXT      SOLUTE CONCENTRATION C(X,T) [M/L**3]
C     XL       LENGTH OF THE FLOW SYSTEM [L]
C     ROOT(N)  ROOTS OF EQ. USED IN INFINITE SERIES SUMMATION
C
C     NBC      SOURCE BOUNDARY CONDITION TYPE (1 OR 3)
C     NX       NUMBER OF X-POSITIONS AT WHICH SOLUTION IS EVALUATED
C     NT       NUMBER OF TIME VALUES AT WHICH SOLUTION IS EVALUATED
C     NROOT    NUMBER OF ROOTS USED IN INFINITE SERIES SUMMATION
C     IPLT     PLOT CONTROL. IF IPLT>0, CONCENTRATION PROFILES ARE PLOTTED
C
      Versn = 
     & '$Id: finite.f,v 1.1 1996/04/03 20:05:22 rsregan Exp rsregan $'
      Versn = 
     &    '@(#)FINITE - Analytical solution for 1-D GW solute transport'
      Versn = '@(#)FINITE - in a finite system with uniform flow'
      Versn = '@(#)FINITE - USGS TWRI, book 3, chap. B7, E.J. Wexler'
      Versn = '@(#)FINITE - Contact: h2osoft@usgs.gov'
      Versn = '@(#)FINITE - Version: 1.1x 1996/04/03'
C
C---DEFINE INPUT/OUTPUT FILES AND PRINT TITLE PAGE
      CALL OFILE
      CALL TITLE
      WRITE (Io,9020)
C
C---READ IN MODEL PARAMETERS
      READ (In,9005) nbc, nx, nt, nroot, iplt
      IF (nbc.EQ.1) WRITE (Io,9025) 'FIRST'
      IF (nbc.EQ.3) WRITE (Io,9025) 'THIRD'
      WRITE (Io,9030) nx, nt, nroot
      READ (In,9010) Cunits, Vunits, Dunits, Kunits, Lunits, Tunits
      READ (In,9015) c0, vx, dx, dk, xl, xsclp
      WRITE (Io,9035) c0, Cunits, vx, Vunits, dx, Dunits, dk, Kunits, 
     &                xl, Lunits, xsclp
      CALL READXT(nx,nt,x,t)
C
C---GET EIGENVALUES (BETA) USED IN SERIES SUMMATION BY SOLVING FOR
C---THE POSITIVE ROOTS OF: BETA*COTAN(BETA)+VX*XL/(2*DX)=0.0
C---FOR A FIRST-TYPE SOURCE BOUNDARY CONDITION,
C---OR: BETA*COTAN(BETA)-BETA**2*DX/(VX*XL)+VX*XL/(4*DX)=0.0
C---FOR A THIRD-TYPE SOURCE BOUNDARY CONDITION.
C
      IF (nbc.EQ.1) THEN
        c = vx*xl/(2.0D0*dx)
        CALL ROOT1(c,root,nroot)
      ELSE
        a = 0.250D0*vx*xl/dx
        c = dx/(xl*vx)
        CALL ROOT3(a,c,root,nroot)
      ENDIF
C
C---BEGIN TIME LOOP
      DO 20 it = 1, nt
C
C---BEGIN X-COORDINATE LOOP
        DO 10 ix = 1, nx
C
C---CALL ROUTINE TO CALCULATE NORMALIZED CONCENTRATION
C---BASED ON TYPE OF BOUNDARY CONDITION SPECIFIED
          IF (nbc.EQ.1) CALL CNORM1(xl,t(it),x(ix),dx,vx,dk,root,cn,
     &                              nroot,ierr(ix,it))
          IF (nbc.EQ.3) CALL CNORM3(xl,t(it),x(ix),dx,vx,dk,root,cn,
     &                              nroot,ierr(ix,it))
          cxt(ix,it) = cn*c0
   10   CONTINUE
C
C---PLOT NORMALIZED CONCENTRATION PROFILE FOR EACH TIME VALUE.
        IF (iplt.GT.0) CALL PLOTNC(nx,nt,it,x,t(it),cxt,c0,xsclp)
   20 CONTINUE
C
C---PRINT OUT TABLES OF CONCENTRATION VALUES
      npage = 1 + (nt-1)/9
      WRITE (Io,9040) Tunits
      ifx = 0
      DO 40 np = 1, npage
        np1 = (np-1)*9
        np2 = 9
        IF ((np1+np2).GT.nt) np2 = nt - np1
        CALL WCONHD(ifx,np1,np2,t)
        ifx = 45
        DO 30 ix = 1, nx
          WRITE (Io,9045) x(ix), (cxt(ix,np1+j),ierr(ix,np1+j),j=1,np2)
          CALL WCONHD(ix,np1,np2,t)
   30   CONTINUE
   40 CONTINUE
C
      CLOSE (In)
      CLOSE (Io)
      STOP
C
C---FORMAT STATEMENTS
 9005 FORMAT (20I4)
 9010 FORMAT (8A10)
 9015 FORMAT (8F10.0)
 9020 FORMAT (/////,31X,'ANALYTICAL SOLUTION TO THE ONE-DIMENSIONAL',/,
     &        29X,'ADVECTIVE-DISPERSIVE SOLUTE-TRANSPORT EQUATION',/,
     &        37X,'FOR A SYSTEM OF FINITE LENGTH',////,41X,
     &        'INPUT DATA',/,41X,10('-'))
 9025 FORMAT (/,26X,A,'-TYPE BOUNDARY CONDITION AT X = 0.0')
 9030 FORMAT (/,26X,'NUMBER OF X-COORDINATES (NX) =',I5,/,26X,
     &        'NUMBER OF TIME VALUES (NT) =',I5,/,26X,
     &     'NUMBER OF ROOTS USED IN INFINITE SERIES SUMMATION (NROOT) ='
     &     ,I5)
 9035 FORMAT (/,26X,'SOLUTE CONCENTRATION ON MODEL BOUNDARY (C0) =',
     &        1P1E13.6,1X,A10,/,26X,
     &        'GROUND-WATER VELOCITY IN X-DIRECTION (VX) =',1P1E13.6,1X,
     &        A10,/,26X,'DISPERSION IN THE X-DIRECTION (DX) =',1P1E13.6,
     &        1X,A10,/,26X,'FIRST-ORDER SOLUTE-DECAY RATE (DK) =',
     &        1P1E13.6,1X,A10,/,26X,
     &        'LENGTH OF FINITE FLOW SYSTEM (XL) =',1P1E13.6,1X,A10,/,
     &        26X,'PLOT SCALING FACTOR (XSCLP) =',1P1E13.6)
 9040 FORMAT ('1',//,16X,'SOLUTE CONCENTRATION AS A FUNCTION OF TIME'
     &        ,15X,'*  INDICATES SOLUTION DID NOT CONVERGE',//,26X,
     &        'TIME VALUES, IN ',A10)
 9045 FORMAT (6X,F12.4,'  ! ',9(F11.5,A1))
      END
C
C
C
      SUBROUTINE CNORM1(XL,T,X,D,V,DK,ROOT,CN,NROOT,IERR)
      INTEGER NROOT
      DOUBLE PRECISION XL, T, X, D, V, DK, ROOT(NROOT), CN
      CHARACTER*1 IERR
C
C       SOLUTION FOR THE ONE-DIMENSIONAL SOLUTE-TRANSPORT EQUATION
C       FOR A SYSTEM OF FINITE LENGTH WITH A FIRST-TYPE SOURCE
C       BOUNDARY CONDITION. VALUE RETURNED IS THE NORMALIZED SOLUTE
C       CONCENTRATION AT A GIVEN X-COORDINATE AND TIME VALUE.
C       FOR NO SOLUTE DECAY, A SIMPLIFIED SOLUTION IS USED.
C
      DOUBLE PRECISION beta, beta2, d2, denom, dkl2d, sigma, u, v2d, 
     &                 vl2d, vl2d2, vmu, vpu, vsqt4d, vupm, vx2d, x1, 
     &                 x2, x3, xl2
      INTEGER n
C
C---INTRINSICS
      DOUBLE PRECISION DEXP, DSIN, DABS, DSQRT
      INTRINSIC DEXP, DSIN, DABS, DSQRT
C
      IERR = ' '
      xl2 = XL*XL
      v2d = V/(2.0D0*D)
      vx2d = v2d*X
      vl2d = v2d*XL
      vl2d2 = vl2d*vl2d
      vsqt4d = V*V*T/(4.0D0*D)
      IF (DK.EQ.0.0D0) THEN
C
C---BEGIN SERIES SUMMATION FOR SOLUTE WITH NO DECAY
        sigma = 0.0D0
        DO 10 n = 1, NROOT
          beta = ROOT(n)
          beta2 = beta*beta
C
C---TERM 1
          denom = beta2 + vl2d2 + vl2d
          x1 = beta*DSIN(beta*X/XL)*DEXP(-beta2*D*T/xl2)
          x1 = x1/denom
          sigma = sigma + x1
          IF (n.GT.25 .AND. DABS(x1).LT.1.0D-14) GOTO 20
   10   CONTINUE
        IERR = '*'
   20   CN = 1.0D0 - 2.0D0*DEXP(vx2d-vsqt4d)*sigma
      ELSE
C
C---BEGIN SERIES SUMMATION FOR SOLUTE WITH DECAY
        sigma = 0.0D0
        DO 30 n = 1, NROOT
          beta = ROOT(n)
          beta2 = beta*beta
C
C---TERM 1
          x1 = (beta2+vl2d2)*DEXP(-beta2*D*T/xl2)
C
C---TERM 2
          dkl2d = DK*XL*XL/D
          denom = (beta2+vl2d2+vl2d)*(beta2+vl2d2+dkl2d)
          x2 = beta*DSIN(beta*X/XL)/denom
          sigma = sigma + x1*x2
C
C---CHECK FOR CONVERGENCE OF SERIES
          IF (n.GT.25 .AND. DABS(x1*x2).LT.1.0D-14) GOTO 40
   30   CONTINUE
        IERR = '*'
C
C---TERM 3
   40   u = DSQRT(V*V+4.0D0*DK*D)
        vmu = V - u
        vpu = V + u
        vupm = (u-V)/vpu
        d2 = D*2.0D0
        x3 = DEXP(vmu*X/d2) + vupm*DEXP((vpu*X-2.0D0*u*XL)/d2)
        x3 = x3/(1.0D0+vupm*DEXP(-u*XL/D))
        CN = x3 - 2.0D0*DEXP(vx2d-vsqt4d-DK*T)*sigma
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE CNORM3(XL,T,X,D,V,DK,ROOT,CN,NROOT,IERR)
      INTEGER NROOT
      DOUBLE PRECISION XL, T, X, D, V, DK, ROOT(NROOT), CN
      CHARACTER*1 IERR
C
C       SOLUTION FOR THE ONE DIMENSIONAL SOLUTE-TRANSPORT EQUATION
C       FOR A SYSTEM OF FINITE LENGTH WITH A THIRD-TYPE SOURCE
C       BOUNDARY CONDITION. VALUE RETURNED IS THE NORMALIZED SOLUTE
C       CONCENTRATION AT A GIVEN X-COORDINATE AND TIME VALUE.
C       FOR NO SOLUTE DECAY, A SIMPLIFIED SOLUTION IS USED.
C
      DOUBLE PRECISION beta, beta2, betaxl, d2, denom, dkl2d, sigma, u, 
     &                 v2d, vl2d, vl2d2, vld, vmu, vpu, vsqt4d, vupm, 
     &                 vx2d, x1, x2, x3, xl2
      INTEGER n
C
C---INTRINSICS
      DOUBLE PRECISION DEXP, DSIN, DABS, DSQRT
      INTRINSIC DEXP, DSIN, DABS, DSQRT
C
      IERR = ' '
      xl2 = XL*XL
      v2d = V/(2.0D0*D)
      vld = V*XL/D
      vx2d = v2d*X
      vl2d = v2d*XL
      vl2d2 = vl2d*vl2d
      vsqt4d = V*V*T/(4.0D0*D)
      IF (DK.EQ.0.0D0) THEN
C
C---BEGIN SERIES SUMMATION FOR SOLUTE WITH NO DECAY
        sigma = 0.0D0
        DO 10 n = 1, NROOT
          beta = ROOT(n)
          beta2 = beta*beta
C
C---TERM 1
          betaxl = beta*X/XL
          x1 = beta*(beta*DCOS(betaxl)+vl2d*DSIN(betaxl))
C
C---TERM 2
          denom = (beta2+vl2d2+vld)*(beta2+vl2d2)
          x2 = DEXP(-beta2*D*T/xl2)/denom
C
          sigma = sigma + x1*x2
          IF (n.GT.25 .AND. DABS(x1*x2).LT.1.0D-14) GOTO 20
   10   CONTINUE
        IERR = '*'
   20   CN = 1.0D0 - 2.0D0*vld*DEXP(vx2d-vsqt4d)*sigma
      ELSE
C
C---BEGIN SERIES SUMMATION FOR SOLUTE WITH DECAY
        sigma = 0.0D0
        DO 30 n = 1, NROOT
          beta = ROOT(n)
          beta2 = beta*beta
C
C---TERM 1
          betaxl = beta*X/XL
          x1 = beta*(beta*DCOS(betaxl)+vl2d*DSIN(betaxl))
C
C---TERM 2
          dkl2d = DK*XL*XL/D
          denom = (beta2+vl2d2+vld)*(beta2+vl2d2+dkl2d)
          x2 = DEXP(-beta2*D*T/xl2)/denom
          sigma = sigma + x1*x2
C
C---CHECK FOR CONVERGENCE OF SERIES
          IF (n.GT.25 .AND. DABS(x1*x2).LT.1.0D-14) GOTO 40
   30   CONTINUE
        IERR = '*'
C
C---TERM 3
   40   u = DSQRT(V*V+4.0D0*DK*D)
        vmu = V - u
        vpu = V + u
        vupm = (u-V)/vpu
        d2 = D*2.0D0
        x3 = DEXP(vmu*X/d2) + vupm*DEXP((vpu*X-2.0D0*u*XL)/d2)
        x3 = 2.0D0*V*x3/(vpu+vmu*vupm*DEXP(-u*XL/D))
        CN = x3 - 2.0D0*vld*DEXP(vx2d-vsqt4d-DK*T)*sigma
      ENDIF
C
      RETURN
      END
C
C
C
      SUBROUTINE ROOT1(C,ROOT,NROOT)
      INTEGER NROOT
      DOUBLE PRECISION ROOT(NROOT), C
C
C       THIS ROUTINE CALCULATES ROOTS OF THE EQUATION: B*COTAN(B)+C=0
C       USING NEWTON'S SECOND-ORDER METHOD.
C
C                      PROGRAM VARIABLES
C       MAXIT    MAXIMUM NUMBER OF ITERATIONS ALLOWED IN ROOT SEARCH
C       EPS      CONVERGENCE CRITERION
C       F1,F2    1ST AND 2ND DERIVATIVES OF THE EQUATION
C       H        SECOND-ORDER CORRECTION FACTOR
C
      DOUBLE PRECISION cotx, eps, f, f1, f2, h, pi, sinx2, x
      INTEGER i, maxit, n
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
C---INTRINSICS
      DOUBLE PRECISION DSIN, DTAN, DABS
      INTRINSIC DSIN, DTAN, DABS
C
      DATA maxit/50/, eps/1.0D-10/
C
C---FIRST ROOT LIES BETWEEN PI/2 AND PI. START WITH .75*PI
      pi = 3.14159265359D0
      ROOT(1) = 0.750D0*pi
C
C---START LOOP FOR EACH ROOT SEARCH
      DO 30 n = 1, NROOT
C
C---BEGIN ITERATIVE LOOP
        DO 10 i = 1, maxit
          x = ROOT(n)
          sinx2 = DSIN(x)*DSIN(x)
          cotx = 1.0D0/DTAN(x)
          f = x*cotx + C
C
C---IF F IS 0.0, EXACT ROOT HAS BEEN FOUND
          IF (f.EQ.0.0D0) GOTO 20
          f1 = cotx - x/sinx2
          f2 = -1.0D0/sinx2 - (sinx2-x*DSIN(x*2.0D0))/(sinx2*sinx2)
          h = (f2/2.0D0)/f1 - f1/f
          h = 1.0D0/h
          ROOT(n) = x + h
C
C---CHECK FOR CONVERGENCE. IF NOT ACHIEVED, RE-ITERATE
          IF (DABS(h).LT.eps) GOTO 20
   10   CONTINUE
        WRITE (Io,9005) maxit, n
        STOP
C
C---NEXT ROOT IS ABOUT PI GREATER THAN LAST ROOT
   20   IF (n.NE.NROOT) ROOT(n+1) = ROOT(n) + pi
   30 CONTINUE
      RETURN
C
C---FORMAT STATEMENTS
 9005 FORMAT (6X,
     &   '**** WARNING ****  ROOT SEARCH ROUTINE DID NOT CONVERGE AFTER'
     &   ,I5,'ITERATIONS WHILE SEARCHING FOR ROOT',I5)
      END
C
C
C
      SUBROUTINE ROOT3(A,C,ROOT,NROOT)
      INTEGER NROOT
      DOUBLE PRECISION ROOT(NROOT), A, C
C
C       THIS ROUTINE CALCULATES ROOTS OF THE EQ: B*COTAN(B)-C*B**2+A=0
C       USING NEWTON'S SECOND-ORDER METHOD.
C
C                      PROGRAM VARIABLES
C       MAXIT    MAXIMUM NUMBER OF ITERATIONS ALLOWED IN ROOT SEARCH
C       EPS      CONVERGENCE CRITERION
C       F1,F2    1ST AND 2ND DERIVATIVES OF THE EQUATION
C       H        SECOND-ORDER CORRECTION FACTOR
C
      DOUBLE PRECISION cotx, eps, f, f1, f2, h, pi, sinx2, x
      INTEGER i, maxit, n
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
C---INTRINSICS
      DOUBLE PRECISION DSIN, DTAN, DABS
      INTRINSIC DSIN, DTAN, DABS
C
      DATA maxit/50/, eps/1.0D-10/
C
C---FIRST ROOT LIES BETWEEN 0.0 AND PI. START WITH 0.5*PI
      pi = 3.14159265359D0
      ROOT(1) = 0.50D0*pi
C
C---START LOOP FOR EACH ROOT SEARCH
      DO 30 n = 1, NROOT
C
C---BEGIN ITERATIVE LOOP
        DO 10 i = 1, maxit
          x = ROOT(n)
          sinx2 = DSIN(x)*DSIN(x)
          cotx = 1.0D0/DTAN(x)
          f = x*cotx - C*x*x + A
C
C---IF F IS 0.0, EXACT ROOT HAS BEEN FOUND
          IF (f.EQ.0.0D0) GOTO 20
          f1 = cotx - x/sinx2 - (2.0D0*C*x)
          f2 = -1.0D0/sinx2 - (sinx2-x*DSIN(x*2.0D0))/(sinx2*sinx2)
     &         - 2.0D0*C
          h = (f2/2.0D0)/f1 - f1/f
          h = 1.0D0/h
          ROOT(n) = x + h
C
C---CHECK FOR CONVERGENCE. IF NOT ACHIEVED, RE-ITERATE
          IF (DABS(h).LT.eps) GOTO 20
   10   CONTINUE
        WRITE (Io,9005) maxit, n
        STOP
C
C---NEXT ROOT IS ABOUT PI GREATER THAN LAST ROOT
   20   IF (n.NE.NROOT) ROOT(n+1) = ROOT(n) + pi
   30 CONTINUE
      RETURN
C
C---FORMAT STATEMENTS
 9005 FORMAT (6X,
     &   '**** WARNING ****  ROOT SEARCH ROUTINE DID NOT CONVERGE AFTER'
     &   ,I5,'ITERATIONS WHILE SEARCHING FOR ROOT',I5)
      END
