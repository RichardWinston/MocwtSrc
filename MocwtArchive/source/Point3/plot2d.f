C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE PLOT2D                   *
C      *                                                      *
C      *            VERSION CURRENT AS OF 11/30/88            *
C      *              CODE CLEANUP 03/27/96 - RSR             *
C      *                                                      *
C       ******************************************************
C
      SUBROUTINE PLOT2D(ZP,TP,NX,NY,IT,NT,IFLG)
      INTEGER NX, NY, IT, NT, IFLG
      REAL TP, ZP
C
C       THIS ROUTINE INITIALIZES A CONTOUR PLOT ON THE RECTANGULAR GRID
C       DEFINED IN THE X-Y PLANE BY THE X AND Y VALUES READ IN. ONE
C       SUBPLOT IS GENERATED FOR EACH Z VALUE AND A NEW PLOT IS
C       GENERATED FOR EACH TIME VALUE.  THE ROUTINE USES DISSPLA (ISCO)
C       SOFTWARE SUBROUTINE CALLS.
C       PLOT SCALING FACTORS (XSCLP,YSCLP) AND CONTOUR INTERVAL (DELTA)
C       ARE SPECIFIED IN THE MAIN PROGRAM.
C       IFLG = 0, means called from POINT2, STRIPF, STRIPI, or GAUSS
C            = 1, means called from POINT3, PATCHF, or PATCHI (3-D pgms)
C
      INCLUDE 'dimens.inc'
      INCLUDE 'pltdat.inc'
      INCLUDE 'units.inc'
C
      REAL hite, x1, xaxis, xmax, xmin, xorig, xpm, xspc, y1, yaxis, 
     &     ymax, ymin, yorig, yp3, ypm, yspc
      INTEGER ipl
      CHARACTER*26 labx, laby
      CHARACTER*36 labx1, laby1
C
      CHARACTER*61 Head
      COMMON /TITLES/ Head(4)
C
C---EXTERNALS
      EXTERNAL DGDISS, PAGE, DONEPL, PHYSOR, AREA2D, HEADIN, INTAXS, 
     &         YAXANG, XREVTK, YREVTK, XNAME, YNAME, GRAF, RESET, 
     &         XNONUM, YNONUM, XGRAXS, YGRAXS, HEIGHT, MESSAG, REALNO, 
     &         CNTOUR, ENDGR, ENDPL
C
C---INTRINSICS
      INTEGER INT
      REAL AMOD
      INTRINSIC INT, AMOD
C
C---CALCULATE PLOT SIZE AND DRAW BORDER
      xspc = 1.5
      IF (IFLG.EQ.1) xspc = 2.0
      yspc = 2.0
      x1 = Xp(NX) - Xp(1)
      xaxis = INT(x1/Xsclp)
      IF (AMOD(x1,Xsclp).GT.0.0) xaxis = xaxis + 1.0
      y1 = Yp(NY) - Yp(1)
      yaxis = INT(y1/Ysclp) + 1.0
      IF (AMOD(y1,Ysclp).GT.0.0) yaxis = yaxis + 1.0
      IF (IT.EQ.1) THEN
        CALL DGDISS
        xpm = (xaxis+xspc)*NT + xspc
        ypm = yaxis + yspc
        CALL PAGE(xpm,ypm)
      ENDIF
C
C---CHOOSE PLOT SIZE BASED ON MAXIMUM COORDINATE VALUES
      xorig = (IT-1)*(xaxis+xspc) + xspc
      yorig = 0.75
      CALL PHYSOR(xorig,yorig)
      CALL AREA2D(xaxis,yaxis)
      IF (IT.EQ.1) THEN
        CALL HEADIN(Head(1),100,1.,4)
        CALL HEADIN(Head(2),100,1.,4)
        CALL HEADIN(Head(3),100,1.,4)
        CALL HEADIN(Head(4),100,1.,4)
      ENDIF
C
C---ROTATE Y VALUES, PUT TICK MARKS ON INSIDE, AND DEFINE AXES LABELS
      CALL INTAXS
      CALL YAXANG(0.)
      CALL XREVTK
      CALL YREVTK
      labx = 'DISTANCE ALONG X-AXIS, IN '
      laby = 'DISTANCE ALONG Y-AXIS, IN '
      labx1 = labx//Lunits
      laby1 = laby//Lunits
C
C---DRAW AND LABEL AXES
      CALL XNAME(labx1,36)
      CALL YNAME(laby1,36)
      xmin = Xp(1)
      ymin = Yp(1)
      xmax = Xsclp*xaxis + xmin
      ymax = Ysclp*yaxis + ymin
      CALL GRAF(xmin,Xsclp,xmax,ymin,Ysclp,ymax)
C
C---DRAW EXTRA AXIS TO CLOSE BOX
      CALL RESET('XREVTK')
      CALL RESET('YREVTK')
      CALL XNONUM
      CALL YNONUM
      CALL XGRAXS(xmin,Xsclp,xmax,xaxis,' ',1,0.0,yaxis)
      CALL YGRAXS(ymin,Ysclp,ymax,yaxis,' ',1,xaxis,0.0)
      CALL RESET('XNONUM')
      CALL RESET('YNONUM')
C
C---PRINT TITLE
      hite = (xaxis-1.0)/(55.*.86)
      IF (hite.GT.0.14) hite = 0.14
      CALL HEIGHT(hite)
      yp3 = yaxis - 0.07 - 1.5*hite
      CALL MESSAG('NORMALIZED CONCENTRATION AT TIME =$',100,0.5,yp3)
      ipl = 3
      IF (AMOD(TP,0.01).EQ.0.0) ipl = 2
      IF (AMOD(TP,0.1).EQ.0.0) ipl = 1
      IF ((TP-INT(TP)).EQ.0.0) ipl = 0
      CALL REALNO(TP,ipl,'ABUT','ABUT')
      CALL MESSAG(Tunits,10,'ABUT','ABUT')
      IF (IFLG.EQ.1) THEN
        yp3 = yp3 - 1.5*hite
        CALL MESSAG('                        AND AT Z =$',100,0.5,yp3)
        ipl = 3
        IF (AMOD(ZP,0.01).EQ.0.0) ipl = 2
        IF (AMOD(ZP,0.1).EQ.0.0) ipl = 1
        IF ((ZP-INT(ZP)).EQ.0.0) ipl = 0
        CALL REALNO(ZP,ipl,'ABUT','ABUT')
        CALL MESSAG(Lunits,10,'ABUT','ABUT')
      ENDIF
C
C---COUNT NUMBER OF DIGITS IN CONTOUR LABEL
      yp3 = yp3 - 1.5*hite
      CALL MESSAG('CONTOUR INTERVAL =$',100,0.5,yp3)
      ipl = 3
      IF (AMOD(Delta,0.01).EQ.0.0) ipl = 2
      IF (AMOD(Delta,0.1).EQ.0.0) ipl = 1
      CALL REALNO(Delta,ipl,'ABUT','ABUT')
      CALL MESSAG('C/Co$',100,'ABUT','ABUT')
C
C---CALL ROUTINE THAT ACTUALLY DOES THE CONTOURING
      CALL CNTOUR(NX,NY,ipl)
C
C---SUBPLOT IS FINISHED
      CALL ENDGR(0)
      IF (IT.EQ.NT) THEN
        CALL ENDPL(0)
        CALL DONEPL
      ENDIF
      RETURN
      END
