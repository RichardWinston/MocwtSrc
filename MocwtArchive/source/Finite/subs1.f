C
C READXT, WCONHD, PLOTNC for one-dimensional programs (seminf and finite)
C
      SUBROUTINE READXT(NX,NT,XCOORD,TIME)
C
C     Read X-coordinates and times at which solute concentrations
C     are computed.
C
      INTEGER NX, NT, i
      DOUBLE PRECISION XCOORD(NX), TIME(NT)
C
      INCLUDE 'units.inc'
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
      READ (In,9005) (XCOORD(i),i=1,NX)
      WRITE (Io,9010) Lunits
      WRITE (Io,9015) (XCOORD(i),i=1,NX)
C
      READ (In,9005) (TIME(i),i=1,NT)
      WRITE (Io,9020) Tunits
      WRITE (Io,9015) (TIME(i),i=1,NT)
C
      RETURN
C
 9005 FORMAT (8F10.0)
 9010 FORMAT (/,26X,'X-COORDINATES AT WHICH SOLUTE CONCENTRATIONS ',
     &        'WILL BE CALCULATED, IN ',A10,/,26X,78('-')/)
 9015 FORMAT (6X,8F12.4)
 9020 FORMAT (/,26X,
     &    'TIMES AT WHICH SOLUTE CONCENTRATIONS WILL BE CALCULATED, IN '
     &    ,A10,/,26X,70('-')/)
      END
C
C
C
      SUBROUTINE WCONHD(IX,NP1,NP2,TIME)
C
C     Print header for solute concentration tables.
C
      INTEGER IX, NP1, NP2
      DOUBLE PRECISION TIME(*)
C
      INCLUDE 'units.inc'
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
      INTEGER j, MOD
      INTRINSIC MOD
C
      IF (MOD(IX,45).EQ.0) THEN
        IF (IX.NE.0) WRITE (Io,9005) Tunits
        WRITE (Io,9010) (TIME(NP1+j),j=1,NP2)
        WRITE (Io,9015) Cunits, Lunits
      ELSE
        IF (MOD(IX,5).EQ.0) WRITE (Io,9020)
      ENDIF
C
      RETURN
C
 9005 FORMAT ('1',//,16X,
     &        'SOLUTE CONCENTRATION AS A FUNCTION OF TIME =',5X,
     &        '(CONTINUED)',//,26X,'TIME VALUES, IN ',A10)
 9010 FORMAT (21X,9F12.4)
 9015 FORMAT (20X,'*',108('-'),/,5X,'X-COORDINATE,  !',44X,
     &        'SOLUTE CONCENTRATION, IN ',A10,/,5X,'IN ',A10,'  !',/,
     &        20X,'!')
 9020 FORMAT (20X,'!')
      END
C
C
C
      SUBROUTINE PLOTNC(NX,NT,IT,XCOORD,TIME,CXT,C0,XSCLP)
C
C     CONVERT X AND C TO SINGLE PRECISION AND DIVIDE BY C0 TO
C     PLOT NORMALIZED CONCENTRATION PROFILE FOR EACH TIME VALUE.
C
      INCLUDE 'dimens.inc'
      INCLUDE 'units.inc'
C
      INTEGER NX, NT, IT
      DOUBLE PRECISION XCOORD(MAXX), TIME, CXT(MAXX,MAXT), C0
      REAL XSCLP
C
      INTEGER i
      REAL xp(MAXX), cp(MAXX), tp
C
C---EXTERNALS
!      EXTERNAL PLOT1D
C
C---INTRINSICS
      REAL SNGL
      INTRINSIC SNGL
C
      DO 10 i = 1, NX
        xp(i) = SNGL(XCOORD(i))
        cp(i) = SNGL(CXT(i,IT)/C0)
   10 CONTINUE
      tp = SNGL(TIME)
!      CALL PLOT1D(xp,cp,NX,tp,IT,NT,Tunits,Lunits,XSCLP)
      RETURN
      END
