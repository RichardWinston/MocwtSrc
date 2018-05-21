C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE GLQPTS                   *
C      *                                                      *
C      *            VERSION CURRENT AS OF 10/01/87            *
C      *              CODE CLEANUP 03/27/96 - RSR             *
C      *                                                      *
C       ******************************************************
C
      SUBROUTINE GLQPTS(N)
      INTEGER N
C
      INTEGER iskip, in2, i, nc, k, j1, j
      CHARACTER*1 skip
C
      DOUBLE PRECISION Wn, Zn
      COMMON /GLPTS/ Wn(256), Zn(256)
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
C       THIS ROUTINE READS THE NORMALIZED ROOTS ZN(I) AND WEIGHTS WN(I)
C       OF THE LEGENDRE POLYNOMIALS FROM THE DATA FILE 'GLQ.PTS'
C
C       N IS THE NUMBER OF INTEGRATION POINTS AND CAN ONLY HAVE A
C       VALUE OF EITHER 4,20,60,104,OR 256
C
C---INTRINSICS
      INTEGER MOD
      INTRINSIC MOD
C
      in2 = 77
      OPEN (in2,FILE='GLQ.PTS',STATUS='OLD')
C
C---SKIP LINES IN FILE UNTIL CORRECT COEFFICIENTS ARE REACHED
      iskip = -1
      IF (N.EQ.4) iskip = 7
      IF (N.EQ.20) iskip = 9
      IF (N.EQ.60) iskip = 15
      IF (N.EQ.104) iskip = 31
      IF (N.EQ.256) iskip = 57
      IF (iskip.EQ.-1) THEN
        WRITE (Io,9015)
        STOP
      ENDIF
      DO 10 i = 1, iskip
        READ (in2,9005) skip
   10 CONTINUE
C
C---READ IN ZN(I) AND WN(I), FOUR VALUES PER LINE
      nc = N/8
      IF (MOD(N,8).NE.0) nc = nc + 1
      DO 20 i = 1, nc
        k = (i-1)*8 - 1
        READ (in2,9010) (Zn(k+j*2),j=1,4)
   20 CONTINUE
      DO 30 i = 1, nc
        k = (i-1)*8 - 1
        READ (in2,9010) (Wn(k+j*2),j=1,4)
   30 CONTINUE
C
C---FILL IN THE SYMMETRIC TERMS
      DO 40 j = 2, N, 2
        j1 = j - 1
        Zn(j) = -Zn(j1)
        Wn(j) = Wn(j1)
   40 CONTINUE
      CLOSE (in2)
      RETURN
C
C---FORMAT STATEMENTS
 9005 FORMAT (A1)
 9010 FORMAT (4D20.0)
 9015 FORMAT (/,21X,'*****  ERROR IN ROUTINE GLQPTS  *****',/,21X,
     &        'NO. OF ROOTS SPECIFIED MUST EQUAL 4,20,60,104 OR 256')
      END
