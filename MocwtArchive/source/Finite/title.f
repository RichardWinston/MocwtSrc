C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE TITLE                    *
C      *                                                      *
C      *            VERSION CURRENT AS OF 10/01/87            *
C      *              CODE CLEANUP 03/27/96 - RSR             *
C      *                                                      *
C       ******************************************************
C
C        THIS ROUTINE CREATES A TITLE BOX ON THE FIRST PAGE OF
C        PROGRAM OUTPUT. THE ROUTINE READS AND PRINTS ALL DATA
C        CARDS UNTIL IT ENCOUNTERS AN '=' IN COLUMN 1. THE FIRST 4
C        LINES READ IN ARE ALSO USED AS TITLES ON PLOTS.
C
      SUBROUTINE TITLE
      INTEGER i, nn, nn1, ns, iyr, mo, ida, ihr, mn, ise
      CHARACTER line*60, t1*61, fmt*20
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
      CHARACTER*61 Head
      COMMON /TITLES/ Head(4)
C
      INTEGER LENCHR
      EXTERNAL SYSTIM, SYSDAT, LENCHR
C
      Head(1) = ' '
      Head(2) = ' '
      Head(3) = ' '
      Head(4) = ' '
      t1 = ' '
      CALL SYSDAT(iyr,mo,ida)
      CALL SYSTIM(ihr,mn,ise)
      WRITE (Io,9010)
      DO 20 i = 1, 60
        nn = 0
        READ (In,9005,END=10) line
        IF (line(1:1).EQ.'=') GOTO 30
        t1 = line
C---STRIP OFF TRAILING BLANKS AND CENTER LINE
        nn = LENCHR(line)
   10   nn1 = nn + 1
        t1(nn1:nn1) = '$'
        IF (i.LT.5) Head(i) = t1
        ns = (66-nn)/2
        WRITE (fmt,9020) ns
        WRITE (Io,fmt) line(1:nn), '*'
   20 CONTINUE
   30 WRITE (Io,9015) iyr, mo, ida, ihr, mn, ise
      RETURN
C
C---FORMAT STATEMENTS
 9005 FORMAT (A60)
 9010 FORMAT ('1',//////////,17X,68('*'))
 9015 FORMAT (17X,'*',66X,'*',/,17X,'*',15X,'PROGRAM RUN ON ',I3.2,
     &        2('/',I2.2),' AT',I3.2,2(':',I2.2),15X,'*',/,17X,'*',66X,
     &        '*',/,17X,68('*'),/,'1')
 9020 FORMAT ('(17X,''*''',I2,'X,A,T85,A)')
      END
