C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE OFILE                    *
C      *                                                      *
C      *            VERSION CURRENT AS OF 10/01/87            *
C      *              CODE CLEANUP 03/27/96 - RSR             *
C      *                                                      *
C       ******************************************************
C
      SUBROUTINE OFILE
      CHARACTER*64 ifname, ofname, temp
      CHARACTER*1 star
      INTEGER i, n, idot
C
      INTEGER In, Io
      COMMON /IOUNIT/ In, Io
C
      CHARACTER Psname*64, Pcolor*16
      COMMON /PLTNAM/ Psname, Pcolor
C
      INTEGER LENCHR
      EXTERNAL LENCHR
C
      DATA star/'*'/
C
      In = 15
      Io = 16
   10 WRITE (*,9005)
      READ (*,9015) ifname
      OPEN (In,FILE=ifname,STATUS='OLD')
      n = LENCHR(ifname)
      IF (n.LT.1) GOTO 10
      DO 20 i = n, 1, -1
        IF (ifname(i:i).EQ.'.') THEN
          idot = i - 1
          GOTO 30
        ENDIF
   20 CONTINUE
      idot = 0
   30 IF (idot.GT.0 .AND. n-idot.LT.5) n = idot
      IF (n.GT.76) n = 76
      ofname = ifname(1:n)//'.prt'
      Psname = ifname(1:n)//'.ps'
      WRITE (*,9010) ofname(1:n+4)
      READ (*,9015) temp
      IF (temp(:4).NE.'    ') ofname = temp
      IF (ofname(1:1).EQ.star) THEN
        Io = 1
      ELSE
        OPEN (Io,FILE=ofname)
      ENDIF
      RETURN
C
C---FORMAT STATEMENTS
 9005 FORMAT (//,5X,'ENTER INPUT FILE NAME')
 9010 FORMAT (//,5X,'ENTER OUTPUT FILE NAME, DEFAULT= ',A)
 9015 FORMAT (A)
      END
