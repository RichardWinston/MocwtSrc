      SUBROUTINE DGDISS
C
C     Routine to designate and initiate graphics output device
C     for use with CA/Disspla software
C
      INTEGER ixarg(10), ipo, ibuf(16), jbuf(16)
      CHARACTER*64 psfile
      EQUIVALENCE (jbuf(1),psfile)
C
      CHARACTER Psname*64, Pcolor*16
      COMMON /PLTNAM/ Psname, Pcolor
C
C     Disspla routines called
      EXTERNAL IOMGR, PSCRPT, TAB, TK41, XWNDOW, COMPRS, SETCLR
C
      INTEGER LENCHR
      EXTERNAL LENCHR
C
      ixarg(1) = 112
      ixarg(2) = 1
      ixarg(3) = 1
      ixarg(4) = 6
      ixarg(5) = 6
      ixarg(6) = 0
      ixarg(7) = 3
      ixarg(8) = -1
      ixarg(9) = 0
      ixarg(10) = 0
C
   10 WRITE (*,9005)
      READ (*,'(I1)',ERR=20) ipo
      IF (ipo.EQ.0) THEN
        ibuf(1) = 5
        CALL IOMGR(ibuf,-102)

        WRITE (*,9015) Psname(:LENCHR(Psname))
        READ (*,'(A)') psfile
        IF (psfile(:4).EQ.'    ') psfile = Psname
        CALL IOMGR(jbuf,-103)
        ibuf(1) = 1
        CALL IOMGR(ibuf,-104)
        CALL PSCRPT(0,0,0)
      ELSEIF (ipo.EQ.1) THEN
        CALL XWNDOW(1,ixarg,9)
      ELSEIF (ipo.EQ.2) THEN
        CALL COMPRS
      ELSEIF (ipo.EQ.3) THEN
        CALL TAB
      ELSEIF (ipo.EQ.4) THEN
        CALL TK41(4114)
      ELSE
        GOTO 20
      ENDIF
      WRITE(*,9020) 'GREEN'
      READ (*,'(A)') Pcolor
      IF (Pcolor(:4).EQ.'    ') Pcolor = 'GREEN'
      CALL SETCLR (Pcolor(1:LENCHR(Pcolor)))
      RETURN
   20 WRITE (*,9010)
      GOTO 10
C
 9005 FORMAT (//,10X,'     SELECT OUTPUT DEVICE',/,10X,
     &        '----------------------------------',/,10X,
     &        '0- POSTSCRIPT FILE'/10X,'1- SCREEN OUTPUT',/,10X,
     &        '2- METAFILE - popfil.dat',/,10X,
     &        '3- TAB 132/15 GRAPHICS',/,10X,'4- TEK 4114',/)
 9010 FORMAT (//,10X,'***Invalid value entered, please try again***')
 9015 FORMAT (//,5X,'ENTER POSTSCRIPT FILE NAME, DEFAULT= ',A)
 9020 FORMAT (//,5X,'ENTER GRAPH COLOR NAME, DEFAULT= ',A)
      END
