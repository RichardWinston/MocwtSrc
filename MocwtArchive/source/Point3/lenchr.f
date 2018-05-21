      INTEGER FUNCTION LENCHR
     I                       ( STRING )
C
C     + + + PURPOSE + + +
C     Get length of a character string
C
C     + + + DUMMY ARGUMENTS + + +
      CHARACTER*(*) STRING
C
C     + + + ARGUMENT DEFINITIONS + + +
C     STRING - Input string to determine length of
C
C     + + + LOCAL VARIABLES + + +
      INTEGER max, i
C
C     + + + INTRINSICS + + +
      INTEGER LEN
      INTRINSIC LEN
C
C     + + + END SPECIFICATIONS + + +
C
      max = LEN(STRING)
      DO 10 i = max, 1, -1
        IF (STRING(i:i).NE.' ') THEN
          LENCHR = i
          GOTO 20
        ENDIF
   10 CONTINUE
      LENCHR = 0
C
   20 RETURN
      END
