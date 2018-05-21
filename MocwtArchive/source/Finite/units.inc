C
C           NOTE:  ANY CONSISTANT SET OF UNITS MAY BE USED IN THE
C           MODEL.  NO FORMAT STATEMENTS NEED TO BE CHANGED AS
C           LABELS FOR ALL VARIABLES ARE SPECIFIED IN MODEL INPUT.
C
C     CHARACTER VARIABLES USED TO SPECIFY UNITS FOR MODEL PARAMETERS
C
C     CUNITS   UNITS OF CONCENTRATION (M/L**3)
C     VUNITS   UNITS OF GROUND-WATER VELOCITY (L/T)
C     DUNITS   UNITS OF DISPERSION COEFFICIENT (L**2/T)
C     KUNITS   UNITS OF SOLUTE DECAY CONSTANT  (1/T)
C     LUNITS   UNITS OF LENGTH (L)
C     TUNITS   UNITS OF TIME (T)
C     QUNITS   UNITS OF FLUID INJECTION RATE (L**3/T)
C
C     VERSN    is used to encode a version number into executables that
C              is accesses with the UNIX what and RCS ident commands
C
      CHARACTER*10 Cunits, Vunits, Dunits, Kunits, Lunits, Tunits,
     &             Qunits
      CHARACTER*80 Versn
      COMMON /UNITS/ Versn, Cunits, Vunits, Dunits, Kunits, Lunits, 
     &               Tunits, Qunits
