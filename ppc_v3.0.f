C
C
      SUBROUTINE PPC(DDPP,ISOLNFLAG,BB,HKR,HKZ,SS,QQ,RW,ZPD,ZPL)

C      ********************************************************
C      *                                                      *
C      *                    **** PPC ****                     *
C      *                                                      *
C      *      COMPUTER PROGRAM FOR CALCULATING DRAWDOWN       *
C      *                                                      *
C      *      IN A CONFINED AQUIFER WITH AXIAL-SYMMETRIC      *
C      *                                                      *
C      *           FLOW TO A PARTIALLY PENETRATING,           *
C      *                                                      *
C      *         INFINITESIMAL-DIAMETER PUMPED WELL           *
C      *                                                      *
C      *          VERSION 3.0 CURRENT AS OF 01/25/08          *
C      *                                                      *
C      ********************************************************
C
C
C      Although this program has been used by the U.S. Geological
C      Survey, no warranty, expressed or implied, is made by the 
C      USGS as to the accuracy and functioning of the program and
C      related program material, nor shall the fact of distribution
C      constitute any such warranty, and no responsibility is
C      assumed by the USGS in connection therewith
C
C---SPECIFICATIONS
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z)
      DIMENSION GAMMA(5)
C
C---COMMON STATEMENTS
      COMMON /PAR1/ IPWD,IRUN,IPWS,NOBWC,IOWS,IDPR
      COMMON /PAR2/ NGAMMA,IDRA,NS,KK,NMAX,NTMS
      COMMON /PAR6/ BETAW,SIGMA,GAMMA
      COMMON /PAR7/ RERRNR,RERRSUM,TDLAST,TLAST
      COMMON /PAR8/ R,ZP,Z1,Z2,WDP
      COMMON /PAR9/ V(20),XLN2,EXPMAX
      COMMON /PAR10/ RD,ZD,ZD1,ZD2
      COMMON /PAR11/ XLD,XDD,WD,SW
C
c      OPEN(UNIT=94,FILE='PPC.out',STATUS='OLD')
C
C
C---THE FOLLOWING PARAMETERS NEED TO BE PASSED FROM MAIN PROGRAM
C      BB=AQUIFER THICKNESS
C      HKR=HORIZONTAL K (L/T)
C      HKZ=VERTICAL K (L/T)
C      SS=SPECIFIC STORAGE (1/L)
C
C      QQ=PUMPING RATE OF WELL (L**3/T)
C      RW=RADIUS OF THE SCREENED INTERVAL OF PUMPED WELL (L)
C      ZPD=DEPTH BELOW TOP OF AQUIFER TO THE TOP OF THE SCREENED 
C          INTERVAL OF THE PUMPED WELL (L)
C      ZPL=DEPTH BELOW TOP OF AQUIFER TO THE BOTTOM OF THE SCREENED
C          INTERVAL OF THE PUMPED WELL (L
C
C      TIMEDD=TIME FOR WHICH DRAWDOWNS WILL BE CALCULATED (T)--MUST BE
C       GREATER THAN 0.0D0
C
c  input for uncoupled code
c      QQ=96000.0
c      RW=0.99
c      ZPD=33.3 
c      ZPL=66.7
C     TIMEDD=100.
c      BB=100.0
c      HKR=140.0
c      HKZ=140.0
c      SS=0.000002
C
C   NOTE, PROGRAM ASSUMES THAT THE PUMPED WELL IS PARTIALLY PENETRATING.
C    THEREFORE, NEED A TEST IN THE MAIN PROGRAM THAT ENSURES THAT THE
C    PUMPED WELL IS PARTIALLY PENETRATING. IF THE PUMPED WELL IS FULLY
C    PENETRATING, DO NOT CALL THIS SUBROUTINE.
C
C1--SET PARAMETERS
C    CONFINED AQUIFER (IAQ=0)
C    DIMENSIONAL ANALYSIS (IFORMAT=1)
C    NO DRAINAGE FROM WT (IDRA=0), NALPHA=0
C    USER-SPECIFIED TIMES (ITS=1) AND NO MEASURED DATA (IMEAS=0)
C    BECAUSE ITS=1, TLAST=0.0D0, NLC=0, AND NOX=0
      IAQ=0
      IFORMAT=1
      IDRA=0
      NALPHA=0
      ITS=1
      TLAST=0.0D0
      NLC=0
      NOX=0
      IMEAS=0
C   PROGRAM SOLUTION VARIABLES
      RERRNR=0.0D0
      RERRSUM=1.D-07
      NMAX=200
      NTMS=0
      NS=8
C
C   PUMPED-WELL INFORMATION
C    WELL IS PARTIALLY PENETRATING (IPWS=0)
C    WELL HAS INFINITESIMAL DIAMETER (IPWD=0); THAT IS, NO WELLBORE 
C      STORAGE
C    WELL-BORE SKIN IS NOT ACCOUNTED FOR HERE (SW=0.0D0)
      IPWS=0
      IPWD=0
      RC=0.0D0
      SW=0.0D0
C    DRAWDOWNS WILL BE CALCULATED FOR A SINGLE TIME (NTSPW=1; IRUN=1)
      NTSPW=1
      IRUN=1
C      
C   OBSERVATION-WELL INFORMATION
C    NO DRAWDOWN CALCULATIONS ARE MADE FOR OBSERVATION WELLS
      NOBWC=0
C
C   CALCULATE AQUIFER PARAMETERS
       AT=HKR*BB
       XKD=HKZ/HKR
       ASC=SS*BB
       SIGMA=0.0D0
C12d-CALCULATE PUMPING-WELL DIMENSIONLESS PARAMETERS 
      RWD=RW/BB
      BETAW=XKD*(RWD*RWD)
      IF(IPWD.EQ.0)WD=0.0D0
      IF(IPWS.EQ.0)THEN
       XDD=ZPD/BB
       XLD=ZPL/BB
      ENDIF
C
C3--DEFINE SELECTED PROGRAM PARAMETERS
      PI=3.141592653589793D0
      IF(IFORMAT.GE.1)THEN
       F1=(SS*RW*RW)/HKR
       F2=QQ/(4.0D0*PI*HKR*BB)
      ENDIF
C
C---EXPMAX IS THE MAXIMUM ALLOWABLE ABSOLUTE VALUE OF EXPONENTIAL ARGUMENTS
       EXPMAX=708.D0
C
C4--CALL LINVST TO CALCULATE COEFFICIENTS USED FOR THE STEHFEST
C     ALGORITHM
       CALL LINVST(V,NS)
       XLN2=DLOG(2.D0)
C
C7--SET KK=1 (NECESSARY FOR SUBROUTINE LTST2):
      KK=1
C
C7a-CACULATE DIMENSIONLESS VARIABLES TO PASS TO LAPLACE TRANSFORM
C    SOLUTION SUBROUTINES 
C
       IOWS=2
       RD=1.0D0
       RDSQ=RD*RD
       ZD=0.0D0
       WDP=0.0D0
C
C7b-DEFINE DIMENSIONLESS TIME (TD
       IF(IFORMAT.GE.1)RDSQ=1.0D0
C
C      DO 30 NT=1,NTS
C       HD=0.0D0
C
C---DETERMINE DIMENSIONLESS TIME OF CURRENT TIME STEP.
C      NTT=NTS-NT+1
C      TD=TIMEDD/F1
C
C---CALCULATE DRAWDOWNS
       DDPPOLD=0.0D0
       EPSILON=0.00001D0
       ISOLNFLAG=0
       TD=1.0D4
C       TD=1.0D-2
       DO 30 NT=1,10
        TD=TD*10.0D0
C---CALCULATE DIMENSIONLESS DRAWDOWN FOR PARTIALLY PENETRATING
C     WELL, CALL LTST2
         CALL LTST2(TD,HD)
         IF(HD.LT.1.D-02)HD=0.0D0
C
         RTD=F1*TD
         RHDPP=F2*HD
C
C---CALCULATE DIMENSIONLESS DRAWDOWN FOR FULLY PENETRATING
C     WELL, CALL LTST1
         CALL LTST1(TD,HD)
         IF(HD.LT.1.D-02)HD=0.0D0
C
         RHDFP=F2*HD
C
C---CALCULATE DRAWDOWN DUE TO PARTIALLY PENETRATING WELL
         DDPP=RHDPP-RHDFP
C
C---WRITE RESULTS
C
         IF(DABS(DDPP).LT.0.001D0)THEN
          DDPPOLD=0.0D0
          GO TO 30
         ENDIF
C
         DDTEST=(DABS(DDPP-DDPPOLD))/DABS(DDPP)
         IF(DDTEST.LT.EPSILON)THEN
          ISOLNFLAG=1
          GO TO 10
         ELSE
          DDPPOLD=DDPP
          ISOLNFLAG=-1
         ENDIF
C---END TIME LOOP FOR DRAWDOWN CALCULATIONS
  30   CONTINUE
C
  10  CONTINUE
C
C---FORMAT STATEMENTS
   14 FORMAT(3X,' DIMENSIONLESS TIME = ',D12.1,' TIME= ',D12.5,
     2' PPHD= ',D12.5,' FPHD= ',D12.5,' DDPP= ',D12.5)
   15 FORMAT(3X,' VALUE RETURNED TO MAIN PROGRAM (DDPP) = ',D12.5)
      RETURN
      END
C
C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE LINVST                   *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C15-SUBROUTINE LINVST CALCULATES COEFFICIENTS USED FOR THE
C   STEHFEST ALGORITHM
C
       SUBROUTINE LINVST(V,NS)
C---SPECIFICATIONS
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
       DIMENSION G(20),V(20),HS(20)
C
       G(1)=1.D0
       NH=NS/2
        DO 1 IS=2,NS
1      G(IS)=G(IS-1)*IS
       HS(1)=2.D0/G(NH-1)
        DO 3 IS=2,NH
       FI=IS
       IF(IS.EQ.NH) GO TO 2
       HS(IS)=FI**(NH)*G(2*IS)/(G(NH-IS)*G(IS)*G(IS-1))
       GO TO 3
2      HS(IS)=FI**(NH)*G(2*IS)/(G(IS)*G(IS-1))
3       CONTINUE
       SN=2*(NH-NH/2*2)-1
        DO 4 IS=1,NS
       V(IS)=0.D0
       K1=(IS+1)/2
       K2=IS
       IF(K2.GT.NH)K2=NH
        DO 5 KS=K1,K2
       IF(2*KS-IS.EQ.0) GO TO 6
       IF(IS.EQ.KS)GO TO 7
       V(IS)=V(IS)+HS(KS)/(G(IS-KS)*G(2*KS-IS))
       GO TO 5
6      V(IS)=V(IS)+HS(KS)/(G(IS-KS))
       GO TO 5
7      V(IS)=V(IS)+HS(KS)/G(2*KS-IS)
5       CONTINUE
       V(IS)=SN*V(IS)
       SN=-SN
4       CONTINUE
       RETURN
       END
C
C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE LTST1                    *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C16-SUBROUTINE LTST1 CALCULATES THE LAPLACE TRANSFORM SOLUTION FOR
C   DRAWDOWN FOR FLOW TO A FULLY PENETRATING WELL OF INFINITESIMAL
C   DIAMETER IN A CONFINED AQUIFER (THEIS SOLUTION).
C
       SUBROUTINE LTST1(TD,HDT)
C---SPECIFICATIONS
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C---COMMON STATEMENTS
      COMMON /PAR2/ NGAMMA,IDRA,NS,KK,NMAX,NTMS
      COMMON /PAR9/ V(20),XLN2,EXPMAX
      COMMON /PAR10/ RD,ZD,ZD1,ZD2
C
       XP=0.D0
      DO 1 I=1,NS
       PP=XLN2*I/TD
C
       CA=RD*DSQRT(PP)
       IF(CA.GT.EXPMAX) CA=EXPMAX
       RE0=BESSK0(CA)
       PDL=RE0/PP
1     XP=XP+V(I)*PDL
       HDT=2.D0*XP*XLN2/TD
C
C---RETURN TO MAIN PROGRAM
       RETURN
       END
C
C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE LTST2                    *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C17-SUBROUTINE LTST2 CALCULATES THE LAPLACE TRANSFORM SOLUTION FOR
C   DRAWDOWN FOR FLOW TO A FINITE DIAMETER, PARTIALLY PENETRATING
C   WELL IN A CONFINED AQUIFER (MODIFIED SOLUTION OF DOUGHERTY AND
C   BABU, 1984). DELAYED DRAWDOWN RESPONSE AT OBSERVATION WELLS
C   IS INCLUDED.
C
       SUBROUTINE LTST2(TD,HD)
C---SPECIFICATIONS
       IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION GAMMA(5)
C---COMMON STATEMENTS
      COMMON /PAR1/ IPWD,IRUN,IPWS,NOBWC,IOWS,IDPR
      COMMON /PAR2/ NGAMMA,IDRA,NS,KK,NMAX,NTMS
      COMMON /PAR6/ BETAW,SIGMA,GAMMA
      COMMON /PAR7/ RERRNR,RERRSUM,TDLAST,TLAST
      COMMON /PAR8/ R,ZP,Z1,Z2,WDP
      COMMON /PAR9/ V(20),XLN2,EXPMAX
      COMMON /PAR10/ RD,ZD,ZD1,ZD2
      COMMON /PAR11/ XLD,XDD,WD,SW
C
       HD=0.D0
       IF(IRUN.EQ.0.AND.KK.EQ.1) RETURN
C
       PI=3.141592653589793D0
C
       IF(IPWS.EQ.1) THEN
        XDD=0.D0
        XLD=1.D0
       ENDIF
C
       XP=0.0D0
C
      DO 1 I=1,NS
       PP=XLN2*I/TD
       Q0=DSQRT(PP)
       Q0RD=Q0*RD
       IF(Q0.GT.EXPMAX) Q0=EXPMAX
       IF(Q0RD.GT.EXPMAX) Q0RD=EXPMAX
       RE0=BESSK0(Q0)
       RE1=BESSK1(Q0)
       RE0X=BESSK0(Q0RD)
       A0=RE0*(XLD-XDD)/(Q0*RE1)
       E0=RE0X*(XLD-XDD)/(Q0*RE1)
       A=0.D0
       E=0.D0
       IF(IPWS.EQ.1) GOTO 30
       IF(IOWS.EQ.1) GOTO 30
       SUMA=0.D0
       SUME=0.D0
C
       NNN=0
C
10     NNN=NNN+1
       IF(NNN.GE.NMAX) GOTO 40
       SUMTA=SUMA
       SUMTE=SUME
       XNPI=NNN*PI
       QN=DSQRT(BETAW*XNPI*XNPI+PP)
       IF(QN.GT.EXPMAX) QN=EXPMAX
       DB=DSIN(XNPI*(1.0D0-XDD))
       DA=DSIN(XNPI*(1.0D0-XLD))
       IF(IPWS.EQ.1) DA=0.D0
       SINES=DB-DA
       RE0=BESSK0(QN)
       RE1=BESSK1(QN)
       XNUM=RE0*SINES*SINES/(XNPI*(XLD-XDD))
       XDEN=0.5D0*QN*RE1*XNPI
       A=XNUM/XDEN
       SUMA=SUMTA+A
C
       IF(KK.GT.1)THEN
        QNRD=QN*RD
        IF(QNRD.GT.EXPMAX) QNRD=EXPMAX
        RE0X=BESSK0(QNRD)
        IF(IOWS.EQ.0) 
     1   XNUM=RE0X*SINES*(DSIN(XNPI*ZD2)
     2   -DSIN(XNPI*ZD1))/(XNPI*(ZD2-ZD1))
        IF(IOWS.EQ.2) 
     1   XNUM=RE0X*SINES*DCOS(XNPI*ZD)
        E=XNUM/XDEN
        SUME=SUMTE+E
       ENDIF
C
       IF(IPWS.EQ.0.AND.NNN.LT.25) GOTO 10
       ERRA=DABS(SUMTA-SUMA)
       IF(KK.EQ.1)THEN
        IF(ERRA.LT.RERRSUM*SUMA) GOTO 40
       ENDIF
       IF(KK.GT.1)THEN
        ERRE=DABS(SUMTE-SUME)
        IF(ERRA.LT.RERRSUM*SUMA.AND.ERRE.LT.RERRSUM*SUME) GOTO 40
       ENDIF
C
       GOTO 10
C
40     CONTINUE
C
       A=SUMA
30     DENOM=(1.D0+WD*PP*(A0+A+SW))
       IF(KK.EQ.1) PDL=(A0+A+SW)/(PP*DENOM)
       IF(KK.GT.1) THEN
         E=SUME
         IF(IDPR.EQ.0) PDL=(E0+E)/(PP*DENOM)
         IF(IDPR.EQ.1) THEN
            SLUGF=1.D0/(1.D0+WDP*PP)
            PDL=SLUGF*(E0+E)/(PP*DENOM)
         ENDIF
       ENDIF
C
      XP=XP+V(I)*PDL
C
1     CONTINUE
C
       HD=2.D0*XP*XLN2/(TD*(XLD-XDD))
C
C      IF(NNN.GE.NMAX) WRITE(IO,100)
C
C---FORMAT STATEMENTS
C100   FORMAT('PROGRAM CONTINUES TO NEXT TIME STEP BECAUSE NNN',
C    2' EXCEEDS NMAX.')
C
C---RETURN TO MAIN PROGRAM
       RETURN
       END
C
C       ******************************************************
C      *                                                      *
C      *                    FUNCTION BESSK0                   *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C19-FUNCTION BESSK0 CALCULATES THE ZERO-ORDER MODIFIED BESSEL FUNCTION
C    OF THE SECOND KIND. SOURCE: PRESS AND OTHERS (1992).
      DOUBLE PRECISION FUNCTION BESSK0(X)
C---SPECIFICATIONS
      DOUBLE PRECISION X,Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,BESSI0
      DATA P1,P2,P3,P4,P5,P6,P7/-0.57721566D0,0.42278420D0,0.23069756D0,
     *    0.3488590D-1,0.262698D-2,0.10750D-3,0.74D-5/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,-0.7832358D-1,0.2189568D-1,
     *    -0.1062446D-1,0.587872D-2,-0.251540D-2,0.53208D-3/
C
      IF (X.LE.2.D0) THEN
        Y=X*X/4.D0
        BESSK0=(-DLOG(X/2.D0)*BESSI0(X))+(P1+Y*(P2+Y*(P3+
     *        Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=(2.D0/X)
        BESSK0=(DEXP(-X)/DSQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *        Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
C
C---RETURN TO MAIN PROGRAM
      RETURN
      END
C
C       ******************************************************
C      *                                                      *
C      *                    FUNCTION BESSI0                   *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C20-FUNCTION BESSI0 CALCULATES THE ZERO-ORDER MODIFIED BESSEL FUNCTION 
C    OF THE FIRST KIND. SOURCE: PRESS AND OTHERS (1992).
      DOUBLE PRECISION FUNCTION BESSI0(X)
C---SPECIFICATIONS
      DOUBLE PRECISION X,Y,AX,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,3.5156229D0,3.0899424D0,1.2067492D
     *0,
     *    0.2659732D0,0.360768D-1,0.45813D-2/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,0.1328592D-1,
     *    0.225319D-2,-0.157565D-2,0.916281D-2,-0.2057706D-1,
     *    0.2635537D-1,-0.1647633D-1,0.392377D-2/
C
      IF (DABS(X).LT.3.75D0) THEN
        Y=(X/3.75D0)**2
        BESSI0=P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7)))))
      ELSE
        AX=DABS(X)
        Y=3.75D0/AX
        BESSI0=(DEXP(AX)/DSQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4
     *      +Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
C
C---RETURN TO MAIN PROGRAM
      RETURN
      END
C
C       ******************************************************
C      *                                                      *
C      *                    FUNCTION BESSK1                   *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C21-FUNCTION BESSK1 CALCULATES THE FIRST-ORDER MODIFIED BESSEL FUNCTION 
C    OF THE SECOND KIND. SOURCE: PRESS AND OTHERS (1992).
      DOUBLE PRECISION FUNCTION BESSK1(X)
C---SPECIFICATIONS
      DOUBLE PRECISION X,Y,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,BESSI1
      DATA P1,P2,P3,P4,P5,P6,P7/1.0D0,0.15443144D0,-0.67278579D0,
     *    -0.18156897D0,-0.1919402D-1,-0.110404D-2,-0.4686D-4/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7/1.25331414D0,0.23498619D0,-0.3655620D-1,
     *    0.1504268D-1,-0.780353D-2,0.325614D-2,-0.68245D-3/
C
      IF (X.LE.2.0) THEN
        Y=X*X/4.0
        BESSK1=(LOG(X/2.0)*BESSI1(X))+(1.0/X)*(P1+Y*(P2+
     *      Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        Y=2.0/X
        BESSK1=(EXP(-X)/SQRT(X))*(Q1+Y*(Q2+Y*(Q3+
     *      Y*(Q4+Y*(Q5+Y*(Q6+Y*Q7))))))
      ENDIF
C
C---RETURN TO MAIN PROGRAM
      RETURN
      END
C
C       ******************************************************
C      *                                                      *
C      *                    FUNCTION BESSI1                   *
C      *                                                      *
C      *          VERSION 1.0 CURRENT AS OF 10/01/99          *
C      *                                                      *
C       ******************************************************
C
C22-FUNCTION BESSI1 CALCULATES THE FIRST-ORDER MODIFIED BESSEL FUNCTION 
C    OF THE FIRST KIND. SOURCE: PRESS AND OTHERS (1992).
      DOUBLE PRECISION FUNCTION BESSI1(X)
C---SPECIFICATIONS
      DOUBLE PRECISION X,Y,AX,P1,P2,P3,P4,P5,P6,P7,
     *    Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9
      DATA P1,P2,P3,P4,P5,P6,P7/0.5D0,0.87890594D0,0.51498869D0,
     *    0.15084934D0,0.2658733D-1,0.301532D-2,0.32411D-3/
      DATA Q1,Q2,Q3,Q4,Q5,Q6,Q7,Q8,Q9/0.39894228D0,-0.3988024D-1,
     *    -0.362018D-2,0.163801D-2,-0.1031555D-1,0.2282967D-1,
     *    -0.2895312D-1,0.1787654D-1,-0.420059D-2/
C
      IF (ABS(X).LT.3.75) THEN
        Y=(X/3.75)**2
        BESSI1=X*(P1+Y*(P2+Y*(P3+Y*(P4+Y*(P5+Y*(P6+Y*P7))))))
      ELSE
        AX=ABS(X)
        Y=3.75/AX
        BESSI1=(EXP(AX)/SQRT(AX))*(Q1+Y*(Q2+Y*(Q3+Y*(Q4+
     *      Y*(Q5+Y*(Q6+Y*(Q7+Y*(Q8+Y*Q9))))))))
      ENDIF
C
C---RETURN TO MAIN PROGRAM
      RETURN
      END
