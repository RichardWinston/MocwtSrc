C
C       ******************************************************
C      *                                                      *
C      *                  SUBROUTINE EXERFC                   *
C      *                                                      *
C      *            VERSION CURRENT AS OF 10/01/87            *
C      *              CODE CLEANUP 03/27/96 - RSR             *
C      *                                                      *
C       ******************************************************
C
      SUBROUTINE EXERFC(X,YY,Z)
      DOUBLE PRECISION X, YY, Z
C
C       THIS ROUTINE USES RATIONAL CHEBYSHEV APPROXIMATIONS
C       FOR EVALUATING THE ERROR FUNCTION AND COMPLEMENTARY
C       ERROR FUNCTION IN ORDER TO EVALUATE THE PRODUCT OF
C       EXP(X) AND ERFC(Y)
C
C---Local variables
      DOUBLE PRECISION erf, erfcy, p1(5), p2(9), p3(6), q1(5), q2(9), 
     &                 q3(6), sqrtpi, sump, sumq, y, y2i, yi
      INTEGER i
C
C---Intrinsics
      DOUBLE PRECISION DABS, DEXP
      INTRINSIC DABS, DEXP
C
      DATA p1/3.209377589138469472562D03, 3.774852376853020208137D02, 
     &     1.138641541510501556495D02, 3.161123743870565596947D0, 
     &     1.857777061846031526730D-01/
      DATA q1/2.844236833439170622273D03, 1.282616526077372275645D03, 
     &     2.440246379344441733056D02, 2.360129095234412093499D01, 
     &     1.0D0/
      DATA p2/1.23033935479799725272D03, 2.05107837782607146532D03, 
     &     1.71204761263407058314D03, 8.81952221241769090411D02, 
     &     2.98635138197400131132D02, 6.61191906371416294775D01, 
     &     8.88314979438837594118D00, 5.64188496988670089180D-01, 
     &     2.15311535474403846343D-08/
      DATA q2/1.23033935480374942043D03, 3.43936767414372163696D03, 
     &     4.36261909014324715820D03, 3.29079923573345962678D03, 
     &     1.62138957456669018874D03, 5.37181101862009857509D02, 
     &     1.17693950891312499305D02, 1.57449261107098347253D01, 1.0D0/
      DATA p3/ - 6.58749161529837803157D-04, 
     &     -1.60837851487422766278D-02, -1.25781726111229246204D-01, 
     &     -3.60344899949804439429D-01, -3.05326634961232344035D-01, 
     &     -1.63153871373020978498D-02/
      DATA q3/2.33520497626869185443D-03, 6.05183413124413191178D-02, 
     &     5.27905102951428412248D-01, 1.87295284992346047209D00, 
     &     2.56852019228982242072D00, 1.0D0/
C
C---FOR  Y = 0.0
      IF (YY.EQ.0.0D0) THEN
        Z = DEXP(X)
      ELSE
        y = DABS(YY)
        sump = 0.0D0
        sumq = 0.0D0
C
C---FOR 0.0 < Y <= .46875
        IF (y.LE.0.46875D0) THEN
          DO 10 i = 1, 5
            y2i = y**(2*(i-1))
            sump = sump + p1(i)*y2i
            sumq = sumq + q1(i)*y2i
   10     CONTINUE
          erf = y*sump/sumq
          IF (YY.LT.0.0D0) erf = -erf
          erfcy = 1.0D0 - erf
          Z = DEXP(X)*erfcy
C
C---FOR .46875 < Y <= 4.0
        ELSEIF (y.LE.4.0D0) THEN
          DO 20 i = 1, 9
            yi = y**(i-1)
            sump = sump + p2(i)*yi
            sumq = sumq + q2(i)*yi
   20     CONTINUE
          Z = DEXP(X-y*y)*sump/sumq
          IF (YY.LT.0.0D0) Z = 2.0D0*DEXP(X) - Z
C
C---FOR  Y > 4.0
        ELSE
          DO 30 i = 1, 6
            y2i = y**(-2*(i-1))
            sump = sump + p3(i)*y2i
            sumq = sumq + q3(i)*y2i
   30     CONTINUE
          sqrtpi = 0.5641895835477562869481D0
          Z = sqrtpi + sump/(y*y*sumq)
          Z = DEXP(X-y*y)*Z/y
          IF (YY.LT.0.0D0) Z = 2.0D0*DEXP(X) - Z
        ENDIF
      ENDIF
      RETURN
      END
