


C     *****************************************************************

      FUNCTION SPENCE(X,Y)

C     ------------------------------------------------------------------
C     COMPUTES THE SPENCE FUNCTION
C
C           PHI(X) = - INTEGRAL FROM 0 TO X OF (LOG(1-T))/T
C
C     AS DEFINED IN STEIN, EQ. A48.
C     Y EQUALS 1-X.  IT IS DEFINED OUTSIDE THE ROUTINE IN CASE IT CAN BE
C     COMPUTED MORE EXACTLY (NECESSARY IF X IS NEAR 1).
C
C     THE SPENCE INTEGRAL FOR N=2 (THE DILOGARITHM) IS DEFINED IN
C     ABRAMOWITZ ANDSTEGUN, HANDBOOK OF MATHEMATICAL FUNCTIONS, AS (EQ.
C     27.7.1)
C
C           F(X) = - INTEGRAL FROM 1 TO X OF (LOG(T))/(T-1).
C
C     COMPARISON SHOWS THAT
C
C           PHI(X) = F(1-X)    FOR 0^X^1.
C
C     SO I CAN ADAPT THE SERIES EXPANSION OF F(X) (AS, EQ. 27.7.2)
C           F(X) = SUM FROM 1 TO INFINITY OF ((1-X)**K)/(K**2))
C     VALID FOR 0 < X < 2.
C
C     THE SERIES EXPANSION FOR PHI(X), VALID FOR 0 < X < 1, IS THEREFORE
C
C           PHI(X) = SUM FROM 1 TO INFINITY OF (X**K)/K**2.
C
C     THIS ROUTINE ALSO USES IDENTITY 27.7.3 FROM AS TO INCREASE THE
C     EFFICIENCY OF COMPUTATION FOR LARGE X (I.E., X NEAR 1.)
C     ------------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      DATA MAXIT,EPS/10000,1.D-10/
      DATA PI/3.1415926536D0/

C     ------------------------------------------------------------------
      IF((X .LT. 0.D0) .OR. (X .GT. 1.D0)) STOP 'SPENCE'

      Z = X
      IF(X .GT. .5D0) Z = Y

      SPENCE = Z
      IF(Z .EQ. 0.D0) GO TO 30

      DO 10 K = 2,MAXIT
      ADD = (Z**K)/(K*K)
      SPENCE = SPENCE + ADD
      CHANGE = DABS(ADD/SPENCE)
      IF(CHANGE .LT. EPS) GO TO 20
10    CONTINUE
20    CONTINUE
      IF(X .LE. .5D0) RETURN

C     -----------------------------------------------------------------
C     USE IDENTITY 27.7.3 FROM AS FOR LARGE X (.5{X^1).
C     -----------------------------------------------------------------
      SPENCE = -SPENCE - DLOG(X)*DLOG(Y) + PI*PI/6.D0
      RETURN

C     -----------------------------------------------------------------
C     SPECIAL CASES_  X=0 OR X=1.
C     -----------------------------------------------------------------
30    IF(X .EQ. 1.D0) SPENCE = PI*PI/6.D0

      RETURN
      END
