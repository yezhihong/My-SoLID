C     *****************************************************************

      FUNCTION TERP(XIN,ICHANN)

C     ----------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON /ROSEXEPL/XH(400),YH(400),XH2(400),YH2(400)
      COMMON /ROSECARDS/NPTS
      DIMENSION X(400),Y(400),DELTA(10),A(10)

C     -----------------------------------------------------------------
C     NTERMS-POINT INTERPOLATION ROUTINE FROM BEVINGTON (PROGRAM 13-2),
C     SLIGHTLY MODIFIED.
C     PLACE PROPER CHANNEL INTO X AND Y.
C     -----------------------------------------------------------------

      NTERMS = 3
      IF(ICHANN .NE. 1) GO TO 20
      DO 10 I = 1,NPTS
      X(I) = XH(I)
      Y(I) = YH(I)
10    CONTINUE
      GO TO 40
20    DO 30 I = 1,NPTS
      X(I) = XH2(I)
      Y(I) = YH2(I)
30    CONTINUE
40    CONTINUE

C     ---------------------------------------------------------------
C     CHECK LIMITS
C     ---------------------------------------------------------------
      IF((XIN .GE. X(1)) .AND. (XIN .LE. X(NPTS))) GO TO 60
      TERP = 0.D0
      RETURN

C     ---------------------------------------------------------------
C     SEARCH FOR APPROPRIATE VALUE OF X1.  AS SOON AS X(I) IS GREATER
C     THAN OR EQUAL TO XIN, SET I1 AND EXIT LOOP.
C     ---------------------------------------------------------------
60    CONTINUE
      DO 90 I = 1,NPTS
      IF(XIN-X(I)) 70,80,90

C     ---------------------------------------------------------------
C     XIN IS LESS THAN X(I).
C     ---------------------------------------------------------------
70    I1 = I - NTERMS/2
      IF(I1 .LE. 0) I1 = 1
      GO TO 100

C     ---------------------------------------------------------------
C     XIN EQUALS X(I)
C     ---------------------------------------------------------------
80    TERP = Y(I)
      RETURN

C     ---------------------------------------------------------------
C     XIN IS GREATER THAN X(I).
C     IF EXIT LOOP WITHOUT SUCCESS, USE THE LAST NTERMS POINTS IN THE
C     ARRAY.  (WITH THE INITIAL TEST, THIS WILL NOT HAPPEN).
C     ---------------------------------------------------------------
90    CONTINUE
      I1 = NPTS - NTERMS + 1
100   I2 = I1 + NTERMS - 1
      IF(I2 .LE. NPTS) GO TO 110
      I2 = NPTS
      I1 = I2 - NTERMS + 1
      IF(I1 .GE. 1) GO TO 110
      I1 = 1
      NTERMS = I2 - I1 + 1

C     ---------------------------------------------------------------
C     EVALUATE DEVIATIONS DELTA
C     ---------------------------------------------------------------
110   DENOM = X(I1+1) - X(I1)
      DELTAX = (XIN-X(I1))/DENOM
      DO 120 I = 1,NTERMS
      IX = I1 + I -1
      DELTA(I) = (X(IX)-X(I1))/DENOM
120   CONTINUE

C     ---------------------------------------------------------------
C     ACCUMULATE COEFFICIENTS A
C     ---------------------------------------------------------------
      A(1) = Y(I1)
      DO 140 K = 2,NTERMS
      PROD = 1.D0
      SUM = 0.D0
      IMAX = K-1
      IXMAX = I1 + K - 1
      DO 130 I = 1,IMAX
      J = K - I
      PROD = PROD*(DELTA(K) - DELTA(J))
      SUM = SUM - A(J)/PROD
130   CONTINUE
      A(K) = SUM + Y(IXMAX)/PROD
140   CONTINUE

C     ------------------------------------------------------------------
C     ACCUMULATE SUM OF EXPANSION
C     ------------------------------------------------------------------
      SUM = A(1)
      DO 160 J = 2,NTERMS
      PROD = 1.D0
      IMAX = J - 1
      DO 150 I = 1,IMAX
150   PROD = PROD*(DELTAX - DELTA(I))
      SUM = SUM + A(J)*PROD
160   CONTINUE

      TERP = SUM
      RETURN
	end
