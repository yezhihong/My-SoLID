      SUBROUTINE FBAR(THETA)

C     ------------------------------------------------------------------
C     COMPUTES TERMS USED TO CALCULATE THE FUNCTION FBAR (SEE STEIN,
C     EQ. A43).  THIS SUBROUTINE MUST BE CALLED EACH TIME THETA CHANGES.
C     ------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,M-Z)
      COMMON /ROSEFBARCONSTANT/FBAR1,FBAR2,FBAR3,BTTOT

      ALPHA = 1.D0/137.03604D0
      PI = 3.1415926536D0
      AP = ALPHA/PI

      SINSQ = DSIN(THETA/2.D0)**2
      COSSQ = DCOS(THETA/2.D0)**2
      PHI = SPENCE(COSSQ,SINSQ)

      FBAR1 = 2.D0*AP*(-14.D0/9.D0) + AP*(PI*PI/6.D0 - PHI)
      FBAR2 = 2.D0*AP*(13.D0/12.D0)
      FBAR3 = - .5D0*AP

      RETURN
      END
