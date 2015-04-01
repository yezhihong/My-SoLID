

      FUNCTION SIGBAR(E)

C     -----------------------------------------------------------------
C     COMPUTES THE ELASTIC CROSS SECTION (AS A FUNCTION OF E = INCIDENT
C     ELECTRON ENERGY) TIMES A FUNCTION FBAR(QSQ) DEFINED IN STEIN, EQ.
C     A44. THE CROSS-SECTION IS RETURNED IN UNITS OF
C             HBARC**2/(STR*MEV**2), 
C     SO THAT MULTIPLICATION BY HBARC**2 = (197.3E-13 MEV-CM)**2 GIVES
C     THE CROSS SECTION IN CM**2/STR.
C     -----------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ME,MT,MTF
      COMMON /ROSEMASS/ME,MT,MTF
      COMMON /ROSERUN/ES,EP,THETA
      COMMON /ROSEFBARCONSTANT/FBAR1,FBAR2,FBAR3,BT

C     -----------------------------------------------------------------
C     THE ELASTIC CROSS-SECTION COMPUTED HERE USES MO AND TSAIS FORM
C     FACTORS FJQ AND GJQ (SEE EQ. B3 IN MT).
C     THESE ARE RELATED TO W1 AND W2 OF STEIN (EQ. A3) BY_
C           FJQ = 4*W2
C           GJQ = 4*MT*MT*W1
C     (COMPARE B3 OF MT WITH A3 OF STEIN).
C     -----------------------------------------------------------------

      ALPHA = 1.D0/137.03604D0
      SINSQ = DSIN(THETA/2.D0)**2
      RECOIL = 1.D0/(1.D0 + 2.D0*E*SINSQ/MT)
      EPRM = E*RECOIL
      QSQ = 4.D0*E*EPRM*SINSQ
      CALL FMFAC(QSQ,FJQ,GJQ)

C     ----------------------------------------------------------------
C     XMOTT HERE IS REALLY THE MOTT CROSS-SECTION DIVIDED BY 4Z**2.
C     ----------------------------------------------------------------
      XMOTT = ALPHA*EPRM*DCOS(THETA/2.D0)/QSQ
      XMOTT = XMOTT*XMOTT
      XELAST = XMOTT*RECOIL*(FJQ + 2.D0*GJQ*(DTAN(THETA/2.D0)/MT)**2)
       
C   ----------------------------------------------------------------
C     COMPUTE FBAR (STEIN, EQ. A44).
C     ----------------------------------------------------------------
      FBAR = 1.D0 + .5772D0*BT + FBAR2*DLOG(QSQ/(ME*ME)) +
     .       FBAR3*DLOG(ES/EP)**2 + FBAR1


C     ----------------------------------------------------------------
C     COMPUTE SIGMA BAR_
C     ----------------------------------------------------------------
      SIGBAR = XELAST*FBAR
*       print *,'elastic:', eprm, qsq, xmott, fjq, gjq, xelast, sigbar

      RETURN
      END
