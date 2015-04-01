

C     ***************************************************************

      SUBROUTINE FMFAC(QSQ,FJQ,GJQ)

C     -----------------------------------------------------------------
C     COMPUTES THE FORM FACTORS FJQ AND GJQ NEEDED BY MO AND TSAI'S
C     INTERNAL BREMSSTRAHLUNG CALCULATION.
C
C     THE FORM FACTORS ARE BEST DEFINED BY EQUATION B.3 IN M & T, WHICH
C     SHOWS THAT THE CROSS-SECTION IN TERMS OF THEM IS:
C
C        (MOTT/4Z**2)(RECOIL FACTOR)(FJQ+(2/M**2)(TAN(THETA/2)**2)GJQ)
C
C     THE RELATIONSHIP BETWEEN THESE FORM FACTORS AND THE ELECTRIC AND
C     MAGNETIC FORM FACTORS GE AND GM (IN MT) IS (EQ. III.2, III.3)_
C
C         FJQ = 4(GE**2 + TAU*GM**2)/(1+TAU)
C         GJQ = -(Q**2)(GM**2)
C     WHERE
C         TAU = -(Q**2)/(4MT**2)
C         MT = MASS OF TARGET (INITIAL MASS)
C         Q = 4-MOMENTUM (MT HAS Q**2 NEGATIVE, Q**2 = E**2 - QVEC**2).
C
C     THE GE AND GM USED THIS WAY IN MT ARE EXACTLY THE SAME AS THE
C     NORMAL ELECTRIC AND MAGNETIC FORM FACTORS ( FRAUENFELDER AND
C     HENLEY EQ.6.22 AND 6.43).  EXCEPT, GE AND GM ARE NORMALIZED TO Z
C     (i.e., GE(0)=GM(0)=Z).
C     -----------------------------------------------------------------
C     THIS SUBROUTINE COMPUTES FORM FACTORS ONE OF FIVE WAYS: BY USING
C     INPUT FORM FACTORS (AND INTERPOLATING), OR BY USING APPROXIMATE
C     ANALYTIC EXPRESSIONS FOR PROTONS, HE3, HE4, OR DEUTERIUM. THE
C     DESIRED FORM FACTOR IS CHOSEN BY SPECIFYING ITARGT (ITARGT = -1 
C     TO 5 CORRESPONDS TO THE ABOVE 5 CHOICES IN ORDER).
C
C     THE CONVENTIONS FOR THE CHOICES ARE:
C     1. CARDS  : THE ROUTINE ASSUMES THAT ARRAYS YH AND YH2 CONTAIN
C                 FJQ AND GJQ.
C     2. PROTON : THE ROUTINE COMPUTES GE AND GM USING EQ. III.4 OF MT.
C                 (Z=1, SO THERE IS NO REASON TO THINK ABOUT NORMALIZA-
C                 TION).
C     3. HE3    : THE ROUTINE COMPUTES THE ELECTRIC AND MAGNETIC FORM
C                 FACTORS FCH AND FMAG USING THE ANALYTICAL FORM IN Mc-
C                 CARTHY, SICK AND WHITNEY (MSW), PHYS. REV. C, VOL 15,
C                 NO. 48 PP. 1396-1414.  FCH AND FMAG ARE EXACTLY THE
C                 SAME AS GE AND GM, EXCEPT FOR NORMALIZATION:
C                       GE = Z*FCH
C                       GM = Z*(1+K)*FMAG
C                 WHERE K IS THE ANOMALOUS MAGNETIC MOMENT OF HE3 RE-
C                 FERRED TO HE3
C                 (i.e., MAG. MOM. = (1+K)*Z*E*HBAR/(2*MHE3*C)).
C                 K = -4.20 FOR HE3.
C     4. HE4    : EXACTLY AS HE3, EXCEPT FMAG = 0.
C     5. D2     : THE ROUTINE COMPUTES THE FORM FACTORS W1 AND W2 USING
C                 THE ANALYTICAL EXPRESSION IN STEIN, EQS. A9-A12.  W1 
C                 AND W2 ARE RELATED TO FJQ AND GJQ BY_
C                       FJQ = 4*W2
C                       GJQ = 4*MT*MT*W1
C                 (COMPARE A3 OF STEIN WITH B3 OF MT.)
C     -----------------------------------------------------------------
C     THIS ROUTINE CAN BE CALLED WITH QSQ = Q**2 EITHER POSITIVE OR
C     NEGATIVE. THE FORM USED IN EACH OF THE FOLLOWING SECTIONS IS THAT
C     USED IN THE REFERENCE PAPER.
C     -----------------------------------------------------------------

      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 ME,MT,MTF
      REAL*8 K,KN,KP,MP
      COMMON /ROSEMASS/ME,MT,MTF
      COMMON /ROSEFMFACCONSTANT/ITARGT
      common /ROSEklgcom/ig,rmod,imod,nklg
      DIMENSION H(6)

C     -----------------------------------------------------------------
C     THE CONSTANTS BELOW ARE USED TO COMPUTE THE HE3 AND HE4 FORM
C     FACTORS K IS THE ANOMALOUS MAGNETIC MOMENT OF HE3.
C     -----------------------------------------------------------------
      DATA K/-4.2D0/

C     -----------------------------------------------------------------
C     SUBSCRIPTS 1, 2, AND 3 REFER, RESPECTIVELY, TO HE3 ELECTRIC, HE3
C     MAGNETIC, AND HE4 ELECTRIC.
C     -----------------------------------------------------------------
      DATA A1,B1,C1,D1,P1,Q01/.675D0,.366D0,.836D0,-6.78D-3,
     .                        .90D0,3.98D0/
      DATA A2,B2,C2/.654D0,.456D0,.821D0/
      DATA A3,B3/.316D0,.675D0/

C     ----------------------------------------------------------------
C     THE FOLLOWING CONSTANTS ARE USED IN THE DEUTERIUM FORM FACTOR
C     COMPUTATION KN AND KP ARE THE ANOMALOUS MAGNETIC MOMENTS OF THE
C     NEUTRON AND PROTON MP IS THE MASS OF THE PROTON (IN GEV).
C     THE CONSTANTS IN ARRAY H ARE USED TO COMPUTE P(Q**2) ON THE WAY
C     TO FIND THE DEUTERIUM FORM FACTORS.
C     ----------------------------------------------------------------
      DATA KN,KP/-1.91348D0,1.7927D0/
      DATA MP/.938256D0/
      DATA H/1.0007D0,1.01807D0,1.05584D0,0.836380D0,0.6864584D0,
     .       0.672830D0/

C     -----------------------------------------------------------------

C      print *,'fmfac begin'

      IF(QSQ.le.0) THEN
        QSQP = DABS(-QSQ)
      ELSE
        QSQP = DABS(QSQ)
      ENDIF
      QSQN = -QSQP

      GO TO (10,20,30,40,80), ITARGT

C     -----------------------------------------------------------------
C     FORM FACTORS WERE READ IN FROM CARDS.  THEY HAVE BEEN TURNED INTO
C     FJQ AND GJQ ALREADY.
C     -----------------------------------------------------------------
10    FJQ = TERP(QSQP,1)
      GJQ = TERP(QSQP,2)
      RETURN

C     -----------------------------------------------------------------
C     PROTON FORM FACTORS FROM MT, EQ. III.4.
C     -----------------------------------------------------------------
20    CONTINUE
      GE = (1.D0 - (QSQN/.71D6))**(-2)
cccklg      
      GM = 2.793D0*GE
ccc
cccccccccc   use donal form factor routine
	qsqd=qsqp*2.568162e-05  ! (1/197.38)**2  converts Mev to fm-2
        
 	call nform(ig,rmod,imod,QSQd,gepd,gend,gmpd,gmnd)

c        write(23,676) QSQd, GE, GM, gepd,gend,gmpd,gmnd
676     FORMAT(1X,7(1X,1PD12.5,1X))

	if(nklg.eq.1)then
	write(6,'(''  day form factor sub '')')
	write(6,'(''   qsq in fermi and in energy='',2e20.5)')qsqd,qsqp
	write(6,'(''   gepdon, dipole='',2e20.5)')gepd,ge
	endif

	ge=gepd
	gm=gmpd
      GO TO 60

C     ------------------------------------------------------------------
C     HE3 FORM FACTORS FROM MSW, P. 1403.
C     NOTICE THAT THIS SECTION REQUIRES Q**2 TO BE POSITIVE AND IN FERMI
C     **(-2) INSTEAD OF MEV**2
C     ------------------------------------------------------------------
30    CONTINUE
      Q = DSQRT(QSQP)/197.32D0
      QSQP = Q*Q

      F0 = DEXP(-A1*A1*QSQP) - B1*B1*QSQP*DEXP(-C1*C1*QSQP)
      DF = D1*DEXP(-((Q-Q01)/P1)**2)
      FCH = F0 + DF

      FMAG = DEXP(-A2*A2*QSQP) - B2*B2*QSQP*DEXP(-C2*C2*QSQP)
      GO TO 50

C     ------------------------------------------------------------------
C     HE4 FORM FACTORS FROM MSW, P. 1403.
C     AS ABOVE, REQUIRE Q**2 POSITIVE, IN FERMI**(-2).
C     ------------------------------------------------------------------
40    CONTINUE    
      QSQP = QSQP/(197.32D0**2)
      FCH = (1.D0 - (A3*A3*QSQP)**6)*DEXP(-B3*B3*QSQP)
      FMAG = 0.D0

C     ------------------------------------------------------------------
C     FOR HE3 AND HE4, NORMALIZE BY Z
C     ------------------------------------------------------------------
50    CONTINUE
      Z = 2.D0
      GE = Z*FCH
      GM = Z*(1.D0+K)*FMAG

C     ------------------------------------------------------------------
C     FOR PROTON, HE3 AND HE4, SQUARE THE FORM FACTORS
C     ------------------------------------------------------------------
60    CONTINUE
      GESQR = GE*GE
      GMSQR = GM*GM

C     ------------------------------------------------------------------
C     COMPUTE FJQ AND GJQ FROM MT (EQ. III.2 AND III.3, MT).
C     ------------------------------------------------------------------
70    TAU = -QSQN/(4.D0*MT*MT)
      FJQ = 4.D0*(GESQR + TAU*GMSQR)/(1.D0+TAU)
      GJQ = -QSQN*GMSQR

      

      RETURN

C     ------------------------------------------------------------------
C     DEUTERIUM FORM FACTORS FROM STEIN, EQ. A10.
C     INITIALLY, Q MUST BE IN GEV.
C     ------------------------------------------------------------------
80    CONTINUE


      GE = (1.D0 - (QSQN/.71D6))**(-2)
cccklg      
      GM = 2.793D0*GE
ccc
ccccccccccuse donal form factor routine
      qsqd=qsqp*2.568162e-05    ! (1/197.38)**2  converts Mev to fm-2
        
      call nform(ig,rmod,imod,QSQd,gepd,gend,gmpd,gmnd)

c        write(23,676) QSQd, GE, GM, gepd,gend,gmpd,gmnd
c676     FORMAT(1X,7(1X,1PD12.5,1X))

      if(nklg.eq.1)then
         write(6,'(''  day form factor sub '')')
         write(6,'(''   qsq in fermi and in energy='',2e20.5)')qsqd,qsqp
         write(6,'(''   gepdon, dipole='',2e20.5)')gepd,ge
      endif
      
      GESQR = (gepd**2+gend**2)/2.
      GMSQR = (gmpd**2+gmnd**2)/2.
      TAU = -QSQN/(4.D0*MT*MT)*4.0
      FJQ = 4.D0*(GESQR + TAU*GMSQR)/(1.D0+TAU)
      GJQ = -QSQN*GMSQR
      

      RETURN




      Q = DSQRT(QSQP)*.001D0
      TAU = Q*Q/(4.D0*MP*MP)

C     ------------------------------------------------------------------
C     COMPUTE THE PROTON TERMS GEP AND GMP (EQ. A7).
C     FIRST FIND P(Q**2) (EQ. A8).
C     ------------------------------------------------------------------
      P = 0.D0
      DO 100 II = 1,6
      I = II - 1
      PROD = 1.
      DO 90 JJ = 1,6
      J = JJ - 1
      IF(J .EQ. I) GO TO 90
      PROD = PROD*(Q-J)/(I-J)
90    CONTINUE
      P = P + H(II)*PROD
100   CONTINUE

      GEP = P/((1.D0+Q*Q/.71D0)**2)
      GMP = (1.D0+KP)*GEP

C     ------------------------------------------------------------------
C     COMPUTE TERMS LEADING TO W1 AND W2 (EQ. A12)
C     ------------------------------------------------------------------
      GMN = KN*GEP
      F1N = TAU*GMN/(1.D0+TAU)
      F1P = (GEP+TAU*GMP)/(1.D0+TAU)
      F2N = GMN/(KN*(1.D0+TAU))
      F2P = (GMP-GEP)/(KP*(1.D0+TAU))
      GP = F1N + F1P
      GS = GP + KN*F2N + KP*F2P

C     ------------------------------------------------------------------
C     COMPUTE FD (EQ. A11).  THIS REQUIRES Q IN INVERSE FERMIS.
C     ------------------------------------------------------------------
      Q = DSQRT(QSQP)/197.3D0
      IF(Q .GT. 1.D-10)
     .FD = (1.580D0/Q)*(DATAN(Q/.930D0) - 2.D0*DATAN(Q/3.19D0)
     .      + DATAN(Q/5.45D0))
      IF(Q .LT. 1.D-10) FD=1.580D0*(1.D0/.930D0-2.D0/3.19D0+1.D0/5.45D0)

C     ------------------------------------------------------------------
C     COMPUTE W1 AND W2 (EQ. A10) AND FJQ GJQ
C     ------------------------------------------------------------------
      W1 = FD*FD*(2.D0/3.D0)*TAU*GS*GS
      W2 = FD*FD*GP*GP + W1
      FJQ = FJQ + 4.D0*W2
      GJQ = GJQ + 4.D0*MT*MT*W1
c      print*,FJQ,GJQ
      RETURN

      
      END

C     *****************************************************************
