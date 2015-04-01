********************************************************************8
*     Xin Qian 2005.10.06
*     Modified from wiser
********************************************************************
*     part: particle type
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p
*     target: target type 0 for hydrogen  1 for deuterium
*     ebeam: beam energy in MeV
*     pmin,pmax: lower (higher) cut of momentum range (MeV)
*     thmin(thmax): lower (higher) edge of theta range (degree)
*     pdel1(thdel1): step size of momentum (angle)
*     sin_h_xs: DXs in nb
*     success: success flag
********************************************************************
      subroutine wiser(part,targ,ebeam,pmin,pdel,pmax,
     c thmin,thdel,thmax,sin_h_xs,success)
      IMPLICIT NONE       
      DOUBLE PRECISION E0MM,PMM,THETA_DEG,RAD_LENgth,SIGMA
      integer part,targ,success
      real ebeam,pmin,pdel,pmax,thmin,thdel,thmax,sin_h_xs
      real rad_len,xsection,p,theta,TotalXSection
      integer ptotal,thtotal,itotalpoint
      integer i,j,num
      
      success=0
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar     
      if (part.eq.1) then
         num=2
      else if (part.eq.2) then
         num=3
      else if (part.eq.3) then
         num=4
      else if (part.eq.4) then
         num=5
      else if (part.eq.5) then
         num=6
      else 
      endif
      
      ptotal=(pmax-pmin)/pdel+1
      thtotal=(thmax-thmin)/thdel+1
      itotalpoint=0
      TotalXsection=0.
      
      if (targ.eq.0) then
         ! Hydrogen target
         rad_len=2.7 + 0.5*(15.*0.0708)/61.28*100.
      else if (targ.eq.1) then
         ! Deuterium target
         rad_len=2.7 + 0.5*(12.*0.169)/122.4*100.
         rad_len=rad_len*2.
      else
      endif
      
      do i=0,ptotal
         do j=0,thtotal
            p=pmin+i*pdel
            theta=thmin+j*thdel
            
            E0MM=dble(ebeam)
            PMM=dble(p)
            theta_deg=dble(theta)
            rad_length=dble(rad_len)
            
            Call WISER_ALL_SIG(E0MM,PMM,theta_deg,RAD_LENgth,part,SIGMA)
            
            xsection=real(SIGMA)
            
            xsection=xsection
            itotalpoint=itotalpoint+1
            TotalXSection=TotalXSection+xsection
         enddo
      enddo
      
      sin_h_xs=TotalXSection/itotalpoint
      
      end






C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
C Xiaochao Zheng, this subroutine added July 07, 2000

      DOUBLE PRECISION FUNCTION QUADMO(FUNCT,PLOWER,PUPPER,EPSLON,NLVL)
      REAL*8 FUNCT,PLOWER,PUPPER,EPSLON                                  
      INTEGER NLVL                                                  
      INTEGER   LEVEL,MINLVL/3/,MAXLVL/24/,IRETRN(50),I                 
      REAL*8 VALINT(50,2), VMX(50), RX(50), FMX(50), FRX(50),           
     1   FMRX(50), ESTRX(50), EPSX(50)                                 
      REAL*8  R, FL, FML, FM, FMR, FR, EST, ESTL, ESTR, ESTINT,VL,      
     1   AREA, ABAREA,   VM, COEF, ROMBRG,   EPS                        

!      IMPLICIT REAL*8 (A-H,O-Z)
         LEVEL = 0                                                     
         NLVL = 0                                                      
         ABAREA = 0.0                                                  
         VL = PLOWER                                                     
         R = PUPPER                                                     
         FL = FUNCT(VL)                                                 
         FM = FUNCT(0.5*(VL+R))                                         
         FR = FUNCT(R)                                                 
         EST = 0.0                                                     
         EPS = EPSLON                                                  
  100 LEVEL = LEVEL+1                                                  
      VM = 0.5*(VL+R)                                                    
      COEF = R-VL                                                       
      IF(COEF.NE.0) GO TO 150                                          
         ROMBRG = EST                                                  
         GO TO 300                                                     
  150 FML = FUNCT(0.5*(VL+VM))                                           
      FMR = FUNCT(0.5*(VM+R))                                           
      ESTL = (FL+4.0*FML+FM)*COEF                                      
      ESTR = (FM+4.0*FMR+FR)*COEF                                      
      ESTINT = ESTL+ESTR                                               
      AREA=DABS(ESTL)+DABS(ESTR)                                       
      ABAREA=AREA+ABAREA-DABS(EST)                                     
      IF(LEVEL.NE.MAXLVL) GO TO 200                                    
         NLVL = NLVL+1                                                 
         ROMBRG = ESTINT                                               
         GO TO 300                                                     
 200  IF((DABS(EST-ESTINT).GT.(EPS*ABAREA)).OR.                        
     1         (LEVEL.LT.MINLVL))  GO TO 400                           
         ROMBRG = (1.6D1*ESTINT-EST)/15.0D0                            
  300    LEVEL = LEVEL-1                                               
         I = IRETRN(LEVEL)                                              
         VALINT(LEVEL, I) = ROMBRG                                     
         GO TO (500, 600), I                                           
  400    IRETRN(LEVEL) = 1                                              
         VMX(LEVEL) = VM                                                 
         RX(LEVEL) = R                                                 
         FMX(LEVEL) = FM                                               
         FMRX(LEVEL) = FMR                                             
         FRX(LEVEL) = FR                                               
         ESTRX(LEVEL) = ESTR                                           
         EPSX(LEVEL) = EPS                                             
         EPS = EPS/1.4                                                 
         R = VM                                                         
         FR = FM                                                       
         FM = FML                                                      
         EST = ESTL                                                    
         GO TO 100                                                     
  500    IRETRN(LEVEL) = 2                                              
         VL = VMX(LEVEL)                                                 
         R = RX(LEVEL)                                                 
         FL = FMX(LEVEL)                                               
         FM = FMRX(LEVEL)                                              
         FR = FRX(LEVEL)                                               
         EST = ESTRX(LEVEL)                                            
         EPS = EPSX(LEVEL)                                             
         GO TO 100                                                     
  600 ROMBRG = VALINT(LEVEL,1)+VALINT(LEVEL,2)                         
      IF(LEVEL.GT.1) GO TO 300                                         
      QUADMO = ROMBRG /12.0D0                                          
      RETURN                                                           
      END                                                              

      Subroutine WISER_ALL_SIG(E0MM,PMM,THETA_DEG,RAD_LEN,TYPE,SIGMA)

!------------------------------------------------------------------------------
! Calculate pi,K,p  cross section for electron beam on a proton target
! IntegrateQs over function WISER_FIT using integration routine QUADMO
! E0         is electron beam energy, OR max of Brem spectra
! P,E       is scattered particle  momentum,energy
! THETA_DEG  is kaon angle in degrees
! RAD_LEN (%)is the radiation length of target, including internal
!                (typically 5%)
!               = .5 *(target radiation length in %) +5.
!       ***  =100. IF BREMSTRULUNG PHOTON BEAM OF 1 EQUIVIVENT QUANTA
***
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! SIGMA      is output cross section in nanobars/GeV-str
!------------------------------------------------------------------------------

      IMPLICIT NONE       
      DOUBLE PRECISION E0MM,PMM,E0,P,THETA_DEG,RAD_LEN,SIGMA
      INTEGER TYPE
      COMMON/WISER_ALL/ E,P_COM,COST,P_T,TYPE_COM,PARTICLE ,M_X,U_MAN
      DOUBLE PRECISION E,P_COM,COST,P_T,M_X,U_MAN
      INTEGER TYPE_COM,PARTICLE
!  Wiser's fit    pi+     pi-    k+     k-     p+      p-   
      REAL*8 A5(6)/-5.49,  -5.23, -5.91, -4.45, -6.77,  -6.53/
      REAL*8 A6(6)/-1.73,  -1.82, -1.74, -3.23,  1.90,  -2.45/
      REAL*8 MASS2(3)/.019488, .2437, .8804/
      REAL*8 MASS(3)/.1396, .4973, .9383/ 
      REAL*8 MP/.9383/,  MP2/.8804/, RADDEG/.0174533/
      REAL*8  M_L,SIG_E
      REAL*8 E_GAMMA_MIN,WISER_ALL_FIT,QUADMO,E08,EPSILON/.003/
      EXTERNAL WISER_ALL_FIT                        
      INTEGER N,CHARGE
                                            
      P=PMM/1000.
      E0=E0MM/1000.

      P_COM = P
      TYPE_COM = TYPE
      PARTICLE = (TYPE+1)/2       ! 1= pi, 2= K, 3 =P
      CHARGE = TYPE -2*PARTICLE +2  ! 1 for + charge, 2 for - charge
      E08 =E0
                
      E =SQRT(MASS2(PARTICLE) + P**2)

      COST = COS(RADDEG * THETA_DEG)
      P_T = P * SIN(RADDEG * THETA_DEG)
      IF(TYPE.LE.4) THEN  !mesons
       IF(CHARGE.EQ.1) THEN   ! K+ n final state
        M_X = MP
       ELSE   ! K- K+ P final state
        M_X = MP+ MASS(PARTICLE)
       ENDIF
      ELSE  ! baryons 
       IF(CHARGE.EQ.1) THEN   ! pi p  final state
        M_X = MASS(1)  ! pion mass
       ELSE   ! P P-bar  P final state
        M_X = 2.*MP
       ENDIF
      ENDIF
      E_GAMMA_MIN = (M_X**2 -MASS2(PARTICLE ) -MP2+2.*MP*E)/
     >  (2.*(MP -E +P*COST))
!      WRITE(10,'(''E_GAMMA_MIN='',F10.2,''  p_t='',F8.2)')
!     >     E_GAMMA_MIN,P_T
!      E_GAMMA_MIN = MP *(E + MASS(PARTILCE))/(MP -P*(1.-COST))
      
*      print *,E_GAMMA_MIN

      IF(E_GAMMA_MIN.GT..1) THEN !Kinematically allowed?
       M_L = SQRT(P_T**2 + MASS2(PARTICLE))    

       IF(TYPE.NE.5) THEN  ! everything but proton
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) *EXP(A6(TYPE) *P_T**2/E)
       ELSE ! proton

        U_MAN = ABS(MP2 + MASS2(PARTICLE) -2.*MP*E)
        SIG_E = QUADMO(WISER_ALL_FIT,E_GAMMA_MIN,E08,EPSILON,N)  *
     >           EXP(A5(TYPE) *M_L) 
       ENDIF
       SIGMA = P**2/E * 1000. * RAD_LEN/100. *SIG_E 
      ELSE   ! Kinematically forbidden
       SIGMA = 0.
      ENDIF

      RETURN
      END


      REAL*8 FUNCTION WISER_ALL_FIT(E_GAMMA)

!---------------------------------------------------------
! Calculates  pi, k, p  cross section for gamma + p -> k
!  It is already divided by E_GAMMA, the bremstrulung spectra
! David Wiser's fit from Thesis, eq. IV-A-2 and Table III.
! Can be called from WISER_SIG using integration routine QUADMO
! E,P are KAON energy and momentum
! P_t is KAON transverse momentum
! P_CM is KAON center of mass momentum
! P_CM_L is KAON center of mass longitudinal momentum
! TYPE:     1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
! E_GAMMA is photon energy.
!             Steve Rock 2/21/96
!---------------------------------------------------------
                           
      IMPLICIT NONE       
      COMMON/WISER_ALL/ E,P,COST,P_T,TYPE,PARTICLE,M_X,U_MAN

      REAL*8  E,P,COST,P_T,M_X,U_MAN
      INTEGER  TYPE  !  1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p, 6=p-bar
      INTEGER PARTICLE   ! 1= pi, 2= K, 3 =P
!  Wiser's fit    pi+     pi-    k+     k-     p+       p- 
      REAL*8 A1(6)/566.,  486.,   368., 18.2,  1.33E5,  1.63E3 / 
      REAL*8 A2(6)/829.,  115.,   1.91, 307.,  5.69E4, -4.30E3 / 
      REAL*8 A3(6)/1.79,  1.77,   1.91, 0.98,  1.41,    1.79 / 
      REAL*8 A4(6)/2.10,  2.18,   1.15, 1.83,   .72,    2.24 /
      REAL*8 A6/1.90/,A7/-.0117/ !proton only
      REAL*8 MASS2(3)/.019488, .2437, .8804/
      REAL*8 MASS(3)/.1396, .4973, .9383/ 
      REAL*8 MP2/.8804/,MP/.9383/, RADDEG/.0174533/
      REAL*8 X_R,S,B_CM, GAM_CM,  P_CM
      REAL*8 P_CM_MAX, P_CM_L
      REAL*8 E_GAMMA
                                            

!Mandlestam variables                                                
      S = MP2 + 2.* E_GAMMA * MP    

!Go to Center of Mass to get X_R
      B_CM = E_GAMMA/(E_GAMMA+MP)
      GAM_CM = 1./SQRT(1.-B_CM**2)
      P_CM_L = -GAM_CM *B_CM *E + 
     >          GAM_CM * P * COST
      P_CM = SQRT(P_CM_L**2 + P_T**2)  


      P_CM_MAX =SQRT (S +(M_X**2-MASS2(PARTICLE))**2/S 
     >    -2.*(M_X**2 +MASS2(PARTICLE)) )/2.
      X_R =  P_CM/P_CM_MAX   
       IF(X_R.GT.1.) THEN  ! Out of kinematic range
        WISER_ALL_FIT = 0.
       ELSEIF(TYPE.NE.5) THEN  ! not the proton
        WISER_ALL_FIT = (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)/E_GAMMA  
       ELSE ! special formula for proton
        WISER_ALL_FIT = ( (A1(TYPE) + A2(TYPE)/SQRT(S)) *
     >   (1. -X_R + A3(TYPE)**2/S)**A4(TYPE)          /
     >   (1.+U_MAN)**(A6+A7*S) )/E_GAMMA  
       ENDIF
      
      RETURN
      END


