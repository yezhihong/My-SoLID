***********************************************************************
*     Xin Qian 2005.10.06  DXs calculation for SIDIS process for pion
*     and kaon  (Modified from Lingyan Zhu's code)
*     Based on Cteq6 pdf and bkk fragmentation function
**********************************************************************
*     targ: 0 for hydrogen and 1 for deuterium
*     ebeam : beam energy in MeV
*     e_pmin,e_pmax,e_thmin,e_thmax,
c e_phmin,e_phmax,p_pmin,p_pmax,p_thmin,p_thmax,p_phmin,
c p_phmax  : momentum and energy range for both arms
*     sigma(2): DXs for positive particle and negative particle
*     type: 0 for pion and 1 for kaon
***********************************************************************

      subroutine SIDIS(targ,type,ebeam,e_pmin,e_pmax,e_thmin,e_thmax,
     c e_phmin,e_phmax,p_pmin,p_pmax,p_thmin,p_thmax,p_phmin,
     c p_phmax,sigma)
      integer targ,type
      real e_pmin,e_pmax,e_thmin,e_thmax,e_phmin,e_phmax
      real p_pmin,p_pmax,p_thmin,p_thmax,p_phmin,p_phmax
      real sigma(2),ebeam
      
      real PI,MP,MN,MPI,MK,ME,MNUC,MH
      parameter (PI=3.14159265)

      real PB,KI_VEC(4),KF_VEC(4),EB,EF_H,EF_E,pf_vec(4)
      real PF_E,Etheta,phi_x,PF_H,ptheta
      
      REAL q2vec,q2,omega,q(3)
      
      real X,S,Y,z,theta_pq,phiq,pt

      real sqrtq2
      integer iset,k
      real dspi,dupi,duk,ddk,pdf(-5:5)
      double precision dspi1,dupi1,duk1,ddk1

      real dm2dp,dpipnl,dpimnl,dkpnl,dkmdl
      real pim_t,pip_t,pim_t_b,pip_t_b,bpt_m,bpt_p
      real rnorm1,rnorm2,rnorm_p,rnorm_m
      real xsection0_m,xsection0_p,xsection_m,xsection_p
      real domega_e,domega_p
      real km_t,km_t_b,kp_t,kp_t_b
      REAL rc_vec(4),p2_rc
      
      SIGMA(1)=0.
      SIGMA(2)=0.

      MP = 938.27231               !Proton mass
      MN = 939.56563               !Neutron mass
      MPI = 139.57018              !Pion mass
      MK = 493.677                 !Kaon mass
      ME = 0.51099906              !Electron mass

      PB=ebeam
      EB=sqrt(pb**2+me**2)
      
      MNUC=MP

      if (type.eq.0) then
         MH=MPI
      else if (type.eq.1) then
         MH=MK
      else
      endif
      
c     initial momentum of beam electron
      KI_VEC(1)=0.
      KI_VEC(2)=0.
      KI_VEC(3)=PB
      ki_vec(4)=eb

c     final electron mometum and angle
      PF_E=(e_pmin+e_pmax)/2.
      Etheta=(e_thmin+e_thmax)/2.*pi/180.
      ephi=0.
      EF_E=sqrt(PF_E**2+ME**2)
      

c     final electron momentum
      kf_vec(3)=pf_e*COS(etheta)
      kf_vec(1)=pf_e*SIN(etheta)*COS(ephi)
      kf_vec(2)=pf_e*SIN(etheta)*SIN(ephi)
      kf_vec(4)=ef_e

c     angle between  hadron and electron planes      
c      Phi_X=0.

c     hadron momentum
      PF_H=(p_pmin+p_pmax)/2.
      Ptheta=(p_thmin+p_thmax)/2.*pi/180.
      pphi=(p_phmin+p_phmax-e_phmin-e_phmax)/2.*pi/180.
      EF_H=sqrt(PF_H**2+MH**2)
      

      pf_vec(3)=pf_h*COS(ptheta)
      pf_vec(1)=pf_h*SIN(ptheta)*COS(pphi)
      pf_vec(2)=pf_h*SIN(ptheta)*SIN(pphi)
      pf_vec(4)=ef_h
      
      p2_rc=0.
      DO i=1,4
         rc_vec(i)=ki_vec(i)-kf_vec(i)-pf_vec(i)
         IF(i.LE.3) p2_rc=p2_rc+rc_vec(i)**2
      ENDDO
      rc_vec(4)=rc_vec(4)+mn
      IF(rc_vec(4)**2.LT.p2_rc+mn**2+1000.) GO TO 999

C
C ---------------------------------------------------------------------
C       Calculate Q, Omega and W.
C ---------------------------------------------------------------------
C
      OMEGA=EB-EF_E
      Q(1)=KI_VEC(1)-KF_VEC(1)
      Q(2)=KI_VEC(2)-KF_VEC(2)
      Q(3)=KI_VEC(3)-KF_VEC(3)
      
      q2vec=q(1)**2+q(2)**2+q(3)**2
      q2 = (q2vec-omega**2)                     !4-mom. transf.^2
      x = q2/(2.d0*mnuc*omega)
      
      theta_pq=acos((q(1)*pf_vec(1)+q(2)*pf_vec(2)+q(3)*pf_vec(3))
     c     /SQRT(q2vec)/pf_h)
      S=2.0*EB*MP+MP*MP
      Y=OMEGA/EB  
      
      

      z=EF_H/OMEGA
      PT=PF_H*SIN(THETA_PQ)
      
      
      

      iset=1
      call setctq6(iset)
      do k=-5,5               
         pdf(k) = real(Ctq6Pdf(k,dble(x),dble(sqrt(q2/1.D6))))
      enddo

      call bkk(dble(z),dble(q2/1.D6),dupi1,dspi1,duk1,ddk1)
      dupi=real(dupi1)
      dspi=real(dspi1)
      duk=real(duk1)
      ddk=real(ddk1)
      
      
      dm2dp=(1-z)**0.083583/(1+z)**1.9838   
      dpipnl=dupi*2.0/(1+dm2dp)
      dpimnl=dpipnl*dm2dp
      
      dkpnl=duk*2.0/(1+dm2dp)
      dkmnl=dkpnl*dm2dp
      
      if (type.eq.0) then
         if (targ.eq.0) then
            pim_t= (dpimnl*4*pdf(1)+dpipnl*pdf(2)
     >           +pdf(-1)*4*dpipnl+pdf(-2)*dpimnl
     >           +pdf(-3)*dspi+pdf(3)*dspi)/9.0
            pip_t= (dpipnl*4*pdf(1)+dpimnl*pdf(2)
     >           +4*pdf(-1)*dpimnl+pdf(-2)*dpipnl
     >           +pdf(-3)*dspi+pdf(3)*dspi)/9.0
         else if (targ.eq.1) then
            pim_t= (dpimnl*4*pdf(1)+dpipnl*pdf(2)
     >           +pdf(-1)*4*dpipnl+pdf(-2)*dpimnl
     >           +pdf(-3)*dspi+pdf(3)*dspi)/9.0
            pip_t= (dpipnl*4*pdf(1)+dpimnl*pdf(2)
     >           +4*pdf(-1)*dpimnl+pdf(-2)*dpipnl
     >           +pdf(-3)*dspi+pdf(3)*dspi)/9.0
            pim_t_b= (dpimnl*4*pdf(2)+dpipnl*pdf(1)
     >           +pdf(-2)*4*dpipnl+pdf(-1)*dpimnl
     >           +pdf(-3)*dspi+pdf(3)*dspi)/9.0
            pip_t_b= (dpipnl*4*pdf(2)+dpimnl*pdf(1)
     >           +4*pdf(-2)*dpimnl+pdf(-1)*dpipnl
     >           +pdf(-3)*dspi+pdf(3)*dspi)/9.0
            pim_t=pim_t+pim_t_b
            pip_t=pip_t+pip_t_b
         else
         endif
         
         
         bpt_m=4.694*1e-6       ! for pi-, 
         bpt_p=4.661*1e-6       ! for pi+ (HERMES GASKELL)
         rnorm1=2.0*PI/137.035**2*s*x*(1.0+(1.0-y)**2)/Q2**2
         rnorm2=rnorm1*EF_E/(Mnuc*omega*2*pi)*(197.3**2/100*1.0D9)
         rnorm_p=rnorm2*bpt_p/PI*exp(-bpt_p*pt**2)
         rnorm_m=rnorm2*bpt_m/PI*exp(-bpt_m*pt**2)
c     ds/domegadedzdpt^2dphi in the unit of b/GeV^3/sr
         xsection0_m=rnorm_m*pim_t*1.0d9/2.0
         xsection0_p=rnorm_p*pip_t*1.0d9/2.0
         
         xsection_m=xsection0_m*Pf_H**3*cos(theta_pq)/omega/EF_H*2.0
         xsection_p=xsection0_p*Pf_H**3*cos(theta_pq)/omega/EF_H*2.0      
         
         xsection_m=xsection_m*1E-3
         xsection_p=xsection_p*1E-3

         decay_pi=1.
         xsection_m=xsection_m*decay_pi
         xsection_p=xsection_p*decay_pi
      else if(type.eq.1) then
         if (targ.eq.0) then
            km_t=( dkmnl*4.*pdf(1)+ddk*pdf(2)
     >           +4.*pdf(-1)*dkpnl+pdf(-2)*ddk
     >           +pdf(-3)*dkmnl+pdf(3)*dkpnl)/9.0
            kp_t=( dkpnl*4.*pdf(1)+ddk*pdf(2)
     >           +4.*pdf(-1)*dkmnl+pdf(-2)*ddk
     >           +pdf(-3)*dkpnl+pdf(3)*dkmnl)/9.0
         else if (targ.eq.1) then
            km_t=( dkmnl*4.*pdf(1)+ddk*pdf(2)
     >           +4.*pdf(-1)*dkpnl+pdf(-2)*ddk
     >           +pdf(-3)*dkmnl+pdf(3)*dkpnl)/9.0
            kp_t=( dkpnl*4.*pdf(1)+ddk*pdf(2)
     >           +4.*pdf(-1)*dkmnl+pdf(-2)*ddk
     >           +pdf(-3)*dkpnl+pdf(3)*dkmnl)/9.0
            km_t_b= (dkmnl*4.*pdf(2)+ddk*pdf(1)
     >           +4.*pdf(-2)*dkpnl+pdf(-1)*ddk
     >           +pdf(-3)*dkmnl+pdf(3)*dkpnl)/9.0
            kp_t_b= (dkpnl*4.*pdf(2)+ddk*pdf(1)
     >           +4.*pdf(-2)*dkmnl+pdf(-1)*ddk
     >           +pdf(-3)*dkpnl+pdf(3)*dkmnl)/9.0
            km_t=km_t+km_t_b
            kp_t=kp_t+kp_t_b
         else
         endif
         bpt_m=4.694*1e-6       ! for pi-, 
         bpt_p=4.661*1e-6       ! for pi+ (HERMES GASKELL)
         rnorm1=2.0*PI/137.035**2*s*x*(1.0+(1.0-y)**2)/Q2**2
         rnorm2=rnorm1*EF_E/(Mnuc*omega*2*pi)*(197.3**2/100*1.0D9)
         rnorm_p=rnorm2*bpt_p/PI*exp(-bpt_p*pt**2)
         rnorm_m=rnorm2*bpt_m/PI*exp(-bpt_m*pt**2)
c     ds/domegadedzdpt^2dphi in the unit of b/GeV^3/sr
         xsection0_m=rnorm_m*km_t*1.0d9/2.0
         xsection0_p=rnorm_p*kp_t*1.0d9/2.0
         xsection_m=xsection0_m*Pf_H**3*cos(theta_pq)/omega/EF_H*2.0
         xsection_p=xsection0_p*Pf_H**3*cos(theta_pq)/omega/EF_H*2.0      
         
         xsection_m=xsection_m
     c        *1E-3
         xsection_p=xsection_p
     c        *1E-3
         
         decay_k=1.
         xsection_m=xsection_m*decay_k
         xsection_p=xsection_p*decay_k
      else
      endif
      
      sigma(1)=xsection_p
      sigma(2)=xsection_m
 999  RETURN
      end

