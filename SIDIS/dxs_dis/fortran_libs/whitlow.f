********************************************************************8
*     Xin Qian 2005. 10. 06 modified from whilow code
*********************************************************************
*     target: target type 0 for hydrogen  1 for deuterium
*     ebeam1: beam energy in MeV
*     pmin1,pmax1: lower (higher) cut of momentum range (MeV)
*     thmin(thmax): lower (higher) edge of theta range (degree)
*     pdel1(thdel1): step size of momentum (angle)
*     sin_em_xs: DXs in nb
*     success: success flag
*********************************************************************
      

      subroutine whitlow(target,ebeam1,pmin1,pdel1,pmax1,
     cthmin,thdel,thmax,sin_em_xs,success)
      integer target ! 0 for hydrogen, 1 for deterium, 2 for He3
      real ebeam1, pmin1, pmax1, pdel1, thmin, thdel, thmax
      real ebeam, pmin, pmax, pdel ! beam energy mometum range and step
      real protonmass,alpha, PI  ! some constatnts
      parameter (protonmass=0.938)
      parameter (Alpha=137.03604)
      parameter (PI=3.14159265)
      integer thtotal,ptotal ! integer of the total steps
      integer i,j,callstats,itotalpoint
      real p,theta,sinsq,qmu2,nu,x,W,xsection ! some kinetic variables
      real result(3),totalXSection,sin_em_xs ! DXs and Total Xs
      real tmp(3)

      itotalpoint=0
      totalXSection=0.

      ebeam=ebeam1/1000.
      pmin=pmin1/1000.
      pmax=pmax1/1000.
      pdel=pdel1/1000.
      
      if (thdel.eq.0.) then
         thtotal=1
      else
         thtotal=(thmax-thmin)/thdel+1
      endif
      
      if (pdel.eq.0.) then
         ptotal=1
      else
         ptotal=(pmax-pmin)/pdel+1
      endif
        
      do i=0,thtotal-1
         do j=0,ptotal-1
            p=pmin+j*pdel
            theta=(thmin+i*thdel)*PI/180.;
            sinsq=sin(theta/2.)**2
            qmu2=4.*ebeam*p*sinsq
            nu=ebeam-p
            x=qmu2/(2.*protonmass*nu)
            W=protonmass**2+2.*protonmass*nu-qmu2
            if (target.lt.1.5) then
               call XSECFT(ebeam,p,sinsq,target,result,success)
            endif
            if (target.eq.2) then
               call XSECFT(ebeam,p,sinsq,0,result,success)
               call XSECFT(ebeam,p,sinsq,1,tmp,success)
               result(3)=result(3)+tmp(3)
            endif
            xsection=result(3)
            xsection=xsection
            p=p*1000.
            theta=theta/PI*180.
            itotalpoint=itotalpoint+1
            totalXSection=totalXSection+xsection
         enddo
      enddo
      
      sin_em_xs=totalXSection/itotalpoint
      
      end
