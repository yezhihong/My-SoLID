      real function output(part,file)
      include ?
      real x,z,pt
      real p_e,p_h,phi_e,phi_h,theta_e,theta_h,weight
      integer part,file

      x=qgetkin(9,0.,part)
      z=qgetkin(4,0.,part)
      pt=qgetkin(6,0.,part)
      
      p_e = getkinp(0,2)
      p_h = getkinp(0,3)
      phi_e = getkinp(4,2)
      phi_h = getkinp(4,3)
      theta_e = getkinp(3,2)
      theta_h = getkinp(3,3)
      
c$$$      open(UNIT=21,FILE='temp.dat',access='DIRECT',
c$$$     c FORM='FORMATTED',RECL=4096,STATUS='NEW')
      weight = wgcros

      write(file,*) x,z,pt,p_e,p_h,phi_e,phi_h,theta_e,theta_h
c$$$      q2=qgetkin(2,0.,part)
c$$$      w=qgetkin(3,0.,part)
c$$$      wp=qgetkin(5,0.,part)
c$$$      theta1=getkinp(3,2,part)
c$$$      theta2=getkinp(3,3,part)
c$$$      hitnum=hitraxy(-8,1,8,3)
c$$$      hitnum1=hitraxy(-8,1,8,2)
c$$$
c$$$      if ((w.ge.2.3).and.(wp.ge.1.6).and.(abs(z-0.5).le.0.2)
c$$$     c     .and.(q2.ge.1.0).and.(x.le.1.).and.(theta1.le.17.)
c$$$     c     .and.(theta2.le.17.).and.(hitnum.ge.1.0).and.
c$$$     c     (hitnum1.ge.1.0)) then
c$$$         call hfnt(9500,wgcros)
c$$$      endif
c$$$      
c$$$      skimed=1.
      end
      
      include 'qgetkin.f'
      include 'hitraxy.f'
