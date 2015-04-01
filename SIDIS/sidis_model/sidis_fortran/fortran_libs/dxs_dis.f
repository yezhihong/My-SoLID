      SUBROUTINE DXS_DIS(targ,part,beam_energy,
     c           p_e,theta_e,phi_e,p_h,theta_h,phi_h,
     c           sdxs_e,sdxs_h,cdxs)
!     targ: 0->proton,1->deuterium,2->he3      
!     part: particle type
!           1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p
!     beam_energy,p_e,p_h(GeV),theta_e,phi_e,theta_h,phi_h (DEG) 
      
      real theta_e,theta_h,p_e,p_h,phi_e,phi_h
      real beam_energy,factor,pe_temp,ph_temp
      real sdxs_e,sdxs_h,cdxs,temp
      real temp1(2),temp2(2)
      integer part,success
             
         pe_temp = p_e*1000.0
         ph_temp = p_h*1000.0
         
         call whitlow(2,beam_energy,pe_temp,0.0,pe_temp
     c        ,theta_e,0.0,theta_e,sdxs_e,success)
         
         if (success.ne.0.or.sdxs_e.le.0.or.sdxs_e.ge.10000.) then
            call qfsrad(beam_energy/1000.,p_e,theta_e/180.*3.1415926
     c           ,0,sdxs_e,success)
            call qfsrad(beam_energy/1000.,p_e,theta_e/180.*3.1415926
     c           ,1,temp,success)
            sdxs_e=sdxs_e + temp
            if (sdxs_e.ge.0.and.sdxs_e.le.1000000) then
            else
               sdxs_e = 0
            endif
         endif
         
         call wiser(part,0,beam_energy,ph_temp,0.1,ph_temp,
     c        theta_h,0.1,theta_h,sdxs_h,success)
         call wiser(part,1,beam_energy,ph_temp,0.1,ph_temp,
     c        theta_h,0.1,theta_h,temp,success)
         sdxs_h = sdxs_h+temp
         
         factor = (2.33369*exp(-0.508963*p_h*
     c        sin(theta_h/180.*3.1415926)*sqrt(0.938*0.938
     c        +2.*0.938*5.892)/sqrt(0.938*0.938+2.*0.938
     c        *beam_energy/1000.)))
         
         if (factor.le.1) factor = 1
         
         sdxs_h = sdxs_h /factor
         if (sdxs_h.le.0.or.success.ne.0) then
            sdxs_h = 0
         endif

         ! proton
         CALL SIDIS(0,0,beam_energy,pe_temp,pe_temp,theta_e,
     c        theta_e,phi_e,phi_e,ph_temp,ph_temp,theta_h
     c        ,theta_h,phi_h,phi_h,temp1)
        ! deutron
         CALL SIDIS(1,0,beam_energy,pe_temp,pe_temp,theta_e,
     c        theta_e,phi_e,phi_e,ph_temp,ph_temp,theta_h
     c        ,theta_h,phi_h,phi_h,temp2)
                     
         if (part.eq.1) then
            cdxs = temp1(1) + temp2(1)
         elseif (part.eq.2) then
            cdxs = temp1(2) + temp2(2)
         elseif (part.eq.3) then
            cdxs = temp1(1) + temp2(1)
         elseif (part.eq.4) then
            cdxs = temp1(2) + temp2(2)   
         endif
      end
