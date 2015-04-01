      SUBROUTINE DXS_DIS(targ_in,part_in,beam_energy_in,
     c           p_e_in,theta_e_in,phi_e_in,p_h_in,theta_h_in,phi_h_in,
     c           sdxs_e,sdxs_h,dxs_hp,dxs_hm)
!     targ: 1->proton,2->deuterium,3->he3      
!     part: particle type
!           1 for pi+;  2 for pi-, 3=k+, 4=k-, 5=p
!     beam_energy,p_e,p_h(GeV),theta_e,phi_e,theta_h,phi_h (DEG) 

      real theta_e_in,theta_h_in,p_e_in,p_h_in,phi_e_in,phi_h_in
      real beam_energy_in
      integer targ_in, part_in

      real theta_e,theta_h,p_e,p_h,phi_e,phi_h
      real beam_energy,factor,pe_temp,ph_temp
      real sdxs_e,sdxs_h,dxs_hp,dxs_hm,temp
      real temp1(2),temp2(2)
      integer targ,part,success
      integer typ ! 0 for pion, 1 for kaon         

      targ = targ_in
      part = part_in
      p_e = p_e_in
      p_h = p_h_in
      theta_e = theta_e_in
      theta_h = theta_h_in
      phi_e = phi_e_in
      phi_h = phi_h_in

      beam_energy = beam_energy_in*1000.0 
      pe_temp = p_e*1000.0
      ph_temp = p_h*1000.0

      if(targ.eq.1) then
          !proton
          call whitlow(0,beam_energy,pe_temp,0.0,pe_temp
     c                  ,theta_e,0.0,theta_e,sdxs_e,success)
      else if(targ.eq.2) then
          !deutron
          call whitlow(1,beam_energy,pe_temp,0.0,pe_temp
     c           ,theta_e,0.0,theta_e,sdxs_e,success)
      else if(targ.eq.3) then
          !he3
          call whitlow(2,beam_energy,pe_temp,0.0,pe_temp
     c            ,theta_e,0.0,theta_e,sdxs_e,success)
      endif   

      if (success.ne.0.or.sdxs_e.le.0.or.sdxs_e.ge.10000.) then
          if(targ.eq.1) then
              !proton
              call qfsrad(beam_energy/1000.,p_e,theta_e/180.*3.1415926
     c           ,0,sdxs_e,success)
          else if(targ.eq.1) then
              !deutron
              call qfsrad(beam_energy/1000.,p_e,theta_e/180.*3.1415926
     c           ,1,sdxs_e,success)
          else if(targ.eq.3) then
              !he3
              call qfsrad(beam_energy/1000.,p_e,theta_e/180.*3.1415926
     c           ,0,sdxs_e,success)
              call qfsrad(beam_energy/1000.,p_e,theta_e/180.*3.1415926
     c           ,1,temp,success)
              sdxs_e=sdxs_e + temp
          endif

          if (sdxs_e.ge.0.and.sdxs_e.le.1000000) then
          else
              sdxs_e = 0
          endif
      endif

      if(targ.eq.1) then
          !proton
          call wiser(part,0,beam_energy,ph_temp,0.1,ph_temp,
     c        theta_h,0.1,theta_h,sdxs_h,success)
          elseif(targ.eq.2) then
          !deutron
          call wiser(part,1,beam_energy,ph_temp,0.1,ph_temp,
     c        theta_h,0.1,theta_h,sdxs_h,success)
          elseif(targ.eq.3) then
          !he3
          call wiser(part,0,beam_energy,ph_temp,0.1,ph_temp,
     c        theta_h,0.1,theta_h,sdxs_h,success)
          call wiser(part,1,beam_energy,ph_temp,0.1,ph_temp,
     c        theta_h,0.1,theta_h,temp,success)
          sdxs_h = sdxs_h+temp
      endif

      factor = (2.33369*exp(-0.508963*p_h*
     c        sin(theta_h/180.*3.1415926)*sqrt(0.938*0.938
     c        +2.*0.938*5.892)/sqrt(0.938*0.938+2.*0.938
     c        *beam_energy/1000.)))

      if (factor.le.1) factor = 1

      sdxs_h = sdxs_h /factor
      if (sdxs_h.le.0.or.success.ne.0) then
          sdxs_h = 0
      endif

      !print*, "1 E0=",beam_energy
      !print*, " Ee=",pe_temp," The=",theta_e," Phe=",phi_e
      !print*, " Eh=",ph_temp," Thh=",theta_h," Phh=",phi_h

      ! proton
      if(part.eq.1 .or. part.eq.2) then
          typ=0
      else if(part.eq.3 .or. part.eq.4) then
          typ=1
      endif

      if(targ.eq.1) then
          ! proton
          CALL SIDIS(0,typ,beam_energy,pe_temp,pe_temp,theta_e,
     c        theta_e,phi_e,phi_e,ph_temp,ph_temp,theta_h
     c        ,theta_h,phi_h,phi_h,temp1)

          dxs_hp = temp1(1)!pi+/k+
          dxs_hm = temp1(2)!pi-/k-
      else if(targ.eq.2) then
          ! deutron
          CALL SIDIS(1,typ,beam_energy,pe_temp,pe_temp,theta_e,
     c        theta_e,phi_e,phi_e,ph_temp,ph_temp,theta_h
     c        ,theta_h,phi_h,phi_h,temp2)

          dxs_hp = temp2(1)!pi+/k+
          dxs_hm = temp2(2)!pi-/k-
      elseif(targ.eq.3) then
          CALL SIDIS(0,typ,beam_energy,pe_temp,pe_temp,theta_e,
     c        theta_e,phi_e,phi_e,ph_temp,ph_temp,theta_h
     c        ,theta_h,phi_h,phi_h,temp1)
          CALL SIDIS(1,typ,beam_energy,pe_temp,pe_temp,theta_e,
     c        theta_e,phi_e,phi_e,ph_temp,ph_temp,theta_h
     c        ,theta_h,phi_h,phi_h,temp2)

          dxs_hp = temp1(1) + temp2(1)!pi+/k+
          dxs_hm = temp1(2) + temp2(2)!pi-/k-
      endif

      if (dxs_hp .lt. 0. .or. dxs_hm .lt. 0.) then
          print*, "OHHHHHH, Wrong SIDIS DXS!"
      endif

      end
