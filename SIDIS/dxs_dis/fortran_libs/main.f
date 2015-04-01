
      PROGRAM main
      PARAMETER (NWPAWC = 1500000)
      PARAMETER (LRECL  = 8192)
      
      COMMON /PAWC/ IPAW(NWPAWC)      

      integer IQUEST(100)
      common/QUEST/IQUEST
      
      real x,q2,w,wp,z,pt
      real hitnum1_1,hitnum1_2,hitnum1_3,hitnum1_4,hitnum1_5
      real hitnum1_6,hitnum1_7,hitnum1_8,hitnum1_9,hitnum1_10
      real hitnum1_11,hitnum1_12,hitnum1_13
      real hitnum2_1,hitnum2_2,hitnum2_3,hitnum2_4,hitnum2_5
      real hitnum2_6,hitnum2_7,hitnum2_8,hitnum2_9,hitnum2_10
      real hitnum2_11,hitnum2_12,hitnum2_13
      real theta_e,theta_h,p_e,p_h,phi_e,phi_h
      real sdxs_e,sdxs_h,cdxs
      real weight
      
      integer neve,part
      common /general/ neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight

      CALL HLIMIT(NWPAWC)   

      IQUEST(10)=256000
      part=1

      call HROPEN(25,'NTUPEL1','./forward_pip_p_11.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/forward_pip_p_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./forward_pim_p_11.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight') 
      call HROPEN(36,'oldfile','./data/forward_pim_p_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./forward_pip_d_11.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/forward_pip_d_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./forward_pim_d_11.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/forward_pim_d_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)



      part=1

      call HROPEN(25,'NTUPEL1','./large_pip_p_11.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/large_pip_p_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./large_pim_p_11.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight') 
      call HROPEN(36,'oldfile','./data/large_pim_p_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./large_pip_d_11.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/large_pip_d_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./large_pim_d_11.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/large_pim_d_11.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)


      part=1

      call HROPEN(25,'NTUPEL1','./forward_pip_p_8p8.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/forward_pip_p_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./forward_pim_p_8p8.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight') 
      call HROPEN(36,'oldfile','./data/forward_pim_p_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./forward_pip_d_8p8.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/forward_pip_d_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./forward_pim_d_8p8.nt','NQE',LRECL
     c     ,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/forward_pim_d_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)



      part=1

      call HROPEN(25,'NTUPEL1','./large_pip_p_8p8.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/large_pip_p_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./large_pim_p_8p8.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight') 
      call HROPEN(36,'oldfile','./data/large_pim_p_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,0)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./large_pip_d_8p8.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/large_pip_d_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,1,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)

      call HROPEN(25,'NTUPEL1','./large_pim_d_8p8.nt','NQE',LRECL,ISTAT)
      call HBNT(9500,'general','D')
      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt,
     c     hitnum1_1,hitnum1_2,
     c     hitnum2_1,hitnum2_2,
     c     theta_e,theta_h,p_e,p_h,phi_e,phi_h,sdxs_e,sdxs_h
     c     ,cdxs,weight')
      call HROPEN(36,'oldfile','./data/large_pim_d_8p8.nt',' ',1024
     c     ,ISTAT)
      call HRIN(2,9999,0)
      call NT20(part,2,1)
      call HROUT(9500,IC,'')
      call HREND('NTUPEL1')
      call HREND('oldfile')
      close(25)
      close(36)








c$$$      part=2
c$$$
c$$$      call HROPEN(25,'NTUPEL1','./convert_kp_p_8p8.nt','NQE',LRECL,ISTAT)
c$$$      call HBNT(9500,'general','D')
c$$$      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt
c$$$     c     ,hitnum,hitnum1,theta_e,theta_h,p_e,p_h,phi_e,phi_h
c$$$     c     ,sdxs_e,sdxs_h,cdxs,weight')
c$$$      call HROPEN(36,'oldfile','../data/csimout_kp_p_8p8.nt',' ',1024
c$$$     c     ,ISTAT)
c$$$      call HRIN(2,9999,0)
c$$$      call NT20(part,3,0)
c$$$      call HROUT(9500,IC,'')
c$$$      call HREND('NTUPEL1')
c$$$      call HREND('oldfile')
c$$$      close(25)
c$$$      close(36)
c$$$
c$$$      call HROPEN(25,'NTUPEL1','./convert_km_p_8p8.nt','NQE',LRECL,ISTAT)
c$$$      call HBNT(9500,'general','D')
c$$$      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt
c$$$     c     ,hitnum,hitnum1,theta_e,theta_h,p_e,p_h,phi_e,phi_h
c$$$     c     ,sdxs_e,sdxs_h,cdxs,weight')
c$$$      call HROPEN(36,'oldfile','../data/csimout_km_p_8p8.nt',' ',1024
c$$$     c     ,ISTAT)
c$$$      call HRIN(2,9999,0)
c$$$      call NT20(part,4,0)
c$$$      call HROUT(9500,IC,'')
c$$$      call HREND('NTUPEL1')
c$$$      call HREND('oldfile')
c$$$      close(25)
c$$$      close(36)
c$$$
c$$$      call HROPEN(25,'NTUPEL1','./convert_kp_d_8p8.nt','NQE',LRECL,ISTAT)
c$$$      call HBNT(9500,'general','D')
c$$$      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt
c$$$     c     ,hitnum,hitnum1,theta_e,theta_h,p_e,p_h,phi_e,phi_h
c$$$     c     ,sdxs_e,sdxs_h,cdxs,weight')
c$$$      call HROPEN(36,'oldfile','../data/csimout_kp_d_8p8.nt',' ',1024
c$$$     c     ,ISTAT)
c$$$      call HRIN(2,9999,0)
c$$$      call NT20(part,3,1)
c$$$      call HROUT(9500,IC,'')
c$$$      call HREND('NTUPEL1')
c$$$      call HREND('oldfile')
c$$$      close(25)
c$$$      close(36)
c$$$
c$$$      call HROPEN(25,'NTUPEL1','./convert_km_d_8p8.nt','NQE',LRECL,ISTAT)
c$$$      call HBNT(9500,'general','D')
c$$$      call HBNAME(9500,'general',neve,'neve,x,q2,w,wp,z,pt
c$$$     c     ,hitnum,hitnum1,theta_e,theta_h,p_e,p_h,phi_e,phi_h
c$$$     c     ,sdxs_e,sdxs_h,cdxs,weight')
c$$$      call HROPEN(36,'oldfile','../data/csimout_km_d_8p8.nt',' ',1024
c$$$     c     ,ISTAT)
c$$$      call HRIN(2,9999,0)
c$$$      call NT20(part,4,1)
c$$$      call HROUT(9500,IC,'')
c$$$      call HREND('NTUPEL1')
c$$$      call HREND('oldfile')
c$$$      close(25)
c$$$      close(36)
      

      end
