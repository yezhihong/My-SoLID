      REAL FUNCTION GETKINP(IFL,ITR)
C
C ===  Gets the parameters of a particle ITR
C===    IFL=0 - full momentum
C===       =1 - enrgy
C===       =2 - TAN(theta)
C===       =3 - THETA (deg)
C===       =4 - PHI (deg)
C===       =5 - X Bjorken
C===       =6 - Q2
C===       =7 - W2
C
      
c$$$      INTEGER ieve,irun,iend
c$$$      COMMON /RUN/ ieve,irun,iend
c$$$*
c$$$      REAL bpara
c$$$      INTEGER ibtyp,ibfla
c$$$      COMMON /BEAM/ ibtyp,ibfla,bpara(6)
c$$$*
c$$$      REAL vert,ptra,xkinc,pkinc,thetcm,phicm,wgcros,anpower
c$$$      INTEGER nver,igev,ltimv,imov,ntdv,itdv,ntra,iget,itra,itvb,itve
c$$$     + ,nkinc,ietype
c$$$      COMMON /KINE/ nver,igev(120),vert(3,120),ltimv(120),imov(120)
c$$$     + ,ntdv(120),itdv(120),ntra,iget(250),ptra(3,250),itra(250)
c$$$     + ,itvb(250),itve(250),nkinc,xkinc(2),pkinc(3,2),ietype,thetcm
c$$$     + ,phicm,wgcros,anpower
c$$$*
c$$$      REAL hit,zhit,hitd,phit,hitsld,hition
c$$$      INTEGER nhit,nhitall,ip1hit,ip2hit
c$$$      COMMON /HIT/ nhit,nhitall,ip1hit(2000),ip2hit(2000),hit(2,2000)
c$$$     + ,zhit(2000),hitd(3,2000),phit(2000),hitsld(2,2000),hition(2000)
c$$$*
c$$$      INTEGER ndig,ndigall,ip1dig,ip2dig,npdig,jpdig
c$$$      COMMON /DIG/ ndig,ndigall,ip1dig(1299),ip2dig(1299),npdig
c$$$     + ,jpdig(1428)

          REAL bpara(6),vert(3,120),ptra(3,250),trlast(3,250),xkinc(2)
     + ,pkinc(3,2),thetcm,phicm,wgcros,anpower,hit(2,2000),zhit(2000)
     + ,hitd(3,2000),phit(2000),hitsld(2,2000),hition(2000)
      INTEGER ieve,irun,iend,ibtyp,ibfla,nver,igev(120),ltimv(120)
     + ,imov(120),ntdv(120),itdv(120),ntra,iget(250),itra(250),itvb(250)
     + ,itve(250),itrstop(250),nkinc,ietype,nhit,nhitall,ip1hit(2000)
     + ,ip2hit(2000),ndig,ndigall,ip1dig(1299),ip2dig(1299),npdig
     + ,jpdig(1428)
*
      COMMON /PAWCR4/ ieve,irun,iend,ibtyp,ibfla,bpara,nver,igev,vert
     + ,ltimv,imov,ntdv,itdv,ntra,iget,ptra,itra,itvb,itve,trlast
     + ,itrstop,nkinc,xkinc,pkinc,ietype,thetcm,phicm,wgcros,anpower
     + ,nhit,nhitall,ip1hit,ip2hit,hit,zhit,hitd,phit,hitsld,hition,ndig
     + ,ndigall,ip1dig,ip2dig,npdig,jpdig

      INTEGER  IFL,ITR
      
C
      real VLIM(2)
C
      INTEGER i,ityp
      REAL ptot,en,am(20),ams,res,thet,sl,pt,phi,eb,q2,xbj,w2
      INTEGER ifirst
      DATA ifirst/1/

      DATA am/0.,0.511E-3,0.511E-3,0.,0.105,0.105,0.134976
     +    ,0.1395700,0.1395700,0.497672,0.493677,0.493677
     +    ,0.93956563,0.9382723,0.9382723,5*0./
C
C     -----------------------------------------------------------------
C
      IF(ifirst.EQ.1) THEN
         VLIM(1)= 9.E20
         VLIM(2)=-9.E20
      ENDIF
      ifirst=0

      GETKINP=-9999.
      IF(IFL.LT.0.OR.IFL.GT.7) GO TO 999
      IF(ITR.LT.1.OR.ITR.GT.ntra) GO TO 999
C
      ptot=SQRT(ptra(1,ITR)**2+ptra(2,ITR)**2+ptra(3,ITR)**2)
      ityp=itra(ITR)
C       ityp=3
      IF(ityp.LT.1.OR.ityp.GT.20) GO TO 999
      ams=am(ityp)
      en=SQRT(ptot**2+ams**2)
      pt=SQRT(ptra(1,ITR)**2+ptra(2,ITR)**2)
      sl=pt/ptra(3,ITR)
      thet=ATAN2(pt,ptra(3,ITR))
      phi=ATAN2(ptra(2,ITR),ptra(1,ITR))
      eb=ABS(ptra(3,1))
      q2=4.*eb*ptot*SIN(thet/2.)**2
      xbj=q2/2./0.938/(eb-ptot) ! for electrons
      w2=0.938**2+2.*0.938*(eb-ptot)-q2
C
      

      IF(IFL.EQ.0) THEN
         res=ptot
      ELSE IF(IFL.EQ.1) THEN
         res=en
      ELSE IF(IFL.EQ.2) THEN
         res=sl
      ELSE IF(IFL.EQ.3) THEN
         res=thet*180./3.1416
      ELSE IF(IFL.EQ.4) THEN
         res=phi*180./3.1416
         IF(res.LT.-90.) res=360.+res
      ELSE IF(IFL.EQ.5) THEN
         res=xbj
      ELSE IF(IFL.EQ.6) THEN
         res=q2
      ELSE IF(IFL.EQ.7) THEN
         res=w2
      ENDIF
      GETKINP=res
      VLIM(1)=MIN(VLIM(1),res)
      VLIM(2)=MAX(VLIM(2),res)
C
 999  RETURN
C
      END


