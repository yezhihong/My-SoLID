      REAL FUNCTION HITRAXY(IXY,IDET,ICH,ITR)
C
C===    Returns the X/Y-coordinate MRS/DRS of the hit originated from
C===            IXY=1,2  - X,Y DRS
C                  -1,-2 - X,Y MRS
C                  -3    - Z MRS
C                  -4    - P mom
C                  -5    - Slope X DRS
C                  -6    - Slope Y DRS
C                  -7    - Full ionization in the detector
C                  -8    - Number of hits in the detector
C                   3    - R   sqrt(X**2+Y**2) in DRS
C                   4    - Phi atan2(Y,X) in DRS
C                   5    - R   sqrt(X**2+Y**2) in MRS
C                   6    - Phi atan2(Y,X) in MRS
C===            the track ITR (the original) 0 - any (last)
C===            in the counter IDET, channel ICH
C===            ICH=0 -  1-st hit
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


c$$$      LOGICAL         CHAIN
c$$$      CHARACTER*128   CFILE
c$$$      INTEGER         IDNEVT,NCHEVT,ICHEVT
c$$$      REAL            OBS(13)
c$$$*
c$$$      COMMON /PAWIDN/ IDNEVT,OBS
c$$$      COMMON /PAWCHN/ CHAIN, NCHEVT, ICHEVT
c$$$      COMMON /PAWCHC/ CFILE
*
*--   Ntuple Variable Declarations
*
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


      INTEGER  IXY,IDET,ICH,ITR
      
C
      COMMON/CHITRAXY/ CHIT(8)
      REAL CHIT
C
      real VHLIM(2)
      INTEGER ifirst
      DATA ifirst/1/
C
      INTEGER i,k,idi,i1,i2,ihi,iddd,itrh,iorh,nfind,idg1,idg2,j,ndhits
      REAL sign_ion
C
C     -----------------------------------------------------------------
C
      IF(ifirst.EQ.1) THEN
         VHLIM(1)= 9.E20
         VHLIM(2)=-9.E20
      ENDIF
      ifirst=0
C
      IF(IXY.EQ.-8) THEN
         HITRAXY=0.
      ELSE
         HITRAXY=-300.
      ENDIF
      IF(nhit.LT.1) GO TO 999
      IF(ndig.LT.1) GO TO 999
      IF(IABS(IXY).LT.1.OR.IABS(IXY).GT.8) GO TO 999
      nfind=0
      sign_ion=0
      ndhits=0
C
      IF(npdig.GT.0) THEN
         DO idi=1,ndig
            k=ip1dig(idi)
            iddd=IBITS(k,0,16)
            IF(iddd.EQ.IDET) THEN
               k=ip2dig(idi)
               idg1=IBITS(k,0,16)
               idg2=IBITS(k,16,16)
C                     write(6,*) idi,iddd,idg1
               IF(idg1.EQ.ICH.OR.ICH.EQ.0) THEN
                  i1=1
                  IF(idi.GT.1) THEN
                     k=ip1dig(idi-1)
                     i1=IBITS(k,16,16)
                  ENDIF
                  k=ip1dig(idi)
                  i2=IBITS(k,16,16)
                  DO j=i1,i2-1
                     ihi=jpdig(j)
                     k=ip1hit(ihi)
                     itrh=IBITS(k,0,16)
                     iorh=IBITS(k,16,16)
C
C                     write(6,*) ' itrh..=',itrh,ITR,iorh
                     IF((itrh.EQ.ITR.AND.iorh.EQ.0).OR.ITR.EQ.0) THEN
                        nfind=nfind+1
C                        WRITE(6,*) 
C     +                        IDNEVT,nfind,idi,ihi,itrh,iorh,i1,i2
                        CHIT(1)=hitd(1,ihi)
                        CHIT(2)=hitd(2,ihi)
                        CHIT(3)=hit (1,ihi)
                        CHIT(4)=hit (2,ihi)
                        CHIT(5)=zhit(ihi)
                        CHIT(6)=phit(ihi)
                        CHIT(7)=hitsld(1,ihi)
                        CHIT(8)=hitsld(2,ihi)
                        ndhits=ndhits+1
C                        write(6,*) 'hit=',ihi,hition(ihi),hit(1,ihi)
C     +                         ,hit(2,ihi)
                        sign_ion=sign_ion+hition(ihi)
                        IF(IXY.GT.0) THEN
                           IF(IXY.LE.2) THEN
                              HITRAXY=hitd(IXY,ihi)
                           ELSE IF(IXY.EQ.3) THEN
                              HITRAXY=SQRT(hitd(1,ihi)**2
     +                                    +hitd(2,ihi)**2)
C                              write(6,*) 'eve,r=',IEVE,HITRAXY
                           ELSE IF(IXY.EQ.4) THEN
                              HITRAXY=ATAN2(hitd(2,ihi),hitd(1,ihi))
     +                                        *180./3.1416
                              IF(HITRAXY.LT.-90.) HITRAXY=HITRAXY+360.
C                              IF(HITRAXY.LT.0.) HITRAXY=HITRAXY+360.
C                              HITRAXY=MOD(HITRAXY,20.)
                           ELSE IF(IXY.EQ.5) THEN
                              HITRAXY=SQRT(hit(1,ihi)**2
     +                                    +hit(2,ihi)**2)
C                              write(6,*) 'eve,r=',IEVE,HITRAXY
                           ELSE IF(IXY.EQ.6) THEN
                              HITRAXY=ATAN2(hit(2,ihi),hit(1,ihi))
     +                                        *180./3.1416
                              IF(HITRAXY.LT.-90.) HITRAXY=HITRAXY+360.
                           ENDIF
                        ELSE 
                           IF(IABS(IXY).LT.3) THEN
                              HITRAXY=hit(-IXY,ihi)
                           ELSE IF(IABS(IXY).EQ.3) THEN
                              HITRAXY=zhit(ihi)
                           ELSE IF(IABS(IXY).EQ.4) THEN
                              HITRAXY=phit(ihi)
                           ELSE IF(IABS(IXY).EQ.5) THEN
                              HITRAXY=hitsld(1,ihi)
                           ELSE IF(IABS(IXY).EQ.6) THEN
                              HITRAXY=hitsld(2,ihi)
                           ENDIF
                        ENDIF
                     ENDIF
C
                  END DO
               ENDIF
            ENDIF
         END DO
      ENDIF
C
      IF(IXY.EQ.-7) THEN
         HITRAXY=sign_ion
      ELSE IF(IXY.EQ.-8) THEN
         HITRAXY=ndhits
      ELSE
         IF(nfind.GT.1) WRITE(6,*) ' == ev=',IDNEVT
     +              ,' Too many hits found=',nfind
     +              ,' IXY,IDET,ICH,ITR=',IXY,IDET,ICH,ITR
      ENDIF
C
      IF(HITRAXY.GT.-299.) THEN
         VHLIM(1)=MIN(VHLIM(1),HITRAXY)
         VHLIM(2)=MAX(VHLIM(2),HITRAXY)
      ENDIF
C
 999  RETURN
C
      END


