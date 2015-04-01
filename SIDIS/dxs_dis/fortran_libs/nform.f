


c this version of nform was copied from [adam.model] feb 24 1988
c by klg

      SUBROUTINE NFORM(IG,RMOD,IMOD,QQ,GEP,GEN,GMP,GMN)   
C 
C  ---------------------------------------------------------- 
C 
C   CALCULATE NUCLEON FORM FACTORS
C 
C 
C   QQ = INPUT Q SQUARED (FM**-2) 
C	RMOD IS MODIFIED NUCLEON RADIUS TO BE USED WITH IG=11
C	IMOD IS FLAG WHETHER TO MODIFY THE FORM FACTOR
C	IMOD = 0 >> NO MODIFICATON
C	IMOD = 1 >> MODIFICATION BY THE TERM FACTOR BELOW 
C  ------------------------------------------------------------ 
      IMPLICIT REAL*8 (A-H,O-Z)
C 
      COMMON/ROSEIO/INPUT,IOUT
C 
C 
      DIMENSION TITLE(15)
      CHARACTER*40 TITLE
      DATA TITLE(1)/'EMPERICAL DIPOLE WITH GEN = 0.0         '/, 
     *     TITLE(2)/'IJL 5 PARAMETER MODEL DIPOLE FIT        '/, 
     *     TITLE(3)/'BEST FIT=GEP+GMP=IJL; GMN=GD;GEN=GALSTER'/, 
     *     TITLE(4)/'BEST FIT EXCEPT GEN= 0.0                '/, 
     *     TITLE(5)/'POINT NUCLEONS                          '/, 
     *     TITLE(6)/'HOEHLER,NPB114,505                      '/, 
     *     TITLE(7)/'BLATNIK + ZOVKO VDM FIT                 '/, 
     *     TITLE(8)/'JANSSENS 1966 STANDARD FIT              '/,
     *     TITLE(9)/'HOEHLER FIT NO. 5.1                     '/,      
     *     TITLE(10)/'HOEHLER FIT NO. 5.3                    '/,
     *	   title(11)/'Av. of BZ and Dipole for Gmn; BZ else  '/,
     *	   title(12)/'Simon GEP,GMN NPA333(1980)381; BZ else '/,
     *	   title(13)/'Simon-Gep,Gmp;Gen BZ;Gmn Av. BZ+dipole '/,
     *     TITLE(14)/'HOEHLER*NUCLEON SIZE CHANGE *GE ONLY*  '/,
     *     TITLE(15)/'DIPOLE WITH MODIFIED MAGNETIC MOMENT   '/
C_______________________________________________________________ 
C   IJL PARAMETERS FOR 5 PARAMETER DIPOLE FIT (IN GEV UNITS)
C   PHYS LETT. 43B, 191(1973) 
      DATA GAM   ,BR    ,BW    ,BF    ,AF   
     +   /0.25  ,0.672 ,1.102 ,0.112 ,-.052 /   
      DATA RMN2  ,RMW2  ,RMF2  ,RMR2  ,GAMR  ,PI
     +   /0.8817,0.6146,1.0384,0.5852,0.112 ,3.14159/   
      DATA RMPI  ,RMPI2 
     +   /.139  ,.019321/   
C     BZ
      DATA TRO   ,TROP  ,TROPP ,TFI   ,TOM   ,TOMP  
     +   /0.585 ,1.30  ,2.10  ,1.039 ,0.614 ,1.40  /
      DATA RMUS  ,RMUV  ,BS    ,BV  
     +   /-.060 ,1.853 ,-.91  ,-1.10 /  
      DATA RMUP  ,RMUN  
     +   /2.79278,-1.91315/ 
!	data for simon parametrization
C  Les masses sont en fm-2
      DATA A1S/0.312/,A2S/1.312/,A3S/-.709/,A4S/0.085/
      DATA A1V/0.694/,A2V/0.719/,A3V/-.418/,A4V/0.005/
      DATA AM1/6./,AM2/15.02/,AM3/44.08/,AM4/154.2/
      DATA AM5/8.5/,AM6/355.4/,AMUP/2.793/
!	end data for simon
!
C     QQ IN FM-2,QQG IN GEV/C**2
C
C      print *,'nform: begin'

	IF (IG .EQ. 15)THEN
	RMUP = 3.631
	RMUN = -2.487
	ENDIF
      QQG=QQ*0.197328**2
C      print *, QQG, 'GeV/c**2', QQ,'fm-2'
      	IF (IMOD .EQ. 0)THEN
	FACTOR = 1.0
	ELSE
	IF (IMOD .EQ. 1) THEN
	FACTOR = (1.+QQ*0.80**2/12.)**2/(1.+QQ*RMOD**2/12.)**2
	ELSE
	STOP 'Meaningless value of IMOD stopped in NFORM'
	ENDIF
	ENDIF
      GOTO(110,120,120,120,150,160,170,180,190,200,110,
     >170,110,210,220),IG  
C     DIPOLE
  110 AAA=0.71  
      GEP=1./(1.+QQG/AAA)**2
      GEN=0.
      GMP=RMUP*GEP  
      GMN=RMUN*GEP  
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
	if(ig .eq. 13)goto 230
	if(ig .eq. 11)goto 230
      GOTO 900  
C     IJL   
  120 TAU=QQG/(4.*RMN2) 
      GT=0.5/(1.+GAM*QQG)**2
      T1=SQRT(QQG+4.*RMPI2) 
      T2=SQRT(QQG)  
      ALPH=2.*T1*LOG((T1+T2)/(2.*RMPI))/(T2*PI)
      TOP=RMR2+8.*GAMR*RMPI/PI  
      BOT=RMR2+QQG+(4.*RMPI2+QQG)*GAMR*ALPH/RMPI
      RHO=TOP/BOT   
      F1S=GT*((1.-BW-BF)+BW/(1.+QQG/RMW2)+BF/(1.+QQG/RMF2)) 
      F1V=GT*((1.-BR)+BR*RHO)   
      F2S=GT*((-0.12-AF)/(1.+QQG/RMW2)+AF/(1.+QQG/RMF2))
      F2V=GT*(3.706*RHO)
      GEP=F1V+F1S-TAU*(F2V+F2S)
      GEN=F1S-F1V-TAU*(F2S-F2V)
      GMP=F1V+F1S+F2V+F2S
      GMN=F1S-F1V+F2S-F2V
      IF(IG.EQ.2) GOTO 900  
      GD=1./(1.+QQG/.71)**2 
      GMN=RMUN*GD   
      GEN=-RMUN*TAU*GD/(1.+5.6*TAU) 
      IF(IG.EQ.3) GOTO 900  
      GEN=0.
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900  
  150 CONTINUE  
      GMN=RMUN  
      GEP=1.
      GEN=0.
      GMP=RMUP  
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900  
  160 CONTINUE  
C     HOEHLER FIT 8.2   
      TAU=QQG/(4.*RMN2) 
      F1S= 0.71/(0.613+QQG)-0.64/(1.040+QQG)-0.13/(3.24+QQG)
      F2S= -0.11/(0.613+QQG)+0.13/(1.040+QQG)-0.02/(3.24+QQG)   
      F1V=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536) 
      F2V=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)   
      F1V=F1V+0.05/(1.464+QQG)-0.52/(6.003+QQG)+0.28/(8.703+QQG)
      F2V=F2V-1.99/(1.464+QQG)+0.20/(6.003+QQG)+0.19/(8.703+QQG)
      GEP=F1V+F1S-TAU*(F2V+F2S)
      GEN=F1S-F1V-TAU*(F2S-F2V)
      GMP=(F1V+F1S+F2V+F2S)
      GMN=(F1S-F1V+F2S-F2V)
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900  
C     BZ
  170 TAU=QQG/(4.*RMN2) 
      RS=1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))   
      RV=1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP)) 
      F1E=(0.5-TAU*(RMUS+2.*RMN2*BS))*RS
      F2E=(0.5-TAU*(RMUV+2.*RMN2*BV))*RV
      F1M=(0.5+RMUS-0.5*BS*QQG)*RS  
      F2M=(0.5+RMUV-0.5*BV*QQG)*RV  
      GEP=F1E+F2E   
      GMP=F1M+F2M   
      GEN=F1E-F2E   
      GMN=F1M-F2M   
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
	if (ig .eq. 12)goto 240	! if simon use gmn,gen from B+Z
      GOTO 900  
C     JANSENS   
  180 F1=1.+QQ/15.7 
      F2=1.+QQ/26.7 
      F3=1.+QQ/8.19 
      GES=0.5*(2.5/F1-1.6/F2+0.1)   
      GMS=0.44*(3.33/F1-2.77/F2+0.44)   
      GEV=0.5*(1.16/F3-0.16)
      GMV=2.353*(1.11/F3-0.11)  
      GEP=GES+GEV   
      GMP=GMS+GMV   
      GEN=GES-GEV   
      GMN=GMS-GMV   
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GO TO 900
C     HOELER FIT 5.1
  190 TAU=QQG/(4.*RMN2) 
      F1=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536) 
      F2=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)   
      F1P=F1+0.63/(0.613+QQG)-0.43/(1.124+QQG)-0.45/(2.723+QQG)
      F2P=F2+0.02/(0.613+QQG)-1.89/(1.323+QQG)+0.27/(7.952+QQG)
      GEP=F1P-TAU*F2P
      GEN=0.0
      GMP=F1P+F2P
      GMN=0.0
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
C     HOEHLER FIT 5.3
  200 TAU=QQG/(4.*RMN2) 
      F1=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536) 
      F2=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)   
      F1P=F1+0.67/(0.613+QQG)-0.39/(0.922+QQG)-0.54/( 2.756+QQG)
      F2P=F2+0.04/(0.613+QQG)-1.88/(1.300+QQG)+0.24/(10.176+QQG)
      GEP=F1P-TAU*F2P
      GEN=0.0
      GMP=F1P+F2P
      GMN=0.0
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900
C     HOEHLER FIT 8.2   
C	Here we modify only the electric form factors leave the
C	magnetic unchanged so this is ig=6 with ge modified
210   TAU=QQG/(4.*RMN2) 
      F1S= 0.71/(0.613+QQG)-0.64/(1.040+QQG)-0.13/(3.24+QQG)
      F2S= -0.11/(0.613+QQG)+0.13/(1.040+QQG)-0.02/(3.24+QQG)   
      F1V=0.5*(0.955+0.09/(1.+QQG/0.355)**2)/(1.+QQG/0.536) 
      F2V=0.5*(5.335+0.962/(1.+QQG/0.268))/(1.+QQG/0.603)   
      F1V=F1V+0.05/(1.464+QQG)-0.52/(6.003+QQG)+0.28/(8.703+QQG)
      F2V=F2V-1.99/(1.464+QQG)+0.20/(6.003+QQG)+0.19/(8.703+QQG)
      GEP=F1V+F1S-TAU*(F2V+F2S)
      GEN=F1S-F1V-TAU*(F2S-F2V)
      GMP=(F1V+F1S+F2V+F2S)
      GMN=(F1S-F1V+F2S-F2V)
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
!     GMP=FACTOR*GMP
!     GMN=FACTOR*GMN
      GOTO 900  
!	here all is modified by we also change the magnetic moment
  220 AAA=0.71  
C      print *,'dipole used for proton'
      GEP=1./(1.+QQG/AAA)**2
      GEN=0.
      GMP=RMUP*GEP  
      GMN=RMUN*GEP  
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
      GOTO 900  
  230 TAU=QQG/(4.*RMN2) 
      RS=1./((1.+QQG/TOM)*(1.+QQG/TFI)*(1.+QQG/TOMP))   
      RV=1./((1.+QQG/TRO)*(1.+QQG/TROP)*(1.+QQG/TROPP)) 
      F1E=(0.5-TAU*(RMUS+2.*RMN2*BS))*RS
      F2E=(0.5-TAU*(RMUV+2.*RMN2*BV))*RV
      F1M=(0.5+RMUS-0.5*BS*QQG)*RS  
      F2M=(0.5+RMUV-0.5*BV*QQG)*RV  
      GEP=F1E+F2E   
      GMP=F1M+F2M   
      GEN=F1E-F2E   
      GMN=(F1M-F2m + gmn)/2.	! using gmn from dipole to av with BZ gmn
      GEP=GEP*FACTOR
      GEN=FACTOR*GEN
      GMP=FACTOR*GMP
      GMN=FACTOR*GMN
c        print *, QQG,GEP,GEN,GMP,GMN

	if(ig .eq. 13) goto 240
      GOTO 900  
240	continue
!	form factor of simon
!	fit of simon NPA 333 (1980) 381
!	parametrizations in terms of 4 poles
!	use gmn,gen from b&z above
	q2 = qq	!momentum transfer in fm-2
      GEP=A1S/(1.D0+Q2/AM1)+A2S/(1.D0+Q2/AM2)+A3S/(1.D0+Q2/AM3)
     ++A4S/(1.D0+Q2/AM4)
      GMP=A1V/(1.D0+Q2/AM5)+A2V/(1.D0+Q2/AM2)+A3V/(1.D0+Q2/AM3)
     ++A4V/(1.D0+Q2/AM6)
      GMP=GMP*AMUp
	goto 900

!
!
C
C      ENTRY NFORMI(IG,RMOD,IMOD) 
C  ---------------------------------------------------------- 
C 
C   ENTRY POINT TO WRITE FORM FACTOR LABEL IN THE OUTPUT
C 
C  ---------------------------------------------------------- 
C 
	iout=6
c        print *, QQG,GEP,GEN,GMP,GMN

*      print *,'nformi: end'

      WRITE(iout,9000)IG,TITLE(IG)
 9000 FORMAT(1H0,'IG = ',I5,4X,'NUCLEON FORM FACTORS USED = ',A40)

C	IF (IMOD .EQ. 1)WRITE(iout,9001)RMOD
C9001	FORMAT(2X,' These form factors have been modified by a ',
C     >	'change in the nucleon size: RMOD = ',F5.2)
C 
c        print *, QQ,GEP,GEN,GMP,GMN

c      print *,'nformi: end'

c      print *,'nform: end'
 900  RETURN
      END 
