
*SIG2N
      REAL FUNCTION SIG2N(E,TH,W,Z,A,PF)
      IMPLICIT REAL (A-H,O-Z)
      PI=ACOS(-1.)
      THR=TH*PI/180.
      DM=1232.
      PIMASS=140.
      PM=940.
      A2=550.
      PFR=60.
      GAM2N=20.
      GAMQFR=40.
      GAMREF=300.
      GAMR=GAMREF
      SIGREF=0.20D-7
      QMSR=4.*596.8*(596.8-380.)*SIN(60.*PI/180./2.)**2
      QVSR=QMSR+380.**2
      SIGKIN=0.5*SIGMOT(596.8,60.*PI/180.)
      SIGKIN=SIGKIN*(QMSR/2./QVSR+TAN(60.*PI/180./2.)**2)
      SIGKIN=SIGKIN*QVSR*FD(QMSR,A2)**2
      SIGKIN=SIGKIN*GAMR/GAMREF
      SIGCON=SIGREF/SIGKIN
      QMS=4.*E*(E-W)*SIN(THR/2.)**2
      QVS=QMS+W**2
      GAMQF=GAMQFR*(PF/PFR)*(SQRT(QVS)/SQRT(QVSR))
      EFFMASS=(PM+DM)/2.
      SIG=(Z*(A-Z)/A)*SIGMOT(E,THR)
      SIG=SIG*(QMS/2./QVS+TAN(THR/2.)**2)
      SIG=SIG*QVS*FD(QMS,A2)**2
      EKAPPA=W-QMS/2./PM
      CMTOT2=PM**2+2.*PM*EKAPPA
C     GAM=SQRT(GAMR**2+GAMQF**2)
      GAM=GAMR
      SIG=SIG*CMTOT2*GAM**2
      SIG=SIG/((CMTOT2-EFFMASS**2)**2+CMTOT2*GAM**2)
      SIG=SIG*(GAMR/GAM)*SIGCON
      SIG2N=SIG
      WTHRESH=QMS/4./PM
      IF(W.GT.WTHRESH)THEN
      THRESH=1.-EXP(-(W-WTHRESH)/GAM2N)
      ELSE
      THRESH=0.
      ENDIF
      SIG2N=SIG2N*THRESH
      RETURN
      END
