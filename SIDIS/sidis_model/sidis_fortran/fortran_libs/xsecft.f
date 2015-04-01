    
      subroutine Xsecft(E0,EP,SINSQ,TG,result,success)
      real E0,EP,SINSQ,result(3)
      integer TG
      real COSSQ,TANSQ,Q2,V,VSQ,TPMV,X,WSQ,tmp,CSMOTT,F2,ST
      real SY,SLP,DSLP,R,DR1,DR2,DR3
      real PMSQ,ALPHA,TPM
      parameter (PMSQ = 8.80369E-1)
      parameter (ALPHA=1.3703604E2)
      parameter (TPM = 1.876559E0)
      real WW1,WW2,WXSECTN
      integer sucess
      success=0
      

      WW2 = 0.0
      WW1 = 0.0
      WXSECTN = 0.0

      if ((TG.ne.0).and.(TG.ne.1)) then
         success=-1
         return
      endif
      
      if ((E0.le.0.).or.(EP.le.0.).or.(SINSQ.le.0.).or.(SINSQ.ge.1.)) 
     c then
         success=-1
         return
      endif
      
      COSSQ = 1.0E0 - SINSQ
      TANSQ = SINSQ/COSSQ
      Q2 = 4.0E0 * E0 * EP * SINSQ
      V = E0 - EP
      VSQ = V * V
      TPMV = TPM * V
      X = Q2/TPMV
      WSQ = PMSQ + TPMV - Q2
      
      if (WSQ.lt.3.) then
         success=-1
         return
      endif
      
      tmp = 1.9732E4 / (2.0E0 * ALPHA * E0 * SINSQ)
      CSMOTT = 1.0E-3 * tmp**2 * COSSQ
      
      call LWW90_F2(X,Q2,TG,12,F2,ST,SY,SLP,DSLP,success)
      if (success.ne.0) then
         return
      endif
      if (TG.eq.1) then
         F2=2.*F2
      endif
      call LWW90_R(X,Q2,R,DR1,DR2,DR3,success)
      if (success.ne.0) then
         return
      endif
      
      WW2 = F2/V
      WW1 = (1.E0 + VSQ/Q2)/(V*(1.E0 + R)) * F2
      WXSECTN = CSMOTT * ( (WW2) + 2.E0*TANSQ* (WW1) )
      
      RESULT(1) = WW2
      RESULT(2) = WW2
      RESULT(3) = WXSECTN
      
      
      
      end
