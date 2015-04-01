      subroutine LWW90_R(X,Q2,R,DR1,DR2,DR3,success)
      real X,Q2,R,DR1,DR2,DR3
      integer success
      real FAC,RLN,Q2thr,R_A,R_B,R_C,tmp
      real A(3) /6.723E-2, 4.6714E-1, 1.89794E0/
      real B(3) /6.347E-2, 5.7468E-1, -3.5342E-1/
      real C(3) /5.992E-2, 5.0885E-1, 2.10807E0/
      success=0
      
      FAC = 1.E0 + 1.2E1 * (Q2/(1.E0 + Q2)) * ((1.25E-1**2)/(X**2 + 
     c (1.25E-1)**2))
      RLN = FAC/log(Q2/4.E-2)
      Q2thr = 5.E0 *(1.E0 - X)**5
      
      tmp = Q2**4 + A(3)**4  
      R_A = A(1) * RLN + A(2)/tmp**(2.5E-1)
  
      R_B = B(1) * RLN + B(2)/Q2 + B(3)/(Q2**2 + 3.0E-1**2);
  
      tmp = (Q2 - Q2thr)**2 + C(3)**2
      R_C = C(1) * RLN + C(2)/tmp**(5.0E-1)
      
      R = (R_A + R_B + R_C)/3.E0
  
      call LWW90_DR1(X,Q2,DR1)
  
      tmp = ((R_A -R)**2 +(R_B -R)**2 +(R_C -R)**2)/2.E0
      DR2 = tmp**(5.0E-1)
      DR3 = 2.3E-2 * (1.E0 + 5.0E-1 * R)
      if((Q2 .lt. 1.E0).or. (X.lt. 1.0E-1)) then
         DR3 = 1.5E0 *DR3
      endif
      if(Q2 .lt. 3.0E-1) then
         success = 1
      endif
      end
      

      subroutine LWW90_DR1(X,Q2,DR1)
      real X,Q2,DR1
      real Q2max,Q2w,S,A,Xl,Xh,tmp,Xmin
      Q2max=6.4E1
      Xmin=1.0E-1
      Q2w=Q2
      if(Q2w .gt. Q2max) then
         Q2w = Q2max
      endif
      S = 6.0E-3 + 3.0E-2 *X**2
      tmp = 8.33E0 * X - 6.6E-1
      A = 5.0E-2
      if(A .lt.tmp) then
         A = tmp
      endif
      Xl = 2.0E-2 + abs(S * log(Q2w/A))
      tmp = X
      if(tmp .lt. Xmin) then
         tmp = Xmin
      endif
      Xh = 1.0E-1 *tmp**(2.0E1)/((8.6E-1)**(2.0E1)+tmp**(2.0E1))
      DR1 = sqrt(Xl**2 + Xh**2)
      end
