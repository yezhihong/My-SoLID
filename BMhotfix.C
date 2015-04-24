//#include "formfactors.h"
#include "TMath.h"
#include "BMhotfix.h"
#include "TArrayD.h"
#include <iostream>

Double_t F1p(Double_t q2) { return ((-1.)*q2/4./TMath::Power(0.938271998,2.)*GMp(q2)+GEp(q2))/(1.-q2/4./TMath::Power(0.938271998,2.)); }
Double_t F2p(Double_t q2) { return ( GMp(q2)-GEp(q2))/(1.-q2/4./TMath::Power(0.938271998,2.)); }

Double_t GMp(Double_t q2) { return KellyM(q2); }
Double_t GEp(Double_t q2) { return KellyE(q2); }


Double_t KellyE(Double_t q2)
  /* JJ Kelly PRC 70, 068202 (2004)
   */
{
  Int_t ia=1, ib=3;
  Double_t a[1]={-0.24};
  Double_t b[3]={10.98, 12.82, 21.97};
  Double_t Mp = 0.938;
  Double_t tau = -q2/(4.*pow(Mp,2));
  Double_t GEKn = 1.0;
  Double_t GEKd = 1.0;
  for (int ja=0; ja<ia; ja++){
    GEKn+=a[ja]*pow(tau,ja+1);
  }
  for (int jb=0; jb<ib; jb++){
    GEKd+=b[jb]*pow(tau,jb+1);
  }
  return GEKn/GEKd;
}

Double_t KellyM(Double_t q2)
  /* JJ Kelly PRC 70, 068202 (2004)
     Magnetic Form Factor fit
     Returned value is ratio to dipole*mu_p
  */
{
  Int_t ia=1, ib=3;
  Double_t a[1]={0.12};
  Double_t b[3]={10.97, 18.86, 6.55};Double_t Mp = 0.938;
  Double_t tau = -q2/(4.*pow(Mp,2));
  Double_t GMKn = 1.0;
  Double_t GMKd = 1.0;
  for (int ja=0; ja<ia; ja++)
    {
      GMKn+=a[ja]*pow(tau,ja+1);
    }
  for (int jb=0; jb<ib; jb++)
    {
      GMKd+=b[jb]*pow(tau,jb+1);
    }
  return 2.79285*GMKn/GMKd;
}


Double_t BH(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi){
  //From BMK, BH unpolarized target case ( 
  yy=q2/2./M/x/kz;
  phi=TMath::Pi()-phi; // conversion phi_trento a phi_BMK
  TArrayD coeff(3);
  F1_2=TMath::Power(F1p(t),2);
  F2_2=TMath::Power(F2p(t),2);
  F_sum_2=TMath::Power(F2p(t)+F1p(t),2);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  J=(1-yy-yy*eps_2/2.)*(1+t/q2)-(1-x)*(2-yy)*t/q2;
  
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  
  K_2=-t*(1-x)*(1-yy-TMath::Power(yy*eps/2.,2))*(1-t_min/t)*(TMath::Sqrt(1+eps_2)+(4*x*(1-x)+eps_2)/(4*(1-x))*(t-t_min)/q2)/q2;
  
  K=TMath::Sqrt(K_2);

  P1=-(J+2*K*TMath::Cos(phi))/(yy*(1+eps_2)); 
  P2=1+t/q2+(J+2*K*TMath::Cos(phi))/(yy*(1+eps_2));

  prefactor=hbar_c*TMath::Power(alpha,3)/(TMath::Pi()*8*q2*TMath::Sqrt(1+eps_2))/x/yy/TMath::Power(1+eps_2,2)/t/P1/P2;
  Jacobian=1/(2*x*M*kz);
  prefactor=Jacobian*prefactor;

  coeff[0]=
    8*K_2*((2+3*eps_2)*q2*(F1_2-t*F2_2/(4*M_2))/t+2*x_2*F_sum_2)
    +TMath::Power(2-yy,2)*((eps_2+2)*(4*x_2*M_2/t*TMath::Power(1+t/q2,2)+4*(1-x)*(1+x*t/q2))*(F1_2-t/4./M_2*F2_2)+4*x_2*(x+(1-x+eps_2/2.)*TMath::Power(1-t/q2,2)-x*(1-2*x)*TMath::Power(t/q2,2))*F_sum_2)
    +8*(1+eps_2)*(1-yy-TMath::Power(eps*yy/2.,2))*(2*eps_2*(1-t/4./M_2)*(F1_2-t/4./M_2*F2_2)-x_2*TMath::Power(1-t/q2,2)*F_sum_2);

  coeff[1]=8*K*(2-yy)*((4*x_2*M_2/t-2*x-eps_2)*(F1_2-t/4./M_2*F2_2)+2*x_2*(1-(1-2*x)*t/q2)*F_sum_2);

  coeff[2]=8*x_2*K_2*(4*M_2/t*(F1_2-t/4./M_2*F2_2)+2*F_sum_2);
  return prefactor*(coeff[0]+TMath::Cos(phi)*coeff[1]+TMath::Cos(2*phi)*coeff[2]);
}

TArrayD DVCS(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi, Double_t hel){
  yy=q2/2./M/x/kz;
  phi=TMath::Pi()-phi; // conversion phi_trento to phi_BMK
  //From BM hotfix, DVCS unpolarized target case
  TArrayD coeff(4);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2);
  K_2=-t*(1-x)*(1-yy-TMath::Power(yy*eps/2.,2))*(1-t_min/t)*(TMath::Sqrt(1+eps_2)+(4*x*(1-x)+eps_2)/(4*(1-x))*(t-t_min)/q2)/q2;
  K=TMath::Sqrt(K_2);
  prefactor=hbar_c*TMath::Power(alpha,3)*x/(TMath::Pi()*8*q2*TMath::Sqrt(1+eps_2))/q2/yy;
  Jacobian=1/(2*x*M*kz);
  prefactor=Jacobian*prefactor;

  coeff[0]=q2*(q2+x*t)/TMath::Power((2-x)*q2+x*t,2)*prefactor*2*(2-2*yy+TMath::Power(yy,2)+eps_2*TMath::Power(yy,2)/2.)/(1+eps_2); //n=0

  coeff[1]=prefactor*8*K/(2-x)/(1+eps_2)*(2-yy)*TMath::Cos(phi); //n=1 cos

  coeff[2]=-prefactor*8*hel*K*yy/(2-x)/(1+eps_2)*TMath::Sqrt(1+eps_2)*TMath::Sin(phi); //n=1 sin

  coeff[3]=prefactor*16*q2*K_2/M_2/TMath::Power(2-x,2)/(1+eps_2)*TMath::Cos(2*phi); //n=2 cos
  return coeff;
}

TArrayD BH_DVCS_interference(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi, Double_t hel){
  yy=q2/2./M/x/kz;
  phi=TMath::Pi()-phi; // conversion phi_trento to phi_BMK
  TArrayD coeff(19);
  TArrayD coeff_CFF(6);
  x_2=TMath::Power(x,2);
  eps=2*x*M/TMath::Sqrt(q2);
  eps_2=TMath::Power(eps,2);
  t_min=-q2*(2*(1-x)*(1-TMath::Sqrt(1+eps_2))+eps_2)/(4*x*(1-x)+eps_2); 
  K_2=-t*(1-x)*(1-yy-TMath::Power(yy*eps/2.,2))*(1-t_min/t)*(TMath::Sqrt(1+eps_2)+(4*x*(1-x)+eps_2)/(4*(1-x))*(t-t_min)/q2)/q2;
  K=TMath::Sqrt(K_2); 
  J=(1-yy-yy*eps_2/2.)*(1+t/q2)-(1-x)*(2-yy)*t/q2;
  P1=-(J+2*K*TMath::Cos(phi))/(yy*(1+eps_2)); 
  P2=1+t/q2+(J+2*K*TMath::Cos(phi))/(yy*(1+eps_2));
  prefactor=hbar_c*TMath::Power(alpha,3)*x*yy/(TMath::Pi()*8*q2*TMath::Sqrt(1+eps_2))/x/TMath::Power(yy,3)/t/P1/P2;
  Jacobian=1/(2*x*M*kz);
  prefactor=Jacobian*prefactor;

  Kp_2=q2*K_2/(1-yy+TMath::Power(eps*yy/2.,2));
  Kp=TMath::Sqrt(Kp_2);
  Double_t p_eff=TMath::Sqrt(2)/(2-x)*Kp/TMath::Sqrt(q2); 

  //Following coefficient are just for H

  coeff[0]=-4*prefactor*(2-yy)*(1+TMath::Sqrt(1+eps_2))/TMath::Power(1+eps_2,2)*(Kp_2/q2*TMath::Power(2-yy,2)/TMath::Sqrt(1+eps_2)+t/q2*(1-yy-TMath::Power(eps*yy/2.,2))*(2-x)*(1+(2*x*(2-x+(TMath::Sqrt(1+eps_2)-1)/2.+eps_2/2./x)*t/q2+eps_2)/((2-x)*(1+TMath::Sqrt(1+eps_2))))); //n=0

  coeff[1]=-16*K*(1-yy-TMath::Power(eps*yy/2.,2))/TMath::Power(1+eps_2,2.5)*(x*t/q2*(1+(1-x)*(TMath::Sqrt(eps_2+1)-1)/2./x+eps_2/4./x)-3*eps_2/4.)-4*K*(2-2*yy+yy*yy+eps_2*yy*yy/2.)*(1+TMath::Sqrt(1+eps_2)-eps_2)/TMath::Power(1+eps_2,2.5)*(1-(1-3*x)*t/q2+x*t/q2*(1-TMath::Sqrt(1+eps_2)+3*eps_2)/(1+TMath::Sqrt(1+eps_2)-eps_2)); coeff[1]=coeff[1]*prefactor*TMath::Cos(phi);// n=1 cos

  coeff[2]=8*hel*K*(2-yy)*yy/(1+eps_2)*(1+(1-x+(TMath::Sqrt(1+eps_2)-1)/2.)/(1+eps_2)*(t-t_min)/q2); coeff[2]=coeff[2]*prefactor*TMath::Sin(phi);// n=1 sin 

  coeff[3]=8*(2-yy)*(1-yy-TMath::Power(eps*yy/2.,2))/TMath::Power(1+eps_2,2)*(Kp_2/q2*2*eps_2/(1+eps_2+TMath::Sqrt(1+eps_2))+x*t*(t-t_min)/q2/q2*(1-x-(TMath::Sqrt(1+eps_2)-1)/2.+eps_2/2./x)); coeff[3]=coeff[3]*prefactor*TMath::Cos(2*phi);// n=2 cos

  coeff[4]=4*hel*yy*(1-yy-TMath::Power(eps*yy/2.,2))/TMath::Power(1+eps_2,1.5)*(1+TMath::Sqrt(1+eps_2)-2*x)*(t-t_min)/q2*((eps_2-x*(TMath::Sqrt(1+eps_2)-1))/(1+TMath::Sqrt(eps_2+1)-2*x)-(2*x+eps_2)/(2*TMath::Sqrt(1+eps_2))*(t-t_min)/q2);    coeff[4]=coeff[4]*prefactor*TMath::Sin(2*phi);// n=2 sin

  coeff[5]=-8*K*(1-yy-TMath::Power(eps*yy/2.,2))*(TMath::Sqrt(1+eps_2)-1)/TMath::Power(1+eps_2,2.5)*(t/q2*(1-x)+(TMath::Sqrt(1+eps_2)-1)/2.*(1+t/q2)); coeff[5]=coeff[5]*prefactor*TMath::Cos(3*phi);// n=3 cos

  //Following ones are for H3_eff

  coeff[6]=prefactor*p_eff*12*TMath::Sqrt(2)*K*(2-yy)*TMath::Sqrt(1-yy-TMath::Power(yy*eps/2.,2))*(eps_2+(2-6*x-eps_2)*t/3./q2)/TMath::Power(1+eps_2,2.5); //n=0

  coeff[7]=prefactor*p_eff*8*TMath::Sqrt(2)*TMath::Sqrt(1-yy-TMath::Power(yy*eps/2.,2))/TMath::Power(1+eps_2,2)*(TMath::Power(2-yy,2)*(t-t_min)/q2*(1-x+((1-x)*x+eps_2/4.)*(t-t_min)/TMath::Sqrt(1+eps_2)/q2)+(1-yy-TMath::Power(yy*eps/2.,2))/TMath::Sqrt(1+eps_2)*(1-(1-2*x)*t/q2)*(eps_2-2*(1+eps_2/2./x)*x*t/q2))*TMath::Cos(phi); //n=1 cos

  coeff[8]=-prefactor*p_eff*8*TMath::Sqrt(2)*K*(2-yy)*TMath::Sqrt(1-yy-TMath::Power(yy*eps/2.,2))*(1+eps_2/2.)/TMath::Power(1+eps_2,2.5)*(1+(1+eps_2/2./x)*x*t/(1+eps_2/2.)/q2)*TMath::Cos(2*phi); //n=2 cos

  coeff[9]=prefactor*hel*p_eff*8*TMath::Sqrt(2)*(2-yy)*yy*TMath::Sqrt(1-yy-TMath::Power(yy*eps/2.,2))/TMath::Power(1+eps_2,2)*Kp_2/q2*TMath::Sin(phi); //n=1 sin

  coeff[10]=prefactor*hel*p_eff*8*TMath::Sqrt(2)*K*yy*TMath::Sqrt(1-yy-TMath::Power(yy*eps/2.,2))*(1+eps_2/2.)/TMath::Power(1+eps_2,2)*(1+(1+eps_2/2./x)*x*t/(1+eps_2/2.)/q2)*TMath::Sin(2*phi); //n=2 sin

  //Here I have implemented two CFF according to hotfix which give delta C when Q goes to infinity
  //C^V
  coeff[11]=prefactor*8*(2-yy)*x*t/TMath::Power(1+eps_2,2)/q2*(TMath::Power(2-yy,2)*Kp_2/q2/TMath::Sqrt(1+eps_2)+(1-yy-yy*yy*eps_2/4.)*(1+TMath::Sqrt(1+eps_2))*(1+t/q2)*(1+(TMath::Sqrt(1+eps_2)-1+2*x)*t/q2/(1+TMath::Sqrt(1+eps_2)))/2.); //n=0

  coeff[12]=prefactor*16*K*x*t/q2/TMath::Sqrt(TMath::Power(1+eps_2,5))*(TMath::Power(2-yy,2)*(1-(1-2*x)*t/q2)+(1-yy-yy*yy*eps_2/4.)*(1-2*x+TMath::Sqrt(1+eps_2))*(t-t_min)/q2/2.)*TMath::Cos(phi); //n=1

  coeff[13]=prefactor*8*(2-yy)*(1-yy-yy*yy*eps_2/4.)*x*t/q2/TMath::Power(1+eps_2,2)*(4*Kp_2/q2/TMath::Sqrt(1+eps_2)+(1+TMath::Sqrt(1+eps_2)-2*x)*(1+t/q2)*(t-t_min)/q2/2.)*TMath::Cos(2*phi); //n=2

  coeff[14]=-prefactor*8*K*(1-yy-yy*yy*eps_2/4.)*x*t/q2/TMath::Sqrt(TMath::Power(1+eps_2,5))*(TMath::Sqrt(1+eps_2)-1+(1+TMath::Sqrt(1+eps_2)-2*x)*t/q2)*TMath::Cos(3*phi); //n=3

  //C^A
  coeff[15]=prefactor*8*(2-yy)*t/q2/TMath::Power(1+eps_2,2)*(TMath::Power(2-yy,2)*Kp_2*(1+TMath::Sqrt(1+eps_2)-2*x)/q2/TMath::Sqrt(1+eps_2)/2.+(1-yy-yy*yy*eps_2/4.)*((1+TMath::Sqrt(1+eps_2))*(1+TMath::Sqrt(1+eps_2)-x+t*(TMath::Sqrt(1+eps_2)-1+x*(3+TMath::Sqrt(1+eps_2)-2*x)/(1+TMath::Sqrt(1+eps_2)))/q2)/2.-2*Kp_2/q2)); //n=0
  
  coeff[16]=-16*prefactor*K*t/q2/TMath::Power(1+eps_2,2)*((1-yy-yy*yy*eps_2/4.)*(1-(1-2*x)*t/q2+(4*x*(1-x)+eps_2)*(t-t_min)/(4.*TMath::Sqrt(1+eps_2))/q2)-TMath::Power(2-yy,2)*(1-x/2.+(1-t/q2)*(1+TMath::Sqrt(1+eps_2)-2*x)/4.+(t-t_min)*(4*x*(1-x)+eps_2)/q2/2./TMath::Sqrt(1+eps_2)))*TMath::Cos(phi); //n=1

  coeff[17]=prefactor*4*(2-yy)*(1-yy-yy*yy*eps_2/4.)*t/q2/TMath::Power(1+eps_2,2)*(4*Kp_2*(1-2*x)/q2/TMath::Sqrt(1+eps_2)-(3-TMath::Sqrt(1+eps_2)-2*x+eps_2/x)*x*(t-t_min)/q2)*TMath::Cos(2*phi); //n=2

  coeff[18]=16*prefactor*K*(1-yy-yy*yy*eps_2/4.)*t*(t-t_min)/q2/q2/TMath::Power(TMath::Sqrt(1+eps_2),5)*(x*(1-x)+eps_2/4.)*TMath::Cos(3*phi); //n=3

  coeff_CFF[0]=coeff[0]+coeff[1]+coeff[3]+coeff[5]; //Re H
  coeff_CFF[1]=coeff[2]+coeff[4]; //Im H
  coeff_CFF[2]=coeff[6]+coeff[7]+coeff[8]; //Re H3_eff
  coeff_CFF[3]=coeff[9]+coeff[10]; //Im H3_eff
  coeff_CFF[4]=coeff[11]+coeff[12]+coeff[13]+coeff[14]; //Re C^V
  coeff_CFF[5]=coeff[15]+coeff[16]+coeff[17]+coeff[18]; //Re C^A
  return coeff_CFF;
}
