#ifndef BMHOTFIX_H
#define BMHOTFIX_H

#include "TArrayD.h"
#include "TMath.h"

Double_t alpha=7.297352533e-3;
Double_t hbar_c=0.38937966e6;
Double_t e=1.6e-19;
Double_t Ebeam=5.75;
Double_t M=0.938271998;
Double_t M_2=M*M;

Double_t prefactor;Double_t Jacobian;
Double_t K;Double_t K_2;
Double_t Kp;Double_t Kp_2;
Double_t J;
Double_t yy;
Double_t P1; Double_t P2;
Double_t t_min;
Double_t F1_2;
Double_t F2_2;
Double_t F_sum_2;

Double_t x_2;
Double_t eps_2;
Double_t eps;

Double_t GMp(Double_t q2);
Double_t GEp(Double_t q2);

Double_t F1p(Double_t q2);
Double_t F2p(Double_t q2);

Double_t KellyE(Double_t q2);
Double_t KellyM(Double_t q2);

Double_t BH(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi);

TArrayD DVCS(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi, Double_t hel);

TArrayD BH_DVCS_interference(Double_t q2, Double_t x, Double_t t, Double_t kz, Double_t phi, Double_t hel);
 
#endif
