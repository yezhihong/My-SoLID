#include <TROOT.h>
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "TFile.h"
#include "TTree.h"
#include "SIDIS.h"

int main(Int_t argc, char *argv[]){

	gRandom->SetSeed(0);// uses time for seed
		/*Inputs&Output{{{*/
		//initialize

		Double_t momentum_ele = 0.0;cerr<<"Incoming Electron Momentum (GeV/c)="; cin>>momentum_ele;
		Double_t momentum_ion = 0.0;cerr<<"Incoming Ion Momentum (GeV/c)="; cin>>momentum_ion;
		Int_t target_flag = 0;cerr<<"Which Target? (1->Proton,2->Deutron,3->He3) = "; cin>> target_flag;
		Int_t particle_flag = 0;cerr<<"Which Hadron? (1->pi+,-1->pi-,2->k+,-2->k-) = ";cin>>particle_flag;

		//define output file
		if (target_flag == 2){
			momentum_ion = momentum_ion/2.;
			/*}}}*/

		/*Define{{{*/
		double Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon, jacoF;
		Double_t theta_gen= 0.0 , phi_gen = 0.0, mom_gen = 0.0;

		Double_t mom_gen_ele,mom_gen_had;
		Double_t theta_gen_ele,theta_gen_had;
		Double_t phi_gen_ele,phi_gen_had;
		Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
		Double_t dxs_hm,dxs_hp,dilute_hp,dilute_hm,weight_hp,weight_hm;
		Int_t nsim = 0;
		/*}}}*/

		//Only initialize once here
		SIDIS *sidis = new SIDIS();
		sidis->SetLAPDF();
		//sidis->Print();


		/*Generator{{{*/
		//For electron
		phi_gen = gRandom->Uniform(0.,2.*PI);//Rad
		theta_gen = acos(gRandom->Uniform(cos(30/DEG),cos(7/DEG)));//Rad
		mom_gen = gRandom->Uniform(0.5,momentum_ele);//GeV/c
		mom_gen_ele = mom_gen; theta_gen_ele = theta_gen*DEG; phi_gen_ele = phi_gen*DEG;

		//For hadron
		theta_gen = acos(gRandom->Uniform(cos(30/DEG),cos(7/DEG)));//Rad
		phi_gen = gRandom->Uniform(0.,2.*PI); //Rad
		mom_gen = gRandom->Uniform(0.5,6);//GeV/c
		mom_gen_had = mom_gen; theta_gen_had = theta_gen*DEG; phi_gen_had = phi_gen*DEG;
		/*}}}*/

		//Initial Kinematics Settings
		sidis->Init(momentum_ele, momentum_ion,          //GeV/c,  GeV/c 
				mom_gen_ele, theta_gen_ele, phi_gen_ele, //GeV/c, Deg, Deg
				mom_gen_had, theta_gen_had, phi_gen_had, //GeV/c, Deg, Deg
				target_flag, particle_flag);             //1.or.2.or.3, 1.or.-1.or.2.or.-2

		//Return Physics Quantities
		x = sidis->x; z = sidis->z; y = sidis->y; Q2 = sidis->Q2; W=sidis->W; Wp = sidis->Wp;
		s = sidis->s; nu = sidis->nu; pt = sidis->pt; jacoF = sidis->jacoF;
		gamma=sidis->gamma; epsilon=sidis->epsilon;

		theta_q = sidis->theta_q;  theta_s = sidis->theta_s; 
		phi_h = sidis->phi_h; phi_s= sidis->phi_s;
		mom_ele = sidis->mom_ele; theta_ele = sidis->theta_ele; phi_ele = sidis->phi_ele;
		mom_had = sidis->mom_had; theta_had = sidis->theta_had; phi_ele = sidis->phi_had;

		/*Get XS in DIS Region{{{*/
		if( Q2 >=1.0&&Q2<=10.&&W>= 2.3&&Wp>=1.6&&z>0.3&&z<0.7){

			sidis->CalcXS();//Start to calculate XS
			dxs_hp = sidis->GetXS_HP(); //Get the hadron+ XS
			dxs_hm = sidis->GetXS_HM(); //Get the hardron- XS
			dilute_hp = sidis->GetDilute_HP(); //Dilution Factors
			dilute_hm = sidis->GetDilute_HM(); //Dilution Factors


			if ((dxs_hp+dxs_hm)!=0){
				count++;
				cerr << count << "\r" ;
				cerr<<Form("N=%d, P_e=%f, Th_e=%f, P_h=%f, Th_h=%f, Q2=%f, W=%f, Wp=%f, z=%f, XS+=%f, XS-=%f",
						count,mom_gen_ele, theta_gen_ele, mom_gen_had, theta_gen_had,
						Q2, W, Wp,z, dxs_hp,dxs_hm)<<"\r";

			}
		}
		/*}}}*/

		//Free memory
		delete sidis;

		return 0;
} 

