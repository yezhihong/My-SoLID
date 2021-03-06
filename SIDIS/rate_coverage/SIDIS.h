/////////////////////////////////////////////////////////////////////////////
//      Semi-Inclusive Deep Inelstic Scattering Cross Section Model        //
// Note:                                                                   //
//   This model is developed by Xin in his "collider" code. I extracted    //
//   the cross section parts and coded it into a C++ class which can be    //
//   easily embeded by other programs.                                     //
//   --- Zhihong Ye, 06/09/2014                                            //
/////////////////////////////////////////////////////////////////////////////
/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
/*}}}*/
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
#include "LHAPDF/LHAPDF.h"
//#include <TMatrix.h>
/*}}}*/

using namespace std;
using namespace LHAPDF;

//const TString LHAPDF_Dir  = "/work/halla/solid/yez/evgen/sidis_model/lhapdf";
char *LHAPDF_Dir = std::getenv("LHAPDF");
const double DEG=180./3.1415926;
const double PI=3.1415926;

class SIDIS
{
	public:
		SIDIS(){}
		virtual ~SIDIS(){}
		SIDIS(SIDIS const&);
		SIDIS& operator=(SIDIS const&){};

		/*void Print(){{{*/
        void Print(){
	       cerr<<"===================================================================================================="<<endl;	
		   cerr<<" --- Useage:"<<endl;
		   cerr<<"      1, Init(energy_ele, energy_ion, P_ele, Th_ele, Ph_ele, P_had, Th_had, Ph_had, tgt_flag, ptc_flg)"<<endl;
		   cerr<<"           where energy_* and P_* are in GeV, and Theta_* and Phi_* are in Deg,"<<endl;
		   cerr<<"                 tgt_flag = 1->Proton, 2->Deutron, 3->He3,"<<endl;
		   cerr<<"                  pct_flag = 1->Pion+, 2->Pion-, 3->Kaon+, 4->Kaon-;"<<endl;
		   cerr<<"      2, CalcXS() --> Calculated dxs_hp and dxs_hm; "<<endl; 
		   cerr<<"      3, GetXS_HP() --> Get SIDIS XS of Pion+ or Kaon+ ,"<<endl;
		   cerr<<"         GetXS_HM() --> Get SIDIS XS of Pion- or Kaon- ,"<<endl;
		   cerr<<"         GetDilute_HP() --> Get Dilution of Pion+ or Kaon+ ,"<<endl;
		   cerr<<"         GetDilute_HM() --> Get Dilution of Pion- or Kaon- ,"<<endl;
		   cerr<<"      4, Other Output Quantities can be calculated after running CalcXS: "<<endl;
		   cerr<<"           Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon, jacoF "<<endl;
	       cerr<<"===================================================================================================="<<endl;	
		
		}
		/*}}}*/

		/*Return Values{{{*/
		double GetXS_HP(){
		  return XS_HP;
		}
		double GetXS_HM(){
			return XS_HM;
		}
		double GetDilute_HP(){
		  return dilute_hp;
		}
		double GetDilute_HM(){
			return dilute_hm;
		}
        /*}}}*/

        /*SetLAPDF{{{*/
		void SetLAPDF(){
			const int SUBSET = 1;
			const string NAME = "cteq6m";

			TString LHAPDF_path=Form("%s/share/lhapdf",LHAPDF_Dir);
			setPDFPath(LHAPDF_path.Data());
			LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
		}
        /*}}}*/
		
		/*void Init( double energy_ele, double energy_ion,...){{{*/
		void Init(double energy_ele,double energy_ion,                        // GeV          GeV   
				double p_ele, double th_ele, double ph_ele,         // GeV/c         DEG            DEG 
				double p_had, double th_had, double ph_had,         // GeV/c         DEG            DEG 
				int tgt_flag, int ptcl_flag){

			/*Define{{{*/
			double mom_pro;
			double energy_pro;
			double mom_ini_ele; 
			double energy_ini_ele;

			target_flag = tgt_flag;
			particle_flag = ptcl_flag;
			double mass_e = 0.511e-3; // electron mass in GeV
			double mass_p = 0.93827;
			double mass_pi = 0.13957;
			//double mass_pi0 = 0.1349766;
			double mass_kaon = 0.493667;
			double mass_n = 0.939566;

			//initialize
			double mass_target= 0.0;;
			if (target_flag ==1){
				mass_target = mass_p;
			}else if (target_flag == 2){
				mass_target = (mass_n+mass_p)/2.;
				energy_ion = energy_ion/2.;
			}else if (target_flag == 3){
				mass_target = (mass_n+2.*mass_p)/3.;
				energy_ion = energy_ion*2./3.;
			}else{
				cout << "target_flag is wrong! 1,2,3 " << endl;
			}

			if (particle_flag == 1){
				mass_hadron = mass_pi;
			}else if (particle_flag == -1){
				mass_hadron = mass_pi;
			}else if (particle_flag == 2){
				mass_hadron = mass_kaon;
			}else if (particle_flag == -2){
				mass_hadron = mass_kaon;
			}else{
				cout << "particle_flag is wrong +-1 and +-2" << endl;
			}	

			// //define the 4-momentum
			//define electron direction as +z assuming a proton/neutron for the ion mass now 
			// approximation for SIDIS
			TLorentzVector *P4_ini_ele = new TLorentzVector(0.,0.,energy_ele,sqrt(energy_ele*energy_ele + mass_e*mass_e));
			TLorentzVector *P4_ini_ion = new TLorentzVector(0.,0.,-energy_ion,sqrt(energy_ion*energy_ion + mass_target*mass_target));
			TLorentzVector *P4_fin_had = new TLorentzVector(0.,0.,0.,1.);
			TLorentzVector *P4_fin_ele = new TLorentzVector(0.,0.,0.,1.);
			TLorentzVector *P4_q = new TLorentzVector(0.,0.,0.,1.);

			TLorentzVector *lrz_P4_q = new TLorentzVector(0.,0.,0.,1.);
			TLorentzVector *lrz_P4_h = new TLorentzVector(0.,0.,0.,1.);
			TLorentzVector *lrz_P4_ef = new TLorentzVector(0.,0.,0.,1.);
			//TLorentzVector *lrz_P4_ion = new TLorentzVector(0.,0.,0.,1.);
			//TLorentzVector *lrz_P4_ei = new TLorentzVector(0.,0.,0.,1.);

			double vn=energy_ion/sqrt(energy_ion*energy_ion + mass_target*mass_target);
			TVector3 vnboost(0.,0.,vn);

			TVector3 p3_q,p3_fin_ele,p3_fin_had;
			TVector3 p3_target_spin(0.,1.,0);
			mom_ini_ele = fabs(energy_ele);
			energy_ini_ele = sqrt(energy_ele*energy_ele + mass_e*mass_e);

			mom_pro = fabs(energy_ion);
			energy_pro = sqrt(energy_ion*energy_ion + mass_target*mass_target);

			/*}}}*/

			/*Kinematics Quantities{{{*/
			//For electron
            px_ele = p_ele*sin(th_ele/180.*3.1415926)*cos(ph_ele/180.*3.1415926);
			py_ele = p_ele*sin(th_ele/180.*3.1415926)*sin(ph_ele/180.*3.1415926);
			pz_ele = p_ele*cos(th_ele/180.*3.1415926);
			E0_ele = sqrt(p_ele*p_ele+mass_e*mass_e);
			P4_fin_ele->SetPxPyPzE(p_ele*sin(th_ele/180.*3.1415926)*cos(ph_ele/180.*3.1415926),
					p_ele*sin(th_ele/180.*3.1415926)*sin(ph_ele/180.*3.1415926),
					p_ele*cos(th_ele/180.*3.1415926)
					,sqrt(p_ele*p_ele+mass_e*mass_e));

			//For hadron
			px_had = p_had*sin(th_had/180.*3.1415926)*cos(ph_had/180.*3.1415926);
		    py_had = p_had*sin(th_had/180.*3.1415926)*sin(ph_had/180.*3.1415926);
			pz_had = p_had*cos(th_had/180.*3.1415926);
			E0_had = sqrt(p_had*p_had+mass_hadron*mass_hadron);
			P4_fin_had->SetPxPyPzE(p_had*sin(th_had/180.*3.1415926)*cos(ph_had/180.*3.1415926),
					p_had*sin(th_had/180.*3.1415926)*sin(ph_had/180.*3.1415926),
					p_had*cos(th_had/180.*3.1415926)
					,sqrt(p_had*p_had+mass_hadron*mass_hadron));


			*P4_q = *P4_ini_ele - *P4_fin_ele;
			Q2 = - (*P4_q)*(*P4_q);
			W = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele);
			Wp = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had);

			s = (*P4_ini_ele + *P4_ini_ion )*(*P4_ini_ele + *P4_ini_ion);
			nu = (*P4_ini_ion) * (*P4_q)/mass_target;
			x = Q2/ 2. / mass_target / nu;
			z = ((*P4_ini_ion)*(*P4_fin_had))/((*P4_ini_ion)*(*P4_q));
			y = ((*P4_ini_ion)*(*P4_q))/((*P4_ini_ion)*(*P4_ini_ele));
			W = sqrt(W);
			Wp = sqrt(Wp);
				
			*lrz_P4_q = *P4_q;
			*lrz_P4_h = *P4_fin_had;
			*lrz_P4_ef = *P4_fin_ele;

			lrz_P4_ef->Boost(vnboost);
			lrz_P4_h->Boost(vnboost); 
			lrz_P4_q->Boost(vnboost);

			gamma = 2*mass_target *x/sqrt(Q2);
			epsilon = (1-y-0.25*gamma*gamma*y*y)/(1-y+0.5*y*y+0.25*gamma*gamma*y*y);

			for(Int_t j=0;j<3;j++){
				p3_q(j)=(*lrz_P4_q)(j);
				p3_fin_ele(j)=(*lrz_P4_ef)(j);
				p3_fin_had(j)=(*lrz_P4_h)(j);
			}

			pt = p3_fin_had.Perp(p3_q);

			theta_q =p3_q.Theta();
			phi_q = p3_q.Phi();
			p3_target_spin.SetXYZ(0.,1.,0.);

			p3_fin_ele.RotateZ(-phi_q);
			p3_fin_ele.RotateY(-theta_q);
			p3_fin_had.RotateZ(-phi_q);
			p3_fin_had.RotateY(-theta_q);
			p3_target_spin.RotateZ(-phi_q);
			p3_target_spin.RotateY(-theta_q);

			phi_s = (Azimuthalphi(p3_target_spin(0),p3_target_spin(1))-Azimuthalphi(p3_fin_ele(0),p3_fin_ele(1)));
			phi_h = (Azimuthalphi(p3_fin_had(0),p3_fin_had(1))-Azimuthalphi(p3_fin_ele(0),p3_fin_ele(1)));
			theta_s = p3_target_spin.Theta();

			mom_ele = P4_fin_ele->P();
			theta_ele = P4_fin_ele->Theta();
			mom_had = P4_fin_had->P();
			theta_had = P4_fin_had->Theta();
			phi_ele = P4_fin_ele->Phi();
			phi_had = P4_fin_had->Phi();

			jacoF = Jacobian(mom_ele,theta_ele,phi_ele,mom_had,theta_had,phi_had,mom_pro,energy_pro,energy_ele);
			//jacoF = 1.0;

			/*}}}*/
			
			delete P4_ini_ele;delete P4_ini_ion; delete P4_q; delete lrz_P4_q; 
			delete lrz_P4_h; delete lrz_P4_ef;
		} 
		/*}}}*/
		
		/*int CalcXS(){{{*/
		int CalcXS(){
			double bpt_m=4.694;
			double bpt_p=4.661;
			double pt_tmp = 0.0;	
			double dxs_all[3][4];
            double dxs_hp = -1000.0;
            double dxs_hm = -1000.0;

			/*Calculate XS{{{*/
			if (pt<0.8){
				//first method 	 
				bpt_p = 1./(0.2+z*z*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
				bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)

				//second method  use the original value
				// bpt_m=4.694; // HERMES parameterization
				// bpt_p=4.661;
                
				pt_tmp = pt;
				dxs(pt_tmp, bpt_m, bpt_p, &dxs_hp, &dxs_hm,dxs_all[0]);

				dilute_hp = dxs_all[0][0]/dxs_all[0][1];
				dilute_hm = dxs_all[0][2]/dxs_all[0][3];

			}
			else{
				// this part is to generate factor K
				//make sure the DXS is the same at PT= 0.8 GeV
				//calculating the TMD part
				double K[2],dxs_temp[10];
				bpt_p = 1./(0.2+z*z*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
				bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)
			
				//second method  use the original value
				// bpt_m=4.694; // HERMES parameterization
				// bpt_p=4.661;

				pt_tmp = 0.8;
				dxs(pt_tmp, bpt_m, bpt_p, &dxs_temp[0], &dxs_temp[1],dxs_all[0]);	   

				//calculating the TMD parts using the new vpt_p values
				bpt_p = 1./(0.25+z*z*0.28);// <pt^2> = 0.25 GeV^2 (quark internal momentum)
				bpt_m = bpt_p; // <kt^2> = 0.28 GeV^2 (struck quark internal momentum)
		
				//taking into account the NLO etc
				pt_tmp = 1.0;
				dxs(pt_tmp, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
				dxs_hpt(pt_tmp, &dxs_temp[4],&dxs_temp[5],dxs_all[2]);

				K[0] = (dxs_temp[0]-dxs_temp[2])/dxs_temp[4];
				K[1] = (dxs_temp[1]-dxs_temp[3])/dxs_temp[5];

				if (K[0]<0.) K[0] = 0.;
				if (K[1]<0.) K[1] = 0.;

				if (pt>1.2){
					pt_tmp = pt;
					dxs(pt_tmp,bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
					dxs_hpt(pt_tmp,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

					dxs_hp = dxs_temp[2] + K[0]*dxs_temp[4];
					dxs_hm = dxs_temp[3] + K[1]*dxs_temp[5];

					dilute_hp = (dxs_all[1][0] + K[0]*dxs_all[2][0])/(dxs_all[1][1] + K[0]*dxs_all[2][1]);
					dilute_hm = (dxs_all[1][2] + K[1]*dxs_all[2][2])/(dxs_all[1][3] + K[1]*dxs_all[2][3]);

				}else{
					pt_tmp = 1.2;
					dxs(pt_tmp, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
					dxs_hpt(pt_tmp,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

					dxs_temp[2] = dxs_temp[2] + K[0]*dxs_temp[4];
					dxs_temp[3] = dxs_temp[3] + K[1]*dxs_temp[5];

					dxs_temp[4] = (dxs_temp[0]-dxs_temp[2])/(1./0.8/0.8-1./1.2/1.2);
					dxs_temp[5] = dxs_temp[0]-dxs_temp[4]/0.8/0.8;
					dxs_temp[6] = (dxs_temp[1]-dxs_temp[3])/(1./0.8/0.8-1./1.2/1.2);
					dxs_temp[7] = dxs_temp[1]-dxs_temp[6]/0.8/0.8;

					dxs_hp = dxs_temp[4]/pt/pt + dxs_temp[5];
					dxs_hm = dxs_temp[6]/pt/pt + dxs_temp[7];

					dilute_hp = ((dxs_all[1][0] + K[0]*dxs_all[2][0])
							/(dxs_all[1][1] + K[0]*dxs_all[2][1])+dxs_all[0][0]/dxs_all[0][1])/2.;
					dilute_hm = ((dxs_all[1][2] + K[1]*dxs_all[2][2])
							/(dxs_all[1][3] + K[1]*dxs_all[2][3])+dxs_all[0][2]/dxs_all[0][3])/2.;
				}
			}
			/*}}}*/

			/*try to take care of the decay{{{*/
			double decay_part = 0.0;
			if (abs(particle_flag)==1){
				if (theta_had>155./180.*3.1415926){
					decay_part = exp(7./cos(theta_had)*mass_hadron/(2.6*mom_had*3.0));
				}else if (theta_had<=155./180.*3.1415926&&theta_had>=140./180.*3.1415926){
					decay_part = exp(4.5/cos(theta_had)*mass_hadron/(2.6*mom_had*3.0));
				}else if (theta_had <140/180.*3.1415926){
					decay_part = exp(-2.5/sin(theta_had)*mass_hadron/(2.6*mom_had*3.0));
				}
			}else{
				if (theta_had>155./180.*3.1415926){
					decay_part = exp(7./cos(theta_had)*mass_hadron/(1.24*mom_had*3.0));
				}else if (theta_had<=155./180.*3.1415926&&theta_had>=140./180.*3.1415926){
					decay_part = exp(4.5/cos(theta_had)*mass_hadron/(1.24*mom_had*3.0));
				}else if (theta_had <140/180.*3.1415926){
					decay_part = exp(-2.5/sin(theta_had)*mass_hadron/(1.24*mom_had*3.0));
				}
			}
			/*}}}*/
				
			//cerr<<Form("++++ pt=%f, jacoF=%f, dxs_hp=%f,dxs_hm=%f, decay=%f, dilu_hp=%f, dilu_hm=%f", pt, jacoF,dxs_hp,dxs_hm,decay_part,dilute_hp,dilute_hm)<<endl;

			/*Save&Return{{{*/
			dxs_hp = decay_part * dxs_hp; 
			dxs_hm = decay_part * dxs_hm;

			int err_hp = -1, err_hm=-1;
			if (dxs_hp>0.&& dxs_hp<10000.){
				err_hp = 0;
			}
			else{
				dxs_hp = 0;
				err_hp = 1;
			}
	
			XS_HP = dxs_hp;
			XS_HM = dxs_hm;

			if (dxs_hm>0.&& dxs_hm<10000.){
				err_hm = 0;
			}
			else{
				dxs_hm = 0;
				err_hm = 1;
			}

			if(err_hp==0&&err_hm==0) 
				return 0;
			else if(err_hp==0&&err_hm==1)
				return 1;
			else if(err_hp==1&&err_hm==0)
				return -1;
			else
				return 2;
			/*}}}*/


		} 
		/*}}}*/

    private:
		/*double Azimuthalphi(double vx, double vy){{{*/
		double Azimuthalphi(double vx, double vy){
			double pmod, phi, cosf;
			pmod=vx*vx+vy*vy;
			if(pmod>0.) {
				pmod=sqrt(pmod);
				cosf=vx/pmod;
			}
			else{
				phi=-100.;
				cosf=10.;
			}
			if(fabs(cosf)<=1.0) phi=acos(cosf);
			if(vy<0.) phi=2*PI-phi;
			return phi;
		}
		/*}}}*/

		/*double Jacobian(double P_ef,double theta_ef,double phi_ef,...){{{*/
		double Jacobian(double P_ef,double theta_ef,double phi_ef,double P_hf,double theta_hf,double phi_hf,double mom_pro,double energy_pro,double P_ei){


			// need to shift the electron phi angle to 0 in order to do calculation for jacobian

			phi_hf -= phi_ef;
			phi_ef = 0.;

			double qbeta = -1.*mom_pro/energy_pro;
			//double qbeta = 0.;
			double qgamma = 1./sqrt(1-qbeta*qbeta);


			double mass_pi = mass_hadron;
			double P_hx,P_hy,P_hz;  
			P_hx = P_hf*sin(theta_hf)*cos(phi_hf);
			P_hy = P_hf*sin(theta_hf)*sin(phi_hf);
			P_hz = qgamma * (P_hf*cos(theta_hf) - qbeta * sqrt(P_hf*P_hf+mass_pi*mass_pi));



			//P_hz = P_hf*cos(theta_hf);

			double qz = qgamma*(P_ei-P_ef*cos(theta_ef)-qbeta*(P_ei-P_ef));

			double alpha,beta,gamma1;
			alpha = atan(P_ef*sin(theta_ef)*cos(phi_ef)/qz);
			beta = atan(-P_ef*sin(theta_ef)*sin(phi_ef)*
					qz/cos(alpha)/
					(pow(P_ef*sin(theta_ef)*sin(phi_ef),2) + 
					 pow(qz,2)));
			gamma1 = atan(sin(beta)/tan(alpha));

			double a2,b2,c2,d2,e2,f2;
			a2 = (-sin(gamma1)*cos(alpha)+cos(gamma1)*sin(beta)*sin(alpha));
			b2 = cos(gamma1)*cos(beta);
			c2 = -(sin(gamma1)*sin(alpha)+cos(gamma1)*sin(beta)*cos(alpha));
			d2 = cos(gamma1)*cos(alpha) + sin(gamma1)*sin(beta)*sin(alpha);
			e2 = sin(gamma1)*cos(beta);
			f2 = cos(gamma1)*sin(alpha) - sin(gamma1)*sin(beta)*cos(alpha);

			double a1,b1,c1,d1,e1,f1;
			a1 = cos(alpha)*cos(alpha) + sin(beta)*sin(beta) * sin(alpha)*sin(alpha);
			b1 = cos(beta)* cos(beta);
			c1 = sin(alpha)*sin(alpha) + sin(beta)*sin(beta) * cos(alpha)*cos(alpha);      
			d1 = 2*cos(beta)*sin(beta)*sin(alpha);
			e1 = -2*cos(beta)*sin(beta)*cos(alpha);
			f1 = 2*cos(beta)*cos(beta)*sin(alpha)*cos(alpha);

			double Px1,Py1,Pz1,Px2,Py2,Pz2;
			Px1 = P_hf*sin(theta_hf)*(-sin(phi_hf));
			Py1 = P_hf*sin(theta_hf)*cos(phi_hf);
			Pz1 = 0.;

			Px2 = -P_hf*cos(phi_hf)/tan(theta_hf);
			Py2 = -P_hf*sin(phi_hf)/tan(theta_hf);
			Pz2 = P_hf*qgamma;

			double Px3 = sin(theta_hf)*cos(phi_hf);
			double Py3 = sin(theta_hf)*sin(phi_hf);
			double Pz3 = qgamma*(cos(theta_hf) -qbeta/sqrt(P_hf*P_hf+mass_pi*mass_pi)*P_hf);

			double P_hyp,P_hxp;
			P_hyp = a2*P_hx + b2*P_hy + c2*P_hz;
			P_hxp = d2*P_hx + e2*P_hy + f2*P_hz;


			//cout << a1*P_hx*P_hx + b1*P_hy*P_hy + c1*P_hz*P_hz + d1*P_hx*P_hy + e1 *P_hy*P_hz + f1 * P_hx*P_hz << "\t" << atan(P_hyp/P_hxp) << endl;

			double jaco_F;
			double jaco2; // py/pa
			double jaco1; // px/pb
			double jaco3; // pz/pd
			double jaco4; // pl/pc
			double jaco5; // pm/pe
			double jaco6; // pm/pf
			double jaco7; // pn/pe
			double jaco8; // pn/pf

			double jaco9;  // px/pa
			double jaco10; // py/pb
			double jaco11; // pz/pe
			double jaco12; // pm/pd
			double jaco13; // pn/pd

			jaco9 = (P_ei*(1.-cos(theta_ef))*(P_ei*(energy_pro+mom_pro)-P_ef*(energy_pro+mom_pro*cos(theta_ef)))+P_ef*P_ei*(1-cos(theta_ef))*(energy_pro+mom_pro*cos(theta_ef)))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef)))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // px / pa

			jaco1 = (-P_ef*P_ei*(P_ei*(mom_pro+energy_pro)-P_ef*(energy_pro+mom_pro*cos(theta_ef)))+P_ei*P_ef*(1-cos(theta_ef))*P_ef*mom_pro)/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef)))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // p x / pb



			jaco2 = - (energy_pro+mom_pro*cos(theta_ef))/(P_ei*(energy_pro+mom_pro)); // py/pa

			jaco10 = (-1*P_ef*mom_pro)/(P_ei*(energy_pro+mom_pro)); // py / pb

			jaco3 = (P_hf*energy_pro/sqrt(P_hf*P_hf + mass_pi*mass_pi) + mom_pro*cos(theta_hf))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // p z / p d
			jaco11 = mom_pro*P_hf/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // p z/ p e

			jaco4 = -1;

			jaco5 = (2*a1*P_hx+d1*P_hy+f1*P_hz)*Px1
				+ (2*b1*P_hy+d1*P_hx+e1*P_hz)*Py1; // p pt^2/ p (phi_pi)   or p m / pf

			jaco6 = (2*a1*P_hx+d1*P_hy+f1*P_hz)*Px2
				+ (2*b1*P_hy+d1*P_hx+e1*P_hz)*Py2
				+ (2*c1*P_hz+e1*P_hy+f1*P_hx)*Pz2; // ppt^2/ p (cos(theta_pi)) or p m / p e

			jaco7 = (P_hxp*(a2*Px1 + b2*Py1) - 
					P_hyp*(d2*Px1 + e2*Py1) )/
				(P_hxp*P_hxp + P_hyp*P_hyp); // p n / pf

			jaco8 = (P_hxp*(a2*Px2 + b2*Py2 + c2*Pz2) - 
					P_hyp*(d2*Px2 + e2*Py2 + f2*Pz2) )/
				(P_hxp*P_hxp + P_hyp*P_hyp); // p n / p e



			jaco12 = (2*a1*P_hx+d1*P_hy+f1*P_hz)*Px3
				+ (2*b1*P_hy+d1*P_hx+e1*P_hz)*Py3
				+ (2*c1*P_hz+e1*P_hy+f1*P_hx)*Pz3; // or p m / p d
			jaco13 = (P_hxp*(a2*Px3 + b2*Py3 + c2*Pz3) - 
					P_hyp*(d2*Px3 + e2*Py3 + f2*Pz3) )/
				(P_hxp*P_hxp + P_hyp*P_hyp); // p n / p d

			//cout << jaco13 << endl;

			//(pz/pd *(p m / p e * p n / pf - p m / pf * // p n / p e) - pz/pe*(p m / p d * p n / pf - p n / p d * p m / pf))

			jaco_F = jaco4*(jaco2*jaco1-jaco9*jaco10)*(jaco3*(jaco6*jaco7-jaco5*jaco8)-jaco11*(jaco12*jaco7-jaco13*jaco5));
			jaco_F = fabs(jaco_F);

			//cout << jaco_F << "\t" << jaco1 << "\t" << jaco2 << "\t" << jaco3 << "\t" << jaco5 << "\t" << jaco6 << "\t" << jaco7 << "\t" << jaco8 << endl;

			return jaco_F;
		}
		/*}}}*/

		/*void dxs(double x, double y, double Q2, ...){{{*/
		void dxs(double pt_tmp,double bpt_m, double bpt_p, double* dxs_hp, double* dxs_hm,double* dxs_all){
			double qu=2./3.;
			double qd=-1./3.;
			double qs=-1./3.;

			*dxs_hp = 1./137.035/137.035/x/y/Q2/1000000. *y*y/(2*(1-epsilon))*(1+gamma*gamma/2./x)*(197.3*197.3/100.*1.0e9);
			*dxs_hm = *dxs_hp;

			*dxs_hp = (*dxs_hp)*bpt_p/PI*exp(-bpt_p*pt_tmp*pt_tmp)*x;
			*dxs_hm = (*dxs_hm)*bpt_m/PI*exp(-bpt_m*pt_tmp*pt_tmp)*x;

			*dxs_hp = *dxs_hp*jacoF;
			*dxs_hm = *dxs_hm*jacoF;

			double df_p_hp=0,df_p_hm=0;
			double df_n_hp=0,df_n_hm=0;

			double uquark,dquark,squark,ubarquark,dbarquark,sbarquark;
			double D_fav,D_unfav,D_s,D_g;
			if (fabs(particle_flag)==1) {
				//pion
				// if (target_flag>=1){
				// hydrogen
				un_ff(1,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);
				uquark = un_pdf(1,x,Q2);
				dquark = un_pdf(2,x,Q2);
				squark = un_pdf(3,x,Q2);
				ubarquark = un_pdf(-1,x,Q2);
				dbarquark = un_pdf(-2,x,Q2);
				sbarquark = un_pdf(-3,x,Q2);

				df_p_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
				df_p_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;


				// }

				//if (target_flag>=2){
				//deuterium add an neutron
				un_ff(1,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);
				uquark = un_pdf(2,x,Q2);
				dquark = un_pdf(1,x,Q2);
				squark = un_pdf(3,x,Q2);
				ubarquark = un_pdf(-2,x,Q2);
				dbarquark = un_pdf(-1,x,Q2);
				sbarquark = un_pdf(-3,x,Q2);

				df_n_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
				df_n_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;



				// }
				//     if (target_flag>=3){
				//       //3he add another proton
				//       un_ff(1,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);
				//       uquark = un_pdf(1,x,Q2);
				//       dquark = un_pdf(2,x,Q2);
				//       squark = un_pdf(3,x,Q2);
				//       ubarquark = un_pdf(-1,x,Q2);
				//       dbarquark = un_pdf(-2,x,Q2);
				//       sbarquark = un_pdf(-3,x,Q2);

				//       df_hp += qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
				//       df_hm += qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;


				//     }//

				// mass_hadron is pion, mom_hadr, 2.6 is lifetime and 3.0 is c


			}else{
				//kaon
				// if (target_flag>=1){
				// hydrogen
				un_ff(2,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);
				uquark = un_pdf(1,x,Q2);
				dquark = un_pdf(2,x,Q2);
				squark = un_pdf(3,x,Q2);
				ubarquark = un_pdf(-1,x,Q2);
				dbarquark = un_pdf(-2,x,Q2);
				sbarquark = un_pdf(-3,x,Q2);

				df_p_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
				df_p_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav;


				//  }
				// if (target_flag>=2){
				//deuterium add an neutron
				un_ff(2,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);
				uquark = un_pdf(2,x,Q2);
				dquark = un_pdf(1,x,Q2);
				squark = un_pdf(3,x,Q2);
				ubarquark = un_pdf(-2,x,Q2);
				dbarquark = un_pdf(-1,x,Q2);
				sbarquark = un_pdf(-3,x,Q2);

				df_n_hp = qu*qu*uquark*D_fav   + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
				df_n_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav   + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav   + qs*qs*sbarquark*D_unfav;
				// }
				//     if (target_flag>=3){
				//3he add another proton
				//       un_ff(2,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);
				//       uquark = un_pdf(1,x,Q2);
				//       dquark = un_pdf(2,x,Q2);
				//       squark = un_pdf(3,x,Q2);
				//       ubarquark = un_pdf(-1,x,Q2);
				//       dbarquark = un_pdf(-2,x,Q2);
				//       sbarquark = un_pdf(-3,x,Q2);

				//       df_hp += qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
				//       df_hm += qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav;
				//     }

			}

			//cout << df_hp/x << jacoF << "\t" << "haha" << endl;

			//xfx is in fact x*pdf

			*(dxs_all+0) = *dxs_hp * df_p_hp/x;
			*(dxs_all+1) = *dxs_hp * df_n_hp/x;
			*(dxs_all+2) = *dxs_hm * df_p_hm/x;
			*(dxs_all+3) = *dxs_hm * df_n_hm/x;


			if (target_flag==1){
				*dxs_hp *= df_p_hp/x;
				*dxs_hm *= df_p_hm/x;
			}else if (target_flag==2){
				*dxs_hp *= (df_p_hp + df_n_hp)/x;
				*dxs_hm *= (df_p_hm + df_n_hm)/x;
			}else if (target_flag==3){
				*dxs_hp *= (df_p_hp*2.+df_n_hp)/x;
				*dxs_hm *= (df_p_hm*2.+df_n_hm)/x;
			}

		}
		/*}}}*/

		/*double un_pdf(int iparton,double x,double Q2){{{*/          
		double un_pdf(int iparton,double x,double Q2)          
		{
			double Q = sqrt(Q2);
			double result = 0;

			if (iparton==1){
				// u quark
				result = xfx(x, Q, 2);
			}else if (iparton==2){
				// d quark
				result = xfx(x, Q, 1);
			}else if (iparton==-1){
				// \bar{u} quark
				result = xfx(x, Q, -2);
			}else if (iparton==-2){
				// \bar{d} quark
				result = xfx(x, Q, -1);
			}else if (iparton==3){
				// strange quark
				result = xfx(x, Q, 3);
			}else if (iparton==-3){
				// \bar{s} quark
				result = xfx(x, Q, -3);
			}else if (iparton==0){
				result = xfx(x,Q,0);
			}

			return(result);
		}
		/*}}}*/

		/*double un_ff(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s, double* D_g){{{*/
		double un_ff(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s, double* D_g){
			//flag ==1 for pion, flag==2 for kaon
			double Q2zero=2.0;
			double lambda=0.227;

			if (z>1) z=1;

			double sv = log(log(Q2/lambda/lambda)/log(Q2zero/lambda/lambda));

			if (flag == 1){
				double N = 1.150 - 1.522*sv + 1.378*pow(sv,2) - 0.527*pow(sv,3);
				double a1 = -0.740 - 1.680*sv + 1.546*pow(sv,2) - 0.596*pow(sv,3);
				double a2 = 1.430 + 0.543*sv - 0.023*pow(sv,2);
				double Ns = 4.250 - 3.147*sv + 0.755*pow(sv,2);
				double a1s = -0.770 -0.573*sv + 0.117*pow(sv,2);
				double a2s = 4.48 + 0.890*sv - 0.138*pow(sv,2);

				double Ng = 5.530 - 9.228*sv + 5.192*sv*sv - 0.966 * sv*sv*sv;
				double a1g = -0.320 + 0.318*sv - 0.561*sv*sv;
				double a2g = 2.7+2.553*sv-0.907*sv*sv;
				double a3g = 0.751*sv+0.496*sv*sv;

				double D_sum = N*pow(z,a1)*pow((1.0-z),a2);
				double D_sum_s = Ns*pow(z,a1s)*pow((1.0-z),a2s);
				double R_D = pow((1.0-z),0.083583)/pow((1.0+z),1.9838);
				*D_fav = D_sum/(1.0+R_D);
				*D_unfav = D_sum/(1.0+1.0/R_D);
				*D_s = D_sum_s/2.0;
				*D_g = Ng * pow(z,a1g) *pow(1-z,a2g)*(1+a3g/z)/2.;

			}else{
				double N = 0.310 - 0.038*sv - 0.042*pow(sv,2);
				double a1 = -0.980 - 0.260*sv + 0.008*pow(sv,2);
				double a2 = 0.970 + 0.978*sv - 0.229*pow(sv,2);
				double Ns =  1.080 - 0.469*sv + 0.003*pow(sv,2);
				double a1s = -0.820 -0.240*sv - 0.035*pow(sv,2);
				double a2s = 2.550 + 1.026*sv - 0.246*pow(sv,2);

				double Ng = 0.310 - 0.325*sv -0.092*sv*sv;
				double a1g = -0.17 - 0.214*sv-0.184*sv*sv;
				double a2g = 0.89 + 2.185*sv-0.471*sv*sv;
				double a3g = 1.154*sv-0.026*sv*sv;

				double D_sum = N*pow(z,a1)*pow((1.0-z),a2);
				double D_sum_s = Ns*pow(z,a1s)*pow((1.0-z),a2s);
				double R_D = pow((1.0-z),0.083583)/pow((1.0+z),1.9838);
				*D_fav = D_sum/(1.0+R_D);
				*D_unfav = D_sum/(1.0+1.0/R_D);
				*D_s = D_sum_s/2.0;
				*D_g = Ng * pow(z,a1g) *pow(1-z,a2g)*(1+a3g/z)/2.;
			}

			return 0;
		}
		/*}}}*/

		/*void dxs_hpt(double pt_tmp, double* dxs_hp,...){{{*/
		void dxs_hpt(double pt_tmp, double* dxs_hp, double* dxs_hm,double* dxs_all){

			double alpha_s;
			alpha_s = alphasPDF(sqrt(Q2));


			//cout << y << "\t" << Q2 << endl;
			//cout << alpha_s << endl;

			*dxs_hp = 1./137.035/137.035/16./PI/PI/Q2/Q2*y/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.;
			*dxs_hm = 1./137.035/137.035/16./PI/PI/Q2/Q2*y/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.;

			double qu=2./3.;
			double qd=-1./3.;
			double qs=-1./3.;

			double Pqq, Pqg, Pgq;
			double zp,xp;
			double df_p_hp = 0,df_p_hm = 0,df_n_hp = 0,df_n_hm = 0;

			double uquark,dquark,squark,ubarquark,dbarquark,sbarquark,gluon;
			double D_fav,D_unfav,D_s,D_g;

			// doing the integral
			for (Int_t i=0;i!=100;i++){
				xp = x + (1-x)/100.*(i+0.5);
				zp =  z / (z + xp*pt_tmp*pt_tmp/z/(1-xp)/Q2);


				if (z/zp>1) continue;


				Pqq = 64*PI/3.*Q2/y/y*((1+(1-y)*(1-y))*((1-xp)*(1-zp) + (1+xp*xp*zp*zp)/(1-xp)/(1-zp)) 
						+ 8*xp*zp*(1-y) //- 4 *sqrt(xp*zp*(1-y)/(1-xp)/(1-zp))*(2-y)*(xp*zp+(1-xp)*(1-zp))*cos(phi_h)
						//+ 4*xp*zp*(1-y)*cos(2. * phi_h)
						) *(1-x)/100. /(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2); //13

				Pqg = 64*PI/3.*Q2/y/y*( (1+(1-y)*(1-y))*((1-xp)*zp + (1+xp*xp*(1-zp)*(1-zp))/(1-xp)/zp) 
						+ 8*xp*(1-y)*(1-zp) 
						//+ 4.*sqrt(xp*(1-y)*(1-zp)/(1-xp)/zp)*(2-y)*(xp*(1-zp)+(1-xp)*zp)*cos(phi_h)
						//+ 4.*xp*(1-y)*(1-zp)*cos(2.*phi_h)
						)*(1-x)/100. /(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2);  //14

				Pgq = 8.*PI*Q2/y/y*((1+(1-y)*(1-y))*(xp*xp+(1-xp)*(1-xp))*(zp*zp+(1-zp)*(1-zp))/zp/(1-zp) 
						+ 16*xp*(1-xp)*(1-y)
						//-4.*sqrt(xp*(1-xp)*(1-y)/zp/(1-zp))*(2-y)*(1-2*xp)*(1-2.*zp)*cos(phi_h)
						//+ 8*xp*(1-xp)*(1-y)*cos(2.*phi_h)
						)*(1-x)/100. /(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2);  //15

				//cout << xp << "\t" << zp << "\t" << z/zp << "\t" << x/xp << "\t" << Pqq << "\t" << Pqg << "\t" << Pgq <<"\t" << (1-x)/100./(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2)<<  endl;
				// Pqq = 1.;
				//     Pqg = 0.;
				//     Pgq = 0.;

				if (fabs(particle_flag)==1){
					//pion
					// if (target_flag >=1){
					//hydrogen add one proton
					un_ff(1,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
					uquark = un_pdf(1,x/xp,Q2);
					dquark = un_pdf(2,x/xp,Q2);
					squark = un_pdf(3,x/xp,Q2);
					ubarquark = un_pdf(-1,x/xp,Q2);
					dbarquark = un_pdf(-2,x/xp,Q2);
					sbarquark = un_pdf(-3,x/xp,Q2);
					gluon = un_pdf(0,x/xp,Q2);

					df_p_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
					df_p_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s)*gluon;
					// }
					//       if (target_flag >=2){
					//deuterium add another neutron
					un_ff(1,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
					uquark = un_pdf(2,x/xp,Q2);
					dquark = un_pdf(1,x/xp,Q2);
					squark = un_pdf(3,x/xp,Q2);
					ubarquark = un_pdf(-2,x/xp,Q2);
					dbarquark = un_pdf(-1,x/xp,Q2);
					sbarquark = un_pdf(-3,x/xp,Q2);
					gluon = un_pdf(0,x/xp,Q2);

					df_n_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
					df_n_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s)*gluon;
					// }
					//       if (target_flag >=3){
					// 	//he3 add another proton
					// 	un_ff(1,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
					// 	uquark = un_pdf(1,x/xp,Q2);
					// 	dquark = un_pdf(2,x/xp,Q2);
					// 	squark = un_pdf(3,x/xp,Q2);
					// 	ubarquark = un_pdf(-1,x/xp,Q2);
					// 	dbarquark = un_pdf(-2,x/xp,Q2);
					// 	sbarquark = un_pdf(-3,x/xp,Q2);
					// 	gluon = un_pdf(0,x/xp,Q2);

					// 	df_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
					// 	df_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s)*gluon;
					//       }c
				}else{
					//kaon
					// if (target_flag >=1){
					//hydrogen add one proton
					un_ff(2,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
					uquark = un_pdf(1,x/xp,Q2);
					dquark = un_pdf(2,x/xp,Q2);
					squark = un_pdf(3,x/xp,Q2);
					ubarquark = un_pdf(-1,x/xp,Q2);
					dbarquark = un_pdf(-2,x/xp,Q2);
					sbarquark = un_pdf(-3,x/xp,Q2);
					gluon = un_pdf(0,x/xp,Q2);
					df_p_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
					df_p_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);

					// }
					// if (target_flag >=2){
					//deuterium add another neutron
					un_ff(2,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
					uquark = un_pdf(2,x/xp,Q2);
					dquark = un_pdf(1,x/xp,Q2);
					squark = un_pdf(3,x/xp,Q2);
					ubarquark = un_pdf(-2,x/xp,Q2);
					dbarquark = un_pdf(-1,x/xp,Q2);
					sbarquark = un_pdf(-3,x/xp,Q2);
					gluon = un_pdf(0,x/xp,Q2);
					df_n_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
					df_n_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
					// }
					// if (target_flag >=3){
					// 	//he3 add another proton
					// 	un_ff(2,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
					// 	uquark = un_pdf(1,x/xp,Q2);
					// 	dquark = un_pdf(2,x/xp,Q2);
					// 	squark = un_pdf(3,x/xp,Q2);
					// 	ubarquark = un_pdf(-1,x/xp,Q2);
					// 	dbarquark = un_pdf(-2,x/xp,Q2);
					// 	sbarquark = un_pdf(-3,x/xp,Q2);
					// 	gluon = un_pdf(0,x/xp,Q2);
					// 	df_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
					// 	df_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
					//       }
				}
			}
			//cout << *dxs_hp << " \t" << df_hp/x << "\t" << jacoF << endl;

			*(dxs_all+0) = *dxs_hp * df_p_hp/x*jacoF;
			*(dxs_all+1) = *dxs_hp * df_n_hp/x*jacoF;
			*(dxs_all+2) = *dxs_hm * df_p_hm/x*jacoF;
			*(dxs_all+3) = *dxs_hm * df_n_hm/x*jacoF;

			if (target_flag==1){
				*dxs_hp = *dxs_hp * df_p_hp*jacoF/x;
				*dxs_hm = *dxs_hm * df_p_hm*jacoF/x;
			}else if (target_flag==2){
				*dxs_hp = *dxs_hp * (df_p_hp+df_n_hp)*jacoF/x;
				*dxs_hm = *dxs_hm * (df_p_hm+df_n_hm)*jacoF/x;
			}else if (target_flag==3){
				*dxs_hp = *dxs_hp * (df_p_hp*2.+df_n_hp)*jacoF/x;
				*dxs_hm = *dxs_hm * (df_p_hm*2.+df_n_hp)*jacoF/x;
			}


		}
		/*}}}*/

		/*Kinematic Quantties{{{*/
	private:
		double mass_hadron;
		int	target_flag;
		int	particle_flag;
			
	public:
		double mom_ele;
		double theta_ele; 
		double phi_ele; 
		double mom_had;
		double theta_had;
		double phi_had;
        
		double theta_q;
		double theta_s;
		double phi_q;
		double phi_s;
		double phi_h;

		double px_ele;
		double py_ele;
		double pz_ele;
		double E0_ele;
		double px_had;
		double py_had;
		double pz_had;
		double E0_had;

		double Q2;
		double W;
		double Wp;
		double x;
		double y;
		double z;
		double pt;
		double nu;
		double s;
		double gamma;
		double epsilon;
		double jacoF;
		double XS_HP;
		double XS_HM;
		double dilute_hp;
		double dilute_hm;
        /*}}}*/
};
