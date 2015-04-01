#include <TROOT.h>
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"
#include "TLorentzVector.h"

using namespace std;
#define PI 3.14159265
#define DEG 180./3.14159265

extern "C"
{
	void dxs_dis__(const int*,const int*,const float* ,      //target,particle,beam_energy
			const float*, const float*,const float*,//p_e,theta_e,phi_e 
			const float*, const float*,const float*,//p_h,theta_h,phi_h 
			float*, float*, float*, float*);//sdxs_e,sdxs_h,dxs_hp, dxs_hm
}

Double_t Azimuthalphi(Double_t vx, Double_t vy);

int main(Int_t argc, char *argv[]){

	gRandom->SetSeed(0);// uses time for seed

	Double_t mass_e = 0.511e-3; // electron mass in GeV
	Double_t mass_p = 0.93827;
	Double_t mass_pi = 0.13957;
	Double_t mass_pi0 = 0.1349766;
	Double_t mass_kaon = 0.493667;
	Double_t mass_n = 0.939566;

	if (argc != 7) {
		cout << "usage: ./GetSIDIS #ele_mom #ion_mom #target #particle #fileno #events" << endl;
		cout << "both mom are in GeV, set ion_mom=0 for fix target" << endl;
		cout << "#target_flag = 1 for proton and 2 for deuteron 3 for 3he" << endl;
		cout << "#particle_flag = 1 for pion+, -1 for pion- and 2 for kaon+, -2 for kaon-" << endl;
		cout << "#fileno is the file number of output, used for batch" << endl;
		cout << "#events is number of event in each file" << endl;
	}
	else{
		/*Inputs&Output{{{*/
		//initialize
		Double_t momentum_ele = atof(argv[1]);
		Double_t momentum_ion = atof(argv[2]);
		Int_t target_flag = atoi(argv[3]);
		Int_t particle_flag = atoi(argv[4]);
		Int_t fileno = atoi(argv[5]);
		Int_t number_of_events = atoi(argv[6]);

		//define output file
		TString filename,prefix;
		prefix="sidis";

		Double_t mass_target;
		if (target_flag ==1){
			prefix += "_p";
			mass_target = mass_p;
		}else if (target_flag == 2){
			prefix += "_d";
			mass_target = (mass_n+mass_p)/2.;
			momentum_ion = momentum_ion/2.;
		}else if (target_flag == 3){
			prefix += "_3he";
			mass_target = (mass_n+2.*mass_p)/3.;
			momentum_ion = momentum_ion*2./3.;
		}else{
			cout << "target_flag is wrong! 1,2,3 " << endl;
			return 0;
		}

		Double_t mass_hadron;

		if (particle_flag == 1){
			prefix += "_pip";
			mass_hadron = mass_pi;
		}else if (particle_flag == 2){
			prefix += "_pim";
			mass_hadron = mass_pi;
		}else if (particle_flag == 3){
			prefix += "_kp";
			mass_hadron = mass_kaon;
		}else if (particle_flag == 4){
			prefix += "_km";
			mass_hadron = mass_kaon;
		}else{
			cout << "particle_flag is wrong 1=pi+,2=pi-,3=k+ and 4=k-" << endl;
			return 0;
		}

		//create filename
		filename.Form("_%d_%d_1_%d.root",Int_t(momentum_ele),Int_t(momentum_ion),Int_t(fileno));
		filename = prefix + filename;
		TFile *file1 = new TFile(filename,"RECREATE");
		TTree *t1 = new TTree("T","T");
		t1->SetDirectory(file1);

		filename.Form("_%d_%d_2_%d.root",Int_t(momentum_ele),Int_t(momentum_ion),Int_t(fileno));
		filename = prefix + filename;
		TFile *file2 = new TFile(filename,"RECREATE");
		TTree *t2 = new TTree("T","T");
		t2->SetDirectory(file2);
		/*}}}*/

		/*Define{{{*/
		// //define the 4-momentum
		//define electron direction as +z assuming a proton/neutron for the ion mass now 
		// approximation for SIDIS
		TLorentzVector *P4_ini_ele = new TLorentzVector(0.,0.,momentum_ele,sqrt(momentum_ele*momentum_ele + mass_e*mass_e));
		TLorentzVector *P4_ini_ion = new TLorentzVector(0.,0.,-momentum_ion,sqrt(momentum_ion*momentum_ion + mass_target*mass_target));
		TLorentzVector *P4_fin_had = new TLorentzVector(0.,0.,0.,1.);
		TLorentzVector *P4_fin_ele = new TLorentzVector(0.,0.,0.,1.);
		TLorentzVector *P4_q = new TLorentzVector(0.,0.,0.,1.);

		TLorentzVector *lrz_P4_q = new TLorentzVector(0.,0.,0.,1.);
		TLorentzVector *lrz_P4_h = new TLorentzVector(0.,0.,0.,1.);
		TLorentzVector *lrz_P4_ef = new TLorentzVector(0.,0.,0.,1.);
		TLorentzVector *lrz_P4_ion = new TLorentzVector(0.,0.,0.,1.);
		TLorentzVector *lrz_P4_ei = new TLorentzVector(0.,0.,0.,1.);

		Double_t vn=momentum_ion/sqrt(momentum_ion*momentum_ion + mass_target*mass_target);
		TVector3 vnboost(0.,0.,vn);

		Double_t theta_gen, phi_gen,mom_gen;
		Double_t theta_gen_ele, phi_gen_ele,mom_gen_ele;
		Double_t theta_gen_had, phi_gen_had,mom_gen_had;
		float beam_energy=0.;
		float theta_ele_temp=0., phi_ele_temp=0.,mom_ele_temp=0.;
		float theta_had_temp=0., phi_had_temp=0.,mom_had_temp=0.;
	    Double_t Q2,W,Wp,x,y,z,pt,nu,s;

		Int_t count[2] = {0,0};

		Double_t gamma,epsilon;

		TVector3 p3_q,p3_fin_ele,p3_fin_had;
		TVector3 p3_target_spin(0.,1.,0);
		Double_t theta_q,phi_q;
		Double_t theta_s,phi_s,phi_h;

		Double_t mom_ele,mom_had;
		Double_t theta_ele,theta_had;
		Double_t phi_ele,phi_had;
		Double_t mom_pro,energy_pro;
		Double_t mom_ini_ele,energy_ini_ele;
		Double_t sdxs_e,sdxs_h,dxs_hp,dxs_hm;
		float sdxs_e_temp,sdxs_h_temp,dxs_hp_temp,dxs_hm_temp;
		Int_t nsim = 0;
		mom_ini_ele = fabs(momentum_ele);
		energy_ini_ele = sqrt(momentum_ele*momentum_ele + mass_e*mass_e);

		mom_pro = fabs(momentum_ion);
		energy_pro = sqrt(momentum_ion*momentum_ion + mass_target*mass_target);
        /*}}}*/

		/*New Branches{{{*/
		t1->Branch("Q2",&Q2,"data/D");
		t1->Branch("W",&W,"data/D");
		t1->Branch("Wp",&Wp,"data/D");
		t1->Branch("x",&x,"data/D");
		t1->Branch("y",&y,"data/D");
		t1->Branch("z",&z,"data/D");
		t1->Branch("nu",&nu,"data/D");
		t1->Branch("s",&s,"data/D");
		t1->Branch("pt",&pt,"data/D");
		t1->Branch("theta_q",&theta_q,"data/D");
		t1->Branch("theta_s",&theta_s,"data/D");
		t1->Branch("phi_h",&phi_h,"data/D");
		t1->Branch("phi_s",&phi_s,"data/D");
		t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t1->Branch("sdxs_e",&sdxs_e,"sdxs_e/D");
		t1->Branch("sdxs_h",&sdxs_h,"sdxs_h/D");
		t1->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t1->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
		t1->Branch("mom_had",&mom_had,"mom_had/D");
		t1->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
		t1->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t1->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
		t1->Branch("theta_had",&theta_had,"theta_had/D");
		t1->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
		t1->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t1->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
		t1->Branch("phi_had",&phi_had,"phi_had/D");
		t1->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
		t1->Branch("nsim",&nsim,"nsim/I");

		t2->Branch("Q2",&Q2,"data/D");
		t2->Branch("W",&W,"data/D");
		t2->Branch("Wp",&Wp,"data/D");
		t2->Branch("x",&x,"data/D");
		t2->Branch("y",&y,"data/D");
		t2->Branch("z",&z,"data/D");
		t2->Branch("nu",&nu,"data/D");
		t2->Branch("s",&s,"data/D");
		t2->Branch("pt",&pt,"data/D");
		t2->Branch("theta_q",&theta_q,"data/D");
		t2->Branch("theta_s",&theta_s,"data/D");
		t2->Branch("phi_h",&phi_h,"data/D");
		t2->Branch("phi_s",&phi_s,"data/D");
		t2->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t2->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t2->Branch("sdxs_e",&sdxs_e,"sdxs_e/D");
		t2->Branch("sdxs_h",&sdxs_h,"sdxs_h/D");
		t2->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t2->Branch("mom_had",&mom_had,"mom_had/D");
		t2->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t2->Branch("theta_had",&theta_had,"theta_had/D");
		t2->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t2->Branch("phi_had",&phi_had,"phi_had/D");
		t2->Branch("nsim",&nsim,"nsim/I");
        /*}}}*/

		bool exitcondition=true;
		while(exitcondition){
			nsim ++;

			/*Generator{{{*/
			//For electron
			phi_gen = gRandom->Uniform(0.,2.*PI);
			//       mom_gen = gRandom->Uniform(0.7,momentum_ele + momentum_ion);
			theta_gen = acos(gRandom->Uniform(cos(30/180.*3.1415926),cos(7/180.*3.1415926)));
			mom_gen = gRandom->Uniform(0.5,momentum_ele);
            
			mom_gen_ele = mom_gen; theta_gen_ele = theta_gen*DEG;  phi_gen_ele = phi_gen*DEG;

			P4_fin_ele->SetPxPyPzE(mom_gen*sin(theta_gen)*cos(phi_gen),
					mom_gen*sin(theta_gen)*sin(phi_gen),
					mom_gen*cos(theta_gen)
					,sqrt(mom_gen*mom_gen+mass_e*mass_e));
			
			//For hadron
			theta_gen = acos(gRandom->Uniform(cos(30/180.*3.1415926),cos(7/180.*3.1415926)));
			phi_gen = gRandom->Uniform(0.,2.*PI);
			//      mom_gen = gRandom->Uniform(0.7,momentum_ele + momentum_ion);
			mom_gen = gRandom->Uniform(0.5,6);
			
			mom_gen_had = mom_gen; theta_gen_had = theta_gen*DEG;  phi_gen_had = phi_gen*DEG;

			P4_fin_had->SetPxPyPzE(mom_gen*sin(theta_gen)*cos(phi_gen),
					mom_gen*sin(theta_gen)*sin(phi_gen),
					mom_gen*cos(theta_gen)
					,sqrt(mom_gen*mom_gen+mass_hadron*mass_hadron));

			*P4_q = *P4_ini_ele - *P4_fin_ele;
			Q2 = - (*P4_q)*(*P4_q);
			W = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele);
			Wp = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had);
            /*}}}*/

			if (Q2 >=1.0 && W>= 2.3*2.3 &&Wp>= 1.6*1.6){
				/*Boost{{{*/
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
                /*}}}*/

				if ( z>0.3&&z<0.7 && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) || (count[1]<number_of_events&&pt>1.0&&Q2<=10.)) ){

					/*Lab Frame{{{*/
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
                    /*}}}*/

					/*Cross Section{{{*/
					sdxs_e = 0.0; sdxs_h = 0.0; dxs_hp = 0.0;dxs_hm = 0.0;

					beam_energy = (float) (momentum_ele);
					mom_ele_temp = (float)(mom_gen_ele);
					mom_had_temp = (float)(mom_gen_had);
					theta_ele_temp = (float)(theta_gen_ele);
					theta_had_temp = (float)(theta_gen_had);
					phi_ele_temp = (float)(phi_gen_ele);
					phi_had_temp = (float)(phi_gen_had);

					//cerr<<Form("--- Ee=%e, The=%e, Phe=%e",mom_ele_temp,theta_ele_temp,phi_ele_temp)<<endl;
					//cerr<<Form("--- Eh=%e, Thh=%e, Phh=%e",mom_had_temp,theta_had_temp,phi_had_temp)<<endl;

					dxs_dis__(&target_flag,&particle_flag,&beam_energy,
							&mom_ele_temp,&theta_ele_temp,&phi_ele_temp,
							&mom_had_temp,&theta_had_temp,&phi_had_temp,
							&sdxs_e_temp,&sdxs_h_temp,&dxs_hp_temp,&dxs_hm_temp);

					sdxs_e = (double)(sdxs_e_temp);
					sdxs_h = (double)(sdxs_h_temp);
					dxs_hp = (double)(dxs_hp_temp);
					dxs_hm = (double)(dxs_hm_temp);
					
					//cerr<<Form("--- Xe=%e, Xh=%e, Xhp=%e, Xhm=%e",sdxs_e,sdxs_h,dxs_hp,dxs_hm)<<endl;

					// try to take care of the decay
/*					
					Double_t decay_part;
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

					dxs_hp = decay_part * dxs_hp; 
					dxs_hm = decay_part * dxs_hm;
					*/
		
					if (dxs_hp>0.&& dxs_hp<10000.){
					}else{
						dxs_hp = 0.;
					}
					if (dxs_hm>0.&& dxs_hm<10000.){
					}else{
						dxs_hm = 0.;
					}
					/*}}}*/

					if ((dxs_hp+dxs_hm)!=0){
						if (Q2<=10.&&pt<=1.0){
							t1->Fill();
							count[0] ++;//cout << 0 << " " << count[0] << endl;
						}
						if (Q2<=10.&&pt>1.0){
							t2->Fill();
							count[1] ++;//cout << 1 << " " << count[1] << endl;
						}
					}
					//cout << nsim << endl;

					cout << count[0] << "\t" << count[1] << "\r";

				}
			}
			///judging exitcondition
			if (count[0] < number_of_events || count[1] < number_of_events) exitcondition=true;
			else exitcondition=false;

		}
		cout << count[0] << "\t" << count[1]  << endl;

		file1->Write();
		file1->Close();
		file2->Write();
		file2->Close();
	}
	return 0;
} 



Double_t Azimuthalphi(Double_t vx, Double_t vy){
	Double_t pmod, phi, cosf;
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
