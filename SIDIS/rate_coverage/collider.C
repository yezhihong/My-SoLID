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

#include "LHAPDF/LHAPDF.h"
using namespace std;
using namespace LHAPDF;

#define PI 3.14159265
#define DEG 180./3.1415926

Double_t Azimuthalphi(Double_t vx, Double_t vy);
Double_t Jacobian(Double_t mom_ele,Double_t theta_ele,Double_t phi_ele,Double_t mom_had,Double_t theta_had,Double_t phi_had,Double_t mom_pro,Double_t energy_pro,Double_t mom_ini_ele,Double_t mass_hadron);
double un_ff(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s,double* D_g);
double un_pdf(int iparton,double x,double Q2);
void dxs(Double_t x, Double_t y, Double_t Q2, Double_t epsilon, Double_t gamma, Double_t z, Double_t pt, Double_t jacoF, Int_t target_flag, Int_t particle_flag,Double_t bpt_m, Double_t bpt_p, Double_t* dxs_hp, Double_t* dxs_hm,Double_t* dxs_all);
void dxs_hpt(Double_t x, Double_t y, Double_t Q2, Double_t z, Double_t pt, Double_t phi_h, Double_t jacoF, Int_t target_flag, Int_t particle_flag, Double_t* dxs_hp, Double_t* dxs_hm,Double_t* dxs_all);

int main(Int_t argc, char *argv[]){

	const int SUBSET = 1;
	const string NAME = "cteq6m";
	gRandom->SetSeed(0);// uses time for seed

	char *LHAPDF_path = getenv("LHAPDF");
	sprintf(LHAPDF_path,"%s/share/lhapdf",LHAPDF_path);
	setPDFPath(LHAPDF_path);
	LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);

	Double_t mass_e = 0.511e-3; // electron mass in GeV
	Double_t mass_p = 0.93827;
	Double_t mass_pi = 0.13957;
	Double_t mass_pi0 = 0.1349766;
	Double_t mass_kaon = 0.493667;
	Double_t mass_n = 0.939566;

	Double_t bpt_m=4.694;
	Double_t bpt_p=4.661;

	if (argc != 8) {
		cout << "usage: ./collider #ele_mom #ion_mom #target #particle #fileno #events config" << endl;
		cout << "both mom are in GeV, set ion_mom=0 for fix target" << endl;
		cout << "#target_flag = 1 for proton and 2 for deuteron 3 for 3he" << endl;
		cout << "#particle_flag = 1 for pion+, -1 for pion- and 2 for kaon+, -2 for kaon-" << endl;
		cout << "#fileno is the file number of output, used for batch" << endl;
		cout << "#events is number of event in each file" << endl;
		cout << "config is 'EIC' or 'SoLID' which needs ion_mom=0" << endl;
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
		string config = argv[7];

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
		}else if (particle_flag == -1){
			prefix += "_pim";
			mass_hadron = mass_pi;
		}else if (particle_flag == 2){
			prefix += "_kp";
			mass_hadron = mass_kaon;
		}else if (particle_flag == -2){
			prefix += "_km";
			mass_hadron = mass_kaon;
		}else{
			cout << "particle_flag is wrong +-1 and +-2" << endl;
			return 0;
		}

		if (config != "EIC" && config != "SoLID"){
			cout << "not supported config" << endl;
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

		filename.Form("_%d_%d_3_%d.root",Int_t(momentum_ele),Int_t(momentum_ion),Int_t(fileno));
		filename = prefix + filename;
		TFile *file3 = new TFile(filename,"RECREATE");
		TTree *t3 = new TTree("T","T");
		t3->SetDirectory(file3);

		filename.Form("_%d_%d_4_%d.root",Int_t(momentum_ele),Int_t(momentum_ion),Int_t(fileno));
		filename = prefix + filename;
		TFile *file4 = new TFile(filename,"RECREATE");
		TTree *t4 = new TTree("T","T");
		t4->SetDirectory(file4);
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

		//  t1->Branch("P4_ini_ele","TLorentzVector",&P4_ini_ele,16000,0);
		//     t1->Branch("P4_ini_ion","TLorentzVector",&P4_ini_ion,16000,0);
		//     t1->Branch("P4_fin_ele","TLorentzVector",&P4_fin_ele,16000,0);
		//     t1->Branch("P4_fin_ion","TLorentzVector",&P4_fin_had,16000,0);
		//     t1->Branch("P4_q","TLorentzVector",&P4_q,16000,0);

		Double_t theta_gen, phi_gen,mom_gen;
		Double_t Q2,W,Wp,x,y,z,pt,nu,s;

		Int_t count[4] = {0,0,0,0};

		Double_t gamma,epsilon;

		TVector3 p3_q,p3_fin_ele,p3_fin_had;
		TVector3 p3_target_spin(0.,1.,0);
		Double_t theta_q,phi_q;
		Double_t theta_s,phi_s,phi_h;


		Double_t mom_ele,mom_had;
		Double_t mom_gen_ele,mom_gen_had;
		Double_t theta_ele,theta_had;
		Double_t theta_gen_ele,theta_gen_had;
		Double_t phi_ele,phi_had;
		Double_t phi_gen_ele,phi_gen_had;
		Double_t mom_pro,energy_pro;
		Double_t mom_ini_ele,energy_ini_ele;
		Double_t jacoF,dxs_hp,dxs_hm,phase_hp,phase_hm,dfff_hp,dfff_hm;
		Int_t nsim = 0;
		Double_t dilute[2];
		mom_ini_ele = fabs(momentum_ele);
		energy_ini_ele = sqrt(momentum_ele*momentum_ele + mass_e*mass_e);

		mom_pro = fabs(momentum_ion);
		energy_pro = sqrt(momentum_ion*momentum_ion + mass_target*mass_target);

		double electron_phase_space= 0.0, hadron_phase_space= 0.0, Phase_space= 0.0;
		if(config=="SoLID" ){
			electron_phase_space=(cos(7/DEG) - cos(30/DEG))*2*PI*(momentum_ele-0.5);   // theta: 7~30 degree,  2pi phi coverage, 0.5~11 GeV Momentum coverage 	
			hadron_phase_space=(cos(7/DEG) - cos(30/DEG))*2*PI*(6-0.5);  //theta, 7~30 degree,  2pi phi coverage, 0.5~6 GeV Momentum coverage
			Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
		}
		cout<<" -- For Config="<<config<<" Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;

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
		t1->Branch("jacoF",&jacoF,"jacoF/D");
		t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t1->Branch("phase_hp",&phase_hp,"phase_hp/D");
		t1->Branch("phase_hm",&phase_hm,"phase_hm/D");
		t1->Branch("dfff_hp",&dfff_hp,"dfff_hp/D");
		t1->Branch("dfff_hm",&dfff_hm,"dfff_hm/D");
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
		t1->Branch("dilute_p",&dilute[0],"data/D");
		t1->Branch("dilute_m",&dilute[1],"data/D");


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
		t2->Branch("jacoF",&jacoF,"jacoF/D");
		t2->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t2->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t2->Branch("phase_hp",&phase_hp,"phase_hp/D");
		t2->Branch("phase_hm",&phase_hm,"phase_hm/D");
		t2->Branch("dfff_hp",&dfff_hp,"dfff_hp/D");
		t2->Branch("dfff_hm",&dfff_hm,"dfff_hm/D");
		t2->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t2->Branch("mom_had",&mom_had,"mom_had/D");
		t2->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t2->Branch("theta_had",&theta_had,"theta_had/D");
		t2->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t2->Branch("phi_had",&phi_had,"phi_had/D");
		t2->Branch("nsim",&nsim,"nsim/I");
		t2->Branch("dilute_p",&dilute[0],"data/D");
		t2->Branch("dilute_m",&dilute[1],"data/D");

		t3->Branch("Q2",&Q2,"data/D");
		t3->Branch("W",&W,"data/D");
		t3->Branch("Wp",&Wp,"data/D");
		t3->Branch("x",&x,"data/D");
		t3->Branch("y",&y,"data/D");
		t3->Branch("z",&z,"data/D");
		t3->Branch("nu",&nu,"data/D");
		t3->Branch("s",&s,"data/D");
		t3->Branch("pt",&pt,"data/D");
		t3->Branch("theta_q",&theta_q,"data/D");
		t3->Branch("theta_s",&theta_s,"data/D");
		t3->Branch("phi_h",&phi_h,"data/D");
		t3->Branch("phi_s",&phi_s,"data/D");
		t3->Branch("jacoF",&jacoF,"jacoF/D");
		t3->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t3->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t3->Branch("phase_hp",&phase_hp,"phase_hp/D");
		t3->Branch("phase_hm",&phase_hm,"phase_hm/D");
		t3->Branch("dfff_hp",&dfff_hp,"dfff_hp/D");
		t3->Branch("dfff_hm",&dfff_hm,"dfff_hm/D");
		t3->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t3->Branch("mom_had",&mom_had,"mom_had/D");
		t3->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t3->Branch("theta_had",&theta_had,"theta_had/D");
		t3->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t3->Branch("phi_had",&phi_had,"phi_had/D");
		t3->Branch("nsim",&nsim,"nsim/I");
		t3->Branch("dilute_p",&dilute[0],"data/D");
		t3->Branch("dilute_m",&dilute[1],"data/D");

		t4->Branch("Q2",&Q2,"data/D");
		t4->Branch("W",&W,"data/D");
		t4->Branch("Wp",&Wp,"data/D");
		t4->Branch("x",&x,"data/D");
		t4->Branch("y",&y,"data/D");
		t4->Branch("z",&z,"data/D");
		t4->Branch("nu",&nu,"data/D");
		t4->Branch("s",&s,"data/D");
		t4->Branch("pt",&pt,"data/D");
		t4->Branch("theta_q",&theta_q,"data/D");
		t4->Branch("theta_s",&theta_s,"data/D");
		t4->Branch("phi_h",&phi_h,"data/D");
		t4->Branch("phi_s",&phi_s,"data/D");
		t4->Branch("jacoF",&jacoF,"jacoF/D");
		t4->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t4->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t4->Branch("phase_hp",&phase_hp,"phase_hp/D");
		t4->Branch("phase_hm",&phase_hm,"phase_hm/D");
		t4->Branch("dfff_hp",&dfff_hp,"dfff_hp/D");
		t4->Branch("dfff_hm",&dfff_hm,"dfff_hm/D");
		t4->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t4->Branch("mom_had",&mom_had,"mom_had/D");
		t4->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t4->Branch("theta_had",&theta_had,"theta_had/D");
		t4->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t4->Branch("phi_had",&phi_had,"phi_had/D");
		t4->Branch("nsim",&nsim,"nsim/I");
		t4->Branch("dilute_p",&dilute[0],"data/D");
		t4->Branch("dilute_m",&dilute[1],"data/D");
        /*}}}*/

		Double_t dxs_all[3][4];

		bool exitcondition=true;
		while(exitcondition){

			nsim ++;

			/*Generator{{{*/
			//For electron
			phi_gen = gRandom->Uniform(0.,2.*PI);
			//       mom_gen = gRandom->Uniform(0.7,momentum_ele + momentum_ion);
			if (config=="EIC") {
				theta_gen = acos(gRandom->Uniform(cos(150/180.*3.1415926),1.));
				mom_gen = gRandom->Uniform(0.7,momentum_ele * 3.);
			}
			else if (config=="SoLID") {
				theta_gen = acos(gRandom->Uniform(cos(30/180.*3.1415926),cos(7/180.*3.1415926)));
				mom_gen = gRandom->Uniform(0.5,momentum_ele);
			}

			mom_gen_ele = mom_gen; theta_gen_ele = theta_gen*DEG; phi_gen_ele = phi_gen*DEG;

			// mom_gen = 10.196474;
			//       theta_gen = 0.7361853;
			//       phi_gen = 0.3896692;

			P4_fin_ele->SetPxPyPzE(mom_gen*sin(theta_gen)*cos(phi_gen),
					mom_gen*sin(theta_gen)*sin(phi_gen),
					mom_gen*cos(theta_gen)
					,sqrt(mom_gen*mom_gen+mass_e*mass_e));
			
			//For hadron
			theta_gen = acos(gRandom->Uniform(cos(30/180.*3.1415926),cos(7/180.*3.1415926)));
			phi_gen = gRandom->Uniform(0.,2.*PI);
			//      mom_gen = gRandom->Uniform(0.7,momentum_ele + momentum_ion);
			mom_gen = gRandom->Uniform(0.5,6);

			mom_gen_had = mom_gen; theta_gen_had = theta_gen*DEG; phi_gen_had = phi_gen*DEG;
			
			// theta_gen = 1.8588846;
			//        phi_gen = 0.8115906;
			//        mom_gen = 2.3559179; 

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
				s = (*P4_ini_ele + *P4_ini_ion )*(*P4_ini_ele + *P4_ini_ion);
				nu = (*P4_ini_ion) * (*P4_q)/mass_target;
				x = Q2/ 2. / mass_target / nu;
				z = ((*P4_ini_ion)*(*P4_fin_had))/((*P4_ini_ion)*(*P4_q));
				y = ((*P4_ini_ion)*(*P4_q))/((*P4_ini_ion)*(*P4_ini_ele));
				W = sqrt(W);
				Wp = sqrt(Wp);

			//lrz_P4_q->SetPxPyPzE(P4_q->Px(),P4_q->Py(),P4_q->Pz(),P4_q->E());
			//lrz_P4_h->SetPxPyPzE(P4_fin_had->Px(),P4_fin_had->Py(),P4_fin_had->Pz(),P4_fin_had->E());
			//lrz_P4_ef->SetPxPyPzE(P4_fin_ele->Px(),P4_fin_ele->Py(),P4_fin_ele->Pz(),P4_fin_ele->E());

				*lrz_P4_q = *P4_q;
				*lrz_P4_h = *P4_fin_had;
				*lrz_P4_ef = *P4_fin_ele;

				lrz_P4_ef->Boost(vnboost);
				lrz_P4_h->Boost(vnboost); 
				lrz_P4_q->Boost(vnboost);

				// 	 *lrz_P4_ei = *P4_ini_ele;
				// *lrz_P4_ion = *P4_ini_ion;
				// 	 lrz_P4_ei->Boost(vnboost); 
				//lrz_P4_ion->Boost(vnboost);

				//cout << lrz_P4_ion->Pz() << endl;
				//cout << x << "\t" << (*lrz_P4_q)*(*lrz_P4_q)/(2.*(*lrz_P4_q)*(*lrz_P4_ion)) << endl;

				gamma = 2*mass_target *x/sqrt(Q2);
				epsilon = (1-y-0.25*gamma*gamma*y*y)/(1-y+0.5*y*y+0.25*gamma*gamma*y*y);

				for(Int_t j=0;j<3;j++){
					p3_q(j)=(*lrz_P4_q)(j);
					p3_fin_ele(j)=(*lrz_P4_ef)(j);
					p3_fin_had(j)=(*lrz_P4_h)(j);
				}

				pt = p3_fin_had.Perp(p3_q);

				if ( (config=="EIC" && z>0.2&&z<0.8&&y>0.05&&y<0.8 && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) 
								|| (count[1]<number_of_events&&pt>1.0&&Q2<=10.)
								|| (count[2]<number_of_events&&pt<=1.0&&Q2>10.)
								|| (count[3]<number_of_events&&pt>1.0&&Q2>10.)) ) || 
						(config=="SoLID" && z>0.3&&z<0.7 && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) || (count[1]<number_of_events&&pt>1.0&&Q2<=10.))) ){

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

					//cout << P4_fin_had->Pz() << "\t" << lrz_P4_h->Pz() << "\t" << lrz_P4_q->Pz() << "\t";
					//cout << pt*pt << "\t" << phi_h << "\t";

					mom_ele = P4_fin_ele->P();
					theta_ele = P4_fin_ele->Theta();
					mom_had = P4_fin_had->P();
					theta_had = P4_fin_had->Theta();
					phi_ele = P4_fin_ele->Phi();
					phi_had = P4_fin_had->Phi();

					jacoF = Jacobian(mom_ele,theta_ele,phi_ele,mom_had,theta_had,phi_had,mom_pro,energy_pro,momentum_ele,mass_hadron);
					//jacoF = 1.0;

					if (pt<0.8){
						/*pt<0.8{{{*/
						//first method 	 
						bpt_p = 1./(0.2+z*z*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
						bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)

						//second method  use the original value
						// bpt_m=4.694; // HERMES parameterization
						// bpt_p=4.661;
						dxs(x,y,Q2,epsilon,gamma,z,pt,jacoF,target_flag,particle_flag, bpt_m, bpt_p, &dxs_hp, &dxs_hm,dxs_all[0]);

						dilute[0] = dxs_all[0][0]/dxs_all[0][1];
						dilute[1] = dxs_all[0][2]/dxs_all[0][3];
                        /*}}}*/
					}
					else{
						/*pt>0.8{{{*/
						// this part is to generate factor K
						//make sure the DXS is the same at PT= 0.8 GeV
						//calculating the TMD part
						Double_t K[2],dxs_temp[10];
						bpt_p = 1./(0.2+z*z*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
						bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)
						//second method  use the original value
						// bpt_m=4.694; // HERMES parameterization
						// bpt_p=4.661;
						dxs(x,y,Q2,epsilon,gamma,z,0.8,jacoF,target_flag,particle_flag, bpt_m, bpt_p, &dxs_temp[0], &dxs_temp[1],dxs_all[0]);	   

						//calculating the TMD parts using the new vpt_p values
						bpt_p = 1./(0.25+z*z*0.28);// <pt^2> = 0.25 GeV^2 (quark internal momentum)
						bpt_m = bpt_p; // <kt^2> = 0.28 GeV^2 (struck quark internal momentum)
						//taking into account the NLO etc
						dxs(x,y,Q2,epsilon,gamma,z,1.0,jacoF,target_flag,particle_flag, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
						dxs_hpt(x,y,Q2,z,1.0,phi_h,jacoF,target_flag,particle_flag,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);



						K[0] = (dxs_temp[0]-dxs_temp[2])/dxs_temp[4];
						K[1] = (dxs_temp[1]-dxs_temp[3])/dxs_temp[5];

						//cout << pt <<"\t"<< Q2 << "\t" << jacoF << "\t" << dxs_temp[0] << "\t" << dxs_temp[2] << "\t" << dxs_temp[4] << "\t" << K[0] << "\t" << K[1] << endl;

						if (K[0]<0.) K[0] = 0.;
						if (K[1]<0.) K[1] = 0.;

						if (pt>1.2){
							dxs(x,y,Q2,epsilon,gamma,z,pt,jacoF,target_flag,particle_flag, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
							dxs_hpt(x,y,Q2,z,pt,phi_h,jacoF,target_flag,particle_flag,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

							dxs_hp = dxs_temp[2] + K[0]*dxs_temp[4];
							dxs_hm = dxs_temp[3] + K[1]*dxs_temp[5];

							dilute[0] = (dxs_all[1][0] + K[0]*dxs_all[2][0])/(dxs_all[1][1] + K[0]*dxs_all[2][1]);
							dilute[1] = (dxs_all[1][2] + K[1]*dxs_all[2][2])/(dxs_all[1][3] + K[1]*dxs_all[2][3]);

						}else{
							dxs(x,y,Q2,epsilon,gamma,z,1.2,jacoF,target_flag,particle_flag, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
							dxs_hpt(x,y,Q2,z,1.2,phi_h,jacoF,target_flag,particle_flag,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

							dxs_temp[2] = dxs_temp[2] + K[0]*dxs_temp[4];
							dxs_temp[3] = dxs_temp[3] + K[1]*dxs_temp[5];

							dxs_temp[4] = (dxs_temp[0]-dxs_temp[2])/(1./0.8/0.8-1./1.2/1.2);
							dxs_temp[5] = dxs_temp[0]-dxs_temp[4]/0.8/0.8;
							dxs_temp[6] = (dxs_temp[1]-dxs_temp[3])/(1./0.8/0.8-1./1.2/1.2);
							dxs_temp[7] = dxs_temp[1]-dxs_temp[6]/0.8/0.8;

							dxs_hp = dxs_temp[4]/pt/pt + dxs_temp[5];
							//dxs_hm = dxs_temp[6]/pt/pt + dxs_temp[6]; //Xin's original
							dxs_hm = dxs_temp[6]/pt/pt + dxs_temp[7]; //ZYE, second term [6]-->[7]

							dilute[0] = ((dxs_all[1][0] + K[0]*dxs_all[2][0])/(dxs_all[1][1] + K[0]*dxs_all[2][1])+dxs_all[0][0]/dxs_all[0][1])/2.;
							dilute[1] = ((dxs_all[1][2] + K[1]*dxs_all[2][2])/(dxs_all[1][3] + K[1]*dxs_all[2][3])+dxs_all[0][2]/dxs_all[0][3])/2.;
						}
						/*}}}*/
					}

					// try to take care of the decay
					/*Decay{{{*/
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

					if (dxs_hp>0.&& dxs_hp<10000.){
					}else{
						dxs_hp = 0.;
					}
					if (dxs_hm>0.&& dxs_hm<10000.){
					}else{
						dxs_hm = 0.;
					}
					//warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2


					if ((dxs_hp+dxs_hm)!=0){
						if (Q2<=10.&&pt<=1.0){
							t1->Fill();
							count[0] ++;//cout << 0 << " " << count[0] << endl;
						}
						if (Q2<=10.&&pt>1.0){
							t2->Fill();
							count[1] ++;//cout << 1 << " " << count[1] << endl;
						}
						if (Q2>10.&&pt<=1.0){
							t3->Fill();
							count[2] ++;//cout << 2 << " " << count[2] << endl;
						}
						if (Q2>10.&&pt>1.0){
							t4->Fill();
							count[3] ++;//cout << 3 << " " << count[3] <<  endl;
						}
					}
					//cout << nsim << endl;
                    /*}}}*/  
					cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << "\r";

				}
			}
			///judging exitcondition
			if (config=="EIC") {
				if (count[0] < number_of_events || count[1] < number_of_events || count[2] < number_of_events || count[3] < number_of_events) exitcondition=true;
				else exitcondition=false;
			} else if (config=="SoLID") {
				if (count[0] < number_of_events || count[1] < number_of_events) exitcondition=true;
				else exitcondition=false;
			} 

		}
		cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << endl;

		file1->Write();
		file1->Close();
		file2->Write();
		file2->Close();
		file3->Write();
		file3->Close();
		file4->Write();
		file4->Close();
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

Double_t Jacobian(Double_t P_ef,Double_t theta_ef,Double_t phi_ef,Double_t P_h,Double_t theta_h,Double_t phi_h,Double_t mom_pro,Double_t energy_pro,Double_t P_ei,Double_t mass_hadron){


	// need to shift the electron phi angle to 0 in order to do calculation for jacobian

	phi_h -= phi_ef;
	phi_ef = 0.;

	Double_t qbeta = -1.*mom_pro/energy_pro;
	//Double_t qbeta = 0.;
	Double_t qgamma = 1./sqrt(1-qbeta*qbeta);


	Double_t mass_pi = mass_hadron;
	Double_t P_hx,P_hy,P_hz;  
	P_hx = P_h*sin(theta_h)*cos(phi_h);
	P_hy = P_h*sin(theta_h)*sin(phi_h);
	P_hz = qgamma * (P_h*cos(theta_h) - qbeta * sqrt(P_h*P_h+mass_pi*mass_pi));



	//P_hz = P_h*cos(theta_h);

	Double_t qz = qgamma*(P_ei-P_ef*cos(theta_ef)-qbeta*(P_ei-P_ef));

	Double_t alpha,beta,gamma1;
	alpha = atan(P_ef*sin(theta_ef)*cos(phi_ef)/qz);
	beta = atan(-P_ef*sin(theta_ef)*sin(phi_ef)*
			qz/cos(alpha)/
			(pow(P_ef*sin(theta_ef)*sin(phi_ef),2) + 
			 pow(qz,2)));
	gamma1 = atan(sin(beta)/tan(alpha));

	Double_t a2,b2,c2,d2,e2,f2;
	a2 = (-sin(gamma1)*cos(alpha)+cos(gamma1)*sin(beta)*sin(alpha));
	b2 = cos(gamma1)*cos(beta);
	c2 = -(sin(gamma1)*sin(alpha)+cos(gamma1)*sin(beta)*cos(alpha));
	d2 = cos(gamma1)*cos(alpha) + sin(gamma1)*sin(beta)*sin(alpha);
	e2 = sin(gamma1)*cos(beta);
	f2 = cos(gamma1)*sin(alpha) - sin(gamma1)*sin(beta)*cos(alpha);

	Double_t a1,b1,c1,d1,e1,f1;
	a1 = cos(alpha)*cos(alpha) + sin(beta)*sin(beta) * sin(alpha)*sin(alpha);
	b1 = cos(beta)* cos(beta);
	c1 = sin(alpha)*sin(alpha) + sin(beta)*sin(beta) * cos(alpha)*cos(alpha);      
	d1 = 2*cos(beta)*sin(beta)*sin(alpha);
	e1 = -2*cos(beta)*sin(beta)*cos(alpha);
	f1 = 2*cos(beta)*cos(beta)*sin(alpha)*cos(alpha);

	Double_t Px1,Py1,Pz1,Px2,Py2,Pz2;
	Px1 = P_h*sin(theta_h)*(-sin(phi_h));
	Py1 = P_h*sin(theta_h)*cos(phi_h);
	Pz1 = 0.;

	Px2 = -P_h*cos(phi_h)/tan(theta_h);
	Py2 = -P_h*sin(phi_h)/tan(theta_h);
	Pz2 = P_h*qgamma;

	Double_t Px3 = sin(theta_h)*cos(phi_h);
	Double_t Py3 = sin(theta_h)*sin(phi_h);
	Double_t Pz3 = qgamma*(cos(theta_h) -qbeta/sqrt(P_h*P_h+mass_pi*mass_pi)*P_h);

	Double_t P_hyp,P_hxp;
	P_hyp = a2*P_hx + b2*P_hy + c2*P_hz;
	P_hxp = d2*P_hx + e2*P_hy + f2*P_hz;


	//cout << a1*P_hx*P_hx + b1*P_hy*P_hy + c1*P_hz*P_hz + d1*P_hx*P_hy + e1 *P_hy*P_hz + f1 * P_hx*P_hz << "\t" << atan(P_hyp/P_hxp) << endl;

	Double_t jaco_F;
	Double_t jaco2; // py/pa
	Double_t jaco1; // px/pb
	Double_t jaco3; // pz/pd
	Double_t jaco4; // pl/pc
	Double_t jaco5; // pm/pe
	Double_t jaco6; // pm/pf
	Double_t jaco7; // pn/pe
	Double_t jaco8; // pn/pf

	Double_t jaco9;  // px/pa
	Double_t jaco10; // py/pb
	Double_t jaco11; // pz/pe
	Double_t jaco12; // pm/pd
	Double_t jaco13; // pn/pd

	jaco9 = (P_ei*(1.-cos(theta_ef))*(P_ei*(energy_pro+mom_pro)-P_ef*(energy_pro+mom_pro*cos(theta_ef)))+P_ef*P_ei*(1-cos(theta_ef))*(energy_pro+mom_pro*cos(theta_ef)))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef)))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // px / pa

	jaco1 = (-P_ef*P_ei*(P_ei*(mom_pro+energy_pro)-P_ef*(energy_pro+mom_pro*cos(theta_ef)))+P_ei*P_ef*(1-cos(theta_ef))*P_ef*mom_pro)/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef)))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // p x / pb



	jaco2 = - (energy_pro+mom_pro*cos(theta_ef))/(P_ei*(energy_pro+mom_pro)); // py/pa

	jaco10 = (-1*P_ef*mom_pro)/(P_ei*(energy_pro+mom_pro)); // py / pb

	jaco3 = (P_h*energy_pro/sqrt(P_h*P_h + mass_pi*mass_pi) + mom_pro*cos(theta_h))/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // p z / p d
	jaco11 = mom_pro*P_h/(P_ei*(energy_pro+mom_pro) - P_ef*(energy_pro+mom_pro*cos(theta_ef))); // p z/ p e

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



void dxs(Double_t x, Double_t y, Double_t Q2, Double_t epsilon, Double_t gamma, Double_t z, Double_t pt, Double_t jacoF, Int_t target_flag, Int_t particle_flag,Double_t bpt_m, Double_t bpt_p, Double_t* dxs_hp, Double_t* dxs_hm,Double_t* dxs_all){
	double qu=2./3.;
	double qd=-1./3.;
	double qs=-1./3.;

	*dxs_hp = 1./137.035/137.035/x/y/Q2/1000000. *y*y/(2*(1-epsilon))*(1+gamma*gamma/2./x)*(197.3*197.3/100.*1.0e9);
	*dxs_hm = *dxs_hp;

	//cout << *dxs_hp << "\t";

	*dxs_hp = (*dxs_hp)*bpt_p/PI*exp(-bpt_p*pt*pt)*x;
	*dxs_hm = (*dxs_hm)*bpt_m/PI*exp(-bpt_m*pt*pt)*x;

	//cout << *dxs_hp << "\t";

	*dxs_hp = *dxs_hp*jacoF;
	*dxs_hm = *dxs_hm*jacoF;

	Double_t df_p_hp=0,df_p_hm=0;
	Double_t df_n_hp=0,df_n_hm=0;

	Double_t uquark,dquark,squark,ubarquark,dbarquark,sbarquark;
	Double_t D_fav,D_unfav,D_s,D_g;
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



void dxs_hpt(Double_t x, Double_t y, Double_t Q2, Double_t z, Double_t pt, Double_t phi_h, Double_t jacoF, Int_t target_flag, Int_t particle_flag, Double_t* dxs_hp, Double_t* dxs_hm,Double_t* dxs_all){

	Double_t alpha_s;
	alpha_s = alphasPDF(sqrt(Q2));


	//cout << y << "\t" << Q2 << endl;
	//cout << alpha_s << endl;

	*dxs_hp = 1./137.035/137.035/16./PI/PI/Q2/Q2*y/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.;
	*dxs_hm = 1./137.035/137.035/16./PI/PI/Q2/Q2*y/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.;

	double qu=2./3.;
	double qd=-1./3.;
	double qs=-1./3.;

	Double_t Pqq, Pqg, Pgq;
	Double_t zp,xp;
	Double_t df_p_hp = 0,df_p_hm = 0,df_n_hp = 0,df_n_hm = 0;

	Double_t uquark,dquark,squark,ubarquark,dbarquark,sbarquark,gluon;
	Double_t D_fav,D_unfav,D_s,D_g;

	// doing the integral
	for (Int_t i=0;i!=100;i++){
		xp = x + (1-x)/100.*(i+0.5);
		zp =  z / (z + xp*pt*pt/z/(1-xp)/Q2);


		if (z/zp>1) continue;


		Pqq = 64*PI/3.*Q2/y/y*((1+(1-y)*(1-y))*((1-xp)*(1-zp) + (1+xp*xp*zp*zp)/(1-xp)/(1-zp)) 
				+ 8*xp*zp*(1-y) //- 4 *sqrt(xp*zp*(1-y)/(1-xp)/(1-zp))*(2-y)*(xp*zp+(1-xp)*(1-zp))*cos(phi_h)
				//+ 4*xp*zp*(1-y)*cos(2. * phi_h)
				) *(1-x)/100. /(xp*pt*pt+z*z*(1-xp)*Q2); //13

		Pqg = 64*PI/3.*Q2/y/y*( (1+(1-y)*(1-y))*((1-xp)*zp + (1+xp*xp*(1-zp)*(1-zp))/(1-xp)/zp) 
				+ 8*xp*(1-y)*(1-zp) 
				//+ 4.*sqrt(xp*(1-y)*(1-zp)/(1-xp)/zp)*(2-y)*(xp*(1-zp)+(1-xp)*zp)*cos(phi_h)
				//+ 4.*xp*(1-y)*(1-zp)*cos(2.*phi_h)
				)*(1-x)/100. /(xp*pt*pt+z*z*(1-xp)*Q2);  //14

		Pgq = 8.*PI*Q2/y/y*((1+(1-y)*(1-y))*(xp*xp+(1-xp)*(1-xp))*(zp*zp+(1-zp)*(1-zp))/zp/(1-zp) 
				+ 16*xp*(1-xp)*(1-y)
				//-4.*sqrt(xp*(1-xp)*(1-y)/zp/(1-zp))*(2-y)*(1-2*xp)*(1-2.*zp)*cos(phi_h)
				//+ 8*xp*(1-xp)*(1-y)*cos(2.*phi_h)
				)*(1-x)/100. /(xp*pt*pt+z*z*(1-xp)*Q2);  //15

		//cout << xp << "\t" << zp << "\t" << z/zp << "\t" << x/xp << "\t" << Pqq << "\t" << Pqg << "\t" << Pgq <<"\t" << (1-x)/100./(xp*pt*pt+z*z*(1-xp)*Q2)<<  endl;
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
