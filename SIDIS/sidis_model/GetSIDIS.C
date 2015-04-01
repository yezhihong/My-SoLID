//////////////////////////////////////////////////////
// SIDIS Events Generators for SoLID or EIC         //
//                                                  //
//Note: Basically the same as Xin Qian's "collider" //
//      but the model is coded in "SIDIS.h"         //
//  -- Zhihong Ye, 06/10/2014                       //
//////////////////////////////////////////////////////
#include <TROOT.h>
#include "TApplication.h"
#include "Rtypes.h"
#include "math.h"
#include "iostream"
#include "fstream"
#include "TFile.h"
#include "TTree.h"
#include "SIDIS.h"

const double mm = 0.001;//default is cm
const double cm = 1.0;//default is cm
const double m = 1000.0;//default is cm
//char *LHAPDF_path = getenv("LHAPDF");
int main(Int_t argc, char *argv[]){

	gRandom->SetSeed(0);// uses time for seed

	if (argc != 8) {
		cout << "usage: ./GetSIDIS #ele_mom #ion_mom #target #particle #fileno #events config" << endl;
		cout << "both mom are in GeV, set ion_mom=0 for fix target" << endl;
		cout << "#target_flag = 1 for proton and 2 for deuteron 3 for 3he" << endl;
		cout << "#particle_flag = 1 for pion+, -1 for pion- and 2 for kaon+, -2 for kaon-" << endl;
		cout << "#fileno is the file number of output, used for batch" << endl;
		cout << "#events is number of event in each file" << endl;
		cout << "config is 'EIC' or 'SoLID' which needs ion_mom=0" << endl;
		return 0;
	}
	else{

		ofstream test_log("test_log.dat");

		/*Inputs&Output{{{*/
		//initialize
		double Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon, jacoF;

		Double_t momentum_ele = atof(argv[1]);
		Double_t momentum_ion = atof(argv[2]);
		Int_t target_flag = atoi(argv[3]);
		Int_t particle_flag = atoi(argv[4]);
		Int_t fileno = atoi(argv[5]);
		Int_t number_of_events = atoi(argv[6]);
		string config = argv[7];

		if (config != "EIC" && config != "SoLID"){
			cout << "not supported config" << endl;
			return 0;
		}    

		//define output file
		TString filename,prefix;
		prefix="sidis";

		if (target_flag ==1){
			prefix += "_p";
		}else if (target_flag == 2){
			prefix += "_d";
			momentum_ion = momentum_ion/2.;
		}else if (target_flag == 3){
			prefix += "_3he";
		}else{
			cout << "target_flag is wrong! 1,2,3 " << endl;
			return 0;
		}

		Double_t charge_had = 0;
		int pid_had = 0;
		double mass_had = 0.0;
		if (particle_flag == 1){
			prefix += "_pip";
			charge_had = 1; pid_had = 211; mass_had = 139.57/1000.;//GeV
		}else if (particle_flag == -1){
			prefix += "_pim";
			charge_had = -1; pid_had = -211; mass_had = 139.57/1000.;//GeV
		}else if (particle_flag == 2){
			prefix += "_kp";
			charge_had = 1; pid_had = 321; mass_had = 493.68/1000.;//GeV
		}else if (particle_flag == -2){
			prefix += "_km";
			charge_had = -1; pid_had = -321; mass_had = 493.68/1000.;//GeV
		}else{
			cout << "particle_flag is wrong +-1 and +-2" << endl;
			return 0;
		}

		ofstream outgemc;
		outgemc.open(Form("%s_%d_%d_%d.LUND",prefix.Data(),Int_t(momentum_ele),Int_t(momentum_ion),Int_t(fileno)));

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

		test_log << "Test A ..."<<endl;

		/*Define{{{*/
		Double_t theta_gen= 0.0 , phi_gen = 0.0, mom_gen = 0.0;

		Int_t count[4] = {0,0,0,0};

		Double_t mom_gen_ele,mom_gen_had;
		Double_t theta_gen_ele,theta_gen_had;
		Double_t phi_gen_ele,phi_gen_had;
		Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
		Double_t cdxs,dxs_hm,dxs_hp,dilute_hp,dilute_hm;
		Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E0_ele,E0_had;
		Int_t nsim = 0;

		//For Beam Position and Vertex info
		const double raster_x_size = 0.5 * cm;
		const double raster_y_size = 0.5 * cm;
		const double target_center = -350.0 * cm;
		const double target_length = 40.0 *cm;
		Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;

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
		t1->Branch("cdxs",&cdxs,"cdxs/I");
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
		t1->Branch("dilute_p",&dilute_hp,"data/D");
		t1->Branch("dilute_m",&dilute_hm ,"data/D");
		t1->Branch("px_ele",&px_ele, "px_ele/D");
		t1->Branch("py_ele",&py_ele, "py_ele/D");
		t1->Branch("pz_ele",&pz_ele, "pz_ele/D");
		t1->Branch("E0_ele",&E0_ele, "E0_ele/D");
		t1->Branch("px_had",&px_had, "px_had/D");
		t1->Branch("py_had",&py_had, "py_had/D");
		t1->Branch("pz_had",&pz_had, "pz_had/D");
		t1->Branch("E0_had",&E0_had, "E0_had/D");
		t1->Branch("vx_ele",&vx_ele, "vx_ele/D");
		t1->Branch("vy_ele",&vy_ele, "vy_ele/D");
		t1->Branch("vz_ele",&vz_ele, "vz_ele/D");
		t1->Branch("vx_had",&vx_had, "vx_had/D");
		t1->Branch("vy_had",&vy_had, "vy_had/D");
		t1->Branch("vz_had",&vz_had, "vz_had/D");


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
		t2->Branch("cdxs",&cdxs,"cdxs/I");
		t2->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t2->Branch("mom_had",&mom_had,"mom_had/D");
		t2->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t2->Branch("theta_had",&theta_had,"theta_had/D");
		t2->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t2->Branch("phi_had",&phi_had,"phi_had/D");
		t2->Branch("nsim",&nsim,"nsim/I");
		t2->Branch("dilute_p",&dilute_hp ,"data/D");
		t2->Branch("dilute_m",&dilute_hm ,"data/D");
		t2->Branch("px_ele",&px_ele, "px_ele/D");
		t2->Branch("py_ele",&py_ele, "py_ele/D");
		t2->Branch("pz_ele",&pz_ele, "pz_ele/D");
		t2->Branch("E0_ele",&E0_ele, "E0_ele/D");
		t2->Branch("px_had",&px_had, "px_had/D");
		t2->Branch("py_had",&py_had, "py_had/D");
		t2->Branch("pz_had",&pz_had, "pz_had/D");
		t2->Branch("E0_had",&E0_had, "E0_had/D");
		t2->Branch("vx_ele",&vx_ele, "vx_ele/D");
		t2->Branch("vy_ele",&vy_ele, "vy_ele/D");
		t2->Branch("vz_ele",&vz_ele, "vz_ele/D");
		t2->Branch("vx_had",&vx_had, "vx_had/D");
		t2->Branch("vy_had",&vy_had, "vy_had/D");
		t2->Branch("vz_had",&vz_had, "vz_had/D");


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
		t3->Branch("cdxs",&cdxs,"cdxs/I");
		t3->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t3->Branch("mom_had",&mom_had,"mom_had/D");
		t3->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t3->Branch("theta_had",&theta_had,"theta_had/D");
		t3->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t3->Branch("phi_had",&phi_had,"phi_had/D");
		t3->Branch("nsim",&nsim,"nsim/I");
		t3->Branch("dilute_p",&dilute_hp ,"data/D");
		t3->Branch("dilute_m",&dilute_hm ,"data/D");
		t3->Branch("px_ele",&px_ele, "px_ele/D");
		t3->Branch("py_ele",&py_ele, "py_ele/D");
		t3->Branch("pz_ele",&pz_ele, "pz_ele/D");
		t3->Branch("E0_ele",&E0_ele, "E0_ele/D");
		t3->Branch("px_had",&px_had, "px_had/D");
		t3->Branch("py_had",&py_had, "py_had/D");
		t3->Branch("pz_had",&pz_had, "pz_had/D");
		t3->Branch("E0_had",&E0_had, "E0_had/D");
		t3->Branch("vx_ele",&vx_ele, "vx_ele/D");
		t3->Branch("vy_ele",&vy_ele, "vy_ele/D");
		t3->Branch("vz_ele",&vz_ele, "vz_ele/D");
		t3->Branch("vx_had",&vx_had, "vx_had/D");
		t3->Branch("vy_had",&vy_had, "vy_had/D");
		t3->Branch("vz_had",&vz_had, "vz_had/D");


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
		t4->Branch("cdxs",&cdxs,"cdxs/I");
		t4->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t4->Branch("mom_had",&mom_had,"mom_had/D");
		t4->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t4->Branch("theta_had",&theta_had,"theta_had/D");
		t4->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t4->Branch("phi_had",&phi_had,"phi_had/D");
		t4->Branch("nsim",&nsim,"nsim/I");
		t4->Branch("dilute_p",&dilute_hp ,"data/D");
		t4->Branch("dilute_m",&dilute_hm ,"data/D");
		t4->Branch("px_ele",&px_ele, "px_ele/D");
		t4->Branch("py_ele",&py_ele, "py_ele/D");
		t4->Branch("pz_ele",&pz_ele, "pz_ele/D");
		t4->Branch("E0_ele",&E0_ele, "E0_ele/D");
		t4->Branch("px_had",&px_had, "px_had/D");
		t4->Branch("py_had",&py_had, "py_had/D");
		t4->Branch("pz_had",&pz_had, "pz_had/D");
		t4->Branch("E0_had",&E0_had, "E0_had/D");
		t4->Branch("vx_ele",&vx_ele, "vx_ele/D");
		t4->Branch("vy_ele",&vy_ele, "vy_ele/D");
		t4->Branch("vz_ele",&vz_ele, "vz_ele/D");
		t4->Branch("vx_had",&vx_had, "vx_had/D");
		t4->Branch("vy_had",&vy_had, "vy_had/D");
		t4->Branch("vz_had",&vz_had, "vz_had/D");
		/*}}}*/
		
		test_log << "Test B ..."<<endl;

		//Only initialize once here
		SIDIS *sidis = new SIDIS();
		sidis->SetLAPDF();
		//sidis->Print();
		
		test_log << "Test C ..."<<endl;

		bool exitcondition=true;	
		while(exitcondition){
			nsim ++;
		
			test_log << "Test D ..."<<endl;

			/*Generator{{{*/
			//For electron
			vx_ele = gRandom->Uniform(-raster_x_size, raster_x_size);
			vy_ele = gRandom->Uniform(-raster_y_size, raster_y_size);
			vz_ele = target_center + gRandom->Uniform(-target_length/2.0, target_length/2.0);

			phi_gen = gRandom->Uniform(0.,2.*PI);
			if (config=="EIC") {
				theta_gen = acos(gRandom->Uniform(cos(150/DEG),1.));
				mom_gen = gRandom->Uniform(0.7,momentum_ele * 3.);
			}
			else if (config=="SoLID") {
				theta_gen = acos(gRandom->Uniform(cos(30/DEG),cos(7/DEG)));
				mom_gen = gRandom->Uniform(0.5,momentum_ele);
			}
			mom_gen_ele = mom_gen; theta_gen_ele = theta_gen*DEG; phi_gen_ele = phi_gen*DEG;

			//For hadron
			//vx_had = gRandom->Uniform(-raster_x_size, raster_x_size);
			//vy_had = gRandom->Uniform(-raster_y_size, raster_y_size);
			//vz_had = target_center + gRandom->Uniform(-target_length/2.0, target_length/2.0);
			vx_had = vx_ele; vy_had = vy_ele; vz_had = vz_ele;

			theta_gen = acos(gRandom->Uniform(cos(30/DEG),cos(7/DEG)));
			phi_gen = gRandom->Uniform(0.,2.*PI);
			mom_gen = gRandom->Uniform(0.5,6);
			mom_gen_had = mom_gen; theta_gen_had = theta_gen*DEG; phi_gen_had = phi_gen*DEG;
			/*}}}*/
		
			test_log << "Test E ..."<<endl;

			sidis->Init(momentum_ele, momentum_ion,
					mom_gen_ele, theta_gen_ele, phi_gen_ele,
					mom_gen_had, theta_gen_had, phi_gen_had,
					target_flag, particle_flag);
			
			test_log << "Test F ..."<<endl;

			mom_ele = sidis->mom_ele; theta_ele = sidis->theta_ele; phi_ele = sidis->phi_ele;
			mom_had = sidis->mom_had; theta_had = sidis->theta_had; phi_ele = sidis->phi_had;
			theta_q=sidis->theta_q;  theta_s=sidis->theta_s; phi_s= sidis->phi_s; phi_h = sidis->phi_h;
			px_ele = sidis->px_ele;	py_ele = sidis->py_ele;	pz_ele = sidis->pz_ele; E0_ele = sidis->E0_ele;
			px_had = sidis->px_had;	py_had = sidis->py_had;	pz_had = sidis->pz_had; E0_had = sidis->E0_had;

			x=sidis->x; y=sidis->y; z=sidis->z; Q2=sidis->Q2; W=sidis->W; Wp=sidis->Wp;
			s=sidis->s; nu=sidis->nu; pt=sidis->pt; gamma=sidis->gamma; epsilon=sidis->epsilon;
			jacoF=sidis->jacoF;
			
			test_log << "Test G ..."<<endl;

			/*Get XS{{{*/
			if (Q2 >=1.0 && W>= 2.3 &&Wp>= 1.6 &&
					( (config=="EIC" && z>0.2&&z<0.8&&y>0.05&&y<0.8
					   && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) 
						   || (count[1]<number_of_events&&pt>1.0&&Q2<=10.)
						   || (count[2]<number_of_events&&pt<=1.0&&Q2>10.)
						   || (count[3]<number_of_events&&pt>1.0&&Q2>10.)) )
					  ||(config=="SoLID" && z>0.3&&z<0.7 
						  && ((count[0]<number_of_events&&pt<=1.0&&Q2<=10.) 
							  || (count[1]<number_of_events&&pt>1.0&&Q2<=10.))) )){


				sidis->CalcXS();
				dxs_hp = sidis->GetXS_HP();
				dxs_hm = sidis->GetXS_HM();
				dilute_hp = sidis->GetDilute_HP();
				dilute_hm = sidis->GetDilute_HM();
				if (particle_flag == 1||particle_flag == 2)
					cdxs = dxs_hp;
				else if (particle_flag ==-1||particle_flag ==-2)
					cdxs = dxs_hm;
				else
					cdxs = 0.0;


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
				cout << count[0] << "\t" << count[1] << "\t" << count[2] << "\t" << count[3] << "\r";
				//cout << nsim << endl;

				/*Output as text format for GEMC{{{*/

				if ((dxs_hp+dxs_hm)!=0){
					//Header:      1#part. 2#x 3#z 4#pt 5#Pol 6#sxs_e 7#sxs_h 8#cxs 9#phi_s 10#phi_h
					outgemc<<Form("    %2d \t %10.4e \t %10.4e \t %10.4e \t %4.3f \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
							2, //ele+had
							x,
							z,
							pt,
							1.0, //pol = 1.0 for now
							Q2,
							W,
							cdxs,
							phi_s,
							phi_h
							)<<endl;

					//electron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
					outgemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
							1, //index
							-1.0,//charge
							1, //=1 for active 
							11,//pid
							0,// parent pid, not in used now
							0,// doughter for decay bookkeeping, not in used now
							px_ele, //GeV
							py_ele, //GeV
							pz_ele, //GeV
							E0_ele, //GeV
							0.0005, //mass not in used	
							vx_ele, //vx,cm
							vy_ele, //vy,cm
							vz_ele  //vz,cm
							)<<endl;
					//hadron info: 1#index. 2#charge 3#type 4#pid 5#mpid 6#daughter 7#px 8#py 9#pz 10#E 11#mass 12#vx 13#vy 14#vz
					outgemc<<Form("%2d \t %4.2f \t %1d \t %8d \t %1d \t %1d \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e \t %10.4e",
							2, //index
							charge_had,//charge
							1, //=1 for active 
							pid_had,//pid
							0,// parent pid, not in used now
							0,// doughter for decay bookkeeping, not in used now
							px_had, //GeV
							py_had, //GeV
							pz_had, //GeV
							E0_had, //GeV
							mass_had, //mass not in used
							vx_had, //vx,cm
							vy_had, //vy,cm
							vz_had  //vz,cm
							)<<endl;
				}
				/*}}}*/

			}
			/*}}}*/
			
			test_log << "Test H ..."<<endl;

			//judging exitcondition
			if (config=="EIC") {
				if (count[0] < number_of_events || count[1] < number_of_events 
						|| count[2] < number_of_events || count[3] < number_of_events) exitcondition=true;
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
		delete sidis;
		test_log.close();
	}

	return 0;
} 
