//////////////////////////////////////////////////////////
// This script works on GEMC 2.0 only, because
// GEMC2.0 and GEMC1.0 have different TTree and branch names
// Zhihong Ye, 09/12/2014
//
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3F.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

//void ExtractGEM(string input_filename)
Int_t main()
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const Double_t DEG=180./3.1415926;
	string input_filename = "eDIS";

	/*Tree and branches{{{*/
	TChain *header = new TChain("header");
	TChain *generated = new TChain("generated");
	TChain *flux = new TChain("flux");
	TString ROOT_Dir = "/v/volatile/halla/solid/yez/sidis/gemc2/3he_pip_11/";
	for(Int_t i=0;i<200;i++){
		header->Add(Form("%s/SIDIS_3he_pip_11_%d.root", ROOT_Dir.Data(), i));
		generated->Add(Form("%s/SIDIS_3he_pip_11_%d.root", ROOT_Dir.Data(), i));
		flux->Add(Form("%s/SIDIS_3he_pip_11_%d.root", ROOT_Dir.Data(), i));
	}
	
	/*
	TChain *header = new TChain("header");
	TChain *generated = new TChain("generated");
	TChain *flux = new TChain("flux");
	header->Add(input_filename.c_str());
	generated->Add(input_filename.c_str());
	flux->Add(input_filename.c_str());
*/

	/*
	TFile *file=new TFile(input_filename.c_str());
	if (file->IsZombie()) {
		cout << "Error opening file" << input_filename << endl;
		//continue;
		exit(-1);
	}
	else cout << "### Open file " << input_filename << endl;
*/
	
	/*Temp Rate Normalization{{{*/
	Double_t filenum=1.;
	if (input_filename.find("_EM_",0) != string::npos) {
		cout << "### EM background from beam on target" <<  endl;
	}else if (input_filename.find("_clean_weighted_",0) != string::npos) {
		cout << "### Background from weighted event generator with no interaction except decay" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "### filenum = " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "### we need filenum for addtional normalization" << endl; return -1;}		
	}	  
	else if (input_filename.find("_dirty_normalized_",0) != string::npos) {
		cout << "### background from normalized event generator with all interaction" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "### filenum = " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "### we need filenum for addtional normalization" << endl; return -2;}
	}	  	
	else {cout << "### Not EM or clean_weighted or dirty_normalized, so I will use filenum=1. " << endl; filenum=1.0;}
    /*}}}*/

	/*Set Branch{{{*/
	//Header Tree:
	// Var#1~#8 are free slots for propogating important info from the "INPUT generator seed"
	// For example, they can be used to store the cross section and other physics quantities
	// In eicRate, we store the following quantities:
	// var1->Wprate, var2->Wmrate, var3->targetPol, var4->x,var5->y, var6->W, var7->Q2, var8->rate 
	//
	//TTree *header = (TTree*) file->Get("header");
	vector <Double_t> *head_evn=0,*head_evn_type=0; //Note: Vectors have to be initialized at first!!!
	vector <Double_t> *head_beamPol=0;
	vector<Double_t> *head_Wmrate=0, *head_Wprate=0, *head_targetPol=0, *head_x=0, *head_Q2=0, *head_W=0, *head_rate=0, *head_y=0;
	header->SetBranchAddress("evn",&head_evn);
	header->SetBranchAddress("evn_type",&head_evn_type);
	header->SetBranchAddress("beamPol",&head_beamPol);
	header->SetBranchAddress("var1",    &head_Wprate);
	header->SetBranchAddress("var2",    &head_Wmrate);
	header->SetBranchAddress("var3",    &head_targetPol);
	header->SetBranchAddress("var4",    &head_x);
	header->SetBranchAddress("var5",    &head_y);
	header->SetBranchAddress("var6",    &head_W);
	header->SetBranchAddress("var7",    &head_Q2);
	header->SetBranchAddress("var8",    &head_rate);

	//TTree *generated = (TTree*) file->Get("generated");
	vector <Int_t> *gen_pid=0;
	vector <Double_t> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	generated->SetBranchAddress("pid",&gen_pid);
	generated->SetBranchAddress("px",&gen_px);
	generated->SetBranchAddress("py",&gen_py);
	generated->SetBranchAddress("pz",&gen_pz);
	generated->SetBranchAddress("vx",&gen_vx);
	generated->SetBranchAddress("vy",&gen_vy);
	generated->SetBranchAddress("vz",&gen_vz);

	//TTree *flux = (TTree*) file->Get("flux");
	vector<Double_t> *flux_id=0,*flux_hitn=0,*flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<Double_t> *flux_trackE=0,*flux_totEdep=0;
	vector<Double_t> *flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0;
	vector<Double_t> *flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
	flux->SetBranchAddress("hitn",&flux_hitn); //number of flux hits in one event
	flux->SetBranchAddress("id",&flux_id); //detector ID
	flux->SetBranchAddress("pid",&flux_pid); //particle ID
	flux->SetBranchAddress("mpid",&flux_mpid);  //mother particle ID
	flux->SetBranchAddress("tid",&flux_tid); //flux track ID
	flux->SetBranchAddress("mtid",&flux_mtid); //mother flux track ID
	flux->SetBranchAddress("otid",&flux_otid);//orginal flux track ID ?
	flux->SetBranchAddress("trackE",&flux_trackE);//energy of the flux particle
	flux->SetBranchAddress("totEdep",&flux_totEdep); //deposited energy of the flux particle
	flux->SetBranchAddress("avg_x",&flux_avg_x);
	flux->SetBranchAddress("avg_y",&flux_avg_y);
	flux->SetBranchAddress("avg_z",&flux_avg_z);
	flux->SetBranchAddress("avg_lx",&flux_avg_lx);
	flux->SetBranchAddress("avg_ly",&flux_avg_ly);
	flux->SetBranchAddress("avg_lz",&flux_avg_lz);
	flux->SetBranchAddress("avg_t",&flux_avg_t);
	flux->SetBranchAddress("px",&flux_px);
	flux->SetBranchAddress("py",&flux_py);
	flux->SetBranchAddress("pz",&flux_pz);
	flux->SetBranchAddress("vx",&flux_vx);
	flux->SetBranchAddress("vy",&flux_vy);
	flux->SetBranchAddress("vz",&flux_vz);
	flux->SetBranchAddress("mvx",&flux_mvx);
	flux->SetBranchAddress("mvy",&flux_mvy);
	flux->SetBranchAddress("mvz",&flux_mvz);
	/*End Set Branch}}}*/
	Int_t nevent = (Int_t)generated->GetEntries();
	cout << "nevent = " << nevent << endl;


	/*End Set Branch}}}*/

	/*Define GEM{{{*/
	//In each segmentation, 40 slides, each slides has 2.5cm width, 0.3cm between two slides.
//	const unsigned int GEM_Layer = 23; //30 module around the circle
	const unsigned int GEM_Plane = 6;// put 6 slides in each module just for check the R-dependence, 
	const Double_t GEM_R_Min[GEM_Plane] = {36,21,25,32,42,55};//cm, GEM inner radius 
	const Double_t GEM_R_Max[GEM_Plane] = {87,98,112,135,100,123};//cm, GEM outer radius
	const Double_t GEM_Z_Set[GEM_Plane] = {-175,-150,-119,-68,5,92}; //cm, GEM Z position
	const Int_t GEM_ID[GEM_Plane] = {1120000,1220000,1320000,1420000,1520000,1620000}; //GEM #3rd layer
//	const Double_t GEM_VP_Z[GEM_Plane] = {-175-0.5,-150-0.5,-119-0.5,-68-0.5,5-0.5,92-0.5}; //GEM virtual plane
//	const Int_t GEM_VP_ID[GEM_Plane] = {1110000,1210000,1310000,1410000,1510000,1610000,}; //GEM virtual plane

	Double_t GEM_X[GEM_Plane], GEM_Y[GEM_Plane], GEM_Z[GEM_Plane];
	Double_t GEM_Px[GEM_Plane], GEM_Py[GEM_Plane], GEM_Pz[GEM_Plane];
	Double_t GEM_Mom[GEM_Plane], GEM_Theta[GEM_Plane], GEM_Phi[GEM_Plane];
	Double_t GEM_E[GEM_Plane], GEM_Edep[GEM_Plane];
	Int_t GEM_ID_Pick[GEM_Plane], GEM_PID_Pick[GEM_Plane], GEM_TID_Pick[GEM_Plane];
	//I just want to initialize everything
	for(unsigned int ip=0;ip<GEM_Plane;ip++){
		GEM_X[ip] = 0.0; GEM_Y[ip] = 0.0; GEM_Z[ip] = 0.0;
		GEM_Px[ip] = 0.0; GEM_Py[ip] = 0.0; GEM_Pz[ip] = 0.0;
		GEM_Mom[ip] = 0.0; GEM_Theta[ip] = 0.0; GEM_Phi[ip] = 0.0;
		GEM_E[ip] = 0.0; GEM_Edep[ip]=0.;
	}

	//Note: Use EC's Virtual Plane to determine the positions and momentum info on EC
	const Int_t LAEC_VP_ID = 3210000;
	const Double_t LAEC_VP_Rmin = 83.0, LAEC_VP_Rmax = 140.0, LAEC_VP_Z = -65.0;//cm
	const Int_t FAEC_VP_ID = 3110000;
	const Double_t FAEC_VP_Rmin = 105.0, FAEC_VP_Rmax = 235.0, FAEC_VP_Z = 427.5;//cm

	//Other Definition	
	/*Particle Definition{{{*/	
	const Int_t Electron = 11;
	//const Int_t Photon = 22;
	const Int_t Proton = 2212;
	//const Int_t Pi0 = 111;
	const Int_t Pion = 211;
	//const Int_t K0 = 311;
	const Int_t Kaon = 321;
	//const Int_t Beam = 0;
	//const Int_t Neutron = 2112;
	//const Int_t Neutrino1 = 12;//Nu_e
	//const Int_t Neutrino2 = 14;//Nu_Mu
	//const Int_t Neutrino3 = 16;//Nu_Tao
	/*}}}*/

	const Double_t MeV2GeV = 1./1000.;
	const Double_t MM2CM = 1./10.;
	/*}}}*/

	/*Tree Define{{{*/
	Double_t rate, mom_gen, theta_gen,phi_gen,mom_flux,theta_flux,phi_flux,r_flux,x_flux,y_flux,z_flux,R_EC;
	Double_t px_flux, py_flux,pz_flux,vx_flux, vy_flux,vz_flux,Edep, E;
	Double_t EC_Cut_G, EC_Cut_E, EC_Cut_P;
	Double_t px_gen, py_gen, pz_gen, vx_gen, vy_gen, vz_gen;
	Int_t evn, evn_flux, ID_flux,PID_flux,ID_Pick,PID_Pick, TID_flux;

	TString pid_temp = "";
	/*Find PID{{{*/
	PID_Pick = -1000;
	if(input_filename.find("_pi0_",0) != string::npos) {
		PID_Pick = Electron; //Pi0 decays into Electron+Position
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pi0_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pi0_down";
		else
			pid_temp = "pi0";
	}
	else if(input_filename.find("_pip_",0) != string::npos) {
		PID_Pick = Pion;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pip_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pip_down";
		else
			pid_temp = "pip";
	}
	else if(input_filename.find("_pim_",0) != string::npos) {
		PID_Pick = -Pion;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pim_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pim_down";
		else
			pid_temp = "pim";
	}
	else if(input_filename.find("_kp_",0) != string::npos) {
		PID_Pick = Kaon;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "kp_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "kp_down";
		else
			pid_temp = "kp";
	}
	else if(input_filename.find("_km_",0) != string::npos) {
		PID_Pick = -Kaon;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "km_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "km_down";
		else
			pid_temp = "km";
	}
	else if(input_filename.find("_p_",0) != string::npos) {
		PID_Pick = Proton;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "p_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "p_down";
		else
			pid_temp = "p";
	}
	else if(input_filename.find("_eDIS_",0) != string::npos) {
		PID_Pick = Electron;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "eDIS_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "eDIS_down";
		else	
			pid_temp = "eDIS";
	}
	else if(input_filename.find("_EM_",0) != string::npos) {
		PID_Pick = Electron;
		pid_temp = "EM";
	}
	else{
		cerr<<"****Attention, I don't understand the file name, so I assume to find electrons"<<endl;	
		PID_Pick = Electron;
		pid_temp = "Real";
	};
	/*}}}*/
	
	TString output_filename = Form("GEM_SIDIS_He3_%s.root",pid_temp.Data() );
	cerr<<"--- I will save the detector info into "<<output_filename.Data()<<endl;
	
	TFile *out_file = new TFile(output_filename.Data(),"recreate");

	/*GEM-Plane-1{{{*/
	TTree *G1 = new TTree("G1","A new Tree for particles going in GEM#1");
	//From Header Tree
	G1->Branch("rate",&rate,"rate/D");
	G1->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	G1->Branch("mom_gen",&mom_gen,"mom_gen/D");
	G1->Branch("theta_gen",&theta_gen,"theta_gen/D");
	G1->Branch("phi_gen",&phi_gen,"phi_gen/D");
	G1->Branch("px_gen",&px_gen,"px_gen/D");
	G1->Branch("py_gen",&py_gen,"py_gen/D");
	G1->Branch("pz_gen",&pz_gen,"pz_gen/D");
	G1->Branch("vx_gen",&vx_gen,"vx_gen/D");
	G1->Branch("vy_gen",&vy_gen,"vy_gen/D");
	G1->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	G1->Branch("evn_flux",&evn_flux,"evn_flux/I");
	G1->Branch("mom_flux",&mom_flux,"mom_flux/D");
	G1->Branch("theta_flux",&theta_flux,"theta_flux/D");
	G1->Branch("phi_flux",&phi_flux,"phi_flux/D");
	G1->Branch("r_flux",&r_flux,"r_flux/D");
	G1->Branch("x_flux",&x_flux,"x_flux/D");
	G1->Branch("y_flux",&y_flux,"y_flux/D");
	G1->Branch("z_flux",&z_flux,"z_flux/D");
	G1->Branch("vx_flux",&vx_flux,"vx_flux/D");
	G1->Branch("vy_flux",&vy_flux,"vy_flux/D");
	G1->Branch("vz_flux",&vz_flux,"vz_flux/D");
	G1->Branch("px_flux",&px_flux,"px_flux/D");
	G1->Branch("py_flux",&py_flux,"py_flux/D");
	G1->Branch("pz_flux",&pz_flux,"pz_flux/D");
	G1->Branch("E",&E,"E/D");
	G1->Branch("Edep",&Edep,"Edep/D");
	G1->Branch("R_EC",&R_EC,"R_EC/D");
	G1->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	G1->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	G1->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	G1->Branch("ID_flux",&ID_flux,"ID_flux/I");
	G1->Branch("PID_flux",&PID_flux,"PID_flux/I");
	G1->Branch("TID_flux",&TID_flux,"TID_flux/I");
	G1->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	G1->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	
	/*GEM-Plane-2{{{*/
	TTree *G2 = new TTree("G2","A new Tree for particles going in GEM#2");
	//From Header Tree
	G2->Branch("rate",&rate,"rate/D");
	G2->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	G2->Branch("mom_gen",&mom_gen,"mom_gen/D");
	G2->Branch("theta_gen",&theta_gen,"theta_gen/D");
	G2->Branch("phi_gen",&phi_gen,"phi_gen/D");
	G2->Branch("px_gen",&px_gen,"px_gen/D");
	G2->Branch("py_gen",&py_gen,"py_gen/D");
	G2->Branch("pz_gen",&pz_gen,"pz_gen/D");
	G2->Branch("vx_gen",&vx_gen,"vx_gen/D");
	G2->Branch("vy_gen",&vy_gen,"vy_gen/D");
	G2->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	G2->Branch("evn_flux",&evn_flux,"evn_flux/I");
	G2->Branch("mom_flux",&mom_flux,"mom_flux/D");
	G2->Branch("theta_flux",&theta_flux,"theta_flux/D");
	G2->Branch("phi_flux",&phi_flux,"phi_flux/D");
	G2->Branch("r_flux",&r_flux,"r_flux/D");
	G2->Branch("x_flux",&x_flux,"x_flux/D");
	G2->Branch("y_flux",&y_flux,"y_flux/D");
	G2->Branch("z_flux",&z_flux,"z_flux/D");
	G2->Branch("vx_flux",&vx_flux,"vx_flux/D");
	G2->Branch("vy_flux",&vy_flux,"vy_flux/D");
	G2->Branch("vz_flux",&vz_flux,"vz_flux/D");
	G2->Branch("px_flux",&px_flux,"px_flux/D");
	G2->Branch("py_flux",&py_flux,"py_flux/D");
	G2->Branch("pz_flux",&pz_flux,"pz_flux/D");
	G2->Branch("E",&E,"E/D");
	G2->Branch("Edep",&Edep,"Edep/D");
	G2->Branch("R_EC",&R_EC,"R_EC/D");
	G2->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	G2->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	G2->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	G2->Branch("ID_flux",&ID_flux,"ID_flux/I");
	G2->Branch("PID_flux",&PID_flux,"PID_flux/I");
	G2->Branch("TID_flux",&TID_flux,"TID_flux/I");
	G2->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	G2->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/

	/*GEM-Plane-3{{{*/
	TTree *G3 = new TTree("G3","A new Tree for particles going in GEM#3");
	//From Header Tree
	G3->Branch("rate",&rate,"rate/D");
	G3->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	G3->Branch("mom_gen",&mom_gen,"mom_gen/D");
	G3->Branch("theta_gen",&theta_gen,"theta_gen/D");
	G3->Branch("phi_gen",&phi_gen,"phi_gen/D");
	G3->Branch("px_gen",&px_gen,"px_gen/D");
	G3->Branch("py_gen",&py_gen,"py_gen/D");
	G3->Branch("pz_gen",&pz_gen,"pz_gen/D");
	G3->Branch("vx_gen",&vx_gen,"vx_gen/D");
	G3->Branch("vy_gen",&vy_gen,"vy_gen/D");
	G3->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	G3->Branch("evn_flux",&evn_flux,"evn_flux/I");
	G3->Branch("mom_flux",&mom_flux,"mom_flux/D");
	G3->Branch("theta_flux",&theta_flux,"theta_flux/D");
	G3->Branch("phi_flux",&phi_flux,"phi_flux/D");
	G3->Branch("r_flux",&r_flux,"r_flux/D");
	G3->Branch("x_flux",&x_flux,"x_flux/D");
	G3->Branch("y_flux",&y_flux,"y_flux/D");
	G3->Branch("z_flux",&z_flux,"z_flux/D");
	G3->Branch("vx_flux",&vx_flux,"vx_flux/D");
	G3->Branch("vy_flux",&vy_flux,"vy_flux/D");
	G3->Branch("vz_flux",&vz_flux,"vz_flux/D");
	G3->Branch("px_flux",&px_flux,"px_flux/D");
	G3->Branch("py_flux",&py_flux,"py_flux/D");
	G3->Branch("pz_flux",&pz_flux,"pz_flux/D");
	G3->Branch("E",&E,"E/D");
	G3->Branch("Edep",&Edep,"Edep/D");
	G3->Branch("R_EC",&R_EC,"R_EC/D");
	G3->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	G3->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	G3->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	G3->Branch("ID_flux",&ID_flux,"ID_flux/I");
	G3->Branch("PID_flux",&PID_flux,"PID_flux/I");
	G3->Branch("TID_flux",&TID_flux,"TID_flux/I");
	G3->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	G3->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	
	/*GEM-Plane-4{{{*/
	TTree *G4 = new TTree("G4","A new Tree for particles going in GEM#4");
	//From Header Tree
	G4->Branch("rate",&rate,"rate/D");
	G4->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	G4->Branch("mom_gen",&mom_gen,"mom_gen/D");
	G4->Branch("theta_gen",&theta_gen,"theta_gen/D");
	G4->Branch("phi_gen",&phi_gen,"phi_gen/D");
	G4->Branch("px_gen",&px_gen,"px_gen/D");
	G4->Branch("py_gen",&py_gen,"py_gen/D");
	G4->Branch("pz_gen",&pz_gen,"pz_gen/D");
	G4->Branch("vx_gen",&vx_gen,"vx_gen/D");
	G4->Branch("vy_gen",&vy_gen,"vy_gen/D");
	G4->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	G4->Branch("evn_flux",&evn_flux,"evn_flux/I");
	G4->Branch("mom_flux",&mom_flux,"mom_flux/D");
	G4->Branch("theta_flux",&theta_flux,"theta_flux/D");
	G4->Branch("phi_flux",&phi_flux,"phi_flux/D");
	G4->Branch("r_flux",&r_flux,"r_flux/D");
	G4->Branch("x_flux",&x_flux,"x_flux/D");
	G4->Branch("y_flux",&y_flux,"y_flux/D");
	G4->Branch("z_flux",&z_flux,"z_flux/D");
	G4->Branch("vx_flux",&vx_flux,"vx_flux/D");
	G4->Branch("vy_flux",&vy_flux,"vy_flux/D");
	G4->Branch("vz_flux",&vz_flux,"vz_flux/D");
	G4->Branch("px_flux",&px_flux,"px_flux/D");
	G4->Branch("py_flux",&py_flux,"py_flux/D");
	G4->Branch("pz_flux",&pz_flux,"pz_flux/D");
	G4->Branch("E",&E,"E/D");
	G4->Branch("Edep",&Edep,"Edep/D");
	G4->Branch("R_EC",&R_EC,"R_EC/D");
	G4->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	G4->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	G4->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	G4->Branch("ID_flux",&ID_flux,"ID_flux/I");
	G4->Branch("PID_flux",&PID_flux,"PID_flux/I");
	G4->Branch("TID_flux",&TID_flux,"TID_flux/I");
	G4->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	G4->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	
	/*GEM-Plane-5{{{*/
	TTree *G5 = new TTree("G5","A new Tree for particles going in GEM#5");
	//From Header Tree
	G5->Branch("rate",&rate,"rate/D");
	G5->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	G5->Branch("mom_gen",&mom_gen,"mom_gen/D");
	G5->Branch("theta_gen",&theta_gen,"theta_gen/D");
	G5->Branch("phi_gen",&phi_gen,"phi_gen/D");
	G5->Branch("px_gen",&px_gen,"px_gen/D");
	G5->Branch("py_gen",&py_gen,"py_gen/D");
	G5->Branch("pz_gen",&pz_gen,"pz_gen/D");
	G5->Branch("vx_gen",&vx_gen,"vx_gen/D");
	G5->Branch("vy_gen",&vy_gen,"vy_gen/D");
	G5->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	G5->Branch("evn_flux",&evn_flux,"evn_flux/I");
	G5->Branch("mom_flux",&mom_flux,"mom_flux/D");
	G5->Branch("theta_flux",&theta_flux,"theta_flux/D");
	G5->Branch("phi_flux",&phi_flux,"phi_flux/D");
	G5->Branch("r_flux",&r_flux,"r_flux/D");
	G5->Branch("x_flux",&x_flux,"x_flux/D");
	G5->Branch("y_flux",&y_flux,"y_flux/D");
	G5->Branch("z_flux",&z_flux,"z_flux/D");
	G5->Branch("vx_flux",&vx_flux,"vx_flux/D");
	G5->Branch("vy_flux",&vy_flux,"vy_flux/D");
	G5->Branch("vz_flux",&vz_flux,"vz_flux/D");
	G5->Branch("px_flux",&px_flux,"px_flux/D");
	G5->Branch("py_flux",&py_flux,"py_flux/D");
	G5->Branch("pz_flux",&pz_flux,"pz_flux/D");
	G5->Branch("E",&E,"E/D");
	G5->Branch("Edep",&Edep,"Edep/D");
	G5->Branch("R_EC",&R_EC,"R_EC/D");
	G5->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	G5->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	G5->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	G5->Branch("ID_flux",&ID_flux,"ID_flux/I");
	G5->Branch("PID_flux",&PID_flux,"PID_flux/I");
	G5->Branch("TID_flux",&TID_flux,"TID_flux/I");
	G5->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	G5->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	
	/*GEM-Plane-6{{{*/
	TTree *G6 = new TTree("G6","A new Tree for particles going in GEM#6");
	//From Header Tree
	G6->Branch("rate",&rate,"rate/D");
	G6->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	G6->Branch("mom_gen",&mom_gen,"mom_gen/D");
	G6->Branch("theta_gen",&theta_gen,"theta_gen/D");
	G6->Branch("phi_gen",&phi_gen,"phi_gen/D");
	G6->Branch("px_gen",&px_gen,"px_gen/D");
	G6->Branch("py_gen",&py_gen,"py_gen/D");
	G6->Branch("pz_gen",&pz_gen,"pz_gen/D");
	G6->Branch("vx_gen",&vx_gen,"vx_gen/D");
	G6->Branch("vy_gen",&vy_gen,"vy_gen/D");
	G6->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	G6->Branch("evn_flux",&evn_flux,"evn_flux/I");
	G6->Branch("mom_flux",&mom_flux,"mom_flux/D");
	G6->Branch("theta_flux",&theta_flux,"theta_flux/D");
	G6->Branch("phi_flux",&phi_flux,"phi_flux/D");
	G6->Branch("r_flux",&r_flux,"r_flux/D");
	G6->Branch("x_flux",&x_flux,"x_flux/D");
	G6->Branch("y_flux",&y_flux,"y_flux/D");
	G6->Branch("z_flux",&z_flux,"z_flux/D");
	G6->Branch("vx_flux",&vx_flux,"vx_flux/D");
	G6->Branch("vy_flux",&vy_flux,"vy_flux/D");
	G6->Branch("vz_flux",&vz_flux,"vz_flux/D");
	G6->Branch("px_flux",&px_flux,"px_flux/D");
	G6->Branch("py_flux",&py_flux,"py_flux/D");
	G6->Branch("pz_flux",&pz_flux,"pz_flux/D");
	G6->Branch("E",&E,"E/D");
	G6->Branch("Edep",&Edep,"Edep/D");
	G6->Branch("R_EC",&R_EC,"R_EC/D");
	G6->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	G6->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	G6->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	G6->Branch("ID_flux",&ID_flux,"ID_flux/I");
	G6->Branch("PID_flux",&PID_flux,"PID_flux/I");
	G6->Branch("TID_flux",&TID_flux,"TID_flux/I");
	G6->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	G6->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	
	/*LAEC{{{*/
	TTree *LAEC = new TTree("LAEC","A new Tree for particles going in LAEC");
	//From Header Tree
	LAEC->Branch("rate",&rate,"rate/D");
	LAEC->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	LAEC->Branch("mom_gen",&mom_gen,"mom_gen/D");
	LAEC->Branch("theta_gen",&theta_gen,"theta_gen/D");
	LAEC->Branch("phi_gen",&phi_gen,"phi_gen/D");
	LAEC->Branch("px_gen",&px_gen,"px_gen/D");
	LAEC->Branch("py_gen",&py_gen,"py_gen/D");
	LAEC->Branch("pz_gen",&pz_gen,"pz_gen/D");
	LAEC->Branch("vx_gen",&vx_gen,"vx_gen/D");
	LAEC->Branch("vy_gen",&vy_gen,"vy_gen/D");
	LAEC->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	LAEC->Branch("evn_flux",&evn_flux,"evn_flux/I");
	LAEC->Branch("mom_flux",&mom_flux,"mom_flux/D");
	LAEC->Branch("theta_flux",&theta_flux,"theta_flux/D");
	LAEC->Branch("phi_flux",&phi_flux,"phi_flux/D");
	LAEC->Branch("r_flux",&r_flux,"r_flux/D");
	LAEC->Branch("x_flux",&x_flux,"x_flux/D");
	LAEC->Branch("y_flux",&y_flux,"y_flux/D");
	LAEC->Branch("z_flux",&z_flux,"z_flux/D");
	LAEC->Branch("vx_flux",&vx_flux,"vx_flux/D");
	LAEC->Branch("vy_flux",&vy_flux,"vy_flux/D");
	LAEC->Branch("vz_flux",&vz_flux,"vz_flux/D");
	LAEC->Branch("px_flux",&px_flux,"px_flux/D");
	LAEC->Branch("py_flux",&py_flux,"py_flux/D");
	LAEC->Branch("pz_flux",&pz_flux,"pz_flux/D");
	LAEC->Branch("E",&E,"E/D");
	LAEC->Branch("Edep",&Edep,"Edep/D");
	LAEC->Branch("R_EC",&R_EC,"R_EC/D");
	LAEC->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	LAEC->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	LAEC->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	LAEC->Branch("ID_flux",&ID_flux,"ID_flux/I");
	LAEC->Branch("PID_flux",&PID_flux,"PID_flux/I");
	LAEC->Branch("TID_flux",&TID_flux,"TID_flux/I");
	LAEC->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	LAEC->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	
	/*FAEC{{{*/
	TTree *FAEC = new TTree("FAEC","A new Tree for particles going in FAEC");
	//From Header Tree
	FAEC->Branch("rate",&rate,"rate/D");
	FAEC->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	FAEC->Branch("mom_gen",&mom_gen,"mom_gen/D");
	FAEC->Branch("theta_gen",&theta_gen,"theta_gen/D");
	FAEC->Branch("phi_gen",&phi_gen,"phi_gen/D");
	FAEC->Branch("px_gen",&px_gen,"px_gen/D");
	FAEC->Branch("py_gen",&py_gen,"py_gen/D");
	FAEC->Branch("pz_gen",&pz_gen,"pz_gen/D");
	FAEC->Branch("vx_gen",&vx_gen,"vx_gen/D");
	FAEC->Branch("vy_gen",&vy_gen,"vy_gen/D");
	FAEC->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	FAEC->Branch("evn_flux",&evn_flux,"evn_flux/I");
	FAEC->Branch("mom_flux",&mom_flux,"mom_flux/D");
	FAEC->Branch("theta_flux",&theta_flux,"theta_flux/D");
	FAEC->Branch("phi_flux",&phi_flux,"phi_flux/D");
	FAEC->Branch("r_flux",&r_flux,"r_flux/D");
	FAEC->Branch("x_flux",&x_flux,"x_flux/D");
	FAEC->Branch("y_flux",&y_flux,"y_flux/D");
	FAEC->Branch("z_flux",&z_flux,"z_flux/D");
	FAEC->Branch("vx_flux",&vx_flux,"vx_flux/D");
	FAEC->Branch("vy_flux",&vy_flux,"vy_flux/D");
	FAEC->Branch("vz_flux",&vz_flux,"vz_flux/D");
	FAEC->Branch("px_flux",&px_flux,"px_flux/D");
	FAEC->Branch("py_flux",&py_flux,"py_flux/D");
	FAEC->Branch("pz_flux",&pz_flux,"pz_flux/D");
	FAEC->Branch("E",&E,"E/D");
	FAEC->Branch("Edep",&Edep,"Edep/D");
	FAEC->Branch("R_EC",&R_EC,"R_EC/D");
	FAEC->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	FAEC->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	FAEC->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	FAEC->Branch("ID_flux",&ID_flux,"ID_flux/I");
	FAEC->Branch("PID_flux",&PID_flux,"PID_flux/I");
	FAEC->Branch("TID_flux",&TID_flux,"TID_flux/I");
	FAEC->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	FAEC->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	/*}}}*/
	/*}}}*/

	/*Read in each event{{{*/
	ofstream laoutfile(Form("SoLID_SIDIS_He3_LA_%s.dat",pid_temp.Data()));
	laoutfile <<Form("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")<<endl;	
	laoutfile <<Form("  Generated: Evt#, PID, Px, Py, Pz, Vx, Vy, Vz, Mom, Theta, Phi")<<endl;	
	laoutfile <<Form("  GEM#1: Plane, PID, TID, Px, Py, Pz, X, Y, Z, Mom, Theta, Phi, E, Edep")<<endl;	
	laoutfile <<"  ..."<<endl;
	laoutfile <<Form("  GEM#4: Plane, PID, TID, Px, Py, Pz, X, Y, Z, Mom, Theta, Phi, E, Edep")<<endl;	
	laoutfile <<Form("  LAEC:  EC-ID, PID, TID, Px, Py, Pz, X, Y, Z, Mom, Theta, Phi, E")<<endl;	
	laoutfile <<Form("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")<<endl;	

	ofstream faoutfile(Form("SoLID_SIDIS_He3_FA_%s.dat",pid_temp.Data()));
	faoutfile <<Form("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")<<endl;	
	faoutfile <<Form("  Gener: Evt#,  PID, Px, Py, Pz, Mom, Theta, Phi")<<endl;	
	faoutfile <<Form("  GEM#1: Plane, PID, TID, Px, Py, Pz, X, Y, Z, Mom, Theta, Phi, E, Edep")<<endl;	
	faoutfile <<"  ..."<<endl;
	faoutfile <<Form("  GEM#6: Plane, PID, TID, Px, Py, Pz, X, Y, Z, Mom, Theta, Phi, E, Edep")<<endl;	
	faoutfile <<Form("  FAEC:  EC-ID, PID, TID, Px, Py, Pz, X, Y, Z, Mom, Theta, Phi, E")<<endl;	
	faoutfile <<Form("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB")<<endl;	

	Int_t nselected = nevent;
	Int_t LA_evn = 0, FA_evn = 0;
	cerr<<"++++++++++++++++ "<<endl;
	for(Long64_t i=0;i<nselected;i++){
		cout<<i<<"\r";
		//  From Header, head_evn, usally only one particle
		header->GetEntry(i);
		if(input_filename.find("_EM_",0) != string::npos) 
			rate = (((1.5e-5)/(1.6e-19))/nevent); //Count to Hz for 15uA electron events;
		else
			rate = head_rate->at(0)/filenum;
		evn = head_evn->at(0);

        /*Fill Generated Info {{{*/		
		//  From Generated, gen_pid,gen_vx/_vy/_vz,gen_px/_py/_pz
		generated->GetEntry(i);

		const Int_t ng = gen_pid->size();//Normally there is only one particle in the gen
		Double_t gen_theta_array[ng], gen_phi_array[ng], gen_mom_array[ng];
		Double_t gen_px_array[ng],gen_py_array[ng],gen_pz_array[ng],gen_vx_array[ng],gen_vy_array[ng],gen_vz_array[ng];
		Int_t Is_Pick = -1;
		for(Int_t ig=0;ig<ng;ig++){
			if((Int_t)gen_pid->at(ig)==PID_Pick)
				Is_Pick = ig;
			gen_mom_array[ig] = sqrt( pow(gen_px->at(ig),2)+pow(gen_py->at(ig),2)+pow(gen_pz->at(ig),2) ); //GeV
			gen_theta_array[ig] = acos(gen_pz->at(ig)/gen_mom_array[ig])*DEG;//Degree
			gen_phi_array[ig] = atan2( gen_py->at(ig), gen_px->at(ig))*DEG;//Degree
			gen_mom_array[ig] *= MeV2GeV;
			gen_px_array[ig] = gen_px->at(ig)*MeV2GeV; //GeV
			gen_py_array[ig] = gen_py->at(ig)*MeV2GeV; //GeV
			gen_pz_array[ig] = gen_pz->at(ig)*MeV2GeV; //GeV
			gen_vx_array[ig] = gen_vx->at(ig)*MM2CM; //cm
			gen_vy_array[ig] = gen_vy->at(ig)*MM2CM; //cm
			gen_vz_array[ig] = gen_vz->at(ig)*MM2CM; //cm
			//cerr<<Form("---#%d@%d Px=%f Py=%f Pz=%f P=%f Theta =%f Phi=%f ",i, ig,gen_px->at(ig),gen_py->at(ig),gen_pz->at(ig),gen_mom_array[ig],gen_theta_array[ig],gen_phi_array[ig])<<endl;
		}
		mom_gen = gen_mom_array[Is_Pick];theta_gen = gen_theta_array[Is_Pick];phi_gen = gen_phi_array[Is_Pick];
		px_gen = gen_px_array[Is_Pick];py_gen = gen_py_array[Is_Pick];pz_gen = gen_pz_array[Is_Pick];
		vx_gen = gen_vx_array[Is_Pick];vy_gen = gen_vy_array[Is_Pick];vz_gen = gen_vz_array[Is_Pick];

		/*}}}End Fill Generated info*/

		/*Fill Flux Info{{{*/
		//   From Flux, flux_avg_x/_y/_z/_vx/_vy/_vz/_px/_py/_pz
		//           flux_hitn,flux_id/_pid/_mpid/_mvx/_mvy/_mvz,
		flux->GetEntry(i);
		Int_t Count_GEM_Hit = 0; //Another 2GEMs after LA and before FA, so at least 2 from LA and 1 from FA
		bool Is_GEM_Hit[GEM_Plane];
		for( unsigned int ip=0;ip<GEM_Plane;ip++){
			GEM_X[ip] = 0.0;   GEM_Y[ip] = 0.0;     GEM_Z[ip] = 0.0;
			GEM_Px[ip] = 0.0;  GEM_Py[ip] = 0.0;    GEM_Pz[ip] = 0.0;
			GEM_Mom[ip] = 0.0; GEM_Theta[ip] = 0.0; GEM_Phi[ip] = 0.0;
			GEM_PID_Pick[ip] = -10;  GEM_ID_Pick[ip] = -99999999; 
			GEM_TID_Pick[ip] = -10;
		   	GEM_E[ip] = 0.0, GEM_Edep[ip]=0.0;
			Is_GEM_Hit[ip] = kFALSE;
		}
		Double_t LAEC_X=0.0,  LAEC_Y=0.0,  LAEC_Z=0.0;
        Double_t LAEC_Px=0.0, LAEC_Py=0.0, LAEC_Pz=0.0;
		Double_t LAEC_Mom=0.0, LAEC_Theta=0.0, LAEC_Phi=0.0, LAEC_E = 0.0;
        Double_t FAEC_X=0.0,  FAEC_Y=0.0,  FAEC_Z=0.0;
		Double_t FAEC_Px=0.0, FAEC_Py=0.0, FAEC_Pz=0.0;
		Double_t FAEC_Mom=0.0, FAEC_Theta=0.0, FAEC_Phi=0.0, FAEC_E = 0.0;
        Int_t LAEC_PID = -10;  Int_t FAEC_PID = -10;
        Int_t LAEC_TID = -10;  Int_t FAEC_TID = -10;
		bool Is_LAEC_Hit = kFALSE, Is_FAEC_Hit = kFALSE;	

		for (unsigned int j=0;j<flux_hitn->size();j++){

			/*Obtains Infos of each flux{{{*/
			evn_flux = j;
			x_flux = flux_avg_x->at(j)*MM2CM; y_flux = flux_avg_y->at(j)*MM2CM; z_flux = flux_avg_z->at(j)*MM2CM;
			vx_flux = flux_vx->at(j)*MM2CM;   vy_flux = flux_vy->at(j)*MM2CM;	vz_flux = flux_vz->at(j)*MM2CM;
			px_flux = flux_px->at(j)*MeV2GeV; py_flux = flux_py->at(j)*MeV2GeV; pz_flux = flux_pz->at(j)*MeV2GeV;
			Edep = flux_totEdep->at(j)*MeV2GeV; E = flux_trackE->at(j)*MeV2GeV;

			r_flux=sqrt(pow(x_flux,2)+pow(y_flux,2));//cm
			mom_flux=sqrt(pow(px_flux,2)+pow(py_flux,2)+pow(pz_flux,2));//GeV
			theta_flux = atan(sqrt(pow(px_flux,2)+pow(py_flux,2))/(pz_flux))*DEG;//Degree
			phi_flux = fabs(atan(y_flux/x_flux))*DEG;
			if(y_flux>0 && x_flux<0) phi_flux += 90.0;
			else if(y_flux<0 && x_flux<0) phi_flux += 180.0;
			else if(y_flux<0 && x_flux>0) phi_flux += 270.0;

			ID_flux = (Int_t) (flux_id->at(j));
			PID_flux = (Int_t) (flux_pid->at(j));
			TID_flux = (Int_t) (flux_tid->at(j));
			if(PID_flux==-Electron) PID_flux = Electron;//Force e+/- to be counted at one time
            /*}}}*/

			//Select the charged particles
//			if(PID_flux!=Electron&&PID_flux!=Proton
//					&&PID_flux!=Pion&&PID_flux!=-Pion
//					&&PID_flux!=Kaon&&PID_flux!=-Kaon) continue;
			if(PID_flux!=PID_Pick) continue;
			//if(flux_pz->at(j)<1e-9)continue;//Cout out backward particles
			if(mom_flux<1e-12) continue; //Cut out Zero-E particles
          
			/*Find Hits on GEMs{{{*/
			for( unsigned int ip=0;ip<GEM_Plane;ip++){
				if(ID_flux!=GEM_ID[ip]) continue; //Only check particles in one GEM layer
				if(r_flux>GEM_R_Max[ip] || r_flux<GEM_R_Min[ip]) continue;//Only inside the GEM Disk

				//Not yet select the true hits but just record everything on the plane
                if(ip==0) G1->Fill();
				else if(ip==1) G2->Fill();
				else if(ip==2) G3->Fill();
				else if(ip==3) G4->Fill();
				else if(ip==4) G5->Fill();
				else if(ip==5) G6->Fill();

				//Now Pick the true hits
				if(mom_flux>=mom_gen*0.8 && mom_flux<mom_gen){//assume the particle can not lose more than 20% of its energy
					GEM_ID_Pick[ip] = ID_flux;
					GEM_PID_Pick[ip] = PID_flux;
					GEM_TID_Pick[ip] = TID_flux;
					GEM_X[ip] = x_flux;
					GEM_Y[ip] = y_flux;
					GEM_Z[ip] = z_flux;
					GEM_Px[ip] = px_flux;
					GEM_Py[ip] = py_flux;
					GEM_Pz[ip] = pz_flux;
					GEM_Mom[ip] = mom_flux; 
					GEM_Theta[ip] = theta_flux; 
					GEM_Phi[ip] = phi_flux;
					GEM_E[ip] = E;
					GEM_Edep[ip] = Edep;
					Is_GEM_Hit[ip] = kTRUE;
					Count_GEM_Hit++;
				}
			}
            /*}}}*/

			/*Find Hits on ECs{{{*/
			if(ID_flux==LAEC_VP_ID && r_flux >= LAEC_VP_Rmin && r_flux < LAEC_VP_Rmax ){
				LAEC->Fill();
				if(mom_flux>=mom_gen*0.8 && mom_flux<mom_gen){//assume the particle can not lose more than 20% of its energy
					Is_LAEC_Hit = kTRUE;
					LAEC_PID = PID_flux;
					LAEC_TID = TID_flux;
					LAEC_X = x_flux;
					LAEC_Y = y_flux;
					LAEC_Z = z_flux;
					LAEC_Px = px_flux;
					LAEC_Py = py_flux;
					LAEC_Pz = pz_flux;
					LAEC_Mom = mom_flux;
					LAEC_Theta = theta_flux;
					LAEC_Phi = phi_flux;
					LAEC_E = E;
				}
			}
			if(ID_flux==FAEC_VP_ID && r_flux >= FAEC_VP_Rmin && r_flux < FAEC_VP_Rmax ){
				FAEC->Fill();
				if(mom_flux>=mom_gen*0.8 && mom_flux<mom_gen){//assume the particle can not lose more than 20% of its energy
					Is_FAEC_Hit = kTRUE;
					FAEC_PID = PID_flux;
					FAEC_TID = TID_flux;
					FAEC_X = x_flux;
					FAEC_Y = y_flux;
					FAEC_Z = z_flux;
					FAEC_Px = px_flux;
					FAEC_Py = py_flux;
					FAEC_Pz = pz_flux;
					FAEC_Mom = mom_flux;
					FAEC_Theta = theta_flux;
					FAEC_Phi = phi_flux;
					FAEC_E = E;
					}
			}
			/*}}}*/
		}
		/*}}}*/

		if(!(Is_FAEC_Hit) && !(Is_LAEC_Hit)) continue;

	    /*Print_t Out{{{*/
		//if(Is_LAEC_Hit&&Count_GEM_Hit_LA>=3 && !(Is_FAEC_Hit)){
		if(Is_LAEC_Hit && Count_GEM_Hit >=3 && !(Is_GEM_Hit[4]) && !(Is_GEM_Hit[5]) ){
            LA_evn ++;
			laoutfile <<Form("%5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
					LA_evn, PID_Pick,
					px_gen, py_gen, pz_gen,
					mom_gen, theta_gen, phi_gen)
				<<endl;
			for(Int_t ip=0;ip<4;ip++){//four GEMCs in front of LAEC
				laoutfile <<Form("%2d %5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
						ip+1, GEM_PID_Pick[ip], GEM_TID_Pick[ip],
						GEM_Px[ip], GEM_Py[ip], GEM_Pz[ip],
						GEM_X[ip], GEM_Y[ip], GEM_Z[ip],
						GEM_Mom[ip],GEM_Theta[ip],GEM_Phi[ip],
						GEM_E[ip], GEM_Edep[ip])
					<<endl;
			}
			laoutfile <<Form("%2d %5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
			        310,  LAEC_PID,LAEC_TID,
					LAEC_Px, LAEC_Py, LAEC_Pz,
					LAEC_X,  LAEC_Y,  LAEC_Z,
					LAEC_Mom,LAEC_Theta,LAEC_Phi,LAEC_E)
					<<endl;
			}
		if(Is_FAEC_Hit&&Count_GEM_Hit>=5){
            FA_evn ++;
			faoutfile <<Form("%5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
					FA_evn, PID_Pick,
					px_gen, py_gen, pz_gen,
					mom_gen, theta_gen, phi_gen)
				<<endl;
			for(Int_t ip=0;ip<6;ip++){//four GEMCs in front of LAEC
				faoutfile <<Form("%2d %5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
						ip+1, GEM_PID_Pick[ip], GEM_TID_Pick[ip],
						GEM_Px[ip], GEM_Py[ip], GEM_Pz[ip],
						GEM_X[ip], GEM_Y[ip], GEM_Z[ip],
						GEM_Mom[ip],GEM_Theta[ip],GEM_Phi[ip],
						GEM_E[ip], GEM_Edep[ip])
					<<endl;
			}
			faoutfile <<Form("%2d %5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
			        320, FAEC_PID, FAEC_TID,
					FAEC_Px, FAEC_Py, FAEC_Pz,
					FAEC_X, FAEC_Y, FAEC_Z,
					FAEC_Mom,FAEC_Theta,FAEC_Phi,FAEC_E)
					<<endl;
			}
	}
    /*}}}End PrInt_t Out*/

	/*End Read in each event}}}*/

	out_file->cd(); 
	G1->Write(); G2->Write(); G3->Write(); G4->Write(); G5->Write(); G6->Write(); 
	LAEC->Write(); FAEC->Write();
    out_file->Close();
//	file->Close();
    header->Delete();
	generated->Delete();
	flux->Delete();
	laoutfile.close();
	faoutfile.close();

	return 0;
}
