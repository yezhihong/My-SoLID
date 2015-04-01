#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TEllipse.h>
#include <TBox.h>
#include <TLorentzVector.h>
#include <TVector3.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>
#include <TColor.h>

using namespace std;

void acceptance(string input_filename)
{
	// gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	// gStyle->SetOptStat(1111111);

	const double DEG=180./3.1415926;
	char the_filename[200];
	sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());

	TFile *file=new TFile(input_filename.c_str());
	if (file->IsZombie()) {
		cout << "Error opening file" << input_filename << endl;
		exit(-1);
	}
	else cout << "open file " << input_filename << endl;

	TTree *header = (TTree*) file->Get("header");
	// vector <int> *evn=0,*evn_type=0;
	vector <double> *evn=0,*evn_type=0;
	vector <double> *beamPol=0;
	vector <double> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0;
	header->SetBranchAddress("evn",&evn);
	header->SetBranchAddress("evn_type",&evn_type);
	header->SetBranchAddress("beamPol",&beamPol);
	header->SetBranchAddress("var1",&var1);
	header->SetBranchAddress("var2",&var2);
	header->SetBranchAddress("var3",&var3);
	header->SetBranchAddress("var4",&var4);
	header->SetBranchAddress("var5",&var5);
	header->SetBranchAddress("var6",&var6);
	header->SetBranchAddress("var7",&var7);
	header->SetBranchAddress("var8",&var8);

	TTree *generated = (TTree*) file->Get("generated");
	vector <int> *gen_pid=0;
	vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	generated->SetBranchAddress("pid",&gen_pid);
	generated->SetBranchAddress("px",&gen_px);
	generated->SetBranchAddress("py",&gen_py);
	generated->SetBranchAddress("pz",&gen_pz);
	generated->SetBranchAddress("vx",&gen_vx);
	generated->SetBranchAddress("vy",&gen_vy);
	generated->SetBranchAddress("vz",&gen_vz);

	TTree *flux = (TTree*) file->Get("flux");
	vector<double> *flux_id=0,*flux_hitn=0;
	vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
	flux->SetBranchAddress("hitn",&flux_hitn);
	flux->SetBranchAddress("id",&flux_id);
	flux->SetBranchAddress("pid",&flux_pid);
	flux->SetBranchAddress("mpid",&flux_mpid);
	flux->SetBranchAddress("tid",&flux_tid);
	flux->SetBranchAddress("mtid",&flux_mtid);
	flux->SetBranchAddress("otid",&flux_otid);
	flux->SetBranchAddress("trackE",&flux_trackE);
	flux->SetBranchAddress("totEdep",&flux_totEdep);
	flux->SetBranchAddress("avg_x",&flux_avg_x);
	flux->SetBranchAddress("avg_y",&flux_avg_y);
	flux->SetBranchAddress("avg_z",&flux_avg_z);
	flux->SetBranchAddress("avg_lx",&flux_avg_lx);
	flux->SetBranchAddress("avg_ly",&flux_avg_ly);
	flux->SetBranchAddress("avg_lz",&flux_avg_lz);
	flux->SetBranchAddress("px",&flux_px);
	flux->SetBranchAddress("py",&flux_py);
	flux->SetBranchAddress("pz",&flux_pz);
	flux->SetBranchAddress("vx",&flux_vx);
	flux->SetBranchAddress("vy",&flux_vy);
	flux->SetBranchAddress("vz",&flux_vz);
	flux->SetBranchAddress("mvx",&flux_mvx);
	flux->SetBranchAddress("mvy",&flux_mvy);
	flux->SetBranchAddress("mvz",&flux_mvz);
	flux->SetBranchAddress("avg_t",&flux_avg_t);

	int nevent = (int)generated->GetEntries();
	int nselected = 0;
	cout << "nevent " << nevent << endl;

	double  weight=1;

	int counter=0;

	for (Long64_t i=0;i<nevent;i++) {
		cerr<<" I = "<<i<<"\r";
		header->GetEntry(i);
		generated->GetEntry(i);

		int pid_gen;
		double p_gen,theta_gen,phi_gen,px_gen,py_gen,pz_gen,vx_gen,vy_gen,vz_gen;
		double t=0,Q2=0;
		double alpha_x,alpha_y;
		double theta_gen_ion,phi_gen_ion;  

		for (int j=0;j<gen_pid->size();j++) {
			pid_gen=gen_pid->at(j);
			px_gen=gen_px->at(j)/1e3;    	//in MeV, convert to GeV
			py_gen=gen_py->at(j)/1e3;		//in MeV, convert to GeV
			pz_gen=gen_pz->at(j)/1e3;      	//in MeV, convert to GeV
			vx_gen=gen_vx->at(j)/1e1;    	//in mm, convert to cm
			vy_gen=gen_vy->at(j)/1e1;		//in mm, convert to cm
			vz_gen=gen_vz->at(j)/1e1;     	//in mm, convert to cm
		}  

		flux->GetEntry(i);    
		double px_flux,py_flux,pz_flux,theta_flux,p_flux,phi_flux;    
		for (int j=0;j<flux_id->size();j++) {
			if (flux_id->at(j) > 30000) continue;
			int beamline = int(flux_id->at(j))%100000/10000-1;
			int side     = int(flux_id->at(j))%10000/1000-1;
			int magnet   = int(flux_id->at(j))%1000/100-1;    
			int number   = int(flux_id->at(j))%100/10-1;
			int window   = int(flux_id->at(j))%10/1-1;      
		}

		for (int j=0;j<flux_id->size();j++) {
			if (flux_id->at(j) < 30000) {
				int beamline = int(flux_id->at(j))%100000/10000-1;
				int side     = int(flux_id->at(j))%10000/1000-1;
				int magnet   = int(flux_id->at(j))%1000/100-1;    
				int number   = int(flux_id->at(j))%100/10-1;
				int window   = int(flux_id->at(j))%10/1-1;

			}
		}
	}
}

