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
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>


struct GEM{
	int pid;
	int plane;
	double px; //GeV
	double py; //GeV
	double pz; //GeV
	double vx; //cm
	double vy; //cm
	double vz; //cm
	double x;  //cm
	double y;  //cm
	double z;  //cm
	double mom;//GeV
	double theta;//degree
	double phi; //degree
};
struct EC{ 
	int id;
	int pid;
	double px;
	double py;
	double pz;
	double vx;
	double vy;
	double vz;
	double x;
	double y;
	double z;
	double mom;
	double theta;
	double phi;
};
struct GENERATED{
	double px;
	double py;
	double pz;
	double vx;
	double vy;
	double vz;
	double mom;
	double theta;
	double phi;
	int pid;
	int evn;
};

void genroot(TString input_filename){

	int N_plane = 0;
	if(input_filename.Contains("FA")) N_plane = 6;
	else if(input_filename.Contains("LA")) N_plane = 4;
	else {
		cerr<<"*** ERROR, I don't know the file name"<<endl;
		exit -1;
	}
    const int N = N_plane;
	//ifstream infile_fa("SoLID_SIDIS_He3_FA.dat");
	GEM gem[N];//4 GEM planes for Large-Angle and 6 GEM planes for Forward-Angle
	GENERATED genr;
	EC ec;

	TString output_filename = input_filename;
	output_filename.ReplaceAll("dat","root");
    TFile *f1 = new TFile(output_filename.Data(),"recreate");
    TTree* T = new TTree("T","A new tree");
	T->Branch("evn", &genr.evn,"evn/I");
	T->Branch("pid", &genr.pid,"pid/I");
	T->Branch("px", &genr.px,"px/D");
	T->Branch("py", &genr.py,"py/D");
	T->Branch("pz", &genr.pz,"pz/D");
	T->Branch("vx", &genr.vx,"vx/D");
	T->Branch("vy", &genr.vy,"vy/D");
	T->Branch("vz", &genr.vz,"vz/D");
	T->Branch("mom", &genr.mom,"mom/D");
	T->Branch("theta", &genr.theta,"theta/D");
	T->Branch("phi", &genr.phi,"phi/D");

	for(int i=0;i<N;i++){
		T->Branch(Form("gem%d_x", i), &gem[i].x,Form("gem%d_x/D",i));
		T->Branch(Form("gem%d_y", i), &gem[i].y,Form("gem%d_y/D",i));
		T->Branch(Form("gem%d_z", i), &gem[i].z,Form("gem%d_z/D",i));
		T->Branch(Form("gem%d_vx", i), &gem[i].vx,Form("gem%d_vx/D",i));
		T->Branch(Form("gem%d_vy", i), &gem[i].vy,Form("gem%d_vy/D",i));
		T->Branch(Form("gem%d_vz", i), &gem[i].vz,Form("gem%d_vz/D",i));
		T->Branch(Form("gem%d_px", i), &gem[i].px,Form("gem%d_px/D",i));
		T->Branch(Form("gem%d_py", i), &gem[i].py,Form("gem%d_py/D",i));
		T->Branch(Form("gem%d_pz", i), &gem[i].pz,Form("gem%d_pz/D",i));
		T->Branch(Form("gem%d_mom", i), &gem[i].mom,Form("gem%d_mom/D",i));
		T->Branch(Form("gem%d_phi", i), &gem[i].phi,Form("gem%d_phi/D",i));
		T->Branch(Form("gem%d_theta", i), &gem[i].theta,Form("gem%d_theta/D",i));
		T->Branch(Form("gem%d_plane", i), &gem[i].plane,Form("gem%d_plane/D",i));
		T->Branch(Form("gem%d_pid", i), &gem[i].pid,Form("gem%d_pid/D",i));
	}
		
	T->Branch("ec_x", &ec.x,"ec_x/D");
	T->Branch("ec_y", &ec.y,"ec_y/D");
	T->Branch("ec_z", &ec.z,"ec_z/D");
	T->Branch("ec_vx", &ec.vx,"ec_vx/D");
	T->Branch("ec_vy", &ec.vy,"ec_vy/D");
	T->Branch("ec_vz", &ec.vz,"ec_vz/D");
	T->Branch("ec_px", &ec.px,"ec_px/D");
	T->Branch("ec_py", &ec.py,"ec_py/D");
	T->Branch("ec_pz", &ec.pz,"ec_pz/D");
	T->Branch("ec_mom", &ec.mom,"ec_mom/D");
	T->Branch("ec_phi", &ec.phi,"ec_phi/D");
	T->Branch("ec_theta", &ec.theta,"ec_theta/D");
	T->Branch("ec_id", &ec.id,"ec_id/I");
	T->Branch("ec_pid", &ec.pid,"ec_pid/I");


	int count=0;
	ifstream infile_la(input_filename.Data());
	while(!(infile_la.eof())){
		infile_la>>genr.evn>>genr.pid>>genr.px>>genr.py>> genr.pz>>genr.vx>>genr.vy>>genr.vz; 
		genr.mom = sqrt( pow(genr.px,2)+pow(genr.py,2)+pow(genr.pz,2) );
		genr.theta = acos(genr.pz/genr.mom) * 180./3.14;
		genr.phi   = atan2(genr.py,genr.px) * 180./3.14;

		for(int j=0;j<N;j++){ 
			infile_la>>gem[j].plane>>gem[j].pid>>gem[j].px>>gem[j].py>> gem[j].pz>>gem[j].vx>>gem[j].vy>>gem[j].vz>>gem[j].x>>gem[j].y>>gem[j].z;
			gem[j].mom = sqrt( pow(gem[j].px,2)+pow(gem[j].py,2)+pow(gem[j].pz,2) );
			gem[j].theta = acos(gem[j].pz/gem[j].mom) * 180./3.14;
			gem[j].phi   = atan2(gem[j].py,gem[j].px) * 180./3.14;
			//cerr<<" GEM Plane = "<<gem[j].plane<<endl;
		}
		infile_la>>ec.id>>ec.pid>>ec.px>>ec.py>>ec.pz>>ec.vx>>ec.vy>>ec.vz>>ec.x>>ec.y>>ec.z; 
		ec.mom = sqrt( pow(ec.px,2)+pow(ec.py,2)+pow(ec.pz,2) );
		ec.theta = acos(ec.pz/ec.mom) * 180./3.14;
		ec.phi   = atan2(ec.py,ec.px) * 180./3.14;
	    T->Fill();

		count++;
		cerr<<count<<"\r";
	}
   T->Write(); f1->Close();
   infile_la.close();
}
