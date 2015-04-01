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
  //#include <TMatrix.h>
  /*}}}*/
using namespace std;

void generate_weight(){
	TString inputname; cerr<<"--- Rootfile please = "; cin >>inputname;
	double ele_energy; cerr<<"--- Electron Energy = (8.8 or 11)"; cin >>ele_energy;

	TFile *file=new TFile(inputname.Data(),"update");  // save Q2<10 pt<1
	TTree *T=(TTree*)file->Get("T");
	int nsim=0;
	double dxs_hp=0,dxs_hm=0;
	T->SetBranchAddress("nsim",&nsim);
	T->SetBranchAddress("dxs_hp",&dxs_hp);
	T->SetBranchAddress("dxs_hm",&dxs_hm);

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(N_entries-1);          //get nsim for this rootfile
	
	double N_simulate=(double)(nsim);
	cout<<"N_simulate: "<<N_simulate<<endl;

	double electron_phase_space=(cos(7/180.*3.1415926) - cos(30/180.*3.1415926))*2*3.14159265*(ele_energy-0.5);   // theta: 7~30 degree,  2pi phi coverage, 0.5~11 GeV Momentum coverage 	
	double hadron_phase_space=(cos(7/180.*3.1415926) - cos(30/180.*3.1415926))*2*3.14159265*(6-0.5);  //theta, 7~30 degree,  2pi phi coverage, 0.5~6 GeV Momentum coverage
	double Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
	cout<<"Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;
	double weight_hp=0;
	double weight_hm=0;

	TBranch *branch_weight_hp=T->Branch("weight_hp",&weight_hp,"weight_hp/D");
	TBranch *branch_weight_hm=T->Branch("weight_hm",&weight_hm,"weight_hm/D");

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		//warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
		weight_hp=dxs_hp*Phase_space/N_simulate;   
		weight_hm=dxs_hm*Phase_space/N_simulate;
		branch_weight_hp->Fill();
		branch_weight_hm->Fill();
	}

	T->Write("",TObject::kOverwrite);
}



