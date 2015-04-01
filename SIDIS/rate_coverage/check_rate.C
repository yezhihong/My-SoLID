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

void check_rate(TString particle){  
	// particle="pip" or "pim"
	TChain *T=new TChain("T");
	T->Add("./collider_3he_pip_11_0_1_0.root");   //Q2<10, pt<1
	T->Add("./collider_3he_pip_11_0_2_0.root");   //Q2<10, pt>1
	double Q2,x,pt,W,Wp,z;   //kinematics
	double mom_ele,mom_had,theta_ele,theta_had;  // acceptance calculation input
	double weight_hp,weight_hm;

	T->SetBranchAddress("Q2",&Q2);
	T->SetBranchAddress("x",&x);
	T->SetBranchAddress("pt",&pt);
	T->SetBranchAddress("W",&W);
	T->SetBranchAddress("Wp",&Wp);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("mom_ele",&mom_ele);
	T->SetBranchAddress("mom_had",&mom_had);
	T->SetBranchAddress("theta_ele",&theta_ele);
	T->SetBranchAddress("theta_had",&theta_had);
	T->SetBranchAddress("weight_hp",&weight_hp);
	T->SetBranchAddress("weight_hm",&weight_hm);

	Long64_t N_entries=T->GetEntries();
	cout<<"total generated events number: "<<N_entries<<endl;

	TFile *file_negative=new TFile("./solid_acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TH2F *accep_negative_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_negative_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance

	double rate_forward=0;
	double rate_large=0;

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		if(mom_ele>1.0 && mom_ele<10.0 && W>=2.3 && Wp>=1.6 && z>0.3 && z<0.7){ //zhihong's cut
			int ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
			int ele_p_bin=int(mom_ele/0.05)+1;      //0.05 GeV per bin
			double ele_forward_acceptance=accep_negative_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			double ele_large_acceptance=accep_negative_large->GetBinContent(ele_theta_bin,ele_p_bin);
			//mandatary 8 degree cut
			if(mom_ele<1.0||theta_ele/3.14159265*180.0<8.0 ||theta_ele/3.14159265*180.>14.8){
				ele_forward_acceptance=0;						                          
			}
			if(mom_ele<3.5||theta_ele/3.14159265*180.0<15.5||theta_ele/3.14159265*180.>24.){
				ele_large_acceptance=0;
			}

			int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
			int hadron_p_bin=int(mom_had/0.05)+1;     //0.05 GeV per bin
			double hadron_acceptance=0;
			hadron_acceptance=accep_negative_forward->GetBinContent(hadron_theta_bin,hadron_p_bin);
			double event_weight=0;   //depend on hadron charge

				event_weight=weight_hp;
//				event_weight=weight_hm;

			//mandatary 8 degree cut
			if(theta_had/3.14159265*180.0<8.0){
				hadron_acceptance=0;
			}

			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

			rate_forward+=event_weight*forward_acceptance;
			rate_large+=event_weight*large_acceptance;
		}// if cut ends here
	}//event loop ends here

	cout<<"rate forward: "<<rate_forward<<endl;
	cout<<"rate large: "<<rate_large<<endl;
}


