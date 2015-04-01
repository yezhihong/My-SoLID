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

/*
//_________________acceptance file__________
acceptance_solid_CLEO_SIDIS_3he_negative_output.root
*/

void my_rate(TString particle,int hadron_flag){
	//particle: pion    kaon
	//hadron_flag   1 for hp      2 for hm
	TString charge = "A";
	if(hadron_flag==1)
		charge = "p";
	else if(hadron_flag==2)
		charge= "m";
	else
		cerr<<"********* Error"<<endl;
	gStyle->SetOptStat(1);
	TChain *T=new TChain("T");
	if(particle=="pion"){
		if(hadron_flag==1){
			T->Add("./collider_3he_pip_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_pip_11_0_2_0.root");   // save Q2<10, pt<1
		}
		else if(hadron_flag==2){
			T->Add("./collider_3he_pim_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_pim_11_0_2_0.root");   // save Q2<10, pt<1
		}
	}else if(particle=="kaon"){
		if(hadron_flag==1){
			T->Add("./collider_3he_kp_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_kp_11_0_2_0.root");   // save Q2<10, pt>1
		}
		if(hadron_flag==2){
			T->Add("./collider_3he_km_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_km_11_0_2_0.root");   // save Q2<10, pt>1
		}
	}

	//useful variables: Q2, W, Wp, x,y,z nu, s, pt, phi_h, phi_s weight_hp, weight_hm   
	//weight is xs in nbarn  weight_hp= dxs_hp*cos_theta_coverage*phi_coverage*energy_coverage/N_simulate
	const double DEG=180./3.1415926;
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
	T->SetBranchAddress("weight_hp",&weight_hp);   //unit in nbar
	T->SetBranchAddress("weight_hm",&weight_hm);

	Long64_t N_entries=T->GetEntries();
	cout<<"total generated events number: "<<N_entries<<endl;

	//deal with acceptance files 
    //CLEO	
	TFile *file_negative=new TFile("./solid_acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
	//CLEO
	//TFile *file_hadron=new TFile("./solid_acceptance/acceptance_solid_CLEO_SIDIS_3he_pionm_output.root");
	TH2F *accep_had=(TH2F*)accep_ele_forward->Clone("acceptance_hadron");//pion is normal acceptance histogram, the same as electron acceptance
	
	double rate_forward=0;
	double rate_large=0;
	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		if(mom_ele>=1.0&&mom_ele<=10.0&&W>=2.3&&Wp>=1.6&&z>0.3&&z<0.7){//any additional cuts should be added in here
			int ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
			int ele_p_bin=int(mom_ele/0.05)+1;      //0.05 GeV per bin for mom
			double ele_forward_acceptance=accep_ele_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			double ele_large_acceptance=accep_ele_large->GetBinContent(ele_theta_bin,ele_p_bin);

			if(mom_ele<1.0||theta_ele*DEG>14.5||theta_ele*DEG<8.0)//GeV, CLEO
				ele_forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(mom_ele<3.5||theta_ele*DEG<15.5||theta_ele*DEG>24)//GeV,CLEO
				ele_large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV

			int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
			int hadron_p_bin=int(mom_had/0.05)+1;     //0.05 GeV per bin for mom
			double hadron_acceptance=accep_had->GetBinContent(hadron_theta_bin,hadron_p_bin);
			if(mom_had<1.0||theta_had*DEG>14.5||theta_had*DEG<8.0)//GeV, CLEO
				hadron_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV

			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

			double event_weight=0;   //depend on hadron charge
				event_weight=weight_hp;
//				event_weight=weight_hm;


			rate_forward+=event_weight*forward_acceptance;
			rate_large+=event_weight*large_acceptance;

		}

	}// events loop ends here


	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;
	cout<<"______forward and large angle integral rate____________________"<<endl;
	cout<<"rate_forward: "<<rate_forward*Lumi*KHz*nBcm2<<endl;
	cout<<"rate_large: "<<rate_large*Lumi*KHz*nBcm2<<endl;
	cout<<"_______________________________________________________________"<<endl;
}

