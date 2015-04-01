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


void my(TString particle,int hadron_flag){
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
	//read in acceptance histograms
	//electron acceptance	
	
	//read in acceptance histograms
	TFile *file_positive=new TFile("./acceptance/acceptance_solid_CLEO_SIDIS_3he_positive_output.root");
	TH2F *accep_positive_forward=(TH2F*)file_positive->Get("acceptance_forwardangle");
	//no need to care about positive large angle accep, since we don't care positron acceptance right now
	
	TFile *file_negative=new TFile("./acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TH2F *accep_negative_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_negative_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance


	double rate_forward=0;
	double rate_large=0;

	double rate_forward_p_g_4[5]={0};
	double rate_large_p_g_4[5]={0};
	double rate_forward_p_l_4=0;
	double rate_large_p_l_4=0;
	int ele_theta_bin= 0;
	int ele_p_bin=0;
	double ele_forward_acceptance= 0.0;
	double ele_large_acceptance= 0.0;

	int N_break = 0;

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		ele_theta_bin= 0;
		ele_p_bin=0;
		ele_forward_acceptance= 0.0;
		ele_large_acceptance= 0.0;

		if(theta_ele*DEG>=7.0&&theta_ele*DEG<=30.0&&mom_ele>=1.0&&mom_ele<=11.0&&W>=2.3&&Wp>=1.6){//any additional cuts should be added in here
		//if(W>=2.3&&Wp>=1.6){//any additional cuts should be added in here
			ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
			ele_p_bin=int(mom_ele/0.1)+1;      //0.05 GeV per bin for mom
			ele_forward_acceptance=accep_negative_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			ele_large_acceptance=accep_negative_large->GetBinContent(ele_theta_bin,ele_p_bin);

			int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
			int hadron_p_bin=int(mom_had/0.1)+1;     //0.05 GeV per bin for mom
			double event_weight=0;   //depend on hadron charge
			
			double hadron_acceptance= 0.0;

			if(hadron_flag==1){  //doing hp
				event_weight=weight_hp;
				hadron_acceptance=accep_positive_forward->GetBinContent(hadron_theta_bin,hadron_p_bin);
			}else if(hadron_flag==2){  //doing hm
				event_weight=weight_hm;
				hadron_acceptance=accep_negative_forward->GetBinContent(hadron_theta_bin,hadron_p_bin);
			}

			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

			if(mom_ele<1.0)//GeV
				forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(mom_ele<3.5)//GeV
				large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV

			/*Rates{{{*/	
			if( mom_ele>=1.0 && mom_had>=0.0){
				rate_forward+=event_weight*forward_acceptance;
			}
			if( mom_ele>=3.5 && mom_had>=0.0){
				rate_large+=event_weight*large_acceptance;
			}

			if(mom_had>=3.5&&mom_had<4.0){
				if(mom_ele>=1.0 && theta_ele*DEG<=18.0)
					rate_forward_p_g_4[4]+=event_weight*forward_acceptance;
				if(mom_ele>=3.5 && theta_ele*DEG>=15.0)
					rate_large_p_g_4[4]+=event_weight*large_acceptance;
			}else if(mom_had>=4.0&&mom_had<4.5){
				if(mom_ele>=1.0 && theta_ele*DEG<=18.0)
					rate_forward_p_g_4[0]+=event_weight*forward_acceptance;
				if(mom_ele>=3.5 && theta_ele*DEG>=15.0)
					rate_large_p_g_4[0]+=event_weight*large_acceptance;
			}else if(mom_had>=4.5&&mom_had<5){
				if(mom_ele>=1.0 && theta_ele*DEG<=18.0)
					rate_forward_p_g_4[1]+=event_weight*forward_acceptance;
				if(mom_ele>=3.5 && theta_ele*DEG>=15.0)
					rate_large_p_g_4[1]+=event_weight*large_acceptance;
			}else if(mom_had>=5&&mom_had<5.5){
				if(mom_ele>=1.0 && theta_ele*DEG<=18.0)
					rate_forward_p_g_4[2]+=event_weight*forward_acceptance;
				if(mom_ele>=3.5 && theta_ele*DEG>=15.0)
					rate_large_p_g_4[2]+=event_weight*large_acceptance;
			}else if(mom_had>=5.5&&mom_had<6){
				if(mom_ele>=1.0 && theta_ele*DEG<=18.0)
					rate_forward_p_g_4[3]+=event_weight*forward_acceptance;
				if(mom_ele>=3.5 && theta_ele*DEG>=15.0)
					rate_large_p_g_4[3]+=event_weight*large_acceptance;
			}else{
				if(mom_ele>=1.0 && theta_ele*DEG<=18.0)
					rate_forward_p_l_4+=event_weight*forward_acceptance;
				if(mom_ele>=3.5 && theta_ele*DEG>=15.0)
					rate_large_p_l_4+=event_weight*large_acceptance;
			}
			/*}}}*/

			cerr<<Form("--- Electron: P=%f,B_P=%d, T=%f, B_T=%d, acc_for=%f, acc_lar=%f ", mom_ele, ele_p_bin, theta_ele, ele_theta_bin, ele_forward_acceptance, ele_large_acceptance)<<endl;
			cerr<<Form("---   Hadron: P=%f,B_P=%d, T=%f, B_T=%d, acc_for=%f,R_F=%f, R_L=%f", mom_had, hadron_p_bin, theta_had, hadron_theta_bin,hadron_acceptance,rate_forward,rate_large)<<endl;


			N_break ++;
			if(N_break>10)
				break;

		}

	}// events loop ends here


	const double Lumi = 3.0e36; // cm-2*s-1
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;

	/*Print&Save{{{*/
	cout<<"______forward and large angle integral rate____________________"<<endl;
	cout<<"rate_forward: "<<rate_forward*Lumi*KHz*nBcm2<<endl;
	cout<<"rate_large: "<<rate_large*Lumi*KHz*nBcm2<<endl;

	cout<<"______________ p~[3.5,4.0) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[4]*Lumi*KHz*nBcm2<<endl;
	cout<<"large:   "<<rate_large_p_g_4[4]*Lumi*KHz*nBcm2<<endl;

	cout<<"______________ p~[4,4.5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[0]*Lumi*KHz*nBcm2<<endl;
	cout<<"large:   "<<rate_large_p_g_4[0]*Lumi*KHz*nBcm2<<endl;

	cout<<"______________ p~[4.5,5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[1]*Lumi*KHz*nBcm2<<endl;
	cout<<"large:   "<<rate_large_p_g_4[1]*Lumi*KHz*nBcm2<<endl;

	cout<<"______________ p~[5,5.5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[2]*Lumi*KHz*nBcm2<<endl;
	cout<<"large:   "<<rate_large_p_g_4[2]*Lumi*KHz*nBcm2<<endl;

	cout<<"______________ p~[5.5, 6) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[3]*Lumi*KHz*nBcm2<<endl;
	cout<<"large:   "<<rate_large_p_g_4[3]*Lumi*KHz*nBcm2<<endl;

	cout<<"______________ 0.8<p<3.5 GeV rate __________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_l_4*Lumi*KHz*nBcm2<<endl;
	cout<<"large:   "<<rate_large_p_l_4*Lumi*KHz*nBcm2<<endl;
	cout<<"_______________________________________________________________"<<endl;
	/*}}}*/
}

