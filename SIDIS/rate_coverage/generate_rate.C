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


void generate_rato(TString particle,int hadron_flag){
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

/*
		if(particle=="pion"){
		if(hadron_flag==1){
			T->Add("./sidis_3he_pip_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./sidis_3he_pip_11_0_2_0.root");   // save Q2<10, pt<1
		}
		else if(hadron_flag==2){
			T->Add("./sidis_3he_pim_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./sidis_3he_pim_11_0_2_0.root");   // save Q2<10, pt<1
		}
	}else if(particle=="kaon"){
		if(hadron_flag==1){
			T->Add("./sidis_3he_kp_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./sidis_3he_kp_11_0_2_0.root");   // save Q2<10, pt>1
		}
		if(hadron_flag==2){
			T->Add("./sidis_3he_km_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./sidis_3he_km_11_0_2_0.root");   // save Q2<10, pt>1
		}
	}
*/
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
		
    //CLEO	
	TFile *file_negative=new TFile("./solid_acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
	//CLEO
	//TFile *file_hadron=new TFile("./solid_acceptance/acceptance_solid_CLEO_SIDIS_3he_pionm_output.root");
	TH2F *accep_had=(TH2F*)accep_ele_forward->Clone("acceptance_hadron");//pion is normal acceptance histogram, the same as electron acceptance
	

    //CDF		
	//TFile *file_negative=new TFile("Xin_CDF_Acceptance_output.root");
//	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
//	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
	//CDF
	//TH2F *accep_had=(TH2F*)file_negative->Get("acceptance");           //negative particle large angle acceptance

	/*Define Histograms{{{*/
	int Nbinx=100;
	int Nbiny=100;
	int Nbin_theta=100;
	int Nbin_p=100;
	//TFile *savefile=new TFile("generate_plot_for_yi_output.root","recreate");
	// __________z vs pt 
	//forward angle
	TH2F *hf_z_pt=new TH2F("hf_z_pt","z VS P_{t} ",Nbinx,0.25,0.75,Nbiny,0,1.8);
	hf_z_pt->GetXaxis()->SetTitle("z");
	hf_z_pt->GetYaxis()->SetTitle("P_{t}(GeV)");
	//large angle
	TH2F *hl_z_pt=new TH2F("hl_z_pt","z VS P_{t} ",Nbinx,0.25,0.75,Nbiny,0,1.8);
	hl_z_pt->GetXaxis()->SetTitle("z");
	hl_z_pt->GetYaxis()->SetTitle("P_{t}(GeV)");
	//total
	TH2F *h_z_pt=new TH2F("h_z_pt","z VS P_{t} ",Nbinx,0.25,0.75,Nbiny,0,1.8);
	h_z_pt->GetXaxis()->SetTitle("z");
	h_z_pt->GetYaxis()->SetTitle("P_{t}(GeV)");

	//_________x vs pt 
	//forward angle
	TH2F *hf_x_pt=new TH2F("hf_x_pt","x VS P_{t} ",Nbinx,0,0.7,Nbiny,0,1.8);
	hf_x_pt->GetXaxis()->SetTitle("x_{bj}");
	hf_x_pt->GetYaxis()->SetTitle("P_{T} (GeV)");
	//large angle
	TH2F *hl_x_pt=new TH2F("hl_x_pt","x VS P_{t} ",Nbinx,0,0.7,Nbiny,0,1.8);
	hl_x_pt->GetXaxis()->SetTitle("x_{bj}");
	hl_x_pt->GetYaxis()->SetTitle("P_{T} (GeV)");
	//total
	TH2F *h_x_pt=new TH2F("h_x_pt","x VS P_{t} ",Nbinx,0,0.7,Nbiny,0,1.8);
	h_x_pt->GetXaxis()->SetTitle("x_{bj}");
	h_x_pt->GetYaxis()->SetTitle("P_{T} (GeV)");

	//_______Q2 vs pt  
	//forward angle
	TH2F *hf_Q2_pt=new TH2F("hf_Q2_pt","Q^{2} VS P_{t} ",Nbiny,0,10,Nbinx,0,1.8);
	hf_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hf_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");
	//large angle
	TH2F *hl_Q2_pt=new TH2F("hl_Q2_pt","Q^{2} VS P_{t} ",Nbinx,0,10,Nbiny,0,1.8);
	hl_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hl_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");
	//total
	TH2F *h_Q2_pt=new TH2F("h_Q2_pt","Q^{2} VS P_{t} ",Nbinx,0,10,Nbiny,0,1.8);
	h_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	h_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");

	//______pt vs p
	//forward angle
	TH2F *hf_pt_p=new TH2F("hf_pt_p","P_{t} VS P ",Nbinx,0,7,Nbiny,0,1.8);
	hf_pt_p->GetXaxis()->SetTitle("P (GeV)");
	hf_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");
	//large angle
	TH2F *hl_pt_p=new TH2F("hl_pt_p","P_{t} VS P ",Nbinx,0,7,Nbiny,0,1.8);
	hl_pt_p->GetXaxis()->SetTitle("P (GeV)");
	hl_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");
	//total
	TH2F *h_pt_p=new TH2F("h_pt_p","P_{t} VS P ",Nbinx,0,7,Nbiny,0,1.8);
	h_pt_p->GetXaxis()->SetTitle("P (GeV)");
	h_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");

	//______x vs p
	//forward angle
	TH2F *hf_x_p=new TH2F("hf_x_p","x VS P ",Nbinx,0,7,Nbiny,0,0.7);
	hf_x_p->GetYaxis()->SetTitle("x_{bj}");
	hf_x_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_x_p=new TH2F("hl_x_p","x VS P ",Nbinx,0,7,Nbiny,0,0.7);
	hl_x_p->GetYaxis()->SetTitle("x_{bj}");
	hl_x_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_x_p=new TH2F("h_x_p","x VS P ",Nbinx,0,7,Nbiny,0,0.7);
	h_x_p->GetYaxis()->SetTitle("x_{bj}");
	h_x_p->GetXaxis()->SetTitle("P (GeV)");

	//______x vs z
	//forward angle
	TH2F *hf_x_z=new TH2F("hf_x_z","x VS z ",Nbinx,0,0.7,Nbiny,0.25,0.75);
	hf_x_z->GetXaxis()->SetTitle("x_{bj}");
	hf_x_z->GetYaxis()->SetTitle("z");
	//large angle
	TH2F *hl_x_z=new TH2F("hl_x_z","x VS z ",Nbinx,0,0.7,Nbiny,0.25,0.75);
	hl_x_z->GetXaxis()->SetTitle("x_{bj}");
	hl_x_z->GetYaxis()->SetTitle("z");
	//total
	TH2F *h_x_z=new TH2F("h_x_z","x VS z ",Nbinx,0,0.7,Nbiny,0.25,0.75);
	h_x_z->GetXaxis()->SetTitle("x_{bj}");
	h_x_z->GetYaxis()->SetTitle("z");

	//______Q2 vs p 
	//forward angle
	TH2F *hf_Q2_p=new TH2F("hf_Q2_p","Q^{2} VS P ",Nbinx,0,7,Nbiny,0,9);
	hf_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hf_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_Q2_p=new TH2F("hl_Q2_p","Q^{2} VS P ",Nbinx,0,7,Nbiny,0,9);
	hl_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hl_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_Q2_p=new TH2F("h_Q2_p","Q^{2} VS P ",Nbinx,0,7,Nbiny,0,9);
	h_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	h_Q2_p->GetXaxis()->SetTitle("P (GeV)");

	//______Q2 vs x 
	//forward angle
	TH2F *hf_Q2_x=new TH2F("hf_Q2_x","Q^{2} VS x ",Nbinx,0,0.7,Nbiny,0,9);
	hf_Q2_x->GetXaxis()->SetTitle("x_{bj}");
	hf_Q2_x->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	//large angle
	TH2F *hl_Q2_x=new TH2F("hl_Q2_x","Q^{2} VS x ",Nbinx,0,0.7,Nbiny,0,9);
	hl_Q2_x->GetXaxis()->SetTitle("x_{bj}");
	hl_Q2_x->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	//total
	TH2F *h_Q2_x=new TH2F("h_Q2_x","Q^{2} VS x ",Nbinx,0,0.7,Nbiny,0,9);
	h_Q2_x->GetXaxis()->SetTitle("x_{bj}");
	h_Q2_x->GetYaxis()->SetTitle("Q^2(GeV^{2})");

	//______W vs x 
	//forward angle
	TH2F *hf_W_x=new TH2F("hf_W_x","W VS x ",Nbinx,0,0.7,Nbiny,2,4.5);
	hf_W_x->GetXaxis()->SetTitle("x_{bj}");
	hf_W_x->GetYaxis()->SetTitle("W (GeV)");
	//large angle
	TH2F *hl_W_x=new TH2F("hl_W_x","W VS x ",Nbinx,0,0.7,Nbiny,2,4.5);
	hl_W_x->GetXaxis()->SetTitle("x_{bj}");
	hl_W_x->GetYaxis()->SetTitle("W (GeV)");
	//total
	TH2F *h_W_x=new TH2F("h_W_x","W VS x ",Nbinx,0,0.7,Nbiny,2,4.5);
	h_W_x->GetXaxis()->SetTitle("x_{bj}");
	h_W_x->GetYaxis()->SetTitle("W (GeV)");

	//______W vs x 
	//forward angle
	TH2F *hf_Wp_x=new TH2F("hf_Wp_x","Wp VS x ",Nbinx,0,0.7,Nbiny,1.5,4);
	hf_Wp_x->GetXaxis()->SetTitle("x_{bj}");
	hf_Wp_x->GetYaxis()->SetTitle("Wp (GeV)");
	//large angle
	TH2F *hl_Wp_x=new TH2F("hl_Wp_x","Wp VS x ",Nbinx,0,0.7,Nbiny,1.5,4);
	hl_Wp_x->GetXaxis()->SetTitle("x_{bj}");
	hl_Wp_x->GetYaxis()->SetTitle("Wp (GeV)");
	//total
	TH2F *h_Wp_x=new TH2F("h_Wp_x","Wp VS x ",Nbinx,0,0.7,Nbiny,1.5,4);
	h_Wp_x->GetXaxis()->SetTitle("x_{bj}");
	h_Wp_x->GetYaxis()->SetTitle("Wp (GeV)");

	//_____p vs theta
	TH2F *hf_p_theta=new TH2F("hf_p_theta","P_{hadron} VS #theta_{hadron}  (forward angle)",Nbin_theta,7.0,30.0,Nbin_p,0,11);
	//forward angle
	hf_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	hf_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	//large angle
	TH2F *hl_p_theta=new TH2F("hl_p_theta","P_{hadron} VS #theta_{hadron} (large angle)",Nbin_theta,7.0,30.0,Nbin_p,0,11);
	hl_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	hl_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	//total
	TH2F *h_p_theta=new TH2F("h_p_theta","P_{hadron} VS #theta_{hadron} (large+forward angle)",Nbin_theta,7.0,30.0,Nbin_p,0,11);
	h_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	h_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");

	//_____TH1F of theta
	//forward angle
	TH1F *hf_theta=new TH1F("hf_theta","rate as function of #theta_{eletron}  (forward angle)",Nbin_theta,7.0,30.0);
	hf_theta->GetYaxis()->SetTitle("rate");
	hf_theta->GetXaxis()->SetTitle("#theta_{eletron}");
	//large angle
	TH1F *hl_theta=new TH1F("hl_theta","rate as function of #theta_{eletron} (large angle)",Nbin_theta,7.0,30.0);
	hl_theta->GetYaxis()->SetTitle("rate");
	hl_theta->GetXaxis()->SetTitle("#theta_{eletron}");
	//total
	TH1F *h_theta=new TH1F("h_theta","rate as function of #theta_{eletron} (large+forward angle)",Nbin_theta,7.0,30.0);
	h_theta->GetYaxis()->SetTitle("rate");
	h_theta->GetXaxis()->SetTitle("#theta_{eletron}");
	/*}}}*/

	double rate_forward=0;
	double rate_large=0;
	double ele_rate_forward=0;
	double ele_rate_large=0;
	double had_rate=0;

	double rate_forward_p_g_4[5]={0};
	double rate_large_p_g_4[5]={0};
	double rate_forward_p_l_4=0;
	double rate_large_p_l_4=0;
	int ele_theta_bin= 0;
	int ele_p_bin=0;
	double ele_forward_acceptance= 0.0;
	double ele_large_acceptance= 0.0;

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		ele_theta_bin= 0;
		ele_p_bin=0;
		ele_forward_acceptance= 0.0;
		ele_large_acceptance= 0.0;

		if(mom_ele>=1.0&&mom_ele<=10.0&&W>=2.3&&Wp>=1.6&&z>0.3&&z<0.7){//any additional cuts should be added in here
			ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
			ele_p_bin=int(mom_ele/0.05)+1;      //0.05 GeV per bin for mom
			ele_forward_acceptance=accep_ele_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			ele_large_acceptance=accep_ele_large->GetBinContent(ele_theta_bin,ele_p_bin);

			int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
			int hadron_p_bin=int(mom_had/0.05)+1;     //0.05 GeV per bin for mom
			double hadron_acceptance=accep_had->GetBinContent(hadron_theta_bin,hadron_p_bin);
			double event_weight=0;   //depend on hadron charge

			if(hadron_flag==1){  //doing hp
				event_weight=weight_hp;
			}else if(hadron_flag==2){  //doing hm
				event_weight=weight_hm;
			}

			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

			if(mom_ele<1.0||theta_ele*DEG>14.5||theta_ele*DEG<8.0)//GeV, CLEO
			//if(mom_ele<1.0||theta_ele*DEG>12.0||theta_ele*DEG<6.6)//GeV, CDF
				forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(mom_ele<3.5||theta_ele*DEG<15.5||theta_ele*DEG>24)//GeV,CLEO
			//if(mom_ele<3.5||theta_ele*DEG<12.5||theta_ele*DEG>22)//GeV,CDF
				large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV

			/*Rates{{{*/	
			if( mom_ele>=1.0){
				ele_rate_forward+=event_weight*ele_forward_acceptance;
			}
			if( mom_ele>=3.5){
				ele_rate_large+=event_weight*ele_large_acceptance;
			}
			if( mom_had>=0.0){
				had_rate+=event_weight*hadron_acceptance;
			}
	
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
			
			/*fill histogram here{{{*/
			//z VS pt
			hf_z_pt->Fill(z,pt,event_weight*forward_acceptance);
			hl_z_pt->Fill(z,pt,event_weight*large_acceptance);
			h_z_pt->Fill(z,pt,event_weight*(forward_acceptance+large_acceptance));

			//x vs pt
			hf_x_pt->Fill(x,pt,event_weight*forward_acceptance);
			hl_x_pt->Fill(x,pt,event_weight*large_acceptance);
			h_x_pt->Fill(x,pt,event_weight*(forward_acceptance+large_acceptance));

			//Q2 vs pt
			hf_Q2_pt->Fill(Q2,pt,event_weight*forward_acceptance);
			hl_Q2_pt->Fill(Q2,pt,event_weight*large_acceptance);
			h_Q2_pt->Fill(Q2,pt,event_weight*(forward_acceptance+large_acceptance));

			//pt vs p
			hf_pt_p->Fill(mom_had,pt,event_weight*forward_acceptance);
			hl_pt_p->Fill(mom_had,pt,event_weight*large_acceptance);
			h_pt_p->Fill(mom_had,pt,event_weight*(forward_acceptance+large_acceptance));

			//x vs p
			hf_x_p->Fill(mom_had,x,event_weight*forward_acceptance);
			hl_x_p->Fill(mom_had,x,event_weight*large_acceptance);
			h_x_p->Fill(mom_had,x,event_weight*(forward_acceptance+large_acceptance));

			//x vs z
			hf_x_z->Fill(x,z,event_weight*forward_acceptance);
			hl_x_z->Fill(x,z,event_weight*large_acceptance);
			h_x_z->Fill(x,z,event_weight*(forward_acceptance+large_acceptance));

			//Q2 vs p
			hf_Q2_p->Fill(mom_had,Q2,event_weight*forward_acceptance);
			hl_Q2_p->Fill(mom_had,Q2,event_weight*large_acceptance);
			h_Q2_p->Fill(mom_had,Q2,event_weight*(forward_acceptance+large_acceptance));

			//Q2 vs x
			hf_Q2_x->Fill(x,Q2,event_weight*forward_acceptance);
			hl_Q2_x->Fill(x,Q2,event_weight*large_acceptance);
			h_Q2_x->Fill(x,Q2,event_weight*(forward_acceptance+large_acceptance));

			//W vs x
			hf_W_x->Fill(x,W,event_weight*forward_acceptance);
			hl_W_x->Fill(x,W,event_weight*large_acceptance);
			h_W_x->Fill(x,W,event_weight*(forward_acceptance+large_acceptance));

			//Wp vs x
			hf_Wp_x->Fill(x,Wp,event_weight*forward_acceptance);
			hl_Wp_x->Fill(x,Wp,event_weight*large_acceptance);
			h_Wp_x->Fill(x,Wp,event_weight*(forward_acceptance+large_acceptance));

			//p vs theta
			hf_p_theta->Fill(theta_ele*DEG,mom_ele,event_weight*forward_acceptance);
			hl_p_theta->Fill(theta_ele*DEG,mom_ele,event_weight*large_acceptance);
			h_p_theta->Fill(theta_ele*DEG,mom_ele, event_weight*(forward_acceptance+large_acceptance));

			//TH1F theta
			hf_theta->Fill(theta_ele*DEG,event_weight*forward_acceptance);
			hl_theta->Fill(theta_ele*DEG,event_weight*large_acceptance);
			h_theta->Fill(theta_ele*DEG,event_weight*(forward_acceptance+large_acceptance));
			/*}}}*/
		}

	}// events loop ends here


	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;

	/*Print&Save{{{*/
	ofstream outf(Form("%s_%s_he3_rate_8.dat", particle.Data(), charge.Data()));

	outf<<"______Single integral rate____________________"<<endl;
	outf<<"eletron rate_forward: "<<ele_rate_forward*Lumi*KHz*nBcm2<<endl;
	outf<<"eletron rate_large: "<<ele_rate_large*Lumi*KHz*nBcm2<<endl;
	outf<<"hadron rate: "<<had_rate*Lumi*KHz*nBcm2<<endl;

	outf<<"______forward and large angle integral rate____________________"<<endl;
	outf<<"rate_forward: "<<rate_forward*Lumi*KHz*nBcm2<<endl;
	outf<<"rate_large: "<<rate_large*Lumi*KHz*nBcm2<<endl;

	outf<<"______________ p~[3.5,4.0) GeV rate _________________________________"<<endl;
	outf<<"forward: "<<rate_forward_p_g_4[4]*Lumi*KHz*nBcm2<<endl;
	outf<<"large:   "<<rate_large_p_g_4[4]*Lumi*KHz*nBcm2<<endl;

	outf<<"______________ p~[4,4.5) GeV rate _________________________________"<<endl;
	outf<<"forward: "<<rate_forward_p_g_4[0]*Lumi*KHz*nBcm2<<endl;
	outf<<"large:   "<<rate_large_p_g_4[0]*Lumi*KHz*nBcm2<<endl;

	outf<<"______________ p~[4.5,5) GeV rate _________________________________"<<endl;
	outf<<"forward: "<<rate_forward_p_g_4[1]*Lumi*KHz*nBcm2<<endl;
	outf<<"large:   "<<rate_large_p_g_4[1]*Lumi*KHz*nBcm2<<endl;

	outf<<"______________ p~[5,5.5) GeV rate _________________________________"<<endl;
	outf<<"forward: "<<rate_forward_p_g_4[2]*Lumi*KHz*nBcm2<<endl;
	outf<<"large:   "<<rate_large_p_g_4[2]*Lumi*KHz*nBcm2<<endl;

	outf<<"______________ p~[5.5,6) GeV rate _________________________________"<<endl;
	outf<<"forward: "<<rate_forward_p_g_4[3]*Lumi*KHz*nBcm2<<endl;
	outf<<"large:   "<<rate_large_p_g_4[3]*Lumi*KHz*nBcm2<<endl;

	outf<<"______________ 0.8<p<3.5 GeV rate __________________________________"<<endl;
	outf<<"forward: "<<rate_forward_p_l_4*Lumi*KHz*nBcm2<<endl;
	outf<<"large:   "<<rate_large_p_l_4*Lumi*KHz*nBcm2<<endl;
	outf<<"_______________________________________________________________"<<endl;


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

	/*Plot-New{{{*/

	TCanvas *can_forward=new TCanvas("can_forward","forward angle", 600,800);
	can_forward->Divide(2,3);
	can_forward->cd(1);
	hf_Q2_x->Draw("colz");
	can_forward->cd(2);
	hf_x_pt->Draw("colz");
	can_forward->cd(3);
	hf_W_x->Draw("colz");
	can_forward->cd(4);
	hf_x_z->Draw("colz");
	can_forward->cd(5);
	hf_Wp_x->Draw("colz");
	can_forward->cd(6);
	hf_z_pt->Draw("COLZ");
	can_forward->Print(Form("CLEO_NewXS_%s_%s_Forward.png", particle.Data(), charge.Data()));

	TCanvas *can_large=new TCanvas("can_large","large angle", 600,800);
	can_large->Divide(2,3);
	can_large->cd(1);
	hl_Q2_x->Draw("colz");
	can_large->cd(2);
	hl_x_pt->Draw("colz");
	can_large->cd(3);
	hl_W_x->Draw("colz");
	can_large->cd(4);
	hl_x_z->Draw("colz");
	can_large->cd(5);
	hl_Wp_x->Draw("colz");
	can_large->cd(6);
	hl_z_pt->Draw("COLZ");
	can_large->Print(Form("CLEO_NewXS_%s_%s_Large.png", particle.Data(), charge.Data()));

	TCanvas *can_total=new TCanvas("can_total","total", 600,800);
	can_total->Divide(2,3);
	can_total->cd(1);
	//	h_Q2_x->Draw("colz");
	hf_Q2_x->Draw();
	hl_Q2_x->SetMarkerColor(8); hl_Q2_x->Draw("same");
	can_total->cd(2);
	//	h_x_pt->Draw("colz");
	hf_x_pt->Draw();
	hl_x_pt->SetMarkerColor(8); hl_x_pt->Draw("same");
	can_total->cd(3);
	//	h_W_x->Draw("colz");
	hf_W_x->Draw();
	hl_W_x->SetMarkerColor(8); hl_W_x->Draw("same");
	can_total->cd(4);
	//	h_x_z->Draw("colz");
	hf_x_z->Draw();
	hl_x_z->SetMarkerColor(8); hl_x_z->Draw("same");
	can_total->cd(5);
	//	h_Wp_x->Draw("colz");
	hf_Wp_x->Draw();
	hl_Wp_x->SetMarkerColor(8); hl_Wp_x->Draw("same");
	can_total->cd(6);
	//	h_z_pt->Draw("COLZ");
	hf_z_pt->Draw();
	hl_z_pt->SetMarkerColor(8); hl_z_pt->Draw("same");
	can_total->Print(Form("CLEO_NewXS_%s_%s_Total.png", particle.Data(), charge.Data()));

	TCanvas *c1=new TCanvas("c1","theta and P_{hadron} plots");
	c1->Divide(2,2);
	c1->cd(1);
	hf_p_theta->Draw("COLZ");
	c1->cd(2);
	hl_p_theta->Draw("COLZ");
	c1->cd(3);
	hf_theta->Draw();
	c1->cd(4);
	hl_theta->Draw();
	c1->Print(Form("CLEO_NewXS_%s_%s_ThetaP.png", particle.Data(), charge.Data()));

	/*}}}*/

	/*Save Histo{{{*/
	TString outputfilename="test.root";
	if(particle=="pion"){
		if(hadron_flag==1){
			outputfilename = "pip_histo.root";
		}
		if(hadron_flag==2){
			outputfilename = "pim_histo.root";
		}
	}

	if(particle=="kaon"){
		if(hadron_flag==1){
			outputfilename = "kp_histo.root";
		}
		if(hadron_flag==2){
			outputfilename = "km_histo.root";
		}
	}

	TFile *fout = new TFile(outputfilename.Data(),"recreate");
	fout->cd();
	h_z_pt->Write();
	h_x_pt->Write();
	h_Q2_pt->Write();
	h_pt_p->Write();
	h_x_p->Write();
	h_x_z->Write();
	h_Q2_p->Write();
	h_Q2_x->Write();
	h_Wp_x->Write();
	h_W_x->Write();
	h_p_theta->Write();
	h_theta->Write();

	hl_z_pt->Write();
	hl_x_pt->Write();
	hl_Q2_pt->Write();
	hl_pt_p->Write();
	hl_x_p->Write();
	hl_x_z->Write();
	hl_Q2_p->Write();
	hl_Q2_x->Write();
	hl_Wp_x->Write();
	hl_W_x->Write();
	hl_p_theta->Write();
	hl_theta->Write();

	hf_z_pt->Write();
	hf_x_pt->Write();
	hf_Q2_pt->Write();
	hf_pt_p->Write();
	hf_x_p->Write();
	hf_x_z->Write();
	hf_Q2_p->Write();
	hf_Q2_x->Write();
	hf_Wp_x->Write();
	hf_W_x->Write();
	hf_p_theta->Write();
	hf_theta->Write();

	fout->Close();
	/*}}}*/
}

