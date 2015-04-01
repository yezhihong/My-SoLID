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


void plot(TString particle,int hadron_flag){
	//particle: pion    kaon
	//hadron_flag   1 for hp      2 for hm

	TString charge = "A";
	if(hadron_flag==1)
		charge = "p";
	else if(hadron_flag==2)
		charge= "m";
	else
		cerr<<"********* Error"<<endl;

	gStyle->SetOptStat(0);
	gStyle->SetLabelSize(0.07,"X");
	gStyle->SetLabelSize(0.07,"Y");
	gStyle->SetLabelSize(0.07,"Z");
	gStyle->SetTitleSize(0.085,"X");
	gStyle->SetTitleSize(0.085,"Y");
	gStyle->SetTitleSize(0.085,"Z");
	gStyle->SetTitleOffset(0.7,"X");
	gStyle->SetTitleOffset(0.7,"Y");
	gStyle->SetTitleSize(0.08,"t");


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
	double mom_ele,mom_had,theta_ele,theta_had,phi_ele,phi_had;  // acceptance calculation input
	double weight_hp,weight_hm;

	T->SetBranchAddress("Q2",&Q2);
	T->SetBranchAddress("x",&x);
	T->SetBranchAddress("pt",&pt);
	T->SetBranchAddress("W",&W);
	T->SetBranchAddress("Wp",&Wp);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("mom_had",&mom_had);
	T->SetBranchAddress("mom_ele",&mom_ele);
	T->SetBranchAddress("phi_had",&phi_had);
	T->SetBranchAddress("phi_ele",&phi_ele);
	T->SetBranchAddress("theta_had",&theta_had);
	T->SetBranchAddress("theta_ele",&theta_ele);
	T->SetBranchAddress("weight_hp",&weight_hp);   //unit in nbar
	T->SetBranchAddress("weight_hm",&weight_hm);

	Long64_t N_entries=T->GetEntries();
	cout<<"total generated events number: "<<N_entries<<endl;

	//deal with acceptance files 
	//read in acceptance histograms
	
	//electron acceptance
		
//	TFile *file_negative=new TFile("./acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
//	TH2F *accep_had=0;
//	accep_had=(TH2F*)accep_ele_forward->Clone("acceptance_hadron");//pion is normal acceptance histogram, the same as electron acceptance

	TFile *file_negative=new TFile("Xin_CDF_Acceptance_output.root");
	TH2F *accep_had=(TH2F*)file_negative->Get("acceptance");           //negative particle large angle acceptance

	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
		
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
	TH2F *hf_p_theta=new TH2F("hf_p_theta","P_{#pi-} VS #theta_{#pi^{-}}  (forward angle)",Nbin_theta,5.0,30.0,Nbin_p,0,11);
	//forward angle
	hf_p_theta->GetYaxis()->SetTitle("P_{#pi-}");
	hf_p_theta->GetXaxis()->SetTitle("#theta_{#pi-}");
	hf_p_theta->GetXaxis()->SetTitleSize(0.085);
	hf_p_theta->GetYaxis()->SetTitleSize(0.085);
	hf_p_theta->GetXaxis()->CenterTitle(1);
	hf_p_theta->GetYaxis()->CenterTitle(1);

	//large angle
	TH2F *hl_p_theta=new TH2F("hl_p_theta","P_{#pi-} VS #theta_{#pi^{-}} (large angle)",Nbin_theta,5.0,30.0,Nbin_p,0,11);
	hl_p_theta->GetYaxis()->SetTitle("P_{#pi-}");
	hl_p_theta->GetXaxis()->SetTitle("#theta_{#pi-}");
	hl_p_theta->GetXaxis()->SetTitleSize(0.085);
	hl_p_theta->GetYaxis()->SetTitleSize(0.085);
	hl_p_theta->GetXaxis()->CenterTitle(1);
	hl_p_theta->GetYaxis()->CenterTitle(1);
		//total
	TH2F *h_p_theta=new TH2F("h_p_theta","P_{#pi-} VS #theta_{#pi^{-}} (large+forward angle)",Nbin_theta,5.0,30.0,Nbin_p,0,11);
	h_p_theta->GetYaxis()->SetTitle("P_{#pi-}");
	h_p_theta->GetXaxis()->SetTitle("#theta_{#pi-}");
	h_p_theta->GetXaxis()->SetTitleSize(0.085);
	h_p_theta->GetYaxis()->SetTitleSize(0.085);
	h_p_theta->GetXaxis()->CenterTitle(1);
	h_p_theta->GetYaxis()->CenterTitle(1);
	

	//_____TH1F of theta
	//forward angle
	TH1F *hf_theta=new TH1F("hf_theta","rate as function of #theta_{#pi^{-}}  (forward angle)",Nbin_theta,5.0,30.0);
	hf_theta->GetYaxis()->SetTitle("rate");
	hf_theta->GetXaxis()->SetTitle("#theta_{#pi^{-}}");
	hf_theta->GetXaxis()->SetTitleSize(0.085);
	hf_theta->GetYaxis()->SetTitleSize(0.085);
	hf_theta->GetXaxis()->CenterTitle(1);
	hf_theta->GetYaxis()->CenterTitle(1);
	
	//large angle
	TH1F *hl_theta=new TH1F("hl_theta","rate as function of #theta_{#pi^{-}} (large angle)",Nbin_theta,5.0,30.0);
	hl_theta->GetYaxis()->SetTitle("rate");
	hl_theta->GetXaxis()->SetTitle("#theta_{#pi^{-}}");
	hl_theta->GetXaxis()->SetTitleSize(0.085);
	hl_theta->GetYaxis()->SetTitleSize(0.085);
	hl_theta->GetXaxis()->CenterTitle(1);
	hl_theta->GetYaxis()->CenterTitle(1);
	
	//total
	TH1F *h_theta=new TH1F("h_theta","rate as function of #theta_{#pi^{-}} (large+forward angle)",Nbin_theta,5.0,30.0);
	h_theta->GetYaxis()->SetTitle("rate");
	h_theta->GetXaxis()->SetTitle("#theta_{#pi^{-}}");
	h_theta->GetXaxis()->SetTitleSize(0.085);
	h_theta->GetYaxis()->SetTitleSize(0.085);
	h_theta->GetXaxis()->CenterTitle(1);
	h_theta->GetYaxis()->CenterTitle(1);
	
	//_____p
	TH1F *hf_p=new TH1F("hf_p","P_{#pi^{-}} (forward angle)",Nbin_p,0,11);
	//forward angle	
	hf_p->GetYaxis()->SetTitle("rate");
	hf_p->GetXaxis()->SetTitle("P_{#pi^{-}}");
	hf_p->GetXaxis()->SetTitleSize(0.085);
	hf_p->GetYaxis()->SetTitleSize(0.085);
	hf_p->GetXaxis()->CenterTitle(1);
	hf_p->GetYaxis()->CenterTitle(1);
	
TH1F *hl_p=new TH1F("hl_p","P_{#pi^{-}} (large angle)",Nbin_p,0,11);
	//forward angle	
	hl_p->GetYaxis()->SetTitle("rate");
	hl_p->GetXaxis()->SetTitle("P_{#pi^{-}}");
	hl_p->GetXaxis()->SetTitleSize(0.085);
	hl_p->GetYaxis()->SetTitleSize(0.085);
	hl_p->GetXaxis()->CenterTitle(1);
	hl_p->GetYaxis()->CenterTitle(1);
	
TH1F *h_p=new TH1F("h_p","P_{#pi^{-}} (both angles)",Nbin_p,0,11);
	//forward angle	
	h_p->GetYaxis()->SetTitle("rate");
	h_p->GetXaxis()->SetTitle("P_{#pi^{-}}");
	h_p->GetXaxis()->SetTitleSize(0.085);
	h_p->GetYaxis()->SetTitleSize(0.085);
	h_p->GetXaxis()->CenterTitle(1);
	h_p->GetYaxis()->CenterTitle(1);
	
	//_____phi
	TH1F *hf_phi=new TH1F("hf_phi","#phi_{#pi^{-}} (forward angle)",Nbin_p,-200,200);
	//forward angle	
	hf_phi->GetYaxis()->SetTitle("rate");
	hf_phi->GetXaxis()->SetTitle("#phi_{#pi^{-}}");
	hf_phi->GetXaxis()->SetTitleSize(0.085);
	hf_phi->GetYaxis()->SetTitleSize(0.085);
	hf_phi->GetXaxis()->CenterTitle(1);
	hf_phi->GetYaxis()->CenterTitle(1);
	
	TH1F *hl_phi=new TH1F("hl_phi","#phi_{#pi^{-}} (large angle)",Nbin_p,-200,200);
	//forward angle	
	hl_phi->GetYaxis()->SetTitle("rate");
	hl_phi->GetXaxis()->SetTitle("#phi_{#pi^{-}}");
	hl_phi->GetXaxis()->SetTitleSize(0.085);
	hl_phi->GetYaxis()->SetTitleSize(0.085);
	hl_phi->GetXaxis()->CenterTitle(1);
	hl_phi->GetYaxis()->CenterTitle(1);
	
	TH1F *h_phi=new TH1F("h_phi","#phi_{#pi^{-}} (both angles)",Nbin_p,-200,200);
	//forward angle	
	h_phi->GetYaxis()->SetTitle("rate");
	h_phi->GetXaxis()->SetTitle("#phi_{#pi^{-}}");
	h_phi->GetXaxis()->SetTitleSize(0.085);
	h_phi->GetYaxis()->SetTitleSize(0.085);
	h_phi->GetXaxis()->CenterTitle(1);
	h_phi->GetYaxis()->CenterTitle(1);
	
	TH1F *ef_phi=new TH1F("ef_phi","#phi_{e^{-}} (forward angle)",Nbin_p,-200,200);
	//forward angle	
	ef_phi->GetYaxis()->SetTitle("rate");
	ef_phi->GetXaxis()->SetTitle("#phi_{e^{-}}");
	ef_phi->GetXaxis()->SetTitleSize(0.085);
	ef_phi->GetYaxis()->SetTitleSize(0.085);
	ef_phi->GetXaxis()->CenterTitle(1);
	ef_phi->GetYaxis()->CenterTitle(1);
	
	TH1F *el_phi=new TH1F("el_phi","#phi_{e^{-}} (large angle)",Nbin_p,-200,200);
	//forward angle	
	el_phi->GetYaxis()->SetTitle("rate");
	el_phi->GetXaxis()->SetTitle("#phi_{e^{-}}");
	el_phi->GetXaxis()->SetTitleSize(0.085);
	el_phi->GetYaxis()->SetTitleSize(0.085);
	el_phi->GetXaxis()->CenterTitle(1);
	el_phi->GetYaxis()->CenterTitle(1);
	
	TH1F *e_phi=new TH1F("e_phi","#phi_{e^{-}} (both angles)",Nbin_p,-200,200);
	//forward angle	
	e_phi->GetYaxis()->SetTitle("rate");
	e_phi->GetXaxis()->SetTitle("#phi_{e^{-}}");
	e_phi->GetXaxis()->SetTitleSize(0.085);
	e_phi->GetYaxis()->SetTitleSize(0.085);
	e_phi->GetXaxis()->CenterTitle(1);
	e_phi->GetYaxis()->CenterTitle(1);
	
	//__Electron___p vs theta
	TH2F *ef_p_theta=new TH2F("ef_p_theta","P_{e^{-}} VS #theta_{e^{-}}  (forward angle)",Nbin_theta,5.0,30.0,Nbin_p,0,11);
	//forward angle
	ef_p_theta->GetYaxis()->SetTitle("P_{e^{-}}");
	ef_p_theta->GetXaxis()->SetTitle("#theta_{e^{-}}");
	ef_p_theta->GetXaxis()->SetTitleSize(0.085);
	ef_p_theta->GetYaxis()->SetTitleSize(0.085);
	ef_p_theta->GetXaxis()->CenterTitle(1);
	ef_p_theta->GetYaxis()->CenterTitle(1);
	
//large angle
	TH2F *el_p_theta=new TH2F("el_p_theta","P_{e^{-}} VS #theta_{e^{-}} (large angle)",Nbin_theta,5.0,30.0,Nbin_p,0,11);
	el_p_theta->GetYaxis()->SetTitle("P_{e^{-}}");
	el_p_theta->GetXaxis()->SetTitle("#theta_{e^{-}}");
	el_p_theta->GetXaxis()->SetTitleSize(0.085);
	el_p_theta->GetYaxis()->SetTitleSize(0.085);
	el_p_theta->GetXaxis()->CenterTitle(1);
	el_p_theta->GetYaxis()->CenterTitle(1);
	
	//total
	TH2F *e_p_theta=new TH2F("e_p_theta","P_{e^{-}} VS #theta_{e^{-}} (large+forward angle)",Nbin_theta,5.0,30.0,Nbin_p,0,11);
	e_p_theta->GetYaxis()->SetTitle("P_{e^{-}}");
	e_p_theta->GetXaxis()->SetTitle("#theta_{e^{-}}");
	e_p_theta->GetXaxis()->SetTitleSize(0.085);
	e_p_theta->GetYaxis()->SetTitleSize(0.085);
	e_p_theta->GetXaxis()->CenterTitle(1);
	e_p_theta->GetYaxis()->CenterTitle(1);
	

	//_____TH1F of theta
	//forward angle
	TH1F *ef_theta=new TH1F("ef_theta","rate as function of #theta_{e^{-}}  (forward angle)",Nbin_theta,5.0,30.0);
	ef_theta->GetYaxis()->SetTitle("rate");
	ef_theta->GetXaxis()->SetTitle("#theta_{e^{-}}");
	ef_theta->GetXaxis()->SetTitleSize(0.085);
	ef_theta->GetYaxis()->SetTitleSize(0.085);
	ef_theta->GetXaxis()->CenterTitle(1);
	ef_theta->GetYaxis()->CenterTitle(1);
	
//large angle
	TH1F *el_theta=new TH1F("el_theta","rate as function of #theta_{e^{-}} (large angle)",Nbin_theta,5.0,30.0);
	el_theta->GetYaxis()->SetTitle("rate");
	el_theta->GetXaxis()->SetTitle("#theta_{e^{-}}");
	el_theta->GetXaxis()->SetTitleSize(0.085);
	el_theta->GetYaxis()->SetTitleSize(0.085);
	el_theta->GetXaxis()->CenterTitle(1);
	el_theta->GetYaxis()->CenterTitle(1);
	
	//total
	TH1F *e_theta=new TH1F("e_theta","rate as function of #theta_{e^{-}} (large+forward angle)",Nbin_theta,5.0,30.0);
	e_theta->GetYaxis()->SetTitle("rate");
	e_theta->GetXaxis()->SetTitle("#theta_{e^{-}}");
	e_theta->GetXaxis()->SetTitleSize(0.085);
	e_theta->GetYaxis()->SetTitleSize(0.085);
	e_theta->GetXaxis()->CenterTitle(1);
	e_theta->GetYaxis()->CenterTitle(1);
	

	//_____p
	TH1F *ef_p=new TH1F("ef_p","P_{e^{-}} (forward angle)",Nbin_p,0,11);
	//forward angle	
	ef_p->GetYaxis()->SetTitle("rate");
	ef_p->GetXaxis()->SetTitle("P_{e^{-}}");
	ef_p->GetXaxis()->SetTitleSize(0.085);
	ef_p->GetYaxis()->SetTitleSize(0.085);
	ef_p->GetXaxis()->CenterTitle(1);
	ef_p->GetYaxis()->CenterTitle(1);
	
	TH1F *el_p=new TH1F("el_p","P_{e^{-}} (large angle)",Nbin_p,0,11);
	//forward angle	
	el_p->GetYaxis()->SetTitle("rate");
	el_p->GetXaxis()->SetTitle("P_{e^{-}}");
	el_p->GetXaxis()->SetTitleSize(0.085);
	el_p->GetYaxis()->SetTitleSize(0.085);
	el_p->GetXaxis()->CenterTitle(1);
	el_p->GetYaxis()->CenterTitle(1);
	
	TH1F *e_p=new TH1F("e_p","P_{e^{-}} (both angles)",Nbin_p,0,11);
	//forward angle	
	e_p->GetYaxis()->SetTitle("rate");
	e_p->GetXaxis()->SetTitle("P_{e^{-}}");
	e_p->GetXaxis()->SetTitleSize(0.085);
	e_p->GetYaxis()->SetTitleSize(0.085);
	e_p->GetXaxis()->CenterTitle(1);
	e_p->GetYaxis()->CenterTitle(1);
	/*}}}*/

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

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		ele_theta_bin= 0;
		ele_p_bin=0;
		ele_forward_acceptance= 0.0;
		ele_large_acceptance= 0.0;

		if(theta_ele*DEG>=4.0&&theta_ele*DEG<=36.0&&mom_ele>=1.0&&mom_ele<=11.0&&W>=2.3&&Wp>=1.6){//any additional cuts should be added in here
		//if(W>=2.3&&Wp>=1.6){//any additional cuts should be added in here
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

	//		if(ele_forward_acceptance>1.)
	//			cerr<<Form("*** Ele Forward Accept = %f for theta=%f, p=%f", ele_forward_acceptance, theta_ele*DEG, mom_ele)<<endl;
	//		if(ele_large_acceptance>1.)
	//			cerr<<Form("*** Ele Large Accept = %f", ele_large_acceptance)<<endl;
	//		if(hadron_acceptance>1.)
	//			cerr<<Form("*** Had Accept = %f", hadron_acceptance)<<endl;

			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

		//	double forward_acceptance= 1.0;
		//	double large_acceptance= 1.0;

			event_weight = 1.0;

			if(mom_ele<1.0)//GeV
				forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
			if(mom_ele<3.5)//GeV
				large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV

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
			hf_p_theta->Fill(theta_had*DEG,mom_had,event_weight*forward_acceptance);
			hl_p_theta->Fill(theta_had*DEG,mom_had,event_weight*large_acceptance);
			h_p_theta->Fill(theta_had*DEG,mom_had, event_weight*(forward_acceptance+large_acceptance));

			//TH1F theta
			hf_theta->Fill(theta_had*DEG,event_weight*forward_acceptance);
			hl_theta->Fill(theta_had*DEG,event_weight*large_acceptance);
			h_theta->Fill(theta_had*DEG,event_weight*(forward_acceptance+large_acceptance));

			//TH1F p 
			h_p->Fill(mom_had, event_weight*(forward_acceptance+large_acceptance));
			hl_p->Fill(mom_had, event_weight*(large_acceptance));
			hf_p->Fill(mom_had, event_weight*(forward_acceptance));

			//TH1F phi
			h_phi->Fill(phi_had*DEG, event_weight*(forward_acceptance+large_acceptance));
			hl_phi->Fill(phi_had*DEG, event_weight*(large_acceptance));
			hf_phi->Fill(phi_had*DEG, event_weight*(forward_acceptance));

			//p vs theta
			e_p_theta->Fill(theta_ele*DEG,mom_ele, event_weight*(forward_acceptance+large_acceptance));
			ef_p_theta->Fill(theta_ele*DEG,mom_ele,event_weight*forward_acceptance);
			el_p_theta->Fill(theta_ele*DEG,mom_ele,event_weight*large_acceptance);

			//TH1F theta
			e_theta->Fill(theta_ele*DEG,event_weight*(forward_acceptance+large_acceptance));
			ef_theta->Fill(theta_ele*DEG,event_weight*forward_acceptance);
			el_theta->Fill(theta_ele*DEG,event_weight*large_acceptance);

			//TH1F p 
			e_p->Fill(mom_ele, event_weight*(forward_acceptance+large_acceptance));
			el_p->Fill(mom_ele, event_weight*(large_acceptance));
			ef_p->Fill(mom_ele, event_weight*(forward_acceptance));

			//TH1F phi
			e_phi->Fill(phi_ele*DEG, event_weight*(forward_acceptance+large_acceptance));
			el_phi->Fill(phi_ele*DEG, event_weight*(large_acceptance));
			ef_phi->Fill(phi_ele*DEG, event_weight*(forward_acceptance));
			/*}}}*/

		}

	}// events loop ends here


	const double Lumi = 3.0e36; // cm-2*s-1
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;

	/*Plot-New{{{*/
/*
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
	can_forward->Print(Form("NewACC_NewXS_%s_%s_Forward_NoP.pdf", particle.Data(), charge.Data()));
	can_forward->Print(Form("NewACC_NewXS_%s_%s_Forward_NoP.png", particle.Data(), charge.Data()));

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
	can_large->Print(Form("NewACC_NewXS_%s_%s_Large_NoP.pdf", particle.Data(), charge.Data()));
	can_large->Print(Form("NewACC_NewXS_%s_%s_Large_NoP.png", particle.Data(), charge.Data()));

	TCanvas *can_total=new TCanvas("can_total","total", 600,800);
	can_total->Divide(2,3);
	can_total->cd(1);
	//	h_Q2_x->Draw("colz");
	hf_Q2_x->SetMarkerColor(4); hf_Q2_x->Draw();
	hl_Q2_x->SetMarkerColor(8); hl_Q2_x->Draw("same");
	can_total->cd(2);
	//	h_x_pt->Draw("colz");
	hf_x_pt->SetMarkerColor(4); hf_x_pt->Draw();
	hl_x_pt->SetMarkerColor(8); hl_x_pt->Draw("same");
	can_total->cd(3);
	//	h_W_x->Draw("colz");
	hf_W_x->SetMarkerColor(4); hf_W_x->Draw();
	hl_W_x->SetMarkerColor(8); hl_W_x->Draw("same");
	can_total->cd(4);
	//	h_x_z->Draw("colz");
	hf_x_z->SetMarkerColor(4); hf_x_z->Draw();
	hl_x_z->SetMarkerColor(8); hl_x_z->Draw("same");
	can_total->cd(5);
	//	h_Wp_x->Draw("colz");
	hf_Wp_x->SetMarkerColor(4); hf_Wp_x->Draw();
	hl_Wp_x->SetMarkerColor(8); hl_Wp_x->Draw("same");
	can_total->cd(6);
	//	h_z_pt->Draw("COLZ");
	hf_z_pt->SetMarkerColor(4);hf_z_pt->Draw();
	hl_z_pt->SetMarkerColor(8); hl_z_pt->Draw("same");
	can_total->Print(Form("NewACC_NewXS_%s_%s_Total_NoP.pdf", particle.Data(), charge.Data()));
	can_total->Print(Form("NewACC_NewXS_%s_%s_Total_NoP.png", particle.Data(), charge.Data()));
*/
	TCanvas *c1=new TCanvas("c1","theta and P_{#pi^{-}} plots", 1000,800);
	c1->Divide(2,4);
	int i=1;
	c1->cd(i++);
	hf_p_theta->Draw("COLZ");
	c1->cd(i++);
	hl_p_theta->Draw("COLZ");
	c1->cd(i++);
	hf_p->Draw();
	c1->cd(i++);
	hl_p->Draw();
	c1->cd(i++);
	hf_theta->Draw();
	c1->cd(i++);
	hl_theta->Draw();
	c1->cd(i++);
	hf_phi->Draw();
	c1->cd(i++);
	hl_phi->Draw();
	//c1->Print("New_Kin_pim_NewAcc.png");
	//c1->Print("New_Kin_pim_NewAcc_NoPhy.png");
	//c1->Print("New_Kin_pim_OldAcc.png");
	c1->Print("New_Kin_pim_OldAcc_NoPhy.png");
	//c1->Print("New_Kin_pim_NoAcc.png");
	
	TCanvas *c2=new TCanvas("c2","theta and P_{e^{-}} plots", 1000,800);
	c2->Divide(2,4);
	i=1;
	c2->cd(i++);
	ef_p_theta->Draw("COLZ");
	c2->cd(i++);
	el_p_theta->Draw("COLZ");
	c2->cd(i++);
	ef_p->Draw();
	c2->cd(i++);
	el_p->Draw();
	c2->cd(i++);
	ef_theta->Draw();
	c2->cd(i++);
	el_theta->Draw();
	c2->cd(i++);
	ef_phi->Draw();
	c2->cd(i++);
	el_phi->Draw();
	
	//c2->Print("New_Kin_ele_NewAcc.png");
	//c2->Print("New_Kin_ele_NewAcc_NoPhy.png");
	//c2->Print("New_Kin_ele_OldAcc.png");
	c2->Print("New_Kin_ele_OldAcc_NoPhy.png");
	//c2->Print("New_Kin_ele_NoAcc.png");

	/*}}}*/
}


