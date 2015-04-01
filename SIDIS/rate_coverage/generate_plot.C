#include<iostream>
#include<fstream>
using namespace std;


void generate_plot(TString particle){
//  change hadron_flag for doing hp and hm plots
//  change add file to the Chain to do pion and kaon

	int hadron_flag=1;  // 1 for hp   2 for hm

	TChain *T=new TChain("T");
	if(particle=="pion"){
		T->Add("./collider_3he_pip_11_0_1_0.root");   //Q2<10, pt<1
		T->Add("./collider_3he_pip_11_0_2_0.root");   //Q2<10, pt>1
	}else if(particle=="kaon"){
		T->Add("./collider_3he_kp_11_0_1_0.root");   //Q2<10, pt<1
		T->Add("./collider_3he_kp_11_0_2_0.root");   //Q2<10, pt>1
	}
	//useful variables: Q2, W, Wp, x,y,z nu, s, pt, phi_h, phi_s weight_hp, weight_hm   
	//weight is xs in nbarn  weight_hp= dxs_hp*cos_theta_coverage*phi_coverage*energy_coverage/N_simulate
	
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

	//read in acceptance histograms
	TFile *file_positive=new TFile("./acceptance/acceptance_solid_CLEO_SIDIS_3he_positive_output.root");
	TH2F *accep_positive_forward=(TH2F*)file_positive->Get("acceptance_forwardangle");
	//no need to care about positive large angle accep, since we don't care positron acceptance right now
	
	TFile *file_negative=new TFile("./acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TH2F *accep_negative_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_negative_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance

	int N_break =0;
	double rate_forward=0;
	double rate_large=0;
	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		
		int ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
		int ele_p_bin=int(mom_ele/0.1)+1;      //0.1 GeV per bin
		double ele_forward_acceptance=accep_negative_forward->GetBinContent(ele_theta_bin,ele_p_bin);
		double ele_large_acceptance=accep_negative_large->GetBinContent(ele_theta_bin,ele_p_bin);

		int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
		int hadron_p_bin=int(mom_had/0.1)+1;     //0.1 GeV per bin
		double hadron_acceptance=0;
		double event_weight=0;   //depend on hadron charge
		if(hadron_flag==1){  //doing hp
			hadron_acceptance=accep_positive_forward->GetBinContent(hadron_theta_bin,hadron_p_bin);
			event_weight=weight_hp;
		}else if(hadron_flag==2){  //doing hm
			hadron_acceptance=accep_negative_forward->GetBinContent(hadron_theta_bin,hadron_p_bin);
			event_weight=weight_hm;
		}
			
		double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
		double large_acceptance=ele_large_acceptance*hadron_acceptance;

		rate_forward+=event_weight*forward_acceptance;
		rate_large+=event_weight*large_acceptance;

	//	cerr<<Form("--- Electron: P=%f,B_P=%d, T=%f, B_T=%d, acc_for=%f, acc_lar=%f ", mom_ele, ele_p_bin, theta_ele, ele_theta_bin, ele_forward_acceptance, ele_large_acceptance)<<endl;
//		cerr<<Form("---   Hadron: P=%f,B_P=%d, T=%f, B_T=%d, acc_for=%f,R_F=%f, R_L=%f", mom_had, hadron_p_bin, theta_had, hadron_theta_bin,hadron_acceptance,rate_forward,rate_large)<<endl;

	//}
	}// events loop ends here

	const double Lumi = 3.0e36; // cm-2*s-1
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;


	cout<<"rate_forward: "<<rate_forward*Lumi*KHz*nBcm2<<endl;
	cout<<"rate_large: "<<rate_large*Lumi*KHz*nBcm2<<endl;
}

