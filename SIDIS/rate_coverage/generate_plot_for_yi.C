#include<iostream>
#include<fstream>
using namespace std;

/*
//_________________acceptance file__________
acceptance_solid_CLEO_SIDIS_3he_negative_output.root
acceptance_solid_CLEO_SIDIS_3he_negative_DIRC_2side_output.root
acceptance_solid_CLEO_SIDIS_3he_negative_DIRC_3side_output.root
acceptance_solid_CLEO_SIDIS_3he_negative_DIRC_4side_output.root
*/


void generate_plot_for_yi(TString particle,int hadron_flag,int DIRC_num){
//particle: pion    kaon
//hadron_flag   1 for hp      2 for hm
//DIRC_num 4: (2 sides)  6: (3 sides)   8: (4 sides)   no meaning for pion case 

	gStyle->SetOptStat(0);

	TChain *T=new TChain("T");
	if(particle=="pion"&&hadron_flag==1){
		T->Add("./collider_3h3_pip_11_0_1_1.root");   // save Q2<10, pt<1
		T->Add("./collider_3he_pip_11_0_2_1.root");   // save Q2<10, pt>1
	}else if(particle=="pion"&&hadron_flag==2){
		T->Add("./collider_3h3_pim_11_0_1_1.root");   // save Q2<10, pt<1
		T->Add("./collider_3he_pim_11_0_2_1.root");   // save Q2<10, pt>1
	}else if(particle=="kaon"){
			T->Add("../collider_3he_kp_11_0_1_1.root");   // save Q2<10, pt<1
		T->Add("../collider_3he_kp_11_0_2_1.root");   // save Q2<10, pt>1
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
	T->SetBranchAddress("weight_hp",&weight_hp);   //unit in nbar
	T->SetBranchAddress("weight_hm",&weight_hm);

	Long64_t N_entries=T->GetEntries();
	cout<<"total generated events number: "<<N_entries<<endl;
	
	//deal with acceptance files 
	//read in acceptance histograms
	//electron acceptance	
	TFile *file_negative=new TFile("./acceptanc/eacceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TFile *file_hadron=0;
	if(particle=="kaon"){
		if(DIRC_num==4){
			file_hadron=new TFile("acceptance_solid_CLEO_SIDIS_3he_negative_DIRC_2side_output.root");
		}else if(DIRC_num==6){
			file_hadron=new TFile("acceptance_solid_CLEO_SIDIS_3he_negative_DIRC_3side_output.root");
		}else if(DIRC_num==8){
			file_hadron=new TFile("acceptance_solid_CLEO_SIDIS_3he_negative_DIRC_4side_output.root");
		}else{
			cout<<"invalid DIRC_num, don't have its acceptance file!!!!!!!"<<endl;
			cout<<"..............now using pion acceptance!!!!!!"<<endl;
			file_hadron=file_negative;
		}
	}else if(particle=="pion"){
		cout<<"...now doing pion calculation,DIRC_num has no meaning, setting it to zero."<<endl;
		DIRC_num=0;
	}else{
		cout<<"not the right input for particle, input is either pion or kaon"<<endl;
	}
	
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
	TH2F *accep_had=0;
	if(particle=="pion"){
		//pion is normal acceptance histogram, the same as electron acceptance
		accep_had=(TH2F*)accep_ele_forward->Clone("acceptance_hadron");
	}else{
		accep_had=(TH2F*)file_hadron->Get("acceptance_forwardangle");
	}
		


	int Nbinx=100;
	int Nbiny=100;
	int Nbin_theta=100;
	int Nbin_p=100;
	TFile *savefile=new TFile("generate_plot_for_yi_output.root","recreate");
	// __________z vs pt 
	//forward angle
	TH2F *hf_z_pt=new TH2F("hf_z_pt","z VS P_{t} ",Nbinx,0,2,Nbiny,0,1);
	hf_z_pt->GetYaxis()->SetTitle("z");
	hf_z_pt->GetXaxis()->SetTitle("P_{t}(GeV)");
	//large angle
	TH2F *hl_z_pt=new TH2F("hl_z_pt","z VS P_{t} ",Nbinx,0,2,Nbiny,0,1);
	hl_z_pt->GetYaxis()->SetTitle("z");
	hl_z_pt->GetXaxis()->SetTitle("P_{t}(GeV)");
	//total
	TH2F *h_z_pt=new TH2F("h_z_pt","z VS P_{t} ",Nbinx,0,2,Nbiny,0,1);
	h_z_pt->GetYaxis()->SetTitle("z");
	h_z_pt->GetXaxis()->SetTitle("P_{t}(GeV)");
	
	//_________x vs pt 
	//forward angle
	TH2F *hf_x_pt=new TH2F("hf_x_pt","x VS P_{t} ",Nbinx,0,2,Nbiny,0,1);
	hf_x_pt->GetYaxis()->SetTitle("x_{bj}");
	hf_x_pt->GetXaxis()->SetTitle("P_{T} (GeV)");
	//large angle
	TH2F *hl_x_pt=new TH2F("hl_x_pt","x VS P_{t} ",Nbinx,0,2,Nbiny,0,1);
	hl_x_pt->GetYaxis()->SetTitle("x_{bj}");
	hl_x_pt->GetXaxis()->SetTitle("P_{T} (GeV)");
	//total
	TH2F *h_x_pt=new TH2F("h_x_pt","x VS P_{t} ",Nbinx,0,2,Nbiny,0,1);
	h_x_pt->GetYaxis()->SetTitle("x_{bj}");
	h_x_pt->GetXaxis()->SetTitle("P_{T} (GeV)");
	
	//_______Q2 vs pt  
	//forward angle
	TH2F *hf_Q2_pt=new TH2F("hf_Q2_pt","Q^{2} VS P_{t} ",Nbinx,0,2,Nbiny,0,10);
	hf_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hf_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");
	//large angle
	TH2F *hl_Q2_pt=new TH2F("hl_Q2_pt","Q^{2} VS P_{t} ",Nbinx,0,2,Nbiny,0,10);
	hl_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hl_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");
	//total
	TH2F *h_Q2_pt=new TH2F("h_Q2_pt","Q^{2} VS P_{t} ",Nbinx,0,2,Nbiny,0,10);
	h_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	h_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");
	
	//______pt vs p
	//forward angle
	TH2F *hf_pt_p=new TH2F("hf_pt_p","P_{t} VS P ",Nbinx,0,7,Nbiny,0,2);
	hf_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");
	hf_pt_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_pt_p=new TH2F("hl_pt_p","P_{t} VS P ",Nbinx,0,7,Nbiny,0,2);
	hl_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");
	hl_pt_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_pt_p=new TH2F("h_pt_p","P_{t} VS P ",Nbinx,0,7,Nbiny,0,2);
	h_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");
	h_pt_p->GetXaxis()->SetTitle("P (GeV)");
	
	//______x vs p
	//forward angle
	TH2F *hf_x_p=new TH2F("hf_x_p","x VS P ",Nbinx,0,7,Nbiny,0,1);
	hf_x_p->GetYaxis()->SetTitle("x_{bj}");
	hf_x_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_x_p=new TH2F("hl_x_p","x VS P ",Nbinx,0,7,Nbiny,0,1);
	hl_x_p->GetYaxis()->SetTitle("x_{bj}");
	hl_x_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_x_p=new TH2F("h_x_p","x VS P ",Nbinx,0,7,Nbiny,0,1);
	h_x_p->GetYaxis()->SetTitle("x_{bj}");
	h_x_p->GetXaxis()->SetTitle("P (GeV)");
	

	//______Q2 vs p 
	//forward angle
	TH2F *hf_Q2_p=new TH2F("hf_Q2_p","Q^{2} VS P ",Nbinx,0,7,Nbiny,0,10);
	hf_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hf_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	//large angle
	TH2F *hl_Q2_p=new TH2F("hl_Q2_p","Q^{2} VS P ",Nbinx,0,7,Nbiny,0,10);
	hl_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	hl_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	//total
	TH2F *h_Q2_p=new TH2F("h_Q2_p","Q^{2} VS P ",Nbinx,0,7,Nbiny,0,10);
	h_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	h_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	
	//_____p vs theta
	TH2F *hf_p_theta=new TH2F("hf_p_theta","P_{hadron} VS #theta_{hadron}  (forward angle)",Nbin_theta,7.0/180.0*3.14159265,30.0/180.0*3.14159265,Nbin_p,0,10);
	//forward angle
	hf_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	hf_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	//large angle
	TH2F *hl_p_theta=new TH2F("hl_p_theta","P_{hadron} VS #theta_{hadron} (large angle)",Nbin_theta,7.0/180.0*3.14159265,30.0/180.0*3.14159265,Nbin_p,0,10);
	hl_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	hl_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	//total
	TH2F *h_p_theta=new TH2F("h_p_theta","P_{hadron} VS #theta_{hadron} (large angle)",Nbin_theta,7.0/180.0*3.14159265,30.0/180.0*3.14159265,Nbin_p,0,10);
	h_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	h_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	
	//_____TH1F of theta
	//forward angle
	TH1F *hf_theta=new TH1F("hf_theta","rate as function of #theta_{hadron}  (forward angle)",Nbin_theta,7.0/180.0*3.14159265,30.0/180.0*3.14159265);
	hf_theta->GetYaxis()->SetTitle("rate");
	hf_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	//large angle
	TH1F *hl_theta=new TH1F("hl_theta","rate as function of #theta_{hadron} (large angle)",Nbin_theta,7.0/180.0*3.14159265,30.0/180.0*3.14159265);
	hl_theta->GetYaxis()->SetTitle("rate");
	hl_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	//total
	TH1F *h_theta=new TH1F("h_theta","rate as function of #theta_{hadron} (large angle)",Nbin_theta,7.0/180.0*3.14159265,30.0/180.0*3.14159265);
	h_theta->GetYaxis()->SetTitle("rate");
	h_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	

	double rate_forward=0;
	double rate_large=0;
	
	double rate_forward_p_g_4[4]={0};
	double rate_large_p_g_4[4]={0};
	double rate_forward_p_l_4=0;
	double rate_large_p_l_4=0;

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		//if(1.0>0.0){//any additional cuts should be added in here
		if(mom_ele>1.0 && mom_ele<10.0 && W>=2.3 && Wp>=1.6 && z>0.3 && z<0.7){//any additional cuts should be added in here
			int ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
			int ele_p_bin=int(mom_ele/0.05)+1;      //0.05 GeV per bin for mom
			double ele_forward_acceptance=accep_ele_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			double ele_large_acceptance=accep_ele_large->GetBinContent(ele_theta_bin,ele_p_bin);

			// mandatary cut 8 degree out
			if(theta_ele/3.14159265*180.0 < 8.0){
				ele_forward_acceptance=0;
				ele_large_acceptance=0;
			}

			int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
			int hadron_p_bin=int(mom_had/0.05)+1;     //0.05 GeV per bin for mom
			double hadron_acceptance=accep_had->GetBinContent(hadron_theta_bin,hadron_p_bin);
			double event_weight=0;   //depend on hadron charge

			//mandatary cut 8 degree out
			if(theta_had/3.14159265*180.0 < 8.0){
				hadron_acceptance=0;
			}
			
			if(hadron_flag==1){  //doing hp
				event_weight=weight_hp;
			}else if(hadron_flag==2){  //doing hm
				event_weight=weight_hm;
			}
			
			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

			//fill histogram here
			//z VS pt
			hf_z_pt->Fill(pt,z,event_weight*forward_acceptance);
			hl_z_pt->Fill(pt,z,event_weight*large_acceptance);
			h_z_pt->Fill(pt,z,event_weight*(forward_acceptance+large_acceptance));

			//x vs pt
			hf_x_pt->Fill(pt,x,event_weight*forward_acceptance);
			hl_x_pt->Fill(pt,x,event_weight*large_acceptance);
			h_x_pt->Fill(pt,x,event_weight*(forward_acceptance+large_acceptance));

			//Q2 vs pt
			hf_Q2_pt->Fill(pt,Q2,event_weight*forward_acceptance);
			hl_Q2_pt->Fill(pt,Q2,event_weight*large_acceptance);
			h_Q2_pt->Fill(pt,Q2,event_weight*(forward_acceptance+large_acceptance));
			
			//pt vs p
			hf_pt_p->Fill(mom_had,pt,event_weight*forward_acceptance);
			hl_pt_p->Fill(mom_had,pt,event_weight*large_acceptance);
			h_pt_p->Fill(mom_had,pt,event_weight*(forward_acceptance+large_acceptance));
			
			//x vs p
			hf_x_p->Fill(mom_had,x,event_weight*forward_acceptance);
			hl_x_p->Fill(mom_had,x,event_weight*large_acceptance);
			h_x_p->Fill(mom_had,x,event_weight*(forward_acceptance+large_acceptance));
			
			//Q2 vs p
			hf_Q2_p->Fill(mom_had,Q2,event_weight*forward_acceptance);
			hl_Q2_p->Fill(mom_had,Q2,event_weight*large_acceptance);
			h_Q2_p->Fill(mom_had,Q2,event_weight*(forward_acceptance+large_acceptance));
			
			//p vs theta
			hf_p_theta->Fill(theta_had,mom_had,event_weight*forward_acceptance);
			hl_p_theta->Fill(theta_had,mom_had,event_weight*large_acceptance);
			h_p_theta->Fill(theta_had,mom_had,event_weight*(forward_acceptance+large_acceptance));
			
			//TH1F theta
			hf_theta->Fill(theta_had,event_weight*forward_acceptance);
			hl_theta->Fill(theta_had,event_weight*large_acceptance);
			h_theta->Fill(theta_had,event_weight*(forward_acceptance+large_acceptance));
			
			rate_forward+=event_weight*forward_acceptance;
			rate_large+=event_weight*large_acceptance;
	
			if(mom_had>=4&&mom_had<4.5){
				rate_forward_p_g_4[0]+=event_weight*forward_acceptance;
				rate_large_p_g_4[0]+=event_weight*large_acceptance;
			}else if(mom_had>=4.5&&mom_had<5){
				rate_forward_p_g_4[1]+=event_weight*forward_acceptance;
				rate_large_p_g_4[1]+=event_weight*large_acceptance;
			}else if(mom_had>=5&&mom_had<5.5){
				rate_forward_p_g_4[2]+=event_weight*forward_acceptance;
				rate_large_p_g_4[2]+=event_weight*large_acceptance;
			}else if(mom_had>=5.5&&mom_had<6){
				rate_forward_p_g_4[3]+=event_weight*forward_acceptance;
				rate_large_p_g_4[3]+=event_weight*large_acceptance;
			}else{
				rate_forward_p_l_4+=event_weight*forward_acceptance;
				rate_large_p_l_4+=event_weight*large_acceptance;
			}
		}

	}// events loop ends here
	cout<<"______forward and large angle integral rate____________________"<<endl;
	cout<<"rate_forward: "<<rate_forward<<endl;
	cout<<"rate_large: "<<rate_large<<endl;
	
	cout<<"______________ p~[4,4.5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[0]<<endl;
	cout<<"large:   "<<rate_large_p_g_4[0]<<endl;
	
	cout<<"______________ p~[4.5,5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[1]<<endl;
	cout<<"large:   "<<rate_large_p_g_4[1]<<endl;
	
	cout<<"______________ p~[5,5.5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[2]<<endl;
	cout<<"large:   "<<rate_large_p_g_4[2]<<endl;
	
	cout<<"______________ p~[5,5.5) GeV rate _________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_g_4[3]<<endl;
	cout<<"large:   "<<rate_large_p_g_4[3]<<endl;
	
	cout<<"______________ p<4 GeV rate __________________________________"<<endl;
	cout<<"forward: "<<rate_forward_p_l_4<<endl;
	cout<<"large:   "<<rate_large_p_l_4<<endl;
	cout<<"_______________________________________________________________"<<endl;

/*
	hl_z_pt->SetMarkerColor(3);
	hl_x_pt->SetMarkerColor(3);
	hl_Q2_pt->SetMarkerColor(3);
	hl_pt_p->SetMarkerColor(3);
	hl_x_p->SetMarkerColor(3);
	hl_Q2_p->SetMarkerColor(3);
*/

/*
	hf_z_pt->SetDirectory(savefile);
	hf_x_pt->SetDirectory(savefile);
	hf_Q2_pt->SetDirectory(savefile);
	hf_pt_p->SetDirectory(savefile);
	hf_x_p->SetDirectory(savefile);
	hf_Q2_p->SetDirectorty(savefile);
	
	hl_z_pt->SetDirectory(savefile);
	hl_x_pt->SetDirectory(savefile);
	hl_Q2_pt->SetDirectory(savefile);
	hl_pt_p->SetDirectory(savefile);
	hl_x_p->SetDirectory(savefile);
	hl_Q2_p->SetDirectorty(savefile);
	
	savefile->Write();
*/


	TCanvas *can_large=new TCanvas("can_large","large angle");
	can_large->Divide(2,3);
	can_large->cd(1);
	hl_z_pt->Draw("COLZ");
	can_large->cd(2);
	hl_pt_p->Draw("colz");
	can_large->cd(3);
	hl_x_pt->Draw("colz");
	can_large->cd(4);
	hl_x_p->Draw("colz");
	can_large->cd(5);
	hl_Q2_pt->Draw("colz");
	can_large->cd(6);
	hl_Q2_p->Draw("colz");
//	can_large->cd(7);
//	hl_mom_theta->Draw("colz");
	
	TCanvas *can_forward=new TCanvas("can_forward","forward angle");
	can_forward->Divide(2,3);
	can_forward->cd(1);
	hf_z_pt->Draw("colz");
	can_forward->cd(2);
	hf_pt_p->Draw("colz");
	can_forward->cd(3);
	hf_x_pt->Draw("colz");
	can_forward->cd(4);
	hf_x_p->Draw("colz");
	can_forward->cd(5);
	hf_Q2_pt->Draw("colz");
	can_forward->cd(6);
	hf_Q2_p->Draw("colz");
//	can_forward->cd(7);
//	hf_mom_theta->Draw("colz");

	TCanvas *can_total=new TCanvas("can_total","total");
	can_total->Divide(2,3);
	can_total->cd(1);
	h_z_pt->Draw("colz");
	can_total->cd(2);
	h_pt_p->Draw("colz");
	can_total->cd(3);
	h_x_pt->Draw("colz");
	can_total->cd(4);
	h_x_p->Draw("colz");
	can_total->cd(5);
	h_Q2_pt->Draw("colz");
	can_total->cd(6);
	h_Q2_p->Draw("colz");

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
}

