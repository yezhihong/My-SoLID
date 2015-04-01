#include<iostream>
#include<fstream>
using namespace std;

const double RADtoDEG=180./3.1415926;

void generate_plot(TString particle,int hadron_flag){
//particle: pion    kaon
//hadron_flag   1 for hp      2 for hm

	gStyle->SetOptStat(0);

	TChain *T=new TChain("T");
	TString outputfile="";
	if(particle=="pion"){
		if(hadron_flag==1){
			T->Add("./collider_3he_pip_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_pip_11_0_2_0.root");   // save Q2<10, pt>1
			outputfile="He3_pip_11_0_histo.root";
		}
		else if(hadron_flag==2){
			T->Add("./collider_3he_pim_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_pim_11_0_2_0.root");   // save Q2<10, pt>1
			outputfile="He3_pim_11_0_histo.root";
		}
	}else if(particle=="kaon"){
		if(hadron_flag==1){
			T->Add("./collider_3he_kp_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_kp_11_0_2_0.root");   // save Q2<10, pt>1
			outputfile="He3_kp_11_0_histo.root";
		}
		else if(hadron_flag==2){
			T->Add("./collider_3he_km_11_0_1_0.root");   // save Q2<10, pt<1
			T->Add("./collider_3he_km_11_0_2_0.root");   // save Q2<10, pt>1
			outputfile="He3_km_11_0_histo.root";
		}
	}
	//useful variables: Q2, W, Wp, x,y,z nu, s, pt, phi_h, phi_s weight_hp, weight_hm   
	//weight is xs in nbarn  weight_hp= dxs_hp*cos_theta_coverage*phi_coverage*energy_coverage/N_simulate

	double Q2,x,pt,W,Wp,z;   //kinematics
	double mom_ele,mom_had,theta_ele,theta_had;  // acceptance calculation input
	double weight_hp,weight_hm,event_weight;
	
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
	/*Define new historgrams{{{*/
	
	int Nbinx=100; //bin_size = (1-0)/200
	int Nbinz=100; //bin_size = (1-0)/200
	int Nbin_q2=100; //bin_size = (10-0)/200
	int Nbin_pt=100; //bin_size = (2.1-0)/200
	double Size_theta = 0.2;  int Nbin_theta=(int)((32.-5.)/Size_theta) + 1;//DEG
	double Size_p = 0.05;     int Nbin_p=(int)((6.5-0.)/Size_p) + 1;//GeV
    int Nbin_dx = 100;

	TFile *outf = new TFile(outputfile.Data(), "recreate");
	
	// __________z vs pt 
	TH2F *h_z_pt=new TH2F("h_z_pt","z VS P_{t} ",Nbin_pt,0,2.1,Nbinz,0,1);
	h_z_pt->GetXaxis()->SetTitle("P_{t}(GeV)");
	h_z_pt->GetYaxis()->SetTitle("z");
	
	//_________x vs pt 
	TH2F *h_x_pt=new TH2F("h_x_pt","x VS P_{t} ",Nbin_pt,0,2.1,Nbinx,0,1);
	h_x_pt->GetXaxis()->SetTitle("P_{T} (GeV)");
	h_x_pt->GetYaxis()->SetTitle("x_{bj}");
	
	//_______Q2 vs pt  
	TH2F *h_Q2_pt=new TH2F("h_Q2_pt","Q^{2} VS P_{t} ",Nbin_pt,0,2.1,Nbin_q2,0,10);
	h_Q2_pt->GetXaxis()->SetTitle("P_{t} (GeV)");
	h_Q2_pt->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	
	//______pt vs p
	TH2F *h_pt_p=new TH2F("h_pt_p","P_{t} VS P ",Nbin_p,0,6.5,Nbin_pt,0,2.1);
	h_pt_p->GetXaxis()->SetTitle("P (GeV)");
	h_pt_p->GetYaxis()->SetTitle("P_{t}(GeV)");
	
	//______x vs p
	TH2F *h_x_p=new TH2F("h_x_p","x VS P ",Nbin_p,0,6.5,Nbinx,0,1);
	h_x_p->GetXaxis()->SetTitle("P (GeV)");
	h_x_p->GetYaxis()->SetTitle("x_{bj}");

	//______Q2 vs p 
	TH2F *h_Q2_p=new TH2F("h_Q2_p","Q^{2} VS P ",Nbin_p,0,6.5,Nbin_q2,0,10);
	h_Q2_p->GetXaxis()->SetTitle("P (GeV)");
	h_Q2_p->GetYaxis()->SetTitle("Q^2(GeV^{2})");
	
	//_____p vs theta
	TH2F *h_p_theta=new TH2F("h_p_theta","P_{hadron} VS #theta_{hadron}",Nbin_theta,5-Size_theta/2.,32+Size_theta/2.,Nbin_p,0-Size_p/2.,6.5+Size_p/2.);
	h_p_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	h_p_theta->GetYaxis()->SetTitle("P_{hadron}");
	
	//_____TH1F of p
	TH1F *h_p=new TH1F("h_p","#p_{hadron}",Nbin_p,0-Size_p/2.,6.5+Size_p/2.);
	h_p->GetXaxis()->SetTitle("#p_{hadron}");
	h_p->GetYaxis()->SetTitle("rate");
	
	//_____TH1F of theta
	TH1F *h_theta=new TH1F("h_theta","#theta_{hadron}",Nbin_theta,5-Size_theta/2.,32+Size_theta/2.);
	h_theta->GetXaxis()->SetTitle("#theta_{hadron}");
	h_theta->GetYaxis()->SetTitle("rate");
    /*}}}*/	

	/*Fill histograms{{{*/
	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		if(1.0>0.0){//any additional cuts should be added in here
			if(hadron_flag==1){  //doing hp
				event_weight=weight_hp;
			}else if(hadron_flag==2){  //doing hm
				event_weight=weight_hm;
			}

			//fill histogram here
			//z VS pt
			h_z_pt->Fill(pt,z,event_weight);

			//x vs pt
			h_x_pt->Fill(pt,x,event_weight);

			//Q2 vs pt
			h_Q2_pt->Fill(pt,Q2,event_weight);
			
			//pt vs p
			h_pt_p->Fill(mom_had,pt,event_weight);
			
			//x vs p
			h_x_p->Fill(mom_had,x,event_weight);
			
			//Q2 vs p
			h_Q2_p->Fill(mom_had,Q2,event_weight);
			
			//p vs theta
			h_p_theta->Fill(theta_had*RADtoDEG,mom_had,event_weight);
			
			//TH1F theta
			h_p->Fill(mom_had,event_weight);

			//TH1F theta
			h_theta->Fill(theta_had*RADtoDEG,event_weight);
		}

	}// events loop ends here
    /*}}}*/

	/*Plot & Sav{{{*/
	/*
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
	c1->Divide(1,2);
	c1->cd(1);
	h_p_theta->Draw("COLZ");
	c1->cd(2);
	h_theta->Draw();
*/
	outf->cd();
	h_z_pt->Write();
	h_x_pt->Write();
	h_Q2_pt->Write();
	h_pt_p->Write();
	h_x_p->Write();
	h_Q2_p->Write();
	h_p_theta->Write();
	h_p->Write();
	h_theta->Write();
	outf->Close();
	/*}}}*/
}

