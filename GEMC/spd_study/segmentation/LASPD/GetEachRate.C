/*Header{{{*/
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TF1.h>
#include <TLatex.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>
/*}}}*/

using namespace std;

void GetEachRate(	TString input_filename){
	gStyle->SetOptStat(0);

	TChain *T = new TChain("T");
	//T->Add(input_filename.Data());
	if(input_filename=="eDIS"||input_filename=="pi0"||input_filename=="pip"||input_filename=="pim"||input_filename=="p"){
		T->Add(Form("LASPD_background_SIDIS_He3_%s.root",input_filename.Data()));
		T->Add(Form("LASPD_background_SIDIS_He3_%s_up.root",input_filename.Data()));
		T->Add(Form("LASPD_background_SIDIS_He3_%s_down.root",input_filename.Data()));
	}
	else if(input_filename=="background"){
		T->Add("LASPD_background_SIDIS_He3_pi0.root");
		T->Add("LASPD_background_SIDIS_He3_pi0_up.root");
		T->Add("LASPD_background_SIDIS_He3_pi0_down.root");
		T->Add("LASPD_background_SIDIS_He3_pip.root");
		T->Add("LASPD_background_SIDIS_He3_pip_up.root");
		T->Add("LASPD_background_SIDIS_He3_pip_down.root");
		T->Add("LASPD_background_SIDIS_He3_pim.root");
		T->Add("LASPD_background_SIDIS_He3_pim_up.root");
		T->Add("LASPD_background_SIDIS_He3_pim_down.root");
		T->Add("LASPD_background_SIDIS_He3_p.root");
		T->Add("LASPD_background_SIDIS_He3_p_up.root");
		T->Add("LASPD_background_SIDIS_He3_p_down.root");
		T->Add("LASPD_background_SIDIS_He3_EM.root");
	}
	else if(input_filename=="all"){
		T->Add("LASPD_background_SIDIS_He3_pi0.root");
		T->Add("LASPD_background_SIDIS_He3_pi0_up.root");
		T->Add("LASPD_background_SIDIS_He3_pi0_down.root");
		T->Add("LASPD_background_SIDIS_He3_pip.root");
		T->Add("LASPD_background_SIDIS_He3_pip_up.root");
		T->Add("LASPD_background_SIDIS_He3_pip_down.root");
		T->Add("LASPD_background_SIDIS_He3_pim.root");
		T->Add("LASPD_background_SIDIS_He3_pim_up.root");
		T->Add("LASPD_background_SIDIS_He3_pim_down.root");
		T->Add("LASPD_background_SIDIS_He3_p.root");
		T->Add("LASPD_background_SIDIS_He3_p_up.root");
		T->Add("LASPD_background_SIDIS_He3_p_down.root");
		T->Add("LASPD_background_SIDIS_He3_eDIS.root");
		T->Add("LASPD_background_SIDIS_He3_eDIS_up.root");
		T->Add("LASPD_background_SIDIS_He3_eDIS_down.root");
		T->Add("LASPD_background_SIDIS_He3_EM.root");
	}
	else
		T->Add(Form("LASPD_background_SIDIS_He3_%s.root",input_filename.Data()));

	/*Tree and Histogram{{{*/
	double rate,mom_flux,theta_flux,phi_flux;
	double r_flux,x_flux,y_flux,z_flux,px_flux,py_flux,pz_flux,vx_flux,vy_flux,vz_flux;
	double R_EC, EDep, E, convert, EC_Cut;
	int evn,evn_flux,nhit_flux,ID_flux, PID_flux,ID_Pick, Module,Slide;

	T->SetBranchAddress("mom_flux", &mom_flux);
	T->SetBranchAddress("theta_flux", &theta_flux);
	T->SetBranchAddress("phi_flux",&phi_flux);
	T->SetBranchAddress("r_flux",&r_flux);
	T->SetBranchAddress("x_flux",&x_flux);
	T->SetBranchAddress("y_flux",&y_flux);
	T->SetBranchAddress("z_flux",&z_flux);
	T->SetBranchAddress("px_flux",&px_flux);
	T->SetBranchAddress("py_flux",&py_flux);
	T->SetBranchAddress("pz_flux",&pz_flux);
	T->SetBranchAddress("vx_flux",&vx_flux);
	T->SetBranchAddress("vy_flux",&vy_flux);
	T->SetBranchAddress("vz_flux",&vz_flux);
	T->SetBranchAddress("R_EC",&R_EC);
	T->SetBranchAddress("EC_Cut",&EC_Cut);
	T->SetBranchAddress("EDep",&EDep);
	T->SetBranchAddress("E",&E);
	T->SetBranchAddress("rate", &rate);
	T->SetBranchAddress("convert",&convert);
	T->SetBranchAddress("ID_flux",&ID_flux);
	T->SetBranchAddress("PID_flux",&PID_flux);
	T->SetBranchAddress("ID_Pick",&ID_Pick);
	T->SetBranchAddress("Module",&Module);
	T->SetBranchAddress("Slide",&Slide);
	T->SetBranchAddress("evn",  &evn);
	T->SetBranchAddress("evn_flux", &evn_flux);
	T->SetBranchAddress("nhit_flux", &nhit_flux);

	TH2F *hxy_e1 = new TH2F("hxy_e1","e^{#pm} (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 100,-145,145,100,-145,145);
	hxy_e1->SetXTitle("x_{flux}");	hxy_e1->SetYTitle("y_{flux}");
	hxy_e1->GetXaxis()->SetTitleSize(0.06); hxy_e1->GetYaxis()->SetTitleSize(0.06);
	hxy_e1->GetXaxis()->CenterTitle(1);	hxy_e1->GetYaxis()->CenterTitle(1);
	TH2F *hxy_e2 = new TH2F("hxy_e2","e^{#pm} (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD (EC triggered)", 100,-145,145,100,-145,145);
	hxy_e2->SetXTitle("x_{flux}");	hxy_e2->SetYTitle("y_{flux}");
	hxy_e2->GetXaxis()->SetTitleSize(0.06); hxy_e2->GetYaxis()->SetTitleSize(0.06);
	hxy_e2->GetXaxis()->CenterTitle(1);	hxy_e2->GetYaxis()->CenterTitle(1);
	TH2F *hxy_g1 = new TH2F("hxy_g1","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 100,-145,145,100,-145,145);
	hxy_g1->SetXTitle("x_{flux}");	hxy_g1->SetYTitle("y_{flux}");
	hxy_g1->GetXaxis()->CenterTitle(1);	hxy_g1->GetYaxis()->CenterTitle(1);
	hxy_g1->GetXaxis()->SetTitleSize(0.06); hxy_g1->GetYaxis()->SetTitleSize(0.06);
	TH2F *hxy_g2 = new TH2F("hxy_g2","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD (EC triggered)", 100,-145,145,100,-145,145);
   	hxy_g2->SetXTitle("x_{flux}");	hxy_g2->SetYTitle("y_{flux}");
	hxy_g2->GetXaxis()->SetTitleSize(0.06); hxy_g2->GetYaxis()->SetTitleSize(0.06);
	hxy_g2->GetXaxis()->CenterTitle(1);	hxy_g2->GetYaxis()->CenterTitle(1);

	TH1F *hr_e1 = new TH1F("hr_e1","e (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD", 70,75,145);
	hr_e1->SetXTitle("r_{flux}");
	hr_e1->GetXaxis()->CenterTitle(1);	
	hr_e1->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr_e2 = new TH1F("hr_e2","e (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD (EC triggered)", 70,75,145);
	hr_e2->SetXTitle("r_{flux}");
	hr_e2->GetXaxis()->CenterTitle(1);	
	hr_e2->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr_g1 = new TH1F("hr_g1","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD", 70,75,145);
	hr_g1->SetXTitle("r_{flux}");
	hr_g1->GetXaxis()->CenterTitle(1);	
	hr_g1->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr_g2 = new TH1F("hr_g2","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD (EC triggered)", 70,75,145);
	hr_g2->SetXTitle("r_{flux}");
	hr_g2->GetXaxis()->CenterTitle(1);	
	hr_g2->GetXaxis()->SetTitleSize(0.06);

	TH2F *hxy1 = new TH2F("hxy1","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 100,-145,145,100,-145,145);
   	hxy1->SetXTitle("x_{flux}");	      hxy1->SetYTitle("y_{flux}");
	hxy1->GetXaxis()->SetTitleSize(0.06); hxy1->GetYaxis()->SetTitleSize(0.06);
	hxy1->GetXaxis()->CenterTitle(1);	  hxy1->GetYaxis()->CenterTitle(1);

	TH2F *hxy2 = new TH2F("hxy2","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD (EC triggered)", 100,-145,145,100,-145,145);
   	hxy2->SetXTitle("x_{flux}");	      hxy2->SetYTitle("y_{flux}");
	hxy2->GetXaxis()->SetTitleSize(0.06); hxy2->GetYaxis()->SetTitleSize(0.06);
	hxy2->GetXaxis()->CenterTitle(1);	  hxy2->GetYaxis()->CenterTitle(1);

	TH1F *hr1 = new TH1F("hr1","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD", 70,75,145);
	hr1->SetXTitle("r_{flux}");
	hr1->GetXaxis()->CenterTitle(1);	
	hr1->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr2 = new TH1F("hr2","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD (EC triggered)", 70,75,145);
	hr2->SetXTitle("r_{flux}");
	hr2->GetXaxis()->CenterTitle(1);	
	hr2->GetXaxis()->SetTitleSize(0.06);

	TH2F *hEDep_bgr = new TH2F("hEDep_bgr","EDep from e^{#pm}+e^{#pm}(#gamma) (#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 70,75,145,100,0,10);
   	hEDep_bgr->SetXTitle("r_{flux}");	      hEDep_bgr->SetYTitle("E_{Dep} (MeV)");
	hEDep_bgr->GetXaxis()->SetTitleSize(0.06); hEDep_bgr->GetYaxis()->SetTitleSize(0.06);
	hEDep_bgr->GetXaxis()->CenterTitle(1);	  hEDep_bgr->GetYaxis()->CenterTitle(1);

	/*}}}*/
	
	int Nevt = T->GetEntries();
    double rate_e1=0, rate_e2=0, rate_g1=0,rate_g2=0;
    double EDep_e1=0, EDep_e2=0, EDep_g1=0,EDep_g2=0;
	double rate_ep=0, rate_em=0, rate_pip=0,rate_pim=0,rate_p=0,rate_g=0;
	for(int i=0;i<Nevt;i++){
		T->GetEntry(i);
		rate /= 1000.0;
		//All Charged Particles
		if(abs(PID_flux)==11||abs(PID_flux)==111||abs(PID_flux)==2212){
			hxy_e1->Fill(x_flux,y_flux,rate*convert);
			hr_e1->Fill(r_flux,rate*convert);
			rate_e1 += rate*convert;
			EDep_e1 += EDep*1000.;
			//if(R_EC>80&&R_EC<=140){
				if(1){
				hxy_e2->Fill(x_flux,y_flux,rate*convert*EC_Cut);
				hr_e2->Fill(r_flux,rate*convert*EC_Cut);
				rate_e2 += rate*convert*EC_Cut;
				EDep_e2 += EDep*1000.*convert*EC_Cut;
			}
		}
		if(PID_flux==22){
				rate *= 2.0;//One photon becomes two electrons

				hxy_g1->Fill(x_flux,y_flux,rate*convert);
				hr_g1->Fill(r_flux,rate*convert);
				rate_g1 += rate*convert;
				EDep_g1 += EDep*1000.;
				//if(R_EC>80&&R_EC<=140){
				if(1){
					hxy_g2->Fill(x_flux,y_flux,rate*convert*EC_Cut);
					hr_g2->Fill(r_flux,rate*convert*EC_Cut);
					rate_g2 += rate*convert*EC_Cut;
					EDep_g2 += EDep*1000.*convert*EC_Cut;
				}
		}
			
		if(abs(PID_flux)==11||abs(PID_flux)==111||abs(PID_flux)==2212||PID_flux==22){
			if(PID_flux==22)
				rate *= 2.0;

			hxy1->Fill(x_flux,y_flux,rate*convert);
			hr1->Fill(r_flux,rate*convert);
			if(PID_flux==11)
				rate_em+= rate*convert;
			else if(PID_flux==-11)
				rate_ep+= rate*convert;
			else if(PID_flux==111)
				rate_pip+= rate*convert;
			else if(PID_flux==-111)
				rate_pim+= rate*convert;
			else if(PID_flux==2212)
				rate_p+= rate*convert;
			else if(PID_flux==22)
				rate_g+= rate*convert;

			if(R_EC>105&&R_EC<=140){
				hxy2->Fill(x_flux,y_flux,rate*convert*EC_Cut);
				hr2->Fill(r_flux,rate*convert*EC_Cut);
			}
			hEDep_bgr->Fill(r_flux, EDep*1000.,rate*convert);
		}
	}

	TLatex* t1 = new TLatex(); t1->SetNDC();t1->SetTextColor(4);

	TCanvas *c1 = new TCanvas("c1","c1",1000,1000);
	c1->Divide(2,2);
	c1->cd(1); hxy_e1->Draw("colz");	
	c1->cd(2); hr_e1->Draw("");	
	hr_e1->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e1));
	c1->cd(3); hxy_e2->Draw("colz");	
	c1->cd(4); hr_e2->Draw("");	
	hr_e2->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e2));
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e2));

	TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
	c2->Divide(2,2);
	c2->cd(1); hxy_g1->Draw("colz");	
	c2->cd(2); hr_g1->Draw("");	
	hr_g1->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_g1));
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_g1));
	c2->cd(3); hxy_g2->Draw("colz");	
	c2->cd(4); hr_g2->Draw("");	
	hr_g2->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_g2));
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_g2));

	TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
	c3->Divide(2,2);
	c3->cd(1); hxy1->Draw("colz");	
	c3->cd(2); hr1->Draw("");	
	hr1->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e1+rate_g1));
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e1+rate_g1));

	t1->DrawLatex(0.25,0.45,Form("#scale[0.9]{e^{#pm}(#pi^{0} #gamma) = %7.3e (KHz)}", rate_g));
	t1->DrawLatex(0.25,0.40,Form("#scale[0.9]{e^{+} = %7.3e (KHz)}", rate_ep));
	t1->DrawLatex(0.25,0.35,Form("#scale[0.9]{e^{-} = %7.3e (KHz)}", rate_em));
	t1->DrawLatex(0.25,0.30,Form("#scale[0.9]{#pi^{+} = %7.3e (KHz)}", rate_pip));
	t1->DrawLatex(0.25,0.25,Form("#scale[0.9]{#pi^{-} = %7.3e (KHz)}", rate_pim));
	t1->DrawLatex(0.25,0.20,Form("#scale[0.9]{p = %7.3e (KHz)}", rate_p));
	
	c3->cd(3); hxy2->Draw("colz");	
	c3->cd(4); hr2->Draw("");	
	hr2->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e2+rate_g2));
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e2+rate_g2));

	TCanvas *c4 = new TCanvas("c4","c4",1000,600);
	c4->Divide(2,1);
	c4->cd(1);
    hEDep_bgr->Draw("colz");
    c4->cd(2);
	TH1F* hEDep_bgr_px = (TH1F*) hEDep_bgr->ProfileX("hEDep_bgr_px");
    hEDep_bgr_px->Draw();

	TCanvas *c5 = new TCanvas("c5","c5",1000,500);
	c5->Divide(2,1);
	c5->cd(1); hxy1->Draw("colz");	
	c5->cd(2); hr1->Draw("");	
    hr1->SetTitle(Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e1+rate_g1));
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %7.3e (KHz)}", rate_e1+rate_g1));
    c5->Print(Form("hist_LASPD_%s.png",input_filename.Data()));

   TFile *new_file = new TFile(Form("histo_%s.root",input_filename.Data()),"recreate");
   new_file->cd();
   hr_e1->Write();hr_e2->Write();hr_g1->Write();hr_g2->Write();
   hxy_e1->Write();hxy_e2->Write();hxy_g1->Write();hxy_g2->Write();
   hxy1->Write(); hxy2->Write(); hr1->Write(); hr2->Write();
   hEDep_bgr->Write();  hEDep_bgr_px->Write();
   new_file->Close();

}

