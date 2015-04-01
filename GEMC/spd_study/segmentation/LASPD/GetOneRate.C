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

void GetOneRate(	TString input_filename){
	gStyle->SetOptStat(0);

	TChain *T = new TChain("T");
	T->Add(input_filename.Data());

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

	TH2F *hxy_e1 = new TH2F("hxy_e1","e^{#pm} (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 100,-150,150,100,-150,150);
	hxy_e1->SetXTitle("x_{flux}");	hxy_e1->SetYTitle("y_{flux}");
	hxy_e1->GetXaxis()->SetTitleSize(0.06); hxy_e1->GetYaxis()->SetTitleSize(0.06);
	hxy_e1->GetXaxis()->CenterTitle(1);	hxy_e1->GetYaxis()->CenterTitle(1);
	TH2F *hxy_e2 = new TH2F("hxy_e2","e^{#pm} (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD (EC triggered)", 100,-150,150,100,-150,150);
	hxy_e2->SetXTitle("x_{flux}");	hxy_e2->SetYTitle("y_{flux}");
	hxy_e2->GetXaxis()->SetTitleSize(0.06); hxy_e2->GetYaxis()->SetTitleSize(0.06);
	hxy_e2->GetXaxis()->CenterTitle(1);	hxy_e2->GetYaxis()->CenterTitle(1);
	TH2F *hxy_g1 = new TH2F("hxy_g1","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 100,-150,150,100,-150,150);
	hxy_g1->SetXTitle("x_{flux}");	hxy_g1->SetYTitle("y_{flux}");
	hxy_g1->GetXaxis()->CenterTitle(1);	hxy_g1->GetYaxis()->CenterTitle(1);
	hxy_g1->GetXaxis()->SetTitleSize(0.06); hxy_g1->GetYaxis()->SetTitleSize(0.06);
	TH2F *hxy_g2 = new TH2F("hxy_g2","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD (EC triggered)", 100,-150,150,100,-150,150);
   	hxy_g2->SetXTitle("x_{flux}");	hxy_g2->SetYTitle("y_{flux}");
	hxy_g2->GetXaxis()->SetTitleSize(0.06); hxy_g2->GetYaxis()->SetTitleSize(0.06);
	hxy_g2->GetXaxis()->CenterTitle(1);	hxy_g2->GetYaxis()->CenterTitle(1);

	TH1F *hr_e1 = new TH1F("hr_e1","e (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD", 80,70,150);
	hr_e1->SetXTitle("r_{flux}");
	hr_e1->GetXaxis()->CenterTitle(1);	
	hr_e1->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr_e2 = new TH1F("hr_e2","e (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD (EC triggered)", 80,70,150);
	hr_e2->SetXTitle("r_{flux}");
	hr_e2->GetXaxis()->CenterTitle(1);	
	hr_e2->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr_g1 = new TH1F("hr_g1","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD", 80,70,150);
	hr_g1->SetXTitle("r_{flux}");
	hr_g1->GetXaxis()->CenterTitle(1);	
	hr_g1->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr_g2 = new TH1F("hr_g2","e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD (EC triggered)", 80,70,150);
	hr_g2->SetXTitle("r_{flux}");
	hr_g2->GetXaxis()->CenterTitle(1);	
	hr_g2->GetXaxis()->SetTitleSize(0.06);

	TH2F *hxy1 = new TH2F("hxy1","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 100,-150,150,100,-150,150);
   	hxy1->SetXTitle("x_{flux}");	      hxy1->SetYTitle("y_{flux}");
	hxy1->GetXaxis()->SetTitleSize(0.06); hxy1->GetYaxis()->SetTitleSize(0.06);
	hxy1->GetXaxis()->CenterTitle(1);	  hxy1->GetYaxis()->CenterTitle(1);

	TH2F *hxy2 = new TH2F("hxy2","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on LASPD (EC triggered)", 100,-150,150,100,-150,150);
   	hxy2->SetXTitle("x_{flux}");	      hxy2->SetYTitle("y_{flux}");
	hxy2->GetXaxis()->SetTitleSize(0.06); hxy2->GetYaxis()->SetTitleSize(0.06);
	hxy2->GetXaxis()->CenterTitle(1);	  hxy2->GetYaxis()->CenterTitle(1);

	TH1F *hr1 = new TH1F("hr1","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD", 80,70,150);
	hr1->SetXTitle("r_{flux}");
	hr1->GetXaxis()->CenterTitle(1);	
	hr1->GetXaxis()->SetTitleSize(0.06);
	TH1F *hr2 = new TH1F("hr2","e^{#pm}+e^{#pm}(#gamma) (DIS+#pi^{0}+#pi^{#pm}+p+EM) on SPD (EC triggered)", 80,70,150);
	hr2->SetXTitle("r_{flux}");
	hr2->GetXaxis()->CenterTitle(1);	
	hr2->GetXaxis()->SetTitleSize(0.06);

	TH2F *hEDep_bgr = new TH2F("hEDep_bgr","EDep from e^{#pm}+e^{#pm}(#gamma) (#pi^{0}+#pi^{#pm}+p+EM) on LASPD", 80,70,150,100,0,10);
   	hEDep_bgr->SetXTitle("r_{flux}");	      hEDep_bgr->SetYTitle("E_{Dep} (MeV)");
	hEDep_bgr->GetXaxis()->SetTitleSize(0.06); hEDep_bgr->GetYaxis()->SetTitleSize(0.06);
	hEDep_bgr->GetXaxis()->CenterTitle(1);	  hEDep_bgr->GetYaxis()->CenterTitle(1);

	/*}}}*/
	
	int Nevt = T->GetEntries();
    double rate_e1=0, rate_e2=0, rate_g1=0,rate_g2=0;
    double EDep_e1=0, EDep_e2=0, EDep_g1=0,EDep_g2=0;
	for(int i=0;i<Nevt;i++){
		T->GetEntry(i);
		rate /= 1000.0;
		if(abs(PID_flux)==11){
			hxy_e1->Fill(x_flux,y_flux,rate*convert);
			hr_e1->Fill(r_flux,rate*convert);
			rate_e1 += rate*convert;
			EDep_e1 += EDep*1000.;
			if(R_EC>84&&R_EC<=140){
				hxy_e2->Fill(x_flux,y_flux,rate*convert*EC_Cut);
				hr_e2->Fill(r_flux,rate*convert*EC_Cut);
				rate_e2 += 2.0*rate*convert*EC_Cut;
				EDep_e2 += EDep*1000.;
			}
		}
		if(PID_flux==22){
			hxy_g1->Fill(x_flux,y_flux,rate*convert);
			hr_g1->Fill(r_flux,rate*convert);
			rate_g1 += 2.0*rate*convert;//One photon becomes two electrons
			EDep_g1 += EDep*1000.;
			if(R_EC>84&&R_EC<=140){
				hxy_g2->Fill(x_flux,y_flux,rate*convert*EC_Cut);
				hr_g2->Fill(r_flux,rate*convert*EC_Cut);
				rate_g2 += 2.0*rate*convert*EC_Cut;//One photon becomes two electrons
				EDep_g2 += EDep*1000.;
			}
		}
			
		if(abs(PID_flux)==11||PID_flux==22){
			hxy1->Fill(x_flux,y_flux,rate*convert);
			hr1->Fill(r_flux,rate*convert);
			if(R_EC>84&&R_EC<=140){
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
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_e1));
	c1->cd(3); hxy_e2->Draw("colz");	
	c1->cd(4); hr_e2->Draw("");	
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_e2));

	TCanvas *c2 = new TCanvas("c2","c2",1000,1000);
	c2->Divide(2,2);
	c2->cd(1); hxy_g1->Draw("colz");	
	c2->cd(2); hr_g1->Draw("");	
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_g1));
	c2->cd(3); hxy_g2->Draw("colz");	
	c2->cd(4); hr_g2->Draw("");	
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_g2));

	TCanvas *c3 = new TCanvas("c3","c3",1000,1000);
	c3->Divide(2,2);
	c3->cd(1); hxy1->Draw("colz");	
	c3->cd(2); hr1->Draw("");	
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_e1+rate_g1));
	c3->cd(3); hxy2->Draw("colz");	
	c3->cd(4); hr2->Draw("");	
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_e2+rate_g2));

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
	t1->DrawLatex(0.25,0.5,Form("#scale[0.9]{Rate = %e (KHz)}", rate_e1+rate_g1));

}

