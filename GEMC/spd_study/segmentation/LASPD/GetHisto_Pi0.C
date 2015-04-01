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

void GetHisto(){
	gStyle->SetOptStat(0);

	const int N = 1;

	TChain *T[N];
	T[0]= new TChain("T");
	T[0]->Add("../LASPD_background_SIDIS_He3_pi0_mpid.root");
	T[0]->Add("../LASPD_background_SIDIS_He3_pi0_up_mpid.root");
	T[0]->Add("../LASPD_background_SIDIS_He3_pi0_down_mpid.root");

	/*Tree and Histogram{{{*/
	double rate,mom_flux,theta_flux,phi_flux;
	double r_flux,x_flux,y_flux,z_flux,px_flux,py_flux,pz_flux,vx_flux,vy_flux,vz_flux;
	double R_EC, EDep, E, convert, EC_Cut;
	int evn,evn_flux,nhit_flux,ID_flux, PID_flux,ID_Pick, Module,Slide;

	/*Define Braches{{{*/
	for(int i=0;i<N;i++){
		T[i]->SetBranchAddress("mom_flux", &mom_flux);
		T[i]->SetBranchAddress("theta_flux", &theta_flux);
		T[i]->SetBranchAddress("phi_flux",&phi_flux);
		T[i]->SetBranchAddress("r_flux",&r_flux);
		T[i]->SetBranchAddress("x_flux",&x_flux);
		T[i]->SetBranchAddress("y_flux",&y_flux);
		T[i]->SetBranchAddress("z_flux",&z_flux);
		T[i]->SetBranchAddress("px_flux",&px_flux);
		T[i]->SetBranchAddress("py_flux",&py_flux);
		T[i]->SetBranchAddress("pz_flux",&pz_flux);
		T[i]->SetBranchAddress("vx_flux",&vx_flux);
		T[i]->SetBranchAddress("vy_flux",&vy_flux);
		T[i]->SetBranchAddress("vz_flux",&vz_flux);
		T[i]->SetBranchAddress("R_EC",&R_EC);
		T[i]->SetBranchAddress("EC_Cut",&EC_Cut);
		T[i]->SetBranchAddress("EDep",&EDep);
		T[i]->SetBranchAddress("E",&E);
		T[i]->SetBranchAddress("rate", &rate);
		T[i]->SetBranchAddress("convert",&convert);
		T[i]->SetBranchAddress("ID_flux",&ID_flux);
		T[i]->SetBranchAddress("PID_flux",&PID_flux);
		T[i]->SetBranchAddress("ID_Pick",&ID_Pick);
		T[i]->SetBranchAddress("Module",&Module);
		T[i]->SetBranchAddress("Slide",&Slide);
		T[i]->SetBranchAddress("evn",  &evn);
		T[i]->SetBranchAddress("evn_flux", &evn_flux);
		T[i]->SetBranchAddress("nhit_flux", &nhit_flux);
	}
	/*}}}*/

	/*R_flux{{{*/
	TH1F *hr_ele = new TH1F("hr_ele","e+/- radius on SPD", 65,75.,140.);
	hr_ele->SetXTitle("r_{flux} (cm)");
	hr_ele->GetXaxis()->CenterTitle(1);	
	hr_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_pho = new TH1F("hr_pho","photon radius on SPD", 65,75.,140.);
	hr_pho->SetXTitle("r_{flux} (cm)");
	hr_pho->GetXaxis()->CenterTitle(1);	
	hr_pho->GetXaxis()->SetTitleSize(0.06);

	/*}}}*/
	
	/*Phi_flux{{{*/
	TH1F *hPhi_ele = new TH1F("hPhi_ele","e+/- Phi on SPD", 360,0.,360.);
	hPhi_ele->SetXTitle("#phi_{flux} (degree)");
	hPhi_ele->GetXaxis()->CenterTitle(1);	
	hPhi_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hPhi_pho = new TH1F("hPhi_pho","photon Phi on SPD",360,0.,360.);
	hPhi_pho->SetXTitle("#phi_{flux} (degree)");
	hPhi_pho->GetXaxis()->CenterTitle(1);	
	hPhi_pho->GetXaxis()->SetTitleSize(0.06);

	/*}}}*/
	
	/*X_flux{{{*/
	TH1F *hx_ele = new TH1F("hx_ele","e+/- x-axis on SPD", 80,-145,145);
	hx_ele->SetXTitle("x_{flux} (cm)");
	hx_ele->GetXaxis()->CenterTitle(1);	
	hx_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hx_pho = new TH1F("hx_pho","photon x-axis on SPD", 80,-145,145);
	hx_pho->SetXTitle("x_{flux} (cm)");
	hx_pho->GetXaxis()->CenterTitle(1);	
	hx_pho->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/

	/*Y_flux{{{*/
	TH1F *hy_ele = new TH1F("hy_ele","e+/- y-axis on SPD", 80,-145,145);
	hy_ele->SetXTitle("y_{flux} (cm)");
	hy_ele->GetXaxis()->CenterTitle(1);	
	hy_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hy_pho = new TH1F("hy_pho","photon y-axis on SPD", 80,-145,145);
	hy_pho->SetXTitle("y_{flux} (cm)");
	hy_pho->GetXaxis()->CenterTitle(1);	
	hy_pho->GetXaxis()->SetTitleSize(0.06);

	/*}}}*/

	/*E_flux log{{{*/
	TH1F *hELog_ele = new TH1F("hELog_ele","e+/- energy (log) on SPD", 100,-6,3.1);
	hELog_ele->SetXTitle("E_{flux} (MeV) (log)");
	hELog_ele->GetXaxis()->CenterTitle(1);	
	hELog_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hELog_pho = new TH1F("hELog_pho","photon energy (log) on SPD", 100,-6,3.1);
	hELog_pho->SetXTitle("E_{flux} (log) (MeV)");
	hELog_pho->GetXaxis()->CenterTitle(1);	
	hELog_pho->GetXaxis()->SetTitleSize(0.06);

	/*}}}*/
	
	/*E_flux log vs r{{{*/
	TH2F *hRELog_ele = new TH2F("hRELog_ele","e+/- energy (log) on SPD", 65,75, 140, 100,-6,3.1);
	hRELog_ele->SetXTitle("R_{flux} (cm)");
	hRELog_ele->SetYTitle("E_{flux} (MeV) (log)");
	hRELog_ele->GetXaxis()->CenterTitle(1);	
	hRELog_ele->GetXaxis()->SetTitleSize(0.06);
	hRELog_ele->GetYaxis()->CenterTitle(1);	
	hRELog_ele->GetYaxis()->SetTitleSize(0.06);
	
	TH2F *hRELog_pho = new TH2F("hRELog_pho","photon energy (log) on SPD",65,75,140, 100,-6,3.1);
	hRELog_pho->SetXTitle("R_{flux} (cm)");
	hRELog_pho->SetYTitle("E_{flux} (MeV) (log)");
	hRELog_pho->GetXaxis()->CenterTitle(1);	
	hRELog_pho->GetXaxis()->SetTitleSize(0.06);
	hRELog_pho->GetYaxis()->CenterTitle(1);	
	hRELog_pho->GetYaxis()->SetTitleSize(0.06);
		/*}}}*/

	/*E_flux{{{*/
	TH1F *hE_ele = new TH1F("hE_ele","e+/- energy on SPD", 100,1e-6,1e4);
	hE_ele->SetXTitle("E_{flux} (MeV)");
	hE_ele->GetXaxis()->CenterTitle(1);	
	hE_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_pho = new TH1F("hE_pho","photon energy on SPD", 100,1e-6,1e4);
	hE_pho->SetXTitle("E_{flux} (MeV)");
	hE_pho->GetXaxis()->CenterTitle(1);	
	hE_pho->GetXaxis()->SetTitleSize(0.06);

	/*}}}*/

	/*EDep_flux (log){{{*/
	TH1F *hEdepLog_ele = new TH1F("hEdepLog_ele","e+/- energy-deposit (log) on SPD", 100,-6,1.0);
	hEdepLog_ele->SetXTitle("Edep_{flux} (log) (MeV)");
	hEdepLog_ele->GetXaxis()->CenterTitle(1);	
	hEdepLog_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdepLog_pho = new TH1F("hEdepLog_pho","photon energy-deposit (log) on SPD", 100,-6,1.0);
	hEdepLog_pho->SetXTitle("Edep_{flux} (log) MeV)");
	hEdepLog_pho->GetXaxis()->CenterTitle(1);	
	hEdepLog_pho->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/
	
	/*EDep_flux (log) vs E(log){{{*/
	TH2F *hELogEdepLog_ele = new TH2F("hELogEdepLog_ele","e+/- energy-deposit (log) vs E(log) on SPD", 100,-6,3.1,100,-6,1.0);
	hELogEdepLog_ele->SetXTitle("E_{flux} (log(MeV))");
	hELogEdepLog_ele->GetXaxis()->CenterTitle(1);	
	hELogEdepLog_ele->GetXaxis()->SetTitleSize(0.06);
	hELogEdepLog_ele->SetYTitle("Edep_{flux} (log) (MeV)");
	hELogEdepLog_ele->GetYaxis()->CenterTitle(1);	
	hELogEdepLog_ele->GetYaxis()->SetTitleSize(0.06);

	TH2F *hELogEdepLog_pho = new TH2F("hELogEdepLog_pho","photon energy-deposit (log) vs E(log) on SPD", 100,-6,3.1,100,-6,1.0);
	hELogEdepLog_pho->SetXTitle("E_{flux} (log(MeV))");
	hELogEdepLog_pho->GetXaxis()->CenterTitle(1);	
	hELogEdepLog_pho->GetXaxis()->SetTitleSize(0.06);
	hELogEdepLog_pho->SetYTitle("Edep_{flux} (log) (MeV)");
	hELogEdepLog_pho->GetYaxis()->CenterTitle(1);	
	hELogEdepLog_pho->GetYaxis()->SetTitleSize(0.06);

	/*}}}*/
	
	/*EDep_flux (log) vs r{{{*/
	TH2F *hREdepLog_ele = new TH2F("hREdepLog_ele","e+/- energy-deposit (log) vs r on SPD", 65,75.,140.,100,-6,1.0);
	hREdepLog_ele->SetXTitle("R_{flux} (cm)");
	hREdepLog_ele->GetXaxis()->CenterTitle(1);	
	hREdepLog_ele->GetXaxis()->SetTitleSize(0.06);
	hREdepLog_ele->SetYTitle("Edep_{flux} (log) (MeV)");
	hREdepLog_ele->GetYaxis()->CenterTitle(1);	
	hREdepLog_ele->GetYaxis()->SetTitleSize(0.06);

	TH2F *hREdepLog_pho = new TH2F("hREdepLog_pho","photon energy-deposit (log) vs r on SPD", 65,75.,140.,100,-6,1.0);
	hREdepLog_pho->SetXTitle("R_{flux} (cm)");
	hREdepLog_pho->GetXaxis()->CenterTitle(1);	
	hREdepLog_pho->GetXaxis()->SetTitleSize(0.06);
	hREdepLog_pho->SetYTitle("Edep_{flux} (log) MeV)");
	hREdepLog_pho->GetYaxis()->CenterTitle(1);	
	hREdepLog_pho->GetYaxis()->SetTitleSize(0.06);
	/*}}}*/
	
	/*EDep_flux (log) vs phi{{{*/
	TH2F *hPhiEdepLog_ele = new TH2F("hPhiEdepLog_ele","e+/e- energy-deposit (log) vs r on SPD",360,0.,360.,100,-6,1.0);
	hPhiEdepLog_ele->SetXTitle("#phi_{flux} (degree)");
	hPhiEdepLog_ele->GetXaxis()->CenterTitle(1);	
	hPhiEdepLog_ele->GetXaxis()->SetTitleSize(0.06);
	hPhiEdepLog_ele->SetYTitle("Edep_{flux} (log) MeV)");
	hPhiEdepLog_ele->GetYaxis()->CenterTitle(1);	
	hPhiEdepLog_ele->GetYaxis()->SetTitleSize(0.06);

	TH2F *hPhiEdepLog_pho = new TH2F("hPhiEdepLog_pho","photon energy-deposit (log) vs r on SPD",360,0.,360.,100,-6,1.0);
	hPhiEdepLog_pho->SetXTitle("#phi_{flux} (degree)");
	hPhiEdepLog_pho->GetXaxis()->CenterTitle(1);	
	hPhiEdepLog_pho->GetXaxis()->SetTitleSize(0.06);
	hPhiEdepLog_pho->SetYTitle("Edep_{flux} (log) MeV)");
	hPhiEdepLog_pho->GetYaxis()->CenterTitle(1);	
	hPhiEdepLog_pho->GetYaxis()->SetTitleSize(0.06);
	/*}}}*/
	
	/*EDep_flux{{{*/
	TH1F *hEdep_ele = new TH1F("hEdep_ele","e+/- energy-deposit on SPD", 100,1e-6,10.0);
	hEdep_ele->SetXTitle("Edep_{flux} (MeV)");
	hEdep_ele->GetXaxis()->CenterTitle(1);	
	hEdep_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_pho = new TH1F("hEdep_pho","photon energy-deposit on SPD", 100,1e-6,10.0);
	hEdep_pho->SetXTitle("Edep_{flux} MeV)");
	hEdep_pho->GetXaxis()->CenterTitle(1);	
	hEdep_pho->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/

	/*}}}*/

	for(int j=0;j<N;j++){	
		int Nevt = T[j]->GetEntries();
		for(int i=0;i<Nevt;i++){
			T[j]->GetEntry(i);
			cerr<<"---- Event = "<< i <<"\r";

			E *= 1000.;
			EDep *= 1000.;
			//All Charged Particles
			if(abs(PID_flux)==11){
				//if(fabs(PID_flux-11)<1.0||fabs(PID_flux+11)<1.0){
				convert = 1.0;
				hx_ele->Fill(x_flux,rate);
				hy_ele->Fill(y_flux,rate);
				hr_ele->Fill(r_flux,rate);
				hPhi_ele->Fill(phi_flux,rate);
				hELog_ele->Fill(log10(E),rate);
				hE_ele->Fill((E),rate);
				hRELog_ele->Fill(r_flux,log10(E),rate);
				hEdepLog_ele->Fill(log10(EDep),rate);
				hELogEdepLog_ele->Fill(log10(E),log10(EDep),rate);
				hREdepLog_ele->Fill(r_flux,log10(EDep),rate);
				hPhiEdepLog_ele->Fill(phi_flux,log10(EDep),rate);
				hEdep_ele->Fill((EDep),rate);
			}
			else if(PID_flux==22){
				//else if(fabs(PID_flux-22)<1.0){
				//rate *= 2.0;//One photon becomes two electrons
				convert = 1.0;
				hx_pho->Fill(x_flux,rate);
				hy_pho->Fill(y_flux,rate);
				hr_pho->Fill(r_flux,rate);
				hPhi_pho->Fill(phi_flux,rate);
				hELog_pho->Fill(log10(E),rate);
				hE_pho->Fill((E),rate);
				hRELog_pho->Fill(r_flux,log10(E),rate);
				hEdepLog_pho->Fill(log10(EDep),rate);
				hELogEdepLog_pho->Fill(log10(E),log10(EDep),rate);
				hREdepLog_pho->Fill(r_flux,log10(EDep),rate);
				hPhiEdepLog_pho->Fill(phi_flux,log10(EDep),rate);
				hEdep_pho->Fill((EDep),rate);
			}
			}
			}

	TFile *new_file = new TFile("background_histo_pi0.root","recreate");
	hr_ele->Write(); hx_ele->Write(); hy_ele->Write(); hELog_ele->Write(); hEdepLog_ele->Write(); hE_ele->Write(); hEdep_ele->Write(); hRELog_ele->Write(); hREdepLog_ele->Write(); 
	hr_pho->Write(); hx_pho->Write(); hy_pho->Write(); hELog_pho->Write(); hEdepLog_pho->Write(); hE_pho->Write(); hEdep_pho->Write(); hRELog_pho->Write(); hREdepLog_pho->Write(); 
	hPhiEdepLog_ele->Write(); hPhiEdepLog_pho->Write();
	hELogEdepLog_ele->Write(); hELogEdepLog_pho->Write();
	hPhi_ele->Write(); hPhi_pho->Write();
	new_file->Close();

	for(int i=0;i<N;i++)
		delete T[i];

}

