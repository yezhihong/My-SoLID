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
	T[0]->Add("LASPD_background_SIDIS_He3_EM.root");
//	T[1]= new TChain("T");
//	T[1]->Add("LASPD_background_SIDIS_He3_pip.root");
//	T[2]= new TChain("T");
//	T[2]->Add("LASPD_background_SIDIS_He3_pip_up.root");
//	T[3]= new TChain("T");
//	T[3]->Add("LASPD_background_SIDIS_He3_pip_down.root");
//	T[4]= new TChain("T");
//	T[4]->Add("LASPD_background_SIDIS_He3_pim.root");
//	T[5]= new TChain("T");
//	T[5]->Add("LASPD_background_SIDIS_He3_pim_up.root");
//	T[6]= new TChain("T");
//	T[6]->Add("LASPD_background_SIDIS_He3_pim_down.root");
//	T[7]= new TChain("T");
//	T[7]->Add("LASPD_background_SIDIS_He3_p.root");
//	T[8]= new TChain("T");
//	T[8]->Add("LASPD_background_SIDIS_He3_p_up.root");
//	T[9]= new TChain("T");
//	T[9]->Add("LASPD_background_SIDIS_He3_p_down.root");
//	T[10]= new TChain("T");
//	T[10]->Add("LASPD_background_SIDIS_He3_pi0.root");
//	T[11]= new TChain("T");
//	T[11]->Add("LASPD_background_SIDIS_He3_pi0_up.root");
//	T[12]= new TChain("T");
//	T[12]->Add("LASPD_background_SIDIS_He3_pi0_down.root");

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
	TH1F *hr_ele = new TH1F("hr_ele","e+/- radius on SPD", 58,82.5,140.5);
	hr_ele->SetXTitle("r_{flux} (cm)");
	hr_ele->GetXaxis()->CenterTitle(1);	
	hr_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_pho = new TH1F("hr_pho","photon radius on SPD", 58,82.5,140.5);
	hr_pho->SetXTitle("r_{flux} (cm)");
	hr_pho->GetXaxis()->CenterTitle(1);	
	hr_pho->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_pip = new TH1F("hr_pip","pion+ radius on SPD", 58,82.5,140.5);
	hr_pip->SetXTitle("r_{flux} (cm)");
	hr_pip->GetXaxis()->CenterTitle(1);	
	hr_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_pim = new TH1F("hr_pim","pion- radius on SPD", 58,82.5,140.5);
	hr_pim->SetXTitle("r_{flux} (cm)");
	hr_pim->GetXaxis()->CenterTitle(1);	
	hr_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_kp = new TH1F("hr_kp","kaon+ radius on SPD", 58,82.5,140.5);
	hr_kp->SetXTitle("r_{flux} (cm)");
	hr_kp->GetXaxis()->CenterTitle(1);	
	hr_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_km = new TH1F("hr_km","kaon- radius on SPD", 58,82.5,140.5);
	hr_km->SetXTitle("r_{flux} (cm)");
	hr_km->GetXaxis()->CenterTitle(1);	
	hr_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *hr_pro = new TH1F("hr_pro","proton radius on SPD", 58,82.5,140.5);
	hr_pro->SetXTitle("r_{flux} (cm)");
	hr_pro->GetXaxis()->CenterTitle(1);	
	hr_pro->GetXaxis()->SetTitleSize(0.06);
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

	TH1F *hx_pip = new TH1F("hx_pip","pion+ x-axis on SPD", 80,-145,145);
	hx_pip->SetXTitle("x_{flux} (cm)");
	hx_pip->GetXaxis()->CenterTitle(1);	
	hx_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *hx_pim = new TH1F("hx_pim","pion- x-axis on SPD", 80,-145,145);
	hx_pim->SetXTitle("x_{flux} (cm)");
	hx_pim->GetXaxis()->CenterTitle(1);	
	hx_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *hx_kp = new TH1F("hx_kp","kaon+ x-axis on SPD", 80,-145,145);
	hx_kp->SetXTitle("x_{flux} (cm)");
	hx_kp->GetXaxis()->CenterTitle(1);	
	hx_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *hx_km = new TH1F("hx_km","kaon- x-axis on SPD", 80,-145,145);
	hx_km->SetXTitle("x_{flux} (cm)");
	hx_km->GetXaxis()->CenterTitle(1);	
	hx_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *hx_pro = new TH1F("hx_pro","proton x-axis on SPD", 80,-145,145);
	hx_pro->SetXTitle("x_{flux} (cm)");
	hx_pro->GetXaxis()->CenterTitle(1);	
	hx_pro->GetXaxis()->SetTitleSize(0.06);
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

	TH1F *hy_pip = new TH1F("hy_pip","pion+ y-axis on SPD", 80,-145,145);
	hy_pip->SetXTitle("y_{flux} (cm)");
	hy_pip->GetXaxis()->CenterTitle(1);	
	hy_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *hy_pim = new TH1F("hy_pim","pion- y-axis on SPD", 80,-145,145);
	hy_pim->SetXTitle("y_{flux} (cm)");
	hy_pim->GetXaxis()->CenterTitle(1);	
	hy_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *hy_kp = new TH1F("hy_kp","kaon+ y-axis on SPD", 80,-145,145);
	hy_kp->SetXTitle("y_{flux} (cm)");
	hy_kp->GetXaxis()->CenterTitle(1);	
	hy_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *hy_km = new TH1F("hy_km","kaon- y-axis on SPD", 80,-145,145);
	hy_km->SetXTitle("y_{flux} (cm)");
	hy_km->GetXaxis()->CenterTitle(1);	
	hy_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *hy_pro = new TH1F("hy_pro","proton y-axis on SPD", 80,-145,145);
	hy_pro->SetXTitle("y_{flux} (cm)");
	hy_pro->GetXaxis()->CenterTitle(1);	
	hy_pro->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/

	/*E_flux log{{{*/
	TH1F *he_ele = new TH1F("he_ele","e+/- energy (log) on SPD", 100,-6,4.5);
	he_ele->SetXTitle("E_{flux} (log) (MeV)");
	he_ele->GetXaxis()->CenterTitle(1);	
	he_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *he_pho = new TH1F("he_pho","photon energy (log) on SPD", 100,-6,4.5);
	he_pho->SetXTitle("E_{flux} (log) (MeV)");
	he_pho->GetXaxis()->CenterTitle(1);	
	he_pho->GetXaxis()->SetTitleSize(0.06);

	TH1F *he_pip = new TH1F("he_pip","pion+ energy (log) on SPD", 100,-6,4.5);
	he_pip->SetXTitle("E_{flux} (log) (MeV)");
	he_pip->GetXaxis()->CenterTitle(1);	
	he_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *he_pim = new TH1F("he_pim","pion- energy (log) on SPD", 100,-6,4.5);
	he_pim->SetXTitle("E_{flux} (log) (MeV)");
	he_pim->GetXaxis()->CenterTitle(1);	
	he_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *he_kp = new TH1F("he_kp","kaon+ energy (log) on SPD", 100,-6,4.5);
	he_kp->SetXTitle("E_{flux} (log) (MeV)");
	he_kp->GetXaxis()->CenterTitle(1);	
	he_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *he_km = new TH1F("he_km","kaon- energy (log) on SPD", 100,-6,4.5);
	he_km->SetXTitle("E_{flux} (log) (MeV)");
	he_km->GetXaxis()->CenterTitle(1);	
	he_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *he_pro = new TH1F("he_pro","proton energy (log) on SPD", 100,-6,4.5);
	he_pro->SetXTitle("E_{flux} (log) (MeV)");
	he_pro->GetXaxis()->CenterTitle(1);	
	he_pro->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/

	/*E_flux{{{*/
	TH1F *hE_ele = new TH1F("hE_ele","e+/- energy on SPD", 100,0,11000);
	hE_ele->SetXTitle("E_{flux} (MeV)");
	hE_ele->GetXaxis()->CenterTitle(1);	
	hE_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_pho = new TH1F("hE_pho","photon energy on SPD", 100,0.,11000.);
	hE_pho->SetXTitle("E_{flux} (MeV)");
	hE_pho->GetXaxis()->CenterTitle(1);	
	hE_pho->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_pip = new TH1F("hE_pip","pion+ energy on SPD", 100,0.,11000.);
	hE_pip->SetXTitle("E_{flux} (MeV)");
	hE_pip->GetXaxis()->CenterTitle(1);	
	hE_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_pim = new TH1F("hE_pim","pion- energy on SPD", 100,0.,11000.);
	hE_pim->SetXTitle("E_{flux} (MeV)");
	hE_pim->GetXaxis()->CenterTitle(1);	
	hE_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_kp = new TH1F("hE_kp","kaon+ energy on SPD", 100,0.,11000.);
	hE_kp->SetXTitle("E_{flux} (MeV)");
	hE_kp->GetXaxis()->CenterTitle(1);	
	hE_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_km = new TH1F("hE_km","kaon- energy on SPD", 100,0.,11000.);
	hE_km->SetXTitle("E_{flux} (MeV)");
	hE_km->GetXaxis()->CenterTitle(1);	
	hE_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *hE_pro = new TH1F("hE_pro","proton energy on SPD", 100,0.,11000.);
	hE_pro->SetXTitle("E_{flux} (MeV)");
	hE_pro->GetXaxis()->CenterTitle(1);	
	hE_pro->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/

	/*EDep_flux (log){{{*/
	TH1F *hedep_ele = new TH1F("hedep_ele","e+/- energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_ele->SetXTitle("Edep_{flux} (log) (MeV)");
	hedep_ele->GetXaxis()->CenterTitle(1);	
	hedep_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hedep_pho = new TH1F("hedep_pho","photon energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_pho->SetXTitle("Edep_{flux} (log) MeV)");
	hedep_pho->GetXaxis()->CenterTitle(1);	
	hedep_pho->GetXaxis()->SetTitleSize(0.06);

	TH1F *hedep_pip = new TH1F("hedep_pip","pion+ energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_pip->SetXTitle("Edep_{flux} (log) (MeV)");
	hedep_pip->GetXaxis()->CenterTitle(1);	
	hedep_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *hedep_pim = new TH1F("hedep_pim","pion- energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_pim->SetXTitle("Edep_{flux} (log) (MeV)");
	hedep_pim->GetXaxis()->CenterTitle(1);	
	hedep_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *hedep_kp = new TH1F("hedep_kp","kaon+ energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_kp->SetXTitle("Edep_{flux} (log) (MeV)");
	hedep_kp->GetXaxis()->CenterTitle(1);	
	hedep_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *hedep_km = new TH1F("hedep_km","kaon- energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_km->SetXTitle("Edep_{flux} (log) (MeV)");
	hedep_km->GetXaxis()->CenterTitle(1);	
	hedep_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *hedep_pro = new TH1F("hedep_pro","proton energy-deposit (log) on SPD", 100,-6,1.5);
	hedep_pro->SetXTitle("Edep_{flux} (log) (MeV)");
	hedep_pro->GetXaxis()->CenterTitle(1);	
	hedep_pro->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/
	
	/*EDep_flux{{{*/
	TH1F *hEdep_ele = new TH1F("hEdep_ele","e+/- energy-deposit on SPD", 100,0.,6.0);
	hEdep_ele->SetXTitle("Edep_{flux} (MeV)");
	hEdep_ele->GetXaxis()->CenterTitle(1);	
	hEdep_ele->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_pho = new TH1F("hEdep_pho","photon energy-deposit on SPD", 100,0.,6.0);
	hEdep_pho->SetXTitle("Edep_{flux} MeV)");
	hEdep_pho->GetXaxis()->CenterTitle(1);	
	hEdep_pho->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_pip = new TH1F("hEdep_pip","pion+ energy-deposit on SPD", 100,0.,6.0);
	hEdep_pip->SetXTitle("Edep_{flux} (MeV)");
	hEdep_pip->GetXaxis()->CenterTitle(1);	
	hEdep_pip->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_pim = new TH1F("hEdep_pim","pion- energy-deposit on SPD", 100,0.,6.0);
	hEdep_pim->SetXTitle("Edep_{flux} (MeV)");
	hEdep_pim->GetXaxis()->CenterTitle(1);	
	hEdep_pim->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_kp = new TH1F("hEdep_kp","kaon+ energy-deposit on SPD", 100,0.,6.0);
	hEdep_kp->SetXTitle("Edep_{flux} (MeV)");
	hEdep_kp->GetXaxis()->CenterTitle(1);	
	hEdep_kp->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_km = new TH1F("hEdep_km","kaon- energy-deposit on SPD", 100,0.,6.0);
	hEdep_km->SetXTitle("Edep_{flux} (MeV)");
	hEdep_km->GetXaxis()->CenterTitle(1);	
	hEdep_km->GetXaxis()->SetTitleSize(0.06);

	TH1F *hEdep_pro = new TH1F("hEdep_pro","proton energy-deposit on SPD", 100,0.,6.0);
	hEdep_pro->SetXTitle("Edep_{flux} (MeV)");
	hEdep_pro->GetXaxis()->CenterTitle(1);	
	hEdep_pro->GetXaxis()->SetTitleSize(0.06);
	/*}}}*/
	/*}}}*/

	for(int j=0;j<N;j++){	
		int Nevt = T[j]->GetEntries();
		for(int i=0;i<Nevt;i++){
			T[j]->GetEntry(i);
			cerr<<"---- Event = "<< i <<"\r";
			//All Charged Particles
			if(abs(PID_flux)==11){
			//if(fabs(PID_flux-11)<1.0||fabs(PID_flux+11)<1.0){
				convert = 1.0;
				hx_ele->Fill(x_flux,rate*convert);
				hy_ele->Fill(y_flux,rate*convert);
				hr_ele->Fill(r_flux,rate*convert);
				he_ele->Fill(log10(E*1000.),rate*convert);
				hedep_ele->Fill(log10(EDep*1000.),rate*convert);
				hE_ele->Fill((E*1000.),rate*convert);
				hEdep_ele->Fill((EDep*1000.),rate*convert);
				}
			else if(PID_flux==22){
			//else if(fabs(PID_flux-22)<1.0){
				//rate *= 2.0;//One photon becomes two electrons
				convert = 1.0;
				hx_pho->Fill(x_flux,rate*convert);
				hy_pho->Fill(y_flux,rate*convert);
				hr_pho->Fill(r_flux,rate*convert);
				he_pho->Fill(log10(E*1000.),rate*convert);
				hedep_pho->Fill(log10(EDep*1000.),rate*convert);
				hE_pho->Fill((E*1000.),rate*convert);
				hEdep_pho->Fill((EDep*1000.),rate*convert);
				}
			else if(PID_flux==211&&j!=9){
			//else if(fabs(PID_flux-211)<1.0){
				convert = 1.0;
				hx_pip->Fill(x_flux,rate*convert);
				hy_pip->Fill(y_flux,rate*convert);
				hr_pip->Fill(r_flux,rate*convert);
				he_pip->Fill(log10(E*1000.),rate*convert);
				hedep_pip->Fill(log10(EDep*1000.),rate*convert);
				hE_pip->Fill((E*1000.),rate*convert);
				hEdep_pip->Fill((EDep*1000.),rate*convert);
				}
			else if(PID_flux==-211&&j!=9){
			//else if(fabs(PID_flux+211)<1.0){
				convert = 1.0;
				hx_pim->Fill(x_flux,rate*convert);
				hy_pim->Fill(y_flux,rate*convert);
				hr_pim->Fill(r_flux,rate*convert);
				he_pim->Fill(log10(E*1000.),rate*convert);
				hedep_pim->Fill(log10(EDep*1000.),rate*convert);
				hE_pim->Fill((E*1000.),rate*convert);
				hEdep_pim->Fill((EDep*1000.),rate*convert);
				}
			else if(PID_flux==321&&j!=9){
			//else if(fabs(PID_flux-321)<1.0){
				convert = 1.0;
				hx_kp->Fill(x_flux,rate*convert);
				hy_kp->Fill(y_flux,rate*convert);
				hr_kp->Fill(r_flux,rate*convert);
				he_kp->Fill(log10(E*1000.),rate*convert);
				hedep_pip->Fill(log10(EDep*1000.),rate*convert);
				hE_kp->Fill((E*1000.),rate*convert);
				hEdep_pip->Fill((EDep*1000.),rate*convert);
				}
			else if(PID_flux==-321&&j!=9){
			//else if(fabs(PID_flux+321)<1.0){
				convert = 1.0;
				hx_km->Fill(x_flux,rate*convert);
				hy_km->Fill(y_flux,rate*convert);
				hr_km->Fill(r_flux,rate*convert);
				he_km->Fill(log10(E*1000.),rate*convert);
				hedep_km->Fill(log10(EDep*1000.),rate*convert);
				hE_km->Fill((E*1000.),rate*convert);
				hEdep_km->Fill((EDep*1000.),rate*convert);
				}
			else if(abs(PID_flux)==2212&&j!=9){
			//else if(fabs(PID_flux-2212)<1.0||fabs(PID_flux+2212)<1.0){
				convert = 1.0;
				hx_pro->Fill(x_flux,rate*convert);
				hy_pro->Fill(y_flux,rate*convert);
				hr_pro->Fill(r_flux,rate*convert);
				he_pro->Fill(log10(E*1000.),rate*convert);
				hedep_pro->Fill(log10(EDep*1000.),rate*convert);
				hE_pro->Fill((E*1000.),rate*convert);
				hEdep_pro->Fill((EDep*1000.),rate*convert);
				}
			else{
				cerr<<"---- None of e+/-, g, pi+/- and p, but PID = "<< PID_flux<<endl;
			}
		}
		}

	TFile *new_file = new TFile("background_histo_pi0.root","recreate");
	hr_ele->Write(); hx_ele->Write(); hy_ele->Write(); he_ele->Write(); hedep_ele->Write(); hE_ele->Write(); hEdep_ele->Write();
	hr_pho->Write(); hx_pho->Write(); hy_pho->Write(); he_pho->Write(); hedep_pho->Write(); hE_pho->Write(); hEdep_pho->Write();
	hr_pip->Write(); hx_pip->Write(); hy_pip->Write(); he_pip->Write(); hedep_pip->Write(); hE_pip->Write(); hEdep_pip->Write();
	hr_pim->Write(); hx_pim->Write(); hy_pim->Write(); he_pim->Write(); hedep_pim->Write(); hE_pim->Write(); hEdep_pim->Write();
	hr_kp->Write(); hx_kp->Write(); hy_kp->Write(); he_kp->Write(); hedep_kp->Write(); hE_kp->Write(); hEdep_kp->Write();
	hr_km->Write(); hx_km->Write(); hy_km->Write(); he_km->Write(); hedep_km->Write(); hE_km->Write(); hEdep_km->Write();
	hr_pro->Write(); hx_pro->Write(); hy_pro->Write(); he_pro->Write(); hedep_pro->Write(); hE_pro->Write(); hEdep_pro->Write();
	new_file->Close();

	for(int i=0;i<N;i++)
		delete T[i];

}

