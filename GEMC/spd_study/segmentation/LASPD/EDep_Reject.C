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

void EDep_Reject(){
	gStyle->SetOptStat(0);
	
	/*Background Tree and Histogram{{{*/
	TChain *T = new TChain("T");
	T->Add("LASPD_background_SIDIS_He3_EM.root");
	T->Add("LASPD_background_SIDIS_He3_pi0_down.root");
	T->Add("LASPD_background_SIDIS_He3_pi0.root");
	T->Add("LASPD_background_SIDIS_He3_pi0_up.root");
	T->Add("LASPD_background_SIDIS_He3_pip_down.root");
	T->Add("LASPD_background_SIDIS_He3_pip.root");
	T->Add("LASPD_background_SIDIS_He3_pip_up.root");
	T->Add("LASPD_background_SIDIS_He3_pim_down.root");
	T->Add("LASPD_background_SIDIS_He3_pim.root");
	T->Add("LASPD_background_SIDIS_He3_pim_up.root");
	T->Add("LASPD_background_SIDIS_He3_p_down.root");
	T->Add("LASPD_background_SIDIS_He3_p.root");
	T->Add("LASPD_background_SIDIS_He3_p_up.root");

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

	TH2F *he_bck = new TH2F("he_bck","EDep (#pi^{0}+#pi^{#pm}+p+EM)) on SPD", 80,70,150, 100, 0.0,10.0);
	he_bck->SetXTitle("r_{flux}");
	he_bck->GetXaxis()->CenterTitle(1);	
	he_bck->GetXaxis()->SetTitleSize(0.06);
	he_bck->SetYTitle("E_{Dep}");
	he_bck->GetYaxis()->CenterTitle(1);	
	he_bck->GetYaxis()->SetTitleSize(0.06);
	/*}}}*/

	const double nSec = 1e-9;
	const double MIP = 5.35;//MeV
    double EDep_Sum[55]= {55*0};
    double area[55] = {55*0};//cm^2

	/*Loop Background Events{{{*/
	int Nevt = T->GetEntries();
	for(int i=0;i<Nevt;i++){
		cerr<<" Evnt# = "<<i<<"\r";
		T->GetEntry(i);
		for(int j=0;j<55;j++){
			double R = 80+j;//cm
		    if(r_flux>=R && r_flux<R+1 && Module==1 && PID_flux == 22){
				area[j] = 3.1415926 * ( (R+1)*(R+1)-R*R );//cm^2
				EDep_Sum[j]+=2.*EDep*convert*rate*1e3*area[j];	//MeV, rate is ~ cm^-2 * s^-1
			}
		}
	}
    /*}}}*/
    
	double EDep_Total = 0.0;
	for(int j=0;j<55;j++){
        EDep_Total += EDep_Sum[j];
	}
	cerr<<Form("---- EDep_Bgd = %f, We need = %d ",EDep_Total*40*nSec,(int)(0.5+EDep_Total*40*nSec/MIP * 2.0 * 30) )<<endl;

}

