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
#include <TClass.h>
#include <TPaletteAxis.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
#include "LHAPDF/LHAPDF.h"
//#include <TMatrix.h>
/*}}}*/

using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

int main(Int_t argc, char *argv[]){
	TString prefix = "./rootfiles/sidis_3he_";
	TString new_prefix = "./skim_rootfiles_11p8/3he_skim_";

	TString particle = "X"; cerr<<"-- What particle (pip,pim,kp,km)? "; cin >> particle;

	Int_t start_num  = 1;
	Int_t end_num  = 2;
	Int_t zflag = 0;
	Int_t Q2flag = 0;

	/*Define old root file{{{*/
	TString posfix,new_filename;
	TChain *T = new TChain("T","T");
	for (Int_t i=start_num; i<=end_num;i++){
		posfix.Form("_11_0_%d_0.root",i);
		new_filename = prefix + particle + posfix;
		cerr<<Form(" @@@ Adding Root File: %s", new_filename.Data())<<endl;
		T->AddFile(new_filename);
		posfix.Form("_8_0_%d_0.root",i);
		new_filename = prefix + particle + posfix;
		cerr<<Form("     Adding Root File: %s", new_filename.Data())<<endl;
		T->AddFile(new_filename);
	}
	cerr<<Form("   Got total number of events = %d", (int)(T->GetEntries()))<<endl;

	Double_t theta_gen, phi_gen,mom_gen;
	Double_t Q2,W,Wp,x,y,z,pt,nu,s;
	Double_t theta_q,phi_q;
	Double_t theta_s,phi_s,phi_h;
	Double_t jacoF,dxs_hp,dxs_hm;
	Double_t mom_ele,mom_had, weight_hp, weight_hm;
	Double_t theta_ele,theta_had;
	Double_t phi_ele,phi_had;
	Double_t mom_pro,energy_pro;
	Double_t mom_ini_ele,energy_ini_ele;
	Double_t dilute[2];
	Int_t nsim;

	T->SetBranchAddress("Q2",&Q2);
	T->SetBranchAddress("W",&W);
	T->SetBranchAddress("Wp",&Wp);
	T->SetBranchAddress("x",&x);
	T->SetBranchAddress("y",&y);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("nu",&nu);
	T->SetBranchAddress("s",&s);
	T->SetBranchAddress("pt",&pt);
	T->SetBranchAddress("theta_q",&theta_q);
	T->SetBranchAddress("theta_s",&theta_s);
	T->SetBranchAddress("phi_h",&phi_h);
	T->SetBranchAddress("phi_s",&phi_s);
	T->SetBranchAddress("jacoF",&jacoF);
	T->SetBranchAddress("dxs_hm",&dxs_hm);
	T->SetBranchAddress("dxs_hp",&dxs_hp);
	T->SetBranchAddress("weight_hp",&weight_hp);
	T->SetBranchAddress("weight_hm",&weight_hm);
	T->SetBranchAddress("mom_ele",&mom_ele);
	T->SetBranchAddress("mom_had",&mom_had);
	T->SetBranchAddress("theta_ele",&theta_ele);
	T->SetBranchAddress("theta_had",&theta_had);
	T->SetBranchAddress("phi_ele",&phi_ele);
	T->SetBranchAddress("phi_had",&phi_had);
	T->SetBranchAddress("nsim",&nsim);
	T->SetBranchAddress("dilute_p",&dilute[0]);
	T->SetBranchAddress("dilute_m",&dilute[1]);
	
	/*}}}*/

	//deal with acceptance files 
    //CLEO	
    TCanvas *dummycavan = new TCanvas();
	TFile *file_negative=new TFile("./solid_acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root","r");
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
	TH2F *accep_had=(TH2F*)accep_ele_forward->Clone("acceptance_hadron");//pion is the same as electron acceptance
	
	Double_t zmin=1000,zmax=-1000;
	Double_t Q2min=1000,Q2max=-1000.;
	Double_t weight = 0.0, ele_energy = 0.0;
	Double_t accp_forward = 0.0, accp_large = 0.0;
	for(zflag=1;zflag<=8;zflag++){
		for(Q2flag=1;Q2flag<=6;Q2flag++){
		//for(Q2flag=1;Q2flag<=6;Q2flag++){
			/*Get Z and Q2 Bin{{{*/
			if (zflag==1){
				zmin = 0.3; zmax = 0.35;
			}else if (zflag==2){
				zmin = 0.35; zmax = 0.4;
			}else if (zflag==3){
				zmin = 0.4; zmax = 0.45;
			}else if (zflag==4){
				zmin = 0.45; zmax = 0.5;
			}else if (zflag==5){
				zmin = 0.5; zmax = 0.55;
			}else if (zflag==6){
				zmin = 0.55; zmax = 0.6;
			}else if (zflag==7){
				zmin = 0.6; zmax = 0.65;
			}else if (zflag==8){
				zmin = 0.65; zmax = 0.7;
			}
			if (Q2flag==1){
				Q2min = 1.; Q2max = 2.0;
			}else if (Q2flag==2){
				Q2min = 2.0; Q2max = 3.0;
			}else if (Q2flag==3){
				Q2min = 3.0; Q2max = 4.0;
			}else if (Q2flag==4){
				Q2min = 4.0; Q2max = 5.0;
			}else if (Q2flag==5){
				Q2min = 5.0; Q2max = 6.0;
			}else if (Q2flag==6){
				Q2min = 6.0; Q2max = 8.0;
			}else if (Q2flag==7){
				Q2min = 8.0; Q2max = 10.0;
			}
			/*}}}*/

			/*Define new rootfile for each bin{{{*/
			posfix.Form("_%d_%d.root",zflag,Q2flag);
			new_filename = new_prefix + particle + posfix;
			TFile *file = new TFile(new_filename,"RECREATE");
			TTree *t1 = new TTree("T","T");
			t1->SetDirectory(file);

			t1->Branch("Q2",&Q2,"data/D");
			t1->Branch("W",&W,"data/D");
			t1->Branch("Wp",&Wp,"data/D");
			t1->Branch("x",&x,"data/D");
			t1->Branch("y",&y,"data/D");
			t1->Branch("z",&z,"data/D");
			t1->Branch("nu",&nu,"data/D");
			t1->Branch("s",&s,"data/D");
			t1->Branch("pt",&pt,"data/D");
			t1->Branch("theta_q",&theta_q,"data/D");
			t1->Branch("theta_s",&theta_s,"data/D");
			t1->Branch("phi_h",&phi_h,"data/D");
			t1->Branch("phi_s",&phi_s,"data/D");
			t1->Branch("jacoF",&jacoF,"jacoF/D");
			t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
			t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
			t1->Branch("weight_hp",&weight_hp,"weight_hp/D");
			t1->Branch("weight_hm",&weight_hm,"weight_hm/D");
			t1->Branch("mom_ele",&mom_ele,"mom_ele/D");
			t1->Branch("mom_had",&mom_had,"mom_had/D");
			t1->Branch("theta_ele",&theta_ele,"theta_ele/D");
			t1->Branch("theta_had",&theta_had,"theta_had/D");
			t1->Branch("phi_ele",&phi_ele,"phi_ele/D");
			t1->Branch("phi_had",&phi_had,"phi_had/D");
			t1->Branch("nsim",&nsim,"nsim/I");
			t1->Branch("dilute_p",&dilute[0],"data/D");
			t1->Branch("dilute_m",&dilute[1],"data/D");
			t1->Branch("ele_energy",&ele_energy,"ele_energy/D");
			t1->Branch("weight",&weight,"data/D");
			t1->Branch("accp_forward",&accp_forward,"data/D");
			t1->Branch("accp_large",&accp_large,"data/D");
			/*}}}*/

			cerr<<Form("--- Working on zflag=%d (min=%3.2f,max=%3.2f), Q2flag=%d (min=%3.2f,max=%3.2f) ...",
					zflag,zmin,zmax,Q2flag,Q2min,Q2max)<<endl;
			for (Int_t i=0;i!=T->GetEntries();i++){
				T->GetEntry(i);
				if (z>=zmin&&z<zmax&&Q2>=Q2min&&Q2<Q2max){

					/*Acceptance Weight{{{*/
					int ele_theta_bin=int(theta_ele/3.14159265*180.0/0.2)+1;    //0.2 degree per bin
					int ele_p_bin=int(mom_ele/0.05)+1;      //0.05 GeV per bin for mom
					double ele_forward_acceptance=accep_ele_forward->GetBinContent(ele_theta_bin,ele_p_bin);
					double ele_large_acceptance=accep_ele_large->GetBinContent(ele_theta_bin,ele_p_bin);

					if(theta_ele/3.1415926*180.<8.0)//A hard angle cut due to the HGS acceptance
                       ele_forward_acceptance = 0.0;

					int hadron_theta_bin=int(theta_had/3.14159265*180.0/0.2)+1;  //0.2 degree per bin
					int hadron_p_bin=int(mom_had/0.05)+1;     //0.05 GeV per b_in
					double hadron_acceptance=accep_had->GetBinContent(hadron_theta_bin,hadron_p_bin);

					if(theta_had/3.1415926*180.<8.0)//A hard angle cut due to the HGS acceptance
                       hadron_acceptance = 0.0;

					double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
					double large_acceptance=ele_large_acceptance*hadron_acceptance;
                    /*}}}*/
		
					ele_energy = nu+mom_ele;

					/*Get Luminosity and Beam Time{{{*/
					//1e36cm^-1*s^-1 for he3, 1e-33 is for nbar->cm
					double luminosity = 1e36 * 1e-33;
					double day = 0; // days
					// nevents = dxs (nbar) * L (nucleons/cm^2/s) * T(days) 
					// 1 day = 24hr*3600s = 86400
					// nbar=10-9 barn = 10^-9*10^-28 m^2=10^-9*10^-28 *10^4 cm^2=10^-33 cm^2
					// so weight  = Lumi * nbar * 86400 *day *acc_ele*acc_had / nsim;
					if(fabs(ele_energy-11.0)<0.1)
						day  = 48; // days
					if(fabs(ele_energy-8.80)<0.1)
						day  = 21; // days
					double time = day * 24 *3600; //sec
					/*}}}*/
					
				   	//weight_hp = dxs_hp * ele_accep *had_accep/Nsim;
					if(particle=="pip"||particle=="kp"){
						weight = weight_hp * luminosity * time * forward_acceptance;
						if(mom_ele>3.5)
							weight += weight_hp * luminosity * time * large_acceptance;
					}
					if(particle=="pim"||particle=="km"){
						weight = weight_hm * luminosity * time * forward_acceptance;
						if(mom_ele>3.5)
							weight += weight_hm * luminosity * time * large_acceptance;
					}
					accp_forward = forward_acceptance;
					accp_large   = large_acceptance;
					t1->Fill();
					if(!(i%10000))
						cerr<<Form("--- Working on zflag=%d, Q2flag=%d, evt=%d",zflag,Q2flag,i)<<"\r";
				}
			}
			file->Write();
			file->Close();
		}
	}

	delete T;
	return 0;
}
