//////////////////////////////////////////////////////////
// This script works on SPDC 2.0 only, because
// SPDC2.0 and SPDC1.0 have different TTree and branch names
// Zhihong Ye, 09/12/2014
//
#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <vector>
#include <TCanvas.h>
#include <TFile.h>
#include <TLine.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TH3F.h>
#include <TF1.h>
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

using namespace std;

double GetSPD(TString input_filename)
	//int main()
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const double DEG=180./3.1415926;

	//	TString input_filename; cerr<<"--- File = "; cin >> input_filename;

	TFile *file=new TFile(input_filename.Data());
	if (file->IsZombie()) {
		cout << "Error opening file" << input_filename << endl;
		exit(-1);
	}
	else cout << "open file " << input_filename << endl;

	/*Set Branch{{{*/
	//Header Tree:
	// Var#1~#8 are free slots for propogating important info from the "INPUT generator seed"
	// For example, they can be used to store the cross section and other physics quantities
	// For SIDIS, we agree to store the following quantities:
	// var1->x, var2->z, var3->pt, var4->Q2,var5->W, var6->cdxs, var7->phi_s, var8->phi_h 
	//
	TTree *header = (TTree*) file->Get("header");
	vector <double> *head_evn=0,*head_evn_type=0; //Note: Vectors have to be initialized at first!!!
	vector <double> *head_beamPol=0;
	vector<double> *head_x=0, *head_z=0, *head_pt=0, *head_Q2=0, *head_W=0, *head_cdxs=0, *head_phis=0, *head_phih=0;
	header->SetBranchAddress("evn",&head_evn);
	header->SetBranchAddress("evn_type",&head_evn_type);
	header->SetBranchAddress("beamPol",&head_beamPol);
	header->SetBranchAddress("var1",    &head_x);
	header->SetBranchAddress("var2",    &head_z);
	header->SetBranchAddress("var3",    &head_pt);
	header->SetBranchAddress("var4",    &head_Q2);
	header->SetBranchAddress("var5",    &head_W);
	header->SetBranchAddress("var6",    &head_cdxs);
	header->SetBranchAddress("var7",    &head_phis);
	header->SetBranchAddress("var8",    &head_phih);

	TTree *generated = (TTree*) file->Get("generated");
	vector <int> *gen_pid=0;
	vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	generated->SetBranchAddress("pid",&gen_pid);
	generated->SetBranchAddress("px",&gen_px);
	generated->SetBranchAddress("py",&gen_py);
	generated->SetBranchAddress("pz",&gen_pz);
	generated->SetBranchAddress("vx",&gen_vx);
	generated->SetBranchAddress("vy",&gen_vy);
	generated->SetBranchAddress("vz",&gen_vz);

	TTree *flux = (TTree*) file->Get("flux");
	vector<double> *flux_id=0,*flux_hitn=0;
	vector<double> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
	flux->SetBranchAddress("hitn",&flux_hitn);
	flux->SetBranchAddress("id",&flux_id);
	flux->SetBranchAddress("pid",&flux_pid);
	flux->SetBranchAddress("mpid",&flux_mpid);
	flux->SetBranchAddress("tid",&flux_tid);
	flux->SetBranchAddress("mtid",&flux_mtid);
	flux->SetBranchAddress("otid",&flux_otid);
	flux->SetBranchAddress("trackE",&flux_trackE);
	flux->SetBranchAddress("totEdep",&flux_totEdep);
	flux->SetBranchAddress("avg_x",&flux_avg_x);
	flux->SetBranchAddress("avg_y",&flux_avg_y);
	flux->SetBranchAddress("avg_z",&flux_avg_z);
	flux->SetBranchAddress("avg_lx",&flux_avg_lx);
	flux->SetBranchAddress("avg_ly",&flux_avg_ly);
	flux->SetBranchAddress("avg_lz",&flux_avg_lz);
	flux->SetBranchAddress("px",&flux_px);
	flux->SetBranchAddress("py",&flux_py);
	flux->SetBranchAddress("pz",&flux_pz);
	flux->SetBranchAddress("vx",&flux_vx);
	flux->SetBranchAddress("vy",&flux_vy);
	flux->SetBranchAddress("vz",&flux_vz);
	flux->SetBranchAddress("mvx",&flux_mvx);
	flux->SetBranchAddress("mvy",&flux_mvy);
	flux->SetBranchAddress("mvz",&flux_mvz);
	flux->SetBranchAddress("avg_t",&flux_avg_t);

	int nevent = (int)generated->GetEntries();
	cout << "nevent " << nevent << endl;

	/*End Set Branch}}}*/

	/*Define SPD{{{*/
	const double LASPD_RMin = 80.;//cm, SPD inner radius 
	const double LASPD_RMax = 135.;//cm, SPD outer radius
	const int LASPD_ID = 5200000; //SPD front virtual plane
	const double LASPD_Z_Set = -67.; //cm, SPD Z position
	const int LASPD_VP_Front_ID = 5210000; //SPD front virtual plane
	const double LASPD_VP_Front_Z = -67.3; //SPD virtual plane
	const int LASPD_VP_Rear_ID = 3210000; //SPD rear virtual plane,is LAEC front
	const double LASPD_VP_Rear_Z = -66.5; //SPD virtual plane

	const double FASPD_RMin = 96.;//cm, SPD inner radius 
	const double FASPD_RMax = 210.;//cm, SPD outer radius
	const int FASPD_ID = 5100000; //SPD front virtual plane
	const double FASPD_Z_Set = 407.; //cm, SPD Z position
	const int FASPD_VP_Front_ID = 5110000; //SPD front virtual plane
	const double FASPD_VP_Front_Z =  406.5; //SPD virtual plane
	const int FASPD_VP_Rear_ID = 4110000; //SPD rear virtual plane,is mrpc VP front
	const double FASPD_V_RearP_Z = 408.0; //SPD virtual plane
	/*}}}End Defined SPD*/

	/*Define EC{{{*/
	//Note: Use EC's Virtual Plane to determine the positions and momentum info on EC
	const int LAEC_ID = 3200000;//Inactive currently
	const double LAEC_RMin = 83.0;
	const double LAEC_RMax = 140.0;
	const double LAEC_Z_Set = 0.0;

	const int LAEC_VP_Front_ID = 3210000;
	const double LAEC_VP_Front_Rmin = 75.0, LAEC_VP_Front_Rmax = 144.0, LAEC_VP_Front_Z = -66.5;//cm
	const int LAEC_VP_Mid_ID = 3220000;
	const double LAEC_VP_Mid_Rmin = 75.0, LAEC_VP_Mid_Rmax = 144.0, LAEC_VP_Mid_Z = -65.4;//cm
	const int LAEC_VP_Inner_ID = 3230000;
	const double LAEC_VP_Inner_Rmin = 75.0, LAEC_VP_Inner_Rmax = 88.0, LAEC_VP_Inner_Z = -40.0;//cm
	const int LAEC_VP_Rear_ID = 3240000;
	const double LAEC_VP_Rear_Rmin = 89.0, LAEC_VP_Rear_Rmax = 144.0, LAEC_VP_Rear_Z = -14.0;//cm

	const int FAEC_ID = 3100000;//Inactive currently
	const double FAEC_RMin = 105.0;
	const double FAEC_RMax = 235.0;
	const double FAEC_Z_Set = 440.0;
	const double FAEC_Dz = 25.0;

	const int FAEC_VP_Front_ID = 3110000;
	const double FAEC_VP_Front_Rmin = 90.0, FAEC_VP_Front_Rmax = 265.0, FAEC_VP_Front_Z = 413.0;//cm
	const int FAEC_VP_Mid_ID = 3120000;
	const double FAEC_VP_Mid_Rmin = 90.0, FAEC_VP_Mid_Rmax = 265.0, FAEC_VP_Mid_Z = 414.6;//cm
	const int FAEC_VP_Inner_ID = 3130000;
	const double FAEC_VP_Inner_Rmin = 96.0, FAEC_VP_Inner_Rmax = 96.1, FAEC_VP_Inner_Z = 440.0;//cm
	const int FAEC_VP_Rear_ID = 3140000;
	const double FAEC_VP_Rear_Rmin = 90.0, FAEC_VP_Rear_Rmax = 265.0, FAEC_VP_Rear_Z = 466.0;//cm
	/*}}}End Defined EC*/

	/*Other Definition{{{*/	
	const int Electron = 11;
	const int Positron =-11;
	const int Photon = 22;
	//	const int Pip = 211;
	//	const int Pim =-211;
	//	const int Beam = 0;
	//	const int Neutron = 2112;
	//	const int Neutrino1 = 12;//Nu_e
	//	const int Neutrino2 = 14;//Nu_Mu
	//	const int Neutrino3 = 16;//Nu_Tao

	const double MM2CM = 1/10.;
	const double MeV2GeV = 1/1000.;
	/*}}}*/

	double Edep = 0.0, Ee= 0.0, Eg= 0.0;
	int PID = 0, Is_Hit=0, Is_Fire=0;

	TString output_filename = input_filename;
	output_filename.ReplaceAll(".root","_out.root");
	TFile *outroot = new TFile(output_filename.Data(),"recreate");
	TTree* T = new TTree("T","T");
	T->Branch("Ee",&Ee,"Ee/D");
	T->Branch("Eg",&Eg,"Eg/D");
	T->Branch("Edep",&Edep,"Edep/D");
	T->Branch("PID",&PID,"PID/I");
	T->Branch("Is_Hit",&Is_Hit,"Is_Hit/I");
	T->Branch("Is_Fire",&Is_Fire,"Is_Fire/I");

	/*Read in each event{{{*/
	Int_t nselected = nevent;
	double Edep_total = 0.0;
	cerr<<"++++++++++++++++ "<<endl;
	int Count_Hit = 0, Count_Fire = 0;
	for(Long64_t i=0;i<nselected;i++){
		cout<<i<<"\r";
		Edep = 0.0; Ee = 0.0; Eg = 0.0;
		Is_Hit=0;
		Is_Fire=0;

		//From Header, head_evn, usally only one particle
		header->GetEntry(i);
		//From Generated, gen_pid,gen_vx/_vy/_vz,gen_px/_py/_pz
		generated->GetEntry(i);

		/*Generated Info {{{*/		
		//Depends on how many particles, for SIDIS, there might be two particles (e and pi/K)
		double gen_mom[gen_pid->size()],gen_theta[gen_pid->size()],gen_phi[gen_pid->size()];
		int Is_Electron=-1, Is_Photon=-1;
		for(unsigned int ig=0;ig<gen_pid->size();ig++){
			int pid = gen_pid->at(ig);

			//Note: in GEMC, momentrum are in MeV, length in mm
			gen_mom[ig] = sqrt( pow(gen_px->at(ig),2)+pow(gen_py->at(ig),2)+pow(gen_pz->at(ig),2) );
			gen_theta[ig] = acos(gen_pz->at(ig)/gen_mom[ig]) * DEG;
			gen_phi[ig] = atan2( gen_py->at(ig), gen_px->at(ig))*DEG; 
			if(pid==Electron){ Is_Electron=ig;
				Ee = gen_mom[Is_Electron];
				Eg = -1000;
				PID = pid;
			}
			else if(pid==Photon){
				Is_Photon=ig;
				Eg = gen_mom[Is_Photon];
				Ee = -1000;
				PID = pid;
			}
		}
		if(Is_Electron<0) continue; //Only see events from electrons
		//if(Is_Photon<0) continue; //Only see events from photons
		/*}}}End Generated info*/

		//From Flux, flux_avg_x/_y/_z/_vx/_vy/_vz/_px/_py/_pz
		//           flux_hitn,flux_id/_pid/_mpid/_mvx/_mvy/_mvz,
		flux->GetEntry(i);

		/*SPDs from Flux{{{*/
		for (unsigned int j=0;j<flux_id->size();j++){
			if((int)(flux_id->at(j)!=FASPD_ID)) continue;//Only look at particle in this detector plane

			double r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2))*MM2CM;//cm
			if(r>FASPD_RMax|| r<FASPD_RMin) continue;//Only inside the SPD Disk

			if(flux_pz->at(j)<-1e-19) continue; //Discard backward particles
			double fmom=sqrt(pow(flux_px->at(j),2)+pow(flux_py->at(j),2)+pow(flux_pz->at(j),2));//MeV
			//if(fmom<1e+1) continue; //Discard low energy particles

			Is_Hit = 1;//Only count particles that pass through the detector
			int PID_flux = (int) (flux_pid->at(j));
			//if((int)flux_pid->at(j)!=Photon) continue; //Only check photons in SPD
			if(abs(PID_flux)!=Electron) continue;//Only check electrons in SPD
			Is_Fire = 1;//Treat the SPD fired when at lease one electron passing through

			//double fE = flux_trackE->at(j)*MeV2GeV;//GeV
			//double flux_theta = acos(flux_pz->at(j)/fmom)*DEG;
			//double flux_phi = atan2(flux_py->at(j),flux_px->at(j))*DEG;
			double fEdep = flux_totEdep->at(j);//MeV
			Edep += fEdep;//MeV
		}//for (unsigned int j=0;j<flux_id->size();j++) 
		/*}}}End SPDCs*/

		if(Is_Hit)		
			Count_Hit ++;//All particles that pass through the detectors

		if(Is_Fire)
			Count_Fire ++;

		Edep_total += Edep;

		T->Fill();
	}
	/*End Read in each event}}}*/

	cerr<<"--- Done: "<<endl;
	cerr<<Form("---- For %s:  Hit = %d, Fire = %d, Edep_avg = %f ", input_filename.Data(), Count_Hit, Count_Fire, Edep_total/nevent)<<endl;

	outroot->cd(); T->Write(); outroot->Close();


	file->Close();
	return Edep_total/nevent;
}

void RunAll(){
 
	double E[21] = {0.3, 0.5, 0.7, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0, 500.0, 1000.0, 3000.0, 5000.0, 7000.0, 9000.0, 11000.0, 0.0, 0.0, 0.0};

	TString runlist;
	cerr << "--- Runlist = "; cin>> runlist;
	ifstream inputfile(runlist.Data());
	TString out_result = runlist;
	out_result.ReplaceAll(".txt","_result.dat");
	ofstream outputfile(out_result.Data());
    TString file_name;
	double Edep,Edep_err;
	int N = 0;
	while(!inputfile.eof()){
	    inputfile >> file_name;
	    Edep = GetSPD(file_name.Data());	
		//Edep_err = rate *( 1./sqrt(rate*10000));
		Edep_err = 0.0;
		outputfile << E[N] <<"      "<<Edep<<"    "<<Edep_err<<endl; 
		N++;
	}

	outputfile.close();
	inputfile.close();


}
