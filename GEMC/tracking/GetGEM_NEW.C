//////////////////////////////////////////////////////////
// This script works on GEMC 2.0 only, because
// GEMC2.0 and GEMC1.0 have different TTree and branch names
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

void GetGEM(string input_filename)
//Int_t main()
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const Double_t DEG=180./3.1415926;

	/*Tree and branches{{{*/

	TChain *header = new TChain("header");
	TChain *generated = new TChain("generated");
	TChain *flux = new TChain("flux");
	header->Add(input_filename.c_str());
	generated->Add(input_filename.c_str());
	flux->Add(input_filename.c_str());

	//Header Tree:
	// Var#1~#8 are free slots for propogating important info from the "INPUT generator seed"
	// For example, they can be used to store the cross section and other physics quantities
	// For SIDIS, we agree to store the following quantities:
	// var1->x, var2->z, var3->pt, var4->Q2,var5->W, var6->cdxs, var7->phi_s, var8->phi_h 
	//
	//TTree *header = (TTree*) file->Get("header");
	vector <Double_t> *head_evn=0,*head_evn_type=0; //Note: Vectors have to be initialized at first!!!
	vector <Double_t> *head_beamPol=0;
	vector<Double_t> *head_x=0, *head_z=0, *head_pt=0, *head_Q2=0, *head_W=0, *head_cdxs=0, *head_phis=0, *head_phih=0;
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

	//TTree *generated = (TTree*) file->Get("generated");
	vector <Int_t> *gen_pid=0;
	vector <Double_t> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	generated->SetBranchAddress("pid",&gen_pid);
	generated->SetBranchAddress("px",&gen_px);
	generated->SetBranchAddress("py",&gen_py);
	generated->SetBranchAddress("pz",&gen_pz);
	generated->SetBranchAddress("vx",&gen_vx);
	generated->SetBranchAddress("vy",&gen_vy);
	generated->SetBranchAddress("vz",&gen_vz);

	//TTree *flux = (TTree*) file->Get("flux");
	vector<Double_t> *flux_id=0,*flux_hitn=0;
	vector<Double_t> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<Double_t> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
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

	Int_t nevent = (Int_t)generated->GetEntries();
	cout << "nevent " << nevent << endl;

	/*End Set Branch}}}*/

	/*Define GEM{{{*/
	//In each segmentation, 40 slides, each slides has 2.5cm width, 0.3cm between two slides.
	const unsigned int GEM_Layer = 23; //30 module around the circle
	const unsigned int GEM_Plane = 6;// put 6 slides in each module just for check the R-dependence, 
	const Double_t GEM_R_Min[GEM_Plane] = {36,21,25,32,42,55};//cm, GEM inner radius 
	const Double_t GEM_R_Max[GEM_Plane] = {87,98,112,135,100,123};//cm, GEM outer radius
	const Double_t GEM_Z_Set[GEM_Plane] = {-175,-150,-119,-68,5,92}; //cm, GEM Z position
	const Int_t GEM_ID[GEM_Plane] = {1100000,1200000,1300000,1400000,1500000,1600000}; //GEM virtual plane
	//const Int_t GEM_ID[GEM_Plane] = {1120000,1220000,1320000,1420000,1520000,1620000}; //GEM virtual plane
	const Double_t GEM_VP_Z[GEM_Plane] = {-175-0.5,-150-0.5,-119-0.5,-68-0.5,5-0.5,92-0.5}; //GEM virtual plane
	const Int_t GEM_VP_ID[GEM_Plane] = {1110000,1210000,1310000,1410000,1510000,1610000,}; //GEM virtual plane

	Double_t GEM_X[GEM_Plane], GEM_Y[GEM_Plane], GEM_Z[GEM_Plane];
	Double_t GEM_Vx[GEM_Plane], GEM_Vy[GEM_Plane], GEM_Vz[GEM_Plane];
	Double_t GEM_Px[GEM_Plane], GEM_Py[GEM_Plane], GEM_Pz[GEM_Plane];
	Double_t GEM_Mom[GEM_Plane], GEM_Theta[GEM_Plane], GEM_Phi[GEM_Plane];
	Double_t GEM_E[GEM_Plane], GEM_Edep[GEM_Plane];
	Int_t GEM_ID_Pick[GEM_Plane], GEM_PID_Pick[GEM_Plane];
	//I just want to initialize everything
	for(unsigned int ip=0;ip<GEM_Plane;ip++){
		GEM_X[ip] = 0.0; GEM_Y[ip] = 0.0; GEM_Z[ip] = 0.0;
		GEM_Vx[ip] = 0.0; GEM_Vy[ip] = 0.0; GEM_Vz[ip] = 0.0;
		GEM_Px[ip] = 0.0; GEM_Py[ip] = 0.0; GEM_Pz[ip] = 0.0;
		GEM_Mom[ip] = 0.0; GEM_Theta[ip] = 0.0; GEM_Phi[ip] = 0.0;
		GEM_E[ip] = 0.0; GEM_Edep[ip]=0.;
	}

	//Note: Use EC's Virtual Plane to determine the positions and momentum info on EC
	const Int_t LAEC_VP_ID = 3210000;
	const Double_t LAEC_VP_Rmin = 83.0, LAEC_VP_Rmax = 140.0, LAEC_VP_Z = -65.0;//cm
	const Int_t FAEC_VP_ID = 3110000;
	const Double_t FAEC_VP_Rmin = 105.0, FAEC_VP_Rmax = 235.0, FAEC_VP_Z = 427.5;//cm

	//Other Definition	
	const Int_t Electron = 11;
	const Int_t Gamma = 22;
	const Int_t Beam = 0;
	const Int_t Neutron = 2112;
	const Int_t Neutrino1 = 12;//Nu_e
	const Int_t Neutrino2 = 14;//Nu_Mu
	const Int_t Neutrino3 = 16;//Nu_Tao

	const Double_t MM2CM = 1/10.;
	const Double_t MeV2GeV = 1/1000.;
	/*}}}*/

	/*Read in each event{{{*/
    TString pid_temp = "";
	/*Find PID{{{*/
	if(input_filename.find("_pi0_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pi0_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pi0_down";
		else
			pid_temp = "pi0";
	}
	else if(input_filename.find("_pip_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pip_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pip_down";
		else
			pid_temp = "pip";
	}
	else if(input_filename.find("_pim_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pim_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pim_down";
		else
			pid_temp = "pim";
	}
	else if(input_filename.find("_kp_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "kp_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "kp_down";
		else
			pid_temp = "kp";
	}
	else if(input_filename.find("_km_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "km_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "km_down";
		else
			pid_temp = "km";
	}
	else if(input_filename.find("_p_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "p_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "p_down";
		else
			pid_temp = "p";
	}
	else if(input_filename.find("_eDIS_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "eDIS_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "eDIS_down";
		else	
			pid_temp = "eDIS";
	}
	else if(input_filename.find("_EM_",0) != string::npos) {
		pid_temp = "EM";
	}
	else
		cerr<<"****ERROR, I don't understand the file name!!! *****"<<endl;
     /*}}}*/

	ofstream laoutfile(Form("SoLID_SIDIS_He3_LA_%s.dat",pid_temp.Data()));
	ofstream faoutfile(Form("SoLID_SIDIS_He3_FA_%s.dat",pid_temp.Data()));


	Int_t nselected = nevent;
	Int_t LA_evn = 0, FA_evn = 0;
	cerr<<"++++++++++++++++ "<<endl;
	for(Long64_t i=0;i<nselected;i++){
		cout<<i<<"\r";
		//From Header, head_evn, usally only one particle
		header->GetEntry(i);
        //Debug
//		cerr<<Form("--- Evn=%d var1=%f var2=%f var3=%f var4=%f var5=%f",
//			   (Int_t)head_evn->at(0),head_x->at(0),head_z->at(0),head_pt->at(0),head_Q2->at(0),head_W->at(0))<<endl;

        //From Generated, gen_pid,gen_vx/_vy/_vz,gen_px/_py/_pz
		generated->GetEntry(i);

        /*Fill Generated Info {{{*/		
		Double_t gen_mom[gen_pid->size()], gen_theta[gen_pid->size()], gen_phi[gen_pid->size()];
		Int_t Is_Electron = -1000;
		//Depends on how many particles, for SIDIS, there might be two particles (e and pi/K)
		for(unsigned int i_gen=0;i_gen<gen_pid->size();i_gen++){
            //Debug
//			cerr<<Form("--- pid=%d px=%e py=%e pz=%e vx=%e vy=%e vz=%e", 
//					gen_pid->at(i_gen),gen_px->at(i_gen),gen_py->at(i_gen),gen_pz->at(i_gen),
//					gen_vx->at(i_gen),gen_vy->at(i_gen),gen_vz->at(i_gen))<<endl;

			gen_mom[i_gen] = 0.0; gen_theta[i_gen]=0.0; gen_phi[i_gen]=0.0;
			if(gen_pid->at(i_gen)!=Electron) continue;
			Is_Electron = i_gen;
			gen_mom[i_gen] = sqrt(pow(gen_px->at(i_gen),2)
					+pow(gen_py->at(i_gen),2)
					+pow(gen_pz->at(i_gen),2));
			gen_theta[i_gen] = acos( gen_pz->at(i_gen)/gen_mom[i_gen] )*DEG;
			gen_phi[i_gen] = atan2(gen_py->at(i_gen),gen_px->at(i_gen))*DEG;
		}
        if(Is_Electron<0) continue;
        /*}}}End Fill Generated info*/

		//From Flux, flux_avg_x/_y/_z/_vx/_vy/_vz/_px/_py/_pz
		//           flux_hitn,flux_id/_pid/_mpid/_mvx/_mvy/_mvz,
		flux->GetEntry(i);

        /*Fill ECs from Flux{{{*/
		Double_t LAEC_X=0.0,  LAEC_Y=0.0,  LAEC_Z=0.0;
        Double_t LAEC_Vx=0.0, LAEC_Vy=0.0, LAEC_Vz=0.0;
        Double_t LAEC_Px=0.0, LAEC_Py=0.0, LAEC_Pz=0.0;
        Double_t FAEC_X=0.0,  FAEC_Y=0.0,  FAEC_Z=0.0;
		Double_t FAEC_Vx=0.0, FAEC_Vy=0.0, FAEC_Vz=0.0;
		Double_t FAEC_Px=0.0, FAEC_Py=0.0, FAEC_Pz=0.0;
		bool Is_LAEC_Hit = kFALSE, Is_FAEC_Hit = kFALSE;	
        Int_t LAEC_PID = -10;  Int_t FAEC_PID = -10; Int_t EC_ID = -10000;
		for (unsigned int j=0;j<flux_id->size();j++) {
            //Debug
//			cerr<<Form("---Evn=%d FID=%d ID=%d PID=%d Gen_Mom=%e E=%e Px=%e Py=%e Pz=%e Vx=%e Vy=%e Vz=%e",
//					(Int_t)head_evn->at(0),j,(Int_t)flux_id->at(j),(Int_t)flux_pid->at(j), gen_mom[Is_Electron],flux_trackE->at(j),
//					flux_px->at(j),flux_py->at(j),flux_pz->at(j),flux_vx->at(j),flux_vy->at(j),flux_vz->at(j)	
//					)<<endl;
			
			EC_ID = (Int_t) flux_id->at(j);
			//if(flux_pz->at(j)<-1e-19) continue; //Discard backward particles
			//if(fabs(flux_pid->at(j)-Electron)>0.5) continue; // Only check electrons

			Double_t EC_Mom=sqrt(pow(flux_px->at(j),2)+pow(flux_py->at(j),2)+pow(flux_pz->at(j),2));//MeV
			//if((gen_mom[Is_Electron]-EC_Mom)<0||(gen_mom[Is_Electron]-EC_Mom)>100) continue;//Discard particles with 100MeV less energy
			//if(EC_Mom<1e+1) continue; //Discard low energy particles
			//cerr<<Form("---  ID = %d, PID=%d, Gen_Mom=%e, EC_Mom=%e ",EC_ID,EC_PID, gen_mom[Is_Electron], EC_Mom)<<endl;
		
			Double_t EC_R = sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2))*MM2CM;//cm
			if(EC_ID == LAEC_VP_ID && (EC_R>LAEC_VP_Rmin && EC_R<LAEC_VP_Rmax)){
				Is_LAEC_Hit = kTRUE;
				LAEC_PID = (Int_t) (flux_pid->at(j));
				LAEC_X = flux_avg_x->at(j)*MM2CM;
			    LAEC_Y = flux_avg_y->at(j)*MM2CM;
			    LAEC_Z = flux_avg_z->at(j)*MM2CM;
			    LAEC_Vx = flux_vx->at(j)*MM2CM;
			    LAEC_Vy = flux_vy->at(j)*MM2CM;
			    LAEC_Vz = flux_vz->at(j)*MM2CM;
			    LAEC_Px = flux_px->at(j)*MeV2GeV;
			    LAEC_Py = flux_py->at(j)*MeV2GeV;
			    LAEC_Pz = flux_pz->at(j)*MeV2GeV;
			//	cerr<<"--- There is a LAEC Hit "<<", ID = "<<EC_ID<<endl;
			}
			else if(EC_ID == FAEC_VP_ID && (EC_R>FAEC_VP_Rmin && EC_R<FAEC_VP_Rmax)){
				Is_FAEC_Hit = kTRUE;
				FAEC_PID = (Int_t) (flux_pid->at(j));
				FAEC_X = flux_avg_x->at(j)*MM2CM;
			    FAEC_Y = flux_avg_y->at(j)*MM2CM;
			    FAEC_Z = flux_avg_z->at(j)*MM2CM;
			    FAEC_Vx = flux_vx->at(j)*MM2CM;
			    FAEC_Vy = flux_vy->at(j)*MM2CM;
			    FAEC_Vz = flux_vz->at(j)*MM2CM;
			    FAEC_Px = flux_px->at(j)*MeV2GeV;
			    FAEC_Py = flux_py->at(j)*MeV2GeV;
			    FAEC_Pz = flux_pz->at(j)*MeV2GeV;
			//	cerr<<"--- There is a FAEC Hit "<<", ID = "<<EC_ID<<endl;
			}
		}
       /*}}}End Fill EC*/
		
		if(!(Is_FAEC_Hit) && !(Is_LAEC_Hit)) continue;
                 
		/*Fill GEMs from Flux{{{*/
		//Count how many GEMs have been hit
		//Int_t Count_GEM_Hit_LA = 0; // 4 GEMs before LA, so at lease 3 planes needed to have hits to make a LA tracks 
		Int_t Count_GEM_Hit = 0; //Another 2GEMs after LA and before FA, so at least 2 from LA and 1 from FA
		bool Is_GEM_Hit[GEM_Plane];
		for( unsigned int ip=0;ip<GEM_Plane;ip++){
			GEM_X[ip] = 0.0;   GEM_Y[ip] = 0.0;     GEM_Z[ip] = 0.0;
			GEM_Vx[ip] = 0.0;  GEM_Vy[ip] = 0.0;    GEM_Vz[ip] = 0.0;
			GEM_Px[ip] = 0.0;  GEM_Py[ip] = 0.0;    GEM_Pz[ip] = 0.0;
			GEM_Mom[ip] = 0.0; GEM_Theta[ip] = 0.0; GEM_Phi[ip] = 0.0;
            GEM_PID_Pick[ip] = -10;  GEM_ID_Pick[ip] = -99999999; 
			GEM_E[ip] = 0.0, GEM_Edep[ip]=0.0;
			Is_GEM_Hit[ip] = kFALSE;

			for (unsigned int j=0;j<flux_id->size();j++) {
				Int_t ID_Pick= (Int_t) (flux_id->at(j));
				if(ID_Pick!=GEM_VP_ID[ip]) continue; //Only check particles in one GEM layer
				Int_t PID_Pick= (Int_t) (flux_pid->at(j));
				if(PID_Pick!=Electron) continue; //Only check particles in one GEM layer

				Double_t r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2))*MM2CM;//cm
				if(r>GEM_R_Max[ip] || r<GEM_R_Min[ip]) continue;//Only inside the GEM Disk
                 
				//if(flux_pz->at(j)<-1e-19) continue; //Discard backward particles
				Double_t fmom=sqrt(pow(flux_px->at(j),2)+pow(flux_py->at(j),2)+pow(flux_pz->at(j),2));//MeV
				if((gen_mom[Is_Electron]-fmom)<0||(gen_mom[Is_Electron]-fmom)>100) continue;//Discard particles with 100MeV less energy
				//if(fmom<1e+1) continue; //Discard low energy particles
				
				Double_t flux_theta = acos(flux_pz->at(j)/fmom)*DEG;
				Double_t flux_phi = atan2(flux_py->at(j),flux_px->at(j))*DEG;
				
				GEM_ID_Pick[ip] = ID_Pick;
				GEM_PID_Pick[ip] = PID_Pick;
				GEM_X[ip] = flux_avg_x->at(j)*MM2CM; 
				GEM_Y[ip] = flux_avg_y->at(j)*MM2CM; 
				GEM_Z[ip] = flux_avg_z->at(j)*MM2CM;
				GEM_Vx[ip] = flux_vx->at(j)*MM2CM;   
				GEM_Vy[ip] = flux_vy->at(j)*MM2CM;   
				GEM_Vz[ip] = flux_vz->at(j)*MM2CM;
				GEM_Px[ip] = flux_px->at(j)*MeV2GeV;   
				GEM_Py[ip] = flux_py->at(j)*MeV2GeV;   
				GEM_Pz[ip] = flux_pz->at(j)*MeV2GeV;
				GEM_Mom[ip] = fmom*MeV2GeV; 
				GEM_Theta[ip] = flux_theta; 
				GEM_Phi[ip] = flux_phi;
				GEM_E[ip] = flux_trackE->at(j)*MeV2GeV;
				GEM_Edep[ip] = flux_totEdep->at(j)*MeV2GeV;
                Is_GEM_Hit[ip] = kTRUE;

				Count_GEM_Hit++;
		}//for (unsigned int j=0;j<flux_id->size();j++) 
		}//for( unsigned int ip=0;ip<GEM_Plane;ip++)
        /*}}}End Fill GEMCs*/

	    /*PrInt_t Out{{{*/
		
		//if(Is_LAEC_Hit&&Count_GEM_Hit_LA>=3 && !(Is_FAEC_Hit)){
		if(Is_LAEC_Hit && Count_GEM_Hit >=3 && !(Is_GEM_Hit[4]) && !(Is_GEM_Hit[5]) ){
            LA_evn ++;
			laoutfile <<Form("%5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
					LA_evn, gen_pid->at(Is_Electron),
//					gen_mom[Is_Electron]*MeV2GeV,gen_theta[Is_Electron],gen_phi[Is_Electron],
					gen_px->at(Is_Electron)*MeV2GeV,
					gen_py->at(Is_Electron)*MeV2GeV,
					gen_pz->at(Is_Electron)*MeV2GeV,
					gen_vx->at(Is_Electron)*MM2CM,
					gen_vy->at(Is_Electron)*MM2CM,
					gen_vz->at(Is_Electron)*MM2CM)
					<<endl;
			for(Int_t ip=0;ip<4;ip++){//four GEMCs in front of LAEC
				laoutfile <<Form("%2d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
						ip+1, GEM_PID_Pick[ip],
//						GEM_Mom[ip],GEM_Theta[ip],GEM_Phi[ip],
						GEM_Px[ip], GEM_Py[ip], GEM_Pz[ip],
						GEM_Vx[ip], GEM_Vy[ip], GEM_Vz[ip],
						GEM_X[ip], GEM_Y[ip], GEM_Z[ip],GEM_E[ip], GEM_Edep[ip])
					<<endl;
			}
			laoutfile <<Form("%2d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
			        310, LAEC_PID,
					LAEC_Px, LAEC_Py, LAEC_Pz,
					LAEC_Vx, LAEC_Vy, LAEC_Vz,
					LAEC_X,  LAEC_Y,  LAEC_Z)
					<<endl;
			}
		if(Is_FAEC_Hit&&Count_GEM_Hit>=5){
            FA_evn ++;
			faoutfile <<Form("%5d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
					FA_evn, gen_pid->at(Is_Electron),
//					gen_mom[Is_Electron]*MeV2GeV, gen_theta[Is_Electron], gen_phi[Is_Electron],
					gen_px->at(Is_Electron)*MeV2GeV,
					gen_py->at(Is_Electron)*MeV2GeV,
					gen_pz->at(Is_Electron)*MeV2GeV,
					gen_vx->at(Is_Electron)*MM2CM,
					gen_vy->at(Is_Electron)*MM2CM,
					gen_vz->at(Is_Electron)*MM2CM)
					<<endl;
			for(Int_t ip=0;ip<6;ip++){//four GEMCs in front of LAEC
				faoutfile <<Form("%2d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
						ip+1, GEM_PID_Pick[ip],
//						GEM_Mom[ip],GEM_Theta[ip],GEM_Phi[ip],
						GEM_Px[ip], GEM_Py[ip], GEM_Pz[ip],
						GEM_Vx[ip], GEM_Vy[ip], GEM_Vz[ip],
						GEM_X[ip], GEM_Y[ip], GEM_Z[ip],GEM_E[ip], GEM_Edep[ip])
					<<endl;
			}
			faoutfile <<Form("%2d %5d %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e %10.2e", 
			        320, FAEC_PID,
					FAEC_Px, FAEC_Py, FAEC_Pz,
					FAEC_Vx, FAEC_Vy, FAEC_Vz,
					FAEC_X, FAEC_Y, FAEC_Z)
					<<endl;
			}
	}
    /*}}}End PrInt_t Out*/

	/*End Read in each event}}}*/

//	file->Close();
    header->Delete();
	generated->Delete();
	flux->Delete();
	laoutfile.close();
	faoutfile.close();
}
