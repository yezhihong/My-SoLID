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
#include <TH3F.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void CheckEvent(string input_filename)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);

	/*Set Branch{{{*/
	TFile *file=new TFile(input_filename.c_str());
	if (file->IsZombie()) {
		cout << "Error opening file" << input_filename << endl;
		exit(-1);
	}
	else cout << "open file " << input_filename << endl;

	TTree *Tgen = (TTree*) file->Get("genT");
	Int_t gen_evn,gen_ngen;
	Int_t gen_id_array[1000];
	Int_t *gen_id=gen_id_array;
	Float_t gen_px_array[1000],gen_py_array[1000],gen_pz_array[1000],gen_p_array[1000],
			gen_phi_array[1000],gen_theta_array[1000],gen_vx_array[1000],gen_vy_array[1000],gen_vz_array[1000];
	Float_t *gen_px=gen_px_array,*gen_py=gen_py_array,*gen_pz=gen_pz_array,*gen_p=gen_p_array,
			*gen_phi=gen_phi_array,*gen_theta=gen_theta_array,*gen_vx=gen_vx_array,
			*gen_vy=gen_vy_array,*gen_vz=gen_vz_array;
	Tgen->SetBranchAddress("evn",&gen_evn);
	Tgen->SetBranchAddress("ngen",&gen_ngen);
	Tgen->SetBranchAddress("id",gen_id);
	Tgen->SetBranchAddress("px",gen_px);
	Tgen->SetBranchAddress("py",gen_py);
	Tgen->SetBranchAddress("pz",gen_pz);
	Tgen->SetBranchAddress("p",gen_p);
	Tgen->SetBranchAddress("phi",gen_phi);
	Tgen->SetBranchAddress("theta",gen_theta);
	Tgen->SetBranchAddress("vx",gen_vx);
	Tgen->SetBranchAddress("vy",gen_vy);
	Tgen->SetBranchAddress("vz",gen_vz);

	TTree *Tflux = (TTree*) file->Get("fluxT");
	Int_t flux_evn,flux_nfluxhit;
	Int_t flux_ID_array[1000],flux_pid_array[1000],flux_mpid_array[1000];
	Int_t *flux_ID=flux_ID_array,*flux_pid=flux_pid_array,*flux_mpid=flux_mpid_array;
	Float_t flux_Edep_array[1000],flux_E_array[1000],flux_x_array[1000],flux_y_array[1000],flux_z_array[1000],
			flux_lx_array[1000],flux_ly_array[1000],flux_lz_array[1000],flux_t_array[1000],flux_px_array[1000],
			flux_py_array[1000],flux_pz_array[1000],flux_vx_array[1000],flux_vy_array[1000],flux_vz_array[1000],
			flux_mvx_array[1000],flux_mvy_array[1000],flux_mvz_array[1000];
	Float_t *flux_Edep=flux_Edep_array,*flux_E=flux_E_array,*flux_x=flux_x_array,*flux_y=flux_y_array,
			*flux_z=flux_z_array,*flux_lx=flux_lx_array,*flux_ly=flux_ly_array,*flux_lz=flux_lz_array,
			*flux_t=flux_t_array,*flux_px=flux_px_array,*flux_py=flux_py_array,*flux_pz=flux_pz_array,
			*flux_vx=flux_vx_array,*flux_vy=flux_vy_array,*flux_vz=flux_vz_array,*flux_mvx=flux_mvx_array,
			*flux_mvy=flux_mvy_array,*flux_mvz=flux_mvz_array;
	Tflux->SetBranchAddress("evn",&flux_evn);
	Tflux->SetBranchAddress("nfluxhit",&flux_nfluxhit);
	Tflux->SetBranchAddress("ID",flux_ID);
	Tflux->SetBranchAddress("Edep",flux_Edep);
	Tflux->SetBranchAddress("E",flux_E);
	Tflux->SetBranchAddress("x",flux_x);
	Tflux->SetBranchAddress("y",flux_y);
	Tflux->SetBranchAddress("z",flux_z);
	Tflux->SetBranchAddress("lx",flux_lx);
	Tflux->SetBranchAddress("ly",flux_ly);
	Tflux->SetBranchAddress("lz",flux_lz);
	Tflux->SetBranchAddress("t",flux_t);
	Tflux->SetBranchAddress("pid",flux_pid);
	Tflux->SetBranchAddress("mpid",flux_mpid);
	Tflux->SetBranchAddress("px",flux_px);
	Tflux->SetBranchAddress("py",flux_py);
	Tflux->SetBranchAddress("pz",flux_pz);
	Tflux->SetBranchAddress("vx",flux_vx);
	Tflux->SetBranchAddress("vy",flux_vy);
	Tflux->SetBranchAddress("vz",flux_vz);
	Tflux->SetBranchAddress("mvx",flux_mvx);
	Tflux->SetBranchAddress("mvy",flux_mvy);
	Tflux->SetBranchAddress("mvz",flux_mvz);

	// Int_t nevent = (Int_t)Tgen->GetEntries();
	Int_t nevent = (Int_t)Tflux->GetEntries();
	cout << nevent << endl;
	/*End Set Branch}}}*/

	/*Read in each event{{{*/
	const int Electron = 11; // -11 for positron
	const int Gamma = 22;
	const int Beam = 0;
	double fmom = -1000.1, r = -1000.1, z = -1000.0;//cm
	int Just_One = 0;
	for(Int_t i=0;i<nevent;i++){
		cout<<i<<"\r";
		Tflux->GetEntry(i);

		if(flux_nfluxhit>=31&&Just_One==0){
			for (Int_t j=0;j<flux_nfluxhit;j++) {
//				if(*(flux_ID+j)>=5100000&&*(flux_ID+j)<=4120000) { //SPD
                if(1){
				//	if(*(flux_pid+j)==Gamma) continue; //SPD doesn't sense photon
	//				if(*(flux_mpid+j)==Beam) continue; //SPD doesn't sense photon

					fmom=sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2)+pow(*(flux_pz+j),2))/1e3;//cm
					r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2))/10.;//cm
					z =*(flux_z+j)/10.; //cm

					if(r > 210.0||r<96.0) continue;//The radius of a mrpc sector is 210cm;
                    Just_One=1;
					if(1){
					//if(*(flux_pid+j)==Electron&&*(flux_mpid+j)==Beam){
					cout << j << endl;
					cout << flux_evn<< " " << flux_nfluxhit << " " << *(flux_ID+j) << " " << *(flux_pid+j) << " " << *(flux_mpid+j)
						<< " " <<  *(flux_Edep+j) << " " << *(flux_E+j) << " " << *(flux_x+j) << " " << *(flux_y+j) << " " << *(flux_z+j) 
						<< " " << *(flux_vx+j) << " " << *(flux_vy+j) << " " << *(flux_vz+j)  << " " << *(flux_px+j) << " " << *(flux_py+j) << " " << *(flux_pz+j) << endl;
					}
					} //Gas & Electron & r>96,r<210

				}//Flux particles in one event
			}
			if(Just_One==1)
				return;

			}
			/*End Read in each event}}}*/
		file->Close(); 
}
