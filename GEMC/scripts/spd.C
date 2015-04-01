#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
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

void spd(string input_filename)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const double DEG=180./3.1415926;

	/*Set Histogram{{{*/
	char the_filename[200];
	sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());

	char output_filename[200];
	sprintf(output_filename, "%s_output.root",the_filename);
	TFile *outputfile=new TFile(output_filename, "recreate");

	const int n=2;

	double P_Min = -0.5;//GeV;
	double P_Max = 11.5;//GeV;
    double T_Min = 0.0; //Deg
    double T_Max = 50.0;//Deg
    double Z_Min = 406.0;//cm
    double Z_Max = 408.0;//cm
    double R_Min = 90.0;//cm
    double R_Max = 220.0;//cm

    double VT_Min = -180.0;//Deg
    double VT_Max = 180.0;//Deg
    double VZ_Min = 300;//cm
    double VZ_Max = 500.0;//cm
    double VR_Min = 90.0;//cm
    double VR_Max = 220.0;//cm

    //Generated 
	TH2F *hgen=new TH2F("gen","gen",250,T_Min,T_Max,220,P_Min,P_Max);   
	hgen->SetTitle("generated particles; Theta (deg); P (GeV)");
	TH2F *hgen_vertexZ=new TH2F("gen_vertexZ","gen_vertexZ",350,VT_Min,VT_Max,50,VZ_Min,VZ_Max);
	hgen_vertexZ->SetTitle("generated particles; vertex Theta (deg); vertex Z (cm)");
	TH2F *hgen_vertexR=new TH2F("gen_vertexR","gen_vertexR",350,VT_Min,VT_Max,50,VR_Min,VR_Min);
	hgen_vertexR->SetTitle("generated particles; vertex Theta (deg); vertex R (cm)");

	TH1F *hflux_mom[n],*hflux_theta[n];
	TH2F *hflux[n],*hflux_vertexZ[n],*hflux_vertexR[n];
	TH2F *hhit_rMom[n],*hhit_xy[n];
	TH2F *hhit_phidiffMom[n],*hhit_thetadiffMom[n];
	TH2F *hhit_rz=new TH2F("hit_rz","hit_rz",1000,Z_Min,Z_Max,300,0,300);  
	for(int i=0;i<n;i++){
		char hstname[100];  
		sprintf(hstname,"flux_mom_%i",i);
		hflux_mom[i]=new TH1F(hstname,hstname,1000,P_Min,P_Max);
		sprintf(hstname,"flux_theta_%i",i);
		hflux_theta[i]=new TH1F(hstname,hstname,1000,T_Min,T_Max); 
		sprintf(hstname,"flux_mom_R_%i",i);

		//Vertex
		sprintf(hstname,"flux_%i",i);
		hflux[i]=new TH2F(hstname,hstname,1000,VT_Min,VT_Max,1000,P_Min,P_Max);     
		hflux[i]->SetTitle("particles detected by MRPC;vertex Theta (deg);P (GeV)");
		sprintf(hstname,"flux_vertexZ_%i",i);
		hflux_vertexZ[i]=new TH2F(hstname,hstname,1000,VT_Min,VT_Max,1000,VZ_Min,VZ_Max);
		hflux_vertexZ[i]->SetTitle("particles detected by MRPC;vertex Theta (deg);vertex Z (cm)");   
		sprintf(hstname,"flux_vertexR_%i",i);
		hflux_vertexR[i]=new TH2F(hstname,hstname,1000,VT_Min,VT_Max,1000,VR_Min,VR_Max);
		hflux_vertexR[i]->SetTitle("particles detected by MRPC;vertex Theta (deg);vertex R (cm)");      

		sprintf(hstname,"hit_rMom_%i",i);
		hhit_rMom[i]=new TH2F(hstname,hstname,1000,0,300,1000,0,P_Max);  
		sprintf(hstname,"hit_phidiffMom_%i",i);
		hhit_phidiffMom[i]=new TH2F(hstname,hstname,7200,-360,360,220,P_Min,P_Max);  
		sprintf(hstname,"hit_thetadiffMom_%i",i);
		hhit_thetadiffMom[i]=new TH2F(hstname,hstname,3600,-180,180,220,P_Min,P_Max);  

		sprintf(hstname,"hit_xy_%i",i);
		hhit_xy[i]=new TH2F(hstname,hstname,1000,-300,300,1000,-300,300);     
	}

	TH3F *hxyz = new TH3F("hxzy","Electron Cascade in MRPC for one electron", 20, -220, 220, 20, -220, 220, 20, 400.,422.);
    hxyz->SetTitle("x (cm); y (cm) ; z (cm)");

	/*End Set Histo}}}*/

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
	Int_t nselected = 0;
	cout << nevent << endl;
    /*End Set Branch}}}*/

	/*Read in each event{{{*/
	double E_Dep = 0.0;
	const int Electron = 11;
	for(Int_t i=0;i<nevent;i++){
		cout << i << "\r";
		double theta=0,mom=0,phi=0,vx=0,vy=0,vz=0;  
		Tgen->GetEntry(i);  
		for(Int_t j=0;j<gen_ngen;j++){
			theta=*(gen_theta+j);
			phi=*(gen_phi+j);      
			mom=*(gen_p+j)/1e3;
			vx=*(gen_vx+j);
			vy=*(gen_vy+j);
			vz=*(gen_vz+j);      

			hgen->Fill(theta,mom);
			hgen_vertexZ->Fill(theta,vz);
			hgen_vertexR->Fill(theta,sqrt(vx*vx+vy*vy));      
		}  

		int counter_hit[n]={0,0};
		Tflux->GetEntry(i);    
		for (Int_t j=0;j<flux_nfluxhit;j++) {

			//if(*(flux_mpid+j)==Electron&&*(flux_pid+j)==Electron){//Only electron
			//if(*(flux_pid+j)==Electron){//Only electron
			if(1){
				int detector_ID=*(flux_ID+j)/1000000;
				int subdetector_ID=(*(flux_ID+j)%1000000)/100000;
				int subsubdetector_ID=((*(flux_ID+j)%1000000)%100000)/10000;

				double fmom=sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2)+pow(*(flux_pz+j),2))/1e3;
				double r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2))/10.;
				double hit_theta=atan((r-sqrt(vx*vx+vy*vy)/10.0)/(*(flux_z+j)-vz))*DEG;

                if(fmom<1e-9) continue;
				if(r > 210.0||r<96.0) continue;//The radius of a mrpc sector is 210cm;

				int hit_id=-1;
				if (*(flux_ID+j)>=5100000) hit_id=0; //SPD
				if (*(flux_ID+j)>=5110000) hit_id=1; //SPD frnt virtual
				if (hit_id==-1) continue;  //skip other subsubdetector

				counter_hit[hit_id]++;

				double hit_y=*(flux_y+j),hit_x=*(flux_x+j);
				double hit_phi=fabs(atan(hit_y/hit_x)*DEG);
				if (hit_y>0 && hit_x>0) hit_phi=hit_phi;
				if (hit_y>0 && hit_x<0) hit_phi=180-hit_phi;
				if (hit_y<0 && hit_x<0) hit_phi=180+hit_phi;
				if (hit_y<0 && hit_x>0) hit_phi=360-hit_phi;    

				hflux_mom[hit_id]->Fill(fmom);
				hhit_rz->Fill(*(flux_z+j)/10,r);
				hhit_phidiffMom[hit_id]->Fill(hit_phi-phi,fmom);
				hhit_thetadiffMom[hit_id]->Fill(hit_theta-theta,fmom);
				hhit_rMom[hit_id]->Fill(r,fmom);   
			
				hhit_xy[hit_id]->Fill(*(flux_x+j)/10,*(flux_y+j)/10);    
				hflux[hit_id]->Fill(theta,fmom);
				hflux_vertexZ[hit_id]->Fill(theta,*(flux_vz+j)/10);
				hflux_vertexR[hit_id]->Fill(theta,sqrt(pow(*(flux_vx+j),2)+pow(*(flux_vy+j),2))/10.);
			}//Only Electron        

			for (Int_t j=0;j<flux_nfluxhit;j++) {
				if(*(flux_mpid+j)==Electron&&*(flux_pid+j)==Electron){//Only electron
					if(counter_hit[0]>1e-33){	
						double r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2));
						hhit_rz->Fill(*(flux_z+j)/10,r);
					}
				}
			}

		}//Flux particles in one events
	}
	file->Close();
    /*End Read in each event}}}*/

	// gStyle->SetOptStat(0);
	outputfile->Write();
	outputfile->Flush();
}
