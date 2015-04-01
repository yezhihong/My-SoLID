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
#include <TGraphErrors.h>
#include <TGraph.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void SPDMRPC(TString input_filename)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const double DEG=180./3.1415926;

	TString input_png = input_filename;
	input_png.ReplaceAll(".root","_zhit.png");

	/*Set Branch{{{*/
	TFile *file=new TFile(input_filename.Data());
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

	/* EC Electron Trigger{{{*/
	const int Ntrigline=6,Ntriglinebin=21;
	int region_index;
	//	if (input_filename.Contains("SIDIS_FA")) region_index=0;
	//	else if (input_filename.Contains("SIDIS_LA")) region_index=1;
	//	else {cout << "need option for FA or LA region" << endl; exit(-1);}
	region_index = 0;

	int det[2]={8,12};  //detecor ID
	double Rmin[2]={90,80};
	double Rmax[2]={230,140};

	TFile *file_trig_cut[2][Ntrigline][2];//[Det_ID][Cut_ID][PID], Det_DC: 0->FA, 1->LA, PID: 0->e, 1->pi
	// Radius(cm)  P Threshold (GeV)
	//   90 - 105      5.0    
	// 105 - 115      4.0    
	// 115 - 130      3.0    
	// 130 - 150      2.0    
	// 150 - 200      1.0 
	// 200 - 230      2.0 
	// 6 point cut, right on Q2=1 line and field bend line
	double trig_cut_range_R[Ntrigline+1]={0,105,115,130,150,200,300};

	file_trig_cut[0][0][0]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH4.4.root");
	file_trig_cut[0][0][1]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH4.4.root");
	file_trig_cut[0][1][0]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH3.5.root");
	file_trig_cut[0][1][1]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH3.5.root");
	file_trig_cut[0][2][0]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH2.6.root");
	file_trig_cut[0][2][1]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH2.6.root");
	file_trig_cut[0][3][0]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root");
	file_trig_cut[0][3][1]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root");
	file_trig_cut[0][4][0]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_Trig0.9.root");
	file_trig_cut[0][4][1]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_Trig0.9.root");
	file_trig_cut[0][5][0]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root");
	file_trig_cut[0][5][1]=new TFile("triggerfile/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root");

	///Large Angle Trigger   has no radial dependence
	for (int i=0;i<Ntrigline;i++){
		file_trig_cut[1][i][0]=new TFile("cut1GeV_innerbackground/Lead2X0PbBlock_Hex.1.SIDIS_Large_RunElectron_GetEfficienciesBackGround_Oct2013_SIDIS_Full_bgd_TrigSH2.0.root");
		file_trig_cut[1][i][1]=new TFile("cut1GeV_innerbackground/Lead2X0PbBlock_Hex.1.SIDIS_Large_RunPion_GetEfficienciesBackGround_Oct2013_SIDIS_Full_bgd_TrigSH2.0.root");
	}


	char *gr_trig_cut_ele_name[2][Ntrigline],*gr_trig_cut_pi_name[2][Ntrigline];
	for(int j=0;j<2;j++){
		for(int i=0;i<Ntrigline;i++){
			gr_trig_cut_ele_name[j][i]="Graph";
			gr_trig_cut_pi_name[j][i]="Graph";    
		}  
	}

	double trig_cut[2][Ntrigline][Ntriglinebin+2][6];
	TGraphErrors *gr_trig_cut_ele[2][Ntrigline],*gr_trig_cut_pi[2][Ntrigline];
	for (int j=0;j<2;j++){
		for (int i=0;i<Ntrigline;i++){  
			gr_trig_cut_ele[j][i]=(TGraphErrors*) file_trig_cut[j][i][0]->Get(gr_trig_cut_ele_name[j][i]);
			gr_trig_cut_pi[j][i]=(TGraphErrors*) file_trig_cut[j][i][1]->Get(gr_trig_cut_pi_name[j][i]);
			double binwidth=gr_trig_cut_ele[j][i]->GetX()[1]-gr_trig_cut_ele[j][i]->GetX()[0];
			if (j==1) { //add one more point for LA to become 21 points like FA
				gr_trig_cut_ele[j][i]->SetPoint(20,gr_trig_cut_ele[j][i]->GetX()[19]+binwidth,gr_trig_cut_ele[j][i]->GetY()[19]);
				gr_trig_cut_pi[j][i]->SetPoint(20,gr_trig_cut_pi[j][i]->GetX()[19]+binwidth,gr_trig_cut_pi[j][i]->GetY()[19]);    
			}
			for (int k=0;k<Ntriglinebin;k++){ //set any point with eff <0.01 as 0
				if (gr_trig_cut_ele[j][i]->GetY()[k]<0.01) gr_trig_cut_ele[j][i]->SetPoint(k,gr_trig_cut_ele[j][i]->GetX()[k],0.);
				if (gr_trig_cut_pi[j][i]->GetY()[k]<0.01) gr_trig_cut_pi[j][i]->SetPoint(k,gr_trig_cut_pi[j][i]->GetX()[k],0.);
			}
			trig_cut[j][i][0][0]=trig_cut_range_R[i];
			trig_cut[j][i][0][1]=trig_cut_range_R[i+1];
			trig_cut[j][i][0][2]=0.;
			trig_cut[j][i][0][3]=gr_trig_cut_ele[j][i]->GetX()[0]-binwidth/2;
			trig_cut[j][i][0][4]=0.;
			trig_cut[j][i][0][5]=0.;
			//	cout << j << " " << i << " " << 0 << " " << gr_trig_cut_ele_name[j][i] << "\t" << gr_trig_cut_ele[j][i]->GetX()[0] << "\t" << gr_trig_cut_ele[j][i]->GetY()[0] << "\t" << gr_trig_cut_pi_name[j][i] << "\t" << gr_trig_cut_pi[j][i]->GetX()[0] << "\t" << gr_trig_cut_pi[j][i]->GetY()[0] << endl;
			for (int k=0;k<Ntriglinebin;k++){
				trig_cut[j][i][k+1][0]=trig_cut_range_R[i];
				trig_cut[j][i][k+1][1]=trig_cut_range_R[i+1];
				trig_cut[j][i][k+1][2]=gr_trig_cut_ele[j][i]->GetX()[k]-binwidth/2;
				trig_cut[j][i][k+1][3]=gr_trig_cut_ele[j][i]->GetX()[k]+binwidth/2;
				trig_cut[j][i][k+1][4]=gr_trig_cut_ele[j][i]->GetY()[k];
				trig_cut[j][i][k+1][5]=gr_trig_cut_pi[j][i]->GetY()[k];
				//		cout << j << " " << i << " " << k+1 << " " << gr_trig_cut_ele_name[j][i] << "\t" << gr_trig_cut_ele[j][i]->GetX()[k] << "\t" << gr_trig_cut_ele[j][i]->GetY()[k] << "\t" << gr_trig_cut_pi_name[j][i] << "\t" << gr_trig_cut_pi[j][i]->GetX()[k] << "\t" << gr_trig_cut_pi[j][i]->GetY()[k] << endl;
			}
			trig_cut[j][i][Ntriglinebin+1][0]=trig_cut_range_R[i];
			trig_cut[j][i][Ntriglinebin+1][1]=trig_cut_range_R[i+1];
			trig_cut[j][i][Ntriglinebin+1][2]=gr_trig_cut_ele[j][i]->GetX()[Ntriglinebin-1]+binwidth/2;
			trig_cut[j][i][Ntriglinebin+1][3]=11.;
			trig_cut[j][i][Ntriglinebin+1][4]=gr_trig_cut_ele[j][i]->GetY()[Ntriglinebin-1];
			trig_cut[j][i][Ntriglinebin+1][5]=gr_trig_cut_pi[j][i]->GetY()[Ntriglinebin-1];
			//	cout << j << " " << i << " " << Ntriglinebin+1 << " " << gr_trig_cut_ele_name[j][i] << "\t" << gr_trig_cut_ele[j][i]->GetX()[Ntriglinebin-1] << "\t" << gr_trig_cut_ele[j][i]->GetY()[Ntriglinebin-1] << "\t" << gr_trig_cut_pi_name[j][i] << "\t" << gr_trig_cut_pi[j][i]->GetX()[Ntriglinebin-1] << "\t" << gr_trig_cut_pi[j][i]->GetY()[Ntriglinebin-1] << endl;
			gr_trig_cut_ele[j][i]->Delete();
			gr_trig_cut_pi[j][i]->Delete();
		}
	}
	/*
	   cout << "here is the trig value" << endl;
	   for (int j=0;j<2;j++){
	   for (int i=0;i<Ntrigline;i++){  
	   for (int k=0;k<Ntriglinebin+2;k++){
	   for (int l=0;l<6;l++){
	//cut[Det_ID][Cut_ID][Data_Point][Cut_Info]
	cout << trig_cut[j][i][k][l] << "\t\t";
	}
	cout << endl;      
	}
	}
	}
	*/
	/*}}}*/

	/*Define MRPC{{{*/
	double MRPC_G = 0, MRPC_E = 0;
	//In each segmentation, 40 slides, each slides has 2.5cm width, 0.3cm between two slides.
	const double MRPC_R_Min = 96.0;
	const double MRPC_R_Max = 205.0;//Originally 210 but consider the size and the gap of each slide
    const int MRPC_Threshold = 4; //The Threshold to deterime how many hits means the slide has been fired
	const int MRPC_Module = 50; //50 module around the circle
	const int MRPC_Slide = 40;// 40 slides in each module, 
	const double MRPC_Width = 2.5; //cm
	const double MRPC_Gap = 0.3; //cm
	double MRPC_R[MRPC_Slide];//Center location of each slide
	int MRPC_Fire[MRPC_Slide];// I sum all 50 modules together, and devide the counts in each slides by 50
	int MRPC_Hit[MRPC_Slide];//each slide has 10 layers, count how many hits in these 10 layers for each event 
	for(int i=0;i<MRPC_Slide;i++){//I just want to initialize everything
		MRPC_R[i] = MRPC_Width/2.0 + i*(MRPC_Width+MRPC_Gap); //Notethat the edge of the first slide is at R=96
        MRPC_Fire[i] =0;
		MRPC_Hit[i] =0;
	}

	//In each player
	const double MRPC_Z = 410; //cm, MRPC Center
	const double MRPC_Gas_Thick =0.025; //cm
	const int MRPC_Gas_Layer =10; //cm
	const double MRPC_Glass_Thick =0.07; //cm const double Hit_Min = -0.5;
	const double Gas[MRPC_Gas_Layer] = {-(2.175+0.95*4)/10., -(2.175+0.95*3)/10., -(2.175+0.95*2)/10., -(2.175+0.95)/10., -2.175/10.,
		2.175/10., (2.175+0.95)/10., (2.175+0.95*2)/10., (2.175+0.95*3)/10., (2.175+0.95*4)/10.};//cm

/*}}}*/

	/*Define SPD{{{*/
	//In each segmentation, 40 slides, each slides has 2.5cm width, 0.3cm between two slides.
	double SPD_Threshold = 0.00001; //SPD is very sensitive to energy deposit
    double EDep_Total=0.0;//MeV
	double SPD_G = 0, SPD_E = 0;
	double SPD_Z = 0;
	if(input_filename.Contains("before")){//before means that spd is before mrpc
		SPD_Z = 412.5; //cm  the real positon of SPD
	}
	if(input_filename.Contains("after")){//after means that spd is before mrpc
		//Before MRPC
		SPD_Z = 407.5; //cm  the real positon of SPD
	}	
	/*}}}*/

	/*Other Definition{{{*/	
	const int Electron = 11;
	const int Gamma = 22;
	const int Beam = 0;
	const int Neutron = 2112;
	const int Neutrino1 = 12;//Nu_e
	const int Neutrino2 = 14;//Nu_Mu
	const int Neutrino3 = 16;//Nu_Tao

	double r = -1000.0, z = -1000.0, R_EC = -1000.0;//cm
	int nselected = nevent;
	int Slide_ID = 0;
	int Module_ID = 0;
	/*}}}*/

	/*Read in each event{{{*/
	//Int_t nselected = 1e5;
	bool Is_Gamma = kFALSE;
	bool Is_Electron = kFALSE;
	cerr<<"++++++++++++++++ "<<endl;
	for(Int_t i=0;i<nselected;i++){
		cout<<i<<"\r";

		Tflux->GetEntry(i);
		Is_Gamma = kFALSE;
		Is_Electron = kFALSE;

		/*SPD{{{*/
			EDep_Total = 0.0;
			for (Int_t j=0;j<flux_nfluxhit;j++) {
				r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2))/10.;//cm
				if(r > 210.0||r<96.0) continue;//The radius of a mrpc sector is 210cm;
				//if(r > R_Max||r<R_Min) continue;//The radius of a mrpc sector is 210cm;
				double fmom=sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2)+pow(*(flux_pz+j),2))/1e3;//GeV
				if(*(flux_pz+j)<-1e-19)continue;

				if(fmom<1e-9) continue;

				if((*(flux_pid+j)==Electron||*(flux_pid+j)==-Electron)&&*(flux_ID+j)==5100000){//#Eelectrons going out 
					EDep_Total += *(flux_Edep+j); //MeV
				}
			}//End Read in each flux particle
			if(EDep_Total>=SPD_Threshold){//the deposit energy convert into light and then into photoelectrons
				SPD_E++;
				Is_Electron=kTRUE; Is_Gamma=kFALSE;
			}
		/*}}}*/

		/*MRPC{{{*/
		if(Is_Electron && !(Is_Gamma)){
		for(int k=0;k<MRPC_Slide;k++){
			MRPC_Hit[k] = 0;
		}

		for (Int_t j=0;j<flux_nfluxhit;j++) {
			if(*(flux_ID+j)>=4000000&&*(flux_ID+j)<4120000) { //ID of 10 gas layers: 4100000 - 41000009
				if(*(flux_pid+j)==Electron||*(flux_pid+j)==-Electron){//MRPC doesn't sense photon
					r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2))/10.;//cm
					if(r > MRPC_R_Max||r<MRPC_R_Min) continue;//The radius of a mrpc sector is 210cm;

					z =*(flux_z+j)/10. - MRPC_Z; //cm

					Slide_ID = (int)((r-96)/(MRPC_Width+MRPC_Gap));

					if(abs(r-(MRPC_R_Min+MRPC_R[Slide_ID]))>MRPC_Width/2.0) continue; //Double Check whether it is still in the Gap
					//cerr<<Form("=== R = %f, Slide_ID = %d, E= %e", r, Slide_ID, *(flux_E+j)/1e3)<<endl;

					for(int k=0;k<MRPC_Gas_Layer;k++){//Check 10layers
						if(abs(z-Gas[k])<=MRPC_Gas_Thick/2.0){ //Charge particle in the gas only
							MRPC_Hit[Slide_ID]++;
						}//Gas Only
					}//Check 10layers
				}// Electron & r>96,r<210, no in the gap
			}//Gass Layers
		}//Flux particles in one event

		for(int k=0;k<MRPC_Slide;k++){
			if(MRPC_Hit[k]>MRPC_Threshold){//If the hits are above the threshold, then I treat this slide has been fired
				MRPC_Fire[k]++; 
				//cerr<<"--- Total#="<<flux_nfluxhit<<", Slide#"<<k<<" has been fired! Count = "<<MRPC_Fire[k]<<endl;
				//cerr<<"--- Total#="<<flux_nfluxhit<<", Slide#"<<k<<" has been fired! Count = "<<MRPC_Hit[k]<<endl;
			}
		}
		}//if(Is_Electron && !(Is_Gamma))
		/*}}}*/
}
/*End Read in each event}}}*/

	/*Output Results{{{*/
	double Count_To_Rate = 0.0;
	if(input_filename.Contains("pi0")){
		Count_To_Rate = 241.0/1e6; //Count to MHz
		if(input_filename.Contains("up")||input_filename.Contains("down")){
			Count_To_Rate = 136.0/1e6;//Count to MHz
		}
	}
	else
		Count_To_Rate = 1.5/1.6; //Count to MHz for 15uA electron events;
	cerr<<" --- Using Converting Factor = "<<Count_To_Rate<<endl;

	TString input_out = input_filename;
	input_out.ReplaceAll(".root",".out");
	ofstream outfile(Form("SPDMRPC_%s",input_out.Data()));

	for(int k=0;k<MRPC_Slide;k++){
	    MRPC_E +=MRPC_Fire[k];
	}
	cerr<<" ====== SPD_E = "<< SPD_E*Count_To_Rate<<endl;
	outfile<<" ====== SPD_E = "<< SPD_E*Count_To_Rate<<endl;

	cerr<<" ====== MRPC_E = "<< MRPC_E*Count_To_Rate<<endl;
	outfile<<" ====== MRPC_E = "<< MRPC_E*Count_To_Rate<<endl;

	/*}}}*/
}
