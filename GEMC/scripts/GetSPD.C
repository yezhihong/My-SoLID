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

void GetSPD(TString input_filename)
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

	/*Define SPD{{{*/
	//In each segmentation, 40 slides, each slides has 2.5cm width, 0.3cm between two slides.
	const double R_Min = 90.0;//cm, VP position, SPD 
	const double R_Max = 230.0;//cm, VP position, SPD
	const int SPD_Module = 30; //30 module around the circle
	const int SPD_Slide = 6;// put 6 slides in each module just for check the R-dependence, 
    const int EC_Mom_Bin = 23;

	double SPD_R[SPD_Slide];//Center location of each slide
	double SPD_In_G[SPD_Module][SPD_Slide];// Number of photons going into the device
	double SPD_In_E[SPD_Module][SPD_Slide];// Number of electrons going into the device
	double SPD_Out_G[SPD_Module][SPD_Slide];// Number of photons going into the device
	double SPD_Out_E[SPD_Module][SPD_Slide];// Number of electrons going into the device

	for(int i=0;i<SPD_Slide;i++){//I just want to initialize everything
		SPD_R[i] = trig_cut_range_R[i]; //{0,105,115,130,150,200,300};
		for(int j=0;j<SPD_Module;j++){
			SPD_In_G[j][i] =0.;
			SPD_In_E[j][i] =0.;
			SPD_Out_G[j][i] =0.;
			SPD_Out_E[j][i] =0.;
		}
	}
	double SPD_Threshold = 1.; //GeV for EC cut, the cut could be tight if using Jin's curves
	//double SPD_Threshold = 0.000001; //MeV for EC cut, the cut could be tight if using Jin's curves
	
	//Other Definition	
	const int Electron = 11;
	const int Gamma = 22;
	const int Beam = 0;
	const int Neutron = 2112;
	const int Neutrino1 = 12;//Nu_e
	const int Neutrino2 = 14;//Nu_Mu
	const int Neutrino3 = 16;//Nu_Tao
   	
	int VP_Before = 0; //SPD virtual plane
	int VP_After  = 0; //EC virtual plane
	double Z_SPD = 0.0; //cm VP position
	double Z_Out = 5.0; //cm VP position
	//Before MRPC
	if(input_filename.Contains("before")){
		VP_Before = 5110000; //SPD virtual plane
		VP_After  = 4110000; //MRPC virtual plane
		Z_SPD = 407.0; //cm VP position, 407.5 is the real positon of SPD
		Z_Out = 409.0; //cm MRPC VP position
	}
	if(input_filename.Contains("after")){
		//After MRPC
		VP_Before = 5110000; //SPD virtual plane
		VP_After  = 3110000; //EC virtual plane
		Z_SPD = 412.0; //cm VP position, 412.5 is the real positon of SPD
		Z_Out = 413.0; //cm EC VP position
	}
	const double Z_EC = 425.5; //cm, EC front, VP position = 413cm
	const double R_Min_Out = 90.0;//cm, HGC front
	const double R_Max_Out = 230.0;//cm,
	 
	double r = -1000.0, z = -1000.0, R_EC = -1000.0;//cm
	int Slide_ID = 0;
	int Module_ID = 0;

	/*}}}*/

	/*Read in each event{{{*/
	//Int_t nselected = 1e5;
	Int_t nselected = nevent;
	double In_G = 0, In_E = 0, Out_G = 0, Out_E = 0;
	double Count_High = 0, Count_Low = 0;
	int Count_Ep =0, Count_Em=0, Count_Both=0;
	double Count_Cut_Ep =0, Count_Cut_Em=0;
	double EC_Cut = 0.0;
	double fmom_old = 0.0;
	int Double_Count = 0;
	int FirstOne0 = 0; //incoming electron with E<1GeV
	int FirstOne1 = 0; //incoming electron with E>1GeV
	int FirstOne2 = 0; //incoming electron with R-dependence E-cut
	int FirstOne5 = 0; //outgoing positron with R-dependence E-cut
	int FirstOne3 = 0; //incoming electron with R-dependence E-cut by slides
	int FirstOne4 = 0; //outgoing electron with R-dependence E-cut by slides
	cerr<<"++++++++++++++++ "<<endl;
	for(Int_t i=0;i<nselected;i++){
		cout<<i<<"\r";

		Tflux->GetEntry(i);
		FirstOne0 = 0;
		FirstOne1 = 0;
		FirstOne2 = 0;
		FirstOne3 = 0;
		FirstOne5 = 0;
		//Double_Count = 0;
		for (Int_t j=0;j<flux_nfluxhit;j++) {
			r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2))/10.;//cm
			if(r > R_Max||r<R_Min) continue;//The radius of a mrpc sector is 210cm;
			double fmom=sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2)+pow(*(flux_pz+j),2))/1e3;//GeV
			if(*(flux_pz+j)<-1e-19)continue;

			if(fmom<1e-9) continue;

			double vertex_theta = atan(sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2))/(*(flux_pz+j)))*DEG;
			R_EC = abs((Z_EC - Z_SPD) * tan(vertex_theta/DEG) + r);
			if(R_EC>230.0 || R_EC<90.00) continue;

			//Selct the right Cut	
			//cut[Det_ID][Cut_ID][Data_Point][Cut_Info]: Cut_Info: R_Min, R_Max, P_Min, P_Max, e_Eff, pi_Eff
			EC_Cut = -1.0;
			for(int k=0;k<SPD_Slide;k++){
				for(int l=0;l<EC_Mom_Bin;l++){
					if(R_EC>trig_cut[0][k][l][0]&&R_EC<=trig_cut[0][k][l][1]){
						if(fmom>trig_cut[0][k][l][2]&&fmom<=trig_cut[0][k][l][3]){
							EC_Cut =trig_cut[0][k][l][4]; 
							Slide_ID = k;
						}
					}
				}
			}
			if(EC_Cut<-1e-9||EC_Cut>1){
				cerr<<"---- I can't find the cut!"<<Form(" --- R= %f,  E =%f", r, fmom)<<endl;
			   	return;
			}

			double hit_y = *(flux_y+j), hit_x = *(flux_x+j);	
			double hit_phi = fabs(atan(hit_y/hit_x)*DEG);
			Module_ID = (int) hit_phi/(360./SPD_Module);

			int ID_Pick = VP_Before;
			//Low Energy Electron <1GeV
			if((*(flux_pid+j)==Electron||*(flux_pid+j)==-Electron)&&*(flux_ID+j)==ID_Pick && fmom<SPD_Threshold){//#Eelectrons going out 
				if(FirstOne0<1)
					Count_Low+=1.0;
				FirstOne0 ++;;
			}
			//High Energy Electron >1GeV
			if((*(flux_pid+j)==Electron||*(flux_pid+j)==-Electron)&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold){//#Eelectrons going out 
				if(FirstOne1<1)
					Count_High+=1.0;
				FirstOne1 ++;;

				if(*(flux_pid+j)==Electron)
					Count_Em++;
				if(*(flux_pid+j)==-Electron)
					Count_Ep++;
			}
			//High Energy Electron with EC R-Cut 
			if((*(flux_pid+j)==Electron)&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold){//#Eelectrons going out 
				//Count_Em++;
				if(FirstOne2<1)
					Count_Cut_Em+=EC_Cut;
				FirstOne2 ++;;
			}
			//High Energy Positron with EC R-Cut 
			if((*(flux_pid+j)==-Electron)&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold){//#Eelectrons going out 
				//Count_Ep++;
				if(FirstOne5<1)
					Count_Cut_Ep+=EC_Cut;
				FirstOne5 ++;;
			}

            //Count by slides
			if(*(flux_pid+j)==Gamma&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold)//#Photons going in
				SPD_In_G[Module_ID][Slide_ID]+=EC_Cut;
			if((*(flux_pid+j)==Electron||*(flux_pid+j)==-Electron)&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold){//#Eelectrons going out 
				if(FirstOne3<1)
					SPD_In_E[Module_ID][Slide_ID]+=EC_Cut;
				FirstOne3 ++;;
			}
		}

		if((FirstOne2>=1) && (FirstOne5>=1))
			cerr<<" Oh! Double Counting !!!  "<< ++Double_Count<<endl;

		FirstOne4 = 0;
		for (Int_t j=0;j<flux_nfluxhit;j++) {
			r=sqrt(pow(*(flux_x+j),2)+pow(*(flux_y+j),2))/10.;//cm
			if(r<R_Min_Out && r>R_Max_Out) continue;
			double fmom=sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2)+pow(*(flux_pz+j),2))/1e3;//GeV
			if(*(flux_pz+j)<-1e-19)continue;
			if(fmom<1e-9) continue;
			double vertex_theta = atan(sqrt(pow(*(flux_px+j),2)+pow(*(flux_py+j),2))/(*(flux_pz+j)))*DEG;
			R_EC = abs((Z_EC - Z_Out) * tan(vertex_theta/DEG) + r);
			if(R_EC>230.0 || R_EC<90.00) continue;

			EC_Cut = -1.0;
			for(int k=0;k<SPD_Slide;k++){
				for(int l=0;l<EC_Mom_Bin;l++){
					if(R_EC>trig_cut[0][k][l][0]&&R_EC<=trig_cut[0][k][l][1]){
						if(fmom>trig_cut[0][k][l][2]&&fmom<=trig_cut[0][k][l][3]){
							EC_Cut =trig_cut[0][k][l][4]; 
							Slide_ID = k;
						}
					}
				}
			}
			if(EC_Cut<-1e-9||EC_Cut>1){
				cerr<<"---- I can't find the cut!"<<Form(" --- R= %f,  E =%f", r, fmom)<<endl;
			   	return;
			}

			double hit_y = *(flux_y+j), hit_x = *(flux_x+j);	
			double hit_phi = fabs(atan(hit_y/hit_x)*DEG);
			Module_ID = (int) hit_phi/(360./SPD_Module);

			int ID_Pick = VP_After;
			if(*(flux_pid+j)==Gamma&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold ){//#Photons going out 
				SPD_Out_G[Module_ID][Slide_ID]+=EC_Cut;
			}
			if((*(flux_pid+j)==Electron||*(flux_pid+j)==-Electron)&&*(flux_ID+j)==ID_Pick && fmom>=SPD_Threshold){//#Eelectrons going out 
				Count_Both++;
				if(FirstOne4<1)
					SPD_Out_E[Module_ID][Slide_ID]+=EC_Cut;
				FirstOne4 ++;;
			}
		}//Flux particles in one event

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
	ofstream outfile(Form("SPD_%s",input_out.Data()));
	for(int k=0;k<SPD_Slide;k++){
		for(int l=0;l<SPD_Module;l++){
			 In_G +=SPD_In_G[l][k];
			 In_E +=SPD_In_E[l][k];
			 Out_E +=SPD_Out_E[l][k];
			 Out_G +=SPD_Out_G[l][k];
		}
	}
	cerr<<" ====== In_E = "<< In_E*Count_To_Rate<<endl;
	cerr<<" ====== In_G = "<< In_G*Count_To_Rate<<endl;
	cerr<<" ====== Out_E = "<< Out_E*Count_To_Rate<<endl;
	cerr<<" ====== Out_G = "<< Out_G*Count_To_Rate<<endl;

	cerr<< Form("&&&& N (P<1GeV) = %d, N (P>1.0GeV) = %d, 1GeV Cut N = (-)%d / (+)%d /(+-)%d",
		   	(int)(Count_Low), (int)(Count_High), (int)(Count_Em), (int)(Count_Ep), (int)(Count_Both))<<endl;
	cerr<< Form("                                        R-Cut Cut N = (-)%d / (+)%d /(+-)%d", (int)(Count_Cut_Em),
		   	(int)(Count_Cut_Ep), (int)(In_E))<<endl;

	outfile<<" ====== In_E = "<< In_E*Count_To_Rate<<endl;
	outfile<<" ====== In_G = "<< In_G*Count_To_Rate<<endl;
	outfile<<" ====== Out_E = "<< Out_E*Count_To_Rate<<endl;
	outfile<<" ====== Out_G = "<< Out_G*Count_To_Rate<<endl;
    /*}}}*/
}
