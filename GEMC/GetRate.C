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
/*}}}*/

using namespace std;

void GetRate(const TString& Detector_Name, const string input_filename)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const Double_t DEG=180./3.1415926;
	const Double_t PI = 3.1415926;

	TFile *file=new TFile(input_filename.c_str());
	if (file->IsZombie()) {
		cout << "Error opening file" << input_filename << endl;
		//continue;
		exit(-1);
	}
	else cout << "### Open file " << input_filename << endl;
	
	/*Temp Rate Normalization{{{*/
	Double_t filenum=1.;
	if (input_filename.find("_EM_",0) != string::npos) {
		cout << "### EM background from beam on target" <<  endl;
	} 
	else{
		cout << "### background from normalized event generator with all interaction" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "### filenum = " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "### we need filenum for addtional normalization" << endl; return;}
	}	  	
    /*}}}*/

	/*Set Branch{{{*/
	//Header Tree:
	// Var#1~#8 are free slots for propogating important info from the "INPUT generator seed"
	// For example, they can be used to store the cross section and other physics quantities
	// In eicRate, we store the following quantities:
	// var1->Wprate, var2->Wmrate, var3->targetPol, var4->x,var5->y, var6->W, var7->Q2, var8->rate 
	//
	TTree *header = (TTree*) file->Get("header");
	vector <Double_t> *head_evn=0,*head_evn_type=0; //Note: Vectors have to be initialized at first!!!
	vector <Double_t> *head_beamPol=0;
	vector<Double_t> *head_Wmrate=0, *head_Wprate=0, *head_targetPol=0, *head_x=0, *head_Q2=0, *head_W=0, *head_rate=0, *head_y=0;
	header->SetBranchAddress("evn",&head_evn);
	header->SetBranchAddress("evn_type",&head_evn_type);
	header->SetBranchAddress("beamPol",&head_beamPol);
	header->SetBranchAddress("var1",    &head_Wprate);
	header->SetBranchAddress("var2",    &head_Wmrate);
	header->SetBranchAddress("var3",    &head_targetPol);
	header->SetBranchAddress("var4",    &head_x);
	header->SetBranchAddress("var5",    &head_y);
	header->SetBranchAddress("var6",    &head_W);
	header->SetBranchAddress("var7",    &head_Q2);
	header->SetBranchAddress("var8",    &head_rate);

	TTree *generated = (TTree*) file->Get("generated");
	vector <Int_t> *gen_pid=0;
	vector <Double_t> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
	generated->SetBranchAddress("pid",&gen_pid);
	generated->SetBranchAddress("px",&gen_px);
	generated->SetBranchAddress("py",&gen_py);
	generated->SetBranchAddress("pz",&gen_pz);
	generated->SetBranchAddress("vx",&gen_vx);
	generated->SetBranchAddress("vy",&gen_vy);
	generated->SetBranchAddress("vz",&gen_vz);

	TTree *flux = (TTree*) file->Get("flux");
	vector<Double_t> *flux_id=0,*flux_hitn=0,*flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<Double_t> *flux_trackE=0,*flux_totEdep=0;
	vector<Double_t> *flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0;
	vector<Double_t> *flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
	flux->SetBranchAddress("hitn",&flux_hitn); //number of flux hits in one event
	flux->SetBranchAddress("id",&flux_id); //detector ID
	flux->SetBranchAddress("pid",&flux_pid); //particle ID
	flux->SetBranchAddress("mpid",&flux_mpid);  //mother particle ID
	flux->SetBranchAddress("tid",&flux_tid); //flux track ID
	flux->SetBranchAddress("mtid",&flux_mtid); //mother flux track ID
	flux->SetBranchAddress("otid",&flux_otid);//orginal flux track ID ?
	flux->SetBranchAddress("trackE",&flux_trackE);//energy of the flux particle
	flux->SetBranchAddress("totEdep",&flux_totEdep); //deposited energy of the flux particle
	flux->SetBranchAddress("avg_x",&flux_avg_x);
	flux->SetBranchAddress("avg_y",&flux_avg_y);
	flux->SetBranchAddress("avg_z",&flux_avg_z);
	flux->SetBranchAddress("avg_lx",&flux_avg_lx);
	flux->SetBranchAddress("avg_ly",&flux_avg_ly);
	flux->SetBranchAddress("avg_lz",&flux_avg_lz);
	flux->SetBranchAddress("avg_t",&flux_avg_t);
	flux->SetBranchAddress("px",&flux_px);
	flux->SetBranchAddress("py",&flux_py);
	flux->SetBranchAddress("pz",&flux_pz);
	flux->SetBranchAddress("vx",&flux_vx);
	flux->SetBranchAddress("vy",&flux_vy);
	flux->SetBranchAddress("vz",&flux_vz);
	flux->SetBranchAddress("mvx",&flux_mvx);
	flux->SetBranchAddress("mvy",&flux_mvy);
	flux->SetBranchAddress("mvz",&flux_mvz);
	/*End Set Branch}}}*/
	Int_t nevent = (Int_t)generated->GetEntries();
	cout << "nevent = " << nevent << endl;

	/*EC Electron Trigger{{{*/
	const Int_t Ntrigline=6,Ntriglinebin=21;
	Int_t region_index;
	//	if (input_filename.Contains("SIDIS_FA")) region_index=0;
	//	else if (input_filename.Contains("SIDIS_LA")) region_index=1;
	//	else {cout << "need option for FA or LA region" << endl; exit(-1);}
	region_index = 0;

	//Int_t det[2]={8,12};  //detecor ID
	//Double_t rmin[2]={105,80};
	//Double_t rmax[2]={235,140};

	TFile *file_trig_cut[2][Ntrigline][2];//[Det_ID][Cut_ID][PID], Det_DC: 0->FA, 1->LA, PID: 0->e, 1->pi
	// Radius(cm)  P Threshold (GeV)
	//   90 - 105      5.0    
	// 105 - 115      4.0    
	// 115 - 130      3.0    
	// 130 - 150      2.0    
	// 150 - 200      1.0 
	// 200 - 230      2.0 
	// 6 point cut, right on Q2=1 line and field bend line
	Double_t trig_cut_range_R[Ntrigline+1]={0,105,115,130,150,200,300};

	TString Trigger_Dir_Rad="./triggerfile/SIDIS_He3_201311/cutRadial_innerbackground/";
	file_trig_cut[0][0][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH4.4.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][0][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH4.4.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][1][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH3.5.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][1][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH3.5.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][2][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH2.6.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][2][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH2.6.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][3][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][3][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][4][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_Trig0.9.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][4][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_Trig0.9.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][5][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunElectron_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root",Trigger_Dir_Rad.Data()));
	file_trig_cut[0][5][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Forward_RunPion_GetEfficiencies_BackGround_Oct2013_SIDIS_TrigSH1.6.root",Trigger_Dir_Rad.Data()));

	///Large Angle Trigger   has no radial dependence
	TString Trigger_Dir_1GeVCut = "./triggerfile/SIDIS_He3_201311/cut1GeV_innerbackground/";
	for (Int_t i=0;i<Ntrigline;i++){
		file_trig_cut[1][i][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Large_RunElectron_GetEfficienciesBackGround_Oct2013_SIDIS_Full_bgd_TrigSH2.0.root",Trigger_Dir_1GeVCut.Data()));
		file_trig_cut[1][i][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Large_RunPion_GetEfficienciesBackGround_Oct2013_SIDIS_Full_bgd_TrigSH2.0.root",Trigger_Dir_1GeVCut.Data()));
	}


	char *gr_trig_cut_ele_name[2][Ntrigline],*gr_trig_cut_pi_name[2][Ntrigline];
	for(Int_t j=0;j<2;j++){
		for(Int_t i=0;i<Ntrigline;i++){
			gr_trig_cut_ele_name[j][i]="Graph";
			gr_trig_cut_pi_name[j][i]="Graph";    
		}  
	}

	Double_t trig_cut[2][Ntrigline][Ntriglinebin+2][6];
	TGraphErrors *gr_trig_cut_ele[2][Ntrigline],*gr_trig_cut_pi[2][Ntrigline];
	for (Int_t j=0;j<2;j++){
		for (Int_t i=0;i<Ntrigline;i++){  
			gr_trig_cut_ele[j][i]=(TGraphErrors*) file_trig_cut[j][i][0]->Get(gr_trig_cut_ele_name[j][i]);
			gr_trig_cut_pi[j][i]=(TGraphErrors*) file_trig_cut[j][i][1]->Get(gr_trig_cut_pi_name[j][i]);
			Double_t binwidth=gr_trig_cut_ele[j][i]->GetX()[1]-gr_trig_cut_ele[j][i]->GetX()[0];
			if (j==1) { //add one more point for LA to become 21 points like FA
				gr_trig_cut_ele[j][i]->SetPoint(20,gr_trig_cut_ele[j][i]->GetX()[19]+binwidth,gr_trig_cut_ele[j][i]->GetY()[19]);
				gr_trig_cut_pi[j][i]->SetPoint(20,gr_trig_cut_pi[j][i]->GetX()[19]+binwidth,gr_trig_cut_pi[j][i]->GetY()[19]);    
			}
			for (Int_t k=0;k<Ntriglinebin;k++){ //set any point with eff <0.01 as 0
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
			for (Int_t k=0;k<Ntriglinebin;k++){
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

	for (Int_t i=0;i<Ntrigline;i++){
		file_trig_cut[0][i][0]->Close();
		file_trig_cut[0][i][1]->Close();
		file_trig_cut[1][i][0]->Close();
		file_trig_cut[1][i][1]->Close();
	}
	
	const Int_t EC_Trigger_Slide = Ntrigline;// put 6 slides in each module just for check the R-dependence, 
	const Int_t EC_Trigger_Mom_Bin = Ntriglinebin+2;
	Double_t EC_R[EC_Trigger_Slide];//Center location of each slide
	for(Int_t i=0;i<EC_Trigger_Slide;i++){
		EC_R[i] = trig_cut_range_R[i]; //{0,105,115,130,150,200,300};
	}
	/*}}}*/

	/*Define Detector{{{*/
	Int_t VP_In = 0; //Front virtual plane
	Double_t Z_In = 0.0; //cm Front VP position
	 Double_t R_Min_In = 0.0;//cm, Front VP position
	 Double_t R_Max_In = 0.0;//cm, Front VP position

	Int_t VP_Out  = 0; //Rear virtual plane
	Double_t Z_Out = 0.0; //cm Rear VP position
	Double_t R_Min_Out = 0.0;//cm, Rear VP position
	Double_t R_Max_Out = 0.0;//cm, Rear VP position

	Int_t VP_EC = 0; //Front virtual plane
	Double_t Z_EC = 0.0; //cm, EC virtual plane front, VP position = 413cm
	Double_t R_Min_EC = 0.0;//cm, EC VP position
	Double_t R_Max_EC = 0.0;//cm, EC VP position
	Int_t EC_ID = 0; //0->FAEC, 1->LAEC
	Double_t EC_Threshold =-1.; //GeV for EC cut, the cut could be tight if using Jin's curves

	//For Segmentation only
	Int_t Module = 30; //30 module around the circle
	Int_t Slide = 6;// put 6 slides in each module just for check the R-dependence, 
	Double_t R_Slide[Slide+1];//edge each slide

	/*FAEC{{{*/
	if(Detector_Name == "FAEC"){
		VP_In = 3110000; //Front virtual plane
		Z_In = 413.0; //cm Front VP position
		R_Min_In = 105.0;//cm, Front VP position
		R_Max_In = 235.0;//cm, Front VP position

		VP_Out  = 3140000; //Rear virtual plane
		Z_Out = 466.0; //cm Rear VP position
		R_Min_Out = 90.0;//cm, Rear VP position
		R_Max_Out = 265.0;//cm, Rear VP position

		VP_EC = 3110000; //Front virtual plane
		Z_EC = 427.5; //cm, EC virtual plane front, VP position = 413cm
		R_Min_EC = 105.0;//cm, EC VP position
		R_Max_EC = 235.0;//cm, EC VP position
		EC_ID = 0; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*FAEC}}}*/
	
	/*MRPC{{{*/
	else if(Detector_Name == "MRPC"){
		VP_In = 4110000; //Front virtual plane
		Z_In = 408.0; //cm Front VP position
		R_Min_In = 96.0;//cm, Front VP position
		R_Max_In = 210.0;//cm, Front VP position

		VP_Out  = 3110000; //Rear virtual plane
		Z_Out = 413.0; //cm Rear VP position
		R_Min_Out = 105.0;//cm, Rear VP position
		R_Max_Out = 235.0;//cm, Rear VP position

		VP_EC = 3110000; //Front virtual plane
		Z_EC = 427.5; //cm, EC virtual plane front, VP position = 413cm
		R_Min_EC = 105.0;//cm, EC VP position
		R_Max_EC = 235.0;//cm, EC VP position
		EC_ID = 0; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*}}}*/
	
	/*FASPD{{{*/
	else if(Detector_Name == "FASPD"){
		VP_In = 5110000; //Front virtual plane
		Z_In = 406.0; //cm Front VP position
		R_Min_In = 96.0;//cm, Front VP position
		R_Max_In = 210.0;//cm, Front VP position

		VP_Out  = 4110000; //Rear virtual plane
		Z_Out = 408.0; //cm Rear VP position
		R_Min_Out = 96.0;//cm, Rear VP position
		R_Max_Out = 210.0;//cm, Rear VP position

		VP_EC = 3110000; //Front virtual plane
		Z_EC = 427.5; //cm, EC virtual plane front, VP position = 413cm
		R_Min_EC = 105.0;//cm, EC VP position
		R_Max_EC = 235.0;//cm, EC VP position
		EC_ID = 0; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*}}}*/
	
	/*HGC{{{*/
	else if(Detector_Name == "HGC"){
		VP_In = 2210000; //Front virtual plane
		Z_In = 305.0; //cm Front VP position
		R_Min_In = 83.1;//cm, Front VP position
		R_Max_In = 230.0;//cm, Front VP position

		VP_Out  = 5110000; //Rear virtual plane
		Z_Out = 406.5; //cm Rear VP position
		R_Min_Out = 96.0;//cm, Rear VP position
		R_Max_Out = 210.0;//cm, Rear VP position

		VP_EC = 3110000; //Front virtual plane
		Z_EC = 427.5; //cm, EC virtual plane front, VP position = 413cm
		R_Min_EC = 105.0;//cm, EC VP position
		R_Max_EC = 235.0;//cm, EC VP position
		EC_ID = 0; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*FAEC}}}*/
	
	/*LGC{{{*/
	else if(Detector_Name == "LGC"){
		VP_In = 2110000; //Front virtual plane
		Z_In = 96.0; //cm Front VP position
		R_Min_In = 58.1;//cm, Front VP position
		R_Max_In = 126.9;//cm, Front VP position

		VP_Out  = 2210000; //Rear virtual plane
		Z_Out = 305.0; //cm Rear VP position
		R_Min_Out = 83.0;//cm, Rear VP position
		R_Max_Out = 235.0;//cm, Rear VP position

		VP_EC = 3110000; //Front virtual plane
		Z_EC = 427.5; //cm, EC virtual plane front, VP position = 413cm
		R_Min_EC = 105.0;//cm, EC VP position
		R_Max_EC = 235.0;//cm, EC VP position
		EC_ID = 0; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*}}}*/
	
	/*LAEC{{{*/
	if(Detector_Name == "LAEC"){
		VP_In = 3210000; //Front virtual plane
		Z_In = -66.5; //cm Front VP position
		R_Min_In = 83.0;//cm, Front VP position
		R_Max_In = 140.0;//cm, Front VP position

		VP_Out  = 3240000; //Rear virtual plane
		Z_Out = -14.0; //cm Rear VP position
		R_Min_Out = 83.0;//cm, Rear VP position
		R_Max_Out = 140.0;//cm, Rear VP position

		VP_EC = 3210000; //Front virtual plane
		Z_EC = -65.0; //cm, EC virtual plane front, VP position = -66.5cm
		R_Min_EC = 83.0;//cm, EC VP position
		R_Max_EC = 140.0;//cm, EC VP position
		EC_ID = 1; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*FAEC}}}*/
	
	/*LASPD{{{*/
	if(Detector_Name == "LASPD"){
		VP_In = 5210000; //Front virtual plane
		Z_In = -67.45; //cm Front VP position
		R_Min_In = 80.0;//cm, Front VP position
		R_Max_In = 135.0;//cm, Front VP position

		VP_Out  = 3210000; //Rear virtual plane
		Z_Out = -66.5; //cm Rear VP position
		R_Min_Out = 83.0;//cm, Rear VP position
		R_Max_Out = 140.0;//cm, Rear VP position

		VP_EC = 3210000; //Front virtual plane
		Z_EC = -66.5; //cm, EC virtual plane front, VP position = 413cm
		R_Min_EC = 83.0;//cm, EC VP position
		R_Max_EC = 140.0;//cm, EC VP position
		EC_ID = 1; //0->FAEC, 1->LAEC
		EC_Threshold =1.; //GeV for EC cut, the cut could be tight if using Jin's curves

		//For Segmentation only
		Module = 30; //30 module around the circle
		Slide = 6;// put 6 slides in each module just for check the R-dependence, 
		for(Int_t i=0;i<Slide+1;i++){
			R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
		}
	}
    /*FAEC}}}*/

	/*}}}*/

	/*Pre-Define{{{*/
	
	/*Particle Definition{{{*/	
	const Int_t Electron = 11;
	const Int_t Photon = 22;
	const Int_t Proton = 2212;
	//const Int_t Pi0 = 111;
	const Int_t Pion = 211;
	//const Int_t K0 = 311;
	const Int_t Kaon = 321;
	//const Int_t Beam = 0;
	//const Int_t Neutron = 2112;
	//const Int_t Neutrino1 = 12;//Nu_e
	//const Int_t Neutrino2 = 14;//Nu_Mu
	//const Int_t Neutrino3 = 16;//Nu_Tao
	/*}}}*/
	
	const Double_t MeV2GeV = 1./1000.;
	const Double_t MM2CM = 1./10.;
	Double_t In_G_All = 0, In_E_All = 0, In_P_All = 0, Out_G_All = 0, Out_E_All = 0, Out_P_All = 0;
	Double_t In_G[Module][Slide];// Number of photons going into the device
	Double_t In_E[Module][Slide];// Number of electrons going into the device
	Double_t In_P[Module][Slide];// Number of electrons going into the device
	Double_t Out_G[Module][Slide];// Number of photons going into the device
	Double_t Out_E[Module][Slide];// Number of electrons going into the device
	Double_t Out_P[Module][Slide];// Number of electrons going into the device
	for(Int_t i=0;i<Slide;i++){//I just want to initialize everything
		EC_R[i] = trig_cut_range_R[i]; //{0,105,115,130,150,200,300};
		for(Int_t j=0;j<Module;j++){
			In_G[j][i] =0.;
			In_E[j][i] =0.;
			In_P[j][i] =0.;
			Out_G[j][i] =0.;
			Out_E[j][i] =0.;
			Out_P[j][i] =0.;
		}
	}
	Double_t Count_Raw_E =0, Count_Raw_G=0, Count_Raw_P=0;
	Double_t Count_Raw_E_Out  =0, Count_Raw_G_Out =0, Count_Raw_P_Out =0;

	Int_t Slide_ID = 0;
	Int_t Module_ID = 0;
	/*End of Pre-Define}}}*/

	/*Tree Define{{{*/
	Double_t rate, mom_gen, theta_gen,phi_gen,mom_flux,theta_flux,phi_flux,r_flux,x_flux,y_flux,z_flux,R_EC;
	Double_t px_flux, py_flux,pz_flux,vx_flux, vy_flux,vz_flux,Edep, E;
	Double_t EC_Cut_G, EC_Cut_E, EC_Cut_P;
	Double_t px_gen, py_gen, pz_gen, vx_gen, vy_gen, vz_gen;
	Int_t evn, evn_flux, ID_flux,PID_flux,ID_Pick,PID_Pick, TID_flux;

	TString pid_temp = "";
	/*Find PID{{{*/
	PID_Pick = -1000;
	if(input_filename.find("_pi0_",0) != string::npos) {
		PID_Pick = Electron; //Pi0 decays into Electron+Position
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pi0_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pi0_down";
		else
			pid_temp = "pi0";
	}
	else if(input_filename.find("_pip_",0) != string::npos) {
		PID_Pick = Pion;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pip_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pip_down";
		else
			pid_temp = "pip";
	}
	else if(input_filename.find("_pim_",0) != string::npos) {
		PID_Pick = -Pion;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "pim_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "pim_down";
		else
			pid_temp = "pim";
	}
	else if(input_filename.find("_kp_",0) != string::npos) {
		PID_Pick = Kaon;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "kp_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "kp_down";
		else
			pid_temp = "kp";
	}
	else if(input_filename.find("_km_",0) != string::npos) {
		PID_Pick = -Kaon;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "km_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "km_down";
		else
			pid_temp = "km";
	}
	else if(input_filename.find("_p_",0) != string::npos||input_filename.find("_proton_",0) != string::npos) {
		PID_Pick = Proton;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "p_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "p_down";
		else
			pid_temp = "p";
	}
	else if(input_filename.find("_eDIS_",0) != string::npos) {
		PID_Pick = Electron;
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "eDIS_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "eDIS_down";
		else	
			pid_temp = "eDIS";
	}
	else if(input_filename.find("_EM_",0) != string::npos) {
		PID_Pick = Electron;
		pid_temp = "EM";
	}
	else
		cerr<<"****ERROR, I don't understand the file name!!! *****"<<endl;
     /*}}}*/
	
	TString geant = "";
	if(input_filename.find("geant4_p96",0) != string::npos)
		geant = "_p96";
	else if(input_filename.find("geant4_p95",0) != string::npos)
		geant = "_p95";
    else
		geant = "";

	TString output_filename = Detector_Name + "_"+"background_SIDIS_He3"+"_"+pid_temp+geant+".root";
	cerr<<"--- I will save the detector info into "<<output_filename.Data()<<endl;
	
	TFile *out_file = new TFile(output_filename.Data(),"recreate");
	TTree *T = new TTree("T","A new Tree for particles going in");
	//From Header Tree
	T->Branch("rate",&rate,"rate/D");
	T->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	T->Branch("mom_gen",&mom_gen,"mom_gen/D");
	T->Branch("theta_gen",&theta_gen,"theta_gen/D");
	T->Branch("phi_gen",&phi_gen,"phi_gen/D");
	T->Branch("px_gen",&px_gen,"px_gen/D");
	T->Branch("py_gen",&py_gen,"py_gen/D");
	T->Branch("pz_gen",&pz_gen,"pz_gen/D");
	T->Branch("vx_gen",&vx_gen,"vx_gen/D");
	T->Branch("vy_gen",&vy_gen,"vy_gen/D");
	T->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	T->Branch("evn_flux",&evn_flux,"evn_flux/I");
	T->Branch("mom_flux",&mom_flux,"mom_flux/D");
	T->Branch("theta_flux",&theta_flux,"theta_flux/D");
	T->Branch("phi_flux",&phi_flux,"phi_flux/D");
	T->Branch("r_flux",&r_flux,"r_flux/D");
	T->Branch("x_flux",&x_flux,"x_flux/D");
	T->Branch("y_flux",&y_flux,"y_flux/D");
	T->Branch("z_flux",&z_flux,"z_flux/D");
	T->Branch("vx_flux",&vx_flux,"vx_flux/D");
	T->Branch("vy_flux",&vy_flux,"vy_flux/D");
	T->Branch("vz_flux",&vz_flux,"vz_flux/D");
	T->Branch("px_flux",&px_flux,"px_flux/D");
	T->Branch("py_flux",&py_flux,"py_flux/D");
	T->Branch("pz_flux",&pz_flux,"pz_flux/D");
	T->Branch("E",&E,"E/D");
	T->Branch("Edep",&Edep,"Edep/D");
	T->Branch("R_EC",&R_EC,"R_EC/D");
	T->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	T->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	T->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	T->Branch("ID_flux",&ID_flux,"ID_flux/I");
	T->Branch("PID_flux",&PID_flux,"PID_flux/I");
	T->Branch("TID_flux",&TID_flux,"TID_flux/I");
	T->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	T->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");
	
	TTree *O = new TTree("O","A new Tree for paritcles going out");
	O->Branch("rate",&rate,"rate/D");
	O->Branch("evn",&evn,"evn/I");

	//From Generate Tree
	O->Branch("mom_gen",&mom_gen,"mom_gen/D");
	O->Branch("theta_gen",&theta_gen,"theta_gen/D");
	O->Branch("phi_gen",&phi_gen,"phi_gen/D");
	O->Branch("px_gen",&px_gen,"px_gen/D");
	O->Branch("py_gen",&py_gen,"py_gen/D");
	O->Branch("pz_gen",&pz_gen,"pz_gen/D");
	O->Branch("vx_gen",&vx_gen,"vx_gen/D");
	O->Branch("vy_gen",&vy_gen,"vy_gen/D");
	O->Branch("vz_gen",&vz_gen,"vz_gen/D");

	//From Flux Tree
	O->Branch("evn_flux",&evn_flux,"evn_flux/I");
	O->Branch("mom_flux",&mom_flux,"mom_flux/D");
	O->Branch("theta_flux",&theta_flux,"theta_flux/D");
	O->Branch("phi_flux",&phi_flux,"phi_flux/D");
	O->Branch("r_flux",&r_flux,"r_flux/D");
	O->Branch("x_flux",&x_flux,"x_flux/D");
	O->Branch("y_flux",&y_flux,"y_flux/D");
	O->Branch("z_flux",&z_flux,"z_flux/D");
	O->Branch("vx_flux",&vx_flux,"vx_flux/D");
	O->Branch("vy_flux",&vy_flux,"vy_flux/D");
	O->Branch("vz_flux",&vz_flux,"vz_flux/D");
	O->Branch("px_flux",&px_flux,"px_flux/D");
	O->Branch("py_flux",&py_flux,"py_flux/D");
	O->Branch("pz_flux",&pz_flux,"pz_flux/D");
	O->Branch("E",&E,"E/D");
	O->Branch("Edep",&Edep,"Edep/D");
	O->Branch("R_EC",&R_EC,"R_EC/D");
	O->Branch("EC_Cut_G",&EC_Cut_G,"EC_Cut_G/D");
	O->Branch("EC_Cut_E",&EC_Cut_E,"EC_Cut_E/D");
	O->Branch("EC_Cut_P",&EC_Cut_P,"EC_Cut_P/D");
	O->Branch("ID_flux",&ID_flux,"ID_flux/I");
	O->Branch("PID_flux",&PID_flux,"PID_flux/I");
	O->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	O->Branch("PID_Pick",&PID_Pick,"PID_Pick/I");

	const double Bin_Size = 0.1;//cm
	const int R_Bin_In = (int) ((R_Max_In - R_Min_In)/Bin_Size);//0.1 cm per bin
	TH1D *hR_In_G = new TH1D("hR_In_G","Photon Rate (Hz/cm^{2}) vs R before the detector", R_Bin_In, R_Min_In, R_Max_In);
	hR_In_G->SetXTitle("R_{flux} (cm)");
	hR_In_G->GetXaxis()->CenterTitle(1);

	TH1D *hR_In_E = new TH1D("hR_In_E","Electron Rate (Hz/cm^{2}) vs R before the detector", R_Bin_In, R_Min_In, R_Max_In);
	hR_In_E->SetXTitle("R_{flux} (cm)");
	hR_In_E->GetXaxis()->CenterTitle(1);

	TH1D *hR_In_P = new TH1D("hR_In_P","Pick particle Rate (Hz/cm^{2}) vs R before the detector", R_Bin_In, R_Min_In, R_Max_In);
	hR_In_P->SetXTitle("R_{flux} (cm)");
	hR_In_P->GetXaxis()->CenterTitle(1);

	const int R_Bin_Out = (int) ((R_Max_Out - R_Min_Out)/Bin_Size);//0.1 cm per bin
	TH1D *hR_Out_G = new TH1D("hR_Out_G","Photon Rate (Hz/cm^{2}) vs R before the detector", R_Bin_Out, R_Min_Out, R_Max_Out);
	hR_Out_G->SetXTitle("R_{flux} (cm)");
	hR_Out_G->GetXaxis()->CenterTitle(1);

	TH1D *hR_Out_E = new TH1D("hR_Out_E","Electron Rate (Hz/cm^{2}) vs R before the detector", R_Bin_Out, R_Min_Out, R_Max_Out);
	hR_Out_E->SetXTitle("R_{flux} (cm)");
	hR_Out_E->GetXaxis()->CenterTitle(1);

	TH1D *hR_Out_P = new TH1D("hR_Out_P","Picked particle Rate (Hz/cm^{2}) vs R before the detector", R_Bin_Out, R_Min_Out, R_Max_Out);
	hR_Out_P->SetXTitle("R_{flux} (cm)");
	hR_Out_P->GetXaxis()->CenterTitle(1);
	/*}}}*/

	/*Read in each event{{{*/
	Int_t FirstOne1 = 0; //incoming photons without cut
	Int_t FirstOne2 = 0; //incoming electrons without cut
	Int_t FirstOne3 = 0; //incoming picked particle without cut
	Int_t FirstOne4 = 0; //incoming photons witht R-Cut
	Int_t FirstOne5 = 0; //incoming electrons witht R-cut
	Int_t FirstOne6 = 0; //incoming picked particle with R-cut

	Int_t FirstOne7 = 0; //outgoing photons without cut
	Int_t FirstOne8 = 0; //outgoing electrons without cut
	Int_t FirstOne9 = 0; //outgoing picked particle without cut
	Int_t FirstOne10 = 0;//outgoing photons witht R-Cut
	Int_t FirstOne11= 0; //outgoing electrons witht R-cut
	Int_t FirstOne12= 0; //outgoing picked particle with R-cut

	cerr<<"++++++++++++++++ "<<endl;
	bool bCheck_Rear = kTRUE;
	Int_t nselected = nevent;
	for(Int_t i=0;i<nselected;i++){
		cout<<i<<"\r";
		header->GetEntry(i);
		if(input_filename.find("_EM_",0) != string::npos) 
			rate = (((1.5e-5)/(1.6e-19))/nevent); //Count to KHz for 15uA electron events;
		else
			rate = head_rate->at(0)/filenum;
		evn = head_evn->at(0);

		generated->GetEntry(i);
		const Int_t ng = gen_pid->size();//Normally there is only one particle in the gen
		Double_t gen_theta_array[ng], gen_phi_array[ng], gen_mom_array[ng];
		Double_t gen_px_array[ng],gen_py_array[ng],gen_pz_array[ng],gen_vx_array[ng],gen_vy_array[ng],gen_vz_array[ng];
		Int_t Is_Pick = -1;
		for(Int_t ig=0;ig<ng;ig++){
			if((Int_t)gen_pid->at(ig)==PID_Pick)
				Is_Pick = ig;
			gen_mom_array[ig] = sqrt( pow(gen_px->at(ig),2)+pow(gen_py->at(ig),2)+pow(gen_pz->at(ig),2) ); //GeV
			gen_theta_array[ig] = acos(gen_pz->at(ig)/gen_mom_array[ig])*DEG;//Degree
			gen_phi_array[ig] = atan2( gen_py->at(ig), gen_px->at(ig))*DEG;//Degree
			gen_mom_array[ig] *= MeV2GeV;
			gen_px_array[ig] = gen_px->at(ig)*MeV2GeV; //GeV
			gen_py_array[ig] = gen_py->at(ig)*MeV2GeV; //GeV
			gen_pz_array[ig] = gen_pz->at(ig)*MeV2GeV; //GeV
			gen_vx_array[ig] = gen_vx->at(ig)*MM2CM; //cm
			gen_vy_array[ig] = gen_vy->at(ig)*MM2CM; //cm
			gen_vz_array[ig] = gen_vz->at(ig)*MM2CM; //cm
			//cerr<<Form("---#%d@%d Px=%f Py=%f Pz=%f P=%f Theta =%f Phi=%f ",i, ig,gen_px->at(ig),gen_py->at(ig),gen_pz->at(ig),gen_mom_array[ig],gen_theta_array[ig],gen_phi_array[ig])<<endl;
		}
		mom_gen = gen_mom_array[Is_Pick];theta_gen = gen_theta_array[Is_Pick];phi_gen = gen_phi_array[Is_Pick];
		px_gen = gen_px_array[Is_Pick];py_gen = gen_py_array[Is_Pick];pz_gen = gen_pz_array[Is_Pick];
		vx_gen = gen_vx_array[Is_Pick];vy_gen = gen_vy_array[Is_Pick];vz_gen = gen_vz_array[Is_Pick];

		if(((pid_temp.Contains("eDIS"))&&theta_gen>=8.0)||!(pid_temp.Contains("eDIS"))){
		//if((!(pid_temp.Contains("EM"))&&theta_gen>=8.0)||(pid_temp.Contains("EM"))){
		//if(1){
			flux->GetEntry(i);
		
			/*Into the Device{{{*/
			FirstOne1 = 0; FirstOne2 = 0; FirstOne3 = 0;
			FirstOne4 = 0; FirstOne5 = 0; FirstOne6 = 0; 
			ID_Pick = VP_In;
			for (unsigned int j=0;j<flux_hitn->size();j++) {
				evn_flux = j;
				x_flux = flux_avg_x->at(j)*MM2CM; y_flux = flux_avg_y->at(j)*MM2CM; z_flux = flux_avg_z->at(j)*MM2CM;
				vx_flux = flux_vx->at(j)*MM2CM;   vy_flux = flux_vy->at(j)*MM2CM;	vz_flux = flux_vz->at(j)*MM2CM;
				px_flux = flux_px->at(j)*MeV2GeV; py_flux = flux_py->at(j)*MeV2GeV; pz_flux = flux_pz->at(j)*MeV2GeV;
				Edep = flux_totEdep->at(j)*MeV2GeV; E = flux_trackE->at(j)*MeV2GeV;

				r_flux=sqrt(pow(x_flux,2)+pow(y_flux,2));//cm
				mom_flux=sqrt(pow(px_flux,2)+pow(py_flux,2)+pow(pz_flux,2));//GeV
				theta_flux = atan(sqrt(pow(px_flux,2)+pow(py_flux,2))/(pz_flux))*DEG;//Degree
				phi_flux = fabs(atan(y_flux/x_flux))*DEG;
				if(y_flux>0 && x_flux<0) phi_flux += 90.0;
				else if(y_flux<0 && x_flux<0) phi_flux += 180.0;
				else if(y_flux<0 && x_flux>0) phi_flux += 270.0;

				Double_t Delta_Z = Z_EC - Z_In;
				Double_t x_EC = x_flux + Delta_Z * px_flux/pz_flux;
				Double_t y_EC = y_flux + Delta_Z * py_flux/pz_flux;
				R_EC = sqrt(x_EC*x_EC+y_EC*y_EC);
				//			R_EC = abs((Z_In - Z_EC) * tan(theta_flux/DEG) + r_flux);//cm

				ID_flux = (Int_t) (flux_id->at(j));
				PID_flux = (Int_t) (flux_pid->at(j));
				TID_flux = (Int_t) (flux_tid->at(j));
				if(PID_flux==-Electron) PID_flux = Electron;//Force e+/- to be counted at one time

				//Select the right Events
				if(ID_flux!=ID_Pick) continue; //Only look at events on the detector plane (VP here) 
				if(PID_flux!=Electron&&PID_flux!=Proton&&PID_flux!=Photon
						&&PID_flux!=Pion&&PID_flux!=-Pion
						&&PID_flux!=Kaon&&PID_flux!=-Kaon) continue;
				if(r_flux > R_Max_In||r_flux<R_Min_In) continue;//Must be within the radius range
				//if(flux_pz->at(j)<1e-9)continue;//Cout out backward particles
				if(mom_flux<1e-9) continue; //Cut out Zero-E particles

				//Photon going in  
				if(PID_flux==Photon){
					if(FirstOne1<1)
						Count_Raw_G+=1.0*rate;
					FirstOne1 ++;;
				}
				//Electron going in  
				if(PID_flux==Electron){
					if(FirstOne2<1)
						Count_Raw_E+=1.0*rate;
					FirstOne2 ++;;
				}
				//Pick Particle going in  
				if(PID_flux==PID_Pick){
					if(FirstOne3<1)
						Count_Raw_P+=1.0*rate;
					FirstOne3 ++;;
				}

				if(R_EC>R_Max_EC || R_EC<R_Min_EC) continue;//Project the particles to the EC and see whether it has been accepted there

				//Selct the right Cut	
				//cut[Det_ID][Cut_ID][Data_Point][Cut_Info]: Cut_Info: R_Min, R_Max, P_Min, P_Max, e_Eff, pi_Eff
				EC_Cut_G = 0;EC_Cut_E = 0;EC_Cut_P = 0;
				if(mom_flux>=0.5*EC_Threshold){//Charged Particles going out 
					/*Look at the EC to check E-Cut{{{*/ 
					Int_t EC_PID = 0;
					for (unsigned int m=0;m<flux_hitn->size();m++) {
						Int_t ID_flux_m = (Int_t) (flux_id->at(m));
						Int_t PID_flux_m = (Int_t) (flux_pid->at(m));
						Int_t TID_flux_m = (Int_t) (flux_tid->at(m));
						if(ID_flux_m==VP_EC){//On EC VP
							if(PID_flux_m==Proton||PID_flux_m==Pion||PID_flux_m==Kaon||PID_flux_m==-Pion||PID_flux_m==-Kaon)
								EC_PID = 5;//Temperately use Pion's EC-Cut for Proton and Kaons
							else if(PID_flux_m==Electron||PID_flux_m==-Electron||PID_flux_m==Photon) 
								EC_PID = 4;
							else
								continue;

							Double_t EC_Cut = 0.0;
							Double_t r_EC=sqrt(pow(flux_avg_x->at(m),2)+pow(flux_avg_y->at(m),2))*MM2CM;//cm
							if(r_EC>R_Max_EC||r_EC<R_Min_EC) continue;//The radius of a mrpc sector is 210cm;

							if(flux_pz->at(m)<1e-9)continue;//Cut out backward particles
							Double_t fmom_EC=sqrt(pow(flux_px->at(m),2)+pow(flux_py->at(m),2)+pow(flux_pz->at(m),2))*MeV2GeV;//GeV
							if(fmom_EC<1e-9) continue; //Cut out Zero-E particles

							/*Find E-Cut{{{*/
							for(Int_t k=0;k<EC_Trigger_Slide;k++){
								for(Int_t l=0;l<EC_Trigger_Mom_Bin;l++){
									if(r_EC>trig_cut[EC_ID][k][l][0]&&r_EC<=trig_cut[EC_ID][k][l][1]){
										if(fmom_EC>trig_cut[EC_ID][k][l][2]&&fmom_EC<=trig_cut[EC_ID][k][l][3]){
											EC_Cut =trig_cut[EC_ID][k][l][EC_PID];
										}
									}
								}
							}
							/*}}}*/
							if(PID_flux_m==Photon && EC_Cut>EC_Cut_G)
								EC_Cut_G = EC_Cut;
							if(PID_flux_m==Electron && EC_Cut>EC_Cut_E)
								EC_Cut_E = EC_Cut;
							if(PID_flux_m==PID_Pick && EC_Cut>EC_Cut_P)
								EC_Cut_P = EC_Cut;
						}
					}//for (Int_t m=0;m<flux_hitn->size();m++) 

					if((EC_Cut_G<-1e-9||EC_Cut_G>1)&&(EC_Cut_E<-1e-9||EC_Cut_E>1)&&(EC_Cut_P<-1e-9||EC_Cut_P>1)){
						cerr<<"----In: I can't find the cut!"<<Form(" --- PID=%d, EC = %d R= %f,  E =%f Cut_G = %f Cut_E = %f Cut_P = %f", EC_PID, EC_ID, r_flux, mom_flux, EC_Cut_G, EC_Cut_E, EC_Cut_P)<<endl;
						return;
					}
					/*}}}*/
				}//if(mom_flux>=0.5*EC_Threshold)

				if(EC_Cut_E<EC_Cut_G)
					EC_Cut_E= EC_Cut_G;

				Module_ID = (Int_t) phi_flux/(360./Module);
				for(Int_t k=0;k<Slide+1;k++){
					if(r_flux>R_Slide[k]&&r_flux<R_Slide[k+1])
						Slide_ID = k;
				}

				if(PID_flux==Photon&& mom_flux>=0.95*EC_Threshold){//#Photons trigger 
					if(FirstOne4<1){
						In_G[Module_ID][Slide_ID]+=EC_Cut_G*rate;
						hR_In_G->Fill(r_flux, EC_Cut_G*rate);
					}
					FirstOne4++;
				}
				if(PID_flux==Electron && mom_flux>=0.95*EC_Threshold){//#Electron trigger 
					if(FirstOne5<1){
						In_E[Module_ID][Slide_ID]+=EC_Cut_E*rate;
						hR_In_E->Fill(r_flux, EC_Cut_E*rate);
					}
					FirstOne5 ++;
				}
				
				if(PID_flux==Proton)	
					EC_Cut_P *= 0.5; //Zhiwen assumes Proton has 50% lower efficiences than pions

				if(PID_flux==PID_Pick && mom_flux>=0.95*EC_Threshold){//#Charged Particles trigger
					if(FirstOne6<1){
						In_P[Module_ID][Slide_ID]+=EC_Cut_P*rate;
						hR_In_P->Fill(r_flux, EC_Cut_P*rate);
					}
					FirstOne6 ++;
				}

				T->Fill();
			}
			/*}}}*/

			/*Out from Device{{{*/
			FirstOne7 = 0; FirstOne8 = 0; FirstOne9 = 0;
			FirstOne10= 0; FirstOne11= 0; FirstOne12= 0; 
			ID_Pick = VP_Out;
			for (unsigned int j=0;j<flux_hitn->size();j++) {
				evn_flux = j;
				x_flux = flux_avg_x->at(j)*MM2CM; y_flux = flux_avg_y->at(j)*MM2CM; z_flux = flux_avg_z->at(j)*MM2CM;
				vx_flux = flux_vx->at(j)*MM2CM;   vy_flux = flux_vy->at(j)*MM2CM;	vz_flux = flux_vz->at(j)*MM2CM;
				px_flux = flux_px->at(j)*MeV2GeV; py_flux = flux_py->at(j)*MeV2GeV; pz_flux = flux_pz->at(j)*MeV2GeV;
				Edep = flux_totEdep->at(j)*MeV2GeV; E = flux_trackE->at(j)*MeV2GeV;

				r_flux=sqrt(pow(x_flux,2)+pow(y_flux,2));//cm
				mom_flux=sqrt(pow(px_flux,2)+pow(py_flux,2)+pow(pz_flux,2));//GeV
				theta_flux = atan(sqrt(pow(px_flux,2)+pow(py_flux,2))/(pz_flux))*DEG;//Degree
				phi_flux = fabs(atan(y_flux/x_flux))*DEG;
				if(y_flux>0 && x_flux<0) phi_flux += 90.0;
				else if(y_flux<0 && x_flux<0) phi_flux += 180.0;
				else if(y_flux<0 && x_flux>0) phi_flux += 270.0;

				Double_t Delta_Z = Z_EC - Z_Out;
				Double_t x_EC = x_flux + Delta_Z * px_flux/pz_flux;
				Double_t y_EC = y_flux + Delta_Z * py_flux/pz_flux;
				R_EC = sqrt(x_EC*x_EC+y_EC*y_EC);
				//			R_EC = abs((Z_In - Z_EC) * tan(theta_flux/DEG) + r_flux);//cm

				ID_flux = (Int_t) (flux_id->at(j));
				PID_flux = (Int_t) (flux_pid->at(j));
				if(PID_flux==-Electron) PID_flux = Electron;//Force e+/- to be counted at one time

				//Select the right Events
				if(ID_flux!=ID_Pick) continue; //Only look at events on the detector plane (VP here) 
				if(PID_flux!=Electron&&PID_flux!=Proton&&PID_flux!=Photon
						&&PID_flux!=Pion&&PID_flux!=-Pion
						&&PID_flux!=Kaon&&PID_flux!=-Kaon) continue;
				if(r_flux > R_Max_Out||r_flux<R_Min_Out) continue;//Must be within the radius range
				//if(flux_pz->at(j)<1e-9)continue;//Cout out backward particles
				if(mom_flux<1e-9) continue; //Cut out Zero-E particles

				//Photon going in  
				if(PID_flux==Photon){
					if(FirstOne7<1)
						Count_Raw_G_Out+=1.0*rate;
					FirstOne7 ++;;
				}
				//Electron going in  
				if(PID_flux==Electron){
					if(FirstOne8<1)
						Count_Raw_E_Out+=1.0*rate;
					FirstOne8 ++;;
				}
				//Pick Particle going in  
				if(PID_flux==PID_Pick){
					if(FirstOne9<1)
						Count_Raw_P_Out+=1.0*rate;
					FirstOne9 ++;;
				}

				if(R_EC>R_Max_EC || R_EC<R_Min_EC) continue;//Project the particles to the EC and see whether it has been accepted there

				//Selct the right Cut	
				//cut[Det_ID][Cut_ID][Data_Point][Cut_Info]: Cut_Info: R_Min, R_Max, P_Min, P_Max, e_Eff, pi_Eff
				EC_Cut_G = 0;EC_Cut_E = 0;EC_Cut_P = 0;
				if(mom_flux>=0.5*EC_Threshold){//Charged Particles going out 
					/*Look at the EC to check E-Cut{{{*/ 
					Int_t EC_PID = 0;
					for (unsigned int m=0;m<flux_hitn->size();m++) {
						Int_t ID_flux_m = (Int_t) (flux_id->at(m));
						Int_t PID_flux_m = (Int_t) (flux_pid->at(m));
						if(ID_flux_m==VP_EC){//On EC VP
							if(PID_flux_m==Proton||PID_flux_m==Pion||PID_flux_m==Kaon||PID_flux_m==-Pion||PID_flux_m==-Kaon)
								EC_PID = 5;//Temperately use Pion's EC-Cut for Proton and Kaons
							else if(PID_flux_m==Electron||PID_flux_m==-Electron||PID_flux_m==Photon) 
								EC_PID = 4;
							else
								continue;

							Double_t EC_Cut = 0.0;
							Double_t r_EC=sqrt(pow(flux_avg_x->at(m),2)+pow(flux_avg_y->at(m),2))*MM2CM;//cm
							if(r_EC>R_Max_EC||r_EC<R_Min_EC) continue;//The radius of a mrpc sector is 210cm;

							if(flux_pz->at(m)<1e-9)continue;//Cut out backward particles
							Double_t fmom_EC=sqrt(pow(flux_px->at(m),2)+pow(flux_py->at(m),2)+pow(flux_pz->at(m),2))*MeV2GeV;//GeV
							if(fmom_EC<1e-9) continue; //Cut out Zero-E particles

							/*Find E-Cut{{{*/
							for(Int_t k=0;k<EC_Trigger_Slide;k++){
								for(Int_t l=0;l<EC_Trigger_Mom_Bin;l++){
									if(r_EC>trig_cut[EC_ID][k][l][0]&&r_EC<=trig_cut[EC_ID][k][l][1]){
										if(fmom_EC>trig_cut[EC_ID][k][l][2]&&fmom_EC<=trig_cut[EC_ID][k][l][3]){
											EC_Cut =trig_cut[EC_ID][k][l][EC_PID]; 
										}
									}
								}
							}
							/*}}}*/
							if(PID_flux_m==Photon && EC_Cut>EC_Cut_G)
								EC_Cut_G = EC_Cut;
							if(PID_flux_m==Electron && EC_Cut>EC_Cut_E)
								EC_Cut_E = EC_Cut;
							if(PID_flux_m==PID_Pick && EC_Cut>EC_Cut_P)
								EC_Cut_P = EC_Cut;
						}
					}//for (Int_t m=0;m<flux_hitn->size();m++) 

					if((EC_Cut_G<-1e-9||EC_Cut_G>1)&&(EC_Cut_E<-1e-9||EC_Cut_E>1)&&(EC_Cut_P<-1e-9||EC_Cut_P>1)){
						cerr<<"----In: I can't find the cut!"<<Form(" --- PID=%d, EC = %d R= %f,  E =%f Cut_G = %f Cut_E = %f Cut_P = %f", EC_PID, EC_ID, r_flux, mom_flux, EC_Cut_G, EC_Cut_E, EC_Cut_P)<<endl;
						return;
					}
					/*}}}*/
				}//if(mom_flux>=0.5*EC_Threshold)

				Module_ID = (Int_t) phi_flux/(360./Module);
				for(Int_t k=0;k<Slide+1;k++){
					if(r_flux>R_Slide[k]&&r_flux<R_Slide[k+1])
						Slide_ID = k;
				}

				if(PID_flux==Photon&& mom_flux>=0.95*EC_Threshold){//#Photons trigger 
					if(FirstOne10<1){
						Out_G[Module_ID][Slide_ID]+=EC_Cut_G*rate;
						hR_Out_G->Fill(r_flux, EC_Cut_G*rate);
					}
					FirstOne10++;
				}
				if(PID_flux==Electron && mom_flux>=0.95*EC_Threshold){//#Electron trigger 
					if(FirstOne11<1){
						Out_E[Module_ID][Slide_ID]+=EC_Cut_E*rate;
						hR_Out_E->Fill(r_flux, EC_Cut_E*rate);
					}
					FirstOne11++;;
				}
	
				if(PID_flux==Proton)	
					EC_Cut_P *= 0.5; //Zhiwen assumes Proton has 50% lower efficiences than pions

				if(PID_flux==PID_Pick && mom_flux>=0.95*EC_Threshold){//#Charged Particles trigger
					if(FirstOne12<1){
						Out_P[Module_ID][Slide_ID]+=EC_Cut_P*rate;
						hR_Out_P->Fill(r_flux, EC_Cut_P*rate);
					}
					FirstOne12++;;
				}

				O->Fill();
			}
			/*}}}*/

		}
	}
		/*End Read in each event}}}*/

	Double_t Count_To_Rate = 1.0/1e3;

	hR_In_G->Write();  hR_In_E->Write();  hR_In_P->Write();
	hR_Out_G->Write(); hR_Out_E->Write(); hR_Out_P->Write();
	T->Write(); out_file->Close(); file->Close();
	
	/*Output Results{{{*/
	cerr<<" --- Using Converting Factor = "<<Count_To_Rate<<endl;

	TString input_out = output_filename;
	input_out.ReplaceAll(".root",".out");
	ofstream outfile(input_out.Data());
	for(Int_t k=0;k<Slide;k++){
		for(Int_t l=0;l<Module;l++){
			In_G_All +=In_G[l][k];
			In_E_All +=In_E[l][k];
			In_P_All +=In_P[l][k];
			Out_G_All +=Out_G[l][k];
			Out_E_All +=Out_E[l][k];
			Out_P_All +=Out_P[l][k];
		}
	}
	cerr<<" ======  In_G = "<< Count_Raw_G*Count_To_Rate <<", after Cut =  "<< In_G_All*Count_To_Rate<<endl;
	cerr<<" ======  In_E = "<< Count_Raw_E*Count_To_Rate <<", after Cut =  "<< In_E_All*Count_To_Rate<<endl;
	cerr<<" ======  In_P = "<< Count_Raw_P*Count_To_Rate <<", after Cut =  "<< In_P_All*Count_To_Rate<<endl;
	cerr<<" ====== Out_G = "<< Count_Raw_G_Out*Count_To_Rate <<", after Cut =  "<< Out_G_All*Count_To_Rate<<endl;
	cerr<<" ====== Out_E = "<< Count_Raw_E_Out*Count_To_Rate <<", after Cut =  "<< Out_E_All*Count_To_Rate<<endl;
	cerr<<" ====== Out_P = "<< Count_Raw_P_Out*Count_To_Rate <<", after Cut =  "<< Out_P_All*Count_To_Rate<<endl;

	outfile<<" ======  In_G = "<< Count_Raw_G*Count_To_Rate <<", after Cut =  "<< In_G_All*Count_To_Rate<<endl;
	outfile<<" ======  In_E = "<< Count_Raw_E*Count_To_Rate <<", after Cut =  "<< In_E_All*Count_To_Rate<<endl;
	outfile<<" ======  In_P = "<< Count_Raw_P*Count_To_Rate <<", after Cut =  "<< In_P_All*Count_To_Rate<<endl;
	outfile<<" ====== Out_G = "<< Count_Raw_G_Out*Count_To_Rate <<", after Cut =  "<< Out_G_All*Count_To_Rate<<endl;
	outfile<<" ====== Out_E = "<< Count_Raw_E_Out*Count_To_Rate <<", after Cut =  "<< Out_E_All*Count_To_Rate<<endl;
	outfile<<" ====== Out_P = "<< Count_Raw_P_Out*Count_To_Rate <<", after Cut =  "<< Out_P_All*Count_To_Rate<<endl;
	/*}}}*/

}

void RunAll(TString detector, TString runlist){
   if(detector!="LASPD"&&detector!="LAEC"&&detector!="LGC"&&detector!="HGC"&&detector!="FASPD"&&detector!="MRPC"&&detector!="FAEC"){
	   cerr<<"--- The detector list is: LASPD, LAEC, LGC, HGC, FASPD, MRPC, FAEC"<<endl;
	   cerr<<"    But you give "<<detector.Data()<<endl;
}
   else{
	   ifstream file(runlist.Data());
       TString file_name;
       while(!file.eof()){
		   file >> file_name;
		   GetRate(detector.Data(), file_name.Data());
	   }
	   file.close();
   }

}
