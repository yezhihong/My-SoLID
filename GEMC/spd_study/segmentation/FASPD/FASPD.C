////////////////////////////////////////////
// This script is modified to read data
// with mpid option turned on during GEMC
// -- Z. Ye, 01/06/2015

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
#include <TRandom.h>
/*}}}*/

using namespace std;

TString FASPD(string input_filename)
{
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	const double DEG=180./3.1415926;

	TFile *file=new TFile(input_filename.c_str());
	if (file->IsZombie()) {
		cout << "*** Error opening file" << input_filename << endl;
		//continue;
		exit(-1);
	}
	else cout << " ==> open file " << input_filename << endl;
	
	/*Temp Rate Normalization{{{*/
	double filenum=1.;

	if (input_filename.find("_EM_",0) != string::npos) {
		cout << "==> EM background from beam on target" <<  endl;
	}
	else if (input_filename.find("background_p",0) != string::npos) {
		cout << "==> background from normalized event generator with all interaction" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "==> filenum=" << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "==> we need filenum for addtional normalization" << endl; return "0";}
	}
	else {cout << "==> not EM or clean_weighted or dirty_normalized " << endl; return "0";}

	/*
	else if (input_filename.find("_clean_weighted_",0) != string::npos) {
		cout << "==> background from weighted event generator with no interaction except decay" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "==> filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "==> we need filenum for addtional normalization" << endl; return "0";}		
	}	  
	else if (input_filename.find("_dirty_normalized_",0) != string::npos) {
		cout << "==> background from normalized event generator with all interaction" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "==> filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "==> we need filenum for addtional normalization" << endl; return "0";}
	}	  	
	else {cout << "==> not EM or clean_weighted or dirty_normalized " << endl; return "0";}
	*/
    /*}}}*/

	/*Set Branch{{{*/
	//Header Tree:
	// Var#1~#8 are free slots for propogating important info from the "INPUT generator seed"
	// For example, they can be used to store the cross section and other physics quantities
	// In eicRate, we store the following quantities:
	// var1->Wprate, var2->Wmrate, var3->targetPol, var4->x,var5->y, var6->W, var7->Q2, var8->rate 
	//
	TTree *header = (TTree*) file->Get("header");
	vector <double> *head_evn=0,*head_evn_type=0; //Note: Vectors have to be initialized at first!!!
	vector <double> *head_beamPol=0;
	vector<double> *head_Wmrate=0, *head_Wprate=0, *head_targetPol=0, *head_x=0, *head_Q2=0, *head_W=0, *head_rate=0, *head_y=0;
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
	vector<double> *flux_id=0,*flux_hitn=0,*flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
	vector<double> *flux_trackE=0,*flux_totEdep=0;
	vector<double> *flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0;
	vector<double> *flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
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
	/*End Set Branch}}}*/
	int nevent = (int)generated->GetEntries();
	cout << "nevent = " << nevent << endl;

	/*EC Electron Trigger{{{*/
	const int Ntrigline=6,Ntriglinebin=21;
	int region_index;
	//	if (input_filename.Contains("SIDIS_FA")) region_index=0;
	//	else if (input_filename.Contains("SIDIS_LA")) region_index=1;
	//	else {cout << "need option for FA or LA region" << endl; exit(-1);}
	region_index = 0;

	//int det[2]={8,12};  //detecor ID
	//double rmin[2]={105,80};
	//double rmax[2]={235,140};

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
	for (int i=0;i<Ntrigline;i++){
		file_trig_cut[1][i][0]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Large_RunElectron_GetEfficienciesBackGround_Oct2013_SIDIS_Full_bgd_TrigSH2.0.root",Trigger_Dir_1GeVCut.Data()));
		file_trig_cut[1][i][1]=new TFile(Form("%s/Lead2X0PbBlock_Hex.1.SIDIS_Large_RunPion_GetEfficienciesBackGround_Oct2013_SIDIS_Full_bgd_TrigSH2.0.root",Trigger_Dir_1GeVCut.Data()));
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
	
	const int EC_Trigger_Slide = Ntrigline;// put 6 slides in each module just for check the R-dependence, 
	const int EC_Trigger_Mom_Bin = Ntriglinebin+2;
	double EC_R[EC_Trigger_Slide];//Center location of each slide
	for(int i=0;i<EC_Trigger_Slide;i++){
		EC_R[i] = trig_cut_range_R[i]; //{0,105,115,130,150,200,300};
	}

	for(int i=0;i<Ntrigline;i++){
	  file_trig_cut[0][i][0]->Close();
	  file_trig_cut[0][i][1]->Close();
	  file_trig_cut[1][i][0]->Close();
	  file_trig_cut[1][i][1]->Close();
	}
	/*}}}*/

	/*Define FASPD (3mm) {{{*/
    TString Detector_Name = "FASPD";   
	const int VP_ID = 5110000; //Front virtual plane
	const double Z_In = 406.50; //cm Front VP position
	const double R_Min_In = 96.0;//cm, Front VP position
	const double R_Max_In = 210.0;//cm, Front VP position
	const int SPD_ID = 5100000; //SPD

    //MRPC virtual plane	
	const int VP_Out  = 4110000; //Rear virtual plane
	const double Z_Out = 408.0; //cm Rear VP position
	const double R_Min_Out = 96.0;//cm, Rear VP position
	const double R_Max_Out = 210.0;//cm, Rear VP position

	const int VP_EC = 3110000; //Front virtual plane
	const double Z_EC = 427.5; //cm, EC virtual plane front, VP position = 413cm
	const double R_Min_EC = 105.0;//cm, EC VP position
	const double R_Max_EC = 235.0;//cm, EC VP position
	double EC_Threshold = 1.; //GeV for EC cut, the cut could be tight if using Jin's curves
    const int EC_ID = 0; //0->FAEC, 1->LAEC

	//For Segmentation only
	const int Module = 30; //30 module around the circle
	const int Slide = 3;// put 3 slides in each module just for check the R-dependence, 
	double R_Slide[Slide+1];//edge each slide
	for(int i=0;i<Slide+1;i++){
		R_Slide[i] = i*(R_Max_In-R_Min_In)/Slide;
	}
		
	//Other Definition
	const double MeV2GeV = 1./1000.;
	const double MM2CM = 1./10.;	
	const int Electron = 11;
	const int Photon = 22;
	const int Proton = 2212;
	//const int Pi0 = 111;
	const int Pion = 211;
	//const int K0 = 311;
	const int Kaon = 321;
	//const int Neutron = 2112;
	//const int Beam = 0;
	//const int Neutron = 2112;
	//const int Neutrino1 = 12;//Nu_e
	//const int Neutrino2 = 14;//Nu_Mu
	//const int Neutrino3 = 16;//Nu_Tao
	/*}}}*/

	/*Tree Define{{{*/
	TString pid_temp = "";
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
	else if(input_filename.find("_EM_",0) != string::npos) {
		pid_temp = "EM";
	}
	else if(input_filename.find("_eDIS_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "eDIS_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "eDIS_down";
		else
			pid_temp = "eDIS";
	}
	else if(input_filename.find("_ePB_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "ePB_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "ePB_down";
		else
			pid_temp = "ePB";
	}
	else if(input_filename.find("_p_",0) != string::npos) {
		if(input_filename.find("_upstream_",0) != string::npos) 
			pid_temp = "p_up";
		else if(input_filename.find("_downstream_",0) != string::npos) 
			pid_temp = "p_down";
		else
			pid_temp = "p";
	}
	else
		cerr<<"****ERROR, I don't understand the file name!!! *****"<<endl;

	TString output_filename = Detector_Name + "_"+"background_SIDIS_He3"+"_"+pid_temp+"_mpid.root";
	cerr<<"--- I will save the detector info into "<<output_filename.Data()<<endl;
	double rate, mom_gen, theta_gen,phi_gen,mom_flux,theta_flux,phi_flux,r_flux,x_flux,y_flux,z_flux,R_EC,EC_Cut,EC_Cut_Max;
	double px_flux, py_flux,pz_flux,vx_flux, vy_flux,vz_flux, E,EDep, convert;
	double px_gen, py_gen, pz_gen, vx_gen, vy_gen, vz_gen;
	int evn, evn_flux,nhit_flux, ID_flux,PID_flux,TID_flux,ID_Pick;
	int Slide_ID, Module_ID;
	
	TFile *out_file = new TFile(output_filename.Data(),"recreate");
	TTree *T = new TTree("T","A new Tree");
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
	T->Branch("nhit_flux",&nhit_flux,"nhit_flux/I");
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
	T->Branch("R_EC",&R_EC,"R_EC/D");
	T->Branch("EDep",&EDep,"EDep/D");
	T->Branch("E",&E,"E/D");
	T->Branch("convert",&convert,"convert/D");
	T->Branch("EC_Cut",&EC_Cut,"EC_Cut/D");
	T->Branch("EC_Cut_Max",&EC_Cut_Max,"EC_Cut_Max/D");
	T->Branch("ID_flux",&ID_flux,"ID_flux/I");
	T->Branch("PID_flux",&PID_flux,"PID_flux/I");
	T->Branch("TID_flux",&TID_flux,"TID_flux/I");
	T->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	T->Branch("Module", &Module_ID,"Module/I");
	T->Branch("Slide", &Slide_ID,"Slide/I");
	/*}}}*/

	/*Histogram Define{{{*/
    TH2F *hist_crg = new TH2F("hist_crg","Triggerred Electrons on FASPD",200,-200,200,200,-200,200);
	hist_crg->SetXTitle("x");
	hist_crg->SetYTitle("y");

	TH2F *hist_pho = new TH2F("hist_pho","Triggerred Photon on FASPD",200,-200,200,200,-200,200);
	hist_pho->SetXTitle("x");
	hist_pho->SetYTitle("y");

	/*}}}*/
	
	/*Read in each event{{{*/
	cerr<<"++++++++++++++++ "<<endl;
	Int_t nselected = nevent;
	double In_G = 0.0, In_E = 0.0;
	int FirstOne0 = 0, FirstOne1 = 0;
    TRandom *tr = new TRandom(); 

	for(Int_t i=0;i<nselected;i++){
		cout<<i<<"\r";
        header->GetEntry(i);
		if(input_filename.find("_EM_",0) != string::npos) 
			rate = ((1.5e-5)/(1.6e-19))/nevent;
		else 
			rate = head_rate->at(0)/filenum;
		evn = head_evn->at(0);

		generated->GetEntry(i);
		const int ng = gen_pid->size();//Normally there is only one particle in the gen
		double gen_theta_array[ng], gen_phi_array[ng], gen_mom_array[ng];
		double gen_px_array[ng],gen_py_array[ng],gen_pz_array[ng],gen_vx_array[ng],gen_vy_array[ng],gen_vz_array[ng];
		int Is_Electron = -1;
				
		for(int ig=0;ig<ng;ig++){
			if((int)gen_pid->at(ig)==Electron)
				Is_Electron = ig;
			gen_mom_array[ig] = sqrt( pow(gen_px->at(ig),2)+pow(gen_py->at(ig),2)+pow(gen_pz->at(ig),2) ); //GeV
			gen_theta_array[ig] = acos(gen_pz->at(ig)/gen_mom_array[ig])*DEG;//Degree
			gen_phi_array[ig] = atan2( gen_py->at(ig), gen_px->at(ig))*DEG;//Degree
			gen_mom_array[ig] *=MeV2GeV ;
			gen_px_array[ig] = gen_px->at(ig)*MeV2GeV; //GeV
			gen_py_array[ig] = gen_py->at(ig)*MeV2GeV; //GeV
			gen_pz_array[ig] = gen_pz->at(ig)*MeV2GeV; //GeV
			gen_vx_array[ig] = gen_vx->at(ig)*MM2CM; //cm
			gen_vy_array[ig] = gen_vy->at(ig)*MM2CM; //cm
			gen_vz_array[ig] = gen_vz->at(ig)*MM2CM; //cm
			//    cerr<<Form("---#%d@%d Px=%f Py=%f Pz=%f P=%f Theta =%f Phi=%f ",i, ig,gen_px->at(ig),gen_py->at(ig),gen_pz->at(ig),gen_mom[ig],gen_theta[ig],gen_phi[ig])<<endl;
		}
		mom_gen = gen_mom_array[Is_Electron];theta_gen = gen_theta_array[Is_Electron];phi_gen = gen_phi_array[Is_Electron];
		px_gen = gen_px_array[Is_Electron];py_gen = gen_py_array[Is_Electron];pz_gen = gen_pz_array[Is_Electron];
		vx_gen = gen_vx_array[Is_Electron];vy_gen = gen_vy_array[Is_Electron];vz_gen = gen_vz_array[Is_Electron];

		//if(gen_theta[Is_Electron]>=5.0){
		if(1){
			FirstOne0 = 0;
			FirstOne1 = 0;

			flux->GetEntry(i);
			EC_Cut_Max = 0.0;
			nhit_flux = flux_hitn->size();
			/*Loop Flux Events{{{*/ 
			for (int j=0;j<nhit_flux;j++) {

				ID_Pick = VP_ID;
				ID_flux = (int) (flux_id->at(j));
				if(ID_flux!=ID_Pick) continue; //Only look at events on the detector plane (VP here) 

				PID_flux = (int) (flux_pid->at(j));
				//Select the right Events
				if(PID_flux!=Electron&&PID_flux!=-Electron
						&&PID_flux!=Pion&&PID_flux!=-Pion
						&&PID_flux!=Kaon&&PID_flux!=-Kaon
						&&PID_flux!=Proton&&PID_flux!=Photon) continue; //only deal with Electron, Proton, Pion,Kaon and Photon
				if(flux_pz->at(j)<1e-9)continue;//Cout out backward particles
				//	if(mom_flux<1e-9) continue; //Cut out Zero-E particles

				evn_flux = j;
				x_flux = flux_avg_x->at(j)*MM2CM; y_flux = flux_avg_y->at(j)*MM2CM; z_flux = flux_avg_z->at(j)*MM2CM;
				vx_flux = flux_vx->at(j)*MM2CM;   vy_flux = flux_vy->at(j)*MM2CM;	vz_flux = flux_vz->at(j)*MM2CM;
				px_flux = flux_px->at(j)*MeV2GeV; py_flux = flux_py->at(j)*MeV2GeV; pz_flux = flux_pz->at(j)*MeV2GeV;
				E = flux_trackE->at(j)*MeV2GeV;

				r_flux=sqrt(pow(x_flux,2)+pow(y_flux,2));//cm
				if(r_flux > R_Max_In||r_flux<R_Min_In) continue;//Must be within the radius range

				mom_flux=sqrt(pow(px_flux,2)+pow(py_flux,2)+pow(pz_flux,2));//GeV
				theta_flux = atan(sqrt(pow(px_flux,2)+pow(py_flux,2))/(pz_flux))*DEG;//Degree
				phi_flux = fabs(atan(y_flux/x_flux)*DEG);
				if(y_flux>0 && x_flux<0) phi_flux += 90.0;
				else if(y_flux<0 && x_flux<0) phi_flux += 180.0;
				else if(y_flux<0 && x_flux>0) phi_flux += 270.0;

				TID_flux = (int) (flux_tid->at(j));

				//Now looking into the SPD with this particular flux at VP
				/*Inside SPD{{{*/
				EDep = 0.00; convert = 0.0;
				for (int k=0;k<nhit_flux;k++) {
					if((int)(flux_id->at(k)!=SPD_ID)) continue;

					if((int)(flux_tid->at(k))!=TID_flux&&(int)(flux_mtid->at(k))!=TID_flux) continue;
					if((int)(flux_pid->at(k))!=PID_flux&&(int)(flux_mpid->at(k))!=PID_flux) continue;
				//	cerr<<"--- Fount one for TID = "<<TID_flux<<"---, new TID = "<<(int)(flux_tid->at(k))<<"---, MTID = "<<(int)(flux_mtid->at(k))<<endl;
					
					double x = flux_avg_x->at(k)*MM2CM, y = flux_avg_y->at(k)*MM2CM;
					double r=sqrt(pow(x,2)+pow(y,2));//cm
					if(r > R_Max_In||r<R_Min_In) continue;//Must be within the radius range

					EDep += flux_totEdep->at(k); //MeV

					int PID_temp = (int) (flux_pid->at(k));
					//Select the right Events
					if(PID_temp==Electron||PID_temp==-Electron
							||PID_temp==Pion||PID_temp==-Pion
							||PID_temp==Kaon||PID_temp==-Kaon
							||PID_temp==Proton) 
						//if there is any charged particles in the SPD that was from this particle at VP,
						//either from photon conversion or already being a charged particle, count as one
						convert = 1.0;
				}
				/*}}}*/

				EDep*=MeV2GeV; //GeV

				Module_ID = (int) phi_flux/(360./Module);
				for(int k=0;k<Slide+1;k++){
					if(r_flux>R_Slide[k]&&r_flux<=R_Slide[k+1])
						Slide_ID = k;
				}
				T->Fill();

				//Selct the right Cut	
				//cut[Det_ID][Cut_ID][Data_Point][Cut_Info]: Cut_Info: R_Min, R_Max, P_Min, P_Max, e_Eff, pi_Eff
				int EC_PID = 0;
				if(PID_flux==Electron||PID_flux==-Electron||PID_flux==Photon)
					EC_PID = 4;
				else if(PID_flux==Pion||PID_flux==-Pion||PID_flux==Kaon||PID_flux==-Kaon||PID_flux==Proton)//Only Pion's EC-R-Cut for Proton and Kaon so far
					EC_PID = 5;
				else
					cerr<<"*** Error, I don't know this particle = "<<PID_flux<<endl;

				double Delta_Z =Z_EC - Z_In;
				double x_EC = x_flux + Delta_Z * px_flux/pz_flux;
				double y_EC = y_flux + Delta_Z * py_flux/pz_flux;
				R_EC = sqrt(x_EC*x_EC+y_EC*y_EC);
				if(R_EC>R_Max_EC || R_EC<R_Min_EC) continue;//Project the particles to the EC and see whether it has been accepted there

				/*Look at the EC to check E-Cut{{{*/ 
				double EC_Cut_Temp = 0.0;
				for (unsigned int m=0;m<flux_hitn->size();m++) {
					int ID_flux_m = (int) (flux_id->at(m));
					//int PID_flux_m = (int) (flux_pid->at(m));

					if(ID_flux_m==VP_EC){//On EC VP
						double r_EC=sqrt(pow(flux_avg_x->at(m),2)+pow(flux_avg_y->at(m),2))/10.;//cm
						if(r_EC>R_Max_EC||r_EC<R_Min_EC) continue;//The radius of a mrpc sector is 210cm;

						//if(flux_pz->at(m)<1e-9)continue;//Cut out backward particles
						double fmom_EC=sqrt(pow(flux_px->at(m),2)+pow(flux_py->at(m),2)+pow(flux_pz->at(m),2))/1e3;//GeV
						if(fmom_EC<1e-9) continue; //Cut out Zero-E particles

						if(fabs(fmom_EC-mom_flux)/mom_flux<0.50){//Alow 50% Eloss, just to make the code faster
							/*Find E-Cut{{{*/
							for(int k=0;k<EC_Trigger_Slide;k++){
								for(int l=0;l<EC_Trigger_Mom_Bin;l++){
									if(r_EC>trig_cut[EC_ID][k][l][0]&&r_EC<=trig_cut[EC_ID][k][l][1]){
										if(fmom_EC>trig_cut[EC_ID][k][l][2]&&fmom_EC<=trig_cut[EC_ID][k][l][3]){
											EC_Cut_Temp =trig_cut[EC_ID][k][l][EC_PID]; 
										}
									}
								}
							}
							/*}}}*/
						}
					}
					if(EC_Cut_Temp>EC_Cut_Max)
						EC_Cut_Max = EC_Cut_Temp;
				}//for (Int_t m=0;m<flux_hitn->size();m++) 

				if(EC_Cut_Max<-1e-9||EC_Cut_Max>1){
					EC_Cut_Max = 0.0;
					cerr<<"----In: I can't find the cut!"<<Form(" --- R= %f,  E =%f", r_flux, mom_flux)<<endl;
					return "0";
				}
				/*}}}*/

				/*Find E-Cut{{{*/
				for(int k=0;k<EC_Trigger_Slide;k++){
					for(int l=0;l<EC_Trigger_Mom_Bin;l++){
						if(R_EC>trig_cut[EC_ID][k][l][0]&&R_EC<=trig_cut[EC_ID][k][l][1]){
							if(mom_flux>trig_cut[EC_ID][k][l][2]&&mom_flux<=trig_cut[EC_ID][k][l][3]){
								EC_Cut =trig_cut[EC_ID][k][l][EC_PID]; 
							}
						}
					}
				}
				/*}}}*/
			
				if(PID_flux==Photon){//#Photons going in
					hist_pho->Fill(y_flux, x_flux, EC_Cut*rate*convert);
				}
				if((PID_flux)==Electron||(PID_flux)==Pion||(PID_flux)==Kaon
						||(PID_flux)==-Electron||(PID_flux)==-Pion||(PID_flux)==-Kaon||PID_flux==Proton){//#Charged Particles going out 
					hist_crg->Fill(y_flux, x_flux, EC_Cut*rate);
				}

				if(PID_flux==Photon&& mom_flux>=0.5*EC_Threshold){//#Photons going in
					if(FirstOne0<1)
						In_G+=EC_Cut_Max*rate;
					FirstOne0++;
				}
				if((PID_flux==Electron||PID_flux==Pion||PID_flux==Kaon
						||PID_flux==-Electron||PID_flux==-Pion||PID_flux==-Kaon||PID_flux==Proton)
						&& mom_flux>=0.5*EC_Threshold){//#Charged Particles going out 
					if(FirstOne1<1)
						In_E+=EC_Cut_Max*rate;
					FirstOne1 ++;;
				}
			}
			/*flux}}}*/

		}//if(1) or if(gen_theta>5.0)
	}//gen
	/*End Read in each event}}}*/

	/*Count_Rate_Ractor{{{*/
	double Count_To_Rate = 1.0/1e3; //Hz to KHz
	cerr<<" ====== In_E = "<< In_E*Count_To_Rate<<endl;
	cerr<<" ====== In_G = "<< In_G*Count_To_Rate<<endl;
	/*}}}*/

	out_file->cd(); T->Write(); hist_crg->Write(); hist_pho->Write(); 
	out_file->Close(); 
	file->Close();
   
	return output_filename;
}

void RunAll(TString input_file){
   ifstream runlist(input_file.Data());
   TString filename = "";

   while(!(runlist.eof())){
      runlist >> filename;
      cerr<<"---- Working on the file = "<<filename.Data()<<endl;
	  TString output_file = FASPD(filename.Data());
//	  if(output_file!="0"){
		 // gSystem->Exec(Form("rm -f file_FASPD/%s", output_file.Data()));
		 // gSystem->Exec(Form("mv %s file_FASPD/", output_file.Data()));
//	  }
   }
   
}
