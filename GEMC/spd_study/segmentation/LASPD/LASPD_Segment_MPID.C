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

const double PI=3.1415926;
const double DEG=180./3.1415926;
const double nSec = 1e-9;//sec
const double Time_Window = 30 * nSec;

void LASPD_Segment(string input_filename){
//int main(){
	//string input_filename= "../background/background_pi0_window_downstream_filenum100.root";
	//string input_filename= "../background/background_pi0_window_upstream_filenum100.root";
//	string input_filename= "../background/background_pi0_filenum100.root";

	gStyle->SetOptStat(0);

	//Photon from Pi0
	gROOT->Reset();
	gStyle->SetPalette(1);
	gStyle->SetOptStat(0);
	
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
	}else if (input_filename.find("_clean_weighted_",0) != string::npos) {
		cout << "==> background from weighted event generator with no interaction except decay" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "==> filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "==> we need filenum for addtional normalization" << endl;}		
	}	  
	else if (input_filename.find("_dirty_normalized_",0) != string::npos) {
		cout << "==> background from normalized event generator with all interaction" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "==> filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "==> we need filenum for addtional normalization" << endl;;}
	}	  	
	else if (input_filename.find("background_",0) != string::npos) {
		cout << "==> background from normalized event generator with all interaction" <<  endl;
		if (input_filename.find("_file",0) != string::npos) {
			filenum=atof(input_filename.substr(input_filename.find("_filenum")+8,input_filename.find("_")).c_str());
			cout << "==> filenum " << filenum << " for addtional normalization, YOU Need to Make Sure It's CORRECT!" <<  endl;
		}
		else {cout << "==> we need filenum for addtional normalization" << endl;;}
	}	  	
	else {cout << "==> not EM or clean_weighted or dirty_normalized " << endl;}
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
	
	/*Define LASPD{{{*/
    TString Detector_Name = "LASPD";   
	const int SPD_VP_ID = 5210000; //Front virtual plane
	const double Z_In = -67.45; //cm Front VP position
	const double R_Min_In = 80.0;//cm, Front VP position
	const double R_Max_In = 135.0;//cm, Front VP position
	const int SPD_ID = 5200000; //Front virtual plane
	
	//const int VP_Out  = 3210000; //Rear virtual plane
	//const double Z_Out = -66.5; //cm Rear VP position
	//const double R_Min_Out = 83.0;//cm, Rear VP position
	//const double R_Max_Out = 140.0;//cm, Rear VP position

	const int VP_EC = 3210000; //Front virtual plane
	const double Z_EC = -66.5; //cm, EC virtual plane front, VP position = 413cm
	const double R_Min_EC = 83.0;//cm, EC VP position
	const double R_Max_EC = 140.0;//cm, EC VP position
	double EC_Threshold = 3.; //GeV for EC cut, the cut could be tight if using Jin's curves
    const int EC_ID = 1; //0->FAEC, 1->LAEC

	//For Segmentation only
	const int Module = 60; //30 module around the circle
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

	TString output_filename = Detector_Name +"_"+pid_temp+Form("_p95_mpid_bck_%dns.root", (int)(Time_Window/nSec));
	cerr<<"--- I will save the detector info into "<<output_filename.Data()<<endl;
	double rate, mom_gen, theta_gen,phi_gen,mom_flux,theta_flux,phi_flux,r_flux,x_flux,y_flux,z_flux,R_EC,EC_Cut,EC_Cut_Max;
	double px_flux, py_flux,pz_flux,vx_flux, vy_flux,vz_flux, E, EDep,EDep_Bck, EDep_Sum, convert;
	double px_gen, py_gen, pz_gen, vx_gen, vy_gen, vz_gen;
	int evn, evn_flux,nhit_flux, ID_flux,PID_flux,ID_Pick,TID_flux, MTID_flux;
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
	T->Branch("EDep_Bck",&EDep_Bck,"EDep_Bck/D");
	T->Branch("EDep_Sum",&EDep_Sum,"EDep_Sum/D");
	T->Branch("E",&E,"E/D");
	T->Branch("convert",&convert,"convert/D");
	T->Branch("EC_Cut",&EC_Cut,"EC_Cut/D");
	T->Branch("EC_Cut_Max",&EC_Cut_Max,"EC_Cut_Max/D");
	T->Branch("ID_flux",&ID_flux,"ID_flux/I");
	T->Branch("PID_flux",&PID_flux,"PID_flux/I");
	T->Branch("ID_Pick",&ID_Pick,"ID_Pick/I");
	T->Branch("TID_flux",&TID_flux,"TID_flux/I");
	T->Branch("MTID_flux",&MTID_flux,"MTID_flux/I");
	T->Branch("Module", &Module_ID,"Module/I");
	T->Branch("Slide", &Slide_ID,"Slide/I");
	/*}}}*/
				
	/*Background Tree and Histogram{{{*/
    TFile *f1 = new TFile("background_histo_EM.root","r");
    //TH1F *hE_ele = (TH1F*) f1->Get("hE_ele");
    //TH1F *hEdep_ele = (TH1F*) f1->Get("hEdep_ele");
    TH1F *hr_ele = (TH1F*) f1->Get("hr_ele");
    //TH1F *hx_ele = (TH1F*) f1->Get("hx_ele");
    //TH1F *hy_ele = (TH1F*) f1->Get("hy_ele");
	//TH2F *hRELog_ele = (TH2F*) f1->Get("hRELog_ele");
	//TH2F *hREdepLog_ele = (TH2F*) f1->Get("hREdepLog_ele");
	TH1F *hEdepLog_ele = (TH1F*) f1->Get("hEdepLog_ele");
	TH1F *hEdep_ele = (TH1F*) f1->Get("hEdep_ele");
	TH1F *hELogEdepLog_ele = (TH1F*) f1->Get("hELogEdepLog_ele");

    //TH1F *hE_pho = (TH1F*) f1->Get("hE_pho");
    //TH1F *hEdep_pho = (TH1F*) f1->Get("hEdep_pho");
    TH1F *hr_pho = (TH1F*) f1->Get("hr_pho");
    //TH1F *hx_pho = (TH1F*) f1->Get("hx_pho");
    //TH1F *hy_pho = (TH1F*) f1->Get("hy_pho");
	//TH2F *hRELog_pho = (TH2F*) f1->Get("hRELog_pho");
	//TH2F *hREdepLog_pho = (TH2F*) f1->Get("hREdepLog_pho");
	TH1F *hEdepLog_pho = (TH1F*) f1->Get("hEdepLog_pho");
	TH1F *hEdep_pho = (TH1F*) f1->Get("hEdep_pho");
	TH1F *hELogEdepLog_pho = (TH1F*) f1->Get("hELogEdepLog_pho");

	//              |--seconds--| |s-1|
	int Count_ele = (int) ( Time_Window * hr_ele->Integral()); 
	int Count_pho = (int) ( Time_Window * hr_pho->Integral());

	cerr<<"==================================================================================="<<endl; 
	cerr<<"--- Overal background events within 50 ns time-window for the entire SPD plane"<<endl;
	cerr<<Form("     electron = %d,  photon = %d",Count_ele, Count_pho)<<endl;  
	cerr<<"==================================================================================="<<endl; 
	/*}}}*/
	
	//TH1F* h_proj;

	/*Read in each event{{{*/
	cerr<<"++++++++++++++++ "<<endl;
	Int_t nselected = nevent;

	for(Int_t i=0;i<nselected;i++){
		cout<<i<<"\r";
        header->GetEntry(i);
		if(input_filename.find("_EM_",0) != string::npos) 
			rate = ((1.5e-5)/(1.6e-19))/nevent; //HZ
		else 
			rate = head_rate->at(0)/filenum; //HZ
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
			flux->GetEntry(i);
			EC_Cut_Max = 0.0;
			nhit_flux = flux_hitn->size();
			/*Loop Flux Events{{{*/ 
			for (int j=0;j<nhit_flux;j++) {
				evn_flux = j;
				x_flux = flux_avg_x->at(j)*MM2CM; y_flux = flux_avg_y->at(j)*MM2CM; z_flux = flux_avg_z->at(j)*MM2CM;
				vx_flux = flux_vx->at(j)*MM2CM;   vy_flux = flux_vy->at(j)*MM2CM;	vz_flux = flux_vz->at(j)*MM2CM;
				px_flux = flux_px->at(j)*MeV2GeV; py_flux = flux_py->at(j)*MeV2GeV; pz_flux = flux_pz->at(j)*MeV2GeV;
				E = flux_trackE->at(j);//MeV

				r_flux=sqrt(pow(x_flux,2)+pow(y_flux,2));//cm
				mom_flux=sqrt(pow(px_flux,2)+pow(py_flux,2)+pow(pz_flux,2));//GeV
				theta_flux = atan(sqrt(pow(px_flux,2)+pow(py_flux,2))/(pz_flux))*DEG;//Degree
				phi_flux = fabs(atan(y_flux/x_flux)*DEG);
				if(y_flux>0 && x_flux<0) phi_flux += 90.0;
				else if(y_flux<0 && x_flux<0) phi_flux += 180.0;
				else if(y_flux<0 && x_flux>0) phi_flux += 270.0;

				Module_ID = (int) phi_flux/(360./Module);
				for(int k=0;k<Slide+1;k++){
					if(r_flux>R_Slide[k]&&r_flux<=R_Slide[k+1])
						Slide_ID = k;
				}
				double Delta_Z =Z_EC - Z_In;
				double x_EC = x_flux + Delta_Z * px_flux/pz_flux;
				double y_EC = y_flux + Delta_Z * py_flux/pz_flux;
				R_EC = sqrt(x_EC*x_EC+y_EC*y_EC);

				ID_Pick = SPD_VP_ID;
				ID_flux = (int) (flux_id->at(j));
				PID_flux = (int) (flux_pid->at(j));
				TID_flux = (int) (flux_tid->at(j));
				MTID_flux = (int) (flux_mtid->at(j));

				//Select the right Events
				if(PID_flux!=Photon) continue; //only deal with Photon
				if(ID_flux!=ID_Pick) continue; //Only look at events on the detector plane (VP here) 
				if(r_flux > R_Max_In||r_flux<R_Min_In) continue;//Must be within the radius range
				//if(flux_pz->at(j)<1e-9)continue;//Cout out backward particles
				//if(mom_flux<1e-9) continue; //Cut out Zero-E particles

				//Now looking into the SPD with this particular flux at VP
				
				/*Inside SPD- forward{{{*/
				EDep = 0.0; convert = 0.0;
				for (int k=0;k<nhit_flux;k++) {
					if((int)(flux_id->at(k)!=SPD_ID)) continue;

					if((int)(flux_tid->at(k))!=TID_flux&&(int)(flux_mtid->at(k))!=TID_flux) continue;
					if((int)(flux_pid->at(k))!=PID_flux&&(int)(flux_mpid->at(k))!=PID_flux) continue;
					//	cerr<<"--- Fount one for TID = "<<TID_flux<<"---, new TID = "<<(int)(flux_tid->at(k))<<"---, MTID = "<<(int)(flux_mtid->at(k))<<endl;

					double x = flux_avg_x->at(k)*MM2CM, y = flux_avg_y->at(k)*MM2CM;
					double r=sqrt(pow(x,2)+pow(y,2));//cm
					if(r > R_Max_In||r_flux<R_Min_In) continue;//Must be within the radius range

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

				/*Inside SPD--backward{{{*/
				EDep_Bck = 0.0;
				for (int k=0;k<nhit_flux;k++) {
					if((int)(flux_id->at(k)!=SPD_ID)) continue;

            //---> 01/19: these condition includes the backward scattering particles from EC. 
			//A disadvantage is that if there are more than one pi0-photons or more high energy particles that also generate backward secondary particles, 
			// they will be wrongly included. 
					if((flux_pz->at(k))<-1e-10) continue;

					double x_temp = flux_avg_x->at(k)*MM2CM, y_temp = flux_avg_y->at(k)*MM2CM;
					double r_temp=sqrt(pow(x_temp,2)+pow(y_temp,2));//cm
					if(r_temp > R_Max_In||r_flux<R_Min_In) continue;//Must be within the radius range

					double phi_temp = fabs(atan(y_temp/x_temp)*DEG);
					if(y_temp>0 && x_temp<0) phi_temp += 90.0;
					else if(y_temp<0 && x_temp<0) phi_temp += 180.0;
					else if(y_temp<0 && x_temp>0) phi_temp += 270.0;

			// From pi0 sources, in general two photons won't come in to the same segmentation, so I apply a roughly big segmentation cut to avoid this situation.
					if(fabs(phi_flux-phi_temp)>(360./Module)) continue;

					EDep_Bck += flux_totEdep->at(k); //MeV

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
	
				/*Electron: Randomly get deposit energy based on its rate{{{*/
				double EDep_Sum_ele = 0.0;	
				for(int k=0;k<Count_ele;k++){
					double EDep_Temp = -0.0;
						
					EDep_Temp = hEdepLog_ele->GetRandom();//log(MeV) 
					if(EDep_Temp>-6.5) //I set EDep=0 to be 1e-7MeV to avoid bug when plotting in Log scale
						EDep_Sum_ele += pow(10,EDep_Temp);
					//EDep_Temp = hEdep_ele->GetRandom();//(MeV) 
					//EDep_Sum_ele += EDep_Temp;
/*
					h_proj = (TH1F*) hREdepLog_ele->ProjectionY("h_proj", RBin,RBin); 
					//assert(h_proj);
					if(h_proj->Integral()>1e-9)
						EDep_Temp = h_proj->GetRandom();//log10(MeV) 
					else
						EDep_Temp = 0.0;
*/
				}
				/*}}}*/

				/*Photon: Randomly get deposit energy based on its rate{{{*/
				//No additional slides along r-direction so I have to sum the entire plane
				double EDep_Sum_pho = 0.0;
				for(int k=0;k<Count_pho;k++){
					double EDep_Temp = 0.0;
					EDep_Temp = hEdepLog_pho->GetRandom();//log(MeV) 
					if(EDep_Temp>-6.5) //I set EDep=0 to be 1e-7MeV to avoid bug when plotting in Log scale
						EDep_Sum_pho += pow(10.,EDep_Temp);//MeV
					//EDep_Temp = hEdep_pho->GetRandom();//(MeV) 
					//EDep_Sum_pho += EDep_Temp;
/*
					h_proj = (TH1F*) hREdepLog_pho->ProjectionY("h_proj", RBin,RBin); 
					//assert(h_proj);
					if(h_proj->Integral()>1e-9)
						EDep_Temp = h_proj->GetRandom();//log(MeV) 
					else
						EDep_Temp = 0.0;
*/
				}
				/*}}}*/

				EDep_Sum = 0.0;
				//	int Bin = (int)((r_flux-75.)/1.0+0.5); //cm, 65-bins from 75. to 140.

				EDep_Sum = EDep_Sum_ele + EDep_Sum_pho;
				if(!(j%5000)){
					cerr<<Form("    Photon: EDep = %f, EDep_Back = %f", EDep, EDep_Bck)<<endl;
					cerr<<Form("Background: EDep = %f, from %f MeV(ele#%d)--%f MeV(pho#%d)", EDep_Sum , EDep_Sum_ele, Count_ele, EDep_Sum_pho,Count_pho)<<endl;
				}
				E*=MeV2GeV;
				T->Fill();
			}
			/*flux}}}*/

		}//if(1) or if(gen_theta>5.0)
	}//gen
	/*End Read in each event}}}*/

	out_file->cd(); T->Write(); out_file->Close(); 
	file->Close();
	f1->Close();
	//file_EDep->Close();
	//file_EDep_200MeV->Close();
	//file_EDep_10MeV->Close();
	//file_EDep_1MeV->Close();
	//delete h_proj;
	
}

int main(){
	string input_filename= "";
	input_filename= "../background/background_pi0_window_downstream_filenum100.root";
	LASPD_Segment(input_filename.c_str());
	input_filename= "../background/background_pi0_window_upstream_filenum100.root";
	LASPD_Segment(input_filename.c_str());
	input_filename= "../background/background_pi0_filenum100.root";
	LASPD_Segment(input_filename.c_str());

}
