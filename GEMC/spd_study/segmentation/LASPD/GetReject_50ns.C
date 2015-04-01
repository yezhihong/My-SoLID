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

void GetReject(){

	const double MIP = 3.82; //MeV

	TChain *T = new TChain("T");
	T->Add("LASPD_pi0_p95_mpid_bck_50ns.root");
	T->Add("LASPD_pi0_up_p95_mpid_bck_50ns.root");
	T->Add("LASPD_pi0_down_p95_mpid_bck_50ns.root");

	double EDep_Sum, EDep,EDep_Bck, convert; 
    T->SetBranchAddress("EDep_Sum", &EDep_Sum);
    T->SetBranchAddress("EDep_Bck", &EDep_Bck);
    T->SetBranchAddress("EDep", &EDep);
    T->SetBranchAddress("convert", &convert);
    int N = T->GetEntries();


    TH2F *h1 = new TH2F("h1","",71, 0.,201., 100,0.0, 200.0);
    h1->SetXTitle("Segments");
    h1->SetYTitle("Rejection (%)");
	h1->SetMarkerStyle(21);
	h1->SetMarkerColor(4);

	ofstream outf("seg_50ns.out");

//	double Segment[70], Reject[70];
	double Segment, Reject;
	double Fire = 0;
	//int k = 0;
	for(int j=1;j<=200;j++){
		cerr<<"--- Work on Modules="<<j<<endl;

		//Segment[k] = j;
        Segment = j;

		Fire = 0.0;
		for(int i=0;i<N;i++){
			T->GetEntry(i);
			if((EDep+EDep_Bck+EDep_Sum/j)>MIP*0.5)
				Fire ++;
			/*
			//This is the part that one photon has the probability of converting into two electrons.
			//If the electrons have high enough energy (>10MeV), it should always fire the detecotr.
			if(EDep>MIP*0.5)
			Fire += convert;

			//This is the part that one phone has the change not to convert into electrons,
			//but rest of background charged particles still can deposit energy which will sum together
			// and pass the threshold
			if((EDep+EDep_Sum)/j>MIP*0.5)
			Fire += (1-convert);
			*/
			if(!(i%100000))
				cerr<<"--- Finished with events = "<<i<<endl;
		}
		Reject = (1.0*N)/(1.0*Fire); 
		outf <<Segment<<"   "<<N<<"   "<<"   "<<Fire<<"   "<<Reject<<endl;
		h1->Fill(Segment, Reject);
		cerr<<"--- Rejection = "<<Reject<<endl;
	}
    outf.close();
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
    h1->Draw();

}
