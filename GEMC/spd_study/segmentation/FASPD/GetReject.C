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

	const double MIP = 0.40; //MeV for 5mm SPD

	int Time_Window = 0; 
	cerr<<"-- Time Window (30 or 50) ns = "; cin>> Time_Window;

	TChain *T = new TChain("T");
	T->Add(Form("FASPD_pi0_p95_mpid_bck_%dns.root",Time_Window));
	//T->Add(Form("FASPD_pi0_up_p95_mpid_bck_%dns.root",Time_Window));
	//T->Add(Form("FASPD_pi0_down_p95_mpid_bck_%dns.root",Time_Window));

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

	ofstream outf(Form("seg_%dns_new.out",Time_Window));

	double  Reject = 0.0;
	int Segment = 0, Fire1= 0, Fire2= 0, Fire3= 0;
	for(int j=10;j<=100;j++){
		cerr<<"--- Work on Modules="<<j<<endl;

        Segment = 4*j;
		Fire1 = 0;
		Fire2 = 0;
		Fire3 = 0;

		for(int i=0;i<N;i++){
			T->GetEntry(i);
			
			//double EDep_Total = EDep+EDep_Bck+EDep_Sum/Segment;
			double EDep_Total = EDep;
			if(EDep_Total>MIP*0.5)
				Fire1++;
	
			EDep_Total = EDep+EDep_Bck;
			if(EDep_Total>MIP*0.5)
				Fire2++;
		
			EDep_Total = EDep+EDep_Bck+EDep_Sum/Segment;
			if(EDep_Total>MIP*0.5)
				Fire3++;
	
	
			if(!(i%100000))
				cerr<<Form("--- Finished with events = %d (EDep_Total = %f (%f), Seg = %d,  Fire = %d)",i, EDep_Total, MIP*0.5, Segment, Fire3)<<endl;
		}

		Reject = (1.0*N)/(1.0*Fire3); 
		outf <<Form("%5d %10d %10d %10d %10d %10.4f",
				      Segment, N, Fire1, Fire2, Fire3, Reject)<<endl;
		h1->Fill(Segment, Reject);
		cerr<<"--- Rejection = "<<Reject<<endl;
	}
    outf.close();
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
    h1->Draw();

}
