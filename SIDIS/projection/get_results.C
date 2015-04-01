/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
//#include "Rtypes.h"
//#include "math.h"

/*}}}*/
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
#include <TMatrixD.h>
#include "LHAPDF/LHAPDF.h"
//#include <TMatrix.h>
/*}}}*/
//#include "solenoid_lib/Sole_LinkDef.h"
#include "solenoid_lib/Sole_inter.h"

using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;
const int N=3;

void getcoef(Double_t,Double_t*,Double_t*);
int get_results(Int_t target_flag,Int_t particle_flag,Int_t z_flag,Int_t Q2_flag);

/*int main{{{*/
int main(){

	Int_t target_flag = 3; cerr<<"--- Target (1->p, 3->he3) = "; cin >> target_flag;
	Int_t particle_flag = 2; cerr<<"--- particle (1->pip, 2->pim) = "; cin >> particle_flag;
	Int_t z_flag= 0;
	Int_t Q2_flag= 0;
	int err = -1000;
	for(int i=1;i<=8;i++){
		for(int j=1;j<=6;j++){
			z_flag = i;
			Q2_flag = j;
			err = get_results( target_flag, particle_flag, z_flag, Q2_flag);
		}
	}
	return err;
}
/*}}}*/

int get_results(Int_t target_flag,Int_t particle_flag,Int_t z_flag,Int_t Q2_flag){
	//   Int_t target_flag; //  1 for proton; 2 for deuteron; 3 for 3he
	//   Int_t Q2_flag; // 1 for Q2_flag <=10; 2 for Q2_flag>10;
	//   Int_t pt_flag; // 1-15
	//   Int_t z_flag; // 1-12
	//   Int_t particle_flag;
	//gSystem->Load("./solenoid_lib/libSole.so");

	/*Define In and Out{{{*/	
	TString filename,log_filename,data_filename;
	TString log_prefix,data_prefix;
	TString header,posfix;

	// loop throug the target_flag
	// for (particle_flag = 1; particle_flag<=2; particle_flag ++ ){
	if (particle_flag==1){
		data_prefix = "./out_rootfiles_p/";
	}else{
		data_prefix = "./out_rootfiles_m/";
	}    
	log_prefix = "./databases_both/";

	TString target = "X";
	if(target_flag==1)
		target ="p";
	else if(target_flag==2)
		target ="d2";
	else if(target_flag==3)
		target ="3he";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	TString particle = "X";
	if(particle_flag==1)
		particle ="pip";
	else if(particle_flag==2)
		particle ="pim";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
	}
	
	posfix.Form("%s_%s_%d_%d.dat",target.Data(),particle.Data(),z_flag,Q2_flag);
	log_filename = log_prefix + posfix;

	cerr<<Form("---- Reading %s", log_filename.Data())<<endl;

	ifstream infile(log_filename);

	ofstream outfile;
	filename.Form("./results/%s_%s_%d_%d.dat",target.Data(),particle.Data(),z_flag,Q2_flag);
	outfile.open(filename);
	cerr<<Form("---- Saving results to %s", filename.Data())<<endl;
	int pt_flag = -1;
	/*}}}*/

	Int_t ncount[2];  
	infile >> ncount[0];
    outfile << ncount[0]<<endl;
	Double_t counter[2][5000],Q2[5000],x[5000],z[5000],y[5000],pt[5000],nevent[5000],nevent_raw[5000],dilute[5000];
	Double_t x1[360],weight[360];
	Double_t angle[4];
	Double_t coef[3];
	Double_t coverage;
    Double_t x_sum = 0;
//	for (Int_t i=0;i<ncount[0];i++){
//		infile >> pt_flag  >> ncount[1];
//		if(pt_flag != i)
//			cerr<<Form("*** ERROR, something wrong here? pt=%d, i=%d",pt_flag,i)<<endl;
    int i=0;
	while(infile >> pt_flag  >> ncount[1]){

		posfix.Form("%s_%s_%d_%d_%d.root",target.Data(),particle.Data(),z_flag,Q2_flag,pt_flag);
		data_filename = data_prefix + posfix;
		//cerr<<"---- OK, let's look at histograms in File = "<< data_filename.Data()<<endl;
		TFile *file1=new TFile(data_filename,"r");

        /*X-bins{{{*/
	    bool JustOne = kTRUE;
		for (Int_t j=0;j<ncount[1];j++){
			infile >> counter[0][i*ncount[1]+j] //pt_bin, from 0->1.6 with 0.2 step
				>> counter[1][i*ncount[1]+j]    //x_bin, from 0.05 to 0.65
				>> z[i*ncount[1]+j] 
				>> Q2[i*ncount[1]+j] 
				>> pt[i*ncount[1]+j] 
				>> x[i*ncount[1]+j] 
				>> y[i*ncount[1]+j] 
				>> dilute[i*ncount[1]+j] 
				>> nevent_raw[i*ncount[1]+j]
				>> nevent[i*ncount[1]+j];

			if( counter[0][i*ncount[1]+j] !=pt_flag &&counter[1][i*ncount[1]+j]!=j) continue;
			if(nevent[i*ncount[1]+j]<1) continue;
	
		    if(JustOne)	
				outfile << pt_flag  << "\t" << ncount[1] << endl;
			JustOne=kFALSE;

			TH1F *h1 = (TH1F*) file1->Get(Form("h_%d_%d",i,j));

			x_sum =0;
			for (Int_t k=0;k<360;k++){
				x1[k] = h1->GetBinContent(k+1);
				x_sum += x1[k];
			}
			h1->Delete();
			if(x_sum <1) continue;

			sole_inter sole;
			sole.init(360,x1);
			sole.connect();
			sole.get_new(x1);
			sole.get_new(x1);
			sole.get_limit(angle);
			for (Int_t k=0;k!=360;k++){
				weight[k] = sole.get_accep(k+0.5);
			}
			coverage = (angle[1]-angle[0]+angle[3]-angle[2])*PI/180.;
			getcoef(coverage,coef,weight);

			// cout << i << "\t" << j << "\t" << coef[0] << endl;

			/*Output{{{*/
			if (nevent[i*ncount[1]+j]>0.){
				outfile << counter[0][i*ncount[1]+j] << " \t" 
					<< counter[1][i*ncount[1]+j] << " \t" 
					<< z[i*ncount[1]+j] << " \t" 
					<< Q2[i*ncount[1]+j] << " \t" 
					<< pt[i*ncount[1]+j] << " \t" 
					<< x[i*ncount[1]+j] << " \t" 
					<< y[i*ncount[1]+j] << " \t" 
					<< 1./sqrt(nevent[i*ncount[1]+j]) << " \t" 
					<< nevent_raw[i*ncount[1]+j] << " \t"
					<< coverage << " \t" 
					<< coef[0] << " \t" 
					<< coef[1] << " \t" 
					<< coef[2] << endl;
			}else{
				outfile << counter[0][i*ncount[1]+j] << " \t" 
					<< counter[1][i*ncount[1]+j] << " \t" 
					<< z[i*ncount[1]+j] << " \t" 
					<< Q2[i*ncount[1]+j] << " \t" 
					<< pt[i*ncount[1]+j] << " \t" 
					<< x[i*ncount[1]+j] << " \t" 
					<< y[i*ncount[1]+j] << " \t" 
					<< 0. << " \t" 
					<< 0. << " \t" 
					<< coverage << " \t" 
					<< coef[0] << " \t" 
					<< coef[1] << " \t" 
					<< coef[2] << endl;
			}
			/*}}}*/
		}
		/*}}}*/
		file1->Close();
        i++;
	}
	infile.close();
	outfile.close();
    return 0;
}

/*getcoef(Double_t angle, Double_t* coef, Double_t* acpt){{{*/
void getcoef(Double_t angle, Double_t* coef, Double_t* acpt){
	//x->phi y->phih 3 terms,w/ acpt
	Double_t a_matrix[N*N]={};
	Double_t x,y;
	for(Int_t i=0;i<180;i++){
		x = PI*(i+0.5)/180.; //phi

		for(Int_t j=0;j<360;j++){
			y = PI*(j+0.5)/180.; //phi_h
			a_matrix[0] += (sin(x)*sin(x))*PI/180.*PI/180.*acpt[j];
			a_matrix[1] += (sin(x)*sin(2*y-x))*PI/180.*PI/180.*acpt[j];
			a_matrix[2] += (sin(x)*sin(2*y+x))*PI/180.*PI/180.*acpt[j];

			a_matrix[N+1] += (sin(2*y-x)*sin(2*y-x))*PI/180.*PI/180.*acpt[j];
			a_matrix[N+2] += (sin(2*y-x)*sin(2*y+x))*PI/180.*PI/180.*acpt[j];

			a_matrix[2*N+2] += (sin(2*y+x)*sin(2*y+x))*PI/180.*PI/180.*acpt[j];
		}
	}
	a_matrix[N] = a_matrix[1];
	a_matrix[2*N] = a_matrix[2];
	a_matrix[2*N+1] = a_matrix[N+2];
	//       cout << "\t Initialize a_matrix done!" << endl;

	TMatrixD a(N,N,a_matrix);
	TMatrixD asav = a;
	a.Invert();
	//   cout << "\t Initialize a invert done!" << endl;

	Double_t qx[N]={};
	for(Int_t i=0;i<180;i++){
		x = PI*(i+0.5)/180.; //phi
		for(Int_t j=0;j<360;j++){
			y = PI*(j+0.5)/180.; //phi_h
			for(Int_t k=0;k<N;k++){
				qx[k] += (pow((a(k,0)*sin(x)+a(k,1)*sin(2*y-x)+a(k,2)*sin(2*y+x)),2))*PI/180.*PI/180.*acpt[j];
			}
		}
	}
	for(Int_t i=0;i<N;i++){
		qx[i] *= PI*angle;
		coef[i] = sqrt(qx[i]/2);
	}
	return;
}
/*}}}*/
