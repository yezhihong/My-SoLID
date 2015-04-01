void plot0(){

	gStyle->SetOptStat(0);

	TFile *fE3 = new TFile("electron/2cm/SIDIS_SPD_E_out.root");
	TFile *fG3 = new TFile("photon/2cm/SIDIS_SPD_G_out.root");

	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(2,2);
	
	c3->cd(1); fE3->cd();gPad->SetLogy(1);
	TH2F *hEDep2cm = new TH2F("hEDep2cm","2cm SPD with Electron (0-11 GeV)",110,1e-9, 11000,110, 1e-1, 20.0);
    T->Draw("Edep:Ee>>hEDep2cm","","colz");

	c3->cd(2); fG3->cd();gPad->SetLogy(1);
	TH2F *hGDep2cm = new TH2F("hGDep2cm","2cm SPD with Photon (0-11 GeV)",110,1e-9, 11000,110, 1e-1, 20.0);
    T->Draw("Edep:Eg>>hGDep2cm","","colz");

	TH2F *h1_2cm = new TH2F("h1_2cm","2cm SPD with Electron (0-11 GeV)",110,1e-9, 11000,110, 1e-1, 5.0);
	TH2F *h2_2cm = new TH2F("h2_2cm","2cm SPD with Photon (0-11 GeV)",110,1e-9, 11000,110, 1e-1, 5.0);
	
	c3->cd(3); fE3->cd();gPad->SetLogy(0);
	TH1F* hEDep2cm_px = (TH1F*) hEDep2cm->ProfileX("hEDep2cm_px"); hEDep2cm_px->Draw();
	h1_2cm->Draw(); hEDep2cm_px->Draw("same");
	hEDep2cm_px->Fit("pol0","","same");
  
	c3->cd(4); fG3->cd();gPad->SetLogy(0);
	TH1F* hGDep2cm_px = (TH1F*) hGDep2cm->ProfileX("hGDep2cm_px"); hGDep2cm_px->Draw();
	h2_2cm->Draw(); hGDep2cm_px->Draw("same");
	hGDep2cm_px->Fit("pol0","","same");
}

void plot1(){
	gStyle->SetOptStat(0);

	TFile *fE3 = new TFile("electron/2cm/SIDIS_SPD_E200MeV_out.root");
	TFile *fG3 = new TFile("photon/2cm/SIDIS_SPD_G200MeV_out.root");

	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(2,2);
	
	c3->cd(1); fE3->cd();gPad->SetLogy(1);
	TH2F *hEDep2cm = new TH2F("hEDep2cm","2cm SPD with Electron (0-200 MeV)",100,1e-9, 200.0,100, 1e-1, 20.0);
    T->Draw("Edep:Ee>>hEDep2cm","","colz");

	c3->cd(2); fG3->cd();gPad->SetLogy(1);
	TH2F *hGDep2cm = new TH2F("hGDep2cm","2cm SPD with Photon (0-200 MeV)",100,1e-9, 200.0,100, 1e-1, 20.0);
    T->Draw("Edep:Eg>>hGDep2cm","","colz");

	c3->cd(3); fE3->cd();gPad->SetLogy(0);
	TH1F* hEDep2cm_px = (TH1F*) hEDep2cm->ProfileX("hEDep2cm_px"); hEDep2cm_px->Draw();
  
	c3->cd(4); fG3->cd();gPad->SetLogy(0);
	TH1F* hGDep2cm_px = (TH1F*) hGDep2cm->ProfileX("hGDep2cm_px"); hGDep2cm_px->Draw();
}

void plot2(){
	gStyle->SetOptStat(0);

	TFile *fE3 = new TFile("electron/2cm/SIDIS_SPD_E10MeV_out.root");
	TFile *fG3 = new TFile("photon/2cm/SIDIS_SPD_G10MeV_out.root");

	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(2,2);
	
	c3->cd(1); fE3->cd();gPad->SetLogy(1);
	TH2F *hEDep2cm = new TH2F("hEDep2cm","2cm SPD with Electron (0-10 MeV)",100,1e-9, 10.0,100, 1e-1, 20.0);
    T->Draw("Edep:Ee>>hEDep2cm","","colz");

	c3->cd(2); fG3->cd();gPad->SetLogy(1);
	TH2F *hGDep2cm = new TH2F("hGDep2cm","2cm SPD with Photon (0-10 MeV)",100,1e-9, 10.0,100, 1e-1, 20.0);
    T->Draw("Edep:Eg>>hGDep2cm","","colz");

	c3->cd(3); fE3->cd();gPad->SetLogy(0);
	TH1F* hEDep2cm_px = (TH1F*) hEDep2cm->ProfileX("hEDep2cm_px"); hEDep2cm_px->Draw();
  
	c3->cd(4); fG3->cd();gPad->SetLogy(0);
	TH1F* hGDep2cm_px = (TH1F*) hGDep2cm->ProfileX("hGDep2cm_px"); hGDep2cm_px->Draw();
}

void plot3(){
	gStyle->SetOptStat(0);

	TFile *fE3 = new TFile("electron/2cm/SIDIS_SPD_E1MeV_out.root");
	TFile *fG3 = new TFile("photon/2cm/SIDIS_SPD_G1MeV_out.root");

	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(2,2);
	
	c3->cd(1); fE3->cd();gPad->SetLogy(1);
	TH2F *hEDep2cm = new TH2F("hEDep2cm","2cm SPD with Electron (0-1 MeV)",100,1e-9, 1.0,100, 1e-2, 20.0);
    T->Draw("Edep:Ee>>hEDep2cm","","colz");

	c3->cd(2); fG3->cd();gPad->SetLogy(1);
	TH2F *hGDep2cm = new TH2F("hGDep2cm","2cm SPD with Photon (0-1 MeV)",100,1e-9, 1.0,100, 1e-2, 20.0);
    T->Draw("Edep:Eg>>hGDep2cm","","colz");

	c3->cd(3); fE3->cd();gPad->SetLogy(0);
	TH1F* hEDep2cm_px = (TH1F*) hEDep2cm->ProfileX("hEDep2cm_px"); hEDep2cm_px->Draw();
  
	c3->cd(4); fG3->cd();gPad->SetLogy(0);
	TH1F* hGDep2cm_px = (TH1F*) hGDep2cm->ProfileX("hGDep2cm_px"); hGDep2cm_px->Draw();
}
