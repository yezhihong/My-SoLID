{
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	
	TFile *output = new TFile("SPD_EDep_200MeV.root","recreate");

	TCanvas *c1 = new TCanvas("c1","c1",1200,800);
	TFile *fE1 = new TFile("electron/5mm/SIDIS_SPD_E200MeV_out.root");
	TFile *fE2 = new TFile("electron/1cm/SIDIS_SPD_E200MeV_out.root");
	TFile *fE3 = new TFile("electron/2cm/SIDIS_SPD_E200MeV_out.root");
//	TFile *fE4 = new TFile("electron/5cm/SIDIS_SPD_E200MeV_out.root");
	
	TFile *fG3 = new TFile("photon/2cm/SIDIS_SPD_G200MeV_out.root");

 
	c1->cd(1); fE1->cd();
	TH1F *hE1 = new TH1F("hE1","5mm LASPD (electron)",100,1e-9, 20.0);
	hE1->SetLineColor(2);
	hE1->SetLineWidth(1.5);
	hE1->SetLineStyle(1);
	T->Draw("Edep>>hE1","");
	hE1->Scale(1/1000.);
	
	c1->cd(1); fE2->cd();
	TH1F *hE2 = new TH1F("hE2","1cm LASPD (electron)",100,1e-9, 20.0);
	hE2->SetLineColor(4);
	hE2->SetLineWidth(1.5);
	hE2->SetLineStyle(1);
	T->Draw("Edep>>hE2","");
	hE2->Scale(1/1000.);

	c1->cd(1); fE3->cd();
	TH1F *hE3 = new TH1F("hE3","2cm LASPD (electron)",100,1e-9, 20.0);
	hE3->SetLineColor(6);
	hE3->SetLineWidth(1.5);
	hE3->SetLineStyle(1);
	T->Draw("Edep>>hE3","");
	hE3->Scale(1/1000.);

	c1->cd(1);
	gPad->SetLogy(1);
	fE1->cd(); hE1->Draw();
	fE2->cd(); hE2->Draw("same");
	fE3->cd(); hE3->Draw("same");
/*
	c1->cd(1); fE4->cd();
	TH1F *hE4 = new TH1F("hE4","5cm LASPD (electron)",400,1e-9, 20.0);
	hE4->SetLineColor(6);
	hE4->SetLineWidth(1.5);
	hE4->SetLineStyle(1);

	T->Draw("Edep>>hE4","");
	c1->cd(1);fE4->cd(); hE4->Draw("same");
*/

	TCanvas *c2 = new TCanvas("c2","c2",1200,600);
	c2->Divide(4,1);
	c2->cd(1); fE1->cd();gPad->SetLogy(1);
	TH2F *hEDep5mm = new TH2F("hEDep5mm","5mm SPD with Electron (0-200 MeV)",100,1e-9, 200,100, 1e-1, 20.0);
    T->Draw("Edep:Ee>>hEDep5mm","","colz");
	c2->cd(2); fE2->cd();gPad->SetLogy(1);
	TH2F *hEDep1cm = new TH2F("hEDep1cm","1cm SPD with Electron (0-200 MeV)",100,1e-9, 200,100, 1e-1, 20.0);
    T->Draw("Edep:Ee>>hEDep1cm","","colz");
	c2->cd(3); fE3->cd();gPad->SetLogy(1);
	TH2F *hEDep2cm = new TH2F("hEDep2cm","2cm SPD with Electron (0-200 MeV)",100,1e-9, 200,100, 1e-1, 20.0);
    T->Draw("Edep:Ee>>hEDep2cm","","colz");
//	c2->cd(4); fE4->cd();gPad->SetLogy(1);
//	TH2F *hEDep5cm = new TH2F("hEDep5cm","50mm SPD with Electron (0-200 MeV)",100,1e-9, 200,100, 1e-1, 20.0);
  //  T->Draw("Edep:Ee>>hEDep5cm","","colz");

	c2->cd(4); fG3->cd();gPad->SetLogy(1);
	TH2F *hGDep2cm = new TH2F("hGDep2cm","2cm SPD with Photon (0-200 MeV)",100,1e-9, 200,100, 1e-1, 20.0);
    T->Draw("Edep:Eg>>hGDep2cm","","colz");
	

	TCanvas *c3 = new TCanvas("c3","c3",1200,600);
	c3->Divide(4,1);
	c3->cd(1); fE1->cd();gPad->SetLogy(0);
	TH1F* hEDep5mm_px = (TH1F*) hEDep5mm->ProfileX("hEDep5mm_px"); hEDep5mm_px->Draw();
  	c3->cd(2); fE2->cd();gPad->SetLogy(0);
	TH1F* hEDep1cm_px = (TH1F*) hEDep1cm->ProfileX("hEDep1cm_px"); hEDep1cm_px->Draw();
	c3->cd(3); fE3->cd();gPad->SetLogx(1);
;gPad->SetLogy(1);
	TH1F* hEDep2cm_px = (TH1F*) hEDep2cm->ProfileX("hEDep2cm_px"); hEDep2cm_px->Draw();
//	c3->cd(4); fE4->cd();gPad->SetLogy(0);
//	TH1F* hEDep5cm_px = (TH1F*) hEDep5cm->ProfileX("hEDep5cm_px"); hEDep2cm_px->Draw();
 
	c3->cd(4); fG3->cd();gPad->SetLogy(0);
	TH1F* hGDep2cm_px = (TH1F*) hGDep2cm->ProfileX("hGDep2cm_px"); hGDep2cm_px->Draw();
	 
   c2->Print("SPD_EDep_2D_200MeV.png");
   c3->Print("SPD_EDep_1D_200MeV.png");

	output->cd();
	hE1->Write();
	hE2->Write();
	hE3->Write();
//	hE4->Write(); 

    hEDep5mm->Write(); hEDep5mm_px->Write();
	hEDep1cm_px->Write(); hEDep1cm->Write();
	hEDep2cm->Write(); hEDep2cm_px->Write();
//	hEDep5cm->Write(); hEDep5cm_px->Write();
	hGDep2cm->Write(); hGDep2cm_px->Write();

	output->Close();

}
