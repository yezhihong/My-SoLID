{
	TChain *T0 = new TChain("T");
	T0->Add("LASPD_background_SIDIS_He3_eDIS.root");
    
	TCanvas *c1 = new TCanvas("c1","c1",800,800);
	TH2F *h0 = new TH2F("h0","e^{#pm} (DIS) distribution on LASPD", 100,-145,145,100,-145,145);
	h0->SetXTitle("x_{flux}");
	h0->SetYTitle("y_{flux}");
    T0->Draw("y_flux:x_flux>>h0","rate","colz");

	TChain *T = new TChain("T");
	T->Add("LASPD_background_SIDIS_He3_pi0_down.root");
	T->Add("LASPD_background_SIDIS_He3_pi0.root");
	T->Add("LASPD_background_SIDIS_He3_pi0_up.root");
	T->Add("LASPD_background_SIDIS_He3_pim.root");
	T->Add("LASPD_background_SIDIS_He3_pim_up.root");
	T->Add("LASPD_background_SIDIS_He3_pim_down.root");
	T->Add("LASPD_background_SIDIS_He3_pip.root");
	T->Add("LASPD_background_SIDIS_He3_pip_up.root");
	T->Add("LASPD_background_SIDIS_He3_pip_down.root");
	T->Add("LASPD_background_SIDIS_He3_p.root");
	T->Add("LASPD_background_SIDIS_He3_p_up.root");
	T->Add("LASPD_background_SIDIS_He3_p_down.root");
	
    TCanvas *c2 = new TCanvas("c2","c2",600,1200);
    c2->Divide(1,2);
	c2->cd(1);	
	TH2F *h1 = new TH2F("h1","e^{#pm} (#pi^{0}+#pi^{#pm}+p+EM) distribution on LASPD", 100,-145,145,100,-145,145);
	h1->SetXTitle("x_{flux}");
	h1->SetYTitle("y_{flux}");
    T->Draw("y_flux:x_flux>>h1","rate*(PID_flux==11)","colz");
	c2->cd(2);	
	TH2F *h2 = new TH2F("h2","#gamma (#pi^{0}+#pi^{#pm}+p+EM) distribution on LASPD", 100,-145,145,100,-145,145);
	h2->SetXTitle("x_{flux}");
	h2->SetYTitle("y_{flux}");
    T->Draw("y_flux:x_flux>>h2","rate*(PID_flux==22)","colz");

}

