{
	TChain *T = new TChain("T");
	T->Add("./FASPD_background_SIDIS_He3_pi0_down_mpid.root");
	T->Add("./FASPD_background_SIDIS_He3_pi0_mpid.root");
	T->Add("./FASPD_background_SIDIS_He3_pi0_up_mpid.root");

    TCanvas *c2 = new TCanvas("c2","c2",1200,1000);
	c2->Divide(1,2);
	c2->cd(1);
	TH1F *h1 = new TH1F("h1","e^{#pm} distribution on FASPD", 100,90,245);
	h1->SetXTitle("r_{flux}");
    T->Draw("r_flux>>h1","rate*(abs(PID_flux)==11)","");

	c2->cd(2);
	TH1F *h2 = new TH1F("h2","#gamma distribution on FASPD", 100,90,245);
	h2->SetXTitle("r_{flux}");
    T->Draw("r_flux>>h2","rate*(PID_flux==22)","");


}

