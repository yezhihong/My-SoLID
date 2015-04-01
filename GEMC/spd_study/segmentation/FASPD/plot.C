{

	gStyle->SetOptStat(0);

    TH2F *h1 = new TH2F("h1","10mm FASPD Rejection Factor with 30ns/50ns Time Window",71, 0.0,801., 100,0.0, 20.0);
    h1->SetXTitle("Segment");
    h1->SetYTitle("Rejection Factor");
	h1->SetMarkerStyle(21);
	h1->SetMarkerColor(4);

    TH2F *h2 = new TH2F("h2","10mm FASPD Rejection Rate with 30ns/50ns Time Window",71, 0.0,801., 100,0.0, 100.0);
    h2->SetXTitle("Segment");
    h2->SetYTitle("Rejection Rate (%)");
	h2->SetMarkerStyle(21);
	h2->SetMarkerColor(4);

 	TH2F *k1 = new TH2F("k1","10mm FASPD Rejection Factor with 30ns/50ns Time Window",71, 0.0,801., 100,0.0, 20.0);
    k1->SetXTitle("Segment");
    k1->SetYTitle("Rejection Factor");
	k1->SetMarkerStyle(22);
	k1->SetMarkerColor(2);

    TH2F *k2 = new TH2F("k2","10mm FASPD Rejection Rate with 30ns/50ns Time Window",71, 0.0,801., 100,0.0, 100.0);
    k2->SetXTitle("Segment");
    k2->SetYTitle("Rejection Rate (%)");
	k2->SetMarkerStyle(22);
	k2->SetMarkerColor(2);

	ifstream inf1("seg_30ns.out");
	ifstream inf2("seg_50ns.out");
	int Seg = 0, N = 0, Fire = 0; 
	double Rej = 0.0;
	double Rej_Rate = 0.0;
	int Dummy;
	while(!(inf1.eof())){
		inf1 >> Seg >> N >> Fire >> Rej;
		cerr<<Form("---- Seg = %d, Rej = %f", Seg, Rej)<<endl;
		h1->Fill(Seg, Rej);
		Rej_Rate = 1.*(N-Fire)/N*100.0;
		h2->Fill(Seg, Rej_Rate);
	}
	inf1.close();

	while(!(inf2.eof())){
		inf2 >> Seg >>N >> Fire >> Rej;
		cerr<<Form("---- Seg = %d, Rej = %f", Seg, Rej)<<endl;
		k1->Fill(Seg, Rej);
		Rej_Rate = 1.*(N-Fire)/N*100.0;
		k2->Fill(Seg, Rej_Rate);
	}
	inf2.close();
	
	TCanvas *c1 = new TCanvas("c1","c1",800,600);
	h1->Draw();
	k1->Draw("same");


	TCanvas *c2 = new TCanvas("c2","c2",800,600);
	h2->Draw();
	k2->Draw("same");
}
