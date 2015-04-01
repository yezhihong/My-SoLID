////////////////////////////////////////////
// Fill the 1D histo with weights
////////////////////////////////////////////
//You need to select electrons that hit on the GEM#i plane, then fill the distributions of x and y for each GEM
//For EM background only 
const double rate = ( (1.5e-5)/(1.6e-19) ) / nevent; //in Hz/cm2
//Fill the histo after events selections, please modify based on GetGEM.C
hx_gem1->Fill(y_gem1, rate);
hx_gem1->Fill(y_gem1, rate);
hr_gem1->Fill(r_gem1, rate);

///////////////////////////////////////////////
//Or easy way (no sure it works)
///
//TChain *flux = new TChain("flux");
//flux->Draw("avg_lx>>hx_gem1","id==1220000 && (pid==11||pid==-11)");//maybe more cuts to select the right event
//hx_gem1->Scale(rate); //You need to make sure it is normalized correctly, 
//TFile *file_histo = new TFile("background_histo.root","recreate");
//hx_gem1->Write(); file_histo->Close();

/////////////////////////////////////////
// Randomly generate background events
/////////////////////////////////////////
double Time_Windows = 75 * 1e-9; //nsec to sec
int Bin = (int)((r_flux-82.5)/1.0+0.5); //r_flux in cm, 58-bins from 82.5 to 140.5, each bin has 1.0cm width
double r_center = hr_gem1->GetBinCenter(Bin);
double Area = (PI*2*r_center*1.0); //PI* ( (r+1)*(r+1)-r*r ) = PI*2*(r+0.5)*1 
int Count = (int) (hr_gem1->GetBinContent(Bin) * Area * Time_Window + 0.5);//rate(Hz/cm2)*area(cm2)*time(sec)
double x_bgrd[Count],y_bgrd[Count];
for(int k=0;k<Count;k++){
  x_bgrd[k] = hx_gem1->GetRandom();
  y_bgrd[k] = hy_gem1->GetRandom();
 }
