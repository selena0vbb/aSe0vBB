// ROOT macro file for plotting example B4 histograms 
// 
// Can be run from ROOT session:
// root[0] .x plotHisto.C
void plotHist(char* filename)
{
  gROOT->Reset();
  gROOT->SetStyle("Plain");

  // Import file
  TFile f(filename);
  f.ls();
  // Create Selenium energy deposit histogram
  TCanvas* c1 = new TCanvas("c1", "", 20, 20, 1000, 1000);
  TH1D* EneDepSe = (TH1D*)f.Get("EneDepSe");
  EneDepSe->SetDirectory(0);
  EneDepSe->SetLineWidth(3);
  string titleSe = string(filename) + " Selenium Energy Deposit";
  EneDepSe->SetTitle(titleSe.c_str());
  EneDepSe->Draw("HIST");

  // Create Gold Electrode energy deposit histogram
  TCanvas* c2 = new TCanvas("c2", "", 20, 20, 1000, 1000);
  TH1D* EneDepAu = (TH1D*)f.Get("EneDepAu");
  EneDepAu->SetDirectory(0);
  EneDepAu->SetLineWidth(3);
  string titleAu = string(filename) + " Electrode Energy Deposit";
  EneDepAu->SetTitle(titleAu.c_str());
  EneDepAu->Draw("HIST");

  
}  
