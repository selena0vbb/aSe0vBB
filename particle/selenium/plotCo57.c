// Root file for plotting the energy deposition for Co57 source (gamma)
// Combines energy deposition of several monoenergetic  based on relative intensity
// Documentation about relative intensity: https://www.nndc.bnl.gov/chart/decaysearchdirect.jsp?nuc=57CO&unc=nds


// run by root[0] .x plotCo57()
void plotCo57()
{
	// Definition of intensities
	double inten122 = 0.856;
	double inten136 = 0.1068;
	double inten14 = 0.0916;
	

	// Location where output data files are stored
	char* file122kev = "./output/122_keV_Simple_Physicslarge.root";
	char* file136kev = "./output/136_keV_Simple_Physicslarge.root";
	char* file14kev = "./output/14_keV_Simple_Physicslarge.root";

	// Load data
	TFile *f122 = new TFile(file122kev);
	TFile *f136 = new TFile(file136kev);
	TFile *f14 = new TFile(file14kev);

	// Create canvas
	TCanvas* c1 = new TCanvas("c1", "", 1000, 1000);

	// Create histograms
	TH1D * h122 = (TH1D*)f122->Get("EneDepSe");
	TH1D * h136 = (TH1D*)f136->Get("EneDepSe");
	TH1D * h14 = (TH1D*)f14->Get("EneDepSe");

	//TH1D* co57h = h122->Scale(inten122);
	h122->Scale(inten122);
	h136->Scale(inten136);
	h14->Scale(inten14);
	h122->Add(h136);
	h122->Add(h14);
	
	// Histogram Propertis
	h122->Draw("HIST");
	h122->SetLineWidth(4);
	h122->SetStats(false);		
	h122->SetTitle("Energy Deposited in Selenium by Gamma Rays of Co57 Source");
	
	// Update Canvas
	c1->Modified();
	c1->Update();
}