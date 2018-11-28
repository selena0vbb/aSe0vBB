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
	const char* file122kev = "./output/122_keV_Simple_Physicslarge.root";
	const char* file136kev = "./output/136_keV_Simple_Physicslarge.root";
	const char* file14kev = "./output/14_keV_Simple_Physicslarge.root";

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

void createCobaltTree(int N, const char* file122, const char* file136, const char* file14, const char* output="/home/apiers/data/particle/cobalt_spectrum.root"){

	TFile *outfile = new TFile(output, "RECREATE");

	// Read data files of particle data
	TFile *f122 = new TFile(file122); 
	TFile *f136 = new TFile(file136); 
	TFile *f14 = new TFile(file14);

	// Get trees from file
	TTree *t122 = (TTree*)f122->Get("aSeData"); 
	TTree *t136 = (TTree*)f136->Get("aSeData"); 
	TTree *t14 = (TTree*)f14->Get("aSeData"); 
	
	// Get number of gamma particles simulated
	int max122Events = (int)t122->GetMaximum("EventID") + 1;
	int max136Events = (int)t136->GetMaximum("EventID") + 1;
	int max14Events = (int)t14->GetMaximum("EventID") + 1;
	
	// Definition of intensities
	double inten122 = 0.856;
	double inten136 = 0.1068;
	double inten14 = 0.0916;

	// Define number of events to use from each energy spectrum 
	int n122 = (int) N*inten122;
	int n136 = (int) N*inten136;
	int n14 = (int) N*inten14;
	
	// Create new trees for the energies filtered by the number of events that come from Co-57 decay
	char *filter122 = new char[50];
	char *filter136 = new char[50];
	char *filter14 = new char[50];
	sprintf(filter122, "EventID < %i", n122);
	sprintf(filter136, "EventID < %i", n136);
	sprintf(filter14, "EventID < %i", n14);
	
	cout << filter122 << "   " << filter136 << "   " << filter14 << endl;	

	outfile->cd();
	TTree* t122Filter = t122->CopyTree(filter122);
	TTree* t136Filter = t136->CopyTree(filter136);
	TTree* t14Filter = t14->CopyTree(filter14);

	// Combine trees and write to file
	TList *list = new TList();
	list->Add(t122Filter); list->Add(t136Filter); list->Add(t14Filter); 
	TTree *coTree = TTree::MergeTrees(list);

	coTree->SetName("aSeData");

	// Create a new file and write the data
	coTree->Write();
	outfile->Close();

	return;

}
